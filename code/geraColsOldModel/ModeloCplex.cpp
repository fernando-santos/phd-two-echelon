#include "ModeloCplex.h"
using namespace std;

double ModeloCplex::limitePrimal;
ModeloCplex::~ModeloCplex(){
	for (int s = 1; s <= nSatellites; ++s)
	{
		for (int r = 0; r < qRotasL2[s]; ++r)
		{
			delete [] A_qr[s][r];
			if ( ptrRotasL2[s][r]->decrNumApontadores() ) delete ptrRotasL2[s][r];
		}
		A_qr[s].clear();
		ptrRotasL2[s].clear();
	}
	delete [] A_qr;
	delete [] ptrRotasL2;
	delete [] qRotasL2;
	env.end();

	delete [] edgesHead;
	delete [] edgesTail;
	delete [] edgesX;
	delete [] demands;
}


ModeloCplex::ModeloCplex( Grafo* g, vector < Rota* > rLambda, vector < Rota* > rGamma, int separateCuts ) : 
						env(), model(env), cplex(model), modelCGH(env), cplexCGH(modelCGH), sepCuts(separateCuts) {
	limitePrimal = MAIS_INFINITO;
	cplex.setOut(env.getNullStream());
	cplexCGH.setOut(env.getNullStream());

	maxVeicL1 = g->getMaxVeiculosL1();
	maxVeicL2 = g->getMaxVeiculosL2();
	maxVeicSat = g->getMaxVeiculosSatellites();

	nSatellites = g->getNumSatellites();
	nCustomers = g->getNumCustomers();
	nVertices = g->getNumVertices();

	setA_qr(rGamma);

	//Montagem do modelo Primal
	initVars( g, rLambda );
	setObjectiveFunction( rLambda );
	setConstraints1( g, rLambda );
	setConstraint2();
	setConstraints3();
	setConstraints4();
	setConstraints5( g );
	setConstraints6( g, rLambda );
	setConstraints7( g );
	capacityCuts = IloRangeArray( env );
}


void ModeloCplex::setA_qr(vector <Rota*> rGamma){
	A_qr = new vector<short int*>[nSatellites+1];
	ptrRotasL2 = new vector<Rota*>[nSatellites+1];
	qRotasL2 = new int[nSatellites+1];
	for ( int s = 0; s <= nSatellites; ++s ) qRotasL2[s] = 0;

	vector<int> v;
	short int* rota;
	int tam;

	//Cada linha de A_qr apontara para um vetor de inteiros de tamanho numCommodities
	//Cada posicao do vetor sera 1 caso o vertice i esteja na rota r, ou 0 caso contrario
	for ( int r = 0; r < rGamma.size(); ++r )
	{
		rota = new short int[nVertices+1];
		memset(rota, 0, (nVertices+1) * sizeof(short int));
		v = rGamma[r]->getVertices();
		tam = v.size()-1;
		for( int j = 1; j < tam; ++j )
		{
			++rota[v[j]]; //coloca-se 1 apenas naquelas posicoes em que o vertice esteja na rota
		}
		A_qr[v[0]].push_back(rota);
		ptrRotasL2[v[0]].push_back(rGamma[r]);
		rGamma[r]->incrNumApontadores();
		++qRotasL2[v[0]];
	}
}


void ModeloCplex::exportModel(const char* arqExport){
	cplex.exportModel(arqExport);
}


double ModeloCplex::solveMaster(){
	cplex.solve();
	return cplex.getObjValue();
}


void ModeloCplex::initVars( Grafo* g, vector < Rota* > rLambda ){
	lambda = IloArray < IloNumVarArray > (env, maxVeicL1);
	lambdaCGH = IloArray < IloIntVarArray > (env, maxVeicL1);
	for ( int k = 0; k < maxVeicL1; ++k )
	{
		lambda[k] = IloNumVarArray( env, rLambda.size(), 0, 1, ILOFLOAT );
		lambdaCGH[k] = IloIntVarArray( env, rLambda.size(), 0, 1 );
		for ( int r = 0; r < rLambda.size(); ++r )
		{
			char nome[20];
			sprintf(nome, "lbd_%01d_%02d", k, r);
			lambda[k][r].setName(nome);
		}
	}

	gamma = IloArray < IloNumVarArray > (env, nSatellites+1);
	gammaCGH = IloArray < IloIntVarArray > (env, nSatellites+1);
	gamma_noCGH = IloArray < IloIntVarArray > (env, nSatellites+1);
	for ( int s = 1; s <= nSatellites; ++s )
	{
		gamma[s] = IloNumVarArray(env, qRotasL2[s], 0, 1, ILOFLOAT );
		gammaCGH[s] = IloIntVarArray(env, qRotasL2[s], 0, 1 );
		gamma_noCGH[s] = IloIntVarArray(env);
		for ( int r = 0; r < qRotasL2[s]; ++r )
		{
			char nome[20];
			sprintf(nome, "gam_%01d_%02d", s, r);
			gamma[s][r].setName(nome);
		}
	}

	int limiteDelta = 0;
	for (int i = 1; i <= nCustomers; ++i) limiteDelta += g->getCargaVertice(nSatellites + i);
	
	delta = IloArray < IloNumVarArray > (env, nSatellites+1);
	deltaCGH = IloArray < IloIntVarArray > (env, nSatellites+1);
	for ( int s = 1; s <= nSatellites; ++s )
	{
		delta[s] = IloNumVarArray( env, maxVeicL1, 0, limiteDelta, ILOFLOAT );
		deltaCGH[s] = IloIntVarArray( env, maxVeicL1, 0, limiteDelta );
		for (int k = 0; k < maxVeicL1; ++k)
		{
			char nome[20];
			sprintf(nome, "del_%01d_%01d", s, k);
			delta[s][k].setName(nome);
		}
	}

	//variaveis necessarias para executar o procedimento de separacao de desigualdades da biblioteca CVRPSEP
	CMGR_CreateCMgr(&newCuts, 100);
	CMGR_CreateCMgr(&oldCuts, 100);
	edgesHead = new int[nCustomers*nCustomers];
	edgesTail = new int[nCustomers*nCustomers];
	edgesX = new double[nCustomers*nCustomers];
	demands = new int[nCustomers+1];
	for ( int i = 1; i <= nCustomers; ++i) demands[i] = g->getCargaVertice(nSatellites + i);
}


void ModeloCplex::setObjectiveFunction( vector < Rota* > rLambda ){
	IloExpr obj(env), objCGH(env);

	for ( int k = 0; k < maxVeicL1; ++k )
	{
		for ( int r = 0; r < rLambda.size(); ++r )
		{
			obj += rLambda[r]->getCusto() * lambda[k][r];
			objCGH += rLambda[r]->getCusto() * lambdaCGH[k][r];
			
		}
	}
		
	for ( int s = 1; s <= nSatellites; ++s )
	{
		for ( int r = 0; r < qRotasL2[s]; ++r )
		{
			obj += ptrRotasL2[s][r]->getCusto() * gamma[s][r];
			objCGH += ptrRotasL2[s][r]->getCusto() * gammaCGH[s][r];
		}
	}

	objCost = IloObjective(env, obj);
	objCostCGH = IloObjective(env, objCGH);

	model.add(objCost);
	modelCGH.add(objCostCGH);
}


void ModeloCplex::setConstraints1( Grafo* g, vector < Rota* > rLambda ){
	IloExpr expCap = IloExpr(env);
	IloExpr expCapCGH = IloExpr(env);
	constraints1 = IloRangeArray(env, maxVeicL1);
	constraints1CGH = IloRangeArray(env, maxVeicL1);

	for ( int k = 0; k < maxVeicL1; ++k )
	{
		IloExpr exp1 = IloExpr(env);
		IloExpr exp1CGH = IloExpr(env);
		for ( int r = 0; r < rLambda.size(); ++r )
		{
			exp1 += lambda[k][r];
			exp1CGH += lambdaCGH[k][r];

			expCap += lambda[k][r];
			expCapCGH += lambdaCGH[k][r];
		}
		char nome[20];
		sprintf(nome, "rest1_k%d", k);
		constraints1[k] = (exp1 <= 1);
		constraints1[k].setName(nome);
		model.add(constraints1[k]);

		constraints1CGH[k] = (exp1CGH <= 1);
		modelCGH.add(constraints1CGH[k]);
		exp1.end(); exp1CGH.end();
	}
	
	//inclui uma restricao para apertar o modelo
	double cargaTotal = 0;
	for ( int i = 1; i <= nCustomers; ++i ) cargaTotal += g->getCargaVertice(nSatellites + i);
	model.add ( expCap >= ceil ( cargaTotal / g->getCapacVeiculosL1() ) );
	modelCGH.add ( expCapCGH >= ceil ( cargaTotal / g->getCapacVeiculosL1() ) );
}


void ModeloCplex::setConstraint2(){
	IloExpr exp2 = IloExpr(env);
	IloExpr exp2CGH = IloExpr(env);
	
	for ( int s = 1; s <= nSatellites; ++s )
	{
		for ( int r = 0; r < qRotasL2[s]; ++r )
		{
			exp2 += gamma[s][r];
			exp2CGH += gammaCGH[s][r];
		}
	}

	constraint2 = (exp2 <= maxVeicL2);
	char nome[20];
	sprintf(nome, "rest2");
	constraint2.setName(nome);
	model.add(constraint2);

	constraint2CGH = (exp2CGH <= maxVeicL2);
	modelCGH.add( constraint2CGH );

	exp2.end();
	exp2CGH.end();
}


void ModeloCplex::setConstraints3(){
	constraints3 = IloRangeArray(env, nSatellites);
	constraints3CGH = IloRangeArray(env, nSatellites);

	for ( int s = 1; s <= nSatellites; ++s)
	{
		IloExpr exp3 = IloExpr(env);
		IloExpr exp3CGH = IloExpr(env);
		for ( int r = 0; r < qRotasL2[s]; ++r )
		{
			exp3 += gamma[s][r];
			exp3CGH += gammaCGH[s][r];
		}

		constraints3[s-1] = (exp3 <= maxVeicSat);
		char nome[20];
		sprintf(nome, "rest3_s%01d", s);
		constraints3[s-1].setName(nome);
		model.add(constraints3[s-1]);

		constraints3CGH[s-1] = (exp3CGH <= maxVeicSat);
		modelCGH.add(constraints3CGH[s-1]);

		exp3.end();
		exp3CGH.end();
	}
}


void ModeloCplex::setConstraints4(){
	constraints4 = IloRangeArray(env, nCustomers);
	constraints4CGH = IloRangeArray(env, nCustomers);
	int count = 0;

	for ( int i = nSatellites + 1; i < nVertices; ++i )
	{
		IloExpr exp4 = IloExpr(env);
		IloExpr exp4CGH = IloExpr(env);
		for ( int s = 1; s <= nSatellites; ++s )
		{
			for ( int r = 0; r < qRotasL2[s]; ++r )
			{
				if ( A_qr[s][r][i] > 0 )
				{
					exp4 += gamma[s][r];
					exp4CGH += gammaCGH[s][r];
				}
			}
		}

		constraints4[count] = (exp4 == 1);
		char nome[20];
		sprintf(nome, "rest4_c%02d", i);
		constraints4[count].setName(nome);
		model.add(constraints4[count]);
		
		constraints4CGH[count] = (exp4CGH >= 1);
		modelCGH.add(constraints4CGH[count]);
		
		exp4.end();
		exp4CGH.end();
		++count;
	}
}


void ModeloCplex::setConstraints5( Grafo* g ){
	constraints5 = IloRangeArray(env, nSatellites);
	constraints5CGH = IloRangeArray(env, nSatellites);
	int coeficiente5;

	for (int s = 1; s <= nSatellites; ++s )
	{
		IloExpr exp5 = IloExpr(env);
		IloExpr exp5CGH = IloExpr(env);
		for ( int r = 0; r < qRotasL2[s]; ++r )
		{
			coeficiente5 = 0;
			for ( int i = nSatellites + 1; i < nVertices; ++i)
			{
				if ( A_qr[s][r][i] > 0 ) coeficiente5 += g->getCargaVertice( i );
			}
			exp5 += coeficiente5 * gamma[s][r];
			exp5CGH += coeficiente5 * gammaCGH[s][r];
		}
		for (int k = 0; k < maxVeicL1; ++k)
		{
			exp5 -= delta[s][k];
			exp5CGH -= deltaCGH[s][k];
		}

		constraints5[s-1] = (exp5 == 0);
		char nome[20];
		sprintf(nome, "rest5_s%01d", s);
		constraints5[s-1].setName(nome);
		model.add(constraints5[s-1]);
		
		constraints5CGH[s-1] = (exp5CGH == 0);
		modelCGH.add(constraints5CGH[s-1]);
		
		exp5.end();
		exp5CGH.end();
	}
}

void ModeloCplex::setConstraints6( Grafo* g, vector < Rota* > rLambda ){
	constraints6 = IloRangeArray(env, nSatellites*maxVeicL1);
	constraints6CGH = IloRangeArray(env, nSatellites*maxVeicL1);
	vector < int > v;
	int count = 0;
	
	for (int s = 1; s <= nSatellites; ++s )
	{
		for (int k = 0; k < maxVeicL1; ++k)
		{
			IloExpr exp6 = delta[s][k];
			IloExpr exp6CGH = deltaCGH[s][k];
			for ( int r = 0; r < rLambda.size(); ++r )
			{
				v = rLambda[r]->getVertices();
				for ( int i = 1; i < v.size(); ++i )
				{
					if ( v[i] == s )
					{
						exp6 -= g->getCapacVeiculosL1() * lambda[k][r];
						exp6CGH -= g->getCapacVeiculosL1() * lambdaCGH[k][r];
						break;
					}
				}
			}

			constraints6[count] = (exp6 <= 0);
			char nome[20];
			sprintf(nome, "rest6_s%01d_k%01d", s, k);
			constraints6[count].setName(nome);
			model.add(constraints6[count]);
			
			constraints6CGH[count] = (exp6CGH <= 0);
			modelCGH.add(constraints6CGH[count]);

			exp6.end();
			exp6CGH.end();
			++count;
		}
	}
}


void ModeloCplex::setConstraints7( Grafo* g ){
	constraints7 = IloRangeArray(env, maxVeicL1);
	constraints7CGH = IloRangeArray(env, maxVeicL1);
	
	for (int k = 0; k < maxVeicL1; ++k)
	{
		IloExpr exp7 = IloExpr(env);
		IloExpr exp7CGH = IloExpr(env);
		for (int s = 1; s <= nSatellites; ++s )
		{
			exp7 += delta[s][k];
			exp7CGH += deltaCGH[s][k];
		}

		constraints7[k] = ( exp7 <= g->getCapacVeiculosL1() );
		char nome[20];
		sprintf(nome, "rest7_k%01d", k);
		constraints7[k].setName(nome);
		model.add(constraints7[k]);
		
		constraints7CGH[k] = ( exp7CGH <= g->getCapacVeiculosL1() );
		modelCGH.add(constraints7CGH[k]);
	}
}


void ModeloCplex::updateDualCosts( Grafo* g, int s ){
	double custoDual;
	IloNumArray muDuals(env), roDuals(env);
	cplex.getDuals(muDuals, constraints4);
	cplex.getDuals(roDuals, capacityCuts);
	double piDual = cplex.getDual(constraints5[s-1]);
	int numEdgesCut, totalCuts = capacityCutsLHS.size();

	for (int i = 0; i < nCustomers; ++i)
	{
		custoDual = -muDuals[i] - g->getCargaVertice(nSatellites + i + 1)*piDual;
		g->setCustoVerticeDual(nSatellites + i + 1, custoDual);
	}
	g->setCustoArestasDual(s);

	for ( int i = 0; i < totalCuts; ++i )
	{
		numEdgesCut = capacityCutsLHS[i].size();
		if ( capacityCutsLHS[i][0] > 0 )
		{
			for ( int j = 0; j < numEdgesCut; j+=2 )
			{
				g->setCustoArestaDual(capacityCutsLHS[i][j]-nSatellites, capacityCutsLHS[i][j+1]-nSatellites, -roDuals[i]);
				g->setCustoArestaDual(capacityCutsLHS[i][j+1]-nSatellites, capacityCutsLHS[i][j]-nSatellites, -roDuals[i]);
			}
		}
		else
		{

			int j, tam = -capacityCutsLHS[i][0];
			for ( j = 1; j <= tam; ++j )
			{
				g->setCustoArestaDual(0, capacityCutsLHS[i][j]-nSatellites, -0.5*roDuals[i]);
				g->setCustoArestaDual(capacityCutsLHS[i][j]-nSatellites, 0, -0.5*roDuals[i]);
			}
			while ( capacityCutsLHS[i][j] >= 0 )
			{
				g->setCustoArestaDual(capacityCutsLHS[i][j]-nSatellites, capacityCutsLHS[i][j+1]-nSatellites, -roDuals[i]);
				g->setCustoArestaDual(capacityCutsLHS[i][j+1]-nSatellites, capacityCutsLHS[i][j]-nSatellites, -roDuals[i]);
				j += 2;
			}
			while ( j < numEdgesCut )
			{
				g->setCustoArestaDual(0, (-capacityCutsLHS[i][j])-nSatellites, 0.5*roDuals[i]);
				g->setCustoArestaDual((-capacityCutsLHS[i][j])-nSatellites, 0, 0.5*roDuals[i]);
				++j;
			}
		}
	}
	muDuals.end(); roDuals.end();
}


double ModeloCplex::getValorSubtrair( Grafo* g, int s ){
	double somaXi = 0, somaPsi = 0;
	for ( int k = 0; k < maxVeicL1; ++k )
	{
		somaXi += cplex.getDual( constraints6[ ( ( s-1 ) * maxVeicL1 ) + k ] );
		somaPsi += cplex.getDual( constraints7[ k ] );
	}
	return cplex.getDual( constraint2 ) + cplex.getDual( constraints3[ s-1 ] ) + (somaXi + somaPsi);
}


void ModeloCplex::insertColumn( Grafo* g, Rota* r, int s ){
	//Armazena a rota na matriz de rotas, que sera usada posteriormente
	short int* rota = new short int[nVertices];
	memset(rota, 0, nVertices*sizeof(short int));
	A_qr[s].push_back(rota);
	ptrRotasL2[s].push_back(r);
	r->incrNumApontadores();

	vector <int> vertRota = r->getVertices();
	int numVertAtualizar = vertRota.size()-2;
	for(int j = 1; j <= numVertAtualizar; ++j)
	{
		++A_qr[s][qRotasL2[s]][vertRota[j]]; //+1 apenas nas posicoes que o vertice esteja na rota
	}

	//CRIA-SE UMA COLUNA ASSOCIADA A ROTA PASSADA COMO PARAMETRO E INSERE NA VARIAVEL ASSOCIADA AO SATELLITE S
	IloNumColumn col = objCost( r->getCusto() );
	IloNumColumn colCGH = objCostCGH( r->getCusto() );

	col += constraint2(1);
	colCGH += constraint2CGH(1);

	col += constraints3[s-1](1);
	colCGH += constraints3CGH[s-1](1);
	int coeficiente5 = 0;

	for ( int i = 0; i < nCustomers; ++i )
	{
		if ( A_qr[s][qRotasL2[s]][nSatellites + i + 1] > 0 )
		{
			col += constraints4[i]( A_qr[s][qRotasL2[s]][nSatellites + i + 1] );
			colCGH += constraints4CGH[i]( A_qr[s][qRotasL2[s]][nSatellites + i + 1] );
			coeficiente5 += ( A_qr[s][qRotasL2[s]][nSatellites + i + 1] * g->getCargaVertice( nSatellites + i + 1 ) );
		}
	}
	col += constraints5[s-1](coeficiente5);
	colCGH += constraints5CGH[s-1](coeficiente5);

	float indexCol;
	int numEdgesCut, numCuts = capacityCutsLHS.size();
	for ( int i = 0; i < numCuts; ++i )
	{
		indexCol = 0;
		numEdgesCut = capacityCutsLHS[i].size();
		if ( capacityCutsLHS[i][0] >= 0 )
		{
			for ( int j = 0; j < numEdgesCut; j+=2 )
			{
				for ( int vR = 1; vR < numVertAtualizar; ++vR )
				{
					if ( ( ( vertRota[vR] == capacityCutsLHS[i][j] ) && ( vertRota[vR+1] == capacityCutsLHS[i][j+1] ) ) || ( ( vertRota[vR] == capacityCutsLHS[i][j+1] ) && ( vertRota[vR+1] == capacityCutsLHS[i][j] ) ) )
					{
						++indexCol;
					}
				}
			}
			if ( indexCol > 0 ) col += capacityCuts[i](indexCol);
		}
		else
		{
			int j, tam = -capacityCutsLHS[i][0];
			for ( j = 1; j <= tam; ++j )
			{
				if ( vertRota[1] == capacityCutsLHS[i][j] ) indexCol += 0.5;
				if ( vertRota[numVertAtualizar] == capacityCutsLHS[i][j] ) indexCol += 0.5;
			}
			while ( capacityCutsLHS[i][j] >= 0 )
			{
				for ( int vR = 1; vR < numVertAtualizar; ++vR )
				{
					if ( ( ( vertRota[vR] == capacityCutsLHS[i][j] ) && ( vertRota[vR+1] == capacityCutsLHS[i][j+1] ) ) || ( ( vertRota[vR] == capacityCutsLHS[i][j+1] ) && ( vertRota[vR+1] == capacityCutsLHS[i][j] ) ) )
					{
						++indexCol;
					}
				}
				j += 2;
			}
			while ( j < numEdgesCut )
			{
				if ( vertRota[1] == -capacityCutsLHS[i][j] ) indexCol -= 0.5;
				if ( vertRota[numVertAtualizar] == -capacityCutsLHS[i][j] ) indexCol -= 0.5;
				++j;
			}
			if ( ( indexCol > 0.0001 ) || ( indexCol < 0.0001 ) ) col += capacityCuts[i](indexCol);
		}
	}

	char nome[20];
	sprintf(nome, "gam_%d_%d", s, qRotasL2[s]);
	gamma[s].add(IloNumVar(col, 0, 1, ILOFLOAT, nome));
	gammaCGH[s].add(IloIntVar(colCGH, 0, 1, nome));

	col.end();
	colCGH.end();
	++qRotasL2[s];
}


bool ModeloCplex::addCapacityCuts( Grafo* g, int* qRotas_no, vector < Rota* >* ptrRotas_no )
{
	int numEdges;
	double maxViolate;
	char solIntegerFeasible;

	buildGraphSolution(numEdges, edgesHead, edgesTail, edgesX, qRotas_no, ptrRotas_no);
	CAPSEP_SeparateCapCuts(nCustomers, demands, g->getCapacVeiculosL2(), numEdges, edgesTail, edgesHead, edgesX, oldCuts, 100, 0.0001, &solIntegerFeasible, &maxViolate, newCuts);

	if ( solIntegerFeasible ) return false; /* Optimal solution found */
	if ( newCuts->Size == 0 ) return false; /* No cuts found */

	double RHS;
	vector<int> vertices;
	short int* sBarra = new short int[1+nCustomers/*/2*/];

	int head, tail, numVertices, indexCCC = capacityCutsLHS.size();
	for ( int i = 0; i < newCuts->Size; ++i )
	{
		vector<int> edgesCut;
		IloExpr expCCC = IloExpr(env);

		if ( newCuts->CPL[i]->IntListSize < nCustomers )
		{
			if ( newCuts->CPL[i]->IntListSize <= nCustomers ) /* <= ( nCustomers / 2 ) */
			{
				for ( int j = 1; j < newCuts->CPL[i]->IntListSize; ++j )
				{
					head = ( newCuts->CPL[i]->IntList[j] + nSatellites );
					for ( int k = j+1; k <= newCuts->CPL[i]->IntListSize; ++k )
					{
						tail = ( newCuts->CPL[i]->IntList[k] + nSatellites );
						edgesCut.push_back( head );
						edgesCut.push_back( tail );

						for ( int s = 1; s <= nSatellites; ++s )
						{
							for ( int r = 0; r < qRotasL2[s]; ++r )
							{
								vertices = ptrRotasL2[s][r]->getVertices();
								numVertices = vertices.size()-2;
								for ( int e = 1; e < numVertices; ++e )
								{
									if ( ( ( vertices[e] == head ) && ( vertices[e+1] == tail ) ) || ( ( vertices[e] == tail ) && ( vertices[e+1] == head ) ) )
									{
										expCCC += gamma[s][r];
									}
								}
							}

							for ( int r = 0; r < qRotas_no[s]; ++r )
							{
								vertices = ptrRotas_no[s][r]->getVertices();
								numVertices = vertices.size()-2;
								for ( int e = 1; e < numVertices; ++e )
								{
									if ( ( ( vertices[e] == head ) && ( vertices[e+1] == tail ) ) || ( ( vertices[e] == tail ) && ( vertices[e+1] == head ) ) )
									{
										expCCC += gamma_no[s][r];
									}
								}
							}
						}
					}
				}
				RHS = newCuts->CPL[i]->RHS;
			}
			else
			{
				//Define o complemento do conjunto S (sBarra) e seu tamanho (tamSBarra)
				int c, tamSBarra = 0;
				for ( int v = 1; v <= nCustomers; ++v )
				{
					for ( c = 1; c <= newCuts->CPL[i]->IntListSize; ++c )
					{
						if ( newCuts->CPL[i]->IntList[c] == v ) break;
					}
					if ( c > newCuts->CPL[i]->IntListSize ) sBarra[tamSBarra++] = v;
				}
				edgesCut.push_back( -tamSBarra );

				//Especifica o corte em funcao dos elementos de S e sBarra: x(sBarra:sBarra) + 0.5*x({0}:sBarra) - 0.5*x({0}:S) <= |sBarra| - \pi(S)
				//Arranjo do vetor edgesCut: O primeiro termo representa a cardinalidade de sBarra. A seguir vem os VERTICES de sBarra para os quais exitem arestas x({0}:sBarra).
				//Os proximos elementos sao as ARESTAS x(sBarra:sBarra) e finalmente tem os VERTICES de S para os quais existem as arestas x({0}:sBarra) - estes elementos
				//tambem sao armazenados com coeficiente negativo
				//OU SEJA, NO MESMO CORTE SAO ARMAZENADOS OS VERTICES DE SBARRA, AS ARESTAS ENTRE OS ELEMENTOS DE SBARRA E OS VERTICES DE S (OS VERTICES DE S TEM COEFICIENTE NEGATIVO)
				for ( int j = 0; j < tamSBarra; ++j ) // + 0.5*x({0}:sBarra)
				{
					tail = sBarra[j] + nSatellites;
					edgesCut.push_back( tail );

					for ( int s = 1; s <= nSatellites; ++s )
					{
						for ( int r = 0; r < qRotasL2[s]; ++r )
						{
							vertices = ptrRotasL2[s][r]->getVertices();
							if ( vertices[1] == tail ) expCCC += 0.5*gamma[s][r];
							if ( vertices[vertices.size()-2] == tail ) expCCC += 0.5*gamma[s][r];
						}

						for ( int r = 0; r < qRotas_no[s]; ++r )
						{
							vertices = ptrRotas_no[s][r]->getVertices();
							if ( vertices[1] == tail ) expCCC += 0.5*gamma_no[s][r];
							if ( vertices[vertices.size()-2] == tail ) expCCC += 0.5*gamma_no[s][r];
						}
					}
				}

				for ( int j = 0; j < tamSBarra-1; ++j ) // x(sBarra:sBarra)
				{
					head = ( sBarra[j] + nSatellites );
					for ( int k = j+1; k < tamSBarra; ++k )
					{
						tail = ( sBarra[k] + nSatellites );
						edgesCut.push_back( head );
						edgesCut.push_back( tail );

						for ( int s = 1; s <= nSatellites; ++s )
						{
							for ( int r = 0; r < qRotasL2[s]; ++r )
							{
								vertices = ptrRotasL2[s][r]->getVertices();
								numVertices = vertices.size()-2;
								for ( int e = 1; e < numVertices; ++e )
								{
									if ( ( ( vertices[e] == head ) && ( vertices[e+1] == tail ) ) || ( ( vertices[e] == tail ) && ( vertices[e+1] == head ) ) )
									{
										expCCC += gamma[s][r];
									}
								}
							}

							for ( int r = 0; r < qRotas_no[s]; ++r )
							{
								vertices = ptrRotas_no[s][r]->getVertices();
								numVertices = vertices.size()-2;
								for ( int e = 1; e < numVertices; ++e )
								{
									if ( ( ( vertices[e] == head ) && ( vertices[e+1] == tail ) ) || ( ( vertices[e] == tail ) && ( vertices[e+1] == head ) ) )
									{
										expCCC += gamma_no[s][r];
									}
								}
							}
						}
					}
				}

				for ( int j = 1; j <= newCuts->CPL[i]->IntListSize; ++j ) // - 0.5*x({0}:S)
				{
					tail = newCuts->CPL[i]->IntList[j] + nSatellites;
					edgesCut.push_back( -tail );
					for ( int s = 1; s <= nSatellites; ++s )
					{
						for ( int r = 0; r < qRotasL2[s]; ++r )
						{
							vertices = ptrRotasL2[s][r]->getVertices();
							if ( vertices[1] == tail ) expCCC -= 0.5*gamma[s][r];
							if ( vertices[vertices.size()-2] == tail ) expCCC -= 0.5*gamma[s][r];
						}

						for ( int r = 0; r < qRotas_no[s]; ++r )
						{
							vertices = ptrRotas_no[s][r]->getVertices();
							if ( vertices[1] == tail ) expCCC -= 0.5*gamma_no[s][r];
							if ( vertices[vertices.size()-2] == tail ) expCCC -= 0.5*gamma_no[s][r];
						}
					}
				}
				RHS = ( tamSBarra - ( newCuts->CPL[i]->IntListSize - newCuts->CPL[i]->RHS ) );
			}

			capacityCutsLHS.push_back(edgesCut);
			capacityCuts.add(expCCC <= RHS);
			model.add(capacityCuts[indexCCC]);
			expCCC.end();
			++indexCCC;
		}
	}

	for ( int i = 0; i < newCuts->Size; ++i )
	{
		CMGR_MoveCnstr(newCuts, oldCuts, i, 0);
	}
	newCuts->Size = 0;
	delete [] sBarra;
	return true;
}


void ModeloCplex::buildGraphSolution( int &edgesSize, int* edgesHead, int* edgesTail, double* edgesX, int* qRotas_no, vector < Rota* >* ptrRotas_no )
{
	edgesSize = 0;
	IloNumArray values(env);
	vector< int > verticesRota;
	int e, head, tail, sizeRota;

	for ( int s = 1; s <= nSatellites; ++s )
	{
		cplex.getValues(values, gamma[s]);
		for (int r = 0; r < qRotasL2[s]; ++r)
		{
			if (values[r] >= 0.0001)
			{
				verticesRota = ptrRotasL2[s][r]->getVertices();
				sizeRota = verticesRota.size()-1;
				for ( int i = 0; i < sizeRota; ++i )
				{
					head = ( verticesRota[i] > nSatellites ) ? (verticesRota[i]-nSatellites) : nCustomers+1;
					tail = ( verticesRota[i+1] > nSatellites ) ? (verticesRota[i+1]-nSatellites) : nCustomers+1;

					for ( e = 1; e <= edgesSize; ++e )
					{
						if ( ( ( edgesHead[e] == head ) && ( edgesTail[e] == tail ) ) || ( (edgesHead[e] == tail) && ( edgesTail[e] == head ) ) )
						{
							edgesX[e] += values[r];
							break;
						}
					}
					if ( e > edgesSize ) //a aresta ainda nao foi registrada no conjunto de arestas
					{
						++edgesSize;
						edgesHead[edgesSize] = head;
						edgesTail[edgesSize] = tail;
						edgesX[edgesSize] = values[r];
					}
				}
			}
		}

		if ( qRotas_no[s] > 0 ) cplex.getValues(values, gamma_no[s]);
		for ( int r = 0; r < qRotas_no[s]; ++r )
		{
			if (values[r] >= 0.0001)
			{
				verticesRota = ptrRotas_no[s][r]->getVertices();
				sizeRota = verticesRota.size()-1;
				for ( int i = 0; i < sizeRota; ++i )
				{
					head = ( verticesRota[i] > nSatellites ) ? (verticesRota[i]-nSatellites) : nCustomers+1;
					tail = ( verticesRota[i+1] > nSatellites ) ? (verticesRota[i+1]-nSatellites) : nCustomers+1;

					for ( e = 1; e <= edgesSize; ++e )
					{
						if ( ( ( edgesHead[e] == head ) && ( edgesTail[e] == tail ) ) || ( (edgesHead[e] == tail) && ( edgesTail[e] == head ) ) )
						{
							edgesX[e] += values[r];
							break;
						}
					}
					if ( e > edgesSize ) //a aresta ainda nao foi registrada no conjunto de arestas
					{
						++edgesSize;
						edgesHead[edgesSize] = head;
						edgesTail[edgesSize] = tail;
						edgesX[edgesSize] = values[r];
					}
				}
			}
		}
	}
	values.end();
}


void ModeloCplex::setLimitePrimal(double lP){
	limitePrimal = lP;
}


double ModeloCplex::getLimitePrimal(){
	return limitePrimal;
}


void ModeloCplex::imprimeRotasBasicas( vector < Rota* > rLambda, Grafo* g ){
	IloNumArray values(env);

	printf("\n");
	double totalLambda = 0;
	for ( int k = 0; k < maxVeicL1; ++k )
	{
		cplex.getValues(values, lambda[k]);
		for (int r = 0; r < values.getSize(); ++r)
		{
			if (values[r] >= 0.0001)
			{
				printf("lambda[%d][%d](%0.03f): ", k, r, values[r]);
				rLambda[r]->imprimir();
				totalLambda += values[r];
			}
		}
	}

	double totalGamma = 0;
	for ( int s = 1; s <= nSatellites; ++s )
	{
		cplex.getValues(values, gamma[s]);
		for (int r = 0; r < qRotasL2[s]; ++r)
		{
			if (values[r] >= 0.0001)
			{
				printf("gamma[%d][%d](%0.03f): ", s, r, values[r]);
				ptrRotasL2[s][r]->imprimir();
				totalGamma += values[r];
			}
		}
	}
	
	for ( int s = 1; s <= nSatellites; ++s )
	{
		cplex.getValues(values, delta[s]);
		for ( int k = 0; k < maxVeicL1; ++k )
		{
			if (values[k] >= 0.0001)
			{
				printf("delta[%d][%d](%0.02f)\n", s, k, values[k]);
			}
		}
	}

	double cargaTotal = 0;
	for ( int i = 1; i <= nCustomers; ++i ) cargaTotal += g->getCargaVertice(nSatellites + i);
	printf( "totalLambda = %0.02f\n", totalLambda );
	printf( "totalGamma = %0.02f\n", totalGamma );
	printf( "cargaTotal = %0.02f\n", cargaTotal );
	printf( "cargaTotal / capacL1 = %0.02f\n", cargaTotal / g->getCapacVeiculosL1() );
	printf( "cargaTotal / capacL2 = %0.02f\n\n", cargaTotal / g->getCapacVeiculosL2() );


	for ( int s = 1; s <= nSatellites; ++s )
	{
		cplex.getValues(values, gamma[s]);
		for (int r = 0; r < qRotasL2[s]; ++r)
		{
			if (values[r] >= 0.0001)
			{
				vector < int > vertices = ptrRotasL2[s][r]->getVertices();
				printf("r = new Rota();\n");
				for ( int i = 0; i < vertices.size(); ++i ) printf("r->inserirVerticeFim( %d );\n", vertices[i]);
				printf("r->setCusto(G);\nrGamma.push_back(r);\n	mCplex->insertColumn(G, r, %d);\n", s);
			}
		}
	}
}
