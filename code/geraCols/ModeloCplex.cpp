#include "ModeloCplex.h"
using namespace std;

int* ModeloCplex::qRotasL2;
int ModeloCplex::maxVeicL2;
int ModeloCplex::maxVeicL1;
int ModeloCplex::maxVeicSat;
int ModeloCplex::nVertices;
int ModeloCplex::nCustomers;
int ModeloCplex::nSatellites;
int ModeloCplex::quantRotasL1;
double ModeloCplex::limitePrimal;
double ModeloCplex::thresholdCapacityCuts;
double ModeloCplex::thresholdCombCuts;
double ModeloCplex::thresholdMultiStarCuts;
double ModeloCplex::thresholdCycleCuts;

vector<short int*>* ModeloCplex::A_qr;
vector<Rota*>* ModeloCplex::ptrRotasL2;

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

	//deleta as matrizes que armazenam os cortes capacityCuts
	for ( int c = 0; c < capacityCutsLHS.size(); ++c )
	{
		for ( int i = 0; i < nVertices; ++i ) delete [] capacityCutsLHS[c][i];
		delete [] capacityCutsLHS[c];
	}

	//deleta as matrizes que armazenam os cortes comb
	for ( int c = 0; c < combCutsLHS.size(); ++c )
	{
		for ( int i = 0; i < nVertices; ++i ) delete [] combCutsLHS[c][i];
		delete [] combCutsLHS[c];
	}

	//deleta as matrizes que armazenam os cortes multiStar
	for ( int m = 0; m < multiStarCutsLHS.size(); ++m )
	{
		for ( int i = 0; i < nVertices; ++i ) delete [] multiStarCutsLHS[m][i];
		delete [] multiStarCutsLHS[m];
	}

	//deleta os vetores que armazenam os cortes de ciclo
	for ( int s = 1; s <= nSatellites; ++s )
	{
		for ( int c = 0; c < cycleCutsLHS[s].size(); ++c )
		{
			delete [] cycleCutsLHS[s][c];
		}
	}

	CMGR_FreeMemCMgr(&oldCuts);
	CMGR_FreeMemCMgr(&newCuts);

	delete [] edgesHead;
	delete [] edgesTail;
	delete [] edgesX;
	delete [] demands;
}


ModeloCplex::ModeloCplex( Grafo* g, vector < Rota* > rLambda, vector < Rota* > rGamma ) : 
						env(), model(env), cplex(model), modelCGH(env), cplexCGH(modelCGH) {
	//abre o arquivo com os parametros de threshold dos cortes e os carrega na memoria
	fstream file("thresholds.txt");
	file >> thresholdCapacityCuts >> thresholdCombCuts >> thresholdMultiStarCuts >> thresholdCycleCuts;
	file.close();

	limitePrimal = MAIS_INFINITO;
	cplex.setOut(env.getNullStream());
	cplexCGH.setOut(env.getNullStream());

	maxVeicL1 = g->getMaxVeiculosL1();
	maxVeicL2 = g->getMaxVeiculosL2();
	maxVeicSat = g->getMaxVeiculosSatellites();

	nVertices = g->getNumVertices();
	nCustomers = g->getNumCustomers();
	nSatellites = g->getNumSatellites();
	quantRotasL1 = rLambda.size();

	setA_qr(rGamma);

	//Montagem do modelo Primal
	initVars( g );
	setObjectiveFunction( rLambda );
	setConstraints1( g );
	setConstraint2();
	setConstraints3();
	setConstraints4();
	setConstraints5( g );
	setConstraints6( g, rLambda );
	combCuts = IloRangeArray( env );
	capacityCuts = IloRangeArray( env );
	multiStarCuts = IloRangeArray( env );
	cycleCuts = IloArray < IloRangeArray > ( env, nSatellites+1 );
	routeCuts = IloArray < IloRangeArray > ( env, nSatellites+1 );
	constraintsArvoreVeic = IloArray < IloRangeArray > ( env, nSatellites+1 );
	cycleCutsLHS = new vector < char* >[nSatellites+1];
	routeCutsLHS = new vector < char* >[nSatellites+1];
	for ( int s = 1; s <= nSatellites; ++s )
	{
		cycleCuts[s] = IloRangeArray( env );
		routeCuts[s] = IloRangeArray( env );
		constraintsArvoreVeic[s] = IloRangeArray( env );
	}
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


double ModeloCplex::solveMaster( bool processaBase )
{
	cplex.solve();
	if ( processaBase )
	{
		IloNumArray redCosts(env);
		for ( int s = 1; s <= nSatellites; ++s )
		{
			cplex.getReducedCosts( redCosts, gamma[s] );
			for ( int r = 0; r < qRotasL2[s]; ++r )
			{
				if ( redCosts[r] < 0.001 ) //a variavel associada esta na base ou tem chances potenciais (custo reduzido 0)
				{
					A_qr[s][r][0] = 0;
				}
				else
				{
					--A_qr[s][r][0];
				}
			}
		}
		redCosts.end();
	}
	return cplex.getObjValue();
}


void ModeloCplex::initVars( Grafo* g ){
	lambda = IloNumVarArray( env, quantRotasL1, 0, maxVeicL1, ILOFLOAT );
	lambdaCGH = IloIntVarArray( env, quantRotasL1, 0, maxVeicL1 );
	for ( int r = 0; r < quantRotasL1; ++r )
	{
		char nome[10];
		sprintf(nome, "lbd_%d", r);
		lambda[r].setName(nome);
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
		delta[s] = IloNumVarArray( env, quantRotasL1, 0, limiteDelta, ILOFLOAT );
		deltaCGH[s] = IloIntVarArray( env, quantRotasL1, 0, limiteDelta );
		for (int r = 0; r < quantRotasL1; ++r)
		{
			char nome[20];
			sprintf(nome, "del_%01d_%01d", s, r);
			delta[s][r].setName(nome);
		}
	}

	//variaveis necessarias para executar o procedimento de separacao de desigualdades da biblioteca CVRPSEP
	CMGR_CreateCMgr(&oldCuts, 1);
	CMGR_CreateCMgr(&newCuts, 100);
	edgesHead = new int[nCustomers*nCustomers];
	edgesTail = new int[nCustomers*nCustomers];
	edgesX = new double[nCustomers*nCustomers];
	demands = new int[nCustomers+1];
	for ( int i = 1; i <= nCustomers; ++i) demands[i] = g->getCargaVertice(nSatellites + i);
	QMin = limiteDelta - ( ( maxVeicL2-1 ) * g->getCapacVeiculosL2() );

	//alocacao da memoria para armazenar a matriz do grafo usado para separar cycleCuts e os valores de y_i para o satellite s
	yCycleCuts = new double[nVertices];
	auxCycleCuts = new char[3*nVertices];
	grafoCycleCuts = new double*[nVertices];
	for ( int i = 0; i < nVertices; ++i) grafoCycleCuts[i] = new double[nVertices];
}


void ModeloCplex::setObjectiveFunction( vector < Rota* > rLambda ){
	IloExpr obj(env), objCGH(env);

	for ( int r = 0; r < rLambda.size(); ++r )
	{
		obj += rLambda[r]->getCusto() * lambda[r];
		objCGH += rLambda[r]->getCusto() * lambdaCGH[r];
		
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


void ModeloCplex::setConstraints1( Grafo* g ){
	double cargaTotal = 0;
	IloExpr exp1(env), exp1CGH(env);
	for ( int i = 1; i <= nCustomers; ++i ) cargaTotal += g->getCargaVertice(nSatellites + i);

	for ( int r = 0; r < quantRotasL1; ++r )
	{
		exp1 += lambda[r];
		exp1CGH += lambdaCGH[r];
	}
	model.add(exp1 <= maxVeicL1);
	modelCGH.add(exp1CGH <= maxVeicL1);

	//inclui uma restricao para apertar o modelo
	model.add ( exp1 >= ceil ( cargaTotal / g->getCapacVeiculosL1() ) );
	modelCGH.add ( exp1CGH >= ceil ( cargaTotal / g->getCapacVeiculosL1() ) );

	exp1.end(); exp1CGH.end();
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
		for ( int r = 0; r < quantRotasL1; ++r )
		{
			exp5 -= delta[s][r];
			exp5CGH -= deltaCGH[s][r];
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
	bool visita;
	vector < int > v;
	constraints6 = IloRangeArray(env, rLambda.size());
	constraints6CGH = IloRangeArray(env, rLambda.size());

	IloExpr exp7(env), exp7CGH(env);
	for ( int r = 0; r < rLambda.size(); ++r )
	{
		IloExpr exp6(env), exp6CGH(env);
		v = rLambda[r]->getVertices();
		for ( int s = 1; s <= nSatellites; ++s )
		{
			visita = false;
			for ( int i = 1; i < v.size(); ++i )
			{
				if ( v[i] == s )
				{
					exp6 += delta[s][r];
					exp6CGH += deltaCGH[s][r];
					visita = true;
					break;
				}
			}
			if ( !visita )
			{
				exp7 += delta[s][r];
				exp7CGH += deltaCGH[s][r];
			}
		}

		constraints6[r] = ( ( exp6 - g->getCapacVeiculosL1()*lambda[r] ) <= 0);
		char nome[20];
		sprintf(nome, "rest6_r%d", r);
		constraints6[r].setName(nome);
		model.add(constraints6[r]);
		
		constraints6CGH[r] = ( ( exp6CGH - g->getCapacVeiculosL1()*lambdaCGH[r] ) <= 0);
		modelCGH.add(constraints6CGH[r]);

		exp6.end();
		exp6CGH.end();
	}

	constraint7 = ( exp7 == 0 );
	constraintCGH7 = ( exp7CGH == 0 );
	model.add( constraint7 );
	modelCGH.add( constraintCGH7 );
	exp7.end(); exp7CGH.end();
}


void ModeloCplex::updateDualCosts( Grafo* g, int s ){
	double custoDual;
	int totalCycleCuts;
	IloNumArray muDuals(env), chi(env), zeta(env);
	int numCycleCuts = cycleCutsLHS[s].size();
	int numRouteCuts = routeCutsLHS[s].size();

	cplex.getDuals(muDuals, constraints4);
	double piDual = cplex.getDual(constraints5[s-1]);
	if ( numCycleCuts > 0 ) cplex.getDuals(chi, cycleCuts[s]);
	if ( numRouteCuts > 0 ) cplex.getDuals(zeta, routeCuts[s]);

	for (int i = nSatellites+1; i < nVertices; ++i)
	{
		totalCycleCuts = 0;
		custoDual = -muDuals[i-nSatellites-1] - g->getCargaVertice(i)*piDual;

		for ( int x = 0; x < numCycleCuts; ++x )
		{
			for ( int y = 1; y <= cycleCutsLHS[s][x][0]; ++y )
			{
				for ( int z = 1; z <= cycleCutsLHS[s][x][0]; ++z )
				{
					if ( ( y != z ) && ( cycleCutsLHS[s][x][z] == i ) ) custoDual += chi[totalCycleCuts];
				}
				++totalCycleCuts;
			}
		}

		for ( int j = 0; j < numRouteCuts; ++j ) //custos duais nos vertices vindos dos routeFeasibilityCuts
		{
			for ( int k = 1; k <= routeCutsLHS[s][j][0]; ++k )
			{
				if ( routeCutsLHS[s][j][k] == i ) custoDual += zeta[j];
			}
		}
		g->setCustoVerticeDual(i, custoDual);
	}

	g->setCustoArestasDual(s);
	muDuals.end();

	//processa os custos duais das combs e multistar juntas (caso existam)
	int totalCombCuts = combCutsLHS.size();
	int totalCapacityCuts = capacityCutsLHS.size();
	int totalMultiStarCuts = multiStarCutsLHS.size();
	if ( ( totalCombCuts > 0 ) || ( totalCapacityCuts > 0 ) || ( totalMultiStarCuts > 0 ) )
	{
		IloNumArray psiDuals(env), roDuals(env), omegaDuals(env);
		if ( totalCombCuts > 0 ) cplex.getDuals(psiDuals, combCuts);
		if ( totalCapacityCuts > 0 ) cplex.getDuals(roDuals, capacityCuts);
		if ( totalMultiStarCuts > 0 ) cplex.getDuals(omegaDuals, multiStarCuts);

		for ( int i = 0; i < nVertices; ++i )
		{
			for ( int j = i+1; j < nVertices; ++j )
			{
				custoDual = 0;
				for ( int c = 0; c < totalCombCuts; ++c )
				{
					if ( combCutsLHS[c][i][j] > 0 ) custoDual += combCutsLHS[c][i][j]*psiDuals[c];
				}
				for ( int c = 0; c < totalCapacityCuts; ++c )
				{
					if ( capacityCutsLHS[c][i][j] != 0 ) custoDual += capacityCutsLHS[c][i][j]*roDuals[c];
				}
				for ( int m = 0; m < totalMultiStarCuts; ++m )
				{
					if ( multiStarCutsLHS[m][i][j] != 0 ) custoDual += multiStarCutsLHS[m][i][j]*omegaDuals[m];
				}
				if ( ( custoDual > 0.00001 ) || ( custoDual < -0.00001 ) )
				{
					if ( i == 0 )
					{
						g->setCustoArestaDual(0, j-nSatellites, -custoDual);
						g->setCustoArestaDual(j-nSatellites, 0, -custoDual);
					}
					else
					{
						g->setCustoArestaDual(i-nSatellites, j-nSatellites, -custoDual);
						g->setCustoArestaDual(j-nSatellites, i-nSatellites, -custoDual);
					}
				}
			}
			if ( i == 0 ) i = nSatellites; //para 'pular' as linhas dos satellites, que tem todos os valores 0
		}
		psiDuals.end(); roDuals.end(); omegaDuals.end();
	}

	//processa os custos duais das cycleCuts inequalities
	totalCycleCuts = 0;
	for ( int i = 0; i < numCycleCuts; ++i )
	{
		for ( int k = 1; k <= cycleCutsLHS[s][i][0]; ++k )
		{
			for ( int j = 1; j < cycleCutsLHS[s][i][0]; ++j )
			{
				g->setCustoArestaDual(cycleCutsLHS[s][i][j]-nSatellites, cycleCutsLHS[s][i][j+1]-nSatellites, -chi[totalCycleCuts]);
				g->setCustoArestaDual(cycleCutsLHS[s][i][j+1]-nSatellites, cycleCutsLHS[s][i][j]-nSatellites, -chi[totalCycleCuts]);
			}

			for ( int j = 3; j <= cycleCutsLHS[s][i][0]; ++j )
			{
				for ( int k = 1; k <= (j-2); ++k )
				{
					g->setCustoArestaDual(cycleCutsLHS[s][i][j]-nSatellites, cycleCutsLHS[s][i][k]-nSatellites, -chi[totalCycleCuts]);
					g->setCustoArestaDual(cycleCutsLHS[s][i][k]-nSatellites, cycleCutsLHS[s][i][j]-nSatellites, -chi[totalCycleCuts]);
				}
			}
			++totalCycleCuts;
		}
	}

	//processa os custos duais das routeFeasibilityCuts inequalities
	for ( int i = 0; i < numRouteCuts; ++i )
	{
		for ( int j = 1; j <= routeCutsLHS[s][i][0]; ++j )
		{
			g->setCustoArestaDual(0, routeCutsLHS[s][i][j]-nSatellites, -zeta[i]);
			g->setCustoArestaDual(routeCutsLHS[s][i][j]-nSatellites, 0, -zeta[i]);

			for ( int k = j+1; k <= routeCutsLHS[s][i][0]; ++k )
			{
				g->setCustoArestaDual(routeCutsLHS[s][i][j]-nSatellites, routeCutsLHS[s][i][k]-nSatellites, -zeta[i]);
				g->setCustoArestaDual(routeCutsLHS[s][i][k]-nSatellites, routeCutsLHS[s][i][j]-nSatellites, -zeta[i]);
			}
		}
	}
	chi.end(); zeta.end();
}


double ModeloCplex::getValorSubtrair( Grafo* g, int s ){
	IloNumArray values(env);
	double somaConstraintsArvoreVeic = 0;
	cplex.getDuals(values, constraintsArvoreVeic[s]);
	for ( int i = 0; i < constraintsArvoreVeic[s].getSize(); ++i ) somaConstraintsArvoreVeic += values[i];
	return ( cplex.getDual( constraint2 ) + cplex.getDual( constraints3[ s-1 ] ) + somaConstraintsArvoreVeic );
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

			if ( A_qr[s][qRotasL2[s]][nSatellites + i + 1] > 1 ) r->ciclo = true;
		}
	}
	col += constraints5[s-1](coeficiente5);
	colCGH += constraints5CGH[s-1](coeficiente5);

	//insere a coluna nas restricoes dos cortes de capacidade (caso existam)
	float aux, indexCol;
	int numCapacityCuts = capacityCutsLHS.size();
	for ( int i = 0; i < numCapacityCuts; ++i )
	{
		indexCol = 0;
		for ( int v = 1; v < numVertAtualizar; ++v )
		{
			if ( capacityCutsLHS[i][vertRota[v]][vertRota[v+1]] != 0 ) indexCol += capacityCutsLHS[i][vertRota[v]][vertRota[v+1]];
			if ( capacityCutsLHS[i][vertRota[v+1]][vertRota[v]] != 0 ) indexCol += capacityCutsLHS[i][vertRota[v+1]][vertRota[v]];
		}
		if ( capacityCutsLHS[i][0][vertRota[1]] != 0 ) indexCol += capacityCutsLHS[i][0][vertRota[1]];
		if ( capacityCutsLHS[i][0][vertRota[numVertAtualizar]] != 0 ) indexCol += capacityCutsLHS[i][0][vertRota[numVertAtualizar]];

		if ( ( indexCol > 0.00001 ) || ( indexCol < -0.00001 ) ) col += capacityCuts[i](indexCol);
	}

	//insere a coluna nas restricoes dos cortes comb (caso existam)
	int numCombCuts = combCutsLHS.size();
	for ( int i = 0; i < numCombCuts; ++i )
	{
		indexCol = 0;
		for ( int v = 1; v < numVertAtualizar; ++v )
		{
			if ( combCutsLHS[i][vertRota[v]][vertRota[v+1]] > 0 ) indexCol += combCutsLHS[i][vertRota[v]][vertRota[v+1]];
			if ( combCutsLHS[i][vertRota[v+1]][vertRota[v]] > 0 ) indexCol += combCutsLHS[i][vertRota[v+1]][vertRota[v]];
		}
		if ( combCutsLHS[i][0][vertRota[1]] > 0 ) indexCol += combCutsLHS[i][0][vertRota[1]];
		if ( combCutsLHS[i][0][vertRota[numVertAtualizar]] > 0 ) indexCol += combCutsLHS[i][0][vertRota[numVertAtualizar]];

		if ( indexCol > 0.00001 ) col += combCuts[i](indexCol);
	}

	//insere a coluna nas restricoes dos cortes multistar (caso existam)
	int numMultiStarCuts = multiStarCutsLHS.size();
	for ( int i = 0; i < numMultiStarCuts; ++i )
	{
		indexCol = 0;
		for ( int v = 1; v < numVertAtualizar; ++v )
		{
			if ( multiStarCutsLHS[i][vertRota[v]][vertRota[v+1]] != 0 ) indexCol += multiStarCutsLHS[i][vertRota[v]][vertRota[v+1]];
			if ( multiStarCutsLHS[i][vertRota[v+1]][vertRota[v]] != 0 ) indexCol += multiStarCutsLHS[i][vertRota[v+1]][vertRota[v]];
		}
		if ( multiStarCutsLHS[i][0][vertRota[1]] != 0 ) indexCol += multiStarCutsLHS[i][0][vertRota[1]];
		if ( multiStarCutsLHS[i][0][vertRota[numVertAtualizar]] != 0 ) indexCol += multiStarCutsLHS[i][0][vertRota[numVertAtualizar]];

		if ( ( indexCol > 0.00001 ) || ( indexCol < -0.00001 ) ) col += multiStarCuts[i](indexCol);
	}

	//insere a coluna nas restricoes dos cortes outfork (caso existam)
	int totalCycleCuts = 0;
	int numCycleCutsLHS = cycleCutsLHS[s].size();
	for ( int i = 0; i < numCycleCutsLHS; ++i )
	{
		indexCol = 0;
		for ( int v = 1; v < numVertAtualizar; ++v )
		{
			for ( int y = 1; y < cycleCutsLHS[s][i][0]; ++y )
			{
				if ( ( ( vertRota[v] == cycleCutsLHS[s][i][y] ) && ( vertRota[v+1] == cycleCutsLHS[s][i][y+1] ) ) ||
					 ( ( vertRota[v] == cycleCutsLHS[s][i][y+1] ) && ( vertRota[v+1] == cycleCutsLHS[s][i][y] ) ) ) ++indexCol;
			}
			for ( int y = 3; y <= cycleCutsLHS[s][i][0]; ++y )
			{
				for ( int z = 1; z <= (y-2); ++z )
				{
					if ( ( ( vertRota[v] == cycleCutsLHS[s][i][y] ) && ( vertRota[v+1] == cycleCutsLHS[s][i][z] ) ) ||
						( ( vertRota[v] == cycleCutsLHS[s][i][z] ) && ( vertRota[v+1] == cycleCutsLHS[s][i][y] ) ) ) ++indexCol;
				}
			}
		}

		for ( int x = 1; x <= cycleCutsLHS[s][i][0]; ++x )
		{
			aux = 0;
			for ( int y = 1; y <= cycleCutsLHS[s][i][0]; ++y )
			{
				if ( ( x != y ) && ( A_qr[s][qRotasL2[s]][cycleCutsLHS[s][i][y]] > 0 ) ) aux += A_qr[s][qRotasL2[s]][cycleCutsLHS[s][i][y]];
			}

			indexCol -= aux;
			if ( ( indexCol > 0.00001 ) || ( indexCol < -0.00001 ) ) col += cycleCuts[s][totalCycleCuts](indexCol);
			++totalCycleCuts;
		}
	}

	//insere a coluna nas restricoes dos routeFeasibilityCuts (caso existam)
	int numRouteCuts = routeCutsLHS[s].size();
	for ( int i = 0; i < numRouteCuts; ++i )
	{
		indexCol = 0;
		for ( int v = 1; v < numVertAtualizar; ++v )
		{
			for ( int y = 1; y < routeCutsLHS[s][i][0]; ++y )
			{
				for ( int z = y+1; z <= routeCutsLHS[s][i][0]; ++z )
				{
					if ( ( ( vertRota[v] == routeCutsLHS[s][i][y] ) && ( vertRota[v+1] == routeCutsLHS[s][i][z] ) ) ||
						( ( vertRota[v] == routeCutsLHS[s][i][z] ) && ( vertRota[v+1] == routeCutsLHS[s][i][y] ) ) ) ++indexCol;
				}
			}
		}
		for ( int y = 1; y <= routeCutsLHS[s][i][0]; ++y )
		{
			if ( vertRota[1] == routeCutsLHS[s][i][y] ) ++indexCol;
			if ( vertRota[numVertAtualizar] == routeCutsLHS[s][i][y] ) ++indexCol;
			if ( A_qr[s][qRotasL2[s]][routeCutsLHS[s][i][y]] > 0 ) indexCol -= A_qr[s][qRotasL2[s]][routeCutsLHS[s][i][y]];
		}

		if ( ( indexCol > 0.00001 ) || ( indexCol < -0.00001 ) ) col += routeCuts[s][i](indexCol);
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
	char solIntegerFeasible;
	bool inseriuCortes = false;
	double maxViolate, RHS, violacao;

	newCuts->Size = 0;
	buildGraphSolution(numEdges, edgesHead, edgesTail, edgesX, qRotas_no, ptrRotas_no);
	CAPSEP_SeparateCapCuts(nCustomers, demands, g->getCapacVeiculosL2(), numEdges, edgesTail, edgesHead, edgesX, oldCuts, 100, 0.0001, &solIntegerFeasible, &maxViolate, newCuts);

	if ( solIntegerFeasible ) return false;
	if ( newCuts->Size == 0 ) return false;

	//procura pela desigualdade com maior violacao entre as retornadas
	int i = mostViolatedCapacityCut( qRotas_no, ptrRotas_no, violacao, 1 );
	if ( ( i < 0 ) || ( violacao < ( thresholdCapacityCuts+0.0001 ) ) ) return false; //sem cortes com a violacao minima

	vector<int> vertices;
	short int* sBarra = new short int[1+nCustomers/2];
	int head, tail, numVertices, indexCCC = capacityCutsLHS.size();
	while ( ( i >= 0 ) && ( violacao > thresholdCapacityCuts ) )
	{
		//aloca a matriz para armazenar o corte e que sera armazenada por capacityCutsLHS
		float** matrizCut = new float*[nVertices];
		for ( int j = 0; j < nVertices; ++j )
		{
			matrizCut[j] = new float[nVertices];
			memset(matrizCut[j], 0, nVertices*sizeof(float));
		}

		if ( newCuts->CPL[i]->IntListSize <= ( nCustomers / 2 ) ) 
		{
			for ( int j = 1; j < newCuts->CPL[i]->IntListSize; ++j )
			{
				head = newCuts->CPL[i]->IntList[j]+nSatellites;
				for ( int k = j+1; k <= newCuts->CPL[i]->IntListSize; ++k )
				{
					tail = newCuts->CPL[i]->IntList[k]+nSatellites;
					if ( head < tail ) ++matrizCut[head][tail];
					else ++matrizCut[tail][head];
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

			for ( int j = 0; j < tamSBarra; ++j ) // + 0.5*x({0}:sBarra)
			{
				matrizCut[0][sBarra[j]+nSatellites] += 0.5;
			}

			for ( int j = 1; j <= newCuts->CPL[i]->IntListSize; ++j ) // - 0.5*x({0}:S)
			{
				matrizCut[0][newCuts->CPL[i]->IntList[j]+nSatellites] -= 0.5;
			}

			for ( int j = 0; j < tamSBarra-1; ++j ) // x(sBarra:sBarra)
			{
				head = ( sBarra[j] + nSatellites );
				for ( int k = j+1; k < tamSBarra; ++k )
				{
					tail = ( sBarra[k] + nSatellites );
					if ( head < tail ) ++matrizCut[head][tail];
					else ++matrizCut[tail][head];
				}
			}
			RHS = ( tamSBarra - ( newCuts->CPL[i]->IntListSize - newCuts->CPL[i]->RHS ) );
		}

		//Monta a restricao a ser incluida no modelo
		IloExpr expCCC = IloExpr(env);
		for ( int s = 1; s <= nSatellites; ++s )
		{
			for ( int r = 0; r < qRotasL2[s]; ++r )
			{
				vertices = ptrRotasL2[s][r]->getVertices();
				numVertices = vertices.size()-2;
				for ( int e = 1; e < numVertices; ++e )
				{
					if ( matrizCut[vertices[e]][vertices[e+1]] != 0 )
					{
						expCCC += matrizCut[vertices[e]][vertices[e+1]]*gamma[s][r];
					}
					if ( matrizCut[vertices[e+1]][vertices[e]] != 0 )
					{
						expCCC += matrizCut[vertices[e+1]][vertices[e]]*gamma[s][r];
					}
				}
				if ( matrizCut[0][vertices[1]] != 0 )
				{
					expCCC += matrizCut[0][vertices[1]]*gamma[s][r];
				}
				if ( matrizCut[0][vertices[numVertices]] != 0 )
				{
					expCCC += matrizCut[0][vertices[numVertices]]*gamma[s][r];
				}
			}

			for ( int r = 0; r < qRotas_no[s]; ++r )
			{
				vertices = ptrRotas_no[s][r]->getVertices();
				numVertices = vertices.size()-2;
				for ( int e = 1; e < numVertices; ++e )
				{
					if ( matrizCut[vertices[e]][vertices[e+1]] != 0 )
					{
						expCCC += matrizCut[vertices[e]][vertices[e+1]]*gamma_no[s][r];
					}
					if ( matrizCut[vertices[e+1]][vertices[e]] != 0 )
					{
						expCCC += matrizCut[vertices[e+1]][vertices[e]]*gamma_no[s][r];
					}
				}
				if ( matrizCut[0][vertices[1]] != 0 )
				{
					expCCC += matrizCut[0][vertices[1]]*gamma_no[s][r];
				}
				if ( matrizCut[0][vertices[numVertices]] != 0 )
				{
					expCCC += matrizCut[0][vertices[numVertices]]*gamma_no[s][r];
				}
			}
		}

		//inclui a restricao no modelo e armazena a matriz do corte no vector
		capacityCutsLHS.push_back( matrizCut );
		capacityCuts.add( expCCC <= RHS );
		model.add(capacityCuts[indexCCC]);
		expCCC.end();
		++indexCCC;

		//verifica o proximo corte a ser inserido (caso satisfaca a violacao minima)
		i = mostViolatedCapacityCut( qRotas_no, ptrRotas_no, violacao, 1 );
		inseriuCortes = true;
	}

	delete [] sBarra;
	return inseriuCortes;
}


bool ModeloCplex::addCombCuts( Grafo* g, int* qRotas_no, vector < Rota* >* ptrRotas_no )
{
	vector<int> vertices;
	bool inseriuCortes = false;
	double maxViolate, violacao;
	int numEdges, nTeeth, minIdx, maxIdx, tmp, tam;
	char **matrizEdgesCut;
	char **inTooth = new char*[nCustomers];
	for ( int i = 0; i < nCustomers; ++i ) inTooth[i] = new char[nCustomers+2]; //armazena se os vertices estao ou nao {0,1} no handle/teeth

	newCuts->Size = 0;
	buildGraphSolution(numEdges, edgesHead, edgesTail, edgesX, qRotas_no, ptrRotas_no);
	COMBSEP_SeparateCombs(nCustomers, demands, g->getCapacVeiculosL2(), QMin, numEdges, edgesTail, edgesHead, edgesX, 100, &maxViolate, newCuts);
	if ( newCuts->Size == 0 ) return false;

	//procura pela desigualdade com maior violacao entre as retornadas
	int indexCombCuts = combCutsLHS.size();
	int i = mostViolatedCapacityCut( qRotas_no, ptrRotas_no, violacao, 2 );
	if ( ( i < 0 ) || ( violacao < ( thresholdCombCuts+0.0001 ) ) ) return false; //sem cortes com a violacao minima
	while ( ( i >= 0 ) && ( violacao > thresholdCombCuts ) )
	{
		//aloca memoria para a matriz q armazena o corte e inicializa com '0' as entradas usadas para monta-lo
		matrizEdgesCut = new char*[nVertices];
		for ( int j = 0; j < nVertices; ++j )
		{
			matrizEdgesCut[j] = new char[nVertices];
			memset(matrizEdgesCut[j], 0, nVertices*sizeof(char));
			if ( j < nCustomers ) memset(inTooth[j], 0, (nCustomers+2)*sizeof(char));
		}

		//marca na matrix auxiliar com 1 os vertices que estao no Handle
		for ( int h = 1; h <= newCuts->CPL[i]->IntListSize; h++ )
		{
			inTooth[0][newCuts->CPL[i]->IntList[h]] = 1;
		}

		//marca na matriz auxiliar com 1 os vertices q estao em cada Tooth
		nTeeth = newCuts->CPL[i]->Key;
		for ( int t = 1; t <= nTeeth; t++ )
		{
			minIdx = newCuts->CPL[i]->ExtList[t];
			if (t == nTeeth) maxIdx = newCuts->CPL[i]->ExtListSize;
			else maxIdx = newCuts->CPL[i]->ExtList[t+1] - 1;

			for ( int k = minIdx; k <= maxIdx; k++ )
			{
				inTooth[t][newCuts->CPL[i]->ExtList[k]] = 1;
			}
		}

		//usa a matriz auxiliar para atribuir agora o valor 1 para aquelas arestas q estao no corte associado ao conjunto H
		for ( int h = 1; h <= newCuts->CPL[i]->IntListSize; h++ )
		{
			tmp = newCuts->CPL[i]->IntList[h];
			for ( int j = 1; j <= nCustomers; ++j )
			{
				if ( inTooth[0][j] == 0 )
				{
					if ( tmp > nCustomers ) ++matrizEdgesCut[0][j+nSatellites]; //caso o vertice de H seja o deposito
					else if ( tmp < j ) ++matrizEdgesCut[tmp+nSatellites][j+nSatellites];
					else ++matrizEdgesCut[j+nSatellites][tmp+nSatellites];
				}
			}
			if ( inTooth[0][nCustomers+1] == 0 ) ++matrizEdgesCut[0][tmp+nSatellites]; //considera o deposito fora de H
		}

		//atribui na matriz valor 1 para aquelas arestas q estao no corte associado aos conjuntos T
		for ( int t = 1; t <= nTeeth; t++ )
		{
			minIdx = newCuts->CPL[i]->ExtList[t];
			if (t == nTeeth) maxIdx = newCuts->CPL[i]->ExtListSize;
			else maxIdx = newCuts->CPL[i]->ExtList[t+1] - 1;

			for ( int j = minIdx; j <= maxIdx; j++ )
			{
				tmp = newCuts->CPL[i]->ExtList[j];
				for ( int k = 1; k <= nCustomers; ++k )
				{
					if ( inTooth[t][k] == 0 )
					{
						if ( tmp > nCustomers ) ++matrizEdgesCut[0][k+nSatellites]; //caso o vertice de T seja o deposito
						else if ( tmp < k ) ++matrizEdgesCut[tmp+nSatellites][k+nSatellites];
						else ++matrizEdgesCut[k+nSatellites][tmp+nSatellites];
					}
				}
				if ( inTooth[t][nCustomers+1] == 0 ) ++matrizEdgesCut[0][tmp+nSatellites]; //considera o deposito fora de T
			}
		}

		//apos a montagem das matrizes, monta-se a restricao associada ao corte para ser incluida no modelo
		IloExpr expComb(env);
		for ( int s = 1; s <= nSatellites; ++s )
		{
			for ( int r = 0; r < qRotasL2[s]; ++r )
			{
				vertices = ptrRotasL2[s][r]->getVertices();
				tam = vertices.size()-2;
				for ( int v = 1; v < tam; ++v )
				{
					if ( matrizEdgesCut[vertices[v]][vertices[v+1]] > 0 ) expComb += matrizEdgesCut[vertices[v]][vertices[v+1]]*gamma[s][r];
					if ( matrizEdgesCut[vertices[v+1]][vertices[v]] > 0 ) expComb += matrizEdgesCut[vertices[v+1]][vertices[v]]*gamma[s][r];
				}
				if ( matrizEdgesCut[0][vertices[1]] > 0 ) expComb += matrizEdgesCut[0][vertices[1]]*gamma[s][r];
				if ( matrizEdgesCut[0][vertices[tam]] > 0 ) expComb += matrizEdgesCut[0][vertices[tam]]*gamma[s][r];
			}

			for ( int r = 0; r < qRotas_no[s]; ++r )
			{
				vertices = ptrRotas_no[s][r]->getVertices();
				tam = vertices.size()-2;
				for ( int v = 1; v < tam; ++v )
				{
					if ( matrizEdgesCut[vertices[v]][vertices[v+1]] > 0 ) expComb += matrizEdgesCut[vertices[v]][vertices[v+1]]*gamma_no[s][r];
					if ( matrizEdgesCut[vertices[v+1]][vertices[v]] > 0 ) expComb += matrizEdgesCut[vertices[v+1]][vertices[v]]*gamma_no[s][r];
				}
				if ( matrizEdgesCut[0][vertices[1]] > 0 ) expComb += matrizEdgesCut[0][vertices[1]]*gamma_no[s][r];
				if ( matrizEdgesCut[0][vertices[tam]] > 0 ) expComb += matrizEdgesCut[0][vertices[tam]]*gamma_no[s][r];
			}
		}

		combCutsLHS.push_back(matrizEdgesCut);
		combCuts.add(expComb >= newCuts->CPL[i]->RHS);
		model.add(combCuts[indexCombCuts]);
		++indexCombCuts;
		expComb.end();

		//verifica o proximo corte a ser inserido (caso satisfaca a violacao minima)
		i = mostViolatedCapacityCut( qRotas_no, ptrRotas_no, violacao, 2 );
		inseriuCortes = true;
	}

	for ( int i = 0; i < nCustomers; ++i ) delete [] inTooth[i]; delete inTooth;
	return inseriuCortes;
}


bool ModeloCplex::addMultiStarCuts( Grafo* g, int* qRotas_no, vector < Rota* >* ptrRotas_no )
{
	vector<int> vertices;
	bool inseriuCortes = false;
	double maxViolate, violacao;
	int numEdges, numVertices, c;
	int indexMSC = multiStarCutsLHS.size();

	newCuts->Size = 0;
	buildGraphSolution(numEdges, edgesHead, edgesTail, edgesX, qRotas_no, ptrRotas_no);
	MSTARSEP_SeparateMultiStarCuts(nCustomers, demands, g->getCapacVeiculosL2(), numEdges, edgesTail, edgesHead, edgesX, oldCuts, 100, &maxViolate, newCuts);
	if ( newCuts->Size == 0 ) return false;

	//procura pela desigualdade com maior violacao entre as retornadas
	int i = mostViolatedCapacityCut( qRotas_no, ptrRotas_no, violacao, 3 );
	if ( ( i < 0 ) || ( violacao < ( thresholdMultiStarCuts+0.0001 ) ) ) return false; //sem cortes com a violacao minima
	while ( ( i >= 0 ) && ( violacao > thresholdCombCuts ) )
	{
		//aloca memoria para a matriz q armazena o corte e inicializa com '0' as entradas usadas para monta-lo
		short int** matrizEdgesCut = new short int*[nVertices];
		for ( int j = 0; j < nVertices; ++j )
		{
			matrizEdgesCut[j] = new short int[nVertices];
			memset(matrizEdgesCut[j], 0, nVertices*sizeof(short int));
		}

		//define as arestas de x(\delta(N))
		for ( int v = 1; v <= ( nCustomers+1 ); ++v )
		{
			for ( c = 1; c <= newCuts->CPL[i]->IntListSize; ++c )
			{
				if ( newCuts->CPL[i]->IntList[c] == v ) break;
			}
			if ( c > newCuts->CPL[i]->IntListSize )
			{
				for ( c = 1; c <= newCuts->CPL[i]->IntListSize; ++c )
				{
					if ( v > nCustomers ) matrizEdgesCut[0][newCuts->CPL[i]->IntList[c]+nSatellites] += newCuts->CPL[i]->B;
					else if ( newCuts->CPL[i]->IntList[c] > nCustomers ) matrizEdgesCut[0][v+nSatellites] += newCuts->CPL[i]->B;
					else if ( v < newCuts->CPL[i]->IntList[c] ) matrizEdgesCut[v+nSatellites][newCuts->CPL[i]->IntList[c]+nSatellites] += newCuts->CPL[i]->B;
					else matrizEdgesCut[newCuts->CPL[i]->IntList[c]+nSatellites][v+nSatellites] += newCuts->CPL[i]->B;
				}
			}
		}

		//define as arestas de x(C:T)
		for ( c = 1; c <= newCuts->CPL[i]->CListSize; ++c )
		{
			for ( int t = 1; t <= newCuts->CPL[i]->ExtListSize; ++t )
			{
				if ( newCuts->CPL[i]->CList[c] > nCustomers ) matrizEdgesCut[0][newCuts->CPL[i]->ExtList[t]+nSatellites] -= newCuts->CPL[i]->A;
				else if ( newCuts->CPL[i]->ExtList[t] > nCustomers ) matrizEdgesCut[0][newCuts->CPL[i]->CList[c]+nSatellites] -= newCuts->CPL[i]->A;
				else if ( newCuts->CPL[i]->CList[c] < newCuts->CPL[i]->ExtList[t] ) matrizEdgesCut[newCuts->CPL[i]->CList[c]+nSatellites][newCuts->CPL[i]->ExtList[t]+nSatellites] -= newCuts->CPL[i]->A;
				else matrizEdgesCut[newCuts->CPL[i]->ExtList[t]+nSatellites][newCuts->CPL[i]->CList[c]+nSatellites] -= newCuts->CPL[i]->A;
			}
		}

		//para cada aresta no corte, verifica quantas vezes cada rota passa e inclui no corte
		IloExpr expMSC = IloExpr(env);
		for ( int s = 1; s <= nSatellites; ++s )
		{
			for ( int r = 0; r < qRotasL2[s]; ++r )
			{
				vertices = ptrRotasL2[s][r]->getVertices();
				numVertices = vertices.size()-2;
				for ( int v = 1; v < numVertices; ++v )
				{
					if ( matrizEdgesCut[vertices[v]][vertices[v+1]] != 0 ) expMSC += matrizEdgesCut[vertices[v]][vertices[v+1]]*gamma[s][r];
					if ( matrizEdgesCut[vertices[v+1]][vertices[v]] != 0 ) expMSC += matrizEdgesCut[vertices[v+1]][vertices[v]]*gamma[s][r];
				}
				if ( matrizEdgesCut[0][vertices[1]] != 0 ) expMSC += matrizEdgesCut[0][vertices[1]]*gamma[s][r];
				if ( matrizEdgesCut[0][vertices[numVertices]] != 0 ) expMSC += matrizEdgesCut[0][vertices[numVertices]]*gamma[s][r];
			}

			for ( int r = 0; r < qRotas_no[s]; ++r )
			{
				vertices = ptrRotas_no[s][r]->getVertices();
				numVertices = vertices.size()-2;
				for ( int v = 1; v < numVertices; ++v )
				{
					if ( matrizEdgesCut[vertices[v]][vertices[v+1]] != 0 ) expMSC += matrizEdgesCut[vertices[v]][vertices[v+1]]*gamma_no[s][r];
					if ( matrizEdgesCut[vertices[v+1]][vertices[v]] != 0 ) expMSC += matrizEdgesCut[vertices[v+1]][vertices[v]]*gamma_no[s][r];
				}
				if ( matrizEdgesCut[0][vertices[1]] != 0 ) expMSC += matrizEdgesCut[0][vertices[1]]*gamma_no[s][r];
				if ( matrizEdgesCut[0][vertices[numVertices]] != 0 ) expMSC += matrizEdgesCut[0][vertices[numVertices]]*gamma_no[s][r];
			}
		}

		multiStarCutsLHS.push_back(matrizEdgesCut);
		multiStarCuts.add(expMSC >= newCuts->CPL[i]->L);
		model.add(multiStarCuts[indexMSC]);
		expMSC.end();
		++indexMSC;

		//verifica o proximo corte a ser inserido (caso satisfaca a violacao minima)
		i = mostViolatedCapacityCut( qRotas_no, ptrRotas_no, violacao, 3 );
		inseriuCortes = true;
	}
	return inseriuCortes;
}


bool ModeloCplex::addCycleCuts( Grafo* g, int* qRotas_no, vector<short int*>* A_ir_no, vector < Rota* >* ptrRotas_no )
{
	bool existeCorte = false;
	double somaLHS, somaY, maiorY;
	vector < int > verticesRota, verticesRotaAux;
	int numVerticesRota, numVerticesRotaAux, indexCycleCuts;

	//pega os valores das variaveis de decisao de todos os satellites, para nao resolver o PL sucessivamente
	IloArray<IloNumArray> values(env, nSatellites+1);
	IloArray<IloNumArray> values_no(env, nSatellites+1);
	for ( int s = 1; s <= nSatellites; ++s )
	{
		values[s] = IloNumArray(env);
		values_no[s] = IloNumArray(env);
		cplex.getValues( values[s], gamma[s] );
		if ( qRotas_no[s] > 0 ) cplex.getValues( values_no[s], gamma_no[s] );
	}

	for ( int s = 1; s <= nSatellites; ++s )
	{
		indexCycleCuts = cycleCuts[s].getSize();
		buildGraphCycleCuts( s, values, values_no, qRotas_no, ptrRotas_no );

		for ( int r = 0; r < qRotasL2[s]; ++r )
		{
			if ( ( values[s][r] > 0.00001 ) && ( ptrRotasL2[s][r]->ciclo ) ) //deve-se incluir no modelo o corte para esta rota
			{
				ptrRotasL2[s][r]->ciclo = false;
				verticesRota = ptrRotasL2[s][r]->getVertices();
				numVerticesRota = verticesRota.size()-2;

				//verifica cada ciclo (que nao contenha subciclos) identificado na rota
				for ( int i = 1; i < numVerticesRota; ++i )
				{
					for ( int j = i+1; j <= numVerticesRota; ++j )
					{
						if ( verticesRota[i] == verticesRota[j] )
						{
							int cont = 0;
							bool corteInvalido = false;
							for ( int n = i; n < j; ++n )
							{
								++cont;
								auxCycleCuts[cont] = verticesRota[n];
							}
							auxCycleCuts[0] = cont;

							//verifica se o corte nao contem subciclos (caso contenha, nao eh um corte valido)
							for ( int x = 1; ( ( x < cont ) && ( !corteInvalido ) ); ++x)
							{
								for ( int y = x+1; ( ( y <= cont ) && ( !corteInvalido ) ); ++y )
								{
									if ( auxCycleCuts[x] == auxCycleCuts[y] ) corteInvalido = true;
								}
							}

							//verifica se o corte viola a solucao atual (ao ser incluido no LP, sobe o bound?)
							if ( !corteInvalido )
							{
								somaY = maiorY = yCycleCuts[auxCycleCuts[1]];
								somaLHS = grafoCycleCuts[auxCycleCuts[1]][auxCycleCuts[2]];
								somaLHS += grafoCycleCuts[auxCycleCuts[2]][auxCycleCuts[1]];
								for ( int x = 2; x < cont; ++x )
								{
									somaY += yCycleCuts[auxCycleCuts[x]];
									somaLHS += grafoCycleCuts[auxCycleCuts[x]][auxCycleCuts[x+1]];
									somaLHS += grafoCycleCuts[auxCycleCuts[x+1]][auxCycleCuts[x]];
									if ( yCycleCuts[auxCycleCuts[x]] > ( maiorY + 0.00001 ) ) maiorY = yCycleCuts[auxCycleCuts[x]];
								}
								somaY += yCycleCuts[auxCycleCuts[cont]];
								if ( yCycleCuts[auxCycleCuts[cont]] > ( maiorY + 0.00001 ) ) maiorY = yCycleCuts[auxCycleCuts[cont]];

								for ( int y = 3; y <= auxCycleCuts[0]; ++y )
								{
									for ( int z = 1; z <= (y-2); ++z )
									{
										somaLHS += grafoCycleCuts[auxCycleCuts[y]][auxCycleCuts[z]];
										somaLHS += grafoCycleCuts[auxCycleCuts[z]][auxCycleCuts[y]];
									}
								}
								if ( somaLHS < ( somaY - maiorY - thresholdCycleCuts ) ) corteInvalido = true;
							}

							//verifica se o corte nao foi incluido anteriormente
							for ( int n = 0; ( ( n < cycleCutsLHS[s].size() ) && ( !corteInvalido ) ); ++n )
							{
								if ( cycleCutsLHS[s][n][0] == auxCycleCuts[0] ) //tem chances do novo corte ja ter sido incluido
								{
									int z = 1;
									while( z <= cycleCutsLHS[s][n][0] )
									{
										if ( cycleCutsLHS[s][n][z] != auxCycleCuts[z] ) break;
										++z;
									}
									if ( z > cycleCutsLHS[s][n][0] ) corteInvalido = true;
								}
							}

							if ( !corteInvalido ) //o corte eh valido e pode ser inserido no modelo
							{
								existeCorte = true;
								char* verticesCiclo = new char[cont+1];
								for ( int x = 0; x <= auxCycleCuts[0]; ++x ) verticesCiclo[x] = auxCycleCuts[x];

								IloExpr expCicloLHS(env);
								IloArray < IloExpr > expCicloRHS (env, auxCycleCuts[0]);
								for ( int a = 0; a < auxCycleCuts[0]; ++a ) expCicloRHS[a] = IloExpr(env);

								for ( int a = 0; a < qRotasL2[s]; ++a )
								{
									verticesRotaAux = ptrRotasL2[s][a]->getVertices();
									numVerticesRotaAux = verticesRotaAux.size()-2;
									for ( int x = 1; x < numVerticesRotaAux; ++x )
									{
										for ( int y = 1; y < verticesCiclo[0]; ++y )
										{
											if ( ( ( verticesRotaAux[x] == verticesCiclo[y] ) && ( verticesRotaAux[x+1] == verticesCiclo[y+1] ) ) ||
												 ( ( verticesRotaAux[x] == verticesCiclo[y+1] ) && ( verticesRotaAux[x+1] == verticesCiclo[y] ) ) )
											{
												expCicloLHS += gamma[s][a];
											}
										}
										for ( int y = 3; y <= verticesCiclo[0]; ++y )
										{
											for ( int z = 1; z <= (y-2); ++z )
											{
												if ( ( ( verticesRotaAux[x] == verticesCiclo[y] ) && ( verticesRotaAux[x+1] == verticesCiclo[z] ) ) ||
													( ( verticesRotaAux[x] == verticesCiclo[z] ) && ( verticesRotaAux[x+1] == verticesCiclo[y] ) ) )
												{
													expCicloLHS += gamma[s][a];
												}
											}
										}
									}

									//inclui os coeficientes da rota quanto ao y_i do corte
									for ( int x = 1; x <= verticesCiclo[0]; ++x )
									{
										for ( int y = 1; y <= verticesCiclo[0]; ++y )
										{
											if ( ( x != y ) && ( A_qr[s][a][verticesCiclo[y]] > 0 ) ) expCicloRHS[x-1] -= A_qr[s][a][verticesCiclo[y]]*gamma[s][a];
										}
									}
								}

								for ( int a = 0; a < qRotas_no[s]; ++a )
								{
									verticesRotaAux = ptrRotas_no[s][a]->getVertices();
									numVerticesRotaAux = verticesRotaAux.size()-2;
									for ( int x = 1; x < numVerticesRotaAux; ++x )
									{
										for ( int y = 1; y < verticesCiclo[0]; ++y )
										{
											if ( ( ( verticesRotaAux[x] == verticesCiclo[y] ) && ( verticesRotaAux[x+1] == verticesCiclo[y+1] ) ) ||
												 ( ( verticesRotaAux[x] == verticesCiclo[y+1] ) && ( verticesRotaAux[x+1] == verticesCiclo[y] ) ) )
											{
												expCicloLHS += gamma_no[s][a];
											}
										}
										for ( int y = 3; y <= verticesCiclo[0]; ++y )
										{
											for ( int z = 1; z <= (y-2); ++z )
											{
												if ( ( ( verticesRotaAux[x] == verticesCiclo[y] ) && ( verticesRotaAux[x+1] == verticesCiclo[z] ) ) ||
													( ( verticesRotaAux[x] == verticesCiclo[z] ) && ( verticesRotaAux[x+1] == verticesCiclo[y] ) ) )
												{
													expCicloLHS += gamma_no[s][a];
												}
											}
										}
									}

									//inclui os coeficientes da rota quanto ao y_i do corte
									for ( int x = 1; x <= verticesCiclo[0]; ++x )
									{
										for ( int y = 1; y <= verticesCiclo[0]; ++y )
										{
											if ( ( x != y ) && ( A_ir_no[s][a][verticesCiclo[y]] > 0 ) ) expCicloRHS[x-1] -= A_ir_no[s][a][verticesCiclo[y]]*gamma_no[s][a];
										}
									}
								}

								//inclui o corte no modelo e as estruturas para a verificacao posterior
								for ( int i = 0; i < verticesCiclo[0]; ++i )
								{
									cycleCuts[s].add( expCicloLHS + expCicloRHS[i] <= 0 );
									model.add(cycleCuts[s][indexCycleCuts]);
									++indexCycleCuts;
								}
								cycleCutsLHS[s].push_back( verticesCiclo );
							}
						}
					}
				}
			}
		}

		//procura por ciclos nas variaveis de decisao associadas ao no que invocou o corte (mesmo procedimento acima, mas para os nos)
		for ( int r = 0; r < qRotas_no[s]; ++r )
		{
			if ( ( values_no[s][r] > 0.00001 ) && ( ptrRotas_no[s][r]->ciclo ) ) //deve-se incluir no modelo o corte para esta rota
			{
				ptrRotas_no[s][r]->ciclo = false;
				verticesRota = ptrRotas_no[s][r]->getVertices();
				numVerticesRota = verticesRota.size()-2;

				//verifica cada ciclo (que nao contenha subciclos) identificado na rota
				for ( int i = 1; i < numVerticesRota; ++i )
				{
					for ( int j = i+1; j <= numVerticesRota; ++j )
					{
						if ( verticesRota[i] == verticesRota[j] )
						{
							int cont = 0;
							bool corteInvalido = false;
							for ( int n = i; n < j; ++n )
							{
								++cont;
								auxCycleCuts[cont] = verticesRota[n];
							}
							auxCycleCuts[0] = cont;

							//verifica se o corte nao contem subciclos (caso contenha, nao eh um corte valido)
							for ( int x = 1; ( ( x < cont ) && ( !corteInvalido ) ); ++x)
							{
								for ( int y = x+1; ( ( y <= cont ) && ( !corteInvalido ) ); ++y )
								{
									if ( auxCycleCuts[x] == auxCycleCuts[y] ) corteInvalido = true;
								}
							}

							//verifica se o corte viola a solucao atual (ao ser incluido no LP, sobe o bound?)
							if ( !corteInvalido )
							{
								somaY = maiorY = yCycleCuts[auxCycleCuts[1]];
								somaLHS = grafoCycleCuts[auxCycleCuts[1]][auxCycleCuts[2]];
								somaLHS += grafoCycleCuts[auxCycleCuts[2]][auxCycleCuts[1]];
								for ( int x = 2; x < cont; ++x )
								{
									somaY += yCycleCuts[auxCycleCuts[x]];
									somaLHS += grafoCycleCuts[auxCycleCuts[x]][auxCycleCuts[x+1]];
									somaLHS += grafoCycleCuts[auxCycleCuts[x+1]][auxCycleCuts[x]];
									if ( yCycleCuts[auxCycleCuts[x]] > ( maiorY + 0.00001 ) ) maiorY = yCycleCuts[auxCycleCuts[x]];
								}
								somaY += yCycleCuts[auxCycleCuts[cont]];
								if ( yCycleCuts[auxCycleCuts[cont]] > ( maiorY + 0.00001 ) ) maiorY = yCycleCuts[auxCycleCuts[cont]];

								for ( int y = 3; y <= auxCycleCuts[0]; ++y )
								{
									for ( int z = 1; z <= (y-2); ++z )
									{
										somaLHS += grafoCycleCuts[auxCycleCuts[y]][auxCycleCuts[z]];
										somaLHS += grafoCycleCuts[auxCycleCuts[z]][auxCycleCuts[y]];
									}
								}
								if ( somaLHS < ( somaY - maiorY - thresholdCycleCuts ) ) corteInvalido = true;
							}

							//verifica se o corte nao foi incluido anteriormente
							for ( int n = 0; ( ( n < cycleCutsLHS[s].size() ) && ( !corteInvalido ) ); ++n )
							{
								if ( cycleCutsLHS[s][n][0] == auxCycleCuts[0] ) //tem chances do novo corte ja ter sido incluido
								{
									int z = 1;
									while( z <= cycleCutsLHS[s][n][0] )
									{
										if ( cycleCutsLHS[s][n][z] != auxCycleCuts[z] ) break;
										++z;
									}
									if ( z > cycleCutsLHS[s][n][0] ) corteInvalido = true;
								}
							}

							if ( !corteInvalido ) //o corte eh valido e pode ser inserido no modelo
							{
								existeCorte = true;
								char* verticesCiclo = new char[cont+1];
								for ( int x = 0; x <= auxCycleCuts[0]; ++x ) verticesCiclo[x] = auxCycleCuts[x];

								IloExpr expCicloLHS(env);
								IloArray < IloExpr > expCicloRHS (env, auxCycleCuts[0]); 
								for ( int a = 0; a < auxCycleCuts[0]; ++a ) expCicloRHS[a] = IloExpr(env);

								for ( int a = 0; a < qRotasL2[s]; ++a )
								{
									verticesRotaAux = ptrRotasL2[s][a]->getVertices();
									numVerticesRotaAux = verticesRotaAux.size()-2;
									for ( int x = 1; x < numVerticesRotaAux; ++x )
									{
										for ( int y = 1; y < verticesCiclo[0]; ++y )
										{
											if ( ( ( verticesRotaAux[x] == verticesCiclo[y] ) && ( verticesRotaAux[x+1] == verticesCiclo[y+1] ) ) ||
												 ( ( verticesRotaAux[x] == verticesCiclo[y+1] ) && ( verticesRotaAux[x+1] == verticesCiclo[y] ) ) )
											{
												expCicloLHS += gamma[s][a];
											}
										}
										for ( int y = 3; y <= verticesCiclo[0]; ++y )
										{
											for ( int z = 1; z <= (y-2); ++z )
											{
												if ( ( ( verticesRotaAux[x] == verticesCiclo[y] ) && ( verticesRotaAux[x+1] == verticesCiclo[z] ) ) ||
													( ( verticesRotaAux[x] == verticesCiclo[z] ) && ( verticesRotaAux[x+1] == verticesCiclo[y] ) ) )
												{
													expCicloLHS += gamma[s][a];
												}
											}
										}
									}

									//inclui os coeficientes da rota quanto ao y_i do corte
									for ( int x = 1; x <= verticesCiclo[0]; ++x )
									{
										for ( int y = 1; y <= verticesCiclo[0]; ++y )
										{
											if ( ( x != y ) && ( A_qr[s][a][verticesCiclo[y]] > 0 ) ) expCicloRHS[x-1] -= A_qr[s][a][verticesCiclo[y]]*gamma[s][a];
										}
									}
								}

								for ( int a = 0; a < qRotas_no[s]; ++a )
								{
									verticesRotaAux = ptrRotas_no[s][a]->getVertices();
									numVerticesRotaAux = verticesRotaAux.size()-2;
									for ( int x = 1; x < numVerticesRotaAux; ++x )
									{
										for ( int y = 1; y < verticesCiclo[0]; ++y )
										{
											if ( ( ( verticesRotaAux[x] == verticesCiclo[y] ) && ( verticesRotaAux[x+1] == verticesCiclo[y+1] ) ) ||
												 ( ( verticesRotaAux[x] == verticesCiclo[y+1] ) && ( verticesRotaAux[x+1] == verticesCiclo[y] ) ) )
											{
												expCicloLHS += gamma_no[s][a];
											}
										}
										for ( int y = 3; y <= verticesCiclo[0]; ++y )
										{
											for ( int z = 1; z <= (y-2); ++z )
											{
												if ( ( ( verticesRotaAux[x] == verticesCiclo[y] ) && ( verticesRotaAux[x+1] == verticesCiclo[z] ) ) ||
													( ( verticesRotaAux[x] == verticesCiclo[z] ) && ( verticesRotaAux[x+1] == verticesCiclo[y] ) ) )
												{
													expCicloLHS += gamma_no[s][a];
												}
											}
										}
									}

									//inclui os coeficientes da rota quanto ao y_i do corte
									for ( int x = 1; x <= verticesCiclo[0]; ++x )
									{
										for ( int y = 1; y <= verticesCiclo[0]; ++y )
										{
											if ( ( x != y ) && ( A_ir_no[s][a][verticesCiclo[y]] > 0 ) ) expCicloRHS[x-1] -= A_ir_no[s][a][verticesCiclo[y]]*gamma_no[s][a];
										}
									}
								}

								//inclui o corte no modelo e as estruturas para a verificacao posterior
								for ( int i = 0; i < verticesCiclo[0]; ++i )
								{
									cycleCuts[s].add( expCicloLHS + expCicloRHS[i] <= 0 );
									model.add(cycleCuts[s][indexCycleCuts]);
									++indexCycleCuts;
								}
								cycleCutsLHS[s].push_back( verticesCiclo );
							}
						}
					}
				}
			}
		}
	}
	return existeCorte;
}


bool ModeloCplex::addRouteCuts( Grafo* g, int* qRotas_no, vector<short int*>* A_ir_no, vector < Rota* >* ptrRotas_no )
{
	bool temCorte = false;
	vector < int > verticesRota;
	float somaLHS, somaY, violacao;
	int numVerticesRota, numCandidatos, numVerticesCorte, v1, v2, v3, v4;

	//pega os valores das variaveis de decisao de todos os satellites, para nao resolver o PL sucessivamente
	IloArray<IloNumArray> values(env, nSatellites+1);
	IloArray<IloNumArray> values_no(env, nSatellites+1);
	for ( int s = 1; s <= nSatellites; ++s )
	{
		values[s] = IloNumArray(env);
		values_no[s] = IloNumArray(env);
		cplex.getValues( values[s], gamma[s] );
		if ( qRotas_no[s] > 0 ) cplex.getValues( values_no[s], gamma_no[s] );
	}

	for ( int s = 1; s <= nSatellites; ++s )
	{
		buildGraphCycleCuts( s, values, values_no, qRotas_no, ptrRotas_no );

		numCandidatos = 0; //auxCycleCuts armazena os vertices candidatos a compor o corte
		for ( int i = nSatellites+1; i < nVertices; ++i )
		{
			if ( ( grafoCycleCuts[0][i] + grafoCycleCuts[i][0] ) > 0.00001 ) auxCycleCuts[numCandidatos++] = i;
		}

		violacao = 0.0001;
		for ( int i = 0; i < numCandidatos; ++i )
		{
			somaY = yCycleCuts[auxCycleCuts[i]];
			somaLHS = grafoCycleCuts[0][auxCycleCuts[i]] + grafoCycleCuts[auxCycleCuts[i]][0];
			if ( ( somaLHS - somaY ) > violacao )
			{
				numVerticesCorte = 1;
				v1 = auxCycleCuts[i];
				violacao = somaLHS - somaY;
			}

			for ( int j = i+1; j < numCandidatos; ++j )
			{
				if ( ( demands[auxCycleCuts[i]-nSatellites] + demands[auxCycleCuts[j]-nSatellites] ) < QMin )
				{
					somaY = yCycleCuts[auxCycleCuts[i]] + yCycleCuts[auxCycleCuts[j]];
					somaLHS = grafoCycleCuts[0][auxCycleCuts[i]] + grafoCycleCuts[auxCycleCuts[i]][0] + grafoCycleCuts[0][auxCycleCuts[j]] + grafoCycleCuts[auxCycleCuts[j]][0] + 
								grafoCycleCuts[auxCycleCuts[i]][auxCycleCuts[j]] + grafoCycleCuts[auxCycleCuts[j]][auxCycleCuts[i]];
					if ( ( somaLHS - somaY  ) > violacao )
					{
						numVerticesCorte = 2;
						v1 = auxCycleCuts[i];
						v2 = auxCycleCuts[j];
						violacao = somaLHS - somaY;
					}

					for ( int k = j+1; k < numCandidatos; ++k )
					{
						if ( ( demands[auxCycleCuts[i]-nSatellites] + demands[auxCycleCuts[j]-nSatellites] + demands[auxCycleCuts[k]-nSatellites] ) < QMin )
						{
							somaY = yCycleCuts[auxCycleCuts[i]] + yCycleCuts[auxCycleCuts[j]] + yCycleCuts[auxCycleCuts[k]];
							somaLHS = grafoCycleCuts[0][auxCycleCuts[i]] + grafoCycleCuts[auxCycleCuts[i]][0] + grafoCycleCuts[0][auxCycleCuts[j]] +
								grafoCycleCuts[auxCycleCuts[j]][0] + grafoCycleCuts[0][auxCycleCuts[k]] + grafoCycleCuts[auxCycleCuts[k]][0] + 
								grafoCycleCuts[auxCycleCuts[i]][auxCycleCuts[j]] + grafoCycleCuts[auxCycleCuts[j]][auxCycleCuts[i]] + grafoCycleCuts[auxCycleCuts[i]][auxCycleCuts[k]] +
								grafoCycleCuts[auxCycleCuts[k]][auxCycleCuts[i]] + grafoCycleCuts[auxCycleCuts[j]][auxCycleCuts[k]] + grafoCycleCuts[auxCycleCuts[k]][auxCycleCuts[j]];
							if ( ( somaLHS - somaY ) > violacao )
							{
								numVerticesCorte = 3;
								v1 = auxCycleCuts[i];
								v2 = auxCycleCuts[j];
								v3 = auxCycleCuts[k];
								violacao = somaLHS - somaY;
							}

							for ( int l = k+1; l < numCandidatos; ++l )
							{
								if ( ( demands[auxCycleCuts[i]-nSatellites] + demands[auxCycleCuts[j]-nSatellites] + demands[auxCycleCuts[k]-nSatellites] + demands[auxCycleCuts[l]-nSatellites] ) < QMin )
								{
									somaY = yCycleCuts[auxCycleCuts[i]] + yCycleCuts[auxCycleCuts[j]] + yCycleCuts[auxCycleCuts[k]] + yCycleCuts[auxCycleCuts[l]];
									somaLHS = grafoCycleCuts[0][auxCycleCuts[i]] + grafoCycleCuts[auxCycleCuts[i]][0] + grafoCycleCuts[0][auxCycleCuts[j]] + grafoCycleCuts[auxCycleCuts[j]][0] +
										grafoCycleCuts[0][auxCycleCuts[k]] + grafoCycleCuts[auxCycleCuts[k]][0] + grafoCycleCuts[0][auxCycleCuts[l]] + grafoCycleCuts[auxCycleCuts[l]][0] +
										grafoCycleCuts[auxCycleCuts[i]][auxCycleCuts[j]] + grafoCycleCuts[auxCycleCuts[j]][auxCycleCuts[i]] + grafoCycleCuts[auxCycleCuts[i]][auxCycleCuts[k]] +
										grafoCycleCuts[auxCycleCuts[k]][auxCycleCuts[i]] + grafoCycleCuts[auxCycleCuts[i]][auxCycleCuts[l]] + grafoCycleCuts[auxCycleCuts[l]][auxCycleCuts[i]] + 
										grafoCycleCuts[auxCycleCuts[j]][auxCycleCuts[k]] + grafoCycleCuts[auxCycleCuts[k]][auxCycleCuts[j]] + grafoCycleCuts[auxCycleCuts[j]][auxCycleCuts[l]] +
										grafoCycleCuts[auxCycleCuts[l]][auxCycleCuts[j]] + grafoCycleCuts[auxCycleCuts[k]][auxCycleCuts[l]] + grafoCycleCuts[auxCycleCuts[l]][auxCycleCuts[k]];
									if ( ( somaLHS - somaY ) > violacao )
									{
										numVerticesCorte = 4;
										v1 = auxCycleCuts[i];
										v2 = auxCycleCuts[j];
										v3 = auxCycleCuts[k];
										v4 = auxCycleCuts[l];
										violacao = somaLHS - somaY;
									}
								}
							}
						}
					}
				}
			}
		}

		//ao terminar os 4 lacos acima, o corte mais violado (caso exista) esta armazenado nas variaveis v1--v4 e sera incluido no modelo
		if ( violacao > 0.0001 )
		{
			char* verticesCorte = new char[numVerticesCorte+1];
			verticesCorte[0] = numVerticesCorte;
			if ( numVerticesCorte == 1 ) //tem 1 vertice no corte
			{
				verticesCorte[1] = v1;
			}
			else if ( numVerticesCorte == 2 ) //tem 2 vertices no corte
			{
				verticesCorte[1] = v1;
				verticesCorte[2] = v2;
			}
			else if ( numVerticesCorte == 3 ) //tem 3 vertices no corte
			{
				verticesCorte[1] = v1;
				verticesCorte[2] = v2;
				verticesCorte[3] = v3;
			}
			else //tem 4 vertices no corte
			{
				verticesCorte[1] = v1;
				verticesCorte[2] = v2;
				verticesCorte[3] = v3;
				verticesCorte[4] = v4;
			}

			IloExpr expCorte(env);
			for ( int r = 0; r < qRotasL2[s]; ++r )
			{
				verticesRota = ptrRotasL2[s][r]->getVertices();
				numVerticesRota = verticesRota.size()-2;
				for ( int x = 1; x < numVerticesRota; ++x ) //inclui os arcos entre os vertices do corte
				{
					for ( int y = 1; y < verticesCorte[0]; ++y )
					{
						for ( int z = y+1; z <= verticesCorte[0]; ++z )
						{
							if ( ( ( verticesRota[x] == verticesCorte[y] ) && ( verticesRota[x+1] == verticesCorte[z] ) ) ||
								( ( verticesRota[x] == verticesCorte[z] ) && ( verticesRota[x+1] == verticesCorte[y] ) ) )
							{
								expCorte += gamma[s][r];
							}
						}
					}
				}

				for ( int y = 1; y <= verticesCorte[0]; ++y ) //inclui os coeficientes do corte quanto aos arcos que saem/entram no deposito e tambem os coeficientes y_i
				{
					if ( verticesRota[1] == verticesCorte[y] ) expCorte += gamma[s][r];
					if ( verticesRota[numVerticesRota] == verticesCorte[y] ) expCorte += gamma[s][r];
					if ( A_qr[s][r][verticesCorte[y]] > 0 ) expCorte -= A_qr[s][r][verticesCorte[y]]*gamma[s][r];
				}
			}

			for ( int r = 0; r < qRotas_no[s]; ++r )
			{
				verticesRota = ptrRotas_no[s][r]->getVertices();
				numVerticesRota = verticesRota.size()-2;
				for ( int x = 1; x < numVerticesRota; ++x ) //inclui os arcos entre os vertices do corte
				{
					for ( int y = 1; y < verticesCorte[0]; ++y )
					{
						for ( int z = y+1; z <= verticesCorte[0]; ++z )
						{
							if ( ( ( verticesRota[x] == verticesCorte[y] ) && ( verticesRota[x+1] == verticesCorte[z] ) ) ||
								( ( verticesRota[x] == verticesCorte[z] ) && ( verticesRota[x+1] == verticesCorte[y] ) ) )
							{
								expCorte += gamma_no[s][r];
							}
						}
					}
				}

				for ( int y = 1; y <= verticesCorte[0]; ++y ) //inclui os coeficientes do corte quanto aos arcos que saem/entram no deposito e tambem os coeficientes y_i
				{
					if ( verticesRota[1] == verticesCorte[y] ) expCorte += gamma_no[s][r];
					if ( verticesRota[numVerticesRota] == verticesCorte[y] ) expCorte += gamma_no[s][r];
					if ( A_ir_no[s][r][verticesCorte[y]] > 0 ) expCorte -= A_ir_no[s][r][verticesCorte[y]]*gamma_no[s][r];
				}
			}

			//inclui o corte no modelo e as estruturas para a verificacao posterior
			routeCuts[s].add( expCorte <= 0 );
			model.add(routeCuts[s][routeCutsLHS[s].size()]);
			routeCutsLHS[s].push_back( verticesCorte );
			expCorte.end();
			temCorte = true;
		}
	}
	return temCorte;
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
			if (values[r] >= 0.00001)
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
			if (values[r] >= 0.00001)
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


void ModeloCplex::buildGraphCycleCuts( int sat, IloArray <IloNumArray>& valuesS,  IloArray <IloNumArray>& valuesS_no, int* qRotas_no, vector < Rota* >* ptrRotas_no )
{
	vector< int > vertRota;
	int numVertRota, inicio, fim;
	memset(yCycleCuts, 0, nVertices*sizeof(double));
	for ( int i = 0; i < nVertices; ++i ) memset(grafoCycleCuts[i], 0, nVertices*sizeof(double));

	inicio = ( sat == 0 ) ? 1 : sat;
	fim = ( sat == 0 ) ? nSatellites : sat;
	for ( int s = inicio; s <= fim; ++s )
	{
		for (int r = 0; r < qRotasL2[s]; ++r) //incluo no grafo os valores das variaveis da raiz
		{
			if (valuesS[s][r] >= 0.00001)
			{
				vertRota = ptrRotasL2[s][r]->getVertices();
				numVertRota = vertRota.size()-2;
				for ( int i = 1; i < numVertRota; ++i )
				{
					grafoCycleCuts[vertRota[i]][vertRota[i+1]] += valuesS[s][r];
					yCycleCuts[vertRota[i]] += valuesS[s][r];
				}
				grafoCycleCuts[0][vertRota[1]] += valuesS[s][r];
				grafoCycleCuts[vertRota[numVertRota]][0] += valuesS[s][r];
				yCycleCuts[vertRota[numVertRota]] += valuesS[s][r];
			}
		}

		for (int r = 0; r < qRotas_no[s]; ++r) //incluo no grafo os valores das variaveis do no
		{
			if (valuesS_no[s][r] >= 0.00001)
			{
				vertRota = ptrRotas_no[s][r]->getVertices();
				numVertRota = vertRota.size()-2;
				for ( int i = 1; i < numVertRota; ++i )
				{
					grafoCycleCuts[vertRota[i]][vertRota[i+1]] += valuesS_no[s][r];
					yCycleCuts[vertRota[i]] += valuesS_no[s][r];
				}
				grafoCycleCuts[0][vertRota[1]] += valuesS_no[s][r];
				grafoCycleCuts[vertRota[numVertRota]][0] += valuesS_no[s][r];
				yCycleCuts[vertRota[numVertRota]] += valuesS_no[s][r];
			}
		}
	}
}


int ModeloCplex::mostViolatedCapacityCut(int* qRotas_no, vector < Rota* >* ptrRotas_no, double& maxViolation, int tipoCorte)
{
	//obtem os valores das variaveis de decisao e monta o grafo associado a solucao atual
	this->solveMaster();
	IloArray<IloNumArray> values(env, nSatellites+1);
	IloArray<IloNumArray> values_no(env, nSatellites+1);
	for ( int s = 1; s <= nSatellites; ++s )
	{
		values[s] = IloNumArray(env);
		values_no[s] = IloNumArray(env);
		cplex.getValues( values[s], gamma[s] );
		if ( qRotas_no[s] > 0 ) cplex.getValues( values_no[s], gamma_no[s] );
	}
	buildGraphCycleCuts(0, values, values_no, qRotas_no, ptrRotas_no);

	//verifica o indice de violacao de cada uma das restricoes e armazena o maior
	float RHS, LHS;
	maxViolation = 0;
	int indexMaxViolation = -1;

	if ( tipoCorte == 1 ) //metodo invocado para retornar capacityCut mais violada
	{
		int* sBarra = new int[1+nCustomers/2];
		for ( int i = 0; i < newCuts->Size; ++i )
		{
			LHS = 0;
			if ( newCuts->CPL[i]->IntListSize < nCustomers )
			{
				if ( newCuts->CPL[i]->IntListSize <= ( nCustomers / 2 ) ) 
				{
					for ( int j = 1; j < newCuts->CPL[i]->IntListSize; ++j )
					{
						for ( int k = j+1; k <= newCuts->CPL[i]->IntListSize; ++k )
						{
							LHS += grafoCycleCuts[newCuts->CPL[i]->IntList[j]+nSatellites][newCuts->CPL[i]->IntList[k]+nSatellites];
							LHS += grafoCycleCuts[newCuts->CPL[i]->IntList[k]+nSatellites][newCuts->CPL[i]->IntList[j]+nSatellites];
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

					//calcula o LHS com os valores da solucao atual para verificar a violacao
					for ( int j = 0; j < tamSBarra; ++j )
					{
						LHS += 0.5*grafoCycleCuts[0][sBarra[j]+nSatellites];
						LHS += 0.5*grafoCycleCuts[sBarra[j]+nSatellites][0];
					}
					for ( int j = 1; j <= newCuts->CPL[i]->IntListSize; ++j )
					{
						LHS -= 0.5*grafoCycleCuts[0][newCuts->CPL[i]->IntList[j]+nSatellites];
						LHS -= 0.5*grafoCycleCuts[newCuts->CPL[i]->IntList[j]+nSatellites][0];
					}
					for ( int j = 0; j < tamSBarra-1; ++j )
					{
						for ( int k = j+1; k < tamSBarra; ++k )
						{
							LHS += grafoCycleCuts[sBarra[j]+nSatellites][sBarra[k]+nSatellites];
							LHS += grafoCycleCuts[sBarra[k]+nSatellites][sBarra[j]+nSatellites];
						}
					}
					RHS = ( tamSBarra - ( newCuts->CPL[i]->IntListSize - newCuts->CPL[i]->RHS ) );
				}

				if ( ( LHS - RHS ) > ( maxViolation + 0.0001 ) )
				{
					indexMaxViolation = i;
					maxViolation = ( LHS - RHS );
				}
			}
		}
		delete [] sBarra;
	}
	else if ( tipoCorte == 2 ) //metodo invocado para retornar combCut mais violada
	{
		bool depositoDentro;
		int tmp, inicioT, fimT, tamComplemento;
		int* conjComplemento = new int[nCustomers];
		for ( int i = 0; i < newCuts->Size; ++i )
		{
			LHS = 0;
			tamComplemento = 0;
			depositoDentro = false;
			for ( int h = 1; h <= newCuts->CPL[i]->IntListSize; h++ )
			{
				for ( int j = 1; ( ( h == 1 ) && ( j <= nCustomers ) ); ++j ) //avalia o complemento do Handle (apenas uma vez)
				{
					for ( tmp = 1; tmp <= newCuts->CPL[i]->IntListSize; ++tmp )
					{
						if ( newCuts->CPL[i]->IntList[tmp] == j ) break;
					}
					if ( tmp > newCuts->CPL[i]->IntListSize )
					{
						conjComplemento[tamComplemento++] = j;
					}
				}
				for ( int j = 0; j < tamComplemento; ++j )
				{
					if ( newCuts->CPL[i]->IntList[h] > nCustomers ) //considera o deposito estando dentro de H
					{
						depositoDentro = true;
						LHS += grafoCycleCuts[0][conjComplemento[j]+nSatellites];
						LHS += grafoCycleCuts[conjComplemento[j]+nSatellites][0];
					}
					else
					{
						LHS += grafoCycleCuts[newCuts->CPL[i]->IntList[h]+nSatellites][conjComplemento[j]+nSatellites];
						LHS += grafoCycleCuts[conjComplemento[j]+nSatellites][newCuts->CPL[i]->IntList[h]+nSatellites];
					}
				}
			}
			//considera separadamente o deposito como estando fora de H
			if ( !depositoDentro )
			{
				for ( int h = 1; h <= newCuts->CPL[i]->IntListSize; h++ )
				{
					LHS += grafoCycleCuts[newCuts->CPL[i]->IntList[h]+nSatellites][0];
					LHS += grafoCycleCuts[0][newCuts->CPL[i]->IntList[h]+nSatellites];
				}
			}

			//atribui na matriz valor 1 para aquelas arestas q estao no corte associado aos conjuntos T
			for ( int t = 1; t <= newCuts->CPL[i]->Key; t++ )
			{
				tamComplemento = 0;
				depositoDentro = false;
				inicioT = newCuts->CPL[i]->ExtList[t];
				if (t == newCuts->CPL[i]->Key) fimT = newCuts->CPL[i]->ExtListSize;
				else fimT = newCuts->CPL[i]->ExtList[t+1] - 1;

				for ( int j = inicioT; j <= fimT; ++j )
				{
					for ( int k = 1; ( ( j == inicioT ) && ( k <= nCustomers ) ); ++k )
					{
						for ( tmp = inicioT; tmp <= fimT; ++tmp )
						{
							if ( newCuts->CPL[i]->ExtList[tmp] == k ) break;
						}
						if ( tmp > fimT )
						{
							conjComplemento[tamComplemento++] = k;
						}
					}
					for ( int k = 0; k < tamComplemento; ++k )
					{
						if ( newCuts->CPL[i]->ExtList[j] > nCustomers ) //considera o deposito estando dentro de T
						{
							depositoDentro = true;
							LHS += grafoCycleCuts[0][conjComplemento[k]+nSatellites];
							LHS += grafoCycleCuts[conjComplemento[k]+nSatellites][0];
						}
						else
						{
							LHS += grafoCycleCuts[newCuts->CPL[i]->ExtList[j]+nSatellites][conjComplemento[k]+nSatellites];
							LHS += grafoCycleCuts[conjComplemento[k]+nSatellites][newCuts->CPL[i]->ExtList[j]+nSatellites];
						}
					}
				}
				//considera separadamente o deposito como estando fora deste Tooth
				if ( !depositoDentro )
				{
					for ( int j = inicioT; j <= fimT; ++j )
					{
						LHS += grafoCycleCuts[newCuts->CPL[i]->ExtList[j]+nSatellites][0];
						LHS += grafoCycleCuts[0][newCuts->CPL[i]->ExtList[j]+nSatellites];
					}
				}
			}

			//verifica o grau de violacao da restricao
			if ( ( newCuts->CPL[i]->RHS - LHS ) > maxViolation + 0.0001 )
			{
				indexMaxViolation = i;
				maxViolation = newCuts->CPL[i]->RHS - LHS;
			}
		}
		delete [] conjComplemento;
	}
	else if ( tipoCorte == 3 ) //metodo invocado para retornar multiStarCut mais violada
	{
		int c;
		float tmp;
		for ( int i = 0; i < newCuts->Size; ++i )
		{
			LHS = tmp = 0;
			//define as arestas de x(\delta(N))
			for ( int v = 1; v <= ( nCustomers+1 ); ++v )
			{
				for ( c = 1; c <= newCuts->CPL[i]->IntListSize; ++c )
				{
					if ( newCuts->CPL[i]->IntList[c] == v ) break;
				}
				if ( c > newCuts->CPL[i]->IntListSize )
				{
					for ( c = 1; c <= newCuts->CPL[i]->IntListSize; ++c )
					{
						if ( v > nCustomers ) //se o deposito estiver fora de N
						{
							LHS += grafoCycleCuts[0][newCuts->CPL[i]->IntList[c]+nSatellites];
							LHS += grafoCycleCuts[newCuts->CPL[i]->IntList[c]+nSatellites][0];
						}
						else if ( newCuts->CPL[i]->IntList[c] > nCustomers ) //se o deposito estiver em N
						{
							LHS += grafoCycleCuts[0][v+nSatellites];
							LHS += grafoCycleCuts[v+nSatellites][0];
						}
						else
						{
							LHS += grafoCycleCuts[v+nSatellites][newCuts->CPL[i]->IntList[c]+nSatellites];
							LHS += grafoCycleCuts[newCuts->CPL[i]->IntList[c]+nSatellites][v+nSatellites];
						}
					}
				}
			}
			LHS *= newCuts->CPL[i]->B;

			//define as arestas de x(C:T)
			for ( c = 1; c <= newCuts->CPL[i]->CListSize; ++c )
			{
				for ( int t = 1; t <= newCuts->CPL[i]->ExtListSize; ++t )
				{
					if ( newCuts->CPL[i]->CList[c] > nCustomers )
					{
						tmp += grafoCycleCuts[0][newCuts->CPL[i]->ExtList[t]+nSatellites];
						tmp += grafoCycleCuts[newCuts->CPL[i]->ExtList[t]+nSatellites][0];
					}
					else if ( newCuts->CPL[i]->ExtList[t] > nCustomers ) 
					{
						tmp += grafoCycleCuts[0][newCuts->CPL[i]->CList[c]+nSatellites];
						tmp += grafoCycleCuts[newCuts->CPL[i]->CList[c]+nSatellites][0];
					}
					else
					{
						tmp += grafoCycleCuts[newCuts->CPL[i]->CList[c]+nSatellites][newCuts->CPL[i]->ExtList[t]+nSatellites];
						tmp += grafoCycleCuts[newCuts->CPL[i]->ExtList[t]+nSatellites][newCuts->CPL[i]->CList[c]+nSatellites];
					}
				}
			}
			LHS -= ( newCuts->CPL[i]->A * tmp );

			if ( ( newCuts->CPL[i]->L - LHS ) > ( maxViolation + 0.0001 ) )
			{
				indexMaxViolation = i;
				maxViolation = newCuts->CPL[i]->L - LHS;
			}
		}
	}
	return indexMaxViolation;
}



void ModeloCplex::setLimitePrimal(double lP){
	limitePrimal = lP;
}


double ModeloCplex::getLimitePrimal(){
	return limitePrimal;
}


void ModeloCplex::processaBase( int numItOpt, int* qRotas_no, vector < short int* >* a_ir_no, vector < Rota* >* ptrRotas_no )
{
	if ( qRotas_no != NULL )
	{
		for ( int s = 1; s <= nSatellites; ++s )
		{
			vector < int > v;
			int n = 0, limite = qRotas_no[s] / 2;
			IloNumVarArray gammaAux = IloNumVarArray(env);
			for ( int r = 0; r < limite; ++r )
			{
				if ( a_ir_no[s][r][0] == 1 )
				{
					v.push_back(r);
					gammaAux.add(gamma_no[s][r]);
				}
				else
				{
					gamma_no[s][r].end();
				}
			}
			for ( int r = limite; r < qRotas_no[s]; ++r ) gammaAux.add( gamma_no[s][r] );

			if ( v.size() > 0 )
			{
				if ( v[v.size()-1] != ( limite-1 ) ) v.push_back( limite );
				for ( int x = 0; x < v.size(); ++x )
				{
					if ( n < v[x] )
					{
						ptrRotas_no[s].erase(ptrRotas_no[s].begin()+n, ptrRotas_no[s].begin()+(v[x]));
						a_ir_no[s].erase(a_ir_no[s].begin()+n, a_ir_no[s].begin()+(v[x]));
						for ( int y = x+1; y < v.size(); ++y ) v[y] -= ( v[x] - n );
					}
					++n;
				}
			}
			else
			{
				ptrRotas_no[s].erase(ptrRotas_no[s].begin(), ptrRotas_no[s].begin()+limite);
				a_ir_no[s].erase(a_ir_no[s].begin(), a_ir_no[s].begin()+limite);
			}

			gamma_no[s].end();
			gamma_no[s] = IloNumVarArray( env );
			gamma_no[s].add( gammaAux );
			qRotas_no[s] = ptrRotas_no[s].size();
		}
	}
	else
	{
		int threshold = numItOpt*0.17;
		for ( int s = 1; s <= nSatellites; ++s )
		{
			int n = 0;
			vector < int > v;
			IloNumVarArray gammaAux = IloNumVarArray(env);
			for ( int r = 0; r < qRotasL2[s]; ++r )
			{
				if ( A_qr[s][r][0] > -threshold )
				{
					v.push_back(r);
					gammaAux.add(gamma[s][r]);
				}
				else
				{
					gamma[s][r].end();
				}
			}

			if ( v.size() > 0 )
			{
				if ( v[v.size()-1] != ( qRotasL2[s]-1 ) ) v.push_back(qRotasL2[s]);
				for ( int x = 0; x < v.size(); ++x )
				{
					if ( n < v[x] )
					{
						ptrRotasL2[s].erase(ptrRotasL2[s].begin()+n, ptrRotasL2[s].begin()+(v[x]));
						A_qr[s].erase(A_qr[s].begin()+n, A_qr[s].begin()+(v[x]));
						for ( int y = x+1; y < v.size(); ++y ) v[y] -= ( v[x] - n );
					}
					++n;
				}
			}
			else
			{
				A_qr[s].clear();
				ptrRotasL2[s].clear();
			}

			gamma[s].end();
			gamma[s] = IloNumVarArray(env);
			gamma[s].add(gammaAux);
			qRotasL2[s] = ptrRotasL2[s].size();
		}
	}
}


void ModeloCplex::imprimeRotasBasicas( vector < Rota* > rLambda, Grafo* g ){
	IloNumArray values(env);

	printf("\n");
	double totalLambda = 0;
	cplex.getValues(values, lambda);
	for (int r = 0; r < values.getSize(); ++r)
	{
		if (values[r] >= 0.0001)
		{
			printf("lambda[%d](%0.03f): ", r, values[r]);
			rLambda[r]->imprimir();
			totalLambda += values[r];
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
		for ( int r = 0; r < quantRotasL1; ++r )
		{
			if (values[r] >= 0.00001)
			{
				printf("delta[%d][%d](%0.02f)\n", s, r, values[r]);
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
}

int ModeloCplex::getNumCuts()
{
	int totalCuts = capacityCutsLHS.size() + combCutsLHS.size() + multiStarCutsLHS.size();
        for ( int s = 1; s <= nSatellites; ++s ) totalCuts += cycleCutsLHS[s].size();
	return totalCuts;
}
