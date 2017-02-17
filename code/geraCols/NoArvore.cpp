#include "NoArvore.h"

Grafo* NoArvore::G;
int NoArvore::nVertices;
int NoArvore::nSatellites;
int NoArvore::totalNosAtivos;
int NoArvore::totalNosCriados;
ModeloCplex* NoArvore::mCplex;
ptrNo NoArvore::listaNosProcessados;

NoArvore::~NoArvore(){
	--totalNosAtivos;
	for (int s = 1; s <= nSatellites; ++s )
	{
		for (int r = 0; r < qRotas_no[s]; ++r)
		{
			delete [] a_ir_no[s][r];
			if ( ptrRotas_no[s][r]->decrNumApontadores() ) delete ptrRotas_no[s][r];
		}
	}
	delete [] a_ir_no;
	delete [] qRotas_no;
	delete [] ptrRotas_no;

	for (int i = 0; i < arcosBranching.size(); ++i)
	{
		delete [] arcosBranching[i];
	}
}

NoArvore::NoArvore(ModeloCplex* modCplex, Grafo* g, vector < ptrNo >& nosCriados, double lD){
	G = g;
	altura = 0;
	indiceNo = 0;
	separarCortes = 0;
	totalNosAtivos = 0;
	totalNosCriados = 0;
	mCplex = modCplex;
	nSatellites = g->getNumSatellites();
	nVertices = g->getNumVertices();
	listaNosProcessados = NULL;
	limiteDual = lD;
	int i, j, s, rhs;

	qRotas_no = new int[nSatellites+1];
	ptrRotas_no = new vector <Rota*>[nSatellites+1];
	a_ir_no = new vector <short int*>[nSatellites+1];
	for ( int s = 0; s <= nSatellites; ++s ) qRotas_no[s] = 0;

	//instancia os objetos do cplex para manipular as restricoes e variaveis do branching
	mCplex->artificiais = IloNumVarArray( mCplex->env );
	mCplex->constraintsArvoreL1 = IloRangeArray( mCplex->env );
	mCplex->constraintsArvoreL2 = IloRangeArray( mCplex->env );
	mCplex->gamma_no = IloArray < IloNumVarArray > ( mCplex->env, nSatellites+1 );
	for ( int s = 0; s <= nSatellites; ++s ) mCplex->gamma_no[s] = IloNumVarArray( mCplex->env );
	//as constraints da arvore da quantidade de veiculos que partem de um satellite sao inicializadas no ModeloCplex


	//variaveis i, j, s e rhs passadas como referencia. armazenarao o arco a ser realizado o branching.
	if ( defineVariavelBranching(i, j, s, rhs ) )
	{
		//nao precisa deletar branch0 e branch1, pois serao armazenados pelos nos e deletados em seus destrutores
		short int* branch0 = new short int[4];
		short int* branch1 = new short int[4];
		if ( s == 0 ) //caso o branching seja em variavel lambda
		{	
			branch0[0] = i; branch0[1] = -1; branch0[2] = 0; branch0[3] = (rhs-1);
			branch1[0] = i; branch1[1] = 1; branch1[2] = 0; branch1[3] = rhs;
		}
		else if ( s < 0 ) //caso o branching seja no numero de veiculos que saem do satellite s = -s
		{
			branch0[0] = 0; branch0[1] = -1; branch0[2] = s; branch0[3] = (rhs-1);
			branch1[0] = 0; branch1[1] = 1; branch1[2] = s; branch1[3] = rhs;
		}
		else //caso o branching seja em variavel gamma
		{
			branch0[0] = i; branch0[1] = j; branch0[2] = s; branch0[3] = 0;
			branch1[0] = i; branch1[1] = j; branch1[2] = s; branch1[3] = 1;
		}
		nosCriados.push_back(new NoArvore(this, branch0));
		nosCriados.push_back(new NoArvore(this, branch1));
	}
}

NoArvore::NoArvore(ptrNo pai, short int* arcoBranching){
	short int *tmp;
	++totalNosAtivos;
	indiceNo = ++totalNosCriados;
	qRotas_no = new int[nSatellites+1];
	ptrRotas_no = new vector <Rota*>[nSatellites+1];
	a_ir_no = new vector <short int*>[nSatellites+1];

	//copia (herda) as informacoes do pai
	this->pai = pai;
	altura = pai->altura+1;
	limiteDual = pai->limiteDual;
	separarCortes = pai->separarCortes;
	for (int s = 0; s <= nSatellites; ++s)
	{
		qRotas_no[s] = pai->qRotas_no[s];
		for (int r = 0; r < qRotas_no[s]; ++r)
		{
			tmp = new short int[nVertices];
			for (int i = 0; i < nVertices; ++i)
			{
				tmp[i] = pai->a_ir_no[s][r][i];
			}
			a_ir_no[s].push_back(tmp);
			ptrRotas_no[s].push_back(pai->ptrRotas_no[s][r]);
			ptrRotas_no[s][r]->incrNumApontadores();
		}
	}

	//herda os arcos todos os arcos fixos do pai
	int arcosBranchingSize = pai->arcosBranching.size();
	for (int i = 0; i < arcosBranchingSize; i++)
	{
		tmp = new short int[4];
		tmp[0] = pai->arcosBranching[i][0];
		tmp[1] = pai->arcosBranching[i][1];
		tmp[2] = pai->arcosBranching[i][2];
		tmp[3] = pai->arcosBranching[i][3];
		arcosBranching.push_back(tmp);
	}
	//por fim, inclue o arco de branching passado como parametro (o ultimo definido para executar o branching)
	arcosBranching.push_back(arcoBranching);
}

double NoArvore::executaBranching(vector<ptrNo>& nosCriados, char opSub){
	//limpa as informacoes do no processado anteriormente no ModeloCplex (variaveis de decisao, artificiais e restricoes da arvore)
	mCplex->artificiais.endElements();
	mCplex->constraintsArvoreL1.endElements();
	mCplex->constraintsArvoreL2.endElements();
	for ( int s = 1; s <= nSatellites; ++s )
	{
		mCplex->gamma_no[s].endElements();
		mCplex->constraintsArvoreVeic[s].endElements();
	}

	//insere no ModeloCplex as informacoes pertinentes ao no
	setVariaveisDecisao_no();
	setRestricoes_no();

	//Uma vez que o modelo deste no esta montado, eh possivel chegar na raiz atraves da geracao de colunas
	//ao fixar os arcos, eh possivel que o modelo se torne inviavel, neste caso, interrompe a geracao de filhos (poda por inviabilidade)
	if ( alcancaRaiz_no( opSub ) )
	{
		//caso o limite dual do no seja maior ou igual ao limite o primal, nao gera os filhos (poda por limite dual)
		if ( limiteDual < ModeloCplex::getLimitePrimal() )
		{
			//variaveis i e j passadas como referencia (arco a ser realizado o branching)
			int i, j, s, rhs;
			if ( defineVariavelBranching( i, j, s, rhs ) )
			{
				//nao precisa deletar branch0 e branch1, pois serao armazenados pelos nos e deletados em seus destrutores
				short int* branch0 = new short int[4];
				short int* branch1 = new short int[4];
				if ( s == 0 ) //caso o branching seja em variavel lambda
				{	
					branch0[0] = i; branch0[1] = -1; branch0[2] = 0; branch0[3] = (rhs-1);
					branch1[0] = i; branch1[1] = 1; branch1[2] = 0; branch1[3] = rhs;
				}
				else if ( s < 0 ) //caso o branching seja no numero de veiculos que saem do satellite s = -s
				{
					branch0[0] = 0; branch0[1] = -1; branch0[2] = s; branch0[3] = (rhs-1);
					branch1[0] = 0; branch1[1] = 1; branch1[2] = s; branch1[3] = rhs;
				}
				else //caso o branching seja em variavel gamma
				{
					branch0[0] = i; branch0[1] = j; branch0[2] = s; branch0[3] = 0;
					branch1[0] = i; branch1[1] = j; branch1[2] = s; branch1[3] = 1;
				}
				nosCriados.push_back(new NoArvore(this, branch0));
				nosCriados.push_back(new NoArvore(this, branch1));
			}
			else return limiteDual; //poda por otimalidade
		}
		return 999999;
	}
	return 999999;
}

void NoArvore::setVariaveisDecisao_no(){
	float indexCol;
	vector < int > vertRota;
	int numCombCuts = mCplex->combCutsLHS.size();
	int numCapacityCuts = mCplex->capacityCutsLHS.size();
	int numMultiStarCuts = mCplex->multiStarCutsLHS.size();
	int coef, count, numVertAtualizar, numEdgesCut, numCycleCuts, numRouteCuts;

	//para cada rota partindo de um satellite inclue uma rota no ModeloCplex
	for ( int s = 1; s <= nSatellites; ++s )
	{
		numCycleCuts = mCplex->cycleCutsLHS[s].size();
		numRouteCuts = mCplex->routeCutsLHS[s].size();

		for (int r = 0; r < qRotas_no[s]; ++r)
		{
			IloNumColumn col = mCplex->objCost(ptrRotas_no[s][r]->getCusto());
			col += mCplex->constraint2(1);
			col += mCplex->constraints3[s-1](1);

			coef = 0; count = 0;
			for ( int i = nSatellites+1; i < nVertices; ++i, ++count )
			{
				if ( a_ir_no[s][r][i] > 0 )
				{
					col += mCplex->constraints4[count]( a_ir_no[s][r][i] );
					coef += ( a_ir_no[s][r][i] * G->getCargaVertice(i) );
				}
			}
			col += mCplex->constraints5[s-1](coef);

			vertRota = ptrRotas_no[s][r]->getVertices();
			numVertAtualizar = vertRota.size()-2;
			//Insere os indices da variavel nas restricoes de corte de capacidade
			for ( int i = 0; i < numCapacityCuts; ++i )
			{
				indexCol = 0;
				for ( int v = 1; v < numVertAtualizar; ++v )
				{
					if ( mCplex->capacityCutsLHS[i][vertRota[v]][vertRota[v+1]] != 0 ) indexCol += mCplex->capacityCutsLHS[i][vertRota[v]][vertRota[v+1]];
					if ( mCplex->capacityCutsLHS[i][vertRota[v+1]][vertRota[v]] != 0 ) indexCol += mCplex->capacityCutsLHS[i][vertRota[v+1]][vertRota[v]];
				}
				if ( mCplex->capacityCutsLHS[i][0][vertRota[1]] != 0 ) indexCol += mCplex->capacityCutsLHS[i][0][vertRota[1]];
				if ( mCplex->capacityCutsLHS[i][0][vertRota[numVertAtualizar]] != 0 ) indexCol += mCplex->capacityCutsLHS[i][0][vertRota[numVertAtualizar]];

				if ( ( indexCol > 0.00001 ) || ( indexCol < -0.00001 ) ) col += mCplex->capacityCuts[i](indexCol);
			}

			//insere a variavel nas restricoes dos cortes comb
			for ( int i = 0; i < numCombCuts; ++i )
			{
				indexCol = 0;
				for ( int v = 1; v < numVertAtualizar; ++v )
				{
					if ( mCplex->combCutsLHS[i][vertRota[v]][vertRota[v+1]] > 0 ) indexCol += mCplex->combCutsLHS[i][vertRota[v]][vertRota[v+1]];
					if ( mCplex->combCutsLHS[i][vertRota[v+1]][vertRota[v]] > 0 ) indexCol += mCplex->combCutsLHS[i][vertRota[v+1]][vertRota[v]];
				}
				if ( mCplex->combCutsLHS[i][0][vertRota[1]] > 0 ) indexCol += mCplex->combCutsLHS[i][0][vertRota[1]];
				if ( mCplex->combCutsLHS[i][0][vertRota[numVertAtualizar]] > 0 ) indexCol += mCplex->combCutsLHS[i][0][vertRota[numVertAtualizar]];

				if ( indexCol > 0.00001 ) col += mCplex->combCuts[i](indexCol);
			}

			//insere a coluna nas restricoes dos cortes multiStar
			for ( int i = 0; i < numMultiStarCuts; ++i )
			{
				indexCol = 0;
				for ( int v = 1; v < numVertAtualizar; ++v )
				{
					if ( mCplex->multiStarCutsLHS[i][vertRota[v]][vertRota[v+1]] != 0 ) indexCol += mCplex->multiStarCutsLHS[i][vertRota[v]][vertRota[v+1]];
					if ( mCplex->multiStarCutsLHS[i][vertRota[v+1]][vertRota[v]] != 0 ) indexCol += mCplex->multiStarCutsLHS[i][vertRota[v+1]][vertRota[v]];
				}
				if ( mCplex->multiStarCutsLHS[i][0][vertRota[1]] != 0 ) indexCol += mCplex->multiStarCutsLHS[i][0][vertRota[1]];
				if ( mCplex->multiStarCutsLHS[i][0][vertRota[numVertAtualizar]] != 0 ) indexCol += mCplex->multiStarCutsLHS[i][0][vertRota[numVertAtualizar]];

				if ( ( indexCol > 0.00001 ) || ( indexCol < -0.00001 ) ) col += mCplex->multiStarCuts[i](indexCol);
			}

			//insere a coluna nas restricoes dos cortes outfork
			for ( int i = 0; i < numCycleCuts; ++i )
			{
				indexCol = 0;
				for ( int v = 1; v < numVertAtualizar; ++v )
				{
					for ( int y = 1; y < mCplex->cycleCutsLHS[s][i][0]; ++y )
					{
						if ( ( ( vertRota[v] == mCplex->cycleCutsLHS[s][i][y] ) && ( vertRota[v+1] == mCplex->cycleCutsLHS[s][i][y+1] ) ) ||
							 ( ( vertRota[v] == mCplex->cycleCutsLHS[s][i][y+1] ) && ( vertRota[v+1] == mCplex->cycleCutsLHS[s][i][y] ) ) ) ++indexCol;
					}
					for ( int y = 3; y <= mCplex->cycleCutsLHS[s][i][0]; ++y )
					{
						for ( int z = 1; z <= (y-2); ++z )
						{
							if ( ( ( vertRota[v] == mCplex->cycleCutsLHS[s][i][y] ) && ( vertRota[v+1] == mCplex->cycleCutsLHS[s][i][z] ) ) ||
								( ( vertRota[v] == mCplex->cycleCutsLHS[s][i][z] ) && ( vertRota[v+1] == mCplex->cycleCutsLHS[s][i][y] ) ) ) ++indexCol;
						}
					}
				}
				for ( int y = 2; y <= mCplex->cycleCutsLHS[s][i][0]; ++y )
				{
					if ( a_ir_no[s][r][mCplex->cycleCutsLHS[s][i][y]] > 0 ) indexCol -= a_ir_no[s][r][mCplex->cycleCutsLHS[s][i][y]];
				}

				if ( ( indexCol > 0.00001 ) || ( indexCol < -0.00001 ) ) col += mCplex->cycleCuts[s][i](indexCol);
			}

			//insere a coluna nas restricoes dos routeFeasibilityCuts
			for ( int i = 0; i < numRouteCuts; ++i )
			{
				indexCol = 0;
				for ( int v = 1; v < numVertAtualizar; ++v )
				{
					for ( int y = 1; y < mCplex->routeCutsLHS[s][i][0]; ++y )
					{
						for ( int z = y+1; z <= mCplex->routeCutsLHS[s][i][0]; ++z )
						{
							if ( ( ( vertRota[v] == mCplex->routeCutsLHS[s][i][y] ) && ( vertRota[v+1] == mCplex->routeCutsLHS[s][i][z] ) ) ||
								( ( vertRota[v] == mCplex->routeCutsLHS[s][i][z] ) && ( vertRota[v+1] == mCplex->routeCutsLHS[s][i][y] ) ) ) ++indexCol;
						}
					}
				}
				for ( int y = 1; y <= mCplex->routeCutsLHS[s][i][0]; ++y )
				{
					if ( vertRota[1] == mCplex->routeCutsLHS[s][i][y] ) ++indexCol;
					if ( vertRota[numVertAtualizar] == mCplex->routeCutsLHS[s][i][y] ) ++indexCol;
					if ( a_ir_no[s][r][mCplex->routeCutsLHS[s][i][y]] > 0 ) indexCol -= a_ir_no[s][r][mCplex->routeCutsLHS[s][i][y]];
				}

				if ( ( indexCol > 0.00001 ) || ( indexCol < -0.00001 ) ) col += mCplex->routeCuts[s][i](indexCol);
			}

			char nome[20];
			sprintf(nome, "gNo_%01d_%02d", s, r);
			mCplex->gamma_no[s].add(IloNumVar(col, 0, 1, ILOFLOAT, nome));
			col.end();
		}
	}
}

void NoArvore::setRestricoes_no() {
	vector<int> verticesRota;
	int sFixo, vFixoI, vFixoJ, quantR, quantV, coefColuna, numVisitasArco, numRestrGamma = 0, numRestrLambda = 0;

	int arcosBranchingSize = arcosBranching.size();
	for ( int i = 0; i < arcosBranchingSize; ++i )
	{
		if ( arcosBranching[i][2] == 0 ) //restricoes de branching associadas as variaveis lambda
		{
			if ( arcosBranching[i][1] < 0 )
			{
				mCplex->constraintsArvoreL1.add(mCplex->lambda[arcosBranching[i][0]] <= arcosBranching[i][3]);
				mCplex->model.add(mCplex->constraintsArvoreL1[numRestrLambda]);

				//insere a variavel artificial associada a esta restricao
				IloNumColumn col = mCplex->objCost(ModeloCplex::limitePrimal);
				col += mCplex->constraintsArvoreL1[numRestrLambda](-1);
				mCplex->artificiais.add(IloNumVar(col, 0, +IloInfinity, ILOFLOAT));
				++numRestrLambda;
				col.end();
			}
			else
			{
				mCplex->constraintsArvoreL1.add(mCplex->lambda[arcosBranching[i][0]] >= arcosBranching[i][3]);
				mCplex->model.add(mCplex->constraintsArvoreL1[numRestrLambda]);

				//insere a variavel artificial associada a esta restricao
				IloNumColumn col = mCplex->objCost(ModeloCplex::limitePrimal);
				col += mCplex->constraintsArvoreL1[numRestrLambda](1);
				mCplex->artificiais.add(IloNumVar(col, 0, +IloInfinity, ILOFLOAT));
				++numRestrLambda;
				col.end();
			}
		}

		if ( arcosBranching[i][2] < 0 ) //restricoes de branching associadas ao numero de veiculos que sai de cada satellite
		{
			sFixo = -arcosBranching[i][2];
			IloExpr exp = IloExpr(mCplex->env);
			quantV = mCplex->constraintsArvoreVeic[sFixo].getSize();

			//obtem as rotas gamma[s] da raiz e do no 
			for ( int r = 0; r < mCplex->qRotasL2[sFixo]; ++r ) exp += mCplex->gamma[sFixo][r];
			for ( int r = 0; r < qRotas_no[sFixo]; ++r ) exp += mCplex->gamma_no[sFixo][r];

			if ( arcosBranching[i][1] < 0 )
			{
				mCplex->constraintsArvoreVeic[sFixo].add( exp <= arcosBranching[i][3] );
				mCplex->model.add(mCplex->constraintsArvoreVeic[sFixo][quantV]);

				//insere a variavel artificial associada a esta restricao
				IloNumColumn col = mCplex->objCost(ModeloCplex::limitePrimal);
				col += mCplex->constraintsArvoreVeic[sFixo][quantV](-1);
				mCplex->artificiais.add(IloNumVar(col, 0, +IloInfinity, ILOFLOAT));
				col.end();
			}
			else
			{
				mCplex->constraintsArvoreVeic[sFixo].add( exp >= arcosBranching[i][3] );
				mCplex->model.add(mCplex->constraintsArvoreVeic[sFixo][quantV]);

				//insere a variavel artificial associada a esta restricao
				IloNumColumn col = mCplex->objCost(ModeloCplex::limitePrimal);
				col += mCplex->constraintsArvoreVeic[sFixo][quantV](1);
				mCplex->artificiais.add(IloNumVar(col, 0, +IloInfinity, ILOFLOAT));
				col.end();
			}
		}

		if ( arcosBranching[i][2] > 0 ) //restricoes de branching associadas as variaveis gamma
		{
			vFixoI = arcosBranching[i][0];
			vFixoJ = arcosBranching[i][1];
			sFixo = arcosBranching[i][2];
			IloExpr exp = IloExpr(mCplex->env);

			//insere na restricao as variaveis de decisao do no raiz (validas para qualquer no da arvore)
			quantR = mCplex->ptrRotasL2[sFixo].size();
			for (int r = 0; r < quantR; ++r)
			{
				verticesRota = mCplex->ptrRotasL2[sFixo][r]->getVertices();
				quantV = verticesRota.size()-1;
				numVisitasArco = 0;
				for (int v = 0; v < quantV; ++v)
				{
					if ( ( verticesRota[v] == vFixoI) && ( verticesRota[v+1] == vFixoJ ) )
					{
						++numVisitasArco;
					}
				}
				if ( numVisitasArco > 0 ) exp += numVisitasArco * mCplex->gamma[sFixo][r];
			}

			//insere na restricao as variaveis de decisao deste no
			quantR = ptrRotas_no[sFixo].size();
			for ( int r = 0; r < quantR; ++r )
			{
				verticesRota = ptrRotas_no[sFixo][r]->getVertices();
				quantV = verticesRota.size()-1;
				numVisitasArco = 0;
				for ( int v = 0; v < quantV; ++v )
				{
					if ( ( verticesRota[v] == vFixoI ) && ( verticesRota[v+1] == vFixoJ ) )
					{
						++numVisitasArco;
					}
				}
				if ( numVisitasArco > 0 ) exp += numVisitasArco * mCplex->gamma_no[sFixo][r];
			}

			//a variavel coefColuna sera usada para controlar o indice da variavel artificial nesta restricao
			if ( arcosBranching[i][3] == 1 )
			{
				mCplex->constraintsArvoreL2.add(exp == 1);
				coefColuna = 1;
			}
			else
			{
				mCplex->constraintsArvoreL2.add(exp == 0);
				coefColuna = -1;
			}

			//apos ter construido a restricao, a insere no modelo
			char nome[20];
			sprintf(nome, "s%d_%d_%d", sFixo, vFixoI, vFixoJ);
			mCplex->constraintsArvoreL2[numRestrGamma].setName(nome);
			mCplex->model.add(mCplex->constraintsArvoreL2[numRestrGamma]);
			exp.end();

			//insere a variavel artificial associada a esta restricao
			sprintf(nome, "artifG%d", numRestrGamma);
			IloNumColumn col = mCplex->objCost(ModeloCplex::limitePrimal);
			col += mCplex->constraintsArvoreL2[numRestrGamma](coefColuna);
			mCplex->artificiais.add(IloNumVar(col, 0, +IloInfinity, ILOFLOAT, nome));
			++numRestrGamma;
			col.end();
		}
	}
}

bool NoArvore::alcancaRaiz_no(char opSub){
	int aux;
	int* pricingSat = new int[nSatellites+1];
	bool alcancouRaiz, bcpCap, bcpComb, bcpMStar, bcpCycle, bcpRoute;

	do
	{
		memset(pricingSat, 0, (nSatellites+1)*sizeof(int));
		do
		{
			alcancouRaiz = true;
			for ( int s = 1; s <= nSatellites; ++s )
			{
				if ( pricingSat[s] == 0 )
				{
					mCplex->solveMaster();
					double valSub = mCplex->getValorSubtrair( G, s );
					atualizaCustosDuais( s );

					if ( opSub == 'G' )
					{
						GRASP grasp( G, 0.6 );
						aux = grasp.run( G->getNumCustomers()*20, valSub );
						if ( aux > 0 )
						{
							alcancouRaiz = false;
							for ( int h = 0; h < aux; ++h ) inserirColuna_no( grasp.getRotaConvertida( h, s ), s );
							for ( int i = 1; i <= nSatellites; ++i ) pricingSat[i] = ( pricingSat[i] == 1 ) ? -1 : pricingSat[i];
						}
						else
						{
							ModeloBC bc(G, s, valSub);
							bc.calculaCaminhoElementar(G);
							aux = bc.rotasNegativas.size();
							if ( aux > 0 )
							{
								alcancouRaiz = false;
								for ( int i = 0; i < aux; ++i ) inserirColuna_no( bc.rotasNegativas[i], s );
								for ( int i = 1; i <= nSatellites; ++i ) pricingSat[i] = ( pricingSat[i] == 1 ) ? -1 : pricingSat[i];
							}
							else
							{
								pricingSat[s] = 1;
							}
						}
					}
					else if ( opSub == 'N' )
					{
						NaoElementar ne( G );
						ne.calculaMatrizCiclo2( G );
						vector <Rota*> rrr = ne.getRotaCustoMinimoCiclo2(G, s, valSub);
						if ( rrr.size() > 0 )
						{
							alcancouRaiz = false;
							inserirColuna_no( rrr[0], s );
							if ( rrr.size() > 1 ) inserirColuna_no( rrr[1], s );
							for ( int i = 1; i <= nSatellites; ++i ) pricingSat[i] = ( pricingSat[i] == 1 ) ? -1 : pricingSat[i];
						}
						else
						{
							pricingSat[s] = 1;
						}
					}
				}
			}

			if ( alcancouRaiz )
			{
				aux = 0;
				for ( int i = 1; i <= nSatellites; ++i ) aux += pricingSat[i];
				if ( aux < nSatellites ) memset(pricingSat, 0, (nSatellites+1)*sizeof(int));
				else break;
			}

		}
		while( true );

		if ( ( this->separarCortes == 1 ) || ( ( this->separarCortes == 2 ) && ( totalNosCriados % 300 == 0 ) ) )
		{
			bcpCap = ( mCplex->sepCuts[0] == '2' ) ? mCplex->addCapacityCuts(G, this->qRotas_no, this->ptrRotas_no) : false;
			if ( ( mCplex->sepCuts[0] != '0' ) && ( ( mCplex->sepCuts[1] != '0' ) || ( mCplex->sepCuts[2] != '0' ) || ( mCplex->sepCuts[3] != '0' ) || ( mCplex->sepCuts[4] != '0' ) ) ) mCplex->solveMaster();
			bcpComb = ( mCplex->sepCuts[1] == '2' ) ? mCplex->addCombCuts(G, this->qRotas_no, this->ptrRotas_no) : false;
			if ( ( mCplex->sepCuts[1] != '0' ) && ( ( mCplex->sepCuts[2] != '0' ) || ( mCplex->sepCuts[3] != '0' ) || ( mCplex->sepCuts[4] != '0' ) ) ) mCplex->solveMaster();
			bcpMStar = ( mCplex->sepCuts[2] == '2' ) ? mCplex->addMultiStarCuts(G, this->qRotas_no, this->ptrRotas_no) : false;
			if ( ( mCplex->sepCuts[2] != '0' ) && ( mCplex->sepCuts[3] != '0' ) || ( mCplex->sepCuts[4] != '0' ) ) mCplex->solveMaster();
			bcpCycle = ( mCplex->sepCuts[3] == '2' ) ? mCplex->addCycleCuts(G, this->qRotas_no, this->a_ir_no, this->ptrRotas_no) : false;
			if ( ( mCplex->sepCuts[3] != '0' ) && ( mCplex->sepCuts[4] != '0' ) ) mCplex->solveMaster();
			bcpRoute = ( mCplex->sepCuts[4] != '0' ) ? mCplex->addRouteCuts(G, this->qRotas_no, this->a_ir_no, this->ptrRotas_no) : false;
			this->separarCortes = 2;
		}
		else
		{
			bcpCap = bcpComb = bcpMStar = bcpCycle = bcpRoute = false;
		}
	}
	while ( bcpCap || bcpComb || bcpMStar || bcpCycle );

	delete [] pricingSat;
	limiteDual = mCplex->solveMaster();
	IloNumArray valoresVarArtif(mCplex->env);
	mCplex->cplex.getValues(valoresVarArtif, mCplex->artificiais);
	int numVarArtificiais = valoresVarArtif.getSize();
	//verifica se o no apresenta uma solucao de LP inviavel (caso ocorra, poda-se o no por inviabilidade)
	for (int i = 0; i < numVarArtificiais; ++i) { if ( valoresVarArtif[i] > 0.00001 ) return false; }

	//apos chegar avaliar o LP deste no, verifica quais sao as variaveis com custo reduzido 0 (possiveis basicas) e marca
	IloNumArray redCosts(mCplex->env);
	for ( int s = 1; s <= nSatellites; ++s )
	{
		mCplex->cplex.getReducedCosts( redCosts, mCplex->gamma_no[s] );
		for ( int r = 0; r < qRotas_no[s]; ++r )
		{
			if ( redCosts[r] < 0.001 ) //a variavel associada esta na base ou tem chances potenciais (custo reduzido 0)
			{
				a_ir_no[s][r][0] = 1;
			}
		}
	}
	redCosts.end();
	if ( ( altura % 8 ) == 0 )
	{
		mCplex->processaBase( 0, qRotas_no, a_ir_no, ptrRotas_no ); //a cada altura, executa a rotina para eliminar rotas
		mCplex->solveMaster();
	}

	return true;
}

void NoArvore::atualizaCustosDuais( int s ){
	double custoDual;
	IloNumArray muDuals(mCplex->env), chiDuals(mCplex->env), zetaDuals(mCplex->env);
	mCplex->cplex.getDuals(muDuals, mCplex->constraints4);
	double piDual = mCplex->cplex.getDual(mCplex->constraints5[s-1]);

	int numCycleCuts = mCplex->cycleCutsLHS[s].size();
	int numRouteCuts = mCplex->routeCutsLHS[s].size();
	if ( numCycleCuts > 0 ) mCplex->cplex.getDuals(chiDuals, mCplex->cycleCuts[s]);
	if ( numRouteCuts > 0 ) mCplex->cplex.getDuals(zetaDuals, mCplex->routeCuts[s]);

	for ( int i = nSatellites+1; i < nVertices; ++i )
	{
		custoDual = -muDuals[i-nSatellites-1] - G->getCargaVertice(i)*piDual;

		for ( int j = 0; j < numCycleCuts; ++j ) //custos duais nos vertices vindos dos cycleCuts
		{
			for ( int k = 2; k <= mCplex->cycleCutsLHS[s][j][0]; ++k )
			{
				if ( mCplex->cycleCutsLHS[s][j][k] == i ) custoDual += chiDuals[j];
			}
		}

		for ( int j = 0; j < numRouteCuts; ++j ) //custos duais nos vertices vindos dos routeFeasibilityCuts
		{
			for ( int k = 1; k <= mCplex->routeCutsLHS[s][j][0]; ++k )
			{
				if ( mCplex->routeCutsLHS[s][j][k] == i ) custoDual += zetaDuals[j];
			}
		}

		G->setCustoVerticeDual(i, custoDual);
	}
	G->setCustoArestasDual(s);
	muDuals.end();

	//custos duais associados aos cortes capacity, comb e multistar
	int totalCombCuts = mCplex->combCutsLHS.size();
	int totalCapacityCuts = mCplex->capacityCutsLHS.size();
	int totalMultiStarCuts = mCplex->multiStarCutsLHS.size();
	if ( ( totalCombCuts > 0 ) || ( totalCapacityCuts > 0 ) || ( totalMultiStarCuts > 0 ) )
	{
		IloNumArray psiDuals(mCplex->env), roDuals(mCplex->env), omegaDuals(mCplex->env);
		if ( totalCombCuts > 0 ) mCplex->cplex.getDuals(psiDuals, mCplex->combCuts);
		if ( totalCapacityCuts > 0 ) mCplex->cplex.getDuals(roDuals, mCplex->capacityCuts);
		if ( totalMultiStarCuts > 0 ) mCplex->cplex.getDuals(omegaDuals, mCplex->multiStarCuts);

		for ( int i = 0; i < nVertices; ++i )
		{
			for ( int j = i+1; j < nVertices; ++j )
			{
				custoDual = 0;
				for ( int c = 0; c < totalCombCuts; ++c )
				{
					if ( mCplex->combCutsLHS[c][i][j] > 0 ) custoDual += mCplex->combCutsLHS[c][i][j]*psiDuals[c];
				}
				for ( int c = 0; c < totalCapacityCuts; ++c )
				{
					if ( mCplex->capacityCutsLHS[c][i][j] != 0 ) custoDual += mCplex->capacityCutsLHS[c][i][j]*roDuals[c];
				}
				for ( int m = 0; m < totalMultiStarCuts; ++m )
				{
					if ( mCplex->multiStarCutsLHS[m][i][j] != 0 ) custoDual += mCplex->multiStarCutsLHS[m][i][j]*omegaDuals[m];
				}
				if ( ( custoDual > 0.00001 ) || ( custoDual < -0.00001 ) )
				{
					if ( i == 0 )
					{
						G->setCustoArestaDual(0, j-nSatellites, -custoDual);
						G->setCustoArestaDual(j-nSatellites, 0, -custoDual);
					}
					else
					{
						G->setCustoArestaDual(i-nSatellites, j-nSatellites, -custoDual);
						G->setCustoArestaDual(j-nSatellites, i-nSatellites, -custoDual);
					}
				}
			}
			if ( i == 0 ) i = nSatellites; //para 'pular' as linhas dos satellites, que tem todos os valores 0
		}
		psiDuals.end(); roDuals.end(); omegaDuals.end();
	}

	//custos duais associados aos cycleCuts
	for ( int i = 0; i < numCycleCuts; ++i )
	{
		for ( int j = 1; j < mCplex->cycleCutsLHS[s][i][0]; ++j )
		{
			G->setCustoArestaDual(mCplex->cycleCutsLHS[s][i][j]-nSatellites, mCplex->cycleCutsLHS[s][i][j+1]-nSatellites, -chiDuals[i]);
			G->setCustoArestaDual(mCplex->cycleCutsLHS[s][i][j+1]-nSatellites, mCplex->cycleCutsLHS[s][i][j]-nSatellites, -chiDuals[i]);
		}
		for ( int j = 3; j <= mCplex->cycleCutsLHS[s][i][0]; ++j )
		{
			for ( int k = 1; k <= (j-2); ++k )
			{
				G->setCustoArestaDual(mCplex->cycleCutsLHS[s][i][j]-nSatellites, mCplex->cycleCutsLHS[s][i][k]-nSatellites, -chiDuals[i]);
				G->setCustoArestaDual(mCplex->cycleCutsLHS[s][i][k]-nSatellites, mCplex->cycleCutsLHS[s][i][j]-nSatellites, -chiDuals[i]);
			}
		}
	}

	//custos duais associados aos routeFeasibilityCuts
	for ( int i = 0; i < numRouteCuts; ++i )
	{
		for ( int j = 1; j <= mCplex->routeCutsLHS[s][i][0]; ++j )
		{
			G->setCustoArestaDual(0, mCplex->routeCutsLHS[s][i][j]-nSatellites, -zetaDuals[i]);
			G->setCustoArestaDual(mCplex->routeCutsLHS[s][i][j]-nSatellites, 0, -zetaDuals[i]);

			for ( int k = j+1; k <= mCplex->routeCutsLHS[s][i][0]; ++k )
			{
				G->setCustoArestaDual(mCplex->routeCutsLHS[s][i][j]-nSatellites, mCplex->routeCutsLHS[s][i][k]-nSatellites, -zetaDuals[i]);
				G->setCustoArestaDual(mCplex->routeCutsLHS[s][i][k]-nSatellites, mCplex->routeCutsLHS[s][i][j]-nSatellites, -zetaDuals[i]);
			}
		}
	}
	chiDuals.end(); zetaDuals.end();

	//insere no grafo os custos duais associados aos arcos que foram fixados para este no
	int auxI, auxJ, numRestrGamma = 0;
	for (int i = 0; i < arcosBranching.size(); ++i)
	{
		if ( arcosBranching[i][2] == s )
		{
			auxI = ( arcosBranching[i][0] <= nSatellites ) ? 0 : ( arcosBranching[i][0]-nSatellites );
			auxJ = ( arcosBranching[i][1] <= nSatellites ) ? 0 : ( arcosBranching[i][1]-nSatellites );
			G->setCustoArestaDual(auxI, auxJ, -mCplex->cplex.getDual(mCplex->constraintsArvoreL2[numRestrGamma]));
		}
		if ( arcosBranching[i][2] > 0 ) ++numRestrGamma;
	}
}

void NoArvore::inserirColuna_no(Rota* r, int s){
	//Armazena a rota na matriz de rotas, que sera usada posteriormente
	short int* rota = new short int[nVertices];
	memset(rota, 0, nVertices*sizeof(short int));
	a_ir_no[s].push_back(rota);
	ptrRotas_no[s].push_back(r);
	r->incrNumApontadores();

	vector <int> vertRota = r->getVertices();
	int numVertAtualizar = vertRota.size()-2;
	for( int j = 1; j <= numVertAtualizar; ++j )
	{
		++a_ir_no[s][qRotas_no[s]][vertRota[j]]; //+1 apenas nas posicoes que o vertice esteja na rota
	}

	//CRIA-SE UMA COLUNA ASSOCIADA A ROTA PASSADA COMO PARAMETRO E INSERE NA VARIAVEL ASSOCIADA AO SATELLITE S
	IloNumColumn col = mCplex->objCost( r->getCusto() );

	col += mCplex->constraint2(1);
	col += mCplex->constraints3[s-1](1);

	int coeficiente5 = 0;
	for ( int i = 0; i < mCplex->nCustomers; ++i )
	{
		if ( a_ir_no[s][qRotas_no[s]][nSatellites+i+1] > 0 )
		{
			col += mCplex->constraints4[i]( a_ir_no[s][qRotas_no[s]][nSatellites+i+1] );
			coeficiente5 += ( a_ir_no[s][qRotas_no[s]][nSatellites+i+1] * G->getCargaVertice( nSatellites+i+1 ) );

			if ( a_ir_no[s][qRotas_no[s]][nSatellites+i+1] > 1 ) r->ciclo = true;
		}
	}
	col += mCplex->constraints5[s-1](coeficiente5);

	//insere a coluna nas restricoes do corte de capacidade
	float indexCol;
	int numCapacityCuts = mCplex->capacityCutsLHS.size();
	for ( int i = 0; i < numCapacityCuts; ++i )
	{
		indexCol = 0;
		for ( int v = 1; v < numVertAtualizar; ++v )
		{
			if ( mCplex->capacityCutsLHS[i][vertRota[v]][vertRota[v+1]] != 0 ) indexCol += mCplex->capacityCutsLHS[i][vertRota[v]][vertRota[v+1]];
			if ( mCplex->capacityCutsLHS[i][vertRota[v+1]][vertRota[v]] != 0 ) indexCol += mCplex->capacityCutsLHS[i][vertRota[v+1]][vertRota[v]];
		}
		if ( mCplex->capacityCutsLHS[i][0][vertRota[1]] != 0 ) indexCol += mCplex->capacityCutsLHS[i][0][vertRota[1]];
		if ( mCplex->capacityCutsLHS[i][0][vertRota[numVertAtualizar]] != 0 ) indexCol += mCplex->capacityCutsLHS[i][0][vertRota[numVertAtualizar]];

		if ( ( indexCol > 0.00001 ) || ( indexCol < -0.00001 ) ) col += mCplex->capacityCuts[i](indexCol);
	}

	//insere a coluna tambem nas restricoes dos cortes comb
	int numCombCuts = mCplex->combCutsLHS.size();
	for ( int i = 0; i < numCombCuts; ++i )
	{
		indexCol = 0;
		for ( int v = 1; v < numVertAtualizar; ++v )
		{
			if ( mCplex->combCutsLHS[i][vertRota[v]][vertRota[v+1]] > 0 ) indexCol += mCplex->combCutsLHS[i][vertRota[v]][vertRota[v+1]];
			if ( mCplex->combCutsLHS[i][vertRota[v+1]][vertRota[v]] > 0 ) indexCol += mCplex->combCutsLHS[i][vertRota[v+1]][vertRota[v]];
		}
		if ( mCplex->combCutsLHS[i][0][vertRota[1]] > 0 ) indexCol += mCplex->combCutsLHS[i][0][vertRota[1]];
		if ( mCplex->combCutsLHS[i][0][vertRota[numVertAtualizar]] > 0 ) indexCol += mCplex->combCutsLHS[i][0][vertRota[numVertAtualizar]];

		if ( indexCol > 0.00001 ) col += mCplex->combCuts[i](indexCol);
	}

	//insere a coluna nas restricoes dos cortes multiStar
	int numMultiStarCuts = mCplex->multiStarCutsLHS.size();
	for ( int i = 0; i < numMultiStarCuts; ++i )
	{
		indexCol = 0;
		for ( int v = 1; v < numVertAtualizar; ++v )
		{
			if ( mCplex->multiStarCutsLHS[i][vertRota[v]][vertRota[v+1]] != 0 ) indexCol += mCplex->multiStarCutsLHS[i][vertRota[v]][vertRota[v+1]];
			if ( mCplex->multiStarCutsLHS[i][vertRota[v+1]][vertRota[v]] != 0 ) indexCol += mCplex->multiStarCutsLHS[i][vertRota[v+1]][vertRota[v]];
		}
		if ( mCplex->multiStarCutsLHS[i][0][vertRota[1]] != 0 ) indexCol += mCplex->multiStarCutsLHS[i][0][vertRota[1]];
		if ( mCplex->multiStarCutsLHS[i][0][vertRota[numVertAtualizar]] != 0 ) indexCol += mCplex->multiStarCutsLHS[i][0][vertRota[numVertAtualizar]];

		if ( ( indexCol > 0.00001 ) || ( indexCol < -0.00001 ) ) col += mCplex->multiStarCuts[i](indexCol);
	}

	//insere a coluna nas restricoes dos cortes outfork (caso existam)
	int numCycleCuts = mCplex->cycleCutsLHS[s].size();
	for ( int i = 0; i < numCycleCuts; ++i )
	{
		indexCol = 0;
		for ( int v = 1; v < numVertAtualizar; ++v )
		{
			for ( int y = 1; y < mCplex->cycleCutsLHS[s][i][0]; ++y )
			{
				if ( ( ( vertRota[v] == mCplex->cycleCutsLHS[s][i][y] ) && ( vertRota[v+1] == mCplex->cycleCutsLHS[s][i][y+1] ) ) ||
					 ( ( vertRota[v] == mCplex->cycleCutsLHS[s][i][y+1] ) && ( vertRota[v+1] == mCplex->cycleCutsLHS[s][i][y] ) ) ) ++indexCol;
			}
			for ( int y = 3; y <= mCplex->cycleCutsLHS[s][i][0]; ++y )
			{
				for ( int z = 1; z <= (y-2); ++z )
				{
					if ( ( ( vertRota[v] == mCplex->cycleCutsLHS[s][i][y] ) && ( vertRota[v+1] == mCplex->cycleCutsLHS[s][i][z] ) ) ||
						( ( vertRota[v] == mCplex->cycleCutsLHS[s][i][z] ) && ( vertRota[v+1] == mCplex->cycleCutsLHS[s][i][y] ) ) ) ++indexCol;
				}
			}
		}
		for ( int y = 2; y <= mCplex->cycleCutsLHS[s][i][0]; ++y )
		{
			if ( a_ir_no[s][qRotas_no[s]][mCplex->cycleCutsLHS[s][i][y]] > 0 ) indexCol -= a_ir_no[s][qRotas_no[s]][mCplex->cycleCutsLHS[s][i][y]];
		}

		if ( ( indexCol > 0.00001 ) || ( indexCol < -0.00001 ) ) col += mCplex->cycleCuts[s][i](indexCol);
	}

	//insere a coluna nas restricoes dos routeFeasibilityCuts (caso existam)
	int numRouteCuts = mCplex->routeCutsLHS[s].size();
	for ( int i = 0; i < numRouteCuts; ++i )
	{
		indexCol = 0;
		for ( int v = 1; v < numVertAtualizar; ++v )
		{
			for ( int y = 1; y < mCplex->routeCutsLHS[s][i][0]; ++y )
			{
				for ( int z = y+1; z <= mCplex->routeCutsLHS[s][i][0]; ++z )
				{
					if ( ( ( vertRota[v] == mCplex->routeCutsLHS[s][i][y] ) && ( vertRota[v+1] == mCplex->routeCutsLHS[s][i][z] ) ) ||
						( ( vertRota[v] == mCplex->routeCutsLHS[s][i][z] ) && ( vertRota[v+1] == mCplex->routeCutsLHS[s][i][y] ) ) ) ++indexCol;
				}
			}
		}
		for ( int y = 1; y <= mCplex->routeCutsLHS[s][i][0]; ++y )
		{
			if ( vertRota[1] == mCplex->routeCutsLHS[s][i][y] ) ++indexCol;
			if ( vertRota[numVertAtualizar] == mCplex->routeCutsLHS[s][i][y] ) ++indexCol;
			if ( a_ir_no[s][qRotas_no[s]][mCplex->routeCutsLHS[s][i][y]] > 0 ) indexCol -= a_ir_no[s][qRotas_no[s]][mCplex->routeCutsLHS[s][i][y]];
		}

		if ( ( indexCol > 0.00001 ) || ( indexCol < -0.00001 ) ) col += mCplex->routeCuts[s][i](indexCol);
	}

	int numVisitasArco, numRestrVeic = 0, numRestrGamma = 0, arcosBranchingSize = arcosBranching.size();
	for ( int i = 0; i < arcosBranchingSize; ++i )
	{
		//inclui na coluna os indices associados as restricoes de branching no numero de veiculos que saem do satellite s (caso existam)
		if ( arcosBranching[i][2] == -s ) col += mCplex->constraintsArvoreVeic[s][numRestrVeic++](1);

		//inclui na coluna os indices associados as restricoes de branching em arcos (caso existam)
		if ( arcosBranching[i][2] == s )
		{
			numVisitasArco = 0;
			for (int j = 0; j <= numVertAtualizar; ++j)
			{
				if ( ( vertRota[j] == arcosBranching[i][0] ) && ( ( vertRota[j+1] == arcosBranching[i][1] ) ) )
				{
					++numVisitasArco;
				}
			}
			if (numVisitasArco > 0)	col += mCplex->constraintsArvoreL2[numRestrGamma](numVisitasArco);
		}
		if ( arcosBranching[i][2] > 0 ) ++numRestrGamma;
	}

	char nome[20];
	sprintf(nome, "gNo_%d_%d", s, qRotas_no[s]);
	mCplex->gamma_no[s].add(IloNumVar(col, 0, 1, ILOFLOAT, nome));
	++qRotas_no[s];
	col.end();
}

bool NoArvore::defineVariavelBranching(int& inicioArco, int& fimArco, int &satellite, int &rhs){
	int valorInt, tamRota;
	vector<int> verticesRota;
	double fracao, numVeicSol, violacao = 0, violacaoS = 0;
	IloNumArray lambdaValues ( mCplex->env );

	//a prioridade 1 eh para fazer branching nas variaveis lambda
	satellite = 0;
	mCplex->cplex.getValues( lambdaValues, mCplex->lambda );
	for ( int r = 0; r < mCplex->quantRotasL1; ++r )
	{
		if ( lambdaValues[r] > 0.00001 )
		{
			valorInt = lambdaValues[r];
			fracao = lambdaValues[r] - valorInt;
			if ( min( fracao, ( 1-fracao ) ) > ( violacao + 0.0001 ) )
			{
				violacao = min( fracao, ( 1-fracao ) );
				inicioArco = r; rhs = ( valorInt+1 );
			}
		}
	}
	lambdaValues.end();
	if ( violacao > 0.00001 ) return true;

	//flag para definir q este no tem lambdas inteiros (os cortes devem ser separados para este no)
	if ( separarCortes == 0 ) separarCortes = 1;

	//a prioridade 2 eh para fazer branching na quantidade de veiculos que partem de um satellite
	IloArray < IloNumArray > gammaValues ( mCplex->env, nSatellites+1 );
	IloArray < IloNumArray > gammaValues_no ( mCplex->env, nSatellites+1 );
	for ( int s = 1; s <= nSatellites; ++s )
	{
		numVeicSol = 0;
		gammaValues[s] = IloNumArray( mCplex->env );
		mCplex->cplex.getValues(gammaValues[s], mCplex->gamma[s]);
		for ( int r = 0; r < mCplex->qRotasL2[s]; ++r )
		{
			if ( gammaValues[s][r] > 0.00001 ) numVeicSol += gammaValues[s][r];
		}

		if ( qRotas_no[s] > 0 )
		{
			gammaValues_no[s] = IloNumArray( mCplex->env );
			mCplex->cplex.getValues(gammaValues_no[s], mCplex->gamma_no[s]);
			for ( int r = 0; r < qRotas_no[s]; ++r )
			{
				if ( gammaValues_no[s][r] > 0.00001 ) numVeicSol += gammaValues_no[s][r];
			}
		}

		//verifica qual satellite tem a maior violacao (caso algum viole a integralidade)
		valorInt = numVeicSol;
		fracao = numVeicSol - valorInt;
		if ( min( fracao, ( 1-fracao ) ) > ( violacao + 0.0001 ) )
		{
			violacao = min( fracao, ( 1-fracao ) );
			satellite = -s; rhs = ( valorInt+1 );
		}
	}
	if ( violacao > 0.00001 )
	{
		for ( int s = 1; s <= nSatellites; ++s )
		{
			gammaValues[s].end();
			gammaValues_no[s].end();
		}
		gammaValues.end();
		gammaValues_no.end();
		return true;
	}

	//a prioridade 3 eh para fazer branching em arcos que partem de um satellite
	//a prioridade 4 eh para fazer branching em quaisquer outros arcos (i, j)[s]
	//Cada celula da matriz representara um arco (i,j) que sera preenchido em funcao da violacao
	//Aquele arco (i,j) que tiver o maior valor de violacao sera fixado no branching de variaveis
	double** matrizViolacao = new double*[nVertices];
	for (int i = 0; i < nVertices; ++i) matrizViolacao[i] = new double[nVertices];

	//A cada iteracao procura-se pela violacao dos arcos de um satellite e armazena.
	//Aquele arco com a maior violacao entre todos os satellites sera selecionado para o branching
	for (int s = 1; s <= nSatellites; ++s)
	{
		for (int i = 0; i < nVertices; ++i) memset(matrizViolacao[i], 0, nVertices*sizeof(double));

		//primeiro busca por variaveis fracionarias entre as variaveis de rotas L2 obtidas na raiz
		for (int r = 0; r < mCplex->qRotasL2[s]; ++r)
		{
			if ( ( gammaValues[s][r] > 0.00001 ) && ( gammaValues[s][r] < 0.99999 ) ) //gamma[s][r] da raiz eh fracionaria
			{
				verticesRota = mCplex->ptrRotasL2[s][r]->getVertices();
				tamRota = verticesRota.size()-1;
				for (int i = 0; i < tamRota; ++i)
				{
					matrizViolacao[verticesRota[i]][verticesRota[i+1]] += gammaValues[s][r];
				}
			}
		}

		//depois computa quais variaveis fracionarias do no
		if ( qRotas_no[s] > 0 )
		{
			for (int r = 0; r < qRotas_no[s]; ++r)
			{
				if ( ( gammaValues_no[s][r] > 0.00001 ) && ( gammaValues_no[s][r] < 0.99999 ) )
				{
					verticesRota = ptrRotas_no[s][r]->getVertices();
					tamRota = verticesRota.size()-1;
					for (int i = 0; i < tamRota; ++i)
					{
						matrizViolacao[verticesRota[i]][verticesRota[i+1]] += gammaValues_no[s][r];
					}
				}
			}
		}

		//obtem o valor do(s) satellite(s) mais violado(s) quanto ao numero de veiculos que partem dele
		for ( int i = 1; i < nVertices; ++i )
		{
			if ( min( matrizViolacao[s][i], ( 1-matrizViolacao[s][i] ) ) > ( violacaoS + 0.00001 ) )
			{
				violacaoS = min( matrizViolacao[s][i], ( 1-matrizViolacao[s][i] ) );
				inicioArco = s; fimArco = i; satellite = s;
			}
			else if ( min( matrizViolacao[i][s], ( 1-matrizViolacao[i][s] ) ) > ( violacaoS + 0.00001 ) )
			{
				violacaoS = min( matrizViolacao[i][s], ( 1-matrizViolacao[i][s] ) );
				inicioArco = i; fimArco = s; satellite = s;
			}
		}

		if ( violacaoS < 0.00001 )
		{
			//obtem o VALOR do(s) arco(s) mais violado(s)
			for ( int i = nSatellites+1; i < nVertices; ++i )
			{
				for ( int j = nSatellites+1; j < nVertices; ++j )
				{
					if ( min( matrizViolacao[i][j], ( 1-matrizViolacao[i][j] ) ) > ( violacao + 0.0001 ) )
					{
						violacao = min( matrizViolacao[i][j], ( 1-matrizViolacao[i][j] ) );
						inicioArco = i; fimArco = j; satellite = s;
					}
				}
			}
		}
	}

	//libera memoria
	for ( int i = 0; i < nVertices; i++ ) delete [] matrizViolacao[i];
	for ( int s = 1; s <= nSatellites; ++s )
	{
		gammaValues[s].end();
		gammaValues_no[s].end();
	}
	gammaValues.end();
	gammaValues_no.end();
	delete [] matrizViolacao;

	if ( ( violacaoS > 0.00001 ) || ( violacao > 0.00001 ) ) return true;//existem arcos para executar branching
	else return false; //solucao inteira, nao executa-se branching e poda por integralidade
}


double NoArvore::runCplexInt()
{
	for ( int s = 1; s <= nSatellites; ++s )
	{
		mCplex->gamma_noCGH[s].endElements();
		for ( int r = 0; r < qRotas_no[s]; ++r )
		{
			IloNumColumn colCGH = mCplex->objCostCGH( ptrRotas_no[s][r]->getCusto() );
			colCGH += mCplex->constraint2CGH(1);
			colCGH += mCplex->constraints3CGH[s-1](1);

			int coeficiente5 = 0;
			for ( int i = 0; i < mCplex->nCustomers; ++i )
			{
				if ( a_ir_no[s][r][nSatellites+i+1] > 0 )
				{
					colCGH += mCplex->constraints4CGH[i]( a_ir_no[s][r][nSatellites+i+1] );
					coeficiente5 += ( a_ir_no[s][r][nSatellites+i+1] * G->getCargaVertice( nSatellites+i+1 ) );
				}
			}
			colCGH += mCplex->constraints5CGH[s-1](coeficiente5);
			mCplex->gamma_noCGH[s].add(IloIntVar(colCGH, 0, 1));
			colCGH.end();
		}
	}

	mCplex->cplexCGH.setParam(IloCplex::Threads, 1);
	mCplex->cplexCGH.setParam(IloCplex::TiLim, 300);
	mCplex->cplexCGH.solve();

	double limitePrimal;
	if ((mCplex->cplexCGH.getStatus() == IloAlgorithm::Feasible) || (mCplex->cplexCGH.getStatus() == IloAlgorithm::Optimal))
	{
		limitePrimal = mCplex->cplexCGH.getObjValue();
	}
	else
	{
		limitePrimal = MAIS_INFINITO;
	}

	return limitePrimal;
}


void NoArvore::imprimir() {
	printf("indice = %d\n", indiceNo);
	printf("valor = %f\n", limiteDual);
	for (int s = 1; s <= nSatellites; ++s)
	{
		printf("qRotas_no[%d] = %d\n", s, qRotas_no[s]);
		for (int r = 0; r < qRotas_no[s]; ++r)
		{
			printf(" ");
			ptrRotas_no[s][r]->imprimir();
		}
	}
	printf("Arcos com restricao:\n");
	for (int i = 0; i < arcosBranching.size(); ++i)
	{
		printf("  (%d, %d)[%d] = %d\n", arcosBranching[i][0], arcosBranching[i][1], arcosBranching[i][2], arcosBranching[i][3]);
	}
	printf("\n");
}

void NoArvore::setProx(ptrNo p){
	prox = p;
}

ptrNo NoArvore::getProx(){
	return prox;
}

double NoArvore::getLimiteDual(){
	return limiteDual;
}

int NoArvore::getIndiceNo(){
	return indiceNo;
}

int NoArvore::getTotalNosCriados(){
	return totalNosCriados;
}

int NoArvore::getTotalNosAtivos(){
	return totalNosAtivos;
}

void NoArvore::addNoProcessado(ptrNo no)
{
	//o procedimento de insercao basicamente verifica se o pai do no a ser inserido esta na lista
	//caso esteja, o pai do no eh trocado pelo filho. Caso contrario, o no eh inserido no inicio da lista
	ptrNo aux, it = listaNosProcessados;

	if ( listaNosProcessados == NULL ) //se a lista estiver vazia, insere o no na primeira posicao
	{
		no->prox = NULL;
		listaNosProcessados = no;
		return;
	}

	if ( no->pai == listaNosProcessados ) //se o primeiro no da lista for pai, o coloca na lista de pais
	{
		aux = listaNosProcessados; 
		no->prox = listaNosProcessados->prox;
		listaNosProcessados = no;
		delete aux;
		return;
	}

	while ( it->prox != NULL ) //percorre a lista verificando se o no pai do no a ser inserido esta na lista
	{
		if ( no->pai == it->prox )
		{
			aux = it->prox;
			no->prox = aux->prox;
			it->prox = no;
			delete aux;
			return;
			
		}
		else if ( it->prox->pai == no )
		{
			delete no;
			return;
		}
		it = it->prox;
	}

	no->prox = listaNosProcessados;
	listaNosProcessados = no;
}

ptrNo NoArvore::getNoProcessado()
{
	if ( listaNosProcessados == NULL ) return NULL;
	ptrNo it = listaNosProcessados;
	listaNosProcessados = it->prox;
	return it;
}

