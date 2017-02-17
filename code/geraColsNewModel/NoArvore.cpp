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
	mCplex = modCplex;
	totalNosAtivos = 0;
	totalNosCriados = 0;
	nSatellites = g->getNumSatellites();
	nVertices = g->getNumVertices();
	listaNosProcessados = NULL;
	limiteDual = lD;
	indiceNo = 0;
	int i, j, s, rhs;

	qRotas_no = new int[nSatellites+1];
	ptrRotas_no = new vector <Rota*>[nSatellites+1];
	a_ir_no = new vector <short int*>[nSatellites+1];
	for ( int s = 0; s <= nSatellites; ++s ) qRotas_no[s] = 0;

	//instancia os objetos do cplex para manipular as restricoes e variaveis do branching
	mCplex->gamma_no = IloArray < IloNumVarArray > ( mCplex->env, nSatellites+1 );
	for ( int s = 0; s <= nSatellites; ++s ) mCplex->gamma_no[s] = IloNumVarArray( mCplex->env );
	mCplex->artificiais = IloNumVarArray( mCplex->env );
	mCplex->constraintsArvoreL1 = IloRangeArray( mCplex->env );
	mCplex->constraintsArvoreL2 = IloRangeArray( mCplex->env );

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
	limiteDual = pai->limiteDual;
	for (int s = 0; s <= nSatellites; ++s)
	{
		qRotas_no[s] = pai->qRotas_no[s];
		for (int r = 0; r < qRotas_no[s]; ++r)
		{
			tmp = new short int[nVertices];
			for (int i = nSatellites+1; i < nVertices; ++i)
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
	for ( int s = 1; s <= nSatellites; ++s ) mCplex->gamma_no[s].endElements();

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
	int coef, count, numVertAtualizar, numEdgesCut, numCuts = mCplex->capacityCutsLHS.size();

	//para cada rota partindo de um satellite inclue uma rota no ModeloCplex
	for ( int s = 1; s <= nSatellites; ++s )
	{
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

			//Insere a variavel do no tambem nas restricoes de corte
			vertRota = ptrRotas_no[s][r]->getVertices();
			numVertAtualizar = vertRota.size()-2;
			for ( int i = 0; i < numCuts; ++i )
			{
				indexCol = 0;
				numEdgesCut = mCplex->capacityCutsLHS[i].size();
				if ( mCplex->capacityCutsLHS[i][0] >= 0 )
				{
					for ( int j = 0; j < numEdgesCut; j+=2 )
					{
						for ( int vR = 1; vR < numVertAtualizar; ++vR )
						{
							if ( ( ( vertRota[vR] == mCplex->capacityCutsLHS[i][j] ) && ( vertRota[vR+1] == mCplex->capacityCutsLHS[i][j+1] ) ) || ( ( vertRota[vR] == mCplex->capacityCutsLHS[i][j+1] ) && ( vertRota[vR+1] == mCplex->capacityCutsLHS[i][j] ) ) )
							{
								++indexCol;
							}
						}
					}
					if ( indexCol > 0 ) col += mCplex->capacityCuts[i](indexCol);
				}
				else
				{
					int j, tam = -mCplex->capacityCutsLHS[i][0];
					for ( j = 1; j <= tam; ++j )
					{
						if ( vertRota[1] == mCplex->capacityCutsLHS[i][j] ) indexCol += 0.5;
						if ( vertRota[numVertAtualizar] == mCplex->capacityCutsLHS[i][j] ) indexCol += 0.5;
					}
					while ( mCplex->capacityCutsLHS[i][j] >= 0 )
					{
						for ( int vR = 1; vR < numVertAtualizar; ++vR )
						{
							if ( ( ( vertRota[vR] == mCplex->capacityCutsLHS[i][j] ) && ( vertRota[vR+1] == mCplex->capacityCutsLHS[i][j+1] ) ) || ( ( vertRota[vR] == mCplex->capacityCutsLHS[i][j+1] ) && ( vertRota[vR+1] == mCplex->capacityCutsLHS[i][j] ) ) )
							{
								++indexCol;
							}
						}
						j += 2;
					}
					while ( j < numEdgesCut )
					{
						if ( vertRota[1] == -mCplex->capacityCutsLHS[i][j] ) indexCol -= 0.5;
						if ( vertRota[numVertAtualizar] == -mCplex->capacityCutsLHS[i][j] ) indexCol -= 0.5;
						++j;
					}
					if ( ( indexCol > 0.0001 ) || ( indexCol < 0.0001 ) ) col += mCplex->capacityCuts[i](indexCol);
				}
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
	int sFixo, vFixoI, vFixoJ, quantR, quantV, coefColuna, numVisitasArco, numRestGamma = 0, numRestrLambda = 0;

	//primeiro insere as restricoes associadas as variaveis lambda
	int arcosBranchingSize = arcosBranching.size();
	for ( int i = 0; i < arcosBranchingSize; ++i )
	{
		if ( arcosBranching[i][2] == 0 )
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
	}

	//depois insere as restricoes associadas as variaveis gamma
	for (int i = 0; i < arcosBranchingSize; ++i)
	{
		if ( arcosBranching[i][2] != 0 )
		{
			IloExpr exp = IloExpr(mCplex->env);
			vFixoI = arcosBranching[i][0];
			vFixoJ = arcosBranching[i][1];
			sFixo = arcosBranching[i][2];

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
			mCplex->constraintsArvoreL2[numRestGamma].setName(nome);
			mCplex->model.add(mCplex->constraintsArvoreL2[numRestGamma]);
			exp.end();

			//insere a variavel artificial associada a esta restricao
			sprintf(nome, "artifG%d", numRestGamma);
			IloNumColumn col = mCplex->objCost(ModeloCplex::limitePrimal);
			col += mCplex->constraintsArvoreL2[numRestGamma](coefColuna);
			mCplex->artificiais.add(IloNumVar(col, 0, +IloInfinity, ILOFLOAT, nome));
			++numRestGamma;
			col.end();
		}
	}
}

bool NoArvore::alcancaRaiz_no(char opSub){
	int aux, tmp;
	bool alcancouRaiz, bcp;

	do
	{
		do
		{
			alcancouRaiz = true;
			for ( int s = 1; s <= nSatellites; ++s )
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
					}
					else
					{
						ModeloBC bc(G, s, valSub);
						bc.calculaCaminhoElementar(G);
						aux = bc.rotasNegativas.size();
						if ( aux > 0 )
						{
							alcancouRaiz = false;
							for ( int i = 0; i < aux; ++i )
							{
								inserirColuna_no( bc.rotasNegativas[i], s );
							}
						}
					}
				}
				else if ( opSub == 'B' )
				{
					ModeloBC bc(G, s, valSub);
					bc.calculaCaminhoElementar(G);
					aux = bc.rotasNegativas.size();
					if ( aux > 0 )
					{
						alcancouRaiz = false;
						for ( int i = 0; i < aux; ++i )
						{
							inserirColuna_no( bc.rotasNegativas[i], s );
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
					}
				}
			}
		}
		while(!alcancouRaiz);

		//bcp = ( mCplex->sepCuts == 2 ) ? mCplex->addCapacityCuts(G, qRotas_no, ptrRotas_no) : false;

	}
	while ( bcp );

	limiteDual = mCplex->solveMaster();
	IloNumArray valoresVarArtif(mCplex->env);
	mCplex->cplex.getValues(valoresVarArtif, mCplex->artificiais);
	int numVarArtificiais = valoresVarArtif.getSize();

	for (int i = 0; i < numVarArtificiais; ++i)
	{
		if (valoresVarArtif[i] > 0.00001) return false;
	}
	return true;
}

void NoArvore::atualizaCustosDuais( int s ){
	double custoDual;
	IloNumArray muDuals(mCplex->env);
	mCplex->cplex.getDuals(muDuals, mCplex->constraints4);
	double piDual = mCplex->cplex.getDual(mCplex->constraints5[s-1]);
	for ( int i = 0; i < mCplex->nCustomers; ++i )
	{
		custoDual = -muDuals[i] - G->getCargaVertice(nSatellites+i+1)*piDual;
		G->setCustoVerticeDual(nSatellites+i+1, custoDual);
	}
	G->setCustoArestasDual(s);
	muDuals.end();

	//Insere no grafo os custos associados as restricoes de corte
	int numEdgesCut;
	int totalCapacityCuts = mCplex->capacityCutsLHS.size();
	if ( totalCapacityCuts > 0 )
	{
		IloNumArray roDuals(mCplex->env);
		mCplex->cplex.getDuals(roDuals, mCplex->capacityCuts);
		for ( int i = 0; i < totalCapacityCuts; ++i )
		{
			numEdgesCut = mCplex->capacityCutsLHS[i].size();
			if ( mCplex->capacityCutsLHS[i][0] >= 0 )
			{
				for ( int j = 0; j < numEdgesCut; j+=2 )
				{
					G->setCustoArestaDual(mCplex->capacityCutsLHS[i][j]-nSatellites, mCplex->capacityCutsLHS[i][j+1]-nSatellites, -roDuals[i]);
					G->setCustoArestaDual(mCplex->capacityCutsLHS[i][j+1]-nSatellites, mCplex->capacityCutsLHS[i][j]-nSatellites, -roDuals[i]);
				}
			}
			else
			{
				int j, tam = -mCplex->capacityCutsLHS[i][0];
				for ( j = 1; j <= tam; ++j )
				{
					G->setCustoArestaDual(0, mCplex->capacityCutsLHS[i][j]-nSatellites, -0.5*roDuals[i]);
					G->setCustoArestaDual(mCplex->capacityCutsLHS[i][j]-nSatellites, 0, -0.5*roDuals[i]);
				}
				while ( mCplex->capacityCutsLHS[i][j] >= 0 )
				{
					G->setCustoArestaDual(mCplex->capacityCutsLHS[i][j]-nSatellites, mCplex->capacityCutsLHS[i][j+1]-nSatellites, -roDuals[i]);
					G->setCustoArestaDual(mCplex->capacityCutsLHS[i][j+1]-nSatellites, mCplex->capacityCutsLHS[i][j]-nSatellites, -roDuals[i]);
					j += 2;
				}
				while ( j < numEdgesCut )
				{
					G->setCustoArestaDual(0, (-mCplex->capacityCutsLHS[i][j])-nSatellites, 0.5*roDuals[i]);
					G->setCustoArestaDual((-mCplex->capacityCutsLHS[i][j])-nSatellites, 0, 0.5*roDuals[i]);
					++j;
				}
			}
		}
		roDuals.end();
	}

	int totalCombCuts = mCplex->combCutsLHS.size();
	if ( totalCombCuts > 0 )
	{
		IloNumArray psiDuals(mCplex->env);
		mCplex->cplex.getDuals(psiDuals, mCplex->combCuts);
		for ( int i = 0; i < nVertices; ++i )
		{
			for ( int j = 0; j < nVertices; ++j )
			{
				custoDual = 0;
				for ( int c = 0; c < totalCombCuts; ++c )
				{
					if ( mCplex->combCutsLHS[c][i][j] > 0 ) custoDual += mCplex->combCutsLHS[c][i][j]*psiDuals[c];
				}
				if ( custoDual > 0.00001 )
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
		psiDuals.end();
	}

	//insere no grafo os custos duais associados aos arcos que foram fixados para este no
	int auxI, auxJ, numRestGamma = 0;
	for (int i = 0; i < arcosBranching.size(); ++i)
	{
		if ( arcosBranching[i][2] == s )
		{
			auxI = ( arcosBranching[i][0] <= nSatellites ) ? 0 : ( arcosBranching[i][0]-nSatellites );
			auxJ = ( arcosBranching[i][1] <= nSatellites ) ? 0 : ( arcosBranching[i][1]-nSatellites );
			G->setCustoArestaDual(auxI, auxJ, -mCplex->cplex.getDual(mCplex->constraintsArvoreL2[numRestGamma]));
		}
		if ( arcosBranching[i][2] != 0 ) ++numRestGamma;
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
		}
	}
	col += mCplex->constraints5[s-1](coeficiente5);

	//insere a coluna nas restricoes dos cortes de capacidade
	float indexCol;
	int numEdgesCut, totalCuts = mCplex->capacityCutsLHS.size();
	for ( int i = 0; i < totalCuts; ++i )
	{
		indexCol = 0;
		numEdgesCut = mCplex->capacityCutsLHS[i].size();
		if ( mCplex->capacityCutsLHS[i][0] >= 0 )
		{
			for ( int j = 0; j < numEdgesCut; j+=2 )
			{
				for ( int vR = 1; vR < numVertAtualizar; ++vR )
				{
					if ( ( ( vertRota[vR] == mCplex->capacityCutsLHS[i][j] ) && ( vertRota[vR+1] == mCplex->capacityCutsLHS[i][j+1] ) ) || ( ( vertRota[vR] == mCplex->capacityCutsLHS[i][j+1] ) && ( vertRota[vR+1] == mCplex->capacityCutsLHS[i][j] ) ) )
					{
						++indexCol;
					}
				}
			}
			if ( indexCol > 0 ) col += mCplex->capacityCuts[i](indexCol);
		}
		else
		{
			int j, tam = -mCplex->capacityCutsLHS[i][0];
			for ( j = 1; j <= tam; ++j )
			{
				if ( vertRota[1] == mCplex->capacityCutsLHS[i][j] ) indexCol += 0.5;
				if ( vertRota[numVertAtualizar] == mCplex->capacityCutsLHS[i][j] ) indexCol += 0.5;
			}
			while ( mCplex->capacityCutsLHS[i][j] >= 0 )
			{
				for ( int vR = 1; vR < numVertAtualizar; ++vR )
				{
					if ( ( ( vertRota[vR] == mCplex->capacityCutsLHS[i][j] ) && ( vertRota[vR+1] == mCplex->capacityCutsLHS[i][j+1] ) ) || ( ( vertRota[vR] == mCplex->capacityCutsLHS[i][j+1] ) && ( vertRota[vR+1] == mCplex->capacityCutsLHS[i][j] ) ) )
					{
						++indexCol;
					}
				}
				j += 2;
			}
			while ( j < numEdgesCut )
			{
				if ( vertRota[1] == -mCplex->capacityCutsLHS[i][j] ) indexCol -= 0.5;
				if ( vertRota[numVertAtualizar] == -mCplex->capacityCutsLHS[i][j] ) indexCol -= 0.5;
				++j;
			}
			if ( ( indexCol > 0.0001 ) || ( indexCol < 0.0001 ) ) col += mCplex->capacityCuts[i](indexCol);
		}
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

	//Inclui na coluna os indices associados as restricoes de branching em arcos, caso existam
	int numVisitasArco, numRestGamma = 0, arcosBranchingSize = arcosBranching.size();
	for ( int i = 0; i < arcosBranchingSize; ++i )
	{
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
			if (numVisitasArco > 0)	col += mCplex->constraintsArvoreL2[numRestGamma](numVisitasArco);
		}
		if ( arcosBranching[i][2] != 0 ) ++numRestGamma;
	}

	char nome[20];
	sprintf(nome, "gNo_%d_%d", s, qRotas_no[s]);
	mCplex->gamma_no[s].add(IloNumVar(col, 0, 1, ILOFLOAT, nome));
	++qRotas_no[s];
	col.end();
}

bool NoArvore::defineVariavelBranching(int& inicioArco, int& fimArco, int &satellite, int &rhs){
	int valorInt, tamRota;
	double fracao, violacao = 0;
	vector<int> verticesRota;
	IloNumArray lambdaValues ( mCplex->env ), gammaValues( mCplex->env );

	satellite = 0;
	mCplex->cplex.getValues(lambdaValues, mCplex->lambda);
	for ( int r = 0; r < mCplex->quantRotasL1; ++r )
	{
		if ( lambdaValues[r] > 0.00001 )
		{
			valorInt = lambdaValues[r];
			fracao = lambdaValues[r] - valorInt;
			if ( min( fracao, ( 1-fracao ) ) > ( violacao + 0.00001 ) )
			{
				violacao = min( fracao, ( 1-fracao ) );
				inicioArco = r; rhs = (valorInt+1);
			}
		}
	}
	lambdaValues.end();
	if ( violacao > 0.00001 )
	{
		gammaValues.end();
		return true;
	}

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
		mCplex->cplex.getValues(gammaValues, mCplex->gamma[s]);
		for (int r = 0; r < mCplex->qRotasL2[s]; ++r)
		{
			if ( ( gammaValues[r] > 0.00001 ) && ( gammaValues[r] < 0.99999 ) ) //gamma[s][r] da raiz eh fracionaria
			{
				verticesRota = mCplex->ptrRotasL2[s][r]->getVertices();
				tamRota = verticesRota.size()-1;
				for (int i = 0; i < tamRota; ++i)
				{
					matrizViolacao[verticesRota[i]][verticesRota[i+1]] += gammaValues[r];
				}
			}
		}

		//depois computa quais variaveis fracionarias do no
		if (qRotas_no[s] > 0)
		{
			mCplex->cplex.getValues(gammaValues, mCplex->gamma_no[s]);
			for (int r = 0; r < qRotas_no[s]; ++r)
			{
				if ( ( gammaValues[r] > 0.00001 ) && ( gammaValues[r] < 0.99999 ) )
				{
					verticesRota = ptrRotas_no[s][r]->getVertices();
					tamRota = verticesRota.size()-1;
					for (int i = 0; i < tamRota; ++i)
					{
						matrizViolacao[verticesRota[i]][verticesRota[i+1]] += gammaValues[r];
					}
				}
			}
		}

		//obtem o VALOR do(s) arco(s) mais violado(s)
		for ( int i = 0; i < nVertices; ++i )
		{
			for ( int j = 0; j < nVertices; ++j )
			{
				if ( min( matrizViolacao[i][j], ( 1-matrizViolacao[i][j] ) ) > ( violacao + 0.00001 ) )
				{
					violacao = min( matrizViolacao[i][j], ( 1-matrizViolacao[i][j] ) );
					inicioArco = i; fimArco = j; satellite = s;
				}
			}
		}
	}

	//libera memoria
	for (int i = 0; i < nVertices; i++) delete [] matrizViolacao[i];
	delete [] matrizViolacao;
	gammaValues.end();

	if ( violacao > 0.00001 ) return true;//existem arcos para executar branching
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

