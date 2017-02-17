#include "ModeloBC.h"
using namespace std;

ILOINCUMBENTCALLBACK2(solIntCallBack, ModeloBC&, modBC, Grafo*, g){
	//caso a solucao seja de custo negativo, a inclui no conjunto de rotas negativas
	if ( getObjValue() < (modBC.valorSubtrair - 0.001 ) )
	{
		IloArray<IloNumArray> x_values(modBC.env, modBC.nVertices+2);
		for ( int i = 0; i <= modBC.nVertices+1; ++i )
		{
			x_values[i] = IloNumArray(modBC.env);
			getValues(x_values[i], modBC.x[i]);
		}

		//cria a rota e a inclui no vector<Rota*>
		int nSatellites = g->getNumSatellites();
		int proxVertice = 0, depositoArtificial = modBC.nVertices+1;
		Rota* r = new Rota();
		do{
			if ( proxVertice == 0 ) r->inserirVerticeFim( modBC.sat );
			else r->inserirVerticeFim( nSatellites + proxVertice );
			for ( int i = 1; i <= depositoArtificial; ++i )
			{
				if ( x_values[proxVertice][i] > 0.5 )
				{
					proxVertice = i;
					break;
				}
			}
		}while(proxVertice != depositoArtificial);		
		r->inserirVerticeFim( modBC.sat );
		r->setCusto(g);
		r->setCustoReduzido(getObjValue());
		modBC.rotasNegativas.push_back(r);

		//libera memoria
		for (int i = 0; i <= modBC.nVertices+1; ++i) x_values[i].end();
		x_values.end();

//		abort(); //interrompe a execucao do cplex apos encontrar a primeira rota de custo reduzido negativo
	}
}


ILOLAZYCONSTRAINTCALLBACK2(cutCallBack, ModeloBC&, modBC, Grafo*, g){		
	//pega os valores das variaveis de decisao x
	IloArray<IloNumArray> x_values(modBC.env, modBC.nVertices+2);
	for ( int i = 0; i <= modBC.nVertices+1; ++i )
	{
		x_values[i] = IloNumArray(modBC.env);
		getValues(x_values[i], modBC.x[i]);
	}

	//pega os valores das variaveis de decisao y
	IloNumArray y_values(modBC.env);
	getValues(y_values, modBC.y);

	//constroi o grafo de dinic baseando-se nos valores das variaveis de decisao
	Fluxo* flow = new Fluxo(x_values, y_values, modBC.nVertices);

	int random, destino, nVDinic = flow->getNVerticesDinic();
	int* ordem = new int[--nVDinic]; //o vertice 0 nao entra como destino, portanto, diminue 1 no numero de vertices Dinic
	for (int i = 0; i < nVDinic; ++i) ordem[i] = ( i+1 );

	bool inW;
	float maiorValorY;
	int indiceMaiorY, i = nVDinic;

	//executa o algoritmo de fluxo maximo de 0 para cada vertice do grafo Dinic (escolhido aleatoriamente)
	while ( i >= 1 )
	{
		random = rand() % i;
		destino = ordem[random];
		ordem [random] = ordem[--i];

		//executa o algoritmo de fluxo maximo, que retornara o conjunto W (ou um conjunto vazio, caso nao seja possivel obter corte)
		vector < int > conjW = flow->calculaFluxoMaximo(y_values, destino);

		//caso o conjunto nao seja vazio (ou seja, existe um conjunto W para esta iteracao), procura pelo corte x(W, V\W)
		if ( conjW.size() > 0 )
		{
			int tam = conjW.size();
			vector < int > arcosCut;
			IloExpr expCorte(modBC.env);

			for ( int j = 0; j < tam; ++j )
			{
				for (int x = 1; x <= modBC.nVertices; ++x)
				{
					inW = false;
					for (int y = 0; y < tam; ++y)
					{
						if (conjW[y] == x)
						{
							inW = true;
							break;
						}
					}
					if ( !inW )
					{
						expCorte += modBC.x[conjW[j]][x];
						arcosCut.push_back(conjW[j]);
						arcosCut.push_back(x);
					}
				}
			}

			if ( arcosCut.size() > 0 )
			{
				//verifica qual o maiorValorY em V-W para compor a restricao
				maiorValorY = 0;
				for ( int x = 1; x <= modBC.nVertices; ++x )
				{
					if ( y_values[x] > maiorValorY )
					{
						inW = false;
						for ( int y = 0; y < tam; ++y )
						{
							if ( conjW[y] == x )
							{
								inW = true;
								break;
							}
						}
						if ( !inW )
						{
							maiorValorY = y_values[x];
							indiceMaiorY = x;
						}
					}
				}

				//adiciona ao vetor de restricoes violadas e ao vetor de arcos dos cortes
				expCorte -= modBC.y[indiceMaiorY];
				IloRange restrCorte  = ( expCorte >= 0 );
				add( restrCorte );
				expCorte.end();
				restrCorte.end();

				//atualiza o grafo para que as arestas incluidas neste corte nao estejam no proximo corte
				flow->atualizaGrafoDinic(g, y_values, arcosCut);
			}
		}
	}
	delete [] ordem;

	//libera a memoria armazenada pelos objetos do cplex e do fluxo
	for (int i = 0; i <= modBC.nVertices+1; ++i) x_values[i].end();
	x_values.end();
	y_values.end();
	delete flow;
}

ModeloBC::ModeloBC( Grafo* g, int s, float valSub ){
	model = IloModel(env);
	cplex = IloCplex(model);
	nVertices = g->getNumCustomers();
	sat = s;
	valorSubtrair = valSub;

	cplex.setOut(env.getNullStream());
	cplex.setWarning(env.getNullStream());

	initVars();
	setObjectiveFunction( g );
	setConstraints1e2();
	setConstraints3e4();
	setConstraint5( g );
}


ModeloBC::~ModeloBC(){
	env.end();
}


void ModeloBC::initVars(){
	x = IloArray<IloIntVarArray>( env, ( nVertices+2 ) );
	for ( int i = 0; i < (nVertices+2); ++i )
	{
		x[i] = IloIntVarArray(env, (nVertices+2), 0, 1);
		for ( int j = 0; j < (nVertices+2); ++j )
		{
			char nome[20];
			sprintf(nome, "x_%02d_%02d", i, j);
			x[i][j].setName(nome);
		}
	}

	y = IloIntVarArray(env, (nVertices+2), 0, 1);
	for ( int j = 0; j < (nVertices+2); ++j )
	{
		char nome[20];
		sprintf(nome, "y_%02d", j);
		y[j].setName(nome);
	}
}


void ModeloBC::setObjectiveFunction( Grafo *g ){
	IloExpr obj(env);
	int depositoArtificial = nVertices+1;

	for ( int i = 1; i <= nVertices; ++i )
	{
		obj += g->getCustoArestaDual( 0 , i ) * x[0][i];
		obj += g->getCustoArestaDual( i , 0 ) * x[i][depositoArtificial];

		obj += MAIS_INFINITO * x[i][0];
		obj += MAIS_INFINITO * x[depositoArtificial][i];
	}
	
	for ( int i = 1; i <= nVertices; ++i )
	{
		for ( int j = 1; j <= nVertices; ++j )
		{
			obj += g->getCustoArestaDual( i , j ) * x[i][j];
		}
	}

	obj += MAIS_INFINITO * x[0][0];
	obj += MAIS_INFINITO * x[0][depositoArtificial];
	obj += MAIS_INFINITO * x[depositoArtificial][0];
	obj += MAIS_INFINITO * x[depositoArtificial][depositoArtificial];
	obj += MAIS_INFINITO * y[depositoArtificial];
	obj += MAIS_INFINITO * y[0];

	model.add(IloObjective(env, obj));
}


void ModeloBC::setConstraints1e2(){
	IloExpr expr = IloExpr(env);
	IloExpr expr2 = IloExpr(env);
	int depositoArtificial = nVertices+1;

	for (int j = 1; j <= nVertices; ++j)
	{
		expr += x[0][j];
	}
	model.add(expr == 1);

	for ( int i = 1; i <= nVertices; ++i )
	{
		expr2 += x[i][depositoArtificial];
	}
	model.add(expr2 == 1);

	expr.end();	expr2.end();
}

void ModeloBC::setConstraints3e4(){
	int depositoArtificial = nVertices+1;

	for ( int i = 1; i <= nVertices; ++i )
	{
		IloExpr expr = IloExpr(env);
		for ( int h = 0; h <= nVertices; ++h )
		{
			expr += x[h][i];
		}
		expr -= y[i];
		model.add(expr == 0);
		expr.end();


		IloExpr expr2 = IloExpr(env);
		for ( int j = 1; j <= depositoArtificial; ++j )
		{
			expr2 += x[i][j];
		}
		expr2 -= y[i];
		model.add(expr2 == 0);
		expr.end();
	}	
}


void ModeloBC::setConstraint5( Grafo* g ){
	IloExpr expr = IloExpr(env);
	int nSatellites = g->getNumSatellites();
	
	for ( int i = 1; i <= nVertices; ++i ) expr += g->getCargaVertice(nSatellites + i) * y[i];
	model.add( expr <= g->getCapacVeiculosL2() );
	expr.end();
}


void ModeloBC::calculaCaminhoElementar(Grafo* g){
	//inclusao dos callbacks para incluir cortes e capturar solucoes inteiras
	cplex.use( cutCallBack( env, *this, g ) );
	cplex.use( solIntCallBack( env, *this, g ) );

	//parametros para trabalhar com modelos sem incluir todas as restricoes
	cplex.setParam( IloCplex::PreInd, 0 );
	cplex.setParam( IloCplex::AggInd, 0 );
    cplex.setParam( IloCplex::HeurFreq, -1 );
   
	//strong branching
	cplex.setParam( IloCplex::VarSel, 3 );
		
	//seta maior prioridade de branching nas variaveis y
	IloNumArray priorities(env, (nVertices+2));
	priorities[0] = 0; priorities[nVertices+1] = 0;
	for (int i = 1; i <= nVertices; ++i) priorities[i] = 1;
	cplex.setPriorities(y, priorities);

	//executa o cplex com o modelo e os parametros definidos acima
	cplex.solve();
}

