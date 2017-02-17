#include "NoArvore.h"
#include "geraRotas.h"
#include "ModeloCplex.h"
#include "ListaNosAtivos.h"

int main(int argc, char** argv){

	if ( ( argc != 3 ) && ( argc != 4 ) )
	{
		printf("Parametros esperados:\n (1)Arquivo de instancias\n (2)Solução do subproblema {G, B, N}\n (3)Semente de números aleatórios(opcional)\n\n");
		exit(0);
	}
	if ( ( argv[2][0] != 'G' ) && (argv[2][0] != 'N') )
	{
		printf("A opção para solução do subproblema deve ser { G, N }\n\n");
		exit(0);
	}

	int tempoInicio = time(0);
	int seed = ( argc > 3 ) ? atoi(argv[3]) : tempoInicio;
	srand( seed );
	printf("instancia = %s\n", argv[1]);
	printf("opcaoSub = %c\n", argv[2][0]);
	printf("seed = %d\n", seed);

	Grafo* G = new Grafo(argv[1]);
	vector < Rota* > rLambda = geraRotasLambda( G );
	vector < Rota* > rGamma = geraRotasGamma( G );
	ModeloCplex *mCplex = new ModeloCplex(G, rLambda, rGamma);
	
	int sat = G->getNumSatellites();
	bool alcancouRaiz = true;
	int aux, tmp;
	do
	{
		alcancouRaiz = true;
		for ( int s = 1; s <= sat; ++s )
		{
			mCplex->solveMaster();
			double valSub = mCplex->getValorSubtrair( G, s );
			mCplex->updateDualCosts(G, s);

			if ( argv[2][0] == 'G' )
			{
				GRASP grasp(G, 0.7);
				aux = grasp.run( G->getNumCustomers()*20, valSub );
				if ( aux > 0 )
				{
					alcancouRaiz = false;
					for (int h = 0; h < aux; ++h) mCplex->insertColumn( G, grasp.getRotaConvertida( h, s ), s );
				}
				else
				{
					int tempoBC = time(0);
					ModeloBC bc(G, s, valSub);
					bc.calculaCaminhoElementar(G);
					aux = bc.rotasNegativas.size();
					if ( aux > 0 )
					{
						alcancouRaiz = false;
						for ( int i = 0; i < aux; ++i )
						{
							mCplex->insertColumn( G, bc.rotasNegativas[i], s );
						}
					}
				}
			}
			else if ( argv[2][0] == 'N' )
			{
				NaoElementar ne( G );
				ne.calculaMatrizCiclo2( G );
				vector <Rota*> rrr = ne.getRotaCustoMinimoCiclo2(G, s, valSub);
				if ( rrr.size() > 0 )
				{
					alcancouRaiz = false;
					mCplex->insertColumn(G, rrr[0], s);
				}
			}
		}
		fflush(stdout);
	}
	while(!alcancouRaiz);

	//A partir deste ponto, a raiz do master foi alcançada e deve ser realizado o branch-and-price
	//A raiz do master representara o primeiro no da arvore, que sera ramificada com a insercao das
	//restricoes de branching. Sera feito um branching em arcos, estendendo a formulacao por rotas
	//Os custos de troca serao uma consequencia das rotas obtidas pelo branching em arcos
	ptrNo noAtual;
	vector < ptrNo > novosNos;
	ListaNosAtivos arvBranching;
	double lP, lD = mCplex->solveMaster();
	NoArvore *raiz = new NoArvore(mCplex, G, novosNos, lD);
	printf("NoCount\t\tNosRestantes\tLowerBound\t\tUpperBound\t\tGAP\t\tTempo\n");
	fflush(stdout);

	//delete raiz; 
	NoArvore::addNoProcessado(raiz);
	for (int i = 0; i < novosNos.size(); ++i) arvBranching.insereNo(novosNos[i]);

	int it = 0;
	while( !arvBranching.vazia() )
	{
		novosNos.clear();
		noAtual = arvBranching.retornaProximo();

		lD = noAtual->getLimiteDual();
		lP = ModeloCplex::getLimitePrimal();
		printf(" %d\t\t  %d\t\t %0.3f\t\t %0.3f\t\t%0.2f\t\t%ld\n", ++it, NoArvore::getTotalNosAtivos(), lD, lP, ((lP - lD) / lP) * 100, time(0)-tempoInicio);
		fflush(stdout);

		lP = noAtual->executaBranching( novosNos, argv[2][0] );

		if ( lP < ModeloCplex::getLimitePrimal() )
		{
			arvBranching.podaPorLimiteDual(lP);
			ModeloCplex::setLimitePrimal(lP);
		}

		//delete noAtual;
		NoArvore::addNoProcessado(noAtual);

		for (int i = 0; i < novosNos.size(); ++i) arvBranching.insereNo(novosNos[i]);
		if ( (time(0) - tempoInicio) > 2000 ) break;
	}

	noAtual = NoArvore::getNoProcessado();
	while ( noAtual != NULL )
	{
		printf("limitePrimal = %f\n", noAtual->runCplexInt());
		delete noAtual;
		noAtual = NoArvore::getNoProcessado();
	}

	while( !arvBranching.vazia() ) delete arvBranching.retornaProximo();
	for ( int i = 0; i < rLambda.size(); ++i ) delete rLambda[i];
	delete mCplex;
	delete G;
	return 0;
}
