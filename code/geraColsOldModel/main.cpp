#include <stdio.h>
#include "NoArvore.h"
#include "geraRotas.h"
#include "ModeloCplex.h"
#include "ListaNosAtivos.h"

int main(int argc, char** argv){

	if ( ( argc < 2 ) || ( argc > 5 ) )
	{
		printf("Parametros esperados:\n (1) Arquivo de instancias\n (2) Solução do subproblema {G (GRASP + BC), *N (Nao Elementar)}\n (3) Separar cortes? {*0 (NAO), 1 (RAIZ), 2 (TODOS)}\n (4) Semente de números aleatórios\n\n");
		exit(0);
	}

	int tempoInicio = time(0);
	char opSub = ( argc > 2 ) ? argv[2][0] : 'N';
	int sepCuts = ( argc > 3 ) ? atoi(argv[3]) : 0;
	int seed = ( argc > 4 ) ? atoi(argv[4]) : tempoInicio;

	if ( ( opSub != 'G' ) && (opSub != 'N') )
	{
		printf("A opção para solução do subproblema deve ser { G, N }\n\n");
		exit(0);
	}
	if ( ( sepCuts < 0 ) || ( sepCuts > 2 ) )
	{
		printf("A opção para separar cortes deve ser { 0, 1, 2 }\n\n");
		exit(0);
	}

	srand( seed );
	printf("instancia = %s\n", argv[1]);
	printf("opcaoSub = %c\n", opSub);
	printf("sepCuts = %d\n", sepCuts);
	printf("seed = %d\n", seed);

	Grafo* G = new Grafo(argv[1]);
	vector < Rota* > rLambda = geraRotasLambda( G );
	vector < Rota* > rGamma = geraRotasGamma( G );
	ModeloCplex *mCplex = new ModeloCplex(G, rLambda, rGamma, sepCuts);

	bool alcancouRaiz, bcp;
	int aux, tmp, sat = G->getNumSatellites();
	int* qRotas_no_temp = new int[sat+1];
	for ( int s = 0; s <= sat; ++s ) qRotas_no_temp[s] = 0;
	do
	{
		do
		{
			alcancouRaiz = true;
			for ( int s = 1; s <= sat; ++s )
			{
				mCplex->solveMaster();
				double valSub = mCplex->getValorSubtrair( G, s );
				mCplex->updateDualCosts(G, s);

				if ( opSub == 'G' )
				{
					GRASP grasp(G, 0.6);
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
				else if ( opSub == 'N' )
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

		bcp = ( mCplex->sepCuts > 0 ) ? mCplex->addCapacityCuts(G, qRotas_no_temp, NULL) : false;

	}
	while ( bcp );

	delete [] qRotas_no_temp;
	double lP, lD = mCplex->solveMaster();
	printf("NoCount\t\tNosRestantes\tLowerBound\t\tUpperBound\t\tGAP\t\tTempo\n");
	//printf(" 1\t\t  0\t\t %0.3f\t\t %0.3f\t\t%0.2f\t\t%ld\n", lD, lP, ((lP - lD) / lP) * 100, time(0)-tempoInicio);
	fflush(stdout);

	//A partir deste ponto, a raiz do master foi alcançada e deve ser realizado o branch-and-price
	//A raiz do master representara o primeiro no da arvore, que sera ramificada com a insercao das
	//restricoes de branching. Sera feito um branching em arcos, estendendo a formulacao por rotas
	//Os custos de troca serao uma consequencia das rotas obtidas pelo branching em arcos
	ptrNo noAtual;
	vector < ptrNo > novosNos;
	ListaNosAtivos arvBranching;
	NoArvore *raiz = new NoArvore(mCplex, G, novosNos, lD);

	ModeloCplex::setLimitePrimal(raiz->runCplexInt());
	for (int i = 0; i < novosNos.size(); ++i) arvBranching.insereNo(novosNos[i]);
	delete raiz;

	int it = 1;
	while( !arvBranching.vazia() )
	{
		novosNos.clear();
		noAtual = arvBranching.retornaProximo();

		lD = noAtual->getLimiteDual();
		lP = ModeloCplex::getLimitePrimal();
		printf(" %d\t\t  %d\t\t %0.3f\t\t %0.3f\t\t%0.2f\t\t%ld\n", ++it, NoArvore::getTotalNosAtivos(), lD, lP, ((lP - lD) / lP) * 100, time(0)-tempoInicio);
		fflush(stdout);

		lP = noAtual->executaBranching( novosNos, opSub );

		if ( lP < ModeloCplex::getLimitePrimal() )
		{
			arvBranching.podaPorLimiteDual(lP);
			ModeloCplex::setLimitePrimal(lP);
		}

		delete noAtual;
		//NoArvore::addNoProcessado(noAtual);

		for (int i = 0; i < novosNos.size(); ++i) arvBranching.insereNo(novosNos[i]);
		if ( (time(0) - tempoInicio) > 10000 ) break;
	}

	noAtual = NoArvore::getNoProcessado();
	if ( noAtual )
	{
		tempoInicio = time(0);
		double menorNoProcessado = 99999;

		while ( noAtual != NULL )
		{
			lP = noAtual->runCplexInt();
			if ( lP < menorNoProcessado ) menorNoProcessado = lP;
			delete noAtual;
			noAtual = NoArvore::getNoProcessado();
		}
		printf("menorLimitePrimalCplex = %f\ntempo = %d\n", menorNoProcessado, (int)time(0)-tempoInicio);
	}

	while( !arvBranching.vazia() ) delete arvBranching.retornaProximo();
	for ( int i = 0; i < rLambda.size(); ++i ) delete rLambda[i];
	delete mCplex;
	delete G;
	return 0;
}
