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
	int seed = ( argc > 4 ) ? atoi(argv[4]) : tempoInicio;
	srand( seed );

	char opSub = ( argc > 2 ) ? argv[2][0] : 'N';
	if ( ( opSub != 'G' ) && (opSub != 'N') )
	{
		printf("A opção para solução do subproblema deve ser { G, N }\n\n");
		exit(0);
	}

	char sepCuts[6];
	if ( argc > 3 )
	{
		if ( strlen(argv[3]) != 5 )
		{
			printf("Use 5 caracteres com valores {0,1,2} para definir o uso de cada um dos cortes!\n");
			exit(0);
		}
		for ( int i = 0; i < 5; ++i )
		{
			if ( ( argv[3][i] != '0' ) && ( argv[3][i] != '1' ) && ( argv[3][i] != '2' ) )
			{
				printf("Use 5 caracteres com valores {0,1,2} para definir o uso de cada um dos cortes!\n");
				exit(0);
			}
		}
		for ( int i = 0; i < 5; ++i ) sepCuts[i] = argv[3][i];
	}
	else
	{
                for ( int i = 0; i < 5; ++i )
                {
                        int r = rand() % 3;
                        if ( r == 0 ) sepCuts[i] = '0';
                        if ( r == 1 ) sepCuts[i] = '1';
                        if ( r == 2 ) sepCuts[i] = '2';
                }

	}
	sepCuts[5] = '\0';

	printf("instancia = %s\n", argv[1]);
	printf("opcaoSub = %c\n", opSub);
	printf("sepCuts = %s\n", sepCuts);
	printf("seed = %d\n", seed);

	Grafo* G = new Grafo(argv[1]);
	vector < Rota* > rLambda = geraRotasLambda( G );
	vector < Rota* > rGamma = geraRotasGamma( G );
	ModeloCplex *mCplex = new ModeloCplex(G, rLambda, rGamma);
	for ( int i = 0; i < 5; ++i ) mCplex->sepCuts[i] = sepCuts[i];

	int aux, tmp = 0;
	int sat = G->getNumSatellites();
	int* pricingSat = new int[sat+1];
	int* qRotas_no_temp = new int[sat+1];
	for ( int s = 0; s <= sat; ++s ) qRotas_no_temp[s] = 0;
	bool alcancouRaiz, bcpCap, bcpComb, bcpMStar, bcpCycle, bcpRoute;

	do
	{
		memset(pricingSat, 0, (sat+1)*sizeof(int));
		do
		{
			alcancouRaiz = true;
			for ( int s = 1; s <= sat; ++s )
			{
				if ( pricingSat[s] == 0 ) //o 'estado' 0 indica que devem ser precificadas rotas para o satellite
				{
					++tmp;
					mCplex->solveMaster( true );
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
							for ( int i = 1; i <= sat; ++i ) pricingSat[i] = ( pricingSat[i] == 1 ) ? -1 : pricingSat[i]; //atualiza o estado do pricing para os satellites
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
								for ( int i = 0; i < aux; ++i ) mCplex->insertColumn( G, bc.rotasNegativas[i], s );
								for ( int i = 1; i <= sat; ++i ) pricingSat[i] = ( pricingSat[i] == 1 ) ? -1 : pricingSat[i]; //atualiza o estado do pricing para os satellites
							}
							else
							{
								pricingSat[s] = 1; //'estado' 1 significa que o satellite nao precisa mais ser precificado
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

							if ( rrr.size() > 1 ) mCplex->insertColumn(G, rrr[1], s);
							for ( int i = 1; i <= sat; ++i ) pricingSat[i] = ( pricingSat[i] == 1 ) ? -1 : pricingSat[i]; //'estado' -1 significa que o satellite i sera precificado no futuro
						}
						else
						{
							pricingSat[s] = 1; //'estado' 1 significa que o satellite nao precisa mais ser precificado
						}
					}
				}
				fflush(stdout);
			}

			if ( alcancouRaiz )
			{
				aux = 0;
				for ( int i = 1; i <= sat; ++i ) aux += pricingSat[i];
				if ( aux < sat ) memset(pricingSat, 0, (sat+1)*sizeof(int));
				else break;
			}
		}
		while( true );

		mCplex->solveMaster();
		bcpCap = ( mCplex->sepCuts[0] != '0' ) ? mCplex->addCapacityCuts(G, qRotas_no_temp, NULL) : false;
		if ( ( mCplex->sepCuts[0] != '0' ) && ( ( mCplex->sepCuts[1] != '0' ) || ( mCplex->sepCuts[2] != '0' ) || ( mCplex->sepCuts[3] != '0' ) || ( mCplex->sepCuts[4] != '0' ) ) ) mCplex->solveMaster();
		bcpComb = ( mCplex->sepCuts[1] != '0' ) ? mCplex->addCombCuts(G, qRotas_no_temp, NULL) : false;
		if ( ( mCplex->sepCuts[1] != '0' ) && ( ( mCplex->sepCuts[2] != '0' ) || ( mCplex->sepCuts[3] != '0' ) || ( mCplex->sepCuts[4] != '0' ) ) ) mCplex->solveMaster();
		bcpMStar = ( mCplex->sepCuts[2] != '0' ) ? mCplex->addMultiStarCuts(G, qRotas_no_temp, NULL) : false;
		if ( ( mCplex->sepCuts[2] != '0' ) && ( mCplex->sepCuts[3] != '0' ) || ( mCplex->sepCuts[4] != '0' ) ) mCplex->solveMaster();
		bcpCycle = ( mCplex->sepCuts[3] != '0' ) ? mCplex->addCycleCuts(G, qRotas_no_temp, NULL, NULL) : false;
		if ( ( mCplex->sepCuts[3] != '0' ) && ( mCplex->sepCuts[4] != '0' ) ) mCplex->solveMaster();
		bcpRoute = ( mCplex->sepCuts[4] != '0' ) ? mCplex->addRouteCuts(G, qRotas_no_temp, NULL, NULL) : false;
	}
	while ( bcpCap || bcpComb || bcpMStar || bcpCycle || bcpRoute );

	mCplex->processaBase( tmp );

	delete [] pricingSat;
	delete [] qRotas_no_temp;
	double lP, lD = mCplex->solveMaster();
	printf("NoCount\t\tNosRestantes\tLowerBound\t\tUpperBound\t\tGAP\t\tTempo\n");
	printf(" 1\t\t  0\t\t %0.3f\t\t %0.3f\t\t%0.2f\t\t%ld\n", lD, lP, ((lP - lD) / lP) * 100, time(0)-tempoInicio);
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
		if ( it % 10 == 1 )
		{
			printf(" %d\t\t  %d\t\t %0.3f\t\t %0.3f\t\t%0.2f\t\t%ld\n", it, NoArvore::getTotalNosAtivos(), lD, lP, ((lP - lD) / lP) * 100, time(0)-tempoInicio);
			fflush(stdout);
		}

		lP = noAtual->executaBranching( novosNos, opSub );

		if ( lP < ModeloCplex::getLimitePrimal() )
		{
			arvBranching.podaPorLimiteDual(lP);
			ModeloCplex::setLimitePrimal(lP);
		}

		//delete noAtual;
		NoArvore::addNoProcessado(noAtual);

		for (int i = 0; i < novosNos.size(); ++i) arvBranching.insereNo(novosNos[i]);
		if ( (time(0) - tempoInicio) > 10000 ) break;
		++it;
	}

	printf("totalCuts = %d\n\n", mCplex->getNumCuts());

	noAtual = NoArvore::getNoProcessado();
	if ( ( noAtual ) && ( (time(0) - tempoInicio) > 10000 ) ) //time limit sem otimo obtido
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
		printf("[%0.2f, %0.2f] (%s)\n", lD, ModeloCplex::getLimitePrimal(), sepCuts);
	}

	while( !arvBranching.vazia() ) delete arvBranching.retornaProximo();
	for ( int i = 0; i < rLambda.size(); ++i ) delete rLambda[i];
	delete mCplex;
	delete G;
	return 0;
}
