#include "NaoElementar.h"

NaoElementar::~NaoElementar(){
	for (int i = 0; i < nLinhas; ++i)
	{
		delete [] matriz[i];
		delete [] matrizArestas[i];
	}
	delete [] linhas;
	delete [] matriz;
	delete [] matrizArestas;
}

NaoElementar::NaoElementar( Grafo* g ){
	nSatellites = g->getNumSatellites();
	mdc = calculaMDC( g );
	nLinhas = g->getNumCustomers()+1;
	nColunas = (g->getCapacVeiculosL2() / mdc)+1;

	matriz = new Celula* [nLinhas]; //a linha 0 representa o deposito
	matrizArestas = new double* [nLinhas];
	for ( int i = 0; i < nLinhas; ++i )
	{
		matriz[i] = new Celula [nColunas]; //a coluna 0 eh usada na inicializacao da PD
		for ( int j = 0; j < nColunas; ++j )
		{
			matriz[i][j].custoMelhorAntec = 0;
			matriz[i][j].melhorAntec = 0;
			matriz[i][j].custoSegMelhorAntec = 0;
			matriz[i][j].segMelhorAntec = 0;
		}
		matrizArestas[i] = new double [nLinhas];
		for ( int j = 0; j < nLinhas; ++j ) matrizArestas[i][j] = g->getCustoArestaDual(i, j);
	}

	//atribui os valores para o cabecalho das linhas
	//os valores do cabecalho sao: indice e peso racionalizado (apos usar MDC)
	linhas = new CabecalhoLinha [nLinhas+1];
	linhas[0].indice = 0;
	linhas[0].pesoRacionalizado = 0;
	for ( int i = 1; i <= g->getNumCustomers(); ++i )
	{
		linhas[i].indice = i;
		linhas[i].pesoRacionalizado = g->getCargaVertice(nSatellites + i) / mdc;
	}
}


void NaoElementar::imprimirMatriz(bool melhorAntecessor){
	if (melhorAntecessor){	
		for (int i = 0; i < nLinhas; ++i){
			for (int j = 0; j < nColunas; ++j){
				if (matriz[i][j].custoMelhorAntec != MAIS_INFINITO){
					printf ("%05.01f(%02d) ", matriz[i][j].custoMelhorAntec, matriz[i][j].melhorAntec);
				}else{
					printf ("-----(-1) ");
				}
			}
			printf ("\n");
		}
	}else{
		for (int i = 0; i < nLinhas; ++i){
			for (int j = 0; j < nColunas; ++j){
				if (matriz[i][j].custoSegMelhorAntec != MAIS_INFINITO){
					printf ("%05.01f(%02d) ", matriz[i][j].custoSegMelhorAntec, matriz[i][j].segMelhorAntec);
				}else{
					printf ("-----(-1) ");
				}
			}
			printf ("\n");
		}
	}
}


void NaoElementar::imprimirCabecalhoLinhas(){
	printf("Linhas:\n");
	for (int i = 0; i < nLinhas; ++i){
		printf ("  (%d) [%d | %d]\n", i, linhas[i].indice, linhas[i].pesoRacionalizado);
	}
}


int NaoElementar::calculaMDC(Grafo* g){
	int inicio = 1, fim = g->getNumCustomers();
	int tamP =  fim+1;
	int* p = new int [tamP];
	p[0] = g->getCapacVeiculosL2();
	int j, menor = MAIS_INFINITO;

	for ( int i = 1; i <= fim; ++i )
	{
		p[i] = g->getCargaVertice(nSatellites + i);
		if ( p[i] < menor )
		{
			menor = p[i]+1;
		}
	}

	//calcula o MDC
	do
	{
		--menor;
		for ( j = 0; j < tamP; ++j )
		{
			if ( p[j] % menor != 0 )
			{
				break;
			}
		}
	}
	while( j < tamP );

	delete [] p;
	return menor;
}


void NaoElementar::calculaMatriz(Grafo* g){
	//No passo base, atribui a todas as celulas da primeira coluna MAIS_INFINITO
	//exceto para o deposito (0), cujo custo para alcancar sera 0
	matriz[0][0].custoMelhorAntec = 0;
	matriz[0][0].melhorAntec = -1;
	for (int i = 1; i < nLinhas; ++i){
		matriz[i][0].custoMelhorAntec = MAIS_INFINITO;
		matriz[i][0].melhorAntec =-1;
	}

	//Nas iteracoes de PD, sera utilizado o peso racionalizado (com o MDC)
	//que eh armazenado no vetor linhas (do tipo CabecalhoLinha)
	//Os estagios de PD serao definidos no limite de carga (a carga eh o recurso neste caso)
	//A cada estagio, aumenta-se uma unidade no limite de carga e verifica-se para cada vertice
	//se eh melhor ir do deposito direto para o vertice ou se eh melhor alcanca-lo por algum outro vertice
	//Ao chegar no ultimo estagio, deve-se verificar qual o vertice eh mais barato para se alcançar 
	//o deposito artificial, mas isto sera feito na reconstrucao do caminho, em outra funcao
	float menorCusto;
	int colAnterior, antecRetorno;
	for ( int q = 1; q < nColunas; ++q )
	{
		matriz[0][q].custoMelhorAntec = 0;
		matriz[0][q].melhorAntec = -1;
		for ( int i = 1; i < nLinhas; ++i )
		{
			colAnterior = q - linhas[i].pesoRacionalizado;
			if ( colAnterior >= 0 )
			{
				menorCusto = getMelhorAntecessor(g, colAnterior, linhas[i].indice, antecRetorno);
				if (menorCusto < ( matriz[i][q-1].custoMelhorAntec - 0.0001 ) )
				{
					matriz[i][q].custoMelhorAntec = menorCusto;
					matriz[i][q].melhorAntec = antecRetorno;
				}
				else
				{
					matriz[i][q].custoMelhorAntec = matriz[i][q-1].custoMelhorAntec;
					matriz[i][q].melhorAntec = matriz[i][q-1].melhorAntec;
				}
			}
			else
			{
				matriz[i][q].custoMelhorAntec = matriz[i][q-1].custoMelhorAntec;
				matriz[i][q].melhorAntec = matriz[i][q-1].melhorAntec;
			}	
		}
	}
	
	//Neste ponto gera-se uma excecao a regra de PD, pois a celula matriz[0][nColunas-1]
	//sera calculada em funcao da ultima coluna (nColunas-1) e nao em funcao das colunas anteriores
	//Esta excecao eh para armazenar na celula matriz[0][nColunas-1] o valor do caminho ate
	//o deposito artificial com a capacidade maxima de recursos, finalizando o calculo da PD
	matriz[0][nColunas-1].custoMelhorAntec = MAIS_INFINITO;
	for ( int i = 1; i < nLinhas; ++i )
	{
		menorCusto = matriz[i][nColunas-1].custoMelhorAntec + matrizArestas[linhas[i].indice][0];
		if (menorCusto < ( matriz[0][nColunas-1].custoMelhorAntec - 0.0001 ) )
		{ 
			matriz[0][nColunas-1].custoMelhorAntec = menorCusto;
			matriz[0][nColunas-1].melhorAntec = linhas[i].indice;
		}
	}
}


float NaoElementar::getMelhorAntecessor(Grafo* g, int coluna, int sucessor, int& antecessorRetorno){
	float menorCusto = MAIS_INFINITO;
	float custoChegada;
	for ( int i = 0; i < nLinhas; ++i )
	{
		custoChegada = matriz[i][coluna].custoMelhorAntec + matrizArestas[linhas[i].indice][sucessor];
		if ( custoChegada < ( menorCusto - 0.0001 ) )
		{
			menorCusto = custoChegada;
			antecessorRetorno = linhas[i].indice;
		}
	}
	return menorCusto;
}


vector<Rota*> NaoElementar::getRotaCustoMinimo(Grafo* g, int s, float valorSub){
	vector < Rota* > rotaNegativa;
	if ( matriz[0][nColunas-1].custoMelhorAntec < ( valorSub - 0.005 ) )
	{
		Rota* r = new Rota();
		r->inserirVerticeFim(s);
		int linhaMat = matriz[0][nColunas-1].melhorAntec;
		int colunaMat = nColunas-1;
		int tmp;

		while( matriz[linhaMat][colunaMat].melhorAntec != -1 )
		{
			if ( matriz[linhaMat][colunaMat].custoMelhorAntec < ( matriz[linhaMat][colunaMat-1].custoMelhorAntec - 0.0001 ) )
			{
				r->inserirVerticeInicio(nSatellites + linhas[linhaMat].indice);
				tmp = colunaMat;
				colunaMat = colunaMat - linhas[linhaMat].pesoRacionalizado;
				linhaMat = ( matriz[linhaMat][tmp].melhorAntec != 0 ) ? ( matriz[linhaMat][tmp].melhorAntec ) : 0;
			}
			else
			{
				--colunaMat;
			}
		}
		r->inserirVerticeInicio(s);
		r->setCustoReduzido( matriz[0][nColunas-1].custoMelhorAntec );
		r->setCusto(g);
		rotaNegativa.push_back(r);
	}
	return rotaNegativa;
}


void NaoElementar::calculaMatrizCiclo2( Grafo* g ){
	//No passo base, atribui a todas as celulas da primeira coluna MAIS_INFINITO
	//exceto para o deposito (0), cujo custo para alcancar sera 0
	matriz[0][0].custoMelhorAntec = 0;
	matriz[0][0].melhorAntec = -1;
	matriz[0][0].custoSegMelhorAntec = MAIS_INFINITO;
	matriz[0][0].segMelhorAntec =-1;		
	for (int i = 1; i < nLinhas; ++i)
	{
		matriz[i][0].custoMelhorAntec = MAIS_INFINITO;
		matriz[i][0].melhorAntec =-1;
		matriz[i][0].custoSegMelhorAntec = MAIS_INFINITO;
		matriz[i][0].segMelhorAntec =-1;
	}

	//A cada iteracao procura-se pelo melhor e pelo segundo melhor antecessor de i
	//O fato de se armazenar tambem o segundo melhor antecessor eh pelo fato de que se o melhor
	//antecessor tiver como antecessor i, este "melhor" caminho nao eh valido, e portanto, deve-se
	//usar o segundo melhor caminho (armazenado pelo segundo melhor antecessor)
	//Obs: O segundo melhor antecessor so sera verificado caso o antecessor do primeiro melhor seja i 
	double custoMelhorAntec, custoSegMelhorAntec;
	int melhorAntec, segMelhorAntec;
	int colAnterior;
	for (int q = 1; q < nColunas; ++q)
	{
		matriz[0][q].custoMelhorAntec = matriz[0][q].custoSegMelhorAntec = 0;
		matriz[0][q].melhorAntec = matriz[0][q].segMelhorAntec = -1;
		for (int i = 1; i < nLinhas; ++i)
		{
			colAnterior = q - linhas[i].pesoRacionalizado;
			if (colAnterior >= 0)
			{
				custoMelhorAntec = getMelhorSegMelhorAntecessor(g, colAnterior, linhas[i].indice, melhorAntec, custoSegMelhorAntec, segMelhorAntec);
				if ((custoMelhorAntec < ( matriz[i][q-1].custoMelhorAntec - 0.00001 ) ) && (melhorAntec != matriz[i][q-1].melhorAntec))
				{
					matriz[i][q].custoMelhorAntec = custoMelhorAntec;
					matriz[i][q].melhorAntec = melhorAntec;
					
					if (custoSegMelhorAntec < ( matriz[i][q-1].custoMelhorAntec - 0.00001 ) )
					{
						matriz[i][q].custoSegMelhorAntec = custoSegMelhorAntec;
						matriz[i][q].segMelhorAntec = segMelhorAntec;	
					}
					else
					{
						matriz[i][q].custoSegMelhorAntec = matriz[i][q-1].custoMelhorAntec;
						matriz[i][q].segMelhorAntec = matriz[i][q-1].melhorAntec;
					}
				}
				else if ( ( custoMelhorAntec < ( matriz[i][q-1].custoMelhorAntec - 0.00001 ) ) && ( melhorAntec == matriz[i][q-1].melhorAntec ) )
				{
					matriz[i][q].custoMelhorAntec = custoMelhorAntec;
					matriz[i][q].melhorAntec = matriz[i][q-1].melhorAntec;

					if (custoSegMelhorAntec < ( matriz[i][q-1].custoSegMelhorAntec - 0.00001 ) )
					{
						matriz[i][q].custoSegMelhorAntec = custoSegMelhorAntec;
						matriz[i][q].segMelhorAntec = segMelhorAntec;
					}
					else
					{
						matriz[i][q].custoSegMelhorAntec = matriz[i][q-1].custoSegMelhorAntec;
						matriz[i][q].segMelhorAntec = matriz[i][q-1].segMelhorAntec;
					}
				}
				else if ( ( custoMelhorAntec < ( matriz[i][q-1].custoSegMelhorAntec - 0.00001 ) ) && ( melhorAntec != matriz[i][q-1].melhorAntec ) )
				{
					matriz[i][q].custoMelhorAntec = matriz[i][q-1].custoMelhorAntec;
					matriz[i][q].melhorAntec = matriz[i][q-1].melhorAntec;
					matriz[i][q].custoSegMelhorAntec = custoMelhorAntec;
					matriz[i][q].segMelhorAntec = melhorAntec;	
				}
				else if ( ( custoSegMelhorAntec < ( matriz[i][q-1].custoSegMelhorAntec - 0.00001 ) ) && ( segMelhorAntec != matriz[i][q-1].melhorAntec ) )
				{
					matriz[i][q].custoMelhorAntec = matriz[i][q-1].custoMelhorAntec;
					matriz[i][q].melhorAntec = matriz[i][q-1].melhorAntec;
					matriz[i][q].custoSegMelhorAntec = custoSegMelhorAntec;
					matriz[i][q].segMelhorAntec = segMelhorAntec;
				}
				else
				{
					matriz[i][q].custoMelhorAntec = matriz[i][q-1].custoMelhorAntec;
					matriz[i][q].melhorAntec = matriz[i][q-1].melhorAntec;
					matriz[i][q].custoSegMelhorAntec = matriz[i][q-1].custoSegMelhorAntec;
					matriz[i][q].segMelhorAntec = matriz[i][q-1].segMelhorAntec;
				}
			}
			else
			{
				matriz[i][q].custoMelhorAntec = matriz[i][q-1].custoMelhorAntec;
				matriz[i][q].melhorAntec = matriz[i][q-1].melhorAntec;
				matriz[i][q].custoSegMelhorAntec = matriz[i][q-1].custoSegMelhorAntec;
				matriz[i][q].segMelhorAntec = matriz[i][q-1].segMelhorAntec;
			}	
		}
	}
	
	//Neste ponto gera-se uma excecao a regra de PD, pois a celula matriz[0][nColunas-1]
	//sera calculada em funcao da ultima coluna (nColunas-1) e nao em funcao das colunas anteriores
	//Esta excecao eh para armazenar na celula matriz[0][nColunas-1] o valor do caminho ate
	//o deposito artificial com a capacidade maxima de recursos, finalizando o calculo da PD
	matriz[0][nColunas-1].custoMelhorAntec = MAIS_INFINITO;
	for ( int i = 1; i < nLinhas; ++i )
	{
		custoMelhorAntec = matriz[i][nColunas-1].custoMelhorAntec + matrizArestas[linhas[i].indice][0];
		if ( custoMelhorAntec < ( matriz[0][nColunas-1].custoMelhorAntec - 0.00001 ) )
		{ 
			matriz[0][nColunas-1].melhorAntec = linhas[i].indice;
			matriz[0][nColunas-1].custoMelhorAntec = custoMelhorAntec;
		}
	}

	matriz[0][nColunas-1].custoSegMelhorAntec = MAIS_INFINITO;
	for ( int i = 1; i < nLinhas; ++i )
	{
		if ( linhas[i].indice != matriz[0][nColunas-1].melhorAntec )
		{
			custoSegMelhorAntec = matriz[i][nColunas-1].custoSegMelhorAntec + matrizArestas[linhas[i].indice][0];
			if ( custoSegMelhorAntec < ( matriz[0][nColunas-1].custoSegMelhorAntec - 0.00001 ) )
			{ 
				matriz[0][nColunas-1].segMelhorAntec = linhas[i].indice;
				matriz[0][nColunas-1].custoSegMelhorAntec = custoSegMelhorAntec;
			}
		}
	}
}


double NaoElementar::getMelhorSegMelhorAntecessor(Grafo* g, int col, int suc, int& melhorAnt, double& custoSegAnt, int& segAnt){
	double custoMelhorAnt = custoSegAnt = MAIS_INFINITO;
	double custoChegada;
	for (int i = 0; i < nLinhas; ++i)
	{
		if ((matriz[i][col].melhorAntec != suc) && (matriz[i][col].custoMelhorAntec < MAIS_INFINITO))
		{
			custoChegada = matriz[i][col].custoMelhorAntec + matrizArestas[linhas[i].indice][suc];
			if ( custoChegada < ( custoMelhorAnt - 0.00001 ) )
			{
				segAnt = melhorAnt;
				custoSegAnt = custoMelhorAnt;
				melhorAnt = linhas[i].indice;
				custoMelhorAnt = custoChegada;
			}
			else if (custoChegada < ( custoSegAnt - 0.00001 ) )
			{
				segAnt = linhas[i].indice;
				custoSegAnt = custoChegada;
			}
		}
		else if (matriz[i][col].custoSegMelhorAntec < MAIS_INFINITO)
		{
			custoChegada = matriz[i][col].custoSegMelhorAntec + matrizArestas[linhas[i].indice][suc];
			if ( custoChegada < ( custoMelhorAnt - 0.00001 ) )
			{
				segAnt = melhorAnt;
				custoSegAnt = custoMelhorAnt;
				melhorAnt = linhas[i].indice;
				custoMelhorAnt = custoChegada;
			}
			else if ( custoChegada < ( custoSegAnt - 0.00001 ) )
			{
				segAnt = linhas[i].indice;
				custoSegAnt = custoChegada;
			}
		}
	}
	
	return custoMelhorAnt;
}


vector < Rota* > NaoElementar::getRotaCustoMinimoCiclo2(Grafo* g, int s, double valorSub){
	vector < Rota* > rotaNegativa;

	if ( matriz[0][nColunas-1].custoMelhorAntec < ( valorSub - 0.0001 ) )
	{
		Rota* r = new Rota();
		r->inserirVerticeFim(s);
		int linhaMat = matriz[0][nColunas-1].melhorAntec;
		int colunaMat = nColunas-1;
		int tmp, melhorAntecPassado = -1000;

		while( matriz[linhaMat][colunaMat].melhorAntec != -1 )
		{
			if ( matriz[linhaMat][colunaMat].melhorAntec != melhorAntecPassado )
			{
				if ( matriz[linhaMat][colunaMat].custoMelhorAntec < ( matriz[linhaMat][colunaMat-1].custoMelhorAntec - 0.00001 ) )
				{
					melhorAntecPassado = linhas[linhaMat].indice;
					r->inserirVerticeInicio(nSatellites + melhorAntecPassado);
					tmp = colunaMat;
					colunaMat = colunaMat - linhas[linhaMat].pesoRacionalizado;
					linhaMat = matriz[linhaMat][tmp].melhorAntec;
				}
				else
				{
					--colunaMat;
				}
			}
			else
			{
				if ( matriz[linhaMat][colunaMat].custoSegMelhorAntec < ( matriz[linhaMat][colunaMat-1].custoSegMelhorAntec - 0.00001 ) )
				{
					melhorAntecPassado = linhas[linhaMat].indice;
					r->inserirVerticeInicio(nSatellites + melhorAntecPassado);
					tmp = colunaMat;
					colunaMat = colunaMat - linhas[linhaMat].pesoRacionalizado;
					linhaMat = matriz[linhaMat][tmp].segMelhorAntec;
				}
				else
				{
					--colunaMat;
				}
			}
		}
		r->inserirVerticeInicio(s);
		r->setCustoReduzido(matriz[0][nColunas-1].custoMelhorAntec);
		r->setCusto( g );
		rotaNegativa.push_back(r);

		if ( matriz[0][nColunas-1].custoSegMelhorAntec < ( valorSub - 0.0001 ) )
		{
			r = new Rota();
			r->inserirVerticeFim(s);
			linhaMat = matriz[0][nColunas-1].segMelhorAntec;
			colunaMat = nColunas-1;
			melhorAntecPassado = -1000;

			while( matriz[linhaMat][colunaMat].melhorAntec != -1 )
			{
				if ( matriz[linhaMat][colunaMat].melhorAntec != melhorAntecPassado )
				{
					if ( matriz[linhaMat][colunaMat].custoMelhorAntec < ( matriz[linhaMat][colunaMat-1].custoMelhorAntec - 0.00001 ) )
					{
						melhorAntecPassado = linhas[linhaMat].indice;
						r->inserirVerticeInicio(nSatellites + melhorAntecPassado);
						tmp = colunaMat;
						colunaMat = colunaMat - linhas[linhaMat].pesoRacionalizado;
						linhaMat = matriz[linhaMat][tmp].melhorAntec;
					}
					else
					{
						--colunaMat;
					}
				}
				else
				{
					if ( matriz[linhaMat][colunaMat].custoSegMelhorAntec < ( matriz[linhaMat][colunaMat-1].custoSegMelhorAntec - 0.00001 ) )
					{
						melhorAntecPassado = linhas[linhaMat].indice;
						r->inserirVerticeInicio(nSatellites + melhorAntecPassado);
						tmp = colunaMat;
						colunaMat = colunaMat - linhas[linhaMat].pesoRacionalizado;
						linhaMat = matriz[linhaMat][tmp].segMelhorAntec;
					}
					else
					{
						--colunaMat;
					}
				}
			}
			r->inserirVerticeInicio(s);
			r->setCustoReduzido(matriz[0][nColunas-1].custoSegMelhorAntec);
			r->setCusto( g );
			rotaNegativa.push_back(r);
		}
	}
	return rotaNegativa;
}
