#include "geraRotas.h"
using namespace std;

vector < vector< int > > arranjoRotasLambda( int* vertices, int tam ){
	vector < vector < int > > retorno;
	if (tam == 1)
	{
		vector < int > r;
		r.push_back(vertices[tam]);
		retorno.push_back( r );
	}
	else
	{
		retorno = arranjoRotasLambda( vertices, tam-1 );
		int tamanhoRetorno = retorno.size();
		for (int i = 0; i < tamanhoRetorno; ++i)
		{
			for (int j = 1; j <= retorno[i].size(); ++j)
			{
				retorno.push_back(retorno[i]);
				retorno[retorno.size()-1].insert(retorno[retorno.size()-1].begin()+j, vertices[tam]);
			}
			retorno[i].insert(retorno[i].begin(), vertices[tam]);
		}
	}
	return retorno;
}

vector < Rota* > geraRotasLambda( Grafo* g ){
	int numSat = g->getNumSatellites();
	int *verticesCoding = new int [numSat+1];
	verticesCoding[0] = 0;
	for (int i = 1; i <= numSat; ++i)
	{
		verticesCoding[i] = pow (2, i-1);
	}

	Rota *aux;
	double custoRota;
	int count, indiceMenorRota;
	vector < Rota* > rotasLambda;
	int *vertices = new int [numSat+1];

	for (int i = 1; i < pow (2, numSat); ++i)
	{
		count = 0;
		for (int j = 1; j <= numSat; ++j)
		{
			if ( (i & verticesCoding[j] ) > 0 ) vertices[++count] = j;
		}

		aux = NULL; custoRota = 999999999;
		vector < vector < int > > vertRota = arranjoRotasLambda(vertices, count);
		for(int i = 0; i < vertRota.size(); ++i)
		{
			Rota* r = new Rota();
			r->inserirVerticeFim(0);
			for(int j = 0; j < vertRota[i].size(); ++j)
			{
				r->inserirVerticeFim(vertRota[i][j]);
			}
			r->inserirVerticeFim(0);
			r->setCusto(g);

			if (r->getCusto() < ( custoRota - 0.00001 ) )
			{
				if ( aux != NULL ) delete aux;
				custoRota = r->getCusto();
				aux = r;
			}
			else
			{	
				delete r;
			}
		}
		rotasLambda.push_back( aux );
	}
	delete [] verticesCoding;
	delete [] vertices;
	return rotasLambda;
}


vector < Rota* > geraRotasGamma( Grafo* g, int satellite ){
	int numSate = g->getNumSatellites();
	int numCust = g->getNumCustomers();
	int *verticesCoding = new int [numCust+1];
	verticesCoding[0] = 0;
	for (int i = 1; i <= numCust; ++i)
	{
		verticesCoding[i] = pow (2, i-1);
	}
	
	int count, carga;
	vector < Rota* > rotasGamma;
	int *vertices = new int [numCust+1];

	for (int i = 1; i < pow (2, numCust); ++i)
	{
		count = 0;
		carga = 0;
		for (int j = 1; j <= numCust; ++j)
		{
			if ( (i & verticesCoding[j] ) > 0 )
			{
				vertices[++count] = j;
				carga += g->getCargaVertice(numSate + j);
			}
		}
		
		if ( carga <= g->getCapacVeiculosL2() )
		{
			vector < vector < int > > vertRota = arranjoRotasLambda(vertices, count);
			for(int i = 0; i < vertRota.size(); ++i)
			{
				Rota* r = new Rota();
				r->inserirVerticeFim( satellite );
				for(int j = 0; j < vertRota[i].size(); ++j)
				{
					r->inserirVerticeFim(numSate + vertRota[i][j]);
				}
				r->inserirVerticeFim( satellite );
				r->setCusto(g);
				rotasGamma.push_back( r );
			}
		}
	}

	delete [] verticesCoding;
	delete [] vertices;
	return rotasGamma;
}


vector < Rota* > geraRotasGamma( Grafo* g ){
	bool inseriuVertice;
	int k, v, cargaTotal;
	int capacidade = g->getCapacVeiculosL2();
	int nSatellites = g->getNumSatellites();
	int nCustomers = g->getNumCustomers();
	int maxVeic = g->getMaxVeiculosL2();

	int numVerticesDescobertos;
	ptrNoListaEnc* rotas = new ptrNoListaEnc[maxVeic];
	int* vertices = new int[nCustomers+1];

	for (int it = 0; it < 1000000; ++it){
		//Inicializa os valores do vetor de vertices a serem cobertos para a iteracao
		numVerticesDescobertos = nCustomers;
		for (int i = 1; i <= nCustomers; ++i){
			vertices[i] = i;
		}

		//A principio, cria uma rota com um unico vertice (escolhido aleatoriamente) para cada veiculo
		for (k = 0; k < maxVeic; ++k){
			do
			{
				v = (rand() % numVerticesDescobertos) + 1;
			}
			while( g->getCargaVertice( vertices[v] + nSatellites ) > capacidade );

			rotas[k] = new NoListaEnc();
			rotas[k]->vertice = vertices[v];
			rotas[k]->prox = NULL;

			//Substitue o vertice coberto pelo ultimo do vetor
			vertices[v] = vertices[numVerticesDescobertos];
			--numVerticesDescobertos;
		}

		//Partindo-se das rotas singulares acima, tenta-se inserir aleatoriamente os vertices nas rotas ate cobrir todos
		k = 0;
		ptrNoListaEnc tmp, aux;
		while ( numVerticesDescobertos > 0 )
		{
			inseriuVertice = false;
			v = (rand() % numVerticesDescobertos) + 1; //escolhe-se o indice de um vertice (vertices[v] representa o vertice)

			//Tenta inserir o vertice escolhido em um veiculo, passando para o proximo, caso nao seja possivel inserir no atual
			for ( int veic = k; veic < (k+maxVeic); ++veic )
			{
				//verifica se a cargaTotal excede a capacidade do veiculo
				cargaTotal = g->getCargaVertice( vertices[v] + nSatellites );
				tmp = rotas[veic % maxVeic];
				while( tmp != NULL )
				{
					cargaTotal += g->getCargaVertice( tmp->vertice + nSatellites );
					tmp = tmp->prox;
				}

				if ( cargaTotal <= capacidade )
				{
					//insere sempre o vertice entre 0 e rota[k]. Esta estrategia pode influenciar no custo da solucao, mas nao na sua viabilidade
					tmp = new NoListaEnc();
					tmp->vertice = vertices[v];
					tmp->prox = rotas[veic % maxVeic];
					rotas[veic % maxVeic] = tmp;

					vertices[v] = vertices[numVerticesDescobertos];
					k = (veic % maxVeic) + 1;
					--numVerticesDescobertos;
					inseriuVertice = true;
					break;
				}
			}

			if ( !inseriuVertice ) break; //Essa tentativa nao deu certo, reinicia os vetores e ponteiros e passa para a proxima tentativa
		}

		if (numVerticesDescobertos == 0) //Significa que obteve as rotas viaveis
		{
			break;
		}
		else //Limpa as rotas sem sucesso obtidas, para iniciar uma nova iteracao procurando por outras rotas
		{
			for (int k = 0; k < maxVeic; ++k)
			{
				tmp = rotas[k];
				while(tmp != NULL)
				{
					aux = tmp;
					tmp = aux->prox;
					delete aux;
				}
			}
		}
	}

	if ( numVerticesDescobertos == 0 )
	{
		ptrNoListaEnc aux;
		vector < Rota* > rotasGamma;
		for ( int k = 0; k < maxVeic; ++k )
		{
			Rota* r = new Rota();
			r->inserirVerticeFim( ( k  % nSatellites ) + 1 );
			while ( rotas[k] != NULL )
			{
				r->inserirVerticeFim(rotas[k]->vertice + nSatellites);
				aux = rotas[k];
				rotas[k] = rotas[k]->prox;
				delete aux;
			}
			r->inserirVerticeFim( ( k  % nSatellites ) + 1 );
			r->setCusto(g);
			rotasGamma.push_back(r);
		}
		delete [] rotas;
		delete [] vertices;
		return rotasGamma;
	}
	else
	{
		printf ( "Não encontrou rotas viaveis para k(L2) = %d\n", maxVeic );
		exit(0);
	}
}
