#include "Grafo.h"
using namespace std;

double Grafo::getCustoVerticeDual(int v){
	return vert[v]->custoDual;
}


double Grafo::getCustoAresta(int i, int j){
	return vert[i]->custoArestas[j];
}


double Grafo::getCustoArestaDual(int i, int j){
	return vert[i]->custoArestasDual[j];
}


int Grafo::getCargaVertice(int i){
	return vert[i]->carga;
}


int Grafo::getCapacVeiculosL1(){
  return capacVeiculosL1;
}


int Grafo::getCapacVeiculosL2(){
  return capacVeiculosL2;
}


int Grafo::getMaxVeiculosL1(){
	return maxVeiculosL1;
}


int Grafo::getMaxVeiculosL2(){
	return maxVeiculosL2;
}


int Grafo::getMaxVeiculosSatellites(){
	return maxVeiculosSatellites;
}


int Grafo::getNumSatellites(){
	return numSatellites;
}


int Grafo::getNumCustomers(){
	return numCustomers;
}


int Grafo::getNumVertices(){
	return numVertices;
}


void Grafo::setCustoArestaDual(int i, int j, double custo){
	vert[i]->custoArestasDual[j] += custo;
}


void Grafo::setCustoVerticeDual(int v, double custoD){
	vert[v]->custoDual = custoD;
}


void Grafo::setCustoArestasDual( int s ){
	for ( int i = 1; i <= numCustomers; ++i )
	{
		vert[0]->custoArestasDual[i] = vert[s]->custoArestas[numSatellites + i] + vert[numSatellites + i]->custoDual;
		vert[i]->custoArestasDual[0] = vert[numSatellites + i]->custoArestas[s];

		for ( int j = 1; j <= numCustomers; ++j )
		{
			if ( i != j )
			{
				vert[i]->custoArestasDual[j] = vert[numSatellites + i]->custoArestas[numSatellites + j] + vert[numSatellites + j]->custoDual;
			}
			else
			{
				vert[i]->custoArestasDual[j] = MAIS_INFINITO;
			}
		}
	}
}


Grafo::~Grafo(){
	for (int i = 0; i < numVertices; ++i){
		delete vert[i];
	}
	delete [] vert;
}


Grafo::Grafo(char* arquivoInstancia){
	int tmpI;
	double tmpF;
	
	//Abre o arquivo de instancia para ser processado e montar o grafo
	ifstream inFile(arquivoInstancia);
	inFile.seekg (0, ios::beg);

	//loop para retirar o cabecalho do arquivo
	char tipoInstancia, buffer[200];
	inFile.getline(buffer, 200, '\n');
	tipoInstancia = buffer[7];
	inFile.getline(buffer, 200, '\n');
	inFile.getline(buffer, 200, '\n');
	while( inFile.get() != ':' );
	inFile >> numVertices;
	while( inFile.get() != ':' );
	inFile >> numSatellites;
	while( inFile.get() != ':' );
	inFile >> numCustomers;

	inFile.getline(buffer, 200, '\n');
	inFile.getline(buffer, 200, '\n');
	while( inFile.get() != ':' );
	inFile >> capacVeiculosL1;
	while( inFile.get() != ':' );
	inFile >> capacVeiculosL2;
	while( inFile.get() != ':' );
	inFile >> maxVeiculosL1;
	while( inFile.get() != ':' );
	inFile >> maxVeiculosL2;
	maxVeiculosSatellites = maxVeiculosL2;
	inFile.getline(buffer, 200, '\n');
	inFile.getline(buffer, 200, '\n');

	//aloca o vetor de ponteiros
	vert = new Vertice*[numVertices];
	for (int i = 0; i < numVertices; ++i){
		vert[i] = new Vertice(numVertices); //um vertice aponta para todos os outros (dele para ele mesmo eh infinito)
	}

	if ( tipoInstancia == 'E' )
	{
		if ( numVertices <= 15)
		{
			for ( int i = 0; i < numVertices; ++i)
			{
				for ( int j = 0; j < numVertices; ++j)
				{
					inFile >> tmpF;
					if ( ! ( ( i == j ) || ( ( i == 0 ) && ( j > numSatellites ) ) || ( ( i > numSatellites ) && ( j == 0 ) ) ) )
					{
						vert[i]->existeArestas[j] = true;
						vert[i]->custoArestas[j] = tmpF;
					}
				}
			}

			inFile.getline(buffer, 200, '\n');
			inFile.getline(buffer, 200, '\n');
			for ( int i = 0; i < numVertices; ++i)
			{
				inFile >> tmpI >> tmpI;
				vert[i]->carga = tmpI;
			}
			inFile.close();
		}
		else
		{
			//Lê todos os dados da instancia e armazena no vetor para as alteracoes necessarias
			double matriz[numVertices][2];
			inFile >> tmpI >> matriz[0][0];
			inFile >> matriz[0][1];
			for (int i = numSatellites+1; i < numVertices; ++i) inFile >> tmpI >> matriz[i][0] >> matriz[i][1];
			
			inFile.getline(buffer, 200, '\n');
			inFile.getline(buffer, 200, '\n');
			for (int i = 1; i <= numSatellites; ++i) inFile >> tmpI >> matriz[i][0] >> matriz[i][1];
			
			for ( int i = 0; i < numVertices; ++i)
			{
				for ( int j = 0; j < numVertices; ++j)
				{
					if ( ! ( ( i == j ) || ( ( i == 0 ) && ( j > numSatellites ) ) || ( ( i > numSatellites ) && ( j == 0 ) ) ) )
					{
						vert[i]->existeArestas[j] = true;
						vert[i]->custoArestas[j] = sqrt (pow ((matriz[i][0] - matriz[j][0]), 2) + pow ((matriz[i][1] - matriz[j][1]), 2));
					}
				}
			}
			
			inFile.getline(buffer, 200, '\n');
			inFile.getline(buffer, 200, '\n');
			inFile.getline(buffer, 200, '\n');
			for (int i = 0; i <= numSatellites; ++i) vert[i]->carga = 0;
			for (int i = numSatellites+1; i < numVertices; ++i) inFile >> tmpI >> vert[i]->carga;
		}
	}
	else
	{
		double matriz[numVertices][2];
		for ( int i = numSatellites+1; i < numVertices; ++i )
		{
			while ( inFile.get() != '\t' );
			inFile >> matriz[i][0] >> matriz[i][1] >> vert[i]->carga >> tmpI;
		}

		for (int i = 1; i <= numSatellites; ++i)
		{
			while ( inFile.get() != '\t' );
			inFile >> matriz[i][0] >> matriz[i][1] >> maxVeiculosSatellites >> tmpI;
			vert[i]->carga = 0;
		}

		while ( inFile.get() != '\t' );
		inFile >> matriz[0][0] >> matriz[0][1];
		vert[0]->carga = 0;

		for ( int i = 0; i < numVertices; ++i)
		{
			for ( int j = 0; j < numVertices; ++j)
			{
				if ( ! ( ( i == j ) || ( ( i == 0 ) && ( j > numSatellites ) ) || ( ( i > numSatellites ) && ( j == 0 ) ) ) )
				{
					vert[i]->existeArestas[j] = true;
					vert[i]->custoArestas[j] = sqrt (pow ((matriz[i][0] - matriz[j][0]), 2) + pow ((matriz[i][1] - matriz[j][1]), 2));
				}
			}
		}
	}
}


void Grafo::imprimeGrafo(){
	for(int i = 0; i < numVertices; ++i)
	{
		printf("%02d} (%02d, %0.02f)   \t[", i, vert[i]->carga, vert[i]->custoDual);
		for (int j = 0; j < numVertices; ++j)
		{
			if (vert[i]->existeArestas[j])
			{
				printf("%02.03f, ", vert[i]->custoArestas[j]);
			}
			else
			{
				printf("- , ");
			}
		}
		printf("]\n");	
	}
	printf("capacidade = %d(L1) %d(L2)\n", capacVeiculosL1, capacVeiculosL2);
	printf("maxVeiculos = %d(L1) %d(L2)\n", maxVeiculosL1, maxVeiculosL2);
	printf("maxVeiculosSatellites = %d\n", maxVeiculosSatellites);
}


void Grafo::imprimeGrafoDual(){
	for(int i = 0; i < numVertices; ++i)
	{
		printf("%02d} (%02d) [", i, vert[i]->carga);
		for (int j = 0; j < numVertices; ++j)
		{
			printf("%02.01f, ", vert[i]->custoArestasDual[j]);
		}
		printf("]\n");	
	}
}
