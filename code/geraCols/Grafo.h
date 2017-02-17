#ifndef GRAFO_H_
#define GRAFO_H_

#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <math.h>
#include "Vertice.h"

class Grafo{
	private:
		Vertice** vert; //um grafo eh composto por um vetor de ponteiros para vertices (nao eh matriz, eh vetor de ponteiros)

		//quantos vertices satellites tem o grafo
		int numSatellites;
		
		//quantos customers tem o grafo
		int numCustomers;
		
		//numSatellites + numCustomers + 1 (deposito {0})
		int numVertices;

		//limite de veiculos para os niveis 1 e 2
		int maxVeiculosL1, maxVeiculosL2;
		int maxVeiculosSatellites;

		//capacidade dos veiculos para os niveis 1 e 2
		int capacVeiculosL1, capacVeiculosL2;

	public:
		int getNumSatellites();
		int getNumCustomers();
		int getNumVertices();
		int getMaxVeiculosL1();
		int getMaxVeiculosL2();
		int getCapacVeiculosL1();
		int getCapacVeiculosL2();
		int getCargaVertice(int);
		int getMaxVeiculosSatellites();
		double getCustoAresta(int, int);
		double getCustoArestaDual(int, int);
		double getCustoVerticeDual(int);
		void setCustoVerticeDual(int, double);
		void setCustoArestasDual(int);
		void setCustoArestaDual(int, int, double);

		void imprimeGrafo();
		void imprimeGrafoDual();

		Grafo(char*);		
		~Grafo();
};

#endif
