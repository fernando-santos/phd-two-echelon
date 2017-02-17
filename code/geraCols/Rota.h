#ifndef ROTA_H_
#define ROTA_H_

#include "Grafo.h"
#include <stdio.h>
#include <vector>
using namespace std;


class Rota{
	private:
		vector< int > vertices;
		double custo;
		double custoReduzido;
		int numApontadores;
		
	public:
		Rota();
		~Rota();
		bool ciclo;
		//true: existe um ciclo na rota e este nao tem corte;
		//false: nao existe ciclo ou existe e o corte ja foi incluido
		
		void setVertices(vector<int>);
		vector<int> getVertices();
		
		void setCusto(Grafo*);
		double getCusto();

		void setCustoReduzido(double);
		void setCustoReduzido(Grafo*);
		double getCustoReduzido();

		void incrNumApontadores();
		bool decrNumApontadores();
		int getNumApontadores();

		void inserirVerticeFim(int);
		void inserirVerticeInicio(int);
		void imprimir();
};

#endif /*ROTA_H_*/
