#ifndef NAOELEMENTAR_H_
#define NAOELEMENTAR_H_

#include "Grafo.h"
#include "Rota.h"

struct CabecalhoLinha{
	int indice;
	int pesoRacionalizado;
};

struct Celula{
	int melhorAntec;
	double custoMelhorAntec;
	int segMelhorAntec;
	double custoSegMelhorAntec;
};

class NaoElementar{
	private:
		Celula** matriz;
		CabecalhoLinha *linhas;
		int nLinhas;
		int nColunas;
		int nSatellites;
		int mdc;

	public:
		NaoElementar( Grafo* );
		~NaoElementar();
		
		//Metodos comum entre o Caminho Nao-Elementar e Nao-Elementar sem Ciclo de tamanho 2
		int calculaMDC(Grafo*);
		void imprimirCabecalhoLinhas();
		void imprimirMatriz(bool);
		
		//Metodos relacionados ao calculo do caminho Nao-Elementar com quaisquer ciclos
		void calculaMatriz(Grafo*);
		float getMelhorAntecessor(Grafo*, int, int, int&);
		vector<Rota*> getRotaCustoMinimo(Grafo*, int, float);

		//Metodos relacionados ao caminho Nao-Elementar sem ciclos de tamanho 2
		void calculaMatrizCiclo2(Grafo*);
		double getMelhorSegMelhorAntecessor(Grafo*, int, int, int&, double&, int&);
		vector < Rota* > getRotaCustoMinimoCiclo2(Grafo*, int, double);
	};

#endif
