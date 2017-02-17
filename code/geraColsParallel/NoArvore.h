#ifndef NOARVOREARCOS_H_
#define NOARVOREARCOS_H_
#include "ModeloCplex.h"
#include <algorithm>

typedef class NoArvore* ptrNo;
class NoArvore{
	private:
		static Grafo* G;
		static ModeloCplex* mCplex;
		static int nSatellites, nVertices;
		static int totalNosCriados, totalNosAtivos;
		static vector < vector < int > > listaLambdaInteiros;
		static ptrNo listaNosProcessados;

		int indiceNo;
		float limiteDual;
		bool lambdaInteiro;
		int alturaNaArvore;
		int *qRotas_no;
		vector<short int*>* a_ir_no;
		vector<Rota*>* ptrRotas_no;
		vector<short int*> arcosBranching; //short int*: ponteiro para um vetor de 4 posicoes: (i, j, s, {0,1}) -> arco fixo (i,j), s = satellite, x_{ij} = 0 ou x_{ij} = 1
		ptrNo prox;

	public:
		~NoArvore();
		NoArvore(ModeloCplex*, Grafo*, vector<ptrNo>&, double);
		NoArvore(ptrNo, short int*);

		bool alcancaRaiz_no(char);
		void setRestricoes_no();
		void setVariaveisDecisao_no();
		double executaBranching(vector<ptrNo>&, char);
		bool defineVariavelBranching(int&, int&, int&);
		void atualizaCustosDuais(int);
		void inserirColuna_no(Rota*, int);
		void imprimir();

		int getIndiceNo();
		void setProx(ptrNo);
		ptrNo getProx();
		double getLimiteDual();
		double runCplexInt();
		static int getTotalNosCriados();
		static int getTotalNosAtivos();
		static void addNoProcessado(ptrNo);
		static ptrNo getNoProcessado();
};

#endif
