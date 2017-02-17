#ifndef NE2C_H_
#define NE2C_H_

#include "Grafo.h"
#include "Rota.h"

typedef struct strLabelNE2C *ptrLabelNE2C;
typedef struct strLabelNE2C{
	int ultimoVertice;
	int verticeAntecessor;
	int cargaAcumulada;
	double custoDual;
	ptrLabelNE2C prox;
	ptrLabelNE2C pai;

	strLabelNE2C(int ult, int antec, int carga, double custo, ptrLabelNE2C p)
				: ultimoVertice(ult), verticeAntecessor(antec), cargaAcumulada(carga), custoDual(custo), prox(p){}
			
	void imprimeLabel()
	{
		printf("c(%f), q(%d), u(%d), a(%d), t(%p) p(%p)\n", custoDual, cargaAcumulada, ultimoVertice, verticeAntecessor, this, pai);
	}

} LabelNE2C;

struct vLabelNE2C{
	ptrLabelNE2C cabeca;
	ptrLabelNE2C posAtual;
	ptrLabelNE2C calda;
};

class NE2C{
	private:
		vLabelNE2C* vetorLabels;
		ptrLabelNE2C menorLabelObtido, lixo;
		double menorCustoAtual, menorCustoObtido, valorSubtrair;
		int numVertices, numSatellites;

	public:
		NE2C(Grafo*, double);
		~NE2C();

		void calculaCaminho(Grafo*);
		bool verificaLabelDominadoEDomina(int, int, double);
		vector <Rota*> getRotaCustoMinimo(Grafo*, int, double);
};

#endif 

