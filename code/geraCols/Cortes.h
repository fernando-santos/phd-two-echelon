#ifndef CORTES_H_
#define CORTES_H_
#include "NoArvore.h"
#include "ModeloCplex.h"

struct conjuntoS
{
	short int tam;
	short int* vertices;
	float violacao;
};

class Cortes{
	private:
		int nSat, nVert;

	public:
		Cortes( ModeloCplex* );
		int capacityRoundedCuts(); 
		void inserirConjuntoS( vector< conjuntoS >&, conjuntoS& );
};

#endif
