#ifndef GERAROTAS_H_
#define GERAROTAS_H_

#include "Grafo.h"
#include "Rota.h"

vector < vector< int > > arranjoRotasLambda( int*, int );
vector <Rota* > geraRotasLambda( Grafo* );
vector <Rota* > geraRotasGamma( Grafo* );
vector <Rota* > geraRotasGamma( Grafo*, int );

typedef struct strNoListaEnc* ptrNoListaEnc;
typedef struct strNoListaEnc{
	int vertice;
	ptrNoListaEnc prox;
}NoListaEnc;

#endif
