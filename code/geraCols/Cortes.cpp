#include "Cortes.h"

Cortes::Cortes(ModeloCplex* mCplex, Grafo* G)
{
	nSat = mCplex->nSatellites;
	nVert = mCplex->nVertices;
	IloNumArray gammaValues( mCplex->env );
	vector < conjuntoS > conjuntosS;
	
	//armazena as demandas dos vertices em um vetor, para evitar sucessivas chamadas ao metodo G->getCargaVertice(i)
	int minVeic;
	int capVeic = G->getCapacVeic();
	int vetorCarga = new int[nVert];
	for ( int i = 0; i < nVert; ++i )
	{
		vetorCarga[i] = G->getCargaVertice(i);
	}

	//celulas matrizArcos[nVert][i] e matrizArcos[i][nVert] armazenam o fluxo que sai e entra (respectivamente) do vertice i
	double somaS;
	double** matrizArcos = new double*[nVert];
	for (int i = 0; i <= nVert; ++i)
	{
		matrizArcos[i] = new double[nVert];
		memset(matrizArcos[i], 0, nVert*sizeof(double));
	}

	//insere os valores dos arcos em cada celula da matriz para saber quais arcos e vertices estarao no grafo residual
	for (int s = 1; s <= nSat; ++s)
	{
		//primeiro busca por variaveis fracionarias entre as variaveis de rotas L2 obtidas na raiz
		mCplex->cplex.getValues(gammaValues, mCplex->gamma[s]);
		for (int r = 0; r < mCplex->qRotasL2[s]; ++r)
		{
			if ( gammaValues[r] > 0.000001 ) //gamma[s][r] > 0
			{
				verticesRota = mCplex->ptrRotasL2[s][r]->getVertices();
				tamRota = verticesRota.size()-1;
				for (int i = 0; i < tamRota; ++i)
				{
					matrizArcos[verticesRota[i]][verticesRota[i+1]] += gammaValues[r];
				}
			}
		}
	}

	//Procura por conjuntos S de cardinalidade 3 que violem a restricao do corte
	for ( int i = (nSat+1); i < numVert; ++i )
	{
		for ( int j = i+1; j < numVert; ++j )
		{
			for ( int k = j+1; k < numVert; ++k )
			{
				somaS = 0;
				for ( int v = 1; v < numVert; ++v )
				{
					if ( ( v != i ) && ( v != j ) && ( v != k ) )
					{
						somaS += matrizArcos[v][i];
						somaS += matrizArcos[v][j];
						somaS += matrizArcos[v][k];
					}
				}

				//avalia quantos veiculos sao necessarios para atender a demanda dos vertices em S
				minVeic = ceil( (float)( vetorCarga[i] + vetorCarga[j] + vetorCarga[k] ) / capVeic );
				
				if ( somaS < minVeic )
				{
					conjuntoS conjS;
					conjS.tamanho = 3;
					conjS.vertices = new short int[3];
					conjS.vertices[0] = i;
					conjS.vertices[1] = j;
					conjS.vertices[2] = k;
					conjS.violacao = minVeic-somaS;
					inserirConjuntoS( conjuntosS, conjS );
				}
			}
		}
	}

	for ( int i = (nSat+1); i < numVert; ++i )
	{
		for ( int j = i+1; j < numVert; ++j )
		{
			for ( int k = j+1; k < numVert; ++k )
			{
				for ( int l = k+1; l < numVert; ++l )
				{
					somaS = 0;
					for ( int v = 1; v < numVert; ++v )
					{
						if ( ( v != i ) && ( v != j ) && ( v != k ) && ( v != l ) )
						{
							somaS += matrizArcos[v][i];
							somaS += matrizArcos[v][j];
							somaS += matrizArcos[v][k];
							somaS += matrizArcos[v][l];
						}
					}

					//avalia quantos veiculos sao necessarios para atender a demanda dos vertices em S
					minVeic = ceil( (float)( vetorCarga[i] + vetorCarga[j] + vetorCarga[k] + vetorCarga[l] ) / capVeic );
					
					if ( somaS < minVeic )
					{
						conjuntoS conjS;
						conjS.tamanho = 4;
						conjS.vertices = new short int[4];
						conjS.vertices[0] = i;
						conjS.vertices[1] = j;
						conjS.vertices[2] = k;
						conjS.vertices[3] = l;
						conjS.violacao = minVeic-somaS;
						inserirConjuntoS( conjuntosS, conjS );
					}
				}
			}
		}
	}

	for ( int i = (nSat+1); i < numVert; ++i )
	{
		for ( int j = i+1; j < numVert; ++j )
		{
			for ( int k = j+1; k < numVert; ++k )
			{
				for ( int l = k+1; l < numVert; ++l )
				{
					for ( int m = l+1; m < numVert; ++m )
					{
						somaS = 0;
						for ( int v = 1; v < numVert; ++v )
						{
							if ( ( v != i ) && ( v != j ) && ( v != k ) && ( v != l ) && ( v != m ) )
							{
								somaS += matrizArcos[v][i];
								somaS += matrizArcos[v][j];
								somaS += matrizArcos[v][k];
								somaS += matrizArcos[v][l];
								somaS += matrizArcos[v][m];
							}
						}

						//avalia quantos veiculos sao necessarios para atender a demanda dos vertices em S
						minVeic = ceil( (float)( vetorCarga[i] + vetorCarga[j] + vetorCarga[k] + vetorCarga[l] + vetorCarga[m] ) / capVeic );
						
						if ( somaS < minVeic )
						{
							conjuntoS conjS;
							conjS.tamanho = 5;
							conjS.vertices = new short int[5];
							conjS.vertices[0] = i;
							conjS.vertices[1] = j;
							conjS.vertices[2] = k;
							conjS.vertices[3] = l;
							conjS.vertices[4] = m;
							conjS.violacao = minVeic-somaS;
							inserirConjuntoS( conjuntosS, conjS );
						}
					}
				}
			}
		}
	}
	
	inserirCortes( mCplex, conjuntosS );
}

void Cortes::inserirConjuntoS( vector< conjuntoS >& setsS, conjuntoS conjS )
{
	int tamSetsS = setsS.size();
	for ( int i = 0; i < tamSetsS; ++i )
	{
		for ( int x = 0; x < setsS[i].tamanho; ++x )
		{
			for ( int y = 0; y < conjS.tamanho; ++y )
			{
				if ( setsS[i].vertices[x] == conjS.vertices[y] )
				{
					if ( setsS[i].violacao < conjS.violacao )
					{
						setsS.erase(setsS.begin()+i);
						setsS.push_back(conjS);
					}
					return;
				}
			}
		}
	}
}

void Cortes::inserirCortes( ModeloCplex* mCplex, vector< conjuntoS >& setsS )
{
	
}
