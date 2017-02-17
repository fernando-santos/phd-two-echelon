#ifndef MODELOCPLEX_H_
#define MODELOCPLEX_H_

#include <string.h>
#include <ilcplex/ilocplex.h>
#include "Rota.h"
#include "Grafo.h"
#include "ModeloBC.h"
#include "Elementar.h"
#include "NaoElementar.h"
#include "GRASP.h"
#include "NE2C.h"
#include "capsep.h"


class ModeloCplex{
	friend class Cortes;
	friend class NoArvore;
	friend class ModeloHeuristica;

	private:
		int* qRotasL2;
		int maxVeicL2;
		int maxVeicL1;
		int maxVeicSat;

		int nSatellites;
		int nCustomers;
		int nVertices;

		vector<short int*>* A_qr;
		vector<Rota*>* ptrRotasL2;
		static double limitePrimal;

		double *edgesX;
		CnstrMgrPointer newCuts, oldCuts;
		int *demands, *edgesHead, *edgesTail;
		vector< vector< int > > capacityCutsLHS;

		IloEnv env;
		IloModel model;
		IloCplex cplex;
		IloModel modelCGH;
		IloCplex cplexCGH;

		//Variaveis de decisao da raiz e arvore BP
		IloArray < IloNumVarArray > lambda;
		IloArray < IloNumVarArray > gamma;
		IloArray < IloNumVarArray > delta;
		IloArray < IloNumVarArray > gamma_no;
		IloNumVarArray artificiais;

		IloArray < IloIntVarArray > lambdaCGH;
		IloArray < IloIntVarArray > deltaCGH;
		IloArray < IloIntVarArray > gammaCGH;
		IloArray < IloIntVarArray > gamma_noCGH;

		//Funcao objetivo e Restricoes Primal
		IloObjective objCost;
		IloRangeArray constraints1;
		IloRange constraint2;
		IloRangeArray constraints3;
		IloRangeArray constraints4;
		IloRangeArray constraints5;
		IloRangeArray constraints6;
		IloRangeArray constraints7;
		IloRangeArray constraintsArvoreL1;
		IloRangeArray constraintsArvoreL2;
		IloRangeArray capacityCuts;

		IloObjective objCostCGH;
		IloRangeArray constraints1CGH;
		IloRange constraint2CGH;
		IloRangeArray constraints3CGH;
		IloRangeArray constraints4CGH;
		IloRangeArray constraints5CGH;
		IloRangeArray constraints6CGH;
		IloRangeArray constraints7CGH;

	public:
		int sepCuts;
		ModeloCplex( Grafo*, vector <Rota*>, vector <Rota*>, int );
		~ModeloCplex();

		//Armazena a matriz A_{qr} no modelo de uma maneira um pouco mais eficiente
		void setA_qr( vector < Rota* > );

		double solveMaster();

		//exporta o modelo para o arquivo passado como parametro
		void exportModel( const char* );
		
		//inicializa as variaveis a serem utilizadas no modelo Primal
		void initVars(Grafo*, vector < Rota* > );

		//define os coeficientes das variaveis na funcao objetivo e nas restricoes do modelo
		void setObjectiveFunction( vector < Rota* > );
		void setConstraints1( Grafo*, vector < Rota* > );
		void setConstraint2();
		void setConstraints3();
		void setConstraints4();
		void setConstraints5( Grafo* );
		void setConstraints6( Grafo*, vector < Rota* > );
		void setConstraints7( Grafo* );

		//atualiza custos Duais no grafo
		void updateDualCosts( Grafo*, int );
		
		//retorna o valor constante a ser subtraida, para saber se uma determinada rota tem ou nao custo reduzido negativo
		double getValorSubtrair( Grafo*, int );

		//insere uma coluna no primal relacionada a rota dos fornecedores/consumidores
		void insertColumn( Grafo*, Rota*, int );
		void buildGraphSolution( int&, int*, int*, double*, int*, vector < Rota* >* );
		bool addCapacityCuts( Grafo*, int*, vector < Rota* >* );

		static double getLimitePrimal();
		static void setLimitePrimal( double );
		void imprimeRotasBasicas( vector < Rota* >, Grafo* );
};

#endif
