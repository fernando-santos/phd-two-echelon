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
#include "combsep.h"
#include "mstarsep.h"


class ModeloCplex{
	friend class NoArvore;
	friend class ModeloHeuristica;

	private:
		static int* qRotasL2;
		static int maxVeicL2, maxVeicL1, maxVeicSat, nVertices, nCustomers, nSatellites, quantRotasL1;
		static double limitePrimal, thresholdCapacityCuts, thresholdCombCuts, thresholdMultiStarCuts, thresholdCycleCuts;

		static vector<short int*>* A_qr;
		static vector<Rota*>* ptrRotasL2;

		double *edgesX;
		double **grafoCycleCuts, *yCycleCuts;
		CnstrMgrPointer newCuts, oldCuts;
		int *demands, *edgesHead, *edgesTail, QMin;
		vector < char** > combCutsLHS;
		vector < float** > capacityCutsLHS;
		vector < short int** > multiStarCutsLHS;
		vector < char* >* cycleCutsLHS;
		vector < char* >* routeCutsLHS;
		char* auxCycleCuts;

		IloEnv env;
		IloModel model;
		IloCplex cplex;
		IloModel modelCGH;
		IloCplex cplexCGH;

		//Variaveis de decisao da raiz e arvore BP
		IloNumVarArray lambda;
		IloArray < IloNumVarArray > gamma;
		IloArray < IloNumVarArray > delta;
		IloArray < IloNumVarArray > gamma_no;
		IloNumVarArray artificiais;

		IloIntVarArray lambdaCGH;
		IloArray < IloIntVarArray > deltaCGH;
		IloArray < IloIntVarArray > gammaCGH;
		IloArray < IloIntVarArray > gamma_noCGH;

		//Funcao objetivo e Restricoes Primal
		IloObjective objCost;
		IloRange constraint2;
		IloRangeArray constraints3;
		IloRangeArray constraints4;
		IloRangeArray constraints5;
		IloRangeArray constraints6;
		IloRange constraint7;

		IloRangeArray constraintsArvoreL1;
		IloRangeArray constraintsArvoreL2;
		IloArray < IloRangeArray > constraintsArvoreVeic;

		IloRangeArray combCuts;
		IloRangeArray capacityCuts;
		IloRangeArray multiStarCuts;
		IloArray < IloRangeArray > cycleCuts;
		IloArray < IloRangeArray > routeCuts;

		IloObjective objCostCGH;
		IloRange constraint2CGH;
		IloRangeArray constraints3CGH;
		IloRangeArray constraints4CGH;
		IloRangeArray constraints5CGH;
		IloRangeArray constraints6CGH;
		IloRange constraintCGH7;

	public:
		char sepCuts[5];
		ModeloCplex( Grafo*, vector < Rota* >, vector < Rota* > );
		~ModeloCplex();

		//Armazena a matriz A_{qr} no modelo de uma maneira um pouco mais eficiente
		void setA_qr( vector < Rota* > );
		double solveMaster( bool processaBase = false );
		void exportModel( const char* );
		
		//inicializa as variaveis a serem utilizadas no modelo Primal
		void initVars( Grafo* );

		//define os coeficientes das variaveis na funcao objetivo e nas restricoes do modelo
		void setObjectiveFunction( vector < Rota* > );
		void setConstraints1( Grafo* );
		void setConstraint2();
		void setConstraints3();
		void setConstraints4();
		void setConstraints5( Grafo* );
		void setConstraints6( Grafo*, vector < Rota* > );

		//atualiza custos Duais no grafo
		void updateDualCosts( Grafo*, int );
		
		//retorna o valor constante a ser subtraida, para saber se uma determinada rota tem ou nao custo reduzido negativo
		double getValorSubtrair( Grafo*, int );

		//insere uma coluna no primal relacionada a rota dos fornecedores/consumidores
		void insertColumn( Grafo*, Rota*, int );

		//separacao dos cortes
		bool addCombCuts( Grafo*, int*, vector < Rota* >* );
		bool addCapacityCuts( Grafo*, int*, vector < Rota* >* );
		bool addMultiStarCuts( Grafo*, int*, vector < Rota* >* );
		bool addCycleCuts( Grafo*, int*, vector<short int*>*, vector < Rota* >* );
		bool addRouteCuts( Grafo*, int*, vector<short int*>*, vector < Rota* >* );

		//montagem dos grafos usados na separacao
		void buildGraphSolution( int&, int*, int*, double*, int*, vector < Rota* >* );
		void buildGraphCycleCuts( int s, IloArray <IloNumArray>&,  IloArray <IloNumArray>&, int*, vector < Rota* >* );
		int mostViolatedCapacityCut(int*, vector < Rota* >*, double&, int);

		static double getLimitePrimal();
		static void setLimitePrimal( double );

		void processaBase( int, int* qRotas_no = NULL, vector < short int* >* a_ir_no = NULL, vector < Rota* >* ptrRotas_no = NULL );

		void imprimeRotasBasicas( vector < Rota* >, Grafo* );

		int getNumCuts();
};

#endif
