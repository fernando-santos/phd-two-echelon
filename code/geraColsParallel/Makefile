OPT = main.o
GRAFO_OBJ = Grafo.o
ROTAS_OBJ = geraRotas.o Rota.o
CPLEX_OBJ = ModeloCplex.o
ARVORE_OBJ = NoArvore.o ListaNosAtivos.o
SUBPROBLEMA_OBJ = NE2C.o GRASP.o Elementar.o NaoElementar.o ModeloBC.o Fluxo.o dinic.o
ALL_OBJS = $(OPT) $(GRAFO_OBJ) $(ROTAS_OBJ) $(CPLEX_OBJ) $(ARVORE_OBJ) $(SUBPROBLEMA_OBJ)

SYSTEM     = x86-64_debian4.0_4.1
LIBFORMAT  = static_pic
CPLEXDIR      = /home/lapo/cplex/cplex121_0.x64linux/cplex121
CONCERTDIR    = /home/lapo/cplex/cplex121_0.x64linux/concert29
CCC = g++
CCOPT = -m64 -O -fPIC -fexceptions -DNDEBUG -DIL_STD
CPLEXBINDIR   = $(CPLEXDIR)/bin/$(BINDIST)
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CCLNFLAGS = -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert -lm -pthread
CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include
CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) 

geraColsVRP2E: $(ALL_OBJS)
	$(CCC) -Wall -W -g -DNOTE $(ALL_OBJS) -o geraColsVRP2E.e $(CCLNFLAGS)

Grafo.o: Grafo.cpp
	$(CCC) $(CCFLAGS) -c -g Grafo.cpp

geraRotas.o: geraRotas.cpp
	$(CCC) $(CCFLAGS) -c -g geraRotas.cpp

Rota.o: Rota.cpp
	$(CCC) $(CCFLAGS) -c -g Rota.cpp

ModeloCplex.o: ModeloCplex.cpp
	$(CCC) $(CCFLAGS) -c -g ModeloCplex.cpp

NoArvore.o: NoArvore.cpp
	$(CCC) $(CCFLAGS) -c -g NoArvore.cpp

ListaNosAtivos.o: ListaNosAtivos.cpp
	$(CCC) $(CCFLAGS) -c -g ListaNosAtivos.cpp

NE2C.o: NE2C.cpp
	$(CCC) $(CCFLAGS) -c -g NE2C.cpp

GRASP.o: GRASP.cpp
	$(CCC) $(CCFLAGS) -c -g GRASP.cpp

Elementar.o: Elementar.cpp
	$(CCC) $(CCFLAGS) -c -g Elementar.cpp

NaoElementar.o: NaoElementar.cpp
	$(CCC) $(CCFLAGS) -c -g NaoElementar.cpp

ModeloBC.o: ModeloBC.cpp
	$(CCC) $(CCFLAGS) -c -g ModeloBC.cpp

Fluxo.o: Fluxo.cpp
	$(CCC) $(CCFLAGS) -c -g Fluxo.cpp

dinic.o: dinic.cpp
	$(CCC) $(CCFLAGS) -c -g dinic.cpp
 
main.o: main.cpp
	$(CCC) $(CCFLAGS) -c -g main.cpp

clean:
	rm *.o -f
	rm *~ -f
	rm geraColsVRP2E.e -f
