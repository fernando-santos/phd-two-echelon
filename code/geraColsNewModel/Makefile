OPT = main.o
GRAFO_OBJ = Grafo.o
ROTAS_OBJ = geraRotas.o Rota.o
CPLEX_OBJ = ModeloCplex.o
ARVORE_OBJ = NoArvore.o ListaNosAtivos.o
SUBPROBLEMA_OBJ = NE2C.o GRASP.o Elementar.o NaoElementar.o ModeloBC.o Fluxo.o dinic.o
CUTS = basegrph.o binpack.o blocks.o brnching.o capsep.o cnstrmgr.o combsep.o compcuts.o compress.o cutbase.o fcapfix.o fcisep.o fcits.o glmsep.o grsearch.o hpmstar.o htoursep.o intap.o memmod.o mstarsep.o mxf.o newhtour.o sort.o strcomb.o strngcmp.o twomatch.o
ALL_OBJS = $(OPT) $(GRAFO_OBJ) $(ROTAS_OBJ) $(CPLEX_OBJ) $(ARVORE_OBJ) $(SUBPROBLEMA_OBJ) $(CUTS)

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

basegrph.o: basegrph.cpp
	$(CCC) $(CCFLAGS) -c -g basegrph.cpp

binpack.o: binpack.cpp
	$(CCC) $(CCFLAGS) -c -g binpack.cpp

blocks.o: blocks.cpp
	$(CCC) $(CCFLAGS) -c -g blocks.cpp

brnching.o: brnching.cpp
	$(CCC) $(CCFLAGS) -c -g brnching.cpp

capsep.o: capsep.cpp
	$(CCC) $(CCFLAGS) -c -g capsep.cpp

cnstrmgr.o: cnstrmgr.cpp
	$(CCC) $(CCFLAGS) -c -g cnstrmgr.cpp

combsep.o: combsep.cpp
	$(CCC) $(CCFLAGS) -c -g combsep.cpp

compcuts.o: compcuts.cpp
	$(CCC) $(CCFLAGS) -c -g compcuts.cpp

compress.o: compress.cpp
	$(CCC) $(CCFLAGS) -c -g compress.cpp

cutbase.o: cutbase.cpp
	$(CCC) $(CCFLAGS) -c -g cutbase.cpp

fcapfix.o: fcapfix.cpp
	$(CCC) $(CCFLAGS) -c -g fcapfix.cpp

fcisep.o: fcisep.cpp
	$(CCC) $(CCFLAGS) -c -g fcisep.cpp

fcits.o: fcits.cpp
	$(CCC) $(CCFLAGS) -c -g fcits.cpp
 
glmsep.o: glmsep.cpp
	$(CCC) $(CCFLAGS) -c -g glmsep.cpp

grsearch.o: grsearch.cpp
	$(CCC) $(CCFLAGS) -c -g grsearch.cpp

hpmstar.o: hpmstar.cpp
	$(CCC) $(CCFLAGS) -c -g hpmstar.cpp

htoursep.o: htoursep.cpp
	$(CCC) $(CCFLAGS) -c -g htoursep.cpp

intap.o: intap.cpp
	$(CCC) $(CCFLAGS) -c -g intap.cpp

memmod.o: memmod.cpp
	$(CCC) $(CCFLAGS) -c -g memmod.cpp

mstarsep.o: mstarsep.cpp
	$(CCC) $(CCFLAGS) -c -g mstarsep.cpp

mxf.o: mxf.cpp
	$(CCC) $(CCFLAGS) -c -g mxf.cpp

newhtour.o: newhtour.cpp
	$(CCC) $(CCFLAGS) -c -g newhtour.cpp

sort.o: sort.cpp
	$(CCC) $(CCFLAGS) -c -g sort.cpp

strcomb.o: strcomb.cpp
	$(CCC) $(CCFLAGS) -c -g strcomb.cpp

strngcmp.o: strngcmp.cpp
	$(CCC) $(CCFLAGS) -c -g strngcmp.cpp

twomatch.o: twomatch.cpp
	$(CCC) $(CCFLAGS) -c -g twomatch.cpp
 
main.o: main.cpp
	$(CCC) $(CCFLAGS) -c -g main.cpp

clean:
	rm *.o -f
	rm *~ -f
	rm geraColsVRP2E.e -f
