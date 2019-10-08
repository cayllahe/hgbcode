SEGSURF = .
# Version LINUX
  XLIB = -L/usr/X11R6/lib -lX11 -lXext
  XINCL = /usr/include/X11R6
  CC = cc
  CCFLAGS = -g -DPC -DUNIXIO -lm
  ODIR = ./linux/obj
  BDIR = ./linux/bin
  LIBS = -lm 
  
SAL=\
$(BDIR)/uprooting \
$(BDIR)/ultramopen\
$(BDIR)/extinctionvalues\
$(BDIR)/pgm2graph\
$(BDIR)/graph2pgm\
$(BDIR)/graphmeshconverter \
$(BDIR)/float2byte

OBJ = $(ODIR)

IDIR = ./include
CDIR = ./src/com
LDIR = ./src/lib

all:  $(SAL) 

clean:	
	rm -f ./makefile
	rm -f ./linux/bin/*
	rm -f ./linux/bin2/*
	rm -f ./linux/obj/*
	rm -f ./hpux/bin/*
	rm -f ./hpux/obj/*
	rm -f $(CDIR)/*~
	rm -f $(LDIR)/*~
	rm -f $(IDIR)/*~
	rm -f ./*~
	rm -f ./*.o




# ===============================================================
# EXECUTABLES
# ===============================================================



$(BDIR)/graphmeshconverter:	$(CDIR)/graphmeshconverter.c $(IDIR)/lgraphmesh.h $(IDIR)/mcweightgraph.h    $(IDIR)/mcmesh.h $(IDIR)/mciomesh.h    $(ODIR)/lgraphmesh.o $(ODIR)/mcweightgraph.o $(ODIR)/mcmesh.o $(ODIR)/mciomesh.o $(ODIR)/mcrbt1.o $(ODIR)/mcprobas.o $(ODIR)/mcgeo.o 
	$(CC) $(CCFLAGS) -I$(IDIR) $(CDIR)/graphmeshconverter.c  $(ODIR)/lgraphmesh.o $(ODIR)/mcweightgraph.o $(ODIR)/mcmesh.o $(ODIR)/mciomesh.o $(ODIR)/mcrbt1.o $(ODIR)/mcprobas.o $(ODIR)/mcgeo.o  -o $(BDIR)/graphmeshconverter $(LIBS)

# $(BDIR)/uprooting:	$(CDIR)/uprooting.c $(IDIR)/mcweightgraph.h $(IDIR)/lwatershedcut.h $(IDIR)/fib.h  $(IDIR)/mcunionfind.h $(IDIR)/llca.h $(IDIR)/lhierarchie.h  $(IDIR)/mcsort.h $(ODIR)/mcweightgraph.o $(ODIR)/mcunionfind.o $(ODIR)/fib.o $(ODIR)/lhierarchie.o $(ODIR)/lcomptree.o $(ODIR)/lminimgraph.o $(ODIR)/mclifo.o  $(ODIR)/llca.o $(ODIR)/mcsort.o
# 	$(CC) $(CCFLAGS) -I$(IDIR) $(CDIR)/uprooting.c $(ODIR)/fib.o $(ODIR)/mcweightgraph.o $(ODIR)/lcomptree.o $(ODIR)/mcunionfind.o $(ODIR)/mclifo.o   $(ODIR)/lhierarchie.o $(ODIR)/llca.o $(ODIR)/lminimgraph.o $(ODIR)/mcsort.o -o $(BDIR)/uprooting $(LIBS)


$(BDIR)/uprooting:	$(CDIR)/uprooting.c $(IDIR)/mcweightgraph.h $(IDIR)/lwatershedcut.h $(IDIR)/fib.h  $(IDIR)/mcunionfind.h $(IDIR)/llca.h $(IDIR)/lhierarchie.h  $(IDIR)/mcsort.h $(ODIR)/mcweightgraph.o $(ODIR)/mcunionfind.o $(ODIR)/fib.o $(ODIR)/lhierarchie.o $(ODIR)/lcomptree.o  $(ODIR)/mclifo.o  $(ODIR)/llca.o $(ODIR)/mcsort.o $(ODIR)/lwatershedcut.o
	$(CC) $(CCFLAGS) -I$(IDIR) $(CDIR)/uprooting.c $(ODIR)/fib.o $(ODIR)/mcweightgraph.o $(ODIR)/lcomptree.o $(ODIR)/mcunionfind.o $(ODIR)/mclifo.o   $(ODIR)/lhierarchie.o $(ODIR)/llca.o  $(ODIR)/mcsort.o $(ODIR)/lwatershedcut.o -o $(BDIR)/uprooting $(LIBS)


$(BDIR)/ultramopen:	$(CDIR)/ultramopen.c $(IDIR)/mcweightgraph.h $(IDIR)/lhierarchie.h $(IDIR)/llca.h  $(IDIR)/lcomptree.h  $(IDIR)/mcunionfind.h $(IDIR)/mcsort.h $(ODIR)/mcweightgraph.o $(ODIR)/lhierarchie.o $(ODIR)/llca.o $(ODIR)/lcomptree.o  $(ODIR)/mcunionfind.o $(ODIR)/mcsort.o
	$(CC) $(CCFLAGS) -I$(IDIR) $(CDIR)/ultramopen.c $(ODIR)/fib.o $(ODIR)/mcweightgraph.o $(ODIR)/mclifo.o  $(ODIR)/lwatershedcut.o  $(ODIR)/mcunionfind.o  $(ODIR)/lhierarchie.o $(ODIR)/llca.o $(ODIR)/lcomptree.o $(ODIR)/mcsort.o -o $(BDIR)/ultramopen $(LIBS)

# $(BDIR)/extinctionvalues:	$(CDIR)/extinctionvalues.c $(IDIR)/mcweightgraph.h $(IDIR)/lhierarchie.h $(IDIR)/llca.h  $(IDIR)/lcomptree.h $(ODIR)/mcweightgraph.o $(ODIR)/lhierarchie.o $(ODIR)/llca.o $(ODIR)/lcomptree.o
# 	$(CC) $(CCFLAGS) -I$(IDIR) $(CDIR)/extinctionvalues.c $(ODIR)/fib.o $(ODIR)/mcweightgraph.o $(ODIR)/mclifo.o  $(ODIR)/lminimgraph.o  $(ODIR)/mcunionfind.o  $(ODIR)/lhierarchie.o $(ODIR)/llca.o $(ODIR)/lcomptree.o -o $(BDIR)/extinctionvalues $(LIBS)

$(BDIR)/extinctionvalues:	$(CDIR)/extinctionvalues.c $(IDIR)/mcweightgraph.h  $(IDIR)/lcomptree.h  $(IDIR)/mcunionfind.h $(IDIR)/mcsort.h $(ODIR)/mcweightgraph.o  $(ODIR)/lcomptree.o $(ODIR)/mcunionfind.o $(ODIR)/mcsort.o
	$(CC) $(CCFLAGS) -I$(IDIR) $(CDIR)/extinctionvalues.c $(ODIR)/mcsort.o $(ODIR)/fib.o $(ODIR)/mcweightgraph.o  $(ODIR)/mcunionfind.o   $(ODIR)/lcomptree.o -o $(BDIR)/extinctionvalues $(LIBS)


$(BDIR)/pgm2graph:	$(CDIR)/pgm2graph.c $(IDIR)/mcweightgraph.h $(IDIR)/mccodimage.h $(IDIR)/mcimage.h $(ODIR)/mcweightgraph.o $(ODIR)/mcimage.o
	$(CC) $(CCFLAGS) -I$(IDIR) $(CDIR)/pgm2graph.c $(ODIR)/mcweightgraph.o  $(ODIR)/mcimage.o  -o $(BDIR)/pgm2graph $(LIBS)

$(BDIR)/graph2pgm:	$(CDIR)/graph2pgm.c $(IDIR)/mcweightgraph.h $(IDIR)/mccodimage.h $(IDIR)/mcimage.h $(ODIR)/mcweightgraph.o $(ODIR)/mcimage.o
	$(CC) $(CCFLAGS) -I$(IDIR) $(CDIR)/graph2pgm.c $(ODIR)/mcweightgraph.o  $(ODIR)/mcimage.o  -o $(BDIR)/graph2pgm $(LIBS)

$(BDIR)/float2byte:	$(CDIR)/float2byte.c  $(IDIR)/mccodimage.h $(IDIR)/mcimage.h $(IDIR)/mcutil.h $(IDIR)/lhisto.h  $(ODIR)/mcimage.o  $(ODIR)/lhisto.o
	$(CC) $(CCFLAGS) -I$(IDIR) $(CDIR)/float2byte.c $(ODIR)/lhisto.o  $(ODIR)/mcimage.o  -o $(BDIR)/float2byte $(LIBS)
 
# ===============================================================
# LIBRAIRIE
# ===============================================================



$(ODIR)/mcgeo.o:	$(LDIR)/mcgeo.c $(IDIR)/mcgeo.h
	$(CC) -c $(CCFLAGS) -I$(IDIR) $(LDIR)/mcgeo.c -o $(ODIR)/mcgeo.o

$(ODIR)/mciomesh.o:	$(LDIR)/mciomesh.c $(IDIR)/mciomesh.h
	$(CC) -c $(CCFLAGS) -I$(IDIR) $(LDIR)/mciomesh.c -o $(ODIR)/mciomesh.o

$(ODIR)/mclifo.o:	$(LDIR)/mclifo.c $(IDIR)/mclifo.h
	$(CC) -c $(CCFLAGS) -I$(IDIR) $(LDIR)/mclifo.c -o $(ODIR)/mclifo.o

$(ODIR)/mcmesh.o:	$(LDIR)/mcmesh.c $(IDIR)/mcmesh.h $(IDIR)/mcprobas.h
	$(CC) -c $(CCFLAGS) -I$(IDIR) $(LDIR)/mcmesh.c -o $(ODIR)/mcmesh.o

$(ODIR)/mcprobas.o:	$(LDIR)/mcprobas.c $(IDIR)/mcprobas.h
	$(CC) -c $(CCFLAGS) -I$(IDIR) $(LDIR)/mcprobas.c -o $(ODIR)/mcprobas.o

$(ODIR)/mcrbt1.o:	$(LDIR)/mcrbt1.c $(IDIR)/mcrbt1.h
	$(CC) -c $(CCFLAGS) -I$(IDIR) $(LDIR)/mcrbt1.c -o $(ODIR)/mcrbt1.o

$(ODIR)/mcsort.o:	$(LDIR)/mcsort.c $(IDIR)/mcsort.h
	$(CC) -c $(CCFLAGS) -I$(IDIR) $(LDIR)/mcsort.c -o $(ODIR)/mcsort.o

$(ODIR)/mcunionfind.o:	$(LDIR)/mcunionfind.c $(IDIR)/mcunionfind.h
	$(CC) -c $(CCFLAGS) -I$(IDIR) $(LDIR)/mcunionfind.c -o $(ODIR)/mcunionfind.o		

$(ODIR)/mcweightgraph.o:	$(LDIR)/mcweightgraph.c $(IDIR)/mcweightgraph.h $(IDIR)/mcsort.h
	$(CC) -c $(CCFLAGS) -I$(IDIR) $(LDIR)/mcweightgraph.c -o $(ODIR)/mcweightgraph.o

$(ODIR)/mcimage.o:	$(LDIR)/mcimage.c $(IDIR)/mcimage.h $(IDIR)/mccodimage.h $(IDIR)/mcutil.h
	$(CC) -c $(CCFLAGS) -I$(IDIR) $(LDIR)/mcimage.c -o $(ODIR)/mcimage.o

$(ODIR)/lwatershedcut.o:	$(LDIR)/lwatershedcut.c  $(IDIR)/mclifo.h $(IDIR)/mcweightgraph.h $(IDIR)/lwatershedcut.h 
	$(CC) -c $(CCFLAGS) -I$(IDIR) $(LDIR)/lwatershedcut.c -o $(ODIR)/lwatershedcut.o



$(ODIR)/lgraphmesh.o:	$(LDIR)/lgraphmesh.c $(IDIR)/mcmesh.h $(IDIR)/mciomesh.h  $(IDIR)/mcweightgraph.h  
	$(CC) -c $(CCFLAGS) -I$(IDIR) $(LDIR)/lgraphmesh.c -o $(ODIR)/lgraphmesh.o 


$(ODIR)/fib.o:	$(LDIR)/fib.c $(IDIR)/fibpriv.h  $(IDIR)/fib.h 
	$(CC) -c $(CCFLAGS) -I$(IDIR) $(LDIR)/fib.c -o $(ODIR)/fib.o 

$(ODIR)/lhierarchie.o:	$(LDIR)/lhierarchie.c $(IDIR)/lhierarchie.h  $(IDIR)/mcweightgraph.h $(IDIR)/mcunionfind.h $(IDIR)/fib.h $(IDIR)/llca.h
	$(CC) -c $(CCFLAGS) -I$(IDIR) $(LDIR)/lhierarchie.c -o $(ODIR)/lhierarchie.o 

$(ODIR)/llca.o:	$(LDIR)/llca.c $(IDIR)/llca.h 
	$(CC) -c $(LDIR)/llca.c -I$(IDIR)  -o $(ODIR)/llca.o

$(ODIR)/lcomptree.o: $(LDIR)/lcomptree.c $(IDIR)/lcomptree.h $(IDIR)/mcunionfind.h $(IDIR)/mcweightgraph.h
	$(CC) -c $(LDIR)/lcomptree.c -I$(IDIR) -o $(ODIR)/lcomptree.o

$(ODIR)/lhisto.o: $(LDIR)/lhisto.c $(IDIR)/lhisto.h $(IDIR)/mcimage.h $(IDIR)/mccodimage.h $(IDIR)/mcutil.h
	$(CC) -c $(LDIR)/lhisto.c -I$(IDIR) -o $(ODIR)/lhisto.o

