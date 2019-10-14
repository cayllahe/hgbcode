SEGSURF = .
# Version LINUX
  XLIB = -L/usr/X11R6/lib -lX11 -lXext
  XINCL = /usr/include/X11R6
  CC = gcc
  CCFLAGS = -O3 -Wall -ggdb -g -DPC -DUNIXIO -lm
  ODIR = ./linux/obj
  BDIR = ./linux/bin
  LIBS = -lm 
  
SAL=\
$(BDIR)/areaFilter

OBJ = $(ODIR)

IDIR = ./include
CDIR = ./src/com
LDIR = ./src/lib

all:  $(SAL) 

clean:	
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

$(BDIR)/areaFilter:	$(CDIR)/areaFilter.c $(IDIR)/mcweightgraph.h $(IDIR)/lhierarchie.h $(IDIR)/llca.h  $(IDIR)/lcomptree.h $(ODIR)/mcweightgraph.o $(ODIR)/lhierarchie.o $(ODIR)/llca.o $(ODIR)/lcomptree.o $(ODIR)/mcimage.o $(ODIR)/fib.o $(ODIR)/mclifo.o $(ODIR)/lwatershedcut.o $(ODIR)/mcunionfind.o $(ODIR)/lhierarchie.o $(ODIR)/llca.o $(ODIR)/mcsort.o 
	$(CC) $(CCFLAGS)  -I$(IDIR) $(CDIR)/areaFilter.c $(ODIR)/fib.o $(ODIR)/mcweightgraph.o $(ODIR)/mclifo.o  $(ODIR)/lwatershedcut.o  $(ODIR)/mcunionfind.o   $(ODIR)/lhierarchie.o $(ODIR)/llca.o $(ODIR)/mcimage.o $(ODIR)/lcomptree.o $(ODIR)/mcsort.o $(LIBS) -o $(BDIR)/areaFilter

# ===============================================================
# LIBRAIRIE
# ===============================================================

$(ODIR)/base.o:	$(LDIR)/base.c $(IDIR)/base.h
	$(CC) -c $(CCFLAGS) -I$(IDIR) $(LDIR)/base.c -o $(ODIR)/base.o

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
	$(CC) -c $(LDIR)/lhisto.c -I$(IDIR) $(LIBS) -o $(ODIR)/lhisto.o

$(ODIR)/bptio.o: $(LDIR)/bptio.c $(IDIR)/bptio.h $(IDIR)/base.h 
	$(CC) -c $(LDIR)/bptio.c -I$(IDIR) $(LIBS) -o $(ODIR)/bptio.o