SAL=\
$(BDIR)/ppm2graph


OBJ = $(ODIR)

IDIR = ./include
CDIR = ./src/com
LDIR = ./src/lib

all:  $(SAL) 

clean:	
	rm -f ./makefile
	rm -r ./linux/bin/*	
	rm -r ./linux/obj/*	
	rm -f $(CDIR)/*~
	rm -f $(LDIR)/*~
	rm -f $(IDIR)/*~
	rm -f ./*~
	rm -f ./*.o




# ===============================================================
# EXECUTABLES
# ===============================================================


$(BDIR)/ppm2graph:	$(CDIR)/ppm2graph.c $(IDIR)/mcweightgraph.h $(IDIR)/mccodimage.h $(IDIR)/mcimage.h $(ODIR)/mcweightgraph.o $(ODIR)/mcimage.o
	$(CC) $(CCFLAGS) -I$(IDIR) $(CDIR)/ppm2graph.c $(ODIR)/mcweightgraph.o  $(ODIR)/mcimage.o  -o $(BDIR)/ppm2graph $(LIBS)



# ===============================================================
# LIBRAIRIE
# ===============================================================


$(ODIR)/mcsort.o:	$(LDIR)/mcsort.c $(IDIR)/mcsort.h
	$(CC) -c $(CCFLAGS) -I$(IDIR) $(LDIR)/mcsort.c -o $(ODIR)/mcsort.o

$(ODIR)/mcunionfind.o:	$(LDIR)/mcunionfind.c $(IDIR)/mcunionfind.h
	$(CC) -c $(CCFLAGS) -I$(IDIR) $(LDIR)/mcunionfind.c -o $(ODIR)/mcunionfind.o		

$(ODIR)/mcweightgraph.o:	$(LDIR)/mcweightgraph.c $(IDIR)/mcweightgraph.h $(IDIR)/mcsort.h
	$(CC) -c $(CCFLAGS) -I$(IDIR) $(LDIR)/mcweightgraph.c -o $(ODIR)/mcweightgraph.o


$(ODIR)/lhierarchie.o:	$(LDIR)/lhierarchie.c $(IDIR)/lhierarchie.h  $(IDIR)/mcweightgraph.h $(IDIR)/mcunionfind.h $(IDIR)/llca.h
	$(CC) -c $(CCFLAGS) -I$(IDIR) $(LDIR)/lhierarchie.c -o $(ODIR)/lhierarchie.o 

$(ODIR)/llca.o:	$(LDIR)/llca.c $(IDIR)/llca.h 
	$(CC) -c $(LDIR)/llca.c -I$(IDIR)  -o $(ODIR)/llca.o

$(ODIR)/lcomptree.o: $(LDIR)/lcomptree.c $(IDIR)/lcomptree.h $(IDIR)/mcunionfind.h $(IDIR)/mcweightgraph.h
	$(CC) -c $(LDIR)/lcomptree.c -I$(IDIR) -o $(ODIR)/lcomptree.o

$(ODIR)/mcimage.o:	$(LDIR)/mcimage.c $(IDIR)/mcimage.h $(IDIR)/mccodimage.h $(IDIR)/mcutil.h
	$(CC) -c $(CCFLAGS) -I$(IDIR) $(LDIR)/mcimage.c -o $(ODIR)/mcimage.o
