SEGSURF = .
# Version LINUX
  XLIB = -L/usr/X11R6/lib -lX11 -lXext
  XINCL = /usr/include/X11R6
  PNGLIB = /usr/local/lib
  CC = cc
#DEBUG
  #CCFLAGS = -pg -rdynamic -DPC -DUNIXIO -lm
 # CCFLAGS = -g 
#  CCFLAGS = -pg -fno-inline
##RELEASE
  #CCFLAGS = -o3 -rdynamic -lm
 CCFLAGS = -O2

  ODIR = ./linux/obj
  BDIR = ./linux/bin
  LIBS = -lm 
  
# $(BDIR)/ultramopen\
# $(BDIR)/segmentImage\
# $(BDIR)/ppm2graph
SAL=\
$(BDIR)/ultramopen\
$(BDIR)/hgbSegmentationInterval

OBJ = $(ODIR)

IDIR = ./include

IDIRPNG = /usr/local/include

CDIR = ./src/com

LDIR = ./src/lib

all:  $(SAL) 

clean:	
	rm -f ./makefile
	rm -f ./linux/obj/*
	rm -f ./linux/bin/*
	rm -f ./linux/bin2/*	
	rm -f ./hpux/bin/*
	rm -f ./hpux/obj/*		
	rm -f ./*.o




# ===============================================================
# EXECUTABLES
# ===============================================================


$(BDIR)/ultramopen: $(CDIR)/ultramopen.c $(IDIR)/mcweightgraph.h $(IDIR)/lhierarchie.h $(IDIR)/llca.h  $(IDIR)/lcomptree.h  $(IDIR)/mcunionfind.h $(IDIR)/mcsort.h $(ODIR)/mcweightgraph.o $(ODIR)/lhierarchie.o $(ODIR)/llca.o $(ODIR)/lcomptree.o  $(ODIR)/mcunionfind.o $(ODIR)/mcsort.o $(ODIR)/MST.o $(ODIR)/list.o $(ODIR)/graphSegmentation.o 
	$(CC) $(CCFLAGS) -I$(IDIR) $(CDIR)/ultramopen.c  $(ODIR)/mcweightgraph.o   $(ODIR)/mcunionfind.o  $(ODIR)/lhierarchie.o $(ODIR)/llca.o $(ODIR)/lcomptree.o $(ODIR)/mcsort.o $(ODIR)/MST.o $(ODIR)/list.o $(ODIR)/graphSegmentation.o -o $(BDIR)/ultramopen $(LIBS)







$(BDIR)/hgbSegmentationInterval: $(CDIR)/HgbSegInterval.c $(IDIR)/mcweightgraph.h $(IDIR)/lhierarchie.h $(IDIR)/llca.h  $(IDIR)/lcomptree.h  $(IDIR)/mcunionfind.h $(IDIR)/mcsort.h $(ODIR)/mcweightgraph.o $(ODIR)/lhierarchie.o $(ODIR)/llca.o $(ODIR)/lcomptree.o  $(ODIR)/mcunionfind.o $(ODIR)/mcsort.o  $(ODIR)/MST.o $(ODIR)/list.o $(ODIR)/graphSegmentation.o 
	$(CC) $(CCFLAGS) -I$(IDIR) $(CDIR)/HgbSegInterval.c  $(ODIR)/mcweightgraph.o   $(ODIR)/mcunionfind.o  $(ODIR)/lhierarchie.o $(ODIR)/llca.o $(ODIR)/lcomptree.o $(ODIR)/mcsort.o $(ODIR)/MST.o $(ODIR)/list.o $(ODIR)/graphSegmentation.o -o $(BDIR)/HgbSegInterval $(LIBS)



# ===============================================================
# LIBRAIRIE
# ===============================================================

$(ODIR)/MST.o:	$(LDIR)/MST.c $(IDIR)/MST.h
	$(CC) -c $(CCFLAGS) -I$(IDIR) $(LDIR)/MST.c  -o $(ODIR)/MST.o

$(ODIR)/list.o:	$(LDIR)/list.c $(IDIR)/list.h
	$(CC) -c $(CCFLAGS) -I$(IDIR) $(LDIR)/list.c  -o $(ODIR)/list.o	

$(ODIR)/graphSegmentation.o:	$(LDIR)/graphSegmentation.c $(IDIR)/graphSegmentation.h
	$(CC) -c $(CCFLAGS) -I$(IDIR) $(LDIR)/graphSegmentation.c  -o $(ODIR)/graphSegmentation.o

$(ODIR)/dealpng.o:	$(LDIR)/dealpng.c $(IDIR)/dealpng.h
	$(CC) -c $(CCFLAGS) -I$(IDIR) -I$(IDIRPNG) -I$(XINCL) -L$(PNGLIB) $(LDIR)/dealpng.c -lpng  -o $(ODIR)/dealpng.o

$(ODIR)/mcsort.o:	$(LDIR)/mcsort.c $(IDIR)/mcsort.h
	$(CC) -c $(CCFLAGS) -I$(IDIR) $(LDIR)/mcsort.c  -o $(ODIR)/mcsort.o

$(ODIR)/mcunionfind.o:	$(LDIR)/mcunionfind.c $(IDIR)/mcunionfind.h
	$(CC) -c $(CCFLAGS) -I$(IDIR) $(LDIR)/mcunionfind.c  -o $(ODIR)/mcunionfind.o		

$(ODIR)/mcweightgraph.o:	$(LDIR)/mcweightgraph.c $(IDIR)/mcweightgraph.h $(IDIR)/mcsort.h
	$(CC) -c $(CCFLAGS) -I$(IDIR) $(LDIR)/mcweightgraph.c  -o $(ODIR)/mcweightgraph.o

$(ODIR)/lhierarchie.o:	$(LDIR)/lhierarchie.c $(IDIR)/lhierarchie.h  $(IDIR)/mcweightgraph.h $(IDIR)/mcunionfind.h $(IDIR)/llca.h
	$(CC) -c $(CCFLAGS) -I$(IDIR) $(LDIR)/lhierarchie.c  -o $(ODIR)/lhierarchie.o 

$(ODIR)/llca.o:	$(LDIR)/llca.c $(IDIR)/llca.h 
	$(CC) -c $(LDIR)/llca.c -I$(IDIR)   -o $(ODIR)/llca.o

$(ODIR)/lcomptree.o: $(LDIR)/lcomptree.c $(IDIR)/lcomptree.h $(IDIR)/mcunionfind.h $(IDIR)/mcweightgraph.h
	$(CC) -c $(LDIR)/lcomptree.c -I$(IDIR)  -o $(ODIR)/lcomptree.o

$(ODIR)/mcimage.o:	$(LDIR)/mcimage.c $(IDIR)/mcimage.h $(IDIR)/mccodimage.h $(IDIR)/mcutil.h
	$(CC) -c $(CCFLAGS) -I$(IDIR) $(LDIR)/mcimage.c  -o $(ODIR)/mcimage.o
