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
  
