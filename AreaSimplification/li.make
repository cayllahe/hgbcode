SEGSURF = .
# Version LINUX
  XLIB = -L/usr/X11R6/lib -lX11 -lXext
  XINCL = /usr/include/X11R6
  CC = gcc
  CCFLAGS = -O3 -Wall -ggdb -g -DPC -DUNIXIO -lm
  ODIR = ./linux/obj
  BDIR = ./linux/bin
  LIBS = -lm 
  
