#!/bin/sh

linux/bin/pgm2graph $1 0 /tmp/_1.graph
linux/bin/extinctionvalues /tmp/_1.graph $2 /tmp/_2.graph /tmp/BM.txt
linux/bin/uprooting /tmp/_2.graph /tmp/_3.graph
linux/bin/ultramopen /tmp/_3.graph /tmp/_4.graph /tmp/BM.txt
linux/bin/graph2pgm /tmp/_4.graph 0 /tmp/_5.pgm
linux/bin/float2byte /tmp/_5.pgm 4 $3
