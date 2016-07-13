#!/bin/bash
lyx -e literate $1.lyx
notangle -t8 -RMakefile $1.nw > Makefile
make all
bash filenamehack.bash
