#!/bin/bash
LDFLAGS="-L/opt/homebrew/opt/openblas/lib"
CPPFLAGS="-I/opt/homebrew/opt/openblas/include"
gcc utility.c main.c -o main $LDFLAGS $CPPFLAGS -llapack -lopenblas
