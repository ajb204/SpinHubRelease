#!/bin/bash

rm ../../bin/specfrq
rm ../../bin/CountMaxThread
rm ../../bin/RelaxFix.out
rm ../../bin/RelaxFix2.out
rm ../../bin/re-shuffle.out
rm ../../bin/re-shuffle2.out
rm ../../bin/seqconv
rm ../../bin/AddFIDvarian.out

#/opt/homebrew/Cellar/llvm/22.1.2/bin/clang++ src/CountMaxThead.cpp -fopenmp -o ../../bin/CountMaxThread
g++ src/specfrq_v2.c -o ../../bin/specfrq
g++ src/CountMaxThread.cpp  -o ../../bin/CountMaxThread
g++ src/RelaxFix.cpp -o ../../bin/RelaxFix.out
g++ src/RelaxFix2.cpp -o ../../bin/RelaxFix2.out
g++ src/re-shuffle.cpp -o ../../bin/re-shuffle.out
g++ src/re-shuffle2.cpp -o ../../bin/re-shuffle2.out
g++ src/seqconv.cpp -o ../../bin/seqconv.out
g++ src/AddFIDvarian.cpp -o ../../bin/AddFIDvarian.out
