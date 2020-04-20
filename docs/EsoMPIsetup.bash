#!/bin/bash
#
# copying all relevant files in the work directory. TAKES THE WORK DIRECTORY AS INPUT
# files are: (in src)
# fitslib.py  workarea.py  thumblib.py runEBM.py tAtmo.py libraryEBM.py
# this one MUST reside in the tree root of the code!

cp src/parallel.py $1
cp src/workarea.py $1
cp src/thumblib.py $1
cp src/fitslib.py $1
cp src/runEBM.py $1
cp src/tAtmo.py $1
cp src/libraryEBM.py $1
cp src/constantsEBM.py $1
cp src/coderoot.py $1

