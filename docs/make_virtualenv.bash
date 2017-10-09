#!/bin/bash
virtualenv local-python
source local-python/bin/activate
pip install --upgrade pip
pip install numpy
pip install scipy
pip install matplotlib
pip install mpy4py
pip install mpi4py
pip install pyfits
pip install ipython

