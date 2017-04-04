#!/usr/bin/env python
#
# This is an example of python script for runnung many times the
# code, varying one parameter existing in one module .h file
# Here is the vegetation albedo in vegetation.h
#
#   GM created
#   GT 29/3/17 some small improvements
#
#

'''
    Input paramters of the computational kernel:
    * simtype,
    * version,
    * simulation number
    
    
    Values of Pressure, Raddi Obliquities and Eccentricities to cycle
    
    WARNING, 0.001 e 0.005 of pressure are not working
    Pressures are epressed in times the Earth value

'''

import sys
import os
import shutil
from posix import system
import logging

from fitslib import create_FITS
from workarea import *

Pressures=[0.01, 0.1, 0.5, 1.0, 3.0, 5.0]
Radii= [0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5]
Obliquities = [0., 15., 23.43929, 30., 45.]
Eccentricities = [ 0.0, 0.01671022, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
simtype="Std"
version="1.1.03"
number=1

# Directories and Files
template_dir="/work/Programming/Esoclimi/Devel"
RisultatiMultipli = "RisultatiMultipli"
Database = "Database"
LogFiles = "LogFiles"
Risultati ="Risultati"
Src = "Src"


fits_out_base = "ESTM1.1.01-"                       # base name of the FITs output
fortran_run_result="year_lat_temp_last1.tlt"        # Data File name
fits_param_file = "NONMELORICORDO"

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s %(message)s',
                    filename='run.log',
                    filemode='w')







for ecc in Eccentricities:
    for obl in Obliquities:
        for dist in Radii:            
            for p in Pressures:
                print p,ecc,obl,dist # BETTER USE A LOGGER
            
                make_work_area()

                os.chdir(Src)
                
                #insert a log here
                logging.info("%s",os.getcwd())
                #executing runEBM.py - this compiles and runs the code
                #GT to include in this code
                str="python runEBM.py PRESSUREScurr.py %d %s %s > log "%(number,version,simtype)
                logging.info("%s",str)
                #system(str)

                #now we have results. Making .fits file (REQUIRES PYFIT)
                #os.chdir("Risultati")
                logging.info("%s","python ../../fits_map_temperature.py")
                
                #create_FITS(fortran_result_file,fits_out_base,fits_param_file)
                
                
                for file in os.listdir("/mydir"):   #maybe not necessary if include the fits conveter
                    if file.endswith(".fits") and file.find("ESTM"):
                        shutil.copy(file,"../../Database")

                os.chdir("..")
                
                #archiving log file
                logging.info("%s",os.getcwd())
                logging.info("%s",os.listdir())
                archive_logs(p,ecc,dist,obl)
                
                #archiving Risults
                archive_results(p,ecc,dist,obl)
                
                #next simulation
                number=number+1

                #back to main dir
                os.chdir("..")
                print "\n\n\n"


