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
from runEBM import *

Pressures=[0.01, 0.1, 0.5, 1.0, 3.0, 5.0]
Radii= [0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5]
Obliquities = [0., 15., 23.43929, 30., 45.]
Eccentricities = [ 0.0, 0.01671022, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
simtype = "Std"
version = "1.1.03"

# Directories and Files
template_dir="/home/LAVORO/Programming/Esoclimi/Templates"
workDir = os.getcwd()

RisultatiMultipli = "RisultatiMultipli"
Database = "Database"
LogFiles = "LogFiles"
Risultati ="Risultati"
Src = "Src"
planet = "EARTH"



fits_out_base = workDir+"/Database/ESTM1.1.01-"     # base name of the FITs output
fortran_run_result="year_lat_temp_last1.tlt"        # Data File name
fits_param_file = "/esopianeti.par"


def make_work_area (_dir):
    '''
        make_work_area (_dir,_p,_ecc,_obl,_dist)
        
        create the working area dircetory and file structure coping from
        template directories
        
        p    = pressure
        ecc  = eccentricity
        obl  =
        dist = distance
        
        '''
    # here we need to catch exeptions (necessary on parallel runs)
    localSrc = _dir+"/"+Src
    os.makedirs(localSrc+"/"+Risultati)

    shutil.copytree(template_dir+"/CCM_RH60", localSrc+"/CCM_RH60")
    
    copyall(template_dir+"/ModulesDef", localSrc) # they are all empty files ????
    
    copyall(template_dir+"/Std",localSrc)
    if simtype=="VegPassive":
        copyall(template_dir+"/VegPassive",localSrc)
        copyall(template_dir+"/VegPassive/Modules/",localSrc)
    elif simtype=="VegAlbedoFB":
        copyall(template_dir+"/VegAlbedoFB", localSrc)
        copyall(template_dir+"/VegAlbedoFB/Modules/",localSrc)

    # varying pressure on an EARTH-LIKE planet, using EARTH template
    #shutil.copy(template_dir+"/Planets/"+planet,workDir+"/PRESSURES.py")
    if planet == "EARTH":
        shutil.copy(template_dir+"/Planets/fo_earth_DMAP.dat",localSrc+"/fo_earth_DMAP.dat")
    return



def esoclimi(numero,ecc,obl,dist,p):
     localWorkDir    = "%s/%d/" % (workDir,numero)
     localSrcDir     = "%s%s/" % (localWorkDir,Src)
     localResultDir = "%s%s/" % (localSrcDir,Risultati)
     os.chdir(workDir)
     os.mkdir(localWorkDir)
     # initilize log file for simulation
     
     logging.info("%d => Begin computation for p=%f ecc=%f obl=%f dits=%s",numero, p,ecc,obl,dist)
     make_work_area(localWorkDir)
     os.chdir(localWorkDir)
     
     logging.info("%s",os.getcwd())
     #Complile and Run
     str="python runEBM.py PRESSUREScurr.py %d %s %s > log " % (numero,version,simtype)
     logging.info("%d => %s",numero, str)
     setupEBM(planet,p,ecc,obl,dist,localSrcDir,numero,version,simtype)
     log_file_local=open(localWorkDir+"/out.log",'w')    #open log file name is out.log
     compileEBM(localSrcDir,log_file_local)
     runEBM(localSrcDir,log_file_local)
     log_file_local.close()
     
     #now we have results. Making .fits file (REQUIRES PYFIT)
     logging.info("%d => %s", numero, "Create FITS file from data")
     create_FITS(localResultDir+fortran_run_result,fits_out_base,localResultDir+fits_param_file)
     
     os.chdir(localWorkDir)
     logging.debug("%d => %s",numero, os.getcwd())
     logging.debug("%d => %s",numero, os.listdir("."))
     
     #archiving Risults
     results_location="%s/Risultati_Press%5.3f_Ecc%4.2f_Dist%3.1f_Obl%5.3f"%(workDir+"/"+RisultatiMultipli,p,ecc,dist,obl)
     logging.debug("%d => Archive results to: %s",numero, results_location)
     archive_results(results_location,Src,localResultDir)
     
     # close log file and archiving it DO We NEED THAT?
     log_dir_name="%s/%s/log_Press%5.3f_Ecc%4.2f_Dist%3.1f_Obl%5.3f"%(workDir,LogFiles,p,ecc,dist,obl)
     logging.debug("Closing log file and archive to: %s",log_dir_name)
     archive_logs(log_dir_name,localWorkDir+"/out.log")
     
     #next simulation
     numero=numero+1
     
     #back to main dir
     os.chdir(workDir)
     CleanAllPartialResults(localWorkDir)
     return


if __name__ == '__main__':
    
    # make directories where final results are stored
    os.makedirs(RisultatiMultipli)
    os.makedirs(Database)
    os.makedirs(LogFiles)
    ##
    # open a logger
    logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s %(message)s',
                    filename=workDir+"/run.log",
                    filemode='w')

    numero=1
    p = Pressures[1]
    dist = Radii[1]
    obl = Obliquities[1]
    ecc = Eccentricities[1]
    esoclimi(numero,ecc,obl,dist,p)
        
