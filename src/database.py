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
simtype = "Std"
version = "1.1.03"

# Directories and Files
template_dir="/home/python/Esoclimi/Templates"
RisultatiMultipli = "RisultatiMultipli"
Database = "Database"
LogFiles = "LogFiles"
Risultati ="Risultati"
Src = "Src"
workDir="/home/python/Esoclimi/Tests"


fits_out_base = "Database/ESTM1.1.01-"                       # base name of the FITs output
fortran_run_result="year_lat_temp_last1.tlt"        # Data File name
fits_param_file = "esopianeti.par"


'''
    make_work_area (_dir,_p,_ecc,_obl,_dist)
    
    create the working area dircetory and file structure coping from
    template directories
    
    p    = pressure
    ecc  = eccentricity
    obl  =
    dist = distance
    
    '''
def make_work_area (_dir,_p,_ecc,_obl,_dist):
    # multiple output directory preparation
    # here we need to catch exeptions (necessary on parallel runs)
    
    
    localSrc = _dir+"/"+Src
    os.makedirs(localSrc+"/"+Risultati)

    shutil.copytree(template_dir+"/CCM_RH60", localSrc+"/CCM_RH60") # GT why here?
    
    copyall(template_dir+"/ModulesDef", localSrc)
    
    copyall(template_dir+"/Std",localSrc)
    if simtype=="VegPassive":
        copyall(template_dir+"/VegPassive",localSrc)
        copyall(template_dir+"/VegPassive/Modules/",localSrc)
    elif simtype=="VegAlbedoFB":
        copyall(template_dir+"/VegAlbedoFB", localSrc)
        copyall(template_dir+"/VegAlbedoFB/Modules/",localSrc)
    # varying pressure on an EARTH-LIKE planet, using EARTH template
    shutil.copy(template_dir+"/Planets/EARTH.py",localSrc+"/PRESSURES.py")
    shutil.copy(template_dir+"/Planets/fo_earth_DMAP.dat",localSrc+"/fo_earth_DMAP.dat")
        ##
        #   configuring planet parameter file
        #
        #   PRESSURES.py is the template (here, EARTH.py)   GT: maybe better to call it PRESSURES.templ ?
        #   PRESSUREScurr.py is what runEBM.py uses
        ###
    fileout=open(localSrc+'/PRESSUREScurr.py',"w")
    with open(localSrc+'/PRESSURES.py', 'r') as f:
            for tl in f:
                pl = tl
                if 'ff=1.0' in tl:
                    pl="ff=%f\n"%_p
                if 'eccP=' in tl:
                    pl='eccP=%f\n'%_ecc
                if 'obliq=' in tl:
                    pl='obliq=%f\n'%_obl
                if 'smaP=' in tl:
                    pl='smaP=%f\n'%_dist
                fileout.write(pl)
    fileout.close()
    return



def esoclimi(numero):
    for ecc in Eccentricities:
        for obl in Obliquities:
            for dist in Radii:
                for p in Pressures:
                    localWorkDir="%s/%d" % (workDir,numero)
                    os.mkdir(localWorkDir)
                    # initilize log file for simulation

                    logging.info("%d => Begin computation for p=%f ecc=%f obl=%f dits=%s",numero, p,ecc,obl,dist)
                    make_work_area(localWorkDir,p,ecc,obl,dist)
                    os.chdir(localWorkDir)
                    
                    #insert a log here
                    logging.info("%s",os.getcwd())
                    #executing runEBM.py - this compiles and runs the code
                    #GT to include in this code
                    str="python runEBM.py PRESSUREScurr.py %d %s %s > log " % (numero,version,simtype)
                    logging.info("%d => %s",numero, str)
                    #setupEBM(PRESSUREScurr,numero,version,simtype, localWorkDir+"/"+Src)
                    #compileEBM(localWorkDir+"/"+Src)
                    #runEBM(localWorkDir+"/"+Src)
                    #now we have results. Making .fits file (REQUIRES PYFIT)
                    #os.chdir("Risultati")
                    logging.info("%d => %s", numero, "Create FITS file from data")
                    
                    #create_FITS(fortran_result_file,fits_out_base,fits_param_file)
                    
                    ## how many files for a run? if one this is not needed
                    #
                    #    for file in os.listdir("/mydir"):
                    #       if file.endswith(".fits") and file.find("ESTM"):
                    #            shutil.copy(file,"../../Database")
                    #
                    os.chdir(localWorkDir)
                    logging.debug("%d => %s",numero, os.getcwd())
                    logging.debug("%d => %s",numero, os.listdir("."))
                    
                    #archiving Risults
                    results_location="%s/Risultati_Press%5.3f_Ecc%4.2f_Dist%3.1f_Obl%5.3f"%(workDir+"/"+RisultatiMultipli,p,ecc,dist,obl)
                    logging.debug("%d => Archive results to: %s",numero, results_location)
                    archive_results(results_location,Src,localWorkDir+"/"+Src+"/"+Risultati)

                    # close log file and archiving it DO We NEED THAT?
                    #log_dir_name="%s/log_Press%5.3f_Ecc%4.2f_Dist%3.1f_Obl%5.3f"%(LogFiles,p,ecc,dist,obl)
                    #logging.debug("Closing log file and archive to: %s",log_dir_name)
                    #archive_logs(log_dir_name,localWorkDir+"/run.log") # WHICH LOGS?

                    #next simulation
                    numero=numero+1
                    
                    #back to main dir
                    os.chdir(workDir)
                    sys.exit()
                    CleanAllPartialResults(localWorkDir)
    
    return


if __name__ == '__main__':

    numero=1

    os.makedirs(RisultatiMultipli)
    os.makedirs(Database)
    os.makedirs(LogFiles)

    logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s %(message)s',
                    filename=workDir+"/run.log",
                    filemode='w')


    esoclimi(numero)
