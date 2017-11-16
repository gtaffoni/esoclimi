#!/usr/bin/env python
#
# This is an example of python script for runnung one time the
# code, varying one parameter existing in one module .h file
# Here is the vegetation albedo in vegetation.h
#
#   GT 20171002 created
#   GT 20171006 fix function esoclimi
#


'''
    With this driver we can setup, compile and run EBM code.
    It requires a set of valid parameters and a configuration file: ConfigAll.ini
    
'''

# system imports
import sys
import os
import shutil
from posix import system
import logging
import time
from time import  time
from time import  sleep
import numpy as np
import random


##############################################################################
# Global parameters
Src               = "Src"           #TIENI


##############################################################################
# program imports
from fitslib import create_FITS
from thumblib import create_THUMBNAILS
from workarea import *
from runEBM import *
from tAtmo import *
import uuid

def make_work_area (_workDir,_code_root_dir,Risultati,Parameters):
    '''
        make_work_area (_dir)
        
        create the working area directory and file structure coping from
        template directories
        code will run in here
        
        '''
    # here we need to catch exeptions (necessary on parallel runs)
    _dir    = "%s/%d/" % (_workDir,Parameters['number'])
    errorcode=0
    template_dir=_code_root_dir+"/Templates/"
    localSrc = _dir+"/"+Src
    localResults=localSrc+"/"+Risultati
    logging.info("Simulation %d => Begin create Local work directories: %s",Parameters['number'], _dir)
    logging.debug("Simulation %d => Create Local Src: %s",Parameters['number'], localSrc)
    if os.path.isdir(localSrc):
        shutil.rmtree(localSrc)
    try:
        os.makedirs(localSrc+"/"+Risultati)
    except:
        logging.error(sys.exc_info()[0])
        raise
    logging.debug("Simulation %d => Create Parameters.txt file",Parameters['number'])
    with open(localSrc+"/"+Risultati+"/Parameters.txt", "w") as text_file:
        text_file.write(Parameters['data'])
    logging.debug("Simulation %d => Copy CCM_RH60",Parameters['number'])
    try:
        shutil.copytree(template_dir+"/CCM_RH60", localSrc+"/CCM_RH60")
    except:
        logging.error(sys.exc_info()[0])
        raise
    logging.debug("Simulation %d => Copy ModulesDef",Parameters['number'])
    copyall(template_dir+"/ModulesDef", localSrc) # they are all empty files ????
    logging.debug("Simulation %d => Copy Std",Parameters['number'])
    copyall(template_dir+"/Std",localSrc)
    if Parameters['simtype']=="VegPassive":
        logging.debug("Simulation %d => Copy VegPassive",Parameters['number'])
        copyall(template_dir+"/VegPassive",localSrc)
        copyall(template_dir+"/VegPassive/Modules/",localSrc)
    elif Parameters['simtype']=="VegAlbedoFB":
        logging.debug("Simulation %d => Copy VegAlbedoFB",Parameters['number'])
        copyall(template_dir+"/VegAlbedoFB", localSrc)
        copyall(template_dir+"/VegAlbedoFB/Modules/",localSrc)
    # varying pressure on an EARTH-LIKE planet, using EARTH template
    logging.debug("Simulation %d => Copy Planets .h",Parameters['number'])
    shutil.copy(template_dir+"/Planets/"+Parameters['planet']+".h", localSrc+"/planet.h")
    if Parameters['planet'] == "EARTH":
        logging.debug("Simulation %d => Copy fo_earth_DMAP.dat for EARTH",Parameters['number'])
        shutil.copy(template_dir+"/Planets/fo_earth_DMAP.dat",localSrc+"/fo_earth_DMAP.dat")
    logging.info("Simulation %d => End create Local work directories: ErrorCode=%d",Parameters['number'],errorcode)
    return errorcode

def emulation(Parameters):
    '''
        Emulate compile and run (only for debug)
    '''
    _return_code=0
    ev=random.randint(0,5)
    if ev==0:
        exitValue=-100
    elif ev==1:
        exitValue=-1.0
    elif ev==2:
        exitValue=-0.5
    elif ev==3:
        exitValue=-2.
    else:
        exitValue=1
    if exitValue< -101.:
        logging.warning("Emulation %d => CATASTROPHIC EXIT VALUE FOUND, ELABORATION STOPPED", Parameters['number'])
        # TODO get this and stop computation
        return(-200)
    elif np.abs(exitValue + 100) < 0.001 :
        logging.warning("Emulation %d => Integration error for case: %s ",  Parameters['number'],Parameters['data'])
        _return_code = -100
    elif np.abs(exitValue + 2.0) < 0.001 : #pressure exceeded (should not happen!)
        logging.info("Emulation %d => Warning, pressure exceeded for case: %s ",  Parameters['number'],Parameters['data'])
        _return_code = -2
    elif np.abs(exitValue + 1.0) < 0.001 : #Runaway GreenHouse
        logging.info("Emulation %d => Runaway GreenHouse: %s ",  Parameters['number'],Parameters['data'])
        _return_code = -1
    elif np.abs(exitValue + 0.5) < 0.001: #SnowBall
        logging.info("Emulation %d => SnowBall: %s ",  Parameters['number'],Parameters['data'])
        _return_code = -0.5
    return(_return_code)


def exoclime(Parameters,workDir,code_work_dir,Risultati,emulate=False):
    '''
        do preparation, compile and run of code
        save results
        check the status
        '''
    # INITIALIZE SOME VARIABLES
    fits_out_base = workDir+"/Database/ESTM1.1.01-"     # base name of the FITs output
    thumbs_out_base = workDir+"/Thumbnails/ESTM1.1.01-" # base name of the FITs output
    fits_param_file = "/esopianeti.par"
    fortran_run_result="year_lat_temp_last1.tlt"        # Data File name
    fortran_value_result="/valori.txt"
    #
    return_code=0
    #
    localWorkDir    = "%s/%d/" % (workDir,Parameters['number'])
    localSrcDir     = "%s%s/" % (localWorkDir,Src)
    localResultDir = "%s%s/" % (localSrcDir,Risultati)
    # TODO: check if dir exists and decide: to remove to stop.
    os.chdir(workDir)

    logging.info("Simulation %d => Begin computation for p=%f ecc=%f obl=%f dist=%s gtype=%d fo=%4.2f",Parameters['number'], Parameters['p'],Parameters['ecc'],Parameters['obl'],Parameters['dist'],Parameters['gg'],Parameters['fo_const'])
    # TODO: non riesce a fare directory: deve uscire con errore 1
    os.chdir(localWorkDir)
    
    #############################
    # Setup EBM                 #
    #############################
    # TODO: non riesce deve uscire con errore 1
    str="SetupEBM %s %s" % (Parameters['version'],Parameters['simtype'])
    logging.info("Simulation %d => %s",Parameters['number'], str)
    setupEBM(Parameters,localSrcDir)
    EBM_Code_log_file=open(localWorkDir+"/out.log",'w')    #open log file name is out.log
    
    ##############################
    # Emulation dfor debugging   #
    ##############################
    if emulate:
        sleep(2)
        exitValue = emulation(Parameters)
        os.chdir(localWorkDir)
        logging.debug("Emulation %d => ExitValue %d",Parameters['number'], exitValue)
        logging.debug("Emulation %d => %s",Parameters['number'], os.getcwd())
        logging.debug("Emulation %d => %s",Parameters['number'], os.listdir("."))
        return(exitValue)

    #############################
    # compiling EBM             #
    #############################
    # TODO: non riesce: deve uscire con errore 2
    str="Compile EBM %s %s" % (Parameters['version'],Parameters['simtype'])
    logging.info("Simulation %d => %s",Parameters['number'], str)
    try:
        exiterror=compileEBM(localSrcDir,EBM_Code_log_file)
    except:
        print("Simulation %d => compile error %s", sys.exc_info()[0])
        logging.error("Simulation %d => compile error %s", sys.exc_info()[0])
        raise
    if not exiterror == 0:
        print("Simulation %d => compile errors %s", sys.exc_info()[0])
        logging.debug("Simulation %d => compile error %s", sys.exc_info()[0])
        raise ValueError('Compilation Error')

    #############################
    # running EBM               #
    #############################
    str="Run EBM %s %s" % (Parameters['version'],Parameters['simtype'])
    logging.info("Simulation %d => %s",Parameters['number'], str)
    try:
        exiterror=runEBM(localSrcDir,EBM_Code_log_file)
    except:
        logging.error("Simulation %d => RUN error %s", Parameters['number'], sys.exc_info()[0])
        raise
    ########################################################
    #calculating atmospheric thickness (Michele Maris code)#
    ########################################################
    # TODO: Change with Mechele help
    #       Check if parfile exists
    try:
        logging.info("Simulation %d => Starting atmospheric thickness code.", Parameters['number'])
        tAtmo(localResultDir+fits_param_file, code_work_dir, EBM_Code_log_file)
    except:
        e = sys.exc_info()[0]
        logging.warning("Simulation %d %s=> It was impossible to run atmospheric thickness code.",Parameters['number'],e)

    EBM_Code_log_file.close()

    #########################################################
    # Run completed       Making .fits file (REQUIRES PYFIT)#
    #########################################################
    logging.info("Simulation %d => %s", Parameters['number'], "Create FITS file from data")
    try:
        date = create_FITS(localResultDir+fortran_run_result,fits_out_base,localResultDir+fits_param_file)
    except:
        e = sys.exc_info()[0]
        logging.warning("Simulation %d => %s Cannot create FITS file",Parameters['number'],e)
        pass
    try:
        create_THUMBNAILS(localResultDir+fortran_run_result,thumbs_out_base, date, Parameters['number'])
    except:
        e = sys.exc_info()[0]
        logging.warning("Simulation %d => %s Cannot create thumbnails",Parameters['number'],e)
        pass


###########################################################################################
# VERY IMPORTANT: CHECKING NON-CONVERGED SNOWBALL/RUNAWAY GREENHOUSE CASES
# no fits produced in that case.
# Fortran code exit value is produce - saving parameters for which we have SB/RG or other weird cases
#
# ExitValue:
#           -2   -> pressure exceeded
#           -0.5 -> snowball (stops) (nTlin)
#           -1   -> Runaway GreenHouse (nSigma)
#           -100 -> integration error
#            1   -> warm
#            2   -> warm-hot
#            3   -> snowball (converged)
#            4   -> waterbelt
#           -200 -> undefined (should not happen)
#
###########################################################################################

    fortran_value_result_file =localResultDir+fortran_value_result
    logging.info("Simulation %d => Open File fortran_value_result: %s", Parameters['number'], fortran_value_result_file)
    exitValue  = np.loadtxt(fortran_value_result_file)
    logging.info("Simulation %d => Fortran Code Exit Value: %d", Parameters['number'], exitValue[25])
    
    
    #
    #   WARNING: TO BE MODIFIED WHEN CHANGING PARAMETERS SPACE EXPLORATION
    #     we definitely must authomatize this!
    
    
    
    if exitValue[25]< -101.: #### TODO: we need to change this, here we are not closing anything and printing anyting
        logging.warning("Simulation %d => CATASTROPHIC EXIT VALUE FOUND, ELABORATION STOPPED", Parameters['number'])
                # TODO get this and stop computation
        return(-200)
    elif np.abs(exitValue[25] + 100) < 0.001 :
        logging.warning("Simulation %d => Integration error for case: %s ",  Parameters['number'],Parameters['data'])
    elif np.abs(exitValue[25] + 2.0) < 0.001 : #pressure exceeded (should not happen!)
        logging.info("Simulation %d => Warning, pressure exceeded for case: %s ",  Parameters['number'],Parameters['data'])
    elif np.abs(exitValue[25] + 1.0) < 0.001 : #Runaway GreenHouse
        logging.info("Simulation %d => Runaway GreenHouse: %s ",  Parameters['number'],Parameters['data'])
    elif np.abs(exitValue[25] + 0.5) < 0.001: #SnowBall
        logging.info("Simulation %d => SnowBall: %s ",  Parameters['number'],Parameters['data'])
    

    os.chdir(localWorkDir)
    logging.debug("Simulation %d => %s",Parameters['number'], os.getcwd())
    logging.debug("Simulation %d => %s",Parameters['number'], os.listdir("."))

    return(exitValue[25])


def archive_exoplanet_data(Parameters, workDir,RisultatiMultipli,LogFiles,Risultati):
    '''
        Archiving data from execution
        '''
    # TODO: this string should be fixed when more params are needed
    results_string="_Press%5.3f_Ecc%4.2f_Dist%3.1f_Obl%5.3f_CO2_%5.3f_GG%d_fo%4.2f"%(Parameters['p'],Parameters['ecc'],Parameters['dist'],Parameters['obl'],Parameters['CO2_Earth_ratio'],Parameters['gg'],Parameters['fo_const'])
    
    localWorkDir    = "%s/%d/" % (workDir,Parameters['number'])
    localSrcDir     = "%s%s/" % (localWorkDir,Src)
    localResultDir = "%s%s/" % (localSrcDir,Risultati)
    
    results_location="%s/Risultati%s"%(workDir+"/"+RisultatiMultipli,results_string)
    logging.debug("Simulation %d => Archive results to: %s",Parameters['number'], results_location)
    #
    # check if dir exists eventually remove (this may happen in case of restart from crash)
    #
    if os.path.isdir(results_location):
        shutil.rmtree(results_location)
    # archive dir
    archive_results(results_location,Src,Parameters['planet'],localResultDir)
    # close log file and archiving it DO We NEED THAT? YES!
    log_dir_name="%s/%s/log%s"%(workDir,LogFiles,results_string)
    #
    # check if file exists (this may happen in case of restart from crash)
    #
    if os.path.isfile(log_dir_name):
        os.remove(log_dir_name)
    # archive file
    logging.debug("Closing log file and archive to: %s",log_dir_name)
    archive_logs(log_dir_name,localWorkDir+"/out.log")
    
    #back to main dir and clean all results
    os.chdir(workDir)
    CleanAllPartialResults(localWorkDir)
    return

def archive_broken_simulations(Parameters, _workDir, Broken):
    '''
        Archiving data from execution
        '''
    import uuid
    
    # TODO: this string should be fixed when more params are needed
    results_string="_Press%5.3f_Ecc%4.2f_Dist%3.1f_Obl%5.3f_CO2_%5.3f_GG%d_fo%4.2f_%d"%(Parameters['p'],Parameters['ecc'],Parameters['dist'],Parameters['obl'],Parameters['CO2_Earth_ratio'],Parameters['gg'],Parameters['fo_const'],Parameters['number'])
    
    localWorkDir    = "%s/%d/" % (_workDir,Parameters['number'])
    
    results_location="%s/Broken%s_%s"%(_workDir+"/"+Broken,results_string,str(uuid.uuid1()))
    logging.debug("Simulation %d => Archive broken results to: %s",Parameters['number'], results_location)
    #
    # check if dir exists eventually remove (this is unlikely to  happen UUID)
    #
    if os.path.isdir(results_location):
        shutil.rmtree(results_location)
    # archive dir
    try:
        shutil.move(localWorkDir,results_location)
    except:
        raise
    os.chdir(_workDir)
    return

