#!/usr/bin/env python
#
# This is an example of python script for runnung many times the
# code, varying one parameter existing in one module .h file
# Here is the vegetation albedo in vegetation.h
#
#   GM created
#   GT 29/3/17 some small improvements
#   GM 11/4 fix small bugs, log non-convergent sim
#   GM 19/7 integrated Michele Maris atmosphere thickness computation
#   GM 7/9 exit values management extended; setup problems fixed
#

'''
    Input paramters of the computational kernel:
    * simtype,
    * version,
    * simulation number
    
    
    Values of Pressure, Radii Obliquities and Eccentricities to cycle
    
    WARNING, 0.001 e 0.005 of pressure are not working
    Pressures are epressed in times the Earth value


'''

import sys
import os
import shutil
from posix import system
import logging
from mpi4py import MPI

# coderoot, templates, src directories
##############################################################################
#
import coderoot as cd
code_root_dir=cd.coderoot()
template_dir=code_root_dir+"Templates/"
src_dir=code_root_dir+"src/"
###############################################################################
# qui potremmo mettere la copia di tutti i files ed evitare la shell bash

from fitslib import create_FITS
from thumblib import create_THUMBNAILS
from workarea import *
from runEBM import *
from tAtmo import *
from time import  time
import numpy as np

def enum(*sequential, **named):
    '''
        simple way to emulate enumerate in python taken from the web
    '''
    enums = dict(zip(sequential, range(len(sequential))), **named)
    return type('Enum', (), enums)




Parameter_set = {'simtype': "Std", 'version': "1.1.03", 'planet': "EARTH" }

# Directories and Files
##############################################################################
#  WARNING: THIS HAS TO BE ADAPTED FOR EACH SETUP
#template_dir="/home/LAVORO/Programming/Esoclimi/Devel/Templates/"
#template_dir="/home/murante/data/EsoClimi/ProvaMpiPiuMichele/Templates/"
#
#
import coderoot as cd
code_root_dir=cd.coderoot()
template_dir=code_root_dir+"Templates/"
src_dir=code_root_dir+"src/"
# sarebbe carino dare al main l'argomento della root dell'albero e farlo fare a lui
###############################################################################
workDir = os.getcwd()
 
RisultatiMultipli = "RisultatiMultipli"
Database = "Database"
LogFiles = "LogFiles"
Risultati ="Risultati"
Thumbnails = "Thumbnails"
Src = "Src"

fits_out_base = workDir+"/Database/ESTM1.1.01-"     # base name of the FITs output
thumbs_out_base = workDir+"/Thumbnails/ESTM1.1.01-"     # base name of the FITs output
fortran_run_result="year_lat_temp_last1.tlt"        # Data File name
fits_param_file = "/esopianeti.par"
fortran_value_result="/valori.txt"



def make_work_area (_dir):
    '''
        make_work_area (_dir,_p,_ecc,_obl,_dist)
        
        create the working area directory and file structure coping from
        template directories  
        code will run in here
        
        p    = pressure
        ecc  = eccentricity
        obl  = obliquity
        dist = semi-major axis
        
        '''
    # here we need to catch exeptions (necessary on parallel runs)
    localSrc = _dir+"/"+Src
    os.makedirs(localSrc+"/"+Risultati)

    shutil.copytree(template_dir+"/CCM_RH60", localSrc+"/CCM_RH60")
    
    copyall(template_dir+"/ModulesDef", localSrc) # they are all empty files ????
    
    copyall(template_dir+"/Std",localSrc)
    if Parameter_set['simtype']=="VegPassive":
        copyall(template_dir+"/VegPassive",localSrc)
        copyall(template_dir+"/VegPassive/Modules/",localSrc) 
    elif Parameter_set['simtype']=="VegAlbedoFB":
        copyall(template_dir+"/VegAlbedoFB", localSrc)
        copyall(template_dir+"/VegAlbedoFB/Modules/",localSrc)

    # varying pressure on an EARTH-LIKE planet, using EARTH template
    shutil.copy(template_dir+"/Planets/"+Parameter_set['planet']+".h", localSrc+"/planet.h")
    if Parameter_set['planet'] == "EARTH":
        shutil.copy(template_dir+"/Planets/fo_earth_DMAP.dat",localSrc+"/fo_earth_DMAP.dat")


    return

def esoclimi_emulate(Parameter_set,nSigmaCrit,nTlim,SigmaCritParams,TlimParams,nPressExceeded,nIntegrationError, PressExceededParams,IntegrationErrorParams):
    '''
       Test fctn for debugging the parallel setup
    '''
    import random
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

    print 'ESOCLIMI_EMULATE, random exit: ',exitValue
    print  Parameter_set['p'], Parameter_set['ecc'], Parameter_set['obl'], Parameter_set['dist'] 


    if exitValue< -101.:
        print 'WARNING, CATASTROPHIC EXIT VALUE FOUND, ELABORATION STOPPED'
        print 'Parameters: ', Parameter_set['ecc'], Parameter_set['obl'], Parameter_set['dist'], Parameter_set['o']
        exit(-200)
    elif np.abs(exitValue + 100) < 0.001 :
        nIntegrationError += 1
        IntegrationErrorParams[0] = np.append(IntegrationErrorParams[0],Parameter_set['ecc'])
        IntegrationErrorParams[1] = np.append(IntegrationErrorParams[1],Parameter_set['obl'])
        IntegrationErrorParams[2] = np.append(IntegrationErrorParams[2],Parameter_set['dist'])
        IntegrationErrorParams[3] = np.append(IntegrationErrorParams[3],Parameter_set['p'])
    elif np.abs(exitValue + 2.0) < 0.001 : #pressure exceeded (should not happen!)
        logging.info("Warning, pressure exceeded ")
        nPressExceeded += 1
        PressExceededParams[0] = np.append(PressExceededParams[0],Parameter_set['ecc'])
        PressExceededParams[1] = np.append(PressExceededParams[1],Parameter_set['obl'])
        PressExceededParams[2] = np.append(PressExceededParams[2],Parameter_set['dist'])
        PressExceededParams[3] = np.append(PressExceededParams[3],Parameter_set['p'])
    elif np.abs(exitValue + 1.0) < 0.001 : #Runaway GreenHouse
        nSigmaCrit += 1
        SigmaCritParams[0] = np.append(SigmaCritParams[0],Parameter_set['ecc'])
        SigmaCritParams[1] = np.append(SigmaCritParams[1],Parameter_set['obl'])
        SigmaCritParams[2] = np.append(SigmaCritParams[2],Parameter_set['dist'])
        SigmaCritParams[3] = np.append(SigmaCritParams[3],Parameter_set['p'])
    elif np.abs(exitValue + 0.5) < 0.001: #SnowBall
        nTlim += 1
        TlimParams[0] = np.append(TlimParams[0],Parameter_set['ecc'])
        TlimParams[1] = np.append(TlimParams[1],Parameter_set['obl'])
        TlimParams[2] = np.append(TlimParams[2],Parameter_set['dist'])
        TlimParams[3] = np.append(TlimParams[3],Parameter_set['p'])

    return(nSigmaCrit, nTlim, SigmaCritParams, TlimParams, nPressExceeded, nIntegrationError, PressExceededParams, IntegrationErrorParams, exitValue)
    
 

def esoclimi(Parameter_set,nSigmaCrit,nTlim,SigmaCritParams,TlimParams,nPressExceeded,nIntegrationError,PresExceededParams,IntegrationErrorParams):
     import numpy as np
     localWorkDir    = "%s/%d/" % (workDir,Parameter_set['number'])
     localSrcDir     = "%s%s/" % (localWorkDir,Src)
     localResultDir = "%s%s/" % (localSrcDir,Risultati)
     os.chdir(workDir)
     os.mkdir(localWorkDir)

     #
     #   WARNING: TO BE MODIFIED WHEN CHANGING PARAMETERS SPACE EXPLORATION
     #
     results_string="_Press%5.3f_Ecc%4.2f_Dist%3.1f_Obl%5.3f_CO2_%5.3f_GG%d"%(Parameter_set['p'],Parameter_set['ecc'],Parameter_set['dist'],Parameter_set['obl'],Parameter_set['p_CO2_P'],Parameter_set['gg'])
     # initilize log file for simulation
     
     logging.info("%d => Begin computation for p=%f ecc=%f obl=%f dits=%s",Parameter_set['number'], Parameter_set['p'],Parameter_set['ecc'],Parameter_set['obl'],Parameter_set['dist'])
     make_work_area(localWorkDir)
     os.chdir(localWorkDir)
     
     logging.info("%s",os.getcwd())

     #Complile and Run
     str="python runEBM.py PRESSUREScurr.py %d %s %s > log " % (Parameter_set['number'],Parameter_set['version'],Parameter_set['simtype'])
     logging.info("%d => %s",Parameter_set['number'], str)
     setupEBM(Parameter_set,localSrcDir)
     log_file_local=open(localWorkDir+"/out.log",'w')    #open log file name is out.log

     ###########################
     # compiling and running   #
     ###########################
     compileEBM(localSrcDir,log_file_local)
     runEBM(localSrcDir,log_file_local)
    
     ########################################################
     #calculating atmospheric thickness (Michele Maris code)#
     ########################################################
     try:
        logging.info("Starting atmospheric thickness code.")
        tAtmo(localResultDir+fits_param_file, workDir, log_file_local)
     except:
        logging.warning("It was impossible to run atmospheric thickness code.")
        print 'Run ',Parameter_set['number'], ' atmo code not running'
        print localResultDir+fits_param_file, workDir

     log_file_local.close()

 
     #########################################################
     #now we have results. Making .fits file (REQUIRES PYFIT)#
     #########################################################
     logging.info("%d => %s", Parameter_set['number'], "Create FITS file from data")
     try:
         date = create_FITS(localResultDir+fortran_run_result,fits_out_base,localResultDir+fits_param_file)
         create_THUMBNAILS(localResultDir+fortran_run_result,thumbs_out_base, date, Parameter_set['number'])
     except:
         logging.warning("Simulation did not converge")

     #VERY IMPORTANT: CHECKING NON-CONVERGED SNOWBALL/RUNAWAY GREENHOUSE CASES 
     # (no fits produced in that case!)

     fortran_value_result_file =localResultDir+fortran_value_result
     logging.info("Open File fortran_value_result: %s", fortran_value_result_file)
     exitValue  = np.loadtxt(fortran_value_result_file)
     logging.info('ExitValue: %d', exitValue[25])

     #
     # ExitValue: -2 -> pressure exceeded, -0.5 -> snowball (stops), -1 -> Runaway GreenHouse -100 integration error
     #   1 -> warm, 2-> warm-hot, 3 -> snowball (converged), 4 -> waterbelt, -200 -> undefined (should not happen)
     # saving parameters for which we have SB/RG or other weird cases


     #
     #   WARNING: TO BE MODIFIED WHEN CHANGING PARAMETERS SPACE EXPLORATION
     #     we definitely must authomatize this!

     if exitValue[25]< -101.:
         print 'WARNING, CATASTROPHIC EXIT VALUE FOUND, ELABORATION STOPPED'
         print 'Parameters: ', Parameter_set['ecc'], Parameter_set['obl'], Parameter_set['dist'], Parameter_set['o']
         exit(-200)
     elif np.abs(exitValue[25] + 100) < 0.001 :
         nIntegrationError += 1
         IntegrationErrorParams[0] = np.append(IntegrationErrorParams[0],Parameter_set['ecc'])
         IntegrationErrorParams[1] = np.append(IntegrationErrorParams[1],Parameter_set['obl'])
         IntegrationErrorParams[2] = np.append(IntegrationErrorParams[2],Parameter_set['dist'])
         IntegrationErrorParams[3] = np.append(IntegrationErrorParams[3],Parameter_set['p'])
     elif np.abs(exitValue[25] + 2.0) < 0.001 : #pressure exceeded (should not happen!)
         logging.info("Warning, pressure exceeded for case: %s %s %s %s\n", Parameter_set['ecc'], Parameter_set['obl'], Parameter_set['dist'], Parameter_set['o']) 
         nPressExceeded += 1
         PressExceededParams[0] = np.append(PressExceededParams[0],Parameter_set['ecc'])
         PressExceededParams[1] = np.append(PressExceededParams[1],Parameter_set['obl'])
         PressExceededParams[2] = np.append(PressExceededParams[2],Parameter_set['dist'])
         PressExceededParams[3] = np.append(PressExceededParams[3],Parameter_set['p'])
     elif np.abs(exitValue[25] + 1.0) < 0.001 : #Runaway GreenHouse
         nSigmaCrit += 1
         SigmaCritParams[0] = np.append(SigmaCritParams[0],Parameter_set['ecc'])
         SigmaCritParams[1] = np.append(SigmaCritParams[1],Parameter_set['obl'])
         SigmaCritParams[2] = np.append(SigmaCritParams[2],Parameter_set['dist'])
         SigmaCritParams[3] = np.append(SigmaCritParams[3],Parameter_set['p'])
     elif np.abs(exitValue[25] + 0.5) < 0.001: #SnowBall
         nTlim += 1
         TlimParams[0] = np.append(TlimParams[0],Parameter_set['ecc'])
         TlimParams[1] = np.append(TlimParams[1],Parameter_set['obl'])
         TlimParams[2] = np.append(TlimParams[2],Parameter_set['dist'])
         TlimParams[3] = np.append(TlimParams[3],Parameter_set['p'])

     os.chdir(localWorkDir)
     logging.debug("%d => %s",Parameter_set['number'], os.getcwd())
     logging.debug("%d => %s",Parameter_set['number'], os.listdir("."))
     
     #archiving Risults
     results_location="%s/Risultati_Press%s"%(workDir+"/"+RisultatiMultipli,results_string)
     logging.debug("%d => Archive results to: %s",Parameter_set['number'], results_location)
     archive_results(results_location,Src,Parameter_set['planet'],localResultDir)
     
     # close log file and archiving it DO We NEED THAT? YES!
     log_dir_name="%s/%s/log%s"%(workDir,LogFiles,results_string)
     logging.debug("Closing log file and archive to: %s",log_dir_name)
     archive_logs(log_dir_name,localWorkDir+"/out.log")

    #back to main dir
     os.chdir(workDir)
     CleanAllPartialResults(localWorkDir)
     return(nSigmaCrit, nTlim, SigmaCritParams, TlimParams, nPressExceeded, nIntegrationError, PressExceededParams, IntegrationErrorParams, exitValue[25])


def make_input_parameters(_data,parameters):
    '''
        Convert input from rank 0 into set of parametes
        WARNING: TO BE MODIFIED WHEN CHANGING PARAMETERS SPACE EXPLORATION
    '''
     #
     #   WARNING: TO BE MODIFIED WHEN CHANGING PARAMETERS SPACE EXPLORATION
     #
    input_params=np.fromstring(_data[1], dtype=float, sep=' ')
    parameters['gg']         = 0      #geography
    parameters['fo_const']   = 0.4    #ocean fraction (only for gg=0)
    parameters['p_CO2_P']    = 3800   #CO2 partial pressure IN PPVM
    parameters['TOAalbfile'] = 'CCM_RH60/ALB_g1_rh60_co2x10.txt'
    parameters['OLRfile']    = 'CCM_RH60/OLR_g1_rh60_co2x10.txt'
    parameters['dist'] = input_params[3]    # semi-major axis of planet orbit
    parameters['obl'] = input_params[2]     # planet axis inclination
    parameters['ecc'] = input_params[1]     # eccentricity of planet orbit
    parameters['p'] = input_params[0]       # pressure
    parameters['number'] = _data[0]
    return(parameters)

def write_restart_file(finput,fcomputed,frestart):
    '''
        write a restart file every X minutes
        
        file1 = input data
        file2 = executed data
        file3 = restart data
        '''
    import difflib
    shutil.copyfile(frestart,frestart+".bak")
    file1=open(finput)
    file2=open(fcomputed)
    file3=open(frestart,"w")
    diff = difflib.ndiff(file1.readlines(), file2.readlines())
    delta = ''.join(x[2:] for x in diff if x.startswith('- '))
    file3.write(delta)
    file3.close()
    return

def write_non_converged_models(nSigmaCrit,nTlim,simulation_index,SigmaCritParams,TlimParams,nPressExceeded,nIntegrationError,PressExceededParams,IntegrationErrorParams):
    '''
    
    Write on files the number and type of non converged models and the corresponding paramters
    
    '''
    print 'nSigmaCrit, nTlim, nPressExceeded, nIntegrationError: ', nSigmaCrit,  nTlim, nPressExceeded, nIntegrationError
    print '\n\n\n'
    print 'nSigmaCrit (Runaway Greenhouse), nTlim (Snowball), nPressExceeded (out of range), nIntegrationError (stepsize too small): ', nSigmaCrit,  nTlim, nPressExceeded, nIntegrationError
    print 'Fractions: ', 1.0*nSigmaCrit/simulation_index, 1.0*nTlim/simulation_index, 1.0*nPressExceeded/simulation_index, 1.0*nIntegrationError/simulation_index
    print '\n Overall number and fraction of non-converged inhabitable cases: ', nSigmaCrit+nTlim+nPressExceeded+nIntegrationError, 1.0*(nSigmaCrit+nTlim+nPressExceeded+nIntegrationError) / simulation_index
    print '\n Total number of runs: ',simulation_index
        
        #recording these on a file!
    f=open('NonConverged.dat','w')
    f.write('nSigmaCrit (Runaway Greenhouse), nTlim (Snowball), nPressExceeded (out of range), nIntegrationError (stepsize too small): %d %d %d %d ' % (nSigmaCrit,  nTlim, nPressExceeded, nIntegrationError) )
    f.write('Fractions: %e %e %e %e\n' % (1.0*nSigmaCrit/simulation_index, 1.0*nTlim/simulation_index, 1.0*nPressExceeded/simulation_index, 1.0*nIntegrationError/simulation_index) )
    f.write('Overall number and fraction of non-converged inhabitable cases: %d %e \n' %(nSigmaCrit+nTlim+nPressExceeded+nIntegrationError, 1.0*(nSigmaCrit+nTlim+nPressExceeded+nIntegrationError) / simulation_index))
    f.write('Total number of runs: %d\n' % simulation_index)
    f.close()
    
#
# WARNING for a miracle this should NOT depend on the number of parameters used but CHECK
#

    with open('SnowBall-Params.dat','w') as f:
        for l in np.matrix(TlimParams).T:
            np.savetxt(f,l,'%e ')
    
    with open('RunawayGreenhouse-Params.dat','w') as f:
        for l in np.matrix(SigmaCritParams).T:
            np.savetxt(f,l,'%e ')

    with open('PressExceeded-Params.dat','w') as f:
        for l in np.matrix(PressExceededParams).T:
            np.savetxt(f,l,'%e ')

    with open('IntegrationError-Params.dat','w') as f:
        for l in np.matrix(IntegrationErrorParams).T:
            np.savetxt(f,l,'%e ')

    return


def collect_non_converged_models_data(nSigmaCrit, nTlim, SigmaCritParams, TlimParams, nPressExceeded, nIntegrationError, PressExceededParams, IntegrationErrorParams, data):
    '''
        Collect and sum up all the non converged models data
        WARNING: THIS HAS TO BE CHANGED WHEN THE NUMBER OF PARAMETERS
        IS VARIED
    '''
#
#  WARNING: THIS HAS TO BE CHANGED WHEN THE NUMBER OF PARAMETERS
#           IS VARIED
#   WE *MUST* STUDY AUTOMATIC WAYS TO DO THESE KIND OF THINGS
#

    print '---- collect non converged models data'
    print '---- data: '
    print  data[0]
    print data[1]
    print data[2][0]
    print data[2][1]
    print data[2][2]
    print data[2][3]
    print data[3][0]
    print data[3][1]
    print data[3][2]
    print data[3][3]
    print data[4]
    print data[5]
    print data[6][0]
    print data[6][1]
    print data[6][2]
    print data[6][3]
    print data[7][0]
    print data[7][1]
    print data[7][2]
    print data[7][3]
    print '------'
    try:
     nSigmaCrit = nSigmaCrit + data[0]
     nTlim = nTlim + data[1]
     SigmaCritParams[0] = np.append(SigmaCritParams[0],data[2][0])
     SigmaCritParams[1] = np.append(SigmaCritParams[1],data[2][1])
     SigmaCritParams[2] = np.append(SigmaCritParams[2],data[2][2])
     SigmaCritParams[3] = np.append(SigmaCritParams[3],data[2][3])
     TlimParams[0] = np.append(TlimParams[0],data[3][0])
     TlimParams[1] = np.append(TlimParams[1],data[3][1])
     TlimParams[2] = np.append(TlimParams[2],data[3][2])
     TlimParams[3] = np.append(TlimParams[3],data[3][3])
     nPressExceeded = nPressExceeded + data[4]
     nIntegrationError  = nIntegrationError + data[5]
     PressExceededParams[0] = np.append(PressExceededParams[0],data[6][0])
     PressExceededParams[1] = np.append(PressExceededParams[1],data[6][1])
     PressExceededParams[2] = np.append(PressExceededParams[2],data[6][2])
     PressExceededParams[3] = np.append(PressExceededParams[3],data[6][3])
     IntegrationErrorParams[0] = np.append(IntegrationErrorParams[0],data[7][0])
     IntegrationErrorParams[1] = np.append(IntegrationErrorParams[1],data[7][1])
     IntegrationErrorParams[2] = np.append(IntegrationErrorParams[2],data[7][2])
     IntegrationErrorParams[3] = np.append(IntegrationErrorParams[3],data[7][3])
    except: #GM emergency hotfix! something went horribly wrong, but the show must go on...
     nSigmaCrit = nSigmaCrit 
     nTlim = nTlim 
     SigmaCritParams[0] = np.append(SigmaCritParams[0],-666.0)
     SigmaCritParams[1] = np.append(SigmaCritParams[1],-666.0)
     SigmaCritParams[2] = np.append(SigmaCritParams[2],-666.0)
     SigmaCritParams[3] = np.append(SigmaCritParams[3],-666.0)
     TlimParams[0] = np.append(TlimParams[0],-666.0)
     TlimParams[1] = np.append(TlimParams[1],-666.0)
     TlimParams[2] = np.append(TlimParams[2],-666.0)
     TlimParams[3] = np.append(TlimParams[3],-666.0)
     nPressExceeded = nPressExceeded
     nIntegrationError  = nIntegrationError 
     PressExceededParams[0] = np.append(PressExceededParams[0],-666.0)
     PressExceededParams[1] = np.append(PressExceededParams[1],-666.0)
     PressExceededParams[2] = np.append(PressExceededParams[2],-666.0)
     PressExceededParams[3] = np.append(PressExceededParams[3],-666.0)
     IntegrationErrorParams[0] = np.append(IntegrationErrorParams[0],-666.0)
     IntegrationErrorParams[1] = np.append(IntegrationErrorParams[1],-666.0)
     IntegrationErrorParams[2] = np.append(IntegrationErrorParams[2],-666.0)
     IntegrationErrorParams[3] = np.append(IntegrationErrorParams[3],-666.0)

    return(nSigmaCrit, nTlim, SigmaCritParams, TlimParams, nPressExceeded, nIntegrationError, PressExceededParams, IntegrationErrorParams)

if __name__ == '__main__':
    
    tags = enum('READY', 'DONE', 'EXIT', 'START','PARAM', 'NCD', 'SendNCD')
    TAGS=('READY', 'DONE', 'EXIT', 'START','PARAM', 'NCD', 'SendNCD')
    comm = MPI.COMM_WORLD # Communicator
    size = comm.size      # Number of processes
    rank = comm.rank      # this process
    status = MPI.Status()
    if size < 2:
        print "ERROR: not enough workers!"
        exit(0)


    #
    # number of non-converged runs
    nSigmaCrit = 0
    nTlim = 0
    nPressExceeded = 0
    nIntegrationError = 0


    #parameter values for non-converged runs
    SigmaCritParams=[ np.empty(shape=0), np.empty(shape=0), np.empty(shape=0), np.empty(shape=0)]
    TlimParams= [ np.empty(shape=0), np.empty(shape=0), np.empty(shape=0), np.empty(shape=0)]
    PressExceededParams=[ np.empty(shape=0), np.empty(shape=0), np.empty(shape=0), np.empty(shape=0)]
    IntegrationErrorParams=[ np.empty(shape=0), np.empty(shape=0), np.empty(shape=0), np.empty(shape=0)]

    # make directories where final results are stored
    if rank == 0:
        os.makedirs(RisultatiMultipli)
        os.makedirs(Database)
        os.makedirs(LogFiles)
        os.makedirs(Thumbnails)
    #
    # open a logger (one each task==rank)
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s %(levelname)s %(message)s',
                        filename=workDir+"/run_"+str(rank)+".log",
                        filemode='w')

    if rank == 0:
        

        # copying Michele stuff in working dir
        shutil.copy(template_dir+"/tAtmo/"+"tAtmo_rc",workDir)
        shutil.copytree(src_dir+'atmosphereGeometry',workDir+'/atmosphereGeometry')
        shutil.copytree(src_dir+'atmosphereLib',workDir+'/atmosphereLib')
        shutil.copytree(src_dir+'pipeline_interface',workDir+'/pipeline_interface')

        #beginning
        restart_interval = 300 # in seconds
        stop_time = 6*3600 #in seconds
        starttime= time()
        oldtime = starttime
        simulation_index = 0
        num_workers = size - 1
        closed_workers = 0
        computed_models_file=workDir+"/executed.out"
        with open(computed_models_file, "w") as myfile:
            myfile.write("# p, ecc, obl, dist\n")
        logging.info("Master starting with %d workers" % num_workers)
        print ("Master starting with %d workers" % num_workers)
        stop_file=workDir+"/stop"
        input_filename=workDir+"/input.txt"
        restart_file=workDir+"/restart.dat"
        running_input=workDir+'/input.%s.dat' % os.getpid()

        try:
            fd = os.open(restart_file, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        except OSError, e:
            if e.errno == 17:
                logging.info("%s" % e)
                shutil.copy(restart_file,running_input)
            else:
                raise
        else:
            shutil.copy(input_filename,running_input)
            shutil.copy(input_filename,restart_file)


        try:
            infile = open(running_input,"r")
        except IOError:
            print "Cannont open Input File"
            exit() ##### verify if it exists a proper way to close MPI

        for line in infile:
            print '\n\n', line, '\n'
            data = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
            source = status.Get_source()
            tag = status.Get_tag()
            print 'Task ', rank, ' received ', TAGS[tag], ' from ', source
            if tag == tags.READY: # Only at first loop
                comm.send([simulation_index,line], dest=source, tag=tags.START)
                logging.info("tags.READY Sending simulation %d to worker %d" % (simulation_index, source))
                print 'Task ',rank,' sent START to ', source
                simulation_index += 1
            elif tag == tags.DONE:
                results = data # collect results from worker
                with open(computed_models_file, "a") as myfile:
                        myfile.write(results)

                data = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)               
                source = status.Get_source()
                tag = status.Get_tag()
            
                print 'Data: ',data
                print data[8]
                print 'Task ', rank, ' received ', TAGS[tag], ' from ', source
                if data[8]<0.:
                    print 'collecting non converged data, task ',rank
                    nSigmaCrit, nTlim, SigmaCritParams, TlimParams, nPressExceeded, nIntegrationError, PressExceededParams, IntegrationErrorParams = collect_non_converged_models_data(nSigmaCrit,nTlim,SigmaCritParams,TlimParams,nPressExceeded,nIntegrationError,PressExceededParams,IntegrationErrorParams,data)

                if time() - oldtime > restart_interval:
                    logging.info("Write restart file ")
                    write_restart_file(input_filename,computed_models_file,restart_file)
                    oldtime = time()
                if os.path.isfile(stop_file):
                    break
                if time() - starttime > stop_time:
                    break
                simulation_index += 1
                logging.info("tags.DONE Got data from worker %d: %s" % (source,results))
                
                comm.send([simulation_index,line], dest=source, tag=tags.START) #send new data to worker
                print 'Task ',rank,' sent START to ', source
            elif tag == tags.EXIT: # ERROR: we should not be here!!!!
                logging.error(" tags.EXIT Worker %d exited." % source)
                closed_workers+=1
    
        #
        ## INPUT file completed. Collect running worker results end ask them to exit
        #
        while closed_workers < num_workers:
            data = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
            source = status.Get_source()
            tag = status.Get_tag()
            print 'Task ', rank, 'received ', TAGS[tag], 'from ', source

            if tag == tags.READY: # Should not happen unless total simulations < numer_of_workers
                comm.send(None, dest=source, tag=tags.EXIT)
                print 'Task ', rank, ' sent EXIT to ', source
            elif tag == tags.DONE: # collect results from running workers and ask to exit
                results = data
                with open(computed_models_file, "a") as myfile:
                    myfile.write(results)
                logging.info("EXIT: Got data from worker %d" % source)                

                data = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)               
                source = status.Get_source()
                tag = status.Get_tag()
            
                print 'Data: ',data
                print data[8]
                print 'Task ', rank, ' received ', TAGS[tag], ' from ', source

                if data[8]<0.:
                    print 'collecting non converged data, task ',rank
                    nSigmaCrit, nTlim, SigmaCritParams, TlimParams, nPressExceeded, nIntegrationError, PressExceededParams, IntegrationErrorParams = collect_non_converged_models_data(nSigmaCrit,nTlim,SigmaCritParams,TlimParams,nPressExceeded,nIntegrationError,PressExceededParams,IntegrationErrorParams,data)



                comm.send(None, dest=source, tag=tags.EXIT)
                print 'Task ', rank, ' sent EXIT to ', source
            elif tag == tags.EXIT: # Collect exit reply from workers
                logging.info("EXIT: Worker %d exited." % source)
                closed_workers+=1



        #
        ## END and close inputfile
        #
        infile.close()
    else: # working tasks
        name = MPI.Get_processor_name()
        logging.info("I am a worker with rank %d on %s." % (rank, name))
        local_simulation_index = 0
        while True:
            if local_simulation_index == 0:
                local_simulation_index += 1
                comm.send(None, dest=0, tag=tags.READY)
                print 'Task ', rank, 'sent READY to 0'
                inputdata = comm.recv(source=0, tag=MPI.ANY_TAG, status=status)
                tag = status.Get_tag()
                print 'Task ', rank, 'received ', TAGS[tag],' from 0'
                logging.info("Receive simulation on worker %d from  0"  % (rank))
                if tag == tags.START:
                    logging.debug("%d,0: begin computation" % (rank))
                    # number of non-converged runs - local data
                    nSigmaCritL = 0
                    nTlimL = 0
                    nPressExceededL = 0
                    nIntegrationErrorL = 0
                    #parameter values for non-converged runs - local data
                    SigmaCritParamsL=[ np.empty(shape=0), np.empty(shape=0), np.empty(shape=0), np.empty(shape=0)]
                    TlimParamsL= [ np.empty(shape=0), np.empty(shape=0), np.empty(shape=0), np.empty(shape=0)]
                    PressExceededParamsL=[ np.empty(shape=0), np.empty(shape=0), np.empty(shape=0), np.empty(shape=0)]
                    IntegrationErrorParamsL=[ np.empty(shape=0), np.empty(shape=0), np.empty(shape=0), np.empty(shape=0)]

                    Parameter_set = make_input_parameters(inputdata,Parameter_set)
                    ####################
                    # running the code:#
                    ####################
                    nSigmaCritL,nTlimL,SigmaCritParamsL,TlimParams,nPressExceededL, nIntegrationErrorL, PressExceededParamsL,IntegrationErrorParamsL,exitValueL=esoclimi_emulate(Parameter_set,nSigmaCritL,nTlimL,SigmaCritParamsL,TlimParamsL,nPressExceededL,nIntegrationErrorL,PressExceededParamsL,IntegrationErrorParamsL)

                    # sending back results:
                    comm.send(inputdata[1], dest=0, tag=tags.DONE)
                    print 'Task ', rank, 'sent DONE to 0'

                    non_converging_data = [nSigmaCritL,nTlimL,SigmaCritParamsL,TlimParamsL,nPressExceededL,nIntegrationErrorL,PressExceededParamsL,IntegrationErrorParamsL,exitValueL]
                    comm.send(non_converging_data, dest=0, tag=tags.NCD)
                    print 'Task ',rank, 'sent NCD to 0'


                elif tag == tags.EXIT:
                    comm.send(None, dest=0, tag=tags.EXIT) # chiudi task mpi e esci
            else:
                local_simulation_index += 1
                inputdata = comm.recv(source=0, tag=MPI.ANY_TAG, status=status)
                tag = status.Get_tag()
                if tag == tags.START:
                    # number of non-converged runs - local data
                    nSigmaCritL = 0
                    nTlimL = 0
                    nPressExceededL = 0
                    nIntegrationErrorL = 0
                    #parameter values for non-converged runs - local data
                    SigmaCritParamsL=[ np.empty(shape=0), np.empty(shape=0), np.empty(shape=0), np.empty(shape=0)]
                    TlimParamsL= [ np.empty(shape=0), np.empty(shape=0), np.empty(shape=0), np.empty(shape=0)]
                    PressExceededParamsL=[ np.empty(shape=0), np.empty(shape=0), np.empty(shape=0), np.empty(shape=0)]
                    IntegrationErrorParamsL=[ np.empty(shape=0), np.empty(shape=0), np.empty(shape=0), np.empty(shape=0)]
                    Parameter_set = make_input_parameters(inputdata,Parameter_set)
                    ####################
                    # running the code:#
                    ####################
                    nSigmaCritL,nTlimL,SigmaCritParamsL,TlimParams,nPressExceededL, nIntegrationErrorL, PressExceededParamsL,IntegrationErrorParamsL,exitValueL=esoclimi_emulate(Parameter_set,nSigmaCritL,nTlimL,SigmaCritParamsL,TlimParamsL,nPressExceededL,nIntegrationErrorL,PressExceededParamsL,IntegrationErrorParamsL)

                    #sending back results:
                    comm.send(inputdata[1], dest=0, tag=tags.DONE)
                    print 'Task ', rank, 'sent DONE to 0'
#get ncd

                    non_converging_data = [nSigmaCritL,nTlimL,SigmaCritParamsL,TlimParamsL,nPressExceededL,nIntegrationErrorL,PressExceededParamsL,IntegrationErrorParamsL,exitValueL]
                    comm.send(non_converging_data, dest=0, tag=tags.NCD)
                    print 'Task ',rank, 'sent NCD to 0'

                elif tag == tags.EXIT:
                        comm.send(None, dest=0, tag=tags.EXIT)
                        print 'Task ', rank, 'sent EXIT to 0'
                        break

#        comm.send(non_converging_data, dest=0, tag=tags.EXIT)
        print 'Task ', rank, ' sent EXIT to 0'

    if rank == 0:
        print 'ORA SCRIVO I SOMMARI'
        write_non_converged_models(nSigmaCrit,nTlim,simulation_index,SigmaCritParams,TlimParams,nPressExceeded,nIntegrationError,PressExceededParams,IntegrationErrorParams)

        exit()    

    exit()
