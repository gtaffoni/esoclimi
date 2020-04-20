#!/usr/bin/env python

#  serial.py
#  Esoclimi
#
#  Created by giuliano taffoni on 09/10/17.
#  Copyright  2017 Guliano Taffoni. All rights reserved.

'''
    serial code to run exoclime driver
'''

import sys
import os
import shutil
from posix import system
import logging
from time import  time
import numpy as np
import time
from exoclime import *
import ConfigParser
#from runEBM import *

Parameters = {'simtype': "Std", 'version': "1.1.03", 'planet': "EARTH" }



def make_input_parameters(_data,parameters):
    '''
        Convert input from rank 0 into set of parametes
        WARNING: TO BE MODIFIED WHEN CHANGING PARAMETERS SPACE EXPLORATION
        Planet flatning da aggiunegre al fortran
        salvare fits anche per confergenza  nulla (non caso numerico)
        '''
    #TODO ADD date from 2000 for files on minutes
    input_params=np.fromstring(_data[1],  sep=' ')
    parameters['TOAalbfile']        = 'CCM_RH60/' + ''.join(_data[1].split()[-2:-1])
    parameters['OLRfile']           = 'CCM_RH60/' + ''.join(_data[1].split()[-1:])
    parameters['p_CO2_P']           = input_params[7]   #CO2 partial pressure IN PPVM
    parameters['CO2_Earth_ratio']   = input_params[6]   #the same, in Earth ratio (for output) TODO: CALCOLIMOCELA!!!!!
    parameters['fo_const']          = input_params[5]
    parameters['gg']                = input_params[4]
    parameters['dist']              = input_params[3]    # semi-major axis of planet orbit
    parameters['obl']               = input_params[2]     # planet axis inclination
    parameters['ecc']               = input_params[1]     # eccentricity of planet orbit
    parameters['p']                 = input_params[0]       # pressure in EARTH units (atm)
    parameters['number']            = _data[0]
    parameters['data']              = _data[1]   #MJD YYYY-MM-DD (ISO)
    return(parameters)



if __name__ == '__main__':
    
    
    #########################################
    # Variables initialization              #
    #########################################
    config = ConfigParser.ConfigParser()
    config.read("Config.ini")
    code_src_dir = config.get("Global", "Code Root Dir")
    # PYTHON PATH
    sys.path.append(code_src_dir)
    # Work directory
    workDir = os.getcwd()
    # Globals
    RisultatiMultipli = "RisultatiMultipli"
    Database          = "Database"
    Thumbnails        = "Thumbnails"
    LogFiles          = "LogFiles"
    Risultati         = "Risultati"
    Broken            = "Broken"
    Src               = "Src"
    #
    # Non convergin models
    #
    N_non_converging=np.zeros(5,dtype=int)
    #   nSigmaCrit (run)  => 0
    #   nTlim (snow)      => 1
    #   nPressExceeded    => 2
    #   nIntegrationError => 3
    #   uncompleted (256) => 4
    #parameter values for non-converged models
    #
    SigmaCritParams=[ np.empty(shape=0) ]
    TlimParams= [ np.empty(shape=0)]
    PressExceededParams=[ np.empty(shape=0)]
    IntegrationErrorParams=[np.empty(shape=0)]
    exitValueL = 0
    #
    # make directories where final results are stored
    #
    os.makedirs(RisultatiMultipli)
    os.makedirs(Database)
    os.makedirs(LogFiles)
    os.makedirs(Thumbnails)
    os.makedirs(Broken)

    #########################################
    # open a logger                         #
    #########################################
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s %(levelname)s %(message)s',
                        filename=workDir+"/run.log",
                        filemode='w')
    #########################################
    # Begin                                 #
    #########################################
    #
    # input data format:
    # p ecc obl dist gg fo_const CO2_Earth_ratio p_CO2_P OLRfile (only name not dir) TOAalbfile (only name not dir)
    #
    #data="0.01 0.0 23.43929 0.8 0 0.1 1.0 380 ALB_g1_rh60_co2x10.txt OLR_g1_rh60_co2x10.txt"   # integration error (exit -100)
    #data="0.018 0.1 23.439 1.5 1.0 0.3 1.0 380 ALB_g1_rh60_co2x10.txt OLR_g1_rh60_co2x10.txt" # Snowball converged (exit 3)
    #data="0.017783 0.6 23.439290 1.0 4 0.70 1.0 380 ALB_g1_rh60_co2x10.txt OLR_g1_rh60_co2x10.txt" # warm-hot (exit 2)
    #data="0.017783 0.00 30.00 0.9 0 0.70 1.0 380 ALB_g1_rh60_co2x10.txt OLR_g1_rh60_co2x10.txt" # warm (exit 1)
    data = "1 0.01671022 23.439290 1.0 4 0.70 1. 380 ALB_g1_rh60_co380_ch1.8.txt OLR_g1_rh60_co380_ch1.8.txt" # Terra!!! (exit 1)
    #data="0.017783 0.70 30.0 0.8 0 0.70 1.0 380 ALB_g1_rh60_co2x10.txt OLR_g1_rh60_co2x10.txt" # Runaway GreenHouse (exit -1)
    #data="0.01 0.8 0.0 0.8 0 0.10 1.0 380 ALB_g1_rh60_co2x10.txt OLR_g1_rh60_co2x10.txt" # pressure exceeded (exit -2)
    
    simulation_index=1
    inputdata=[simulation_index,data]
    Parameters = make_input_parameters(inputdata,Parameters)
    Parameters['number']=simulation_index
    logging.info("Master starting....")

    ######################################################
    # prepare, compile and execute code, store results   #
    ######################################################
    #
    # Make directory structure for code execution
    #
    try:
        make_work_area(workDir,code_src_dir,Risultati,Parameters)
    except:
        logging.error(sys.exc_info()[0])
        exitValueL = 256
        pass

#
# prepare compile and execute
#
    if not exitValueL == 256:
        try:
            exitValueL = exoclime(Parameters, workDir, code_src_dir, Risultati, emulate=False)
        except:
            logging.error(sys.exc_info()[0])
            exitValueL = 256
            pass

    if exitValueL == -200:
        print 'WARNING, CATASTROPHIC EXIT VALUE FOUND, ELABORATION STOPPED'
        print 'Parameters: ', Parameters['data']
        comm.free()
        exit(-200)
    #
    # Archive results
    #
    if exitValueL == 256:
        try:
            archive_broken_simulations(Parameters, workDir, Broken)
        except:
            logging.error("unable to archive broken simulations: %s",sys.exc_info()[0])
            pass
    else:
        try:
            archive_exoplanet_data(Parameters, workDir,RisultatiMultipli,LogFiles,Risultati)
        except:
            logging.error("unable to archive  simulations: %s",sys.exc_info()[0])
            pass


    print ("Exit Value: %d" % exitValueL)


