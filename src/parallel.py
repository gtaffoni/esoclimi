#!/usr/bin/env python

#  parallel.py
#  Esoclimi
#
#  Created by giuliano taffoni on 01/10/17.
#  Copyright  2017 Guliano Taffoni. All rights reserved.

'''
    parallel code to run exoclime driver
    '''

import sys 
import os
import shutil
from posix import system
import logging
from random import randint
from time import sleep
from time import  time
import numpy as np
import time
from exoclime import *
from mpi4py import MPI
import ConfigParser
import difflib


Parameters = {'simtype': "Std", 'version': "1.1.03", 'planet': "EARTH" }

def enum(*sequential, **named):
    '''
        simple way to emulate enumerate in python taken from the web
        
    '''
    enums = dict(zip(sequential, range(len(sequential))), **named)
    return type('Enum', (), enums)


def write_intermediate_results(sim_index,input,computed,tmp_computed,snowball, tmp_snowball, runaway, tmp_runaway, pressure, tmp_pressure, integration, tmp_integration, nonconverging_file, uncompleted, tmp_uncompleted, n,out):
    '''
               write a restart file  and consolidate the outputs
               n = non converginf array
    '''
    
    #
    # Backup: Necessary in case of crash.
    #
    shutil.copyfile(computed,computed+".bak")
    shutil.copyfile(out,out+".bak")
    #
    # Update computed models file with data in temporary_computed (the computed file of actual run)
    #
    merge_files(computed,tmp_computed)
    #
    # Update list of non converging models
    #
    for _file in [tmp_snowball,tmp_runaway,tmp_pressure,tmp_pressure]:
        with open(nonconverging_file, 'a') as outfile:
            with open(_file) as infile:
                for line in infile:
                    outfile.write(line)
    #
    # Update other files.
    #
    merge_files(snowball,tmp_snowball)
    merge_files(runaway, tmp_runaway)
    merge_files(pressure, tmp_pressure)
    merge_files(integration, tmp_integration)
    merge_files(uncompleted, tmp_uncompleted)
    #
    # Update list of non converging models
    #
    #   nSigmaCrit (snow) => n0
    #   nTlim  (run) => n1
    #   nPressExceeded => n2
    #   nIntegrationError =>  n3
    o = ConfigParser.SafeConfigParser()
    o.read(out)
    o.set("Main program","Simulation index",str(sim_index))
    o.set("Divergent Models","Snowball",str(n[1]))
    o.set("Divergent Models","Runaway greenhouse",str(n[0]))
    o.set("Divergent Models","Pressure execeded",str(n[2]))
    o.set("Divergent Models","Integration error",str(n[3]))
    o.set("Uncompleted Models","Number",str(n[4]))
    with open(out, 'wb') as configfile:
        o.write(configfile)
    return

#(data, tmp_computed_models_file, tmp_snowball_file,tmp_runaway_greenhouse_file,
# tmp_press_exceeded_file, tmp_integration_error_file,tmp_uncompleted_file,
# N_non_converging)

def collect_results(_data, cmf, sf, rgf, pef, ief, ucf, n):
    '''
        collect the results from a worker and save temporary data
        '''
    results = _data[1] # results from worker
    #
    # Computed models file
    #
    with open(cmf, "a") as tmp_file:
        tmp_file.write(results)
    #
    # Divergence File
    #
    #   nSigmaCrit (run)  => 0
    #   nTlim (snow)      => 1
    #   nPressExceeded    => 2
    #   nIntegrationError => 3
    #   nUncompleted      => 4
    if np.abs(_data[2] + 100) < 0.001 :   # Integration Error
        with open(ief, "a") as tmp_file:
            tmp_file.write(results)
        n[3] += 1
    elif np.abs(_data[2] + 2.0) < 0.001 : #pressure exceeded (should not happen!)
        with open(pef, "a") as tmp_file:
            tmp_file.write(results)
        n[2] += 1
    elif np.abs(_data[2] + 1.0) < 0.001 : #Runaway GreenHouse
        with open(rgf, "a") as tmp_file:
            tmp_file.write(results)
        n[0] += 1
    elif np.abs(_data[2] + 0.5) < 0.001: #SnowBall
        with open(sf, "a") as tmp_file:
            tmp_file.write(results)
        n[1] += 1
    elif np.abs(_data[2] - 256) < 0.001: #Uncompleted (256)
        with open(ucf, "a") as tmp_file:
            tmp_file.write(results)
        n[4] += 1
    return n

def write_restart_file(_input,_computed,_restart):
    file1=open(_computed)
    file2=open(_input)
    file3=open(_restart,"w")
    # Read lines from files
    lines1 = file1.read().strip().splitlines()
    lines2 = file2.read().strip().splitlines()
    diff = difflib.unified_diff(lines1, lines2, fromfile='file1', tofile='file2', lineterm='', n=0)
    lines = list(diff)[2:]
    added   = [line[1:] for line in lines if line[0] == '+']
    removed = [line[1:] for line in lines if line[0] == '-']

    number_of_models_to_compute=0
    for line in added:
        if line not in removed:
            file3.writelines(line + "\n")
            number_of_models_to_compute += 1
    if  number_of_models_to_compute != (len(lines2)- len(lines1)):
        return number_of_models_to_compute
    if  number_of_models_to_compute == 0:
        return 0
    return 1

def merge_files(f1,f2):
    '''
        Append the content of f2 to f1 and empty f2
        f1 and f2 are file names
        '''
    with open(f1, 'a') as outfile:
        with open(f2) as infile:
            for line in infile:
                outfile.write(line)
    open(f2, 'w').close()
    return

def make_output_file_data_structure(filename):
    o = ConfigParser.SafeConfigParser()
    o.add_section("Main program")
    o.add_section("Divergent Models")
    o.add_section("Uncompleted Models")
    o.set("Main program","Simulation index","0")
    o.set("Divergent Models","Snowball","0")
    o.set("Divergent Models","Runaway greenhouse","0")
    o.set("Divergent Models","Pressure execeded","0")
    o.set("Divergent Models","Integration error","0")
    o.set("Uncompleted Models","Number","0")
    file = open(filename,"w")
    o.write(file)
    file.close()
    del o
    return

    #def make_input_parameters(_data,parameters):
    #'''
    #    Convert input from rank 0 into set of parametes
    #    WARNING: TO BE MODIFIED WHEN CHANGING PARAMETERS SPACE EXPLORATION
    #    '''
def make_input_parameters(_data,parameters):
    '''
        Convert input from rank 0 into set of parametes
        WARNING: TO BE MODIFIED WHEN CHANGING PARAMETERS SPACE EXPLORATION
        '''
    input_params=np.fromstring(_data[1],  sep=' ')
    parameters['TOAalbfile']        = 'CCM_RH60/' + ''.join(_data[1].split()[-2:-1])
    parameters['OLRfile']           = 'CCM_RH60/' + ''.join(_data[1].split()[-1:])
    parameters['p_CO2_P']           = input_params[7]   #CO2 partial pressure IN PPVM
    parameters['CO2_Earth_ratio']   = input_params[6]   #the same, in Earth ratio (for output)
    parameters['fo_const']          = input_params[5]
    parameters['gg']                = input_params[4]
    parameters['dist']              = input_params[3]    # semi-major axis of planet orbit
    parameters['obl']               = input_params[2]     # planet axis inclination
    parameters['ecc']               = input_params[1]     # eccentricity of planet orbit
    parameters['p']                 = input_params[0]       # pressure
    parameters['number']            = _data[0]
    parameters['data']              = _data[1]
    return(parameters)

if __name__ == '__main__':
    
    tags = enum('READY', 'DONE', 'EXIT', 'START','PARAM', 'NCD', 'NCDS', 'GO')
    TAGS=('READY', 'DONE', 'EXIT', 'START','PARAM', 'NCD', 'NCDS', 'GO')
    comm = MPI.COMM_WORLD # Communicator
    size = comm.size      # Number of processes
    rank = comm.rank      # this process
    status = MPI.Status()
    if size < 2:
        print "ERROR: not enough workers!"
        exit(0)
    #######################################
    #   Initializations                   #
    #######################################
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
    runs              = "Runs"
    #
    # Non convergin models
    #
    N_non_converging=np.zeros(5,dtype=int)
    #   nSigmaCrit (run)  => 0
    #   nTlim (snow)      => 1
    #   nPressExceeded    => 2
    #   nIntegrationError => 3
    #   uncompleted (256) => 4
    #
    # open a logger (one each task==rank)
    logging_file_name=workDir+"/run_"+str(rank)+".log"
    logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s %(message)s',
                    filename=logging_file_name,
                    filemode='w')

    #########################################
    #   Parallel distribution               #
    #########################################
    if rank == 0:  ####      MASTER     #####
        # Initializations
        restart_interval = int(config.get("Checkpoint", "Restart Interval"))
        n_of_runs_before_restart = int(config.get("Checkpoint", "Number of runs before restart"))
        stop_time = int(config.get("Checkpoint", "Stop Time"))
        #
        stop_file = workDir+"/"+config.get("Checkpoint", "Stop File")
        cpu_file = workDir+"/"+config.get("Code", "CPU File")
        computed_models_file = workDir+"/"+config.get("Code", "Executed runs")
        input_filename = workDir+"/"+config.get("Code", "Input File")
        output_filename = workDir+"/"+config.get("Code", "Output File")
        restart_file = workDir+"/"+config.get("Code", "Restart File")
        nonconverging_file = workDir+"/"+config.get("Code", "Non Converging Models File")
        uncompleted_file = workDir+"/"+config.get("Code", "Uncompleted Models File")
        snowball_file = workDir+"/"+config.get("Results", "SnowBall File")
        runaway_greenhouse_file = workDir+"/"+config.get("Results", "RunawayGreenhouse File")
        press_exceeded_file = workDir+"/"+config.get("Results", "PressExceeded File")
        integration_error_file = workDir+"/"+config.get("Results", "IntegrationError File")
        #
        running_input=workDir+'/input.%s.dat' % os.getpid() # TODO: use tmp random files
        #
        tmp_snowball_file = snowball_file+"_tmp" # TODO: use tmp random files ( maybe hidden)
        open(tmp_snowball_file, 'w').close()
        #
        tmp_runaway_greenhouse_file = runaway_greenhouse_file+"_tmp" # TODO: use tmp random files ( maybe hidden)
        open(tmp_runaway_greenhouse_file, 'w').close()
        #
        tmp_press_exceeded_file = press_exceeded_file+"_tmp" # TODO: use tmp random files ( maybe hidden)
        open(tmp_press_exceeded_file, 'w').close()
        #
        tmp_integration_error_file = integration_error_file+"_tmp" # TODO: use tmp random files ( maybe hidden)
        open(tmp_integration_error_file, 'w').close()
        #
        tmp_computed_models_file = computed_models_file+"_tmp"  # TODO: use tmp random files ( maybe hidden)
        open(tmp_computed_models_file, 'w').close()
        #
        # DEBUG: maybe it is redundant but avoid to copy an open file at checkpoint
        tmp_running_input=running_input+"_tmp"  # TODO: use tmp random files ( maybe hidden)
        
        open(tmp_running_input, 'w').close()
        #
        tmp_uncompleted_file = uncompleted_file+"_tmp"
        open(tmp_uncompleted_file, 'w').close()
        #
        starttime= time()
        oldtime = starttime
        simulation_index = 0
        num_workers = size - 1
        num_computed = 0
        closed_workers = 0
        #
        create_directory_structure = False
        output = ConfigParser.ConfigParser() # to write and read results
        #
        #   Open file with computed models and append header only if file does not exists.
        #
        #if not (os.path.isfile(computed_models_file)):
        #    with open(computed_models_file, "w") as tmp_file:
        #        tmp_file.write("# type, p, ecc, obl, dist, gg, fo_const\n") # NOTE: to update when adding parameters
        #
        #   Open  executed model file if there is one, create restart file and read simulation index, and other counters from output file
        #
        if os.path.isfile(stop_file):
            print "WARNING: I found a stop file, I assume that either"
            print "the execution is completed (no more model to evaluate) or"
            print ("the stop file %s has been created by the user" % stop_file)
            comm.Abort()
            exit(0)
                
        try:
            fd = os.open(computed_models_file, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        except OSError, e:
            if e.errno == 17: #EXISTS
                write_restart_file_err = write_restart_file(input_filename, computed_models_file,restart_file)
                if write_restart_file_err == 1:
                    shutil.copy(restart_file,running_input)
                    shutil.copy(restart_file,tmp_running_input)
                    output.read(output_filename)
                    simulation_index = int(output.get("Main program", "simulation index"))
                    N_non_converging[0]=int(output.get("Divergent Models","Snowball"))
                    N_non_converging[1]=int(output.get("Divergent Models","Runaway greenhouse"))
                    N_non_converging[2]=int(output.get("Divergent Models","Pressure execeded"))
                    N_non_converging[3]=int(output.get("Divergent Models","Integration error"))
                    N_non_converging[4]=int(output.get("Uncompleted Models","Number"))
                elif write_restart_file_err == 0 :
                    print "WARNING: Execution completed, no more model to evaluate"
                    open(stop_file, 'a').close()
                    comm.Abort()
                    exit(0)
                else:
                    print "ERROR: Cannot Create Restart file"
                    comm.Abort()
                    exit()
            else:
                raise
        else:
            shutil.copy(input_filename,running_input)
            shutil.copy(running_input,tmp_running_input)
            shutil.copy(input_filename,restart_file)
            make_output_file_data_structure(output_filename)
            create_directory_structure = True

        try:
            infile = open(running_input,"r")
        except IOError:
            print "Cannont open Input File"
            comm.Abort(1)
            exit(0)


        logging.info("Master starting with %d workers from simulation number %d" % (num_workers, simulation_index))
        print ("Master starting with %d workers from simulation number %d" % (num_workers, simulation_index))
        if create_directory_structure :
        #
        # create directories where final results are stored
        #
            os.makedirs(RisultatiMultipli)
            os.makedirs(Database)
            os.makedirs(Broken)
            os.makedirs(LogFiles)
            os.makedirs(Thumbnails)
            logging.debug("RANK 0 Create dir strucrture")
        #
        #   BEGIN COMPUTATION LOOP
        #
        for line in infile:
            data = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
            source = status.Get_source()
            tag = status.Get_tag()
            if tag == tags.READY: # Only at first loop
                logging.info("tags.READY Sending simulation %d to worker %d" % (simulation_index, source))
                comm.send([simulation_index,line], dest=source, tag=tags.START)
                simulation_index += 1
            elif tag == tags.DONE:
                N_non_converging = collect_results(data, tmp_computed_models_file, tmp_snowball_file,tmp_runaway_greenhouse_file,
                                                   tmp_press_exceeded_file, tmp_integration_error_file,tmp_uncompleted_file,
                                                   N_non_converging)
                logging.info("tags.DONE Got data from worker %d. Simulation index: %d" % (source,data[0]))
                logging.debug("tags.DONE Got data from worker %d. Simulation index: %d, input: %s, err value: %d" % (source,data[0],data[1],data[2]))
            #
            #       CHECK IF IT IS NECESSARY TO MAKE A CHECKPOINT
            #
                if time() - oldtime > restart_interval or simulation_index%n_of_runs_before_restart==0:
                    logging.info("tags.DONE Master Write restart file")
                    #
                    # WARNING: simulation index is always greater than the number of completed runs of a factor num_workers-1.
                    #
                    write_intermediate_results(simulation_index-(num_workers-1), input_filename, computed_models_file,tmp_computed_models_file,
                                        snowball_file, tmp_snowball_file, runaway_greenhouse_file, tmp_runaway_greenhouse_file,
                                        press_exceeded_file, tmp_press_exceeded_file, integration_error_file, tmp_integration_error_file,
                                        nonconverging_file,  uncompleted_file, tmp_uncompleted_file, N_non_converging,output_filename)
                    oldtime = time()
                if os.path.isfile(stop_file):
                    logging.info("tags.DONE Master stop file found: closing the simulation.")
                    comm.send(None, dest=source, tag=tags.EXIT)
                    break
                if time() - starttime > stop_time:
                    logging.info("Stopping time limit exceeded ")
                    comm.send(None, dest=source, tag=tags.EXIT)
                    break
            #
            #    RUN NEW SIMULATION
            #
                logging.info("tags.READY Sending simulation %d to worker %d" % (simulation_index, source))
                comm.send([simulation_index,line], dest=source, tag=tags.START) #send new data to worker
                simulation_index += 1
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
            logging.info("Closed worker number: %d (worker=%d)" % (closed_workers,source))

            if tag == tags.READY: # Should not happen unless total simulations < numer_of_workers
                logging.info("Closing ==> tag.READY: close worker %d" % source)
                comm.send(None, dest=source, tag=tags.EXIT)
            elif tag == tags.DONE: # collect results from running workers and ask to exit
                N_non_converging = collect_results(data, tmp_computed_models_file, tmp_snowball_file,tmp_runaway_greenhouse_file,
                                                    tmp_press_exceeded_file, tmp_integration_error_file, tmp_uncompleted_file,
                                                    N_non_converging)
                #
                # MAKE A CHECKPOINT EACH TIME WE FINISH a TASK
                #
                write_intermediate_results(simulation_index, input_filename, computed_models_file,tmp_computed_models_file,
                           snowball_file, tmp_snowball_file, runaway_greenhouse_file, tmp_runaway_greenhouse_file,
                           press_exceeded_file, tmp_press_exceeded_file, integration_error_file, tmp_integration_error_file,
                           nonconverging_file,  uncompleted_file, tmp_uncompleted_file, N_non_converging,output_filename)
                logging.info("Closing ==> tag.DONE: Got data from worker %d. Simulation index: %d" % (source, data[0]))
                comm.send(None, dest=source, tag=tags.EXIT)
            elif tag == tags.EXIT: # Collect exit reply from workers
                logging.info("Closing ==> tag.EXIT Worker %d exited (n_closed=%d)." % (source,closed_workers))
                closed_workers+=1
                logging.info("Closing ==> tags.EXIT Master Write restart file")

#
## END and close inputfile
#
        infile.close()
        #
        # Check if we need to write final  data: paranoid check
        #
        tmp_computed_models_file_index=0
        with open(tmp_computed_models_file) as tmp_computed_models_file_fd:
            for line in tmp_computed_models_file_fd:
                tmp_computed_models_file_index += 1
        if tmp_computed_models_file_index != 0:
            logging.info("Closing: Master Write  %d results  to files",tmp_computed_models_file_index)
            write_intermediate_results(simulation_index, input_filename, computed_models_file,tmp_computed_models_file,
                            snowball_file, tmp_snowball_file, runaway_greenhouse_file, tmp_runaway_greenhouse_file,
                            press_exceeded_file, tmp_press_exceeded_file, integration_error_file, tmp_integration_error_file,
                            nonconverging_file, uncompleted_file, tmp_uncompleted_file, N_non_converging,output_filename)
        #
        # Print results
        #
        print 'nSigmaCrit (Runaway Greenhouse), nTlim (Snowball), nPressExceeded (out of range), nIntegrationError (stepsize too small): ', N_non_converging[0],  N_non_converging[1], N_non_converging[2], N_non_converging[3]
        print 'Fractions: ', 1.0*N_non_converging[0]/simulation_index, 1.0*N_non_converging[1]/simulation_index, 1.0*N_non_converging[2]/(simulation_index), 1.0*N_non_converging[3]/(simulation_index)
        print '\n Overall number and fraction of non-converged inhabitable cases: ', sum(N_non_converging), 1.0* sum(N_non_converging) / (simulation_index+1)
        print '\n Number of uncompleted runs: ', N_non_converging[4]
        print '\n Total number of runs: ',simulation_index




        

    else: # working tasks (slaves)
        name = MPI.Get_processor_name()
        logging.info("I am a worker with rank %d on %s." % (rank, name))
        local_simulation_index = 0
        while True:
            if local_simulation_index == 0:
                local_simulation_index += 1
                exitValueL = 0
                comm.send(None, dest=0, tag=tags.READY)
                inputdata = comm.recv(source=0, tag=MPI.ANY_TAG, status=status)
                tag = status.Get_tag()
                logging.info("Receive message %s on worker %d from  0"  % (tag,rank))
                if tag == tags.START:
                    logging.debug("START on rank %d, local sim %d: begin computation" % (rank,local_simulation_index))
                    
                    Parameters = make_input_parameters(inputdata,Parameters)
                    logging.debug("Simulations %d => Check exitValueL=%d" % (Parameters['number'] ,exitValueL))
                    ######################################################
                    # prepare, compile and execute code, store results   #
                    #
                    #
                    # Make directory structure for code execution
                    #
                    try:
                        make_work_area(workDir,code_src_dir,Risultati,Parameters)
                    except:
                        logging.error("Make work area: %s" % sys.exc_info()[0])
                        exitValueL = 256
                        pass
                    #
                    # prepare compile and execute
                    #
                    if not exitValueL == 256:
                        try:
                            logging.debug("Simulation %d => Run exoplates",Parameters['number'] )
                            exitValueL = exoclime(Parameters, workDir, code_src_dir, Risultati, emulate=True)
                        except:
                            logging.error("Exoclime %s" % sys.exc_info()[0])
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
                        logging.debug("Simulation %d => Archive Broken",Parameters['number'] )
                        try:
                            archive_broken_simulations(Parameters, workDir, Broken)
                        except:
                            logging.error("Unable to archive broken simulations: %s",sys.exc_info()[0])
                            pass
                    else:
                        logging.debug("Simulation %d => Archive data",Parameters['number'] )
                        try:
                            archive_exoplanet_data(Parameters, workDir,RisultatiMultipli,LogFiles,Risultati)
                        except:
                            logging.error("Unable to archive  simulations: %s",sys.exc_info()[0])
                            pass
                    logging.debug("START %d,0: End computation, ExitValue %d" % (rank,exitValueL))
                    # END                                               #
                    #####################################################
                    #sending back results:
                    inputdata.append(exitValueL)
                    comm.send(inputdata, dest=0, tag=tags.DONE)
                    logging.debug("START %d,0: Data sent to Master" % (rank))
                elif tag == tags.EXIT:
                    logging.info("EXIT and Break from simulation %d"  % (rank))
                    comm.send(None, dest=0, tag=tags.EXIT) # chiudi task mpi e esci
                    break
            else:
                local_simulation_index += 1
                exitValueL = 0
                inputdata = comm.recv(source=0, tag=MPI.ANY_TAG, status=status)
                tag = status.Get_tag()
                if tag == tags.START:
                    logging.debug("START on rank %d, local sim %d: begin computation" % (rank,local_simulation_index))
                    Parameters = make_input_parameters(inputdata,Parameters)
                    logging.debug("Simulation %d => Check exitValueL=%d" % (Parameters['number'] ,exitValueL))

                    ######################################################
                    # prepare, compile and execute code, store results   #
                    #
                    #
                    # Make directory structure for code execution
                    #
                    try:
                        make_work_area(workDir,code_src_dir,Risultati,Parameters)
                    except:
                        logging.error("Make work area: %s" % sys.exc_info()[0])
                        exitValueL = 256
                        pass
                    #
                    if not exitValueL == 256:
                        try:
                            logging.debug("Simulation %d => Run exoplates",Parameters['number'] )
                            exitValueL = exoclime(Parameters, workDir, code_src_dir, Risultati, emulate=True)
                        except:
                            logging.error("Exoclime %s" % sys.exc_info()[0])
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
                            logging.debug("Simulation %d => Archive Broken",Parameters['number'] )
                            archive_broken_simulations(Parameters, workDir, Broken)
                        except:
                            logging.error("unable to archive broken simulations: %s",sys.exc_info()[0])
                            pass
                    else:
                        try:
                            logging.debug("Simulation %d => Archive data",Parameters['number'] )
                            archive_exoplanet_data(Parameters, workDir,RisultatiMultipli,LogFiles,Risultati)
                        except:
                            logging.error("unable to archive  simulations: %s",sys.exc_info()[0])
                            pass
                    
                    # END                                               #
                    #####################################################
                    #sending back results:
                    inputdata.append(exitValueL)
                    comm.send(inputdata, dest=0, tag=tags.DONE)
                elif tag == tags.EXIT:
                    logging.info("EXIT and break from simulation %d"  % (rank))
                    comm.send(None, dest=0, tag=tags.EXIT)
                    break
        logging.info("Rank %d ends" % (rank))

    MPI.Finalize()

