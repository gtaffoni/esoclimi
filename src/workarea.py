#!/usr/bin/env python

#  work_env.py
#  Esoclimi
#
#  Created by giuliano taffoni on 04/04/17.
#  Copyright 2017 Guliano Taffoni. All rights reserved.



def copyall (src, dest):
    '''
        copyall(src,dest)
        
        copy all files from source to destination one by one
        
    '''
    import os
    import shutil
    src_files = os.listdir(src)
    for file_name in src_files:
        full_file_name = os.path.join(src, file_name)
        if (os.path.isfile(full_file_name)):
            shutil.copy(full_file_name, dest)
    return


def collective_move (src, dest):
    '''
        collective_move (src, dest)
        
        move all files from src directory  to dest directory
        one by one
    '''
    import os
    import shutil
    src_files = os.listdir(src)
    for file_name in src_files:
        full_file_name = os.path.join(src, file_name)
        if (os.path.isfile(full_file_name)):
            shutil.move(full_file_name, dest)
    return





def archive_results(str, src, planet, Risultati):
    '''
        archive_results(str,src,Risultati)
        
        archive results from run directories
        
        str   =  multiple results directory
        src  =  directory where PRESSUREScurr.py and parEBM.h are located
        Risultati = results dir
        
    '''
    import os
    import shutil
    shutil.copy(src+"/parEBM.h",Risultati)
    shutil.copy(src+"/"+planet+".h",Risultati)
    os.mkdir(str)
    collective_move(Risultati,str)
    return



def archive_logs(str,log):
    '''
        archive_logs(src,logs)
        
        archive logs collected during execution
        
    '''
    import shutil
    shutil.move(log,str)


def CleanAllPartialResults(localDir):
    '''
        CleanAllPartialResults(localDir)
        
        remove all data from temporary running directories before
        a new run is executed
    '''
    import shutil
    shutil.rmtree(localDir)
    return
