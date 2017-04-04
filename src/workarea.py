#!/usr/bin/env python

#  work_env.py
#  Esoclimi
#
#  Created by giuliano taffoni on 04/04/17.
#  Copyright © 2017 Guliano Taffoni. All rights reserved.


'''
    copyall(src,dest) 
    
        copy all files from source to destination one by one
    
'''
def copyall (src, dest):
    src_files = os.listdir(src)
    for file_name in src_files:
        full_file_name = os.path.join(src, file_name)
        if (os.path.isfile(full_file_name)):
            shutil.copy(full_file_name, dest)
    return

'''
   collective_move (src, dest)
   
        move all files from src directory  to dest directory
        one by one
'''

def collective_move (src, dest):
    src_files = os.listdir(src)
    for file_name in src_files:
        full_file_name = os.path.join(src, file_name)
        if (os.path.isfile(full_file_name)):
            shutil.move(full_file_name, dest)
    return


'''
    make_work_area (_p,_ecc,_obl,_dist)
    
        create the working area dircetory and file structure coping from
        template directories
        
        p    = pressure
        ecc  = eccentricity
        obl  =
        dist = distance
    
'''
def make_work_area (_p,_ecc,_obl,_dist):
    # multiple output directory preparation
    # here we need to catch exeptions (necessary on parallel runs)
    os.makedirs(RisultatiMultipli)
    os.makedirs(Database)
    os.makedirs(LogFiles)
    os.makedirs(Risultati)
    
    shutil.copytree(template_dir+"/CCM_RH60", Src) # GT why here?

    copyall("ModulesDef", Src)

    copyall("Std",Src)
    
    if simtype=="VegPassive":
        copyall(template_dir+"/VegPassive",Src)
        copyall(template_dir+"/VegPassive/Modules/",Src)
    elif simtype=="VegAlbedoFB":
        copyall(template_dir+"/VegAlbedoFB", Src)
        copyall(template_dir+"/VegAlbedoFB/Modules/",Src)
    # varying pressure on an EARTH-LIKE planet, using EARTH template
    shutil.copy("Planets/EARTH.py",Src+"/PRESSURES.py")
    shutil.copy("Planets/fo_earth_DMAP.dat",Src+"/fo_earth_DMAP.dat")
##
#   configuring planet parameter file
#
#   PRESSURES.py is the template (here, EARTH.py)   GT: maybe better to call it PRESSURES.templ ?
#   PRESSUREScurr.py is what runEBM.py uses
###
        fileout=open(Src+'/PRESSUREScurr.py',"w")
        with open(Src+'/PRESSURES.py', 'r') as f:
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


'''
    archive_results(p,ecc,dist,obl)
    
        archive results from run directory

        p    = pressure
        ecc  = eccentricity
        obl  =
        dist = distance
        


'''
def archive_results(p,ecc,dist,obl):
    shutil.copy("PRESSUREScurr.py",Risultati)
    shutil.copy("parEBM.h",Risultati)
    str="%s/Risultati_Press%5.3f_Ecc%4.2f_Dist%3.1f_Obl%5.3f"%(RisultatiMultipli,p,ecc,dist,obl)
    os.mkdir(str)
    str="%s/Risultati_Press%5.3f_Ecc%4.2f_Dist%3.1f_Obl%5.3f"%(RisultatiMultipli,p,ecc,dist,obl)
    collective_move(Risultati,str)
    return

'''
    archive_logs(p,ecc,dist,obl)
    
        archive logs collected during execution

    p    = pressure
    ecc  = eccentricity
    obl  =
    dist = distance
'''
def archive_logs(p,ecc,dist,obl):
    str="%s/log_Press%5.3f_Ecc%4.2f_Dist%3.1f_Obl%5.3f"%(LogFiles,p,ecc,dist,obl)
    shutil.move(log,str)

