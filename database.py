#!/usr/bin/python
"""
This is an example of python script for runnung many times the 
code, varying one parameter existing in one module .h file 
Here is the vegetation albedo in vegetation.h

GT 29/3/17 some small improvements
"""

import sys
import os
import shutil
from posix import system


def copyall (src, dest):
    src_files = os.listdir(src)
    for file_name in src_files:
        full_file_name = os.path.join(src, file_name)
        if (os.path.isfile(full_file_name)):
            shutil.copy(full_file_name, dest)
    return

def collective_move (src, dest):
    src_files = os.listdir(src)
    for file_name in src_files:
        full_file_name = os.path.join(src, file_name)
        if (os.path.isfile(full_file_name)):
            shutil.move(full_file_name, dest)
    return


#values to cycle
Pressures=[0.01, 0.1, 0.5, 1.0, 3.0, 5.0]  #times the Earth value
#nota, 0.001 e 0.005 non funzionavano
#Pressures=[1.0, 3.0, 5.0]  #times the Earth value
Radii= [0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5]
Obliquities = [0., 15., 23.43929, 30., 45.]
Eccentricities = [ 0.0, 0.01671022, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]


#simtype, version, simulation number
simtype="Std"
version="1.1.03"
number=1

# multiple output directory preparation
# here we need to catch exeptions (necessary on parallel runs)
os.makedirs("RisultatiMultipli")
os.makedirs("Database")
os.makedirs("LogFiles")

shutil.copytree("CCM_RH60", "Src")

for ecc in Eccentricities:
    for obl in Obliquities:
        for dist in Radii:            
            for p in Pressures:
                    print p,ecc,obl,dist

                    copyall("ModulesDef", "Src")
                    copyall("Std","Src")
                    if simtype=="VegPassive":
                        copyall("VegPassive","Src")
                        copyall("VegPassive/Modules/","Src")
                    elif simtype=="VegAlbedoFB":
                        copyall("VegAlbedoFB", "Src")
                        copyall("VegAlbedoFB/Modules/","Src")
                
    # varying pressure on an EARTH-LIKE planet, using EARTH template
                    shutil.copy("Planets/EARTH.py","Src/PRESSURES.py")
                    shutil.copy("Planets/fo_earth_DMAP.dat","Src/fo_earth_DMAP.dat")
    #set working dir
                    os.chdir("Src") # (a call to system would open a shell, chd and return )
                    system("pwd") # GT a che serve?
    

    #configuring planet parameter file
    #PRESSURES.py is the template (here, EARTH.py)
    #PRESSUREScurr.py is what runEBM.py uses
                    fileout=open('PRESSUREScurr.py',"w")
                    with open('PRESSURES.py', 'r') as f:
                        for tl in f:
                            pl = tl
                            if 'ff=1.0' in tl:
                                pl="ff=%f\n"%p
                            if 'eccP=' in tl:
                                pl='eccP=%f\n'%ecc
                            if 'obliq=' in tl:
                                pl='obliq=%f\n'%obl
                            if 'smaP=' in tl:
                                pl='smaP=%f\n'%dist
                            fileout.write(pl)
                    fileout.close()                    

                    os.makedirs("Risultati")
    #executing runEBM.py - this compiles and runs the code

# GT Parliamone e' proprio necessario? non converrebbe fare una funzione?
                    str="python runEBM.py PRESSUREScurr.py %d %s %s > log "%(number,version,simtype)
                    print str
                    system(str)

    #now we have results. Making .fits file (REQUIRES PYFIT)
                    os.chdir("Risultati")
                    system("python ../../fits_map_temperature.py")

                    for file in os.listdir("/mydir"):
                        if file.endswith(".fits") and file.find(ESTM):
                            shutil.copy(file,"../../Database")

                    os.chdir("..")
    #archiving log file
                    system("pwd")
                    system("ls PR*")
                    str="../LogFiles/log_Press%5.3f_Ecc%4.2f_Dist%3.1f_Obl%5.3f"%(p,ecc,dist,obl)
                    shutil.move("log",str)

    #archiving Risultati
                    shutil.copy("PRESSUREScurr.py","Risultati")
                    shutil.copy("parEBM.h","Risultati")
                    str="../RisultatiMultipli/Risultati_Press%5.3f_Ecc%4.2f_Dist%3.1f_Obl%5.3f"%(p,ecc,dist,obl)
                    os.mkdir(str)
                    str="../RisultatiMultipli/Risultati_Press%5.3f_Ecc%4.2f_Dist%3.1f_Obl%5.3f"%(p,ecc,dist,obl)
                    collective_move("Risultati",str)
    #next simulation
                    number=number+1

    #back to main dir
                    os.chdir("..")
                    print "\n\n\n"
