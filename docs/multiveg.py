#!/usr/bin/python
"""
This is an example of python script for runnung many times the 
code, varying one parameter written in a .h file
Here is the vegetation albedo
"""
import sys
import os
from posix import system

#values to cycle
VegetationAlbedo=[0.15, 0.13, 0.1] 

#simtype, version, simulation number
simtype="VegAlbedoFB"
version="1.1.02"
number=1

# multiple output directory preparation
system("mkdir -p RisultatiMultipli")
system("mkdir -p Database")
system("cp -r CCM_RH60 Src")


for va in VegetationAlbedo:

    # preparing Src (usually done by run.bash)
    system("cp ModulesDef/* Src")
    system("cp Std/* Src")
    if simtype=="VegPassive":
        system("cp VegPassive/* Src")
        system("cp= VegPassive/Modules/* Src")
    elif simtype=="VegAlbedoFB":
        system("cp VegPassive/* Src")
        system("cp VegPassive/Modules/* Src")
        system("cp VegAlbedoFB/* Src")
        system("cp VegAlbedoFB/Modules/* Src")
                
    # varying pressure on an EARTH-LIKE planet, using EARTH template
    system("cp Planets/EARTH.py Src/VEGALBEDOES.py")
    system("cp Planets/fo_earth_DMAP.dat Src")

    #set working dir
    os.chdir("Src") # (a call to system would open a shell, chd and return )
    system("pwd") 
    

    #changing vegetation albedo in  vegetation.h


    #configuring planet parameter file
    #PRESSURES.py is the template (here, EARTH.py)
    #PRESSUREScurr.py is what runEBM.py uses
    # note that here you can change any parameter in EARTH.py in the same way
    # note that you could also decide to change startpar.h or other files
    # (see vegetati=on example multipleveg.py) 
    system("cp vegetation.h vegetation_template.h")
    fileout=open('vegetation.h',"w")
    with open('vegetation_template.h', 'r') as f:
        for tl in f:
            pl=tl
            if 'parameter(aveg=' in tl  and 'character' not in tl:
                pl="        parameter(aveg=%4.2f)\n"%va
            fileout.write(pl)
    fileout.close()

    system("mkdir Risultati")
    #executing runEBM.py - this compiles and runs the code

    str="python runEBM.py VEGALBEDOES.py %d %s %s"%(number,version,simtype)
    print str
    system(str)


    #now we have results. Making .fits file (REQUIRES PYFIT)
    os.chdir("Risultati")
    system("python ../../fits_map_temperature.py")
    system("cp ESTM*fits ../../Database")

    #archiving Risultati
    os.chdir("..")
    system("cp vegetation.h Risultati")
    system("cp parEBM.h Risultati");
    str="mkdir ../RisultatiMultipli/RisultatiVA%4.2f"%va
    system(str)
    str="mv Risultati/* ../RisultatiMultipli/RisultatiVA%4.2f"%va
    system(str)

    #next simulation
    number=number+1

    #back to main dir
    os.chdir("..")
