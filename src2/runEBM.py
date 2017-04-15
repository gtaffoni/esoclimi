#!/usr/bin/python
#  runEBM.py
#  Esoclimi
#
#  Created by Giuseppe Murante.
#  Modified by Guliano Taffoni.

"""
     
    setupEBM(planet,NUMBER=0,VERSION="Unknwn",SIMTYPE="---",_dir):
    
    setup the code for compilation
    
    _dir == location of the source files to setup and compile
    
"""


def setupEBM(planet,_p,_ecc,_obl,_dist, _gg, _p_CO2_P, _fo_const, _TOAalbfile, _OLRfile, _dir,NUMBER=0,VERSION="Unknwn",SIMTYPE="---"):
# REMEMBER TO ADD CALCULUS OF PLANET MASS AND CHANGE BELOW


     import sys
     from libraryEBM import modulationPAR
     from constantsEBM import *
     import logging
     import shutil

     logging.info("reading %s parameters", planet)
     ### assign parameter values from input better in the future
     ff            = _p
     eccP          = _ecc
     obliq         = _obl
     smaP          = _dist
     gg            = _gg
     p_CO2_P       = _p_CO2_P
     fo_const      = _fo_const
     tabTOAalbfile = _TOAalbfile
     tabOLRfile    = _OLRfile

     ####
     startpar_file = _dir+"startpar.h"
     parEBM_file   = _dir+"/parEBM.h"
     planet_file   = _dir+planet+".h"
     planet_orig   = _dir+"/planet.h" 


     #EARTH VALUE! warning it's in PASCAL - this is why I divide by 10, in the calling cycle is in ppvm
     OLRmodel='ccm3tab0'   # ['ccm3tab0','CCMcal00', 'CCMcalCF'] 
                           # ccm3tab0: OLR calculated with CCM3 taken at face value (no cloud forcing)
                           # CCMcal00: OLR calibrated with CCM, WITHOUT correction factors
                           # CCMcalCF: OLR calibrated with CCM, WITH correction factors 

     p_CO2_P /= 10.
     p_CH4_P=p_CH4_E       # planet CH4 partial pressure [bar]
     press_E  = 101325.0 + p_CO2_P 
     R=2.2        # ratio max(modulation term)/min(modulation term) of the diffusion coefficient
     pmass    = 1.0 #planet mass in Earth masses  

     # calculate parameters c0 and c1 of modulation term, Eqs. (A11), (A12) 
     C0par,C1par=modulationPAR(obliq,R)

     # calculate cp_P and molwtP given pressure and chemical composition
     pressP= pressE*ff
     cp_P,molwtP=AtmPar(pressP,p_CO2_P,p_CH4_P)
    
     
     # parsing "planet.h", producing planet+".h" (e.g., EARTH.h)  
     tf=open(planet_orig,'r') # input template file
     pf=open(planet_file,'w')        # output parameter file
     
     # updates the EBM parameters file (e.g. parEBM.h) with current values of the planets parameter file (e.g. EARTH.py)
     # note: here I will only change the parameter that we are exploring!
     for tl in tf:
     
        pl = tl
        
        if 'parameter(planet' in tl and 'character' not in tl:
           pl="        parameter(planet='%s')\n" % planet

        if 'parameter(Pmass' in tl and 'real' not in tl:
           pl="        parameter(Pmass=%f)\n" % pmass   #WARNING this is SET not calculated
            
        if 'parameter(smaP' in tl and 'real' not in tl:
           pl='        parameter(smaP=%7.5f) ' % smaP      
           pl=pl+'! semi-major axis of planet orbit [AU]\n'
     
        if 'parameter(eccP' in tl and 'real' not in tl:
           pl='        parameter(eccP=%7.5f) ' % eccP      
           pl=pl+'! eccentricity of planet orbit\n'
                      
        if 'parameter(obliq' in tl and 'real' not in tl:
             pl='        parameter(obliq=%8.4f)   ' % obliq
             pl=pl+'! [deg] planet axis inclination\n'

        if 'parameter(gg' in tl and 'integer' not in tl:
           pl='        parameter(gg=%i)   ' % gg
           pl=pl+'! planet geography \n'
           
        if 'parameter(fo_const' in tl and 'real' not in tl:
           pl='        parameter(fo_const=%5.2f)   ' % fo_const
           pl=pl+'! constant ocean fraction (if gg=0) \n'
           
        if 'parameter(TOAalbFile' in tl and 'character' not in tl:
           pl="        parameter(TOAalbFile='%s')   " % tabTOAalbfile
           pl=pl+'! file with TOA albedo matrix \n'     
           
        if 'parameter(p_CO2_P' in tl and 'real' not in tl:
           pl='        parameter(p_CO2_P=%9.2e)       ' % p_CO2_P
           pl=pl+'! [Pa] CO2 partial pressure \n' 
           
        if 'parameter(pressP' in tl and 'real' not in tl:
           pl='        parameter(pressP=%11.5e)   ' % pressP
           pl=pl+'! [Pa] planet total pressure\n'  
     
        # these are calculated and must be modified
        if 'parameter(C0par' in tl and 'real' not in tl:
             pl='        parameter(C0par=%8.5f)         ' % C0par
             pl=pl+'! diffusion parameter C0\n'
             
        if 'parameter(C1par' in tl and 'real' not in tl:
             pl='        parameter(C1par=%8.5f)         ' % C1par
             pl=pl+'! diffusion parameter C1\n'

        if 'parameter(molwtP' in tl and 'real 'not in tl:
           pl='        parameter(molwtP=%5.2f)   ' % molwtP
           pl=pl+'! Planet atmosphere mean molecular weight\n'
     
        if 'parameter(cp_P' in tl and 'real 'not in tl:
           pl='        parameter(cp_P=%6.3f)   ' % cp_P
           pl=pl+'! Planet atmosphere specific heat capacity at T=0 C\n' 
             
        # addendum GM, 8/2/2017: data for the generation of the fit file   
        if 'parameter(VERSION' in tl and 'character' not in tl:
           pl='        parameter(VERSION="%s")         ' % VERSION
           pl=pl+'!code version\n' 
        if 'parameter(SIMTYPE' in tl and 'character' not in tl:
           pl='        parameter(SIMTYPE="%s")         ' % SIMTYPE
           pl=pl+'!simulation type\n' 
        if 'parameter(NUMBER' in tl and 'integer' not in tl:
           pl='        parameter(NUMBER=%s)         ' % int(NUMBER)
           pl=pl+'!simulation number for current date\n' 
           
     
        pf.write(pl) 
     
     
     ###############   OLR taken from ccm3 calculations  
     
     if OLRmodel=='ccm3tab0':
     
        # EXTRACT AN ARRAY OF OLR values VERSUS TEMPERATURE
        # APPROPRIATE FOR THE ADOPTED PRESSURE
        
        import numpy as np 
        tab=open(_dir+"/"+tabOLRfile,'r')  # CHECK HERE!!!!
        
        # read pressure array from the header
        parr=[]
        tabHeader=tab.readline()
        e=tabHeader.split()
        for k in range(1,len(e)): 
           parr.append( float(e[k]) )
        pOLR=np.array(parr) 
        Np=len(parr) 
        
        # read temperature array from the first column
        Tarr=[]
        NT=0
        for line in tab:
           e=line.split()
           Tarr.append( float(e[0]) )
           NT += 1
        TOLR=np.array(Tarr)	 
     
        # read matrix of OLR values
        OLRmat=np.zeros((NT,Np),float) 
        tab=open(_dir+"/"+tabOLRfile,'r')
        iT=0
        for line in tab:
           if '#' not in line:
        	 e=line.split() 
        	 for ip in range(Np): 
        	    OLRmat[iT,ip]=float( e[ip+1] ) 
        	 iT += 1 
        	
        # search index of pressure columns to be interpolated
        for ip in range(Np-1):
           if pOLR[ip] <= pressP/1.e5 <= pOLR[ip+1]:
        	 ip0=ip 
        
        x0=pressP/1.e5
        xi=pOLR[ip0]
        xp=pOLR[ip0+1]
        
        vOLR=np.zeros((NT),float)
        for iT in range(NT):  
        
           yi=OLRmat[iT,ip0]
           yp=OLRmat[iT,ip0+1]  
        
           vOLR[iT] = yi+(yp-yi)*(x0-xi)/(xp-xi) 
        
        # store data in parameters file
        pf.write('\n')
        pf.write('c       Array of OLR values\n')
        pf.write('	     integer NOLR\n')
        pl='        parameter(NOLR=%i)\n' % NT
        pf.write(pl)
        pl='        real*8 TOLR(NOLR),vOLR(NOLR)\n'
        pf.write(pl)
        
        NDATA=6 # number of data stored in one line of the fortran statement 'data'
        
        pl='        data TOLR/' 
        exit=0
        iCF=0
        while (iCF < NT-1):
           k=0
           while (k < NDATA): 
        	 if iCF == NT-1: 
        	    pl=pl+'%5.1f/ \n'%TOLR[iCF]
        	    exit=1
        	 if iCF < NT-1:
        	    pl=pl+'%5.1f,'%TOLR[iCF]  
        	 k   += 1  
        	 iCF += 1
           if exit==0:
        	 pl=pl+'\n'
        	 pl=pl+'     >     ' 
        pf.write(pl)
        
        pl='        data vOLR/' 
        exit=0
        iCF=0
        while (iCF < NT-1):
           k=0
           while (k < NDATA): 
        	 if iCF == NT-1: 
        	    pl=pl+'%5.1f/ \n'%vOLR[iCF]
        	    exit=1
        	 if iCF < NT-1:
        	    pl=pl+'%5.1f,'%vOLR[iCF]  
        	 k   += 1  
        	 iCF += 1
           if exit==0:
        	 pl=pl+'\n'
        	 pl=pl+'     >     ' 
        pf.write(pl)
        
     tf.close()
     pf.close()
        
     #now, adds the include at the end of startpar.h
     shutil.copy(startpar_file,parEBM_file)
     with open(parEBM_file,"a") as f:
          f.write('\n')
          f.write('        include \'%s\'         \n' % (planet+'.h'))
          f.write('\n')


     return


def compileEBM(runDir,logfile):
    import logging
    #from posix import system
    import subprocess
    import os
    logging.info("Compile codes")
    origin = os.getcwd()
    # move to the proper directory
    os.chdir(runDir)
    logging.info("")
    p = subprocess.call("make", stdout=logfile,stderr=subprocess.STDOUT,shell=True)
    # p is the return code of Make so we can make some check
    #system('make')
    os.chdir(origin)


def runEBM(runDir,logfile):
    import logging
    import subprocess
    import os
    origin = os.getcwd()
    # move to the proper directory
    logging.info("Run program")
    os.chdir(runDir)
    p = subprocess.call("./codeEBM.x", stdout=logfile,stderr=subprocess.STDOUT,shell=True)
    os.chdir(origin)
