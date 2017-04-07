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


def setupEBM(planet,_p,_ecc,_obl,_dist,_dir,NUMBER=0,VERSION="Unknwn",SIMTYPE="---"):
     import sys
     from libraryEBM import modulationPAR
     from constantsEBM import *
     import logging
     logging.info("reading %s parameters", planet)
     ### assign parameter values from inpunt better in the future
     ff    = _p
     eccP  = _ecc
     obliq = _obl
     smaP  = _dist
     ####
     startpar_file = _dir+"/startpar.h"
     parEBM_file   = _dir+"/parEBM.h"
     # calculate parameters c0 and c1 of modulation term, Eqs. (A11), (A12) 
     C0par,C1par=modulationPAR(obliq,R)
     
     tf=open(startpar_file,'r') # input template file
     pf=open(parEBM_file,'w')        # output parameter file
     
     # updates the EBM parameters file (e.g. parEBM.h) with current values of the planets parameter file (e.g. EARTH.py)
     for tl in tf:
     
        pl = tl
        
        if 'parameter(planet' in tl and 'character' not in tl:
           pl="        parameter(planet='%s')\n" % planet
        
        if 'parameter(N=' in tl and 'integer' not in tl:
           pl='	parameter(N=%i)' % N
           pl=pl+'! latitude grid points\n'
           
        if 'parameter(Tstart' in tl and 'real' not in tl:
           pl='        parameter(Tstart=%5.0f) ' % Tstart      
           pl=pl+'! initial temperature [HOT START: 350 K]\n'
     
        if 'parameter(Mstar' in tl and 'real' not in tl:
           pl='        parameter(Mstar=%7.3f*Msun)   ' % Mstar
           pl=pl+'! stellar mass [Msun]\n'
           
        if 'parameter(LumStar' in tl and 'real' not in tl:
           pl='        parameter(LumStar=%8.4f*LumSun)   ' % LumStar
           pl=pl+'! stellar luminosity [Lsun]\n' 
           
        if 'parameter(LumEvol' in tl and 'real' not in tl:
           pl='        parameter(LumEvol=%7.3f)   ' % LumEvol
           pl=pl+'! evolving stellar luminosity\n'   
     
        if 'parameter(smaP' in tl and 'real' not in tl:
           pl='        parameter(smaP=%7.5f) ' % smaP      
           pl=pl+'! semi-major axis of planet orbit [AU]\n'
     
        if 'parameter(eccP' in tl and 'real' not in tl:
           pl='        parameter(eccP=%7.5f) ' % eccP      
           pl=pl+'! eccentricity of planet orbit\n'
           
        if 'parameter(omegaPERI' in tl and 'real' not in tl:
           pl='        parameter(omegaPERI=%7.5f) ' % omegaPERI      
           pl=pl+'! argument of pericenter (from ascending node in the orbital plane)\n'   
        
        if 'parameter(Ls0' in tl and 'real' not in tl:
           pl='        parameter(Ls0=%s)   ' % Ls0_string
           if 'pi' in Ls0_string: msg='simulation starts at spring equinox'
           if '0.' in Ls0_string: msg='simulation starts at winter solstice' 
           pl=pl+'! '+msg+'\n' 
           
        if 'parameter(obliq' in tl and 'real' not in tl:
             pl='        parameter(obliq=%8.4f)   ' % obliq
             pl=pl+'! [deg] planet axis inclination\n'
     
        if 'parameter(Prot' in tl and 'real' not in tl:
           pl='        parameter(Prot=%7.3f)   ' % Prot
           pl=pl+'! planet rotation period\n'
           
        if 'parameter(Rplanet' in tl and 'real' not in tl:
           pl='        parameter(Rplanet=%7.4f*Rearth)   ' % Rplan
           pl=pl+'! planet radius [m] \n'
           
        if 'parameter(gP' in tl:
           pl='        parameter(gP=%6.2f)   ' % gP
           pl=pl+'! planet gravitational acceleration\n'
           
        if 'parameter(gg' in tl and 'integer' not in tl:
           pl='        parameter(gg=%i)   ' % gg
           pl=pl+'! planet geography \n'
           
        if 'parameter(fo_const' in tl and 'real' not in tl:
           pl='        parameter(fo_const=%5.2f)   ' % fo_const
           pl=pl+'! constant ocean fraction (if gg=0) \n'
           
        if 'parameter(albedoFile' in tl and 'character' not in tl:
           pl="        parameter(albedoFile='%s')   " % albedoFile
           pl=pl+'! file with albedo data \n'    
     
        if 'parameter(TOAalbFile' in tl and 'character' not in tl:
           pl="        parameter(TOAalbFile='%s')   " % tabTOAalbfile
           pl=pl+'! file with TOA albedo matrix \n'     
           
        if 'parameter(OLRmodel' in tl and 'character' not in tl:
           pl="	parameter(OLRmodel='%s')\n" % OLRmodel 
        
        if 'parameter(p_CO2_P' in tl and 'real' not in tl:
           pl='        parameter(p_CO2_P=%9.2e)       ' % p_CO2_P
           pl=pl+'! [Pa] CO2 partial pressure \n' 
           
     #   if 'parameter(tauIR0' in tl and 'real' not in tl and OLRmodel=='SP':
     #      pl='        parameter(tauIR0=%7.4f)    ' %  tauIR0
     #      pl=pl+'! Spiegel et al value=0.79\n'  
     #      
     #   if 'parameter(SPindex' in tl and 'real' not in tl and OLRmodel=='SP':
     #      pl='        parameter(SPindex=%7.4f)    ' %  SPindex
     #      pl=pl+'! exponent in Spiegel tauIR: (T/273)^SPindex\n'
          
        if 'parameter(pressP' in tl and 'real' not in tl:
           pl='        parameter(pressP=%11.5e)   ' % pressP
           pl=pl+'! [Pa] planet total pressure\n'  
     
        if 'parameter(molwtP' in tl and 'real 'not in tl:
           pl='        parameter(molwtP=%5.2f)   ' % molwtP
           pl=pl+'! Planet atmosphere mean molecular weight\n'
     
        if 'parameter(cp_P' in tl and 'real 'not in tl:
           pl='        parameter(cp_P=%6.3f)   ' % cp_P
           pl=pl+'! Planet atmosphere specific heat capacity at T=0 C\n' 
        
        if 'parameter(D0par' in tl and 'real' not in tl:
           pl='        parameter(D0par=%5.3f)         ' % D0par
           pl=pl+'! [W/m2] diffusion parameter D0\n'
           
        if 'parameter(Rpar' in tl and 'real' not in tl:
           pl='        parameter(Rpar=%5.1f)         ' % R
           pl=pl+'! diffusion parameter R\n'      
           
        if 'parameter(C0par' in tl and 'real' not in tl:
           pl='        parameter(C0par=%8.5f)         ' % C0par
           pl=pl+'! diffusion parameter C0\n'
           
        if 'parameter(C1par' in tl and 'real' not in tl:
           pl='        parameter(C1par=%8.5f)         ' % C1par
           pl=pl+'! diffusion parameter C1\n'
           
        if 'parameter(asl' in tl and 'real' not in tl:
           pl='        parameter(asl=%5.2f)         ' % asl
           pl=pl+'! surface albedo of lands (=0.2 in WK97)\n' 
           
        if 'parameter(asil' in tl and 'real' not in tl:
           pl='        parameter(asil=%5.2f)         ' % asil
           pl=pl+'! surface albedo of ice on lands (=0.85 Pierrehumbert)\n'
        
        if 'parameter(asio' in tl and 'real' not in tl:
           pl='        parameter(asio=%5.2f)         ' % asio
           pl=pl+'! surface albedo of ice on ocean (=0.50 Pierrehumbert)\n'
     
        if 'parameter(cloudIR' in tl and 'real' not in tl and OLRmodel=='WK':
           pl='        parameter(cloudIR=%4.1f)         ' % cloudIR
           pl=pl+'! cloud IR absorption in the outgo [=14 W in WK97]\n'
     
        if 'parameter(fcw' in tl and 'real' not in tl:
           pl='        parameter(fcw=%5.2f)         ' % fcw
           pl=pl+'! cloud coverage on water (mean from Sanroma & Palle 2011)\n'
      
        if 'parameter(fcl' in tl and 'real' not in tl:
           pl='        parameter(fcl=%5.2f)         ' % fcl
           pl=pl+'! cloud coverage on land (mean from Sanroma & Palle 2011)\n'
     
        if 'parameter(fci' in tl and 'real' not in tl:
           pl='        parameter(fci=%5.2f)         ' % fci
           pl=pl+'! cloud coverage on ice (mean from Sanroma & Palle 2011)\n' 
     
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
        
        #print TOLR
        #print vOLR
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

