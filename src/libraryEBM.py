#  LIBRARY OF PYTHON MODULES USED TO RUN EBM CODES
 

#########################################################################################

def modulationPAR(obliquity,R): 

   """
   compute parameters c0 and c1 of the modulation term of the diffusion coefficient

   the parameters area calculated as a function of obliquity and R

   the results are independent of the orbital period which appears as an internal parameter
   
   the calculations are based on Eqs. (A7), (A8), (A11), (A12)

   stricitly speaking THE RESULTS ARE EXACT FOR CIRCULAR ORBITS, but we do not
   expect significant variations for elliptical orbits
   """

   import numpy as np
   import math as m
   
   if R==1.: 
      c0=1.
      c1=0.
      return c0,c1
   
   print 'Starting calculation of C0 and C1'
   delta0=obliquity*np.pi/180.  #  OBLIQUITY
   Pdays=365                    #  ORBITAL PERIOD (days)
   P=86400.*Pdays               #  ORBITAL PERIOD (sec)
   omega=2.*np.pi/P             #  planet angular velocity

   degstep=0.5
   latitudes=np.arange(-90.,90.+degstep,degstep)

   muTOTAL=0.
   sumDAYS=0. 
   maxMU=0.
   minMU=1.

   for day in range(Pdays):
      t = 86400.*day  
      Ls=omega*t 
      delta=m.asin(-m.sin(delta0)*m.cos(Ls+np.pi/2)) 
   
      muAREA=0.
      sumWT=0. 
      for lat in latitudes:
      
         theta=lat*np.pi/180.
         cosH=-m.tan(theta)*m.tan(delta)
         if cosH > +1.: cosH=+1.
         if cosH < -1.: cosH=-1.
         H=m.acos(cosH)
         if H > 0:
            muDAY = m.sin(theta)*m.sin(delta)+m.cos(theta)*m.cos(delta)*m.sin(H)/H
         else:
            muDAY = 0.
      
         weight = m.cos(theta)
         muAREA=muAREA+weight*muDAY 
         sumWT = sumWT+ weight
      
         maxMU=max(maxMU,muDAY)
         minMU=min(minMU,muDAY)

      muAREA=muAREA/sumWT 
   
      muTOTAL=muTOTAL+muAREA
      sumDAYS=sumDAYS+1 
   

   muTOTAL=muTOTAL/sumDAYS  # final result of quantity defined in Eq. (A8)

   c1=maxMU/(R-1.)+muTOTAL   # OK CONTROLLATO
   c1=1./c1                 # Eq. (A11)

   c0=1.-c1*muTOTAL          # OK CONTROLLATO, Eq. (A7)

   return c0,c1


#########################################################################################


def LWTR(P):

#     input: total pressure (Pa)
#     output: liquid water temperature range (K)

#     data from CRC Handbook of Chemistry and Physics 
#     or from engineeringtoolbox.com
      
   P1=[5.,10.,1000.,10000.]        
   T1=[273.15,273.15,273.15,273.09]
      
   P2=[5.,10.,20.,30.,50.,100.,300.,1000.,2000.,3000.,5000.,7000.,8300.,10000.]
   T2=[273.15,273.15,288.,299.55,306.03,318.97,342.26,372.78,393.39,407.,425.,438.,445.,454.]
     
   P=P*1.e2 # convert to millibar
   
   for i in range(len(P1)-1):
      if P1[i] <= P <= P1[i+1]: 
         TWmin=T1[i]+(T1[i+1]-T1[i])*(P-P1[i])/(P1[i+1]-P1[i]) 
	 
   for i in range(len(P2)-1):
      if P2[i] <= P <= P2[i+1]: 
         TWmax=T2[i]+(T2[i+1]-T2[i])*(P-P2[i])/(P2[i+1]-P2[i]) 	 
      
   return TWmin,TWmax

#########################################################################################



def AtmPar(pDRY,p_CO2_P,p_CH4_P):

   """
   Calculate mean specific heat and molecular weight of the dry atmosphere
   resulting from variations of CO2 and CH4 partial pressure.
   The atmosphere is assumed to be of current terrestrial type (N2 and O2).
   All pressure units are in bar.

   """   
   from constantsEBM import *  
   
   # specific heat capacity of Earth air without greenhouse gases
   cp_Emgh=cp_E*pressEdry-cp_CO2*p_CO2_E-cp_CH4*p_CH4_E
   cp_Emgh=cp_Emgh/(pressEdry-p_CO2_E-p_CH4_E)
   
   cp_P=cp_Emgh*(pDRY-p_CO2_P-p_CH4_P)+cp_CO2*p_CO2_P+cp_CH4*p_CH4_P
   cp_P=cp_P/pDRY
   
   # molecular weight of Earth air without greenhouse gases
   molwtEmgh=molwtE*pressEdry-molwtCO2*p_CO2_E-molwtCH4*p_CH4_E
   molwtEmgh=molwtEmgh/(pressEdry-p_CO2_E-p_CH4_E)
   
   molwtP=molwtEmgh*(pDRY-p_CO2_P-p_CH4_P)+molwtCO2*p_CO2_P+molwtCH4*p_CH4_P
   molwtP=molwtP/pDRY  
   
   return cp_P,molwtP



def LinInterp(x,y,x0): # x MUST BE INCREASING
    "Linear interpolation at position x0 in arrays x,y"
    n=len(x)
    for i in range(n-1):
        if x[i] <= x0 <= x[i+1] :
            y0 = y[i] + (y[i+1]-y[i])*(x0-x[i])/(x[i+1]-x[i]) 
    if x0 < x[0]:
        y0 = 0.5*(y[0]+y[1])
        print 'LinInterp: EXTRAPOLATION AT START OF RANGE'
    if x0 > x[-1]: 
        y0 = 0.5*(y[-1]+y[-2])
        print 'LinInterp: EXTRAPOLATION AT END OF RANGE' 
    return y0
    
