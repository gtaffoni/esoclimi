from libraryEBM import *
from constantsEBM import *

N= 54              # latitude grid points 
Tstart=275.        # initial temperature

## PARAMETERS from JENKINS ET AL. 2015, ApJ

Mstar=1.037             # stellar mass [Msun]
#Mstar=1.0             # stellar mass [Msun]
#LumStar=1.21 # stellar luminosity [Lsun] from Teff and R*, tuned with Jenkins estimates of Seff 
LumStar=1.1 # stellar luminosity [Lsun] from Teff and R*, tuned with Jenkins estimates of Seff 
LumEvol=1.00          # evolving stellar luminosity [fraction of present luminosity]

smaP=1.048 # semi-major axis of planet orbit;  tuned with Jenkins measurements of orbital period 
eccP=0.0   # eccentricity from (e cos omega) and (e sin omega), with omega=293 of planet orbit  [earthECC=0.01671022]
omegaPERI=0. # fixed, should be 293.  # argument of pericenter (angle from line of nodes to pericenter in the orbital plane) 


Ls0_string='pi2/4.'   # simulation starts at spring equinox
#Ls0_string='0.'      # simulation starts at winter solstice 

#obliq=22.5 #23.43929        # obliquity of planet axis of rotation [deg] 
obliq=0.0       # obliquity of planet axis of rotation [deg] 

# PLANET "SOLAR" DAY (FROM NOON TO NOON) 
# estimated from equation: Prot = 1 / { 1/(P[rot,sid]) - 1/(P[orbit,sid])
Prot=1.0            # planet rotation period [days]

#Rjup=10.9733 # Jupiter mean radius in units of Earth mean radii
Rplan=1.63           # planet radius [R(earth)]

gg=0          # ocean fractions: 0: constant fo
	      # 1: present Earth WK97; 2: equatorial cont.WK97; 3: polar cont.WK97
	      # 4: present Earth, file 'fo_earth_DMAP.dat' 

fo_const=0.0          # constant fraction of ocean (only used when gg=0)

OLRmodel='ccm3tab0'   # ['ccm3tab0','CCMcal00', 'CCMcalCF'] 
                      # ccm3tab0: OLR calculated with CCM3 taken at face value (no cloud forcing)
                      # CCMcal00: OLR calibrated with CCM, WITHOUT correction factors
		      # CCMcalCF: OLR calibrated with CCM, WITH correction factors 

if OLRmodel=='ccm3tab0':  
   tabOLRfile='CCM_RH60/OLR_g15_rh60_co380_ch1.8.txt' 
tabTOAalbfile='CCM_RH60/ALB_g15_rh60_co380_ch1.8.txt'
 
ff=3.0
pressP=ff*pressEdry   # total dry pressure of the planet (N2/O2 atmosphere, greenhouse gases included, water excluded)
p_CO2_P=p_CO2_E       # planet CO2 partial pressure [bar]
p_CH4_P=p_CH4_E       # planet CH4 partial pressure [bar]

cp_P,molwtP=AtmPar(pressP,p_CO2_P,p_CH4_P)

gP = 1.5*gE # corresponds to rh0=5.075 g/cm3 and M(K452b)=4.0 M(earth)
# gP=gE # corresponds to rho=3.38 g/cm3 and M(K452b)=2.7 M(earth)

albedoFile='none'

#asl=0.18             # surface albedo of lands (=0.2 in WK97)
asl=0.35             # surface albedo of lands (=0.2 in WK97)
asil=0.70            # surface albedo of ice on lands (=0.85 Pierrehumbert)
asio=asil            # surface albedo of ice on ocean (=0.50 Pierrehumbert) 

#fcw=0.70              # cloud coverage on water (mean from Sanroma & Palle 2011=0.70)
fcw=0.40              # cloud coverage on water (mean from Sanroma & Palle 2011=0.70)
#fcl=0.60              # cloud coverage on land (mean from Sanroma & Palle 2011=0.50)
fcl=0.30              # cloud coverage on land (mean from Sanroma & Palle 2011=0.50)
fci=fcl               # cloud coverage on ice (mean da Sanroma & Palle 2011=0.50) 
 
R=2.2        # ratio max(modulation term)/min(modulation term) of the diffusion coefficient
D0par=0.66   # [W/m2]
    
