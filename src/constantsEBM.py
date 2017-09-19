from libraryEBM import *

pi2=6.28318530717958647 # 2*pi (USED IN THE EBM CODE)

# EARTH DATA
pressE=1.01325e5   # [Pa] Earth surface pressure   1.01325 bar
pressEdry=1.0031e5  # [Pa] Earth surface pressure minus pH2O at rh=0.5 and T=<T_surf,earth> 
                   # pH2O = rh x p_sat,w(T_m,E)  
p_CO2_E=3.8e2  # [PPMV] Earth CO2 partial pressure   NOTE, in PPMV
p_CH4_E=1.8    # [PPMV] Earth CH4 partial pressure   NOTE, in PPMV
cp_E=1.005e3   # [J/(kg K)] Earth air specific heat capacity at T=0 C
molwtE=28.97   # Earth air mean molecular weight
 
gE=9.8         # [m/s2] Gravitational acceleration

# SPECIFIC HEAT CAPACITIES  FROM phys.py Pierrehumbert
cp_N2= 1.037e3   # N2 specific heat capacity at T=0 C p=1 bar [J/(kg K)]
cp_O2= 0.916e3   # O2 specific heat capacity at T=0 C p=1 bar [J/(kg K)]
cp_CO2=0.820e3   # CO2 specific heat capacity at T=0 C p=1 bar [J/(kg K)]
cp_CH4=2.195e3   # CH4 specific heat capacity at T=0 C p=1 bar [J/(kg K)]
cp_H2O=1.847e3   # H2O specific heat capacity at T=0 C p=1 bar [J/(kg K)]

# MEAN MOLECULAR WEIGHTS
molwtN2=28.    # Mean molecular weight of N2
molwtO2=32.    # Mean molecular weight of O2
molwtCO2=44.01 # Mean molecular weight of CO2
molwtCH4=16.04 # Mean molecular weight of CH4
molwtH2O=18.01 # Mean molecular weight of H2O

