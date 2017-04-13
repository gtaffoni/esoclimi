*
*       Constants and parameters for codeEBM
* 
        character*16 planet
        character*50 TOAalbFile
        integer N, Ns
        real*8 latbc1, latbc2
	integer gg ! type of geography  
        real*8 obliq, OBLIQUITY, Pmass
        real*8 Prot, OmegaP
        real*8 Rplanet, gP
        real*8 fo_const
        character*16 OLRmodel
        real*8 p_CO2_P
        real*8 Mstar, LumStar, LumEvol
        real*8 smaP, eccP, pressP
        real*8 Porb, omegaORB, omegaPERI, LSP, nu0, Einiz, Miniz, Ls0
        real*8 q0, cporb, cAUPorb	
        real*8 molwtP, cp_P
        real*8 asl, asil, asio
        real*8 fcw, fcl, fci, snowball_LIMIT		
      
	integer NUMBER
        character*7 VERSION   
        character*11 SIMTYPE 
	
        parameter(planet='K452b')
        parameter(TOAalbFile='CCM_RH60/ALB_g16_rh60_co380_ch1.8.txt')   ! file with TOA albedo matrix 
c        OLR file is computed in python script!
 	
c       SIMULATION PARAMETERS

c       LATITUDE GRID
	parameter(N=54)! latitude grid points
        parameter(latbc1=28.) ! latitude limits of the baroclinic zone (Barry et al. 2002)
        parameter(latbc2=68.)
*       CHOOSE ONE OF THE TWO FOLLOWING LINES; COMMENT THE OTHER	
*       parameter(xmin=-1, xmax=1)           ! EQUISPACED IN SIN(LATITUDE)
        parameter(xmin=-pi2/4,xmax=pi2/4)    ! EQUISPACED in LATITUDE
	
        parameter(dx = (xmax-xmin)/N)   
		
	
c       TIME INTEGRATION PARAMETERS
	parameter(Ns=48)               ! number of data in one orbital period
	

c       PLANETARY PARAMETERS 
        parameter(Pmass=1.0) [Earth masses]  !Mass of the planet. Warning, here is WRONG for K452b!

        parameter(obliq=  22.5)   ! [deg] planet axis inclination
	parameter(OBLIQUITY=pi2*obliq/360.)   ! [rad] planet axis inclination  
	 
        parameter(Prot=  1.000)   ! planet rotation period
	parameter(omegaP=pi2/(Prot*86400.))  ! [rad/s] planet rotational angular velocity  
	
        parameter(Rplanet= 1.63*Rearth)   ! planet radius [m] 
        parameter(gP=  1.6*gE)   ! planet gravitational acceleration
	                                    
c       GEOGRAPHY  
        parameter(gg=0)   ! planet geography 
	                ! 0: constant fo; 1: present Earth; 2: equatorial cont.; 3: polar cont. 
        parameter(fo_const= 0.70)   ! constant ocean fraction (if gg=0) 

		
c       OLR MODEL  ['CCMcal00', 'CCMcalCF']
	parameter(OLRmodel='ccm3tab0')
	                                ! CCMcal00: OLR calibrated with CCM, WITHOUT correction factors
	                                ! CCMcalCF: OLR calibrated with CCM, WITH correction factors   
	
c       OUTGOING THERMAL RADIATION  
c       planet CO2 partial pressure     ! ONLY USED if OLRmodel='WK'
        parameter(p_CO2_P= p_CO2_E)       ! [Pa] CO2 partial pressure 
	parameter(pressP=3.0*(pressE+p_CO2_P)) 	 
	
	
c       STELLAR PARAMETERS   
        parameter(Mstar=  1.037*Msun)   ! stellar mass [Msun]
        parameter(LumStar=  1.21*LumSun)   ! stellar luminosity [Lsun]
        parameter(LumEvol=  1.000)   ! evolving stellar luminosity

c       ORBITAL PARAMETERS
        parameter(smaP=1.048) ! semi-major axis of planet orbit [AU]
        parameter(eccP=0.0*earthECC) ! eccentricity of planet orbit
*       parameter(Porb=86400.0d0*365.259635)       ! [s] orbital period
	parameter(cPorb=(pi2/(cgrav*Mstar)**0.5))
        parameter(cAUPorb=cPorb*(AU**1.5))

*       parameter(Porb=86400.0d0*365.259635)       ! [s] orbital period
	parameter(Porb=cAUPorb*(smaP**1.5))   
	
	parameter(omegaORB=pi2/Porb)     ! [rad/s] orbital angular velocity SI
	
        parameter(omegaPERI=0.0) ! argument of pericenter (from ascending node in the orbital plane)

	parameter(LSP=omegaPERI*pi2/360.)  ! [rad] PLANETOCENTRIC ORBITAL LONGITUDE OF THE STAR AT PERICENTER 
        parameter(nu0=-LSP)                ! [rad] initial value of true anomaly (t=0)

c       calculate initial value of Eccentric Anomaly from Eq. (A27)	
        parameter(Einiz=2.*datan(((1.+eccP)/(1.-eccP))**(-0.5)  
     >            *dtan(nu0/2.)))

c       calculate initial value of Mean Anomaly from Eq. (A28)
        parameter(Miniz=Einiz-eccP*dsin(Einiz))  	

c       initial value of the planetocentric stellar longitude	
        parameter(Ls0=pi2/4.)   ! simulation starts at spring equinox

c       INCOMING STELLAR RADIATION at distance equal to the semi-major axis
	parameter(q0=(LumStar*LumEvol/LumSun)*SolarConstant/(smaP**2))  
	
c       ATMOSPHERE 
	parameter(molwtP=28.97) ! Planet atmosphere mean molecular weight
        parameter(cp_P=1.005d3) ! [J/(kg K)] Planet atmosphere specific heat capacity at T=0 C 


c       ALBEDO PARAMETERS
        parameter(asl=0.18)   ! surface albedo of lands (=0.2 in WK97)
	parameter(asil=0.70)   ! surface albedo of ice on lands (=0.85 Pierrehumbert)	
	parameter(asio=asil)   ! surface albedo of ice on ocean (=0.50 Pierrehumbert) 

c       CLOUD COVERAGE		 
	parameter(fcw=0.70) ! cloud coverage on water (mean from Sanroma & Palle 2011)
	parameter(fcl=0.60) ! cloud coverage on land (mean from Sanroma & Palle 2011)
	parameter(fci=fcl)  ! cloud coverage on ice (mean da Sanroma & Palle 2011) 
	parameter(snowball_LIMIT=1.00)  ! when global ice coverage > snowball_LIMIT we set cloud coverage=0



c       DATA FOR Top-of-atmosphere albedo       
        INTEGER NT_TOA,Np_TOA,Nz_TOA,Nas_TOA
        parameter(NT_TOA=19, Np_TOA=16, Nz_TOA=21, Nas_TOA=21) 
        REAL*8 Matrix_TOA(NT_TOA,Np_TOA,Nz_TOA,Nas_TOA)
        REAL*8 T_TOA(NT_TOA),p_TOA(Np_TOA),z_TOA(Nz_TOA),as_TOA(Nas_TOA)




c       DATA FOR THE PRODUCTION OF FIT FILES
        parameter(VERSION="1.1.03")         !code version
        parameter(SIMTYPE="Std")         !simulation type
        parameter(NUMBER=1)         !simulation number for current date

