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
	
        parameter(planet='EARTH')
        parameter(TOAalbFile='CCM_RH60/ALB_g1_rh60_co380_ch1.8.txt')   ! file with TOA albedo matrix 
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
        parameter(Pmass=1.0)  ! [Earth masses ]Mass of planet

        parameter(obliq=  earthOBLIQ)   ! [deg] planet axis inclination
	parameter(OBLIQUITY=pi2*obliq/360.)   ! [rad] planet axis inclination  
	 
        parameter(Prot=  1.000)   ! planet rotation period
	parameter(omegaP=pi2/(Prot*86400.))  ! [rad/s] planet rotational angular velocity  
	
        parameter(Rplanet= 1.0000*Rearth)   ! planet radius [m] 
        parameter(gP=  9.80)   ! planet gravitational acceleration
	                                    
c       GEOGRAPHY  
        parameter(gg=4)   ! planet geography 
	                ! 0: constant fo; 1: present Earth; 2: equatorial cont.; 3: polar cont. 
        parameter(fo_const= 0.70)   ! constant ocean fraction (if gg=0) 

		
c       OLR MODEL  ['CCMcal00', 'CCMcalCF']
	parameter(OLRmodel='ccm3tab0')
	                                ! CCMcal00: OLR calibrated with CCM, WITHOUT correction factors
	                                ! CCMcalCF: OLR calibrated with CCM, WITH correction factors   
	
c       OUTGOING THERMAL RADIATION  
c       planet CO2 partial pressure     ! ONLY USED if OLRmodel='WK'
        parameter(p_CO2_P= p_CO2_E)       ! [Pa] CO2 partial pressure 
	parameter(pressP=pressE+p_CO2_P) 	 
	
	
c       STELLAR PARAMETERS   
        parameter(Mstar=  1.000*Msun)   ! stellar mass [Msun]
        parameter(LumStar=  1.0000*LumSun)   ! stellar luminosity [Lsun]
        parameter(LumEvol=  1.000)   ! evolving stellar luminosity

c       ORBITAL PARAMETERS
        parameter(smaP=1.00000) ! semi-major axis of planet orbit [AU]
        parameter(eccP=earthECC) ! eccentricity of planet orbit
*       parameter(Porb=86400.0d0*365.259635)       ! [s] orbital period
	parameter(cPorb=(pi2/(cgrav*Mstar)**0.5))
        parameter(cAUPorb=cPorb*(AU**1.5))

*       parameter(Porb=86400.0d0*365.259635)       ! [s] orbital period
	parameter(Porb=cAUPorb*(smaP**1.5))   
	
	parameter(omegaORB=pi2/Porb)     ! [rad/s] orbital angular velocity SI
	
        parameter(omegaPERI=-77.06300) ! argument of pericenter (from ascending node in the orbital plane)

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
        parameter(asl=0.20)   ! surface albedo of lands (=0.2 in WK97)
	parameter(asil=0.85)   ! surface albedo of ice on lands (=0.85 Pierrehumbert)	
	parameter(asio=0.50)   ! surface albedo of ice on ocean (=0.50 Pierrehumbert) 

c       CLOUD COVERAGE		 
	parameter(fcw=0.70) ! cloud coverage on water (mean from Sanroma & Palle 2011)
	parameter(fcl=0.50) ! cloud coverage on land (mean from Sanroma & Palle 2011)
	parameter(fci=fcl)  ! cloud coverage on ice (mean da Sanroma & Palle 2011) 
	parameter(snowball_LIMIT=1.00)  ! when global ice coverage > snowball_LIMIT we set cloud coverage=0



c       DATA FOR THE PRODUCTION OF FIT FILES
        parameter(VERSION="1.1.03")         !code version
        parameter(SIMTYPE="Std")         !simulation type
        parameter(NUMBER=1)         !simulation number for current date

c       DATA FOR Top-of-atmosphere albedo       
        INTEGER NT_TOA,Np_TOA,Nz_TOA,Nas_TOA
        parameter(NT_TOA=19, Np_TOA=16, Nz_TOA=21, Nas_TOA=21) 
        REAL*8 Matrix_TOA(NT_TOA,Np_TOA,Nz_TOA,Nas_TOA)
        REAL*8 T_TOA(NT_TOA),p_TOA(Np_TOA),z_TOA(Nz_TOA),as_TOA(Nas_TOA)



c       Array of OLR values
c             integer NOLR
c        parameter(NOLR=71)
c        real*8 TOLR(NOLR),vOLR(NOLR)
c        data TOLR/100.0,105.0,110.0,115.0,120.0,125.0,
c     >     130.0,135.0,140.0,145.0,150.0,155.0,
c     >     160.0,165.0,170.0,175.0,180.0,185.0,
c     >     190.0,195.0,200.0,205.0,210.0,215.0,
c     >     220.0,225.0,230.0,235.0,240.0,245.0,
c     >     250.0,255.0,260.0,265.0,270.0,275.0,
c     >     280.0,285.0,290.0,295.0,300.0,305.0,
c     >     310.0,315.0,320.0,325.0,330.0,335.0,
c     >     340.0,345.0,350.0,355.0,360.0,365.0,
c     >     370.0,375.0,380.0,385.0,390.0,395.0,
c     >     400.0,405.0,410.0,415.0,420.0,425.0,
c     >     430.0,435.0,440.0,445.0,450.0/ 
c        data vOLR/  5.6,  6.8,  8.1,  9.7, 11.4, 13.3,
c     >      15.5, 17.9, 20.6, 23.6, 26.8, 30.4,
c     >      34.3, 38.5, 43.1, 48.1, 53.6, 59.4,
c     >      65.7, 72.4, 79.6, 87.3, 95.4,103.9,
c     >     112.8,122.2,131.9,142.1,152.5,163.3,
c     >     174.4,185.7,197.2,209.0,220.9,232.9,
c     >     244.9,256.5,268.0,278.8,288.6,297.2,
c     >     304.9,312.4,320.0,327.7,335.2,340.2,
c     >     349.2,355.5,361.3,367.1,373.2,379.6,
c     >     386.5,394.0,402.4,411.7,422.2,434.0,
c     >     447.4,462.5,479.9,499.4,521.8,547.1,
c     >     576.0,610.6,622.4,684.0,688.3/ 
