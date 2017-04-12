*
*       Constants and parameters for codeEBM
* 
        integer N, idebug, Ns
	integer nouti,nprompti,maxNorbits
	integer gg ! type of geography  
	
        real*8 pi2, dx, h, tauir, sigma, T0, a, b, q0
        real*8 omegaORB, omegaE, omegaP
	real*8 earthOBLIQ,earthECC
	real*8 OBLIQUITY,obliq
	real*8 D0par,C0par,C1par,Rpar ! diffusion parameters
        real*8 Porb,cPorb,cAUPorb, AU, cq0
	real*8 pressE, pressEdry, pressP
	real*8 smaP,eccP ! semimajor axis, eccentricity
        real*8 Mstar, Msun
	real*8 cgrav
	real*8 fo_const ! value of constant ocean fraction
        real*8 xmin, xmax
        real*8 Tstart  ! initial temperature 
        real*8 h1, hmin, eps, deltaTconv
	real*8 anglePN,PN,eoc0
	real*8 omtilde,omgrande,meanL
	real*8 omegaPERI,LSP
	real*8 Einiz,Miniz,ompi,nu0
	real*8 SolarConstant,LumSun,LumStar,LumEvol
	real*8 Prot
	real*8 P_CO2_E, P_CO2_P 
        real*8 Ls0
	real*8 fcw,fcl,fci    ! cloud coverage  
	real*8 molwtE,molwtP ! mean molecular weight of the atmosphere
	real*8 cp_P,cp_E     ! specific heat capacity of the atmosphere 
	real*8 column,tauIR0,tauIRf,snowball_LIMIT	
	real*8 molwtCO2,cp_CO2,molwtH2O,cp_H2O
	real*8 asl,asil,asio 
	real*8 fcGLOBAL_E
	real*8 t0_melt  ! minimum time before ice melts at T > 0 C
	real*8 SPindex  ! exponent in Spiegel tauIR: (T/273)^SPindex
	real*8 RH_E,RH ! relative humidity
	real*8 gE, gP, gPnor ! gravitational acceleration
        real*8 Rearth, Rplanet
        real*8 T1bcE ! mean temperature baroclinic zone northern Earth
        real*8 SbcE  ! mean radiative heating Earth baroclinic zone 
        real*8 vTgradE ! mean d(vapour)/dT gradient on Earth
        real*8 DTbcE ! mean temperature gradient north hemisphere
        real*8 latbc1,latbc2 
        real*8 Tcxl1,Tcxl2 ! limits for complex life
        real*8 CML50,COCEAN,CSOLID,CATM_E ! thermal capacity per unit area  
        real*8 cloudOLRforcingE
        real*8 calpha,cbeta  ! cloud albedo parameters (function CloudAlbCess)
        real*8 lambdaE ! ratio of moist over dry atm. flux on Earth
 	character*16 planet,zendistType,albedoType,iceType
	character*16 OLRmodel
        real*8 fixAlbedo
        real*8 Tlim1 ! simulation interrupted when Tmax < Tlim1
	character*50 albedoFile,TOAalbFile
	
        parameter(planet='EARTH')
        parameter(TOAalbFile='CCM_RH60/ALB_g1_rh60_co380_ch1.8.txt')   ! file with TOA albedo matrix 

c       MATHEMATICAL/PHYSICAL CONSTANTS 
	parameter(pi2=6.28318530717958647d0)  ! 2*pi
	parameter(sigma=5.67040d-8)	      ! [W/m2 K4] Stefan-Boltzsmann constant SI
        parameter(cgrav=6.67259d-11)          ! gravitational constant SI

c       ASTRONOMICAL DATA	 
	parameter(AU=1.4959787066d11)               ! [m] astronomical unit SI 
        parameter(SolarConstant=1.360d3)            ! [W/m2] solar constant SI ROUNDED VALUE
c       parameter(SolarConstant=1.3616d3)           ! [W/m2] solar constant SI [Kopp & Lean 2010: 1360.8+0.8]
        parameter(LumSun=SolarConstant*2.*pi2*AU*AU)! Sun luminosity calculated from Solar Constant SI
c       parameter(LumSun=3.845d26)                  ! [W] Sun luminosity SI
	parameter(MSun=1.9891d30)                   ! [kg] Sun mass SI

c       CHEMICAL/PHYSICAL DATA
        parameter(molwtCO2=44.01)       ! CO2 molecular weight 
	parameter(molwtH2O=18.01)       ! H2O molecular weight
	parameter(molwtE=28.97)         ! Earth air mean molecular weight
	parameter(cp_CO2=0.820d3)       ! [J/(kg K)] CO2 specific heat capacity at T=0 C 
	parameter(cp_H2O=1.847d3)       ! [J/(kg K)] H2O specific heat capacity at T=0 C  	
	parameter(cp_E=1.005d3)         ! [J/(kg K)] Earth air specific heat capacity at T=0 C  

c       THERMAL CAPACITY PER UNIT AREA OF DIFFERENT TYPES OF SURFACE
c       see Pierrehumbert(2010,p.445)  and  WK97
        parameter(CSOLID=1.d6)           ! [J/(m2 K)]  ~0.5m rock/ice  	 
        parameter(CML50=2.10d8)          ! [J/(m2 K)] 50m mixed ocean layer   
        parameter(CATM_E=CML50*(2.4/50)) ! [J/(m2 K)] Earth atmosph.=2.4m mixed ocean layer 

c       EARTH DATA
        parameter(earthECC=0.01671022)  ! Current eccentricity
	parameter(earthOBLIQ=23.43929)  ! Current obliquity of the axis  NOT USED
	parameter(omegaE=pi2/86400.d0)  ! [rad/s] Rotational angular velocity SI
        parameter(pressE=1.01325d5)     ! [Pa] Earth surface pressure SI
	parameter(pressEdry=1.0031d5)   ! [Pa] Earth dry surface pressure SI
        parameter(p_CO2_E=3.8d1)        ! [Pa] CO2 partial pressure SI 
	parameter(Rearth=6.371e6)       ! [m] Earth mean radius SI
        parameter(gE=9.80)              ! [m/s2] Gravitational acceleration SI

        parameter(DTbcE=26.6835)        ! [K] Temperature difference (lat=+28,+68) in Earth 
        parameter(T1bcE=293.3641)       ! [K] Mean T (lat=+28) in Earth
        parameter(SbcE=2.06450d2)       ! [W/m2] Mean absorbed radiation in the baroclinic zone 
        parameter(vTgradE=74.43897)     ! [Pa K-1] Mean d(vapour)/dT gradient	
        parameter(cloudOLRforcingE=26.4) ! [W/m2] Mean global cloud OLR forcing Pierrehumbert (2010)
	parameter(fcGLOBAL_E=0.667)     ! Global cloud coverage  
        parameter(lambdaE=0.7)          ! ratio of moist over dry atm. flux on Earth

        parameter(calpha=-0.11)         ! cloud albedo parameter (function CloudAlbCess)
        parameter(cbeta=7.98d-3)         ! cloud albedo parameter (function CloudAlbCess) 
        parameter(zendistType='orbital')  ! mean orbital zenith distance
c       parameter(zendistType='instant')  ! instantaneous zenith distance
 	



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
	parameter(maxNorbits=100)      ! max integration time (number of orbital periods)     
	parameter(deltaTconv=1.0d-2)   ! simulation stops when annual global temperature converges within this accuracy
	
	parameter(h1 = 80000.0)
	parameter(hmin=60.) 
	parameter(eps=1e-6) 
	
c       PRINT OPTIONS	 
	parameter(nouti=2)             ! numero di output in totale
	parameter(nprompti=10)         ! numero di calcoli della temp media a schermo
	parameter(idebug=0)            ! set 1 for print debugging 
	
c       TEMPERATURES 
        parameter(Tstart=  275) ! initial temperature [HOT START: 350 K]
        parameter(Tlim1=200.d0)      ! [K] SIMULATION INTERRUPTED WHEN Tmax < Tlim1

c       STELLAR PARAMETERS   
        parameter(Mstar=  1.000*Msun)   ! stellar mass [Msun]
        parameter(LumStar=  1.0000*LumSun)   ! stellar luminosity [Lsun]
        parameter(LumEvol=  1.000)   ! evolving stellar luminosity

c       ORBITAL PARAMETERS
        parameter(smaP=1.00000) ! semi-major axis of planet orbit [AU]
        parameter(eccP=0.01671) ! eccentricity of planet orbit

*       parameter(Porb=86400.0d0*365.259635)       ! [s] orbital period
	parameter(cPorb=(pi2/(cgrav*Mstar)**0.5))
        parameter(cAUPorb=cPorb*(AU**1.5))
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

c       PLANETARY PARAMETERS 

        parameter(obliq= 23.4393)   ! [deg] planet axis inclination
	parameter(OBLIQUITY=pi2*obliq/360.)   ! [rad] planet axis inclination  
	
        parameter(Prot=  1.000)   ! planet rotation period
	parameter(omegaP=pi2/(Prot*86400.))  ! [rad/s] planet rotational angular velocity  
	
        parameter(Rplanet= 1.0000*Rearth)   ! planet radius [m] 
        parameter(gP=  9.80)   ! planet gravitational acceleration
	                                    
c       GEOGRAPHY  
        parameter(gg=4)   ! planet geography 
	                ! 0: constant fo; 1: present Earth; 2: equatorial cont.; 3: polar cont. 
        parameter(fo_const= 0.70)   ! constant ocean fraction (if gg=0) 

c       THERMAL CAPACITY PER UNIT AREA OF THE OCEANS
        parameter(COCEAN=CML50) ! [J/(m2 K)] thermal capacity per unit area of the ocean
	
c       ALBEDO
        parameter(albedoFile='none')   ! file with albedo data 
        parameter(albedoType='aut')  ! albedo is calculated from the surface and atmospheric properties
c       parameter(albedoType='fix')  ! the planet albedo is fixed
        parameter(fixAlbedo=0.35)    ! fixed value of albedo

c       ICE
        parameter(iceType='aut')   ! ice coverage is calculated from the surface temperature
c       parameter(iceType='none')  ! ice coverage is set to zero

	
c       OLR MODEL  ['CCMcal00', 'CCMcalCF']
	parameter(OLRmodel='ccm3tab0')
	                                ! CCMcal00: OLR calibrated with CCM, WITHOUT correction factors
	                                ! CCMcalCF: OLR calibrated with CCM, WITH correction factors   
	
c       OUTGOING THERMAL RADIATION  
c       planet CO2 partial pressure     ! ONLY USED if OLRmodel='WK'
        parameter(p_CO2_P= 3.80e+01)       ! [Pa] CO2 partial pressure 
	
c       ATMOSPHERE 
	parameter(column=1.00)  ! vertical column density of the atmosphere [Earth=1.0] 
        parameter(pressP=1.00310e+05)   ! [Pa] planet total pressure
        parameter(molwtP=28.97)   ! Planet atmosphere mean molecular weight
        parameter(cp_P=1005.000)   ! Planet atmosphere specific heat capacity at T=0 C
	parameter(RH_E=0.6)     ! relative humidity of earth model
        parameter(RH=RH_E)      ! relative humidity of the planet   
		 
c       DIFFUSION PARAMETERS
        parameter(D0par=0.660)         ! [W/m2] diffusion parameter D0
        parameter(C0par= 0.52598)         ! diffusion parameter C0
        parameter(C1par= 0.99144)         ! diffusion parameter C1
        parameter(Rpar=  2.2)         ! diffusion parameter R
	
c       ALBEDO PARAMETERS
        parameter(asl= 0.18)         ! surface albedo of lands (=0.2 in WK97)
        parameter(asil= 0.70)         ! surface albedo of ice on lands (=0.85 Pierrehumbert)
        parameter(asio= 0.70)         ! surface albedo of ice on ocean (=0.50 Pierrehumbert)

c       CLOUD COVERAGE		 
        parameter(fcw= 0.70)         ! cloud coverage on water (mean from Sanroma & Palle 2011)
        parameter(fcl= 0.60)         ! cloud coverage on land (mean from Sanroma & Palle 2011)
        parameter(fci= 0.60)         ! cloud coverage on ice (mean from Sanroma & Palle 2011)
	parameter(snowball_LIMIT=1.00)  ! when global ice coverage > snowball_LIMIT we set cloud coverage=0

c       ASTROBIOLOGICAL PARAMETERS
        parameter(Tcxl1=273.15)      ! [K] minimum temperature for complex life
        parameter(Tcxl2=Tcxl1+50.)   ! [K] maximum temperature for complex life


c       DATA FOR Top-of-atmosphere albedo       
        INTEGER NT_TOA,Np_TOA,Nz_TOA,Nas_TOA
        parameter(NT_TOA=19, Np_TOA=16, Nz_TOA=21, Nas_TOA=21) 
        REAL*8 Matrix_TOA(NT_TOA,Np_TOA,Nz_TOA,Nas_TOA)
        REAL*8 T_TOA(NT_TOA),p_TOA(Np_TOA),z_TOA(Nz_TOA),as_TOA(Nas_TOA)

c       Array of OLR values
	     integer NOLR
        parameter(NOLR=71)
        real*8 TOLR(NOLR),vOLR(NOLR)
        data TOLR/100.0,105.0,110.0,115.0,120.0,125.0,
     >     130.0,135.0,140.0,145.0,150.0,155.0,
     >     160.0,165.0,170.0,175.0,180.0,185.0,
     >     190.0,195.0,200.0,205.0,210.0,215.0,
     >     220.0,225.0,230.0,235.0,240.0,245.0,
     >     250.0,255.0,260.0,265.0,270.0,275.0,
     >     280.0,285.0,290.0,295.0,300.0,305.0,
     >     310.0,315.0,320.0,325.0,330.0,335.0,
     >     340.0,345.0,350.0,355.0,360.0,365.0,
     >     370.0,375.0,380.0,385.0,390.0,395.0,
     >     400.0,405.0,410.0,415.0,420.0,425.0,
     >     430.0,435.0,440.0,445.0,450.0/ 
        data vOLR/  5.6,  6.8,  8.1,  9.7, 11.4, 13.3,
     >      15.5, 17.9, 20.6, 23.6, 26.8, 30.4,
     >      34.3, 38.5, 43.1, 48.1, 53.6, 59.4,
     >      65.7, 72.4, 79.6, 87.3, 95.4,103.9,
     >     112.8,122.2,131.9,142.1,152.5,163.3,
     >     174.4,185.7,197.2,209.0,220.9,232.9,
     >     244.9,256.5,268.0,278.8,288.6,297.2,
     >     304.9,312.4,320.0,327.7,335.2,340.2,
     >     349.2,355.5,361.3,367.1,373.2,379.6,
     >     386.5,394.0,402.4,411.7,422.2,434.0,
     >     447.4,462.5,479.9,499.4,521.8,547.1,
     >     576.0,610.6,622.4,684.0,688.3/ 
