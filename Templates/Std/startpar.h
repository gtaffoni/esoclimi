*
*       Constants and parameters for codeEBM
* 
        character*50 albedoFile
    
        real*8 pi2, sigma, cgrav, AU, SolarConstant, LumSun, Msun
        real*8 molwtCO2, molwtH2O, molwtE, cp_CO2, cp_H2O, cp_E
        real*8 CSOLID, CML50, CATM_E
        real*8 earthECC, earthOBLIQ, omegaE
        real*8 pressE, pressEdry, p_CO2_E, Rearth, gE
        real*8 DTbce, T1bce, SbcE, vTgradE, cloudOLRforcingE
        real*8 fcGLOBAL_E, lambdaE
        real*8 calpha, cbeta
        real*8 xmin, xmax, dx
        integer maxNorbits, nouti, nprompti, idebug
        real*8 deltaTconv, ddeltaTconv, h1, hmin,eps
        real*8 Tstart, Tlim1
        real*8 COCEAN
        character*16 albedoType,iceType, zenDistType
        real*8 fixAlbedo, column, RH_E, RH
        real*8 D0par, C0par, C1par, Rpar
        real*8 Tcxl1, Tcxl2
	
c       MAT/maxNorHEMATICAL/PHYSICAL CONSTANTS 
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
        parameter(p_CO2_E=38.0)         ! [Pa] CO2 partial pressure SI 
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

	
c       TIME INTEGRATION PARAMETERS
        parameter(maxNorbits=1000)      ! max integration time (number of orbital periods)
        parameter(deltaTconv=1.0d-4)   ! simulation stops when DELTA(annual global temperature)/annual global temperature
*                                        converges within this accuracy
        parameter(ddeltaTconv=1.d-4)   ! also the derivative of DeltaT/T is checked
	
	parameter(h1 = 80000.0)
	parameter(hmin=60.) 
	parameter(eps=1e-6) 
	
c       PRINT OPTIONS	 
	parameter(nouti=2)             ! numero di output in totale
	parameter(nprompti=10)         ! numero di calcoli della temp media a schermo
	parameter(idebug=0)            ! set 1 for print debugging 
	
c       TEMPERATURES 
        parameter(Tstart=275.d0)     ! [K] start T; T=350 for HOT START (Spiegel09 page 598)
        parameter(Tlim1=200.d0)      ! [K] SIMULATION INTERRUPTED WHEN Tmax < Tlim1

	

c       THERMAL CAPACITY PER UNIT AREA OF THE OCEANS
        parameter(COCEAN=CML50) ! [J/(m2 K)] thermal capacity per unit area of the ocean
	 
c       ALBEDO
        parameter(albedoFile='none')
        parameter(albedoType='aut')  ! albedo is calculated from the surface and atmospheric properties
c       parameter(albedoType='fix')  ! the planet albedo is fixed
        parameter(fixAlbedo=0.35)    ! fixed value of albedo

c       ICE
        parameter(iceType='aut')   ! ice coverage is calculated from the surface temperature
c       parameter(iceType='none')  ! ice coverage is set to zero

		
	
c       ATMOSPHERE 
	parameter(column=1.00)  ! vertical column density of the atmosphere [Earth=1.0] 
	parameter(RH_E=0.6)     ! relative humidity of earth model
        parameter(RH=RH_E)      ! relative humidity of the planet   
		 
c       DIFFUSION PARAMETERS
        parameter(D0par=0.58d0)  ! [W/m2] 
	parameter(C0par=0.39972) 
	parameter(C1par=1.25577)  
	parameter(Rpar=5.0)
	

c       ASTROBIOLOGICAL PARAMETERS
        parameter(Tcxl1=273.15)      ! [K] minimum temperature for complex life
        parameter(Tcxl2=Tcxl1+50.)   ! [K] maximum temperature for complex life



