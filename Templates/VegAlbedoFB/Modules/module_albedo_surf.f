*******************************************************************
*       calculate mean zonal surface albedo
*******************************************************************
      subroutine albsup(fo,fio,fil,aso,asc,as) 
      include 'parEBM.h'  ! read fcw, fcl, fci from parameters file 
      real*8 fo,fio,fil,aso,asc  
      real*8 asc_ice ! albedo of clouds over ice
      real*8 as ! output: total surface albedo as 

      include 'vegetation.h' 
      common/vegetation/veg,fveg

      fc_ice=fci 
      asc_ice=asil ! albedo of clouds over ices = albedo of ice
      
      as=  fo*(1.-fio)*                ! OCEANI non ghiacciati
     >    (aso*(1.-fcw) + asc*fcw) +   ! senza nubi e con nubi
     >     fo*fio*                     ! OCEANI GHIACCIATI
     >    (asio*(1.-fc_ice) + asc_ice*fc_ice) +  ! senza e con nubi
     >     (1.-fo)*(1.-fil)*           ! TERRE non ghiacciate
     >    ( asl*(1.-fcl)*(1-fveg) +    ! senza nubi, senza vegetazione
     >    aveg*(1.-fcl)*fveg +        ! senza nubi, con vegetazione
     >    asc*fcl) +                   ! con nubi (fcl puo' variare sulla vegetazione??)
     >     (1.-fo)*fil*                ! TERRE GHIACCIATE
     >    (asil*(1.-fc_ice) + asc_ice*fc_ice)    ! senza e con nubi   


! EARTH.py
!asl=0.18             # surface albedo of lands (=0.2 in WK97)
!asil=0.70            # surface albedo of ice on lands (=0.85 Pierrehumbert)
!asio=asil            # surface albedo of ice on ocean (=0.50 Pierrehumbert)

!fcw=0.70              # cloud coverage on water (mean from Sanroma & Palle 2011=0.70)
!fcl=0.60              # cloud coverage on land (mean from Sanroma & Palle 2011=0.50)
!fci=fcl               # cloud coverage on ice (mean da Sanroma & Palle 2011=0.50
            
!fo ocean fraction
!fio fraction of ice on ocean
!fil fraction of ice on lands
!aso  ocean albedo - computed in EBMfunctions.f
!asc  cloud albedo - computed in EBMfunctions.f
!aveg    vegetation albedo (in vegetation.h under VegAlbedoFB)
!fveg    vegetaton fraction (in module_vegalbedo_local.f)

!asc_ice, fc_ice - see above 

     
 
      return
      end 

***************** the following is for the computation of surface albedo, INCLUDING ices but not clouds***


      
      subroutine albsup_noclouds(fo,fio,fil,aso,asc,as) 
      include 'parEBM.h'  ! read fcw, fcl, fci from parameters file 
      real*8 fo,fio,fil,aso,asc  
      real*8 asc_ice ! albedo of clouds over ice
      real*8 as ! output: total surface albedo as 

      real*8 fcwn, fcln, fcin

      fcwn=0.    !NOTE: this eliminates the clouds
      fcln=0.
      fcin=0.

      
      fc_ice=fcin 
      asc_ice=asil ! albedo of clouds over ices = albedo of ice
      


      as=  fo*(1.-fio)*                ! OCEANI non ghiacciati
     >    (aso*(1.-fcw) + asc*fcw) +   ! senza nubi e con nubi
     >     fo*fio*                     ! OCEANI GHIACCIATI
     >    (asio*(1.-fc_ice) + asc_ice*fc_ice) +  ! senza e con nubi
     >     (1.-fo)*(1.-fil)*           ! TERRE non ghiacciate
     >    ( asl*(1.-fcl)*(1-fveg) +    ! senza nubi, senza vegetazione
     >    aveg*(1.-fcl)*fveg +        ! senza nubi, con vegetazione
     >    asc*fcl) +                   ! con nubi (fcl puo' variare sulla vegetazione??)
     >     (1.-fo)*fil*                ! TERRE GHIACCIATE
     >    (asil*(1.-fc_ice) + asc_ice*fc_ice)    ! senza e con nubi   
      
            
      return
      end 


      
      real*8 function alb_Zsurf(hour) 

      include 'parEBM.h' ! read albedo parameters asl,asio,asil
      real*8 hour
      real*8 phi,delta,T,time
      real*8 mu, ZZ
      real*8 f_ice, fi, as, aso, asi, asc 
      integer ilat
      real*8 fo(N) 
      real*8 inc,nfr,rfr  ! Fresnel parameters 
      real*8 alfa, beta  
      real*8 TOA_albedo,TOA_diff,CloudAlbCess  
      real*8 fio,fil
      real*8 ToaQuadLinInterp,pi,Zdeg 
      real*8 mzd(N)  ! MEAN ORBITAL ZENITH DISTANCE
      
      common /albparam/ time,phi,delta,fo,T,ilat 
      common /interp_TOA/ Matrix_TOA,T_TOA,p_TOA,z_TOA,as_TOA 
      common /meanZenDist/ mzd 

      include 'module_vegalbedo_local.f'

*       calculate instantaneous ZENITH DISTANCE   Eq. (A20) 
      mu= dsin(phi)*dsin(delta) + dcos(phi)*dcos(delta)*dcos(hour)  
      
      if (abs(mu) .lt. 1d-10) mu=0.
      
      ZZ=dacos(mu)  

      pi=3.1415926535
      Zdeg=ZZ*180./pi ! convert to degrees 

*       calculate OCEAN albedo at current zenith distance  Eq. (A14)
*       Briegleb et al,J. Clim. Appl. Meteorol. 25, 214–224 (1986) 
       aso = 0.026/(1.1*mu**1.7 + 0.065) 
     >      + 0.15*(mu-0.1)*(mu-0.5)*(mu-1.0)  ! OK CONTROLLATA  

*        ice cover
      fio=f_ice(T,time,ilat) ! ocean ice cover
      fil=fio                ! land ice cover

*        cloud albedo at mean zonal annual zenith distance 
*        based on Cess 1976 trend 
      if(zendistType.eq.'instant') asc=dmax1(0.,calpha+cbeta*Zdeg) 
      if(zendistType.eq.'orbital') asc=calpha+cbeta*mzd(ilat)  
      
*       calculate mean zonal surface albedo Eq. (A13)
*       values of land and ice albedos are read from parEBM.h
*       cloud albedo is treated as part of the surface albedo
      call albsup_noclouds(fo(ilat),fio,fil,aso,asc,as) 

*       calculate Top-of-Atmosphere albedo as a function of:
*       T, p(DRY), zenith distance and surface albedo  
      alb_Zsurf=as 
      return
      end

      
