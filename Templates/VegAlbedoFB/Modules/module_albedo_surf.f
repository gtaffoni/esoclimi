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

