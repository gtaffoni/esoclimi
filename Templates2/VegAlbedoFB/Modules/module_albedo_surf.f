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
     >    aveg*(1.-fvc)*fveg +        ! senza nubi, con vegetazione
     >    asc*fcl) +                   ! con nubi (fcl puo' variare sulla vegetazione??)
     >     (1.-fo)*fil*                ! TERRE GHIACCIATE
     >    (asil*(1.-fc_ice) + asc_ice*fc_ice)    ! senza e con nubi   

            
      return
      end 

