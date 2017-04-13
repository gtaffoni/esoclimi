c       inizialization of parameters in barry's transport coefficient
c       should be set equal to the EARTH's VALUES 

      DTbc=DTbcE
      T1bc=T1bcE
      Sbc=SbcE

c       LATITUDE GRID 
c       dx, xmin are read from parameters file parEBM.h

      do i=1,N 
         x(i)= (xmin+dx*(i-1)+dx/2.0) 
         if(xmin < -1.5) then !integration in phi
            phiDEG(i)=x(i)*360./pi2
	    phi(i)=x(i)
         else                 !integration in sin(phi)
            phiDEG(i)=dasin(x(i))*360./pi2
	    phi(i)=dasin(x(i))
         endif 
      end do   

c       Mean orbital zonal zenith distance 

      do i=1,N ! loop latitudini 
         HHS(I)=0.
         OrbitMeanZD=0.d0 
         orbitdays=Porb/(Prot*86400.)
         if(orbitdays.lt.3) then
            print *,'WARNING: LESS THAN 3 PLANET ROTATIONS IN ONE ORBIT'
         end if  
         do j = 1,orbitdays     ! loop anno
            time=dfloat(j)*Prot*86400. 
            DMM=DiurnalMeanMu(time,x(i),HourD)
            HHS(I)=HHS(I)+HourD
            OrbitMeanZD=OrbitMeanZD+
     >           HourD*dacos(DMM)*360./pi2  
         end do                 ! fine loop sull'anno   
         OrbitMeanZD=OrbitMeanZD/HHS(I)
c     MEDIA SU TUTTI I GIORNI DELL'ANNO, PESANDO CON LE ORE DI LUCE
         mzd(i)= OrbitMeanZD 
c     print *,phiDEG(i),mzd(i), HHS(I) 
      end do                    ! fine loop latitudini 
      
c     SET INITIAL VALUES   
      sigmaLIM=0.5              !1.d-4   ! maximum fraction of planet surface in critical conditions od water vapor
      exitFLAG=1.0              ! flag that describes the type of exit: 
                                !  1=regular exit
                                ! -0.5=Tmax below minimum tamperature (in practice outer edge); 
                                ! -1=critical conditions of water vapor (at inner edge)
                                ! -100=simulation interrupted because stepsize too small

      cmH2Olim=0.1d5/1.         ! VALORE MASSIMO DI COLUMNAR MASS DI VAPORE ACQUEO
                                ! unita' di misura: Pa / gravita' normalizzata al valore terrestre
      
      time0=0.0d0               ! initial time 
      t1=0.0 

      do i=1,N 
         f(i)=Tstart            ! zonal temperatures
      enddo  
      tsum_old=Tstart
      
      do j=1,N   
         Htime(j)=0.
      enddo 
      

      do j=1,N
         do is=1,Ns
            tempmat(is,j)=Tstart ! initialize temperature matrix
            habmat(is,j)=0.     ! initialize habitability function matrix
            cxlmat(is,j)=0.
            fcmat(is,j)=0.5     ! initialize cloud cover matrix
            asrmat(is,j)=q0/4   ! initialize insolation matrix 
            atmfluxmat(is,j)=0. ! initialize atmospheric flux matrix
            dryfluxmat(is,j)=0.
         end do
      end do
      
      nprompt = nprompti        ! number of orbits between printouts 
      
      pressPtot=pressP          ! total pressure without water vapor
      cp_Ptot=cp_P              ! specific heat capacity without water vapor
      molwtPtot=molwtP          ! molecular weight without water vapor
      
      do i=1,N
         zonalfc(i)=0.5         ! annual zonal cloud cover
      enddo  
      

c     PLANET GEOGRAPHY  
c     parameter gg is read from parameters file parEBM.h

      if(gg.eq.0) then          !  constant ocean fraction in all zones
         do i=1,N  
            fo(i)=fo_const
         enddo
         fom=fo_const
         write(*,500) fo_const,N
      end if 
      
      if(gg.ge.1.and.gg.le.3) then !  zonal ocean fractions WK97 
         call geography(x,fo,fom)
         if(gg.eq.1) write(*,501) N,fom
         if(gg.eq.2) write(*,502) N,fom 
         if(gg.eq.3) write(*,503) N,fom
      end if 
      
      if(gg.eq.4) then 
         OPEN(65,file='fo_earth_DMAP.dat',status='old') ! ocean fraction file   
         do i=1,nlatDMAP        !  nlat=46 in file  'fo_earth_DMAP.dat'
            READ(65,6565) latDMAP(i),foDMAP(i) 
         end do
         CLOSE(65)
         fom=0.7095 
         do i=1,N  
            fo(i)=lininterp(latDMAP,foDMAP,nlatDMAP,phiDEG(i))
c     print *,phiDEG(i),fo(i)
         end do
         write(*,504) N,fom
      end if
 6565 format(f5.1,2x,f5.3)   

       
 500  format('GEOGRAPHY: constant fo=',f6.4,' Nlat=',i3)
      
 501  format('GEOGRAPHY: present Earth         Nlat=',
     >     i3,5x,'<fo>=',f6.4)
 502  format('GEOGRAPHY: equatorial continent  Nlat=',
     >     i3,5x,'<fo>=',f6.4)
 503  format('GEOGRAPHY: polar continent       Nlat=',
     >     i3,5x,'<fo>=',f6.4)  
 504  format('GEOGRAPHY: present Earth DMAP    Nlat=',
     >     i3,5x,'<fo>=',f6.4)
      
