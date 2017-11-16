************************************************************
*   this updates seasonal vector and matrix                 *
************************************************************
      subroutine update_seasonal_qties(is, t2, pressPtot, 
     >     cp_Ptot, molwtPtot, tclx1, tclx2,
     >     cmH2O, cmH2Olim)
      

      implicit none
      integer is
      real*8 t2, PressPtot, cp_Ptot, 
     >     molwtPtot, tclx1, 
     >     tclx2,cmH2O, cmH2Olim


      real*8 Tice, Tvapor, dfi, Ddry
      integer j

      include 'parEBM.h'
      
      include 'matrices.h'
      include 'functions.h'
      common/matrices/f, fo,  fcmat, olrmat, albmat, asrmat,
     >     phi, dryfluxmat, atmfluxmat, tts, habmat, cxlmat, 
     >     boilmat, RGmat
      common /tempmatrix/ tempmat


	  do j=1,N
             tempmat(is,j)=f(j) ! update temperature matrix used in f_ice
             call subfc(f(j),t2,j,fo(j),fcmat(is,j)) ! update matrix of cloud cover fcmat
             olrmat(is,j)=Iterm(j,f(j))
             albmat(is,j)=Aterm(t2,dsin(phi(j)),f(j),j,fo)
             asrmat(is,j)=Sterm(t2,phi(j))
     >            *(1.-Aterm(t2,dsin(phi(j)),f(j),j,fo))
	  end do 
	            

c             Poleward atmospheric enery flux [W]
c             based on the definition of diffusion coefficient D
          do j=1,N-1
             atmfluxmat(is,j)=-pi2*Rplanet**2
     >            *Dterm(t2,phi(j),f(j),Ddry)
     >            *(f(j+1)-f(j))/(phi(j+1)-phi(j)) ! local latitude temperature gradient
     >            *dcos(phi(j))
             dryfluxmat(is,j)=-pi2*Rplanet**2
     >            *Ddry
     >            *(f(j+1)-f(j))/(phi(j+1)-phi(j)) ! local latitude temperature gradient
     >            *dcos(phi(j))
          end do
          atmfluxmat(is,N)=0.
          dryfluxmat(is,N)=0.
          
c             global temperature at current season
          do j=1,N 
             dfi=dlat(j)
             tts(is) = tts(is) + dfi*f(j)/2.0 ! area weighted temperature of current season
          end do

c     update atmospheric quantities with global water vapor content at current season
      
          pressPtot = pressP + RH*P_H2O(tts(is)) ! total pressure
    
          cp_Ptot =(pressP*cp_P+RH*P_H2O(tts(is))*cp_H2O) ! thermal capacity 
     >         /pressPtot
      
          molwtPtot=(pressP*molwtP+RH*P_H2O(tts(is))*molwtH2O) ! molecular weight 
     >         /pressPtot    
      
          call LWTR(pressPtot,Tice,Tvapor) ! liquid water temperature range 

c---------------------------------------------------------------------------	  
	  do j=1,N 
             habmat(is,j)=0.    ! update habitability function
             cxlmat(is,j)=0.
             boilmat(is,j)=0.   !  boiling point matrix
             RGmat(is,j)=0.     !   water vapour column limit matrix
             if( f(j) .ge. Tice .and. f(j) .le. Tvapor )then
                habmat(is,j)=1.0
             end if
             if(f(j).ge.Tcxl1 .and. f(j).le.dmin1(Tcxl2,Tvapor)) then
                cxlmat(is,j)=1.0
             end if
             if( f(j) .gt. Tvapor ) then !  20 dic 2013
                boilmat(is,j)=1.0
             endif
             cmH2O=RH*P_H2O(f(j))/(gP/gE) ! COLUMNAR MASS of WATER VAPOUR
                                          ! units:   Pa / (gravity in earth units)
             if( cmH2O .gt. cmH2Olim ) then !  22 dic 2013
                RGmat(is,j)=1.0
             endif
          end do 


      return
      end



************************************************************
*   this writes on file needed  seasonal data              *
************************************************************
      subroutine write_seasons(i, is, nlastorbits, t2, f, fo, phi, 
     >     phiDEG)
      
      implicit none
      integer i, is, nlastorbits
      include 'parEBM.h'
      real*8 f(N), fo(N), t2
      real*8 phi(N) ! latitude in radians 
      real*8 phiDEG(N) ! latitude in degrees 

      integer j

      real*8 Aterm, f_ice

      do j=1,N  
         write(19,30) (i+is/dfloat(Ns)),phiDEG(j),f(j) 
 30      format(f14.8,1x,f14.8,1x,f14.8)
      enddo	   
      
c     store season-latitude-temperature of the last 5 orbits after convergence
      if( nlastorbits .ge.1 .and. nlastorbits .le. 5) then   
         do j=1,N
            write(23,30) (i+is/dfloat(Ns)),phiDEG(j),f(j) 
         enddo
      endif
      
c     store season-latitude-temperature of the last orbit of the simulations (5th orbit after convergence)
      if( nlastorbits .eq. 5 ) then 
         do j=1,N 
            write(27,30) (i+is/dfloat(Ns)),phiDEG(j),f(j)
         enddo
      endif 
      
*     store season-latitude-albedo and other results of the lst two orbits
      if(albedoFile.eq.'none') then
         if( nlastorbits .eq. 5) then ! solo l'ultimo anno 
            do j=1,N
               write(40,9430) (i+is/dfloat(Ns)),phiDEG(j),
     >              Aterm(t2,dsin(phi(j)),f(j),j,fo),
     >              f_ice(f(j),t2,j)   
 9430          format(f10.6,1x,f14.8,1x,f14.8,1x,f14.8)          
            enddo
         endif	   
      endif
      
      
      return
      end
      

************************************************************
*   this evaluates and writes zonal quantities             *
*     for current season                                   *
************************************************************
      subroutine zonal_qties_current_season(i, is, Tice, Tvapor, f, 
     >     fcmese, nlastorbits, t2, zonalT,
     >     fo, fc, Htime, Harea, phiDEG, tts)

      implicit none
      integer i, is, nlastorbits,k
      include 'parEBM.h'
      include 'functions.h'

      real*8 tempmat(Ns,N)  ! temperature matrix 
      common /tempmatrix/ tempmat

      real*8 f(N)
      real*8 phiDEG(N)          ! latitude in degrees 
      real*8 fcmese(Ns), fcTOT, fc(N)
      real*8 Htime(N), Harea(Ns)
      real*8 tts(Ns)
      real*8 fo(N) ! ocean coverage
      real*8 Tice, Tvapor, t2

      real*8 zonalT(N)

      integer j, nmese
      real*8 dfi, lathab, nonhab, fi1, fi2, LongS
      character*256 bfr, fname

*        -----------------------------------------------------------------
*        START CALCULATE AREA WEIGHTED GLOBAL QUANTITIES FOR CURRENT SEASON
*        ------------------------------------------------------------------
 	  
        fcmese(is)=0.	 
        do j=1,N                ! LOOP ON LATITUDES  
           dfi=dlat(j)  
           call subfc(f(j),t2,j,fo(j),fc(j))  
           fcmese(is)=fcmese(is)+dfi*fc(j)/2.  
        end do                  ! END LOOP ON LATITUDES  
        
*        ----------------------------------------------------------------
*        END CALCULATE AREA WEIGHTED GLOBAL QUANTITIES FOR CURRENT SEASON
*        ----------------------------------------------------------------

          if(nlastorbits.eq.5) then ! ONLY THE LAST ORBIT 
             do j=1,N
                if( f(j) .ge. Tice .and. f(j) .le. Tvapor )then
                   Htime(j)=Htime(j)+1./Ns
                endif
             enddo
          endif  
          
          Harea(is)=0.d0
          lathab=0.
          nonhab=0.
          
          if(nlastorbits.eq.5) then ! ONLY THE LAST ORBIT  
             do j=1,N
                if(xmin < -1.5) then 
                   fi1 = dsin(xmin+dx*(j-1)) !integration in latitude
                   fi2 = dsin(xmin+dx*j)
                else
                   fi1 = xmin+dx*(j-1) !integration in sin(latitude)
                   fi2 = xmin+dx*j
                endif
                dfi=fi2-fi1
                
                if( f(j) .ge. Tice .and. f(j) .le. Tvapor )then
                   Harea(is) = Harea(is) + dfi/2.0
                endif
                
                if( Htime(j) .gt. 0.999) then
                   lathab=lathab+dfi/2. !  1./N
                endif   
                
                if( Htime(j) .eq. 0.) then
                   nonhab=nonhab+dfi/2.	  
                endif   
                
             enddo
          endif 
          
          if(nlastorbits.eq.5) then ! ONLY THE LAST ORBIT  
             call ZonalMean(tempmat,zonalT)   !!perche' dentro il ciclo delle stagioni??? cosi' lo fa 48 volte!
          endif 
          
          if(nlastorbits.eq.5) then ! ONLY THE LAST ORBIT   
             write(25,33) (i+is/dfloat(Ns)), Harea(is)
          endif 

          if(nlastorbits.eq.5) then ! ONLY THE LAST ORBIT  
             nmese=is*(360./dfloat(Ns))
             bfr='Risultati/ZonalTemp'
             write(fname,'(a19,i0.3,a4)') bfr,nmese,".txt"
             
             open(unit=24,file=fname)
             do j=1,N
                write(24,32) phiDEG(j),f(j) 
             enddo
             close(24)
	 
c      store stellar longitude and mean global monthly temperature
c      WARNING: stellar longitude is correct only for circular orbits 
             LongS=is*(360./Ns)	 
             if ( LongS .gt. 310. ) then
                LongS=LongS-360.
             endif
             write(29,37) LongS,(tts(is)-273.15)
             
 37          format(f5.0,2x,f7.2)
          endif        


 32       format(f8.3,2x,6(f12.6,2x),2(1pe12.5,2x))
 33       format(f12.6,2x,f12.6)


      return
      end


************************************************************
*   this calculates and updates mean annual quantities     *
*       for the current orbit                              *
************************************************************
      subroutine calculate_annual_means(i, Tice, Tvapor, phiDEG, 
     >     annualglobalT, sigmaCRIT, fctot, printfhab, printchab,
     >     hcxl, TotOLR, TotASR, sigmaRG, sigmaBoil, DelT_NE, DTbc,
     >     meanTNorthHem, fhabN, meanAlbedoNorth, TotOLRN,
     >     T1bc, Sbc, vTgrad, nhab, fhab, chab, DelT_EP)

      implicit none
      include 'parEBM.h'
      include 'functions.h'
      common /tempmatrix/ tempmat

      include 'matrices.h'
      common/matrices/f, fo,  fcmat, olrmat, albmat, asrmat,
     >     phi, dryfluxmat, atmfluxmat, tts, habmat, cxlmat, 
     >     boilmat, RGmat
      
      integer i, is, j
      real*8 phiDEG(N)          ! latitude in degrees 
      real*8 fcmese(Ns), fcTOT, fc(N)
      real*8 Htime(N), Harea(Ns)
      real*8 Tice, Tvapor
      real*8 zonalT(N)


      real*8 sumcoslat, annualglobalT, fhab
      real*8 printfhab, printchab, printnhab
      real*8 TotOLR, TotASR, TotOLRAlb
      real*8 meanTNorthHem,TotOLRN, fhabN
      real*8 habtime(N), chab, nhab, sigmaBoil, sigmaRG
      real*8 fi1, fi2, dfi, ftme, ftmw
      real*8 meanTequat, wbc, meanTNorthPole, meanTSouthPole
      real*8 T1bc, T2bc, DTbc, sigmaCRIT,hcxl, meanAlbedoNorth
      real*8 DelT_NE, DelT_SE, DelT_EP, vTgrad, Sbc


      integer k, ks, kl

c       sum of the weights of planet area
       sumcoslat=0.
       do j=1,N
          sumcoslat=sumcoslat+dcos(phi(j)) 
       end do 

c        area-weighted global mean
       call GlobalMean(tempmat,phi,annualglobalT) ! temperature
       call GlobalMean(habmat,phi,fhab) ! habitability 
       printfhab=fhab    !!!!! PERCHE??
       call GlobalMean(cxlmat,phi,hcxl) ! complex life habitability
       call GlobalMean(olrmat,phi,TotOLR) ! OLR 
       call GlobalMean(asrmat,phi,TotASR) ! OLR 
       call GlobalMean(fcmat,phi,fcTot) 

c        area-weight global mean of the Northern hemisphere 
       call HemisphereMean(tempmat,phi,meanTNorthHem) 
       call HemisphereMean(habmat,phi,fhabN) 
       call HemisphereMean(olrmat,phi,TotOLRN)

c        ftime in Eq. (5)      
       do j=1,N
          habtime(j)=0.
          do is=1,Ns
             habtime(j)=habtime(j)+habmat(is,j)/Ns
          end do
       end do
       
      
       chab=0.                  !  CONTINUOUSLY HABITABLE LATITUDES    [ Eq. (8) in paper ]
       nhab=0.                  !  CONTINUOSLY NON-HABITABLE LATITUDES
       do j=1,N
          if(habtime(j).gt.0.999) then
             chab=chab+(dcos(phi(j))/sumcoslat)
          end if
          if(habtime(j).eq.0.) then
             nhab=nhab+(dcos(phi(j))/sumcoslat)
          end if
       end do
       printchab=chab !!!PERCHE?
       printnhab=nhab 
       
       
       sigmaBoil=0.0
       sigmaRG=0.0 
       do j=1,N
          
          if(xmin < -1.5) then 
             fi1 = dsin(xmin+dx*(j-1)) !integration in latitude
             fi2 = dsin(xmin+dx*j)
          else
             fi1 = xmin+dx*(j-1) !integration in sin(latitude)
             fi2 = xmin+dx*j
          endif
          dfi=fi2-fi1
          
          ftme=0.0
          ftmw=0.0
          do is=1,Ns
*     find out if this latitude zone is continuosly over the evaporation limit
             ftme=ftme+boilmat(is,j)
             ftmw=ftmw+RGmat(is,j) 
          end do
*     if zone is continuously hot, add it %area to sigmaBoil
          if( dabs(ftme-1.0*Ns).le.0.001) then
             sigmaBoil = sigmaBoil + dfi/2
          endif
          if( dabs(ftmw-1.0*Ns).le.0.001) then
             sigmaRG = sigmaRG + dfi/2
          endif
          
       end do
       
       sigmaCRIT=dmax1(sigmaBoil,sigmaRG) 
       
       meanTequat=0.
       wbc=0.
       do ks=1,Ns
          do kl=1,N
             if(phiDEG(kl).ge.-5..and.phiDEG(kl).le.+5.) then
                meanTequat=meanTequat+tempmat(ks,kl)*dcos(phi(kl)) 
                wbc=wbc+dcos(phi(kl))
             end if
          end do
       end do  
       meanTequat=meanTequat/wbc 
       
       meanTNorthPole=0.
       wbc=0.
       do ks=1,Ns 
          do kl=1,N
             if(phiDEG(kl).ge.85.) then
                meanTNorthPole=meanTNorthPole+tempmat(ks,kl)*
     >               dcos(phi(kl))
                wbc=wbc+dcos(phi(kl))
             end if
          end do 
       end do  
       meanTNorthPole=meanTNorthPole/wbc 
       
       meanTSouthPole=0.
       wbc=0.
       do ks=1,Ns 
          do kl=1,N
             if(phiDEG(kl).le.-85.) then
                meanTSouthPole=meanTSouthPole+tempmat(ks,kl)*
     >               dcos(phi(kl))
                wbc=wbc+dcos(phi(kl))
             end if
          end do 
       end do  
       meanTSouthPole=meanTSouthPole/wbc 
       
       T1bc=0.
       wbc=0.
       do ks=1,Ns
          do kl=1,N
             if(phiDEG(kl).ge.latbc1-5.and.
     >            phiDEG(kl).le.latbc1+5.) then
                T1bc=T1bc+tempmat(ks,kl)*dcos(phi(kl)) 
                wbc=wbc+dcos(phi(kl))
             end if
          end do
       end do  
       T1bc=T1bc/wbc 
       
       T2bc=0.
       wbc=0.
       do ks=1,Ns
          do kl=1,N
             if(phiDEG(kl).ge.latbc2-5.and.
     >            phiDEG(kl).le.latbc2+5.) then
                T2bc=T2bc+tempmat(ks,kl)*dcos(phi(kl)) 
                wbc=wbc+dcos(phi(kl))
             end if
          end do
       end do  
       T2bc=T2bc/wbc  
       
       DTbc=T1bc-T2bc
*     DTbc=(meanTNorthPole-meanTequat)/(pi2/4)  
*     if(DTbc.ge.0.) then
       if(DTbc.le.0.) then
          if(i.lt.10) DTbc=DTbcE
          if(i.ge.10) then
             print *,'delta T < 0'
             print *,'TRY TO REDUCE THE PLANET OBLIQUITY'
             stop
          end if
       end if   
       
       DelT_NE=dabs(meanTequat-meanTNorthPole)
       DelT_SE=dabs(meanTequat-meanTSouthPole)
       DelT_EP=0.5*(DelT_NE+DelT_SE)
       
       vTgrad=(P_H2O(T1bc)-P_H2O(T2bc))/DTbc
*     vTgrad=(P_H2O(meanTequat)-P_H2O(meanTNorthPole))/DelT_NE
       
*     Tbc=0.
       Sbc=0.
       wbc=0.
       do ks=1,Ns 
          do kl=1,N
             if(phiDEG(kl).ge.latbc1.and.phiDEG(kl).le.latbc2) then
*     Tbc=Tbc+tempmat(ks,kl)*dcos(phi(kl))
                Sbc=Sbc+asrmat(ks,kl)*dcos(phi(kl))
                wbc=wbc+dcos(phi(kl))
             end if 
          end do
       end do   
*     Tbc=Tbc/wbc
       Sbc=Sbc/wbc 


      return
      end
          


************************************************************
*   this calculates the zonal means                        *
*       for the current orbit                              *
************************************************************
      subroutine calculate_zonal_means(zonalALB, 
     >     zonalOLR, zonalASR, zonalATFX, zonalDRFX,
     >     zonalfc, annualglobalA, meanAlbedoNorth,
     >     Tmin, Tmax)


      implicit none
      include 'parEBM.h'
      include 'functions.h'
      common /tempmatrix/ tempmat

      include 'matrices.h'
      common/matrices/f, fo,  fcmat, olrmat, albmat, asrmat,
     >     phi, dryfluxmat, atmfluxmat, tts, habmat, cxlmat, 
     >     boilmat, RGmat
      
      integer i, ks, j, is, kl
      real*8 phiDEG(N)          ! latitude in degrees 
      real*8 fcmese(Ns), fcTOT, fc(N)
      real*8 Htime(N), Harea(Ns)

      real*8 zonalT(N)
      real*8 zonalALB(N)
      real*8 Tmin, Tmax
      real*8 annualglobalA, wbc, meanAlbedoNorth
      real*8 zonalATFX(N), zonalDRFX(N), zonalOLF(N), zonalASR(N)
      real*8 zonalOLR(N), zonalfc(N)

       call ZonalMean(fcmat,zonalfc) ! cloud coverage
       call ZonalMean(atmfluxmat,zonalATFX) ! horizontal atmopsheric flux
       call ZonalMean(dryfluxmat,zonalDRFX) ! horizontal atmopsheric flux
       call ZonalMean(olrmat,zonalOLR) ! OLR
       call ZonalMean(asrmat,zonalASR) ! Absorbed Stellar Radiation S(1-A)
       
c     mean annual zonal ALBEDO      
       do j=1,N  
          zonalALB(j)=0.
          ks=0
          do is=1,Ns
             if(albmat(is,j).gt.0.) then
                zonalALB(j)=zonalALB(j)+albmat(is,j)
                ks=ks+1
             end if
          end do
          zonalALB(j)=zonalALB(j)/ks
       end do
       
       Tmin=1.d4
       Tmax=0.
       do ks=1,Ns
          do kl=1,N
             Tmin=dmin1(Tmin,tempmat(ks,kl))
             Tmax=dmax1(Tmax,tempmat(ks,kl))
          end do
       end do
       
c     mean annual global albedo
       annualglobalA=0.
       wbc=0.
       do ks=1,Ns
          do kl=1,N
             if(albmat(ks,kl).gt.0.) then
                annualglobalA=annualglobalA+albmat(ks,kl)*dcos(phi(kl)) 
                wbc=wbc+dcos(phi(kl))
             end if
          end do
       end do
       annualglobalA=annualglobalA/wbc  
       
c     mean annual albedo northern hemisphere
       meanAlbedoNorth=0.
       wbc=0.
       do ks=1,Ns
          do kl=1,N
             if(phi(kl).ge.0.) then
                if(albmat(ks,kl).gt.0.) then
                   meanAlbedoNorth=meanAlbedoNorth+albmat(ks,kl)*
     >                  cos(phi(kl)) 
                   wbc=wbc+dcos(phi(kl))
                end if
             end if
          end do
       end do
       meanAlbedoNorth=meanAlbedoNorth/wbc  


      return
      end


************************************************************
*   this calculates applies convergence criteria           *
************************************************************
      subroutine calculate_convergence(i, nprompt, 
     >     convergence, annualglobalT,tsum, tsum_old,
     >     annualglobalA, fcTOT, iceTOT, pressPtot, sigmaCRIT, 
     >     sigmaLIM, Tmax, f, tts, t2, exitFLAG)

      implicit none
      include 'parEBM.h'


      integer i, nprompt, convergence, is
      real*8 f(N), t2
      real*8 tts(Ns)
      real*8 annualglobalT, annualglobalA, fcTOT
      real*8 pressPtot, sigmaCRIT, sigmaLIM, Tmax
      real*8 tsum, tsum_old, ICEtot
      real*8 exitFLAG
      
      real*8 icecover

      if( (i/nprompt)*nprompt .eq. i) then  
         tsum=0.  
         do is=1,Ns
            tsum=tsum+tts(is)/Ns  
         enddo
         
         iceTOT=icecover(f,t2) 
         
         write(*,702) i,annualglobalT,annualglobalT-273.15,
     >        annualglobalA,fcTOT,
     >        iceTOT,pressPtot
         
         if(i.ge.30) then       ! search for convergence only after 30 orbits
                                ! (20 with ice of WK97 and 10 with our ice)
            if(dabs(tsum-tsum_old).lt.deltaTconv) then
               if(sigmaCRIT.le.sigmaLIM.and.Tmax.ge.Tlim1) then
                  write(*,953) tsum,tsum_old
 953              format('Simulation converged: <T>=',f8.4,
     >                 '  <T>old=',f8.4) 	   
                  convergence=1  
               endif
            endif 
         endif 
         
         tsum_old=tsum 
      endif
      
 702  format('orbit',i4,2x,'<T>=',f9.4,'K ',f9.4,'C',
     >       3x,'<A>=',f6.4,' <clouds>=',f6.4,' <ice>=',f6.4,
     >       '  pTOT=',1pe10.4,'Pa') 

* ---------------------------------------------------------------------------------------
*     if convergence has not been achieved and the simulation has run at least two orbits,
*     test condition of FORCED EXIT 
* ---------------------------------------------------------------------------------------
      if(convergence.eq.0.and.i.ge.2) then

c       total pressure exceeds maximum allowed value (here 10 bars)      
         if(pressPtot.gt.10.d5) then  
            write(*,2954) pressPtot,annualglobalT 
            exitFLAG=-2.
            return
*     go to 955 using exitFLAG for this
         end if

 2954    format('Simulation interrupted: ',
     >        'Total pressure ',1pe10.4,' out of allowed range',
     >        3x,'<T>=',0pf9.4,'K')  

c        next forced criteria are applied only after 30 orbits
c        in this way we give time to the ice cover routine to make permanent ices, if any             
         if(i.ge.50) then  
            if(sigmaCRIT.gt.sigmaLIM) then
               exitFLAG=-1.
               write(*,9541) annualglobalT,sigmaCRIT
               return
*     go to 955        using exitFLAG for this
            endif
            
c       mean planet temperature lower than half the minimum value  
            if(Tmax.lt.Tlim1/2.) then  
               exitFLAG=-0.5
               write(*,9542) Tmax,Tlim1 
               return
*     go to 955 using exitFLAG for this
            endif 

         end if                 ! end if(i.ge.30)  
       
 9541    format('Simulation interrupted: <T>=',f9.4,
     >        ' sigmaCRIT=',f8.6)  
 9542    format('Simulation interrupted: Tmax=',f6.1,
     >             'K < ',f6.1,'K')       
      end if                    ! end  if(convergence.eq.0)
           
            
      return
      end


************************************************************
*   final output and closing of files                     *
************************************************************
      subroutine make_output(i)

      implicit none
      include 'parEBM.h'
      
      common /tempmatrix/ tempmat

      include 'matrices.h'
      common/matrices/f, fo,  fcmat, olrmat, albmat, asrmat,
     >     phi, dryfluxmat, atmfluxmat, tts, habmat, cxlmat, 
     >     boilmat, RGmat

      common/output1/Tmin, Tmax,
     >     printfhab, printchab, hcxl, TotOLR, sigmaRG,
     >     sigmaBoil, DelT_NE, meanTNorthHem, fhabN,
     >     meanAlbedoNorth, TotOLRN,
     >     fhab, chab,
     >     fom, t2
      common/output2/
     >     annualglobalT, DelT_EP, nhab,
     >     annualglobalA, TotASR, exitFLAG,
     >     fcTOT, iceTOT, phiDEG, zonalT, Htime, zonalOLR,
     >     zonalASR, zonalALB, zonalATFX, zonalDRFX

      common /gasplusvapor_pars/ pressPtot,cp_Ptot,molwtPtot

      common /transportpar/ DTbc,T1bc,Sbc,vTgrad 


      integer i, j
      real*8 Tmin, Tmax
      real*8  printfhab, printchab, hcxl, TotOLR, sigmaRG
      real*8     sigmaBoil, DelT_NE, DTbc, meanTNorthHem, fhabN
      real*8     meanAlbedoNorth, TotOLRN, T1bc, Sbc
      real*8     vTgrad, fhab, chab, fom, t2
      real*8     annualglobalT, DelT_EP, nhab
      real*8     annualglobalA, TotASR, exitFLAG
      real*8     fcTOT, iceTOT
      real*8 phiDEG(N), zonalT(N), Htime(N), zonalOLR(N)
      real*8 zonalASR(N), zonalALB(N), zonalATFX(N), zonalDRFX(N)
      real*8 f_ice
      real*8 pressPtot, cp_Ptot, molwtPtot, Tice, Tvapor

      integer date_time(8)
      character*10 bb(3)
      character*256 bfr


      write(*,1955) i,Tmin,Tmax
1955  format('Exit simulation after ',i3,' orbits',2x,
     >       'Tmin= ',f5.1,2x,'Tmax= ',f5.1)

      close(19) 
      close(23) 
      close(25)
      close(27)
      close(29)
      close(75) 
      
      if(albedoFile.eq.'none') close(40)  
      
      write(*,701) printfhab,printchab,hcxl
 701  format('h=',f6.4,3x,'chl=',f6.4,3x,
     >     'h(complex life)=',f6.4)

      write(*,7016)  TotOLR,sigmaRG,sigmaBoil
 7016 format('TotOLR=',f6.1,'W/m2',3x,
     >     'sigmaRG=',f8.6,3x,'sigmaBoil=',f8.6)   

      write(*,7012) DelT_NE,DTbc
 7012 format('<DT>ENP=',f6.2,'K   DTbc=',f8.4,'K')
      
      write(*,7015) meanTNorthHem,fhabN,meanAlbedoNorth,TotOLRN
 7015 format('<T>north=',f8.4,3x,'<h>north=',f6.4,3x,
     >     '<A>north=',f8.4,3x,'<OLR>north=',f6.1,'W/m2') 
      
      write(*,7013) T1bc,Sbc
 7013 format('baroclinic zone: <T1bc>north=',f8.4,
     >     'K  <Sbc>=',1pe12.5,'W/m2')

      write(*,7014) vTgrad,vTgrad/vTgradE
 7014 format('vTgrad=',1pe13.6,2x,'vTgrad/vTgradE=',0pf6.3) 
      
c     write mean annual zonal values: latitude, temperature, ftime Eq.(5), OLR 
      open(unit=216,file='Risultati/ZonalData.txt',status='unknown')
      do j=1,N
         write(216,32) phiDEG(j),zonalT(j),Htime(j),zonalOLR(j),
     >        zonalASR(j),zonalALB(j),f_ice(zonalT(j),t2,j),
     >        zonalATFX(j),zonalDRFX(j)
      enddo
 32   format(f8.3,2x,6(f12.6,2x),2(1pe12.5,2x))
      close(216)   
      close(55)
      close(56)
      close(57)
      close(58)
      close(59) 
      
      close(77)                 ! warnings file 
      
      open(unit=26,file='Risultati/label.txt',status='unknown') 
      write(26,34) "a=",smaP,"e=",eccP,"h=",fhab,"hl=",chab     
 34   format(a2,f4.2,1x,a2,f5.3,1x,a2,f5.3,1x,a3,f5.3)             
      close(26)
      
      open(unit=61,file='labeltab.txt',status='unknown') 
      write(61,35) LumStar/LumSun,omegaPERI,obliq,Prot,fom
 35   format(f7.3,2x,f6.1,2x,f5.1,2x,f4.0,2x,f4.2)
      close(61) 


* class of the solution (only for CONVERGED sims):
      call LWTR(pressPtot,Tice,Tvapor) 
      if (exitFlag .ge. 0.0) then
         if (annualGlobalT .gt. Tice .and. Tmax .lt. Tvapor ) then
            exitFlag=1.0        !WARM
         else if (Tmax .ge. Tvapor ) then
            exitFlag=2.0        !WARM_HOT
         else if (iceTOT .gt. 0.99) then
            exitFlag=3.0        !SNOWBALL
         else if (annualGlobalT .le. Tice .and. iceTOT .le. 0.99) then
            exitFlag=4.0        !WATERBELT
         else
            exitFlag=-200.0     !UNDEFINED
         endif
      endif
      
c     write a summary of input and output parameters 
      open(unit=28,file='Risultati/valori.txt',status='unknown')
      write(28,39) Mstar,LumStar,smaP,eccP,omegaPERI,obliq,
     >     Prot,fom,pressP,q0,Porb/86400,annualglobalT,DelT_EP,
     >     fhab,chab,nhab,annualglobalA,float(i),Tmin,Tmax,
     >     asl,Rpar,TotOLR,sigmaRG,sigmaBoil,exitFlag
      close(28)
 39   format(26(1p e12.5,1x))  
      
      open(unit=29,file='Risultati/GlobalData.txt',status='unknown')
      write(29,*) '# planet    <T>    <A>    <OLR>    <ASR>   ',
     >     '<clouds>    <ice>    <h>    <hcxl>    ',
     >     '(DT)ep    Tmin    Tmax    dTdphi    dqdT '
      write(29,395) planet,annualglobalT,annualglobalA,TotOLR,TotASR,
     >     fcTOT,iceTOT,fhab,hcxl,
     >     DelT_EP,Tmin,Tmax,DTbc,vTgrad
 395  format(1x,a16,26(1x,1pe12.5)) 
      close(29)



c  This produces the file "esopianeti.par", needed by fits_map_temperature.py to create
c  the .fits file from the run. The above python script also need year_lat_temp_last1.tlt

      open(unit=28,file="Risultati/esopianeti.par",status='unknown')
      write(28,'("NAME!",A,"!planet name!STR")') planet
      write(28,'("MSTAR!",E16.8,"!solar masses!F")') mstar
      write(28,'("LUMSTAR!",E16.8,"!luminosity in solar units!F")') 
     >     Lumstar
      write(28,'("SMA!",F16.8,"!semi-major axis [AU]!F")') smaP
      write(28,'("ECC!",F16.8,"!eccentricity!F")') eccP
      write(28,'("OMEGAPER!",F16.8,"!argument of pericenter!F")') 
     >omegaPERI
      write(28,'("GRAV!",F16.8,
     >"!surface gravity [gEARTH=9.8 m/s**2]!F")') gP
      write(28,'("OBLIQ!",F16.8,
     >"!obliquity of planet rotation axis [deg]!F")') obliq
      write(28,'("PROT!",F16.8,"!planet rotation period [days]!F")') 
     >     Prot
      write(28,'("MPLAN!",E16.8,
     >"!planet mass [Earth masses]!F")') Pmass
      write(28,'("RPLAN!",E16.8,
     >"!planet radius [R(earth)=6.371e6 m]!F")') Rplanet
      bfr='("GEO!",I2,
     >"!planet geography (ocean fract. in lat zones)!I")'
      write(28,bfr) gg
      write(28,'("FO_CONST!",F16.8,"!const fraction of oceans!F")'), 
     >     fo_const
      write(28,'("PRESS!",E16.8 
     >"!total dry pressure at planet surface [Pa]!F")') pressP
      write(28,'("PO2!209460.0!O2 partial pressure [Pa]!F")')    !!!!WARNING, FOR FUTURE USE
      write(28,'("PN2!780840.0!N2 partial pressure [Pa]!F")')    !!!!WARNING, FOR FUTURE USE
      write(28,'("P_CO2!",E16.8,"!CO2 partial pressure [ppvm]!F")') 
     >     p_CO2_P*10 !!!WARNING, from Pascal to PPMV
      write(28,'("P_CH4!1.8!CH4 partial pressure [ppvm]!F")') !!!!WARNING, FOR FUTURE USE
      write(28,'("P_O3!0.0!PO3 partial pressure [ppvm]!F")')  !!!!WARNING, FOR FUTURE USE
      write(28,'("RH!",F16.8,"!relative humidity!F")') RH
      write(28,'("TMGLOB!",E16.8,"!mean orbital globaltemperature!F")'), 
     > annualGlobalT
      write(28,'("HLW!",F16.8,
     >"!mean orbital global liquid-water habitability!F")') fhab
      write(28,'("H050!",F16.8,
     >"!mean orbital globalcomplex-life habitability!F")') hcxl
      write(28,'("ALB!",F16.8,
     >"!mean orbital albedo!F")') annualglobalA
      write(28,'("CLOUDS!",F16.8,
     >"!mean orbital cloud covering!F")') fctot
      write(28,'("NORBIT!",I3,"!number of orbit before convergence!I")')
     >     , i

      write(28,'("CONTHAB!",F16.8,
     >"!index of continuous liquid-water habitability!F")') chab
      write(28,'("DTEP!",F16.8,
     >"!equator-pole temperature difference!F")') DelT_EP
      write(28,'("MOLR!",F16.8,
     >"!mean OLR!F")') TotOLR
      write(28,'("MARS!",F16.8,
     >"!mean ASR !F")') TotASR
      write(28,'("ICE!",F16.8,
     >"!mean global ice coverage!F")') iceTOT

* class of the solution:
      if (annualGlobalT .gt. Tice .and. Tmax .lt. Tvapor ) then
         write(28,'("CLASS! WARM!class of the solution!STR")') 
      else if (Tmax .ge. Tvapor ) then
         write(28,'("CLASS! WARM_HOT!class of the solution!STR")') 
      else if (iceTOT .gt. 0.99) then
         write(28,'("CLASS! SNOWBALL!class of the solution!STR")') 
      else if (annualGlobalT .le. Tice .and. iceTOT .le. 0.99) then
         write(28,'("CLASS! WATERBELT!class of the solution!STR")') 
      else
         write(28,'("CLASS! UNDEFINED!class of the solution!STR")') 
      endif


      call date_and_time(bb(1),bb(2),bb(3),date_time)
      write(28,'("DATE!",I0.2,I0.2,I4,"!date of the run!STR")') , 
     >     date_time(3),date_time(2),date_time(1)
      write(28,'("NUMBER!",I0.4,
     >"!Ordered map number for the above date!I")') NUMBER   
      write(28,'("VERSION!",A,"!code version!STR")')  VERSION 
      write(28,'("SIMTYPE!",A,"!type of run!STR")') SIMTYPE 
      write(28,
     >'("PAPER1!Vladilo+2015,ApJ,804,50!Reference paper (model)!STR")')
      bfr='("PAPER2!Silva+2016,Int.J.Astrobiology,pp.1-22 doi:10.1017/S1
     >473550416000215"!Ref. paper (habitability)!STR")'
      write(28,bfr)      
      write(28,'("PRJNAME!EXOCLIMATES!project name!STR")')
      write(28,'("COMMENT!! The planet geography expressed in 
     >terms of")')
      write(28,'("COMMENT!! ocean fractions in each latitude zone")')
      write(28,'("COMMENT!! 0: constant freaction of oceans")')
      write(28,'("COMMENT!! 1: present Earth WK97; 2: equatorial cont.WK97; 
     >3: polar cont.WK97")')
      write(28,'("COMMENT!! 4: present Earth, file fo_earth_DMAP.dat")')
      
      close(28)



 98   print *,"finito!"
      
      return
      end
      
