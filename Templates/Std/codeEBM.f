
      program codeEBM

      implicit none
      integer k
      include 'parEBM.h'

      include 'var-incl.h'
      include 'matrices.h'
      include 'functions.h'

c       MODULE INCLUDES
      include 'module_include.f'

      include 'commons.h'

c        SPLASH MESSAGE
      include 'splash_output.f'

c       INITIALIZE VARIABLES
      include 'initialize.f'
c       MODULE INITIALIZATION
      include 'module_initialize.f'


c       READ FILE WITH TOP-OF_ATMOSPHERE DATA
      include 'read_toa.f'

c        calculate LIQUID WATER TEMPERATURE RANGE at total pressure pressPtot
      call LWTR(pressPtot,Tice,Tvapor) 

c
c  STUDY OF BISTABILITY
c
c  OVERWRITING TEMPERATURE INITIALIZATION
c
c  Here I evaluate DeltaT in Nt=20 intervals and use the value CurrentT*DeltaT
c  THIS VALUE MUST BE CHANGED BY-HANDS FOR EACH RUN
c
c

      print *, '**************************'
      print *, 'Tstart:'
      print *,Tstart
      print *, Tice, Tvapor
      print *, '**************************'
c
c
c
 
      do i=1,N
         f(i)=Tstart            ! zonal temperatures
      enddo
      tsum_old=Tstart
*      ddeltaTold=0.0
      

c       PRINT INPUT MODEL PARAMETERS
      include 'print_input.f'
         

c        OPEN OUTPUT FILES
      include 'open_output.f'
c        MODULE OPEN OUTPUT FILES
      include 'module_open_output.f'


      nlastorbits=0      ! number of orbits after convergence or after forced exit  
      convergence=0      ! flag for convergence
      nconv=0            ! flag for continous convergence, nprompt orbits
      
*************************************************************************************
*     START LOOP ON ORBITS
*************************************************************************************
      do i=0, maxNorbits+5  

         if(convergence.eq.2.or.i.gt.maxNorbits) then
            nlastorbits=nlastorbits+1 
         endif

* --------------------------------------------------
*     START LOOP ON SEASONS   
* --------------------------------------------------
         do is=1,Ns 
          
            tts(is)=0.d0        ! initialize mean global temperature of current season
            t2 = t1 +  Porb/Ns 

c         Module integrator	  
            include 'integrator_call.f'

*         This updates all seasonal quantities
            call update_seasonal_qties(is, t2, pressPtot, 
     >           cp_Ptot, molwtPtot, tclx1, tclx2,
     >           cmH2O, cmH2Olim)
*         as above, for module-dependent quantities
            include 'module_update_season.f'


*         This makes outputs using the above quantities
            call write_seasons(i, is, nlastorbits, t2, f, fo, 
     >           phi, phiDEG)
*         as above, for module-dependent quantities
            include 'module_write_seasons.f'

            t1=t2

            call zonal_qties_current_season(i, is, Tice, Tvapor, f, 
     >           fcmese, nlastorbits, t2, zonalT,
     >           fo, fc, Htime, Harea, phiDEG, tts)
          
         end do                 !  END LOOP is=1,Ns  
* --------------------------------------------------
*     END LOOP ON SEASONS  - we are still inside the loops on orbits
* --------------------------------------------------  


c     calculate and update mean annual quantities for the current orbit
         call calculate_annual_means(i, Tice, Tvapor, phiDEG, 
     >        annualglobalT, sigmaCRIT, fctot, printfhab, printchab,
     >        hcxl, TotOLR, TotASR, sigmaRG, sigmaBoil, DelT_NE, DTbc,
     >        meanTNorthHem, fhabN, meanAlbedoNorth, TotOLRN,
     >        T1bc, Sbc, vTgrad, nhab, fhab, chab, delT_EP)
*     as above, for module-dependent variables
         include 'module_annual_mean.f'



c     mean annual zonal arrays
         call calculate_zonal_means(zonalALB, 
     >        zonalOLR, zonalASR, zonalATFX, zonalDRFX,
     >        zonalfc, annualglobalA, meanAlbedoNorth,
     >        Tmin, Tmax)
       
       
* -----------------------------------------------------------------------------
*     APPLY CONVERGENCE CRITERIA every nprompt orbits, starting after 30 orbits
* -----------------------------------------------------------------------------

*     This calls the correct calculate_convergence subroutine (different for
*         vegetation and no vegetation)
         include 'module_convergence.f'

         
*      write annual module output and splash output if needed
         include 'module_annual_output.f'

* simulation interrupted:
         if(exitFLAG.eq.-2.0.or.exitFLAG.eq.-1.0.or.exitFLAG.eq.-0.5) 
     >        go to 955
     
* ---------------------------------------------------
*     REGULAR EXIT
*     the regular exit is allowed 5 orbits after convergence has been achieved
*     in this way, the data of the last 5 orbits are stored for analysis of temporal evolution
* ---------------------------------------------------   
         if(nlastorbits.eq.5) go to 955     
      
      end do  

*************************************************************************************
*     END LOOP ON ORBITS
*************************************************************************************

c      data for output is in a common
      
 955  call make_output(i)



*       final output for module, and close output files of module
      include 'module_final_output.f'

      stop
      end

 
*************************************************************************************
*     END MAIN PROGRAM
************************************************************************************* 



