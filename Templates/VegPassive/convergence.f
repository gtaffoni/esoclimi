
************************************************************
*   this calculates applies convergence criteria           *
************************************************************
      subroutine calculate_convergence(i, nprompt, nconv,
     >     convergence, annualglobalT,tsum, tsum_old, ddeltaTold,
     >     annualGlobalVeg, annualGlobalVeg_old,
     >     annualglobalA, fcTOT, iceTOT, pressPtot, sigmaCRIT, 
     >     sigmaLIM, Tmax, f, tts, t2, exitFLAG)

      implicit none
      include 'parEBM.h'
      include 'vegetation.h' !this contains annualglobalVeg[_old] and deltaVconv

      integer i, nprompt, convergence, is, nconv
      real*8 f(N), t2
      real*8 tts(Ns)
      real*8 annualglobalT, annualglobalA, fcTOT
      real*8 pressPtot, sigmaCRIT, sigmaLIM, Tmax
      real*8 tsum, tsum_old, ICEtot, ddeltaTold
      real*8 exitFLAG
      
      real*8 icecover

      tsum=0.  
      do is=1,Ns
         tsum=tsum+tts(is)/Ns  
      enddo
         
      iceTOT=icecover(f,t2) 

      if( (i/nprompt)*nprompt .eq. i) then  
         
         write(*,702) i,annualglobalT,annualglobalT-273.15,
     >        annualglobalA,fcTOT,
     >        iceTOT,pressPtot
         
      endif

      if(i.ge.30) then          ! search for convergence only after 30 orbits
                                ! (20 with ice of WK97 and 10 with our ice)
         if(dabs(tsum-tsum_old)/tsum_old.lt.deltaTconv) then

            if(sigmaCRIT.le.sigmaLIM.and.Tmax.ge.Tlim1) then

               if(convergence.eq.0) then
                  convergence=1 !this means: check next npropt orbits wth the SAME tsum_old
                  
               else if(convergence.eq.1.and.nconv.lt.nprompt) then   

* rimettere anche il check sulla vegetazione                  
                  write(*,9553) nconv,dabs((tsum-tsum_old)/tsum_old),
     >                 deltaTconv, dabs(tsum-tsum_old),
     >                 deltaTconv*tsum_old

                  if(nconv.eq.0) then
                     ddeltaTold = dabs(tsum-tsum_old)/tsum_old
                  else

                     if( dabs( ( (tsum-tsum_old)/tsum_old - ddeltaTold))
     >                .gt. ddeltaTconv) then
                         write(*,9554) dabs( (tsum-tsum_old)/tsum_old -
     >                      ddeltaTold), ddeltaTconv
 9554                    format(
     >         '   Convergence NOT acheived: variation in dT/T is ',
     >              1pe10.4, ' wanted: ', 1pe10.4)
                         print *, dabs( (tsum-tsum_old)/tsum_old -
     >                      ddeltaTold), ddeltaTconv, ddeltaTold
                         nconv=-1
                         convergence=0
                      endif
                     
                   endif
                   
                  
                  nconv = nconv+1
               else if(convergence.eq.1.and.nconv.eq.nprompt) then
                  convergence=2                                   
                  write(*,953) tsum,tsum_old
               endif

 9553          format('   Checking achieved convergence, orbit ',i4,
     >              ' deltaT/T=',1pe10.4, ' wanted: ',1pe10.4,
     >          '(deltaT:',
     >              f8.4, ' wanted: ',f8.4,')')
 953           format('Simulation converged: <T>=',f8.4,
     >              '  <T>old=',f8.4) 	   
            endif               ! sigmaCRIT, Tmax
         else                   !converce condiction failed. Reset convergence and timer
            convergence=0
            nconv=0
         endif ! dabs*tsum-tsum_old).le.deltaTconv

      endif 
         
      if(convergence.eq.0) then
         tsum_old=tsum 
         annualglobalVeg_old =  annualglobalVeg
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
            if(Tmax.lt.Tlim1) then  
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

