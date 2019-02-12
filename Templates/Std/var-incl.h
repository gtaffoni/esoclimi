c general variable, vector and matrix definition

      real*8 df(N),fout(N)
      real*8 Htime(N), Harea(Ns)
      real*8 fyear, fhab, lathab, nonhab, hcxl
      real*8 time,time0,tsum,tsum_old,ddeltaT_old
      integer i,is,j,nmese,kl,ks
      real*8 LongS, tclx1, tclx2
      integer nok, nbad, nprompt, nout  
      real*8 t2 , t1 !, h1, hmin, eps
      real*8 fi1, fi2, dfi 
      real*8 x(N) 
      real*8 phiDEG(N) ! latitude in degrees 
      real*8 fclouds    ! cloud coverage
      real*8 fcmese(Ns), fcTOT, fc(N)
      real*8 iceTOT, fom
      real*8 sumcoslat,habtime(N)
      real*8 chab,nhab,fhabN
      real*8 printfhab,printchab,printnhab
      real*8 TotOLR, TotASR, TotOLRAlb

      real*8 DTbc,Sbc,wbc ! BAROCLINIC ZONE PARAMETERs FOR BARRY TRANSPORT 
      real*8 vTgrad            ! gradient FOR MOIST VERSION
      real*8 meanTpoles,meanTequat
      real*8 meanTNorthPole,meanTSouthPole
      real*8 DelT_NE,DelT_SE,DelT_EP
      real*8 meanTNorthHem,TotOLRN
      real*8 T1bc,T2bc
      real*8 Ddry

      real*8 Tice,Tvapor ! liquid water temperature range 
      real*8 Tmin,Tmax   ! minimum and maximum planet temperature
      real*8 pressPtot,cp_Ptot,molwtPtot
      real*8 annualglobalT,annualglobalA
      real*8 meanAlbedoNorth
      real*8 zonalfc(N),zonalOLR(N),zonalATFX(N),zonalDRFX(N)
      real*8 zonalT(N),zonalASR(N),zonalALB(N)
    
      real*8 ftme,ftmw,sigmaBoil,sigmaRG,sigmaCRIT
      real*8 sigmaLIM,exitFLAG,cmH2Olim,cmH2O
      integer convergence,nlastorbits,nconv

      integer nlatDMAP
      parameter(nlatDMAP=46)
      real*8 latDMAP(nlatDMAP),foDMAP(nlatDMAP) 

      INTEGER iT,ip,iz,ias  
      INTEGER MT,Mp,Mz,Mas 

      character*50 fname
      character*50 bfr  
 
      integer orbitdays
      real*8 mzd(N)  ! mean orbital zenith distance
      real*8 DMM ! function used to calculate the mean zenith distance
      real*8 OrbitMeanZD,HourD,HHS(N)   ! orbit mean zenith distance
