c auxiliary file containing the definitions of matrices used
c in common/matrices/, used for passing data to subroutines


      real*8 f(N)
      real*8 tempmat(Ns,N)  ! temperature matrix 
      real*8 fo(N) ! ocean coverage
      real*8 fcmat(Ns,N)    ! cloud cover matrix
      real*8 cxlmat(Ns,N)   ! complex life habitability
      real*8 olrmat(Ns,N),asrmat(Ns,N)  ! OLR and ASR matrix 
      real*8 albmat(Ns,N)   ! TOA albedo matrix
      real*8 phi(N) ! latitude in radians 
      real*8 atmfluxmat(Ns,N),dryfluxmat(Ns,N) ! meridional atm.  flux
      real*8 boilmat(Ns,N),RGmat(Ns,N)                   
      real*8 habmat(Ns,N)   ! habitability function matrix
      real*8 tts(Ns)



