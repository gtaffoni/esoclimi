c      write(*,'(a)') TOAalbFile
  
      OPEN(55,file=TOAalbFile,status='old') !TOAalbFile
c       skip first 4 lines with arrays       
      read(55,8120) MT,(T_TOA(iT),iT=1,NT_TOA) 
      read(55,8121) Mp,(p_TOA(iT),iT=1,Np_TOA) 
      read(55,8122) Mz,(z_TOA(iT),iT=1,Nz_TOA)  
      read(55,8123) Mas,(as_TOA(iT),iT=1,Nas_TOA) 

 8120 format(i2,1x,19f6.1) 
 8121 format(i2,1x,16f9.5)
 8122 format(i2,1x,21f6.2)
 8123 format(i2,1x,21f6.3) 
      
       DO ias=1,Nas_TOA         ! surface albedo
          DO iz=1,Nz_TOA        ! zenith distance
             DO ip=1,Np_TOA     ! pressure
                READ(55,8130) (Matrix_TOA(iT,ip,iz,ias),iT=1,NT_TOA) ! temperature
             END DO
          END DO
       END DO
       
 8130  format(21F7.3)
