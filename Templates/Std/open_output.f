      fname='Risultati/init.dat'
      open(unit=20,file=fname,status='unknown')
      do i=1,N
      write(20,20) xmin+(i-1)*dx+dx/2,f(i),fo(i)
20    format(f8.5,1x,f14.8,1x,f14.8)
      enddo
      close(20) 
       
      open(unit=55,file='Risultati/Sterm_lat000.txt',status='unknown')
      open(unit=56,file='Risultati/Sterm_lat+45.txt',status='unknown')
      open(unit=57,file='Risultati/Sterm_lat-45.txt',status='unknown')
      open(unit=58,file='Risultati/Sterm_lat+90.txt',status='unknown')
      open(unit=59,file='Risultati/Sterm_lat-90.txt',status='unknown') 

      open(unit=77,file='Risultati/warnings.txt',status='unknown')
      open(unit=75,file='Risultati/cloudalbedo.txt',status='unknown')
 

c       MOST INPUT AND OUTPUT DATA ARE STORED IN THE FILE "valori.txt"
c       here we write an initial version of the file to be used when the simulation is interrupted 
c       because the stepsize is too small  
      open(unit=28,file='Risultati/valori.txt',status='unknown')
      write(28,39) Mstar,LumStar,smaP,eccP,omegaPERI,obliq,
     >  Prot,fom,pressP,q0,Porb/86400,-100.,-100.,
     >  -100.,-100.,-100.,-100.,-100.,-100.,-100.,
     >  asl,Rpar,-100.,-100.,-100.,-100.
      close(28)
39    format(26(1p e12.5,1x))  

      
      if(albedoFile.eq.'none') 
     >  open(unit=40,file='Risultati/year_lat_alb.tla',status='unknown')
     
      open(unit=19,file='Risultati/year_lat_temp.tlt',status='unknown')      
      open(unit=23,file='Risultati/year_lat_temp_last5.tlt',
     >              status='unknown')
      open(unit=27,file='Risultati/year_lat_temp_last1.tlt',
     >              status='unknown')
      
      if(albedoFile.eq.'none') write(40,31) Ns, N ! display albedo-latitude, last two orbits 
      
      write(19,31) Ns, N
      write(23,31) Ns, N
      write(27,31) Ns, N
31    format(i14,1x,i14)      

      open(unit=29,file='Risultati/longOrb_temp.ebm',status='unknown')
      open(unit=25,file='Risultati/year_farea.ebm',status='unknown')
33    format(f12.6,2x,f12.6)

