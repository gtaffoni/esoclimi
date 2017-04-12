      open(unit=24,file='Risultati/year_lat_veg_last5.tlt',
     >              status='unknown')
      open(unit=28,file='Risultati/year_lat_veg_last1.tlt',
     >              status='unknown')
      open(unit=30,file='Risultati/year_veg_average.dat',
     >              status='unknown')
      
      
      write(24,319) Ns, N
      write(28,319) Ns, N
 319  format(i14,1x,i14)      

