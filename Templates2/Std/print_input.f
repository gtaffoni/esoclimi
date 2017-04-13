      write(*,4962)  Mstar/Msun ,q0   
4962  format('STELLAR MASS=',f5.3,' M_sun',3x,
     >       'INSOLATION=',f7.1,' W/m2 at r=a')  

      write(*,496) planet,smaP,eccP,Porb/86400.
496   format('PLANET: ',a16,'a=',f6.4,' AU',3x,'e=',f6.4,
     >    3x,'ORBITAL PERIOD: ',f6.2,' terrestrial days')     

      write(*,498) Rplanet/Rearth,Prot,obliq 
498   format('RADIUS: ',f7.4,' R_earth',3x,
     >       'ROTATION PERIOD: ',f6.3,' d',3x,
     >       'AXIS OBLIQUITY: ',f7.3,' deg')  
     
      write(*,4982) D0par,Rpar,C0par,C1par
4982  format('DIFFUSION PARs.: D0=',f6.3,'W/m2',3x,
     >       'R=',f5.2,3x,'C0=',f8.5,3x,'C1=',f8.5)  

      write(*,4985) pressP,Tice,Tvapor,Tlim1
4985  format('pDRY=',1pe11.4,'Pa',2x,
     >       'T(ice)=',0pf6.2,'K - T(vap)=',f6.2,'K',2x,
     >       'Tlim1=',f5.1,'K')


      write(*,505) OLRmodel,TOAalbFile,RH
 505  format('OLR model: ',a16,3x,'Tab TOA albedo: ',a20,3x,
     >     'Relative humidity: ',0pf4.2)  
      
