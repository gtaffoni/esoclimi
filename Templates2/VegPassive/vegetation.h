* include this AFTER parEBM.h
        real*8 Tveg,width,dead,aveg,gromax,seedfall,Vegstart

  
	parameter(Tveg=293.3)
        parameter(Vegstart=0.1)
 	parameter(width=1/20.0)
        parameter(dead=0.1/Porb)
	parameter(gromax=0.3/Porb)
        parameter(aveg= 0.15)
        parameter(seedfall=0.1)

        real*8 annualGlobalVegZonal, annualGlobalVeg
        real*8 veg(N),vegmat(Ns,N)
