* include this AFTER parEBM.h
        real*8 Tveg,width,dead,aveg,gromax,seedfall,Vegstart
        real*8 dTdead

  
	parameter(Tveg=293.3)
        parameter(dTdead=50.0)
        parameter(Vegstart=0.3)
 	parameter(width=1/20.0)
        parameter(dead=0.01/Porb)
	parameter(gromax=0.3/Porb)
        parameter(aveg= 0.15)
        parameter(seedfall=0.01)

        real*8 annualGlobalVegZonal, annualGlobalVeg
        real*8 veg(N),vegmat(Ns,N), fveg
 

