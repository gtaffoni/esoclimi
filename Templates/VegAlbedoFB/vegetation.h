<<<<<<< HEAD
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
 

=======
* include this AFTER parEBM.h
        real*8 Tveg,dTdead,width,dead,aveg,gromax,seedfall,Vegstart
        real*8 annualGlobalVegZonal, annualGlobalVeg,
     >        annualGlobalVeg_old, deltaVconv  
        real*8 veg(N),vegmat(Ns,N),fveg

        parameter(deltaVconv=1.0d-4)
* forse questo dovrebbe essere deltaTconv * seedfall?
  
        parameter(Tveg=293.3)
*       parameter(Tveg=303.3)
        parameter(dTdead=30.0)
        parameter(Vegstart=0.3)
        parameter(width=1/20.0)
*        parameter(dead=0.1/Porb)
        parameter(dead=0.01/Porb)
        parameter(gromax=0.3/Porb)
        parameter(aveg= 0.15)
*        parameter(seedfall=0.1)
        parameter(seedfall=0.01)


>>>>>>> EsoMPI
