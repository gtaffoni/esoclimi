* include this AFTER parEBM.h
        real*8 Tveg,dTdead,width,dead,aveg,gromax,seedfall,Vegstart
        real*8 annualGlobalVegZonal, annualGlobalVeg,
     >        annualGlobalVeg_old, deltaVconv  
        real*8 veg(N),vegmat(Ns,N),fveg

        parameter(deltaVconv=1.0d-4)
* forse questo dovrebbe essere deltaTconv * seedfall?
  
* parametri tundra
        parameter(Tveg=278.15)
        parameter(dTdead=30.0)
        parameter(Vegstart=0.7)
        parameter(width=1.0/5.0)
        parameter(dead=1.0/Porb)
        parameter(gromax=2.5/Porb)
        parameter(aveg=0.2)
        parameter(seedfall=0.01)


