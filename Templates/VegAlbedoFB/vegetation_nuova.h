* include this AFTER parEBM.h
        real*8 Tveg,dTdead,width,dead,aveg,gromax,seedfall,Vegstart
        real*8 annualGlobalVegZonal, annualGlobalVeg,
     >        annualGlobalVeg_old, deltaVconv  
        real*8 veg(N),vegmat(Ns,N),fveg

        parameter(deltaVconv=1.0d-4)
* forse questo dovrebbe essere deltaTconv * seedfall?
  
        parameter(Tveg=283.15)
        parameter(dTdead=30.0)
        parameter(Vegstart=0.9)
        parameter(width=1.0/25.0)
        parameter(dead=0.03/Porb)
        parameter(gromax=0.15/Porb)
        parameter(aveg=0.1)
        parameter(seedfall=0.01)

