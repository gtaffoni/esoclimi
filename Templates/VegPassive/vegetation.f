
	real*8 function grow(Tcur)
        include 'parEBM.h'
        include 'vegetation.h'
        real*8 Tcur

        grow=gromax*(1.0-(width*(Tveg-Tcur))**2)
	if (grow.lt.0.) grow=0.

        return
        end


	real*8 function death(Tcur)
        include 'parEBM.h'
        include 'vegetation.h'
        real*8 Tcur

* vegetation can't live if temp is too far from Tveg OR negative
        death=dead
	if (Tcur.lt.(Tveg-dTdead).or.Tcur.gt.(Tveg+dTdead).or.Tcur<273.5)
     >     death=dead*100 

        return
        end

