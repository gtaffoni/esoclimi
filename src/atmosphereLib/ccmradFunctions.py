import climt_lite as climt
from ClimateUtilities import *
import phys

#Set default values of globals
g = 9.8 #Gravity



#ToDo: See warnings in the comments on the OLR function
#
#      General questions: Does NCAR radiation us mass mixing ratio
#        of water, or specific humidity? This affects high pressure
#        and very warm calculations.
#        Does ClimT use the same water specification as NCAR, or does
#        it do a translation?
#
#
#ToDo: Generalize the specification of the profile, so we
#      can easily switch to a dry adiabat or a nonideal gas
#      adiabat. It can be handled by using a plug-compatible
#      adiabat function.
#
#
#      The argument list for the radiation object r(...) has
#      gotten unwieldy. Make up an object to simplify this
#
#      The setting of the surface albedos to zero has not yet been
#      implemented, nor is there a provision to set the albedos.
#      The default values are .07.  We should set default
#      to zero, and allow user to set albedos
#
#      The ccmrad package uses mb for pressure and ppmv for greenhouse
#      gas concentrations.  This ought to be harmonized with the
#      miniClimt radiation package, which uses Pascals for pressure
#      and mass specific concentrations.
#
#      Make a radcomp function analogous to miniClimt, so that
#      this package can be used to do radiative-convective models
#      in a way analogous to miniClimt, without calling climt directly.
#      radcompCCM can return, [flux,heatIR,heatSolar], and can
#      take a few extra composition arguments (co2, qH2O,ch4,o3 ...).
#      Currently this isn't used in Chapter 4, but it could be useful for
#      later chapters.
#
#      The call to setTstrat, which is supposed to
#      set the globals, doesn't work reliably.  Check this out.
#      What is more reliable is to "import ccmradFunctions as ccm"
#      then do things like ccm.Tstrat = 200.
#      This may also be a problem with setGravity
#
#      Changing the liquid particle size doesn't seem to affect the OLR
#      Is this a bug, or something that ccmrad assumes? The OLR
#      does respond correctly to the ice particle size, though,
#      which is what was used in the previous version of climt
#      to make the plots of particle size effect in Chapter 5.
#
#      
#
#Remarks:  Note  that the ground temperature Ts specified in ClimT and
#NCAR doesn't just affect the radiation through sigma Ts**4. Evidently
#it does an estimate of the boundary layer (subgrid) opacity and
#uses this to compute the net upward flux.  The resulting flux gives
#the net cooling of the surface itself correctly. It's probably
#necessary to do this because even just the water in 2 meters
#of air has substantial opacity.

#Changes: Surface pressure now incorporates pressure of CO2.
#         This makes the NCAR radiation calculation of OLR accurate
#         up to around 2 bars of CO2. The surface radiation at
#         CO2 over .2 bar seems strange in the presence of water vapor,
#         though, and so the low level fluxes may not be as accurate.
#         This needs to be checked against the KPA radiation code.
#         (Total pressure computation in GetMoistProfile updated
#          6/29/2009 to be accurate at higher mole fractions of CO2.
#          Note that pAir/g is NOT the mass path of air when
#          CO2 concentration is high!)
#
#         I have not yet incorporated the effect of water vapor
#         or CH4 on surface pressure.  This is not so important
#         for CH4, which rarely gets that massive, but it is
#         relevant for water vapor in hot atmospheres. 
#
#         I have zeroed out the ozone. We could allow the user
#         to specify an ozone profile, if desired
#
#         I have replaced the moist adiabat computation by
#         the climate book version.
#
#         5/2010 I have modified this to work with the new
#         version of climt_lite.  Not yet fully checked

#General notes on the radiation computation done in
#this script:
#
#Note  that unless the vertical resolution is increased
#a great deal, this calculation may not get the OLR right
#for temperatures above 320K, when there is a lot of water
#(hence opacity) in the upper atmosphere.

#-----------------Function definitions-----------------------------
#
#
#These two routines set globals used in the
#other computations: the non-CO2 surface pressure,
#the acceleration of gravity, and the stratospheric
#temperature.
#They are used so that these infrequently changed
#parameters don't have to be included as function arguments.
def setGravity(gval):
    global g
    g = gval# Gravity in m/s**2


def setTstrat(StratosphereTemperature):
    global Tstrat
    Tstrat = StratosphereTemperature
    
#**ToDo: Modify the OLR function to allow the optional
#input of a cloud profile parameter. Do that as a keyword
#argument. (Note that AllRad essentially does all of that,
#and returns shortwave albedo information as well). 

#Utility function to set the modified moist adiabat
#temperature and moisture profiles, and set up pressure array
#Note that the co2 concentration (ppmv) is interpreted as
#the mole fraction in the total atmosphere, when the concentration
#is high. (Not the most convenient way to specify things at
#high concentrations. Would be better to specify mass path
#of CO2 and mass path of air.
def getMoistProfile(psAir,Ts,rh,co2=300.):
    ps = psAir/(1.-1.e-6*co2) #Corrected computation of ps)
    ps = ps/100. #Convert to mb.
    p = setpLin(ps)
#    p = setpLog(ps,1.,r.nlev) #Top at 1 mb
    m = phys.MoistAdiabat()
    p1,T,molecon,q = m(100.*ps,Ts,100.*p)#Convert pressure to Pa
    q = 1000.*q*rh #convert to g/kg, and set relative humidity
    #----
    #Simple isothermal Stratosphere model:
    T = numpy.where(T<Tstrat,Tstrat,T)#Stratospheric temperature
    #Set stratospheric water vapor to constant tropopause concentration
    #The following three statements are a trick to find the
    #concentration at the tropopause (where T=Tstrat) and set the
    #whole stratosphere to that value
    q = numpy.where(T>Tstrat,q,1.e30)
    qmin = min(q)
    q = numpy.where(q>1000.,qmin,q)
    return ps,p,T,q

#Define OLR function.  Greenhouse gas amounts
# (e.g. CO2 and CH4) are molar concentrations expressed
#in ppm. (essentially the mixing ratio, for small
#concentratins). Note that the effect of water vapor on surface
#pressure is not currently included, so this won't
#work properly at temperatures much above 310K.  The effect
#of added CO2 on surface pressure is included, but not the
#effect of methane; also though the pressure broadening
#from increased surface pressure is taken into account, the
#ccm radiation model does not incorporate the effect of
#the difference between self-broadening and foreign broadening
#of CO2.
#
#psAir is specified in Pa
#
#The cloud water profile can now be passed to
#the function as a keyword argument, cloudWater. .
#If the argument is present, the cloud fraction is set to 1
#at all levels
#
def OLR(psAir,Ts,rh = .5,co2=300.,ch4=1.e-30,**kwargs):
    #Check keyword arguments for cloud water
    if 'cloudWater' in kwargs.keys():
        clwp = kwargs['cloudWater']        
        #Set up cloud fraction profile
        cldf = numpy.zeros(r.nlev,numpy.Float) + 1.
    else:
        clwp = None
    #
    #Allow specification of other trace gases as keyword args
    #(so far just N2O)
    n2o =small
    if 'N2O' in kwargs.keys():
        n2o = kwargs['N2O']
    #
    #Set up the profile
    ps,p,T,q = getMoistProfile(psAir,Ts,rh,co2)
    #
    #Note that in the old version of climt, things
    #like r.params.value['co2'] = 300. would set the
    #value. Now this no longer works; the r.Params is
    #now used only for output.  The input values must
    #be set as arguments to r(...)
    #---------
    #Do the computation:
    if clwp == None:
        r(p=p,ps=ps,T=T,Ts=Ts,q=q,co2 = co2,n2o = n2o, ch4=ch4,g=g,o3 = small) #Last arg zeroes out ozone
    else:
        r(p=p,ps=ps,T=T,Ts=Ts,q=q,co2=co2,ch4=ch4,n2o=n2o,clwp=clwp,cldf=cldf,g=g,o3 = small) #Cloud case
    return -r.lwflx[0][0][0] #Extra indices make sure we return a number not array

#This is just like the OLR function, except it returns
#the shortwave albedo  as well as the OLR
#**ToDo: Make this optionally return longwave and shortwave
#heating rates, for use in radiative-convective models.
#Also profiles of atmospheric solar absorption? Also allow
#for changing ozone?
#
#ToDo: Put particle size in the call to the radiation
#method.  It's not there now, because the old version of climT
#left r_liq out of the kwarg list.
#
def AllRad(psAir,Ts,rh = .5,co2=300.,ch4=1.e-30,**kwargs):
    #Check if we want to return profiles of flux and heating
    #If doProfile = True, we will return flux and heating profiles
    #instead of just the albedo and OLR.  This is useful for doing
    #radiative-convective models.
    #
    #(Note: to be useful for radiative-convective models,
    #we have to have an option to input T and q as well. )
    doProfile = False
    if 'doProfile' in kwargs.keys():
        doProfile = kwargs['doProfile']
    #
    #Check for specification of ozone profile
    o3 = small
    if 'o3' in kwargs.keys():
        o3 = kwargs['o3']
    #Allow specification of other trace gases as keyword args
    #(so far just N2O)
    n2o =small
    if 'N2O' in kwargs.keys():
        n2o = kwargs['N2O']
    #Check keyword arguments for cloud water
    if 'cloudWater' in kwargs.keys():
        clwp = kwargs['cloudWater']        
        #Set up cloud fraction profile
        cldf = numpy.zeros(r.nlev,numpy.float) + 1.
    else:
        clwp = None
    #
    #Check for particle size
    if 'r_liq' in  kwargs.keys():
        r_liq = kwargs['r_liq']
    else:
        r_liq = 10.*numpy.ones(r.nlev,numpy.float)
    if 'r_ice' in  kwargs.keys():
        r_ice = kwargs['r_ice']
    else:
        r_ice = 30.*numpy.ones(r.nlev,numpy.float)

    #Set up the profile
    ps,p,T,q = getMoistProfile(psAir,Ts,rh,co2)  
    #Do the computation:
    if clwp == None:
        r(p=p,ps=ps,T=T,Ts=Ts,q=q,co2 = co2,n2o = n2o, ch4=ch4\
          ,solin = S,zen=0.,g=g,o3 = o3) #Last arg zeroes out ozone       
    else:
        r(p=p,ps=ps,T=T,Ts=Ts,q=q,clwp=clwp,cldf=cldf,r_ice = r_ice,co2 = co2,\
          n2o = n2o, ch4=ch4,solin = S,zen=0.,g=g,o3 = o3) #Cloud case
    if doProfile:
        #**Need to append some indices?
        return p,r.lwflx,r.swflx,r.lwhr,r.swhr
    else:
        #Extra indices make sure we return a number not array
        return -r.lwflx[0][0][0],1.-r.swflx[0][0][0]/S,\
               (r.swflx[0][0][0]-r.swflx[-1][0][0])/S

#"Canonical" dry air OLR function
#psAir is specified in Pascal
#
#Note that the ccm radiation code describes the input greenhouse
#gas concentrations as "volumetric mixing ratios," but the actual
#computes the path by multiplying total pressure by the mixing ratio
#Hence, if we take into account the co2 when computing the pressure,
#we are really treating the value as a molar concentration, not as
#a mixing ratio.  The maximum would thus be 1000000 ppmv.  The ccm
#radiation code is not really supposed to be valid at concentrations
#where the distinction between mixing ratio and concentration matters,
#but treating things this way allows us to push the envelope of validity
#and at least crudely incorporate co2 self-broadening.
#
#The surface pressure adjustment makes little difference out
#to 10% CO2, but as one goes to 20% and above it starts to become
#important
def OLRDryAir(psAir,Ts,co2 = 300.,ch4 = 1.e-30,**kwargs):
    #Allow specification of other trace gases as keyword args
    #(so far just N2O)
    n2o =small
    if 'N2O' in kwargs.keys():
        n2o = kwargs['N2O']
    #
    ps = (psAir/100.)/(1.-1.e-6*co2) #Total surface pressure includes CO2
    p = setpLin(ps)
    T = Ts*(p/ps)**(2./7.) #Ignores effect of co2 on r/cp
    T = numpy.where(T<Tstrat,Tstrat,T)#Isothermal stratosphere
    q = small
    r(p=p,ps=ps,T=T,Ts=Ts,q=q,co2 = co2,n2o = n2o, ch4=ch4,g=g,o3 = small) #Last arg zeroes out ozone
    #Extra indices make sure we return a number not array
    return -r.lwflx[0][0][0]

#This is like the OLR function, except that
#it returns both the OLR and the net surface cooling.
#It also allows for the boundary layer rh to be different
#from the rh aloft. The boundary layer rh is defined as the
#rh in the lowest 50mb of the atmosphere. Finally, it
#allows the ground temperature to be specified separately
#from the surface air temperature (necessary in the computation
#of cooling, because of the very high opacity of water vapor
#
#psAir is air pressure in Pa
#
#FIXED: 11/7/2011
#       Fixed bug in climt_call r(...) where I had Ts=Ts instead
#       of Ts = Tg.  The keyword argument Ts in climt_lite is the
#       ground temperature, which I usually call Tg.
def BddRad(psAir,Ts,Tg,rh = .5,rhbdd=.8,co2=300.,ch4=1.e-30):
    #ToDo: Allow for other trace gases, as in OLR
    #
    ps = (psAir/100.)*(1.+1.e-6*co2)
    p = setpLin(ps)
    #Set the moist adiabat
    m = phys.MoistAdiabat()
    p1,T,molecon,q = m(100.*ps,Ts,100.*p)#Convert pressure to Pa
    q = 1000.*q*rh #convert to g/kg, and set relative humidity
    #---Simple isothermal Stratosphere model:
    T = numpy.where(T<Tstrat,Tstrat,T)#Stratospheric temperature
    #Set stratospheric water vapor to constant tropopause concentration
    #The following three statements are a trick to find the
    #concentration at the tropopause (where T=Tstrat) and set the
    #whole stratosphere to that value
    q = numpy.where(T>Tstrat,q,1.e30)
    qmin = min(q)
    q = numpy.where(q>1000.,qmin,q)
    #----
    #Reset the boundarylayer humidity
    for i in range(len(q)):
        if p[i]>p[-1]-50.:
            q[i] = q[i]*rhbdd/rh
    r(p=p,ps=ps,T=T,Ts=Tg,q=q,co2 = co2 , ch4=ch4,g=g,o3 = small)
    #Extra indices make sure we return a number not array
    return (-r.lwflx[0][0][0],-r.lwflx[-1][0][0])

def SurfCoeffs(psAir,Ts,rh = .5,rhbdd=.8,co2=300.,ch4=1.e-30):
    olr,IRcool1 = BddRad(psAir,Ts,Ts,rh,rhbdd,co2,ch4)
    olr,IRcool2 = BddRad(psAir,Ts,Ts+1.,rh,rhbdd,co2,ch4)
    g = IRcool1/(phys.sigma*Ts**4)
    gs = (IRcool2-IRcool1)/(4.*phys.sigma*Ts**3)
    return g,gs

#------------------------------------------------------------------------
#This function computes the coefficient determining the incremental
#atmospheric transmission caused by nonzero air/ground temperature
#temperature difference Tg-Ts.  It is used in energy balance models
#for atmospheres which are not optically thick, in which case the
#OLR depends (somewhat) on Tg-Ts.  The effect of Tg-Ts on
#OLR is TopCoeff*(Tg-Ts).  Note that this is NOT the same as the
#coefficient (1-a_+) in Eq. (6.6) in the text, which refers to
#the transmission of the net upwelling from the surface not the
#incremental upwelling. TopCoeff is the coefficient giving the linearization
#of OLR as a function of Tg, about the state Tg = Tsa.  It is mostly
#the coefficient one wants when doing computations solutions of a
#radiative-convective model for an atmosphere that is not optically
#thick (in which case the OLR depends on Tg-Tsa).
#
#Note: There was a bug in BddRad which didn't properly use the
#ground temperature, and caused TopCoeff to return zero. This
#was fixed 11/7/2011
def TopCoeff(psAir,Ts,rh = .5,rhbdd=.8,co2=300.,ch4=1.e-30):
    olr1,IRcool = BddRad(psAir,Ts,Ts,rh,rhbdd,co2,ch4)
    olr2,IRcool = BddRad(psAir,Ts,Ts+1.,rh,rhbdd,co2,ch4)
    gtop = (olr2-olr1)/(4.*phys.sigma*Ts**3)
    return gtop

#Function to get linear fit a and b coefficient,
#for OLR = a + b*(T-T0)
def OLRcoeffs(psAir,T0,T1,rh = .5,co2=300.,ch4=1.e-30):
    a = OLR(psAir,T0,rh,co2,ch4)
    b = (OLR(psAir,T1,rh,co2,ch4)-OLR(psAir,T0,rh,co2,ch4))/(T1-T0)
    return a,b

#CLIMT moist adiabat function no longer needed, since
#we now have our own. But remember to convert q to g/kg!

#Utility function to set up a logarithmic pressure grid
def setpLog(ps,ptop,n):
	pfact = (ps/ptop)**(1./(n-1))
	p = [ptop]
	for i in range(1,n):
		p.append( pfact*p[i-1])
	return numpy.array(p)

#Utility function to set up a linear pressure grid
#  (Fix this so that the lowest level is pbot? But
#   NCAR wants pressures on half-levels, staggered from surface
def setpLin(pbot):
    return  (numpy.arange(r.nlev, dtype ='d')+ 0.5 ) * pbot/r.nlev


#To Do: Incorporate a simple nth order polynomial fit.
#Alternately, could use routine polint from Numerical Recipes
#to generate OLR from a table.

#Initialize globals
Tstrat = 150. #Stratospheric temperature
S = 300. #Incoming solar, for shortwave calculation
         #This is now passed as an argument to the call
         #(used in AllRad only, to compute albedos)

#---Instantiate the radiation model---
#REMINDER: climt wants moisture mass mixing ratio in g/kg !
#          Please keep that in mind if you modify the script to
#          set the q profile yourself.
r = climt.radiation()
#
#Set albedos (for shortwave calculation
#**ToDo: These no longer work. Need to set them as arguments
#The default value of the albedos is .07
#r.Params.value['aldir']= 0.
#r.Params.value['aldif']= 0.
#r.Params.value['asdir'] = 0.
#r.Params.value['asdif'] = 0.

#Define array of near-zeros for zeroing out ozone, water vapor, etc
small = numpy.zeros(r.nlev,numpy.float)+1.e-30
