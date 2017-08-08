__DESCRIPTION__="""
modelAtmophere class
"""

class atmosphereGeometry :
   def __init__(self,p,T,molarCon,Mc,Mnc,Rstar,Name,Radius,Gravity) :
      import numpy as np
      self.Name=Name
      self.Radius=Radius
      self.Gravity=Gravity
      self.p=p.copy()
      self.T=T.copy()
      self.molarCon=molarCon
      self.Mc=Mc
      self.Mnc=Mnc
      self.Rstar=Rstar
      self.Nsamp=len(p)
      #
      # moles per unit volume
      self.n=p/T/Rstar
      #
      # molecular weight
      self.MolW=molarCon*Mc+(1-molarCon)*Mnc
      #
      self.calcZ()
      #
      self.calcMolecularDepth()
      #
      self.calcMassDepth()
   def copy(self) :
      import copy
      return copy.deepcopy(self)
   def __len__(self) : 
      return self.Nsamp
   def calcZ(self):
      """computes Z(p,T) assuming constant molecular weight between samples"""
      import numpy as np
      self.eff=np.log(self.p)
      self.Delta_eff=self.eff[1:]-self.eff[:-1]
      self.Delta_T=self.T[1:]-self.T[:-1]
      self.Mean_T=(self.T[1:]+self.T[:-1])/2
      self.Mean_MW=(self.MolW[1:]+self.MolW[:-1])/2
      #
      self.Delta_z=-self.Rstar/self.Gravity*self.Mean_T/self.Mean_MW*self.Delta_eff
      self.z=np.concatenate([[0.],self.Delta_z]).cumsum()
      self.r=self.z+self.Radius
      self.atmosphereArea=np.pi*(self.r.max()**2-self.r.min()**2)
      #
      self.Delta_z_cond=-self.Rstar/self.Gravity*self.Mean_T/self.Mc*self.Delta_eff
      self.z_cond=np.concatenate([[0.],self.Delta_z_cond]).cumsum()
      #
      self.Delta_z_ncond=-self.Rstar/self.Gravity*self.Mean_T/self.Mnc*self.Delta_eff
      self.z_ncond=np.concatenate([[0.],self.Delta_z_ncond]).cumsum()
      #
      # those quantities are used to define uncertanty from M approximation
      self.Delta_z_inf=-self.Rstar/self.Gravity*self.Mean_T/self.MolW[:-1]*self.Delta_eff
      self.z_inf=np.concatenate([[0.],self.Delta_z_inf]).cumsum()
      #
      self.Delta_z_sup=-self.Rstar/self.Gravity*self.Mean_T/self.MolW[1:]*self.Delta_eff
      self.z_sup=np.concatenate([[0.],self.Delta_z_sup]).cumsum()
      #
      self.Delta_z_unc=abs(np.array([self.Delta_z_sup-self.Delta_z,self.Delta_z_inf-self.Delta_z])).max(axis=0)
      self.z_unc=abs(np.array([self.z_sup-self.z,self.z_inf-self.z])).max(axis=0)
      #
      self.zenithal_column_density=self.calc_zenithal_column_density()
      #
      self.zenithal_column_mass=self.calc_zenithal_column_mass()
      #
   def ztropo(self) : return self.z[-1]
      #
   def calc_zenithal_column_density(self) :
      """ trapezoidal integration for line integral on the vertical"""
      import numpy as np
      dz=self.z[1:]-self.z[:-1]
      nn=0.5*(self.n[1:]+self.n[:-1])
      return np.sort(dz*nn).sum()
      #
   def calc_zenithal_column_mass(self) :
      """ trapezoidal integration for line integral on the vertical"""
      import numpy as np
      dz=self.z[1:]-self.z[:-1]
      mm=self.MolW*self.n
      nn=0.5*(mm[1:]+mm[:-1])
      return np.sort(dz*nn).sum()
      #
   def los_integral(self,i) :
      """ trapezoidal integration for line integral along a specific z[i]"""
      import numpy as np
      if i > len(self)-2 : return None
      if i < 0 : return None
      #for the ith sample computes the z along each sampled position
      xqi=(self.r[i:]**2-self.r[i]**2)**0.5
      yqi=self.n[i:]
      #trapezoidal integration
      lositg=2*np.sort((xqi[1:]-xqi[:-1])*(yqi[1:]+yqi[:-1])/2.).sum()
      return lositg,xqi,yqi
      #
   def list_los_integral(self) :
      """ splits the sampling intervall in z and for each z assumed as impact factor computes the line integral"""
      import numpy as np
      #samples integrals along each possible sample position
      LOSitg=np.zeros(len(self))
      for i in range(len(self)-1) :
         LOSitg[i]=self.los_integral(i)[0]
      return LOSitg
      #
   def calcMolecularDepth(self) :
      """total amount of kmoles seen in projection"""
      import numpy as np
      #samples trapzodial integration along each possible sample position
      self.LOS_integral=self.list_los_integral()
      #computes trapezoidal integration over all the possible sampling directions
      y=self.r*self.LOS_integral
      itg=np.sort((self.r[1:]-self.r[:-1])*(y[1:]+y[:-1])/2.).sum()
      self.MolecularDepth=2*np.pi*itg
      self.MolecularDepth_Relative_Contribution=np.concatenate([[0.],2*np.pi*((self.r[1:]-self.r[:-1])*(y[1:]+y[:-1])/2.).cumsum()])/self.MolecularDepth
      return self.MolecularDepth
      #
   def los_mass_integral(self,i) :
      """ trapezoidal integration for line integral along a specific z[i] in mass"""
      import numpy as np
      if i > len(self)-2 : return None
      if i < 0 : return None
      #for the ith sample computes the z along each sampled position
      xqi=(self.r[i:]**2-self.r[i]**2)**0.5
      yqi=self.n[i:]*self.MolW[i:]
      #trapezoidal integration
      lositg=2*np.sort((xqi[1:]-xqi[:-1])*(yqi[1:]+yqi[:-1])/2.).sum()
      return lositg,xqi,yqi
      #
   def list_los_integral_mass(self) :
      """ splits the sampling intervall in z and for each z assumed as impact factor computes the line integral in mass"""
      import numpy as np
      #samples integrals along each possible sample position
      LOSitg=np.zeros(len(self))
      for i in range(len(self)-1) :
         LOSitg[i]=self.los_mass_integral(i)[0]
      return LOSitg
      #
   def calcMassDepth(self) :
      """total amount of mass seen in projection"""
      import numpy as np
      #samples trapzodial integration along each possible sample position
      self.LOS_mass_integral=self.list_los_integral_mass()
      #computes trapezoidal integration over all the possible sampling directions
      y=self.r*self.LOS_mass_integral
      itg=np.sort((self.r[1:]-self.r[:-1])*(y[1:]+y[:-1])/2.).sum()
      self.MassDepth=2*np.pi*itg
      self.MassDepth_Relative_Contribution=np.concatenate([[0.],2*np.pi*((self.r[1:]-self.r[:-1])*(y[1:]+y[:-1])/2.).cumsum()])/self.MassDepth
      return self.MassDepth

     