#!/usr/bin/python
def tAtmo(exop_par_file_name, workdir, logfile):
   """
   output a csv file to be handled by using a TabulatedMoistAdiabat

   kwywords

   #nonCondensible_name = name of non condesible
   #nonCondensible_MolecularWeight = Molercular weight of non condensible
   #Condensible_name = Molercular weight of condensible
   #Condensible_MolecularWeight = Molercular weight of condensible
   #
   # Surface_Radius=radius of planet at surface [m]
   # Surface_Gravity=gravity at surface [m/s^2]
   # Surface_Pressure = pressure [Pascal] at surface
   # Surface_T = T [K] at surface
   # Top_Pressure = minimum pressure at top of atmosphere
   #

   #columns:
   # T [k]
   # p [Pascal] total pressure
   # molarCond [0-1 fraction] fraction of condensable (0 no condensable, 1 all condensable)
   # massCond [mass] fraction of condensable
   # MolW = molecular weight
   #

   """

   from pipeline_interface import CONFIG
   CONFIG(workdir+'/src/tAtmo_rc') #reads the rc_file
   #
   #####################################

   import sys
   import os
   import time
   import shutil
   import numpy as np
   
   from atmosphereLib import MoistAdiabat as PhysMoistAdiabat
   from atmosphereLib import Rstar as Rstar
   from atmosphereGeometry import atmosphereGeometry
   from pipeline_interface import exop_par_file


   # parameters processing - paramfile passed as argument
   #if len(sys.argv) < 1+1 :
   #   print "\n \n %s exop_par_file \n"%sys.argv[0]
   #   sys.exit()

   #default values from CONFIG
   overWriteParFile=CONFIG.OverWriteParFile()
   ptop=CONFIG.top_pressure() 

   #from input - now passed as argument
   #exop_par_file_name=sys.argv[1]

   #gets the par file 
   EXOP=exop_par_file(exop_par_file_name)

   #gets relevant information
   Surface_Name=EXOP['NAME'][2]
   Surface_Radius=float(EXOP['RPLAN'][2])
   Surface_Gravity=float(EXOP['GRAV'][2])
   Surface_DryPressure=float(EXOP['PRESS'][2])
   Surface_RelativeHumidity=float(EXOP['RH'][2])
   Surface_T=float(EXOP['TMGLOB'][2])

#   print
#   print "outPfx                   : ",exop_par_file_name
#   print "ptop                     : ",ptop
#   print "Surface_Name             : ",Surface_Name
#   print "Surface_Radius           : ",Surface_Radius
#   print "Surface_Gravity          : ",Surface_Gravity
#   print "Surface_DryPressure      : ",Surface_DryPressure
#   print "Surface_RelativeHumidity : ",Surface_RelativeHumidity
#   print "Surface_T                : ",Surface_T
#   print

   logfile.write(" \n\n")
   logfile.write("   Running tAtmo code...\n\n");
   logfile.write( "outPfx                   : %s\n"%exop_par_file_name)
   logfile.write( "ptop                     : %f\n"%ptop)
   logfile.write( "Surface_Name             : %s\n"%Surface_Name)
   logfile.write( "Surface_Radius           : %f\n"%Surface_Radius)
   logfile.write( "Surface_Gravity          : %f\n"%Surface_Gravity)
   logfile.write( "Surface_DryPressure      : %f\n"%Surface_DryPressure)
   logfile.write( "Surface_RelativeHumidity : %f\n"%Surface_RelativeHumidity)
   logfile.write( "Surface_T                : %f\n"%Surface_T)
   logfile.write(" ")


   # some simple assumptions
   Ts=Surface_T
   ps=Surface_DryPressure

   # model calculation
   tic=time.time()
   MA=PhysMoistAdiabat()
   MA.ptop=ptop
   p,T,molarCon,massCon = MA(ps,Ts)

   # derived quantities
   Mc=MA.condensible.MolecularWeight
   Mnc=MA.noncon.MolecularWeight
   MolW=molarCon*Mc+(1-molarCon)*Mnc
   qL=Mc/MolW*molarCon

   maG=atmosphereGeometry(p,T,molarCon,Mc,Mnc,Rstar,Surface_Name,Surface_Radius,Surface_Gravity)

   logfile.write( " " )
   logfile.write( "ztop = %f m\n"%maG.z.max())
   logfile.write( "Natmo = %f 1e15 Kmoles\n"%(maG.MolecularDepth/1e15))
   logfile.write( " ")
   logfile.write( "Zenithal Mass Integral = %f Kgr/m2\n"%maG.zenithal_column_mass)
   logfile.write( " ")

   tic=time.time()-tic

   outPfx=exop_par_file_name.split('.par')[0]
   outPath=outPfx.split('/')[:-1]
   outPfx=outPfx.split('/')[-1]
#   outPath.append(outPfx)
   outPath='/'.join(outPath)

   ofileName=outPath+'/'+outPfx+'.agm'
   oparName=outPath+'/'+outPfx+'.par'

   oparBackupName=outPath+'/'+outPfx+'.par.bak'

#   try :
#      os.mkdir(outPath)
#   except :
#      print 
#      print "Warning: impossible to create directory '%s', may be it is already present"%outPfx
#      print
#
#   os.system('cp -rauvp %s %s'%(exop_par_file_name,oparBackupName))
   shutil.copy(exop_par_file_name,oparBackupName)
   # creates csv file
   logfile.write( "Writing %s\n"%ofileName)

   # replacement of original file name
#   os.system('cp -rauvp %s %s'%(oparName,exop_par_file_name))

   Creation_date=time.asctime()
   Creation_code=sys.argv[0]
   Creation_code_version='0.7-2017 may 10-2017 jul 3-'
   Creation_code_author='M.Maris'

   ofile=open(ofileName,'w')
   ofile.write("""#
   #
   #!File=%s
   #
   #!Creation_date=%s
   #!Creation_code=%s
   #!Creation_code_version=%s
   #!Creation_code_author=%s
   #
   """%(ofileName,Creation_date,Creation_code,Creation_code_version,Creation_code_author))

   ofile.write("""#
   #=====================================
   #!ParFile_Input=%s
   #!ParFile_Output=%s
   #!ParFile_Backup=%s
   #!Output_Path=%s
   """%(exop_par_file_name,oparName,oparBackupName,outPfx))

   ofile.write("""#
   #=====================================
   #Units are MKS:
   #Length/height in m
   #time in s
   #Mass in Kgr
   #pressure in Pascal
   #gravity in m/s^2
   #T in K
   """)

   ofile.write("""#s
   #=====================================
   #NON CONDENSIBLE DATA
   #!nonCondensible_name = %s
   #!nonCondensible_MolecularWeight = %e
   """%(MA.noncon.name,MA.noncon.MolecularWeight))
   
   ofile.write("""#
   #=====================================
   #CONDENSIBLE DATA
   #!Condensible_name = %s
   #!Condensible_MolecularWeight = %e
   """%(MA.noncon.name,MA.noncon.MolecularWeight))

   ofile.write("""#
   #=====================================
   #SURFACE DATA
   #!Surface_Name=%s
   #!Surface_Radius=%e
   #!Surface_Gravity= %e
   #!Surface_DryPressure = %e
   #!Surface_RelativeHumidity = %e
   #!Surface_T = %e
   """%(Surface_Name,Surface_Radius,Surface_Gravity,Surface_DryPressure,Surface_RelativeHumidity,Surface_T))
     
   ofile.write("""#
   #=====================================
   #Used values (some transformation if needed)
   #!ps= %e
   #!Ts= %e
   #!Top_Pressure = %e
   """%(ps,Ts,ptop))

   ofile.write("""#
   #=====================================
   #INTEGRATED DATA
   #!ztop=%20.5e
   #!zenithal_column_density=%20.5e
   #!zenithal_column_mass=%20.5e
   #!MolecularDepth15=%20.5e
   #
   #ztop = maximum value of z in m
   #zenithal_column_density=column density along the zenith, Kmoles/m^2
   #zenithal_column_mass=column mass along the zenith, Kgr/m^2
   #MolecularDepth15 = number of molecules observable in units of 10^15 Kmoles
   """%(maG.z.max(),maG.zenithal_column_density,maG.zenithal_column_mass,maG.MolecularDepth/1e15))

   ofile.write("""#
   #=====================================
   #Nsamples=%d
   #
   #Universal gas constant [mks]
   #Rstar=%e
   """%(len(p),Rstar))
   
   ofile.write("""#
   #=====================================
   #COLUMNS
   # ismp sample number
   # p [Pascal] total pressure
   # T [k]
   # molarCond [0-1 fraction] fraction of condensable 
   #           (0 not condensable, 1 all condensable)
   # massCond [mass] mass of condensable
   # MolW = molecular weight
   # z = height above the surface [m]
   # z_unc = z [m] for an atmosphere without condensable gas
   # z_inf = z [m] for Molecular weight of stratum at minimum
   # z_sup = z [m] for Molecular weight of stratum at maximum
   # LOS_integral = Kmoles/m^2 integral of molecules along LOS(z)
   # RelativeIntegral = relative contribution to MolecularDepth15 from each strata
   #
   # MolW = molarCond*Condensible_MolecularWeight + 
   #        (1-molarCond)*nonCondensible_MolecularWeight
   """)

   ofile.write("""#
   #=====================================
   #BEGIN
   ismp,p,T,molarCond,massCond,MolW,z,z_unc,z_inf,z_sup,LOS_integral,RelativeIntegral
   """)

   for ik in range(len(p)) :
      line="%i"%ik 
      line+=",%20.16e"%p[ik]
      line+=",%20.16e"%T[ik]
      line+=",%20.16e"%molarCon[ik]
      line+=",%20.16e"%massCon[ik]
      line+=",%20.16e"%MolW[ik]
      line+=",%20.16e"%maG.z[ik]
      line+=",%20.16e"%maG.z_unc[ik]
      line+=",%20.16e"%maG.z_inf[ik]
      line+=",%20.16e"%maG.z_sup[ik]
      line+=",%20.16e"%maG.z_sup[ik]
      line+=",%20.16e"%maG.LOS_integral[ik]
      line+=",%20.16e"%maG.MolecularDepth_Relative_Contribution[ik]
      line+="\n"
      ofile.write(line)

   ofile.write("""#
   #END
   #=====================================
   """)

   if not 'G_PTOP' in EXOP.keys() :
      EXOP.append('G_PTOP',"%e"%ptop,"Pressure at top for atmospheric geometry calculation [pa]","F")
      EXOP.append('G_ZTOP',"%e"%maG.z.max(),"Maximum atmosphere height [m]","F")
      EXOP.append('G_ZCD',"%e"%maG.zenithal_column_density,"Zenithal Column Density [Kmoles/m2]","F")
      EXOP.append('G_ZCM',"%e"%maG.zenithal_column_mass,"Zenithal Column Mass [Kgr/m2]","F")
      EXOP.append('G_MDP',"%e"%maG.MolecularDepth,"Projected Molecular Depth [Kmoles]","F")
      EXOP.append('COMMENT','',"G_... are parameters from atmospheric geometry calculation","")
      EXOP.append('COMMENT','',"other results in: "+outPfx+"/","")
   else :
      EXOP.update('G_PTOP',"%e"%ptop)
      EXOP.update('G_ZTOP',"%e"%maG.z.max())
      EXOP.update('G_ZCD',"%e"%maG.zenithal_column_density)
      EXOP.update('G_ZCM',"%e"%maG.zenithal_column_mass)
      EXOP.update('G_MDP',"%e"%maG.MolecularDepth)

   if overWriteParFile :
      logfile.write( "overWriting %s\n"%exop_par_file_name)
      EXOP.tofile(exop_par_file_name)
   else :
      logfile.write( "Writing %s\n"%oparName)
      EXOP.tofile(oparName)
   

   return
 
 
 
 
