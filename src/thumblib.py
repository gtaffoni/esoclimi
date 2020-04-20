#!/usr/bin/python

import sys
import os
from posix import system
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm  
from mpl_toolkits.axes_grid1.inset_locator import inset_axes



def create_THUMBNAILS(inpfile, out_file, date, number):


    out_file=out_file+date[0:2]+"."+date[2:4]+"."+date[4:8]+"-%04d"%number+".png"

    f=open(inpfile,'r')

    line=f.readline()
    e=line.split()
    Nt=int(e[0])   # number of pixels in time
    Nl=int(e[1])   # number of pixels in latitude

    time=np.zeros((Nt),float)      # inizialization of empty vector for time
    lat =np.zeros((Nl),float)      # inizialization of empty vector for latitudes
    tempmat=np.zeros((Nl,Nt),float)   # inizialization of empty images for temperatures

    for i in range(Nt):
        for j in range(Nl):
            line=f.readline()
            e=line.split() 
            try:
                time[i]=float(e[0]) #%int(float(e[0]))
                lat[j] =float(e[1])
                tempmat[j,i]=float(e[2])  # image is transposed: [j,i] instead of [i,j]
            except:
                print 'Not converged'
                return

    time=time-int(time[0])
#plt.figure(figsize=(10,8)) 

    Tmin = np.min(tempmat) - 0.05*np.abs(np.min(tempmat))
    Tmax = np.max(tempmat) + 0.05*np.abs(np.max(tempmat))


#fig,ax=plt.subplots(1,1,1) 

#    Tmin=200 #215.
#    Tmax=330 #315.
    Tstep=5.
    TempLevels=np.arange(Tmin,Tmax+0.5*Tstep,Tstep)

    plt.title('Temperature seasonal cycle ')

    plt.xlabel('Time (yr)',fontsize=20) 
    plt.ylabel('Latitude ($^\circ$N)',fontsize=20)  

    plt.xticks([0.,0.2,0.4,0.6,0.8,1.0],fontsize=18)
    plt.yticks([-80.,-60.,-40,-20.,0.,20.,40,60.,80.],fontsize=18)

    plt.axis([0,1,-90,90])

    plt.contourf(time,lat,tempmat,levels=TempLevels,cmap=plt.cm.rainbow
                 ,norm=plt.Normalize(vmax=Tmax,vmin=Tmin)) #,alpha=0.5) 
    cb=plt.colorbar(shrink=1.0,ticks=np.arange(Tmin,Tmax+1,10))
    cb.set_label('Temperature [K]') 
    plt.savefig(out_file)
    plt.clf()

    return



