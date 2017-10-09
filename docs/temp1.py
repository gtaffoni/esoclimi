import numpy as np
import pylab as plt

def temp(i):
    y,x,t=np.loadtxt('Risultati/year_lat_temp_last1.tlt',skiprows=1,unpack=True)
    xstd=x[i*54:i*54+54]
    tstd=t[i*54:i*54+54]

    bfr="season %d"%i

    plt.clf()
    plt.title(bfr)
    plt.scatter(xstd,tstd,s=3.0,color='black')
    plt.plot(xstd,tstd,label='Std',color='black')
    plt.legend(loc='best')
    plt.show()

    bfr="temp-season%02d.png"%i
    plt.title(bfr)
    plt.scatter(xstd,tstd,s=3.0,color='black')
    plt.plot(xstd,tstd,label='Std',color='black')
    plt.legend(loc='best')
    bfr="temp-season%02d.png"%i
    plt.savefig(bfr)

    
    return
