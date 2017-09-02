import numpy as np
import pylab as plt

def veg(i):
    y,x,t,t2=np.loadtxt('Risultati/year_lat_veg_last1.tlt',skiprows=1,unpack=True)
    xvpass=x[i*54:i*54+54]
    tvpass=t[i*54:i*54+54]
    
    bfr="season %d"%i

    plt.clf()
    plt.title(bfr)
    plt.scatter(xvpass,tvpass,s=3.0,color='red')
    plt.plot(xvpass,tvpass,label='Veg ',color='red')
    plt.legend(loc='best')
    plt.show()

    bfr="veg-season%02d.png"%i
    plt.title(bfr)
    plt.scatter(xvpass,tvpass,s=3.0,color='red')
    plt.plot(xvpass,tvpass,label='Veg ',color='red')
    plt.legend(loc='best')
    bfr="veg-season%02d.png"%i
    plt.savefig(bfr)

    
    return


def vegfraz(i):
    y,x,t2,t=np.loadtxt('Risultati/year_lat_veg_last1.tlt',skiprows=1,unpack=True)
    xvpass=x[i*54:i*54+54]
    tvpass=t[i*54:i*54+54]
    
    bfr="season %d"%i

    plt.clf()
    plt.title(bfr)
    plt.scatter(xvpass,tvpass,s=3.0,color='red')
    plt.plot(xvpass,tvpass,label='Veg ',color='red')
    plt.legend(loc='best')
    plt.show()

    bfr="vegfraz-season%02d.png"%i
    plt.title(bfr)
    plt.scatter(xvpass,tvpass,s=3.0,color='red')
    plt.plot(xvpass,tvpass,label='Veg ',color='red')
    plt.legend(loc='best')
    bfr="vegfraz-season%02d.png"%i
    plt.savefig(bfr)

    
    return
