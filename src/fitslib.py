#!/usr/bin/env python

#  fitslib.py
#  Esoclimi
#
#  Created by giuliano taffoni on 04/04/17.
#  Copyright 2017 Guliano Taffoni. All rights reserved.

"""
    
    
    create_FITS(in, out, param)
    
    
    Create a FITS file [out] from output file [in] of the Fortran Code Esopianeti
    The stracture of the FITS is given by the parameter file [param]
    It requires pyfits, math and numpy

   
"""

## this function parse the output file from fortran code
def create_FITS(in_file,out_file,param_file):
    import numpy as np
    import pyfits
    import math
    
    prihdr = pyfits.Header() # Create Header
    with open(in_file, 'r') as f:
        first_line = f.readline().rstrip().split()
    #print "Creating fits with Ns, N=",first_line
    N=int(math.floor(float(first_line[1])))
    Ns=int(math.floor(float(first_line[0])))
    prihdr['N']  = (N, 'number of latitude zones')
    prihdr['Ns']  = (Ns, 'number of time steps per orbit')

    lines = [line.rstrip('\n') for line in open(pfile)]
    for line in lines:
        if line:
            words = line.split('!')
            if words[0] == "COMMENT":
                prihdr['COMMENT'] = words[2]
            else:
                if words[3]=="F":
                    prihdr[words[0]]  = (float(words[1]),words[2])
                elif words[3]=="I":
                    prihdr[words[0]]  = (int(words[1]),words[2])
                else:
                    prihdr[words[0]]  = (words[1],words[2])
            if words[0] == "DATE":
                date=words[1]
            if words[0] == "NUMBER":
                number=words[1]

    prihdu = pyfits.PrimaryHDU(header=prihdr)

#print "Date of creation of map: ",date
    out_file=out_file+date[0:2]+"."+date[2:4]+"."+date[4:8]+"-"+number+".fits"
#print "Creating ",out_file

    data=np.loadtxt(in_file, skiprows=1)    # Read data from input file
    year, lat, temp = data.T
    year = year - math.floor(year[0])       # This subtracts the integer part of "year"

    c1=pyfits.Column(name='year', format='E', array=year)
    c2=pyfits.Column(name='lat', format='E', array=lat)
    c3=pyfits.Column(name='temp', format='E', array=temp)

    cols=pyfits.ColDefs([c1,c2,c3])

    tbhdu=pyfits.BinTableHDU.from_columns(cols)

    thdulist = pyfits.HDUList([prihdu, tbhdu])

    thdulist.writeto(out_file)
