#Obliquita': 0., 15., 23.43929, 30., 45.
#Eccentricita':  0.0, 0.01671022, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8
#Semi-major axis: 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5
#Pressioni: 0.01, 0.0177827941004, 0.0316227766017, 0.056234132519,
#0.1, 0.177827941004,, 0.316227766017, 0.56234132519, 1.0,
#1.77827941004, 3.16227766017, 5.6234132519, 10.0
#CO2: 0.1,1,10, 100 volte il valore terrestre. SOSTITUIAMO PER ORA CON LA GEOGRAFIA
#
# (da aggiungere - non metto ancora la geografia)

import numpy as np

press=[0.01, 0.0177827941004, 0.0316227766017, 0.056234132519,1.77827941004, 3.16227766017, 5.6234132519, 8.]
ecc=[ 0.0, 0.01671022, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
sma=[0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5]
obl=[0., 15., 23.43929, 30., 45.]
gtype=[0,2,3,4] #note, only for fixed ocean fraction we use geogr
geogr=[0.1, 0.3, 0.5, 0.7, 0.9]

f=open("input-moreparams.txt","w")
f.write('#  press ecc sma obl gtype oceanfrac  \n')
for p in press:
 for e in ecc:
  for s in sma:
   for o in obl: 
    for gg in gtype:
     if gg == 0:
      for g in geogr:
       line = str(p) + ' ' + str(e) + ' ' + str(o) + ' ' + str(s) + ' ' + str(gg) + ' ' + str(g) + '\n'
       f.write(line)
     else:
      g=0.7
      line = str(p) + ' ' + str(e) + ' ' + str(o) + ' ' + str(s) + ' ' + str(gg) + ' ' + str(g) + '\n'
      f.write(line)

f.close()



