#!/bin/bash

#removing old Risultati.bak
if [ -d Risultati.bak ]; then
 rm -rf Risultati.bak
fi
#saving Risultati
if [ -d Risultati ]; then
  mv Risultati Risultati.bak
fi

cp -r CCM_RH60 Src
cd Src
rm *.f
rm Risultati/*
rm *.py
rm Makefile
rm *.dat
cd ..

# command line input (you must know what you're doing)
if [ $1 ]; then
 ibuild=$1
 iplanet=$2
 number=$3
 version=$4
 simtype=$5
else #shell input

 # choosing module
 ibuild=0
 echo "Choose your module:"
 echo "(0) Standart ESTM"
 echo "(1) ESTM + Vegetation, passive"
 echo "(2) ESTM + Vegetation, Albedo FeedBack"
 read ibuild

 #choosing planet
 iplanet=0
 echo "Choose your planet:"
 echo "(0) EARTH"
 echo "(1) Kepler 452b"
 read iplanet

 #number of simulation
 number=1
 echo "Enter simulation number (ENTER to skip): "
 read number

 #code version
 version="1.1.01"
 echo "Enter code version (e.g., 1.1.01, ENTER to skip) "
 read version

fi

if [ ! $3 ]; then
 number=0
fi
if [ ! $number ]; then
 number=0
fi

#setting build
 if [ $ibuild -eq 0 ]; then
  build='Std'
 fi
 if [ $ibuild -eq 1 ]; then
  build='VegPassive'
 fi
 if [ $ibuild -eq 2 ]; then
  build='VegAlbedoFB'
 fi

#setting planet
 if [ $iplanet -eq 0 ]; then
  planet='EARTH'
 fi 
 if [ $iplanet -eq 1 ]; then
  planet='K452b'
 fi 

echo " configuring  " $build $planet

# copying standard files
mkdir -p Src
cp -r CCM_RH60 Src
cp ModulesDef/* Src
cp Std/* Src

#comping appropriate module files
if [ $build = "VegPassive" ]; then
 cp VegPassive/* Src
 cp VegPassive/Modules/* Src
fi
if [ $build = "VegAlbedoFB" ]; then
 cp VegPassive/* Src
 cp VegPassive/Modules/* Src
 cp VegAlbedoFB/* Src
 cp VegAlbedoFB/Modules/* Src
fi

#copying planet data
if [ $planet = "EARTH"  ]; then
 cp Planets/EARTH.py Src
 cp Planets/fo_earth_DMAP.dat Src
 PLANET="EARTH.py"
fi
if [ $planet = "K452b"  ]; then
 cp Planets/K452b.py Src
 PLANET="K452b.py"
fi

# compiling and executing code (in directory Src)
cd Src
mkdir -p Risultati
rm Risultati/*

if [ $number -eq 0  ]
 then
 python runEBM.py $PLANET 
 else
 python runEBM.py $PLANET $number $version $build
fi
# moving risultati in the main dir
mv Risultati ../
rm *.o
exit
