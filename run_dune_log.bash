#!/bin/bash

#removing old Risultati.bak
if [ -d Risultati.bak ]; then
 rm -rf Risultati.bak
fi
#saving Risultati
if [ -d Risultati ]; then
  mv Risultati Risultati.bak
fi
#saving Src
if [ -d Src ]; then
  mv Src Src.bak
fi

mkdir -p Src

cp -r Templates/CCM_RH60 Src
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
 echo "(2) Dune"
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
 if [ $iplanet -eq 2 ]; then
  planet='Dune'
 fi

echo " configuring  " $build $planet

# copying standard files
cp Templates/ModulesDef/* Src
cp Templates/Std/* Src
cp Templates/Std/startpar_standalone.h Src/startpar.h
cp source/runEBM_standalone.py Src/runEBM.py
cp source/libraryEBM.py Src
cp source/constantsEBM.py Src

#comping appropriate module files
if [ $build = "VegPassive" ]; then
 cp Templates/VegPassive/* Src
 cp Templates/VegPassive/Modules/* Src
fi
if [ $build = "VegAlbedoFB" ]; then
 cp Templates/VegPassive/* Src
 cp Templates/VegPassive/Modules/* Src
 cp Templates/VegAlbedoFB/* Src
 cp Templates/VegAlbedoFB/Modules/* Src
fi

#copying planet data
if [ $planet = "EARTH"  ]; then
 cp source/EARTH.py Src
 cp Templates/Planets/fo_earth_DMAP.dat Src
 PLANET="EARTH.py"
fi
if [ $planet = "K452b"  ]; then
 cp source/K452b.py Src
 PLANET="K452b.py"
fi
if [ $planet = "Dune"  ]; then
 cp source/Dune.py Src
 PLANET="Dune.py"
fi

# compiling and executing code (in directory Src)
cd Src
mkdir -p Risultati
rm Risultati/*

if [ $number -eq 0  ]
 then
 python runEBM.py $PLANET >& log_estm
 else
 python runEBM.py $PLANET $number $version $build >& log_estm
fi
cp log_estm Risultati
# moving risultati in the main dir
mv Risultati ../
rm *.o
exit

