	print *,' '
	print *,' '
	print *,'******************************************************'
<<<<<<< HEAD
	print *,'*  Welcome to ESTM 1.1.1                             *'
=======
	print *,'*  Welcome to ESTM 1.1.04                            *'
>>>>>>> EsoMPI
        print *,'*  estimating latitude dependent surface temperature *' 
	print *,'*  and vegetation fraction                           *' 
	print *,'*  of an (exo) planet                                *' 
	print *,'*  In this version vegetation reduces the albedo     *'
        print *,'*  of continents                                     *'
        write(*,6524) " *  Albedo of vegetation is ",aveg, 
     >    "                     *"
	print *,'******************************************************'
	print *,' '
	print *,' '


 6524   format(a28,f5.3,a22)
