CC     = gfortran
LNKCMD = gfortran



#OPTIMIZE     = -O3 -Wall
OPTIMIZE     = -O5 
#OPTIMIZE = -O0 -ggdb
OBJS   = codeEBM.o  core_fctns.o  integrate.o  EBMfunctions.o other_functions.o
INCL   = commons.h  functions.h  matrices.h  parEBM.h  startpar.h  var-incl.h initialize.f Makefile print_input.f open_output.f  splash_output.f module_annual_mean.f  module_annual_output.f  module_final_output.f  module_include.f  integrator_call.f  module_open_output.f  module_update_season.f  module_write_seasons.f module_albedo_surf.f

CFLAGS = $(OPTIMIZE)
LIBS   = 

codeEBM.x: $(OBJS)
	@$(LNKCMD) $(OBJS) $(LIBS) -o codeEBM.x
	@echo "===============  ::  codeEBM.x  built"


$(OBJS): $(INCL)

%.o : %.f
	@echo "...............[$(CC)]" $< "->" $@
	@$(CC) -c $(CFLAGS) $< -o $@


clean:
	@echo "cleaning object files: "	
	@echo $(OBJS) 
	@rm -f $(OBJS) 
	@rm -f codeEBM.x

cleanmain:
	@echo "cleaning source files from main directory: "	
	@echo *.f *.h *.py fo_earth_DMAP.dat Risultati/
	@rm -f *.f *.h *.py fo_earth_DMAP.dat Risultati/*

