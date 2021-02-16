ifneq ("$(wildcard .git)","")
	GIT_VERSION := "$(shell git describe --dirty --always --tags)"
else
	GIT_VERSION := "no git version info available"
endif

CFLAGS = -O2
FPPFLAGS = -fpp -D__GIT_VERSION__=\"$(GIT_VERSION)\"
DFLAGS = -warn all -check all

FC = ifort

SOURCE = open_file.f90 constants.f90 \
	utilities.f90 temperature_laws.f90 \
	control_parameters_class.f90 \
	energy_parameters_class.f90 \
	mc_lat_class.f90 \
	energy_mod.f90 \
	mmc_mod.f90  \
	rates_hopping_class.f90 rates_desorption_class.f90 rates_dissociation_class.f90 \
	rates_association_class.f90  rates_bimolecular_class.f90 \
	reaction_class.f90  \
	kmc_mod.f90 \
	kmc.f90

kmc_tian: $(subst .f90,.o,$(SOURCE))
	$(FC) -o $@ $^
	
.PHONY: clean
clean: 
	rm -f *.o *.mod kmc_tian 
	
%.o: %.f90
	$(FC) -c $(DFLAGS) $(FPPFLAGS) -o $@ $^
