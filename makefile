# git version enviromental variable
ifneq ("$(wildcard .git)","")
	GIT_VERSION := "$(shell git describe --dirty --always --tags)"
else
	GIT_VERSION := "no git version info available"
endif

# preprocessor flags
FPPFLAGS = -fpp -D__GIT_VERSION__=\"$(GIT_VERSION)\"

# compiler
FC = ifort

# compiler flags
CFLAGS = 

# files
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
OBJS = $(subst .f90,.o,$(SOURCE))
EXEC = kmc_tian

# debug build
DDIR   = debug
DEXEC  = $(DDIR)/$(EXEC)
DOBJS  = $(addprefix $(DDIR)/, $(OBJS))
DFLAGS = -warn all -check all -module $(DDIR)

# release build
RDIR   = release
REXEC  = $(RDIR)/$(EXEC)
ROBJS  = $(addprefix $(RDIR)/, $(OBJS))
RFLAGS = -O2 -module $(RDIR)

.PHONY: all debug release prepare remake clean

# Default rules
all: prepare release

# Debug rules
debug: $(DEXEC)

$(DEXEC): $(DOBJS) 
	$(FC) $(CFLAGS) $(DFLAGS) -o $(DEXEC) $^

$(DDIR)/%.o: %.f90
	$(FC) -c $(CFLAGS) $(DFLAGS) $(FPPFLAGS) -o $@ $^
	
# Release rules
release: $(REXEC)

$(REXEC): $(ROBJS) 
	$(FC) $(CFLAGS) $(RFLAGS) -o $(REXEC) $^

$(RDIR)/%.o: %.f90
	$(FC) -c $(CFLAGS) $(RFLAGS) $(FPPFLAGS) -o $@ $^
	
# additional rules

prepare:
	mkdir -p $(DDIR) $(RDIR)
	
remake: clean all
	
clean: 
	rm -f $(DDIR)/* $(RDIR)/*
