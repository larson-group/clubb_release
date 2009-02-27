################################################################################
# $Id$
#
# This version should work with GNU make only.
#
# For CLUBB we edit this by hand rather than implementing imake or automake.
#
# Include your machine specific file here. Select the one appropriate for your
# platform or create a customized one.
#
################################################################################

#include config/linux_ia32_g95_debug.in
include config/linux_ia32_g95_optimize.in
#include config/linux_ia32_gfortran.in
#include config/linux_ia32_absoft.in
#include config/linux_ia32_pg.in
#include config/linux_ia32_ss12_optimize.in
#include config/linux_ia32_ss12_debug.in
#include config/linux_ia32_intel.in
#include config/linux_ia64_intel.in
#include config/solaris_ia32_g95.in
#include config/solaris_generic_ss12_debug.in
#include config/solaris_generic_ss12_optimize.in
#include config/osf1_alpha_digital.in

# == Object files ==

OBJDIR = $(PREFIX)/obj
SRCDIR = $(PREFIX)/src

COAMPS_MICRO_DIR = COAMPS_micro
NR_DIR = Numerical_recipes
UNREL_DIR = Benchmark_cases/unreleased_cases
#UNREL_DIR = Unreleased_cases

# == Object files with licensing restrictions ==
# We cannot freely distribute COAMPS microphysics or Numerical Recipes code.
#   For those w/o access to COAMPS and NR, the -DCOAMPS compiler flag 
# should *not* be used, and a checkout of clubb_dist will omit
# the COAMPS and NR directories, and so COAMPS_MICRO = null.
#   For those w/ access to COAMPS and NR, the -DCOAMPS flag can be used,
# and a checkout of clubb_minimal, e.g. will contain the COAMPS and NR
# directories, and COAMPS_MICRO = list of COAMPS .o files.
# In either case, CLUBB should compile.  
COAMPS_MICRO = $(patsubst %.F,%.o, $(notdir $(wildcard $(COAMPS_MICRO_DIR)/*.F)))
ifeq ($(COAMPS_MICRO),)
# Append nothing
else
FCFLAGS := $(FCFLAGS) -DCOAMPS_MICRO
endif

NUM_REC = $(patsubst %.f90,%.o, $(notdir $(wildcard $(NR_DIR)/*.f90)))
ifeq ($(NUM_REC),)
# Append nothing
else
FCFLAGS := $(FCFLAGS) -DTUNER
endif

UNRELEASED = $(patsubst %.f90,%.o, $(notdir $(wildcard $(UNREL_DIR)/*.F90)))
ifeq ($(UNRELEASED),)
# Append nothing
else
FCFLAGS := $(FCFLAGS) -DUNRELEASED_CODE 
endif

# == Other object files ==
 
LATIN_HYP = math_utilities.o \
	lh_sampler.o \
	micro_calcs.o permute_height_time.o \
	random.o \
	matrix_operations.o 

BUGS_RAD = kinds.o std_atmosphere_mod.o bugsrad_planck.o bugsrad_physconst.o \
	gases_ckd.o \
	cloud_correlate.o bugs_lwr.o bugs_swr.o bugs_rad.o \
	cloudg.o comscp1.o comscp2.o gascon.o \
	newexp.o rayle.o two_rt_lw.o two_rt_sw.o \
	two_rt_sw_ocastrndm.o two_rt_lw_ocastrndm.o \
	two_rt_lw_sel.o two_rt_lw_iter.o cos_solar_zen.o \
	two_rt_sw_bs.o bugsrad_driver.o

MODEL = $(COAMPS_MICRO) \
	clubb_driver.o sounding.o rad_lwsw.o cloud_rad.o astex.o \
	atex.o fire.o bomex.o \
	arm.o arm_97.o arm_0003.o arm_3year.o dycoms2_rf01.o \
	dycoms2_rf02.o cloud_sed_mod.o mpace_a.o mpace_b.o \
	nov11.o jun25.o clex9_nov02.o clex9_oct14.o wangara.o \
	diag_ustar.o cobra.o rico.o gabls2.o gabls3.o \
	KK_microphys_module.o microphys_driver.o \
	Parabolic_constants.o AiryFunction.o Parabolic.o \
	lba.o ice_dfsn_mod.o coamps_micro_driver.o soil_vegetation.o

PARAMETERIZATION = \
	$(LATIN_HYP) \
	constants.o endian.o output_file_module.o output_grads.o \
	output_netcdf.o grid.o model_flags.o file_functions.o \
	stats_type.o stats_variables.o stats_subs.o \
	stats_zt.o stats_zm.o stats_sfc.o stats_precision.o sfc.o \
	erf.o diffusion.o mean_adv.o hyper_diffusion_4th_ord.o \
        T_in_K_mod.o pdf_closure_module.o advance_xp2_xpyp_module.o \
	hydrostatic.o mixing_length.o advance_xm_wpxp_module.o saturation.o \
	advance_windm_edsclrm_module.o advance_wp2_wp3_mod.o Skw_module.o \
	error_code.o numer_check.o variables_diagnostic_module.o \
	variables_prognostic_module.o clubb_core.o lapack_wrap.o \
	parameters_tunable.o pos_definite_mod.o calendar.o fill_holes.o \
	mono_flux_limiter.o clip_semi_implicit.o clip_explicit.o \
	array_index.o interpolation.o parameter_indices.o \
	parameters_model.o parameters_microphys.o


INTERFACING = grads_common.o inputfields.o inputgrads.o 

TUNING = $(MODEL) $(INTERFACING) $(NUM_REC) error.o

INT2TXT		  = int2txt.o
CLUBB_TUNER	  = $(MODEL) $(INTERFACING) $(NUM_REC) error.o clubb_tuner.o
CLUBB_STANDALONE  = $(MODEL) $(INTERFACING) clubb_standalone.o
JACOBIAN	  = $(MODEL) $(INTERFACING) jacobian.o
CLUBB_INPUTFIELDS = $(MODEL) $(INTERFACING) clubb_inputfields.o
#CLUBB_BUDGET_TERMS= $(MODEL) $(INTERFACING) $(NUM_REC) error.o budget_terms.o clubb_tuner_budget_terms.o

# List of files to check for compiler warnings
WARNABLE = $(filter std_atmosphere_mod.o calendar.o cos_solar_zen.o bugsrad_driver.o, $(BUGS_RAD) ) \
	$(filter coamps_micro_driver.o saturation.o, $(COAMPS_MICRO) ) \
	$(filter-out AiryFunction.o Parabolic.o Parabolic_Constants.o $(COAMPS_MICRO), $(MODEL) ) \
	$(PARAMETERIZATION) 

# == Targets == ( clubb_tuner and jacobian assume statistics )
all:	clubb_standalone clubb_tuner int2txt clubb_inputfields  \
	jacobian
# clubb_tuner_budget_terms
# Note:
# The budget terms tuner is currently non-working
# -dschanen 15 March 2007
	perl $(PREFIX)/utilities/CLUBBStandardsCheck.pl $(PREFIX)/src/*.F90
	perl $(PREFIX)/utilities/CLUBBStandardsCheck.pl $(PREFIX)/src/CLUBB_core/*.F90
	perl $(PREFIX)/utilities/CLUBBStandardsCheck.pl $(PREFIX)/src/Benchmark_cases/*.F90
	perl $(PREFIX)/utilities/CLUBBStandardsCheck.pl $(PREFIX)/src/Latin_hypercube/*.F90


# == Rules ==
# Note: since we're using a virtual path the make command is going to use the
# first file it finds with that name.  Be careful when making copies of files.
 
VPATH=.:Latin_hypercube:Benchmark_cases:Benchmark_cases/unreleased_cases:BUGSrad:CLUBB_core:$(OBJDIR):$(NR_DIR):$(COAMPS_MICRO_DIR)

# Cancel implicit rules (e.g for modula-2 `m2c -c $@.mod')
.SUFFIXES:
# Define valid Fortran suffixes for CLUBB
.SUFFIXES: .a .o .f .F .f90 .F90 .h

.f.o:
	cd $(OBJDIR); $(FC) $(FCFLAGS) $(if $(filter $@, $(WARNABLE) ), $(WARNINGS) ) -o $(OBJDIR)/$@ -c $(SRCDIR)/$<  
.f90.o:
	cd $(OBJDIR); $(FC) $(FCFLAGS) $(if $(filter $@, $(WARNABLE) ), $(WARNINGS) ) -o $(OBJDIR)/$@ -c $(SRCDIR)/$<  
.F.o:
	cd $(OBJDIR); $(FC) $(FCFLAGS) $(if $(filter $@, $(WARNABLE) ), $(WARNINGS) ) -o $(OBJDIR)/$@ -c $(SRCDIR)/$<  
.F90.o:
	cd $(OBJDIR); $(FC) $(FCFLAGS) $(if $(filter $@, $(WARNABLE) ), $(WARNINGS) ) -o $(OBJDIR)/$@ -c $(SRCDIR)/$<  
	

# == Static libraries (default) ==

libclubb_param.a:	$(PARAMETERIZATION)
	cd $(OBJDIR); $(AR) $(ARFLAGS) $(LIBDIR)/$@ $(PARAMETERIZATION)
	$(RANLIB) $(LIBDIR)/$@

libclubb_bugsrad.a:	$(BUGS_RAD)
	cd $(OBJDIR); $(AR) $(ARFLAGS) $(LIBDIR)/$@ $(OBJDIR)/$(BUGS_RAD)
	$(RANLIB) $(LIBDIR)/$@

# == Shared libraries ==

libclubb_param.so:	$(PARAMETERIZATION)
	cd $(OBJDIR); $(SHARED) $(SHAREDFLAGS) -o $(LIBDIR)/$@ $(OBJDIR)/$(PARAMETERIZATION)

libclubb_bugsrad.so:	$(BUGS_RAD)
	cd $(OBJDIR); $(SHARED) $(SHAREDFLAGS) -o $(LIBDIR)/$@ $(OBJDIR)/$(BUGS_RAD)

# Remove comment to use shared libs (make sure you know what you're doing)
CLUBBLIBS = libclubb_param.a libclubb_bugsrad.a
#CLUBBLIBS = libclubb_param.so libclubb_bugsrad.so

int2txt: $(INT2TXT)
		$(FC) $(FCFLAGS) $(OBJDIR)/$@.o -o $(BINDIR)/int2txt

clubb_tuner: $(CLUBBLIBS) $(CLUBB_TUNER)
	cd $(OBJDIR); $(FC) $(FCFLAGS) -o $(BINDIR)/clubb_tuner \
	$(CLUBB_TUNER) $(LDFLAGS)

clubb_standalone: $(CLUBBLIBS) $(CLUBB_STANDALONE)
	cd $(OBJDIR); $(FC) $(FCFLAGS) -o $(BINDIR)/clubb_standalone \
	$(CLUBB_STANDALONE) $(LDFLAGS)  

clubb_inputfields: $(CLUBBLIBS) $(CLUBB_INPUTFIELDS)
	cd $(OBJDIR); $(FC) $(FCFLAGS) -o $(BINDIR)/clubb_inputfields \
	$(CLUBB_INPUTFIELDS) $(LDFLAGS)

#clubb_tuner_budget_terms:	$(CLUBBLIBS) $(CLUBB_BUDGET_TERMS)
#	cd $(OBJDIR); $(FC) $(FCFLAGS) -o $(BINDIR)/clubb_tuner_budget_terms \
#	$(OBJDIR)/$(CLUBB_BUDGET_TERMS) $(LDFLAGS) 

jacobian: $(CLUBBLIBS) $(JACOBIAN)
	cd $(OBJDIR); $(FC) $(FCFLAGS) -o $(BINDIR)/jacobian \
	$(JACOBIAN) $(LDFLAGS)

clean: 
	-rm -f $(OBJDIR)/*.o $(OBJDIR)/*.mod $(OBJDIR)/core*

distclean: 
	-rm -f $(OBJDIR)/*.o $(OBJDIR)/*.mod $(OBJDIR)/core* \
	$(BINDIR)/clubb_tuner \
	$(BINDIR)/int2txt \
	$(BINDIR)/clubb_standalone \
	$(BINDIR)/clubb_inputfields \
	$(BINDIR)/jacobian \
	$(LIBDIR)/lib*

# Dependencies
include make.deps
