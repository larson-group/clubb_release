################################################################################
# $Id: Makefile,v 1.128 2008-06-23 18:11:26 dschanen Exp $
#
# This version should work with GNU make and Sun's dmake.
#
# For HOC we edit this by hand rather than implementing imake or automake.
#
# Include your machine specific file here. Select the one appropriate for your
# platform or create a customized one.
#
################################################################################

include config/linux_ia32_g95_debug.in
#include config/linux_ia32_g95_optimize.in
#include config/linux_ia32_gfortran.in
#include config/linux_ia32_absoft.in
#include config/linux_ia32_pg.in
#include config/linux_ia32_ss12.in
#include config/linux_ia64_intel.in
#include config/solaris_ia32_g95.in
#include config/solaris_generic_ss12_debug.in
#include config/solaris_generic_ss12_optimize.in
#include config/osf1_alpha_digital.in

# == Object files ==

OBJDIR = $(PREFIX)/obj
SRCDIR = $(PREFIX)/src

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
	two_rt_lw_sel.o two_rt_lw_iter.o calendar.o cos_solar_zen.o \
	two_rt_sw_bs.o bugsrad_hoc.o

COAMPS_MICRO = adjmlt.o adjtq.o coamps_micro_driver.o conice.o \
        eqa10g.o eqa11g.o eqa12.o eqa12g.o \
        eqa13g.o eqa14g.o eqa15.o eqa17g.o eqa18.o \
        eqa18g.o eqa19.o eqa19g.o eqa20g.o eqa21.o \
        eqa21g.o eqa22.o eqa22g.o eqa25.o eqa26.o eqa27.o \
        eqa27r.o eqa28.o eqa5g.o eqa6.o eqa7.o eqa7g.o \
        eqa8g.o eqa9.o eqa9g.o esati.o esat_new.o \
        esatv.o frzh.o nrmcol.o nrmtqi.o \
        nrmtqw.o qsatvi.o qtadj.o slope.o tgqg.o \
        tgqi.o tgqr.o tgqs.o saturation.o gamma.o


MODEL = $(COAMPS_MICRO) \
	hoc.o sounding.o rad_lwsw.o cloud_rad.o astex.o atex.o fire.o bomex.o \
	arm.o arm_97.o arm_0003.o arm_3year.o dycoms2_rf01.o \
	dycoms2_rf02.o cloud_sed_mod.o mpace_a.o mpace_b.o \
	nov11.o jun25.o clex9_nov02.o clex9_oct14.o wangara.o \
	diag_ustar.o cobra.o rico.o gabls2.o \
	KK_microphys.o microphys_driver.o cons.o air.o parab.o \
	lba.o interpolation.o ice_dfsn_mod.o

PARAMETERIZATION = \
	$(LATIN_HYP) \
	constants.o endian.o outputfile.o outputgrads.o \
	outputnc.o grid.o model_flags.o file_functions.o \
        stats_type.o stats_hoc.o stats_subs.o \
        stats_zt.o stats_zm.o stats_sfc.o stats_prec.o sfc.o \
	erf.o diffusion.o mean_adv.o temp_in_K_mod.o \
	closure_new.o diag_var.o interpolation.o \
        hydrostatic.o mixing_length.o mixing.o saturation.o \
	compute_um_edsclrm_mod.o wp23.o Skw_func.o error_code.o \
	numer_check.o diag_variables.o prog_variables.o \
	parameterization_interface.o lapack_wrap.o  \
	parameters.o pos_definite_mod.o calendar.o \
	fill_holes_mult.o


INTERFACING = gradsaverage.o inputfields.o inputgrads.o 

NUM_REC	= nrtype.o nr.o nrutil.o ran_state.o ran1.o amebsa.o amoeba.o

TUNING = $(MODEL) $(INTERFACING) $(NUM_REC) error.o

INT2TXT		= int2txt.o
HOC_TUNER	= $(MODEL) $(INTERFACING) $(NUM_REC) error.o hoc_tuner.o
HOC_STANDALONE	= $(MODEL) $(INTERFACING) hoc_standalone.o
COMPARE_RUNS	= $(INTERFACING) compare_runs.o
JACOBIAN	= $(MODEL) $(INTERFACING) jacobian.o
HOC_INPUTFIELDS	= $(MODEL) $(INTERFACING) hoc_inputfields.o
HOC_BUDGET_TERMS= $(MODEL) $(INTERFACING) $(NUM_REC) error.o budget_terms.o hoc_tuner_budget_terms.o

# == Targets == ( hoc_tuner and jacobian assume statistics )
all:	hoc_standalone hoc_tuner int2txt hoc_inputfields  \
	compare_runs jacobian
# hoc_tuner_budget_terms
# Note:
# The budget terms tuner is currently non-working
# -dschanen 15 March 2007

# == Rules ==
# Note: since we're using a virtual path the make command is going to use the
# first file it finds with that name.  Be careful when making copies of files.
 
VPATH=.:LH:GCSS:BUGSrad:HOC_parameterization:coamps_micro:$(OBJDIR)

# Cancel implicit rules (e.g for modula-2 `m2c -c $@.mod')
.SUFFIXES:
# Define valid Fortran suffixes for HOC
.SUFFIXES: .a .o .f .F .f90 .F90 .h

.f.o:
	cd $(OBJDIR); $(FC) $(FCFLAGS) -o $(OBJDIR)/$@ -c $(SRCDIR)/$<
.f90.o:
	cd $(OBJDIR); $(FC) $(FCFLAGS) -o $(OBJDIR)/$@ -c $(SRCDIR)/$<
.F.o:
	cd $(OBJDIR); $(FC) $(FCFLAGS) -o $(OBJDIR)/$@ -c $(SRCDIR)/$<
.F90.o:
	cd $(OBJDIR); $(FC) $(FCFLAGS) -o $(OBJDIR)/$@ -c $(SRCDIR)/$<

# == Static libraries (default) ==

libhoc_param.a:	$(PARAMETERIZATION)
	cd $(OBJDIR); $(AR) $(ARFLAGS) $(LIBDIR)/$@ $(PARAMETERIZATION)
	$(RANLIB) $(LIBDIR)/$@

libhoc_bugsrad.a:	$(BUGS_RAD)
	cd $(OBJDIR); $(AR) $(ARFLAGS) $(LIBDIR)/$@ $(OBJDIR)/$(BUGS_RAD)
	$(RANLIB) $(LIBDIR)/$@

# == Shared libraries ==

libhoc_param.so:	$(PARAMETERIZATION)
	cd $(OBJDIR); $(SHARED) $(SHAREDFLAGS) -o $(LIBDIR)/$@ $(OBJDIR)/$(PARAMETERIZATION)

libhoc_bugsrad.so:	$(BUGS_RAD)
	cd $(OBJDIR); $(SHARED) $(SHAREDFLAGS) -o $(LIBDIR)/$@ $(OBJDIR)/$(BUGS_RAD)

# Remove comment to use shared libs (make sure you know what you're doing)
HOCLIBS = libhoc_param.a libhoc_bugsrad.a
#HOCLIBS = libhoc_param.so libhoc_bugsrad.so

int2txt: $(INT2TXT)
		$(FC) $(FCFLAGS) $(OBJDIR)/$@.o -o $(BINDIR)/int2txt

hoc_tuner: $(HOCLIBS) $(HOC_TUNER)
	cd $(OBJDIR); $(FC) $(FCFLAGS) -o $(BINDIR)/hoc_tuner \
	$(HOC_TUNER) $(LDFLAGS)

hoc_standalone:	$(HOCLIBS) $(HOC_STANDALONE)
	cd $(OBJDIR); $(FC) $(FCFLAGS) -o $(BINDIR)/hoc_standalone \
	$(HOC_STANDALONE) $(LDFLAGS)  

hoc_inputfields: $(HOCLIBS) $(HOC_INPUTFIELDS)
	cd $(OBJDIR); $(FC) $(FCFLAGS) -o $(BINDIR)/hoc_inputfields \
	$(HOC_INPUTFIELDS) $(LDFLAGS)

hoc_tuner_budget_terms:	$(HOCLIBS) $(HOC_BUDGET_TERMS)
	cd $(OBJDIR); $(FC) $(FCFLAGS) -o $(BINDIR)/hoc_tuner_budget_terms \
	$(OBJDIR)/$(HOC_BUDGET_TERMS) $(LDFLAGS) 

jacobian: $(HOCLIBS) $(JACOBIAN)
	cd $(OBJDIR); $(FC) $(FCFLAGS) -o $(BINDIR)/jacobian \
	$(JACOBIAN) $(LDFLAGS)

compare_runs:	$(HOCLIBS) $(COMPARE_RUNS)
	cd $(OBJDIR); $(FC) $(FCFLAGS) -o $(BINDIR)/compare_runs \
	$(COMPARE_RUNS) $(LDFLAGS)

clean: 
	-rm -f $(OBJDIR)/*.o $(OBJDIR)/*.mod $(OBJDIR)/core*

distclean: 
	-rm -f $(OBJDIR)/*.o $(OBJDIR)/*.mod $(OBJDIR)/core* \
	$(BINDIR)/hoc_tuner \
	$(BINDIR)/int2txt \
	$(BINDIR)/hoc_tuner_budget_terms \
	$(BINDIR)/hoc_standalone \
	$(BINDIR)/hoc_inputfields \
	$(BINDIR)/compare_runs \
	$(BINDIR)/jacobian \
	$(LIBDIR)/lib*

# Dependencies
include make.deps
