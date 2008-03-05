################################################################################
# $Id: Makefile,v 1.94 2008-03-05 17:22:52 faschinj Exp $
#
# This version should work with BSD, GNU make, smake and dmake.
#
# For HOC we edit this by hand rather than implementing imake or automake.
#
# Include your machine specific file here. Select the one appropriate for your
# machine or create a custom one.
#
################################################################################

#include config.linux_ia32_pg.in
#include config.linux_ia32_g95.in
#include config.linux_ia64_intel.in
#include config.solaris_generic_ss11.in
#include config.osf1_alpha_digital.in
include config.messner.in

# == Object files ==

LATIN_HYP = autoconv_driver.o corrcoef.o gaus_condt.o \
	gaus_mixt_points.o gaussj.o generate_k_order.o \
	latin_hyper_sample.o lh_sampler.o ltqnorm.o \
	matmult.o micro_calcs.o permute_height_time.o \
	ql_estimate.o ran2.o rand_permute.o \
	rtpthlp_2_sptp.o sample_points.o st_2_rtthl.o \
	truncate_gaus_mixt.o

BUGS_RAD = kinds.o bugsrad_planck.o bugsrad_physconst.o gases_ckd.o \
	cloud_correlate.o bugs_lwr.o bugs_swr.o bugs_rad.o \
	cloudg.o comscp1.o comscp2.o gascon.o \
	newexp.o rayle.o two_rt_lw.o two_rt_sw.o \
	two_rt_sw_ocastrndm.o two_rt_lw_ocastrndm.o \
	two_rt_lw_sel.o two_rt_lw_iter.o cos_solar_zen.o \
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
        tgqi.o tgqr.o tgqs.o rsati.o gamma.o

MODEL = hoc.o sounding.o rad_lwsw.o cloud_rad.o astex.o atex.o fire.o bomex.o \
	arm.o arm_97.o arm_0003.o arm_3year.o dycoms2_rf01.o dycoms2_rf02.o mpace_a.o mpace_b.o \
	nov11.o jun25.o wangara.o \
	linear_interpolation.o diag_ustar.o icedfs.o cobra.o rico.o gabls2.o \
	lba.o zlinterp.o file_functions.o

PARAMETERIZATION = \
	$(LATIN_HYP) \
	$(COAMPS_MICRO) \
	constants.o endian.o outputfile.o outputgrads.o \
	outputnc.o grid.o model_flags.o \
        stats_type.o stats_hoc.o stats_subs.o\
        stats_zt.o stats_zm.o stats_sfc.o \
	erf.o diffusion.o KK_microphys.o microphysics.o cons.o air.o parab.o \
	closure_new.o diag_var.o \
        hydrostatic.o mixing_length.o mixing.o rsat.o sat_rcm.o sfc.o \
	update.o wp23.o tau.o Skw_func.o \
	numer_check.o diag_variables.o prog_variables.o \
	parameterization_interface.o lapack_wrap.o rsati.o \
	parameters.o


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
# The budget terms tuner is currently non-working, and will
# be fixed soon -dschanen 15 March 2007

# == Rules ==
# Note: since we're using a virtual path the make command is going to use the
# first file it finds with that name.  Be careful when making copies of files.
 
VPATH=.:LH:GCSS:BUGSrad:HOC_parameterization:coamps_micro

# Cancel implicit rules (e.g for modula-2 `m2c -c $@.mod')
.SUFFIXES:
# Define valid Fortran suffixes for HOC
.SUFFIXES: .a .o .f .F .f90 .F90 .h

.f.o:
	$(FC) $(FCFLAGS) -c $<
.f90.o:
	$(FC) $(FCFLAGS) -c $<
.F.o:
	$(FC) $(FCFLAGS) -c $<
.F90.o:
	$(FC) $(FCFLAGS) -c $<

# == Static libraries (default) ==

libhoc_param.a:	$(PARAMETERIZATION)
	$(AR) $(ARFLAGS) $(LIBDIR)/$@ $(PARAMETERIZATION)
	$(RANLIB) $(LIBDIR)/$@

libhoc_bugsrad.a:	$(BUGS_RAD)
	$(AR) $(ARFLAGS) $(LIBDIR)/$@ $(BUGS_RAD)
	$(RANLIB) $(LIBDIR)/$@

# == Shared libraries ==

libhoc_param.so:	$(PARAMETERIZATION)
	$(SHARED) $(SHAREDFLAGS) -o $(LIBDIR)/$@ $(PARAMETERIZATION)

libhoc_bugsrad.so:	$(BUGS_RAD)
	$(SHARED) $(SHAREDFLAGS) -o $(LIBDIR)/$@ $(BUGS_RAD)

# Remove comment to use shared libs (make sure you know what you're doing)
HOCLIBS = libhoc_param.a libhoc_bugsrad.a
#HOCLIBS = libhoc_param.so libhoc_bugsrad.so

int2txt:	$(INT2TXT)
		$(FC) $(FCFLAGS) $@.o -o $(BINDIR)/int2txt

hoc_tuner:	$(HOCLIBS) $(HOC_TUNER)
		$(FC) $(FCFLAGS) -o $(BINDIR)/hoc_tuner \
		$(HOC_TUNER) $(LDFLAGS)

hoc_standalone:	$(HOCLIBS) $(HOC_STANDALONE)
		$(FC) $(FCFLAGS) -o $(BINDIR)/hoc_standalone \
		$(HOC_STANDALONE) $(LDFLAGS)  

hoc_inputfields:	$(HOCLIBS) $(HOC_INPUTFIELDS)
			$(FC) $(FCFLAGS) -o $(BINDIR)/hoc_inputfields \
			$(HOC_INPUTFIELDS) $(LDFLAGS)

hoc_tuner_budget_terms:	$(HOCLIBS) $(HOC_BUDGET_TERMS)
			$(FC) $(FCFLAGS) -o $(BINDIR)/hoc_tuner_budget_terms \
			$(HOC_BUDGET_TERMS) $(LDFLAGS) 

jacobian:	$(HOCLIBS) $(JACOBIAN)
		$(FC) $(FCFLAGS) -o $(BINDIR)/jacobian \
		$(JACOBIAN) $(LDFLAGS)

compare_runs:	$(HOCLIBS) $(COMPARE_RUNS)
		$(FC) $(FCFLAGS) -o $(BINDIR)/compare_runs \
		$(COMPARE_RUNS) $(LDFLAGS)

clean: 
	-rm -f *.o *.a *.mod core*

distclean: 
	-rm -f *.o *.mod core* \
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
