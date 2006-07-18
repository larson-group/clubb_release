# $Id: Makefile,v 1.58 2006-07-18 01:32:20 ajsmith4 Exp $
# This version should work with BSD, GNU make, smake and dmake

# For HOC we edit this by hand rather than implementing imake or automake

# Include machine specific file. Select the one appropriate for your
# machine or create a custom one.

#include config.linux_ia32_pg.in
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
	two_rt_lw_sel.o two_rt_lw_iter.o \
	two_rt_sw_bs.o bugsrad_hoc.o

PARAMETERIZATION = \
	$(LATIN_HYP) \
	polpak.o constants.o endian.o outputfile.o outputgrads.o \
	outputnc.o grid.o model_flags.o statistics.o \
	banbks.o bandec.o erf.o tridag.o rain.o \
	closure_new.o diag_var.o gcss.o rad_lwsw.o linearinterpolation.o \
        hydrostatic.o \
	length.o mixing.o rsat.o sat_rcm.o sfc.o \
	sounding.o update.o wp23.o tau.o \
	numer_check.o diag_variables.o prog_variables.o \
	parameterization_interface.o 

INTERFACING = gradsaverage.o inputfields.o inputgrads.o hoc.o 

NUM_REC	= nrtype.o nr.o nrutil.o ran_state.o ran1.o amebsa.o amoeba.o

TUNING = $(INTERFACING) $(NUM_REC) error.o

INT2TXT		= int2txt.o
HOC_TUNER	= $(TUNING) hoc_tuner.o
HOC_STANDALONE	= $(INTERFACING) hoc_standalone.o
COMPARE_RUNS	= $(INTERFACING) compare_runs.o
JACOBIAN	= $(INTERFACING) jacobian.o
HOC_INPUTFIELDS	= $(INTERFACING) hoc_inputfields.o
HOC_BUDGET_TERMS= $(TUNING) budget_terms.o hoc_tuner_budget_terms.o

# == Targets == ( hoc_tuner and jacobian assume statistics )
all:	hoc_tuner int2txt hoc_tuner_budget_terms \
	hoc_standalone hoc_inputfields compare_runs jacobian

# == Rules ==
# Note: since we're using a virtual path the make command is going to use the
# first file it finds with that name.  Be careful when making copies of files.
VPATH=.:LH:BUGSrad:HOC_parameterization

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

# == Dependencies ==
polpak.o	: polpak.f90
endian.o	: endian.F
constants.o	: constants.F
model_flags.o	: model_flags.F
outputfile.o	: outputfile.F
rain.o		: constants.o model_flags.o statistics.o tridag.o \
	polpak.o rain.F
outputgrads.o	: endian.o outputfile.o outputgrads.F
outputnc.o	: outputfile.o outputnc.F
diag_variables.o: constants.o grid.o diag_variables.F 
prog_variables.o: constants.o grid.o prog_variables.F 
numer_check.o	: prog_variables.o diag_variables.o numer_check.F
inputgrads.o	: endian.o inputgrads.F
grid.o		: grid.F
statistics.o	: outputgrads.o outputnc.o prog_variables.o diag_variables.o \
	statistics.F 

banbks.o	: banbks.F
bandec.o	: bandec.F
erf.o		: erf.F
tridag.o	: tridag.F

closure_new.o	: constants.o closure_new.F 
diag_var.o	: constants.o tridag.o statistics.o diag_var.F 
gcss.o		: constants.o statistics.o rad_lwsw.o linearinterpolation.F \
                  gcss.F
rad_lwsw.o      : rad_lwsw.F
linearinterpolation.o : linearinterpolation.F
mixing.o	: constants.o statistics.o banbks.o bandec.o mixing.F 
parameterization_interface.o : constants.o gcss.o hydrostatic.o mixing.o \
	model_flags.o numer_check.o \
	length.o statistics.o sounding.o diag_var.o rsat.o \
	grid.o sat_rcm.o update.o prog_variables.o diag_variables.o \
	sfc.o closure_new.o rain.o parameterization_interface.F 
hoc.o		: parameterization_interface.o hoc.F
hydrostatic.o	: constants.o hydrostatic.F 
length.o	: constants.o length.F 
rsat.o		: constants.o rsat.F 
sat_rcm.o	: constants.o sat_rcm.F 
sfc.o		: constants.o sfc.F 
sounding.o	: constants.o sounding.F
update.o	: tridag.o update.F
wp23.o		: constants.o statistics.o banbks.o bandec.o wp23.F 
tau.o		: constants.o tau.F 
inputfields.o	: grid.o inputgrads.o inputfields.F 

nrtype.o	: nrtype.f90
nr.o		: nrtype.o nr.f90 
nrutil.o	: nrtype.o nr.o nrutil.f90
ran1.o		: nrtype.o nr.o nrutil.o ran1.f90
ran_state.o	: nrtype.o nr.o nrutil.o ran_state.f90
amebsa.o	: nrtype.o nr.o nrutil.o ran_state.o ran1.o amebsa.f90
amoeba.o	: nrtype.o nr.o nrutil.o amoeba.f90

gradsaverage.o	: inputgrads.o gradsaverage.F 
error.o		: nrtype.o nr.o gradsaverage.o hoc.o error.F
budget_terms.o  : amoeba.o inputfields.o gradsaverage.o hoc.o \
	error.o budget_terms.F

autoconv_driver.o	: autoconv_driver.F
corrcoef.o		: corrcoef.F
gaus_condt.o		: gaus_condt.F
gaus_mixt_points.o	: gaus_mixt_points.F
gaussj.o		: gaussj.F
generate_k_order.o	: generate_k_order.F
latin_hyper_sample.o	: latin_hyper_sample.F
lh_sampler.o		: lh_sampler.F
ltqnorm.o		: constants.o ltqnorm.F
matmult.o		: matmult.F
micro_calcs.o		: constants.o micro_calcs.F 
permute_height_time.o	: permute_height_time.F
ql_estimate.o		: ql_estimate.F
ran2.o			: ran2.F
rand_permute.o		: rand_permute.F
rtpthlp_2_sptp.o	: rtpthlp_2_sptp.F
sample_points.o		: sample_points.F
st_2_rtthl.o		: st_2_rtthl.F
truncate_gaus_mixt.o	: truncate_gaus_mixt.F

kinds.o			: kinds.F90
bugsrad_planck.o	: kinds.o bugsrad_planck.F90
cloud_correlate.o	: kinds.o cloud_correlate.F90
bugsrad_physconst.o	: kinds.o bugsrad_physconst.F90
bugs_lwr.o		: kinds.o bugs_lwr.F bugsrad_planck.o gases_ckd.o 
bugs_rad.o 		: kinds.o bugs_rad.F bugs_swr.o bugs_lwr.o gascon.o \
	bugsrad_planck.o cloud_correlate.o bugsrad_physconst.o \
	gases_ckd.o cloudg.o newexp.o rayle.o
bugs_swr.o 		: kinds.o gases_ckd.o bugs_swr.F 
cloudg.o 		: kinds.o cloudg.F
comscp1.o 		: kinds.o comscp1.F
comscp2.o		: kinds.o comscp2.F
gascon.o 		: kinds.o bugsrad_physconst.o gascon.F
gases_ckd.o 		: kinds.o bugsrad_physconst.o gases_ckd.F90
newexp.o 		: kinds.o newexp.F90
rayle.o 		: kinds.o rayle.F
two_rt_lw.o 		: kinds.o two_rt_lw.F
two_rt_lw_sel.o 	: kinds.o two_rt_lw_sel.F
two_rt_lw_iter.o 	: kinds.o two_rt_lw_iter.F
two_rt_lw_ocastrndm.o	: kinds.o two_rt_lw_ocastrndm.F90
two_rt_sw.o		: kinds.o two_rt_sw.F
two_rt_sw_bs.o		: kinds.o bandec.o banbks.o two_rt_sw_bs.F
two_rt_sw_ocastrndm.o	: kinds.o two_rt_sw_ocastrndm.F90
bugsrad_hoc.o 		: kinds.o constants.o bugs_rad.o bugsrad_hoc.F90 

# == Executables ==
hoc_tuner.o		: nrtype.o nr.o hoc.o error.o hoc_tuner.F
hoc_tuner_budget_terms.o: nrtype.o nr.o hoc.o error.o \
	budget_terms.o hoc_tuner_budget_terms.F
hoc_standalone.o	: hoc.o hoc_standalone.F
hoc_inputfields.o	: hoc.o inputfields.o hoc_inputfields.F
int2txt.o		: int2txt.f90
compare_runs.o		: endian.o inputgrads.o gradsaverage.o compare_runs.F
jacobian.o		: gradsaverage.o hoc.o jacobian.F
