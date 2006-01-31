# $Id: Makefile,v 1.44 2006-01-31 22:26:31 golaz Exp $
# This version should work with BSD, GNU make, smake and dmake

# For HOC we edit this by hand rather than implementing imake or automake

# Include machine specific file. Select the one appropriate for your
# machine or create a custom one.

#include config.linux_ia32_pg.in
#include config.linux_ia64_intel.in
 include config.messner.in

# == Object files ==

LATIN_HYP	= autoconv_driver.o corrcoef.o gaus_condt.o \
		gaus_mixt_points.o gaussj.o generate_k_order.o \
		latin_hyper_sample.o lh_sampler.o ltqnorm.o \
		matmult.o micro_calcs.o permute_height_time.o \
		ql_estimate.o ran2.o rand_permute.o \
		rtpthlp_2_sptp.o sample_points.o st_2_rtthl.o \
		truncate_gaus_mixt.o

BUGS_RAD	= kinds.o bugsrad_planck.o bugsrad_physconst.o gases_ckd.o \
		cloud_correlate.o bugs_lwr.o bugs_swr.o bugs_rad.o \
		cloudg.o comscp1.o comscp2.o gascon.o \
		newexp.o rayle.o two_rt_lw.o two_rt_sw.o \
		two_rt_sw_ocastrndm.o two_rt_lw_ocastrndm.o \
		two_rt_lw_sel.o two_rt_lw_iter.o \
		two_rt_sw_bs.o bugsrad_hoc.o

COMMON	= $(LATIN_HYP) $(BUGS_RAD) \
	polpak.o constants.o arrays.o endian.o outputgrads.o \
	grid.o statistics.o inputgrads.o \
	banbks.o bandec.o erf.o tridag.o rain.o \
	closure_new.o diag_var.o gcss.o hydrostatic.o \
	length.o mixing.o rsat.o sat_rcm.o sfc.o \
	sounding.o update.o wp23.o tau.o \
	inputfields.o gradsaverage.o hoc.o error.o 

NUM_REC	= nrtype.o nr.o nrutil.o ran_state.o ran1.o amebsa.o amoeba.o

INT2TXT		= int2txt.o 
HOC_TUNER	= $(COMMON) $(NUM_REC) hoc_tuner.o 
HOC_STANDALONE	= $(COMMON) hoc_standalone.o 
COMPARE_RUNS	= $(COMMON) compare_runs.o 
JACOBIAN	= $(COMMON) jacobian.o 
HOC_INPUTFIELDS	= $(COMMON) hoc_inputfields.o
HOC_BUDGET_TERMS= $(COMMON) $(NUM_REC) budget_terms.o hoc_tuner_budget_terms.o 

# == Targets == ( hoc_tuner and jacobian assume statistics )
all:	hoc_tuner int2txt hoc_tuner_budget_terms \
	hoc_standalone hoc_inputfields compare_runs jacobian 

# == Rules ==
VPATH=.:LH:BUGSrad

.SUFFIXES:	.o .f .F .f90 .F90 .h

.f.o:
	$(FC) $(FCFLAGS) -c $<
.f90.o:
	$(FC) $(FCFLAGS) -c $<
.F.o:
	$(FC) $(FCFLAGS) -c $<
.F90.o:
	$(FC) $(FCFLAGS) -c $<

int2txt:	$(INT2TXT)
		$(FC) $(FCFLAGS) $(LDFLAGS) $@.o -o int2txt

hoc_tuner:	$(HOC_TUNER)
		$(FC) $(FCFLAGS) $(LDFLAGS) -o hoc_tuner \
		$(HOC_TUNER)

hoc_standalone:	$(HOC_STANDALONE)
		$(FC) $(FCFLAGS) $(LDFLAGS) -o hoc_standalone \
		$(HOC_STANDALONE)

hoc_inputfields:	$(HOC_INPUTFIELDS)
			$(FC) $(FCFLAGS) $(LDFLAGS) -o hoc_inputfields \
			$(HOC_INPUTFIELDS)

hoc_tuner_budget_terms:	$(HOC_BUDGET_TERMS)
			$(FC) $(FCFLAGS) $(LDFLAGS) -o hoc_tuner_budget_terms \
			$(HOC_BUDGET_TERMS)

jacobian:	$(JACOBIAN)
		$(FC) $(FCFLAGS) $(LDFLAGS) -o jacobian \
		$(JACOBIAN)

compare_runs:	$(COMPARE_RUNS)
		$(FC) $(FCFLAGS) $(LDFLAGS) -o compare_runs \
		$(COMPARE_RUNS)

install: 
	-cp int2txt ../tune/ 
	-cp hoc_tuner ../tune/
	-cp hoc_inputfields ../inputfields/
	-cp hoc_standalone ../standalone/
	-cp compare_runs ../compare_runs/
	-cp jacobian ../jacobian/
	-cp hoc_tuner_budget_terms ../tune_budgets/
clean: 
	-rm -f *.o *.mod core*

# == Dependencies ==
polpak.o	: polpak.f90
endian.o	: endian.F
constants.o	: constants.F
rain.o		: constants.o tridag.o polpak.o rain.F 
outputgrads.o	: endian.o outputgrads.F
arrays.o	: constants.o grid.o arrays.F 
inputgrads.o	: endian.o inputgrads.F
grid.o		: grid.F grid.h
statistics.o	: outputgrads.o arrays.o statistics.F 

banbks.o	: banbks.F
bandec.o	: bandec.F
erf.o		: erf.F
tridag.o	: tridag.F

closure_new.o	: constants.o closure_new.F 
diag_var.o	: constants.o tridag.o statistics.o diag_var.F 
gcss.o		: constants.o statistics.o gcss.F
mixing.o	: constants.o statistics.o banbks.o bandec.o mixing.F 
hoc.o		: constants.o gcss.o hydrostatic.o mixing.o rain.o \
		length.o statistics.o sounding.o diag_var.o rsat.o \
		grid.o sat_rcm.o update.o arrays.o sfc.o closure_new.o \
		hoc.F
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
cloud_correlate.o	: cloud_correlate.F90
bugsrad_physconst.o	: kinds.o bugsrad_physconst.F90
bugs_lwr.o		: bugs_lwr.F bugsrad_planck.o gases_ckd.o kinds.o 
bugs_rad.o 		: bugs_rad.F bugs_swr.o bugs_lwr.o gascon.o kinds.o \
			bugsrad_planck.o cloud_correlate.o bugsrad_physconst.o \
			gases_ckd.o cloudg.o newexp.o rayle.o
bugs_swr.o 		: bugs_swr.F gases_ckd.o kinds.o
cloudg.o 		: kinds.o cloudg.F
comscp1.o 		: comscp1.F
comscp2.o		: comscp2.F
gascon.o 		: bugsrad_physconst.o gascon.F
gases_ckd.o 		: bugsrad_physconst.o gases_ckd.F90
newexp.o 		: kinds.o newexp.F90
rayle.o 		: rayle.F
two_rt_lw.o 		: two_rt_lw.F
two_rt_lw_sel.o 	: two_rt_lw_sel.F
two_rt_lw_iter.o 	: two_rt_lw_iter.F
two_rt_lw_ocastrndm.o	: two_rt_lw_ocastrndm.F90
two_rt_sw.o		: two_rt_sw.F
two_rt_sw_bs.o		: bandec.o banbks.o two_rt_sw_bs.F
two_rt_sw_ocastrndm.o	: two_rt_sw_ocastrndm.F90
bugsrad_hoc.o 		: bugsrad_hoc.F90 bugs_rad.o constants.o

# == Executables ==
hoc_tuner.o		: nrtype.o nr.o hoc.o error.o hoc_tuner.F
hoc_tuner_budget_terms.o: nrtype.o nr.o hoc.o error.o \
			budget_terms.o hoc_tuner_budget_terms.F
hoc_standalone.o	: hoc.o hoc_standalone.F
hoc_inputfields.o	: hoc.o inputfields.o hoc_inputfields.F
int2txt.o		: int2txt.f90
compare_runs.o		: endian.o inputgrads.o gradsaverage.o compare_runs.F
jacobian.o		: gradsaverage.o hoc.o jacobian.F
