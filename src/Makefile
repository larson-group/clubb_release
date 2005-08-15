# $Id: Makefile,v 1.33 2005-08-15 21:17:07 dschanen Exp $
# Uses GNU make.  BSD make doesn't like the LH/*.F syntax.

# Compiler invocation

FC=pgf90 # Portland Group PGF90
#FC=g95 # G95
#FC=f95 # Sun Studio / Forte
#FC=ifort # Intel Fortran


DEBUG=-g -C#-Ktrap=inv,unf -Minform=inform 

# Optimization (assuming pentium4)
# Sun fortran currently appears to crash if statistics.F is done with 
#   optimization ( due to a bug in ube I think)

#OPTIMIZE=-O2 -msse2 -march=pentium4# g95
#OPTIMIZE=-xtarget=opteron -ftrap=%none# Sun Studio x86/Opteron
#OPTIMIZE=-O2 -xN #Intel complier, P4 w/out SSE3
#OPTIMIZE=-O2 -tp p7 -Kieee -r8# PGF90 w/ double precision reals.
OPTIMIZE=-O2 -tp p7 -Kieee# PGF90, P4

# Linking Flags
LDFLAGS = # use -s to strip

# Compiler flags
#   -DSTATS enables statistics
#   -DSCALARS enables passive scalars
#   -DIFORT -DPGF the non-standard record lengths in binary I/O require the 
#    specification of the compiler;  PGF90 defines records to 1 byte, Intel's ifort to 4 bytes.

FCFLAGS = $(OPTIMIZE) $(DEBUG) -DSTATS -DPGF
#FCFLAGS = $(OPTIMIZE) $(DEBUG) -DSTATS -DSCALARS


# Object files
OBJS	= polpak.o \
	constants.o arrays.o endian.o outputgrads.o grid.o statistics.o \
	banbks.o bandec.o erf.o tridag.o rain.o \
	closure_new.o diag_var.o gcss.o hydrostatic.o \
	length.o mixing.o rsat.o sat_rcm.o sfc.o \
	sounding.o update.o wp23.o tau.o \
	inputgrads.o gradsaverage.o inputfields.o hoc.o error.o \
	\
	nrtype.o nr.o nrutil.o ran_state.o ran1.o amebsa.o amoeba.o \
	budget_terms.o 

LH_OBJS	= autoconv_driver.o corrcoef.o gaus_condt.o \
	gaus_mixt_points.o gaussj.o generate_k_order.o \
	latin_hyper_sample.o lh_sampler.o ltqnorm.o \
	matmult.o micro_calcs.o permute_height_time.o \
	ql_estimate.o ran2.o rand_permute.o \
	rtpthlp_2_sptp.o sample_points.o st_2_rtthl.o \
	truncate_gaus_mixt.o

TARGETS	= int2txt.o hoc_tuner.o hoc_standalone.o compare_runs.o \
	jacobian.o hoc_tuner_budget_terms.o hoc_inputfields.o

# Targets ( hoc_tuner and jacobian assume statistics )
all:	hoc_tuner int2txt hoc_tuner_budget_terms \
	hoc_standalone hoc_inputfields compare_runs jacobian 

# Rules
.SUFFIXES:	.o .f .F .f90

.f.o:
	$(FC) $(FCFLAGS) -c $<
.f90.o:
	$(FC) $(FCFLAGS) -c $<
.F.o:
	$(FC) $(FCFLAGS) -c $<
$(LH_OBJS):
	$(FC) $(FCFLAGS) -c $<

hoc_tuner:	$(OBJS) $(LH_OBJS) $(TARGETS)
		$(FC) $(FCFLAGS) $(LDFLAGS) -o hoc_tuner \
		$(OBJS) $(LH_OBJS) $@.o

hoc_standalone: $(OBJS) $(LH_OBJS) $(TARGETS)
		$(FC) $(FCFLAGS) $(LDFLAGS) -o hoc_standalone \
		$(OBJS) $(LH_OBJS) $@.o 

hoc_inputfields:$(OBJS) $(LH_OBJS) $(TARGETS)
		$(FC) $(FCFLAGS) $(LDFLAGS) -o hoc_inputfields \
		$(OBJS) $(LH_OBJS) $@.o

int2txt:	$(TARGETS)
		$(FC)$(FCFLAGS)  $(LDFLAGS) $@.o -o int2txt

compare_runs:	$(OBJS) $(LH_OBJS) $(TARGETS)
		$(FC) $(FCFLAGS) $(LDFLAGS) -o compare_runs \
		$(OBJS) $(LH_OBJS) $@.o 

jacobian:	$(OBJS) $(LH_OBJS) $(TARGETS)
		$(FC) $(FCFLAGS) $(LDFLAGS) -o jacobian \
		$(OBJS) $(LH_OBJS) $@.o 

hoc_tuner_budget_terms:	$(OBJS) $(LH_OBJS) $(TARGETS)
			$(FC) $(FCFLAGS) $(LDFLAGS) -o hoc_tuner_budget_terms \
			$(OBJS) $(LH_OBJS) $@.o 

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

# Dependencies
polpak.o	: polpak.f90
endian.o	: endian.F
constants.o	: constants.F
rain.o		: constants.F polpak.f90 rain.F 
outputgrads.o	: endian.F outputgrads.F
arrays.o	: arrays.F constants.F
inputgrads.o	: endian.F inputgrads.F
grid.o		: grid.F grid.h
statistics.o	: gcss.F outputgrads.F statistics.F arrays.F 

banbks.o	: banbks.F
bandec.o	: bandec.F
erf.o		: erf.F
tridag.o	: tridag.F

closure_new.o	: closure_new.F constants.F
diag_var.o	: constants.F tridag.F diag_var.F 
gcss.o		: constants.F gcss.F 
mixing.o	: constants.F statistics.F banbks.F bandec.F mixing.F 
hoc.o		: constants.F gcss.F hydrostatic.F mixing.F rain.F \
		  length.F statistics.F sounding.F diag_var.F rsat.F \
		  grid.F sat_rcm.F update.F arrays.F hoc.F
hydrostatic.o	: constants.F hydrostatic.F 
length.o	: constants.F length.F 
rsat.o		: constants.F rsat.F 
sat_rcm.o	: constants.F sat_rcm.F 
sfc.o		: constants.F sfc.F 
sounding.o	: constants.F sounding.F
update.o	: tridag.F update.F
wp23.o		: constants.F statistics.F banbks.F bandec.F wp23.F 
tau.o		: constants.F tau.F 
inputfields.o	: grid.F inputgrads.F inputfields.F 

nrtype.o	: nrtype.f90
nr.o		: nr.f90 nrtype.f90
nrutil.o	: nr.f90 nrtype.f90 nrutil.f90
ran1.o		: nr.f90 nrtype.f90 nrutil.f90 ran1.f90
ran_state.o	: nr.f90 nrtype.f90 nrutil.f90 ran_state.f90
amebsa.o	: nr.f90 nrtype.f90 nrutil.f90 ran_state.f90 ran1.f90 amebsa.f90
amoeba.o	: nr.f90 nrtype.f90 nrutil.f90 amoeba.f90
gradsaverage.o	: inputgrads.F gradsaverage.F 
error.o		: nr.f90 nrtype.f90 gradsaverage.F error.F
budget_terms.o  : nr.f90 inputfields.F gradsaverage.F hoc.F \
		  error.F budget_terms.F

autoconv_driver.o	: LH/autoconv_driver.F
corrcoef.o		: LH/corrcoef.F
gaus_condt.o		: LH/gaus_condt.F
gaus_mixt_points.o	: LH/gaus_mixt_points.F
gaussj.o		: LH/gaussj.F
generate_k_order.o	: LH/generate_k_order.F
latin_hyper_sample.o	: LH/latin_hyper_sample.F
lh_sampler.o		: LH/lh_sampler.F
ltqnorm.o		: LH/ltqnorm.F
matmult.o		: LH/matmult.F
micro_calcs.o		: LH/micro_calcs.F constants.F
permute_height_time.o	: LH/permute_height_time.F
ql_estimate.o		: LH/ql_estimate.F
ran2.o			: LH/ran2.F
rand_permute.o		: LH/rand_permute.F
rtpthlp_2_sptp.o	: LH/rtpthlp_2_sptp.F
sample_points.o		: LH/sample_points.F
st_2_rtthl.o		: LH/st_2_rtthl.F
truncate_gaus_mixt.o	: LH/truncate_gaus_mixt.F

# Executables
hoc_tuner.o		: nr.f90 nrtype.f90 hoc.F error.F hoc_tuner.F
hoc_tuner_budget_terms.o: nr.f90 nrtype.f90 hoc.F error.F \
			  budget_terms.F hoc_tuner_budget_terms.F
hoc_standalone.o	: hoc.F hoc_standalone.F
hoc_inputfields.o	: hoc.F inputfields.F hoc_inputfields.F
int2txt.o		: int2txt.f90
compare_runs.o		: endian.F inputgrads.F gradsaverage.F compare_runs.F
jacobian.o		: gradsaverage.F hoc.F jacobian.F
