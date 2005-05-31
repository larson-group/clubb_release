# $Id: Makefile,v 1.8 2005-05-31 15:55:06 dschanen Exp $

# Compiler invocation

FC=pgf90
#FC=g95
#FC=f95 

# Compiler flags
DEBUG = -g -C

#G95_P3_OPTIMIZE = -O3 -ffast-math -msse -march=pentium3
#SUNPRO_X86_OPTIMIZE = -fast -ftrap=%none
PGF_P4_OPTIMIZE = -O3 -fastsse -Mconcur -tp p7

#FCFLAGS = $(PGF_P4_OPTIMIZE) -DSTATS 
FCFLAGS = $(PGF_P4_OPTIMIZE) $(DEBUG) -DSTATS

# Linking Flags
LDFLAGS = $(FCFLAGS)

# Object files
OBJS	= polpak.o \
	constants.o arrays.o endian.o outputgrads.o grid.o statistics.o \
	banbks.o bandec.o erf.o tridag.o rain.o \
	closure_new.o diag_var.o gcss.o hoc.o hydrostatic.o \
	length.o mixing.o rsat.o sfc.o sounding.o update.o wp23.o tau.o \
	inputgrads.o gradsaverage.o error.o \
        autoconv_driver.o corrcoef.o gaus_condt.o \
        gaus_mixt_points.o gaussj.o generate_k_order.o \
        latin_hyper_sample.o lh_sampler.o ltqnorm.o \
        matmult.o micro_calcs.o permute_height_time.o \
        ql_estimate.o ran2.o rand_permute.o \
        rtpthlp_2_sptp.o sample_points.o st_2_rtthl.o \
        truncate_gaus_mixt.o \
	nrtype.o nr.o nrutil.o ran_state.o ran1.o amebsa.o amoeba.o \

TARGETS = int2txt.o hoc_tuner.o hoc_standalone.o compare_runs.o

# Targets
all:	hoc_tuner hoc_standalone int2txt compare_runs

# Rules
.SUFFIXES:	.o .f .F .f90

.f.o:
	$(FC) $(FCFLAGS) -c $<
.f90.o:
	$(FC) $(FCFLAGS) -c $<
.F.o:
	$(FC) $(FCFLAGS) -c $<

hoc_tuner:	$(OBJS) $(TARGETS)
		$(FC) $(LDFLAGS) -o hoc_tuner \
		$(OBJS) $@.o

hoc_standalone: $(OBJS) $(TARGETS)
	        $(FC) $(LDFLAGS) -o hoc_standalone \
                $(OBJS) $@.o 

int2txt:	$(TARGETS)
		$(FC) $(LDFLAGS) $@.o -o int2txt

compare_runs:	$(OBJS) $(TARGETS)
	        $(FC) $(LDFLAGS) -o compare_runs \
                $(OBJS) $@.o 

install: 
	-cp int2txt ../tune/ 
	-cp hoc_tuner ../tune/
	-cp hoc_standalone ../standalone/
	-cp compare_runs ../compare_runs/
clean: 
	-rm -f *.o *.mod core.*


# Dependencies
polpak.o      : polpak.f90
endian.o      : endian.F
constants.o   : constants.F
rain.o        : rain.F polpak.f90
outputgrads.o : endian.F outputgrads.F
arrays.o      : arrays.F constants.F
inputgrads.o  : endian.F inputgrads.F
grid.o        : grid.F grid.h
statistics.o  : outputgrads.F statistics.F

banbks.o      : banbks.F
bandec.o      : bandec.F
erf.o         : erf.F
tridag.o      : tridag.F

closure_new.o : closure_new.F constants.F
diag_var.o    : diag_var.F constants.F
gcss.o        : gcss.F constants.F
hoc.o         : constants.F gcss.F rain.F hoc.F 
hydrostatic.o : hydrostatic.F constants.F
length.o      : length.F constants.F
mixing.o      : mixing.F constants.F
rsat.o        : rsat.F constants.F
sfc.o         : sfc.F constants.F
sounding.o    : sounding.F
update.o      : update.F
wp23.o        : wp23.F constants.F
tau.o         : tau.F constants.F

nrtype.o      : nrtype.f90
nr.o          : nr.f90 nrtype.f90
nrutil.o      : nr.f90 nrtype.f90 nrutil.f90
ran1.o        : nr.f90 nrtype.f90 nrutil.f90 ran1.f90
ran_state.o   : nr.f90 nrtype.f90 nrutil.f90 ran_state.f90
amebsa.o      : nr.f90 nrtype.f90 nrutil.f90 ran_state.f90 ran1.f90 amebsa.f90
amoeba.o      : nr.f90 nrtype.f90 nrutil.f90 amoeba.f90
gradsaverage.o: inputgrads.F gradsaverage.F 
error.o       : nr.f90 nrtype.f90 gradsaverage.F error.F

autoconv_driver.o	: autoconv_driver.F
corrcoef.o		: corrcoef.F
gaus_condt.o		: gaus_condt.F
gaus_mixt_points.o	: gaus_mixt_points.F
gaussj.o		: gaussj.F
generate_k_order.o	: generate_k_order.F
latin_hyper_sample.o	: latin_hyper_sample.F
lh_sampler.o		: lh_sampler.F
ltqnorm.o		: ltqnorm.F
matmult.o		: matmult.F
micro_calcs.o		: micro_calcs.F constants.F
permute_height_time.o	: permute_height_time.F
ql_estimate.o		: ql_estimate.F
ran2.o			: ran2.F
rand_permute.o		: rand_permute.F
rtpthlp_2_sptp.o	: rtpthlp_2_sptp.F
sample_points.o		: sample_points.F
st_2_rtthl.o		: st_2_rtthl.F
truncate_gaus_mixt.o	: truncate_gaus_mixt.F

# Executables
hoc_tuner.o	: nr.f90 nrtype.f90 hoc.F error.o hoc_tuner.F
hoc_standalone.o: hoc.F hoc_standalone.F
int2txt.o	: int2txt.f90
compare_runs.o	: endian.F inputgrads.F gradsaverage.F compare_runs.F
