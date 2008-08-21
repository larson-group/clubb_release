# Compiler invocation

FC=pgf90
#FC=g95
#FC=f95 

# Compiler flags
DEBUG = -g -C

#G95_P3_OPTIMIZE = -O3 -ffast-math -msse -march=pentium3
#SUNPRO_X86_OPTIMIZE = -fast -ftrap=%none
PGF_P4_OPTIMIZE = -O3 -fastsse -Mconcur -tp p7

FCFLAGS = $(PGF_P4_OPTIMIZE)
#FCFLAGS = $(DEBUG) 

# Linking Flags
LDFLAGS = -s $(FCFLAGS)

# Object files
OBJS90	= nrtype.o nr.o nrutil.o ran_state.o ran1.o bessj0.o amebsa.o amoeba.o

TARGETS = xamebsa.o xamoeba.o

# Target
all:	xamebsa xamoeba

$(OBJS90):
	$(FC) $(FCFLAGS) -c $*.f90 

$(TARGETS):
	$(FC) $(FCFLAGS) -c $*.f90


xamebsa:	$(OBJS90) $(TARGETS)
		$(FC) $(LDFLAGS) -o xamebsa \
		$(OBJS90) $@.o

xamoeba:	$(OBJS90) $(TARGETS)
		$(FC) $(LDFLAGS) -o xamoeba \
		$(OBJS90) $@.o
clean: 
	-rm -f *.o *.mod xamebsa xamoeba

# Dependencies

nrtype.o      : nrtype.f90
nr.o          : nr.f90 nrtype.f90
nrutil.o      : nr.f90 nrtype.f90 nrutil.f90
ran1.o        : nr.f90 nrtype.f90 nrutil.f90 ran1.f90
ran_state.o   : nr.f90 nrtype.f90 nrutil.f90 ran_state.f90
amebsa.o      : nr.f90 nrtype.f90 nrutil.f90 ran_state.f90 ran1.f90 amebsa.f90
amoeba.o      : nr.f90 nrtype.f90 nrutil.f90 amoeba.f90
bessj0.o      : bessj0.f90
xamebsa.o     : nr.f90 nrtype.f90 nrutil.f90 ran_state.f90 ran1.f90 amebsa.f90 \
		xamebsa.f90	
xamoeba.o     : nr.f90 nrtype.f90 nrutil.f90 amoeba.f90 bessj0.f90 xamoeba.f90
