# Makefile created by mkmf $Id$ 


.DEFAULT:
	-touch $@
all: a.out
SRC =
OBJ =
clean: neat
	-rm -f .cppdefs $(OBJ) a.out
neat:
	-rm -f $(TMPFILES)
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
a.out: $(OBJ) 
	$(LD) $(OBJ) -o a.out  $(LDFLAGS)
