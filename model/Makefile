#------------------------------------------------------------
# <next few lines under version control, D O  N O T  E D I T>
# $Date$
# $Author$
# $Revision$
# $Id$
#------------------------------------------------------------

#------------------------------------------------------------
# Makefile for sqg program
# always link appropriate makefile.conf and do "make clean" first
#------------------------------------------------------------

#------------------------------------------------------------
# Get F90, FFLAGS, NCFLAGS, RM from makefile.conf
include makefile.conf
#------------------------------------------------------------

#------------------------------------------------------------
SRC    = .
BIN    = .
#------------------------------------------------------------

#------------------------------------------------------------
common_OBJ = $(SRC)/fft.f90

sqg_OBJ = $(common_OBJ) \
          $(SRC)/sqg_mod.f90 \
          $(SRC)/sqg.f90

pvinv_OBJ = $(common_OBJ) \
            $(SRC)/pvinv_mod.f90 \
            $(SRC)/pvinv.f90

arBase_OBJ = $(SRC)/add_or_remove_base.f90

sqg_stats_OBJ = $(SRC)/sqg_stats.f90

sqg_diag_OBJ= $(SRC)/sqg_diag.f90
#------------------------------------------------------------

#------------------------------------------------------------
.f90.o :
	$(F90) $(FFLAGS) -c $< -o $@

all : sqg pvinv \
      add_or_remove_base \
      sqg_stats \
      sqg_diag

sqg : $(sqg_OBJ)
	$(F90) $(FFLAGS) -o $(BIN)/sqg $(sqg_OBJ) $(NCFLAGS)

pvinv : $(pvinv_OBJ)
	$(F90) $(FFLAGS) -o $(BIN)/pvinv $(pvinv_OBJ) $(NCFLAGS)

add_or_remove_base : $(arBase_OBJ)
	$(F90) $(FFLAGS) -o $(BIN)/add_or_remove_base $(arBase_OBJ) $(NCFLAGS)

sqg_stats : $(sqg_stats_OBJ)
	$(F90) $(FFLAGS) -o $(BIN)/sqg_stats $(sqg_stats_OBJ) $(NCFLAGS)

sqg_diag : $(sqg_diag_OBJ)
	$(F90) $(FFLAGS) -o $(BIN)/sqg_diag $(sqg_diag_OBJ) $(NCFLAGS)

clean :
	$(RM) $(SRC)/*.o $(RM) $(SRC)/*.mod \
$(RM) $(BIN)/sqg \
$(RM) $(BIN)/pvinv \
$(RM) $(BIN)/add_or_remove_base \
$(RM) $(BIN)/sqg_stats \
$(RM) $(BIN)/sqg_diag
#------------------------------------------------------------
