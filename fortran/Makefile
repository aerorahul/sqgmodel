#------------------------------------------------------------
# <next few lines under version control, D O  N O T  E D I T>
# $Date$
# $Author$
# $Revision$
# $Id$
#------------------------------------------------------------

#------------------------------------------------------------
# Makefile for sQG program
# always do "make clean" first.
#------------------------------------------------------------

#------------------------------------------------------------
machine := ${shell uname}
ifeq (${machine}, Linux)
NETCDF  = /usr
NCFLAGS = -L${NETCDF}/lib -lnetcdf -lnetcdff -I${NETCDF}/include
mFFLAGS =
endif
ifeq (${machine}, Darwin)
NETCDF  = /opt/netcdf
NCFLAGS = -L${NETCDF}/lib -lnetcdf -I${NETCDF}/include
mFFLAGS = -shared-intel
endif
#------------------------------------------------------------

#------------------------------------------------------------
F90    = ifort
FFLAGS = -g -check -traceback $(mFFLAGS) # Use these flags for development / testing
FFLAGS = -w -O3 -ftz -align all -fno-alias -fp-model precise -convert big_endian -unroll5 $(mFFLAGS) 
RM     = rm -f
SRC    = .
BIN    = .
#------------------------------------------------------------

#------------------------------------------------------------
sQG_OBJ = $(SRC)/SPECTRAL.f90 \
          $(SRC)/sqg_2b.f90   \
          $(SRC)/fft.f90

arBase_OBJ = $(SRC)/add_or_remove_base.f90

sQG_diag_OBJ = $(SRC)/sqg_diag.f90
#------------------------------------------------------------

#------------------------------------------------------------
.f90.o :
	$(F90) $(FFLAGS) -c $< -o $@

sQG : $(sQG_OBJ)
	$(F90) $(FFLAGS) -o $(BIN)/sQG $(sQG_OBJ) $(NCFLAGS)

add_or_remove_base : $(arBase_OBJ)
	$(F90) $(FFLAGS) -o $(BIN)/add_or_remove_base $(arBase_OBJ) $(NCFLAGS)

sQG_diag: $(sQG_diag_OBJ)
	$(F90) $(FFLAGS) -o $(BIN)/sQG_diag $(sQG_diag_OBJ) $(NCFLAGS)

all: sQG \
     add_or_remove_base \
     sQG_diag

clean:
	$(RM) $(SRC)/*.o $(RM) $(SRC)/*.mod \
$(RM) $(BIN)/sQG   \
$(RM) $(BIN)/add_or_remove_base \
$(RM) $(BIN)/sQG_diag
#------------------------------------------------------------
