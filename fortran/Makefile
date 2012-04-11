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
machine := ${shell uname -m}
host := ${shell hostname -s}

ifeq (${machine}, x86_64)
NETCDF = /usr/local/other/netcdf/4.1.1/serial_Intel-11.1.069 
NCFLAGS = -L${NETCDF}/lib -lnetcdff -lnetcdf -I${NETCDF}/include
mFFLAGS =
endif

ifeq (${host}, stratopause)
NETCDF = /usr
NCFLAGS = -L${NETCDF}/lib -lnetcdff -lnetcdf -I${NETCDF}/include
mFFLAGS =
endif

ifeq (${host}, trident)
NETCDF  = /ford1/local/EL6/netcdf4
NCFLAGS = -L${NETCDF}/lib -lnetcdff -lnetcdf -I${NETCDF}/include
mFFLAGS =
endif

ifeq (${host}, blue)
NETCDF  = /ford1/local/EL5/netcdf4
NCFLAGS = -L${NETCDF}/lib -lnetcdff -lnetcdf -I${NETCDF}/include
mFFLAGS =
endif
#------------------------------------------------------------

#------------------------------------------------------------
F90    = ifort
FFLAGS = -g -check -traceback $(mFFLAGS) # Use these flags for development / testing
FFLAGS = -w -O3 -ftz -align all -fno-alias -fp-model precise -convert big_endian -unroll5 -parallel $(mFFLAGS) 
RM     = rm -f
SRC    = .
BIN    = .
#------------------------------------------------------------

#------------------------------------------------------------
sQG_OBJ = $(SRC)/SPECTRAL.f90 \
          $(SRC)/sqg_2b.f90   \
          $(SRC)/fft.f90

arBase_OBJ = $(SRC)/add_or_remove_base.f90

sqg_stats_OBJ = $(SRC)/sqg_stats.f90

sqg_diag_OBJ= $(SRC)/sqg_diag.f90
#------------------------------------------------------------

#------------------------------------------------------------
.f90.o :
	$(F90) $(FFLAGS) -c $< -o $@

sQG : $(sQG_OBJ)
	$(F90) $(FFLAGS) -o $(BIN)/sQG $(sQG_OBJ) $(NCFLAGS)

add_or_remove_base : $(arBase_OBJ)
	$(F90) $(FFLAGS) -o $(BIN)/add_or_remove_base $(arBase_OBJ) $(NCFLAGS)

sqg_stats: $(sqg_stats_OBJ)
	$(F90) $(FFLAGS) -o $(BIN)/sqg_stats $(sqg_stats_OBJ) $(NCFLAGS)

sqg_diag: $(sqg_diag_OBJ)
	$(F90) $(FFLAGS) -o $(BIN)/sqg_diag $(sqg_diag_OBJ) $(NCFLAGS)

all: sQG \
     add_or_remove_base \
     sqg_stats \
     sqg_diag

clean:
	$(RM) $(SRC)/*.o $(RM) $(SRC)/*.mod \
$(RM) $(BIN)/sQG   \
$(RM) $(BIN)/add_or_remove_base \
$(RM) $(BIN)/sqg_stats \
$(RM) $(BIN)/sqg_diag
#------------------------------------------------------------
