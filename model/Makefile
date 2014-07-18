#------------------------------------------------------------
# <next few lines under version control, D O  N O T  E D I T>
# $Date$
# $Author$
# $Revision$
# $Id$
#------------------------------------------------------------

#------------------------------------------------------------
# Makefile for sqg program
# always do "make clean" first.
#------------------------------------------------------------

#------------------------------------------------------------
machine := ${shell uname -m}
host := ${shell hostname -s}

# Discover
ifeq (${machine}, x86_64)
NETCDF = /usr/local/other/netcdf/4.1.1/serial_Intel-11.1.069
HDF    = /usr/local/other/hdf5/1.8.4/serial_intel-11.1.069
ZLIB   = /usr/local/other/zlib/1.2.3-intel-11.1.069
SZIP   = /usr/local/other/szip/2.1/Intel-11.1.069
INCS = -I$(NETCDF)/include
LIBS = -L$(NETCDF)/lib -lnetcdf \
       -L$(HDF)/lib -lhdf5_hl -lhdf5  \
       -L$(ZLIB)/lib -lz -lm  \
       -L$(SZIP)/lib -lsz
mFFLAGS =
endif

#stratopause, dorne, valyria
ifeq (${host}, stratopause)
NETCDF = /usr
LIBS = ${shell $(NETCDF)/bin/nc-config --flibs}
INCS = -I$(NETCDF)/include/netcdf
mFFLAGS =
else
ifeq (${host}, dorne)
NETCDF  = /opt/Baselibs/v3.2.1/i686-apple-darwin12.0.0/ifort/Darwin
LIBS = ${shell $(NETCDF)/bin/nc-config --flibs}
INCS = -I$(NETCDF)/include/netcdf
mFFLAGS = -shared-intel
else
ifeq (${host}, valyria)
NETCDF  = /opt/Baselibs/v3.2.1/i686-apple-darwin12.0.0/ifort/Darwin
LIBS = ${shell $(NETCDF)/bin/nc-config --flibs}
INCS = -I$(NETCDF)/include/netcdf
mFFLAGS = -shared-intel
endif
endif
endif

NCFLAGS = $(LIBS) $(INCS)
#------------------------------------------------------------

#------------------------------------------------------------
F90    = ifort
FFLAGS = -w -O3 -ftz -align all -fno-alias -fp-model precise -cm -axT -u -convert big_endian -parallel $(mFFLAGS)
# Use these flags for development / testing
#FFLAGS = -O0 -g -check all -ftz -vec-report0 -fp-model precise -ftrapuv -traceback $(mFFLAGS)
RM     = rm -f
SRC    = .
BIN    = .
#------------------------------------------------------------

#------------------------------------------------------------
common_OBJ = $(SRC)/spectral_mod.f90 \
             $(SRC)/fft.f90

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
