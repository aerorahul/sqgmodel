#!/usr/bin/env python

#############################################################
# <next few lines under version control, D O  N O T  E D I T>
# $Date$
# $Author$
# $Revision$
# $Id$
#############################################################

#############################################################
#  ic_tracer - create tracer field of unit magnitude
#              will write out tr_init.nc
#############################################################

import sys
import numpy as np
from netCDF4 import Dataset
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

parser = ArgumentParser(description = 'Generate initial conditions for tracer', formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('-f','--filename',help='theta filename',type=str,required=False,default='th_init.nc')
args = parser.parse_args()

fname = args.filename

nc = Dataset(fname,mode='r')
theta = nc.variables['theta']
thB = theta[0,:]
thT = theta[1,:]

Nx,Ny,Nz = np.shape(theta)

maxval = np.max(np.max(np.abs(thB)))
trB = np.abs(thB / maxval)
maxval = np.max(np.max(np.abs(thT)))
trT = np.abs(thT / maxval)

tracer = np.zeros((Nz,Ny,Nx))
tracer[0,:] = trB
tracer[1,:] = trT
nc = Dataset('tr_init.nc',mode='w',clobber=True)
Dim = nc.createDimension('nx',size=Nx)
Dim = nc.createDimension('ny',size=Ny)
Dim = nc.createDimension('nz',size=2)
Var = nc.createVariable('tracer','f8',('nz','ny','nx',))
nc.close()

nc = Dataset('tr_init.nc',mode='a',clobber=True)
nc.variables['tracer'][:] = tracer
nc.close()
