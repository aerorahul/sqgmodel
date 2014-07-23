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
theta = nc.variables['theta'][:]
nc.close()
Nt,Nz,Ny,Nx = theta.shape
tracer = np.zeros((1,Nz,Ny,Nx))

for i in range(Nz):
    tracer[0,i,:] = np.abs(theta[0,i,:] / np.max(np.max(np.abs(theta[0,i,:]))))

nc = Dataset('tr_init.nc',mode='w',clobber=True)
Dim = nc.createDimension('nx',  size=Nx)
Dim = nc.createDimension('ny',  size=Ny)
Dim = nc.createDimension('nz',  size=2)
Dim = nc.createDimension('time',size=1)
Var = nc.createVariable('tracer','f8',('nz','ny','nx',))
nc.close()

nc = Dataset('tr_init.nc',mode='a',clobber=True)
nc.variables['tracer'][:] = tracer
nc.close()
