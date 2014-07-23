#!/usr/bin/env python

#############################################################
# <next few lines under version control, D O  N O T  E D I T>
# $Date$
# $Author$
# $Revision$
# $Id$
#############################################################

import sys
import numpy as np
from netCDF4    import Dataset
from matplotlib import pyplot,cm
from argparse   import ArgumentParser, ArgumentDefaultsHelpFormatter

parser = ArgumentParser(description = 'Generate initial conditions', formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('-l','--localization',help='localization function',type=int,choices=[1,2],required=False, default=1)
parser.add_argument('-i','--initial_condition',help='initial condition',type=int,choices=[1,2],required=False, default=1)
args = parser.parse_args()

# Localization function
# 1 = Gaussian
# 2 = Cosine Bell
# 3 = Gaspari-Cohn
localization = args.localization

# Initial condition
# 1 = Neumann IC (theta)
# 2 = Dirchilet IC (phi)
init_condition = args.initial_condition

# Offset the initial disturbance from y = 0
offset = 0

# Domain size and Grid points
N  = 128                       # Number of gridpoints
Nx = N                         # Grid points in X
Ny = N/2                       # Grid points in Y
Lx =  20.0                     # Domain length in X
Ly = 2*5.539118                # Domain length in Y
H = 1.0                        # baroclinic instability
Ro = 0.0                       # Rossby number
Hnorm = 1                      # Norm to normalize (0:None, 1:Velocity, 2:Vorticity)

#  set coordinates (x,y)
xx = np.linspace(0.0,Lx,Nx+1)
yy = np.linspace(0.0,Ly,Ny+1)
x = xx[1:] - Lx/2.0
y = yy[1:] - Ly/2.0
xg, yg = np.meshgrid(x,y)

#  fourier wavenumber operators
facx = 2*np.pi/Lx
facy = 2*np.pi/Ly

dx = np.linspace(-Nx/2,Nx/2,Nx,endpoint=False) * facx
dy = np.linspace(-Ny/2,Ny/2,Ny,endpoint=False) * facy
DX, DY = np.meshgrid(dx,dy)
DX = np.fft.fftshift(DX)
DY = np.fft.fftshift(DY)

#	Inversion Parameters
m        = np.sqrt(DX**2+DY**2)
m[0,0]   = 1.0
IZ       = 1.0 / (m*np.tanh(m*H))
IZo      = 1.0 / (m*np.sinh(m*H))
m[0,0]   = 0.0
IZ[0,0]  = 1.0
IZo[0,0] = 1.0
DZ       = 1.0 / IZ
DZo      = 1.0 / IZo
DZ[0,0]  = 0.0
DZo[0,0] = 0.0
IZ[0,0]  = 0.0
IZo[0,0] = 0.0

if ( localization == 1 ):

    print 'Gaussian bell upper vort ICs'

    amp = -5.0                        # vortex amplitude
    a  = 0.5 ; b = a                  # vortex scale
    x0 = 0.0 ; y0 = y[Ny/2-1-offset]  # vortex origin

    arg = (((xg-x0)/a)**2 + ((yg-y0)/b)**2) / 2.0
    fxy = amp * np.exp(-arg)

elif ( localization == 2 ):

    print 'Cosine bell upper vort ICs'

    amp = 0.175                         # vortex amplitude
    amp = -0.51375                      # vortex amplitude
    amp = -1.0                          # vortex amplitude
    a  = 1.75  ;  b  = a;               # vortex scale
    a  = 3.5   ;  b  = a;               # vortex scale
    x0 = 0.0 ; y0 = y[Ny/2-1-offset]    # vortex origin

    arg = np.sqrt(((xg-x0)/4.0*a)**2 + ((yg-y0)/4.0*b)**2)
    fxy = ((amp/16.0)*(1.0+np.cos(np.pi*arg))**4)*(arg<=1)

if   ( init_condition == 1 ):

    print 'Neumann ICs'

    thetaT = np.fft.fft2(fxy)     # spectral theta at Tropopause
    thetaB = np.zeros((Ny,Nx))    # spectral theta at Ground

elif ( init_condition == 2 ):

    print 'Dirichlet ICs'

    thetaT = DZ*np.fft.fft2(fxy)  # spectral theta at Tropopause
    thetaB = np.zeros((Ny,Nx))    # spectral theta at Ground

###################################################
# normalize so that max leading-order wind is O(1).
###################################################
# lower boundary
Phi0B = -IZ*thetaB + IZo*thetaT             # spectral phi
Phi0Bx = np.real(np.fft.ifft2(1j*DX*Phi0B))
Phi0By = np.real(np.fft.ifft2(1j*DY*Phi0B))
vB = -Phi0Bx
uB =  Phi0By
magVB = ((vB**2) + (uB**2))**0.5

# upper boundary
Phi0T = -IZo*thetaB + IZ*thetaT             # spectral phi
Phi0Tx = np.real(np.fft.ifft2(1j*DX*Phi0T))
Phi0Ty = np.real(np.fft.ifft2(1j*DY*Phi0T))
vT = -Phi0Tx
uT =  Phi0Ty
magVT = ((vT**2) + (uT**2))**0.5

#   velocity norm
Vnorm = np.max([np.max(np.max(magVB)), np.max(np.max(magVT))])

#   vorticity norm
Zeta0T = np.real(np.fft.ifft2(-(m**2)*(Phi0T)))
Zeta0B = np.real(np.fft.ifft2(-(m**2)*(Phi0B)))
Znorm  = np.max([np.max(Zeta0T), np.max(Zeta0B)])

if   ( Hnorm == 0 ):
    th_B = np.real(np.fft.ifft2(thetaB))
    th_T = np.real(np.fft.ifft2(thetaT))
elif ( Hnorm == 1 ):
    th_B = np.real(np.fft.ifft2(thetaB/Vnorm))
    th_T = np.real(np.fft.ifft2(thetaT/Vnorm))
elif ( Hnorm == 2 ):
    th_B = np.abs(amp) * np.real(np.fft.ifft2(thetaB/Znorm))
    th_T = np.abs(amp) * np.real(np.fft.ifft2(thetaT/Znorm))

# zero out theta beyond tolerance
tol = 1.0e-6
th_B = th_B * (np.abs(th_B) > tol)
th_T = th_T * (np.abs(th_T) > tol)

#   plot figures
cmap=cm.get_cmap(name='PuOr_r',lut=16)
f1 = pyplot.figure(1)
pyplot.clf()
pyplot.pcolor(x,y,th_B,cmap=cmap)
t = 'Lower boundary potential temperature'
pyplot.title(t,fontsize=14,fontweight='bold')
pyplot.colorbar()
cmax = np.max(np.max(np.abs(th_B)))
pyplot.clim(-cmax,cmax)
pyplot.axis('image')
f1.canvas.set_window_title(t)

f2 = pyplot.figure(2)
pyplot.clf()
pyplot.pcolor(x,y,th_T,cmap=cmap)
t = 'Upper boundary potential temperature'
pyplot.title(t,fontsize=14,fontweight='bold')
pyplot.colorbar()
cmax = np.max(np.max(np.abs(th_T)))
pyplot.clim(-cmax,cmax)
pyplot.axis('image')
f2.canvas.set_window_title(t)

f3 = pyplot.figure(3)
pyplot.clf()
pyplot.pcolor(x,y,th_T-th_B,cmap=cmap)
t = 'Upper -Lower boundary potential temperature'
pyplot.title(t,fontsize=14,fontweight='bold')
pyplot.colorbar()
cmax = np.max(np.max(np.abs(th_T-th_B)))
pyplot.clim(-cmax,cmax)
pyplot.axis('image')
f3.canvas.set_window_title(t)

f4 = pyplot.figure(4)
pyplot.clf()
pyplot.plot(x,th_T[Ny/2-1,:],'-r.');
t = 'X - Cross Section'
pyplot.title(t,fontsize=14,fontweight='bold')
f4.canvas.set_window_title(t)

f5 = pyplot.figure(5)
pyplot.clf()
pyplot.plot(y,th_T[:,Nx/2-1],'-r.');
t = 'Y - Cross Section'
pyplot.title(t,fontsize=14,fontweight='bold')
f5.canvas.set_window_title(t)

#	dump to disk:
theta = np.zeros((1,2,Ny,Nx))
theta[0,0,:] = th_B
theta[0,1,:] = th_T
nc = Dataset('th_init.nc',mode='w',clobber=True)
Dim = nc.createDimension('nx',  size=Nx)
Dim = nc.createDimension('ny',  size=Ny)
Dim = nc.createDimension('nz',  size=2)
Dim = nc.createDimension('time',size=1)
Var = nc.createVariable('theta','f8',('time','nz','ny','nx',))
nc.close()

nc = Dataset('th_init.nc',mode='a',clobber=True)
nc.variables['theta'][:] = theta
nc.close()

# determine max flow induced on home & opposing boundaries
thetaB = np.fft.fft2(th_B)
thetaT = np.fft.fft2(th_T)

# lower boundary
Phi0B  = -IZ*thetaB # spectral phi
Phi0Bx = np.real(np.fft.ifft2(1j*DX*Phi0B))
Phi0By = np.real(np.fft.ifft2(1j*DY*Phi0B))
magVBb = ((Phi0Bx**2) + (Phi0By**2))**0.5; # bottom on bottom

Phi0B  = IZo*thetaT # spectral phi
Phi0Bx = np.real(np.fft.ifft2(1j*DX*Phi0B));
Phi0By = np.real(np.fft.ifft2(1j*DY*Phi0B));
magVBt = ((Phi0Bx**2) + (Phi0By**2))**0.5; # bottom on top

VBb = np.max(np.max(magVBb))
VBt = np.max(np.max(magVBt))
print 'Max lower flow on lower, upper, and ratio: %6.3f %6.3f %6.3f' % (VBb, VBt, VBt/VBb)

# upper boundary
Phi0T  = IZ*thetaT # spectral phi
Phi0Tx = np.real(np.fft.ifft2(1j*DX*Phi0T))
Phi0Ty = np.real(np.fft.ifft2(1j*DY*Phi0T))
magVTt = ((Phi0Tx**2) + (Phi0Ty**2))**0.5 # top on top

Phi0T  = -IZo*thetaB # spectral phi
Phi0Tx = np.real(np.fft.ifft2(1j*DX*Phi0T))
Phi0Ty = np.real(np.fft.ifft2(1j*DY*Phi0T))
magVTb = ((Phi0Tx**2) + (Phi0Ty**2))**0.5 # top on top

VTt = np.max(np.max(magVTt))
VTb = np.max(np.max(magVTb))
print 'Max upper flow on lower, upper, and ratio: %6.3f %6.3f %6.3f' % (VTt, VTb, VTb/VTt)

pyplot.show()

sys.exit(0)
