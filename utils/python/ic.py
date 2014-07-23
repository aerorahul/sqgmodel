#!/usr/bin/env python

#############################################################
# <next few lines under version control, D O  N O T  E D I T>
# $Date$
# $Author$
# $Revision$
# $Id$
#############################################################
#
#############################################################
#  2 boundary initial conditions
#############################################################

import sys
import numpy as np
from netCDF4    import Dataset
from matplotlib import pyplot,cm
from argparse   import ArgumentParser, ArgumentDefaultsHelpFormatter

parser = ArgumentParser(description = 'Generate initial conditions', formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('-i','--initial_condition',help='initial condition, {%(choices)s}',type=int,choices=np.arange(1,12),required=False,default=1,metavar='1 = turb; 2 = tilted vorts; 3 = baroclinic instability; 4 = Gaussian; 5 = 2 vorts; 6 = edge wave; 7 = topo wave; 8 = cosine bell; 11 = turb w/ no baro flow!')
parser.add_argument('-nx','--Nx',help='no. of X grid-points',type=int,choices=[32,64,128,256,512],required=False, default=128)
parser.add_argument('-ny','--Ny',help='no. of Y grid-points',type=int,choices=[32,64,128,256,512],required=False, default=64)
parser.add_argument('-Ro','--Rossby',help='Rossby number',type=float,choices=[0.0, 0.1, 0.001],required=False,default=0.0)
args = parser.parse_args()

# Initial condition
Nic = args.initial_condition
# 1 = turb; 2 = tilted vorts; 3 = baroclinic instability; 4 =
# Gaussian, 11 = turb w/ no baro flow!, 5 = 2 vorts, 6 =
# edge wave, 7 = topo wave, 8 = cosine bell

# number of gridpoints
Nx = args.Nx
Ny = args.Ny

# Rossby number
Ro = args.Rossby

# domain size
Lx =  20; Ly = 2*5.539118; H = 1.; # baroclinic instability
#Lx = 5*8*np.pi; Ly = 5*4*np.pi; H = 6; # Hsqg waves
#Lx = 4*np.pi; Ly = 2*np.pi; H = 1.0; # Hsqg waves
#Lx = 7*np.pi; Ly = 7*np.pi; H = 1; # original turbulence runs.
#Lx = 7*np.pi; Ly = 7*np.pi; H = 2; # original turbulence runs.
#Lx = 2*np.pi*14; Ly = 2*np.pi*14; H = 1; # original turbulence runs.
#Lx =    4*np.pi; Ly =   4*np.pi; H = 1.; # 128 turbulence (dx=dy=.0982 ~100km).
#Lx =  3.9118; Ly = 5.539118; H = 1.; # baroclinic instability
#Lx =  4; Ly = 5.539118; H = 1.; # baroclinic instability
#Lx =    2*np.pi; Ly =   2*np.pi; H = 1.; #  tilted vortex.

# vortex amplitude
amp = -1

# set coordinates (x,y)
xx = np.linspace(0.0,Lx,Nx+1)
yy = np.linspace(0.0,Ly,Ny+1)
x = xx[1:] - Lx/2.0
y = yy[1:] - Ly/2.0
xg, yg = np.meshgrid(x,y)

# fourier wavenumber operators
facx = 2*np.pi/Lx
facy = 2*np.pi/Ly

dx = np.linspace(-Nx/2,Nx/2,Nx,endpoint=False) * facx
dy = np.linspace(-Ny/2,Ny/2,Ny,endpoint=False) * facy
DX, DY = np.meshgrid(dx,dy)
DX = np.fft.fftshift(DX)
DY = np.fft.fftshift(DY)

#	Inversion Parameters
if ( H > 10 ):
    print '***  sQG inversion ***'
    DZ       = -np.sqrt(DX**2+DY**2)
    DZi      = DZ
    DZi[0,0] = 1
    DZi      = 1.0 / DZi
    DZi[0,0] = 0
else:
    print '*** 2sQG inversion ***'
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

if   ( Nic == 1 or Nic == 11 ):

    print 'turbulence ICs'

    # specify the peak wavenumber
    #mk0 = 14;     mm = 25  # original
    mk0 = 28;     mm = 25  # original
    sigma = 0.5 # gaussian spread factor
    k0 = mk0 * (np.pi/Lx)

    # magnitude of horizontal wavenumber:
    MK = np.sqrt(DX**2+DY**2)

    # define phi^0 (zero nyquist modes):
    # polvani formula:
    theta0Te = (MK**(mm/4 - 1)) / ((MK + k0)**(mm/2))

    # upper boundary
    np.random.seed()

    a1 = np.random.rand(Ny/2-1,Nx/2-1)
    b1 = np.random.rand(Ny/2-1,1)
    a2 = np.random.rand(Ny/2-1,Nx/2-1)
    b2 = np.random.rand(1,Nx/2-1)

    t1 =  np.hstack([a1,b1,a2])
    t2 =  np.hstack([b2, np.zeros((1,1)), -np.fliplr(b2)])
    t3 = np.hstack([-np.flipud(np.fliplr(a2)), -np.flipud(b1), -np.flipud(np.fliplr(a1))])
    c1 = np.vstack([t1,t2,t3])
    c2 = np.fft.fftshift(np.vstack([np.zeros((1,Nx)), np.hstack([np.zeros((Ny-1,1)),c1])]))

    thetaB    = theta0Te * np.exp(1j*2.0*np.pi*c2) # spectral theta
    thetaB_xy = np.real(np.fft.ifft2(thetaB))

    # upper boundary
    np.random.seed()

    a1 = np.random.rand(Ny/2-1,Nx/2-1)
    b1 = np.random.rand(Ny/2-1,1)
    a2 = np.random.rand(Ny/2-1,Nx/2-1)
    b2 = np.random.rand(1,Nx/2-1)

    t1 =  np.hstack([a1,b1,a2])
    t2 =  np.hstack([b2, np.zeros((1,1)), -np.fliplr(b2)])
    t3 = np.hstack([-np.flipud(np.fliplr(a2)), -np.flipud(b1), -np.flipud(np.fliplr(a1))])
    c1 = np.vstack([t1,t2,t3])
    c2 = np.fft.fftshift(np.vstack([np.zeros((1,Nx)), np.hstack([np.zeros((Ny-1,1)),c1])]))

    Hnorm = 1

    # option for no initial barotropic flow:
    if   ( Nic == 1  ):
        thetaT = thetaB
    elif ( Nic == 11 ):
        thetaT  = theta0Te * np.exp(1j*2.0*np.pi*c2) # spectral theta

    thetaT_xy = np.real(np.fft.ifft2(thetaT))

elif ( Nic == 2 ):

    print 'tilted vortex ICs'

    #### BOTTOM VORTEX
    #  snyder parameters
    a = 2 ; b = a # scale control
    c = 2; l = np.sqrt(2)
    #x0 = -.5; y0 = 0 #vortex origin
    x0 = 0; y0 = 0 #vortex origin

    # "Plateau" vortex:
    arg   = (b*(xg-x0))**2 + (a*(yg-y0))**2
    theta  = amp*(1.0 - np.tanh( (arg -c)/l**2 ))
    thetaB = 0.0*np.fft.fft2(theta) # spectral theta

    thetaB_xy = 0.0*theta

    #### TOPPOM VORTEX
    x0 = 1/a; y0 = 0 # vortex origin

    # "Plateau" vortex:
    arg   = (b*(xg-x0))**2 + (a*(yg-y0))**2
    theta  = amp*(1.0 - np.tanh( (arg -c)/l**2 ))
    thetaT = np.fft.fft2(theta) # spectral theta

    thetaT_xy = theta

    Hnorm = 1

elif ( Nic == 3 ):

    print 'baroclinic instability ICs'

    ang = 0.83775
    amp = -2.2

    thetaB = np.fft.fft2(amp*np.sin((2.0*np.pi*xg/Lx) + ang))
    thetaT = np.fft.fft2(amp*np.sin((2.0*np.pi*xg/Lx)))

    thetaB_xy = amp*np.sin((2.0*np.pi*xg/Lx) + ang)
    thetaT_xy = amp*np.sin((2.0*np.pi*xg/Lx))

    Hnorm = 0

elif ( Nic == 4 ):

    print 'Gaussian upper vort ICs'

    amp = -5.0
    a = 0.5; b = a; # scale control
    x0 = 0; y0 = 0; # vortex origin

    arg = (((xg-x0)/a)**2 + ((yg-y0)/b)**2)/2

    theta  = amp*np.exp(-arg)
    thetaT = np.fft.fft2(theta) # spectral theta
    thetaB = np.zeros((Ny,Nx))

    zetaBT = np.sum(np.abs(thetaT)**2)

    zetaT  = -( (H/2*m**2)*(1+np.cosh(H*m))*IZo ) * thetaT
    zetaT  = np.sum(np.abs(zetaT)**2)

    Bind = zetaBT/zetaT

    thetaB_xy = thetaB
    thetaT_xy = theta

    Hnorm = 2

elif ( Nic == 5 ):

    print '2 vorticies put side by side'

    #amp = -1.; # cyclone (cool air mass in middle)
    amp = 1.; # anticyclone (warm air mass in middle)
    a= 2; b= a; # scale control
    x0 = 1; y0 = 0; # vortex to the right (play with x0)
    x1= -1; # vortex to the left (play with x0)

    arg  = (b*(xg-x0))**2 + (a*(yg-y0))**2
    arg1 = (b*(xg-x1))**2 + (a*(yg-y0))**2

    theta1 = amp*np.exp(-arg1)
    theta  = amp*np.exp(-arg) + theta1

    thetaT = np.fft.fft2(theta) # spectral theta

    thetaB = np.zeros((Ny,Nx))

    thetaB_xy = thetaB
    thetaT_xy = theta

    Hnorm = 0

elif ( Nic == 6 ):

    print 'edge wave w/ vort-norm'

    theta = np.cos(xg) * np.cos(yg)

    thetaT = np.fft.fft2(theta)
    thetaB = np.zeros((Ny,Nx))

    thetaB_xy = thetaB
    thetaT_xy = theta

    Hnorm = 2

elif ( Nic == 7 ):

    print 'topo wave w/ vort-norm'

    #  arg   = (xg**2 + yg**2)/2
    #  theta = -np.exp(-arg)
    theta  = np.cos(xg) * np.sin(yg)
    thetaB = np.fft.fft2(theta)

    coeff  = np.cosh(H*m) - (H*m*np.sinh(H*m))

    thetaT = thetaB / coeff

    thetaB_xy = theta
    thetaT_xy = theta / coeff

    Hnorm = 2

elif ( Nic == 8 ):

    print 'Cosine bell upper vort ICs'

    amp = 0.175
    a= 1.75; b= a  # vortex scale control
    x0 = 0; y0 = 0 # vortex origin

    arg   = np.sqrt(((xg-x0)/4*a)**2 + ((yg-y0)/4*b)**2)
    theta  = ((amp/16)*(1.0+np.cos(np.pi*arg))**4)*(arg<=1)

    thetaT = np.fft.fft2(theta) # spectral theta
    thetaB = np.zeros((Ny,Nx))

    zetaBT = np.sum(np.abs(thetaT)**2)

    zetaT  = -( (H/2*m**2)*(1.0+np.cosh(H*m))*IZo )*thetaT
    zetaT  = np.sum(np.abs(zetaT)**2)

    Bind = zetaBT/zetaT

    thetaB_xy = thetaB
    thetaT_xy = theta

    Hnorm = 2

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
    th_B = thetaB_xy
    th_T = thetaT_xy
elif ( Hnorm == 1 ):
    th_B = np.real(np.fft.ifft2(thetaB/Vnorm))
    th_T = np.real(np.fft.ifft2(thetaT/Vnorm))
    th_B = thetaB_xy / Vnorm
    th_T = thetaT_xy / Vnorm
elif ( Hnorm == 2 ):
    th_B = np.abs(amp) * np.real(np.fft.ifft2(thetaB/Znorm))
    th_T = np.abs(amp) * np.real(np.fft.ifft2(thetaT/Znorm))
    th_B = (np.abs(amp)/Znorm) * thetaB_xy
    th_T = (np.abs(amp)/Znorm) * thetaT_xy

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
