%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% <next few lines under version control, D O  N O T  E D I T>
% $Date$
% $Author$
% $Revision$
% $Id$
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% function compareprofiles - compares localization profiles 
%                            in theta and phi.
%                            A cosine bell is used to i
%                            specify the function.
%
%   compareprofiles
% 
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function compareprofiles

Nx = 128; Ny = 64;
Lx =  20; Ly = 2*5.539118; H=1.0;

%  set coordinates (x,y)
xx = 0:Lx/Nx:Lx;   x = xx(2:Nx+1) - Lx/2;
yy = 0:Ly/Ny:Ly;   y = yy(2:Ny+1) - Ly/2;
[xg, yg] = meshgrid(x,y);

%  fourier wavenumber operators
facx = 2*pi/Lx; 
facy = 2*pi/Ly;

dx = [-Nx/2:Nx/2-1] * (2*pi/Lx);
dy = [-Ny/2:Ny/2-1] * (2*pi/Ly);
[DX DY] = meshgrid(dx,dy);
DX = fftshift(DX);   DY = fftshift(DY);

%	Inversion Parameters
m = sqrt(DX.^2+DY.^2);
m(1,1) = 1; 
IZ = real(1./(m.*tanh(m*H))); 
IZo = real(1./(m.*sinh(m*H))); 
m(1,1) = 0;
IZ(1,1) = 1; IZo(1,1) = 1; 
DZ = 1./IZ; DZo = 1./IZo; 
DZ(1,1) = 0; DZo(1,1) = 0; IZ(1,1) = 0; IZo(1,1) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LW = 'linewidth';
lw = 2.0;
FS = 'fontsize';
fs = 14;
FW = 'fontweight';
fw = 'bold';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x0=0.0; y0=0.0;
a = 2.0; b = a;
amp=16.0;
arg = sqrt(((xg-x0)/4*a).^2 + ((yg-y0)/4*b).^2);
fxy = ((amp/16).*(1+cos(pi.*arg)).^4).*(arg<=1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Localized theta
thT_xy = fxy;
thB_xy = zeros(Ny,Nx);
thT_sp = fft2(thT_xy);
thB_sp = fft2(thB_xy);

% lower boundary
phiB_sp = -IZ.*thB_sp + IZo.*thT_sp; % spectral phi
uB      = -real(ifft2(i*DY.*phiB_sp));
vB      =  real(ifft2(i*DX.*phiB_sp));
phiB_xy =  real(ifft2(phiB_sp));

% upper boundary
phiT_sp = -IZo.*thB_sp + IZ.*thT_sp; % spectral phi
uT   = -real(ifft2(i*DY.*phiT_sp));
vT   =  real(ifft2(i*DX.*phiT_sp));
phiT_xy =  real(ifft2(phiT_sp));

maxval = max(max(abs(thT_xy)));

figure(1); clf;
plot(x,thT_xy(32,:)./maxval,'r',  ...
     x,phiT_xy(32,:)./maxval,'b', ... 
		 x,vT(32,:)./maxval,'k',      ...
		 LW,lw);
xlabel('X',FS,fs,FW,fw)
ylabel('Amplitude',FS,fs,FW,fw)

l = legend('\theta','v','\phi');

axis([-Lx/2 Lx/2 -0.6 1.1]);
h2 = gca;
set(h2,'Fontweight','bold','Fontsize',12);

t = title('Localized Theta');

set(l,FS,fs,FW,fw);
set(t,FS,fs,FW,fw);

savefigure(gcf,'loctheta');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Localized phi 
phiT_xy = fxy;
phiT_sp = fft2(phiT_xy);
thT_sp = DZ.*phiT_sp;
thT_xy = real(ifft2(thT_sp));

uT = -real(ifft2(i*DY.*phiT_sp));
vT =  real(ifft2(i*DX.*phiT_sp));

maxval = max(max(abs(thT_xy)));

figure(2); clf;
plot(x,thT_xy(32,:)./maxval,'r',  ...
     x,phiT_xy(32,:)./maxval,'b', ... 
		 x,vT(32,:)./maxval,'k',      ...
		 LW,lw);
xlabel('X',FS,fs,FW,fw)
ylabel('Amplitude',FS,fs,FW,fw)

l = legend('\theta','v','\phi');

axis([-Lx/2 Lx/2 -0.6 1.1]);
h2 = gca;
set(h2,'Fontweight','bold','Fontsize',12);

t = title('Localized Geopotential');

set(l,FS,fs,FW,fw);
set(t,FS,fs,FW,fw);

savefigure(gcf,'locphi');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
