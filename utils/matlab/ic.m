%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% <next few lines under version control, D O  N O T  E D I T>
% $Date$
% $Author$
% $Revision$
% $Id$
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
%
%  2 boundary initial conditions;  17 july 02
%     - vorticity norm - 08 feb 2005 -- djm
%
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

Nic = 8; %(1 = turb; 2 = tilted vorts; 3 = baroclinic instability; 4 =
         %Gaussian, 11 = turb w/ no baro flow!, 5 = 2 vorts, 6 =
         %edge wave, 7 = topo wave, 8 = cosine bell)

N  = 128; % number of gridpoints
%N  = 32; % number of gridpoints
%Nx = N; Ny = N;
%Nx = N; Ny = N/2;
%Nx = 512; Ny = 128;
%Nx = 256; Ny = 64;
%Nx = 128; Ny = 64;
Nx = N; Ny = N/2;

Ro = 0.0; % Rossby number

% domain size
Lx =  20; Ly = 2*5.539118; H = 1.; % baroclinic instability
%Lx = 5*8*pi; Ly = 5*4*pi; H = 6; % Hsqg waves
%Lx = 4*pi; Ly = 2*pi; H = 1.0; % Hsqg waves
%Lx = 7*pi; Ly = 7*pi; H = 1; % original turbulence runs.
%Lx = 7*pi; Ly = 7*pi; H = 2; % original turbulence runs.
%Lx = 2*pi*14; Ly = 2*pi*14; H = 1; % original turbulence runs.
%Lx =    4*pi; Ly =   4*pi; H = 1.; % 128 turbulence (dx=dy=.0982 ~100km).
%Lx =  3.9118; Ly = 5.539118; H = 1.; % baroclinic instability
%Lx =  4; Ly = 5.539118; H = 1.; % baroclinic instability
%Lx =    2*pi; Ly =   2*pi; H = 1.; %  tilted vortex.

amp = -1; % vortex amplitude
%amp = 0.1; % vortex amplitude

facx = 2*pi/Lx; facy = 2*pi/Ly;

%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

%  set coordinates (x,y)

xx = 0:Lx/Nx:Lx;   x = xx(2:Nx+1) - Lx/2;
yy = 0:Ly/Ny:Ly;   y = yy(2:Ny+1) - Ly/2;

[xg, yg] = meshgrid(x,y);

%  fourier wavenumber operators

% fixed 03/03
dx = [-Nx/2:Nx/2-1] * (2*pi/Lx);
dy = [-Ny/2:Ny/2-1] * (2*pi/Ly);

[DX DY] = meshgrid(dx,dy);

DX = fftshift(DX);   DY = fftshift(DY);

if H > 10
  disp('*** sQG inversion ***');
  DZ = -sqrt(DX.^2+DY.^2);
  DZi = DZ;  DZi(1,1) = 1;  DZi = 1./DZi;  DZi(1,1) = 0;
else
  disp('*** 2sQG inversion ***');
  m = sqrt(DX.^2+DY.^2);
  m(1,1) = 1; 
  IZ = real(1./(m.*tanh(m*H))); 
  IZo = real(1./(m.*sinh(m*H))); 
  m(1,1) = 0;
  IZ(1,1) = 1; IZo(1,1) = 1; 
  DZ = 1./IZ; DZo = 1./IZo; 
  DZ(1,1) = 0; DZo(1,1) = 0; IZ(1,1) = 0; IZo(1,1) = 0;
end

switch(Nic)
  case{1,11}
   disp('  turbulence ICs')

   % specify the peak wavenumber
   %mk0 = 14;     mm = 25;  % original
   mk0 = 28;     mm = 25;  % original
   sigma = .5; % gaussian spread factor
   k0 = mk0 * (pi/Lx);

   % magnitude of horizontal wavenumber:
   MK = sqrt(DX.^2+DY.^2);

   % define phi^0 (zero nyquist modes):
   % polvani formula:
   theta0Te = (MK.^(mm/4 - 1))./((MK + k0).^(mm/2));

   % upper boundary
   rand('state',sum(100*clock));

   a1 = rand(Ny/2-1,Nx/2-1);  b1 = rand(Ny/2-1,     1);
   a2 = rand(Ny/2-1,Nx/2-1);  b2 = rand(     1,Nx/2-1);

   c1 = [a1 b1 a2 ; b2 0. -b2(1,end:-1:1);
	 -a2(end:-1:1,end:-1:1) -b1(end:-1:1,1) -a1(end:-1:1,end:-1:1)];
   c2 = fftshift([0 zeros(1,Nx-1) ; zeros(Ny-1,1) c1]);

   thetaB  = theta0Te .* exp(i*2*pi*c2); % spectral theta
					 % upper boundary
   rand('state',sum(1000*clock));

   a1 = rand(Ny/2-1,Nx/2-1);  b1 = rand(Ny/2-1,     1);
   a2 = rand(Ny/2-1,Nx/2-1);  b2 = rand(     1,Nx/2-1);

   c1 = [a1 b1 a2 ; b2 0. -b2(1,end:-1:1);
	 -a2(end:-1:1,end:-1:1) -b1(end:-1:1,1) -a1(end:-1:1,end:-1:1)];
   c2 = fftshift([0 zeros(1,Nx-1) ; zeros(Ny-1,1) c1]);

   Hnorm = 1;

   % option for no initial barotropic flow:
   switch(Nic)
    case{11}
     thetaT = thetaB;
    case{11}
     thetaT  = theta0Te .* exp(i*2*pi*c2); % spectral theta
   end
end

switch(Nic)
 case{2}
  disp('  tilted vortex ICs')

  %%%% BOTTOM VORTEX
  %  snyder parameters (see ic.m for details)
  a= 2; b= a; % scale control
  c = 2; l = sqrt(2);
  %x0 = -.5; y0 = 0; %vortex origin
  x0 = 0; y0 = 0; %vortex origin

  % "Plateau" vortex:
  arg   = (b*(xg-x0)).^2 + (a*(yg-y0)).^2;
  theta  = amp*(1 -tanh( (arg -c)/l^2 ));
  thetaB = 0.*fft2(theta); % spectral theta

	thetaB_xy = 0.*theta;

  %%%% TOPPOM VORTEX
  x0 = 1/a; y0 = 0; %vortex origin

  % "Plateau" vortex:
  arg   = (b*(xg-x0)).^2 + (a*(yg-y0)).^2;
  theta  = amp*(1 -tanh( (arg -c)/l^2 ));
  thetaT = fft2(theta); % spectral theta	 

	thetaT_xy = theta;

  Hnorm = 1;
  
 case{3}
  disp(' baroclinic instability ICs')

  ang = 0.83775;
  amp = -2.2;
  
  thetaB = fft2(amp*sin((2.*pi*xg/Lx) + ang));
  thetaT = fft2(amp*sin((2.*pi*xg/Lx)));

  thetaB_xy = amp*sin((2.*pi*xg/Lx) + ang);
  thetaT_xy = amp*sin((2.*pi*xg/Lx));

  Hnorm = 0;

 case{4}
  disp(' Gaussian upper vort ICs')

  amp = -5.;
  a= .5; b= a; % scale control
  x0 = 0; y0 = 0; %vortex origin
  
  arg   = (((xg-x0)/a).^2 + ((yg-y0)/b).^2)/2;
  
  theta  = amp*exp(-arg);
  thetaT = fft2(theta); % spectral theta
  thetaB = zeros(Ny,Nx);

  zetaBT = sum(abs(thetaT(:)).^2);
  
  zetaT  = -( (H/2*m.^2).*(1+cosh(H*m)).*IZo ).*thetaT;
  zetaT  = sum(abs(zetaT(:)).^2);
  
  Bind = zetaBT/zetaT

	thetaB_xy = thetaB;
	thetaT_xy = theta;

  Hnorm = 2;

 case{5}
  disp('2 vorticies put side by side')

  %amp = -1.; % cyclone (cool air mass in middle)
  amp = 1.; % anticyclone (warm air mass in middle)
  a= 2; b= a; % scale control
  x0 = 1; y0 = 0; % vortex to the right (play with x0)
  x1= -1; % vortex to the left (play with x0)

  arg = (b*(xg-x0)).^2 + (a*(yg-y0)).^2;
  arg1 = (b*(xg-x1)).^2 + (a*(yg-y0)).^2;

  theta1  = amp*exp(-arg1);
  theta = amp*exp(-arg) + theta1;

  thetaT = fft2(theta); % spectral theta

  thetaB = zeros(Ny,Nx);

	thetaB_xy = thetaB;
	thetaT_xy = theta;
  
	Hnorm = 0;

 case{6}
  disp(' edge wave w/ vort-norm')

  theta = cos(xg).*cos(yg);

  thetaT = fft2(theta);
  thetaB = zeros(Ny,Nx);

	thetaB_xy = thetaB;
	thetaT_xy = theta;
  
  Hnorm = 2;

 case{7}
  disp(' topo wave w/ vort-norm')

%  arg   = (xg.^2 + yg.^2)/2;
%  theta  = -exp(-arg);
  theta  = cos(xg).*sin(yg);

  thetaB = fft2(theta);
  
  coeff  = cosh(H*m) - (H*m.*sinh(H*m));
  
  thetaT = thetaB ./ coeff;

	thetaB_xy = theta;
	thetaT_xy = theta ./ coeff;
  
  Hnorm = 2;

 case{8}
  disp(' Cosine bell upper vort ICs')

  amp = 0.175;
  a= 1.75; b= a; % scale control
  x0 = 0; y0 = 0; %vortex origin
  
  arg   = sqrt(((xg-x0)/4*a).^2 + ((yg-y0)/4*b).^2);
 
 	theta  = ((amp/16).*(1+cos(pi.*arg)).^4).*(arg<=1);
  
	thetaT = fft2(theta); % spectral theta
  thetaB = zeros(Ny,Nx);

  zetaBT = sum(abs(thetaT(:)).^2);
  
  zetaT  = -( (H/2*m.^2).*(1+cosh(H*m)).*IZo ).*thetaT;
  zetaT  = sum(abs(zetaT(:)).^2);
  
  Bind = zetaBT/zetaT

	thetaB_xy = thetaB;
	thetaT_xy = theta;
  
  Hnorm = 2;

end % switch(Nic)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% normalize so that max leading-order wind is O(1).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lower boundary
Phi0B = -IZ.*thetaB + IZo.*thetaT; % spectral phi
Phi0Bx   = real(ifft2(    i*DX.*Phi0B));
Phi0By   = real(ifft2(    i*DY.*Phi0B));
magVB = ((Phi0Bx.^2) + (Phi0By.^2)).^.5;
% upper boundary
Phi0T = -IZo.*thetaB + IZ.*thetaT; % spectral phi
Phi0Tx   = real(ifft2(    i*DX.*Phi0T));
Phi0Ty   = real(ifft2(    i*DY.*Phi0T));
magVT = ((Phi0Tx.^2) + (Phi0Ty.^2)).^.5;

Vnorm = max([max(max(magVB)) max(max(magVT))]);

%  vorticity
Zeta0T = real(ifft2(-(m.^2).*(Phi0T)));
Zeta0B = real(ifft2(-(m.^2).*(Phi0B)));

Znorm  = max([max(Zeta0T(:)) max(Zeta0B(:))]);

switch(Hnorm)
 case{0}
  th_B = real(ifft2(thetaB));
  th_T = real(ifft2(thetaT));
  th_B = thetaB_xy;
  th_T = thetaT_xy;
 case{1}
  th_B = real(ifft2(thetaB/Vnorm));
  th_T = real(ifft2(thetaT/Vnorm));
  th_B = thetaB_xy/Vnorm;
  th_T = thetaT_xy/Vnorm;
 case{2}  % vorticity norm
  th_B = abs(amp)*real(ifft2(thetaB/Znorm));
  th_T = abs(amp)*real(ifft2(thetaT/Znorm));
  th_B = (abs(amp)/Znorm)*thetaB_xy;
  th_T = (abs(amp)/Znorm)*thetaT_xy;
end

%  plot here

%  colormap (cold -> hot)
setcolor; hj = hj32; PCA = PCA32;

figure(1);  clf
pcolor(x,y,th_B);shading flat
caxis([-1 1]*PCA * max(max(abs(th_B))));
title('Lower boundary potential temperature')
hold on;  axis image
colormap(hj); colorbar

figure(2);  clf
pcolor(x,y,th_T);shading flat
caxis([-1 1]*PCA * max(max(abs(th_T))));
title('Upper boundary potential temperature')
hold on;  axis image
colormap(hj); colorbar

figure(3);  clf
baro = (th_T-th_B)/H;
pcolor(x,y,baro);shading flat
caxis([-1 1]*PCA * max(max(abs(baro))));
title('\theta_T - \theta_B / H')
hold on;  axis image
colormap(hj); colorbar

%  th_init.nc
nc              = netcdf('th_init.nc','clobber');
nc('nx')        = Nx;
nc('ny')        = Ny;
nc{'thetaB'}    = ncfloat('ny','nx');
nc{'thetaT'}    = ncfloat('ny','nx');
nc{'thetaB'}(:) = th_B;
nc{'thetaT'}(:) = th_T;
close(nc);

% determine max flow induced on home & opposing boundaries
thetaB = fft2(th_B);  thetaT = fft2(th_T);

% lower boundary
Phi0B = -IZ.*thetaB;
Phi0Bx   = real(ifft2(    i*DX.*Phi0B));
Phi0By   = real(ifft2(    i*DY.*Phi0B));
magVBb = ((Phi0Bx.^2) + (Phi0By.^2)).^.5; % bottom on bottom

Phi0B = IZo.*thetaT; % spectral phi
Phi0Bx   = real(ifft2(    i*DX.*Phi0B));
Phi0By   = real(ifft2(    i*DY.*Phi0B));
magVBt = ((Phi0Bx.^2) + (Phi0By.^2)).^.5; % bottom on top

VBb = max(max(magVBb)); VBt = max(max(magVBt));
disp(['Max lower flow on lower, upper, and ratio: ' num2str(VBb,3) ...
		'  ' num2str(VBt,3) '  ' num2str(VBt/VBb,3)]);

% upper boundary
Phi0T = IZ.*thetaT; % spectral phi
Phi0Tx   = real(ifft2(    i*DX.*Phi0T));
Phi0Ty   = real(ifft2(    i*DY.*Phi0T));
magVTt = ((Phi0Tx.^2) + (Phi0Ty.^2)).^.5; % top on top

Phi0T = -IZo.*thetaB; % spectral phi
Phi0Tx   = real(ifft2(    i*DX.*Phi0T));
Phi0Ty   = real(ifft2(    i*DY.*Phi0T));
magVTb = ((Phi0Tx.^2) + (Phi0Ty.^2)).^.5; % top on top

VTt = max(max(magVTt)); VTb = max(max(magVTb));
disp(['Max upper flow on lower, upper, and ratio: ' num2str(VTt,3) ...
		'  ' num2str(VTb,3) '  ' num2str(VTb/VTt,3)]);




