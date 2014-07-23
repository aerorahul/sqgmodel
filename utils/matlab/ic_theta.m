%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% <next few lines under version control, D O  N O T  E D I T>
% $Date$
% $Author$
% $Revision$
% $Id$
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% SPECIFY potential temperature (theta) at the tropopause.
% OBTAIN, None needed
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clear all; close all; clc;
warning off all;

% Localization function
% 1 = Gaussian
% 2 = Cosine Bell	
Nic = 2;

% Offset the initial disturbance from y = 0
offset = 0;

save_netcdf  = true;
save_figures = ~true;

% Domain size and Grid points
N  = 128;                      % Number of gridpoints
Nx = N; Ny = N/2;              % Grid points in X and Y
Lx =  20.0; Ly = 2*5.539118;   % Domain length
H = 1.0;                       % baroclinic instability
Ro = 0.0;                      % Rossby number
Hnorm = 1;                     % Norm to normalize (0:None, 1:Velocity, 2:Vorticity)

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
IZ  = real(1./(m.*tanh(m*H))); 
IZo = real(1./(m.*sinh(m*H))); 
m(1,1)  = 0;
IZ(1,1) = 1; IZo(1,1) = 1; 
DZ = 1./IZ;  DZo = 1./IZo; 
DZ(1,1) = 0; DZo(1,1) = 0; IZ(1,1) = 0; IZo(1,1) = 0;

switch(Nic)
	case{1}
  	disp(' Gaussian Formulation')

		amp = -5.0;
		a  = 0.5   ;  b = a;     % scale control
		x0 = 0     ;  y0 = y(Ny/2-offset);    % vortex origin

		arg = (((xg-x0)/a).^2 + ((yg-y0)/b).^2)/2;
		fxy = amp*exp(-arg);
 
		thetaT = fft2(fxy);           % spectral theta at Tropopause
		thetaB = zeros(Ny,Nx);        % spectral theta at Ground
		
	case{2}
		disp(' Cosine bell upper vort ICs')

		amp = 0.175;
		amp = -0.51375;
		amp = -1.0;
		a  = 1.75  ;  b  = a;    % scale control
		a  = 3.5   ;  b  = a;    % scale control
		x0 = 0     ;  y0 = y(Ny/2-offset);    % vortex origin
			
		arg = sqrt(((xg-x0)/4*a).^2 + ((yg-y0)/4*b).^2);
		fxy = ((amp/16).*(1+cos(pi.*arg)).^4).*(arg<=1);
		 
		thetaT = fft2(fxy);           % spectral theta at Tropopause
		thetaB = zeros(Ny,Nx);        % spectral theta at Ground
		
end % switch(Nic)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% normalize so that max leading-order wind is O(1).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lower boundary
Phi0B = -IZ.*thetaB + IZo.*thetaT; % spectral phi
Phi0Bx = real(ifft2(i*DX.*Phi0B));
Phi0By = real(ifft2(i*DY.*Phi0B));
vB = -Phi0Bx;
uB =  Phi0By;
magVB = ((vB.^2) + (uB.^2)) .^ 0.5;

% upper boundary
Phi0T = -IZo.*thetaB + IZ.*thetaT; % spectral phi
Phi0Tx = real(ifft2(i*DX.*Phi0T));
Phi0Ty = real(ifft2(i*DY.*Phi0T));
vT = -Phi0Tx;
uT =  Phi0Ty;
magVT = ((vT.^2) + (uT.^2)) .^ 0.5;

%	velocity norm
Vnorm = max([max(max(magVB)) max(max(magVT))]);

%  vorticity norm
Zeta0T = real(ifft2(-(m.^2).*(Phi0T)));
Zeta0B = real(ifft2(-(m.^2).*(Phi0B)));
Znorm  = max([max(Zeta0T(:)) max(Zeta0B(:))]);

switch(Hnorm)
	case{0}
		th_B = real(ifft2(thetaB));
		th_T = real(ifft2(thetaT));
 	case{1}	% velocity norm
		th_B = real(ifft2(thetaB/Vnorm));
		th_T = real(ifft2(thetaT/Vnorm));
 	case{2}	% vorticity norm
		th_B = abs(amp)*real(ifft2(thetaB/Znorm));
		th_T = abs(amp)*real(ifft2(thetaT/Znorm));
end

tol = 10^-6;
th_T = th_T.*(abs(th_T) > tol);

%	Plot
%  colormap (cold -> hot)
[hj, PCA] = setcolor(128);

f(1) = figure(1);  clf
pcolor(x,y,th_B);shading flat
caxis([-1 1]*PCA * max(max(abs(th_B))));
title('Lower boundary potential temperature')
hold on;  axis image
colormap(hj); colorbar
t{1} = ['Bottom.png'];

f(2) = figure(2);  clf
pcolor(x,y,th_T);shading flat
caxis([-1 1]*PCA * max(max(abs(th_T))));
title('Upper boundary potential temperature')
hold on;  axis image
colormap(hj); colorbar
t{2} = ['Toppom.png'];

f(3) = figure(3);  clf
baro = (th_T-th_B)/H;
pcolor(x,y,baro);shading flat
caxis([-1 1]*PCA * max(max(abs(baro))));
title('\theta_T - \theta_B / H')
hold on;  axis image
colormap(hj); colorbar
t{3} = ['Toppom-Bottom.png'];

f(4) = figure(4); clf
surf(x,y,th_T); shading flat;
caxis([-1 1]*PCA * max(max(abs(th_T))));
hold on; axis image 
title('Surface plot')
colormap(hj); colorbar
t{4} = ['Surface.png'];

f(5) = figure(5); clf
plot(x,th_T(Ny/2,:),'-r.');
title('X - Cross Section')
t{5} = ['x-sec.png'];

f(6) = figure(6); clf
plot(y,th_T(:,Nx/2),'-r.');
title('Y - Cross Section')
t{6} = ['y-sec.png'];

%	Dump to disk:
if( save_netcdf )
	nc = netcdf('th_init.nc','clobber');
	nc('nx')=Nx;
	nc('ny')=Ny;
	nc{'thetaB'}=ncfloat('ny','nx');
	nc{'thetaT'}=ncfloat('ny','nx');
	nc{'thetaB'}(:)=th_B;
	nc{'thetaT'}(:)=th_T;
	nc.('Localization') = 'theta';
	close(nc);
end
	
if( save_figures )
	for i = 1:6
		savefigure(f(i),t{i});
	end
end

% determine max flow induced on home & opposing boundaries
thetaB = fft2(th_B);  thetaT = fft2(th_T);

% lower boundary
Phi0B = -IZ.*thetaB; % spectral phi
Phi0Bx   = real(ifft2(i*DX.*Phi0B));
Phi0By   = real(ifft2(i*DY.*Phi0B));
magVBb = ((Phi0Bx.^2) + (Phi0By.^2)).^.5; % bottom on bottom

Phi0B = IZo.*thetaT; % spectral phi
Phi0Bx   = real(ifft2(i*DX.*Phi0B));
Phi0By   = real(ifft2(i*DY.*Phi0B));
magVBt = ((Phi0Bx.^2) + (Phi0By.^2)).^.5; % bottom on top

VBb = max(max(magVBb)); VBt = max(max(magVBt));
disp(['Max lower flow on lower, upper, and ratio: ' num2str(VBb,3) ...
		'  ' num2str(VBt,3) '  ' num2str(VBt/VBb,3)]);

% upper boundary
Phi0T = IZ.*thetaT; % spectral phi
Phi0Tx   = real(ifft2(i*DX.*Phi0T));
Phi0Ty   = real(ifft2(i*DY.*Phi0T));
magVTt = ((Phi0Tx.^2) + (Phi0Ty.^2)).^.5; % top on top

Phi0T = -IZo.*thetaB; % spectral phi
Phi0Tx   = real(ifft2(i*DX.*Phi0T));
Phi0Ty   = real(ifft2(i*DY.*Phi0T));
magVTb = ((Phi0Tx.^2) + (Phi0Ty.^2)).^.5; % top on top

VTt = max(max(magVTt)); VTb = max(max(magVTb));
disp(['Max upper flow on lower, upper, and ratio: ' num2str(VTt,3) ...
		'  ' num2str(VTb,3) '  ' num2str(VTb/VTt,3)]);
