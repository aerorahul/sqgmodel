%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% <next few lines under version control, D O  N O T  E D I T>
% $Date$
% $Author$
% $Revision$
% $Id$
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%  ic_tracer - create tracer field of unit magnitude 
%              will write out tr_init.nc 
%
%  fname - name of the th_init.nc file to read from
%
%  ic_tracer('th_init.nc')
%
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function ic_tracer(fname)

nc = netcdf(fname,'nowrite');
thB = nc{'thetaB'}(:);
thT = nc{'thetaT'}(:);
close(nc);

Nx = size(thB,2);
Ny = size(thT,1);

maxval = max(max(abs(thB)));
trB = abs(thB ./ maxval);
maxval = max(max(abs(thT)));
trT = abs(thT ./ maxval);

nc = netcdf('tr_init.nc','clobber');
nc('nx') = Nx;
nc('ny') = Ny;
nc{'tracerB'} = ncfloat('ny','nx');
nc{'tracerT'} = ncfloat('ny','nx');
nc{'tracerB'}(:) = trB;
nc{'tracerT'}(:) = trT;
close(nc);

return;
