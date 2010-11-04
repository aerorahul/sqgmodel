%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% <next few lines under version control, D O  N O T  E D I T>
% $Date$
% $Author$
% $Revision$
% $Id$
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%  normalize_ic - normalize and scale ic
%
%  fname - name of the file to read/write from/to
%  scale - scale to factor after normalizing by max. value
%
%    normalize_ic(fname,scale)
%
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function normalize_ic(fname,scale)

nc = netcdf(fname,'write');

th = nc{'thetaT'}(:);
maxval = max(max(abs(th)));
disp(['max(max(abs(thetaT))) = ' num2str(maxval)])
if (maxval ~= 0)
	th = th ./ (maxval * scale);
end
nc{'thetaT'}(:) = th;

th = nc{'thetaB'}(:);
maxval = max(max(abs(th)));
disp(['max(max(abs(thetaB))) = ' num2str(maxval)])
if (maxval ~= 0)
	th = th ./ (maxval * scale);
end
nc{'thetaB'}(:) = th;

close(nc);

return;
