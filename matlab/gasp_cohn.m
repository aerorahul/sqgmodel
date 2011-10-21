%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   gasp_cohn - function that calculates the Gaspari-Cohn localization 
%               function value between two different points.  Given 
%               by equation 4.10 of Gaspari and Cohn QJRMS 1999.,
%				volume 125
%
%   local = gasp_cohn(dist, r_max)
%
%     dist - distance between two grid points
%    r_max - radius where all covariances become zero
%    local - localization factor
%
%     created May 2004 Ryan Torn, U. Washington
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function local = gasp_cohn(dist, r_max);

b = r_max; a = 0.5 * b;
gcv = dist ./ a;

if  dist <= a % less than 1/2 length 

  local = -0.25*(gcv).^5 + 0.5*(gcv).^4 + ...
             (5/8)*(gcv).^3 - (5./3.)*(gcv).^2 + 1;

elseif ( dist > a ) & ( dist <= b ) % greater than one half length           

  local = (1./12.)*(gcv).^5 - 0.5*(gcv).^4 + (5./8.)*(gcv).^3 + ...
          (5./3.)*(gcv).^2 - 5*(gcv) + 4 - (2./3.)*(gcv).^-1;

else

  local = 0.0;

end;
