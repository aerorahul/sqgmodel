%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% <next few lines under version control, D O  N O T  E D I T>
% $Date$
% $Author$
% $Revision$
% $Id$
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
%
% Set up the colormap
%
% D. Muraki's:
% colormap (PCA adjusts caxis limits)
% This gives 128 colors:
z32 = zeros(1,32);  o64 = 64*ones(1,64);
hj128 = [[z32,(0:2:64),o64];[(0:1:64),(63:-1:0)];[o64,(64:-2:0),z32]]'/64;
PCA128 = 64/63;
% This gives 64 colors:
z16 = zeros(1,16);  o32 = 32*ones(1,32);
hj64 = [[z16,(0:2:32),o32];[(0:1:32),(31:-1:0)];[o32,(32:-2:0),z16]]'/32;
PCA64 = 32/31;
% This gives 32 colors:
z8 = zeros(1,8);  o16 = 16*ones(1,16);
hj32  = [[z8,(0:2:16),o16];[(0:1:16),(15:-1:0)];[o16,(16:-2:0),z8]]'/16;
PCA32 = 16/15;
% This gives 16 colors:
z4 = zeros(1,4); o8 = 8*ones(1,8);
hj16 = [[z4,(0:2:4),8,8,o8];[(0:1:6),8,8,8,(6:-1:0)];...
         [o8,8,8,(4:-2:0),z4]]'/8;
PCA16 = 8/7;
% This gives 8 colors:
z2 = zeros(1,2);  o4 = 4*ones(1,4);
hj4  = [[z2,(0:2:4),o4];[(0:1:4),(3:-1:0)];[o4,(8:-2:0),z2]]'/4;
PCA4 = 4/3;

