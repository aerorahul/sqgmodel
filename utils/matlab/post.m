%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% <next few lines under version control, D O  N O T  E D I T>
% $Date$
% $Author$
% $Revision$
% $Id$
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
%
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
%  sqg post-proc;  01 sept 99 (happy bastille day)
%						 08 july 00 (greg in BoCo)
%                  15 july 05 (netcdf option works; R. Mahajan & GJH)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
%
%  le menu

% relative path to smat files
%pth = '../run';
pth = './';
wd = cd;

% set file type here by uncommenting the appropriate line
ftype = 'netcdf'; fsuf = 'nc';
%ftype = ' ascii'; fsuf = 'dat';

disp(['-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='])
disp(['currently set to read ' ftype ' files...']);
disp(['-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='])

clf;  go = 0;
disp(' ')
disp('  sqg post-processor ')
disp('    R)  exit')
disp('    0)  read smat.dat & initialize (-100 = BOTTOM; 100 = TOP)')
disp('    1)  read file & initialize (11 = _specific_ filename')
disp('    2)  toggle jet on')
disp('    3)  toggle jet off (default)')
disp('    4)  contour slice')
disp('    5)  pcolor slice')
disp('    6)  contour movie')
disp('    7)  pcolor movie')
disp('    8)  read th_init.dat & initialize')
disp('    9)  final -> th_init')
disp('   10)  mean theta')
disp('   11)  turbulent vortices')
disp('   12)  no-mean theta')
disp('   13)  no-mean theta movie')

while(go >=0)

go = input('  choice (return to exit):  ');
if isempty(go);  go = -1;  end

switch go
%  read input file & initialize
	case {0,100,-100,1,8,11}
		if (go==0)
			head = 1; 
			fname = strcat('smat.',fsuf);
		elseif (go==-100)
			head = 1; 
			fname = strcat([pth '/smat.'],fsuf);
			varname = 'thetaB';
			disp('Reading smat.nc ... ')
		elseif (go==100)
			head = 1; 
			fname = strcat([pth '/smat.'],fsuf);
			varname = 'thetaT';
			disp('Reading smat.nc ... ')
		elseif (go==1)
%		  flen = 24; % length of filenames
		  flen = 10; % length of filenames
		  wd = cd;  cd (pth);
		  fchk = ls; % load directory info
		  fi = findstr('smat_',fchk);
		  disp(' ');disp([pth ' contains the following smat output files...'])
		  disp(' ')
		  for j = 1:1:length(fi)
			 fnam(j) = cellstr(fchk(fi(j):fi(j)+flen));
			 disp([' ' num2str(j) ' = ' char(fnam(j))]);
		  end
		  disp(' '); go = input(' Enter a file number:  ');
		  if isempty(go);  go = -1;  end
		  disp([' You selected: ' char(fnam(go))]);
		  fname = char(fnam(go));
%		  fid  = fopen(char(fname(go)));
		  head = 1; 
		elseif (go==11)
		  head = 1; 
		  disp(' '); fname = input(' Enter a file name:  ','s');
%		  fid     = fopen(fname);		  
      else
			head = 0; 
			fid     = fopen('th_init.dat');
%			fid     = fopen('/mmmtmp/muraki/sqg/th_init.dat');
		end

		if ftype == 'netcdf'
		  disp(' ')
		  disp(' Reading data from netCDF...')
		  nc = netcdf(fname);
		  timevar=nc{'time'}(:); 
		  avar=nc{char(varname)}(:);
		  Nt=size(avar,1);
		  Ny=size(avar,2);
		  Nx=size(avar,3);
		  Lx=nc.('XL')(1);
		  Ly=nc.('YL')(1);
		  H =nc.('H')(1);
	
		  disp(' ')
		  disp(['  [Nx Ny Nt] = ' num2str(Nx) ' , ' ...
				  num2str(Ny) ' , ' num2str(Nt)])
		  disp(['  [Lx Ly H ] = ' num2str(Lx) ' , ' ...
				  num2str(Ly) ' , ' num2str(H)])
		  
		  avar2=reshape(avar,[Nt,Ny,Nx]);
		  
		  th_plot = permute(avar2,[2,3,1]);

		else
		  fid = fopen(fname);		  
		  th_read = fscanf(fid,'%g',inf);
		  fclose(fid);

		  if(head==1)
			 Nx = th_read(1);  Ny  = th_read(2); H = th_read(3);
			 Lx = th_read(4);  Ly  = th_read(5);  Ro = th_read(6);
			 nD = th_read(7);  tau = th_read(8);
			 ntims = th_read(9);  iplot = th_read(10);  dt = th_read(11);
			 
			 Nt = floor((size(th_read,1)-11)/(Nx*Ny));
			 th_plot = reshape(th_read(12:Ny*Nx*Nt+11),Ny,Nx,Nt);
			 
			 disp(' ')
			 disp(['  [Nx Ny Nt] = ' num2str(Nx) ' , ' ...
					 num2str(Ny) ' , ' num2str(Nt)])
			 disp(['  [Lx Ly H ] = ' num2str(Lx) ' , ' ...
					 num2str(Ly) ' , ' num2str(H)])
			 disp(['  [Ro tau  ] = ' num2str(Ro) ' , ' num2str(tau)])
		  else
			 Nx = th_read(1);  Ny  = th_read(2);
			 Nt = floor((size(th_read,1)-2)/(Nx*Ny));
			 th_plot = reshape(th_read(3:Ny*Nx*Nt+2),Ny,Nx,Nt);
			 
			 disp(' ')
			 disp('  initial condition analysis')
			 disp(' ')
			 disp(['  [Lx Ly Ro] = ' num2str(Lx) ' , ' ...
					 num2str(Ly) ' , ' num2str(Ro)])
		  end
		end

		%  set coordinate vectors
		xx = 0:Lx/Nx:Lx;  x = xx(1:Nx) - (Lx/2);
		yy = 0:Ly/Ny:Ly;  y = yy(1:Ny) - (Ly/2);
		
		[xg, yg] = meshgrid(x,y);
		
		Jet  = 0;  
		cxx  = max(max(abs(th_plot(:,:,1)- Jet*yg)));
		tcon = ([-1:0.1:1] *cxx);
		
		% colormap (PCA adjusts caxis limits)
    cxx = max(max(abs(th_plot(:,:,1)- Jet*yg)));
    z16 = zeros(1,16);  o32 = 32*ones(1,32);
    hj  = [[z16,(0:2:32),o32]; ...
					[(0:1:32),(31:-1:0)];[o32,(32:-2:0),z16]]'/32;
    hj  = [0 0 0 ; hj ; 0 0 0 ]; % blacken out of range
    PCA = 32/31;

%  toggle Jet ON
  case {2}
    Jet = 1;    
    cxx  = max(max(abs(th_plot(:,:,1)- Jet*yg)));
    tcon = ([-1:0.1:1] *cxx);

%  toggle Jet OFF
  case {3}
    Jet  = 0;  
    cxx  = max(max(abs(th_plot(:,:,1)- Jet*yg)));
    tcon = ([-1:0.1:1] *cxx);

%  contour slice
  case {4}
    slice = Nt;  disp(' ')
    while (slice > 0)
		clf;
		%		Jet = 0.5;
		Jet = 0.0;
		cxx  = max(max(abs(th_plot(:,:,slice)- Jet*yg))); % BI
%		tcon = ([-1.5:0.1:1.5]);
		tcon = ([-1:0.1:-.1 .1:.1:1] *cxx);
%		cint=0.553; tcon = -[0:1.:5]*cint -.3; % RR
      mask = th_plot(:,:,slice).*(th_plot(:,:,slice) > 0);
%      [c h] = contour(x,y,mask,tcon,'k-'); set(h,'Linewidth',2); hold on;
      [c h] = contour(x,y,th_plot(:,:,slice),tcon); 
%		set(h,'Linewidth',2); hold on;
      mask = th_plot(:,:,slice).*(th_plot(:,:,slice) < 0);
%      [c h] = contour(x,y,mask,tcon,'k-'); set(h,'Linewidth',2);
%      [c h] = contour(x,y,mask,tcon,'k--'); set(h,'Linewidth',2);
		axis equal; % axis square;  
		axis([-Lx Lx -Ly Ly]/2); % axis square;  
		h2 = gca; 
%		set(h2,'Xlim',[-2 2]); set(h2,'Linewidth',2); set(h2,'Fontweight','bold');
% special set for 2X domain:
%set(h2,'Ylim',[-5.5391/2 5.5391/2]);
      caxis([-1 1]*cxx);  colorbar
      cc = text(x(1)+Lx/20,y(1)+Ly/20,[num2str(slice)]);
      set(cc,'fontsize',18)
      slice = input('  slice #,R:  ');
		hold off;
      if isempty(slice);  slice = -1;  end
    end

%  pcolor slice
  case {5}
    slice = 1;  disp(' ')
    while (slice > 0)
		if (Jet == 1)
		  [ubase,tbase] = HoskinsWestJet(1,Nx,Ny,Lx,Ly);
		  tplot = th_plot(:,:,slice) - tbase;
		else
		  tplot = th_plot(:,:,slice);		  
		end
		clf;
		cxx  = max(max(abs(tplot))); 
      pcolor(x,y,tplot);
      axis image;
      colormap(hj);  
      colormap(jet);  
%		shading flat;
		shading interp;
%      caxis([-1 1]*PCA * cxx);  colorbar % perts only
      caxis([-1 0]*PCA * cxx);  colorbar % + basic state
      cc = text(x(1)+Lx/20,y(1)+Ly/20,[num2str(slice)]);
      set(cc,'fontsize',18)
      STRslice = input('  slice #,f,b,R:  ','s');
      if isempty(STRslice)  
	  slice = -1;
	else
	  switch STRslice
	    case {'f'}
	      slice = min([Nt slice+1]);
	    case {'b'}
	      slice = max([ 1 slice-1]);
	    otherwise
	      slice = min([Nt max([ 1 str2num(STRslice) ])]);
          end
      end
    end

%  contour movies
  case {6}
    clear mov; mov = moviein(Nt);

    for ii = 1:Nt
      contour(x,y,th_plot(:,:,ii) - Jet*yg,tcon)
		axis equal; % axis square;  
		axis([-Lx Lx -Ly Ly]/2); % axis square;  
      cc = text(x(1)+Lx/20,y(1)+Ly/20,[num2str(ii)]);
      set(cc,'fontsize',18)
      mov(:,ii) = getframe;
    end
    disp(' ')
    disp('  to play: movie(mov,5,1)')

%  pcolor movies
  case {7}
    clear mov; mov = moviein(Nt);

	 ii = 1;
    for jj = 1:1:Nt
		cxx = max(max(max(th_plot(:,:,:))));
		ii = ii + 1;
      pcolor(x,y,th_plot(:,:,jj)- Jet*yg)
%      disp([ii max(max(th_plot(:,:,ii))) min(min(th_plot(:,:,ii)))])
      caxis([-1 1]*PCA * cxx)
%      caxis([-1 0]*PCA * cxx);  colorbar % + basic state
%      colormap(hj);  shading interp
      colormap(jet);  shading interp
%      cc = text(x(1)+Lx/20,y(1)+Ly/20,[num2str(ii)]);
%      set(cc,'fontsize',18)
      axis image;  colorbar
      mov(:,ii) = getframe;
		if 1 == 0
		filn = 'dd_qgp';
		if (ii < 10)
		  outf = strcat(filn,'0',num2str(ii),'.jpg');
		  print('-djpeg40','-r100','-zbuffer',outf);
		else
		  outf = strcat(filn,num2str(ii),'.jpg');
		  print('-djpeg40','-r100','-zbuffer',outf);
		end
		hold off;
    end % if 1 == 1
	 end 

    disp(' ')
    disp('  to play: movie(mov,5,1)')
    disp('  to show: pcolor(x,y,th_plot(:,:,1));shading interp')
    disp('           caxis([-1 1]/2*cxx);axis image;colorbar')

%  final -> th_init.dat
  case {9}
    fid = fopen('th_init.dat','w');
    fprintf(fid,'%i\n',Nx);
    fprintf(fid,'%i\n',Ny);
    fprintf(fid,'%f\n',th_plot(:,:,end));
    fclose(fid);

%  mean theta
  case {10}
    mtheta = [];
    for jj=1:Nt
      mtheta = [mtheta ; sum(sum(th_plot(:,:,jj)))/(Nx*Ny)];
    end

    plot(0:Nt-1,mtheta,'r.')

%  turbulent vortices
  case {11}
    disp(' ')
    slice = input('  slice #:  ');
    pcolor(x,y,0.5*(1+sign(abs(th_plot(:,:,slice))-2)).*th_plot(:,:,slice))
    caxis([-1 1]*PCA * cxx)
    colormap(hj); shading flat
    cc = text(x(1)+Lx/20,y(1)+Ly/20,[num2str(slice)]);
    set(cc,'fontsize',18)
    axis image;  colorbar

%  pcolor slice
  case {12}
    slice = Nt;  disp(' ')
    while (slice > 0)
      nomean = th_plot(:,:,slice) - sum(sum(th_plot(:,:,slice))/(Nx*Ny));
      nxx    = max(max(abs(nomean)));
      pcolor(x,y,nomean);
      axis image;
      colormap(hj);  shading flat
      caxis([-1 1]*PCA * nxx/1);  colorbar
      cc = text(x(1)+Lx/20,y(1)+Ly/20,[num2str(slice)]);
      set(cc,'fontsize',18)
      STRslice = input('  slice #,f,b,R:  ','s');
      if isempty(STRslice)  
	  slice = -1;
	else
	  switch STRslice
	    case {'f'}
	      slice = min([Nt slice+1]);
	    case {'b'}
	      slice = max([ 1 slice-1]);
	    otherwise
	      slice = min([Nt max([ 1 str2num(STRslice) ])]);
          end
      end
    end

%  pcolor movies
  case {13}
    clear mov; mov = moviein(Nt);

    for ii = 1:Nt
      meanth = sum(sum(th_plot(:,:,ii)))/(Nx*Ny);
      pcolor(x,y,th_plot(:,:,ii)- Jet*yg - meanth)
%      disp([ii max(max(th_plot(:,:,ii))) min(min(th_plot(:,:,ii)))])
      caxis([-1 1]*PCA*(cxx-meanth))
      colormap(hj);  shading interp
      cc = text(x(1)+Lx/20,y(1)+Ly/20,[num2str(ii)]);
      set(cc,'fontsize',18)
      axis image;  colorbar
      mov(:,ii) = getframe;
    end
    colorbar
    disp(' ')
    disp('  to play: movie(mov,5,1)')
    disp('  to show: pcolor(x,y,th_plot(:,:,1));shading interp')
    disp('           caxis([-1 1]/2*cxx);axis image;colorbar')

  otherwise
    go = -1;
end %switch

if (go ~= -1) 
  go = 0;
end %if

disp(' ')
end %while

% return to original directory
cd (wd);
