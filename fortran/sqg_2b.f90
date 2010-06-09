!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! <next few lines under version control, D O  N O T  E D I T>
! $Date$
! $Author$
! $Revision$
! $Id$
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This block under version control: DO NO EDIT
!
! Date: 2009-05-01 11:00:35 -0700 (Fri, 01 May 2009) 
! Revision: 45 
! Author: hakim 
! Id: sqg_2b.f90 45 2009-05-01 18:00:35Z hakim 
!
! new svn version started 12 september 2005 from...
! sqg-2b.f version 2.7 Tue Sep  6 17:29:41 PDT 2005
!
! Originators: G. J. Hakim,  University of Washington
!              D. J. Muraki, Simon Fraser University
!
! Fortran90 code for a uniform-PV two-surface QG+1 model in spectral form.
!
! Compile: ifort sqg_2b.f90 fft_90.f90 (make 2sqg)
!
! sQG dynamics are recovered for: Ross=0.0, H-->\infty
!  2D dynamics are recovered for: Ross=0.0, H-->0; \theta->\theta*H

PROGRAM sqg_spectral

  call main
  stop

END PROGRAM sqg_spectral

SUBROUTINE main

!PROGRAM sqg_spectral
  ! Originator: G. Hakim, University of Washington
  USE spectral

  IMPLICIT NONE

  ! spectral values
  complex, dimension(2*kmax,2*lmax) :: tspB,toldB,sspB,tbB,tbT
  complex, dimension(2*kmax,2*lmax) :: tspT,toldT,sspT,sbold,sB
  ! grid point values
  real, dimension(mmax,nmax) :: txB,tyB,vB,uB
  real, dimension(mmax,nmax) :: txT,tyT,vT,uT
  ! spectral work arrays
  real, dimension(2*kmax,2*lmax) :: oxyB,oxyT
  ! tendencies
  real, dimension(mmax,nmax) :: ttB,ttT,lap,blank,hx,hy,hu,hv 
  complex, dimension(mmax,nmax) :: ttspB,ttspT
  ! basic state
  real, dimension(mmax,nmax) :: ulinB,ulinT,tbyB,tbyT
  ! misc
  real dx,dx2,dy,dy2,x,asx,asy,asig,rr,ctemp,rtemp,ak,bl,dco,btc
  real hamp,amp,ran1,noise,time,cxB,cyB,cxT,cyT,lam
  integer i,j,k,l,irec,irecl,LL,MM,idu,itime,kk,iseed
  logical first,bot,top

  ! Pertaining to netCDF files
  logical head
  integer cdfid(2)
  character(len=64) cdf_file(2) 

!!!!!!!!!!!!!!!!!!!!!!!!
! Start of executable: !
!!!!!!!!!!!!!!!!!!!!!!!!

  blank = 0.

  ! initialize diffusion (n,kmax,lmax,tau,dco):
  call diff(dco); ! dco = 0.
  ! dco = 2.2204460E-16 ! fixed across all resolution

  if (verbose .gt. 0) print *,'Running with Rossby number: ',Ross

! netCDF file output
  if(cdf) then
      if (verbose .gt. 0) print *, '...writing to netCDF files'
      
      head=.TRUE.
      cdf_file(1)= path//'smat_B.nc'; cdf_file(2)= path//'smat_T.nc'
      
      call cdf_dump(cdf_file(1),head,oxyB,0)
      call cdf_dump(cdf_file(2),head,oxyT,0)
  else
     if (verbose .gt. 0) print *, '...writing to ASCII files'

      ! matlab output file ([Nx XL YL Ross n tau ; Ny ntims iplot dt]):
      open(12,file=path//'smat_B.dat') 
      write(12,'(2i5,f16.8)') 2*kmax, 2*lmax, H
      write(12,'(3f16.8)') XL, YL, Ross
      write(12,'(i5,f16.8)') n, tau
      write(12,'(i6,i8,f16.8)') ntims, iplot, dt  
      ! matlab output file ([Nx XL YL Ross n tau ; Ny ntims iplot dt]):
      open(13,file=path//'smat_T.dat') 
      write(13,'(2i5,f16.8)') 2*kmax, 2*lmax, H
      write(13,'(3f16.8)') XL, YL, Ross
      write(13,'(i5,f16.8)') n, tau
      write(13,'(i6,i8,f16.8)') ntims, iplot, dt  
  endif

  ! initialize (read from Matlab file):
  call init(oxyB,oxyT)
  !oxyB = 0.; 
  !oxyT = 0. ! override

  ! advection flags
  top = .TRUE.; bot = .TRUE.

 ! if (maxval(abs(oxyT)) .lt. 1.e-5) top = .FALSE.
 ! if (maxval(abs(oxyB)) .lt. 1.e-5) bot = .FALSE.
  if (model .eq. 2) then ! tropo sQG
     top = .TRUE.; bot = .FALSE.
  elseif (model .eq. 3 .or. model .eq. 0) then ! surface sQG
     top = .FALSE.; bot = .TRUE.
  elseif (model .eq. 4) then ! tropo HsQG
     top = .TRUE.; bot = .FALSE.
  endif

  if (verbose .gt. 0) print*,'max BOTTOM initial value=',maxval(abs(oxyB)),bot
  if (verbose .gt. 0) print*,'max TOP initial value=',maxval(abs(oxyT)),top
  ! map into spectral space at the same resolution:
  call xy_to_sp(cmplx(oxyB,0.),tspB,2*kmax,2*lmax,kmax,lmax)
  call xy_to_sp(cmplx(oxyT,0.),tspT,2*kmax,2*lmax,kmax,lmax)

  ! zero k = 1
  !      tspB(2,:) = 0.; tspB(2*kmax,:) = 0.
  !      tspT(2,:) = 0.; tspT(2*kmax,:) = 0.

  ! init jet
  if (hw) then 
     call init_jet(tbB,tbT,tbyB,tbyT,ulinB,ulinT,lam)
     ulinB = ulinB + 0; ulinT = ulinT - 0*H*lam! barotropic wind (Ross = 0!)
  else
     tbB = 0;tbT = 0;tbyB=0.;tbyT=0.;ulinB=0.;ulinT=0.;lam=0.
  endif
  irec = 0; first = .TRUE.
  if (verbose .gt. 1) print*,'lam = ',lam
  if (verbose .gt. 1) print*,'extrema ulinT = ',maxval(ulinT),minval(ulinT)

  open(21,file=path//'mean.dat',status='replace')
  write(21,*) 'Mean theta, time:'; close(21)

! initialize terrain:
  hamp = 0.6
  hu = 0.; hv = 0.; hx = 0.; hy = 0.
  if (verbose .gt. 1)   print*,'iterr = ',iterr
  if (iterr .and. Ross .eq. 0) then 
     ! initialize terrain
     call terrain(hamp,hx,hy,hu,hv)
     call invert(tspB,0*tspT,txB,txT,tyB,tyT,vB,vT,uB,uT,tbB,tbT, &
          &        tbyB,tbyT,ulinB,ulinT,first, &
          &        .TRUE.,.TRUE.,0*lam,sB,sbold,lap)
     hu = uT; hv = vT
     if (verbose .gt. 1) print*,'extrema hx = ',maxval(hx),minval(hx)
     if (verbose .gt. 1) print*,'extrema hy = ',maxval(hy),minval(hy)
     if (verbose .gt. 1) print*,'extrema hu = ',maxval(hu),minval(hu)
     if (verbose .gt. 1) print*,'extrema hv = ',maxval(hv),minval(hv)
  endif

!!!!!!!!!!!!!!!!!!!!
! BEGIN: Time loop !
!!!!!!!!!!!!!!!!!!!!
  do itime = 1, ntims

! option to zero l = 0 modes, which are nonlinear solutions
!     tspB(:,1) = 0; tspT(:,1) = 0;

     time= real(itime-1)*dt
     if (mod((itime-1),imean) .eq. 0) then
        open(21,file=path//'mean.dat',position='append')
        if (verbose .gt. 0) write(6,*) '...writing to mean.dat'
        write(21,*)  real(tspB(1,1)),real(tspT(1,1))
        write(21,*) time
        close(21)
     endif

     ! this is for growing the most unstable mode
     if (grow .and. cyT .gt. .001) then ! rescale 
        tspB = tspB/2.; tspT = tspT/2.
        cyT = 0.
        if (verbose .gt. 1) print*,'RESCALING'
     endif

     ! save old streamfunction for Ekman calculation.
     sbold = sB; if (first) sbold = 0.

     ! Invert theta for streamfunction; compute gradients for advection:
     call invert(tspB,tspT,txB,txT,tyB,tyT,vB,vT,uB,uT,tbB,tbT,tbyB,tbyT,ulinB,ulinT,first,bot,top,lam,sB,sbold,lap)
     ! option to compute potential enstrophy norm and growth rate:
     ! call norm(tspB,tspT,itime)

     ! write data to da file for plots (last entry for linear field +):
     if (mod((itime-1),iplot) .eq. 0) then
        if (hw) then ! add basic state to output
           call dump(tspB+tbB,tspT+tbT,1,lam,itime,cdf_file,cdfid)
        else ! no basic state
           call dump(tspB,tspT,0,lam,itime,cdf_file,cdfid)
        endif
     endif

     ! spectral advection:
     if (bot) call advect(uB,vB,txB,tyB,tbyB,hx,hy,ulinB,ttB,lam,lap)
     if (top) call advect(uT+hu,vT+hv,txT,tyT,tbyT,blank,blank,ulinT,ttT,lam,blank)

     ! Write Courant numbers to stdout
     if (mod((itime-1),10) .eq. 0) then 
        cxB = maxval(abs(uB+ulinB))*dt/(XL/real(2*kmax))
        cyB = maxval(abs(vB))*dt/(YL/real(2*lmax))
        cxT = maxval(abs(uT+ulinT))*dt/(XL/real(2*kmax))
        cyT = maxval(abs(vT))*dt/(YL/real(2*lmax))
        write(*,'(A13,F10.3,4F8.3)') 'time,cx,cy = ', real(itime-1)*dt, &
             &           cxB,cyB,cxT,cyT
     endif

     ! FFT back to spectral space:
     if (bot) then 
        ttspB = cmplx(ttB,0.); call ft_2d(ttspB,mmax,nmax,-1)
     endif
     if (top) then 
        ttspT = cmplx(ttT,0.); call ft_2d(ttspT,mmax,nmax,-1)
     endif

     ! advance one time step with explicit (hyper-)diffusion:
     if (bot) call tadvB(tspB,ttspB,dco,first)
     if (top) call tadvT(tspT,ttspT,dco,first) 
     first = .FALSE.
     tspB(kmax,:) = 0.; tspB(lmax,:) = 0.
  end do ! itime

!!!!!!!!!!!!!!!!!!
! END: Time loop !
!!!!!!!!!!!!!!!!!!

  close(12); close(13)

  ! write final time to disk for future restart:
  irecl = 2*kmax*2*lmax*32
  open(8,file=path//'final+1.dat',access='direct', &
       &     form='unformatted', recl = irecl,status='unknown')
  call sp_to_xy(tspB,oxyB,kmax,lmax,2*kmax,2*lmax)
  call sp_to_xy(tspT,oxyT,kmax,lmax,2*kmax,2*lmax)
  write(8,rec=1) oxyB
  write(8,rec=2) oxyT

  close(8)

  !      open(9,file=path//'linear+1.dat',status='unknown')
  !      write(9,*) btc
  !      do j=1,nmax; write(9,*) tbby(j),ulin(j); enddo
  !      do j=1,2*lmax; write(9,*) tbb(j); enddo
  !      close(9)

end SUBROUTINE main

!!!!!!!!!!!!!!!
! SUBPROGRAMS !
!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE invert(itspB,itspT,txBr,txTr,tyBr,tyTr,vBr,vTr,uBr,uTr, &
     &    tbB,tbT,tbyB,tbyT,ulinB,ulinT,first,bot,top,lam,sB,sold,sblre)

  USE spectral
  IMPLICIT NONE

  ! Invert PV and transform to a larger grid for de-aliasing.

  complex, intent(in), dimension(2*kmax,2*lmax) :: itspB,itspT, &
       &                                               tbB,tbT,sold
  complex, intent(out), dimension(2*kmax,2*lmax) :: sB
  real, intent(in) :: lam
  real, intent(in), dimension(mmax,nmax) :: tbyB,tbyT,ulinB,ulinT
  logical, intent(in) :: first,bot,top
  complex, dimension(mmax,nmax) :: txB,txT,tyB,tyT, &
       &                                 uB,uT,vB,vT,copy
  real, intent(out), dimension(mmax,nmax) :: txBr,txTr,tyBr,tyTr, &
       &                                           uBr,uTr,vBr,vTr,sblre

  complex, dimension(2*kmax,2*lmax) :: tspB,tspT,temps
  complex, dimension(mmax,nmax) :: szspB,szspT,szzspB,szzspT, &
       &                                 u1B,v1B,u1T,v1T,temp

  complex, dimension(2*kmax,2*lmax) :: u1spB,u1spT,v1spB,v1spT, &
       &                                     tempspB,tempspT,htempB, &
       &                                     htempT,tempxyB,tempxyT

  complex,dimension(2*kmax,2*lmax),save:: dx,dy,dz,dzo,iz,izo,dzi,Id
  real, dimension(mmax,nmax) :: tx,ty,u,v
  integer, save :: pf,gf,bf
  integer :: i,j,k,l
  logical, save :: correct

!!!!!!!!!!!!!!!!! Set-up on first time step !!!!!!!!!!!!!!!!!!!!!!!!!!
  if (first) then 
     if (verbose .gt. 1) print*,'...calling setup'
     call d_setup(dx,dy,dz,dzo,iz,izo,Id)    ! derivative operators
     !         print*,maxval(abs(dx))
     !         print*,maxval(abs(dy))
     !         print*,maxval(abs(dz))
     !         print*,maxval(abs(iz))
     !         print*,maxval(abs(dzo))
     !         print*,maxval(abs(izo))
     correct = .FALSE. ! flag for next-order corrections
     if (Ross .gt. 1.e-5) correct = .TRUE. 
     ! switches for removing rotation or divergence (sQG only)
     !         pf = 1; gf = 1; bf = 0  ! control: full sqg+1 corrections
     !         pf = 0; gf = 1; bf = 1  ! no rotational correction
     !         pf = 1; gf = 0; bf = -1 ! no divergent correction
  endif

  tspB = itspB; tspT = itspT ! local copies of spectral boundary theta
!!!!!!!!!!!!!!!!!!!!!!!!!!! Grad Theta !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  txT = 0.; tyT = 0.; vT = 0.; uT = 0.
  txB = 0.; tyB = 0.; vB = 0.; uB = 0.

  call d_s2b(tspB,txB,1,dx) ! x derivative: small to big domain.
  call d_s2b(tspB,tyB,1,dy) ! y derivative: small to big domain.
  call d_s2b(tspT,txT,1,dx) ! x derivative: small to big domain.
  call d_s2b(tspT,tyT,1,dy) ! y derivative: small to big domain.

!!!!!!!!!!!!!!!!!!!!!!!!!!! sQG Bottom !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (bot .or. correct) then ! velocity contribution from lower boundary
     call d_b2b(txB,vB,1,-iz) 
     call d_b2b(-tyB,uB,1,-iz)
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!! sQG Toppom !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (top .or. correct) then ! velocity contribution from upper boundary
     call d_b2b(txT,vT,1,iz) ! z integral: big to big domain.
     call d_b2b(-tyT,uT,1,iz) ! z integral: big to big domain.
  endif
!!!!!!!!!!!!!!!!!!!!!!!! sQG Cross Boundary !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (bot .and. top .or. correct) then! lower velocity from upper boundary
     call d_b2b(txT,temp,1,izo);  vB = vB + temp;
     call d_b2b(-tyT,temp,1,izo); uB = uB + temp;
  endif
  if (bot .and. top .or. correct) then! upper velocity from lower boundary
     call d_b2b(txB,temp,1,-izo);  vT = vT + temp;
     call d_b2b(-tyB,temp,1,-izo); uT = uT + temp;
  endif
  ! FFT back to xy grid.
  if (bot .or. correct) then 
     call ft_2d(txB,mmax,nmax,1) ! t_x
     call ft_2d(vB,mmax,nmax,1) ! v0
     call ft_2d(uB,mmax,nmax,1) ! u0
     call ft_2d(tyB,mmax,nmax,1) ! t_y
  endif
  if (top .or. correct) then 
     call ft_2d(txT,mmax,nmax,1) ! t_x
     call ft_2d(tyT,mmax,nmax,1) ! t_y
     call ft_2d(vT,mmax,nmax,1) ! v0
     call ft_2d(uT,mmax,nmax,1) ! u0         
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!! sQG +1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (correct) then
     szzspB=0.;szspB=0.;u1spB=0.;v1spB=0. ! these _MUST_ be zeroed.
     szzspT=0.;szspT=0.;u1spT=0.;v1spT=0. 
     u1spB=0.;u1spT=0.;v1spB=0.;v1spT=0.
     u1B=0.;u1T=0.;v1B=0.;v1T=0.

     ! basic-state terms: add periodic parts here; linear parts are handled last.
     if (hw) then 
        tspB = tspB + tbB; tspT = tspT + tbT
        uB = uB + ulinB; uT = uT + ulinT - (lam*H) ! don't double count!
        tyB = tyB + tbyB; tyT = tyT + tbyT
     endif

     ! big-grid theta (szsp) and theta_z (szzsp):
     call d_s2b(tspB,szspB,1,Id);  call ft_2d(szspB,mmax,nmax,1)
     call d_s2b(tspT,szspT,1,Id);  call ft_2d(szspT,mmax,nmax,1)

     call d_s2b(tspB,szzspB,1,-dz); call d_s2b(tspT,temp,1,dzo); 
     szzspB = szzspB + temp
     call ft_2d(szzspB,mmax,nmax,1)
     call d_s2b(tspT,szzspT,1,dz); call d_s2b(tspB,temp,1,-dzo); 
     szzspT = szzspT + temp
     call ft_2d(szzspT,mmax,nmax,1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Homogeneous contributions (Note: xy_to_sp dealiases; no d_b2s)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Compute u1 = -F1_z contributions (nonlinearities)
     temp = -uB*szspB
     call xy_to_sp(temp,tempspB,mmax,nmax,kmax,lmax)
     temp = -uT*szspT
     call xy_to_sp(temp,tempspT,mmax,nmax,kmax,lmax)
     ! u1 bottom
     call d_s2s(tempspB,htempB,1,dz)
     call d_s2s(tempspT,htempT,1,-dzo)
     u1spB = - (htempB + htempT) ! -F1_z
     ! u1 toppom
     call d_s2s(tempspB,htempB,1,dzo)
     call d_s2s(tempspT,htempT,1,-dz)
     u1spT = - (htempB + htempT) ! -F1_z

     ! Compute v1 = -G1_z contributions (nonlinearities)
     temp = vB*szspB
     call xy_to_sp(temp,tempspB,mmax,nmax,kmax,lmax)
     temp = vT*szspT
     call xy_to_sp(temp,tempspT,mmax,nmax,kmax,lmax)
     ! v1 bottom
     call d_s2s(tempspB,htempB,1,-dz)
     call d_s2s(tempspT,htempT,1,dzo)
     v1spB = - (htempB + htempT) ! -G1_z
     ! v1 toppom
     call d_s2s(tempspB,htempB,1,-dzo)
     call d_s2s(tempspT,htempT,1,dz)
     v1spT = - (htempB + htempT) ! -G1_z

     ! Compute Phi1 contributions (nonlinearities)
     temp = szzspB*szspB
     call xy_to_sp(temp,tempspB,mmax,nmax,kmax,lmax)
     temp = szzspT*szspT
     call xy_to_sp(temp,tempspT,mmax,nmax,kmax,lmax)
     ! u1 bottom
     call d_s2s(tempspB,htempB,1,iz*dy)
     call d_s2s(tempspT,htempT,1,-izo*dy)
     u1spB = u1spB - (htempB + htempT) ! add -P1_y
     ! u1 toppom
     call d_s2s(tempspB,htempB,1,izo*dy)
     call d_s2s(tempspT,htempT,1,-iz*dy)
     u1spT = u1spT - (htempB + htempT) ! add -P1_y
     ! v1 bottom
     call d_s2s(tempspB,htempB,1,iz*dx)
     call d_s2s(tempspT,htempT,1,-izo*dx)
     v1spB = v1spB + (htempB + htempT) ! add P1_x
     ! v1 toppom
     call d_s2s(tempspB,htempB,1,izo*dx)
     call d_s2s(tempspT,htempT,1,-iz*dx)
     v1spT = v1spT + (htempB + htempT) ! add P1_x

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Particular contributions 
     ! Notes: \Phi contrib gives 2 factor!!!
     !        Dealiasing is crucial here (no shortcuts, DJM!)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! u1 bottom
     temp = (-2.*tyB*szspB) + (uB*szzspB)
     call xy_to_sp(temp,tempspB,mmax,nmax,kmax,lmax)
     u1spB = u1spB + tempspB !add (-P1_y - F1_z)
     ! u1 toppom
     temp = (-2.*tyT*szspT) + (uT*szzspT)
     call xy_to_sp(temp,tempspT,mmax,nmax,kmax,lmax)
     u1spT = u1spT + tempspT !add (-P1_y - F1_z)
     ! v1 bottom
     temp = (2.*txB*szspB) + (vB*szzspB)
     call xy_to_sp(temp,tempspB,mmax,nmax,kmax,lmax)
     v1spB = v1spB + tempspB !add (P1_x - G1_z)
     ! v1 toppom
     temp = (2.*txT*szspT) + (vT*szzspT)
     call xy_to_sp(temp,tempspT,mmax,nmax,kmax,lmax)
     v1spT = v1spT + tempspT !add (P1_x - G1_z)

     ! linear-shear contributions (note: recycling var names)
     if (hw) then 
        call d_s2s(tspB,tempspB,1,-izo*dx*dx) !phi_xx^B
        call d_s2s(tspT,tempspT,1,izo*dx*dx) !phi_xx^T
        call d_s2s(tspB,htempB,1,-izo*dy*dy)  !phi_yy^B
        call d_s2s(tspT,htempT,1,izo*dy*dy)  !phi_yy^T
        call d_s2s(tspB,tempxyB,1,-izo*dx*dy) !phi_xy^B
        call d_s2s(tspT,tempxyT,1,izo*dx*dy) !phi_xy^T

        !         u1spB = u1spB - lam*(H*(htempT - tempspT) + (tspT - tspB) )
        !         u1spT = u1spT + lam*(H*(htempB - tempspB) + (tspT - tspB) )
        ! as per DM email 12/24/02:
        u1spB = u1spB - lam*(H*(htempT - tempspT) - tspB )
        u1spT = u1spT + lam*(H*(htempB - tempspB) + tspT )
        v1spB = v1spB + lam*2*H*tempxyT
        v1spT = v1spT - lam*2*H*tempxyB
     endif

     ! map from small to large spectral array
     call d_s2b(u1spB,u1B,1,Id); call d_s2b(v1spB,v1B,1,Id)
     call d_s2b(u1spT,u1T,1,Id); call d_s2b(v1spT,v1T,1,Id)

     ! FFT u1 and v1 to grid point space
     call ft_2d(u1B,mmax,nmax,1); call ft_2d(v1B,mmax,nmax,1)
     call ft_2d(u1T,mmax,nmax,1); call ft_2d(v1T,mmax,nmax,1)

     ! remove periodic base state (added above):
     uB = uB - ulinB; uT = uT - ulinT + (lam*H)
     tyB = tyB - tbyB; tyT = tyT - tbyT
  endif ! correct

  ! return u = u_0 + Ro * u_1; v = v_0 + Ro * v_1
  ! (look into using the f90 "transfer" function to return real parts...)

  txBr = real(txB); tyBr = real(tyB);
  txTr = real(txT); tyTr = real(tyT);      
  uBr = real(uB + Ross*u1B); uTr = real(uT + Ross*u1T)
  vBr = real(vB + Ross*v1B); vTr = real(vT + Ross*v1T)

  ! Ekman layer calculations (need streamfunction and Laplacian.
  if (gamma .gt. 0) then
     ! reset these since they were changed for next-order calculations
     tspB = itspB; tspT = itspT ! local copies of spectral boundary theta
     sb = 0. ! surface O(1) streamfunction
     call d_s2s(tspB,sB,1,-iz) ! bottom contribution
     call d_s2s(tspT,temps,1,izo) ! toppom contribution
     sB = sB + temps;
     ! now compute Laplacian from previous time-step's sB
     temp = 0.
     call d_s2b(sold,temp,1,dx*dx + dy*dy)
     call ft_2d(temp,mmax,nmax,1)
     sblre = real(temp)
  endif

  return
end SUBROUTINE invert

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE xy_to_sp(xy,sp,mx,ny,km,lm)
  !     Originator: G. Hakim, NCAR/MMM

  ! Map an (x,y) array onto a _smaller_ spectral array.
  ! Input: xy(mx,ny) --- a grid point array.
  ! Output: sp(2*km,2*lm) --- a spectral array.

  IMPLICIT NONE

  integer, intent(in) :: km,lm,mx,ny
  !      real, intent(in) :: xy(mx,ny)
  complex, intent(in) :: xy(mx,ny)
  complex, dimension(2*km,2*lm), intent(out) :: sp
  complex, dimension(mx,ny) :: cop
  real :: rtemp,ctemp
  integer :: i,j,k,l,kk,kmp1,lmp1,k2,l2,LL,MM,stat

  ! initialize arrays:
  sp = 0.; cop = xy

  kmp1 = km + 1; lmp1 = lm + 1; k2 = mx - km; l2 = ny - lm

  call ft_2d(cop,mx,ny,-1)

  !      do k=1,kmp1; do l=1,lmp1
  do k=1,km; do l=1,lm
     kk = km + k; ll = lm + l

     ! waves: 0 <= k <= +/-kmax; 0 <= l <= +/-lmax:
     sp(k,l) = cop(k,l)/real(mx*ny)

     ! +/-kmax <= k <= -1; 0 <= l <= +/-lmax:
     if (k .gt. 1) then
        sp(kk,l) = cop(k2+k,l)/real(mx*ny)
     endif

     ! 0 <= k <= +/-kmax; +/-lmax <= l <= -1:
     if (l .gt. 1) then
        sp(k,ll) = cop(k,l2+l)/real(mx*ny)
     endif

     ! +/-kmax <= k <= -1; +/-lmax <= l <= -1:
     if ((k .gt. 1) .and. (l .gt. 1)) then
        sp(kk,ll) = cop(k2+k,l2+l)/real(mx*ny)
     endif
  enddo; enddo

  return
end SUBROUTINE xy_to_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE sp_to_xy(sp,xy2,km,lm,mx,ny)
  !     Originator: G. Hakim, NCAR/MMM

  ! Map an (km,lm) spectral array onto a _bigger_ grid point array.
  ! Input: sp(2*km,2*lm) --- a spectral array.
  ! Output: xy(mx,ny) --- a grid point array.

  IMPLICIT NONE

  integer, intent(in) :: km,lm,mx,ny
  real, dimension(mx,ny), intent(out) :: xy2
  complex :: xy(mx,ny)
  complex, intent(in) :: sp(2*km,2*lm)
  real rtemp,ctemp
  integer k,l,kk,ll,kmp1,lmp1,k2,l2,MM,stat

  xy = 0.

  kmp1 = km + 1; lmp1 = lm + 1; k2 = mx - km; l2 = ny - lm

  !      do 10 k = 1,kmp1; do 10 l = 1,lmp1
  do k = 1,km
     do l = 1,lm

        ! waves: 0 <= k <= +/-kmax; 0 <= l <= +/-lmax:
        xy(k,l) = sp(k,l)

        kk = km + k; ll = lm + l

        ! +/-kmax <= k <= -1; 0 <= l <= +/-lmax:
        if (k .ge. 1) then
           xy(k2+k,l) = sp(kk,l)
        endif

        ! 0 <= k <= +/-kmax; +/-lmax <= l <= -1:
        if (l .ge. 1) then
           xy(k,l2+l) = sp(k,ll)
        endif

        ! +/-kmax <= k <= -1; +/-lmax <= l <= -1:
        if ((k .ge. 1) .and. (l .ge. 1)) then
           xy(k2+k,l2+l) = sp(kk,ll)
        endif
     enddo !k
  enddo !l

  call ft_2d(xy,mx,ny,1)
  xy2 = real(xy)

  return
end SUBROUTINE sp_to_xy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE diff(dco)
  !     Originator: G. Hakim, NCAR/MMM

  USE spectral
  IMPLICIT NONE

  ! Compute the diffusion coefficient for del^n diffusion (n even>0).
  ! Tau gives the e-folding time scale for damping modes at the 
  ! Nyquist frequency.

  real, intent(out) :: dco
  real in

  ! 10/21/99 gjh mods...
  !      in = (cmplx(0.,1.)**n)
  in = 1
  !      dco = (((real(facx*kmax)**n) + (real(facy*lmax)**n))*tau)**(-1)
  dco = 1.0/(((real(facx*kmax)**2) &
       &          + (real(facy*lmax)**2))**(n/2)*tau)
  dco = in*dco
  if (dco .lt. 0.) dco = -1.*dco

  if (verbose .gt. 1) print*,'diffusion coefficient:',dco

  return
end SUBROUTINE diff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE init(tB,tT)
  !     Originator: G. Hakim, NCAR/MMM

  USE spectral
  IMPLICIT NONE

  ! Initialize basic state and perturbation fields.

  real, intent(out), dimension(2*kmax,2*lmax) :: tB,tT
  integer :: i,j,twokmax,twolmax,irecl
  real :: fac
  logical, parameter :: matlab = .TRUE. ! matlab input files

  if (verbose .gt. 0)  print*,'initializing theta fields...'
  tB = 0.; tT = 0.

!!!!!!!!!!!!!!!!!!!!!!!!!!! Bottom !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  open(9,file=path//'th_init_B.dat', status='old')
  read(9,*) twokmax; read(9,*) twolmax
  if (grow .or. matlab) then !matlab file
     do i = 1, 2*kmax; do j=1, 2*lmax
        read(9,*) tB(i,j)      
     enddo; enddo
  else ! fortran file
     read(9,*) tB
  endif
  close(9)

!!!!!!!!!!!!!!!!!!!!!!!!!!! Toppom !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  open(9,file=path//'th_init_T.dat', status='old')
  read(9,*) twokmax; read(9,*) twolmax
  if (twokmax .ne. 2*kmax .and. twolmax .ne. 2*lmax) then  
     print*,'x and y resolution from th_init_T.dat do not match SPECTRAL.f90!!! ',twokmax,twolmax
     STOP
  endif
  if (grow .or. matlab) then !matlab file
     do i = 1, 2*kmax; do j=1, 2*lmax
        read(9,*) tT(i,j)      
     enddo; enddo
  else ! fortran file
     read(9,*) tT
  endif
  close(9)

  ! normalize initial amplitude (RMS00)
  !      fac = .15/maxval(abs(tB)); tB = fac*tB; tT = fac*tT
  !      tB = .15*tB/maxval(abs(tB)); tT = .15*tT/maxval(abs(tB))

  !      tT = -tT
  !      tB = -tB
  !      tB = 0.
  !      tT = 0.
  !      tT = tT/1000.
  !      tT = tT/10.
  !      tB = tB/10.

  ! load from a restart file (will need to save told, told2, eventually)
  if (restart) then 
     if (verbose .gt. 1) print*,'!!! USING RESTART FILE !!!'
     irecl = 2*kmax*2*lmax*32
     open(9,file=path//'final+1.dat',access='direct', &
          &        form='unformatted', recl = irecl,status='old')
     read(9,rec=1) tB
     read(9,rec=2) tT
     close(9)
  endif

  return
end SUBROUTINE init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE init_jet(tbB,tbT,tbyB,tbyT,ulinB,ulinT,lam)
  !     Originator: G. Hakim, University of Washington

  USE spectral
  IMPLICIT NONE

  ! Initialize basic state jet.
  ! tbB = periodic basic state theta; bottom boundary (spectral grid).
  ! tbT = periodic basic state theta; top boundary (spectral grid).
  ! tbyB = y derivative of periodic basic state theta; bottom boundary (grid point).
  ! tbyT = y derivative of periodic basic state theta; top boundary (grid point).
  ! ulinB = basic state wind; bottom boundary (grid point).
  ! ulinT = basic state wind; top boundary (grid point).

  complex, intent(out), dimension(2*kmax,2*lmax) :: tbB,tbT
  real, intent(out), dimension(mmax,nmax) :: tbyB,tbyT,ulinB,ulinT
  complex, dimension(mmax,nmax) :: tyB,tyT,uB,uT,temp
  real, intent(out) :: lam
  real, dimension(2*kmax,2*lmax) :: tbxyB,tbxyT
  real :: HW_ubar,HW_theta,HW_thetay
  real :: y,yp,dyy
  integer :: j,icho
  complex,dimension(2*kmax,2*lmax) :: dx,dy,dz,dzo,iz,izo,dzi,Id

  if (verbose .gt. 0) print*,'initializing basic state...'

  ! first determine lamda, the linear shear parameter
  lam = (HW_theta(0.,0.) - HW_theta(YL,0.)) / YL
  if (verbose .gt. 1) print*,'lam = ',lam

  ! spectral variables
  dyy = YL/real(2*lmax)
  do j=1, 2*lmax
     y = real(j-1)*dyy; yp = y - (0.5*(YL - hwp))
     tbxyB(:,j) = HW_theta(y,0.) + lam*yp
     tbxyT(:,j) = HW_theta(y,H) + lam*yp
     !         print*,'periodic theta grid point:',y,yp,tbxyT(1,j)
  enddo
  ! map into spectral space at the same resolution:
  call xy_to_sp(cmplx(tbxyB,0.),tbB,2*kmax,2*lmax,kmax,lmax)
  call xy_to_sp(cmplx(tbxyT,0.),tbT,2*kmax,2*lmax,kmax,lmax)

  ! grid point variables
  dyy = YL/real(nmax)
  do j=1, nmax
     y = real(j-1)*dyy; yp = y - (0.5*(YL - hwp))
     tbyB(:,j) = HW_thetay(y,0.) + lam
     tbyT(:,j) = HW_thetay(y,H) + lam
     ! old U
     !         ulinB(:,j) = HW_ubar(y,0.)
     !         ulinT(:,j) = HW_ubar(y,H)
  enddo

  ! new U: solve numerically given theta
  uB = 0.; uT = 0.
  call d_setup(dx,dy,dz,dzo,iz,izo,Id) ! derivative operators
  call d_s2b(tbB,tyB,1,dy); call d_s2b(tbT,tyT,1,dy)
  call d_b2b(-tyB,uB,1,-iz); call d_b2b(-tyT,uT,1,iz)
  call d_b2b(-tyT,temp,1,izo); uB = uB + temp;
  call d_b2b(-tyB,temp,1,-izo); uT = uT + temp;
  call ft_2d(uB,mmax,nmax,1); call ft_2d(uT,mmax,nmax,1)
  ulinB = real(uB); ulinT = real(uT)
  ulinT = ulinT + lam*H

  return
end SUBROUTINE init_jet

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Numerical Recipes random number generator:
FUNCTION ran1(idum)
  INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
  REAL ran1,AM,EPS,RNMX
  PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836, &
       & NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
  INTEGER j,k,iv(NTAB),iy
  SAVE iv,iy
  DATA iv /NTAB*0/, iy /0/
  if (idum.le.0.or.iy.eq.0) then
     idum=max(-idum,1)
     do j=NTAB+8,1,-1
        k=idum/IQ
        idum=IA*(idum-k*IQ)-IR*k
        if (idum.lt.0) idum=idum+IM
        if (j.le.NTAB) iv(j)=idum
     enddo
     iy=iv(1)
  end if
  k=idum/IQ
  idum=IA*(idum-k*IQ)-IR*k
  if (idum.lt.0) idum=idum+IM
  j=1+iy/NDIV
  iy=iv(j)
  iv(j)=idum
  ran1=min(AM*iy,RNMX)
  return
END FUNCTION ran1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE dump(tspB,tspT,ilam,lam,it,cdf_file,cdfid)
  !     Originator: G. Hakim, NCAR/MMM
  !     Modifier  : Rahul Mahajan, UW July 2005. (netcdf output)

  !     Write theta or streamfunction to disk.
  USE spectral
  IMPLICIT NONE

  complex, dimension(2*kmax,2*lmax) :: tspB,tspT,copy,tbB,tbT
  real, dimension(2*kmax,2*lmax) :: oxyB,oxyT
  real, intent(in) :: lam
  real :: y,yp,dy
  integer :: ilam
  integer, save :: irec
  integer :: i,j
  character(len=16) :: form
  character(len=8) :: size

! Pertaining to netCDF files
  logical head
  integer :: it, cdfid(2)
  character(len=64) cdf_file(2) 

! set up the format form
  write(size,'(i8)') (2*kmax*2*lmax); size = adjustr(size)
  form = '('//size//'F10.6'//')'!; print*,'FORMAT',form

  copy = tspB
  call sp_to_xy(copy,oxyB,kmax,lmax,2*kmax,2*lmax)
  copy = tspT
  call sp_to_xy(copy,oxyT,kmax,lmax,2*kmax,2*lmax)

  ! add in linear shear
  if (ilam .eq. 1) then 
     dy = YL/real(2*lmax)
     do j=1, 2*lmax
! fixed 08/10/2005 GJH & RBM
        y = real(j-1)*dy; 
        oxyB(:,j) = oxyB(:,j) - lam*y
        oxyT(:,j) = oxyT(:,j) - lam*y
     enddo
  endif

  if (verbose .gt. 0) print*,'Writing to disk...'

  if(cdf) then
     ! write to netCDF files
     head = .FALSE. 
     call cdf_dump(cdf_file(1),head,oxyB,it)
     call cdf_dump(cdf_file(2),head,oxyT,it)
  else
     !  write to matlab
     do i=1,2*kmax
        do j=1,2*lmax
           write(12,*) oxyB(i,j)
           write(13,*) oxyT(i,j)
        enddo
     enddo
  endif
  
  ! write the init files if this is a grow run
  if (grow) then 
     open(31,file=path//'th_init_B.dat') 
     write(31,'(i5)') 2*kmax
     write(31,'(i5)') 2*lmax
     write(31,*) oxyB; close(31)
     open(31,file=path//'th_init_T.dat') 
     write(31,'(i5)') 2*kmax
     write(31,'(i5)') 2*lmax
     write(31,*) oxyT; close(31)
  endif

  return
end SUBROUTINE dump

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE advect(u,v,f_x,f_y,fb_y,h_x,h_y,ub,tf,lam,lap)
  !     Originator: G. Hakim, U. Washington

  ! Spectral advection on the (x,y) grid.
  USE spectral
  IMPLICIT NONE

  real, intent(in), dimension(mmax,nmax) :: u,v,f_x,f_y,fb_y,ub,lap,h_x,h_y
  real, intent(in) :: lam
  real, intent(out), dimension(mmax,nmax) :: tf
  real :: terr,ekman,dx,x,dy,y,ran1,famp,rw
  integer :: i,j
  integer, save :: iseed
  real :: rphasex,rphasey

  dx = XL/real((mmax)); dy = YL/real((nmax))

! random wave-one forcing; random walk in phase (30 August 2006)
!  rw = 32. ! random wavenumber
  rw = 4. ! random wavenumber
  famp = 0.1*H ! forcing amplitude
  famp = 0. ! forcing amplitude
  rphasex = 2.*pi*(ran1(iseed)-0.5)
  rphasey = 2.*pi*(ran1(iseed)-0.5)

  do i=1,mmax; do j=1,nmax
     x = real(i-1)*dx; y = real(j-1)*dy

     if (linear) then 
        tf(i,j) = -((v(i,j) * (fb_y(i,j) - lam)) + &
             &                 (ub(i,j) * f_x(i,j)))
     else
        tf(i,j) = -((v(i,j) * (f_y(i,j) + fb_y(i,j) - lam)) + &
             &                 ((u(i,j) + ub(i,j)) * f_x(i,j)))
     endif

     tf(i,j) = tf(i,j) - (ekman(x,y)*lap(i,j)) ! Ekman layer

! terrain:
     if (iterr) then
        terr = -( (v(i,j) * h_y(i,j)) + ((u(i,j) + ub(i,j)) * h_x(i,j)) )
        tf(i,j) = tf(i,j) + terr 
     endif

! random wave-one forcing; random walk in phase (30 August 2006)
     tf(i,j) = tf(i,j) - famp*sin((rw*x*2*pi/XL)-rphasex)*sin((rw*y*2*pi/YL)-rphasey)

  enddo; enddo


  return
end SUBROUTINE advect

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE tadvB(dat,tend,dco,first)
  !     Originator: G. Hakim, NCAR/MMM

  ! Time-advance subroutine.
  USE spectral
  IMPLICIT NONE

  complex, intent(inout), dimension(2*kmax,2*lmax) :: dat
  !f2py intent(in,out) :: dat
  complex, dimension(2*kmax,2*lmax), save :: told,told2
  !      real, intent(in), dimension(mmax,nmax) :: rtend,itend
  complex, intent(in), dimension(mmax,nmax) :: tend
  real, intent(in) :: dco
  real :: ak,bl
  real, save :: dts
  integer :: k,l,kk,ll
  complex :: ttsp,temp
  logical :: first

  if (first) then
     !  adams-bash
     dts = dt*(12./23.); told = 0.; told2 = 0.
  else
     dts = dt
  endif

  do k=1,kmax; do l=1,lmax
     ! 0 <= k <= +/-kmax; 0 <= l <= +/-lmax:
     ak = facx*real(k - 1); bl = facy*real(l - 1)
     ttsp = tend(k,l)/real(mmax*nmax)
     call tstep_ab(ttsp,dat(k,l),told(k,l),told2(k,l), &
          &        ak,bl,dco,dts)

     kk = kmax + k; ll = lmax + l
     ! +/-kmax <= k <= -1; 0 <= l <= +/-lmax:
     !         if (k .gt. 1 .and. k .lt. kmax) then
     if (k .gt. 1) then
        ak = -1.*facx*real(kmaxp1 - k); bl = facy*real(l - 1)
        ttsp = tend(k2+k,l)/real(mmax*nmax)
        call tstep_ab(ttsp,dat(kk,l),told(kk,l), &
             &           told2(kk,l),ak,bl,dco,dts)
     endif
     ! 0 <= k <= +/-kmax; +/-lmax <= l <= -1:
     !         if (l .le. lmax) then
     if (l .gt. 1) then
        ak = facx*real(k-1); bl = -1.*facy*real(lmaxp1 - l)
        ttsp = tend(k,l2+l)/real(mmax*nmax)
        call tstep_ab(ttsp,dat(k,ll),told(k,ll), &
             &           told2(k,ll),ak,bl,dco,dts)
     endif
     ! +/-kmax <= k <= -1; +/-lmax <= l <= -1:
     !         if (k .le. kmax .and. l .le. lmax) then
     if ((k .gt. 1) .and. (l .gt. 1)) then
        ak=-1.*facx*real(kmaxp1 - k); bl=-1.*facy*real(lmaxp1 - l)
        ttsp = tend(k2+k,l2+l)/real(mmax*nmax)
        call tstep_ab(ttsp,dat(kk,ll),told(kk,ll), &
             &         told2(kk,ll),ak,bl,dco,dts)
     endif
  enddo; enddo

  return
end SUBROUTINE tadvB

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE tadvT(dat,tend,dco,first)
  !     Originator: G. Hakim, NCAR/MMM

  ! Time-advance subroutine.
  USE spectral
  IMPLICIT NONE

  complex, intent(inout), dimension(2*kmax,2*lmax) :: dat
  !f2py intent(in,out) :: dat
  complex, dimension(2*kmax,2*lmax), save :: told,told2
  !      real, intent(in), dimension(mmax,nmax) :: rtend,itend
  complex, intent(in), dimension(mmax,nmax) :: tend
  real, intent(in) :: dco
  real :: ak,bl
  integer :: k,l,kk,ll
  complex :: ttsp,temp
  logical :: first
  real, save :: dts

  if (first) then
     !  adams-bash
     dts = dt*(12./23.); told = 0.; told2 = 0.
  else
     dts = dt
  endif

  do k=1,kmax; do l=1,lmax
     ! 0 <= k <= +/-kmax; 0 <= l <= +/-lmax:
     ak = facx*real(k - 1); bl = facy*real(l - 1)
     ttsp = tend(k,l)/real(mmax*nmax)
     call tstep_ab(ttsp,dat(k,l),told(k,l),told2(k,l), &
          &        ak,bl,dco,dts)

     kk = kmax + k; ll = lmax + l
     ! +/-kmax <= k <= -1; 0 <= l <= +/-lmax:
     !         if (k .gt. 1 .and. k .lt. kmax) then
     if (k .gt. 1) then
        ak = -1.*facx*real(kmaxp1 - k); bl = facy*real(l - 1)
        ttsp = tend(k2+k,l)/real(mmax*nmax)
        call tstep_ab(ttsp,dat(kk,l),told(kk,l), &
             &           told2(kk,l),ak,bl,dco,dts)
     endif
     ! 0 <= k <= +/-kmax; +/-lmax <= l <= -1:
     !         if (l .le. lmax) then
     if (l .gt. 1) then
        ak = facx*real(k-1); bl = -1.*facy*real(lmaxp1 - l)
        ttsp = tend(k,l2+l)/real(mmax*nmax)
        call tstep_ab(ttsp,dat(k,ll),told(k,ll), &
             &           told2(k,ll),ak,bl,dco,dts)
     endif
     ! +/-kmax <= k <= -1; +/-lmax <= l <= -1:
     !         if (k .le. kmax .and. l .le. lmax) then
     if ((k .gt. 1) .and. (l .gt. 1)) then
        ak=-1.*facx*real(kmaxp1 - k); bl=-1.*facy*real(lmaxp1 - l)
        ttsp = tend(k2+k,l2+l)/real(mmax*nmax)
        call tstep_ab(ttsp,dat(kk,ll),told(kk,ll), &
             &         told2(kk,ll),ak,bl,dco,dts)
     endif
  enddo; enddo

  return
end SUBROUTINE tadvT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine tstep_ab(tend,dat,told,told2,ak,bl,dco,dts)
  ! 3rd-order Adams-Bashforth time step (Durran 1991).
  USE spectral
  IMPLICIT NONE

  complex, intent(inout) :: dat, told, told2
  complex, intent(in) :: tend
  real, intent(in) :: ak,bl,dco,dts
  complex :: new, temp
  real :: tfac,relax

  ! changed 02/07/03 to include relaxation to jet
  relax = 0.
  if (trl .lt. 1.e3 .and. trl .gt. 0) relax = 1./trl 
  ! 12/17/2008: damps l=0 faster for stability (used to zero in main block)
  if (bl .eq. 0) relax = relax*4

  tfac = dts*((dco*(((ak**2)+(bl**2))**(n/2))) + relax)

  ! new t-step:
  !      tfac = dco*dts*((ak**n)+(bl**n))
  !      tfac = dco*dts*((ak**2)+(bl**2))**(n/2)

  told = told*exp(-1.*tfac)
  told2 = told2*exp(-tfac)

  temp = (23.*tend) - (16.*told) + (5.*told2)
  new = dat + ( (dts/12.)*(temp) )

  !      dat=new*exp(-1.*tfac); told2=told*exp(tfac); told=tend
  dat=new*exp(-1.*tfac); told2=told; told=tend

  return
end subroutine tstep_ab

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ft_2d(f,ni,nj,isign)
  !     Originator: G. J. Hakim, University of Washington.

  ! FFT-calling subroutine. Plug in appropriate calls.

  IMPLICIT NONE

  !      include '/usr/include/DXMLDEF.FOR' ! for the dec fft routines
  !      logical, parameter :: dec = .TRUE. ! for the dec fft routines
  logical, parameter :: dec = .FALSE.; real :: cfft_2d !netlib

  integer, intent(in) :: ni, nj, isign
  real, dimension(ni,nj) :: re, im
  complex, dimension(ni,nj) :: f
  !f2py intent(in,out) :: f
  character(len=1) :: csign
  !  integer*4 :: stat
  integer :: stat

  if (dec) then 
     if (isign .eq. -1) then 
        csign = 'F'
     elseif (isign .eq. 1) then 
        csign = 'B'
        f = f*ni*nj
     else
        print*,'Error in ft_2d call'
        stop
     endif
!     stat = CFFT_2D('C','C',csign,f,f,ni,nj,ni,1,1)
     if (stat .ne. 0) then 
        print*,'Error in dec fft call:',stat
        stop
     endif

  elseif (.not. dec) then 

     re = real(f); im = aimag(f)
     call fft(re,im,ni*nj,ni,ni,isign)
     call fft(re,im,ni*nj,nj,nj*ni,isign)
     f = cmplx(re,im)

  endif

  return
end subroutine ft_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine d_setup(dx,dy,dz,dzo,iz,izo,Id)
  !     Originator: G. Hakim, University of Washington.

  ! Set up matrices for derivatives and integrals.

  USE spectral
  IMPLICIT NONE

  complex, dimension(2*kmax,2*lmax) :: dx,dy,dz,dzo,iz,izo,Id
  real :: m
  integer :: k,l

  ! dx
  dx(1,:) = 0.; dx(kmax+1,:) = 0.
  do k=2,kmax
     dx(k,:) = facx*real(k-1)*cmplx(0.,1.)
     dx(2*kmax-k+2,:) = -1.*dx(k,:)
  enddo

  ! dy
  dy(:,1) = 0.; dy(:,lmax+1) = 0.
  do l=2,lmax
     dy(:,l) = facy*real(l-1)*cmplx(0.,1.)
     dy(:,2*lmax-l+2) = -1.*dy(:,l)
  enddo

  ! dz,dzi
  do k=1,2*kmax; do l=1,2*lmax
     if (k .eq. 1 .and. l .eq. 1) then 
        dz(k,l) = 0.; dzo(k,l) = 0.; iz(k,l) = 0.; izo(k,l) = 0.
     elseif (k .eq. 1 .and. l .eq. lmax+1) then
        dz(k,l) = 0.; dzo(k,l) = 0.; iz(k,l) = 0.; izo(k,l) = 0.
     elseif (k .eq. kmax+1 .and. l .eq. lmax+1) then
        dz(k,l) = 0.; dzo(k,l) = 0.; iz(k,l) = 0.; izo(k,l) = 0.
     elseif (l .eq. 1 .and. k .eq. kmax+1) then
        dz(k,l) = 0.; dzo(k,l) = 0.; iz(k,l) = 0.; izo(k,l) = 0.
     else
        ! must force real, since small complex residual may remain:
        m = real(((-(dx(k,1)**2)) - (dy(1,l)**2))**0.5)
        dz(k,l) = m / tanh(m*H) ! h->\infty: tanh->1, dz->m
        iz(k,l) = 1./ (m * tanh(m*H)) ! h->\infty: iz-> 1/m
        ! operators for opposing boundary
        if (m*H .lt. 50 .and. model .ne. 0) then
           dzo(k,l) = m / sinh(m*H) ! h->\infty: dzo->0
           izo(k,l) = 1./ (m * sinh(m*H)) ! h->\infty: izo->0
        elseif (model .eq. 0) then 
           print *,'SELECTED 2D EULER DYNAMICS....'
           dz(k,l) = m**2
           iz(k,l) = 1./dz(k,l)
        else
           dzo(k,l) = 0.
           izo(k,l) = 0.
        endif
     endif
  enddo; enddo

  ! Id
  Id = 1.0

  return
end subroutine d_setup

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine d_s2b(in,out,dflag,dn)
  !     Originator: G. J. Hakim, University of Washington.

  ! small array to big array.
  ! dflag =  n: n derivatives. dflag = -n: n integrations.

  USE spectral
  IMPLICIT NONE

  complex, intent(in), dimension(2*kmax,2*lmax) :: in,dn
  integer, intent(in) :: dflag
  complex, intent(out), dimension(mmax,nmax) :: out
  integer :: k,l,kk,ll

  do k = 1,kmax; do l = 1,lmax
     ! waves: 0 <= k <= +/-kmax; 0 <= l <= +/-lmax:

     out(k,l) = (dn(k,l)**dflag)*in(k,l)

     kk = kmax + k; ll = lmax + l

     ! +/-kmax <= k <= -1; 0 <= l <= +/-lmax:
     if (k .gt. 1) then
        out(k2+k,l) = (dn(kk,l)**dflag)*in(kk,l)
     endif

     ! 0 <= k <= +/-kmax; +/-lmax <= l <= -1:
     if (l .gt. 1) then
        out(k,l2+l) = (dn(k,ll)**dflag)*in(k,ll)
     endif

     ! +/-kmax <= k <= -1; +/-lmax <= l <= -1:
     if ((k .gt. 1) .and. (l .gt. 1)) then
        out(k2+k,l2+l) = (dn(kk,ll)**dflag)*in(kk,ll)
     endif
  enddo; enddo

  return
end subroutine d_s2b

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine d_s2s(in,out,dflag,dn)
  !     Originator: G. J. Hakim, University of Washington.

  ! small array to small array.
  ! dflag =  n: n derivatives. dflag = -n: n integrations.

  USE spectral
  IMPLICIT NONE

  complex, intent(in), dimension(2*kmax,2*lmax) :: in,dn
  integer, intent(in) :: dflag
  complex, intent(out), dimension(2*kmax,2*lmax) :: out
  integer :: k,l

  out = 0.

  do k = 1,2*kmax; do l = 1,2*lmax
     if (dn(k,l) .ne. 0) out(k,l) = (dn(k,l)**dflag)*in(k,l)
  enddo; enddo

  return
end subroutine d_s2s

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine d_b2b(in,out,dflag,dn)
  !     Originator: G. J. Hakim, University of Washington.

  ! big array to big array.
  ! dflag =  n: n derivatives. dflag = -n: n integrations.

  USE spectral
  IMPLICIT NONE

  complex, intent(in), dimension(mmax,nmax) :: in
  complex, intent(in), dimension(2*kmax,2*lmax) :: dn
  integer, intent(in) :: dflag
  complex, intent(out), dimension(mmax,nmax) :: out
  integer :: k,l,kk,ll

  out = 0.

  do k = 1,kmax; do l = 1,lmax
     ! waves: 0 <= k <= +/-kmax; 0 <= l <= +/-lmax:

     if (dn(k,l) .ne. 0) out(k,l) = (dn(k,l)**dflag)*in(k,l)

     kk = kmax + k; ll = lmax + l

     ! +/-kmax <= k <= -1; 0 <= l <= +/-lmax:
     if (k .gt. 1) then
        if (dn(kk,l) .ne. 0)  &
             &           out(k2+k,l) = (dn(kk,l)**dflag)*in(k2+k,l)
     endif

     ! 0 <= k <= +/-kmax; +/-lmax <= l <= -1:
     if (l .gt. 1) then
        if (dn(k,ll) .ne. 0) &
             &           out(k,l2+l) = (dn(k,ll)**dflag)*in(k,l2+l)
     endif

     ! +/-kmax <= k <= -1; +/-lmax <= l <= -1:
     if ((k .gt. 1) .and. (l .gt. 1)) then
        if (dn(kk,ll) .ne. 0)  &
             &           out(k2+k,l2+l) = (dn(kk,ll)**dflag)*in(k2+k,l2+l)
     endif
  enddo; enddo

  return
end subroutine d_b2b

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine d_b2s(in,out,dflag,dn)
  !     Originator: G. J. Hakim, University of Washington.

  ! big array to small array.
  ! dflag =  n: n derivatives. dflag = -n: n integrations.

  USE spectral
  IMPLICIT NONE

  complex, intent(in), dimension(mmax,nmax) :: in
  complex, intent(in), dimension(2*kmax,2*lmax) :: dn
  integer, intent(in) :: dflag
  complex, intent(out), dimension(2*kmax,2*lmax) :: out
  integer :: k,l,kk,ll

  do k = 1,kmax; do l = 1,lmax
     ! waves: 0 <= k <= +/-kmax; 0 <= l <= +/-lmax:

     out(k,l) = (dn(k,l)**dflag)*in(k,l)

     kk = kmax + k; ll = lmax + l

     ! +/-kmax <= k <= -1; 0 <= l <= +/-lmax:
     if (k .ge. 1) then
        out(kk,l) = (dn(kk,l)**dflag)*in(k2+k,l)
     endif

     ! 0 <= k <= +/-kmax; +/-lmax <= l <= -1:
     if (l .ge. 1) then
        out(k,ll) = (dn(k,ll)**dflag)*in(k,l2+l)
     endif

     ! +/-kmax <= k <= -1; +/-lmax <= l <= -1:
     if ((k .ge. 1) .and. (l .ge. 1)) then
        out(kk,ll) = (dn(kk,ll)**dflag)*in(k2+k,l2+l)
     endif
  enddo; enddo

  return
end subroutine d_b2s

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! BEGIN: Hoskins-West base state function calls
REAL FUNCTION HW_ubar(y,z)
  USE spectral
  IMPLICIT NONE
  real y,yp,z

  yp = y - (0.5*(YL - hwp))
  if (yp .lt. 0. .and. amu .ne. 0) yp = 0
  if (yp .gt. hwp .and. amu .ne. 0) yp = hwp
  HW_UBAR = shear*z - &
       &     (amu/2.)*(z + ((sinh(ryl*z)/sinh(ryl))*cos(yp*ryl)))

  return
end FUNCTION HW_ubar

REAL FUNCTION HW_theta(y,z)
  ! 9/2006: added amp and correction factor for HW 'hiccup'
  ! NOW INCOMPATIBLE WITH HW_ubar!!! MUST solve for U numerically!!!
  USE spectral
  IMPLICIT NONE
  real y,yp,z
  real amp,hiccup

  amp = 1.5; ! control jet strength
  hiccup = 0.75;
  yp = y - (0.5*(YL - hwp))
  if (yp .lt. 0. .and. amu .ne. 0) yp = 0
  if (yp .gt. hwp .and. amu .ne. 0) yp = hwp
  HW_theta = (-shear*yp)+(amu/2.)* &
       &     (yp+(hiccup*(cosh(ryl*z)/sinh(ryl))*sin(yp*ryl)))
  HW_theta = amp*HW_theta

  return
end FUNCTION HW_theta

real function HW_THETAY(y,z)
  USE spectral
  IMPLICIT NONE
  real y,yp,z

  yp = y - (0.5*(YL - hwp))
  ! fixed this block 02/20/03

  if (amu .eq. 0) then
     HW_THETAY = -shear
  else
     if (yp .lt. 0 .or. yp .gt. hwp) then 
        HW_THETAY = 0.
     else
        HW_THETAY = -shear + (amu/2.) + &
          &        ((amu*ryl/2.)*(cosh(ryl*z)*cos(yp*ryl))/sinh(ryl))
     endif
  endif
  return
end function HW_THETAY
! END: Hoskins-West base state function calls


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! hacked from qg code to check mode growth rates
SUBROUTINE norm(tb,tt,itime)
  USE spectral

  IMPLICIT NONE

  complex, intent(in), dimension(2*kmax,2*lmax) :: tb,tt
  integer, intent(in) :: itime
  real :: V,ak,bl,Vgr,ttau,cmax
  real, save :: V_old,V_0
  integer :: i,j,k,l,m,kk
  logical, parameter :: mesg = .TRUE.
  !      logical, parameter :: mesg = .FALSE.

  ! enstrophy norm (for normal modes, any norm will do)
  V = (sum(abs(tb)**2) + sum(abs(tt)**2))
  V = 0.5*(V**0.5)

  ttau = real(itime -5)*dt
  if (itime .eq. 5) then ! wait 5 time steps for norms to stabilize.
     V_0 = V
  endif

  if (itime .gt. 2) then
     !         Vgr = ((V - V_old)/(dt*0.5*(V+V_old)))
     Vgr = (log(V) - log(V_old))/dt
  endif

  if (mesg) print*,'V inst growth rates:',itime*dt,Vgr
  !      if (itime .gt. 5) then 
  !         if (mesg) print*,'V mean growth rates:',(alog(V/V_0))/ttau
  !      endif

  V_old = V

  return
end SUBROUTINE norm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION ekman(x,y)
  USE spectral

  real x,y

  if (x .ge. 2.*XL/3.) then
     ekman = gamma
  else
     ekman = 0.
     ekman = gamma
     !         ekman = gamma/4.
  endif

  ekman = gamma
  !      if (abs(y-(YL/2)) .ge. 3) ekman = 0.5

  return
end FUNCTION ekman

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cdf_dump(output_file,head,a,it)
  !     Originator: R. Mahajan, University of Washington.

  use spectral

  implicit none
  
	include 'netcdf.inc'
  
  character(len=64) output_file
  integer nx, ny, cdfid, vardim(4), avar, timevar, domvar, domdim
  integer it, k, count(4), start(4), ierr
  real time
  real, dimension(2*kmax,2*lmax) :: a
  logical head   ! to write the header/data part of the cdf file
  
  !     Main Program starts
  
  nx=2*kmax; ny=2*lmax
  
  if(head) then
     
     !           Create a new NetCDF file
     cdfid = nccre(output_file, NCCLOB, ierr)
     
     !           Define dimensions
     vardim(1) = ncddef(cdfid, 'nx', nx, ierr)
     vardim(2) = ncddef(cdfid, 'ny', ny, ierr)
     vardim(3) = ncddef(cdfid, 'nz', 1, ierr)
     vardim(4) = ncddef(cdfid, 'time', NCUNLIM, ierr)
     
     domdim = ncddef(cdfid, 'xy', 2, ierr)
     
     !           Define variables
     avar = ncvdef(cdfid, 'theta', NCFLOAT, 4, vardim, ierr)
     timevar = ncvdef(cdfid, 'time', NCFLOAT, 1, vardim(4), ierr)
     
     domvar = ncvdef(cdfid, 'domain', NCFLOAT, 1, domdim, ierr)
     
     call ncendf(cdfid, ierr) ! End definitions
     
     ! Feed in the dimensions of the domain
     call ncvpt1(cdfid, domvar, 1, XL, ierr)
     call ncvpt1(cdfid, domvar, 2, YL, ierr)
     
  else 
     
     ierr = nf_open(output_file, 1, cdfid)
     time = (it-1)*dt
     k = 1 + (it-1)/iplot
     
     !           Get back variable ID's from the netCDF file
     timevar = ncvid(cdfid, 'time', ierr) 
     avar = ncvid(cdfid, 'theta', ierr) 
     
     call ncvpt1(cdfid, timevar, k, time, ierr)
     
     count(1) = nx;  start(1) = 1
     count(2) = ny;  start(2) = 1
     count(3) = 1;   start(3) = 1
     count(4) = 1;   start(4) = k
     
     call ncvpt(cdfid, avar, start, count, a, ierr)
     
  endif
  
  call ncclos(cdfid,ierr)

  return
end subroutine cdf_dump

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Terrain and barotropic wind on the _advection_ grid:
SUBROUTINE terrain(hamp,hx,hy,hu,hv)

!     Originator: G. Hakim, NCAR/MMM
! feb 2005: added tropo winds from z = 0 for Hsqg

  USE spectral
  IMPLICIT NONE
  
  complex, dimension(mmax,nmax) :: hxsp,hysp
  real, dimension(mmax,nmax), intent(out) :: hx,hy,hu,hv
  complex, dimension(mmax,nmax) :: hh,hxs,hys
  real, dimension(2*kmax,2*lmax) :: hspR
  complex, dimension(2*kmax,2*lmax) :: hspC
  real :: ubaro(nmax),hamp,x,y,ddx,ddy,xcen,ycen,amdx,amdy,asx,asy
  real :: asig,rr,c1,c2,ak,bl,lam
  integer :: i,j,k,l,LL,MM,kk,mmd2,nmd2,mmaxp1,nmaxp1
  real, dimension(mmax,nmax) :: Rblank
  complex, dimension(2*kmax,2*lmax) :: Cblank
  ! fool's day 2009 for angie
  complex,dimension(2*kmax,2*lmax),save:: dx,dy,dz,dzo,iz,izo,dzi,Id

  ! barotropic wind
  hu = 0.; hv = 0.
  
  ddx = XL/real((mmax)); ddy = YL/real((nmax))
  xcen = mmax*ddx/2; ycen = nmax*ddy/2
  
  print*,'xcen,ycen:',xcen,ycen
  
  do i=1,mmax; do j=1,nmax
     x = (i-1)*ddx; y = (j-1)*ddy
     amdx=min(abs(x-xcen),abs(XL+x-xcen),abs(XL-x+xcen))
     amdy=min(abs(y-ycen),abs(YL+y-ycen),abs(YL-y+ycen))
     asx = 1.0; asy = 1.0; asig = 1.0
     rr = (((amdx/asx)**2) + ((amdy/asy)**2))**.5
     ! gaussian topography:
     hh(i,j) = hamp*exp(-1.*((rr/asig)**2))
  enddo; enddo
  
  print*,real(hh(:,nmax/2))

  call ft_2d(hh,mmax,nmax,-1)
  hh = hh / (mmax*nmax)

  ! form spectral h_x, h_y:
  ! fool's day 2009 for angie
  call d_setup(dx,dy,dz,dzo,iz,izo,Id)    ! derivative operators
  call d_b2b(hh,hxs,1,dx) ! x derivative: small to big domain.
  call d_b2b(hh,hys,1,dy) ! y derivative: small to big domain.

  ! back to grid point space
  call ft_2d(hxs,mmax,nmax,1)
  call ft_2d(hys,mmax,nmax,1)
  hx = real(hxs)
  hy = real(hys)

  print*,hx(:,nmax/2)

  if (model .eq. 4) then ! HsQG needs topo winds on the tropo
     Rblank = 0.; Cblank = 0. ! for HsQG inversion call
  
     !first make spectral topo height on small grid
     dx = XL/real((2*kmax)); dy = YL/real((2*lmax))
     do i=1,2*kmax; do j=1,2*lmax
        x = (i-1)*ddx; y = (j-1)*ddy
        amdx=min(abs(x-xcen),abs(XL+x-xcen),abs(XL-x+xcen))
        amdy=min(abs(y-ycen),abs(YL+y-ycen),abs(YL-y+ycen))
        asx = 1.0; asy = 1.0; asig = 1.0
        rr = (((amdx/asx)**2) + ((amdy/asy)**2))**.5
        hspr(i,j) = hamp*exp(-1.*((rr/asig)**2))
     enddo; enddo
     print*,'max topo=',maxval(abs(hspR))
     call xy_to_sp(cmplx(hspR,0.),hspC,2*kmax,2*lmax,kmax,lmax)
     print*,'max topo spectral=',maxval(abs(hspC))
     call invert(-hspC,Cblank,Rblank,Rblank,Rblank,Rblank,Rblank, & 
          & hv,Rblank,hu,Cblank,Cblank,Rblank,Rblank,Rblank, & 
          & Rblank,.TRUE.,.TRUE.,.TRUE.,lam,Cblank,Cblank,Rblank)
     ! hu and hv have the tropo winds due to topography
     print*,'max tropo winds due to topography: ',maxval(abs(hu)),maxval(abs(hv))
  endif
  
  return
end SUBROUTINE terrain




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2sQG+1 model version log
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! 15 February 2002
!
! Log:
!
! v.2.0:
! 1) We cannot prove the linear shear terms are correct, but the comparison 
!    with RMS is favorable. Nonlinear wave tests were inconclusive.
!
! 2) Errors in the inversion routine have been corrected. Cross-boundary
!    terms were handled incorrectly and went undetected in tests because 
!    they go to zero in the limiting cases.
!
! 13 March 2003
!
! v.2.1:
! 1) Added option for wider domain.
!
! 03 December 2003
!
! v.2.2:
! 1) Added an Ekman layer; formally valid for leading-order only.
!
! 2) Removed all multiple loops sharing continuation lines to conform 
!    to f90/95 standards.
!
! 3) Activated Rayleigh damping.
!
! 19 August 2004
!
! v.2.3:
! 1) renamed the code with dashes instead of underscores, and .f90 
!    instead of .f. For one thing, emacs, recognizes the difference.
! 2) extracted SPECTRAL module; it is not compiled separately.
! 3) the code was reformatted to "standard" f90 form, rather than 
!    f77 form for compliance with f2py. 
!
! v.2.4:
! 8 February 2005
!
! 1) added terrain
!
! v.2.5
! 17 February 2005
!
! 1) new terrain call. new ic2 (v2.0) supplies terrain through oxyB.
!
! v.2.6
! 14 July 2005
! 1) netcdf output files enabled (Rahul Mahajan)
!
! v2.6.1
! 19 July 2005
! 1) fixed minor bug in subroutine dump to catch cdf or matlab output.
!
! v2.6.2
! 10 August 2005
! 1) minor bug fix in dump for linear gradient.
!
! v2.7
! 6 September 2005
! 1) convert main source into a subroutine for rolv
! 2) pull out diffusion time scale (tau) to SPECTRAL module.
!
! v2.7.1
! 1) minor revs for netcdf include and print of Ross #.
!
! **** check svn logs beyond this point. ****
