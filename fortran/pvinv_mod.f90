!============================================================
! <next few lines under version control, D O  N O T  E D I T>
! $Date$
! $Author$
! $Revision$
! $Id$
!============================================================

!
! source: svn+ssh://modon.atmos.washington.edu/home/disk/modon/hakim/svn/pvinv
!
! pvinv.f --- a QG+1 pv inversion code.
!
! Originator: G. J. Hakim, University of Washington; hakim@atmos.washington.edu
!
! revision log:
!
! 17 November 2005: revision 44 was tagged as release v2.0.
! 23 Septmeber 2005: revision 12 was tagged as release v1.0.
! 20 September 2005: boundary fields returned; memory issues fixed.
! ----subversion control used beyond this time
! 14 September 2005: v0.7.1 Making output data of dump routine 3 dimensional,
! 09 September 2005: v0.7 calc u,v on boundaries to next order; write to files.
! 08 Septmeber 2005: v. 0.6.2 this is Rolv's version---mod this for v.0.7
! 20 August 2004: v. 0.6 extracted spectral.f90; f90 format; epv routine.
! 3 August 2004: v. 0.5 next-order corrected boundary winds
! 30 July 2004: v. 0.4 option to read ic2.m output for ICs.
! 29 July 2004: v. 0.3.1 write out pv.
! 26 July 2004: v. 0.3 minor bug fix.
! 25 July 2004: v. 0.2 (dimensional values; theta output) 
! 14 July 2004: v. 0.1
!
! Compile: f90 pvinv.f fft_90.f (or, using Makefile: make pvinv)
!

PROGRAM pvinv
  call program()
  stop
END PROGRAM pvinv

subroutine program()

  USE spectral
  IMPLICIT NONE
  
  complex, dimension(2*kmax,2*lmax) :: tb,tt,phib,phit,ub,vb,ut,vt
  complex, dimension(2*kmax,2*lmax,pmax) :: pvsp
  complex, dimension(2*kmax,2*lmax,pmax) :: u,v,w,phi,theta
  real, dimension(2*kmax,2*lmax) :: itbxy,ittxy
  real, dimension(2*kmax,2*lmax,pmax) :: ipvxy
  ! dimensional real fields
  real, dimension(pmax) :: pheight
  real, dimension(2*kmax,2*lmax,pmax) :: ophi,ou,ov,ow,otheta,opv
  real, dimension(2*kmax,2*lmax) :: otb,ott,ophib,ophit,oub,ovb,out,uvt

  ! set up initial condition in (x,y) space
  if (verbose .gt. 0) then
     print*,' '; print*,'initializing pv and theta fields...'
  endif
  ! second entry to init: 0 = define IC in subroutine; 1 = read matlab file
  call init(pmax,1,itbxy,ittxy,ipvxy)

  ! now send to the main inversion routine
  call main(pmax,itbxy,ittxy,ipvxy,phi,theta,u,v,w, &
     &                        tb,tt,phib,phit,ub,vb,ut,vt,pvsp)

  if (verbose .gt. 0) then 
     print*,' '; print*,'writing results to file...'
  endif
  call dump(1,pmax,phi,theta,u,v,w,tb,tt,phib,phit,ub,vb,ut,vt,pvsp, &
  &ophi,otheta,ou,ov,ow,otb,ott,ophib,ophit,oub,ovb,out,uvt,opv,pheight)

  stop

end subroutine program

subroutine main(nlevs,itbxy,ittxy,ipvxy,phi,theta,u,v,w, &
     &                           tbsp,ttsp,phib,phit,ub,vb,ut,vt,pvsp)
  
  USE spectral

  IMPLICIT NONE

  integer, intent(in) :: nlevs
  real, intent(in), dimension(2*kmax,2*lmax) :: itbxy,ittxy
  real, intent(in), dimension(2*kmax,2*lmax,nlevs) :: ipvxy
  complex, intent(out), dimension(2*kmax,2*lmax,nlevs) :: phi,u,v,w,theta
  complex, intent(out), dimension(2*kmax,2*lmax) ::phib,phit
  complex, intent(out), dimension(2*kmax,2*lmax,nlevs) :: pvsp
  complex, intent(out), dimension(2*kmax,2*lmax) :: tbsp,ttsp,ub,ut,vb,vt
  complex, dimension(2*kmax,2*lmax) :: phi0b,phi0t
  complex, dimension(2*kmax,2*lmax,nlevs) :: F1,G1,P1,phi1
  complex, dimension(2*kmax,2*lmax) :: phi1b,phi1t
  complex, dimension(2*kmax,2*lmax,nlevs) :: phi0,wcont
  real, dimension(2*kmax,2*lmax,nlevs) :: epv
  real, dimension(2*kmax,2*lmax,nlevs) :: vort,div
  integer :: k
  real :: dxs,dys,dzs,dts ! for rolv's python call to get dimensional values
  
!!!!!!!!!!!!!!!!!!!!!!!!
 ! Start of executable: !
!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! test for echo
  call dx_echo(dxs,dys,dzs,dts)
  if (verbose .gt. 1)   print*,'dxs,dys,dzs,dts',dxs,dys,dzs,dts

  ! max values
  if (verbose .gt. 0)  then
     print*,'max |theta(z=0)| = ',maxval(abs(itbxy))
     print*,'max |theta(z=1)| = ',maxval(abs(ittxy))
     print*,'max |pv| = ',maxval(abs(ipvxy))
  endif

  ! open output files:
  call openf(nlevs)

  ! map into spectral space at the same resolution
  call xy_to_sp(itbxy,tbsp,2*kmax,2*lmax,kmax,lmax)
  call xy_to_sp(ittxy,ttsp,2*kmax,2*lmax,kmax,lmax)

  do k=1,nlevs
     call xy_to_sp(ipvxy(:,:,k),pvsp(:,:,k),2*kmax,2*lmax,kmax,lmax)
  enddo

  if (verbose .gt. 0)  then 
     print*,' ';print*,'inverting for leading order geopotential...'
  endif
  call invert(nlevs,tbsp,ttsp,pvsp,phi0b,phi0t,phi0)

  phi = phi0; phib = phi0b; phit = phi0t
  if (verbose .gt. 1) then 
     print*,' ';print*,'max nondimensional leading-order phi = ',maxval(abs(phi))
  endif

  if (order .eq. 1 .or. order .eq. 2 .and. Ross .ne. 0.) then 
     if (verbose .gt. 0) then 
        print*,' '; print*,'inverting for O(R) fields...'
     endif
     call invert_R(nlevs,tbsp,ttsp,phi0b,phi0t,phi0,F1,G1,P1,phi1,phi1b,phi1t,pvsp)
     if (order .eq. 1) then 
        phi = Ross*phi1; phib = Ross*phi1b; phit = Ross*phi1t
     else
        phi = phi + Ross*phi1
        phib = phib + Ross*phi1b
        phit = phit + Ross*phi1t
     endif
  endif

  if (verbose .gt. 0)   print*,'calculating u and v...'
  call uvwtp(nlevs,tbsp,ttsp,phi0,phi0b,phi0t,F1,G1,P1,phi1,u,ub,ut, &
       &     v,vb,vt,w,theta)
  if (verbose .gt. 1)   print*,'max |u| = ',maxval(abs(u))
  if (verbose .gt. 1)   print*,'max |v| = ',maxval(abs(v))

  ! check w by vertically integrating the continuity equation.
  call w_check(nlevs,u,v,w,wcont)

  ! check Ertel PV
  call epv_pseudoheight(nlevs,u,v,theta,ub,vb,tbsp,ut,vt,ttsp,epv)
  !call epv_pseudoheight_new(nlevs,u,v,theta,tbsp,ttsp,phi0,phi0b,phi0t,epv)
  !open(41,file='pv.dat'); write(41,*) epv; close(41)
  do k = 1,nlevs
!     if (verbose .gt. 1) print*,'pv check at level ',k,': ',maxval(abs(epv(:,:,k)-ipvxy(:,:,k)))
     if (verbose .gt. 1) print*,'pv check at level ',k,': ',maxval(ipvxy(:,:,k))
  enddo

  ! compute vorticity and divergence
  call vortdiv(nlevs,u,v,vort,div)
  print *,'max surface vort & div:',maxval(vort(:,:,2)),maxval(div(:,:,2))
  ! write out fields for checking
  !open(41,file='vort.dat'); write(41,*) vort; close(41)
  !open(41,file='div.dat'); write(41,*) div; close(41)

  !Added by rolv to print out values
  call dx_echo(dxs,dys,dzs,dts)
  if (verbose .gt. 1) then
     print *,'Added by rolv 19 Aug 2004'
     print *,'Choose these values as default for wrf'
     print *,'dx:',dxs
     print *,'dy:',dys
     print *,'dz:',dzs
     print *,'dt:',dts
  endif

end subroutine main

!!!!!!!!!!!!!!!
! SUBPROGRAMS !
!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE invert(nlevs,tbsp,ttsp,pvsp,phi0b,phi0t,phi0)
  USE spectral

  ! Originators: G. J. Hakim,  University of Washington
  !
  ! 14 July 2004; version 0.1

  IMPLICIT NONE

  ! Invert PV for streamfunction; compute spectral derivatives on 
  ! the transform (advection) grid.

  integer, intent(in) :: nlevs
  complex, intent(in), dimension(2*kmax,2*lmax) :: tbsp,ttsp
  complex, intent(in), dimension(2*kmax,2*lmax,nlevs) :: pvsp
  complex, intent(out), dimension(2*kmax,2*lmax) :: phi0b,phi0t
  complex, intent(out), dimension(2*kmax,2*lmax,nlevs) :: phi0
  complex, dimension(2*kmax,2*lmax,nlevs) :: bpvsp,lap     
  complex, dimension(2*kmax,2*lmax,nlevs-1) :: pz      
  complex, dimension(2*kmax,2*lmax) :: pxx,pyy,pzz,errsp      
  real, dimension(2*kmax,2*lmax) :: errxy,tmpxy
  real :: ak,bl,dz
  integer :: j,k,l,kk,ll
  complex, dimension(nlevs) :: psi
  complex,dimension(2*kmax,2*lmax),save:: dx,dy,Id
  ! filler for setup_tanh
  real, dimension(nlevs) :: facz,faczo,zl,zhl

  dz = ZH/real(nlevs)

  ! derivative operators
  call d_setup(dx,dy,Id)

  ! factors for tanh vertical grid
  call setup_tanh(nlevs,facz,faczo,zl,zhl)

  ! initialize arrays to zero:
  phi0b=0.;phi0t=0.;phi0=0.

  bpvsp = pvsp
  do k=1,2*kmax; do l=1,2*lmax

     !     add boundary theta to a copy of spectral pv:
     bpvsp(k,l,1) = bpvsp(k,l,1) + (tbsp(k,l)/(facz(1)*dz))
     bpvsp(k,l,nlevs) = bpvsp(k,l,nlevs) - (ttsp(k,l)/(facz(nlevs)*dz))

     !     get wavenumbers
     call get_waves(k,l,ak,bl)

     !     spectral inversion (zero mean)
     if (k .eq. 1 .and. l .eq. 1) then
        phi0b(k,l) = 0.; phi0t(k,l) = 0.; phi0(k,l,:) = 0.
     else
        call matinv(nlevs,bpvsp(k,l,:),ak,bl,psi,1)
        phi0b(k,l) = psi(1) - (0.5*dz*facz(1)*tbsp(k,l))
        phi0t(k,l) = psi(nlevs) + (0.5*dz*facz(nlevs)*ttsp(k,l))
        phi0(k,l,:) = psi(:)
     endif
  enddo; enddo

  ! BEGIN leading-order PV check
  ! phi_xx + phi_yy + phi_zz = q

  call laplacian(nlevs,phi0,tbsp,ttsp,lap)

  ! loop and print our errors
  do k = 1,nlevs

     errsp = pvsp(:,:,k) - lap(:,:,k)
     call sp_to_xy(errsp,errxy,kmax,lmax,2*kmax,2*lmax)
     if (verbose .gt. 1)  print*,'max leading-order solution error at level ',k,' : ',maxval(abs(errxy))

     call sp_to_xy(lap,errxy,kmax,lmax,2*kmax,2*lmax)
     if (verbose .gt. 1) print*,'leading-order lap level ',k,' : ',maxval(abs(errxy))
  enddo
  ! END leading-order PV check

  return
end SUBROUTINE invert

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_waves(k,l,ak,bl)

  !     Originator: G. Hakim, University of Washington.

  ! compute x and y wavenumbers for use in spectral calculations.

  USE spectral
  IMPLICIT NONE

  integer, intent(in) :: k,l
  real, intent(out) :: ak, bl

  ak = facx*real(k - 1); bl = facy*real(l - 1)

  !     other spectral quadrants
  if (k .ge. kmax .and. l .le. lmax) then 
     ak = -1.*facx*real(2*kmax - k + 1)
  elseif (l .ge. lmax .and. k .le. kmax) then 
     bl = -1.*facy*real(2*lmax -l + 1)
  elseif (k .ge. kmax .and. l .ge. lmax) then 
     ak = -1.*facx*real(2*kmax - k + 1)
     bl = -1.*facy*real(2*lmax - l + 1)
  endif

  return
end subroutine get_waves

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine d_setup(dx,dy,Id)

  !     Originator: G. Hakim, University of Washington.

  ! Set up matrices for derivatives and integrals.

  USE spectral
  IMPLICIT NONE

  complex, dimension(2*kmax,2*lmax) :: dx,dy,Id
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

  ! Id
  Id = 1.0

  return
end subroutine d_setup

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

  out = 0.

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
SUBROUTINE matinv(nlevs,pv,ak,bl,psi,idn)
  USE spectral

  ! Solve for psi(z) given pv(z) for a given (k,l) in ** uniform ** N^2:

  IMPLICIT NONE

  integer, intent(in) :: nlevs
  complex, intent(in) :: pv(nlevs)
  real, intent(in) :: ak,bl
  integer, intent(in) :: idn
  complex, intent(out) :: psi(nlevs)
  complex, dimension(nlevs) :: f
  real, dimension(nlevs) :: e
  real :: A,B,C,dz
  integer :: j

  psi = 0.

  ! set up coeficients A, B, and C (note they are set for uniform n^2!!!):
  dz = ZH/real(nlevs)
  A = -1./(dz*dz); B = -1.*((ak*ak) + (bl*bl) + (2./(dz*dz)))
  C = -1./(dz*dz)

  ! idn determines inversion BCs (1.=Neumann, -1.=Dirichlet)
  !      idn = 1

  ! first pass:
  e(1) = A / (B - (C*real(idn)))
  f(1) = pv(1) / (B-(C*real(idn)))

  do j = 2, nlevs 
     e(j) = A / (B - (C*e(j-1)))
     f(j) = (pv(j) + (C*f(j-1))) / (B - (C*e(j-1)))
     !        print*,'deep check=',j,f(j),nlevs
  enddo

  !  print*,'f = ',f(nlevs)
  ! second pass
  if(e(nlevs).eq.idn) then
     psi(nlevs) = 0.0
  else
     psi(nlevs) = f(nlevs) / (1. - (real(idn)*e(nlevs)))  

     !        print*,'deep check=',f(nlevs),real(idn),e(nlevs)
  end if

  do j = 1, nlevs-1                      
     psi(nlevs-j) = (e(nlevs-j)*psi(nlevs+1-j)) + f(nlevs-j)
  enddo

  return
end SUBROUTINE matinv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE xy_to_sp(xy,sp,mx,ny,km,lm)
  !     Originator: G. Hakim, NCAR/MMM

  ! Map an (x,y) array onto a _smaller_ spectral array.
  ! Input: xy(mx,ny) --- a grid point array.
  ! Output: sp(2*km,2*lm) --- a spectral array.

  IMPLICIT NONE

  integer, intent(in) :: km,lm,mx,ny
  real, intent(in) :: xy(mx,ny)
  !      complex, intent(in) :: xy(mx,ny)
  complex, intent(out) :: sp(2*km,2*lm)
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
     sp(k,l) = cop(k,l)

     ! +/-kmax <= k <= -1; 0 <= l <= +/-lmax:
     if (k .gt. 1) then
        sp(kk,l) = cop(k2+k,l)
     endif

     ! 0 <= k <= +/-kmax; +/-lmax <= l <= -1:
     if (l .gt. 1) then
        sp(k,ll) = cop(k,l2+l)
     endif

     ! +/-kmax <= k <= -1; +/-lmax <= l <= -1:
     if ((k .gt. 1) .and. (l .gt. 1)) then
        sp(kk,ll) = cop(k2+k,l2+l)
     endif
  enddo; enddo

  sp = sp/real(mx*ny)

  return
end SUBROUTINE xy_to_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE sp_to_xy(sp,xy2,km,lm,mx,ny)
  !     Originator: G. Hakim, NCAR/MMM

  ! Map an (km,lm) spectral array onto a _bigger_ grid point array.
  ! Input: sp(2*km,2*lm) --- a spectral array.
  ! Output: xy(mx,ny) --- a grid point array.

  IMPLICIT NONE

  complex, intent(in), dimension(2*km,2*lm) :: sp
  integer, intent(in) :: km,lm,mx,ny
  real, intent(out), dimension(mx,ny) :: xy2
  complex, dimension(mx,ny) :: xy
  real rtemp,ctemp
  integer k,l,kk,ll,kmp1,lmp1,k2,l2,MM,stat

  xy = 0.; xy2 = 0.

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
SUBROUTINE init(nlevs,icho,tbxy,ttxy,pvxy)
  USE spectral
  !
  ! initialize pv and theta fields
  !
  IMPLICIT NONE

  integer, intent(in) :: nlevs
  real, intent(out), dimension(2*kmax,2*lmax) :: tbxy,ttxy
  real, intent(out), dimension(2*kmax,2*lmax,nlevs) :: pvxy
  integer, intent(in) :: icho
  real, dimension(2*kmax,2*lmax) :: itbxy,ittxy
  real :: dx,dy,xcen,ycen,x,y,linb,lint,pva,amdx,amdy
  real :: dz,C1,C2,ubaro,asx,asy,rr,ang,z,zlev,asz,rr2,xn,yb,yy
  real :: tnotb,tnott,xnot,ynot,RAN1,top,bot,dth,ath
  integer :: i,j,k,ilev,iold,irecl,iseed,itrans,jtrans,it,jt
  integer :: twokmax,twolmax
  logical :: trans
  real, parameter :: amp = -2., asig=0.5


  dx=XL/real((2*kmax)); dy=YL/real((2*lmax)); dz=ZH/real(nlevs)

  ! make IC right here
  if (icho .eq. 0) then 

     xcen = kmax*dx; ycen = lmax*dy

     do i=1,2*kmax; do j=1,2*lmax
        x = (i-1)*dx; y = (j-1)*dy
        xnot = ycen; ynot = ycen + 0.    
        amdx=min( abs(x-xnot),abs(XL+x-xnot),abs(XL-x+xnot) )
        amdy=min( abs(y-ynot),abs(YL+y-ynot),abs(YL-y+ynot) )
        asx = 2.0; asy = 0.5! elliptical vortex
        !         asx = 1.0; asy = 1.0! symmetrical vortex
        rr = (((amdx/asx)**2) + ((amdy/asy)**2))**.5
        !
        ! plane waves
        !            ttxy(i,j) = amp*cos((2.*pi*x/XL) ) ! + (2.*pi*y/YL))
        !            tbxy(i,j) = 0.
        !            ttxy(i,j) = amp*cos((2.*pi*x/XL) ) ! + (2.*pi*y/YL))
        ! gaussian:
        !            tbxy(i,j) = amp*exp(-1.*((rr/asig)**2))
        !            ttxy(i,j) = 0.
        ! second deriv of gaussian:
        ttxy(i,j) = amp*(1.-(rr*rr))*exp(-1.*((rr/asig)**2))
        tbxy(i,j) = -ttxy(i,j)
        ttxy(i,j) = 0.
        !            tbxy(i,j) = 0.
        ! random:
        !            tbxy(i,j) = amp*(ran1(iseed)-0.5)
        !            ttxy(i,j) = amp*(ran1(iseed)-0.5)
        do k = 1, nlevs
           !               pva = 3.0
           pva = 0.0
           Z = (REAL(k)-.5)*DZ 
           !               zlev = 1.0 ; asz = .35
           zlev = .5 ; asz = .1
           xn = xnot !- (z - 1.) ! this term adds vertical tilt
           amdx=min( abs(x-xn),abs(XL+x-xn),abs(XL-x+xn) )
           rr = (((amdx/asx)**2) + ((amdy/asy)**2))**.5 
           rr2 = abs(zlev-z)
           pvxy(i,j,k) = pva*(1.-(rr*rr))* &
                &              exp(-1.*((rr/asig)**2))* &
                &              exp(-1.*((rr2/(asz))**2))
           if (i.eq.kmax.and.j.eq.lmax) then 
              if (verbose .gt. 1) print*,'pv: ',k,pvxy(i,j,k)
           endif
        enddo

     enddo; enddo

  else
     !     read the IC from a matlab file
     !     for now, just use lower boundary theta

     open(9,file=path//'th_init_B.dat', status='old')
     read(9,*) twokmax; read(9,*) twolmax
     if (twokmax .ne. 2*kmax .or. twolmax .ne. 2*lmax) then 
        print*,'INIT ERROR IN DIMENSIONS...',twokmax,twolmax,2*kmax,2*lmax
        stop
     endif
     do i = 1, 2*kmax; do j=1, 2*lmax
        read(9,*) tbxy(i,j)      
     enddo; enddo
     close(9)

     ! ttxy = 0.

     open(9,file=path//'th_init_T.dat', status='old')
     read(9,*) twokmax; read(9,*) twolmax
     if (twokmax .ne. 2*kmax .or. twolmax .ne. 2*lmax) then 
        print*,'INIT ERROR IN DIMENSIONS...',twokmax,twolmax,2*kmax,2*lmax
        stop
     endif
     do i = 1, 2*kmax; do j=1, 2*lmax
        read(9,*) ttxy(i,j)      
     enddo; enddo
     close(9)

     !ttxy = tbxy ! TESTING!
     ! **** In future, put PV file read here
     pvxy = 0.

  endif

  return
end SUBROUTINE init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Numerical Recipes random number generator:
REAL FUNCTION ran1(idum)
  INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
  REAL AM,EPS,RNMX
  PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836, &
       &NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
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
     end do
     iy=iv(1)
  endif
  k=idum/IQ
  idum=IA*(idum-k*IQ)-IR*k
  if (idum.lt.0) idum=idum+IM
  j=1+iy/NDIV
  iy=iv(j)
  iv(j)=idum
  ran1=min(AM*iy,RNMX)
  return
end FUNCTION ran1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!completely new....
SUBROUTINE dump(iwrite,nlevs,phi,theta,u,v,w, &
     &          tb,tt,phib,phit,ub,vb,ut,vt,pv, &
     &          ophi,otheta,ou,ov,ow, &
     &          otb,ott,ophib,ophit,oub,ovb,out,ovt,opv,pheight)

  USE spectral

  IMPLICIT NONE

  ! input variables
  integer, intent(in) :: iwrite,nlevs
  complex, intent(in), dimension(2*kmax,2*lmax,nlevs) :: &
       & phi,u,v,w,theta,pv
  complex, intent (in), dimension(2*kmax,2*lmax) :: &
       & tb,tt,phib,phit,ub,vb,ut,vt
  ! output variables
  real, intent(out), dimension(nlevs) :: pheight
  real, intent(out), dimension(2*kmax,2*lmax,nlevs) :: & 
     & ophi,ou,ov,ow,otheta,opv
  real, intent (out), dimension(2*kmax,2*lmax) :: & 
     & otb,ott,ophib,ophit,oub,ovb,out,ovt
  integer :: k
  real :: dz
  logical :: go
  real :: grav,tnot,Ns,km,Cor,Ls,Hs,Us,Ws,Ps,Ts,PVs

  call scaling(grav,tnot,Ns,km,Cor,Ls,Hs,Us,Ws,Ps,Ts)
  call PVscaling(PVs)

  dz = ZH/real(nlevs)
  if (verbose .gt. 1) then 
     print*,'Scaling parameters...'
     print*,'grav = ',grav,' m s^-2'
     print*,'Cor = ',Cor,' s^-1'
     print*,'Ls = ',Ls, ' m'
     print*,'Hs = ',Hs, ' m'
     print*,'Us = ',Us, ' m s^-1'
     print*,'Ws = ',Ws*100, ' cm s^-1'
     print*,'Ps = ',Ps, ' m^2 s^-2'
     print*,'Ts = ',Ts, ' K'
     print*,'Ns = ',Ns, ' s^-1'
     print*,'PVs = ',PVs, ' PVU'

     print*,'nondim t(h) =',maxval(abs(tt))
     print*,'nondim pv =',maxval(abs(pv))
  endif

  ! keeps syntex clean below & lets python send an integer
  go = .FALSE.
  if (iwrite .eq. 1) go = .TRUE.

  ! lower boundary theta
  call sp_to_xy(tb,otb,kmax,lmax,2*kmax,2*lmax); otb = Ts*otb
  if (go) write(18,*) otb
  if (verbose .gt. 1) print*,'max z=0 dim |t| = ',maxval(abs(otb))

  ! upper boundary theta
  call sp_to_xy(tt,ott,kmax,lmax,2*kmax,2*lmax); ott = Ts*ott
  if (go) write(19,*) ott
  if (verbose .gt. 1) print*,'max z=H dim |t| = ',maxval(abs(ott))

  ! lower boundary geopotential
  call sp_to_xy(phib,ophib,kmax,lmax,2*kmax,2*lmax); ophib = Ps*ophib
  if (go) write(12,*) ophib
  if (verbose .gt. 1) print*,'max z=0 dim |p| = ',maxval(abs(ophib))

  ! upper boundary geopotential
  call sp_to_xy(phit,ophit,kmax,lmax,2*kmax,2*lmax); ophit = Ps*ophit
  if (go) write(13,*) ophit
    if (verbose .gt. 1) print*,'max z=H dim |p| = ',maxval(abs(ophit))
  
  ! lower boundary u & v
  call sp_to_xy(ub,oub,kmax,lmax,2*kmax,2*lmax); oub = Us*oub
  if (go) write(22,*) oub
  if (verbose .gt. 1) print*,'max z=0 dim |u| = ',maxval(abs(oub))
  call sp_to_xy(vb,ovb,kmax,lmax,2*kmax,2*lmax); ovb = Us*ovb
  if (go) write(23,*) ovb
  if (verbose .gt. 1) print*,'max z=0 dim |v| = ',maxval(abs(ovb))
  
  ! upper boundary u & v
  call sp_to_xy(ut,out,kmax,lmax,2*kmax,2*lmax); out = Us*out
  if (go) write(24,*) out
  if (verbose .gt. 1) print*,'max z=H dim |u| = ',maxval(abs(out))
  call sp_to_xy(vt,ovt,kmax,lmax,2*kmax,2*lmax); ovt = Us*ovt
  if (go) write(25,*) ovt
  if (verbose .gt. 1) print*,'max z=H dim |v| = ',maxval(abs(ovt))

  ! interior geopotential, u, v, and w
  do k = 1, nlevs
     if (verbose .gt. 1) print*,'level ',k

     call sp_to_xy(phi(:,:,k),ophi(:,:,k),kmax,lmax,2*kmax,2*lmax) 
     ophi(:,:,k) = Ps*ophi(:,:,k)
     if (go) then
        write(14,*) ophi(:,:,k);
        if (verbose .gt. 1) print*,'max dim |p| = ',maxval(abs(ophi(:,:,k)))
     endif

     call sp_to_xy(u(:,:,k),ou(:,:,k),kmax,lmax,2*kmax,2*lmax)
     ou(:,:,k) = Us*ou(:,:,k)
     if (go) then
        write(15,*) ou(:,:,k)
        if (verbose .gt. 1) print*,'max dim |u| = ',maxval(abs(ou(:,:,k)))
     endif

     call sp_to_xy(v(:,:,k),ov(:,:,k),kmax,lmax,2*kmax,2*lmax)
     ov(:,:,k) = Us*ov(:,:,k)
     if (go) then
        write(16,*) ov(:,:,k)
        if (verbose .gt. 1) print*,'max dim |v| = ',maxval(abs(ov(:,:,k)))
     endif

     call sp_to_xy(w(:,:,k),ow(:,:,k),kmax,lmax,2*kmax,2*lmax)
     ow(:,:,k) = Ws*ow(:,:,k)
     if (go) then
        write(17,*) ow(:,:,k)
        if (verbose .gt. 1) print*,'max dim |w| = ',maxval(abs(ow(:,:,k)))
     endif
     
     call sp_to_xy(theta(:,:,k),otheta(:,:,k),kmax,lmax,2*kmax,2*lmax)
     otheta(:,:,k) = Ts*otheta(:,:,k)
     if (go) then
        write(20,*) otheta(:,:,k) 
        if (verbose .gt. 1) print*,'max dim |t| = ',maxval(abs(otheta(:,:,k)))
     endif

     call sp_to_xy(pv(:,:,k),opv(:,:,k),kmax,lmax,2*kmax,2*lmax); 
     opv(:,:,k) = PVs*opv(:,:,k)
     if (go) then
        write(21,*) opv(:,:,k); 
        if (verbose .gt. 1) print*,'max dim |pv| = ',maxval(abs(opv(:,:,k)))
     endif
     pheight(k) = Hs*(real(k) - 0.5)*dz
     
  enddo
 
  if (go) then
     close(12); close(13); close(14); close(15); close(16); close(17)
     close(18); close(19)
     close(20); close(21); close(22); close(23); close(24); close(25)
  endif
  
  ! NEW---write dimensional pseudoheight levels (in meters):
  if (go) then
     open(12,file=path//'data/pvinv_z.dat') 
     write(12,*) 2*kmax; write(12,*) 2*lmax; write(12,*) nlevs
     write(12,*) XL; write(12,*) YL; write(12,*) ZH; write(12,*) Ross
     do k = 1, nlevs
        write(12,*) Hs*(real(k) - 0.5)*dz
     enddo
     close(12)
  endif

  return
end SUBROUTINE dump

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE dump_old(nlevs,tb,tt,sb,st,s,pv,ut,vt)
  USE spectral

  IMPLICIT NONE

  integer, intent(in) :: nlevs
  real, dimension(2*kmax,2*lmax), intent(in) :: ut,vt
  complex, intent (in), dimension(2*kmax,2*lmax) :: tb,tt,sb,st
  complex, dimension(2*kmax,2*lmax,nlevs) :: s,pv
  real, dimension(2*kmax,2*lmax,nlevs+1) :: ost,oth,opv,stemp,ptemp
  complex, dimension(2*kmax,2*lmax) :: copy
  real :: dz,sfmin,y,z,dy,pvmax
  integer, save :: irec
  integer :: i,j,k
  logical :: ilin
  character(len=18) :: form
  character(len=10) :: size

  oth = 0.; ost = 0.; opv = 0.

  dz = ZH/real(nlevs); dy = YL/real(2.*lmax)

  ! boundary streamfunction
  copy = sb
  call sp_to_xy(copy,ost(:,:,1),kmax,lmax,2*kmax,2*lmax)
  copy = st
  call sp_to_xy(copy,ost(:,:,nlevs+1),kmax,lmax,2*kmax,2*lmax)

  ! boundary theta
  copy = tb
  call sp_to_xy(copy,oth(:,:,1),kmax,lmax,2*kmax,2*lmax)
  copy = tt
  call sp_to_xy(copy,oth(:,:,nlevs+1),kmax,lmax,2*kmax,2*lmax)

  ! boundary pv
  opv(:,:,1) = amiss
  opv(:,:,nlevs+1) = amiss

  ! interior fields
  do j = 1, nlevs
     copy(:,:) = s(:,:,j)
     call sp_to_xy(copy,stemp(:,:,j),kmax,lmax,2*kmax,2*lmax)
     copy(:,:) = pv(:,:,j)
     call sp_to_xy(copy,ptemp(:,:,j),kmax,lmax,2*kmax,2*lmax)
  enddo

  do i = 1, 2*kmax; do j = 1, 2*lmax; do k = 1, nlevs+1 
     y = (j-1)*dy; z = (k-1)*dz
     if (k .ge. 2 .and. k .le. nlevs) then
        oth(i,j,k) = oth(i,j,k) + ((stemp(i,j,k)-stemp(i,j,k-1))/dz) 
        ost(i,j,k) = (stemp(i,j,k) + stemp(i,j,k-1)) / 2.
        opv(i,j,k) = (ptemp(i,j,k) + ptemp(i,j,k-1)) / 2.
     endif
  enddo; enddo; enddo

   if (verbose .gt. 1) print*,'max |phi(z=0)| = ',maxval(abs(ost(:,:,1)))
   if (verbose .gt. 1) print*,'max |phi(z=1)| = ',maxval(abs(ost(:,:,nlevs+1)))
   if (verbose .gt. 1) print*,'max |phi(0 < z < 1)| = ',maxval(abs(ost(:,:,2:nlevs)))

  !  write to matlab
  do i=1,2*kmax
     do j=1,2*lmax
        write(12,*) oth(i,j,1)
        write(13,*) oth(i,j,nlevs+1)
        write(14,*) ost(i,j,1)
        write(15,*) ost(i,j,nlevs+1)
        write(16,*) ut(i,j)
        write(17,*) vt(i,j)
     enddo
  enddo

  return
end SUBROUTINE dump_old

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE openf(nlevs)

  USE spectral

  IMPLICIT NONE

  integer, intent(in) :: nlevs
  integer ::irecl

  irecl = 2*kmax*2*lmax*(nlevs+1)*32

  ! matlab interface files:
  open(12,file=path//'data/pvinv_pb.dat') 
  write(12,*) 2*kmax; write(12,*) 2*lmax; write(12,*) nlevs
  write(12,*) XL; write(12,*) YL; write(12,*) ZH; write(12,*) Ross

  open(13,file=path//'data/pvinv_pt.dat') 
  write(13,*) 2*kmax; write(13,*) 2*lmax; write(13,*) nlevs
  write(13,*) XL; write(13,*) YL; write(13,*) ZH; write(13,*) Ross

  open(14,file=path//'data/pvinv_phi.dat') 
  write(14,*) 2*kmax; write(14,*) 2*lmax; write(14,*) nlevs
  write(14,*) XL; write(14,*) YL; write(14,*) ZH; write(14,*) Ross

  open(15,file=path//'data/pvinv_u.dat') 
  write(15,*) 2*kmax; write(15,*) 2*lmax; write(15,*) nlevs
  write(15,*) XL; write(15,*) YL; write(15,*) ZH; write(15,*) Ross

  open(16,file=path//'data/pvinv_v.dat') 
  write(16,*) 2*kmax; write(16,*) 2*lmax; write(16,*) nlevs
  write(16,*) XL; write(16,*) YL; write(16,*) ZH; write(16,*) Ross

  open(17,file=path//'data/pvinv_w.dat') 
  write(17,*) 2*kmax; write(17,*) 2*lmax; write(17,*) nlevs
  write(17,*) XL; write(17,*) YL; write(17,*) ZH; write(17,*) Ross

  open(18,file=path//'data/pvinv_tb.dat') 
  write(18,*) 2*kmax; write(18,*) 2*lmax; write(18,*) nlevs
  write(18,*) XL; write(18,*) YL; write(18,*) ZH; write(18,*) Ross

  open(19,file=path//'data/pvinv_tt.dat') 
  write(19,*) 2*kmax; write(19,*) 2*lmax; write(19,*) nlevs
  write(19,*) XL; write(19,*) YL; write(19,*) ZH; write(19,*) Ross

  open(20,file=path//'data/pvinv_t.dat') 
  write(20,*) 2*kmax; write(20,*) 2*lmax; write(20,*) nlevs
  write(20,*) XL; write(20,*) YL; write(20,*) ZH; write(20,*) Ross

  open(21,file=path//'data/pvinv_pv.dat') 
  write(21,*) 2*kmax; write(21,*) 2*lmax; write(21,*) nlevs
  write(21,*) XL; write(21,*) YL; write(21,*) ZH; write(21,*) Ross

  !these are new (09/05)
  open(22,file=path//'data/pvinv_ub.dat') 
  write(22,*) 2*kmax; write(22,*) 2*lmax; write(22,*) nlevs
  write(22,*) XL; write(22,*) YL; write(22,*) ZH; write(22,*) Ross

  open(23,file=path//'data/pvinv_vb.dat') 
  write(23,*) 2*kmax; write(23,*) 2*lmax; write(23,*) nlevs
  write(23,*) XL; write(23,*) YL; write(23,*) ZH; write(23,*) Ross

  open(24,file=path//'data/pvinv_ut.dat') 
  write(24,*) 2*kmax; write(24,*) 2*lmax; write(24,*) nlevs
  write(24,*) XL; write(24,*) YL; write(24,*) ZH; write(24,*) Ross

  open(25,file=path//'data/pvinv_vt.dat') 
  write(25,*) 2*kmax; write(25,*) 2*lmax; write(25,*) nlevs
  write(25,*) XL; write(25,*) YL; write(25,*) ZH; write(25,*) Ross

  return
end SUBROUTINE openf

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
  character(len=1) :: csign
  !     integer*4 :: stat
  integer :: stat

  if (dec) then 
     if (isign .eq. -1) then 
        csign = 'F'
     elseif (isign .eq. 1) then 
        csign = 'B'
        f = f*ni*nj
     else
        print*,'ERROR IN FT_2D CALl'
        stop
     endif
     !         stat = CFFT_2D('C','C',csign,f,f,ni,nj,ni,1,1)
     if (stat .ne. 0) then 
        print*,'ERROR IN DEC FFT CALL:',stat
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
!
! this is the old inversion routine. the new one has the quadrant work 
! hidden in operators, so it scales better for +1 code.
!
SUBROUTINE invert_old(nlevs,tbsp,ttsp,pvsp,phi0b,phi0t,phi0)

  USE spectral

  IMPLICIT NONE

  ! Invert PV for streamfunction; compute spectral derivatives on 
  ! the transform (advection) grid.

  integer, intent(in) :: nlevs
  complex, intent(in), dimension(2*kmax,2*lmax) :: tbsp,ttsp
  complex, intent(in), dimension(2*kmax,2*lmax,nlevs) :: pvsp
  complex, intent(out), dimension(2*kmax,2*lmax) :: phi0b,phi0t
  complex, intent(out), dimension(2*kmax,2*lmax,nlevs) :: phi0
  complex, dimension(2*kmax,2*lmax,nlevs) :: bpvsp      
  real :: ak,bl,kap,sqg,amod,dz
  integer :: j,k,l,kk,ll
  complex, dimension(nlevs) :: psi

  dz = ZH/real(nlevs)

  ! initialize arrays to zero:
  phi0b=0.;phi0t=0.;phi0=0.

  ! add boundary theta to a copy of spectral pv:
  bpvsp = pvsp
  do k=1,2*kmax; do l=1,2*lmax
     bpvsp(k,l,1) = bpvsp(k,l,1) + (tbsp(k,l)/dz)
     bpvsp(k,l,nlevs) = bpvsp(k,l,nlevs) - (ttsp(k,l)/dz)
  enddo; enddo

  ! inversion of (kmax,lmax) waves into (mmax,nmax) arrays:
  do k = 1,kmaxp1
     do l = 1,lmaxp1

        ! waves: 0 <= k <= +/-kmax; 0 <= l <= +/-lmax:
        ak = facx*real(k - 1); bl = facy*real(l - 1)

        if (k .eq. 1 .and. l .eq. 1) then
           phi0(k,l,:) = 0.; phi0(k,l,:) = 0.
        else
           call matinv(nlevs,bpvsp(k,l,:),ak,bl,psi,1)
           phi0b(k,l) = psi(1) - (0.5*dz*tbsp(k,l))
           phi0t(k,l) = psi(nlevs) + (0.5*dz*ttsp(k,l))
           phi0(k,l,:) = psi
        endif

        kk = kmax + k; ll = lmax + l

        ! +/-kmax <= k <= -1; 0 <= l <= +/-lmax:
        if (k .le. kmax) then
           ak = -1.*facx*real(kmaxp1 - k); bl = facy*real(l - 1)
           call matinv(nlevs,bpvsp(kk,l,:),ak,bl,psi,1)
           phi0b(kk,l) = psi(1) - (0.5*dz*tbsp(kk,l))
           phi0t(kk,l) = psi(nlevs) + (0.5*dz*ttsp(kk,l))
           phi0(kk,l,:) = psi
        endif

        ! 0 <= k <= +/-kmax; +/-lmax <= l <= -1:
        if (l .le. lmax) then
           ak = facx*real(k-1); bl = -1.*facy*real(lmaxp1 - l)
           call matinv(nlevs,bpvsp(k,ll,:),ak,bl,psi,1)
           phi0b(k,ll) = psi(1) - (0.5*dz*tbsp(k,ll))
           phi0t(k,ll) = psi(nlevs) + (0.5*dz*ttsp(k,ll))
           phi0(k,ll,:) = psi
        endif

        ! +/-kmax <= k <= -1; +/-lmax <= l <= -1:
        if (k .le. kmax .and. l .le. lmax) then
           ak=-1.*facx*real(kmaxp1 - k); bl=-1.*facy*real(lmaxp1 - l)
           call matinv(nlevs,bpvsp(kk,ll,:),ak,bl,psi,1)
           phi0b(kk,ll) = psi(1) - (0.5*dz*tbsp(kk,ll))
           phi0t(kk,ll) = psi(nlevs) + (0.5*dz*ttsp(kk,ll))
           phi0(kk,ll,:) = psi
        endif
     end do
  end do

  return
end SUBROUTINE invert_old

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! recover winds, theta, and geopotential from potential data.
!
subroutine uvwtp(nlevs,tb,tt,phi0,phi0b,phi0t,F1,G1,P1,phi1, &
     &                 u,ub,ut,v,vb,vt,w,thF)

  USE spectral

  IMPLICIT NONE

  integer, intent(in) :: nlevs
  complex, dimension(2*kmax,2*lmax,nlevs), intent(in) :: phi0, &
       &                   F1,G1,P1,phi1
  complex, dimension(2*kmax,2*lmax), intent(in) :: phi0b,phi0t,tb,tt
  complex, dimension(2*kmax,2*lmax,nlevs), intent(out) :: u,v,w,thF
  complex, dimension(2*kmax,2*lmax), intent(out) :: ub,ut,vb,vt
  complex, dimension(2*kmax,2*lmax) :: dx,dy,Id,uF,vF,temp,ftemp,oftemp
  complex, dimension(2*kmax,2*lmax,nlevs) :: u1,v1
  complex, dimension(2*kmax,2*lmax,nlevs) :: phi
  complex, dimension(2*kmax,2*lmax) :: ug,vg
  real :: dz
  integer :: k,l,p
  complex :: gtemp,ogtemp
  ! testing...
  real, dimension(2*kmax,2*lmax) :: rtemp,rutemp,rvtemp,ws,wsg

  dz = ZH/real(nlevs)

  ! derivative operators
  call d_setup(dx,dy,Id)

  ! 3D leading-order winds
  if (order .eq. 0 .or. order .eq. 2) then 

     do k = 1,nlevs
        call d_s2s(phi0(:,:,k),uF,1,-dy); u(:,:,k) = uF
        call d_s2s(phi0(:,:,k),vF,1,dx); v(:,:,k) = vF
     enddo

     ! leading-order boundary winds
     call d_s2s(phi0b,uF,1,-dy); ub = uF
     call d_s2s(phi0t,uF,1,-dy); ut = uF

     call d_s2s(phi0b,vF,1,dx); vb = vF
     call d_s2s(phi0t,vF,1,dx); vt = vF

  else
     u = 0.; v = 0.
  endif ! order 0

  if (order .eq. 1 .or. order .eq. 2 .and. Ross .ne. 0.) then 

  if (verbose .gt. 0) print*,'computing next-order corrections to u and v...'

     ! next-order corrections
     ! v1 =  phi1_x - G1_z
     ! u1 = -phi1_y - F1_z
     ! w  =  F1_x + G1_y

     ! compute vertical derivatives by finite differences
     do p=1,nlevs
        if (p .eq. 1) then ! linear interp; zero BC
           u1(:,:,1) = -(F1(:,:,1) + F1(:,:,2))/(2*dz)
           v1(:,:,1) = -(G1(:,:,1) + G1(:,:,2))/(2*dz)
           ub(:,:) = ub(:,:) - Ross*(2*F1(:,:,1)/dz)
           vb(:,:) = vb(:,:) - Ross*(2*G1(:,:,1)/dz)  
        elseif (p .eq. nlevs) then ! linear interp; zero BC
           u1(:,:,nlevs) = (F1(:,:,nlevs) + F1(:,:,nlevs-1))/(2*dz)
           v1(:,:,nlevs) = (G1(:,:,nlevs) + G1(:,:,nlevs-1))/(2*dz)
           ut(:,:) = ut(:,:) + Ross*(2*F1(:,:,nlevs)/dz)
           vt(:,:) = vt(:,:) + Ross*(2*G1(:,:,nlevs)/dz)
        else
           u1(:,:,p) = -(F1(:,:,p+1) - F1(:,:,p-1))/(2*dz)
           v1(:,:,p) = -(G1(:,:,p+1) - G1(:,:,p-1))/(2*dz)
        endif
     enddo ! p 

     !      u1 = 0; v1 = 0 ! test of phi contrib
     ! compute horizonal derivatives spectrally
     do p = 1, nlevs
        call d_s2s(P1(:,:,p),temp,1,dx); v1(:,:,p) = v1(:,:,p) + temp
        call d_s2s(P1(:,:,p),temp,1,-dy); u1(:,:,p) = u1(:,:,p) + temp
        if (p .eq. 1) then ! recall phi1_z = 0
           call d_s2s(P1(:,:,p),temp,1,dx); vb = vb + Ross*temp
           call d_s2s(P1(:,:,p),temp,1,-dy); ub = ub + Ross*temp
        endif
        if (p .eq. nlevs) then ! recall phi1_z = 0
           call d_s2s(P1(:,:,p),temp,1,dx); vt = vt + Ross*temp
           call d_s2s(P1(:,:,p),temp,1,-dy); ut = ut + Ross*temp
        endif
        call d_s2s(F1(:,:,p),w(:,:,p),1,dx)
        call d_s2s(G1(:,:,p),temp,1,dy)
        w(:,:,p) = Ross*(w(:,:,p) + temp)
     enddo
     
  else
     u1 = 0.; v1 = 0.; w = 0.
  endif ! order 1

  u = u + Ross*u1
  v = v + Ross*v1

  ! OLD
  ! compute theta and get it on the regular, staggered vertical grid
  !if (order .eq. 0) then 
  !   phi = phi0
  !else
  !   phi = phi0 + (Ross*phi1) ! total geopotential
  !endif
  !thF = 0.
  ! compute vertical derivatives by finite differences
  ! on staggered grid starting dz off boundary
  !do p=1,nlevs-1
  !   thF(:,:,p) = (phi(:,:,p+1) - phi(:,:,p))/dz
  !enddo ! p 

  ! NEW---try computing theta from potential functions
  ! theta = Phi_z + G_x - F_y

  if (1 .eq. 1) then
  thF = 0.
  if (order .eq. 0) then 
     phi = phi0             ! total phi potential -- different from above!
  else
     phi = phi0 + (Ross*P1) ! total phi potential -- different from above!     
  endif
  do p=1,nlevs
     !P1_z:
     if (p .eq. 1) then 
        thF(:,:,p) = (phi(:,:,p+1) - (phi(:,:,p)-tb*dz))/(2.*dz)
     elseif (p .eq. nlevs) then 
        thF(:,:,p) = ((phi(:,:,p)+tt*dz) - phi(:,:,p-1))/(2.*dz)
     else
        thF(:,:,p) = (phi(:,:,p+1) - phi(:,:,p-1))/(2.*dz)
     endif
     if (order .ne. 0) then 
        !G_x:
        call d_s2s(G1(:,:,p),temp,1,dx); thF(:,:,p) = thF(:,:,p) + Ross*temp
        !F_y
        call d_s2s(F1(:,:,p),temp,1,dy); thF(:,:,p) = thF(:,:,p) - Ross*temp
     endif
  enddo ! p 
  endif

! NEW NEW NEW--theta at half-levels
  if (1 .eq. 0) then 
  thF = 0.
  phi = phi0   
  do p=1,nlevs-1
     thF(:,:,p) = (phi(:,:,p+1) - phi(:,:,p))/(dz)
     thF(:,:,p) = thF(:,:,p) + Ross*((P1(:,:,p+1) - P1(:,:,p))/(dz))
     if (order .ne. 0) then 
        !G_x:
        call d_s2s(G1(:,:,p),temp,1,dx); call d_s2s(G1(:,:,p+1),ftemp,1,dx) 
        thF(:,:,p) = thF(:,:,p) + Ross*(temp + ftemp)/2.
        !F_y
        call d_s2s(F1(:,:,p),temp,1,dy); call d_s2s(F1(:,:,p+1),ftemp,1,dy); 
        thF(:,:,p) = thF(:,:,p) - Ross*(temp + ftemp)/2.
     endif
  enddo ! p 
  endif

  ! compute the geostrophic wind from the _total_ geopotential 
  do p = 1, nlevs
     rutemp = 0.; rvtemp = 0.
     call d_s2s(phi(:,:,p),vg,1,dx);  call d_s2s(phi(:,:,p),ug,1,-dy);
     ! move full wind onto physical grid
     call d_s2s(vg,rvtemp,1,Id); call d_s2s(ug,rutemp,1,Id)
     wsg = sqrt(((rutemp**2) + (rvtemp**2)))
     ! move full wind onto physical grid
     rutemp = 0.; rvtemp = 0.
     call d_s2s(v(:,:,p),rvtemp,1,Id); call d_s2s(u(:,:,p),rutemp,1,Id)
     ws = sqrt(((rutemp**2) + (rvtemp**2)))
     if (verbose .gt. 1) print*,'max wind speed for full and geostrophic winds (total phi): ', &
          & maxval(ws),maxval(wsg)
  enddo

  return
end subroutine uvwtp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine invert_R(nlevs,tbsp,ttsp,phi0b,phi0t,phi0,F1,G1,P1,phi1,phi1b,phi1t,pvsp)

  USE spectral

  ! Originator: G. J. Hakim,  University of Washington

  IMPLICIT NONE

  ! Invert for O(R) potentials.
  !
  ! phi0 : leading order geopotential.
  ! phi1 : O(Ross) geopotential correction.
  ! phi1b, phi1t : O(Ross) boundary geopotential corrections.
  ! F1, G1, P1 : 3D O(Ross) potentials.
  ! pz, pzz, pxy, etc. : derivatives of phi0.
  ! bpvsp : boundary condition trick (make homogeneous BCs; correct int).
  !
  integer, intent(in) :: nlevs
  complex, intent(in), dimension(2*kmax,2*lmax) :: tbsp,ttsp,phi0b,phi0t
  complex, intent(in), dimension(2*kmax,2*lmax,nlevs) :: phi0,pvsp
  !
  complex, intent(out), dimension(2*kmax,2*lmax,nlevs) :: F1,G1,P1,phi1
  complex, intent(out), dimension(2*kmax,2*lmax) :: phi1b,phi1t
  !
  complex, dimension(2*kmax,2*lmax,nlevs) :: F1s,G1s,P1s,phi1s,lap
  real :: ak,bl,dz,solv
  integer :: j,k,l,kk,ll,p
  complex,dimension(2*kmax,2*lmax),save:: dx,dy,Id
  complex, dimension(2*kmax,2*lmax,nlevs) :: pzS,pzzS
  complex, dimension(2*kmax,2*lmax) :: bcb,bct,errsp,tb,tt
  real, dimension(mmax,nmax) :: pzx,pzy,pzz,pxx,pyy,pxy,Rtemp,pvxy,tbx,tby,pzxold,pzyold

  complex, dimension(mmax,nmax) :: Ctemp
  complex, dimension(nlevs) :: psi
  real, dimension(2*kmax,2*lmax) :: Ptemp,errxy,tmpxy
  complex, dimension(2*kmax,2*lmax) :: solvB,solvT
  ! filler for setup_tanh
  real, dimension(nlevs) :: facz,faczo,zl,zhl

  dz = ZH/real(nlevs)

  ! derivative operators
  call d_setup(dx,dy,Id)

  ! factors for tanh vertical grid
  call setup_tanh(nlevs,facz,faczo,zl,zhl)
  
  tb = tbsp; tt = ttsp;
  tb(1,1) = 0.; tt(1,1) = 0.

  ! solvability condition [\int q dV = \int \theta^1 + \phi0_z(\phi0_xx + \phi0_yy)dxdy]_0^1
  ! \theta^1 = 0 ==> q = C/H where C = avg[\phi0_z(\phi0_xx + \phi0_yy)].
  ! bottom
  call d_s2b(phi0b,Ctemp,1,dx*dx); call ft_2d(Ctemp,mmax,nmax,1); pxx = real(Ctemp)
  call d_s2b(phi0b,Ctemp,1,dy*dy); call ft_2d(Ctemp,mmax,nmax,1); pyy = real(Ctemp)
  ! this is an estimate of pzz 1/8 dz off the surface; pre-use solvB
  solvB = (((phi0(:,:,1)-phi0b)*2./(facz(1)*dz)) - tb)*4./(facz(1)*dz) 
  call d_s2b(solvB,Ctemp,1,Id); call ft_2d(Ctemp,mmax,nmax,1); pzz = real(Ctemp)
  Rtemp = pzz*(pxx + pyy); call xy_to_sp(Rtemp,solvB,mmax,nmax,kmax,lmax)
  ! top
  call d_s2b(phi0t,Ctemp,1,dx*dx); call ft_2d(Ctemp,mmax,nmax,1); pxx = real(Ctemp)
  call d_s2b(phi0t,Ctemp,1,dy*dy); call ft_2d(Ctemp,mmax,nmax,1); pyy = real(Ctemp)
  ! this is an estimate of pzz 1/8 dz off the surface; pre-use solvT
  solvT = (tt - ((phi0t - phi0(:,:,nlevs))*2./(facz(nlevs)*dz)))*4./(facz(nlevs)*dz)
  call d_s2b(solvT,Ctemp,1,Id); call ft_2d(Ctemp,mmax,nmax,1); pzz = real(Ctemp)
  Rtemp = pzz*(pxx + pyy); call xy_to_sp(Rtemp,solvT,mmax,nmax,kmax,lmax)
  solv = (solvT(1,1) - solvB(1,1))/H
  if (verbose .gt. 1)  print*,'*-*-*-* O(R) solvability constant : ',solv

  ! phi_z at half levels; phi_zz on grid levels
  do k = 1,nlevs
     if (k .le. nlevs-1) then
        pzS(:,:,k) = (phi0(:,:,k+1) - phi0(:,:,k)) / (faczo(k)*dz)
     endif
     if (k .eq. 1) then
        pzzS(:,:,k) = (pzS(:,:,k) - tb) / (facz(k)*dz)
     elseif (k .eq. nlevs) then
        pzzS(:,:,k) = (tt - pzS(:,:,k-1)) / (facz(k)*dz)
     else
        pzzS(:,:,k) = (pzS(:,:,k) - pzS(:,:,k-1)) / (facz(k)*dz)
     endif
  enddo
  pzzS(1,1,:) = 0.
  
  ! terms with mixed derivatives involving d/dz at half levels on big grid
  do k=1,nlevs

     pzz=0.;pzx=0.;pzy=0.;pxx=0.;pyy=0.;pxy=0.

     if (k .gt. 2 .and. k .lt. nlevs) then
        pzxold = pzx; pzyold = pzy
     endif
     if (k .le. nlevs-1) then
        call d_s2b(pzS(:,:,k),Ctemp,1,dx)
        call ft_2d(Ctemp,mmax,nmax,1); pzx = real(Ctemp)
        call d_s2b(pzS(:,:,k),Ctemp,1,dy)
        call ft_2d(Ctemp,mmax,nmax,1); pzy = real(Ctemp)
     endif
     ! start averaging back to grid levels
     if (k .eq. 1) then
        pzxold = pzx; pzyold = pzy
        call d_s2b(tb,Ctemp,1,dx)
        call ft_2d(Ctemp,mmax,nmax,1); tbx = real(Ctemp)
        call d_s2b(tb,Ctemp,1,dy)
        call ft_2d(Ctemp,mmax,nmax,1); tby = real(Ctemp)
        tbx = -tbx; tby = -tby
        pzx = (pzx + tbx)/2.; pzy = (pzy + tby)/2.
     elseif (k .eq. nlevs) then ! recycle tbx & tby names here
        call d_s2b(tt,Ctemp,1,dx)
        call ft_2d(Ctemp,mmax,nmax,1); tbx = real(Ctemp)
        call d_s2b(tt,Ctemp,1,dy)
        call ft_2d(Ctemp,mmax,nmax,1); tby = real(Ctemp)
        pzx = (pzxold + tbx)/2.; pzy = (pzyold + tby)/2.
     else
        if (k .eq. nlevs-1) then 
           pzxold = pzx; pzyold = pzy
        endif
        pzx = (pzx + pzxold)/2.; pzy = (pzy + pzyold)/2.
     endif
     ! end averaging back to grid levels
     
     call d_s2b(pzzS(:,:,k),Ctemp,1,Id)
     call ft_2d(Ctemp,mmax,nmax,1); pzz = real(Ctemp)
     call d_s2b(phi0(:,:,k),Ctemp,1,dx*dx)
     call ft_2d(Ctemp,mmax,nmax,1); pxx = real(Ctemp)
     call d_s2b(phi0(:,:,k),Ctemp,1,dy*dy)
     call ft_2d(Ctemp,mmax,nmax,1); pyy = real(Ctemp)
     call d_s2b(phi0(:,:,k),Ctemp,1,dx*dy)
     call ft_2d(Ctemp,mmax,nmax,1); pxy = real(Ctemp)

     ! compute source terms
     Rtemp = 2*(pzx*pxy - pzy*pxx)
     call xy_to_sp(Rtemp,F1s(:,:,k),mmax,nmax,kmax,lmax)
     Rtemp = 2*(pzx*pyy - pzy*pxy)
     call xy_to_sp(Rtemp,G1s(:,:,k),mmax,nmax,kmax,lmax)
     Rtemp = (pzx*pzx) + (pzy*pzy) - ((pxx + pyy)*pzz)
!     Rtemp = (pzx*pzx) + (pzy*pzy) + (pzz*pzz) ! test for zero pv
     call xy_to_sp(Rtemp,P1s(:,:,k),mmax,nmax,kmax,lmax)
     P1s(1,1,k) = P1s(1,1,k) + solv
     Rtemp = 2*(pxx*pyy - pxy*pxy)
     call xy_to_sp(Rtemp,phi1s(:,:,k),mmax,nmax,kmax,lmax)

  enddo

  ! inversion (mod 19 July to correct F,G calls to matinv)
  do k=1,2*kmax; do l=1,2*lmax

     ! get wavenumbers
     call get_waves(k,l,ak,bl)

     ! spectral inversion (zero mean)
     if (k .eq. 1 .and. l .eq. 1) then
        F1s(k,l,:)=0.; G1s(k,l,:)=0.; P1s(k,l,:)=0.; phi1s(k,l,:)=0.
!        F1s(k,l,:)=0.; G1s(k,l,:)=0.; phi1s(k,l,:)=0.
        F1(k,l,:) = 0.; G1(k,l,:)=0.; P1(k,l,:) = 0.
!     endif
     else
        ! phi (geopotential) 7/20/04 mods here
        ! phi1 BC turns out to be homog. Neumann
        call matinv(nlevs,phi1s(k,l,:),ak,bl,psi,1)
        phi1b(k,l) = psi(1)
        phi1t(k,l) = psi(nlevs)
        phi1(k,l,:) = psi(:)
        ! F
        call matinv(nlevs,F1s(k,l,:),ak,bl,psi,-1); F1(k,l,:) = psi(:)
        ! G
        call matinv(nlevs,G1s(k,l,:),ak,bl,psi,-1); G1(k,l,:) = psi(:)
        ! Phi
        call matinv(nlevs,P1s(k,l,:),ak,bl,psi,1); P1(k,l,:) = psi(:)
     endif
  enddo; enddo

  ! finally, correct phi from Phi
  phi1 = phi1 + P1
  phi1b = phi1b + P1(:,:,1) ! P1_z = 0--> P1(z=0) = P1(z=dz/2)
  phi1t = phi1t + P1(:,:,nlevs)

  !---------- Phi1 solution check
  bcb = 0.; bct = 0.
  call laplacian(nlevs,P1,bcb,bct,lap)

  do k = 1,nlevs
     errsp = P1s(:,:,k) - lap(:,:,k); errsp(1,1) = 0.
     call sp_to_xy(errsp,errxy,kmax,lmax,2*kmax,2*lmax)
     if (verbose .gt. 1) print*,'max P1 solution error at level ',k,' : ',maxval(abs(errxy))

  enddo
  !---------- F1 solution check
  bcb = 2*F1(:,:,1)/dz; bct = -2*F1(:,:,nlevs)/dz
  call laplacian(nlevs,F1,bcb,bct,lap)
  do k = 1,nlevs
     errsp = F1s(:,:,k) - lap(:,:,k)
     call sp_to_xy(errsp,errxy,kmax,lmax,2*kmax,2*lmax)
     if (verbose .gt. 1) print*,'max F1 solution error at level ',k,' : ',maxval(abs(errxy))
  enddo

  !---------- G1 solution check
  bcb = 2*G1(:,:,1)/dz; bct = -2*G1(:,:,nlevs)/dz
  call laplacian(nlevs,G1,bcb,bct,lap)
  do k = 1,nlevs
     errsp = G1s(:,:,k) - lap(:,:,k)
     call sp_to_xy(errsp,errxy,kmax,lmax,2*kmax,2*lmax)
     if (verbose .gt. 1) print*,'max G1 solution error at level ',k,' : ',maxval(abs(errxy))
  enddo

  return
end subroutine invert_R


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine w_check(nlevs,u,v,win,wout)

  USE spectral

  ! Originator: G. J. Hakim,  University of Washington

  IMPLICIT NONE

  integer, intent(in) :: nlevs
  complex, dimension(2*kmax,2*lmax,nlevs), intent(in) :: u,v,win
  complex, dimension(2*kmax,2*lmax,nlevs), intent(out) :: wout
  complex, dimension(2*kmax,2*lmax) :: dx,dy,Id
  complex, dimension(2*kmax,2*lmax,nlevs) :: w
  complex, dimension(2*kmax,2*lmax) :: temp,wnew,wold,ux,vy,wF
  real, dimension(2*kmax,2*lmax) :: wxo,wxn      
  real :: dz
  integer:: p

  ! vertical grid distance
  dz = ZH/real(nlevs)

  ! derivative operators
  call d_setup(dx,dy,Id)

  ! integrate u_x + v_y = -Ross w_z
  ! note that u and v are next-order, so Ross on RHS cancels w/ those that 
  ! are on the LHS.
  wold = 0.
  do p = 1, nlevs
     call d_s2s(u(:,:,p),ux,1,dx) 
     call d_s2s(v(:,:,p),vy,1,dy) 
     wnew = wold - (ux + vy)*dz
     !         print*,'max ux = ',maxval(abs(ux)),maxval(abs(vy)),dz, &
     !     &        maxval(abs(ux+vy))
     !wF = (wnew + wold)/2. !why???
     wF = wnew
     wout(:,:,p) = wnew
     wold = wnew; wxn = 0.
     call sp_to_xy(win(:,:,p),wxo,kmax,lmax,2*kmax,2*lmax)
     call sp_to_xy(wF,wxn,kmax,lmax,2*kmax,2*lmax)
     if (verbose .gt. 1) print*,'max w comparison (old,new): ',maxval(abs(wxo)), &
          &        maxval(abs(wxn))
  enddo

  return
end subroutine w_check

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine laplacian(nlevs,func,fb,ft,lap)

  USE spectral

  ! compute an isotropic 3D laplacian given a staggered field & BCs.
  ! Originator: G. J. Hakim,  University of Washington

  IMPLICIT NONE

  integer, intent(in) :: nlevs
  complex, dimension(2*kmax,2*lmax,nlevs), intent(in) :: func
  complex, dimension(2*kmax,2*lmax), intent(in) :: fb,ft
  complex, dimension(2*kmax,2*lmax,nlevs), intent(out) :: lap
  integer :: k
  complex, dimension(2*kmax,2*lmax,nlevs-1) :: fz      
  complex, dimension(2*kmax,2*lmax) :: fxx,fyy,fzz      
  complex,dimension(2*kmax,2*lmax) :: dx,dy,Id
  real :: dz

  dz = ZH/real(nlevs)

  ! derivative operators
  call d_setup(dx,dy,Id)

  dz = ZH/real(nlevs)

  ! first z-loop for f_z at at intermediate levels
  do k = 1,nlevs-1
     fz(:,:,k) = (func(:,:,k+1) - func(:,:,k)) / dz
  enddo

  ! second z-loop for laplacian at grid levels and f_zz; lap calc
  do k = 1,nlevs
     
     call d_s2s(func(:,:,k),fxx,1,dx*dx)
     call d_s2s(func(:,:,k),fyy,1,dy*dy)

    if (k .eq. 1) then
        fzz = (fz(:,:,k) - fb) / dz
     elseif (k .eq. nlevs) then
        fzz = (ft - fz(:,:,k-1)) / dz
     else
        fzz = (fz(:,:,k) - fz(:,:,k-1)) / dz
     endif
     ! this is critical to get the inversion to check; the mean is zero in the inversion
     fzz(1,1) = 0.

     lap(:,:,k) = fxx + fyy + fzz

  enddo
! END leading-order PV check

return
end subroutine laplacian

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine vortdiv(nlevs,u,v,vort,div)

  USE spectral

  ! given horizontal wind (u,v) compute vertical component of 
  ! vorticity and horizontal divergence. NONDIMENSIONAL I/O.
  ! make these values dimensional by multiplying by U/L
  ! Originator: G. J. Hakim,  University of Washington

  IMPLICIT NONE

  integer, intent(in) :: nlevs
  complex, dimension(2*kmax,2*lmax,nlevs), intent(in) :: u,v
  real, dimension(2*kmax,2*lmax,nlevs), intent(out) :: vort,div
  complex, dimension(2*kmax,2*lmax) :: ux,uy,vx,vy
  complex,dimension(2*kmax,2*lmax) :: dx,dy,Id,vs,ds
  integer :: k

  ! derivative operators
  call d_setup(dx,dy,Id)

  do k = 1,nlevs
     call d_s2s(v(:,:,k),vx,1,dx); call d_s2s(v(:,:,k),vy,1,dy) 
     call d_s2s(u(:,:,k),ux,1,dx); call d_s2s(u(:,:,k),uy,1,dy) 
     ! spectral vorticity and divergence
     vs = vx - uy
     ds = ux + vy
     ! grid point
     call sp_to_xy(vs,vort(:,:,k),kmax,lmax,2*kmax,2*lmax)
     call sp_to_xy(ds,div(:,:,k),kmax,lmax,2*kmax,2*lmax)
  enddo

  return
end subroutine vortdiv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine setup_tanh(nlevs,facz,faczo,zl,zhl)

  USE spectral

  ! setup Jacobian factors for tanh vertical grid

  IMPLICIT NONE

  integer, intent(in) :: nlevs
  real, dimension(nlevs), intent(out) :: facz, faczo,zl,zhl
  integer :: k
  real :: fac,alpha,acosh,znot,z,arg,dz

  ! if we're using the regular grid, set jacobian values to one:
  if (.not. tanhgrid) then 
     facz = 1.; faczo = 1.
     return
  endif

  dz = ZH/real(nlevs)

  ! resolution ratio factor (near-boundary res is fac * mid level)
  fac = 10; 

  ! nondimensional midlevel in *both* coordinates
  znot = .5; 

  ! grid stretch factor
  alpha = 2.*acosh(sqrt(fac));

  ! first do Jacobian factors for interior grid levels that have the 
  ! first level located one-half grid level off the boundary.
  do k = 1,nlevs     
     z = (real(k) - .5)*dz     
     arg = alpha*(z-.5)
     facz(k) = 2.*tanh(alpha/2.)*(1./(1. - (tanh(arg))**2))/alpha
     zl(k) = znot + (tanh(alpha*(z-.5))) / (2.*tanh(alpha/2.))
     print*,'jacobian:',k,z,zl(k),facz(k)
  enddo

  !now do the "half levels", which are actually spaced evenly with no staggering.
  do k = 1,nlevs-1
     z = real(k)*dz     
     arg = alpha*(z-.5)
     faczo(k) = 2.*tanh(alpha/2.)*(1./(1. - (tanh(arg))**2))/alpha
     zhl(k) = znot + (tanh(alpha*(z-.5))) / (2.*tanh(alpha/2.))
     print*,'jacobian:',k,z,zhl(k),faczo(k)
  enddo

  return
end subroutine setup_tanh

real function acosh(x)

  real, intent(in) :: x
  
  if (x .lt. 1) then 
     print*,'error in acosh routine...bad argument'
     return
  else
     acosh = log(x + sqrt(x*x - 1))
  endif
  return
     
end function acosh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine epv_pseudoheight(nlevs,u,v,t,ub,vb,tb,ut,vt,tt,out)

  USE spectral

  IMPLICIT NONE

  ! source: svn+ssh://modon.atmos.washington.edu/home/disk/modon/hakim/svn/pvinv
  !
  ! Originator: G. J. Hakim (University of Washington); hakim@atmos.washington.edu
  !
  !
  ! inputs: horizontal winds on grid levels: u,v.
  !         horizontal winds on grid boundaries: ub,vb.
  !         potential temperature on grid levels: t
  !         potential temperature on grid boundaries: tb
  ! output: Ertel potential vorticity: out 
  !
  ! epv = v_x - u_y + t_z + Ross * [(v_x - u_y) * t_z - v_z * t_x + u_z * t_y]
  !
  ! all variables are SPECTRAL and NON-DIMENSIONAL, _except_ for out.!        

  integer, intent(in) :: nlevs
  integer, parameter :: nx = 2*kmax, ny = 2*lmax
  complex, dimension(nx,ny,nlevs), intent(inout) :: u,v,t
  complex, dimension(nx,ny), intent(in) :: ub,vb,tb,ut,vt,tt
  real, dimension(nx,ny,nlevs), intent(out) :: out
  !
  complex, dimension(nx,ny,nlevs) :: uz,vz,tz
  complex, dimension(2*kmax,2*lmax) :: dx,dy,Id,vx,uy,tx,ty
  complex, dimension(nx,ny) :: tmp
  real :: dz
  integer :: i,j,k

  ! vertical grid distance
  dz = ZH/real(nlevs); if (verbose .gt. 1) print*,'dz = ',dz

  ! derivative operators
  call d_setup(dx,dy,Id)

  ! zero mean theta
  !tb(1,1) = 0.; tt(1,1) = 0.; t(1,1,:) = 0.

  ! vertical derivatives by centered finite differences, except near boundaries.
  ! the near-boundary points use centered differencing and averaging to get the 
  ! mid-grid points (e.g. grid level 1.5).
  do k=1,nlevs
     if (k .eq. 1) then 
        vz(:,:,k) = (v(:,:,k+1) + v(:,:,k) - 2.*vb(:,:)) / (2.*dz)
        uz(:,:,k) = (u(:,:,k+1) + u(:,:,k) - 2.*ub(:,:)) / (2.*dz)
        tz(:,:,k) = (t(:,:,k) - tb(:,:)) / dz
     elseif (k .eq. nlevs) then 
        vz(:,:,k) = (2.*vt(:,:) - v(:,:,k-1) - v(:,:,k)) / (2.*dz)
        uz(:,:,k) = (2.*ut(:,:) - u(:,:,k-1) - u(:,:,k)) / (2.*dz)
        tz(:,:,k) = (tt(:,:) - t(:,:,k-1)) / dz
     else
        vz(:,:,k) = (v(:,:,k+1)-v(:,:,k-1)) / (2.*dz)
        uz(:,:,k) = (u(:,:,k+1)-u(:,:,k-1)) / (2.*dz)
        tz(:,:,k) = (t(:,:,k) - t(:,:,k-1)) / dz
     endif
     if (verbose .gt. 1) print*,'max sp vz = ',maxval(abs(vz(:,:,k)))
  enddo

  ! spectral horizontal derivatives and epv at each level
  out = 0.
  do k=1,nlevs
     call d_s2s(v(:,:,k),vx,1,dx) 
     call d_s2s(u(:,:,k),uy,1,dy) 
     call d_s2s(t(:,:,k),tx,1,dx) 
     call d_s2s(t(:,:,k),ty,1,dy) 

     tmp = vx - uy + tz(:,:,k) + &
          & Ross*( (vx - uy)*tz(:,:,k) - vz(:,:,k)*tx + uz(:,:,k)*ty )
     ! map back to physical space
     call sp_to_xy(tmp,out(:,:,k),kmax,lmax,nx,ny)
     !call sp_to_xy(ty(:,:),out(:,:,k),kmax,lmax,nx,ny)
     if (verbose .gt. 1) print*,'max recovered pv at level ',k,' : ',maxval(abs(out(:,:,k)))

  enddo

  return
end subroutine epv_pseudoheight
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine epv_pseudoheight_new(nlevs,u,v,t,tb,tt,phi0,phi0b,phi0t,out)

  USE spectral

  IMPLICIT NONE

  ! this version should be asymptotically consistent...

  ! source: svn+ssh://modon.atmos.washington.edu/home/disk/modon/hakim/svn/pvinv
  !
  ! Originator: G. J. Hakim (University of Washington); hakim@atmos.washington.edu
  !
  !
  ! inputs: horizontal winds on grid levels: u,v.
  !         horizontal winds on grid boundaries: ub,vb.
  !         potential temperature on grid levels: t
  !         potential temperature on grid boundaries: tb
  ! output: Ertel potential vorticity: out 
  !
  ! epv = v_x - u_y + t_z + Ross * [(v_x - u_y) * t_z - v_z * t_x + u_z * t_y]
  !
  ! all variables are SPECTRAL and NON-DIMENSIONAL, _except_ for out.!        

  integer, intent(in) :: nlevs
  integer, parameter :: nx = 2*kmax, ny = 2*lmax
  complex, dimension(nx,ny,nlevs), intent(inout) :: u,v,t,phi0
  complex, dimension(nx,ny), intent(inout) :: tb,tt,phi0b,phi0t
  real, dimension(nx,ny,nlevs), intent(out) :: out
  !
  complex, dimension(nx,ny,nlevs) :: uz,vz,tz,u0,v0
  complex, dimension(2*kmax,2*lmax) :: dx,dy,Id,vx,uy,tx,ty,v0x,u0y,u0b,u0t,v0b,v0t
  complex, dimension(nx,ny) :: tmp
  real :: dz
  integer :: i,j,k

  ! vertical grid distance
  dz = ZH/real(nlevs); if (verbose .gt. 1) print*,'dz = ',dz

  ! derivative operators
  call d_setup(dx,dy,Id)

  ! zero mean theta
  tb(1,1) = 0.; tt(1,1) = 0.; t(1,1,:) = 0.

  ! leading order winds
  do k = 1,nlevs
     call d_s2s(phi0(:,:,k),u0(:,:,k),1,-dy)
     call d_s2s(phi0(:,:,k),v0(:,:,k),1,dx)
  enddo
  
  ! leading-order boundary winds
  call d_s2s(phi0b,u0b,1,-dy); call d_s2s(phi0t,u0t,1,-dy)
  call d_s2s(phi0b,v0b,1,dx);  call d_s2s(phi0t,v0t,1,dx)
  
  ! vertical derivatives
  do k=1,nlevs
     if (k .eq. 1) then 
        vz(:,:,k) = (v0(:,:,k+1) + v0(:,:,k) - 2.*v0b(:,:)) / (2.*dz)
        uz(:,:,k) = (u0(:,:,k+1) + u0(:,:,k) - 2.*u0b(:,:)) / (2.*dz)
        tz(:,:,k) = (t(:,:,k) - tb(:,:)) / dz
     elseif (k .eq. nlevs) then 
        vz(:,:,k) = (2.*v0t(:,:) - v0(:,:,k-1) - v0(:,:,k)) / (2.*dz)
        uz(:,:,k) = (2.*u0t(:,:) - u0(:,:,k-1) - u0(:,:,k)) / (2.*dz)
        tz(:,:,k) = (tt(:,:) - t(:,:,k-1)) / dz
     else
        vz(:,:,k) = (v0(:,:,k+1)-v0(:,:,k-1)) / (2.*dz)
        uz(:,:,k) = (u0(:,:,k+1)-u0(:,:,k-1)) / (2.*dz)
        tz(:,:,k) = (t(:,:,k) - t(:,:,k-1)) / dz
     endif
  enddo

  ! spectral horizontal derivatives and epv at each level
  out = 0.
  do k=1,nlevs
     call d_s2s(v(:,:,k),vx,1,dx) 
     call d_s2s(u(:,:,k),uy,1,dy) 
     call d_s2s(v0(:,:,k),v0x,1,dx) 
     call d_s2s(u0(:,:,k),u0y,1,dy) 
     call d_s2s(t(:,:,k),tx,1,dx) 
     call d_s2s(t(:,:,k),ty,1,dy) 
     ! linear piece
     tmp = vx - uy + tz(:,:,k) 
     ! nonlinear piece
     !tmp = tmp + Ross*( (v0x - u0y)*tz(:,:,k) - vz(:,:,k)*tx + uz(:,:,k)*ty )
     ! map back to physical space
     call sp_to_xy(tmp,out(:,:,k),kmax,lmax,nx,ny)
     !call sp_to_xy(vx-uy,out(:,:,k),kmax,lmax,nx,ny)
     !call sp_to_xy(tz(:,:,k),out(:,:,k),kmax,lmax,nx,ny)
     print*,'--max recovered pv at level ',k,' : ',maxval(abs(out(:,:,k)))
     !print*,'max recovered pv at level ',k,' : ',minval(out(:,:,k))
  enddo

  return
end subroutine epv_pseudoheight_new
