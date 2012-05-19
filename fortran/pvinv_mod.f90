!============================================================
! <next few lines under version control, D O  N O T  E D I T>
! $Date$
! $Author$
! $Revision$
! $Id$
!============================================================

!========================================================================
! QG+1 PV inversion module
!========================================================================

!========================================================================
MODULE pvinv_mod

    USE netcdf
    USE spectral_mod

CONTAINS

!========================================================================
SUBROUTINE pvinv_main()

    implicit none

    real, dimension(2*kmax,2*lmax)      :: itbxy,ittxy
    real, dimension(2*kmax,2*lmax,pmax) :: ipvxy

    complex, dimension(2*kmax,2*lmax)      :: phib,phit,ub,ut,vb,vt,tbsp,ttsp
    complex, dimension(2*kmax,2*lmax,pmax) :: phi,u,v,theta,w,pvsp

    complex, dimension(2*kmax,2*lmax)      :: phi0b,phi0t,phi1b,phi1t
    complex, dimension(2*kmax,2*lmax,pmax) :: F1,G1,P1
    complex, dimension(2*kmax,2*lmax,pmax) :: phi0,phi1
    complex, dimension(2*kmax,2*lmax,pmax) :: wcont

    real, dimension(2*kmax,2*lmax,pmax) :: epv,vort,div

    integer :: k, itim, ntims

    character(len=64), parameter :: inpfile = 'pvInput.nc'
    character(len=64), parameter :: outfile = 'pvOutput.nc'

    real, dimension(2*kmax,2*lmax,2),    parameter :: blank2 = 0.0
    real, dimension(2*kmax,2*lmax,pmax), parameter :: blank3 = 0.0

    ntims = get_ntimes(inpfile)

    call write_diag(outfile,0,blank2,blank2,blank2,blank2,blank3,blank3,blank3,blank3,blank3,blank3)

    do itim = 1, ntims

    call init(inpfile,itim,itbxy,ittxy,ipvxy)

    ! max values
    if (verbose .gt. 0)  then
    print*,'max |theta(z=0)| = ',maxval(abs(itbxy))
    print*,'max |theta(z=1)| = ',maxval(abs(ittxy))
    print*,'max |pv| = ',maxval(abs(ipvxy))
    endif

    ! map into spectral space at the same resolution
    call xy_to_sp(cmplx(itbxy,0.),tbsp,2*kmax,2*lmax,kmax,lmax)
    call xy_to_sp(cmplx(ittxy,0.),ttsp,2*kmax,2*lmax,kmax,lmax)
    do k=1,pmax
        call xy_to_sp(cmplx(ipvxy(:,:,k),0.),pvsp(:,:,k),2*kmax,2*lmax,kmax,lmax)
    enddo

    if (verbose .gt. 0) print*,'inverting for leading order geopotential...'
    call invert(tbsp,ttsp,pvsp,phi0b,phi0t,phi0)

    phi = phi0; phib = phi0b; phit = phi0t
    if (verbose .gt. 1) print*,'max nondimensional leading-order phi = ',maxval(abs(phi))

    if (order .eq. 1 .or. order .eq. 2 .and. Ross .ne. 0.) then

    if (verbose .gt. 0) print*,'inverting for O(R) fields...'
    call invert_R(tbsp,ttsp,phi0b,phi0t,phi0,F1,G1,P1,phi1,phi1b,phi1t,pvsp)

    if (order .eq. 1) then 
        phi = Ross*phi1; phib = Ross*phi1b; phit = Ross*phi1t
    else
        phi = phi + Ross*phi1
        phib = phib + Ross*phi1b
        phit = phit + Ross*phi1t
    endif

    endif

    if (verbose .gt. 0) print*,'calculating u and v...'
    call uvwtp(tbsp,ttsp,phi0,phi0b,phi0t,F1,G1,P1,phi1,u,ub,ut,v,vb,vt,w,theta)
    if (verbose .gt. 1) print*,'max |u| = ',maxval(abs(u))
    if (verbose .gt. 1) print*,'max |v| = ',maxval(abs(v))

    ! check w by vertically integrating the continuity equation.
    call w_check(u,v,w,wcont)

    ! check Ertel PV
    call epv_pseudoheight(u,v,theta,ub,vb,tbsp,ut,vt,ttsp,epv)
    !call epv_pseudoheight_new(u,v,theta,tbsp,ttsp,phi0,phi0b,phi0t,epv)
    do k = 1,pmax
    !    if (verbose .gt. 1) print*,'pv check at level ',k,': ',maxval(abs(epv(:,:,k)-ipvxy(:,:,k)))
        if (verbose .gt. 1) print*,'pv check at level ',k,': ',maxval(ipvxy(:,:,k))
    enddo

    ! compute vorticity and divergence
    call vortdiv(u,v,vort,div)
    if (verbose .gt. 1) print *,'max surface vorticity, divergence:',maxval(vort(:,:,2)),maxval(div(:,:,2))

    call dump(u,v,theta,phi,w,pvsp,ub,vb,ut,vt,tbsp,ttsp,phib,phit,itim,outfile)

enddo

END SUBROUTINE pvinv_main

SUBROUTINE invert(tbsp,ttsp,pvsp,phi0b,phi0t,phi0)
! Invert PV for streamfunction; compute spectral derivatives on 
! the transform (advection) grid.

    implicit none

    complex, dimension(2*kmax,2*lmax),      intent(in)  :: tbsp,ttsp
    complex, dimension(2*kmax,2*lmax,pmax), intent(in)  :: pvsp
    complex, dimension(2*kmax,2*lmax),      intent(out) :: phi0b,phi0t
    complex, dimension(2*kmax,2*lmax,pmax), intent(out) :: phi0

    complex, dimension(2*kmax,2*lmax)        :: pxx,pyy,pzz,errsp
    complex, dimension(2*kmax,2*lmax,pmax)   :: bpvsp,lap
    complex, dimension(2*kmax,2*lmax,pmax-1) :: pz

    real :: ak,bl,dz
    integer :: j,k,l,kk,ll

    real,    dimension(pmax) :: facz,faczo,zl,zhl
    complex, dimension(pmax) :: psi
    real,    dimension(2*kmax,2*lmax)       :: errxy
    complex, dimension(2*kmax,2*lmax), save :: dx,dy,Id

    dz = ZH/real(pmax)

    ! derivative operators
    call d_setup(dx,dy,Id)

    ! factors for tanh vertical grid
    call setup_tanh(facz,faczo,zl,zhl)

    ! initialize arrays to zero:
    phi0b=0.;phi0t=0.;phi0=0.

    bpvsp = pvsp
    do k=1,2*kmax; do l=1,2*lmax

        !     add boundary theta to a copy of spectral pv:
         bpvsp(k,l,1) = bpvsp(k,l,1) + (tbsp(k,l)/(facz(1)*dz))
         bpvsp(k,l,pmax) = bpvsp(k,l,pmax) - (ttsp(k,l)/(facz(pmax)*dz))

         !     get wavenumbers
         call get_waves(k,l,ak,bl)

         !     spectral inversion (zero mean)
         if (k .eq. 1 .and. l .eq. 1) then
            phi0b(k,l) = 0.; phi0t(k,l) = 0.; phi0(k,l,:) = 0.
         else
            call matinv(bpvsp(k,l,:),ak,bl,psi,1)
            phi0b(k,l) = psi(1) - (0.5*dz*facz(1)*tbsp(k,l))
            phi0t(k,l) = psi(pmax) + (0.5*dz*facz(pmax)*ttsp(k,l))
            phi0(k,l,:) = psi(:)
         endif

    enddo; enddo

    ! BEGIN leading-order PV check
    ! phi_xx + phi_yy + phi_zz = q

    call laplacian(phi0,tbsp,ttsp,lap)

    ! loop and print our errors
    do k = 1,pmax

        errsp = pvsp(:,:,k) - lap(:,:,k)
        call sp_to_xy(errsp,errxy,kmax,lmax,2*kmax,2*lmax)
        if (verbose .gt. 1)  print*,'max leading-order solution error at level ',k,' : ',maxval(abs(errxy))

        call sp_to_xy(lap(:,:,k),errxy,kmax,lmax,2*kmax,2*lmax)
        if (verbose .gt. 1) print*,'leading-order lap level ',k,' : ',maxval(abs(errxy))

    enddo
    ! END leading-order PV check

    return
END SUBROUTINE invert

SUBROUTINE get_waves(k,l,ak,bl)
! compute x and y wavenumbers for use in spectral calculations.

    implicit none

    integer, intent(in)  :: k,l
    real,    intent(out) :: ak, bl

    ak = facx*real(k - 1); bl = facy*real(l - 1)

    ! other spectral quadrants
    if (k .ge. kmax .and. l .le. lmax) then 
        ak = -1.*facx*real(2*kmax - k + 1)
    elseif (l .ge. lmax .and. k .le. kmax) then 
        bl = -1.*facy*real(2*lmax -l + 1)
    elseif (k .ge. kmax .and. l .ge. lmax) then 
        ak = -1.*facx*real(2*kmax - k + 1)
        bl = -1.*facy*real(2*lmax - l + 1)
    endif

    return
END SUBROUTINE get_waves

SUBROUTINE d_setup(dx,dy,Id)
! Set up matrices for derivatives and integrals.

    implicit none

    complex, dimension(2*kmax,2*lmax), intent(out) :: dx,dy,Id
    real    :: m
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
END SUBROUTINE d_setup

SUBROUTINE d_s2s(temp_in,temp_out,dflag,dn)
! small array to small array.
! dflag =  n: n derivatives. dflag = -n: n integrations.

    implicit none

    complex, dimension(:,:), intent(in)  :: temp_in,dn
    integer,                 intent(in)  :: dflag
    complex, dimension(:,:), intent(out) :: temp_out

    integer :: k,l

    temp_out = 0.

    do k = 1,2*kmax; do l = 1,2*lmax
        if (dn(k,l) .ne. 0) temp_out(k,l) = (dn(k,l)**dflag)*temp_in(k,l)
    enddo; enddo

    return
END SUBROUTINE d_s2s

SUBROUTINE d_s2b(temp_in,temp_out,dflag,dn)
! small array to big array.
! dflag =  n: n derivatives. dflag = -n: n integrations.

    implicit none

    complex, dimension(:,:), intent(in)  :: temp_in,dn
    integer,                 intent(in)  :: dflag
    complex, dimension(:,:), intent(out) :: temp_out

    integer :: k,l,kk,ll

    temp_out = 0.

    do k = 1,kmax; do l = 1,lmax

        ! waves: 0 <= k <= +/-kmax; 0 <= l <= +/-lmax:
        temp_out(k,l) = (dn(k,l)**dflag)*temp_in(k,l)

        kk = kmax + k; ll = lmax + l

        ! +/-kmax <= k <= -1; 0 <= l <= +/-lmax:
        if (k .gt. 1) temp_out(k2+k,l) = (dn(kk,l)**dflag)*temp_in(kk,l)

        ! 0 <= k <= +/-kmax; +/-lmax <= l <= -1:
        if (l .gt. 1) temp_out(k,l2+l) = (dn(k,ll)**dflag)*temp_in(k,ll)

        ! +/-kmax <= k <= -1; +/-lmax <= l <= -1:
        if ((k .gt. 1) .and. (l .gt. 1)) temp_out(k2+k,l2+l) = (dn(kk,ll)**dflag)*temp_in(kk,ll)

    enddo; enddo

    return
END SUBROUTINE d_s2b

SUBROUTINE matinv(pv,ak,bl,psi,idn)
! Solve for psi(z) given pv(z) for a given (k,l) in ** uniform ** N^2:

    implicit none

    complex, intent(in)  :: pv(pmax)
    real,    intent(in)  :: ak,bl
    integer, intent(in)  :: idn
    complex, intent(out) :: psi(pmax)

    complex, dimension(pmax) :: f
    real,    dimension(pmax) :: e

    real    :: A,B,C,dz
    integer :: j

    psi = 0.

    ! set up coeficients A, B, and C (note they are set for uniform n^2!!!):
    dz = ZH/real(pmax)
    A = -1./(dz*dz); B = -1.*((ak*ak) + (bl*bl) + (2./(dz*dz)))
    C = -1./(dz*dz)

    ! idn determines inversion BCs (1.=Neumann, -1.=Dirichlet)
    !      idn = 1

    ! first pass:
    e(1) = A / (B - (C*real(idn)))
    f(1) = pv(1) / (B-(C*real(idn)))

    do j = 2, pmax 
        e(j) = A / (B - (C*e(j-1)))
        f(j) = (pv(j) + (C*f(j-1))) / (B - (C*e(j-1)))
        !print*,'deep check=',j,f(j),pmax
    enddo

    !print*,'f = ',f(pmax)
    !second pass
    if(e(pmax).eq.idn) then
        psi(pmax) = 0.0
    else
        psi(pmax) = f(pmax) / (1. - (real(idn)*e(pmax)))  
        !print*,'deep check=',f(pmax),real(idn),e(pmax)
    endif

    do j = 1, pmax-1                      
        psi(pmax-j) = (e(pmax-j)*psi(pmax+1-j)) + f(pmax-j)
    enddo

    return
END SUBROUTINE matinv
!========================================================================

!========================================================================
SUBROUTINE xy_to_sp(xy,sp,mx,ny,km,lm)
! Map an (x,y) array onto a _smaller_ spectral array.
! Input: xy(mx,ny) --- a grid point array.
! Output: sp(2*km,2*lm) --- a spectral array.

    implicit none

    integer,                       intent(in)  :: km,lm,mx,ny
    complex, dimension(:,:),       intent(in)  :: xy
    complex, dimension(2*km,2*lm), intent(out) :: sp

    complex, dimension(mx,ny) :: copy
    integer                   :: kmp1,lmp1,k,l,k2,l2,kk,ll

    ! initialize arrays:
    sp = 0.0 ; copy = xy

    kmp1 = km + 1; lmp1 = lm + 1; k2 = mx - km; l2 = ny - lm

    call ft_2d(copy,mx,ny,-1)

    !do k=1,kmp1; do l=1,lmp1
    do k=1,km; do l=1,lm

        kk = km + k; ll = lm + l

        ! waves: 0 <= k <= +/-kmax; 0 <= l <= +/-lmax:
        sp(k,l) = copy(k,l)/real(mx*ny)

        ! +/-kmax <= k <= -1; 0 <= l <= +/-lmax:
        if (k .gt. 1) then
            sp(kk,l) = copy(k2+k,l)/real(mx*ny)
        endif

        ! 0 <= k <= +/-kmax; +/-lmax <= l <= -1:
        if (l .gt. 1) then
            sp(k,ll) = copy(k,l2+l)/real(mx*ny)
        endif

        ! +/-kmax <= k <= -1; +/-lmax <= l <= -1:
        if ((k .gt. 1) .and. (l .gt. 1)) then
            sp(kk,ll) = copy(k2+k,l2+l)/real(mx*ny)
        endif

    enddo; enddo

    return
END SUBROUTINE xy_to_sp
!========================================================================

!========================================================================
SUBROUTINE sp_to_xy(sp,xy,km,lm,mx,ny)
! Map an (km,lm) spectral array onto a _bigger_ grid point array.
! Input: sp(2*km,2*lm) --- a spectral array.
! Output: xy(mx,ny) --- a grid point array.

    implicit none

    integer,                       intent(in)  :: km,lm,mx,ny
    complex, dimension(:,:),       intent(in)  :: sp
    real,    dimension(mx,ny),     intent(out) :: xy

    complex :: copy(mx,ny)
    integer :: kmp1,lmp1,k,l,k2,l2,kk,ll

    copy = 0.0

    kmp1 = km + 1; lmp1 = lm + 1; k2 = mx - km; l2 = ny - lm

    !do 10 k = 1,kmp1; do 10 l = 1,lmp1
    do k = 1,km ; do l = 1,lm

        ! waves: 0 <= k <= +/-kmax; 0 <= l <= +/-lmax:
        copy(k,l) = sp(k,l)

        kk = km + k; ll = lm + l

        ! +/-kmax <= k <= -1; 0 <= l <= +/-lmax:
        if (k .ge. 1) then
           copy(k2+k,l) = sp(kk,l)
        endif

        ! 0 <= k <= +/-kmax; +/-lmax <= l <= -1:
        if (l .ge. 1) then
           copy(k,l2+l) = sp(k,ll)
        endif

        ! +/-kmax <= k <= -1; +/-lmax <= l <= -1:
        if ((k .ge. 1) .and. (l .ge. 1)) then
           copy(k2+k,l2+l) = sp(kk,ll)
        endif

     enddo ; enddo

    call ft_2d(copy,mx,ny,1)
    xy = real(copy)

    return
END SUBROUTINE sp_to_xy
!========================================================================

!========================================================================
SUBROUTINE init(infile,itime,tbxy,ttxy,pvxy)
! Initialize pv and theta fields

  implicit none

  character(len=*),                    intent(in)  :: infile
  integer,                             intent(in)  :: itime
  real, dimension(2*kmax,2*lmax),      intent(out) :: tbxy,ttxy
  real, dimension(2*kmax,2*lmax,pmax), intent(out) :: pvxy

  real, dimension(2*kmax,2*lmax,2) :: theta

  real :: x,y,z,dx,dy,dz
  real :: xcen,ycen,amdx,amdy,rr,rrz
  real :: xn
  integer :: i,j,k
  integer :: iseed
  integer :: ncid, varid, start(4), count(4)

  real, parameter :: pva = 0.0
  real, parameter :: amp = -2.0, asig = 0.5
  real, parameter :: asx = 1.0,  asy  = 1.0, asz = 0.1
  real, parameter :: zlev = 0.5
 
  if (itime == 0) then ! make IC
  
     dx = XL / real(2*kmax)
     dy = YL / real(2*lmax)
     dz = ZH / real(pmax)

     do i=1,2*kmax; do j=1,2*lmax

        x = (i-1)*dx; y = (j-1)*dy

        xcen = kmax*dx; ycen = lmax*dy
        amdx=min( abs(x-xcen),abs(XL+x-xcen),abs(XL-x+xcen) )
        amdy=min( abs(y-ycen),abs(YL+y-ycen),abs(YL-y+ycen) )

        rr = (((amdx/asx)**2) + ((amdy/asy)**2))**0.5
        
        ! gaussian:
        tbxy(i,j) = 0.
        ttxy(i,j) = amp*exp(-1.0 * ((rr/asig)**2))

        ! second deriv of gaussian:
        !tbxy(i,j) = 0.
        !ttxy(i,j) = -amp*(1.-(rr*rr))*exp(-1.*((rr/asig)**2))

        ! random:
        !tbxy(i,j) = amp*(ran1(iseed)-0.5)
        !ttxy(i,j) = amp*(ran1(iseed)-0.5)

        ! plane waves
        !tbxy(i,j) = 0.
        !ttxy(i,j) = amp*cos((2.*pi*x/XL) ) ! + (2.*pi*y/YL))

        do k = 1, pmax

           Z = (real(k) - 0.5)*dz

           xn = xcen !- (z - 1.) ! this term adds vertical tilt

           amdx=min( abs(x-xn),abs(XL+x-xn),abs(XL-x+xn) )
           rr = (((amdx/asx)**2) + ((amdy/asy)**2))**0.5 
           rrz = abs(zlev-z)
           pvxy(i,j,k) = pva * (1.0 - (rr*rr))        * &
                           exp(-1.0 * ((rr/asig)**2)) * &
                           exp(-1.0 * ((rrz/asz)**2))

           if ( ( i.eq.kmax) .and. (j.eq.lmax) ) then 
              if (verbose .gt. 1) print*,'pv: ',k,pvxy(i,j,k)
           endif

        enddo

     enddo; enddo

  else ! read from a file

     tbxy = 0.0 ; ttxy = 0.0 ; pvxy = 0.0

     start(1) = 1     ; count(1) = 2*kmax
     start(2) = 1     ; count(2) = 2*lmax
     start(3) = 1     ; count(3) = 2
     start(4) = itime ; count(4) = 1

     call nc_check( nf90_open(trim(infile), NF90_NOWRITE, ncid), 'init', 'open, ' // trim(infile) )
     call nc_check( nf90_inq_varid(ncid, 'theta', varid), 'init', 'inq_varid theta, ' // trim(infile) )
     call nc_check( nf90_get_var(ncid, varid, theta, start, count), 'init', 'get_var theta, ' // trim(infile) )
     call nc_check( nf90_close(ncid), 'init', 'close, ' // trim(infile) )

     tbxy = theta(:,:,1) ; ttxy = theta(:,:,2)

  endif

  return
END SUBROUTINE init
!========================================================================

!========================================================================
FUNCTION ran1(idum)
! Numerical Recipes random number generator:

    implicit none

    REAL    ran1
    INTEGER idum

    INTEGER IA,IM,IQ,IR,NTAB,NDIV
    REAL    AM,EPS,RNMX
    PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836, &
               NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
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
!========================================================================

!========================================================================
SUBROUTINE dump(u,v,theta,phi,w,pv,ub,vb,ut,vt,tb,tt,phib,phit,it,outfile)
! Write to disk.

    implicit none

    complex, dimension(2*kmax,2*lmax,pmax), intent(in) :: u,v,theta,phi,w,pv
    complex, dimension(2*kmax,2*lmax),      intent(in) :: ub,vb,ut,vt,tb,tt,phib,phit
    integer,                                intent(in) :: it
    character(len=*),                       intent(in) :: outfile

    real, dimension(pmax)               :: pheight
    real, dimension(2*kmax,2*lmax,pmax) :: otheta,ophi,ou,ov,ow,opv
    real, dimension(2*kmax,2*lmax,2)    :: otheta_s,ophi_s,ou_s,ov_s
    integer :: k
    real :: dz
    real :: grav,tnot,Rhos,Ns,km,Cor,Ls,Hs,Us,Ws,Ps,Ts,PVs

    call scaling(grav,tnot,Rhos,Ns,km,Cor,Ls,Hs,Us,Ws,Ps,Ts,PVs)

    dz = ZH/real(pmax)
    if (verbose .gt. 1) then 
    print*,'Scaling parameters...'
    print*,'grav = ',grav,' m s^-2'
    print*,'Cor  = ',Cor,' s^-1'
    print*,'Ls   = ',Ls, ' m'
    print*,'Hs   = ',Hs, ' m'
    print*,'Us   = ',Us, ' m s^-1'
    print*,'Ws   = ',Ws*100, ' cm s^-1'
    print*,'Ps   = ',Ps, ' m^2 s^-2'
    print*,'Ts   = ',Ts, ' K'
    print*,'Ns   = ',Ns, ' s^-1'
    print*,'PVs  = ',PVs, ' PVU'
    endif

    ! lower boundary theta
    call sp_to_xy(tb,otheta_s(:,:,1),kmax,lmax,2*kmax,2*lmax); otheta_s(:,:,1) = Ts*otheta_s(:,:,1)
    if (verbose .gt. 1) print*,'max z=0 dim |t| = ',maxval(abs(otheta_s(:,:,1)))

    ! upper boundary theta
    call sp_to_xy(tt,otheta_s(:,:,2),kmax,lmax,2*kmax,2*lmax); otheta_s(:,:,2) = Ts*otheta_s(:,:,2)
    if (verbose .gt. 1) print*,'max z=H dim |t| = ',maxval(abs(otheta_s(:,:,2)))

    ! lower boundary geopotential
    call sp_to_xy(phib,ophi_s(:,:,1),kmax,lmax,2*kmax,2*lmax); ophi_s(:,:,1) = Ps*ophi_s(:,:,1)
    if (verbose .gt. 1) print*,'max z=0 dim |p| = ',maxval(abs(ophi_s(:,:,1)))

    ! upper boundary geopotential
    call sp_to_xy(phit,ophi_s(:,:,2),kmax,lmax,2*kmax,2*lmax); ophi_s(:,:,2) = Ps*ophi_s(:,:,2)
    if (verbose .gt. 1) print*,'max z=0 dim |p| = ',maxval(abs(ophi_s(:,:,2)))

    ! lower boundary u & v
    call sp_to_xy(ub,ou_s(:,:,1),kmax,lmax,2*kmax,2*lmax); ou_s(:,:,1) = Us*ou_s(:,:,1)
    if (verbose .gt. 1) print*,'max z=0 dim |u| = ',maxval(abs(ou_s(:,:,1)))
    call sp_to_xy(vb,ov_s(:,:,1),kmax,lmax,2*kmax,2*lmax); ov_s(:,:,1) = Us*ov_s(:,:,1)
    if (verbose .gt. 1) print*,'max z=0 dim |v| = ',maxval(abs(ov_s(:,:,1)))

    ! upper boundary u & v
    call sp_to_xy(ut,ou_s(:,:,2),kmax,lmax,2*kmax,2*lmax); ou_s(:,:,2) = Us*ou_s(:,:,2)
    if (verbose .gt. 1) print*,'max z=H dim |u| = ',maxval(abs(ou_s(:,:,2)))
    call sp_to_xy(vt,ov_s(:,:,2),kmax,lmax,2*kmax,2*lmax); ov_s(:,:,2) = Us*ov_s(:,:,2)
    if (verbose .gt. 1) print*,'max z=H dim |v| = ',maxval(abs(ov_s(:,:,2)))

    ! interior theta, geopotential, u, v, w and pv
    do k = 1, pmax

        if (verbose .gt. 1) print*,'level ',k

        call sp_to_xy(theta(:,:,k),otheta(:,:,k),kmax,lmax,2*kmax,2*lmax)
        otheta(:,:,k) = Ts*otheta(:,:,k)
        if (verbose .gt. 1) print*,'max dim |t| = ',maxval(abs(otheta(:,:,k)))

        call sp_to_xy(phi(:,:,k),ophi(:,:,k),kmax,lmax,2*kmax,2*lmax) 
        ophi(:,:,k) = Ps*ophi(:,:,k)
        if (verbose .gt. 1) print*,'max dim |p| = ',maxval(abs(ophi(:,:,k)))

        call sp_to_xy(u(:,:,k),ou(:,:,k),kmax,lmax,2*kmax,2*lmax)
        ou(:,:,k) = Us*ou(:,:,k)
        if (verbose .gt. 1) print*,'max dim |u| = ',maxval(abs(ou(:,:,k)))

        call sp_to_xy(v(:,:,k),ov(:,:,k),kmax,lmax,2*kmax,2*lmax)
        ov(:,:,k) = Us*ov(:,:,k)
        if (verbose .gt. 1) print*,'max dim |v| = ',maxval(abs(ov(:,:,k)))

        call sp_to_xy(w(:,:,k),ow(:,:,k),kmax,lmax,2*kmax,2*lmax)
        ow(:,:,k) = Ws*ow(:,:,k)
        if (verbose .gt. 1) print*,'max dim |w| = ',maxval(abs(ow(:,:,k)))

        call sp_to_xy(pv(:,:,k),opv(:,:,k),kmax,lmax,2*kmax,2*lmax); 
        opv(:,:,k) = PVs*opv(:,:,k)
        if (verbose .gt. 1) print*,'max dim |pv| = ',maxval(abs(opv(:,:,k)))

        pheight(k) = Hs*(real(k) - 0.5)*dz

    enddo

    if (verbose .gt. 0) print*,'Writing to disk...'
    call write_diag(outfile,it,ou_s,ov_s,otheta_s,ophi_s,ou,ov,otheta,ophi,ow,opv)

    return
END SUBROUTINE dump
!========================================================================

!========================================================================
SUBROUTINE ft_2d(f,ni,nj,isign)
! FFT-calling SUBROUTINE.

    implicit none

    integer,                   intent(in)    :: ni, nj, isign
    complex, dimension(ni,nj), intent(inout) :: f

    real, dimension(ni,nj) :: re, im

    re = real(f); im = aimag(f)
    call fft(re,im,ni*nj,ni,ni,   isign)
    call fft(re,im,ni*nj,nj,nj*ni,isign)
    f = cmplx(re,im)

    return
END SUBROUTINE ft_2d
!========================================================================

!========================================================================
SUBROUTINE invert_old(tbsp,ttsp,pvsp,phi0b,phi0t,phi0)
! this is the old inversion routine. the new one has the quadrant work 
! hidden in operators, so it scales better for +1 code.

  implicit none

  ! Invert PV for streamfunction; compute spectral derivatives on 
  ! the transform (advection) grid.

  complex, intent(in), dimension(2*kmax,2*lmax) :: tbsp,ttsp
  complex, intent(in), dimension(2*kmax,2*lmax,pmax) :: pvsp
  complex, intent(out), dimension(2*kmax,2*lmax) :: phi0b,phi0t
  complex, intent(out), dimension(2*kmax,2*lmax,pmax) :: phi0
  complex, dimension(2*kmax,2*lmax,pmax) :: bpvsp      
  real :: ak,bl,kap,sqg,amod,dz
  integer :: j,k,l,kk,ll
  complex, dimension(pmax) :: psi

  dz = ZH/real(pmax)

  ! initialize arrays to zero:
  phi0b=0.;phi0t=0.;phi0=0.

  ! add boundary theta to a copy of spectral pv:
  bpvsp = pvsp
  do k=1,2*kmax; do l=1,2*lmax
     bpvsp(k,l,1) = bpvsp(k,l,1) + (tbsp(k,l)/dz)
     bpvsp(k,l,pmax) = bpvsp(k,l,pmax) - (ttsp(k,l)/dz)
  enddo; enddo

  ! inversion of (kmax,lmax) waves into (mmax,nmax) arrays:
  do k = 1,kmaxp1
     do l = 1,lmaxp1

        ! waves: 0 <= k <= +/-kmax; 0 <= l <= +/-lmax:
        ak = facx*real(k - 1); bl = facy*real(l - 1)

        if (k .eq. 1 .and. l .eq. 1) then
           phi0(k,l,:) = 0.; phi0(k,l,:) = 0.
        else
           call matinv(bpvsp(k,l,:),ak,bl,psi,1)
           phi0b(k,l) = psi(1) - (0.5*dz*tbsp(k,l))
           phi0t(k,l) = psi(pmax) + (0.5*dz*ttsp(k,l))
           phi0(k,l,:) = psi
        endif

        kk = kmax + k; ll = lmax + l

        ! +/-kmax <= k <= -1; 0 <= l <= +/-lmax:
        if (k .le. kmax) then
           ak = -1.*facx*real(kmaxp1 - k); bl = facy*real(l - 1)
           call matinv(bpvsp(kk,l,:),ak,bl,psi,1)
           phi0b(kk,l) = psi(1) - (0.5*dz*tbsp(kk,l))
           phi0t(kk,l) = psi(pmax) + (0.5*dz*ttsp(kk,l))
           phi0(kk,l,:) = psi
        endif

        ! 0 <= k <= +/-kmax; +/-lmax <= l <= -1:
        if (l .le. lmax) then
           ak = facx*real(k-1); bl = -1.*facy*real(lmaxp1 - l)
           call matinv(bpvsp(k,ll,:),ak,bl,psi,1)
           phi0b(k,ll) = psi(1) - (0.5*dz*tbsp(k,ll))
           phi0t(k,ll) = psi(pmax) + (0.5*dz*ttsp(k,ll))
           phi0(k,ll,:) = psi
        endif

        ! +/-kmax <= k <= -1; +/-lmax <= l <= -1:
        if (k .le. kmax .and. l .le. lmax) then
           ak=-1.*facx*real(kmaxp1 - k); bl=-1.*facy*real(lmaxp1 - l)
           call matinv(bpvsp(kk,ll,:),ak,bl,psi,1)
           phi0b(kk,ll) = psi(1) - (0.5*dz*tbsp(kk,ll))
           phi0t(kk,ll) = psi(pmax) + (0.5*dz*ttsp(kk,ll))
           phi0(kk,ll,:) = psi
        endif
     end do
  end do

  return
END SUBROUTINE invert_old
!========================================================================

!========================================================================
SUBROUTINE uvwtp(tb,tt,phi0,phi0b,phi0t,F1,G1,P1,phi1, &
                 u,ub,ut,v,vb,vt,w,thF)
! recover winds, theta, and geopotential from potential data.

  implicit none

  complex, dimension(2*kmax,2*lmax,pmax), intent(in) :: phi0, &
       &                   F1,G1,P1,phi1
  complex, dimension(2*kmax,2*lmax), intent(in) :: phi0b,phi0t,tb,tt
  complex, dimension(2*kmax,2*lmax,pmax), intent(out) :: u,v,w,thF
  complex, dimension(2*kmax,2*lmax), intent(out) :: ub,ut,vb,vt
  complex, dimension(2*kmax,2*lmax) :: dx,dy,Id,uF,vF,temp,ftemp,oftemp
  complex, dimension(2*kmax,2*lmax,pmax) :: u1,v1
  complex, dimension(2*kmax,2*lmax,pmax) :: phi
  complex, dimension(2*kmax,2*lmax) :: ug,vg
  real :: dz
  integer :: k,l,p
  complex :: gtemp,ogtemp
  ! testing...
  real, dimension(2*kmax,2*lmax) :: rtemp,rutemp,rvtemp,ws,wsg

  dz = ZH/real(pmax)

  ! derivative operators
  call d_setup(dx,dy,Id)

  ! 3D leading-order winds
  if (order .eq. 0 .or. order .eq. 2) then 

     do k = 1,pmax
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
     do p=1,pmax
        if (p .eq. 1) then ! linear interp; zero BC
           u1(:,:,1) = -(F1(:,:,1) + F1(:,:,2))/(2*dz)
           v1(:,:,1) = -(G1(:,:,1) + G1(:,:,2))/(2*dz)
           ub(:,:) = ub(:,:) - Ross*(2*F1(:,:,1)/dz)
           vb(:,:) = vb(:,:) - Ross*(2*G1(:,:,1)/dz)  
        elseif (p .eq. pmax) then ! linear interp; zero BC
           u1(:,:,pmax) = (F1(:,:,pmax) + F1(:,:,pmax-1))/(2*dz)
           v1(:,:,pmax) = (G1(:,:,pmax) + G1(:,:,pmax-1))/(2*dz)
           ut(:,:) = ut(:,:) + Ross*(2*F1(:,:,pmax)/dz)
           vt(:,:) = vt(:,:) + Ross*(2*G1(:,:,pmax)/dz)
        else
           u1(:,:,p) = -(F1(:,:,p+1) - F1(:,:,p-1))/(2*dz)
           v1(:,:,p) = -(G1(:,:,p+1) - G1(:,:,p-1))/(2*dz)
        endif
     enddo ! p 

     !      u1 = 0; v1 = 0 ! test of phi contrib
     ! compute horizonal derivatives spectrally
     do p = 1, pmax
        call d_s2s(P1(:,:,p),temp,1,dx); v1(:,:,p) = v1(:,:,p) + temp
        call d_s2s(P1(:,:,p),temp,1,-dy); u1(:,:,p) = u1(:,:,p) + temp
        if (p .eq. 1) then ! recall phi1_z = 0
           call d_s2s(P1(:,:,p),temp,1,dx); vb = vb + Ross*temp
           call d_s2s(P1(:,:,p),temp,1,-dy); ub = ub + Ross*temp
        endif
        if (p .eq. pmax) then ! recall phi1_z = 0
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
  !do p=1,pmax-1
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
  do p=1,pmax
     !P1_z:
     if (p .eq. 1) then 
        thF(:,:,p) = (phi(:,:,p+1) - (phi(:,:,p)-tb*dz))/(2.*dz)
     elseif (p .eq. pmax) then 
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
  do p=1,pmax-1
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
  do p = 1, pmax
     call d_s2s(phi(:,:,p),vg,1,dx);  call d_s2s(phi(:,:,p),ug,1,-dy);
     ! move geostrophic wind onto physical grid
     rutemp = 0.; rvtemp = 0.
     call sp_to_xy(vg,rvtemp,kmax,lmax,2*kmax,2*lmax)
     call sp_to_xy(ug,rutemp,kmax,lmax,2*kmax,2*lmax)
     !call d_s2s(vg,rvtemp,1,Id); call d_s2s(ug,rutemp,1,Id)
     wsg = sqrt(((rutemp**2) + (rvtemp**2)))
     ! move full wind onto physical grid
     rutemp = 0.; rvtemp = 0.
     call sp_to_xy(v(:,:,p),rvtemp,kmax,lmax,2*kmax,2*lmax)
     call sp_to_xy(u(:,:,p),rutemp,kmax,lmax,2*kmax,2*lmax)
     !call d_s2s(v(:,:,p),rvtemp,1,Id); call d_s2s(u(:,:,p),rutemp,1,Id)
     ws = sqrt(((rutemp**2) + (rvtemp**2)))
     if (verbose .gt. 1) print*,'max wind speed for full and geostrophic winds (total phi) at level: ', p, maxval(ws),maxval(wsg)
  enddo

  return
END SUBROUTINE uvwtp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE invert_R(tbsp,ttsp,phi0b,phi0t,phi0,F1,G1,P1,phi1,phi1b,phi1t,pvsp)
! Invert for O(R) potentials.
! phi0 : leading order geopotential.
! phi1 : O(Ross) geopotential correction.
! phi1b, phi1t : O(Ross) boundary geopotential corrections.
! F1, G1, P1 : 3D O(Ross) potentials.
! pz, pzz, pxy, etc. : derivatives of phi0.
! bpvsp : boundary condition trick (make homogeneous BCs; correct int).

  implicit none

  complex, intent(in), dimension(2*kmax,2*lmax) :: tbsp,ttsp,phi0b,phi0t
  complex, intent(in), dimension(2*kmax,2*lmax,pmax) :: phi0,pvsp
  !
  complex, intent(out), dimension(2*kmax,2*lmax,pmax) :: F1,G1,P1,phi1
  complex, intent(out), dimension(2*kmax,2*lmax) :: phi1b,phi1t
  !
  complex, dimension(2*kmax,2*lmax,pmax) :: F1s,G1s,P1s,phi1s,lap
  real :: ak,bl,dz,solv
  integer :: j,k,l,kk,ll,p
  complex,dimension(2*kmax,2*lmax),save:: dx,dy,Id
  complex, dimension(2*kmax,2*lmax,pmax) :: pzS,pzzS
  complex, dimension(2*kmax,2*lmax) :: bcb,bct,errsp,tb,tt
  real, dimension(mmax,nmax) :: pzx,pzy,pzz,pxx,pyy,pxy,Rtemp,pvxy,tbx,tby,pzxold,pzyold

  complex, dimension(mmax,nmax) :: Ctemp
  complex, dimension(pmax) :: psi
  real, dimension(2*kmax,2*lmax) :: Ptemp,errxy
  complex, dimension(2*kmax,2*lmax) :: solvB,solvT
  ! filler for setup_tanh
  real, dimension(pmax) :: facz,faczo,zl,zhl

  dz = ZH/real(pmax)

  ! derivative operators
  call d_setup(dx,dy,Id)

  ! factors for tanh vertical grid
  call setup_tanh(facz,faczo,zl,zhl)
  
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
  Rtemp = pzz*(pxx + pyy); call xy_to_sp(cmplx(Rtemp,0.),solvB,mmax,nmax,kmax,lmax)
  ! top
  call d_s2b(phi0t,Ctemp,1,dx*dx); call ft_2d(Ctemp,mmax,nmax,1); pxx = real(Ctemp)
  call d_s2b(phi0t,Ctemp,1,dy*dy); call ft_2d(Ctemp,mmax,nmax,1); pyy = real(Ctemp)
  ! this is an estimate of pzz 1/8 dz off the surface; pre-use solvT
  solvT = (tt - ((phi0t - phi0(:,:,pmax))*2./(facz(pmax)*dz)))*4./(facz(pmax)*dz)
  call d_s2b(solvT,Ctemp,1,Id); call ft_2d(Ctemp,mmax,nmax,1); pzz = real(Ctemp)
  Rtemp = pzz*(pxx + pyy); call xy_to_sp(cmplx(Rtemp,0.),solvT,mmax,nmax,kmax,lmax)
  solv = (solvT(1,1) - solvB(1,1))/H
  if (verbose .gt. 1)  print*,'*-*-*-* O(R) solvability constant : ',solv

  ! phi_z at half levels; phi_zz on grid levels
  do k = 1,pmax
     if (k .le. pmax-1) then
        pzS(:,:,k) = (phi0(:,:,k+1) - phi0(:,:,k)) / (faczo(k)*dz)
     endif
     if (k .eq. 1) then
        pzzS(:,:,k) = (pzS(:,:,k) - tb) / (facz(k)*dz)
     elseif (k .eq. pmax) then
        pzzS(:,:,k) = (tt - pzS(:,:,k-1)) / (facz(k)*dz)
     else
        pzzS(:,:,k) = (pzS(:,:,k) - pzS(:,:,k-1)) / (facz(k)*dz)
     endif
  enddo
  pzzS(1,1,:) = 0.
  
  ! terms with mixed derivatives involving d/dz at half levels on big grid
  do k=1,pmax

     pzz=0.;pzx=0.;pzy=0.;pxx=0.;pyy=0.;pxy=0.

     if (k .gt. 2 .and. k .lt. pmax) then
        pzxold = pzx; pzyold = pzy
     endif
     if (k .le. pmax-1) then
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
     elseif (k .eq. pmax) then ! recycle tbx & tby names here
        call d_s2b(tt,Ctemp,1,dx)
        call ft_2d(Ctemp,mmax,nmax,1); tbx = real(Ctemp)
        call d_s2b(tt,Ctemp,1,dy)
        call ft_2d(Ctemp,mmax,nmax,1); tby = real(Ctemp)
        pzx = (pzxold + tbx)/2.; pzy = (pzyold + tby)/2.
     else
        if (k .eq. pmax-1) then 
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
     call xy_to_sp(cmplx(Rtemp,0.),F1s(:,:,k),mmax,nmax,kmax,lmax)
     Rtemp = 2*(pzx*pyy - pzy*pxy)
     call xy_to_sp(cmplx(Rtemp,0.),G1s(:,:,k),mmax,nmax,kmax,lmax)
     Rtemp = (pzx*pzx) + (pzy*pzy) - ((pxx + pyy)*pzz)
!     Rtemp = (pzx*pzx) + (pzy*pzy) + (pzz*pzz) ! test for zero pv
     call xy_to_sp(cmplx(Rtemp,0.),P1s(:,:,k),mmax,nmax,kmax,lmax)
     P1s(1,1,k) = P1s(1,1,k) + solv
     Rtemp = 2*(pxx*pyy - pxy*pxy)
     call xy_to_sp(cmplx(Rtemp,0.),phi1s(:,:,k),mmax,nmax,kmax,lmax)

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
        call matinv(phi1s(k,l,:),ak,bl,psi,1)
        phi1b(k,l) = psi(1)
        phi1t(k,l) = psi(pmax)
        phi1(k,l,:) = psi(:)
        ! F
        call matinv(F1s(k,l,:),ak,bl,psi,-1); F1(k,l,:) = psi(:)
        ! G
        call matinv(G1s(k,l,:),ak,bl,psi,-1); G1(k,l,:) = psi(:)
        ! Phi
        call matinv(P1s(k,l,:),ak,bl,psi,1); P1(k,l,:) = psi(:)
     endif
  enddo; enddo

  ! finally, correct phi from Phi
  phi1 = phi1 + P1
  phi1b = phi1b + P1(:,:,1) ! P1_z = 0--> P1(z=0) = P1(z=dz/2)
  phi1t = phi1t + P1(:,:,pmax)

  !---------- Phi1 solution check
  bcb = 0.; bct = 0.
  call laplacian(P1,bcb,bct,lap)

  do k = 1,pmax
     errsp = P1s(:,:,k) - lap(:,:,k); errsp(1,1) = 0.
     call sp_to_xy(errsp,errxy,kmax,lmax,2*kmax,2*lmax)
     if (verbose .gt. 1) print*,'max P1 solution error at level ',k,' : ',maxval(abs(errxy))

  enddo
  !---------- F1 solution check
  bcb = 2*F1(:,:,1)/dz; bct = -2*F1(:,:,pmax)/dz
  call laplacian(F1,bcb,bct,lap)
  do k = 1,pmax
     errsp = F1s(:,:,k) - lap(:,:,k)
     call sp_to_xy(errsp,errxy,kmax,lmax,2*kmax,2*lmax)
     if (verbose .gt. 1) print*,'max F1 solution error at level ',k,' : ',maxval(abs(errxy))
  enddo

  !---------- G1 solution check
  bcb = 2*G1(:,:,1)/dz; bct = -2*G1(:,:,pmax)/dz
  call laplacian(G1,bcb,bct,lap)
  do k = 1,pmax
     errsp = G1s(:,:,k) - lap(:,:,k)
     call sp_to_xy(errsp,errxy,kmax,lmax,2*kmax,2*lmax)
     if (verbose .gt. 1) print*,'max G1 solution error at level ',k,' : ',maxval(abs(errxy))
  enddo

  return
END SUBROUTINE invert_R


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE w_check(u,v,win,wout)

  ! Originator: G. J. Hakim,  University of Washington

  implicit none

  complex, dimension(2*kmax,2*lmax,pmax), intent(in) :: u,v,win
  complex, dimension(2*kmax,2*lmax,pmax), intent(out) :: wout
  complex, dimension(2*kmax,2*lmax) :: dx,dy,Id
  complex, dimension(2*kmax,2*lmax,pmax) :: w
  complex, dimension(2*kmax,2*lmax) :: temp,wnew,wold,ux,vy,wF
  real, dimension(2*kmax,2*lmax) :: wxo,wxn      
  real :: dz
  integer:: p

  ! vertical grid distance
  dz = ZH/real(pmax)

  ! derivative operators
  call d_setup(dx,dy,Id)

  ! integrate u_x + v_y = -Ross w_z
  ! note that u and v are next-order, so Ross on RHS cancels w/ those that 
  ! are on the LHS.
  wold = 0.
  do p = 1, pmax
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
END SUBROUTINE w_check

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE laplacian(func,fb,ft,lap)

  ! compute an isotropic 3D laplacian given a staggered field & BCs.
  ! Originator: G. J. Hakim,  University of Washington

  implicit none

  complex, dimension(2*kmax,2*lmax,pmax), intent(in) :: func
  complex, dimension(2*kmax,2*lmax), intent(in) :: fb,ft
  complex, dimension(2*kmax,2*lmax,pmax), intent(out) :: lap
  integer :: k
  complex, dimension(2*kmax,2*lmax,pmax-1) :: fz      
  complex, dimension(2*kmax,2*lmax) :: fxx,fyy,fzz      
  complex,dimension(2*kmax,2*lmax) :: dx,dy,Id
  real :: dz

  dz = ZH/real(pmax)

  ! derivative operators
  call d_setup(dx,dy,Id)

  dz = ZH/real(pmax)

  ! first z-loop for f_z at at intermediate levels
  do k = 1,pmax-1
     fz(:,:,k) = (func(:,:,k+1) - func(:,:,k)) / dz
  enddo

  ! second z-loop for laplacian at grid levels and f_zz; lap calc
  do k = 1,pmax
     
     call d_s2s(func(:,:,k),fxx,1,dx*dx)
     call d_s2s(func(:,:,k),fyy,1,dy*dy)

    if (k .eq. 1) then
        fzz = (fz(:,:,k) - fb) / dz
     elseif (k .eq. pmax) then
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
END SUBROUTINE laplacian

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE vortdiv(u,v,vort,div)

  ! given horizontal wind (u,v) compute vertical component of 
  ! vorticity and horizontal divergence. NONDIMENSIONAL I/O.
  ! make these values dimensional by multiplying by U/L
  ! Originator: G. J. Hakim,  University of Washington

  implicit none

  complex, dimension(2*kmax,2*lmax,pmax), intent(in) :: u,v
  real, dimension(2*kmax,2*lmax,pmax), intent(out) :: vort,div
  complex, dimension(2*kmax,2*lmax) :: ux,uy,vx,vy
  complex,dimension(2*kmax,2*lmax) :: dx,dy,Id,vs,ds
  integer :: k

  ! derivative operators
  call d_setup(dx,dy,Id)

  do k = 1,pmax
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
END SUBROUTINE vortdiv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE setup_tanh(facz,faczo,zl,zhl)

  ! setup Jacobian factors for tanh vertical grid

  implicit none

  real, dimension(pmax), intent(out) :: facz, faczo,zl,zhl
  integer :: k
  real :: fac,alpha,acosh,znot,z,arg,dz

  ! if we're using the regular grid, set jacobian values to one:
  if (.not. tanhgrid) then 
     facz = 1.; faczo = 1.
     return
  endif

  dz = ZH/real(pmax)

  ! resolution ratio factor (near-boundary res is fac * mid level)
  fac = 10; 

  ! nondimensional midlevel in *both* coordinates
  znot = .5; 

  ! grid stretch factor
  alpha = 2.*acosh(sqrt(fac));

  ! first do Jacobian factors for interior grid levels that have the 
  ! first level located one-half grid level off the boundary.
  do k = 1,pmax     
     z = (real(k) - .5)*dz     
     arg = alpha*(z-.5)
     facz(k) = 2.*tanh(alpha/2.)*(1./(1. - (tanh(arg))**2))/alpha
     zl(k) = znot + (tanh(alpha*(z-.5))) / (2.*tanh(alpha/2.))
     print*,'jacobian:',k,z,zl(k),facz(k)
  enddo

  !now do the "half levels", which are actually spaced evenly with no staggering.
  do k = 1,pmax-1
     z = real(k)*dz     
     arg = alpha*(z-.5)
     faczo(k) = 2.*tanh(alpha/2.)*(1./(1. - (tanh(arg))**2))/alpha
     zhl(k) = znot + (tanh(alpha*(z-.5))) / (2.*tanh(alpha/2.))
     print*,'jacobian:',k,z,zhl(k),faczo(k)
  enddo

  return
END SUBROUTINE setup_tanh

FUNCTION acosh(x)

  real, intent(in) :: x
  real :: acosh
  
  if (x .lt. 1) then 
     print*,'error in acosh routine...bad argument'
     return
  else
     acosh = log(x + sqrt(x*x - 1))
  endif
  return
     
END FUNCTION acosh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE epv_pseudoheight(u,v,t,ub,vb,tb,ut,vt,tt,out)

  implicit none

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

  integer, parameter :: nx = 2*kmax, ny = 2*lmax
  complex, dimension(nx,ny,pmax), intent(inout) :: u,v,t
  complex, dimension(nx,ny), intent(in) :: ub,vb,tb,ut,vt,tt
  real, dimension(nx,ny,pmax), intent(out) :: out
  !
  complex, dimension(nx,ny,pmax) :: uz,vz,tz
  complex, dimension(2*kmax,2*lmax) :: dx,dy,Id,vx,uy,tx,ty
  complex, dimension(nx,ny) :: tmp
  real :: dz
  integer :: i,j,k

  ! vertical grid distance
  dz = ZH/real(pmax); if (verbose .gt. 1) print*,'dz = ',dz

  ! derivative operators
  call d_setup(dx,dy,Id)

  ! zero mean theta
  !tb(1,1) = 0.; tt(1,1) = 0.; t(1,1,:) = 0.

  ! vertical derivatives by centered finite differences, except near boundaries.
  ! the near-boundary points use centered differencing and averaging to get the 
  ! mid-grid points (e.g. grid level 1.5).
  do k=1,pmax
     if (k .eq. 1) then 
        vz(:,:,k) = (v(:,:,k+1) + v(:,:,k) - 2.*vb(:,:)) / (2.*dz)
        uz(:,:,k) = (u(:,:,k+1) + u(:,:,k) - 2.*ub(:,:)) / (2.*dz)
        tz(:,:,k) = (t(:,:,k) - tb(:,:)) / dz
     elseif (k .eq. pmax) then 
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
  do k=1,pmax
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
END SUBROUTINE epv_pseudoheight
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE epv_pseudoheight_new(u,v,t,tb,tt,phi0,phi0b,phi0t,out)

  implicit none

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

  integer, parameter :: nx = 2*kmax, ny = 2*lmax
  complex, dimension(nx,ny,pmax), intent(inout) :: u,v,t,phi0
  complex, dimension(nx,ny), intent(inout) :: tb,tt,phi0b,phi0t
  real, dimension(nx,ny,pmax), intent(out) :: out
  !
  complex, dimension(nx,ny,pmax) :: uz,vz,tz,u0,v0
  complex, dimension(2*kmax,2*lmax) :: dx,dy,Id,vx,uy,tx,ty,v0x,u0y,u0b,u0t,v0b,v0t
  complex, dimension(nx,ny) :: tmp
  real :: dz
  integer :: i,j,k

  ! vertical grid distance
  dz = ZH/real(pmax); if (verbose .gt. 1) print*,'dz = ',dz

  ! derivative operators
  call d_setup(dx,dy,Id)

  ! zero mean theta
  tb(1,1) = 0.; tt(1,1) = 0.; t(1,1,:) = 0.

  ! leading order winds
  do k = 1,pmax
     call d_s2s(phi0(:,:,k),u0(:,:,k),1,-dy)
     call d_s2s(phi0(:,:,k),v0(:,:,k),1,dx)
  enddo
  
  ! leading-order boundary winds
  call d_s2s(phi0b,u0b,1,-dy); call d_s2s(phi0t,u0t,1,-dy)
  call d_s2s(phi0b,v0b,1,dx);  call d_s2s(phi0t,v0t,1,dx)
  
  ! vertical derivatives
  do k=1,pmax
     if (k .eq. 1) then 
        vz(:,:,k) = (v0(:,:,k+1) + v0(:,:,k) - 2.*v0b(:,:)) / (2.*dz)
        uz(:,:,k) = (u0(:,:,k+1) + u0(:,:,k) - 2.*u0b(:,:)) / (2.*dz)
        tz(:,:,k) = (t(:,:,k) - tb(:,:)) / dz
     elseif (k .eq. pmax) then 
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
  do k=1,pmax
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
END SUBROUTINE epv_pseudoheight_new
!========================================================================

!========================================================================
SUBROUTINE write_diag(output_file,it,u_s,v_s,theta_s,phi_s,u,v,theta,phi,w,pv)
! write to disk

    implicit none

    character(len=*),                    intent(in) :: output_file
    integer,                             intent(in) :: it
    real, dimension(2*kmax,2*lmax,2),    intent(in) :: u_s,v_s,theta_s,phi_s
    real, dimension(2*kmax,2*lmax,pmax), intent(in) :: u,v,theta,phi,w,pv

    integer :: ncid, vardim(5), svardim(4), pvardim(4), varid
    integer :: count(4), start(4), ierr
    real    :: time

    if ( it .eq. 0 ) then

        ! Create a new NetCDF file
        call nc_check( nf90_create(output_file, NF90_CLOBBER .or. NF90_64BIT_OFFSET, ncid), 'write_diag', 'create ' // trim(output_file) )

        ! Define dimensions
        call nc_check( nf90_def_dim(ncid, "nx",   2*kmax,         vardim(1)), 'write_diag', 'def_dim, nx '   // trim(output_file) )
        call nc_check( nf90_def_dim(ncid, "ny",   2*lmax,         vardim(2)), 'write_diag', 'def_dim, ny '   // trim(output_file) )
        call nc_check( nf90_def_dim(ncid, "nz",   2,              vardim(3)), 'write_diag', 'def_dim, nz '   // trim(output_file) )
        call nc_check( nf90_def_dim(ncid, "np",   pmax,           vardim(4)), 'write_diag', 'def_dim, np '   // trim(output_file) )
        call nc_check( nf90_def_dim(ncid, "time", NF90_UNLIMITED, vardim(5)), 'write_diag', 'def_dim, time ' // trim(output_file) )

        ! Define variable dimensions
        svardim(1) = vardim(1) ; pvardim(1) = vardim(1)
        svardim(2) = vardim(2) ; pvardim(2) = vardim(2)
        svardim(3) = vardim(3) ; pvardim(3) = vardim(4)
        svardim(4) = vardim(5) ; pvardim(4) = vardim(5)

        ! Define variables
        call nc_check( nf90_def_var(ncid, "time",    NF90_FLOAT, vardim(5), varid),  'write_diag', 'def_var, time '    // trim(output_file) )
        call nc_check( nf90_def_var(ncid, "u_s",     NF90_FLOAT, svardim,   varid),  'write_diag', 'def_var, u_s '     // trim(output_file) )
        call nc_check( nf90_put_att(ncid, varid, "description", "surface zonal wind"), 'write_diag', 'put_att, description ' // trim(output_file) )
        call nc_check( nf90_def_var(ncid, "v_s",     NF90_FLOAT, svardim,   varid),  'write_diag', 'def_var, v_s '     // trim(output_file) )
        call nc_check( nf90_put_att(ncid, varid, "description", "surface meridional wind"), 'write_diag', 'put_att, description ' // trim(output_file) )
        call nc_check( nf90_def_var(ncid, "theta_s", NF90_FLOAT, svardim,   varid),  'write_diag', 'def_var, theta_s ' // trim(output_file) )
        call nc_check( nf90_put_att(ncid, varid, "description", "surface potential temperature"), 'write_diag', 'put_att, description ' // trim(output_file) )
        call nc_check( nf90_def_var(ncid, "phi_s",   NF90_FLOAT, svardim,   varid),  'write_diag', 'def_var, phi_s '   // trim(output_file) )
        call nc_check( nf90_put_att(ncid, varid, "description", "surface geopotential"), 'write_diag', 'put_att, description ' // trim(output_file) )
        call nc_check( nf90_def_var(ncid, "u",       NF90_FLOAT, pvardim,   varid),  'write_diag', 'def_var, u '       // trim(output_file) )
        call nc_check( nf90_put_att(ncid, varid, "description", "interior zonal wind"), 'write_diag', 'put_att, description ' // trim(output_file) )
        call nc_check( nf90_def_var(ncid, "v",       NF90_FLOAT, pvardim,   varid),  'write_diag', 'def_var, v '       // trim(output_file) )
        call nc_check( nf90_put_att(ncid, varid, "description", "interior meridional wind"), 'write_diag', 'put_att, description ' // trim(output_file) )
        call nc_check( nf90_def_var(ncid, "theta",   NF90_FLOAT, pvardim,   varid),  'write_diag', 'def_var, theta '   // trim(output_file) )
        call nc_check( nf90_put_att(ncid, varid, "description", "interior potential temperature"), 'write_diag', 'put_att, description ' // trim(output_file) )
        call nc_check( nf90_def_var(ncid, "phi",     NF90_FLOAT, pvardim,   varid),  'write_diag', 'def_var, phi '     // trim(output_file) )
        call nc_check( nf90_put_att(ncid, varid, "description", "interior geopotential"), 'write_diag', 'put_att, description ' // trim(output_file) )
        call nc_check( nf90_def_var(ncid, "w",       NF90_FLOAT, pvardim,   varid),  'write_diag', 'def_var, w '       // trim(output_file) )
        call nc_check( nf90_put_att(ncid, varid, "description", "interior vertical velocity"), 'write_diag', 'put_att, description ' // trim(output_file) )
        call nc_check( nf90_def_var(ncid, "pv",      NF90_FLOAT, pvardim,   varid),  'write_diag', 'def_var, pv '      // trim(output_file) )
        call nc_check( nf90_put_att(ncid, varid, "description", "interior potential vorticity"), 'write_diag', 'put_att, description ' // trim(output_file) )

        ! Put global attributes
        call nc_check( nf90_put_att(ncid, NF90_GLOBAL, "model", model), 'write_diag', 'put_att, model ' // trim(output_file) )
        call nc_check( nf90_put_att(ncid, NF90_GLOBAL, "ntims", ntims), 'write_diag', 'put_att, ntims ' // trim(output_file) )
        call nc_check( nf90_put_att(ncid, NF90_GLOBAL, "dt",    dt),    'write_diag', 'put_att, dt '    // trim(output_file) )
        call nc_check( nf90_put_att(ncid, NF90_GLOBAL, "iplot", iplot), 'write_diag', 'put_att, iplot ' // trim(output_file) )
        call nc_check( nf90_put_att(ncid, NF90_GLOBAL, "XL",    XL),    'write_diag', 'put_att, XL '    // trim(output_file) )
        call nc_check( nf90_put_att(ncid, NF90_GLOBAL, "YL",    YL),    'write_diag', 'put_att, YL '    // trim(output_file) )
        call nc_check( nf90_put_att(ncid, NF90_GLOBAL, "H",     H),     'write_diag', 'put_att, H '     // trim(output_file) )
        call nc_check( nf90_put_att(ncid, NF90_GLOBAL, "Ross",  Ross),  'write_diag', 'put_att, Ross '  // trim(output_file) )
        call nc_check( nf90_put_att(ncid, NF90_GLOBAL, "gamma", gamma), 'write_diag', 'put_att, gamma ' // trim(output_file) )
        call nc_check( nf90_put_att(ncid, NF90_GLOBAL, "n",     n),     'write_diag', 'put_att, n '     // trim(output_file) )
        call nc_check( nf90_put_att(ncid, NF90_GLOBAL, "tau",   tau),   'write_diag', 'put_att, tau '   // trim(output_file) )
        call nc_check( nf90_put_att(ncid, NF90_GLOBAL, "trl",   trl),   'write_diag', 'put_att, trl '   // trim(output_file) )
        call nc_check( nf90_put_att(ncid, NF90_GLOBAL, "amu",   amu),   'write_diag', 'put_att, amu '   // trim(output_file) )
        call nc_check( nf90_put_att(ncid, NF90_GLOBAL, "shear", shear), 'write_diag', 'put_att, shear ' // trim(output_file) )

        call nc_check( nf90_enddef(ncid), 'write_diag', 'enddef, ' // trim(output_file) )
        call nc_check( nf90_close(ncid),  'write_diag', 'close, '  // trim(output_file) )

    else 

        time = (it-1)*dt

        ! Open the netCDF file
        call nc_check( nf90_open(output_file, NF90_WRITE, ncid), 'write_diag', 'open, '           // trim(output_file) )

        ! Write time variable
        call nc_check( nf90_inq_varid(ncid, "time",    varid) ,  'write_diag', 'inq_varid, time ' // trim(output_file) )

        count(1) = 2*kmax;  start(1) = 1
        count(2) = 2*lmax;  start(2) = 1
        count(3) = 2;       start(3) = 1
        count(4) = 1;       start(4) = it

        ! Write surface variables
        call nc_check( nf90_put_var(ncid, varid, time, (/it/)),          'write_diag', 'put_var, time '      // trim(output_file) )
        call nc_check( nf90_inq_varid(ncid, "u_s",     varid) ,          'write_diag', 'inq_varid, u_s '     // trim(output_file) )
        call nc_check( nf90_put_var(ncid, varid, u_s,     start, count), 'write_diag', 'put_var, u_s '       // trim(output_file) )
        call nc_check( nf90_inq_varid(ncid, "v_s",     varid) ,          'write_diag', 'inq_varid, v_s '     // trim(output_file) )
        call nc_check( nf90_put_var(ncid, varid, v_s,     start, count), 'write_diag', 'put_var, v_s '       // trim(output_file) )
        call nc_check( nf90_inq_varid(ncid, "theta_s", varid),           'write_diag', 'inq_varid, theta_s ' // trim(output_file) )
        call nc_check( nf90_put_var(ncid, varid, theta_s, start, count), 'write_diag', 'put_var, theta_s '   // trim(output_file) )
        call nc_check( nf90_inq_varid(ncid, "phi_s",   varid),           'write_diag', 'inq_varid, phi_s '   // trim(output_file) )
        call nc_check( nf90_put_var(ncid, varid, phi_s,   start, count), 'write_diag', 'put_var, phi_s '     // trim(output_file) )

        count(1) = 2*kmax;  start(1) = 1
        count(2) = 2*lmax;  start(2) = 1
        count(3) = pmax;    start(3) = 1
        count(4) = 1;       start(4) = it

        ! Write interior variables
        call nc_check( nf90_inq_varid(ncid, "u",     varid) ,          'write_diag', 'inq_varid, u '     // trim(output_file) )
        call nc_check( nf90_put_var(ncid, varid, u,     start, count), 'write_diag', 'put_var, u '       // trim(output_file) )
        call nc_check( nf90_inq_varid(ncid, "v",     varid) ,          'write_diag', 'inq_varid, v '     // trim(output_file) )
        call nc_check( nf90_put_var(ncid, varid, v,     start, count), 'write_diag', 'put_var, v '       // trim(output_file) )
        call nc_check( nf90_inq_varid(ncid, "theta", varid),           'write_diag', 'inq_varid, theta ' // trim(output_file) )
        call nc_check( nf90_put_var(ncid, varid, theta, start, count), 'write_diag', 'put_var, theta '   // trim(output_file) )
        call nc_check( nf90_inq_varid(ncid, "phi",   varid),           'write_diag', 'inq_varid, phi '   // trim(output_file) )
        call nc_check( nf90_put_var(ncid, varid, phi,   start, count), 'write_diag', 'put_var, phi '     // trim(output_file) )
        call nc_check( nf90_inq_varid(ncid, "w",     varid) ,          'write_diag', 'inq_varid, w '     // trim(output_file) )
        call nc_check( nf90_put_var(ncid, varid, w,     start, count), 'write_diag', 'put_var, w '       // trim(output_file) )
        call nc_check( nf90_inq_varid(ncid, "pv",    varid) ,          'write_diag', 'inq_varid, pv '    // trim(output_file) )
        call nc_check( nf90_put_var(ncid, varid, pv,    start, count), 'write_diag', 'put_var, pv '      // trim(output_file) )

        ! Close the netCDF file
        call nc_check( nf90_close(ncid),                               'write_diag', 'close, '           // trim(output_file) )

    endif

    return
END SUBROUTINE write_diag
!========================================================================

!========================================================================
SUBROUTINE scaling(grav,tnot,Rhos,Ns,km,Cor,Ls,Hs,Us,Ws,Ps,Ts,PVs)
! scaling parameters:

    implicit none

    real, intent(out) :: grav,tnot,Rhos,Ns,km,Cor,Ls,Hs,Us,Ws,Ps,Ts,PVs

! The _independent_ scaling parameters are:
! - Rossby number.
! - horizontal length scale.
! - vertical length scale.
! - surface temperature (tnot).
! - Coriolis parameter.
! - gravitational constant.
! - density.

! The _dependent_ scaling parameters are:
! - N ~ f * L / H
! - U ~ R * f * L
! - W ~ R * U * H / L
! - P ~ U * f * L
! - T ~ U * f * L * tnot / (g * H)
! - PV ~ U * f^2 * L * tnot / (rho * g * H^2)

    grav = 9.81        ! gravitational constant (m s^-2)
    tnot = 300.        ! theta constant at Z = 0 (K)
    km   = 1000.       ! 1000 m in a kilometer
    Cor  = 1.e-4       ! Coriolis parameter (s^-1)
    Hs   = 10.*km      ! vertical length scale (m)
    Ls   = 1000.*km    ! horizontal length scale (m)
    Rhos = 1.0         ! density (kg m^-3)
    Ns   = Cor*Ls/Hs   ! bouyancy frequency (s-^-1)

    if (Ross .ne. 0) then 
        Us = Ross*Cor*Ls   ! Horizontal wind scale (m s^-1)
        if (verbose .gt. 1)  print*,'Us = ',Us
    else
        Us = 1.0
    endif

    Ws  = Ross*Hs*Us/Ls                 ! Vertical wind scale (m s^-1)
    Ps  = Us*Cor*Ls                     ! geopotential scale (m^2 s^-2)
    Ts  = Ps*tnot/(grav*Hs)             ! potential temperature (K)
    PVs = 1.e6 * Cor * Ts / (Rhos * Hs) ! PV scaling (PVU)

    return
END SUBROUTINE scaling
!========================================================================

!========================================================================
SUBROUTINE dx_echo(dxs,dys,dzs,dts)
! send back dimensional values for dx, dy, dz, and dt.

    implicit none

    real, intent(out) :: dxs,dys,dzs,dts
    real :: grav,tnot,Rhos,Ns,km,Cor,Ls,Hs,Us,Ws,Ps,Ts,PVs

    call scaling(grav,tnot,Rhos,Ns,km,Cor,Ls,Hs,Us,Ws,Ps,Ts,PVs)

    dxs = Ls*XL/real(2*kmax) ! meters
    dys = Ls*YL/real(2*kmax) ! meters
    dzs = Hs*ZH/real(pmax)   ! meters
    dts = Ls*dt/Us           ! seconds

    return
END SUBROUTINE dx_echo
!========================================================================

!========================================================================
SUBROUTINE nc_check(ierr,subr_name,context)

  ! check for netcdf errors

  implicit none

  integer,          intent(in) :: ierr
  character(len=*), intent(in) :: subr_name, context

  character(len=129) :: error_msg

  if (ierr /= nf90_noerr) then
    error_msg = trim(subr_name) // ': ' // trim(context) // ': ' // trim(nf90_strerror(ierr))
    print*,trim(adjustl(error_msg))
    stop
  end if

  return
END SUBROUTINE nc_check
!========================================================================

!========================================================================
FUNCTION get_ntimes(infile)

    implicit none

    character(len=*), intent(in) :: infile
    integer :: get_ntimes
    integer :: ncid, dimid

    call nc_check( nf90_open(trim(infile), NF90_NOWRITE, ncid), 'get_ntimes', 'open, ' // trim(infile) )
    call nc_check( nf90_inq_dimid(ncid, 'time', dimid), 'get_ntimes', 'inq_dimid time, ' // trim(infile) )
    call nc_check( nf90_inquire_dimension(ncid, dimid, len = get_ntimes), 'get_ntimes', 'inquire_dimension time,' // trim(infile) )
    call nc_check( nf90_close(ncid), 'get_ntimes', 'close, ' // trim(infile) )

    return

END FUNCTION get_ntimes
!========================================================================

END MODULE pvinv_mod
!========================================================================
