!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! <next few lines under version control, D O  N O T  E D I T>
! $Date$
! $Author$
! $Revision$
! $Id$
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! v.3  6 September 2005
!
! a module to be used with spectral surface-based models.
!
MODULE spectral

  integer, parameter :: model = 1  !0: 2D; 1:2sQG; 2:sQG trop; 3:sQG sfc; 4:HsQG
  real,    parameter :: dt =  0.01        ! model time step
  integer, parameter :: ntims =  1001     ! number of model time steps
  integer, parameter :: iplot =  100     ! plots every iplot time steps
  integer, parameter :: imean =  int(10.e3/dt) ! compute area mean every imean

  logical, parameter :: iterr = .FALSE.   ! flag for terrain
  logical, parameter :: linear = .FALSE.  ! flag for linear advection
  logical, parameter :: trop = .FALSE.    ! flag for tropopause geometry

  real, parameter :: pi=3.14159265358979323846264338327950288419

  integer, parameter:: kmax = 128/2 ! number of x waves
  integer, parameter:: lmax = 64/2  ! number of y waves

  integer, parameter:: n = 8  ! diffusion parameter
  real, parameter :: tau = 20*dt ! diffusion time scale 

!  real, parameter :: XL = 2*pi   ! x domain length
!  real, parameter :: YL = 2*pi   ! y domain length
!  real, parameter :: XL = 20.0   ! x domain length
!  real, parameter :: XL = 30.0   ! x domain length
  real, parameter :: XL = 20.0   ! x domain length
  real, parameter :: YL = 2*5.539118  ! y domain length
  real, parameter :: H  = 1.0   ! z domain length
!  real, parameter :: H  = .1   ! z domain length
!  real, parameter :: H  = 0.01   ! z domain length
!  real, parameter :: H  = 10.   ! z domain length

!  real, parameter :: Ross = 0.1 ! Rossby number
  real, parameter :: Ross = 0. ! Rossby number

!  real, parameter :: gamma=0.0  ! ekman parameter
  real, parameter :: gamma=0.0  ! ekman parameter
!  real, parameter :: gamma=0.075*1.5  ! ekman parameter

!  logical, parameter :: hw = .FALSE.!  jet off
  logical, parameter :: hw = .TRUE.!  jet on
  logical, parameter :: restart = .FALSE. ! use restart file for ic
!  logical, parameter :: restart = .TRUE. ! use restart file for ic
  logical, parameter :: grow = .FALSE. ! grow a normal mode
  real, parameter :: amu = 1.0 ! HW jet parameter (0->1)
  real, parameter :: shear = 1.0 ! shear parameter (1 for HW jet)

!  real, parameter :: trl = 40. ! jet relaxation parameter (e-folding time)
  real, parameter :: trl = 15. ! jet relaxation parameter (e-folding time)
!  real, parameter :: trl = 1000 ! jet relaxation parameter (OFF)
  logical, parameter :: cdf = .TRUE.      ! write to netcdf files
!  logical, parameter :: cdf = .FALSE.     ! write to text ASCII files

  integer, parameter:: verbose = 2  ! (0: no mesgs; 1: imp only; 2: all)
!  character path*2; parameter (path='./') ! file location:
!  character path*7; parameter (path='../run/') ! file location:
  character path*2; parameter (path='./') ! file location:

  ! do not modify below this point:
  integer, parameter :: kmid=kmax/2,lmid=lmax/2,mmax=3.125*kmax, &
       &        nmax=3.125*lmax,l2=nmax-lmax,k2=mmax-kmax, &
       &        lmaxp1=lmax+1,kmaxp1=kmax+1
  real, parameter :: facx=2.*pi/XL,facy=2.*pi/YL,eps=0.1
  real, parameter :: amiss=-99.,hwp=5.539118,ryl=2.*pi/hwp

! START pvinv specific-------------------------------------------------------
  real, parameter :: ZH = H ! z domain length
  integer, parameter :: pmax = 21 ! number of vertical levels

!  integer, parameter :: order = 0 ! leading-order only
!  integer, parameter :: order = 1 ! first-order only
  integer, parameter :: order = 2 ! use this for full inversion

CONTAINS

  ! scaling parameters:
  SUBROUTINE scaling(grav,tnot,Ns,km,Cor,Ls,Hs,Us,Ws,Ps,Ts)

    real, intent(out) :: grav,tnot,Ns,km,Cor,Ls,Hs,Us,Ws,Ps,Ts

    grav = 9.81            ! gravitational constant (m s^-2)
    tnot = 300.            ! theta constant at Z = 0 (K)
    Ns = 1.e-2             ! bouyancy frequency
    km = 1000.             ! 1000 m in a kilometer
    Cor = 1.e-4            ! Coriolis parameter (s^-1)
    Hs = 10.*km            ! vertical length scale (m)
    Ls = Ns*Hs/Cor         ! horizontal length scale (m)

    if (Ross .ne. 0) then 
       Us = Ross*Cor*Ls       ! Horizontal wind scale (m s^-1)
       if (verbose .gt. 1)        print*,'Us = ',Us
    else
       Us = 1.
    endif

    Ws = Ross*Hs*Us/Ls     ! Vertical wind scale (m s^-1)
    Ps = Us*Cor*Ls         ! geopotential scale (m^2 s^-2)
    Ts = Ps*tnot/(grav*Hs) ! potential temperature (K)

    return

  end SUBROUTINE scaling

  ! rolv's routine to get dimensional dx, dy, dz, and dt
  SUBROUTINE dx_echo(dxs,dys,dzs,dts)
    ! send back dimensional values for dx, dy, dz, and dt.

    real, intent(out) :: dxs,dys,dzs,dts
    real :: grav,tnot,Ns,km,Cor,Ls,Hs,Us,Ws,Ps,Ts

    call scaling(grav,tnot,Ns,km,Cor,Ls,Hs,Us,Ws,Ps,Ts)

    dxs = Ls*XL/real(2*kmax) ! meters
    dys = Ls*YL/real(2*kmax) ! meters
    dzs = Hs*ZH/real(pmax)   ! meters
    dts = Ls*dt/Us           ! seconds

    return
  end SUBROUTINE dx_echo
  ! END pvinv specific-------------------------------------------------------
ENDMODULE spectral
