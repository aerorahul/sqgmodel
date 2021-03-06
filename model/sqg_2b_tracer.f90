!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! <next few lines under version control, D O  N O T  E D I T>
! $Date$
! $Author$
! $Revision$
! $Id$
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Fortran90 code for a uniform-PV two-surface QG+1 model in spectral form.
!
! sQG dynamics are recovered for: Ross=0.0, H-->\infty
!  2D dynamics are recovered for: Ross=0.0, H-->0; \theta->\theta*H

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM sqg_spectral

    use spectral

    call main

    stop

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE main

    implicit none

    ! spectral values
    complex, dimension(2*kmax,2*lmax) :: thspB,thbB,thspT,thbT
    complex, dimension(2*kmax,2*lmax) :: trspB,trspT
    complex, dimension(2*kmax,2*lmax) :: sbold,sB
    ! spectral work arrays
    real,    dimension(2*kmax,2*lmax) :: thxyB,thpxyB,thxyT,thpxyT
    real,    dimension(2*kmax,2*lmax) :: trxyB,trxyT
    ! grid point values
    real,    dimension(mmax,nmax) :: thxB,thyB,thxT,thyT
    real,    dimension(mmax,nmax) :: uB,vB,uT,vT
    real,    dimension(mmax,nmax) :: trxB,tryB,trxT,tryT
    ! tendencies
    real,    dimension(mmax,nmax) :: lap,hx,hy,hu,hv 
    real,    dimension(mmax,nmax) :: tthB,tthT
    real,    dimension(mmax,nmax) :: ttrB,ttrT
    complex, dimension(mmax,nmax) :: tthspB,tthspT
    complex, dimension(mmax,nmax) :: ttrspB,ttrspT
    ! basic state
    real,    dimension(mmax,nmax) :: ulinB,ulinT
    real,    dimension(mmax,nmax) :: thbyB,thbyT
    ! misc
    complex, dimension(2*kmax,2*lmax) :: Cblank
    real,    dimension(2*kmax,2*lmax) :: Rblanks
    real,    dimension(mmax,nmax)     :: Rblankb
    real                              :: dco,lam,time
    real                              :: cxB,cyB,cxT,cyT
    integer                           :: itime
    logical                           :: first,bot,top
    ! output
    character(len=64) :: smatfile,finalfile

    !!!!!!!!!!!!!!!!!!!!!!!!
    ! Start of executable: !
    !!!!!!!!!!!!!!!!!!!!!!!!

    Rblanks = 0.; Rblankb = 0.
    Cblank  = 0.

    if (verbose .gt. 0) print *,'Running with Rossby number: ', Ross

    ! initialize diffusion: 
    call diff(dco)

    ! initialize (read from Matlab file or restart file, option to add pert. to restart file):
    if (restart) then
        call init_restart(thxyB,thxyT)
        if (add_pert) then
            call init(thpxyB,thpxyT)
            thxyB = thxyB + thpxyB
            thxyT = thxyT + thpxyT
        endif
    else
        call init(thxyB,thxyT)
    endif

    ! initialize tracer:
    if (verbose .gt. 1) print*,'itracer = ',itracer
    if (itracer) call init_tracer(trxyB,trxyT)

    if (verbose .gt. 0) print *,'...writing to netCDF files'
    smatfile = trim(adjustl(path)) // 'smat.nc'
    call write_diag(smatfile,0,Rblanks,Rblanks,Rblanks,Rblanks)

    ! advection flags
    top = .TRUE.; bot = .TRUE.
    ! if (maxval(abs(thxyT)) .lt. 1.e-5) top = .FALSE.
    ! if (maxval(abs(thxyB)) .lt. 1.e-5) bot = .FALSE.
    if (model .eq. 2) then                       ! tropo sQG
        top = .TRUE.; bot = .FALSE.
    elseif (model .eq. 3 .or. model .eq. 0) then ! surface sQG
        top = .FALSE.; bot = .TRUE.
    elseif (model .eq. 4) then                   ! tropo HsQG
        top = .TRUE.; bot = .FALSE.
    endif

    if (verbose .gt. 0) print*,'max BOTTOM initial value=',maxval(abs(thxyB)),bot
    if (verbose .gt. 0) print*,'max TOPPOM initial value=',maxval(abs(thxyT)),top

    ! map into spectral space at the same resolution:
    call xy_to_sp(cmplx(thxyB,0.),thspB,2*kmax,2*lmax,kmax,lmax)
    call xy_to_sp(cmplx(thxyT,0.),thspT,2*kmax,2*lmax,kmax,lmax)
    if (itracer) then
        call xy_to_sp(cmplx(trxyB,0.),trspB,2*kmax,2*lmax,kmax,lmax)
        call xy_to_sp(cmplx(trxyT,0.),trspT,2*kmax,2*lmax,kmax,lmax)
    endif

    ! option to zero k = 1 modes
    !thspB(2,:) = 0.; thspB(2*kmax,:) = 0.
    !thspT(2,:) = 0.; thspT(2*kmax,:) = 0.

    ! init jet
    if (hw) then

        call init_jet(thbB,thbT,thbyB,thbyT,ulinB,ulinT,lam)
        ulinB = ulinB + 0; ulinT = ulinT - 0*H*lam     ! barotropic wind (Ross = 0!)

        ! write basic state to disk (linear shear:=  ON : .TRUE., OFF : .FALSE.)
        call write_diag(trim(adjustl(path)) // 'base.nc',0,Rblanks,Rblanks,Rblanks,Rblanks)
        call dump(thbB,thbT,Cblank,Cblank,.TRUE.,lam,1,trim(adjustl(path)) // 'base.nc')

    else

        thbB = 0;thbT = 0;thbyB=0.;thbyT=0.;ulinB=0.;ulinT=0.;lam=0.

    endif

    first = .TRUE.
    if (verbose .gt. 1) print*,'lam = ',lam
    if (verbose .gt. 1) print*,'extrema ulinT = ',maxval(ulinT),minval(ulinT)

    open(21,file = trim(adjustl(path)) // 'mean.dat',status='replace')
    write(21,*) 'Mean theta, time:'; close(21)

    ! initialize terrain:
    hu = 0.; hv = 0.; hx = 0.; hy = 0.
    if (verbose .gt. 1)   print*,'iterr = ',iterr
    if (iterr .and. Ross .eq. 0) then 
        call terrain(hx,hy,hu,hv)
        call invert(thspB,0*thspT,thxB,thxT,thyB,thyT,vB,vT,uB,uT,thbB,thbT, &
          &         thbyB,thbyT,ulinB,ulinT,first, .TRUE.,.TRUE.,0*lam,sB,sbold,lap)
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
        !thspB(:,1) = 0; thspT(:,1) = 0;

        time = real(itime-1)*dt
        if (mod((itime-1),imean) .eq. 0) then
            open(21,file=trim(adjustl(path)) // 'mean.dat',position='append')
            if (verbose .gt. 0) write(6,*) '...writing to mean.dat'
            write(21,*)  real(thspB(1,1)),real(thspT(1,1))
            write(21,*) time
            close(21)
        endif

        ! this is for growing the most unstable mode
        if (grow .and. cyT .gt. .001) then ! rescale 
            if (verbose .gt. 1) print*,'RESCALING'
            thspB = thspB/2.; thspT = thspT/2.
            cyT = 0.
        endif

        ! save old streamFUNCTION for Ekman calculation.
        sbold = sB; if (first) sbold = 0.

        ! Invert theta for streamFUNCTION; compute gradients for advection:
        call invert(thspB,thspT,thxB,thxT,thyB,thyT,vB,vT,uB,uT,thbB,thbT,thbyB,thbyT,ulinB,ulinT,first,bot,top,lam,sB,sbold,lap)

        ! Invert tracer to compute gradients for advection:
        if (itracer) then
        call invert_tracer(trspB,trspT,trxB,trxT,tryB,tryT,first,bot,top)
        endif
 
        ! option to compute potential enstrophy norm and growth rate:
        if (inorm) call norm(thspB,thspT,itime)

        ! write data to da file for plots (last entry for linear field +):
        if (mod((itime-1),iplot) .eq. 0) then
            if (hw) then ! add basic state to output
                call dump(thspB+thbB,thspT+thbT,trspB,trspT,.TRUE.,lam,itime,smatfile)
            else ! no basic state
                call dump(thspB,thspT,trspB,trspT,.FALSE.,lam,itime,smatfile)
            endif
        endif

        ! spectral advection:
        if (bot) call advect(uB,vB,thxB,thyB,thbyB,hx,hy,ulinB,tthB,lam,lap)
        if (top) call advect(uT+hu,vT+hv,thxT,thyT,thbyT,Rblankb,Rblankb,ulinT,tthT,lam,Rblankb)
        if (itracer) then
            if (bot) call advect(uB,vB,trxB,tryB,Rblankb,hx,hy,ulinB,ttrB,0.0,lap)
            if (top) call advect(uT+hu,vT+hv,trxT,tryT,Rblankb,Rblankb,Rblankb,ulinT,ttrT,0.0,Rblankb)
        endif

        ! Write Courant numbers to stdout
        if (mod((itime-1),10) .eq. 0) then 
            cxB = maxval(abs(uB+ulinB))*dt/(XL/real(2*kmax))
            cyB = maxval(abs(vB))*dt/(YL/real(2*lmax))
            cxT = maxval(abs(uT+ulinT))*dt/(XL/real(2*kmax))
            cyT = maxval(abs(vT))*dt/(YL/real(2*lmax))
            write(*,'(A23,F10.3,4F8.3)') 'time,cxB,cyB,cxT,cyT = ', &
            &            real(itime-1)*dt,cxB,cyB,cxT,cyT
        endif

        ! FFT back to spectral space:
        if (bot) then 
            tthspB = cmplx(tthB,0.); call ft_2d(tthspB,mmax,nmax,-1)
            if (itracer) then
                ttrspB = cmplx(ttrB,0.); call ft_2d(ttrspB,mmax,nmax,-1)
            endif
        endif
        if (top) then 
            tthspT = cmplx(tthT,0.); call ft_2d(tthspT,mmax,nmax,-1)
            if (itracer) then
                ttrspT = cmplx(ttrT,0.); call ft_2d(ttrspT,mmax,nmax,-1)
            endif
        endif

        ! advance one time step with explicit (hyper-)diffusion:
        if (bot) call thadvB(thspB,tthspB,dco,first)
        if (top) call thadvT(thspT,tthspT,dco,first)
        if (itracer) then
            if (bot) call tradvB(trspB,ttrspB,dco,first)
            if (top) call tradvT(trspT,ttrspT,dco,first)
        endif

        thspB(kmax,:) = 0.; thspB(lmax,:) = 0.
        if (itracer) then
            trspB(kmax,:) = 0.; trspB(lmax,:) = 0.
        endif
    
        first = .FALSE.

    end do ! itime
    !!!!!!!!!!!!!!!!!!
    ! END: Time loop !
    !!!!!!!!!!!!!!!!!!

    ! write (final + 1) time to disk for future restart:
    finalfile = trim(adjustl(path)) // 'final+1.nc'
    call write_diag(finalfile,0,Rblanks,Rblanks,Rblanks,Rblanks)
    call dump(thspB,thspT,Cblank,Cblank,.FALSE.,lam,1,finalfile)

END SUBROUTINE main
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE invert(ithspB,ithspT,thxBr,thxTr,thyBr,thyTr,vBr,vTr,uBr,uTr, &
     &    thbB,thbT,thbyB,thbyT,ulinB,ulinT,first,bot,top,lam,sB,sold,sblre)

  ! Invert PV and transform to a larger grid for de-aliasing.

  implicit none

  complex, dimension(2*kmax,2*lmax), intent(in)  :: ithspB,ithspT
  complex, dimension(2*kmax,2*lmax), intent(in)  :: thbB,thbT
  complex, dimension(2*kmax,2*lmax), intent(in)  :: sold
  real,    dimension(mmax,nmax),     intent(in)  :: thbyB,thbyT
  real,    dimension(mmax,nmax),     intent(in)  :: ulinB,ulinT
  real,                              intent(in)  :: lam
  logical,                           intent(in)  :: first,bot,top
  complex, dimension(2*kmax,2*lmax), intent(out) :: sB
  real,    dimension(mmax,nmax),     intent(out) :: thxBr,thxTr,thyBr,thyTr
  real,    dimension(mmax,nmax),     intent(out) :: uBr,uTr,vBr,vTr
  real,    dimension(mmax,nmax),     intent(out) :: sblre

  complex, dimension(mmax,nmax) :: thxB,thxT,thyB,thyT
  complex, dimension(mmax,nmax) :: uB,uT,vB,vT
  complex, dimension(mmax,nmax) :: copy

  complex, dimension(2*kmax,2*lmax) :: thspB,thspT
  complex, dimension(mmax,nmax)     :: szspB,szspT,szzspB,szzspT
  complex, dimension(mmax,nmax)     :: u1B,v1B,u1T,v1T

  complex, dimension(2*kmax,2*lmax) :: temps
  complex, dimension(2*kmax,2*lmax) :: u1spB,u1spT,v1spB,v1spT
  complex, dimension(2*kmax,2*lmax) :: tempspB,tempspT,tempxyB,tempxyT
  complex, dimension(2*kmax,2*lmax) :: htempB,htempT
  complex, dimension(mmax,nmax)     :: temp

  real, dimension(mmax,nmax) :: thx,thy
  real, dimension(mmax,nmax) :: u,v

  complex, dimension(2*kmax,2*lmax), save :: dx,dy,dz,dzo,iz,izo,dzi,Id
  integer,                           save :: pf,gf,bf
  logical,                           save :: correct

  !integer :: i,j,k,l

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

  thspB = ithspB; thspT = ithspT ! local copies of spectral boundary theta
!!!!!!!!!!!!!!!!!!!!!!!!!!! Grad Theta !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  thxT = 0.; thyT = 0.; vT = 0.; uT = 0.
  thxB = 0.; thyB = 0.; vB = 0.; uB = 0.

  call d_s2b(thspB,thxB,1,dx) ! x derivative: small to big domain.
  call d_s2b(thspB,thyB,1,dy) ! y derivative: small to big domain.
  call d_s2b(thspT,thxT,1,dx) ! x derivative: small to big domain.
  call d_s2b(thspT,thyT,1,dy) ! y derivative: small to big domain.

!!!!!!!!!!!!!!!!!!!!!!!!!!! sQG Bottom !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (bot .or. correct) then ! velocity contribution from lower boundary
     call d_b2b(thxB,vB,1,-iz) 
     call d_b2b(-thyB,uB,1,-iz)
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!! sQG Toppom !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (top .or. correct) then ! velocity contribution from upper boundary
     call d_b2b(thxT,vT,1,iz) ! z integral: big to big domain.
     call d_b2b(-thyT,uT,1,iz) ! z integral: big to big domain.
  endif
!!!!!!!!!!!!!!!!!!!!!!!! sQG Cross Boundary !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (bot .and. top .or. correct) then! lower velocity from upper boundary
     call d_b2b(thxT,temp,1,izo);  vB = vB + temp;
     call d_b2b(-thyT,temp,1,izo); uB = uB + temp;
  endif
  if (bot .and. top .or. correct) then! upper velocity from lower boundary
     call d_b2b(thxB,temp,1,-izo);  vT = vT + temp;
     call d_b2b(-thyB,temp,1,-izo); uT = uT + temp;
  endif
  ! FFT back to xy grid.
  if (bot .or. correct) then 
     call ft_2d(thxB,mmax,nmax,1) ! t_x
     call ft_2d(vB,mmax,nmax,1) ! v0
     call ft_2d(uB,mmax,nmax,1) ! u0
     call ft_2d(thyB,mmax,nmax,1) ! t_y
  endif
  if (top .or. correct) then 
     call ft_2d(thxT,mmax,nmax,1) ! t_x
     call ft_2d(thyT,mmax,nmax,1) ! t_y
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
        thspB = thspB + thbB; thspT = thspT + thbT
        uB = uB + ulinB; uT = uT + ulinT - (lam*H) ! don't double count!
        thyB = thyB + thbyB; thyT = thyT + thbyT
     endif

     ! big-grid theta (szsp) and theta_z (szzsp):
     call d_s2b(thspB,szspB,1,Id);  call ft_2d(szspB,mmax,nmax,1)
     call d_s2b(thspT,szspT,1,Id);  call ft_2d(szspT,mmax,nmax,1)

     call d_s2b(thspB,szzspB,1,-dz); call d_s2b(thspT,temp,1,dzo); 
     szzspB = szzspB + temp
     call ft_2d(szzspB,mmax,nmax,1)
     call d_s2b(thspT,szzspT,1,dz); call d_s2b(thspB,temp,1,-dzo); 
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
     temp = (-2.*thyB*szspB) + (uB*szzspB)
     call xy_to_sp(temp,tempspB,mmax,nmax,kmax,lmax)
     u1spB = u1spB + tempspB !add (-P1_y - F1_z)
     ! u1 toppom
     temp = (-2.*thyT*szspT) + (uT*szzspT)
     call xy_to_sp(temp,tempspT,mmax,nmax,kmax,lmax)
     u1spT = u1spT + tempspT !add (-P1_y - F1_z)
     ! v1 bottom
     temp = (2.*thxB*szspB) + (vB*szzspB)
     call xy_to_sp(temp,tempspB,mmax,nmax,kmax,lmax)
     v1spB = v1spB + tempspB !add (P1_x - G1_z)
     ! v1 toppom
     temp = (2.*thxT*szspT) + (vT*szzspT)
     call xy_to_sp(temp,tempspT,mmax,nmax,kmax,lmax)
     v1spT = v1spT + tempspT !add (P1_x - G1_z)

     ! linear-shear contributions (note: recycling var names)
     if (hw) then 
        call d_s2s(thspB,tempspB,1,-izo*dx*dx) !phi_xx^B
        call d_s2s(thspT,tempspT,1,izo*dx*dx) !phi_xx^T
        call d_s2s(thspB,htempB,1,-izo*dy*dy)  !phi_yy^B
        call d_s2s(thspT,htempT,1,izo*dy*dy)  !phi_yy^T
        call d_s2s(thspB,tempxyB,1,-izo*dx*dy) !phi_xy^B
        call d_s2s(thspT,tempxyT,1,izo*dx*dy) !phi_xy^T

        !         u1spB = u1spB - lam*(H*(htempT - tempspT) + (thspT - thspB) )
        !         u1spT = u1spT + lam*(H*(htempB - tempspB) + (thspT - thspB) )
        ! as per DM email 12/24/02:
        u1spB = u1spB - lam*(H*(htempT - tempspT) - thspB )
        u1spT = u1spT + lam*(H*(htempB - tempspB) + thspT )
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
     thyB = thyB - thbyB; thyT = thyT - thbyT
  endif ! correct

  ! return u = u_0 + Ro * u_1; v = v_0 + Ro * v_1
  ! (look into using the f90 "transfer" FUNCTION to return real parts...)

  thxBr = real(thxB); thyBr = real(thyB);
  thxTr = real(thxT); thyTr = real(thyT);      
  uBr = real(uB + Ross*u1B); uTr = real(uT + Ross*u1T)
  vBr = real(vB + Ross*v1B); vTr = real(vT + Ross*v1T)

  ! Ekman layer calculations (need streamFUNCTION and Laplacian.
  if (gamma .gt. 0) then
     ! reset these since they were changed for next-order calculations
     thspB = ithspB; thspT = ithspT ! local copies of spectral boundary theta
     sb = 0. ! surface O(1) streamFUNCTION
     call d_s2s(thspB,sB,1,-iz) ! bottom contribution
     call d_s2s(thspT,temps,1,izo) ! toppom contribution
     sB = sB + temps;
     ! now compute Laplacian from previous time-step's sB
     temp = 0.
     call d_s2b(sold,temp,1,dx*dx + dy*dy)
     call ft_2d(temp,mmax,nmax,1)
     sblre = real(temp)
  endif

  return
END SUBROUTINE invert
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE invert_tracer(itrspB,itrspT,trxBr,trxTr,tryBr,tryTr,first,bot,top)

  ! Invert tracer feild and compute gradients for advection.

  implicit none

  complex,  dimension(2*kmax,2*lmax), intent(in)  :: itrspB,itrspT
  logical,                            intent(in)  :: first,bot,top
  real,     dimension(mmax,nmax),     intent(out) :: trxBr,trxTr,tryBr,tryTr

  complex, dimension(mmax,nmax)     :: trxB,trxT,tryB,tryT
  complex, dimension(2*kmax,2*lmax) :: trspB,trspT
  
  complex, dimension(2*kmax,2*lmax), save :: dx1,dy1,dz1,dzo1,iz1,izo1,dzi1,Id1
  integer,                           save :: pf1,gf1,bf1
  logical,                           save :: correct
  
!!!!!!!!!!!!!!!!! Set-up on first time step !!!!!!!!!!!!!!!!!!!!!!!!!!
  if (first) then 
     print*,'calling setup'
     call d_setup(dx1,dy1,dz1,dzo1,iz1,izo1,Id1)    ! derivative operators
     !         print*,maxval(abs(dx1))
     !         print*,maxval(abs(dy1))
     !         print*,maxval(abs(dz1))
     !         print*,maxval(abs(iz1))
     !         print*,maxval(abs(dzo1))
     !         print*,maxval(abs(izo1))
     correct = .FALSE. ! flag for next-order corrections
     if (Ross .gt. 1.e-5) correct = .TRUE. 
     ! switches for removing rotation or divergence (sQG only)
     !         pf1 = 1; gf1 = 1; bf1 = 0  ! control: full sqg+1 corrections
     !         pf1 = 0; gf1 = 1; bf1 = 1  ! no rotational correction
     !         pf1 = 1; gf1 = 0; bf1 = -1 ! no divergent correction
  endif

  trspB = itrspB; trspT = itrspT ! local copies of spectral boundary tracer
  !!!!!!!!!!!!!!!!!!!!!!!!!!! Grad Tracer !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  trxT = 0.; tryT = 0.
  trxB = 0.; tryB = 0.

  call d_s2b(trspB,trxB,1,dx1) ! x derivative: small to big domain.
  call d_s2b(trspB,tryB,1,dy1) ! y derivative: small to big domain.
  call d_s2b(trspT,trxT,1,dx1) ! x derivative: small to big domain.
  call d_s2b(trspT,tryT,1,dy1) ! y derivative: small to big domain.
  
  if (bot .or. correct) then 
     call ft_2d(trxB,mmax,nmax,1) ! t_x
     call ft_2d(tryB,mmax,nmax,1) ! t_y
  endif
  if (top .or. correct) then 
     call ft_2d(trxT,mmax,nmax,1) ! t_x
     call ft_2d(tryT,mmax,nmax,1) ! t_y
  endif

  trxBr = real(trxB); tryBr = real(tryB);
  trxTr = real(trxT); tryTr = real(tryT);      

  return 
END SUBROUTINE invert_tracer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE xy_to_sp(xy,sp,mx,ny,km,lm)

  ! Map an (x,y) array onto a _smaller_ spectral array.
  ! Input: xy(mx,ny) --- a grid point array.
  ! Output: sp(2*km,2*lm) --- a spectral array.

  implicit none

  integer,                       intent(in)  :: mx,ny,km,lm
  complex, dimension(mx,ny),     intent(in)  :: xy
  complex, dimension(2*km,2*lm), intent(out) :: sp

  complex, dimension(mx,ny) :: copy
  integer                   :: kmp1,lmp1,k,l,k2,l2,kk,ll

  ! initialize arrays:
  sp = 0.; copy = xy

  kmp1 = km + 1; lmp1 = lm + 1; k2 = mx - km; l2 = ny - lm

  call ft_2d(copy,mx,ny,-1)

  !      do k=1,kmp1; do l=1,lmp1
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE sp_to_xy(sp,xy,km,lm,mx,ny)

  ! Map an (km,lm) spectral array onto a _bigger_ grid point array.
  ! Input: sp(2*km,2*lm) --- a spectral array.
  ! Output: xy(mx,ny) --- a grid point array.

  implicit none

  integer,                       intent(in)  :: km,lm,mx,ny
  complex, dimension(2*km,2*lm), intent(in)  :: sp
  real,    dimension(mx,ny),     intent(out) :: xy

  complex, dimension(mx,ny) :: copy
  integer                   :: kmp1,lmp1,k,l,k2,l2,kk,ll

  ! initialize arrays:
  copy = 0.

  kmp1 = km + 1; lmp1 = lm + 1; k2 = mx - km; l2 = ny - lm

  !      do 10 k = 1,kmp1; do 10 l = 1,lmp1
  do k = 1,km
     do l = 1,lm

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
     enddo !k
  enddo !l

  call ft_2d(copy,mx,ny,1)
  xy = real(copy)

  return
END SUBROUTINE sp_to_xy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE diff(dco)

  ! Compute the diffusion coefficient for del^n diffusion (n even>0).
  ! Tau gives the e-folding time scale for damping modes at the 
  ! Nyquist frequency.

  implicit none

  real, intent(out) :: dco
  real              :: temp

  ! 10/21/99 gjh mods...
  !      temp = (cmplx(0.,1.)**n)
  temp = 1
  !      dco = (((real(facx*kmax)**n) + (real(facy*lmax)**n))*tau)**(-1)
  dco = 1.0/(((real(facx*kmax)**2) &
       &          + (real(facy*lmax)**2))**(n/2)*tau)
  dco = temp*dco
  if (dco .lt. 0.) dco = -1.*dco

  if (verbose .gt. 1) print*,'diffusion coefficient:',dco

  return
END SUBROUTINE diff
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE init(thB,thT)

  ! Initialize basic state and perturbation fields.

  implicit none

  real, dimension(2*kmax,2*lmax), intent(out) :: thB,thT

  real    :: fac
  integer :: twokmax,twolmax
  integer :: ncid, dimid, varid

  if (verbose .gt. 0)  print*,'initializing theta fields...'
  thB = 0.; thT = 0.

  call nc_check( nf90_open(trim(adjustl(path)) // 'th_init.nc', NF90_NOWRITE, ncid) )
  call nc_check( nf90_inq_dimid(ncid, 'nx', dimid) )
  call nc_check( nf90_inquire_dimension(ncid, dimid, len = twokmax) )
  call nc_check( nf90_inq_dimid(ncid, 'ny', dimid) )
  call nc_check( nf90_inquire_dimension(ncid, dimid, len = twolmax) )
  call nc_check( nf90_inq_varid(ncid, 'thetaB', varid) )
  call nc_check( nf90_get_var(ncid, varid, thB) )
  call nc_check( nf90_inq_varid(ncid, 'thetaT', varid) )
  call nc_check( nf90_get_var(ncid, varid, thT) )
  call nc_check( nf90_close(ncid) )

  if (twokmax .ne. 2*kmax .and. twolmax .ne. 2*lmax) then  
     print*,'x and y resolution from th_init.nc do not match SPECTRAL.f90!!!'
     print*,'x and y from SPECTRAL.f90 : ', 2*kmax, 2*lmax
     print*,'x and y from th_init.nc   : ', twokmax,twolmax
     stop
  endif

  ! normalize initial amplitude (RMS00)
  !      fac = .15/maxval(abs(thB)); thB = fac*thB; thT = fac*thT
  !      thB = .15*thB/maxval(abs(thB)); thT = .15*thT/maxval(abs(thB))

  !      thT = -thT
  !      thB = -thB
  !      thB = 0.
  !      thT = 0.
  !      thT = thT/1000.
  !      thT = thT/10.
  !      thB = thB/10.

  return
END SUBROUTINE init
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE init_restart(thB,thT)

  ! Initialize theta fields (from restart file).

  implicit none

  real, dimension(2*kmax,2*lmax), intent(out) :: thB,thT

  integer :: twokmax, twolmax
  integer :: ncid, dimid, varid

  if (verbose .gt. 0) print*,'!!! USING RESTART FILE !!!'
  thB = 0.; thT = 0.

  ! load from a restart file (will need to save told, told2, eventually)
  call nc_check( nf90_open(trim(adjustl(path)) // 'restart.nc', NF90_NOWRITE, ncid) )
  call nc_check( nf90_inq_dimid(ncid, 'nx', dimid) )
  call nc_check( nf90_inquire_dimension(ncid, dimid, len = twokmax) )
  call nc_check( nf90_inq_dimid(ncid, 'ny', dimid) )
  call nc_check( nf90_inquire_dimension(ncid, dimid, len = twolmax) )
  call nc_check( nf90_inq_varid(ncid, 'thetaB', varid) )
  call nc_check( nf90_get_var(ncid, varid, thB) )
  call nc_check( nf90_inq_varid(ncid, 'thetaT', varid) )
  call nc_check( nf90_get_var(ncid, varid, thT) )
  call nc_check (nf90_close(ncid) )

  if (twokmax .ne. 2*kmax .and. twolmax .ne. 2*lmax) then  
     print*,'x and y resolution from restart.nc do not match SPECTRAL.f90!!!'
     print*,'x and y from SPECTRAL.f90 : ', 2*kmax, 2*lmax
     print*,'x and y from restart.nc   : ', twokmax,twolmax
     stop
  endif

  return
END SUBROUTINE init_restart
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE init_tracer(trB,trT)

  ! Initialize tracer fields

  implicit none

  real, dimension(2*kmax,2*lmax), intent(out) :: trB,trT

  integer :: twokmax, twolmax
  integer :: ncid, dimid, varid

  if (verbose .gt. 0)  print*,'initializing tracer fields...'
  trB = 0.; trT = 0.

  call nc_check( nf90_open(trim(adjustl(path)) //  'tr_init.nc', NF90_NOWRITE, ncid) )
  call nc_check( nf90_inq_dimid(ncid, 'nx', dimid) )
  call nc_check( nf90_inquire_dimension(ncid, dimid, len = twokmax) )
  call nc_check( nf90_inq_dimid(ncid, 'ny', dimid) )
  call nc_check( nf90_inquire_dimension(ncid, dimid, len = twolmax) )
  call nc_check( nf90_inq_varid(ncid, 'tracerB', varid) )
  call nc_check( nf90_get_var(ncid, varid, trB) )
  call nc_check( nf90_inq_varid(ncid, 'tracerT', varid) )
  call nc_check( nf90_get_var(ncid, varid, trT) )
  call nc_check (nf90_close(ncid) )

  if (twokmax .ne. 2*kmax .and. twolmax .ne. 2*lmax) then  
     print*,'x and y resolution from tr_init.nc do not match SPECTRAL.f90!!!'
     print*,'x and y from SPECTRAL.f90 : ', 2*kmax, 2*lmax
     print*,'x and y from tr_init.nc   : ', twokmax,twolmax
     stop
  endif

  return
END SUBROUTINE init_tracer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE init_jet(thbB,thbT,thbyB,thbyT,ulinB,ulinT,lam)

  ! Initialize basic state jet.
  ! thbB = periodic basic state theta; bottom boundary (spectral grid).
  ! thbT = periodic basic state theta; top boundary (spectral grid).
  ! thbyB = y derivative of periodic basic state theta; bottom boundary (grid point).
  ! thbyT = y derivative of periodic basic state theta; top boundary (grid point).
  ! ulinB = basic state wind; bottom boundary (grid point).
  ! ulinT = basic state wind; top boundary (grid point).

  implicit none

  complex, dimension(2*kmax,2*lmax), intent(out) :: thbB,thbT
  real,    dimension(mmax,nmax),     intent(out) :: thbyB,thbyT,ulinB,ulinT
  real,                              intent(out) :: lam

  complex, dimension(mmax,nmax)     :: thyB,thyT,uB,uT,temp
  complex, dimension(2*kmax,2*lmax) :: dx,dy,dz,dzo,iz,izo,dzi,Id
  real,    dimension(2*kmax,2*lmax) :: thbxyB,thbxyT
  real                              :: y,yp,dyy
  integer                           :: j

  if (verbose .gt. 0) print*,'initializing basic state...'

  ! first determine lamda, the linear shear parameter
  lam = (HW_theta(0.,0.) - HW_theta(YL,0.)) / YL
  if (verbose .gt. 1) print*,'lam = ',lam

  ! spectral variables
  dyy = YL/real(2*lmax)
  do j=1, 2*lmax
     y = real(j-1)*dyy; yp = y - (0.5*(YL - hwp))
     thbxyB(:,j) = HW_theta(y,0.) + lam*yp
     thbxyT(:,j) = HW_theta(y,H) + lam*yp
     !print*,'periodic theta grid point:',y,yp,thbxyT(1,j)
  enddo

  ! map into spectral space at the same resolution:
  call xy_to_sp(cmplx(thbxyB,0.),thbB,2*kmax,2*lmax,kmax,lmax)
  call xy_to_sp(cmplx(thbxyT,0.),thbT,2*kmax,2*lmax,kmax,lmax)

  ! grid point variables
  dyy = YL/real(nmax)
  do j=1, nmax
     y = real(j-1)*dyy; yp = y - (0.5*(YL - hwp))
     thbyB(:,j) = HW_thetay(y,0.) + lam
     thbyT(:,j) = HW_thetay(y,H) + lam
     ! old U
     !ulinB(:,j) = HW_ubar(y,0.)
     !ulinT(:,j) = HW_ubar(y,H)
  enddo

  ! new U: solve numerically given theta
  uB = 0.; uT = 0.
  call d_setup(dx,dy,dz,dzo,iz,izo,Id) ! derivative operators
  call d_s2b(thbB,thyB,1,dy); call d_s2b(thbT,thyT,1,dy)
  call d_b2b(-thyB,uB,1,-iz); call d_b2b(-thyT,uT,1,iz)
  call d_b2b(-thyT,temp,1,izo); uB = uB + temp;
  call d_b2b(-thyB,temp,1,-izo); uT = uT + temp;
  call ft_2d(uB,mmax,nmax,1); call ft_2d(uT,mmax,nmax,1)
  ulinB = real(uB); ulinT = real(uT)
  ulinT = ulinT + lam*H

  return
END SUBROUTINE init_jet
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION ran1(idum)

! Numerical Recipes random number generator:

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE dump(thspB,thspT,trspB,trspT,ilam,lam,it,outfile)

  !     Write theta or streamFUNCTION to disk.

  implicit none

  complex, dimension(2*kmax,2*lmax), intent(in) :: thspB,thspT
  complex, dimension(2*kmax,2*lmax), intent(in) :: trspB,trspT
  logical,                           intent(in) :: ilam
  real,                              intent(in) :: lam
  integer,                           intent(in) :: it
  character(len=*),                  intent(in) :: outfile

  complex, dimension(2*kmax,2*lmax) :: copy
  complex, dimension(2*kmax,2*lmax) :: thbB,thbT
  real,    dimension(2*kmax,2*lmax) :: thxyB,thxyT
  real,    dimension(2*kmax,2*lmax) :: trxyB,trxyT
  real                              :: y,dy
  integer                           :: j

  copy = thspB
  call sp_to_xy(copy,thxyB,kmax,lmax,2*kmax,2*lmax)
  copy = thspT
  call sp_to_xy(copy,thxyT,kmax,lmax,2*kmax,2*lmax)
  copy = trspB
  call sp_to_xy(copy,trxyB,kmax,lmax,2*kmax,2*lmax)
  copy = trspT
  call sp_to_xy(copy,trxyT,kmax,lmax,2*kmax,2*lmax)

  ! add in linear shear
  if ( ilam ) then 
     dy = YL/real(2*lmax)
     do j=1, 2*lmax
        ! fixed 08/10/2005 GJH & RBM
        y = real(j-1)*dy; 
        thxyB(:,j) = thxyB(:,j) - lam*y
        thxyT(:,j) = thxyT(:,j) - lam*y
     enddo
  endif

  if (verbose .gt. 0) print*,'Writing to disk...'
  call write_diag(outfile,it,thxyB,thxyT,trxyB,trxyT)
  
  ! write the init files if this is a grow run
  if (grow) then 
     open(31,file=trim(adjustl(path)) // 'th_init_B.dat') 
     write(31,'(i5)') 2*kmax
     write(31,'(i5)') 2*lmax
     write(31,*) thxyB; close(31)
     open(31,file=trim(adjustl(path)) // 'th_init_T.dat') 
     write(31,'(i5)') 2*kmax
     write(31,'(i5)') 2*lmax
     write(31,*) thxyT; close(31)
  endif

  return
END SUBROUTINE dump
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE advect(u,v,f_x,f_y,fb_y,h_x,h_y,ub,tf,lam,lap)

  ! Spectral advection on the (x,y) grid.

  implicit none

  real, dimension(mmax,nmax), intent(in)  :: u,v,f_x,f_y,fb_y,ub,h_x,h_y,lap
  real,                       intent(in)  :: lam
  real, dimension(mmax,nmax), intent(out) :: tf

  real          :: terr,dx,x,dy,y,famp,rw
  real          :: rphasex,rphasey
  integer       :: i,j
  integer, save :: iseed

  dx = XL/real((mmax)); dy = YL/real((nmax))

! random wave-one forcing; random walk in phase (30 August 2006)
  !rw = 32. ! random wavenumber
  rw = 4.   ! random wavenumber
  famp = 0.1*H ! forcing amplitude
  famp = 0.    ! forcing amplitude
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
END SUBROUTINE advect
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE thadvB(dat,tend,dco,first)

  ! Time-advance SUBROUTINE.

  implicit none

  complex, dimension(2*kmax,2*lmax), intent(inout) :: dat
  complex, dimension(mmax,nmax),     intent(in)    :: tend
  real,                              intent(in)    :: dco
  logical,                           intent(in)    :: first

  real    :: ak,bl
  integer :: k,l,kk,ll
  complex :: ttsp

  complex, dimension(2*kmax,2*lmax), save :: told,told2
  real,                              save :: dts

  logical, parameter :: flag = .TRUE.

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
          &        ak,bl,dco,dts,flag)

     kk = kmax + k; ll = lmax + l
     ! +/-kmax <= k <= -1; 0 <= l <= +/-lmax:
     !         if (k .gt. 1 .and. k .lt. kmax) then
     if (k .gt. 1) then
        ak = -1.*facx*real(kmaxp1 - k); bl = facy*real(l - 1)
        ttsp = tend(k2+k,l)/real(mmax*nmax)
        call tstep_ab(ttsp,dat(kk,l),told(kk,l), &
             &           told2(kk,l),ak,bl,dco,dts,flag)
     endif
     ! 0 <= k <= +/-kmax; +/-lmax <= l <= -1:
     !         if (l .le. lmax) then
     if (l .gt. 1) then
        ak = facx*real(k-1); bl = -1.*facy*real(lmaxp1 - l)
        ttsp = tend(k,l2+l)/real(mmax*nmax)
        call tstep_ab(ttsp,dat(k,ll),told(k,ll), &
             &           told2(k,ll),ak,bl,dco,dts,flag)
     endif
     ! +/-kmax <= k <= -1; +/-lmax <= l <= -1:
     !         if (k .le. kmax .and. l .le. lmax) then
     if ((k .gt. 1) .and. (l .gt. 1)) then
        ak=-1.*facx*real(kmaxp1 - k); bl=-1.*facy*real(lmaxp1 - l)
        ttsp = tend(k2+k,l2+l)/real(mmax*nmax)
        call tstep_ab(ttsp,dat(kk,ll),told(kk,ll), &
             &         told2(kk,ll),ak,bl,dco,dts,flag)
     endif
  enddo; enddo

  return
END SUBROUTINE thadvB
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE thadvT(dat,tend,dco,first)

  ! Time-advance SUBROUTINE.

  implicit none

  complex, dimension(2*kmax,2*lmax), intent(inout) :: dat
  complex, dimension(mmax,nmax),     intent(in)    :: tend
  real,                              intent(in)    :: dco
  logical,                           intent(in)    :: first

  real    :: ak,bl
  integer :: k,l,kk,ll
  complex :: ttsp

  complex, dimension(2*kmax,2*lmax), save :: told,told2
  real,                              save :: dts

  logical, parameter :: flag = .TRUE.

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
          &        ak,bl,dco,dts,flag)

     kk = kmax + k; ll = lmax + l
     ! +/-kmax <= k <= -1; 0 <= l <= +/-lmax:
     !         if (k .gt. 1 .and. k .lt. kmax) then
     if (k .gt. 1) then
        ak = -1.*facx*real(kmaxp1 - k); bl = facy*real(l - 1)
        ttsp = tend(k2+k,l)/real(mmax*nmax)
        call tstep_ab(ttsp,dat(kk,l),told(kk,l), &
             &           told2(kk,l),ak,bl,dco,dts,flag)
     endif
     ! 0 <= k <= +/-kmax; +/-lmax <= l <= -1:
     !         if (l .le. lmax) then
     if (l .gt. 1) then
        ak = facx*real(k-1); bl = -1.*facy*real(lmaxp1 - l)
        ttsp = tend(k,l2+l)/real(mmax*nmax)
        call tstep_ab(ttsp,dat(k,ll),told(k,ll), &
             &           told2(k,ll),ak,bl,dco,dts,flag)
     endif
     ! +/-kmax <= k <= -1; +/-lmax <= l <= -1:
     !         if (k .le. kmax .and. l .le. lmax) then
     if ((k .gt. 1) .and. (l .gt. 1)) then
        ak=-1.*facx*real(kmaxp1 - k); bl=-1.*facy*real(lmaxp1 - l)
        ttsp = tend(k2+k,l2+l)/real(mmax*nmax)
        call tstep_ab(ttsp,dat(kk,ll),told(kk,ll), &
             &         told2(kk,ll),ak,bl,dco,dts,flag)
     endif
  enddo; enddo

  return
END SUBROUTINE thadvT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE tradvB(dat,tend,dco,first)

  ! Time-advance SUBROUTINE.

  implicit none

  complex, dimension(2*kmax,2*lmax), intent(inout) :: dat
  complex, dimension(mmax,nmax),     intent(in)    :: tend
  real,                              intent(in)    :: dco
  logical,                           intent(in)    :: first

  real    :: ak,bl
  integer :: k,l,kk,ll
  complex :: ttsp

  complex, dimension(2*kmax,2*lmax), save :: told,told2
  real,                              save :: dts

  logical, parameter :: flag = .FALSE.

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
          &        ak,bl,dco,dts,flag)

     kk = kmax + k; ll = lmax + l
     ! +/-kmax <= k <= -1; 0 <= l <= +/-lmax:
     !         if (k .gt. 1 .and. k .lt. kmax) then
     if (k .gt. 1) then
        ak = -1.*facx*real(kmaxp1 - k); bl = facy*real(l - 1)
        ttsp = tend(k2+k,l)/real(mmax*nmax)
        call tstep_ab(ttsp,dat(kk,l),told(kk,l), &
             &           told2(kk,l),ak,bl,dco,dts,flag)
     endif
     ! 0 <= k <= +/-kmax; +/-lmax <= l <= -1:
     !         if (l .le. lmax) then
     if (l .gt. 1) then
        ak = facx*real(k-1); bl = -1.*facy*real(lmaxp1 - l)
        ttsp = tend(k,l2+l)/real(mmax*nmax)
        call tstep_ab(ttsp,dat(k,ll),told(k,ll), &
             &           told2(k,ll),ak,bl,dco,dts,flag)
     endif
     ! +/-kmax <= k <= -1; +/-lmax <= l <= -1:
     !         if (k .le. kmax .and. l .le. lmax) then
     if ((k .gt. 1) .and. (l .gt. 1)) then
        ak=-1.*facx*real(kmaxp1 - k); bl=-1.*facy*real(lmaxp1 - l)
        ttsp = tend(k2+k,l2+l)/real(mmax*nmax)
        call tstep_ab(ttsp,dat(kk,ll),told(kk,ll), &
             &         told2(kk,ll),ak,bl,dco,dts,flag)
     endif
  enddo; enddo

  return
END SUBROUTINE tradvB
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE tradvT(dat,tend,dco,first)

  ! Time-advance SUBROUTINE.

  implicit none

  complex, dimension(2*kmax,2*lmax), intent(inout) :: dat
  complex, dimension(mmax,nmax) ,    intent(in)    :: tend
  real,                              intent(in)    :: dco
  logical,                           intent(in)    :: first

  real    :: ak,bl
  integer :: k,l,kk,ll
  complex :: ttsp

  complex, dimension(2*kmax,2*lmax), save :: told,told2
  real,                              save :: dts

  logical, parameter :: flag = .FALSE.

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
          &        ak,bl,dco,dts,flag)

     kk = kmax + k; ll = lmax + l
     ! +/-kmax <= k <= -1; 0 <= l <= +/-lmax:
     !         if (k .gt. 1 .and. k .lt. kmax) then
     if (k .gt. 1) then
        ak = -1.*facx*real(kmaxp1 - k); bl = facy*real(l - 1)
        ttsp = tend(k2+k,l)/real(mmax*nmax)
        call tstep_ab(ttsp,dat(kk,l),told(kk,l), &
             &           told2(kk,l),ak,bl,dco,dts,flag)
     endif
     ! 0 <= k <= +/-kmax; +/-lmax <= l <= -1:
     !         if (l .le. lmax) then
     if (l .gt. 1) then
        ak = facx*real(k-1); bl = -1.*facy*real(lmaxp1 - l)
        ttsp = tend(k,l2+l)/real(mmax*nmax)
        call tstep_ab(ttsp,dat(k,ll),told(k,ll), &
             &           told2(k,ll),ak,bl,dco,dts,flag)
     endif
     ! +/-kmax <= k <= -1; +/-lmax <= l <= -1:
     !         if (k .le. kmax .and. l .le. lmax) then
     if ((k .gt. 1) .and. (l .gt. 1)) then
        ak=-1.*facx*real(kmaxp1 - k); bl=-1.*facy*real(lmaxp1 - l)
        ttsp = tend(k2+k,l2+l)/real(mmax*nmax)
        call tstep_ab(ttsp,dat(kk,ll),told(kk,ll), &
             &         told2(kk,ll),ak,bl,dco,dts,flag)
     endif
  enddo; enddo

  return
END SUBROUTINE tradvT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE tstep_ab(tend,dat,told,told2,ak,bl,dco,dts,flag)

  ! 3rd-order Adams-Bashforth time step (Durran 1991).

  implicit none

  complex, intent(inout) :: dat, told, told2
  complex, intent(in)    :: tend
  real,    intent(in)    :: ak,bl,dco,dts
  logical, intent(in)    :: flag

  complex :: new, temp
  real    :: tfac,relax

  ! changed 02/07/03 to include relaxation to jet
  relax = 0.

  ! set relaxation only if theta and not tracer
  if (flag) then
    if (trl .lt. 1.e3 .and. trl .gt. 0) relax = 1./trl 
  else
    relax = 0.
  endif
  
  ! 12/17/2008: damps l=0 faster for stability (used to zero in main block)
  !if (bl .eq. 0) relax = relax*4

  tfac = dts*((dco*(((ak**2)+(bl**2))**(n/2))) + relax)

  ! new t-step:
  !tfac = dco*dts*((ak**n)+(bl**n))
  !tfac = dco*dts*((ak**2)+(bl**2))**(n/2)

  told = told*exp(-1.*tfac)
  told2 = told2*exp(-tfac)

  temp = (23.*tend) - (16.*told) + (5.*told2)
  new = dat + ( (dts/12.)*(temp) )

  !dat=new*exp(-1.*tfac); told2=told*exp(tfac); told=tend
  dat=new*exp(-1.*tfac); told2=told; told=tend

  return
END SUBROUTINE tstep_ab
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ft_2d(f,ni,nj,isign)

  ! FFT-calling SUBROUTINE. Plug in appropriate calls.

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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE d_setup(dx,dy,dz,dzo,iz,izo,Id)

  ! Set up matrices for derivatives and integrals.

  implicit none

  complex, dimension(2*kmax,2*lmax), intent(out) :: dx,dy,dz,dzo,iz,izo,Id

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
        dz(k,l) = m / tanh(m*H)       ! h->\infty: tanh->1, dz->m
        iz(k,l) = 1./ (m * tanh(m*H)) ! h->\infty: iz-> 1/m
        ! operators for opposing boundary
        if (m*H .lt. 50 .and. model .ne. 0) then
           dzo(k,l) = m / sinh(m*H)       ! h->\infty: dzo->0
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
END SUBROUTINE d_setup
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE d_s2b(temp_in,temp_out,dflag,dn)

  ! small array to big array.
  ! dflag =  n: n derivatives. dflag = -n: n integrations.

  implicit none

  complex, dimension(2*kmax,2*lmax) , intent(in)  :: temp_in,dn
  integer,                            intent(in)  :: dflag
  complex, dimension(mmax,nmax),      intent(out) :: temp_out

  integer :: k,l,kk,ll

  do k = 1,kmax; do l = 1,lmax
     ! waves: 0 <= k <= +/-kmax; 0 <= l <= +/-lmax:

     temp_out(k,l) = (dn(k,l)**dflag)*temp_in(k,l)

     kk = kmax + k; ll = lmax + l

     ! +/-kmax <= k <= -1; 0 <= l <= +/-lmax:
     if (k .gt. 1) then
        temp_out(k2+k,l) = (dn(kk,l)**dflag)*temp_in(kk,l)
     endif

     ! 0 <= k <= +/-kmax; +/-lmax <= l <= -1:
     if (l .gt. 1) then
        temp_out(k,l2+l) = (dn(k,ll)**dflag)*temp_in(k,ll)
     endif

     ! +/-kmax <= k <= -1; +/-lmax <= l <= -1:
     if ((k .gt. 1) .and. (l .gt. 1)) then
        temp_out(k2+k,l2+l) = (dn(kk,ll)**dflag)*temp_in(kk,ll)
     endif
  enddo; enddo

  return
END SUBROUTINE d_s2b
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE d_s2s(temp_in,temp_out,dflag,dn)

  ! small array to small array.
  ! dflag =  n: n derivatives. dflag = -n: n integrations.

  implicit none

  complex, dimension(2*kmax,2*lmax), intent(in)  :: temp_in,dn
  integer,                           intent(in)  :: dflag
  complex, dimension(2*kmax,2*lmax), intent(out) :: temp_out

  integer :: k,l

  temp_out = 0.

  do k = 1,2*kmax; do l = 1,2*lmax
     if (dn(k,l) .ne. 0) temp_out(k,l) = (dn(k,l)**dflag)*temp_in(k,l)
  enddo; enddo

  return
END SUBROUTINE d_s2s
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE d_b2b(temp_in,temp_out,dflag,dn)

  ! big array to big array.
  ! dflag =  n: n derivatives. dflag = -n: n integrations.

  implicit none

  complex, dimension(mmax,nmax),     intent(in)  :: temp_in
  complex, dimension(2*kmax,2*lmax), intent(in)  :: dn
  integer,                           intent(in)  :: dflag
  complex, dimension(mmax,nmax),     intent(out) :: temp_out

  integer :: k,l,kk,ll

  temp_out = 0.

  do k = 1,kmax; do l = 1,lmax
     ! waves: 0 <= k <= +/-kmax; 0 <= l <= +/-lmax:

     if (dn(k,l) .ne. 0) temp_out(k,l) = (dn(k,l)**dflag)*temp_in(k,l)

     kk = kmax + k; ll = lmax + l

     ! +/-kmax <= k <= -1; 0 <= l <= +/-lmax:
     if (k .gt. 1) then
        if (dn(kk,l) .ne. 0)  &
             &           temp_out(k2+k,l) = (dn(kk,l)**dflag)*temp_in(k2+k,l)
     endif

     ! 0 <= k <= +/-kmax; +/-lmax <= l <= -1:
     if (l .gt. 1) then
        if (dn(k,ll) .ne. 0) &
             &           temp_out(k,l2+l) = (dn(k,ll)**dflag)*temp_in(k,l2+l)
     endif

     ! +/-kmax <= k <= -1; +/-lmax <= l <= -1:
     if ((k .gt. 1) .and. (l .gt. 1)) then
        if (dn(kk,ll) .ne. 0)  &
             &           temp_out(k2+k,l2+l) = (dn(kk,ll)**dflag)*temp_in(k2+k,l2+l)
     endif
  enddo; enddo

  return
END SUBROUTINE d_b2b
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE d_b2s(temp_in,temp_out,dflag,dn)

  ! big array to small array.
  ! dflag =  n: n derivatives. dflag = -n: n integrations.

  implicit none

  complex, dimension(mmax,nmax),     intent(in)  :: temp_in
  complex, dimension(2*kmax,2*lmax), intent(in)  :: dn
  integer,                           intent(in)  :: dflag
  complex, dimension(2*kmax,2*lmax), intent(out) :: temp_out
  integer :: k,l,kk,ll

  do k = 1,kmax; do l = 1,lmax
     ! waves: 0 <= k <= +/-kmax; 0 <= l <= +/-lmax:

     temp_out(k,l) = (dn(k,l)**dflag)*temp_in(k,l)

     kk = kmax + k; ll = lmax + l

     ! +/-kmax <= k <= -1; 0 <= l <= +/-lmax:
     if (k .ge. 1) then
        temp_out(kk,l) = (dn(kk,l)**dflag)*temp_in(k2+k,l)
     endif

     ! 0 <= k <= +/-kmax; +/-lmax <= l <= -1:
     if (l .ge. 1) then
        temp_out(k,ll) = (dn(k,ll)**dflag)*temp_in(k,l2+l)
     endif

     ! +/-kmax <= k <= -1; +/-lmax <= l <= -1:
     if ((k .ge. 1) .and. (l .ge. 1)) then
        temp_out(kk,ll) = (dn(kk,ll)**dflag)*temp_in(k2+k,l2+l)
     endif
  enddo; enddo

  return
END SUBROUTINE d_b2s
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Hoskins-West base state FUNCTION calls
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION HW_ubar(y,z)

  implicit none

  real, intent(in)  :: y,z

  real :: HW_ubar
  real :: yp

  yp = y - (0.5*(YL - hwp))
  if (yp .lt. 0.  .and. amu .ne. 0) yp = 0.
  if (yp .gt. hwp .and. amu .ne. 0) yp = hwp
  HW_ubar = shear*z - &
       &     (amu/2.)*(z + ((sinh(ryl*z)/sinh(ryl))*cos(yp*ryl)))

  return
END FUNCTION HW_ubar
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION HW_theta(y,z)

  ! 9/2006: added amp and correction factor for HW 'hiccup'
  ! NOW INCOMPATIBLE WITH HW_ubar!!! MUST solve for U numerically!!!

  implicit none

  real, intent(in)  :: y,z

  real :: HW_theta
  real :: yp

  real, parameter :: amp    = 1.0 !control jet strength (GJH's value 1.50)
  real, parameter :: hiccup = 1.0 !                     (GJH's value 0.75)

  yp = y - (0.5*(YL - hwp))
  if (yp .lt. 0. .and. amu .ne. 0) yp = 0
  if (yp .gt. hwp .and. amu .ne. 0) yp = hwp
  HW_theta = (-shear*yp)+(amu/2.)* &
       &     (yp+(hiccup*(cosh(ryl*z)/sinh(ryl))*sin(yp*ryl)))
  HW_theta = amp*HW_theta

  return
END FUNCTION HW_theta
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION HW_thetay(y,z)

  implicit none

  real, intent(in)  :: y,z

  real :: HW_thetay
  real :: yp

  yp = y - (0.5*(YL - hwp))
  ! fixed this block 02/20/03

  if (amu .eq. 0) then
     HW_thetay = -shear
  else
     if (yp .lt. 0 .or. yp .gt. hwp) then 
        HW_thetay = 0.
     else
        HW_thetay = -shear + (amu/2.) + &
          &        ((amu*ryl/2.)*(cosh(ryl*z)*cos(yp*ryl))/sinh(ryl))
     endif
  endif
  return
END FUNCTION HW_thetay
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! END: Hoskins-West base state FUNCTION calls

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE norm(thB,thT,itime)

  ! hacked from qg code to check mode growth rates

  implicit none

  complex, dimension(2*kmax,2*lmax), intent(in) :: thB,thT
  integer,                           intent(in) :: itime

  real :: V,Vgr,ttau

  real, save :: V_old,V_0

  ! enstrophy norm (for normal modes, any norm will do)
  V = (sum(abs(thB)**2) + sum(abs(thT)**2))
  V = 0.5*(V**0.5)

  ttau = real(itime -5)*dt
  if (itime .eq. 5) then ! wait 5 time steps for norms to stabilize.
     V_0 = V
  endif

  if (itime .gt. 2) then
     !Vgr = ((V - V_old)/(dt*0.5*(V+V_old)))
     Vgr = (log(V) - log(V_old))/dt
  endif

  if (verbose .gt. 0) then
    print*,'V inst growth rates:',itime*dt,Vgr
    !if (itime .gt. 5) print*,'V mean growth rates:',(alog(V/V_0))/ttau
  endif

  V_old = V

  return
END SUBROUTINE norm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION ekman(x,y)

  ! Calculate Ekman pumping

  implicit none

  real             :: ekman
  real, intent(in) :: x,y

  if (x .ge. 2.*XL/3.) then
     ekman = gamma
  else
     ekman = 0.
     ekman = gamma
     !ekman = gamma/4.
  endif

  ekman = gamma
  !if (abs(y-(YL/2)) .ge. 3) ekman = 0.5

  return
END FUNCTION ekman
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE write_diag(output_file,it,thB,thT,trB,trT)

  ! write to disk
  
  implicit none
  
  character(len=*),               intent(in) :: output_file
  integer,                        intent(in) :: it
  real, dimension(2*kmax,2*lmax), intent(in) :: thB, thT, trB, trT

  integer :: ncid, vardim(3), varid
  integer :: k, count(3), start(3), ierr
  real    :: time
 
  if ( it .eq. 0 ) then
     
    !           Create a new NetCDF file
    call nc_check( nf90_create(output_file, ior(NF90_CLOBBER,NF90_64BIT_OFFSET), ncid) )
     
    !           Define dimensions
    call nc_check( nf90_def_dim(ncid, "nx",   2*kmax,         vardim(1)) )
    call nc_check( nf90_def_dim(ncid, "ny",   2*lmax,         vardim(2)) )
    call nc_check( nf90_def_dim(ncid, "time", NF90_UNLIMITED, vardim(3)) )
     
    !           Define variables
    call nc_check( nf90_def_var(ncid, "time",   NF90_FLOAT, vardim(3), varid) )
    call nc_check( nf90_def_var(ncid, "thetaB", NF90_FLOAT, vardim,    varid) )
    call nc_check( nf90_put_att(ncid, varid, "description", "bottom potential temperature") )
    call nc_check( nf90_def_var(ncid, "thetaT", NF90_FLOAT, vardim,    varid) )
    call nc_check( nf90_put_att(ncid, varid, "description", "toppom potential temperature") )
    if (itracer) then
      call nc_check( nf90_def_var(ncid, "tracerB", NF90_FLOAT, vardim,    varid) )
      call nc_check( nf90_put_att(ncid, varid, "description", "bottom tracer") )
      call nc_check( nf90_def_var(ncid, "tracerT", NF90_FLOAT, vardim,    varid) )
      call nc_check( nf90_put_att(ncid, varid, "description", "toppom tracer") )
    endif
     
    call nc_check( nf90_put_att(ncid, NF90_GLOBAL, "model", model) )
    call nc_check( nf90_put_att(ncid, NF90_GLOBAL, "ntims", ntims) )
    call nc_check( nf90_put_att(ncid, NF90_GLOBAL, "dt", dt)       )
    call nc_check( nf90_put_att(ncid, NF90_GLOBAL, "iplot", iplot) )
    call nc_check( nf90_put_att(ncid, NF90_GLOBAL, "XL", XL)       )
    call nc_check( nf90_put_att(ncid, NF90_GLOBAL, "YL", YL)       )
    call nc_check( nf90_put_att(ncid, NF90_GLOBAL, "H", H)         )
    call nc_check( nf90_put_att(ncid, NF90_GLOBAL, "Ross", Ross)   )
    call nc_check( nf90_put_att(ncid, NF90_GLOBAL, "gamma", gamma) )
    call nc_check( nf90_put_att(ncid, NF90_GLOBAL, "n", n)         )
    call nc_check( nf90_put_att(ncid, NF90_GLOBAL, "tau", tau)     )
    call nc_check( nf90_put_att(ncid, NF90_GLOBAL, "trl", trl)     )
    call nc_check( nf90_put_att(ncid, NF90_GLOBAL, "amu", amu)     )
    call nc_check( nf90_put_att(ncid, NF90_GLOBAL, "shear", shear) )

    call nc_check( nf90_enddef(ncid) )
    call nc_check( nf90_close(ncid)  )

  else 
     
    time = (it-1)*dt
    k = 1 + (it-1)/iplot
     
    count(1) = 2*kmax;  start(1) = 1
    count(2) = 2*lmax;  start(2) = 1
    count(3) = 1;       start(3) = k
     
    !           Open the netCDF file, write variables and close
    call nc_check( nf90_open(output_file, NF90_WRITE, ncid) )
    call nc_check( nf90_inq_varid(ncid, "time",   varid) )
    call nc_check( nf90_put_var(ncid, varid, time, (/k/))       )
    call nc_check( nf90_inq_varid(ncid, "thetaB", varid) )
    call nc_check( nf90_put_var(ncid, varid, thB, start, count) )
    call nc_check( nf90_inq_varid(ncid, "thetaT", varid) )
    call nc_check( nf90_put_var(ncid, varid, thT, start, count) )
    if (itracer) then
      call nc_check( nf90_inq_varid(ncid, "tracerB", varid) )
      call nc_check( nf90_put_var(ncid, varid, trB, start, count) )
      call nc_check( nf90_inq_varid(ncid, "tracerT", varid) )
      call nc_check( nf90_put_var(ncid, varid, trT, start, count) )
    endif
    call nc_check( nf90_close(ncid) )
     
  endif
  
  return
END SUBROUTINE write_diag
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE terrain(hx,hy,hu,hv)

! Terrain and barotropic wind on the _advection_ grid:
! feb 2005: added tropo winds from z = 0 for Hsqg

  implicit none
  
  real, intent(out), dimension(mmax,nmax) :: hx,hy,hu,hv

  complex, dimension(mmax,nmax)     :: hxsp,hysp
  complex, dimension(mmax,nmax)     :: hh,hxs,hys

  complex, dimension(2*kmax,2*lmax) :: hspC
  real,    dimension(2*kmax,2*lmax) :: hspR
  complex, dimension(2*kmax,2*lmax) :: Cblank
  real,    dimension(mmax,nmax)     :: Rblank

  real    :: x,y,ddx,ddy,xcen,ycen,amdx,amdy,asx,asy
  real    :: asig,rr,lam
  integer :: i,j

  ! fool's day 2009 for angie
  complex, dimension(2*kmax,2*lmax), save :: dx,dy,dz,dzo,iz,izo,dzi,Id

  ! barotropic wind
  hu = 0.; hv = 0.
  
  ddx = XL/real((mmax)); ddy = YL/real((nmax))
  xcen = mmax*ddx/2; ycen = nmax*ddy/2
  
  if ( verbose .gt. 1 ) print*,'xcen,ycen:',xcen,ycen
  
  do i=1,mmax; do j=1,nmax
     x = (i-1)*ddx; y = (j-1)*ddy
     amdx=min(abs(x-xcen),abs(XL+x-xcen),abs(XL-x+xcen))
     amdy=min(abs(y-ycen),abs(YL+y-ycen),abs(YL-y+ycen))
     asx = 1.0; asy = 1.0; asig = 1.0
     rr = (((amdx/asx)**2) + ((amdy/asy)**2))**.5
     ! gaussian topography:
     hh(i,j) = hamp*exp(-1.*((rr/asig)**2))
  enddo; enddo
  
  if ( verbose .gt. 1 ) print*,real(hh(:,nmax/2))

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

  if ( verbose .gt. 1 ) print*,hx(:,nmax/2)

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
     if ( verbose .gt. 1 ) print*,'max topo=',maxval(abs(hspR))
     call xy_to_sp(cmplx(hspR,0.),hspC,2*kmax,2*lmax,kmax,lmax)
     if ( verbose .gt. 1 ) print*,'max topo spectral=',maxval(abs(hspC))
     call invert(-hspC,Cblank,Rblank,Rblank,Rblank,Rblank,Rblank, & 
          & hv,Rblank,hu,Cblank,Cblank,Rblank,Rblank,Rblank, & 
          & Rblank,.TRUE.,.TRUE.,.TRUE.,lam,Cblank,Cblank,Rblank)
     ! hu and hv have the tropo winds due to topography
     if ( verbose .gt. 1 ) print*,'max tropo winds due to topography: ',maxval(abs(hu)),maxval(abs(hv))
  endif
  
  return
END SUBROUTINE terrain
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE nc_check(ierr)

  ! check for netcdf errors

  implicit none

  integer, intent(in) :: ierr

  if (ierr /= nf90_noerr) then
    print*,'Netcdf error: ', trim( nf90_strerror(ierr) )
    stop
  end if

  return
END SUBROUTINE nc_check
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END PROGRAM sqg_spectral
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
! 1) new terrain call. new ic2 (v2.0) supplies terrain through thxyB.
!
! v.2.6
! 14 July 2005
! 1) netcdf output files enabled (Rahul Mahajan)
!
! v2.6.1
! 19 July 2005
! 1) fixed minor bug in SUBROUTINE dump to catch cdf or matlab output.
!
! v2.6.2
! 10 August 2005
! 1) minor bug fix in dump for linear gradient.
!
! v2.7
! 6 September 2005
! 1) convert main source into a SUBROUTINE for rolv
! 2) pull out diffusion time scale (tau) to SPECTRAL module.
!
! v2.7.1
! 1) minor revs for netcdf include and print of Ross #.
!
! **** check svn logs beyond this point. ****
