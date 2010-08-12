!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! <next few lines under version control, D O  N O T  E D I T>
! $Date: 2010-08-12 12:03:57 -0700 (Thu, 12 Aug 2010) $
! $Author$
! $Revision$
! $Id$
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM sqg_diag

	USE NETCDF

  IMPLICIT NONE

	logical, parameter                  :: debug = .false.

	integer                             :: kmax, lmax, tmax, t, n, ens_size
  real, allocatable, dimension(:,:,:) :: thxyB_c, thxyT_c, trxyB_c, trxyT_c
  real, allocatable, dimension(:,:,:) :: thxyB_p, thxyT_p, trxyB_p, trxyT_p
  real, allocatable, dimension(:,:,:) :: thxyB, thxyT, trxyB, trxyT
  real, allocatable, dimension(:,:,:) :: thxyB_sum, thxyT_sum, trxyB_sum, trxyT_sum
  real, allocatable, dimension(:,:,:) :: thxyB_mean, thxyT_mean, trxyB_mean, trxyT_mean
  real, allocatable, dimension(:,:,:) :: thxyB_spread, thxyT_spread, trxyB_spread, trxyT_spread
	character(len=64)                   :: charbuf, smatCfile, smatPfile
	character(len=5)                    :: nchar
	
	call getarg(1,charbuf)
	read(charbuf,'(I5)') ens_size
	
	print*,'calculating ensemble mean ...'
	do n = 1, ens_size
		write( nchar, '(i5.5)' ) n
		
		if ( mod(n,100) == 0 ) print*,' ... member ', trim(adjustl(nchar))

		smatCfile = 'smat_C_' // trim(adjustl(nchar)) // '.nc'
		smatPfile = 'smat_P_' // trim(adjustl(nchar)) // '.nc'

		if ( (.not. file_exist(smatCfile))  .and. ( .not. file_exist(smatPfile)) ) then
			print*, 'Both ', trim(adjustl(smatCfile)), ' and ', trim(adjustl(smatPfile)), ' must exist'
			stop
		endif

		if ( n == 1 ) then
			
			call read_dimensions(smatCfile,kmax,lmax,tmax)
			
			allocate( thxyB_c(  2*kmax,2*lmax,tmax) ) ; allocate( thxyT_c(  2*kmax,2*lmax,tmax) )
			allocate( trxyB_c(  2*kmax,2*lmax,tmax) ) ; allocate( trxyT_c(  2*kmax,2*lmax,tmax) )
			allocate( thxyB_p(  2*kmax,2*lmax,tmax) ) ; allocate( thxyT_p(  2*kmax,2*lmax,tmax) )
			allocate( trxyB_p(  2*kmax,2*lmax,tmax) ) ; allocate( trxyT_p(  2*kmax,2*lmax,tmax) )
			allocate( thxyB(    2*kmax,2*lmax,tmax) ) ; allocate( thxyT(    2*kmax,2*lmax,tmax) )
			allocate( trxyB(    2*kmax,2*lmax,tmax) ) ; allocate( trxyT(    2*kmax,2*lmax,tmax) )
			allocate( thxyB_sum(2*kmax,2*lmax,tmax) ) ; allocate( thxyT_sum(2*kmax,2*lmax,tmax) )
			allocate( trxyB_sum(2*kmax,2*lmax,tmax) ) ; allocate( trxyT_sum(2*kmax,2*lmax,tmax) )

			thxyB_sum = 0. ; thxyT_sum = 0.
			trxyB_sum = 0. ; trxyT_sum = 0.

		endif

		call read_smat(smatCfile,thxyB_c,thxyT_c,trxyB_c,trxyT_c)
		call read_smat(smatPfile,thxyB_p,thxyT_p,trxyB_p,trxyT_p)

		thxyB = thxyB_p - thxyB_c ; thxyT = thxyT_p - thxyT_c
		trxyB = trxyB_p           ; trxyT = trxyT_p

		thxyB_sum = thxyB_sum + thxyB ; thxyT_sum = thxyT_sum + thxyT
		trxyB_sum = trxyB_sum + trxyB ; trxyT_sum = trxyT_sum + trxyT

		if ( n == ens_size ) then
			deallocate(thxyB_c) ; deallocate(thxyT_c)
			deallocate(trxyB_c) ; deallocate(trxyT_c)
			deallocate(thxyB_p) ; deallocate(thxyT_p)
			deallocate(trxyB_p) ; deallocate(trxyT_p)
			deallocate(thxyB)   ; deallocate(thxyT)
			deallocate(trxyB)   ; deallocate(trxyT)
		endif

	enddo
	
	allocate( thxyB_mean(2*kmax,2*lmax,tmax) ) ; allocate( thxyT_mean(2*kmax,2*lmax,tmax) )
	allocate( trxyB_mean(2*kmax,2*lmax,tmax) ) ; allocate( trxyT_mean(2*kmax,2*lmax,tmax) )

	thxyB_mean = thxyB_sum / ens_size	; thxyT_mean = thxyT_sum / ens_size
	trxyB_mean = trxyB_sum / ens_size	; trxyT_mean = trxyT_sum / ens_size
			
	deallocate(thxyB_sum) ; deallocate(thxyT_sum)
	deallocate(trxyB_sum) ; deallocate(trxyT_sum)

	print*,'calculating ensemble spread ...'
	do n = 1, ens_size
		write( nchar, '(i5.5)' ) n
		
		if ( mod(n,100) == 0 ) print*,' ... member ', trim(adjustl(nchar))

		smatCfile = 'smat_C_' // trim(adjustl(nchar)) // '.nc'
		smatPfile = 'smat_P_' // trim(adjustl(nchar)) // '.nc'

		if ( n == 1 ) then
			
			allocate( thxyB_c(  2*kmax,2*lmax,tmax) ) ; allocate( thxyT_c(  2*kmax,2*lmax,tmax) )
			allocate( trxyB_c(  2*kmax,2*lmax,tmax) ) ; allocate( trxyT_c(  2*kmax,2*lmax,tmax) )
			allocate( thxyB_p(  2*kmax,2*lmax,tmax) ) ; allocate( thxyT_p(  2*kmax,2*lmax,tmax) )
			allocate( trxyB_p(  2*kmax,2*lmax,tmax) ) ; allocate( trxyT_p(  2*kmax,2*lmax,tmax) )
			allocate( thxyB(    2*kmax,2*lmax,tmax) ) ; allocate( thxyT(    2*kmax,2*lmax,tmax) )
			allocate( trxyB(    2*kmax,2*lmax,tmax) ) ; allocate( trxyT(    2*kmax,2*lmax,tmax) )
			allocate( thxyB_sum(2*kmax,2*lmax,tmax) ) ; allocate( thxyT_sum(2*kmax,2*lmax,tmax) )
			allocate( trxyB_sum(2*kmax,2*lmax,tmax) ) ; allocate( trxyT_sum(2*kmax,2*lmax,tmax) )

			thxyB_sum  = 0. ; thxyT_sum  = 0.
			trxyB_sum  = 0. ; trxyT_sum  = 0.

		endif

		call read_smat(smatCfile,thxyB_c,thxyT_c,trxyB_c,trxyT_c)
		call read_smat(smatPfile,thxyB_p,thxyT_p,trxyB_p,trxyT_p)

		thxyB = thxyB_p - thxyB_c ; thxyT = thxyT_p - thxyT_c
		trxyB = trxyB_p           ; trxyT = trxyT_p

		thxyB_sum = thxyB_sum + ( thxyB - thxyB_mean )**2 ; thxyT_sum = thxyT_sum + ( thxyT -	thxyT_mean )**2
		trxyB_sum = trxyB_sum + ( trxyB - trxyB_mean )**2 ; trxyT_sum = trxyT_sum + ( trxyT - trxyT_mean )**2

		if ( n == ens_size ) then
			deallocate(thxyB_c) ; deallocate(thxyT_c)
			deallocate(trxyB_c) ; deallocate(trxyT_c)
			deallocate(thxyB_p) ; deallocate(thxyT_p)
			deallocate(trxyB_p) ; deallocate(trxyT_p)
			deallocate(thxyB)   ; deallocate(thxyT)
			deallocate(trxyB)   ; deallocate(trxyT)
		endif

	enddo
	
	allocate( thxyB_spread(2*kmax,2*lmax,tmax) ) ; allocate( thxyT_spread(2*kmax,2*lmax,tmax) )
	allocate( trxyB_spread(2*kmax,2*lmax,tmax) ) ; allocate( trxyT_spread(2*kmax,2*lmax,tmax) )

	thxyB_spread = sqrt(thxyB_sum / (ens_size - 1)) ; thxyT_spread = sqrt(thxyT_sum / (ens_size - 1))
	trxyB_spread = sqrt(trxyB_sum / (ens_size - 1)) ; trxyT_spread = sqrt(trxyT_sum / (ens_size - 1))
			
	deallocate(thxyB_sum)  ; deallocate(thxyT_sum)
	deallocate(trxyB_sum)  ; deallocate(trxyT_sum)
	
	call write_diag(thxyB_mean,   thxyT_mean,   trxyB_mean,   trxyT_mean, &
	                thxyB_spread, thxyT_spread, trxyB_spread, trxyT_spread)

	deallocate(thxyB_mean)   ; deallocate(thxyT_mean)
	deallocate(trxyB_mean)   ; deallocate(trxyT_mean)
	deallocate(thxyB_spread) ; deallocate(thxyT_spread)
	deallocate(trxyB_spread) ; deallocate(trxyT_spread)
	
	stop

CONTAINS

FUNCTION file_exist(filename)

	implicit none
	character(len=*), intent (in) :: filename
	logical                       :: file_exist
	
	inquire(file=trim(adjustl(filename)), exist=file_exist)

  return
end FUNCTION file_exist

SUBROUTINE read_dimensions(filename,kmax_out,lmax_out,tmax_out)
	
	implicit none
	character(len=*), intent (in)  :: filename
	integer,          intent (out) :: kmax_out, lmax_out, tmax_out
	integer                        :: ncid, dimid, twokmax, twolmax
	character(len=64), parameter   :: routine_name = 'read_dimensions'

	if(debug) print*,'reading data from ... ',  trim(adjustl(filename))
	call check( nf90_open(trim(adjustl(filename)), NF90_WRITE, ncid), routine_name )
	call check( nf90_inq_dimid(ncid, 'time', dimid), routine_name )
	call check( nf90_inquire_dimension(ncid, dimid, len = tmax_out), routine_name )
	call check( nf90_inq_dimid(ncid, 'nx', dimid), routine_name )
	call check( nf90_inquire_dimension(ncid, dimid, len = twokmax), routine_name )
	call check( nf90_inq_dimid(ncid, 'ny', dimid), routine_name )
	call check( nf90_inquire_dimension(ncid, dimid, len = twolmax), routine_name )
	call check (nf90_close(ncid), routine_name )
	kmax_out = twokmax / 2
	lmax_out = twolmax / 2

	return
end SUBROUTINE read_dimensions

SUBROUTINE read_smat(filename,thB,thT,trB,trT)

  IMPLICIT NONE
	
	character(len=*),                    intent (in)  :: filename
  real, dimension(2*kmax,2*lmax,tmax), intent (out) :: thB,thT,trB,trT
	integer                      :: ncid, varid, ierr
	character(len=64), parameter :: routine_name = 'read_smat'

  thB = 0.; thT = 0.
  trB = 0.; trT = 0.

	if(debug) print*,'reading data from ... ', trim(adjustl(filename))
	call check( nf90_open(trim(adjustl(filename)), NF90_NOWRITE, ncid), routine_name )
	call check( nf90_inq_varid(ncid, 'thetaB', varid), routine_name )
	call check( nf90_get_var(ncid, varid, thB), routine_name )
	call check( nf90_inq_varid(ncid, 'thetaT', varid), routine_name )
	call check( nf90_get_var(ncid, varid, thT), routine_name )
	
	ierr = nf90_inq_varid(ncid, 'tracerB', varid)
	if ( ierr == nf90_noerr ) then
		call check( nf90_inq_varid(ncid, 'tracerB', varid), routine_name )
		call check( nf90_get_var(ncid, varid, trB), routine_name )
		call check( nf90_inq_varid(ncid, 'tracerT', varid), routine_name )
		call check( nf90_get_var(ncid, varid, trT), routine_name )
	endif

	call check (nf90_close(ncid), routine_name )

  return
end SUBROUTINE read_smat

SUBROUTINE write_diag(thBm,thTm,trBm,trTm,thBs,thTs,trBs,trTs)

  IMPLICIT NONE

  real, dimension(2*kmax,2*lmax, tmax), intent (in) :: thBm,thTm,thBs,thTs
  real, dimension(2*kmax,2*lmax, tmax), intent (in) :: trBm,trTm,trBs,trTs
	integer                      :: ncid, varid, vardim(3)
	character(len=64)            :: diag_file = 'smat_diag.nc'
	character(len=64), parameter :: routine_name = 'write_diag'

	if(debug) print*,'writing data to   ... ',  trim(adjustl(diag_file))

	call check( nf90_create(trim(adjustl(diag_file)), NF90_CLOBBER .or. NF90_64BIT_OFFSET, ncid), routine_name )

	call check( nf90_def_dim(ncid, "nx",   2*kmax, vardim(1)), routine_name )
	call check( nf90_def_dim(ncid, "ny",   2*lmax, vardim(2)), routine_name )
	call check( nf90_def_dim(ncid, "time",   tmax, vardim(3)), routine_name )

	call check( nf90_def_var(ncid, "thetaB_mean", NF90_FLOAT, vardim,    varid), routine_name )
	call check( nf90_put_att(ncid, varid, "description", "ensemble mean bottom potential temperature"), routine_name )
	call check( nf90_def_var(ncid, "thetaT_mean", NF90_FLOAT, vardim,    varid), routine_name )
	call check( nf90_put_att(ncid, varid, "description", "ensemble mean toppom potential temperature"), routine_name )
	call check( nf90_def_var(ncid, "tracerB_mean", NF90_FLOAT, vardim,    varid), routine_name )
	call check( nf90_put_att(ncid, varid, "description", "ensemble mean bottom tracer"), routine_name )
	call check( nf90_def_var(ncid, "tracerT_mean", NF90_FLOAT, vardim,    varid), routine_name )
	call check( nf90_put_att(ncid, varid, "description", "ensemble mean toppom tracer"), routine_name )
	call check( nf90_def_var(ncid, "thetaB_spread", NF90_FLOAT, vardim,    varid), routine_name )
	call check( nf90_put_att(ncid, varid, "description", "ensemble spread bottom potential temperature"), routine_name )
	call check( nf90_def_var(ncid, "thetaT_spread", NF90_FLOAT, vardim,    varid), routine_name )
	call check( nf90_put_att(ncid, varid, "description", "ensemble spread toppom potential temperature"), routine_name )
	call check( nf90_def_var(ncid, "tracerB_spread", NF90_FLOAT, vardim,    varid), routine_name )
	call check( nf90_put_att(ncid, varid, "description", "ensemble spread bottom tracer"), routine_name )
	call check( nf90_def_var(ncid, "tracerT_spread", NF90_FLOAT, vardim,    varid), routine_name )
	call check( nf90_put_att(ncid, varid, "description", "ensemble spread toppom tracer"), routine_name )
	
	call check( nf90_put_att(ncid, NF90_GLOBAL, "ens_size", ens_size), routine_name )

	call check( nf90_enddef(ncid), routine_name )
	call check( nf90_close(ncid), routine_name  )

	call check( nf90_open(trim(adjustl(diag_file)), NF90_WRITE, ncid), routine_name )

	call check( nf90_inq_varid(ncid, "thetaB_mean",    varid), routine_name )
	call check( nf90_put_var(ncid, varid, thBm), routine_name )
	call check( nf90_inq_varid(ncid, "thetaT_mean",    varid), routine_name )
	call check( nf90_put_var(ncid, varid, thTm), routine_name )
	call check( nf90_inq_varid(ncid, "tracerB_mean",   varid), routine_name )
	call check( nf90_put_var(ncid, varid, trBm), routine_name )
	call check( nf90_inq_varid(ncid, "tracerT_mean",   varid), routine_name )
	call check( nf90_put_var(ncid, varid, trTm), routine_name )
	call check( nf90_inq_varid(ncid, "thetaB_spread",  varid), routine_name )
	call check( nf90_put_var(ncid, varid, thBs), routine_name )
	call check( nf90_inq_varid(ncid, "thetaT_spread",  varid), routine_name )
	call check( nf90_put_var(ncid, varid, thTs), routine_name )
	call check( nf90_inq_varid(ncid, "tracerB_spread", varid), routine_name )
	call check( nf90_put_var(ncid, varid, trBs), routine_name )
	call check( nf90_inq_varid(ncid, "tracerT_spread", varid), routine_name )
	call check( nf90_put_var(ncid, varid, trTs), routine_name )
	
	call check( nf90_close(ncid), routine_name  )

  return
end SUBROUTINE write_diag

subroutine check(ierr, routine_name)

	implicit none
	integer,          intent (in) :: ierr
	character(len=*), intent (in) :: routine_name

	if (ierr /= nf90_noerr) then
  	print*,'Netcdf error: ', trim(adjustl(nf90_strerror(ierr)))
  	print*,'in ', trim(adjustl(routine_name))
	  stop
	endif

  return
end subroutine check
  
END PROGRAM sqg_diag
