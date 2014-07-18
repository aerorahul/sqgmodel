!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! <next few lines under version control, D O  N O T  E D I T>
! $Date$
! $Author$
! $Revision$
! $Id$
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


PROGRAM add_or_remove_base

	USE NETCDF

  IMPLICIT NONE

	logical, parameter :: debug = .true.

	integer                             :: kmax, lmax, tmax, t
  real, allocatable, dimension(:,:)   :: thbxyB,thbxyT
  real, allocatable, dimension(:,:,:) :: thxyB, thxyT
	character(len=64)                   :: basefile, smatfile, action_to_perform
	logical                             :: add_base, remove_base

	! Files to manipulate
	basefile = 'base.nc'
	smatfile = 'smat.nc'
	
	! Actions to perform
	action_to_perform = ''
	add_base          = .FALSE.
	remove_base       = .FALSE.
	
	print*,'add / remove base state to / from the smat file'
	print*,'ENTER, add OR remove'
	call getarg(1,action_to_perform)
	if ( trim(adjustl(action_to_perform)) == 'add' ) then
		print*,' ... adding base state'
		add_base    = .TRUE.
		remove_base = .FALSE.
	elseif ( trim(adjustl(action_to_perform)) == 'remove' ) then
		print*,' ... removing base state:'
		add_base    = .FALSE.
		remove_base = .TRUE.
	else
		print*,trim(adjustl(action_to_perform)), ' operation is not defined'
		print*,'try : add or remove'
		stop
	endif
	
	! check if the necessary files exist
	if ( .not. file_exist(basefile) ) then
		print*,trim(adjustl(basefile)), ' does not exist'
		stop
	endif
	if ( .not. file_exist(smatfile) ) then
		print*,trim(adjustl(smatfile)), ' does not exist'
		stop
	endif

	! read dimensions from smatfile
	call read_dimensions(smatfile,kmax,lmax,tmax)
	if (debug) print*,'nx, ny, nt = ', 2*kmax, 2*lmax, tmax

	! allocate
	allocate(thbxyB(2*kmax,2*lmax));      allocate(thbxyT(2*kmax,2*lmax))
	allocate( thxyB(2*kmax,2*lmax,tmax)); allocate( thxyT(2*kmax,2*lmax,tmax))
	
	! read from basefile
	call read_base(basefile,thbxyB,thbxyT)

	! read from smatfile
	call read_smat(smatfile,thxyB,thxyT)

	! add / remove  base state
	if ( add_base ) then
		do t = 1, tmax
			thxyB(:,:,t) = thxyB(:,:,t) + thbxyB(:,:)
			thxyT(:,:,t) = thxyT(:,:,t) + thbxyT(:,:)
		enddo
	elseif ( remove_base ) then
		do t = 1, tmax
			thxyB(:,:,t) = thxyB(:,:,t) - thbxyB(:,:)
			thxyT(:,:,t) = thxyT(:,:,t) - thbxyT(:,:)
		enddo
	endif

	! write to smatfile
	call write_smat(smatfile,thxyB,thxyT)

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

	if(debug) print*,'reading data from ... ',  trim(adjustl(filename))
	call check( nf90_open(trim(adjustl(filename)), NF90_WRITE, ncid) )
	call check( nf90_inq_dimid(ncid, 'time', dimid) )
	call check( nf90_inquire_dimension(ncid, dimid, len = tmax_out) )
	call check( nf90_inq_dimid(ncid, 'nx', dimid) )
	call check( nf90_inquire_dimension(ncid, dimid, len = twokmax) )
	call check( nf90_inq_dimid(ncid, 'ny', dimid) )
	call check( nf90_inquire_dimension(ncid, dimid, len = twolmax) )
	call check (nf90_close(ncid) )
	kmax_out = twokmax / 2
	lmax_out = twolmax / 2

	return
end SUBROUTINE read_dimensions

SUBROUTINE read_base(filename,thB,thT)

  IMPLICIT NONE
	
	character(len=*),               intent (in)  :: filename
  real, dimension(2*kmax,2*lmax), intent (out) :: thB,thT
	integer :: ncid, varid

  thB = 0.; thT = 0.
	
	if (debug) print*,'reading data from ... ', trim(adjustl(filename))
	call check( nf90_open(trim(adjustl(filename)), NF90_NOWRITE, ncid) )
	call check( nf90_inq_varid(ncid, 'thetaB', varid) )
	call check( nf90_get_var(ncid, varid, thB) )
	call check( nf90_inq_varid(ncid, 'thetaT', varid) )
	call check( nf90_get_var(ncid, varid, thT) )
	call check (nf90_close(ncid) )

  return
end SUBROUTINE read_base

SUBROUTINE read_smat(filename,thB,thT)

  IMPLICIT NONE
	
	character(len=*),                    intent (in)  :: filename
  real, dimension(2*kmax,2*lmax,tmax), intent (out) :: thB,thT
	integer :: ncid, varid

  thB = 0.; thT = 0.

	if(debug) print*,'reading data from ... ', trim(adjustl(filename))
	call check( nf90_open(trim(adjustl(filename)), NF90_NOWRITE, ncid) )
	call check( nf90_inq_varid(ncid, 'thetaB', varid) )
	call check( nf90_get_var(ncid, varid, thB) )
	call check( nf90_inq_varid(ncid, 'thetaT', varid) )
	call check( nf90_get_var(ncid, varid, thT) )
	call check (nf90_close(ncid) )

  return
end SUBROUTINE read_smat

SUBROUTINE write_smat(filename,thB,thT)

  IMPLICIT NONE

	character(len=*),                     intent (in) :: filename
  real, dimension(2*kmax,2*lmax, tmax), intent (in) :: thB,thT
	integer :: ncid, varid

	if(debug) print*,'writing data to   ... ',  trim(adjustl(filename))
	call check( nf90_open(trim(adjustl(filename)), NF90_WRITE, ncid) )
	call check( nf90_inq_varid(ncid, 'thetaB', varid) )
	call check( nf90_put_var(ncid, varid, thB) )
	call check( nf90_inq_varid(ncid, 'thetaT', varid) )
	call check( nf90_put_var(ncid, varid, thT) )
	call check (nf90_close(ncid) )

  return
end SUBROUTINE write_smat 

subroutine check(ierr)

	implicit none
	integer, intent (in) :: ierr

	if (ierr /= nf90_noerr) then
  	print*,'Netcdf error: ', trim( nf90_strerror(ierr) )
	  stop
	endif

  return
end subroutine check
  
END PROGRAM add_or_remove_base
