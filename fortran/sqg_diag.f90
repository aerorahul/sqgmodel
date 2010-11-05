!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! <next few lines under version control, D O  N O T  E D I T>
! $Date: 2010-08-12 12:03:57 -0700 (Thu, 12 Aug 2010) $
! $Author: rahulm $
! $Revision: 27 $
! $Id: sqg_diag.f90 27 2010-11-05 16:49:46Z rahulm $
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM sqg_diag

	USE NETCDF

  IMPLICIT NONE

	logical, parameter                  :: debug = .false.

	integer                             :: kmax, lmax, tmax, t, n, ens_size
  real, allocatable, dimension(:,:,:) :: thxyB, thxyT, trxyB, trxyT
	character(len=1024)                 :: fpth, etype, smatfile, diagfile
	character(len=5)                    :: nchar
	
	print*, 'ensemble size:'
	read*,  ens_size
	print*, 'path to experiment files:'
	read*,  fpth
	print*, 'type of experiment [ C / P ]:'
	read*,  etype

	diagfile = trim(adjustl(fpth)) // 'smat_' // trim(adjustl(etype)) // '_diag.nc'
	
	do n = 1, ens_size
		write( nchar, '(i5.5)' ) n
		
		if ( mod(n,100) == 0 ) print*,' ... member ', trim(adjustl(nchar))

		smatfile = trim(adjustl(fpth)) // 'smat_' // trim(adjustl(etype)) // '_' // trim(adjustl(nchar)) // '.nc'

		if ( (.not. file_exist(smatfile)) ) then
			print*, trim(adjustl(smatfile)), ' must exist'
			stop
		endif

		if ( n == 1 ) then

			call read_dimensions(smatfile,kmax,lmax,tmax)

			allocate( thxyB(2*kmax,2*lmax,tmax) ) ; allocate( thxyT(2*kmax,2*lmax,tmax) )
			allocate( trxyB(2*kmax,2*lmax,tmax) ) ; allocate( trxyT(2*kmax,2*lmax,tmax) )

		endif

		call read_smat(smatfile,thxyB,thxyT,trxyB,trxyT)

		call write_diag(n,thxyB,thxyT,trxyB,trxyT)

		if ( n == ens_size ) then
			deallocate(thxyB)   ; deallocate(thxyT)
			deallocate(trxyB)   ; deallocate(trxyT)
		endif

	enddo

	call transfer_ga(smatfile,diagfile)
	
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

SUBROUTINE write_diag(n,thB,thT,trB,trT)
  
	implicit none

	integer,                              intent(in)  :: n
	real, dimension(2*kmax,2*lmax, tmax), intent (in) :: thB,thT,trB,trT
	integer                      :: ncid, varid, vardim(4), start(4), count(4), mddim
	character(len=64)            :: copy_string
	character(len=64), parameter :: routine_name = 'write_diag'

	if ( n == 1 ) then
	
		if(debug) print*,'creating file ... ', trim(adjustl(diagfile))

		call check( nf90_create(trim(adjustl(diagfile)), NF90_CLOBBER .or. NF90_64BIT_OFFSET, ncid), routine_name )
	
		call check( nf90_def_dim(ncid, "nx",           2*kmax, vardim(1)), routine_name )
		call check( nf90_def_dim(ncid, "ny",           2*lmax, vardim(2)), routine_name )
		call check( nf90_def_dim(ncid, "time",           tmax, vardim(3)), routine_name )
		call check( nf90_def_dim(ncid, "copy", NF90_UNLIMITED, vardim(4)), routine_name )
		call check( nf90_def_dim(ncid, "metadatalength",   64,     mddim), routine_name )

		call check( nf90_def_var(ncid, "copy", NF90_INT, vardim(4),     varid), routine_name )
		call check( nf90_def_var(ncid, "CopyMetaData", NF90_CHAR, (/mddim, vardim(4)/), varid), routine_name )
		call check( nf90_def_var(ncid, "thetaB", NF90_FLOAT, vardim,    varid), routine_name )
		call check( nf90_put_att(ncid, varid, "description", "bottom potential temperature"), routine_name )
		call check( nf90_def_var(ncid, "thetaT", NF90_FLOAT, vardim,    varid), routine_name )
		call check( nf90_put_att(ncid, varid, "description", "toppom potential temperature"), routine_name )

		if ( trim(adjustl(etype)) == 'P' ) then
			call check( nf90_def_var(ncid, "tracerB", NF90_FLOAT, vardim,    varid), routine_name )
			call check( nf90_put_att(ncid, varid, "description", "bottom tracer"), routine_name )
			call check( nf90_def_var(ncid, "tracerT", NF90_FLOAT, vardim,    varid), routine_name )
			call check( nf90_put_att(ncid, varid, "description", "toppom tracer"), routine_name )
		endif
		
		call check( nf90_put_att(ncid, NF90_GLOBAL, "ens_size", ens_size), routine_name )

		call check( nf90_enddef(ncid), routine_name )
		call check (nf90_close(ncid),  routine_name )

	endif

	if ( n <= ens_size) then
		write(copy_string,'(a15, 1X, i6)') 'ensemble member', n
	elseif ( n == ens_size + 1 ) then
		write(copy_string,'(a13)') 'ensemble mean'
	elseif ( n == ens_size + 2 ) then
		write(copy_string,'(a15)') 'ensemble spread'
	endif

	start(1) = 1      ; start(2) = 1      ; start(3) = 1    ; start(4) = n
	count(1) = 2*kmax ; count(2) = 2*lmax ; count(3) = tmax ; count(4) = 1
	
	call check( nf90_open(trim(adjustl(diagfile)), NF90_WRITE, ncid), routine_name )

	call check( nf90_inq_varid(ncid, "copy",      varid), routine_name )
	call check( nf90_put_var(ncid, varid, n, (/start(4)/)), routine_name )
	call check( nf90_inq_varid(ncid, "CopyMetaData", varid), routine_name )
	call check( nf90_put_var(ncid, varid, copy_string, (/start(1), start(4)/)), routine_name )
	call check( nf90_inq_varid(ncid, "thetaB",    varid ), routine_name )
	call check( nf90_put_var(ncid, varid, thB, start, count), routine_name )
	call check( nf90_inq_varid(ncid, "thetaT",    varid), routine_name )
	call check( nf90_put_var(ncid, varid, thT, start, count), routine_name )
	
	if ( trim(adjustl(etype)) == 'P' ) then
		call check( nf90_inq_varid(ncid, "tracerB",   varid), routine_name )
		call check( nf90_put_var(ncid, varid, trB, start, count), routine_name )
		call check( nf90_inq_varid(ncid, "tracerT",   varid), routine_name )
		call check( nf90_put_var(ncid, varid, trT, start, count), routine_name )
	endif
	
	call check( nf90_close(ncid), routine_name )
  
	return
end SUBROUTINE write_diag

SUBROUTINE transfer_ga(fromfile,tofile)

	implicit none
	character(len=*), intent (in)  :: fromfile, tofile
	character(len=64), parameter   :: routine_name = 'transfer_ga'
	integer                        :: ncidf, ncidt

	call check( nf90_open(trim(adjustl(fromfile)), NF90_NOWRITE, ncidf), routine_name )
	call check( nf90_open(trim(adjustl(tofile)),   NF90_WRITE,   ncidt), routine_name )
	
	call check( nf90_redef(ncidt), routine_name )

	call check( nf90_copy_att(ncidf,NF90_GLOBAL,'model',ncidt,NF90_GLOBAL), routine_name )
	call check( nf90_copy_att(ncidf,NF90_GLOBAL,'dt',   ncidt,NF90_GLOBAL), routine_name )
	call check( nf90_copy_att(ncidf,NF90_GLOBAL,'iplot',ncidt,NF90_GLOBAL), routine_name )
	call check( nf90_copy_att(ncidf,NF90_GLOBAL,'XL',   ncidt,NF90_GLOBAL), routine_name )
	call check( nf90_copy_att(ncidf,NF90_GLOBAL,'YL',   ncidt,NF90_GLOBAL), routine_name )
	call check( nf90_copy_att(ncidf,NF90_GLOBAL,'H',    ncidt,NF90_GLOBAL), routine_name )
	call check( nf90_copy_att(ncidf,NF90_GLOBAL,'Ross', ncidt,NF90_GLOBAL), routine_name )
	call check( nf90_copy_att(ncidf,NF90_GLOBAL,'gamma',ncidt,NF90_GLOBAL), routine_name )
	call check( nf90_copy_att(ncidf,NF90_GLOBAL,'n',    ncidt,NF90_GLOBAL), routine_name )
	call check( nf90_copy_att(ncidf,NF90_GLOBAL,'tau',  ncidt,NF90_GLOBAL), routine_name )
	call check( nf90_copy_att(ncidf,NF90_GLOBAL,'trl',  ncidt,NF90_GLOBAL), routine_name )
	call check( nf90_copy_att(ncidf,NF90_GLOBAL,'amu',  ncidt,NF90_GLOBAL), routine_name )
	call check( nf90_copy_att(ncidf,NF90_GLOBAL,'shear',ncidt,NF90_GLOBAL), routine_name )

	call check( nf90_enddef(ncidt), routine_name )
	
	call check( nf90_close(ncidf), routine_name )
	call check( nf90_close(ncidt), routine_name )
  
	return
end SUBROUTINE transfer_ga 
  
subroutine check(ierr, routine_name)

	implicit none
	integer,          intent (in) :: ierr
	character(len=*), intent (in) :: routine_name

	if (ierr /= nf90_noerr) then
  	print*,'Netcdf error: ', trim(adjustl(nf90_strerror(ierr)))
  	print*,'in subroutine ', trim(adjustl(routine_name))
	  stop
	endif

  return
end subroutine check
  
END PROGRAM sqg_diag
