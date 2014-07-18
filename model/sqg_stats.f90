!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! <next few lines under version control, D O  N O T  E D I T>
! $Date: 2010-08-12 12:03:57 -0700 (Thu, 12 Aug 2010) $
! $Author$
! $Revision$
! $Id$
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM sqg_stats

	USE NETCDF

  IMPLICIT NONE

	logical, parameter                  :: debug = .false.

	integer                             :: kmax, lmax, tmax, t, n, ens_size
	integer                             :: ncidC, ncidP, ncidD
  real, allocatable, dimension(:,:,:) :: thxyB_c, thxyT_c, trxyB_c, trxyT_c
  real, allocatable, dimension(:,:,:) :: thxyB_p, thxyT_p, trxyB_p, trxyT_p
  real, allocatable, dimension(:,:,:) :: thxyB, thxyT, trxyB, trxyT
  real, allocatable, dimension(:,:,:) :: thxyB_sum, thxyT_sum, trxyB_sum, trxyT_sum
  real, allocatable, dimension(:,:,:) :: thxyB_mean, thxyT_mean, trxyB_mean, trxyT_mean
  real, allocatable, dimension(:,:,:) :: thxyB_spread, thxyT_spread, trxyB_spread, trxyT_spread
	character(len=1024)                 :: fpth_c, fpth_p, smatCfile, smatPfile, smatDfile
	
	print*, 'path to control experiment diag file:'
	read*,  fpth_c
	print*, 'path to perturbed experiment diag file:'
	read*,  fpth_p
	
	smatCfile = trim(adjustl(fpth_c)) // 'smat_C_diag.nc'
	smatPfile = trim(adjustl(fpth_p)) // 'smat_P_diag.nc'
	smatDfile = trim(adjustl(fpth_p)) // 'smat_D_diag.nc'

	if ( (.not. file_exist(smatCfile))  .and. ( .not. file_exist(smatPfile)) ) then
		print*, 'Both ', trim(adjustl(smatCfile)), ' and ', trim(adjustl(smatPfile)), ' must exist'
		stop
	endif

	call check_smat()
	
	call read_dimensions(smatCfile,ens_size,kmax,lmax,tmax)

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

	call check( nf90_open(trim(adjustl(smatCfile)), NF90_NOWRITE, ncidC), 'sqg_stats')
	call check( nf90_open(trim(adjustl(smatPfile)), NF90_NOWRITE, ncidP), 'sqg_stats')
	
	print*,'calculating ensemble mean ...'
	do n = 1, ens_size

		if ( mod(n,100) == 0 ) write(6,'(a11, 1x, i6)') ' ... member', n 

		call read_smat(ncidC,n,thxyB_c,thxyT_c,trxyB_c,trxyT_c)
		call read_smat(ncidP,n,thxyB_p,thxyT_p,trxyB_p,trxyT_p)

		thxyB = thxyB_p - thxyB_c ; thxyT = thxyT_p - thxyT_c
		trxyB = trxyB_p           ; trxyT = trxyT_p

		call write_diag(n,thxyB,thxyT,trxyB,trxyT)

		thxyB_sum = thxyB_sum + thxyB ; thxyT_sum = thxyT_sum + thxyT
		trxyB_sum = trxyB_sum + trxyB ; trxyT_sum = trxyT_sum + trxyT

	enddo
	
	call check( nf90_close(ncidC), 'sqg_stats')
	call check( nf90_close(ncidP), 'sqg_stats')

	deallocate(thxyB_c) ; deallocate(thxyT_c)
	deallocate(trxyB_c) ; deallocate(trxyT_c)
	deallocate(thxyB_p) ; deallocate(thxyT_p)
	deallocate(trxyB_p) ; deallocate(trxyT_p)
	deallocate(thxyB)   ; deallocate(thxyT)
	deallocate(trxyB)   ; deallocate(trxyT)

	allocate( thxyB_mean(2*kmax,2*lmax,tmax) ) ; allocate( thxyT_mean(2*kmax,2*lmax,tmax) )
	allocate( trxyB_mean(2*kmax,2*lmax,tmax) ) ; allocate( trxyT_mean(2*kmax,2*lmax,tmax) )

	thxyB_mean = thxyB_sum / ens_size	; thxyT_mean = thxyT_sum / ens_size
	trxyB_mean = trxyB_sum / ens_size	; trxyT_mean = trxyT_sum / ens_size
			
	deallocate(thxyB_sum) ; deallocate(thxyT_sum)
	deallocate(trxyB_sum) ; deallocate(trxyT_sum)

	allocate( thxyB(    2*kmax,2*lmax,tmax) ) ; allocate( thxyT(    2*kmax,2*lmax,tmax) )
	allocate( trxyB(    2*kmax,2*lmax,tmax) ) ; allocate( trxyT(    2*kmax,2*lmax,tmax) )
	allocate( thxyB_sum(2*kmax,2*lmax,tmax) ) ; allocate( thxyT_sum(2*kmax,2*lmax,tmax) )
	allocate( trxyB_sum(2*kmax,2*lmax,tmax) ) ; allocate( trxyT_sum(2*kmax,2*lmax,tmax) )

	thxyB_sum  = 0. ; thxyT_sum  = 0.
	trxyB_sum  = 0. ; trxyT_sum  = 0.

	call check( nf90_open(trim(adjustl(smatDfile)), NF90_NOWRITE, ncidD), 'sqg_stats')

	print*,'calculating ensemble spread ...'
	do n = 1, ens_size

		if ( mod(n,100) == 0 ) write(6,'(a11, 1x, i6)') ' ... member', n 

		call read_smat(ncidD,n,thxyB,thxyT,trxyB,trxyT)

		thxyB_sum = thxyB_sum + ( thxyB - thxyB_mean )**2 ; thxyT_sum = thxyT_sum + ( thxyT -	thxyT_mean )**2
		trxyB_sum = trxyB_sum + ( trxyB - trxyB_mean )**2 ; trxyT_sum = trxyT_sum + ( trxyT - trxyT_mean )**2

	enddo
	
	call check( nf90_close(ncidD), 'sqg_stats')
	
	deallocate(thxyB)   ; deallocate(thxyT)
	deallocate(trxyB)   ; deallocate(trxyT)

	allocate( thxyB_spread(2*kmax,2*lmax,tmax) ) ; allocate( thxyT_spread(2*kmax,2*lmax,tmax) )
	allocate( trxyB_spread(2*kmax,2*lmax,tmax) ) ; allocate( trxyT_spread(2*kmax,2*lmax,tmax) )

	thxyB_spread = sqrt(thxyB_sum / (ens_size - 1)) ; thxyT_spread = sqrt(thxyT_sum / (ens_size - 1))
	trxyB_spread = sqrt(trxyB_sum / (ens_size - 1)) ; trxyT_spread = sqrt(trxyT_sum / (ens_size - 1))
			
	deallocate(thxyB_sum)  ; deallocate(thxyT_sum)
	deallocate(trxyB_sum)  ; deallocate(trxyT_sum)
	
	call write_diag(ens_size+1,thxyB_mean,  thxyT_mean,  trxyB_mean,  trxyT_mean)
	call write_diag(ens_size+2,thxyB_spread,thxyT_spread,trxyB_spread,trxyT_spread)

	deallocate(thxyB_mean)   ; deallocate(thxyT_mean)
	deallocate(trxyB_mean)   ; deallocate(trxyT_mean)
	deallocate(thxyB_spread) ; deallocate(thxyT_spread)
	deallocate(trxyB_spread) ; deallocate(trxyT_spread)
	
	write(6,'(a14)') 'sqg_stats done'

	stop

CONTAINS

FUNCTION file_exist(filename)

	implicit none
	character(len=*), intent (in) :: filename
	logical                       :: file_exist
	
	inquire(file=trim(adjustl(filename)), exist=file_exist)

  return
end FUNCTION file_exist

SUBROUTINE read_dimensions(filename,copy_out,kmax_out,lmax_out,tmax_out)
	
	implicit none
	character(len=*), intent (in)  :: filename
	integer,          intent (out) :: copy_out,kmax_out, lmax_out, tmax_out
	integer                        :: ncid, dimid, twokmax, twolmax
	character(len=64), parameter   :: routine_name = 'read_dimensions'

	if(debug) print*,'reading data from ... ',  trim(adjustl(filename))
	call check( nf90_open(trim(adjustl(filename)), NF90_WRITE, ncid), routine_name )
	call check( nf90_inq_dimid(ncid, 'copy', dimid), routine_name )
	call check( nf90_inquire_dimension(ncid, dimid, len = copy_out), routine_name )
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

SUBROUTINE read_smat(ncid,n,thB,thT,trB,trT)

  IMPLICIT NONE
	
	integer,intent(in)                                :: ncid, n
  real, dimension(2*kmax,2*lmax,tmax), intent (out) :: thB, thT, trB, trT
  real, dimension(2*kmax,2*lmax,tmax,1)             :: val
	integer                      :: varid, ierr, start(4), count(4)
	character(len=64), parameter :: routine_name = 'read_smat'
	
	val = 0.
  thB = 0.; thT = 0.
  trB = 0.; trT = 0.
	
	start(1) = 1      ; start(2) = 1      ; start(3) = 1    ; start(4) = n
	count(1) = 2*kmax ; count(2) = 2*lmax ; count(3) = tmax ; count(4) = 1

	call check( nf90_inq_varid(ncid, 'thetaB', varid), routine_name )
	call check( nf90_get_var(ncid, varid, val, start, count), routine_name )
	thB(:,:,:) = val(:,:,:,1);
	call check( nf90_inq_varid(ncid, 'thetaT', varid), routine_name )
	call check( nf90_get_var(ncid, varid, val, start, count), routine_name )
	thT(:,:,:) = val(:,:,:,1);
	
	ierr = nf90_inq_varid(ncid, 'tracerB', varid)
	if ( ierr == nf90_noerr ) then
		call check( nf90_inq_varid(ncid, 'tracerB', varid), routine_name )
		call check( nf90_get_var(ncid, varid, val, start, count), routine_name )
		trB(:,:,:) = val(:,:,:,1);
		call check( nf90_inq_varid(ncid, 'tracerT', varid), routine_name )
		call check( nf90_get_var(ncid, varid, val, start, count), routine_name )
		trT(:,:,:) = val(:,:,:,1);
	endif

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
	
		if(debug) print*,'creating file ... ', trim(adjustl(smatDfile))

		call check( nf90_create(trim(adjustl(smatDfile)), NF90_CLOBBER .or. NF90_64BIT_OFFSET, ncid), routine_name )
	
		call check( nf90_def_dim(ncid, "nx",           2*kmax, vardim(1)), routine_name )
		call check( nf90_def_dim(ncid, "ny",           2*lmax, vardim(2)), routine_name )
		call check( nf90_def_dim(ncid, "time",           tmax, vardim(3)), routine_name )
		call check( nf90_def_dim(ncid, "copy", NF90_UNLIMITED, vardim(4)), routine_name )
		call check( nf90_def_dim(ncid, "metadatalength",   64,     mddim), routine_name )

		call check( nf90_def_var(ncid, "copy", NF90_INT, vardim(4),     varid),  routine_name )
		call check( nf90_def_var(ncid, "CopyMetaData", NF90_CHAR, (/mddim, vardim(4)/), varid), routine_name )
		call check( nf90_def_var(ncid, "thetaB", NF90_FLOAT, vardim,    varid),  routine_name )
		call check( nf90_put_att(ncid, varid, "description", "bottom potential temperature"),   routine_name )
		call check( nf90_def_var(ncid, "thetaT", NF90_FLOAT, vardim,    varid),  routine_name )
		call check( nf90_put_att(ncid, varid, "description", "toppom potential temperature"),   routine_name )
		call check( nf90_def_var(ncid, "tracerB", NF90_FLOAT, vardim,    varid), routine_name )
		call check( nf90_put_att(ncid, varid, "description", "bottom tracer"),   routine_name )
		call check( nf90_def_var(ncid, "tracerT", NF90_FLOAT, vardim,    varid), routine_name )
		call check( nf90_put_att(ncid, varid, "description", "toppom tracer"),   routine_name )
		
		call check( nf90_put_att(ncid, NF90_GLOBAL, "title", 'Perturbation - Control'), routine_name )
		call check( nf90_put_att(ncid, NF90_GLOBAL, "ens_size", ens_size),              routine_name )

		call check( nf90_enddef(ncid), routine_name )
		call check( nf90_close(ncid),  routine_name )

		call transfer_ga(smatCfile,smatDfile)

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
	
	call check( nf90_open(trim(adjustl(smatDfile)), NF90_WRITE, ncid), routine_name )

	call check( nf90_inq_varid(ncid, "copy",          varid), routine_name )
	call check( nf90_put_var(ncid, varid, n,   (/start(4)/)), routine_name )
	call check( nf90_inq_varid(ncid, "CopyMetaData",  varid), routine_name )
	call check( nf90_put_var(ncid, varid, copy_string, (/start(1), start(4)/)), routine_name )
	call check( nf90_inq_varid(ncid, "thetaB",        varid), routine_name )
	call check( nf90_put_var(ncid, varid, thB, start, count), routine_name )
	call check( nf90_inq_varid(ncid, "thetaT",        varid), routine_name )
	call check( nf90_put_var(ncid, varid, thT, start, count), routine_name )
	call check( nf90_inq_varid(ncid, "tracerB",       varid), routine_name )
	call check( nf90_put_var(ncid, varid, trB, start, count), routine_name )
	call check( nf90_inq_varid(ncid, "tracerT",       varid), routine_name )
	call check( nf90_put_var(ncid, varid, trT, start, count), routine_name )
	
	call check( nf90_close(ncid), routine_name )
  
	return
end SUBROUTINE write_diag

SUBROUTINE check_smat()
	
	implicit none
	integer :: ncidC, ncidP, dimid
	integer :: ivalC, ivalP
	real    :: rvalC, rvalP
	character(len=64)            :: dim_name
	character(len=64), parameter :: routine_name = 'check_smat'
	
	call check( nf90_open(trim(adjustl(smatCfile)), NF90_NOWRITE, ncidC), routine_name )
	call check( nf90_open(trim(adjustl(smatPfile)), NF90_NOWRITE, ncidP), routine_name )

	call check( nf90_inq_dimid(ncidC,'nx',dimid), routine_name )  
	call check( nf90_inquire_dimension(ncidC,dimid,dim_name,ivalC), routine_name )  
	call check( nf90_inq_dimid(ncidP,'nx',dimid), routine_name )  
	call check( nf90_inquire_dimension(ncidP,dimid,dim_name,ivalP), routine_name )  
	call compare_integers(ivalC, ivalP, 'nx')
	
	call check( nf90_inq_dimid(ncidC,'ny',dimid), routine_name )  
	call check( nf90_inquire_dimension(ncidC,dimid,dim_name,ivalC), routine_name )  
	call check( nf90_inq_dimid(ncidP,'ny',dimid), routine_name )  
	call check( nf90_inquire_dimension(ncidP,dimid,dim_name,ivalP), routine_name )  
	call compare_integers(ivalC, ivalP, 'ny')
	
	call check( nf90_inq_dimid(ncidC,'time',dimid), routine_name )  
	call check( nf90_inquire_dimension(ncidC,dimid,dim_name,ivalC), routine_name )  
	call check( nf90_inq_dimid(ncidP,'time',dimid), routine_name )  
	call check( nf90_inquire_dimension(ncidP,dimid,dim_name,ivalP), routine_name )  
	call compare_integers(ivalC, ivalP, 'time')
	
	call check( nf90_inq_dimid(ncidC,'copy',dimid), routine_name )  
	call check( nf90_inquire_dimension(ncidC,dimid,dim_name,ivalC), routine_name )  
	call check( nf90_inq_dimid(ncidP,'copy',dimid), routine_name )  
	call check( nf90_inquire_dimension(ncidP,dimid,dim_name,ivalP), routine_name )  
	call compare_integers(ivalC, ivalP, 'copy')
	
	call check( nf90_get_att(ncidC,NF90_GLOBAL,'ens_size',ivalC), routine_name )  
	call check( nf90_get_att(ncidP,NF90_GLOBAL,'ens_size',ivalP), routine_name )
	call compare_integers(ivalC, ivalP, 'ens_size')
	
	call check( nf90_get_att(ncidC,NF90_GLOBAL,'model',ivalC), routine_name )  
	call check( nf90_get_att(ncidP,NF90_GLOBAL,'model',ivalP), routine_name )
	call compare_integers(ivalC, ivalP, 'model')
	
	call check( nf90_get_att(ncidC,NF90_GLOBAL,'dt',rvalC), routine_name )  
	call check( nf90_get_att(ncidP,NF90_GLOBAL,'dt',rvalP), routine_name )
	call compare_reals(rvalC, rvalP, 'dt')
	
	call check( nf90_get_att(ncidC,NF90_GLOBAL,'iplot',ivalC), routine_name )  
	call check( nf90_get_att(ncidP,NF90_GLOBAL,'iplot',ivalP), routine_name )
	call compare_integers(ivalC, ivalP, 'iplot')
	
	call check( nf90_get_att(ncidC,NF90_GLOBAL,'XL',rvalC), routine_name )  
	call check( nf90_get_att(ncidP,NF90_GLOBAL,'XL',rvalP), routine_name )
	call compare_reals(rvalC, rvalP, 'XL')
	
	call check( nf90_get_att(ncidC,NF90_GLOBAL,'YL',rvalC), routine_name )  
	call check( nf90_get_att(ncidP,NF90_GLOBAL,'YL',rvalP), routine_name )
	call compare_reals(rvalC, rvalP, 'YL')
	
	call check( nf90_get_att(ncidC,NF90_GLOBAL,'H',rvalC), routine_name )  
	call check( nf90_get_att(ncidP,NF90_GLOBAL,'H',rvalP), routine_name )
	call compare_reals(rvalC, rvalP, 'H')
	
	call check( nf90_get_att(ncidC,NF90_GLOBAL,'Ross',rvalC), routine_name )  
	call check( nf90_get_att(ncidP,NF90_GLOBAL,'Ross',rvalP), routine_name )
	call compare_reals(rvalC, rvalP, 'Ross')
	
	call check( nf90_get_att(ncidC,NF90_GLOBAL,'gamma',rvalC), routine_name )  
	call check( nf90_get_att(ncidP,NF90_GLOBAL,'gamma',rvalP), routine_name )
	call compare_reals(rvalC, rvalP, 'gamma')
	
	call check( nf90_get_att(ncidC,NF90_GLOBAL,'tau',rvalC), routine_name )  
	call check( nf90_get_att(ncidP,NF90_GLOBAL,'tau',rvalP), routine_name )
	call compare_reals(rvalC, rvalP, 'tau')
	
	call check( nf90_get_att(ncidC,NF90_GLOBAL,'trl',rvalC), routine_name )  
	call check( nf90_get_att(ncidP,NF90_GLOBAL,'trl',rvalP), routine_name )
	call compare_reals(rvalC, rvalP, 'trl')
	
	call check( nf90_get_att(ncidC,NF90_GLOBAL,'amu',rvalC), routine_name )  
	call check( nf90_get_att(ncidP,NF90_GLOBAL,'amu',rvalP), routine_name )
	call compare_reals(rvalC, rvalP, 'amu')
	
	call check( nf90_get_att(ncidC,NF90_GLOBAL,'shear',rvalC), routine_name )  
	call check( nf90_get_att(ncidP,NF90_GLOBAL,'shear',rvalP), routine_name )
	call compare_reals(rvalC, rvalP, 'shear')
	
	call check( nf90_close(ncidC), routine_name )
	call check( nf90_close(ncidP), routine_name )

	return
end SUBROUTINE check_smat 

SUBROUTINE compare_integers(ival1, ival2, varname)

	implicit none
	integer, intent(in)          :: ival1, ival2
	character(len=*), intent(in) :: varname
	character(len=64), parameter :: routine_name = 'compare_ints'

	if ( ival1 /= ival2 ) then
		write(6,'(a18, 1x, a20)') 'variable mismatch:', trim(adjustl(varname))
		write(6,'(a1, 1x, i6, 1x, a1, 1x, i6)') 'C', ival1, 'P', ival2
		stop
	endif

	return
end SUBROUTINE compare_integers

SUBROUTINE compare_reals(rval1, rval2, varname)

	implicit none
	real, intent(in)             :: rval1, rval2
	character(len=*), intent(in) :: varname
	character(len=64), parameter :: routine_name = 'compare_reals'

	if ( rval1 /= rval2 ) then
		write(6,'(a18, 1x, a20)') 'variable mismatch:', trim(adjustl(varname))
		write(6,'(a1, 1x, f6.6, 1x, a1, 1x, f6.6)') 'C', rval1, 'P', rval2
		stop
	endif

	return
end SUBROUTINE compare_reals

SUBROUTINE transfer_ga(fromfile,tofile)

	implicit none
	character(len=*), intent (in)  :: fromfile, tofile
	integer                        :: ncidf, ncidt
	character(len=64), parameter   :: routine_name = 'transfer_ga'

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
  
END PROGRAM sqg_stats
