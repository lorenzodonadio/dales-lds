!> \file modgenstat_lite.f90
!!  Genstat calculates slab averages of several variables and saves them to NetCDF.

module modgenstat_lite

  use modprecision
  implicit none
  PUBLIC :: initgenstat_lite, writestat_lite, exitgenstat_lite
  save

  !NetCDF variables
  ! integer :: nvar = 10  ! Only saving 10 variables as specified
  integer :: nvar = 4  ! Only saving 4 variables as specified
  integer :: ncid, nrec = 0
  character(80) :: fname = 'profiles_lite.xxx.nc'  ! Updated file name
  character(80), allocatable, dimension(:,:) :: ncname
  character(80), dimension(1,4) :: tncname

  real :: timeav
  integer(kind=longint) ::  itimeav, tnextwrite
  logical :: lstat = .false.
contains  

  subroutine initgenstat_lite
    use modmpi, only : myid
    use modglobal, only : kmax, k1, ifnamopt, fname_options, dtav_glob, timeav_glob, btime, tres, checknamelisterror,cexpnr
    use modstat_nc, only : open_nc, define_nc, ncinfo, nctiminfo, writestat_dims_nc

    implicit none
    integer ierr

    
    namelist/NAMGENSTATLITE/ timeav,lstat
    
    timeav = timeav_glob

    if (myid == 0) then
      open(ifnamopt, file=fname_options, status='old', iostat=ierr)
      read(ifnamopt, NAMGENSTATLITE, iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'NAMGENSTATLITE')
      write(*, NAMGENSTATLITE)
      close(ifnamopt)
      
      itimeav = timeav / tres
      tnextwrite = itimeav + btime
      if(.not. lstat) return
      
      fname(15:17) = cexpnr  ! Update filename with experiment number
      allocate(ncname(nvar, 4))  ! Only 10 variables to allocate
      call nctiminfo(tncname(1,:))
      
      ! Define the variables to be written
      call ncinfo(ncname( 1,:),'rhof','Full level slab averaged density','kg/m^3','tt')
      call ncinfo(ncname( 2,:),'rhobf','Full level base-state density','kg/m^3','tt')
      call ncinfo(ncname( 3,:),'rhobh','Half level base-state density','kg/m^3','mt')
      call ncinfo(ncname( 4,:),'presh','Pressure at cell center','Pa','tt')
      ! call ncinfo(ncname( 5,:),'u','West-East velocity','m/s','tt')
      ! call ncinfo(ncname( 6,:),'v','South-North velocity','m/s','tt')
      ! call ncinfo(ncname( 7,:),'thl','Liquid water potential temperature','K','tt')
      ! call ncinfo(ncname( 8,:),'thv','Virtual potential temperature','K','tt')
      ! call ncinfo(ncname( 9,:),'qt','Total water specific humidity','kg/kg','tt')
      ! call ncinfo(ncname(10,:),'ql','Liquid water specific humidity','kg/kg','tt')
      
      ! Open the NetCDF file and write initial dimensions
      call open_nc(fname, ncid, nrec, n3=kmax)
      if (nrec == 0) then
        call define_nc(ncid, 1, tncname)
        call writestat_dims_nc(ncid)
      end if
      call define_nc(ncid, nvar, ncname)
      write(*,*) 'Initialized succesfully netcdf genstatlite', fname, ncid
    end if
  end subroutine initgenstat_lite

  subroutine writestat_lite
    use modglobal, only : kmax, k1, rk3step, timee,rtimee
    use modfields, only : presh, rhof, rhobf, rhobh
    use modmpi, only : myid
    use modstat_nc, only : writestat_nc
    implicit none

    real, dimension(k1, nvar) :: vars

    ! Assign values to the variables to be written
    if (myid /= 0) return
    if (.not. lstat) return
    if (rk3step/=3) return

      if (timee>=tnextwrite) then
        tnextwrite = tnextwrite+itimeav
        write(*,*) 'Attempt to write to', ncid
        
        vars(:, 1) = rhof
        vars(:, 2) = rhobf
        vars(:, 3) = rhobh
        vars(:, 4) = presh
        ! vars(:, 5) = umn
        ! vars(:, 6) = vmn
        ! vars(:, 7) = thlmn
        ! vars(:, 8) = thvmn
        ! vars(:, 9) = qtmn
        ! vars(:,10) = qlmn
        
        write(*,*) Time:", ncid, 1, tncname, (/rtimee/), nrec
        ! Write the time information to NetCDF
        call writestat_nc(ncid, 1, tncname, (/rtimee/), nrec, .true.)
        
        write(*,*) Vars:", ncid, nvar, ncname, vars(1:kmax,:), nrec, kmax
        ! Write the slab-averaged variables to NetCDF
        call writestat_nc(ncid, nvar, ncname, vars(1:kmax,:), nrec, kmax)
      end if
  end subroutine writestat_lite

  subroutine exitgenstat_lite
    use modmpi, only : myid
    use modstat_nc, only : exitstat_nc
    implicit none

    if (myid == 0) call exitstat_nc(ncid)
  end subroutine exitgenstat_lite

end module modgenstat_lite
