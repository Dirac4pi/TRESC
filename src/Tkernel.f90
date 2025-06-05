!> @file Tkernel.f90
!!
!! @brief TRESC startup procedure
!!
!! @syntax Fortran 2008 free format
!!
!! @code UTF-8
!!
!! @author dirac4pi
program Tkernel
  use Fundamentals
  implicit none
  integer            :: trloop_i, trloop_j, trloop_k
  integer            :: cmd_count, stat1
  character(len=200) :: cmd1
  cmd_count = command_argument_count()
  if (cmd_count == 0) then
    call terminate('0 argument received')
!------------------------------------------------------------------------------
  else if (cmd_count == 1) then
    call get_command_argument(1, cmd1, status=stat1)
    if (stat1 > 0) then
      call terminate('1st argument retrieval fails')
    else if (stat1 < 0) then
      call terminate('1st argument too long')
    end if

    if (cmd1(1:1) /= '-') then
      call sconf_calc(.true.)
    else if ((trim(cmd1) == '-v' .or. trim(cmd1) == '-V')) then
      call printversion()
    else if ((trim(cmd1) == '--version')) then
      call printversion()
    else
      call terminate('unknown argument '//trim(cmd1))
    end if
    !do trloop_i = 1, 10
    !  do batch calculation, plot with batch_plot.py
    !end do
!------------------------------------------------------------------------------
  else if (cmd_count == 2) then
    call get_command_argument(1, cmd1, status=stat1)
    if (stat1 > 0) then
      call terminate('1st argument retrieval fails')
    else if (stat1 < 0) then
      call terminate('1st argument too long')
    end if

    if ((trim(cmd1) == '-2c' .or. trim(cmd1) == '-2C')) then
      call mo2cgrid(.true.)
    else if ((trim(cmd1) == '-1c' .or. trim(cmd1) == '-1C')) then
      call mo1cgrid(.true.)
    else
      call terminate('unknown argument '//trim(cmd1))
    end if
  else
    call terminate('too many arguments')
  end if
end program Tkernel
  
!-----------------------------------------------------------------------
!> print current version of TRESC
subroutine printversion()
  use Fundamentals
  implicit none
  character(len=200) :: env_TRESC = ''
  write(*,*) 'Thomas Relativistic Electronic Structure Calculation'
  write(*,*) 'version: '//trim(version)
  call getenv('TRESC',env_TRESC)
  if (env_TRESC == '') then
    write(*,*) '$TRESC = None'
  else
    write(*,*) '$TRESC = '//trim(env_TRESC)
  end if
end subroutine printversion
    
!-----------------------------------------------------------------------
!> standard single-configuration calculation
subroutine sconf_calc(kill_)
  use SCF
  use Fundamentals
  implicit none
  logical,intent(in) :: kill_  ! whether to kill the process after SCF
  integer            :: molstat
  call get_command_argument(1, address_molecule, status=molstat)
  if (molstat > 0) then
    call terminate('molecule adress retrieval fails')
  else if (molstat < 0) then
    call terminate('molecule adress too long')
  end if
  ! get information from .gbs and .xyz files
  call generate_output()
  call read_keywords()
  call read_gbs()
  call read_geometry()
  call input_check()
  call input_print()
  ! SCF process
  if (DKH_order == 0 .or. DKH_order == 2) then
    call DKH_Hamiltonian()
    call DKH_SCF()
    call outprint()
    call glob_init(kill_)
  else if (DKH_order == 1) then
    write(60,*)
    write(60,'(A)') 'Module Hamiltonian:'
    write(60,'(A)') '  The SRTP effect requires at least a second-order DKH'
    write(60,'(A)') '  Hamiltonian to be accurately described, since the lead'
    write(60,'(A)') '  term of SRTP effect expanded at c=0 is in order c^-4.'
    write(60,'(A)') '  Only DKH2+ Hamiltonian guarantees the completeness of'
    write(60,'(A)') '  the expansion term at order c^-4.'
    call terminate('fpFW Hamiltonian is not supported')
  else
    call terminate('DKH_order is set incorrectly')
  end if
end subroutine sconf_calc

!-----------------------------------------------------------------------
!> dump grid value of 2-component molecular orbitals for visualization
subroutine mo2cgrid(kill_)
  use Representation
  use LAPACK95
  use Fundamentals
  implicit none
  logical,intent(in)      :: kill_
  integer                 :: channel, count
  complex(dp),allocatable :: arr(:), supparr(:)
  character(len=30)       :: cmd2
  integer                 :: stat2
  call get_command_argument(2, cmd2, status=stat2)
  if (stat2 > 0) then
    call terminate('2nd argument retrieval fails')
  else if (stat2 < 0) then
    call terminate('2nd argument too long')
  end if

  address_basis = '.gbs'
  call read_gbs()
  address_molecule = '.xyz'
  call read_geometry()
  call assign_cs(.true.)
  open(newunit=channel, file=trim(cmd2)//'.mo', action='read', iostat=ios)
  if (ios /= 0) call terminate('error when opening .mo file')
  allocate(supparr(2*cbdm))
  read(channel,*,iostat=ios) supparr(1)
  count = 1
  do while(ios == 0)
    count = count + 1
    read(channel,*,iostat=ios) supparr(count)
  end do
  count = count - 1
  if (count == 2*cbdm) then
    call mogrid_becke(.true., supparr, cmd2)
  else if (count == 2*sbdm) then
    allocate(arr(2*sbdm))
    arr(:) = supparr(1:2*sbdm)
    call zgemm( 'N', 'N', 2*cbdm, 1, 2*sbdm, c1, &
    exc2s, 2*cbdm, arr, 2*sbdm, c0, supparr, 2*cbdm)
    call mogrid_becke(.true., supparr, cmd2)
  end if
  close(channel)
  if (kill_) call terminate('normal')
end subroutine mo2cgrid

!-----------------------------------------------------------------------
!> dump grid value of scalar molecular orbitals for visualization
subroutine mo1cgrid(kill_)
  use Representation
  use LAPACK95
  use Fundamentals
  implicit none
  logical,intent(in)      :: kill_
  integer                 :: channel, count
  real(dp),allocatable    :: arr(:), supparr(:)
  character(len=30)       :: cmd2
  integer                 :: stat2
  call get_command_argument(2, cmd2, status=stat2)
  if (stat2 > 0) then
    call terminate('2nd argument retrieval fails')
  else if (stat2 < 0) then
    call terminate('2nd argument too long')
  end if

  address_basis = '.gbs'
  call read_gbs()
  address_molecule = '.xyz'
  call read_geometry()
  call assign_cs(.false.)
  open(newunit=channel, file=trim(cmd2)//'.mo', action='read', iostat=ios)
  if (ios /= 0) call terminate('error when opening .mo file')
  allocate(supparr(cbdm))
  read(channel,*,iostat=ios) supparr(1)
  count = 1
  do while(ios == 0)
    count = count + 1
    read(channel,*,iostat=ios) supparr(count)
  end do
  count = count - 1
  if (count == cbdm) then
    call mogrid_becke(.false., supparr*c1, cmd2)
  else if (count == sbdm) then
    allocate(arr(sbdm))
    arr(:) = supparr(1:sbdm)
    call dgemm( 'N', 'N', cbdm, 1, sbdm, 1.0_dp, &
    c2s, cbdm, arr, sbdm, 0.0_dp, supparr, cbdm)
    call mogrid_becke(.false., supparr*c1, cmd2)
  end if
  close(channel)
  if (kill_) call terminate('normal')
end subroutine mo1cgrid