!> @file Tkernel.f90
!!
!! @brief TRESC kernel program
!!
!! @syntax Fortran 2008 free format
!!
!! @code UTF-8
!!
!! @author dirac4pi
program Tkernel
  use Fundamentals
  implicit none
  integer            :: tri, trj, trk
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
    !do tri = 1, 10
    !  do batch calculation, plot with batch_plot.py
    !end do
!------------------------------------------------------------------------------
  else if (cmd_count == 3) then
    call get_command_argument(1, cmd1, status=stat1)
    if (stat1 > 0) then
      call terminate('1st argument retrieval fails')
    else if (stat1 < 0) then
      call terminate('1st argument too long')
    end if

    call lowercase(cmd1)
    if (trim(cmd1) == '-cub2c') then
      call mo2cgrid(.true., .true.)
    else if (trim(cmd1) == '-mog2c') then
      call mo2cgrid(.false., .true.)
    else if (trim(cmd1) == '-cub1c') then
      call mo1cgrid(.true., .true.)
    else if (trim(cmd1) == '-mog1c') then
      call mo1cgrid(.false., .true.)
    else
      call terminate('unknown argument '//trim(cmd1))
    end if
  else
    call terminate('invalid arguments')
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
  call load_b_gbs()
  call load_g_xyz()
  call input_check()
  call input_print()
  ! SCF process
  call DKH_Hamiltonian()
  call DKH_SCF()
  call info_print()
  call glob_init(kill_)
end subroutine sconf_calc

!-----------------------------------------------------------------------
!> dump grid value of 2-component spherical molecular orbitals for visualization
subroutine mo2cgrid(cube, kill_)
  use Representation
  use Fundamentals
  implicit none
  logical,intent(in)      :: cube
  logical,intent(in)      :: kill_
  complex(dp),allocatable :: arr(:)
  complex(dp),allocatable :: arr2(:)
  real(dp),allocatable    :: supparr(:)
  character(len=60)       :: cmd2, cmd3
  integer                 :: imo
  call get_command_argument(2, cmd2, status=ios)
  if (ios > 0) then
    call terminate('2nd argument retrieval fails')
  else if (ios < 0) then
    call terminate('2nd argument too long')
  end if
  jobname = trim(cmd2(:index(cmd2,'.molden.d')-1))
  address_molden = trim(cmd2)//'/'//trim(jobname)//'-realpart1.molden.input'
  call load_gb_molden(tomol = .false.)
  cbdm = m_cbdm
  sbdm = m_sbdm
  basis_inf = m_basis_inf
  mol = m_mol
  atom_count = m_atom_count
  allocate(arr(2*sbdm), arr2(2*cbdm), supparr(2*sbdm))
  call assign_csf()
  call get_command_argument(3, cmd3, status=ios)
  read(cmd3,'(I)',iostat=ios) imo
  if (ios /= 0) call terminate('unable to recognize orbital index')
  if (imo <= sbdm) then
    address_molden = trim(cmd2)//'/'//trim(jobname)//'-realpart1.molden.input'
  else
    address_molden = trim(cmd2)//'/'//trim(jobname)//'-realpart2.molden.input'
  end if
  call load_MO_molden(imo, supparr(1:sbdm), imo, supparr(sbdm+1:2*sbdm))
  arr = supparr
  if (imo <= sbdm) then
    address_molden = trim(cmd2)//'/'//trim(jobname)//'-imgpart1.molden.input'
  else
    address_molden = trim(cmd2)//'/'//trim(jobname)//'-imgpart2.molden.input'
  end if
  call load_MO_molden(imo, supparr(1:sbdm), imo, supparr(sbdm+1:2*sbdm))
  arr = arr + ci * supparr
  call matmul('N', exc2s, arr, arr2)
  if (cube) then
    call dump_cube(.true., arr2, cmd3)
  else
    call mogrid_becke(.true., arr2, cmd3)
  end if
  deallocate(arr, arr2, supparr)
  if (kill_) call terminate('normal')
end subroutine mo2cgrid

!-----------------------------------------------------------------------
!> dump grid value of scalar spherical molecular orbitals for visualization
subroutine mo1cgrid(cube, kill_)
  use Representation
  use Fundamentals
  implicit none
  logical,intent(in)      :: cube
  logical,intent(in)      :: kill_
  integer                 :: spin
  real(dp),allocatable    :: arr(:)
  complex(dp),allocatable :: arr2(:)
  character(len=60)       :: cmd2, cmd3, fimo
  integer                 :: imo
  integer                 :: ii
  call get_command_argument(2, cmd2, status=ios)
  if (ios > 0) then
    call terminate('2nd argument retrieval fails')
  else if (ios < 0) then
    call terminate('2nd argument too long')
  end if
  address_molden = trim(cmd2)
  call load_gb_molden(tomol = .false.)
  cbdm = m_cbdm
  sbdm = m_sbdm
  basis_inf = m_basis_inf
  mol = m_mol
  atom_count = m_atom_count
  allocate(arr(2*sbdm), arr2(2*cbdm))
  call assign_csf()
  call get_command_argument(3, cmd3, status=ios)
  do ii = 1, len(cmd3)
    if (is_alpha(cmd3(ii:ii))) then
      if (cmd3(ii:ii)=='A' .or. cmd3(ii:ii)=='a') then
        spin = 1
      else if(cmd3(ii:ii)=='B' .or. cmd3(ii:ii)=='b') then
        spin = 2
      else
        call terminate('unable to recognize spin')
      end if
      exit
    end if
  end do
  if (ii == len(cmd3)) call terminate('no spin info in orbital index')
  fimo = cmd3(:ii-1)//cmd3(ii+1:)
  read(fimo,'(I)',iostat=ios) imo
  if (ios /= 0) call terminate('unable to recognize orbital index')
  if (spin == 1) then      ! alpha
    call load_MO_molden(imo, arr(1:sbdm), 0, arr(sbdm+1:2*sbdm))
  else if(spin == 2) then  ! beta
    call load_MO_molden(0, arr(1:sbdm), imo, arr(sbdm+1:2*sbdm))
  end if
  call matmul('N', exc2s, arr*c1, arr2)
  if (cube) then
    call dump_cube(.false., arr2, cmd3, spin)
  else
    call mogrid_becke(.false., arr2, cmd3)
  end if
  deallocate(arr, arr2)
  if (kill_) call terminate('normal')
end subroutine mo1cgrid