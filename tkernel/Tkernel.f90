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
  !===========================================================
  else if (cmd_count == 1) then
    call get_command_argument(1, cmd1, status=stat1)
    if (stat1 > 0) then
      call terminate('1st argument retrieval fails')
    else if (stat1 < 0) then
      call terminate('1st argument too long')
    end if
    if (cmd1(1:1) /= '-') then
      call Elec_Struc_Calc(.true.)
    else if ((trim(cmd1) == '-v' .or. trim(cmd1) == '-V')) then
      call Printversion()
    else if ((trim(cmd1) == '--version')) then
      call Printversion()
    else
      call terminate('unknown argument '//trim(cmd1))
    end if
    !do tri = 1, 10
    !  do batch calculation, plot with batch_plot.py
    !end do
    !===========================================================
  else if (cmd_count == 2) then
    call get_command_argument(1, cmd1, status=stat1)
    if (stat1 > 0) then
      call terminate('1st argument retrieval fails')
    else if (stat1 < 0) then
      call terminate('1st argument too long')
    end if
    if (trim(cmd1) == '-moana' .or. trim(cmd1) == '--moana') then
      call MOANA(.true.)
    else if (trim(cmd1) == '-MOANA' .or. trim(cmd1) == '--MOANA') then
      call MOANA(.true.)
    else
      call terminate('unknown argument '//trim(cmd1))
    end if
  else if (cmd_count == 3) then
    call get_command_argument(1, cmd1, status=stat1)
    if (stat1 > 0) then
      call terminate('1st argument retrieval fails')
    else if (stat1 < 0) then
      call terminate('1st argument too long')
    end if

    call lowercase(cmd1)
    if (trim(cmd1) == '-cub2c') then
      call Grid2c('c', .true.)
    else if (trim(cmd1) == '-mog2c') then
      call Grid2c('m', .true.)
    else if (trim(cmd1) == '-pro2c') then
      call Grid2c('p', .true.)
    else if (trim(cmd1) == '-cub1c') then
      call Grid1c('c', .true.)
    else if (trim(cmd1) == '-mog1c') then
      call Grid1c('m', .true.)
    else if (trim(cmd1) == '-pro1c') then
      call Grid1c('p', .true.)
    else
      call terminate('unknown argument '//trim(cmd1))
    end if
  else
    call terminate('invalid arguments')
  end if
end program Tkernel
  
!-----------------------------------------------------------------------
!> print current version of TRESC
subroutine Printversion()
  use Fundamentals
  use omp_lib, only: omp_get_max_threads
  use xc_f03_lib_m, only: xc_f03_version
  implicit none
  character(len=200) :: env_TRESC = ''
  integer            :: max_threads
  integer            :: vmajor, vminor, vmicro
  write(*,*) '***Thomas Relativistic Electronic Structure Calculation***'
  write(*,*) '***         https://github.com/dirac4pi/TRESC          ***'
  write(*,*) 'tkernel version: '//trim(version)
  call getenv('TRESC',env_TRESC)
  if (env_TRESC == '') then
    write(*,*) '$TRESC = None'
  else
    write(*,*) '$TRESC = '//trim(env_TRESC)
  end if
  write(*,*) 'check dynamic link library compatibility'
  write(*,*) 'OpenMP checking ...'
  max_threads = omp_get_max_threads()
  write(*,'(A,I0)') &
  " [checked] OpenMP link succeeded, OpenMP max threads = ",max_threads
  write(*,*) 'LibXC checking ...'
  call xc_f03_version(vmajor, vminor, vmicro)
  write(*,'(A,I0,A,I0,A,I0)') &
  " [checked] LibXC link succeeded, LibXC version = ",&
  vmajor, ".", vminor, ".", vmicro
end subroutine Printversion

!-----------------------------------------------------------------------
!> standard electronic structure calculation
subroutine Elec_Struc_Calc(kill_)
  use SCF
  use Fundamentals
  implicit none
  logical,intent(in) :: kill_  ! whether to kill the process after SCF
  integer            :: trestat
  call get_command_argument(1, address_tre, status=trestat)
  if (trestat > 0) call terminate('input file (.tre) retrieval fails')
  if (trestat < 0) call terminate('input file name (.tre) too long')
  call generate_output()
  call Load_keywords_tre()
  call Load_basis_gbs()
  call Load_geom_xyz()
  call Check_all_loads()
  if (pVp1e) then
    call Hess_Hamiltonian()
  else
    call Hartree_Hamiltonian()
  end if
  if (.not.pVp1e) then
    call SCF_V2e()
  else if (pVp1e .and. .not.pVp2e) then
    call SCF_AVA2e()
  else if (pVp1e .and. pVp2e) then
    call SCF_ARVRA2e()
  end if
  call Globinit(kill_)
end subroutine Elec_Struc_Calc

!-----------------------------------------------------------------------
!> dump grid value of 2-component spherical molecular orbitals for visualization
!!
!! cube = 'c': generate cube files for real space
!! cube = 'p': generate cube files for both real space and momentum space
!! cube = 'm': generate mog files for real space
subroutine Grid2c(cube, kill_)
  use Representation
  use Fundamentals
  use Hamiltonian
  implicit none
  character(len=1),intent(in)      :: cube
  logical,intent(in)               :: kill_
  complex(dp),allocatable          :: arr(:)
  complex(dp),allocatable          :: arr2(:)
  real(dp),allocatable             :: supparr(:)
  character(len=60)                :: cmd2, cmd3
  integer                          :: imo
  integer                          :: ii
  call get_command_argument(2, cmd2, status=ios)
  if (ios > 0) then
    call terminate('2nd argument retrieval fails')
  else if (ios < 0) then
    call terminate('2nd argument too long')
  end if
  jobname = trim(cmd2(:index(cmd2,'.molden.d')-1))
  address_MOLDEN = trim(cmd2)//'/'//trim(jobname)//'-realpart1.molden.input'
  cbdm = 5000
  allocate(cbdata(3*cbdm))
  call Load_geombasis_MOLDEN(.false.)
  call Normalize_Molden_Basis(M_cbdata, M_cbdm, isorca)
  cbdm = M_cbdm
  sbdm = M_sbdm
  cbdata = M_cbdata
  mol = M_mol
  atom_count = M_atom_count
  allocate(arr(2*sbdm), arr2(2*cbdm), supparr(2*sbdm))
  call Assign_csf()
  call get_command_argument(3, cmd3, status=ios)
  read(cmd3,'(I)',iostat=ios) imo
  if (ios /= 0) call terminate('unable to recognize orbital index')
  if (imo <= sbdm) then
    address_MOLDEN = trim(cmd2)//'/'//trim(jobname)//'-realpart1.molden.input'
  else
    address_MOLDEN = trim(cmd2)//'/'//trim(jobname)//'-realpart2.molden.input'
  end if
  call Load_MO_MOLDEN(imo, supparr(1:sbdm), imo, supparr(sbdm+1:2*sbdm))
  arr = supparr
  if (imo <= sbdm) then
    address_MOLDEN = trim(cmd2)//'/'//trim(jobname)//'-imgpart1.molden.input'
  else
    address_MOLDEN = trim(cmd2)//'/'//trim(jobname)//'-imgpart2.molden.input'
  end if
  call Load_MO_MOLDEN(imo, supparr(1:sbdm), imo, supparr(sbdm+1:2*sbdm))
  arr = arr + ci * supparr
  call matmul('N', exc2s, arr, arr2)
  if (cube == 'c' .or. cube == 'C') then
    write(*,*) ' generating real-space cube files...'
    call Basis2real_cube(.true., arr2, cmd3)
  else if (cube == 'p' .or. cube == 'P') then
    write(*,*) ' generating real-space cube files...'
    call Basis2real_cube(.true., arr2, cmd3)
    write(*,*) ' generating momentum-space cube files...'
    call Basis2momentum_cube(.true., arr2, cmd3)
  else if (cube == 'm' .or. cube == 'M') then
    write(*,*) ' generating real-space mog files...'
    call Basis2real_Becke_mog(.true., arr2, cmd3)
  else
    call terminate('Grid2c: unknown grid type')
  end if
  deallocate(arr, arr2, supparr)
  if (kill_) call terminate('normal')
end subroutine Grid2c

!-----------------------------------------------------------------------
!> dump grid value of scalar spherical molecular orbitals for visualization
!!
!! cube = 'c': generate cube files for real space
!! cube = 'p': generate cube files for both real space and momentum space
!! cube = 'm': generate mog files for real space
subroutine Grid1c(cube, kill_)
  use Representation
  use Fundamentals
  use Hamiltonian
  implicit none
  character(len=1),intent(in)      :: cube
  logical,intent(in)               :: kill_
  integer                          :: spin
  real(dp),allocatable             :: arr(:)
  complex(dp),allocatable          :: arr2(:)
  character(len=60)                :: cmd2, cmd3, fimo
  integer                          :: imo
  integer                          :: ii
  call get_command_argument(2, cmd2, status=ios)
  if (ios > 0) then
    call terminate('2nd argument retrieval fails')
  else if (ios < 0) then
    call terminate('2nd argument too long')
  end if
  address_MOLDEN = trim(cmd2)
  cbdm = 5000
  allocate(cbdata(3*cbdm))
  call Load_geombasis_MOLDEN(.false.)
  call Normalize_Molden_Basis(M_cbdata, M_cbdm, isorca)
  cbdm = M_cbdm
  sbdm = M_sbdm
  cbdata = M_cbdata
  mol = M_mol
  atom_count = M_atom_count
  allocate(arr(2*sbdm), arr2(2*cbdm))
  call Assign_csf()
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
    call Load_MO_MOLDEN(imo, arr(1:sbdm), 0, arr(sbdm+1:2*sbdm))
  else if(spin == 2) then  ! beta
    call Load_MO_MOLDEN(0, arr(1:sbdm), imo, arr(sbdm+1:2*sbdm))
  end if
  call matmul('N', exc2s, arr*c1, arr2)
  if (cube == 'c' .or. cube == 'C') then
    write(*,*) ' generating real-space cube file...'
    call Basis2real_cube(.false., arr2, cmd3, spin)
  else if (cube == 'p' .or. cube == 'P') then
    write(*,*) ' generating real-space cube file...'
    call Basis2real_cube(.false., arr2, cmd3, spin)
    write(*,*) ' generating momentum-space cube file...'
    call Basis2momentum_cube(.false., arr2, cmd3, spin)
  else if (cube == 'm' .or. cube == 'M') then
    write(*,*) ' generating real-space mog file...'
    call Basis2real_Becke_mog(.false., arr2, cmd3)
  else
    call terminate('Grid1c: unknown grid type')
  end if
  deallocate(arr, arr2)
  if (kill_) call terminate('normal')
end subroutine Grid1c

!-----------------------------------------------------------------------
!> Molecular Orbital ANAlysis
!!
!! spin pure states decomposition and TRS-breaking analysis for
!! molden input files
subroutine MOANA(kill_)
  use Fundamentals
  use Atoms
  use Hamiltonian
  use SCF
  implicit none
  logical,intent(in)        :: kill_
  real(dp),allocatable      :: moa(:,:), mob(:,:)
  real(dp),allocatable      :: enea(:), eneb(:)
  real(dp),allocatable      :: occene(:)
  real(dp),allocatable      :: sortedoccene(:)
  real(dp),allocatable      :: occcoe(:,:)
  integer,allocatable       :: sortocc(:)
  character(len=60)         :: cmd2
  integer                   :: MoldenNa, MOldenNb
  integer                   :: ii, jj, kk, ll
  integer                   :: contri      ! contr of atoMi, shelLi
  integer                   :: faci(3)     ! xyz factor of |AOi>
  integer                   :: contrj      ! contr of atoMj, shelLj
  integer                   :: facj(3)     ! xyz factor of |AOj>
  real(dp)                  :: expi(16)    ! expo of |AOi>
  real(dp)                  :: expj(16)    ! expo of |AOj>
  real(dp)                  :: coei(16)    ! coefficient of |AOi>
  real(dp)                  :: coej(16)    ! coefficient of |AOj>
  real(dp)                  :: codi(3)     ! coordinate of center of |AOi>
  real(dp)                  :: codj(3)     ! coordinate of center of |AOj>

  call get_command_argument(2, cmd2, status=ios)
  if (ios > 0) then
    call terminate('2nd argument retrieval fails')
  else if (ios < 0) then
    call terminate('2nd argument too long')
  end if
  address_MOLDEN = trim(cmd2)
  cbdm = 5000
  allocate(cbdata(3*cbdm))
  call Load_geombasis_MOLDEN(.false.)
  call Normalize_Molden_Basis(M_cbdata, M_cbdm, isorca)
  cbdm = M_cbdm
  sbdm = M_sbdm
  cbdata = M_cbdata
  mol = M_mol
  atom_count = M_atom_count
  !-------------<basis overlap integral>-------------
  allocate(i_j(cbdm, cbdm), source = 0.0_dp)
  do ii = 1, cbdm
    contri = cbdata(ii) % contr
    faci   = cbdata(ii) % fac
    expi(1:contri) = cbdata(ii) % expo(1:contri)
    coei(1:contri) = cbdata(ii) % Ncoe(1:contri)
    codi   = cbdata(ii) % pos
    do jj = 1, cbdm
      contrj = cbdata(jj) % contr
      facj   = cbdata(jj) % fac
      expj(1:contrj) = cbdata(jj) % expo(1:contrj)
      coej(1:contrj) = cbdata(jj) % Ncoe(1:contrj)
      codj   = cbdata(jj) % pos
      do kk = 1, contri
        do ll = 1, contrj
          i_j(ii,jj) = i_j(ii,jj) + Integral_S_1e(&
          coei(kk)*coej(ll),faci,facj,expi(kk),expj(ll),codi,codj)
        end do
      end do
    end do
  end do
  call Assign_csf(i_j)
  allocate(moa(M_sbdm,M_sbdm), mob(M_sbdm, M_sbdm))
  allocate(enea(M_sbdm), eneb(M_sbdm))
  call Load_MO_MOLDEN(moa, mob, MoldenNa, MoldenNb, enea, eneb)
  electron_count = MoldenNa + MoldenNb
  allocate(occindex(electron_count), occene(electron_count))
  allocate(sortedoccene(electron_count))
  allocate(sortocc(electron_count))
  allocate(occcoe(2*M_sbdm, electron_count), source=0.0_dp)
  forall (ii=1:electron_count) occindex(ii) = ii
  allocate(AO2MO(2*sbdm,electron_count), oper3(2*fbdm,electron_count), &
           source=c0)
  occene(1:MoldenNa) = enea(1:MoldenNa)
  occcoe(1:M_sbdm, 1:MoldenNa) = moa(1:M_sbdm, 1:MoldenNa)
  occene(MoldenNa+1:electron_count) = eneb(1:MoldenNb)
  occcoe(M_sbdm+1:2*M_sbdm, MoldenNa+1:electron_count) = mob(1:M_sbdm, 1:MoldenNb)
  call sort(occene, sortedoccene, sortocc)
  forall (ii=1:electron_count) AO2MO(1:2*sbdm, sortocc(ii)) = &
      occcoe(1:2*sbdm, sortocc(ii)) * c1
  address_tre = 'MOANA.tre'
  call generate_output()
  write(60,*)
  write(60,'(A)') "  Note: MOANA does not read unoccupied orbitals;"
  write(60,'(A)') "        if you need to select for active space,"
  write(60,'(A)') "        better choose RO orbitals, natural orbitals,"
  write(60,'(A)') "        or biorthogonalized orbitals, etc."
  write(60,*)
  call matmul('N', 'N', exf2s, AO2MO, oper3)
  call Calc_S2HF()
  call Kramers(moana=.true.)
  write(60,*)
  write(60,*)
  call RGI()
  deallocate(i_j, AO2MO, oper3, moa, mob, enea, eneb, occindex)
  deallocate(occene, occcoe, sortocc, sortedoccene)
  if (kill_) call terminate('normal')
end subroutine MOANA

