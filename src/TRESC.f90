!> @file TRESC.f90
!!
!! @brief TRESC default startup procedure
!!
!! @author dirac4pi
program TRESC
  implicit none
  integer :: keep
  logical :: kill
  integer :: trloop_i, trloop_j, trloop_k
  call sconf_calc(0, .true.)
  !do trloop_i = 1, 10
  !    do batch calculation, plot with batch_plot.py
  !end do
end program TRESC
  
  
!-----------------------------------------------------------------------
!> test of subroutine V_Integral_1e from module Hamiltonian
subroutine V_Integral_1e_test()
  use Fundamentals
  use Hamiltonian
  implicit none
  real(dp) :: Z = 6.0
  real(dp) :: coe_i = 1.00000000
  real(dp) :: coe_j = -0.00014450
  real(dp) :: expo_i = 0.76100000
  real(dp) :: expo_j = 8.00000000
  character(len = 3) :: fac_i = 'zzx'
  character(len = 3) :: fac_j = 'yyz'
  real(dp) :: cod_i_inp(3), cod_j_inp(3)
  cod_i_inp(1) = 0.0256
  cod_i_inp(2) = -3.68
  cod_i_inp(3) = 5.395
  cod_j_inp(1) = 2.46
  cod_j_inp(2) = -4.198
  cod_j_inp(3) = 3.59
  write(*,*) V_Integral_1e(&
  Z,coe_i,coe_j,fac_i,fac_j,expo_i,expo_j,cod_i_inp(1),cod_i_inp(2),&
  cod_i_inp(3),cod_j_inp(1),cod_j_inp(2),cod_j_inp(3),0.0_dp)
end subroutine V_Integral_1e_test
  
!-----------------------------------------------------------------------
!> test of subroutine V_Integral_2e from module Hamiltonian
subroutine V_Integral_2e_test()
  use Fundamentals
  use Hamiltonian
  implicit none
  real(dp) :: ai = 1.05700000000000
  real(dp) :: aj = 0.761000000000000
  real(dp) :: ak = 33.8700000000000
  real(dp) :: al = 0.761000000000000
  character(len = 3) :: fac_i = 'xz '
  character(len = 3) :: fac_j = 'yzz'
  character(len = 3) :: fac_k = '   '
  character(len = 3) :: fac_l = 'yzz'
  real(dp) :: cod_i(3), cod_j(3), cod_k(3), cod_l(3)
  cod_i(1) = -9.05312365384779
  cod_i(2) = 1.79901241370262
  cod_i(3) = -1.65096194451096
    
  cod_j(1) = -9.72713758901604
  cod_j(2) = 0.845838045467446
  cod_j(3) = 0.000000000000000E+000
    
  cod_k(1) = -9.05312365384779
  cod_k(2) = 1.79901241370262
  cod_k(3) = 1.65096194451096
    
  cod_l(1) = -9.72713758901604
  cod_l(2) = 0.845838045467446
  cod_l(3) = 0.000000000000000E+000
    
  write(*,*) &
  3.141241390*1.7805530990*6.081100279-002*1.7805530990*V_Integral_2e(&
  fac_i,fac_j,fac_k,fac_l,ai,aj,ak,al,cod_i(1),cod_i(2),cod_i(3),cod_j(1),&
  cod_j(2),cod_j(3),cod_k(1),cod_k(2),cod_k(3),cod_l(1),cod_l(2),cod_l(3))
end subroutine V_Integral_2e_test
    
!-----------------------------------------------------------------------
!> standard single-configuration calculation
subroutine sconf_calc(keep_, kill_)
  use SCF
  use Fundamentals
  implicit none
  integer,intent(in) :: keep_  ! initialization level after SCF
  logical,intent(in) :: kill_  ! whether to kill the process after SCF
  integer :: cmd_count, molstat
  ! get .xyz file address
  cmd_count = command_argument_count()
  if (cmd_count == 0) then
    write(*,*) 'TRESC error: no argument received'
    stop 1
  else if (cmd_count >= 2) then
    write(*,*) 'TRESC error: too many arguments'
    stop 1
  end if
  call get_command_argument(1, address_molecular, status=molstat)
  if (molstat > 0) then
    write(*,*) 'TRESC error: molecular adress retrieval fails'
    stop 1
  else if (molstat < 0) then
    write(*,*) 'TRESC error: molecular adress too long'
    stop 1
  end if
  ! get information from .gbs and .xyz files
  call generate_output()
  call read_keywords()
  call read_gbs()
  call read_geometry()
  call input_check()
  call input_print()
  ! Hartree-Fock SCF process
  if (DKH_order == 0 .or. DKH_order == 2) then
    call DKH_Hamiltonian()
    call DKH_Hartree_Fock(keep_,kill_)
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