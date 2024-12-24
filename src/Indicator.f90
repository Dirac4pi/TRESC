! |-------<THOMAS RELATIVISTIC ELECTRONIC STRUCTURE CALCULATIONS>------|
! |-----------------------------<DIRAC4PI>-----------------------------|
! |-----------------------------<INDICATOR>----------------------------|
! |--------------------------<PROCESS CONTROL>-------------------------|
    
module Indicator
    use Atoms
    use Hamiltonian
    use SCF
    contains
    
!-----------------------------------------------------------------------
! test of subroutine V_Integral_1e from module Hamiltonian
    
    subroutine V_Integral_1e_test()
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
        write(*,*) V_Integral_1e(Z,coe_i,coe_j,fac_i,fac_j,expo_i,expo_j,cod_i_inp,cod_j_inp,0.0_dp)
    end subroutine V_Integral_1e_test
    
!-----------------------------------------------------------------------
! test of subroutine V_Integral_2e from module Hamiltonian
    
    subroutine V_Integral_2e_test()
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
        
        write(*,*) 3.14124139088628 * 1.78055309907994 * 6.081100279080574E-002 * 1.78055309907994 *& 
            V_Integral_2e(fac_i,fac_j,fac_k,fac_l,ai,aj,ak,al,cod_i,cod_j,cod_k,cod_l)
    end subroutine V_Integral_2e_test
        
!-----------------------------------------------------------------------
! standard single-configuration calculation
! do NOT use universal loop variables
    
    subroutine sconf_calc(keep, kill)
        implicit none
        integer,intent(in) :: keep                                    ! initialization level
        logical,intent(in) :: kill                                    ! whether to kill the whole program after SCF process
        integer :: cmd_count, molstat
        ! get .xyz file address
        cmd_count = command_argument_count()
        if (cmd_count == 0) then
            write(*,*) 'TRESC error: no argument received'
            pause
            stop
        else if (cmd_count >= 2) then
            write(*,*) 'TRESC error: too many arguments'
            pause
            stop
        end if
        call get_command_argument(1, address_molecular, status=molstat)
        if (molstat > 0) then
            write(*,*) 'TRESC error: molecular adress retrieval fails'
            pause
            stop
        else if (molstat < 0) then
            write(*,*) 'TRESC error: molecular adress too long'
            pause
            stop
        end if
        ! get information from .gbs and .xyz files
        call generate_tot()
        call read_keywords()
        call read_gbs()
        call read_geometry()
        call input_check()
        call input_print()
        ! Hartree-Fock SCF process
        if (DKH_order == 0 .or. DKH_order == 2) then
            call DKH_Hamiltonian()
            if (molden) then
                call DKH_Hartree_Fock(keep = 1,kill = .false.)
                call dump_molden()
                pause
                stop
            else
                call DKH_Hartree_Fock(keep,kill)
            end if
        else if (DKH_order == 1) then
            write(60,*)
            write(60,'(A)') 'Module Hamiltonian: The SRTP effect requires at least a second-order'
            write(60,'(A)') '   DKH Hamiltonian to be accurately described, since the lead term of'
            write(60,'(A)') '   SRTP effect expanded at c=0 is in order c^-4. Only DKH2+ Hamiltonian'
            write(60,'(A)') '   guarantees the completeness of the expansion term at order c^-4.'
            call terminate('fpFW Hamiltonian is not supported')
        else
            call terminate('DKH_order is set incorrectly')
        end if
    end subroutine sconf_calc
    
    
end module Indicator