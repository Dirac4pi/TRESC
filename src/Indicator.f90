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
        
        write(*,*) 3.14124139088628 * 1.78055309907994 * 6.081100279080574E-002 * 1.78055309907994 * V_Integral_2e(fac_i,fac_j,fac_k,fac_l,ai,aj,ak,al,cod_i,cod_j,cod_k,cod_l)
    end subroutine V_Integral_2e_test
        
        
!-----------------------------------------------------------------------
! standard input process control
! DO NOT use universal loop variables
    
    subroutine elec_stc_calc()
        implicit none
        integer :: cmd_count, molstat, basstat
        character(len = 50) :: address_molecular = 'E:\TRESC\examples\carbene-3et\carbene-3et.tip'
        character(len = 50) :: address_basis = 'E:\TRESC\examples\carbene-3et\dkh-def2-svp-o.gbs'
        cmd_count = command_argument_count()
        if (cmd_count /= 2 .and. cmd_count /= 0) then
            write(*,*) 'Error: number of command line arguments is greater or less than 2'
            stop
        end if
        if (cmd_count == 2) then
        call get_command_argument(1, address_molecular, status=molstat)
            if (molstat > 0) then
                write(*,*) 'Error: molecular adress argument retrieval fails'
                stop
            else if (molstat < 0) then
                write(*,*) 'Error: molecular adress argument too long'
                stop
            end if
            call get_command_argument(2, address_basis, status=basstat)
            if (basstat > 0) then
                write(*,*) 'Error: basis set adress argument retrieval fails'
                stop
            else if (basstat < 0) then
                write(*,*) 'Error: basis set adress argument too long'
                stop
            end if
        end if
        call generate_tot(address_molecular)                                         ! generate .tot file before doing anything
        call read_gbs(address_basis)                                                 ! read basis set file before read geometry file
        call read_geometry(address_molecular)
        call input_print(address_basis, address_molecular)
        call read_keywords(address_molecular)
        if (DKH_order == 0) then
            call DKH_Hamiltonian()
            call DKH_Hartree_Fock()
        else if (DKH_order == 1) then
            write(60,*)
            write(60,'(A)') 'Module Hamiltonian: The SRTP effect requires at least a second-order'
            write(60,'(A)') '   DKH Hamiltonian to be accurately described, since the lead term of'
            write(60,'(A)') '   SRTP effect expanded at c=0 is in order c^-4. Only DKH2+ Hamiltonian'
            write(60,'(A)') '   guarantees the completeness of the expansion term at order c^-4.'
            call terminate('fpFW Hamiltonian is not supported by current version')
        else if (DKH_order == 2) then
            call DKH_Hamiltonian()
            call DKH_Hartree_Fock()
        else
            call terminate('DKH_order is set incorrectly')
        end if
        ! initialization
        deallocate(shell_in_element)
        deallocate(atom_basis)
        deallocate(molecular)
        deallocate(AO2MO)
        deallocate(rou_m)
        deallocate(i_j)
        deallocate(i_V_j)
        deallocate(i_p2_j)
        deallocate(S0_5)
        deallocate(exS0_5)
        deallocate(Fock1)
        deallocate(Fock2)
        deallocate(Fock)
        deallocate(exi_T_j)
        deallocate(exi_V_j)
        deallocate(oper)
        deallocate(oper3)
        deallocate(oper4)
        deallocate(oper6)
        deallocate(orbE)
        deallocate(isupp_ev_f)
        deallocate(Rsd)
        deallocate(DIISmat)
        deallocate(ipiv)
        deallocate(rou_history)
        deallocate(rou_pre)
        deallocate(Fock2_mic)
        deallocate(Fock2_assigned)
        deallocate(swintegral)
        deallocate(basis_inf)
        ndschwarz = .true.
        ini_rou = .true.
        if (DKH_order == 2) then
            deallocate(exSOC)
            deallocate(AO2p2)
            deallocate(evl_p2)
            if (SRTP_type) then
                deallocate(exSR)
            end if
        end if
        call terminate('normal')
    end subroutine elec_stc_calc
    
    
end module Indicator