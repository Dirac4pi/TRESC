! |-------<THOMAS RELATIVISTIC ELECTRONIC STRUCTURE CALCULATIONS>------|
! |-----------------------------<DIRAC4PI>-----------------------------|
! |---------------------------<FUNDAMENTALS>---------------------------|
! |---------------<BASIC CONTROL & PARAMETER DEFINITION>---------------|

module Fundamentals
!-----------------------------------------------------------------------
! definitions for programming
    
    integer :: loop_i, loop_j, loop_k, loop_m, loop_n, loop_l     ! universal loop variables
    integer :: ios                                                ! universal operating status
    integer,parameter :: dp = SELECTED_REAL_KIND(12)
    integer,parameter :: sp = kind(0.0)
    integer :: exists
    character(len = 50) :: job_name
    integer,private :: clock_time1, clock_time2                   ! job clock time
    real(sp),private :: cpu_time1, cpu_time2                      ! job cpu time
    integer :: threads_use = 8                                    ! Number of threads be used
    integer :: cpu_threads                                        ! Number of threads in CPU
    integer :: stacksize = 2000                                   ! stack size of each thread (MB)
    character(len = 6) :: stackchar
    
    
!-----------------------------------------------------------------------
! definitions of physical and mathematical parameters

    real(dp),parameter :: pi = 3.14159265358979323_dp
    real(dp),parameter :: speed_light = 137.035999074_dp
    real(dp),parameter :: Ang2Bohr = 0.529177249_dp
    real(dp),parameter :: fm2Bohr = 52917.7249_dp
    real(dp),parameter :: QED_rad = sqrt(1.0/(2.0*pi*speed_light) - 0.328/(pi*pi*speed_light*speed_light))
    
!-----------------------------------------------------------------------
! calculation settings (with default)
    
    !-----------------<module Hamiltonian>-----------------
    integer :: DKH_order = 2                                      ! 0: nonrelativistic; 1: fpFW; 2: DKH2
    logical(kind=4) :: Breit_type = .false.                       ! Breit term (electron-electron magnetic interaction and retarded effect)
    logical(kind=4) :: SRTP_type = .false.                        ! Second Relativized Thomas Precession
    logical(kind=4) :: STTP_type = .false.                        ! Spin Tensor Thomas Precession
    logical(kind=4) :: mDCB_type = .false.                        ! modified Dirac Coulomb Breit
    logical :: finitenuc = .false.
    !--------------------<module Atoms>--------------------
    integer :: charge = 0                                         ! charge of the system
    integer :: spin_mult = 675                                    ! spin multiplicity of the system, its default value depends on the number of electrons
    integer :: electron_count                                     ! number of electrons in molecular
    !---------------------<module SCF>---------------------
    real(dp) :: schwarz_VT = 1E-10                                ! Schwarz screening cutoff, default as 1E-10
    integer :: maxiter = 128                                      ! upper limit of convergence loops
    real(dp) :: conver_tol = 1E-8                                 ! convergence tolerence of energy
    integer :: nodiis = 8                                         ! initial iteration steps without any mixing
	integer :: subsp = 5                                          ! dimension of suboptimal subspace
    real(dp) :: damp = 0.7                                        ! mix the generated coefficient and original coefiicient by residual
    real(dp) :: prtlev = 0.1                                      ! print level (minimun coefficient of AO will be printed in output file)
    real(dp) :: cutdiis = 1E-5                                    ! cut DIIS when threshold is reached
    logical :: keepspin = .false.                                 ! whether to avoid spin multiplicity mutations in case of degenerate frontier orbitals
    !------------------<module Radiation>------------------
    
    contains
    
!-----------------------------------------------------------------------
! standard generation of .tot file
    
    subroutine generate_tot(address_molecular)
        implicit none
        character(len = *),intent(in) :: address_molecular
        character(len = len(address_molecular)) :: address_output
        character(8)  :: date
        character(10) :: time
        character(5)  :: zone
        integer,dimension(8) :: values
        if (index(address_molecular,".tip") /= 0) then
            address_output = address_molecular(1:index(address_molecular,".tip") - 1)
            job_name = address_molecular(1 : index(address_molecular,".tip") - 1 )
            open(60, file = trim(address_output)//".tot", status = "replace", action = "write")
            open(64, file = trim(address_output)//".ao2mo",Access = 'direct', form = 'unformatted', RecL = 16, status = "replace", action = "write")
        else
            address_output = address_molecular(1:index(address_molecular,'\',back = .true.))
            job_name = 'untitled'
            open(60, file = trim(address_output)//"untitled.tot", status = "replace", action = "write")
        end if
        call system_clock(clock_time1)
        call cpu_time(cpu_time1)
        call date_and_time(date,time,zone,values)
        write(60,"(A21 I4 '.' I2.2 '.' I2.2 '  ' I2 ':' I2.2)") 'Job execution start: ', values(1), values(2), values(3), values(5), values(6)
        write(60,"('This file is the output file of job ' A20)") job_name
        write(60,"(a)") "Program/version: TRESC/development"
        write(60,*)
        write(60,"(a)") "Acknowledge: "
        write(60,"(a)") "   TRESC is a 2-component DKH2 Hartree-Fock molecular self-consistent field"
        write(60,"(a)") "   calculation program, use Cartesian basis functions for all basis."
        write(60,"(a)") "   Make sure you have installed Gaussian and set the environment variables."
        write(60,"(a)") "   During the calculation, Converged orbital coefficients from Gaussian will"
        write(60,"(a)") "   be used to assign the initial density matrix, see 'Gaujob.gjf'."
        write(60,*)
        write(60,"(a)") "Packages used: "
        write(60,"(a)") "   LAPACK: linear algebra routines"
        write(60,"(a)") "   BLAS: basic vector and matrix operations"
        write(60,*)
    end subroutine generate_tot
    
!-----------------------------------------------------------------------
! standard termination of current job
    
    subroutine terminate(terminate_message)
        implicit none
        character(len = *),intent(in) :: terminate_message
        if (terminate_message == 'normal') then
            write(60,*)
            write(60,"(a)") '---------<NORMAL TERMINATION>---------'
        else if (terminate_message /= '') then
            write(60,*)
            write(60,"(a)") '---------<ERROR TERMINATION>---------'
            write(60,"(a)") "Error message: "//terminate_message
        else
            write(60,*)
            write(60,"(a)") '---------<ERROR TERMINATION>---------'
            write(60,"(a)") "Cause of error unknown."
        end if
        call system_clock(clock_time2)
        call cpu_time(cpu_time2)
        write(60,"(A,I2,A,I2,A,I2,A)") 'Job  cpu  time: ',floor((cpu_time2 - cpu_time1)/3600_sp),' h' &
        ,mod(floor((cpu_time2 - cpu_time1)/60_sp),60),' min',mod(int(cpu_time2 - cpu_time1),60),' s'
        write(60,"(A,I2,A,I2,A,I2,A)") 'Job clock time: ',floor(real(clock_time2 - clock_time1)/36000000_sp),' h' &
        ,mod(floor(real(clock_time2 - clock_time1)/600000_sp),60),' min',mod((clock_time2 - clock_time1)/10000,60),' s'
        close(60)
        if (terminate_message == 'normal') then
            write(*,"('TRESC: Normal termination of job ' A20)") job_name
            ! normal termination won't stop the program
        else
            write(*,"('TRESC: Error termination of job ' A20)") job_name
            pause
            stop
        end if
    end subroutine terminate
    
!-----------------------------------------------------------------------
! bionomial coefficient of given N, M (C_N^M)
    
    real(dp) function binomialcoe(N,M)
        integer,intent(in) :: N,M
        integer :: bloop_i, numerator, denominator
        if (N < M .or. N < 0 .or. M < 0) then
            call terminate('binomialcoe is called incorrectly')
        else if(M == 0 .or. N == 0 .or. M == N) then
            binomialcoe = 1.0
            return
        end if
        binomialcoe = 1.0
        numerator = 1
        bloop_i = N
        do while(bloop_i > M)
            numerator = numerator * bloop_i
            bloop_i = bloop_i - 1
        end do
        denominator = 1
        bloop_i = N - M
        do while(bloop_i > 1)
            denominator = denominator * bloop_i
            bloop_i = bloop_i - 1
        end do
        binomialcoe = real(numerator) / real(denominator)
        return
    end function binomialcoe
    
!-----------------------------------------------------------------------
! normalization factor of each basis
    
    real(dp) function AON(a,l,m,n)
        real(dp),intent(in) :: a
        integer,intent(in) :: l,m,n
        integer :: s,hl,hm,hn
        if (l == 0) then
            hl = 1
        else
            hl = 1
            do s=1, 2*l-1, 2
                hl = hl * s
            end do
        end if
        if (m == 0) then
            hm = 1
        else
            hm = 1
            do s=1, 2*m-1, 2
                hm = hm * s
            end do
        end if
        if (n == 0) then
            hn = 1
        else
            hn = 1
            do s=1, 2*n-1, 2
                hn = hn * s
            end do
        end if
        AON = (2.0_dp*a/pi)**(0.75)*sqrt((4.0_dp*a)**(l+m+n)/(real(hl)*real(hm)*real(hn)))
        return
    end function AON
    
!-----------------------------------------------------------------------
! factorial of given a
    
    real(dp) function factorial(a)
        integer,intent(in) :: a
        integer :: na
        factorial = 1.0_dp
        if (a == 0) then
            return
        else if (a < 0) then
            call terminate('factorial is called incorrectly')
        else
            do na=1, a
                factorial = factorial * real(na)
            end do
        end if
        return
    end function factorial
    
end module Fundamentals