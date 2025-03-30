!> @file Fundamentals.f90
!!
!! @brief basic processes and definitions
!!
!! @syntax Fortran 2008 free format
!!
!! @code UTF-8
!!
!! @author dirac4pi
module Fundamentals
  use ifport
  
!-----------------------------------------------------------------------
! definitions for programming
  
  integer :: loop_i, loop_j, loop_k                ! universal loop variables
  integer :: loop_m, loop_n, loop_l
  integer :: ios                                   ! universal operating status
  integer :: pid                                   ! process ID
  integer,parameter :: dp = selected_real_kind(12) ! double precesion
  integer,parameter :: sp = kind(0.0)              ! single precesion
  ! safe minimal that 1/safmin does not overflow
  real(dp),parameter :: safmin = 1E-14
  logical :: exists                                ! whether the file exists
  integer,private :: clock_time1, clock_time2      ! job clock time (wall time)
  real(sp),private :: cpu_time1, cpu_time2         ! job cpu time
  integer :: threads_use = 8                       ! Number of threads be used
  integer :: cpu_threads                           ! Number of threads in CPU
  character(len = 50) :: address_molecular
  character(len = 50) :: address_job
  character(len = 50) :: address_basis
  character(len = 50) :: wd                        ! working directory
  character(len = 20) :: usrname
  
!-----------------------------------------------------------------------
! definitions of physical and mathematical parameters

  complex(dp),parameter :: c0 = cmplx(0.0,0.0,dp)  ! 0 + 0i
  complex(dp),parameter :: c1 = cmplx(1.0,0.0,dp)  ! 1 + 0i
  complex(dp),parameter :: ci = cmplx(0.0,1.0,dp)  ! 0 + 1i
  real(dp),parameter :: pi = 3.14159265358979323_dp
  real(dp),parameter :: speedc = 137.035999074_dp  ! speed of light in vacuum
  real(dp),parameter :: Ang2Bohr = 0.529177249_dp
  real(dp),parameter :: fm2Bohr = 52917.7249_dp
  ! Correction factor for the QED radiation effect on the Bohr magnetic moment
  real(dp),parameter :: QED_rad = sqrt(1.0/(2.0*pi*speedc) - &
  0.328/(pi*pi*speedc*speedc))
  
!-----------------------------------------------------------------------
! calculation settings (with default)
  
  !-----------------<module Hamiltonian>-----------------
  integer :: DKH_order = 2               ! 0: nonrelativistic; 1: fpFW; 2: DKH2
  logical(kind=4) :: SRTP_type = .false. ! Second Relativized Thomas Precession
  logical(kind=4) :: STTP_type = .false. ! Spin Tensor Thomas Precession
  logical :: finitenuc = .false.
  real(dp) :: cutS = 1E-5                ! threshold of evl(i_j)
  !--------------------<module Atoms>--------------------
  integer :: charge = 0                  ! charge of the system
  integer :: spin_mult = 675             ! spin multiplicity of the system
  integer :: electron_count              ! number of electrons of the system
  logical :: s_h = .true.                ! Cartesian / spher-harmo basis
  !---------------------<module SCF>---------------------
  real(dp) :: schwarz_VT = 1E-9          ! Schwarz screening cutoff
  integer :: maxiter = 128               ! upper limit of convergence loops
  real(dp) :: conver_tol = 1E-6          ! convergence tolerence of energy
  real(dp) :: damp = 0.0                 ! dynamical damp (-(dE)^damp+1)
  integer :: nodiis = 8                  ! initial iteration steps without DIIS
  integer :: subsp = 5                   ! dimension of suboptimal subspace
  real(dp) :: diisdamp = 0.7             ! damp coefficient in DIIS
  real(dp) :: prtlev = 0.1               ! minimun AO coefficient print output
  real(dp) :: cutdiis = 0.0              ! cut DIIS when threshold is reached
  real(dp) :: cutdamp = 0.01             ! cut damp when threshold is reached
  logical :: keepspin = .false.          ! avoid spin mutations
  logical :: d4 = .false.                ! use DFT-D4 dispersion correction
  character(len=8) :: guess_type = 'gaussian' ! initial guess for SCF
  logical :: molden = .false.            ! save MOs to .molden file
  !------------------<module Functional>------------------
  ! https://libxc.gitlab.io/functionals/
  integer :: fx_id = -1                  ! exchange functional id, default HF
  integer :: fc_id = -1                  ! correlation functional id, default HF
  real(dp) :: x_HF = 1.0_dp              ! HF(exact exchange) componnet
  character(len=20) :: funcemd4 = 'hf'   ! functional name in emd4

  contains
  
!-----------------------------------------------------------------------
!> generation of output file
  subroutine generate_output()
    implicit none
    character(8)  :: date
    character(10) :: time
    character(5)  :: zone
    integer,dimension(8) :: values
    if (index(address_molecular,".xyz") /= 0) then
      address_job = address_molecular(1:index(address_molecular,".xyz") - 1)
      open(60, file = trim(address_job)//".esc", &
      status = "replace", action = "write")
    else
      address_job = address_molecular(1:index(address_molecular,&
      '\',back = .true.))
      open(60, file = trim(address_job)//"untitled.esc",&
      status = "replace", action = "write")
    end if
    call system_clock(clock_time1)
    call cpu_time(cpu_time1)
    call date_and_time(date,time,zone,values)
    write(60,"(A21,I4,A1,I2.2,A1,I2.2,A2,I2,A1,I2.2)") 'Job execution start: ',&
    values(1),'.',values(2),'.',values(3),'  ',values(5),':',values(6)
    ios = getcwd(wd)
    if (ios == 0) write(60,"(A)") "WD = "//wd
    pid = getpid()
    call getlog(usrname)
    write(60,"(A,I6,A,A)") "PID = ", pid, ";  username = ", usrname
    write(60,"(A60)") 'This file is the output file of job '//address_job
    write(60,"(A)") "Program/version: TRESC/development"
    write(60,*)
    write(60,"(A)") "Acknowledge: "
    write(60,"(A)") "  TRESC is a molecular 2-component DKH2 HF/KS SCF"
    write(60,"(A)") "  calculation program. For more:"
    write(60,"(A)") "  https://github.com/Dirac4pi/TRESC.git"
    write(60,"(A)") "  All results default to Atomic Units (A.U.)."
    write(60,*)
    write(60,"(A)") "External libs:"
    write(60,"(A)") "  IFPORT: cross-platform development"
    write(60,"(A)") "  Intel MKL: linear algebra routines"
    write(60,"(A)") "  OpenMP: multithreaded parallel operations"
    write(60,"(A)") "  Libxc: exchange, correlation functionals"
    write(60,*)
  end subroutine generate_output
  
!-----------------------------------------------------------------------
!> standard termination of current job
  subroutine terminate(terminate_message)
    implicit none
    logical :: isopen
    character(len = *),intent(in) :: terminate_message
    call lowercase(terminate_message)
    inquire(unit = 60, opened = isopen)
    if (isopen) then
      if (terminate_message == 'normal' .or. terminate_message == 'keep') then
        write(60,*)
        write(60,"(A)") '----------<NORMAL TERMINATION>----------'
      else if (terminate_message /= '') then
        write(60,*)
        write(60,"(A)") '----------<ERROR TERMINATION>----------'
        write(60,"(A)") "Error message: "//terminate_message
      else
        write(60,*)
        write(60,"(A)") '----------<ERROR TERMINATION>----------'
        write(60,"(A)") "Unknown error"
      end if
      call system_clock(clock_time2)
      call cpu_time(cpu_time2)
      write(60,"(A,I2,A,I2,A,I2,A)") 'Job  cpu  time: ',&
      floor((cpu_time2 - cpu_time1)/3600_sp),' h' &
      ,mod(floor((cpu_time2 - cpu_time1)/60_sp),60),' min',&
      mod(int(cpu_time2 - cpu_time1),60),' s'
      write(60,"(A,I2,A,I2,A,I2,A)") 'Job clock time: ',&
      floor(real(clock_time2 - clock_time1)/36000000_sp),' h' &
      ,mod(floor(real(clock_time2 - clock_time1)/600000_sp),60),' min',&
      mod((clock_time2 - clock_time1)/10000,60),' s'
      if (terminate_message == 'keep') write(60,'(A)') 'process holding...'
      close(60)
      if (terminate_message == 'normal') then
        write(*,"(A)") 'TRESC: Normal termination of job '//address_job
        stop 0
      else if (terminate_message == 'keep') then
        write(*,"(A)") 'TRESC: Normal termination of job '//address_job
        write(*,'(A)') 'holding...'
      else
        write(0,"(A)") 'TRESC Error: '//terminate_message
        stop 1
      end if
    else
      if (terminate_message == 'normal') then
        write(*,"(A)") 'TRESC: Normal termination'
        stop 0
      else if (terminate_message == 'keep') then
        write(*,"(A)") 'TRESC: holding...'
      else
        write(0,"(A)") 'TRESC Error: '//terminate_message
        stop 1
      end if
    end if
  end subroutine terminate
  
!-----------------------------------------------------------------------
!> dump complex matrix (double precesion) to address_job.name binary file
  subroutine dump_matrix_cmplx(name, m, dmi, dmj)
    character(len = *),intent(in) :: name
    integer,intent(in) :: dmi, dmj
    integer :: channel, a, b
    complex(dp),intent(in) :: m(dmi, dmj)
    if (size(m) < dmi*dmj) then
      call terminate('dump matrix failed, dm too large')
    else if (size(m) > dmi*dmj) then
      call terminate('dump matrix failed, dm too small')
    end if
    call lowercase(name)
    open(newunit=channel, file=trim(address_job)//'.'//trim(name), &
    access='direct', form='unformatted', &
    recl=16, status='replace', action='write', iostat=ios)
    if (ios /= 0) call terminate('dump binary file .'//trim(name)//' failed')
    write(channel,rec=1) dmi
    write(channel,rec=2) dmj
    do a = 1, dmi
      do b = 1, dmj
        write(channel,rec=2+dmj*(a-1)+b)  m(a,b)
      end do
    end do
    close(channel)
  end subroutine dump_matrix_cmplx
  
!-----------------------------------------------------------------------
!> dump real matrix (double precesion) to address_job.name binary file
  subroutine dump_matrix_real(name, m, dmi, dmj)
    character(len = *),intent(in) :: name
    integer,intent(in) :: dmi, dmj
    integer :: channel, a, b
    real(dp),intent(in) :: m(dmi, dmj)
    if (size(m) < dmi*dmj) then
      call terminate('dump matrix failed, dm too large')
    else if (size(m) > dmi*dmj) then
      call terminate('dump matrix failed, dm too small')
    end if
    call lowercase(name)
    open(newunit=channel, file=trim(address_job)//'.'//trim(name), &
    access='direct', form='unformatted', &
    recl=8, status='replace', action='write', iostat=ios)
    if (ios /= 0) call terminate('dump binary file .'//trim(name)//' failed')
    write(channel,rec=1) dmi
    write(channel,rec=2) dmj
    do a = 1, dmi
      do b = 1, dmj
        write(channel,rec=2+dmj*(a-1)+b)  m(a,b)
      end do
    end do
    close(channel)
  end subroutine dump_matrix_real
!-----------------------------------------------------------------------
!> load complex matrix (double precesion) from address_job.name binary file
  subroutine load_matrix_cmplx(name, m, dmi, dmj)
    character(len = *),intent(in) :: name
    integer :: channel, a, b
    integer,intent(out) :: dmi, dmj
    complex(dp),allocatable :: m(:,:)
    if (allocated(m)) then
      dmi = -1
      dmj = -1
      return
    end if
    call lowercase(name)
    open(newunit=channel, file=trim(address_job)//'.'//trim(name), &
    access='direct', form='unformatted',&
    recl=16, status='old', action='read', iostat=ios)
    if (ios /= 0) call terminate('read binary file .'//trim(name)//' failed')
    read(channel,rec=1) dmi
    read(channel,rec=2) dmj
    allocate(m(dmi,dmj))
    m = c0
    do a = 1, dmi
      do b = 1, dmj
        read(channel,rec=2+dmj*(a-1)+b)  m(a,b)
      end do
    end do
    close(channel)
  end subroutine load_matrix_cmplx
  
!-----------------------------------------------------------------------
!> load real matrix (double precesion) from address_job.name binary file
  subroutine load_matrix_real(name, m, dmi, dmj)
    character(len = *),intent(in) :: name
    integer :: channel, a, b
    integer,intent(out) :: dmi, dmj
    real(dp),allocatable :: m(:,:)
    if (allocated(m)) then
      dmi = -1
      dmj = -1
      return
    end if
    call lowercase(name)
    open(newunit=channel, file=trim(address_job)//'.'//trim(name), &
    access='direct', form='unformatted',&
    recl=8, status='old', action='read', iostat=ios)
    if (ios /= 0) call terminate('read binary file .'//trim(name)//' failed')
    read(channel,rec=1) dmi
    read(channel,rec=2) dmj
    allocate(m(dmi,dmj))
    m = 0.0_dp
    do a = 1, dmi
      do b = 1, dmj
        read(channel,rec=2+dmj*(a-1)+b)  m(a,b)
      end do
    end do
    close(channel)
  end subroutine load_matrix_real
  
!-----------------------------------------------------------------------
!> convert a string to lowercase
  subroutine lowercase(inpstr)
    character(len = *) :: inpstr
    integer :: lenstr, a, asc
    lenstr = len(inpstr)
    do a = 1, lenstr, 1
      asc = iachar(inpstr(a:a))
      if (asc >= 65 .and. asc <= 90) inpstr(a:a) = achar(asc + 32)
    end do
  end subroutine lowercase
  
end module Fundamentals