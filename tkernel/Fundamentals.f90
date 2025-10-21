!> @file Fundamentals.f90
!!
!! @brief definitions, fundamental processes and matrix operations
!!
!! @syntax Fortran 2008 free format
!!
!! @code UTF-8
!!
!! @author dirac4pi

module Fundamentals
  use IFPORT
  use LAPACK95
  
!-----------------------------------------------------------------------
! definitions for programming
  
  integer             :: ios                        ! universal operating status
  integer             :: pid                        ! process ID
  integer,parameter   :: dp = selected_real_kind(12)! double precesion
  integer,parameter   :: sp = kind(0.0)             ! single precesion
  ! safe minimal that 1/safmin does not overflow
  real(dp),parameter  :: safmin = 1E-14
  integer             :: info                       ! info from lapack routines
  logical             :: exists                     ! whether the file exists
  integer,private     :: clock_time1, clock_time2   ! job clock time (wall time)
  real(sp),private    :: cpu_time1, cpu_time2       ! job cpu time
  integer             :: threads = 8                ! number of threads
  integer             :: new_threads
  integer,parameter   :: align_size = 32            ! AVX2:32, AVX512:64
  !!DIR$ ATTRIBUTES ALIGN:align_size :: mat
  integer             :: nproc                      ! number of processors
  character(len=200)  :: address_molecule
  character(len=200)  :: address_job
  character(len=100)  :: jobname
  character(len=200)  :: address_basis
  character(len=200)  :: address_molden
  character(len=200)  :: address_fch
  character(len=200)  :: wd = ''                    ! working directory
  character(len=50)   :: usrname = ''
  character(len=10),parameter   :: version = 'dev'  ! TRESC version
  
!-----------------------------------------------------------------------
! definitions of physical and mathematical parameters

  complex(dp),parameter :: c0       = cmplx(0.0,0.0,dp)     ! 0 + 0i
  complex(dp),parameter :: c1       = cmplx(1.0,0.0,dp)     ! 1 + 0i
  complex(dp),parameter :: ci       = cmplx(0.0,1.0,dp)     ! 0 + 1i
  real(dp),parameter    :: pi       = 3.14159265358979324_dp
  real(dp),parameter    :: rpi      = 1.7724538509055159_dp ! dsqrt(pi)
  real(dp),parameter    :: c        = 137.035999074_dp      ! vacuum light speed
  real(dp),parameter    :: c2       = c**2                  ! c^2
  real(dp),parameter    :: Ang2Bohr = 0.529177249_dp
  real(dp),parameter    :: fm2Bohr  = 52917.7249_dp
  ! Correction factor for the QED radiation effect on the Bohr magnetic moment
  real(dp),parameter    :: QED_rad  = sqrt(1.0/(2.0*pi*c) - 0.328/(pi*pi*c2))
  
!-----------------------------------------------------------------------
! calculation settings (with default)
  
  !-----------------<module Hamiltonian>-----------------
  logical(kind=4) :: pVp1e       = .false.! one-electron pVp potetial (spinor)
  logical(kind=4) :: pVp2e       = .false.! two-electron pVp potetial (spinor)
  logical(kind=4) :: pppVp       = .false.! Second Relativized Thomas Precession
  real(dp)        :: cutS        = 1E-5   ! threshold of evl(i_j)
  !--------------------<module Atoms>--------------------
  integer         :: charge      = 0      ! charge of the system
  integer         :: spin_mult   = 675    ! spin multiplicity of the system
  integer         :: electron_count       ! number of electrons of the system
  !---------------------<module SCF>---------------------
  real(dp)        :: schwarz_VT  = 1E-9   ! Schwarz screening cutoff
  integer         :: maxiter     = 128    ! upper limit of convergence loops
  real(dp)        :: conver_tol  = 1E-6   ! convergence tolerence of energy
  real(dp)        :: damp        = 0.0    ! dynamical damp (-(dE)^damp+1)
  real(dp)        :: damp_new
  integer         :: nodiis      = 8      ! initial iteration steps without DIIS
  integer         :: subsp       = 5      ! dimension of suboptimal subspace
  real(dp)        :: diisdamp    = 0.7    ! damp coefficient in DIIS
  real(dp)        :: diisdamp_new
  real(dp)        :: prtlev      = 0.1    ! minimun AO coefficient print output
  real(dp)        :: cutdiis     = 0.0    ! cut DIIS when threshold is reached
  real(dp)        :: cutdiis_new
  real(dp)        :: cutdamp     = 0.01   ! cut damp when threshold is reached
  real(dp)        :: cutdamp_new
  character       :: cspin       = 'n'    ! constrained spin multiplicity mode
                                          ! n: normal  d: degenerate  f: force
  logical         :: d4          = .false.! use DFT-D4 dispersion correction
  character(len=8):: guess_type  = 'molden'  ! initial guess for SCF
  logical         :: molden      = .false.! save MOs to .molden file
  !------------------<module Representation>------------------
  real(dp)        :: beta(3)     = 0.0_dp ! v/c of frame of motion
  real(dp)        :: beta2       = 0.0_dp ! |beta|**2
  real(dp)        :: gamma       = 1.0_dp ! Lorentz contraction
  !------------------<module Functional>------------------
  ! https://libxc.gitlab.io/functionals/
  integer         :: fx_id       = -1     ! exchange functional ID(default HF
  integer         :: fc_id       = -1     ! correlation functional ID(default HF
  real(dp)        :: x_HF        = 1.0_dp ! HF(exact exchange) componnet
  character(len=20) :: funcemd4  = 'hf'   ! functional name in emd4

  interface dump_matrix
    module procedure dump_matrix_cmplx
    module procedure dump_matrix_real
  end interface dump_matrix

  interface load_matrix
    module procedure load_matrix_cmplx
    module procedure load_matrix_real
  end interface load_matrix

  interface inverse
    module procedure real_symm_inverse
    module procedure cmplx_inverse
  end interface inverse

  interface diag
    module procedure real_symm_diag
    module procedure cmplx_diag
  end interface diag

  interface matmul
    module procedure real_matmul
    module procedure cmplx_matmul
    module procedure real_matvec
    module procedure cmplx_matvec
  end interface matmul

  interface norm
    module procedure real_fbnorm
    module procedure cmplx_fbnorm
  end interface norm

  interface det
    module procedure real_det
    module procedure cmplx_det
  end interface det

  private :: dump_matrix_real, dump_matrix_cmplx, load_matrix_real
  private :: load_matrix_cmplx, real_symm_inverse, cmplx_inverse
  private :: real_symm_diag, cmplx_diag, real_matmul, cmplx_matmul
  private :: real_matvec, cmplx_matvec, real_fbnorm, cmplx_fbnorm
  private :: real_det, cmplx_det
  public  :: dump_matrix, load_matrix, generate_output, terminate, matmul, norm
  public  :: lowercase, is_alpha, diag, inverse, can_orth, symm_orth, atnz2block
  public  :: wigner_d, binomialcoe, factorial, det

  contains
  
!-----------------------------------------------------------------------
!> generation of output file
  subroutine generate_output()
    implicit none
    character(8)         :: date
    character(10)        :: time
    character(5)         :: zone
    integer              :: values(8)
    if (index(address_molecule,".xyz") /= 0) then
      address_job = address_molecule(1:index(address_molecule,".xyz") - 1)
      jobname = address_job(index(address_molecule,'/',back=.true.)+1 :)
      open(60, file=trim(address_job)//".esc", status="replace", action="write")
    else
      address_job = address_molecule(1:index(address_molecule,'/',back=.true.))
      jobname = 'untitled'
      open(60, file=trim(address_job)//"untitled.esc", &
      status="replace", action="write")
    end if
    call system_clock(clock_time1)
    call cpu_time(cpu_time1)
    call date_and_time(date,time,zone,values)
    write(60,"(A21,I4,A1,I2.2,A1,I2.2,A2,I2,A1,I2.2)") 'Job execution start: ',&
    values(1),'.',values(2),'.',values(3),'  ',values(5),':',values(6)
    ios = getcwd(wd)
    if (ios == 0) write(60,"(A)") "WD = "//trim(wd)
    pid = getpid()
    call getlog(usrname)
    write(60,"(A,I7,A,A)") "PID = ", pid, ";  username = ", trim(usrname)
    write(60,"(A)") 'This file is the output file of job '//trim(address_job)
    write(60,"(A)") "Program/version: TRESC/"//trim(version)
    write(60,*)
    write(60,"(A)") "Acknowledge: "
    write(60,"(A)") "  TRESC is a molecular 2-component DKH2 HF/KS-SCF"
    write(60,"(A)") "  calculation program. For more:"
    write(60,"(A)") "  https://github.com/Dirac4pi/TRESC"
    write(60,"(A)") "  All results default to Atomic Units (A.U.)."
    write(60,*)
    write(60,"(A)") "External libs:"
    write(60,"(A)") "  IFPORT: cross-platform Fortran development"
    write(60,"(A)") "  Intel MKL: linear algebra routines"
    write(60,"(A)") "  OpenMP: multithreaded parallel operations"
    write(60,"(A)") "  Libxc: exchange & correlation functionals"
    write(60,*)
  end subroutine generate_output
  
!-----------------------------------------------------------------------
!> standard termination of current job
  subroutine terminate(terminate_message)
    implicit none
    logical                       :: isopen
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
        write(*,"(A)") 'tkernel: Normal termination of job '//trim(address_job)
        stop 0
      else if (terminate_message == 'keep') then
        write(*,"(A)") 'tkernel: Normal termination of job '//trim(address_job)
        write(*,'(A)') 'holding...'
      else
        write(0,"(A)") 'tkernel error: '//terminate_message
        stop 1
      end if
    else
      if (terminate_message == 'normal') then
        write(*,"(A)") 'tkernel: Normal termination'
        stop 0
      else if (terminate_message == 'keep') then
        write(*,"(A)") 'tkernel: holding...'
      else
        write(0,"(A)") 'tkernel error: '//terminate_message
        stop 1
      end if
    end if
  end subroutine terminate
  
!-----------------------------------------------------------------------
!> dump complex matrix (double precesion) to address_job.name binary file
  subroutine dump_matrix_cmplx(name, m, dmi, dmj)
    character(len = *),intent(in) :: name
    integer,intent(in)            :: dmi, dmj
    integer                       :: channel, a, b
    complex(dp),intent(in)        :: m(dmi, dmj)
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
    integer,intent(in)            :: dmi, dmj
    integer                       :: channel, a, b
    real(dp),intent(in)           :: m(dmi, dmj)
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
    integer                       :: channel, a, b
    integer,intent(out)           :: dmi, dmj
    complex(dp),allocatable       :: m(:,:)
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
    integer                       :: channel, a, b
    integer,intent(out)           :: dmi, dmj
    real(dp),allocatable          :: m(:,:)
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
!> determine if a character is an alphabet
  logical function is_alpha(ch)
    implicit none
    character, intent(in) :: ch
    integer               :: ascii
    ascii = iachar(ch)
    is_alpha = (ascii>=65 .and. ascii<=90) .or. (ascii>=97 .and. ascii<=122)
  end function is_alpha

!-----------------------------------------------------------------------
!> convert a string to lowercase
  subroutine lowercase(inpstr)
    implicit none
    character(len = *) :: inpstr
    integer            :: lenstr, a, asc
    lenstr = len(inpstr)
    do a = 1, lenstr, 1
      asc = iachar(inpstr(a:a))
      if (asc >= 65 .and. asc <= 90) inpstr(a:a) = achar(asc + 32)
    end do
  end subroutine lowercase

!-----------------------------------------------------------------------
!> diagonalisation of given real symmetric matrix
  subroutine real_symm_diag(mat, dm, U, evl)
    implicit none
    integer,intent(in)   :: dm
    real(dp),intent(in)  :: mat(dm,dm)
    real(dp),intent(out) :: U(dm,dm)
    real(dp),intent(out) :: evl(dm)
    integer              :: evl_count   ! number of eigenvalues found by dsyevr
    integer              :: isupev(2*dm)! indices of the nonzero elements
    integer              :: lwork       ! lwork for input of dsyevr
    integer              :: liwork      ! liwork for input of dsyevr
    integer,allocatable  :: iwork(:)    ! iwork for input of dsyevr
    real(dp)             :: mat_u(dm,dm)! upper triangular part
    real(dp),allocatable :: work(:)     ! work for input of dsyevr
    integer              :: fi, fj      ! loop variables for real_symm_diag
    allocate(work(1))
    allocate(iwork(1))
    mat_u = 0.0_dp
    do fi = 1, dm
      do fj = fi, dm
        mat_u(fi,fj) = mat(fi,fj)
      end do
    end do
    call dsyevr('V','A','U',dm,mat_u,dm,0.0_dp,0.0_dp,0,0,safmin,&
    evl_count,evl,U,dm,isupev,work,-1,iwork,-1,info)
    liwork = iwork(1)
    lwork = nint(work(1))
    deallocate(work)
    deallocate(iwork)
    allocate(work(lwork))
    allocate(iwork(liwork))
    call dsyevr('V','A','U',dm,mat_u,dm,0.0_dp,0.0_dp,0,0,safmin,&
    evl_count,evl,U,dm,isupev,work,lwork,iwork,liwork,info)
    if(info < 0) then
      call terminate('real_symm_diag: illegal input of dsyevr')
    else if(info > 0) then
      call terminate('real_symm_diag: internal error of dsyevr')
    end if
    if (evl_count /= dm) &
    call terminate('real_symm_diag: diagonalisation failed, evl_count /= dm')
    deallocate(work)
    deallocate(iwork)
  end subroutine real_symm_diag

!-----------------------------------------------------------------------
!> diagonalisation of given complex Hermitian matrix
  subroutine cmplx_diag(mat, dm, U, evl)
    implicit none
    integer,intent(in)      :: dm
    complex(dp),intent(in)  :: mat(dm,dm)
    complex(dp),intent(out) :: U(dm,dm)
    real(dp),intent(out)    :: evl(dm)
    integer                 :: evl_count   ! number of eigenvalues
    integer                 :: isupev(2*dm)! indices of the nonzero elements
    integer                 :: lwork       ! lwork for input of zheevr
    integer                 :: lrwork      ! lrwork for input of zheevr
    integer                 :: liwork      ! liwork for input of zheevr
    complex(dp),allocatable :: work(:)     ! work for input of zheevr
    real(dp),allocatable    :: rwork(:)    ! rwork for input of zheevr
    integer,allocatable     :: iwork(:)    ! iwork for input of zheevr
    allocate(work(1))
    allocate(rwork(1))
    allocate(iwork(1))
    call zheevr('V','A','U',dm,mat,dm,0.0_dp,0.0_dp,0,0,safmin,&
    evl_count,evl,U,dm,isupev,work,-1,rwork,-1,iwork,-1,info)
    lwork = nint(real(work(1)))
    lrwork = nint(rwork(1))
    liwork = iwork(1)
    deallocate(work)
    deallocate(rwork)
    deallocate(iwork)
    allocate(work(lwork))
    allocate(rwork(lrwork))
    allocate(iwork(liwork))
    call zheevr('V','A','U',dm,mat,dm,0.0_dp,0.0_dp,0,0,safmin,&
    evl_count,evl,U,dm,isupev,work,lwork,rwork,lrwork,iwork,liwork,info)
    if(info < 0) then
      call terminate('cmplx_diag: illegal input of zheevr')
    else if(info > 0) then
      call terminate('cmplx_diag: internal error of zheevr')
    end if
    if (evl_count /= dm) &
    call terminate('cmplx_diag: diagonalisation failed, evl_count /= dm')
    deallocate(work)
    deallocate(rwork)
    deallocate(iwork)
  end subroutine cmplx_diag

!-----------------------------------------------------------------------
!> inverse of given real symmetric matrix
  subroutine real_symm_inverse(mat, dm)
    implicit none
    integer,intent(in)      :: dm
    real(dp)                :: mat(dm,dm)
    real(dp)                :: mat_u(dm,dm)
    real(dp),allocatable    :: work(:)
    integer                 :: lwork
    integer                 :: ipiv(dm)
    integer                 :: fi, fj
    mat_u = 0.0_dp
    do fi = 1, dm
      do fj = fi, dm
        mat_u(fi,fj) = mat(fi,fj)
      end do
    end do
    allocate(work(1))
    call dsytrf('U', dm, mat_u, dm, ipiv, work, -1, info)
    lwork = nint(work(1))
    deallocate(work)
    allocate(work(lwork))
    call dsytrf('U', dm, mat_u, dm, ipiv, work, lwork, info)
    if (info < 0) then
      call terminate('real_symm_inverse: illegal input of dsytrf')
    else if (info > 0) then
      call terminate('real_symm_inverse: dsytrf error: matrix is singular')
    end if
    call dsytri('U', dm, mat_u, dm, ipiv, work, info)
    deallocate(work)
    if (info < 0) then
      call terminate('real_symm_inverse: illegal input of dsytri')
    else if (info > 0) then
      call terminate('real_symm_inverse: dsytri error: internel error')
    end if
    mat = mat_u
    do fi = 1, dm
      do fj = fi, dm
        mat(fj,fi) = mat_u(fi,fj)
      end do
    end do
  end subroutine real_symm_inverse
  
!-----------------------------------------------------------------------
!> inverse of given complex matrix
  subroutine cmplx_inverse(mat, dm)
    implicit none
    integer,intent(in)      :: dm
    complex(dp)             :: mat(dm,dm)
    complex(dp),allocatable :: work(:)
    integer                 :: lwork
    integer                 :: ipiv(dm)
    call zgetrf(dm, dm, mat, dm, ipiv, info)
    if (info < 0) then
      call terminate('cmplx_inverse: illegal input of zgetrf')
    else if (info > 0) then
      call terminate('cmplx_inverse: zgetrf error: matrix is singular')
    end if
    allocate(work(1))
    call zgetri(dm, mat, dm, ipiv, work, -1, info)
    lwork = nint(real(work(1)))
    deallocate(work)
    allocate(work(lwork))
    call zgetri(dm, mat, dm, ipiv, work, lwork, info)
    deallocate(work)
    if (info < 0) then
      call terminate('cmplx_inverse: illegal input of zgetri')
    else if (info > 0) then
      call terminate('cmplx_inverse: zgetri error: internel error')
    end if
  end subroutine cmplx_inverse

!-----------------------------------------------------------------------
!> symmetric orthogonalisation of given real symmetric matrix
!!
!! X mat XT = I
  subroutine symm_orth(mat, dm, X, min_evl)
    implicit none
    integer,intent(in)   :: dm
    real(dp),intent(in)  :: mat(dm,dm)
    real(dp),intent(out) :: X(dm,dm)
    real(dp),intent(out) :: min_evl
    real(dp)             :: U(dm,dm)
    real(dp)             :: supU(dm,dm)
    real(dp)             :: evl(dm)
    real(dp)             :: mat_u(dm,dm)
    integer              :: fi
    call real_symm_diag(mat, dm, U, evl)
    min_evl = evl(1)
    mat_u = 0.0_dp
    do fi = 1, dm
      if (evl(fi) < 0.0) call terminate(&
      'symm_orth: eigenvalue less than zero, may due to code error')
      if (abs(evl(fi)) < safmin) then
        write(60,'(A)') &
        '  --matrix is not invertible (not full rank), call can_orth directly'
        min_evl = 0.0_dp
        return
      end if
      if (min_evl > evl(fi)) min_evl = evl(fi)
      mat_u(fi,fi) = 1.0_dp / dsqrt(evl(fi))
    end do
    call real_matmul('N', 'N', U, mat_u, supU)
    call real_matmul('N', 'T', supU, U, X)
  end subroutine symm_orth
  
!-----------------------------------------------------------------------
!> canonical orthogonalisation of given real symmetric matrix
!!
!! X mat XT = I
  subroutine can_orth(mat, dm, X, dm2)
    implicit none
    integer,intent(in)               :: dm
    real(dp),intent(in)              :: mat(dm,dm)
    real(dp),allocatable,intent(out) :: X(:,:)
    integer,intent(out)              :: dm2
    real(dp)                         :: U(dm,dm)
    real(dp)                         :: evl(dm)
    integer                          :: supU(dm)
    integer                          :: fi, fj, fk
    call real_symm_diag(mat, dm, U, evl)
    supU = 0
    dm2 = dm
    do fi = 1, dm
      if (evl(fi) < cutS) then
        supU(fi) = 1
        dm2 = dm2 - 1
      end if
    end do
    allocate(X(dm,dm2))
    fk = 1
    do fi = 1, dm
      if (supU(fi) == 1) cycle
      do fj = 1, dm
        X(fj,fk) = U(fj,fi) / dsqrt(evl(fi))
      end do
      fk = fk + 1
    end do
  end subroutine can_orth

!------------------------------------------------------------
!> transfer alternating zero matrix to block matrix
  subroutine atnz2block(atnz, dm)
    implicit none
    integer, intent(in) :: dm
    integer             :: fi, fj
    complex(dp)         :: atnz(dm,dm), aoper(dm,dm)
    if (dm <= 1 .or. mod(dm,2) /= 0) &
    call terminate('atnz2block called incorrectly')
    aoper = atnz
    do fi = 1, dm
      do fj = 1, dm
        if (fi <= dm/2) then
          if (fj <= dm/2) then
            atnz(fi,fj) = aoper(2*fi-1,2*fj-1)
          else
            atnz(fi,fj) = aoper(2*fi-1,2*(fj-dm/2))
          end if
        else
          if (fj <= dm/2) then
            atnz(fi,fj) = aoper(2*(fi-dm/2)-1,2*fj)
          else
            atnz(fi,fj) = aoper(2*(fi-dm/2),2*(fj-dm/2))
          end if
        end if
      end do
    end do
  end subroutine atnz2block

!------------------------------------------------------------
!> real matrix and matrix multiplicity
!!
!! makesure: dima(2) = dimb(1)  dimc(1) = dima(1)  dimc(2) = dimb(2)
  subroutine real_matmul(transa, transb, mata, matb, matc)
    implicit none
    character(len=1),intent(in) :: transa
    character(len=1),intent(in) :: transb
    real(dp),intent(in)         :: mata(:,:)
    real(dp),intent(in)         :: matb(:,:)
    real(dp),intent(out)        :: matc(:,:)
    integer                     :: dima(2)
    integer                     :: dimb(2)
    integer                     :: dimc(2)
    integer                     :: lda, ldb, swt
    dima = shape(mata)
    lda = dima(1)
    if (transa == 'T' .or. transa == 't') then
      swt = dima(1)
      dima(1) = dima(2)
      dima(2) = swt
    end if
    dimb = shape(matb)
    ldb = dimb(1)
    if (transb == 'T' .or. transb == 't') then
      swt = dimb(1)
      dimb(1) = dimb(2)
      dimb(2) = swt
    end if
    if (dima(2) /= dimb(1)) call terminate('real_matmul: dima(2) /= dimb(1)')
    dimc = shape(matc)
    if (dimc(1) /= dima(1)) call terminate('real_matmul: dimc(1) /= dima(1)')
    if (dimc(2) /= dimb(2)) call terminate('real_matmul: dimc(2) /= dimb(2)')
    call dgemm(transa, transb, dima(1), dimb(2), dima(2), 1.0_dp, mata, lda, &
    matb, ldb, 0.0_dp, matc, dima(1))
  end subroutine real_matmul

!------------------------------------------------------------
!> real matrix and vector multiplicity
  subroutine real_matvec(transa, mata, vecb, vecc)
    implicit none
    character(len=1),intent(in) :: transa
    real(dp),intent(in)         :: mata(:,:)
    real(dp),intent(in)         :: vecb(:)
    real(dp),intent(out)        :: vecc(:)
    integer                     :: dima(2)
    integer                     :: dimb(1)
    integer                     :: dimc(1)
    integer                     :: lda, ldb, swt
    dima = shape(mata)
    lda = dima(1)
    if (transa == 'N' .or. transa == 'n') then
      ! do nothing
    else if (transa == 'T' .or. transa == 't') then
      swt = dima(1)
      dima(1) = dima(2)
      dima(2) = swt
    else
      call terminate("real_matvec: transa should be 'N' or 'T'")
    end if
    dimb = shape(vecb)
    ldb = dimb(1)
    if (dima(2) /= dimb(1)) call terminate('real_matvec: dima(2) /= dimb(1)')
    dimc = shape(vecc)
    if (dimc(1) /= dima(1)) call terminate('real_matvec: dimc(1) /= dima(1)')
    call dgemm(transa, 'N', dima(1), 1, dima(2), 1.0_dp, mata, lda, &
    vecb, ldb, 0.0_dp, vecc, dima(1))
  end subroutine real_matvec

!------------------------------------------------------------
!> real matrix and vector multiplicity
  subroutine cmplx_matvec(transa, mata, vecb, vecc)
    implicit none
    character(len=1),intent(in) :: transa
    complex(dp),intent(in)      :: mata(:,:)
    complex(dp),intent(in)      :: vecb(:)
    complex(dp),intent(out)     :: vecc(:)
    integer                     :: dima(2)
    integer                     :: dimb(1)
    integer                     :: dimc(1)
    integer                     :: lda, ldb, swt
    dima = shape(mata)
    lda = dima(1)
    if (transa == 'N' .or. transa == 'n') then
      ! do nothing
    else if (transa == 'T' .or. transa == 't') then
      swt = dima(1)
      dima(1) = dima(2)
      dima(2) = swt
    else if (transa == 'C' .or. transa == 'c') then
      swt = dima(1)
      dima(1) = dima(2)
      dima(2) = swt
    else
      call terminate("cmplx_matvec: transa should be 'N' or 'T' or 'C'")
    end if
    dimb = shape(vecb)
    ldb = dimb(1)
    if (dima(2) /= dimb(1)) call terminate('cmplx_matvec: dima(2) /= dimb(1)')
    dimc = shape(vecc)
    if (dimc(1) /= dima(1)) call terminate('cmplx_matvec: dimc(1) /= dima(1)')
    call zgemm(transa, 'N', dima(1), 1, dima(2), c1, mata, lda, &
    vecb, ldb, c0, vecc, dima(1))
  end subroutine cmplx_matvec

!------------------------------------------------------------
!> complex matrix and matrix multiplicity
!!
!! makesure: dima(2) = dimb(1)  dimc(1) = dima(1)  dimc(2) = dimb(2)
  subroutine cmplx_matmul(transa, transb, mata, matb, matc)
    implicit none
    character(len=1),intent(in) :: transa
    character(len=1),intent(in) :: transb
    complex(dp),intent(in)      :: mata(:,:)
    complex(dp),intent(in)      :: matb(:,:)
    complex(dp),intent(out)     :: matc(:,:)
    integer                     :: dima(2)
    integer                     :: dimb(2)
    integer                     :: dimc(2)
    integer                     :: lda, ldb, swt
    dima = shape(mata)
    lda = dima(1)
    if (transa == 'T' .or. transa == 't') then
      swt = dima(1)
      dima(1) = dima(2)
      dima(2) = swt
    else if (transa == 'C' .or. transa == 'c') then
      swt = dima(1)
      dima(1) = dima(2)
      dima(2) = swt
    end if
    dimb = shape(matb)
    ldb = dimb(1)
    if (transb == 'T' .or. transb == 't') then
      swt = dimb(1)
      dimb(1) = dimb(2)
      dimb(2) = swt
    else if (transb == 'C' .or. transb == 'c') then
      swt = dimb(1)
      dimb(1) = dimb(2)
      dimb(2) = swt
    end if
    if (dima(2) /= dimb(1)) call terminate('cmplx_matmul: dima(2) /= dimb(1)')
    dimc = shape(matc)
    if (dimc(1) /= dima(1)) call terminate('cmplx_matmul: dimc(1) /= dima(1)')
    if (dimc(2) /= dimb(2)) call terminate('cmplx_matmul: dimc(2) /= dimb(2)')
    call zgemm( transa, transb, dima(1), dimb(2), dima(2), c1, mata, lda, &
    matb, ldb, c0, matc, dima(1))
  end subroutine cmplx_matmul

!------------------------------------------------------------
!> Frobenius norm of a real matrix
  real(dp) function real_fbnorm(mat) result(fbnorm)
    implicit none
    real(dp),intent(in)          :: mat(:,:)
    integer                      :: ii, jj
    integer                      :: dm(2)
    fbnorm = 0.0_dp
    dm = shape(mat)
    do ii = 1, dm(1)
      do jj = 1, dm(2)
        fbnorm = fbnorm + mat(ii,jj)**2
      end do
    end do
    fbnorm = dsqrt(fbnorm)
  end function real_fbnorm

!------------------------------------------------------------
!> Frobenius norm of a complex matrix
  real(dp) function cmplx_fbnorm(mat) result(fbnorm)
    implicit none
    complex(dp),intent(in)       :: mat(:,:)
    integer                      :: ii, jj
    integer                      :: dm(2)
    fbnorm = 0.0_dp
    dm = shape(mat)
    do ii = 1, dm(1)
      do jj = 1, dm(2)
        fbnorm = fbnorm + real(mat(ii,jj)*conjg(mat(ii,jj)),dp)
      end do
    end do
    fbnorm = dsqrt(fbnorm)
  end function cmplx_fbnorm

!------------------------------------------------------------
!> calculate determinant of a real matrix
  real(dp) function real_det(mat) result(val)
    implicit none
    real(dp), intent(in)    :: mat(:,:)
    real(dp), allocatable   :: LU(:,:)
    integer, allocatable    :: ipiv(:)
    integer                 :: n, i
    real(dp)                :: sign
    
    n = size(mat, 1)
    if (size(mat, 2) /= n) call terminate("real_det: mat must be square")
    allocate(LU(n, n), ipiv(n))
    LU = mat ! have to copy because dgetrf will replace mat
    ! det(A)=det(P^−1LU)=det(P^−1)det(L)det(U)
    !                   =(−1)^k*U_11*U_22*...*U_nn
    call dgetrf(n, n, LU, n, ipiv, info)  ! LU factorization
    if (info > 0) then      ! mat is singular matrix
      val = 0.0_dp
    else if (info < 0) then ! invalid argument
      call terminate("real_det: invalid argument in dgetrf")
    else
      val = 1.0_dp
      sign = 1.0_dp
      do i = 1, n
        val = val * LU(i, i)
        if (ipiv(i) /= i) sign = -sign
      end do
      val = sign * val
    end if
  end function real_det

!------------------------------------------------------------
!> calculate determinant of a matrix matrix
  complex(dp) function cmplx_det(mat) result(val)
    implicit none
    complex(dp), intent(in)    :: mat(:,:)
    complex(dp), allocatable   :: LU(:,:)
    integer, allocatable       :: ipiv(:)
    integer                    :: n, i
    real(dp)                   :: sign
    
    n = size(mat, 1)
    if (size(mat, 2) /= n) call terminate("cmplx_det: mat must be square")
    allocate(LU(n, n), ipiv(n))
    LU = mat ! have to copy because zgetrf will replace mat
    ! det(A)=det(P^−1LU)=det(P^−1)det(L)det(U)
    !                   =(−1)^k*U_11*U_22*...*U_nn
    call zgetrf(n, n, LU, n, ipiv, info)  ! LU factorization
    if (info > 0) then      ! mat is singular matrix
      val = c0
    else if (info < 0) then ! invalid argument
      call terminate("cmplx_det: invalid argument in zgetrf")
    else
      val = c1
      sign = 1.0_dp
      do i = 1, n
        val = val * LU(i, i)
        if (ipiv(i) /= i) sign = -sign
      end do
      val = sign * val
    end if
  end function cmplx_det

!------------------------------------------------------------
!> calculate (reduced) Wigner d-matrix d^S_MK(theta) for specific quantum number
!!
!! |S,M> = exp(-i*phi*M)*d^S_MK*exp(-i*chi*K) |S,K>
!!
!! S,M,K and specific angle theta based on z-y-z convention of Euler angles.
!!
!! !!Note: the input S,M,K is actually 2S,2M,2K respectively
  pure recursive real(dp) function wigner_d(S, M, K, theta) result(smalld)
    implicit none
    integer,intent(in)  ::  S          ! spin quantum number S/2
    integer,intent(in)  ::  M          ! spin magnetic quantum number M/2
    integer,intent(in)  ::  K          ! spin magnetic quantum number K/2
    real(dp),intent(in) ::  theta      ! input angle
    integer             ::  SpM        ! S+M
    integer             ::  SmM        ! S-M
    integer             ::  SpK        ! S+K
    integer             ::  SmK        ! S-K
    integer             ::  MmK        ! M-K
    integer             ::  vs, vsmin, vsmax
    real(dp)            ::  coe, vsd

    if (S < abs(M) .or. S < abs(K)) then
      smalld = 0.0_dp
      return
    end if
    ! general formula
    if (S > 4) then
      SpM = nint(0.5_dp*(real(S,dp)+real(M,dp)))
      SmM = nint(0.5_dp*(real(S,dp)-real(M,dp)))
      SpK = nint(0.5_dp*(real(S,dp)+real(K,dp)))
      SmK = nint(0.5_dp*(real(S,dp)-real(K,dp)))
      MmK = nint(0.5_dp*(real(M,dp)-real(K,dp)))
      coe = dsqrt(factorial(SpM)*factorial(SmM)*factorial(SpK)*factorial(SmK))
      vsmax = min(SpK, SmM)
      vsmin = max(0, -MmK)
      smalld = 0.0_dp
      do vs = vsmin, vsmax
        vsd = (-1.0_dp)**(MmK+vs) * dcos(0.5_dp*theta)**(SpK+SmM-2*vs) * &
        dsin(0.5_dp*theta)**(MmK+2*vs)
        vsd = vsd / &
        (factorial(SpK-vs)*factorial(vs)*factorial(MmK+vs)*factorial(SmM-vs))
        smalld = smalld + vsd
      end do
      smalld = smalld * coe
      return
    end if
    ! d^S_MK = (-1)^(K-M)*d^S_KM = d^S_-K-M
    if (.not.(M >= 0 .and. M >= abs(K))) then
      if (M >= 0) then
        if (K > 0) then        ! K > M >= 0
          smalld = (-1.0_dp)**((K-M)/2)*wigner_d(S, K, M, theta)
        else if (K < 0) then   ! M >= 0 > K  (e.g. M=1, K=-2)
          smalld = wigner_d(S, -K, -M, theta)
        end if
      else if (M < 0) then
        if (K > 0) then        ! K > 0 > M
          smalld = (-1.0_dp)**((K-M)/2)*wigner_d(S, K, M, theta)
        else if (K <= 0) then  ! M < 0  k <= 0
          smalld = wigner_d(S, -K, -M, theta)
        end if
      end if
      return
    end if
    select case(S)
      case(0) ! S = 0
        smalld = 1.0_dp                           ! M = K = 0
      case(1) ! S = 1/2
        if (M == 1 .and. K == 1) then             ! M = 1/2, K = 1/2
          smalld = dcos(0.5_dp*theta)
        else if (M == 1 .and. K == -1) then       ! M = 1/2, K = -1/2
          smalld = -dsin(0.5_dp*theta)
        end if
      case(2) ! S = 1
        if (M == 2 .and. K == 2) then             ! M = 1, K = 1
          smalld = 0.5_dp*(1.0_dp+dcos(theta))
        else if (M == 2 .and. K == 0) then        ! M = 1, K = 0
          smalld = -dsin(theta)/dsqrt(2.0_dp)
        else if (M == 2 .and. K == -2) then       ! M = 1, K = -1
          smalld = 0.5_dp*(1.0_dp-dcos(theta))
        else if (M == 0 .and. K == 0) then        ! M = 0, K = 0
          smalld = dcos(theta)
        end if
      case(3) ! S = 3/2
        if (M == 3 .and. K == 3) then             ! M = 3/2, K = 3/2
          smalld = 0.5_dp*(1.0_dp+dcos(theta))*dcos(0.5_dp*theta)
        else if (M == 3 .and. K == 1) then        ! M = 3/2, K = 1/2
          smalld = -0.5_dp*dsqrt(3.0_dp)*(1.0_dp+dcos(theta))*dsin(0.5_dp*theta)
        else if (M == 3 .and. K == -1) then       ! M = 3/2, K = -1/2
          smalld = 0.5_dp*dsqrt(3.0_dp)*(1.0_dp-dcos(theta))*dcos(0.5_dp*theta)
        else if (M == 3 .and. K == -3) then       ! M = 3/2, K = -3/2
          smalld = -0.5_dp*(1.0_dp-dcos(theta))*dsin(0.5_dp*theta)
        else if (M == 1 .and. K == 1) then        ! M = 1/2, K = 1/2
          smalld = 0.5_dp*(3.0_dp*dcos(theta)-1.0_dp)*dcos(0.5_dp*theta)
        else if (M == 1 .and. K == -1) then       ! M = 1/2, K = -1/2
          smalld = -0.5_dp*(3.0_dp*dcos(theta)+1.0_dp)*dsin(0.5_dp*theta)
        end if
      case(4) ! S = 2
        if (M == 4 .and. K == 4) then             ! M = 2, K = 2
          smalld = 0.25_dp*(1.0_dp+dcos(theta))**2
        else if (M == 4 .and. K == 2) then        ! M = 2, K = 1
          smalld = -0.5_dp*dsin(theta)*(1.0_dp+dcos(theta))
        else if (M == 4 .and. K == 0) then        ! M = 2, K = 0
          smalld = dsqrt(0.375_dp)*dsin(theta)**2
        else if (M == 4 .and. K == -2) then       ! M = 2, K = -1
          smalld = -0.5_dp*dsin(theta)*(1.0_dp-dcos(theta))
        else if (M == 4 .and. K == -4) then       ! M = 2, K = -2
          smalld = 0.25_dp*(1.0_dp-dcos(theta))**2
        else if (M == 2 .and. K == 2) then        ! M = 1, K = 1
          smalld = 0.5_dp*(2.0_dp*dcos(theta)**2+dcos(theta)-1.0_dp)
        else if (M == 2 .and. K == 0) then        ! M = 1, K = 0
          smalld = -dsqrt(0.375_dp)*dsin(2.0_dp*theta)
        else if (M == 2 .and. K == -2) then       ! M = 1, K = -1
          smalld = 0.5_dp*(-2.0_dp*dcos(theta)**2+dcos(theta)+1.0_dp)
        else if (M == 0 .and. K == 0) then        ! M = 0, K = 0
          smalld = 0.5_dp*(3.0_dp*dcos(theta)**2-1.0_dp)
        end if
    end select
  end function wigner_d

!-----------------------------------------------------------------------
!> fast factorial of given a
  pure elemental real(dp) function factorial(a)
    implicit none
    integer,intent(in) :: a
    integer :: na
    select case (a)
    case (0)
      factorial = 1.0_dp
    case (1)
      factorial = 1.0_dp
    case (2)
      factorial = 2.0_dp
    case (3)
      factorial = 6.0_dp
    case (4)
      factorial = 24.0_dp
    case (5)
      factorial = 120.0_dp
    case (6)
      factorial = 720.0_dp
    case (7)
      factorial = 5054.0_dp
    case (8)
      factorial = 40320.0_dp
    case (9)
      factorial = 362880.0_dp
    case (10)
      factorial = 3628800.0_dp
    case (11)
      factorial = 39916800.0_dp
    case (12)
      factorial = 479001600.0_dp
    case (13)
      factorial = 6227020800.0_dp
    case (14)
      factorial = 87178291200.0_dp
    case (15)
      factorial = 1307674368000.0_dp
    case (16)
      factorial = 20922789888000.0_dp
    end select
    if (a <= 16) then
      return
    else
      factorial = 20922789888000.0_dp
      do na = 17, a
        factorial = factorial * real(na)
      end do
    end if
    return
  end function factorial
  
!-----------------------------------------------------------------------
!> fast bionomial coefficient of given N, M (C_N^M)
!!
!! M <= 9: return
!!
!! 9 < M <= 13: recursive
!!
!! M > 13: direct
  pure recursive real(dp) function binomialcoe(N,M) result(bicoe)
    implicit none
    integer,intent(in) :: N,M
    integer :: fi, numerator, denominator
    if (M == 0 .or. N == 0 .or. M == N) then
      bicoe = 1.0
      return
    else if (M == 1) then
      bicoe = real(N,dp)
      return
    else if (2*M > N) then
      bicoe = binomialcoe(N,N-M)
      return
    end if
    select case (N)
    case (4)
      bicoe = 6.0_dp
    case (5)
      bicoe = 10.0_dp
    case (6)
      select case (M)
      case (2)
        bicoe = 15.0_dp
      case (3)
        bicoe = 20.0_dp
      end select
    case (7)
      select case (M)
      case (2)
        bicoe = 21.0_dp
      case (3)
        bicoe = 35.0_dp
      end select
    case (8)
      select case (M)
      case (2)
        bicoe = 28.0_dp
      case (3)
        bicoe = 56.0_dp
      case (4)
        bicoe = 70.0_dp
      end select
    case (9)
      select case (M)
      case (2)
        bicoe = 36.0_dp
      case (3)
        bicoe = 84.0_dp
      case (4)
        bicoe = 126.0_dp
      end select
    end select
    if (N <= 9) then
      return
    else if (N <= 13) then
      bicoe = binomialcoe(N-1,M) + binomialcoe(N-1,M-1)
      return
    else
      bicoe = 1.0
      numerator = 1
      fi = N
      do while(fi > M)
        numerator = numerator * fi
        fi = fi - 1
      end do
      denominator = factorial(N - M)
      bicoe = real(numerator) / real(denominator)
      return
    end if
  end function binomialcoe
end module Fundamentals