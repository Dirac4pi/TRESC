!> @file Fundamentals.f90
!!
!! @brief definitions, fundamental processes and matrix operations
!!
!! @syntax Fortran 2008 free format
!!
!! @code UTF-8
!!
!! @author dirac4pi


!! 去掉s_h，只支持spherical-harmonic基组，把Xm跟c2s组合成f2s减少矩阵操作
!! 进一步提升单双电子积分精度，尤其是高角动量部分，提升Taylor级数跟截断阈值
module Fundamentals
  use IFPORT
  use LAPACK95
  
!-----------------------------------------------------------------------
! definitions for programming
  
  integer             :: loop_i, loop_j, loop_k     ! universal loop variables
  integer             :: loop_m, loop_n, loop_l
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
  integer         :: DKH_order   = 2      ! 0: nonrelativistic; 1: fpFW; 2: DKH2
  logical(kind=4) :: SRTP_type   = .false.! Second Relativized Thomas Precession
  logical(kind=4) :: STTP_type   = .false.! Spin Tensor Thomas Precession
  logical         :: finitenuc   = .false.
  real(dp)        :: cutS        = 1E-5   ! threshold of evl(i_j)
  !--------------------<module Atoms>--------------------
  integer         :: charge      = 0      ! charge of the system
  integer         :: spin_mult   = 675    ! spin multiplicity of the system
  integer         :: electron_count       ! number of electrons of the system
  logical         :: s_h         = .true. ! Cartesian / spher-harmo basis
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
  logical         :: keepspin    = .false.! avoid spin mutations
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


  private :: dump_matrix_real, dump_matrix_cmplx, load_matrix_real
  private :: load_matrix_cmplx, real_symm_inverse, cmplx_inverse
  private :: real_symm_diag, cmplx_diag, real_matmul, cmplx_matmul
  private :: real_matvec, cmplx_matvec
  public  :: dump_matrix, load_matrix, generate_output, terminate, matmul
  public  :: lowercase, is_alpha, diag, inverse, can_orth, symm_orth, atnz2block

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
    integer              :: dialoop_i, dialoop_j
    allocate(work(1))
    allocate(iwork(1))
    mat_u = 0.0_dp
    do dialoop_i = 1, dm
      do dialoop_j = dialoop_i, dm
        mat_u(dialoop_i,dialoop_j) = mat(dialoop_i,dialoop_j)
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
    integer                 :: dialoop_i, dialoop_j
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
    integer                 :: invloop_i, invloop_j
    mat_u = 0.0_dp
    do invloop_i = 1, dm
      do invloop_j = invloop_i, dm
        mat_u(invloop_i,invloop_j) = mat(invloop_i,invloop_j)
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
    do invloop_i = 1, dm
      do invloop_j = invloop_i, dm
        mat(invloop_j,invloop_i) = mat_u(invloop_i,invloop_j)
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
! symmetric orthogonalisation of given real symmetric matrix
! X mat XT = I
  
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
    integer              :: symloop_i
    call real_symm_diag(mat, dm, U, evl)
    min_evl = evl(1)
    mat_u = 0.0_dp
    do symloop_i = 1, dm
      if (evl(symloop_i) < 0.0) call terminate(&
      'symm_orth: eigenvalue less than zero, may due to code error')
      if (abs(evl(symloop_i)) < safmin) then
        write(60,'(A)') &
        '  --matrix is not invertible (not full rank), call can_orth directly'
        min_evl = 0.0_dp
        return
      end if
      if (min_evl > evl(symloop_i)) min_evl = evl(symloop_i)
      mat_u(symloop_i,symloop_i) = 1.0_dp / dsqrt(evl(symloop_i))
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
    integer                          :: canloop_i, canloop_j, canloop_k
    call real_symm_diag(mat, dm, U, evl)
    supU = 0
    dm2 = dm
    do canloop_i = 1, dm
      if (evl(canloop_i) < cutS) then
        supU(canloop_i) = 1
        dm2 = dm2 - 1
      end if
    end do
    allocate(X(dm,dm2))
    canloop_k = 1
    do canloop_i = 1, dm
      if (supU(canloop_i) == 1) cycle
      do canloop_j = 1, dm
        X(canloop_j,canloop_k) = U(canloop_j,canloop_i) / dsqrt(evl(canloop_i))
      end do
      canloop_k = canloop_k + 1
    end do
  end subroutine can_orth

!------------------------------------------------------------
!> transfer alternating zero matrix to block matrix
  subroutine atnz2block(atnz, dm)
    implicit none
    integer, intent(in) :: dm
    integer             :: aloop_i, aloop_j
    complex(dp)         :: atnz(dm,dm), aoper(dm,dm)
    if (dm <= 1 .or. mod(dm,2) /= 0) &
    call terminate('atnz2block called incorrectly')
    aoper = atnz
    do aloop_i = 1, dm
      do aloop_j = 1, dm
        if (aloop_i <= dm/2) then
          if (aloop_j <= dm/2) then
            atnz(aloop_i,aloop_j) = aoper(2*aloop_i-1,2*aloop_j-1)
          else
            atnz(aloop_i,aloop_j) = aoper(2*aloop_i-1,2*(aloop_j-dm/2))
          end if
        else
          if (aloop_j <= dm/2) then
            atnz(aloop_i,aloop_j) = aoper(2*(aloop_i-dm/2)-1,2*aloop_j)
          else
            atnz(aloop_i,aloop_j) = aoper(2*(aloop_i-dm/2),2*(aloop_j-dm/2))
          end if
        end if
      end do
    end do
  end subroutine atnz2block

!------------------------------------------------------------
!> real matrix and matrix multiplicity
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
    if (transa == 'N' .or. transa == 'n') then
      ! do nothing
    else if (transa == 'T' .or. transa == 't') then
      swt = dima(1)
      dima(1) = dima(2)
      dima(2) = swt
    else
      call terminate("real_matmul: transa should be 'N' or 'T'")
    end if
    dimb = shape(matb)
    ldb = dimb(1)
    if (transb == 'N' .or. transb == 'n') then
      ! do nothing
    else if (transb == 'T' .or. transb == 't') then
      swt = dimb(1)
      dimb(1) = dimb(2)
      dimb(2) = swt
    else
      call terminate("real_matmul: transb should be 'N' or 'T'")
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
      call terminate("cmplx_matmul: transa should be 'N' or 'T' or 'C'")
    end if
    dimb = shape(matb)
    ldb = dimb(1)
    if (transb == 'N' .or. transb == 'n') then
      ! do nothing
    else if (transb == 'T' .or. transb == 't') then
      swt = dimb(1)
      dimb(1) = dimb(2)
      dimb(2) = swt
    else if (transb == 'C' .or. transb == 'c') then
      swt = dimb(1)
      dimb(1) = dimb(2)
      dimb(2) = swt
    else
      call terminate("cmplx_matmul: transb should be 'N' or 'T' or 'C'")
    end if
    if (dima(2) /= dimb(1)) call terminate('cmplx_matmul: dima(2) /= dimb(1)')
    dimc = shape(matc)
    if (dimc(1) /= dima(1)) call terminate('cmplx_matmul: dimc(1) /= dima(1)')
    if (dimc(2) /= dimb(2)) call terminate('cmplx_matmul: dimc(2) /= dimb(2)')
    call zgemm( transa, transb, dima(1), dimb(2), dima(2), c1, mata, lda, &
    matb, ldb, c0, matc, dima(1))
  end subroutine cmplx_matmul

end module Fundamentals