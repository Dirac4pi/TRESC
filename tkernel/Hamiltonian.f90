!> @file Hamiltonian.f90
!!
!! @brief static electronic Hamiltonian integrals
!!
!! @syntax Fortran 2008 free format
!!
!! @code UTF-8
!!
!! @author dirac4pi
module Hamiltonian
  use Atoms
  use Fundamentals
  use LAPACK95
  use Rys
  use GRysroot
  use OMP_LIB

! <AOi|AOj> related
  !DIR$ ATTRIBUTES ALIGN:align_size :: i_j, i_j_s, M_i_j
  real(dp),allocatable    :: i_j(:,:)      ! <AOi|AOj>
  real(dp),allocatable    :: i_j_s(:,:)    ! <AOi|AOj> (sphe-har, reverse then)
  real(dp),allocatable    :: M_i_j(:,:)    ! <AOi|AOjm>

! <AOi|p^2|AOj> related
  !DIR$ ATTRIBUTES ALIGN:align_size :: i_p2_j, AO2p2, evl_p2, exi_T_j
  real(dp),allocatable    :: i_p2_j(:,:)   ! <AOi|p^2|AOj>
  real(dp),allocatable    :: i_T_j(:,:)    ! <AOi|T|AOj>
  real(dp),allocatable    :: i_T_j_read(:,:)    ! <AOi|T|AOj>
  ! unitary transformation from AO basis to p^2 eigenstate (validated)
  real(dp),allocatable    :: AO2p2(:,:)    ! transformation matrix from AO to p2
  ! all eigenvalues of <AOi|p^2|AOj> found by dsyevr
  real(dp),allocatable    :: evl_p2(:)     ! eigenvalue of i_p2_j
  complex(dp),allocatable :: exi_T_j(:,:)  ! extended i_p2_j matrix
  
! <AOi|V|AOj> related
  !DIR$ ATTRIBUTES ALIGN:align_size :: i_V_j, exi_V_j
  real(dp),allocatable    :: i_V_j(:,:)    ! <AOi|V|AOj>
  complex(dp),allocatable :: exi_V_j(:,:)  ! extended i_V_j matrix
  
! pVp related
  !DIR$ ATTRIBUTES ALIGN:align_size :: pxVpx, pyVpy, pzVpz
  !DIR$ ATTRIBUTES ALIGN:align_size :: pxVpy, pyVpx, pxVpz
  !DIR$ ATTRIBUTES ALIGN:align_size :: pzVpx, pyVpz, pzVpy, expVp
  real(dp),allocatable    :: pxVpx(:,:)    ! <AOi|pVp|AOj>
  real(dp),allocatable    :: pyVpy(:,:)
  real(dp),allocatable    :: pzVpz(:,:)
  real(dp),allocatable    :: pxVpy(:,:)
  real(dp),allocatable    :: pyVpx(:,:)
  real(dp),allocatable    :: pxVpz(:,:)
  real(dp),allocatable    :: pzVpx(:,:)
  real(dp),allocatable    :: pyVpz(:,:)
  real(dp),allocatable    :: pzVpy(:,:)
  complex(dp),allocatable :: expVp(:,:)    ! extended pVp-related matrix
  
! pppVp related
  !DIR$ ATTRIBUTES ALIGN:align_size :: px3Vpx, py3Vpy, pz3Vpz
  !DIR$ ATTRIBUTES ALIGN:align_size :: px3Vpy, py3Vpx, px3Vpz
  !DIR$ ATTRIBUTES ALIGN:align_size :: pz3Vpx, py3Vpz, pz3Vpy
  real(dp),allocatable    :: px3Vpx(:,:)   ! <AOi|pppVp|AOj>
  real(dp),allocatable    :: py3Vpy(:,:)
  real(dp),allocatable    :: pz3Vpz(:,:)
  real(dp),allocatable    :: px3Vpy(:,:)
  real(dp),allocatable    :: py3Vpx(:,:)
  real(dp),allocatable    :: px3Vpz(:,:)
  real(dp),allocatable    :: pz3Vpx(:,:)
  real(dp),allocatable    :: py3Vpz(:,:)
  real(dp),allocatable    :: pz3Vpy(:,:)
  ! <AOi|pxVpy3|AOj> = Trans(<AOi|py3Vpx|AOj>)

! SRTP related (momentum tensor)
  !DIR$ ATTRIBUTES ALIGN:align_size :: pxpx, pypy, pzpz
  !DIR$ ATTRIBUTES ALIGN:align_size :: pxpy, pypx, pxpz
  !DIR$ ATTRIBUTES ALIGN:align_size :: pzpx, pypz, pzpy
  real(dp),allocatable    :: pxpx(:,:)   ! <AOi|pipj|AOj>
  real(dp),allocatable    :: pxpy(:,:)
  real(dp),allocatable    :: pxpz(:,:)
  real(dp),allocatable    :: pypx(:,:)
  real(dp),allocatable    :: pypy(:,:)
  real(dp),allocatable    :: pypz(:,:)
  real(dp),allocatable    :: pzpx(:,:)
  real(dp),allocatable    :: pzpy(:,:)
  real(dp),allocatable    :: pzpz(:,:)


  ! interation Taylor expansion coefficients for Integral_V_1e and Integral_V_2e_OS
  !DIR$ ATTRIBUTES ALIGN:align_size :: intTaycoe
  real(dp)                :: intTaycoe(64,43)

  ! Per-thread workspace for batches of Obara-Saika two-electron integrals.
  ! It is resized only when a shell quartet requires a larger recurrence range.
  type :: V2eVecScratch
    real(dp),allocatable :: Gnm(:,:,:,:,:)
    real(dp),allocatable :: Itrans(:,:,:,:)
    real(dp),allocatable :: I(:,:,:)
    real(dp),allocatable :: PL(:,:)
  end type V2eVecScratch

  contains

!-----------------------------------------------------------------------
!> ensure that a vectorised Obara-Saika work array has the required shape
  recursive subroutine Reserve_V2eVecScratch(scratch, nlanes, nij, nkl, nrec, nl, npl)
    implicit none
    type(V2eVecScratch),intent(inout) :: scratch
    integer,intent(in) :: nlanes, nij, nkl, nrec, nl, npl
    logical :: resize

    resize = .not. allocated(scratch%Gnm)
    if (.not. resize) then
      resize = size(scratch%Gnm, 1) < nlanes .or. size(scratch%Gnm, 3) < nij .or. &
               size(scratch%Gnm, 4) < nkl .or. size(scratch%Gnm, 5) < nrec .or. &
               size(scratch%Itrans, 3) < nl .or. size(scratch%Itrans, 4) < nrec .or. &
               size(scratch%PL, 2) < npl
    end if
    if (.not. resize) return

    if (allocated(scratch%Gnm)) deallocate(scratch%Gnm)
    if (allocated(scratch%Itrans)) deallocate(scratch%Itrans)
    if (allocated(scratch%I)) deallocate(scratch%I)
    if (allocated(scratch%PL)) deallocate(scratch%PL)
    allocate(scratch%Gnm(nlanes, 3, nij, nkl, nrec))
    allocate(scratch%Itrans(nlanes, 3, nl, nrec))
    allocate(scratch%I(nlanes, 3, nrec))
    allocate(scratch%PL(nlanes, npl))
  end subroutine Reserve_V2eVecScratch

!-----------------------------------------------------------------------
!> release a thread-local vectorised Obara-Saika work array
  subroutine Release_V2eVecScratch(scratch)
    implicit none
    type(V2eVecScratch),intent(inout) :: scratch

    if (allocated(scratch%Gnm)) deallocate(scratch%Gnm)
    if (allocated(scratch%Itrans)) deallocate(scratch%Itrans)
    if (allocated(scratch%I)) deallocate(scratch%I)
    if (allocated(scratch%PL)) deallocate(scratch%PL)
  end subroutine Release_V2eVecScratch

!-----------------------------------------------------------------------
!> add one floating-point term using Neumaier compensated summation
  pure subroutine Neumaier_Add(sum, correction, term)
    implicit none
    real(dp),intent(inout) :: sum, correction
    real(dp),intent(in)    :: term
    real(dp)               :: updated

    updated = sum + term
    if (abs(sum) >= abs(term)) then
      correction = correction + (sum-updated) + term
    else
      correction = correction + (term-updated) + sum
    end if
    sum = updated
  end subroutine Neumaier_Add

!-------------------------------------------------------------------
!> calculate non-relativistic scalar one-electron integrals
!!
!! and prepare for SCF process
  subroutine Hartree_Hamiltonian()
    implicit none
    integer           :: si, sj            ! loop variable DKH_Hamiltonian
    character(len=30) :: ch30
    write(60,'(A)') 'Module Hamiltonian:'
    ! OpenMP set up
    write(60,'(A)') '  ----------<PARALLEL>----------'
    nproc = omp_get_num_procs()
    new_threads = min(max(1,threads), nproc)
    write(60,'(A,I3,A,I3,A,I3)') &
    '  threads requested:',threads,'; threads using:',new_threads,'; node nproc:',nproc
    if (new_threads == 1) then
      write(60,'(A)') '  Calculation will be performed serially.'
    else
      call getenv('KMP_AFFINITY',ch30)
      if (ch30 == '') then
        write(60,'(A)') '  KMP_AFFINITY = None'
      else
        write(60,'(A)') '  KMP_AFFINITY = '//trim(ch30)
      end if
    end if
    write(60,'(A)') '  ----------<HAMILTONIAN>----------'
    write(60,'(A)') '  spinor basis causes additional cost in scalar SCF.'
    !-----------------------------------------------
    ! 1e integral calculation
    write(60,'(A)') '  normalization of GTOs'
    call Calc_Ncoe(cbdata, cbdm)
    write(60,'(A)') '  complete!'
    write(60,'(A)') '  initialize Gaussian integration'
    call Integral_init()
    write(60,'(A)') '  complete!'
    write(60,'(A)') '  one-electron integral calculation'
    call Assign_matrices_1e()
    write(60,'(A)') '  complete! stored in:'
    write(60,'(A)') '  i_j, i_p2_j, i_V_j'
    ! cbdm -> sbdm -> fbdm
    write(60,"(A)") '  perform basis transformation'
    allocate(i_j_s(cbdm,cbdm))
    i_j_s = i_j
    call Assign_csf(i_j)
    call csgo(i_j_s)
    call cfgo(i_V_j)
    call cfgo(i_p2_j)
    if (fbdm == sbdm) then
      write(60,"(A)") '  complete! by symm_orth.'
    else
      write(60,"(A,I4)") '  complete! by can_orth. fbdm = ', fbdm
    end if
    write(60,'(A)') 'exit module Hamiltonian'
  end subroutine Hartree_Hamiltonian

!-------------------------------------------------------------------
!> calculate relativistic spinor one-electron integrals proposed by Hess
!!
!! (doi:10.1063/1.1515314, include pVp-related integrals)
!!
!! and prepare for SCF process
  subroutine Hess_Hamiltonian()
    implicit none
    integer           :: si            ! loop variable DKH_Hamiltonian
    character(len=30) :: ch30
    write(60,'(A)') 'Module Hamiltonian:'
    ! OpenMP set up
    write(60,'(A)') '  ----------<PARALLEL>----------'
    nproc = omp_get_num_procs()
    new_threads = min(max(1,threads), nproc)
    write(60,'(A,I3,A,I3,A,I3)') &
    '  threads requested:',threads,'; threads using:',new_threads,'; node nproc:',nproc
    if (new_threads == 1) then
      write(60,'(A)') '  Calculation will be performed serially.'
    else
      call getenv('KMP_AFFINITY',ch30)
      if (ch30 == '') then
        write(60,'(A)') '  KMP_AFFINITY = None'
      else
        write(60,'(A)') '  KMP_AFFINITY = '//trim(ch30)
      end if
    end if
    write(60,'(A)') '  ----------<HAMILTONIAN>----------'
    write(60,'(A)') &
    '  QED effect: radiative correction(c^-3, c^-4, spin-dependent).'
    write(60,'(A)') &
    '  1e DKH transformation: scalar terms up to c^-2 order, spin-'
    write(60,'(A)') &
    '  dependent terms up to c^-4 order.'
    write(60,'(A)') &
    '  Incompleteness of basis introduces error in 1e Fock.'
    !-----------------------------------------------
    ! 1e integral calculation
    write(60,'(A)') '  normalization of GTOs'
    call Calc_Ncoe(cbdata, cbdm)
    write(60,'(A)') '  complete!'
    write(60,'(A)') '  initialize Gaussian integration'
    call Integral_init()
    write(60,'(A)') '  complete!'
    write(60,'(A)') '  one-electron integral calculation'
    call Assign_matrices_1e()
    write(60,'(A)') '  complete! stored in:'
    write(60,'(A)') '  i_j, i_p2_j, i_V_j, i_pVp_j (9 matrices)'
    if (srtp) write(60,'(A)') '  i_pp_j (9 matrices)'
    ! cbdm -> sbdm -> fbdm
    write(60,"(A)") '  perform basis transformation'
    allocate(i_j_s(cbdm,cbdm))
    i_j_s = i_j
    call Assign_csf(i_j)
    call csgo(i_j_s)
    call cfgo(i_V_j)
    call cfgo(i_p2_j)
    call cfgo(pxVpx)
    call cfgo(pyVpy)
    call cfgo(pzVpz)
    call cfgo(pxVpy)
    call cfgo(pyVpx)
    call cfgo(pxVpz)
    call cfgo(pzVpx)
    call cfgo(pyVpz)
    call cfgo(pzVpy)
    if (srtp) then
      call cfgo(pxpx)
      call cfgo(pypy)
      call cfgo(pzpz)
      call cfgo(pxpy)
      call cfgo(pypx)
      call cfgo(pxpz)
      call cfgo(pzpx)
      call cfgo(pypz)
      call cfgo(pzpy)
    end if
    if (fbdm == sbdm) then
      write(60,"(A)") '  complete! by symm_orth.'
    else
      write(60,"(A,I4)") '  complete! by can_orth. fbdm = ', fbdm
    end if
    ! <AOi|p^2|AOj> diagonalization
    write(60,'(A)') '  <AOi|p^2|AOj> diagonalization'
    allocate(AO2p2(fbdm,fbdm))
    allocate(evl_p2(fbdm))
    call diag(i_p2_j, fbdm, AO2p2, evl_p2)
    do si = 1, fbdm
      if (evl_p2(si) < 0.0) &
      call terminate('evl(T) less than zero, may due to code error')
    end do
    ! (AO2p2)^T(i_p2_j)(AO2p2)=evl_p2
    write(60,'(A,I4,A)') '  complete!', fbdm, ' eigenvalues found.'
    write(60,'(A)') 'exit module Hamiltonian'
  end subroutine Hess_Hamiltonian
  
!-----------------------------------------------------------------------
!> Assign value to one-electron integral matrices
  subroutine Assign_matrices_1e()
    implicit none
    integer          :: i, j                 ! openMP parallel variable
    integer          :: contri               ! contr of atoMi, shelLi
    integer          :: faci(3)              ! xyz factor of |AOi>
    integer          :: contrj               ! contr of atoMj, shelLj
    integer          :: facj(3)              ! xyz factor of |AOj>
    real(dp)         :: expi(16)             ! expo of |AOi>
    real(dp)         :: expj(16)             ! expo of |AOj>
    real(dp)         :: coei(16)             ! coefficient of |AOi>
    real(dp)         :: coej(16)             ! coefficient of |AOj>
    real(dp)         :: codi(3)              ! coordinate of center of |AOi>
    real(dp)         :: codj(3)              ! coordinate of center of |AOj>
    real(dp)         :: coedx_i(32)          ! coefficient.derivative x.|AOi>
    real(dp)         :: coedy_i(32)          ! coefficient.derivative y.|AOi>
    real(dp)         :: coedz_i(32)          ! coefficient.derivative z.|AOi>
    real(dp)         :: coedx_j(32)          ! coefficient.derivative x.|AOj>
    real(dp)         :: coedy_j(32)          ! coefficient.derivative y.|AOj>
    real(dp)         :: coedz_j(32)          ! coefficient.derivative z.|AOj>
    integer          :: facdx_i(3,2)         ! x,y,z factor.derivative x.|AOi>
    integer          :: facdy_i(3,2)         ! x,y,z factor.derivative y.|AOi>
    integer          :: facdz_i(3,2)         ! x,y,z factor.derivative z.|AOi>
    integer          :: facdx_j(3,2)         ! x,y,z factor.derivative x.|AOj>
    integer          :: facdy_j(3,2)         ! x,y,z factor.derivative y.|AOj>
    integer          :: facdz_j(3,2)         ! x,y,z factor.derivative z.|AOj>
    integer          :: si, sj, sk, sl       ! loop variables Assign_matrices_1e
    type threadlocal    ! thread-local storage to avoid thread-sync overhead
      real(dp), allocatable :: i_j(:,:), i_V_j(:,:), i_p2_j(:,:), pxVpx(:,:),&
      pyVpy(:,:), pzVpz(:,:), pxVpy(:,:), pyVpx(:,:), pyVpz(:,:), pzVpy(:,:),&
      pxVpz(:,:), pzVpx(:,:), pxpx(:,:), pypy(:,:), pzpz(:,:),&
      pxpy(:,:), pypx(:,:), pxpz(:,:), pzpx(:,:), pypz(:,:), pzpy(:,:)
    end type
    type(threadlocal) :: tl
    allocate(i_j(cbdm,cbdm), i_V_j(cbdm,cbdm), i_p2_j(cbdm,cbdm), source=0.0_dp)
    if (pVp1e) then
      allocate(pxVpx(cbdm,cbdm),pyVpy(cbdm,cbdm),pzVpz(cbdm,cbdm),source=0.0_dp)
      allocate(pxVpy(cbdm,cbdm),pyVpx(cbdm,cbdm),pyVpz(cbdm,cbdm),source=0.0_dp)
      allocate(pzVpy(cbdm,cbdm),pxVpz(cbdm,cbdm),pzVpx(cbdm,cbdm),source=0.0_dp)
      if (srtp) then
        allocate(pxpx(cbdm,cbdm),pypy(cbdm,cbdm),source=0.0_dp)
        allocate(pzpz(cbdm,cbdm),pxpy(cbdm,cbdm),source=0.0_dp)
        allocate(pypx(cbdm,cbdm),pxpz(cbdm,cbdm),source=0.0_dp)
        allocate(pzpx(cbdm,cbdm),pypz(cbdm,cbdm),source=0.0_dp)
        allocate(pzpy(cbdm,cbdm),source=0.0_dp)
      end if
    end if
    ! parallel zone, running results consistent with serial
    !$omp parallel num_threads(new_threads) default(shared) private(i,j,si,sj,sk,&
    !$omp& sl,contri,faci,contrj,facj,expi,expj,coei,coej,codi,codj,&
    !$omp& facdx_i,facdy_i,facdz_i,coedx_i,coedy_i,coedz_i,facdx_j,&
    !$omp& facdy_j,facdz_j,coedx_j,coedy_j,coedz_j,tl) if(new_threads > 1)
    allocate(tl%i_j(cbdm,cbdm),tl%i_V_j(cbdm,cbdm), source=0.0_dp)
    allocate(tl%i_p2_j(cbdm,cbdm), source=0.0_dp)
    if (pVp1e) then
      allocate(tl%pxVpx(cbdm,cbdm),tl%pyVpy(cbdm,cbdm),source=0.0_dp)
      allocate(tl%pzVpz(cbdm,cbdm),source=0.0_dp)
      allocate(tl%pxVpy(cbdm,cbdm),tl%pyVpx(cbdm,cbdm),source=0.0_dp)
      allocate(tl%pyVpz(cbdm,cbdm),tl%pzVpy(cbdm,cbdm),source=0.0_dp)
      allocate(tl%pxVpz(cbdm,cbdm),tl%pzVpx(cbdm,cbdm),source=0.0_dp)
      if (srtp) then
        allocate(tl%pxpx(cbdm,cbdm),tl%pypy(cbdm,cbdm),source=0.0_dp)
        allocate(tl%pzpz(cbdm,cbdm),source=0.0_dp)
        allocate(tl%pzpy(cbdm,cbdm),tl%pxpy(cbdm,cbdm),source=0.0_dp)
        allocate(tl%pypx(cbdm,cbdm),tl%pxpz(cbdm,cbdm),source=0.0_dp)
        allocate(tl%pzpx(cbdm,cbdm),tl%pypz(cbdm,cbdm),source=0.0_dp)
      end if
    end if
    !$omp do schedule(dynamic, 5) collapse(2)
    do i = 1, cbdm
      do j = 1, cbdm
        si = i
        contri = cbdata(si) % contr
        faci   = cbdata(si) % fac
        expi(1:contri) = cbdata(si) % expo(1:contri)
        coei(1:contri) = cbdata(si) % Ncoe(1:contri)
        codi    = cbdata(si) % pos
        facdx_i = cbdata(si) % facdx
        facdy_i = cbdata(si) % facdy
        facdz_i = cbdata(si) % facdz
        coedx_i(1:2*contri) = cbdata(si) % coedx(1:2*contri)
        coedy_i(1:2*contri) = cbdata(si) % coedy(1:2*contri)
        coedz_i(1:2*contri) = cbdata(si) % coedz(1:2*contri)
        !---------------------------------------
        sj = j
        contrj = cbdata(sj) % contr
        facj   = cbdata(sj) % fac
        expj(1:contrj) = cbdata(sj) % expo(1:contrj)
        coej(1:contrj) = cbdata(sj) % Ncoe(1:contrj)
        codj    = cbdata(sj) % pos
        facdx_j = cbdata(sj) % facdx
        facdy_j = cbdata(sj) % facdy
        facdz_j = cbdata(sj) % facdz
        coedx_j(1:2*contrj) = cbdata(sj) % coedx(1:2*contrj)
        coedy_j(1:2*contrj) = cbdata(sj) % coedy(1:2*contrj)
        coedz_j(1:2*contrj) = cbdata(sj) % coedz(1:2*contrj)
        ! calc <AOi|V|AOj>
        tl%i_V_j(si,sj) = Calc_V_1e(contri,contrj,coei(1:contri),coej(1:contrj),&
        faci,facj,expi(1:contri),expj(1:contrj),codi,codj)
        if (pVp1e) then
          ! calc <AOi|pVp|AOj>,totally 9 matrices,6 matrices will be calculated
          tl%pxVpx(si,sj)=Calc_pVp_1e(faci(1),facj(1),contri,contrj, &
          coedx_i(1:2*contri),coedx_j(1:2*contrj),facdx_i,facdx_j,   &
          expi(1:contri),expj(1:contrj),codi,codj)
          tl%pyVpy(si,sj)=Calc_pVp_1e(faci(2),facj(2),contri,contrj, &
          coedy_i(1:2*contri),coedy_j(1:2*contrj),facdy_i,facdy_j,   &
          expi(1:contri),expj(1:contrj),codi,codj)
          tl%pzVpz(si,sj)=Calc_pVp_1e(faci(3),facj(3),contri,contrj, &
          coedz_i(1:2*contri),coedz_j(1:2*contrj),facdz_i,facdz_j,   &
          expi(1:contri),expj(1:contrj),codi,codj)
          tl%pxVpy(si,sj)=Calc_pVp_1e(faci(1),facj(2),contri,contrj, &
          coedx_i(1:2*contri),coedy_j(1:2*contrj),facdx_i,facdy_j,   &
          expi(1:contri),expj(1:contrj),codi,codj)
          tl%pyVpx(sj,si) = tl%pxVpy(si,sj)
          tl%pxVpz(si,sj)=Calc_pVp_1e(faci(1),facj(3),contri,contrj, &
          coedx_i(1:2*contri),coedz_j(1:2*contrj),facdx_i,facdz_j,   &
          expi(1:contri),expj(1:contrj),codi,codj)
          tl%pzVpx(sj,si) = tl%pxVpz(si,sj)
          tl%pyVpz(si,sj)=Calc_pVp_1e(faci(2),facj(3),contri,contrj, &
          coedy_i(1:2*contri),coedz_j(1:2*contrj),facdy_i,facdz_j,   &
          expi(1:contri),expj(1:contrj),codi,codj)
          tl%pzVpy(sj,si) = tl%pyVpz(si,sj)
        end if
        !-------------<basis overlap integral>-------------
        do sk = 1, contri
          do sl = 1, contrj
            tl%i_j(si,sj) = tl%i_j(si,sj) + Integral_S_1e(&
            coei(sk)*coej(sl),faci,facj,expi(sk),expj(sl),codi,codj)
          end do
        end do
        !-------------<momentum tensor integral>-------------
        if (srtp) then
          tl%pxpx(si,sj) = Calc_pp_1e(faci(1),facj(1),contri,contrj, &
          coedx_i(1:2*contri),coedx_j(1:2*contrj),facdx_i,facdx_j, &
          expi(1:contri),expj(1:contrj),codi,codj)
          tl%pypy(si,sj) = Calc_pp_1e(faci(2),facj(2),contri,contrj, &
          coedy_i(1:2*contri),coedy_j(1:2*contrj),facdy_i,facdy_j, &
          expi(1:contri),expj(1:contrj),codi,codj)
          tl%pzpz(si,sj) = Calc_pp_1e(faci(3),facj(3),contri,contrj, &
          coedz_i(1:2*contri),coedz_j(1:2*contrj),facdz_i,facdz_j, &
          expi(1:contri),expj(1:contrj),codi,codj)
          tl%pxpy(si,sj) = tl%pxpy(si,sj) + Calc_pp_1e(                 &
          faci(1),facj(2),contri,contrj,coedx_i(1:2*contri),coedy_j(1:2*contrj),&
          facdx_i,facdy_j,expi(1:contri),expj(1:contrj),codi,codj)
          tl%pxpz(si,sj) = tl%pxpz(si,sj) + Calc_pp_1e(                 &
          faci(1),facj(3),contri,contrj,coedx_i(1:2*contri),coedz_j(1:2*contrj),&
          facdx_i,facdz_j,expi(1:contri),expj(1:contrj),codi,codj)
          tl%pypz(si,sj) = tl%pypz(si,sj) + Calc_pp_1e(                 &
          faci(2),facj(3),contri,contrj,coedy_i(1:2*contri),coedz_j(1:2*contrj),&
          facdy_i,facdz_j,expi(1:contri),expj(1:contrj),codi,codj)
          tl%pypx(sj,si) = tl%pxpy(si,sj)
          tl%pzpx(sj,si) = tl%pxpz(si,sj)
          tl%pzpy(sj,si) = tl%pypz(si,sj)
          tl%i_p2_j(si,sj) = tl%pxpx(si,sj) + tl%pypy(si,sj) + tl%pzpz(si,sj)
        else
          tl%i_p2_j(si,sj) = Calc_pp_1e(faci(1),facj(1),contri,contrj, &
          coedx_i(1:2*contri),coedx_j(1:2*contrj),facdx_i,facdx_j, &
          expi(1:contri),expj(1:contrj),codi,codj) + &
          Calc_pp_1e(faci(2),facj(2),contri,contrj,coedy_i(1:2*contri), &
          coedy_j(1:2*contrj),facdy_i,facdy_j,expi(1:contri),expj(1:contrj),codi,codj) + &
          Calc_pp_1e(faci(3),facj(3),contri,contrj,coedz_i(1:2*contri), &
          coedz_j(1:2*contrj),facdz_i,facdz_j,expi(1:contri),expj(1:contrj),codi,codj)
        end if
      end do
    end do
    !$omp end do
    !-------------<thread sync>-------------
    !$omp critical
    i_j = i_j + tl%i_j
    i_V_j = i_V_j + tl%i_V_j
    i_p2_j = i_p2_j + tl%i_p2_j
    if (pVp1e) then
      pxVpx = pxVpx + tl%pxVpx
      pyVpy = pyVpy + tl%pyVpy
      pzVpz = pzVpz + tl%pzVpz
      pxVpy = pxVpy + tl%pxVpy
      pyVpx = pyVpx + tl%pyVpx
      pxVpz = pxVpz + tl%pxVpz
      pzVpx = pzVpx + tl%pzVpx
      pyVpz = pyVpz + tl%pyVpz
      pzVpy = pzVpy + tl%pzVpy
      if (srtp) then
        pxpx = pxpx + tl%pxpx
        pypy = pypy + tl%pypy
        pzpz = pzpz + tl%pzpz
        pxpy = pxpy + tl%pxpy
        pxpz = pxpz + tl%pxpz
        pypz = pypz + tl%pypz
        pypx = pypx + tl%pypx
        pzpx = pzpx + tl%pzpx
        pzpy = pzpy + tl%pzpy
      end if
    end if
    !$omp end critical
    ! free memory for parallel threads explicitly
    deallocate(tl%i_j, tl%i_V_j, tl%i_p2_j)
    if (pVp1e) then
      deallocate(tl%pxVpx, tl%pyVpy, tl%pzVpz)
      deallocate(tl%pxVpy, tl%pyVpx, tl%pxVpz)
      deallocate(tl%pzVpx, tl%pyVpz, tl%pzVpy)
      if (srtp) then
        deallocate(tl%pxpx, tl%pypy, tl%pzpz)
        deallocate(tl%pxpy, tl%pypx, tl%pxpz)
        deallocate(tl%pzpx, tl%pypz, tl%pzpy)
      end if
    end if
    !$omp end parallel
  end subroutine Assign_matrices_1e

!-----------------------------------------------------------------------
!> calculate <AOi|p^2|AOj>
  real(dp) pure function Calc_pp_1e(&
  ni,nj,contri,contrj,coei,coej,faci,facj,expi,expj,codi,codj) result(val)
    implicit none
    integer,intent(in)  :: ni              ! number of x/y/z in <i|
    integer,intent(in)  :: nj              ! number of x/y/z in <j|
    integer,intent(in)  :: contri          ! contr of |AOi>
    integer,intent(in)  :: contrj          ! contr of |AOj>
    real(dp),intent(in) :: coei(2*contri)  ! coefficients of |AOi>
    real(dp),intent(in) :: coej(2*contrj)  ! coefficients of |AOj>
    integer,intent(in)  :: faci(3,2)       ! xyz factor of |AOi>
    integer,intent(in)  :: facj(3,2)       ! xyz factor of |AOj>
    real(dp),intent(in) :: expi(2*contri)  ! expo of |AOi>
    real(dp),intent(in) :: expj(2*contrj)  ! expo of |AOj>
    real(dp),intent(in) :: codi(3)         ! centrol coordintates.|AOi>
    real(dp),intent(in) :: codj(3)         ! centrol coordintates.|AOj>
    integer             :: numi, numj
    integer             :: ti,tj,tk,tl     ! loop variables for Calc_pp_1e
    if (ni == 0) then
      numi = 1
    else
      numi = 2
    end if
    if (nj == 0) then
      numj = 1
    else
      numj = 2
    end if
    val = 0.0_dp
    ! center of potential atom set to zero
    do tk = 1, numi
      do ti = 1, contri
        do tl = 1, numj
          do tj = 1, contrj
            val = val + Integral_S_1e(     &
            coei((tk-1)*contri+ti)*        &
            coej((tl-1)*contrj+tj),        &
            faci(:,tk),                    &
            facj(:,tl),                    &
            expi(ti),                      &
            expj(tj),                      &
            codi,                          &
            codj)
          end do
        end do
      end do
    end do
  end function Calc_pp_1e

!-----------------------------------------------------------------------
!> calculate <AOi|V|AOj>
  real(dp) pure function Calc_V_1e(&
  contri,contrj,coei,coej,faci,facj,expi,expj,codi,codj) result(val)
    implicit none
    integer,intent(in)  :: contri          ! contr of |AOi>
    integer,intent(in)  :: contrj          ! contr of |AOj>
    real(dp),intent(in) :: coei(contri)    ! coefficients of |AOi>
    real(dp),intent(in) :: coej(contrj)    ! coefficients of |AOj>
    integer,intent(in)  :: faci(3)         ! xyz factor of |AOi>
    integer,intent(in)  :: facj(3)         ! xyz factor of |AOj>
    real(dp),intent(in) :: expi(contri)    ! expo of |AOi>
    real(dp),intent(in) :: expj(contrj)    ! expo of |AOj>
    real(dp),intent(in) :: codi(3)         ! centrol coordintates.|AOi>
    real(dp),intent(in) :: codj(3)         ! centrol coordintates.|AOj>
    real(dp)            :: bcodi(3)        ! extended coordintates.|AOi>
    real(dp)            :: bcodj(3)        ! extended coordintates.|AOj>
    real(dp)            :: Z_pot, R_pot
    real(dp)            :: codpot(3)
    integer             :: ti,bj,tpot      ! loop variables for Calc_V_1e
    val = 0.0_dp
    do tpot = 1, atom_count
      codpot = mol(tpot) % pos
      Z_pot = real(mol(tpot) % atom_number)
      R_pot = mol(tpot) % rad / fm2Bohr
      ! center of potential atom set to zero
      bcodi(:) = codi(:) - codpot(:)
      bcodj(:) = codj(:) - codpot(:)
      do ti = 1, contri
        do bj = 1, contrj
          val = val + Integral_V_1e(     &
          Z_pot,                         &
          coei(ti),                      &
          coej(bj),                      &
          faci,                          &
          facj,                          &
          expi(ti),                      &
          expj(bj),                      &
          bcodi,                         &
          bcodj,                         &
          R_pot)
        end do
      end do
    end do
  end function Calc_V_1e

!-----------------------------------------------------------------------
!> calculate (AOiAOj|V|AOkAOl)
  recursive real(dp) function Calc_V_2e(&
  contri,contrj,contrk,contrl,coei,coej,coek,coel,codi,codj,codk,codl,&
  codA,codB,A,B,Gij,Gkl,Gimij,Gimkl,faci,facj,fack,facl) result(val)
    implicit none
    integer,intent(in)  :: contri, contrj, contrk, contrl ! contraction of |AO>
    real(dp),intent(in) :: coei(contri), coej(contrj)     ! coefficients of |AO>
    real(dp),intent(in) :: coek(contrk), coel(contrl)     ! coefficients of |AO>
    real(dp),intent(in) :: codi(3),codj(3),codk(3),codl(3)! centrol of |AO>
    real(dp),intent(in) :: codA(3,contri,contrj)          ! PRISM parameters
    real(dp),intent(in) :: codB(3,contrk,contrl)          ! PRISM parameters
    real(dp),intent(in) :: A(contri,contrj)               ! PRISM parameters
    real(dp),intent(in) :: B(contrk,contrl)               ! PRISM parameters
    real(dp),intent(in) :: Gij(3,contri,contrj)           ! PRISM parameters
    real(dp),intent(in) :: Gkl(3,contrk,contrl)           ! PRISM parameters
    real(dp),intent(in) :: Gimij(3,contri,contrj)         ! PRISM parameters
    real(dp),intent(in) :: Gimkl(3,contrk,contrl)         ! PRISM parameters
    integer,intent(in)  :: faci(3),facj(3),fack(3),facl(3)! xyz factor of |AO>
    integer,parameter   :: v2e_vec_batch = 256
    type(V2eVecScratch) :: scratch
    real(dp)            :: coe(v2e_vec_batch), Avec(v2e_vec_batch), Bvec(v2e_vec_batch)
    real(dp)            :: codAvec(v2e_vec_batch,3), codBvec(v2e_vec_batch,3)
    real(dp)            :: Gvec(v2e_vec_batch,3), Gimvec(v2e_vec_batch,3)
    real(dp)            :: correction
    integer             :: i, j, uo, up, contr ! loop variables
    val = 0.0_dp
    correction = 0.0_dp
    contr = 0
    do i = 1, contri
      do j = 1, contrj
        do uo = 1, contrk
          do up = 1, contrl
            contr = contr + 1
            coe(contr) = coei(i)*coej(j)*coek(uo)*coel(up)
            Avec(contr) = A(i,j)
            Bvec(contr) = B(uo,up)
            codAvec(contr,:) = codA(:,i,j)
            codBvec(contr,:) = codB(:,uo,up)
            Gvec(contr,:) = Gij(:,i,j) + Gkl(:,uo,up)
            Gimvec(contr,:) = Gimij(:,i,j) * Gimkl(:,uo,up)
            if (contr == v2e_vec_batch) then
              call Integral_V_2e_OS_PRISM_vec(scratch, contr, coe, codi, codj, codk, codl, &
                codAvec, codBvec, Avec, Bvec, Gvec, Gimvec, faci, facj, fack, facl, val, correction)
              contr = 0
            end if
          end do
        end do
      end do
    end do
    if (contr > 0) then
      call Integral_V_2e_OS_PRISM_vec(scratch, contr, coe, codi, codj, codk, codl, &
        codAvec(1:contr,:), codBvec(1:contr,:), Avec(1:contr), Bvec(1:contr), Gvec(1:contr,:), &
        Gimvec(1:contr,:), faci, facj, fack, facl, val, correction)
    end if
    val = val + correction
  end function Calc_V_2e
  
!-----------------------------------------------------------------------
!> calculate <AOi|pVp|AOj>
  real(dp) pure function Calc_pVp_1e(&
  ni,nj,contri,contrj,coei,coej,faci,facj,expi,expj,codi,codj) result(val)
    implicit none
    integer,intent(in)    :: ni            ! number of x/y/z in |AOi>
    integer,intent(in)    :: nj            ! number of x/y/z in |AOj>
    integer,intent(in)    :: contri        ! contr of |AOi>
    integer,intent(in)    :: contrj        ! contr of |AOj>
    ! coefficients of first-order derivative of |AOi>
    real(dp),intent(in)   :: coei(2*contri)
    ! coefficients of first-order derivative of |AOj>
    real(dp),intent(in)   :: coej(2*contrj)
    integer,intent(in)    :: faci(3,2)     ! xyz factor of |AOi>
    integer,intent(in)    :: facj(3,2)     ! xyz factor of |AOj>
    real(dp),intent(in)   :: expi(contri)  ! expo of |AOi>
    real(dp),intent(in)   :: expj(contrj)  ! expo of |AOj>
    real(dp),intent(in)   :: codi(3)       ! centrol coordintates.|AOi>
    real(dp),intent(in)   :: codj(3)       ! centrol coordintates.|AOj>
    real(dp)              :: scodi(3)      ! extended coordintates.|AOi>
    real(dp)              :: scodj(3)      ! extended coordintates.|AOj>
    real(dp)              :: Z_pot, R_pot
    real(dp)              :: codpot(3)
    integer               :: numi, numj
    integer               :: ti,tj,tk,tl,tpot ! loop variables for Calc_pVp_1e
    if (ni == 0) then
      numi = 1
    else
      numi = 2
    end if
    if (nj == 0) then
      numj = 1
    else
      numj = 2
    end if
    val = 0.0_dp
    do tpot = 1, atom_count
      codpot = mol(tpot) % pos
      Z_pot = real(mol(tpot) % atom_number)
      R_pot = mol(tpot) % rad / fm2Bohr
      ! center of potential atom set to zero
      scodi(:) = codi(:) - codpot(:)
      scodj(:) = codj(:) - codpot(:)
      do tk = 1, numi
      do ti = 1, contri
        do tl = 1, numj
        do tj = 1, contrj
          val = val + Integral_V_1e(      &
          Z_pot,                          &
          coei((tk-1)*contri+ti),         &
          coej((tl-1)*contrj+tj),         &
          faci(:,tk),                     &
          facj(:,tl),                     &
          expi(ti),                       &
          expj(tj),                       &
          scodi,                          &
          scodj,                          &
          R_pot)
        end do
        end do
      end do
      end do
    end do
  end function Calc_pVp_1e

!-----------------------------------------------------------------------
!> calculate (AOiAOj|pVp|AOkAOl) = (pAOipAOj|V|AOkAOl)
  recursive real(dp) function Calc_pVp_2eij(&
  ni,nj,contri,contrj,contrk,contrl,coei,coej,coek,coel,codi,codj,codk,codl,&
  codA,codB,A,B,Gij,Gkl,Gimij,Gimkl,faci,facj,fack,facl) result(val)
    implicit none
    integer,intent(in)  :: ni                         ! number of x/y/z in |AOi>
    integer,intent(in)  :: nj                         ! number of x/y/z in |AOj>
    integer,intent(in)  :: contri, contrj, contrk, contrl ! contraction of |AO>
    ! coefficients of first-order derivative of |AOi> and |AOj>
    real(dp),intent(in) :: coei(2*contri), coej(2*contrj)
    real(dp),intent(in) :: coek(contrk), coel(contrl)     ! coefficients of |AO>
    real(dp),intent(in) :: codi(3),codj(3),codk(3),codl(3)! centrol of |AO>
    real(dp),intent(in) :: codA(3,contri,contrj)          ! PRISM parameters
    real(dp),intent(in) :: codB(3,contrk,contrl)          ! PRISM parameters
    real(dp),intent(in) :: A(contri,contrj)               ! PRISM parameters
    real(dp),intent(in) :: B(contrk,contrl)               ! PRISM parameters
    real(dp),intent(in) :: Gij(3,contri,contrj)           ! PRISM parameters
    real(dp),intent(in) :: Gkl(3,contrk,contrl)           ! PRISM parameters
    real(dp),intent(in) :: Gimij(3,contri,contrj)         ! PRISM parameters
    real(dp),intent(in) :: Gimkl(3,contrk,contrl)         ! PRISM parameters
    ! xyz factor of first-order derivative of |AOi> and |AOj>
    integer,intent(in)  :: faci(3,2),facj(3,2)
    integer,intent(in)  :: fack(3),facl(3)                ! xyz factor of |AO>
    integer,parameter   :: v2e_vec_batch = 256
    type(V2eVecScratch) :: scratch
    real(dp)            :: coe(v2e_vec_batch), Avec(v2e_vec_batch), Bvec(v2e_vec_batch)
    real(dp)            :: codAvec(v2e_vec_batch,3), codBvec(v2e_vec_batch,3)
    real(dp)            :: Gvec(v2e_vec_batch,3), Gimvec(v2e_vec_batch,3)
    real(dp)            :: correction
    integer             :: numi, numj
    integer             :: um, un, uo, up, contr          ! loop variables
    integer             :: ii, jj ! loop variables
    if (ni == 0) then
      numi = 1
    else
      numi = 2
    end if
    if (nj == 0) then
      numj = 1
    else
      numj = 2
    end if
    val = 0.0_dp
    correction = 0.0_dp
    contr = 0
    do ii = 1, numi
    do um = 1, contri
      do jj = 1, numj
      do un = 1, contrj
        do uo = 1, contrk
          do up = 1, contrl
            contr = contr + 1
            coe(contr) = coei((ii-1)*contri+um)*coej((jj-1)*contrj+un)*coek(uo)*coel(up)
            Avec(contr) = A(um,un)
            Bvec(contr) = B(uo,up)
            codAvec(contr,:) = codA(:,um,un)
            codBvec(contr,:) = codB(:,uo,up)
            Gvec(contr,:) = Gij(:,um,un) + Gkl(:,uo,up)
            Gimvec(contr,:) = Gimij(:,um,un) * Gimkl(:,uo,up)
            if (contr == v2e_vec_batch) then
              call Integral_V_2e_OS_PRISM_vec(scratch, contr, coe, codi, codj, codk, codl, &
                codAvec, codBvec, Avec, Bvec, Gvec, Gimvec, faci(:,ii), facj(:,jj), fack, facl, val, correction)
              contr = 0
            end if
          end do
        end do
      end do
      if (contr > 0) then
        call Integral_V_2e_OS_PRISM_vec(scratch, contr, coe, codi, codj, codk, codl, &
          codAvec(1:contr,:), codBvec(1:contr,:), Avec(1:contr), Bvec(1:contr), Gvec(1:contr,:), &
          Gimvec(1:contr,:), faci(:,ii), facj(:,jj), fack, facl, val, correction)
        contr = 0
      end if
      end do
      end do
    end do
    val = val + correction
  end function Calc_pVp_2eij

!-----------------------------------------------------------------------
!> calculate (AOiAOj|pVp|AOkAOl) = (pAOiAOj|V|pAOkAOl)
  recursive real(dp) function Calc_pVp_2eik(&
  ni,nk,contri,contrj,contrk,contrl,coei,coej,coek,coel,codi,codj,codk,codl,&
  codA,codB,A,B,Gij,Gkl,Gimij,Gimkl,faci,facj,fack,facl) result(val)
    implicit none
    integer,intent(in)  :: ni                         ! number of x/y/z in |AOi>
    integer,intent(in)  :: nk                         ! number of x/y/z in |AOk>
    integer,intent(in)  :: contri, contrj, contrk, contrl ! contraction of |AO>
    ! coefficients of first-order derivative of |AOi> and |AOk>
    real(dp),intent(in) :: coei(2*contri), coek(2*contrk)
    real(dp),intent(in) :: coej(contrj), coel(contrl)     ! coefficients of |AO>
    real(dp),intent(in) :: codi(3),codj(3),codk(3),codl(3)! centrol of |AO>
    real(dp),intent(in) :: codA(3,contri,contrj)          ! PRISM parameters
    real(dp),intent(in) :: codB(3,contrk,contrl)          ! PRISM parameters
    real(dp),intent(in) :: A(contri,contrj)               ! PRISM parameters
    real(dp),intent(in) :: B(contrk,contrl)               ! PRISM parameters
    real(dp),intent(in) :: Gij(3,contri,contrj)           ! PRISM parameters
    real(dp),intent(in) :: Gkl(3,contrk,contrl)           ! PRISM parameters
    real(dp),intent(in) :: Gimij(3,contri,contrj)         ! PRISM parameters
    real(dp),intent(in) :: Gimkl(3,contrk,contrl)         ! PRISM parameters
    ! xyz factor of first-order derivative of |AOi> and |AOk>
    integer,intent(in)  :: faci(3,2),fack(3,2)
    integer,intent(in)  :: facj(3),facl(3)                ! xyz factor of |AO>
    integer,parameter   :: v2e_vec_batch = 256
    type(V2eVecScratch) :: scratch
    real(dp)            :: coe(v2e_vec_batch), Avec(v2e_vec_batch), Bvec(v2e_vec_batch)
    real(dp)            :: codAvec(v2e_vec_batch,3), codBvec(v2e_vec_batch,3)
    real(dp)            :: Gvec(v2e_vec_batch,3), Gimvec(v2e_vec_batch,3)
    real(dp)            :: correction
    integer             :: numi, numk
    integer             :: um, un, uo, up, contr          ! loop variables
    integer             :: ii, kk ! loop variables
    if (ni == 0) then
      numi = 1
    else
      numi = 2
    end if
    if (nk == 0) then
      numk = 1
    else
      numk = 2
    end if
    val = 0.0_dp
    correction = 0.0_dp
    contr = 0
    do ii = 1, numi
    do um = 1, contri
      do un = 1, contrj
        do kk = 1, numk
        do uo = 1, contrk
          do up = 1, contrl
            contr = contr + 1
            coe(contr) = coei((ii-1)*contri+um)*coej(un)*coek((kk-1)*contrk+uo)*coel(up)
            Avec(contr) = A(um,un)
            Bvec(contr) = B(uo,up)
            codAvec(contr,:) = codA(:,um,un)
            codBvec(contr,:) = codB(:,uo,up)
            Gvec(contr,:) = Gij(:,um,un) + Gkl(:,uo,up)
            Gimvec(contr,:) = Gimij(:,um,un) * Gimkl(:,uo,up)
            if (contr == v2e_vec_batch) then
              call Integral_V_2e_OS_PRISM_vec(scratch, contr, coe, codi, codj, codk, codl, &
                codAvec, codBvec, Avec, Bvec, Gvec, Gimvec, faci(:,ii), facj, fack(:,kk), facl, val, correction)
              contr = 0
            end if
          end do
        end do
        if (contr > 0) then
          call Integral_V_2e_OS_PRISM_vec(scratch, contr, coe, codi, codj, codk, codl, &
            codAvec(1:contr,:), codBvec(1:contr,:), Avec(1:contr), Bvec(1:contr), Gvec(1:contr,:), &
            Gimvec(1:contr,:), faci(:,ii), facj, fack(:,kk), facl, val, correction)
          contr = 0
        end if
        end do
      end do
      end do
    end do
    val = val + correction
  end function Calc_pVp_2eik
  
!-----------------------------------------------------------------------
!> calculate <AOi|p^3Vp|AOj> (9 matrices)
!!
!! px3Vpx py3Vpy pz3Vpz px3Vpy py3Vpx px3Vpz pz3Vpx py3Vpz pz3Vpy
  real(dp) pure function Calc_pppVp_1e(&
  di,ni,nj,contri,contrj,coei,coej,faci,facj,expi,expj,codi,codj) result(val)
    implicit none
    integer,intent(in)  :: di              ! 1:d|AOi>/dx, 1:d|AOi>/dy, 1:d|AOi>/dz
    integer,intent(in)  :: ni              ! number of x/y/z in <i|
    integer,intent(in)  :: nj              ! number of x/y/z in <j|
    integer,intent(in)  :: contri          ! contr of |AOi>
    integer,intent(in)  :: contrj          ! contr of |AOj>
    real(dp),intent(in) :: coei(2*contri)  ! coefficients of |AOi>
    real(dp),intent(in) :: coej(2*contrj)  ! coefficients of |AOj>
    integer,intent(in)  :: faci(3,2)       ! xyz factor of |AOi>
    integer,intent(in)  :: facj(3,2)       ! xyz factor of |AOj>
    real(dp),intent(in) :: expi(contri)    ! expo of |AOi>
    real(dp),intent(in) :: expj(contrj)    ! expo of |AOj>
    real(dp),intent(in) :: codi(3)         ! centrol coordintates.|AOi>
    real(dp),intent(in) :: codj(3)         ! centrol coordintates.|AOj>
    ! coordintates.|AOi> after 3 derivative actions
    real(dp)            :: tcodi(3)
    ! coordintates.|AOj> after 3 derivative actions
    real(dp)            :: tcodj(3)
    real(dp)            :: Z_pot, R_pot
    real(dp)            :: codpot(3)
    ! coefficients of |AOi> after 3 derivative actions
    real(dp)            :: tcoei(8*contri)
    ! xyz factor.|AOi> after 3 derivative actions
    integer             :: tfaci(3,8)      ! di**max, di**max-2 ... di**0
    integer             :: numi, numj      ! number of calls to Integral_V_1e
    integer             :: ti,tj,tk,tl,tpot! loop variables for Calc_pppVp_1e

    if (nj == 0) then
      numj = 1
    else
      numj = 2
    end if
  ! combine polynomials and ignore terms with coefficient 0
    select case(ni)
      case(0)
        tfaci(:,[1,2]) = faci(:,[2,2])         ! 1(1); 2(2,3);
        tfaci(di,1) = 3
        tfaci(di,2) = 1
        tcoei(1:contri) = 8*expi(:)**2*coei(1:contri)
        tcoei(contri+1:2*contri) = -6*expi(:)*coei(1:contri)
        numi = 2
      case(1)
        tfaci(:,[1,2,3]) = faci(:,[2,2,2])     ! 1(1); 2(2,3,5); 3(4,6);
        tfaci(di,1) = 4
        tfaci(di,2) = 2
        tfaci(di,3) = 0
        tcoei(1:contri) = 12*expi(:)**2*coei(1:contri)
        tcoei(contri+1:2*contri) = -10*expi(:)*coei(1:contri) + &
        4*expi(:)**2*coei(contri+1:2*contri)
        tcoei(2*contri+1:3*contri) = 2*coei(1:contri) - &
        2*expi(:)*coei(contri+1:2*contri)
        numi = 3
      case(2)
        tfaci(:,[1,2,3]) = faci(:,[2,2,2])     ! 1(1); 2(2,3,5); 3(4,6,7);
        tfaci(di,1) = 5
        tfaci(di,2) = 3
        tfaci(di,3) = 1
        tcoei(1:contri) = 16*expi(:)**2*coei(1:contri)
        tcoei(contri+1:2*contri) = -14*expi(:)*coei(1:contri) + &
        8*expi(:)**2*coei(contri+1:2*contri)
        tcoei(2*contri+1:3*contri) = 6*coei(1:contri) - &
        6*expi(:)*coei(contri+1:2*contri)
        numi = 3
      case(3)
        tfaci(:,[1,2,3,4]) = faci(:,[2,2,2,2]) ! 1(1); 2(2,3,5); 3(4,6,7); 4(8)
        tfaci(di,1) = 6
        tfaci(di,2) = 4
        tfaci(di,3) = 2
        tfaci(di,4) = 0
        tcoei(1:contri) = 20*expi(:)**2*coei(1:contri)
        tcoei(contri+1:2*contri) = -18*expi(:)*coei(1:contri) + &
        12*expi(:)**2*coei(contri+1:2*contri)
        tcoei(2*contri+1:3*contri) = 12*coei(1:contri) - &
        10*expi(:)*coei(contri+1:2*contri)
        tcoei(3*contri+1:4*contri) = 2*coei(contri+1:2*contri)
        numi = 4
      case(4)
        tfaci(:,[1,2,3,4]) = faci(:,[2,2,2,2]) ! 1(1); 2(2,3,5); 3(4,6,7); 4(8)
        tfaci(di,1) = 7
        tfaci(di,2) = 5
        tfaci(di,3) = 3
        tfaci(di,4) = 1
        tcoei(1:contri) = 24*expi(:)**2*coei(1:contri)
        tcoei(contri+1:2*contri) = -22*expi(:)*coei(1:contri) + &
        16*expi(:)**2*coei(contri+1:2*contri)
        tcoei(2*contri+1:3*contri) = 20*coei(1:contri) - &
        14*expi(:)*coei(contri+1:2*contri)
        tcoei(3*contri+1:4*contri) = 6*coei(contri+1:2*contri)
        numi = 4
    end select
    val = 0.0_dp
    do tpot = 1, atom_count
      codpot = mol(tpot) % pos
      Z_pot = real(mol(tpot) % atom_number)
      R_pot = mol(tpot) % rad / fm2Bohr
      ! center of potential atom set to zero
      tcodi(:) = codi(:) - codpot(:)
      tcodj(:) = codj(:) - codpot(:)
      do tk = 1, numi
        do ti = 1, contri
          do tl = 1, numj
            do tj = 1, contrj
              val = val + Integral_V_1e(     &
              Z_pot,                         &
              tcoei((tk-1)*contri+ti),       &
              coej((tl-1)*contrj+tj),        &
              tfaci(:,tk),                   &
              facj(:,tl),                    &
              expi(ti),                      &
              expj(tj),                      &
              tcodi,                         &
              tcodj,                         &
              R_pot)
            end do
          end do
        end do
      end do
    end do
  end function Calc_pppVp_1e
  
!-----------------------------------------------------------------------
!> normalization of primitive shell (GTFs) and contracted shell (GTOs)
!!
!! since all GTFs in a GTO have the same normalization coefficient, normalizing
!!
!! of GTO is equivalent to normalizing the GTO after normalizing the all GTFs.
  pure subroutine Calc_Ncoe(incbdata, incbdm)
    implicit none
    type(basis_data),intent(inout) :: incbdata(:)! input cbdata
    integer, intent(in)            :: incbdm     ! input cbdm
    integer                        :: contr    ! contr of atom, shell
    real(dp)                       :: expo(16) ! expo of |AO>
    real(dp)                       :: coe(16)  ! coe of |AO>
    integer                        :: fac(3)   ! xyz factor of |AO>
    integer                        :: L        ! angular quantum number of |AO>
    integer                        :: M        ! magnetic quantum number of |AO>
    real(dp)                       :: cod(3)   ! central coordinate of |AO>
    real(dp)                       :: i_i      ! full space integral
    integer                        :: oi,oj,ok ! loop variables for Calc_Ncoe
    ! normalization of primitive shell (GTFs)
    do oi = 1, incbdm
      contr = incbdata(oi) % contr
      expo(1:contr) = incbdata(oi) % expo(1:contr)
      coe(1:contr)  = incbdata(oi) % coe(1:contr)
      L = incbdata(oi) % L
      M = incbdata(oi) % M
      do oj = 1, contr
        incbdata(oi)%Ncoe(oj) = coe(oj) * &
        AON(expo(oj),AO_fac(1,L,M),AO_fac(2,L,M),AO_fac(3,L,M))
      end do
    end do
    ! normalization of contracted shell (GTOs)
    do oi = 1, incbdm
      contr = incbdata(oi) % contr
      expo(1:contr) = incbdata(oi) % expo(1:contr)
      coe(1:contr) = incbdata(oi)%Ncoe(1:contr)
      cod  = incbdata(oi) % pos
      L = incbdata(oi) % L
      M = incbdata(oi) % M
      i_i = 0.0_dp
      do oj = 1, contr
        do ok = 1, contr
          i_i = i_i + Integral_S_1e(         &
          coe(oj)*coe(ok),                   &
          AO_fac(:,L,M),                     &
          AO_fac(:,L,M),                     &
          expo(oj),                          &
          expo(ok),                          &
          cod,                               &
          cod)
        end do
      end do
      incbdata(oi)%Ncoe(1:contr)=incbdata(oi)%Ncoe(1:contr)/i_i**0.5_dp
    end do
  end subroutine Calc_Ncoe

!------------------------------------------------------------
!> initialising Integral_V_1e and Integral_V_2e_OS and initialising Pairwise
!!
!! Recursive Integral Strategy for Multipoles (PRISM), assign Gaussian
!!
!! production of GTOs and \hat{p} GTOs, ref 10.1002/jcc.540040206
!!
!! must be called after Calc_Ncoe and before Gaussian Integration
  subroutine Integral_init()
    implicit none
    integer  :: ii, jj, kk, ll
    real(dp) :: expTaycoe(64)
    real(dp) :: erfTaycoe(64)
    real(dp) :: ai(16), aj(16)
    real(dp) :: codi(3), codj(3)
    real(dp) :: coedx_i(32), coedy_i(32), coedz_i(32)
    real(dp) :: coedx_j(32), coedy_j(32), coedz_j(32)
    integer  :: cti, ctj
    ! initialising Integral_V_1e and Integral_V_2e_OS
    expTaycoe = 0.0_dp     ! Taylor expansion coefficients for exp(-x)
    erfTaycoe = 0.0_dp     ! Taylor expansion coefficients for I_0
    do jj = 1, 64
      expTaycoe(jj) = (-1.0_dp)**(jj-1)*(1.0_dp/factorial(jj-1))
      erfTaycoe(jj) = (-1.0_dp)**(jj-1)*(1.0_dp/(real(2*jj-1)*factorial(jj-1)))
    end do
    intTaycoe = 0.0_dp
    do ii = 0, 42
      if (ii == 0) then     ! Taycoes for I_0, GNC^1*[1 + (GNC*br2)^2 + ...]
        do jj = 1, 64
          intTaycoe(jj,1) = erfTaycoe(jj)
        end do
      else if(ii == 1) then ! Taycoes for I_1, GNC^2*[1 + (GNC*br2)^2 + ...]
        do jj = 1, 63
          intTaycoe(jj,2) = -0.5_dp*expTaycoe(jj+1)
        end do
      else                  ! Taycoes for I_n, GNC^(n+1)*[1 + (GNC*br2)^2 + ...]
        do jj = 1, 64-int(real(ii)/1.99)-mod(ii,2)
          intTaycoe(jj,ii+1) = 0.5_dp*(real(ii-1,dp)*intTaycoe(jj+1,ii-1) - &
                                       expTaycoe(jj+1))
        end do
      end if
    end do
    ! initialising PRISM
    ! data in AOpair can also be used for <AOi|\hat{p}|\hat{p}|AOj>
    allocate(AOpair(cbdm,cbdm))
    do ii = 1, cbdm
      cti = cbdata(ii) % contr
      ai = cbdata(ii) % expo
      codi = cbdata(ii) % pos
      do jj = 1, cbdm
        ctj = cbdata(jj) % contr
        aj = cbdata(jj) % expo
        codj = cbdata(jj) % pos
        !------------------------------
        AOpair(ii,jj) % contri = cti
        AOpair(ii,jj) % contrj = ctj
        allocate(AOpair(ii,jj) % sumexpo(cti, ctj),   &
                 AOpair(ii,jj) % cod(3, cti, ctj),    &
                 AOpair(ii,jj) % Gij(3, cti, ctj),    &
                 AOpair(ii,jj) % Gimij(3, cti, ctj),  &
                 source = 0.0_dp)
        do kk = 1, cti
          do ll = 1, ctj
            AOpair(ii,jj) % cod(:,kk,ll) = &
            (ai(kk)*codi(:)+aj(ll)*codj(:)) / (ai(kk)+aj(ll))
            AOpair(ii,jj) % sumexpo(kk,ll) = ai(kk)+aj(ll)
            AOpair(ii,jj) % Gij(:,kk,ll) = &
            (ai(kk)*aj(ll)/(ai(kk)+aj(ll))) * (codi(:)-codj(:))**2
            AOpair(ii,jj) % Gimij(:,kk,ll) = &
            rpi / dsqrt(AOpair(ii,jj)%sumexpo(kk,ll)) * &
            dexp(-AOpair(ii,jj)%Gij(:,kk,ll))
          end do
        end do
      end do
    end do
    ! assign facdx, coedx, facdy, coedy, facdz, coedz in cbdata
    do ii = 1, cbdm
      cti = cbdata(ii) % contr
      ai = cbdata(ii) % expo
      !---------------------------------------
      ! factor of x^m*d(exp)
      cbdata(ii)%facdx(:,1) = cbdata(ii)%fac
      cbdata(ii)%facdx(1,1) = cbdata(ii)%facdx(1,1) + 1
      cbdata(ii)%facdy(:,1) = cbdata(ii)%fac
      cbdata(ii)%facdy(2,1) = cbdata(ii)%facdy(2,1) + 1
      cbdata(ii)%facdz(:,1) = cbdata(ii)%fac
      cbdata(ii)%facdz(3,1) = cbdata(ii)%facdz(3,1) + 1
      !-------------------------------
      ! factor of d(x^m)*exp
      cbdata(ii)%facdx(:,2) = cbdata(ii)%fac
      cbdata(ii)%facdx(1,2) = max(cbdata(ii)%facdx(1,2)-1, 0)
      cbdata(ii)%facdy(:,2) = cbdata(ii)%fac
      cbdata(ii)%facdy(2,2) = max(cbdata(ii)%facdy(2,2)-1, 0)
      cbdata(ii)%facdz(:,2) = cbdata(ii)%fac
      cbdata(ii)%facdz(3,2) = max(cbdata(ii)%facdz(3,2)-1, 0)
      !---------------------------------------
      ! coefficient of x^m*d(exp)
      cbdata(ii)%coedx(1:cti) = -2.0_dp * cbdata(ii)%Ncoe(1:cti) * ai(1:cti)
      cbdata(ii)%coedy(1:cti) = -2.0_dp * cbdata(ii)%Ncoe(1:cti) * ai(1:cti)
      cbdata(ii)%coedz(1:cti) = -2.0_dp * cbdata(ii)%Ncoe(1:cti) * ai(1:cti)
      !---------------------------------
      ! coefficient of d(x^m)*exp
      cbdata(ii)%coedx(cti+1:2*cti) = cbdata(ii)%fac(1) * cbdata(ii)%Ncoe(1:cti)
      cbdata(ii)%coedy(cti+1:2*cti) = cbdata(ii)%fac(2) * cbdata(ii)%Ncoe(1:cti)
      cbdata(ii)%coedz(cti+1:2*cti) = cbdata(ii)%fac(3) * cbdata(ii)%Ncoe(1:cti)
    end do
  end subroutine Integral_init

!-----------------------------------------------------------------------
!> full space integration of product of 2 Gaussian functions in Cartesian
!!
!! it's not compatible with PRISM
  real(dp) pure function Integral_S_1e(&
  coe,faci,facj,expi,expj,codi,codj) result(val)
    implicit none
    real(dp),intent(in) :: coe           ! coefficient before Gaussian product
    integer,intent(in)  :: faci(3)       ! xyz factor before Gaussian product
    integer,intent(in)  :: facj(3)       ! xyz factor before Gaussian product
    real(dp),intent(in) :: expi          ! exponent before Gaussian product
    real(dp),intent(in) :: expj          ! exponent before Gaussian product
    real(dp),intent(in) :: codi(3)       ! coordinate before Gaussian product
    real(dp),intent(in) :: codj(3)       ! coordinate before Gaussian product
    real(dp)            :: pcoe          ! coefficient after Gaussian product
    real(dp)            :: invexpo       ! inverse of Gaussian product exponent
    real(dp)            :: cod(3)        ! coordinate of Gaussian product
    real(dp)            :: codpos(3)     ! linear transformation of cod_i
    real(dp)            :: integral(3)   ! integral of x,y,z polinomial
    real(dp)            :: itm
    real(dp)            :: mic, mic_, mic__
    integer             :: gi,gj,gk      ! loop variables for Integral_S_1e
    integral = 0.0_dp
    ! Gaussian function produntion
    invexpo = 1.0_dp / (expi+expj)
    pcoe = coe * exp(-(sum((codi(:)-codj(:))**2)) * (expi*expj)*invexpo)
    cod(:) = (codi(:)*expi+codj(:)*expj) * invexpo
    ! linear transformation of integral variables
    codpos = codi - codj
    cod = cod - codj
    itm = dsqrt(pi*invexpo)
    ! binomial expansion of x,y,z
    do gk = 1, 3
      do gi = 0, faci(gk)
        ! integral x^m*exp(-b*(x-x0)^2)
        do gj = 0, faci(gk)-gi+facj(gk)
          if (gj == 0) then
            mic = itm
          else if(gj == 1) then
            mic_ = mic
            mic = cod(gk) * mic_
          else
            mic__ = mic_
            mic_ = mic
            mic = cod(gk)*mic_ + 0.5_dp*(real(gj-1,dp)*mic__)*invexpo
          end if
        end do
        integral(gk) = integral(gk) + &
        binom(faci(gk),gi)*(-codpos(gk))**(gi)*mic
      end do
    end do
    val = pcoe * integral(1) * integral(2) * integral(3)
    return
  end function Integral_S_1e
  
!-----------------------------------------------------------------------
!> integration of electron-nuclear attraction potential in Cartesian coordinate
!!
!! scheme: Obara-Saika
!!
!! inategral transformation: (x-xi)^m*(x-xj)^n*expo(-b*x^2)*expo(x^2*t^2)dxdt
!!
!! it's not compatible with PRISM
  real(dp) pure function Integral_V_1e(&
  Z,coei,coej,faci,facj,expi,expj,codi,codj,rn) result(val)
    implicit none
    real(dp),intent(in) :: Z                  ! atomic number of potential atom
    real(dp),intent(in) :: coei               ! coeffcient of |AOi>
    real(dp),intent(in) :: coej               ! coeffcient of |AOj>
    integer,intent(in)  :: faci(3)            ! xyz factor of |AOi>
    integer,intent(in)  :: facj(3)            ! xyz factor of |AOj>
    real(dp),intent(in) :: expi               ! exponent of |AOi>
    real(dp),intent(in) :: expj               ! exponent of |AOj>
    real(dp)            :: expo               ! exponnet of product shell
    real(dp)            :: invexpo            ! inverse of expo
    real(dp)            :: coe                ! coefficient of product shell
    real(dp),intent(in) :: codi(3)            ! coordination of i
    real(dp),intent(in) :: codj(3)            ! coordination of j
    real(dp),intent(in) :: rn                 ! nuclear radius in Bohr
    real(dp)            :: invrn, GNC         ! Gauss finite nuclear correction
    real(dp)            :: cod(3)
    real(dp)            :: R2
    ! t-containing coefficients, t2pb_xyz(1) = coeff.(t^2+b)^(1/2).x integration
    !DIR$ ATTRIBUTES ALIGN:align_size :: t2pb_xyz
    real(dp)            :: t2pb_xyz(16,3)
    ! t-containing coefficients, t2pb(1) = coeff.(t^2+b)^(1/2).x,y,z integration
    !DIR$ ATTRIBUTES ALIGN:align_size :: t2pb, mic, mic_, mic__
    real(dp)            :: t2pb(16), mic(16), mic_(16), mic__(16)
    real(dp)            :: int, int_mic, int_mic_, int_mic__
    ! number of Taylor expansion series of integration at X=0
    integer             :: tayeps
    ! direct integration (X > xts); Taylor expansion integration (X <= xts)
    real(dp)            :: xts
    integer             :: vi,vj,vk,vo    ! loop variables for Integral_V_1e
    integer             :: vmic,vmic_
    integer             :: max1, max2, max3
    real(dp)            :: tmp, br2, expbr2, invbr2, prec

    select case (2*(sum(faci)+sum(facj)+4))
      case(0:5)
        tayeps = 5
        xts = 0.01
      case(6:10)
        tayeps = 12
        xts = 0.2
      case(11:15)
        tayeps = 15
        xts = 0.5
      case(16:20)
        tayeps = 20
        xts = 1.0
      case(21:25)
        tayeps = 25
        xts = 2.0
      case(26:30)
        tayeps = 30
        xts = 3.0
      case(31:40)
        tayeps = 35
        xts = 4.0
      case(41:50)
        tayeps = 45
        xts = 5.0
    end select
    !--------------------------
    ! Gaussian production
    coe = -Z * coei * coej
    expo = expi + expj
    invexpo = 1.0_dp / (expi + expj)
    coe = coe * exp(-(sum((codi(:)-codj(:))**2)) * (expi*expj)*invexpo)
    cod(:) = (codi(:)*expi+codj(:)*expj) * invexpo
    R2 = sum(cod(:)**2)
    br2 = expo*R2
    ! Gauss finite nuclear correction
    invrn = 1.0_dp / rn
    GNC = invrn / dsqrt(invrn*invrn+expo)
    expbr2 = exp(-br2*GNC*GNC)
    ! use binomial expansion for easy storage of t-containing coefficients.
    !--------------------------
    ! integral of x,y,z, generate coefficient exp(-b*((cod(1))^2*t^2)/(t^2+b))
    do vk = 1, 3
      t2pb_xyz(:,vk) = 0.0_dp
      do vi = 0, faci(vk)
        vj = 0
        do vj = 0, facj(vk)
          mic = 0.0_dp
          mic_ = 0.0_dp
          mic__ = 0.0_dp
          max1 = faci(vk)+facj(vk)-vi-vj
          do vmic = 0, max1
            if (vmic == 0) then
              mic(1) = rpi
            else if(vmic == 1) then
              mic_(1) = rpi
              mic(1) = 0.0_dp
              mic(2) = cod(vk) * expo * rpi
            else
              mic__ = mic_
              mic_ = mic
              mic = 0.0_dp
              do vmic_ = 0, vmic
                mic(vmic_+2) = mic(vmic_+2) + &
                0.5_dp * real(vmic-1,dp) * mic__(vmic_+1) + &
                cod(vk) * expo * mic_(vmic_+1)
              end do
            end if
          end do
          t2pb_xyz(:,vk) = t2pb_xyz(:,vk) + &
          binom(faci(vk),vi) * &
          binom(facj(vk),vj) * &
          (-codi(vk))**(vi) * (-codj(vk))**(vj) * mic
        end do
      end do
    end do
    !--------------------------
    ! integral of t
    ! 2*(-1)^k*b^((1-m)/2)*coe*(u-1)^k*(u+1)^k*exp(-b*R2*u^2)
    t2pb = 0.0_dp
    max1 = faci(1)+facj(1)+1
    max2 = faci(2)+facj(2)+1
    max3 = faci(3)+facj(3)+1
    do vi = 0, max1
      do vj = 0, max2
        tmp = t2pb_xyz(vi+1,1) * t2pb_xyz(vj+1,2)
        do vk = 0, max3
          t2pb(vi+vj+vk+2) = t2pb(vi+vj+vk+2) + tmp * t2pb_xyz(vk+1,3)
        end do
      end do
    end do
    val = 0.0_dp
    max1 = sum(faci) + sum(facj) + 4
    if (abs(br2) < 1E-13) then
      ! k = vi, t2pb(vi+2) = coe
      do vi = 0, max1
        int = 0.0_dp
        do vj = 0, vi
          tmp = (-1.0_dp)**(vj)*binom(vi,vj)
          do vk = 0, vi
            int_mic = GNC**(2*vi-vj-vk+1) / &
            (real(2*vi-vj-vk,dp) + 1.0_dp)
            int = int + tmp * binom(vi,vk) * int_mic
          end do
        end do
        val = val + 2.0_dp*(-1.0_dp)**(vi) * expo**(-vi-1) * t2pb(vi+2) * int
      end do
    else if (abs(br2) <= xts) then
      ! k = vi, t2pb(vi+2) = coe
      prec = GNC*GNC*br2
      do vi = 0, max1
        int = 0.0_dp
        do vj = 0, vi
          tmp = (-1.0_dp)**(vj) * binom(vi,vj)
          do vk = 0, vi
            int_mic = 0.0_dp
            if (2*vi-vj-vk == 0) then
              do vo = 1, tayeps
                int_mic = int_mic + prec**(vo-1)*intTaycoe(vo,1)
              end do
              int_mic = int_mic * GNC
            else if (2*vi-vj-vk == 1) then
              do vo = 1, tayeps
                int_mic = int_mic + prec**(vo-1)*intTaycoe(vo,2)
              end do
              int_mic = int_mic * GNC**2
            else
              do vo = 1, tayeps
                int_mic = int_mic + prec**(vo-1)*intTaycoe(vo,2*vi-vj-vk+1)
              end do
              int_mic = int_mic * GNC**(2*vi-vj-vk+1)
            end if
            int = int + tmp * binom(vi,vk) * int_mic
          end do
        end do
        val = val + 2.0_dp*(-1.0_dp)**(vi) * expo**(-vi-1) * t2pb(vi+2) * int
      end do
    else
      !k = vi, t2pb(vi+2) = coe
      invbr2 = 0.5_dp / br2
      prec = 0.5_dp * dsqrt(pi/br2) * erf(dsqrt(br2)*GNC)
      do vi = 0, max1
        int = 0.0_dp
        do vj = 0, vi
          tmp = (-1.0_dp)**(vj) * binom(vi,vj)
          do vk = 0, vi
            do vmic = 0, 2*vi-vj-vk
              if (vmic == 0) then
                int_mic = prec
              else if(vmic == 1) then
                int_mic_ = int_mic
                int_mic = (1.0_dp - expbr2) * invbr2
              else
                int_mic__ = int_mic_
                int_mic_ = int_mic
                int_mic = (real(vmic-1,dp)*int_mic__ - &
                GNC**(vmic-1)*expbr2) * invbr2
              end if
            end do
            int = int + tmp * binom(vi,vk) * int_mic
          end do
        end do
        val = val + 2.0_dp*(-1.0_dp)**(vi) * expo**(-vi-1) * t2pb(vi+2) * int
      end do
    end if
    val = val / rpi
    val = val * coe
    return
  end function Integral_V_1e
  
!-----------------------------------------------------------------------
!> integration of two-electron repulsion potential in Cartesian coordinate
!!
!! scheme: Obara-Saika
!!
!! EXPRESS: |AOi>:(x1 - xi), |AOj>:(x1 - xj), |AOk>:(x2 - xk), |AOl>:(x2 - xl)
!! xA  xB
!! Li > Lj, Lk > Ll
!! 
!!
!! The recurrence is vectorised across primitive quartets.  The contraction is
!! performed here in lexical primitive order with Neumaier compensation.
  recursive subroutine Integral_V_2e_OS_PRISM_vec(scratch, contr, coe, codi, codj, codk, codl, &
  codA, codB, A, B, G, Gim, faci, facj, fack, facl, sumint, correction)
    implicit none
    type(V2eVecScratch),intent(inout) :: scratch
    integer,intent(in)  :: contr
    real(dp),intent(in) :: codi(3), codj(3), codk(3), codl(3)
    real(dp),intent(in) :: coe(:), codA(:,:), codB(:,:)
    real(dp),intent(in) :: A(:), B(:), G(:,:), Gim(:,:)
    integer,intent(in)  :: faci(3), facj(3), fack(3), facl(3)
    real(dp),intent(inout) :: sumint, correction
    real(dp)            :: integral(contr)
    ! composite parameters, ref 10.1002/jcc.540040206
    real(dp)            :: rou(contr), D(contr,3), X(contr)
    real(dp)            :: codAi(contr,3), codBk(contr,3)
    real(dp)            :: coe1A(contr,3),coe1B(contr,3),coe2A(contr)
    real(dp)            :: coe2B(contr),coe2AB(contr),coe3A(contr),coe3B(contr)
    integer             :: facij1(3), fackl1(3), nt
    real(dp)            :: int_mic(contr), int_mic_(contr), int_mic__(contr)
    real(dp)            :: supp1(contr), supp2(contr), supp3(contr), supp4(contr)
    integer             :: tayeps  ! Taylor expansion series of integral at X=0
    ! direct integration (X > xts); Taylor expansion integration (X <= xts)
    real(dp)            :: xts
    integer             :: ri, rj, rk, ii, lane  ! loop variables for Integral_V_2e_OS
    integer             :: nij, nkl, nrec, nl, npl
    logical             :: X_zero(contr), X_taylor(contr), X_direct(contr)
    rou = A*B / (A+B)
    D(:,1) = rou(:) * (codA(:,1)-codB(:,1))**2
    D(:,2) = rou(:) * (codA(:,2)-codB(:,2))**2
    D(:,3) = rou(:) * (codA(:,3)-codB(:,3))**2
    X(:)   = sum(D(:,:),dim=2)
    ! for non-normalized inputs, do not consider Gx, Gy, and Gz.
    codAi(:,1)  = codA(:,1) - codi(1)
    codAi(:,2)  = codA(:,2) - codi(2)
    codAi(:,3)  = codA(:,3) - codi(3)
    codBk(:,1)  = codB(:,1) - codk(1)
    codBk(:,2)  = codB(:,2) - codk(2)
    codBk(:,3)  = codB(:,3) - codk(3)
    coe1A(:,1) = (A(:)*(codA(:,1)-codB(:,1))/(A(:)+B(:)))
    coe1A(:,2) = (A(:)*(codA(:,2)-codB(:,2))/(A(:)+B(:)))
    coe1A(:,3) = (A(:)*(codA(:,3)-codB(:,3))/(A(:)+B(:)))
    coe1B(:,1) = (B(:)*(codB(:,1)-codA(:,1))/(A(:)+B(:)))
    coe1B(:,2) = (B(:)*(codB(:,2)-codA(:,2))/(A(:)+B(:)))
    coe1B(:,3) = (B(:)*(codB(:,3)-codA(:,3))/(A(:)+B(:)))
    coe2A    = 1.0_dp/(2.0_dp*A)
    coe2B    = 1.0_dp/(2.0_dp*B)
    coe2AB   = 1.0_dp/(2.0_dp*(A+B))
    coe3A    = A/(2.0_dp*B*(A+B))
    coe3B    = B/(2.0_dp*A*(A+B))
    nt       = sum(faci) + sum(facj) + sum(fack) + sum(facl)
    facij1   = faci + facj + 1
    fackl1   = fack + facl + 1
    nij      = maxval(facij1)
    nkl      = maxval(fackl1)
    nrec     = maxval(facij1 + fackl1)
    nl       = maxval(facl) + 1
    npl      = nt + 6
    call Reserve_V2eVecScratch(scratch, contr, nij, nkl, nrec, nl, npl)

    associate (Gnm => scratch%Gnm(1:contr, :, 1:nij, 1:nkl, 1:nrec), &
               Itrans => scratch%Itrans(1:contr, :, 1:nl, 1:nrec), &
               I => scratch%I(1:contr, :, 1:nrec), PL => scratch%PL(1:contr, 1:npl))

    !=============================================================
    ! Obara-Saika scheme for high angular momentum Gaussian functions
    ! Ix(ni+nj,0,nk+nl,0,u)
    select case (2*(nt+6)-2)
    case(0:5)
      tayeps = 5
      xts = 0.01
    case(6:10)
      tayeps = 12
      xts = 0.2
    case(11:15)
      tayeps = 15
      xts = 0.5
    case(16:20)
      tayeps = 20
      xts = 1.0
    case(21:25)
      tayeps = 25
      xts = 2.0
    case(26:30)
      tayeps = 30
      xts = 3.0
    case(31:40)
      tayeps = 35
      xts = 4.0
    case(41:50)
      tayeps = 45
      xts = 5.0
    end select
    Gnm = 0.0_dp
    !---------------------------------------------------------------------
    ! reduce the factor (1-t^2)^(1/2)*exp(-Dx*t^2)
    do ii = 1, 3
      Gnm(:,ii,1,1,1) = Gim(:,ii)
      if (facij1(ii) >= 2) then
        Gnm(:,ii,2,1,1) = Gnm(:,ii,2,1,1) + Gnm(:,ii,1,1,1)*codAi(:,ii)
        Gnm(:,ii,2,1,2) = Gnm(:,ii,2,1,2) + Gnm(:,ii,1,1,1)*coe1B(:,ii)
      end if
      do ri = 3, facij1(ii)
        do rj = 1, ri - 1
          Gnm(:,ii,ri  ,1,rj  ) = Gnm(:,ii,ri,1,rj)     + &
          Gnm(:,ii,ri-2,1,rj  ) * real(ri-2,dp)*coe2A(:) + &
          Gnm(:,ii,ri-1,1,rj  ) * codAi(:,ii)
          Gnm(:,ii,ri  ,1,rj+1) = Gnm(:,ii,ri,1,rj+1)   - &
          Gnm(:,ii,ri-2,1,rj  ) * real(ri-2,dp)*coe3B(:) + &
          Gnm(:,ii,ri-1,1,rj  ) * coe1B(:,ii)
        end do
      end do
      if (fackl1(ii) >= 2) then
        Gnm(:,ii,1,2,1) = Gnm(:,ii,1,2,1) + Gnm(:,ii,1,1,1)*codBk(:,ii)
        Gnm(:,ii,1,2,2) = Gnm(:,ii,1,2,2) + Gnm(:,ii,1,1,1)*coe1A(:,ii)
      end if
      do ri = 3, fackl1(ii)
        do rj = 1, ri - 1
          Gnm(:,ii,1,ri  ,rj  ) = Gnm(:,ii,1,ri,rj)     + &
          Gnm(:,ii,1,ri-2,rj  ) * real(ri-2,dp)*coe2B(:) + &
          Gnm(:,ii,1,ri-1,rj  ) * codBk(:,ii)
          Gnm(:,ii,1,ri  ,rj+1) = Gnm(:,ii,1,ri,rj+1)   - &
          Gnm(:,ii,1,ri-2,rj  ) * real(ri-2,dp)*coe3A(:) + &
          Gnm(:,ii,1,ri-1,rj  ) * coe1A(:,ii)
        end do
      end do
      !---------------------------------------------------------------------
      ! use G(n+1,m) recursion only, codAi and codBk are asymmetric
      do rk = 2, fackl1(ii)
        if (facij1(ii) >= 2) then
          do rj = 1, rk
            Gnm(:,ii,2,rk  ,rj  ) = Gnm(:,ii,2,rk,rj)     + &
            Gnm(:,ii,1,rk  ,rj  ) * codAi(:,ii)
            Gnm(:,ii,2,rk  ,rj+1) = Gnm(:,ii,2,rk,rj+1)   + &
            Gnm(:,ii,1,rk  ,rj  ) * coe1B(:,ii)           + &
            Gnm(:,ii,1,rk-1,rj  ) * real(rk-1,dp)*coe2AB(:)
          end do
        end if
        do ri = 3, facij1(ii)
          do rj = 1, ri + rk - 2
            Gnm(:,ii,ri  ,rk  ,rj  ) = &
            Gnm(:,ii,ri  ,rk  ,rj  ) + &
            Gnm(:,ii,ri-2,rk  ,rj  ) * real(ri-2,dp)*coe2A(:) +&
            Gnm(:,ii,ri-1,rk  ,rj  ) * codAi(:,ii)
            Gnm(:,ii,ri  ,rk  ,rj+1) = &
            Gnm(:,ii,ri  ,rk  ,rj+1) - &
            Gnm(:,ii,ri-2,rk  ,rj  ) * real(ri-2,dp)*coe3B(:) +&
            Gnm(:,ii,ri-1,rk  ,rj  ) * coe1B(:,ii) + &
            Gnm(:,ii,ri-1,rk-1,rj  ) * real(rk-1,dp)*coe2AB(:)
          end do
        end do
      end do
    end do
    !---------------------------------------------------------------------
    ! transfer from codi to codj, codk to codl
    Itrans = 0.0_dp
    I = 0.0_dp
    do ii = 1, 3
      do rk = 1, facl(ii) + 1
        do ri = 0, facj(ii)
          do rj = 1, facij1(ii) + fackl1(ii)
            Itrans(:,ii,rk,rj) = Itrans(:,ii,rk,rj) + &
            binom(facj(ii),ri)*(codi(ii)-codj(ii))**(ri) * &
            Gnm(:,ii,facij1(ii)-ri,fackl1(ii)-(rk-1),rj)
          end do
        end do
      end do
      do ri = 0, facl(ii)
        do rj = 1, facij1(ii) + fackl1(ii)
          I(:,ii,rj) = I(:,ii,rj) + &
          binom(facl(ii),ri)*(codk(ii)-codl(ii))**(ri) * &
          Itrans(:,ii,ri+1,rj)
        end do
      end do
    end do
    !---------------------------------------------------------------------
    ! product to PL
    PL = 0.0_dp
    do ri = 0, facij1(1)+fackl1(1)-1
      do rj = 0, facij1(2)+fackl1(2)-1
        do rk = 0, facij1(3)+fackl1(3)-1
          PL(:,ri+rj+rk+1) = PL(:,ri+rj+rk+1) + &
          I(:,1,ri+1) * I(:,2,rj+1) * I(:,3,rk+1)
        end do
      end do
    end do
    do ri = 1, nt+6
      PL(:,ri) = PL(:,ri) * 2.0_dp * dsqrt(rou(:)/pi)
    end do
    !---------------------------------------------------------------------
    ! integral of t exp(-X*t^2)*PL(t^2), 0 -> 1
    X_zero = abs(X) < 1.0E-13_dp
    X_taylor = .not. X_zero .and. abs(X) <= xts
    X_direct = .not. (X_zero .or. X_taylor)
    integral = 0.0_dp
    do ri = 1, nt+6
      where (X_zero)
        integral = integral + PL(:,ri) / real(2*ri-1,dp)
      end where
    end do
    int_mic = 0.0_dp
    do rk = 1, tayeps
      where (X_taylor)
        int_mic = int_mic + X**(rk-1)*intTaycoe(rk,1)
      end where
    end do
    where (X_taylor)
      integral = PL(:,1) * int_mic
    end where
    do ri = 2, nt+6
      int_mic = 0.0_dp
      do rk = 1, tayeps
        where (X_taylor)
          int_mic = int_mic + X**(rk-1)*intTaycoe(rk,2*ri-1)
        end where
      end do
      where (X_taylor)
        integral = integral + PL(:,ri) * int_mic
      end where
    end do
    where (X_direct)
      supp1 = dsqrt(pi/X) * erf(dsqrt(X)) / 2.0_dp
      supp3 = 1.0_dp / (2.0_dp*X)
      supp4 = -exp(-X)
      supp2 = (1.0_dp+supp4) * supp3
    end where
    do ri = 1, nt+6
      where (X_direct)
        int_mic = supp1
      end where
      if (ri >= 2) then
        where (X_direct)
          int_mic_ = int_mic
          int_mic = supp2
        end where
      end if
      do rj = 2, 2*ri - 2
        where (X_direct)
          int_mic__ = int_mic_
          int_mic_ = int_mic
          int_mic = (supp4+real(rj-1,dp)*int_mic__) * supp3
        end where
      end do
      where (X_direct)
        integral = integral + PL(:,ri) * int_mic
      end where
    end do
    !DIR$ NOVECTOR
    do lane = 1, contr
      call Neumaier_Add(sumint, correction, coe(lane)*integral(lane))
    end do
    end associate
  end subroutine Integral_V_2e_OS_PRISM_vec

!-----------------------------------------------------------------------
!> integration of two-electron repulsion potential in Cartesian coordinate
!!
!! scheme: Rys quadrature
!!
!! EXPRESS: |AOi>:(x1 - xi), |AOj>:(x1 - xj), |AOk>:(x2 - xk), |AOl>:(x2 - xl)
!! xA  xB
!! Li > Lj, Lk > Ll
  real(dp) pure function Integral_V_2e_Rys(&
  faci,facj,fack,facl,ai,aj,ak,al,codi,codj,codk,codl) result(int)
    implicit none
    integer,intent(in)  :: faci(3), facj(3), fack(3), facl(3)
    real(dp),intent(in) :: ai, aj, ak, al
    real(dp),intent(in) :: codi(3), codj(3), codk(3), codl(3)
    ! composite parameters, ref 10.1002/jcc.540040206
    real(dp)            :: codA(3), codB(3), A, B, rou
    real(dp)            :: D(3), X, G(3)
    integer             :: nt             ! number of polyfactor
    integer             :: nroots         ! Rys quadrature: number of roots
    ! Rys quadrature: roots and weights ref 10.1016/0021-9991(76)90008-5
    real(dp)            :: Gim(3)
    real(dp)            :: f0, f1, f2, f3, f4
    real(dp)            :: coe1A(3),coe1B(3),coe2A,coe2B,coe2AB,coe3A,coe3B
    real(dp)            :: u(6), w(6), t2
    !DIR$ ATTRIBUTES ALIGN:align_size :: Gnm
    real(dp)            :: Gnm(9,9,20,3)
    ! ni transfer to nj, nk tranfer to nl
    !DIR$ ATTRIBUTES ALIGN:align_size :: Itrans,I,PL
    real(dp)            :: Itrans(5,20,3), I(20,3), PL(24)
    integer             :: facij1(3), fackl1(3)
    integer             :: ri, rj, rk, ii  ! loop variables for Integral_V_2e_OS
    ! Gaussian product
    codA(:) = (ai*codi(:)+aj*codj(:)) / (ai+aj)
    codB(:) = (ak*codk(:)+al*codl(:)) / (ak+al)
    A = ai + aj
    B = ak + al
    rou = A*B / (A+B)
    D(:) = rou * (codA(:)-codB(:))**2
    X = sum(D)
    ! for non-normalized inputs, do not consider Gx, Gy, and Gz.
    G(:) = (ai*aj/(ai+aj)) * (codi(:)-codj(:))**2 + &
    (ak*al/(ak+al)) * (codk(:)-codl(:))**2
    Gim(:) = pi / dsqrt(A*B) * exp(-G(:))
    coe1A(:) = (A*(codA(:)-codB(:))/(A+B))
    coe1B(:) = (B*(codB(:)-codA(:))/(A+B))
    coe2A = 1.0_dp/(2.0_dp*A)
    coe2B = 1.0_dp/(2.0_dp*B)
    coe2AB = 1.0_dp/(2.0_dp*(A+B))
    coe3A = A/(2.0_dp*B*(A+B))
    coe3B = B/(2.0_dp*A*(A+B))
    ! change the definitions of codA and codB for ease of computation
    codA = codA - codi
    codB = codB - codk
    ! Rys quadrature scheme for low angular momentum Gaussian functions
    nt = sum(faci) + sum(facj) + sum(fack) + sum(facl)
    facij1 = faci + facj + 1
    fackl1 = fack + facl + 1
    !=============================================================
    ! Rys quadrature for low angular momentum Gaussian functions
    nroots = int(real(nt)/1.99) + 1
    ! call rys_roots(from LibCInt), GRysroots(from GAMESS)
    call Grysroots(nroots, X, u, w)
    int = 0.0_dp
    do ri = 1, nroots
      !t2 = u(ri) / (rou+u(ri))             ! for rys_roots
      t2 = u(ri)                           ! for GRysroots
      f0 = coe2AB*t2
      f3 = coe2A-coe3B*t2
      f4 = coe2B-coe3A*t2
      !---------------------------Gnm---------------------------
      Gnm = 0.0_dp
      do ii = 1, 3
        f1 = codA(ii)+coe1B(ii)*t2
        f2 = codB(ii)+coe1A(ii)*t2
        Gnm(1,1,1,ii) = Gim(ii)
        do rj = 2, facij1(ii)
          if (rj == 2) then
            Gnm(2,1,1,ii) = Gnm(2,1,1,ii) + &
            Gnm(1,1,1,ii)*f1
          else
            Gnm(rj,1,1,ii) = Gnm(rj,1,1,ii) + &
            Gnm(rj-2,1,1,ii)*real(rj-2)*f3 + &
            Gnm(rj-1,1,1,ii)*f1
          end if
        end do
        do rj = 2, fackl1(ii)
          if (rj == 2) then
            Gnm(1,2,1,ii) = Gnm(1,2,1,ii) + &
            Gnm(1,1,1,ii)*f2
          else
            Gnm(1,rj,1,ii) = Gnm(1,rj,1,ii) + &
            Gnm(1,rj-2,1,ii)*real(rj-2)*f4 + &
            Gnm(1,rj-1,1,ii)*f2
          end if
        end do
        do rj = 2, fackl1(ii)
          do rk = 2, facij1(ii)
            if (rk == 2) then
              Gnm(2,rj,1,ii) = Gnm(2,rj,1,ii) + &
              Gnm(1,rj,1,ii)*f1 + &
              Gnm(1,rj-1,1,ii)*real(rj-1)*f0
            else
              Gnm(rk,rj,1,ii) = Gnm(rk,rj,1,ii) + &
              Gnm(rk-2,rj,1,ii)*real(rk-2)*f3 + &
              Gnm(rk-1,rj,1,ii)*f1 + &
              Gnm(rk-1,rj-1,1,ii)*real(rj-1)*f0
            end if
          end do
        end do
      end do
      !---------------------------I---------------------------
      Itrans = 0.0_dp
      I = 0.0_dp
      do ii = 1, 3
        do rk = 1, facl(ii)+1
          do rj = 0, facj(ii)
            Itrans(rk,1,ii) = Itrans(rk,1,ii) + &
            binom(facj(ii),rj)*(codi(ii)-codj(ii))**rj*&
            Gnm(facij1(ii)-rj,fackl1(ii)+1-rk,1,ii)
          end do
        end do
        do rj = 0, facl(ii)
          I(1,ii) = I(1,ii) + &
          binom(facl(ii),rj)*(codk(ii)-codl(ii))**(rj)*&
          Itrans(rj+1,1,ii)
        end do
      end do
      !---------------------------PL---------------------------
      PL(1) = I(1,1) * I(1,2) * I(1,3)
      int = int + w(ri)*PL(1)
    end do
    int = int * 2.0_dp * dsqrt(rou/pi)
  end function Integral_V_2e_Rys
  
end module Hamiltonian
