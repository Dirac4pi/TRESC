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
  !DIR$ ATTRIBUTES ALIGN:align_size :: pz3Vpx, py3Vpz, pz3Vpy, exSR
  real(dp),allocatable    :: px3Vpx(:,:)   ! <AOi|pppVp|AOj>
  real(dp),allocatable    :: py3Vpy(:,:)
  real(dp),allocatable    :: pz3Vpz(:,:)
  real(dp),allocatable    :: px3Vpy(:,:)
  real(dp),allocatable    :: py3Vpx(:,:)
  real(dp),allocatable    :: px3Vpz(:,:)
  real(dp),allocatable    :: pz3Vpx(:,:)
  real(dp),allocatable    :: py3Vpz(:,:)
  real(dp),allocatable    :: pz3Vpy(:,:)
  complex(dp),allocatable :: exSR(:,:)     ! extended SR matrix
  ! <AOi|pxVpy3|AOj> = Trans(<AOi|py3Vpx|AOj>)

! interation Taylor expansion coefficients for Integral_V_1e and Integral_V_2e
  !DIR$ ATTRIBUTES ALIGN:align_size :: intTaycoe
  real(dp)                :: intTaycoe(64,43)

  contains

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
    write(60,'(A,I3,A,I3)') &
    '  threads using:',threads,'; node nproc:',nproc
    if (nproc <= threads) then
      write(*,'(A,I3)') &
      'tkernel: Calculation will perform serially! node nproc: ',nproc
      write(60,'(A)') '  Warning: calculation will be performed SERIALLY!'
    else
      call getenv('KMP_AFFINITY',ch30)
      if (ch30 == '') then
        write(60,'(A)') '  KMP_AFFINITY = None'
      else
        write(60,'(A)') '  KMP_AFFINITY = '//trim(ch30)
      end if
    end if
    write(60,'(A)') '  Integral_V initialize'
    call Integral_V_init()
    write(60,'(A)') '  complete! stored in: intTaycoe'
    write(60,'(A)') '  ----------<HAMILTONIAN>----------'
    write(60,'(A)') '  spinor basis causes additional cost in scalar SCF.'
    !-----------------------------------------------
    ! 1e integral calculation
    write(60,'(A)') '  normalization of GTOs'
    call Calc_Ncoe(cbdata, cbdm)
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



    ! allocate(i_T_j_read(sbdm,sbdm))
    ! open(14, file='ao_overlap.mat', status="old", action="read")
    ! do si = 1, sbdm
    !   read(14, *) i_T_j_read(si,:)
    ! end do
    ! do si = 1, sbdm
    !   do sj = 1, sbdm
    !     if (abs(i_j_s(si,sj)-i_T_j_read(si,sj)) > 1.0d-6) write(60,'(I6,A,I1,I6,E12.5,E12.5,I5,I5)') si,' L=',sbdata(si)%L, sj, i_j_s(si,sj), abs(i_j_s(si,sj)-i_T_j_read(si,sj)), sbdata(sj)%L, sbdata(sj)%M
    !   end do
    ! end do
    ! close(14)
    ! write(60,*) 'complete!!!!!!'
    ! stop






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
    write(60,'(A,I3,A,I3)') &
    '  threads using:',threads,'; node nproc:',nproc
    if (nproc <= threads) then
      write(*,'(A,I3)') &
      'tkernel: Calculation will perform serially! node nproc: ',nproc
      write(60,'(A)') '  Warning: calculation will be performed SERIALLY!'
    else
      call getenv('KMP_AFFINITY',ch30)
      if (ch30 == '') then
        write(60,'(A)') '  KMP_AFFINITY = None'
      else
        write(60,'(A)') '  KMP_AFFINITY = '//trim(ch30)
      end if
    end if
    write(60,'(A)') '  Integral_V initialize'
    call Integral_V_init()
    write(60,'(A)') '  complete! stored in: intTaycoe'
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
    write(60,'(A)') '  one-electron integral calculation'
    call Assign_matrices_1e()
    write(60,'(A)') '  complete! stored in:'
    write(60,'(A)') '  i_j, i_p2_j, i_V_j, i_pVp_j (9 matrices)'
    if (pppVp) write(60,'(A)') '  i_pppVp_j (9 matrices)'
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
    if (pppVp) then
      call cfgo(px3Vpx)
      call cfgo(py3Vpy)
      call cfgo(pz3Vpz)
      call cfgo(px3Vpy)
      call cfgo(py3Vpx)
      call cfgo(px3Vpz)
      call cfgo(pz3Vpx)
      call cfgo(py3Vpz)
      call cfgo(pz3Vpy)
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
    integer          :: Li                   ! angular quantum number of |AOi>
    integer          :: Mi                   ! magnetic quantum number of |AOi>
    integer          :: contrj               ! contr of atoMj, shelLj
    integer          :: Lj                   ! angular quantum number of |AOj>
    integer          :: Mj                   ! magnetic quantum number of |AOj>
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
      real(dp), allocatable :: i_j(:,:), i_p2_j(:,:), i_V_j(:,:), pxVpx(:,:),&
      pyVpy(:,:), pzVpz(:,:), pxVpy(:,:), pyVpx(:,:), pyVpz(:,:), pzVpy(:,:),&
      pxVpz(:,:), pzVpx(:,:), px3Vpx(:,:), py3Vpy(:,:), pz3Vpz(:,:),&
      px3Vpy(:,:), py3Vpx(:,:), px3Vpz(:,:), pz3Vpx(:,:),py3Vpz(:,:),pz3Vpy(:,:)
    end type
    type(threadlocal) :: tl
    allocate(i_j(cbdm,cbdm),i_V_j(cbdm,cbdm),i_p2_j(cbdm,cbdm), source=0.0_dp)
    if(pVp1e) then
      allocate(pxVpx(cbdm,cbdm),pyVpy(cbdm,cbdm),pzVpz(cbdm,cbdm),source=0.0_dp)
      allocate(pxVpy(cbdm,cbdm),pyVpx(cbdm,cbdm),pyVpz(cbdm,cbdm),source=0.0_dp)
      allocate(pzVpy(cbdm,cbdm),pxVpz(cbdm,cbdm),pzVpx(cbdm,cbdm),source=0.0_dp)
      if(pppVp) then
        allocate (px3Vpx(cbdm,cbdm),py3Vpy(cbdm,cbdm),source=0.0_dp)
        allocate (pz3Vpz(cbdm,cbdm),px3Vpy(cbdm,cbdm),source=0.0_dp)
        allocate (py3Vpx(cbdm,cbdm),px3Vpz(cbdm,cbdm),source=0.0_dp)
        allocate (pz3Vpx(cbdm,cbdm),py3Vpz(cbdm,cbdm),source=0.0_dp)
        allocate (pz3Vpy(cbdm,cbdm),source=0.0_dp)
      end if
    end if
    ! parallel zone, running results consistent with serial
    !$omp parallel num_threads(threads) default(shared) private(i,j,si,sj,sk,&
    !$omp& sl,contri,Li,Mi,contrj,Lj,Mj,expi,expj,coei,coej,codi,codj,&
    !$omp& facdx_i,facdy_i,facdz_i,coedx_i,coedy_i,coedz_i,facdx_j,&
    !$omp& facdy_j,facdz_j,coedx_j,coedy_j,coedz_j,tl) if(threads < nproc)
    allocate(&
    tl%i_j(cbdm,cbdm),tl%i_V_j(cbdm,cbdm),tl%i_p2_j(cbdm,cbdm), source=0.0_dp)
    if (pVp1e) then
      allocate(tl%pxVpx(cbdm,cbdm),tl%pyVpy(cbdm,cbdm),source=0.0_dp)
      allocate(tl%pzVpz(cbdm,cbdm),source=0.0_dp)
      allocate(tl%pxVpy(cbdm,cbdm),tl%pyVpx(cbdm,cbdm),source=0.0_dp)
      allocate(tl%pyVpz(cbdm,cbdm),tl%pzVpy(cbdm,cbdm),source=0.0_dp)
      allocate(tl%pxVpz(cbdm,cbdm),tl%pzVpx(cbdm,cbdm),source=0.0_dp)
      if (pppVp) then
        allocate (tl%px3Vpx(cbdm,cbdm),tl%py3Vpy(cbdm,cbdm),source=0.0_dp)
        allocate (tl%pz3Vpz(cbdm,cbdm),tl%px3Vpy(cbdm,cbdm),source=0.0_dp)
        allocate (tl%py3Vpx(cbdm,cbdm),tl%px3Vpz(cbdm,cbdm),source=0.0_dp)
        allocate (tl%pz3Vpx(cbdm,cbdm),tl%py3Vpz(cbdm,cbdm),source=0.0_dp)
        allocate (tl%pz3Vpy(cbdm,cbdm),source=0.0_dp)
      end if
    end if
    !$omp do schedule(dynamic, 5) collapse(2)
    do i = 1, cbdm
      do j = 1, cbdm
        si = i
        contri = cbdata(si) % contr
        Li     = cbdata(si) % L
        Mi     = cbdata(si) % M
        expi(1:contri) = cbdata(si) % expo(1:contri)
        coei(1:contri) = cbdata(si) % Ncoe(1:contri)
        codi = cbdata(si) % pos
        !---------------------------------------
        ! factor of x^m*d(exp)
        facdx_i(:,1) = AO_fac(:,Li,Mi)
        facdx_i(1,1) = facdx_i(1,1) + 1
        facdy_i(:,1) = AO_fac(:,Li,Mi)
        facdy_i(2,1) = facdy_i(2,1) + 1
        facdz_i(:,1) = AO_fac(:,Li,Mi)
        facdz_i(3,1) = facdz_i(3,1) + 1
        !-------------------------------
        ! factor of d(x^m)*exp
        facdx_i(:,2) = AO_fac(:,Li,Mi)
        facdx_i(1,2) = max(facdx_i(1,2)-1,0)
        facdy_i(:,2) = AO_fac(:,Li,Mi)
        facdy_i(2,2) = max(facdy_i(2,2)-1,0)
        facdz_i(:,2) = AO_fac(:,Li,Mi)
        facdz_i(3,2) = max(facdz_i(3,2)-1,0)
        !---------------------------------------
        ! coefficient of x^m*d(exp)
        coedx_i(1:contri) = -2.0_dp * coei(1:contri) * expi(1:contri)
        coedy_i(1:contri) = -2.0_dp * coei(1:contri) * expi(1:contri)
        coedz_i(1:contri) = -2.0_dp * coei(1:contri) * expi(1:contri)
        !---------------------------------
        ! coefficient of d(x^m)*exp
        coedx_i(contri+1:2*contri) = AO_fac(1,Li,Mi) * coei(1:contri)
        coedy_i(contri+1:2*contri) = AO_fac(2,Li,Mi) * coei(1:contri)
        coedz_i(contri+1:2*contri) = AO_fac(3,Li,Mi) * coei(1:contri)
        sj = j
        contrj = cbdata(sj) % contr
        Lj     = cbdata(sj) % L
        Mj     = cbdata(sj) % M
        expj(1:contrj) = cbdata(sj) % expo(1:contrj)
        coej(1:contrj) = cbdata(sj) % Ncoe(1:contrj)
        codj = cbdata(sj) % pos
        !---------------------------------------
        ! factor of x^m*d(exp)
        facdx_j(:,1) = AO_fac(:,Lj,Mj)
        facdx_j(1,1) = facdx_j(1,1) + 1
        facdy_j(:,1) = AO_fac(:,Lj,Mj)
        facdy_j(2,1) = facdy_j(2,1) + 1
        facdz_j(:,1) = AO_fac(:,Lj,Mj)
        facdz_j(3,1) = facdz_j(3,1) + 1
        !-------------------------------
        ! factor of d(x^m)*exp
        facdx_j(:,2) = AO_fac(:,Lj,Mj)
        facdx_j(1,2) = max(facdx_j(1,2)-1,0)
        facdy_j(:,2) = AO_fac(:,Lj,Mj)
        facdy_j(2,2) = max(facdy_j(2,2)-1,0)
        facdz_j(:,2) = AO_fac(:,Lj,Mj)
        facdz_j(3,2) = max(facdz_j(3,2)-1,0)
        !---------------------------------------
        ! coefficient of x^m*d(exp)
        coedx_j(1:contrj) = -2.0_dp * coej(1:contrj) * expj(1:contrj)
        coedy_j(1:contrj) = -2.0_dp * coej(1:contrj) * expj(1:contrj)
        coedz_j(1:contrj) = -2.0_dp * coej(1:contrj) * expj(1:contrj)
        !---------------------------------
        ! coefficient of d(x^m)*exp
        coedx_j(contrj+1:2*contrj) = AO_fac(1,Lj,Mj) * coej(1:contrj)
        coedy_j(contrj+1:2*contrj) = AO_fac(2,Lj,Mj) * coej(1:contrj)
        coedz_j(contrj+1:2*contrj) = AO_fac(3,Lj,Mj) * coej(1:contrj)
        ! calc <AOi|V|AOj>
        tl%i_V_j(si,sj) = Calc_V_1e(&
        contri,contrj,coei(1:contri),coej(1:contrj),&
        AO_fac(:,Li,Mi),AO_fac(:,Lj,Mj),expi(1:contri),expj(1:contrj),codi,codj)
        
        if (pVp1e) then
          ! calc <AOi|pVp|AOj>,totally 9 matrices,6 matrices will be calculated
          tl%pxVpx(si,sj)=Calc_pVp_1e(AO_fac(1,Li,Mi),AO_fac(1,Lj,Mj),&
          contri,contrj,coedx_i(1:2*contri),coedx_j(1:2*contrj),&
          facdx_i,facdx_j,expi(1:contri),expj(1:contrj),codi,codj)
          tl%pyVpy(si,sj)=Calc_pVp_1e(AO_fac(2,Li,Mi),AO_fac(2,Lj,Mj),&
          contri,contrj,coedy_i(1:2*contri),coedy_j(1:2*contrj),&
          facdy_i,facdy_j,expi(1:contri),expj(1:contrj),codi,codj)
          tl%pzVpz(si,sj)=Calc_pVp_1e(AO_fac(3,Li,Mi),AO_fac(3,Lj,Mj),&
          contri,contrj,coedz_i(1:2*contri),coedz_j(1:2*contrj),&
          facdz_i,facdz_j,expi(1:contri),expj(1:contrj),codi,codj)
          tl%pxVpy(si,sj)=Calc_pVp_1e(AO_fac(1,Li,Mi),AO_fac(2,Lj,Mj),&
          contri,contrj,coedx_i(1:2*contri),coedy_j(1:2*contrj),&
          facdx_i,facdy_j,expi(1:contri),expj(1:contrj),codi,codj)
          tl%pyVpx(sj,si) = tl%pxVpy(si,sj)
          tl%pxVpz(si,sj)=Calc_pVp_1e(AO_fac(1,Li,Mi),AO_fac(3,Lj,Mj),&
          contri,contrj,coedx_i(1:2*contri),coedz_j(1:2*contrj),&
          facdx_i,facdz_j,expi(1:contri),expj(1:contrj),codi,codj)
          tl%pzVpx(sj,si) = tl%pxVpz(si,sj)
          tl%pyVpz(si,sj)=Calc_pVp_1e(AO_fac(2,Li,Mi),AO_fac(3,Lj,Mj),&
          contri,contrj,coedy_i(1:2*contri),coedz_j(1:2*contrj),&
          facdy_i,facdz_j,expi(1:contri),expj(1:contrj),codi,codj)
          tl%pzVpy(sj,si) = tl%pyVpz(si,sj)
          ! calc <AOi|p^3Vp|AOj> and <AOi|pVp^3|AOj>, totally 18 matrices,
          ! 9 matrices will be calculated
          if (pppVp) then
            tl%px3Vpx(si,sj) = Calc_pppVp_1e(1,AO_fac(1,Li,Mi),AO_fac(1,Lj,Mj),&
            contri,contrj,coedx_i(1:2*contri),coedx_j(1:2*contrj),&
            facdx_i,facdx_j,expi(1:contri),expj(1:contrj),codi,codj)
            tl%py3Vpy(si,sj) = Calc_pppVp_1e(2,AO_fac(2,Li,Mi),AO_fac(2,Lj,Mj),&
            contri,contrj,coedy_i(1:2*contri),coedy_j(1:2*contrj),&
            facdy_i,facdy_j,expi(1:contri),expj(1:contrj),codi,codj)
            tl%pz3Vpz(si,sj) = Calc_pppVp_1e(3,AO_fac(3,Li,Mi),AO_fac(3,Lj,Mj),&
            contri,contrj,coedz_i(1:2*contri),coedz_j(1:2*contrj),&
            facdz_i,facdz_j,expi(1:contri),expj(1:contrj),codi,codj)
            tl%px3Vpy(si,sj) = Calc_pppVp_1e(1,AO_fac(1,Li,Mi),AO_fac(2,Lj,Mj),&
            contri,contrj,coedx_i(1:2*contri),coedy_j(1:2*contrj),&
            facdx_i,facdy_j,expi(1:contri),expj(1:contrj),codi,codj)
            tl%py3Vpx(si,sj) = Calc_pppVp_1e(2,AO_fac(2,Li,Mi),AO_fac(1,Lj,Mj),&
            contri,contrj,coedy_i(1:2*contri),coedx_j(1:2*contrj),&
            facdy_i,facdx_j,expi(1:contri),expj(1:contrj),codi,codj)
            tl%px3Vpz(si,sj) = Calc_pppVp_1e(1,AO_fac(1,Li,Mi),AO_fac(3,Lj,Mj),&
            contri,contrj,coedx_i(1:2*contri),coedz_j(1:2*contrj),&
            facdx_i,facdz_j,expi(1:contri),expj(1:contrj),codi,codj)
            tl%pz3Vpx(si,sj) = Calc_pppVp_1e(3,AO_fac(3,Li,Mi),AO_fac(1,Lj,Mj),&
            contri,contrj,coedz_i(1:2*contri),coedx_j(1:2*contrj),&
            facdz_i,facdx_j,expi(1:contri),expj(1:contrj),codi,codj)
            tl%py3Vpz(si,sj) = Calc_pppVp_1e(2,AO_fac(2,Li,Mi),AO_fac(3,Lj,Mj),&
            contri,contrj,coedy_i(1:2*contri),coedz_j(1:2*contrj),&
            facdy_i,facdz_j,expi(1:contri),expj(1:contrj),codi,codj)
            tl%pz3Vpy(si,sj) = Calc_pppVp_1e(3,AO_fac(3,Li,Mi),AO_fac(2,Lj,Mj),&
            contri,contrj,coedz_i(1:2*contri),coedy_j(1:2*contrj),&
            facdz_i,facdy_j,expi(1:contri),expj(1:contrj),codi,codj)
          end if
        end if
        !-------------<basis overlap integral>-------------
        do sk = 1, contri
          do sl = 1, contrj
            tl%i_j(si,sj) = tl%i_j(si,sj) + &
            Integral_S_1e(coei(sk)*coej(sl),&
            AO_fac(:,Li,Mi),AO_fac(:,Lj,Mj),expi(sk),expj(sl),codi,codj)
          end do
        end do
        !-------------<P2 integral>-------------
        tl%i_p2_j(si,sj) = tl%i_p2_j(si,sj) + &
        Calc_pp_1e(AO_fac(1,Li,Mi),AO_fac(1,Lj,Mj),contri,contrj,&
        coedx_i(1:2*contri),coedx_j(1:2*contrj),facdx_i,facdx_j,&
        expi(1:contri),expj(1:contrj),codi,codj)
        tl%i_p2_j(si,sj) = tl%i_p2_j(si,sj) + &
        Calc_pp_1e(AO_fac(2,Li,Mi),AO_fac(2,Lj,Mj),contri,contrj,&
        coedy_i(1:2*contri),coedy_j(1:2*contrj),facdy_i,facdy_j,&
        expi(1:contri),expj(1:contrj),codi,codj)
        tl%i_p2_j(si,sj) = tl%i_p2_j(si,sj) + &
        Calc_pp_1e(AO_fac(3,Li,Mi),AO_fac(3,Lj,Mj),contri,contrj,&
        coedz_i(1:2*contri),coedz_j(1:2*contrj),facdz_i,facdz_j,&
        expi(1:contri),expj(1:contrj),codi,codj)
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
      ! consider the relation p = -iD,
      ! p3Vp and pVp3 change sign, p2 and pVp nochange
      if (pppVp) then
        px3Vpx = px3Vpx - tl%px3Vpx
        py3Vpy = py3Vpy - tl%py3Vpy
        pz3Vpz = pz3Vpz - tl%pz3Vpz
        px3Vpy = px3Vpy - tl%px3Vpy
        py3Vpx = py3Vpx - tl%py3Vpx
        px3Vpz = px3Vpz - tl%px3Vpz
        pz3Vpx = pz3Vpx - tl%pz3Vpx
        py3Vpz = py3Vpz - tl%py3Vpz
        pz3Vpy = pz3Vpy - tl%pz3Vpy
      end if
    end if
    !$omp end critical
    ! free memory for parallel threads explicitly
    deallocate(tl%i_j, tl%i_p2_j, tl%i_V_j)
    if (pVp1e) then
      deallocate(tl%pxVpx, tl%pyVpy, tl%pzVpz)
      deallocate(tl%pxVpy, tl%pyVpx, tl%pxVpz)
      deallocate(tl%pzVpx, tl%pyVpz, tl%pzVpy)
      if (pppVp) then
        deallocate(tl%px3Vpx, tl%py3Vpy, tl%pz3Vpz)
        deallocate(tl%px3Vpy, tl%py3Vpx, tl%px3Vpz)
        deallocate(tl%pz3Vpx, tl%py3Vpz, tl%pz3Vpy)
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
    return
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
    return
  end function Calc_V_1e

!-----------------------------------------------------------------------
!> calculate (AOiAOj|V|AOkAOl)
  real(dp) pure function Calc_V_2e(&
  contri,contrj,contrk,contrl,coei,coej,coek,coel,faci,facj,fack,facl,&
  expi,expj,expk,expl,codi,codj,codk,codl) result(val)
    implicit none
    integer,intent(in)  :: contri          ! contr of |AOi>
    integer,intent(in)  :: contrj          ! contr of |AOj>
    integer,intent(in)  :: contrk          ! contr of |AOk>
    integer,intent(in)  :: contrl          ! contr of |AOl>
    real(dp),intent(in) :: coei(contri)    ! coefficients of |AOi>
    real(dp),intent(in) :: coej(contrj)    ! coefficients of |AOj>
    real(dp),intent(in) :: coek(contrk)    ! coefficients of |AOk>
    real(dp),intent(in) :: coel(contrl)    ! coefficients of |AOl>
    integer,intent(in)  :: faci(3)         ! xyz factor of |AOi>
    integer,intent(in)  :: facj(3)         ! xyz factor of |AOj>
    integer,intent(in)  :: fack(3)         ! xyz factor of |AOk>
    integer,intent(in)  :: facl(3)         ! xyz factor of |AOl>
    real(dp),intent(in) :: expi(contri)    ! expo of |AOi>
    real(dp),intent(in) :: expj(contrj)    ! expo of |AOj>
    real(dp),intent(in) :: expk(contrk)    ! expo of |AOk>
    real(dp),intent(in) :: expl(contrl)    ! expo of |AOl>
    real(dp),intent(in) :: codi(3)         ! centrol coordintates.|AOi>
    real(dp),intent(in) :: codj(3)         ! centrol coordintates.|AOj>
    real(dp),intent(in) :: codk(3)         ! centrol coordintates.|AOk>
    real(dp),intent(in) :: codl(3)         ! centrol coordintates.|AOl>
    integer             :: um, un, uo, up  ! loop variables for Calc_V_2e
    val = 0.0_dp
    do um = 1, contri
      do un = 1, contrj
        do uo = 1, contrk
          do up = 1, contrl
            val = val +                              &
            coei(um)*coej(un)*coek(uo)*coel(up)*     &
            Integral_V_2e(                           &
            faci,                                    &
            facj,                                    &
            fack,                                    &
            facl,                                    &
            expi(um),                                &
            expj(un),                                &
            expk(uo),                                &
            expl(up),                                &
            codi,                                    &
            codj,                                    &
            codk,                                    &
            codl)
          end do
        end do
      end do
    end do
    return
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
    return
  end function Calc_pVp_1e

!-----------------------------------------------------------------------
!> calculate (AOiAOj|pVp|AOkAOl) = (pAOipAOj|V|AOkAOl)
  real(dp) pure function Calc_pVp_2eij(&
  ni,nj,contri,contrj,contrk,contrl,&
  coei,coej,coek,coel,&
  faci,facj,fack,facl,&
  expi,expj,expk,expl,&
  codi,codj,codk,codl) result(val)
    implicit none
    integer,intent(in)    :: ni            ! number of x/y/z in |AOi>
    integer,intent(in)    :: nj            ! number of x/y/z in |AOj>
    integer,intent(in)    :: contri        ! contr of |AOi>
    integer,intent(in)    :: contrj        ! contr of |AOj>
    integer,intent(in)    :: contrk        ! contr of |AOk>
    integer,intent(in)    :: contrl        ! contr of |AOl>
    ! coefficients of first-order derivative of |AOi>
    real(dp),intent(in)   :: coei(2*contri)
    ! coefficients of first-order derivative of |AOj>
    real(dp),intent(in)   :: coej(2*contrj)
    real(dp),intent(in)   :: coek(contrk)  ! coefficients of |AOk>
    real(dp),intent(in)   :: coel(contrl)  ! coefficients of |AOl>
    integer,intent(in)    :: faci(3,2)     ! xyz factor of |AOi>
    integer,intent(in)    :: facj(3,2)     ! xyz factor of |AOj>
    integer,intent(in)    :: fack(3)       ! xyz factor of |AOk>
    integer,intent(in)    :: facl(3)       ! xyz factor of |AOl>
    real(dp),intent(in)   :: expi(contri)  ! expo of |AOi>
    real(dp),intent(in)   :: expj(contrj)  ! expo of |AOj>
    real(dp),intent(in)   :: expk(contrk)  ! expo of |AOk>
    real(dp),intent(in)   :: expl(contrl)  ! expo of |AOl>
    real(dp),intent(in)   :: codi(3)       ! centrol coordintates.|AOi>
    real(dp),intent(in)   :: codj(3)       ! centrol coordintates.|AOj>
    real(dp),intent(in)   :: codk(3)       ! centrol coordintates.|AOk>
    real(dp),intent(in)   :: codl(3)       ! centrol coordintates.|AOl>
    integer               :: numi, numj
    integer               :: um, un, uo, up! loop variables for Calc_pVp_2eij
    integer               :: ii, jj        ! loop variables for Calc_pVp_2eij

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
    do ii = 1, numi
    do um = 1, contri
      do jj = 1, numj
      do un = 1, contrj
        do uo = 1, contrk
          do up = 1, contrl
            val = val +                                         &
            coei((ii-1)*contri+um)*coej((jj-1)*contrj+un)*      &
            coek(uo)*coel(up)*                                  &
            Integral_V_2e(                                      &
            faci(:,ii),                                         &
            facj(:,jj),                                         &
            fack,                                               &
            facl,                                               &
            expi(um),                                           &
            expj(un),                                           &
            expk(uo),                                           &
            expl(up),                                           &
            codi,                                               &
            codj,                                               &
            codk,                                               &
            codl)
          end do
        end do
      end do
      end do
    end do
    end do
    return
  end function Calc_pVp_2eij

!-----------------------------------------------------------------------
!> calculate (AOiAOj|pVp|AOkAOl) = (pAOiAOj|V|pAOkAOl)
  real(dp) pure function Calc_pVp_2eik(&
  ni,nk,contri,contrj,contrk,contrl,&
  coei,coej,coek,coel,&
  faci,facj,fack,facl,&
  expi,expj,expk,expl,&
  codi,codj,codk,codl) result(val)
    implicit none
    integer,intent(in)    :: ni            ! number of x/y/z in |AOi>
    integer,intent(in)    :: nk            ! number of x/y/z in |AOk>
    integer,intent(in)    :: contri        ! contr of |AOi>
    integer,intent(in)    :: contrj        ! contr of |AOj>
    integer,intent(in)    :: contrk        ! contr of |AOk>
    integer,intent(in)    :: contrl        ! contr of |AOl>
    ! coefficients of first-order derivative of |AOi>
    real(dp),intent(in)   :: coei(2*contri)
    real(dp),intent(in)   :: coej(contrj)  ! coefficients of |AOj>
    ! coefficients of first-order derivative of |AOk>
    real(dp),intent(in)   :: coek(2*contrk)
    real(dp),intent(in)   :: coel(contrl)  ! coefficients of |AOl>
    integer,intent(in)    :: faci(3,2)     ! xyz factor of |AOi>
    integer,intent(in)    :: facj(3)       ! xyz factor of |AOj>
    integer,intent(in)    :: fack(3,2)     ! xyz factor of |AOk>
    integer,intent(in)    :: facl(3)       ! xyz factor of |AOl>
    real(dp),intent(in)   :: expi(contri)  ! expo of |AOi>
    real(dp),intent(in)   :: expj(contrj)  ! expo of |AOj>
    real(dp),intent(in)   :: expk(contrk)  ! expo of |AOk>
    real(dp),intent(in)   :: expl(contrl)  ! expo of |AOl>
    real(dp),intent(in)   :: codi(3)       ! centrol coordintates.|AOi>
    real(dp),intent(in)   :: codj(3)       ! centrol coordintates.|AOj>
    real(dp),intent(in)   :: codk(3)       ! centrol coordintates.|AOk>
    real(dp),intent(in)   :: codl(3)       ! centrol coordintates.|AOl>
    integer               :: numi, numk
    integer               :: um, un, uo, up! loop variables for Calc_pVp_2eik
    integer               :: ii, kk        ! loop variables for Calc_pVp_2eik

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
    do ii = 1, numi
    do um = 1, contri
      do un = 1, contrj
        do kk = 1, numk
        do uo = 1, contrk
          do up = 1, contrl
            val = val +                                         &
            coei((ii-1)*contri+um)*coej(un)*                    &
            coek((kk-1)*contrk+uo)*coel(up)*                    &
            Integral_V_2e(                                      &
            faci(:,ii),                                         &
            facj,                                               &
            fack(:,kk),                                         &
            facl,                                               &
            expi(um),                                           &
            expj(un),                                           &
            expk(uo),                                           &
            expl(up),                                           &
            codi,                                               &
            codj,                                               &
            codk,                                               &
            codl)
          end do
        end do
        end do
      end do
    end do
    end do
    return
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
    return
  end function Calc_pppVp_1e
  
!-----------------------------------------------------------------------
!> normalization of primitive shell (GTFs) and contracted shell (GTOs)
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

!-----------------------------------------------------------------------
!> full space integration of product of 2 Gaussian functions in Cartesian
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
        binomialcoe(faci(gk),gi)*(-codpos(gk))**(gi)*mic
      end do
    end do
    val = pcoe * integral(1) * integral(2) * integral(3)
    return
  end function Integral_S_1e
  
!-----------------------------------------------------------------------
!> integration of electron-nuclear attraction potential in Cartesian coordinate:
!!
!! inategral transformation: (x-xi)^m*(x-xj)^n*expo(-b*x^2)*expo(x^2*t^2)dxdt
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
    ! Gaussian shell production
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
          binomialcoe(faci(vk),vi) * &
          binomialcoe(facj(vk),vj) * &
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
          tmp = (-1.0_dp)**(vj)*binomialcoe(vi,vj)
          do vk = 0, vi
            int_mic = GNC**(2*vi-vj-vk+1) / &
            (real(2*vi-vj-vk,dp) + 1.0_dp)
            int = int + tmp * binomialcoe(vi,vk) * int_mic
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
          tmp = (-1.0_dp)**(vj) * binomialcoe(vi,vj)
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
            int = int + tmp * binomialcoe(vi,vk) * int_mic
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
          tmp = (-1.0_dp)**(vj) * binomialcoe(vi,vj)
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
            int = int + tmp * binomialcoe(vi,vk) * int_mic
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
!! EXPRESS: |AOi>:(x1 - xi), |AOj>:(x1 - xj), |AOk>:(x2 - xk), |AOl>:(x2 - xl)
!! xA  xB
!! Li > Lj, Lk > Ll
  real(dp) pure function Integral_V_2e(&
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
    real(dp)            :: int_mic, int_mic_, int_mic__
    real(dp)            :: supp1, supp2, supp3, supp4
    integer             :: facij1(3), fackl1(3)
    integer             :: tayeps  ! Taylor expansion series of integral at X=0
    ! direct integration (X > xts); Taylor expansion integration (X <= xts)
    real(dp)            :: xts
    integer             :: ri, rj, rk, ii  ! loop variables for Integral_V_2e
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
    if (.false.) then       !int(real(nt)/1.99) + 1 <= 5
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
              binomialcoe(facj(ii),rj)*(codi(ii)-codj(ii))**rj*&
              Gnm(facij1(ii)-rj,fackl1(ii)+1-rk,1,ii)
            end do
          end do
          do rj = 0, facl(ii)
            I(1,ii) = I(1,ii) + &
            binomialcoe(facl(ii),rj)*(codk(ii)-codl(ii))**(rj)*&
            Itrans(rj+1,1,ii)
          end do
        end do
        !---------------------------PL---------------------------
        PL(1) = I(1,1) * I(1,2) * I(1,3)
        int = int + w(ri)*PL(1)
      end do
      int = int * 2.0_dp * dsqrt(rou/pi)
    !=============================================================
    ! direct integral for high angular momentum Gaussian functions
    else
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
        Gnm(1,1,1,ii) = Gim(ii)
        do ri = 2, facij1(ii)
          if (ri == 2) then
            Gnm(2,1,1,ii) = Gnm(2,1,1,ii) + &
            Gnm(1,1,1,ii)*codA(ii)
            Gnm(2,1,2,ii) = Gnm(2,1,2,ii) + &
            Gnm(1,1,1,ii)*coe1B(ii)
          else
            do rj = 1, ri - 1
              Gnm(ri  ,1,rj  ,ii) = Gnm(ri,1,rj,ii) + &
              Gnm(ri-2,1,rj  ,ii) * real(ri-2)*coe2A + &
              Gnm(ri-1,1,rj  ,ii) * codA(ii)
              Gnm(ri  ,1,rj+1,ii) = Gnm(ri,1,rj+1,ii) - &
              Gnm(ri-2,1,rj  ,ii) * real(ri-2)*coe3B + &
              Gnm(ri-1,1,rj  ,ii) * coe1B(ii)
            end do
          end if
        end do
        do ri = 2, fackl1(ii)
          if (ri == 2) then
            Gnm(1,2,1,ii) = Gnm(1,2,1,ii) + &
            Gnm(1,1,1,ii)*codB(ii)
            Gnm(1,2,2,ii) = Gnm(1,2,2,ii) + &
            Gnm(1,1,1,ii)*coe1A(ii)
          else
            do rj = 1, ri - 1
              Gnm(1,ri  ,rj  ,ii) = Gnm(1,ri,rj,ii) + &
              Gnm(1,ri-2,rj  ,ii) * real(ri-2)*coe2B + &
              Gnm(1,ri-1,rj  ,ii) * codB(ii)
              Gnm(1,ri  ,rj+1,ii) = Gnm(1,ri,rj+1,ii) - &
              Gnm(1,ri-2,rj  ,ii) * real(ri-2)*coe3A + &
              Gnm(1,ri-1,rj  ,ii) * coe1A(ii)
            end do
          end if
        end do
        !---------------------------------------------------------------------
        ! use G(n+1,m) recursion only, codA and codB are asymmetric
        do rk = 2, fackl1(ii)
          do ri = 2, facij1(ii)
            if (ri == 2) then
              do rj = 1, rk
                Gnm(2,rk  ,rj  ,ii) = Gnm(2,rk,rj,ii) + &
                Gnm(1,rk  ,rj  ,ii) * codA(ii)
                Gnm(2,rk  ,rj+1,ii) = Gnm(2,rk,rj+1,ii) + &
                Gnm(1,rk  ,rj  ,ii) * coe1B(ii) + &
                Gnm(1,rk-1,rj  ,ii) * real(rk-1)*coe2AB
              end do
            else
              do rj = 1, ri + rk - 2
                Gnm(ri  ,rk  ,rj  ,ii) = &
                Gnm(ri  ,rk  ,rj  ,ii) + &
                Gnm(ri-2,rk  ,rj  ,ii) * real(ri-2)*coe2A +&
                Gnm(ri-1,rk  ,rj  ,ii) * codA(ii)
                Gnm(ri  ,rk  ,rj+1,ii) = &
                Gnm(ri  ,rk  ,rj+1,ii) - &
                Gnm(ri-2,rk  ,rj  ,ii) * real(ri-2)*coe3B +&
                Gnm(ri-1,rk  ,rj  ,ii) * coe1B(ii) + &
                Gnm(ri-1,rk-1,rj  ,ii) * real(rk-1)*coe2AB
              end do
            end if
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
              Itrans(rk,rj,ii) = Itrans(rk,rj,ii) + &
              binomialcoe(facj(ii),ri) * (codi(ii)-codj(ii))**(ri) * &
              Gnm(facij1(ii)-ri, fackl1(ii)-(rk-1),rj,ii)
            end do
          end do
        end do
        do ri = 0, facl(ii)
          do rj = 1, facij1(ii) + fackl1(ii)
            I(rj,ii) = I(rj,ii) + binomialcoe(facl(ii),ri) * &
            (codk(ii)-codl(ii))**(ri) * Itrans(ri+1,rj,ii)
          end do
        end do
      end do
      !---------------------------------------------------------------------
      ! product to PL
      PL = 0.0_dp
      do ri = 0, facij1(1)+fackl1(1)-1
        do rj = 0, facij1(2)+fackl1(2)-1
          do rk = 0, facij1(3)+fackl1(3)-1
            PL(ri+rj+rk+1) = PL(ri+rj+rk+1) + &
            I(ri+1,1) * I(rj+1,2) * I(rk+1,3)
          end do
        end do
      end do
      PL = PL * 2.0_dp * dsqrt(rou/pi)
      !---------------------------------------------------------------------
      ! integral of t exp(-X*t^2)*PL(t^2), 0 -> 1
      int = 0.0_dp
      if (abs(X) < 1E-13) then 
        do ri = 1, nt+6
          int_mic = 1.0_dp / (real(2*ri-2)+1.0_dp)
          int = int + PL(ri) * int_mic
        end do
      else if (abs(X) <= xts) then
        do ri = 1, nt+6
          int_mic = 0.0_dp
          if (2*ri-2 == 0) then
            do rk = 1, tayeps
              int_mic = int_mic + X**(rk-1)*intTaycoe(rk,1)
            end do
          else if (2*ri-2 == 1) then
            do rk = 1, tayeps
              int_mic = int_mic + X**(rk-1)*intTaycoe(rk,2)
            end do
          else
            do rk = 1, tayeps
              int_mic = int_mic + X**(rk-1)*intTaycoe(rk,2*ri-1)
            end do
          end if
          int = int + PL(ri) * int_mic
        end do
      else
        supp1 = dsqrt(pi/X) * erf(dsqrt(X)) / 2.0_dp
        supp3 = 1.0_dp / (2.0_dp*X)
        supp4 = -exp(-X)
        supp2 = (1.0_dp+supp4) * supp3
        do ri = 1, nt+6
          do rj = 0, 2*ri - 2
            if (rj == 0) then
              int_mic = supp1
            else if(rj == 1) then
              int_mic_ = int_mic
              int_mic = supp2
            else
              int_mic__ = int_mic_
              int_mic_ = int_mic
              int_mic = (supp4+real(rj-1)*int_mic__) * supp3
            end if
          end do
          int = int + PL(ri) * int_mic
        end do
      end if
    end if
    return
  end function Integral_V_2e

!------------------------------------------------------------
!> initialising Integral_V_1e and Integral_V_2e
  subroutine Integral_V_init()
    implicit none
    integer             :: ii, jj
    real(dp)            :: expTaycoe(64)
    real(dp)            :: erfTaycoe(64)
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
  end subroutine Integral_V_init
  
end module Hamiltonian