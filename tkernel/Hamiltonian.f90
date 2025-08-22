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
  !DIR$ ATTRIBUTES ALIGN:align_size :: i_j, i_j_s, m_i_j
  real(dp),allocatable    :: i_j(:,:)      ! <AOi|AOj>
  real(dp),allocatable    :: i_j_s(:,:)    ! <AOi|AOj> (sphe-har, reverse then)
  real(dp),allocatable    :: m_i_j(:,:)    ! <AOi|AOjm>

! <AOi|p^2|AOj> related
  !DIR$ ATTRIBUTES ALIGN:align_size :: i_p2_j, AO2p2, evl_p2, exi_T_j
  real(dp),allocatable    :: i_p2_j(:,:)   ! <AOi|p^2|AOj>
  ! unitary transformation from AO basis to p^2 eigenstate (validated)
  real(dp),allocatable    :: AO2p2(:,:)    ! transformation matrix from AO to p2
  ! all eigenvalues of <AOi|p^2|AOj> found by dsyevr
  real(dp),allocatable    :: evl_p2(:)     ! eigenvalue of i_p2_j
  complex(dp),allocatable :: exi_T_j(:,:)  ! extended i_p2_j matrix
  
! <AOi|V|AOj> related
  !DIR$ ATTRIBUTES ALIGN:align_size :: i_V_j, exi_V_j
  real(dp),allocatable    :: i_V_j(:,:)    ! <AOi|V|AOj>
  complex(dp),allocatable :: exi_V_j(:,:)  ! extended i_V_j matrix
  
  !DIR$ ATTRIBUTES ALIGN:align_size :: pxVpx, pyVpy, pzVpz
  !DIR$ ATTRIBUTES ALIGN:align_size :: pxVpy, pyVpx, pxVpz
  !DIR$ ATTRIBUTES ALIGN:align_size :: pzVpx, pyVpz, pzVpy, exSOC
! DKH2 related
  real(dp),allocatable    :: pxVpx(:,:)    ! <AOi|pVp|AOj>
  real(dp),allocatable    :: pyVpy(:,:)
  real(dp),allocatable    :: pzVpz(:,:)
  real(dp),allocatable    :: pxVpy(:,:)
  real(dp),allocatable    :: pyVpx(:,:)
  real(dp),allocatable    :: pxVpz(:,:)
  real(dp),allocatable    :: pzVpx(:,:)
  real(dp),allocatable    :: pyVpz(:,:)
  real(dp),allocatable    :: pzVpy(:,:)
  complex(dp),allocatable :: exSOC(:,:)    ! extended SOC matrix
  
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

! V_integral_1e and V_integral_2e related
  !DIR$ ATTRIBUTES ALIGN:align_size :: NX_mic, NXint, NXint_
  real(dp)                :: NX_mic(32,43)
  real(dp)                :: NXint(32)
  real(dp)                :: NXint_(32)


  contains

!-------------------------------------------------------------------
!> calculate one-electron integrals in single-conf calculations
  subroutine DKH_Hamiltonian()
    implicit none
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
    write(60,'(A)') '  integration initialize'
    call int_init()
    write(60,'(A)') '  complete! stored in: NX_mic, NXint, NXint_'
    write(60,'(A)') '  ----------<1E INTEGRALS>----------'
    
!------------------<NONRELATIVISTIC HAMILTONIAN>------------------
    if (DKH_order == 0) then
      write(60,'(A)') &
      '  spinor basis causes additional cost in scalar SCF.'
      !-----------------------------------------------
      ! 1e integral calculation 
      write(60,'(A)') '  one electron integral calculation'
      call assign_matrices_1e()
      write(60,'(A)') '  complete! stored in:'
      write(60,'(A)') '  i_j, i_p2_j, i_V_j'
      ! cbdm -> sbdm -> fbdm
      write(60,"(A)") '  perform basis transformation'
      if (guess_type == 'molden') then
        allocate(i_j_s(cbdm,cbdm))
        i_j_s = i_j
      end if
      call assign_csf(i_j)
      if (guess_type == 'molden') call csgo(i_j_s)
      call cfgo(i_V_j)
      call cfgo(i_p2_j)
      if (fbdm == sbdm) then
        write(60,"(A)") '  complete! by symm_orth.'
      else
        write(60,"(A,I4)") '  complete! by can_orth. fbdm = ', fbdm
      end if
      
!------------------<DKH2 HAMILTONIAN>------------------
    else if (DKH_order == 2) then
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
      write(60,'(A)') '  one electron integral calculation'
      call assign_matrices_1e()
      write(60,'(A)') '  complete! stored in:'
      write(60,'(A)') '  i_j, i_p2_j, i_V_j, i_pVp_j (9 matrices)'
      if (SRTP_type) write(60,'(A)') '  i_pppVp_j (9 matrices)'
      ! cbdm -> sbdm -> fbdm
      write(60,"(A)") '  perform basis transformation'
      if (guess_type == 'molden') then
        allocate(i_j_s(cbdm,cbdm))
        i_j_s = i_j
      end if
      call assign_csf(i_j)
      if (guess_type == 'molden') call csgo(i_j_s)
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
      if (SRTP_type) then
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
      do loop_i = 1, fbdm
        if (evl_p2(loop_i) < 0.0) &
        call terminate('evl(T) less than zero, may due to code error')
      end do
      ! (AO2p2)^T(i_p2_j)(AO2p2)=evl_p2
      write(60,'(A,I4,A)') '  complete!', fbdm, ' eigenvalues found.'
    end if
    write(60,'(A)') 'exit module Hamiltonian'
  end subroutine DKH_Hamiltonian
  
!-----------------------------------------------------------------------
!> assign value to one electron integral matrices
  subroutine assign_matrices_1e()
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
    type threadlocal    ! thread-local storage to avoid thread-sync overhead
      real(dp), allocatable :: i_j(:,:), i_p2_j(:,:), i_V_j(:,:), pxVpx(:,:),&
      pyVpy(:,:), pzVpz(:,:), pxVpy(:,:), pyVpx(:,:), pyVpz(:,:), pzVpy(:,:),&
      pxVpz(:,:), pzVpx(:,:), px3Vpx(:,:), py3Vpy(:,:), pz3Vpz(:,:),&
      px3Vpy(:,:), py3Vpx(:,:), px3Vpz(:,:), pz3Vpx(:,:),py3Vpz(:,:),pz3Vpy(:,:)
    end type
    type(threadlocal) :: tl
    allocate(i_j(cbdm,cbdm),i_V_j(cbdm,cbdm),i_p2_j(cbdm,cbdm), source=0.0_dp)
    if(DKH_order == 2) then
      allocate(pxVpx(cbdm,cbdm),pyVpy(cbdm,cbdm),pzVpz(cbdm,cbdm),source=0.0_dp)
      allocate(pxVpy(cbdm,cbdm),pyVpx(cbdm,cbdm),pyVpz(cbdm,cbdm),source=0.0_dp)
      allocate(pzVpy(cbdm,cbdm),pxVpz(cbdm,cbdm),pzVpx(cbdm,cbdm),source=0.0_dp)
      if(SRTP_type) then
        allocate (px3Vpx(cbdm,cbdm),py3Vpy(cbdm,cbdm),source=0.0_dp)
        allocate (pz3Vpz(cbdm,cbdm),px3Vpy(cbdm,cbdm),source=0.0_dp)
        allocate (py3Vpx(cbdm,cbdm),px3Vpz(cbdm,cbdm),source=0.0_dp)
        allocate (pz3Vpx(cbdm,cbdm),py3Vpz(cbdm,cbdm),source=0.0_dp)
        allocate (pz3Vpy(cbdm,cbdm),source=0.0_dp)
      end if
    end if
    ! parallel zone, running results consistent with serial
    !$omp parallel num_threads(threads) default(shared) private(i,j,loop_i,  &
    !$omp& loop_j,loop_k,loop_m,contri,Li,Mi,contrj,Lj,Mj,expi,expj,coei,coej, &
    !$omp& codi,codj,facdx_i,facdy_i,facdz_i,coedx_i,coedy_i,coedz_i,facdx_j,  &
    !$omp& facdy_j,facdz_j,coedx_j,coedy_j,coedz_j,tl) if(threads < nproc)
    allocate(&
    tl%i_j(cbdm,cbdm),tl%i_V_j(cbdm,cbdm),tl%i_p2_j(cbdm,cbdm), source=0.0_dp)
    if (DKH_order == 2) then
      allocate(&
      tl%pxVpx(cbdm,cbdm),tl%pyVpy(cbdm,cbdm),tl%pzVpz(cbdm,cbdm),source=0.0_dp)
      allocate(&
      tl%pxVpy(cbdm,cbdm),tl%pyVpx(cbdm,cbdm),tl%pyVpz(cbdm,cbdm),source=0.0_dp)
      allocate(&
      tl%pzVpy(cbdm,cbdm),tl%pxVpz(cbdm,cbdm),tl%pzVpx(cbdm,cbdm),source=0.0_dp)
      if (SRTP_type) then
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
        loop_i = i
        contri = basis_inf(loop_i) % contr
        Li     = basis_inf(loop_i) % L
        Mi     = basis_inf(loop_i) % M
        expi(1:contri) = basis_inf(loop_i) % expo(1:contri)
        coei(1:contri) = basis_inf(loop_i) % Ncoe(1:contri,Mi)
        codi = basis_inf(loop_i) % pos
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
        loop_j = j
        contrj = basis_inf(loop_j) % contr
        Lj     = basis_inf(loop_j) % L
        Mj     = basis_inf(loop_j) % M
        expj(1:contrj) = basis_inf(loop_j) % expo(1:contrj)
        coej(1:contrj) = basis_inf(loop_j) % Ncoe(1:contrj,Mj)
        codj = basis_inf(loop_j) % pos
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
        ! coeffjcjent of x^m*d(exp)
        coedx_j(1:contrj) = -2.0_dp * coej(1:contrj) * expj(1:contrj)
        coedy_j(1:contrj) = -2.0_dp * coej(1:contrj) * expj(1:contrj)
        coedz_j(1:contrj) = -2.0_dp * coej(1:contrj) * expj(1:contrj)
        !---------------------------------
        ! coeffjcjent of d(x^m)*exp
        coedx_j(contrj+1:2*contrj) = AO_fac(1,Lj,Mj) * coej(1:contrj)
        coedy_j(contrj+1:2*contrj) = AO_fac(2,Lj,Mj) * coej(1:contrj)
        coedz_j(contrj+1:2*contrj) = AO_fac(3,Lj,Mj) * coej(1:contrj)
        ! calc <AOi|V|AOj>
        tl%i_V_j(loop_i,loop_j) = calc_1e_V(&
        contri,contrj,coei(1:contri),coej(1:contrj),&
        AO_fac(:,Li,Mi),AO_fac(:,Lj,Mj),expi(1:contri),expj(1:contrj),codi,codj)
        
        if (DKH_order == 2) then
          ! calc <AOi|pVp|AOj>,totally 9 matrices,6 matrices will be calculated
          tl%pxVpx(loop_i,loop_j)=calc_1e_pVp(AO_fac(1,Li,Mi),AO_fac(1,Lj,Mj),&
          contri,contrj,coedx_i(1:2*contri),coedx_j(1:2*contrj),&
          facdx_i,facdx_j,expi(1:contri),expj(1:contrj),codi,codj)
          tl%pyVpy(loop_i,loop_j)=calc_1e_pVp(AO_fac(2,Li,Mi),AO_fac(2,Lj,Mj),&
          contri,contrj,coedy_i(1:2*contri),coedy_j(1:2*contrj),&
          facdy_i,facdy_j,expi(1:contri),expj(1:contrj),codi,codj)
          tl%pzVpz(loop_i,loop_j)=calc_1e_pVp(AO_fac(3,Li,Mi),AO_fac(3,Lj,Mj),&
          contri,contrj,coedz_i(1:2*contri),coedz_j(1:2*contrj),&
          facdz_i,facdz_j,expi(1:contri),expj(1:contrj),codi,codj)
          tl%pxVpy(loop_i,loop_j)=calc_1e_pVp(AO_fac(1,Li,Mi),AO_fac(2,Lj,Mj),&
          contri,contrj,coedx_i(1:2*contri),coedy_j(1:2*contrj),&
          facdx_i,facdy_j,expi(1:contri),expj(1:contrj),codi,codj)
          tl%pyVpx(loop_j,loop_i) = tl%pxVpy(loop_i,loop_j)
          tl%pxVpz(loop_i,loop_j)=calc_1e_pVp(AO_fac(1,Li,Mi),AO_fac(3,Lj,Mj),&
          contri,contrj,coedx_i(1:2*contri),coedz_j(1:2*contrj),&
          facdx_i,facdz_j,expi(1:contri),expj(1:contrj),codi,codj)
          tl%pzVpx(loop_j,loop_i) = tl%pxVpz(loop_i,loop_j)
          tl%pyVpz(loop_i,loop_j)=calc_1e_pVp(AO_fac(2,Li,Mi),AO_fac(3,Lj,Mj),&
          contri,contrj,coedy_i(1:2*contri),coedz_j(1:2*contrj),&
          facdy_i,facdz_j,expi(1:contri),expj(1:contrj),codi,codj)
          tl%pzVpy(loop_j,loop_i) = tl%pyVpz(loop_i,loop_j)
          ! calc <AOi|p^3Vp|AOj> and <AOi|pVp^3|AOj>, totally 18 matrices,
          ! 9 matrices will be calculated
          if (SRTP_type) then
            tl%px3Vpx(loop_i,loop_j) = calc_1e_pppVp(1,AO_fac(1,Li,Mi),&
            AO_fac(1,Lj,Mj),contri,contrj,coedx_i(1:2*contri),&
            coedx_j(1:2*contrj),facdx_i,facdx_j,expi(1:contri),&
            expj(1:contrj),codi,codj)
            tl%py3Vpy(loop_i,loop_j) = calc_1e_pppVp(2,AO_fac(2,Li,Mi),&
            AO_fac(2,Lj,Mj),contri,contrj,coedy_i(1:2*contri),&
            coedy_j(1:2*contrj),facdy_i,facdy_j,expi(1:contri),&
            expj(1:contrj),codi,codj)
            tl%pz3Vpz(loop_i,loop_j) = calc_1e_pppVp(3,AO_fac(3,Li,Mi),&
            AO_fac(3,Lj,Mj),contri,contrj,coedz_i(1:2*contri),&
            coedz_j(1:2*contrj),facdz_i,facdz_j,expi(1:contri),&
            expj(1:contrj),codi,codj)
            tl%px3Vpy(loop_i,loop_j) = calc_1e_pppVp(1,AO_fac(1,Li,Mi),&
            AO_fac(2,Lj,Mj),contri,contrj,coedx_i(1:2*contri),&
            coedy_j(1:2*contrj),facdx_i,facdy_j,expi(1:contri),&
            expj(1:contrj),codi,codj)
            tl%py3Vpx(loop_i,loop_j) = calc_1e_pppVp(2,AO_fac(2,Li,Mi),&
            AO_fac(1,Lj,Mj),contri,contrj,coedy_i(1:2*contri),&
            coedx_j(1:2*contrj),facdy_i,facdx_j,expi(1:contri),&
            expj(1:contrj),codi,codj)
            tl%px3Vpz(loop_i,loop_j) = calc_1e_pppVp(1,AO_fac(1,Li,Mi),&
            AO_fac(3,Lj,Mj),contri,contrj,coedx_i(1:2*contri),&
            coedz_j(1:2*contrj),facdx_i,facdz_j,expi(1:contri),&
            expj(1:contrj),codi,codj)
            tl%pz3Vpx(loop_i,loop_j) = calc_1e_pppVp(3,AO_fac(3,Li,Mi),&
            AO_fac(1,Lj,Mj),contri,contrj,coedz_i(1:2*contri),&
            coedx_j(1:2*contrj),facdz_i,facdx_j,expi(1:contri),&
            expj(1:contrj),codi,codj)
            tl%py3Vpz(loop_i,loop_j) = calc_1e_pppVp(2,AO_fac(2,Li,Mi),&
            AO_fac(3,Lj,Mj),contri,contrj,coedy_i(1:2*contri),&
            coedz_j(1:2*contrj),facdy_i,facdz_j,expi(1:contri),&
            expj(1:contrj),codi,codj)
            tl%pz3Vpy(loop_i,loop_j) = calc_1e_pppVp(3,AO_fac(3,Li,Mi),&
            AO_fac(2,Lj,Mj),contri,contrj,coedz_i(1:2*contri),&
            coedy_j(1:2*contrj),facdz_i,facdy_j,expi(1:contri),&
            expj(1:contrj),codi,codj)
          end if
        end if
        !-------------<basis overlap integral>-------------
        do loop_k = 1, contri
          do loop_m = 1, contrj
            tl%i_j(loop_i,loop_j) = tl%i_j(loop_i,loop_j) + &
            Gaussian_Product_Integral(coei(loop_k)*coej(loop_m),&
            AO_fac(:,Li,Mi),AO_fac(:,Lj,Mj),expi(loop_k),expj(loop_m),codi,codj)
          end do
        end do
        !-------------<P2 integral>-------------
        tl%i_p2_j(loop_i,loop_j) = tl%i_p2_j(loop_i,loop_j) + &
        calc_1e_pp(AO_fac(1,Li,Mi),AO_fac(1,Lj,Mj),contri,contrj,&
        coedx_i(1:2*contri),coedx_j(1:2*contrj),facdx_i,facdx_j,&
        expi(1:contri),expj(1:contrj),codi,codj)
        tl%i_p2_j(loop_i,loop_j) = tl%i_p2_j(loop_i,loop_j) + &
        calc_1e_pp(AO_fac(2,Li,Mi),AO_fac(2,Lj,Mj),contri,contrj,&
        coedy_i(1:2*contri),coedy_j(1:2*contrj),facdy_i,facdy_j,&
        expi(1:contri),expj(1:contrj),codi,codj)
        tl%i_p2_j(loop_i,loop_j) = tl%i_p2_j(loop_i,loop_j) + &
        calc_1e_pp(AO_fac(3,Li,Mi),AO_fac(3,Lj,Mj),contri,contrj,&
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
    if (DKH_order == 2) then
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
      if (SRTP_type) then
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
    if (DKH_order == 2) then
      deallocate(tl%pxVpx, tl%pyVpy, tl%pzVpz)
      deallocate(tl%pxVpy, tl%pyVpx, tl%pxVpz)
      deallocate(tl%pzVpx, tl%pyVpz, tl%pzVpy)
      if (SRTP_type) then
        deallocate(tl%px3Vpx, tl%py3Vpy, tl%pz3Vpz)
        deallocate(tl%px3Vpy, tl%py3Vpx, tl%px3Vpz)
        deallocate(tl%pz3Vpx, tl%py3Vpz, tl%pz3Vpy)
      end if
    end if
    !$omp end parallel
  end subroutine assign_matrices_1e

!-----------------------------------------------------------------------
!> calc values to matrix <AOi|V|AOj>
  real(dp) pure function calc_1e_pp(&
  ni,nj,contri,contrj,coei,coej,faci,facj,expi,expj,codi,codj) result(val)
    implicit none
    integer,intent(in)  :: ni              ! number of x/y/z in <i|
    integer,intent(in)  :: nj              ! number of x/y/z in <j|
    integer,intent(in)  :: contri          ! contr of |i>
    integer,intent(in)  :: contrj          ! contr of |j>
    real(dp),intent(in) :: coei(2*contri)  ! coefficients of i^th orbital
    real(dp),intent(in) :: coej(2*contrj)  ! coefficients of j^th orbital
    integer,intent(in)  :: faci(3,2)       ! xyz factor of i^th orbital
    integer,intent(in)  :: facj(3,2)       ! xyz factor of j^th orbital
    real(dp),intent(in) :: expi(2*contri)  ! expo of i^th orbital
    real(dp),intent(in) :: expj(2*contrj)  ! expo of j^th orbital
    real(dp),intent(in) :: codi(3)         ! centrol coordintates.i^th orbital
    real(dp),intent(in) :: codj(3)         ! centrol coordintates.j^th orbital
    integer             :: numi, numj
    integer :: pploop_i,pploop_j,pploop_k,pploop_l ! loop variables

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
    do pploop_k = 1, numi
      do pploop_i = 1, contri
        do pploop_l = 1, numj
          do pploop_j = 1, contrj
            val = val + Gaussian_Product_Integral(    &
            coei((pploop_k-1)*contri+pploop_i)*       &
            coej((pploop_l-1)*contrj+pploop_j),       &
            faci(:,pploop_k),                         &
            facj(:,pploop_l),                         &
            expi(pploop_i),                           &
            expj(pploop_j),                           &
            codi,                                     &
            codj)
          end do
        end do
      end do
    end do
    return
  end function calc_1e_pp

!-----------------------------------------------------------------------
!> calc values to matrix <AOi|V|AOj>
  real(dp) pure function calc_1e_V(&
  contri,contrj,coei,coej,faci,facj,expi,expj,codi,codj) result(val)
    implicit none
    integer,intent(in)  :: contri          ! contr of |i>
    integer,intent(in)  :: contrj          ! contr of |j>
    real(dp),intent(in) :: coei(contri)    ! coefficients of i^th orbital
    real(dp),intent(in) :: coej(contrj)    ! coefficients of j^th orbital
    integer,intent(in)  :: faci(3)         ! xyz factor of i^th orbital
    integer,intent(in)  :: facj(3)         ! xyz factor of j^th orbital
    real(dp),intent(in) :: expi(contri)    ! expo of i^th orbital
    real(dp),intent(in) :: expj(contrj)    ! expo of j^th orbital
    real(dp),intent(in) :: codi(3)         ! centrol coordintates.i^th orbital
    real(dp),intent(in) :: codj(3)         ! centrol coordintates.j^th orbital
    real(dp)            :: bcodi(3)        ! extended coordintates.i^th orbital
    real(dp)            :: bcodj(3)        ! extended coordintates.j^th orbital
    real(dp)            :: Z_pot, R_pot
    real(dp)            :: codpot(3)
    integer :: bloop_i,bloop_j,bloop_pot   ! loop variables only for calc_1e_V

    val = 0.0_dp
    do bloop_pot = 1, atom_count
      codpot = mol(bloop_pot) % pos
      Z_pot = real(mol(bloop_pot) % atom_number)
      R_pot = mol(bloop_pot) % rad / fm2Bohr
      ! center of potential atom set to zero
      bcodi(:) = codi(:) - codpot(:)
      bcodj(:) = codj(:) - codpot(:)
      do bloop_i = 1, contri
        do bloop_j = 1, contrj
          val = val + V_Integral_1e(        &
          Z_pot,                            &
          coei(bloop_i),                    &
          coej(bloop_j),                    &
          faci,                             &
          facj,                             &
          expi(bloop_i),                    &
          expj(bloop_j),                    &
          bcodi,                            &
          bcodj,                            &
          R_pot)
        end do
      end do
    end do
    return
  end function calc_1e_V
  
!-----------------------------------------------------------------------
!> calc values to matrix <AOi|pVp|AOj> (9 matrices)
!!
!! pxVpx pyVpy pzVpz pxVpy pyVpx pxVpz pzVpx pyVpz pzVpy
  real(dp) pure function calc_1e_pVp(&
  ni,nj,contri,contrj,coei,coej,faci,facj,expi,expj,codi,codj) result(val)
    implicit none
    integer,intent(in)    :: ni            ! number of x/y/z in <i|
    integer,intent(in)    :: nj            ! number of x/y/z in <j|
    integer,intent(in)    :: contri        ! contr of |i>
    integer,intent(in)    :: contrj        ! contr of |j>
    ! coefficients of first-order derivative of i^th orbital
    real(dp),intent(in)   :: coei(2*contri)
    ! coefficients of first-order derivative of j^th orbital
    real(dp),intent(in)   :: coej(2*contrj)
    integer,intent(in)    :: faci(3,2)     ! xyz factor of i^th orbital
    integer,intent(in)    :: facj(3,2)     ! xyz factor of j^th orbital
    real(dp),intent(in)   :: expi(contri)  ! expo of i^th orbital
    real(dp),intent(in)   :: expj(contrj)  ! expo of j^th orbital
    real(dp),intent(in)   :: codi(3)       ! centrol coordintates.i^th orbital
    real(dp),intent(in)   :: codj(3)       ! centrol coordintates.j^th orbital
    real(dp)              :: scodi(3)      ! extended coordintates.i^th orbital
    real(dp)              :: scodj(3)      ! extended coordintates.j^th orbital
    real(dp)              :: Z_pot, R_pot
    real(dp)              :: codpot(3)
    integer               :: numi, numj
    integer :: sloop_i,sloop_j,sloop_k,sloop_l,sloop_pot   ! loop variables

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
    do sloop_pot = 1, atom_count
      codpot = mol(sloop_pot) % pos
      Z_pot = real(mol(sloop_pot) % atom_number)
      R_pot = mol(sloop_pot) % rad / fm2Bohr
      ! center of potential atom set to zero
      scodi(:) = codi(:) - codpot(:)
      scodj(:) = codj(:) - codpot(:)
      do sloop_k = 1, numi
        do sloop_i = 1, contri
          do sloop_l = 1, numj
            do sloop_j = 1, contrj
              val = val + V_Integral_1e(                &
              Z_pot,                                    &
              coei((sloop_k-1)*contri+sloop_i),         &
              coej((sloop_l-1)*contrj+sloop_j),         &
              faci(:,sloop_k),                          &
              facj(:,sloop_l),                          &
              expi(sloop_i),                            &
              expj(sloop_j),                            &
              scodi,                                    &
              scodj,                                    &
              R_pot)
            end do
          end do
        end do
      end do
    end do
    return
  end function calc_1e_pVp
  
!-----------------------------------------------------------------------
!> calc values to matrix <AOi|p^3Vp|AOj> (9 matrices)
!!
!! px3Vpx py3Vpy pz3Vpz px3Vpy py3Vpx px3Vpz pz3Vpx py3Vpz pz3Vpy
  real(dp) pure function calc_1e_pppVp(&
  di,ni,nj,contri,contrj,coei,coej,faci,facj,expi,expj,codi,codj) result(val)
    implicit none
    integer,intent(in)  :: di              ! 1:d|i>/dx, 1:d|i>/dy, 1:d|i>/dz
    integer,intent(in)  :: ni              ! number of x/y/z in <i|
    integer,intent(in)  :: nj              ! number of x/y/z in <j|
    integer,intent(in)  :: contri          ! contr of |i>
    integer,intent(in)  :: contrj          ! contr of |j>
    real(dp),intent(in) :: coei(2*contri)  ! coefficients of i^th orbital
    real(dp),intent(in) :: coej(2*contrj)  ! coefficients of j^th orbital
    integer,intent(in)  :: faci(3,2)       ! xyz factor of i^th orbital
    integer,intent(in)  :: facj(3,2)       ! xyz factor of j^th orbital
    real(dp),intent(in) :: expi(contri)    ! expo of i^th orbital
    real(dp),intent(in) :: expj(contrj)    ! expo of j^th orbital
    real(dp),intent(in) :: codi(3)         ! centrol coordintates.i^th orbital
    real(dp),intent(in) :: codj(3)         ! centrol coordintates.j^th orbital
    ! coordintates.i^th orbital after 3 derivative actions
    real(dp)            :: tcodi(3)
    ! coordintates.j^th orbital after 3 derivative actions
    real(dp)            :: tcodj(3)
    real(dp)            :: Z_pot, R_pot
    real(dp)            :: codpot(3)
    ! coefficients of i^th orbital after 3 derivative actions
    real(dp)            :: tcoei(8*contri)
    ! xyz factor.i^th orbital after 3 derivative actions
    integer             :: tfaci(3,8)      ! di**max, di**max-2 ... di**0
    integer             :: numi, numj      ! number of calls to  V_integral_1e
    integer :: tloop_i,tloop_j,tloop_k,tloop_l,tloop_pot ! loop variables

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
    do tloop_pot = 1, atom_count
      codpot = mol(tloop_pot) % pos
      Z_pot = real(mol(tloop_pot) % atom_number)
      R_pot = mol(tloop_pot) % rad / fm2Bohr
      ! center of potential atom set to zero
      tcodi(:) = codi(:) - codpot(:)
      tcodj(:) = codj(:) - codpot(:)
      do tloop_k = 1, numi
        do tloop_i = 1, contri
          do tloop_l = 1, numj
            do tloop_j = 1, contrj
              val = val + V_Integral_1e(                             &
              Z_pot,                                                 &
              tcoei((tloop_k-1)*contri+tloop_i),                     &
              coej((tloop_l-1)*contrj+tloop_j),                      &
              tfaci(:,tloop_k),                                      &
              facj(:,tloop_l),                                       &
              expi(tloop_i),                                         &
              expj(tloop_j),                                         &
              tcodi,                                                 &
              tcodj,                                                 &
              R_pot)
            end do
          end do
        end do
      end do
    end do
    return
  end function calc_1e_pppVp
  
  
!-----------------------------------------------------------------------
!> full space integration of product of 2 Gaussian functions in Cartesian
  real(dp) pure function Gaussian_Product_Integral(&
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
    integer             :: gloop_i,gloop_j,gloop_k! loop variables
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
    do gloop_k = 1, 3
      do gloop_i = 0, faci(gloop_k)
        ! integral x^m*exp(-b*(x-x0)^2)
        do gloop_j = 0, faci(gloop_k)-gloop_i+facj(gloop_k)
          if (gloop_j == 0) then
            mic = itm
          else if(gloop_j == 1) then
            mic_ = mic
            mic = cod(gloop_k) * mic_
          else
            mic__ = mic_
            mic_ = mic
            mic = cod(gloop_k)*mic_ + 0.5_dp*(real(gloop_j-1,dp)*mic__)*invexpo
          end if
        end do
        integral(gloop_k) = integral(gloop_k) + &
        binomialcoe(faci(gloop_k),gloop_i)*(-codpos(gloop_k))**(gloop_i)*mic
      end do
    end do
    val = pcoe * integral(1) * integral(2) * integral(3)
    return
  end function Gaussian_Product_Integral
  
!-----------------------------------------------------------------------
!> integration of electron-nuclear attraction potential in Cartesian coordinate:
!!
!! inategral transformation: (x-xi)^m*(x-xj)^n*expo(-b*x^2)*expo(x^2*t^2)dxdt
  real(dp) pure function V_Integral_1e(&
  Z,coei,coej,faci,facj,expi,expj,codi,codj,rn) result(val)
    implicit none
    real(dp),intent(in) :: Z                  ! atomic number of potential atom
    real(dp),intent(in) :: coei               ! coeffcient of |i>
    real(dp),intent(in) :: coej               ! coeffcient of |j>
    integer,intent(in)  :: faci(3)            ! xyz factor of |i>
    integer,intent(in)  :: facj(3)            ! xyz factor of |j>
    real(dp),intent(in) :: expi               ! exponent of |i>
    real(dp),intent(in) :: expj               ! exponent of |j>
    real(dp)            :: expo               ! exponnet of product shell
    real(dp)            :: invexpo            ! inverse of expo
    real(dp)            :: coe                ! coefficient of product shell
    real(dp),intent(in) :: codi(3)            ! coordination of i
    real(dp),intent(in) :: codj(3)            ! coordination of j
    ! finite nuclear potential will not be considered if rn = 0.
    real(dp),intent(in) :: rn
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
    integer             :: vloop_i,vloop_j,vloop_k,vloop_o ! loop variables
    integer             :: vloop_mic,vloop_mic_
    integer             :: max1, max2, max3
    real(dp)            :: tmp, br2, expbr2, invbr2, prec

    select case (2*(sum(faci)+sum(facj) + 4))
      case(0:10)
        tayeps = 5
        xts = 0.01
      case(11:15)
        tayeps = 8
        xts = 0.1
      case(16:25)
        tayeps = 15
        xts = 1.0
      case(26:35)
        tayeps = 20
        xts = 1.0
      case(36:40)
        tayeps = 25
        xts = 3.0
      case(41:50)
        tayeps = 30
        xts = 4.0
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
    expbr2 = exp(-br2)
    ! use binomial expansion for easy storage of t-containing coefficients.
    !--------------------------
    ! integral of x,y,z, generate coefficient exp(-b*((cod(1))^2*t^2)/(t^2+b))
    do vloop_k = 1, 3
      t2pb_xyz(:,vloop_k) = 0.0_dp
      do vloop_i = 0, faci(vloop_k)
        vloop_j = 0
        do vloop_j = 0, facj(vloop_k)
          mic = 0.0_dp
          mic_ = 0.0_dp
          mic__ = 0.0_dp
          max1 = faci(vloop_k)+facj(vloop_k)-vloop_i-vloop_j
          do vloop_mic = 0, max1
            if (vloop_mic == 0) then
              mic(1) = rpi
            else if(vloop_mic == 1) then
              mic_(1) = rpi
              mic(1) = 0.0_dp
              mic(2) = cod(vloop_k) * expo * rpi
            else
              mic__ = mic_
              mic_ = mic
              mic = 0.0_dp
              do vloop_mic_ = 0, vloop_mic
                mic(vloop_mic_+2) = mic(vloop_mic_+2) + &
                0.5_dp * real(vloop_mic-1) * mic__(vloop_mic_+1) + &
                cod(vloop_k) * expo * mic_(vloop_mic_+1)
              end do
            end if
          end do
          t2pb_xyz(:,vloop_k) = t2pb_xyz(:,vloop_k) + &
          binomialcoe(faci(vloop_k),vloop_i) * &
          binomialcoe(facj(vloop_k),vloop_j) * &
          (-codi(vloop_k))**(vloop_i) * (-codj(vloop_k))**(vloop_j) * mic
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
    do vloop_i = 0, max1
      do vloop_j = 0, max2
        tmp = t2pb_xyz(vloop_i+1,1) * t2pb_xyz(vloop_j+1,2)
        do vloop_k = 0, max3
          t2pb(vloop_i+vloop_j+vloop_k+2) = t2pb(vloop_i+vloop_j+vloop_k+2) + &
          tmp * t2pb_xyz(vloop_k+1,3)
        end do
      end do
    end do
    val = 0.0_dp
    max1 = sum(faci) + sum(facj) + 4
    if (abs(br2) < 1E-13) then
      ! k = vloop_i, t2pb(vloop_i + 3) = coe
      do vloop_i = 0, max1
        int = 0.0_dp
        do vloop_j = 0, vloop_i
          tmp = (-1.0_dp)**(vloop_j)*binomialcoe(vloop_i,vloop_j)
          do vloop_k = 0, vloop_i
            int_mic = 1.0_dp / (real(2*vloop_i-vloop_j-vloop_k) + 1.0_dp)
            int = int + tmp * binomialcoe(vloop_i,vloop_k) * int_mic
          end do
        end do
        val = val + 2.0_dp*(-1.0_dp)**(vloop_i) * &
        expo**(-vloop_i-1) * t2pb(vloop_i+2) * int
      end do
    else if (abs(br2) <= xts) then
      ! k = vloop_i, t2pb(vloop_i + 3) = coe
      do vloop_i = 0, max1
        int = 0.0_dp
        do vloop_j = 0, vloop_i
          tmp = (-1.0_dp)**(vloop_j) * binomialcoe(vloop_i,vloop_j)
          do vloop_k = 0, vloop_i
            int_mic = 0.0_dp
            if (2*vloop_i-vloop_j-vloop_k == 0) then
              do vloop_o = 1, tayeps
                int_mic = int_mic + br2**(vloop_o-1)*NXint_(vloop_o)
              end do
            else if (2*vloop_i-vloop_j-vloop_k == 1) then
              do vloop_o = 1, tayeps
                int_mic = int_mic + br2**(vloop_o-1)*NXint(vloop_o)
              end do
            else
              do vloop_o = 1, tayeps
                int_mic = int_mic + br2**(vloop_o-1)*(NXint(vloop_o) + &
                0.5_dp*NX_mic(vloop_o,2*vloop_i-vloop_j-vloop_k-1))
              end do
            end if
            int = int + tmp * binomialcoe(vloop_i,vloop_k) * int_mic
          end do
        end do
        val = val + 2.0_dp*(-1.0_dp)**(vloop_i) * &
        expo**(-vloop_i-1) * t2pb(vloop_i+2) * int
      end do
    else
      !k = vloop_i, t2pb(vloop_i + 3) = coe
      invbr2 = 0.5_dp / br2
      prec = 0.5_dp * dsqrt(pi/br2) * erf(dsqrt(br2))
      do vloop_i = 0, max1
        int = 0.0_dp
        do vloop_j = 0, vloop_i
          tmp = (-1.0_dp)**(vloop_j) * binomialcoe(vloop_i,vloop_j)
          do vloop_k = 0, vloop_i
            do vloop_mic = 0, 2*vloop_i - vloop_j - vloop_k
              if (vloop_mic == 0) then
                int_mic = prec
              else if(vloop_mic == 1) then
                int_mic_ = int_mic
                int_mic = (1.0_dp - expbr2) * invbr2
              else
                int_mic__ = int_mic_
                int_mic_ = int_mic
                int_mic = (-expbr2 + real(vloop_mic - 1) * int_mic__) * invbr2
              end if
            end do
            int = int + tmp * binomialcoe(vloop_i,vloop_k) * int_mic
          end do
        end do
        val = val + 2.0_dp*(-1.0_dp)**(vloop_i) * &
        expo**(-vloop_i-1) * t2pb(vloop_i+2) * int
      end do
    end if
    val = val / rpi
    if (finitenuc) then
      ! Finite nuclear model correction, spherical electric charge
      val = val - &
      2.0_dp/3.0_dp * pi * rn * rn * expbr2 * (-codi(1))**(faci(1)) * &
      (-codi(2))**(faci(2)) * (-codi(3))**(faci(3)) * (-codj(1))**(facj(1)) * &
      (-codj(2))**(facj(2)) * (-codj(3))**(facj(3))
    end if
    val = val * coe
    return
  end function V_Integral_1e
  
!-----------------------------------------------------------------------
!> integration of two electron repulsion potential in Cartesian coordinate
!!
!! EXPRESS: |AOi>:(x1 - xi), |AOj>:(x1 - xj), |AOk>:(x2 - xk), |AOl>:(x2 - xl)
!! xA  xB
!! Li > Lj, Lk > Ll
  pure real(dp) function V_Integral_2e(&
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
    integer             :: rloop_i, rloop_j, rloop_k, ii
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
    if (int(real(nt)/2.0) + 1 <= 5) then
      nroots = int(real(nt)/2.0) + 1
      ! call rys_roots(libcint-like), GRysroots(GAMESS-like)
      call GRysroots(nroots, X, u, w)
      int = 0.0_dp
      do rloop_i = 1, nroots
        !t2 = u(rloop_i) / (rou+u(rloop_i))       ! for rys_roots
        t2 = u(rloop_i)                           ! for GRysroots
        f0 = coe2AB*t2
        f3 = coe2A-coe3B*t2
        f4 = coe2B-coe3A*t2
        !---------------------------Gnm---------------------------
        Gnm = 0.0_dp
        do ii = 1, 3
          f1 = codA(ii)+coe1B(ii)*t2
          f2 = codB(ii)+coe1A(ii)*t2
          Gnm(1,1,1,ii) = Gim(ii)
          do rloop_j = 2, facij1(ii)
            if (rloop_j == 2) then
              Gnm(2,1,1,ii) = Gnm(2,1,1,ii) + &
              Gnm(1,1,1,ii)*f1
            else
              Gnm(rloop_j,1,1,ii) = Gnm(rloop_j,1,1,ii) + &
              Gnm(rloop_j-2,1,1,ii)*real(rloop_j-2)*f3 + &
              Gnm(rloop_j-1,1,1,ii)*f1
            end if
          end do
          do rloop_j = 2, fackl1(ii)
            if (rloop_j == 2) then
              Gnm(1,2,1,ii) = Gnm(1,2,1,ii) + &
              Gnm(1,1,1,ii)*f2
            else
              Gnm(1,rloop_j,1,ii) = Gnm(1,rloop_j,1,ii) + &
              Gnm(1,rloop_j-2,1,ii)*real(rloop_j-2)*f4 + &
              Gnm(1,rloop_j-1,1,ii)*f2
            end if
          end do
          do rloop_j = 2, fackl1(ii)
            do rloop_k = 2, facij1(ii)
              if (rloop_k == 2) then
                Gnm(2,rloop_j,1,ii) = Gnm(2,rloop_j,1,ii) + &
                Gnm(1,rloop_j,1,ii)*f1 + &
                Gnm(1,rloop_j-1,1,ii)*real(rloop_j-1)*f0
              else
                Gnm(rloop_k,rloop_j,1,ii) = Gnm(rloop_k,rloop_j,1,ii) + &
                Gnm(rloop_k-2,rloop_j,1,ii)*real(rloop_k-2)*f3 + &
                Gnm(rloop_k-1,rloop_j,1,ii)*f1 + &
                Gnm(rloop_k-1,rloop_j-1,1,ii)*real(rloop_j-1)*f0
              end if
            end do
          end do
        end do
        !---------------------------I---------------------------
        Itrans = 0.0_dp
        I = 0.0_dp
        do ii = 1, 3
          do rloop_k = 1, facl(ii)+1
            do rloop_j = 0, facj(ii)
              Itrans(rloop_k,1,ii) = Itrans(rloop_k,1,ii) + &
              binomialcoe(facj(ii),rloop_j)*(codi(ii)-codj(ii))**rloop_j*&
              Gnm(facij1(ii)-rloop_j,fackl1(ii)+1-rloop_k,1,ii)
            end do
          end do
          do rloop_j = 0, facl(ii)
            I(1,ii) = I(1,ii) + &
            binomialcoe(facl(ii),rloop_j)*(codk(ii)-codl(ii))**(rloop_j)*&
            Itrans(rloop_j+1,1,ii)
          end do
        end do
        !---------------------------PL---------------------------
        PL(1) = I(1,1) * I(1,2) * I(1,3)
        int = int + w(rloop_i)*PL(1)
      end do
      int = int * 2.0_dp * dsqrt(rou/pi)
    !=============================================================
    ! direct integral for high angular momentum Gaussian functions
    else
      ! Ix(ni+nj,0,nk+nl,0,u)
      select case (2*(nt+6)-2)
        case(0:15)
          tayeps = 8
          xts = 0.01
        case(16:25)
          tayeps = 15
          xts = 1.0
        case(26:35)
          tayeps = 20
          xts = 1.0
        case(36:40)
          tayeps = 25
          xts = 3.0
        case(41:50)
          tayeps = 30
          xts = 4.0
      end select
      Gnm = 0.0_dp
      !---------------------------------------------------------------------
      ! reduce the factor (1-t^2)^(1/2)*exp(-Dx*t^2)
      do ii = 1, 3
        Gnm(1,1,1,ii) = Gim(ii)
        do rloop_i = 2, facij1(ii)
          if (rloop_i == 2) then
            Gnm(2,1,1,ii) = Gnm(2,1,1,ii) + &
            Gnm(1,1,1,ii)*codA(ii)
            Gnm(2,1,2,ii) = Gnm(2,1,2,ii) + &
            Gnm(1,1,1,ii)*coe1B(ii)
          else
            do rloop_j = 1, rloop_i - 1
              Gnm(rloop_i  ,1,rloop_j  ,ii) = Gnm(rloop_i,1,rloop_j,ii) + &
              Gnm(rloop_i-2,1,rloop_j  ,ii) * real(rloop_i-2)*coe2A + &
              Gnm(rloop_i-1,1,rloop_j  ,ii) * codA(ii)
              Gnm(rloop_i  ,1,rloop_j+1,ii) = Gnm(rloop_i,1,rloop_j+1,ii) - &
              Gnm(rloop_i-2,1,rloop_j  ,ii) * real(rloop_i-2)*coe3B + &
              Gnm(rloop_i-1,1,rloop_j  ,ii) * coe1B(ii)
            end do
          end if
        end do
        do rloop_i = 2, fackl1(ii)
          if (rloop_i == 2) then
            Gnm(1,2,1,ii) = Gnm(1,2,1,ii) + &
            Gnm(1,1,1,ii)*codB(ii)
            Gnm(1,2,2,ii) = Gnm(1,2,2,ii) + &
            Gnm(1,1,1,ii)*coe1A(ii)
          else
            do rloop_j = 1, rloop_i - 1
              Gnm(1,rloop_i  ,rloop_j  ,ii) = Gnm(1,rloop_i,rloop_j,ii) + &
              Gnm(1,rloop_i-2,rloop_j  ,ii) * real(rloop_i-2)*coe2B + &
              Gnm(1,rloop_i-1,rloop_j  ,ii) * codB(ii)
              Gnm(1,rloop_i  ,rloop_j+1,ii) = Gnm(1,rloop_i,rloop_j+1,ii) - &
              Gnm(1,rloop_i-2,rloop_j  ,ii) * real(rloop_i-2)*coe3A + &
              Gnm(1,rloop_i-1,rloop_j  ,ii) * coe1A(ii)
            end do
          end if
        end do
        !---------------------------------------------------------------------
        ! use G(n+1,m) recursion only, codA and codB are asymmetric
        do rloop_k = 2, fackl1(ii)
          do rloop_i = 2, facij1(ii)
            if (rloop_i == 2) then
              do rloop_j = 1, rloop_k
                Gnm(2,rloop_k  ,rloop_j  ,ii) = Gnm(2,rloop_k,rloop_j,ii) + &
                Gnm(1,rloop_k  ,rloop_j  ,ii) * codA(ii)
                Gnm(2,rloop_k  ,rloop_j+1,ii) = Gnm(2,rloop_k,rloop_j+1,ii) + &
                Gnm(1,rloop_k  ,rloop_j  ,ii) * coe1B(ii) + &
                Gnm(1,rloop_k-1,rloop_j  ,ii) * real(rloop_k-1)*coe2AB
              end do
            else
              do rloop_j = 1, rloop_i + rloop_k - 2
                Gnm(rloop_i  ,rloop_k  ,rloop_j  ,ii) = &
                Gnm(rloop_i  ,rloop_k  ,rloop_j  ,ii) + &
                Gnm(rloop_i-2,rloop_k  ,rloop_j  ,ii) * real(rloop_i-2)*coe2A +&
                Gnm(rloop_i-1,rloop_k  ,rloop_j  ,ii) * codA(ii)
                Gnm(rloop_i  ,rloop_k  ,rloop_j+1,ii) = &
                Gnm(rloop_i  ,rloop_k  ,rloop_j+1,ii) - &
                Gnm(rloop_i-2,rloop_k  ,rloop_j  ,ii) * real(rloop_i-2)*coe3B +&
                Gnm(rloop_i-1,rloop_k  ,rloop_j  ,ii) * coe1B(ii) + &
                Gnm(rloop_i-1,rloop_k-1,rloop_j  ,ii) * real(rloop_k-1)*coe2AB
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
        do rloop_k = 1, facl(ii) + 1
          do rloop_i = 0, facj(ii)
            do rloop_j = 1, facij1(ii) + fackl1(ii)
              Itrans(rloop_k,rloop_j,ii) = Itrans(rloop_k,rloop_j,ii) + &
              binomialcoe(facj(ii),rloop_i) * (codi(ii)-codj(ii))**(rloop_i) * &
              Gnm(facij1(ii)-rloop_i, fackl1(ii)-(rloop_k-1),rloop_j,ii)
            end do
          end do
        end do
        do rloop_i = 0, facl(ii)
          do rloop_j = 1, facij1(ii) + fackl1(ii)
            I(rloop_j,ii) = I(rloop_j,ii) + binomialcoe(facl(ii),rloop_i) * &
            (codk(ii)-codl(ii))**(rloop_i) * Itrans(rloop_i+1,rloop_j,ii)
          end do
        end do
      end do
      !---------------------------------------------------------------------
      ! product to PL
      PL = 0.0_dp
      do rloop_i = 0, facij1(1)+fackl1(1)-1
        do rloop_j = 0, facij1(2)+fackl1(2)-1
          do rloop_k = 0, facij1(3)+fackl1(3)-1
            PL(rloop_i+rloop_j+rloop_k+1) = PL(rloop_i+rloop_j+rloop_k+1) + &
            I(rloop_i+1,1) * I(rloop_j+1,2) * I(rloop_k+1,3)
          end do
        end do
      end do
      PL = PL * 2.0_dp * dsqrt(rou/pi)
      !---------------------------------------------------------------------
      ! integral of t exp(-X*t^2)*PL(t^2), 0 -> 1
      int = 0.0_dp
      if (abs(X) < 1E-13) then 
        do rloop_i = 1, nt+6
          int_mic = 1.0_dp / (real(2*rloop_i-2)+1.0_dp)
          int = int + PL(rloop_i) * int_mic
        end do
      else if (abs(X) <= xts) then
        do rloop_i = 1, nt+6
          int_mic = 0.0_dp
          if (2*rloop_i-2 == 0) then
            do rloop_k = 1, tayeps
              int_mic = int_mic + X**(rloop_k-1)*NXint_(rloop_k)
            end do
          else if (2*rloop_i-2 == 1) then
            do rloop_k = 1, tayeps
              int_mic = int_mic + X**(rloop_k-1)*NXint(rloop_k)
            end do
          else
            do rloop_k = 1, tayeps
              int_mic = int_mic + X**(rloop_k-1)*(NXint(rloop_k) + &
              0.5_dp*NX_mic(rloop_k,2*rloop_i-3))
            end do
          end if
          int = int + PL(rloop_i) * int_mic
        end do
      else
        supp1 = dsqrt(pi/X) * erf(dsqrt(X)) / 2.0_dp
        supp3 = 1.0_dp / (2.0_dp*X)
        supp4 = -exp(-X)
        supp2 = (1.0_dp+supp4) * supp3
        do rloop_i = 1, nt+6
          do rloop_j = 0, 2*rloop_i - 2
            if (rloop_j == 0) then
              int_mic = supp1
            else if(rloop_j == 1) then
              int_mic_ = int_mic
              int_mic = supp2
            else
              int_mic__ = int_mic_
              int_mic_ = int_mic
              int_mic = (supp4+real(rloop_j-1)*int_mic__) * supp3
            end if
          end do
          int = int + PL(rloop_i) * int_mic
        end do
      end if
    end if
    return
  end function V_Integral_2e

!------------------------------------------------------------
!> initialising V_integral_1e and V_integral_2e
  subroutine int_init()
    implicit none
    integer :: ii, jj
    NX_mic = 0.0_dp
    do ii = 0, 42
      if (ii == 0) then ! dsqrt(pi/X) * erf(dsqrt(X)) / 2.0_dp
        do jj = 1, 32
          NX_mic(jj,ii+1)=(-1.0_dp)**jj*(1.0_dp/(real(2*jj+1)*factorial(jj)))
        end do
      else if(ii == 1) then ! (1.0_dp - exp(-X)) / (2.0_dp * X)
        do jj = 1, 32
          NX_mic(jj,ii+1)=(-1.0_dp)**jj*(1.0_dp/(factorial(jj+1)))
        end do
      else ! (-exp(-X)+real(ii-1)*int_mic__) / (2.0_dp*X)
        do jj = 1, 31
          NX_mic(jj,ii+1)=((-1.0_dp)**jj*(0.5_dp/(factorial(jj+1))) + &
          0.5_dp*(NX_mic(jj+1,ii-1)))*real(ii+1)
        end do
        NX_mic(32,ii+1)=((0.5_dp/(factorial(33))))*real(ii+1)
      end if
    end do
    NXint = 0.0_dp
    NXint_ = 0.0_dp
    do ii = 1, 32
      NXint(ii) = (-1.0_dp)**(ii+1)*(0.5_dp/factorial(ii))
      NXint_(ii) = (-1.0_dp)**(ii+1)*(1.0_dp/(real(2*ii-1)*factorial(ii-1)))
    end do
  end subroutine int_init
  
end module Hamiltonian