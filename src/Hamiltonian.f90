!> @file Hamiltonian.f90
!!
!! @brief electronic static Hamiltonian integrals
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

  integer                 :: info          ! info of calling lapack routines
! <AOi|AOj> related
  !DIR$ ATTRIBUTES ALIGN:align_size :: i_j, Xm, exXm
  real(dp),allocatable    :: i_j(:,:)      ! <AOi|AOj>
  real(dp),allocatable    :: Xm(:,:)       ! orthogonal transform unitary matrix
  real(dp)                :: smallest_evl  ! smallest eigenvalue of i_j
  complex(dp),allocatable :: exXm(:,:)     ! extended Xm matrix

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
  
  !DIR$ ATTRIBUTES ALIGN:align_size :: i_pxVpx_j, i_pyVpy_j, i_pzVpz_j
  !DIR$ ATTRIBUTES ALIGN:align_size :: i_pxVpy_j, i_pyVpx_j, i_pxVpz_j
  !DIR$ ATTRIBUTES ALIGN:align_size :: i_pzVpx_j, i_pyVpz_j, i_pzVpy_j, exSOC
! DKH2 related
  real(dp),allocatable    :: i_pxVpx_j(:,:)! <AOi|pVp|AOj>
  real(dp),allocatable    :: i_pyVpy_j(:,:)
  real(dp),allocatable    :: i_pzVpz_j(:,:)
  real(dp),allocatable    :: i_pxVpy_j(:,:)
  real(dp),allocatable    :: i_pyVpx_j(:,:)
  real(dp),allocatable    :: i_pxVpz_j(:,:)
  real(dp),allocatable    :: i_pzVpx_j(:,:)
  real(dp),allocatable    :: i_pyVpz_j(:,:)
  real(dp),allocatable    :: i_pzVpy_j(:,:)
  complex(dp),allocatable :: exSOC(:,:)    ! extended SOC matrix
  
  !DIR$ ATTRIBUTES ALIGN:align_size :: i_px3Vpx_j, i_py3Vpy_j, i_pz3Vpz_j
  !DIR$ ATTRIBUTES ALIGN:align_size :: i_px3Vpy_j, i_py3Vpx_j, i_px3Vpz_j
  !DIR$ ATTRIBUTES ALIGN:align_size :: i_pz3Vpx_j, i_py3Vpz_j, i_pz3Vpy_j, exSR
  real(dp),allocatable    :: i_px3Vpx_j(:,:)! <AOi|pppVp|AOj>
  real(dp),allocatable    :: i_py3Vpy_j(:,:)
  real(dp),allocatable    :: i_pz3Vpz_j(:,:)
  real(dp),allocatable    :: i_px3Vpy_j(:,:)
  real(dp),allocatable    :: i_py3Vpx_j(:,:)
  real(dp),allocatable    :: i_px3Vpz_j(:,:)
  real(dp),allocatable    :: i_pz3Vpx_j(:,:)
  real(dp),allocatable    :: i_py3Vpz_j(:,:)
  real(dp),allocatable    :: i_pz3Vpy_j(:,:)
  complex(dp),allocatable :: exSR(:,:)    ! extended SR matrix
  
  ! <AOi|pxVpy3|AOj> = Trans(<AOi|py3Vpx|AOj>)
 
  contains

!-------------------------------------------------------------------
!> calculate one-electron integrals in single-conf calculations
  subroutine DKH_Hamiltonian()
    implicit none
    write(60,'(A)') 'Module Hamiltonian:'
    ! OpenMP set up
    write(60,'(A)') '  ----------<PARALLEL>----------'
    cpu_threads = omp_get_num_procs()
    write(60,'(A,I3,A,I3)') &
    '  threads using:',threads_use,' CPU threads:',cpu_threads
    if (cpu_threads <= threads_use) then
      write(*,'(A,I2)') &
      'TRESC: Calculation will perform serially! CPU threads is',cpu_threads
      write(60,'(A)') '  Warning: calculation will be performed SERIALLY!'
    else
      loop_i = omp_get_proc_bind()
      if (loop_i == 0) then
        write(60,'(A)') '  processor bind OFF'
      else
        loop_i = omp_get_num_places()
        write(60,'(A,I2,A)') '  processor bind ON, ',loop_i,' places'
      end if
    end if
    write(60,'(A)') '  ----------<1E INTEGRALS>----------'
    
    
!------------------<NONRELATIVISTIC HAMILTONIAN>------------------
    if (DKH_order == 0) then
      write(60,'(A)') &
      '  spinor basis causes additional cost in scalar SCF.'
      !-----------------------------------------------
      ! 1e integral calculation 
      write(60,'(A)') '  one electron integral calculation'
      call assign_matrices_1e()
      write(60,'(A)') &
      '  complete! stored in:'
      write(60,'(A)') '  i_j, i_p2_j, i_V_j'
      if (s_h) then
        ! transform to spherical-harmonic basis
        write(60,"(A)") '  perform spherical-harmonic basis transformation'
        call assign_cs()
        call sphehar(i_j)
        call sphehar(i_V_j)
        call sphehar(i_p2_j)
        write(60,"(A)") &
        '  complete! spherical-harmonic basis (5D,7F,9G) involved.'
      end if
      ! symmetric orthogonalisation
      write(60,'(A,E10.3)') &
      '  symmetric orthogonalisation of 1e integrals, cutS = ', cutS
      allocate(Xm(sbdm,sbdm))
      call symm_orth(i_j, sbdm, Xm, smallest_evl)
      if (smallest_evl < 0.0) &
      call terminate('evl(i_j) less than zero, may due to code error')
      if (smallest_evl < cutS) then
        write(60,'(A,F9.6,A)') &
        '  complete! smallest eigenvalue = ', smallest_evl, ', less than cutS'
        write(60,'(A)') '  canonical orthogonalisation of 1e integrals'
        deallocate(Xm)
        call can_orth(i_j, sbdm, Xm, fbdm)
      else
        fbdm = sbdm
      end if
      allocate(exXm(2*sbdm,2*fbdm))
      exXm = c0
      do loop_i = 1, sbdm
        do loop_j = 1, fbdm
          exXm(2*loop_i-1,2*loop_j-1) = cmplx(Xm(loop_i,loop_j),0.0_dp,dp)
          exXm(2*loop_i,2*loop_j) = cmplx(Xm(loop_i,loop_j),0.0_dp,dp)
        end do
      end do
      ! transpose(Xm)(Xm)=I
      call detach(i_j)
      call detach(i_p2_j)
      call detach(i_V_j)
      if (smallest_evl < cutS) then
        write(60,'(A,I4)') '  complete! basis dimension reduce to ', fbdm
      else
        write(60,'(A,F9.6)') '  complete! smallest eigenvalue = ', smallest_evl
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
      if (s_h) then
        ! transform to spherical-harmonic basis
        write(60,"(A)") '  perform spherical-harmonic basis transformation'
        call assign_cs()
        call sphehar(i_j)
        call sphehar(i_V_j)
        call sphehar(i_p2_j)
        call sphehar(i_pxVpx_j)
        call sphehar(i_pyVpy_j)
        call sphehar(i_pzVpz_j)
        call sphehar(i_pxVpy_j)
        call sphehar(i_pyVpx_j)
        call sphehar(i_pxVpz_j)
        call sphehar(i_pzVpx_j)
        call sphehar(i_pyVpz_j)
        call sphehar(i_pzVpy_j)
        if (SRTP_type) then
          call sphehar(i_px3Vpx_j)
          call sphehar(i_py3Vpy_j)
          call sphehar(i_pz3Vpz_j)
          call sphehar(i_px3Vpy_j)
          call sphehar(i_py3Vpx_j)
          call sphehar(i_px3Vpz_j)
          call sphehar(i_pz3Vpx_j)
          call sphehar(i_py3Vpz_j)
          call sphehar(i_pz3Vpy_j)
        end if
        write(60,"(A)") &
        '  complete! spherical-harmonic basis (5D,7F,9G) involved.'
      end if
      ! symmetric orthogonalisation
      write(60,'(A)') '  symmetric orthogonalisation of 1e integrals'
      allocate(Xm(sbdm, sbdm))
      call symm_orth(i_j, sbdm, Xm, smallest_evl)
      if (smallest_evl < 0.0) &
      call terminate('evl(i_j) less than zero, may due to code error')
      if (smallest_evl < cutS) then
        if (smallest_evl > safmin) then
          write(60,'(A,F9.6,A)') &
          '  complete! smallest eigenvalue = ', smallest_evl, ', less than cutS'
        end if
        write(60,'(A)') '  canonical orthogonalisation of 1e integrals'
        deallocate(Xm)
        call can_orth(i_j, sbdm, Xm, fbdm)
      else
        fbdm = sbdm
      end if
      allocate(exXm(2*sbdm,2*fbdm))
      exXm = c0
      do loop_i = 1, sbdm
        do loop_j = 1, fbdm
          exXm(2*loop_i-1,2*loop_j-1) = cmplx(Xm(loop_i,loop_j),0.0_dp,dp)
          exXm(2*loop_i,2*loop_j) = cmplx(Xm(loop_i,loop_j),0.0_dp,dp)
        end do
      end do
      ! transpose(Xm)(Xm)=I
      call detach(i_j)
      call detach(i_p2_j)
      call detach(i_V_j)
      !-------------------------------------------
      call detach(i_pxVpx_j)
      call detach(i_pyVpy_j)
      call detach(i_pzVpz_j)
      call detach(i_pxVpy_j)
      call detach(i_pyVpx_j)
      call detach(i_pxVpz_j)
      call detach(i_pzVpx_j)
      call detach(i_pyVpz_j)
      call detach(i_pzVpy_j)
      if (SRTP_type) then
        call detach(i_px3Vpx_j)
        call detach(i_py3Vpy_j)
        call detach(i_pz3Vpz_j)
        call detach(i_px3Vpy_j)
        call detach(i_py3Vpx_j)
        call detach(i_px3Vpz_j)
        call detach(i_pz3Vpx_j)
        call detach(i_py3Vpz_j)
        call detach(i_pz3Vpy_j)
      end if
      if (smallest_evl < cutS) then
        write(60,'(A,I4)') '  complete! basis dimension reduce to ', fbdm
      else
        write(60,'(A,F9.6)') '  complete! smallest eigenvalue = ', smallest_evl
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
    integer :: i                             ! openMP parallel variable
    integer :: contraction_i                 ! contraction of atom_i, shell_i
    integer :: L_i                           ! angular quantum number of |AOi>
    integer :: M_i                           ! magnetic quantum number of |AOi>
    integer :: contraction_j                 ! contraction of atom_j, shell_j
    integer :: L_j                           ! angular quantum number of |AOj>
    integer :: M_j                           ! magnetic quantum number of |AOj>
    real(dp) :: exponents_i(20)              ! exponents of |AOi>
    real(dp) :: exponents_j(20)              ! exponents of |AOj>
    real(dp) :: coefficient_i(20)            ! coefficient of |AOi>
    real(dp) :: coefficient_j(20)            ! coefficient of |AOj>
    real(dp) :: coordinate_i(3)              ! coordinate of center of |AOi>
    real(dp) :: coordinate_j(3)              ! coordinate of center of |AOj>
    real(dp) :: coedx_i(40)                  ! coefficient.derivative x.|AOi>
    real(dp) :: coedy_i(40)                  ! coefficient.derivative y.|AOi>
    real(dp) :: coedz_i(40)                  ! coefficient.derivative z.|AOi>
    real(dp) :: coedx_j(40)                  ! coefficient.derivative x.|AOj>
    real(dp) :: coedy_j(40)                  ! coefficient.derivative y.|AOj>
    real(dp) :: coedz_j(40)                  ! coefficient.derivative z.|AOj>
    character(len=5) :: chardx_i(2)          ! x,y,z factor.derivative x.|AOi>
    character(len=5) :: chardy_i(2)          ! x,y,z factor.derivative y.|AOi>
    character(len=5) :: chardz_i(2)          ! x,y,z factor.derivative z.|AOi>
    character(len=5) :: chardx_j(2)          ! x,y,z factor.derivative x.|AOj>
    character(len=5) :: chardy_j(2)          ! x,y,z factor.derivative y.|AOj>
    character(len=5) :: chardz_j(2)          ! x,y,z factor.derivative z.|AOj>
    real(dp) :: integral_overlap
    real(dp) :: integral_dx
    real(dp) :: integral_dy
    real(dp) :: integral_dz
    real(dp) :: integral_V
    real(dp) :: integral_DKH2(9)
    real(dp) :: integral_SRTP(9)
    allocate (i_j(cbdm,cbdm))
    i_j = 0.0_dp
    allocate (i_V_j(cbdm,cbdm))
    i_V_j = 0.0_dp
    allocate (i_p2_j(cbdm,cbdm))
    i_p2_j = 0.0_dp
    if(DKH_order == 2) then
      allocate (i_pxVpx_j(cbdm,cbdm))
      i_pxVpx_j = 0.0_dp
      allocate (i_pyVpy_j(cbdm,cbdm))
      i_pyVpy_j = 0.0_dp
      allocate (i_pzVpz_j(cbdm,cbdm))
      i_pzVpz_j = 0.0_dp
      allocate (i_pxVpy_j(cbdm,cbdm))
      i_pxVpy_j = 0.0_dp
      allocate (i_pyVpx_j(cbdm,cbdm))
      i_pyVpx_j = 0.0_dp
      allocate (i_pyVpz_j(cbdm,cbdm))
      i_pyVpz_j = 0.0_dp
      allocate (i_pzVpy_j(cbdm,cbdm))
      i_pzVpy_j = 0.0_dp
      allocate (i_pxVpz_j(cbdm,cbdm))
      i_pxVpz_j = 0.0_dp
      allocate (i_pzVpx_j(cbdm,cbdm))
      i_pzVpx_j = 0.0_dp
      if(SRTP_type) then
        allocate (i_px3Vpx_j(cbdm,cbdm))
        i_px3Vpx_j = 0.0_dp
        allocate (i_py3Vpy_j(cbdm,cbdm))
        i_py3Vpy_j = 0.0_dp
        allocate (i_pz3Vpz_j(cbdm,cbdm))
        i_pz3Vpz_j = 0.0_dp
        allocate (i_px3Vpy_j(cbdm,cbdm))
        i_px3Vpy_j = 0.0_dp
        allocate (i_py3Vpx_j(cbdm,cbdm))
        i_py3Vpx_j = 0.0_dp
        allocate (i_px3Vpz_j(cbdm,cbdm))
        i_px3Vpz_j = 0.0_dp
        allocate (i_pz3Vpx_j(cbdm,cbdm))
        i_pz3Vpx_j = 0.0_dp
        allocate (i_py3Vpz_j(cbdm,cbdm))
        i_py3Vpz_j = 0.0_dp
        allocate (i_pz3Vpy_j(cbdm,cbdm))
        i_pz3Vpy_j = 0.0_dp
      end if
    end if
    ! parallel zone, running results consistent with serial
    !$omp parallel num_threads(threads_use) default(shared) private(i,loop_i,&
    !$omp& loop_j,loop_k,loop_m,contraction_i,L_i,M_i,contraction_j,L_j,M_j, &
    !$omp& exponents_i,exponents_j,coefficient_i,coefficient_j,coordinate_i,&
    !$omp& coordinate_j,chardx_i,chardy_i,chardz_i,coedx_i,coedy_i,coedz_i, &
    !$omp& chardx_j,chardy_j,chardz_j,coedx_j,coedy_j,coedz_j,integral_dx,&
    !$omp& integral_dy,integral_dz,integral_overlap,integral_V,integral_DKH2,&
    !$omp& integral_SRTP) if(threads_use < cpu_threads)
    !$omp do schedule(static)
    do i = 1, cbdm
      loop_i = i
      exponents_i = 0.0_dp
      coefficient_i = 0.0_dp
      coedx_i = 0.0_dp
      coedy_i = 0.0_dp
      coedz_i = 0.0_dp
      contraction_i = atom_basis(molecule(basis_inf(loop_i) % atom) % &
      basis_number + basis_inf(loop_i) % shell - 1) % contraction
      L_i = basis_inf(loop_i) % L
      M_i = basis_inf(loop_i) % M
      do loop_k = 1, contraction_i
        exponents_i(loop_k) = atom_basis(molecule(basis_inf(loop_i) % atom) % &
        basis_number + basis_inf(loop_i) % shell - 1) % exponents(loop_k)
        coefficient_i(loop_k) = atom_basis(molecule(basis_inf(loop_i) % atom) &
        % basis_number + basis_inf(loop_i) % shell - 1) &
        % Ncoefficient(loop_k,M_i)
      end do
      coordinate_i = molecule(basis_inf(loop_i) % atom) % nucleus_position
      !-------------------------------
      ! factor of d(x^m)*exp
      if (index(AO_xyz_factor(L_i,M_i),'x') == 4) then
        chardx_i(1) = adjustl('  '//AO_xyz_factor(L_i,M_i)(1:3))
      else if(index(AO_xyz_factor(L_i,M_i),'x') == 3) then
        chardx_i(1) = adjustl('  '//AO_xyz_factor(L_i,M_i)(1:2)//&
        AO_xyz_factor(L_i,M_i)(4:4))
      else if(index(AO_xyz_factor(L_i,M_i),'x') == 2) then
        chardx_i(1) = adjustl('  '//AO_xyz_factor(L_i,M_i)(1:1)//&
        AO_xyz_factor(L_i,M_i)(3:4))
      else if(index(AO_xyz_factor(L_i,M_i),'x') == 1) then
        chardx_i(1) = adjustl('  '//AO_xyz_factor(L_i,M_i)(2:4))
      else
        chardx_i(1) = '     '
      end if
      if (index(AO_xyz_factor(L_i,M_i),'y') == 4) then
        chardy_i(1) = adjustl('  '//AO_xyz_factor(L_i,M_i)(1:3))
      else if(index(AO_xyz_factor(L_i,M_i),'y') == 3) then
        chardy_i(1) = adjustl('  '//AO_xyz_factor(L_i,M_i)(1:2)//&
        AO_xyz_factor(L_i,M_i)(4:4))
      else if(index(AO_xyz_factor(L_i,M_i),'y') == 2) then
        chardy_i(1) = adjustl('  '//AO_xyz_factor(L_i,M_i)(1:1)//&
        AO_xyz_factor(L_i,M_i)(3:4))
      else if(index(AO_xyz_factor(L_i,M_i),'y') == 1) then
        chardy_i(1) = adjustl('  '//AO_xyz_factor(L_i,M_i)(2:4))
      else
        chardy_i(1) = '     '
      end if
      if (index(AO_xyz_factor(L_i,M_i),'z') == 4) then
        chardz_i(1) = adjustl('  '//AO_xyz_factor(L_i,M_i)(1:3))
      else if(index(AO_xyz_factor(L_i,M_i),'z') == 3) then
        chardz_i(1) = adjustl('  '//AO_xyz_factor(L_i,M_i)(1:2)//&
        AO_xyz_factor(L_i,M_i)(4:4))
      else if(index(AO_xyz_factor(L_i,M_i),'z') == 2) then
        chardz_i(1) = adjustl('  '//AO_xyz_factor(L_i,M_i)(1:1)//&
        AO_xyz_factor(L_i,M_i)(3:4))
      else if(index(AO_xyz_factor(L_i,M_i),'z') == 1) then
        chardz_i(1) = adjustl('  '//AO_xyz_factor(L_i,M_i)(2:4))
      else
        chardz_i(1) = '     '
      end if
      !---------------------------------------
      ! factor of x^m*d(exp)
      chardx_i(2) = 'x'//AO_xyz_factor(L_i,M_i)
      chardy_i(2) = 'y'//AO_xyz_factor(L_i,M_i)
      chardz_i(2) = 'z'//AO_xyz_factor(L_i,M_i)
      !---------------------------------------
      ! coefficient of x^m*d(exp)
      loop_k = 2
      do while(loop_k <= 2*contraction_i)
        coedx_i(loop_k) = -2.0_dp * coefficient_i(loop_k/2) * &
        exponents_i(loop_k/2)
        coedy_i(loop_k) = -2.0_dp * coefficient_i(loop_k/2) * &
        exponents_i(loop_k/2)
        coedz_i(loop_k) = -2.0_dp * coefficient_i(loop_k/2) * &
        exponents_i(loop_k/2)
        loop_k = loop_k + 2
      end do
      !---------------------------------
      ! coefficient of d(x^m)*exp
      loop_k = 1
      do while(loop_k < 2*contraction_i)
        coedx_i(loop_k) = 0.0_dp
        coedy_i(loop_k) = 0.0_dp
        coedz_i(loop_k) = 0.0_dp
        do loop_m = 1, len(AO_xyz_factor(L_i,M_i))
          if (AO_xyz_factor(L_i,M_i)(loop_m:loop_m) == 'x') &
          coedx_i(loop_k) = coedx_i(loop_k) + coefficient_i((loop_k + 1) / 2)
          if (AO_xyz_factor(L_i,M_i)(loop_m:loop_m) == 'y') &
          coedy_i(loop_k) = coedy_i(loop_k) + coefficient_i((loop_k + 1) / 2)
          if (AO_xyz_factor(L_i,M_i)(loop_m:loop_m) == 'z') &
          coedz_i(loop_k) = coedz_i(loop_k) + coefficient_i((loop_k + 1) / 2)
        end do
        loop_k = loop_k + 2
      end do
      do loop_j = 1, cbdm
        exponents_j = 0.0_dp
        coefficient_j = 0.0_dp
        coedx_j = 0.0_dp
        coedy_j = 0.0_dp
        coedz_j = 0.0_dp
        contraction_j = atom_basis(molecule(basis_inf(loop_j) % atom) % &
        basis_number + basis_inf(loop_j) % shell - 1) % contraction
        L_j = basis_inf(loop_j) % L
        M_j = basis_inf(loop_j) % M
        do loop_k = 1, contraction_j
          exponents_j(loop_k) = atom_basis(molecule(basis_inf(loop_j) % atom)&
          % basis_number + basis_inf(loop_j) % shell - 1) % exponents(loop_k)
          coefficient_j(loop_k) = &
          atom_basis(molecule(basis_inf(loop_j) % atom) % basis_number + &
          basis_inf(loop_j) % shell - 1) % Ncoefficient(loop_k,M_j)
        end do
        coordinate_j = molecule(basis_inf(loop_j) % atom) % nucleus_position
        if (index(AO_xyz_factor(L_j,M_j),'x') == 4) then
          chardx_j(1) = adjustl('  '//AO_xyz_factor(L_j,M_j)(1:3))
        else if(index(AO_xyz_factor(L_j,M_j),'x') == 3) then
          chardx_j(1) = adjustl('  '//AO_xyz_factor(L_j,M_j)(1:2)//&
          AO_xyz_factor(L_j,M_j)(4:4))
        else if(index(AO_xyz_factor(L_j,M_j),'x') == 2) then
          chardx_j(1) = adjustl('  '//AO_xyz_factor(L_j,M_j)(1:1)//&
          AO_xyz_factor(L_j,M_j)(3:4))
        else if(index(AO_xyz_factor(L_j,M_j),'x') == 1) then
          chardx_j(1) = adjustl('  '//AO_xyz_factor(L_j,M_j)(2:4))
        else
          chardx_j(1) = '     '
        end if
        if (index(AO_xyz_factor(L_j,M_j),'y') == 4) then
          chardy_j(1) = adjustl('  '//AO_xyz_factor(L_j,M_j)(1:3))
        else if(index(AO_xyz_factor(L_j,M_j),'y') == 3) then
          chardy_j(1) = adjustl('  '//AO_xyz_factor(L_j,M_j)(1:2)//&
          AO_xyz_factor(L_j,M_j)(4:4))
        else if(index(AO_xyz_factor(L_j,M_j),'y') == 2) then
          chardy_j(1) = adjustl('  '//AO_xyz_factor(L_j,M_j)(1:1)//&
          AO_xyz_factor(L_j,M_j)(3:4))
        else if(index(AO_xyz_factor(L_j,M_j),'y') == 1) then
          chardy_j(1) = adjustl('  '//AO_xyz_factor(L_j,M_j)(2:4))
        else
          chardy_j(1) = '     '
        end if
        if (index(AO_xyz_factor(L_j,M_j),'z') == 4) then
          chardz_j(1) = adjustl('  '//AO_xyz_factor(L_j,M_j)(1:3))
        else if(index(AO_xyz_factor(L_j,M_j),'z') == 3) then
          chardz_j(1) = adjustl('  '//AO_xyz_factor(L_j,M_j)(1:2)//&
          AO_xyz_factor(L_j,M_j)(4:4))
        else if(index(AO_xyz_factor(L_j,M_j),'z') == 2) then
          chardz_j(1) = adjustl('  '//AO_xyz_factor(L_j,M_j)(1:1)//&
          AO_xyz_factor(L_j,M_j)(3:4))
        else if(index(AO_xyz_factor(L_j,M_j),'z') == 1) then
          chardz_j(1) = adjustl('  '//AO_xyz_factor(L_j,M_j)(2:4))
        else
          chardz_j(1) = '     '
        end if
        chardx_j(2) = 'x'//AO_xyz_factor(L_j,M_j)
        chardy_j(2) = 'y'//AO_xyz_factor(L_j,M_j)
        chardz_j(2) = 'z'//AO_xyz_factor(L_j,M_j)
        loop_k = 2
        do while(loop_k <= 2*contraction_j)
          coedx_j(loop_k) = -2.0_dp * coefficient_j(loop_k/2) * &
          exponents_j(loop_k/2)
          coedy_j(loop_k) = -2.0_dp * coefficient_j(loop_k/2) * &
          exponents_j(loop_k/2)
          coedz_j(loop_k) = -2.0_dp * coefficient_j(loop_k/2) * &
          exponents_j(loop_k/2)
          loop_k = loop_k + 2
        end do
        loop_k = 1
        do while(loop_k < 2*contraction_j)
          coedx_j(loop_k) = 0.0_dp
          coedy_j(loop_k) = 0.0_dp
          coedz_j(loop_k) = 0.0_dp
          do loop_m = 1, len(AO_xyz_factor(L_j,M_j))
            if (AO_xyz_factor(L_j,M_j)(loop_m:loop_m) == 'x') &
            coedx_j(loop_k) = coedx_j(loop_k) + coefficient_j((loop_k + 1) / 2)
            if (AO_xyz_factor(L_j,M_j)(loop_m:loop_m) == 'y') &
            coedy_j(loop_k) = coedy_j(loop_k) + coefficient_j((loop_k + 1) / 2)
            if (AO_xyz_factor(L_j,M_j)(loop_m:loop_m) == 'z') &
            coedz_j(loop_k) = coedz_j(loop_k) + coefficient_j((loop_k + 1) / 2)
          end do
          loop_k = loop_k + 2
        end do
        ! calc <AOi|V|AOj>
        integral_V = calc_1e_V(contraction_i,contraction_j,coefficient_i,&
        coefficient_j,AO_xyz_factor(L_i,M_i),AO_xyz_factor(L_j,M_j),&
        exponents_i,exponents_j,coordinate_i,coordinate_j)
        
        if (DKH_order == 2) then
          ! calc <AOi|pVp|AOj>,totally 9 matrices,6 matrices will be calculated
          integral_DKH2(1) = calc_1e_pVp(contraction_i,contraction_j,&
          coedx_i,coedx_j,chardx_i,chardx_j,exponents_i,exponents_j,&
          coordinate_i,coordinate_j)
          integral_DKH2(2) = calc_1e_pVp(contraction_i,contraction_j,&
          coedy_i,coedy_j,chardy_i,chardy_j,exponents_i,exponents_j,&
          coordinate_i,coordinate_j)
          integral_DKH2(3) = calc_1e_pVp(contraction_i,contraction_j,&
          coedz_i,coedz_j,chardz_i,chardz_j,exponents_i,exponents_j,&
          coordinate_i,coordinate_j)
          integral_DKH2(4) = calc_1e_pVp(contraction_i,contraction_j,&
          coedx_i,coedy_j,chardx_i,chardy_j,exponents_i,exponents_j,&
          coordinate_i,coordinate_j)
          integral_DKH2(5) = integral_DKH2(4)
          integral_DKH2(6) = calc_1e_pVp(contraction_i,contraction_j,&
          coedx_i,coedz_j,chardx_i,chardz_j,exponents_i,exponents_j,&
          coordinate_i,coordinate_j)
          integral_DKH2(7) = integral_DKH2(6)
          integral_DKH2(8) = calc_1e_pVp(contraction_i,contraction_j,&
          coedy_i,coedz_j,chardy_i,chardz_j,exponents_i,exponents_j,&
          coordinate_i,coordinate_j)
          integral_DKH2(9) = integral_DKH2(8)
          ! calc <AOi|p^3Vp|AOj> and <AOi|pVp^3|AOj>, totally 18 matrices,
          ! 9 matrices will be calculated
          if (SRTP_type) then
            integral_SRTP(1) = calc_1e_pppVp('xx',contraction_i,contraction_j,&
            coedx_i,coedx_j,chardx_i,chardx_j,exponents_i,exponents_j,&
            coordinate_i,coordinate_j)
            integral_SRTP(2) = calc_1e_pppVp('yy',contraction_i,contraction_j,&
            coedy_i,coedy_j,chardy_i,chardy_j,exponents_i,exponents_j,&
            coordinate_i,coordinate_j)
            integral_SRTP(3) = calc_1e_pppVp('zz',contraction_i,contraction_j,&
            coedz_i,coedz_j,chardz_i,chardz_j,exponents_i,exponents_j,&
            coordinate_i,coordinate_j)
            integral_SRTP(4) = calc_1e_pppVp('xy',contraction_i,contraction_j,&
            coedx_i,coedy_j,chardx_i,chardy_j,exponents_i,exponents_j,&
            coordinate_i,coordinate_j)
            integral_SRTP(5) = calc_1e_pppVp('yx',contraction_i,contraction_j,&
            coedy_i,coedx_j,chardy_i,chardx_j,exponents_i,exponents_j,&
            coordinate_i,coordinate_j)
            integral_SRTP(6) = calc_1e_pppVp('xz',contraction_i,contraction_j,&
            coedx_i,coedz_j,chardx_i,chardz_j,exponents_i,exponents_j,&
            coordinate_i,coordinate_j)
            integral_SRTP(7) = calc_1e_pppVp('zx',contraction_i,contraction_j,&
            coedz_i,coedx_j,chardz_i,chardx_j,exponents_i,exponents_j,&
            coordinate_i,coordinate_j)
            integral_SRTP(8) = calc_1e_pppVp('yz',contraction_i,contraction_j,&
            coedy_i,coedz_j,chardy_i,chardz_j,exponents_i,exponents_j,&
            coordinate_i,coordinate_j)
            integral_SRTP(9) = calc_1e_pppVp('zy',contraction_i,contraction_j,&
            coedz_i,coedy_j,chardz_i,chardy_j,exponents_i,exponents_j,&
            coordinate_i,coordinate_j)
          end if
        end if
        !-------------<basis overlap integral>-------------
        integral_overlap = 0.0_dp
        do loop_k = 1, contraction_i
          do loop_m = 1, contraction_j
            integral_overlap = integral_overlap + Gaussian_Product_Integral(&
            coefficient_i(loop_k)*coefficient_j(loop_m),AO_xyz_factor(L_i,M_i),&
            AO_xyz_factor(L_j,M_j),exponents_i(loop_k),exponents_j(loop_m),&
            coordinate_i(1),coordinate_i(2),coordinate_i(3),&
            coordinate_j(1),coordinate_j(2),coordinate_j(3))
          end do
        end do
        !-------------<P2 integral>-------------
        integral_dx = 0.0_dp
        do loop_k = 1, 2*contraction_i
          do loop_m = 1, 2*contraction_j
            if (mod(loop_k,2) == 0 .and. mod(loop_m,2) == 0) then
              integral_dx = integral_dx + Gaussian_Product_Integral(&
              coedx_i(loop_k)*coedx_j(loop_m),chardx_i(2),chardx_j(2),&
              exponents_i(loop_k/2),exponents_j(loop_m/2),&
              coordinate_i(1),coordinate_i(2),coordinate_i(3),&
              coordinate_j(1),coordinate_j(2),coordinate_j(3))
            else if (mod(loop_k,2) == 1 .and. mod(loop_m,2) == 0) then
              integral_dx = integral_dx + Gaussian_Product_Integral(&
              coedx_i(loop_k)*coedx_j(loop_m),chardx_i(1),chardx_j(2),&
              exponents_i((loop_k+1)/2),exponents_j(loop_m/2),&
              coordinate_i(1),coordinate_i(2),coordinate_i(3),&
              coordinate_j(1),coordinate_j(2),coordinate_j(3))
            else if (mod(loop_k,2) == 0 .and. mod(loop_m,2) == 1) then
              integral_dx = integral_dx + Gaussian_Product_Integral(&
              coedx_i(loop_k)*coedx_j(loop_m),chardx_i(2),chardx_j(1),&
              exponents_i(loop_k/2),exponents_j((loop_m+1)/2),&
              coordinate_i(1),coordinate_i(2),coordinate_i(3),&
              coordinate_j(1),coordinate_j(2),coordinate_j(3))
            else if (mod(loop_k,2) == 1 .and. mod(loop_m,2) == 1) then
              integral_dx = integral_dx + Gaussian_Product_Integral(&
              coedx_i(loop_k)*coedx_j(loop_m),chardx_i(1),chardx_j(1),&
              exponents_i((loop_k+1)/2),exponents_j((loop_m+1)/2),&
              coordinate_i(1),coordinate_i(2),coordinate_i(3),&
              coordinate_j(1),coordinate_j(2),coordinate_j(3))
            end if
          end do
        end do
        integral_dy = 0.0_dp
        do loop_k = 1, 2*contraction_i
          do loop_m = 1, 2*contraction_j
            if (mod(loop_k,2) == 0 .and. mod(loop_m,2) == 0) then
              integral_dy = integral_dy + Gaussian_Product_Integral(&
              coedy_i(loop_k)*coedy_j(loop_m),chardy_i(2),chardy_j(2),&
              exponents_i(loop_k/2),exponents_j(loop_m/2),&
              coordinate_i(1),coordinate_i(2),coordinate_i(3),&
              coordinate_j(1),coordinate_j(2),coordinate_j(3))
            else if (mod(loop_k,2) == 1 .and. mod(loop_m,2) == 0) then
              integral_dy = integral_dy + Gaussian_Product_Integral(&
              coedy_i(loop_k)*coedy_j(loop_m),chardy_i(1),chardy_j(2),&
              exponents_i((loop_k+1)/2),exponents_j(loop_m/2),&
              coordinate_i(1),coordinate_i(2),coordinate_i(3),&
              coordinate_j(1),coordinate_j(2),coordinate_j(3))
            else if (mod(loop_k,2) == 0 .and. mod(loop_m,2) == 1) then
              integral_dy = integral_dy + Gaussian_Product_Integral(&
              coedy_i(loop_k)*coedy_j(loop_m),chardy_i(2),chardy_j(1),&
              exponents_i(loop_k/2),exponents_j((loop_m+1)/2),&
              coordinate_i(1),coordinate_i(2),coordinate_i(3),&
              coordinate_j(1),coordinate_j(2),coordinate_j(3))
            else if (mod(loop_k,2) == 1 .and. mod(loop_m,2) == 1) then
              integral_dy = integral_dy + Gaussian_Product_Integral(&
              coedy_i(loop_k)*coedy_j(loop_m),chardy_i(1),chardy_j(1),&
              exponents_i((loop_k+1)/2),exponents_j((loop_m+1)/2),&
              coordinate_i(1),coordinate_i(2),coordinate_i(3),&
              coordinate_j(1),coordinate_j(2),coordinate_j(3))
            end if
          end do
        end do
        integral_dz = 0.0_dp
        do loop_k = 1, 2*contraction_i
          do loop_m = 1, 2*contraction_j
            if (mod(loop_k,2) == 0 .and. mod(loop_m,2) == 0) then
              integral_dz = integral_dz + Gaussian_Product_Integral(&
              coedz_i(loop_k)*coedz_j(loop_m),chardz_i(2),chardz_j(2),&
              exponents_i(loop_k/2),exponents_j(loop_m/2),&
              coordinate_i(1),coordinate_i(2),coordinate_i(3),&
              coordinate_j(1),coordinate_j(2),coordinate_j(3))
            else if (mod(loop_k,2) == 1 .and. mod(loop_m,2) == 0) then
              integral_dz = integral_dz + Gaussian_Product_Integral(&
              coedz_i(loop_k)*coedz_j(loop_m),chardz_i(1),chardz_j(2),&
              exponents_i((loop_k+1)/2),exponents_j(loop_m/2),&
              coordinate_i(1),coordinate_i(2),coordinate_i(3),&
              coordinate_j(1),coordinate_j(2),coordinate_j(3))
            else if (mod(loop_k,2) == 0 .and. mod(loop_m,2) == 1) then
              integral_dz = integral_dz + Gaussian_Product_Integral(&
              coedz_i(loop_k)*coedz_j(loop_m),chardz_i(2),chardz_j(1),&
              exponents_i(loop_k/2),exponents_j((loop_m+1)/2),&
              coordinate_i(1),coordinate_i(2),coordinate_i(3),&
              coordinate_j(1),coordinate_j(2),coordinate_j(3))
            else if (mod(loop_k,2) == 1 .and. mod(loop_m,2) == 1) then
              integral_dz = integral_dz + Gaussian_Product_Integral(&
              coedz_i(loop_k)*coedz_j(loop_m),chardz_i(1),chardz_j(1),&
              exponents_i((loop_k+1)/2),exponents_j((loop_m+1)/2),&
              coordinate_i(1),coordinate_i(2),coordinate_i(3),&
              coordinate_j(1),coordinate_j(2),coordinate_j(3))
            end if
          end do
        end do
        !$omp atomic
        i_j(loop_i,loop_j) = i_j(loop_i,loop_j) + integral_overlap
        !$omp atomic
        i_p2_j(loop_i,loop_j) = i_p2_j(loop_i,loop_j) + &
        integral_dx + integral_dy + integral_dz
        !$omp atomic
        i_V_j(loop_i,loop_j) = i_V_j(loop_i,loop_j) + integral_V
        if (DKH_order == 2) then
          !$omp atomic
          i_pxVpx_j(loop_i,loop_j) = i_pxVpx_j(loop_i,loop_j) + integral_DKH2(1)
          !$omp atomic
          i_pyVpy_j(loop_i,loop_j) = i_pyVpy_j(loop_i,loop_j) + integral_DKH2(2)
          !$omp atomic
          i_pzVpz_j(loop_i,loop_j) = i_pzVpz_j(loop_i,loop_j) + integral_DKH2(3)
          !$omp atomic
          i_pxVpy_j(loop_i,loop_j) = i_pxVpy_j(loop_i,loop_j) + integral_DKH2(4)
          !$omp atomic
          i_pyVpx_j(loop_j,loop_i) = i_pyVpx_j(loop_j,loop_i) + integral_DKH2(5)
          !$omp atomic
          i_pxVpz_j(loop_i,loop_j) = i_pxVpz_j(loop_i,loop_j) + integral_DKH2(6)
          !$omp atomic
          i_pzVpx_j(loop_j,loop_i) = i_pzVpx_j(loop_j,loop_i) + integral_DKH2(7)
          !$omp atomic
          i_pyVpz_j(loop_i,loop_j) = i_pyVpz_j(loop_i,loop_j) + integral_DKH2(8)
          !$omp atomic
          i_pzVpy_j(loop_j,loop_i) = i_pzVpy_j(loop_j,loop_i) + integral_DKH2(9)
          if (SRTP_type) then
            !$omp atomic
            i_px3Vpx_j(loop_i,loop_j) = &
            i_px3Vpx_j(loop_i,loop_j) + integral_SRTP(1)
            !$omp atomic
            i_py3Vpy_j(loop_i,loop_j) = &
            i_py3Vpy_j(loop_i,loop_j) + integral_SRTP(2)
            !$omp atomic
            i_pz3Vpz_j(loop_i,loop_j) = &
            i_pz3Vpz_j(loop_i,loop_j) + integral_SRTP(3)
            !$omp atomic
            i_px3Vpy_j(loop_i,loop_j) = &
            i_px3Vpy_j(loop_i,loop_j) + integral_SRTP(4)
            !$omp atomic
            i_py3Vpx_j(loop_i,loop_j) = &
            i_py3Vpx_j(loop_i,loop_j) + integral_SRTP(5)
            !$omp atomic
            i_px3Vpz_j(loop_i,loop_j) = &
            i_px3Vpz_j(loop_i,loop_j) + integral_SRTP(6)
            !$omp atomic
            i_pz3Vpx_j(loop_i,loop_j) = &
            i_pz3Vpx_j(loop_i,loop_j) + integral_SRTP(7)
            !$omp atomic
            i_py3Vpz_j(loop_i,loop_j) = &
            i_py3Vpz_j(loop_i,loop_j) + integral_SRTP(8)
            !$omp atomic
            i_pz3Vpy_j(loop_i,loop_j) = &
            i_pz3Vpy_j(loop_i,loop_j) + integral_SRTP(9)
          end if
        end if
      end do
    end do
    !$omp end do
    !$omp end parallel
    ! consider the relation p = -iD, 
    ! p3Vp and pVp3 change sign, p2 and pVp nochange
    if (SRTP_type) then
      i_px3Vpx_j = - i_px3Vpx_j
      i_py3Vpy_j = - i_py3Vpy_j
      i_pz3Vpz_j = - i_pz3Vpz_j
      i_px3Vpy_j = - i_px3Vpy_j
      i_py3Vpx_j = - i_py3Vpx_j
      i_px3Vpz_j = - i_px3Vpz_j
      i_pz3Vpx_j = - i_pz3Vpx_j
      i_py3Vpz_j = - i_py3Vpz_j
      i_pz3Vpy_j = - i_pz3Vpy_j
    end if
  end subroutine assign_matrices_1e

!-----------------------------------------------------------------------
!> calc values to matrix <AOi|V|AOj>
  pure function calc_1e_V(&
  contr_i,contr_j,coe_i,coe_j,fac_i,fac_j,expo_i,expo_j,cod_i,cod_j) &
  result(integral_pot)
    implicit none
    integer,intent(in) :: contr_i
    integer,intent(in) :: contr_j
    integer :: bloop_i,bloop_j,bloop_pot   ! loop variables only for calc_1e_V
    real(dp),intent(in) :: coe_i(contr_i)  ! coefficients of i^th orbital
    real(dp),intent(in) :: coe_j(contr_j)  ! coefficients of j^th orbital
    real(dp),intent(in) :: expo_i(contr_i) ! exponents of i^th orbital
    real(dp),intent(in) :: expo_j(contr_j) ! exponents of j^th orbital
    real(dp),intent(in) :: cod_i(3)        ! centrol coordintates.i^th orbital
    real(dp),intent(in) :: cod_j(3)        ! centrol coordintates.j^th orbital
    real(dp) :: bcod_i(3)                  ! extended coordintates.i^th orbital
    real(dp) :: bcod_j(3)                  ! extended coordintates.j^th orbital
    real(dp) :: Z_pot, R_pot
    real(dp) :: integral_pot
    real(dp) :: cod_pot(3)
    character(len = *),intent(in) :: fac_i ! xyz factor of i^th orbital
    character(len = *),intent(in) :: fac_j ! xyz factor of j^th orbital
    integral_pot = 0.0_dp
    do bloop_pot = 1, atom_count
      cod_pot = molecule(bloop_pot) % nucleus_position
      Z_pot = real(molecule(bloop_pot) % atom_number)
      R_pot = molecule(bloop_pot) % nucleus_radius / fm2Bohr
      ! center of potential atom set to zero
      bcod_i(1) = cod_i(1) - cod_pot(1)
      bcod_i(2) = cod_i(2) - cod_pot(2)
      bcod_i(3) = cod_i(3) - cod_pot(3)
      bcod_j(1) = cod_j(1) - cod_pot(1)
      bcod_j(2) = cod_j(2) - cod_pot(2)
      bcod_j(3) = cod_j(3) - cod_pot(3)
      do bloop_i = 1, contr_i
        do bloop_j = 1, contr_j
          integral_pot = integral_pot + V_Integral_1e(&
          Z_pot,coe_i(bloop_i),coe_j(bloop_j),fac_i,fac_j,&
          expo_i(bloop_i),expo_j(bloop_j),bcod_i(1),bcod_i(2),bcod_i(3),&
          bcod_j(1),bcod_j(2),bcod_j(3),R_pot)
        end do
      end do
    end do
    return 
  end function calc_1e_V
  
!-----------------------------------------------------------------------
!> calc values to matrix <AOi|pVp|AOj> (9 matrices)
!!
!! pxVpx pyVpy pzVpz pxVpy pyVpx pxVpz pzVpx pyVpz pzVpy
  pure function calc_1e_pVp(&
  contr_i,contr_j,coe_i,coe_j,fac_i,fac_j,expo_i,expo_j,cod_i,cod_j) &
  result(integral_pot)
    implicit none
    integer,intent(in) :: contr_i
    integer,intent(in) :: contr_j
    integer :: sloop_i,sloop_j,sloop_pot    ! loop variables for calc_1e_pVp
    ! coefficients of first-order derivative of i^th orbital
    real(dp),intent(in) :: coe_i(2*contr_i)
    ! coefficients of first-order derivative of j^th orbital
    real(dp),intent(in) :: coe_j(2*contr_j)
    real(dp),intent(in) :: expo_i(contr_i)  ! exponents of i^th orbital
    real(dp),intent(in) :: expo_j(contr_j)  ! exponents of j^th orbital
    real(dp),intent(in) :: cod_i(3)         ! centrol coordintates.i^th orbital
    real(dp),intent(in) :: cod_j(3)         ! centrol coordintates.j^th orbital
    real(dp) :: scod_i(3)                   ! extended coordintates.i^th orbital
    real(dp) :: scod_j(3)                   ! extended coordintates.j^th orbital
    real(dp) :: Z_pot, R_pot
    real(dp) :: integral_pot
    real(dp) :: cod_pot(3)
    character(len = *),intent(in) :: fac_i(2)  ! xyz factor of i^th orbital
    character(len = *),intent(in) :: fac_j(2)  ! xyz factor of j^th orbital
    integral_pot = 0.0_dp
    do sloop_pot = 1, atom_count
      cod_pot = molecule(sloop_pot) % nucleus_position
      Z_pot = real(molecule(sloop_pot) % atom_number)
      R_pot = molecule(sloop_pot) % nucleus_radius / fm2Bohr
      ! center of potential atom set to zero
      scod_i(1) = cod_i(1) - cod_pot(1)
      scod_i(2) = cod_i(2) - cod_pot(2)
      scod_i(3) = cod_i(3) - cod_pot(3)
      scod_j(1) = cod_j(1) - cod_pot(1)
      scod_j(2) = cod_j(2) - cod_pot(2)
      scod_j(3) = cod_j(3) - cod_pot(3)
      do sloop_i = 1, 2*contr_i
        do sloop_j = 1, 2*contr_j
          integral_pot = integral_pot + V_Integral_1e(&
          Z_pot,coe_i(sloop_i),coe_j(sloop_j),fac_i(2-mod(sloop_i,2)),&
          fac_j(2-mod(sloop_j,2)),expo_i((sloop_i+mod(sloop_i,2))/2),&
          expo_j((sloop_j+mod(sloop_j,2))/2),scod_i(1),scod_i(2),scod_i(3),&
          scod_j(1),scod_j(2),scod_j(3),R_pot)
        end do
      end do
    end do
    return
  end function calc_1e_pVp
  
!-----------------------------------------------------------------------
!> calc values to matrix <AOi|p^3Vp|AOj> (9 matrices)
!!
!! px3Vpx py3Vpy pz3Vpz px3Vpy py3Vpx px3Vpz pz3Vpx py3Vpz pz3Vpy
  pure function calc_1e_pppVp(&
  char_der,contr_i,contr_j,coe_i,coe_j,fac_i,fac_j,expo_i,expo_j,cod_i,cod_j) &
  result(integral_pot)
    implicit none
    integer,intent(in) :: contr_i
    integer,intent(in) :: contr_j
    integer :: tloop_i,tloop_j,tloop_pot    ! loop variables for calc_1e_pppVp
    real(dp),intent(in) :: coe_i(2*contr_i) ! coefficients of i^th orbital
    real(dp),intent(in) :: coe_j(2*contr_j) ! coefficients of j^th orbital
    ! coefficients of i^th orbital after 2 derivative actions
    real(dp) :: mcoe_i(4*contr_i)
    ! coefficients of i^th orbital after 3 derivative actions
    real(dp) :: tcoe_i(8*contr_i)
    real(dp),intent(in) :: expo_i(contr_i)  ! exponents of i^th orbital
    real(dp),intent(in) :: expo_j(contr_j)  ! exponents of j^th orbital
    real(dp),intent(in) :: cod_i(3)         ! centrol coordintates.i^th orbital
    real(dp),intent(in) :: cod_j(3)         ! centrol coordintates.j^th orbital
    real(dp) :: tcod_i(3) ! coordintates.i^th orbital after 3 derivative actions
    real(dp) :: tcod_j(3) ! coordintates.j^th orbital after 3 derivative actions
    real(dp) :: Z_pot, R_pot
    real(dp) :: integral_pot
    real(dp) :: cod_pot(3)
    ! which pVp matrix is assignning (e.g. yy for py3Vpy)
    character(len = 2),intent(in) :: char_der
    character(len = *),intent(in) :: fac_i(2) ! xyz factor of i^th orbital
    ! xyz factor.i^th orbital after 2 derivative actions
    character(len = 8) :: mfac_i(4)
    ! xyz factor.i^th orbital after 3 derivative actions
    character(len = 16) :: tfac_i(8)
    character(len = *),intent(in) :: fac_j(2) ! xyz factor of j^th orbital
    ! find the derivative of p|AOi> respect to x/y/z
    character(len = 1) :: ele_der
    if (char_der == 'xx' .or. char_der == 'xy' .or. char_der == 'xz') then
      ele_der = 'x'
    else if (char_der == 'yy' .or. char_der == 'yx' .or. char_der == 'yz') then
      ele_der = 'y'
    else if (char_der == 'zz' .or. char_der == 'zx' .or. char_der == 'zy') then
      ele_der = 'z'
    end if
    ! derivative of xyz factor of fac_i(1,2)
    do tloop_i = 1, 2
      do tloop_j = 1, len(fac_i(tloop_i))
        if (index(fac_i(tloop_i),ele_der) == tloop_j) then
          if (tloop_j == 0) then
            mfac_i(2*tloop_i-1) = '  '
          else if (tloop_j == 1) then
            mfac_i(2*tloop_i-1) = adjustl('  '//&
            fac_i(tloop_i)(2:len(fac_i(tloop_i))))
          else if (tloop_j == len(fac_i(tloop_i))) then
            mfac_i(2*tloop_i-1) = '  '//&
            fac_i(tloop_i)(1:len(fac_i(tloop_i))-1)
          else
            mfac_i(2*tloop_i-1) = adjustl('  '//fac_i(tloop_i)(1:tloop_j-1)//&
            fac_i(tloop_i)(tloop_j+1:len(fac_i(tloop_i))))
          end if
          exit
        end if
      end do
    end do
    ! derivative of exponent of fac_i(1,2)
    mfac_i(2) = ele_der//fac_i(1)
    mfac_i(4) = ele_der//fac_i(2)
    ! derivative of xyz factor of mfac_i(1,2,3,4)
    do tloop_i = 1, 4
      do tloop_j = 1, len(mfac_i(tloop_i))
        if (index(mfac_i(tloop_i),ele_der) == tloop_j) then
          if (tloop_j == 0) then
            tfac_i(2*tloop_i-1) = '  '
          else if (tloop_j == 1) then
            tfac_i(2*tloop_i-1) = adjustl('  '//&
            mfac_i(tloop_i)(2:len(mfac_i(tloop_i))))
          else if (tloop_j == len(mfac_i(tloop_i))) then
            tfac_i(2*tloop_i-1) = '  '//&
            mfac_i(tloop_i)(1:len(mfac_i(tloop_i))-1)
          else
            tfac_i(2*tloop_i-1) = adjustl('  '//mfac_i(tloop_i)(1:tloop_j-1)//&
            mfac_i(tloop_i)(tloop_j+1:len(mfac_i(tloop_i))))
          end if
          exit
        end if
      end do
    end do
    ! derivative of exponent of fac_i(1,2,3,4)
    tfac_i(2) = ele_der//mfac_i(1)
    tfac_i(4) = ele_der//mfac_i(2)
    tfac_i(6) = ele_der//mfac_i(3)
    tfac_i(8) = ele_der//mfac_i(4)
    ! derivative of xyz factor of coe_i
    do tloop_i = 1, 2*contr_i
      mcoe_i(2*tloop_i-1) = 0.0_dp
      do tloop_j = 1, len(fac_i(2-mod(tloop_i,2)))
        if (index(fac_i(2-mod(tloop_i,2))(tloop_j:tloop_j),ele_der) /= 0) &
        mcoe_i(2*tloop_i-1) = mcoe_i(2*tloop_i-1) + coe_i(tloop_i)
      end do
    end do
    ! derivative of exponent of coe_i
    do tloop_i = 1, 2*contr_i
      mcoe_i(2*tloop_i) = -2.0_dp * coe_i(tloop_i) * &
      expo_i((tloop_i+mod(tloop_i,2))/2)
    end do
    ! derivative of xyz factor of mcoe_i
    do tloop_i = 1, 4*contr_i
      tcoe_i(2*tloop_i-1) = 0.0_dp
      do tloop_j = 1, len(mfac_i(4-mod(4-mod(tloop_i,4),4)))
        if (index(mfac_i(4-mod(4-mod(tloop_i,4),4))(tloop_j:tloop_j),ele_der) &
        /= 0) tcoe_i(2*tloop_i-1) = tcoe_i(2*tloop_i-1) + mcoe_i(tloop_i)
      end do
    end do
    ! derivative of exponent of mcoe_i
    do tloop_i = 1, 4*contr_i
      tcoe_i(2*tloop_i) = -2.0_dp * mcoe_i(tloop_i) * &
      expo_i((tloop_i+mod(4-mod(tloop_i,4),4))/4)
    end do
    integral_pot = 0.0_dp
    do tloop_pot = 1, atom_count
      cod_pot = molecule(tloop_pot) % nucleus_position
      Z_pot = real(molecule(tloop_pot) % atom_number)
      R_pot = molecule(tloop_pot) % nucleus_radius / fm2Bohr
      ! center of potential atom set to zero
      tcod_i(1) = cod_i(1) - cod_pot(1)
      tcod_i(2) = cod_i(2) - cod_pot(2)
      tcod_i(3) = cod_i(3) - cod_pot(3)
      tcod_j(1) = cod_j(1) - cod_pot(1)
      tcod_j(2) = cod_j(2) - cod_pot(2)
      tcod_j(3) = cod_j(3) - cod_pot(3)
      do tloop_i = 1, 8*contr_i
        do tloop_j = 1, 2*contr_j
          integral_pot = integral_pot + &
          V_Integral_1e(Z_pot,tcoe_i(tloop_i),coe_j(tloop_j),&
          tfac_i(8-mod(8-mod(tloop_i,8),8)),fac_j(2-mod(tloop_j,2)),&
          expo_i((tloop_i+mod(8-mod(tloop_i,8),8))/8),&
          expo_j((tloop_j+mod(tloop_j,2))/2),tcod_i(1),tcod_i(2),tcod_i(3),&
          tcod_j(1),tcod_j(2),tcod_j(3),R_pot)
        end do
      end do
    end do
    return
  end function calc_1e_pppVp
  
  
!-----------------------------------------------------------------------
!> full space integration of product of 2 Gaussian functions in Cartesian: 
  pure elemental real(dp) function Gaussian_Product_Integral(&
  coef,fac_i,fac_j,expo_i,expo_j,cod_ix,cod_iy,cod_iz,&
  cod_jx,cod_jy,cod_jz)
    implicit none
    real(dp),intent(in) :: coef       ! coefficient before Gaussian product
    real(dp) :: coeff                 ! coefficient after Gaussian product
    real(dp) :: expo                  ! exponent of Gaussian product
    real(dp) :: cod(3)                ! coordinate of Gaussian product
    real(dp),intent(in) :: expo_i
    real(dp),intent(in) :: expo_j
    real(dp),intent(in) :: cod_ix, cod_iy, cod_iz
    real(dp),intent(in) :: cod_jx, cod_jy, cod_jz
    real(dp) :: cod_i_pos(3)          ! linear transformation of cod_i
    real(dp) :: integral_x
    real(dp) :: integral_x_mic,integral_x_mic_pre,integral_x_mic_pre_pre
    real(dp) :: integral_y
    real(dp) :: integral_y_mic,integral_y_mic_pre,integral_y_mic_pre_pre
    real(dp) :: integral_z
    real(dp) :: integral_z_mic,integral_z_mic_pre,integral_z_mic_pre_pre
    integer :: m_x_i, m_y_i, m_z_i, m_x_j, m_y_j, m_z_j ! largest power of x,y,z
    integer :: gloop_i,gloop_j        ! loop variables only for GPI function
    character(len = *),intent(in) :: fac_i
    character(len = *),intent(in) :: fac_j
    m_x_i = 0
    m_y_i = 0
    m_z_i = 0
    m_x_j = 0
    m_y_j = 0
    m_z_j = 0
    integral_x = 0.0
    integral_y = 0.0
    integral_z = 0.0
    ! Gaussian function produntion
    coeff = coef * exp(-((cod_ix-cod_jx)**2) / (1.0_dp/expo_i+1.0_dp/expo_j))
    coeff = coeff * exp(-((cod_iy-cod_jy)**2) / (1.0_dp/expo_i+1.0_dp/expo_j))
    coeff = coeff * exp(-((cod_iz-cod_jz)**2) / (1.0_dp/expo_i+1.0_dp/expo_j))
    expo = expo_i + expo_j
    cod(1) = (cod_ix * expo_i + cod_jx * expo_j)/(expo_i + expo_j)
    cod(2) = (cod_iy * expo_i + cod_jy * expo_j)/(expo_i + expo_j)
    cod(3) = (cod_iz * expo_i + cod_jz * expo_j)/(expo_i + expo_j)
    do gloop_i = 1, len(fac_i)
      if(fac_i(gloop_i:gloop_i) == 'x') m_x_i = m_x_i + 1
      if(fac_i(gloop_i:gloop_i) == 'y') m_y_i = m_y_i + 1
      if(fac_i(gloop_i:gloop_i) == 'z') m_z_i = m_z_i + 1
    end do
    do gloop_i = 1, len(fac_j)
      if(fac_j(gloop_i:gloop_i) == 'x') m_x_j = m_x_j + 1
      if(fac_j(gloop_i:gloop_i) == 'y') m_y_j = m_y_j + 1
      if(fac_j(gloop_i:gloop_i) == 'z') m_z_j = m_z_j + 1
    end do
    ! linear transformation of integral variables
    cod_i_pos(1) = cod_ix - cod_jx
    cod_i_pos(2) = cod_iy - cod_jy
    cod_i_pos(3) = cod_iz - cod_jz
    cod(1) = cod(1) - cod_jx
    cod(2) = cod(2) - cod_jy
    cod(3) = cod(3) - cod_jz
    ! binomial expansion -x
    do gloop_i = 0, m_x_i
      do gloop_j = 0, m_x_i - gloop_i + m_x_j   ! integral x^m*exp(-b*(x-x0)^2)
        if (gloop_j == 0) then
          integral_x_mic = dsqrt(pi/expo)
        else if(gloop_j == 1) then
          integral_x_mic_pre = integral_x_mic
          integral_x_mic = cod(1) * integral_x_mic_pre
        else
          integral_x_mic_pre_pre = integral_x_mic_pre
          integral_x_mic_pre = integral_x_mic
          integral_x_mic = cod(1) * integral_x_mic_pre + &
          (real(gloop_j-1)*integral_x_mic_pre_pre)/(2.0_dp*expo)
        end if
      end do
      integral_x = integral_x + binomialcoe(m_x_i,gloop_i) * &
      (-cod_i_pos(1))**(gloop_i) * integral_x_mic
    end do
    ! binomial expansion -y
    do gloop_i = 0, m_y_i
      do gloop_j = 0, m_y_i - gloop_i + m_y_j
        if (gloop_j == 0) then
          integral_y_mic = dsqrt(pi/expo)
        else if(gloop_j == 1) then
          integral_y_mic_pre = integral_y_mic
          integral_y_mic = cod(2) * integral_y_mic_pre
        else
          integral_y_mic_pre_pre = integral_y_mic_pre
          integral_y_mic_pre = integral_y_mic
          integral_y_mic = cod(2) * integral_y_mic_pre + &
          (real(gloop_j-1)*integral_y_mic_pre_pre)/(2.0_dp*expo)
        end if
      end do
      integral_y = integral_y + binomialcoe(m_y_i,gloop_i) * &
      (-cod_i_pos(2))**(gloop_i) * integral_y_mic
    end do
    ! binomial expansion -z
    do gloop_i = 0, m_z_i
      do gloop_j = 0, m_z_i - gloop_i + m_z_j
        if (gloop_j == 0) then
          integral_z_mic = dsqrt(pi/expo)
        else if(gloop_j == 1) then
          integral_z_mic_pre = integral_z_mic
          integral_z_mic = cod(3) * integral_z_mic_pre
        else
          integral_z_mic_pre_pre = integral_z_mic_pre
          integral_z_mic_pre = integral_z_mic
          integral_z_mic = cod(3) * integral_z_mic_pre + &
          (real(gloop_j-1)*integral_z_mic_pre_pre)/(2.0_dp*expo)
        end if
      end do
      integral_z = integral_z + binomialcoe(m_z_i,gloop_i) * &
      (-cod_i_pos(3))**(gloop_i) * integral_z_mic
    end do
    Gaussian_Product_Integral = coeff * integral_x * integral_y * integral_z
    return
  end function Gaussian_Product_Integral
  
!-----------------------------------------------------------------------
!> integration of electron-nuclear attraction potential in Cartesian coordinate:
!!
!! inategral transformation: (x-xi)^m*(x-xj)^n*expo(-b*x^2)*expo(x^2*t^2)dxdt
  pure elemental real(dp) function V_Integral_1e(&
  Z,coe_i,coe_j,fac_i,fac_j,expo_i,expo_j,cod_ix,cod_iy,cod_iz,cod_jx,cod_jy,&
  cod_jz,rn)
    implicit none
    real(dp),intent(in) :: expo_i, expo_j
    real(dp) :: expo
    real(dp),intent(in) :: coe_i, coe_j
    real(dp) :: coe
    real(dp),intent(in) :: Z
    ! finite nuclear potential will not be considered if rn = 0.
    real(dp),intent(in) :: rn
    real(dp),intent(in) :: cod_ix, cod_iy, cod_iz   ! coordination of i
    real(dp),intent(in) :: cod_jx, cod_jy, cod_jz   ! coordination of j
    real(dp) :: cod(3)
    real(dp) :: R2
    ! t-containing coefficients, t2pb_x(1) = coeff.(t^2+b)^(1/2).x integration
    real(dp) :: t2pb_x(10)
    ! t-containing coefficients, t2pb_y(1) = coeff.(t^2+b)^(1/2).y integration
    real(dp) :: t2pb_y(10)
    ! t-containing coefficients, t2pb_z(1) = coeff.(t^2+b)^(1/2).z integration
    real(dp) :: t2pb_z(10)
    ! t-containing coefficients, t2pb(1) = coeff.(t^2+b)^(1/2).x,y,z integration
    real(dp) :: t2pb(14)
    real(dp) :: t2pb_mic(10)
    real(dp) :: t2pb_mic_pre(10)
    real(dp) :: t2pb_mic_pre_pre(10)
    real(dp) :: integral, integral_mic, integral_mic_pre, integral_mic_pre_pre
    real(dp) :: NX_mic(30), NX_mic_pre(30), NX_mic_pre_pre(30)
    integer :: tayeps  ! number of Taylor expansion series of integration at X=0
    ! direct integration (X > xts); Taylor expansion integration (X <= xts)
    real(dp) :: xts
    integer :: m_x_i, m_y_i, m_z_i, m_x_j, m_y_j, m_z_j ! largest power of x,y,z
    integer :: vloop_i,vloop_j,vloop_k,vloop_o ! loop variable for V_Integral_1e
    integer :: vloop_mic,vloop_mic_pre
    character(len = *),intent(in) :: fac_i
    character(len = *),intent(in) :: fac_j
    m_x_i = 0
    m_y_i = 0
    m_z_i = 0
    m_x_j = 0
    m_y_j = 0
    m_z_j = 0
    do vloop_i = 1, len(fac_i)
      if(fac_i(vloop_i:vloop_i) == 'x') m_x_i = m_x_i + 1
      if(fac_i(vloop_i:vloop_i) == 'y') m_y_i = m_y_i + 1
      if(fac_i(vloop_i:vloop_i) == 'z') m_z_i = m_z_i + 1
    end do
    do vloop_i = 1, len(fac_j)
      if(fac_j(vloop_i:vloop_i) == 'x') m_x_j = m_x_j + 1
      if(fac_j(vloop_i:vloop_i) == 'y') m_y_j = m_y_j + 1
      if(fac_j(vloop_i:vloop_i) == 'z') m_z_j = m_z_j + 1
    end do
    if (2*(m_x_i + m_x_j + m_y_i + m_y_j + m_z_i + m_z_j + 4) <= 10) then
      tayeps = 5
      xts = 0.01
    else if (2*(m_x_i + m_x_j + m_y_i + m_y_j + m_z_i + m_z_j + 4) <= 15) then
      tayeps = 8
      xts = 0.1
    else if (2*(m_x_i + m_x_j + m_y_i + m_y_j + m_z_i + m_z_j + 4) <= 25) then
      tayeps = 15
      xts = 1.0
    else if (2*(m_x_i + m_x_j + m_y_i + m_y_j + m_z_i + m_z_j + 4) <= 35) then
      tayeps = 20
      xts = 1.0
    else if (2*(m_x_i + m_x_j + m_y_i + m_y_j + m_z_i + m_z_j + 4) <= 40) then
      tayeps = 25
      xts = 3.0
    else if (2*(m_x_i + m_x_j + m_y_i + m_y_j + m_z_i + m_z_j + 4) <= 50) then
      tayeps = 30
      xts = 4.0
    end if
    NX_mic = 0.0_dp
    NX_mic_pre = 0.0_dp
    NX_mic_pre_pre = 0.0_dp
    
    !--------------------------
    ! Gaussian function production
    coe = -Z * coe_i * coe_j
    coe = coe * exp(-((cod_ix - cod_jx)**2.0_dp)&
    /(1.0_dp/expo_i + 1.0_dp/expo_j))
    coe = coe * exp(-((cod_iy - cod_jy)**2.0_dp)&
    /(1.0_dp/expo_i + 1.0_dp/expo_j))
    coe = coe * exp(-((cod_iz - cod_jz)**2.0_dp)&
    /(1.0_dp/expo_i + 1.0_dp/expo_j))
    !------
    expo = expo_i + expo_j
    !------
    cod(1) = (cod_ix * expo_i + cod_jx * expo_j)/(expo_i + expo_j)
    cod(2) = (cod_iy * expo_i + cod_jy * expo_j)/(expo_i + expo_j)
    cod(3) = (cod_iz * expo_i + cod_jz * expo_j)/(expo_i + expo_j)
    R2 = cod(1)*cod(1) + cod(2)*cod(2) + cod(3)*cod(3)
    ! use binomial expansion for easy storage of t-containing coefficients.
    !--------------------------
    ! integral of x, generate a coefficient exp(-b*((cod(1))^2*t^2)/(t^2+b))
    t2pb_x = 0.0_dp
    vloop_i = 0
    do while(vloop_i <= m_x_i)
      vloop_j = 0
      do while(vloop_j <= m_x_j)
        t2pb_mic = 0.0_dp
        t2pb_mic_pre = 0.0_dp
        t2pb_mic_pre_pre = 0.0_dp
        do vloop_mic = 0, m_x_i + m_x_j - vloop_i - vloop_j
          if (vloop_mic == 0) then
            t2pb_mic(1) = dsqrt(pi)
          else if(vloop_mic == 1) then
            t2pb_mic_pre(1) = dsqrt(pi)
            t2pb_mic = 0.0_dp
            t2pb_mic(2) = cod(1) * expo * t2pb_mic_pre(1)
          else
            t2pb_mic_pre_pre = t2pb_mic_pre
            t2pb_mic_pre = t2pb_mic
            t2pb_mic = 0.0_dp
            do vloop_mic_pre = 0, vloop_mic
              t2pb_mic(vloop_mic_pre + 2) = t2pb_mic(vloop_mic_pre + 2) + &
              0.5_dp * real(vloop_mic-1) * t2pb_mic_pre_pre(vloop_mic_pre+1)
              t2pb_mic(vloop_mic_pre + 2) = t2pb_mic(vloop_mic_pre + 2) + &
              cod(1) * expo * t2pb_mic_pre(vloop_mic_pre+1)
            end do
          end if
        end do
        t2pb_x = t2pb_x + &
        binomialcoe(m_x_i,vloop_i) * binomialcoe(m_x_j,vloop_j) * &
        (-cod_ix)**(vloop_i) * (-cod_jx)**(vloop_j) * t2pb_mic
        vloop_j = vloop_j + 1
      end do
      vloop_i = vloop_i + 1 
    end do
    !--------------------------
    ! integral of x, generate a coefficient exp(-b*((cod(2))^2*t^2)/(t^2+b))
    t2pb_y = 0.0_dp
    vloop_i = 0
    do while(vloop_i <= m_y_i)
      vloop_j = 0
      do while(vloop_j <= m_y_j)
        t2pb_mic = 0.0_dp
        t2pb_mic_pre = 0.0_dp
        t2pb_mic_pre_pre = 0.0_dp
        do vloop_mic = 0, m_y_i + m_y_j - vloop_i - vloop_j
          if (vloop_mic == 0) then
            t2pb_mic(1) = dsqrt(pi)
          else if(vloop_mic == 1) then
            t2pb_mic_pre(1) = dsqrt(pi)
            t2pb_mic = 0.0_dp
            t2pb_mic(2) = cod(2) * expo * t2pb_mic_pre(1)
          else
            t2pb_mic_pre_pre = t2pb_mic_pre
            t2pb_mic_pre = t2pb_mic
            t2pb_mic = 0.0_dp
            do vloop_mic_pre = 0, vloop_mic
              t2pb_mic(vloop_mic_pre + 2) = t2pb_mic(vloop_mic_pre + 2) + &
              0.5_dp * real(vloop_mic-1) * t2pb_mic_pre_pre(vloop_mic_pre+1)
              t2pb_mic(vloop_mic_pre + 2) = t2pb_mic(vloop_mic_pre + 2) + &
              cod(2) * expo * t2pb_mic_pre(vloop_mic_pre+1)
            end do
          end if
        end do
        t2pb_y = t2pb_y + &
        binomialcoe(m_y_i,vloop_i) * binomialcoe(m_y_j,vloop_j) * &
        (-cod_iy)**(vloop_i) * (-cod_jy)**(vloop_j) * t2pb_mic
        vloop_j = vloop_j + 1
      end do
      vloop_i = vloop_i + 1 
    end do
    !--------------------------
    ! integral of z, generate a coefficient exp(-b*((cod(3))^2*t^2)/(t^2+b))
    t2pb_z = 0.0_dp
    vloop_i = 0
    do while(vloop_i <= m_z_i)
      vloop_j = 0
      do while(vloop_j <= m_z_j)
        t2pb_mic = 0.0_dp
        t2pb_mic_pre = 0.0_dp
        t2pb_mic_pre_pre = 0.0_dp
        do vloop_mic = 0, m_z_i + m_z_j - vloop_i - vloop_j
          if (vloop_mic == 0) then
            t2pb_mic(1) = dsqrt(pi)
          else if(vloop_mic == 1) then
            t2pb_mic_pre(1) = dsqrt(pi)
            t2pb_mic = 0.0_dp
            t2pb_mic(2) = cod(3) * expo * t2pb_mic_pre(1)
          else
            t2pb_mic_pre_pre = t2pb_mic_pre
            t2pb_mic_pre = t2pb_mic
            t2pb_mic = 0.0_dp
            do vloop_mic_pre = 0, vloop_mic
              t2pb_mic(vloop_mic_pre + 2) = t2pb_mic(vloop_mic_pre + 2) + &
              0.5_dp * real(vloop_mic-1) * t2pb_mic_pre_pre(vloop_mic_pre+1)
              t2pb_mic(vloop_mic_pre + 2) = t2pb_mic(vloop_mic_pre + 2) + &
              cod(3) * expo * t2pb_mic_pre(vloop_mic_pre+1)
            end do
          end if
        end do
        t2pb_z = t2pb_z + &
        binomialcoe(m_z_i,vloop_i) * binomialcoe(m_z_j,vloop_j) * &
        (-cod_iz)**(vloop_i) * (-cod_jz)**(vloop_j) * t2pb_mic
        vloop_j = vloop_j + 1
      end do
      vloop_i = vloop_i + 1 
    end do
    !--------------------------
    ! integral of t
    ! 2*(-1)^k*b^((1-m)/2)*coe*(u-1)^k*(u+1)^k*exp(-b*R2*u^2)
    t2pb = 0.0_dp
    do vloop_i = 0, m_x_i + m_x_j + 1
      do vloop_j = 0, m_y_i + m_y_j + 1
        do vloop_k = 0, m_z_i + m_z_j + 1
          t2pb(vloop_i+vloop_j+vloop_k+2) = t2pb(vloop_i+vloop_j+vloop_k+2)+ &
          t2pb_x(vloop_i+1) * t2pb_y(vloop_j+1) * t2pb_z(vloop_k+1)
        end do
      end do
    end do
    V_Integral_1e = 0.0_dp
    if (abs(expo*R2) < 1E-13) then
      vloop_i = 0 ! k = vloop_i, t2pb(vloop_i + 3) = coe
      do while(vloop_i <= m_x_i + m_x_j + m_y_i + m_y_j + m_z_i + m_z_j + 4)
        integral = 0.0_dp
        do vloop_j = 0, vloop_i
          do vloop_k = 0, vloop_i
            integral_mic = 1.0_dp / (real(2*vloop_i-vloop_j-vloop_k) + 1.0_dp)
            integral = integral + &
            binomialcoe(vloop_i,vloop_j) * binomialcoe(vloop_i,vloop_k) * &
            (-1.0_dp)**(vloop_j) * integral_mic
          end do
        end do
        V_Integral_1e = V_Integral_1e + 2.0_dp * &
        (-1.0_dp)**(vloop_i) * expo**(-vloop_i-1) * t2pb(vloop_i+2) * integral
        vloop_i = vloop_i + 1
      end do
    else if (abs(expo*R2) <= xts) then
      vloop_i = 0 ! k = vloop_i, t2pb(vloop_i + 3) = coe
      do while(vloop_i <= m_x_i + m_x_j + m_y_i + m_y_j + m_z_i + m_z_j + 4)
        integral = 0.0_dp
        vloop_j = 0
        do while(vloop_j <= vloop_i)
          vloop_k = 0
          do while(vloop_k <= vloop_i)
            vloop_mic = 0
            do while(vloop_mic <= 2 * vloop_i - vloop_j - vloop_k)
              if (vloop_mic == 0) then 
                do vloop_o=1,tayeps
                  NX_mic(vloop_o)=(-1.0_dp)**(vloop_o)*&
                  (1.0_dp/(real(2*vloop_o+1)*factorial(vloop_o)))
                end do
              else if(vloop_mic == 1) then
                NX_mic_pre = NX_mic
                do vloop_o=1,tayeps
                  NX_mic(vloop_o)=(-1.0_dp)**(vloop_o)*&
                  (1.0_dp/(factorial(vloop_o+1)))
                end do
              else
                NX_mic_pre_pre = NX_mic_pre
                NX_mic_pre = NX_mic
                do vloop_o=1,tayeps-1
                  NX_mic(vloop_o)=((-1.0_dp)**(vloop_o)*&
                  (0.5_dp/(factorial(vloop_o+1)))+(NX_mic_pre_pre(vloop_o+1))&
                  /2.0_dp) * real(vloop_mic+1)
                end do
                NX_mic(tayeps)=((-1.0_dp)**tayeps*&
                (0.5_dp/(factorial(tayeps+1)))) * real(vloop_mic+1)
              end if
              vloop_mic = vloop_mic + 1
            end do
            if (2 * vloop_i - vloop_j - vloop_k == 0) then
              integral_mic = 0.0_dp
              do vloop_o=1,tayeps
                integral_mic = integral_mic + &
                (-1.0_dp)**(vloop_o+1)*(expo*R2)**(vloop_o-1)*&
                (1.0_dp/(real(2*(vloop_o)-1)*factorial(vloop_o-1)))
              end do
            else if (2 * vloop_i - vloop_j - vloop_k == 1) then
              integral_mic = 0.0_dp
              do vloop_o=1,tayeps
                integral_mic = integral_mic + &
                (-1.0_dp)**(vloop_o+1)*(expo*R2)**(vloop_o-1)*&
                (0.5_dp/factorial(vloop_o))
              end do
            else
              integral_mic = 0.0_dp
              do vloop_o=1,tayeps
                integral_mic = integral_mic + &
                (-1.0_dp)**(vloop_o+1)*(expo*R2)**(vloop_o-1)*&
                (0.5_dp/factorial(vloop_o))&
                  + (expo*R2)**(vloop_o-1)*NX_mic_pre_pre(vloop_o)/2.0_dp
              end do
            end if
            integral = integral + &
            binomialcoe(vloop_i,vloop_j) * binomialcoe(vloop_i,vloop_k) * &
            (-1.0_dp)**(vloop_j) * integral_mic
            vloop_k = vloop_k + 1
          end do
          vloop_j = vloop_j + 1
        end do
        V_Integral_1e = V_Integral_1e + 2.0_dp * &
        (-1.0_dp)**(vloop_i) * expo**(-vloop_i-1) * t2pb(vloop_i+2) * integral
        vloop_i = vloop_i + 1
      end do
    else
      vloop_i = 0 ! k = vloop_i, t2pb(vloop_i + 3) = coe
      do while(vloop_i <= m_x_i + m_x_j + m_y_i + m_y_j + m_z_i + m_z_j + 4)
        integral = 0.0_dp
        do vloop_j = 0, vloop_i
          do vloop_k = 0, vloop_i
            do vloop_mic = 0, 2*vloop_i - vloop_j - vloop_k
              if (vloop_mic == 0) then
                integral_mic = 0.5_dp*dsqrt(pi/(expo*R2)) * erf(dsqrt(expo*R2))
              else if(vloop_mic == 1) then
                integral_mic_pre = integral_mic
                integral_mic = (1.0_dp - exp(-expo*R2)) / (2.0_dp*expo*R2)
              else
                integral_mic_pre_pre = integral_mic_pre
                integral_mic_pre = integral_mic
                integral_mic = (-(exp(-expo*R2)) + real(vloop_mic - 1) * &
                integral_mic_pre_pre) / (2.0_dp*expo*R2)
              end if
            end do
            integral = integral + &
            binomialcoe(vloop_i,vloop_j) * binomialcoe(vloop_i,vloop_k) * &
            (-1.0_dp)**(vloop_j) * integral_mic
          end do
        end do
        V_Integral_1e = V_Integral_1e + 2.0_dp * &
        (-1.0_dp)**(vloop_i) * expo**(-vloop_i-1) * t2pb(vloop_i+2) * integral
        vloop_i = vloop_i + 1
      end do
    end if
    V_Integral_1e = V_Integral_1e / dsqrt(pi)
    if (finitenuc) then
      ! Finite nuclear model correction, spherical electric charge
      V_Integral_1e = V_Integral_1e - &
      2.0_dp/3.0_dp * pi * rn * rn * exp(-expo*R2) * (-cod_ix)**(m_x_i) * &
      (-cod_iy)**(m_y_i) * (-cod_iz)**(m_z_i) * (-cod_jx)**(m_x_j) * &
      (-cod_jy)**(m_y_j) * (-cod_jz)**(m_z_j)
    end if
    V_Integral_1e = V_Integral_1e * coe
    return
  end function V_Integral_1e
  
!-----------------------------------------------------------------------
!> integration of two electron repulsion potential in Cartesian coordinate
!!
!! EXPRESS: |AOi>:(x1 - xi), |AOj>:(x1 - xj), |AOk>:(x2 - xk), |AOl>:(x2 - xl)
!!
!! Li > Lj, Lk > Ll
  pure elemental real(dp) function V_Integral_2e(&
  fac_i,fac_j,fac_k,fac_l,ai,aj,ak,al,cod_ix,cod_iy,cod_iz,&
  cod_jx,cod_jy,cod_jz,cod_kx,cod_ky,cod_kz,cod_lx,cod_ly,cod_lz)
    implicit none
    real(dp),intent(in) :: ai, aj, ak, al
    real(dp),intent(in) :: cod_ix,cod_iy,cod_iz,cod_jx,cod_jy,cod_jz
    real(dp),intent(in) :: cod_kx,cod_ky,cod_kz,cod_lx,cod_ly,cod_lz
    real(dp) :: xi, xj, xk, xl, yi, yj, yk, yl, zi, zj, zk, zl
    ! composite parameters, ref 10.1002/jcc.540040206
    real(dp) :: xA, xB, yA, yB, zA, zB, A, B, rou, Dx, Dy, Dz, X, Gx, Gy, Gz
    integer :: nroots                   ! Rys quadrature: number of roots
    ! Rys quadrature: roots and weights ref 10.1016/0021-9991(76)90008-5
    real(dp) :: u(6), w(6), t2
    real(dp) :: Gnmx(9,9,18),Gnmy(9,9,18),Gnmz(9,9,18)
    ! ni transfer to nj, nk tranfer to nl
    real(dp) :: Ixtrans(5,18),Ix(18),Iytrans(5,18),Iy(18)
    real(dp) :: Iztrans(5,18),Iz(18),PL(22)
    real(dp) :: integral, integral_mic, integral_mic_pre, integral_mic_pre_pre
    real(dp) :: NX_mic(30), NX_mic_pre(30), NX_mic_pre_pre(30)
    character(len = *),intent(in) :: fac_i,fac_j,fac_k,fac_l
    integer :: tayeps  ! number of Taylor expansion series of integration at X=0
    ! direct integration (X > xts); Taylor expansion integration (X <= xts)
    real(dp) :: xts
    integer :: nxi,nxj,nxk,nxl,nyi,nyj,nyk,nyl,nzi,nzj,nzk,nzl
    integer :: rloop_i,rloop_j,rloop_k  ! loop variables only for V_Integral_2e
    xi = cod_ix
    xj = cod_jx
    xk = cod_kx
    xl = cod_lx
    yi = cod_iy
    yj = cod_jy
    yk = cod_ky
    yl = cod_ly
    zi = cod_iz
    zj = cod_jz
    zk = cod_kz
    zl = cod_lz
    nxi = 0
    nxj = 0
    nxk = 0
    nxl = 0
    nyi = 0
    nyj = 0
    nyk = 0
    nyl = 0
    nzi = 0
    nzj = 0
    nzk = 0
    nzl = 0
    do rloop_i = 1, len(fac_i)
      if(fac_i(rloop_i:rloop_i) == 'x') nxi = nxi + 1
      if(fac_i(rloop_i:rloop_i) == 'y') nyi = nyi + 1
      if(fac_i(rloop_i:rloop_i) == 'z') nzi = nzi + 1
    end do
    do rloop_i = 1, len(fac_j)
      if(fac_j(rloop_i:rloop_i) == 'x') nxj = nxj + 1
      if(fac_j(rloop_i:rloop_i) == 'y') nyj = nyj + 1
      if(fac_j(rloop_i:rloop_i) == 'z') nzj = nzj + 1
    end do
    do rloop_i = 1, len(fac_k)
      if(fac_k(rloop_i:rloop_i) == 'x') nxk = nxk + 1
      if(fac_k(rloop_i:rloop_i) == 'y') nyk = nyk + 1
      if(fac_k(rloop_i:rloop_i) == 'z') nzk = nzk + 1
    end do
    do rloop_i = 1, len(fac_l)
      if(fac_l(rloop_i:rloop_i) == 'x') nxl = nxl + 1
      if(fac_l(rloop_i:rloop_i) == 'y') nyl = nyl + 1
      if(fac_l(rloop_i:rloop_i) == 'z') nzl = nzl + 1
    end do
    ! Gaussian product
    xA = (ai * xi + aj * xj) / (ai + aj)
    xB = (ak * xk + al * xl) / (ak + al)
    yA = (ai * yi + aj * yj) / (ai + aj)
    yB = (ak * yk + al * yl) / (ak + al)
    zA = (ai * zi + aj * zj) / (ai + aj)
    zB = (ak * zk + al * zl) / (ak + al)
    A = ai + aj
    B = ak + al
    rou = A * B / (A + B)
    Dx = rou * (xA - xB)**2
    Dy = rou * (yA - yB)**2
    Dz = rou * (zA - zB)**2
    X = Dx + Dy + Dz
    ! for non-normalized inputs, don't consider Gx, Gy, and Gz.
    Gx = (ai * aj / (ai + aj)) * (xi - xj)**2 + &
    (ak * al / (ak + al)) * (xk - xl)**2
    Gy = (ai * aj / (ai + aj)) * (yi - yj)**2 + &
    (ak * al / (ak + al)) * (yk - yl)**2
    Gz = (ai * aj / (ai + aj)) * (zi - zj)**2 + &
    (ak * al / (ak + al)) * (zk - zl)**2
    ! Rys quadrature scheme for low angular momentum Gaussian functions
    if (int((nxi+nxj+nxk+nxl+nyi+nyj+nyk+nyl+nzi+nzj+nzk+nzl)/2.0)+1 <= 5) then
      nroots = int((nxi+nxj+nxk+nxl+nyi+nyj+nyk+nyl+nzi+nzj+nzk+nzl)/2.0)+1
      call GRysroots(nroots, X, u, w)  ! rys_roots (libcint), GRysroots (GAMESS)
      integral = 0.0_dp
      do rloop_i = 1, nroots
        !t2 = u(rloop_i) / (rou + u(rloop_i))     ! for rys_roots
        t2 = u(rloop_i)                           ! for GRysroots
        !---------------------------Gnmx---------------------------
        Gnmx = 0.0_dp
        Gnmx(1,1,1) = pi / dsqrt(A*B) * exp(-Gx)
        do rloop_j = 2, nxi + nxj + 1
          if (rloop_j == 2) then
            Gnmx(2,1,1) = Gnmx(2,1,1) + &
            Gnmx(1,1,1)*(xA-xi) + Gnmx(1,1,1)*(B*(xB-xA)/(A+B)) * t2
          else
            Gnmx(rloop_j,1,1) = Gnmx(rloop_j,1,1) + &
            Gnmx(rloop_j-2,1,1) * real(rloop_j-2)/(2.0_dp*A) + &
            Gnmx(rloop_j-1,1,1) * (xA - xi)
            Gnmx(rloop_j,1,1) = Gnmx(rloop_j,1,1) - &
            Gnmx(rloop_j-2,1,1) * real(rloop_j-2)*B/(2.0_dp*A*(A+B)) * t2 + &
            Gnmx(rloop_j-1,1,1) * (B*(xB-xA)/(A+B)) * t2
          end if
        end do
        do rloop_j = 2, nxk + nxl + 1
          if (rloop_j == 2) then
            Gnmx(1,2,1) = Gnmx(1,2,1) + Gnmx(1,1,1)*(xB-xk) + &
            Gnmx(1,1,1)*(A*(xA-xB)/(A+B)) * t2
          else
            Gnmx(1,rloop_j,1) = Gnmx(1,rloop_j,1) + &
            Gnmx(1,rloop_j-2,1) * real(rloop_j-2)/(2.0_dp*B) + &
            Gnmx(1,rloop_j-1,1) * (xB - xk)
            Gnmx(1,rloop_j,1) = Gnmx(1,rloop_j,1) - &
            Gnmx(1,rloop_j-2,1) * real(rloop_j-2)*A/(2.0_dp*B*(A+B)) * t2 + &
            Gnmx(1,rloop_j-1,1) * (A*(xA-xB)/(A+B)) * t2
          end if
        end do
        do rloop_j = 2, nxk + nxl + 1
          do rloop_k = 2, nxi + nxj + 1
            if (rloop_k == 2) then
              Gnmx(2,rloop_j,1) = Gnmx(2,rloop_j,1) + &
              Gnmx(1,rloop_j,1) * (xA - xi)
              Gnmx(2,rloop_j,1) = Gnmx(2,rloop_j,1) + &
              Gnmx(1,rloop_j,1) * (B*(xB-xA)/(A+B)) * t2 + &
              Gnmx(1,rloop_j-1,1) * real(rloop_j-1)/(2.0_dp*(A+B)) * t2
            else
              Gnmx(rloop_k,rloop_j,1) = Gnmx(rloop_k,rloop_j,1) + &
              Gnmx(rloop_k-2,rloop_j,1) * real(rloop_k-2)/(2.0_dp*A) + &
              Gnmx(rloop_k-1,rloop_j,1) * (xA - xi)
              Gnmx(rloop_k,rloop_j,1) = Gnmx(rloop_k,rloop_j,1) - &
              Gnmx(rloop_k-2,rloop_j,1) * real(rloop_k-2)*B/(2.0_dp*A*(A+B)) * &
              t2 + Gnmx(rloop_k-1,rloop_j,1) * (B*(xB-xA)/(A+B)) * t2 + &
              Gnmx(rloop_k-1,rloop_j-1,1) * real(rloop_j-1)/(2.0_dp*(A+B)) * t2
            end if
          end do
        end do
        !---------------------------Gnmy---------------------------
        Gnmy = 0.0_dp
        Gnmy(1,1,1) = pi / dsqrt(A*B) * exp(-Gy)
        do rloop_j = 2, nyi + nyj + 1
          if (rloop_j == 2) then
            Gnmy(2,1,1) = Gnmy(2,1,1) + &
            Gnmy(1,1,1)*(yA-yi) + Gnmy(1,1,1)*(B*(yB-yA)/(A+B)) * t2
          else
            Gnmy(rloop_j,1,1) = Gnmy(rloop_j,1,1) + &
            Gnmy(rloop_j-2,1,1) * real(rloop_j-2)/(2.0_dp*A) + &
            Gnmy(rloop_j-1,1,1) * (yA - yi)
            Gnmy(rloop_j,1,1) = Gnmy(rloop_j,1,1) - &
            Gnmy(rloop_j-2,1,1) * real(rloop_j-2)*B/(2.0_dp*A*(A+B)) * t2 + &
            Gnmy(rloop_j-1,1,1) * (B*(yB-yA)/(A+B)) * t2
          end if
        end do
        do rloop_j = 2, nyk + nyl + 1
          if (rloop_j == 2) then
            Gnmy(1,2,1) = Gnmy(1,2,1) + &
            Gnmy(1,1,1)*(yB-yk) + Gnmy(1,1,1)*(A*(yA-yB)/(A+B)) * t2
          else
            Gnmy(1,rloop_j,1) = Gnmy(1,rloop_j,1) + &
            Gnmy(1,rloop_j-2,1) * real(rloop_j-2)/(2.0_dp*B) + &
            Gnmy(1,rloop_j-1,1) * (yB - yk)
            Gnmy(1,rloop_j,1) = Gnmy(1,rloop_j,1) - &
            Gnmy(1,rloop_j-2,1) * real(rloop_j-2)*A/(2.0_dp*B*(A+B)) * t2 + &
            Gnmy(1,rloop_j-1,1) * (A*(yA-yB)/(A+B)) * t2
          end if
        end do
        do rloop_j = 2, nyk + nyl + 1
          do rloop_k = 2, nyi + nyj + 1
            if (rloop_k == 2) then
              Gnmy(2,rloop_j,1) = Gnmy(2,rloop_j,1) + &
              Gnmy(1,rloop_j,1) * (yA - yi)
              Gnmy(2,rloop_j,1) = Gnmy(2,rloop_j,1) + &
              Gnmy(1,rloop_j,1) * (B*(yB-yA)/(A+B)) * t2 + &
              Gnmy(1,rloop_j-1,1) * real(rloop_j-1)/(2.0_dp*(A+B)) * t2
            else
              Gnmy(rloop_k,rloop_j,1) = Gnmy(rloop_k,rloop_j,1) + &
              Gnmy(rloop_k-2,rloop_j,1) * real(rloop_k-2)/(2.0_dp*A) + &
              Gnmy(rloop_k-1,rloop_j,1) * (yA - yi)
              Gnmy(rloop_k,rloop_j,1) = Gnmy(rloop_k,rloop_j,1) - &
              Gnmy(rloop_k-2,rloop_j,1) * real(rloop_k-2)*B/(2.0_dp*A*(A+B)) * &
              t2 + Gnmy(rloop_k-1,rloop_j,1) * (B*(yB-yA)/(A+B)) * t2 + &
              Gnmy(rloop_k-1,rloop_j-1,1) * real(rloop_j-1)/(2.0_dp*(A+B)) * t2
            end if
          end do
        end do
        !---------------------------Gnmz---------------------------
        Gnmz = 0.0_dp
        Gnmz(1,1,1) = pi / dsqrt(A*B) * exp(-Gz)
        do rloop_j = 2, nzi + nzj + 1
          if (rloop_j == 2) then
            Gnmz(2,1,1) = Gnmz(2,1,1) + &
            Gnmz(1,1,1)*(zA-zi) + Gnmz(1,1,1)*(B*(zB-zA)/(A+B)) * t2
          else
            Gnmz(rloop_j,1,1) = Gnmz(rloop_j,1,1) + &
            Gnmz(rloop_j-2,1,1) * real(rloop_j-2)/(2.0_dp*A) + &
            Gnmz(rloop_j-1,1,1) * (zA - zi)
            Gnmz(rloop_j,1,1) = Gnmz(rloop_j,1,1) - &
            Gnmz(rloop_j-2,1,1) * real(rloop_j-2)*B/(2.0_dp*A*(A+B)) * t2 + &
            Gnmz(rloop_j-1,1,1) * (B*(zB-zA)/(A+B)) * t2
          end if
        end do
        do rloop_j = 2, nzk + nzl + 1
          if (rloop_j == 2) then
            Gnmz(1,2,1) = Gnmz(1,2,1) + &
            Gnmz(1,1,1)*(zB-zk) + Gnmz(1,1,1)*(A*(zA-zB)/(A+B)) * t2
          else
            Gnmz(1,rloop_j,1) = Gnmz(1,rloop_j,1) + &
            Gnmz(1,rloop_j-2,1) * real(rloop_j-2)/(2.0_dp*B) + &
            Gnmz(1,rloop_j-1,1) * (zB - zk)
            Gnmz(1,rloop_j,1) = Gnmz(1,rloop_j,1) - &
            Gnmz(1,rloop_j-2,1) * real(rloop_j-2)*A/(2.0_dp*B*(A+B)) * t2 + &
            Gnmz(1,rloop_j-1,1) * (A*(zA-zB)/(A+B)) * t2
          end if
        end do
        do rloop_j = 2, nzk + nzl + 1
          do rloop_k = 2, nzi + nzj + 1
            if (rloop_k == 2) then
              Gnmz(2,rloop_j,1) = Gnmz(2,rloop_j,1) + &
              Gnmz(1,rloop_j,1) * (zA - zi)
              Gnmz(2,rloop_j,1) = Gnmz(2,rloop_j,1) + &
              Gnmz(1,rloop_j,1) * (B*(zB-zA)/(A+B)) * t2 + &
              Gnmz(1,rloop_j-1,1) * real(rloop_j-1)/(2.0_dp*(A+B)) * t2
            else
              Gnmz(rloop_k,rloop_j,1) = Gnmz(rloop_k,rloop_j,1) + &
              Gnmz(rloop_k-2,rloop_j,1) * real(rloop_k-2)/(2.0_dp*A) + &
              Gnmz(rloop_k-1,rloop_j,1) * (zA - zi)
              Gnmz(rloop_k,rloop_j,1) = Gnmz(rloop_k,rloop_j,1) - &
              Gnmz(rloop_k-2,rloop_j,1) * real(rloop_k-2)*B/(2.0_dp*A*(A+B)) * &
              t2 + Gnmz(rloop_k-1,rloop_j,1) * (B*(zB-zA)/(A+B)) * t2 + &
              Gnmz(rloop_k-1,rloop_j-1,1) * real(rloop_j-1)/(2.0_dp*(A+B)) * t2
            end if
          end do
        end do
        !---------------------------Ix---------------------------
        Ixtrans = 0.0_dp
        Ix = 0.0_dp
        do rloop_k = 1, nxl + 1
          do rloop_j = 0, nxj
            Ixtrans(rloop_k,1) = Ixtrans(rloop_k,1) + &
            binomialcoe(nxj, rloop_j) * (xi-xj)**(rloop_j) * &
            Gnmx(nxi+nxj+1-rloop_j, nxk+nxl+1-(rloop_k-1), 1)
          end do
        end do
        do rloop_j = 0, nxl
          Ix(1) = Ix(1) + &
          binomialcoe(nxl,rloop_j) * (xk-xl)**(rloop_j) * Ixtrans(rloop_j+1, 1)
        end do
        !---------------------------Iy---------------------------
        Iytrans = 0.0_dp
        Iy = 0.0_dp
        do rloop_k = 1, nyl + 1
          do rloop_j = 0, nyj
            Iytrans(rloop_k,1) = Iytrans(rloop_k,1) + &
            binomialcoe(nyj, rloop_j) * (yi-yj)**(rloop_j) * &
            Gnmy(nyi+nyj+1-rloop_j, nyk+nyl+1-(rloop_k-1), 1)
          end do
        end do
        do rloop_j = 0, nyl
          Iy(1) = Iy(1) + &
          binomialcoe(nyl,rloop_j) * (yk-yl)**(rloop_j) * Iytrans(rloop_j+1, 1)
        end do
        !---------------------------Iz---------------------------
        Iztrans = 0.0_dp
        Iz = 0.0_dp
        do rloop_k = 1, nzl + 1
          do rloop_j = 0, nzj
            Iztrans(rloop_k,1) = Iztrans(rloop_k,1) + &
            binomialcoe(nzj, rloop_j) * (zi-zj)**(rloop_j) * &
            Gnmz(nzi+nzj+1-rloop_j, nzk+nzl+1-(rloop_k-1), 1)
          end do
        end do
        do rloop_j = 0, nzl
          Iz(1) = Iz(1) + &
          binomialcoe(nzl,rloop_j) * (zk-zl)**(rloop_j) * Iztrans(rloop_j+1, 1)
        end do
        !---------------------------PL---------------------------
        PL(1) = Ix(1) * Iy(1) * Iz(1) * 2.0_dp * dsqrt(rou/pi)
        integral = integral + w(rloop_i)*PL(1)
      end do
      V_Integral_2e = integral
      return
    ! direct integral for high angular momentum Gaussian functions
    else
      ! Ix(ni+nj,0,nk+nl,0,u)
      if (2*(nxi+nxj+1+nxk+nxl+1+nyi+nyj+1+nyk+nyl+1+nzi+nzj+1+nzk+&
      nzl+1)-2 <= 15) then
        tayeps = 8
        xts = 0.1
      else if (2*(nxi+nxj+1+nxk+nxl+1+nyi+nyj+1+nyk+nyl+1+nzi+nzj+1+&
      nzk+nzl+1)-2 <= 25) then
        tayeps = 15
        xts = 1.0
      else if (2*(nxi+nxj+1+nxk+nxl+1+nyi+nyj+1+nyk+nyl+1+nzi+nzj+1+&
      nzk+nzl+1)-2 <= 35) then
        tayeps = 20
        xts = 1.0
      else if (2*(nxi+nxj+1+nxk+nxl+1+nyi+nyj+1+nyk+nyl+1+nzi+nzj+1+&
      nzk+nzl+1)-2 <= 40) then
        tayeps = 25
        xts = 3.0
      else if (2*(nxi+nxj+1+nxk+nxl+1+nyi+nyj+1+nyk+nyl+1+nzi+nzj+1+&
      nzk+nzl+1)-2 <= 50) then
        tayeps = 30
        xts = 4.0
      end if
      NX_mic = 0.0_dp
      NX_mic_pre = 0.0_dp
      NX_mic_pre_pre = 0.0_dp

      Gnmx = 0.0_dp
      Gnmy = 0.0_dp
      Gnmz = 0.0_dp
      ! reduce the factor (1-t^2)^(1/2)*exp(-Dx*t^2)
      Gnmx(1,1,1) = pi / dsqrt(A*B) * exp(-Gx)
      do rloop_i = 2, nxi + nxj + 1
        if (rloop_i == 2) then
          Gnmx(2,1,1) = Gnmx(2,1,1) + Gnmx(1,1,1) * (xA - xi)
          Gnmx(2,1,2) = Gnmx(2,1,2) + Gnmx(1,1,1) * (B*(xB-xA)/(A+B))
        else
          do rloop_j = 1, rloop_i - 1
            Gnmx(rloop_i,1,rloop_j) = Gnmx(rloop_i,1,rloop_j) + &
            Gnmx(rloop_i-2,1,rloop_j) * real(rloop_i-2)/(2.0_dp*A) + &
            Gnmx(rloop_i-1,1,rloop_j) * (xA - xi)
            Gnmx(rloop_i,1,rloop_j + 1) = Gnmx(rloop_i,1,rloop_j + 1) - &
            Gnmx(rloop_i-2,1,rloop_j) * real(rloop_i-2)*B/(2.0_dp*A*(A+B)) + &
            Gnmx(rloop_i-1,1,rloop_j) * (B*(xB-xA)/(A+B))
          end do
        end if
      end do
      do rloop_i = 2, nxk + nxl + 1
        if (rloop_i == 2) then
          Gnmx(1,2,1) = Gnmx(1,2,1) + Gnmx(1,1,1) * (xB - xk)
          Gnmx(1,2,2) = Gnmx(1,2,2) + Gnmx(1,1,1) * (A*(xA-xB)/(A+B))
        else
          do rloop_j = 1, rloop_i - 1
            Gnmx(1,rloop_i,rloop_j) = Gnmx(1,rloop_i,rloop_j) + &
            Gnmx(1,rloop_i-2,rloop_j) * real(rloop_i-2)/(2.0_dp*B) + &
            Gnmx(1,rloop_i-1,rloop_j) * (xB - xk)
            Gnmx(1,rloop_i,rloop_j + 1) = Gnmx(1,rloop_i,rloop_j + 1) - &
            Gnmx(1,rloop_i-2,rloop_j) * real(rloop_i-2)*A/(2.0_dp*B*(A+B)) + &
            Gnmx(1,rloop_i-1,rloop_j) * (A*(xA-xB)/(A+B))
          end do
        end if
      end do
      ! use G(n+1,m) recursion only, xA and xB are asymmetric
      ! ---------------------------------------------------------------------
      do rloop_k = 2, nxk + nxl + 1
        do rloop_i = 2, nxi + nxj + 1
          if (rloop_i == 2) then
            do rloop_j = 1, rloop_k
              Gnmx(2,rloop_k,rloop_j) = Gnmx(2,rloop_k,rloop_j) + &
              Gnmx(1,rloop_k,rloop_j) * (xA - xi)
              Gnmx(2,rloop_k,rloop_j + 1) = Gnmx(2,rloop_k,rloop_j + 1) + &
              Gnmx(1,rloop_k,rloop_j) * (B*(xB-xA)/(A+B)) + &
              Gnmx(1,rloop_k-1,rloop_j) * real(rloop_k-1)/(2.0_dp*(A+B))
            end do
          else
            do rloop_j = 1, rloop_i + rloop_k - 2
              Gnmx(rloop_i,rloop_k,rloop_j) = Gnmx(rloop_i,rloop_k,rloop_j) + &
              Gnmx(rloop_i-2,rloop_k,rloop_j) * real(rloop_i-2)/(2.0_dp*A) + &
              Gnmx(rloop_i-1,rloop_k,rloop_j) * (xA - xi)
              Gnmx(rloop_i,rloop_k,rloop_j+1) = Gnmx(rloop_i,rloop_k,rloop_j+1)&
              - Gnmx(rloop_i-2,rloop_k,rloop_j) * real(rloop_i-2)*B/(2.0_dp*A*&
              (A+B)) + Gnmx(rloop_i-1,rloop_k,rloop_j) * (B*(xB-xA)/(A+B)) + &
              Gnmx(rloop_i-1,rloop_k-1,rloop_j) * real(rloop_k-1)/(2.0_dp*(A+B))
            end do
          end if
        end do
      end do
      ! ---------------------------------------------------------------------
      ! Iy(ni+nj,0,nk+nl,0,u)
      ! reduce the factor (1-t^2)^(1/2)*eyp(-Dy*t^2)
      Gnmy(1,1,1) = pi / dsqrt(A*B) * exp(-Gy)
      do rloop_i = 2, nyi + nyj + 1
        if (rloop_i == 2) then
          Gnmy(2,1,1) = Gnmy(2,1,1) + Gnmy(1,1,1) * (yA - yi)
          Gnmy(2,1,2) = Gnmy(2,1,2) + Gnmy(1,1,1) * (B*(yB-yA)/(A+B))
        else
          do rloop_j = 1, rloop_i - 1
            Gnmy(rloop_i,1,rloop_j) = Gnmy(rloop_i,1,rloop_j) + &
            Gnmy(rloop_i-2,1,rloop_j) * real(rloop_i-2)/(2.0_dp*A) + &
            Gnmy(rloop_i-1,1,rloop_j) * (yA - yi)
            Gnmy(rloop_i,1,rloop_j + 1) = Gnmy(rloop_i,1,rloop_j + 1) - &
            Gnmy(rloop_i-2,1,rloop_j) * real(rloop_i-2)*B/(2.0_dp*A*(A+B)) + &
            Gnmy(rloop_i-1,1,rloop_j) * (B*(yB-yA)/(A+B))
          end do
        end if
      end do
      do rloop_i = 2, nyk + nyl + 1
        if (rloop_i == 2) then
          Gnmy(1,2,1) = Gnmy(1,2,1) + Gnmy(1,1,1) * (yB - yk)
          Gnmy(1,2,2) = Gnmy(1,2,2) + Gnmy(1,1,1) * (A*(yA-yB)/(A+B))
        else
          do rloop_j = 1, rloop_i - 1
            Gnmy(1,rloop_i,rloop_j) = Gnmy(1,rloop_i,rloop_j) + &
            Gnmy(1,rloop_i-2,rloop_j) * real(rloop_i-2)/(2.0_dp*B) + &
            Gnmy(1,rloop_i-1,rloop_j) * (yB - yk)
            Gnmy(1,rloop_i,rloop_j + 1) = Gnmy(1,rloop_i,rloop_j + 1) - &
            Gnmy(1,rloop_i-2,rloop_j) * real(rloop_i-2)*A/(2.0_dp*B*(A+B)) + &
            Gnmy(1,rloop_i-1,rloop_j) * (A*(yA-yB)/(A+B))
          end do
        end if
      end do
      do rloop_k = 2, nyk + nyl + 1
        do rloop_i = 2, nyi + nyj + 1
          if (rloop_i == 2) then
            do rloop_j = 1, rloop_k
              Gnmy(2,rloop_k,rloop_j) = Gnmy(2,rloop_k,rloop_j) + &
              Gnmy(1,rloop_k,rloop_j) * (yA - yi)
              Gnmy(2,rloop_k,rloop_j + 1) = Gnmy(2,rloop_k,rloop_j + 1) + &
              Gnmy(1,rloop_k,rloop_j) * (B*(yB-yA)/(A+B)) + &
              Gnmy(1,rloop_k-1,rloop_j) * real(rloop_k-1)/(2.0_dp*(A+B))
            end do
          else
            do rloop_j = 1, rloop_i + rloop_k - 2
              Gnmy(rloop_i,rloop_k,rloop_j) = Gnmy(rloop_i,rloop_k,rloop_j) + &
              Gnmy(rloop_i-2,rloop_k,rloop_j) * real(rloop_i-2)/(2.0_dp*A) + &
              Gnmy(rloop_i-1,rloop_k,rloop_j) * (yA - yi)
              Gnmy(rloop_i,rloop_k,rloop_j+1) = Gnmy(rloop_i,rloop_k,rloop_j+1)&
              - Gnmy(rloop_i-2,rloop_k,rloop_j) * real(rloop_i-2)*B/(2.0_dp*A*&
              (A+B)) + Gnmy(rloop_i-1,rloop_k,rloop_j) * (B*(yB-yA)/(A+B)) + &
              Gnmy(rloop_i-1,rloop_k-1,rloop_j) * real(rloop_k-1)/(2.0_dp*(A+B))
            end do
          end if
        end do
      end do
      ! Iz(ni+nj,0,nk+nl,0,u)
      ! reduce the factor (1-t^2)^(1/2)*ezp(-Dz*t^2)
      Gnmz(1,1,1) = pi / dsqrt(A*B) * exp(-Gz)
      do rloop_i = 2, nzi + nzj + 1
        if (rloop_i == 2) then
          Gnmz(2,1,1) = Gnmz(2,1,1) + Gnmz(1,1,1) * (zA - zi)
          Gnmz(2,1,2) = Gnmz(2,1,2) + Gnmz(1,1,1) * (B*(zB-zA)/(A+B))
        else
          do rloop_j = 1, rloop_i - 1
            Gnmz(rloop_i,1,rloop_j) = Gnmz(rloop_i,1,rloop_j) + &
            Gnmz(rloop_i-2,1,rloop_j) * real(rloop_i-2)/(2.0_dp*A) + &
            Gnmz(rloop_i-1,1,rloop_j) * (zA - zi)
            Gnmz(rloop_i,1,rloop_j + 1) = Gnmz(rloop_i,1,rloop_j + 1) - &
            Gnmz(rloop_i-2,1,rloop_j) * real(rloop_i-2)*B/(2.0_dp*A*(A+B)) + &
            Gnmz(rloop_i-1,1,rloop_j) * (B*(zB-zA)/(A+B))
          end do
        end if
      end do
      do rloop_i = 2, nzk + nzl + 1
        if (rloop_i == 2) then
          Gnmz(1,2,1) = Gnmz(1,2,1) + Gnmz(1,1,1) * (zB - zk)
          Gnmz(1,2,2) = Gnmz(1,2,2) + Gnmz(1,1,1) * (A*(zA-zB)/(A+B))
        else
          do rloop_j = 1, rloop_i - 1
            Gnmz(1,rloop_i,rloop_j) = Gnmz(1,rloop_i,rloop_j) + &
            Gnmz(1,rloop_i-2,rloop_j) * real(rloop_i-2)/(2.0_dp*B) + &
            Gnmz(1,rloop_i-1,rloop_j) * (zB - zk)
            Gnmz(1,rloop_i,rloop_j + 1) = Gnmz(1,rloop_i,rloop_j + 1) - &
            Gnmz(1,rloop_i-2,rloop_j) * real(rloop_i-2)*A/(2.0_dp*B*(A+B)) + &
            Gnmz(1,rloop_i-1,rloop_j) * (A*(zA-zB)/(A+B))
          end do
        end if
      end do
      do rloop_k = 2, nzk + nzl + 1
        do rloop_i = 2, nzi + nzj + 1
          if (rloop_i == 2) then
            do rloop_j = 1, rloop_k
              Gnmz(2,rloop_k,rloop_j) = Gnmz(2,rloop_k,rloop_j) + &
              Gnmz(1,rloop_k,rloop_j) * (zA - zi)
              Gnmz(2,rloop_k,rloop_j + 1) = Gnmz(2,rloop_k,rloop_j + 1) + &
              Gnmz(1,rloop_k,rloop_j) * (B*(zB-zA)/(A+B)) + &
              Gnmz(1,rloop_k-1,rloop_j) * real(rloop_k-1)/(2.0_dp*(A+B))
            end do
          else
            do rloop_j = 1, rloop_i + rloop_k - 2
              Gnmz(rloop_i,rloop_k,rloop_j) = Gnmz(rloop_i,rloop_k,rloop_j) + &
              Gnmz(rloop_i-2,rloop_k,rloop_j) * real(rloop_i-2)/(2.0_dp*A) + &
              Gnmz(rloop_i-1,rloop_k,rloop_j) * (zA - zi)
              Gnmz(rloop_i,rloop_k,rloop_j+1) = Gnmz(rloop_i,rloop_k,rloop_j+1)&
              - Gnmz(rloop_i-2,rloop_k,rloop_j) * real(rloop_i-2)*B/(2.0_dp*A*&
              (A+B)) + Gnmz(rloop_i-1,rloop_k,rloop_j) * (B*(zB-zA)/(A+B)) + &
              Gnmz(rloop_i-1,rloop_k-1,rloop_j) * real(rloop_k-1)/(2.0_dp*(A+B))
            end do
          end if
        end do
      end do
      ! transfer from xi to xj, xk to xl
      Ixtrans = 0.0_dp
      Ix = 0.0_dp
      do rloop_k = 1, nxl + 1
        do rloop_i = 0, nxj
          do rloop_j = 1, nxi + nxj + 1 + nxk + nxl + 1
            Ixtrans(rloop_k,rloop_j) = Ixtrans(rloop_k,rloop_j) + &
            binomialcoe(nxj,rloop_i) * (xi-xj)**(rloop_i) * &
            Gnmx(nxi+nxj+1-rloop_i, nxk+nxl+1-(rloop_k-1),rloop_j)
          end do
        end do
      end do
      do rloop_i = 0, nxl
        do rloop_j = 1, nxi + nxj + 1 + nxk + nxl + 1
          Ix(rloop_j) = Ix(rloop_j) + binomialcoe(nxl,rloop_i) * &
          (xk-xl)**(rloop_i) * Ixtrans(rloop_i + 1,rloop_j)
        end do
      end do
      ! transfer from yi to yj, yk to yl
      Iytrans = 0.0_dp
      Iy = 0.0_dp
      do rloop_k = 1, nyl + 1
        do rloop_i = 0, nyj
          do rloop_j = 1, nyi + nyj + 1 + nyk + nyl + 1
            Iytrans(rloop_k,rloop_j) = Iytrans(rloop_k,rloop_j) + &
            binomialcoe(nyj,rloop_i) * (yi-yj)**(rloop_i) * &
            Gnmy(nyi+nyj+1-rloop_i, nyk+nyl+1-(rloop_k-1),rloop_j)
          end do
        end do
      end do
      do rloop_i = 0, nyl
        do rloop_j = 1, nyi + nyj + 1 + nyk + nyl + 1
          Iy(rloop_j) = Iy(rloop_j) + binomialcoe(nyl,rloop_i) * &
          (yk-yl)**(rloop_i) * Iytrans(rloop_i + 1,rloop_j)
        end do
      end do
      ! transfer from zi to zj, zk to zl
      Iztrans = 0.0_dp
      Iz = 0.0_dp
      do rloop_k = 1, nzl + 1
        do rloop_i = 0, nzj
          do rloop_j = 1, nzi + nzj + 1 + nzk + nzl + 1
            Iztrans(rloop_k,rloop_j) = Iztrans(rloop_k,rloop_j) + &
            binomialcoe(nzj,rloop_i) * (zi-zj)**(rloop_i) * &
            Gnmz(nzi+nzj+1-rloop_i, nzk+nzl+1-(rloop_k-1),rloop_j)
          end do
        end do
      end do
      do rloop_i = 0, nzl
        do rloop_j = 1, nzi + nzj + 1 + nzk + nzl + 1
          Iz(rloop_j) = Iz(rloop_j) + binomialcoe(nzl,rloop_i) * &
          (zk-zl)**(rloop_i) * Iztrans(rloop_i + 1,rloop_j)
        end do
      end do
      ! product of Ix, Iy, Iz
      PL = 0.0_dp
      !--------------------------------------------------------------
      do rloop_i = 0, nxi + nxj + 1 + nxk + nxl
        do rloop_j = 0, nyi + nyj + 1 + nyk + nyl
          do rloop_k = 0, nzi + nzj + 1 + nzk + nzl
            PL(rloop_i+rloop_j+rloop_k+1) = PL(rloop_i+rloop_j+rloop_k+1) + &
            Ix(rloop_i+1) * Iy(rloop_j+1) * Iz(rloop_k+1)
          end do
        end do
      end do
      PL = PL * 2.0_dp * dsqrt(rou/pi)
      ! integral of t exp(-X*t^2)*PL(t^2), 0 -> 1
      integral = 0.0_dp
      if (abs(X) < 1E-13) then 
        do rloop_i=1,nxi+nxj+1+nxk+nxl+1+nyi+nyj+1+nyk+nyl+1+nzi+nzj+1+nzk+nzl+1
          integral_mic = 1.0_dp / (real(2*rloop_i-2) + 1.0_dp)
          integral = integral + PL(rloop_i) * integral_mic
        end do
      else if (abs(X) <= xts) then
        rloop_i = 1
        do while(rloop_i <= &
        nxi+nxj+1+nxk+nxl+1+nyi+nyj+1+nyk+nyl+1+nzi+nzj+1+nzk+nzl+1)
          rloop_j = 0
          do while(rloop_j <= 2 * rloop_i - 2)
            if (rloop_j == 0) then ! dsqrt(pi/X) * erf(dsqrt(X)) / 2.0_dp
              do rloop_k=1,tayeps
                NX_mic(rloop_k)=(-1.0_dp)**(rloop_k)*&
                (1.0_dp/(real(2*rloop_k+1)*factorial(rloop_k)))
              end do
            else if(rloop_j == 1) then ! (1.0_dp - exp(-X)) / (2.0_dp * X)
              NX_mic_pre = NX_mic
              do rloop_k=1,tayeps
                NX_mic(rloop_k)=(-1.0_dp)**(rloop_k)*&
                (1.0_dp/(factorial(rloop_k+1)))
              end do
            else ! (-exp(-X)+real(rloop_j-1)*integral_mic_pre_pre) / (2.0_dp*X)
              NX_mic_pre_pre = NX_mic_pre
              NX_mic_pre = NX_mic
              do rloop_k=1,tayeps-1
                NX_mic(rloop_k)=((-1.0_dp)**(rloop_k)*&
                (0.5_dp/(factorial(rloop_k+1)))+&
                (NX_mic_pre_pre(rloop_k+1))/2.0_dp) * real(rloop_j+1)
              end do
              NX_mic(tayeps)=((-1.0_dp)**tayeps*&
              (0.5_dp/(factorial(tayeps+1)))) * real(rloop_j+1)
            end if
            rloop_j = rloop_j + 1
          end do
          if (2 * rloop_i - 2 == 0) then
            integral_mic = 0.0_dp
            do rloop_k=1,tayeps
              integral_mic = integral_mic + (-1.0_dp)**(rloop_k+1)*&
              X**(rloop_k-1)*(1.0_dp/(real(2*(rloop_k)-1)*factorial(rloop_k-1)))
            end do
          else if (2 * rloop_i - 2 == 1) then
            integral_mic = 0.0_dp
            do rloop_k=1,tayeps
              integral_mic = integral_mic + (-1.0_dp)**(rloop_k+1)*&
              X**(rloop_k-1)*(0.5_dp/factorial(rloop_k))
            end do
          else
            integral_mic = 0.0_dp
            do rloop_k=1,tayeps
              integral_mic = integral_mic + (-1.0_dp)**(rloop_k+1)*&
              X**(rloop_k-1)*(0.5_dp/factorial(rloop_k))&
              + X**(rloop_k-1)*NX_mic_pre_pre(rloop_k)/2.0_dp
            end do
          end if
          integral = integral + PL(rloop_i) * integral_mic
          rloop_i = rloop_i + 1
        end do
      else
        do rloop_i=1,nxi+nxj+1+nxk+nxl+1+nyi+nyj+1+nyk+nyl+1+nzi+nzj+1+nzk+nzl+1
          do rloop_j = 0, 2*rloop_i - 2
            if (rloop_j == 0) then
              integral_mic = dsqrt(pi/X) * erf(dsqrt(X)) / 2.0_dp
            else if(rloop_j == 1) then
              integral_mic_pre = integral_mic
              integral_mic = (1.0_dp - exp(-X)) / (2.0_dp * X)
            else
              integral_mic_pre_pre = integral_mic_pre
              integral_mic_pre = integral_mic
              integral_mic = (-exp(-X) + real(rloop_j - 1) * &
              integral_mic_pre_pre) / (2.0_dp * X)
            end if
          end do
          integral = integral + PL(rloop_i) * integral_mic
        end do
      end if
      V_Integral_2e = integral
      return
    end if
  end function V_Integral_2e
  
!-----------------------------------------------------------------------
!> diagonalisation of given real symmetric matrix
  subroutine diag(mat, dm, U, evl)
    implicit none
    integer,intent(in) :: dm
    real(dp),intent(in) :: mat(dm,dm)
    real(dp),intent(out) :: U(dm,dm)
    real(dp),intent(out) :: evl(dm)
    integer :: evl_count             ! number of eigenvalues found by dsyevr
    integer :: isupp_ev(2*dm)        ! indices indicating the nonzero elements
    integer :: lwork                 ! lwork for input of dsyevr
    integer :: liwork                ! liwork for input of dsyevr
    integer,allocatable :: iwork(:)  ! iwork for input of dsyevr
    real(dp) :: mat_u(dm,dm)         ! upper triangular part
    real(dp),allocatable :: work(:)  ! work for input of dsyevr
    integer :: dialoop_i, dialoop_j
    allocate(work(1))
    allocate(iwork(1))
    mat_u = 0.0_dp
    do dialoop_i = 1, dm
      do dialoop_j = dialoop_i, dm
        mat_u(dialoop_i,dialoop_j) = mat(dialoop_i,dialoop_j)
      end do
    end do
    call dsyevr('V','A','U',dm,mat_u,dm,0.0_dp,0.0_dp,0,0,safmin,&
    evl_count,evl,U,dm,isupp_ev,work,-1,iwork,-1,info)
    liwork = iwork(1)
    lwork = nint(work(1))
    deallocate(work)
    deallocate(iwork)
    allocate(work(lwork))
    allocate(iwork(liwork))
    call dsyevr('V','A','U',dm,mat_u,dm,0.0_dp,0.0_dp,0,0,safmin,&
    evl_count,evl,U,dm,isupp_ev,work,lwork,iwork,liwork,info)
    if(info < 0) then
      call terminate('illegal input of dsyevr')
    else if(info > 0) then
      call terminate('internal error of dsyevr')
    end if
    if (evl_count /= dm) call terminate('diagonalisation failed')
    deallocate(work)
    deallocate(iwork)
  end subroutine diag
  
!-----------------------------------------------------------------------
! symmetric orthogonalisation of given real symmetric matrix
! X mat XT = I
  
  subroutine symm_orth(mat, dm, X, min_evl)
    implicit none
    integer,intent(in) :: dm
    real(dp),intent(in) :: mat(dm,dm)
    real(dp),intent(out) :: X(:,:)
    real(dp),intent(out) :: min_evl
    real(dp) :: U(dm,dm)
    real(dp) :: supU(dm,dm)
    real(dp) :: evl(dm)
    real(dp) :: mat_u(dm,dm)
    integer :: symloop_i
    call diag(mat, dm, U, evl)
    min_evl = evl(1)
    mat_u = 0.0_dp
    do symloop_i = 1, dm
      if (evl(symloop_i) < 0.0) call terminate(&
      'symm_orth: eigenvalue less than zero, may due to code error')
      if (abs(evl(symloop_i)) < safmin) then
        write(60,'(A)') &
        '  --matrix is not invertible (not full rank), call can_orth directly'
        min_evl = 0.0_dp
        X = 0.0_dp
        return
      end if
      if (min_evl > evl(symloop_i)) min_evl = evl(symloop_i)
      mat_u(symloop_i,symloop_i) = 1.0_dp / dsqrt(evl(symloop_i))
    end do
    call dgemm('N', 'N', dm, dm, dm, 1.0_dp, U, dm, mat_u, dm, 0.0_dp, supU, dm)
    call dgemm('N', 'T', dm, dm, dm, 1.0_dp, supU, dm, U, dm, 0.0_dp, X, dm)
  end subroutine symm_orth
  
!-----------------------------------------------------------------------
!> canonical orthogonalisation of given real symmetric matrix
!!
!! X mat XT = I
  subroutine can_orth(mat, dm, X, dm2)
    implicit none
    integer,intent(in) :: dm
    real(dp),intent(in) :: mat(dm,dm)
    real(dp),allocatable,intent(out) :: X(:,:)
    integer,intent(out) :: dm2
    real(dp) :: U(dm,dm)
    real(dp) :: evl(dm)
    integer :: supU(dm)
    integer :: canloop_i, canloop_j, canloop_k
    call diag(mat, dm, U, evl)
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
  
!-----------------------------------------------------------------------
!> transfer matrix on Cartesian basis to matrix on spher-harmo basis
  subroutine sphehar(m)
    implicit none
    real(dp),allocatable :: m(:, :) ! target matrix
    real(dp),allocatable :: sphoper(:,:)
    if (.not. allocated(c2s)) &
    call terminate('call sphehar error, c2s not allocated')
    if (.not. allocated(m)) &
    call terminate('call sphehar error, matrix not allocated')
    allocate(sphoper(cbdm,sbdm))
    call dgemm('N', 'N', cbdm, sbdm, cbdm, 1.0_dp, m, cbdm, c2s, &
    cbdm, 0.0_dp, sphoper, cbdm)
    deallocate(m)
    allocate(m(sbdm,sbdm))
    call dgemm('T', 'N', sbdm, sbdm, cbdm, 1.0_dp, c2s, cbdm, sphoper, &
    cbdm, 0.0_dp, m, sbdm)
    deallocate(sphoper)
  end subroutine sphehar
  
!-----------------------------------------------------------------------
!> transfer matrix on spher-harmo basis to elimination linear dependence
  subroutine detach(m)
    implicit none
    real(dp),allocatable :: m(:, :) ! target matrix
    real(dp),allocatable :: detoper(:,:)
    if (.not. allocated(Xm)) &
    call terminate('call detach error, Xm not allocated')
    if (.not. allocated(m)) &
    call terminate('call detach error, matrix not allocated')
    allocate(detoper(fbdm, sbdm))
    call dgemm( 'T', 'N', fbdm, sbdm, sbdm, 1.0_dp, Xm, sbdm, m, &
    sbdm, 0.0_dp, detoper, fbdm)
    deallocate(m)
    allocate(m(fbdm, fbdm))
    call dgemm( 'N', 'N', fbdm, fbdm, sbdm, 1.0_dp, detoper, fbdm, Xm, &
    sbdm, 0.0_dp, m, fbdm)
    deallocate(detoper)
  end subroutine detach
  
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
    integer :: bloop_i, numerator, denominator
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
      bloop_i = N
      do while(bloop_i > M)
        numerator = numerator * bloop_i
        bloop_i = bloop_i - 1
      end do
      denominator = factorial(N - M)
      bicoe = real(numerator) / real(denominator)
      return
    end if
  end function binomialcoe
  
end module Hamiltonian