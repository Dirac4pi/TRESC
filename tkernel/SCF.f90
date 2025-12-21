!> @file SCF.f90
!!
!! @brief single-configuration self-consistent-field procedure
!!
!! @syntax Fortran 2008 free format
!!
!! @code UTF-8
!!
!! @author dirac4pi
module SCF
  use Hamiltonian
  use Atoms
  use Fundamentals
  use Representation
  
  ! spinor MO coefficients, in order of (AO1,0), (0,AO1), ... (AOn,0), (0,AOn)
  !DIR$ ATTRIBUTES ALIGN:align_size :: AO2MO, AO2MOalpha, AO2MObeta
  !DIR$ ATTRIBUTES ALIGN:align_size :: rho_m, AOsupp, Fock
  complex(dp),allocatable :: AO2MO(:,:)     ! MO coeff
  integer,allocatable     :: occindex(:)    ! occupied orbital number of MO
  real(dp),allocatable    :: M_AO2MO_a(:,:) ! alpha MO coeff read form MOLDEN
  real(dp),allocatable    :: T_AO2MO_a(:,:) ! transfered alpha MO coeff
  real(dp),allocatable    :: M_AO2MO_b(:,:) ! beta MO coeff read form MOLDEN
  real(dp),allocatable    :: T_AO2MO_b(:,:) ! transfered beta MO coeff
  complex(dp),allocatable :: rho_m(:,:)     ! density matrix, complex Hermitian
  real(dp)                :: RMSDP, maxDP
  real(dp),allocatable    :: AOsupp(:,:)
  complex(dp),allocatable :: Fock(:,:)      ! Fock matrix
  integer                 :: iter           ! SCF iteration number
  
  logical                 :: ini_rho =.true.! initial density matrix loaded
  
  !--------------<one-electron Fock>--------------
  !DIR$ ATTRIBUTES ALIGN:align_size :: Fock1, oper1, oper2, oper3
  !DIR$ ATTRIBUTES ALIGN:align_size :: oper4, oper5, oper6, Ve
  !DIR$ ATTRIBUTES ALIGN:align_size :: pxVepx, pyVepy, pzVepz, pxVepy
  !DIR$ ATTRIBUTES ALIGN:align_size :: pyVepx, pxVepz, pzVepx, pyVepz
  !DIR$ ATTRIBUTES ALIGN:align_size :: pzVepy, Ap, ApRp, SRp, ARVRA
  !DIR$ ATTRIBUTES ALIGN:align_size :: ARVeRA, AVA, AVeA, exAO2p2
  complex(dp),allocatable :: Fock1(:,:)     ! one-electron Fock matrix
  real(dp),allocatable    :: oper1(:,:)     ! operator matrix
  real(dp),allocatable    :: oper2(:,:)     ! operator matrix
  complex(dp),allocatable :: oper3(:,:)     ! operator matrix
  complex(dp),allocatable :: oper4(:,:)     ! operator matrix
  complex(dp),allocatable :: oper5(:,:)     ! operator matrix
  complex(dp),allocatable :: oper6(:,:)     ! operator matrix
  real(dp),allocatable    :: Ve(:,:)        ! V/(E+E')
  real(dp),allocatable    :: pxVepx(:,:)    ! pxVpx/(E+E')
  real(dp),allocatable    :: pyVepy(:,:)    ! pyVpy/(E+E')
  real(dp),allocatable    :: pzVepz(:,:)    ! pzVpz/(E+E')
  real(dp),allocatable    :: pxVepy(:,:)    ! pxVpy/(E+E')
  real(dp),allocatable    :: pyVepx(:,:)    ! pyVpx/(E+E')
  real(dp),allocatable    :: pxVepz(:,:)    ! pxVpz/(E+E')
  real(dp),allocatable    :: pzVepx(:,:)    ! pzVpx/(E+E')
  real(dp),allocatable    :: pyVepz(:,:)    ! pyVpz/(E+E')
  real(dp),allocatable    :: pzVepy(:,:)    ! pzVpy/(E+E')
  real(dp),allocatable    :: Ap(:,:)        ! Ap
  real(dp),allocatable    :: ApRp(:,:)      ! ApRp
  complex(dp),allocatable :: SRp(:,:)       ! ApRp(1+p^2/4c^2)
  complex(dp),allocatable :: ARVRA(:,:)     ! ApRpVRpAp
  complex(dp),allocatable :: ARVeRA(:,:)    ! ApRp(V/(E+E'))RpAp
  complex(dp),allocatable :: AVA(:,:)       ! ApVAp
  complex(dp),allocatable :: AVeA(:,:)      ! Ap(V/(E+E'))Ap
  complex(dp),allocatable :: exAO2p2(:,:)   ! extended AO2p2 matrix
  
  !--------------<two-electron Fock>--------------
  !DIR$ ATTRIBUTES ALIGN:align_size :: mHFcol, mHFexc, mKSexc, mKScor,mpVpcol_11
  !DIR$ ATTRIBUTES ALIGN:align_size :: mpVpexc_11, mpVpcol_22, mpVpexc_22
  !DIR$ ATTRIBUTES ALIGN:align_size :: iijj_V, iijj_pxVpx,iijj_pyVpy, iijj_pzVpz
  !DIR$ ATTRIBUTES ALIGN:align_size :: ijij_pxVpx, ijij_pyVpy, ijij_pzVpz
  complex(dp),allocatable :: mHFcol(:,:)    ! 2e HF Coulomb matrix
  complex(dp),allocatable :: mHFexc(:,:)    ! 2e HF Exchange matrix
  complex(dp),allocatable :: mKSexc(:,:)    ! 2e KS Exchange matrix
  complex(dp),allocatable :: mKScor(:,:)    ! 2e KS correlation matrix
  complex(dp),allocatable :: mpVpcol_11(:,:)! 2e pVp Coulomb matrix(one with ii)
  complex(dp),allocatable :: mpVpexc_11(:,:)! 2e pVp Exchange matrix
  complex(dp),allocatable :: mpVpcol_22(:,:)! 2e pVp Coulomb matrix(one with ii)
  complex(dp),allocatable :: mpVpexc_22(:,:)! 2e pVp Exchange matrix
  real(dp),allocatable    :: iijj_V(:,:)    ! <ij|V|ij>
  real(dp),allocatable    :: iijj_pxVpx(:,:)! <ij|pxVpx|ij> = (pxipxi|V|jj)
  real(dp),allocatable    :: iijj_pyVpy(:,:)! <ij|pyVpy|ij> = (pyipyi|V|jj)
  real(dp),allocatable    :: iijj_pzVpz(:,:)! <ij|pzVpz|ij> = (pzipzi|V|jj)
  real(dp),allocatable    :: ijij_pxVpx(:,:)! <ii|pxVpx|jj> = (pxij|V|pxij)
  real(dp),allocatable    :: ijij_pyVpy(:,:)! <ii|pyVpy|jj> = (pyij|V|pyij)
  real(dp),allocatable    :: ijij_pzVpz(:,:)! <ii|pzVpz|jj> = (pzij|V|pzij)
  
  ! orbital energy and mol energy
  !DIR$ ATTRIBUTES ALIGN:align_size :: orbE
  real(dp),allocatable    :: orbE(:)        ! orbital energy
  real(dp)                :: molE_pre, molE ! mol energy
  real(dp)                :: Virial         ! Virial ratio
  real(dp)                :: nucE           ! nuclear repulsion energy
  real(dp)                :: HFCol          ! Hartree-Fock Coulomb energy
  real(dp)                :: HFexc          ! Hartree-Fock exchange energy
  real(dp)                :: KSexc          ! Kohn-Shanm exchange energy
  real(dp)                :: KScor          ! Kohn-Shanm correlation energy
  real(dp)                :: Ecore          ! one-electron (core) energy
  real(dp)                :: E2e            ! two-electron energy
  real(dp)                :: T              ! kinetic energy
  real(dp)                :: V              ! electron-nuclear attraction energy
  real(dp)                :: EpVp           ! 1e pVp-related energy
  real(dp)                :: ESR            ! SRTP and radiation energy
  real(dp)                :: EpVpcol        ! 2e pvp-related Coulomb energy
  real(dp)                :: EpVpexc        ! 2e pvp-related Exchange energy
  real(dp)                :: emd4           ! dispersion energy calc by DFT-D4
  real(dp)                :: scf_kappa      ! deviation parameter from TRS

  !DIR$ ATTRIBUTES ALIGN:align_size :: rho_pre, rho_history, Rsd
  ! damping & DIIS(AX=B)
  ! privious rho_m of current iteration
  complex(dp),allocatable :: rho_pre(:,:)
  complex(dp),allocatable :: rho_history(:,:,:)! coeff of subsp iteration
  complex(dp),allocatable :: Rsd(:,:,:)        ! residuals of rho_m
  complex(dp),allocatable :: DIISmat(:,:)      ! A
  real(dp)                :: damp_coe          ! damp coeff of direct/DIIS SCF
  
  contains

!------------------------------------------------------------
!> HF/KS-SCF procedure for 2e scalar V interaction
!!
!! call-up condition: pVp1e = .False. pVp2e = .False.
  subroutine SCF_V2e()
    implicit none
    integer             :: ii, jj, kk, ll, mm   ! loop variable SCF_V2e
    write(60,'(A)') 'Module SCF:'
    write(60,'(A)') '  construct one-electron Fock matrix'
    call Assign_Fock_1e()
    write(60,'(A)') '  complete! stored in Fock1'
    write(60,'(A)') '  calculate nuclear repulsion energy'
    nucE = 0.0_dp
    do ii = 1, atom_count
      do jj = ii + 1, atom_count
        nucE = nucE + (real(mol(ii) % atom_number) * &
        real(mol(jj) % atom_number)) / dsqrt(&
        (mol(ii) % pos(1) - mol(jj) % pos(1))**2 + &
        (mol(ii) % pos(2) - mol(jj) % pos(2))**2 + &
        (mol(ii) % pos(3) - mol(jj) % pos(3))**2)
      end do
    end do
    if (nucE >= 1E12) call terminate(&
    'nuclear repulsive energy anomaly, may due to overlap atomic coordinates')
    write(60,'(A,F12.7,A)') &
    '  complete! nuclear repulsive energy = ', nucE, ' Eh'
    call Assign_Schwarz_V2e()
    write(60,'(A)') '  initialize Libxc'
    if (fx_id /= -1) call Fockxc_init()
    write(60,'(A)') '  complete!'
    allocate(Fock(2*fbdm,2*fbdm))
    allocate(orbE(2*fbdm))
    allocate(oper6(2*fbdm,2*fbdm))
    allocate(oper3(2*fbdm,2*fbdm))
    allocate(oper4(2*fbdm,2*fbdm))
    ! ----------<DIIS>----------
    ! rho_m(new) = sigma(i,subsp) DIIScoe(i)*(rho_history(i)+diisdamp*Rsd(i))
    allocate(Rsd(subsp,2*sbdm,2*sbdm))
    allocate(DIISmat(subsp+1,subsp+1))
    allocate(rho_history(subsp, 2*sbdm, 2*sbdm))
    allocate(rho_pre(2*sbdm, 2*sbdm))
    ! ----------<DIIS>----------
    write(60,'(A)') '  ----------<SCF>----------'
    write(60,'(A)') '  SCF settings:'
    write(60,'(A,I3.3)')  '  -- maxiter   = ', maxiter
    write(60,'(A,E10.3)') '  -- conv_tol  =',  conver_tol
    write(60,'(A,F6.3)')  '  -- damp      =',  damp
    write(60,'(A,E10.3)') '  -- cutdamp   =',  cutdamp
    write(60,'(A,I3.3)')  '  -- nodiis    = ', nodiis
    write(60,'(A,I3.3)')  '  -- subsp     = ', subsp
    write(60,'(A,F6.3)')  '  -- diisdamp  =',  diisdamp
    write(60,'(A,E10.3)') '  -- cutdiis   =',  cutdiis
    do iter = 1, maxiter
      write(60,*)
      write(60,*)
      write(60,'(A,I3)') '  SCF iter ',iter
      if (iter == 1) then
        write(60,'(A)') '  load density matrix'
        call Assign_rho()
        write(60,'(A)') '  complete! stored in rho_m'
      end if
      ! construct 2e Fock matrices
      write(60,'(A)') '  construct scalar 2e Fock matrices'
      call Assign_Fock_V2e()
      write(60,'(A)') '  complete! stored in mHFcol, mHFexc'
      if (fx_id /= -1) then
        if (fc_id /= -1) then
          write(60,'(A)') '  mKSexc, mKScor'
          Fock = Fock1 + mHFcol + x_HF*mHFexc + (1.0_dp-x_HF)*mKSexc + mKScor
        else
          x_HF = 0.0_dp
          KScor = 0.0_dp
          write(60,'(A)') '  mKSexc'
          Fock = Fock1 + mHFcol + mKSexc
        end if
      else    ! pure Hartree-Fock
        x_HF = 1.0_dp
        KScor = 0.0_dp
        KSexc = 0.0_dp
        Fock = Fock1 + mHFcol + mHFexc
      end if
      ! diagonalization of Fock matrix
      write(60,'(A)') '  diagonalization of Fock matrix'
      call diag(Fock, 2*fbdm, oper3, orbE)
      write(60,'(A,I5,A)') '  complete!',2*fbdm,' eigenvectors found'
      ! de-orthogonalization
      write(60,'(A)') '  generate and dump AO2MO to .ao2mo file'
      call matmul('N', 'N', exs2f, oper3, AO2MO)
      call Dump_matrix('ao2mo', AO2MO, 2*sbdm, 2*fbdm)
      ! construct new density matrix
      write(60,'(A)') '  construct new density matrix'
      call Assign_rho()
      write(60,'(A)') '  complete! stored in rho_m'

      call SCFprint()
      ! convergence check
      if (iter == 1) then
        write(60,'(A,F12.6)') '  SCF energy (A.U.) = ',molE
      else
        write(60,'(A,F12.6,A,E10.3)') '  SCF energy (A.U.) = ',molE,&
        '; delta E = ', molE-molE_pre
        if (abs(molE-molE_pre) < conver_tol   .and. &
        (abs(RMSDP) < 0.5*abs(molE-molE_pre)  .or.  &
        ! orbital degeneracy induces density matrix oscillations
        abs(RMSDP) > 10.0*abs(molE-molE_pre)) .and. &
        (damp_coe < 0.01                      .or.  &
        damp_coe < 0.1*(log10(conver_tol)-log10(abs(molE-molE_pre))))) then
          write(60,'(A)') '  convergence tolerance met, SCF done!'
          exit
        else
          write(60,'(A)') '  convergence tolerance not met'
        end if
      end if
      write(60,'(A)') '  DIIS information'
      call DIIS()
      rho_pre = rho_m
    end do
    if (abs(molE-molE_pre) < conver_tol .and. iter < maxiter) then
      write(60,'(A)') '  SCF succeed!'
    else
      write(60,'(A)') '  SCF failed!'
    end if
    if (d4) then
      emd4 = DFTD4()
      molE = molE + emd4
    end if
  end subroutine SCF_V2e

!------------------------------------------------------------
!> HF/KS-SCF procedure for 2e scalar AVA interaction
!!
!! call-up condition: pVp1e = .True. pVp2e = .False.
  subroutine SCF_AVA2e()
    implicit none
    integer             :: ii, jj, kk, ll, mm   ! loop variable SCF_AVA2e
    write(60,'(A)') 'Module SCF:'
    write(60,'(A)') '  construct one-electron Fock matrix'
    call Assign_Fock_1e()
    write(60,'(A)') '  complete! stored in Fock1'
    write(60,'(A)') '  calculate nuclear repulsion energy'
    nucE = 0.0_dp
    do ii = 1, atom_count
      do jj = ii + 1, atom_count
        nucE = nucE + (real(mol(ii) % atom_number) * &
        real(mol(jj) % atom_number)) / dsqrt(&
        (mol(ii) % pos(1) - mol(jj) % pos(1))**2 + &
        (mol(ii) % pos(2) - mol(jj) % pos(2))**2 + &
        (mol(ii) % pos(3) - mol(jj) % pos(3))**2)
      end do
    end do
    if (nucE >= 1E12) call terminate(&
    'nuclear repulsive energy anomaly, may due to overlap atomic coordinates')
    write(60,'(A,F12.7,A)') &
    '  complete! nuclear repulsive energy = ', nucE, ' Eh'
    call Assign_Schwarz_V2e()
    write(60,'(A)') '  initialize Libxc'
    if (fx_id /= -1) call Fockxc_init()
    write(60,'(A)') '  complete!'
    allocate(Fock(2*fbdm,2*fbdm))
    allocate(orbE(2*fbdm))
    allocate(oper6(2*fbdm,2*fbdm))
    allocate(oper3(2*fbdm,2*fbdm))
    allocate(oper4(2*fbdm,2*fbdm))
    ! ----------<DIIS>----------
    ! rho_m(new) = sigma(i,subsp) DIIScoe(i)*(rho_history(i)+diisdamp*Rsd(i))
    allocate(Rsd(subsp,2*sbdm,2*sbdm))
    allocate(DIISmat(subsp+1,subsp+1))
    allocate(rho_history(subsp, 2*sbdm, 2*sbdm))
    allocate(rho_pre(2*sbdm, 2*sbdm))
    ! ----------<DIIS>----------
    write(60,'(A)') '  ----------<SCF>----------'
    write(60,'(A)') '  SCF settings:'
    write(60,'(A,I3.3)')  '  -- maxiter   = ', maxiter
    write(60,'(A,E10.3)') '  -- conv_tol  =',  conver_tol
    write(60,'(A,F6.3)')  '  -- damp      =',  damp
    write(60,'(A,E10.3)') '  -- cutdamp   =',  cutdamp
    write(60,'(A,I3.3)')  '  -- nodiis    = ', nodiis
    write(60,'(A,I3.3)')  '  -- subsp     = ', subsp
    write(60,'(A,F6.3)')  '  -- diisdamp  =',  diisdamp
    write(60,'(A,E10.3)') '  -- cutdiis   =',  cutdiis
    do iter = 1, maxiter
      write(60,*)
      write(60,*)
      write(60,'(A,I3)') '  SCF iter ',iter
      if (iter == 1) then
        write(60,'(A)') '  load density matrix'
        call Assign_rho()
        write(60,'(A)') '  complete! stored in rho_m'
      end if
      ! construct 2e Fock matrices
      write(60,'(A)') '  construct relativistic scalar 2e Fock matrices'
      call Assign_Fock_AVA2e()
      write(60,'(A)') '  complete! stored in mHFcol, mHFexc'
      if (fx_id /= -1) then
        if (fc_id /= -1) then
          write(60,'(A)') '  mKSexc, mKScor'
          Fock = Fock1 + mHFcol + x_HF*mHFexc + (1.0_dp-x_HF)*mKSexc + mKScor
        else
          x_HF = 0.0_dp
          KScor = 0.0_dp
          write(60,'(A)') '  mKSexc'
          Fock = Fock1 + mHFcol + mKSexc
        end if
      else    ! pure Hartree-Fock
        x_HF = 1.0_dp
        KScor = 0.0_dp
        KSexc = 0.0_dp
        Fock = Fock1 + mHFcol + mHFexc
      end if

      ! diagonalization of Fock matrix
      write(60,'(A)') '  diagonalization of Fock matrix'
      call diag(Fock, 2*fbdm, oper3, orbE)
      write(60,'(A,I5,A)') '  complete!',2*fbdm,' eigenvectors found'
      ! de-orthogonalization
      write(60,'(A)') '  generate and dump AO2MO to .ao2mo file'
      call matmul('N', 'N', exs2f, oper3, AO2MO)
      call Dump_matrix('ao2mo', AO2MO, 2*sbdm, 2*fbdm)
      ! construct new density matrix
      write(60,'(A)') '  construct new density matrix'
      call Assign_rho()
      write(60,'(A)') '  complete! stored in rho_m'
      call SCFprint()
      ! convergence check
      if (iter == 1) then
        write(60,'(A,F12.6)') '  SCF energy (A.U.) = ',molE
      else
        write(60,'(A,F12.6,A,E10.3)') '  SCF energy (A.U.) = ',molE,&
        '; delta E = ', molE-molE_pre
        if (abs(molE-molE_pre) < conver_tol   .and. &
        (abs(RMSDP) < 0.5*abs(molE-molE_pre)  .or.  &
        ! orbital degeneracy induces density matrix oscillations
        abs(RMSDP) > 10.0*abs(molE-molE_pre)) .and. &
        (damp_coe < 0.01                      .or.  &
        damp_coe < 0.1*(log10(conver_tol)-log10(abs(molE-molE_pre))))) then
          write(60,'(A)') '  convergence tolerance met, SCF done!'
          exit
        else
          write(60,'(A)') '  convergence tolerance not met'
        end if
      end if
      write(60,'(A)') '  DIIS information'
      call DIIS()
      rho_pre = rho_m
    end do
    if (abs(molE-molE_pre) < conver_tol .and. iter < maxiter) then
      write(60,'(A)') '  SCF succeed!'
    else
      write(60,'(A)') '  SCF failed!'
    end if
    if (d4) then
      emd4 = DFTD4()
      molE = molE + emd4
    end if
  end subroutine SCF_AVA2e

!------------------------------------------------------------
!> HF/KS-SCF procedure for 2e spinor AVA and ARVRA interaction
!!
!! call-up condition: pVp1e = .True. pVp2e = .True.
  subroutine SCF_ARVRA2e()
    implicit none
    integer             :: ii, jj, kk, ll, mm   ! loop variable SCF_ARVRA2e
    write(60,'(A)') 'Module SCF:'
    write(60,'(A)') '  construct one-electron Fock matrix'
    call Assign_Fock_1e()
    write(60,'(A)') '  complete! stored in Fock1'
    write(60,'(A)') '  calculate nuclear repulsion energy'
    nucE = 0.0_dp
    do ii = 1, atom_count
      do jj = ii + 1, atom_count
        nucE = nucE + (real(mol(ii) % atom_number) * &
        real(mol(jj) % atom_number)) / dsqrt(&
        (mol(ii) % pos(1) - mol(jj) % pos(1))**2 + &
        (mol(ii) % pos(2) - mol(jj) % pos(2))**2 + &
        (mol(ii) % pos(3) - mol(jj) % pos(3))**2)
      end do
    end do
    if (nucE >= 1E12) call terminate(&
    'nuclear repulsive energy anomaly, may due to overlap atomic coordinates')
    write(60,'(A,F12.7,A)') &
    '  complete! nuclear repulsive energy = ', nucE, ' Eh'
    call Assign_Schwarz_pVp2e()
    write(60,'(A)') '  initialize Libxc'
    if (fx_id /= -1) call Fockxc_init()
    write(60,'(A)') '  complete!'
    allocate(Fock(2*fbdm,2*fbdm))
    allocate(orbE(2*fbdm))
    allocate(oper6(2*fbdm,2*fbdm))
    allocate(oper3(2*fbdm,2*fbdm))
    allocate(oper4(2*fbdm,2*fbdm))
    ! ----------<DIIS>----------
    ! rho_m(new) = sigma(i,subsp) DIIScoe(i)*(rho_history(i)+diisdamp*Rsd(i))
    allocate(Rsd(subsp,2*sbdm,2*sbdm))
    allocate(DIISmat(subsp+1,subsp+1))
    allocate(rho_history(subsp, 2*sbdm, 2*sbdm))
    allocate(rho_pre(2*sbdm, 2*sbdm))
    ! ----------<DIIS>----------
    write(60,'(A)') '  ===============<SCF>==============='
    pVp2e = .False.
    damp_ = damp
    diisdamp_ = diisdamp
    write(60,*)
    write(60,*)
    write(60,*)
    write(60,*)
    write(60,'(A)') '        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(60,'(A)') '        !!!!!!<SCALAR CONVERGENCE STAGE>!!!!!!'
    write(60,'(A)') '        !!SCF settings:                     !!'
    write(60,'(A,I3.3,A)')  '        !!-- maxiter   = ', maxiter, &
    '                !!'
    write(60,'(A,E10.3,A)') '        !!-- conv_tol  =',  conver_tol, &
    '          !!'
    write(60,'(A,F6.3,A)')  '        !!-- damp      =',  damp, &
    '              !!'
    write(60,'(A,E10.3,A)') '        !!-- cutdamp   =',  cutdamp, &
    '          !!'
    write(60,'(A,I3.3,A)')  '        !!-- nodiis    = ', nodiis, &
    '                !!'
    write(60,'(A,I3.3,A)')  '        !!-- subsp     = ', subsp, &
    '                !!'
    write(60,'(A,F6.3,A)')  '        !!-- diisdamp  =',  diisdamp, &
    '              !!'
    write(60,'(A,E10.3,A)') '        !!-- cutdiis   =',  cutdiis, &
    '          !!'
    write(60,'(A)') '        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(60,'(A)') '        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(60,*)
    write(60,*)
    write(60,*)
    do iter = 1, maxiter
      write(60,*)
      write(60,*)
      write(60,'(A,I3)') '  scalar SCF iter ',iter
      if (iter == 1) then
        write(60,'(A)') '  load density matrix'
        call Assign_rho()
        write(60,'(A)') '  complete! stored in rho_m'
      end if
      ! construct 2e Fock matrices
      write(60,'(A)') '  construct relativistic scalar 2e Fock matrices'
      call Assign_Fock_AVA2e()
      write(60,'(A)') '  complete! stored in mHFcol, mHFexc'
      if (fx_id /= -1) then
        if (fc_id /= -1) then
          write(60,'(A)') '  mKSexc, mKScor'
          Fock = Fock1 + mHFcol + x_HF*mHFexc + (1.0_dp-x_HF)*mKSexc + mKScor
        else
          x_HF = 0.0_dp
          KScor = 0.0_dp
          write(60,'(A)') '  mKSexc'
          Fock = Fock1 + mHFcol + mKSexc
        end if
      else    ! pure Hartree-Fock
        x_HF = 1.0_dp
        KScor = 0.0_dp
        KSexc = 0.0_dp
        Fock = Fock1 + mHFcol + mHFexc
      end if
      ! diagonalization of Fock matrix
      write(60,'(A)') '  diagonalization of Fock matrix'
      call diag(Fock, 2*fbdm, oper3, orbE)
      write(60,'(A,I5,A)') '  complete!',2*fbdm,' eigenvectors found'
      ! de-orthogonalization
      write(60,'(A)') '  generate and dump AO2MO to .ao2mo file'
      call matmul('N', 'N', exs2f, oper3, AO2MO)
      call Dump_matrix('ao2mo', AO2MO, 2*sbdm, 2*fbdm)
      ! construct new density matrix
      write(60,'(A)') '  construct new density matrix'
      call Assign_rho()
      write(60,'(A)') '  complete! stored in rho_m'
      call SCFprint()
      ! convergence check
      if (iter == 1) then
        write(60,'(A,F12.6)') '  SCF energy (A.U.) = ',molE
      else
        write(60,'(A,F12.6,A,E10.3)') '  SCF energy (A.U.) = ',molE,&
        '; delta E = ', molE-molE_pre
        if (abs(molE-molE_pre) < conver_tol   .and. &
        (abs(RMSDP) < 0.5*abs(molE-molE_pre)  .or.  &
        ! orbital degeneracy induces density matrix oscillations
        abs(RMSDP) > 10.0*abs(molE-molE_pre)) .and. &
        (damp_coe < 0.01                      .or.  &
        damp_coe < 0.1*(log10(conver_tol)-log10(abs(molE-molE_pre))))) then
          write(60,'(A)') '  convergence tolerance met, SCF done!'
          exit
        else
          write(60,'(A)') '  convergence tolerance not met'
        end if
      end if
      write(60,'(A)') '  DIIS information'
      call DIIS()
      rho_pre = rho_m
    end do
    if (abs(molE-molE_pre) < conver_tol .and. iter+1 < maxiter) then
      write(60,'(A)') '  SCF succeed!'
    else
      write(60,'(A)') '  SCF failed! Skip spinor convergence stage!'
      return
    end if
    write(60,*)
    write(60,*)
    write(60,*)
    write(60,*)
    write(60,'(A)') '        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(60,'(A)') '        !!!!!!<SPINOR CONVERGENCE STAGE>!!!!!!'
    pVp2e = .True.
    if (subsp <= 8) then
      nodiis = subsp
    else
      nodiis = 5
      subsp = 5
    end if
    write(60,'(A)') '        !!SCF settings:                     !!'
    write(60,'(A,I3.3,A)')  '        !!-- maxiter  =            ', maxiter,   &
    '      !!'
    write(60,'(A,E10.3,A)') '        !!-- conv_tol =           ',  conver_tol,&
    '!!'
    write(60,'(A,F6.3,A)')  '        !!-- damp     = dynamic to',  damp,      &
    '    !!'
    write(60,'(A,E10.3,A)') '        !!-- cutdamp  =           ',  cutdamp,   &
    '!!'
    write(60,'(A,I3.3,A)')  '        !!-- nodiis   =            ', nodiis,    &
    '      !!'
    write(60,'(A,I3.3,A)')  '        !!-- subsp    =            ', subsp,     &
    '      !!'
    write(60,'(A,F6.3,A)')  '        !!-- diisdamp = dynamic to',  diisdamp,  &
    '    !!'
    write(60,'(A,E10.3,A)') '        !!-- cutdiis  =           ',  cutdiis,   &
    '!!'
    write(60,'(A)') '        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(60,'(A)') '        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(60,*)
    write(60,*)
    write(60,*)
    do iter = 1, maxiter
      write(60,*)
      write(60,*)
      write(60,'(A,I3)') '  spinor SCF iter ', iter
      if (iter == 1) then
        write(60,'(A)') '  load density matrix from scalar convergence stage'
        call Assign_rho()
        write(60,'(A)') '  complete! stored in rho_m'
      end if
      ! construct 2e Fock matrices
      write(60,'(A)') '  construct relativistic spinor 2e Fock matrices'
      call Assign_Fock_ARVRA2e()
      write(60,'(A)') '  complete! stored in mHFcol, mHFexc'
      write(60,'(A)') '  mpVpcol_11, mpVpexc_11, mpVpcol_22, mpVpexc_22'
      ! assign Fock matrix
      if (fx_id /= -1) then
        if (fc_id /= -1) then
          write(60,'(A)') '  mKSexc, mKScor'
          Fock = Fock1 + mHFcol + x_HF*mHFexc + (1.0_dp-x_HF)*mKSexc + &
          mKScor + mpVpcol_11 + mpVpexc_11 + mpVpcol_22 + mpVpexc_22
        else
          x_HF = 0.0_dp
          KScor = 0.0_dp
          write(60,'(A)') '  mKSexc'
          Fock = Fock1 + mHFcol + mKSexc + &
          mpVpcol_11 + mpVpexc_11 + mpVpcol_22 + mpVpexc_22
        end if
      else    ! pure Hartree-Fock
        x_HF = 1.0_dp
        KScor = 0.0_dp
        KSexc = 0.0_dp
        Fock = Fock1 + mHFcol + mHFexc + &
        mpVpcol_11 + mpVpexc_11 + mpVpcol_22 + mpVpexc_22
      end if
      ! diagonalization of Fock matrix
      write(60,'(A)') '  diagonalization of Fock matrix'
      call diag(Fock, 2*fbdm, oper3, orbE)
      write(60,'(A,I5,A)') '  complete!',2*fbdm,' eigenvectors found'
      ! de-orthogonalization
      write(60,'(A)') '  generate and dump AO2MO to .ao2mo file'
      call matmul('N', 'N', exs2f, oper3, AO2MO)
      call Dump_matrix('ao2mo', AO2MO, 2*sbdm, 2*fbdm)
      ! construct new density matrix
      write(60,'(A)') '  construct new density matrix'
      call Assign_rho()
      write(60,'(A)') '  complete! stored in rho_m'
      call SCFprint()
      ! convergence check
      if (iter == 1) then
        write(60,'(A,F12.6)') '  SCF energy (A.U.) = ',molE
      else
        write(60,'(A,F12.6,A,E10.3)') '  SCF energy (A.U.) = ',molE,&
        '; delta E = ', molE-molE_pre
        if (abs(molE-molE_pre) < conver_tol   .and. &
        (abs(RMSDP) < 0.5*abs(molE-molE_pre)  .or.  &
        ! orbital degeneracy induces density matrix oscillations
        abs(RMSDP) > 10.0*abs(molE-molE_pre)) .and. &
        (damp_coe < 0.01                      .or.  &
        damp_coe < 0.1*(log10(conver_tol)-log10(abs(molE-molE_pre))))) then
          write(60,'(A)') '  convergence tolerance met, SCF done!'
          exit
        else
          write(60,'(A)') '  convergence tolerance not met'
        end if
      end if
      write(60,'(A)') '  DIIS information'
      ! dynamic damping control
      if (iter /= 1) then
        if (abs(molE - molE_pre) > 1E-3) then
          if (damp_ < 0.7) damp = 0.7_dp
          if (diisdamp_ < 0.7) diisdamp = 0.7_dp
        else if (abs(molE - molE_pre) > 1E-4) then
          if (damp_ < 0.5) damp = 0.5_dp
          if (diisdamp_ < 0.5) diisdamp = 0.5_dp
        else if (abs(molE - molE_pre) > 1E-5) then
          if (damp_ < 0.3) damp = 0.3_dp
          if (diisdamp_ < 0.3) diisdamp = 0.3_dp
        else
          damp = damp_
          diisdamp = diisdamp_
        end if
      else
        if (damp_ < 0.7) damp = 0.7_dp
        if (diisdamp_ < 0.7) diisdamp = 0.7_dp
      end if
      call DIIS()
      rho_pre = rho_m
    end do
    if (abs(molE-molE_pre) < conver_tol .and. iter < maxiter) then
      write(60,'(A)') '  SCF succeed!'
    else
      write(60,'(A)') '  SCF failed!'
    end if
    if (d4) then
      emd4 = DFTD4()
      molE = molE + emd4
    end if
  end subroutine SCF_ARVRA2e

!------------------------------------------------------------
!> generate next rho_m by DIIS method
  subroutine DIIS()
    implicit none
    integer             :: jj, kk, ll, mm   ! loop variable DIIS
    if (iter <= nodiis) then
      !--------------<damping>-----------------
      if (abs(molE - molE_pre) >= cutdamp) then
        damp_coe = damp
        rho_m = (1.0_dp-damp_coe) * rho_m + damp_coe * rho_pre
        write(60,'(A,F6.3)') '  -- damped ', damp_coe
      else if (abs(molE - molE_pre) >= 0.01*cutdamp) then
        damp_coe = max(0.5_dp*damp*log10(abs(molE-molE_pre)) + &
        damp*(1.0_dp-0.5_dp*log10(cutdamp)), damp_coe-0.05_dp)
        rho_m = (1.0_dp-damp_coe) * rho_m + damp_coe * rho_pre
        write(60,'(A,F6.3)') '  -- damped ', damp_coe
      else
        damp_coe = 0.0_dp
        write(60,'(A)') '  -- undamped '
      end if
      !--------------<damping>-----------------
      if (iter <= nodiis - subsp) then
        write(60,'(A)') '  -- no DIIS acceleration'
      else
        ! update Rsd
        do jj = 1, 2*sbdm
          do kk = 1, 2*sbdm
            Rsd(iter-(nodiis-subsp), jj, kk) = rho_m(jj, kk) - rho_pre(jj, kk)
          end do
        end do
        ! update rho_history
        do jj = 1, 2*sbdm
          do kk = 1, 2*sbdm
            rho_history(iter-(nodiis-subsp), jj, kk) = rho_m(jj, kk)
          end do
        end do
        write(60,'(A,I2,A,I2)') &
        '  -- DIIS subspace filling ',iter-(nodiis-subsp),'/',subsp
      end if
    else
      !--------------<DIIS damping>-----------------
      if (abs(molE - molE_pre) >= cutdiis) then
        damp_coe = diisdamp
        write(60,'(A,F6.3)') '  -- DIIS damped ', damp_coe
      else if (abs(molE - molE_pre) >= cutdiis/100.0) then
        damp_coe = max(0.5_dp*diisdamp*log10(abs(molE-molE_pre)) + &
        diisdamp*(1.0_dp-0.5_dp*log10(cutdiis)), damp_coe-0.05_dp)
        write(60,'(A,F6.3)') '  -- DIIS damped ', damp_coe
      else
        damp_coe = 0.0_dp
        write(60,'(A)') '  -- DIIS undamped'
      end if
      !--------------<DIIS damping>-----------------
      ! update Rsd
      do jj = 2, subsp
        do kk = 1, 2*sbdm
          do ll = 1, 2*sbdm
            Rsd(jj-1, kk, ll) = Rsd(jj, kk, ll)
          end do
        end do
      end do
      do jj = 1, 2*sbdm
        do kk = 1, 2*sbdm
          Rsd(subsp, jj, kk) = rho_m(jj, kk) - rho_pre(jj, kk)
        end do
      end do
      ! update rho_history
      do jj = 2, subsp
        do kk = 1, 2*sbdm
          do ll = 1, 2*sbdm
            rho_history(jj-1, kk, ll) = rho_history(jj, kk, ll)
          end do
        end do
      end do
      do jj = 1, 2*sbdm
        do kk = 1, 2*sbdm
          rho_history(subsp, jj, kk) = rho_m(jj, kk)
        end do
      end do
      ! construct DIISmat
      DIISmat = c0
      do jj = 1, subsp
        do kk = 1, subsp
          do ll = 1, 2*sbdm
            do mm = 1, 2*sbdm
              DIISmat(jj, kk) = DIISmat(jj, kk) + &
              conjg(Rsd(jj, ll, mm))*Rsd(kk, ll, mm)
            end do
          end do
        end do
      end do
      do jj = 1, subsp
        DIISmat(subsp+1, jj) = c1
        DIISmat(jj, subsp+1) = c1
      end do
      ! solve residual equation
      ! dgesv and dspsv will cause Integral_V_2e conflict for unknown reason
      ! since DIISmat (and its inverse) is real symmetric, plus the column 
      ! vector is simple, use the inverse of DIISmat to solve directly
      call inverse(DIISmat, subsp+1)
      ! generate new rho_m
      rho_m = c0
      do jj = 1, subsp
        do kk = 1, 2*sbdm
          do ll = 1, 2*sbdm
            rho_m(kk,ll) = rho_m(kk,ll) + DIISmat(jj,subsp+1) * &
            (rho_history(jj,kk,ll)+damp_coe*Rsd(jj,kk,ll))
          end do
        end do
      end do
      write(60,'(A,E10.3,A,E10.3,A)') '  -- predicted residual (', &
      -real(DIISmat(subsp+1,subsp+1)), ',', -aimag(DIISmat(subsp+1,subsp+1)),')'
      do jj = 1, subsp
        write(60,'(A,I2,A,E10.3,A,E10.3,A)') '  -- subsp coeff', jj,'      (', &
        real(DIISmat(jj,subsp+1)), ',', aimag(DIISmat(jj,subsp+1)), ')'
      end do
    end if
  end subroutine DIIS

!------------------------------------------------------------
!> print wave-function information during SCF process
!!
!! will calculate molE and molE_pre
  subroutine SCFprint()
    implicit none
    integer             :: jj   ! loop variable for SCFprint
    ! frontier orbital energy
    write(60,'(A)') '  frontier orbital energy (A.U.)'
    call Calc_S2HForb(occindex(electron_count))
    write(60,'(A,I3,F12.6,A,F6.3)') &
    '  -- HOMO ', occindex(electron_count), orbE(occindex(electron_count)), &
    ' <Sz> = ',Szorb
    call Calc_S2HForb(occindex(electron_count)+1)
    write(60,'(A,I3,F12.6,A,F6.3)') &
    '  -- LUMO ',occindex(electron_count)+1,orbE(occindex(electron_count)+1),&
    ' <Sz> = ',Szorb
    write(60,'(A,F12.6)') '  -- gap ',&
    orbE(occindex(electron_count)+1)-orbE(occindex(electron_count))
    ! energy components calculation
    write(60,'(A)') '  calculate energy components (A.U.)'
    Ecore = 0.0_dp
    T = 0.0_dp
    V = 0.0_dp
    EpVp = 0.0_dp
    ESR = 0.0_dp
    call matmul('C', 'N', oper3, Fock1, oper6)
    call matmul('N', 'N', oper6, oper3, oper4)
    do jj = 1, electron_count
      Ecore = Ecore + real(oper4(occindex(jj),occindex(jj)))
    end do
    ! kinetic energy
    call matmul('C', 'N', oper3, exi_T_j, oper6)
    call matmul('N', 'N', oper6, oper3, oper4)
    do jj = 1, electron_count
      T = T + real(oper4(occindex(jj),occindex(jj)))
    end do
    ! electron-nuclear attraction energy
    call matmul('C', 'N', oper3, AVA, oper6)
    call matmul('N', 'N', oper6, oper3, oper4)
    do jj = 1, electron_count
      V = V + real(oper4(occindex(jj),occindex(jj)))
    end do
    if (pVp1e) then
      call matmul('C', 'N', oper3, expVp, oper6)
      call matmul('N', 'N', oper6, oper3, oper4)
      do jj = 1, electron_count
        EpVp = EpVp + real(oper4(occindex(jj),occindex(jj)))
      end do
      if (pppVp) then
        call matmul('C', 'N', oper3, exSR, oper6)
        call matmul('N', 'N', oper6, oper3, oper4)
        do jj = 1, electron_count
          ESR = ESR + real(oper4(occindex(jj),occindex(jj)))
        end do
      end if
    end if
    write(60,'(A,F12.6)') '  -- One-electron (core) energy           ', Ecore
    write(60,'(A,F12.6)') '  -- -- Kinetic                           ', T
    write(60,'(A,F12.6)') '  -- -- Electron-nuclear attraction       ', V
    write(60,'(A,F12.6)') '  -- -- pVp-related                       ', EpVp
    write(60,'(A,F12.6)') '  -- -- pppVp-related                     ', ESR
    E2e = 0.0_dp
    HFCol = 0.0_dp
    HFexc = 0.0_dp
    EpVpcol = 0.0_dp
    EpVpexc = 0.0_dp
    ! HF Coulomb energy
    call matmul('C', 'N', oper3, mHFcol, oper6)
    call matmul('N', 'N', oper6, oper3, oper4)
    do jj = 1, electron_count
      HFCol = HFCol + real(oper4(occindex(jj),occindex(jj)))
    end do
    ! HF exchange energy
    call matmul('C', 'N', oper3, mHFexc, oper6)
    call matmul('N', 'N', oper6, oper3, oper4)
    do jj = 1, electron_count
      HFexc = HFexc + real(oper4(occindex(jj),occindex(jj)))
    end do
    if (pVp2e) then
      call matmul('C', 'N', oper3, mpVpcol_11+mpVpcol_22, oper6)
      call matmul('N', 'N', oper6, oper3, oper4)
      do jj = 1, electron_count
        EpVpcol = EpVpcol + real(oper4(occindex(jj),occindex(jj)))
      end do
      call matmul('C', 'N', oper3, mpVpexc_11+mpVpexc_22, oper6)
      call matmul('N', 'N', oper6, oper3, oper4)
      do jj = 1, electron_count
        EpVpexc = EpVpexc + real(oper4(occindex(jj),occindex(jj)))
      end do
    end if
    E2e = 0.5_dp*(HFCol+x_HF*HFexc+EpVpcol+EpVpexc)
    write(60,'(A,F12.6)') &
    '  -- Two-electron energy                  ', E2e
    write(60,'(A,F12.6)') &
    '  -- -- Coulomb                           ', 0.5_dp*HFCol
    write(60,'(A,F12.6)') &
    '  -- -- HF Exchange                       ', 0.5_dp*x_HF*HFexc
    write(60,'(A,F12.6)') &
    '  -- -- pVp-related Coulomb               ', 0.5_dp*EpVpcol
    write(60,'(A,F12.6)') &
    '  -- -- pVp-related Exchange              ', 0.5_dp*EpVpexc
    write(60,'(A,F12.6)') &
    '  -- KS Exchange                          ', (1-x_HF)*KSexc
    write(60,'(A,F12.6)') &
    '  -- KS Correlation                       ', KScor

    ! electronic energy
    if (iter /= 1) molE_pre = molE
    molE = nucE + Ecore + E2e + (1.0_dp-x_HF)*KSexc + KScor

    ! (non-relativistic) Virial ratio
    Virial = -(nucE+Ecore-T-EpVp-ESR+0.5_dp*(HFCol+x_HF*HFexc)+&
    (1.0_dp-x_HF)*KSexc+KScor) / T
    write(60,'(A,F12.6)') &
    '  -- -<V>/<T>                             ', Virial
  end subroutine SCFprint

!------------------------------------------------------------
!> print wave-function information after SCF process
  subroutine Outputprint()
    implicit none
    integer :: iatom, ishell
    integer :: ii, jj              ! loop variable Outputprint
    close(80)
    write(60,*)
    write(60,*)
    write(60,'(A)') &
    '  ============================================================='
    write(60,'(A)') '                           MOL INFO'
    write(60,'(A)') &
    '  ============================================================='
    write(60,'(A,F12.6)') &
    '  Total electronic energy / Eh                  ...',molE
    write(60,'(A,F12.6)') &
    '  Nuclear repulsive energy / Eh                 ...',nucE
    write(60,'(A,F12.6)') &
    '  One-electron energy / Eh                      ...',Ecore
    write(60,'(A,F12.6)') &
    '  -- Kinetic energy / Eh                        ...',T
    write(60,'(A,F12.6)') &
    '  -- Electron-nuclear attraction energy / Eh    ...',V
    write(60,'(A,F12.6)') &
    '  -- pVp-related energy / Eh                    ...',EpVp
    write(60,'(A,F12.6)') &
    '  -- pppVp-related energy / Eh                  ...',ESR
    write(60,'(A,F12.6)') &
    '  Two-electron energy / Eh                      ...',E2e
    write(60,'(A,F12.6)') &
    '  -- HF Coulomb energy / Eh                     ...',0.5_dp*HFCol
    write(60,'(A,F12.6)') &
    '  -- HF Exchange energy / Eh                    ...',0.5_dp*x_HF*HFexc
    write(60,'(A,F12.6)') &
    '  -- pVp-related Coulomb energy / Eh            ...',0.5_dp*EpVpcol
    write(60,'(A,F12.6)') &
    '  -- pVp-related Exchange energy / Eh           ...',0.5_dp*EpVpexc
    write(60,'(A,F12.6)') &
    '  KS Exchange energy / Eh                       ...',(1.0_dp-x_HF)*KSexc
    write(60,'(A,F12.6)') &
    '  KS Correlation energy / Eh                    ...',KScor
    if (d4) write(60,'(A,F12.6)') &
    '  Dispersion energy (DFT-D4) / Eh               ...',emd4
    write(60,'(A,F12.6)') &
    '  Virial ratio                                  ...',Virial
    if (pVp1e) then
      write(60,'(A)') '  -- Note: relativistic calculation causes the'
      write(60,'(A)') '  -- Virial ratio to deviate (usually below) 2.0'
    end if
    write(60,'(A,F12.6)') &
    '  Total Alpha electron                          ...',totalpha
    write(60,'(A,F12.6)') &
    '  Total Beta electron                           ...',totbeta
    write(60,'(A,F12.6)') &
    '  <Sz*(Sz+1)> / hbar**2                         ...',&
    ((totalpha-totbeta)/2.0)*((totalpha-totbeta)/2.0+1.0_dp)
    write(60,'(A,F12.6)') &
    '  <S**2> / hbar**2                              ...',S__2
    if (fx_id /= 0) then
      write(60,'(A)') '  -- Note: there is little theoretical justification'
      write(60,'(A)') '  -- to calculate <S**2> in a DFT calculation.'
    end if
    scf_kappa = Krammers()
    write(60,'(A)') '  Time Reversal Symmetry (TRS) deviation parameter'
    write(60,'(A,E12.5)') '  -- kappa     =', scf_kappa
    write(60,'(A,E12.5)') '  -- ref kappa =', dsqrt(real(Nalpha-Nbeta,dp))
    write(60,'(A,E12.5)') '  -- SOC kappa =', &
    scf_kappa - dsqrt(real(Nalpha-Nbeta,dp))
    write(60,'(A)') &
    '  ============================================================='
    write(60,*)
    write(60,*)
    write(60,*)
    call RGI()
    write(60,*)
    write(60,*)
    write(60,*)
    write(60,'(A)') &
    '  ============================================================='
    write(60,'(A)') &
    '                       CANONICAL ORB INFO'
    write(60,'(A)') &
    '  ============================================================='
    do ii = 1, electron_count
      if (ii < electron_count) then
        call Calc_S2HForb(ii)
        write(60,'(A,I3.3,A,F12.6)') &
        '  HOMO-', electron_count-ii, &
        '   energy / Eh                       ... ', orbE(ii)
        write(60,'(A3,I3.3,A,F12.6)') &
        '  #',ii,'       <S**2> / hbar**2                  ... ', S__2orb
        write(60,'(A,F12.6)') &
        '             <Sz> / hbar                       ... ', Szorb
        write(60,'(A)') &
        '                 orb-in-atom  atom-in-mol      RE        IM'
        iatom = 1
        ishell = 1
        do jj = 1, sbdm
          if (sbdata(jj)%atom /= iatom) then
            ishell = 1
            iatom = iatom + 1
          end if
          if (abs(AO2MO(jj,ii)) >= prtlev) then
            write(60,'(A,I2.2,A,I1,A,A2,A,I2.2,A,F9.6,A1,F9.6,A1)') &
            '             --   A #',ishell,' L=',sbdata(jj)%L-1,'     ',&
            element_list(mol(sbdata(jj)%atom)%atom_number),' #',sbdata(jj)%atom&
            , '    (', real(AO2MO(jj,ii)), ',', aimag(AO2MO(jj,ii)), ')'
          end if
          if (abs(AO2MO(sbdm+jj,ii)) >= prtlev) then
            write(60,'(A,I2.2,A,I1,A,A2,A,I2.2,A,F9.6,A1,F9.6,A1)') &
            '             --   B #',ishell,' L=',sbdata(jj)%L-1,'     ',&
            element_list(mol(sbdata(jj)%atom)%atom_number),' #',sbdata(jj)%atom&
            , '    (', real(AO2MO(sbdm+jj,ii)),',',aimag(AO2MO(sbdm+jj,ii)), ')'
          end if
          ishell = ishell + 1
        end do
      else
        call Calc_S2HForb(ii)
        write(60,'(A,A,F12.6)') &
        '  HOMO    ', '   energy / Eh                       ... ', orbE(ii)
        write(60,'(A3,I3.3,A,F12.6)') &
        '  #',ii,'       <S**2> / hbar**2                  ... ', S__2orb
        write(60,'(A,F12.6)') &
        '             <Sz> / hbar                       ... ', Szorb
        write(60,'(A)') &
        '                 orb-in-atom  atom-in-mol      RE        IM'
        iatom = 1
        ishell = 1
        do jj = 1, sbdm
          if (sbdata(jj)%atom /= iatom) then
            ishell = 1
            iatom = iatom + 1
          end if
          if (abs(AO2MO(jj,ii)) >= prtlev) then
            write(60,'(A,I2.2,A,I1,A,A2,A,I2.2,A,F9.6,A1,F9.6,A1)') &
            '             --   A #',ishell,' L=',sbdata(jj)%L-1,'     ',&
            element_list(mol(sbdata(jj)%atom)%atom_number),' #',sbdata(jj)%atom&
            , '    (', real(AO2MO(jj,ii)), ',', aimag(AO2MO(jj,ii)), ')'
          end if
          if (abs(AO2MO(sbdm+jj,ii)) >= prtlev) then
            write(60,'(A,I2.2,A,I1,A,A2,A,I2.2,A,F9.6,A1,F9.6,A1)') &
            '             --   B #',ishell,' L=',sbdata(jj)%L-1,'     ',&
            element_list(mol(sbdata(jj)%atom)%atom_number),' #',sbdata(jj)%atom&
            , '    (', real(AO2MO(sbdm+jj,ii)),',',aimag(AO2MO(sbdm+jj,ii)), ')'
          end if
          ishell = ishell + 1
        end do
      end if
    end do
    do ii = electron_count+1, electron_count+5
      if (ii == electron_count + 1) then
        call Calc_S2HForb(ii)
        write(60,'(A,A,F12.6)') &
        '  LUMO    ', '   energy / Eh                       ... ', orbE(ii)
        write(60,'(A3,I3.3,A,F12.6)') &
        '  #',ii,'       <S**2> / hbar**2                  ... ', S__2orb
        write(60,'(A,F12.6)') &
        '             <Sz> / hbar                       ... ', Szorb
        write(60,'(A)') &
        '                 orb-in-atom  atom-in-mol      RE        IM'
        iatom = 1
        ishell = 1
        do jj = 1, sbdm
          if (sbdata(jj)%atom /= iatom) then
            ishell = 1
            iatom = iatom + 1
          end if
          if (abs(AO2MO(jj,ii)) >= prtlev) then
            write(60,'(A,I2.2,A,I1,A,A2,A,I2.2,A,F9.6,A1,F9.6,A1)') &
            '             --   A #',ishell,' L=',sbdata(jj)%L-1,'     ',&
            element_list(mol(sbdata(jj)%atom)%atom_number),' #',sbdata(jj)%atom&
            , '    (', real(AO2MO(jj,ii)), ',', aimag(AO2MO(jj,ii)), ')'
          end if
          if (abs(AO2MO(sbdm+jj,ii)) >= prtlev) then
            write(60,'(A,I2.2,A,I1,A,A2,A,I2.2,A,F9.6,A1,F9.6,A1)') &
            '             --   B #',ishell,' L=',sbdata(jj)%L-1,'     ',&
            element_list(mol(sbdata(jj)%atom)%atom_number),' #',sbdata(jj)%atom&
            , '    (', real(AO2MO(sbdm+jj,ii)),',',aimag(AO2MO(sbdm+jj,ii)), ')'
          end if
          ishell = ishell + 1
        end do
      else
        call Calc_S2HForb(ii)
        write(60,'(A,I3.3,A,F12.6)') &
        '  LUMO+', ii-electron_count-1, &
        '   energy / Eh                       ... ', orbE(ii)
        write(60,'(A3,I3.3,A,F12.6)') &
        '  #',ii,'       <S**2> / hbar**2                  ... ', S__2orb
        write(60,'(A,F12.6)') &
        '             <Sz> / hbar                       ... ', Szorb
        write(60,'(A)') &
        '                 orb-in-atom  atom-in-mol      RE        IM'
        iatom = 1
        ishell = 1
        do jj = 1, sbdm
          if (sbdata(jj)%atom /= iatom) then
            ishell = 1
            iatom = iatom + 1
          end if
          if (abs(AO2MO(jj,ii)) >= prtlev) then
            write(60,'(A,I2.2,A,I1,A,A2,A,I2.2,A,F9.6,A1,F9.6,A1)') &
            '             --   A #',ishell,' L=',sbdata(jj)%L-1,'     ',&
            element_list(mol(sbdata(jj)%atom)%atom_number),' #',sbdata(jj)%atom&
            , '    (', real(AO2MO(jj,ii)), ',', aimag(AO2MO(jj,ii)), ')'
          end if
          if (abs(AO2MO(sbdm+jj,ii)) >= prtlev) then
            write(60,'(A,I2.2,A,I1,A,A2,A,I2.2,A,F9.6,A1,F9.6,A1)') &
            '             --   B #',ishell,' L=',sbdata(jj)%L-1,'     ',&
            element_list(mol(sbdata(jj)%atom)%atom_number),' #',sbdata(jj)%atom&
            , '    (', real(AO2MO(sbdm+jj,ii)),',',aimag(AO2MO(sbdm+jj,ii)), ')'
          end if
          ishell = ishell + 1
        end do
      end if
    end do
    write(60,'(A)') &
    '  ============================================================='
    write(60,*)
    if (molden) then
      write(60,'(A)') '  dumping AO2MO to '//trim(address_job)//'.molden.d'
      if (.not. pVp1e) then
        write(60,'(A)') &
        '  -- Note: for scalar MOs, only realpart will be generated.'
      end if
      call Dump_MOLDEN()
      write(60,'(A)') '  complete!'
    end if
    write(60,'(A)') 'exit module SCF'
  end subroutine Outputprint

!------------------------------------------------------------
!> initialization of global variables
  subroutine Globinit(kill)
    implicit none
    logical,intent(in) :: kill   ! kill the process after SCF
    deallocate(exi_T_j, exi_V_j)
    deallocate(oper4, oper6)
    deallocate(Rsd, DIISmat)
    deallocate(rho_history, rho_pre)
    deallocate(iijj_V)
    if (pVp2e) then
      deallocate(iijj_pxVpx, iijj_pyVpy, iijj_pzVpz)
      deallocate(ijij_pxVpx, ijij_pyVpy, ijij_pzVpz)
    end if
    deallocate(Fock1, mHFexc, mHFcol)
    if (pVp2e) then
      deallocate(mpVpexc_11, mpVpcol_11, mpVpexc_22, mpVpcol_22)
    end if
    if (allocated(mKScor)) deallocate(mKScor)
    if (allocated(mKSexc)) deallocate(mKSexc)
    deallocate(i_j, i_p2_j, i_V_j, c2s, exc2s)
    deallocate(s2f, exs2f, f2s, exf2s, c2f, exc2f)
    deallocate(c2soper, exc2soper, s2coper, exs2coper, f2soper, exf2soper)
    deallocate(s2foper, exs2foper, c2foper, exc2foper)
    if (allocated(M_c2s)) deallocate(M_c2s, M_exc2s)
    if (kill) then
      if (fx_id /= -1) call Fockxc_end()
      deallocate(AO2MO, rho_m, Fock, orbE, oper3, occindex)
    end if
    ini_rho = .true.
    deallocate(AVA)
    if (pVp1e) then
      deallocate(expVp, AO2p2, evl_p2, Ap, ApRp, SRp, ARVRA, exAO2p2)
      if (pppVp) deallocate(exSR)
    end if
    deallocate(cbdata)
    deallocate(sbdata)
    if (allocated(M_cbdata)) deallocate(M_cbdata)
    if (kill) then
      call terminate('normal')
    else
      call terminate('keep')
    end if
  end subroutine Globinit
  
!------------------------------------------------------------
!> initial guess and generate density matix
  subroutine Assign_rho()
    implicit none
    integer              :: ii, jj, kk    ! loop variables for Assign_rho
    integer              :: Na, Nb
    integer              :: degenlow, degenhigh, load  ! degenerat region
    real(dp)             :: maxval
    integer              :: imax
    if (ini_rho) then
      if (.not. allocated(rho_m)) allocate(rho_m(2*sbdm,2*sbdm), source=c0)
      if (.not. allocated(occindex)) allocate(occindex(electron_count))
      ! load MO coefficient
      !-------------------------load from MOLDEN-------------------------
      if (guess_type == 'molden') then
        if (.not. allocated(AO2MO)) allocate(AO2MO(2*sbdm,2*fbdm))
        call load_geombasis_MOLDEN(.true.)
        if (M_atom_count /= atom_count) call terminate(&
        'n_atoms in MOLDEN is not consistent with n_atoms in .xyz')
        if (isorca) then
          ! The GTF coefficients in the ORCA MOLDEN file are
          ! already normalized at the GTO level
          do ii = 1, M_cbdm
            M_cbdata(ii)%Ncoe = M_cbdata(ii)%coe
          end do
        else
          call Calc_Ncoe(M_cbdata, M_cbdm)
        end if
        write(60,'(A)') '  -- basis and geometry in MOLDEN were loaded'
        write(60,'(A,I4,A,I4)') '  -- M_cbdm / M_sbdm: ',M_cbdm, ' /', M_sbdm
        allocate(M_AO2MO_a(M_sbdm,M_sbdm))
        allocate(M_AO2MO_b(M_sbdm,M_sbdm))
        allocate(T_AO2MO_a(sbdm,M_sbdm))
        allocate(T_AO2MO_b(sbdm,M_sbdm))
        call load_MO_MOLDEN(M_AO2MO_a, M_AO2MO_b)
        write(60,'(A)') '  -- MO coeffs in MOLDEN were loaded'
        call M_basis_proj(M_AO2MO_a, T_AO2MO_a)
        call M_basis_proj(M_AO2MO_b, T_AO2MO_b)
        write(60,'(A)') '  -- MO coeffs were projected to job basis'

        do ii = 1,sbdm
          do jj = 1,sbdm
            do kk = 1,Nalpha
              rho_m(ii,jj) = rho_m(ii,jj) + &
              T_AO2MO_a(ii,kk)*T_AO2MO_a(jj,kk)
            end do
            do kk = 1,Nbeta
              rho_m(sbdm+ii,sbdm+jj) = rho_m(sbdm+ii,sbdm+jj) + &
              T_AO2MO_b(ii,kk)*T_AO2MO_b(jj,kk)
            end do
          end do
        end do
        deallocate(M_AO2MO_a, M_AO2MO_b, T_AO2MO_a, T_AO2MO_b)
      !-------------------------load from .ao2mo-------------------------
      else if (guess_type == 'ao2mo') then
        if (allocated(AO2MO)) deallocate(AO2MO)
        call load_matrix('ao2mo', AO2MO, ii, jj)
        if (ii /= 2*sbdm .or. jj /= 2*fbdm) call terminate(&
        'basis dimension in .ao2mo file mismatch with current job')
        do ii = 1, 2*sbdm
          do jj = 1, 2*sbdm
            do kk = 1, electron_count
              rho_m(ii,jj) = rho_m(ii,jj) + AO2MO(ii,kk)*conjg(AO2MO(jj,kk))
            end do
          end do
        end do
      end if
      rho_pre = rho_m
    ! keep spin multiplicity if frontier MOs are degenerate
    ! should keep oper3 and AO2MO unchanged, only change the occindex
    else if (.not.ini_rho .and. (cspin=='n' .or. cspin=='d')) then
      forall (ii=1:electron_count) occindex(ii) = ii
      if (allocated(rho_m)) deallocate(rho_m)
      allocate(rho_m(2*sbdm,2*sbdm), source=c0)
      do ii = 1, 2*sbdm
        do jj = 1, 2*sbdm
          do kk = 1, electron_count
            rho_m(ii,jj) = rho_m(ii,jj) + AO2MO(ii,kk)*conjg(AO2MO(jj,kk))
          end do
        end do
      end do
      call Calc_S2HF()
      if (cspin=='d' &
      .and. abs((totalpha-totbeta) - real(Nalpha-Nbeta,dp)) > 1.0 &
      .and. abs((totalpha-real(Nalpha,dp))-(real(Nbeta,dp)-totbeta)) < 0.1) then
        write(60,'(A)') &
        '  -- order of degenerate frontier alpha/beta orbitals changed'
        do ii = electron_count-1, 1, -1
          if (abs(orbE(ii)-orbE(electron_count)) /&
          abs(orbE(electron_count)) > 0.04) then
            degenlow = ii + 1
            exit
          end if
        end do
        do ii = electron_count+1, 6*electron_count
          if (abs(orbE(ii)-orbE(electron_count)) /&
          abs(orbE(electron_count)) > 0.04) then
            degenhigh = ii - 1
            exit
          end if
        end do
        write(60,'(A,I3,A,I3)') &
        '  -- -- degenerate space: ',degenlow,' - ',degenhigh
        Na = 0
        Nb = 0
        do kk = 1, degenlow-1
          call Calc_S2HForb(kk)
          if (Szorb > 0.0) then
            Na = Na + 1
          else
            Nb = Nb + 1
          end if
        end do
        write(60,'(A,I3,A3,I3)') &
        '  -- -- alpha and beta in nondegenerate region: ', Na, ' / ', Nb
        load = degenlow
        do kk = degenlow, degenhigh
          call Calc_S2HForb(kk)
          if (Szorb > 0.0 .and. Na < Nalpha) then
            write(60,'(A,I3,A,I3)') '  -- -- load alpha ',kk,' on ',load
            occindex(load) = kk
            load = load + 1
            Na = Na + 1
          else if(Szorb<0.0 .and. Nb<Nbeta) then
            write(60,'(A,I3,A,I3)') '  -- -- load beta ',kk,' on ',load
            occindex(load) = kk
            load = load + 1
            Nb = Nb + 1
          end if
          if (load == electron_count + 1) exit
        end do
        if (load < electron_count + 1) call terminate(&
        "cspin error, can't keep spin in degenerate space")
        rho_m = c0
        do ii = 1, 2*sbdm
          do jj = 1, 2*sbdm
            do kk = 1, electron_count
              rho_m(ii,jj) = rho_m(ii,jj) + &
              AO2MO(ii,occindex(kk)) * conjg(AO2MO(jj,occindex(kk)))
            end do
          end do
        end do
      end if
    ! force keep spin multiplicity
    ! Can't use constrained-DFT because the number of alpha and beta electrons
    ! is not strictly integer.
    else if(.not.ini_rho .and. cspin=='f') then
      Na = 0
      Nb = 0
      load = 1
      do kk = 1, 2*fbdm
        call Calc_S2HForb(kk)
        if (Szorb > 0.0 .and. Na < Nalpha) then
          occindex(load) = kk
          load = load + 1
          Na = Na + 1
        else if (Szorb < 0.0 .and. Nb < Nbeta) then
          occindex(load) = kk
          load = load + 1
          Nb = Nb + 1
        end if
        if (Na == Nalpha .and. Nb == Nbeta) exit
      end do
      if (allocated(rho_m)) deallocate(rho_m)
      allocate(rho_m(2*sbdm,2*sbdm), source=c0)
      do ii = 1, 2*sbdm
        do jj = 1, 2*sbdm
          do kk = 1, electron_count
            rho_m(ii,jj) = rho_m(ii,jj) + &
            AO2MO(ii,occindex(kk)) * conjg(AO2MO(jj,occindex(kk)))
          end do
        end do
      end do
    end if
    if (.not. ini_rho) then
      if (iter /= 1) then
        maxDP = (real(rho_m(1,1))-real(rho_pre(1,1)))**2 + &
        (aimag(rho_m(1,1))-aimag(rho_pre(1,1)))**2
        RMSDP = 0.0_dp
        do ii = 1, 2*sbdm
          do jj = 1, 2*sbdm
            if ((real(rho_m(ii,jj))-&
            real(rho_pre(ii,jj)))**2 + (aimag(rho_m(ii,jj))&
            -aimag(rho_pre(ii,jj)))**2 > maxDP) then
              maxDP = (real(rho_m(ii,jj))-real(rho_pre(ii,jj)))**2 + &
              (aimag(rho_m(ii,jj))-aimag(rho_pre(ii,jj)))**2
            end if
            RMSDP = RMSDP + &
            (real(rho_m(ii,jj))-real(rho_pre(ii,jj)))**2 + &
            (aimag(rho_m(ii,jj))-aimag(rho_pre(ii,jj)))**2
          end do
        end do
        maxDP = dsqrt(maxDP)
        RMSDP = RMSDP / (4*sbdm*sbdm)
        RMSDP = dsqrt(RMSDP)
        write(60,'(A,E10.3)') '  -- maxDP                  ', maxDP
        write(60,'(A,E10.3)') '  -- RMSDP                  ', RMSDP
      end if
      ! calc <S**2>
      call Calc_S2HF()
      write(60,'(A,F9.5)') '  -- <S**2>                ',S__2
      write(60,'(A,F9.5)') '  -- total alpha electron  ',totalpha
      write(60,'(A,F9.5)') '  -- total beta electron   ',totbeta
    else
      ini_rho = .false.
    end if
  end subroutine Assign_rho
  
!------------------------------------------------------------
!> assign Coulomb matrix for Assign_Fock_2e
!!
!! ref: page 261 of "Quantum Chemistry: Basic Principles and
!!
!! ab-initio Calculations, Volume 2 | 2nd Edition"
  pure subroutine Assign_Coulomb(Col, rho, int, ui, uj, uk, ul)
    implicit none
    complex(dp), intent(inout) :: Col(2*cbdm, 2*cbdm) ! Coulomb matrix
    complex(dp), intent(in)    :: rho(2*cbdm, 2*cbdm) ! (reduced) density matrix
    real(dp),    intent(in)    :: int                 ! scalar integral (ij|kl)
    integer,     intent(in)    :: ui, uj, uk, ul
    ! avoid duplicate assignment of Col
    complex(dp)                :: ospin               ! oalpha = obeta = ospin
    ! matrix element change of sigma operators
    integer                    :: Col_assigned(2,8)
    integer                    :: assigned
    integer                    :: iassign
    logical                    :: carry

    Col_assigned = 0
    carry = .true.
    assigned = 1
    Col_assigned(1,1) = ui
    Col_assigned(2,1) = uj
    ospin = int*rho(uk,ul) + int*rho(cbdm+uk,cbdm+ul)
    if (uk /= ul) ospin = ospin + int*rho(ul,uk) + int*rho(cbdm+ul,cbdm+uk)
    Col(ui,uj) = Col(ui,uj) + ospin
    Col(cbdm+ui,cbdm+uj) = Col(cbdm+ui,cbdm+uj) + ospin
    do iassign = 1, assigned
      if (Col_assigned(1,iassign)==uj .and. Col_assigned(2,iassign)==ui) then
        carry = .false.
        exit
      end if
    end do
    if (carry) then
      assigned = assigned + 1
      Col_assigned(1,assigned) = uj
      Col_assigned(2,assigned) = ui
      Col(uj,ui) = Col(uj,ui) + ospin
      Col(cbdm+uj,cbdm+ui) = Col(cbdm+uj,cbdm+ui) + ospin
    end if
    carry = .true.
    do iassign = 1, assigned
      if (Col_assigned(1,iassign)==uk .and. Col_assigned(2,iassign)==ul) then
        carry = .false.
        exit
      end if
    end do
    if (carry) then
      assigned = assigned + 1
      Col_assigned(1,assigned) = uk
      Col_assigned(2,assigned) = ul
      ospin = int*rho(ui,uj) + int*rho(cbdm+ui,cbdm+uj)
      if (ui /= uj) ospin = ospin + int*rho(uj,ui) + int*rho(cbdm+uj,cbdm+ui)
      Col(uk,ul) = Col(uk,ul) + ospin
      Col(cbdm+uk,cbdm+ul) = Col(cbdm+uk,cbdm+ul) + ospin
    end if
    carry = .true.
    do iassign = 1, assigned
      if (Col_assigned(1,iassign)==ul .and. Col_assigned(2,iassign)==uk) then
        carry = .false.
        exit
      end if
    end do
    if (carry) then
      assigned = assigned + 1
      Col_assigned(1,assigned) = ul
      Col_assigned(2,assigned) = uk
      Col(ul,uk) = Col(ul,uk) + ospin
      Col(cbdm+ul,cbdm+uk) = Col(cbdm+ul,cbdm+uk) + ospin
    end if
  end subroutine Assign_Coulomb

!------------------------------------------------------------
!> assign Exchange matrix for Assign_Fock_2e
!!
!! ref: page 261 of "Quantum Chemistry: Basic Principles and
!!
!! ab-initio Calculations, Volume 2 | 2nd Edition"
  pure subroutine Assign_Exchange(Exc, rho, int, ui, uj, uk, ul)
    implicit none
    complex(dp), intent(inout) :: Exc(2*cbdm, 2*cbdm) ! Exchange matrix
    complex(dp), intent(in)    :: rho(2*cbdm, 2*cbdm) ! (reduced) density matrix
    real(dp),    intent(in)    :: int                 ! scalar integral (ij|kl)
    integer,     intent(in)    :: ui, uj, uk, ul
    ! avoid duplicate assignment of Exc
    integer                    :: Exc_assigned(2,8)
    integer                    :: assigned
    integer                    :: iassign
    logical                    :: carry
    Exc_assigned = 0
    carry = .true.
    assigned = 1
    Exc_assigned(1,1) = ui
    Exc_assigned(2,1) = uk
    Exc(ui,uk) = Exc(ui,uk) - int*rho(uj,ul)
    Exc(cbdm+ui,cbdm+uk) = Exc(cbdm+ui,cbdm+uk) - int*rho(cbdm+uj,cbdm+ul)
    if (ui == uk .and. ul /= uj) then
      Exc(ui,uk) = Exc(ui,uk) - int*rho(ul,uj)
      Exc(cbdm+ui,cbdm+uk) = Exc(cbdm+ui,cbdm+uk) - int*rho(cbdm+ul,cbdm+uj)
    end if
    do iassign = 1, assigned
      if (Exc_assigned(1,iassign)==uk .and. Exc_assigned(2,iassign)==ui) then
        carry = .false.
        exit
      end if
    end do
    if (carry) then
      assigned = assigned + 1
      Exc_assigned(1,assigned) = uk
      Exc_assigned(2,assigned) = ui
      Exc(uk,ui) = Exc(uk,ui) - int*rho(ul,uj)
      Exc(cbdm+uk,cbdm+ui) = Exc(cbdm+uk,cbdm+ui) - int*rho(cbdm+ul,cbdm+uj)
      if (uk == ui .and. ul /= uj) then
        Exc(uk,ui) = Exc(uk,ui) - int*rho(uj,ul)
        Exc(cbdm+uk,cbdm+ui) = Exc(cbdm+uk,cbdm+ui) - int*rho(cbdm+uj,cbdm+ul)
      end if
    end if
    carry = .true.
    do iassign = 1, assigned
      if (Exc_assigned(1,iassign)==ui .and. Exc_assigned(2,iassign)==ul) then
        carry = .false.
        exit
      end if
    end do
    if (carry) then
      assigned = assigned + 1
      Exc_assigned(1,assigned) = ui
      Exc_assigned(2,assigned) = ul
      Exc(ui,ul) = Exc(ui,ul) - int*rho(uj,uk)
      Exc(cbdm+ui,cbdm+ul) = Exc(cbdm+ui,cbdm+ul) - int*rho(cbdm+uj,cbdm+uk)
      if (ui == ul .and. uk /= uj) then
        Exc(ui,ul) = Exc(ui,ul) - int*rho(uk,uj)
        Exc(cbdm+ui,cbdm+ul) = Exc(cbdm+ui,cbdm+ul) - int*rho(cbdm+uk,cbdm+uj)
      end if
    end if
    carry = .true.
    do iassign = 1, assigned
      if (Exc_assigned(1,iassign)==ul .and. Exc_assigned(2,iassign)==ui) then
        carry = .false.
        exit
      end if
    end do
    if (carry) then
      assigned = assigned + 1
      Exc_assigned(1,assigned) = ul
      Exc_assigned(2,assigned) = ui
      Exc(ul,ui) = Exc(ul,ui) - int*rho(uk,uj)
      Exc(cbdm+ul,cbdm+ui) = Exc(cbdm+ul,cbdm+ui) - int*rho(cbdm+uk,cbdm+uj)
      if (ul == ui .and. uk /= uj) then
        Exc(ul,ui) = Exc(ul,ui) - int*rho(uj,uk)
        Exc(cbdm+ul,cbdm+ui) = Exc(cbdm+ul,cbdm+ui) - int*rho(cbdm+uj,cbdm+uk)
      end if
    end if
    carry = .true.
    do iassign = 1, assigned
      if (Exc_assigned(1,iassign)==uj .and. Exc_assigned(2,iassign)==uk) then
        carry = .false.
        exit
      end if
    end do
    if (carry) then
      assigned = assigned + 1
      Exc_assigned(1,assigned) = uj
      Exc_assigned(2,assigned) = uk
      Exc(uj,uk) = Exc(uj,uk) - int*rho(ui,ul)
      Exc(cbdm+uj,cbdm+uk) = Exc(cbdm+uj,cbdm+uk) - int*rho(cbdm+ui,cbdm+ul)
      if (uj == uk .and. ul /= ui) then
        Exc(uj,uk) = Exc(uj,uk) - int*rho(ul,ui)
        Exc(cbdm+uj,cbdm+uk) = Exc(cbdm+uj,cbdm+uk) - int*rho(cbdm+ul,cbdm+ui)
      end if
    end if
    carry = .true.
    do iassign = 1, assigned
      if (Exc_assigned(1,iassign)==uk .and. Exc_assigned(2,iassign)==uj) then
        carry = .false.
        exit
      end if
    end do
    if (carry) then
      assigned = assigned + 1
      Exc_assigned(1,assigned) = uk
      Exc_assigned(2,assigned) = uj
      Exc(uk,uj) = Exc(uk,uj) - int*rho(ul,ui)
      Exc(cbdm+uk,cbdm+uj) = Exc(cbdm+uk,cbdm+uj) - int*rho(cbdm+ul,cbdm+ui)
      if (uk == uj .and. ui /= ul) then
        Exc(uk,uj) = Exc(uk,uj) - int*rho(ui,ul)
        Exc(cbdm+uk,cbdm+uj) = Exc(cbdm+uk,cbdm+uj) - int*rho(cbdm+ui,cbdm+ul)
      end if
    end if
    carry = .true.
    do iassign = 1, assigned
      if (Exc_assigned(1,iassign)==uj .and. Exc_assigned(2,iassign)==ul) then
        carry = .false.
        exit
      end if
    end do
    if (carry) then
      assigned = assigned + 1
      Exc_assigned(1,assigned) = uj
      Exc_assigned(2,assigned) = ul
      Exc(uj,ul) = Exc(uj,ul) - int*rho(ui,uk)
      Exc(cbdm+uj,cbdm+ul) = Exc(cbdm+uj,cbdm+ul) - int*rho(cbdm+ui,cbdm+uk)
      if (uj == ul .and. uk /= ui) then
        Exc(uj,ul) = Exc(uj,ul) - int*rho(uk,ui)
        Exc(cbdm+uj,cbdm+ul) = Exc(cbdm+uj,cbdm+ul) - int*rho(cbdm+uk,cbdm+ui)
      end if
    end if
    carry = .true.
    do iassign = 1, assigned
      if (Exc_assigned(1,iassign)==ul .and. Exc_assigned(2,iassign)==uj) then
        carry = .false.
        exit
      end if
    end do
    if (carry) then
      assigned = assigned + 1
      Exc_assigned(1,assigned) = ul
      Exc_assigned(2,assigned) = uj
      Exc(ul,uj) = Exc(ul,uj) - int*rho(uk,ui)
      Exc(cbdm+ul,cbdm+uj) = Exc(cbdm+ul,cbdm+uj) - int*rho(cbdm+uk,cbdm+ui)
      if (ul == uj .and. ui /= uk) then
        Exc(ul,uj) = Exc(ul,uj) - int*rho(ui,uk)
        Exc(cbdm+ul,cbdm+uj) = Exc(cbdm+ul,cbdm+uj) - int*rho(cbdm+ui,cbdm+uk)
      end if
    end if
  end subroutine Assign_Exchange

!------------------------------------------------------------
!> assign Coulomb matrix with Sigma for Assign_Fock_2e
!!
!! sigma = 1:x, 2:y, 3:z
!!
!! (sigmax_i,sigmay_j|px1Vpy1|k,l) = (sigmax_i,sigmay_j|px1Vpy1|l,k)
!!
!!=(sigmay_j,sigmax_i|py1Vpx1|k,l) = (sigmay_j,sigmax_i|py1Vpx1|l,k)
!!
!!=(k,l|px2Vpy2|sigmax_i,sigmay_j) = (l,k|px2Vpy2|sigmax_i,sigmay_j)
!!
!!=(k,l|py2Vpx2|sigmay_j,sigmax_i) = (l,k|py2Vpx2|sigmay_j,sigmax_i)
!!
!! if sigmai = sigmaj, i and j should be half-traversed
!!
!! if sigmai /= sigmaj, i and j should be traversed
  pure subroutine Assign_Sigma_Col(&
  sigmai, sigmaj, Col1, Col2, rho1, rho2, int, ui, uj, uk, ul)
    implicit none
    integer,intent(in)         :: sigmai, sigmaj     ! sigma operator on i/j
    complex(dp), intent(inout) :: Col1(2*cbdm,2*cbdm)! electron1 Coulomb matrix
    complex(dp), intent(in)    :: rho1(2*cbdm,2*cbdm)! (reduced) density matrix
    complex(dp), intent(inout) :: Col2(2*cbdm,2*cbdm)! electron2 Coulomb matrix
    complex(dp), intent(in)    :: rho2(2*cbdm,2*cbdm)! (reduced) density matrix
    real(dp),    intent(in)    :: int                ! scalar integral (ij|kl)
    integer,     intent(in)    :: ui, uj, uk, ul
    complex(dp)                :: ospin              ! oalpha = obeta = ospin
    complex(dp)                :: ci2a,ci2b,cj2a,cj2b! coeff of sigma operators
    ! matrix element change of sigma operators
    integer                    :: di2a,di2b,dj2a,dj2b
    select case (sigmai)
    case (1)           ! oalpha -> beta, obeta -> alpha
      ci2a = c1
      ci2b = c1
      di2a = cbdm
      di2b = 0
    case (2)           ! oalpha -> i*beta, obeta -> -i*alpha
      ci2a = -ci
      ci2b = ci
      di2a = cbdm
      di2b = 0
    case (3)           ! oalpha -> alpha, obeta -> -beta
      ci2a = c1
      ci2b = -c1
      di2a = 0
      di2b = cbdm
    end select
    select case (sigmaj)
    case (1)           ! oalpha -> beta, obeta -> alpha
      cj2a = c1
      cj2b = c1
      dj2a = cbdm
      dj2b = 0
    case (2)           ! oalpha -> i*beta, obeta -> -i*alpha
      cj2a = -ci
      cj2b = ci
      dj2a = cbdm
      dj2b = 0
    case (3)           ! oalpha -> alpha, obeta -> -beta
      cj2a = c1
      cj2b = -c1
      dj2a = 0
      dj2b = cbdm
    end select
    !-------------------------------------------------------------------
    ! (sigmax_i,sigmay_j|px1Vpy1|k,l) and (sigmax_i,sigmay_j|px1Vpy1|l,k)
    ! conjg is needed because it operate on <i|
    ospin = int*(rho1(uk,ul) + rho1(cbdm+uk,cbdm+ul))
    if (uk /= ul) ospin = ospin + int*(rho1(ul,uk) + rho1(cbdm+ul,cbdm+uk))
    Col1(di2a+ui,dj2a+uj) = Col1(di2a+ui,dj2a+uj) + conjg(ci2a)*cj2a*ospin
    Col1(di2b+ui,dj2b+uj) = Col1(di2b+ui,dj2b+uj) + conjg(ci2b)*cj2b*ospin
    !-------------------------------------------------------------------
    ! (sigmay_j,sigmax_i|py1Vpx1|k,l) and (sigmay_j,sigmax_i|py1Vpx1|l,k)
    ! conjg is needed because it operate on <j|
    if (sigmai /= sigmaj .or. ui /= uj) then
      Col1(dj2a+uj,di2a+ui) = Col1(dj2a+uj,di2a+ui) + ci2a*conjg(cj2a)*ospin
      Col1(dj2b+uj,di2b+ui) = Col1(dj2b+uj,di2b+ui) + ci2b*conjg(cj2b)*ospin
    end if
    !-------------------------------------------------------------------
    ! (k,l|px2Vpy2|sigmax_i,sigmay_j) and (l,k|px2Vpy2|sigmax_i,sigmay_j)
    ospin = int*(conjg(ci2a)*cj2a*rho2(di2a+ui,dj2a+uj) + &
    conjg(ci2b)*cj2b*rho2(di2b+ui,dj2b+uj))
    Col2(uk,ul) = Col2(uk,ul) + ospin
    Col2(cbdm+uk,cbdm+ul) = Col2(cbdm+uk,cbdm+ul) + ospin
    if (uk /= ul) then
      Col2(ul,uk) = Col2(ul,uk) + ospin
      Col2(cbdm+ul,cbdm+uk) = Col2(cbdm+ul,cbdm+uk) + ospin
    end if
    !-------------------------------------------------------------------
    ! (k,l|py2Vpx2|sigmay_j,sigmax_i) and (l,k|py2Vpx2|sigmay_j,sigmax_i)
    if (sigmai /= sigmaj .or. ui /= uj) then
      ospin = int*(ci2a*conjg(cj2a)*rho2(dj2a+uj,di2a+ui) + &
      ci2b*conjg(cj2b)*rho2(dj2b+uj,di2b+ui))
      Col2(uk,ul) = Col2(uk,ul) + ospin
      Col2(cbdm+uk,cbdm+ul) = Col2(cbdm+uk,cbdm+ul) + ospin
      if (uk /= ul) then
        Col2(ul,uk) = Col2(ul,uk) + ospin
        Col2(cbdm+ul,cbdm+uk) = Col2(cbdm+ul,cbdm+uk) + ospin
      end if
    end if
  end subroutine Assign_Sigma_Col

!------------------------------------------------------------
!> assign Exchange matrix with Sigma for Assign_Fock_2e
!!
!! sigma = 1:x, 2:y, 3:z
!!
!! (sigmax_i,j|px1Vpy2|sigmay_k,l) = (j,sigmax_i|px1Vpy2|l,sigmay_k)
!!
!!=(sigmay_k,l|px2Vpy1|sigmax_i,j) = (l,sigmay_k|px2Vpy1|j,sigmax_i)
!!
!! if sigmai = sigmak, i and k should be half-traversed
!!
!! if sigmai /= sigmak, i and k should be traversed
  pure subroutine Assign_Sigma_Exc(&
  sigmai, sigmak, Exc1, Exc2, rho1, rho2, int, ui, uj, uk, ul)
    implicit none
    integer,intent(in)         :: sigmai, sigmak     ! sigma operator on i/k
    complex(dp), intent(inout) :: Exc1(2*cbdm,2*cbdm)! electron1 Coulomb matrix
    complex(dp), intent(in)    :: rho1(2*cbdm,2*cbdm)! (reduced) density matrix
    complex(dp), intent(inout) :: Exc2(2*cbdm,2*cbdm)! electron2 Coulomb matrix
    complex(dp), intent(in)    :: rho2(2*cbdm,2*cbdm)! (reduced) density matrix
    real(dp),    intent(in)    :: int                ! scalar integral (ij|kl)
    integer,     intent(in)    :: ui, uj, uk, ul
    ! coefficients to get alpha/beta component after sigma operation
    complex(dp)                :: ci2a, ci2b, ck2a, ck2b
    ! matrix element change to get alpha/beta component after sigma operation
    integer                    :: di2a, di2b, dk2a, dk2b
    select case (sigmai)
    case (1)                 ! oalpha -> beta, obeta -> alpha
      ci2a = c1
      ci2b = c1
      di2a = cbdm
      di2b = 0
    case (2)                 ! oalpha -> i*beta, obeta -> -i*alpha
      ci2a = -ci
      ci2b = ci
      di2a = cbdm
      di2b = 0
    case (3)                 ! oalpha -> alpha, obeta -> -beta
      ci2a = c1
      ci2b = -c1
      di2a = 0
      di2b = cbdm
    end select
    select case (sigmak)
    case (1)                 ! oalpha -> beta, obeta -> alpha
      ck2a = c1
      ck2b = c1
      dk2a = cbdm
      dk2b = 0
    case (2)                 ! oalpha -> i*beta, obeta -> -i*alpha
      ck2a = -ci
      ck2b = ci
      dk2a = cbdm
      dk2b = 0
    case (3)                 ! oalpha -> alpha, obeta -> -beta
      ck2a = c1
      ck2b = -c1
      dk2a = 0
      dk2b = cbdm
    end select
    !-------------------------------------------------------------------
    ! (sigmax_i,j|px1Vpy2|sigmay_k,l)
    Exc1(di2a+ui,dk2a+uk) = Exc1(di2a+ui,dk2a+uk) - &
    conjg(ci2a)*ck2a*int*rho1(uj,ul)
    Exc1(di2b+ui,dk2b+uk) = Exc1(di2b+ui,dk2b+uk) - &
    conjg(ci2b)*ck2b*int*rho1(cbdm+uj,cbdm+ul)
    !-------------------------------------------------------------------
    ! (sigmay_k,l|px2Vpy1|sigmax_i,j)
    if (sigmai /= sigmak .or. ui /= uk .or. uj /= ul) then
      Exc1(dk2a+uk,di2a+ui) = Exc1(dk2a+uk,di2a+ui) - &
      ci2a*conjg(ck2a)*int*rho1(ul,uj)
      Exc1(dk2b+uk,di2b+ui) = Exc1(dk2b+uk,di2b+ui) - &
      ci2b*conjg(ck2b)*int*rho1(cbdm+ul,cbdm+uj)
    end if
    !-------------------------------------------------------------------
    ! (j,sigmax_i|px1Vpy2|l,sigmay_k)
    Exc2(uj,ul) = Exc2(uj,ul) - &
    ci2a*conjg(ck2a)*int*rho2(di2a+ui,dk2a+uk)
    Exc2(cbdm+uj,cbdm+ul) = Exc2(cbdm+uj,cbdm+ul) - &
    ci2b*conjg(ck2b)*int*rho2(di2b+ui,dk2b+uk)
    !-------------------------------------------------------------------
    ! (l,sigmay_k|px2Vpy1|j,sigmax_i)
    if (sigmai /= sigmak .or. ui /= uk .or. uj /= ul) then
      Exc2(ul,uj) = Exc2(ul,uj) - &
      conjg(ci2a)*ck2a*int*rho2(dk2a+uk,di2a+ui)
      Exc2(cbdm+ul,cbdm+uj) = Exc2(cbdm+ul,cbdm+uj) - &
      conjg(ci2b)*ck2b*int*rho2(dk2b+uk,di2b+ui)
    end if

  end subroutine Assign_Sigma_Exc

!------------------------------------------------------------
!> apply Hess's RI technique to DKH2 2e integrals
!!
!! in'_p2 = in_p2(i,j)*mA(i,i)*mA(j,j), in(2*cbdm,2*cbdm) -> in(2*fbdm,2*fbdm)
  subroutine Assign_matrices_2e(in, mA)
    implicit none
    complex(dp),allocatable,intent(inout) :: in(:,:)
    real(dp), intent(in)                  :: mA(fbdm,fbdm)
    integer                               :: ii, jj
    complex(dp)                           :: supp(2*fbdm,2*cbdm)
    complex(dp)                           :: supp2(2*fbdm,2*fbdm)
    complex(dp)                           :: suppA(2*fbdm,2*fbdm)
    ! transform 'in' to orthogonal AO basis with suppA
    call matmul('C', 'N', exc2f, in, supp)
    call matmul('N', 'N', supp, exc2f, suppA)
    ! transform 'in' to p^2 eigenbasis with suppA
    call matmul('T', 'N', exAO2p2, suppA, supp2)
    call matmul('N', 'N', supp2, exAO2p2, suppA)
    ! construct new 'in' with suppA
    do ii = 1, fbdm
      do jj = 1, fbdm
        suppA(ii,jj) = suppA(ii,jj) * mA(ii,ii)*mA(jj,jj)
        suppA(fbdm+ii,jj) = suppA(fbdm+ii,jj) * mA(ii,ii)*mA(jj,jj)
        suppA(ii,fbdm+jj) = suppA(ii,fbdm+jj) * mA(ii,ii)*mA(jj,jj)
        suppA(fbdm+ii,fbdm+jj) = suppA(fbdm+ii,fbdm+jj) * mA(ii,ii)*mA(jj,jj)
      end do
    end do
    ! transform 'in' to orthogonal AO basis
    deallocate(in)
    allocate(in(2*fbdm,2*fbdm))
    call matmul('N', 'N', exAO2p2, suppA, supp2)
    call matmul('N', 'T', supp2, exAO2p2, in)
  end subroutine Assign_matrices_2e

!------------------------------------------------------------
!> construct non-relativistic/DKH2 spinor 1-electron Fock matrix (Fock1)
!!
!! should be called before SCF for once
  subroutine Assign_Fock_1e()
    implicit none
    real(dp)    :: temp_pool(fbdm, fbdm)
    integer     :: ii, jj               ! loop variables for Assign_Fock_1e
    complex(dp) :: itm(fbdm)
    real(dp)    :: edc(fbdm)            ! (p2+c2)^0.5
    real(dp)    :: coe                  ! coefficient for pVp-related terms
    real(dp)    :: i_V_j_p2(fbdm, fbdm) ! i_V_j in p^2 eigenbasis
    allocate(Fock1(2*fbdm,2*fbdm), exi_T_j(2*fbdm,2*fbdm), &
    exi_V_j(2*fbdm,2*fbdm), source = c0)
    if (.not. pVp1e) then
      Fock1(1:fbdm,1:fbdm) = 0.5_dp*i_p2_j + i_V_j
      Fock1(fbdm+1:2*fbdm,fbdm+1:2*fbdm) = 0.5_dp*i_p2_j + i_V_j
      exi_T_j(1:fbdm,1:fbdm) = 0.5_dp*i_p2_j
      exi_T_j(fbdm+1:2*fbdm,fbdm+1:2*fbdm) = 0.5_dp*i_p2_j
      exi_V_j(1:fbdm,1:fbdm) = i_V_j
      exi_V_j(fbdm+1:2*fbdm,fbdm+1:2*fbdm) = i_V_j
      allocate(AVA(2*fbdm,2*fbdm), source = c0)
      AVA = exi_V_j
    else
      if (.not. pVp2e) then
        write(60,'(A)') '  -- only 1e pVp enabled, we typically approximate'
        write(60,'(A)') '  -- the equivalent pVp-related matrix elements as'
        write(60,'(A)') '  -- half of the 1e pVp-related matrix elements.'
        coe = 0.5_dp
      else
        coe = 1.0_dp
      end if
      allocate(oper1(fbdm,fbdm), oper2(fbdm,fbdm), oper3(2*fbdm,2*fbdm))
      allocate(oper4(2*fbdm,2*fbdm), oper5(2*fbdm,2*fbdm), source = c0)
      allocate(Ap(fbdm,fbdm), ApRp(fbdm,fbdm), source = 0.0_dp) ! ApRp = RpAp
      allocate(SRp(2*fbdm,2*fbdm), source = c0)
      allocate(ARVRA(2*fbdm,2*fbdm), AVA(2*fbdm,2*fbdm), source = c0)
      allocate(ARVeRA(2*fbdm,2*fbdm), AVeA(2*fbdm,2*fbdm), source = c0)
      allocate(Ve(fbdm,fbdm))
      allocate(pxVepx(fbdm,fbdm), pyVepy(fbdm,fbdm), pzVepz(fbdm,fbdm))
      allocate(pxVepy(fbdm,fbdm), pyVepx(fbdm,fbdm), pxVepz(fbdm,fbdm))
      allocate(pzVepx(fbdm,fbdm), pyVepz(fbdm,fbdm), pzVepy(fbdm,fbdm))
      allocate(exAO2p2(2*fbdm,2*fbdm), expVp(2*fbdm,2*fbdm), source = c0)
      edc = dsqrt(evl_p2 + c2)
      !----------------------
      forall (ii=1:fbdm) Ap(ii, ii) = &
      dsqrt( (dsqrt(evl_p2(ii)/c2 + 1.0_dp) + 1.0_dp) / &
      (2.0_dp * dsqrt(evl_p2(ii)/c2 + 1.0_dp)) )
      !----------------------
      forall (ii=1:fbdm) ApRp(ii,ii) = Ap(ii,ii) / (edc(ii)+c)
      !----------------------
      itm = (1.0_dp + evl_p2/(4.0_dp*c2) + QED_rad) * c1
      forall (ii = 1:fbdm)
        SRp(ii, ii) = itm(ii)
        SRp(fbdm+ii, fbdm+ii)   = itm(ii)
      end forall
      !----------------------
      ! Ap V Ap
      exi_V_j(1:fbdm,1:fbdm) = i_V_j
      exi_V_j(fbdm+1:2*fbdm,fbdm+1:2*fbdm) = i_V_j
      call matmul('T', 'N', AO2p2, i_V_j, oper2)
      call matmul('N', 'N', oper2, AO2p2, i_V_j_p2)
      call matmul('N', 'N', Ap, i_V_j_p2, oper2)
      call matmul('N', 'N', oper2, Ap, oper1)
      AVA(1:fbdm,1:fbdm) = AVA(1:fbdm,1:fbdm) + oper1 * c1
      AVA(fbdm+1:2*fbdm,fbdm+1:2*fbdm) = &
      AVA(fbdm+1:2*fbdm,fbdm+1:2*fbdm) + oper1 * c1
      !----------------------
      ! ApRp pxVpx+pyVpy+pzVpz ApRp
      temp_pool = pxVpx+pyVpy+pzVpz
      call matmul('N', 'T', ApRp, AO2p2, oper2)
      call matmul('N', 'N', oper2, temp_pool, oper1)
      call matmul('N', 'N', oper1, AO2p2, oper2)
      call matmul('N', 'N', oper2, ApRp, oper1)
      ARVRA(1:fbdm,1:fbdm) = ARVRA(1:fbdm,1:fbdm) + oper1 * c1
      ARVRA(fbdm+1:2*fbdm,fbdm+1:2*fbdm) = &
      ARVRA(fbdm+1:2*fbdm,fbdm+1:2*fbdm) + oper1 * c1
      !----------------------
      ! ApRp pxVpy-pyVpx ApRp
      temp_pool = pxVpy-pyVpx
      call matmul('N', 'T', ApRp, AO2p2, oper2)
      call matmul('N', 'N', oper2, temp_pool, oper1)
      call matmul('N', 'N', oper1, AO2p2, oper2)
      call matmul('N', 'N', oper2, ApRp, oper1)
      ARVRA(1:fbdm,1:fbdm) = ARVRA(1:fbdm,1:fbdm) + oper1 * ci
      ARVRA(fbdm+1:2*fbdm,fbdm+1:2*fbdm) = &
      ARVRA(fbdm+1:2*fbdm,fbdm+1:2*fbdm) - oper1 * ci
      !----------------------
      ! ApRp pzVpx-pxVpz ApRp
      temp_pool = pzVpx-pxVpz
      call matmul('N', 'T', ApRp, AO2p2, oper2)
      call matmul('N', 'N', oper2, temp_pool, oper1)
      call matmul('N', 'N', oper1, AO2p2, oper2)
      call matmul('N', 'N', oper2, ApRp, oper1)
      ARVRA(fbdm+1:2*fbdm,1:fbdm) = ARVRA(fbdm+1:2*fbdm,1:fbdm) - oper1 * c1
      ARVRA(1:fbdm,fbdm+1:2*fbdm) = ARVRA(1:fbdm,fbdm+1:2*fbdm) + oper1 * c1
      !----------------------
      ! ApRp pyVpz-pzVpy ApRp
      temp_pool = pyVpz-pzVpy
      call matmul('N', 'T', ApRp, AO2p2, oper2)
      call matmul('N', 'N', oper2, temp_pool, oper1)
      call matmul('N', 'N', oper1, AO2p2, oper2)
      call matmul('N', 'N', oper2, ApRp, oper1)
      ARVRA(fbdm+1:2*fbdm,1:fbdm) = ARVRA(fbdm+1:2*fbdm,1:fbdm) + oper1 * ci
      ARVRA(1:fbdm,fbdm+1:2*fbdm) = ARVRA(1:fbdm,fbdm+1:2*fbdm) + oper1 * ci
      !----------------------
      Fock1 = Fock1 + AVA
      if (pppVp) then
        call matmul('N', 'N', SRp, ARVRA, oper3)
        call matmul('N', 'N', oper3, SRp, oper4)
        Fock1 = Fock1 + coe*oper4
        expVp = coe*oper4
      else
        Fock1 = Fock1 + coe*ARVRA
        expVp = coe*ARVRA
      end if
      if (pppVp) then
        allocate(exSR(2*fbdm,2*fbdm), source = c0)
        !----------------------
        ! ApRp (px3Vpx+py3Vpy+pz3Vpz+pxVpx3+pyVpy3+pzVpz3) ApRp
        temp_pool = px3Vpx+py3Vpy+pz3Vpz+transpose(px3Vpx+py3Vpy+pz3Vpz)
        call matmul('N', 'T', ApRp, AO2p2, oper2)
        call matmul('N', 'N', oper2, temp_pool, oper1)
        call matmul('N', 'N', oper1, AO2p2, oper2)
        call matmul('N', 'N', oper2, ApRp, oper1)
        Fock1(1:fbdm,1:fbdm) = &
        Fock1(1:fbdm,1:fbdm) - coe*oper1/(2.0_dp*c2) * c1
        Fock1(fbdm+1:2*fbdm,fbdm+1:2*fbdm) = &
        Fock1(fbdm+1:2*fbdm,fbdm+1:2*fbdm) - coe*oper1/(2.0_dp*c2) * c1
        exSR(1:fbdm,1:fbdm) = &
        exSR(1:fbdm,1:fbdm) - coe*oper1/(2.0_dp*c2) * c1
        exSR(fbdm+1:2*fbdm,fbdm+1:2*fbdm) = &
        exSR(fbdm+1:2*fbdm,fbdm+1:2*fbdm) - coe*oper1/(2.0_dp*c2) * c1
        !----------------------
        ! ApRp (px3Vpy-py3Vpx+pxVpy3-pyVpx3) ApRp
        temp_pool = px3Vpy-py3Vpx+transpose(py3Vpx-px3Vpy)
        call matmul('N', 'T', ApRp, AO2p2, oper2)
        call matmul('N', 'N', oper2, temp_pool, oper1)
        call matmul('N', 'N', oper1, AO2p2, oper2)
        call matmul('N', 'N', oper2, ApRp, oper1)
        Fock1(1:fbdm,1:fbdm) = &
        Fock1(1:fbdm,1:fbdm) - coe*oper1/(2.0_dp*c2) * ci
        Fock1(fbdm+1:2*fbdm,fbdm+1:2*fbdm) = &
        Fock1(fbdm+1:2*fbdm,fbdm+1:2*fbdm) + coe*oper1/(2.0_dp*c2) * ci
        exSR(1:fbdm,1:fbdm) = &
        exSR(1:fbdm,1:fbdm) - coe*oper1/(2.0_dp*c2) * ci
        exSR(fbdm+1:2*fbdm,fbdm+1:2*fbdm) = &
        exSR(fbdm+1:2*fbdm,fbdm+1:2*fbdm) + coe*oper1/(2.0_dp*c2) * ci
        !----------------------
        ! ApRp (pz3Vpx-px3Vpz+pzVpx3-pxVpz3) ApRp
        temp_pool = pz3Vpx-px3Vpz+transpose(px3Vpz-pz3Vpx)
        call matmul('N', 'T', ApRp, AO2p2, oper2)
        call matmul('N', 'N', oper2, temp_pool, oper1)
        call matmul('N', 'N', oper1, AO2p2, oper2)
        call matmul('N', 'N', oper2, ApRp, oper1)
        Fock1(fbdm+1:2*fbdm,1:fbdm) = &
        Fock1(fbdm+1:2*fbdm,1:fbdm) + coe*oper1/(2.0_dp*c2) * c1
        Fock1(1:fbdm,fbdm+1:2*fbdm) = &
        Fock1(1:fbdm,fbdm+1:2*fbdm) - coe*oper1/(2.0_dp*c2) * c1
        exSR(fbdm+1:2*fbdm,1:fbdm) = &
        exSR(fbdm+1:2*fbdm,1:fbdm) + coe*oper1/(2.0_dp*c2) * c1
        exSR(1:fbdm,fbdm+1:2*fbdm) = &
        exSR(1:fbdm,fbdm+1:2*fbdm) - coe*oper1/(2.0_dp*c2) * c1
        !----------------------
        ! ApRp (py3Vpz-pz3Vpy+pyVpz3-pzVpy3) ApRp
        temp_pool = py3Vpz-pz3Vpy+transpose(pz3Vpy-py3Vpz)
        call matmul('N', 'T', ApRp, AO2p2, oper2)
        call matmul('N', 'N', oper2, temp_pool, oper1)
        call matmul('N', 'N', oper1, AO2p2, oper2)
        call matmul('N', 'N', oper2, ApRp, oper1)
        Fock1(fbdm+1:2*fbdm,1:fbdm) = &
        Fock1(fbdm+1:2*fbdm,1:fbdm) - coe*oper1/(2.0_dp*c2) * ci
        Fock1(1:fbdm,fbdm+1:2*fbdm) = &
        Fock1(1:fbdm,fbdm+1:2*fbdm) - coe*oper1/(2.0_dp*c2) * ci
        exSR(fbdm+1:2*fbdm,1:fbdm) = &
        exSR(fbdm+1:2*fbdm,1:fbdm) - coe*oper1/(2.0_dp*c2) * ci
        exSR(1:fbdm,fbdm+1:2*fbdm) = &
        exSR(1:fbdm,fbdm+1:2*fbdm) - coe*oper1/(2.0_dp*c2) * ci
      end if
      !----------------------
      ! start building Fock1
      do jj = 1, fbdm
        do ii = 1, fbdm
          Ve(ii,jj) = i_V_j_p2(ii,jj) / (c*(edc(ii) + edc(jj)))
          pxVepx(ii,jj) = pxVpx(ii,jj) / (c*(edc(ii) + edc(jj)))
          pyVepy(ii,jj) = pyVpy(ii,jj) / (c*(edc(ii) + edc(jj)))
          pzVepz(ii,jj) = pzVpz(ii,jj) / (c*(edc(ii) + edc(jj)))
          pxVepy(ii,jj) = pxVpy(ii,jj) / (c*(edc(ii) + edc(jj)))
          pyVepx(ii,jj) = pyVpx(ii,jj) / (c*(edc(ii) + edc(jj)))
          pxVepz(ii,jj) = pxVpz(ii,jj) / (c*(edc(ii) + edc(jj)))
          pzVepx(ii,jj) = pzVpx(ii,jj) / (c*(edc(ii) + edc(jj)))
          pyVepz(ii,jj) = pyVpz(ii,jj) / (c*(edc(ii) + edc(jj)))
          pzVepy(ii,jj) = pzVpy(ii,jj) / (c*(edc(ii) + edc(jj)))
        end do
      end do
      !----------------------
      ! Ap Ve Ap
      call matmul('N', 'T', Ap, AO2p2, oper2)
      call matmul('N', 'N', oper2, Ve, oper1)
      call matmul('N', 'N', oper1, AO2p2, oper2)
      call matmul('N', 'N', oper2, Ap, oper1)
      AVeA(1:fbdm,1:fbdm) = oper1 * c1
      AVeA(fbdm+1:2*fbdm,fbdm+1:2*fbdm) = oper1 * c1
      !----------------------
      ! ApRp pxVepx+pyVepy+pzVepz ApRp
      temp_pool = pxVepx+pyVepy+pzVepz
      call matmul('N', 'T', ApRp, AO2p2, oper2)
      call matmul('N', 'N', oper2, temp_pool, oper1)
      call matmul('N', 'N', oper1, AO2p2, oper2)
      call matmul('N', 'N', oper2, ApRp, oper1)
      ARVeRA(1:fbdm,1:fbdm) = ARVeRA(1:fbdm,1:fbdm) + oper1 * c1
      ARVeRA(fbdm+1:2*fbdm,fbdm+1:2*fbdm) = &
      ARVeRA(fbdm+1:2*fbdm,fbdm+1:2*fbdm) + oper1 * c1
      !----------------------
      ! ApRp pxVepy-pyVepx ApRp
      temp_pool = pxVepy-pyVepx
      call matmul('N', 'T', ApRp, AO2p2, oper2)
      call matmul('N', 'N', oper2, temp_pool, oper1)
      call matmul('N', 'N', oper1, AO2p2, oper2)
      call matmul('N', 'N', oper2, ApRp, oper1)
      ARVeRA(1:fbdm,1:fbdm) = ARVeRA(1:fbdm,1:fbdm) + oper1 * ci
      ARVeRA(fbdm+1:2*fbdm,fbdm+1:2*fbdm) = &
      ARVeRA(fbdm+1:2*fbdm,fbdm+1:2*fbdm) - oper1 * ci
      !----------------------
      ! ApRp pzVepx-pxVepz ApRp
      temp_pool = pzVepx-pxVepz
      call matmul('N', 'T', ApRp, AO2p2, oper2)
      call matmul('N', 'N', oper2, temp_pool, oper1)
      call matmul('N', 'N', oper1, AO2p2, oper2)
      call matmul('N', 'N', oper2, ApRp, oper1)
      ARVeRA(fbdm+1:2*fbdm,1:fbdm) = ARVeRA(fbdm+1:2*fbdm,1:fbdm) - oper1 * c1
      ARVeRA(1:fbdm,fbdm+1:2*fbdm) = ARVeRA(1:fbdm,fbdm+1:2*fbdm) + oper1 * c1
      !----------------------
      ! ApRp pyVepz-pzVepy ApRp
      temp_pool = pyVepz-pzVepy
      call matmul('N', 'T', ApRp, AO2p2, oper2)
      call matmul('N', 'N', oper2, temp_pool, oper1)
      call matmul('N', 'N', oper1, AO2p2, oper2)
      call matmul('N', 'N', oper2, ApRp, oper1)
      ARVeRA(fbdm+1:2*fbdm,1:fbdm) = ARVeRA(fbdm+1:2*fbdm,1:fbdm) + oper1 * ci
      ARVeRA(1:fbdm,fbdm+1:2*fbdm) = ARVeRA(1:fbdm,fbdm+1:2*fbdm) + oper1 * ci
      !----------------------
      ! term of order c^-4, negative terms, no RI insertion
      forall (ii = 1:2*fbdm)
        oper5(ii, ii) = 0.5_dp * c1
      end forall
      call matmul('N', 'N', ARVeRA, oper5, oper3)
      call matmul('N', 'N', oper3, AVA, oper4)
      Fock1 = Fock1 - coe*oper4
      call matmul('N', 'N', ARVRA, oper5, oper3)
      call matmul('N', 'N', oper3, AVeA, oper4)
      Fock1 = Fock1 - coe*oper4
      call matmul('N', 'N', AVeA, oper5, oper3)
      call matmul('N', 'N', oper3, ARVRA, oper4)
      Fock1 = Fock1 - coe*oper4
      call matmul('N', 'N', AVA, oper5, oper3)
      call matmul('N', 'N', oper3, ARVeRA, oper4)
      Fock1 = Fock1 - coe*oper4
      !----------------------
      ! term of order c^-4, positive terms
      ! RI insertion: I = sigma.Pi(1/P^2)sigma.Pj
      ! oper5 is half of 1/P^2 or P^2, depends on whether RI insert or extract
      forall (ii = 1:fbdm) ! RI insert
        oper5(ii,ii) = 0.5_dp * ((edc(ii)+c)**2/evl_p2(ii)) * c1
        oper5(fbdm+ii,fbdm+ii) = 0.5_dp * ((edc(ii)+c)**2/evl_p2(ii)) * c1
      end forall
      call matmul('N', 'N', ARVeRA, oper5, oper3)
      call matmul('N', 'N', oper3, ARVRA, oper4)
      Fock1 = Fock1 + coe*oper4
      call matmul('N', 'N', ARVRA, oper5, oper3)
      call matmul('N', 'N', oper3, ARVeRA, oper4)
      Fock1 = Fock1 + coe*oper4
      forall (ii = 1:fbdm) ! RI extract
        oper5(ii,ii) = 0.5_dp * (evl_p2(ii)/(edc(ii)+c)**2) * c1
        oper5(fbdm+ii,fbdm+ii) = 0.5_dp * (evl_p2(ii)/(edc(ii)+c)**2) * c1
      end forall
      call matmul('N', 'N', AVeA, oper5, oper3)
      call matmul('N', 'N', oper3, AVA, oper4)
      Fock1 = Fock1 + coe*oper4
      call matmul('N', 'N', AVA, oper5, oper3)
      call matmul('N', 'N', oper3, AVeA, oper4)
      Fock1 = Fock1 + coe*oper4
      !----------------------
      ! kinetic energy
      forall (ii = 1:fbdm) ! RI extract
        Fock1(ii,ii) = &
        Fock1(ii,ii) + (c*edc(ii)-c2) * c1
        Fock1(fbdm+ii,fbdm+ii) = &
        Fock1(fbdm+ii,fbdm+ii) + (c*edc(ii)-c2) * c1
        exi_T_j(ii,ii) = &
        exi_T_j(ii,ii) + (c*edc(ii)-c2) * c1
        exi_T_j(fbdm+ii,fbdm+ii) = &
        exi_T_j(fbdm+ii,fbdm+ii) + (c*edc(ii)-c2) * c1
      end forall
      !----------------------
      ! transform from p^2 eigenbasis to orthogonal normalized AO basis
      exAO2p2(1:fbdm,1:fbdm) = AO2p2 * c1
      exAO2p2(fbdm+1:2*fbdm,fbdm+1:2*fbdm) = AO2p2 * c1
      call matmul('N', 'N', exAO2p2, Fock1, oper3)
      call matmul('N', 'T', oper3, exAO2p2, Fock1)
      call matmul('N', 'N', exAO2p2, exi_T_j, oper3)
      call matmul('N', 'T', oper3, exAO2p2, exi_T_j)
      call matmul('N', 'N', exAO2p2, expVp, oper3)
      call matmul('N', 'T', oper3, exAO2p2, expVp)
      call matmul('N', 'N', exAO2p2, AVA, oper3)
      call matmul('N', 'T', oper3, exAO2p2, AVA)
      if (pppVp) then
        call matmul('N', 'N', exAO2p2, exSR, oper3)
        call matmul('N', 'T', oper3, exAO2p2, exSR)
      end if
      deallocate(Ve, pxVepx, pyVepy, pzVepz, pxVepy, pyVepx)
      deallocate(pxVepz, pzVepx, pyVepz, pzVepy, ARVeRA, AVeA)
      deallocate(oper1, oper2, oper3, oper4, oper5)
      ! deallocate one-electron ints to reduce memory usage
      deallocate(pxVpx, pyVpy, pzVpz, pxVpy, pyVpx)
      deallocate(pxVpz, pzVpx, pyVpz, pzVpy)
      if (pppVp) then
        deallocate(px3Vpx, py3Vpy, pz3Vpz, px3Vpy, py3Vpx)
        deallocate(px3Vpz, pz3Vpx, py3Vpz, pz3Vpy)
      end if
    end if
  end subroutine Assign_Fock_1e

!------------------------------------------------------------
!> screen of non-relativistic 4-indicator integrals based on the Cauchy–Schwarz
!!
!! inequality, allocate iijj_V(cbdm,cbdm)
!!
!! <ij|V|kl>^2 <= <ij|V|ij>*<kl|V|kl>
!!
!! should be called before Assign_Fock_V2e and Assign_Fock_AVA2e
  subroutine Assign_Schwarz_V2e()
    implicit none
    integer  :: ui, uj              ! loop variables for Assign_Schwarz_V2e
    !----------------------------------
    integer  :: contri              ! contr of atomi, shelli
    integer  :: Li                  ! angular quantum number of |AOi>
    integer  :: Mi                  ! magnetic quantum number of |AOi>
    !----------------------------------
    integer  :: contrj              ! contr of atomj, shellj
    integer  :: Lj                  ! angular quantum number of |AOj>
    integer  :: Mj                  ! magnetic quantum number of |AOj>
    !DIR$ ATTRIBUTES ALIGN:align_size :: expi, expj, coei, coej
    real(dp) :: expi(20)            ! expo of |AOi>
    real(dp) :: expj(20)            ! expo of |AOj>
    real(dp) :: coei(20)            ! coefficient of |AOi>
    real(dp) :: coej(20)            ! coefficient of |AOj>
    real(dp) :: codi(3)             ! coordinate of center of |AOi>
    real(dp) :: codj(3)             ! coordinate of center of |AOj>
    if (allocated(iijj_V)) return
    write(60,'(A)') '  -- calculate Schwarz matrix <ij|V|ij>'
    allocate(iijj_V(cbdm,cbdm), source=0.0_dp)
    do ui = 1, cbdm
      contri = cbdata(ui) % contr
      Li     = cbdata(ui) % L
      Mi     = cbdata(ui) % M
      expi(1:contri) = cbdata(ui) % expo(1:contri)
      coei(1:contri) = cbdata(ui) % Ncoe(1:contri)
      codi = cbdata(ui) % pos
      do uj = ui, cbdm
        contrj = cbdata(uj) % contr
        Lj     = cbdata(uj) % L
        Mj     = cbdata(uj) % M
        expj(1:contrj) = cbdata(uj) % expo(1:contrj)
        coej(1:contrj) = cbdata(uj) % Ncoe(1:contrj)
        codj = cbdata(uj) % pos
        !---------------------------------
        ! <ij|V|ij> = (ii|V|jj)
        iijj_V(ui,uj) = iijj_V(ui,uj) + Calc_V_2e(&
        contri,contri,contrj,contrj,&
        coei,coei,coej,coej,&
        AO_fac(:,Li,Mi),AO_fac(:,Li,Mi),AO_fac(:,Lj,Mj),AO_fac(:,Lj,Mj),&
        expi,expi,expj,expj,&
        codi,codi,codj,codj)
        iijj_V(uj,ui) = iijj_V(ui,uj)
      end do
    end do
    write(60,'(A,E10.2,A)') &
    '  -- complete! cutoff:',Schwarz_VT,'; stored in iijj_V'
  end subroutine Assign_Schwarz_V2e

!------------------------------------------------------------
!> screen of DKH2 spinor 4-indicator integrals based on the Cauchy–Schwarz ineq-
!!
!! uality, allocate iijj_V(cbdm,cbdm), iijj_pVp(cbdm,cbdm), ijij_pVp(cbdm,cbdm)
!!
!! <ij|pxVpy|kl>^2 = <pxij|V|pykl>^2 <= <pxij|V|pxij>*<pykl|V|pykl>
!!
!! for pVp, we store as (pxipxi|jj), (pyipyi|jj), (pzipzi|jj)
!!
!! should be called before Assign_Fock_ARVRA2e
  subroutine Assign_Schwarz_pVp2e()
    implicit none
    integer          :: ui, uj               ! loop var for Assign_Schwarz_pVp2e
    integer          :: um, un
    integer          :: uo, up
    integer          :: ii, ii2, numi
    !----------------------------------
    integer          :: contri               ! contr of atomi, shelli
    integer          :: Li                   ! angular quantum number of |AOi>
    integer          :: Mi                   ! magnetic quantum number of |AOi>
    !----------------------------------
    integer          :: contrj               ! contr of atomj, shellj
    integer          :: Lj                   ! angular quantum number of |AOj>
    integer          :: Mj                   ! magnetic quantum number of |AOj>
    !DIR$ ATTRIBUTES ALIGN:align_size :: expi, expj, coei, coej
    real(dp)         :: expi(20)             ! expo of |AOi>
    real(dp)         :: expj(20)             ! expo of |AOj>
    real(dp)         :: coei(20)             ! coefficient of |AOi>
    real(dp)         :: coej(20)             ! coefficient of |AOj>
    real(dp)         :: codi(3)              ! coordinate of center of |AOi>
    real(dp)         :: codj(3)              ! coordinate of center of |AOj>
    integer          :: i                    ! openMP parallel variable
    !DIR$ ATTRIBUTES ALIGN:align_size :: coedx_i, coedy_i, coedz_i
    real(dp)         :: coedx_i(32)          ! coefficient.derivative x.|AOi>
    real(dp)         :: coedy_i(32)          ! coefficient.derivative y.|AOi>
    real(dp)         :: coedz_i(32)          ! coefficient.derivative z.|AOi>
    integer          :: facdx_i(3,2)         ! x,y,z factor.derivative x.|AOi>
    integer          :: facdy_i(3,2)         ! x,y,z factor.derivative y.|AOi>
    integer          :: facdz_i(3,2)         ! x,y,z factor.derivative z.|AOi>
    type threadlocal    ! thread-local storage to avoid thread-sync overhead
      real(dp),allocatable :: iijj_V(:,:), iijj_pxVpx(:,:)
      real(dp),allocatable :: iijj_pyVpy(:,:), iijj_pzVpz(:,:)
      real(dp),allocatable :: ijij_pxVpx(:,:), ijij_pyVpy(:,:)
      real(dp),allocatable :: ijij_pzVpz(:,:)
    end type
    type(threadlocal) :: tl
    if (allocated(iijj_V)) return
    if (allocated(iijj_pxVpx)) return
    if (allocated(ijij_pxVpx)) return
    write(60,'(A)') '  -- calculate Schwarz matrices <ij|V|ij>, <ij|pVp|ij>'
    allocate(iijj_V(cbdm,cbdm), iijj_pxVpx(cbdm,cbdm), source=0.0_dp)
    allocate(iijj_pyVpy(cbdm,cbdm), iijj_pzVpz(cbdm,cbdm), source=0.0_dp)
    allocate(ijij_pxVpx(cbdm,cbdm), ijij_pyVpy(cbdm,cbdm), source=0.0_dp)
    allocate(ijij_pzVpz(cbdm,cbdm), source=0.0_dp)
    !$omp parallel num_threads(threads) default(shared) private(i,ui,numi,  &
    !$omp& uj,um,un,uo,up,ii,ii2,contri,Li,Mi,contrj,Lj,Mj,expi,expj,coei,coej,&
    !$omp& codi,codj,facdx_i,facdy_i,facdz_i,coedx_i,coedy_i,coedz_i,tl)  &
    !$omp& if(threads < nproc)
    allocate(tl%iijj_V(cbdm,cbdm), tl%iijj_pxVpx(cbdm,cbdm), source=0.0_dp)
    allocate(tl%iijj_pyVpy(cbdm,cbdm), tl%iijj_pzVpz(cbdm,cbdm), source=0.0_dp)
    allocate(tl%ijij_pxVpx(cbdm,cbdm), tl%ijij_pyVpy(cbdm,cbdm), source=0.0_dp)
    allocate(tl%ijij_pzVpz(cbdm,cbdm), source=0.0_dp)
    !$omp do schedule(dynamic, 5)
    do i = 1, cbdm
      ui = i
      contri = cbdata(ui) % contr
      Li     = cbdata(ui) % L
      Mi     = cbdata(ui) % M
      expi(1:contri) = cbdata(ui) % expo(1:contri)
      coei(1:contri) = cbdata(ui) % Ncoe(1:contri)
      codi = cbdata(ui) % pos
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
      do uj = 1, cbdm
        ! no need to differentiate |AOj>
        contrj = cbdata(uj) % contr
        Lj     = cbdata(uj) % L
        Mj     = cbdata(uj) % M
        expj(1:contrj) = cbdata(uj) % expo(1:contrj)
        coej(1:contrj) = cbdata(uj) % Ncoe(1:contrj)
        codj = cbdata(uj) % pos
        !---------------------------------
        ! <ij|V|ij> = (ii|V|jj)
        tl%iijj_V(ui,uj) = Calc_V_2e(&
        contri,contri,contrj,contrj,&
        coei(1:contri),coei(1:contri),coej(1:contrj),coej(1:contrj),&
        AO_fac(:,Li,Mi),AO_fac(:,Li,Mi),AO_fac(:,Lj,Mj),AO_fac(:,Lj,Mj),&
        expi(1:contri),expi(1:contri),expj(1:contrj),expj(1:contrj),&
        codi,codi,codj,codj)
        !---------------------------------
        ! <ij|pxVpx|ij> = (pxipxi|V|jj)
        tl%iijj_pxVpx(ui,uj) = Calc_pVp_2eij(&
        AO_fac(1,Li,Mi),AO_fac(1,Li,Mi),contri,contri,contrj,contrj,&
        coedx_i(1:2*contri),coedx_i(1:2*contri),coej(1:contrj),coej(1:contrj),&
        facdx_i,facdx_i,AO_fac(:,Lj,Mj),AO_fac(:,Lj,Mj),&
        expi(1:contri),expi(1:contri),expj(1:contrj),expj(1:contrj),&
        codi,codi,codj,codj)
        !---------------------------------
        ! <ij|pyVpy|ij> = (pyipyi|V|jj)
        tl%iijj_pyVpy(ui,uj) = Calc_pVp_2eij(&
        AO_fac(2,Li,Mi),AO_fac(2,Li,Mi),contri,contri,contrj,contrj,&
        coedy_i(1:2*contri),coedy_i(1:2*contri),coej(1:contrj),coej(1:contrj),&
        facdy_i,facdy_i,AO_fac(:,Lj,Mj),AO_fac(:,Lj,Mj),&
        expi(1:contri),expi(1:contri),expj(1:contrj),expj(1:contrj),&
        codi,codi,codj,codj)
        !---------------------------------
        ! <ij|pzVpz|ij> = (pzipzi|V|jj)
        tl%iijj_pzVpz(ui,uj) = Calc_pVp_2eij(&
        AO_fac(3,Li,Mi),AO_fac(3,Li,Mi),contri,contri,contrj,contrj,&
        coedz_i(1:2*contri),coedz_i(1:2*contri),coej(1:contrj),coej(1:contrj),&
        facdz_i,facdz_i,AO_fac(:,Lj,Mj),AO_fac(:,Lj,Mj),&
        expi(1:contri),expi(1:contri),expj(1:contrj),expj(1:contrj),&
        codi,codi,codj,codj)
        !---------------------------------
        ! <ii|pxVpx|jj> = (pxij|V|pxij)
        tl%ijij_pxVpx(ui,uj) = Calc_pVp_2eik(&
        AO_fac(1,Li,Mi),AO_fac(1,Li,Mi),contri,contrj,contri,contrj,&
        coedx_i(1:2*contri),coej(1:contrj),coedx_i(1:2*contri),coej(1:contrj),&
        facdx_i,AO_fac(:,Lj,Mj),facdx_i,AO_fac(:,Lj,Mj),&
        expi(1:contri),expj(1:contrj),expi(1:contri),expj(1:contrj),&
        codi,codj,codi,codj)
        !---------------------------------
        ! <ii|pyVpy|jj> = (pyij|V|pyij)
        tl%ijij_pyVpy(ui,uj) = Calc_pVp_2eik(&
        AO_fac(2,Li,Mi),AO_fac(2,Li,Mi),contri,contrj,contri,contrj,&
        coedy_i(1:2*contri),coej(1:contrj),coedy_i(1:2*contri),coej(1:contrj),&
        facdy_i,AO_fac(:,Lj,Mj),facdy_i,AO_fac(:,Lj,Mj),&
        expi(1:contri),expj(1:contrj),expi(1:contri),expj(1:contrj),&
        codi,codj,codi,codj)
        !---------------------------------
        ! <ii|pzVpz|jj> = (pzij|V|pzij)
        tl%ijij_pzVpz(ui,uj) = Calc_pVp_2eik(&
        AO_fac(3,Li,Mi),AO_fac(3,Li,Mi),contri,contrj,contri,contrj,&
        coedz_i(1:2*contri),coej(1:contrj),coedz_i(1:2*contri),coej(1:contrj),&
        facdz_i,AO_fac(:,Lj,Mj),facdz_i,AO_fac(:,Lj,Mj),&
        expi(1:contri),expj(1:contrj),expi(1:contri),expj(1:contrj),&
        codi,codj,codi,codj)
      end do
    end do
    !$omp end do
    !-------------<thread sync>-------------
    !$omp critical
    iijj_V = iijj_V + tl%iijj_V
    ! consider the relation p = -iD, no change sign
    iijj_pxVpx = iijj_pxVpx + tl%iijj_pxVpx
    iijj_pyVpy = iijj_pyVpy + tl%iijj_pyVpy
    iijj_pzVpz = iijj_pzVpz + tl%iijj_pzVpz
    ijij_pxVpx = ijij_pxVpx + tl%ijij_pxVpx
    ijij_pyVpy = ijij_pyVpy + tl%ijij_pyVpy
    ijij_pzVpz = ijij_pzVpz + tl%ijij_pzVpz
    !$omp end critical
    ! free memory for parallel threads explicitly
    deallocate(tl%iijj_V, tl%iijj_pxVpx, tl%iijj_pyVpy, tl%iijj_pzVpz)
    deallocate(tl%ijij_pxVpx, tl%ijij_pyVpy, tl%ijij_pzVpz)
    !$omp end parallel
    write(60,'(A,E10.2,A)') &
    '  -- complete! cutoff:',Schwarz_VT,&
    '; stored in iijj_V, iijj_pVp(3 matrices), ijij_pVp(3 matrices)'
  end subroutine Assign_Schwarz_pVp2e

!------------------------------------------------------------
!> construct non-relativistic 2-electron Fock matrices
!!
!! (mHFcol, mHFexc, mKSexc, mKScor)
!!
!! "direct" calculation to avoid memory overflow and disk r&w
!!
!! should be called in each SCF iteration
  subroutine Assign_Fock_V2e()
    implicit none
    integer     :: i, j                ! for parallel computation, ui=i, uk=j
    integer     :: ui, uj              ! loop variables for Assign_Fock_V2e
    integer     :: uk, ul
    real(dp)    :: intV                ! scalar integral
    !----------------------------------
    integer     :: contri              ! contr of atom_i, shell_i
    integer     :: Li                  ! angular quantum number of |AOi>
    integer     :: Mi                  ! magnetic quantum number of |AOi>
    !----------------------------------
    integer     :: contrj              ! contr of atom_j, shell_j
    integer     :: Lj                  ! angular quantum number of |AOj>
    integer     :: Mj                  ! magnetic quantum number of |AOj>
    !----------------------------------
    integer     :: contrk              ! contr of atom_k, shell_k
    integer     :: Lk                  ! angular quantum number of |AOk>
    integer     :: Mk                  ! magnetic quantum number of |AOk>
    !----------------------------------
    integer     :: contrl              ! contr of atom_l, shell_l
    integer     :: Ll                  ! angular quantum number of |AOl>
    integer     :: Ml                  ! magnetic quantum number of |AOl>
    !DIR$ ATTRIBUTES ALIGN:align_size :: expi, expj, expk, expl
    !DIR$ ATTRIBUTES ALIGN:align_size :: coei, coej, coek, coel
    real(dp)    :: expi(20)            ! expo of |AOi>
    real(dp)    :: expj(20)            ! expo of |AOj>
    real(dp)    :: expk(20)            ! expo of |AOk>
    real(dp)    :: expl(20)            ! expo of |AOl>
    real(dp)    :: coei(20)            ! coefficient of |AOi>
    real(dp)    :: coej(20)            ! coefficient of |AOj>
    real(dp)    :: coek(20)            ! coefficient of |AOk>
    real(dp)    :: coel(20)            ! coefficient of |AOl>
    real(dp)    :: codi(3)             ! coordinate of center of |AOi>
    real(dp)    :: codj(3)             ! coordinate of center of |AOj>
    real(dp)    :: codk(3)             ! coordinate of center of |AOk>
    real(dp)    :: codl(3)             ! coordinate of center of |AOl>
    !DIR$ ATTRIBUTES ALIGN:align_size :: HFcol_mic, HFexc_mic
    complex(dp) :: HFcol_mic(2*cbdm,2*cbdm)  ! micro 2e Fock matrix
    complex(dp) :: HFexc_mic(2*cbdm,2*cbdm)  ! micro 2e Fock matrix
    complex(dp) :: supp(2*sbdm,2*cbdm)
    ! transform rho_m to Cartesian basis
    call scgo(rho_m)
    if (allocated(mHFcol)) deallocate(mHFcol)
    allocate(mHFcol(2*cbdm,2*cbdm), source=c0)
    if (allocated(mHFexc)) deallocate(mHFexc)
    allocate(mHFexc(2*cbdm,2*cbdm), source=c0)
    ! parallel zone, running results consistent with serial
    !$omp parallel num_threads(threads) default(shared) private(i,j,ui,uj,uk,&
    !$omp& ul,intV,contri,Li,Mi,contrj,Lj,Mj,contrk,Lk,Mk,contrl,Ll,Ml,&
    !$omp& expi,expj,expk,expl,coei,coej,coek,coel,codi,codj,codk,codl,&
    !$omp& HFcol_mic,HFexc_mic) if(threads < nproc)
    HFcol_mic = c0
    HFexc_mic = c0
    !$omp do schedule(dynamic,5) collapse(2)
    ! utilizing permutation symmetry:(11|22)
    ! (ij|kl)=(ji|kl)=(ij|lk)=(ji|lk)=(kl|ij)=(kl|ji)=(lk|ij)=(lk|ji)
    do i = cbdm, 1, -1
      do j = cbdm, 1, -1
        !------------------------------<ui>---------------------------------
        ui     = i
        contri = cbdata(ui) % contr
        Li     = cbdata(ui) % L
        Mi     = cbdata(ui) % M
        expi(1:contri) = cbdata(ui) % expo(1:contri)
        coei(1:contri) = cbdata(ui) % Ncoe(1:contri)
        codi = cbdata(ui) % pos
        !------------------------------<uk>---------------------------------
        uk     = j
        contrk = cbdata(uk) % contr
        Lk     = cbdata(uk) % L
        Mk     = cbdata(uk) % M
        expk(1:contrk) = cbdata(uk) % expo(1:contrk)
        coek(1:contrk) = cbdata(uk) % Ncoe(1:contrk)
        codk = cbdata(uk) % pos
        !------------------------------<uj>---------------------------------
        do uj = ui, 1, -1
          contrj = cbdata(uj) % contr
          Lj     = cbdata(uj) % L
          Mj     = cbdata(uj) % M
          expj(1:contrj) = cbdata(uj) % expo(1:contrj)
          coej(1:contrj) = cbdata(uj) % Ncoe(1:contrj)
          codj = cbdata(uj) % pos
          !----------------------------<ul>---------------------------------
          do ul = min(uk,uj+(ui*(ui-1)-uk*(uk-1))/2), 1, -1
            ! Schwarz screening, |(ij|kl)| <= dsqrt( (ii|kk)*(jj|ll) )
            if (dsqrt(iijj_V(ui,uk)*iijj_V(uj,ul)) < Schwarz_VT) cycle
            contrl = cbdata(ul) % contr
            Ll     = cbdata(ul) % L
            Ml     = cbdata(ul) % M
            expl(1:contrl) = cbdata(ul) % expo(1:contrl)
            coel(1:contrl) = cbdata(ul) % Ncoe(1:contrl)
            codl = cbdata(ul) % pos
            !===========================(ij|kl)===============================
            intV = calc_V_2e(&
            contri,contrj,contrk,contrl,&
            coei(1:contri),coej(1:contrj),coek(1:contrk),coel(1:contrl),&
            AO_fac(:,Li,Mi),AO_fac(:,Lj,Mj),AO_fac(:,Lk,Mk),AO_fac(:,Ll,Ml),&
            expi(1:contri),expj(1:contrj),expk(1:contrk),expl(1:contrl),&
            codi,codj,codk,codl)
            ! assign two-electron Fock matrices
            !------------------------<COULOMB>------------------------
            call Assign_Coulomb(HFcol_mic,rho_m,intV,ui,uj,uk,ul)
            !------------------------<EXCHANGE>------------------------
            call Assign_Exchange(HFexc_mic,rho_m,intV,ui,uj,uk,ul)
          end do
        end do
      end do
    end do
    !$omp end do
    !-------------<thread sync>-------------
    !$omp critical
    mHFcol = mHFcol + HFcol_mic
    mHFexc = mHFexc + HFexc_mic
    !$omp end critical
    !$omp end parallel
    ! assign Kohn-Sham matrices
    if (fx_id /= -1) then
      if (allocated(mKSexc)) deallocate(mKSexc)
      allocate(mKSexc(2*cbdm, 2*cbdm))
      if (fc_id /= -1) then
        if (allocated(mKScor)) deallocate(mKScor)
        allocate(mKScor(2*cbdm, 2*cbdm))
        call basis2grid_Becke(rho_m, KSexc, KScor, mKSexc, mKScor)
        call cfgo(mKSexc)
        call cfgo(mKScor)
      else
        call basis2grid_Becke(rho_m=rho_m, ex=KSexc, Fockx=mKSexc)
        call cfgo(mKSexc)
      end if
    end if
    call cfgo(mHFcol)
    call cfgo(mHFexc)
  end subroutine Assign_Fock_V2e
  
!------------------------------------------------------------
!> construct DKH2 scalar 2-electron Fock matrices
!!
!! (mHFcol, mHFexc, mKSexc, mKScor)
!!
!! "direct" calculation to avoid memory overflow and disk r&w
!!
!! should be called in each SCF iteration
  subroutine Assign_Fock_AVA2e()
    implicit none
    integer     :: i, j                ! for parallel computation, ui=i, uk=j
    integer     :: ui, uj              ! loop variables for Assign_Fock_ARVRA2e
    integer     :: uk, ul
    real(dp)    :: intV
    !----------------------------------
    integer     :: contri              ! contr of atom_i, shell_i
    integer     :: Li                  ! angular quantum number of |AOi>
    integer     :: Mi                  ! magnetic quantum number of |AOi>
    !----------------------------------
    integer     :: contrj              ! contr of atom_j, shell_j
    integer     :: Lj                  ! angular quantum number of |AOj>
    integer     :: Mj                  ! magnetic quantum number of |AOj>
    !----------------------------------
    integer     :: contrk              ! contr of atom_k, shell_k
    integer     :: Lk                  ! angular quantum number of |AOk>
    integer     :: Mk                  ! magnetic quantum number of |AOk>
    !----------------------------------
    integer     :: contrl              ! contr of atom_l, shell_l
    integer     :: Ll                  ! angular quantum number of |AOl>
    integer     :: Ml                  ! magnetic quantum number of |AOl>
    !DIR$ ATTRIBUTES ALIGN:align_size :: expi, expj, expk, expl
    !DIR$ ATTRIBUTES ALIGN:align_size :: coei, coej, coek, coel
    real(dp)    :: expi(20)            ! expo of |AOi>
    real(dp)    :: expj(20)            ! expo of |AOj>
    real(dp)    :: expk(20)            ! expo of |AOk>
    real(dp)    :: expl(20)            ! expo of |AOl>
    real(dp)    :: coei(20)            ! coefficient of |AOi>
    real(dp)    :: coej(20)            ! coefficient of |AOj>
    real(dp)    :: coek(20)            ! coefficient of |AOk>
    real(dp)    :: coel(20)            ! coefficient of |AOl>
    real(dp)    :: codi(3)             ! coordinate of center of |AOi>
    real(dp)    :: codj(3)             ! coordinate of center of |AOj>
    real(dp)    :: codk(3)             ! coordinate of center of |AOk>
    real(dp)    :: codl(3)             ! coordinate of center of |AOl>
    !DIR$ ATTRIBUTES ALIGN:align_size :: Vcol_mic, Vexc_mic
    complex(dp) :: Vcol_mic(2*cbdm,2*cbdm)     ! micro ApVAp Coulomb matrix
    complex(dp) :: Vexc_mic(2*cbdm,2*cbdm)     ! micro ApVAp Exchange matrix
    !DIR$ ATTRIBUTES ALIGN:align_size :: rho_Ap, suppff
    complex(dp),allocatable :: rho_Ap(:,:)     ! D' = D * Ap * Ap
    complex(dp) :: suppff(2*fbdm,2*fbdm)
    allocate(rho_Ap(2*sbdm,2*sbdm))
    rho_Ap = rho_m
    ! transform rho_Ap to orthogonal AO basis
    call sfgo(rho_Ap)
    ! transform rho_Ap to p^2 eigenbasis
    call matmul('T', 'N', exAO2p2, rho_Ap, suppff)
    call matmul('N', 'N', suppff, exAO2p2, rho_Ap)
    ! construct reduced density matrix
    do ui = 1, fbdm
      do uj = 1, fbdm
        rho_Ap(ui,uj) = rho_Ap(ui,uj) * Ap(ui,ui)*Ap(uj,uj)
        rho_Ap(fbdm+ui,uj) = rho_Ap(fbdm+ui,uj) * Ap(ui,ui)*Ap(uj,uj)
        rho_Ap(ui,fbdm+uj) = rho_Ap(ui,fbdm+uj) * Ap(ui,ui)*Ap(uj,uj)
        rho_Ap(fbdm+ui,fbdm+uj) = rho_Ap(fbdm+ui,fbdm+uj) * Ap(ui,ui)*Ap(uj,uj)
      end do
    end do
    ! transform rho_Ap to orthogonal AO basis
    call matmul('N', 'N', exAO2p2, rho_Ap, suppff)
    call matmul('N', 'T', suppff, exAO2p2, rho_Ap)
    ! transform rho_Ap to Cartesian basis with rho_Ap
    call fsgo(rho_Ap)
    call scgo(rho_Ap)
    ! transform rho_m to Cartesian basis with rho_m
    call scgo(rho_m)

    ! initialize 2e Fock matrices
    if (allocated(mHFcol)) deallocate(mHFcol)
    allocate(mHFcol(2*cbdm,2*cbdm), source=c0)
    if (allocated(mHFexc)) deallocate(mHFexc)
    allocate(mHFexc(2*cbdm,2*cbdm), source=c0)

    ! parallel zone, running results consistent with serial
    !$omp parallel num_threads(threads) default(shared) private(i,j,ui,uj,uk,&
    !$omp& ul,intV,contri,Li,Mi,contrj,Lj,Mj,contrk,Lk,Mk,contrl,Ll,Ml,&
    !$omp& expi,expj,expk,expl,coei,coej,coek,coel,codi,codj,codk,codl,&
    !$omp& Vcol_mic,Vexc_mic) if(threads < nproc)
    Vcol_mic = c0
    Vexc_mic = c0
    !$omp do schedule(dynamic,5) collapse(2)
    ! utilizing permutation symmetry:(11|22)
    ! (ij|kl)=(ji|kl)=(ij|lk)=(ji|lk)=(kl|ij)=(kl|ji)=(lk|ij)=(lk|ji)
    do i = cbdm, 1, -1
      do j = cbdm, 1, -1
        !------------------------------<ui>---------------------------------
        ui     = i
        contri = cbdata(ui) % contr
        Li     = cbdata(ui) % L
        Mi     = cbdata(ui) % M
        expi(1:contri) = cbdata(ui) % expo(1:contri)
        coei(1:contri) = cbdata(ui) % Ncoe(1:contri)
        codi = cbdata(ui) % pos
        !------------------------------<uk>---------------------------------
        uk     = j
        contrk = cbdata(uk) % contr
        Lk     = cbdata(uk) % L
        Mk     = cbdata(uk) % M
        expk(1:contrk) = cbdata(uk) % expo(1:contrk)
        coek(1:contrk) = cbdata(uk) % Ncoe(1:contrk)
        codk = cbdata(uk) % pos
        !------------------------------<uj>---------------------------------
        do uj = ui, 1, -1
          contrj = cbdata(uj) % contr
          Lj     = cbdata(uj) % L
          Mj     = cbdata(uj) % M
          expj(1:contrj) = cbdata(uj) % expo(1:contrj)
          coej(1:contrj) = cbdata(uj) % Ncoe(1:contrj)
          codj = cbdata(uj) % pos
          !----------------------------<ul>---------------------------------
          do ul = min(uk,uj+(ui*(ui-1)-uk*(uk-1))/2), 1, -1
            ! Schwarz screening, |(ij|kl)| <= dsqrt( (ii|kk)*(jj|ll) )
            if (dsqrt(iijj_V(ui,uk)*iijj_V(uj,ul)) < Schwarz_VT) cycle
            contrl = cbdata(ul) % contr
            Ll     = cbdata(ul) % L
            Ml     = cbdata(ul) % M
            expl(1:contrl) = cbdata(ul) % expo(1:contrl)
            coel(1:contrl) = cbdata(ul) % Ncoe(1:contrl)
            codl = cbdata(ul) % pos
            !===========================(ij|kl)===============================
            intV = calc_V_2e(&
            contri,contrj,contrk,contrl,&
            coei(1:contri),coej(1:contrj),coek(1:contrk),coel(1:contrl),&
            AO_fac(:,Li,Mi),AO_fac(:,Lj,Mj),AO_fac(:,Lk,Mk),AO_fac(:,Ll,Ml),&
            expi(1:contri),expj(1:contrj),expk(1:contrk),expl(1:contrl),&
            codi,codj,codk,codl)
            ! assign two-electron Fock matrices
            !------------------------<COULOMB>------------------------
            call Assign_Coulomb(Vcol_mic,rho_Ap,intV,ui,uj,uk,ul)
            !------------------------<EXCHANGE>------------------------
            call Assign_Exchange(Vexc_mic,rho_Ap,intV,ui,uj,uk,ul)
          end do
        end do
      end do
    end do
    !$omp end do
    !-------------<thread sync>-------------
    !$omp critical
    mHFcol = mHFcol + Vcol_mic
    mHFexc = mHFexc + Vexc_mic
    !$omp end critical
    !$omp end parallel
    call Assign_matrices_2e(mHFcol,Ap)
    call Assign_matrices_2e(mHFexc,Ap)
    ! assign Kohn-Sham matrices
    if (fx_id /= -1) then
      if (allocated(mKSexc)) deallocate(mKSexc)
      allocate(mKSexc(2*cbdm, 2*cbdm))
      if (fc_id /= -1) then
        if (allocated(mKScor)) deallocate(mKScor)
        allocate(mKScor(2*cbdm, 2*cbdm))
        call basis2grid_Becke(rho_m, KSexc, KScor, mKSexc, mKScor)
        call cfgo(mKSexc)
        call cfgo(mKScor)
      else
        call basis2grid_Becke(rho_m=rho_m, ex=KSexc, Fockx=mKSexc)
        call cfgo(mKSexc)
      end if
    end if
  end subroutine Assign_Fock_AVA2e
!------------------------------------------------------------
!> construct DKH2 spinor 2-electron Fock matrices
!!
!! (mHFcol, mHFexc, mKSexc, mKScor, mpVpexc_11,mpVpcol_11,mpVpexc_22,mpVpcol_22)
!!
!! "direct" calculation to avoid memory overflow and disk r&w
!!
!! should be called in each SCF iteration
  subroutine Assign_Fock_ARVRA2e()
    implicit none
    integer     :: i, j                ! for parallel computation, ui=i, uk=j
    integer     :: ui, uj              ! loop variables for Assign_Fock_ARVRA2e
    integer     :: uk, ul
    real(dp)    :: intV, intCpxVpx     ! scalar integrals
    real(dp)    :: intCpyVpy, intCpzVpz
    real(dp)    :: intCpxVpy, intCpxVpz! intC contributes to Coulomb matrix
    real(dp)    :: intCpyVpz, intXpxVpx! intX contributes to Exchange matrix
    real(dp)    :: intXpyVpy, intXpzVpz
    real(dp)    :: intXpxVpy, intXpxVpz
    real(dp)    :: intXpyVpz
    !----------------------------------
    integer     :: contri              ! contr of atom_i, shell_i
    integer     :: Li                  ! angular quantum number of |AOi>
    integer     :: Mi                  ! magnetic quantum number of |AOi>
    !DIR$ ATTRIBUTES ALIGN:align_size :: coedx_i, coedy_i, coedz_i
    real(dp)    :: coedx_i(32)         ! coefficient.derivative x.|AOi>
    real(dp)    :: coedy_i(32)         ! coefficient.derivative y.|AOi>
    real(dp)    :: coedz_i(32)         ! coefficient.derivative z.|AOi>
    integer     :: facdx_i(3,2)        ! x,y,z factor.derivative x.|AOi>
    integer     :: facdy_i(3,2)        ! x,y,z factor.derivative y.|AOi>
    integer     :: facdz_i(3,2)        ! x,y,z factor.derivative z.|AOi>
    !----------------------------------
    integer     :: contrj              ! contr of atom_j, shell_j
    integer     :: Lj                  ! angular quantum number of |AOj>
    integer     :: Mj                  ! magnetic quantum number of |AOj>
    !DIR$ ATTRIBUTES ALIGN:align_size :: coedx_j, coedy_j, coedz_j
    real(dp)    :: coedx_j(32)         ! coefficient.derivative x.|AOj>
    real(dp)    :: coedy_j(32)         ! coefficient.derivative y.|AOj>
    real(dp)    :: coedz_j(32)         ! coefficient.derivative z.|AOj>
    integer     :: facdx_j(3,2)        ! x,y,z factor.derivative x.|AOj>
    integer     :: facdy_j(3,2)        ! x,y,z factor.derivative y.|AOj>
    integer     :: facdz_j(3,2)        ! x,y,z factor.derivative z.|AOj>
    !----------------------------------
    integer     :: contrk              ! contr of atom_k, shell_k
    integer     :: Lk                  ! angular quantum number of |AOk>
    integer     :: Mk                  ! magnetic quantum number of |AOk>
    !DIR$ ATTRIBUTES ALIGN:align_size :: coedx_k, coedy_k, coedz_k
    real(dp)    :: coedx_k(32)         ! coefficient.derivative x.|AOk>
    real(dp)    :: coedy_k(32)         ! coefficient.derivative y.|AOk>
    real(dp)    :: coedz_k(32)         ! coefficient.derivative z.|AOk>
    integer     :: facdx_k(3,2)        ! x,y,z factor.derivative x.|AOk>
    integer     :: facdy_k(3,2)        ! x,y,z factor.derivative y.|AOk>
    integer     :: facdz_k(3,2)        ! x,y,z factor.derivative z.|AOk>
    !----------------------------------
    integer     :: contrl              ! contr of atom_l, shell_l
    integer     :: Ll                  ! angular quantum number of |AOl>
    integer     :: Ml                  ! magnetic quantum number of |AOl>
    !DIR$ ATTRIBUTES ALIGN:align_size :: expi, expj, expk, expl
    !DIR$ ATTRIBUTES ALIGN:align_size :: coei, coej, coek, coel
    real(dp)    :: expi(20)            ! expo of |AOi>
    real(dp)    :: expj(20)            ! expo of |AOj>
    real(dp)    :: expk(20)            ! expo of |AOk>
    real(dp)    :: expl(20)            ! expo of |AOl>
    real(dp)    :: coei(20)            ! coefficient of |AOi>
    real(dp)    :: coej(20)            ! coefficient of |AOj>
    real(dp)    :: coek(20)            ! coefficient of |AOk>
    real(dp)    :: coel(20)            ! coefficient of |AOl>
    real(dp)    :: codi(3)             ! coordinate of center of |AOi>
    real(dp)    :: codj(3)             ! coordinate of center of |AOj>
    real(dp)    :: codk(3)             ! coordinate of center of |AOk>
    real(dp)    :: codl(3)             ! coordinate of center of |AOl>
    !DIR$ ATTRIBUTES ALIGN:align_size :: Vcol_mic, Vexc_mic
    !DIR$ ATTRIBUTES ALIGN:align_size :: R1VR1col_mic, R1VR1exc_mic
    !DIR$ ATTRIBUTES ALIGN:align_size :: R2VR2col_mic, R2VR2exc_mic
    !DIR$ ATTRIBUTES ALIGN:align_size :: rho_Ap, rho_ApRp, suppff
    complex(dp) :: Vcol_mic(2*cbdm,2*cbdm)     ! micro ApVAp Coulomb matrix
    complex(dp) :: Vexc_mic(2*cbdm,2*cbdm)     ! micro ApVAp Exchange matrix
    complex(dp) :: R1VR1col_mic(2*cbdm,2*cbdm) ! micro ApRiVRiAp Coulomb matrix
    complex(dp) :: R1VR1exc_mic(2*cbdm,2*cbdm) ! micro ApRiVRiAp Exchange matrix
    complex(dp) :: R2VR2col_mic(2*cbdm,2*cbdm) ! micro ApRjVRjAp Coulomb matrix
    complex(dp) :: R2VR2exc_mic(2*cbdm,2*cbdm) ! micro ApRjVRjAp Exchange matrix
    complex(dp) :: suppff(2*fbdm,2*fbdm)
    complex(dp),allocatable :: rho_Ap(:,:)       ! D' = D * Ap * Ap
    complex(dp),allocatable :: rho_ApRp(:,:)     ! D' = D * ApRp * ApRp
    allocate(rho_Ap(2*sbdm,2*sbdm))
    allocate(rho_ApRp(2*fbdm,2*fbdm))
    rho_Ap = rho_m
    ! transform rho_Ap to orthogonal AO basis
    call sfgo(rho_Ap)
    ! transform rho_Ap to p^2 eigenbasis
    call matmul('T', 'N', exAO2p2, rho_Ap, suppff)
    call matmul('N', 'N', suppff, exAO2p2, rho_Ap)
    ! construct reduced density matrix
    rho_ApRp = rho_Ap
    do ui = 1, fbdm
      do uj = 1, fbdm
        rho_ApRp(ui,uj) = rho_ApRp(ui,uj) * ApRp(ui,ui)*ApRp(uj,uj)
        rho_ApRp(fbdm+ui,uj) = rho_ApRp(fbdm+ui,uj) * ApRp(ui,ui)*ApRp(uj,uj)
        rho_ApRp(ui,fbdm+uj) = rho_ApRp(ui,fbdm+uj) * ApRp(ui,ui)*ApRp(uj,uj)
        rho_ApRp(fbdm+ui,fbdm+uj) = rho_ApRp(fbdm+ui,fbdm+uj) * &
        ApRp(ui,ui)*ApRp(uj,uj)
        rho_Ap(ui,uj) = rho_Ap(ui,uj) * Ap(ui,ui)*Ap(uj,uj)
        rho_Ap(fbdm+ui,uj) = rho_Ap(fbdm+ui,uj) * Ap(ui,ui)*Ap(uj,uj)
        rho_Ap(ui,fbdm+uj) = rho_Ap(ui,fbdm+uj) * Ap(ui,ui)*Ap(uj,uj)
        rho_Ap(fbdm+ui,fbdm+uj) = rho_Ap(fbdm+ui,fbdm+uj) * Ap(ui,ui)*Ap(uj,uj)
      end do
    end do
    ! transform rho_Ap and rho_ApRp to orthogonal AO basis
    call matmul('N', 'N', exAO2p2, rho_ApRp, suppff)
    call matmul('N', 'T', suppff, exAO2p2, rho_ApRp)
    call matmul('N', 'N', exAO2p2, rho_Ap, suppff)
    call matmul('N', 'T', suppff, exAO2p2, rho_Ap)
    ! transform rho_Ap and rho_ApRp to Cartesian basis
    call fsgo(rho_Ap)
    call scgo(rho_Ap)
    call fsgo(rho_ApRp)
    call scgo(rho_ApRp)
    ! transform rho_m to Cartesian basis
    call scgo(rho_m)

    ! initialize 2e Fock matrices
    if (allocated(mHFcol)) deallocate(mHFcol)
    allocate(mHFcol(2*cbdm,2*cbdm), source=c0)
    if (allocated(mHFexc)) deallocate(mHFexc)
    allocate(mHFexc(2*cbdm,2*cbdm), source=c0)
    if (allocated(mpVpcol_11)) deallocate(mpVpcol_11)
    allocate(mpVpcol_11(2*cbdm,2*cbdm), source=c0)
    if (allocated(mpVpexc_11)) deallocate(mpVpexc_11)
    allocate(mpVpexc_11(2*cbdm,2*cbdm), source=c0)
    if (allocated(mpVpcol_22)) deallocate(mpVpcol_22)
    allocate(mpVpcol_22(2*cbdm,2*cbdm), source=c0)
    if (allocated(mpVpexc_22)) deallocate(mpVpexc_22)
    allocate(mpVpexc_22(2*cbdm,2*cbdm), source=c0)

    ! parallel zone, running results consistent with serial
    !$omp parallel num_threads(threads) default(shared) private(i,j,ui,uj,uk,&
    !$omp& ul,intV,intCpxVpx,intCpyVpy,intCpzVpz,intCpxVpy,intCpxVpz,intCpyVpz,&
    !$omp& intXpxVpx,intXpyVpy,intXpzVpz,intXpxVpy,intXpxVpz,intXpyVpz,&
    !$omp& contri,Li,Mi,contrj,Lj,Mj,contrk,Lk,Mk,contrl,Ll,Ml,&
    !$omp& expi,expj,expk,expl,coei,coej,coek,coel,codi,codj,codk,codl,&
    !$omp& coedx_i,coedy_i,coedz_i,facdx_i,facdy_i,facdz_i,coedx_j,coedy_j,&
    !$omp& coedz_j,facdx_j,facdy_j,facdz_j,coedx_k,coedy_k,coedz_k,facdx_k,&
    !$omp& facdy_k,facdz_k,Vcol_mic,Vexc_mic,R1VR1col_mic,R1VR1exc_mic,&
    !$omp& R2VR2col_mic,R2VR2exc_mic) if(threads < nproc)
    Vcol_mic = c0
    Vexc_mic = c0
    R1VR1col_mic = c0
    R1VR1exc_mic = c0
    R2VR2col_mic = c0
    R2VR2exc_mic = c0
    !$omp do schedule(dynamic,5) collapse(2)
    ! have to calculate all (ij|kl) to get SOC terms correctly
    do i = 1, cbdm
      do j = 1, cbdm
        !------------------------------<ui>---------------------------------
        ui     = i
        contri = cbdata(ui) % contr
        Li     = cbdata(ui) % L
        Mi     = cbdata(ui) % M
        expi(1:contri) = cbdata(ui) % expo(1:contri)
        coei(1:contri) = cbdata(ui) % Ncoe(1:contri)
        codi = cbdata(ui) % pos
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
        !------------------------------<uk>---------------------------------
        uk     = j
        contrk = cbdata(uk) % contr
        Lk     = cbdata(uk) % L
        Mk     = cbdata(uk) % M
        expk(1:contrk) = cbdata(uk) % expo(1:contrk)
        coek(1:contrk) = cbdata(uk) % Ncoe(1:contrk)
        codk = cbdata(uk) % pos
        !---------------------------------------
        ! factor of x^m*d(exp)
        facdx_k(:,1) = AO_fac(:,Lk,Mk)
        facdx_k(1,1) = facdx_k(1,1) + 1
        facdy_k(:,1) = AO_fac(:,Lk,Mk)
        facdy_k(2,1) = facdy_k(2,1) + 1
        facdz_k(:,1) = AO_fac(:,Lk,Mk)
        facdz_k(3,1) = facdz_k(3,1) + 1
        !-------------------------------
        ! factor of d(x^m)*exp
        facdx_k(:,2) = AO_fac(:,Lk,Mk)
        facdx_k(1,2) = max(facdx_k(1,2)-1,0)
        facdy_k(:,2) = AO_fac(:,Lk,Mk)
        facdy_k(2,2) = max(facdy_k(2,2)-1,0)
        facdz_k(:,2) = AO_fac(:,Lk,Mk)
        facdz_k(3,2) = max(facdz_k(3,2)-1,0)
        !---------------------------------------
        ! coefficient of x^m*d(exp)
        coedx_k(1:contrk) = -2.0_dp * coek(1:contrk) * expk(1:contrk)
        coedy_k(1:contrk) = -2.0_dp * coek(1:contrk) * expk(1:contrk)
        coedz_k(1:contrk) = -2.0_dp * coek(1:contrk) * expk(1:contrk)
        !---------------------------------
        ! coefficient of d(x^m)*exp
        coedx_k(contrk+1:2*contrk) = AO_fac(1,Lk,Mk) * coek(1:contrk)
        coedy_k(contrk+1:2*contrk) = AO_fac(2,Lk,Mk) * coek(1:contrk)
        coedz_k(contrk+1:2*contrk) = AO_fac(3,Lk,Mk) * coek(1:contrk)
        !------------------------------<uj>---------------------------------
        do uj = 1, cbdm
          contrj = cbdata(uj) % contr
          Lj     = cbdata(uj) % L
          Mj     = cbdata(uj) % M
          expj(1:contrj) = cbdata(uj) % expo(1:contrj)
          coej(1:contrj) = cbdata(uj) % Ncoe(1:contrj)
          codj = cbdata(uj) % pos
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
          !----------------------------<ul>---------------------------------
          do ul = 1, cbdm
            contrl = cbdata(ul) % contr
            Ll     = cbdata(ul) % L
            Ml     = cbdata(ul) % M
            expl(1:contrl) = cbdata(ul) % expo(1:contrl)
            coel(1:contrl) = cbdata(ul) % Ncoe(1:contrl)
            codl = cbdata(ul) % pos
            !===========================(ij|kl)===========================
            if (uj<=ui .and. ul<=uk .and. ul-uj<=(ui*(ui-1)-uk*(uk-1))/2) then
              if (dsqrt(iijj_V(ui,uk)*iijj_V(uj,ul)) > Schwarz_VT) then
                intV = calc_V_2e(&
                contri,contrj,contrk,contrl,&
                coei(1:contri),coej(1:contrj),coek(1:contrk),coel(1:contrl),&
                AO_fac(:,Li,Mi),AO_fac(:,Lj,Mj),&
                AO_fac(:,Lk,Mk),AO_fac(:,Ll,Ml),&
                expi(1:contri),expj(1:contrj),expk(1:contrk),expl(1:contrl),&
                codi,codj,codk,codl)
                !------------------------<V-COULOMB>------------------------
                call Assign_Coulomb(Vcol_mic,rho_Ap,intV,ui,uj,uk,ul)
                !------------------------<V-EXCHANGE>------------------------
                call Assign_Exchange(Vexc_mic,rho_Ap,intV,ui,uj,uk,ul)
              end if
            end if
            if (ul <= uk) then
              if (uj <= ui) then
                !==========================(pxipxj|kl)==========================
                if (dsqrt(iijj_pxVpx(ui,uk)*iijj_pxVpx(uj,ul)) > Schwarz_VT)then
                  intCpxVpx = Calc_pVp_2eij(&
                  AO_fac(1,Li,Mi),AO_fac(1,Lj,Mj),contri,contrj,contrk,contrl,&
                  coedx_i(1:2*contri),coedx_j(1:2*contrj),&
                  coek(1:contrk),coel(1:contrl),&
                  facdx_i,facdx_j,AO_fac(:,Lk,Mk),AO_fac(:,Ll,Ml),&
                  expi(1:contri),expj(1:contrj),expk(1:contrk),expl(1:contrl),&
                  codi,codj,codk,codl)
                  !----------------------<RVR-COULOMB>----------------------
                  call Assign_Sigma_Col(1,1,R1VR1col_mic,R2VR2col_mic,rho_Ap,&
                  rho_ApRp,intCpxVpx,ui,uj,uk,ul)
                end if
                !==========================(pyipyj|kl)==========================
                if (dsqrt(iijj_pyVpy(ui,uk)*iijj_pyVpy(uj,ul)) > Schwarz_VT)then
                  intCpyVpy = Calc_pVp_2eij(&
                  AO_fac(2,Li,Mi),AO_fac(2,Lj,Mj),contri,contrj,contrk,contrl,&
                  coedy_i(1:2*contri),coedy_j(1:2*contrj),&
                  coek(1:contrk),coel(1:contrl),&
                  facdy_i,facdy_j,AO_fac(:,Lk,Mk),AO_fac(:,Ll,Ml),&
                  expi(1:contri),expj(1:contrj),expk(1:contrk),expl(1:contrl),&
                  codi,codj,codk,codl)
                  !----------------------<RVR-COULOMB>----------------------
                  call Assign_Sigma_Col(2,2,R1VR1col_mic,R2VR2col_mic,rho_Ap,&
                  rho_ApRp,intCpyVpy,ui,uj,uk,ul)
                end if
                !==========================(pzipzj|kl)==========================
                if (dsqrt(iijj_pzVpz(ui,uk)*iijj_pzVpz(uj,ul)) > Schwarz_VT)then
                  intCpzVpz = Calc_pVp_2eij(&
                  AO_fac(3,Li,Mi),AO_fac(3,Lj,Mj),contri,contrj,contrk,contrl,&
                  coedz_i(1:2*contri),coedz_j(1:2*contrj),&
                  coek(1:contrk),coel(1:contrl),&
                  facdz_i,facdz_j,AO_fac(:,Lk,Mk),AO_fac(:,Ll,Ml),&
                  expi(1:contri),expj(1:contrj),expk(1:contrk),expl(1:contrl),&
                  codi,codj,codk,codl)
                  !----------------------<RVR-COULOMB>----------------------
                  call Assign_Sigma_Col(3,3,R1VR1col_mic,R2VR2col_mic,rho_Ap,&
                  rho_ApRp,intCpzVpz,ui,uj,uk,ul)
                end if
              end if
              !===========================(pxipyj|kl)===========================
              if (dsqrt(iijj_pxVpx(ui,uk)*iijj_pyVpy(uj,ul)) > Schwarz_VT) then
                intCpxVpy = Calc_pVp_2eij(&
                AO_fac(1,Li,Mi),AO_fac(2,Lj,Mj),contri,contrj,contrk,contrl,&
                coedx_i(1:2*contri),coedy_j(1:2*contrj),&
                coek(1:contrk),coel(1:contrl),&
                facdx_i,facdy_j,AO_fac(:,Lk,Mk),AO_fac(:,Ll,Ml),&
                expi(1:contri),expj(1:contrj),expk(1:contrk),expl(1:contrl),&
                codi,codj,codk,codl)
                !----------------------<RVR-COULOMB>----------------------
                call Assign_Sigma_Col(1,2,R1VR1col_mic,R2VR2col_mic,rho_Ap,&
                rho_ApRp,intCpxVpy,ui,uj,uk,ul)
              end if
              !===========================(pxipzj|kl)===========================
              if (dsqrt(iijj_pxVpx(ui,uk)*iijj_pzVpz(uj,ul)) > Schwarz_VT) then
                intCpxVpz = Calc_pVp_2eij(&
                AO_fac(1,Li,Mi),AO_fac(3,Lj,Mj),contri,contrj,contrk,contrl,&
                coedx_i(1:2*contri),coedz_j(1:2*contrj),&
                coek(1:contrk),coel(1:contrl),&
                facdx_i,facdz_j,AO_fac(:,Lk,Mk),AO_fac(:,Ll,Ml),&
                expi(1:contri),expj(1:contrj),expk(1:contrk),expl(1:contrl),&
                codi,codj,codk,codl)
                !----------------------<RVR-COULOMB>----------------------
                call Assign_Sigma_Col(1,3,R1VR1col_mic,R2VR2col_mic,rho_Ap,&
                rho_ApRp,intCpxVpz,ui,uj,uk,ul)
              end if
              !===========================(pyipzj|kl)===========================
              if (dsqrt(iijj_pyVpy(ui,uk)*iijj_pzVpz(uj,ul)) > Schwarz_VT) then
                intCpyVpz = Calc_pVp_2eij(&
                AO_fac(2,Li,Mi),AO_fac(3,Lj,Mj),contri,contrj,contrk,contrl,&
                coedy_i(1:2*contri),coedz_j(1:2*contrj),&
                coek(1:contrk),coel(1:contrl),&
                facdy_i,facdz_j,AO_fac(:,Lk,Mk),AO_fac(:,Ll,Ml),&
                expi(1:contri),expj(1:contrj),expk(1:contrk),expl(1:contrl),&
                codi,codj,codk,codl)
                !----------------------<RVR-COULOMB>----------------------
                call Assign_Sigma_Col(2,3,R1VR1col_mic,R2VR2col_mic,rho_Ap,&
                rho_ApRp,intCpyVpz,ui,uj,uk,ul)
              end if
            end if
            if (uk <= ui) then
              !===========================(pxij|pxkl)===========================
              if (dsqrt(ijij_pxVpx(ui,uj)*ijij_pxVpx(uk,ul)) > Schwarz_VT) then
                intXpxVpx = Calc_pVp_2eik(&
                AO_fac(1,Li,Mi),AO_fac(1,Lk,Mk),contri,contrj,contrk,contrl,&
                coedx_i(1:2*contri),coej(1:contrj),&
                coedx_k(1:2*contrk),coel(1:contrl),&
                facdx_i,AO_fac(:,Lj,Mj),facdx_k,AO_fac(:,Ll,Ml),&
                expi(1:contri),expj(1:contrj),expk(1:contrk),expl(1:contrl),&
                codi,codj,codk,codl)
                !----------------------<RVR-EXCHANGE>----------------------
                call Assign_Sigma_Exc(1,1,R1VR1exc_mic,R2VR2exc_mic,rho_Ap,&
                rho_ApRp,intXpxVpx,ui,uj,uk,ul)
              end if
              !===========================(pyij|pykl)===========================
              if (dsqrt(ijij_pyVpy(ui,uj)*ijij_pyVpy(uk,ul)) > Schwarz_VT) then
                intXpyVpy = Calc_pVp_2eik(&
                AO_fac(2,Li,Mi),AO_fac(2,Lk,Mk),contri,contrj,contrk,contrl,&
                coedy_i(1:2*contri),coej(1:contrj),&
                coedy_k(1:2*contrk),coel(1:contrl),&
                facdy_i,AO_fac(:,Lj,Mj),facdy_k,AO_fac(:,Ll,Ml),&
                expi(1:contri),expj(1:contrj),expk(1:contrk),expl(1:contrl),&
                codi,codj,codk,codl)
                !----------------------<RVR-EXCHANGE>----------------------
                call Assign_Sigma_Exc(2,2,R1VR1exc_mic,R2VR2exc_mic,rho_Ap,&
                rho_ApRp,intXpyVpy,ui,uj,uk,ul)
              end if
              !===========================(pzij|pzkl)===========================
              if (dsqrt(ijij_pzVpz(ui,uj)*ijij_pzVpz(uk,ul)) > Schwarz_VT) then
                intXpzVpz = Calc_pVp_2eik(&
                AO_fac(3,Li,Mi),AO_fac(3,Lk,Mk),contri,contrj,contrk,contrl,&
                coedz_i(1:2*contri),coej(1:contrj),&
                coedz_k(1:2*contrk),coel(1:contrl),&
                facdz_i,AO_fac(:,Lj,Mj),facdz_k,AO_fac(:,Ll,Ml),&
                expi(1:contri),expj(1:contrj),expk(1:contrk),expl(1:contrl),&
                codi,codj,codk,codl)
                !----------------------<RVR-EXCHANGE>----------------------
                call Assign_Sigma_Exc(3,3,R1VR1exc_mic,R2VR2exc_mic,rho_Ap,&
                rho_ApRp,intXpzVpz,ui,uj,uk,ul)
              end if
            end if
            !===========================(pxij|pykl)===========================
            if (dsqrt(ijij_pxVpx(ui,uj)*ijij_pyVpy(uk,ul)) > Schwarz_VT) then
              intXpxVpy = Calc_pVp_2eik(&
              AO_fac(1,Li,Mi),AO_fac(2,Lk,Mk),contri,contrj,contrk,contrl,&
              coedx_i(1:2*contri),coej(1:contrj),&
              coedy_k(1:2*contrk),coel(1:contrl),&
              facdx_i,AO_fac(:,Lj,Mj),facdy_k,AO_fac(:,Ll,Ml),&
              expi(1:contri),expj(1:contrj),expk(1:contrk),expl(1:contrl),&
              codi,codj,codk,codl)
              !----------------------<RVR-EXCHANGE>----------------------
              call Assign_Sigma_Exc(1,2,R1VR1exc_mic,R2VR2exc_mic,rho_Ap,&
              rho_ApRp,intXpxVpy,ui,uj,uk,ul)
            end if
            !===========================(pxij|pzkl)===========================
            if (dsqrt(ijij_pxVpx(ui,uj)*ijij_pzVpz(uk,ul)) > Schwarz_VT) then
              intXpxVpz = Calc_pVp_2eik(&
              AO_fac(1,Li,Mi),AO_fac(3,Lk,Mk),contri,contrj,contrk,contrl,&
              coedx_i(1:2*contri),coej(1:contrj),&
              coedz_k(1:2*contrk),coel(1:contrl),&
              facdx_i,AO_fac(:,Lj,Mj),facdz_k,AO_fac(:,Ll,Ml),&
              expi(1:contri),expj(1:contrj),expk(1:contrk),expl(1:contrl),&
              codi,codj,codk,codl)
              !----------------------<RVR-EXCHANGE>----------------------
              call Assign_Sigma_Exc(1,3,R1VR1exc_mic,R2VR2exc_mic,rho_Ap,&
              rho_ApRp,intXpxVpz,ui,uj,uk,ul)
            end if
            !===========================(pyij|pzkl)===========================
            if (dsqrt(ijij_pyVpy(ui,uj)*ijij_pzVpz(uk,ul)) > Schwarz_VT) then
              intXpyVpz = Calc_pVp_2eik(&
              AO_fac(2,Li,Mi),AO_fac(3,Lk,Mk),contri,contrj,contrk,contrl,&
              coedy_i(1:2*contri),coej(1:contrj),&
              coedz_k(1:2*contrk),coel(1:contrl),&
              facdy_i,AO_fac(:,Lj,Mj),facdz_k,AO_fac(:,Ll,Ml),&
              expi(1:contri),expj(1:contrj),expk(1:contrk),expl(1:contrl),&
              codi,codj,codk,codl)
              !----------------------<RVR-EXCHANGE>----------------------
              call Assign_Sigma_Exc(2,3,R1VR1exc_mic,R2VR2exc_mic,rho_Ap,&
              rho_ApRp,intXpyVpz,ui,uj,uk,ul)
            end if
          end do
        end do
      end do
    end do
    !$omp end do
    !-------------<thread sync>-------------
    !$omp critical
    mHFcol = mHFcol + Vcol_mic
    mHFexc = mHFexc + Vexc_mic
    mpVpcol_11 = mpVpcol_11 + R1VR1col_mic
    mpVpexc_11 = mpVpexc_11 + R1VR1exc_mic
    mpVpcol_22 = mpVpcol_22 + R2VR2col_mic
    mpVpexc_22 = mpVpexc_22 + R2VR2exc_mic
    !$omp end critical
    !$omp end parallel
    call Assign_matrices_2e(mHFcol,Ap)
    call Assign_matrices_2e(mHFexc,Ap)
    call Assign_matrices_2e(mpVpcol_11,ApRp)
    call Assign_matrices_2e(mpVpexc_11,ApRp)
    call Assign_matrices_2e(mpVpcol_22,Ap)
    call Assign_matrices_2e(mpVpexc_22,Ap)
    ! assign Kohn-Sham matrices
    if (fx_id /= -1) then
      if (allocated(mKSexc)) deallocate(mKSexc)
      allocate(mKSexc(2*cbdm, 2*cbdm))
      if (fc_id /= -1) then
        if (allocated(mKScor)) deallocate(mKScor)
        allocate(mKScor(2*cbdm, 2*cbdm))
        call basis2grid_Becke(rho_m, KSexc, KScor, mKSexc, mKScor)
        call cfgo(mKSexc)
        call cfgo(mKScor)
      else
        call basis2grid_Becke(rho_m=rho_m, ex=KSexc, Fockx=mKSexc)
        call cfgo(mKSexc)
      end if
    end if
  end subroutine Assign_Fock_ARVRA2e

!------------------------------------------------------------
!> calculate <S**2> based on oper3 generated by SCF
!!
!! divide all occupied orbitals into two sets of alpha, beta orbitals with
!!
!! their original coefficients, then use UHF method to solve Lowdin <S**2>.
  subroutine Calc_S2HF()
    implicit none
    integer     :: ii, jj, kk                ! loop variables for Calc_S2HF
    complex(dp) :: alal, albe, beal, bebe    ! components of pair density matrix
    complex(dp) :: suppa, suppb
    ! oper3 and i_j are orthogonally normalized
    totalpha = 0.0_dp
    totbeta = 0.0_dp
    do ii = 1, electron_count
      do kk = 1, fbdm
        totalpha = totalpha + &
        real(conjg(oper3(kk,occindex(ii))) * &
        oper3(kk,occindex(ii)),dp)
      end do
      do kk = fbdm+1, 2*fbdm
        totbeta = totbeta + &
        real(conjg(oper3(kk,occindex(ii))) * &
        oper3(kk,occindex(ii)),dp)
      end do
    end do
    alal = totalpha * (totalpha-1.0_dp)
    bebe = totbeta * (totbeta-1.0_dp)
    albe = c0
    beal = c0
    do ii = 1, electron_count
      do jj = 1, electron_count
        suppa = c0
        suppb = c0
        do kk = 1, fbdm
          suppa = suppa + &
          conjg(oper3(kk,occindex(ii))) * &
          oper3(fbdm+kk,occindex(jj))
          suppb = suppb + &
          conjg(oper3(fbdm+kk,occindex(jj))) * &
          oper3(kk,occindex(ii))
        end do
        albe = albe + suppa*suppb
        suppa = c0
        suppb = c0
        do kk = fbdm+1, 2*fbdm
          suppa = suppa + &
          conjg(oper3(kk,occindex(ii))) * &
          oper3(kk-fbdm,occindex(jj))
          suppb = suppb + &
          conjg(oper3(kk-fbdm,occindex(jj))) * &
          oper3(kk,occindex(ii))
        end do
        beal = beal + suppa*suppb
      end do
    end do
    S__2 = -real(electron_count*(electron_count-4),dp)/4.0_dp + &
    real(alal/2.0_dp + bebe/2.0_dp) - real(albe/2.0_dp + beal/2.0_dp)
  end subroutine Calc_S2HF
  
!------------------------------------------------------------
!> calculate <S**2> for orbitals based on oper3 generated by SCF
  subroutine Calc_S2HForb(orbnum)
    implicit none
    integer,intent(in) :: orbnum
    integer            :: ii       !loop variables for Calc_S2HForb
    real(dp)           :: suppa, suppb
    totalphaorb = 0.0_dp
    do ii = 1, fbdm
      totalphaorb = totalphaorb + &
      real(conjg(oper3(ii,orbnum))*oper3(ii,orbnum),dp)
    end do
    totbetaorb = 0.0_dp
    do ii = fbdm+1, 2*fbdm
      totbetaorb = totbetaorb + &
      real(conjg(oper3(ii,orbnum))*oper3(ii,orbnum),dp)
    end do
    suppa = c0
    suppb = c0
    do ii = fbdm+1, 2*fbdm
      suppa = suppa + conjg(oper3(ii,orbnum))*oper3(ii-fbdm,orbnum)
      suppb = suppb + conjg(oper3(ii-fbdm,orbnum))*oper3(ii,orbnum)
    end do
    S__2orb = 3.0/4.0 + totalphaorb*(totalphaorb-1.0)/2.0 + &
    totbetaorb*(totbetaorb-1.0)/2.0 - real(suppa*suppb)
    Szorb = (totalphaorb-totbetaorb)/2.0_dp
  end subroutine Calc_S2HForb

!------------------------------------------------------------
!> spin projection via Rotation Group Integration -- 
!! to get spin pure states form 2-componnet states
!!
!! projection operator: P^S_MK = (2S+1)/(8pi^2)*Int(D^S_MK*(omega)*R(omega))
!!
!! ref: Percus, J. K., and A. Rotenberg. J Math Phys, 3.5 (1962): 928-932.
!!
!! SU(2) rotation is based on z-y-z convention of Euler angles
subroutine RGI()
  use Lebedev
  implicit none

  ! SU(2) rotation matrix(z-y-z convention):
  ! R(phi,theta,chi) = exp(-(i/2)*phi*sigma_z) * 
  ! exp(-(i/2)*theta*sigma_y) * exp(-(i/2)*chi*sigma_z)

  ! Wigner D-matrix
  ! D^S_MK = exp(-i*phi*M)*d^S_MK*exp(-i*chi*K)

  real(dp)                :: mcs2           ! mc: basical configuration
  integer                 :: BCSmult        ! 1: singlet, 2: doublet ...
  integer                 :: BCSz           ! 0: sz=0, 1: sz=1/2 ...
  integer                 :: ii, jj         ! loop variables
  real(dp)                :: phi, theta, chi! Euler angles
  real(dp)                :: wsmalld(6)     ! Wigner d-matrix
  complex(dp)             :: wbigd(6)       ! Wigner D-matrix
  complex(dp)             :: spinproj(6)    ! <psi|P^S_MK|psi>
  complex(dp)             :: orb_i(2*fbdm,electron_count)
  complex(dp)             :: orb_o(2*fbdm,electron_count)
  complex(dp)        :: overlap(electron_count,electron_count)! <MO|R(omega)|MO>
  complex(dp)             :: rotproj        ! <psi|R(omega)|psi>
  real(dp)           :: codx(434),cody(434),codz(434) ! Cartesian Lebedev points
  real(dp)                :: weight(434)    ! Lebedev weights
  real(dp)           :: phis(434), thetas(434)! Spherical Lebedev points
  integer                 :: n
  real(dp)                :: R

  write(60,'(A)') &
  '  ============================================================='
  write(60,'(A)') '                   Rotation Group Integration'
  write(60,'(A)') &
  '  ============================================================='
  if (.not. allocated(oper3)) then
    write(60,'(A)') '  RGI: oper3 not allocated, unable to perform'
    return
  end if

  ! designation of basical configuration, BCSmult and BCSz
  if (S__2 < 0.0_dp) S__2 = -S__2  ! negative S__2 arises from numerical error
  if (mod(electron_count,2) == 0) then
    do BCSmult = 0, 50, 2
      mcs2 = real(BCSmult)/2.0*(real(BCSmult)/2.0+1.0)
      if (S__2 < mcs2) exit
    end do
  else
    do BCSmult = 1, 51, 2
      mcs2 = real(BCSmult)/2.0*(real(BCSmult)/2.0+1.0)
      if (S__2 < mcs2) exit
    end do
  end if
  BCSmult = BCSmult - 1
  if (mod(electron_count,2)/=0 .and. S__2 < 0.75) BCSmult = 2
  BCSz = nint(totalpha) - nint(totbeta)
  if (mod(BCSmult,2) == 0 .and. mod(BCSz,2) == 0) then
    write(60,'(A)') '  RGI: BCSmult & BCSz are contradictory, unable to perform'
    return
  else if (mod(BCSmult,2) /= 0 .and. mod(BCSz,2) /= 0) then
    write(60,'(A)') '  RGI: BCSmult & BCSz are contradictory, unable to perform'
    return
  end if
  write(60,'(A,I2,A,I2,A)') '  basical S_mult =', bcsmult,', S_z =',bcsz,'/2'

  ! Lebedev gird points
  call LD0434(codx, cody, codz, weight, n)
  weight = weight * 4.0_dp*pi ! sphere weight
  do ii = 1, n
    ! polar angle theta (0 ≤ theta ≤ pi)
    thetas(ii) = acos(codz(ii))
    ! azimuth angle phi (0 ≤ phi < 2pi)
    R = hypot(codx(ii), cody(ii))
    if (R > 1.0E-12_dp) then
      ! not on the z-axis
      phis(ii) = atan2(cody(ii), codx(ii))
      if (phis(ii) < 0.0_dp) phis(ii) = phis(ii) + 2.0_dp*pi
    else
      ! on the z-axis（theta=0 or theta=pi），phi set to 0 to avoid NaN
      phis(ii) = 0.0_dp
    end if
  end do

  ! Consider configuration components satisfies abs(S_z^a - S_z^b) <= 1 and 
  ! abs(S^a - S^b) <= 1. According to Wigner–Eckart theorem, the CG coefficients
  ! (transition matrix elements) between states are not zero only if they
  ! satisfies the 2 conditions. For more see RGI_Configurations
  spinproj = c0
  orb_i = oper3(:,1:electron_count)
  ! To construct the pure states of S^2 and S_z, it is necessary to first obtain
  ! the components with the same S_z from the projected 2-component state, and
  ! then perform rotation projection. 
  ! For components with the same S^2 but different S_z, converting them to the
  ! same S_z using the Wigner D-matrix will cause mutual contamination between
  ! sub-states.
  ! Therefore, we need to set K = M for all possible M.
  if (BCSmult == 1 .and. BCSz == 0) then ! basical configuration is closed shell
    !---------------
    ! S     M
    ! 0     0
    ! 2     2
    !       0
    !      -2
    !---------------
    do ii = 1, n
      phi = phis(ii)
      theta = thetas(ii)
      wsmalld = 0.0_dp
      ! projection to (0,0)
      wsmalld(1) = wigner_d(0, 0, 0, theta)
      ! projection to (1,1)
      wsmalld(2) = wigner_d(2, 2, 2, theta)
      ! projection to (1,0)
      wsmalld(3) = wigner_d(2, 0, 0, theta)
      ! projection to (1,-1)
      wsmalld(4) = wigner_d(2, -2, -2, theta)
      ! 100 chi points are verified to provide sufficient integration accuracy
      do jj = 0, 99
        wbigd = c0
        chi = 2.0_dp*pi*real(jj,dp)/(100.0_dp)
        wbigd(1) = wsmalld(1)
        wbigd(2) = wsmalld(2) * exp(ci*chi) * exp(ci*phi)
        wbigd(3) = wsmalld(3)
        wbigd(4) = wsmalld(4) * exp(-ci*chi) * exp(-ci*phi)
        call orb_SU2trafo(orb_i, orb_o, phi, theta, chi)
        call matmul('C', 'N', orb_i, orb_o, overlap)
        ! <psi|R(omega)|psi> = det(<MO|R(omega)|MO>) = det(overlap)
        rotproj = det(overlap)
        spinproj(1:4) = spinproj(1:4) + &
        weight(ii)*2.0_dp*pi/(100.0_dp)*wbigd(1:4)*rotproj
      end do
    end do
    spinproj(1) = spinproj(1) * (1.0_dp/(8.0_dp*pi**2))
    spinproj(2:4) = spinproj(2:4) * (3.0_dp/(8.0_dp*pi**2))
    write(60,'(A,E12.5)') &
    '  |S^0| (basical configuration)                 ...', &
    dsqrt(real(spinproj(1)*conjg(spinproj(1))))
    write(60,'(A,E12.5)') &
    '  |T^1|                                         ...', &
    dsqrt(real(spinproj(2)*conjg(spinproj(2))))
    write(60,'(A,E12.5)') &
    '  |T^0|                                         ...', &
    dsqrt(real(spinproj(3)*conjg(spinproj(3))))
    write(60,'(A,E12.5)') &
    '  |T^-1|                                        ...', &
    dsqrt(real(spinproj(4)*conjg(spinproj(4))))
  else if (BCSmult == 2 .and. BCSz == 1) then ! basical configuration is doublet
    !---------------
    ! 2S    2M
    !  1     1
    !       -1
    !  3     3
    !        1
    !       -1
    !---------------
    do ii = 1, n
      phi = phis(ii)
      theta = thetas(ii)
      wsmalld = 0.0_dp
      ! projection to (1/2,1/2)
      wsmalld(1) = wigner_d(1, 1, 1, theta)
      ! projection to (1/2,-1/2)
      wsmalld(2) = wigner_d(1, -1, -1, theta)
      ! projection to (3/2,3/2)
      wsmalld(3) = wigner_d(3, 3, 3, theta)
      ! projection to (3/2,1/2)
      wsmalld(4) = wigner_d(3, 1, 1, theta)
      ! projection to (3/2,-1/2)
      wsmalld(5) = wigner_d(3, -1, -1, theta)
      ! 100 chi points are verified to provide sufficient integration accuracy
      do jj = 0, 99
        wbigd = c0
        chi = 2.0_dp*pi*real(jj,dp)/(100.0_dp)
        wbigd(1) = wsmalld(1) * exp(0.5_dp*ci*chi) * exp(0.5_dp*ci*phi)
        wbigd(2) = wsmalld(2) * exp(-0.5_dp*ci*chi) * exp(-0.5_dp*ci*phi)
        wbigd(3) = wsmalld(3) * exp(1.5_dp*ci*chi) * exp(1.5_dp*ci*phi)
        wbigd(4) = wsmalld(4) * exp(0.5_dp*ci*chi) * exp(0.5_dp*ci*phi)
        wbigd(5) = wsmalld(5) * exp(-0.5_dp*ci*chi) * exp(-0.5_dp*ci*phi)
        call orb_SU2trafo(orb_i, orb_o, phi, theta, chi)
        call matmul('C', 'N', orb_i, orb_o, overlap)
        ! <psi|R(omega)|psi> = det(<MO|R(omega)|MO>) = det(overlap)
        rotproj = det(overlap)
        spinproj(1:5) = spinproj(1:5) + &
        weight(ii)*2.0_dp*pi/(100.0_dp)*wbigd(1:5)*rotproj
      end do
    end do
    spinproj(1:2) = spinproj(1:2) * (2.0_dp/(8.0_dp*pi**2))
    spinproj(3:5) = spinproj(3:5) * (4.0_dp/(8.0_dp*pi**2))
    write(60,'(A,E12.5)') &
    '  |D^1/2| (basical configuration)               ...', &
    dsqrt(real(spinproj(1)*conjg(spinproj(1))))
    write(60,'(A,E12.5)') &
    '  |D^-1/2|                                      ...', &
    dsqrt(real(spinproj(2)*conjg(spinproj(2))))
    write(60,'(A,E12.5)') &
    '  |4^3/2|                                       ...', &
    dsqrt(real(spinproj(3)*conjg(spinproj(3))))
    write(60,'(A,E12.5)') &
    '  |4^1/2|                                       ...', &
    dsqrt(real(spinproj(4)*conjg(spinproj(4))))
    write(60,'(A,E12.5)') &
    '  |4^-1/2|                                      ...', &
    dsqrt(real(spinproj(5)*conjg(spinproj(5))))
  else if (BCSmult == 3 .and. BCSz == 2) then ! basical configuration is triplet
    !---------------
    ! 2S    2M
    !  0     0
    !  2     2
    !        0
    !  4     4
    !        2
    !        0
    !---------------
    do ii = 1, n
      phi = phis(ii)
      theta = thetas(ii)
      wsmalld = 0.0_dp
      ! projection to (0,0)
      wsmalld(1) = wigner_d(0, 0, 0, theta)
      ! projection to (1,1)
      wsmalld(2) = wigner_d(2, 2, 2, theta)
      ! projection to (1,0)
      wsmalld(3) = wigner_d(2, 0, 0, theta)
      ! projection to (2,2)
      wsmalld(4) = wigner_d(4, 4, 4, theta)
      ! projection to (2,1)
      wsmalld(5) = wigner_d(4, 2, 2, theta)
      ! projection to (2,0)
      wsmalld(6) = wigner_d(4, 0, 0, theta)
      ! 100 chi points are verified to provide sufficient integration accuracy
      do jj = 0, 99
        wbigd = c0
        chi = 2.0_dp*pi*real(jj,dp)/(100.0_dp)
        wbigd(1) = wsmalld(1)
        wbigd(2) = wsmalld(2) * exp(ci*chi) * exp(ci*phi)
        wbigd(3) = wsmalld(3)
        wbigd(4) = wsmalld(4) * exp(2.0_dp*ci*chi) * exp(2.0_dp*ci*phi)
        wbigd(5) = wsmalld(5) * exp(ci*chi) * exp(ci*phi)
        wbigd(6) = wsmalld(6)
        call orb_SU2trafo(orb_i, orb_o, phi, theta, chi)
        call matmul('C', 'N', orb_i, orb_o, overlap)
        ! <psi|R(omega)|psi> = det(<MO|R(omega)|MO>) = det(overlap)
        rotproj = det(overlap)
        spinproj = spinproj + weight(ii)*2.0_dp*pi/(100.0_dp)*wbigd*rotproj
      end do
    end do
    spinproj(1) = spinproj(1) * (1.0_dp/(8.0_dp*pi**2))
    spinproj(2:3) = spinproj(2:3) * (3.0_dp/(8.0_dp*pi**2))
    spinproj(4:6) = spinproj(4:6) * (5.0_dp/(8.0_dp*pi**2))
    write(60,'(A,E12.5)') &
    '  |S^0|                                         ...', &
    dsqrt(real(spinproj(1)*conjg(spinproj(1))))
    write(60,'(A,E12.5)') &
    '  |T^1|  (basical configuration)                ...', &
    dsqrt(real(spinproj(2)*conjg(spinproj(2))))
    write(60,'(A,E12.5)') &
    '  |T^0|                                         ...', &
    dsqrt(real(spinproj(3)*conjg(spinproj(3))))
    write(60,'(A,E12.5)') &
    '  |5^2|                                         ...', &
    dsqrt(real(spinproj(4)*conjg(spinproj(4))))
    write(60,'(A,E12.5)') &
    '  |5^1|                                         ...', &
    dsqrt(real(spinproj(5)*conjg(spinproj(5))))
    write(60,'(A,E12.5)') &
    '  |5^0|                                         ...', &
    dsqrt(real(spinproj(6)*conjg(spinproj(6))))
  else if (BCSmult == 4 .and. BCSz == 3) then ! basical configuration is quartet
    !---------------
    ! 2S    2M
    !  1     1
    !  3     3
    !        1
    !  5     5
    !        3
    !        1
    !---------------
    do ii = 1, n
      phi = phis(ii)
      theta = thetas(ii)
      wsmalld = 0.0_dp
      ! projection to (1/2,1/2)
      wsmalld(1) = wigner_d(1, 1, 1, theta)
      ! projection to (3/2,3/2)
      wsmalld(2) = wigner_d(3, 3, 3, theta)
      ! projection to (3/2,1/2)
      wsmalld(3) = wigner_d(3, 1, 1, theta)
      ! projection to (5/2,5/2)
      wsmalld(4) = wigner_d(5, 5, 5, theta)
      ! projection to (5/2,3/2)
      wsmalld(5) = wigner_d(5, 3, 3, theta)
      ! projection to (5/2,1/2)
      wsmalld(6) = wigner_d(5, 1, 1, theta)
      ! 100 chi points are verified to provide sufficient integration accuracy
      do jj = 0, 99
        wbigd = c0
        chi = 2.0_dp*pi*real(jj,dp)/(100.0_dp)
        wbigd(1) = wsmalld(1) * exp(0.5_dp*ci*chi) * exp(0.5_dp*ci*phi)
        wbigd(2) = wsmalld(2) * exp(1.5_dp*ci*chi) * exp(1.5_dp*ci*phi)
        wbigd(3) = wsmalld(3) * exp(0.5_dp*ci*chi) * exp(0.5_dp*ci*phi)
        wbigd(4) = wsmalld(4) * exp(2.5_dp*ci*chi) * exp(2.5_dp*ci*phi)
        wbigd(5) = wsmalld(5) * exp(1.5_dp*ci*chi) * exp(1.5_dp*ci*phi)
        wbigd(6) = wsmalld(6) * exp(0.5_dp*ci*chi) * exp(0.5_dp*ci*phi)
        call orb_SU2trafo(orb_i, orb_o, phi, theta, chi)
        call matmul('C', 'N', orb_i, orb_o, overlap)
        ! <psi|R(omega)|psi> = det(<MO|R(omega)|MO>) = det(overlap)
        rotproj = det(overlap)
        spinproj = spinproj + weight(ii)*2.0_dp*pi/(100.0_dp)*wbigd*rotproj
      end do
    end do
    spinproj(1) = spinproj(1) * (2.0_dp/(8.0_dp*pi**2))
    spinproj(2:3) = spinproj(2:3) * (4.0_dp/(8.0_dp*pi**2))
    spinproj(4:6) = spinproj(4:6) * (6.0_dp/(8.0_dp*pi**2))
    write(60,'(A,E12.5)') &
    '  |D^1/2|                                       ...', &
    dsqrt(real(spinproj(1)*conjg(spinproj(1))))
    write(60,'(A,E12.5)') &
    '  |4^3/2|  (basical configuration)              ...', &
    dsqrt(real(spinproj(2)*conjg(spinproj(2))))
    write(60,'(A,E12.5)') &
    '  |4^1/2|                                       ...', &
    dsqrt(real(spinproj(3)*conjg(spinproj(3))))
    write(60,'(A,E12.5)') &
    '  |6^5/2|                                       ...', &
    dsqrt(real(spinproj(4)*conjg(spinproj(4))))
    write(60,'(A,E12.5)') &
    '  |6^3/2|                                       ...', &
    dsqrt(real(spinproj(5)*conjg(spinproj(5))))
    write(60,'(A,E12.5)') &
    '  |6^1/2|                                       ...', &
    dsqrt(real(spinproj(6)*conjg(spinproj(6))))
  else if (BCSmult == 5 .and. BCSz == 4) then ! basical configuration is quintet
    !---------------
    ! 2S    2M
    !  2     2
    !  4     4
    !        2
    !  6     6
    !        4
    !        2
    !---------------
    do ii = 1, n
      phi = phis(ii)
      theta = thetas(ii)
      wsmalld = 0.0_dp
      ! projection to (1,1)
      wsmalld(1) = wigner_d(2, 2, 2, theta)
      ! projection to (2,2)
      wsmalld(2) = wigner_d(4, 4, 4, theta)
      ! projection to (2,1)
      wsmalld(3) = wigner_d(4, 2, 2, theta)
      ! projection to (3,3)
      wsmalld(4) = wigner_d(6, 6, 6, theta)
      ! projection to (3,2)
      wsmalld(5) = wigner_d(6, 4, 4, theta)
      ! projection to (3,1)
      wsmalld(6) = wigner_d(6, 2, 2, theta)
      ! 100 chi points are verified to provide sufficient integration accuracy
      do jj = 0, 99
        wbigd = c0
        chi = 2.0_dp*pi*real(jj,dp)/(100.0_dp)
        wbigd(1) = wsmalld(1) * exp(ci*chi) * exp(ci*phi)
        wbigd(2) = wsmalld(2) * exp(2.0_dp*ci*chi) * exp(2.0_dp*ci*phi)
        wbigd(3) = wsmalld(3) * exp(ci*chi) * exp(ci*phi)
        wbigd(4) = wsmalld(4) * exp(3.0_dp*ci*chi) * exp(3.0_dp*ci*phi)
        wbigd(5) = wsmalld(5) * exp(2.0_dp*ci*chi) * exp(2.0_dp*ci*phi)
        wbigd(6) = wsmalld(6) * exp(ci*chi) * exp(ci*phi)
        call orb_SU2trafo(orb_i, orb_o, phi, theta, chi)
        call matmul('C', 'N', orb_i, orb_o, overlap)
        ! <psi|R(omega)|psi> = det(<MO|R(omega)|MO>) = det(overlap)
        rotproj = det(overlap)
        spinproj = spinproj + weight(ii)*2.0_dp*pi/(100.0_dp)*wbigd*rotproj
      end do
    end do
    spinproj(1) = spinproj(1) * (3.0_dp/(8.0_dp*pi**2))
    spinproj(2:3) = spinproj(2:3) * (5.0_dp/(8.0_dp*pi**2))
    spinproj(4:6) = spinproj(4:6) * (7.0_dp/(8.0_dp*pi**2))
    write(60,'(A,E12.5)') &
    '  |T^1|                                         ...', &
    dsqrt(real(spinproj(1)*conjg(spinproj(1))))
    write(60,'(A,E12.5)') &
    '  |5^2|  (basical configuration)                ...', &
    dsqrt(real(spinproj(2)*conjg(spinproj(2))))
    write(60,'(A,E12.5)') &
    '  |5^1|                                         ...', &
    dsqrt(real(spinproj(3)*conjg(spinproj(3))))
    write(60,'(A,E12.5)') &
    '  |7^3|                                         ...', &
    dsqrt(real(spinproj(4)*conjg(spinproj(4))))
    write(60,'(A,E12.5)') &
    '  |7^2|                                         ...', &
    dsqrt(real(spinproj(5)*conjg(spinproj(5))))
    write(60,'(A,E12.5)') &
    '  |7^1|                                         ...', &
    dsqrt(real(spinproj(6)*conjg(spinproj(6))))
  else if (BCSmult == 6 .and. BCSz == 5) then ! basical configuration is sextet
    !---------------
    ! 2S    2M
    !  3     3
    !  5     5
    !        3
    !  7     7
    !        5
    !        3
    !---------------
    do ii = 1, n
      phi = phis(ii)
      theta = thetas(ii)
      wsmalld = 0.0_dp
      ! projection to (3/2,3/2)
      wsmalld(1) = wigner_d(3, 3, 3, theta)
      ! projection to (5/2,5/2)
      wsmalld(2) = wigner_d(5, 5, 5, theta)
      ! projection to (5/2,3/2)
      wsmalld(3) = wigner_d(5, 3, 3, theta)
      ! projection to (7/2,7/2)
      wsmalld(4) = wigner_d(7, 7, 7, theta)
      ! projection to (7/2,5/2)
      wsmalld(5) = wigner_d(7, 5, 5, theta)
      ! projection to (7/2,3/2)
      wsmalld(6) = wigner_d(7, 3, 3, theta)
      ! 100 chi points are verified to provide sufficient integration accuracy
      do jj = 0, 99
        wbigd = c0
        chi = 2.0_dp*pi*real(jj,dp)/(100.0_dp)
        wbigd(1) = wsmalld(1) * exp(1.5_dp*ci*chi) * exp(1.5_dp*ci*phi)
        wbigd(2) = wsmalld(2) * exp(2.5_dp*ci*chi) * exp(2.5_dp*ci*phi)
        wbigd(3) = wsmalld(3) * exp(1.5_dp*ci*chi) * exp(1.5_dp*ci*phi)
        wbigd(4) = wsmalld(4) * exp(3.5_dp*ci*chi) * exp(3.5_dp*ci*phi)
        wbigd(5) = wsmalld(5) * exp(2.5_dp*ci*chi) * exp(2.5_dp*ci*phi)
        wbigd(6) = wsmalld(6) * exp(1.5_dp*ci*chi) * exp(1.5_dp*ci*phi)
        call orb_SU2trafo(orb_i, orb_o, phi, theta, chi)
        call matmul('C', 'N', orb_i, orb_o, overlap)
        ! <psi|R(omega)|psi> = det(<MO|R(omega)|MO>) = det(overlap)
        rotproj = det(overlap)
        spinproj = spinproj + weight(ii)*2.0_dp*pi/(100.0_dp)*wbigd*rotproj
      end do
    end do
    spinproj(1) = spinproj(1) * (4.0_dp/(8.0_dp*pi**2))
    spinproj(2:3) = spinproj(2:3) * (6.0_dp/(8.0_dp*pi**2))
    spinproj(4:6) = spinproj(4:6) * (8.0_dp/(8.0_dp*pi**2))
    write(60,'(A,E12.5)') &
    '  |4^3/2|                                       ...', &
    dsqrt(real(spinproj(1)*conjg(spinproj(1))))
    write(60,'(A,E12.5)') &
    '  |6^5/2|  (basical configuration)              ...', &
    dsqrt(real(spinproj(2)*conjg(spinproj(2))))
    write(60,'(A,E12.5)') &
    '  |6^3/2|                                       ...', &
    dsqrt(real(spinproj(3)*conjg(spinproj(3))))
    write(60,'(A,E12.5)') &
    '  |8^7/2|                                       ...', &
    dsqrt(real(spinproj(4)*conjg(spinproj(4))))
    write(60,'(A,E12.5)') &
    '  |8^5/2|                                       ...', &
    dsqrt(real(spinproj(5)*conjg(spinproj(5))))
    write(60,'(A,E12.5)') &
    '  |8^3/2|                                       ...', &
    dsqrt(real(spinproj(6)*conjg(spinproj(6))))
  else if (BCSmult > BCSz+1) then  ! spin polarization open-shell
    write(60,'(A)') &
    '  RGI: basical configuration is spin polarization, in this'
    write(60,'(A)') &
    '       case, decompose single-determinant into spin pure'
    write(60,'(A)') &
    '       states is meaningless.'
    return
  else
    write(60,'(A)') '  RGI: BCSmult < BCSz+1, unable to perform'
    return
  end if
  write(60,'(A)') &
  '  ============================================================='
end subroutine RGI

!-----------------------------------------------------------------------
!> calculate molecule deviations from time-reversal symmetry
!!
!! (Krammers degeneration)
!!
!! kappa = norm(k2.conjg(k2)+I)
!!
!! k2 = <MO_i|-i*sigma_y|MO_j>
real(dp) function Krammers() result(kappa)
  implicit none
  integer           :: ii            ! loop variables
  complex(dp)       :: ci_j(fbdm, fbdm)
  complex(dp)       :: k1(fbdm, electron_count)
  complex(dp)       :: k1TR(fbdm, electron_count)
  complex(dp)       :: supp(electron_count, electron_count)
  complex(dp)       :: k2(electron_count, electron_count)
  complex(dp)       :: conjgk2(electron_count, electron_count)
  complex(dp)       :: k2k2pI(electron_count, electron_count)
  kappa = 0.0_dp
  ci_j = i_j * c1
  k1 = oper3(1:fbdm,1:electron_count)
  k1TR = -conjg(oper3(fbdm+1:2*fbdm,1:electron_count))
  call matmul('C', 'N', k1, k1TR, supp)
  k2 = supp

  k1 = oper3(fbdm+1:2*fbdm,1:electron_count)
  k1TR = conjg(oper3(1:fbdm,1:electron_count))
  call matmul('C', 'N', k1, k1TR, supp)
  k2 = k2 + supp

  conjgk2 = conjg(k2)
  call matmul('N', 'N', k2, conjgk2, k2k2pI)

  forall(ii=1:electron_count) k2k2pI(ii,ii) = k2k2pI(ii,ii) + c1
  kappa = norm(k2k2pI)
end function Krammers

!-----------------------------------------------------------------------
!> calculate dispersion correction by Grimme's DFT-D4
real(dp) function DFTD4()
  implicit none
  character(len=1)  :: ch1
  character(len=10) :: ch10
  character(len=30) :: ch30
  write(ch1,"(I1)") charge
  write(60,'(A)') '  calling DFT-D4 for dispersion correction'
  ! take .xyz file as input file of DFT-D4
  ios = system('dftd4 '//trim(address_xyz)//' -f '//trim(funcemd4)//&
  ' -c '//ch1//' --noedisp --json '//trim(address_job)//'.emd4 -s -s')
  if (ios == -1) then
    write(60,'(A)') '  DFT-D4 failed with error calling, emd4 set zero'
    DFTD4 = 0.0
    return
  end if
  open(20,file=trim(address_job)//'.emd4',&
  status='old',action='read',iostat=ios)
  if (ios /= 0) then
    write(60,'(A)') '  DFT-D4 failed with empty result, emd4 set zero'
    DFTD4 = 0.0
    return
  else
    read(20,*) ! read '{'
    do
      read(20, *, iostat = ios) ch10, ch1, ch30  ! ch1 for ':'
      if (ios /= 0) then
        write(60,'(A)') '  no energy print in json, emd4 set zero'
        DFTD4 = 0.0
        close(20)
        return
      else if (index(ch10, 'energy') /= 0) then
        read(ch30, '(E30.20)') DFTD4
        if (DFTD4 > 0.0) then
          write(60,'(A)') &
          '  get non-negative d4 correction, dispersion energy set zero'
          DFTD4 = 0.0
        else
          write(60,'(A,F10.5,A)') '  complete! emd4 = ',DFTD4,' Eh'
        end if
        close(20)
        return
      end if
    end do
  end if
end function DFTD4

!-----------------------------------------------------------------------
!> project MO coefficient from m_basis to basis
!!
!! cB = i_j_s^(-1).M_i_j.cA
!!
!! cA(M_sbdm, M_sbdm) -> cB(sbdm, M_sbdm) number of project MOs is M_sbdm
  subroutine M_basis_proj(cA, cB)
    implicit none
    real(dp),intent(in)  :: cA(M_sbdm, M_sbdm)
    real(dp),intent(out) :: cB(sbdm, M_sbdm)
    integer              :: ii, jj, kk, mm
    integer              :: contri, contrj
    integer              :: Li, Lj
    integer              :: Mi, Mj
    real(dp)             :: expi(16), expj(16)
    real(dp)             :: coei(16), coej(16)
    real(dp)             :: codi(3), codj(3)
    real(dp)             :: spp(sbdm,M_sbdm)
    real(dp)             :: i_j_inv(sbdm,sbdm)
    real(dp)             :: M_min_evl
    if (.not. allocated(M_i_j)) then
      ! assign M_i_j
      allocate(M_i_j(cbdm,M_cbdm), source=0.0_dp)
      do ii = 1, cbdm
        contri = cbdata(ii) % contr
        Li     = cbdata(ii) % L
        Mi     = cbdata(ii) % M
        expi(1:contri) = cbdata(ii) % expo(1:contri)
        coei(1:contri) = cbdata(ii) % Ncoe(1:contri)
        codi = cbdata(ii) % pos
        do jj = 1, M_cbdm
          contrj = M_cbdata(jj) % contr
          Lj     = M_cbdata(jj) % L
          Mj     = M_cbdata(jj) % M
          expj(1:contrj) = M_cbdata(jj) % expo(1:contrj)
          coej(1:contrj) = M_cbdata(jj) % Ncoe(1:contrj)
          codj = M_cbdata(jj) % pos
          do kk = 1, contri
            do mm = 1, contrj
              M_i_j(ii,jj) = M_i_j(ii,jj) + &
              Integral_S_1e(coei(kk)*coej(mm),&
              AO_fac(:,Li,Mi),M_AO_fac(:,Lj,Mj),&
              expi(kk),expj(mm),codi,codj)
            end do
          end do
        end do
      end do
      if (.not. allocated(M_c2s)) call M_assign_cs()
      call M_csgo(M_i_j)
    end if
    i_j_inv = i_j_s
    call inverse(i_j_inv, sbdm)
    ! cB = i_j_s^(-1).M_i_j.cA, cB is not orthogonal, and no need to
    call matmul('N', 'N', M_i_j, cA, spp)
    call matmul('N', 'N', i_j_inv, spp, cB)
  end subroutine M_basis_proj

!-----------------------------------------------------------------------
!> dump spinor orbitals to .molden.d
  subroutine Dump_MOLDEN()
    implicit none
    character(len=200)     :: dir
    integer                :: channel
    integer                :: dmi, dmj, dmk     ! loop variables
    dir = trim(address_job)//'.molden.d'
    call execute_command_line('mkdir -p '//trim(dir), wait=.true., exitstat=ios)
    if (ios /= 0) call terminate(&
    "dump to MOLDEN failed, molden.d can't be created")
    if (.not. allocated(AO2MO)) &
    call terminate('dump to MOLDEN failed, AO2MO is empty')
    write(60,'(A)') "  AO2MO is like:"
    write(60,'(A)') "  MO 1 2 3      ...       fbdm         ...       2*fbdm"
    write(60,'(A)') "  AO1                       |                       |"
    write(60,'(A)') "  AO2                       |                       |"
    write(60,'(A)') "    |                       |                       |"
    write(60,'(A)') "  ...        part1          |         part2         |"
    write(60,'(A)') "    |    alpha real&img     |    alpha real&img     |"
    write(60,'(A)') "    |                       |                       |"
    write(60,'(A)') "    |                       |                       |"
    write(60,'(A)') "  sbdm______________________|_______________________|"
    write(60,'(A)') "    |                       |                       |"
    write(60,'(A)') "    |                       |                       |"
    write(60,'(A)') "    |                       |                       |"
    write(60,'(A)') "  ...         part1         |         part2         |"
    write(60,'(A)') "    |     beta real&img     |     beta real&img     |"
    write(60,'(A)') "    |                       |                       |"
    write(60,'(A)') "    |                       |                       |"
    write(60,'(A)') "  2*sbdm____________________|_______________________|"
    
    ! MOLDEN file contains the real part1 of AO2MO
    open(newunit=channel, file=trim(dir)//'/'//trim(jobname)//&
    '-realpart1.molden.input', status='replace', action='write', iostat=ios)
    if (ios /= 0) call terminate('creat .molden.input failed')
    write(channel, '(A)') '[Molden Format]'
    write(channel, '(A)') '[Title]'
    write(channel, '(A)') &
    'generated by TRESC, real part1 of MOs of job '//trim(address_job)
    write(channel, *)
    ! mol geometry
    write(channel, '(A)') '[Atoms] AU'
    do dmi = 1, atom_count
      write(channel, '(A,I4,I4,F16.10,F16.10,F16.10)') &
      element_list(mol(dmi)%atom_number), dmi, &
      mol(dmi)%atom_number, mol(dmi)%pos(1), &
      mol(dmi)%pos(2), mol(dmi)%pos(3)
    end do
    ! basis of each atom
    write(channel, '(A)') '[GTO]'
    do dmi = 1, atom_count
      write(channel, '(I3,I2)') dmi, 0
      do dmj = 0, shell_in_element(mol(dmi)%atom_number)-1
        if (atom_basis(mol(dmi)%basis_number+dmj)%L == 0) then
          write(channel, '(A2,I3,A)') 's', &
          atom_basis(mol(dmi)%basis_number+dmj)%contr,' 1.0'
        else if (atom_basis(mol(dmi)%basis_number+dmj)%L == 1) then
          write(channel, '(A2,I3,A)') 'p', &
          atom_basis(mol(dmi)%basis_number+dmj)%contr,' 1.0'
        else if (atom_basis(mol(dmi)%basis_number+dmj)%L == 2) then
          write(channel, '(A2,I3,A)') 'd', &
          atom_basis(mol(dmi)%basis_number+dmj)%contr,' 1.0'
        else if (atom_basis(mol(dmi)%basis_number+dmj)%L == 3) then
          write(channel, '(A2,I3,A)') 'f', &
          atom_basis(mol(dmi)%basis_number+dmj)%contr,' 1.0'
        else if (atom_basis(mol(dmi)%basis_number+dmj)%L == 4) then
          write(channel, '(A2,I3,A)') 'g', &
          atom_basis(mol(dmi)%basis_number+dmj)%contr,' 1.0'
        end if
        do dmk = 1, atom_basis(mol(dmi)%basis_number+dmj)%contr
          write(channel, '(F20.10,F20.10)') &
          atom_basis(mol(dmi)%basis_number+dmj)%expo(dmk), &
          atom_basis(mol(dmi)%basis_number+dmj)%coe(dmk)
        end do
      end do
      write(channel, *)
    end do
    write(channel, '(A)') '[5D]'
    write(channel, '(A)') '[7F]'
    write(channel, '(A)') '[9G]'
    ! MO coefficient
    write(channel, '(A)') '[MO]'
    do dmi = 1, fbdm
      write(channel, '(A5,E23.14)') ' Ene=', orbE(dmi)
      write(channel, '(A)') ' Spin= Alpha'
      call Calc_S2HForb(dmi)
      if (dmi <= electron_count) then
        write(channel, '(A7,F12.6)') ' Occup=', 0.5+Szorb
      else
        write(channel, '(A)') ' Occup=    0.000000'
      end if
      do dmj = 1, sbdm
        write(channel, '(I4,F20.12)') dmj, real(AO2MO(dmj, dmi))
      end do
    end do
    do dmi = 1, fbdm
      write(channel, '(A5,E23.14)') ' Ene=', orbE(dmi)
      write(channel, '(A)') ' Spin= Beta'
      call Calc_S2HForb(dmi)
      if (dmi <= electron_count) then
        write(channel, '(A7,F12.6)') ' Occup=', 0.5-Szorb
      else
        write(channel, '(A)') ' Occup=    0.000000'
      end if
      do dmj = 1, sbdm
        write(channel, '(I4,F20.12)') dmj, real(AO2MO(sbdm+dmj, dmi))
      end do
    end do
    close(channel)

    ! MOLDEN file contains the real part2 of AO2MO
    open(newunit=channel, file=trim(dir)//'/'//trim(jobname)//&
    '-realpart2.molden.input', status='replace', action='write', iostat=ios)
    if (ios /= 0) call terminate('creat .molden.input failed')
    write(channel, '(A)') '[Molden Format]'
    write(channel, '(A)') '[Title]'
    write(channel, '(A)') &
    'generated by TRESC, real part2 of MOs of job '//trim(address_job)
    write(channel, *)
    ! mol geometry
    write(channel, '(A)') '[Atoms] AU'
    do dmi = 1, atom_count
      write(channel, '(A,I4,I4,F16.10,F16.10,F16.10)') &
      element_list(mol(dmi)%atom_number), dmi, &
      mol(dmi)%atom_number, mol(dmi)%pos(1), &
      mol(dmi)%pos(2), mol(dmi)%pos(3)
    end do
    ! basis of each atom
    write(channel, '(A)') '[GTO]'
    do dmi = 1, atom_count
      write(channel, '(I3,I2)') dmi, 0
      do dmj = 0, shell_in_element(mol(dmi)%atom_number)-1
        if (atom_basis(mol(dmi)%basis_number+dmj)%L == 0) then
          write(channel, '(A2,I3,A)') 's', &
          atom_basis(mol(dmi)%basis_number+dmj)%contr,' 1.0'
        else if (atom_basis(mol(dmi)%basis_number+dmj)%L == 1) then
          write(channel, '(A2,I3,A)') 'p', &
          atom_basis(mol(dmi)%basis_number+dmj)%contr,' 1.0'
        else if (atom_basis(mol(dmi)%basis_number+dmj)%L == 2) then
          write(channel, '(A2,I3,A)') 'd', &
          atom_basis(mol(dmi)%basis_number+dmj)%contr,' 1.0'
        else if (atom_basis(mol(dmi)%basis_number+dmj)%L == 3) then
          write(channel, '(A2,I3,A)') 'f', &
          atom_basis(mol(dmi)%basis_number+dmj)%contr,' 1.0'
        else if (atom_basis(mol(dmi)%basis_number+dmj)%L == 4) then
          write(channel, '(A2,I3,A)') 'g', &
          atom_basis(mol(dmi)%basis_number+dmj)%contr,' 1.0'
        end if
        do dmk = 1, atom_basis(mol(dmi)%basis_number+dmj)%contr
          write(channel, '(F20.10,F20.10)') &
          atom_basis(mol(dmi)%basis_number+dmj)%expo(dmk), &
          atom_basis(mol(dmi)%basis_number+dmj)%coe(dmk)
        end do
      end do
      write(channel, *)
    end do
    write(channel, '(A)') '[5D]'
    write(channel, '(A)') '[7F]'
    write(channel, '(A)') '[9G]'
    ! MO coefficient
    write(channel, '(A)') '[MO]'
    do dmi = fbdm+1, 2*fbdm
      write(channel, '(A5,E23.14)') ' Ene=', orbE(dmi)
      write(channel, '(A)') ' Spin= Alpha'
      call Calc_S2HForb(dmi)
      if (dmi <= electron_count) then
        write(channel, '(A7,F12.6)') ' Occup=', 0.5+Szorb
      else
        write(channel, '(A)') ' Occup=    0.000000'
      end if
      do dmj = 1, sbdm
        write(channel, '(I4,F20.12)') dmj, real(AO2MO(dmj, dmi))
      end do
    end do
    do dmi = fbdm+1, 2*fbdm
      write(channel, '(A5,E23.14)') ' Ene=', orbE(dmi)
      write(channel, '(A)') ' Spin= Beta'
      call Calc_S2HForb(dmi)
      if (dmi <= electron_count) then
        write(channel, '(A7,F12.6)') ' Occup=', 0.5-Szorb
      else
        write(channel, '(A)') ' Occup=    0.000000'
      end if
      do dmj = 1, sbdm
        write(channel, '(I4,F20.12)') dmj, real(AO2MO(sbdm+dmj, dmi))
      end do
    end do
    close(channel)
    
    if (pVp1e) then
      ! MOLDEN file contains the imaginary part1 of AO2MO
      open(newunit=channel, file=trim(dir)//'/'//trim(jobname)//&
      '-imgpart1.molden.input', status='replace', action='write', iostat=ios)
      if (ios /= 0) call terminate('creat .molden.input failed')
      write(channel, '(A)') '[Molden Format]'
      write(channel, '(A)') '[Title]'
      write(channel, '(A)') &
      'generated by TRESC, imaginary part1 of MOs of job '//trim(address_job)
      write(channel, *)
      ! mol geometry
      write(channel, '(A)') '[Atoms] AU'
      do dmi = 1, atom_count
        write(channel, '(A,I4,I4,F16.10,F16.10,F16.10)') &
        element_list(mol(dmi)%atom_number), dmi, &
        mol(dmi)%atom_number, mol(dmi)%pos(1), &
        mol(dmi)%pos(2), mol(dmi)%pos(3)
      end do
      ! basis of each atom
      write(channel, '(A)') '[GTO]'
      do dmi = 1, atom_count
        write(channel, '(I3,I2)') dmi, 0
        do dmj = 0, shell_in_element(mol(dmi)%atom_number)-1
          if (atom_basis(mol(dmi)%basis_number+dmj)%L == 0) then
            write(channel, '(A2,I3,A)') 's', &
            atom_basis(mol(dmi)%basis_number+dmj)%contr,' 1.0'
          else if (atom_basis(mol(dmi)%basis_number+dmj)%L == 1) then
            write(channel, '(A2,I3,A)') 'p', &
            atom_basis(mol(dmi)%basis_number+dmj)%contr,' 1.0'
          else if (atom_basis(mol(dmi)%basis_number+dmj)%L == 2) then
            write(channel, '(A2,I3,A)') 'd', &
            atom_basis(mol(dmi)%basis_number+dmj)%contr,' 1.0'
          else if (atom_basis(mol(dmi)%basis_number+dmj)%L == 3) then
            write(channel, '(A2,I3,A)') 'f', &
            atom_basis(mol(dmi)%basis_number+dmj)%contr,' 1.0'
          else if (atom_basis(mol(dmi)%basis_number+dmj)%L == 4) then
            write(channel, '(A2,I3,A)') 'g', &
            atom_basis(mol(dmi)%basis_number+dmj)%contr,' 1.0'
          end if
          do dmk = 1, atom_basis(mol(dmi)%basis_number+dmj)%contr
            write(channel, '(F20.10,F20.10)') &
            atom_basis(mol(dmi)%basis_number+dmj)%expo(dmk), &
            atom_basis(mol(dmi)%basis_number+dmj)%coe(dmk)
          end do
        end do
        write(channel, *)
      end do
      write(channel, '(A)') '[5D]'
      write(channel, '(A)') '[7F]'
      write(channel, '(A)') '[9G]'
      ! MO coefficient
      write(channel, '(A)') '[MO]'
      do dmi = 1, fbdm
        write(channel, '(A5,E23.14)') ' Ene=', orbE(dmi)
        write(channel, '(A)') ' Spin= Alpha'
        call Calc_S2HForb(dmi)
        if (dmi <= electron_count) then
          write(channel, '(A7,F12.6)') ' Occup=', 0.5+Szorb
        else
          write(channel, '(A)') ' Occup=    0.000000'
        end if
        do dmj = 1, sbdm
          write(channel, '(I4,F20.12)') dmj, aimag(AO2MO(dmj, dmi))
        end do
      end do
      do dmi = 1, fbdm
        write(channel, '(A5,E23.14)') ' Ene=', orbE(dmi)
        write(channel, '(A)') ' Spin= Beta'
        call Calc_S2HForb(dmi)
        if (dmi <= electron_count) then
          write(channel, '(A7,F12.6)') ' Occup=', 0.5-Szorb
        else
          write(channel, '(A)') ' Occup=    0.000000'
        end if
        do dmj = 1, sbdm
          write(channel, '(I4,F20.12)') dmj, aimag(AO2MO(sbdm+dmj, dmi))
        end do
      end do
      close(channel)

      ! MOLDEN file contains the imaginary part2 of AO2MO
      open(newunit=channel, file=trim(dir)//'/'//trim(jobname)//&
      '-imgpart2.molden.input', status='replace', action='write', iostat=ios)
      if (ios /= 0) call terminate('creat .molden.input failed')
      write(channel, '(A)') '[Molden Format]'
      write(channel, '(A)') '[Title]'
      write(channel, '(A)') &
      'generated by TRESC, imaginary part2 of MOs of job '//trim(address_job)
      write(channel, *)
      ! mol geometry
      write(channel, '(A)') '[Atoms] AU'
      do dmi = 1, atom_count
        write(channel, '(A,I4,I4,F16.10,F16.10,F16.10)') &
        element_list(mol(dmi)%atom_number), dmi, &
        mol(dmi)%atom_number, mol(dmi)%pos(1), &
        mol(dmi)%pos(2), mol(dmi)%pos(3)
      end do
      ! basis of each atom
      write(channel, '(A)') '[GTO]'
      do dmi = 1, atom_count
        write(channel, '(I3,I2)') dmi, 0
        do dmj = 0, shell_in_element(mol(dmi)%atom_number)-1
          if (atom_basis(mol(dmi)%basis_number+dmj)%L == 0) then
            write(channel, '(A2,I3,A)') 's', &
            atom_basis(mol(dmi)%basis_number+dmj)%contr,' 1.0'
          else if (atom_basis(mol(dmi)%basis_number+dmj)%L == 1) then
            write(channel, '(A2,I3,A)') 'p', &
            atom_basis(mol(dmi)%basis_number+dmj)%contr,' 1.0'
          else if (atom_basis(mol(dmi)%basis_number+dmj)%L == 2) then
            write(channel, '(A2,I3,A)') 'd', &
            atom_basis(mol(dmi)%basis_number+dmj)%contr,' 1.0'
          else if (atom_basis(mol(dmi)%basis_number+dmj)%L == 3) then
            write(channel, '(A2,I3,A)') 'f', &
            atom_basis(mol(dmi)%basis_number+dmj)%contr,' 1.0'
          else if (atom_basis(mol(dmi)%basis_number+dmj)%L == 4) then
            write(channel, '(A2,I3,A)') 'g', &
            atom_basis(mol(dmi)%basis_number+dmj)%contr,' 1.0'
          end if
          do dmk = 1, atom_basis(mol(dmi)%basis_number+dmj)%contr
            write(channel, '(F20.10,F20.10)') &
            atom_basis(mol(dmi)%basis_number+dmj)%expo(dmk), &
            atom_basis(mol(dmi)%basis_number+dmj)%coe(dmk)
          end do
        end do
        write(channel, *)
      end do
      write(channel, '(A)') '[5D]'
      write(channel, '(A)') '[7F]'
      write(channel, '(A)') '[9G]'
      ! MO coefficient
      write(channel, '(A)') '[MO]'
      do dmi = fbdm+1, 2*fbdm
        write(channel, '(A5,E23.14)') ' Ene=', orbE(dmi)
        write(channel, '(A)') ' Spin= Alpha'
        call Calc_S2HForb(dmi)
        if (dmi <= electron_count) then
          write(channel, '(A7,F12.6)') ' Occup=', 0.5+Szorb
        else
          write(channel, '(A)') ' Occup=    0.000000'
        end if
        do dmj = 1, sbdm
          write(channel, '(I4,F20.12)') dmj, aimag(AO2MO(dmj, dmi))
        end do
      end do
      do dmi = fbdm+1, 2*fbdm
        write(channel, '(A5,E23.14)') ' Ene=', orbE(dmi)
        write(channel, '(A)') ' Spin= Beta'
        call Calc_S2HForb(dmi)
        if (dmi <= electron_count) then
          write(channel, '(A7,F12.6)') ' Occup=', 0.5-Szorb
        else
          write(channel, '(A)') ' Occup=    0.000000'
        end if
        do dmj = 1, sbdm
          write(channel, '(I4,F20.12)') dmj, aimag(AO2MO(sbdm+dmj, dmi))
        end do
      end do
      close(channel)
    end if
  end subroutine Dump_MOLDEN
  
end module SCF