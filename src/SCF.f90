!> @file SCF.f90
!!
!! @brief single-configuration self-consistent field calculation
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
  
  ! density matrices and mol orbital coefficients
  integer                 :: Nalpha, Nbeta  ! number of alpha and beta elctron
  ! spinor MO coefficients, in order of (AO1,0), (0,AO1), ... (AOn,0), (0,AOn)
  !DIR$ ATTRIBUTES ALIGN:align_size :: AO2MO, AO2MOalpha, AO2MObeta
  !DIR$ ATTRIBUTES ALIGN:align_size :: rho_m, rotation, AOsupp, Fock
  complex(dp),allocatable :: AO2MO(:,:)
  real(dp),allocatable    :: Gaualpha(:)    ! Gaussian alpha orbital coefficient
  real(dp),allocatable    :: AO2MOalpha(:,:)! Gaussian alpha orbital coefficient
  real(dp),allocatable    :: Gaubeta(:)     ! Gaussian beta orbital coefficient
  real(dp),allocatable    :: AO2MObeta(:,:) ! Gaussian beta orbital coefficient
  complex(dp),allocatable :: rho_m(:,:)     ! density matrix, complex Hermitian
  complex(dp),allocatable :: rotation(:,:)  ! rotate one orbital with another
  real(dp)                :: RMSDP, maxDP
  real(dp),allocatable    :: AOsupp(:,:)
  integer                 :: evl_count_f    ! number of eigenvalues in Fock diag
  complex(dp),allocatable :: Fock(:,:)      ! Fock matrix
  integer,allocatable     :: isupp_ev_f(:)
  integer                 :: fliwork
  integer                 :: flwork
  integer                 :: lrwork
  complex(dp),allocatable :: fwork(:)       ! work of Fock for zheevr input
  real(dp),allocatable    :: rwork(:)       ! rwork of Fock for zheevr input
  integer,allocatable     :: fiwork(:)      ! iwork of Fock for zheevr input
  integer                 :: finfo          ! info of calling lapack functions
  
  logical                 :: ini_rou =.true.! initial density matrix loaded
  
  !--------------<one electron Fock>--------------
  !DIR$ ATTRIBUTES ALIGN:align_size :: Fock1, oper1, oper2, oper3
  !DIR$ ATTRIBUTES ALIGN:align_size :: oper4, oper5, oper6, Ve
  !DIR$ ATTRIBUTES ALIGN:align_size :: pxVepx, pyVepy, pzVepz, pxVepy
  !DIR$ ATTRIBUTES ALIGN:align_size :: pyVepx, pxVepz, pzVepx, pyVepz
  !DIR$ ATTRIBUTES ALIGN:align_size :: pzVepy, Ap, ApRp, SRp, ARVRA
  !DIR$ ATTRIBUTES ALIGN:align_size :: ARVeRA, AVA, AVeA, exAO2p2
  complex(dp),allocatable :: Fock1(:,:)     ! one electron Fock matrix
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
  
  !--------------<two electron Fock>--------------
  !DIR$ ATTRIBUTES ALIGN:align_size :: Fock2HFcol, Fock2HFexc, Fock2KSexc
  !DIR$ ATTRIBUTES ALIGN:align_size :: Fock2KScor, swint
  complex(dp),allocatable :: Fock2HFcol(:,:)! HF Coulomb matrix
  complex(dp),allocatable :: Fock2HFexc(:,:)! HF Exchange matrix
  complex(dp),allocatable :: Fock2KSexc(:,:)! KS Exchange matrix
  complex(dp),allocatable :: Fock2KScor(:,:)! KS correlation matrix
  real(dp),allocatable    :: swint(:,:)! <ij||ij> as well as <kl||kl>
  logical                 :: ndschwarz = .true.
  
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
  real(dp)                :: Ecore          ! one selectron (core) energy
  real(dp)                :: T              ! kinetic energy
  real(dp)                :: V              ! electron-nuclear attraction energy
  real(dp)                :: ESOC           ! SOC energy
  real(dp)                :: ESR            ! SRTP energy
  real(dp)                :: emd4           ! dispersion energy calc by DFT-D4

  !DIR$ ATTRIBUTES ALIGN:align_size :: rho_pre, rho_pre_pre, rho_history, Rsd
  ! damping & DIIS(AX=B)
  ! privious rho_m of current iteration
  complex(dp),allocatable :: rho_pre(:,:)
  ! privious rho_m of privious iteration
  complex(dp),allocatable :: rho_pre_pre(:,:)
  complex(dp),allocatable :: rho_history(:,:,:)! coeff of subsp iteration
  complex(dp),allocatable :: Rsd(:,:,:)        ! residuals of rho_m
  complex(dp),allocatable :: DIISmat(:,:)      ! A
  complex(dp),allocatable :: DIIuork(:)
  integer                 :: lDIIuork
  integer,allocatable     :: ipiv(:)
  integer                 :: DIIsinfo
  real(dp)                :: damp_coe          ! damp coeff of direct/DIIS SCF
  real(dp)                :: dE_pre
  logical                 :: forward = .true.
  
  contains

!------------------------------------------------------------
!> HF/KS-SCF procedure for DKH0/DKH2 Hamiltonian
  subroutine DKH_SCF()
    implicit none
    character(len = 40) :: keyword
    write(60,'(A)') 'Module SCF:'
    write(60,'(A)') '  construct one electron Fock matrix'
    call Fock1e()
    write(60,'(A)') '  complete! stored in Fock1'
    write(60,'(A)') '  calculate nuclear repulsion energy'
    nucE = 0.0_dp
    do loop_i = 1, atom_count
      do loop_j = loop_i + 1, atom_count
        nucE = nucE + (real(mol(loop_i) % atom_number) * &
        real(mol(loop_j) % atom_number))/dsqrt((mol(loop_i) % &
        pos(1) - mol(loop_j) % pos(1))**2 + &
        (mol(loop_i) % pos(2) - mol(loop_j) % &
        pos(2))**2 + (mol(loop_i) % pos(3) - &
        mol(loop_j) % pos(3))**2)
      end do
    end do
    if (nucE >= 1E12) call terminate(&
    'nuclear repulsive energy anomaly, may due to overlap atomic coordinates')
    write(60,'(A,F12.7,A)') &
    '  complete! nuclear repulsive energy = ', nucE, ' Eh'
    write(60,'(A)') '  initialize the functional'
    if (fx_id /= -1) call Fockxc_init()
    write(60,'(A)') '  complete!'
    write(60,'(A)') '  SCF settings:'
    write(60,'(A,I4)') '  -- maxiter =',maxiter
    write(60,'(A,E10.3)') '  -- conv_tol = ',conver_tol
    write(60,'(A,F6.3)') '  -- damp =',damp
    write(60,'(A,E10.3)') '  -- cutdamp = ',cutdamp
    write(60,'(A,I4)') '  -- nodiis =',nodiis
    write(60,'(A,I3)') '  -- subsp =',subsp
    write(60,'(A,F6.3)') '  -- diisdamp =',diisdamp
    write(60,'(A,E10.3)') '  -- cutdiis = ',cutdiis
    write(60,'(A)') '  ----------<SCF>----------'
    allocate(Fock(2*fbdm,2*fbdm))
    allocate(orbE(2*fbdm))
    allocate(oper6(2*fbdm,2*fbdm))
    allocate(isupp_ev_f(4*fbdm))
    allocate(oper3(2*fbdm,2*fbdm))
    allocate(oper4(2*fbdm,2*fbdm))
    
    ! DIIS
    ! AO2MO(new) = ��(i,subsp) DIIScoe(i)*(rho_history(i)+diisdamp*Rsd(i))
    allocate(Rsd(subsp,2*sbdm,2*sbdm))
    allocate(DIISmat(subsp+1,subsp+1))
    allocate(ipiv(subsp+1))
    allocate(rho_history(subsp, 2*sbdm, 2*sbdm))
    allocate(rho_pre(2*sbdm, 2*sbdm))
    allocate(rho_pre_pre(2*sbdm, 2*sbdm))
    do loop_i = 1, maxiter
      if (loop_i /= 1) then
        write(60,*)
        write(60,*)
        write(60,'(A,I3)') '  SCF iter ',loop_i
        open(12, file=address_molecule, status="old", action="read")
        do
          read(12,*,iostat = ios) keyword
          if (ios /= 0) exit
          call lowercase(keyword)
          if (index(keyword,'threads') == 1) then
            if (index(keyword,'=') /= 0) then
              read(keyword(index(keyword,'=') + 1 : len(trim(keyword))),&
              "(I2)",iostat = ios) threads_new
              if (threads_new /= threads) then
                threads = threads_new
                write(60,'(A,I2)') '  number of threads change to ',threads
              end if
            end if
          else if (index(keyword,'damp') == 1) then
            if (index(keyword,'=') /= 0) then
              read(keyword(index(keyword,'=') + 1 : len(trim(keyword))),&
              "(F20.12)",iostat = ios) damp_new
              if (damp_new /= damp) then
                damp = damp_new
                write(60,'(A,F20.12)') '  dynamical damp change to ',damp
              end if
            end if
          else if (index(keyword,'diisdamp') == 1) then
            if (index(keyword,'=') /= 0) then
              read(keyword(index(keyword,'=') + 1 : len(trim(keyword))),&
              "(F20.12)",iostat = ios) diisdamp_new
              if (diisdamp_new /= diisdamp) then
                diisdamp = diisdamp_new
                write(60,'(A,F20.12)') &
                '  dynamical damp in DIIS change to ',diisdamp
              end if
            end if
          else if (index(keyword,'cutdiis') == 1) then
            if (index(keyword,'=') /= 0) then
              read(keyword(index(keyword,'=') + 1 : len(trim(keyword))),&
              "(F20.12)",iostat = ios) cutdiis_new
              if (cutdiis_new /= cutdiis) then
                cutdiis = cutdiis_new
                write(60,'(A,F20.12)') '  DIIS threshold change to ',cutdiis
              end if
            end if
          else if (index(keyword,'cutdamp') == 1) then
            if (index(keyword,'=') /= 0) then
              read(keyword(index(keyword,'=') + 1 : len(trim(keyword))),&
              "(F20.12)",iostat = ios) cutdamp_new
              if (cutdamp_new /= cutdamp) then
                cutdamp = cutdamp_new
                write(60,'(A,F20.12)') '  damp threshold change to ',cutdamp
              end if
            end if
          end if
        end do
        close(12)
      else
        write(60,'(A)') '  read density matrix'
        call assign_rou()
        write(60,'(A)') '  complete! stored in rho_m'
      end if

      write(60,'(A)') '  construct two electron Fock matrix'
      call Fock2e()
      if (fx_id /= -1) then
        if (fc_id /= -1) then
          Fock = Fock1 + Fock2HFcol + x_HF*Fock2HFexc + &
          (1.0_dp-x_HF)*Fock2KSexc + Fock2KScor
          write(60,'(A)') '  complete! stored in Fock2HFcol, Fock2HFexc'
          write(60,'(A)') '  Fock2KSexc, Fock2KScor'
        else
          x_HF = 0.0_dp
          KScor = 0.0_dp
          Fock = Fock1 + Fock2HFcol + Fock2KSexc
          write(60,'(A)') '  complete! stored in Fock2HFcol, Fock2HFexc'
          write(60,'(A)') '  Fock2KSexc'
        end if
      else
        ! pure Hartree-Fock
        x_HF = 1.0_dp
        KScor = 0.0_dp
        KSexc = 0.0_dp
        Fock = Fock1 + Fock2HFcol + Fock2HFexc
        write(60,'(A)') '  complete! stored in Fock2HFcol, Fock2HFexc'
      end if

      write(60,'(A)') '  diagonalization of Fock matrix'
      allocate(fwork(1))
      allocate(fiwork(1))
      allocate(rwork(1))
      call zheevr(&
      'V','A','U',2*fbdm,Fock,2*fbdm,0.0_dp,0.0_dp,0,0,safmin,evl_count_f,orbE,&
      oper3,2*fbdm,isupp_ev_f,fwork,-1,rwork,-1,fiwork,-1,finfo)
      flwork = nint(real(fwork(1)))
      fliwork = fiwork(1)
      lrwork = nint(rwork(1))
      deallocate(fwork)
      deallocate(fiwork)
      deallocate(rwork)
      allocate(fwork(flwork))
      allocate(fiwork(fliwork))
      allocate(rwork(lrwork))
      call zheevr(&
      'V','A','U',2*fbdm,Fock,2*fbdm,0.0_dp,0.0_dp,0,0,safmin,evl_count_f,orbE,&
      oper3,2*fbdm,isupp_ev_f,fwork,flwork,rwork,lrwork,fiwork,fliwork,finfo)
      if (finfo < 0) call terminate(&
      'Fock matrix diagonalization failure, illegal input of zheevr')
      if (finfo > 0) call terminate(&
      'Fock matrix diagonalization failure, internal error of zheevr')
      if (evl_count_f < 2*fbdm) then
        call terminate('number of MO less than 2*fbdm')
      else
        write(60,'(A,I5,A)') '  complete!',evl_count_f,' eigenvectors found'
      end if
      ! frontier orbital energy
      write(60,'(A)') '  frontier orbital energy (A.U.)'
      call calc_S2HForb(electron_count)
      write(60,'(A,I3,F12.6,A,F6.3)') &
      '  -- HOMO ', electron_count, orbE(electron_count), ' <Sz> = ',Szorb
      call calc_S2HForb(electron_count+1)
      write(60,'(A,I3,F12.6,A,F6.3)') &
      '  -- LUMO ', electron_count+1, orbE(electron_count+1), ' <Sz> = ',Szorb
      write(60,'(A,F12.6)') &
      '  -- gap  ', orbE(electron_count+1) - orbE(electron_count)
      deallocate(fwork)
      deallocate(fiwork)
      deallocate(rwork)
      ! energy components calculation
      write(60,'(A)') '  calculate energy components (A.U.)'
      ! HF Coulomb energy
      HFCol = 0.0_dp
      call zgemm('C', 'N', 2*fbdm, 2*fbdm, 2*fbdm, c1, oper3, 2*fbdm, &
      Fock2HFcol, 2*fbdm, c0, oper6, 2*fbdm)
      call zgemm( 'N', 'N', 2*fbdm, 2*fbdm, 2*fbdm, c1, oper6, 2*fbdm, &
      oper3, 2*fbdm, c0, oper4, 2*fbdm)
      do loop_j = 1, electron_count
        HFCol = HFCol + real(oper4(loop_j,loop_j))
      end do
      write(60,'(A,F12.6)') '  -- HF Coulomb energy                    ', HFCol

      ! HF exchange energy
      HFexc = 0.0_dp
      call zgemm('C', 'N', 2*fbdm, 2*fbdm, 2*fbdm, c1, oper3, 2*fbdm, &
      Fock2HFexc, 2*fbdm, c0, oper6, 2*fbdm)
      call zgemm( 'N', 'N', 2*fbdm, 2*fbdm, 2*fbdm, c1, oper6, 2*fbdm, &
      oper3, 2*fbdm, c0, oper4, 2*fbdm)
      do loop_j = 1, electron_count
        HFexc = HFexc + real(oper4(loop_j,loop_j))
      end do
      write(60,'(A,F12.6)') &
      '  -- HF exchange energy                   ', x_HF*HFexc
      write(60,'(A,F12.6)') &
      '  -- KS exchange energy                   ', (1-x_HF)*KSexc
      write(60,'(A,F12.6)') '  -- KS correlation energy                ', KScor

      ! core energy
      Ecore = 0.0_dp
      call zgemm('C', 'N', 2*fbdm, 2*fbdm, 2*fbdm, c1, oper3, 2*fbdm, &
      Fock1, 2*fbdm, c0, oper6, 2*fbdm)
      call zgemm( 'N', 'N', 2*fbdm, 2*fbdm, 2*fbdm, c1, oper6, 2*fbdm, &
      oper3, 2*fbdm, c0, oper4, 2*fbdm)
      do loop_j = 1, electron_count
        Ecore = Ecore + real(oper4(loop_j,loop_j))
      end do
      write(60,'(A,F12.6)') '  -- core energy                          ', Ecore

      ! kinetic energy
      T = 0.0_dp
      call zgemm( 'C', 'N', 2*fbdm, 2*fbdm, 2*fbdm, c1, oper3, 2*fbdm, &
      exi_T_j, 2*fbdm, c0, oper6, 2*fbdm)
      call zgemm( 'N', 'N', 2*fbdm, 2*fbdm, 2*fbdm, c1, oper6, 2*fbdm, &
      oper3, 2*fbdm, c0, oper4, 2*fbdm)
      do loop_j = 1, electron_count
        T = T + real(oper4(loop_j,loop_j))
      end do
      write(60,'(A,F12.6)') '  -- electron kinetic energy              ', T

      ! electron-nuclear attraction energy
      V = 0.0_dp
      call zgemm( 'C', 'N', 2*fbdm, 2*fbdm, 2*fbdm, c1, oper3, 2*fbdm, &
      exi_V_j, 2*fbdm, c0, oper6, 2*fbdm)
      call zgemm( 'N', 'N', 2*fbdm, 2*fbdm, 2*fbdm, c1, oper6, 2*fbdm, &
      oper3, 2*fbdm, c0, oper4, 2*fbdm)
      do loop_j = 1, electron_count
        V = V + real(oper4(loop_j,loop_j))
      end do
      write(60,'(A,F12.6)') '  -- electron-nuclear attraction energy   ', V

      if (DKH_order == 2) then
        ESOC = 0.0_dp
        call zgemm( 'C', 'N', 2*fbdm, 2*fbdm, 2*fbdm, c1, oper3, 2*fbdm, &
        exSOC, 2*fbdm, c0, oper6, 2*fbdm)
        call zgemm( 'N', 'N', 2*fbdm, 2*fbdm, 2*fbdm, c1, oper6, 2*fbdm, &
        oper3, 2*fbdm, c0, oper4, 2*fbdm)
        do loop_j = 1, electron_count
          ESOC = ESOC + real(oper4(loop_j,loop_j))
        end do
        write(60,'(A,F12.6)') '  -- spin-orbital coupling energy         ', ESOC
        if (SRTP_type) then
          ESR = 0.0_dp
          call zgemm( 'C', 'N', 2*fbdm, 2*fbdm, 2*fbdm, c1, oper3, 2*fbdm, &
          exSR, 2*fbdm, c0, oper6, 2*fbdm)
          call zgemm( 'N', 'N', 2*fbdm, 2*fbdm, 2*fbdm, c1, oper6, 2*fbdm, &
          oper3, 2*fbdm, c0, oper4, 2*fbdm)
          do loop_j = 1, electron_count
            ESR = ESR + real(oper4(loop_j,loop_j))
          end do
          write(60,'(A,F12.6)') &
          '  -- SRTP & radiative correction energy   ', ESR
        end if
      end if

      ! electronic energy
      if (loop_i /= 1) molE_pre = molE
      molE = nucE + Ecore + 0.5_dp*(HFCol+x_HF*HFexc)+(1.0_dp-x_HF)*KSexc +KScor
      
      ! Virial ratio
      if (DKH_order == 0) then
        Virial = -(nucE+V+0.5_dp*(HFCol+x_HF*HFexc)+&
        (1.0_dp-x_HF)*KSexc +KScor)/T
      else if (DKH_order == 2) then
        if (SRTP_type) then
          Virial = -(nucE+Ecore-T-Esoc-ESR+0.5_dp*(HFCol+x_HF*HFexc)+&
          (1.0_dp-x_HF)*KSexc +KScor)/T
        else
          Virial = -(nucE+Ecore-T-Esoc+0.5_dp*(HFCol+x_HF*HFexc)+&
          (1.0_dp-x_HF)*KSexc +KScor)/T
        end if
      end if
      write(60,'(A,F12.6)') '  -- -<V>/<T>                             ', Virial
      ! de-orthogonalization
      call zgemm( 'N', 'N', 2*sbdm, 2*fbdm, 2*fbdm, c1, exXm, 2*sbdm, &
      oper3, 2*fbdm, c0, AO2MO, 2*sbdm)
      write(60,'(A)') '  AO2MO dump to .ao2mo file'
      call dump_matrix('ao2mo', AO2MO, 2*sbdm, 2*fbdm)
      ! convergence check
      if (loop_i == 1) then
        write(60,'(A,F12.6)') '  SCF energy (A.U.) = ',molE
      else
        write(60,'(A,F12.6,A,E10.3)') '  SCF energy (A.U.) = ',molE,&
        '; delta E = ',molE - molE_pre
        if (abs(molE - molE_pre) < conver_tol .and. abs(RMSDP) < 5*conver_tol &
        .and. forward .and. (damp_coe < 0.01 .or. &
        damp_coe < 0.1*(log10(conver_tol)-log10(abs(molE - molE_pre))))) then
          write(60,'(A)') '  convergence tolerance met, SCF done!'
          exit
        else
          write(60,'(A)') '  convergence tolerance not met'
        end if
      end if
      write(60,'(A)') '  construct density matrix'
      call assign_rou()
      write(60,'(A)') '  complete! stored in rho_m'
      write(60,'(A)') '  DIIS information'
      ! generate next rho_m by DIIS method
      if (loop_i <= nodiis) then
        
        !--------------<damping>-----------------
        if (loop_i <= nodiis - subsp) then
          if (loop_i > 2 .and. abs(molE - molE_pre) >= 1.5*abs(dE_pre)) then
            damp_coe = damp_coe + (1.0_dp-damp_coe)*0.5_dp
            if (damp_coe < 0.91) then
              rho_m = (1.0_dp - damp_coe) * rho_pre + damp_coe * rho_pre_pre
              write(60,'(A)') '  -- fallback '
              forward = .false.
              cycle
            end if
          end if
        end if
        if (abs(molE - molE_pre) >= cutdamp) then
          damp_coe = damp
          rho_m = (1.0 - damp_coe) * rho_m + damp_coe * rho_pre
          write(60,'(A,F6.3)') '  -- damped ', damp_coe
        else if (abs(molE - molE_pre) >= cutdamp/100.0) then
          damp_coe = 0.5_dp*damp*log10(abs(molE - molE_pre)) + &
          damp*(1.0_dp-0.5_dp*log10(cutdamp))
          rho_m = (1.0 - damp_coe) * rho_m + damp_coe * rho_pre
          write(60,'(A,F6.3)') '  -- damped ', damp_coe
        else
          damp_coe = 0.0_dp
          write(60,'(A)') '  -- undamped '
        end if
        dE_pre = molE - molE_pre
        !--------------<damping>-----------------

        if (loop_i <= nodiis - subsp) then
          write(60,'(A)') '  -- no DIIS acceleration'
        else
          ! update Rsd
          do loop_j = 1, 2*sbdm
            do loop_k = 1, 2*sbdm
              Rsd(loop_i-(nodiis-subsp), loop_j, loop_k) = &
              rho_m(loop_j, loop_k) - rho_pre(loop_j, loop_k)
            end do
          end do
          ! update rho_history
          do loop_j = 1, 2*sbdm
            do loop_k = 1, 2*sbdm
              rho_history(loop_i-(nodiis-subsp), loop_j, loop_k) = &
              rho_m(loop_j, loop_k)
            end do
          end do
          write(60,'(A,I2,A,I2)') &
          '  -- DIIS subspace filling ',loop_i-(nodiis-subsp),'/',subsp
        end if
      else

        !--------------<DIIS damping>-----------------
        if (abs(molE - molE_pre) >= cutdiis) then
          damp_coe = diisdamp
          write(60,'(A,F6.3)') '  -- DIIS damped ', damp_coe
        else if (abs(molE - molE_pre) >= cutdiis/100.0) then
          damp_coe = 0.5_dp*diisdamp*log10(abs(molE - molE_pre)) + &
          diisdamp*(1.0_dp-0.5_dp*log10(cutdiis))
          write(60,'(A,F6.3)') '  -- DIIS damped ', damp_coe
        else
          damp_coe = 0.0_dp
          write(60,'(A)') '  -- DIIS undamped'
        end if
        !--------------<DIIS damping>-----------------

        ! update Rsd
        do loop_j = 2, subsp
          do loop_k = 1, 2*sbdm
            do loop_l = 1, 2*sbdm
              Rsd(loop_j - 1, loop_k, loop_l) = Rsd(loop_j, loop_k, loop_l)
            end do
          end do
        end do
        do loop_j = 1, 2*sbdm
          do loop_k = 1, 2*sbdm
            Rsd(subsp, loop_j, loop_k) = rho_m(loop_j, loop_k) - &
            rho_pre(loop_j, loop_k)
          end do
        end do
        ! update rho_history
        do loop_j = 2, subsp
          do loop_k = 1, 2*sbdm
            do loop_l = 1, 2*sbdm
              rho_history(loop_j - 1, loop_k, loop_l) = &
              rho_history(loop_j, loop_k, loop_l)
            end do
          end do
        end do
        do loop_j = 1, 2*sbdm
          do loop_k = 1, 2*sbdm
            rho_history(subsp, loop_j, loop_k) = rho_m(loop_j, loop_k)
          end do
        end do
        ! construct DIISmat
        DIISmat = c0
        do loop_j = 1, subsp
          do loop_k = 1, subsp
            do loop_l = 1, 2*sbdm
              do loop_m = 1, 2*sbdm
                DIISmat(loop_j, loop_k) = DIISmat(loop_j, loop_k) + &
                conjg(Rsd(loop_j, loop_l, loop_m))*Rsd(loop_k, loop_l, loop_m)
              end do
            end do
          end do
        end do
        do loop_j = 1, subsp
          DIISmat(subsp+1, loop_j) = c1
          DIISmat(loop_j, subsp+1) = c1
        end do
        ! solveg residual equation
        ! dgesv and dspsv will cause V_integral_2e conflict for unknown reason
        ! since DIISmat (and its inverse) is real symmetric, plus the column 
        ! vector is simple, use the inverse of DIISmat to solve directly
        call zgetrf(subsp+1, subsp+1, DIISmat, subsp+1, ipiv, DIISinfo)
        if (DIISinfo < 0) then
          call terminate('DIIS solution failure, illegal input of dgetrf')
        else if (DIISinfo > 0) then
          call terminate('DIIS solution failure, internal error of dgetrf')
        end if
        allocate(DIIuork(1))
        call zgetri(subsp+1, DIISmat, subsp+1, ipiv, DIIuork, -1, DIISinfo)
        lDIIuork = nint(real(DIIuork(1)))
        deallocate(DIIuork)
        allocate(DIIuork(lDIIuork))
        call zgetri(subsp+1,DIISmat,subsp+1,ipiv,DIIuork,lDIIuork,DIISinfo)
        deallocate(DIIuork)
        if (DIISinfo < 0) then
          call terminate('DIIS solution failure, illegal input of dgetri')
        else if (DIISinfo > 0) then
          call terminate('DIIS solution failure, internal error of dgetri')
        end if
        ! generate new rho_m
        rho_m = c0
        do loop_j = 1, subsp
          do loop_k = 1, 2*sbdm
            do loop_l = 1, 2*sbdm
              rho_m(loop_k,loop_l) = rho_m(loop_k,loop_l) + &
              DIISmat(loop_j,subsp+1) * (rho_history(loop_j,loop_k,loop_l) + &
              damp_coe*Rsd(loop_j,loop_k,loop_l))
            end do
          end do
        end do
        write(60,'(A,F10.6,A,F10.6)') &
        '  -- predicted residual', -real(DIISmat(subsp+1,subsp+1)), &
        ',',-aimag(DIISmat(subsp+1,subsp+1))
        do loop_j = 1, subsp
          write(60,'(A,I2,F10.6,A,F10.6)') &
          '  -- subsp coe',loop_j, real(DIISmat(loop_j,subsp+1)), &
          ',', aimag(DIISmat(loop_j,subsp+1))
        end do
      end if
      if (forward) then
        if (loop_i > 1) rho_pre_pre = rho_pre
        rho_pre = rho_m
      end if
      forward = .true.
    end do
    if (abs(molE - molE_pre) < conver_tol .and. loop_i < maxiter) then
      if (DKH_order == 0) write(60,'(A)') '  DKH0 SCF succeed!'
      if (DKH_order == 2) write(60,'(A)') '  DKH2 SCF succeed!'
    else
      if (DKH_order == 0) write(60,'(A)') '  DKH0 SCF failed!'
      if (DKH_order == 2) write(60,'(A)') '  DKH2 SCF failed!'
    end if
    if (d4) then
      emd4 = dftd4()
      molE = molE + emd4
    end if
  end subroutine DKH_SCF

  !------------------------------------------------------------
  !> print final wave function information
  subroutine outprint()
    implicit none
    integer :: iatom, ishell
    close(80)
    write(60,*)
    write(60,*)
    write(60,'(A)') &
    '  ============================================================='
    write(60,'(A)') '                        mol INFO'
    write(60,'(A)') &
    '  ============================================================='
    write(60,'(A,F12.6)') &
    '  total electronic energy / Eh                  ...',molE
    write(60,'(A,F12.6)') &
    '  nuclear repulsive energy / Eh                 ...',nucE
    write(60,'(A,F12.6)') &
    '  HF Coulomb energy / Eh                        ...',HFCol
    write(60,'(A,F12.6)') &
    '  HF exchange energy / Eh                       ...',x_HF*HFexc
    write(60,'(A,F12.6)') &
    '  KS exchange energy / Eh                       ...',(1.0_dp-x_HF)*KSexc
    write(60,'(A,F12.6)') &
    '  KS correlation energy / Eh                    ...',KScor
    write(60,'(A,F12.6)') &
    '  core energy / Eh                              ...',Ecore
    write(60,'(A,F12.6)') &
    '  electron kinetic energy / Eh                  ...',T
    write(60,'(A,F12.6)') &
    '  electron-nuclear attraction energy / Eh       ...',V
    if (DKH_order == 2) then
      write(60,'(A,F12.6)') &
      '  spin-orbital coupling energy / Eh             ...',ESOC
      if (SRTP_type) write(60,'(A,F12.6)') &
      '  SRTP & radiative correction energy / Eh       ...',ESR
    end if
    if (d4) write(60,'(A,F12.6)') &
    '  dispersion energy (DFT-D4) / Eh               ...',emd4
    write(60,'(A,F12.6)') &
    '  Virial ratio                                  ...',Virial
    if (DKH_order /= 0) then
      write(60,*)
      write(60,'(A)') '  Note: relativistic calculation causes the Virial ratio'
      write(60,'(A)') '  to deviate (usually below) 2.0'
      write(60,*)
    end if
    write(60,'(A,F12.6)') &
    '  total alpha electron                          ...',totalpha
    write(60,'(A,F12.6)') &
    '  total beta electron                           ...',totbeta
    write(60,'(A,F12.6)') &
    '  <Sz*(Sz+1)> / hbar**2                         ...',&
    ((totalpha-totbeta)/2.0)*((totalpha-totbeta)/2.0+1.0_dp)
    write(60,'(A,F12.6)') &
    '  <S**2> / hbar**2                              ...',S__2
    if (DKH_order == 2) then
      write(60,*)
      write(60,'(A)') &
      '  Note: <S**2> may be contaminated by electron correlation,'
      write(60,'(A)') &
      '  discussion of SOC suggests <S**2> data of 1e orbitals.'
    end if
    write(60,'(A)') &
    '  ============================================================='
    write(60,*)
    write(60,*)
    write(60,*)
    write(60,'(A)') &
    '  ============================================================='
    write(60,'(A)') &
    '                      CANONICAL MO INFO'
    write(60,'(A)') &
    '  ============================================================='
    do loop_i = 1, electron_count
      if (loop_i < electron_count) then
        call calc_S2HForb(loop_i)
        write(60,'(A,I3.3,A,F12.6)') &
        '  HOMO-', electron_count-loop_i, &
        '   energy / Eh                       ... ', orbE(loop_i)
        write(60,'(A,F12.6)') &
        '             <S**2> / hbar**2                  ... ', S__2orb
        write(60,'(A,F12.6)') &
        '             <Sz> / hbar                       ... ', Szorb
        write(60,'(A)') &
        '                 orb-in-atom  atom-in-mol      RE        IM'
        iatom = 1
        ishell = 1
        do loop_j = 1, sbdm
          if (basis_inf(loop_j)%atom /= iatom) then
            ishell = 1
            iatom = iatom + 1
          end if
          if (abs(AO2MO(loop_j,loop_i)) >= prtlev) then
            write(60,'(A,I2.2,A,I1,A,A2,A,I2.2,A,F9.6,A1,F9.6,A1)') &
            '             --   A #',ishell,' L=',basis_inf(loop_j)%L-1,'     ',&
            element_list(mol(basis_inf(loop_j)%atom)%atom_number),' #',&
            basis_inf(loop_j)%atom,'    (', real(AO2MO(loop_j,&
            loop_i)), ',', aimag(AO2MO(loop_j,loop_i)), ')'
          end if
          if (abs(AO2MO(sbdm+loop_j,loop_i)) >= prtlev) then
            write(60,'(A,I2.2,A,I1,A,A2,A,I2.2,A,F9.6,A1,F9.6,A1)') &
            '             --   B #',ishell,' L=',basis_inf(loop_j)%L-1,'     ',&
            element_list(mol(basis_inf(loop_j)%atom)%atom_number),' #',&
            basis_inf(loop_j)%atom,'    (', real(AO2MO(sbdm+loop_j,&
            loop_i)), ',', aimag(AO2MO(sbdm+loop_j,loop_i)), ')'
          end if
          ishell = ishell + 1
        end do
      else
        call calc_S2HForb(loop_i)
        write(60,'(A,A,F12.6)') &
        '  HOMO    ', '   energy / Eh                       ... ', orbE(loop_i)
        write(60,'(A,F12.6)') &
        '             <S**2> / hbar**2                  ... ', S__2orb
        write(60,'(A,F12.6)') &
        '             <Sz> / hbar                       ... ', Szorb
        write(60,'(A)') &
        '                 orb-in-atom  atom-in-mol      RE        IM'
        iatom = 1
        ishell = 1
        do loop_j = 1, sbdm
          if (basis_inf(loop_j)%atom /= iatom) then
            ishell = 1
            iatom = iatom + 1
          end if
          if (abs(AO2MO(loop_j,loop_i)) >= prtlev) then
            write(60,'(A,I2.2,A,I1,A,A2,A,I2.2,A,F9.6,A1,F9.6,A1)') &
            '             --   A #',ishell,' L=',basis_inf(loop_j)%L-1,'     ',&
            element_list(mol(basis_inf(loop_j)%atom)%atom_number),' #',&
            basis_inf(loop_j)%atom,'    (', real(AO2MO(loop_j,&
            loop_i)), ',', aimag(AO2MO(loop_j,loop_i)), ')'
          end if
          if (abs(AO2MO(sbdm+loop_j,loop_i)) >= prtlev) then
            write(60,'(A,I2.2,A,I1,A,A2,A,I2.2,A,F9.6,A1,F9.6,A1)') &
            '             --   B #',ishell,' L=',basis_inf(loop_j)%L-1,'     ',&
            element_list(mol(basis_inf(loop_j)%atom)%atom_number),' #',&
            basis_inf(loop_j)%atom,'    (', real(AO2MO(sbdm+loop_j,&
            loop_i)), ',', aimag(AO2MO(sbdm+loop_j,loop_i)), ')'
          end if
          ishell = ishell + 1
        end do
      end if
    end do
    do loop_i = electron_count+1, electron_count+5
      if (loop_i == electron_count + 1) then
        call calc_S2HForb(loop_i)
        write(60,'(A,A,F12.6)') &
        '  LUMO    ', '   energy / Eh                       ... ', orbE(loop_i)
        write(60,'(A,F12.6)') &
        '             <S**2> / hbar**2                  ... ', S__2orb
        write(60,'(A,F12.6)') &
        '             <Sz> / hbar                       ... ', Szorb
        write(60,'(A)') &
        '                 orb-in-atom  atom-in-mol      RE        IM'
        iatom = 1
        ishell = 1
        do loop_j = 1, sbdm
          if (basis_inf(loop_j)%atom /= iatom) then
            ishell = 1
            iatom = iatom + 1
          end if
          if (abs(AO2MO(loop_j,loop_i)) >= prtlev) then
            write(60,'(A,I2.2,A,I1,A,A2,A,I2.2,A,F9.6,A1,F9.6,A1)') &
            '             --   A #',ishell,' L=',basis_inf(loop_j)%L-1,'     ',&
            element_list(mol(basis_inf(loop_j)%atom)%atom_number),' #',&
            basis_inf(loop_j)%atom,'    (', real(AO2MO(loop_j,&
            loop_i)), ',', aimag(AO2MO(loop_j,loop_i)), ')'
          end if
          if (abs(AO2MO(sbdm+loop_j,loop_i)) >= prtlev) then
            write(60,'(A,I2.2,A,I1,A,A2,A,I2.2,A,F9.6,A1,F9.6,A1)') &
            '             --   B #',ishell,' L=',basis_inf(loop_j)%L-1,'     ',&
            element_list(mol(basis_inf(loop_j)%atom)%atom_number),' #',&
            basis_inf(loop_j)%atom,'    (', real(AO2MO(sbdm+loop_j,&
            loop_i)), ',', aimag(AO2MO(sbdm+loop_j,loop_i)), ')'
          end if
          ishell = ishell + 1
        end do
      else
        call calc_S2HForb(loop_i)
        write(60,'(A,I3.3,A,F12.6)') &
        '  LUMO+', loop_i-electron_count-1, &
        '   energy / Eh                       ... ', orbE(loop_i)
        write(60,'(A,F12.6)') &
        '             <S**2> / hbar**2                  ... ', S__2orb
        write(60,'(A,F12.6)') &
        '             <Sz> / hbar                       ... ', Szorb
        write(60,'(A)') &
        '                 orb-in-atom  atom-in-mol      RE        IM'
        iatom = 1
        ishell = 1
        do loop_j = 1, sbdm
          if (basis_inf(loop_j)%atom /= iatom) then
            ishell = 1
            iatom = iatom + 1
          end if
          if (abs(AO2MO(loop_j,loop_i)) >= prtlev) then
            write(60,'(A,I2.2,A,I1,A,A2,A,I2.2,A,F9.6,A1,F9.6,A1)') &
            '             --   A #',ishell,' L=',basis_inf(loop_j)%L-1,'     ',&
            element_list(mol(basis_inf(loop_j)%atom)%atom_number),' #',&
            basis_inf(loop_j)%atom,'    (', real(AO2MO(loop_j,&
            loop_i)), ',', aimag(AO2MO(loop_j,loop_i)), ')'
          end if
          if (abs(AO2MO(sbdm+loop_j,loop_i)) >= prtlev) then
            write(60,'(A,I2.2,A,I1,A,A2,A,I2.2,A,F9.6,A1,F9.6,A1)') &
            '             --   B #',ishell,' L=',basis_inf(loop_j)%L-1,'     ',&
            element_list(mol(basis_inf(loop_j)%atom)%atom_number),' #',&
            basis_inf(loop_j)%atom,'    (', real(AO2MO(sbdm+loop_j,&
            loop_i)), ',', aimag(AO2MO(sbdm+loop_j,loop_i)), ')'
          end if
          ishell = ishell + 1
        end do
      end if
    end do
    write(60,'(A)') &
    '  ============================================================='
    write(60,*)
    if (molden) then
      write(60,'(A)') '  dumping AO2MO to '//trim(address_job)//'.molden.input'
      call dump_molden()
      write(60,'(A)') '  complete!'
    end if
    write(60,'(A)') 'exit module SCF'
  end subroutine outprint
  !------------------------------------------------------------
  !> initialization of global variables
  subroutine glob_init(kill)
    implicit none
    logical,intent(in) :: kill   ! kill the process after SCF
    deallocate(exi_T_j, exi_V_j)
    deallocate(oper4, oper6)
    deallocate(isupp_ev_f, Rsd, DIISmat, ipiv)
    deallocate(rho_history, rho_pre, rho_pre_pre)
    deallocate(swint)
    deallocate(Fock1, Fock2HFexc, Fock2HFcol)
    if (allocated(Fock2KScor)) deallocate(Fock2KScor)
    if (allocated(Fock2KSexc)) deallocate(Fock2KSexc)
    deallocate(i_j, i_p2_j, i_V_j, Xm, exXm)
    if (s_h) then
      deallocate(c2s, exc2s)
    end if
    if (kill) then
      if (fx_id /= -1) call Fockxc_end()
      deallocate(AO2MO, rho_m, Fock, orbE, oper3)
    end if
    ndschwarz = .true.
    ini_rou = .true.
    if (DKH_order == 2) then
      deallocate(exSOC, AO2p2, evl_p2)
      if (SRTP_type) deallocate(exSR)
    end if
    if (kill) then
      call terminate('normal')
    else
      call terminate('keep')
    end if
  end subroutine glob_init
  
!------------------------------------------------------------
!> initial guess and generate density matix
  subroutine assign_rou()
    implicit none
    integer :: ploop_i,ploop_j,ploop_k! loop variables for subroutine assign_rou
    integer :: mat_dimension
    integer :: Na, Nb, degenlow, degenhigh, load, unload ! degenerat region
    real(dp) :: rdMO(5)
    character(len = 512) :: line_str
    if (ini_rou) then
      ini_rou = .false.
      allocate(rho_m(2*sbdm,2*sbdm))
      allocate(AO2MO(2*sbdm,2*fbdm))
      Nalpha = (electron_count-(spin_mult-1))/2 + (spin_mult-1)
      Nbeta = (electron_count-(spin_mult-1))/2
      
      ! read MO coefficient
      if (guess_type == 'gaussian') then
        allocate(Gaualpha(sbdm*sbdm))
        allocate(Gaubeta(sbdm*sbdm))
        allocate(AO2MOalpha(sbdm,sbdm))
        allocate(AO2MObeta(sbdm,sbdm))
        open(61, file = trim(address_job)//'.gaualpha', &
        status = 'old', action = 'read', iostat = ios)
        if (ios /= 0) then
          ios = system('rwfdump '//trim(address_job)//'.chk '&
          //trim(address_job)//'.gaualpha 524R')
          if (ios == -1) &
          call terminate('Cannot generate .gaualpha for initial density matrix')
          open(61, file = trim(address_job)//'.gaualpha', status = 'old', &
          action = 'read', iostat = ios)
          if (ios /= 0) &
          call terminate('Cannot open .gaualpha for initial density matrix')
        end if
        open(62, file = trim(address_job)//'.gaubeta', &
        status = 'old', action = 'read', iostat = ios)
        if (ios /= 0) then
          ios = system('rwfdump '//trim(address_job)//'.chk '&
          //trim(address_job)//'.gaubeta 526R')
          if (ios == -1) &
          call terminate('Cannot generate .gaubeta for initial density matrix')
          open(62, file = trim(address_job)//'.gaubeta', status = 'old', &
          action = 'read', iostat = ios)
          if (ios /= 0) &
          call terminate('Cannot open .gaubeta for initial density matrix')
        end if
        do
          read(61,'(a512)') line_str
          if (index(line_str,'Dump of file') /= 0) exit
        end do
        if (index(line_str(index(line_Str,'length')-4:&
        index(line_Str,'length')-2),'524') == 0) &
        call terminate('Gaussian alpha orbital coefficient should be RWF 524')
        read(line_str(index(line_Str,'(')-9:index(line_Str,'(')-2),&
        '(I8)',iostat = ios) mat_dimension
        if (ios /= 0) call terminate('Unmatched Gaualpha content')
        if (mat_dimension /= sbdm*sbdm) &
        call terminate('Gaussian alpha orbital dimension not match')
        do ploop_i=1,sbdm*sbdm/5
          read(61,*) rdMO
          Gaualpha((ploop_i-1)*5 + 1) = rdMO(1)
          Gaualpha((ploop_i-1)*5 + 2) = rdMO(2)
          Gaualpha((ploop_i-1)*5 + 3) = rdMO(3)
          Gaualpha((ploop_i-1)*5 + 4) = rdMO(4)
          Gaualpha((ploop_i-1)*5 + 5) = rdMO(5)
        end do
        ploop_i = ploop_i - 5
        if (mod(sbdm*sbdm,5) /= 0) then
          if (mod(sbdm*sbdm,5) == 1) read(61,*) Gaualpha(ploop_i*5 + 1)
          if (mod(sbdm*sbdm,5) == 2) read(61,*) Gaualpha(ploop_i*5 + 1),&
          Gaualpha(ploop_i*5 + 2)
          if (mod(sbdm*sbdm,5) == 3) read(61,*) Gaualpha(ploop_i*5 + 1),&
          Gaualpha(ploop_i*5 + 2),Gaualpha(ploop_i*5 + 3)
          if (mod(sbdm*sbdm,5) == 4) read(61,*) Gaualpha(ploop_i*5 + 1),&
          Gaualpha(ploop_i*5 + 2),Gaualpha(ploop_i*5 + 3),Gaualpha(ploop_i*5+4)
        end if
        AO2MOalpha = transpose(reshape(Gaualpha,[sbdm,sbdm]))
        close(61)
        do
          read(62,'(a512)') line_str
          if (index(line_str,'Dump of file') /= 0) exit
        end do
        if (index(line_str(index(line_Str,'length')-4:&
        index(line_Str,'length')-2),'526') == 0) &
        call terminate('Gaussian beta orbital coefficient should be RWF 526')
        read(line_str(index(line_Str,'(')-9:index(line_Str,'(')-2),&
        '(I8)',iostat = ios) mat_dimension
        if (ios /= 0) call terminate('Unmatched Gaubeta content')
        if (mat_dimension /= sbdm*sbdm) &
        call terminate('Gaussian beta orbital dimension not match')
        do ploop_i=1,sbdm*sbdm/5
          read(62,*) rdMO
          Gaubeta((ploop_i-1)*5 + 1) = rdMO(1)
          Gaubeta((ploop_i-1)*5 + 2) = rdMO(2)
          Gaubeta((ploop_i-1)*5 + 3) = rdMO(3)
          Gaubeta((ploop_i-1)*5 + 4) = rdMO(4)
          Gaubeta((ploop_i-1)*5 + 5) = rdMO(5)
        end do
        ploop_i = ploop_i - 5
        if (mod(sbdm*sbdm,5) /= 0) then
          if (mod(sbdm*sbdm,5) == 1) read(62,*) Gaubeta(ploop_i*5 + 1)
          if (mod(sbdm*sbdm,5) == 2) read(62,*) Gaubeta(ploop_i*5 + 1),&
          Gaubeta(ploop_i*5 + 2)
          if (mod(sbdm*sbdm,5) == 3) read(62,*) Gaubeta(ploop_i*5 + 1),&
          Gaubeta(ploop_i*5 + 2),Gaubeta(ploop_i*5 + 3)
          if (mod(sbdm*sbdm,5) == 4) read(62,*) Gaubeta(ploop_i*5 + 1),&
          Gaubeta(ploop_i*5 + 2),Gaubeta(ploop_i*5 + 3),Gaubeta(ploop_i*5 + 4)
        end if
        AO2MObeta = transpose(reshape(Gaubeta,[sbdm,sbdm]))
        close(62)
        deallocate(Gaualpha)
        deallocate(Gaubeta)
        rho_m = c0
        do ploop_i = 1,sbdm
          do ploop_j = 1,sbdm
            do ploop_k = 1,Nalpha
              rho_m(ploop_i,ploop_j) = rho_m(ploop_i,ploop_j) +&
              AO2MOalpha(ploop_k,ploop_i)*AO2MOalpha(ploop_k,ploop_j)
            end do
            do ploop_k = 1,Nbeta
              rho_m(sbdm+ploop_i,sbdm+ploop_j) = &
              rho_m(sbdm+ploop_i,sbdm+ploop_j) +&
              AO2MObeta(ploop_k,ploop_i)*AO2MObeta(ploop_k,ploop_j)
            end do
          end do
        end do
        deallocate(AO2MOalpha)
        deallocate(AO2MObeta)
      else if (guess_type == 'read') then
        call load_matrix('ao2mo', AO2MO, ploop_i, ploop_j)
        if (ploop_i /= 2*sbdm .or. ploop_j /= 2*fbdm) call terminate(&
        'basis dimension in .ao2mo file mismatch with current job')
        rho_m = c0
        do ploop_i = 1, 2*sbdm
          do ploop_j = 1, 2*sbdm
            do ploop_k = 1, electron_count
              rho_m(ploop_i,ploop_j) = rho_m(ploop_i,ploop_j) + &
              AO2MO(ploop_i,ploop_k)*conjg(AO2MO(ploop_j,ploop_k))
            end do
          end do
        end do
      end if
      rho_pre = rho_m
    else
      rho_m = c0
      do ploop_i = 1, 2*sbdm
        do ploop_j = 1, 2*sbdm
          do ploop_k = 1, electron_count
            rho_m(ploop_i,ploop_j) = rho_m(ploop_i,ploop_j) + &
            AO2MO(ploop_i,ploop_k)*conjg(AO2MO(ploop_j,ploop_k))
          end do
        end do
      end do
      call calc_S2HF()
      if (keepspin .and. abs((totalpha-totbeta) - &
      real(Nalpha-Nbeta,dp)) > 1.0.and. abs((totalpha-real(Nalpha,dp)) - &
      (real(Nbeta,dp)-totbeta)) < 0.1) then
        write(60,'(A)') &
        '  -- order of degenerate frontier alpha/beta orbitals changed'
        do ploop_i = electron_count-1, 1, -1
          if (abs(orbE(ploop_i)-orbE(electron_count)) &
          / abs(orbE(electron_count)) > 0.04) then
            degenlow = ploop_i + 1
            exit
          end if
        end do
        do ploop_i = electron_count+1, 6*electron_count
          if (abs(orbE(ploop_i)-orbE(electron_count)) &
          / abs(orbE(electron_count)) > 0.04) then
            degenhigh = ploop_i - 1
            exit
          end if
        end do
        write(60,'(A,I3,A,I3)') &
        '  -- -- degenerate space: ',degenlow,' - ',degenhigh
        allocate(rotation(2*fbdm,degenhigh-degenlow+1))
        Na = 0
        Nb = 0
        do ploop_k = 1, degenlow-1
          call calc_S2HForb(ploop_k)
          if (Szorb > 0.0) then
            Na = Na + 1
          else
            Nb = Nb + 1
          end if
        end do
        load = degenlow
        unload = degenhigh
        do ploop_k = degenlow, degenhigh
          call calc_S2HForb(ploop_k)
          if (Szorb > 0.0) then
            if (Na < Nalpha) then
              do ploop_i = 1, 2*fbdm
                rotation(ploop_i,load-degenlow+1) = oper3(ploop_i,ploop_k)
              end do
              write(60,'(A,I3,A,I3,A)') &
              '  -- -- load ',ploop_k,' -> ',load,' (alpha)'
              load = load + 1
              Na = Na + 1
            else
              do ploop_i = 1, 2*fbdm
                rotation(ploop_i,unload-degenlow+1) = oper3(ploop_i,ploop_k)
              end do
              write(60,'(A,I3,A,I3,A)') &
              '  -- -- unload ',ploop_k,' -> ',unload,' (alpha)'
              unload = unload - 1
            end if                            
          else
              if (Nb < Nbeta) then
              do ploop_i = 1, 2*fbdm
                rotation(ploop_i,load-degenlow+1) = oper3(ploop_i,ploop_k)
              end do
              write(60,'(A,I3,A,I3,A)') &
              '  -- -- load ',ploop_k,' -> ',load,' (beta)'
              load = load + 1
              Nb = Nb + 1
            else
              do ploop_i = 1, 2*fbdm
                rotation(ploop_i,unload-degenlow+1) = oper3(ploop_i,ploop_k)
              end do
              write(60,'(A,I3,A,I3,A)') &
              '  -- -- unload ',ploop_k,' -> ',unload,' (beta)'
              unload = unload - 1
            end if 
          end if
        end do
        do ploop_k = degenlow, degenhigh
          do ploop_i = 1, 2*fbdm
            oper3(ploop_i,ploop_k) = rotation(ploop_i,ploop_k-degenlow+1)
          end do
        end do
        deallocate(rotation)
        call zgemm( 'N', 'N', 2*sbdm, 2*fbdm, 2*fbdm, c1, exXm, 2*sbdm, &
        oper3, 2*fbdm, c0, AO2MO, 2*sbdm)
        rho_m = c0
        do ploop_i = 1, 2*sbdm
          do ploop_j = 1, 2*sbdm
            do ploop_k = 1, electron_count
              rho_m(ploop_i,ploop_j) = rho_m(ploop_i,ploop_j) + &
              AO2MO(ploop_i,ploop_k)*conjg(AO2MO(ploop_j,ploop_k))
            end do
          end do
        end do
      end if
      if (loop_i /= 1) then
        maxDP = (real(rho_m(1,1))-real(rho_pre(1,1)))**2 + &
        (aimag(rho_m(1,1))-aimag(rho_pre(1,1)))**2
        RMSDP = 0.0_dp
        do ploop_i = 1, 2*sbdm
          do ploop_j = 1, 2*sbdm
            if ((real(rho_m(ploop_i,ploop_j))-&
            real(rho_pre(ploop_i,ploop_j)))**2 + (aimag(rho_m(ploop_i,ploop_j))&
            -aimag(rho_pre(ploop_i,ploop_j)))**2 > maxDP) then
              maxDP = (real(rho_m(ploop_i,ploop_j))-&
              real(rho_pre(ploop_i,ploop_j)))**2 + &
              (aimag(rho_m(ploop_i,ploop_j))-aimag(rho_pre(ploop_i,ploop_j)))**2
            end if
            RMSDP = RMSDP + &
            (real(rho_m(ploop_i,ploop_j))-real(rho_pre(ploop_i,ploop_j)))**2 + &
            (aimag(rho_m(ploop_i,ploop_j))-aimag(rho_pre(ploop_i,ploop_j)))**2
          end do
        end do
        maxDP = dsqrt(maxDP)
        RMSDP = RMSDP / (4*sbdm*sbdm)
        RMSDP = dsqrt(RMSDP)
        write(60,'(A,E10.3)') '  -- maxDP                  ', maxDP
        write(60,'(A,E10.3)') '  -- RMSDP                  ', RMSDP
      end if
      ! calc <S**2>
      call calc_S2HF()
      write(60,'(A,F9.5)') '  -- <S**2>                ',S__2
      write(60,'(A,F9.5)') '  -- total alpha electron  ',totalpha
      write(60,'(A,F9.5)') '  -- total beta electron   ',totbeta
    end if
  end subroutine assign_rou
  
!------------------------------------------------------------
!> construct one electron Fock matrix
  subroutine Fock1e()
    implicit none
    real(dp)    :: temp_pool(fbdm, fbdm)
    integer     :: floop_i, floop_j   ! loop variables for subroutine Fock1e
    complex(dp) :: itm(fbdm)
    real(dp)    :: edc(fbdm)          ! (p2+c2)^0.5
    allocate(Fock1(2*fbdm,2*fbdm), exi_T_j(2*fbdm,2*fbdm), &
    exi_V_j(2*fbdm,2*fbdm), source = c0)
    if (DKH_order == 0) then
      Fock1(1:fbdm,1:fbdm) = 0.5_dp*i_p2_j + i_V_j
      Fock1(fbdm+1:2*fbdm,fbdm+1:2*fbdm) = 0.5_dp*i_p2_j + i_V_j
      exi_T_j(1:fbdm,1:fbdm) = 0.5_dp*i_p2_j
      exi_T_j(fbdm+1:2*fbdm,fbdm+1:2*fbdm) = 0.5_dp*i_p2_j
      exi_V_j(1:fbdm,1:fbdm) = i_V_j
      exi_V_j(fbdm+1:2*fbdm,fbdm+1:2*fbdm) = i_V_j
    else if (DKH_order == 2) then
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
      allocate(exAO2p2(2*fbdm,2*fbdm), exSOC(2*fbdm,2*fbdm), source = c0)
      edc = dsqrt(evl_p2 + c2)
      !----------------------
      forall (floop_i = 1:fbdm) Ap(floop_i, floop_i) = &
      dsqrt( (dsqrt(evl_p2(floop_i)/c2 + 1.0_dp) + 1.0_dp) / &
      (2.0_dp * dsqrt(evl_p2(floop_i)/c2 + 1.0_dp)) )
      !----------------------
      forall (floop_i = 1:fbdm) ApRp(floop_i, floop_i) = &
      Ap(floop_i,floop_i) / (edc(floop_i)+c)
      !----------------------
      itm = (1.0_dp + evl_p2 / (4.0_dp * c2) + QED_rad) * c1
      forall (floop_i = 1:fbdm)
        SRp(floop_i, floop_i) = itm(floop_i)
        SRp(fbdm+floop_i, fbdm+floop_i)   = itm(floop_i)
      end forall
      !----------------------
      ! Ap V Ap
      call dgemm( 'N', 'T', fbdm, fbdm, fbdm, 1.0_dp, Ap, fbdm, &
      AO2p2, fbdm, 0.0_dp, oper2, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
      i_V_j, fbdm, 0.0_dp, oper1, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper1, fbdm, &
      AO2p2, fbdm, 0.0_dp, oper2, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
      Ap, fbdm, 0.0_dp, oper1, fbdm)
      AVA(1:fbdm,1:fbdm) = AVA(1:fbdm,1:fbdm) + oper1 * c1
      AVA(fbdm+1:2*fbdm,fbdm+1:2*fbdm) = &
      AVA(fbdm+1:2*fbdm,fbdm+1:2*fbdm) + oper1 * c1
      exi_V_j(1:fbdm,1:fbdm) = i_V_j
      exi_V_j(fbdm+1:2*fbdm,fbdm+1:2*fbdm) = i_V_j
      !----------------------
      ! ApRp pxVpx+pyVpy+pzVpz ApRp
      temp_pool = pxVpx+pyVpy+pzVpz
      call dgemm( 'N', 'T', fbdm, fbdm, fbdm, 1.0_dp, ApRp, fbdm, &
      AO2p2, fbdm, 0.0_dp, oper2, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
      temp_pool, fbdm, 0.0_dp, oper1, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper1, fbdm, &
      AO2p2, fbdm, 0.0_dp, oper2, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
      ApRp, fbdm, 0.0_dp, oper1, fbdm)
      ARVRA(1:fbdm,1:fbdm) = ARVRA(1:fbdm,1:fbdm) + oper1 * c1
      ARVRA(fbdm+1:2*fbdm,fbdm+1:2*fbdm) = &
      ARVRA(fbdm+1:2*fbdm,fbdm+1:2*fbdm) + oper1 * c1
      !----------------------
      ! ApRp pxVpy-pyVpx ApRp
      temp_pool = pxVpy-pyVpx
      call dgemm( 'N', 'T', fbdm, fbdm, fbdm, 1.0_dp, ApRp, fbdm, &
      AO2p2, fbdm, 0.0_dp, oper2, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
      temp_pool, fbdm, 0.0_dp, oper1, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper1, fbdm, &
      AO2p2, fbdm, 0.0_dp, oper2, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
      ApRp, fbdm, 0.0_dp, oper1, fbdm)
      ARVRA(1:fbdm,1:fbdm) = ARVRA(1:fbdm,1:fbdm) + oper1 * ci
      ARVRA(fbdm+1:2*fbdm,fbdm+1:2*fbdm) = &
      ARVRA(fbdm+1:2*fbdm,fbdm+1:2*fbdm) - oper1 * ci
      !----------------------
      ! ApRp pzVpx-pxVpz ApRp
      temp_pool = pzVpx-pxVpz
      call dgemm( 'N', 'T', fbdm, fbdm, fbdm, 1.0_dp, ApRp, fbdm, &
      AO2p2, fbdm, 0.0_dp, oper2, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
      temp_pool, fbdm, 0.0_dp, oper1, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper1, fbdm, &
      AO2p2, fbdm, 0.0_dp, oper2, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
      ApRp, fbdm, 0.0_dp, oper1, fbdm)
      ARVRA(fbdm+1:2*fbdm,1:fbdm) = ARVRA(fbdm+1:2*fbdm,1:fbdm) - oper1 * c1
      ARVRA(1:fbdm,fbdm+1:2*fbdm) = ARVRA(1:fbdm,fbdm+1:2*fbdm) + oper1 * c1
      !----------------------
      ! ApRp pyVpz-pzVpy ApRp
      temp_pool = pyVpz-pzVpy
      call dgemm( 'N', 'T', fbdm, fbdm, fbdm, 1.0_dp, ApRp, fbdm, &
      AO2p2, fbdm, 0.0_dp, oper2, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
      temp_pool, fbdm, 0.0_dp, oper1, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper1, fbdm, &
      AO2p2, fbdm, 0.0_dp, oper2, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
      ApRp, fbdm, 0.0_dp, oper1, fbdm)
      ARVRA(fbdm+1:2*fbdm,1:fbdm) = ARVRA(fbdm+1:2*fbdm,1:fbdm) + oper1 * ci
      ARVRA(1:fbdm,fbdm+1:2*fbdm) = ARVRA(1:fbdm,fbdm+1:2*fbdm) + oper1 * ci
      !----------------------
      Fock1 = Fock1 + AVA
      if (SRTP_type) then
        call zgemm( 'N', 'N', 2*fbdm, 2*fbdm, 2*fbdm, c1, SRp, 2*fbdm, &
        ARVRA, 2*fbdm, c0, oper3, 2*fbdm)
        call zgemm( 'N', 'N', 2*fbdm, 2*fbdm, 2*fbdm, c1, oper3, 2*fbdm, &
        SRp, 2*fbdm, c0, oper4, 2*fbdm)
        Fock1 = Fock1 + oper4
        exSOC = ARVRA
      else
        Fock1 = Fock1 + ARVRA
        exSOC = ARVRA
      end if
      if (SRTP_type) then
        allocate(exSR(2*fbdm,2*fbdm), source = c0)
        !----------------------
        ! ApRp (px3Vpx+py3Vpy+pz3Vpz+pxVpx3+pyVpy3+pzVpz3) ApRp
        temp_pool = px3Vpx+py3Vpy+pz3Vpz+transpose(px3Vpx+py3Vpy+pz3Vpz)
        call dgemm( 'N', 'T', fbdm, fbdm, fbdm, 1.0_dp, ApRp, fbdm, &
        AO2p2, fbdm, 0.0_dp, oper2, fbdm)
        call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
        temp_pool, fbdm, 0.0_dp, oper1, fbdm)
        call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper1, fbdm, &
        AO2p2, fbdm, 0.0_dp, oper2, fbdm)
        call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
        ApRp, fbdm, 0.0_dp, oper1, fbdm)
        Fock1(1:fbdm,1:fbdm) = Fock1(1:fbdm,1:fbdm) - oper1/(2.0_dp*c2) * c1
        Fock1(fbdm+1:2*fbdm,fbdm+1:2*fbdm) = &
        Fock1(fbdm+1:2*fbdm,fbdm+1:2*fbdm) - oper1/(2.0_dp*c2) * c1
        exSR(1:fbdm,1:fbdm) = exSR(1:fbdm,1:fbdm) - oper1/(2.0_dp*c2) * c1
        exSR(fbdm+1:2*fbdm,fbdm+1:2*fbdm) = &
        exSR(fbdm+1:2*fbdm,fbdm+1:2*fbdm) - oper1/(2.0_dp*c2) * c1
        !----------------------
        ! ApRp (px3Vpy-py3Vpx+pxVpy3-pyVpx3) ApRp
        temp_pool = px3Vpy-py3Vpx+transpose(py3Vpx-px3Vpy)
        call dgemm( 'N', 'T', fbdm, fbdm, fbdm, 1.0_dp, ApRp, fbdm, &
        AO2p2, fbdm, 0.0_dp, oper2, fbdm)
        call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
        temp_pool, fbdm, 0.0_dp, oper1, fbdm)
        call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper1, fbdm, &
        AO2p2, fbdm, 0.0_dp, oper2, fbdm)
        call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
        ApRp, fbdm, 0.0_dp, oper1, fbdm)
        Fock1(1:fbdm,1:fbdm) = Fock1(1:fbdm,1:fbdm) - oper1/(2.0_dp*c2) * ci
        Fock1(fbdm+1:2*fbdm,fbdm+1:2*fbdm) = &
        Fock1(fbdm+1:2*fbdm,fbdm+1:2*fbdm) + oper1/(2.0_dp*c2) * ci
        exSR(1:fbdm,1:fbdm) = exSR(1:fbdm,1:fbdm) - oper1/(2.0_dp*c2) * ci
        exSR(fbdm+1:2*fbdm,fbdm+1:2*fbdm) = &
        exSR(fbdm+1:2*fbdm,fbdm+1:2*fbdm) + oper1/(2.0_dp*c2) * ci
        !----------------------
        ! ApRp (pz3Vpx-px3Vpz+pzVpx3-pxVpz3) ApRp
        temp_pool = pz3Vpx-px3Vpz+transpose(px3Vpz-pz3Vpx)
        call dgemm( 'N', 'T', fbdm, fbdm, fbdm, 1.0_dp, ApRp, fbdm, &
        AO2p2, fbdm, 0.0_dp, oper2, fbdm)
        call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
        temp_pool, fbdm, 0.0_dp, oper1, fbdm)
        call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper1, fbdm, &
        AO2p2, fbdm, 0.0_dp, oper2, fbdm)
        call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
        ApRp, fbdm, 0.0_dp, oper1, fbdm)
        Fock1(fbdm+1:2*fbdm,1:fbdm) = &
        Fock1(fbdm+1:2*fbdm,1:fbdm) + oper1/(2.0_dp*c2) * c1
        Fock1(1:fbdm,fbdm+1:2*fbdm) = &
        Fock1(1:fbdm,fbdm+1:2*fbdm) - oper1/(2.0_dp*c2) * c1
        exSR(fbdm+1:2*fbdm,1:fbdm) = &
        exSR(fbdm+1:2*fbdm,1:fbdm) + oper1/(2.0_dp*c2) * c1
        exSR(1:fbdm,fbdm+1:2*fbdm) = &
        exSR(1:fbdm,fbdm+1:2*fbdm) - oper1/(2.0_dp*c2) * c1
        !----------------------
        ! ApRp (py3Vpz-pz3Vpy+pyVpz3-pzVpy3) ApRp
        temp_pool = py3Vpz-pz3Vpy+transpose(pz3Vpy-py3Vpz)
        call dgemm( 'N', 'T', fbdm, fbdm, fbdm, 1.0_dp, ApRp, fbdm, &
        AO2p2, fbdm, 0.0_dp, oper2, fbdm)
        call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
        temp_pool, fbdm, 0.0_dp, oper1, fbdm)
        call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper1, fbdm, &
        AO2p2, fbdm, 0.0_dp, oper2, fbdm)
        call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
        ApRp, fbdm, 0.0_dp, oper1, fbdm)
        Fock1(fbdm+1:2*fbdm,1:fbdm) = &
        Fock1(fbdm+1:2*fbdm,1:fbdm) - oper1/(2.0_dp*c2) * ci
        Fock1(1:fbdm,fbdm+1:2*fbdm) = &
        Fock1(1:fbdm,fbdm+1:2*fbdm) - oper1/(2.0_dp*c2) * ci
        exSR(fbdm+1:2*fbdm,1:fbdm) = &
        exSR(fbdm+1:2*fbdm,1:fbdm) - oper1/(2.0_dp*c2) * ci
        exSR(1:fbdm,fbdm+1:2*fbdm) = &
        exSR(1:fbdm,fbdm+1:2*fbdm) - oper1/(2.0_dp*c2) * ci
      end if
      !----------------------
      ! start building Fock1
      do floop_j = 1, fbdm
        do floop_i = 1, fbdm
          Ve(floop_i,floop_j) = i_V_j(floop_i,floop_j)/&
          (c * (edc(floop_i) + edc(floop_j)))
        end do
      end do
      do floop_j = 1, fbdm
        do floop_i = 1, fbdm
          pxVepx(floop_i,floop_j) = pxVpx(floop_i,floop_j)/&
          (c * (edc(floop_i) + edc(floop_j)))
        end do
      end do
      do floop_j = 1, fbdm
        do floop_i = 1, fbdm
          pyVepy(floop_i,floop_j) = pyVpy(floop_i,floop_j)/&
          (c * (edc(floop_i) + edc(floop_j)))
        end do
      end do
      do floop_j = 1, fbdm
        do floop_i = 1, fbdm
          pzVepz(floop_i,floop_j) = pzVpz(floop_i,floop_j)/&
          (c * (edc(floop_i) + edc(floop_j)))
        end do
      end do
      do floop_j = 1, fbdm
        do floop_i = 1, fbdm
          pxVepy(floop_i,floop_j) = pxVpy(floop_i,floop_j)/&
          (c * (edc(floop_i) + edc(floop_j)))
        end do
      end do
      do floop_j = 1, fbdm
        do floop_i = 1, fbdm
          pyVepx(floop_i,floop_j) = pyVpx(floop_i,floop_j)/&
          (c * (edc(floop_i) + edc(floop_j)))
        end do
      end do
      do floop_j = 1, fbdm
        do floop_i = 1, fbdm
          pxVepz(floop_i,floop_j) = pxVpz(floop_i,floop_j)/&
          (c * (edc(floop_i) + edc(floop_j)))
        end do
      end do
      do floop_j = 1, fbdm
        do floop_i = 1, fbdm
          pzVepx(floop_i,floop_j) = pzVpx(floop_i,floop_j)/&
          (c * (edc(floop_i) + edc(floop_j)))
        end do
      end do
      do floop_j = 1, fbdm
        do floop_i = 1, fbdm
          pyVepz(floop_i,floop_j) = pyVpz(floop_i,floop_j)/&
          (c * (edc(floop_i) + edc(floop_j)))
        end do
      end do
      do floop_j = 1, fbdm
        do floop_i = 1, fbdm
          pzVepy(floop_i,floop_j) = pzVpy(floop_i,floop_j)/&
          (c * (edc(floop_i) + edc(floop_j)))
        end do
      end do
      !----------------------
      ! Ap Ve Ap
      call dgemm( 'N', 'T', fbdm, fbdm, fbdm, 1.0_dp, Ap, fbdm, &
      AO2p2, fbdm, 0.0_dp, oper2, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
      Ve, fbdm, 0.0_dp, oper1, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper1, fbdm, &
      AO2p2, fbdm, 0.0_dp, oper2, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
      Ap, fbdm, 0.0_dp, oper1, fbdm)
      AVeA(1:fbdm,1:fbdm) = oper1 * c1
      AVeA(fbdm+1:2*fbdm,fbdm+1:2*fbdm) = oper1 * c1
      !----------------------
      ! ApRp pxVepx+pyVepy+pzVepz ApRp
      temp_pool = pxVepx+pyVepy+pzVepz
      call dgemm( 'N', 'T', fbdm, fbdm, fbdm, 1.0_dp, ApRp, fbdm, &
      AO2p2, fbdm, 0.0_dp, oper2, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
      temp_pool, fbdm, 0.0_dp, oper1, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper1, fbdm, &
      AO2p2, fbdm, 0.0_dp, oper2, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
      ApRp, fbdm, 0.0_dp, oper1, fbdm)
      ARVeRA(1:fbdm,1:fbdm) = ARVeRA(1:fbdm,1:fbdm) + oper1 * c1
      ARVeRA(fbdm+1:2*fbdm,fbdm+1:2*fbdm) = &
      ARVeRA(fbdm+1:2*fbdm,fbdm+1:2*fbdm) + oper1 * c1
      !----------------------
      ! ApRp pxVepy-pyVepx ApRp
      temp_pool = pxVepy-pyVepx
      call dgemm( 'N', 'T', fbdm, fbdm, fbdm, 1.0_dp, ApRp, fbdm, &
      AO2p2, fbdm, 0.0_dp, oper2, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
      temp_pool, fbdm, 0.0_dp, oper1, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper1, fbdm, &
      AO2p2, fbdm, 0.0_dp, oper2, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
      ApRp, fbdm, 0.0_dp, oper1, fbdm)
      ARVeRA(1:fbdm,1:fbdm) = ARVeRA(1:fbdm,1:fbdm) + oper1 * ci
      ARVeRA(fbdm+1:2*fbdm,fbdm+1:2*fbdm) = &
      ARVeRA(fbdm+1:2*fbdm,fbdm+1:2*fbdm) - oper1 * ci
      !----------------------
      ! ApRp pzVepx-pxVepz ApRp
      temp_pool = pzVepx-pxVepz
      call dgemm( 'N', 'T', fbdm, fbdm, fbdm, 1.0_dp, ApRp, fbdm, &
      AO2p2, fbdm, 0.0_dp, oper2, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
      temp_pool, fbdm, 0.0_dp, oper1, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper1, fbdm, &
      AO2p2, fbdm, 0.0_dp, oper2, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
      ApRp, fbdm, 0.0_dp, oper1, fbdm)
      ARVeRA(fbdm+1:2*fbdm,1:fbdm) = ARVeRA(fbdm+1:2*fbdm,1:fbdm) - oper1 * c1
      ARVeRA(1:fbdm,fbdm+1:2*fbdm) = ARVeRA(1:fbdm,fbdm+1:2*fbdm) + oper1 * c1
      !----------------------
      ! ApRp pyVepz-pzVepy ApRp
      temp_pool = pyVepz-pzVepy
      call dgemm( 'N', 'T', fbdm, fbdm, fbdm, 1.0_dp, ApRp, fbdm, &
      AO2p2, fbdm, 0.0_dp, oper2, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
      temp_pool, fbdm, 0.0_dp, oper1, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper1, fbdm, &
      AO2p2, fbdm, 0.0_dp, oper2, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
      ApRp, fbdm, 0.0_dp, oper1, fbdm)
      ARVeRA(fbdm+1:2*fbdm,1:fbdm) = ARVeRA(fbdm+1:2*fbdm,1:fbdm) + oper1 * ci
      ARVeRA(1:fbdm,fbdm+1:2*fbdm) = ARVeRA(1:fbdm,fbdm+1:2*fbdm) + oper1 * ci
      !----------------------
      ! term of order c^-4, negative terms, no RI insertion
      forall (floop_i = 1:2*fbdm)
        SRp(floop_i, floop_i) = 0.5_dp * c1
      end forall
      call zgemm( 'N', 'N', 2*fbdm, 2*fbdm, 2*fbdm, c1, ARVeRA, 2*fbdm, &
      oper5, 2*fbdm, c0, oper3, 2*fbdm)
      call zgemm( 'N', 'N', 2*fbdm, 2*fbdm, 2*fbdm, c1, oper3, 2*fbdm, &
      AVA, 2*fbdm, c0, oper4, 2*fbdm)
      Fock1 = Fock1 - oper4
      call zgemm( 'N', 'N', 2*fbdm, 2*fbdm, 2*fbdm, c1, ARVRA, 2*fbdm, &
      oper5, 2*fbdm, c0, oper3, 2*fbdm)
      call zgemm( 'N', 'N', 2*fbdm, 2*fbdm, 2*fbdm, c1, oper3, 2*fbdm, &
      AVeA, 2*fbdm, c0, oper4, 2*fbdm)
      Fock1 = Fock1 - oper4
      call zgemm( 'N', 'N', 2*fbdm, 2*fbdm, 2*fbdm, c1, AVeA, 2*fbdm, &
      oper5, 2*fbdm, c0, oper3, 2*fbdm)
      call zgemm( 'N', 'N', 2*fbdm, 2*fbdm, 2*fbdm, c1, oper3, 2*fbdm, &
      ARVRA, 2*fbdm, c0, oper4, 2*fbdm)
      Fock1 = Fock1 - oper4
      call zgemm( 'N', 'N', 2*fbdm, 2*fbdm, 2*fbdm, c1, AVA, 2*fbdm, &
      oper5, 2*fbdm, c0, oper3, 2*fbdm)
      call zgemm( 'N', 'N', 2*fbdm, 2*fbdm, 2*fbdm, c1, oper3, 2*fbdm, &
      ARVeRA, 2*fbdm, c0, oper4, 2*fbdm)
      Fock1 = Fock1 - oper4
      !----------------------
      ! term of order c^-4, positive terms
      ! RI insertion: I = sigma.Pi(1/P^2)sigma.Pj
      ! oper5 is half of 1/P^2 or P^2, depends on whether RI insert or extract
      forall (floop_i = 1:fbdm) ! RI insert
        oper5(floop_i,floop_i) = &
        0.5_dp * ((edc(floop_i)+c)**2/evl_p2(floop_i)) * c1
        oper5(fbdm+floop_i,fbdm+floop_i) = &
        0.5_dp * ((edc(floop_i)+c)**2/evl_p2(floop_i)) * c1
      end forall
      call zgemm( 'N', 'N', 2*fbdm, 2*fbdm, 2*fbdm, c1, ARVeRA, 2*fbdm, &
      oper5, 2*fbdm, c0, oper3, 2*fbdm)
      call zgemm( 'N', 'N', 2*fbdm, 2*fbdm, 2*fbdm, c1, oper3, 2*fbdm, &
      ARVRA, 2*fbdm, c0, oper4, 2*fbdm)
      Fock1 = Fock1 + oper4
      call zgemm( 'N', 'N', 2*fbdm, 2*fbdm, 2*fbdm, c1, ARVRA, 2*fbdm, &
      oper5, 2*fbdm, c0, oper3, 2*fbdm)
      call zgemm( 'N', 'N', 2*fbdm, 2*fbdm, 2*fbdm, c1, oper3, 2*fbdm, &
      ARVeRA, 2*fbdm, c0, oper4, 2*fbdm)
      Fock1 = Fock1 + oper4
      forall (floop_i = 1:fbdm) ! RI extract
        oper5(floop_i,floop_i) = &
        0.5_dp * (evl_p2(floop_i)/(edc(floop_i)+c)**2) * c1
        oper5(fbdm+floop_i,fbdm+floop_i) = &
        0.5_dp * (evl_p2(floop_i)/(edc(floop_i)+c)**2) * c1
      end forall
      call zgemm( 'N', 'N', 2*fbdm, 2*fbdm, 2*fbdm, c1, AVeA, 2*fbdm, &
      oper5, 2*fbdm, c0, oper3, 2*fbdm)
      call zgemm( 'N', 'N', 2*fbdm, 2*fbdm, 2*fbdm, c1, oper3, 2*fbdm, &
      AVA, 2*fbdm, c0, oper4, 2*fbdm)
      Fock1 = Fock1 + oper4
      call zgemm( 'N', 'N', 2*fbdm, 2*fbdm, 2*fbdm, c1, AVA, 2*fbdm, &
      oper5, 2*fbdm, c0, oper3, 2*fbdm)
      call zgemm( 'N', 'N', 2*fbdm, 2*fbdm, 2*fbdm, c1, oper3, 2*fbdm, &
      AVeA, 2*fbdm, c0, oper4, 2*fbdm)
      Fock1 = Fock1 + oper4
      !----------------------
      ! kinetic energy
      forall (floop_i = 1:fbdm) ! RI extract
        Fock1(floop_i,floop_i) = &
        Fock1(floop_i,floop_i) + (c*edc(floop_i)-c2) * c1
        Fock1(fbdm+floop_i,fbdm+floop_i) = &
        Fock1(fbdm+floop_i,fbdm+floop_i) + (c*edc(floop_i)-c2) * c1
        exi_T_j(floop_i,floop_i) = &
        exi_T_j(floop_i,floop_i) + (c*edc(floop_i)-c2) * c1
        exi_T_j(fbdm+floop_i,fbdm+floop_i) = &
        exi_T_j(fbdm+floop_i,fbdm+floop_i) + (c*edc(floop_i)-c2) * c1
      end forall
      !----------------------
      ! transform from p^2 eigenbasis to orthogonal normalized AO basis
      exAO2p2(1:fbdm,1:fbdm) = AO2p2 * c1
      exAO2p2(fbdm+1:2*fbdm,fbdm+1:2*fbdm) = AO2p2 * c1
      call zgemm( 'N', 'N', 2*fbdm, 2*fbdm, 2*fbdm, c1, &
      exAO2p2, 2*fbdm, Fock1, 2*fbdm, c0, oper3, 2*fbdm)
      call zgemm( 'N', 'T', 2*fbdm, 2*fbdm, 2*fbdm, c1, &
      oper3, 2*fbdm, exAO2p2, 2*fbdm, c0, Fock1, 2*fbdm)
      call zgemm( 'N', 'N', 2*fbdm, 2*fbdm, 2*fbdm, c1, &
      exAO2p2, 2*fbdm, exi_T_j, 2*fbdm, c0, oper3, 2*fbdm)
      call zgemm( 'N', 'T', 2*fbdm, 2*fbdm, 2*fbdm, c1, &
      oper3, 2*fbdm, exAO2p2, 2*fbdm, c0, exi_T_j, 2*fbdm)
      call zgemm( 'N', 'N', 2*fbdm, 2*fbdm, 2*fbdm, c1, &
      exAO2p2, 2*fbdm, exi_V_j, 2*fbdm, c0, oper3, 2*fbdm)
      call zgemm( 'N', 'T', 2*fbdm, 2*fbdm, 2*fbdm, c1, &
      oper3, 2*fbdm, exAO2p2, 2*fbdm, c0, exi_V_j, 2*fbdm)
      call zgemm( 'N', 'N', 2*fbdm, 2*fbdm, 2*fbdm, c1, &
      exAO2p2, 2*fbdm, exSOC, 2*fbdm, c0, oper3, 2*fbdm)
      call zgemm( 'N', 'T', 2*fbdm, 2*fbdm, 2*fbdm, c1, &
      oper3, 2*fbdm, exAO2p2, 2*fbdm, c0, exSOC, 2*fbdm)
      if (SRTP_type) then
        call zgemm( 'N', 'N', 2*fbdm, 2*fbdm, 2*fbdm, c1, &
        exAO2p2, 2*fbdm, exSR, 2*fbdm, c0, oper3, 2*fbdm)
        call zgemm( 'N', 'T', 2*fbdm, 2*fbdm, 2*fbdm, c1, &
        oper3, 2*fbdm, exAO2p2, 2*fbdm, c0, exSR, 2*fbdm)
      end if
      deallocate(Ve, pxVepx, pyVepy, pzVepz, pxVepy, pyVepx)
      deallocate(pxVepz, pzVepx, pyVepz, pzVepy, ARVeRA, AVeA)
      deallocate(oper1, oper2, oper3, oper4, oper5)
      deallocate(Ap, ApRp, SRp, ARVRA, AVA, exAO2p2)
      ! deallocate one-electron ints to reduce memory usage
      deallocate(pxVpx, pyVpy, pzVpz, pxVpy, pyVpx)
      deallocate(pxVpz, pzVpx, pyVpz, pzVpy)
      if (SRTP_type) then
        deallocate(px3Vpx, py3Vpy, pz3Vpz, px3Vpy, py3Vpx)
        deallocate(px3Vpz, pz3Vpx, py3Vpz, pz3Vpy)
      end if
    end if
  end subroutine Fock1e
  
!------------------------------------------------------------
!> construct Fock2 matrices (Fock2HFcol, Fock2HFexc, Fock2KSexc, Fock2KScor)
!!
!! "direct" calculation to avoid memory overflow and large amount of disk r&w
!!
!! but 2 electron int should be calculated in each SCF iteration
  subroutine Fock2e()
    implicit none
    integer  :: i, j                ! for parallel computation, ui = i, uk = j
    integer  :: ui, uj              ! loop variables for Fock2e routine
    integer  :: uk, ul
    integer  :: um, un
    integer  :: uo, up
    integer  :: addi, addj          ! serial scalar for atomic operation
    real(dp) :: int
    !----------------------------------
    integer  :: contri              ! contr of atomi, shelli
    integer  :: atomi               ! which atom is the i^th component of |AOi>
    integer  :: shelli              ! which shell is the i^th component of |AOi>
    integer  :: Li                  ! angular quantum number of |AOi>
    integer  :: Mi                  ! magnetic quantum number of |AOi>
    !----------------------------------
    integer  :: contrj              ! contr of atomj, shellj
    integer  :: atomj               ! which atom is the j_th component of |AOj>
    integer  :: shellj              ! which shell is the j_th component of |AOj>
    integer  :: Lj                  ! angular quantum number of |AOj>
    integer  :: Mj                  ! magnetic quantum number of |AOj>
    !----------------------------------
    integer  :: contrk              ! contr of atoMk, shelLk
    integer  :: Lk                  ! angular quantum number of |AOk>
    integer  :: Mk                  ! magnetic quantum number of |AOk>
    !----------------------------------
    integer  :: contrl              ! contr of atoMl, shelLl
    integer  :: Ll                  ! angular quantum number of |AOl>
    integer  :: Ml                  ! magnetic quantum number of |AOl>
    !DIR$ ATTRIBUTES ALIGN:align_size :: expi, expj, expk, expl
    !DIR$ ATTRIBUTES ALIGN:align_size :: coei, coej, coek, coel
    real(dp) :: expi(20)            ! expo of |AOi>
    real(dp) :: expj(20)            ! expo of |AOj>
    real(dp) :: expk(20)            ! expo of |AOk>
    real(dp) :: expl(20)            ! expo of |AOl>
    real(dp) :: coei(20)            ! coefficient of |AOi>
    real(dp) :: coej(20)            ! coefficient of |AOj>
    real(dp) :: coek(20)            ! coefficient of |AOk>
    real(dp) :: coel(20)            ! coefficient of |AOl>
    real(dp) :: codi(3)             ! coordinate of center of |AOi>
    real(dp) :: codj(3)             ! coordinate of center of |AOj>
    real(dp) :: codk(3)             ! coordinate of center of |AOk>
    real(dp) :: codl(3)             ! coordinate of center of |AOl>
    real(dp) :: swint_mic
    !DIR$ ATTRIBUTES ALIGN:align_size :: HFcol_mic, HFexc_mic
    complex(dp) :: HFcol_mic(2*cbdm,2*cbdm)  ! micro 2e Fock matrix
    complex(dp) :: HFexc_mic(2*cbdm,2*cbdm)  ! micro 2e Fock matrix
    complex(dp),allocatable :: supp1(:,:), supp2(:,:), supp3(:,:)
    !DIR$ ATTRIBUTES ALIGN:align_size :: Fock2_assigned
    integer :: Fock2_assigned(2,8)         ! avoid duplicate assignment of Fock2
    integer :: assigned
    integer :: iassign
    logical :: carry
    if (ndschwarz) then
      ndschwarz = .false.
      write(60,'(a)') '  -- Schwarz screening of <ij||kl>'
      allocate(Fock2HFcol(2*fbdm,2*fbdm))
      allocate(Fock2HFexc(2*fbdm,2*fbdm))
      allocate(swint(cbdm,cbdm))
      swint = 0.0_dp
      do ui = 1, cbdm
        contri = atom_basis(mol(basis_inf(ui)%atom)%&
        basis_number+basis_inf(ui)%shell-1) % contr
        Li = basis_inf(ui) % L
        Mi = basis_inf(ui) % M
        expi(1:contri) = atom_basis(mol(basis_inf(ui)%atom)%&
        basis_number+basis_inf(ui)%shell-1) % expo(1:contri)
        coei(1:contri) = atom_basis(mol(basis_inf(ui)%atom)%&
        basis_number+basis_inf(ui)%shell-1) % Ncoe(1:contri,Mi)
        codi = mol(basis_inf(ui)%atom) % pos
        do uj = ui, cbdm
          contrj = atom_basis(mol(basis_inf(uj)%atom)%&
          basis_number+basis_inf(uj)%shell-1) % contr
          Lj = basis_inf(uj) % L
          Mj = basis_inf(uj) % M
          expj(1:contrj) = atom_basis(mol(basis_inf(uj)%atom)%&
          basis_number+basis_inf(uj)%shell-1) % expo(1:contrj)
          coej(1:contrj) = atom_basis(mol(basis_inf(uj)%atom)%&
          basis_number+basis_inf(uj)%shell-1) % Ncoe(1:contrj,Mj)
          codj = mol(basis_inf(uj)%atom) % pos
          if (Li >= Lj) then
            do um = 1, contri
              do un = 1, contrj
                do uo = 1, contri
                  do up = 1, contrj
                    swint_mic=coei(um)*coej(un)*coei(uo)*coej(up)*&
                    V_integral_2e(&
                    AO_fac(:,Li,Mi), AO_fac(:,Lj,Mj), &
                    AO_fac(:,Li,Mi), AO_fac(:,Lj,Mj), &
                    expi(um),expj(un),expi(uo),expj(up),&
                    codi,codj,codi,codj)
                    swint(ui,uj) = swint(ui,uj) + swint_mic
                  end do
                end do
              end do
            end do
          else
            do um = 1, contri
              do un = 1, contrj
                do uo = 1, contri
                  do up = 1, contrj
                    swint_mic=coei(um)*coej(un)*coei(uo)*coej(up)*&
                    V_integral_2e(&
                    AO_fac(:,Lj,Mj), AO_fac(:,Li,Mi), &
                    AO_fac(:,Lj,Mj), AO_fac(:,Li,Mi), &
                    expj(un),expi(um),expj(up),expi(uo),&
                    codj,codi,codj,codi)
                    swint(ui,uj) = swint(ui,uj) + swint_mic
                  end do
                end do
              end do
            end do
          end if
          swint(uj,ui) = swint(ui,uj)
        end do
      end do
      write(60,'(A,E10.2,A)') &
      '  -- complete! cutoff:',schwarz_VT,'; stored in swint'
    end if
    if (s_h) then
      allocate(supp1(2*sbdm,2*cbdm))
      call zgemm( 'N', 'C', 2*sbdm, 2*cbdm, 2*sbdm, c1, rho_m, &
      2*sbdm, exc2s, 2*cbdm, c0, supp1, 2*sbdm)
      deallocate(rho_m)
      allocate(rho_m(2*cbdm,2*cbdm))
      call zgemm( 'N', 'N', 2*cbdm, 2*cbdm, 2*sbdm, c1, exc2s, &
      2*cbdm, supp1, 2*sbdm, c0, rho_m, 2*cbdm)
      deallocate(supp1)
    end if
    allocate(supp1(2*cbdm,2*cbdm))
    allocate(supp3(2*cbdm,2*cbdm))
    supp1 = c0
    supp3 = c0
    ! parallel zone, running results consistent with serial
    !$omp parallel num_threads(threads) default(shared) private(i,j,ui,uj,uk,&
    !$omp& ul,um,un,uo,up,int,contri,Li,Mi,contrj,Lj,Mj,contrk,Lk,Mk,contrl,Ll,&
    !$omp& Ml,expi,expj,expk,expl,coei,coej,coek,coel,codi,codj,codk,codl,addi,&
    !$omp& addj,HFcol_mic,HFexc_mic,Fock2_assigned,assigned,iassign,carry)&
    !$omp& if(threads < nproc)
    HFcol_mic = c0
    HFexc_mic = c0
    !$omp do schedule(dynamic,5) collapse(2)
    ! (ij|kl)=(ji|kl)=(ij|lk)=(ji|lk)=(kl|ij)=(kl|ji)=(lk|ij)=(lk|ji)
    do i = cbdm, 1, -1
      do j = cbdm, 1, -1
        !----------------------------<ui>-------------------------------------
        ui = i
        contri = atom_basis(mol(basis_inf(ui)%atom)%&
        basis_number+basis_inf(ui)%shell-1) % contr
        Li = basis_inf(ui) % L
        Mi = basis_inf(ui) % M
        expi(1:contri) = atom_basis(mol(basis_inf(ui)%atom)%&
        basis_number+basis_inf(ui)%shell-1) % expo(1:contri)
        coei(1:contri) = atom_basis(mol(basis_inf(ui)%atom)%&
        basis_number+basis_inf(ui)%shell-1) % Ncoe(1:contri,Mi)
        codi = mol(basis_inf(ui)%atom) % pos
        !----------------------------<uk>-----------------------------------
        uk = j
        contrk = atom_basis(mol(basis_inf(uk)%atom)%&
        basis_number+basis_inf(uk)%shell-1) % contr
        Lk = basis_inf(uk) % L
        Mk = basis_inf(uk) % M
        expk(1:contrk) = atom_basis(mol(basis_inf(uk)%atom)%&
        basis_number+basis_inf(uk)%shell-1) % expo(1:contrk)
        coek(1:contrk) = atom_basis(mol(basis_inf(uk)%atom)%&
        basis_number+basis_inf(uk)%shell-1) % Ncoe(1:contrk,Mk)
        codk = mol(basis_inf(uk)%atom) % pos
        !----------------------------<uj>---------------------------------
        do uj = ui, 1, -1
          contrj = atom_basis(mol(basis_inf(uj)%atom)%&
          basis_number+basis_inf(uj)%shell-1) % contr
          Lj = basis_inf(uj) % L
          Mj = basis_inf(uj) % M
          expj(1:contrj) = atom_basis(mol(basis_inf(uj)%atom)%&
          basis_number+basis_inf(uj)%shell-1) % expo(1:contrj)
          coej(1:contrj) = atom_basis(mol(basis_inf(uj)%atom)%&
          basis_number+basis_inf(uj)%shell-1) % Ncoe(1:contrj,Mj)
          codj = mol(basis_inf(uj)%atom) % pos
          !----------------------------<ul>-------------------------------
          do ul = min(uk,uj+(ui*(ui-1)-uk*(uk-1))/2), 1, -1
            ! Schwarz screening, |<ij||kl>| <= dsqrt(<ij||ij>) * dsqrt(<kl||kl>)
            if (dsqrt(swint(ui,uj)*swint(uk,ul)) < schwarz_VT) cycle
            contrl = atom_basis(mol(basis_inf(ul)%atom)%&
            basis_number+basis_inf(ul)%shell-1) % contr
            Ll = basis_inf(ul) % L
            Ml = basis_inf(ul) % M
            expl(1:contrl) = atom_basis(mol(basis_inf(ul)%atom)%&
            basis_number+basis_inf(ul)%shell-1) % expo(1:contrl)
            coel(1:contrl) = atom_basis(mol(basis_inf(ul)%atom)%&
            basis_number+basis_inf(ul)%shell-1) % Ncoe(1:contrl,Ml)
            codl = mol(basis_inf(ul)%atom) % pos
            int = 0.0_dp
            !===========================<ij||kl>===============================
            if (Li >= Lj .and. Lk >= Ll) then
              do um = 1, contri
                do un = 1, contrj
                  do uo = 1, contrk
                    do up = 1, contrl
                      int = int + coei(um) * &
                      coej(un) * coek(uo) * &
                      coel(up) * V_integral_2e(&
                      AO_fac(:,Li,Mi), AO_fac(:,Lj,Mj),&
                      AO_fac(:,Lk,Mk), AO_fac(:,Ll,Ml),&
                      expi(um),expj(un),&
                      expk(uo),expl(up),&
                      codi,codj,codk,codl)
                    end do
                  end do
                end do
              end do
            else if (Li >= Lj .and. Lk < Ll) then
              do um = 1, contri
                do un = 1, contrj
                  do uo = 1, contrk
                    do up = 1, contrl
                      int = int + coei(um) * &
                      coej(un) * coek(uo) * &
                      coel(up) * V_integral_2e(&
                      AO_fac(:,Li,Mi), AO_fac(:,Lj,Mj),&
                      AO_fac(:,Ll,Ml), AO_fac(:,Lk,Mk),&
                      expi(um),expj(un),&
                      expl(up),expk(uo),&
                      codi,codj,codl,codk)
                    end do
                  end do
                end do
              end do
            else if (Li < Lj .and. Lk >= Ll) then
              do um = 1, contri
                do un = 1, contrj
                  do uo = 1, contrk
                    do up = 1, contrl
                      int = int + coei(um) * &
                      coej(un) * coek(uo) * &
                      coel(up) * V_integral_2e(&
                      AO_fac(:,Lj,Mj), AO_fac(:,Li,Mi),&
                      AO_fac(:,Lk,Mk), AO_fac(:,Ll,Ml),&
                      expj(un),expi(um),&
                      expk(uo),expl(up),&
                      codj,codi,codk,codl)
                    end do
                  end do
                end do
              end do
            else
              do um = 1, contri
                do un = 1, contrj
                  do uo = 1, contrk
                    do up = 1, contrl
                      int = int + coei(um) *&
                      coej(un) * coek(uo) * &
                      coel(up) * V_integral_2e(&
                      AO_fac(:,Lj,Mj), AO_fac(:,Li,Mi),&
                      AO_fac(:,Ll,Ml), AO_fac(:,Lk,Mk),&
                      expj(un),expi(um),&
                      expl(up),expk(uo),&
                      codj,codi,codl,codk)
                    end do
                  end do
                end do
              end do
            end if
            ! assign values to two-electron Fock matrix
            ! ref page 261 of Quantum Chemistry: Basic Principles and
            ! ab-initio Calculations, Volume 2 | 2nd Edition
            !------------------------<COULOMB>------------------------
            Fock2_assigned = 0
            carry = .true.
            assigned = 1
            Fock2_assigned(1,1) = ui
            Fock2_assigned(2,1) = uj
            HFcol_mic(ui,uj) = HFcol_mic(ui,uj) + int*rho_m(uk,ul)
            HFcol_mic(ui,uj) = HFcol_mic(ui,uj) + int*rho_m(cbdm+uk,cbdm+ul)
            HFcol_mic(cbdm+ui,cbdm+uj) = &
            HFcol_mic(cbdm+ui,cbdm+uj) + int*rho_m(uk,ul)
            HFcol_mic(cbdm+ui,cbdm+uj) = &
            HFcol_mic(cbdm+ui,cbdm+uj) + int*rho_m(cbdm+uk,cbdm+ul)
            if (uk /= ul) then
              HFcol_mic(ui,uj) = HFcol_mic(ui,uj) + int*rho_m(ul,uk)
              HFcol_mic(ui,uj) = HFcol_mic(ui,uj) + int*rho_m(cbdm+ul,cbdm+uk)
              HFcol_mic(cbdm+ui,cbdm+uj) = &
              HFcol_mic(cbdm+ui,cbdm+uj) + int*rho_m(ul,uk)
              HFcol_mic(cbdm+ui,cbdm+uj) = &
              HFcol_mic(cbdm+ui,cbdm+uj) + int*rho_m(cbdm+ul,cbdm+uk)
            end if
            do iassign = 1, assigned
              if (Fock2_assigned(1,iassign) == uj .and. &
              Fock2_assigned(2,iassign) == ui) then
                carry = .false.
                exit
              end if
            end do
            if (carry) then
              assigned = assigned + 1
              Fock2_assigned(1,assigned) = uj
              Fock2_assigned(2,assigned) = ui
              HFcol_mic(uj,ui) = HFcol_mic(uj,ui) + int*rho_m(uk,ul)
              HFcol_mic(uj,ui) = HFcol_mic(uj,ui) + int*rho_m(cbdm+uk,cbdm+ul)
              HFcol_mic(cbdm+uj,cbdm+ui) = &
              HFcol_mic(cbdm+uj,cbdm+ui) + int*rho_m(uk,ul)
              HFcol_mic(cbdm+uj,cbdm+ui) = &
              HFcol_mic(cbdm+uj,cbdm+ui) + int*rho_m(cbdm+uk,cbdm+ul)
              if (uk /= ul) then
                HFcol_mic(uj,ui) = HFcol_mic(uj,ui) + int*rho_m(ul,uk)
                HFcol_mic(uj,ui) = HFcol_mic(uj,ui) + int*rho_m(cbdm+ul,cbdm+uk)
                HFcol_mic(cbdm+uj,cbdm+ui) = &
                HFcol_mic(cbdm+uj,cbdm+ui) + int*rho_m(ul,uk)
                HFcol_mic(cbdm+uj,cbdm+ui) = &
                HFcol_mic(cbdm+uj,cbdm+ui) + int*rho_m(cbdm+ul,cbdm+uk)
              end if
            end if
            carry = .true.
            do iassign = 1, assigned
              if (Fock2_assigned(1,iassign) == uk .and. &
              Fock2_assigned(2,iassign) == ul) then
                carry = .false.
                exit
              end if
            end do
            if (carry) then
              assigned = assigned + 1
              Fock2_assigned(1,assigned) = uk
              Fock2_assigned(2,assigned) = ul
              HFcol_mic(uk,ul) = HFcol_mic(uk,ul)+int*rho_m(ui,uj)
              HFcol_mic(uk,ul) = HFcol_mic(uk,ul)+int*rho_m(cbdm+ui,cbdm+uj)
              HFcol_mic(cbdm+uk,cbdm+ul) = &
              HFcol_mic(cbdm+uk,cbdm+ul) + int*rho_m(ui,uj)
              HFcol_mic(cbdm+uk,cbdm+ul) = &
              HFcol_mic(cbdm+uk,cbdm+ul) + int*rho_m(cbdm+ui,cbdm+uj)
              if (ui /= uj) then
                HFcol_mic(uk,ul) = HFcol_mic(uk,ul) + int*rho_m(uj,ui)
                HFcol_mic(uk,ul) = HFcol_mic(uk,ul) + int*rho_m(cbdm+uj,cbdm+ui)
                HFcol_mic(cbdm+uk,cbdm+ul) = &
                HFcol_mic(cbdm+uk,cbdm+ul) + int*rho_m(uj,ui)
                HFcol_mic(cbdm+uk,cbdm+ul) = &
                HFcol_mic(cbdm+uk,cbdm+ul) + int*rho_m(cbdm+uj,cbdm+ui)
              end if
            end if
            carry = .true.
            do iassign = 1, assigned
              if (Fock2_assigned(1,iassign) == ul .and. &
              Fock2_assigned(2,iassign) == uk) then
                carry = .false.
                exit
              end if
            end do
            if (carry) then
              assigned = assigned + 1
              Fock2_assigned(1,assigned) = ul
              Fock2_assigned(2,assigned) = uk
              HFcol_mic(ul,uk) = HFcol_mic(ul,uk) + int*rho_m(ui,uj)
              HFcol_mic(ul,uk) = HFcol_mic(ul,uk) + int*rho_m(cbdm+ui,cbdm+uj)
              HFcol_mic(cbdm+ul,cbdm+uk) = &
              HFcol_mic(cbdm+ul,cbdm+uk) + int*rho_m(ui,uj)
              HFcol_mic(cbdm+ul,cbdm+uk) = &
              HFcol_mic(cbdm+ul,cbdm+uk) + int*rho_m(cbdm+ui,cbdm+uj)
              if (ui /= uj) then
                HFcol_mic(ul,uk) = HFcol_mic(ul,uk) + int*rho_m(uj,ui)
                HFcol_mic(ul,uk) = HFcol_mic(ul,uk) + int*rho_m(cbdm+uj,cbdm+ui)
                HFcol_mic(cbdm+ul,cbdm+uk) = &
                HFcol_mic(cbdm+ul,cbdm+uk) + int*rho_m(uj,ui)
                HFcol_mic(cbdm+ul,cbdm+uk) = &
                HFcol_mic(cbdm+ul,cbdm+uk) + int*rho_m(cbdm+uj,cbdm+ui)
              end if
            end if
            !------------------------<EXCHANGE int>------------------------
            Fock2_assigned = 0
              carry = .true.
              assigned = 1
              Fock2_assigned(1,1) = ui
              Fock2_assigned(2,1) = uk
              HFexc_mic(ui,uk) = HFexc_mic(ui,uk) - int*rho_m(uj,ul)
              HFexc_mic(cbdm+ui,cbdm+uk) = &
              HFexc_mic(cbdm+ui,cbdm+uk) - int*rho_m(cbdm+uj,cbdm+ul)
              if (ui == uk .and. ul /= uj) then
                HFexc_mic(ui,uk) = HFexc_mic(ui,uk) - int*rho_m(ul,uj)
                HFexc_mic(cbdm+ui,cbdm+uk) = &
                HFexc_mic(cbdm+ui,cbdm+uk) - int*rho_m(cbdm+ul,cbdm+uj)
              end if
            do iassign = 1, assigned
              if (Fock2_assigned(1,iassign) == uk .and. &
              Fock2_assigned(2,iassign) == ui) then
                carry = .false.
                exit
              end if
            end do
            if (carry) then
              assigned = assigned + 1
              Fock2_assigned(1,assigned) = uk
              Fock2_assigned(2,assigned) = ui
              HFexc_mic(uk,ui) = HFexc_mic(uk,ui) - int*rho_m(ul,uj)
              HFexc_mic(cbdm+uk,cbdm+ui) = &
              HFexc_mic(cbdm+uk,cbdm+ui) - int*rho_m(cbdm+ul,cbdm+uj)
              if (uk == ui .and. ul /= uj) then
                HFexc_mic(uk,ui) = HFexc_mic(uk,ui) - int*rho_m(uj,ul)
                HFexc_mic(cbdm+uk,cbdm+ui) = &
                HFexc_mic(cbdm+uk,cbdm+ui) - int*rho_m(cbdm+uj,cbdm+ul)
              end if
            end if
            carry = .true.
            do iassign = 1, assigned
              if (Fock2_assigned(1,iassign) == ui .and. &
              Fock2_assigned(2,iassign) == ul) then
                carry = .false.
                exit
              end if
            end do
            if (carry) then
              assigned = assigned + 1
              Fock2_assigned(1,assigned) = ui
              Fock2_assigned(2,assigned) = ul
              HFexc_mic(ui,ul) = HFexc_mic(ui,ul) - int*rho_m(uj,uk)
              HFexc_mic(cbdm+ui,cbdm+ul) = &
              HFexc_mic(cbdm+ui,cbdm+ul) - int*rho_m(cbdm+uj,cbdm+uk)
              if (ui == ul .and. uk /= uj) then
                HFexc_mic(ui,ul) = HFexc_mic(ui,ul) - int*rho_m(uk,uj)
                HFexc_mic(cbdm+ui,cbdm+ul) = &
                HFexc_mic(cbdm+ui,cbdm+ul) - int*rho_m(cbdm+uk,cbdm+uj)
              end if
            end if
            carry = .true.
            do iassign = 1, assigned
              if (Fock2_assigned(1,iassign) == ul .and. &
               Fock2_assigned(2,iassign) == ui) then
                carry = .false.
                exit
              end if
            end do
            if (carry) then
              assigned = assigned + 1
              Fock2_assigned(1,assigned) = ul
              Fock2_assigned(2,assigned) = ui
              HFexc_mic(ul,ui) = HFexc_mic(ul,ui) - int*rho_m(uk,uj)
              HFexc_mic(cbdm+ul,cbdm+ui) = &
              HFexc_mic(cbdm+ul,cbdm+ui) - int*rho_m(cbdm+uk,cbdm+uj)
              if (ul == ui .and. uk /= uj) then
                HFexc_mic(ul,ui) = HFexc_mic(ul,ui) - int*rho_m(uj,uk)
                HFexc_mic(cbdm+ul,cbdm+ui) = &
                HFexc_mic(cbdm+ul,cbdm+ui) - int*rho_m(cbdm+uj,cbdm+uk)
              end if
            end if
            carry = .true.
            do iassign = 1, assigned
              if (Fock2_assigned(1,iassign) == uj .and. &
              Fock2_assigned(2,iassign) == uk) then
                carry = .false.
                exit
              end if
            end do
            if (carry) then
              assigned = assigned + 1
              Fock2_assigned(1,assigned) = uj
              Fock2_assigned(2,assigned) = uk
              HFexc_mic(uj,uk) = HFexc_mic(uj,uk) - int*rho_m(ui,ul)
              HFexc_mic(cbdm+uj,cbdm+uk) = &
              HFexc_mic(cbdm+uj,cbdm+uk) - int*rho_m(cbdm+ui,cbdm+ul)
              if (uj == uk .and. ul /= ui) then
                HFexc_mic(uj,uk) = HFexc_mic(uj,uk) - int*rho_m(ul,ui)
                HFexc_mic(cbdm+uj,cbdm+uk) = &
                HFexc_mic(cbdm+uj,cbdm+uk) - int*rho_m(cbdm+ul,cbdm+ui)
              end if
            end if
            carry = .true.
            do iassign = 1, assigned
              if (Fock2_assigned(1,iassign) == uk .and. &
              Fock2_assigned(2,iassign) == uj) then
                carry = .false.
                exit
              end if
            end do
            if (carry) then
              assigned = assigned + 1
              Fock2_assigned(1,assigned) = uk
              Fock2_assigned(2,assigned) = uj
              HFexc_mic(uk,uj) = HFexc_mic(uk,uj) - int*rho_m(ul,ui)
              HFexc_mic(cbdm+uk,cbdm+uj) = &
              HFexc_mic(cbdm+uk,cbdm+uj) - int*rho_m(cbdm+ul,cbdm+ui)
              if (uk == uj .and. ui /= ul) then
                HFexc_mic(uk,uj) = HFexc_mic(uk,uj) - int*rho_m(ui,ul)
                HFexc_mic(cbdm+uk,cbdm+uj) = &
                HFexc_mic(cbdm+uk,cbdm+uj) - int*rho_m(cbdm+ui,cbdm+ul)
              end if
            end if
            carry = .true.
            do iassign = 1, assigned
              if (Fock2_assigned(1,iassign) == uj .and. &
              Fock2_assigned(2,iassign) == ul) then
                carry = .false.
                exit
              end if
            end do
            if (carry) then
              assigned = assigned + 1
              Fock2_assigned(1,assigned) = uj
              Fock2_assigned(2,assigned) = ul
              HFexc_mic(uj,ul) = HFexc_mic(uj,ul) - int*rho_m(ui,uk)
              HFexc_mic(cbdm+uj,cbdm+ul) = &
              HFexc_mic(cbdm+uj,cbdm+ul) - int*rho_m(cbdm+ui,cbdm+uk)
              if (uj == ul .and. uk /= ui) then
                HFexc_mic(uj,ul) = HFexc_mic(uj,ul) - int*rho_m(uk,ui)
                HFexc_mic(cbdm+uj,cbdm+ul) = &
                HFexc_mic(cbdm+uj,cbdm+ul) - int*rho_m(cbdm+uk,cbdm+ui)
              end if
            end if
            carry = .true.
            do iassign = 1, assigned
              if (Fock2_assigned(1,iassign) == ul .and. &
              Fock2_assigned(2,iassign) == uj) then
                carry = .false.
                exit
              end if
            end do
            if (carry) then
              assigned = assigned + 1
              Fock2_assigned(1,assigned) = ul
              Fock2_assigned(2,assigned) = uj
              HFexc_mic(ul,uj) = HFexc_mic(ul,uj) - int*rho_m(uk,ui)
              HFexc_mic(cbdm+ul,cbdm+uj) = &
              HFexc_mic(cbdm+ul,cbdm+uj) - int*rho_m(cbdm+uk,cbdm+ui)
              if (ul == uj .and. ui /= uk) then
                HFexc_mic(ul,uj) = HFexc_mic(ul,uj) - int*rho_m(ui,uk)
                HFexc_mic(cbdm+ul,cbdm+uj) = &
                HFexc_mic(cbdm+ul,cbdm+uj) - int*rho_m(cbdm+ui,cbdm+uk)
              end if
            end if
          end do
        end do
      end do
    end do
    !$omp end do
    !-------------<thread sync>-------------
    !$omp critical
    supp1 = supp1 + HFcol_mic
    supp3 = supp3 + HFexc_mic
    !$omp end critical
    !$omp end parallel
    ! assign Kohn-Sham matrices
    if (fx_id /= -1) then
      if (allocated(Fock2KSexc)) deallocate(Fock2KSexc)
      allocate(Fock2KSexc(2*cbdm, 2*cbdm))
      if (allocated(oper5)) deallocate(oper5)
      if (fc_id /= -1) then
        if (allocated(Fock2KScor)) deallocate(Fock2KScor)
        allocate(Fock2KScor(2*cbdm, 2*cbdm))
        call basis2grid_Becke(rho_m, KSexc, KScor, Fock2KSexc, Fock2KScor)
        if (s_h) then
          allocate(oper5(2*cbdm,2*sbdm))
          ! exchange
          call zgemm( 'N', 'N', 2*cbdm, 2*sbdm, 2*cbdm, c1, &
          Fock2KSexc, 2*cbdm, exc2s, 2*cbdm, c0, oper5, 2*cbdm)
          deallocate(Fock2KSexc)
          allocate(Fock2KSexc(2*sbdm,2*sbdm))
          call zgemm( 'C', 'N', 2*sbdm, 2*sbdm, 2*cbdm, c1, &
          exc2s, 2*cbdm, oper5, 2*cbdm, c0, Fock2KSexc, 2*sbdm)
          ! correlation
          call zgemm( 'N', 'N', 2*cbdm, 2*sbdm, 2*cbdm, c1, &
          Fock2KScor, 2*cbdm, exc2s, 2*cbdm, c0, oper5, 2*cbdm)
          deallocate(Fock2KScor)
          allocate(Fock2KScor(2*sbdm,2*sbdm))
          call zgemm( 'C', 'N', 2*sbdm, 2*sbdm, 2*cbdm, c1, &
          exc2s, 2*cbdm, oper5, 2*cbdm, c0, Fock2KScor, 2*sbdm)
          deallocate(oper5)
        end if
        allocate(oper5(2*fbdm,2*sbdm))
        ! exchange
        call zgemm( 'C', 'N', 2*fbdm, 2*sbdm, 2*sbdm, c1, &
        exXm, 2*sbdm, Fock2KSexc, 2*sbdm, c0, oper5, 2*fbdm)
        deallocate(Fock2KSexc)
        allocate(Fock2KSexc(2*fbdm, 2*fbdm))
        call zgemm( 'N', 'N', 2*fbdm, 2*fbdm, 2*sbdm, c1, &
        oper5, 2*fbdm, exXm, 2*sbdm, c0, Fock2KSexc, 2*fbdm)
        ! correlation
        call zgemm( 'C', 'N', 2*fbdm, 2*sbdm, 2*sbdm, c1, &
        exXm, 2*sbdm, Fock2KScor, 2*sbdm, c0, oper5, 2*fbdm)
        deallocate(Fock2KScor)
        allocate(Fock2KScor(2*fbdm, 2*fbdm))
        call zgemm( 'N', 'N', 2*fbdm, 2*fbdm, 2*sbdm, c1, &
        oper5, 2*fbdm, exXm, 2*sbdm, c0, Fock2KScor, 2*fbdm)
      else
        call basis2grid_Becke(rho_m=rho_m, ex=KSexc, Fockx=Fock2KSexc)
        if (s_h) then
          allocate(oper5(2*cbdm,2*sbdm))
          ! exchange-correlation
          call zgemm( 'N', 'N', 2*cbdm, 2*sbdm, 2*cbdm, c1, &
          Fock2KSexc, 2*cbdm, exc2s, 2*cbdm, c0, oper5, 2*cbdm)
          deallocate(Fock2KSexc)
          allocate(Fock2KSexc(2*sbdm,2*sbdm))
          call zgemm( 'C', 'N', 2*sbdm, 2*sbdm, 2*cbdm, c1, &
          exc2s, 2*cbdm, oper5, 2*cbdm, c0, Fock2KSexc, 2*sbdm)
        end if
        allocate(oper5(2*fbdm,2*sbdm))
        ! exchange-correlation
        call zgemm( 'C', 'N', 2*fbdm, 2*sbdm, 2*sbdm, c1, &
        exXm, 2*sbdm, Fock2KSexc, 2*sbdm, c0, oper5, 2*fbdm)
        deallocate(Fock2KSexc)
        allocate(Fock2KSexc(2*fbdm, 2*fbdm))
        call zgemm( 'N', 'N', 2*fbdm, 2*fbdm, 2*sbdm, c1, &
        oper5, 2*fbdm, exXm, 2*sbdm, c0, Fock2KSexc, 2*fbdm)
      end if
    end if



    if (s_h) then
      ! transform to spherical-harmonic basis
      allocate(supp2(2*cbdm,2*sbdm))
      call zgemm( 'N', 'N', 2*cbdm, 2*sbdm, 2*cbdm, c1, &
      rho_m, 2*cbdm, exc2s, 2*cbdm, c0, supp2, 2*cbdm)
      deallocate(rho_m)
      allocate(rho_m(2*sbdm,2*sbdm))
      call zgemm( 'C', 'N', 2*sbdm, 2*sbdm, 2*cbdm, c1, &
      exc2s, 2*cbdm, supp2, 2*cbdm, c0, rho_m, 2*sbdm)
    
      ! Coulomb
      call zgemm( 'N', 'N', 2*cbdm, 2*sbdm, 2*cbdm, c1, &
      supp1, 2*cbdm, exc2s, 2*cbdm, c0, supp2, 2*cbdm)
      deallocate(supp1)
      allocate(supp1(2*sbdm,2*sbdm))
      call zgemm( 'C', 'N', 2*sbdm, 2*sbdm, 2*cbdm, c1, &
      exc2s, 2*cbdm, supp2, 2*cbdm, c0, supp1, 2*sbdm)
      ! Exchange
      call zgemm( 'N', 'N', 2*cbdm, 2*sbdm, 2*cbdm, c1, &
      supp3, 2*cbdm, exc2s, 2*cbdm, c0, supp2, 2*cbdm)
      deallocate(supp3)
      allocate(supp3(2*sbdm,2*sbdm))
      call zgemm( 'C', 'N', 2*sbdm, 2*sbdm, 2*cbdm, c1, &
      exc2s, 2*cbdm, supp2, 2*cbdm, c0, supp3, 2*sbdm)
      deallocate(supp2)
    end if
    
    ! orth to Fock2
    allocate(supp2(2*fbdm,2*sbdm))
    ! Coulomb
    call zgemm( 'C', 'N', 2*fbdm, 2*sbdm, 2*sbdm, c1, &
    exXm, 2*sbdm, supp1, 2*sbdm, c0, supp2, 2*fbdm)
    call zgemm( 'N', 'N', 2*fbdm, 2*fbdm, 2*sbdm, c1, &
    supp2, 2*fbdm, exXm, 2*sbdm, c0, Fock2HFcol, 2*fbdm)
    ! Exchange
    call zgemm( 'C', 'N', 2*fbdm, 2*sbdm, 2*sbdm, c1, &
    exXm, 2*sbdm, supp3, 2*sbdm, c0, supp2, 2*fbdm)
    call zgemm( 'N', 'N', 2*fbdm, 2*fbdm, 2*sbdm, c1, &
    supp2, 2*fbdm, exXm, 2*sbdm, c0, Fock2HFexc, 2*fbdm)
    deallocate(supp1)
    deallocate(supp2)
    deallocate(supp3)
  end subroutine Fock2e

!------------------------------------------------------------
!> transfer alternating zero matrix to block matrix
  subroutine atnz2block(atnz, dm)
    implicit none
    integer, intent(in) :: dm
    integer :: aloop_i, aloop_j
    complex(dp) :: atnz(dm,dm), aoper(dm,dm)
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
!> calculate <S**2> based on oper3 generated by SCF
!!
!! divide all occupied orbitals into two sets of alpha, beta orbitals with
!!
!! their original coefficients, then use UHF method to solve Lowdin <S**2>.
  subroutine calc_S2HF()
    implicit none
    integer cloop_i, cloop_j, cloop_k   ! loop variables for calc_S2HF
    complex(dp) :: alal, albe, beal, bebe  ! components of pair density matrix
    complex(dp) :: suppa, suppb
    ! oper3 and i_j are orthogonally normalized
    totalpha = 0.0_dp
    totbeta = 0.0_dp
    do cloop_i = 1, electron_count
      do cloop_k = 1, fbdm
        totalpha = totalpha + &
        real(conjg(oper3(cloop_k,cloop_i))*oper3(cloop_k,cloop_i),dp)
      end do
      do cloop_k = fbdm+1, 2*fbdm
        totbeta = totbeta + &
        real(conjg(oper3(cloop_k,cloop_i))*oper3(cloop_k,cloop_i),dp)
      end do
    end do
    alal = totalpha * (totalpha - 1.0_dp)
    bebe = totbeta * (totbeta - 1.0_dp)
    albe = c0
    do cloop_i = 1, electron_count
      do cloop_j = 1, electron_count
        suppa = c0
        suppb = c0
        do cloop_k = 1, fbdm
          suppa = suppa + &
          conjg(oper3(cloop_k,cloop_i))*oper3(fbdm+cloop_k,cloop_j)
          suppb = suppb + &
          conjg(oper3(fbdm+cloop_k,cloop_j))*oper3(cloop_k,cloop_i)
        end do
        albe = albe + suppa*suppb
      end do
    end do
    beal = c0
    do cloop_i = 1, electron_count
      do cloop_j = 1, electron_count
        suppa = c0
        suppb = c0
        do cloop_k = fbdm+1, 2*fbdm
          suppa = suppa + &
          conjg(oper3(cloop_k,cloop_i))*oper3(cloop_k-fbdm,cloop_j)
          suppb = suppb + &
          conjg(oper3(cloop_k-fbdm,cloop_j))*oper3(cloop_k,cloop_i)
        end do
        beal = beal + suppa*suppb
      end do
    end do
    S__2 = -real(electron_count*(electron_count-4),dp)/4.0_dp + &
    real(alal/2.0_dp + bebe/2.0_dp) - real(albe/2.0_dp + beal/2.0_dp)
  end subroutine calc_S2HF
  
!------------------------------------------------------------
!> calculate <S**2> for orbitals based on oper3 generated by SCF
  subroutine calc_S2HForb(orbnum)
    implicit none
    integer,intent(in) :: orbnum
    integer :: zloop_i     !loop variables for calc_S2HForb
    real(dp) :: suppa, suppb
    totalphaorb = 0.0_dp
    do zloop_i = 1, fbdm
      totalphaorb = totalphaorb + &
      real(conjg(oper3(zloop_i,orbnum))*oper3(zloop_i,orbnum),dp)
    end do
    totbetaorb = 0.0_dp
    do zloop_i = fbdm+1, 2*fbdm
      totbetaorb = totbetaorb + &
      real(conjg(oper3(zloop_i,orbnum))*oper3(zloop_i,orbnum),dp)
    end do
    suppa = c0
    suppb = c0
    do zloop_i = fbdm+1, 2*fbdm
      suppa = suppa + conjg(oper3(zloop_i,orbnum))*oper3(zloop_i-fbdm,orbnum)
      suppb = suppb + conjg(oper3(zloop_i-fbdm,orbnum))*oper3(zloop_i,orbnum)
    end do
    S__2orb = 3.0/4.0 + totalphaorb*(totalphaorb-1.0)/2.0 + &
    totbetaorb*(totbetaorb-1.0)/2.0 - real(suppa*suppb)
    Szorb = (totalphaorb-totbetaorb)/2.0_dp
  end subroutine calc_S2HForb
  
!-----------------------------------------------------------------------
!> dump mol orbital information to .molden.input file
  subroutine dump_molden()
    implicit none
    integer :: channel, dmi, dmj, dmk
    if (.not. allocated(AO2MO)) &
    call terminate('dump mol orbital failed, AO2MO is empty')
    
    ! molden file contains the real part of mol orbital
    open(newunit=channel, file=trim(address_job)//'-real.molden.input', &
    status='replace', action='write', iostat=ios)
    if (ios /= 0) call terminate('dump .molden.input failed')
    write(channel, '(A)') '[Molden Format]'
    write(channel, '(A)') '[Title]'
    write(channel, '(A)') &
    'The real part of molecular orbitals of job '//trim(address_job)
    write(channel, *)
    ! mol geometry
    write(channel, '(A)') '[Atoms] AU'
    do dmi = 1, atom_count
      write(channel, '(A,I3,I3,F13.7,F13.7,F13.7)') &
      element_list(mol(dmi)%atom_number), dmi, &
      mol(dmi)%atom_number, mol(dmi)%pos(1), &
      mol(dmi)%pos(2), mol(dmi)%pos(3)
    end do
    ! basis of each atom
    write(channel, '(A)') '[GTO]'
    do dmi = 1, atom_count
      write(channel, '(I3,I2)') dmi, 0
      do dmj = 0, shell_in_element(mol(dmi) % atom_number)-1
        if (atom_basis(mol(dmi)%basis_number+dmj)%&
        L == 0) then
          write(channel, '(A2,I2,A)') 's', &
          atom_basis(mol(dmi)%basis_number+dmj)%contr,' 1.0'
        else if (atom_basis(mol(dmi)%basis_number+dmj)%&
        L == 1) then
          write(channel, '(A2,I2,A)') 'p', &
          atom_basis(mol(dmi)%basis_number+dmj)%contr,' 1.0'
        else if (atom_basis(mol(dmi)%basis_number+dmj)%&
        L == 2) then
          write(channel, '(A2,I2,A)') 'd', &
          atom_basis(mol(dmi)%basis_number+dmj)%contr,' 1.0'
        else if (atom_basis(mol(dmi)%basis_number+dmj)%&
        L == 3) then
          write(channel, '(A2,I2,A)') 'f', &
          atom_basis(mol(dmi)%basis_number+dmj)%contr,' 1.0'
        else if (atom_basis(mol(dmi)%basis_number+dmj)%&
        L == 4) then
          write(channel, '(A2,I2,A)') 'g', &
          atom_basis(mol(dmi)%basis_number+dmj)%contr,' 1.0'
        end if
        do dmk = 1, atom_basis(mol(dmi)%basis_number+dmj)%contr
          write(channel, '(F13.7,F13.7)') &
          atom_basis(mol(dmi)%basis_number+dmj)%expo(dmk), &
          atom_basis(mol(dmi)%basis_number+dmj)%coe(dmk)
        end do
      end do
      write(channel, *)
    end do
    if (s_h) then
      write(channel, '(A)') '[5D]'
      write(channel, '(A)') '[7F]'
      write(channel, '(A)') '[9G]'
    else
      write(channel, '(A)') '[6D]'
      write(channel, '(A)') '[10F]'
      write(channel, '(A)') '[15G]'
    end if
    ! MO coefficient
    write(channel, '(A)') '[MO]'
    do dmi = 1, fbdm
      write(channel, '(A4,E23.14)') 'Ene=', orbE(dmi)
      write(channel, '(A)') 'Spin= Alpha'
      call calc_S2HForb(dmi)
      if (dmi <= electron_count) then
        write(channel, '(A7,F12.6)') 'Occup=', 0.5 + Szorb
      else
        write(channel, '(A)') 'Occup= 0.00'
      end if
      do dmj = 1, sbdm
        write(channel, '(I4,F20.12)') dmj, real(AO2MO(dmj, dmi))
      end do
    end do
    do dmi = 1, fbdm
      write(channel, '(A4,E23.14)') 'Ene=', orbE(dmi)
      write(channel, '(A)') 'Spin= Beta'
      call calc_S2HForb(dmi)
      if (dmi <= electron_count) then
        write(channel, '(A7,F12.6)') 'Occup=', 0.5 - Szorb
      else
        write(channel, '(A)') 'Occup= 0.00'
      end if
      do dmj = 1, sbdm
        write(channel, '(I4,F20.12)') dmj, real(AO2MO(sbdm+dmj, dmi))
      end do
    end do
    close(channel)
    
    if (DKH_order /= 0) then
      ! molden file contains the maginary part of mol orbital
      open(newunit=channel, file=trim(address_job)//'-img.molden.input', &
      status='replace', action='write', iostat=ios)
      if (ios /= 0) call terminate('dump .molden.input failed')
      write(channel, '(A)') '[Molden Format]'
      write(channel, '(A)') '[Title]'
      write(channel, '(A)') &
      'The imaginary part of molecular orbitals of job '//trim(address_job)
      write(channel, *)
      ! mol geometry
      write(channel, '(A)') '[Atoms] AU'
      do dmi = 1, atom_count
        write(channel, '(A,I3,I3,F13.7,F13.7,F13.7)') &
        element_list(mol(dmi)%atom_number), dmi, &
        mol(dmi)%atom_number, mol(dmi)%pos(1), &
        mol(dmi)%pos(2), mol(dmi)%pos(3)
      end do
      ! basis of each atom
      write(channel, '(A)') '[GTO]'
      do dmi = 1, atom_count
        write(channel, '(I3,I2)') dmi, 0
        do dmj = 0, shell_in_element(mol(dmi) % atom_number)-1
          if (atom_basis(mol(dmi)%basis_number+dmj)%&
          L == 0) then
            write(channel, '(A2,I2,A)') 's', &
            atom_basis(mol(dmi)%basis_number+dmj)%contr,' 1.0'
          else if (atom_basis(mol(dmi)%basis_number+dmj)%&
          L == 1) then
            write(channel, '(A2,I2,A)') 'p', &
            atom_basis(mol(dmi)%basis_number+dmj)%contr,' 1.0'
          else if (atom_basis(mol(dmi)%basis_number+dmj)%&
          L == 2) then
            write(channel, '(A2,I2,A)') 'd', &
            atom_basis(mol(dmi)%basis_number+dmj)%contr,' 1.0'
          else if (atom_basis(mol(dmi)%basis_number+dmj)%&
          L == 3) then
            write(channel, '(A2,I2,A)') 'f', &
            atom_basis(mol(dmi)%basis_number+dmj)%contr,' 1.0'
          else if (atom_basis(mol(dmi)%basis_number+dmj)%&
          L == 4) then
            write(channel, '(A2,I2,A)') 'g', &
            atom_basis(mol(dmi)%basis_number+dmj)%contr,' 1.0'
          end if
          do dmk = 1, atom_basis(mol(dmi)%basis_number+dmj)%contr
            write(channel, '(F13.7,F13.7)') &
            atom_basis(mol(dmi)%basis_number+dmj)%expo(dmk), &
            atom_basis(mol(dmi)%basis_number+dmj)%coe(dmk)
          end do
        end do
        write(channel, *)
      end do
      if (s_h) then
        write(channel, '(A)') '[5D]'
        write(channel, '(A)') '[7F]'
        write(channel, '(A)') '[9G]'
      else
        write(channel, '(A)') '[6D]'
        write(channel, '(A)') '[10F]'
        write(channel, '(A)') '[15G]'
      end if
      ! MO coefficient
      write(channel, '(A)') '[MO]'
      do dmi = 1, fbdm
        write(channel, '(A4,E23.14)') 'Ene=', orbE(dmi)
        write(channel, '(A)') 'Spin= Alpha'
        call calc_S2HForb(dmi)
        if (dmi <= electron_count) then
          write(channel, '(A7,F12.6)') 'Occup=', 0.5 + Szorb
        else
          write(channel, '(A)') 'Occup= 0.000'
        end if
        do dmj = 1, sbdm
          write(channel, '(I4,F20.12)') dmj, aimag(AO2MO(dmj, dmi))
        end do
      end do
      do dmi = 1, fbdm
        write(channel, '(A4,E23.14)') 'Ene=', orbE(dmi)
        write(channel, '(A)') 'Spin= Beta'
        call calc_S2HForb(dmi)
        if (dmi <= electron_count) then
          write(channel, '(A7,F12.6)') 'Occup=', 0.5 - Szorb
        else
          write(channel, '(A)') 'Occup= 0.000'
        end if
        do dmj = 1, sbdm
          write(channel, '(I4,F20.12)') dmj, aimag(AO2MO(sbdm+dmj, dmi))
        end do
      end do
      close(channel)
    end if
  end subroutine dump_molden
  
end module SCF