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
  
  ! density matrices and mol orbital coefficients
  integer                 :: Nalpha, Nbeta  ! number of alpha and beta elctron
  ! spinor MO coefficients, in order of (AO1,0), (0,AO1), ... (AOn,0), (0,AOn)
  !DIR$ ATTRIBUTES ALIGN:align_size :: AO2MO, AO2MOalpha, AO2MObeta
  !DIR$ ATTRIBUTES ALIGN:align_size :: rho_m, AOsupp, Fock
  complex(dp),allocatable :: AO2MO(:,:)     ! MO coeff
  integer,allocatable     :: occindex(:)    ! occupid orbital number of MO
  real(dp),allocatable    :: m_AO2MO_a(:,:) ! alpha MO coeff read form molden
  real(dp),allocatable    :: t_AO2MO_a(:,:) ! transfered alpha MO coeff
  real(dp),allocatable    :: m_AO2MO_b(:,:) ! beta MO coeff read form molden
  real(dp),allocatable    :: t_AO2MO_b(:,:) ! transfered beta MO coeff
  complex(dp),allocatable :: rho_m(:,:)     ! density matrix, complex Hermitian
  real(dp)                :: RMSDP, maxDP
  real(dp),allocatable    :: AOsupp(:,:)
  complex(dp),allocatable :: Fock(:,:)      ! Fock matrix
  
  logical                 :: ini_rho =.true.! initial density matrix loaded
  
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
  real(dp)                :: scf_kappa      ! deviation parameter from TRS

  !DIR$ ATTRIBUTES ALIGN:align_size :: rho_pre, rho_pre_pre, rho_history, Rsd
  ! damping & DIIS(AX=B)
  ! privious rho_m of current iteration
  complex(dp),allocatable :: rho_pre(:,:)
  ! privious rho_m of privious iteration
  complex(dp),allocatable :: rho_pre_pre(:,:)
  complex(dp),allocatable :: rho_history(:,:,:)! coeff of subsp iteration
  complex(dp),allocatable :: Rsd(:,:,:)        ! residuals of rho_m
  complex(dp),allocatable :: DIISmat(:,:)      ! A
  real(dp)                :: damp_coe          ! damp coeff of direct/DIIS SCF
  real(dp)                :: dE_pre
  logical                 :: forward = .true.
  
  contains

!------------------------------------------------------------
!> HF/KS-SCF procedure for NR/DKH2 Hamiltonian
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
    write(60,'(A,I3.3)')  '  -- maxiter   = ',maxiter
    write(60,'(A,E10.3)') '  -- conv_tol  =',conver_tol
    write(60,'(A,F6.3)')  '  -- damp      =',damp
    write(60,'(A,E10.3)') '  -- cutdamp   =',cutdamp
    write(60,'(A,I3.3)')  '  -- nodiis    = ',nodiis
    write(60,'(A,I3.3)')  '  -- subsp     = ',subsp
    write(60,'(A,F6.3)')  '  -- diisdamp  =',diisdamp
    write(60,'(A,E10.3)') '  -- cutdiis   =',cutdiis
    write(60,'(A)') '  ----------<SCF>----------'
    allocate(Fock(2*fbdm,2*fbdm))
    allocate(orbE(2*fbdm))
    allocate(oper6(2*fbdm,2*fbdm))
    allocate(oper3(2*fbdm,2*fbdm))
    allocate(oper4(2*fbdm,2*fbdm))
    
    ! DIIS
    ! AO2MO(new) = ��(i,subsp) DIIScoe(i)*(rho_history(i)+diisdamp*Rsd(i))
    allocate(Rsd(subsp,2*sbdm,2*sbdm))
    allocate(DIISmat(subsp+1,subsp+1))
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
              "(I2)",iostat = ios) new_threads
              if (new_threads /= threads) then
                threads = new_threads
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
        write(60,'(A)') '  load density matrix'
        call assign_rho()
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

      ! diagonalization of Fock matrix
      write(60,'(A)') '  diagonalization of Fock matrix'
      call diag(Fock, 2*fbdm, oper3, orbE)
      write(60,'(A,I5,A)') '  complete!',2*fbdm,' eigenvectors found'
      ! de-orthogonalization
      write(60,'(A)') '  generate and dump AO2MO to .ao2mo file'
      call matmul('N', 'N', exs2f, oper3, AO2MO)
      call dump_matrix('ao2mo', AO2MO, 2*sbdm, 2*fbdm)
      ! construct new density matrix
      write(60,'(A)') '  construct new density matrix'
      call assign_rho()
      write(60,'(A)') '  complete! stored in rho_m'


      ! frontier orbital energy
      write(60,'(A)') '  frontier orbital energy (A.U.)'
      call calc_S2HForb(occindex(electron_count))
      write(60,'(A,I3,F12.6,A,F6.3)') &
      '  -- HOMO ', occindex(electron_count), orbE(occindex(electron_count)), &
      ' <Sz> = ',Szorb
      call calc_S2HForb(occindex(electron_count)+1)
      write(60,'(A,I3,F12.6,A,F6.3)') &
      '  -- LUMO ',occindex(electron_count)+1,orbE(occindex(electron_count)+1),&
      ' <Sz> = ',Szorb
      write(60,'(A,F12.6)') '  -- gap ',&
      orbE(occindex(electron_count)+1)-orbE(occindex(electron_count))
      ! energy components calculation
      write(60,'(A)') '  calculate energy components (A.U.)'
      ! HF Coulomb energy
      HFCol = 0.0_dp
      call matmul('C', 'N', oper3, Fock2HFcol, oper6)
      call matmul('N', 'N', oper6, oper3, oper4)
      do loop_j = 1, electron_count
        HFCol = HFCol + real(oper4(occindex(loop_j),occindex(loop_j)))
      end do
      write(60,'(A,F12.6)') '  -- HF Coulomb energy                    ', HFCol

      ! HF exchange energy
      HFexc = 0.0_dp
      call matmul('C', 'N', oper3, Fock2HFexc, oper6)
      call matmul('N', 'N', oper6, oper3, oper4)
      do loop_j = 1, electron_count
        HFexc = HFexc + real(oper4(occindex(loop_j),occindex(loop_j)))
      end do
      write(60,'(A,F12.6)') &
      '  -- HF exchange energy                   ', x_HF*HFexc
      write(60,'(A,F12.6)') &
      '  -- KS exchange energy                   ', (1-x_HF)*KSexc
      write(60,'(A,F12.6)') '  -- KS correlation energy                ', KScor

      ! core energy
      Ecore = 0.0_dp
      call matmul('C', 'N', oper3, Fock1, oper6)
      call matmul('N', 'N', oper6, oper3, oper4)
      do loop_j = 1, electron_count
        Ecore = Ecore + real(oper4(occindex(loop_j),occindex(loop_j)))
      end do
      write(60,'(A,F12.6)') '  -- core energy                          ', Ecore

      ! kinetic energy
      T = 0.0_dp
      call matmul('C', 'N', oper3, exi_T_j, oper6)
      call matmul('N', 'N', oper6, oper3, oper4)
      do loop_j = 1, electron_count
        T = T + real(oper4(occindex(loop_j),occindex(loop_j)))
      end do
      write(60,'(A,F12.6)') '  -- electron kinetic energy              ', T

      ! electron-nuclear attraction energy
      V = 0.0_dp
      call matmul('C', 'N', oper3, exi_V_j, oper6)
      call matmul('N', 'N', oper6, oper3, oper4)
      do loop_j = 1, electron_count
        V = V + real(oper4(occindex(loop_j),occindex(loop_j)))
      end do
      write(60,'(A,F12.6)') '  -- electron-nuclear attraction energy   ', V

      if (DKH_order == 2) then
        ESOC = 0.0_dp
        call matmul('C', 'N', oper3, exSOC, oper6)
        call matmul('N', 'N', oper6, oper3, oper4)
        do loop_j = 1, electron_count
          ESOC = ESOC + real(oper4(occindex(loop_j),occindex(loop_j)))
        end do
        write(60,'(A,F12.6)') '  -- spin-orbital coupling energy         ', ESOC
        if (SRTP_type) then
          ESR = 0.0_dp
          call matmul('C', 'N', oper3, exSR, oper6)
          call matmul('N', 'N', oper6, oper3, oper4)
          do loop_j = 1, electron_count
            ESR = ESR + real(oper4(occindex(loop_j),occindex(loop_j)))
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
        ! solve residual equation
        ! dgesv and dspsv will cause V_integral_2e conflict for unknown reason
        ! since DIISmat (and its inverse) is real symmetric, plus the column 
        ! vector is simple, use the inverse of DIISmat to solve directly
        call inverse(DIISmat, subsp+1)
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
    if (abs(molE-molE_pre) < conver_tol .and. loop_i < maxiter) then
      if (DKH_order == 0) write(60,'(A)') '  NR SCF succeed!'
      if (DKH_order == 2) write(60,'(A)') '  DKH2 SCF succeed!'
    else
      if (DKH_order == 0) write(60,'(A)') '  NR SCF failed!'
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
    write(60,'(A)') '                           MOL INFO'
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
      write(60,'(A)') '  -- Note: relativistic calculation causes the'
      write(60,'(A)') '  -- Virial ratio to deviate (usually below) 2.0'
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
    if (fx_id /= 0) then
      write(60,'(A)') &
      '  -- Note: there is little theoretical justification'
      write(60,'(A)') &
      '  -- to calculate <S**2> in a DFT calculation.'
    end if
    scf_kappa = Krammers()
    write(60,'(A)') '  deviation parameters from time reversal symmetry (TRS)'
    write(60,'(A,E12.5)') '  -- kappa     =', scf_kappa
    write(60,'(A,E12.5)') '  -- ref kappa =', dsqrt(real(Nalpha-Nbeta,dp))
    write(60,'(A,E12.5)') '  -- SOC kappa =', &
    scf_kappa - dsqrt(real(Nalpha-Nbeta,dp))
    write(60,'(A)') &
    '  ============================================================='
    write(60,*)
    write(60,*)
    write(60,*)
    if (DKH_order /= 0) call RGI()
    write(60,*)
    write(60,*)
    write(60,*)
    write(60,'(A)') &
    '  ============================================================='
    write(60,'(A)') &
    '                       CANONICAL ORB INFO'
    write(60,'(A)') &
    '  ============================================================='
    do loop_i = 1, electron_count
      if (loop_i < electron_count) then
        call calc_S2HForb(loop_i)
        write(60,'(A,I3.3,A,F12.6)') &
        '  HOMO-', electron_count-loop_i, &
        '   energy / Eh                       ... ', orbE(loop_i)
        write(60,'(A3,I3.3,A,F12.6)') &
        '  #',loop_i,'       <S**2> / hbar**2                  ... ', S__2orb
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
        write(60,'(A3,I3.3,A,F12.6)') &
        '  #',loop_i,'       <S**2> / hbar**2                  ... ', S__2orb
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
        write(60,'(A3,I3.3,A,F12.6)') &
        '  #',loop_i,'       <S**2> / hbar**2                  ... ', S__2orb
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
        write(60,'(A3,I3.3,A,F12.6)') &
        '  #',loop_i,'       <S**2> / hbar**2                  ... ', S__2orb
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
      write(60,'(A)') '  dumping AO2MO to '//trim(address_job)//'.molden.d'
      if (DKH_order == 0) then
        write(60,'(A)') &
        '  -- Note: for scalar MOs, only realpart will be generated.'
      end if
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
    deallocate(Rsd, DIISmat)
    deallocate(rho_history, rho_pre, rho_pre_pre)
    deallocate(swint)
    deallocate(Fock1, Fock2HFexc, Fock2HFcol)
    if (allocated(Fock2KScor)) deallocate(Fock2KScor)
    if (allocated(Fock2KSexc)) deallocate(Fock2KSexc)
    deallocate(i_j, i_p2_j, i_V_j, c2s, exc2s, s2f, exs2f, c2f, exc2f)
    if (allocated(m_c2s)) deallocate(m_c2s, m_exc2s)
    if (kill) then
      if (fx_id /= -1) call Fockxc_end()
      deallocate(AO2MO, rho_m, Fock, orbE, oper3, occindex)
    end if
    ndschwarz = .true.
    ini_rho = .true.
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
  subroutine assign_rho()
    implicit none
    integer              :: ploop_i,ploop_j,ploop_k    ! loop variables
    integer              :: Na, Nb
    integer              :: degenlow, degenhigh, load  ! degenerat region
    real(dp)             :: rdMO(5)
    character(len = 512) :: line_str
    if (ini_rho) then
      allocate(rho_m(2*sbdm,2*sbdm))
      allocate(AO2MO(2*sbdm,2*fbdm))
      allocate(occindex(electron_count))
      Nalpha = (electron_count-(spin_mult-1))/2 + (spin_mult-1)
      Nbeta = (electron_count-(spin_mult-1))/2
      ! load MO coefficient
      !-------------------------load from molden-------------------------
      if (guess_type == 'molden') then
        call load_gb_molden(.true.)
        if (m_atom_count /= atom_count) call terminate(&
        'n_atoms in .molden is not consistent with n_atoms in .xyz')
        write(60,'(A)') &
        '  -- basis and geometry in molden were loaded'
        write(60,'(A,I4,A,I4)') '  -- m_cbdm / m_sbdm: ',m_cbdm, ' /', m_sbdm
        allocate(m_AO2MO_a(m_sbdm,m_sbdm))
        allocate(m_AO2MO_b(m_sbdm,m_sbdm))
        allocate(t_AO2MO_a(sbdm,m_sbdm))
        allocate(t_AO2MO_b(sbdm,m_sbdm))
        call load_MO_molden(m_AO2MO_a, m_AO2MO_b)
        write(60,'(A)') '  -- MO coeffs in molden were loaded'
        call m_basis_proj(m_AO2MO_a, t_AO2MO_a)
        call m_basis_proj(m_AO2MO_b, t_AO2MO_b)
        write(60,'(A)') '  -- MO coeffs were projected to job basis'
        rho_m = c0
        do ploop_i = 1,sbdm
          do ploop_j = 1,sbdm
            do ploop_k = 1,Nalpha
              rho_m(ploop_i,ploop_j) = &
              rho_m(ploop_i,ploop_j) + &
              t_AO2MO_a(ploop_i,ploop_k)*t_AO2MO_a(ploop_j,ploop_k)
            end do
            do ploop_k = 1,Nbeta
              rho_m(sbdm+ploop_i,sbdm+ploop_j) = &
              rho_m(sbdm+ploop_i,sbdm+ploop_j) + &
              t_AO2MO_b(ploop_i,ploop_k)*t_AO2MO_b(ploop_j,ploop_k)
            end do
          end do
        end do
        deallocate(m_AO2MO_a, m_AO2MO_b, t_AO2MO_a, t_AO2MO_b)
      !-------------------------load from .ao2mo-------------------------
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
    ! keep spin multiplicity if frontier MOs are degenerate
    ! should keep oper3 and AO2MO unchanged, only change the occindex
    else if (.not.ini_rho .and. (cspin=='n' .or. cspin=='d')) then
      forall (ploop_i=1:electron_count) occindex(ploop_i) = ploop_i
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
      if (cspin=='d' .and. abs((totalpha-totbeta) - &
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
        write(60,'(A,I3,A3,I3)') &
        '  -- -- alpha and beta in nondegenerate region: ', Na, ' / ', Nb
        load = degenlow
        do ploop_k = degenlow, degenhigh
          call calc_S2HForb(ploop_k)
          if (Szorb>0.0 .and. Na<Nalpha) then
            write(60,'(A,I3,A,I3)') '  -- -- load alpha ',ploop_k,' on ',load
            occindex(load) = ploop_k
            load = load + 1
            Na = Na + 1
          else if(Szorb<0.0 .and. Nb<Nbeta) then
            write(60,'(A,I3,A,I3)') '  -- -- load beta ',ploop_k,' on ',load
            occindex(load) = ploop_k
            load = load + 1
            Nb = Nb + 1
          end if
          if (load == electron_count + 1) exit
        end do
        if (load < electron_count + 1) call terminate(&
        "cspin error, can't keep spin in degenerate space")
        rho_m = c0
        do ploop_i = 1, 2*sbdm
          do ploop_j = 1, 2*sbdm
            do ploop_k = 1, electron_count
              rho_m(ploop_i,ploop_j) = rho_m(ploop_i,ploop_j) + &
              AO2MO(ploop_i,occindex(ploop_k)) * &
              conjg(AO2MO(ploop_j,occindex(ploop_k)))
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
      do ploop_k = 1, 2*fbdm
        call calc_S2HForb(ploop_k)
        if (Szorb > 0.0 .and. Na < Nalpha) then
          occindex(load) = ploop_k
          load = load + 1
          Na = Na + 1
        else if (Szorb < 0.0 .and. Nb < Nbeta) then
          occindex(load) = ploop_k
          load = load + 1
          Nb = Nb + 1
        end if
        if (Na==Nalpha .and. Nb==Nbeta) exit
      end do
      rho_m = c0
      do ploop_i = 1, 2*sbdm
        do ploop_j = 1, 2*sbdm
          do ploop_k = 1, electron_count
            rho_m(ploop_i,ploop_j) = rho_m(ploop_i,ploop_j) + &
            AO2MO(ploop_i,occindex(ploop_k)) * &
            conjg(AO2MO(ploop_j,occindex(ploop_k)))
          end do
        end do
      end do
    end if
    if (.not. ini_rho) then
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
    else
      ini_rho = .false.
    end if
  end subroutine assign_rho
  
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
      call matmul('N', 'T', Ap, AO2p2, oper2)
      call matmul('N', 'N', oper2, i_V_j, oper1)
      call matmul('N', 'N', oper1, AO2p2, oper2)
      call matmul('N', 'N', oper2, Ap, oper1)
      AVA(1:fbdm,1:fbdm) = AVA(1:fbdm,1:fbdm) + oper1 * c1
      AVA(fbdm+1:2*fbdm,fbdm+1:2*fbdm) = &
      AVA(fbdm+1:2*fbdm,fbdm+1:2*fbdm) + oper1 * c1
      exi_V_j(1:fbdm,1:fbdm) = i_V_j
      exi_V_j(fbdm+1:2*fbdm,fbdm+1:2*fbdm) = i_V_j
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
      if (SRTP_type) then
        call matmul('N', 'N', SRp, ARVRA, oper3)
        call matmul('N', 'N', oper3, SRp, oper4)
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
        call matmul('N', 'T', ApRp, AO2p2, oper2)
        call matmul('N', 'N', oper2, temp_pool, oper1)
        call matmul('N', 'N', oper1, AO2p2, oper2)
        call matmul('N', 'N', oper2, ApRp, oper1)
        Fock1(1:fbdm,1:fbdm) = Fock1(1:fbdm,1:fbdm) - oper1/(2.0_dp*c2) * c1
        Fock1(fbdm+1:2*fbdm,fbdm+1:2*fbdm) = &
        Fock1(fbdm+1:2*fbdm,fbdm+1:2*fbdm) - oper1/(2.0_dp*c2) * c1
        exSR(1:fbdm,1:fbdm) = exSR(1:fbdm,1:fbdm) - oper1/(2.0_dp*c2) * c1
        exSR(fbdm+1:2*fbdm,fbdm+1:2*fbdm) = &
        exSR(fbdm+1:2*fbdm,fbdm+1:2*fbdm) - oper1/(2.0_dp*c2) * c1
        !----------------------
        ! ApRp (px3Vpy-py3Vpx+pxVpy3-pyVpx3) ApRp
        temp_pool = px3Vpy-py3Vpx+transpose(py3Vpx-px3Vpy)
        call matmul('N', 'T', ApRp, AO2p2, oper2)
        call matmul('N', 'N', oper2, temp_pool, oper1)
        call matmul('N', 'N', oper1, AO2p2, oper2)
        call matmul('N', 'N', oper2, ApRp, oper1)
        Fock1(1:fbdm,1:fbdm) = Fock1(1:fbdm,1:fbdm) - oper1/(2.0_dp*c2) * ci
        Fock1(fbdm+1:2*fbdm,fbdm+1:2*fbdm) = &
        Fock1(fbdm+1:2*fbdm,fbdm+1:2*fbdm) + oper1/(2.0_dp*c2) * ci
        exSR(1:fbdm,1:fbdm) = exSR(1:fbdm,1:fbdm) - oper1/(2.0_dp*c2) * ci
        exSR(fbdm+1:2*fbdm,fbdm+1:2*fbdm) = &
        exSR(fbdm+1:2*fbdm,fbdm+1:2*fbdm) + oper1/(2.0_dp*c2) * ci
        !----------------------
        ! ApRp (pz3Vpx-px3Vpz+pzVpx3-pxVpz3) ApRp
        temp_pool = pz3Vpx-px3Vpz+transpose(px3Vpz-pz3Vpx)
        call matmul('N', 'T', ApRp, AO2p2, oper2)
        call matmul('N', 'N', oper2, temp_pool, oper1)
        call matmul('N', 'N', oper1, AO2p2, oper2)
        call matmul('N', 'N', oper2, ApRp, oper1)
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
        call matmul('N', 'T', ApRp, AO2p2, oper2)
        call matmul('N', 'N', oper2, temp_pool, oper1)
        call matmul('N', 'N', oper1, AO2p2, oper2)
        call matmul('N', 'N', oper2, ApRp, oper1)
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
          Ve(floop_i,floop_j) = i_V_j(floop_i,floop_j) / &
          (c*(edc(floop_i) + edc(floop_j)))
          pxVepx(floop_i,floop_j) = pxVpx(floop_i,floop_j) / &
          (c*(edc(floop_i) + edc(floop_j)))
          pyVepy(floop_i,floop_j) = pyVpy(floop_i,floop_j) / &
          (c*(edc(floop_i) + edc(floop_j)))
          pzVepz(floop_i,floop_j) = pzVpz(floop_i,floop_j) / &
          (c*(edc(floop_i) + edc(floop_j)))
          pxVepy(floop_i,floop_j) = pxVpy(floop_i,floop_j) / &
          (c*(edc(floop_i) + edc(floop_j)))
          pyVepx(floop_i,floop_j) = pyVpx(floop_i,floop_j) / &
          (c*(edc(floop_i) + edc(floop_j)))
          pxVepz(floop_i,floop_j) = pxVpz(floop_i,floop_j) / &
          (c*(edc(floop_i) + edc(floop_j)))
          pzVepx(floop_i,floop_j) = pzVpx(floop_i,floop_j) / &
          (c*(edc(floop_i) + edc(floop_j)))
          pyVepz(floop_i,floop_j) = pyVpz(floop_i,floop_j) / &
          (c*(edc(floop_i) + edc(floop_j)))
          pzVepy(floop_i,floop_j) = pzVpy(floop_i,floop_j) / &
          (c*(edc(floop_i) + edc(floop_j)))
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
      forall (floop_i = 1:2*fbdm)
        SRp(floop_i, floop_i) = 0.5_dp * c1
      end forall
      call matmul('N', 'N', ARVeRA, oper5, oper3)
      call matmul('N', 'N', oper3, AVA, oper4)
      Fock1 = Fock1 - oper4
      call matmul('N', 'N', ARVRA, oper5, oper3)
      call matmul('N', 'N', oper3, AVeA, oper4)
      Fock1 = Fock1 - oper4
      call matmul('N', 'N', AVeA, oper5, oper3)
      call matmul('N', 'N', oper3, ARVRA, oper4)
      Fock1 = Fock1 - oper4
      call matmul('N', 'N', AVA, oper5, oper3)
      call matmul('N', 'N', oper3, ARVeRA, oper4)
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
      call matmul('N', 'N', ARVeRA, oper5, oper3)
      call matmul('N', 'N', oper3, ARVRA, oper4)
      Fock1 = Fock1 + oper4
      call matmul('N', 'N', ARVRA, oper5, oper3)
      call matmul('N', 'N', oper3, ARVeRA, oper4)
      Fock1 = Fock1 + oper4
      forall (floop_i = 1:fbdm) ! RI extract
        oper5(floop_i,floop_i) = &
        0.5_dp * (evl_p2(floop_i)/(edc(floop_i)+c)**2) * c1
        oper5(fbdm+floop_i,fbdm+floop_i) = &
        0.5_dp * (evl_p2(floop_i)/(edc(floop_i)+c)**2) * c1
      end forall
      call matmul('N', 'N', AVeA, oper5, oper3)
      call matmul('N', 'N', oper3, AVA, oper4)
      Fock1 = Fock1 + oper4
      call matmul('N', 'N', AVA, oper5, oper3)
      call matmul('N', 'N', oper3, AVeA, oper4)
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
      call matmul('N', 'N', exAO2p2, Fock1, oper3)
      call matmul('N', 'T', oper3, exAO2p2, Fock1)
      call matmul('N', 'N', exAO2p2, exi_T_j, oper3)
      call matmul('N', 'T', oper3, exAO2p2, exi_T_j)
      call matmul('N', 'N', exAO2p2, exi_V_j, oper3)
      call matmul('N', 'T', oper3, exAO2p2, exi_V_j)
      call matmul('N', 'N', exAO2p2, exSOC, oper3)
      call matmul('N', 'T', oper3, exAO2p2, exSOC)
      if (SRTP_type) then
        call matmul('N', 'N', exAO2p2, exSR, oper3)
        call matmul('N', 'T', oper3, exAO2p2, exSR)
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
    complex(dp) :: supp(2*sbdm,2*cbdm)
    !DIR$ ATTRIBUTES ALIGN:align_size :: Fock2_assigned
    integer :: Fock2_assigned(2,8)         ! avoid duplicate assignment of Fock2
    integer :: assigned
    integer :: iassign
    logical :: carry
    if (ndschwarz) then
      ndschwarz = .false.
      write(60,'(a)') '  -- Schwarz screening of <ij||kl>'
      allocate(swint(cbdm,cbdm))
      swint = 0.0_dp
      do ui = 1, cbdm
        contri = basis_inf(ui) % contr
        Li     = basis_inf(ui) % L
        Mi     = basis_inf(ui) % M
        expi(1:contri) = basis_inf(ui) % expo(1:contri)
        coei(1:contri) = basis_inf(ui) % Ncoe(1:contri,Mi)
        codi = basis_inf(ui) % pos
        do uj = ui, cbdm
          contrj = basis_inf(uj) % contr
          Lj     = basis_inf(uj) % L
          Mj     = basis_inf(uj) % M
          expj(1:contrj) = basis_inf(uj) % expo(1:contrj)
          coej(1:contrj) = basis_inf(uj) % Ncoe(1:contrj,Mj)
          codj = basis_inf(uj) % pos
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
    ! transform rho_m to Cartesian basis
    call matmul('N', 'C', rho_m, exc2s, supp)
    deallocate(rho_m)
    allocate(rho_m(2*cbdm,2*cbdm))
    call matmul('N', 'N', exc2s, supp, rho_m)
    if (allocated(Fock2HFcol)) deallocate(Fock2HFcol)
    allocate(Fock2HFcol(2*cbdm,2*cbdm), source=c0)
    if (allocated(Fock2HFexc)) deallocate(Fock2HFexc)
    allocate(Fock2HFexc(2*cbdm,2*cbdm), source=c0)
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
        ui     = i
        contri = basis_inf(ui) % contr
        Li     = basis_inf(ui) % L
        Mi     = basis_inf(ui) % M
        expi(1:contri) = basis_inf(ui) % expo(1:contri)
        coei(1:contri) = basis_inf(ui) % Ncoe(1:contri,Mi)
        codi = basis_inf(ui) % pos
        !----------------------------<uk>-----------------------------------
        uk     = j
        contrk = basis_inf(uk) % contr
        Lk     = basis_inf(uk) % L
        Mk     = basis_inf(uk) % M
        expk(1:contrk) = basis_inf(uk) % expo(1:contrk)
        coek(1:contrk) = basis_inf(uk) % Ncoe(1:contrk,Mk)
        codk = basis_inf(uk) % pos
        !----------------------------<uj>---------------------------------
        do uj = ui, 1, -1
          contrj = basis_inf(uj) % contr
          Lj     = basis_inf(uj) % L
          Mj     = basis_inf(uj) % M
          expj(1:contrj) = basis_inf(uj) % expo(1:contrj)
          coej(1:contrj) = basis_inf(uj) % Ncoe(1:contrj,Mj)
          codj = mol(basis_inf(uj)%atom) % pos
          !----------------------------<ul>-------------------------------
          do ul = min(uk,uj+(ui*(ui-1)-uk*(uk-1))/2), 1, -1
            ! Schwarz screening, |<ij||kl>| <= dsqrt(<ij||ij>) * dsqrt(<kl||kl>)
            if (dsqrt(swint(ui,uj)*swint(uk,ul)) < schwarz_VT) cycle
            contrl = basis_inf(ul) % contr
            Ll     = basis_inf(ul) % L
            Ml     = basis_inf(ul) % M
            expl(1:contrl) = basis_inf(ul) % expo(1:contrl)
            coel(1:contrl) = basis_inf(ul) % Ncoe(1:contrl,Ml)
            codl = basis_inf(ul) % pos
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
    Fock2HFcol = Fock2HFcol + HFcol_mic
    Fock2HFexc = Fock2HFexc + HFexc_mic
    !$omp end critical
    !$omp end parallel
    ! assign Kohn-Sham matrices
    if (fx_id /= -1) then
      if (allocated(Fock2KSexc)) deallocate(Fock2KSexc)
      allocate(Fock2KSexc(2*cbdm, 2*cbdm))
      if (fc_id /= -1) then
        if (allocated(Fock2KScor)) deallocate(Fock2KScor)
        allocate(Fock2KScor(2*cbdm, 2*cbdm))
        call basis2grid_Becke(rho_m, KSexc, KScor, Fock2KSexc, Fock2KScor)
        call cfgo(Fock2KSexc)
        call cfgo(Fock2KScor)
      else
        call basis2grid_Becke(rho_m=rho_m, ex=KSexc, Fockx=Fock2KSexc)
        call cfgo(Fock2KSexc)
      end if
    end if
    call csgo(rho_m) ! transform rho_m to spherical-harmonic basis
    call cfgo(Fock2HFcol)
    call cfgo(Fock2HFexc)
  end subroutine Fock2e
  
!------------------------------------------------------------
!> calculate <S**2> based on oper3 generated by SCF
!!
!! divide all occupied orbitals into two sets of alpha, beta orbitals with
!!
!! their original coefficients, then use UHF method to solve Lowdin <S**2>.
  subroutine calc_S2HF()
    implicit none
    integer     :: cloop_i, cloop_j, cloop_k ! loop variables for calc_S2HF
    complex(dp) :: alal, albe, beal, bebe    ! components of pair density matrix
    complex(dp) :: suppa, suppb
    ! oper3 and i_j are orthogonally normalized
    totalpha = 0.0_dp
    totbeta = 0.0_dp
    do cloop_i = 1, electron_count
      do cloop_k = 1, fbdm
        totalpha = totalpha + &
        real(conjg(oper3(cloop_k,occindex(cloop_i))) * &
        oper3(cloop_k,occindex(cloop_i)),dp)
      end do
      do cloop_k = fbdm+1, 2*fbdm
        totbeta = totbeta + &
        real(conjg(oper3(cloop_k,occindex(cloop_i))) * &
        oper3(cloop_k,occindex(cloop_i)),dp)
      end do
    end do
    alal = totalpha * (totalpha-1.0_dp)
    bebe = totbeta * (totbeta-1.0_dp)
    albe = c0
    beal = c0
    do cloop_i = 1, electron_count
      do cloop_j = 1, electron_count
        suppa = c0
        suppb = c0
        do cloop_k = 1, fbdm
          suppa = suppa + &
          conjg(oper3(cloop_k,occindex(cloop_i))) * &
          oper3(fbdm+cloop_k,occindex(cloop_j))
          suppb = suppb + &
          conjg(oper3(fbdm+cloop_k,occindex(cloop_j))) * &
          oper3(cloop_k,occindex(cloop_i))
        end do
        albe = albe + suppa*suppb
        suppa = c0
        suppb = c0
        do cloop_k = fbdm+1, 2*fbdm
          suppa = suppa + &
          conjg(oper3(cloop_k,occindex(cloop_i))) * &
          oper3(cloop_k-fbdm,occindex(cloop_j))
          suppb = suppb + &
          conjg(oper3(cloop_k-fbdm,occindex(cloop_j))) * &
          oper3(cloop_k,occindex(cloop_i))
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
    integer            :: zloop_i     !loop variables for calc_S2HForb
    real(dp)           :: suppa, suppb
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
  integer                 :: ii, jj, kk, mm ! loop variables
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
  if (DKH_order == 0) then
    write(60,'(A)') '  RGI: DKH_order = 0, unable to perform'
    return
  end if

  ! designation of basical configuration, BCSmult and BCSz
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
  integer           :: ii, jj          ! loop variables
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
!> project AO2MO coefficient from m_basis to basis
!!
!! cB = i_j_s^(-1) . m_i_j . cA . m_s2fT
!!
!! cA(m_sbdm, m_sbdm), cB(sbdm, m_sbdm) number of project MOs is m_sbdm
  subroutine m_basis_proj(cA, cB)
    implicit none
    real(dp),intent(in)  :: cA(m_sbdm, m_sbdm)
    real(dp),intent(out) :: cB(sbdm, m_sbdm)
    integer              :: pjloop_i, pjloop_j, pjloop_k, pjloop_m
    integer              :: contri, contrj
    integer              :: Li, Lj
    integer              :: Mi, Mj
    real(dp)             :: expi(16), expj(16)
    real(dp)             :: coei(16), coej(16)
    real(dp)             :: codi(3), codj(3)
    real(dp),allocatable :: spp(:,:), spp2(:,:)
    real(dp)             :: m_s2f(m_sbdm,m_sbdm)
    real(dp)             :: i_j_inv(sbdm,sbdm)
    real(dp)             :: m_min_evl
    real(dp)             :: itmat(sbdm, m_sbdm)
    real(dp)             :: sum
    if (.not. allocated(m_i_j)) then
      ! assign m_i_j
      allocate(m_i_j(cbdm,m_cbdm), source=0.0_dp)
      do pjloop_i = 1, cbdm
        contri = basis_inf(pjloop_i) % contr
        Li     = basis_inf(pjloop_i) % L
        Mi     = basis_inf(pjloop_i) % M
        expi(1:contri) = basis_inf(pjloop_i) % expo(1:contri)
        coei(1:contri) = basis_inf(pjloop_i) % Ncoe(1:contri,Mi)
        codi = basis_inf(pjloop_i) % pos
        do pjloop_j = 1, m_cbdm
          contrj = m_basis_inf(pjloop_j) % contr
          Lj     = m_basis_inf(pjloop_j) % L
          Mj     = m_basis_inf(pjloop_j) % M
          expj(1:contrj) = m_basis_inf(pjloop_j) % expo(1:contrj)
          coej(1:contrj) = m_basis_inf(pjloop_j) % Ncoe(1:contrj,Mj)
          codj = m_basis_inf(pjloop_j) % pos
          do pjloop_k = 1, contri
            do pjloop_m = 1, contrj
              m_i_j(pjloop_i,pjloop_j) = m_i_j(pjloop_i,pjloop_j) + &
              Gaussian_Product_Integral(coei(pjloop_k)*coej(pjloop_m),&
              AO_fac(:,Li,Mi),AO_fac(:,Lj,Mj),& 
              expi(pjloop_k),expj(pjloop_m),codi,codj)
            end do
          end do
        end do
      end do
      call m_assign_cs()
      call m_csgo(m_i_j)
    end if
    i_j_inv = i_j_s
    call inverse(i_j_inv, sbdm)
    ! itmat = i_j_s^(-1) . m_i_j . cA, itmat is not orthogonal
    allocate(spp(sbdm,m_sbdm))
    call matmul('N', 'N', m_i_j, cA, spp)
    call matmul('N', 'N', i_j_inv, spp, itmat)
    deallocate(spp)
    ! symmetric orthogonalisation
    ! m_s2fT.cBT.i_j_s.cB.m_s2f = I
    ! note that this m_s2f is not same as s2f which satisfy s2f.i_j_s.s2fT = I
    ! cB is not orthogonal so we have to make whole MO(cBT.i_j_s.cB)
    ! orthogonal instead of just AO(i_j_s)
    allocate(spp(m_sbdm, sbdm))
    allocate(spp2(m_sbdm, m_sbdm))
    call matmul('T', 'N', itmat, i_j_s, spp)
    call matmul('N', 'N', spp, itmat, spp2)
    call symm_orth(spp2, m_sbdm, m_s2f, m_min_evl)
    ! unlikely to meet linear dependency, because MO satisfy Pauli exclusion
    ! principle so will not repeat each other
    call matmul('N', 'N', itmat, m_s2f, cB)
  end subroutine m_basis_proj

!-----------------------------------------------------------------------
!> dump mol orbital information to .molden.d file
! AO2MO is like:
!                  |                  |
!                  |                  |
!                  |                  |
!                  |                  |
! a part1 real&img | a part2 real&img |
!                  |                  |
!                  |                  |
!__________________|__________________|
!                  |                  |
!                  |                  |
!                  |                  |
!                  |                  |
! b part1 real&img | b part1 real&img |
!                  |                  |
!                  |                  |
!__________________|__________________|
  subroutine dump_molden()
    implicit none
    character(len=200)     :: dir
    integer                :: channel
    integer                :: dmi, dmj, dmk     ! loop variables
    dir = trim(address_job)//'.molden.d'
    call execute_command_line('mkdir -p '//trim(dir), wait=.true., exitstat=ios)
    if (ios /= 0) call terminate(&
    "dump to molden failed, molden.d can't be created")
    if (.not. allocated(AO2MO)) &
    call terminate('dump to molden failed, AO2MO is empty')
    
    ! molden file contains the real part1 of AO2MO
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
      call calc_S2HForb(dmi)
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
      call calc_S2HForb(dmi)
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

    ! molden file contains the real part2 of AO2MO
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
      call calc_S2HForb(dmi)
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
      call calc_S2HForb(dmi)
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
    
    if (DKH_order == 2) then
      ! molden file contains the imaginary part1 of AO2MO
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
        call calc_S2HForb(dmi)
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
        call calc_S2HForb(dmi)
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

      ! molden file contains the imaginary part2 of AO2MO
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
        call calc_S2HForb(dmi)
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
        call calc_S2HForb(dmi)
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
  end subroutine dump_molden
  
end module SCF