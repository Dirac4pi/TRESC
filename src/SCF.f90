!> @file SCF.f90
!!
!! @brief single-configuration self-consistent field calculation
!!
!! @author dirac4pi
module SCF
  use Hamiltonian
  use Atoms
  use Fundamentals
  
  ! density matrices and molecular orbital coefficients
  integer Nalpha, Nbeta                    ! number of alpha and beta elctron
  ! spinor MO coefficients, in order of (AO1,0), (0,AO1), ... (AOn,0), (0,AOn)
  complex(dp),allocatable :: AO2MO(:,:)
  real(dp),allocatable :: Gaualpha(:)      ! Gaussian alpha orbital coefficient
  real(dp),allocatable :: AO2MOalpha(:,:)  ! Gaussian alpha orbital coefficient
  real(dp),allocatable :: Gaubeta(:)       ! Gaussian beta orbital coefficient
  real(dp),allocatable :: AO2MObeta(:,:)   ! Gaussian beta orbital coefficient
  complex(dp),allocatable :: rou_m(:,:)    ! density matrix, complex Hermitian
  complex(dp),allocatable :: rotation(:,:) ! rotate one orbital with another
  real(dp) :: RMSDP, maxDP
  real(dp),allocatable :: AOsupp(:,:)
  integer :: evl_count_f                   ! number of eigenvalues in Fock diag
  complex(dp),allocatable :: Fock(:,:)     ! Fock matrix
  integer,allocatable :: isupp_ev_f(:)
  integer :: fliwork
  integer :: flwork
  integer :: lrwork
  complex(dp),allocatable :: fwork(:)      ! work of Fock for zheevr input
  real(dp),allocatable :: rwork(:)         ! rwork of Fock for zheevr input
  integer,allocatable :: fiwork(:)         ! iwork of Fock for zheevr input
  integer :: finfo                         ! info of calling lapack functions
  
  logical :: ini_rou = .true.              ! initial density matrix loaded
  logical :: diiscuted = .false.           ! whether cutdiis ever reached
  
  !--------------<one electron Fock>--------------
  complex(dp),allocatable :: Fock1(:,:)    ! one electron Fock matrix
  real(dp),allocatable :: oper1(:,:)       ! operator matrix
  real(dp),allocatable :: oper2(:,:)       ! operator matrix
  complex(dp),allocatable :: oper3(:,:)    ! operator matrix
  complex(dp),allocatable :: oper4(:,:)    ! operator matrix
  complex(dp),allocatable :: oper5(:,:)    ! operator matrix
  complex(dp),allocatable :: oper6(:,:)    ! operator matrix
  real(dp),allocatable :: Ve(:,:)          ! V/(E+E')
  real(dp),allocatable :: pxVepx(:,:)      ! pxVpx/(E+E')
  real(dp),allocatable :: pyVepy(:,:)      ! pyVpy/(E+E')
  real(dp),allocatable :: pzVepz(:,:)      ! pzVpz/(E+E')
  real(dp),allocatable :: pxVepy(:,:)      ! pxVpy/(E+E')
  real(dp),allocatable :: pyVepx(:,:)      ! pyVpx/(E+E')
  real(dp),allocatable :: pxVepz(:,:)      ! pxVpz/(E+E')
  real(dp),allocatable :: pzVepx(:,:)      ! pzVpx/(E+E')
  real(dp),allocatable :: pyVepz(:,:)      ! pyVpz/(E+E')
  real(dp),allocatable :: pzVepy(:,:)      ! pzVpy/(E+E')
  real(dp),allocatable :: Ap(:,:)          ! Ap
  real(dp),allocatable :: ApRp(:,:)        ! ApRp
  complex(dp),allocatable :: SRp(:,:)      ! ApRp(1+p^2/4c^2)
  complex(dp),allocatable :: ARVRA(:,:)    ! ApRpVRpAp
  complex(dp),allocatable :: ARVeRA(:,:)   ! ApRp(V/(E+E'))RpAp
  complex(dp),allocatable :: AVA(:,:)      ! ApVAp
  complex(dp),allocatable :: AVeA(:,:)     ! Ap(V/(E+E'))Ap
  complex(dp),allocatable :: exAO2p2(:,:)  ! extended AO2p2 matrix
  
  !--------------<two electron Fock>--------------
  complex(dp),allocatable :: Fock2(:,:)      ! two electron Fock matrix
  complex(dp),allocatable :: Fock2_mic(:,:)  ! two electron Fock matrix
  
  real(dp),allocatable :: swintegral(:,:)    ! <ij||ij> as well as <kl||kl>
  logical :: ndschwarz = .true.
  integer :: Fock2_assigned(2,8)          ! avoid duplicate assignment of Fock2
  integer :: assigned
  integer :: iassign
  logical :: carry
  
  complex(dp),allocatable :: supp1(:,:), supp2(:,:)
  
  ! orbital energy and molecular energy
  real(dp),allocatable :: orbE(:)  ! orbital energy
  real(dp) :: molE_pre, molE       ! molecular energy
  real(dp) :: nucE                 ! nuclear repulsion energy
  real(dp) :: eleE                 ! electron repulsion energy
  real(dp) :: T                    ! kinetic energy
  real(dp) :: V                    ! electron-nuclear attraction energy
  real(dp) :: ESOC                 ! SOC energy
  real(dp) :: ESR                  ! SRTP energy
  real(dp) :: emd4                 ! dispersion energy calculated by DFT-D4

  ! damping & DIIS(AX=B)
  ! privious rou_m of current iteration
  complex(dp),allocatable :: rou_pre(:,:)
  ! privious rou_m of privious iteration
  complex(dp),allocatable :: rou_pre_pre(:,:)
  complex(dp),allocatable :: rou_history(:,:,:)  ! coeff of subsp iteration
  complex(dp),allocatable :: Rsd(:,:,:)          ! residuals of rou_m
  complex(dp),allocatable :: DIISmat(:,:)        ! A
  complex(dp),allocatable :: DIISwork(:)
  integer :: lDIISwork
  integer,allocatable :: ipiv(:)
  integer :: DIIsinfo
  real(dp) :: damp_coe
  real(dp) :: dE_pre
  logical :: forward = .true.
  
  contains

!------------------------------------------------------------
!> Hartree-Fock SCF procedure for DKH0 / DKH2 Hamiltonian
  subroutine DKH_Hartree_Fock(keep, kill)
    implicit none
    integer,intent(in) :: keep   ! initialization level
    logical,intent(in) :: kill   ! kill the process after SCF
    integer :: iatom, ishell
    write(60,'(A)') 'Module SCF:'
    write(60,'(A)') '  construct one electron Fock matrix'
    call Fock1e()
    write(60,'(A)') '  complete! stored in Fock1'
    write(60,'(A)') '  calculate nuclear repulsion energy'
    nucE = 0.0_dp
    do loop_i = 1, atom_count
      do loop_j = loop_i + 1, atom_count
        nucE = nucE + (real(molecular(loop_i) % atom_number) * &
        real(molecular(loop_j) % atom_number))/sqrt((molecular(loop_i) % &
        nucleus_position(1) - molecular(loop_j) % nucleus_position(1))**2 + &
        (molecular(loop_i) % nucleus_position(2) - molecular(loop_j) % &
        nucleus_position(2))**2 + (molecular(loop_i) % nucleus_position(3) - &
        molecular(loop_j) % nucleus_position(3))**2)
      end do
    end do
    if (isnan(nucE)) call terminate(&
    'nuclear repulsive energy anomaly, may due to overlap atomic coordinates')
    write(60,'(A,F12.7,A)') &
    '  complete! nuclear repulsive energy = ', nucE, ' Eh'
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
    ! AO2MO(new) = ��(i,subsp) DIIScoe(i)*(rou_history(i)+diisdamp*Rsd(i))
    allocate(Rsd(subsp,2*sbdm,2*sbdm))
    allocate(DIISmat(subsp+1,subsp+1))
    allocate(ipiv(subsp+1))
    allocate(rou_history(subsp, 2*sbdm, 2*sbdm))
    allocate(rou_pre(2*sbdm, 2*sbdm))
    allocate(rou_pre_pre(2*sbdm, 2*sbdm))
    do loop_i = 1, maxiter
      if (loop_i /= 1) &
      write(60,*)
      write(60,*)
      write(60,'(A,I3)') '  SCF iter ',loop_i
      if (loop_i == 1) then
        write(60,'(A)') '  read density matrix'
        call assign_rou()
        write(60,'(A)') '  complete! stored in rou_m'
      end if
      write(60,'(A)') '  construct two electron Fock matrix'
      call Fock2e()
      write(60,'(A)') '  complete! stored in Fock2'
      write(60,'(A)') '  diagonalization of Fock matrix'
      
      Fock = Fock1 + Fock2
      
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
      eleE = 0.0_dp
      call zgemm('C', 'N', 2*fbdm, 2*fbdm, 2*fbdm, c1, oper3, 2*fbdm, &
      Fock2, 2*fbdm, c0, oper6, 2*fbdm)
      call zgemm( 'N', 'N', 2*fbdm, 2*fbdm, 2*fbdm, c1, oper6, 2*fbdm, &
      oper3, 2*fbdm, c0, oper4, 2*fbdm)
      do loop_j = 1, electron_count
        eleE = eleE + real(oper4(loop_j,loop_j))
      end do
      write(60,'(A,F12.6)') '  -- electron repulsive energy            ', eleE
      ! calculate molE
      if (loop_i /= 1) molE_pre = molE
      molE = nucE - eleE/2.0_dp
      do loop_j = 1, electron_count
        molE = molE + orbE(loop_j)
      end do
      T = 0.0_dp
      call zgemm( 'C', 'N', 2*fbdm, 2*fbdm, 2*fbdm, c1, oper3, 2*fbdm, &
      exi_T_j, 2*fbdm, c0, oper6, 2*fbdm)
      call zgemm( 'N', 'N', 2*fbdm, 2*fbdm, 2*fbdm, c1, oper6, 2*fbdm, &
      oper3, 2*fbdm, c0, oper4, 2*fbdm)
      do loop_j = 1, electron_count
        T = T + real(oper4(loop_j,loop_j))
      end do
      write(60,'(A,F12.6)') '  -- electron kinetic energy              ', T
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
      if (DKH_order == 0) then
        write(60,'(A,F12.6)') &
        '  -- -<V>/<T>                             ', -(molE-T)/T
      else if (DKH_order == 2) then
        if (SRTP_type) then
          write(60,'(A,F12.6)') &
          '  -- -<V>/<T>                             ', -(molE-T-ESOC-ESR)/T
        else
          write(60,'(A,F12.6)') &
          '  -- -<V>/<T>                             ', -(molE-T-ESOC)/T
        end if
      end if
      ! de-orthogonalization
      call zgemm( 'N', 'N', 2*sbdm, 2*fbdm, 2*fbdm, c1, exXm, 2*sbdm, &
      oper3, 2*fbdm, c0, AO2MO, 2*sbdm)
      write(60,'(A)') '  AO2MO dump to .ao2mo file'
      call dump_matrix_cmplx(name='ao2mo', m=AO2MO, dmi=2*sbdm, dmj=2*fbdm)
      ! convergence check
      if (loop_i == 1) then
        write(60,'(A,F12.6)') '  SCF energy (A.U.) = ',molE
      else
        write(60,'(A,F12.6,A,E10.3)') '  SCF energy (A.U.) = ',molE,&
        '; delta E = ',molE - molE_pre
        if (abs(molE - molE_pre) < conver_tol .and. &
        forward .and. damp_coe < 0.01) then
          write(60,'(A)') '  convergence tolerance met, SCF done!'
          exit
        else
          write(60,'(A)') '  convergence tolerance not met'
        end if
      end if
      write(60,'(A)') '  construct density matrix'
      call assign_rou()
      write(60,'(A)') '  complete! stored in rou_m'
      write(60,'(A)') '  DIIS information'
      if ((abs(molE - molE_pre) < cutdiis .or. diiscuted) .and. forward) then
        diiscuted = .true.
        write(60,*) '  -- DIIS cutted'
        cycle
      end if
      ! generate next AO2MO by DIIS method
      if (loop_i <= nodiis - subsp) then
        write(60,'(A)') '  -- no DIIS acceleration'
        
        
        !--------------<damping>-----------------
        if (loop_i > 2 .and. abs(molE - molE_pre) >= 1.5*abs(dE_pre)) then
          damp_coe = damp_coe + (1.0_dp-damp_coe)*0.5_dp
          if (damp_coe < 0.91) then
            rou_m = (1.0_dp - damp_coe) * rou_pre + damp_coe * rou_pre_pre
            write(60,'(A)') '  -- fallback '
            forward = .false.
            cycle
          end if
        end if
        if (abs(molE - molE_pre) >= cutdamp) then
          damp_coe = damp
          rou_m = (1.0 - damp_coe) * rou_m + damp_coe * rou_pre
          write(60,'(A,F6.3)') '  -- damped ', damp_coe
        else if (abs(molE - molE_pre) >= cutdamp/100.0) then
          damp_coe = 0.5_dp*damp*log10(abs(molE - molE_pre)) + &
          damp*(1.0_dp-0.5_dp*log10(cutdamp))
          rou_m = (1.0 - damp_coe) * rou_m + damp_coe * rou_pre
          write(60,'(A,F6.3)') '  -- damped ', damp_coe
        else
          damp_coe = 0.0_dp
          write(60,'(A)') '  -- undamped '
        end if
        dE_pre = molE - molE_pre
        !--------------<damping>-----------------
        
        
      else if (nodiis - subsp < loop_i .and. loop_i <= nodiis) then
        ! update Rsd
        do loop_j = 1, 2*sbdm
          do loop_k = 1, 2*sbdm
            Rsd(loop_i-(nodiis-subsp), loop_j, loop_k) = &
            rou_m(loop_j, loop_k) - rou_pre(loop_j, loop_k)
          end do
        end do
        ! update rou_history
        do loop_j = 1, 2*sbdm
          do loop_k = 1, 2*sbdm
            rou_history(loop_i-(nodiis-subsp), loop_j, loop_k) = &
            rou_m(loop_j, loop_k)
          end do
        end do
        write(60,'(A,I2,A,I2)') &
        '  -- DIIS subspace filling ',loop_i-(nodiis-subsp),'/',subsp
      else
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
            Rsd(subsp, loop_j, loop_k) = rou_m(loop_j, loop_k) - &
            rou_pre(loop_j, loop_k)
          end do
        end do
        ! update rou_history
        do loop_j = 2, subsp
          do loop_k = 1, 2*sbdm
            do loop_l = 1, 2*sbdm
              rou_history(loop_j - 1, loop_k, loop_l) = &
              rou_history(loop_j, loop_k, loop_l)
            end do
          end do
        end do
        do loop_j = 1, 2*sbdm
          do loop_k = 1, 2*sbdm
            rou_history(subsp, loop_j, loop_k) = rou_m(loop_j, loop_k)
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
        ! dgesv and dspsv will cause V_Integral_2e conflict for unknown reason
        ! since DIISmat (and its inverse) is real symmetric, plus the column 
        ! vector is simple, use the inverse of DIISmat to solve directly
        call zgetrf(subsp+1, subsp+1, DIISmat, subsp+1, ipiv, DIISinfo)
        if (DIISinfo < 0) then
          call terminate('DIIS solution failure, illegal input of dgetrf')
        else if (DIISinfo > 0) then
          call terminate('DIIS solution failure, internal error of dgetrf')
        end if
        allocate(DIISwork(1))
        call zgetri(subsp+1, DIISmat, subsp+1, ipiv, DIISwork, -1, DIISinfo)
        lDIISwork = nint(real(DIISwork(1)))
        deallocate(DIISwork)
        allocate(DIISwork(lDIISwork))
        call zgetri(subsp+1,DIISmat,subsp+1,ipiv,DIISwork,lDIISwork,DIISinfo)
        deallocate(DIISwork)
        if (DIISinfo < 0) then
          call terminate('DIIS solution failure, illegal input of dgetri')
        else if (DIISinfo > 0) then
          call terminate('DIIS solution failure, internal error of dgetri')
        end if
        ! generate new rou_m
        rou_m = c0
        do loop_j = 1, subsp
          do loop_k = 1, 2*sbdm
            do loop_l = 1, 2*sbdm
              rou_m(loop_k,loop_l) = rou_m(loop_k,loop_l) + &
              DIISmat(loop_j,subsp+1) * (rou_history(loop_j,loop_k,loop_l) + &
              diisdamp*Rsd(loop_j,loop_k,loop_l))
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
        if (loop_i > 1) rou_pre_pre = rou_pre
        rou_pre = rou_m
      end if
      forward = .true.
    end do
    if (abs(molE - molE_pre) < conver_tol .and. loop_i < maxiter) then
      if (DKH_order == 0) write(60,'(A)') '  DKH0 Hartree-Fock SCF succeed!'
      if (DKH_order == 2) write(60,'(A)') '  DKH2 Hartree-Fock SCF succeed!'
    else
      if (DKH_order == 0) write(60,'(A)') '  DKH0 Hartree-Fock SCF failed!'
      if (DKH_order == 2) write(60,'(A)') '  DKH2 Hartree-Fock SCF failed!'
    end if
    if (d4) then
      emd4 = dftd4('hf')
      molE = molE + emd4
    end if
    ! print final wave function information
    write(60,*)
    write(60,*)
    write(60,'(A)') &
    '  ============================================================='
    write(60,'(A)') '                        MOLECULAR INFO'
    write(60,'(A)') &
    '  ============================================================='
    write(60,'(A,F12.6)') &
    '  total electronic energy / Eh                  ...',molE
    write(60,'(A,F12.6)') &
    '  nuclear repulsive energy / Eh                 ...',nucE
    write(60,'(A,F12.6)') &
    '  electron repulsive energy / Eh                ...',eleE
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
    if (DKH_order == 0) then
      write(60,'(A,F12.6)') &
      '  Virial ratio                                  ...',-(molE-T)/T
    else if (DKH_order == 2) then
      if (SRTP_type) then
        write(60,'(A,F12.6)') &
        '  Virial ratio                                  ...',&
        -(molE-T-ESOC-ESR)/T
      else
        write(60,'(A,F12.6)') &
        '  Virial ratio                                  ...',-(molE-T-ESOC)/T
      end if
      write(60,*)
      write(60,'(A)') '  Note: relativistic calculation causes the Virial ratio'
      write(60,'(A)') '  to deviate (usually below) 2.0'
      write(60,*)
    end if
    write(60,'(A,F12.6)') &
    '  total alpha electron                          ...',totalpha
    write(60,'(A,F12.6)') &
    '  total alpha electron                          ...',totbeta
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
          if (abs(AO2MO(2*loop_j-1,loop_i)) >= prtlev) then
            write(60,'(A,I2.2,A,I1,A,A2,A,I2.2,A,F9.6,A1,F9.6,A1)') &
            '             --   A #',ishell,' L=',basis_inf(loop_j)%L-1,'     ',&
            element_list(molecular(basis_inf(loop_j)%atom)%atom_number),' #',&
            basis_inf(loop_j)%atom,'    (', real(AO2MO(2*loop_j-1,&
            loop_i)), ',', aimag(AO2MO(2*loop_j-1,loop_i)), ')'
          end if
          if (abs(AO2MO(2*loop_j,loop_i)) >= prtlev) then
            write(60,'(A,I2.2,A,I1,A,A2,A,I2.2,A,F9.6,A1,F9.6,A1)') &
            '             --   B #',ishell,' L=',basis_inf(loop_j)%L-1,'     ',&
            element_list(molecular(basis_inf(loop_j)%atom)%atom_number),' #',&
            basis_inf(loop_j)%atom,'    (', real(AO2MO(2*loop_j,&
            loop_i)), ',', aimag(AO2MO(2*loop_j,loop_i)), ')'
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
          if (abs(AO2MO(2*loop_j-1,loop_i)) >= prtlev) then
            write(60,'(A,I2.2,A,I1,A,A2,A,I2.2,A,F9.6,A1,F9.6,A1)') &
            '             --   A #',ishell,' L=',basis_inf(loop_j)%L-1,'     ',&
            element_list(molecular(basis_inf(loop_j)%atom)%atom_number),' #',&
            basis_inf(loop_j)%atom,'    (', real(AO2MO(2*loop_j-1,&
            loop_i)), ',', aimag(AO2MO(2*loop_j-1,loop_i)), ')'
          end if
          if (abs(AO2MO(2*loop_j,loop_i)) >= prtlev) then
            write(60,'(A,I2.2,A,I1,A,A2,A,I2.2,A,F9.6,A1,F9.6,A1)') &
            '             --   B #',ishell,' L=',basis_inf(loop_j)%L-1,'     ',&
            element_list(molecular(basis_inf(loop_j)%atom)%atom_number),' #',&
            basis_inf(loop_j)%atom,'    (', real(AO2MO(2*loop_j,&
            loop_i)), ',', aimag(AO2MO(2*loop_j,loop_i)), ')'
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
          if (abs(AO2MO(2*loop_j-1,loop_i)) >= prtlev) then
            write(60,'(A,I2.2,A,I1,A,A2,A,I2.2,A,F9.6,A1,F9.6,A1)') &
            '             --   A #',ishell,' L=',basis_inf(loop_j)%L-1,'     ',&
            element_list(molecular(basis_inf(loop_j)%atom)%atom_number),' #',&
            basis_inf(loop_j)%atom,'    (', real(AO2MO(2*loop_j-1,&
            loop_i)), ',', aimag(AO2MO(2*loop_j-1,loop_i)), ')'
          end if
          if (abs(AO2MO(2*loop_j,loop_i)) >= prtlev) then
            write(60,'(A,I2.2,A,I1,A,A2,A,I2.2,A,F9.6,A1,F9.6,A1)') &
            '             --   B #',ishell,' L=',basis_inf(loop_j)%L-1,'     ',&
            element_list(molecular(basis_inf(loop_j)%atom)%atom_number),' #',&
            basis_inf(loop_j)%atom,'    (', real(AO2MO(2*loop_j,&
            loop_i)), ',', aimag(AO2MO(2*loop_j,loop_i)), ')'
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
          if (abs(AO2MO(2*loop_j-1,loop_i)) >= prtlev) then
            write(60,'(A,I2.2,A,I1,A,A2,A,I2.2,A,F9.6,A1,F9.6,A1)') &
            '             --   A #',ishell,' L=',basis_inf(loop_j)%L-1,'     ',&
            element_list(molecular(basis_inf(loop_j)%atom)%atom_number),' #',&
            basis_inf(loop_j)%atom,'    (', real(AO2MO(2*loop_j-1,&
            loop_i)), ',', aimag(AO2MO(2*loop_j-1,loop_i)), ')'
          end if
          if (abs(AO2MO(2*loop_j,loop_i)) >= prtlev) then
            write(60,'(A,I2.2,A,I1,A,A2,A,I2.2,A,F9.6,A1,F9.6,A1)') &
            '             --   B #',ishell,' L=',basis_inf(loop_j)%L-1,'     ',&
            element_list(molecular(basis_inf(loop_j)%atom)%atom_number),' #',&
            basis_inf(loop_j)%atom,'    (', real(AO2MO(2*loop_j,&
            loop_i)), ',', aimag(AO2MO(2*loop_j,loop_i)), ')'
          end if
          ishell = ishell + 1
        end do
      end if
    end do
    write(60,'(A)') &
    '  ============================================================='
    write(60,*)
    if (molden) then
      write(60,'(A)') '  dumping AO2MO to '//trim(address_job)//'.molden'
      call dump_molden()
      write(60,'(A)') '  complete!'
    end if
    ! initialization
    deallocate(exi_T_j)
    deallocate(exi_V_j)
    deallocate(oper4)
    deallocate(oper6)
    deallocate(isupp_ev_f)
    deallocate(Rsd)
    deallocate(DIISmat)
    deallocate(ipiv)
    deallocate(rou_history)
    deallocate(rou_pre)
    deallocate(rou_pre_pre)
    deallocate(Fock2_mic)
    deallocate(swintegral)
    deallocate(Fock1)
    deallocate(Fock2)
    deallocate(i_V_j)
    deallocate(i_j)
    deallocate(i_p2_j)
    deallocate(Xm)
    deallocate(exXm)
    if (s_h) then
      deallocate(c2s)
      deallocate(exc2s)
    end if
    if (keep == 0) then
      deallocate(shell_in_element)
      deallocate(atom_basis)
      deallocate(basis_inf)
      deallocate(molecular)
      deallocate(AO2MO)
      deallocate(rou_m)
      deallocate(Fock)
      deallocate(orbE)
      deallocate(oper3)
    end if
    ndschwarz = .true.
    ini_rou = .true.
    if (DKH_order == 2) then
      deallocate(exSOC)
      deallocate(AO2p2)
      deallocate(evl_p2)
      if (SRTP_type) deallocate(exSR)
    end if
    write(60,'(A)') 'exit module SCF'
    if (kill) then
      call terminate('normal')
    else
      call terminate('keep')
    end if
  end subroutine DKH_Hartree_Fock
  
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
      allocate(rou_m(2*sbdm,2*sbdm))
      Nalpha = (electron_count - (spin_mult - 1)) / 2 + (spin_mult - 1)
      Nbeta = (electron_count - (spin_mult - 1)) / 2
      
      ! read MO coefficient
      if (guess_type == 'gaussian') then
        allocate(AO2MO(2*sbdm,2*fbdm))
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
        '(I)',iostat = ios) mat_dimension
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
        '(I)',iostat = ios) mat_dimension
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
        rou_m = c0
        do ploop_i = 1,sbdm
          do ploop_j = 1,sbdm
            do ploop_k = 1,Nalpha
              rou_m(2*ploop_i-1,2*ploop_j-1) = rou_m(2*ploop_i-1,2*ploop_j-1) +&
              AO2MOalpha(ploop_k,ploop_i)*AO2MOalpha(ploop_k,ploop_j)
            end do
            do ploop_k = 1,Nbeta
              rou_m(2*ploop_i,2*ploop_j) = rou_m(2*ploop_i,2*ploop_j) +&
              AO2MObeta(ploop_k,ploop_i)*AO2MObeta(ploop_k,ploop_j)
            end do
          end do
        end do
        deallocate(AO2MOalpha)
        deallocate(AO2MObeta)
      else if (guess_type == 'read') then
        call load_matrix_cmplx(name='ao2mo', m=AO2MO, dmi=ploop_i, dmj=ploop_j)
        if (ploop_i /= 2*sbdm .or. ploop_j /= 2*fbdm) call terminate(&
        'basis dimension in .ao2mo file mismatch with current job')
        rou_m = c0
        do ploop_i = 1, 2*sbdm
          do ploop_j = 1, 2*sbdm
            do ploop_k = 1, electron_count
              rou_m(ploop_i,ploop_j) = rou_m(ploop_i,ploop_j) + &
              AO2MO(ploop_i,ploop_k)*conjg(AO2MO(ploop_j,ploop_k))
            end do
          end do
        end do
      end if
      rou_pre = rou_m
    else
      rou_m = c0
      do ploop_i = 1, 2*sbdm
        do ploop_j = 1, 2*sbdm
          do ploop_k = 1, electron_count
            rou_m(ploop_i,ploop_j) = rou_m(ploop_i,ploop_j) + &
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
        rou_m = c0
        do ploop_i = 1, 2*sbdm
          do ploop_j = 1, 2*sbdm
            do ploop_k = 1, electron_count
              rou_m(ploop_i,ploop_j) = rou_m(ploop_i,ploop_j) + &
              AO2MO(ploop_i,ploop_k)*conjg(AO2MO(ploop_j,ploop_k))
            end do
          end do
        end do
      end if
      if (loop_i /= 1) then
        maxDP = (real(rou_m(1,1))-real(rou_pre(1,1)))**2 + &
        (aimag(rou_m(1,1))-aimag(rou_pre(1,1)))**2
        RMSDP = 0.0_dp
        do ploop_i = 1, 2*sbdm
          do ploop_j = 1, 2*sbdm
            if ((real(rou_m(ploop_i,ploop_j))-&
            real(rou_pre(ploop_i,ploop_j)))**2 + (aimag(rou_m(ploop_i,ploop_j))&
            -aimag(rou_pre(ploop_i,ploop_j)))**2 > maxDP) then
              maxDP = (real(rou_m(ploop_i,ploop_j))-&
              real(rou_pre(ploop_i,ploop_j)))**2 + &
              (aimag(rou_m(ploop_i,ploop_j))-aimag(rou_pre(ploop_i,ploop_j)))**2
            end if
            RMSDP = RMSDP + &
            (real(rou_m(ploop_i,ploop_j))-real(rou_pre(ploop_i,ploop_j)))**2 + &
            (aimag(rou_m(ploop_i,ploop_j))-aimag(rou_pre(ploop_i,ploop_j)))**2
          end do
        end do
        maxDP = sqrt(maxDP)
        RMSDP = RMSDP / (4*sbdm*sbdm)
        RMSDP = sqrt(RMSDP)
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
    integer :: floop_i, floop_j   ! loop variables for subroutine Fock1e
    allocate(Fock1(2*fbdm,2*fbdm))
    allocate(exi_T_j(2*fbdm,2*fbdm))
    allocate(exi_V_j(2*fbdm,2*fbdm))
    if (DKH_order == 0) then
      Fock1 = c0
      exi_T_j = c0
      exi_V_j = c0
      do floop_i = 1, fbdm
        do floop_j = 1, fbdm
          Fock1(2*floop_i-1,2*floop_j-1) = &
          i_p2_j(floop_i,floop_j)/2.0_dp + i_V_j(floop_i,floop_j)
          Fock1(2*floop_i,2*floop_j) = &
          i_p2_j(floop_i,floop_j)/2.0_dp + i_V_j(floop_i,floop_j)
          exi_T_j(2*floop_i-1,2*floop_j-1) = i_p2_j(floop_i,floop_j)/2.0_dp
          exi_T_j(2*floop_i,2*floop_j) = i_p2_j(floop_i,floop_j)/2.0_dp
          exi_V_j(2*floop_i-1,2*floop_j-1) = i_V_j(floop_i,floop_j)
          exi_V_j(2*floop_i,2*floop_j) = i_V_j(floop_i,floop_j)
        end do
      end do
    else if (DKH_order == 2) then
      allocate(oper1(fbdm,fbdm))
      allocate(oper2(fbdm,fbdm))
      allocate(oper3(2*fbdm,2*fbdm))
      allocate(oper4(2*fbdm,2*fbdm))
      allocate(oper5(2*fbdm,2*fbdm))
      allocate(Ap(fbdm,fbdm))
      allocate(ApRp(fbdm,fbdm))  ! ApRp = RpAp
      allocate(SRp(2*fbdm,2*fbdm))
      allocate(ARVRA(2*fbdm,2*fbdm))
      allocate(AVA(2*fbdm,2*fbdm))
      allocate(ARVeRA(2*fbdm,2*fbdm))
      allocate(AVeA(2*fbdm,2*fbdm))
      allocate(Ve(fbdm,fbdm))
      allocate(pxVepx(fbdm,fbdm))
      allocate(pyVepy(fbdm,fbdm))
      allocate(pzVepz(fbdm,fbdm))
      allocate(pxVepy(fbdm,fbdm))
      allocate(pyVepx(fbdm,fbdm))
      allocate(pxVepz(fbdm,fbdm))
      allocate(pzVepx(fbdm,fbdm))
      allocate(pyVepz(fbdm,fbdm))
      allocate(pzVepy(fbdm,fbdm))
      allocate(exAO2p2(2*fbdm,2*fbdm))
      allocate(exSOC(2*fbdm,2*fbdm))
      ARVRA = c0
      AVA = c0
      ARVeRA = c0
      AVeA = c0
      Fock1 = c0
      exi_T_j = c0
      exi_V_j = c0
      exSOC = c0
      Ap = 0.0_dp
      do floop_i=1,fbdm
        Ap(floop_i,floop_i)= sqrt((sqrt(evl_p2(floop_i)/(speedc*speedc)+1.0_dp)&
        + 1.0_dp)/(2.0_dp*sqrt(evl_p2(floop_i)/(speedc*speedc) + 1.0_dp)))
      end do
      ApRp = 0.0_dp
      do floop_i=1,fbdm
        ApRp(floop_i,floop_i) = Ap(floop_i,floop_i)/(sqrt(evl_p2(floop_i) + &
        speedc*speedc) + speedc)
      end do
      SRp = c0
      do floop_i=1,fbdm
        SRp(2*floop_i-1,2*floop_i-1) = cmplx(1.0_dp + &
        evl_p2(floop_i)/(4.0_dp*speedc*speedc) + QED_rad,0.0_dp,dp)
        SRp(2*floop_i,2*floop_i) = cmplx(1.0_dp + &
        evl_p2(floop_i)/(4.0_dp*speedc*speedc) + QED_rad,0.0_dp,dp)
      end do
      ! Ap V Ap
      call dgemm( 'N', 'T', fbdm, fbdm, fbdm, 1.0_dp, Ap, fbdm, &
      AO2p2, fbdm, 0.0_dp, oper2, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
      i_V_j, fbdm, 0.0_dp, oper1, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper1, fbdm, &
      AO2p2, fbdm, 0.0_dp, oper2, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
      Ap, fbdm, 0.0_dp, oper1, fbdm)
      do floop_i=1, fbdm
        do floop_j=1, fbdm
          AVA(2*floop_i-1,2*floop_j-1) = AVA(2*floop_i-1,2*floop_j-1) + &
          cmplx(oper1(floop_i,floop_j),0.0_dp,dp)
          AVA(2*floop_i,2*floop_j) = AVA(2*floop_i,2*floop_j) + &
          cmplx(oper1(floop_i,floop_j),0.0_dp,dp)
          exi_V_j(2*floop_i-1,2*floop_j-1) = i_V_j(floop_i,floop_j)
          exi_V_j(2*floop_i,2*floop_j) = i_V_j(floop_i,floop_j)
        end do
      end do
      ! ApRp pxVpx+pyVpy+pzVpz ApRp
      call dgemm( 'N', 'T', fbdm, fbdm, fbdm, 1.0_dp, ApRp, fbdm, &
      AO2p2, fbdm, 0.0_dp, oper2, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
      i_pxVpx_j+i_pyVpy_j+i_pzVpz_j, fbdm, 0.0_dp, oper1, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper1, fbdm, &
      AO2p2, fbdm, 0.0_dp, oper2, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
      ApRp, fbdm, 0.0_dp, oper1, fbdm)
      do floop_i=1, fbdm
        do floop_j=1, fbdm
          ARVRA(2*floop_i-1,2*floop_j-1) = ARVRA(2*floop_i-1,2*floop_j-1) + &
          cmplx(oper1(floop_i,floop_j),0.0_dp,dp)
          ARVRA(2*floop_i,2*floop_j) = ARVRA(2*floop_i,2*floop_j) + &
          cmplx(oper1(floop_i,floop_j),0.0_dp,dp)
        end do
      end do
      ! ApRp pxVpy-pyVpx ApRp
      call dgemm( 'N', 'T', fbdm, fbdm, fbdm, 1.0_dp, ApRp, fbdm, &
      AO2p2, fbdm, 0.0_dp, oper2, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
      i_pxVpy_j-i_pyVpx_j, fbdm, 0.0_dp, oper1, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper1, fbdm, &
      AO2p2, fbdm, 0.0_dp, oper2, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
      ApRp, fbdm, 0.0_dp, oper1, fbdm)
      do floop_i=1, fbdm
        do floop_j=1, fbdm
          ARVRA(2*floop_i-1,2*floop_j-1) = ARVRA(2*floop_i-1,2*floop_j-1) + &
          cmplx(0.0_dp,oper1(floop_i,floop_j),dp)
          ARVRA(2*floop_i,2*floop_j) = ARVRA(2*floop_i,2*floop_j) - &
          cmplx(0.0_dp,oper1(floop_i,floop_j),dp)
        end do
      end do
      ! ApRp pzVpx-pxVpz ApRp
      call dgemm( 'N', 'T', fbdm, fbdm, fbdm, 1.0_dp, ApRp, fbdm, &
      AO2p2, fbdm, 0.0_dp, oper2, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
      i_pzVpx_j-i_pxVpz_j, fbdm, 0.0_dp, oper1, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper1, fbdm, &
      AO2p2, fbdm, 0.0_dp, oper2, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
      ApRp, fbdm, 0.0_dp, oper1, fbdm)
      do floop_i=1, fbdm
        do floop_j=1, fbdm
          ARVRA(2*floop_i,2*floop_j-1) = ARVRA(2*floop_i,2*floop_j-1) - &
          cmplx(oper1(floop_i,floop_j),0.0_dp,dp)
          ARVRA(2*floop_i-1,2*floop_j) = ARVRA(2*floop_i-1,2*floop_j) + &
          cmplx(oper1(floop_i,floop_j),0.0_dp,dp)
        end do
      end do
      ! ApRp pyVpz-pzVpy ApRp
      call dgemm( 'N', 'T', fbdm, fbdm, fbdm, 1.0_dp, ApRp, fbdm, &
      AO2p2, fbdm, 0.0_dp, oper2, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
      i_pyVpz_j-i_pzVpy_j, fbdm, 0.0_dp, oper1, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper1, fbdm, &
      AO2p2, fbdm, 0.0_dp, oper2, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
      ApRp, fbdm, 0.0_dp, oper1, fbdm)
      do floop_i=1, fbdm
        do floop_j=1, fbdm
          ARVRA(2*floop_i,2*floop_j-1) = ARVRA(2*floop_i,2*floop_j-1) + &
          cmplx(0.0_dp,oper1(floop_i,floop_j),dp)
          ARVRA(2*floop_i-1,2*floop_j) = ARVRA(2*floop_i-1,2*floop_j) + &
          cmplx(0.0_dp,oper1(floop_i,floop_j),dp)
        end do
      end do
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
        allocate(exSR(2*fbdm,2*fbdm))
        exSR = c0
        ! ApRp px3Vpx+py3Vpy+pz3Vpz ApRp
        call dgemm( 'N', 'T', fbdm, fbdm, fbdm, 1.0_dp, ApRp, fbdm, &
        AO2p2, fbdm, 0.0_dp, oper2, fbdm)
        call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
        i_px3Vpx_j+i_py3Vpy_j+i_pz3Vpz_j, fbdm, 0.0_dp, oper1, fbdm)
        call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper1, fbdm, &
        AO2p2, fbdm, 0.0_dp, oper2, fbdm)
        call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
        ApRp, fbdm, 0.0_dp, oper1, fbdm)
        do floop_i=1, fbdm
          do floop_j=1, fbdm
            Fock1(2*floop_i-1,2*floop_j-1) = Fock1(2*floop_i-1,2*floop_j-1) - &
            cmplx(oper1(floop_i,floop_j) / (2.0_dp*speedc*speedc),0.0_dp,dp)
            Fock1(2*floop_i,2*floop_j) = Fock1(2*floop_i,2*floop_j) - &
            cmplx(oper1(floop_i,floop_j) / (2.0_dp*speedc*speedc),0.0_dp,dp)
            exSR(2*floop_i-1,2*floop_j-1) = exSR(2*floop_i-1,2*floop_j-1) - &
            cmplx(oper1(floop_i,floop_j) / (2.0_dp*speedc*speedc),0.0_dp,dp)
            exSR(2*floop_i,2*floop_j) = exSR(2*floop_i,2*floop_j) - &
            cmplx(oper1(floop_i,floop_j) / (2.0_dp*speedc*speedc),0.0_dp,dp)
          end do
        end do
        ! ApRp px3Vpy-py3Vpx ApRp
        call dgemm( 'N', 'T', fbdm, fbdm, fbdm, 1.0_dp, ApRp, fbdm, &
        AO2p2, fbdm, 0.0_dp, oper2, fbdm)
        call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
        i_px3Vpy_j-i_py3Vpx_j, fbdm, 0.0_dp, oper1, fbdm)
        call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper1, fbdm, &
        AO2p2, fbdm, 0.0_dp, oper2, fbdm)
        call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
        ApRp, fbdm, 0.0_dp, oper1, fbdm)
        do floop_i=1, fbdm
          do floop_j=1, fbdm
            Fock1(2*floop_i-1,2*floop_j-1) = Fock1(2*floop_i-1,2*floop_j-1) - &
            cmplx(0.0_dp,oper1(floop_i,floop_j) / (2.0_dp*speedc*speedc),dp)
            Fock1(2*floop_i,2*floop_j) = Fock1(2*floop_i,2*floop_j) + &
            cmplx(0.0_dp,oper1(floop_i,floop_j) / (2.0_dp*speedc*speedc),dp)
            exSR(2*floop_i-1,2*floop_j-1) = exSR(2*floop_i-1,2*floop_j-1) - &
            cmplx(0.0_dp,oper1(floop_i,floop_j) / (2.0_dp*speedc*speedc),dp)
            exSR(2*floop_i,2*floop_j) = exSR(2*floop_i,2*floop_j) + &
            cmplx(0.0_dp,oper1(floop_i,floop_j) / (2.0_dp*speedc*speedc),dp)
          end do
        end do
        ! ApRp pz3Vpx-px3Vpz ApRp
        call dgemm( 'N', 'T', fbdm, fbdm, fbdm, 1.0_dp, ApRp, fbdm, &
        AO2p2, fbdm, 0.0_dp, oper2, fbdm)
        call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
        i_pz3Vpx_j-i_px3Vpz_j, fbdm, 0.0_dp, oper1, fbdm)
        call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper1, fbdm, &
        AO2p2, fbdm, 0.0_dp, oper2, fbdm)
        call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
        ApRp, fbdm, 0.0_dp, oper1, fbdm)
        do floop_i=1, fbdm
          do floop_j=1, fbdm
            Fock1(2*floop_i,2*floop_j-1) = Fock1(2*floop_i,2*floop_j-1) + &
            cmplx(oper1(floop_i,floop_j) / (2.0_dp*speedc*speedc),0.0_dp,dp)
            Fock1(2*floop_i-1,2*floop_j) = Fock1(2*floop_i-1,2*floop_j) - &
            cmplx(oper1(floop_i,floop_j) / (2.0_dp*speedc*speedc),0.0_dp,dp)
            exSR(2*floop_i,2*floop_j-1) = exSR(2*floop_i,2*floop_j-1) + &
            cmplx(oper1(floop_i,floop_j) / (2.0_dp*speedc*speedc),0.0_dp,dp)
            exSR(2*floop_i-1,2*floop_j) = exSR(2*floop_i-1,2*floop_j) - &
            cmplx(oper1(floop_i,floop_j) / (2.0_dp*speedc*speedc),0.0_dp,dp)
          end do
        end do
        ! ApRp py3Vpz-pz3Vpy ApRp
        call dgemm( 'N', 'T', fbdm, fbdm, fbdm, 1.0_dp, ApRp, fbdm, &
        AO2p2, fbdm, 0.0_dp, oper2, fbdm)
        call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
        i_py3Vpz_j-i_pz3Vpy_j, fbdm, 0.0_dp, oper1, fbdm)
        call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper1, fbdm, &
        AO2p2, fbdm, 0.0_dp, oper2, fbdm)
        call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
        ApRp, fbdm, 0.0_dp, oper1, fbdm)
        do floop_i=1, fbdm
          do floop_j=1, fbdm
            Fock1(2*floop_i,2*floop_j-1) = Fock1(2*floop_i,2*floop_j-1) - &
            cmplx(0.0_dp,oper1(floop_i,floop_j) / (2.0_dp*speedc*speedc),dp)
            Fock1(2*floop_i-1,2*floop_j) = Fock1(2*floop_i-1,2*floop_j) - &
            cmplx(0.0_dp,oper1(floop_i,floop_j) / (2.0_dp*speedc*speedc),dp)
            exSR(2*floop_i,2*floop_j-1) = exSR(2*floop_i,2*floop_j-1) - &
            cmplx(0.0_dp,oper1(floop_i,floop_j) / (2.0_dp*speedc*speedc),dp)
            exSR(2*floop_i-1,2*floop_j) = exSR(2*floop_i-1,2*floop_j) - &
            cmplx(0.0_dp,oper1(floop_i,floop_j) / (2.0_dp*speedc*speedc),dp)
          end do
        end do
        ! ApRp pxVpx3+pyVpy3+pzVpz3 ApRp
        call dgemm( 'N', 'T', fbdm, fbdm, fbdm, 1.0_dp, ApRp, fbdm, &
        AO2p2, fbdm, 0.0_dp, oper2, fbdm)
        call dgemm( 'N', 'T', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
        i_px3Vpx_j+i_py3Vpy_j+i_pz3Vpz_j, fbdm, 0.0_dp, oper1, fbdm)
        call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper1, fbdm, &
        AO2p2, fbdm, 0.0_dp, oper2, fbdm)
        call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
        ApRp, fbdm, 0.0_dp, oper1, fbdm)
        do floop_i=1, fbdm
          do floop_j=1, fbdm
            Fock1(2*floop_i-1,2*floop_j-1) = Fock1(2*floop_i-1,2*floop_j-1) - &
            cmplx(oper1(floop_i,floop_j) / (2.0_dp*speedc*speedc),0.0_dp,dp)
            Fock1(2*floop_i,2*floop_j) = Fock1(2*floop_i,2*floop_j) - &
            cmplx(oper1(floop_i,floop_j) / (2.0_dp*speedc*speedc),0.0_dp,dp)
            exSR(2*floop_i-1,2*floop_j-1) = exSR(2*floop_i-1,2*floop_j-1) - &
            cmplx(oper1(floop_i,floop_j) / (2.0_dp*speedc*speedc),0.0_dp,dp)
            exSR(2*floop_i,2*floop_j) = exSR(2*floop_i,2*floop_j) - &
            cmplx(oper1(floop_i,floop_j) / (2.0_dp*speedc*speedc),0.0_dp,dp)
          end do
        end do
        ! ApRp pxVpy3-pyVpx3 ApRp
        call dgemm( 'N', 'T', fbdm, fbdm, fbdm, 1.0_dp, ApRp, fbdm, &
        AO2p2, fbdm, 0.0_dp, oper2, fbdm)
        call dgemm( 'N', 'T', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
        i_py3Vpx_j-i_px3Vpy_j, fbdm, 0.0_dp, oper1, fbdm)
        call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper1, fbdm, &
        AO2p2, fbdm, 0.0_dp, oper2, fbdm)
        call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
        ApRp, fbdm, 0.0_dp, oper1, fbdm)
        do floop_i=1, fbdm
          do floop_j=1, fbdm
            Fock1(2*floop_i-1,2*floop_j-1) = Fock1(2*floop_i-1,2*floop_j-1) - &
            cmplx(0.0_dp,oper1(floop_i,floop_j) / (2.0_dp*speedc*speedc),dp)
            Fock1(2*floop_i,2*floop_j) = Fock1(2*floop_i,2*floop_j) + &
            cmplx(0.0_dp,oper1(floop_i,floop_j) / (2.0_dp*speedc*speedc),dp)
            exSR(2*floop_i-1,2*floop_j-1) = exSR(2*floop_i-1,2*floop_j-1) - &
            cmplx(0.0_dp,oper1(floop_i,floop_j) / (2.0_dp*speedc*speedc),dp)
            exSR(2*floop_i,2*floop_j) = exSR(2*floop_i,2*floop_j) + &
            cmplx(0.0_dp,oper1(floop_i,floop_j) / (2.0_dp*speedc*speedc),dp)
          end do
        end do
        ! ApRp pxVpz3-pzVpx3 ApRp
        call dgemm( 'N', 'T', fbdm, fbdm, fbdm, 1.0_dp, ApRp, fbdm, &
        AO2p2, fbdm, 0.0_dp, oper2, fbdm)
        call dgemm( 'N', 'T', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
        i_px3Vpz_j-i_pz3Vpx_j, fbdm, 0.0_dp, oper1, fbdm)
        call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper1, fbdm, &
        AO2p2, fbdm, 0.0_dp, oper2, fbdm)
        call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
        ApRp, fbdm, 0.0_dp, oper1, fbdm)
        do floop_i=1, fbdm
          do floop_j=1, fbdm
            Fock1(2*floop_i,2*floop_j-1) = Fock1(2*floop_i,2*floop_j-1) + &
            cmplx(oper1(floop_i,floop_j) / (2.0_dp*speedc*speedc),0.0_dp,dp)
            Fock1(2*floop_i-1,2*floop_j) = Fock1(2*floop_i-1,2*floop_j) - &
            cmplx(oper1(floop_i,floop_j) / (2.0_dp*speedc*speedc),0.0_dp,dp)
            exSR(2*floop_i,2*floop_j-1) = exSR(2*floop_i,2*floop_j-1) + &
            cmplx(oper1(floop_i,floop_j) / (2.0_dp*speedc*speedc),0.0_dp,dp)
            exSR(2*floop_i-1,2*floop_j) = exSR(2*floop_i-1,2*floop_j) - &
            cmplx(oper1(floop_i,floop_j) / (2.0_dp*speedc*speedc),0.0_dp,dp)
          end do
        end do
        ! ApRp pyVpz3-pzVpy3 ApRp
        call dgemm( 'N', 'T', fbdm, fbdm, fbdm, 1.0_dp, ApRp, fbdm, &
        AO2p2, fbdm, 0.0_dp, oper2, fbdm)
        call dgemm( 'N', 'T', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
        i_pz3Vpy_j-i_py3Vpz_j, fbdm, 0.0_dp, oper1, fbdm)
        call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper1, fbdm, &
        AO2p2, fbdm, 0.0_dp, oper2, fbdm)
        call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
        ApRp, fbdm, 0.0_dp, oper1, fbdm)
        do floop_i=1, fbdm
          do floop_j=1, fbdm
            Fock1(2*floop_i,2*floop_j-1) = Fock1(2*floop_i,2*floop_j-1) - &
            cmplx(0.0_dp,oper1(floop_i,floop_j) / (2.0_dp*speedc*speedc),dp)
            Fock1(2*floop_i-1,2*floop_j) = Fock1(2*floop_i-1,2*floop_j) - &
            cmplx(0.0_dp,oper1(floop_i,floop_j) / (2.0_dp*speedc*speedc),dp)
            exSR(2*floop_i,2*floop_j-1) = exSR(2*floop_i,2*floop_j-1) - &
            cmplx(0.0_dp,oper1(floop_i,floop_j) / (2.0_dp*speedc*speedc),dp)
            exSR(2*floop_i-1,2*floop_j) = exSR(2*floop_i-1,2*floop_j) - &
            cmplx(0.0_dp,oper1(floop_i,floop_j) / (2.0_dp*speedc*speedc),dp)
          end do
        end do
      end if
      do floop_i=1,fbdm
        do floop_j=1,fbdm
          Ve(floop_i,floop_j) = i_V_j(floop_i,floop_j)/&
          (speedc * (sqrt(evl_p2(floop_i)+speedc*speedc) + &
          sqrt(evl_p2(floop_j)+speedc*speedc)))
        end do
      end do
      do floop_i=1,fbdm
        do floop_j=1,fbdm
          pxVepx(floop_i,floop_j) = i_pxVpx_j(floop_i,floop_j)/&
          (speedc * (sqrt(evl_p2(floop_i)+speedc*speedc) + &
          sqrt(evl_p2(floop_j)+speedc*speedc)))
        end do
      end do
      do floop_i=1,fbdm
        do floop_j=1,fbdm
          pyVepy(floop_i,floop_j) = i_pyVpy_j(floop_i,floop_j)/&
          (speedc * (sqrt(evl_p2(floop_i)+speedc*speedc) + &
          sqrt(evl_p2(floop_j)+speedc*speedc)))
        end do
      end do
      do floop_i=1,fbdm
        do floop_j=1,fbdm
          pzVepz(floop_i,floop_j) = i_pzVpz_j(floop_i,floop_j)/&
          (speedc * (sqrt(evl_p2(floop_i)+speedc*speedc) + &
          sqrt(evl_p2(floop_j)+speedc*speedc)))
        end do
      end do
      do floop_i=1,fbdm
        do floop_j=1,fbdm
          pxVepy(floop_i,floop_j) = i_pxVpy_j(floop_i,floop_j)/&
          (speedc * (sqrt(evl_p2(floop_i)+speedc*speedc) + &
          sqrt(evl_p2(floop_j)+speedc*speedc)))
        end do
      end do
      do floop_i=1,fbdm
        do floop_j=1,fbdm
          pyVepx(floop_i,floop_j) = i_pyVpx_j(floop_i,floop_j)/&
          (speedc * (sqrt(evl_p2(floop_i)+speedc*speedc) + &
          sqrt(evl_p2(floop_j)+speedc*speedc)))
        end do
      end do
      do floop_i=1,fbdm
        do floop_j=1,fbdm
          pxVepz(floop_i,floop_j) = i_pxVpz_j(floop_i,floop_j)/&
          (speedc * (sqrt(evl_p2(floop_i)+speedc*speedc) + &
          sqrt(evl_p2(floop_j)+speedc*speedc)))
        end do
      end do
      do floop_i=1,fbdm
        do floop_j=1,fbdm
          pzVepx(floop_i,floop_j) = i_pzVpx_j(floop_i,floop_j)/&
          (speedc * (sqrt(evl_p2(floop_i)+speedc*speedc) + &
          sqrt(evl_p2(floop_j)+speedc*speedc)))
        end do
      end do
      do floop_i=1,fbdm
        do floop_j=1,fbdm
          pyVepz(floop_i,floop_j) = i_pyVpz_j(floop_i,floop_j)/&
          (speedc * (sqrt(evl_p2(floop_i)+speedc*speedc) + &
          sqrt(evl_p2(floop_j)+speedc*speedc)))
        end do
      end do
      do floop_i=1,fbdm
        do floop_j=1,fbdm
          pzVepy(floop_i,floop_j) = i_pzVpy_j(floop_i,floop_j)/&
          (speedc * (sqrt(evl_p2(floop_i)+speedc*speedc) + &
          sqrt(evl_p2(floop_j)+speedc*speedc)))
        end do
      end do
      ! Ap Ve Ap
      call dgemm( 'N', 'T', fbdm, fbdm, fbdm, 1.0_dp, Ap, fbdm, &
      AO2p2, fbdm, 0.0_dp, oper2, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
      Ve, fbdm, 0.0_dp, oper1, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper1, fbdm, &
      AO2p2, fbdm, 0.0_dp, oper2, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
      Ap, fbdm, 0.0_dp, oper1, fbdm)
      do floop_i=1, fbdm
        do floop_j=1, fbdm
          AVeA(2*floop_i-1,2*floop_j-1) = AVeA(2*floop_i-1,2*floop_j-1) + &
          cmplx(oper1(floop_i,floop_j),0.0_dp,dp)
          AVeA(2*floop_i,2*floop_j) = AVeA(2*floop_i,2*floop_j) + &
          cmplx(oper1(floop_i,floop_j),0.0_dp,dp)
        end do
      end do
      ! ApRp pxVepx+pyVepy+pzVepz ApRp
      call dgemm( 'N', 'T', fbdm, fbdm, fbdm, 1.0_dp, ApRp, fbdm, &
      AO2p2, fbdm, 0.0_dp, oper2, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
      pxVepx+pyVepy+pzVepz, fbdm, 0.0_dp, oper1, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper1, fbdm, &
      AO2p2, fbdm, 0.0_dp, oper2, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
      ApRp, fbdm, 0.0_dp, oper1, fbdm)
      do floop_i=1, fbdm
        do floop_j=1, fbdm
          ARVeRA(2*floop_i-1,2*floop_j-1) = ARVeRA(2*floop_i-1,2*floop_j-1) + &
          cmplx(oper1(floop_i,floop_j),0.0_dp,dp)
          ARVeRA(2*floop_i,2*floop_j) = ARVeRA(2*floop_i,2*floop_j) + &
          cmplx(oper1(floop_i,floop_j),0.0_dp,dp)
        end do
      end do
      ! ApRp pxVepy-pyVepx ApRp
      call dgemm( 'N', 'T', fbdm, fbdm, fbdm, 1.0_dp, ApRp, fbdm, &
      AO2p2, fbdm, 0.0_dp, oper2, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
      pxVepy-pyVepx, fbdm, 0.0_dp, oper1, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper1, fbdm, &
      AO2p2, fbdm, 0.0_dp, oper2, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
      ApRp, fbdm, 0.0_dp, oper1, fbdm)
      do floop_i=1, fbdm
        do floop_j=1, fbdm
          ARVeRA(2*floop_i-1,2*floop_j-1) = ARVeRA(2*floop_i-1,2*floop_j-1) + &
          cmplx(0.0_dp,oper1(floop_i,floop_j),dp)
          ARVeRA(2*floop_i,2*floop_j) = ARVeRA(2*floop_i,2*floop_j) - &
          cmplx(0.0_dp,oper1(floop_i,floop_j),dp)
        end do
      end do
      ! ApRp pzVepx-pxVepz ApRp
      call dgemm( 'N', 'T', fbdm, fbdm, fbdm, 1.0_dp, ApRp, fbdm, &
      AO2p2, fbdm, 0.0_dp, oper2, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
      pzVepx-pxVepz, fbdm, 0.0_dp, oper1, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper1, fbdm, &
      AO2p2, fbdm, 0.0_dp, oper2, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
      ApRp, fbdm, 0.0_dp, oper1, fbdm)
      do floop_i=1, fbdm
        do floop_j=1, fbdm
          ARVeRA(2*floop_i,2*floop_j-1) = ARVeRA(2*floop_i,2*floop_j-1) - &
          cmplx(oper1(floop_i,floop_j),0.0_dp,dp)
          ARVeRA(2*floop_i-1,2*floop_j) = ARVeRA(2*floop_i-1,2*floop_j) + &
          cmplx(oper1(floop_i,floop_j),0.0_dp,dp)
        end do
      end do
      ! ApRp pyVepz-pzVepy ApRp
      call dgemm( 'N', 'T', fbdm, fbdm, fbdm, 1.0_dp, ApRp, fbdm, &
      AO2p2, fbdm, 0.0_dp, oper2, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
      pyVepz-pzVepy, fbdm, 0.0_dp, oper1, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper1, fbdm, &
      AO2p2, fbdm, 0.0_dp, oper2, fbdm)
      call dgemm( 'N', 'N', fbdm, fbdm, fbdm, 1.0_dp, oper2, fbdm, &
      ApRp, fbdm, 0.0_dp, oper1, fbdm)
      do floop_i=1, fbdm
        do floop_j=1, fbdm
          ARVeRA(2*floop_i,2*floop_j-1) = ARVeRA(2*floop_i,2*floop_j-1) + &
          cmplx(0.0_dp,oper1(floop_i,floop_j),dp)
          ARVeRA(2*floop_i-1,2*floop_j) = ARVeRA(2*floop_i-1,2*floop_j) + &
          cmplx(0.0_dp,oper1(floop_i,floop_j),dp)
        end do
      end do
    
      oper5 = c0
      do floop_i=1, fbdm
        oper5(2*floop_i-1,2*floop_i-1) = cmplx(0.5_dp,0.0_dp,dp)
        oper5(2*floop_i,2*floop_i) = cmplx(0.5_dp,0.0_dp,dp)
      end do
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
      do floop_i=1, fbdm
        oper5(2*floop_i-1,2*floop_i-1) = cmplx(0.5_dp,0.0_dp,dp) * &
        cmplx((sqrt(evl_p2(floop_i) + speedc*speedc) + &
        speedc)**2 / evl_p2(floop_i),0.0_dp,dp)
        oper5(2*floop_i,2*floop_i) = cmplx(0.5_dp,0.0_dp,dp) * &
        cmplx((sqrt(evl_p2(floop_i) + speedc*speedc) + &
        speedc)**2 / evl_p2(floop_i),0.0_dp,dp)
      end do
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
      do floop_i=1, fbdm
        oper5(2*floop_i-1,2*floop_i-1) = cmplx(evl_p2(floop_i) / &
        (2.0_dp*(sqrt(evl_p2(floop_i)+speedc*speedc) + speedc)**2),0.0_dp,dp)
        oper5(2*floop_i,2*floop_i) = cmplx(evl_p2(floop_i) / &
        (2.0_dp * speedc * sqrt(evl_p2(floop_i) + speedc*speedc)),0.0_dp,dp)
      end do
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
      do floop_i=1,fbdm
        Fock1(2*floop_i-1,2*floop_i-1) = &
        Fock1(2*floop_i-1,2*floop_i-1) + cmplx(speedc*sqrt(evl_p2(floop_i) + &
        speedc*speedc),0.0_dp,dp) - speedc*speedc
        Fock1(2*floop_i,2*floop_i) = &
        Fock1(2*floop_i,2*floop_i) + cmplx(speedc*sqrt(evl_p2(floop_i) + &
        speedc*speedc),0.0_dp,dp) - speedc*speedc
        exi_T_j(2*floop_i-1,2*floop_i-1) = cmplx(speedc*sqrt(evl_p2(floop_i) + &
        speedc*speedc),0.0_dp,dp) - speedc*speedc
        exi_T_j(2*floop_i,2*floop_i) = cmplx(speedc*sqrt(evl_p2(floop_i) + &
        speedc*speedc),0.0_dp,dp) - speedc*speedc
      end do
      ! Transform from p^2 eigenbasis to orthogonal normalized AO basis
      exAO2p2 = c0
      do floop_i=1,fbdm
        do floop_j=1,fbdm
          exAO2p2(2*floop_i-1,2*floop_j-1) = &
          cmplx(AO2p2(floop_i,floop_j),0.0_dp,dp)
          exAO2p2(2*floop_i,2*floop_j) = &
          cmplx(AO2p2(floop_i,floop_j),0.0_dp,dp)
        end do
      end do
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
      deallocate(Ve)
      deallocate(pxVepx)
      deallocate(pyVepy)
      deallocate(pzVepz)
      deallocate(pxVepy)
      deallocate(pyVepx)
      deallocate(pxVepz)
      deallocate(pzVepx)
      deallocate(pyVepz)
      deallocate(pzVepy)
      deallocate(ARVeRA)
      deallocate(AVeA)
      deallocate(oper1)
      deallocate(oper2)
      deallocate(oper3)
      deallocate(oper4)
      deallocate(oper5)
      deallocate(Ap)
      deallocate(ApRp)
      deallocate(SRp)
      deallocate(ARVRA)
      deallocate(AVA)
      deallocate(exAO2p2)
      ! deallocate one-electron integrals to reduce memory usage
      deallocate(i_pxVpx_j)
      deallocate(i_pyVpy_j)
      deallocate(i_pzVpz_j)
      deallocate(i_pxVpy_j)
      deallocate(i_pyVpx_j)
      deallocate(i_pxVpz_j)
      deallocate(i_pzVpx_j)
      deallocate(i_pyVpz_j)
      deallocate(i_pzVpy_j)
      if (SRTP_type) then
        deallocate(i_px3Vpx_j)
        deallocate(i_py3Vpy_j)
        deallocate(i_pz3Vpz_j)
        deallocate(i_px3Vpy_j)
        deallocate(i_py3Vpx_j)
        deallocate(i_px3Vpz_j)
        deallocate(i_pz3Vpx_j)
        deallocate(i_py3Vpz_j)
        deallocate(i_pz3Vpy_j)
      end if
    end if
  end subroutine Fock1e
  
!------------------------------------------------------------
!> construct two electron Fock operator
!!
!! "direct" calculation to avoid memory overflow and large amount of disk r&w
!!
!! but 2 electron integral should be calculated in each SCF iteration
  subroutine Fock2e()
    implicit none
    integer :: i                   ! for parallel computation, dloop_i = i
    integer :: dloop_i, dloop_j    ! loop variables for Fock2e routine
    integer :: dloop_k, dloop_l
    integer :: dloop_m, dloop_n
    integer :: dloop_o, dloop_p
    real(dp) :: integral
    !----------------------------------
    integer :: contraction_i       ! contraction of atom_i, shell_i
    integer :: atom_i              ! which atom is the i^th component of |AOi>
    integer :: shell_i             ! which shell is the i^th component of |AOi>
    integer :: L_i                 ! angular quantum number of |AOi>
    integer :: M_i                 ! magnetic quantum number of |AOi>
    !----------------------------------
    integer :: contraction_j       ! contraction of atom_j, shell_j
    integer :: atom_j              ! which atom is the j_th component of |AOj>
    integer :: shell_j             ! which shell is the j_th component of |AOj>
    integer :: L_j                 ! angular quantum number of |AOj>
    integer :: M_j                 ! magnetic quantum number of |AOj>
    !----------------------------------
    integer :: contraction_k       ! contraction of atom_k, shell_k
    integer :: L_k                 ! angular quantum number of |AOk>
    integer :: M_k                 ! magnetic quantum number of |AOk>
    !----------------------------------
    integer :: contraction_l       ! contraction of atom_l, shell_l
    integer :: L_l                 ! angular quantum number of |AOl>
    integer :: M_l                 ! magnetic quantum number of |AOl>
    real(dp) :: exponents_i(20)    ! exponents of |AOi>
    real(dp) :: exponents_j(20)    ! exponents of |AOj>
    real(dp) :: exponents_k(20)    ! exponents of |AOk>
    real(dp) :: exponents_l(20)    ! exponents of |AOl>
    real(dp) :: coefficient_i(20)  ! coefficient of |AOi>
    real(dp) :: coefficient_j(20)  ! coefficient of |AOj>
    real(dp) :: coefficient_k(20)  ! coefficient of |AOk>
    real(dp) :: coefficient_l(20)  ! coefficient of |AOl>
    real(dp) :: coordinate_i(3)    ! coordinate of center of |AOi>
    real(dp) :: coordinate_j(3)    ! coordinate of center of |AOj>
    real(dp) :: coordinate_k(3)    ! coordinate of center of |AOk>
    real(dp) :: coordinate_l(3)    ! coordinate of center of |AOl>
    real(dp) :: swintegral_mic
    integer :: swloop_i, swloop_j  ! loop variables only for schwarz screening
    integer :: swloop_m, swloop_n
    integer :: swloop_o, swloop_p
    integer :: shell_start_i                    ! start point of an shell
    integer :: shell_start_j                    ! start point of an shell
    real(dp),allocatable :: swexponents_i(:)    ! expo of |AOi> for schwarz
    real(dp),allocatable :: swexponents_j(:)    ! expo of |AOj> for schwarz
    real(dp),allocatable :: swcoefficient_i(:)  ! coeff of |AOi> for schwarz
    real(dp),allocatable :: swcoefficient_j(:)  ! coeff of |AOj> for schwarz
    if (ndschwarz) then
      ndschwarz = .false.
      write(60,'(a)') '  -- Schwarz screening of <ij||kl>'
      allocate(Fock2(2*fbdm,2*fbdm))
      allocate(swintegral(cbdm,cbdm))
      swintegral = 0.0_dp
      swloop_i = 1
      atom_i = 1
      shell_i = 1
      shell_start_i = 1
      do while(swloop_i <= cbdm)
        if (shell_i > shell_in_element(molecular(atom_i) % atom_number)) then
          shell_i = 1
          atom_i = atom_i + 1
        end if
        contraction_i = atom_basis(molecular(atom_i) % &
        basis_number + shell_i - 1) % contraction
        L_i = atom_basis(molecular(atom_i) % basis_number + shell_i - 1) % &
        angular_quantum_number + 1
        M_i = swloop_i - shell_start_i + 1
        allocate(swexponents_i(contraction_i))
        allocate(swcoefficient_i(contraction_i))
        swexponents_i = &
        atom_basis(molecular(atom_i) % basis_number + shell_i - 1) % exponents
        do swloop_m=1,contraction_i
          swcoefficient_i(swloop_m) = atom_basis(molecular(atom_i) % &
          basis_number + shell_i - 1) % Ncoefficient(swloop_m,M_i)
        end do
        coordinate_i = molecular(atom_i) % nucleus_position
        swloop_j = swloop_i
        atom_j = atom_i
        shell_j = shell_i
        shell_start_j = shell_start_i
        do while(swloop_j <= cbdm)
          if (shell_j > shell_in_element(molecular(atom_j) % atom_number)) then
            shell_j = 1
            atom_j = atom_j + 1
          end if
          contraction_j = &
          atom_basis(molecular(atom_j) % basis_number+shell_j-1) % contraction
          L_j = atom_basis(molecular(atom_j) % basis_number + shell_j - 1) % &
          angular_quantum_number + 1
          M_j = swloop_j - shell_start_j + 1
          allocate(swexponents_j(contraction_j))
          allocate(swcoefficient_j(contraction_j))
          swexponents_j = &
          atom_basis(molecular(atom_j) % basis_number + shell_j - 1) % exponents
          do swloop_m=1,contraction_j
            swcoefficient_j(swloop_m) = atom_basis(molecular(atom_j) % &
            basis_number + shell_j - 1) % Ncoefficient(swloop_m,M_j)
          end do
          coordinate_j = molecular(atom_j) % nucleus_position
          do swloop_m = 1, contraction_i
            do swloop_n = 1, contraction_j
              do swloop_o = 1, contraction_i
                do swloop_p = 1, contraction_j
                  if (L_i >= L_j) then
                    swintegral_mic = swcoefficient_i(swloop_m) * &
                    swcoefficient_j(swloop_n) * swcoefficient_i(swloop_o) * &
                    swcoefficient_j(swloop_p) * V_Integral_2e(&
                    AO_xyz_factor(L_i,M_i), AO_xyz_factor(L_j,M_j), &
                    AO_xyz_factor(L_i,M_i), AO_xyz_factor(L_j,M_j), &
                    swexponents_i(swloop_m),swexponents_j(swloop_n),&
                    swexponents_i(swloop_o),swexponents_j(swloop_p),&
                    coordinate_i(1),coordinate_i(2),coordinate_i(3),&
                    coordinate_j(1),coordinate_j(2),coordinate_j(3),&
                    coordinate_i(1),coordinate_i(2),coordinate_i(3),&
                    coordinate_j(1),coordinate_j(2),coordinate_j(3))
                  else
                    swintegral_mic = swcoefficient_i(swloop_m) * &
                    swcoefficient_j(swloop_n) * swcoefficient_i(swloop_o) * &
                    swcoefficient_j(swloop_p) * V_Integral_2e(&
                    AO_xyz_factor(L_j,M_j), AO_xyz_factor(L_i,M_i), &
                    AO_xyz_factor(L_j,M_j), AO_xyz_factor(L_i,M_i), &
                    swexponents_j(swloop_n),swexponents_i(swloop_m),&
                    swexponents_j(swloop_p),swexponents_i(swloop_o),&
                    coordinate_j(1),coordinate_j(2),coordinate_j(3),&
                    coordinate_i(1),coordinate_i(2),coordinate_i(3),&
                    coordinate_j(1),coordinate_j(2),coordinate_j(3),&
                    coordinate_i(1),coordinate_i(2),coordinate_i(3))
                  end if
                  swintegral(swloop_i,swloop_j) = &
                  swintegral(swloop_i,swloop_j) + swintegral_mic
                end do
              end do
            end do
          end do
          swintegral(swloop_j,swloop_i) = swintegral(swloop_i,swloop_j)
          deallocate(swexponents_j)
          deallocate(swcoefficient_j)
          swloop_j = swloop_j + 1
          if (swloop_j - shell_start_j >= (L_j + 1) * L_j / 2) then
            shell_j = shell_j + 1
            shell_start_j = swloop_j
          end if
        end do
        deallocate(swexponents_i)
        deallocate(swcoefficient_i)
        swloop_i = swloop_i + 1
        if (swloop_i - shell_start_i >= (L_i + 1) * L_i / 2) then
          shell_i = shell_i + 1
          shell_start_i = swloop_i
        end if
      end do
      write(60,'(A,E10.2,A)') &
      '  -- complete! cutoff:',schwarz_VT,'; stored in swintegral'
      allocate(Fock2_mic(2*cbdm,2*cbdm))
    end if
    if (s_h) then
      allocate(supp1(2*sbdm,2*cbdm))
      call zgemm( 'N', 'C', 2*sbdm, 2*cbdm, 2*sbdm, c1, rou_m, &
      2*sbdm, exc2s, 2*cbdm, c0, supp1, 2*sbdm)
      deallocate(rou_m)
      allocate(rou_m(2*cbdm,2*cbdm))
      call zgemm( 'N', 'N', 2*cbdm, 2*cbdm, 2*sbdm, c1, exc2s, &
      2*cbdm, supp1, 2*sbdm, c0, rou_m, 2*cbdm)
      deallocate(supp1)
    end if
    allocate(supp1(2*cbdm,2*cbdm))
    supp1 = c0
    ! parallel zone, running results consistent with serial
    !$omp parallel num_threads(threads_use) default(shared) private(i,dloop_i,&
    !$omp& dloop_j,dloop_k,dloop_l,dloop_m,dloop_n,dloop_o,dloop_p,integral,&
    !$omp& contraction_i,L_i,M_i,contraction_j,L_j,M_j,contraction_k,L_k,M_k,&
    !$omp& contraction_l,L_l,M_l,exponents_i,exponents_j,exponents_k,&
    !$omp& exponents_l,coefficient_i,coefficient_j,coefficient_k,coefficient_l,&
    !$omp& coordinate_i,coordinate_j,coordinate_k,coordinate_l,Fock2_mic,&
    !$omp& Fock2_assigned,assigned,iassign,carry) if(threads_use < cpu_threads)
    !$omp do schedule(dynamic,5)
    ! (ij|kl)=(ji|kl)=(ij|lk)=(ji|lk)=(kl|ij)=(kl|ji)=(lk|ij)=(lk|ji)
    !----------------------------<dloop_i>-------------------------------------
    do i = cbdm, 1, -1
      Fock2_mic = c0
      dloop_i = i
      contraction_i = atom_basis(molecular(basis_inf(dloop_i) % atom) % &
      basis_number + basis_inf(dloop_i) % shell - 1) % contraction
      L_i = basis_inf(dloop_i) % L
      M_i = basis_inf(dloop_i) % M
      do dloop_m = 1, contraction_i
        exponents_i(dloop_m) = atom_basis(molecular(basis_inf(dloop_i) % atom) &
        % basis_number + basis_inf(dloop_i) % shell - 1) % exponents(dloop_m)
        coefficient_i(dloop_m) = atom_basis(molecular(basis_inf(dloop_i)%atom) &
        % basis_number + basis_inf(dloop_i)%shell-1) % Ncoefficient(dloop_m,M_i)
      end do
      coordinate_i = molecular(basis_inf(dloop_i) % atom) % nucleus_position
      !----------------------------<dloop_k>-----------------------------------
      do dloop_k = cbdm, 1, -1
        contraction_k = atom_basis(molecular(basis_inf(dloop_k) % atom) % &
        basis_number + basis_inf(dloop_k) % shell - 1) % contraction
        L_k = basis_inf(dloop_k) % L
        M_k = basis_inf(dloop_k) % M
        do dloop_m = 1, contraction_k
          exponents_k(dloop_m) = atom_basis(molecular(basis_inf(dloop_k) % &
          atom) % basis_number + basis_inf(dloop_k) % shell - 1) &
          % exponents(dloop_m)
          coefficient_k(dloop_m) = atom_basis(molecular(basis_inf(dloop_k) % &
          atom) % basis_number + basis_inf(dloop_k) % shell - 1) &
          % Ncoefficient(dloop_m,M_k)
        end do
        coordinate_k = molecular(basis_inf(dloop_k) % atom) % nucleus_position
        !----------------------------<dloop_j>---------------------------------
        do dloop_j = dloop_i, 1, -1
          contraction_j = atom_basis(molecular(basis_inf(dloop_j) % atom) % &
          basis_number + basis_inf(dloop_j) % shell - 1) % contraction
          L_j = basis_inf(dloop_j) % L
          M_j = basis_inf(dloop_j) % M
          do dloop_m = 1, contraction_j
            exponents_j(dloop_m) = atom_basis(molecular(basis_inf(dloop_j) % &
            atom) % basis_number + basis_inf(dloop_j) % shell - 1) % &
            exponents(dloop_m)
            coefficient_j(dloop_m) = atom_basis(molecular(basis_inf(dloop_j) % &
            atom) % basis_number + basis_inf(dloop_j) % shell - 1) % &
            Ncoefficient(dloop_m,M_j)
          end do
          coordinate_j = molecular(basis_inf(dloop_j) % atom) % nucleus_position
          !----------------------------<dloop_l>-------------------------------
          do dloop_l = min(dloop_k,dloop_j+&
          (dloop_i*(dloop_i-1)-dloop_k*(dloop_k-1))/2), 1, -1
            ! Schwarz screening, |<ij||kl>| <= sqrt(<ij||ij>) * sqrt(<kl||kl>)
            if (sqrt(swintegral(dloop_i,dloop_j)*swintegral(dloop_k,dloop_l))&
            < schwarz_VT) cycle
            contraction_l = atom_basis(molecular(basis_inf(dloop_l) % atom) % &
            basis_number + basis_inf(dloop_l) % shell - 1) % contraction
            L_l = basis_inf(dloop_l) % L
            M_l = basis_inf(dloop_l) % M
            do dloop_m = 1, contraction_l
              exponents_l(dloop_m) = atom_basis(molecular(basis_inf(dloop_l) % &
              atom) % basis_number + basis_inf(dloop_l) % shell - 1) &
              % exponents(dloop_m)
              coefficient_l(dloop_m) = atom_basis(molecular(basis_inf(dloop_l)%&
              atom) % basis_number + basis_inf(dloop_l) % shell - 1) &
              % Ncoefficient(dloop_m,M_l)
            end do
            coordinate_l = molecular(basis_inf(dloop_l)%atom) % nucleus_position
            integral = 0.0_dp
            !===========================<ij||kl>===============================
            do dloop_m = 1, contraction_i
              do dloop_n = 1, contraction_j
                do dloop_o = 1, contraction_k
                  do dloop_p = 1, contraction_l
                    if (L_i >= L_j .and. L_k >= L_l) then
                      integral = integral + coefficient_i(dloop_m) * &
                      coefficient_j(dloop_n) * coefficient_k(dloop_o) * &
                      coefficient_l(dloop_p) * V_Integral_2e(&
                      AO_xyz_factor(L_i,M_i), AO_xyz_factor(L_j,M_j),&
                      AO_xyz_factor(L_k,M_k), AO_xyz_factor(L_l,M_l),&
                      exponents_i(dloop_m),exponents_j(dloop_n),&
                      exponents_k(dloop_o),exponents_l(dloop_p),&
                      coordinate_i(1),coordinate_i(2),coordinate_i(3),&
                      coordinate_j(1),coordinate_j(2),coordinate_j(3),&
                      coordinate_k(1),coordinate_k(2),coordinate_k(3),&
                      coordinate_l(1),coordinate_l(2),coordinate_l(3))
                    else if (L_i >= L_j .and. L_k < L_l) then
                      integral = integral + coefficient_i(dloop_m) * &
                      coefficient_j(dloop_n) * coefficient_k(dloop_o) * &
                      coefficient_l(dloop_p) * V_Integral_2e(&
                      AO_xyz_factor(L_i,M_i), AO_xyz_factor(L_j,M_j),&
                      AO_xyz_factor(L_l,M_l), AO_xyz_factor(L_k,M_k),&
                      exponents_i(dloop_m),exponents_j(dloop_n),&
                      exponents_l(dloop_p),exponents_k(dloop_o),&
                      coordinate_i(1),coordinate_i(2),coordinate_i(3),&
                      coordinate_j(1),coordinate_j(2),coordinate_j(3),&
                      coordinate_l(1),coordinate_l(2),coordinate_l(3),&
                      coordinate_k(1),coordinate_k(2),coordinate_k(3))
                    else if (L_i < L_j .and. L_k >= L_l) then
                      integral = integral + coefficient_i(dloop_m) * &
                      coefficient_j(dloop_n) * coefficient_k(dloop_o) * &
                      coefficient_l(dloop_p) * V_Integral_2e(&
                      AO_xyz_factor(L_j,M_j), AO_xyz_factor(L_i,M_i),&
                      AO_xyz_factor(L_k,M_k), AO_xyz_factor(L_l,M_l),&
                      exponents_j(dloop_n),exponents_i(dloop_m),&
                      exponents_k(dloop_o),exponents_l(dloop_p),&
                      coordinate_j(1),coordinate_j(2),coordinate_j(3),&
                      coordinate_i(1),coordinate_i(2),coordinate_i(3),&
                      coordinate_k(1),coordinate_k(2),coordinate_k(3),&
                      coordinate_l(1),coordinate_l(2),coordinate_l(3))
                    else
                      integral = integral + coefficient_i(dloop_m) *&
                       coefficient_j(dloop_n) * coefficient_k(dloop_o) * &
                       coefficient_l(dloop_p) * V_Integral_2e(&
                       AO_xyz_factor(L_j,M_j), AO_xyz_factor(L_i,M_i),&
                       AO_xyz_factor(L_l,M_l), AO_xyz_factor(L_k,M_k),&
                      exponents_j(dloop_n),exponents_i(dloop_m),&
                      exponents_l(dloop_p),exponents_k(dloop_o),&
                      coordinate_j(1),coordinate_j(2),coordinate_j(3),&
                      coordinate_i(1),coordinate_i(2),coordinate_i(3),&
                      coordinate_l(1),coordinate_l(2),coordinate_l(3),&
                      coordinate_k(1),coordinate_k(2),coordinate_k(3))
                    end if
                    !----------<DEBUG>----------
                    !if (dloop_i == 79 .and. dloop_j == 34 .and. dloop_k == 52 &
                    !.and. dloop_l == 34) then
                    !  write(60,*) 'integral',integral
                    !  write(60,*) &
                    !  AO_xyz_factor(L_i,M_i), AO_xyz_factor(L_j,M_j), &
                    !  AO_xyz_factor(L_k,M_k), AO_xyz_factor(L_l,M_l)
                    !  write(60,*) &
                    !  exponents_i(dloop_m),exponents_j(dloop_n),&
                    !  exponents_k(dloop_o),exponents_l(dloop_p)
                    !  write(60,*) &
                    !  coordinate_i,coordinate_j,&
                    !  coordinate_k,coordinate_l
                    !  write(60,*)
                    !  write(60,*) &
                    !  coefficient_i(dloop_m), coefficient_j(dloop_n), &
                    !  coefficient_k(dloop_o), coefficient_l(dloop_p)
                    !end if
                    !----------<DEBUG>----------
                  end do
                end do
              end do
            end do
            !----------<DEBUG>----------
            !if (inin) write(60,'(I3,I3,I3,I3,F)') &
            !dloop_i, dloop_j, dloop_k, dloop_l,integral
            !----------<DEBUG>----------

            ! assign values to two-electron Fock matrix
            ! ref page 261 of Quantum Chemistry: Basic Principles and
            ! Ab-initio Calculations, Volume 2 | 2nd Edition
            !------------------------<COULOMB INTEGRAL>------------------------
            Fock2_assigned = 0
              carry = .true.
              assigned = 1
              Fock2_assigned(1,1) = dloop_i
              Fock2_assigned(2,1) = dloop_j
              Fock2_mic(2*dloop_i-1,2*dloop_j-1) = &
              Fock2_mic(2*dloop_i-1,2*dloop_j-1) + &
              integral*rou_m(2*dloop_k-1,2*dloop_l-1)     ! alpha->alpha Coulomb
              Fock2_mic(2*dloop_i-1,2*dloop_j-1) = &
              Fock2_mic(2*dloop_i-1,2*dloop_j-1) + &
              integral*rou_m(2*dloop_k,2*dloop_l)         ! alpha->beta Coulomb
              Fock2_mic(2*dloop_i,2*dloop_j) = &
              Fock2_mic(2*dloop_i,2*dloop_j) + &
              integral*rou_m(2*dloop_k-1,2*dloop_l-1)     ! beta->alpha Coulomb
              Fock2_mic(2*dloop_i,2*dloop_j) = &
              Fock2_mic(2*dloop_i,2*dloop_j) + &
              integral*rou_m(2*dloop_k,2*dloop_l)         ! beta->beta Coulomb
              if (dloop_k /= dloop_l) then
                Fock2_mic(2*dloop_i-1,2*dloop_j-1) = &
                Fock2_mic(2*dloop_i-1,2*dloop_j-1) + &
                integral*rou_m(2*dloop_l-1,2*dloop_k-1)   ! alpha->alpha Coulomb
                Fock2_mic(2*dloop_i-1,2*dloop_j-1) = &
                Fock2_mic(2*dloop_i-1,2*dloop_j-1) + &
                integral*rou_m(2*dloop_l,2*dloop_k)       ! alpha->beta Coulomb
                Fock2_mic(2*dloop_i,2*dloop_j) = &
                Fock2_mic(2*dloop_i,2*dloop_j) + &
                integral*rou_m(2*dloop_l-1,2*dloop_k-1)   ! beta->alpha Coulomb
                Fock2_mic(2*dloop_i,2*dloop_j) = &
                Fock2_mic(2*dloop_i,2*dloop_j) + &
                integral*rou_m(2*dloop_l,2*dloop_k)       ! beta->beta Coulomb
              end if
            do iassign = 1, assigned
              if (Fock2_assigned(1,iassign) == dloop_j .and. &
              Fock2_assigned(2,iassign) == dloop_i) then
                carry = .false.
                exit
              end if
            end do
            if (carry) then
              assigned = assigned + 1
              Fock2_assigned(1,assigned) = dloop_j
              Fock2_assigned(2,assigned) = dloop_i
              Fock2_mic(2*dloop_j-1,2*dloop_i-1) = &
              Fock2_mic(2*dloop_j-1,2*dloop_i-1) + &
              integral*rou_m(2*dloop_k-1,2*dloop_l-1)     ! alpha->alpha Coulomb
              Fock2_mic(2*dloop_j-1,2*dloop_i-1) = &
              Fock2_mic(2*dloop_j-1,2*dloop_i-1) + &
              integral*rou_m(2*dloop_k,2*dloop_l)         ! alpha->beta Coulomb
              Fock2_mic(2*dloop_j,2*dloop_i) = &
              Fock2_mic(2*dloop_j,2*dloop_i) + &
              integral*rou_m(2*dloop_k-1,2*dloop_l-1)     ! beta->alpha Coulomb
              Fock2_mic(2*dloop_j,2*dloop_i) = &
              Fock2_mic(2*dloop_j,2*dloop_i) + &
              integral*rou_m(2*dloop_k,2*dloop_l)         ! beta->beta Coulomb
              if (dloop_k /= dloop_l) then
                Fock2_mic(2*dloop_j-1,2*dloop_i-1) = &
                Fock2_mic(2*dloop_j-1,2*dloop_i-1) + &
                integral*rou_m(2*dloop_l-1,2*dloop_k-1)   ! alpha->alpha Coulomb
                Fock2_mic(2*dloop_j-1,2*dloop_i-1) = &
                Fock2_mic(2*dloop_j-1,2*dloop_i-1) + &
                integral*rou_m(2*dloop_l,2*dloop_k)       ! alpha->beta Coulomb
                Fock2_mic(2*dloop_j,2*dloop_i) = &
                Fock2_mic(2*dloop_j,2*dloop_i) + &
                integral*rou_m(2*dloop_l-1,2*dloop_k-1)   ! beta->alpha Coulomb
                Fock2_mic(2*dloop_j,2*dloop_i) = &
                Fock2_mic(2*dloop_j,2*dloop_i) + &
                integral*rou_m(2*dloop_l,2*dloop_k)       ! beta->beta Coulomb
              end if
            end if
            carry = .true.
            do iassign = 1, assigned
              if (Fock2_assigned(1,iassign) == dloop_k .and. &
              Fock2_assigned(2,iassign) == dloop_l) then
                carry = .false.
                exit
              end if
            end do
            if (carry) then
              assigned = assigned + 1
              Fock2_assigned(1,assigned) = dloop_k
              Fock2_assigned(2,assigned) = dloop_l
              Fock2_mic(2*dloop_k-1,2*dloop_l-1) = &
              Fock2_mic(2*dloop_k-1,2*dloop_l-1) + &
              integral*rou_m(2*dloop_i-1,2*dloop_j-1)     ! alpha->alpha Coulomb
              Fock2_mic(2*dloop_k-1,2*dloop_l-1) = &
              Fock2_mic(2*dloop_k-1,2*dloop_l-1) + &
              integral*rou_m(2*dloop_i,2*dloop_j)         ! alpha->beta Coulomb
              Fock2_mic(2*dloop_k,2*dloop_l) = &
              Fock2_mic(2*dloop_k,2*dloop_l) + &
              integral*rou_m(2*dloop_i-1,2*dloop_j-1)     ! beta->alpha Coulomb
              Fock2_mic(2*dloop_k,2*dloop_l) = &
              Fock2_mic(2*dloop_k,2*dloop_l) + &
              integral*rou_m(2*dloop_i,2*dloop_j)         ! beta->beta Coulomb
              if (dloop_i /= dloop_j) then
                Fock2_mic(2*dloop_k-1,2*dloop_l-1) = &
                Fock2_mic(2*dloop_k-1,2*dloop_l-1) + &
                integral*rou_m(2*dloop_j-1,2*dloop_i-1)   ! alpha->alpha Coulomb
                Fock2_mic(2*dloop_k-1,2*dloop_l-1) = &
                Fock2_mic(2*dloop_k-1,2*dloop_l-1) + &
                integral*rou_m(2*dloop_j,2*dloop_i)       ! alpha->beta Coulomb
                Fock2_mic(2*dloop_k,2*dloop_l) = &
                Fock2_mic(2*dloop_k,2*dloop_l) + &
                integral*rou_m(2*dloop_j-1,2*dloop_i-1)   ! beta->alpha Coulomb
                Fock2_mic(2*dloop_k,2*dloop_l) = &
                Fock2_mic(2*dloop_k,2*dloop_l) + &
                integral*rou_m(2*dloop_j,2*dloop_i)       ! beta->beta Coulomb
              end if
            end if
            carry = .true.
            do iassign = 1, assigned
              if (Fock2_assigned(1,iassign) == dloop_l .and. &
              Fock2_assigned(2,iassign) == dloop_k) then
                carry = .false.
                exit
              end if
            end do
            if (carry) then
              assigned = assigned + 1
              Fock2_assigned(1,assigned) = dloop_l
              Fock2_assigned(2,assigned) = dloop_k
              Fock2_mic(2*dloop_l-1,2*dloop_k-1) = &
              Fock2_mic(2*dloop_l-1,2*dloop_k-1) + &
              integral*rou_m(2*dloop_i-1,2*dloop_j-1)     ! alpha->alpha Coulomb
              Fock2_mic(2*dloop_l-1,2*dloop_k-1) = &
              Fock2_mic(2*dloop_l-1,2*dloop_k-1) + &
              integral*rou_m(2*dloop_i,2*dloop_j)         ! alpha->beta Coulomb
              Fock2_mic(2*dloop_l,2*dloop_k) = &
              Fock2_mic(2*dloop_l,2*dloop_k) + &
              integral*rou_m(2*dloop_i-1,2*dloop_j-1)     ! beta->alpha Coulomb
              Fock2_mic(2*dloop_l,2*dloop_k) = &
              Fock2_mic(2*dloop_l,2*dloop_k) + &
              integral*rou_m(2*dloop_i,2*dloop_j)         ! beta->beta Coulomb
              if (dloop_i /= dloop_j) then
                Fock2_mic(2*dloop_l-1,2*dloop_k-1) = &
                Fock2_mic(2*dloop_l-1,2*dloop_k-1) + &
                integral*rou_m(2*dloop_j-1,2*dloop_i-1)   ! alpha->alpha Coulomb
                Fock2_mic(2*dloop_l-1,2*dloop_k-1) = &
                Fock2_mic(2*dloop_l-1,2*dloop_k-1) + &
                integral*rou_m(2*dloop_j,2*dloop_i)       ! alpha->beta Coulomb
                Fock2_mic(2*dloop_l,2*dloop_k) = &
                Fock2_mic(2*dloop_l,2*dloop_k) + &
                integral*rou_m(2*dloop_j-1,2*dloop_i-1)   ! beta->alpha Coulomb
                Fock2_mic(2*dloop_l,2*dloop_k) = &
                Fock2_mic(2*dloop_l,2*dloop_k) + &
                integral*rou_m(2*dloop_j,2*dloop_i)       ! beta->beta Coulomb
              end if
            end if
            !------------------------<EXCHANGE INTEGRAL>------------------------
            Fock2_assigned = 0
              carry = .true.
              assigned = 1
              Fock2_assigned(1,1) = dloop_i
              Fock2_assigned(2,1) = dloop_k
              Fock2_mic(2*dloop_i-1,2*dloop_k-1) = &
              Fock2_mic(2*dloop_i-1,2*dloop_k-1) - &
              integral*rou_m(2*dloop_j-1,2*dloop_l-1)    ! alpha->alpha Exchange
              Fock2_mic(2*dloop_i,2*dloop_k) = &
              Fock2_mic(2*dloop_i,2*dloop_k) - &
              integral*rou_m(2*dloop_j,2*dloop_l)        ! beta->beta Exchange
              if (dloop_i == dloop_k .and. dloop_l /= dloop_j) then
                Fock2_mic(2*dloop_i-1,2*dloop_k-1) = &
                Fock2_mic(2*dloop_i-1,2*dloop_k-1) - &
                integral*rou_m(2*dloop_l-1,2*dloop_j-1)  ! alpha->alpha Exchange
                Fock2_mic(2*dloop_i,2*dloop_k) = &
                Fock2_mic(2*dloop_i,2*dloop_k) - &
                integral*rou_m(2*dloop_l,2*dloop_j)      ! beta->beta Exchange
              end if
            do iassign = 1, assigned
              if (Fock2_assigned(1,iassign) == dloop_k .and. &
              Fock2_assigned(2,iassign) == dloop_i) then
                carry = .false.
                exit
              end if
            end do
            if (carry) then
              assigned = assigned + 1
              Fock2_assigned(1,assigned) = dloop_k
              Fock2_assigned(2,assigned) = dloop_i
              Fock2_mic(2*dloop_k-1,2*dloop_i-1) = &
              Fock2_mic(2*dloop_k-1,2*dloop_i-1) - &
              integral*rou_m(2*dloop_l-1,2*dloop_j-1)    ! alpha->alpha Exchange
              Fock2_mic(2*dloop_k,2*dloop_i) = &
              Fock2_mic(2*dloop_k,2*dloop_i) - &
              integral*rou_m(2*dloop_l,2*dloop_j)        ! beta->beta Exchange
              if (dloop_k == dloop_i .and. dloop_l /= dloop_j) then
                Fock2_mic(2*dloop_k-1,2*dloop_i-1) = &
                Fock2_mic(2*dloop_k-1,2*dloop_i-1) - &
                integral*rou_m(2*dloop_j-1,2*dloop_l-1)  ! alpha->alpha Exchange
                Fock2_mic(2*dloop_k,2*dloop_i) = &
                Fock2_mic(2*dloop_k,2*dloop_i) - &
                integral*rou_m(2*dloop_j,2*dloop_l)      ! beta->beta Exchange
              end if
            end if
            carry = .true.
            do iassign = 1, assigned
              if (Fock2_assigned(1,iassign) == dloop_i .and. &
              Fock2_assigned(2,iassign) == dloop_l) then
                carry = .false.
                exit
              end if
            end do
            if (carry) then
              assigned = assigned + 1
              Fock2_assigned(1,assigned) = dloop_i
              Fock2_assigned(2,assigned) = dloop_l
              Fock2_mic(2*dloop_i-1,2*dloop_l-1) = &
              Fock2_mic(2*dloop_i-1,2*dloop_l-1) - &
              integral*rou_m(2*dloop_j-1,2*dloop_k-1)    ! alpha->alpha Exchange
              Fock2_mic(2*dloop_i,2*dloop_l) = &
              Fock2_mic(2*dloop_i,2*dloop_l) - &
              integral*rou_m(2*dloop_j,2*dloop_k)        ! beta->beta Exchange
              if (dloop_i == dloop_l .and. dloop_k /= dloop_j) then
                Fock2_mic(2*dloop_i-1,2*dloop_l-1) = &
                Fock2_mic(2*dloop_i-1,2*dloop_l-1) - &
                integral*rou_m(2*dloop_k-1,2*dloop_j-1)  ! alpha->alpha Exchange
                Fock2_mic(2*dloop_i,2*dloop_l) = &
                Fock2_mic(2*dloop_i,2*dloop_l) - &
                integral*rou_m(2*dloop_k,2*dloop_j)      ! beta->beta Exchange
              end if
            end if
            carry = .true.
            do iassign = 1, assigned
              if (Fock2_assigned(1,iassign) == dloop_l .and. &
               Fock2_assigned(2,iassign) == dloop_i) then
                carry = .false.
                exit
              end if
            end do
            if (carry) then
              assigned = assigned + 1
              Fock2_assigned(1,assigned) = dloop_l
              Fock2_assigned(2,assigned) = dloop_i
              Fock2_mic(2*dloop_l-1,2*dloop_i-1) = &
              Fock2_mic(2*dloop_l-1,2*dloop_i-1) - &
              integral*rou_m(2*dloop_k-1,2*dloop_j-1)    ! alpha->alpha Exchange
              Fock2_mic(2*dloop_l,2*dloop_i) = &
              Fock2_mic(2*dloop_l,2*dloop_i) - &
              integral*rou_m(2*dloop_k,2*dloop_j)        ! beta->beta Exchange
              if (dloop_l == dloop_i .and. dloop_k /= dloop_j) then
                Fock2_mic(2*dloop_l-1,2*dloop_i-1) = &
                Fock2_mic(2*dloop_l-1,2*dloop_i-1) - &
                integral*rou_m(2*dloop_j-1,2*dloop_k-1)  ! alpha->alpha Exchange
                Fock2_mic(2*dloop_l,2*dloop_i) = &
                Fock2_mic(2*dloop_l,2*dloop_i) - &
                integral*rou_m(2*dloop_j,2*dloop_k)      ! beta->beta Exchange
              end if
            end if
            carry = .true.
            do iassign = 1, assigned
              if (Fock2_assigned(1,iassign) == dloop_j .and. &
              Fock2_assigned(2,iassign) == dloop_k) then
                carry = .false.
                exit
              end if
            end do
            if (carry) then
              assigned = assigned + 1
              Fock2_assigned(1,assigned) = dloop_j
              Fock2_assigned(2,assigned) = dloop_k
              Fock2_mic(2*dloop_j-1,2*dloop_k-1) = &
              Fock2_mic(2*dloop_j-1,2*dloop_k-1) - &
              integral*rou_m(2*dloop_i-1,2*dloop_l-1)    ! alpha->alpha Exchange
              Fock2_mic(2*dloop_j,2*dloop_k) = &
              Fock2_mic(2*dloop_j,2*dloop_k) - &
              integral*rou_m(2*dloop_i,2*dloop_l)        ! beta->beta Exchange
              if (dloop_j == dloop_k .and. dloop_l /= dloop_i) then
                Fock2_mic(2*dloop_j-1,2*dloop_k-1) = &
                Fock2_mic(2*dloop_j-1,2*dloop_k-1) - &
                integral*rou_m(2*dloop_l-1,2*dloop_i-1)  ! alpha->alpha Exchange
                Fock2_mic(2*dloop_j,2*dloop_k) = &
                Fock2_mic(2*dloop_j,2*dloop_k) - &
                integral*rou_m(2*dloop_l,2*dloop_i)      ! beta->beta Exchange
              end if
            end if
            carry = .true.
            do iassign = 1, assigned
              if (Fock2_assigned(1,iassign) == dloop_k .and. &
              Fock2_assigned(2,iassign) == dloop_j) then
                carry = .false.
                exit
              end if
            end do
            if (carry) then
              assigned = assigned + 1
              Fock2_assigned(1,assigned) = dloop_k
              Fock2_assigned(2,assigned) = dloop_j
              Fock2_mic(2*dloop_k-1,2*dloop_j-1) = &
              Fock2_mic(2*dloop_k-1,2*dloop_j-1) - &
              integral*rou_m(2*dloop_l-1,2*dloop_i-1)    ! alpha->alpha Exchange
              Fock2_mic(2*dloop_k,2*dloop_j) = &
              Fock2_mic(2*dloop_k,2*dloop_j) - &
              integral*rou_m(2*dloop_l,2*dloop_i)        ! beta->beta Exchange
              if (dloop_k == dloop_j .and. dloop_i /= dloop_l) then
                Fock2_mic(2*dloop_k-1,2*dloop_j-1) = &
                Fock2_mic(2*dloop_k-1,2*dloop_j-1) - &
                integral*rou_m(2*dloop_i-1,2*dloop_l-1)  ! alpha->alpha Exchange
                Fock2_mic(2*dloop_k,2*dloop_j) = &
                Fock2_mic(2*dloop_k,2*dloop_j) - &
                integral*rou_m(2*dloop_i,2*dloop_l)      ! beta->beta Exchange
              end if
            end if
            carry = .true.
            do iassign = 1, assigned
              if (Fock2_assigned(1,iassign) == dloop_j .and. &
              Fock2_assigned(2,iassign) == dloop_l) then
                carry = .false.
                exit
              end if
            end do
            if (carry) then
              assigned = assigned + 1
              Fock2_assigned(1,assigned) = dloop_j
              Fock2_assigned(2,assigned) = dloop_l
              Fock2_mic(2*dloop_j-1,2*dloop_l-1) = &
              Fock2_mic(2*dloop_j-1,2*dloop_l-1) - &
              integral*rou_m(2*dloop_i-1,2*dloop_k-1)    ! alpha->alpha Exchange
              Fock2_mic(2*dloop_j,2*dloop_l) = &
              Fock2_mic(2*dloop_j,2*dloop_l) - &
              integral*rou_m(2*dloop_i,2*dloop_k)        ! beta->beta Exchange
              if (dloop_j == dloop_l .and. dloop_k /= dloop_i) then
                Fock2_mic(2*dloop_j-1,2*dloop_l-1) = &
                Fock2_mic(2*dloop_j-1,2*dloop_l-1) - &
                integral*rou_m(2*dloop_k-1,2*dloop_i-1)  ! alpha->alpha Exchange
                Fock2_mic(2*dloop_j,2*dloop_l) = &
                Fock2_mic(2*dloop_j,2*dloop_l) - &
                integral*rou_m(2*dloop_k,2*dloop_i)      ! beta->beta Exchange
              end if
            end if
            carry = .true.
            do iassign = 1, assigned
              if (Fock2_assigned(1,iassign) == dloop_l .and. &
              Fock2_assigned(2,iassign) == dloop_j) then
                carry = .false.
                exit
              end if
            end do
            if (carry) then
              assigned = assigned + 1
              Fock2_assigned(1,assigned) = dloop_l
              Fock2_assigned(2,assigned) = dloop_j
              Fock2_mic(2*dloop_l-1,2*dloop_j-1) = &
              Fock2_mic(2*dloop_l-1,2*dloop_j-1) - &
              integral*rou_m(2*dloop_k-1,2*dloop_i-1)    ! alpha->alpha Exchange
              Fock2_mic(2*dloop_l,2*dloop_j) = &
              Fock2_mic(2*dloop_l,2*dloop_j) - &
              integral*rou_m(2*dloop_k,2*dloop_i)        ! beta->beta Exchange
              if (dloop_l == dloop_j .and. dloop_i /= dloop_k) then
                Fock2_mic(2*dloop_l-1,2*dloop_j-1) = &
                Fock2_mic(2*dloop_l-1,2*dloop_j-1) - &
                integral*rou_m(2*dloop_i-1,2*dloop_k-1)  ! alpha->alpha Exchange
                Fock2_mic(2*dloop_l,2*dloop_j) = &
                Fock2_mic(2*dloop_l,2*dloop_j) - &
                integral*rou_m(2*dloop_i,2*dloop_k)      ! beta->beta Exchange
              end if
            end if
          end do
        end do
      end do
      !$omp critical
      supp1 = supp1 + Fock2_mic
      !$omp end critical
    end do
    !$omp end do
    !$omp end parallel
    if (s_h) then
      ! transform to spherical-harmonic basis
      allocate(supp2(2*cbdm,2*sbdm))
      call zgemm( 'N', 'N', 2*cbdm, 2*sbdm, 2*cbdm, c1, &
      rou_m, 2*cbdm, exc2s, 2*cbdm, c0, supp2, 2*cbdm)
      deallocate(rou_m)
      allocate(rou_m(2*sbdm,2*sbdm))
      call zgemm( 'C', 'N', 2*sbdm, 2*sbdm, 2*cbdm, c1, &
      exc2s, 2*cbdm, supp2, 2*cbdm, c0, rou_m, 2*sbdm)
      deallocate(supp2)
    
      allocate(supp2(2*cbdm,2*sbdm))
      call zgemm( 'N', 'N', 2*cbdm, 2*sbdm, 2*cbdm, c1, &
      supp1, 2*cbdm, exc2s, 2*cbdm, c0, supp2, 2*cbdm)
      deallocate(supp1)
      allocate(supp1(2*sbdm,2*sbdm))
      call zgemm( 'C', 'N', 2*sbdm, 2*sbdm, 2*cbdm, c1, &
      exc2s, 2*cbdm, supp2, 2*cbdm, c0, supp1, 2*sbdm)
      deallocate(supp2)
    end if
    
    ! orth to Fock2
    allocate(supp2(2*fbdm,2*sbdm))
    call zgemm( 'C', 'N', 2*fbdm, 2*sbdm, 2*sbdm, c1, &
    exXm, 2*sbdm, supp1, 2*sbdm, c0, supp2, 2*fbdm)
    call zgemm( 'N', 'N', 2*fbdm, 2*fbdm, 2*sbdm, c1, &
    supp2, 2*fbdm, exXm, 2*sbdm, c0, Fock2, 2*fbdm)
    deallocate(supp1)
    deallocate(supp2)
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
!> calculate <S**2> based on oper3 generated by Hartree-Fock SCF
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
      do cloop_k = 1, 2*fbdm-1, 2
        totalpha = totalpha + &
        real(conjg(oper3(cloop_k,cloop_i))*oper3(cloop_k,cloop_i),dp)
      end do
      do cloop_k = 2, 2*fbdm, 2
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
        do cloop_k = 1, 2*fbdm-1, 2
          suppa = suppa + conjg(oper3(cloop_k,cloop_i))*oper3(cloop_k+1,cloop_j)
          suppb = suppb + conjg(oper3(cloop_k+1,cloop_j))*oper3(cloop_k,cloop_i)
        end do
        albe = albe + suppa*suppb
      end do
    end do
    beal = c0
    do cloop_i = 1, electron_count
      do cloop_j = 1, electron_count
        suppa = c0
        suppb = c0
        do cloop_k = 2, 2*fbdm, 2
          suppa = suppa + conjg(oper3(cloop_k,cloop_i))*oper3(cloop_k-1,cloop_j)
          suppb = suppb + conjg(oper3(cloop_k-1,cloop_j))*oper3(cloop_k,cloop_i)
        end do
        beal = beal + suppa*suppb
      end do
    end do
    S__2 = -real(electron_count*(electron_count-4),dp)/4.0_dp + &
    real(alal/2.0_dp + bebe/2.0_dp) - real(albe/2.0_dp + beal/2.0_dp)
  end subroutine calc_S2HF
  
!------------------------------------------------------------
!> calculate <S**2> for orbitals based on oper3 generated by Hartree-Fock SCF
  subroutine calc_S2HForb(orbnum)
    implicit none
    integer,intent(in) :: orbnum
    integer :: zloop_i     !loop variables for calc_S2HForb
    real(dp) :: suppa, suppb
    totalphaorb = 0.0_dp
    do zloop_i = 1, 2*fbdm-1, 2
      totalphaorb = totalphaorb + &
      real(conjg(oper3(zloop_i,orbnum))*oper3(zloop_i,orbnum),dp)
    end do
    totbetaorb = 0.0_dp
    do zloop_i = 2, 2*fbdm, 2
      totbetaorb = totbetaorb + &
      real(conjg(oper3(zloop_i,orbnum))*oper3(zloop_i,orbnum),dp)
    end do
    suppa = c0
    suppb = c0
    do zloop_i = 2, 2*fbdm, 2
      suppa = suppa + conjg(oper3(zloop_i,orbnum))*oper3(zloop_i-1,orbnum)
      suppb = suppb + conjg(oper3(zloop_i-1,orbnum))*oper3(zloop_i,orbnum)
    end do
    S__2orb = 3.0/4.0 + totalphaorb*(totalphaorb-1.0)/2.0 + &
    totbetaorb*(totbetaorb-1.0)/2.0 - real(suppa*suppb)
    Szorb = (totalphaorb-totbetaorb)/2.0_dp
  end subroutine calc_S2HForb
  
!-----------------------------------------------------------------------
!> dump molecular orbital information to .molden file
  subroutine dump_molden()
    implicit none
    integer :: channel, dmi, dmj, dmk
    if (.not. allocated(AO2MO)) &
    call terminate('dump molecular orbital failed, AO2MO is empty')
    
    ! molden file contains the real part of molecular orbital
    open(newunit=channel, file=trim(address_job)//'-real.molden', &
    status='replace', action='write', iostat=ios)
    if (ios /= 0) call terminate('dump .molden failed')
    write(channel, '(A)') '[Molden Format]'
    write(channel, '(A)') '[Title]'
    write(channel, '(A)') &
    'The real part of molecular orbitals of job '//trim(address_job)
    write(channel, *)
    ! molecular geometry
    write(channel, '(A)') '[Atoms] AU'
    do dmi = 1, atom_count
      write(channel, '(A,I3,I3,F13.7,F13.7,F13.7)') &
      element_list(molecular(dmi)%atom_number), dmi, &
      molecular(dmi)%atom_number, molecular(dmi)%nucleus_position(1), &
      molecular(dmi)%nucleus_position(2), molecular(dmi)%nucleus_position(3)
    end do
    ! basis of each atom
    write(channel, '(A)') '[GTO]'
    do dmi = 1, atom_count
      write(channel, '(I3,I2)') dmi, 0
      do dmj = 0, shell_in_element(molecular(dmi) % atom_number)-1
        if (atom_basis(molecular(dmi)%basis_number+dmj)%&
        angular_quantum_number == 0) then
          write(channel, '(A2,I2,A)') 's', &
          atom_basis(molecular(dmi)%basis_number+dmj)%contraction,' 1.0'
        else if (atom_basis(molecular(dmi)%basis_number+dmj)%&
        angular_quantum_number == 1) then
          write(channel, '(A2,I2,A)') 'p', &
          atom_basis(molecular(dmi)%basis_number+dmj)%contraction,' 1.0'
        else if (atom_basis(molecular(dmi)%basis_number+dmj)%&
        angular_quantum_number == 2) then
          write(channel, '(A2,I2,A)') 'd', &
          atom_basis(molecular(dmi)%basis_number+dmj)%contraction,' 1.0'
        else if (atom_basis(molecular(dmi)%basis_number+dmj)%&
        angular_quantum_number == 3) then
          write(channel, '(A2,I2,A)') 'f', &
          atom_basis(molecular(dmi)%basis_number+dmj)%contraction,' 1.0'
        else if (atom_basis(molecular(dmi)%basis_number+dmj)%&
        angular_quantum_number == 4) then
          write(channel, '(A2,I2,A)') 'g', &
          atom_basis(molecular(dmi)%basis_number+dmj)%contraction,' 1.0'
        end if
        do dmk = 1, atom_basis(molecular(dmi)%basis_number+dmj)%contraction
          write(channel, '(F13.7,F13.7)') &
          atom_basis(molecular(dmi)%basis_number+dmj)%exponents(dmk), &
          atom_basis(molecular(dmi)%basis_number+dmj)%coefficient(dmk)
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
        write(channel, '(I4,F20.12)') dmj, real(AO2MO(2*dmj-1, dmi))
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
        write(channel, '(I4,F20.12)') dmj, real(AO2MO(2*dmj, dmi))
      end do
    end do
    close(channel)
    
    if (DKH_order /= 0) then
      ! molden file contains the maginary part of molecular orbital
      open(newunit=channel, file=trim(address_job)//'-img.molden', &
      status='replace', action='write', iostat=ios)
      if (ios /= 0) call terminate('dump .molden failed')
      write(channel, '(A)') '[Molden Format]'
      write(channel, '(A)') '[Title]'
      write(channel, '(A)') &
      'The imaginary part of molecular orbitals of job '//trim(address_job)
      write(channel, *)
      ! molecular geometry
      write(channel, '(A)') '[Atoms] AU'
      do dmi = 1, atom_count
        write(channel, '(A,I3,I3,F13.7,F13.7,F13.7)') &
        element_list(molecular(dmi)%atom_number), dmi, &
        molecular(dmi)%atom_number, molecular(dmi)%nucleus_position(1), &
        molecular(dmi)%nucleus_position(2), molecular(dmi)%nucleus_position(3)
      end do
      ! basis of each atom
      write(channel, '(A)') '[GTO]'
      do dmi = 1, atom_count
        write(channel, '(I3,I2)') dmi, 0
        do dmj = 0, shell_in_element(molecular(dmi) % atom_number)-1
          if (atom_basis(molecular(dmi)%basis_number+dmj)%&
          angular_quantum_number == 0) then
            write(channel, '(A2,I2,A)') 's', &
            atom_basis(molecular(dmi)%basis_number+dmj)%contraction,' 1.0'
          else if (atom_basis(molecular(dmi)%basis_number+dmj)%&
          angular_quantum_number == 1) then
            write(channel, '(A2,I2,A)') 'p', &
            atom_basis(molecular(dmi)%basis_number+dmj)%contraction,' 1.0'
          else if (atom_basis(molecular(dmi)%basis_number+dmj)%&
          angular_quantum_number == 2) then
            write(channel, '(A2,I2,A)') 'd', &
            atom_basis(molecular(dmi)%basis_number+dmj)%contraction,' 1.0'
          else if (atom_basis(molecular(dmi)%basis_number+dmj)%&
          angular_quantum_number == 3) then
            write(channel, '(A2,I2,A)') 'f', &
            atom_basis(molecular(dmi)%basis_number+dmj)%contraction,' 1.0'
          else if (atom_basis(molecular(dmi)%basis_number+dmj)%&
          angular_quantum_number == 4) then
            write(channel, '(A2,I2,A)') 'g', &
            atom_basis(molecular(dmi)%basis_number+dmj)%contraction,' 1.0'
          end if
          do dmk = 1, atom_basis(molecular(dmi)%basis_number+dmj)%contraction
            write(channel, '(F13.7,F13.7)') &
            atom_basis(molecular(dmi)%basis_number+dmj)%exponents(dmk), &
            atom_basis(molecular(dmi)%basis_number+dmj)%coefficient(dmk)
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
          write(channel, '(I4,F20.12)') dmj, aimag(AO2MO(2*dmj-1, dmi))
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
          write(channel, '(I4,F20.12)') dmj, aimag(AO2MO(2*dmj, dmi))
        end do
      end do
      close(channel)
    end if
  end subroutine dump_molden
  
end module SCF