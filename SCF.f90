! |-------<THOMAS RELATIVISTIC ELECTRONIC STRUCTURE CALCULATIONS>------|
! |------------------------------<DIRAC4Π>-----------------------------|
! |--------------------------------<SCF>-------------------------------|
! |----------------<SELF CONSISTENT FIELD CALCULATIONS>----------------|

module SCF
    use Fundamentals
    use Hamiltonian
    use Atoms
    use lapack95
    use omp_lib
    
    ! density matrices and molecular orbital coefficients
    integer Nalpha, Nbeta                                            ! number of alpha and beta elctron for initial wave function
    complex(dp),allocatable :: AO2MO(:,:)                            ! linear combination coefficients of MOs, basis spinor represented by (χ1,0), (0,χ1), (χ2,0), (0,χ2) ... (χn,0), (0,χn)
    real(dp),allocatable :: Gaualpha(:)                              ! Gaussian alpha orbital coefficient
    real(dp),allocatable :: AO2MOalpha(:,:)                          ! Gaussian alpha orbital coefficient
    real(dp),allocatable :: Gaubeta(:)                               ! Gaussian beta orbital coefficient
    real(dp),allocatable :: AO2MObeta(:,:)                           ! Gaussian beta orbital coefficient
    complex(dp),allocatable :: rou_m(:,:)                            ! density matrix, complex Hermitian matrix
    real(dp),allocatable :: AOsupp(:,:)
    integer :: evl_count_f                                           ! number of eigenvalues found in Fock diagonalization
    complex(dp),allocatable :: Fock(:,:)                             ! Fock matrix
    integer,allocatable :: isupp_ev_f(:)
    integer :: fliwork
    integer :: flwork
    integer :: lrwork
    complex(dp),allocatable :: fwork(:)                              ! work of Fock matrix for input of zheevr
    real(dp),allocatable :: rwork(:)                                 ! rwork of Fock matrix for input of zheevr
    integer,allocatable :: fiwork(:)                                 ! iwork of Fock matrix for input of zheevr
    integer :: finfo                                                 ! information of calling lapack functions
    
    logical :: ini_rou = .true.                                      ! whether to read the initial density matrix
    logical :: inin = .true.
    
    
    !--------------<one electron Fock>--------------
    complex(dp),allocatable :: Fock1(:,:)                            ! one electron Fock matrix
    real(dp),allocatable :: oper1(:,:)                               ! operator matrix used to calculate one electron Fock matrix
    real(dp),allocatable :: oper2(:,:)                               ! operator matrix used to calculate one electron Fock matrix
    complex(dp),allocatable :: oper3(:,:)                            ! operator matrix used to calculate one electron Fock matrix
    complex(dp),allocatable :: oper4(:,:)                            ! operator matrix used to calculate one electron Fock matrix
    complex(dp),allocatable :: oper5(:,:)                            ! operator matrix used to calculate one electron Fock matrix
    complex(dp),allocatable :: oper6(:,:)                               ! operator matrix used to calculate one electron Fock matrix
    real(dp),allocatable :: Ve(:,:)                                  ! V/(E+E')
    real(dp),allocatable :: pxVepx(:,:)                              ! pxVpx/(E+E')
    real(dp),allocatable :: pyVepy(:,:)                              ! pyVpy/(E+E')
    real(dp),allocatable :: pzVepz(:,:)                              ! pzVpz/(E+E')
    real(dp),allocatable :: pxVepy(:,:)                              ! pxVpy/(E+E')
    real(dp),allocatable :: pyVepx(:,:)                              ! pyVpx/(E+E')
    real(dp),allocatable :: pxVepz(:,:)                              ! pxVpz/(E+E')
    real(dp),allocatable :: pzVepx(:,:)                              ! pzVpx/(E+E')
    real(dp),allocatable :: pyVepz(:,:)                              ! pyVpz/(E+E')
    real(dp),allocatable :: pzVepy(:,:)                              ! pzVpy/(E+E')
    real(dp),allocatable :: Ap(:,:)                                  ! Ap
    real(dp),allocatable :: ApRp(:,:)                                ! ApRp
    complex(dp),allocatable :: SRp(:,:)                              ! ApRp(1+p^2/4c^2)
    complex(dp),allocatable :: ARVRA(:,:)                            ! ApRpVRpAp
    complex(dp),allocatable :: ARVeRA(:,:)                           ! ApRp(V/(E+E'))RpAp
    complex(dp),allocatable :: AVA(:,:)                              ! ApVAp
    complex(dp),allocatable :: AVeA(:,:)                             ! Ap(V/(E+E'))Ap
    complex(dp),allocatable :: exAO2p2(:,:)                          ! extended AO2p2 matrix
    
    
    !--------------<two electron Fock>--------------
    complex(dp),allocatable :: Fock2(:,:)                            ! two electron Fock matrix
    complex(dp),allocatable :: Fock2_mic(:,:)                        ! two electron Fock matrix
    type basis_inf_type                                              ! for openMP parallel computation
        integer :: atom                                              ! atom number
        integer :: shell                                             ! shell number
        integer :: L                                                 ! angular quantum number
        integer :: M                                                 ! magnetic quantum number
    end type basis_inf_type
    
    type(basis_inf_type),allocatable :: basis_inf(:)
    
    real(dp),allocatable :: swintegral(:,:)                          ! <ij||ij> as well as <kl||kl>
    logical :: ndschwarz = .true.
    logical,allocatable :: Fock2_assigned(:,:)                       ! avoid duplicate assignment of Fock2
    complex(dp),allocatable :: oper(:,:)
    
    ! orbital energy and molecular energy
    real(dp),allocatable :: orbE(:)                                  ! orbital energy
    real(dp),allocatable :: oper3_pre(:,:)                           ! privious oper3 of current iteration
    complex(dp),allocatable :: oper3_history(:,:,:)                  ! coefficients of subsp iteration
    real(dp) :: molE_pre, molE                                    	 ! molecular energy
    real(dp) :: nucE   												 ! nuclear repulsion energy
    real(dp) :: eleE                                                 ! electron repulsion energy
    real(dp) :: T                                                    ! kinetic energy
    real(dp) :: V                                                    ! electron-nuclear attraction energy
    real(dp) :: ESOC                                                 ! SOC energy
    real(dp) :: ESR                                                  ! second relativistic Thomas precession energy

    ! DIIS Ax=B
	real(dp),allocatable :: Rsd(:,:,:)                               ! residuals of oper3
	real(dp),allocatable :: DIISmat(:,:)                             ! A
    real(dp),allocatable :: DIISwork(:)
    integer :: lDIISwork
	integer,allocatable :: ipiv(:)
	integer :: DIIsinfo
	
    
    contains

!------------------------------------------------------------
! standard Hartree Fock SCF procedure based on the Pulay (DIIS) mixing method for DKH0 & DKH2 Hamiltonian 
    subroutine DKH_Hartree_Fock()
        write(60,'(a)') 'Module SCF:'
        write(60,'(a)') '   construct one electron Fock matrix'
        call Fock1e()
        write(60,'(a)') '   complete! stored in Fock1'
        write(60,'(a)') '   calculate nuclear repulsion energy'
        nucE = 0.0_dp
        do loop_i = 1, atom_count
        	do loop_j = loop_i + 1, atom_count
        		nucE = nucE + (real(molecular(loop_i) % atom_number) * real(molecular(loop_j) % atom_number))/sqrt(&
				(molecular(loop_i) % nucleus_position(1) - molecular(loop_j) % nucleus_position(1))**2 + &
				(molecular(loop_i) % nucleus_position(2) - molecular(loop_j) % nucleus_position(2))**2 + &
				(molecular(loop_i) % nucleus_position(3) - molecular(loop_j) % nucleus_position(3))**2)
			end do
        end do
        if (isnan(nucE)) call terminate('nuclear repulsive energy anomaly, possibly due to overlapping atomic coordinates')
        write(60,'(a45,e15.7)') '   complete! nuclear repulsive energy: (A.U.)', nucE
        cpu_threads = omp_get_num_procs()
        write(60,'(a,i3,a,i3)') '   threads using:',threads_use,' CPU threads:',cpu_threads
        if (cpu_threads <= threads_use) then
            write(*,*) 'Warning! Calculation will be performed serially, CPU threads is',cpu_threads
            write(60,'(a)') '   calculation will be performed SERIALLY!'
        end if
        write(60,'(a,i4,a,e15.7)') '   SCF maxloop:',maxloop,';  convergence tolerance:',conver_tol
        write(60,'(a,i3,a,i3)') '   nodiis =',nodiis,';  subsp =',subsp
        write(60,'(a)') '   -----<SCF>-----'
        allocate(Fock(2*basis_dimension,2*basis_dimension))
        allocate(orbE(2*basis_dimension))
        allocate(oper6(2*basis_dimension,2*basis_dimension))
        allocate(isupp_ev_f(4*basis_dimension))
        allocate(oper3(2*basis_dimension,2*basis_dimension))
        allocate(oper4(2*basis_dimension,2*basis_dimension))
        

        ! DIIS
        ! AO2MO(new) = Σ(i,subsp) DIIScoe(i)*(oper3history(i)+nudge*Rsd(i))
        allocate(Rsd(subsp,2*basis_dimension,2*basis_dimension))
        allocate(DIISmat(subsp+1,subsp+1))
        allocate(ipiv(subsp+1))
        allocate(oper3_history(subsp, 2*basis_dimension, 2*basis_dimension))
        allocate(oper3_pre(2*basis_dimension, 2*basis_dimension))
        do loop_i = 1, maxloop
            write(60,*)
            write(60,'(a,i3)') '   SCF loop ',loop_i
            if (loop_i == 1) then
                write(60,'(a)') '   read density matrix from Gaussian checkpoint file'
            else
                write(60,'(a)') '   construct density matrix'
            end if
            call assign_rou()
            write(60,'(a)') '   complete! stored in rou_m'
            write(60,'(a)') '   construct two electron Fock matrix'
            call Fock2e()
            write(60,'(a)') '   complete! stored in Fock2'
            write(60,'(a)') '   diagonalization of Fock matrix'
            Fock = cmplx(0.0_dp,0.0_dp,dp)
            do loop_j = 1, 2*basis_dimension
                do loop_k = loop_j, 2*basis_dimension
                    Fock(loop_j,loop_k) = Fock1(loop_j,loop_k) + Fock2(loop_j,loop_k)
                end do
            end do
            allocate(fwork(1))
            allocate(fiwork(1))
            allocate(rwork(1))
            call zheevr('V','A','U',2*basis_dimension,Fock,2*basis_dimension,0.0_dp,0.0_dp,0,0,dlamch('S'),&
                evl_count_f,orbE,oper3,2*basis_dimension,isupp_ev_f,fwork,-1,rwork,-1,fiwork,-1,finfo)
            flwork = int(fwork(1))
            fliwork = fiwork(1)
            lrwork = int(rwork(1))
            deallocate(fwork)
            deallocate(fiwork)
            deallocate(rwork)
            allocate(fwork(flwork))
            allocate(fiwork(fliwork))
            allocate(rwork(lrwork))
            call zheevr('V','A','U',2*basis_dimension,Fock,2*basis_dimension,0.0_dp,0.0_dp,0,0,dlamch('S'),&
                evl_count_f,orbE,oper3,2*basis_dimension,isupp_ev_f,fwork,flwork,rwork,lrwork,fiwork,fliwork,finfo)
            if (finfo < 0) then
                call terminate('Fock matrix diagonalization failure, illegal input of zheevr')
            else if (finfo > 0) then
                call terminate('Fock matrix diagonalization failure, internal error of zheevr')
            end if
            if (evl_count_f < 2*basis_dimension) then
				call terminate('number of MO less than 2*basis_dimension')
			else
				write(60,'(a,i5,a)') '   complete!',evl_count_f,' eigenvectors found'
            end if
            
            ! print orbital energy
            write(60,'(a)') '   print orbital energy to LUMO+2: (A.U.)'
			do loop_j = 1, electron_count
                if (loop_j < electron_count) then
				    write(60,'(a,i3.3,e20.6)') '      HOMO-',electron_count - loop_j,orbE(loop_j)
                else
                    write(60,'(a,e20.6)') '      HOMO    ',orbE(loop_j)
                end if
            end do
            do loop_j = electron_count+1, electron_count+3
                if (loop_j == electron_count + 1) then
                    write(60,'(a,e20.6)') '      LUMO    ',orbE(loop_j)
                else
                    write(60,'(a,i3.3,e20.6)') '      LUMO+',loop_j - electron_count - 1,orbE(loop_j)
                end if
            end do
            deallocate(fwork)
            deallocate(fiwork)
            deallocate(rwork)
            ! energy components calculation
            write(60,'(a)') '   calculate energy components'
            eleE = 0.0_dp
            call zgemm( 'C', 'N', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), oper3, 2*basis_dimension, Fock2, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), oper6, 2*basis_dimension)
            call zgemm( 'N', 'N', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), oper6, 2*basis_dimension, oper3, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), oper4, 2*basis_dimension)
            do loop_j = 1, electron_count
				eleE = eleE + real(oper4(loop_j,loop_j))
            end do
            write(60,'(a,e15.7)') '      electron repulsive energy: (A.U.)', eleE
                ! calculate molE
                if (loop_i /= 1) molE_pre = molE
                molE = nucE - eleE/2.0_dp
			    do loop_j = 1, electron_count
				    molE = molE + orbE(loop_j)
                end do
            T = 0.0_dp
            call zgemm( 'C', 'N', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), oper3, 2*basis_dimension, exi_T_j, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), oper6, 2*basis_dimension)
            call zgemm( 'N', 'N', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), oper6, 2*basis_dimension, oper3, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), oper4, 2*basis_dimension)
            do loop_j = 1, electron_count
				T = T + real(oper4(loop_j,loop_j))
            end do
            write(60,'(a,e15.7)') '      electron kinetic energy: (A.U.)', T
            V = 0.0_dp
            call zgemm( 'C', 'N', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), oper3, 2*basis_dimension, exi_V_j, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), oper6, 2*basis_dimension)
            call zgemm( 'N', 'N', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), oper6, 2*basis_dimension, oper3, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), oper4, 2*basis_dimension)
            do loop_j = 1, electron_count
				V = V + real(oper4(loop_j,loop_j))
            end do
            write(60,'(a,e15.7)') '      electron-nuclear attraction energy: (A.U.)', V
            if (DKH_order == 2) then
                ESOC = 0.0_dp
                call zgemm( 'C', 'N', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), oper3, 2*basis_dimension, exSOC, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), oper6, 2*basis_dimension)
                call zgemm( 'N', 'N', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), oper6, 2*basis_dimension, oper3, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), oper4, 2*basis_dimension)
                do loop_j = 1, electron_count
				    ESOC = ESOC + real(oper4(loop_j,loop_j))
                end do
                write(60,'(a,e15.7)') '      spin-orbital coupling energy: (A.U.)', ESOC
                if (SRTP_type) then
                    ESR = 0.0_dp
                    call zgemm( 'C', 'N', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), oper3, 2*basis_dimension, exSR, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), oper6, 2*basis_dimension)
                    call zgemm( 'N', 'N', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), oper6, 2*basis_dimension, oper3, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), oper4, 2*basis_dimension)
                    do loop_j = 1, electron_count
				        ESR = ESR + real(oper4(loop_j,loop_j))
                    end do
                    write(60,'(a,e15.7)') '      second relativistic Thomas precession energy: (A.U.)', ESR
                end if
            end if
            write(60,'(a,e15.7)') '      -V/T ', -(molE-T)/T
            ! convergence check
            if (loop_i == 1) then
                write(60,'(a,f12.6)') '   molecular energy (A.U.) = ',molE
                ! de-Lowdin orthogonalization
                AO2MO = oper3
				!call zgemm( 'N', 'N', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), exS0_5, 2*basis_dimension, oper3, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), AO2MO, 2*basis_dimension)
            else 
                if (abs(molE - molE_pre) < conver_tol) exit
                write(60,'(a,f12.6,a,e15.7)') '   molecular energy (A.U.) = ',molE,'; delta E = ',molE - molE_pre
                write(60,'(a)') '   convergence tolerance not met'
            end if

            write(60,'(a)') '   DIIS information'
            ! generate next AO2MO by DIIS method
			if (loop_i <= nodiis - subsp) then
	            if (loop_i == nodiis - subsp) oper3_pre = oper3
                write(60,'(a)') '      no DIIS acceleration'
			else if (nodiis - subsp < loop_i .and. loop_i <= nodiis) then
				! 更新Rsd
				do loop_j = 1, 2*basis_dimension
                    do loop_k = 1, 2*basis_dimension
					    Rsd(loop_i-(nodiis-subsp), loop_j, loop_k) = oper3(loop_j, loop_k) - oper3_pre(loop_j, loop_k) ! Rsd存储迭代产生的轨道能量残差
                    end do
				end do
				oper3_pre = oper3
				! 更新oper3_history
				do loop_j = 1, 2*basis_dimension
					do loop_k = 1, 2*basis_dimension
						oper3_history(loop_i-(nodiis-subsp), loop_j, loop_k) = oper3(loop_j, loop_k)
					end do
                end do
                write(60,'(a)') '      DIIS subspace filling'
			else
				! 更新Rsd
				do loop_j = 2, subsp
					do loop_k = 1, 2*basis_dimension
                        do loop_l = 1, 2*basis_dimension
						    Rsd(loop_j - 1, loop_k, loop_l) = Rsd(loop_j, loop_k, loop_l)
                        end do
					end do
				end do
				do loop_j = 1, 2*basis_dimension
                    do loop_k = 1, 2*basis_dimension
					    Rsd(subsp, loop_j, loop_k) = oper3(loop_j, loop_k) - oper3_pre(loop_j, loop_k) ! Rsd存储迭代产生的轨道能量残差
                    end do
				end do
				oper3_pre = oper3
				! 更新oper3_history
				do loop_j = 2, subsp
					do loop_k = 1, 2*basis_dimension
						do loop_l = 1, 2*basis_dimension
							oper3_history(loop_j - 1, loop_k, loop_l) = oper3_history(loop_j, loop_k, loop_l)
						end do
					end do
				end do
				do loop_j = 1, 2*basis_dimension
					do loop_k = 1, 2*basis_dimension
						oper3_history(subsp, loop_j, loop_k) = oper3(loop_j, loop_k)
					end do
				end do
				! 根据当下Rsd产生DIISmat
				DIISmat = 0.0_dp
				do loop_j = 1, subsp
					do loop_k = 1, subsp
						do loop_l = 1, 2*basis_dimension
                            do loop_m = 1, 2*basis_dimension
							    DIISmat(loop_j, loop_k) = DIISmat(loop_j, loop_k) + Rsd(loop_j, loop_l, loop_m)*Rsd(loop_k, loop_l, loop_m)
                            end do
						end do
					end do
				end do
				do loop_j = 1, subsp
					DIISmat(subsp+1, loop_j) = 1.0_dp
					DIISmat(loop_j, subsp+1) = 1.0_dp
				end do
				! 求解残差方程
				! dgesv and dspsv Will cause V_Integral_2e conflict for unknown reason
				! since DIISmat (and its inverse) is real symmetric, plus the column vector is simple, use the inverse of DIISmat to solve directly
                call dgetrf( subsp+1, subsp+1, DIISmat, subsp+1, ipiv, DIISinfo )
                if (DIISinfo < 0) then
	                call terminate('DIIS solution failure, illegal input of dgetrf')
	            else if (DIISinfo > 0) then
	                call terminate('DIIS solution failure, internal error of dgetrf')
                end if
                allocate(DIISwork(1))
                call dgetri( subsp+1, DIISmat, subsp+1, ipiv, DIISwork, -1, DIISinfo )
                lDIISwork = int(DIISwork(1))
                deallocate(DIISwork)
                allocate(DIISwork(lDIISwork))
                call dgetri( subsp+1, DIISmat, subsp+1, ipiv, DIISwork, lDIISwork, DIISinfo )
                deallocate(DIISwork)
				if (DIISinfo < 0) then
	                call terminate('DIIS solution failure, illegal input of dgetri')
	            else if (DIISinfo > 0) then
	                call terminate('DIIS solution failure, internal error of dgetri')
	            end if
	            ! 根据DIIScoe和oper3_history产生新的oper3
	            oper3 = cmplx(0.0_dp,0.0_dp,dp)
	            do loop_j = 1, subsp
	            	do loop_k = 1, 2*basis_dimension
		            	do loop_l = 1, 2*basis_dimension
							oper3(loop_k,loop_l) = oper3(loop_k,loop_l) + DIISmat(loop_j,subsp+1) * (oper3_history(loop_j,loop_k,loop_l) + nudge*Rsd(loop_j,loop_k,loop_l))
						end do
					end do
                end do
                write(60,'(a26,f10.6)') '      predicted residual =', -DIISmat(subsp+1,subsp+1)
                do loop_j = 1, subsp
                    write(60,'(a16,i2,a2,f10.6)') '      subsp coe ',loop_j,' =', DIISmat(loop_j,subsp+1)
                end do
            end if
            ! de-Lowdin orthogonalization
            AO2MO = oper3
            !call zgemm( 'N', 'N', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), exS0_5, 2*basis_dimension, oper3, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), AO2MO, 2*basis_dimension)
	  		write(60,'(a)') '   MO coefficient stored in AO2MO, orbital energy stored in orbE'
        end do
        if (abs(molE - molE_pre) < conver_tol) then
            write(60,'(a)') 'SCF done! exit module SCF'
            if (DKH_order == 0) then
            	write(60,'(a,f12.6)') 'DKH0 Hatree-Fock SCF energy (A.U.): ',molE
            else if (DKH_order == 2) then
				write(60,'(a,f12.6)') 'DKH2 Hatree-Fock SCF energy (A.U.): ',molE
			end if
            write(60,'(a)') 'MO coefficients and orbital energy were printed in AO2MO.ditm and orbE.ditm'
        else
            write(60,'(a)') 'SCF convergence failed, exit module SCF'
        end if
    end subroutine DKH_Hartree_Fock
    
!------------------------------------------------------------
! gives an initial density matrix based on charge and spin multiplicity
    subroutine assign_rou()
        integer :: ploop_i,ploop_j,ploop_k,ploop_l                      ! ploop is loop variables only for subroutine assign_rou
        integer :: mat_dimension
        real(dp) :: rdMO(5)
        character(len = 512) :: line_str
        if (ini_rou) then
            ini_rou = .false.
            allocate(AO2MO(2*basis_dimension,2*basis_dimension))
            allocate(Gaualpha(basis_dimension*basis_dimension))
            allocate(Gaubeta(basis_dimension*basis_dimension))
            allocate(AO2MOalpha(basis_dimension,basis_dimension))
            allocate(AO2MObeta(basis_dimension,basis_dimension))
            allocate(rou_m(2*basis_dimension,2*basis_dimension))
            Nalpha = (electron_count - (spin_mult - 1)) / 2 + (spin_mult - 1)
            Nbeta = (electron_count - (spin_mult - 1)) / 2
            open(61, file = 'E:\TRESC\H2\Gaualpha.itm', status = 'old', action = 'read', iostat = ios)
            if (ios /= 0) then
                call system('rwfdump E:\TRESC\H2\Gaujob.chk E:\TRESC\H2\Gaualpha.itm 524R')
                open(61, file = 'E:\TRESC\H2\Gaualpha.itm', status = 'old', action = 'read', iostat = ios)
                if (ios /= 0) call terminate('Cannot generate Gaussian alpha orbital for initial density matrix')
            end if
            open(62, file = 'E:\TRESC\H2\Gaubeta.itm', status = 'old', action = 'read', iostat = ios)
            if (ios /= 0) then
                call system('rwfdump E:\TRESC\H2\Gaujob.chk E:\TRESC\H2\Gaubeta.itm 526R')
                open(62, file = 'E:\TRESC\H2\Gaubeta.itm', status = 'old', action = 'read', iostat = ios)
                if (ios /= 0) call terminate('Cannot generate Gaussian beta orbital for initial density matrix')
            end if
            do while(.true.)
                read(61,'(a512)') line_str
                if (index(line_str,'Dump of file') /= 0) exit
            end do
            if (index(line_str(index(line_Str,'length')-4:index(line_Str,'length')-2),'524') == 0) call terminate('Gaussian alpha orbital coefficient should be RWF 524')
            read(line_str(index(line_Str,'(')-9:index(line_Str,'(')-2),'(I)',iostat = ios) mat_dimension
            if (ios /= 0) call terminate('Unmatched Gaualpha content')
            if (mat_dimension /= basis_dimension*basis_dimension) call terminate('Gaussian alpha orbital dimension not match')
            do ploop_i=1,basis_dimension*basis_dimension/5
                read(61,*) rdMO
                Gaualpha((ploop_i-1)*5 + 1) = rdMO(1)
                Gaualpha((ploop_i-1)*5 + 2) = rdMO(2)
                Gaualpha((ploop_i-1)*5 + 3) = rdMO(3)
                Gaualpha((ploop_i-1)*5 + 4) = rdMO(4)
                Gaualpha((ploop_i-1)*5 + 5) = rdMO(5)
            end do
            if (mod(basis_dimension*basis_dimension,5) /= 0) then
                if (mod(basis_dimension*basis_dimension,5) == 1) read(61,*) Gaualpha(ploop_i*5 + 1)
                if (mod(basis_dimension*basis_dimension,5) == 2) read(61,*) Gaualpha(ploop_i*5 + 1),Gaualpha(ploop_i*5 + 2)
                if (mod(basis_dimension*basis_dimension,5) == 3) read(61,*) Gaualpha(ploop_i*5 + 1),Gaualpha(ploop_i*5 + 2),Gaualpha(ploop_i*5 + 3)
                if (mod(basis_dimension*basis_dimension,5) == 4) read(61,*) Gaualpha(ploop_i*5 + 1),Gaualpha(ploop_i*5 + 2),Gaualpha(ploop_i*5 + 3),Gaualpha(ploop_i*5 + 4)
            end if
            AO2MOalpha = transpose(reshape(Gaualpha,[basis_dimension,basis_dimension]))
            !call dgemm( 'T', 'N', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, Sp0_5, basis_dimension, AOsupp, basis_dimension, 0.0_dp, AO2MOalpha, basis_dimension)
            close(61)
            do while(.true.)
                read(62,'(a512)') line_str
                if (index(line_str,'Dump of file') /= 0) exit
            end do
            if (index(line_str(index(line_Str,'length')-4:index(line_Str,'length')-2),'526') == 0) call terminate('Gaussian beta orbital coefficient should be RWF 526')
            read(line_str(index(line_Str,'(')-9:index(line_Str,'(')-2),'(I)',iostat = ios) mat_dimension
            if (ios /= 0) call terminate('Unmatched Gaubeta content')
            if (mat_dimension /= basis_dimension*basis_dimension) call terminate('Gaussian beta orbital dimension not match')
            do ploop_i=1,basis_dimension*basis_dimension/5
                read(62,*) rdMO
                Gaubeta((ploop_i-1)*5 + 1) = rdMO(1)
                Gaubeta((ploop_i-1)*5 + 2) = rdMO(2)
                Gaubeta((ploop_i-1)*5 + 3) = rdMO(3)
                Gaubeta((ploop_i-1)*5 + 4) = rdMO(4)
                Gaubeta((ploop_i-1)*5 + 5) = rdMO(5)
            end do
            if (mod(basis_dimension*basis_dimension,5) /= 0) then
                if (mod(basis_dimension*basis_dimension,5) == 1) read(62,*) Gaubeta(ploop_i*5 + 1)
                if (mod(basis_dimension*basis_dimension,5) == 2) read(62,*) Gaubeta(ploop_i*5 + 1),Gaubeta(ploop_i*5 + 2)
                if (mod(basis_dimension*basis_dimension,5) == 3) read(62,*) Gaubeta(ploop_i*5 + 1),Gaubeta(ploop_i*5 + 2),Gaubeta(ploop_i*5 + 3)
                if (mod(basis_dimension*basis_dimension,5) == 4) read(62,*) Gaubeta(ploop_i*5 + 1),Gaubeta(ploop_i*5 + 2),Gaubeta(ploop_i*5 + 3),Gaubeta(ploop_i*5 + 4)
            end if
            AO2MObeta = transpose(reshape(Gaubeta,[basis_dimension,basis_dimension]))
            !call dgemm( 'T', 'N', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, Sp0_5, basis_dimension, AOsupp, basis_dimension, 0.0_dp, AO2MObeta, basis_dimension)
            close(62)
            deallocate(Gaualpha)
            deallocate(Gaubeta)
            rou_m = cmplx(0.0_dp,0.0_dp,dp)
            do ploop_i = 1,basis_dimension
                do ploop_j = 1,basis_dimension
                    do ploop_k = 1,Nalpha
                        rou_m(2*ploop_i-1,2*ploop_j-1) = rou_m(2*ploop_i-1,2*ploop_j-1) + AO2MOalpha(ploop_k,ploop_i)*AO2MOalpha(ploop_k,ploop_j)
                    end do
                    do ploop_k = 1,Nbeta
                        rou_m(2*ploop_i,2*ploop_j) = rou_m(2*ploop_i,2*ploop_j) + AO2MObeta(ploop_k,ploop_i)*AO2MObeta(ploop_k,ploop_j)
                    end do
                end do
            end do
            AO2MO = cmplx(0.0_dp,0.0_dp,dp)
            do ploop_i = 1, basis_dimension
                do ploop_j = 1, basis_dimension
                    AO2MO(2*ploop_i-1,2*ploop_j-1) = AO2MOalpha(ploop_i,ploop_j)
                    AO2MO(2*ploop_i,2*ploop_j) = AO2MObeta(ploop_i,ploop_j)
                end do
            end do
            deallocate(AO2MOalpha)
            deallocate(AO2MObeta)
        else
        	rou_m = cmplx(0.0_dp,0.0_dp,dp)
            do ploop_i = 1, 2*basis_dimension
                do ploop_j = 1, 2*basis_dimension
                    do ploop_k = 1, electron_count
                        rou_m(ploop_i,ploop_j) = rou_m(ploop_i,ploop_j) + AO2MO(ploop_k,ploop_i)*conjg(AO2MO(ploop_k,ploop_j))
                    end do
                end do
            end do
            call zgemm( 'N', 'N', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), exS0_5, 2*basis_dimension, rou_m, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), AO2MO, 2*basis_dimension)
            call zgemm( 'N', 'N', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), AO2MO, 2*basis_dimension, exS0_5, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), rou_m, 2*basis_dimension)
            !--------------------------<debug>--------------------------
            !do loop_j = 1, 2*basis_dimension
            !    write(60,*) real(rou_m(loop_j, 1)), real(rou_m(loop_j, 2))
            !end do
            !--------------------------<debug>--------------------------
        end if
    end subroutine assign_rou
    
!------------------------------------------------------------
! construct one electron Fock operator
    subroutine Fock1e()
        integer :: floop_i, floop_j, floop_k, floop_l                   ! floop is loop variables only for subroutine Fock1e
        allocate(Fock1(2*basis_dimension,2*basis_dimension))
        allocate(exi_T_j(2*basis_dimension,2*basis_dimension))
        allocate(exi_V_j(2*basis_dimension,2*basis_dimension))
        if (DKH_order == 0) then
            Fock1 = cmplx(0.0_dp,0.0_dp,dp)
            exi_T_j = cmplx(0.0_dp,0.0_dp,dp)
            exi_V_j = cmplx(0.0_dp,0.0_dp,dp)
            do floop_i = 1, basis_dimension
                do floop_j = 1, basis_dimension
                    Fock1(2*floop_i-1,2*floop_j-1) = i_p2_j(floop_i,floop_j)/2.0_dp + i_V_j(floop_i,floop_j)
                    Fock1(2*floop_i,2*floop_j) = i_p2_j(floop_i,floop_j)/2.0_dp + i_V_j(floop_i,floop_j)
                    exi_T_j(2*floop_i-1,2*floop_j-1) = i_p2_j(floop_i,floop_j)/2.0_dp
                    exi_T_j(2*floop_i,2*floop_j) = i_p2_j(floop_i,floop_j)/2.0_dp
                    exi_V_j(2*floop_i-1,2*floop_j-1) = i_V_j(floop_i,floop_j)
                    exi_V_j(2*floop_i,2*floop_j) = i_V_j(floop_i,floop_j)
                end do
            end do
        else if (DKH_order == 2) then
            allocate(oper1(basis_dimension,basis_dimension))
            allocate(oper2(basis_dimension,basis_dimension))
            allocate(oper3(2*basis_dimension,2*basis_dimension))
            allocate(oper4(2*basis_dimension,2*basis_dimension))
            allocate(oper5(2*basis_dimension,2*basis_dimension))
            allocate(Ap(basis_dimension,basis_dimension))
            allocate(ApRp(basis_dimension,basis_dimension))                 ! ApRp = RpAp
            allocate(SRp(2*basis_dimension,2*basis_dimension))
            allocate(ARVRA(2*basis_dimension,2*basis_dimension))
            allocate(AVA(2*basis_dimension,2*basis_dimension))
            allocate(ARVeRA(2*basis_dimension,2*basis_dimension))
            allocate(AVeA(2*basis_dimension,2*basis_dimension))
            allocate(Ve(basis_dimension,basis_dimension))
            allocate(pxVepx(basis_dimension,basis_dimension))
            allocate(pyVepy(basis_dimension,basis_dimension))
            allocate(pzVepz(basis_dimension,basis_dimension))
            allocate(pxVepy(basis_dimension,basis_dimension))
            allocate(pyVepx(basis_dimension,basis_dimension))
            allocate(pxVepz(basis_dimension,basis_dimension))
            allocate(pzVepx(basis_dimension,basis_dimension))
            allocate(pyVepz(basis_dimension,basis_dimension))
            allocate(pzVepy(basis_dimension,basis_dimension))
            allocate(exAO2p2(2*basis_dimension,2*basis_dimension))
            allocate(exSOC(2*basis_dimension,2*basis_dimension))
            ARVRA = cmplx(0.0_dp,0.0_dp,dp)
            AVA = cmplx(0.0_dp,0.0_dp,dp)
            ARVeRA = cmplx(0.0_dp,0.0_dp,dp)
            AVeA = cmplx(0.0_dp,0.0_dp,dp)
            Fock1 = cmplx(0.0_dp,0.0_dp,dp)
            exi_T_j = cmplx(0.0_dp,0.0_dp,dp)
            exi_V_j = cmplx(0.0_dp,0.0_dp,dp)
            exSOC = cmplx(0.0_dp,0.0_dp,dp)
            Ap = 0.0_dp
            do floop_i=1,basis_dimension
                Ap(floop_i,floop_i) = sqrt((sqrt(evl_p2(floop_i)/(speed_light*speed_light) + 1.0_dp) + 1.0_dp)/(2.0_dp*sqrt(evl_p2(floop_i)/(speed_light*speed_light) + 1.0_dp)))
            end do
            ApRp = 0.0_dp
            do floop_i=1,basis_dimension
                ApRp(floop_i,floop_i) = Ap(floop_i,floop_i)/(sqrt(evl_p2(floop_i) + speed_light*speed_light) + speed_light)
            end do
            SRp = cmplx(0.0_dp,0.0_dp,dp)
            do floop_i=1,basis_dimension
                SRp(2*floop_i-1,2*floop_i-1) = cmplx(1.0_dp + evl_p2(floop_i)/(4.0_dp*speed_light*speed_light) + QED_rad,0.0_dp,dp)
                SRp(2*floop_i,2*floop_i) = cmplx(1.0_dp + evl_p2(floop_i)/(4.0_dp*speed_light*speed_light) + QED_rad,0.0_dp,dp)
            end do
            ! Ap V Ap
            call dgemm( 'N', 'T', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, Ap, basis_dimension, AO2p2, basis_dimension, 0.0_dp, oper2, basis_dimension)
            call dgemm( 'N', 'N', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, oper2, basis_dimension, i_V_j, basis_dimension, 0.0_dp, oper1, basis_dimension)
            call dgemm( 'N', 'N', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, oper1, basis_dimension, AO2p2, basis_dimension, 0.0_dp, oper2, basis_dimension)
            call dgemm( 'N', 'N', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, oper2, basis_dimension, Ap, basis_dimension, 0.0_dp, oper1, basis_dimension)
            do floop_i=1, basis_dimension
                do floop_j=1, basis_dimension
                    AVA(2*floop_i-1,2*floop_j-1) = AVA(2*floop_i-1,2*floop_j-1) + cmplx(oper1(floop_i,floop_j),0.0_dp,dp)
                    AVA(2*floop_i,2*floop_j) = AVA(2*floop_i,2*floop_j) + cmplx(oper1(floop_i,floop_j),0.0_dp,dp)
                    exi_V_j(2*floop_i-1,2*floop_j-1) = i_V_j(floop_i,floop_j)
                    exi_V_j(2*floop_i,2*floop_j) = i_V_j(floop_i,floop_j)
                end do
            end do
            ! ApRp pxVpx+pyVpy+pzVpz ApRp
            call dgemm( 'N', 'T', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, ApRp, basis_dimension, AO2p2, basis_dimension, 0.0_dp, oper2, basis_dimension)
            call dgemm( 'N', 'N', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, oper2, basis_dimension, i_pxVpx_j+i_pyVpy_j+i_pzVpz_j, basis_dimension, 0.0_dp, oper1, basis_dimension)
            call dgemm( 'N', 'N', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, oper1, basis_dimension, AO2p2, basis_dimension, 0.0_dp, oper2, basis_dimension)
            call dgemm( 'N', 'N', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, oper2, basis_dimension, ApRp, basis_dimension, 0.0_dp, oper1, basis_dimension)
            do floop_i=1, basis_dimension
                do floop_j=1, basis_dimension
                    ARVRA(2*floop_i-1,2*floop_j-1) = ARVRA(2*floop_i-1,2*floop_j-1) + cmplx(oper1(floop_i,floop_j),0.0_dp,dp)
                    ARVRA(2*floop_i,2*floop_j) = ARVRA(2*floop_i,2*floop_j) + cmplx(oper1(floop_i,floop_j),0.0_dp,dp)
                end do
            end do
            ! ApRp pxVpy-pyVpx ApRp
            call dgemm( 'N', 'T', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, ApRp, basis_dimension, AO2p2, basis_dimension, 0.0_dp, oper2, basis_dimension)
            call dgemm( 'N', 'N', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, oper2, basis_dimension, i_pxVpy_j-i_pyVpx_j, basis_dimension, 0.0_dp, oper1, basis_dimension)
            call dgemm( 'N', 'N', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, oper1, basis_dimension, AO2p2, basis_dimension, 0.0_dp, oper2, basis_dimension)
            call dgemm( 'N', 'N', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, oper2, basis_dimension, ApRp, basis_dimension, 0.0_dp, oper1, basis_dimension)
            do floop_i=1, basis_dimension
                do floop_j=1, basis_dimension
                    ARVRA(2*floop_i-1,2*floop_j-1) = ARVRA(2*floop_i-1,2*floop_j-1) + cmplx(0.0_dp,oper1(floop_i,floop_j),dp)
                    ARVRA(2*floop_i,2*floop_j) = ARVRA(2*floop_i,2*floop_j) - cmplx(0.0_dp,oper1(floop_i,floop_j),dp)
                end do
            end do
            ! ApRp pzVpx-pxVpz ApRp
            call dgemm( 'N', 'T', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, ApRp, basis_dimension, AO2p2, basis_dimension, 0.0_dp, oper2, basis_dimension)
            call dgemm( 'N', 'N', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, oper2, basis_dimension, i_pzVpx_j-i_pxVpz_j, basis_dimension, 0.0_dp, oper1, basis_dimension)
            call dgemm( 'N', 'N', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, oper1, basis_dimension, AO2p2, basis_dimension, 0.0_dp, oper2, basis_dimension)
            call dgemm( 'N', 'N', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, oper2, basis_dimension, ApRp, basis_dimension, 0.0_dp, oper1, basis_dimension)
            do floop_i=1, basis_dimension
                do floop_j=1, basis_dimension
                    ARVRA(2*floop_i,2*floop_j-1) = ARVRA(2*floop_i,2*floop_j-1) - cmplx(oper1(floop_i,floop_j),0.0_dp,dp)
                    ARVRA(2*floop_i-1,2*floop_j) = ARVRA(2*floop_i-1,2*floop_j) + cmplx(oper1(floop_i,floop_j),0.0_dp,dp)
                end do
            end do
            ! ApRp pyVpz-pzVpy ApRp
            call dgemm( 'N', 'T', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, ApRp, basis_dimension, AO2p2, basis_dimension, 0.0_dp, oper2, basis_dimension)
            call dgemm( 'N', 'N', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, oper2, basis_dimension, i_pyVpz_j-i_pzVpy_j, basis_dimension, 0.0_dp, oper1, basis_dimension)
            call dgemm( 'N', 'N', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, oper1, basis_dimension, AO2p2, basis_dimension, 0.0_dp, oper2, basis_dimension)
            call dgemm( 'N', 'N', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, oper2, basis_dimension, ApRp, basis_dimension, 0.0_dp, oper1, basis_dimension)
            do floop_i=1, basis_dimension
                do floop_j=1, basis_dimension
                    ARVRA(2*floop_i,2*floop_j-1) = ARVRA(2*floop_i,2*floop_j-1) + cmplx(0.0_dp,oper1(floop_i,floop_j),dp)
                    ARVRA(2*floop_i-1,2*floop_j) = ARVRA(2*floop_i-1,2*floop_j) + cmplx(0.0_dp,oper1(floop_i,floop_j),dp)
                end do
            end do
            Fock1 = Fock1 + AVA
            if (SRTP_type) then
                call zgemm( 'N', 'N', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), SRp, 2*basis_dimension, ARVRA, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), oper3, 2*basis_dimension)
                call zgemm( 'N', 'N', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), oper3, 2*basis_dimension, SRp, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), oper4, 2*basis_dimension)
                Fock1 = Fock1 + oper4
                exSOC = oper4
            else
                Fock1 = Fock1 + ARVRA
                exSOC = ARVRA
            end if
            if (SRTP_type) then
                allocate(exSR(2*basis_dimension,2*basis_dimension))
                exSR = cmplx(0.0_dp,0.0_dp,dp)
                ! ApRp px3Vpx+py3Vpy+pz3Vpz ApRp
                call dgemm( 'N', 'T', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, ApRp, basis_dimension, AO2p2, basis_dimension, 0.0_dp, oper2, basis_dimension)
                call dgemm( 'N', 'N', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, oper2, basis_dimension, i_px3Vpx_j+i_py3Vpy_j+i_pz3Vpz_j, basis_dimension, 0.0_dp, oper1, basis_dimension)
                call dgemm( 'N', 'N', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, oper1, basis_dimension, AO2p2, basis_dimension, 0.0_dp, oper2, basis_dimension)
                call dgemm( 'N', 'N', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, oper2, basis_dimension, ApRp, basis_dimension, 0.0_dp, oper1, basis_dimension)
                do floop_i=1, basis_dimension
                    do floop_j=1, basis_dimension
                        Fock1(2*floop_i-1,2*floop_j-1) = Fock1(2*floop_i-1,2*floop_j-1) - cmplx(oper1(floop_i,floop_j) / (2.0_dp*speed_light*speed_light),0.0_dp,dp)
                        Fock1(2*floop_i,2*floop_j) = Fock1(2*floop_i,2*floop_j) - cmplx(oper1(floop_i,floop_j) / (2.0_dp*speed_light*speed_light),0.0_dp,dp)
                        exSR(2*floop_i-1,2*floop_j-1) = exSR(2*floop_i-1,2*floop_j-1) - cmplx(oper1(floop_i,floop_j) / (2.0_dp*speed_light*speed_light),0.0_dp,dp)
                        exSR(2*floop_i,2*floop_j) = exSR(2*floop_i,2*floop_j) - cmplx(oper1(floop_i,floop_j) / (2.0_dp*speed_light*speed_light),0.0_dp,dp)
                    end do
                end do
                ! ApRp px3Vpy-py3Vpx ApRp
                call dgemm( 'N', 'T', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, ApRp, basis_dimension, AO2p2, basis_dimension, 0.0_dp, oper2, basis_dimension)
                call dgemm( 'N', 'N', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, oper2, basis_dimension, i_px3Vpy_j-i_py3Vpx_j, basis_dimension, 0.0_dp, oper1, basis_dimension)
                call dgemm( 'N', 'N', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, oper1, basis_dimension, AO2p2, basis_dimension, 0.0_dp, oper2, basis_dimension)
                call dgemm( 'N', 'N', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, oper2, basis_dimension, ApRp, basis_dimension, 0.0_dp, oper1, basis_dimension)
                do floop_i=1, basis_dimension
                    do floop_j=1, basis_dimension
                        Fock1(2*floop_i-1,2*floop_j-1) = Fock1(2*floop_i-1,2*floop_j-1) - cmplx(0.0_dp,oper1(floop_i,floop_j) / (2.0_dp*speed_light*speed_light),dp)
                        Fock1(2*floop_i,2*floop_j) = Fock1(2*floop_i,2*floop_j) + cmplx(0.0_dp,oper1(floop_i,floop_j) / (2.0_dp*speed_light*speed_light),dp)
                        exSR(2*floop_i-1,2*floop_j-1) = exSR(2*floop_i-1,2*floop_j-1) - cmplx(0.0_dp,oper1(floop_i,floop_j) / (2.0_dp*speed_light*speed_light),dp)
                        exSR(2*floop_i,2*floop_j) = exSR(2*floop_i,2*floop_j) + cmplx(0.0_dp,oper1(floop_i,floop_j) / (2.0_dp*speed_light*speed_light),dp)
                    end do
                end do
                ! ApRp pz3Vpx-px3Vpz ApRp
                call dgemm( 'N', 'T', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, ApRp, basis_dimension, AO2p2, basis_dimension, 0.0_dp, oper2, basis_dimension)
                call dgemm( 'N', 'N', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, oper2, basis_dimension, i_pz3Vpx_j-i_px3Vpz_j, basis_dimension, 0.0_dp, oper1, basis_dimension)
                call dgemm( 'N', 'N', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, oper1, basis_dimension, AO2p2, basis_dimension, 0.0_dp, oper2, basis_dimension)
                call dgemm( 'N', 'N', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, oper2, basis_dimension, ApRp, basis_dimension, 0.0_dp, oper1, basis_dimension)
                do floop_i=1, basis_dimension
                    do floop_j=1, basis_dimension
                        Fock1(2*floop_i,2*floop_j-1) = Fock1(2*floop_i,2*floop_j-1) + cmplx(oper1(floop_i,floop_j) / (2.0_dp*speed_light*speed_light),0.0_dp,dp)
                        Fock1(2*floop_i-1,2*floop_j) = Fock1(2*floop_i-1,2*floop_j) - cmplx(oper1(floop_i,floop_j) / (2.0_dp*speed_light*speed_light),0.0_dp,dp)
                        exSR(2*floop_i,2*floop_j-1) = exSR(2*floop_i,2*floop_j-1) + cmplx(oper1(floop_i,floop_j) / (2.0_dp*speed_light*speed_light),0.0_dp,dp)
                        exSR(2*floop_i-1,2*floop_j) = exSR(2*floop_i-1,2*floop_j) - cmplx(oper1(floop_i,floop_j) / (2.0_dp*speed_light*speed_light),0.0_dp,dp)
                    end do
                end do
                ! ApRp py3Vpz-pz3Vpy ApRp
                call dgemm( 'N', 'T', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, ApRp, basis_dimension, AO2p2, basis_dimension, 0.0_dp, oper2, basis_dimension)
                call dgemm( 'N', 'N', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, oper2, basis_dimension, i_py3Vpz_j-i_pz3Vpy_j, basis_dimension, 0.0_dp, oper1, basis_dimension)
                call dgemm( 'N', 'N', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, oper1, basis_dimension, AO2p2, basis_dimension, 0.0_dp, oper2, basis_dimension)
                call dgemm( 'N', 'N', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, oper2, basis_dimension, ApRp, basis_dimension, 0.0_dp, oper1, basis_dimension)
                do floop_i=1, basis_dimension
                    do floop_j=1, basis_dimension
                        Fock1(2*floop_i,2*floop_j-1) = Fock1(2*floop_i,2*floop_j-1) - cmplx(0.0_dp,oper1(floop_i,floop_j) / (2.0_dp*speed_light*speed_light),dp)
                        Fock1(2*floop_i-1,2*floop_j) = Fock1(2*floop_i-1,2*floop_j) - cmplx(0.0_dp,oper1(floop_i,floop_j) / (2.0_dp*speed_light*speed_light),dp)
                        exSR(2*floop_i,2*floop_j-1) = exSR(2*floop_i,2*floop_j-1) - cmplx(0.0_dp,oper1(floop_i,floop_j) / (2.0_dp*speed_light*speed_light),dp)
                        exSR(2*floop_i-1,2*floop_j) = exSR(2*floop_i-1,2*floop_j) - cmplx(0.0_dp,oper1(floop_i,floop_j) / (2.0_dp*speed_light*speed_light),dp)
                    end do
                end do
                ! ApRp pxVpx3+pyVpy3+pzVpz3 ApRp
                call dgemm( 'N', 'T', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, ApRp, basis_dimension, AO2p2, basis_dimension, 0.0_dp, oper2, basis_dimension)
                call dgemm( 'N', 'T', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, oper2, basis_dimension, i_px3Vpx_j+i_py3Vpy_j+i_pz3Vpz_j, basis_dimension, 0.0_dp, oper1, basis_dimension)
                call dgemm( 'N', 'N', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, oper1, basis_dimension, AO2p2, basis_dimension, 0.0_dp, oper2, basis_dimension)
                call dgemm( 'N', 'N', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, oper2, basis_dimension, ApRp, basis_dimension, 0.0_dp, oper1, basis_dimension)
                do floop_i=1, basis_dimension
                    do floop_j=1, basis_dimension
                        Fock1(2*floop_i-1,2*floop_j-1) = Fock1(2*floop_i-1,2*floop_j-1) - cmplx(oper1(floop_i,floop_j) / (2.0_dp*speed_light*speed_light),0.0_dp,dp)
                        Fock1(2*floop_i,2*floop_j) = Fock1(2*floop_i,2*floop_j) - cmplx(oper1(floop_i,floop_j) / (2.0_dp*speed_light*speed_light),0.0_dp,dp)
                        exSR(2*floop_i-1,2*floop_j-1) = exSR(2*floop_i-1,2*floop_j-1) - cmplx(oper1(floop_i,floop_j) / (2.0_dp*speed_light*speed_light),0.0_dp,dp)
                        exSR(2*floop_i,2*floop_j) = exSR(2*floop_i,2*floop_j) - cmplx(oper1(floop_i,floop_j) / (2.0_dp*speed_light*speed_light),0.0_dp,dp)
                    end do
                end do
                ! ApRp pxVpy3-pyVpx3 ApRp
                call dgemm( 'N', 'T', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, ApRp, basis_dimension, AO2p2, basis_dimension, 0.0_dp, oper2, basis_dimension)
                call dgemm( 'N', 'T', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, oper2, basis_dimension, i_py3Vpx_j-i_px3Vpy_j, basis_dimension, 0.0_dp, oper1, basis_dimension)
                call dgemm( 'N', 'N', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, oper1, basis_dimension, AO2p2, basis_dimension, 0.0_dp, oper2, basis_dimension)
                call dgemm( 'N', 'N', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, oper2, basis_dimension, ApRp, basis_dimension, 0.0_dp, oper1, basis_dimension)
                do floop_i=1, basis_dimension
                    do floop_j=1, basis_dimension
                        Fock1(2*floop_i-1,2*floop_j-1) = Fock1(2*floop_i-1,2*floop_j-1) - cmplx(0.0_dp,oper1(floop_i,floop_j) / (2.0_dp*speed_light*speed_light),dp)
                        Fock1(2*floop_i,2*floop_j) = Fock1(2*floop_i,2*floop_j) + cmplx(0.0_dp,oper1(floop_i,floop_j) / (2.0_dp*speed_light*speed_light),dp)
                        exSR(2*floop_i-1,2*floop_j-1) = exSR(2*floop_i-1,2*floop_j-1) - cmplx(0.0_dp,oper1(floop_i,floop_j) / (2.0_dp*speed_light*speed_light),dp)
                        exSR(2*floop_i,2*floop_j) = exSR(2*floop_i,2*floop_j) + cmplx(0.0_dp,oper1(floop_i,floop_j) / (2.0_dp*speed_light*speed_light),dp)
                    end do
                end do
                ! ApRp pxVpz3-pzVpx3 ApRp
                call dgemm( 'N', 'T', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, ApRp, basis_dimension, AO2p2, basis_dimension, 0.0_dp, oper2, basis_dimension)
                call dgemm( 'N', 'T', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, oper2, basis_dimension, i_px3Vpz_j-i_pz3Vpx_j, basis_dimension, 0.0_dp, oper1, basis_dimension)
                call dgemm( 'N', 'N', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, oper1, basis_dimension, AO2p2, basis_dimension, 0.0_dp, oper2, basis_dimension)
                call dgemm( 'N', 'N', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, oper2, basis_dimension, ApRp, basis_dimension, 0.0_dp, oper1, basis_dimension)
                do floop_i=1, basis_dimension
                    do floop_j=1, basis_dimension
                        Fock1(2*floop_i,2*floop_j-1) = Fock1(2*floop_i,2*floop_j-1) + cmplx(oper1(floop_i,floop_j) / (2.0_dp*speed_light*speed_light),0.0_dp,dp)
                        Fock1(2*floop_i-1,2*floop_j) = Fock1(2*floop_i-1,2*floop_j) - cmplx(oper1(floop_i,floop_j) / (2.0_dp*speed_light*speed_light),0.0_dp,dp)
                        exSR(2*floop_i,2*floop_j-1) = exSR(2*floop_i,2*floop_j-1) + cmplx(oper1(floop_i,floop_j) / (2.0_dp*speed_light*speed_light),0.0_dp,dp)
                        exSR(2*floop_i-1,2*floop_j) = exSR(2*floop_i-1,2*floop_j) - cmplx(oper1(floop_i,floop_j) / (2.0_dp*speed_light*speed_light),0.0_dp,dp)
                    end do
                end do
                ! ApRp pyVpz3-pzVpy3 ApRp
                call dgemm( 'N', 'T', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, ApRp, basis_dimension, AO2p2, basis_dimension, 0.0_dp, oper2, basis_dimension)
                call dgemm( 'N', 'T', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, oper2, basis_dimension, i_pz3Vpy_j-i_py3Vpz_j, basis_dimension, 0.0_dp, oper1, basis_dimension)
                call dgemm( 'N', 'N', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, oper1, basis_dimension, AO2p2, basis_dimension, 0.0_dp, oper2, basis_dimension)
                call dgemm( 'N', 'N', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, oper2, basis_dimension, ApRp, basis_dimension, 0.0_dp, oper1, basis_dimension)
                do floop_i=1, basis_dimension
                    do floop_j=1, basis_dimension
                        Fock1(2*floop_i,2*floop_j-1) = Fock1(2*floop_i,2*floop_j-1) - cmplx(0.0_dp,oper1(floop_i,floop_j) / (2.0_dp*speed_light*speed_light),dp)
                        Fock1(2*floop_i-1,2*floop_j) = Fock1(2*floop_i-1,2*floop_j) - cmplx(0.0_dp,oper1(floop_i,floop_j) / (2.0_dp*speed_light*speed_light),dp)
                        exSR(2*floop_i,2*floop_j-1) = exSR(2*floop_i,2*floop_j-1) - cmplx(0.0_dp,oper1(floop_i,floop_j) / (2.0_dp*speed_light*speed_light),dp)
                        exSR(2*floop_i-1,2*floop_j) = exSR(2*floop_i-1,2*floop_j) - cmplx(0.0_dp,oper1(floop_i,floop_j) / (2.0_dp*speed_light*speed_light),dp)
                    end do
                end do
            end if
            do floop_i=1,basis_dimension
                do floop_j=1,basis_dimension
                    Ve(floop_i,floop_j) = i_V_j(floop_i,floop_j)/(speed_light * (sqrt(evl_p2(floop_i)+speed_light*speed_light) + sqrt(evl_p2(floop_j)+speed_light*speed_light)))
                end do
            end do
            do floop_i=1,basis_dimension
                do floop_j=1,basis_dimension
                    pxVepx(floop_i,floop_j) = i_pxVpx_j(floop_i,floop_j)/(speed_light * (sqrt(evl_p2(floop_i)+speed_light*speed_light) + sqrt(evl_p2(floop_j)+speed_light*speed_light)))
                end do
            end do
            do floop_i=1,basis_dimension
                do floop_j=1,basis_dimension
                    pyVepy(floop_i,floop_j) = i_pyVpy_j(floop_i,floop_j)/(speed_light * (sqrt(evl_p2(floop_i)+speed_light*speed_light) + sqrt(evl_p2(floop_j)+speed_light*speed_light)))
                end do
            end do
            do floop_i=1,basis_dimension
                do floop_j=1,basis_dimension
                    pzVepz(floop_i,floop_j) = i_pzVpz_j(floop_i,floop_j)/(speed_light * (sqrt(evl_p2(floop_i)+speed_light*speed_light) + sqrt(evl_p2(floop_j)+speed_light*speed_light)))
                end do
            end do
            do floop_i=1,basis_dimension
                do floop_j=1,basis_dimension
                    pxVepy(floop_i,floop_j) = i_pxVpy_j(floop_i,floop_j)/(speed_light * (sqrt(evl_p2(floop_i)+speed_light*speed_light) + sqrt(evl_p2(floop_j)+speed_light*speed_light)))
                end do
            end do
            do floop_i=1,basis_dimension
                do floop_j=1,basis_dimension
                    pyVepx(floop_i,floop_j) = i_pyVpx_j(floop_i,floop_j)/(speed_light * (sqrt(evl_p2(floop_i)+speed_light*speed_light) + sqrt(evl_p2(floop_j)+speed_light*speed_light)))
                end do
            end do
            do floop_i=1,basis_dimension
                do floop_j=1,basis_dimension
                    pxVepz(floop_i,floop_j) = i_pxVpz_j(floop_i,floop_j)/(speed_light * (sqrt(evl_p2(floop_i)+speed_light*speed_light) + sqrt(evl_p2(floop_j)+speed_light*speed_light)))
                end do
            end do
            do floop_i=1,basis_dimension
                do floop_j=1,basis_dimension
                    pzVepx(floop_i,floop_j) = i_pzVpx_j(floop_i,floop_j)/(speed_light * (sqrt(evl_p2(floop_i)+speed_light*speed_light) + sqrt(evl_p2(floop_j)+speed_light*speed_light)))
                end do
            end do
            do floop_i=1,basis_dimension
                do floop_j=1,basis_dimension
                    pyVepz(floop_i,floop_j) = i_pyVpz_j(floop_i,floop_j)/(speed_light * (sqrt(evl_p2(floop_i)+speed_light*speed_light) + sqrt(evl_p2(floop_j)+speed_light*speed_light)))
                end do
            end do
            do floop_i=1,basis_dimension
                do floop_j=1,basis_dimension
                    pzVepy(floop_i,floop_j) = i_pzVpy_j(floop_i,floop_j)/(speed_light * (sqrt(evl_p2(floop_i)+speed_light*speed_light) + sqrt(evl_p2(floop_j)+speed_light*speed_light)))
                end do
            end do
            ! Ap Ve Ap
            call dgemm( 'N', 'T', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, Ap, basis_dimension, AO2p2, basis_dimension, 0.0_dp, oper2, basis_dimension)
            call dgemm( 'N', 'N', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, oper2, basis_dimension, Ve, basis_dimension, 0.0_dp, oper1, basis_dimension)
            call dgemm( 'N', 'N', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, oper1, basis_dimension, AO2p2, basis_dimension, 0.0_dp, oper2, basis_dimension)
            call dgemm( 'N', 'N', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, oper2, basis_dimension, Ap, basis_dimension, 0.0_dp, oper1, basis_dimension)
            do floop_i=1, basis_dimension
                do floop_j=1, basis_dimension
                    AVeA(2*floop_i-1,2*floop_j-1) = AVeA(2*floop_i-1,2*floop_j-1) + cmplx(oper1(floop_i,floop_j),0.0_dp,dp)
                    AVeA(2*floop_i,2*floop_j) = AVeA(2*floop_i,2*floop_j) + cmplx(oper1(floop_i,floop_j),0.0_dp,dp)
                end do
            end do
            ! ApRp pxVepx+pyVepy+pzVepz ApRp
            call dgemm( 'N', 'T', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, ApRp, basis_dimension, AO2p2, basis_dimension, 0.0_dp, oper2, basis_dimension)
            call dgemm( 'N', 'N', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, oper2, basis_dimension, pxVepx+pyVepy+pzVepz, basis_dimension, 0.0_dp, oper1, basis_dimension)
            call dgemm( 'N', 'N', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, oper1, basis_dimension, AO2p2, basis_dimension, 0.0_dp, oper2, basis_dimension)
            call dgemm( 'N', 'N', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, oper2, basis_dimension, ApRp, basis_dimension, 0.0_dp, oper1, basis_dimension)
            do floop_i=1, basis_dimension
                do floop_j=1, basis_dimension
                    ARVeRA(2*floop_i-1,2*floop_j-1) = ARVeRA(2*floop_i-1,2*floop_j-1) + cmplx(oper1(floop_i,floop_j),0.0_dp,dp)
                    ARVeRA(2*floop_i,2*floop_j) = ARVeRA(2*floop_i,2*floop_j) + cmplx(oper1(floop_i,floop_j),0.0_dp,dp)
                end do
            end do
            ! ApRp pxVepy-pyVepx ApRp
            call dgemm( 'N', 'T', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, ApRp, basis_dimension, AO2p2, basis_dimension, 0.0_dp, oper2, basis_dimension)
            call dgemm( 'N', 'N', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, oper2, basis_dimension, pxVepy-pyVepx, basis_dimension, 0.0_dp, oper1, basis_dimension)
            call dgemm( 'N', 'N', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, oper1, basis_dimension, AO2p2, basis_dimension, 0.0_dp, oper2, basis_dimension)
            call dgemm( 'N', 'N', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, oper2, basis_dimension, ApRp, basis_dimension, 0.0_dp, oper1, basis_dimension)
            do floop_i=1, basis_dimension
                do floop_j=1, basis_dimension
                    ARVeRA(2*floop_i-1,2*floop_j-1) = ARVeRA(2*floop_i-1,2*floop_j-1) + cmplx(0.0_dp,oper1(floop_i,floop_j),dp)
                    ARVeRA(2*floop_i,2*floop_j) = ARVeRA(2*floop_i,2*floop_j) - cmplx(0.0_dp,oper1(floop_i,floop_j),dp)
                end do
            end do
            ! ApRp pzVepx-pxVepz ApRp
            call dgemm( 'N', 'T', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, ApRp, basis_dimension, AO2p2, basis_dimension, 0.0_dp, oper2, basis_dimension)
            call dgemm( 'N', 'N', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, oper2, basis_dimension, pzVepx-pxVepz, basis_dimension, 0.0_dp, oper1, basis_dimension)
            call dgemm( 'N', 'N', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, oper1, basis_dimension, AO2p2, basis_dimension, 0.0_dp, oper2, basis_dimension)
            call dgemm( 'N', 'N', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, oper2, basis_dimension, ApRp, basis_dimension, 0.0_dp, oper1, basis_dimension)
            do floop_i=1, basis_dimension
                do floop_j=1, basis_dimension
                    ARVeRA(2*floop_i,2*floop_j-1) = ARVeRA(2*floop_i,2*floop_j-1) - cmplx(oper1(floop_i,floop_j),0.0_dp,dp)
                    ARVeRA(2*floop_i-1,2*floop_j) = ARVeRA(2*floop_i-1,2*floop_j) + cmplx(oper1(floop_i,floop_j),0.0_dp,dp)
                end do
            end do
            ! ApRp pyVepz-pzVepy ApRp
            call dgemm( 'N', 'T', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, ApRp, basis_dimension, AO2p2, basis_dimension, 0.0_dp, oper2, basis_dimension)
            call dgemm( 'N', 'N', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, oper2, basis_dimension, pyVepz-pzVepy, basis_dimension, 0.0_dp, oper1, basis_dimension)
            call dgemm( 'N', 'N', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, oper1, basis_dimension, AO2p2, basis_dimension, 0.0_dp, oper2, basis_dimension)
            call dgemm( 'N', 'N', basis_dimension, basis_dimension, basis_dimension, 1.0_dp, oper2, basis_dimension, ApRp, basis_dimension, 0.0_dp, oper1, basis_dimension)
            do floop_i=1, basis_dimension
                do floop_j=1, basis_dimension
                    ARVeRA(2*floop_i,2*floop_j-1) = ARVeRA(2*floop_i,2*floop_j-1) + cmplx(0.0_dp,oper1(floop_i,floop_j),dp)
                    ARVeRA(2*floop_i-1,2*floop_j) = ARVeRA(2*floop_i-1,2*floop_j) + cmplx(0.0_dp,oper1(floop_i,floop_j),dp)
                end do
            end do
        
            oper5 = cmplx(0.0_dp,0.0_dp,dp)
            do floop_i=1, basis_dimension
                oper5(2*floop_i-1,2*floop_i-1) = cmplx(0.5_dp,0.0_dp,dp)
                oper5(2*floop_i,2*floop_i) = cmplx(0.5_dp,0.0_dp,dp)
            end do
            call zgemm( 'N', 'N', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), ARVeRA, 2*basis_dimension, oper5, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), oper3, 2*basis_dimension)
            call zgemm( 'N', 'N', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), oper3, 2*basis_dimension, AVA, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), oper4, 2*basis_dimension)
            Fock1 = Fock1 - oper4
            call zgemm( 'N', 'N', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), ARVRA, 2*basis_dimension, oper5, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), oper3, 2*basis_dimension)
            call zgemm( 'N', 'N', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), oper3, 2*basis_dimension, AVeA, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), oper4, 2*basis_dimension)
            Fock1 = Fock1 - oper4
            call zgemm( 'N', 'N', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), AVeA, 2*basis_dimension, oper5, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), oper3, 2*basis_dimension)
            call zgemm( 'N', 'N', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), oper3, 2*basis_dimension, ARVRA, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), oper4, 2*basis_dimension)
            Fock1 = Fock1 - oper4
            call zgemm( 'N', 'N', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), AVA, 2*basis_dimension, oper5, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), oper3, 2*basis_dimension)
            call zgemm( 'N', 'N', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), oper3, 2*basis_dimension, ARVeRA, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), oper4, 2*basis_dimension)
            Fock1 = Fock1 - oper4
            do floop_i=1, basis_dimension
                oper5(2*floop_i-1,2*floop_i-1) = cmplx(0.5_dp,0.0_dp,dp) * cmplx((sqrt(evl_p2(floop_i) + speed_light*speed_light) + speed_light)**2 / evl_p2(floop_i),0.0_dp,dp)
                oper5(2*floop_i,2*floop_i) = cmplx(0.5_dp,0.0_dp,dp) * cmplx((sqrt(evl_p2(floop_i) + speed_light*speed_light) + speed_light)**2 / evl_p2(floop_i),0.0_dp,dp)
            end do
            call zgemm( 'N', 'N', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), ARVeRA, 2*basis_dimension, oper5, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), oper3, 2*basis_dimension)
            call zgemm( 'N', 'N', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), oper3, 2*basis_dimension, ARVRA, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), oper4, 2*basis_dimension)
            Fock1 = Fock1 + oper4
            call zgemm( 'N', 'N', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), ARVRA, 2*basis_dimension, oper5, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), oper3, 2*basis_dimension)
            call zgemm( 'N', 'N', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), oper3, 2*basis_dimension, ARVeRA, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), oper4, 2*basis_dimension)
            Fock1 = Fock1 + oper4
            do floop_i=1, basis_dimension
                oper5(2*floop_i-1,2*floop_i-1) = cmplx(evl_p2(floop_i) / (2.0_dp*(sqrt(evl_p2(floop_i)+speed_light*speed_light) + speed_light)**2),0.0_dp,dp)
                oper5(2*floop_i,2*floop_i) = cmplx(evl_p2(floop_i) / (2.0_dp * speed_light * sqrt(evl_p2(floop_i) + speed_light*speed_light)),0.0_dp,dp)
            end do
            call zgemm( 'N', 'N', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), AVeA, 2*basis_dimension, oper5, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), oper3, 2*basis_dimension)
            call zgemm( 'N', 'N', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), oper3, 2*basis_dimension, AVA, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), oper4, 2*basis_dimension)
            Fock1 = Fock1 + oper4
            call zgemm( 'N', 'N', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), AVA, 2*basis_dimension, oper5, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), oper3, 2*basis_dimension)
            call zgemm( 'N', 'N', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), oper3, 2*basis_dimension, AVeA, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), oper4, 2*basis_dimension)
            Fock1 = Fock1 + oper4
            do floop_i=1,basis_dimension
                Fock1(2*floop_i-1,2*floop_i-1) = Fock1(2*floop_i-1,2*floop_i-1) + cmplx(speed_light*sqrt(evl_p2(floop_i) + speed_light*speed_light),0.0_dp,dp) - speed_light*speed_light
                Fock1(2*floop_i,2*floop_i) = Fock1(2*floop_i,2*floop_i) + cmplx(speed_light*sqrt(evl_p2(floop_i) + speed_light*speed_light),0.0_dp,dp) - speed_light*speed_light
                exi_T_j(2*floop_i-1,2*floop_i-1) = cmplx(speed_light*sqrt(evl_p2(floop_i) + speed_light*speed_light),0.0_dp,dp) - speed_light*speed_light
                exi_T_j(2*floop_i,2*floop_i) = cmplx(speed_light*sqrt(evl_p2(floop_i) + speed_light*speed_light),0.0_dp,dp) - speed_light*speed_light
            end do
            ! Transform from p^2 eigenbasis to orthogonal normalized AO basis
            exAO2p2 = cmplx(0.0_dp,0.0_dp,dp)
            do floop_i=1,basis_dimension
                do floop_j=1,basis_dimension
                    exAO2p2(2*floop_i-1,2*floop_j-1) = cmplx(AO2p2(floop_i,floop_j),0.0_dp,dp)
                    exAO2p2(2*floop_i,2*floop_j) = cmplx(AO2p2(floop_i,floop_j),0.0_dp,dp)
                end do
            end do
            call zgemm( 'N', 'N', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), exAO2p2, 2*basis_dimension, Fock1, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), oper3, 2*basis_dimension)
            call zgemm( 'N', 'T', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), oper3, 2*basis_dimension, exAO2p2, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), Fock1, 2*basis_dimension)
            call zgemm( 'N', 'N', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), exAO2p2, 2*basis_dimension, exi_T_j, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), oper3, 2*basis_dimension)
            call zgemm( 'N', 'T', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), oper3, 2*basis_dimension, exAO2p2, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), exi_T_j, 2*basis_dimension)
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
! construct two electron Fock operator,
! "direct" calculation for 2 electron Fock, which will avoid memory overflow and large amount of disk R&W
! but 2 electron integral will be calculated in each round of SCF
    subroutine Fock2e()
        integer :: i                                                ! for parallel computation, dloop_i is replaced by i
        integer :: dloop_i, dloop_j, dloop_k, dloop_l               ! dloop is loop variables only for Fock2e routine
        integer :: dloop_m, dloop_n, dloop_o, dloop_p               ! dloop is loop variables only for Fock2e routine
        real(dp) :: integral
        !----------------------------------
        integer :: contraction_i                                    ! contraction of atom_i, shell_i
        integer :: atom_i                                           ! which atom is the i^th component of |χi>
        integer :: shell_i                                          ! which shell is the i^th component of |χi>
        integer :: L_i                                              ! angular quantum number of |χi>
        integer :: M_i                                              ! magnetic quantum number of |χi>
        !----------------------------------
        integer :: contraction_j                                    ! contraction of atom_j, shell_j
        integer :: atom_j                                           ! which atom is the j_th component of |χj>
        integer :: shell_j                                          ! which shell is the j_th component of |χj>
        integer :: L_j                                              ! angular quantum number of |χj>
        integer :: M_j                                              ! magnetic quantum number of |χj>
        !----------------------------------
        integer :: contraction_k                                    ! contraction of atom_k, shell_k
        integer :: atom_k                                           ! which atom is the i^th component of |χk>
        integer :: shell_k                                          ! which shell is the i^th component of |χk>
        integer :: L_k                                              ! angular quantum number of |χk>
        integer :: M_k                                              ! magnetic quantum number of |χk>
        !----------------------------------
        integer :: contraction_l                                    ! contraction of atom_l, shell_l
        integer :: atom_l                                           ! which atom is the i^th component of |χl>
        integer :: shell_l                                          ! which shell is the i^th component of |χl>
        integer :: L_l                                              ! angular quantum number of |χl>
        integer :: M_l                                              ! magnetic quantum number of |χl>
        real(dp) :: exponents_i(20)                                 ! exponents of |χi>
        real(dp) :: exponents_j(20)                                 ! exponents of |χj>
        real(dp) :: exponents_k(20)                                 ! exponents of |χk>
        real(dp) :: exponents_l(20)                                 ! exponents of |χl>
        real(dp) :: coefficient_i(20)                               ! coefficient of |χi>
        real(dp) :: coefficient_j(20)                               ! coefficient of |χj>
        real(dp) :: coefficient_k(20)                               ! coefficient of |χk>
        real(dp) :: coefficient_l(20)                               ! coefficient of |χl>
        real(dp) :: coordinate_i(3)                                 ! coordinate of center of |χi>
        real(dp) :: coordinate_j(3)                                 ! coordinate of center of |χj>
        real(dp) :: coordinate_k(3)                                 ! coordinate of center of |χk>
        real(dp) :: coordinate_l(3)                                 ! coordinate of center of |χl>
        real(dp) :: swintegral_mic
        integer :: swloop_i, swloop_j                               ! swloop is loop variables only for schwarz screening
        integer :: swloop_m, swloop_n, swloop_o, swloop_p           ! swloop is loop variables only for schwarz screening
        integer :: shell_start_i                                    ! start point of same angular momentum shells
        integer :: shell_start_j                                    ! start point of same angular momentum shells
        real(dp),allocatable :: swexponents_i(:)                    ! exponents of |χi> for schwarz screening
        real(dp),allocatable :: swexponents_j(:)                    ! exponents of |χj> for schwarz screening
        real(dp),allocatable :: swcoefficient_i(:)                  ! coefficient of |χi> for schwarz screening
        real(dp),allocatable :: swcoefficient_j(:)                  ! coefficient of |χj> for schwarz screening
        if (ndschwarz) then
        	ndschwarz = .false.
            write(60,'(a)') '      Schwarz screening of <ij||kl>'
            allocate(Fock2(2*basis_dimension,2*basis_dimension))
            allocate(swintegral(basis_dimension,basis_dimension))
            allocate(basis_inf(basis_dimension))
            swintegral = 0.0_dp
            swloop_i = 1
            atom_i = 1
            shell_i = 1
            shell_start_i = 1
            do while(swloop_i <= basis_dimension)
                if (shell_i > shell_in_element(molecular(atom_i) % atom_number)) then
                    shell_i = 1
                    atom_i = atom_i + 1
                end if
                contraction_i = atom_basis(molecular(atom_i) % basis_number + shell_i - 1) % contraction
                L_i = atom_basis(molecular(atom_i) % basis_number + shell_i - 1) % angular_quantum_number + 1
                M_i = swloop_i - shell_start_i + 1
				! prepare for openMP parallel computation
                basis_inf(swloop_i) % atom = atom_i
                basis_inf(swloop_i) % shell = shell_i
                basis_inf(swloop_i) % L = L_i
                basis_inf(swloop_i) % M = M_i
                allocate(swexponents_i(contraction_i))
                allocate(swcoefficient_i(contraction_i))
                swexponents_i = atom_basis(molecular(atom_i) % basis_number + shell_i - 1) % exponents
                do swloop_m=1,contraction_i
                    swcoefficient_i(swloop_m) = atom_basis(molecular(atom_i) % basis_number + shell_i - 1) % Ncoefficient(swloop_m,M_i)
                end do
                coordinate_i = molecular(atom_i) % nucleus_position
                swloop_j = swloop_i
                atom_j = atom_i
                shell_j = shell_i
                shell_start_j = shell_start_i
                do while(swloop_j <= basis_dimension)
                    if (shell_j > shell_in_element(molecular(atom_j) % atom_number)) then
                        shell_j = 1
                        atom_j = atom_j + 1
                    end if
                    contraction_j = atom_basis(molecular(atom_j) % basis_number + shell_j - 1) % contraction
                    L_j = atom_basis(molecular(atom_j) % basis_number + shell_j - 1) % angular_quantum_number + 1
                    M_j = swloop_j - shell_start_j + 1
                    allocate(swexponents_j(contraction_j))
                    allocate(swcoefficient_j(contraction_j))
                    swexponents_j = atom_basis(molecular(atom_j) % basis_number + shell_j - 1) % exponents
                    do swloop_m=1,contraction_j
                        swcoefficient_j(swloop_m) = atom_basis(molecular(atom_j) % basis_number + shell_j - 1) % Ncoefficient(swloop_m,M_j)
                    end do
                    coordinate_j = molecular(atom_j) % nucleus_position                       
                    do swloop_m = 1, contraction_i
                        do swloop_n = 1, contraction_j
                            do swloop_o = 1, contraction_i
                                do swloop_p = 1, contraction_j
							        if (L_i >= L_j) then
								        swintegral_mic = swcoefficient_i(swloop_m) * swcoefficient_j(swloop_n) * swcoefficient_i(swloop_o) * swcoefficient_j(swloop_p) * &
								        V_Integral_2e(AO_xyz_factor(L_i,M_i), AO_xyz_factor(L_j,M_j), AO_xyz_factor(L_i,M_i), AO_xyz_factor(L_j,M_j), &
									        swexponents_i(swloop_m),swexponents_j(swloop_n),swexponents_i(swloop_o),swexponents_j(swloop_p),coordinate_i,coordinate_j,coordinate_i,coordinate_j)
							        else
								        swintegral_mic = swcoefficient_i(swloop_m) * swcoefficient_j(swloop_n) * swcoefficient_i(swloop_o) * swcoefficient_j(swloop_p) * &
								        V_Integral_2e(AO_xyz_factor(L_j,M_j), AO_xyz_factor(L_i,M_i), AO_xyz_factor(L_j,M_j), AO_xyz_factor(L_i,M_i), &
									        swexponents_j(swloop_n),swexponents_i(swloop_m),swexponents_j(swloop_p),swexponents_i(swloop_o),coordinate_j,coordinate_i,coordinate_j,coordinate_i)
                                    end if
                                    swintegral(swloop_i,swloop_j) = swintegral(swloop_i,swloop_j) + swintegral_mic
                                end do
                            end do
                        end do
                    end do
                    if (isnan(swintegral(swloop_i,swloop_j))) call terminate('calculation of swintegral (Schwarz screening) failure, NaN detected.')
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
            write(60,'(a,e10.2,a)') '      complete! cutoff:',schwarz_VT,'; stored in swintegral'
            allocate(Fock2_assigned(basis_dimension,basis_dimension))
            allocate(Fock2_mic(2*basis_dimension,2*basis_dimension))
            allocate(oper(2*basis_dimension,2*basis_dimension))
        end if
        Fock2 = cmplx(0.0_dp,0.0_dp,dp)
        ! parallel zone, running results consistent with serial
        !$omp parallel num_threads(threads_use) default(shared) private(i,dloop_i,dloop_j,dloop_k,dloop_l,dloop_m,dloop_n,dloop_o,dloop_p,integral, &
        !$omp& contraction_i,atom_i,shell_i,L_i,M_i,contraction_j,atom_j,shell_j,L_j,M_j, &
        !$omp& contraction_k,atom_k,shell_k,L_k,M_k,contraction_l,atom_l,shell_l,L_l,M_l, &
        !$omp& exponents_i,exponents_j,exponents_k,exponents_l,coefficient_i,coefficient_j,coefficient_k,coefficient_l,coordinate_i,coordinate_j, &
        !$omp& coordinate_k,coordinate_l,Fock2_mic,Fock2_assigned) if(threads_use < cpu_threads)
        !$omp do schedule(dynamic,5)
        ! (ij|kl)=(ji|kl)=(ij|lk)=(ji|lk)=(kl|ij)=(kl|ji)=(lk|ij)=(lk|ji)
        !--------------------------------------------------------<dloop_i>--------------------------------------------------------
        do i = basis_dimension, 1, -1
            Fock2_mic = cmplx(0.0_dp,0.0_dp,dp)
            dloop_i = i
            contraction_i = atom_basis(molecular(basis_inf(dloop_i) % atom) % basis_number + basis_inf(dloop_i) % shell - 1) % contraction
            L_i = basis_inf(dloop_i) % L
            M_i = basis_inf(dloop_i) % M
            do dloop_m = 1, contraction_i
                exponents_i(dloop_m) = atom_basis(molecular(basis_inf(dloop_i) % atom) % basis_number + basis_inf(dloop_i) % shell - 1) % exponents(dloop_m)
                coefficient_i(dloop_m) = atom_basis(molecular(basis_inf(dloop_i) % atom) % basis_number + basis_inf(dloop_i) % shell - 1) % Ncoefficient(dloop_m,M_i)
            end do
            coordinate_i = molecular(basis_inf(dloop_i) % atom) % nucleus_position
        !--------------------------------------------------------<dloop_k>--------------------------------------------------------
            do dloop_k = basis_dimension, 1, -1
                contraction_k = atom_basis(molecular(basis_inf(dloop_k) % atom) % basis_number + basis_inf(dloop_k) % shell - 1) % contraction
                L_k = basis_inf(dloop_k) % L
                M_k = basis_inf(dloop_k) % M
                do dloop_m = 1, contraction_k
                    exponents_k(dloop_m) = atom_basis(molecular(basis_inf(dloop_k) % atom) % basis_number + basis_inf(dloop_k) % shell - 1) % exponents(dloop_m)
                    coefficient_k(dloop_m) = atom_basis(molecular(basis_inf(dloop_k) % atom) % basis_number + basis_inf(dloop_k) % shell - 1) % Ncoefficient(dloop_m,M_k)
                end do
                coordinate_k = molecular(basis_inf(dloop_k) % atom) % nucleus_position
        !--------------------------------------------------------<dloop_j>--------------------------------------------------------
                do dloop_j = dloop_i, 1, -1
                    contraction_j = atom_basis(molecular(basis_inf(dloop_j) % atom) % basis_number + basis_inf(dloop_j) % shell - 1) % contraction
                    L_j = basis_inf(dloop_j) % L
                    M_j = basis_inf(dloop_j) % M
                    do dloop_m = 1, contraction_j
                        exponents_j(dloop_m) = atom_basis(molecular(basis_inf(dloop_j) % atom) % basis_number + basis_inf(dloop_j) % shell - 1) % exponents(dloop_m)
                        coefficient_j(dloop_m) = atom_basis(molecular(basis_inf(dloop_j) % atom) % basis_number + basis_inf(dloop_j) % shell - 1) % Ncoefficient(dloop_m,M_j)
                    end do
                    coordinate_j = molecular(basis_inf(dloop_j) % atom) % nucleus_position
        !--------------------------------------------------------<dloop_l>--------------------------------------------------------
                    do dloop_l = min(dloop_k,dloop_j+(dloop_i*(dloop_i-1)-dloop_k*(dloop_k-1))/2), 1, -1
                        ! Schwarz screening of <ij||kl>, |<ij||kl>| <= sqrt(<ij||ij>) * sqrt(<kl||kl>)
                        if (sqrt(swintegral(dloop_i,dloop_j)*swintegral(dloop_k,dloop_l)) < schwarz_VT) cycle
                        contraction_l = atom_basis(molecular(basis_inf(dloop_l) % atom) % basis_number + basis_inf(dloop_l) % shell - 1) % contraction
                        L_l = basis_inf(dloop_l) % L
                        M_l = basis_inf(dloop_l) % M
                        do dloop_m = 1, contraction_l
                            exponents_l(dloop_m) = atom_basis(molecular(basis_inf(dloop_l) % atom) % basis_number + basis_inf(dloop_l) % shell - 1) % exponents(dloop_m)
                            coefficient_l(dloop_m) = atom_basis(molecular(basis_inf(dloop_l) % atom) % basis_number + basis_inf(dloop_l) % shell - 1) % Ncoefficient(dloop_m,M_l)
                        end do
                        coordinate_l = molecular(basis_inf(dloop_l) % atom) % nucleus_position
                        integral = 0.0_dp
        				!================================================<ij||kl>================================================
                        do dloop_m = 1, contraction_i
                            do dloop_n = 1, contraction_j
                                do dloop_o = 1, contraction_k
                                    do dloop_p = 1, contraction_l
	                                    if (L_i >= L_j .and. L_k >= L_l) then
	                                        integral = integral + coefficient_i(dloop_m) * coefficient_j(dloop_n) * coefficient_k(dloop_o) * coefficient_l(dloop_p) * &
	                                        	V_Integral_2e(AO_xyz_factor(L_i,M_i), AO_xyz_factor(L_j,M_j), AO_xyz_factor(L_k,M_k), AO_xyz_factor(L_l,M_l), &
	                                            exponents_i(dloop_m),exponents_j(dloop_n),exponents_k(dloop_o),exponents_l(dloop_p),coordinate_i,coordinate_j,coordinate_k,coordinate_l)
	                                    else if (L_i >= L_j .and. L_k < L_l) then
											integral = integral + coefficient_i(dloop_m) * coefficient_j(dloop_n) * coefficient_k(dloop_o) * coefficient_l(dloop_p) * &
                                            	V_Integral_2e(AO_xyz_factor(L_i,M_i), AO_xyz_factor(L_j,M_j), AO_xyz_factor(L_l,M_l), AO_xyz_factor(L_k,M_k), &
                                                exponents_i(dloop_m),exponents_j(dloop_n),exponents_l(dloop_p),exponents_k(dloop_o),coordinate_i,coordinate_j,coordinate_l,coordinate_k)
	                                    else if (L_i < L_j .and. L_k >= L_l) then
											integral = integral + coefficient_i(dloop_m) * coefficient_j(dloop_n) * coefficient_k(dloop_o) * coefficient_l(dloop_p) * &
                                            	V_Integral_2e(AO_xyz_factor(L_j,M_j), AO_xyz_factor(L_i,M_i), AO_xyz_factor(L_k,M_k), AO_xyz_factor(L_l,M_l), &
                                                exponents_j(dloop_n),exponents_i(dloop_m),exponents_k(dloop_o),exponents_l(dloop_p),coordinate_j,coordinate_i,coordinate_k,coordinate_l)
	                                    else
											integral = integral + coefficient_i(dloop_m) * coefficient_j(dloop_n) * coefficient_k(dloop_o) * coefficient_l(dloop_p) * &
                                            	V_Integral_2e(AO_xyz_factor(L_j,M_j), AO_xyz_factor(L_i,M_i), AO_xyz_factor(L_l,M_l), AO_xyz_factor(L_k,M_k), &
                                                exponents_j(dloop_n),exponents_i(dloop_m),exponents_l(dloop_p),exponents_k(dloop_o),coordinate_j,coordinate_i,coordinate_l,coordinate_k)
                                        end if
                                        !----------<DEBUG>----------
                                        !if (dloop_i == 79 .and. dloop_j == 34 .and. dloop_k == 52 .and. dloop_l == 34) then
                                        !    write(60,*) 'integral',integral
                                        !    write(60,*) AO_xyz_factor(L_i,M_i), AO_xyz_factor(L_j,M_j), AO_xyz_factor(L_k,M_k), AO_xyz_factor(L_l,M_l)
                                        !    write(60,*) exponents_i(dloop_m),exponents_j(dloop_n),exponents_k(dloop_o),exponents_l(dloop_p)
                                        !    write(60,*) coordinate_i,coordinate_j,coordinate_k,coordinate_l
                                        !    write(60,*)
                                        !    write(60,*) coefficient_i(dloop_m), coefficient_j(dloop_n), coefficient_k(dloop_o), coefficient_l(dloop_p)
                                        !end if
                                        !----------<DEBUG>----------
                                    end do
                                end do
                            end do
                        end do
                        !----------<DEBUG>----------
                        !if (inin) write(60,'(i3,i3,i3,i3,f)')  dloop_i, dloop_j, dloop_k, dloop_l, integral
                        !----------<DEBUG>----------
                        if (isnan(integral)) call terminate('calculation of <ij||kl> failure, NaN detected.')
                        ! assign values to two-electron Fock matrix
                        ! ref page 261 of Quantum Chemistry: Basic Principles and Ab Initio Calculations, Volume 2 | 2nd Edition
                        !------------------------<COULOMB INTEGRAL>------------------------
                        Fock2_assigned = .true.
                            Fock2_assigned(dloop_i,dloop_j) = .false.
                            Fock2_mic(2*dloop_i-1,2*dloop_j-1) = Fock2_mic(2*dloop_i-1,2*dloop_j-1) + integral*rou_m(2*dloop_k-1,2*dloop_l-1)       ! alpha->alpha Coulomb
                            Fock2_mic(2*dloop_i-1,2*dloop_j-1) = Fock2_mic(2*dloop_i-1,2*dloop_j-1) + integral*rou_m(2*dloop_k,2*dloop_l)           ! alpha->beta Coulomb
                            Fock2_mic(2*dloop_i,2*dloop_j) = Fock2_mic(2*dloop_i,2*dloop_j) + integral*rou_m(2*dloop_k-1,2*dloop_l-1)               ! beta->alpha Coulomb
                            Fock2_mic(2*dloop_i,2*dloop_j) = Fock2_mic(2*dloop_i,2*dloop_j) + integral*rou_m(2*dloop_k,2*dloop_l)                   ! beta->beta Coulomb
                            if (dloop_k /= dloop_l) then
                                Fock2_mic(2*dloop_i-1,2*dloop_j-1) = Fock2_mic(2*dloop_i-1,2*dloop_j-1) + integral*rou_m(2*dloop_l-1,2*dloop_k-1)   ! alpha->alpha Coulomb
                                Fock2_mic(2*dloop_i-1,2*dloop_j-1) = Fock2_mic(2*dloop_i-1,2*dloop_j-1) + integral*rou_m(2*dloop_l,2*dloop_k)       ! alpha->beta Coulomb
                                Fock2_mic(2*dloop_i,2*dloop_j) = Fock2_mic(2*dloop_i,2*dloop_j) + integral*rou_m(2*dloop_l-1,2*dloop_k-1)           ! beta->alpha Coulomb
                                Fock2_mic(2*dloop_i,2*dloop_j) = Fock2_mic(2*dloop_i,2*dloop_j) + integral*rou_m(2*dloop_l,2*dloop_k)               ! beta->beta Coulomb
                            end if
                        if (Fock2_assigned(dloop_j,dloop_i)) then
                            Fock2_assigned(dloop_j,dloop_i) = .false.
                            Fock2_mic(2*dloop_j-1,2*dloop_i-1) = Fock2_mic(2*dloop_j-1,2*dloop_i-1) + integral*rou_m(2*dloop_k-1,2*dloop_l-1)       ! alpha->alpha Coulomb
                            Fock2_mic(2*dloop_j-1,2*dloop_i-1) = Fock2_mic(2*dloop_j-1,2*dloop_i-1) + integral*rou_m(2*dloop_k,2*dloop_l)           ! alpha->beta Coulomb
                            Fock2_mic(2*dloop_j,2*dloop_i) = Fock2_mic(2*dloop_j,2*dloop_i) + integral*rou_m(2*dloop_k-1,2*dloop_l-1)               ! beta->alpha Coulomb
                            Fock2_mic(2*dloop_j,2*dloop_i) = Fock2_mic(2*dloop_j,2*dloop_i) + integral*rou_m(2*dloop_k,2*dloop_l)                   ! beta->beta Coulomb
                            if (dloop_k /= dloop_l) then
                                Fock2_mic(2*dloop_j-1,2*dloop_i-1) = Fock2_mic(2*dloop_j-1,2*dloop_i-1) + integral*rou_m(2*dloop_l-1,2*dloop_k-1)   ! alpha->alpha Coulomb
                                Fock2_mic(2*dloop_j-1,2*dloop_i-1) = Fock2_mic(2*dloop_j-1,2*dloop_i-1) + integral*rou_m(2*dloop_l,2*dloop_k)       ! alpha->beta Coulomb
                                Fock2_mic(2*dloop_j,2*dloop_i) = Fock2_mic(2*dloop_j,2*dloop_i) + integral*rou_m(2*dloop_l-1,2*dloop_k-1)           ! beta->alpha Coulomb
                                Fock2_mic(2*dloop_j,2*dloop_i) = Fock2_mic(2*dloop_j,2*dloop_i) + integral*rou_m(2*dloop_l,2*dloop_k)               ! beta->beta Coulomb
                            end if
                        end if
                        if (Fock2_assigned(dloop_k,dloop_l)) then
                            Fock2_assigned(dloop_k,dloop_l) = .false.
                            Fock2_mic(2*dloop_k-1,2*dloop_l-1) = Fock2_mic(2*dloop_k-1,2*dloop_l-1) + integral*rou_m(2*dloop_i-1,2*dloop_j-1)       ! alpha->alpha Coulomb
                            Fock2_mic(2*dloop_k-1,2*dloop_l-1) = Fock2_mic(2*dloop_k-1,2*dloop_l-1) + integral*rou_m(2*dloop_i,2*dloop_j)           ! alpha->beta Coulomb
                            Fock2_mic(2*dloop_k,2*dloop_l) = Fock2_mic(2*dloop_k,2*dloop_l) + integral*rou_m(2*dloop_i-1,2*dloop_j-1)               ! beta->alpha Coulomb
                            Fock2_mic(2*dloop_k,2*dloop_l) = Fock2_mic(2*dloop_k,2*dloop_l) + integral*rou_m(2*dloop_i,2*dloop_j)                   ! beta->beta Coulomb
                            if (dloop_i /= dloop_j) then
                                Fock2_mic(2*dloop_k-1,2*dloop_l-1) = Fock2_mic(2*dloop_k-1,2*dloop_l-1) + integral*rou_m(2*dloop_j-1,2*dloop_i-1)   ! alpha->alpha Coulomb
                                Fock2_mic(2*dloop_k-1,2*dloop_l-1) = Fock2_mic(2*dloop_k-1,2*dloop_l-1) + integral*rou_m(2*dloop_j,2*dloop_i)       ! alpha->beta Coulomb
                                Fock2_mic(2*dloop_k,2*dloop_l) = Fock2_mic(2*dloop_k,2*dloop_l) + integral*rou_m(2*dloop_j-1,2*dloop_i-1)           ! beta->alpha Coulomb
                                Fock2_mic(2*dloop_k,2*dloop_l) = Fock2_mic(2*dloop_k,2*dloop_l) + integral*rou_m(2*dloop_j,2*dloop_i)               ! beta->beta Coulomb
                            end if
                        end if
                        if (Fock2_assigned(dloop_l,dloop_k)) then
                            Fock2_assigned(dloop_l,dloop_k) = .false.
                            Fock2_mic(2*dloop_l-1,2*dloop_k-1) = Fock2_mic(2*dloop_l-1,2*dloop_k-1) + integral*rou_m(2*dloop_i-1,2*dloop_j-1)       ! alpha->alpha Coulomb
                            Fock2_mic(2*dloop_l-1,2*dloop_k-1) = Fock2_mic(2*dloop_l-1,2*dloop_k-1) + integral*rou_m(2*dloop_i,2*dloop_j)           ! alpha->beta Coulomb
                            Fock2_mic(2*dloop_l,2*dloop_k) = Fock2_mic(2*dloop_l,2*dloop_k) + integral*rou_m(2*dloop_i-1,2*dloop_j-1)               ! beta->alpha Coulomb
                            Fock2_mic(2*dloop_l,2*dloop_k) = Fock2_mic(2*dloop_l,2*dloop_k) + integral*rou_m(2*dloop_i,2*dloop_j)                   ! beta->beta Coulomb
                            if (dloop_i /= dloop_j) then
                                Fock2_mic(2*dloop_l-1,2*dloop_k-1) = Fock2_mic(2*dloop_l-1,2*dloop_k-1) + integral*rou_m(2*dloop_j-1,2*dloop_i-1)   ! alpha->alpha Coulomb
                                Fock2_mic(2*dloop_l-1,2*dloop_k-1) = Fock2_mic(2*dloop_l-1,2*dloop_k-1) + integral*rou_m(2*dloop_j,2*dloop_i)       ! alpha->beta Coulomb
                                Fock2_mic(2*dloop_l,2*dloop_k) = Fock2_mic(2*dloop_l,2*dloop_k) + integral*rou_m(2*dloop_j-1,2*dloop_i-1)           ! beta->alpha Coulomb
                                Fock2_mic(2*dloop_l,2*dloop_k) = Fock2_mic(2*dloop_l,2*dloop_k) + integral*rou_m(2*dloop_j,2*dloop_i)               ! beta->beta Coulomb
                            end if
                        end if
                        !------------------------<EXCHANGE INTEGRAL>------------------------
                        Fock2_assigned = .true.
                            Fock2_assigned(dloop_i,dloop_k) = .false.
                            Fock2_mic(2*dloop_i-1,2*dloop_k-1) = Fock2_mic(2*dloop_i-1,2*dloop_k-1) - integral*rou_m(2*dloop_j-1,2*dloop_l-1)       ! alpha->alpha Exchange
                            Fock2_mic(2*dloop_i,2*dloop_k) = Fock2_mic(2*dloop_i,2*dloop_k) - integral*rou_m(2*dloop_j,2*dloop_l)                   ! beta->beta Exchange
                            if (dloop_i == dloop_k .and. dloop_l /= dloop_j) then
                                Fock2_mic(2*dloop_i-1,2*dloop_k-1) = Fock2_mic(2*dloop_i-1,2*dloop_k-1) - integral*rou_m(2*dloop_l-1,2*dloop_j-1)   ! alpha->alpha Exchange
                                Fock2_mic(2*dloop_i,2*dloop_k) = Fock2_mic(2*dloop_i,2*dloop_k) - integral*rou_m(2*dloop_l,2*dloop_j)               ! beta->beta Exchange
                            end if
                        if (Fock2_assigned(dloop_k,dloop_i)) then
                            Fock2_assigned(dloop_k,dloop_i) = .false.
                            Fock2_mic(2*dloop_k-1,2*dloop_i-1) = Fock2_mic(2*dloop_k-1,2*dloop_i-1) - integral*rou_m(2*dloop_l-1,2*dloop_j-1)       ! alpha->alpha Exchange
                            Fock2_mic(2*dloop_k,2*dloop_i) = Fock2_mic(2*dloop_k,2*dloop_i) - integral*rou_m(2*dloop_l,2*dloop_j)                   ! beta->beta Exchange
                            if (dloop_k == dloop_i .and. dloop_l /= dloop_j) then
                                Fock2_mic(2*dloop_k-1,2*dloop_i-1) = Fock2_mic(2*dloop_k-1,2*dloop_i-1) - integral*rou_m(2*dloop_j-1,2*dloop_l-1)   ! alpha->alpha Exchange
                                Fock2_mic(2*dloop_k,2*dloop_i) = Fock2_mic(2*dloop_k,2*dloop_i) - integral*rou_m(2*dloop_j,2*dloop_l)               ! beta->beta Exchange
                            end if
                        end if
                        if (Fock2_assigned(dloop_i,dloop_l)) then
                            Fock2_assigned(dloop_i,dloop_l) = .false.
                            Fock2_mic(2*dloop_i-1,2*dloop_l-1) = Fock2_mic(2*dloop_i-1,2*dloop_l-1) - integral*rou_m(2*dloop_j-1,2*dloop_k-1)       ! alpha->alpha Exchange
                            Fock2_mic(2*dloop_i,2*dloop_l) = Fock2_mic(2*dloop_i,2*dloop_l) - integral*rou_m(2*dloop_j,2*dloop_k)                   ! beta->beta Exchange
                            if (dloop_i == dloop_l .and. dloop_k /= dloop_j) then
                                Fock2_mic(2*dloop_i-1,2*dloop_l-1) = Fock2_mic(2*dloop_i-1,2*dloop_l-1) - integral*rou_m(2*dloop_k-1,2*dloop_j-1)   ! alpha->alpha Exchange
                                Fock2_mic(2*dloop_i,2*dloop_l) = Fock2_mic(2*dloop_i,2*dloop_l) - integral*rou_m(2*dloop_k,2*dloop_j)               ! beta->beta Exchange
                            end if
                        end if
                        if (Fock2_assigned(dloop_l,dloop_i)) then
                            Fock2_assigned(dloop_l,dloop_i) = .false.
                            Fock2_mic(2*dloop_l-1,2*dloop_i-1) = Fock2_mic(2*dloop_l-1,2*dloop_i-1) - integral*rou_m(2*dloop_k-1,2*dloop_j-1)       ! alpha->alpha Exchange
                            Fock2_mic(2*dloop_l,2*dloop_i) = Fock2_mic(2*dloop_l,2*dloop_i) - integral*rou_m(2*dloop_k,2*dloop_j)                   ! beta->beta Exchange
                            if (dloop_l == dloop_i .and. dloop_k /= dloop_j) then
                                Fock2_mic(2*dloop_l-1,2*dloop_i-1) = Fock2_mic(2*dloop_l-1,2*dloop_i-1) - integral*rou_m(2*dloop_j-1,2*dloop_k-1)   ! alpha->alpha Exchange
                                Fock2_mic(2*dloop_l,2*dloop_i) = Fock2_mic(2*dloop_l,2*dloop_i) - integral*rou_m(2*dloop_j,2*dloop_k)               ! beta->beta Exchange
                            end if
                        end if
                        if (Fock2_assigned(dloop_j,dloop_k)) then
                            Fock2_assigned(dloop_j,dloop_k) = .false.
                            Fock2_mic(2*dloop_j-1,2*dloop_k-1) = Fock2_mic(2*dloop_j-1,2*dloop_k-1) - integral*rou_m(2*dloop_i-1,2*dloop_l-1)       ! alpha->alpha Exchange
                            Fock2_mic(2*dloop_j,2*dloop_k) = Fock2_mic(2*dloop_j,2*dloop_k) - integral*rou_m(2*dloop_i,2*dloop_l)                   ! beta->beta Exchange
                            if (dloop_j == dloop_k .and. dloop_l /= dloop_i) then
                                Fock2_mic(2*dloop_j-1,2*dloop_k-1) = Fock2_mic(2*dloop_j-1,2*dloop_k-1) - integral*rou_m(2*dloop_l-1,2*dloop_i-1)   ! alpha->alpha Exchange
                                Fock2_mic(2*dloop_j,2*dloop_k) = Fock2_mic(2*dloop_j,2*dloop_k) - integral*rou_m(2*dloop_l,2*dloop_i)               ! beta->beta Exchange
                            end if
                        end if
                        if (Fock2_assigned(dloop_k,dloop_j)) then
                            Fock2_assigned(dloop_k,dloop_j) = .false.
                            Fock2_mic(2*dloop_k-1,2*dloop_j-1) = Fock2_mic(2*dloop_k-1,2*dloop_j-1) - integral*rou_m(2*dloop_l-1,2*dloop_i-1)       ! alpha->alpha Exchange
                            Fock2_mic(2*dloop_k,2*dloop_j) = Fock2_mic(2*dloop_k,2*dloop_j) - integral*rou_m(2*dloop_l,2*dloop_i)                   ! beta->beta Exchange
                            if (dloop_k == dloop_j .and. dloop_i /= dloop_l) then
                                Fock2_mic(2*dloop_k-1,2*dloop_j-1) = Fock2_mic(2*dloop_k-1,2*dloop_j-1) - integral*rou_m(2*dloop_i-1,2*dloop_l-1)   ! alpha->alpha Exchange
                                Fock2_mic(2*dloop_k,2*dloop_j) = Fock2_mic(2*dloop_k,2*dloop_j) - integral*rou_m(2*dloop_i,2*dloop_l)               ! beta->beta Exchange
                            end if
                        end if
                        if (Fock2_assigned(dloop_j,dloop_l)) then
                            Fock2_assigned(dloop_j,dloop_l) = .false.
                            Fock2_mic(2*dloop_j-1,2*dloop_l-1) = Fock2_mic(2*dloop_j-1,2*dloop_l-1) - integral*rou_m(2*dloop_i-1,2*dloop_k-1)       ! alpha->alpha Exchange
                            Fock2_mic(2*dloop_j,2*dloop_l) = Fock2_mic(2*dloop_j,2*dloop_l) - integral*rou_m(2*dloop_i,2*dloop_k)                   ! beta->beta Exchange
                            if (dloop_j == dloop_l .and. dloop_k /= dloop_i) then
                                Fock2_mic(2*dloop_j-1,2*dloop_l-1) = Fock2_mic(2*dloop_j-1,2*dloop_l-1) - integral*rou_m(2*dloop_k-1,2*dloop_i-1)   ! alpha->alpha Exchange
                                Fock2_mic(2*dloop_j,2*dloop_l) = Fock2_mic(2*dloop_j,2*dloop_l) - integral*rou_m(2*dloop_k,2*dloop_i)               ! beta->beta Exchange
                            end if
                        end if
                        if (Fock2_assigned(dloop_l,dloop_j)) then
                            Fock2_assigned(dloop_l,dloop_j) = .false.
                            Fock2_mic(2*dloop_l-1,2*dloop_j-1) = Fock2_mic(2*dloop_l-1,2*dloop_j-1) - integral*rou_m(2*dloop_k-1,2*dloop_i-1)       ! alpha->alpha Exchange
                            Fock2_mic(2*dloop_l,2*dloop_j) = Fock2_mic(2*dloop_l,2*dloop_j) - integral*rou_m(2*dloop_k,2*dloop_i)                   ! beta->beta Exchange
                            if (dloop_l == dloop_j .and. dloop_i /= dloop_k) then
                                Fock2_mic(2*dloop_l-1,2*dloop_j-1) = Fock2_mic(2*dloop_l-1,2*dloop_j-1) - integral*rou_m(2*dloop_i-1,2*dloop_k-1)   ! alpha->alpha Exchange
                                Fock2_mic(2*dloop_l,2*dloop_j) = Fock2_mic(2*dloop_l,2*dloop_j) - integral*rou_m(2*dloop_i,2*dloop_k)               ! beta->beta Exchange
                            end if
                        end if
                    end do
                end do
            end do
            !$omp critical
            Fock2 = Fock2 + Fock2_mic
            !$omp end critical
        end do
        !$omp end do
        !$omp end parallel
        ! Lowdin orth to Fock2
        call zgemm( 'T', 'N', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), exS0_5, 2*basis_dimension, Fock2, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), oper, 2*basis_dimension)
        call zgemm( 'N', 'N', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), oper, 2*basis_dimension, exS0_5, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), Fock2, 2*basis_dimension)
    end subroutine Fock2e

    !------------------------------------------------------------
	! transfer alternating zero matrix to block matrix
    subroutine atnz2block(atnz, dm)
    	integer, intent(in) :: dm
    	integer :: aloop_i, aloop_j
        complex(dp) :: atnz(dm,dm), aoper(dm,dm)
    	if (dm <= 1 .or. mod(dm,2) /= 0) call terminate('atnz2block called incorrectly')
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
end module SCF