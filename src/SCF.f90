! |-------<THOMAS RELATIVISTIC ELECTRONIC STRUCTURE CALCULATIONS>------|
! |------------------------ -----<DIRAC4PI>----------------------------|
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
    complex(dp),allocatable :: AO2MO(:,:)                            ! linear combination coefficients of MOs, basis spinor represented by (AO1,0), (0,AO1), (AO2,0), (0,AO2) ... (AOn,0), (0,AOn)
    real(dp),allocatable :: Gaualpha(:)                              ! Gaussian alpha orbital coefficient
    real(dp),allocatable :: AO2MOalpha(:,:)                          ! Gaussian alpha orbital coefficient
    real(dp),allocatable :: Gaubeta(:)                               ! Gaussian beta orbital coefficient
    real(dp),allocatable :: AO2MObeta(:,:)                           ! Gaussian beta orbital coefficient
    complex(dp),allocatable :: rou_m(:,:)                            ! density matrix, complex Hermitian matrix
    complex(dp),allocatable :: rotation(:,:)                         ! rotate one orbital with another
    real(dp) :: RMSDP, maxDP
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
    logical :: diiscuted = .false.                                   ! whether cutdiis ever reached
    
    
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
    
    real(dp),allocatable :: swintegral(:,:)                          ! <ij||ij> as well as <kl||kl>
    logical :: ndschwarz = .true.
    logical,allocatable :: Fock2_assigned(:,:)                       ! avoid duplicate assignment of Fock2
    complex(dp),allocatable :: oper(:,:)
    
    ! orbital energy and molecular energy
    real(dp),allocatable :: orbE(:)                                  ! orbital energy
    real(dp) :: molE_pre, molE                                       ! molecular energy
    real(dp) :: nucE                                                 ! nuclear repulsion energy
    real(dp) :: eleE                                                 ! electron repulsion energy
    real(dp) :: T                                                    ! kinetic energy
    real(dp) :: V                                                    ! electron-nuclear attraction energy
    real(dp) :: ESOC                                                 ! SOC energy
    real(dp) :: ESR                                                  ! second relativistic Thomas precession energy
    real(dp) :: emd4                                                 ! dispersion energy calculated by DFT-D4

    ! DIIS Ax=B
    complex(dp),allocatable :: rou_pre(:,:)                          ! privious rou_m of current iteration
    complex(dp),allocatable :: rou_history(:,:,:)                    ! coefficients of subsp iteration
    complex(dp),allocatable :: Rsd(:,:,:)                            ! residuals of rou_m
    complex(dp),allocatable :: DIISmat(:,:)                          ! A
    complex(dp),allocatable :: DIISwork(:)
    integer :: lDIISwork
    integer,allocatable :: ipiv(:)
    integer :: DIIsinfo
    
    
    !--------------<DEBUG>--------------
    real(dp),allocatable :: fockerr(:)
    real(dp),allocatable :: focksum(:)
    real(dp) :: gaufock(5)
    integer :: rdd
        
    
    contains

!------------------------------------------------------------
! standard Hartree Fock SCF procedure based on the Pulay (DIIS) mixing method for DKH0 & DKH2 Hamiltonian 
    subroutine DKH_Hartree_Fock(keep, kill)
        integer,intent(in) :: keep                                    ! initialization level
        logical,intent(in) :: kill                                    ! whether to kill the whole program after SCF process
        integer :: iatom, ishell
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
        write(60,'(a45,f12.7)') '   complete! nuclear repulsive energy: (A.U.)', nucE
        write(60,'(a,i4,a,e15.7)') '   SCF maxiter:',maxiter,';  convergence tolerance:',conver_tol
        write(60,'(a,i3,a,i3,a,f12.6,a,e15.7)') '   nodiis =',nodiis,';  subsp =',subsp,';  damp =',damp,';  cutdiis =',cutdiis
        write(60,'(a)') '   -------------------------<SCF>-------------------------'
        allocate(Fock(2*basis_dimension,2*basis_dimension))
        allocate(orbE(2*basis_dimension))
        allocate(oper6(2*basis_dimension,2*basis_dimension))
        allocate(isupp_ev_f(4*basis_dimension))
        allocate(oper3(2*basis_dimension,2*basis_dimension))
        allocate(oper4(2*basis_dimension,2*basis_dimension))
        

        ! DIIS
        ! AO2MO(new) = ¦²(i,subsp) DIIScoe(i)*(rou_history(i)+damp*Rsd(i))
        allocate(Rsd(subsp,2*basis_dimension,2*basis_dimension))
        allocate(DIISmat(subsp+1,subsp+1))
        allocate(ipiv(subsp+1))
        allocate(rou_history(subsp, 2*basis_dimension, 2*basis_dimension))
        allocate(rou_pre(2*basis_dimension, 2*basis_dimension))
        do loop_i = 1, maxiter
            if (loop_i /= 1) write(60,'(a)') '   -------------------------------------------------------------'
            write(60,*)
            write(60,*)
            write(60,'(a,i3)') '   SCF iter ',loop_i
            if (loop_i == 1) then
                write(60,'(a)') '   read density matrix from Gaussian checkpoint file'
                call assign_rou()
                write(60,'(a)') '   complete! stored in rou_m'
            end if
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
            
            
            !---------------------<DEBUG>---------------------
            ! print Fock and check with Gaussian Fock
            !call zgemm( 'T', 'N', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), exSp0_5, 2*basis_dimension, Fock1, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), oper3, 2*basis_dimension)
            !call zgemm( 'N', 'N', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), oper3, 2*basis_dimension, exSp0_5, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), Fock1, 2*basis_dimension)
            !call zgemm( 'T', 'N', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), exSp0_5, 2*basis_dimension, Fock2, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), oper3, 2*basis_dimension)
            !call zgemm( 'N', 'N', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), oper3, 2*basis_dimension, exSp0_5, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), Fock2, 2*basis_dimension)
            !allocate(fockerr(basis_dimension))
            !allocate(focksum(basis_dimension))
            !fockerr = 0.0_dp
            !focksum = 0.0_dp
            !write(60,*) 'Fock2(alpha)'
            !open(111,file='E:\TRESC\Fe(CN)6\Fock2alpha.out',action='read')
            !outeralpha: do loop_j=1, basis_dimension, 5
            !	read(111,*)
            !	do loop_k=loop_j, basis_dimension
            !		if (loop_j+loop_k > 2*basis_dimension) exit outeralpha
            !		if (loop_k == loop_j) then
            !			read(111,*) rdd, gaufock(1)
            !		else if (loop_k == loop_j + 1) then
            !			read(111,*) rdd, gaufock(1), gaufock(2)
            !		else if (loop_k == loop_j + 2) then
            !			read(111,*) rdd, gaufock(1), gaufock(2), gaufock(3)
            !		else if (loop_k == loop_j + 3) then
            !			read(111,*) rdd, gaufock(1), gaufock(2), gaufock(3), gaufock(4)
            !		else
            !			read(111,*) rdd, gaufock
            !		end if
            !		do loop_m=1, min(loop_k-loop_j+1,5)
            !			if (abs(real(Fock2(2*loop_k-1,2*(loop_j+loop_m-1)-1))-gaufock(loop_m))&
            !				/abs(real(Fock2(2*loop_k-1,2*(loop_j+loop_m-1)-1))) >= 0.01 .and. &
            !				abs(real(Fock2(2*loop_k-1,2*(loop_j+loop_m-1)-1))) >= 1E-7) then
            !				if (abs(real(Fock2(2*loop_k-1,2*(loop_j+loop_m-1)-1))-gaufock(loop_m)) >= 0.1) then
            !					write(60,'(i3,a,i1,a,i3,a,i1,a,f12.7,a,f12.7)') loop_k, ' L=',basis_inf(loop_k)%L-1,'  ',loop_j+loop_m-1, ' L=',basis_inf(loop_j+loop_m-1)%L-1,&
            !						'  ',gaufock(loop_m), '             ',real(Fock2(2*loop_k-1,2*(loop_j+loop_m-1)-1))-gaufock(loop_m)
            !				else if (abs(real(Fock2(2*loop_k-1,2*(loop_j+loop_m-1)-1))-gaufock(loop_m)) >= 0.01) then
            !					write(60,'(i3,a,i1,a,i3,a,i1,a,f12.7,a,f12.7)') loop_k, ' L=',basis_inf(loop_k)%L-1,'  ',loop_j+loop_m-1, ' L=',basis_inf(loop_j+loop_m-1)%L-1,&
            !						'  ',gaufock(loop_m), '         ',real(Fock2(2*loop_k-1,2*(loop_j+loop_m-1)-1))-gaufock(loop_m)
            !				else if (abs(real(Fock2(2*loop_k-1,2*(loop_j+loop_m-1)-1))-gaufock(loop_m)) >= 0.001) then
            !					write(60,'(i3,a,i1,a,i3,a,i1,a,f12.7,a,f12.7)') loop_k, ' L=',basis_inf(loop_k)%L-1,'  ',loop_j+loop_m-1, ' L=',basis_inf(loop_j+loop_m-1)%L-1,&
            !						'  ',gaufock(loop_m), '     ',real(Fock2(2*loop_k-1,2*(loop_j+loop_m-1)-1))-gaufock(loop_m)
            !				else
            !					write(60,'(i3,a,i1,a,i3,a,i1,a,f12.7,a,f12.7)') loop_k, ' L=',basis_inf(loop_k)%L-1,'  ',loop_j+loop_m-1, ' L=',basis_inf(loop_j+loop_m-1)%L-1,&
            !						'  ',gaufock(loop_m), ' ',real(Fock2(2*loop_k-1,2*(loop_j+loop_m-1)-1))-gaufock(loop_m)
            !				end if
            !			end if
            !			Fock(2*loop_k-1,2*(loop_j+loop_m-1)-1) = cmplx(gaufock(loop_m),0.0,dp)
            !			Fock(2*(loop_j+loop_m-1)-1,2*loop_k-1) = cmplx(gaufock(loop_m),0.0,dp)
            !			fockerr(loop_k) = fockerr(loop_k) + real(Fock2(2*loop_k-1,2*(loop_j+loop_m-1)-1))-gaufock(loop_m)
            !			focksum(loop_k) = focksum(loop_k) + gaufock(loop_m)
            !		end do
            !	end do
            !end do outeralpha
            !write(60,*)
            !write(60,*)
            !write(60,*) '        Focksum           Fockerr'
            !do loop_j=1, basis_dimension
            !	write(60,'(i3,f12.7,a,f12.7)') loop_j, Focksum(loop_j), '   ', Fockerr(loop_j)
            !end do
            !close(111)
            !fockerr = 0.0_dp
            !focksum = 0.0_dp
            !write(60,*) 'Fock2(beta)'
            !open(111,file='E:\TRESC\Fe(CN)6\Fock2beta.out',action='read')
            !outerbeta: do loop_j=1, basis_dimension, 5
            !	read(111,*)
            !	do loop_k=loop_j, basis_dimension
            !		if (loop_j+loop_k > 2*basis_dimension) exit outerbeta
            !		if (loop_k == loop_j) then
            !			read(111,*) rdd, gaufock(1)
            !		else if (loop_k == loop_j + 1) then
            !			read(111,*) rdd, gaufock(1), gaufock(2)
            !		else if (loop_k == loop_j + 2) then
            !			read(111,*) rdd, gaufock(1), gaufock(2), gaufock(3)
            !		else if (loop_k == loop_j + 3) then
            !			read(111,*) rdd, gaufock(1), gaufock(2), gaufock(3), gaufock(4)
            !		else
            !			read(111,*) rdd, gaufock
            !		end if
            !		do loop_m=1, min(loop_k-loop_j+1,5)
            !			if (abs(real(Fock2(2*loop_k,2*(loop_j+loop_m-1)))-gaufock(loop_m))&
            !				/abs(real(Fock2(2*loop_k,2*(loop_j+loop_m-1)))) >= 0.01 .and. &
            !				abs(real(Fock2(2*loop_k,2*(loop_j+loop_m-1)))) >= 1E-7) then
            !				if (abs(real(Fock2(2*loop_k,2*(loop_j+loop_m-1)))-gaufock(loop_m)) >= 0.1) then
            !					write(60,'(i3,a,i1,a,i3,a,i1,a,f12.7,a,f12.7)') loop_k, ' L=',basis_inf(loop_k)%L-1,'  ',loop_j+loop_m-1, ' L=',basis_inf(loop_j+loop_m-1)%L-1,&
            !						'  ',gaufock(loop_m), '             ',real(Fock2(2*loop_k,2*(loop_j+loop_m-1)))-gaufock(loop_m)
            !				else if (abs(real(Fock2(2*loop_k,2*(loop_j+loop_m-1)))-gaufock(loop_m)) >= 0.01) then
            !					write(60,'(i3,a,i1,a,i3,a,i1,a,f12.7,a,f12.7)') loop_k, ' L=',basis_inf(loop_k)%L-1,'  ',loop_j+loop_m-1, ' L=',basis_inf(loop_j+loop_m-1)%L-1,&
            !						'  ',gaufock(loop_m), '         ',real(Fock2(2*loop_k,2*(loop_j+loop_m-1)))-gaufock(loop_m)
            !				else if (abs(real(Fock2(2*loop_k,2*(loop_j+loop_m-1)))-gaufock(loop_m)) >= 0.001) then
            !					write(60,'(i3,a,i1,a,i3,a,i1,a,f12.7,a,f12.7)') loop_k, ' L=',basis_inf(loop_k)%L-1,'  ',loop_j+loop_m-1, ' L=',basis_inf(loop_j+loop_m-1)%L-1,&
            !						'  ',gaufock(loop_m), '     ',real(Fock2(2*loop_k,2*(loop_j+loop_m-1)))-gaufock(loop_m)
            !				else
            !					write(60,'(i3,a,i1,a,i3,a,i1,a,f12.7,a,f12.7)') loop_k, ' L=',basis_inf(loop_k)%L-1,'  ',loop_j+loop_m-1, ' L=',basis_inf(loop_j+loop_m-1)%L-1,&
            !						'  ',gaufock(loop_m), ' ',real(Fock2(2*loop_k,2*(loop_j+loop_m-1)))-gaufock(loop_m)
            !				end if
            !			end if
            !			Fock(2*loop_k,2*(loop_j+loop_m-1)) = cmplx(gaufock(loop_m),0.0,dp)
            !			Fock(2*(loop_j+loop_m-1),2*loop_k) = cmplx(gaufock(loop_m),0.0,dp)
            !			fockerr(loop_k) = fockerr(loop_k) + real(Fock2(2*loop_k,2*(loop_j+loop_m-1)))-gaufock(loop_m)
            !			focksum(loop_k) = focksum(loop_k) + gaufock(loop_m)
            !		end do
            !	end do
            !end do outerbeta
            !write(60,*)
            !write(60,*)
            !write(60,*) '   Focksum           Fockerr'
            !do loop_j=1, basis_dimension
            !	if (abs(Fockerr(loop_j)/Focksum(loop_j)) >= 0.1) then
            !		write(60,'(i3,f12.7,a,f12.7,a,i1,a,i2)') loop_j, Focksum(loop_j), '             ', Fockerr(loop_j), '  L=', basis_inf(loop_j)%L-1, '  A=', basis_inf(loop_j)%atom
            !	else if (abs(Fockerr(loop_j)/Focksum(loop_j)) >= 0.01) then
            !		write(60,'(i3,f12.7,a,f12.7,a,i1,a,i2)') loop_j, Focksum(loop_j), '         ', Fockerr(loop_j), '  L=', basis_inf(loop_j)%L-1, '  A=', basis_inf(loop_j)%atom
            !	else if (abs(Fockerr(loop_j)/Focksum(loop_j)) >= 0.001) then
            !		write(60,'(i3,f12.7,a,f12.7,a,i1,a,i2)') loop_j, Focksum(loop_j), '     ', Fockerr(loop_j), '  L=', basis_inf(loop_j)%L-1, '  A=', basis_inf(loop_j)%atom
            !	else
            !		write(60,'(i3,f12.7,a,f12.7,a,i1,a,i2)') loop_j, Focksum(loop_j), ' ', Fockerr(loop_j), '  L=', basis_inf(loop_j)%L-1, '  A=', basis_inf(loop_j)%atom
            !	end if
            !end do
            !close(111)
            !deallocate(fockerr)
            !deallocate(focksum)
            !call zgemm( 'T', 'N', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), exS0_5, 2*basis_dimension, Fock, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), oper3, 2*basis_dimension)
            !call zgemm( 'N', 'N', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), oper3, 2*basis_dimension, exS0_5, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), Fock, 2*basis_dimension)
            !do loop_j=1, 2*basis_dimension
            !	do loop_k=1, loop_j-1
            !		Fock(loop_j,loop_k) = cmplx(0.0,0.0,dp)
            !	end do
            !end do
            !---------------------<DEBUG>---------------------
            
            
            allocate(fwork(1))
            allocate(fiwork(1))
            allocate(rwork(1))
            call zheevr('V','A','U',2*basis_dimension,Fock,2*basis_dimension,0.0_dp,0.0_dp,0,0,dlamch('S'),&
                evl_count_f,orbE,oper3,2*basis_dimension,isupp_ev_f,fwork,-1,rwork,-1,fiwork,-1,finfo)
            flwork = nint(real(fwork(1)))
            fliwork = fiwork(1)
            lrwork = nint(rwork(1))
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
            ! frontier orbital energy
            write(60,'(a)') '   frontier orbital energy (A.U.)'
            call calc_S2HForb(electron_count)
            write(60,'(a,i3,f12.6,a,f6.3)') '   --- HOMO ', electron_count, orbE(electron_count), ' <Sz> = ',Szorb
            call calc_S2HForb(electron_count+1)
            write(60,'(a,i3,f12.6,a,f6.3)') '   --- LUMO ', electron_count+1, orbE(electron_count+1), ' <Sz> = ',Szorb
            write(60,'(a,f12.6)') '   --- gap  ', orbE(electron_count+1) - orbE(electron_count)
            deallocate(fwork)
            deallocate(fiwork)
            deallocate(rwork)
            ! energy components calculation
            write(60,'(a)') '   calculate energy components (A.U.)'
            eleE = 0.0_dp
            call zgemm( 'C', 'N', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), oper3, 2*basis_dimension, Fock2, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), oper6, 2*basis_dimension)
            call zgemm( 'N', 'N', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), oper6, 2*basis_dimension, oper3, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), oper4, 2*basis_dimension)
            do loop_j = 1, electron_count
                eleE = eleE + real(oper4(loop_j,loop_j))
            end do
            write(60,'(a,f12.6)') '   --- electron repulsive energy', eleE
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
            write(60,'(a,f12.6)') '   --- electron kinetic energy', T
            V = 0.0_dp
            call zgemm( 'C', 'N', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), oper3, 2*basis_dimension, exi_V_j, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), oper6, 2*basis_dimension)
            call zgemm( 'N', 'N', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), oper6, 2*basis_dimension, oper3, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), oper4, 2*basis_dimension)
            do loop_j = 1, electron_count
                V = V + real(oper4(loop_j,loop_j))
            end do
            write(60,'(a,f12.6)') '   --- electron-nuclear attraction energy', V
            if (DKH_order == 2) then
                ESOC = 0.0_dp
                call zgemm( 'C', 'N', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), oper3, 2*basis_dimension, exSOC, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), oper6, 2*basis_dimension)
                call zgemm( 'N', 'N', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), oper6, 2*basis_dimension, oper3, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), oper4, 2*basis_dimension)
                do loop_j = 1, electron_count
                    ESOC = ESOC + real(oper4(loop_j,loop_j))
                end do
                write(60,'(a,f12.6)') '   --- spin-orbital coupling energy', ESOC
                if (SRTP_type) then
                    ESR = 0.0_dp
                    call zgemm( 'C', 'N', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), oper3, 2*basis_dimension, exSR, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), oper6, 2*basis_dimension)
                    call zgemm( 'N', 'N', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), oper6, 2*basis_dimension, oper3, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), oper4, 2*basis_dimension)
                    do loop_j = 1, electron_count
                        ESR = ESR + real(oper4(loop_j,loop_j))
                    end do
                    write(60,'(a,f12.6)') '   --- second relativistic Thomas precession energy', ESR
                end if
            end if
            if (DKH_order == 0) then
                write(60,'(a,f12.6)') '   --- -<V/T> ', -(molE-T)/T
            else if (DKH_order == 2) then
                write(60,'(a,f12.6)') '   --- -<V/T> ', -(molE-T-ESOC)/T
            end if
            ! de-Lowdin orthogonalization
            call zgemm( 'T', 'N', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), oper3, 2*basis_dimension, exS0_5, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), AO2MO, 2*basis_dimension)
            write(60,'(a)') '   MO coefficient dump to .ao2mo file'
            call dump_matrix_cmplx(name='ao2mo', m=AO2MO, dm=2*basis_dimension)
            ! convergence check
            if (loop_i == 1) then
                write(60,'(a,f12.6)') '   SCF energy (A.U.) = ',molE
            else 
                if (abs(molE - molE_pre) < conver_tol) then
                    write(60,'(a,f12.6,a,e15.7)') '   SCF energy (A.U.) = ',molE,'; delta E = ',molE - molE_pre
                    write(60,'(a)') '   convergence tolerance met, SCF done!'
                    write(60,'(a)') '   -------------------------------------------------------------'
                    exit
                end if
                write(60,'(a,f12.6,a,e15.7)') '   SCF energy (A.U.) = ',molE,'; delta E = ',molE - molE_pre
                write(60,'(a)') '   convergence tolerance not met'
            end if
            write(60,'(a)') '   construct density matrix'
            call assign_rou()
            write(60,'(a)') '   complete! stored in rou_m'
            write(60,'(a)') '   DIIS information'
            if (abs(molE - molE_pre) < cutdiis .or. diiscuted) then
                diiscuted = .true.
                write(60,*) '   --- DIIS cutted'
                cycle
            end if
            ! generate next AO2MO by DIIS method
            if (loop_i <= nodiis - subsp) then
                write(60,'(a)') '   --- no DIIS acceleration'
            else if (nodiis - subsp < loop_i .and. loop_i <= nodiis) then
                ! update Rsd
                do loop_j = 1, 2*basis_dimension
                    do loop_k = 1, 2*basis_dimension
                        Rsd(loop_i-(nodiis-subsp), loop_j, loop_k) = rou_m(loop_j, loop_k) - rou_pre(loop_j, loop_k)
                    end do
                end do
                ! update rou_history
                do loop_j = 1, 2*basis_dimension
                    do loop_k = 1, 2*basis_dimension
                        rou_history(loop_i-(nodiis-subsp), loop_j, loop_k) = rou_m(loop_j, loop_k)
                    end do
                end do
                write(60,'(a,i2,a,i2)') '   --- DIIS subspace filling ',loop_i-(nodiis-subsp),'/',subsp
            else
                ! update Rsd
                do loop_j = 2, subsp
                    do loop_k = 1, 2*basis_dimension
                        do loop_l = 1, 2*basis_dimension
                            Rsd(loop_j - 1, loop_k, loop_l) = Rsd(loop_j, loop_k, loop_l)
                        end do
                    end do
                end do
                do loop_j = 1, 2*basis_dimension
                    do loop_k = 1, 2*basis_dimension
                        Rsd(subsp, loop_j, loop_k) = rou_m(loop_j, loop_k) - rou_pre(loop_j, loop_k)
                    end do
                end do
                ! update rou_history
                do loop_j = 2, subsp
                    do loop_k = 1, 2*basis_dimension
                        do loop_l = 1, 2*basis_dimension
                            rou_history(loop_j - 1, loop_k, loop_l) = rou_history(loop_j, loop_k, loop_l)
                        end do
                    end do
                end do
                do loop_j = 1, 2*basis_dimension
                    do loop_k = 1, 2*basis_dimension
                        rou_history(subsp, loop_j, loop_k) = rou_m(loop_j, loop_k)
                    end do
                end do
                ! construct DIISmat
                DIISmat = cmplx(0.0,0.0,dp)
                do loop_j = 1, subsp
                    do loop_k = 1, subsp
                        do loop_l = 1, 2*basis_dimension
                            do loop_m = 1, 2*basis_dimension
                                DIISmat(loop_j, loop_k) = DIISmat(loop_j, loop_k) + conjg(Rsd(loop_j, loop_l, loop_m))*Rsd(loop_k, loop_l, loop_m)
                            end do
                        end do
                    end do
                end do
                do loop_j = 1, subsp
                    DIISmat(subsp+1, loop_j) = cmplx(1.0,0.0,dp)
                    DIISmat(loop_j, subsp+1) = cmplx(1.0,0.0,dp)
                end do
                ! solveg residual equation
                ! dgesv and dspsv will cause V_Integral_2e conflict for unknown reason
                ! since DIISmat (and its inverse) is real symmetric, plus the column vector is simple, use the inverse of DIISmat to solve directly
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
                call zgetri(subsp+1, DIISmat, subsp+1, ipiv, DIISwork, lDIISwork, DIISinfo)
                deallocate(DIISwork)
                if (DIISinfo < 0) then
                    call terminate('DIIS solution failure, illegal input of dgetri')
                else if (DIISinfo > 0) then
                    call terminate('DIIS solution failure, internal error of dgetri')
                end if
                ! generate new rou_m
                rou_m = cmplx(0.0_dp,0.0_dp,dp)
                do loop_j = 1, subsp
                    do loop_k = 1, 2*basis_dimension
                        do loop_l = 1, 2*basis_dimension
                            rou_m(loop_k,loop_l) = rou_m(loop_k,loop_l) + DIISmat(loop_j,subsp+1) * (rou_history(loop_j,loop_k,loop_l) + damp*Rsd(loop_j,loop_k,loop_l))
                        end do
                    end do
                end do
                write(60,'(a,f10.6,a,f10.6)') '   --- predicted residual', -real(DIISmat(subsp+1,subsp+1)), ',',-aimag(DIISmat(subsp+1,subsp+1))
                do loop_j = 1, subsp
                    write(60,'(a,i2,f10.6,a,f10.6)') '   --- subsp coe',loop_j, real(DIISmat(loop_j,subsp+1)), ',', aimag(DIISmat(loop_j,subsp+1))
                end do
            end if
        end do
        write(60,*)
        write(60,*)
        if (abs(molE - molE_pre) < conver_tol .and. loop_i < maxiter) then
            if (DKH_order == 0) write(60,'(a)') '   DKH0 Hartree-Fock SCF succeed!'
            if (DKH_order == 2) write(60,'(a)') '   DKH2 Hartree-Fock SCF succeed!'
        else
            if (DKH_order == 0) write(60,'(a)') '   DKH0 Hartree-Fock SCF failed!'
            if (DKH_order == 2) write(60,'(a)') '   DKH2 Hartree-Fock SCF failed!'
        end if
        if (d4) then
            emd4 = dftd4('hf')
            molE = molE + emd4
        end if
        ! print final wave function information
        write(60,*)
        write(60,*)
        write(60,'(a)') '   -------------------------------------------------------------'
        write(60,'(a)') '                       MOLECULAR INFO'
        write(60,'(a)') '   -------------------------------------------------------------'
        write(60,'(a,f12.6)') '   total electronic energy / Eh                  ...',molE
        write(60,'(a,f12.6)') '   nuclear repulsive energy / Eh                 ...',nucE
        write(60,'(a,f12.6)') '   electron repulsive energy / Eh                ...',eleE
        write(60,'(a,f12.6)') '   electron kinetic energy / Eh                  ...',T
        write(60,'(a,f12.6)') '   electron-nuclear attraction energy / Eh       ...',V
        if (DKH_order == 2) then
            write(60,'(a,f12.6)') '   spin-orbital coupling energy / Eh             ...',ESOC
            if (SRTP_type) write(60,'(a,f12.6)') '   SRTP energy / Eh                              ...',ESR
        end if
        if (d4) write(60,'(a,f12.6)') '   dispersion energy (DFT-D4) / Eh               ...',emd4
        if (DKH_order == 0) then
            write(60,'(a,f12.6)') '   Virial ratio                                  ...',-(molE-T)/T
        else if (DKH_order == 2) then
            write(60,'(a,f12.6)') '   Virial ratio                                  ...',-(molE-T-ESOC)/T
            write(60,*)
            write(60,'(a)') '   Note: relativistic calculation cause the Virial ratio'
            write(60,'(a)') '   to deviate (usually below) 2.0.'
            write(60,*)
        end if
        write(60,'(a,f12.6)') '   total alpha electron                          ...',totalpha
        write(60,'(a,f12.6)') '   total alpha electron                          ...',totbeta
        write(60,'(a,f12.6)') '   <Sz*(Sz+1)> / hbar**2                         ...',((totalpha-totbeta)/2.0)*((totalpha-totbeta)/2.0+1.0_dp)
        write(60,'(a,f12.6)') '   <S**2> / hbar**2                              ...',S__2
        if (DKH_order == 2) then
            write(60,*)
            write(60,'(a)') '   Note: <S**2> may be contaminated by electron correlation,'
            write(60,'(a)') '   discussion of SOC suggests <S**2> data of 1e orbitals.'
        end if
        write(60,*)
        write(60,*)
        write(60,'(a)') '   -------------------------------------------------------------'
        write(60,'(a)') '                     CANONICAL MO INFO'
        write(60,'(a)') '   -------------------------------------------------------------'
        do loop_i = 1, electron_count
            if (loop_i < electron_count) then
                call calc_S2HForb(loop_i)
                write(60,'(a,i3.3,a,f12.6)') '   HOMO-', electron_count-loop_i, '   energy / Eh                       ... ', orbE(loop_i)
                write(60,'(a,f12.6)') '              <S**2> / hbar**2                  ... ', S__2orb
                write(60,'(a,f12.6)') '              <Sz> / hbar                       ... ', Szorb
                write(60,'(a)') '                 --- orb-in-atom  atom-in-mol      RE        IM'
                iatom = 1
                ishell = 1
                do loop_j = 1, basis_dimension
                    if (basis_inf(loop_j)%atom /= iatom) then
                        ishell = 1
                        iatom = iatom + 1
                    end if
                    if (abs(oper3(2*loop_j-1,loop_i)) >= prtlev) then
                        write(60,'(a,i2.2,a,i1,a,a2,a,i2.2,a,f9.6,a1,f9.6,a1)') '              ---  A #',ishell,' L=',basis_inf(loop_j)%L-1,'     ',&
                            element_list(molecular(basis_inf(loop_j)%atom)%atom_number),' #',basis_inf(loop_j)%atom,'    (', real(oper3(2*loop_j-1,loop_i)), ',', aimag(oper3(2*loop_j-1,loop_i)), ')'
                    end if
                    if (abs(oper3(2*loop_j,loop_i)) >= prtlev) then
                        write(60,'(a,i2.2,a,i1,a,a2,a,i2.2,a,f9.6,a1,f9.6,a1)') '              ---  B #',ishell,' L=',basis_inf(loop_j)%L-1,'     ',&
                            element_list(molecular(basis_inf(loop_j)%atom)%atom_number),' #',basis_inf(loop_j)%atom,'    (', real(oper3(2*loop_j,loop_i)), ',', aimag(oper3(2*loop_j,loop_i)), ')'
                    end if
                    ishell = ishell + 1
                end do
            else
                call calc_S2HForb(loop_i)
                write(60,'(a,a,f12.6)') '   HOMO    ', '   energy / Eh                       ... ', orbE(loop_i)
                write(60,'(a,f12.6)') '              <S**2> / hbar**2                  ... ', S__2orb
                write(60,'(a,f12.6)') '              <Sz> / hbar                       ... ', Szorb
                write(60,'(a)') '                 --- orb-in-atom  atom-in-mol      RE        IM'
                iatom = 1
                ishell = 1
                do loop_j = 1, basis_dimension
                    if (basis_inf(loop_j)%atom /= iatom) then
                        ishell = 1
                        iatom = iatom + 1
                    end if
                    if (abs(oper3(2*loop_j-1,loop_i)) >= prtlev) then
                        write(60,'(a,i2.2,a,i1,a,a2,a,i2.2,a,f9.6,a1,f9.6,a1)') '              ---  A #',ishell,' L=',basis_inf(loop_j)%L-1,'     ',&
                            element_list(molecular(basis_inf(loop_j)%atom)%atom_number),' #',basis_inf(loop_j)%atom,'    (', real(oper3(2*loop_j-1,loop_i)), ',', aimag(oper3(2*loop_j-1,loop_i)), ')'
                    end if
                    if (abs(oper3(2*loop_j,loop_i)) >= prtlev) then
                        write(60,'(a,i2.2,a,i1,a,a2,a,i2.2,a,f9.6,a1,f9.6,a1)') '              ---  B #',ishell,' L=',basis_inf(loop_j)%L-1,'     ',&
                            element_list(molecular(basis_inf(loop_j)%atom)%atom_number),' #',basis_inf(loop_j)%atom,'    (', real(oper3(2*loop_j,loop_i)), ',', aimag(oper3(2*loop_j,loop_i)), ')'
                    end if
                    ishell = ishell + 1
                end do
            end if
        end do
        do loop_i = electron_count+1, electron_count+5
            if (loop_i == electron_count + 1) then
                call calc_S2HForb(loop_i)
                write(60,'(a,a,f12.6)') '   LUMO    ', '   energy / Eh                       ... ', orbE(loop_i)
                write(60,'(a,f12.6)') '              <S**2> / hbar**2                  ... ', S__2orb
                write(60,'(a,f12.6)') '              <Sz> / hbar                       ... ', Szorb
                write(60,'(a)') '                 --- orb-in-atom  atom-in-mol      RE        IM'
                iatom = 1
                ishell = 1
                do loop_j = 1, basis_dimension
                    if (basis_inf(loop_j)%atom /= iatom) then
                        ishell = 1
                        iatom = iatom + 1
                    end if
                    if (abs(oper3(2*loop_j-1,loop_i)) >= prtlev) then
                        write(60,'(a,i2.2,a,i1,a,a2,a,i2.2,a,f9.6,a1,f9.6,a1)') '              ---  A #',ishell,' L=',basis_inf(loop_j)%L-1,'     ',&
                            element_list(molecular(basis_inf(loop_j)%atom)%atom_number),' #',basis_inf(loop_j)%atom,'    (', real(oper3(2*loop_j-1,loop_i)), ',', aimag(oper3(2*loop_j-1,loop_i)), ')'
                    end if
                    if (abs(oper3(2*loop_j,loop_i)) >= prtlev) then
                        write(60,'(a,i2.2,a,i1,a,a2,a,i2.2,a,f9.6,a1,f9.6,a1)') '              ---  B #',ishell,' L=',basis_inf(loop_j)%L-1,'     ',&
                            element_list(molecular(basis_inf(loop_j)%atom)%atom_number),' #',basis_inf(loop_j)%atom,'    (', real(oper3(2*loop_j,loop_i)), ',', aimag(oper3(2*loop_j,loop_i)), ')'
                    end if
                    ishell = ishell + 1
                end do
            else
                call calc_S2HForb(loop_i)
                write(60,'(a,i3.3,a,f12.6)') '   LUMO+', loop_i-electron_count-1, '   energy / Eh                       ... ', orbE(loop_i)
                write(60,'(a,f12.6)') '              <S**2> / hbar**2                  ... ', S__2orb
                write(60,'(a,f12.6)') '              <Sz> / hbar                       ... ', Szorb
                write(60,'(a)') '                 --- orb-in-atom  atom-in-mol      RE        IM'
                iatom = 1
                ishell = 1
                do loop_j = 1, basis_dimension
                    if (basis_inf(loop_j)%atom /= iatom) then
                        ishell = 1
                        iatom = iatom + 1
                    end if
                    if (abs(oper3(2*loop_j-1,loop_i)) >= prtlev) then
                        write(60,'(a,i2.2,a,i1,a,a2,a,i2.2,a,f9.6,a1,f9.6,a1)') '              ---  A #',ishell,' L=',basis_inf(loop_j)%L-1,'     ',&
                            element_list(molecular(basis_inf(loop_j)%atom)%atom_number),' #',basis_inf(loop_j)%atom,'    (', real(oper3(2*loop_j-1,loop_i)), ',', aimag(oper3(2*loop_j-1,loop_i)), ')'
                    end if
                    if (abs(oper3(2*loop_j,loop_i)) >= prtlev) then
                        write(60,'(a,i2.2,a,i1,a,a2,a,i2.2,a,f9.6,a1,f9.6,a1)') '              ---  B #',ishell,' L=',basis_inf(loop_j)%L-1,'     ',&
                            element_list(molecular(basis_inf(loop_j)%atom)%atom_number),' #',basis_inf(loop_j)%atom,'    (', real(oper3(2*loop_j,loop_i)), ',', aimag(oper3(2*loop_j,loop_i)), ')'
                    end if
                    ishell = ishell + 1
                end do
            end if
        end do
        AO2MO = oper3
        write(60,*)
        write(60,'(a)') 'exit module SCF'
        
        ! initialization
        deallocate(exi_T_j)
        deallocate(exi_V_j)
        deallocate(oper)
        deallocate(oper4)
        deallocate(oper6)
        deallocate(isupp_ev_f)
        deallocate(Rsd)
        deallocate(DIISmat)
        deallocate(ipiv)
        deallocate(rou_history)
        deallocate(rou_pre)
        deallocate(Fock2_mic)
        deallocate(Fock2_assigned)
        deallocate(swintegral)
        deallocate(Fock1)
        deallocate(Fock2)
        deallocate(i_V_j)
        deallocate(i_j)
        deallocate(i_p2_j)
        deallocate(S0_5)
        deallocate(exS0_5)
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
        if (kill) then
            call terminate('normal')
        else
            call terminate('keep')
        end if
    end subroutine DKH_Hartree_Fock
    
!------------------------------------------------------------
! gives an initial density matrix based on charge and spin multiplicity
    subroutine assign_rou()
        integer :: ploop_i,ploop_j,ploop_k,ploop_l                      ! ploop is loop variables only for subroutine assign_rou
        integer :: mat_dimension
        integer :: Na, Nb, degenlow, degenhigh, load, unload            ! degenerat area
        real(dp) :: rdMO(5)
        character(len = 512) :: line_str
        if (ini_rou) then
            ini_rou = .false.
            allocate(rou_m(2*basis_dimension,2*basis_dimension))
            Nalpha = (electron_count - (spin_mult - 1)) / 2 + (spin_mult - 1)
            Nbeta = (electron_count - (spin_mult - 1)) / 2
            
            ! read MO coefficient
            if (guess_type == 'gaussian') then
                allocate(AO2MO(2*basis_dimension,2*basis_dimension))
                allocate(Gaualpha(basis_dimension*basis_dimension))
                allocate(Gaubeta(basis_dimension*basis_dimension))
                allocate(AO2MOalpha(basis_dimension,basis_dimension))
                allocate(AO2MObeta(basis_dimension,basis_dimension))
                open(61, file = trim(address_job)//'.gaualpha', status = 'old', action = 'read', iostat = ios)
                if (ios /= 0) then
                    call system('rwfdump '//trim(address_job)//'.chk '//trim(address_job)//'.gaualpha 524R')
                    open(61, file = trim(address_job)//'.gaualpha', status = 'old', action = 'read', iostat = ios)
                    if (ios /= 0) call terminate('Cannot generate Gaussian alpha orbital for initial density matrix')
                end if
                open(62, file = trim(address_job)//'.gaubeta', status = 'old', action = 'read', iostat = ios)
                if (ios /= 0) then
                    call system('rwfdump '//trim(address_job)//'.chk '//trim(address_job)//'.gaubeta 526R')
                    open(62, file = trim(address_job)//'.gaubeta', status = 'old', action = 'read', iostat = ios)
                    if (ios /= 0) call terminate('Cannot generate Gaussian beta orbital for initial density matrix')
                end if
                do
                    read(61,'(a512)') line_str
                    if (index(line_str,'Dump of file') /= 0) exit
                end do
                if (index(line_str(index(line_Str,'length')-4:index(line_Str,'length')-2),'524') == 0) call terminate('Gaussian alpha orbital coefficient should be RWF 524')
                read(line_str(index(line_Str,'(')-9:index(line_Str,'(')-2),'(I)',iostat = ios) mat_dimension
                if (ios /= 0) call terminate('Unmatched Gaualpha content')
                if (mat_dimension /= basis_dimension*basis_dimension) then
                    write(*,*) mat_dimension, basis_dimension
                    call terminate('Gaussian alpha orbital dimension not match')
                end if
                do ploop_i=1,basis_dimension*basis_dimension/5
                    read(61,*) rdMO
                    Gaualpha((ploop_i-1)*5 + 1) = rdMO(1)
                    Gaualpha((ploop_i-1)*5 + 2) = rdMO(2)
                    Gaualpha((ploop_i-1)*5 + 3) = rdMO(3)
                    Gaualpha((ploop_i-1)*5 + 4) = rdMO(4)
                    Gaualpha((ploop_i-1)*5 + 5) = rdMO(5)
                end do
                ploop_i = ploop_i - 5
                if (mod(basis_dimension*basis_dimension,5) /= 0) then
                    if (mod(basis_dimension*basis_dimension,5) == 1) read(61,*) Gaualpha(ploop_i*5 + 1)
                    if (mod(basis_dimension*basis_dimension,5) == 2) read(61,*) Gaualpha(ploop_i*5 + 1),Gaualpha(ploop_i*5 + 2)
                    if (mod(basis_dimension*basis_dimension,5) == 3) read(61,*) Gaualpha(ploop_i*5 + 1),Gaualpha(ploop_i*5 + 2),Gaualpha(ploop_i*5 + 3)
                    if (mod(basis_dimension*basis_dimension,5) == 4) read(61,*) Gaualpha(ploop_i*5 + 1),Gaualpha(ploop_i*5 + 2),Gaualpha(ploop_i*5 + 3),Gaualpha(ploop_i*5 + 4)
                end if
                AO2MOalpha = transpose(reshape(Gaualpha,[basis_dimension,basis_dimension]))
                close(61)
                do
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
                ploop_i = ploop_i - 5
                if (mod(basis_dimension*basis_dimension,5) /= 0) then
                    if (mod(basis_dimension*basis_dimension,5) == 1) read(62,*) Gaubeta(ploop_i*5 + 1)
                    if (mod(basis_dimension*basis_dimension,5) == 2) read(62,*) Gaubeta(ploop_i*5 + 1),Gaubeta(ploop_i*5 + 2)
                    if (mod(basis_dimension*basis_dimension,5) == 3) read(62,*) Gaubeta(ploop_i*5 + 1),Gaubeta(ploop_i*5 + 2),Gaubeta(ploop_i*5 + 3)
                    if (mod(basis_dimension*basis_dimension,5) == 4) read(62,*) Gaubeta(ploop_i*5 + 1),Gaubeta(ploop_i*5 + 2),Gaubeta(ploop_i*5 + 3),Gaubeta(ploop_i*5 + 4)
                end if
                AO2MObeta = transpose(reshape(Gaubeta,[basis_dimension,basis_dimension]))
                close(62)
                deallocate(Gaualpha)
                deallocate(Gaubeta)
                AO2MO = cmplx(0.0_dp,0.0_dp,dp)
                do ploop_i = 1, basis_dimension
                    do ploop_j = 1, basis_dimension
                        AO2MO(2*ploop_i-1,2*ploop_j-1) = AO2MOalpha(ploop_i,ploop_j)
                        AO2MO(2*ploop_i,2*ploop_j) = AO2MObeta(ploop_i,ploop_j)
                    end do
                end do
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
                deallocate(AO2MOalpha)
                deallocate(AO2MObeta)
            else if (guess_type == 'read') then
                call read_matrix_cmplx(name='ao2mo', m=AO2MO, dm=mat_dimension)
                if (mat_dimension /= 2*basis_dimension) call terminate('basis dimension in .ao2mo file mismatch with current job')
                rou_m = cmplx(0.0_dp,0.0_dp,dp)
                do ploop_i = 1, 2*basis_dimension
                    do ploop_j = 1, 2*basis_dimension
                        do ploop_k = 1, electron_count
                            rou_m(ploop_i,ploop_j) = rou_m(ploop_i,ploop_j) + AO2MO(ploop_k,ploop_i)*conjg(AO2MO(ploop_k,ploop_j))
                        end do
                    end do
                end do
            end if
        else
            rou_pre = rou_m
            rou_m = cmplx(0.0_dp,0.0_dp,dp)
            do ploop_i = 1, 2*basis_dimension
                do ploop_j = 1, 2*basis_dimension
                    do ploop_k = 1, electron_count
                        rou_m(ploop_i,ploop_j) = rou_m(ploop_i,ploop_j) + AO2MO(ploop_k,ploop_i)*conjg(AO2MO(ploop_k,ploop_j))
                    end do
                end do
            end do
            call calc_S2HF()
            if (keepspin .and. abs((totalpha-totbeta) - real(Nalpha-Nbeta,dp)) > 1.0 .and. abs((totalpha-real(Nalpha,dp)) - (real(Nbeta,dp)-totbeta)) < 0.1) then
                write(60,'(a)') '   --- order of degenerate frontier alpha/beta orbitals changed'
                do ploop_i = electron_count-1, 1, -1
                    if (abs(orbE(ploop_i)-orbE(electron_count)) / abs(orbE(electron_count)) > 0.04) then
                        degenlow = ploop_i + 1
                        exit
                    end if
                end do
                do ploop_i = electron_count+1, 6*electron_count
                    if (abs(orbE(ploop_i)-orbE(electron_count)) / abs(orbE(electron_count)) > 0.04) then
                        degenhigh = ploop_i - 1
                        exit
                    end if
                end do
                write(60,'(a,i3,a,i3)') '   --- --- degenerate space: ',degenlow,' - ',degenhigh
                allocate(rotation(2*basis_dimension,degenhigh-degenlow+1))
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
                            do ploop_i = 1, 2*basis_dimension
                                rotation(ploop_i,load-degenlow+1) = oper3(ploop_i,ploop_k)
                            end do
                            write(60,'(a,i3,a,i3,a)') '   --- --- load ',ploop_k,' -> ',load,' (alpha)'
                            load = load + 1
                            Na = Na + 1
                        else
                            do ploop_i = 1, 2*basis_dimension
                                rotation(ploop_i,unload-degenlow+1) = oper3(ploop_i,ploop_k)
                            end do
                            write(60,'(a,i3,a,i3,a)') '   --- --- unload ',ploop_k,' -> ',unload,' (alpha)'
                            unload = unload - 1
                        end if                            
                    else
                            if (Nb < Nbeta) then
                            do ploop_i = 1, 2*basis_dimension
                                rotation(ploop_i,load-degenlow+1) = oper3(ploop_i,ploop_k)
                            end do
                            write(60,'(a,i3,a,i3,a)') '   --- --- load ',ploop_k,' -> ',load,' (beta)'
                            load = load + 1
                            Nb = Nb + 1
                        else
                            do ploop_i = 1, 2*basis_dimension
                                rotation(ploop_i,unload-degenlow+1) = oper3(ploop_i,ploop_k)
                            end do
                            write(60,'(a,i3,a,i3,a)') '   --- --- unload ',ploop_k,' -> ',unload,' (beta)'
                            unload = unload - 1
                        end if 
                    end if
                end do
                do ploop_k = degenlow, degenhigh
                    do ploop_i = 1, 2*basis_dimension
                        oper3(ploop_i,ploop_k) = rotation(ploop_i,ploop_k-degenlow+1)
                    end do
                end do
                deallocate(rotation)
                call zgemm( 'T', 'N', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), oper3, 2*basis_dimension, exS0_5, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), AO2MO, 2*basis_dimension)
                rou_m = cmplx(0.0_dp,0.0_dp,dp)
                do ploop_i = 1, 2*basis_dimension
                    do ploop_j = 1, 2*basis_dimension
                        do ploop_k = 1, electron_count
                            rou_m(ploop_i,ploop_j) = rou_m(ploop_i,ploop_j) + AO2MO(ploop_k,ploop_i)*conjg(AO2MO(ploop_k,ploop_j))
                        end do
                    end do
                end do
            end if
            if (loop_i /= 1) then
                maxDP = (real(rou_m(1,1))-real(rou_pre(1,1)))**2 + (aimag(rou_m(1,1))-aimag(rou_pre(1,1)))**2
                RMSDP = 0.0_dp
                do ploop_i = 1, 2*basis_dimension
                    do ploop_j = 1, 2*basis_dimension
                        if ((real(rou_m(ploop_i,ploop_j))-real(rou_pre(ploop_i,ploop_j)))**2 + (aimag(rou_m(ploop_i,ploop_j))-aimag(rou_pre(ploop_i,ploop_j)))**2 > maxDP) then
                            maxDP = (real(rou_m(ploop_i,ploop_j))-real(rou_pre(ploop_i,ploop_j)))**2 + (aimag(rou_m(ploop_i,ploop_j))-aimag(rou_pre(ploop_i,ploop_j)))**2
                        end if
                        RMSDP = RMSDP + (real(rou_m(ploop_i,ploop_j))-real(rou_pre(ploop_i,ploop_j)))**2 + (aimag(rou_m(ploop_i,ploop_j))-aimag(rou_pre(ploop_i,ploop_j)))**2
                    end do
                end do
                maxDP = sqrt(maxDP)
                RMSDP = RMSDP / (4*basis_dimension*basis_dimension)
                RMSDP = sqrt(RMSDP)
                write(60,'(a,e12.6)') '   --- maxDP ', maxDP
                write(60,'(a,e12.6)') '   --- RMSDP ', RMSDP
            end if
            ! calc <S**2>
            call calc_S2HF()
            write(60,'(a,f8.5)') '   --- <S**2> ',S__2
            write(60,'(a,f9.5)') '   --- total alpha electron ',totalpha
            write(60,'(a,f9.5)') '   --- total beta electron ',totbeta
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
                exSOC = ARVRA
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
            call zgemm( 'N', 'N', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), exAO2p2, 2*basis_dimension, exi_V_j, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), oper3, 2*basis_dimension)
            call zgemm( 'N', 'T', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), oper3, 2*basis_dimension, exAO2p2, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), exi_V_j, 2*basis_dimension)
            call zgemm( 'N', 'N', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), exAO2p2, 2*basis_dimension, exSOC, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), oper3, 2*basis_dimension)
            call zgemm( 'N', 'T', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), oper3, 2*basis_dimension, exAO2p2, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), exSOC, 2*basis_dimension)
            if (SRTP_type) then
                call zgemm( 'N', 'N', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), exAO2p2, 2*basis_dimension, exSR, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), oper3, 2*basis_dimension)
                call zgemm( 'N', 'T', 2*basis_dimension, 2*basis_dimension, 2*basis_dimension, cmplx(1.0_dp,0.0_dp,dp), oper3, 2*basis_dimension, exAO2p2, 2*basis_dimension, cmplx(0.0_dp,0.0_dp,dp), exSR, 2*basis_dimension)
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
        integer :: atom_i                                           ! which atom is the i^th component of |AOi>
        integer :: shell_i                                          ! which shell is the i^th component of |AOi>
        integer :: L_i                                              ! angular quantum number of |AOi>
        integer :: M_i                                              ! magnetic quantum number of |AOi>
        !----------------------------------
        integer :: contraction_j                                    ! contraction of atom_j, shell_j
        integer :: atom_j                                           ! which atom is the j_th component of |AOj>
        integer :: shell_j                                          ! which shell is the j_th component of |AOj>
        integer :: L_j                                              ! angular quantum number of |AOj>
        integer :: M_j                                              ! magnetic quantum number of |AOj>
        !----------------------------------
        integer :: contraction_k                                    ! contraction of atom_k, shell_k
        integer :: atom_k                                           ! which atom is the i^th component of |AOk>
        integer :: shell_k                                          ! which shell is the i^th component of |AOk>
        integer :: L_k                                              ! angular quantum number of |AOk>
        integer :: M_k                                              ! magnetic quantum number of |AOk>
        !----------------------------------
        integer :: contraction_l                                    ! contraction of atom_l, shell_l
        integer :: atom_l                                           ! which atom is the i^th component of |AOl>
        integer :: shell_l                                          ! which shell is the i^th component of |AOl>
        integer :: L_l                                              ! angular quantum number of |AOl>
        integer :: M_l                                              ! magnetic quantum number of |AOl>
        real(dp) :: exponents_i(20)                                 ! exponents of |AOi>
        real(dp) :: exponents_j(20)                                 ! exponents of |AOj>
        real(dp) :: exponents_k(20)                                 ! exponents of |AOk>
        real(dp) :: exponents_l(20)                                 ! exponents of |AOl>
        real(dp) :: coefficient_i(20)                               ! coefficient of |AOi>
        real(dp) :: coefficient_j(20)                               ! coefficient of |AOj>
        real(dp) :: coefficient_k(20)                               ! coefficient of |AOk>
        real(dp) :: coefficient_l(20)                               ! coefficient of |AOl>
        real(dp) :: coordinate_i(3)                                 ! coordinate of center of |AOi>
        real(dp) :: coordinate_j(3)                                 ! coordinate of center of |AOj>
        real(dp) :: coordinate_k(3)                                 ! coordinate of center of |AOk>
        real(dp) :: coordinate_l(3)                                 ! coordinate of center of |AOl>
        real(dp) :: swintegral_mic
        integer :: swloop_i, swloop_j                               ! swloop is loop variables only for schwarz screening
        integer :: swloop_m, swloop_n, swloop_o, swloop_p           ! swloop is loop variables only for schwarz screening
        integer :: shell_start_i                                    ! start point of same angular momentum shells
        integer :: shell_start_j                                    ! start point of same angular momentum shells
        real(dp),allocatable :: swexponents_i(:)                    ! exponents of |AOi> for schwarz screening
        real(dp),allocatable :: swexponents_j(:)                    ! exponents of |AOj> for schwarz screening
        real(dp),allocatable :: swcoefficient_i(:)                  ! coefficient of |AOi> for schwarz screening
        real(dp),allocatable :: swcoefficient_j(:)                  ! coefficient of |AOj> for schwarz screening
        if (ndschwarz) then
            ndschwarz = .false.
            write(60,'(a)') '   --- Schwarz screening of <ij||kl>'
            allocate(Fock2(2*basis_dimension,2*basis_dimension))
            allocate(swintegral(basis_dimension,basis_dimension))
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
            write(60,'(a,e10.2,a)') '   --- complete! cutoff:',schwarz_VT,'; stored in swintegral'
            allocate(Fock2_assigned(basis_dimension,basis_dimension))
            allocate(Fock2_mic(2*basis_dimension,2*basis_dimension))
            allocate(oper(2*basis_dimension,2*basis_dimension))
        end if
        Fock2 = cmplx(0.0_dp,0.0_dp,dp)
        ! parallel zone, running results consistent with serial
        !$omp parallel num_threads(threads_use) default(shared) private(i,dloop_i,dloop_j,dloop_k,dloop_l,dloop_m,dloop_n,dloop_o,dloop_p,integral, &
        !$omp& contraction_i,L_i,M_i,contraction_j,L_j,M_j, &
        !$omp& contraction_k,L_k,M_k,contraction_l,L_l,M_l, &
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
    
    !------------------------------------------------------------
    ! calculate <S**2> based on oper3 generated by Hartree-Fock SCF
    ! divide all occupied orbitals into two sets of alpha, beta orbitals with their original coefficients, then use the UHF method to solve Lowdin <S**2>.
    subroutine calc_S2HF()
        integer cloop_i, cloop_j, cloop_k                                           ! cloop is loop variables only for calc_S2HF routine
        complex(dp) :: alal, albe, beal, bebe                                       ! components of pair density matrix
        complex(dp) :: supp1, supp2
        ! oper3 and i_j are orthogonally normalized
        totalpha = 0.0_dp
        totbeta = 0.0_dp
        do cloop_i = 1, electron_count
            do cloop_k = 1, 2*basis_dimension-1, 2
                totalpha = totalpha + real(conjg(oper3(cloop_k,cloop_i))*oper3(cloop_k,cloop_i),dp)
            end do
            do cloop_k = 2, 2*basis_dimension, 2
                totbeta = totbeta + real(conjg(oper3(cloop_k,cloop_i))*oper3(cloop_k,cloop_i),dp)
            end do
        end do
        alal = totalpha * (totalpha - 1.0_dp)
        bebe = totbeta * (totbeta - 1.0_dp)
        albe = cmplx(0.0,0.0,dp)
        do cloop_i = 1, electron_count
            do cloop_j = 1, electron_count
                supp1 = cmplx(0.0,0.0,dp)
                supp2 = cmplx(0.0,0.0,dp)
                do cloop_k = 1, 2*basis_dimension-1, 2
                    supp1 = supp1 + conjg(oper3(cloop_k,cloop_i))*oper3(cloop_k+1,cloop_j)
                    supp2 = supp2 + conjg(oper3(cloop_k+1,cloop_j))*oper3(cloop_k,cloop_i)
                end do
                albe = albe + supp1*supp2
            end do
        end do
        beal = cmplx(0.0,0.0,dp)
        do cloop_i = 1, electron_count
            do cloop_j = 1, electron_count
                supp1 = cmplx(0.0,0.0,dp)
                supp2 = cmplx(0.0,0.0,dp)
                do cloop_k = 2, 2*basis_dimension, 2
                    supp1 = supp1 + conjg(oper3(cloop_k,cloop_i))*oper3(cloop_k-1,cloop_j)
                    supp2 = supp2 + conjg(oper3(cloop_k-1,cloop_j))*oper3(cloop_k,cloop_i)
                end do
                beal = beal + supp1*supp2
            end do
        end do
        S__2 = -real(electron_count*(electron_count-4),dp)/4.0_dp + real(alal/2.0_dp + bebe/2.0_dp) - real(albe/2.0_dp + beal/2.0_dp)
    end subroutine calc_S2HF
    
    !------------------------------------------------------------
    ! calculate <S**2> of certain orbital based on oper3 generated by Hartree-Fock SCF
    subroutine calc_S2HForb(orbnum)
        integer,intent(in) :: orbnum
        integer :: zloop_i, zloop_j                                  ! zloop is loop variables only for calc_S2HForb routine
        real(dp) :: supp1, supp2
        totalphaorb = 0.0_dp
        do zloop_i = 1, 2*basis_dimension-1, 2
            totalphaorb = totalphaorb + real(conjg(oper3(zloop_i,orbnum))*oper3(zloop_i,orbnum),dp)
        end do
        totbetaorb = 0.0_dp
        do zloop_i = 2, 2*basis_dimension, 2
            totbetaorb = totbetaorb + real(conjg(oper3(zloop_i,orbnum))*oper3(zloop_i,orbnum),dp)
        end do
        supp1 = cmplx(0.0,0.0,dp)
        supp2 = cmplx(0.0,0.0,dp)
        do zloop_i = 2, 2*basis_dimension, 2
            supp1 = supp1 + conjg(oper3(zloop_i,orbnum))*oper3(zloop_i-1,orbnum)
            supp2 = supp2 + conjg(oper3(zloop_i-1,orbnum))*oper3(zloop_i,orbnum)
        end do
        S__2orb = 3.0/4.0 + totalphaorb*(totalphaorb-1.0)/2.0 + totbetaorb*(totbetaorb-1.0)/2.0 - real(supp1*supp2)
        Szorb = (totalphaorb-totbetaorb)/2.0_dp
    end subroutine calc_S2HForb
    
    !-----------------------------------------------------------------------
    ! dump molecular orbital information to .molden file
    subroutine dump_molden()
        integer :: channel, dmi, dmj, dmk
        if (.not. allocated(AO2MO)) call terminate('dump molecular orbital failed, AO2MO = NULL')
        
        
        ! molden file contains the real part of molecular orbital
        open(newunit=channel, file=trim(address_job)//'-real.molden', status='replace', action='write', iostat=ios)
        if (ios /= 0) call terminate('dump .molden failed')
        write(channel, '(a)') '[Molden Format]'
        write(channel, '(a)') '[Title]'
        write(channel, '(a)') 'molden file contains the real of molecular orbital of '//trim(address_molecular)
        write(channel, *)
        ! molecular geometry
        write(channel, '(a)') '[Atoms] AU'
        do dmi = 1, atom_count
            write(channel, '(a,i3,i2,f13.7,f13.7,f13.7)') element_list(molecular(dmi)%atom_number), &
                dmi, molecular(dmi)%atom_number, molecular(dmi)%nucleus_position(1), &
                molecular(dmi)%nucleus_position(2), molecular(dmi)%nucleus_position(3)
        end do
        ! basis of each atom
        write(channel, '(a)') '[GTO]'
        do dmi = 1, atom_count
            write(channel, '(i3,i2)') dmi, 0
            do dmj = 0, shell_in_element(molecular(dmi) % atom_number)-1
                if (atom_basis(molecular(dmi)%basis_number+dmj)%angular_quantum_number == 0) then
                    write(channel, '(a2,i2,a)') 's', atom_basis(molecular(dmi)%basis_number+dmj)%contraction,' 1.0'
                else if (atom_basis(molecular(dmi)%basis_number+dmj)%angular_quantum_number == 1) then
                    write(channel, '(a2,i2,a)') 'p', atom_basis(molecular(dmi)%basis_number+dmj)%contraction,' 1.0'
                else if (atom_basis(molecular(dmi)%basis_number+dmj)%angular_quantum_number == 2) then
                    write(channel, '(a2,i2,a)') 'd', atom_basis(molecular(dmi)%basis_number+dmj)%contraction,' 1.0'
                else if (atom_basis(molecular(dmi)%basis_number+dmj)%angular_quantum_number == 3) then
                    write(channel, '(a2,i2,a)') 'f', atom_basis(molecular(dmi)%basis_number+dmj)%contraction,' 1.0'
                else if (atom_basis(molecular(dmi)%basis_number+dmj)%angular_quantum_number == 4) then
                    write(channel, '(a2,i2,a)') 'g', atom_basis(molecular(dmi)%basis_number+dmj)%contraction,' 1.0'
                end if
                do dmk = 1, atom_basis(molecular(dmi)%basis_number+dmj)%contraction
                    write(channel, '(f13.7,f13.7)') atom_basis(molecular(dmi)%basis_number+dmj)%exponents(dmk), &
                        atom_basis(molecular(dmi)%basis_number+dmj)%coefficient(dmk)
                end do
            end do
            write(channel, *)
        end do
        write(channel, '(a)') '[6D]'
        write(channel, '(a)') '[10F]'
        write(channel, '(a)') '[15G]'
        ! MO coefficient
        write(channel, '(a)') '[MO]'
        do dmi = 1, basis_dimension
            write(channel, '(a4,e23.14)') 'Ene=', orbE(dmi)
            write(channel, '(a)') 'Spin= Alpha'
            call calc_S2HForb(dmi)
            if (dmi <= electron_count) then
                write(channel, '(a7,f12.6)') 'Occup=', 0.5 + Szorb
            else
                write(channel, '(a)') 'Occup= 0.00'
            end if
            do dmj = 1, basis_dimension
                write(channel, '(i4,f20.12)') dmj, real(AO2MO(2*dmj-1, dmi))
            end do
        end do
        do dmi = 1, basis_dimension
            write(channel, '(a4,e23.14)') 'Ene=', orbE(dmi)
            write(channel, '(a)') 'Spin= Beta'
            call calc_S2HForb(dmi)
            if (dmi <= electron_count) then
                write(channel, '(a7,f12.6)') 'Occup=', 0.5 - Szorb
            else
                write(channel, '(a)') 'Occup= 0.00'
            end if
            do dmj = 1, basis_dimension
                write(channel, '(i4,f20.12)') dmj, real(AO2MO(2*dmj, dmi))
            end do
        end do
        close(channel)
        
        ! molden file contains the maginary part of molecular orbital
        open(newunit=channel, file=trim(address_job)//'-img.molden', status='replace', action='write', iostat=ios)
        if (ios /= 0) call terminate('dump .molden failed')
        write(channel, '(a)') '[Molden Format]'
        write(channel, '(a)') '[Title]'
        write(channel, '(a)') 'molden file contains the maginary part of molecular orbital of '//trim(address_molecular)
        write(channel, *)
        ! molecular geometry
        write(channel, '(a)') '[Atoms] AU'
        do dmi = 1, atom_count
            write(channel, '(a,i3,i2,f13.7,f13.7,f13.7)') element_list(molecular(dmi)%atom_number), &
                dmi, molecular(dmi)%atom_number, molecular(dmi)%nucleus_position(1), &
                molecular(dmi)%nucleus_position(2), molecular(dmi)%nucleus_position(3)
        end do
        ! basis of each atom
        write(channel, '(a)') '[GTO]'
        do dmi = 1, atom_count
            write(channel, '(i3,i2)') dmi, 0
            do dmj = 0, shell_in_element(molecular(dmi) % atom_number)-1
                if (atom_basis(molecular(dmi)%basis_number+dmj)%angular_quantum_number == 0) then
                    write(channel, '(a2,i2,a)') 's', atom_basis(molecular(dmi)%basis_number+dmj)%contraction,' 1.0'
                else if (atom_basis(molecular(dmi)%basis_number+dmj)%angular_quantum_number == 1) then
                    write(channel, '(a2,i2,a)') 'p', atom_basis(molecular(dmi)%basis_number+dmj)%contraction,' 1.0'
                else if (atom_basis(molecular(dmi)%basis_number+dmj)%angular_quantum_number == 2) then
                    write(channel, '(a2,i2,a)') 'd', atom_basis(molecular(dmi)%basis_number+dmj)%contraction,' 1.0'
                else if (atom_basis(molecular(dmi)%basis_number+dmj)%angular_quantum_number == 3) then
                    write(channel, '(a2,i2,a)') 'f', atom_basis(molecular(dmi)%basis_number+dmj)%contraction,' 1.0'
                else if (atom_basis(molecular(dmi)%basis_number+dmj)%angular_quantum_number == 4) then
                    write(channel, '(a2,i2,a)') 'g', atom_basis(molecular(dmi)%basis_number+dmj)%contraction,' 1.0'
                end if
                do dmk = 1, atom_basis(molecular(dmi)%basis_number+dmj)%contraction
                    write(channel, '(f13.7,f13.7)') atom_basis(molecular(dmi)%basis_number+dmj)%exponents(dmk), &
                        atom_basis(molecular(dmi)%basis_number+dmj)%coefficient(dmk)
                end do
            end do
            write(channel, *)
        end do
        write(channel, '(a)') '[6D]'
        write(channel, '(a)') '[10F]'
        write(channel, '(a)') '[15G]'
        ! MO coefficient
        write(channel, '(a)') '[MO]'
        do dmi = 1, basis_dimension
            write(channel, '(a4,e23.14)') 'Ene=', orbE(dmi)
            write(channel, '(a)') 'Spin= Alpha'
            call calc_S2HForb(dmi)
            if (dmi <= electron_count) then
                write(channel, '(a7,f12.6)') 'Occup=', 0.5 + Szorb
            else
                write(channel, '(a)') 'Occup= 0.000'
            end if
            do dmj = 1, basis_dimension
                write(channel, '(i4,f20.12)') dmj, aimag(AO2MO(2*dmj-1, dmi))
            end do
        end do
        do dmi = 1, basis_dimension
            write(channel, '(a4,e23.14)') 'Ene=', orbE(dmi)
            write(channel, '(a)') 'Spin= Beta'
            call calc_S2HForb(dmi)
            if (dmi <= electron_count) then
                write(channel, '(a7,f12.6)') 'Occup=', 0.5 - Szorb
            else
                write(channel, '(a)') 'Occup= 0.000'
            end if
            do dmj = 1, basis_dimension
                write(channel, '(i4,f20.12)') dmj, aimag(AO2MO(2*dmj, dmi))
            end do
        end do
        close(channel)
    end subroutine dump_molden
    
end module SCF