! |-------<THOMAS RELATIVISTIC ELECTRONIC STRUCTURE CALCULATIONS>------|
! |---------------- -------------<DIRAC4PI>----------------------------|
! |-------------------------------<ATOMS>------------------------------|
! |--------<ATOMIC INFORMATION & CORRESPONDING ATOMIC BASIS SET>-------|

module Atoms
    use Fundamentals
    
	integer,parameter :: element_count = 30
    
	character(len = 2),parameter :: element_list(element_count) =                             &
  	(/ ' H',                                                                                'He', &  ! 1-2
       'Li','Be',                                                  ' B',' C',' N',' O',' F','Ne', &  ! 3-10
       'Na','Mg',                                                  'Al','Si',' P',' S','Cl','Ar', &  ! 11-18
       ' K','Ca','Sc','Ti',' V','Cr','Mn','Fe','Co','Ni','Cu','Zn'                               /)  ! 19-30
    
    integer,parameter :: element_massnumber(element_count) =                                  &
  	(/    2,                                                                                   4, &  ! 1-2
          7,   9,                                                    11,  12,  14,  16,  19,  20, &  ! 3-10
         23,  24,                                                    27,  28,  31,  32,  35,  40, &  ! 11-18
         39,  40,  45,  48,  51,  52,  55,  56,  59,  59,  64,  65                               /)  ! 19-30
    
    integer :: basis_dimension                                    ! number of basis in molecular
    integer :: basis_count                                        ! number of basis functions in basis set
    integer,allocatable :: shell_in_element(:)                    ! number of shells in each element
    type atom_basis_type
        integer :: atom_number                                    ! atomic number of basis atom
        integer :: angular_quantum_number                         ! angular quantum number of basis shell
        integer :: contraction                                    ! contraction of basis shell
        real(dp),allocatable :: exponents(:)                      ! exponents of primitive basis shell
        real(dp),allocatable :: coefficient(:)                    ! contraction coefficient of primitive basis shell
        real(dp),allocatable :: Ncoefficient(:,:)                 ! normalized coefficient of primitive shell
    end type atom_basis_type
    type(atom_basis_type),allocatable :: atom_basis(:)            ! basis functions CGTOs
    
    integer :: atom_count                                         ! number of atoms in molecular
    type molecular_type     			              
        integer :: basis_number                                   ! serial number of each atom in basis set
        integer :: atom_number                                    ! atomic number of each atom in molecular
        integer :: atom_charge                                    ! charge of the atom (for initial guess of density matrix)
        integer :: atom_spin                                      ! spin of the atom (for initial guess of density matrix)
        real(dp) :: nucleus_position(3)                           ! position of each atom, (x,y,z)
        real(dp) :: nucleus_radius                                ! nucleus radius of each atom in molecular 0.836*A^(1/3)+0.57 (fm)
    end type molecular_type
    type(molecular_type),allocatable :: molecular(:)              ! system composed only of atoms
    
    real(dp) :: S__2                                              ! <S**2> of molecular
    real(dp) :: totalpha, totbeta                                 ! total alpha and beta electron of molecular
    real(dp) :: S__2orb, Szorb                                    ! <S**2>, Sz of orbital
    real(dp) :: totalphaorb, totbetaorb                           ! total alpha and beta electron of orbital
    
    contains
    
!-----------------------------------------------------------------------
! Get CGTOs from .gbs file, sharing of exponents is not allowed 
! default: 6d cartesian type segment contraction basis set
! Spherical harmonic basis functions should be used for most basis sets, but too many code changes are required, author may make improvements in the future.
! Minimal basis sets are not available
    
    subroutine read_gbs(address_basis)
        implicit none
        integer :: try
        character(len = *),intent(in) :: address_basis
        character(len = 512) :: line_str
        character(len = 2) :: basis_element_name
        character(len = 1) :: basis_angular_name
        if(index(address_basis,'.gbs') == 0) then
            call terminate('input basis set file is not .gbs file')
        end if
        inquire(file = address_basis, exist = exists)
        if (.not. exists) then
            call terminate('basis set file address does not exist')
        end if
        open(11,file = address_basis,status="old",action="read")
        do
            read(11,'(a512)') line_str
            if (index(line_str,'!') == 0) then
                exit
            end if
        end do
        read(11,*)
        
        call get_basis_count(address_basis)
        allocate(atom_basis(basis_count))
        
        read(11,*) basis_element_name                                                              ! default for neutral atom basis set
        loop_i = 1
        loop_j = 1
        loop_m = 0                                                                                 ! count for basis function in each element
        element_loop: do
            if (ADJUSTL(basis_element_name) /= ADJUSTL(element_list(loop_i))) then
                loop_i = loop_i + 1
                if (loop_i >= element_count + 1) then
                    call terminate('basis set file: incorrect element symbol or element not included in current version of program')
                end if
                cycle
            end if
            atom_basis(loop_j) % atom_number = loop_i
            angular_loop: do
                read(11,*,iostat = ios) basis_angular_name, atom_basis(loop_j) % contraction       ! scale factor defaults to 1.00
                if (ios /= 0) then
                    call terminate('read basis set file failed, may caused by incorrect formatting')
                end if 
                if (basis_angular_name == 'S') then
                    atom_basis(loop_j) % angular_quantum_number = 0
                else if (basis_angular_name == 'P') then
                    atom_basis(loop_j) % angular_quantum_number = 1
                else if (basis_angular_name == 'D') then
                    atom_basis(loop_j) % angular_quantum_number = 2
                else if (basis_angular_name == 'F') then
                    atom_basis(loop_j) % angular_quantum_number = 3
                else if (basis_angular_name == 'G') then
                    atom_basis(loop_j) % angular_quantum_number = 4
                end if
                allocate(atom_basis(loop_j) % exponents(atom_basis(loop_j) % contraction))                         ! allocate memory according to contraction
                allocate(atom_basis(loop_j) % coefficient(atom_basis(loop_j) % contraction))                       ! allocate memory according to contraction
                allocate(atom_basis(loop_j) % Ncoefficient(atom_basis(loop_j) % contraction,(atom_basis&
                    (loop_j) % angular_quantum_number+2)*(atom_basis(loop_j) % angular_quantum_number+1)/2))       ! allocate memory according to contraction
                atom_basis(loop_j) % atom_number = loop_i
                loop_k = 1
                contraction_loop: do
                    if (loop_k <= atom_basis(loop_j) % contraction) then
                        read(11,*,iostat = ios) atom_basis(loop_j) % exponents(loop_k), atom_basis(loop_j) % coefficient(loop_k)
                        if (ios /= 0) then
                            call terminate('read basis set file failed, try to denote the scientific notation in basis file with E')
                        end if 
                        loop_k = loop_k + 1
                        cycle
                    end if
                    loop_j = loop_j + 1
                    exit
                end do contraction_loop 
                loop_m = loop_m + 1
                if (loop_m == shell_in_element(loop_i)) then
                    exit
                end if
            end do angular_loop
            read(11,*)                                                         ! skip the separator between elements
            loop_m = 0
            loop_i = 1
            read(11,*,iostat = ios) basis_element_name
            if (ios < 0) then
                exit
            end if
        end do element_loop
        close(11)
    end subroutine read_gbs
    
!-----------------------------------------------------------------------
! Get the number of basis functions in .gbs file
    
    subroutine get_basis_count(address_basis)
        implicit none
        character(len = *),intent(in) :: address_basis
        character(len = 512) :: line_str
        allocate(shell_in_element(50))
        loop_i = 0
        loop_j = 0
        open(13,file = address_basis)
        do
            read(13,'(a512)') line_str
            if (index(line_str,'!') == 0) then
                exit
            end if
        end do
        read(13,*)
        do
            read(13,'(a512)',iostat = ios) line_str
            if (ios < 0) then
                exit
            end if
            if (index(line_str,'S') /= 0 .or. index(line_str,'P') /= 0 .or. index(line_str,'D') /= 0 &
                .or. index(line_str,'F') /= 0 .or. index(line_str,'G') /= 0) then
                if (index(line_str,'1.00') /= 0 .or. index(line_str,'1.0') /= 0) then
                    loop_i = loop_i + 1
                    cycle
                end if
                end if
            if (index(line_str,'0') /= 0 .and. index(line_str,'.') == 0) then
                do loop_l = 1, element_count
                    if (ADJUSTL(line_str(1:2)) == ADJUSTL(element_list(loop_l))) then
                        loop_k = loop_l
                        exit
                    end if
                end do
                if (loop_l == element_count+1) then
                    call terminate('basis set file: incorrect element symbol or element not included in current version of program')
                end if
            end if
            if (index(line_str,'**') /= 0) then
                shell_in_element(loop_k) = loop_i - loop_j
                loop_j = loop_i
            end if  
        end do
        close(13)
        basis_count = loop_i
    end subroutine get_basis_count

    
!-----------------------------------------------------------------------
! Get static geometry of molecule from .tip file
    
    subroutine read_geometry(address_molecular)
        implicit none
        character(len = *),intent(in) :: address_molecular
        character(len = 2) :: molecular_element_name, title_note
        if(index(address_molecular,'.tip') == 0) then
            call terminate('input geometry file is not .tip file')
        end if
        inquire(file = address_molecular, exist = exists)
        if (.not. exists) then
            call terminate('geometry file address does not exist')
        end if
        open(12,file = address_molecular,status="old",action="read")
        read(12,*) title_note
        do
            read(title_note,"(I)",iostat = ios) atom_count
            if (ios == 0) then
                exit
            else if (ios /= 0 .and. index(title_note,'!') == 1) then
                read(12,*) title_note
                cycle
            else
                call terminate("use '!' tag for comments in .tip file")
            end if
        end do
        read(12,*)
        allocate(molecular(atom_count + 1))                                    ! allocate 1 additional group for inspection
        loop_i = 1
        loop_j = 1
        do while (loop_j <= atom_count)
            read(12,*,iostat = ios) molecular_element_name, molecular(loop_j) % &
            nucleus_position(1), molecular(loop_j) % nucleus_position(2), molecular(loop_j) % nucleus_position(3)
            molecular(loop_j) % nucleus_position(1) = molecular(loop_j) % nucleus_position(1) / Ang2Bohr ! change the unit from Angstorm to Bohr
            molecular(loop_j) % nucleus_position(2) = molecular(loop_j) % nucleus_position(2) / Ang2Bohr
            molecular(loop_j) % nucleus_position(3) = molecular(loop_j) % nucleus_position(3) / Ang2Bohr
            if (ios /= 0) then
                call terminate('read geometry file failed, may caused by the absence of atoms actually recorded')
            end if
            do
                if (ADJUSTL(molecular_element_name) /= ADJUSTL(element_list(loop_i))) then
                    loop_i = loop_i + 1
                    if (loop_i >= element_count + 1) then
                        call terminate('geometry file: incorrect element symbol or element not contained in current version of program')
                    end if
                    cycle
                end if
                exit
            end do
            molecular(loop_j) % atom_number = loop_i
            loop_k = 1
            do
                if (atom_basis(loop_k) % atom_number /= loop_i) then
                    loop_k = loop_k + 1
                    if (loop_k >= basis_count) then
                        call terminate('there are elements in geometry file not contained in basis set file')
                    end if
                    cycle
                end if
                exit
            end do
            molecular(loop_j) % basis_number = loop_k
            molecular(loop_j) % nucleus_radius = 0.836 * real(element_massnumber(loop_i)) ** (1/3) + 0.57
            loop_i = 1
            loop_j = loop_j + 1
        end do
        read(12,*,iostat = ios) molecular_element_name, molecular(atom_count+1) % &
        nucleus_position(1), molecular(atom_count+1) % nucleus_position(2), molecular(atom_count+1) % nucleus_position(3)
        if (ios == 0) then
            call terminate('read geometry file failed, may caused by the redundancy of atoms actually recorded')
        end if
        ! calculation settings will be read in module Fundamentals
        close(12)
        loop_i = 1
        electron_count = 0
        do while (loop_i <= atom_count)                                        ! calculate number of electrons in neutral system
            electron_count = electron_count + molecular(loop_i) % atom_number
            loop_i = loop_i + 1
        end do
        molecular % atom_charge = molecular % atom_number
        molecular % atom_spin = 0
    end subroutine read_geometry

!-----------------------------------------------------------------------
! print basis set information, geometry information and calculation settings
    
    subroutine input_print (address_basis, address_molecular)
        character(len = *),intent(in) :: address_basis
        character(len = *),intent(in) :: address_molecular
        write(60,"(a)") 'Module Atoms:'
        write(60,"(a)") '   input basis set file path: '//address_basis
        write(60,"(a)") '   -------------------------<BASIS>-------------------------'
        loop_i = 1
        do while(loop_i <= basis_count)                                                                      ! print basis set information, print only the basis set of atoms in the molecule
            loop_j = 1
            do while(loop_j <= atom_count)
                if (molecular(loop_j) % basis_number == loop_i) then
                    exit
                end if
                loop_j = loop_j + 1
            end do
            if (loop_j == atom_count + 1) then
                loop_i = loop_i + 1
                cycle
            end if
            write(60,"('   ' A2 '     0')") element_list(atom_basis(loop_i) % atom_number)                   ! default for neutral atom basis set
            loop_k = 0
            do while(loop_k <= shell_in_element(molecular(loop_j) % atom_number) - 1)
                write(60,"('                      l = ' I1)")  atom_basis(loop_i + loop_k) % angular_quantum_number
                loop_m = 1
                do while(loop_m <= atom_basis(loop_i + loop_k) % contraction)
                    write(60,"('   exp: ' F15.8 '      coe: ' F15.8)") atom_basis(loop_i + loop_k) % &
                    exponents(loop_m), atom_basis(loop_i + loop_k) % coefficient(loop_m)
                    loop_m = loop_m + 1
                end do
                loop_k = loop_k + 1
            end do
            loop_i = loop_i + 1
            write(60,*)
        end do
        write(60,*)
        write(60,"(a)") '   input geometry file path: '//address_molecular
        write(60,"(a)") '   -------------------------<GEOMETRY>-------------------------'
        loop_i = 1
        do while(loop_i <= atom_count)                                                                       ! print geometry information
            write(60,"('   ' A2 '    ' F9.4 '    ' F9.4 '    ' F9.4)") element_list(molecular(loop_i) % &
            atom_number), molecular(loop_i) % nucleus_position(1), molecular(loop_i) % &
            nucleus_position(2), molecular(loop_i) % nucleus_position(3)
            loop_i = loop_i + 1
        end do
        ! get the dimension of basis of molecular, use 6d, 10f, 18g basis for convenience
        basis_dimension = 0
        loop_i = 1
        do while(loop_i <= atom_count)
            loop_j = 1
            do while(loop_j <= shell_in_element(molecular(loop_i) % atom_number))
                basis_dimension = basis_dimension + &
                (atom_basis(molecular(loop_i) % basis_number + loop_j - 1) % angular_quantum_number + 2) * &
                (atom_basis(molecular(loop_i) % basis_number + loop_j - 1) % angular_quantum_number + 1) / 2
                loop_j = loop_j + 1
            end do
            loop_i = loop_i + 1
        end do
        write(60,"(A34,I4,A2,I4)") '   scalar/spinor basis dimension: ', basis_dimension, ' /', 2*basis_dimension
        write(64,rec = 1) 2*basis_dimension     ! first record in .ao2mo file is the dimension of AO2MO
        write(60,"(A)") '   Cartesian type basis (6D,10F,18G) used.'
    end subroutine input_print
    
    !-----------------------------------------------------------------------
    ! read keywords from .tip file
    
    subroutine read_keywords(address_molecular)
        implicit none
        character(len = *),intent(in) :: address_molecular
        character(len = 40) :: module_name, keyword
        character(len = 2) :: title_note
        integer :: atom_count
        integer :: atomcharge, atomspin, iatom
        write(60,"(a)") "   -------------------------<KEYWORDS>-------------------------"
        open(12, file=address_molecular, status="old", action="read")
        read(12,*) title_note
        do
            read(title_note,"(I)",iostat = ios) atom_count
            if (ios == 0) exit
            read(12,*) title_note
        end do
        read(12,*)
        loop_j = 1
        do while (loop_j <= atom_count)
            read(12,*)
            loop_j= loop_j + 1
        end do
        entire: do
            read(12,*,iostat = ios) module_name
            interval: do                                           ! Allow multiple blank lines between each module
                if (ios /= 0) then
                    exit entire
                else if (module_name == '' .or. index(module_name,'!') == 1) then
                    read(12,*,iostat = ios) module_name
                    cycle interval
                end if
                exit interval
            end do interval
            ! ----------------------<module Atoms>----------------------
            if (trim(module_name) == '%Atoms') then
            	module_Atoms: do
                	read(12,*,iostat = ios) keyword
                	if (ios /= 0 .or. keyword == '') then
                    	call terminate('empty line detected in module Atoms.')
                    else if (index(keyword,'!') == 1) then
                        cycle module_Atoms
                	else if (trim(keyword) == 'endAtoms') then
                		exit module_Atoms
                	else if (trim(keyword) == 'end') then
                		call terminate("end of module should be written as 'endmodulename'")
                	else if (index(keyword,'charge') == 1) then
                		if (index(keyword,'=') /= 0) then
                			read(keyword(index(keyword,'=') + 1 : len(trim(keyword))),"(I)",iostat = ios) charge
                			if (ios /= 0) call terminate("charge setting should be written as 'charge=n'")
                		else
                			call terminate("charge setting should be written as 'charge=n'")
                		end if
                		write(60,"(A,I2)") "   Charge is changed to ",charge
                	else if (index(keyword,'spin') == 1) then
                		if (index(keyword,'=') /= 0) then
                			read(keyword(index(keyword,'=') + 1 : len(trim(keyword))),"(I)",iostat = ios) spin_mult
                			if (ios /= 0) call terminate("spin multiplicity setting should be written as 'spin=n'")
                		else
                			call terminate("spin multiplicity setting should be written as 'spin=n'")
                		end if
                		write(60,"(A,I2)") "   Spin multiplicity is changed to ",spin_mult
                	else
	                	call terminate('unknown keyword detected in module Atoms')
                    end if
                    
            	end do module_Atoms 
            ! -------------------<module Hamiltonian>-------------------
            else if (trim(module_name) == '%Hamiltonian') then
            	module_Hamiltonian: do
	                read(12,*,iostat = ios) keyword
	                if (ios /= 0 .or. keyword == '') then
	                	call terminate('empty line detected in module Hamiltonian.')
                    else if (index(keyword,'!') == 1) then
                        cycle module_Hamiltonian
	                else if (trim(keyword) == 'endHamiltonian') then
	                	exit module_Hamiltonian
	                else if (trim(keyword) == 'end') then
	                	call terminate("end of module should be written as 'endmodulename'")
                    else if (index(keyword,'threads') == 1) then
                        if (index(keyword,'=') /= 0) then
                			read(keyword(index(keyword,'=') + 1 : len(trim(keyword))),"(I)",iostat = ios) threads_use
                			if (ios /= 0) call terminate("threads setting should be written as 'threads=n'")
                		else
                			call terminate("threads setting should be written as 'threads=n'")
                        end if
                    else if (index(keyword,'stack') == 1) then
                        if (index(keyword,'=') /= 0) then
                			read(keyword(index(keyword,'=') + 1 : len(trim(keyword))),"(I)",iostat = ios) stacksize
                			if (ios /= 0) call terminate("stack size setting should be written as 'stack=n'(MB)")
                		else
                			call terminate("stack size setting should be written as 'stack=n'(MB)")
                        end if
                    else if (trim(keyword) == 'dkh0') then
				        DKH_order = 0
				        write(60,"(a)") "   nonrelativistic Hamiltonian will be considered."
                    else if (trim(keyword) == 'dkh1') then
				        DKH_order = 1
				        write(60,"(a)") "   fpFW Hamiltonian will be considered."
                    else if (trim(keyword) == 'dkh2') then
				        DKH_order = 2
				        write(60,"(a)") "   DKH2 Hamiltonian will be considered."
                    else if (trim(keyword) == 'breit') then
				        Breit_type = .true.
				        write(60,"(a)") "   Breit term will be considered."
			        else if (trim(keyword) == 'srtp') then
				        SRTP_type = .true.
				        write(60,"(a)") "   second Relativized Thomas Precession will be considered."
	                else if (trim(keyword) == 'sttp') then
	                	STTP_type = .true.
	                	write(60,"(a)") "   spin Tensor Thomas Precession will be considered."
	                else if (trim(keyword) == 'mdcb') then
	                	mDCB_type = .true.
	                	write(60,"(a)") "   high speed modification of Dirac Coulomb Breit equation will be considered"
                    else if (trim(keyword) == 'finitenuc') then
	                	finitenuc = .true.
	                	write(60,"(a)") "   finite nuclear effect will be considered"
	                else
	                	call terminate('unknown keyword detected in module Hamiltonian')
	                end if
                end do module_Hamiltonian
            ! --------------------<module Radiation>--------------------
            else if (trim(module_name) == '%Radiation') then
            	module_Radiation: do
	                read(12,*,iostat = ios) keyword
	                if (ios /= 0 .or. keyword == '') then
	                	call terminate('empty line detected in module Radiation.')
                    else if (index(keyword,'!') == 1) then
                        cycle module_Radiation
	                else if (trim(keyword) == 'endRadiation') then
	                	exit module_Radiation
	                else if (trim(keyword) == 'end') then
	                	call terminate("end of module should be written as 'endmodulename'")
	                else
	                	call terminate('unknown keyword detected in module Radiation')
	                end if
                end do module_Radiation
            ! -----------------------<module SCF>-----------------------
            else if (trim(module_name) == '%SCF') then
            	module_SCF: do
	                read(12,*,iostat = ios) keyword
	                if (ios /= 0 .or. keyword == '') then
	                	call terminate('empty line detected in module SCF.')
                    else if (index(keyword,'!') == 1) then
                        cycle module_SCF
	                else if (trim(keyword) == 'endSCF') then
	                	exit module_SCF
	                else if (trim(keyword) == 'end') then
	                	call terminate("end of module should be written as 'endmodulename'")
                    else if (index(keyword,'maxiter') == 1) then
                        if (index(keyword,'=') /= 0) then
                			read(keyword(index(keyword,'=') + 1 : len(trim(keyword))),"(I)",iostat = ios) maxiter
                			if (ios /= 0) call terminate("maxiter setting should be written as 'maxiter=n'")
                		else
                			call terminate("maxiter setting should be written as 'maxiter=n'")
                        end if
                    else if (index(keyword,'schwarz') == 1) then
                        if (index(keyword,'=') /= 0) then
                			read(keyword(index(keyword,'=') + 1 : len(trim(keyword))),"(F)",iostat = ios) schwarz_VT
                			if (ios /= 0) call terminate("Schwarz screening setting should be written as 'schwarz=n', no scientific notation")
                		else
                			call terminate("Schwarz screening setting should be written as 'schwarz=n', no scientific notation")
                        end if
                    else if (index(keyword,'convertol') == 1) then
                        if (index(keyword,'=') /= 0) then
                			read(keyword(index(keyword,'=') + 1 : len(trim(keyword))),"(F)",iostat = ios) conver_tol
                			if (ios /= 0) call terminate("Convergence tolerence setting should be written as 'convertol=n', no scientific notation")
                		else
                			call terminate("Convergence tolerence setting should be written as 'convertol=n', no scientific notation")
                        end if
                    else if (index(keyword,'nodiis') == 1) then
                        if (index(keyword,'=') /= 0) then
                			read(keyword(index(keyword,'=') + 1 : len(trim(keyword))),"(I)",iostat = ios) nodiis
                			if (ios /= 0) call terminate("Initial iterations with no DIIS should be written as 'nodiis=n'")
                		else
                			call terminate("Initial iterations with no DIIS should be written as 'nodiis=n'")
                        end if
                    else if (index(keyword,'subsp') == 1) then
                        if (index(keyword,'=') /= 0) then
                			read(keyword(index(keyword,'=') + 1 : len(trim(keyword))),"(I)",iostat = ios) subsp
                			if (ios /= 0) call terminate("Dimension of suboptimal subspace should be written as 'subsp=n'")
                		else
                			call terminate("Dimension of suboptimal subspace should be written as 'subsp=n'")
                        end if
                    else if (index(keyword,'damp') == 1) then
                        if (index(keyword,'=') /= 0) then
                			read(keyword(index(keyword,'=') + 1 : len(trim(keyword))),"(F)",iostat = ios) damp
                			if (ios /= 0) call terminate("Mix parameter should be written as 'damp=n'")
                		else
                			call terminate("Mix parameter should be written as 'damp=n'")
                        end if
                    else if (index(keyword,'prtlev') == 1) then
                        if (index(keyword,'=') /= 0) then
                			read(keyword(index(keyword,'=') + 1 : len(trim(keyword))),"(F)",iostat = ios) prtlev
                			if (ios /= 0) call terminate("Print level setting should be written as 'prtlev=n'")
                		else
                			call terminate("Print level setting should be written as 'prtlev=n'")
                        end if
                    else if (index(keyword,'cutdiis') == 1) then
                        if (index(keyword,'=') /= 0) then
                			read(keyword(index(keyword,'=') + 1 : len(trim(keyword))),"(F)",iostat = ios) cutdiis
                			if (ios /= 0) call terminate("Cutting DIIS threshold should be written as 'cutdiis=n'")
                		else
                			call terminate("Cutting DIIS threshold should be written as 'cutdiis=n'")
                        end if
                    else if (trim(keyword) == 'keepspin') then
				        keepspin = .true.
				        write(60,"(a)") "   Avoid spin multiplicity mutations."
	                else
	                	call terminate('unknown keyword detected in module SCF')
	                end if
                end do module_SCF
            else if (trim(module_name) == '%') then
                call terminate("module name should be written as '%modulename'.")
            else
                call terminate('unknown module name detected')
            end if
        end do entire
        write(60,"(a)") "All calculation settings are default except for printing above, exit module Atoms"
        close(12)
        
!-----------------------------------------------------------------------
! judge the self-consistency of computational settings
        
        electron_count = electron_count - charge
        if (electron_count <= 0) then
            call terminate('charge is set incorrectly.')
        end if
        if (spin_mult == 675) then                                  ! default lowest spin
            if (mod(electron_count,2) == 0) then
                spin_mult = 1
            else
                spin_mult = 2
            end if
        else
            if (mod(electron_count,2) == 0) then
                if (mod(spin_mult,2) == 0 .or. spin_mult <= 0 .or. spin_mult > electron_count+1) then
                    call terminate('spin multiplicity is set incorrectly.')
                end if
            else
                if (mod(spin_mult,2) == 1 .or. spin_mult <= 0 .or. spin_mult > electron_count+1) then
                    call terminate('spin multiplicity is set incorrectly.')
                end if
            end if
        end if
        if (stacksize <= 500 .or. stacksize >= 100000) call terminate('stack size setting incorrectly, its magnitude is MB.')
        if (DKH_order == 0 .and. SRTP_type) call terminate('SRTP is not suitable for non-relativistic calculation')
        if (DKH_order == 0 .and. STTP_type) call terminate('STTP is not suitable for non-relativistic calculation')
        if (DKH_order == 0 .and. mDCB_type) call terminate('mDCB is not suitable for non-relativistic calculation')
        if (SRTP_type .and. STTP_type) call terminate('cannot set both SRTP and STTP.')
        if (SRTP_type .and. mDCB_type) call terminate('cannot set both SRTP and mDCB.')
        if (STTP_type .and. mDCB_type) call terminate('cannot set both STTP and mDCB.')
        if (subsp <= 1) call terminate('subsp two small')
        if (nodiis - subsp < 2) call terminate('nodiis should be set larger')
        if (damp > 1.0 .or. damp < 0.0) call terminate('damp shall be in range [0,1]')
    end subroutine read_keywords
    
    
end module Atoms