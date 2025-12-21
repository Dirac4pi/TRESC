!> @file Atoms.f90
!!
!! @brief construct polyatomic model
!!
!! @syntax Fortran 2008 free format
!!
!! @code UTF-8
!!
!! @author dirac4pi
module Atoms
  use Fundamentals
  
!-----------------------------------------------------------------------
! element definations & atom properties
  integer,parameter :: element_count = 30
  character(len = 2),parameter :: element_list(element_count) =       &
  (/ ' H',                                                      'He', &  ! 1-2
     'Li','Be',                        ' B',' C',' N',' O',' F','Ne', &  ! 3-10
     'Na','Mg',                        'Al','Si',' P',' S','Cl','Ar', &  ! 11-18
     ' K','Ca','Sc','Ti',' V','Cr','Mn','Fe','Co','Ni','Cu','Zn'     /)  ! 19-30
  
  integer,parameter :: element_massnumber(element_count) =            &
  (/  2,                                                           4, &  ! 1-2
      7,   9,                            11,  12,  14,  16,  19,  20, &  ! 3-10
     23,  24,                            27,  28,  31,  32,  35,  40, &  ! 11-18
     39,  40,  45,  48,  51,  52,  55,  56,  59,  59,  64,  65       /)  ! 19-30

  real(dp),parameter :: CSD_CovR(element_count) =                     &
  (/ 0.59,                                                      0.53, &  ! 1-2
     2.42,1.81,                        1.59,1.38,1.34,1.25,1.08,1.10, &  ! 3-10
     3.14,2.66,                        2.29,2.10,2.02,1.98,1.93,2.00, &  ! 11-18
     3.84,3.33,3.21,3.02,2.89,2.63,2.83,2.68,2.61,2.34,2.49,2.31     /)  ! 19-30

!-----------------------------------------------------------------------
! AO basis definations
integer, parameter :: AO_fac(3,5,15) = reshape( [                 &! (L,M)
!(x,y,z)   (x,y,z)   (x,y,z)   (x,y,z)   (x,y,z)
  0,0,0,    0,0,1,    2,0,0,    3,0,0,    4,0,0,                      &! M=1  <S
  0,0,0,    1,0,0,    0,2,0,    0,3,0,    0,4,0,                      &! M=2
  0,0,0,    0,1,0,    0,0,2,    0,0,3,    0,0,4,                      &! M=3  <P
  0,0,0,    0,0,0,    1,1,0,    1,2,0,    3,1,0,                      &! M=4
  0,0,0,    0,0,0,    1,0,1,    2,1,0,    3,0,1,                      &! M=5
  0,0,0,    0,0,0,    0,1,1,    2,0,1,    1,3,0,                      &! M=6  <D
  0,0,0,    0,0,0,    0,0,0,    1,0,2,    0,3,1,                      &! M=7
  0,0,0,    0,0,0,    0,0,0,    0,1,2,    1,0,3,                      &! M=8
  0,0,0,    0,0,0,    0,0,0,    0,2,1,    0,1,3,                      &! M=9
  0,0,0,    0,0,0,    0,0,0,    1,1,1,    2,2,0,                      &! M=10 <F
  0,0,0,    0,0,0,    0,0,0,    0,0,0,    2,0,2,                      &! M=11
  0,0,0,    0,0,0,    0,0,0,    0,0,0,    0,2,2,                      &! M=12
  0,0,0,    0,0,0,    0,0,0,    0,0,0,    2,1,1,                      &! M=13
  0,0,0,    0,0,0,    0,0,0,    0,0,0,    1,2,1,                      &! M=14
  0,0,0,    0,0,0,    0,0,0,    0,0,0,    1,1,2                       &! M=15 <G
  ], [3,5,15])

!-----------------------------------------------------------------------
! AO basis definations, sequence consistent with MOLDEN format
integer, parameter :: M_AO_fac(3,5,15) = reshape( [               &! (L,M)
!(x,y,z)   (x,y,z)   (x,y,z)   (x,y,z)   (x,y,z)
  0,0,0,    1,0,0,    2,0,0,    3,0,0,    4,0,0,                      &! M=1  <S
  0,0,0,    0,1,0,    0,2,0,    0,3,0,    0,4,0,                      &! M=2
  0,0,0,    0,0,1,    0,0,2,    0,0,3,    0,0,4,                      &! M=3  <P
  0,0,0,    0,0,0,    1,1,0,    1,2,0,    3,1,0,                      &! M=4
  0,0,0,    0,0,0,    1,0,1,    2,1,0,    3,0,1,                      &! M=5
  0,0,0,    0,0,0,    0,1,1,    2,0,1,    1,3,0,                      &! M=6  <D
  0,0,0,    0,0,0,    0,0,0,    1,0,2,    0,3,1,                      &! M=7
  0,0,0,    0,0,0,    0,0,0,    0,1,2,    1,0,3,                      &! M=8
  0,0,0,    0,0,0,    0,0,0,    0,2,1,    0,1,3,                      &! M=9
  0,0,0,    0,0,0,    0,0,0,    1,1,1,    2,2,0,                      &! M=10 <F
  0,0,0,    0,0,0,    0,0,0,    0,0,0,    2,0,2,                      &! M=11
  0,0,0,    0,0,0,    0,0,0,    0,0,0,    0,2,2,                      &! M=12
  0,0,0,    0,0,0,    0,0,0,    0,0,0,    2,1,1,                      &! M=13
  0,0,0,    0,0,0,    0,0,0,    0,0,0,    1,2,1,                      &! M=14
  0,0,0,    0,0,0,    0,0,0,    0,0,0,    1,1,2                       &! M=15 <G
  ], [3,5,15])

!-----------------------------------------------------------------------
! AO basis definations, sequence consistent with Gaussian09(.fch)
  integer, parameter :: g_AO_fac(3,5,15) = reshape( [                   &! (L,M)
!(x,y,z)   (x,y,z)   (x,y,z)   (x,y,z)   (x,y,z)
  0,0,0,    1,0,0,    2,0,0,    3,0,0,    0,0,4,                      &! M=1  <S
  0,0,0,    0,1,0,    0,2,0,    0,3,0,    0,1,3,                      &! M=2
  0,0,0,    0,0,1,    0,0,2,    0,0,3,    0,2,2,                      &! M=3  <P
  0,0,0,    0,0,0,    1,1,0,    1,2,0,    0,3,1,                      &! M=4
  0,0,0,    0,0,0,    1,0,1,    2,1,0,    0,4,0,                      &! M=5
  0,0,0,    0,0,0,    0,1,1,    2,0,1,    1,0,3,                      &! M=6  <D
  0,0,0,    0,0,0,    0,0,0,    1,0,2,    1,1,2,                      &! M=7
  0,0,0,    0,0,0,    0,0,0,    0,1,2,    1,2,1,                      &! M=8
  0,0,0,    0,0,0,    0,0,0,    0,2,1,    1,3,0,                      &! M=9
  0,0,0,    0,0,0,    0,0,0,    1,1,1,    2,0,2,                      &! M=10 <F
  0,0,0,    0,0,0,    0,0,0,    0,0,0,    2,1,1,                      &! M=11
  0,0,0,    0,0,0,    0,0,0,    0,0,0,    2,2,0,                      &! M=12
  0,0,0,    0,0,0,    0,0,0,    0,0,0,    3,0,1,                      &! M=13
  0,0,0,    0,0,0,    0,0,0,    0,0,0,    3,1,0,                      &! M=14
  0,0,0,    0,0,0,    0,0,0,    0,0,0,    4,0,0                       &! M=15 <G
  ], [3,5,15])

  ! transformation matrix from Cartesian basis to spher-harmo basis
  real(dp),parameter :: c2sd(6,5) = (/                           &! (6D,5D)
  -0.5_dp, -0.5_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                       &! D0
  0.0_dp,  0.0_dp,  0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp,                       &! D+1
  0.0_dp,  0.0_dp,  0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp,                       &! D-1
  dsqrt(3.0_dp)/2.0_dp, -dsqrt(3.0_dp)/2.0_dp,                            &
  0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                                         &! D+2
  0.0_dp,  0.0_dp,  0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp                       /)! D-2
  
  real(dp),parameter :: c2sf(10,7) = (/                          &! (10F,7F)
  0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, -1.5_dp/dsqrt(5.0_dp), 0.0_dp,  &
  0.0_dp, -1.5_dp/dsqrt(5.0_dp), 0.0_dp,                                  &! F0
  -dsqrt(3.0_dp/8.0_dp), 0.0_dp, 0.0_dp, -dsqrt(3.0_dp/40.0_dp),          &
  0.0_dp, 0.0_dp, dsqrt(6.0_dp/5.0_dp), 0.0_dp, 0.0_dp, 0.0_dp,           &! F+1
  0.0_dp, -dsqrt(3.0_dp/8.0_dp), 0.0_dp, 0.0_dp, -dsqrt(3.0_dp/40.0_dp),  &
  0.0_dp, 0.0_dp, dsqrt(6.0_dp/5.0_dp), 0.0_dp, 0.0_dp,                   &! F-1
  0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, dsqrt(3.0_dp)/2.0_dp, 0.0_dp,   &
  0.0_dp, -dsqrt(3.0_dp)/2.0_dp, 0.0_dp,                                  &! F+2
  0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,         &
  0.0_dp, 1.0_dp,                                                         &! F-2
  dsqrt(5.0_dp/8.0_dp), 0.0_dp, 0.0_dp, -3.0_dp/dsqrt(8.0_dp), 0.0_dp,    &
  0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                                 &! F+3
  0.0_dp, -dsqrt(5.0_dp/8.0_dp), 0.0_dp, 0.0_dp, 3.0_dp/dsqrt(8.0_dp),    &
  0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp                                 /)! F-3
  
  real(dp),parameter :: c2sg(15,9) = (/                          &! (15G,9G)
  0.375_dp, 0.375_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,     &
  0.0_dp, 3.0_dp/4.0_dp*dsqrt(3.0_dp/35.0_dp),                            &
  -3.0_dp*dsqrt(3.0_dp/35.0_dp), -3.0_dp*dsqrt(3.0_dp/35.0_dp), 0.0_dp,   &
  0.0_dp, 0.0_dp,                                                         &! G0
  0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -1.5_dp*dsqrt(5.0_dp/14.0_dp), 0.0_dp,  &
  0.0_dp, 2.0_dp*dsqrt(5.0_dp/14.0_dp), 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,   &
  0.0_dp, -1.5_dp/dsqrt(14.0_dp), 0.0_dp,                                 &! G+1
  0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                         &
  -1.5_dp*dsqrt(5.0_dp/14.0_dp), 0.0_dp, 2.0_dp*dsqrt(5.0_dp/14.0_dp),    &
  0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -1.5_dp/dsqrt(14.0_dp), 0.0_dp,         &! G-1
  -dsqrt(5.0_dp)/4.0_dp, dsqrt(5.0_dp)/4.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,    &
  0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 3.0_dp*dsqrt(3.0_dp/28.0_dp),   &
  -3.0_dp*dsqrt(3.0_dp/28.0_dp), 0.0_dp, 0.0_dp, 0.0_dp,                  &! G+2
  0.0_dp, 0.0_dp, 0.0_dp, -dsqrt(5.0_dp/28.0_dp), 0.0_dp,                 &
  -dsqrt(5.0_dp/28.0_dp), 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
  0.0_dp, 0.0_dp, 3.0_dp/dsqrt(7.0_dp),                                   &! G-2
  0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, dsqrt(5.0_dp/8.0_dp), 0.0_dp, 0.0_dp,   &
  0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -3.0_dp/dsqrt(8.0_dp),  &
  0.0_dp,                                                                 &! G+3
  0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -dsqrt(5.0_dp/8.0_dp),  &
  0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 3.0_dp/dsqrt(8.0_dp), 0.0_dp,   &
  0.0_dp,                                                                 &! G-3
  dsqrt(35.0_dp)/8.0_dp, dsqrt(35.0_dp)/8.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,   &
  0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
  0.0_dp,                                                                 &! G+4
  0.0_dp, 0.0_dp, 0.0_dp, dsqrt(5.0_dp)/2.0_dp, 0.0_dp,                   &
  -dsqrt(5.0_dp)/2.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,  &
  0.0_dp, 0.0_dp, 0.0_dp                                                 /)! G-4

  real(dp),parameter :: g_c2sg(15,9) = (/                        &! (15G,9G)
  1.0_dp, 0.0_dp, -3.0_dp*dsqrt(3.0_dp/35.0_dp), 0.0_dp, 0.375_dp, 0.0_dp,&
  0.0_dp, 0.0_dp, 0.0_dp, -3.0_dp*dsqrt(3.0_dp/35.0_dp), 0.0_dp,          &
  3.0_dp/4.0_dp*dsqrt(3.0_dp/35.0_dp), 0.0_dp, 0.0_dp, 0.375_dp,          &! G0
  0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 2.0_dp*dsqrt(5.0_dp/14.0_dp),   &
  0.0_dp, -1.5_dp/dsqrt(14.0_dp), 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,         &
  -1.5_dp*dsqrt(5.0_dp/14.0_dp), 0.0_dp, 0.0_dp,                          &! G+1
  0.0_dp, 2.0_dp*dsqrt(5.0_dp/14.0_dp), 0.0_dp,                           &
  -1.5_dp*dsqrt(5.0_dp/14.0_dp), 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,  &
  0.0_dp, -1.5_dp/dsqrt(14.0_dp), 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,         &! G-1
  0.0_dp, 0.0_dp, -3.0_dp*dsqrt(3.0_dp/28.0_dp), 0.0_dp,                  &
  dsqrt(5.0_dp)/4.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                   &
  3.0_dp*dsqrt(3.0_dp/28.0_dp), 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,           &
  -dsqrt(5.0_dp)/4.0_dp,                                                  &! G+2
  0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 3.0_dp/dsqrt(7.0_dp),   &
  0.0_dp, -dsqrt(5.0_dp/28.0_dp), 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,         &
  -dsqrt(5.0_dp/28.0_dp), 0.0_dp,                                         &! G-2
  0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                 &
  -3.0_dp/dsqrt(8.0_dp), 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                  &
  dsqrt(5.0_dp/8.0_dp), 0.0_dp, 0.0_dp,                                   &! G+3
  0.0_dp, 0.0_dp, 0.0_dp, -dsqrt(5.0_dp/8.0_dp), 0.0_dp, 0.0_dp, 0.0_dp,  &
  0.0_dp, 0.0_dp, 0.0_dp, 3.0_dp/dsqrt(8.0_dp), 0.0_dp, 0.0_dp, 0.0_dp,   &
  0.0_dp,                                                                 &! G-3
  0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, dsqrt(35.0_dp)/8.0_dp, 0.0_dp, 0.0_dp,  &
  0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -3.0_dp/4.0_dp*dsqrt(3.0_dp), 0.0_dp,   &
  0.0_dp, dsqrt(35.0_dp)/8.0_dp,                                          &! G+4
  0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,         &
  -dsqrt(5.0_dp)/2.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                  &
  dsqrt(5.0_dp)/2.0_dp, 0.0_dp                                           /)! G-4

  !DIR$ ATTRIBUTES ALIGN:align_size :: c2s, exc2s, s2c, exs2c, f2s, exf2s
  !DIR$ ATTRIBUTES ALIGN:align_size :: M_c2s, M_exc2s, s2f, exs2f, c2f, exc2f
  real(dp),allocatable    :: c2s(:,:)  ! convert Cartesian to spher-harmo
  complex(dp),allocatable :: exc2s(:,:)
  real(dp),allocatable    :: M_c2s(:,:)  ! convert Cartesian to spher-harmo
  complex(dp),allocatable :: M_exc2s(:,:)
  real(dp),allocatable    :: s2f(:,:)    ! orthogonal transform matrix
  complex(dp),allocatable :: exs2f(:,:)
  real(dp),allocatable    :: f2s(:,:)    ! reverse orthogonal transform matrix
  complex(dp),allocatable :: exf2s(:,:)
  real(dp),allocatable    :: c2f(:,:)    ! total transformation matrix
  complex(dp),allocatable :: exc2f(:,:)

  ! operator in basis transformation
  real(dp),allocatable    :: c2soper(:,:)
  complex(dp),allocatable :: exc2soper(:,:)
  real(dp),allocatable    :: s2coper(:,:)
  complex(dp),allocatable :: exs2coper(:,:)
  real(dp),allocatable    :: s2foper(:,:)
  complex(dp),allocatable :: exs2foper(:,:)
  real(dp),allocatable    :: f2soper(:,:)
  complex(dp),allocatable :: exf2soper(:,:)
  real(dp),allocatable    :: c2foper(:,:)
  complex(dp),allocatable :: exc2foper(:,:)
  real(dp),allocatable    :: sphoper(:,:)
  
  integer    :: cbdm                   ! Cartesian basis dimension
  integer    :: M_cbdm                 ! Cartesian basis dimension (MOLDEN)
  integer    :: sbdm                   ! spher-harmo basis dimension
  integer    :: M_sbdm                 ! spher-harmo basis dimension (MOLDEN)
  ! final basis dimension (after canonical orthogonalisation)
  integer    :: fbdm                   ! final basis dimension
  integer    :: basis_count            ! number of basis in basis set
  ! number of shells in each element
  integer    :: shell_in_element(50)

  ! The arrays in mol, atom_basis and cbdata will not involved in
  ! arithmetic, but are read frequently, so they are stored as static arrays
  ! and are not memory aligned.

  type atom_basis_type                 ! dimension basis_count
    integer  :: atom_number            ! atomic number of basis atom
    integer  :: L                      ! angular quantum number of shell
    integer  :: contr                  ! contr of basis shell
    real(dp) :: expo(16)               ! expo of primitive shell
    real(dp) :: coe(16)                ! contr coeff of prim shell
  end type atom_basis_type
  
  type(atom_basis_type) :: atom_basis(5000)   ! basis CGTOs
  type(atom_basis_type) :: M_atom_basis(5000) ! basis CGTOs from MOLDEN

  ! info of each Cartesian basis(directly related to atoms, dimension cbdm)
  type basis_data
    integer  :: atom                   ! which atom center in mol
    integer  :: shell                  ! shell number
    integer  :: L                      ! angular quantum number, S:1, P:2 ...
    integer  :: M                      ! magnetic quantum number
    real(dp) :: pos(3)                 ! position of each atom, (x,y,z)
    integer  :: contr                  ! contr of basis shell
    real(dp) :: expo(16)               ! expo of primitive shell
    real(dp) :: coe(16)                ! coeff of primitive shell
    real(dp) :: Ncoe(16)               ! normalized coeff of primitive shell
  end type basis_data
  type(basis_data), allocatable :: cbdata(:)  ! Cartesian basis data
  type(basis_data), allocatable :: sbdata(:)  ! Spherical-harmonic basis data
  type(basis_data), allocatable :: M_cbdata(:)! Cartesian basis data (MOLDEN)

!-----------------------------------------------------------------------
! molecular definations
  type mol_type
    integer  :: basis_number           ! serial number of each atom in basis set
    integer  :: atom_number            ! atomic number of each atom in mol
    real(dp) :: pos(3)                 ! position of each atom, (x,y,z)
    ! nucleus radius of each atom in mol 0.836*A^(1/3)+0.57 (fm)
    real(dp) :: rad
  end type mol_type
  integer        :: atom_count         ! number of atoms in mol
  type(mol_type) :: mol(100)           ! system composed only of atom
  integer        :: M_atom_count       ! number of atoms in M_mol
  type(mol_type) :: M_mol(100)         ! system composed only of atom (MOLDEN)
  type(mol_type) :: trafomol(100)      ! mol in frame of motion
  
  real(dp)   :: S__2                   ! <S**2> of mol
  real(dp)   :: totalpha, totbeta      ! total alpha, beta electron of mol
  real(dp)   :: S__2orb, Szorb         ! <S**2>, Sz of orbital
  real(dp)   :: totalphaorb, totbetaorb! total alpha, beta electron of orbital
  
  private :: Get_basis_count, Load_1MO_MOLDEN, Load_MOs_MOLDEN
  private :: csgo_1c, csgo_2c, sfgo_1c, sfgo_2c, cfgo_1c, cfgo_2c
  public  :: Load_basis_gbs, Load_geom_xyz, Load_geombasis_MOLDEN
  public  :: Load_keywords_tre, Inputcheck, Inputprint, Assign_csf, M_Assign_cs
  public  :: Load_MO_MOLDEN, csgo, sfgo, cfgo, M_csgo

  interface Load_MO_MOLDEN
    module procedure Load_1MO_MOLDEN
    module procedure Load_MOs_MOLDEN
  end interface Load_MO_MOLDEN
  ! basis transformation is irreversible, only involve:
  ! A_s = U^T A_c U
  !     = c2s^T A_c c2s
  ! C_c = U C_s
  !     = c2s A_s
  ! R_c = U R_s U^T
  !     = c2s R_s c2s^T
  ! where A is operator matrix, C is coefficient matrix, R is density matrix
  interface csgo                     ! operator matrix only
    module procedure csgo_1c
    module procedure csgo_2c
  end interface csgo
  interface scgo                     ! desity matrix only
    module procedure scgo_1c
    module procedure scgo_2c
  end interface scgo
  interface sfgo                     ! operator matrix only
    module procedure sfgo_1c
    module procedure sfgo_2c
  end interface sfgo
  interface fsgo                     ! desity matrix only
    module procedure fsgo_1c
    module procedure fsgo_2c
  end interface fsgo
  interface cfgo                     ! operator matrix only
    module procedure cfgo_1c
    module procedure cfgo_2c
  end interface cfgo

  contains
  
!-----------------------------------------------------------------------
!> get CGTOs from .gbs file, sharing of expo is not allowed
!!
!> default and only: spher-harmo segment contr basis
!!
!> minimal basis set are not available
  subroutine Load_basis_gbs()
    implicit none
    character(len=3) :: basis_element_name
    character(len=1) :: basis_angular_name
    integer          :: ai, aj, ak, am        ! loop variables Load_basis_gbs
    if(index(address_basis,'.gbs') == 0) call terminate(&
    'input basis set file is not .gbs file')
    inquire(file=address_basis, exist=exists, name=address_basis)
    if (.not. exists) call terminate(&
    'basis file '//address_basis//' do not exist')
    open(11, file=address_basis, status="old", action="read", iostat=ios)
    if (ios /= 0) call terminate(&
    'basis set file '//address_basis//' could not be read')
    call Get_basis_count()

    do
      read(11,*) basis_element_name
      if (index(basis_element_name,'!')==0 .and. basis_element_name/='') exit
    end do
    ai = 1
    aj = 1
    am = 0
    do
      if (adjustl(basis_element_name(index(basis_element_name,'-')+1:&
      index(basis_element_name,'-')+2)) /= adjustl(element_list(ai))) then
        ai = ai + 1
        if (ai >= element_count + 1) call terminate(&
        'basis file: element not included in current program')
        cycle
      end if
      atom_basis(aj) % atom_number = ai
      do
        read(11,*,iostat = ios) basis_angular_name, &
        atom_basis(aj) % contr  ! scale factor defaults to 1.00
        if (ios /= 0) call terminate(&
        'read basis file failed, may caused by incorrect formatting')
        if (basis_angular_name == 'S') then
          atom_basis(aj) % L = 0
        else if (basis_angular_name == 'P') then
          atom_basis(aj) % L = 1
        else if (basis_angular_name == 'D') then
          atom_basis(aj) % L = 2
        else if (basis_angular_name == 'F') then
          atom_basis(aj) % L = 3
        else if (basis_angular_name == 'G') then
          atom_basis(aj) % L = 4
        end if
        atom_basis(aj) % atom_number = ai
        ak = 1
        do
          if (ak <= atom_basis(aj) % contr) then
            read(11,*,iostat = ios) atom_basis(aj) % &
            expo(ak), atom_basis(aj) % coe(ak)
            if (ios /= 0) call terminate(&
            'read basis file failed, try denote scientific notation with E')
            ak = ak + 1
            cycle
          end if
          aj = aj + 1
          exit
        end do
        am = am + 1
        if (am == shell_in_element(ai)) exit
      end do
      read(11,*)  ! skip the separator between elements
      am = 0
      ai = 1
      read(11,*,iostat = ios) basis_element_name
      if (ios < 0) exit
      if (ios > 0) call terminate('Error encountered while reading .gbs file.')
    end do
    close(11)
  end subroutine Load_basis_gbs

!-----------------------------------------------------------------------
!> load geometry and CGTOs from MOLDEN file
!!
!! default and only: spher-harmo segment contr basis
!!
!! tomol: whether the wave-function is projected to the current molecule (geom)
  subroutine Load_geombasis_MOLDEN(tomol)
    implicit none
    logical,intent(in) :: tomol
    character(len=200) :: line
    character(len=2)   :: e
    integer            :: a, b, ix, iy, iz
    real(dp)           :: f1, f2, f3
    integer            :: contr           ! contr of atom, shell
    real(dp)           :: expo(16)        ! expo of |AO>
    real(dp)           :: coe(16)         ! coe of |AO>
    integer            :: ai, aj , ak, al ! loop variables Load_geombasis_MOLDEN
    open(14, file=address_MOLDEN, status="old", action="read", iostat=ios)
    if (ios /= 0) call terminate(&
    'MOLDEN file '//address_MOLDEN//' could not be read')
    do
      read(14,'(A)') line
      if (index(line,'[Title]')/=0) exit
      if (index(line,'[GTO]')/=0) exit    ! sometimes no [Title] in MOLDEN file
    end do
    read(14,'(A)') line
    call lowercase(line)
    if (index(line,'orca') /= 0) then
      isorca = .true.
      write(60,'(A)') '  -- MOLDEN file is generated by ORCA program'
    else
      isorca = .false.
      write(60,'(A)') '  -- MOLDEN file is NOT generated by ORCA program'
    end if
    close(14)
    ! load geometry
    open(14, file=address_MOLDEN, status="old", action="read")
    do
      read(14,*) line
      if (index(line,'[Atoms]')/=0) exit
    end do
    ai = 1
    do
      read(14,*,iostat=ios) e, a, b, f1, f2, f3
      if (ios /= 0) exit
      M_mol(ai)%atom_number = b
      M_mol(ai)%pos(1) = f1
      M_mol(ai)%pos(2) = f2
      M_mol(ai)%pos(3) = f3
      ai = ai + 1
    end do
    M_atom_count = ai - 1
    close(14)
    ! load basis
    open(14, file=address_MOLDEN, status="old", action="read")
    do
      read(14,*) line
      if (index(line,'[GTO]') /= 0) exit
    end do
    ai = 1
    aj = 1
    read(14,*) line
    M_cbdm = 0
    M_sbdm = 0
    allocate(M_cbdata(3*cbdm))
    do
      read(14,*,iostat=ios) e, a
      if (ios /= 0) exit
      call lowercase(e)
      select case(trim(e))
        case('s')
          b = 1
        case('p')
          b = 2
        case('d')
          b = 3
        case('f')
          b = 4
        case('g')
          b = 5
        case default
          aj = aj + 1
          cycle
      end select
      M_cbdm = M_cbdm + (b*(b+1))/2
      M_sbdm = M_sbdm + 2*b-1
      do al = 1, a
        read(14,*) f1, f2
        M_cbdata(ai:ai+(b*(b+1))/2-1)%expo(al) = f1
        M_cbdata(ai:ai+(b*(b+1))/2-1)%coe(al)  = f2
      end do
      do ak = 1, (b*(b+1))/2
        M_cbdata(ai)%L     = b
        M_cbdata(ai)%M     = ak
        M_cbdata(ai)%contr = a
        M_cbdata(ai)%atom  = aj
        !!!!!! this is very important !!!!!!
        ! mol: the wave-function is projected to the current molecule (geometry)
        ! M_mol: keep the original wave-function
        if (tomol) then
          M_cbdata(ai)%pos = mol(aj)%pos
        else
          M_cbdata(ai)%pos = M_mol(aj)%pos
        end if
        ai = ai + 1
      end do
    end do
    close(14)
  end subroutine Load_geombasis_MOLDEN

!-----------------------------------------------------------------------
!> load all MO coeffs from .MOLDEN file
!!
!! default and only: spher-harmo segment contr basis
  subroutine Load_MOs_MOLDEN(moa, mob)
    implicit none
    real(dp),intent(out)  ::  moa(:,:)   ! alpha orbitals
    real(dp),intent(out)  ::  mob(:,:)   ! beta orbitals
    integer               ::  ii, jj, kk ! loop variables
    character(len=200)    ::  line
    integer               ::  a

    open(14, file=address_MOLDEN, status="old", action="read", iostat=ios)
    if (ios /= 0) call terminate(&
    'MOLDEN file '//address_MOLDEN//' could not be read')
    ii = 1
    jj = 1
    do
      read(14,'(A)',iostat=ios) line
      if (ios /= 0) exit
      if (index(line,'=')/=0 .and. index(line,'Alpha')/=0) then
        read(14,'(A)') line
        do kk = 1, M_sbdm
          read(14,*) a, moa(kk,ii)
        end do
        ii = ii + 1
      else if (index(line,'=')/=0 .and. index(line,'Beta')/=0) then
        read(14,'(A)') line
        do kk = 1, M_sbdm
          read(14,*) a, mob(kk,jj)
        end do
        jj = jj + 1
      end if
    end do
    close(14)
    ! some of high angular momentum ORCA MOs are normalized to -1 rather than 1
    ! so convert them
    if (isorca) then
      if (.not. allocated(M_c2s)) call M_Assign_cs()
      ii = 1
      jj = 1
      do while(ii <= M_cbdm)
        if (M_cbdata(ii)%L == 4) then
          moa(jj+5,:) = -moa(jj+5,:)
          mob(jj+5,:) = -mob(jj+5,:)
          moa(jj+6,:) = -moa(jj+6,:)
          mob(jj+6,:) = -mob(jj+6,:)
        else if (M_cbdata(ii)%L == 5) then
          moa(jj+5,:) = -moa(jj+5,:)
          mob(jj+5,:) = -mob(jj+5,:)
          moa(jj+6,:) = -moa(jj+6,:)
          mob(jj+6,:) = -mob(jj+6,:)
          moa(jj+7,:) = -moa(jj+7,:)
          mob(jj+7,:) = -mob(jj+7,:)
          moa(jj+8,:) = -moa(jj+8,:)
          mob(jj+8,:) = -mob(jj+8,:)
        end if
        ii = ii + (M_cbdata(ii)%L*(M_cbdata(ii)%L+1))/2
        jj = jj + 2*M_cbdata(ii)%L-1
      end do
    end if
  end subroutine Load_MOs_MOLDEN

!-----------------------------------------------------------------------
!> load MO coeffs from .MOLDEN file
!!
!! default and only: spher-harmo segment contr basis
!!
!! imo: which mo to be loaded, if imo<=0, will ignore it
  subroutine Load_1MO_MOLDEN(imoa, moa, imob, mob)
    implicit none
    integer,intent(in)    ::  imoa       ! which alpha mo to be loaded
    integer,intent(in)    ::  imob       ! which beta mo to be loaded
    real(dp),intent(out)  ::  moa(:)     ! alpha orbital
    real(dp),intent(out)  ::  mob(:)     ! beta orbital
    integer               ::  ii, jj, kk ! loop variables
    character(len=200)    ::  line
    integer               ::  a
    open(14, file=address_MOLDEN, status="old", action="read", iostat=ios)
    if (ios /= 0) call terminate(&
    'MOLDEN file '//address_MOLDEN//' could not be read')
    ii = 1
    jj = 1
    do
      read(14,'(A)',iostat=ios) line
      if (ios /= 0) exit
      if (index(line,'=')/=0 .and. index(line,'Alpha')/=0) then
        if (ii == imoa) then
          read(14,'(A)') line
          do kk = 1, M_sbdm
            read(14,*) a, moa(kk)
          end do
        else
          read(14,'(A)') line
          do kk = 1, M_sbdm
            read(14,'(A)') line
          end do
        end if
        ii = ii + 1
      else if (index(line,'=')/=0 .and. index(line,'Beta')/=0) then
        if (jj == imob) then
          read(14,'(A)') line
          do kk = 1, M_sbdm
            read(14,*) a, mob(kk)
          end do
        else
          read(14,'(A)') line
          do kk = 1, M_sbdm
            read(14,'(A)') line
          end do
        end if
        jj = jj + 1
      end if
      if (ii == M_sbdm+1 .and. jj == M_sbdm+1) exit
      if ((ii>imoa .or. imoa<=0) .and. (jj>imob .or. imob<=0)) exit
    end do
    close(14)
    if (imoa <= 0) moa = 0.0_dp
    if (imob <= 0) mob = 0.0_dp
    ! some of high angular momentum ORCA MOs are normalized to -1 rather than 1
    ! so convert them
    if (isorca) then
      if (.not. allocated(M_c2s)) call M_Assign_cs()
      ii = 1
      jj = 1
      do while(ii <= M_cbdm)
        if (M_cbdata(ii)%L == 4) then
          moa(jj+5) = -moa(jj+5)
          mob(jj+5) = -mob(jj+5)
          moa(jj+6) = -moa(jj+6)
          mob(jj+6) = -mob(jj+6)
        else if (M_cbdata(ii)%L == 5) then
          moa(jj+5) = -moa(jj+5)
          mob(jj+5) = -mob(jj+5)
          moa(jj+6) = -moa(jj+6)
          mob(jj+6) = -mob(jj+6)
          moa(jj+7) = -moa(jj+7)
          mob(jj+7) = -mob(jj+7)
          moa(jj+8) = -moa(jj+8)
          mob(jj+8) = -mob(jj+8)
        end if
        ii = ii + (M_cbdata(ii)%L*(M_cbdata(ii)%L+1))/2
        jj = jj + 2*M_cbdata(ii)%L-1
      end do
    end if
  end subroutine Load_1MO_MOLDEN

!-----------------------------------------------------------------------
!> get the number of basis functions in .gbs file
  subroutine Get_basis_count()
    implicit none
    character(len = 200) :: line
    integer              :: ai, aj, ak, al ! loop variables Get_basis_count
    ai = 0
    aj = 0
    open(13,file = address_basis,status = "old",action = "read")
    do
      read(13,'(A200)') line
      if (index(line,'!') == 0 .and. line /= '') exit
    end do
    do
      if (ios < 0) exit
      if (ios > 0) call terminate('Error encountered while reading .gbs file.')
      if (index(line,'S') /= 0 .or. &
      index(line,'P') /= 0 .or. index(line,'D') /= 0 &
      .or. index(line,'F') /= 0 .or. index(line,'G') /= 0) then
        if (index(line,'1.00') /= 0 .or. index(line,'1.0') /= 0) then
          ai = ai + 1
          read(13,'(A200)',iostat = ios) line
          cycle
        end if
      end if
      if (index(line,'0') /= 0 .and. index(line,'.') == 0) then
        do al = 1, element_count
          if (adjustl(line(index(line,'-')+1:&
          index(line,'-')+2)) == adjustl(element_list(al))) then
            ak = al
            exit
          end if
        end do
        if (al == element_count+1) then
          call terminate(&
          'basis file: element not included in current program')
        end if
      end if
      if (index(line,'***') /= 0) then
        shell_in_element(ak) = ai - aj
        aj = ai
      end if  
      read(13,'(A200)',iostat = ios) line
    end do
    close(13)
    basis_count = ai
  end subroutine Get_basis_count

!-----------------------------------------------------------------------
!> load static geometry of mol from .xyz file
!!
!! will also calculate cbdm, sbdm and allocate cbdata, sbdata
!!
!! must be called after Load_basis_gbs
  subroutine Load_geom_xyz()
    implicit none
    character(len = 2) :: mol_element_name, title_note
    integer            :: atom                 ! atom of |AO>
    integer            :: shell                ! shell of |AO>
    integer            :: shell_start          ! start point of an shell
    integer            :: L                    ! angular quantum number of |AO>
    integer            :: M                    ! magnetic quantum number of |AO>
    integer            :: ai, aj, ak           ! loop variables Load_geom_xyz
    if(index(address_xyz,'.xyz') == 0) address_xyz = trim(address_xyz)//'.xyz'
    inquire(file=address_xyz, exist=exists)
    if (.not. exists) &
    call terminate('geometry file '//address_xyz//' does not exist')
    open(12, file=address_xyz, status="old", action="read")
    read(12,*) title_note
    do
      read(title_note, "(I3)", iostat=ios) atom_count
      if (ios == 0) then
        exit
      else if (ios /= 0 .and. index(title_note,'!') == 1) then
        read(12,*) title_note
        cycle
      else
        call terminate("use '!' tag for comments in .xyz file")
      end if
    end do
    read(12,*)
    do aj = 1, atom_count, 1
      read(12,*,iostat=ios) mol_element_name, mol(aj)%pos(1:3)
      mol(aj)%pos(1:3) = mol(aj)%pos(1:3) / Ang2Bohr
      if (ios /= 0) call terminate(&
      'read geometry file failed, may caused by the absence of atom')
      ai = 1
      do while(adjustl(mol_element_name) /= adjustl(element_list(ai)))
        ai = ai + 1
        if (ai > element_count) call terminate(&
        'geometry file: element not contained in current program')
      end do
      ak = 1
      do while(atom_basis(ak) % atom_number /= ai)
        ak = ak + 1
        if (ak >= basis_count) call terminate(&
        'there are elements in geometry file not contained in basis set file')
      end do
      mol(aj) % basis_number = ak
      mol(aj) % atom_number = ai
      mol(aj) % rad = 0.836_dp*real(element_massnumber(ai))**(1.0/3.0) + 0.57_dp
    end do
    read(12,*,iostat=ios) mol_element_name, &
    mol(atom_count+1)%pos(1:3)
    if (ios == 0) call terminate(&
    'read geometry file failed, may caused by redundancy of atoms')
    close(12)
    electron_count = sum(mol(1:atom_count) % atom_number)

    ! get the dimension of Cartesian basis
    cbdm = 0
    sbdm = 0
    do ai = 1, atom_count
      do aj = 1, shell_in_element(mol(ai) % atom_number)
        cbdm = cbdm + &
        (atom_basis(mol(ai)%basis_number+aj-1)%L+2) * &
        (atom_basis(mol(ai)%basis_number+aj-1)%L+1) / 2
        sbdm = sbdm + 2 * (atom_basis(mol(ai)%basis_number+aj-1)%L) + 1
      end do
    end do
    ! generate cbdata
    allocate(cbdata(cbdm))
    ai = 1
    atom = 1
    shell = 1
    shell_start = 1
    do while(ai <= cbdm)
      if (shell > shell_in_element(mol(atom) % atom_number)) then
        shell = 1
        atom = atom + 1
      end if
      L = atom_basis(mol(atom)%basis_number+shell-1)%L + 1
      M = ai - shell_start + 1
      cbdata(ai)%atom  = atom
      cbdata(ai)%shell = shell
      cbdata(ai)%L     = L
      cbdata(ai)%M     = M
      cbdata(ai)%pos   = mol(atom)%pos
      cbdata(ai)%contr = atom_basis(mol(atom)%basis_number+shell-1)%contr
      cbdata(ai)%expo  = atom_basis(mol(atom)%basis_number+shell-1)%expo
      cbdata(ai)%coe   = atom_basis(mol(atom)%basis_number+shell-1)%coe
      ai = ai + 1
      if (ai-shell_start >= (L+1)*L/2) then
        shell = shell + 1
        shell_start = ai
      end if
    end do
    ! generate sbdata
    allocate(sbdata(sbdm))
    ai = 1
    atom = 1
    shell = 1
    shell_start = 1
    do while(ai <= sbdm)
      if (shell > shell_in_element(mol(atom) % atom_number)) then
        shell = 1
        atom = atom + 1
      end if
      L = atom_basis(mol(atom)%basis_number+shell-1)%L + 1
      M = nint(real(int((ai-shell_start+1)/2)) * (-1.0)**(ai-shell_start+1))
      sbdata(ai)%atom  = atom
      sbdata(ai)%shell = shell
      sbdata(ai)%L     = L
      sbdata(ai)%M     = M
      sbdata(ai)%pos   = mol(atom)%pos
      ai = ai + 1
      if (ai-shell_start >= 2*L-1) then
        shell = shell + 1
        shell_start = ai
      end if
    end do
  end subroutine Load_geom_xyz

!-----------------------------------------------------------------------
!> print basis set information, geometry information and calculation settings
  subroutine Inputprint()
    implicit none
    integer :: ai, aj, ak, am  ! loop variables Inputprint
    write(60,"(A)") 'Module Atoms:'
    write(60,"(A)") '  input basis set file path: '
    write(60,"(A)") '  '//trim(address_basis)
    write(60,"(A)") '  ----------<BASIS>----------'
    ai = 1
    do while(ai <= basis_count)  ! print only basis of atoms in mol
      do aj = 1, atom_count
        if (mol(aj) % basis_number == ai) exit
      end do
      if (aj == atom_count + 1) then
        ai = ai + 1
        cycle
      end if
      write(60,"(A10,A2,A)") &
      '  element ',element_list(atom_basis(ai) % atom_number),':'
      ak = 0
      do while(ak <= shell_in_element(mol(aj) % atom_number) - 1)
        write(60,"(A15,I1,A1)") '  -- shell l = ',&
        atom_basis(ai + ak) % L,':'
        am = 1
        do while(am <= atom_basis(ai + ak) % contr)
          write(60,"(A7,E12.5E2,A6,F9.5)") &
          '  exp: ',atom_basis(ai + ak) % &
          expo(am), ', coe:',atom_basis(ai + ak) % coe(am)
          am = am + 1
        end do
        ak = ak + 1
      end do
      ai = ai + 1
      write(60,*)
    end do
    write(60,*)
    write(60,"(A)") '  input geometry file path: '
    if (index(address_xyz,'/') == 0) then
      write(60,"(A)") '  WD/'//trim(address_xyz)
    else
      write(60,"(A)") '  '//trim(address_xyz)
    end if
    write(60,"(A)") '  ----------<GEOMETRY>----------'
    do ai = 1, atom_count
      write(60,"(A2, A2, A3, F9.5, A3, F9.5, A3, F9.5)") &
      '  ',element_list(mol(ai) % atom_number), '   ',&
      mol(ai) % pos(1), '   ', mol(ai) % pos(2), '   ', mol(ai) % pos(3)
    end do
    write(60,"(A22,I4,A2,I4)") '  scalar/spinor cbdm: ', cbdm, ' /', 2*cbdm
    write(60,"(A22,I4,A2,I4)") '  scalar/spinor sbdm: ', sbdm, ' /', 2*sbdm
  end subroutine Inputprint
  
!-----------------------------------------------------------------------
!> read keywords from .tre file
  subroutine Load_keywords_tre()
    implicit none
    integer              :: aj
    character(len = 40)  :: module_name, keyword, pathword
    write(60,"(A)") "  ----------<KEYWORDS>----------"
    open(12, file=address_tre, status="old", action="read")
    entire: do
      read(12,*,iostat = ios) module_name
      call lowercase(module_name)
      interval: do  ! Allow multiple blank lines between each module
        if (ios /= 0) then
          exit entire
        else if (module_name == '' .or. index(module_name,'!') == 1) then
          read(12,*,iostat = ios) module_name
          call lowercase(module_name)
          cycle interval
        end if
        exit interval
      end do interval
      ! ----------------------<module Atoms>----------------------
      if (trim(module_name) == '%atoms') then
        module_Atoms: do
          read(12,*,iostat = ios) keyword
          pathword = keyword
          call lowercase(keyword)
          if (ios /= 0 .or. keyword == '') then
            call terminate('empty line detected in module Atoms.')
          else if (index(keyword,'!') == 1) then
            cycle module_Atoms
          else if (trim(keyword) == 'endatoms') then
            exit module_Atoms
          else if (trim(keyword) == 'end') then
            call terminate("end of module should be written as 'endmodulename'")
          else if (index(keyword,'charge') == 1) then
            if (index(keyword,'=') /= 0) then
              read(keyword(index(keyword,'=') + 1 : len(trim(keyword))),&
              "(I2)",iostat = ios) charge
              if (ios /= 0) call terminate(&
              "charge setting should be written as 'charge=n'")
            else
              call terminate("charge setting should be written as 'charge=n'")
            end if
            write(60,"(A,I2)") "  charge set to ",charge
          else if (index(keyword,'spin') == 1) then
            if (index(keyword,'=') /= 0) then
              read(keyword(index(keyword,'=') + 1 : len(trim(keyword))),&
              "(I2)",iostat = ios) spin_mult
              if (ios /= 0) call terminate(&
              "spin multiplicity setting should be written as 'spin=n'")
            else
              call terminate(&
              "spin multiplicity setting should be written as 'spin=n'")
            end if
            write(60,"(A,I2)") "  spin multiplicity set to ",spin_mult
          else if (index(keyword,'basis') == 1) then
            if (index(keyword,'=') /= 0) then
              call getenv('TRESC',address_basis)
              if (address_basis == '') call terminate (&
              "Environment variable 'TRESC' could not be found")
              if (index(keyword,'.gbs') /= 0) then
                address_basis = trim(address_basis) // '/basis/'&
                // pathword(index(pathword,'=')+1 : len(trim(pathword)))
              else
                address_basis = trim(address_basis) // '/basis/'&
                // pathword(index(pathword,'=')+1 : len(trim(pathword)))//'.gbs'
              end if
            else
              call terminate(&
              "basis setting should be written as 'basis=***.gbs'")
            end if
          else if (index(keyword,'geom') == 1) then
            if (index(keyword,'=') /= 0) then
              address_xyz = pathword(index(pathword,'=')+1:len(trim(pathword)))
            else
              call terminate(&
              "geometry setting should be written as 'geom=***.xyz'")
            end if
          else
            call terminate('unknown keyword detected in module Atoms')
          end if
        end do module_Atoms 
      ! -------------------<module Hamiltonian>-------------------
      else if (trim(module_name) == '%hamiltonian') then
        module_Hamiltonian: do
          read(12,*,iostat = ios) keyword
          call lowercase(keyword)
          if (ios /= 0 .or. keyword == '') then
            call terminate('empty line detected in module Hamiltonian.')
          else if (index(keyword,'!') == 1) then
            cycle module_Hamiltonian
          else if (trim(keyword) == 'endhamiltonian') then
            exit module_Hamiltonian
          else if (trim(keyword) == 'end') then
            call terminate("end of module should be written as 'endmodulename'")
          else if (index(keyword,'threads') == 1) then
            if (index(keyword,'=') /= 0) then
              read(keyword(index(keyword,'=') + 1 : len(trim(keyword))),&
              "(I2)",iostat = ios) threads
              if (ios /= 0) call terminate(&
              "threads setting should be written as 'threads=n'")
              write(60,'(A2,I3,A)') &
              '  ',threads,' threads will be used'
            else
              call terminate("threads setting should be written as 'threads=n'")
            end if
          else if (trim(keyword) == 'pvp1e') then
            pVp1e = .true.
            write(60,"(A)") "  one-electron pVp potential will be considered"
          else if (trim(keyword) == 'pvp2e') then
            pVp2e = .true.
            write(60,"(A)") "  two-electron pVp potential will be considered"
          else if (trim(keyword) == 'pvp') then
            pVp1e = .true.
            pVp2e = .true.
            write(60,"(A)") "  one-electron pVp potential will be considered"
            write(60,"(A)") "  two-electron pVp potential will be considered"
          else if (trim(keyword) == 'pppvp') then
            pppVp = .true.
            write(60,"(A)") &
            "  Second Relativized Thomas Precession will be considered"
          else if (index(keyword,'cuts') == 1) then
            if (index(keyword,'=') /= 0) then
              read(keyword(index(keyword,'=') + 1 : len(trim(keyword))),&
              "(F20.12)",iostat = ios) cutS
              if (ios /= 0) call terminate(&
              "cutS should be written as 'cutS=n', no scientific notation")
              write(60,"(A,E10.3)") &
              "  linear dependence threshold changed to ", cutS
            else
              call terminate(&
              "cutS should be written as 'cutS=n', no scientific notation")
            end if
          else
            call terminate('unknown keyword detected in module Hamiltonian')
          end if
        end do module_Hamiltonian
      ! --------------------<module Representation>--------------------
      else if (trim(module_name) == '%representation') then
        module_Representation: do
          read(12,*,iostat = ios) keyword
          call lowercase(keyword)
          if (ios /= 0 .or. keyword == '') then
            call terminate('empty line detected in module Representation.')
          else if (index(keyword,'!') == 1) then
            cycle module_Representation
          else if (trim(keyword) == 'endrepresentation') then
            exit module_Representation
          else if (trim(keyword) == 'end') then
            call terminate("end of module should be written as 'endmodulename'")
          else if (index(keyword,'betax=') == 1) then
            if (index(keyword,'=') /= 0) then
              read(keyword(index(keyword,'=')+1 : len(trim(keyword))),&
              "(E12.5)",iostat = ios) beta(1)
              if (ios /= 0) call terminate(&
              "v_x/c should be written as 'betax=n'")
            else
              call terminate("v_x/c should be written as 'betax=n'")
            end if
          else if (index(keyword,'betay=') == 1) then
            if (index(keyword,'=') /= 0) then
              read(keyword(index(keyword,'=')+1 : len(trim(keyword))),&
              "(E12.5)",iostat = ios) beta(2)
              if (ios /= 0) call terminate(&
              "v_y/c should be written as 'betay=n'")
            else
              call terminate("v_y/c should be written as 'betay=n'")
            end if
          else if (index(keyword,'betaz=') == 1) then
            if (index(keyword,'=') /= 0) then
              read(keyword(index(keyword,'=')+1 : len(trim(keyword))),&
              "(E12.5)",iostat = ios) beta(3)
              if (ios /= 0) call terminate(&
              "v_z/c should be written as 'betaz=n'")
            else
              call terminate("v_z/c should be written as 'betaz=n'")
            end if
          else
            call terminate('unknown keyword detected in module Representation')
          end if
        end do module_Representation
      ! -----------------------<module SCF>-----------------------
      else if (trim(module_name) == '%scf') then
        module_SCF: do
          read(12,*,iostat = ios) keyword
          call lowercase(keyword)
          if (ios /= 0 .or. keyword == '') then
            call terminate('empty line detected in module SCF.')
          else if (index(keyword,'!') == 1) then
            cycle module_SCF
          else if (trim(keyword) == 'endscf') then
            exit module_SCF
          else if (trim(keyword) == 'end') then
            call terminate("end of module should be written as 'endmodulename'")
          else if (index(keyword,'maxiter') == 1) then
            if (index(keyword,'=') /= 0) then
              read(keyword(index(keyword,'=') + 1 : len(trim(keyword))),&
              "(I3)",iostat = ios) maxiter
              if (ios /= 0) call terminate(&
              "maxiter setting should be written as 'maxiter=n'")
              write(60,'(A,I4)') '  SCF max iteration changed to ',maxiter
            else
              call terminate("maxiter setting should be written as 'maxiter=n'")
            end if
          else if (index(keyword,'schwarz') == 1) then
            if (index(keyword,'=') /= 0) then
              read(keyword(index(keyword,'=') + 1 : len(trim(keyword))),&
              "(F20.12)",iostat = ios) schwarz_VT
              if (ios /= 0) call terminate(&
              "Schwarz screening setting should be written as 'schwarz=n',"//&
              " no scientific notation")
              write(60,'(A,E10.3)') &
              '  Schwarz screen threshold changed to ',schwarz_VT
            else
              call terminate(&
              "Schwarz screening setting should be written as 'schwarz=n',"//&
              " no scientific notation")
            end if
          else if (index(keyword,'convertol') == 1) then
            if (index(keyword,'=') /= 0) then
              read(keyword(index(keyword,'=') + 1 : len(trim(keyword))),&
              "(F20.12)",iostat = ios) conver_tol
              if (ios /= 0) call terminate(&
              "Convergence tolerence setting should be written as "//&
              "'convertol=n', no scientific notation")
              write(60,'(A,E10.3)') &
              '  convergence tolerance changed to ',conver_tol
            else
              call terminate(&
              "Convergence tolerence setting should be written as "//&
              "'convertol=n', no scientific notation")
            end if
          else if (index(keyword,'nodiis') == 1) then
            if (index(keyword,'=') /= 0) then
              read(keyword(index(keyword,'=') + 1 : len(trim(keyword))),&
              "(I3)",iostat = ios) nodiis
              if (ios /= 0) call terminate(&
              "Initial iterations with no DIIS should be written as 'nodiis=n'")
              write(60,'(A,I3)') '  nodiis changed to ',nodiis
            else
              call terminate(&
              "Initial iterations with no DIIS should be written as 'nodiis=n'")
            end if
          else if (index(keyword,'subsp') == 1) then
            if (index(keyword,'=') /= 0) then
              read(keyword(index(keyword,'=') + 1 : len(trim(keyword))),&
              "(I3)",iostat = ios) subsp
              if (ios /= 0) call terminate(&
              "Dimension of suboptimal subspace should be written as 'subsp=n'")
              write(60,'(A,I3)') '  subsp changed to ',subsp
            else
              call terminate(&
              "Dimension of suboptimal subspace should be written as 'subsp=n'")
            end if
          else if (index(keyword,'damp') == 1) then
            if (index(keyword,'=') /= 0) then
              read(keyword(index(keyword,'=') + 1 : len(trim(keyword))),&
              "(F20.12)",iostat = ios) damp
              if (ios /= 0) call terminate(&
              "Mix parameter should be written as 'damp=n'")
              write(60,'(A,E10.3)') '  damp changed to ',damp
            else
              call terminate("Mix parameter should be written as 'diisdamp=n'")
            end if
          else if (index(keyword,'diisdamp') == 1) then
            if (index(keyword,'=') /= 0) then
              read(keyword(index(keyword,'=') + 1 : len(trim(keyword))),&
              "(F20.12)",iostat = ios) diisdamp
              if (ios /= 0) call terminate(&
              "DIIS mix parameter should be written as 'diisdamp=n'")
              write(60,'(A,E10.3)') '  diisdamp changed to ',diisdamp
            else
              call terminate("Mix parameter should be written as 'diisdamp=n'")
            end if
          else if (index(keyword,'prtlev') == 1) then
            if (index(keyword,'=') /= 0) then
              read(keyword(index(keyword,'=') + 1 : len(trim(keyword))),&
              "(F20.12)",iostat = ios) prtlev
              if (ios /= 0) call terminate(&
              "Print level setting should be written as 'prtlev=n'")
              write(60,'(A,E10.3)') '  print level changed to ',prtlev
            else
              call terminate(&
              "Print level setting should be written as 'prtlev=n'")
            end if
          else if (index(keyword,'molden') == 1) then
            MOLDEN = .true.
            write(60,'(A)') '  canonical orbitals will be dumped to .molden.d'
          else if (index(keyword,'cutdiis') == 1) then
            if (index(keyword,'=') /= 0) then
              read(keyword(index(keyword,'=') + 1 : len(trim(keyword))),&
              "(F20.12)",iostat = ios) cutdiis
              if (ios /= 0) call terminate(&
              "Cutting DIIS threshold should be written as 'cutdiis=n'")
            else
              call terminate(&
              "Cutting DIIS threshold should be written as 'cutdiis=n'")
            end if
          else if (index(keyword,'cutdamp') == 1) then
            if (index(keyword,'=') /= 0) then
              read(keyword(index(keyword,'=') + 1 : len(trim(keyword))),&
              "(F20.12)",iostat = ios) cutdamp
              if (ios /= 0) call terminate(&
              "Cutting damp threshold should be written as 'cutdamp=n'")
              write(60,'(A,E10.3)') '  cut damp threshold changed to ',cutdamp
            else
              call terminate&
              ("Cutting damp threshold should be written as 'cutdamp=n'")
            end if
          else if (index(keyword,'cspin') == 1) then
            if (index(keyword,'=') /= 0) then
              read(keyword(index(keyword,'=')+1 : index(keyword,'=')+1),&
              "(A1)",iostat = ios) cspin
              if (ios /= 0) call terminate(&
            "Constrained spin multiplicity mode should be written as 'cspin=*'")
            else
              call terminate(&
            "Constrained spin multiplicity mode should be written as 'cspin=*'")
            end if
            write(60,"(A)")"  constrained spin multiplicity mode set to "//cspin
            if (cspin == 'f' .or. cspin == 'd') then
              write(60,*)
              write(60,'(A)') "  WARNNING:"
              write(60,'(A)') "  -- when cspin = f/d, the selective occupation"
              write(60,'(A)') "  -- of the MOs will likely lead to a violation"
              write(60,'(A)') "  -- of the Aufbau principle, further leading to"
              write(60,'(A)') "  -- variational instability! Make sure to check"
              write(60,'(A)') "  -- the orbital occupation during SCF to"
              write(60,'(A)') "  -- determine if the convergence are reliable!"
              write(60,*)
              write(*,'(A)') "  WARNNING:"
              write(*,'(A)') "  -- when cspin = f/d, the selective occupation"
              write(*,'(A)') "  -- of the MOs will likely lead to a violation"
              write(*,'(A)') "  -- of the Aufbau principle, further leading to"
              write(*,'(A)') "  -- variational instability! Make sure to check"
              write(*,'(A)') "  -- the orbital occupation during SCF to"
              write(*,'(A)') "  -- determine if the convergence are reliable!"
            end if
          else if (index(keyword,'emd4') == 1) then
            d4 = .true.
            if (index(keyword,'=') /= 0) then
              read(keyword(index(keyword,'=') + 1 : len(trim(keyword))),&
              "(A)",iostat = ios) funcemd4
              if (ios /= 0) call terminate(&
              "dispersion correction should be written as 'emd4=func_name'")
            else
              call terminate(&
              "dispersion correction should be written as 'emd4=func_name'")
            end if
            write(60,"(A)") "  DFT-D4 empirical dispersion will be considered"
            write(60,"(A)") "  -- functional name: "//trim(funcemd4)
          else if (index(keyword,'guess') == 1) then
            if (index(keyword,'=') /= 0) then
              guess_type = keyword(index(keyword,'=') + 1 : len(trim(keyword)))
              if (ios /= 0) call terminate(&
              "Initial guess setting should be written as 'guess=***'")
              if (guess_type == 'ao2mo') then
                write(60,"(A)") "  initial guess load from .ao2mo file"
              else if (guess_type == 'molden') then
                write(60,"(A)") "  initial guess load from MOLDEN file"
              end if
            else
              call terminate(&
              "Initial guess setting should be written as 'guess=***'")
            end if
          else
            call terminate('unknown keyword detected in module SCF')
          end if
        end do module_SCF
      else if (trim(module_name) == '%functional') then
        module_functional: do
          read(12,*,iostat = ios) keyword
          call lowercase(keyword)
          if (ios /= 0 .or. keyword == '') then
            call terminate('empty line detected in module Functional.')
          else if (index(keyword,'!') == 1) then
            cycle module_functional
          else if (trim(keyword) == 'endfunctional') then
            exit module_functional
          else if (trim(keyword) == 'end') then
            call terminate("end of module should be written as 'endmodulename'")
          else if (index(keyword,'xid') == 1) then
            if (index(keyword,'=') /= 0) then
              read(keyword(index(keyword,'=') + 1 : len(trim(keyword))),&
              "(I6)",iostat = ios) fx_id
              if (ios /= 0) call terminate("xid should be written as 'xid=n'")
            else
              call terminate("xid should be written as 'xid=n'")
            end if
            write(60,'(A,I6)') '  exchange functional set to ', fx_id
          else if (index(keyword,'cid') == 1) then
            if (index(keyword,'=') /= 0) then
              read(keyword(index(keyword,'=') + 1 : len(trim(keyword))),&
              "(I6)",iostat = ios) fc_id
              if (ios /= 0) call terminate("cid should be written as 'cid=n'")
            else
              call terminate("cid should be written as 'cid=n'")
            end if
            write(60,'(A,I6)') '  correlation functional set to ', fc_id
          else if (index(keyword,'xhf') == 1) then
            if (index(keyword,'=') /= 0) then
              read(keyword(index(keyword,'=') + 1 : len(trim(keyword))),&
              "(F20.12)",iostat = ios) x_HF
              if (ios /= 0) call terminate(&
              "HF exchange component should be written as 'xhf=n'")
            else
              call terminate(&
              "HF exchange component should be written as 'xhf=n'")
            end if
            write(60,'(A,E10.3)') '  HF components set to ', x_HF
          else
            call terminate('unknown keyword detected in module Functional')
          end if
        end do module_functional
      else if (trim(module_name) == '%') then
        call terminate("module name should be written as '%modulename'")
      else
        call terminate('unknown module name detected')
      end if
    end do entire
    write(60,"(A)") "Other settings are default, exit module Atoms"
    close(12)
  end subroutine Load_keywords_tre
  
!-----------------------------------------------------------------------
!> check the self-consistency of computational settings
  subroutine Inputcheck()
    implicit none
    if (address_basis == '') call terminate('No basis set specified.')
    if (guess_type == 'molden') then
      address_MOLDEN = trim(address_job)//'.molden.input'
      inquire(file=address_MOLDEN, exist=exists, name=address_MOLDEN)
      if (.not. exists) then
        address_MOLDEN = trim(address_job)//'.molden'
        inquire(file=address_MOLDEN, exist=exists, name=address_MOLDEN)
        if (.not. exists) call terminate(&
        'MOLDEN file '//address_MOLDEN//' do not exist')
      end if
    else if (guess_type == 'fch') then

    else if (guess_type == 'ao2mo') then

    else
      call terminate('Unrecognisable initial guess setting.')
    end if
    electron_count = electron_count - charge
    if (electron_count <= 0) call terminate('charge is set incorrectly.')
    if (spin_mult == 675) then  ! default lowest spin
      if (mod(electron_count,2) == 0) then
        spin_mult = 1
      else
        spin_mult = 2
      end if
    else
      if (mod(electron_count,2) == 0) then
        if (mod(spin_mult,2) == 0 .or. spin_mult <= 0 .or. &
        spin_mult > electron_count+1) then
          call terminate('spin multiplicity is set incorrectly.')
        end if
      else
        if (mod(spin_mult,2) == 1 .or. spin_mult <= 0 .or. &
        spin_mult > electron_count+1) then
          call terminate('spin multiplicity is set incorrectly.')
        end if
      end if
    end if
    Nalpha = (electron_count-(spin_mult-1))/2 + (spin_mult-1)
    Nbeta  = (electron_count-(spin_mult-1))/2
    if (cspin/='n' .and. cspin/='d' .and. cspin/='f') call terminate(&
    'unrecognized constrained spin mode')
    if (.not. pVp1e .and. pppVp) call terminate(&
    'pppVp must be used in conjunction with pVp1e')
    if (.not. pVp1e .and. pVp2e) call terminate(&
    'pVp2e must be used in conjunction with pVp1e')
    if (subsp <= 1) call terminate('subsp two small')
    if (nodiis - subsp < 2) call terminate('nodiis should be set larger')
    if (damp > 0.91 .or. damp < 0.0) call terminate(&
    'damp shall be in range [0.05,0.9]')
    if (cutdamp < 0.0) call terminate('cutdamp should large than 0')
    if (diisdamp > 1.0 .or. diisdamp < 0.0) call terminate(&
    'diisdamp shall be in range [0,1]')
    if (cutdiis < 0.0) call terminate('cutdiis should large than 0')
    if (cutS < 0.0) call terminate('cutS should large than 0')
    if (x_HF > 1.000001 .or. x_HF < -0.000001) &
    call terminate('xhf should be in range [0,1]')
    if (beta(1) < 0.0 .or. beta(1) >= 1.0) call terminate("&
    betax should be in range [0,1)")
    if (beta(2) < 0.0 .or. beta(2) >= 1.0) call terminate("&
    betay should be in range [0,1)")
    if (beta(3) < 0.0 .or. beta(3) >= 1.0) call terminate("&
    betaz should be in range [0,1)")
    beta2 = sum(beta(:)**2)
    if (beta2 < 0.0 .or. beta2 >= 1.0) call terminate("&
    beta2 should be in range [0,1)")
    gamma = (1.0_dp-beta2)**(-0.5_dp)
  end subroutine Inputcheck

!-----------------------------------------------------------------------
!> normalization factor of GTFs
  pure elemental real(dp) function AON(a,l,m,n)
    implicit none
    real(dp),parameter  :: d(0:4) = &
    [1.0_dp, 1.0_dp, 1.0_dp/3.0_dp, 1.0_dp/15.0_dp, 1.0_dp/105.0_dp]
    real(dp),parameter  :: cfc = 2.0_dp/pi
    real(dp),intent(in) :: a
    integer,intent(in)  :: l,m,n
    AON = (cfc*a)**(0.75_dp)*dsqrt((4.0_dp*a)**(l+m+n) * d(l)*d(m)*d(n))
    return
  end function AON

!-----------------------------------------------------------------------
!> Assign matrix c2s, s2f and c2f
  subroutine Assign_csf(i_j)
    implicit none
    real(dp),allocatable,optional  :: i_j(:,:)      ! overlap matrix
    real(dp),allocatable           :: sph(:,:)
    integer                        :: ii, jj, kk, ll
    real(dp)                       :: min_evl       ! smallest eigenvalue of i_j
    allocate(c2soper(cbdm,sbdm), exc2soper(2*cbdm,2*sbdm))
    allocate(s2coper(sbdm,cbdm), exs2coper(2*sbdm,2*cbdm))
    allocate(c2s(cbdm,sbdm))
    c2s = 0.0_dp
    ii = 1 ! row
    kk = 1 ! column
    do while(ii <= cbdm)
      if (cbdata(ii) % L == 1) then ! S
        c2s(ii, kk) = 1.0_dp
        kk = kk + 1
        ii = ii + 1
      else if(cbdata(ii) % L == 2) then ! P
        do jj = 0, 2
          c2s(ii+jj, kk+jj) = 1.0_dp
        end do
        kk = kk + 3
        ii = ii + 3
      else if(cbdata(ii) % L == 3) then ! D
        do jj = 0, 5
          do ll = 0, 4
            c2s(ii+jj, kk+ll) = c2sd(1+jj, 1+ll)
          end do
        end do
        kk = kk + 5
        ii = ii + 6
      else if(cbdata(ii) % L == 4) then ! F
        do jj = 0, 9
          do ll = 0, 6
            c2s(ii+jj, kk+ll) = c2sf(1+jj, 1+ll)
          end do
        end do
        kk = kk + 7
        ii = ii + 10
      else if(cbdata(ii) % L == 5) then ! G
        do jj = 0, 14
          do ll = 0, 8
            c2s(ii+jj, kk+ll) = c2sg(1+jj, 1+ll)
          end do
        end do
        kk = kk + 9
        ii = ii + 15
      end if
    end do
    allocate(exc2s(2*cbdm,2*sbdm),source=c0)
    exc2s(1:cbdm,1:sbdm) = c2s * c1
    exc2s(cbdm+1:2*cbdm,sbdm+1:2*sbdm) = c2s * c1

    if (present(i_j)) then
      ! assign s2f and c2f
      allocate(sph(sbdm,sbdm))
      call csgo_1c(i_j)
      sph = i_j
      allocate(s2f(sbdm,sbdm))  ! be careful s2f and exs2f are not unitary
      call symm_orth(i_j, sbdm, s2f, min_evl)
      if (min_evl < 0.0_dp) &
      call terminate('evl(i_j) less than zero, may due to code error')
      if (min_evl < cutS) then
        deallocate(s2f)
        call can_orth(i_j, sbdm, s2f, fbdm)
      else
        fbdm = sbdm
      end if
      allocate(f2soper(fbdm,sbdm), exf2soper(2*fbdm,2*sbdm))
      allocate(s2foper(sbdm,fbdm), exs2foper(2*sbdm,2*fbdm))
      allocate(c2foper(cbdm,fbdm), exc2foper(2*cbdm,2*fbdm))
      !call sfgo_1c(i_j)
      allocate(c2f(cbdm,fbdm))
      call matmul('N', 'N', c2s, s2f, c2f)
      ! assign exs2f, exc2f
      allocate(exs2f(2*sbdm,2*fbdm),source=c0)
      exs2f(1:sbdm,1:fbdm) = s2f * c1
      exs2f(sbdm+1:2*sbdm,fbdm+1:2*fbdm) = s2f * c1
      allocate(exc2f(2*cbdm,2*fbdm),source=c0)
      exc2f(1:cbdm,1:fbdm) = c2f * c1
      exc2f(cbdm+1:2*cbdm,fbdm+1:2*fbdm) = c2f * c1
    end if

    if (present(i_j)) then
      ! assign f2s
      allocate(f2s(fbdm,sbdm), source=0.0_dp)
      allocate(exf2s(2*fbdm,2*sbdm), source=c0)
      call matmul('T', 'N', s2f, sph, f2s)
      exf2s(1:fbdm,1:sbdm) = f2s * c1
      exf2s(fbdm+1:2*fbdm,sbdm+1:2*sbdm) = f2s * c1
      deallocate(sph)
    end if
  end subroutine Assign_csf

!-----------------------------------------------------------------------
!> Assign matrix M_c2s and M_exc2s
  subroutine M_Assign_cs()
    implicit none
    integer            :: ii, jj, kk, ll
    if (allocated(M_c2s)) deallocate(M_c2s)
    allocate(M_c2s(M_cbdm,M_sbdm))
    M_c2s = 0.0_dp
    ii = 1 ! row
    kk = 1 ! column
    do while(ii <= M_cbdm)
      if (M_cbdata(ii) % L == 1) then ! S
        M_c2s(ii, kk) = 1.0_dp
        kk = kk + 1
        ii = ii + 1
      else if(M_cbdata(ii) % L == 2) then ! P
        do jj = 0, 2
          M_c2s(ii+jj, kk+jj) = 1.0_dp
        end do
        kk = kk + 3
        ii = ii + 3
      else if(M_cbdata(ii) % L == 3) then ! D
        do jj = 0, 5
          do ll = 0, 4
            M_c2s(ii+jj, kk+ll) = c2sd(1+jj, 1+ll)
          end do
        end do
        kk = kk + 5
        ii = ii + 6
      else if(M_cbdata(ii) % L == 4) then ! F
        do jj = 0, 9
          do ll = 0, 6
            M_c2s(ii+jj, kk+ll) = c2sf(1+jj, 1+ll)
          end do
        end do
        kk = kk + 7
        ii = ii + 10
      else if(M_cbdata(ii) % L == 5) then ! G
        do jj = 0, 14
          do ll = 0, 8
            M_c2s(ii+jj, kk+ll) = c2sg(1+jj, 1+ll)
          end do
        end do
        kk = kk + 9
        ii = ii + 15
      end if
    end do
    if (allocated(M_exc2s)) deallocate(M_exc2s)
    allocate(M_exc2s(2*M_cbdm,2*M_sbdm),source=c0)
    allocate(sphoper(cbdm,M_sbdm))
    M_exc2s(1:M_cbdm,1:M_sbdm) = M_c2s * c1
    M_exc2s(M_cbdm+1:2*M_cbdm,M_sbdm+1:2*M_sbdm) = M_c2s * c1
  end subroutine M_Assign_cs

!-----------------------------------------------------------------------
!> transfer real scalar matrix on Cartesian basis to matrix on spher-harmo basis
!!
!! m must be operator matrix
!!
!! m(cbdm,cbdm) -> m(sbdm,sbdm)
  subroutine csgo_1c(m)
    implicit none
    real(dp),allocatable :: m(:,:) ! target matrix
    if (.not. allocated(c2s)) call terminate('csgo_1c: c2s not allocated')
    if (.not. allocated(m)) call terminate('csgo_1c: matrix not allocated')
    call matmul('N', 'N', m, c2s, c2soper)
    deallocate(m)
    allocate(m(sbdm,sbdm))
    call matmul('T', 'N', c2s, c2soper, m)
  end subroutine csgo_1c

!-----------------------------------------------------------------------
!> transfer complex 2c matrix on Cartesian basis to matrix on spher-harmo basis
!!
!! m must be operator matrix
!!
!! m(2*cbdm,2*cbdm) -> m(2*sbdm,2*sbdm)
  subroutine csgo_2c(m)
    implicit none
    complex(dp),allocatable :: m(:,:) ! target matrix
    if (.not. allocated(exc2s)) call terminate('csgo_2c: exc2s not allocated')
    if (.not. allocated(m)) call terminate('csgo_2c: matrix not allocated')
    call matmul('N', 'N', m, exc2s, exc2soper)
    deallocate(m)
    allocate(m(2*sbdm,2*sbdm))
    call matmul('C', 'N', exc2s, exc2soper, m)
  end subroutine csgo_2c

!-----------------------------------------------------------------------
!> transfer real scalar matrix on spher-harmo basis to matrix on Cartesian basis
!!
!! m must be density matrix
!!
!! m(sbdm,sbdm) -> m(cbdm,cbdm)
  subroutine scgo_1c(m)
    implicit none
    real(dp),allocatable :: m(:,:) ! target matrix
    if (.not. allocated(c2s)) call terminate('scgo_1c: c2s not allocated')
    if (.not. allocated(m)) call terminate('scgo_1c: matrix not allocated')
    call matmul('N', 'T', m, c2s, s2coper)
    deallocate(m)
    allocate(m(cbdm,cbdm))
    call matmul('N', 'N', c2s, s2coper, m)
  end subroutine scgo_1c

!-----------------------------------------------------------------------
!> transfer complex 2c matrix on spher-harmo basis to matrix on Cartesian basis
!!
!! m must be density matrix
!!
!! m(2*sbdm,2*sbdm) -> m(2*cbdm,2*cbdm)
  subroutine scgo_2c(m)
    implicit none
    complex(dp),allocatable :: m(:,:) ! target matrix
    if (.not. allocated(exc2s)) call terminate('scgo_2c: exc2s not allocated')
    if (.not. allocated(m)) call terminate('scgo_2c: matrix not allocated')
    call matmul('N', 'C', m, exc2s, exs2coper)
    deallocate(m)
    allocate(m(2*cbdm,2*cbdm))
    call matmul('N', 'N', exc2s, exs2coper, m)
  end subroutine scgo_2c

!-----------------------------------------------------------------------
!> transfer real scalar matrix on spher-harmo basis to matrix on final SCF basis
!!
!! m must be operator matrix when use can_orth
!!
!! m(sbdm,sbdm) -> m(fbdm,fbdm)
  subroutine sfgo_1c(m)
    implicit none
    real(dp),allocatable :: m(:,:) ! target matrix
    if (.not. allocated(s2f)) call terminate('sfgo_1c: s2f not allocated')
    if (.not. allocated(m)) call terminate('sfgo_1c: matrix not allocated')
    call matmul('N', 'N', m, s2f, s2foper)
    deallocate(m)
    allocate(m(fbdm,fbdm))
    call matmul('T', 'N', s2f, s2foper, m)
  end subroutine sfgo_1c

!-----------------------------------------------------------------------
!> transfer complex 2c matrix on spher-harmo basis to matrix on final SCF basis
!!
!! m must be operator matrix when use can_orth
!!
!! m(2*sbdm,2*sbdm) -> m(2*fbdm,2*fbdm)
  subroutine sfgo_2c(m)
    implicit none
    complex(dp),allocatable :: m(:,:) ! target matrix
    if (.not. allocated(exs2f)) call terminate('sfgo_2c: exs2f not allocated')
    if (.not. allocated(m)) call terminate('sfgo_2c: matrix not allocated')
    call matmul('N', 'N', m, exs2f, exs2foper)
    deallocate(m)
    allocate(m(2*fbdm,2*fbdm))
    call matmul('C', 'N', exs2f, exs2foper, m)
  end subroutine sfgo_2c

!-----------------------------------------------------------------------
!> transfer real scalar matrix on final SCF basis to matrix on spher-harmo basis
!!
!! m must be density matrix when use can_orth
!!
!! m(fbdm,fbdm) -> m(sbdm,sbdm)
  subroutine fsgo_1c(m)
    implicit none
    real(dp),allocatable :: m(:,:) ! target matrix
    if (.not. allocated(f2s)) call terminate('fsgo_1c: f2s not allocated')
    if (.not. allocated(m)) call terminate('fsgo_1c: matrix not allocated')
    call matmul('N', 'N', m, f2s, f2soper)
    deallocate(m)
    allocate(m(sbdm,sbdm))
    call matmul('T', 'N', f2s, f2soper, m)
  end subroutine fsgo_1c

!-----------------------------------------------------------------------
!> transfer complex 2c matrix on final SCF basis to matrix on spher-harmo basis
!!
!! m must be density matrix when use can_orth
!!
!! m(2*fbdm,2*fbdm) -> m(2*sbdm,2*sbdm)
  subroutine fsgo_2c(m)
    implicit none
    complex(dp),allocatable :: m(:,:) ! target matrix
    if (.not. allocated(exf2s)) call terminate('fsgo_2c: exf2s not allocated')
    if (.not. allocated(m)) call terminate('fsgo_2c: matrix not allocated')
    call matmul('N', 'N', m, exf2s, exf2soper)
    deallocate(m)
    allocate(m(2*sbdm,2*sbdm))
    call matmul('C', 'N', exf2s, exf2soper, m)
  end subroutine fsgo_2c

!-----------------------------------------------------------------------
!> transfer real scalar matrix on Cartesian basis to matrix on final SCF basis
!!
!! m must be operator matrix
!!
!! m(cbdm,cbdm) -> m(fbdm,fbdm)
  subroutine cfgo_1c(m)
    implicit none
    real(dp),allocatable :: m(:,:) ! target matrix
    if (.not. allocated(c2f)) call terminate('cfgo_1c: c2f not allocated')
    if (.not. allocated(m)) call terminate('cfgo_1c: matrix not allocated')
    call matmul('N', 'N', m, c2f, c2foper)
    deallocate(m)
    allocate(m(fbdm,fbdm))
    call matmul('T', 'N', c2f, c2foper, m)
  end subroutine cfgo_1c

!-----------------------------------------------------------------------
!> transfer complex 2c matrix on Cartesian basis to matrix on final SCF basis
!!
!! m must be operator matrix
!!
!! m(2*cbdm,2*cbdm) -> m(2*fbdm,2*fbdm)
  subroutine cfgo_2c(m)
    implicit none
    complex(dp),allocatable :: m(:,:) ! target matrix
    if (.not. allocated(exc2f)) call terminate('cfgo_2c: exc2f not allocated')
    if (.not. allocated(m)) call terminate('cfgo_2c: matrix not allocated')
    call matmul('N', 'N', m, exc2f, exc2foper)
    deallocate(m)
    allocate(m(2*fbdm,2*fbdm))
    call matmul('C', 'N', exc2f, exc2foper, m)
  end subroutine cfgo_2c

!-----------------------------------------------------------------------
!> transfer matrix on Cartesian basis to matrix on spher-harmo basis
!!
!! m(cbdm,M_cbdm) -> m(sbdm,M_sbdm)
  subroutine M_csgo(m)
    implicit none
    real(dp),allocatable :: m(:,:) ! target matrix
    if (.not. allocated(M_c2s)) call terminate('M_csgo: M_c2s not allocated')
    if (.not. allocated(c2s)) call terminate('M_csgo: c2s not allocated')
    if (.not. allocated(m)) call terminate('M_csgo: matrix not allocated')
    call matmul('N', 'N', m, M_c2s, sphoper)
    deallocate(m)
    allocate(m(sbdm,M_sbdm))
    call matmul('T', 'N', c2s, sphoper, m)
  end subroutine M_csgo
  
end module Atoms