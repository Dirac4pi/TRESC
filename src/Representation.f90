!> @file Representation.f90
!!
!! @brief representation transformation of states
!!
!! @syntax Fortran 2008 free format
!!
!! @code UTF-8
!!
!! @author dirac4pi
module Representation
  use Atoms
  use Functional
  use Fundamentals
  use Lebedev

  contains

!------------------------------------------------------------
!> transform states from basis repre to real space repre (to generate .mog file)
!!
!! grad=0:grid =1:gradient =2:Laplacian =a:alpha =b:beta
  subroutine mogrid_becke(is2c, arr, title)
    implicit none
    logical,intent(in)          :: is2c       ! true: 2c, false: 1c
    complex(dp),intent(in)      :: arr(:)     ! input array (basis repre)
    character(len=*),intent(in) :: title      ! title of .mog file
    integer(8)                  :: i, ii, jj, kk, ll  ! loop varables
    integer(8)                  :: nl, nr
    integer                     :: xloop_k, xloop_l
    integer                     :: Archannel, Aichannel, Brchannel, Bichannel
    integer                     :: xchannel, ychannel ,zchannel
    integer                     :: r          ! record position
    integer(8)                  :: npoints    ! number of grid points
    !DIR$ ATTRIBUTES ALIGN:align_size :: datsa, datsb, pos3
    complex(dp)                 :: datsa(43400), datsb(43400)
    real(dp)                    :: pos3(434,3,100)
    inquire(file='.mogx',exist=exists)
    if (.not. exists) then
      open(newunit=xchannel, file='.mogx', access='direct', &
      form='unformatted', recl=8, status='replace', action='write', iostat=ios)
      if (ios /= 0) call terminate('dump binary .mog failed')
      open(newunit=ychannel, file='.mogy', access='direct', &
      form='unformatted', recl=8, status='replace', action='write', iostat=ios)
      if (ios /= 0) call terminate('dump binary .mog failed')
      open(newunit=zchannel, file='.mogz', access='direct', &
      form='unformatted', recl=8, status='replace', action='write', iostat=ios)
      if (ios /= 0) call terminate('dump binary .mog failed')
    end if
    if (is2c) then
      open(newunit=Archannel, file=trim(title)//'.mogar', access='direct', &
      form='unformatted', recl=8, status='replace', action='write', iostat=ios)
      if (ios /= 0) call terminate('dump binary .mog failed')
      open(newunit=Aichannel, file=trim(title)//'.mogai', access='direct', &
      form='unformatted', recl=8, status='replace', action='write', iostat=ios)
      if (ios /= 0) call terminate('dump binary .mog failed')
      open(newunit=Brchannel, file=trim(title)//'.mogbr', access='direct', &
      form='unformatted', recl=8, status='replace', action='write', iostat=ios)
      if (ios /= 0) call terminate('dump binary .mog failed')
      open(newunit=Bichannel, file=trim(title)//'.mogbi', access='direct', &
      form='unformatted', recl=8, status='replace', action='write', iostat=ios)
      if (ios /= 0) call terminate('dump binary .mog failed')
    else
      open(newunit=Archannel, file=trim(title)//'.mogv', access='direct', &
      form='unformatted', recl=8, status='replace', action='write', iostat=ios)
      if (ios /= 0) call terminate('dump binary .mog failed')
    end if
    npoints = 0
    do ii = 1, atom_count
      if (ii == 1 .or. ii == 3 .or. ii == 4 .or. (ii >= 67 .and. ii <= 75)) then
        if (mol(ii)%atom_number <= 2) then
          nr = 35
          nl = 230
        else if (mol(ii)%atom_number <= 10) then
          nr = 65
          nl = 434
        else if (mol(ii)%atom_number <= 18) then
          nr = 80
          nl = 434
        else
          nr = 100
          nl = 434
        end if
      else
        nr = 0
        nl = 0
      end if
      npoints = npoints + nr*nl
    end do
    if (is2c) then
      write(Archannel,rec=1) npoints
      write(Aichannel,rec=1) npoints
      write(Brchannel,rec=1) npoints
      write(Bichannel,rec=1) npoints
    else
      write(Archannel,rec=1) npoints
    end if
    if (.not. exists) then
      write(xchannel,rec=1) npoints
      write(ychannel,rec=1) npoints
      write(zchannel,rec=1) npoints
    end if
    r = 2
    do ii = 1, atom_count
      if (ii == 1 .or. ii == 3 .or. ii == 4 .or. (ii >= 67 .and. ii <= 75)) then
        if (mol(ii)%atom_number <= 2) then
          nr = 35
          nl = 230
        else if (mol(ii)%atom_number <= 10) then
          nr = 65
          nl = 434
        else if (mol(ii)%atom_number <= 18) then
          nr = 80
          nl = 434
        else
          nr = 100
          nl = 434
        end if
      else
        nr = 0
        nl = 0
      end if
      call Chebyshev2_Lebedev_noweight(ii, nr, nl, pos3)
      !$omp parallel num_threads(14) default(shared) private(i)&
      !$omp& if(atom_count > 20)
      !$omp do schedule(static)
      do i = 1, nr
        if (is2c) then
          call grid(arr, nl, pos3(1:nl,:,i), 1, datsa((i-1)*nl+1:i*nl))
          call grid(arr, nl, pos3(1:nl,:,i), 2, datsb((i-1)*nl+1:i*nl))
        else
          ! always alpha
          call grid(arr, nl, pos3(1:nl,:,i), 1, datsa((i-1)*nl+1:i*nl))
        end if
      end do
      !$omp end do
      !$omp end parallel
      do jj = 1, nr
        do kk = 1, nl
          if (is2c) then
            write(Archannel, rec=r) real(datsa((jj-1)*nl+kk))
            write(Aichannel, rec=r) aimag(datsa((jj-1)*nl+kk))
            write(Brchannel, rec=r) real(datsb((jj-1)*nl+kk))
            write(Bichannel, rec=r) aimag(datsb((jj-1)*nl+kk))
          else
            write(Archannel, rec=r) real(datsa((jj-1)*nl+kk))
          end if
          if (.not. exists) then
            write(xchannel, rec=r) pos3(kk,1,jj)
            write(ychannel, rec=r) pos3(kk,2,jj)
            write(zchannel, rec=r) pos3(kk,3,jj)
          end if
          r = r + 1
        end do
      end do
    end do
    if (is2c) then
      close(Archannel)
      close(Aichannel)
      close(Brchannel)
      close(Bichannel)
    else
      close(Archannel)
    end if
    if (.not. exists) then
      close(xchannel)
      close(ychannel)
      close(zchannel)
    end if
  end subroutine mogrid_becke

!------------------------------------------------------------
!> transform states from basis repre to real space repre (Becke's fuzzy grid)
!!
!! directly used to generate xc energies for xc functional
!!
!! returns only the integral value
  subroutine basis2grid_Becke(rho_m, ex, ec, Fockx, Fockc)
    implicit none
    real(dp),intent(out)   :: ex         ! exchange functional
    real(dp), optional     :: ec         ! correlation functional
    complex(dp),intent(in) :: rho_m(2*cbdm, 2*cbdm)! Cartesian density matrix
    integer(8)             :: i          ! loop variables for basis2grid_Becke
    integer(8)             :: nr, nl     ! number of grid points
    !DIR$ ATTRIBUTES ALIGN:align_size :: pos3w1
    real(dp)               :: pos3w1(434,4,100)
    complex(dp),intent(out):: Fockx(2*cbdm, 2*cbdm)
    complex(dp),optional   :: Fockc(2*cbdm, 2*cbdm)
    !DIR$ ATTRIBUTES ALIGN:align_size :: Fockxmic, Fockcmic
    real(dp)               :: Fockxmic(2*cbdm, 2*cbdm)
    real(dp)               :: Fockcmic(2*cbdm, 2*cbdm)
    real(dp)               :: exm, ecm
    Fockx = c0
    ex = 0.0_dp
    if (present(Fockc)) Fockc = c0
    if (present(ec)) ec = 0.0_dp
    ! single-centre Gaussian weight integral
    ! Processes all the gird lattice points of an atom center at once,
    ! improving efficiency and avoiding stack overflows

    !$omp parallel num_threads(min(threads,atom_count)) default(shared)&
    !$omp& private(i, nr, nl, Fockxmic, Fockcmic, exm, ecm, pos3w1)&
    !$omp& if(threads < nproc)
    !$omp do schedule(dynamic, 1)
    do i = 1, atom_count
      if (mol(i)%atom_number <= 2) then
        nr = 35
        nl = 434
      else if (mol(i)%atom_number <= 10) then
        nr = 65
        nl = 434
      else if (mol(i)%atom_number <= 18) then
        nr = 80
        nl = 434
      else
        nr = 100
        nl = 434
      end if
      call Chebyshev2_Lebedev(i, nr, nl, pos3w1)
      if (fc_id /= -1) then
        call Fockxc_main(rho_m, nr, nl, pos3w1, exm, ecm, Fockxmic, Fockcmic)
      else
        call Fockxc_main(rho_m=rho_m, nr=nr, nl=nl, pos3w1=pos3w1, &
        ex=exm, matx=Fockxmic)
      end if
    end do
    !$omp end do
    !$omp critical
    Fockx = Fockx + Fockxmic
    ex = ex + exm
    if (present(Fockc)) Fockc = Fockc + Fockcmic
    if (present(ec)) ec = ec + ecm
    !$omp end critical
    !$omp end parallel
  end subroutine basis2grid_Becke

!------------------------------------------------------------
!> assign radial grid points (Gauss weight) using 2nd Chebyshev method
!!
!! and combined with the Lebedev sphere grid to calculate single-centre integral
  pure subroutine Chebyshev2_Lebedev(count, nr, nl, pos3w1)
    implicit none
    integer(8),intent(in) :: count, nr, nl
    integer :: ii, jj
    integer :: nl2
    !DIR$ ATTRIBUTES ALIGN:align_size :: wr, R
    real(dp) :: wr(nr), R(nr)
    real(dp),intent(out) :: pos3w1(434,4,100) ! x, y, z, weight
    !DIR$ ATTRIBUTES ALIGN:align_size :: lx, ly, lz, wl
    real(dp) :: lx(nl), ly(nl), lz(nl), wl(nl)
    real(dp) :: p, x_i, point(3), spacew
    ! assign radial grid points (Gauss weight) using 2nd Chebyshev method
    if (mol(count)%atom_number == 1) then
      p = CSD_CovR(mol(count)%atom_number)
    else
      p = CSD_CovR(mol(count)%atom_number) * 0.5_dp
    end if
    do ii = 1, nr
      x_i = cos(ii*pi / real(nr+1))
      R(ii) = (1+x_i) / (1-x_i) * p
      wr(ii) = (2.0_dp*pi / real(nr+1)) * p**3 * ((x_i+1)**2.5 / (1-x_i)**3.5)
    end do
    if (nl == 230) then
      call LD0230(lx, ly, lz, wl, nl2)
    else if (nl == 434) then
      call LD0434(lx, ly, lz, wl, nl2)
    end if
    wl = wl * 4.0_dp * pi ! sphere weight
    do ii = 1, nr
      pos3w1(1:nl,1,ii) = R(ii) * lx(:)
      pos3w1(1:nl,2,ii) = R(ii) * ly(:)
      pos3w1(1:nl,3,ii) = R(ii) * lz(:)
      pos3w1(1:nl,4,ii) = wr(ii) * wl(:)
    end do
    pos3w1(1:nl,1,1:nr) = pos3w1(1:nl,1,1:nr) + mol(count)%pos(1)
    pos3w1(1:nl,2,1:nr) = pos3w1(1:nl,2,1:nr) + mol(count)%pos(2)
    pos3w1(1:nl,3,1:nr) = pos3w1(1:nl,3,1:nr) + mol(count)%pos(3)
    do ii = 1, nr
      do jj = 1, nl
        point(:) = pos3w1(jj,1:3,ii)
        call Becke_weight(point, count, spacew)
        pos3w1(jj,4,ii) = pos3w1(jj,4,ii) * spacew
      end do
    end do
  end subroutine Chebyshev2_Lebedev

!------------------------------------------------------------
!> assign radial grid points (no weight) using 2nd Chebyshev method
!!
!! and combined with the Lebedev sphere grid to calculate single-centre integral
  pure subroutine Chebyshev2_Lebedev_noweight(count, nr, nl, pos3)
    implicit none
    integer(8),intent(in) :: count, nr, nl
    integer :: ii, jj
    integer :: nl2
    !DIR$ ATTRIBUTES ALIGN:align_size :: R
    real(dp) :: R(nr)
    real(dp),intent(out) :: pos3(434,3,100) ! x, y, z
    !DIR$ ATTRIBUTES ALIGN:align_size :: lx, ly, lz, wl
    real(dp) :: lx(nl), ly(nl), lz(nl), wl(nl)
    real(dp) :: p, x_i
    ! assign radial grid points (Gauss weight) using 2nd Chebyshev method
    if (mol(count)%atom_number == 1) then
      p = CSD_CovR(mol(count)%atom_number)
    else
      p = CSD_CovR(mol(count)%atom_number) * 0.5_dp
    end if
    do ii = 1, nr
      x_i = cos(ii*pi / real(nr+1))
      R(ii) = (1+x_i) / (1-x_i) * p
    end do
    if (nl == 230) then
      call LD0230(lx, ly, lz, wl, nl2)
    else if (nl == 434) then
      call LD0434(lx, ly, lz, wl, nl2)
    end if
    do ii = 1, nr
      pos3(1:nl,1,ii) = R(ii) * lx(:)
      pos3(1:nl,2,ii) = R(ii) * ly(:)
      pos3(1:nl,3,ii) = R(ii) * lz(:)
    end do
    pos3(1:nl,1,1:nr) = pos3(1:nl,1,1:nr) + &
    mol(count)%pos(1)
    pos3(1:nl,2,1:nr) = pos3(1:nl,2,1:nr) + &
    mol(count)%pos(2)
    pos3(1:nl,3,1:nr) = pos3(1:nl,3,1:nr) + &
    mol(count)%pos(3)
  end subroutine Chebyshev2_Lebedev_noweight

!------------------------------------------------------------
!> calculate weight of a given coordinate in Becke's fuzzy grid
  pure subroutine Becke_weight(point, count, weight)
    implicit none
    real(dp), intent(in) :: point(3)
    integer(8), intent(in) :: count
    real(dp), intent(out) :: weight
    real(dp) :: weights(atom_count)
    real(dp) :: ri, rj, R, chi, miu_, miu, a, niu, s
    integer :: ii, jj
    weights = 1.0_dp
    do ii = 1, atom_count
      ri = dsqrt(sum((point(:)-mol(ii)%pos(:))**2))
      do jj = 1, atom_count
        if (jj == ii) cycle
        rj = dsqrt(sum((point(:)-mol(jj)%pos(:))**2))
        R = dsqrt(sum((mol(ii)%pos(:)-&
        mol(jj)%pos(:))**2))
        miu = (ri - rj) / R
        chi = CSD_CovR(mol(ii)%atom_number) / &
        CSD_CovR(mol(jj)%atom_number)
        miu_ = (chi - 1.0_dp) / (chi + 1.0_dp)
        a = miu_ / (miu_**2 - 1.0_dp)
        if (a > 0.5_dp) a = 0.5_dp
        if (a < -0.5_dp) a = -0.5_dp
        niu = miu + a*(1.0_dp-miu**2)
        s = 0.5_dp*(1.0_dp-Becke_miu(niu, 3))
        weights(ii) = weights(ii) * s
      end do
    end do
    weight = weights(count) / sum(weights)
  end subroutine Becke_weight

!------------------------------------------------------------
!> parameter iterations for weighting curve, iteration, minimum N of 1
  pure recursive real(dp) function Becke_miu(x, p) result(y)
    implicit none
    real(dp),intent(in) :: x
    integer,intent(in) :: p ! number of remaining iterations
    if (p <= 1) then
      y = 1.5_dp*x - 0.5_dp*x**3
    else
      y = 1.5_dp*Becke_miu(x, p-1) - 0.5_dp*Becke_miu(x, p-1)**3
    end if
  end function Becke_miu
end module Representation