!> @file Representation.f90
!!
!! @brief representation transformation and frame transformation
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
    integer(8)                  :: nl(atom_count), nr(atom_count), cnr, cnl
    integer                     :: xloop_k, xloop_l
    integer                     :: Archannel, Aichannel, Brchannel, Bichannel
    integer                     :: xchannel, ychannel ,zchannel
    integer                     :: r          ! record position
    character(len=200)          :: env_TRESC
    character(len=100)          :: line
    logical                     :: mog_element
    !DIR$ ATTRIBUTES ALIGN:align_size :: datsa, datsb, pos3
    complex(dp)                 :: datsa(43400), datsb(43400)
    real(dp)                    :: pos3(434,3,100), trafopos3(434,3,100)
    character(len=100)          :: input_line
    inquire(file='.mogx',exist=exists)
    if (.not. exists) then
      open(newunit=xchannel, file=trim(jobname)//'.mogx', access='direct', &
      form='unformatted', recl=8, status='replace', action='write', iostat=ios)
      if (ios /= 0) call terminate('dump binary .mog failed')
      open(newunit=ychannel, file=trim(jobname)//'.mogy', access='direct', &
      form='unformatted', recl=8, status='replace', action='write', iostat=ios)
      if (ios /= 0) call terminate('dump binary .mog failed')
      open(newunit=zchannel, file=trim(jobname)//'.mogz', access='direct', &
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
    call getenv('TRESC',env_TRESC)
    open(30, file=trim(env_TRESC)//'/settings.ini', status='old', &
    action='read', iostat=ios)
    if (ios /= 0) then
      call terminate("can't open file "//trim(env_TRESC)//'/settings.ini')
    end if
    do
      read(30, '(A)', iostat=ios) line
      if (ios /= 0) then
        call terminate("no mogtype in "//trim(env_TRESC)//'/settings.ini')
      else if(index(line,'mogtype=') == 1) then
        call lowercase(line)
        if (line(index(line,'=')+1:index(line,'=')+1) == 'e') then
          mog_element = .true.
        else if(line(index(line,'=')+1:index(line,'=')+1) == 'i') then
          mog_element = .false.
        end if
        exit
      end if
    end do
    if (ios /= 0) then
      call terminate("can't read file "//trim(env_TRESC)//'/settings.ini')
    end if
    nr = 0
    nl = 0
    if (mog_element) then
      do
        read(30, '(A)', iostat=ios) line
        if (index(line,'//') == 1) cycle
        if (index(line,'nr[') /= 1) exit
        read(line(4:index(line,']=')-1), '(I)') r
        read(line(index(line,']=')+2:index(line,',')-1), '(I)') cnr
        line = line(index(line,']=')+2:)
        read(line(index(line,']=')+2:index(line,';')-1), '(I)') cnl
        do ii = 1, atom_count
          if (mol(ii)%atom_number == r) then
            nr(ii) = cnr
            nl(ii) = cnl
          end if
        end do
      end do
    else
      do
        read(30, '(A)', iostat=ios) line
        if (index(line,'//') == 1) cycle
        if (index(line,'nr{') /= 1) exit
        read(line(4:index(line,'}=')-1), '(I)') r
        read(line(index(line,'}=')+2:index(line,',')-1), '(I)') cnr
        line = line(index(line,'}=')+2:)
        read(line(index(line,'}=')+2:index(line,';')-1), '(I)') cnl
        nr(r) = cnr
        nl(r) = cnl
      end do
    end if
    close(30)
    if (is2c) then
      write(*,*) "Is this molecule in motion?"
      write(*,*) "if yes, input x,y,z components of velocity in NATRUAL UNIT"
      write(*,*) "if not, press ENTER"
      read(*, '(A)', iostat=ios) input_line
      if (ios == 0) read(input_line, *, iostat=ios) beta(1), beta(2), beta(3)
      if (ios /= 0 .or. len_trim(input_line) == 0) beta = 0.0
      write(*,*) "recieved..."
      beta2 = sum(beta(:)**2)
      if (beta2 < 0.0 .or. beta2 >= 1.0) call terminate("&
      beta2 should be in range [0,1)")
      gamma = (1.0_dp-beta2)**(-0.5_dp)
      call mol_frametrafo()
      if (beta2 > 1E-6) then
      write(*,*)
      write(*,*) "------------------------------------------------------------"
      write(*,*) "Electronic structure in motion frame differs from rest frame"
      write(*,*) "due to changes in measurement;"
      write(*,*) "Electronic structure in motion frame is not self-consistent."
      write(*,*) "------------------------------------------------------------"
      write(*,*)
      end if
    end if
    r = 1
    do ii = 1, atom_count
      call Chebyshev2_Lebedev_noweight(ii, nr(ii), nl(ii), pos3)
      !$omp parallel num_threads(16) default(shared) private(i)&
      !$omp& if(atom_count > 5)
      !$omp do schedule(static)
      do i = 1, nr(ii)
        if (is2c) then
          call grid_frametrafo_mo(&
          arr,nl(ii),pos3(1:nl(ii),:,i),trafopos3(1:nl(ii),:,i),&
          datsa((i-1)*nl(ii)+1:i*nl(ii)),datsb((i-1)*nl(ii)+1:i*nl(ii)))
        else
          ! always alpha
          call grid(arr,nl(ii),pos3(1:nl(ii),:,i),1,&
          datsa((i-1)*nl(ii)+1:i*nl(ii)))
        end if
      end do
      !$omp end do
      !$omp end parallel
      pos3 = trafopos3
      do jj = 1, nr(ii)
        do kk = 1, nl(ii)
          if (is2c) then
            write(Archannel, rec=r) real(datsa((jj-1)*nl(ii)+kk))
            write(Aichannel, rec=r) aimag(datsa((jj-1)*nl(ii)+kk))
            write(Brchannel, rec=r) real(datsb((jj-1)*nl(ii)+kk))
            write(Bichannel, rec=r) aimag(datsb((jj-1)*nl(ii)+kk))
            if (.not. exists) then
              write(xchannel, rec=r) trafopos3(kk,1,jj)
              write(ychannel, rec=r) trafopos3(kk,2,jj)
              write(zchannel, rec=r) trafopos3(kk,3,jj)
            end if
          else
            write(Archannel, rec=r) real(datsa((jj-1)*nl(ii)+kk))
            if (.not. exists) then
              write(xchannel, rec=r) pos3(kk,1,jj)
              write(ychannel, rec=r) pos3(kk,2,jj)
              write(zchannel, rec=r) pos3(kk,3,jj)
            end if
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
!> generate Gaussian cube file for input MO
  subroutine dump_cube(is2c, arr, title, spin)
    implicit none
    logical,intent(in)          :: is2c
    complex(dp),intent(in)      :: arr(2*cbdm)
    character(len=*),intent(in) :: title
    character(len=200)          :: env_TRESC
    character(len=100)          :: line
    integer,optional            :: spin
    integer                     :: end
    integer                     :: i, ii, jj, kk
    real(dp)                    :: maxx, minx, maxy, miny, maxz, minz
    real(dp)                    :: deltax, deltay, deltaz
    integer(8)                  :: nx, ny, nz
    !DIR$ ATTRIBUTES ALIGN:align_size :: points, datsa, datsb
    real(dp),allocatable        :: points(:,:)
    complex(dp),allocatable     :: datsa(:), datsb(:)
    !DIR$ ATTRIBUTES ALIGN:align_size :: dats
    complex(dp),allocatable     :: dats(:,:,:)

    call getenv('TRESC',env_TRESC)
    open(30, file=trim(env_TRESC)//'/settings.ini', status='old', &
    action='read', iostat=ios)
    if (ios /= 0) then
      call terminate("can't open file "//trim(env_TRESC)//'/settings.ini')
    end if
    do
      read(30, '(A)', iostat=ios) line
      if (ios /= 0) then
        call terminate("no nx,ny,nz in "//trim(env_TRESC)//'/settings.ini')
      else if(index(line,'nx=') == 1) then
        read(line(index(line,'=')+1:index(line,';')-1), '(I)') nx
      else if(index(line,'ny=') == 1) then
        read(line(index(line,'=')+1:index(line,';')-1), '(I)') ny
      else if(index(line,'nz=') == 1) then
        read(line(index(line,'=')+1:index(line,';')-1), '(I)') nz
        exit
      end if
    end do
    close(30)
    allocate(dats(2*nz,ny,nx))
    allocate(points(nz,3), datsa(nz), datsb(nz))
    if (is2c) then
      open(99, file=trim(title)//"-real.cub", status="replace", action="write")
      write(99,'(A)') 'generate by TRESC: real part of 2-component orb '&
      //trim(title)
      write(99,'(A,I,A)') 'total ',nx*ny*nz,' grid points'
      open(100, file=trim(title)//"-img.cub", status="replace", action="write")
      write(100,'(A)') 'generated by TRESC: imaginary part of 2-component orb '&
      //trim(title)
      write(100,'(A,I,A)') 'total ',nx*ny*nz,' grid points'
    else
      open(99, file=trim(title)//".cub", status="replace", action="write")
      write(99,'(A)') 'generated by TRESC'
      write(99,'(A,I,A)') 'total ',nx*ny*nz,' points'
    end if
    ! generate grid points
    maxx = maxval(mol(:)%pos(1)) + 3.0
    minx = minval(mol(:)%pos(1)) - 3.0
    deltax = (maxx-minx) / nx
    maxy = maxval(mol(:)%pos(2)) + 3.0
    miny = minval(mol(:)%pos(2)) - 3.0
    deltay = (maxy-miny) / ny
    maxz = maxval(mol(:)%pos(3)) + 3.0
    minz = minval(mol(:)%pos(3)) - 3.0
    deltaz = (maxz-minz) / nz
    if (is2c) then
      write(99,'(I4,3F14.8)') -atom_count, minx, miny, minz
      write(99,'(I3,3F14.8)') nx, deltax, 0.0, 0.0
      write(99,'(I3,3F14.8)') ny, 0.0, deltay, 0.0
      write(99,'(I3,3F14.8)') nz, 0.0, 0.0, deltaz
      write(100,'(I4,3F14.8)') -atom_count, minx, miny, minz
      write(100,'(I3,3F14.8)') nx, deltax, 0.0, 0.0
      write(100,'(I3,3F14.8)') ny, 0.0, deltay, 0.0
      write(100,'(I3,3F14.8)') nz, 0.0, 0.0, deltaz
      do ii = 1, atom_count
        write(99,'(I4,4F14.8)') mol(ii)%atom_number, real(mol(ii)%atom_number),&
        mol(ii)%pos(1), mol(ii)%pos(2), mol(ii)%pos(3)
        write(100,'(I4,4F14.8)')mol(ii)%atom_number, real(mol(ii)%atom_number),&
        mol(ii)%pos(1), mol(ii)%pos(2), mol(ii)%pos(3)
      end do
      write(99,'(I3,A4,A4)') 2, trim(title), trim(title)
      write(100,'(I3,A4,A4)') 2, trim(title), trim(title)
    else
      write(99,'(I4,3F14.8)') -atom_count, minx, miny, minz
      write(99,'(I3,3F14.8)') nx, deltax, 0.0, 0.0
      write(99,'(I3,3F14.8)') ny, 0.0, deltay, 0.0
      write(99,'(I3,3F14.8)') nz, 0.0, 0.0, deltaz
      do ii = 1, atom_count
        write(99,'(I4,4F14.8)') mol(ii)%atom_number, real(mol(ii)%atom_number),&
        mol(ii)%pos(1), mol(ii)%pos(2), mol(ii)%pos(3)
      end do
      do ii = 1, len(title)
        if (is_alpha(title(ii:ii))) exit
      end do
      write(99,'(I3,A4)') 1, trim(title(:ii-1)//title(ii+1:))
    end if
    !$omp parallel num_threads(10) default(shared) private(i, ii, jj, kk, &
    !$omp& points, datsa, datsb) if(atom_count > 5)
    !$omp do schedule(static)
    do i = 1, nx
      ii = i
      points(:,1) = minx + (ii-1) * deltax
      do jj = 1, ny
        points(:,2) = miny + (jj-1) * deltay
        do kk = 1, nz
          points(kk,3) = minz + (kk-1) * deltaz
        end do
        if (is2c) then
          call grid(arr, nz, points, 1, datsa)
          call grid(arr, nz, points, 2, datsb)
          do kk = 1, nz
            dats(2*kk-1,jj,ii) = datsa(kk)
            dats(2*kk,jj,ii) = datsb(kk)
          end do
        else
          call grid(arr, nz, points, spin, dats(1:nz,jj,ii))
        end if
      end do
    end do
    !$omp end do
    !$omp end parallel
    do ii = 1, nx
      do jj = 1, ny
        if (is2c) then
          do kk = 1, 2*nz, 6
            end = min(kk+5, 2*nz)
            write(99, '(6ES15.6)') real(dats(kk:end,jj,ii))
            write(100, '(6ES15.6)') aimag(dats(kk:end,jj,ii))
          end do
        else
          do kk = 1, nz, 6
            end = min(kk+5, nz)
            write(99, '(6ES15.6)') real(dats(kk:end,jj,ii))
          end do
        end if
      end do
    end do
    close(99)
    if (is2c) close(100)
    deallocate(dats, points, datsa, datsb)
  end subroutine dump_cube

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
    integer               :: ii, jj
    integer               :: nl2
    !DIR$ ATTRIBUTES ALIGN:align_size :: wr, R
    real(dp)              :: wr(nr), R(nr)
    real(dp),intent(out)  :: pos3w1(434,4,100) ! x, y, z, weight
    !DIR$ ATTRIBUTES ALIGN:align_size :: lx, ly, lz, wl
    real(dp)              :: lx(nl), ly(nl), lz(nl), wl(nl)
    real(dp)              :: p, x_i, point(3), spacew
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
    integer               :: ii, jj
    integer               :: nl2
    !DIR$ ATTRIBUTES ALIGN:align_size :: R
    real(dp)              :: R(nr)
    real(dp),intent(out)  :: pos3(434,3,100) ! x, y, z
    !DIR$ ATTRIBUTES ALIGN:align_size :: lx, ly, lz, wl
    real(dp)              :: lx(nl), ly(nl), lz(nl), wl(nl)
    real(dp)              :: p, x_i
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
    real(dp), intent(in)   :: point(3)
    integer(8), intent(in) :: count
    real(dp), intent(out)  :: weight
    real(dp)               :: weights(atom_count)
    real(dp)               :: ri, rj, R, chi, miu_, miu, a, niu, s
    integer                :: ii, jj
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
    integer,intent(in)  :: p ! number of remaining iterations
    if (p <= 1) then
      y = 1.5_dp*x - 0.5_dp*x**3
    else
      y = 1.5_dp*Becke_miu(x, p-1) - 0.5_dp*Becke_miu(x, p-1)**3
    end if
  end function Becke_miu

!------------------------------------------------------------
!> assign molecule in frame of motion (trafomol)
  subroutine mol_frametrafo()
    implicit none
    integer :: ii, jj, kk
    real(dp) :: trafocoe, dot
    ! Lorentz transformation of x, leads to length contraction
    ! x' = x - gamma/(gamma+1.0_dp)(x.beta) * beta
    trafomol = mol
    trafocoe = -gamma/(gamma+1.0_dp)
    do ii = 1, atom_count
      dot = trafocoe*sum(beta(:)*mol(ii)%pos(:))
      trafomol(ii)%pos(:) = mol(ii)%pos(:) + dot*beta(:)
    end do
  end subroutine mol_frametrafo

!------------------------------------------------------------
!> reference frame transformation of a (basis repre) vector
!!
!! input "points" is the coordinates in rest frame
!!
!! “simultaneous” in motion frame: x_rest = L(x_motion)|t_motion=0,
!! in this case, x' = x_rest; x = x_motion
  pure subroutine grid_frametrafo_mo(arr, n, points, trafopoints, datsa, datsb)
    implicit none
    complex(dp),intent(in)  :: arr(:)
    integer(8),intent(in)   :: n
    real(dp),intent(in)     :: points(n, 3)
    real(dp),intent(out)    :: trafopoints(n, 3)
    complex(dp),intent(out) :: datsa(n), datsb(n)
    real(dp)                :: sz(n), trafosz(n), amp(n)
    real(dp)                :: a2, b2
    integer                 :: ii, jj, kk
    integer                 :: contr
    real(dp)                :: vec(n, 3)
    integer                 :: L, M
    integer                 :: fac(3)
    real(dp)                :: coeff
    real(dp)                :: expo
    real(dp)                :: val(n)
    real(dp)                :: codcoe, coddot, spincoe
    ! Lorentz transformation of x, leads to length contraction
    ! x = x' + gamma**2/(gamma+1.0_dp)(x'.beta) * beta
    codcoe = -gamma/(gamma+1.0_dp)
    do ii = 1, n
      coddot = codcoe*sum(beta(:)*points(ii,:))
      trafopoints(ii,:) = points(ii,:) + coddot*beta(:)
    end do
    
    datsa = c0
    datsb = c0
    do kk = 1, cbdm
      val = 0.0_dp
      vec(:,1) = points(:,1) - basis_inf(kk) % pos(1)
      vec(:,2) = points(:,2) - basis_inf(kk) % pos(2)
      vec(:,3) = points(:,3) - basis_inf(kk) % pos(3)
      contr    = basis_inf(kk) % contr
      L        = basis_inf(kk) % L
      M        = basis_inf(kk) % M
      fac(:)   = AO_fac(:,L,M)
      do jj = 1, contr
        expo   = basis_inf(kk) % expo(jj)
        coeff  = basis_inf(kk) % Ncoe(jj,M)
        val(:) = val(:) + coeff * exp(-expo*(sum(vec(:,:)**2,dim=2)))
      end do
      val(:) = val(:) * (vec(:,1)**fac(1))
      val(:) = val(:) * (vec(:,2)**fac(3))
      val(:) = val(:) * (vec(:,3)**fac(3))
      datsa  = datsa + val * arr(kk)
      datsb  = datsb + val * arr(kk+cbdm)
    end do

    ! Imitative transformation of s, causes s_z to away from +-1/2,
    ! also means mixing of alpha and beta components
    ! s' = (1-gamma/(gamma+1)*beta(3)**2)/(1-beta(3)**2)**0.5 * s
    spincoe = (1.0_dp+codcoe*beta(3)**2) * (1.0_dp-beta(3)**2)**(-0.5_dp)
    sz = real(datsa*conjg(datsa) - datsb*conjg(datsb))
    amp = real(datsa*conjg(datsa) + datsb*conjg(datsb))
    trafosz = spincoe * sz
    do ii = 1, n
      a2 = 0.5_dp*(trafosz(ii)+amp(ii))
      b2 = 0.5_dp*(amp(ii)-trafosz(ii))
      datsa(ii) = datsa(ii) * a2/(datsa(ii)*conjg(datsa(ii)))
      datsb(ii) = datsb(ii) * b2/(datsb(ii)*conjg(datsb(ii)))
    end do
  end subroutine grid_frametrafo_mo

end module Representation