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
!> transform states from basis repre to real space repre (Becke's fuzzy grid)
!!
!! generate .mog file can be used to visualize MOs via vis2c
  subroutine Basis2real_Becke_mog(is2c, arr, title)
    implicit none
    logical,intent(in)          :: is2c       ! true: 2c, false: 1c
    complex(dp),intent(in)      :: arr(:)     ! input array (basis repre)
    character(len=*),intent(in) :: title      ! title of .mog file
    integer(8)                  :: i, ii, jj, kk, ll  ! loop variables
    integer(8)                  :: nl(atom_count), nr(atom_count), cnr, cnl
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
    open(30, file=trim(env_TRESC)//'/gridsettings.ini', status='old', &
    action='read', iostat=ios)
    if (ios /= 0) then
      call terminate("can't open file "//trim(env_TRESC)//'/gridsettings.ini')
    end if
    do
      read(30, '(A)', iostat=ios) line
      if (ios /= 0) then
        call terminate("no mogtype in "//trim(env_TRESC)//'/gridsettings.ini')
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
    nr = 0
    nl = 0
    if (mog_element) then
      do
        read(30, '(A)', iostat=ios) line
        if (ios /= 0) exit
        if (index(line,'#') == 1) cycle
        if (index(line,'nr[') == 1) then
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
        end if
      end do
    else
      do
        read(30, '(A)', iostat=ios) line
        if (ios /= 0) exit
        if (index(line,'#') == 1) cycle
        if (index(line,'nr{') == 1) then
          read(line(4:index(line,'}=')-1), '(I)') r
          read(line(index(line,'}=')+2:index(line,',')-1), '(I)') cnr
          line = line(index(line,'}=')+2:)
          read(line(index(line,'}=')+2:index(line,';')-1), '(I)') cnl
          nr(r) = cnr
          nl(r) = cnl
        end if
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
      call Mol_real_Lorentztrafo()
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
          call Grid_real_Lorentztrafo(&
          arr,nl(ii),pos3(1:nl(ii),:,i),trafopos3(1:nl(ii),:,i),&
          datsa((i-1)*nl(ii)+1:i*nl(ii)),datsb((i-1)*nl(ii)+1:i*nl(ii)))
        else
          ! always alpha
          call Grid_real(arr,nl(ii),pos3(1:nl(ii),:,i),1,&
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
  end subroutine Basis2real_Becke_mog

!------------------------------------------------------------
!> transform states from basis repre to real space repre (uniform cube grid)
!!
!! generate .cub file can be used to visualize MOs via vis2c
!!
!! spin(optional)=0:total =1:alpha =2:beta
  subroutine Basis2real_cube(is2c, arr, title, spin)
    implicit none
    logical,intent(in)          :: is2c
    complex(dp),intent(in)      :: arr(2*cbdm)
    character(len=*),intent(in) :: title
    character(len=200)          :: env_TRESC
    character(len=100)          :: line
    integer,optional            :: spin
    integer                     :: end
    integer                     :: i, ii, jj, kk
    real(dp)                    :: maxx, minx, maxy, miny, maxz, minz, edge
    real(dp)                    :: deltax, deltay, deltaz
    integer(8)                  :: nx, ny, nz
    !DIR$ ATTRIBUTES ALIGN:align_size :: points, datsa, datsb
    real(dp),allocatable        :: points(:,:)
    complex(dp),allocatable     :: datsa(:), datsb(:)
    !DIR$ ATTRIBUTES ALIGN:align_size :: dats
    complex(dp),allocatable     :: dats(:,:,:)

    call getenv('TRESC',env_TRESC)
    open(30, file=trim(env_TRESC)//'/gridsettings.ini', status='old', &
    action='read', iostat=ios)
    if (ios /= 0) call terminate(&
    "can't open file "//trim(env_TRESC)//'/gridsettings.ini')
    nx = 64        ! default settings
    ny = 64
    nz = 64
    edge = 3.0_dp
    do
      read(30, '(A)', iostat=ios) line
      if (ios /= 0) then
        exit
      else if(index(line,'nx=') == 1) then
        read(line(index(line,'=')+1:index(line,';')-1), '(I)') nx
        if (nx <= 0) call terminate("nx should be positive integer")
      else if(index(line,'ny=') == 1) then
        read(line(index(line,'=')+1:index(line,';')-1), '(I)') ny
        if (ny <= 0) call terminate("ny should be positive integer")
      else if(index(line,'nz=') == 1) then
        read(line(index(line,'=')+1:index(line,';')-1), '(I)') nz
        if (nz <= 0) call terminate("nz should be positive integer")
      else if(index(line,'edge=') == 1) then
        read(line(index(line,'=')+1:index(line,';')-1), '(F)') edge
        if (edge <= 0.0) call terminate("edge should be positive float")
      end if
    end do
    close(30)
    allocate(dats(2*nz,ny,nx))
    allocate(points(nz,3), datsa(nz), datsb(nz))
    if (is2c) then
      open(99,file=trim(title)//"-realreal.cub",status="replace",action="write")
      write(99,'(A)') 'generate by TRESC: real part of 2-component orb '&
      //trim(title)//' in real space'
      write(99,'(A,I,A)') 'total ',nx*ny*nz,' grid points'
      open(100,file=trim(title)//"-realimg.cub",status="replace",action="write")
      write(100,'(A)') 'generated by TRESC: imaginary part of 2-component orb '&
      //trim(title)//' in real space'
      write(100,'(A,I,A)') 'total ',nx*ny*nz,' grid points'
    else
      open(99, file=trim(title)//"-real.cub", status="replace", action="write")
      write(99,'(A)') 'generated by TRESC: scalar orb '&
      //trim(title)//' in real space'
      write(99,'(A,I,A)') 'total ',nx*ny*nz,' points'
    end if
    ! generate grid points
    maxx = maxval(mol(:)%pos(1)) + edge
    minx = minval(mol(:)%pos(1)) - edge
    deltax = (maxx-minx) / nx
    maxy = maxval(mol(:)%pos(2)) + edge
    miny = minval(mol(:)%pos(2)) - edge
    deltay = (maxy-miny) / ny
    maxz = maxval(mol(:)%pos(3)) + edge
    minz = minval(mol(:)%pos(3)) - edge
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
    !$omp parallel num_threads(16) default(shared) private(i, ii, jj, kk, &
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
          call Grid_real(arr, nz, points, 1, datsa)
          call Grid_real(arr, nz, points, 2, datsb)
          do kk = 1, nz
            dats(2*kk-1,jj,ii) = datsa(kk)
            dats(2*kk,jj,ii) = datsb(kk)
          end do
        else
          call Grid_real(arr, nz, points, spin, dats(1:nz,jj,ii))
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
  end subroutine Basis2real_cube

!------------------------------------------------------------
!> transform states from basis repre to momentum space repre (uniform cube grid)
!!
!! generate .cub file can be used to visualize MOs via vis2c
!!
!! spin(optional)=0:total =1:alpha =2:beta
  subroutine Basis2momentum_cube(is2c, arr, title, spin)
    implicit none
    logical,intent(in)          :: is2c
    complex(dp),intent(in)      :: arr(2*cbdm)
    character(len=*),intent(in) :: title
    character(len=200)          :: env_TRESC
    character(len=100)          :: line
    integer,optional            :: spin
    integer                     :: end
    integer                     :: i, ii, jj, kk
    real(dp)                    :: maxexpo, maxcoeinshell, maxexpoinshell
    real(dp)                    :: maxx, minx, maxy, miny, maxz, minz, sigma
    real(dp)                    :: deltax, deltay, deltaz
    integer(8)                  :: nx, ny, nz
    !DIR$ ATTRIBUTES ALIGN:align_size :: mmts, datsa, datsb
    real(dp),allocatable        :: mmts(:,:)
    complex(dp),allocatable     :: datsa(:), datsb(:)
    !DIR$ ATTRIBUTES ALIGN:align_size :: dats
    complex(dp),allocatable     :: dats(:,:,:)

    call getenv('TRESC',env_TRESC)
    open(30, file=trim(env_TRESC)//'/gridsettings.ini', status='old', &
    action='read', iostat=ios)
    if (ios /= 0) call terminate(&
    "can't open file "//trim(env_TRESC)//'/gridsettings.ini')
    nx = 64        ! default settings
    ny = 64
    nz = 64
    sigma = 3.0_dp
    do
      read(30, '(A)', iostat=ios) line
      if (ios /= 0) then
        exit
      else if(index(line,'nx=') == 1) then
        read(line(index(line,'=')+1:index(line,';')-1), '(I)') nx
        if (nx <= 0) call terminate("nx should be positive integer")
      else if(index(line,'ny=') == 1) then
        read(line(index(line,'=')+1:index(line,';')-1), '(I)') ny
        if (ny <= 0) call terminate("ny should be positive integer")
      else if(index(line,'nz=') == 1) then
        read(line(index(line,'=')+1:index(line,';')-1), '(I)') nz
        if (nz <= 0) call terminate("nz should be positive integer")
      else if(index(line,'sigma=') == 1) then
        read(line(index(line,'=')+1:index(line,';')-1), '(F)') sigma
        if (sigma <= 0.0) call terminate("sigma should be positive float")
      end if
    end do
    close(30)
    allocate(dats(2*nz,ny,nx))
    allocate(mmts(nz,3), datsa(nz), datsb(nz))
    if (is2c) then
      open(99,file=trim(title)//"-mmtreal.cub",status="replace",action="write")
      write(99,'(A)') 'generate by TRESC: real part of 2-component orb '&
      //trim(title)//' in momentum space'
      write(99,'(A,I,A)') 'total ',nx*ny*nz,' grid points'
      open(100,file=trim(title)//"-mmtimg.cub",status="replace",action="write")
      write(100,'(A)') 'generated by TRESC: imaginary part of 2-component orb '&
      //trim(title)//' in momentum space'
      write(100,'(A,I,A)') 'total ',nx*ny*nz,' grid points'
    else
      open(99, file=trim(title)//"-mmt.cub", status="replace", action="write")
      write(99,'(A)') 'generated by TRESC: scalar orb '&
      //trim(title)//' in momentum space'
      write(99,'(A,I,A)') 'total ',nx*ny*nz,' points'
    end if
    ! generate grid points
    ! three-sigma boundary of the primitive GTO contributing most significantly
    ! to the GTO contributing more than 5% to the MO
    maxexpo = 0.0_dp
    do ii = 1, cbdm
      if (abs(arr(ii)) < 0.05) cycle
      maxcoeinshell = cbdata(ii)%Ncoe(1)
      maxexpoinshell = cbdata(ii)%expo(1)
      do jj = 1, cbdata(ii)%contr
        if (cbdata(ii)%Ncoe(jj) > maxcoeinshell) then
          maxcoeinshell = cbdata(ii)%Ncoe(jj)
          maxexpoinshell = cbdata(ii)%expo(jj)
        end if
      end do
      if (maxexpoinshell > maxexpo) maxexpo = maxexpoinshell
    end do
    ! amplitude of GTO under the momentum representation exhibits inversion
    ! symmetry, abs(psi(p)) = abs(psi(-p)), therefore, the square box is centred
    ! on (0,0,0).
    maxx = sigma*dsqrt(2.0_dp*maxexpo)
    minx = -maxx
    deltax = (maxx-minx) / real(nx,dp)
    maxy = maxx
    miny = minx
    deltay = deltax
    maxz = maxx
    minz = minx
    deltaz = deltax
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
    !$omp parallel num_threads(16) default(shared) private(i, ii, jj, kk, &
    !$omp& mmts, datsa, datsb) if(atom_count > 5)
    !$omp do schedule(static)
    do i = 1, nx
      ii = i
      mmts(:,1) = minx + (ii-1) * deltax
      do jj = 1, ny
        mmts(:,2) = miny + (jj-1) * deltay
        do kk = 1, nz
          mmts(kk,3) = minz + (kk-1) * deltaz
        end do
        if (is2c) then
          call Grid_momentum(arr, nz, mmts, 1, datsa)
          call Grid_momentum(arr, nz, mmts, 2, datsb)
          do kk = 1, nz
            dats(2*kk-1,jj,ii) = datsa(kk)
            dats(2*kk,jj,ii) = datsb(kk)
          end do
        else
          call Grid_momentum(arr, nz, mmts, spin, dats(1:nz,jj,ii))
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
    deallocate(dats, mmts, datsa, datsb)
  end subroutine Basis2momentum_cube

!------------------------------------------------------------
!> transform states from basis repre to real space repre (Becke's fuzzy grid)
!!
!! it's used to generate xc energies for xc functional
!!
!! returns only the integral value
  subroutine Basis2real_Becke_XCintegral(rho_m, ex, ec, Fockx, Fockc)
    implicit none
    real(dp),intent(out)   :: ex         ! exchange functional
    real(dp), optional     :: ec         ! correlation functional
    complex(dp),intent(in) :: rho_m(2*cbdm, 2*cbdm)! Cartesian density matrix
    integer(8)             :: i          ! loop variable
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
  end subroutine Basis2real_Becke_XCintegral

!------------------------------------------------------------
!> assign radial grid points (Gauss weight) using 2nd Chebyshev method
!!
!! and combined with the Lebedev sphere grid to calculate single-centre integral
  pure subroutine Chebyshev2_Lebedev(count, nr, nl, pos3w1)
    implicit none
    integer(8),intent(in) :: count, nr, nl
    integer               :: ii, jj
    integer               :: nl2   ! nl2 = nl
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
!> grid points in real space of input (basis repre) vector
!!
!! spin=0:total =1:alpha =2:beta
  pure subroutine Grid_real(arr, n, points, spin, dats)
    implicit none
    complex(dp),intent(in)  :: arr(:)
    integer(8),intent(in)   :: n
    real(dp),intent(in)     :: points(n, 3)
    integer,intent(in)      :: spin
    complex(dp),intent(out) :: dats(n)
    integer                 :: ii, jj, kk
    integer                 :: contr
    real(dp)                :: vec(n, 3)
    integer                 :: L, M
    integer                 :: fac(3)
    real(dp)                :: coeff
    real(dp)                :: expo
    real(dp)                :: val(n)
    integer                 :: init, final, step
    dats = c0
    if (spin == 0) then
      init = 1
      step = 1
      final = 2*cbdm
    else if (spin == 1) then
      init = 1
      step = 1
      final = cbdm
    else if (spin == 2) then
      init = cbdm+1
      step = 1
      final = 2*cbdm
    end if
    do ii = init, final, step
      if (ii <= cbdm) then
        kk = ii
      else
        kk = ii - cbdm
      end if
      val = 0.0_dp
      vec(:,1) = points(:,1) - cbdata(kk) % pos(1)
      vec(:,2) = points(:,2) - cbdata(kk) % pos(2)
      vec(:,3) = points(:,3) - cbdata(kk) % pos(3)
      contr    = cbdata(kk) % contr
      L        = cbdata(kk) % L
      M        = cbdata(kk) % M
      fac(:)   = AO_fac(:,L,M)
      do jj = 1, contr
        expo   = cbdata(kk) % expo(jj)
        coeff  = cbdata(kk) % Ncoe(jj)
        val(:) = val(:) + coeff * exp(-expo*(sum(vec(:,:)**2,dim=2)))
      end do
      val(:) = val(:) * (vec(:,1)**fac(1))
      val(:) = val(:) * (vec(:,2)**fac(2))
      val(:) = val(:) * (vec(:,3)**fac(3))
      dats   = dats + val * arr(ii)
    end do
  end subroutine Grid_real

!------------------------------------------------------------
!> grid points in momentum space of input (basis repre) vector
!!
!! spin=0:total =1:alpha =2:beta
  pure subroutine Grid_momentum(arr, n, mmts, spin, dats)
    implicit none
    complex(dp),intent(in)  :: arr(:)
    integer(8),intent(in)   :: n
    real(dp), intent(in)    :: mmts(n,3)
    integer,intent(in)      :: spin
    complex(dp),intent(out) :: dats(n)
    integer                 :: ii, jj, kk, mm
    integer                 :: contr
    real(dp)                :: rx, ry, rz
    real(dp)                :: px, py, pz
    integer                 :: L, M
    integer                 :: fac(3)
    complex(dp)             :: coeff(cbdm,16)
    real(dp)                :: expo(16)
    complex(dp)             :: val
    integer                 :: init, final, step
    real(dp)                :: p2, inv2sqrtb(cbdm,16)

    do ii = 1, cbdm
      inv2sqrtb(ii,:) = 1.0_dp / (2.0_dp*sqrt(cbdata(ii)%expo))
      L = cbdata(ii)%L - 1
      expo = cbdata(ii)%expo
      coeff(ii,:) = cbdata(ii)%Ncoe(:)*(-ci)**L/((2.0_dp**(real(L,dp)+1.5_dp))*&
      (expo(:)**(0.5_dp*real(L,dp)+1.5_dp)))
    end do
    if (spin == 0) then
      init = 1
      step = 1
      final = 2*cbdm
    else if (spin == 1) then
      init = 1
      step = 1
      final = cbdm
    else if (spin == 2) then
      init = cbdm+1
      step = 1
      final = 2*cbdm
    end if
    dats = c0
    do mm = 1, n
      px = mmts(mm,1)
      py = mmts(mm,2)
      pz = mmts(mm,3)
      p2 = px**2 + py**2 + pz**2
      do ii = init, final, step
        if (ii <= cbdm) then
          kk = ii
        else
          kk = ii - cbdm
        end if
        val = c0
        contr  = cbdata(kk)%contr
        expo   = cbdata(kk)%expo
        fac(:) = AO_fac(:,cbdata(kk)%L,cbdata(kk)%M)
        rx     = cbdata(kk)%pos(1)
        ry     = cbdata(kk)%pos(2)
        rz     = cbdata(kk)%pos(3)
        do jj = 1, contr
          val = val + coeff(kk,jj)                       *&
          exp(-p2/(4.0_dp*expo(jj)))                     *&
          Hermite_poly(fac(1),px*inv2sqrtb(kk,jj))       *&
          Hermite_poly(fac(2),py*inv2sqrtb(kk,jj))       *&
          Hermite_poly(fac(3),pz*inv2sqrtb(kk,jj))
        end do
        val = val * exp(-ci*(px*rx+py*ry+pz*rz))
        dats(mm) = dats(mm) + val * arr(ii)
      end do
    end do
  end subroutine Grid_momentum

!------------------------------------------------------------
!> reference frame transformation in real space of mol
!!
!! assign molecule in frame of motion (trafomol)
  subroutine Mol_real_Lorentztrafo()
    implicit none
    integer  :: ii, jj, kk
    real(dp) :: trafocoe, dot
    ! Lorentz transformation of x, leads to length contraction
    ! x' = x - gamma/(gamma+1.0_dp)(x.beta) * beta
    trafomol = mol
    trafocoe = -gamma/(gamma+1.0_dp)
    do ii = 1, atom_count
      dot = trafocoe*sum(beta(:)*mol(ii)%pos(:))
      trafomol(ii)%pos(:) = mol(ii)%pos(:) + dot*beta(:)
    end do
  end subroutine Mol_real_Lorentztrafo

!------------------------------------------------------------
!> reference frame transformation in real space of input (basis repre) vector
!!
!! input "points" is the coordinates in rest frame
!!
!! “simultaneous” in motion frame: x_rest = L(x_motion)|t_motion=0,
!! in this case, x' = x_rest; x = x_motion
  pure subroutine Grid_real_Lorentztrafo(arr,n,points,trafopoints,datsa,datsb)
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
      vec(:,1) = points(:,1) - cbdata(kk) % pos(1)
      vec(:,2) = points(:,2) - cbdata(kk) % pos(2)
      vec(:,3) = points(:,3) - cbdata(kk) % pos(3)
      contr    = cbdata(kk) % contr
      L        = cbdata(kk) % L
      M        = cbdata(kk) % M
      fac(:)   = AO_fac(:,L,M)
      do jj = 1, contr
        expo   = cbdata(kk) % expo(jj)
        coeff  = cbdata(kk) % Ncoe(jj)
        val(:) = val(:) + coeff * exp(-expo*(sum(vec(:,:)**2,dim=2)))
      end do
      val(:) = val(:) * (vec(:,1)**fac(1))
      val(:) = val(:) * (vec(:,2)**fac(2))
      val(:) = val(:) * (vec(:,3)**fac(3))
      datsa  = datsa + val * arr(kk)
      datsb  = datsb + val * arr(kk+cbdm)
    end do

    ! Imitative transformation of s, causes s_z to away from +/- 1/2,
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
  end subroutine Grid_real_Lorentztrafo

!------------------------------------------------------------
!> perform SU(2) Euler transformation of input (spin part of) orbital(s)
!!
!! in spinor space, Euler angle: phi [0,2pi), theta [0,pi], chi [0,2pi)
!!
!! based on z-y-z convention: 
!! exp(-(i/2)*phi*sigma_z) * exp(-(i/2)*theta*sigma_y) * exp(-(i/2)*chi*sigma_z)
  pure subroutine MO_spinor_SU2trafo(orb_i, orb_o, phi, theta, chi)
    implicit none
    complex(dp),intent(in)                :: orb_i(:,:) ! input orbitals
    complex(dp),intent(out)               :: orb_o(:,:) ! output orbitals
    real(dp),intent(in)                   :: phi        ! 1st Euler angle
    real(dp),intent(in)                   :: theta      ! 2nd Euler angle
    real(dp),intent(in)                   :: chi        ! 3rd Euler angle
    integer                               :: dm(2)
    integer                               :: ii, jj     ! loop variables
    complex(dp)                           :: SU2rot(2,2)! SU(2) rotation matirx
    dm = shape(orb_i)
    orb_o = c0
    ! S_phi_theta_chi = 
    !  _                                                                  _
    ! |                                                                    |
    ! | exp(-i/2*(phi+chi))cos(theta/2),  -exp(-i/2*(phi-chi))sin(theta/2) |
    ! | exp(i/2*(phi-chi))sin(theta/2) ,   exp(i/2*(phi+chi))cos(theta/2)  |
    ! |_                                                                  _|
    SU2rot(1,1) = exp(-0.5_dp*ci*(phi+chi))*dcos(0.5_dp*theta)
    SU2rot(1,2) = -exp(-0.5_dp*ci*(phi-chi))*dsin(0.5_dp*theta)
    SU2rot(2,1) = exp(0.5_dp*ci*(phi-chi))*dsin(0.5_dp*theta)
    SU2rot(2,2) = exp(0.5_dp*ci*(phi+chi))*dcos(0.5_dp*theta)
    do ii = 1, dm(2)
      do jj = 1, dm(1)/2
        ! alpha component after rotation
        orb_o(jj,ii) = SU2rot(1,1)*orb_i(jj,ii) + &
        SU2rot(1,2)*orb_i(dm(1)/2+jj,ii)
        ! beta componnet after rotation
        orb_o(dm(1)/2+jj,ii) = SU2rot(2,1)*orb_i(jj,ii) + &
        SU2rot(2,2)*orb_i(dm(1)/2+jj,ii)
      end do
    end do
  end subroutine MO_spinor_SU2trafo

end module Representation