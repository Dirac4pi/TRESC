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
    complex(dp)                 :: datsa(59000), datsb(59000)
    real(dp)                    :: pos3(590,3,100), trafopos3(590,3,100)
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
    integer                     :: i, ii, jj, kk, ll
    real(dp)                    :: maxx, minx, maxy, miny, maxz, minz
    real(dp)                    :: grid, edge
    logical                     :: autoedge
    real(dp)                    :: autoedgepos(54,3), temp(3)
    complex(dp)                 :: autoedgedat(54)
    real(dp)                    :: autoedgeavr
    integer(8)                  :: nx, ny, nz
    !DIR$ ATTRIBUTES ALIGN:align_size :: expx, expy, expz, val
    real(dp),allocatable        :: expx(:,:,:), expy(:,:,:), expz(:,:,:)
    real(dp),allocatable        :: facx(:,:), facy(:,:), facz(:,:)
    real(dp),allocatable        :: val(:)
    real(dp)                    :: coe, vec(cbdm), fac(3)
    real(dp)                    :: xx, yy, zz
    !DIR$ ATTRIBUTES ALIGN:align_size :: datsa, datsb, tldats
    complex(dp),allocatable     :: datsa(:,:,:), datsb(:,:,:), tldats(:,:,:)

    call getenv('TRESC',env_TRESC)
    open(30, file=trim(env_TRESC)//'/gridsettings.ini', status='old', &
    action='read', iostat=ios)
    if (ios /= 0) call terminate(&
    "can't open file "//trim(env_TRESC)//'/gridsettings.ini')
    grid = 0.15_dp
    edge = 5.0_dp
    autoedge = .false.
    do
      read(30, '(A)', iostat=ios) line
      if (ios /= 0) then
        exit
      else if(index(line,'rgrid=') == 1) then
        read(line(index(line,'=')+1:index(line,';')-1), '(F)') grid
        if (grid <= 0.0) call terminate("grid should be positive float")
      else if(index(line,'edge=') == 1) then
        if (line(index(line,'=')+1:index(line,';')-1) == 'auto') then
          autoedge = .true.
        else
          read(line(index(line,'=')+1:index(line,';')-1), '(F)') edge
          if (edge <= 0.0) call terminate("edge should be positive float")
        end if
      end if
    end do
    close(30)
    if (autoedge) then
      edge = 5.0
      temp = (/0.25, 0.5, 0.75/)
      do ii = 1, 100
        maxx = maxval(mol(:)%pos(1)) + edge
        minx = minval(mol(:)%pos(1)) - edge
        nx = int((maxx-minx) / grid) + 1
        maxy = maxval(mol(:)%pos(2)) + edge
        miny = minval(mol(:)%pos(2)) - edge
        ny = int((maxy-miny) / grid) + 1
        maxz = maxval(mol(:)%pos(3)) + edge
        minz = minval(mol(:)%pos(3)) - edge
        nz = int((maxz-minz) / grid) + 1
        autoedgepos(1:9,1) = maxx
          autoedgepos(1:3,2) = miny + 0.25_dp*(maxy-miny)
          autoedgepos(1:3,3) = minz + (maxz-minz)*temp(:)
          autoedgepos(4:6,2) = miny + 0.5_dp*(maxy-miny)
          autoedgepos(4:6,3) = minz + (maxz-minz)*temp(:)
          autoedgepos(7:9,2) = miny + 0.75_dp*(maxy-miny)
          autoedgepos(7:9,3) = minz + (maxz-minz)*temp(:)
        autoedgepos(10:18,1) = minx
          autoedgepos(10:12,2) = miny + 0.25_dp*(maxy-miny)
          autoedgepos(10:12,3) = minz + (maxz-minz)*temp(:)
          autoedgepos(13:15,2) = miny + 0.5_dp*(maxy-miny)
          autoedgepos(13:15,3) = minz + (maxz-minz)*temp(:)
          autoedgepos(16:18,2) = miny + 0.75_dp*(maxy-miny)
          autoedgepos(16:18,3) = minz + (maxz-minz)*temp(:)
        autoedgepos(19:27,2) = maxy
          autoedgepos(19:21,1) = minx + 0.25_dp*(maxx-minx)
          autoedgepos(19:21,3) = minz + (maxz-minz)*temp(:)
          autoedgepos(22:24,1) = minx + 0.5_dp*(maxx-minx)
          autoedgepos(22:24,3) = minz + (maxz-minz)*temp(:)
          autoedgepos(25:27,1) = minx + 0.75_dp*(maxx-minx)
          autoedgepos(25:27,3) = minz + (maxz-minz)*temp(:)
        autoedgepos(28:36,2) = miny
          autoedgepos(28:30,1) = minx + 0.25_dp*(maxx-minx)
          autoedgepos(28:30,3) = minz + (maxz-minz)*temp(:)
          autoedgepos(31:33,1) = minx + 0.5_dp*(maxx-minx)
          autoedgepos(31:33,3) = minz + (maxz-minz)*temp(:)
          autoedgepos(34:36,1) = minx + 0.75_dp*(maxx-minx)
          autoedgepos(34:36,3) = minz + (maxz-minz)*temp(:)
        autoedgepos(37:45,3) = maxz
          autoedgepos(37:39,1) = minx + 0.25_dp*(maxx-minx)
          autoedgepos(37:39,2) = miny + (maxy-miny)*temp(:)
          autoedgepos(40:42,1) = minx + 0.5_dp*(maxx-minx)
          autoedgepos(40:42,2) = miny + (maxy-miny)*temp(:)
          autoedgepos(43:45,1) = minx + 0.75_dp*(maxx-minx)
          autoedgepos(43:45,2) = miny + (maxy-miny)*temp(:)
        autoedgepos(46:54,3) = minz
          autoedgepos(46:48,1) = minx + 0.25_dp*(maxx-minx)
          autoedgepos(46:48,2) = miny + (maxy-miny)*temp(:)
          autoedgepos(49:51,1) = minx + 0.5_dp*(maxx-minx)
          autoedgepos(49:51,2) = miny + (maxy-miny)*temp(:)
          autoedgepos(52:54,1) = minx + 0.75_dp*(maxx-minx)
          autoedgepos(52:54,2) = miny + (maxy-miny)*temp(:)
        call Grid_real(arr, 54, autoedgepos, 0, autoedgedat)
        autoedgeavr = sum(abs(autoedgedat(:))) / 54.0
        ! 1E-3 is enough for 0.05 isosurface
        if (autoedgeavr > 1E-5 .and. autoedgeavr < 1E-3) then
          exit
        else if (autoedgeavr >= 1E-3) then
          edge = edge * 1.2
        else
          edge = edge * 0.8
        end if
      end do
      if (ii == 100) then
        call terminate(&
        'Basis2real_cube: autoedge failed, please set edge manualy')
      else
        write(*,*) ' edge tuned to ', edge
      end if
    end if
    ! generate grid points
    maxx = maxval(mol(:)%pos(1)) + edge
    minx = minval(mol(:)%pos(1)) - edge
    nx = int((maxx-minx) / grid) + 1
    maxy = maxval(mol(:)%pos(2)) + edge
    miny = minval(mol(:)%pos(2)) - edge
    ny = int((maxy-miny) / grid) + 1
    maxz = maxval(mol(:)%pos(3)) + edge
    minz = minval(mol(:)%pos(3)) - edge
    nz = int((maxz-minz) / grid) + 1
    write(*,'(A,I4,A,I4,A,I4)') '  grid ', nx, '*', ny, '*', nz
    if (nx*ny*nz > 1E7) then
      call terminate("nx*ny*nz too large, increase `rgrid` or decrease `edge`")
    else if (nx*ny*nz < 1E3) then
      write(*,*) "Warning: nx*ny*nz too small, may lead to inaccurate plotting."
    end if
    allocate(datsa(nz,ny,nx), datsb(nz,ny,nx), source=c0)
    allocate(val(nz))
    allocate(expx(nx,cbdm,16),expy(ny,cbdm,16),expz(nz,cbdm,16))
    allocate(facx(nx,cbdm),facy(ny,cbdm),facz(nz,cbdm))
    if (is2c) then
      open(99,file=trim(title)//"-rar.cub",status="replace",action="write")
      write(99,'(A)') &
      'Real Space Uniform Grid:Alpha Real part of Spinor Orb '//trim(title)
      write(99,'(A,I,A)') 'total ',nx*ny*nz,' grid points'
      open(100,file=trim(title)//"-rai.cub",status="replace",action="write")
      write(100,'(A)') &
      'Real Space Uniform Grid:Alpha Imaginary part of Spinor Orb '//trim(title)
      write(100,'(A,I,A)') 'total ',nx*ny*nz,' grid points'
      open(101,file=trim(title)//"-rbr.cub",status="replace",action="write")
      write(101,'(A)') &
      'Real Space Uniform Grid:Beta Real part of Spinor Orb '//trim(title)
      write(101,'(A,I,A)') 'total ',nx*ny*nz,' grid points'
      open(102,file=trim(title)//"-rbi.cub",status="replace",action="write")
      write(102,'(A)') &
      'Real Space Uniform Grid:Beta Imaginary part of Spinor Orb '//trim(title)
      write(102,'(A,I,A)') 'total ',nx*ny*nz,' grid points'
    else
      open(99, file=trim(title)//"-real.cub", status="replace", action="write")
      write(99,'(A)') &
      'Real Space Uniform Grid:Real part of Scalar Orb '//trim(title)
      write(99,'(A,I,A)') 'total ',nx*ny*nz,' points'
    end if
    if (is2c) then
      write(99,'(I4,3F14.8)') atom_count, minx, miny, minz
      write(99,'(I3,3F14.8)') nx, grid, 0.0, 0.0
      write(99,'(I3,3F14.8)') ny, 0.0, grid, 0.0
      write(99,'(I3,3F14.8)') nz, 0.0, 0.0, grid
      write(100,'(I4,3F14.8)') atom_count, minx, miny, minz
      write(100,'(I3,3F14.8)') nx, grid, 0.0, 0.0
      write(100,'(I3,3F14.8)') ny, 0.0, grid, 0.0
      write(100,'(I3,3F14.8)') nz, 0.0, 0.0, grid
      write(101,'(I4,3F14.8)') atom_count, minx, miny, minz
      write(101,'(I3,3F14.8)') nx, grid, 0.0, 0.0
      write(101,'(I3,3F14.8)') ny, 0.0, grid, 0.0
      write(101,'(I3,3F14.8)') nz, 0.0, 0.0, grid
      write(102,'(I4,3F14.8)') atom_count, minx, miny, minz
      write(102,'(I3,3F14.8)') nx, grid, 0.0, 0.0
      write(102,'(I3,3F14.8)') ny, 0.0, grid, 0.0
      write(102,'(I3,3F14.8)') nz, 0.0, 0.0, grid
      do ii = 1, atom_count
        write(99,'(I4,4F14.8)') mol(ii)%atom_number, real(mol(ii)%atom_number),&
        mol(ii)%pos(1), mol(ii)%pos(2), mol(ii)%pos(3)
        write(100,'(I4,4F14.8)')mol(ii)%atom_number, real(mol(ii)%atom_number),&
        mol(ii)%pos(1), mol(ii)%pos(2), mol(ii)%pos(3)
        write(101,'(I4,4F14.8)')mol(ii)%atom_number, real(mol(ii)%atom_number),&
        mol(ii)%pos(1), mol(ii)%pos(2), mol(ii)%pos(3)
        write(102,'(I4,4F14.8)')mol(ii)%atom_number, real(mol(ii)%atom_number),&
        mol(ii)%pos(1), mol(ii)%pos(2), mol(ii)%pos(3)
      end do
    else
      write(99,'(I4,3F14.8)') atom_count, minx, miny, minz
      write(99,'(I3,3F14.8)') nx, grid, 0.0, 0.0
      write(99,'(I3,3F14.8)') ny, 0.0, grid, 0.0
      write(99,'(I3,3F14.8)') nz, 0.0, 0.0, grid
      do ii = 1, atom_count
        write(99,'(I4,4F14.8)') mol(ii)%atom_number, real(mol(ii)%atom_number),&
        mol(ii)%pos(1), mol(ii)%pos(2), mol(ii)%pos(3)
      end do
      do ii = 1, len(title)
        if (is_alpha(title(ii:ii))) exit
      end do
    end if
    do ii = 1, nx
      xx = minx + (ii-1) * grid
      do jj = 1, cbdm
        vec(jj) = xx - cbdata(jj)%pos(1)
        fac(:) = cbdata(jj)%fac
        facx(ii,jj) = vec(jj) ** fac(1)
        do kk = 1, cbdata(jj)%contr
          expx(ii,jj,kk) = dexp(-cbdata(jj)%expo(kk)*vec(jj)*vec(jj))
        end do
      end do
    end do
    do ii = 1, ny
      yy = miny + (ii-1) * grid
      do jj = 1, cbdm
        vec(jj) = yy - cbdata(jj)%pos(2)
        fac(:) = cbdata(jj)%fac
        facy(ii,jj) = vec(jj) ** fac(2)
        do kk = 1, cbdata(jj)%contr
          expy(ii,jj,kk) = dexp(-cbdata(jj)%expo(kk)*vec(jj)*vec(jj))
        end do
      end do
    end do
    do ii = 1, nz
      zz = minz + (ii-1) * grid
      do jj = 1, cbdm
        vec(jj) = zz - cbdata(jj)%pos(3)
        fac(:) = cbdata(jj)%fac
        facz(ii,jj) = vec(jj) ** fac(3)
        do kk = 1, cbdata(jj)%contr
          expz(ii,jj,kk) = dexp(-cbdata(jj)%expo(kk)*vec(jj)*vec(jj))
        end do
      end do
    end do
    !$omp parallel num_threads(16) default(shared) private(i, ii, jj, kk, ll, &
    !$omp& coe, val, tldats)
    allocate(tldats(2*nz,ny,nx), source=c0)
    !$omp do schedule(static)
    do i = 1, nx
      ii = i
      do jj = 1, ny
        if (is2c) then
          do kk = 1, cbdm
            val = 0.0_dp
            do ll = 1, cbdata(kk)%contr
              coe = cbdata(kk) % Ncoe(ll)
              val(:) = val(:) + coe*expx(ii,kk,ll)*expy(jj,kk,ll)*expz(:,kk,ll)
            end do
            val(:) = val(:) * facx(ii,kk)*facy(jj,kk)*facz(:,kk)
            tldats(1:nz,jj,ii) = tldats(1:nz,jj,ii) + val*arr(kk)
            tldats(nz+1:2*nz,jj,ii) = tldats(nz+1:2*nz,jj,ii) + val*arr(cbdm+kk)
          end do
        else
          do kk = 1, cbdm
            val = 0.0_dp
            do ll = 1, cbdata(kk)%contr
              coe = cbdata(kk) % Ncoe(ll)
              val(:) = val(:) + coe*expx(ii,kk,ll)*expy(jj,kk,ll)*expz(:,kk,ll)
            end do
            val(:) = val(:) * facx(ii,kk)*facy(jj,kk)*facz(:,kk)
            if (spin==1) tldats(1:nz,jj,ii)=tldats(1:nz,jj,ii)+val*arr(kk)
            if (spin==2) tldats(1:nz,jj,ii)=tldats(1:nz,jj,ii)+val*arr(cbdm+kk)
          end do
        end if
      end do
    end do
    !$omp end do
    !$omp critical
    if (is2c) then
      datsa(:,:,:) = datsa(:,:,:) + tldats(1:nz,:,:)
      datsb(:,:,:) = datsb(:,:,:) + tldats(nz+1:2*nz,:,:)
    else
      datsa(:,:,:) = datsa(:,:,:) + tldats(1:nz,:,:)
    end if
    !$omp end critical
    deallocate(tldats)
    !$omp end parallel
    do ii = 1, nx
      do jj = 1, ny
        if (is2c) then
          do kk = 1, nz, 6
            end = min(kk+5, 2*nz)
            write(99,'(6ES15.6)') real(datsa(kk:end,jj,ii))
            write(100,'(6ES15.6)') aimag(datsa(kk:end,jj,ii))
            write(101,'(6ES15.6)') real(datsb(kk:end,jj,ii))
            write(102,'(6ES15.6)') aimag(datsb(kk:end,jj,ii))
          end do
        else
          do kk = 1, nz, 6
            end = min(kk+5, nz)
            write(99,'(6ES15.6)') real(datsa(kk:end,jj,ii))
          end do
        end if
      end do
    end do
    close(99)
    if (is2c) then
      close(100)
      close(101)
      close(102)
    end if
    deallocate(datsa,datsb,val)
    deallocate(expx,expy,expz)
    deallocate(facx,facy,facz)
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
    real(dp)                    :: arrabs(2*cbdm)
    integer                     :: arrindex(2*cbdm)
    integer                     :: arr97(cbdm), n95
    real(dp)                    :: sum95, sum100
    character(len=*),intent(in) :: title
    character(len=200)          :: env_TRESC
    character(len=100)          :: line
    integer,optional            :: spin
    integer                     :: end
    integer                     :: i, ii, jj, kk, ll, mm
    real(dp)                    :: maxexpo, maxcoeinshell, maxexpoinshell
    real(dp)                    :: maxx, minx, maxy, miny, maxz, minz
    real(dp)                    :: grid, sigma
    logical                     :: autosigma
    real(dp)                    :: autosigmammt(54,3), temp(3)
    complex(dp)                 :: autosigmadat(54)
    real(dp)                    :: autosigmaavr
    integer(8)                  :: nx, ny, nz
    !DIR$ ATTRIBUTES ALIGN:align_size :: Hermitex, Hermitey, Hermitez
    real(dp),allocatable        :: Hermitex(:,:,:)
    real(dp),allocatable        :: Hermitey(:,:,:)
    real(dp),allocatable        :: Hermitez(:,:,:)
    real(dp)                    :: rx, ry, rz
    real(dp)                    :: px, py, pz
    integer                     :: L, M
    integer                     :: fac(3)
    real(dp)                    :: expo(16)
    real(dp)                    :: p2, coe(cbdm,16), inv2sqrtb(cbdm,16)
    complex(dp)                 :: val
    !DIR$ ATTRIBUTES ALIGN:align_size :: datsa, datsb, tldats
    complex(dp),allocatable     :: datsa(:,:,:), datsb(:,:,:), tldats(:,:,:)

    call getenv('TRESC',env_TRESC)
    open(30, file=trim(env_TRESC)//'/gridsettings.ini', status='old', &
    action='read', iostat=ios)
    if (ios /= 0) call terminate(&
    "can't open file "//trim(env_TRESC)//'/gridsettings.ini')
    grid = 5.0_dp
    sigma = 3.0_dp
    autosigma = .false.
    do
      read(30, '(A)', iostat=ios) line
      if (ios /= 0) then
        exit
      else if(index(line,'mgrid=') == 1) then
        read(line(index(line,'=')+1:index(line,';')-1), '(F)') grid
        if (grid <= 0) call terminate("grid should be positive float")
      else if(index(line,'sigma=') == 1) then
        if (line(index(line,'=')+1:index(line,';')-1) == 'auto') then
          autosigma = .true.
        else
          read(line(index(line,'=')+1:index(line,';')-1), '(F)') sigma
          if (sigma <= 0.0) call terminate("sigma should be positive float")
        end if
      end if
    end do
    close(30)
    ! generate grid points
    ! three-sigma boundary of the primitive GTO with the largest exponent in the
    ! smallest set of GTOs that total contribution to the MO is more than 95%.
    sum100 = sum(abs(arr(:))**2)
    call sort(arr, arrabs, arrindex)
    sum95 = 0.0_dp
    n95 = 0
    do while (sum95/sum100 < 0.99)
      n95 = n95 + 1
      sum95 = sum95 + arrabs(2*cbdm-n95+1)**2
    end do
    maxexpo = 0.0_dp
    do ii = 1, n95
      ll = arrindex(2*cbdm-ii+1)
      if (ll > cbdm) ll = ll - cbdm
      maxcoeinshell = abs(cbdata(ll)%Ncoe(1))
      maxexpoinshell = cbdata(ll)%expo(1)
      do jj = 1, cbdata(ll)%contr
        if (abs(cbdata(ll)%Ncoe(jj)) > maxcoeinshell) then
          maxcoeinshell = abs(cbdata(ll)%Ncoe(jj))
          maxexpoinshell = cbdata(ll)%expo(jj)
        end if
      end do
      if (maxexpoinshell > maxexpo) maxexpo = maxexpoinshell
    end do
    ! amplitude of GTO under the momentum representation exhibits inversion
    ! symmetry, abs(psi(p)) = abs(psi(-p)), therefore, the square box is centred
    ! on (0,0,0).
    if (autosigma) then
      sigma = 1.5
      temp = (/0.25, 0.5, 0.75/)
      do ii = 1, 100
        maxx = sigma*dsqrt(2.0_dp*maxexpo)
        minx = -maxx
        nx = int((maxx-minx) / grid) + 1
        maxy = maxx
        miny = minx
        ny = nx
        maxz = maxx
        minz = minx
        nz = nx
        autosigmammt(1:9,1) = maxx
          autosigmammt(1:3,2) = miny + 0.25_dp*(maxy-miny)
          autosigmammt(1:3,3) = minz + (maxz-minz)*temp(:)
          autosigmammt(4:6,2) = miny + 0.5_dp*(maxy-miny)
          autosigmammt(4:6,3) = minz + (maxz-minz)*temp(:)
          autosigmammt(7:9,2) = miny + 0.75_dp*(maxy-miny)
          autosigmammt(7:9,3) = minz + (maxz-minz)*temp(:)
        autosigmammt(10:18,1) = minx
          autosigmammt(10:12,2) = miny + 0.25_dp*(maxy-miny)
          autosigmammt(10:12,3) = minz + (maxz-minz)*temp(:)
          autosigmammt(13:15,2) = miny + 0.5_dp*(maxy-miny)
          autosigmammt(13:15,3) = minz + (maxz-minz)*temp(:)
          autosigmammt(16:18,2) = miny + 0.75_dp*(maxy-miny)
          autosigmammt(16:18,3) = minz + (maxz-minz)*temp(:)
        autosigmammt(19:27,2) = maxy
          autosigmammt(19:21,1) = minx + 0.25_dp*(maxx-minx)
          autosigmammt(19:21,3) = minz + (maxz-minz)*temp(:)
          autosigmammt(22:24,1) = minx + 0.5_dp*(maxx-minx)
          autosigmammt(22:24,3) = minz + (maxz-minz)*temp(:)
          autosigmammt(25:27,1) = minx + 0.75_dp*(maxx-minx)
          autosigmammt(25:27,3) = minz + (maxz-minz)*temp(:)
        autosigmammt(28:36,2) = miny
          autosigmammt(28:30,1) = minx + 0.25_dp*(maxx-minx)
          autosigmammt(28:30,3) = minz + (maxz-minz)*temp(:)
          autosigmammt(31:33,1) = minx + 0.5_dp*(maxx-minx)
          autosigmammt(31:33,3) = minz + (maxz-minz)*temp(:)
          autosigmammt(34:36,1) = minx + 0.75_dp*(maxx-minx)
          autosigmammt(34:36,3) = minz + (maxz-minz)*temp(:)
        autosigmammt(37:45,3) = maxz
          autosigmammt(37:39,1) = minx + 0.25_dp*(maxx-minx)
          autosigmammt(37:39,2) = miny + (maxy-miny)*temp(:)
          autosigmammt(40:42,1) = minx + 0.5_dp*(maxx-minx)
          autosigmammt(40:42,2) = miny + (maxy-miny)*temp(:)
          autosigmammt(43:45,1) = minx + 0.75_dp*(maxx-minx)
          autosigmammt(43:45,2) = miny + (maxy-miny)*temp(:)
        autosigmammt(46:54,3) = minz
          autosigmammt(46:48,1) = minx + 0.25_dp*(maxx-minx)
          autosigmammt(46:48,2) = miny + (maxy-miny)*temp(:)
          autosigmammt(49:51,1) = minx + 0.5_dp*(maxx-minx)
          autosigmammt(49:51,2) = miny + (maxy-miny)*temp(:)
          autosigmammt(52:54,1) = minx + 0.75_dp*(maxx-minx)
          autosigmammt(52:54,2) = miny + (maxy-miny)*temp(:)
        call Grid_momentum(arr, 54, autosigmammt, 0, autosigmadat)
        autosigmaavr = sum(abs(autosigmadat(:))) / 54.0
        ! 1E-3 is enough for 0.05 isosurface
        if (autosigmaavr > 1E-5 .and. autosigmaavr < 1E-3) then
          exit
        else if (autosigmaavr >= 1E-3) then
          sigma = sigma * 1.2
        else
          sigma = sigma * 0.8
        end if
      end do
      if (ii == 100) then
        call terminate(&
        'Basis2momentum_cube: autosigma failed, please set sigma manualy')
      else
        write(*,*) ' sigma tuned to ', sigma
      end if
    end if
    maxx = sigma*dsqrt(2.0_dp*maxexpo)
    minx = -maxx
    nx = int((maxx-minx) / grid) + 1
    maxy = maxx
    miny = minx
    ny = nx
    maxz = maxx
    minz = minx
    nz = nx
    write(*,'(A,I4,A,I4,A,I4)') '  grid ', nx, '*', ny, '*', nz
    if (nx*ny*nz > 1E7) then
      call terminate("nx*ny*nz too large: increase `mgrid` or decrease `sigma`")
    else if (nx*ny*nz < 5E3) then
      write(*,*) "Warning: nx*ny*nz too small, may lead to inaccurate plotting."
    end if
    allocate(datsa(nz,ny,nx), datsb(nz,ny,nx))
    allocate(Hermitex(nx,cbdm,16),Hermitey(ny,cbdm,16),Hermitez(nz,cbdm,16))
    if (is2c) then
      open(99,file=trim(title)//"-mar.cub",status="replace",action="write")
      write(99,'(A)') &
      'Momentum Uniform Grid:Alpha Real part of Spinor Orb '//trim(title)
      write(99,'(A,I,A)') 'total ',nx*ny*nz,' grid points'
      open(100,file=trim(title)//"-mai.cub",status="replace",action="write")
      write(100,'(A)') &
      'Momentum Uniform Grid:Alpha Imaginary part of Spinor Orb '//trim(title)
      write(100,'(A,I,A)') 'total ',nx*ny*nz,' grid points'
      open(101,file=trim(title)//"-mbr.cub",status="replace",action="write")
      write(101,'(A)') &
      'Momentum Uniform Grid:Beta Real part of Spinor Orb '//trim(title)
      write(101,'(A,I,A)') 'total ',nx*ny*nz,' grid points'
      open(102,file=trim(title)//"-mbi.cub",status="replace",action="write")
      write(102,'(A)') &
      'Momentum Uniform Grid:Beta Imaginary part of Spinor Orb '//trim(title)
      write(102,'(A,I,A)') 'total ',nx*ny*nz,' grid points'
    else
      open(99, file=trim(title)//"-mmt.cub", status="replace", action="write")
      write(99,'(A)') &
      'Momentum Uniform Grid:Real part of Scalar Orb '//trim(title)
      write(99,'(A,I,A)') 'total ',nx*ny*nz,' points'
    end if
    if (is2c) then
      write(99,'(I4,3F14.8)') atom_count, minx, miny, minz
      write(99,'(I3,3F14.8)') nx, grid, 0.0, 0.0
      write(99,'(I3,3F14.8)') ny, 0.0, grid, 0.0
      write(99,'(I3,3F14.8)') nz, 0.0, 0.0, grid
      write(100,'(I4,3F14.8)') atom_count, minx, miny, minz
      write(100,'(I3,3F14.8)') nx, grid, 0.0, 0.0
      write(100,'(I3,3F14.8)') ny, 0.0, grid, 0.0
      write(100,'(I3,3F14.8)') nz, 0.0, 0.0, grid
      write(101,'(I4,3F14.8)') atom_count, minx, miny, minz
      write(101,'(I3,3F14.8)') nx, grid, 0.0, 0.0
      write(101,'(I3,3F14.8)') ny, 0.0, grid, 0.0
      write(101,'(I3,3F14.8)') nz, 0.0, 0.0, grid
      write(102,'(I4,3F14.8)') atom_count, minx, miny, minz
      write(102,'(I3,3F14.8)') nx, grid, 0.0, 0.0
      write(102,'(I3,3F14.8)') ny, 0.0, grid, 0.0
      write(102,'(I3,3F14.8)') nz, 0.0, 0.0, grid
      do ii = 1, atom_count
        write(99,'(I4,4F14.8)') mol(ii)%atom_number, real(mol(ii)%atom_number),&
        mol(ii)%pos(1), mol(ii)%pos(2), mol(ii)%pos(3)
        write(100,'(I4,4F14.8)')mol(ii)%atom_number, real(mol(ii)%atom_number),&
        mol(ii)%pos(1), mol(ii)%pos(2), mol(ii)%pos(3)
        write(101,'(I4,4F14.8)')mol(ii)%atom_number, real(mol(ii)%atom_number),&
        mol(ii)%pos(1), mol(ii)%pos(2), mol(ii)%pos(3)
        write(102,'(I4,4F14.8)')mol(ii)%atom_number, real(mol(ii)%atom_number),&
        mol(ii)%pos(1), mol(ii)%pos(2), mol(ii)%pos(3)
      end do
    else
      write(99,'(I4,3F14.8)') atom_count, minx, miny, minz
      write(99,'(I3,3F14.8)') nx, grid, 0.0, 0.0
      write(99,'(I3,3F14.8)') ny, 0.0, grid, 0.0
      write(99,'(I3,3F14.8)') nz, 0.0, 0.0, grid
      do ii = 1, atom_count
        write(99,'(I4,4F14.8)') mol(ii)%atom_number, real(mol(ii)%atom_number),&
        mol(ii)%pos(1), mol(ii)%pos(2), mol(ii)%pos(3)
      end do
      do ii = 1, len(title)
        if (is_alpha(title(ii:ii))) exit
      end do
    end if
    do ii = 1, cbdm
      inv2sqrtb(ii,:) = 1.0_dp / (2.0_dp*sqrt(cbdata(ii)%expo(:)))
      L = cbdata(ii)%L - 1
      expo = cbdata(ii)%expo
      coe(ii,:) = cbdata(ii)%Ncoe(:)*(-ci)**L/((2.0_dp**(real(L,dp)+1.5_dp))*&
      (expo(:)**(0.5_dp*real(L,dp)+1.5_dp)))
    end do
    do ii = 1, nx
      px = minx + (ii-1) * grid
      do jj = 1, cbdm
        expo   = cbdata(jj)%expo
        fac(:) = cbdata(jj)%fac
        do kk = 1, cbdata(jj)%contr
          Hermitex(ii,jj,kk) = exp(-px*px/(4.0_dp*expo(kk))) * &
                                Hermite_poly(fac(1),px*inv2sqrtb(jj,kk))
        end do
      end do
    end do
    do ii = 1, ny
      py = miny + (ii-1) * grid
      do jj = 1, cbdm
        expo   = cbdata(jj)%expo
        fac(:) = cbdata(jj)%fac
        do kk = 1, cbdata(jj)%contr
          Hermitey(ii,jj,kk) = exp(-py*py/(4.0_dp*expo(kk))) * &
                                Hermite_poly(fac(2),py*inv2sqrtb(jj,kk))
        end do
      end do
    end do
    do ii = 1, nz
      pz = minz + (ii-1) * grid
      do jj = 1, cbdm
        expo   = cbdata(jj)%expo
        fac(:) = cbdata(jj)%fac
        do kk = 1, cbdata(jj)%contr
          Hermitez(ii,jj,kk) = exp(-pz*pz/(4.0_dp*expo(kk))) * &
                                Hermite_poly(fac(3),pz*inv2sqrtb(jj,kk))
        end do
      end do
    end do
    !$omp parallel num_threads(16) default(shared) private(i, ii, jj, kk, ll, &
    !$omp& mm, px, py, pz, rx, ry, rz, val, tldats)
    allocate(tldats(2*nz,ny,nx), source=c0)
    !$omp do schedule(static)
    do i = 1, nx
      ii = i
      px = minx + (ii-1) * grid
      do jj = 1, ny
        py = miny + (jj-1) * grid
        if (is2c) then
          do kk = 1, nz
            pz = minz + (kk-1) * grid
            do ll = 1, cbdm
              val = c0
              rx  = cbdata(ll)%pos(1)
              ry  = cbdata(ll)%pos(2)
              rz  = cbdata(ll)%pos(3)
              do mm = 1, cbdata(ll)%contr
                val = val + coe(ll,mm) * &
                Hermitex(ii,ll,mm)*Hermitey(jj,ll,mm)*Hermitez(kk,ll,mm)
              end do
              val = val * exp(-ci*(px*rx+py*ry+pz*rz))
              tldats(kk,jj,ii) = tldats(kk,jj,ii) + val*arr(ll)
              tldats(kk+nz,jj,ii) = tldats(kk+nz,jj,ii) + val*arr(cbdm+ll)
            end do
          end do
        else
          do kk = 1, nz
            pz = minz + (kk-1) * grid
            do ll = 1, cbdm
              val = c0
              rx  = cbdata(ll)%pos(1)
              ry  = cbdata(ll)%pos(2)
              rz  = cbdata(ll)%pos(3)
              do mm = 1, cbdata(ll)%contr
                val = val + coe(ll,mm) * &
                Hermitex(ii,ll,mm)*Hermitey(jj,ll,mm)*Hermitez(kk,ll,mm)
              end do
              val = val * exp(-ci*(px*rx+py*ry+pz*rz))
              if (spin==1) tldats(kk,jj,ii)=tldats(kk,jj,ii)+val*arr(ll)
              if (spin==2) tldats(kk,jj,ii)=tldats(kk,jj,ii)+val*arr(cbdm+ll)
            end do
          end do
        end if
      end do
    end do
    !$omp end do
    !$omp critical
    if (is2c) then
      datsa(:,:,:) = datsa(:,:,:) + tldats(1:nz,:,:)
      datsb(:,:,:) = datsb(:,:,:) + tldats(nz+1:2*nz,:,:)
    else
      datsa(:,:,:) = datsa(:,:,:) + tldats(1:nz,:,:)
    end if
    !$omp end critical
    deallocate(tldats)
    !$omp end parallel
    do ii = 1, nx
      do jj = 1, ny
        if (is2c) then
          do kk = 1, nz, 6
            end = min(kk+5, 2*nz)
            write(99,'(6ES15.6)') real(datsa(kk:end,jj,ii))
            write(100,'(6ES15.6)') aimag(datsa(kk:end,jj,ii))
            write(101,'(6ES15.6)') real(datsb(kk:end,jj,ii))
            write(102,'(6ES15.6)') aimag(datsb(kk:end,jj,ii))
          end do
        else
          do kk = 1, nz, 6
            end = min(kk+5, nz)
            write(99,'(6ES15.6)') real(datsa(kk:end,jj,ii))
          end do
        end if
      end do
    end do
    close(99)
    if (is2c) then
      close(100)
      close(101)
      close(102)
    end if
    deallocate(datsa, datsb)
    deallocate(Hermitex, Hermitey, Hermitez)
  end subroutine Basis2momentum_cube

!------------------------------------------------------------
!> transform states from basis repre to real space repre (Becke's fuzzy grid)
!!
!! it's used to generate exchange energies, exchange Fock matrix
!!
!! correlation energy, correlation Fock matrix
  subroutine Basis2real_Becke_XandC(rho_m, cAO2MO, ex, ec, Fockx, Fockc)
    implicit none
    real(dp),intent(out)   :: ex         ! exchange energy
    real(dp),intent(out)   :: ec         ! correlation energy
    complex(dp),intent(in) :: rho_m(2*cbdm, 2*cbdm)! density matrix
    complex(dp),intent(in) :: cAO2MO(2*cbdm, 2*fbdm)! coefficient matrix
    integer(8)             :: i,j          ! loop variable
    integer(8)             :: nr, nl     ! number of grid points
    !DIR$ ATTRIBUTES ALIGN:align_size :: pos3w1
    real(dp)               :: pos3w1(590,4,100)
    complex(dp),intent(out):: Fockx(2*cbdm, 2*cbdm)
    complex(dp),intent(out):: Fockc(2*cbdm, 2*cbdm)
    !DIR$ ATTRIBUTES ALIGN:align_size :: Fockxmic, Fockcmic
    real(dp)               :: Fockxmic(2*cbdm, 2*cbdm)
    real(dp)               :: Fockcmic(2*cbdm, 2*cbdm)
    real(dp)               :: exm, ecm
    Fockx = c0
    ex = 0.0_dp
    Fockc = c0
    ec = 0.0_dp
    ! single-centre Gaussian weight integral
    ! Processes all the gird lattice points of an atom center at once,
    ! improving efficiency and avoiding stack overflows

    !$omp parallel num_threads(min(threads,atom_count)) default(shared)&
    !$omp& private(i, nr, nl, Fockxmic, Fockcmic, exm, ecm, pos3w1)&
    !$omp& if(threads < nproc)
    !$omp do schedule(dynamic, 1)
    do i = 1, atom_count
      ! grid precision equal to 'int=ultrafine' in Gaussian program
      if (mol(i)%atom_number <= 2) then
        nr = 50
        nl = 590
      else if (mol(i)%atom_number <= 10) then
        nr = 70
        nl = 590
      else if (mol(i)%atom_number <= 18) then
        nr = 90
        nl = 590
      else
        nr = 100
        nl = 590
      end if
      call Chebyshev2_Lebedev(i, nr, nl, pos3w1)
      select case (xc_f03_func_info_get_family(x_info))
      case (XC_FAMILY_LDA, XC_FAMILY_HYB_LDA)
        call Fock_XandC_LDA(rho_m, nr, nl, pos3w1, exm, ecm, Fockxmic, Fockcmic)
      case (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
        call Fock_XandC_GGA(rho_m, nr, nl, pos3w1, exm, ecm, Fockxmic, Fockcmic)
      case (XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA)
        call Fock_XandC_metaGGA(rho_m, cAO2MO, nr, nl, pos3w1, exm, ecm, &
                                Fockxmic, Fockcmic)
      end select
    end do
    !$omp end do
    !$omp critical
    Fockx = Fockx + Fockxmic
    ex = ex + exm
    Fockc = Fockc + Fockcmic
    ec = ec + ecm
    !$omp end critical
    !$omp end parallel
  end subroutine Basis2real_Becke_XandC

!------------------------------------------------------------
!> transform states from basis repre to real space repre (Becke's fuzzy grid)
!!
!! it's used to generate exchange-correlation energies and Fock matrix
  subroutine Basis2real_Becke_XC(rho_m, cAO2MO, exc, Fockxc)
    implicit none
    real(dp),intent(out)   :: exc        ! exchange-correlation energy
    complex(dp),intent(in) :: rho_m(2*cbdm, 2*cbdm)! density matrix
    complex(dp),intent(in) :: cAO2MO(2*cbdm, 2*fbdm)! coefficient matrix
    integer(8)             :: i          ! loop variable
    integer(8)             :: nr, nl     ! number of grid points
    !DIR$ ATTRIBUTES ALIGN:align_size :: pos3w1
    real(dp)               :: pos3w1(590,4,100)
    complex(dp),intent(out):: Fockxc(2*cbdm, 2*cbdm)
    !DIR$ ATTRIBUTES ALIGN:align_size :: Fockxcmic
    real(dp)               :: Fockxcmic(2*cbdm, 2*cbdm)
    real(dp)               :: excm
    Fockxc = c0
    exc = 0.0_dp
    ! single-centre Gaussian weight integral
    ! Processes all the gird lattice points of an atom center at once,
    ! improving efficiency and avoiding stack overflows

    !$omp parallel num_threads(min(threads,atom_count)) default(shared)&
    !$omp& private(i, nr, nl, Fockxcmic, excm, pos3w1)&
    !$omp& if(threads < nproc)
    !$omp do schedule(dynamic, 1)
    do i = 1, atom_count
      ! integral precision equivalent to 'int=ultrafine' in Gaussian program
      if (mol(i)%atom_number <= 2) then
        nr = 50
        nl = 590
      else if (mol(i)%atom_number <= 10) then
        nr = 70
        nl = 590
      else if (mol(i)%atom_number <= 18) then
        nr = 90
        nl = 590
      else
        nr = 100
        nl = 590
      end if
      call Chebyshev2_Lebedev(i, nr, nl, pos3w1)
      select case (xc_f03_func_info_get_family(xc_info))
      case (XC_FAMILY_LDA, XC_FAMILY_HYB_LDA)
        call Fock_XC_LDA(rho_m, nr, nl, pos3w1, excm, Fockxcmic)
      case (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
        call Fock_XC_GGA(rho_m, nr, nl, pos3w1, excm, Fockxcmic)
      case (XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA)
        call Fock_XC_metaGGA(rho_m, cAO2MO, nr, nl, pos3w1, excm, Fockxcmic)
      end select
    end do
    !$omp end do
    !$omp critical
    Fockxc = Fockxc + Fockxcmic
    exc = exc + excm
    !$omp end critical
    !$omp end parallel
  end subroutine Basis2real_Becke_XC

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
    real(dp),intent(out)  :: pos3w1(590,4,100) ! x, y, z, weight
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
    else if (nl == 590) then
      call LD0590(lx, ly, lz, wl, nl2)
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
    real(dp),intent(out)  :: pos3(590,3,100) ! x, y, z
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
    else if (nl == 590) then
      call LD0590(lx, ly, lz, wl, nl2)
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
    real(dp)                :: vecsum(n)
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
      vecsum(:) = sum(vec(:,:)**2,dim=2)
      contr    = cbdata(kk) % contr
      fac(:)   = cbdata(kk) % fac
      do jj = 1, contr
        expo   = cbdata(kk) % expo(jj)
        coeff  = cbdata(kk) % Ncoe(jj)
        val(:) = val(:) + coeff * exp(-expo*vecsum(:))
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
        fac(:) = cbdata(kk)%fac
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
      fac(:)   = cbdata(kk) % fac
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