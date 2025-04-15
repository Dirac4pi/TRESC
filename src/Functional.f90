!> @file Functional.f90
!!
!! @brief exchange-correlation functionals
!!
!! @syntax Fortran 2008 free format
!!
!! @code UTF-8
!!
!! @author dirac4pi
module functional
  use xc_f03_lib_m
  use Fundamentals
  use Atoms
  use OMP_LIB

  ! exchange and correlation functional infomations
  type(xc_f03_func_info_t) :: x_info
  type(xc_f03_func_info_t) :: c_info
  type(xc_f03_func_t) :: x_func
  type(xc_f03_func_t) :: c_func

  contains

!------------------------------------------------------------
!> grid points, spin=0:total =1:alpha =2:beta
  pure subroutine grid(arr, n, points, spin, dats)
    implicit none
    complex(dp),intent(in) :: arr(2*cbdm)
    integer(8),intent(in) :: n
    real(dp),intent(in) :: points(n, 3)
    integer,intent(in) :: spin
    complex(dp),intent(out) :: dats(n)
    integer :: ii, jj, kk
    integer :: contraction
    real(dp) :: vec(n, 3)
    integer :: L, M
    character(len=4) :: fac
    real(dp) :: coeff
    real(dp) :: expo
    real(dp) :: val(n)
    integer :: init, final, step
    dats = c0
    if (spin == 0) then
      init = 1
      step = 1
      final = 2*cbdm
    else if (spin == 1) then
      init = 1
      step = 2
      final = 2*cbdm-1
    else if (spin == 2) then
      init = 2
      step = 2
      final = 2*cbdm
    end if
    do ii = init, final, step
      kk = int((ii+1)/2)
      val = 0.0_dp
      vec(:,1) = points(:,1) - molecule(basis_inf(kk)%atom)%nucleus_position(1)
      vec(:,2) = points(:,2) - molecule(basis_inf(kk)%atom)%nucleus_position(2)
      vec(:,3) = points(:,3) - molecule(basis_inf(kk)%atom)%nucleus_position(3)
      contraction = atom_basis(molecule(basis_inf(kk)%atom) % &
      basis_number + basis_inf(kk) % shell - 1)%contraction
      L = basis_inf(kk) % L
      M = basis_inf(kk) % M
      fac = AO_xyz_factor(L,M)
      do jj = 1, contraction
        expo = atom_basis(molecule(basis_inf(kk)%atom) % &
        basis_number + basis_inf(kk) % shell - 1)%exponents(jj)
        coeff = atom_basis(molecule(basis_inf(kk)%atom) % &
        basis_number + basis_inf(kk) % shell - 1)%Ncoefficient(jj,M)
        val(:) = val(:) + coeff * exp(-expo*(sum(vec(:,:)**2,dim=2)))
      end do
      do jj = 1, 4
        if (fac(jj:jj) == 'x') then
          val(:) = val(:) * vec(:,1)
        else if(fac(jj:jj) == 'y') then
          val(:) = val(:) * vec(:,2)
        else if(fac(jj:jj) == 'z') then
          val(:) = val(:) * vec(:,3)
        end if
      end do
      dats = dats + val * arr(ii)
    end do
  end subroutine grid

!------------------------------------------------------------
!> norm of gradient of grid points, spin=0:total =1:alpha =2:beta
  pure complex(dp) function gradient(arr, point, spin) result(pointvl)
    implicit none
    complex(dp),intent(in) :: arr(2*cbdm)
    real(dp),intent(in) :: point(3)
    integer,intent(in) :: spin
    integer :: ii, jj, kk
    integer :: contraction
    real(dp) :: vec(3)
    integer :: L, M
    integer :: mx, my, mz
    character(len=4) :: fac
    real(dp) :: coeff
    real(dp) :: expo
    complex(dp) :: val(3), valdx, valdy, valdz
    integer :: init, final, step
    pointvl = c0
    if (spin == 0) then
      init = 1
      step = 1
      final = 2*cbdm
    else if (spin == 1) then
      init = 1
      step = 2
      final = 2*cbdm-1
    else if (spin == 2) then
      init = 2
      step = 2
      final = 2*cbdm
    end if
    val = c0
    do ii = init, final, step
      kk = int((ii+1)/2)
      vec = point - molecule(basis_inf(kk)%atom)%nucleus_position
      contraction = atom_basis(molecule(basis_inf(kk)%atom)%&
      basis_number)%contraction
      L = basis_inf(kk) % L
      M = basis_inf(kk) % M
      fac = AO_xyz_factor(L,M)
      valdx = c0
      valdy = c0
      valdz = c0
      mx = 0
      my = 0
      mz = 0
      do jj = 1, 4
        if (fac(jj:jj) == 'x') then
          mx = mx + 1
        else if(fac(jj:jj) == 'y') then
          my = my + 1
        else if(fac(jj:jj) == 'z') then
          mz = mz + 1
        end if
      end do
      do jj = 1, contraction
        expo = atom_basis(molecule(basis_inf(kk)%atom)%&
        basis_number)%exponents(jj)
        coeff = atom_basis(molecule(basis_inf(kk)%atom)%&
        basis_number)%Ncoefficient(jj,M)
        ! d/dx
        if (mx == 0) then
          valdx = valdx + &
          (-2*expo*vec(1)**(mx+1))*vec(2)**my*vec(3)**mz&
          * coeff * exp(-expo*(vec(1)**2+vec(2)**2+vec(3)**2))
        else
          valdx = valdx + &
          (mx*vec(1)**(mx-1)-2*expo*vec(1)**(mx+1))*vec(2)**my*vec(3)**mz&
          * coeff * exp(-expo*(vec(1)**2+vec(2)**2+vec(3)**2))
        end if
        ! d/dy
        if (my == 0) then
          valdy = valdy + &
          (-2*expo*vec(2)**(my+1))*vec(1)**mx*vec(3)**mz&
          * coeff * exp(-expo*(vec(1)**2+vec(2)**2+vec(3)**2))
        else
          valdy = valdy + &
          (my*vec(2)**(my-1)-2*expo*vec(2)**(my+1))*vec(1)**mx*vec(3)**mz&
          * coeff * exp(-expo*(vec(1)**2+vec(2)**2+vec(3)**2))
        end if
        ! d/dz
        if (mz == 0) then
          valdz = valdz + &
          (-2*expo*vec(3)**(mz+1))*vec(1)**mx*vec(2)**my&
          * coeff * exp(-expo*(vec(1)**2+vec(2)**2+vec(3)**2))
        else
          valdz = valdz + &
          (mz*vec(3)**(mz-1)-2*expo*vec(3)**(mz+1))*vec(1)**mx*vec(2)**my&
          * coeff * exp(-expo*(vec(1)**2+vec(2)**2+vec(3)**2))
        end if
      end do
      val(1) = val(1) + valdx * arr(ii)
      val(2) = val(2) + valdy * arr(ii)
      val(3) = val(3) + valdz * arr(ii)
    end do
    pointvl = sqrt(sum(val**2))
  end function gradient

!------------------------------------------------------------
!> Laplacian of grid points, spin=0:total =1:alpha =2:beta
  pure complex(dp) function Laplacian(arr, point, spin) result(pointvl)
    implicit none
    complex(dp),intent(in) :: arr(2*cbdm)
    real(dp),intent(in) :: point(3)
    integer,intent(in) :: spin
    integer :: ii, jj, kk
    integer :: contraction
    real(dp) :: vec(3)
    integer :: L, M
    integer :: mx, my, mz
    character(len=4) :: fac
    real(dp) :: coeff
    real(dp) :: expo
    complex(dp) :: val, valdx, valdy, valdz
    integer :: init, final, step
    pointvl = c0
    if (spin == 0) then
      init = 1
      step = 1
      final = 2*cbdm
    else if (spin == 1) then
      init = 1
      step = 2
      final = 2*cbdm-1
    else if (spin == 2) then
      init = 2
      step = 2
      final = 2*cbdm
    end if
    do ii = init, final, step
      kk = int((ii+1)/2)
      val = c0
      vec = point - molecule(basis_inf(kk)%atom)%nucleus_position
      contraction = atom_basis(molecule(basis_inf(kk)%atom)%&
      basis_number)%contraction
      L = basis_inf(kk) % L
      M = basis_inf(kk) % M
      fac = AO_xyz_factor(L,M)
      mx = 0
      my = 0
      mz = 0
      do jj = 1, 4
        if (fac(jj:jj) == 'x') then
          mx = mx + 1
        else if(fac(jj:jj) == 'y') then
          my = my + 1
        else if(fac(jj:jj) == 'z') then
          mz = mz + 1
        end if
      end do
      do jj = 1, contraction
        expo = atom_basis(molecule(basis_inf(kk)%atom)%&
        basis_number)%exponents(jj)
        coeff = atom_basis(molecule(basis_inf(kk)%atom)%&
        basis_number)%Ncoefficient(jj,M)
        ! d2/dx2
        if (mx == 0) then
          valdx = (4*expo**2*vec(1)**2 - 2*expo)&
          * vec(2)**my * vec(3)**mz&
          * coeff * exp(-expo*(vec(1)**2+vec(2)**2+vec(3)**2))
        else if (mx == 1) then
          valdx = (4*expo**2*vec(1)**3 - 4*expo*vec(1))&
          * vec(2)**my * vec(3)**mz&
          * coeff * exp(-expo*(vec(1)**2+vec(2)**2+vec(3)**2))
        else
          valdx = ((mx-1)*mx*vec(1)**(mx-2) - 2*expo*(2*mx+1)*vec(1)**mx + &
          4*expo**2*vec(1)**(mx+2))*vec(2)**my*vec(3)**mz&
          * coeff * exp(-expo*(vec(1)**2+vec(2)**2+vec(3)**2))
        end if
        ! d2/dy2
        if (my == 0) then
          valdy = (4*expo**2*vec(2)**2 - 2*expo)&
          * vec(1)**mx * vec(3)**mz&
          * coeff * exp(-expo*(vec(1)**2+vec(2)**2+vec(3)**2))
        else if (my == 1) then
          valdy = (4*expo**2*vec(2)**3 - 4*expo*vec(2))&
          * vec(1)**mx * vec(3)**mz&
          * coeff * exp(-expo*(vec(1)**2+vec(2)**2+vec(3)**2))
        else
          valdy = ((my-1)*my*vec(2)**(my-2) - 2*expo*(2*my+1)*vec(2)**my + &
          4*expo**2*vec(2)**(my+2))*vec(1)**mx*vec(3)**mz&
          * coeff * exp(-expo*(vec(1)**2+vec(2)**2+vec(3)**2))
        end if
        ! d2/dz2
        if (mz == 0) then
          valdz = (4*expo**2*vec(3)**2 - 2*expo)&
          * vec(1)**mx * vec(2)**my&
          * coeff * exp(-expo*(vec(1)**2+vec(2)**2+vec(3)**2))
        else if (mz == 1) then
          valdz = (4*expo**2*vec(3)**3 - 4*expo*vec(3))&
          * vec(1)**mx * vec(2)**my&
          * coeff * exp(-expo*(vec(1)**2+vec(2)**2+vec(3)**2))
        else
          valdz = ((mz-1)*mz*vec(3)**(mz-2) - 2*expo*(2*mz+1)*vec(3)**mz + &
          4*expo**2*vec(3)**(mz+2))*vec(1)**mx*vec(2)**my&
          * coeff * exp(-expo*(vec(1)**2+vec(2)**2+vec(3)**2))
        end if
        val = val + valdx + valdy + valdz ! Laplace operator
      end do
      val = val * arr(ii)
      pointvl = pointvl + val
    end do
  end function Laplacian

!------------------------------------------------------------
!> calculate a series of grid data for xc functional
  pure subroutine grid_xc(&
  mat, n, point, rhoa, rhob, drhoa, drhob, AOAO, dxAOAO, dyAOAO, dzAOAO)
    implicit none
    complex(dp),intent(in) :: mat(2*cbdm, 2*cbdm)
    integer(8),intent(in) :: n
    real(dp),intent(in) :: point(n, 3)
    real(dp),intent(out) :: rhoa(n), rhob(n)
    real(dp),intent(out) :: drhoa(n, 3), drhob(n, 3)
    real(dp),intent(out) :: AOAO(n, cbdm), dxAOAO(n, cbdm)
    real(dp),intent(out) :: dyAOAO(n, cbdm), dzAOAO(n, cbdm)
    complex(dp) :: crhoa(n), crhob(n)
    complex(dp) :: cdrhoa(n, 3), cdrhob(n, 3)
    integer :: ii, jj
    integer :: contraction
    real(dp) :: vec(n, 3)
    integer :: L, M
    integer :: mx, my, mz
    character(len=4) :: fac
    real(dp) :: coeff
    real(dp) :: b
    real(dp) :: val(n), valdx(n), valdy(n), valdz(n)
    real(dp) :: d(n)      ! basis function values
    do ii = 1, cbdm
      val = 0.0_dp
      valdx = 0.0_dp
      valdy = 0.0_dp
      valdz = 0.0_dp
      contraction = atom_basis(molecule(basis_inf(ii)%atom) % &
      basis_number + basis_inf(ii) % shell - 1) % contraction
      L = basis_inf(ii) % L
      M = basis_inf(ii) % M
      fac = AO_xyz_factor(L,M)
      mx = 0
      my = 0
      mz = 0
      do jj = 1, 4
        if (fac(jj:jj) == 'x') then
          mx = mx + 1
        else if(fac(jj:jj) == 'y') then
          my = my + 1
        else if(fac(jj:jj) == 'z') then
          mz = mz + 1
        end if
      end do
      vec(:,1) = point(:,1) - molecule(basis_inf(ii)%atom)%nucleus_position(1)
      vec(:,2) = point(:,2) - molecule(basis_inf(ii)%atom)%nucleus_position(2)
      vec(:,3) = point(:,3) - molecule(basis_inf(ii)%atom)%nucleus_position(3)
      do jj = 1, contraction
        b = atom_basis(molecule(basis_inf(ii)%atom) % &
        basis_number + basis_inf(ii) % shell - 1) % exponents(jj)
        coeff = atom_basis(molecule(basis_inf(ii)%atom) % &
        basis_number + basis_inf(ii) % shell - 1) % Ncoefficient(jj,M)
        d(:) = coeff*exp(-b*(sum(vec(:,:)**2,dim=2)))

        ! calculate the grid value
        ! grid value
        val(:) = val(:) + d(:)

        ! d/dx
        if (mx == 0) then
          valdx(:) = valdx(:) + (-2*b*vec(:,1))*d(:)
        else
          valdx(:) = valdx(:) + (mx*vec(:,1)**(mx-1)-2*b*vec(:,1)**(mx+1))*d(:)
        end if

        ! d/dy
        if (my == 0) then
          valdy(:) = valdy(:) + (-2*b*vec(:,2))*d(:)
        else
          valdy(:) = valdy(:) + (my*vec(:,2)**(my-1)-2*b*vec(:,2)**(my+1))*d(:)
        end if
        
        ! d/dz
        if (mz == 0) then
          valdz(:) = valdz(:) + (-2*b*vec(:,3))*d(:)
        else
          valdz(:) = valdz(:) + (mz*vec(:,3)**(mz-1)-2*b*vec(:,3)**(mz+1))*d(:)
        end if
      end do
      val(:) = val(:) * vec(:,1)**mx*vec(:,2)**my*vec(:,3)**mz
      valdx(:) = valdx(:) * vec(:,2)**my*vec(:,3)**mz
      valdy(:) = valdy(:) * vec(:,1)**mx*vec(:,3)**mz
      valdz(:) = valdz(:) * vec(:,1)**mx*vec(:,2)**my
      AOAO(:,ii) = val(:)
      dxAOAO(:,ii) = valdx(:)
      dyAOAO(:,ii) = valdy(:)
      dzAOAO(:,ii) = valdz(:)
    end do
    crhoa = c0
    crhob = c0
    cdrhoa = c0
    cdrhob = c0
    do ii = 1, cbdm
      do jj = 1, cbdm
        crhoa(:) = crhoa(:) + mat(2*ii-1,2*jj-1)*AOAO(:,ii)*AOAO(:,jj)
        crhob(:) = crhob(:) + mat(2*ii,2*jj)*AOAO(:,ii)*AOAO(:,jj)

        cdrhoa(:,1) = cdrhoa(:,1) + mat(2*ii-1,2*jj-1)*(&
        dxAOAO(:,ii)*AOAO(:,jj) + AOAO(:,ii)*dxAOAO(:,jj))
        cdrhoa(:,2) = cdrhoa(:,2) + mat(2*ii-1,2*jj-1)*(&
        dyAOAO(:,ii)*AOAO(:,jj) + AOAO(:,ii)*dyAOAO(:,jj))
        cdrhoa(:,3) = cdrhoa(:,3) + mat(2*ii-1,2*jj-1)*(&
        dzAOAO(:,ii)*AOAO(:,jj) + AOAO(:,ii)*dzAOAO(:,jj))

        cdrhob(:,1) = cdrhob(:,1) + mat(2*ii,2*jj)*(&
        dxAOAO(:,ii)*AOAO(:,jj) + AOAO(:,ii)*dxAOAO(:,jj))
        cdrhob(:,2) = cdrhob(:,2) + mat(2*ii,2*jj)*(&
        dyAOAO(:,ii)*AOAO(:,jj) + AOAO(:,ii)*dyAOAO(:,jj))
        cdrhob(:,3) = cdrhob(:,3) + mat(2*ii,2*jj)*(&
        dzAOAO(:,ii)*AOAO(:,jj) + AOAO(:,ii)*dzAOAO(:,jj))
      end do
    end do
    rhoa = real(crhoa,dp)
    rhob = real(crhob,dp)
    drhoa = real(cdrhoa,dp)
    drhob = real(cdrhob,dp)
  end subroutine grid_xc

!------------------------------------------------------------
!> assign a group of points value of exchange-correlation Fock
!! matrix element Fockxc(i,j)
!! 
!! variables ref: https://libxc.gitlab.io/manual/previous/libxc-5.0.x/
  subroutine Fockxc_main(rho_m, nr, nl, pos3w1, ex, ec, matx, matc)
    implicit none
    complex(dp),intent(in) :: rho_m(2*cbdm, 2*cbdm)
    integer :: ii, jj, kk, ll
    integer(8),intent(in) :: nr, nl
    real(dp),intent(in) :: pos3w1(100,434,4)
    real(dp),intent(out) :: ex
    real(dp),optional :: ec
    real(dp),intent(out) :: matx(2*cbdm, 2*cbdm)
    real(dp),optional :: matc(2*cbdm, 2*cbdm)
    real(dp) :: pos(nl, 3)
    real(dp) :: weight(nl)
    real(dp) :: rhoa(nl), rhob(nl)
    real(dp) :: xvrho(2*nl), cvrho(2*nl), xvsigma(3*nl), cvsigma(3*nl)
    real(dp) :: rho(2*nl), drhoa(nl, 3), drhob(nl, 3), sigma(3*nl)
    real(dp) :: exin(nl), ecin(nl), exout, ecout
    real(dp) :: matxmic(2*cbdm, 2*cbdm)
    real(dp) :: matcmic(2*cbdm, 2*cbdm)
    real(dp) :: AOAO(nl, cbdm), dxAOAO(nl, cbdm)
    real(dp) :: dyAOAO(nl, cbdm), dzAOAO(nl, cbdm)
    matx = 0.0_dp
    if (present(matc)) matc = 0.0_dp
    ex = 0.0_dp
    if (present(ec)) ec = 0.0_dp
    matxmic = 0.0_dp
    matcmic = 0.0_dp
    exout = 0.0_dp
    ecout = 0.0_dp
    do ii = 1, nr
      pos(:,1) = pos3w1(ii,1:nl,1)
      pos(:,2) = pos3w1(ii,1:nl,2)
      pos(:,3) = pos3w1(ii,1:nl,3)
      weight(:) = pos3w1(ii,1:nl,4)

      call grid_xc(&
      rho_m, nl, pos, rhoa, rhob, drhoa, drhob, AOAO, dxAOAO, dyAOAO, dzAOAO)
      rho(1:2*nl-1:2) = rhoa(:)
      rho(2:2*nl:2) = rhob(:)
      sigma(1:3*nl-2:3) = sum(drhoa(:,:)**2,dim=2)
      sigma(2:3*nl-1:3) = drhoa(:,1)*drhob(:,1) + drhoa(:,2)*drhob(:,2) + &
      drhoa(:,3)*drhob(:,3)
      sigma(3:3*nl:3) = sum(drhob(:,:)**2,dim=2)
      if (fc_id /= -1) then
        select case (xc_f03_func_info_get_family(x_info))
          case(XC_FAMILY_LDA)
            call xc_f03_lda_exc_vxc(x_func, nl, rho, exin, xvrho)
            call xc_f03_lda_exc_vxc(c_func, nl, rho, ecin, cvrho)
            do kk = 1, cbdm
              do ll = kk, cbdm
                ! exchange | alpha
                matxmic(2*kk-1, 2*ll-1) = sum(weight(:)*(&
                xvrho(1:2*nl-1:2)*AOAO(:,kk)*AOAO(:,ll)))
                matxmic(2*ll-1, 2*kk-1) = matxmic(2*kk-1, 2*ll-1)
                ! exchange | beta
                matxmic(2*kk, 2*ll) = sum(weight(:)*(&
                xvrho(2:2*nl:2)*AOAO(:,kk)*AOAO(:,ll)))
                matxmic(2*ll, 2*kk) = matxmic(2*kk, 2*ll)

                ! correlation | alpha
                matcmic(2*kk-1, 2*ll-1) = sum(weight(:)*(&
                cvrho(1:2*nl-1:2)*AOAO(:,kk)*AOAO(:,ll)))
                matcmic(2*ll-1, 2*kk-1) = matcmic(2*kk-1, 2*ll-1)
                ! correlation | beta
                matcmic(2*kk, 2*ll) = sum(weight(:)*(&
                cvrho(2:2*nl:2)*AOAO(:,kk)*AOAO(:,ll)))
                matcmic(2*ll, 2*kk) = matcmic(2*kk, 2*ll)
              end do
            end do
            exout = sum(weight(:)*(rhoa(:)+rhob(:))*exin(:))
            ecout = sum(weight(:)*(rhoa(:)+rhob(:))*ecin(:))
          case(XC_FAMILY_GGA)
            call xc_f03_gga_exc_vxc(x_func, nl, rho, sigma, exin, xvrho,xvsigma)
            call xc_f03_gga_exc_vxc(c_func, nl, rho, sigma, ecin, cvrho,cvsigma)
            do kk = 1, cbdm
              do ll = kk, cbdm
                ! exchange | alpha
                matxmic(2*kk-1, 2*ll-1) = sum(weight(:)*(&
                xvrho(1:2*nl-1:2)*AOAO(:,kk)*AOAO(:,ll) + &
                2.0_dp*xvsigma(1:3*nl-2:3)*(&
                drhoa(:,1)*(dxAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dxAOAO(:,ll)) + &
                drhoa(:,2)*(dyAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dyAOAO(:,ll)) + &
                drhoa(:,3)*(dzAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dzAOAO(:,ll)))+ &
                xvsigma(2:3*nl-1:3)*(&
                drhob(:,1)*(dxAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dxAOAO(:,ll)) + &
                drhob(:,2)*(dyAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dyAOAO(:,ll)) + &
                drhob(:,3)*(dzAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dzAOAO(:,ll)))))
                matxmic(2*ll-1, 2*kk-1) = matxmic(2*kk-1, 2*ll-1)
                ! exchange | beta
                matxmic(2*kk, 2*ll) = sum(weight(:)*(&
                xvrho(2:2*nl:2)*AOAO(:,kk)*AOAO(:,ll) + &
                2.0_dp*xvsigma(3:3*nl:3)*(&
                drhob(:,1)*(dxAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dxAOAO(:,ll)) + &
                drhob(:,2)*(dyAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dyAOAO(:,ll)) + &
                drhob(:,3)*(dzAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dzAOAO(:,ll)))+ &
                xvsigma(2:3*nl-1:3)*(&
                drhoa(:,1)*(dxAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dxAOAO(:,ll)) + &
                drhoa(:,2)*(dyAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dyAOAO(:,ll)) + &
                drhoa(:,3)*(dzAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dzAOAO(:,ll)))))
                matxmic(2*ll, 2*kk) = matxmic(2*kk, 2*ll)

                ! correlation | alpha
                matcmic(2*kk-1, 2*ll-1) = sum(weight(:)*(&
                cvrho(1:2*nl-1:2)*AOAO(:,kk)*AOAO(:,ll) + &
                2.0_dp*cvsigma(1:3*nl-2:3)*(&
                drhoa(:,1)*(dxAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dxAOAO(:,ll)) + &
                drhoa(:,2)*(dyAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dyAOAO(:,ll)) + &
                drhoa(:,3)*(dzAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dzAOAO(:,ll)))+ &
                cvsigma(2:3*nl-1:3)*(&
                drhob(:,1)*(dxAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dxAOAO(:,ll)) + &
                drhob(:,2)*(dyAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dyAOAO(:,ll)) + &
                drhob(:,3)*(dzAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dzAOAO(:,ll)))))
                matcmic(2*ll-1, 2*kk-1) = matcmic(2*kk-1, 2*ll-1)
                ! correlation | beta
                matcmic(2*kk, 2*ll) = sum(weight(:)*(&
                cvrho(2:2*nl:2)*AOAO(:,kk)*AOAO(:,ll) + &
                2.0_dp*cvsigma(3:3*nl:3)*(&
                drhob(:,1)*(dxAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dxAOAO(:,ll)) + &
                drhob(:,2)*(dyAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dyAOAO(:,ll)) + &
                drhob(:,3)*(dzAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dzAOAO(:,ll)))+ &
                cvsigma(2:3*nl-1:3)*(&
                drhoa(:,1)*(dxAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dxAOAO(:,ll)) + &
                drhoa(:,2)*(dyAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dyAOAO(:,ll)) + &
                drhoa(:,3)*(dzAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dzAOAO(:,ll)))))
                matcmic(2*ll, 2*kk) = matcmic(2*kk, 2*ll)
              end do
            end do
            exout = sum(weight(:)*(rhoa(:)+rhob(:))*exin(:))
            ecout = sum(weight(:)*(rhoa(:)+rhob(:))*ecin(:))
        end select
      else
        select case (xc_f03_func_info_get_family(x_info))
          case(XC_FAMILY_LDA)
            call xc_f03_lda_exc_vxc(x_func, nl, rho, exin, xvrho)
            do kk = 1, cbdm
              do ll = kk, cbdm
                ! exchange | alpha
                matxmic(2*kk-1, 2*ll-1) = sum(weight(:)*(&
                xvrho(1:2*nl-1:2)*AOAO(:,kk)*AOAO(:,ll)))
                matxmic(2*ll-1, 2*kk-1) = matxmic(2*kk-1, 2*ll-1)
                ! exchange | beta
                matxmic(2*kk, 2*ll) = sum(weight(:)*(&
                xvrho(2:2*nl:2)*AOAO(:,kk)*AOAO(:,ll)))
                matxmic(2*ll, 2*kk) = matxmic(2*kk, 2*ll)
              end do
            end do
            exout = sum(weight(:)*(rhoa(:)+rhob(:))*exin(:))
          case(XC_FAMILY_GGA)
            call xc_f03_gga_exc_vxc(x_func, nl, rho, sigma, exin, xvrho,xvsigma)
            do kk = 1, cbdm
              do ll = kk, cbdm
                ! exchange | alpha
                matxmic(2*kk-1, 2*ll-1) = sum(weight(:)*(&
                xvrho(1:2*nl-1:2)*AOAO(:,kk)*AOAO(:,ll) + &
                2.0_dp*xvsigma(1:3*nl-2:3)*(&
                drhoa(:,1)*(dxAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dxAOAO(:,ll)) + &
                drhoa(:,2)*(dyAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dyAOAO(:,ll)) + &
                drhoa(:,3)*(dzAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dzAOAO(:,ll)))+ &
                xvsigma(2:3*nl-1:3)*(&
                drhob(:,1)*(dxAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dxAOAO(:,ll)) + &
                drhob(:,2)*(dyAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dyAOAO(:,ll)) + &
                drhob(:,3)*(dzAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dzAOAO(:,ll)))))
                matxmic(2*ll-1, 2*kk-1) = matxmic(2*kk-1, 2*ll-1)
                ! exchange | beta
                matxmic(2*kk, 2*ll) = sum(weight(:)*(&
                xvrho(2:2*nl:2)*AOAO(:,kk)*AOAO(:,ll) + &
                2.0_dp*xvsigma(3:3*nl:3)*(&
                drhob(:,1)*(dxAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dxAOAO(:,ll)) + &
                drhob(:,2)*(dyAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dyAOAO(:,ll)) + &
                drhob(:,3)*(dzAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dzAOAO(:,ll)))+ &
                xvsigma(2:3*nl-1:3)*(&
                drhoa(:,1)*(dxAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dxAOAO(:,ll)) + &
                drhoa(:,2)*(dyAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dyAOAO(:,ll)) + &
                drhoa(:,3)*(dzAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dzAOAO(:,ll)))))
                matxmic(2*ll, 2*kk) = matxmic(2*kk, 2*ll)
              end do
            end do
            exout = sum(weight(:)*(rhoa(:)+rhob(:))*exin(:))
        end select
      end if
      matx = matx + matxmic
      ex = ex + exout
      if (present(matc)) matc = matc + matcmic
      if (present(ec)) ec = ec + ecout
    end do
  end subroutine Fockxc_main

!------------------------------------------------------------
!> initialising Libxc, if x is a xc functional,
!! you can leave out c_id
  subroutine Fockxc_init()
    implicit none
    integer :: ii, vmajor, vminor, vmicro
    character :: x_name, c_name
    ! Print out the version  
    call xc_f03_version(vmajor, vminor, vmicro)
    write(60,'("  -- Libxc version: ",I1,".",I1,".",I1)') vmajor, vminor, vmicro
    call xc_f03_func_init(x_func, fx_id, XC_POLARIZED)
    x_info = xc_f03_func_get_info(x_func)
    if (xc_f03_func_info_get_kind(x_info) == XC_EXCHANGE) then
      if (fc_id /= -1) then
        call xc_f03_func_init(c_func, fc_id, XC_POLARIZED)
        c_info = xc_f03_func_get_info(c_func)
      else
        call terminate('Fockxc_init: missing data for correlation functional')
      end if
    end if
    write(60,'(A13,I5,A12,A)') &
    "  -- fx_id = ",fx_id,';func name: ',trim(xc_f03_func_info_get_name(x_info))
    if (fc_id /= -1) then
    write(60,'(A13,I5,A12,A)') &
    "  -- fc_id = ",fc_id,';func name: ',trim(xc_f03_func_info_get_name(c_info))
    end if
  end subroutine Fockxc_init

!------------------------------------------------------------
!> deallocate memory of Libxc, if x is a xc functional,
!! you can leave out c_func
  subroutine Fockxc_end()
    implicit none
    call xc_f03_func_end(x_func)
    if (fc_id /= -1) call xc_f03_func_end(c_func)
  end subroutine Fockxc_end

end module functional