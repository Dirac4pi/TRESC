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
  type(xc_f03_func_info_t) :: xc_info
  type(xc_f03_func_t) :: x_func
  type(xc_f03_func_t) :: c_func
  type(xc_f03_func_t) :: xc_func

  contains

!------------------------------------------------------------
!> calculate a series of density grid data for LDA functionals
  pure subroutine Grid_XC_LDA(mat, n, point, rhoa, rhob, AOAO)
    implicit none
    complex(dp),intent(in) :: mat(2*cbdm, 2*cbdm)
    integer(8),intent(in)  :: n
    real(dp),intent(in)    :: point(n, 3)
    real(dp),intent(out)   :: rhoa(n), rhob(n)
    real(dp),intent(out)   :: AOAO(n, cbdm)
    complex(dp)            :: crhoa(n), crhob(n)
    integer                :: ii, jj
    integer                :: contr
    real(dp)               :: vec(n, 3)
    integer                :: L, M
    integer                :: fac(3)
    real(dp)               :: coeff
    real(dp)               :: b
    real(dp)               :: val(n)
    real(dp)               :: d(n)      ! Gaussian core function value
    do ii = 1, cbdm
      val = 0.0_dp
      contr = cbdata(ii) % contr
      L     = cbdata(ii) % L
      M     = cbdata(ii) % M
      fac = AO_fac(:,L,M)
      vec(:,1) = point(:,1) - cbdata(ii) % pos(1)
      vec(:,2) = point(:,2) - cbdata(ii) % pos(2)
      vec(:,3) = point(:,3) - cbdata(ii) % pos(3)
      do jj = 1, contr
        b     = cbdata(ii) % expo(jj)
        coeff = cbdata(ii) % Ncoe(jj)
        d(:)  = coeff*exp(-b*(sum(vec(:,:)**2,dim=2)))
        ! calculate the Grid value
        ! Grid value
        val(:) = val(:) + d(:)
      end do
      AOAO(:,ii) = val(:) * vec(:,1)**fac(1)*vec(:,2)**fac(2)*vec(:,3)**fac(3)
    end do
    crhoa = c0
    crhob = c0
    do ii = 1, cbdm
      do jj = 1, cbdm
        crhoa(:) = crhoa(:) + mat(ii,jj)*AOAO(:,ii)*AOAO(:,jj)
        crhob(:) = crhob(:) + mat(cbdm+ii,cbdm+jj)*AOAO(:,ii)*AOAO(:,jj)
      end do
    end do
    rhoa = real(crhoa,dp)
    rhob = real(crhob,dp)
  end subroutine Grid_XC_LDA


!------------------------------------------------------------
!> calculate a series of density grid data for GGA functionals
  pure subroutine Grid_XC_GGA(&
  mat, n, point, rhoa, rhob, drhoa, drhob, AOAO, dxAOAO, dyAOAO, dzAOAO)
    implicit none
    complex(dp),intent(in)   :: mat(2*cbdm, 2*cbdm)
    integer(8),intent(in)    :: n
    real(dp),intent(in)      :: point(n, 3)
    real(dp),intent(out)     :: rhoa(n), rhob(n)
    real(dp),intent(out)     :: drhoa(n, 3), drhob(n, 3)
    real(dp),intent(out)     :: AOAO(n, cbdm), dxAOAO(n, cbdm)
    real(dp),intent(out)     :: dyAOAO(n, cbdm), dzAOAO(n, cbdm)
    complex(dp)              :: crhoa(n), crhob(n)
    complex(dp)              :: cdrhoa(n, 3), cdrhob(n, 3)
    integer                  :: ii, jj
    integer                  :: contr
    real(dp)                 :: vec(n, 3)
    integer                  :: L, M
    integer                  :: fac(3)
    real(dp)                 :: coeff
    real(dp)                 :: b
    real(dp)                 :: val(n), valdx(n), valdy(n), valdz(n)
    real(dp)                 :: d(n)      ! Gaussian core function value
    do ii = 1, cbdm
      val = 0.0_dp
      valdx = 0.0_dp
      valdy = 0.0_dp
      valdz = 0.0_dp
      contr = cbdata(ii) % contr
      L     = cbdata(ii) % L
      M     = cbdata(ii) % M
      fac = AO_fac(:,L,M)
      vec(:,1) = point(:,1) - cbdata(ii) % pos(1)
      vec(:,2) = point(:,2) - cbdata(ii) % pos(2)
      vec(:,3) = point(:,3) - cbdata(ii) % pos(3)
      do jj = 1, contr
        b     = cbdata(ii) % expo(jj)
        coeff = cbdata(ii) % Ncoe(jj)
        d(:)  = coeff*exp(-b*(sum(vec(:,:)**2,dim=2)))

        ! calculate the Grid value
        ! Grid value
        val(:) = val(:) + d(:)

        ! d/dx
        if (fac(1) == 0) then
          valdx(:) = valdx(:) + (-2*b*vec(:,1))*d(:)
        else
          valdx(:) = valdx(:) + &
          (fac(1)*vec(:,1)**(fac(1)-1)-2*b*vec(:,1)**(fac(1)+1))*d(:)
        end if

        ! d/dy
        if (fac(2) == 0) then
          valdy(:) = valdy(:) + (-2*b*vec(:,2))*d(:)
        else
          valdy(:) = valdy(:) + &
          (fac(2)*vec(:,2)**(fac(2)-1)-2*b*vec(:,2)**(fac(2)+1))*d(:)
        end if
        
        ! d/dz
        if (fac(3) == 0) then
          valdz(:) = valdz(:) + (-2*b*vec(:,3))*d(:)
        else
          valdz(:) = valdz(:) + &
          (fac(3)*vec(:,3)**(fac(3)-1)-2*b*vec(:,3)**(fac(3)+1))*d(:)
        end if
      end do
      AOAO(:,ii) = val(:) * vec(:,1)**fac(1)*vec(:,2)**fac(2)*vec(:,3)**fac(3)
      dxAOAO(:,ii) = valdx(:) * vec(:,2)**fac(2)*vec(:,3)**fac(3)
      dyAOAO(:,ii) = valdy(:) * vec(:,1)**fac(1)*vec(:,3)**fac(3)
      dzAOAO(:,ii) = valdz(:) * vec(:,1)**fac(1)*vec(:,2)**fac(2)
    end do
    crhoa = c0
    crhob = c0
    cdrhoa = c0
    cdrhob = c0
    do ii = 1, cbdm
      do jj = 1, cbdm
        crhoa(:) = crhoa(:) + mat(ii,jj)*AOAO(:,ii)*AOAO(:,jj)
        crhob(:) = crhob(:) + mat(cbdm+ii,cbdm+jj)*AOAO(:,ii)*AOAO(:,jj)

        cdrhoa(:,1) = cdrhoa(:,1) + mat(ii,jj)*&
                                    (dxAOAO(:,ii)*AOAO(:,jj) + &
                                    AOAO(:,ii)*dxAOAO(:,jj))
        cdrhoa(:,2) = cdrhoa(:,2) + mat(ii,jj)*&
                                    (dyAOAO(:,ii)*AOAO(:,jj) + &
                                    AOAO(:,ii)*dyAOAO(:,jj))
        cdrhoa(:,3) = cdrhoa(:,3) + mat(ii,jj)*&
                                    (dzAOAO(:,ii)*AOAO(:,jj) + &
                                    AOAO(:,ii)*dzAOAO(:,jj))

        cdrhob(:,1) = cdrhob(:,1) + mat(cbdm+ii,cbdm+jj)*&
                                    (dxAOAO(:,ii)*AOAO(:,jj) + &
                                    AOAO(:,ii)*dxAOAO(:,jj))
        cdrhob(:,2) = cdrhob(:,2) + mat(cbdm+ii,cbdm+jj)*&
                                    (dyAOAO(:,ii)*AOAO(:,jj) + &
                                    AOAO(:,ii)*dyAOAO(:,jj))
        cdrhob(:,3) = cdrhob(:,3) + mat(cbdm+ii,cbdm+jj)*&
                                    (dzAOAO(:,ii)*AOAO(:,jj) + &
                                    AOAO(:,ii)*dzAOAO(:,jj))
      end do
    end do
    rhoa = real(crhoa,dp)
    rhob = real(crhob,dp)
    drhoa = real(cdrhoa,dp)
    drhob = real(cdrhob,dp)
  end subroutine Grid_XC_GGA

!------------------------------------------------------------
!> calculate a series of density grid data for meta-GGA functionals
  pure subroutine Grid_XC_metaGGA(&
  mat, n, point, rhoa, rhob, drhoa, drhob, AOAO, dxAOAO, dyAOAO, dzAOAO, &
  dx2AOAO, dy2AOAO, dz2AOAO, lapla, laplb)
    implicit none
    complex(dp),intent(in)   :: mat(2*cbdm, 2*cbdm)
    integer(8),intent(in)    :: n
    real(dp),intent(in)      :: point(n, 3)
    real(dp),intent(out)     :: rhoa(n), rhob(n)
    real(dp),intent(out)     :: drhoa(n, 3), drhob(n, 3)
    real(dp),intent(out)     :: AOAO(n, cbdm)
    real(dp),intent(out)     :: dxAOAO(n,cbdm), dyAOAO(n,cbdm), dzAOAO(n,cbdm)
    real(dp),intent(out)     :: dx2AOAO(n,cbdm),dy2AOAO(n,cbdm),dz2AOAO(n,cbdm)
    real(dp),intent(out)     :: lapla(n), laplb(n)
    complex(dp)              :: crhoa(n), crhob(n)
    complex(dp)              :: cdrhoa(n, 3), cdrhob(n, 3)
    complex(dp)              :: clapla(n, 3), claplb(n, 3)
    integer                  :: ii, jj
    integer                  :: contr
    real(dp)                 :: vec(n, 3)
    integer                  :: L, M
    integer                  :: fac(3)
    real(dp)                 :: coeff
    real(dp)                 :: b
    real(dp)                 :: val(n), valdx(n), valdy(n), valdz(n)
    real(dp)                 :: vlaplx(n), vlaply(n), vlaplz(n)
    real(dp)                 :: d(n)      ! Gaussian core function value
    do ii = 1, cbdm
      val = 0.0_dp
      valdx = 0.0_dp
      valdy = 0.0_dp
      valdz = 0.0_dp
      vlaplx = 0.0_dp
      vlaply = 0.0_dp
      vlaplz = 0.0_dp
      contr = cbdata(ii) % contr
      L     = cbdata(ii) % L
      M     = cbdata(ii) % M
      fac = AO_fac(:,L,M)
      vec(:,1) = point(:,1) - cbdata(ii) % pos(1)
      vec(:,2) = point(:,2) - cbdata(ii) % pos(2)
      vec(:,3) = point(:,3) - cbdata(ii) % pos(3)
      do jj = 1, contr
        b     = cbdata(ii) % expo(jj)
        coeff = cbdata(ii) % Ncoe(jj)
        d(:)  = coeff*exp(-b*(sum(vec(:,:)**2,dim=2)))

        ! calculate the Grid value
        ! Grid value
        val(:) = val(:) + d(:)

        ! d/dx
        if (fac(1) == 0) then
          valdx(:) = valdx(:) + (-2*b*vec(:,1))*d(:)
        else
          valdx(:) = valdx(:) + &
          (fac(1)*vec(:,1)**(fac(1)-1)-2*b*vec(:,1)**(fac(1)+1))*d(:)
        end if

        ! d/dy
        if (fac(2) == 0) then
          valdy(:) = valdy(:) + (-2*b*vec(:,2))*d(:)
        else
          valdy(:) = valdy(:) + &
          (fac(2)*vec(:,2)**(fac(2)-1)-2*b*vec(:,2)**(fac(2)+1))*d(:)
        end if
        
        ! d/dz
        if (fac(3) == 0) then
          valdz(:) = valdz(:) + (-2*b*vec(:,3))*d(:)
        else
          valdz(:) = valdz(:) + &
          (fac(3)*vec(:,3)**(fac(3)-1)-2*b*vec(:,3)**(fac(3)+1))*d(:)
        end if

        ! d2/dx2
        if (fac(1) == 0) then
          vlaplx(:) = vlaplx(:) + 2*b*(2*b*vec(:,1)**2 - 1)*d(:)
        else if (fac(1) == 1) then
          vlaplx(:) = vlaplx(:) + 2*b*(2*b*vec(:,1)**3 - 3*vec(:,1))*d(:)
        else
          vlaplx(:) = vlaplx(:) + &
          (fac(1)*(fac(1)-1)*vec(:,1)**(fac(1)-2) - &
          2*b*(2*fac(1)+1)*vec(:,1)**fac(1) + &
          4*b*b*vec(:,1)**(fac(1)+2))*d(:)
        end if

        ! d2/dy2
        if (fac(2) == 0) then
          vlaply(:) = vlaply(:) + 2*b*(2*b*vec(:,2)**2 - 1)*d(:)
        else if (fac(2) == 1) then
          vlaply(:) = vlaply(:) + 2*b*(2*b*vec(:,2)**3 - 3*vec(:,2))*d(:)
        else
          vlaply(:) = vlaply(:) + &
          (fac(2)*(fac(2)-1)*vec(:,2)**(fac(2)-2) - &
          2*b*(2*fac(2)+1)*vec(:,2)**fac(2) + &
          4*b*b*vec(:,2)**(fac(2)+2))*d(:)
        end if

        ! d2/dz2
        if (fac(3) == 0) then
          vlaplz(:) = vlaplz(:) + 2*b*(2*b*vec(:,3)**2 - 1)*d(:)
        else if (fac(3) == 1) then
          vlaplz(:) = vlaplz(:) + 2*b*(2*b*vec(:,3)**3 - 3*vec(:,3))*d(:)
        else
          vlaplz(:) = vlaplz(:) + &
          (fac(3)*(fac(3)-1)*vec(:,3)**(fac(3)-2) - &
          2*b*(2*fac(3)+1)*vec(:,3)**fac(3) + &
          4*b*b*vec(:,3)**(fac(3)+2))*d(:)
        end if
      end do
      AOAO(:,ii) = val(:) * vec(:,1)**fac(1)*vec(:,2)**fac(2)*vec(:,3)**fac(3)
      dxAOAO(:,ii) = valdx(:) * vec(:,2)**fac(2)*vec(:,3)**fac(3)
      dyAOAO(:,ii) = valdy(:) * vec(:,1)**fac(1)*vec(:,3)**fac(3)
      dzAOAO(:,ii) = valdz(:) * vec(:,1)**fac(1)*vec(:,2)**fac(2)
      dx2AOAO(:,ii) = vlaplx(:) * vec(:,2)**fac(2)*vec(:,3)**fac(3)
      dy2AOAO(:,ii) = vlaply(:) * vec(:,1)**fac(1)*vec(:,3)**fac(3)
      dz2AOAO(:,ii) = vlaplz(:) * vec(:,1)**fac(1)*vec(:,2)**fac(2)
    end do
    crhoa = c0
    crhob = c0
    cdrhoa = c0
    cdrhob = c0
    clapla = c0
    claplb = c0
    do ii = 1, cbdm
      do jj = 1, cbdm
        crhoa(:) = crhoa(:) + mat(ii,jj)*AOAO(:,ii)*AOAO(:,jj)
        crhob(:) = crhob(:) + mat(cbdm+ii,cbdm+jj)*AOAO(:,ii)*AOAO(:,jj)

        cdrhoa(:,1) = cdrhoa(:,1) + mat(ii,jj)*&
                                    (dxAOAO(:,ii)*AOAO(:,jj) + &
                                    AOAO(:,ii)*dxAOAO(:,jj))
        cdrhoa(:,2) = cdrhoa(:,2) + mat(ii,jj)*&
                                    (dyAOAO(:,ii)*AOAO(:,jj) + &
                                    AOAO(:,ii)*dyAOAO(:,jj))
        cdrhoa(:,3) = cdrhoa(:,3) + mat(ii,jj)*&
                                    (dzAOAO(:,ii)*AOAO(:,jj) + &
                                    AOAO(:,ii)*dzAOAO(:,jj))

        cdrhob(:,1) = cdrhob(:,1) + mat(cbdm+ii,cbdm+jj)*&
                                    (dxAOAO(:,ii)*AOAO(:,jj) + &
                                    AOAO(:,ii)*dxAOAO(:,jj))
        cdrhob(:,2) = cdrhob(:,2) + mat(cbdm+ii,cbdm+jj)*&
                                    (dyAOAO(:,ii)*AOAO(:,jj) + &
                                    AOAO(:,ii)*dyAOAO(:,jj))
        cdrhob(:,3) = cdrhob(:,3) + mat(cbdm+ii,cbdm+jj)*&
                                    (dzAOAO(:,ii)*AOAO(:,jj) + &
                                    AOAO(:,ii)*dzAOAO(:,jj))

        clapla(:,1) = clapla(:,1) + mat(ii,jj)*&
                                    (dx2AOAO(:,ii)*AOAO(:,jj) + &
                                    2*dxAOAO(:,ii)*dxAOAO(:,jj) + &
                                    AOAO(:,ii)*dx2AOAO(:,jj))
        clapla(:,2) = clapla(:,2) + mat(ii,jj)*&
                                    (dy2AOAO(:,ii)*AOAO(:,jj) + &
                                    2*dyAOAO(:,ii)*dyAOAO(:,jj) + &
                                    AOAO(:,ii)*dy2AOAO(:,jj))
        clapla(:,3) = clapla(:,3) + mat(ii,jj)*&
                                    (dz2AOAO(:,ii)*AOAO(:,jj) + &
                                    2*dzAOAO(:,ii)*dzAOAO(:,jj) + &
                                    AOAO(:,ii)*dz2AOAO(:,jj))

        claplb(:,1) = claplb(:,1) + mat(cbdm+ii,cbdm+jj)*&
                                    (dx2AOAO(:,ii)*AOAO(:,jj) + &
                                    2*dxAOAO(:,ii)*dxAOAO(:,jj) + &
                                    AOAO(:,ii)*dx2AOAO(:,jj))
        claplb(:,2) = claplb(:,2) + mat(cbdm+ii,cbdm+jj)*&
                                    (dy2AOAO(:,ii)*AOAO(:,jj) + &
                                    2*dyAOAO(:,ii)*dyAOAO(:,jj) + &
                                    AOAO(:,ii)*dy2AOAO(:,jj))
        claplb(:,3) = claplb(:,3) + mat(cbdm+ii,cbdm+jj)*&
                                    (dz2AOAO(:,ii)*AOAO(:,jj) + &
                                    2*dzAOAO(:,ii)*dzAOAO(:,jj) + &
                                    AOAO(:,ii)*dz2AOAO(:,jj))
      end do
    end do
    rhoa = real(crhoa,dp)
    rhob = real(crhob,dp)
    drhoa = real(cdrhoa,dp)
    drhob = real(cdrhob,dp)
    lapla = real(clapla(:,1)+clapla(:,2)+clapla(:,3),dp)
    laplb = real(claplb(:,1)+claplb(:,2)+claplb(:,3),dp)
  end subroutine Grid_XC_metaGGA

!------------------------------------------------------------
!> calculate a series of non-relativistic Lagrangian kinetic energy density grid
!! data for meta-GGA functionals
  pure subroutine Lag_Ek_density(mat, Nele, n, point, taua, taub)
  implicit none
    complex(dp),intent(in)   :: mat(2*cbdm, 2*fbdm)
    integer,intent(in)       :: Nele
    integer(8),intent(in)    :: n
    real(dp),intent(in)      :: point(n, 3)
    real(dp),intent(out)     :: taua(n), taub(n)
    real(dp)                 :: p2a(n), p2b(n)
    real(dp)                 :: dxAOAO(n,cbdm), dyAOAO(n,cbdm), dzAOAO(n,cbdm)
    integer                  :: ii, jj
    integer                  :: contr
    real(dp)                 :: vec(n, 3)
    integer                  :: L, M
    integer                  :: fac(3)
    real(dp)                 :: coeff
    real(dp)                 :: b
    real(dp)                 :: valdx(n), valdy(n), valdz(n)
    real(dp)                 :: d(n)      ! Gaussian core function value
    complex(dp)              :: Pxpsi(n), Pypsi(n), Pzpsi(n)
    do ii = 1, cbdm
      valdx = 0.0_dp
      valdy = 0.0_dp
      valdz = 0.0_dp
      contr = cbdata(ii) % contr
      L     = cbdata(ii) % L
      M     = cbdata(ii) % M
      fac   = AO_fac(:,L,M)
      vec(:,1) = point(:,1) - cbdata(ii) % pos(1)
      vec(:,2) = point(:,2) - cbdata(ii) % pos(2)
      vec(:,3) = point(:,3) - cbdata(ii) % pos(3)
      do jj = 1, contr
        b     = cbdata(ii) % expo(jj)
        coeff = cbdata(ii) % Ncoe(jj)
        d(:)  = coeff*exp(-b*(sum(vec(:,:)**2,dim=2)))

        ! d/dx
        if (fac(1) == 0) then
          valdx(:) = valdx(:) + (-2*b*vec(:,1))*d(:)
        else
          valdx(:) = valdx(:) + &
          (fac(1)*vec(:,1)**(fac(1)-1)-2*b*vec(:,1)**(fac(1)+1))*d(:)
        end if

        ! d/dy
        if (fac(2) == 0) then
          valdy(:) = valdy(:) + (-2*b*vec(:,2))*d(:)
        else
          valdy(:) = valdy(:) + &
          (fac(2)*vec(:,2)**(fac(2)-1)-2*b*vec(:,2)**(fac(2)+1))*d(:)
        end if

        ! d/dz
        if (fac(3) == 0) then
          valdz(:) = valdz(:) + (-2*b*vec(:,3))*d(:)
        else
          valdz(:) = valdz(:) + &
          (fac(3)*vec(:,3)**(fac(3)-1)-2*b*vec(:,3)**(fac(3)+1))*d(:)
        end if
      end do
      dxAOAO(:,ii) = valdx(:) * vec(:,2)**fac(2)*vec(:,3)**fac(3)
      dyAOAO(:,ii) = valdy(:) * vec(:,1)**fac(1)*vec(:,3)**fac(3)
      dzAOAO(:,ii) = valdz(:) * vec(:,1)**fac(1)*vec(:,2)**fac(2)
    end do
    p2a = 0.0_dp
    p2b = 0.0_dp
    do ii = 1, Nele
      Pxpsi = c0
      Pypsi = c0
      Pzpsi = c0
      do jj = 1, cbdm
        Pxpsi(:) = Pxpsi(:) + mat(jj,ii) * dxAOAO(:,jj)
        Pypsi(:) = Pypsi(:) + mat(jj,ii) * dyAOAO(:,jj)
        Pzpsi(:) = Pzpsi(:) + mat(jj,ii) * dzAOAO(:,jj)
      end do
      p2a(:) = p2a(:) + real(Pxpsi(:)*conjg(Pxpsi(:)) + &
                             Pypsi(:)*conjg(Pypsi(:)) + &
                             Pzpsi(:)*conjg(Pzpsi(:)),dp)
      Pxpsi = c0
      Pypsi = c0
      Pzpsi = c0
      do jj = 1, cbdm
        Pxpsi(:) = Pxpsi(:) + mat(cbdm+jj,ii) * dxAOAO(:,jj)
        Pypsi(:) = Pypsi(:) + mat(cbdm+jj,ii) * dyAOAO(:,jj)
        Pzpsi(:) = Pzpsi(:) + mat(cbdm+jj,ii) * dzAOAO(:,jj)
      end do
      p2b(:) = p2b(:) + real(Pxpsi(:)*conjg(Pxpsi(:)) + &
                             Pypsi(:)*conjg(Pypsi(:)) + &
                             Pzpsi(:)*conjg(Pzpsi(:)),dp)
    end do
    do ii = 1, n
      taua(ii) = dsqrt(p2a(ii)+c2)*c - c2
      taub(ii) = dsqrt(p2b(ii)+c2)*c - c2
    end do
  end subroutine Lag_Ek_density
!------------------------------------------------------------
!> assign a group of points value of LDA exchange Fock matrix Fockx(i_j)
!! and correlation Fock matrix Fockc(i_j)
!! 
!! variables ref: https://libxc.gitlab.io/manual/previous/libxc-5.0.x/
  subroutine Fock_XandC_LDA(rho_m, nr, nl, pos3w1, ex, ec, matx, matc)
    implicit none
    complex(dp),intent(in) :: rho_m(2*cbdm, 2*cbdm)
    integer                :: ii, jj, kk, ll
    integer(8),intent(in)  :: nr, nl
    real(dp),intent(in)    :: pos3w1(590,4,100)
    real(dp),intent(out)   :: ex
    real(dp),intent(out)   :: ec
    real(dp),intent(out)   :: matx(2*cbdm, 2*cbdm)
    real(dp),intent(out)   :: matc(2*cbdm, 2*cbdm)
    !DIR$ ATTRIBUTES ALIGN:align_size :: pos, weight, rhoa, rhob
    real(dp)               :: pos(nl, 3)
    real(dp)               :: weight(nl)
    real(dp)               :: rhoa(nl), rhob(nl)
    !DIR$ ATTRIBUTES ALIGN:align_size :: xvrho, cvrho
    real(dp)               :: xvrho(2*nl), cvrho(2*nl)
    !DIR$ ATTRIBUTES ALIGN:align_size :: rho
    real(dp)               :: rho(2*nl)
    !DIR$ ATTRIBUTES ALIGN:align_size :: exin, ecin
    real(dp)               :: exin(nl), ecin(nl)
    !DIR$ ATTRIBUTES ALIGN:align_size :: matxmic, matcmic
    real(dp)               :: matxmic(2*cbdm, 2*cbdm)
    real(dp)               :: matcmic(2*cbdm, 2*cbdm)
    !DIR$ ATTRIBUTES ALIGN:align_size :: AOAO
    real(dp)               :: AOAO(nl, cbdm)
    matx    = 0.0_dp
    matc    = 0.0_dp
    ex      = 0.0_dp
    ec      = 0.0_dp
    do ii = 1, nr
      pos(:,1) = pos3w1(1:nl,1,ii)
      pos(:,2) = pos3w1(1:nl,2,ii)
      pos(:,3) = pos3w1(1:nl,3,ii)
      weight(:) = pos3w1(1:nl,4,ii)
      call Grid_XC_LDA(rho_m, nl, pos, rhoa, rhob, AOAO)
      rho(1:2*nl-1:2) = rhoa(:)
      rho(2:2*nl:2) = rhob(:)
      call xc_f03_lda_exc_vxc(x_func, nl, rho, exin, xvrho)
      call xc_f03_lda_exc_vxc(c_func, nl, rho, ecin, cvrho)
      matxmic = 0.0_dp
      matcmic = 0.0_dp
      do kk = 1, cbdm
        do ll = kk, cbdm
          ! exchange | alpha
          matxmic(kk, ll) = sum(weight(:)*(&
          xvrho(1:2*nl-1:2)*AOAO(:,kk)*AOAO(:,ll)))
          matxmic(ll, kk) = matxmic(kk, ll)
          ! exchange | beta
          matxmic(cbdm+kk, cbdm+ll) = sum(weight(:)*(&
          xvrho(2:2*nl:2)*AOAO(:,kk)*AOAO(:,ll)))
          matxmic(cbdm+ll, cbdm+kk) = matxmic(cbdm+kk, cbdm+ll)
          ! correlation | alpha
          matcmic(kk, ll) = sum(weight(:)*(&
          cvrho(1:2*nl-1:2)*AOAO(:,kk)*AOAO(:,ll)))
          matcmic(ll, kk) = matcmic(kk, ll)
          ! correlation | beta
          matcmic(cbdm+kk, cbdm+ll) = sum(weight(:)*(&
          cvrho(2:2*nl:2)*AOAO(:,kk)*AOAO(:,ll)))
          matcmic(cbdm+ll, cbdm+kk) = matcmic(cbdm+kk, cbdm+ll)
        end do
      end do
      matx  = matx + matxmic
      ex    = ex + sum(weight(:)*(rhoa(:)+rhob(:))*exin(:))
      matc  = matc + matcmic
      ec    = ec + sum(weight(:)*(rhoa(:)+rhob(:))*ecin(:))
    end do
  end subroutine Fock_XandC_LDA

!------------------------------------------------------------
!> assign a group of points value of GGA exchange Fock matrix Fockx(i_j)
!! and correlation Fock matrix Fockc(i_j)
!! 
!! variables ref: https://libxc.gitlab.io/manual/previous/libxc-5.0.x/
  subroutine Fock_XandC_GGA(rho_m, nr, nl, pos3w1, ex, ec, matx, matc)
    implicit none
    complex(dp),intent(in) :: rho_m(2*cbdm, 2*cbdm)
    integer                :: ii, jj, kk, ll
    integer(8),intent(in)  :: nr, nl
    real(dp),intent(in)    :: pos3w1(590,4,100)
    real(dp),intent(out)   :: ex
    real(dp),intent(out)   :: ec
    real(dp),intent(out)   :: matx(2*cbdm, 2*cbdm)
    real(dp),intent(out)   :: matc(2*cbdm, 2*cbdm)
    !DIR$ ATTRIBUTES ALIGN:align_size :: pos, weight, rhoa, rhob
    real(dp)               :: pos(nl, 3)
    real(dp)         :: weight(nl)
    real(dp)         :: rhoa(nl), rhob(nl)
    !DIR$ ATTRIBUTES ALIGN:align_size :: xvrho, cvrho, xvsigma, cvsigma
    real(dp)         :: xvrho(2*nl), cvrho(2*nl), xvsigma(3*nl), cvsigma(3*nl)
    !DIR$ ATTRIBUTES ALIGN:align_size :: rho, drhoa, drhob, sigma
    real(dp)         :: rho(2*nl), drhoa(nl, 3), drhob(nl, 3), sigma(3*nl)
    !DIR$ ATTRIBUTES ALIGN:align_size :: exin, ecin
    real(dp)         :: exin(nl), ecin(nl)
    !DIR$ ATTRIBUTES ALIGN:align_size :: matxmic, matcmic
    real(dp)         :: matxmic(2*cbdm, 2*cbdm)
    real(dp)         :: matcmic(2*cbdm, 2*cbdm)
    !DIR$ ATTRIBUTES ALIGN:align_size :: AOAO, dxAOAO, dyAOAO, dzAOAO
    real(dp)         :: AOAO(nl, cbdm), dxAOAO(nl, cbdm)
    real(dp)         :: dyAOAO(nl, cbdm), dzAOAO(nl, cbdm)
    matx    = 0.0_dp
    matc    = 0.0_dp
    ex      = 0.0_dp
    ec      = 0.0_dp
    do ii = 1, nr
      pos(:,1) = pos3w1(1:nl,1,ii)
      pos(:,2) = pos3w1(1:nl,2,ii)
      pos(:,3) = pos3w1(1:nl,3,ii)
      weight(:) = pos3w1(1:nl,4,ii)
      call Grid_XC_GGA(&
      rho_m, nl, pos, rhoa, rhob, drhoa, drhob, AOAO, dxAOAO, dyAOAO, dzAOAO)
      rho(1:2*nl-1:2) = rhoa(:)
      rho(2:2*nl:2) = rhob(:)
      sigma(1:3*nl-2:3) = sum(drhoa(:,:)**2,dim=2)
      sigma(2:3*nl-1:3) = drhoa(:,1)*drhob(:,1) + drhoa(:,2)*drhob(:,2) + &
      drhoa(:,3)*drhob(:,3)
      sigma(3:3*nl:3) = sum(drhob(:,:)**2,dim=2)
      call xc_f03_gga_exc_vxc(x_func, nl, rho, sigma, exin, xvrho,xvsigma)
      call xc_f03_gga_exc_vxc(c_func, nl, rho, sigma, ecin, cvrho,cvsigma)
      matxmic = 0.0_dp
      matcmic = 0.0_dp
      do kk = 1, cbdm
        do ll = kk, cbdm
          ! exchange | alpha
          matxmic(kk, ll) = sum(weight(:)*(&
          xvrho(1:2*nl-1:2)*AOAO(:,kk)*AOAO(:,ll) + &
          2.0_dp*xvsigma(1:3*nl-2:3)*(&
          drhoa(:,1)*(dxAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dxAOAO(:,ll)) + &
          drhoa(:,2)*(dyAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dyAOAO(:,ll)) + &
          drhoa(:,3)*(dzAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dzAOAO(:,ll)))+ &
          xvsigma(2:3*nl-1:3)*(&
          drhob(:,1)*(dxAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dxAOAO(:,ll)) + &
          drhob(:,2)*(dyAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dyAOAO(:,ll)) + &
          drhob(:,3)*(dzAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dzAOAO(:,ll)))))
          matxmic(ll, kk) = matxmic(kk, ll)
          ! exchange | beta
          matxmic(cbdm+kk, cbdm+ll) = sum(weight(:)*(&
          xvrho(2:2*nl:2)*AOAO(:,kk)*AOAO(:,ll) + &
          2.0_dp*xvsigma(3:3*nl:3)*(&
          drhob(:,1)*(dxAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dxAOAO(:,ll)) + &
          drhob(:,2)*(dyAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dyAOAO(:,ll)) + &
          drhob(:,3)*(dzAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dzAOAO(:,ll)))+ &
          xvsigma(2:3*nl-1:3)*(&
          drhoa(:,1)*(dxAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dxAOAO(:,ll)) + &
          drhoa(:,2)*(dyAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dyAOAO(:,ll)) + &
          drhoa(:,3)*(dzAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dzAOAO(:,ll)))))
          matxmic(cbdm+ll, cbdm+kk) = matxmic(cbdm+kk, cbdm+ll)
          ! correlation | alpha
          matcmic(kk, ll) = sum(weight(:)*(&
          cvrho(1:2*nl-1:2)*AOAO(:,kk)*AOAO(:,ll) + &
          2.0_dp*cvsigma(1:3*nl-2:3)*(&
          drhoa(:,1)*(dxAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dxAOAO(:,ll)) + &
          drhoa(:,2)*(dyAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dyAOAO(:,ll)) + &
          drhoa(:,3)*(dzAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dzAOAO(:,ll)))+ &
          cvsigma(2:3*nl-1:3)*(&
          drhob(:,1)*(dxAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dxAOAO(:,ll)) + &
          drhob(:,2)*(dyAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dyAOAO(:,ll)) + &
          drhob(:,3)*(dzAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dzAOAO(:,ll)))))
          matcmic(ll, kk) = matcmic(kk, ll)
          ! correlation | beta
          matcmic(cbdm+kk, cbdm+ll) = sum(weight(:)*(&
          cvrho(2:2*nl:2)*AOAO(:,kk)*AOAO(:,ll) + &
          2.0_dp*cvsigma(3:3*nl:3)*(&
          drhob(:,1)*(dxAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dxAOAO(:,ll)) + &
          drhob(:,2)*(dyAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dyAOAO(:,ll)) + &
          drhob(:,3)*(dzAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dzAOAO(:,ll)))+ &
          cvsigma(2:3*nl-1:3)*(&
          drhoa(:,1)*(dxAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dxAOAO(:,ll)) + &
          drhoa(:,2)*(dyAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dyAOAO(:,ll)) + &
          drhoa(:,3)*(dzAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dzAOAO(:,ll)))))
          matcmic(cbdm+ll, cbdm+kk) = matcmic(cbdm+kk, cbdm+ll)
        end do
      end do
      matx  = matx + matxmic
      ex    = ex + sum(weight(:)*(rhoa(:)+rhob(:))*exin(:))
      matc  = matc + matcmic
      ec    = ec + sum(weight(:)*(rhoa(:)+rhob(:))*ecin(:))
    end do
  end subroutine Fock_XandC_GGA

!------------------------------------------------------------
!> assign a group of points value of meta-GGA exchange Fock matrix Fockx(i_j)
!! and correlation Fock matrix Fockc(i_j)
!! 
!! variables ref: https://libxc.gitlab.io/manual/previous/libxc-5.0.x/
  subroutine Fock_XandC_metaGGA(rho_m,cAO2MO,nr,nl,pos3w1,ex,ec,matx,matc)
    implicit none
    complex(dp),intent(in) :: rho_m(2*cbdm, 2*cbdm)
    complex(dp),intent(in) :: cAO2MO(2*cbdm, 2*fbdm)
    integer                :: ii, jj, kk, ll
    integer(8),intent(in)  :: nr, nl
    real(dp),intent(in)    :: pos3w1(590,4,100)
    real(dp),intent(out)   :: ex
    real(dp),intent(out)   :: ec
    real(dp),intent(out)   :: matx(2*cbdm, 2*cbdm)
    real(dp),intent(out)   :: matc(2*cbdm, 2*cbdm)
    !DIR$ ATTRIBUTES ALIGN:align_size :: pos, weight, rhoa, rhob
    real(dp)               :: pos(nl, 3)
    real(dp)         :: weight(nl)
    real(dp)         :: rhoa(nl), rhob(nl)
    !DIR$ ATTRIBUTES ALIGN:align_size :: xvrho, cvrho, xvsigma, cvsigma
    real(dp)         :: xvrho(2*nl), cvrho(2*nl), xvsigma(3*nl), cvsigma(3*nl)
    !DIR$ ATTRIBUTES ALIGN:align_size :: xvlapl, xvtau, cvlapl, cvtau
    real(dp)         :: xvlapl(2*nl), xvtau(2*nl), cvlapl(2*nl), cvtau(2*nl)
    !DIR$ ATTRIBUTES ALIGN:align_size :: rho, drhoa, drhob, sigma
    real(dp)         :: rho(2*nl), drhoa(nl, 3), drhob(nl, 3), sigma(3*nl)
    !DIR$ ATTRIBUTES ALIGN:align_size :: exin, ecin
    real(dp)         :: exin(nl), ecin(nl)
    !DIR$ ATTRIBUTES ALIGN:align_size :: matxmic, matcmic
    real(dp)         :: matxmic(2*cbdm, 2*cbdm)
    real(dp)         :: matcmic(2*cbdm, 2*cbdm)
    !DIR$ ATTRIBUTES ALIGN:align_size :: AOAO, dxAOAO, dyAOAO, dzAOAO
    real(dp)         :: AOAO(nl, cbdm)
    real(dp)         :: dxAOAO(nl, cbdm), dyAOAO(nl, cbdm), dzAOAO(nl, cbdm)
    !DIR$ ATTRIBUTES ALIGN:align_size :: dx2AOAO, dy2AOAO, dz2AOAO
    real(dp)         :: dx2AOAO(nl,cbdm),dy2AOAO(nl,cbdm),dz2AOAO(nl,cbdm)
    !DIR$ ATTRIBUTES ALIGN:align_size :: lapla, laplb, lapl
    real(dp)         :: lapla(nl), laplb(nl), lapl(2*nl)
    !DIR$ ATTRIBUTES ALIGN:align_size :: taua, taub, tau
    real(dp)         :: taua(nl), taub(nl), tau(2*nl)
    matx    = 0.0_dp
    matc    = 0.0_dp
    ex      = 0.0_dp
    ec      = 0.0_dp
    do ii = 1, nr
      pos(:,1) = pos3w1(1:nl,1,ii)
      pos(:,2) = pos3w1(1:nl,2,ii)
      pos(:,3) = pos3w1(1:nl,3,ii)
      weight(:) = pos3w1(1:nl,4,ii)
      call Grid_XC_metaGGA(&
      rho_m, nl, pos, rhoa, rhob, drhoa, drhob, AOAO, dxAOAO, &
      dyAOAO, dzAOAO, dx2AOAO, dy2AOAO, dz2AOAO, lapla, laplb)
      call Lag_Ek_density(cAO2MO, electron_count, nl, pos, taua, taub)
      rho(1:2*nl-1:2) = rhoa(:)
      rho(2:2*nl:2) = rhob(:)
      lapl(1:2*nl-1:2) = lapla(:)
      lapl(2:2*nl:2) = laplb(:)
      tau(1:2*nl-1:2) = taua(:)
      tau(2:2*nl:2) = taub(:)
      sigma(1:3*nl-2:3) = sum(drhoa(:,:)**2,dim=2)
      sigma(2:3*nl-1:3) = drhoa(:,1)*drhob(:,1) + drhoa(:,2)*drhob(:,2) + &
      drhoa(:,3)*drhob(:,3)
      sigma(3:3*nl:3) = sum(drhob(:,:)**2,dim=2)
      call xc_f03_mgga_exc_vxc(x_func, nl, rho, sigma, lapl, tau, exin, &
                               xvrho, xvsigma, xvlapl, xvtau)
      call xc_f03_mgga_exc_vxc(c_func, nl, rho, sigma, lapl, tau, ecin, &
                               cvrho, cvsigma, cvlapl, cvtau)
      matxmic = 0.0_dp
      matcmic = 0.0_dp
      do kk = 1, cbdm
        do ll = kk, cbdm
          ! exchange | alpha
          matxmic(kk, ll) = sum(weight(:)*(&
          ! GGA part: density, gradient-----------------------------------
          xvrho(1:2*nl-1:2)*AOAO(:,kk)*AOAO(:,ll)                      + &
          2.0_dp*xvsigma(1:3*nl-2:3)*(&
          drhoa(:,1)*(dxAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dxAOAO(:,ll)) + &
          drhoa(:,2)*(dyAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dyAOAO(:,ll)) + &
          drhoa(:,3)*(dzAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dzAOAO(:,ll)))+ &
          xvsigma(2:3*nl-1:3)*(&
          drhob(:,1)*(dxAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dxAOAO(:,ll)) + &
          drhob(:,2)*(dyAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dyAOAO(:,ll)) + &
          drhob(:,3)*(dzAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dzAOAO(:,ll)))+ &
          ! meta-GGA part: Laplacian, kinetic energy----------------------
          xvlapl(1:2*nl-1:2)*(&
          dx2AOAO(:,kk)*AOAO(:,ll)                                     + &
          2.0_dp*dxAOAO(:,kk)*dxAOAO(:,ll)                             + &
          AOAO(:,kk)*dx2AOAO(:,ll)                                     + &
          dy2AOAO(:,kk)*AOAO(:,ll)                                     + &
          2.0_dp*dyAOAO(:,kk)*dyAOAO(:,ll)                             + &
          AOAO(:,kk)*dy2AOAO(:,ll)                                     + &
          dz2AOAO(:,kk)*AOAO(:,ll)                                     + &
          2.0_dp*dzAOAO(:,kk)*dzAOAO(:,ll)                             + &
          AOAO(:,kk)*dz2AOAO(:,ll))                                    + &
          0.5_dp*xvtau(1:2*nl-1:2)*(&
          dxAOAO(:,kk)*dxAOAO(:,ll)                                    + &
          dyAOAO(:,kk)*dyAOAO(:,ll)                                    + &
          dzAOAO(:,kk)*dzAOAO(:,ll))))
          matxmic(ll, kk) = matxmic(kk, ll)
          ! exchange | beta
          matxmic(cbdm+kk, cbdm+ll) = sum(weight(:)*(&
          ! GGA part: density, gradient-----------------------------------
          xvrho(2:2*nl:2)*AOAO(:,kk)*AOAO(:,ll) + &
          2.0_dp*xvsigma(3:3*nl:3)*(&
          drhob(:,1)*(dxAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dxAOAO(:,ll)) + &
          drhob(:,2)*(dyAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dyAOAO(:,ll)) + &
          drhob(:,3)*(dzAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dzAOAO(:,ll)))+ &
          xvsigma(2:3*nl-1:3)*(&
          drhoa(:,1)*(dxAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dxAOAO(:,ll)) + &
          drhoa(:,2)*(dyAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dyAOAO(:,ll)) + &
          drhoa(:,3)*(dzAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dzAOAO(:,ll)))+ &
          ! meta-GGA part: Laplacian, kinetic energy----------------------
          xvlapl(2:2*nl:2)*(&
          dx2AOAO(:,kk)*AOAO(:,ll)                                     + &
          2.0_dp*dxAOAO(:,kk)*dxAOAO(:,ll)                             + &
          AOAO(:,kk)*dx2AOAO(:,ll)                                     + &
          dy2AOAO(:,kk)*AOAO(:,ll)                                     + &
          2.0_dp*dyAOAO(:,kk)*dyAOAO(:,ll)                             + &
          AOAO(:,kk)*dy2AOAO(:,ll)                                     + &
          dz2AOAO(:,kk)*AOAO(:,ll)                                     + &
          2.0_dp*dzAOAO(:,kk)*dzAOAO(:,ll)                             + &
          AOAO(:,kk)*dz2AOAO(:,ll))                                    + &
          0.5_dp*xvtau(2:2*nl:2)*(&
          dxAOAO(:,kk)*dxAOAO(:,ll)                                    + &
          dyAOAO(:,kk)*dyAOAO(:,ll)                                    + &
          dzAOAO(:,kk)*dzAOAO(:,ll))))
          matxmic(cbdm+ll, cbdm+kk) = matxmic(cbdm+kk, cbdm+ll)
          ! correlation | alpha
          matcmic(kk, ll) = sum(weight(:)*(&
          ! GGA part: density, gradient-----------------------------------
          cvrho(1:2*nl-1:2)*AOAO(:,kk)*AOAO(:,ll)                     + &
          2.0_dp*cvsigma(1:3*nl-2:3)*(&
          drhoa(:,1)*(dxAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dxAOAO(:,ll)) + &
          drhoa(:,2)*(dyAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dyAOAO(:,ll)) + &
          drhoa(:,3)*(dzAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dzAOAO(:,ll)))+ &
          cvsigma(2:3*nl-1:3)*(&
          drhob(:,1)*(dxAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dxAOAO(:,ll)) + &
          drhob(:,2)*(dyAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dyAOAO(:,ll)) + &
          drhob(:,3)*(dzAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dzAOAO(:,ll)))+ &
          ! meta-GGA part: Laplacian, kinetic energy----------------------
          cvlapl(1:2*nl-1:2)*(&
          dx2AOAO(:,kk)*AOAO(:,ll)                                     + &
          2.0_dp*dxAOAO(:,kk)*dxAOAO(:,ll)                             + &
          AOAO(:,kk)*dx2AOAO(:,ll)                                     + &
          dy2AOAO(:,kk)*AOAO(:,ll)                                     + &
          2.0_dp*dyAOAO(:,kk)*dyAOAO(:,ll)                             + &
          AOAO(:,kk)*dy2AOAO(:,ll)                                     + &
          dz2AOAO(:,kk)*AOAO(:,ll)                                     + &
          2.0_dp*dzAOAO(:,kk)*dzAOAO(:,ll)                             + &
          AOAO(:,kk)*dz2AOAO(:,ll))                                    + &
          0.5_dp*cvtau(1:2*nl-1:2)*(&
          dxAOAO(:,kk)*dxAOAO(:,ll)                                    + &
          dyAOAO(:,kk)*dyAOAO(:,ll)                                    + &
          dzAOAO(:,kk)*dzAOAO(:,ll))))
          matcmic(ll, kk) = matcmic(kk, ll)
          ! correlation | beta
          matcmic(cbdm+kk, cbdm+ll) = sum(weight(:)*(&
          ! GGA part: density, gradient-----------------------------------
          cvrho(2:2*nl:2)*AOAO(:,kk)*AOAO(:,ll) + &
          2.0_dp*cvsigma(3:3*nl:3)*(&
          drhob(:,1)*(dxAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dxAOAO(:,ll)) + &
          drhob(:,2)*(dyAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dyAOAO(:,ll)) + &
          drhob(:,3)*(dzAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dzAOAO(:,ll)))+ &
          cvsigma(2:3*nl-1:3)*(&
          drhoa(:,1)*(dxAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dxAOAO(:,ll)) + &
          drhoa(:,2)*(dyAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dyAOAO(:,ll)) + &
          drhoa(:,3)*(dzAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dzAOAO(:,ll)))+ &
          ! meta-GGA part: Laplacian, kinetic energy----------------------
          cvlapl(2:2*nl:2)*(&
          dx2AOAO(:,kk)*AOAO(:,ll)                                     + &
          2.0_dp*dxAOAO(:,kk)*dxAOAO(:,ll)                             + &
          AOAO(:,kk)*dx2AOAO(:,ll)                                     + &
          dy2AOAO(:,kk)*AOAO(:,ll)                                     + &
          2.0_dp*dyAOAO(:,kk)*dyAOAO(:,ll)                             + &
          AOAO(:,kk)*dy2AOAO(:,ll)                                     + &
          dz2AOAO(:,kk)*AOAO(:,ll)                                     + &
          2.0_dp*dzAOAO(:,kk)*dzAOAO(:,ll)                             + &
          AOAO(:,kk)*dz2AOAO(:,ll))                                    + &
          0.5_dp*cvtau(2:2*nl:2)*(&
          dxAOAO(:,kk)*dxAOAO(:,ll)                                    + &
          dyAOAO(:,kk)*dyAOAO(:,ll)                                    + &
          dzAOAO(:,kk)*dzAOAO(:,ll))))
          matcmic(cbdm+ll, cbdm+kk) = matcmic(cbdm+kk, cbdm+ll)
        end do
      end do
      matx  = matx + matxmic
      ex    = ex + sum(weight(:)*(rhoa(:)+rhob(:))*exin(:))
      matc  = matc + matcmic
      ec    = ec + sum(weight(:)*(rhoa(:)+rhob(:))*ecin(:))
    end do
  end subroutine Fock_XandC_metaGGA

!------------------------------------------------------------
!> assign a group of points value of LDA exchange-correlation Fock
!! matrix element Fockxc(i,j)
!! 
!! variables ref: https://libxc.gitlab.io/manual/previous/libxc-5.0.x/
  subroutine Fock_XC_LDA(rho_m, nr, nl, pos3w1, exc, matxc)
    implicit none
    complex(dp),intent(in) :: rho_m(2*cbdm, 2*cbdm)
    integer                :: ii, jj, kk, ll
    integer(8),intent(in)  :: nr, nl
    real(dp),intent(in)    :: pos3w1(590,4,100)
    real(dp),intent(out)   :: exc
    real(dp),intent(out)   :: matxc(2*cbdm, 2*cbdm)
    !DIR$ ATTRIBUTES ALIGN:align_size :: pos, weight, rhoa, rhob
    real(dp)               :: pos(nl, 3)
    real(dp)               :: weight(nl)
    real(dp)               :: rhoa(nl), rhob(nl)
    !DIR$ ATTRIBUTES ALIGN:align_size :: xcvrho
    real(dp)               :: xcvrho(2*nl)
    !DIR$ ATTRIBUTES ALIGN:align_size :: rho
    real(dp)               :: rho(2*nl)
    !DIR$ ATTRIBUTES ALIGN:align_size :: excin
    real(dp)               :: excin(nl)
    !DIR$ ATTRIBUTES ALIGN:align_size :: matxcmic
    real(dp)               :: matxcmic(2*cbdm, 2*cbdm)
    !DIR$ ATTRIBUTES ALIGN:align_size :: AOAO
    real(dp)               :: AOAO(nl, cbdm)
    matxc    = 0.0_dp
    exc      = 0.0_dp
    do ii = 1, nr
      pos(:,1) = pos3w1(1:nl,1,ii)
      pos(:,2) = pos3w1(1:nl,2,ii)
      pos(:,3) = pos3w1(1:nl,3,ii)
      weight(:) = pos3w1(1:nl,4,ii)
      call Grid_XC_LDA(rho_m, nl, pos, rhoa, rhob, AOAO)
      rho(1:2*nl-1:2) = rhoa(:)
      rho(2:2*nl:2) = rhob(:)
      call xc_f03_lda_exc_vxc(xc_func, nl, rho, excin, xcvrho)
      matxcmic = 0.0_dp
      do kk = 1, cbdm
        do ll = kk, cbdm
          ! alpha
          matxcmic(kk, ll) = sum(weight(:)*(&
          xcvrho(1:2*nl-1:2)*AOAO(:,kk)*AOAO(:,ll)))
          matxcmic(ll, kk) = matxcmic(kk, ll)
          ! beta
          matxcmic(cbdm+kk, cbdm+ll) = sum(weight(:)*(&
          xcvrho(2:2*nl:2)*AOAO(:,kk)*AOAO(:,ll)))
          matxcmic(cbdm+ll, cbdm+kk) = matxcmic(cbdm+kk, cbdm+ll)
        end do
      end do
      matxc  = matxc + matxcmic
      exc    = exc + sum(weight(:)*(rhoa(:)+rhob(:))*excin(:))
    end do
  end subroutine Fock_XC_LDA

!------------------------------------------------------------
!> assign a group of points value of GGA exchange-correlation Fock
!! matrix element Fockxc(i,j)
!! 
!! variables ref: https://libxc.gitlab.io/manual/previous/libxc-5.0.x/
  subroutine Fock_XC_GGA(rho_m, nr, nl, pos3w1, exc, matxc)
    implicit none
    complex(dp),intent(in) :: rho_m(2*cbdm, 2*cbdm)
    integer                :: ii, jj, kk, ll
    integer(8),intent(in)  :: nr, nl
    real(dp),intent(in)    :: pos3w1(590,4,100)
    real(dp),intent(out)   :: exc
    real(dp),intent(out)   :: matxc(2*cbdm, 2*cbdm)
    !DIR$ ATTRIBUTES ALIGN:align_size :: pos, weight, rhoa, rhob
    real(dp)               :: pos(nl, 3)
    real(dp)               :: weight(nl)
    real(dp)               :: rhoa(nl), rhob(nl)
    !DIR$ ATTRIBUTES ALIGN:align_size :: xcvrho, xcvsigma
    real(dp)               :: xcvrho(2*nl), xcvsigma(3*nl)
    !DIR$ ATTRIBUTES ALIGN:align_size :: rho, drhoa, drhob, sigma
    real(dp)               :: rho(2*nl), drhoa(nl, 3), drhob(nl, 3), sigma(3*nl)
    !DIR$ ATTRIBUTES ALIGN:align_size :: excin
    real(dp)               :: excin(nl)
    !DIR$ ATTRIBUTES ALIGN:align_size :: matxcmic
    real(dp)               :: matxcmic(2*cbdm, 2*cbdm)
    !DIR$ ATTRIBUTES ALIGN:align_size :: AOAO, dxAOAO, dyAOAO, dzAOAO
    real(dp)               :: AOAO(nl, cbdm), dxAOAO(nl, cbdm)
    real(dp)               :: dyAOAO(nl, cbdm), dzAOAO(nl, cbdm)
    matxc    = 0.0_dp
    exc      = 0.0_dp
    do ii = 1, nr
      pos(:,1) = pos3w1(1:nl,1,ii)
      pos(:,2) = pos3w1(1:nl,2,ii)
      pos(:,3) = pos3w1(1:nl,3,ii)
      weight(:) = pos3w1(1:nl,4,ii)
      call Grid_XC_GGA(&
      rho_m, nl, pos, rhoa, rhob, drhoa, drhob, AOAO, dxAOAO, dyAOAO, dzAOAO)
      rho(1:2*nl-1:2) = rhoa(:)
      rho(2:2*nl:2) = rhob(:)
      sigma(1:3*nl-2:3) = sum(drhoa(:,:)**2,dim=2)
      sigma(2:3*nl-1:3) = drhoa(:,1)*drhob(:,1) + drhoa(:,2)*drhob(:,2) + &
      drhoa(:,3)*drhob(:,3)
      sigma(3:3*nl:3) = sum(drhob(:,:)**2,dim=2)
      call xc_f03_gga_exc_vxc(xc_func,nl,rho,sigma,excin,xcvrho,xcvsigma)
      matxcmic = 0.0_dp
      do kk = 1, cbdm
        do ll = kk, cbdm
          ! alpha
          matxcmic(kk, ll) = sum(weight(:)*(&
          xcvrho(1:2*nl-1:2)*AOAO(:,kk)*AOAO(:,ll) + &
          2.0_dp*xcvsigma(1:3*nl-2:3)*(&
          drhoa(:,1)*(dxAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dxAOAO(:,ll)) + &
          drhoa(:,2)*(dyAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dyAOAO(:,ll)) + &
          drhoa(:,3)*(dzAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dzAOAO(:,ll)))+ &
          xcvsigma(2:3*nl-1:3)*(&
          drhob(:,1)*(dxAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dxAOAO(:,ll)) + &
          drhob(:,2)*(dyAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dyAOAO(:,ll)) + &
          drhob(:,3)*(dzAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dzAOAO(:,ll)))))
          matxcmic(ll, kk) = matxcmic(kk, ll)
          ! beta
          matxcmic(cbdm+kk, cbdm+ll) = sum(weight(:)*(&
          xcvrho(2:2*nl:2)*AOAO(:,kk)*AOAO(:,ll) + &
          2.0_dp*xcvsigma(3:3*nl:3)*(&
          drhob(:,1)*(dxAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dxAOAO(:,ll)) + &
          drhob(:,2)*(dyAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dyAOAO(:,ll)) + &
          drhob(:,3)*(dzAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dzAOAO(:,ll)))+ &
          xcvsigma(2:3*nl-1:3)*(&
          drhoa(:,1)*(dxAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dxAOAO(:,ll)) + &
          drhoa(:,2)*(dyAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dyAOAO(:,ll)) + &
          drhoa(:,3)*(dzAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dzAOAO(:,ll)))))
          matxcmic(cbdm+ll, cbdm+kk) = matxcmic(cbdm+kk, cbdm+ll)
        end do
      end do
      matxc  = matxc + matxcmic
      exc    = exc + sum(weight(:)*(rhoa(:)+rhob(:))*excin(:))
    end do
  end subroutine Fock_XC_GGA

!------------------------------------------------------------
!> assign a group of points value of meta-GGA exchange-correlation Fock
!! matrix element Fockxc(i,j)
!! 
!! variables ref: https://libxc.gitlab.io/manual/previous/libxc-5.0.x/
  subroutine Fock_XC_metaGGA(rho_m, cAO2MO, nr, nl, pos3w1, exc, matxc)
    implicit none
    complex(dp),intent(in) :: rho_m(2*cbdm, 2*cbdm)
    complex(dp),intent(in) :: cAO2MO(2*cbdm, 2*fbdm)
    integer                :: ii, jj, kk, ll
    integer(8),intent(in)  :: nr, nl
    real(dp),intent(in)    :: pos3w1(590,4,100)
    real(dp),intent(out)   :: exc
    real(dp),intent(out)   :: matxc(2*cbdm, 2*cbdm)
    !DIR$ ATTRIBUTES ALIGN:align_size :: pos, weight, rhoa, rhob
    real(dp)               :: pos(nl, 3)
    real(dp)               :: weight(nl)
    real(dp)               :: rhoa(nl), rhob(nl)
    !DIR$ ATTRIBUTES ALIGN:align_size :: xcvrho, xcvsigma
    real(dp)               :: xcvrho(2*nl), xcvsigma(3*nl)
    !DIR$ ATTRIBUTES ALIGN:align_size :: xcvlapl, xcvtau
    real(dp)               :: xcvlapl(2*nl), xcvtau(2*nl)
    !DIR$ ATTRIBUTES ALIGN:align_size :: rho, drhoa, drhob, sigma
    real(dp)               :: rho(2*nl), drhoa(nl, 3), drhob(nl, 3), sigma(3*nl)
    !DIR$ ATTRIBUTES ALIGN:align_size :: excin
    real(dp)               :: excin(nl)
    !DIR$ ATTRIBUTES ALIGN:align_size :: matxcmic
    real(dp)               :: matxcmic(2*cbdm, 2*cbdm)
    !DIR$ ATTRIBUTES ALIGN:align_size :: AOAO, dxAOAO, dyAOAO, dzAOAO
    real(dp)               :: AOAO(nl, cbdm)
    real(dp)               :: dxAOAO(nl,cbdm), dyAOAO(nl,cbdm), dzAOAO(nl,cbdm)
    !DIR$ ATTRIBUTES ALIGN:align_size :: dx2AOAO, dy2AOAO, dz2AOAO
    real(dp)               :: dx2AOAO(nl,cbdm),dy2AOAO(nl,cbdm),dz2AOAO(nl,cbdm)
    !DIR$ ATTRIBUTES ALIGN:align_size :: lapla, laplb, lapl
    real(dp)               :: lapla(nl), laplb(nl), lapl(2*nl)
    !DIR$ ATTRIBUTES ALIGN:align_size :: taua, taub, tau
    real(dp)               :: taua(nl), taub(nl), tau(2*nl)
    matxc    = 0.0_dp
    exc      = 0.0_dp
    do ii = 1, nr
      pos(:,1) = pos3w1(1:nl,1,ii)
      pos(:,2) = pos3w1(1:nl,2,ii)
      pos(:,3) = pos3w1(1:nl,3,ii)
      weight(:) = pos3w1(1:nl,4,ii)
      call Grid_XC_metaGGA(&
      rho_m, nl, pos, rhoa, rhob, drhoa, drhob, AOAO, dxAOAO, &
      dyAOAO, dzAOAO, dx2AOAO, dy2AOAO, dz2AOAO, lapla, laplb)
      call Lag_Ek_density(cAO2MO, electron_count, nl, pos, taua, taub)
      rho(1:2*nl-1:2) = rhoa(:)
      rho(2:2*nl:2) = rhob(:)
      lapl(1:2*nl-1:2) = lapla(:)
      lapl(2:2*nl:2) = laplb(:)
      tau(1:2*nl-1:2) = taua(:)
      tau(2:2*nl:2) = taub(:)
      sigma(1:3*nl-2:3) = sum(drhoa(:,:)**2,dim=2)
      sigma(2:3*nl-1:3) = drhoa(:,1)*drhob(:,1) + drhoa(:,2)*drhob(:,2) + &
      drhoa(:,3)*drhob(:,3)
      sigma(3:3*nl:3) = sum(drhob(:,:)**2,dim=2)
      call xc_f03_mgga_exc_vxc(xc_func, nl, rho, sigma, lapl, tau, excin, &
                               xcvrho, xcvsigma, xcvlapl, xcvtau)
      matxcmic = 0.0_dp
      do kk = 1, cbdm
        do ll = kk, cbdm
          ! alpha
          matxcmic(kk, ll) = sum(weight(:)*(&
          ! GGA part: density, gradient-----------------------------------
          xcvrho(1:2*nl-1:2)*AOAO(:,kk)*AOAO(:,ll)                     + &
          2.0_dp*xcvsigma(1:3*nl-2:3)*(&
          drhoa(:,1)*(dxAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dxAOAO(:,ll)) + &
          drhoa(:,2)*(dyAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dyAOAO(:,ll)) + &
          drhoa(:,3)*(dzAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dzAOAO(:,ll)))+ &
          xcvsigma(2:3*nl-1:3)*(&
          drhob(:,1)*(dxAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dxAOAO(:,ll)) + &
          drhob(:,2)*(dyAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dyAOAO(:,ll)) + &
          drhob(:,3)*(dzAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dzAOAO(:,ll)))+ &
          ! meta-GGA part: Laplacian, kinetic energy----------------------
          xcvlapl(1:2*nl-1:2)*(&
          dx2AOAO(:,kk)*AOAO(:,ll)                                     + &
          2.0_dp*dxAOAO(:,kk)*dxAOAO(:,ll)                             + &
          AOAO(:,kk)*dx2AOAO(:,ll)                                     + &
          dy2AOAO(:,kk)*AOAO(:,ll)                                     + &
          2.0_dp*dyAOAO(:,kk)*dyAOAO(:,ll)                             + &
          AOAO(:,kk)*dy2AOAO(:,ll)                                     + &
          dz2AOAO(:,kk)*AOAO(:,ll)                                     + &
          2.0_dp*dzAOAO(:,kk)*dzAOAO(:,ll)                             + &
          AOAO(:,kk)*dz2AOAO(:,ll))                                    + &
          0.5_dp*xcvtau(1:2*nl-1:2)*(&
          dxAOAO(:,kk)*dxAOAO(:,ll)                                    + &
          dyAOAO(:,kk)*dyAOAO(:,ll)                                    + &
          dzAOAO(:,kk)*dzAOAO(:,ll))))
          matxcmic(ll, kk) = matxcmic(kk, ll)
          ! beta
          matxcmic(cbdm+kk, cbdm+ll) = sum(weight(:)*(&
          ! GGA part: density, gradient-----------------------------------
          xcvrho(2:2*nl:2)*AOAO(:,kk)*AOAO(:,ll) + &
          2.0_dp*xcvsigma(3:3*nl:3)*(&
          drhob(:,1)*(dxAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dxAOAO(:,ll)) + &
          drhob(:,2)*(dyAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dyAOAO(:,ll)) + &
          drhob(:,3)*(dzAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dzAOAO(:,ll)))+ &
          xcvsigma(2:3*nl-1:3)*(&
          drhoa(:,1)*(dxAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dxAOAO(:,ll)) + &
          drhoa(:,2)*(dyAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dyAOAO(:,ll)) + &
          drhoa(:,3)*(dzAOAO(:,kk)*AOAO(:,ll)+AOAO(:,kk)*dzAOAO(:,ll)))+ &
          ! meta-GGA part: Laplacian, kinetic energy----------------------
          xcvlapl(2:2*nl:2)*(&
          dx2AOAO(:,kk)*AOAO(:,ll)                                     + &
          2.0_dp*dxAOAO(:,kk)*dxAOAO(:,ll)                             + &
          AOAO(:,kk)*dx2AOAO(:,ll)                                     + &
          dy2AOAO(:,kk)*AOAO(:,ll)                                     + &
          2.0_dp*dyAOAO(:,kk)*dyAOAO(:,ll)                             + &
          AOAO(:,kk)*dy2AOAO(:,ll)                                     + &
          dz2AOAO(:,kk)*AOAO(:,ll)                                     + &
          2.0_dp*dzAOAO(:,kk)*dzAOAO(:,ll)                             + &
          AOAO(:,kk)*dz2AOAO(:,ll))                                    + &
          0.5_dp*xcvtau(2:2*nl:2)*(&
          dxAOAO(:,kk)*dxAOAO(:,ll)                                    + &
          dyAOAO(:,kk)*dyAOAO(:,ll)                                    + &
          dzAOAO(:,kk)*dzAOAO(:,ll))))
          matxcmic(cbdm+ll, cbdm+kk) = matxcmic(cbdm+kk, cbdm+ll)
        end do
      end do
      matxc  = matxc + matxcmic
      exc    = exc + sum(weight(:)*(rhoa(:)+rhob(:))*excin(:))
    end do
  end subroutine Fock_XC_metaGGA

!------------------------------------------------------------
!> standard initialize Libxc
  subroutine Libxc_init()
    implicit none
    integer :: ii, vmajor, vminor, vmicro
    ! print out the version  
    call xc_f03_version(vmajor, vminor, vmicro)
    write(60,'("  -- Libxc version: ",I1,".",I1,".",I1)') vmajor, vminor, vmicro
    if (fx_id /= -1) then
      ! check the kind of fx_id and fc_id
      call xc_f03_func_init(x_func, fx_id, XC_POLARIZED)
      x_info = xc_f03_func_get_info(x_func)
      if (xc_f03_func_info_get_kind(x_info) == XC_EXCHANGE) then
        if (fc_id /= -1) then
          call xc_f03_func_init(c_func, fc_id, XC_POLARIZED)
          c_info = xc_f03_func_get_info(c_func)
          if (xc_f03_func_info_get_kind(c_info) /= XC_CORRELATION) call &
          terminate('Libxc_init: fc_id seems not a correlation functional')
        else
          call terminate('Libxc_init: missing data for correlation functional')
        end if
      else
        call terminate('Libxc_init: fx_id seems not an exchange functional')
      end if
      write(60,'(A13,I3,A13,A)') "  -- fx_id = ",fx_id,'; func name: ',&
      trim(xc_f03_func_info_get_name(x_info))
      write(60,'(A13,I3,A13,A)') "  -- fc_id = ",fc_id,'; func name: ',&
      trim(xc_f03_func_info_get_name(c_info))
      if (xc_f03_func_info_get_family(x_info) == XC_FAMILY_HYB_GGA .or. &
          xc_f03_func_info_get_family(x_info) == XC_FAMILY_HYB_LDA .or. &
          xc_f03_func_info_get_family(x_info) == XC_FAMILY_HYB_MGGA) then
        x_HF = xc_f03_hyb_exx_coef(x_func)
      else
        x_HF = 0.0_dp  ! pure functional
      end if
      write(60,'(A,E10.3)') "  -- x_HF = ", x_HF
    else if (fxc_id /= -1) then
      ! check the kind of fxc_id
      call xc_f03_func_init(xc_func, fxc_id, XC_POLARIZED)
      xc_info = xc_f03_func_get_info(xc_func)
      if (xc_f03_func_info_get_kind(xc_info) /= XC_EXCHANGE_CORRELATION) then
        call terminate(&
        'Libxc_init: fxc_id seems not an exchange-correlation functional')
      end if
      write(60,'(A13,I3,A13,A)') "  -- fxc_id = ",fxc_id,'; func name: ',&
      trim(xc_f03_func_info_get_name(xc_info))
      if (xc_f03_func_info_get_family(xc_info) == XC_FAMILY_HYB_GGA .or. &
          xc_f03_func_info_get_family(xc_info) == XC_FAMILY_HYB_LDA .or. &
          xc_f03_func_info_get_family(xc_info) == XC_FAMILY_HYB_MGGA) then
        x_HF = xc_f03_hyb_exx_coef(xc_func)
      else
        x_HF = 0.0_dp  ! pure functional
      end if
      write(60,'(A,E10.3)') "  -- x_HF = ", x_HF
    end if
  end subroutine Libxc_init

!------------------------------------------------------------
!> standard exit of Libxc
  subroutine Libxc_exit()
    implicit none
    if (fx_id /= -1) call xc_f03_func_end(x_func)
    if (fc_id /= -1) call xc_f03_func_end(c_func)
    if (fxc_id /= -1) call xc_f03_func_end(xc_func)
  end subroutine Libxc_exit

end module functional