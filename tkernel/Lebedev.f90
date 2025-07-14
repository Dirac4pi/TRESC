!> @file Lebedev.f90
!!
!! @brief Generate Lebedev integral grid coordinates and weights
!!
!! @syntax Fortran 2008 free format
!!
!! @code UTF-8
!!
!! @author Christoph van Wuellen (modified by Dirac4pi)
module Lebedev
  implicit none

  contains

!> core grid generation subroutine
  pure subroutine gen_oh(code, num, x, y, z, w, a, b, v)
    implicit none
    integer, intent(in) :: code
    integer, intent(inout) :: num
    double precision, intent(inout) :: x(:), y(:), z(:), w(:)
    double precision, intent(inout) :: a, b, v
    double precision :: c

    select case(code)
    case(1)
      a = 1.0d0
      x(1) =  a; y(1) = 0.0d0; z(1) = 0.0d0; w(1) = v
      x(2) = -a; y(2) = 0.0d0; z(2) = 0.0d0; w(2) = v
      x(3) = 0.0d0; y(3) =  a; z(3) = 0.0d0; w(3) = v
      x(4) = 0.0d0; y(4) = -a; z(4) = 0.0d0; w(4) = v
      x(5) = 0.0d0; y(5) = 0.0d0; z(5) =  a; w(5) = v
      x(6) = 0.0d0; y(6) = 0.0d0; z(6) = -a; w(6) = v
      num = num + 6

    case(2)
      a = sqrt(0.5d0)
      x( 1) = 0.0d0; y( 1) =  a; z( 1) =  a; w( 1) = v
      x( 2) = 0.0d0; y( 2) = -a; z( 2) =  a; w( 2) = v
      x( 3) = 0.0d0; y( 3) =  a; z( 3) = -a; w( 3) = v
      x( 4) = 0.0d0; y( 4) = -a; z( 4) = -a; w( 4) = v
      x( 5) =  a; y( 5) = 0.0d0; z( 5) =  a; w( 5) = v
      x( 6) = -a; y( 6) = 0.0d0; z( 6) =  a; w( 6) = v
      x( 7) =  a; y( 7) = 0.0d0; z( 7) = -a; w( 7) = v
      x( 8) = -a; y( 8) = 0.0d0; z( 8) = -a; w( 8) = v
      x( 9) =  a; y( 9) =  a; z( 9) = 0.0d0; w( 9) = v
      x(10) = -a; y(10) =  a; z(10) = 0.0d0; w(10) = v
      x(11) =  a; y(11) = -a; z(11) = 0.0d0; w(11) = v
      x(12) = -a; y(12) = -a; z(12) = 0.0d0; w(12) = v
      num = num + 12

    case(3)
      a = sqrt(1.0d0 / 3.0d0)
      x(1) =  a; y(1) =  a; z(1) =  a; w(1) = v
      x(2) = -a; y(2) =  a; z(2) =  a; w(2) = v
      x(3) =  a; y(3) = -a; z(3) =  a; w(3) = v
      x(4) = -a; y(4) = -a; z(4) =  a; w(4) = v
      x(5) =  a; y(5) =  a; z(5) = -a; w(5) = v
      x(6) = -a; y(6) =  a; z(6) = -a; w(6) = v
      x(7) =  a; y(7) = -a; z(7) = -a; w(7) = v
      x(8) = -a; y(8) = -a; z(8) = -a; w(8) = v
      num = num + 8

    case(4)
      b = sqrt(1.0d0 - 2.0d0 * a * a)
      x( 1) =  a; y( 1) =  a; z( 1) =  b; w( 1) = v
      x( 2) = -a; y( 2) =  a; z( 2) =  b; w( 2) = v
      x( 3) =  a; y( 3) = -a; z( 3) =  b; w( 3) = v
      x( 4) = -a; y( 4) = -a; z( 4) =  b; w( 4) = v
      x( 5) =  a; y( 5) =  a; z( 5) = -b; w( 5) = v
      x( 6) = -a; y( 6) =  a; z( 6) = -b; w( 6) = v
      x( 7) =  a; y( 7) = -a; z( 7) = -b; w( 7) = v
      x( 8) = -a; y( 8) = -a; z( 8) = -b; w( 8) = v
      x( 9) =  a; y( 9) =  b; z( 9) =  a; w( 9) = v
      x(10) = -a; y(10) =  b; z(10) =  a; w(10) = v
      x(11) =  a; y(11) = -b; z(11) =  a; w(11) = v
      x(12) = -a; y(12) = -b; z(12) =  a; w(12) = v
      x(13) =  a; y(13) =  b; z(13) = -a; w(13) = v
      x(14) = -a; y(14) =  b; z(14) = -a; w(14) = v
      x(15) =  a; y(15) = -b; z(15) = -a; w(15) = v
      x(16) = -a; y(16) = -b; z(16) = -a; w(16) = v
      x(17) =  b; y(17) =  a; z(17) =  a; w(17) = v
      x(18) = -b; y(18) =  a; z(18) =  a; w(18) = v
      x(19) =  b; y(19) = -a; z(19) =  a; w(19) = v
      x(20) = -b; y(20) = -a; z(20) =  a; w(20) = v
      x(21) =  b; y(21) =  a; z(21) = -a; w(21) = v
      x(22) = -b; y(22) =  a; z(22) = -a; w(22) = v
      x(23) =  b; y(23) = -a; z(23) = -a; w(23) = v
      x(24) = -b; y(24) = -a; z(24) = -a; w(24) = v
      num = num + 24

    case(5)
      b = sqrt(1.0d0 - a * a)
      x( 1) =  a; y( 1) =  b; z( 1) = 0.0d0; w( 1) = v
      x( 2) = -a; y( 2) =  b; z( 2) = 0.0d0; w( 2) = v
      x( 3) =  a; y( 3) = -b; z( 3) = 0.0d0; w( 3) = v
      x( 4) = -a; y( 4) = -b; z( 4) = 0.0d0; w( 4) = v
      x( 5) =  b; y( 5) =  a; z( 5) = 0.0d0; w( 5) = v
      x( 6) = -b; y( 6) =  a; z( 6) = 0.0d0; w( 6) = v
      x( 7) =  b; y( 7) = -a; z( 7) = 0.0d0; w( 7) = v
      x( 8) = -b; y( 8) = -a; z( 8) = 0.0d0; w( 8) = v
      x( 9) =  a; y( 9) = 0.0d0; z( 9) =  b; w( 9) = v
      x(10) = -a; y(10) = 0.0d0; z(10) =  b; w(10) = v
      x(11) =  a; y(11) = 0.0d0; z(11) = -b; w(11) = v
      x(12) = -a; y(12) = 0.0d0; z(12) = -b; w(12) = v
      x(13) =  b; y(13) = 0.0d0; z(13) =  a; w(13) = v
      x(14) = -b; y(14) = 0.0d0; z(14) =  a; w(14) = v
      x(15) =  b; y(15) = 0.0d0; z(15) = -a; w(15) = v
      x(16) = -b; y(16) = 0.0d0; z(16) = -a; w(16) = v
      x(17) = 0.0d0; y(17) =  a; z(17) =  b; w(17) = v
      x(18) = 0.0d0; y(18) = -a; z(18) =  b; w(18) = v
      x(19) = 0.0d0; y(19) =  a; z(19) = -b; w(19) = v
      x(20) = 0.0d0; y(20) = -a; z(20) = -b; w(20) = v
      x(21) = 0.0d0; y(21) =  b; z(21) =  a; w(21) = v
      x(22) = 0.0d0; y(22) = -b; z(22) =  a; w(22) = v
      x(23) = 0.0d0; y(23) =  b; z(23) = -a; w(23) = v
      x(24) = 0.0d0; y(24) = -b; z(24) = -a; w(24) = v
      num = num + 24

    case(6)
      c = sqrt(1.0d0 - a * a - b * b)
      x( 1) =  a; y( 1) =  b; z( 1) =  c; w( 1) = v
      x( 2) = -a; y( 2) =  b; z( 2) =  c; w( 2) = v
      x( 3) =  a; y( 3) = -b; z( 3) =  c; w( 3) = v
      x( 4) = -a; y( 4) = -b; z( 4) =  c; w( 4) = v
      x( 5) =  a; y( 5) =  b; z( 5) = -c; w( 5) = v
      x( 6) = -a; y( 6) =  b; z( 6) = -c; w( 6) = v
      x( 7) =  a; y( 7) = -b; z( 7) = -c; w( 7) = v
      x( 8) = -a; y( 8) = -b; z( 8) = -c; w( 8) = v
      x( 9) =  a; y( 9) =  c; z( 9) =  b; w( 9) = v
      x(10) = -a; y(10) =  c; z(10) =  b; w(10) = v
      x(11) =  a; y(11) = -c; z(11) =  b; w(11) = v
      x(12) = -a; y(12) = -c; z(12) =  b; w(12) = v
      x(13) =  a; y(13) =  c; z(13) = -b; w(13) = v
      x(14) = -a; y(14) =  c; z(14) = -b; w(14) = v
      x(15) =  a; y(15) = -c; z(15) = -b; w(15) = v
      x(16) = -a; y(16) = -c; z(16) = -b; w(16) = v
      x(17) =  b; y(17) =  a; z(17) =  c; w(17) = v
      x(18) = -b; y(18) =  a; z(18) =  c; w(18) = v
      x(19) =  b; y(19) = -a; z(19) =  c; w(19) = v
      x(20) = -b; y(20) = -a; z(20) =  c; w(20) = v
      x(21) =  b; y(21) =  a; z(21) = -c; w(21) = v
      x(22) = -b; y(22) =  a; z(22) = -c; w(22) = v
      x(23) =  b; y(23) = -a; z(23) = -c; w(23) = v
      x(24) = -b; y(24) = -a; z(24) = -c; w(24) = v
      x(25) =  b; y(25) =  c; z(25) =  a; w(25) = v
      x(26) = -b; y(26) =  c; z(26) =  a; w(26) = v
      x(27) =  b; y(27) = -c; z(27) =  a; w(27) = v
      x(28) = -b; y(28) = -c; z(28) =  a; w(28) = v
      x(29) =  b; y(29) =  c; z(29) = -a; w(29) = v
      x(30) = -b; y(30) =  c; z(30) = -a; w(30) = v
      x(31) =  b; y(31) = -c; z(31) = -a; w(31) = v
      x(32) = -b; y(32) = -c; z(32) = -a; w(32) = v
      x(33) =  c; y(33) =  a; z(33) =  b; w(33) = v
      x(34) = -c; y(34) =  a; z(34) =  b; w(34) = v
      x(35) =  c; y(35) = -a; z(35) =  b; w(35) = v
      x(36) = -c; y(36) = -a; z(36) =  b; w(36) = v
      x(37) =  c; y(37) =  a; z(37) = -b; w(37) = v
      x(38) = -c; y(38) =  a; z(38) = -b; w(38) = v
      x(39) =  c; y(39) = -a; z(39) = -b; w(39) = v
      x(40) = -c; y(40) = -a; z(40) = -b; w(40) = v
      x(41) =  c; y(41) =  b; z(41) =  a; w(41) = v
      x(42) = -c; y(42) =  b; z(42) =  a; w(42) = v
      x(43) =  c; y(43) = -b; z(43) =  a; w(43) = v
      x(44) = -c; y(44) = -b; z(44) =  a; w(44) = v
      x(45) =  c; y(45) =  b; z(45) = -a; w(45) = v
      x(46) = -c; y(46) =  b; z(46) = -a; w(46) = v
      x(47) =  c; y(47) = -b; z(47) = -a; w(47) = v
      x(48) = -c; y(48) = -b; z(48) = -a; w(48) = v
      num = num + 48
    end select
  end subroutine gen_oh

!> generate 230 grid points in Lebedev scheme
  pure subroutine LD0230(x, y, z, w, n)
    implicit none
    double precision, intent(out) :: x(:), y(:), z(:), w(:)
    integer, intent(out) :: n
    double precision :: a, b, v

    n = 1
    v = -0.5522639919727325d-1
    call gen_oh(1, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    v = 0.4450274607445226d-2
    call gen_oh(3, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.4492044687397611d+0
    v = 0.4496841067921404d-2
    call gen_oh(4, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.2520419490210201d+0
    v = 0.5049153450478750d-2
    call gen_oh(4, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.6981906658447242d+0
    v = 0.3976408018051883d-2
    call gen_oh(4, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.6587405243460960d+0
    v = 0.4401400650381014d-2
    call gen_oh(4, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.4038544050097660d-1
    v = 0.1724544350544401d-1
    call gen_oh(4, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.5823842309715585d+0
    v = 0.4231083095357343d-2
    call gen_oh(5, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.3545877390518688d+0
    v = 0.5198069864064399d-2
    call gen_oh(5, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.2272181808998187d+0
    b = 0.4864661535886647d+0
    v = 0.4695720972568883d-2
    call gen_oh(6, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    n = n - 1
  end subroutine LD0230

!> generate 302 grid points in Lebedev scheme
  pure subroutine LD0302(x, y, z, w, n)
    implicit none
    double precision, intent(out) :: x(:), y(:), z(:), w(:)
    integer, intent(out) :: n
    double precision :: a, b, v

    n = 1
    v = 0.8545911725128148d-3
    call gen_oh(1, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    v = 0.3599119285025571d-2
    call gen_oh(3, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.3515640345570105d+0
    v = 0.3449788424305883d-2
    call gen_oh(4, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.6566329410219612d+0
    v = 0.3604822601419882d-2
    call gen_oh(4, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.4729054132581005d+0
    v = 0.3576729661743367d-2
    call gen_oh(4, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.9618308522614784d-1
    v = 0.2352101413689164d-2
    call gen_oh(4, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.2219645236294178d+0
    v = 0.3108953122413675d-2
    call gen_oh(4, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.7011766416089545d+0
    v = 0.3650045807677255d-2
    call gen_oh(4, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.2644152887060663d+0
    v = 0.2982344963171804d-2
    call gen_oh(5, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.5718955891878961d+0
    v = 0.3600820932216460d-2
    call gen_oh(5, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.2510034751770465d+0
    b = 0.8000727494073952d+0
    v = 0.3571540554273387d-2
    call gen_oh(6, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.1233548532583327d+0
    b = 0.4127724083168531d+0
    v = 0.3392312205006170d-2
    call gen_oh(6, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    n = n - 1
  end subroutine LD0302

!> generate 434 grid points in Lebedev scheme
  pure subroutine LD0434(x, y, z, w, n)
    implicit none
    double precision, intent(out) :: x(:), y(:), z(:), w(:)
    integer, intent(out) :: n
    double precision :: a, b, v

    n = 1
    v = 0.5265897968224436d-3
    call gen_oh(1, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    v = 0.2548219972002607d-2
    call gen_oh(2, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    v = 0.2512317418927307d-2
    call gen_oh(3, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.6909346307509111d+0
    v = 0.2530403801186355d-2
    call gen_oh(4, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.1774836054609158d+0
    v = 0.2014279020918528d-2
    call gen_oh(4, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.4914342637784746d+0
    v = 0.2501725168402936d-2
    call gen_oh(4, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.6456664707424256d+0
    v = 0.2513267174597564d-2
    call gen_oh(4, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.2861289010307638d+0
    v = 0.2302694782227416d-2
    call gen_oh(4, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.7568084367178018d-1
    v = 0.1462495621594614d-2
    call gen_oh(4, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.3927259763368002d+0
    v = 0.2445373437312980d-2
    call gen_oh(4, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.8818132877794288d+0
    v = 0.2417442375638981d-2
    call gen_oh(5, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.9776428111182649d+0
    v = 0.1910951282179532d-2
    call gen_oh(5, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.2054823696403044d+0
    b = 0.8689460322872412d+0
    v = 0.2416930044324775d-2
    call gen_oh(6, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.5905157048925271d+0
    b = 0.7999278543857286d+0
    v = 0.2512236854563495d-2
    call gen_oh(6, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.5550152361076807d+0
    b = 0.7717462626915901d+0
    v = 0.2496644054553086d-2
    call gen_oh(6, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.9371809858553722d+0
    b = 0.3344363145343455d+0
    v = 0.2236607760437849d-2
    call gen_oh(6, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    n = n - 1
  end subroutine LD0434

!> generate 590 grid points in Lebedev scheme
  pure subroutine LD0590(x, y, z, w, n)
    implicit none
    double precision, intent(out) :: x(:), y(:), z(:), w(:)
    integer, intent(out) :: n
    double precision :: a, b, v

    n = 1
    v = 0.3095121295306187d-3
    call gen_oh(1, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    v = 0.1852379698597489d-2
    call gen_oh(3, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.7040954938227469d+0
    v = 0.1871790639277744d-2
    call gen_oh(4, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.6807744066455243d+0
    v = 0.1858812585438317d-2
    call gen_oh(4, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.6372546939258752d+0
    v = 0.1852028828296213d-2
    call gen_oh(4, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.5044419707800358d+0
    v = 0.1846715956151242d-2
    call gen_oh(4, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.4215761784010967d+0
    v = 0.1818471778162769d-2
    call gen_oh(4, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.3317920736472123d+0
    v = 0.1749564657281154d-2
    call gen_oh(4, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.2384736701421887d+0
    v = 0.1617210647254411d-2
    call gen_oh(4, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.1459036449157763d+0
    v = 0.1384737234851692d-2
    call gen_oh(4, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.6095034115507196d-1
    v = 0.9764331165051050d-3
    call gen_oh(4, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.6116843442009876d+0
    v = 0.1857161196774078d-2
    call gen_oh(5, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.3964755348199858d+0
    v = 0.1705153996395864d-2
    call gen_oh(5, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.1724782009907724d+0
    v = 0.1300321685886048d-2
    call gen_oh(5, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.5610263808622060d+0
    b = 0.3518280927733519d+0
    v = 0.1842866472905286d-2
    call gen_oh(6, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.4742392842551980d+0
    b = 0.2634716655937950d+0
    v = 0.1802658934377451d-2
    call gen_oh(6, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.5984126497885380d+0
    b = 0.1816640840360209d+0
    v = 0.1849830560443660d-2
    call gen_oh(6, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.3791035407695563d+0
    b = 0.1720795225656878d+0
    v = 0.1713904507106709d-2
    call gen_oh(6, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.2778673190586244d+0
    b = 0.8213021581932511d-1
    v = 0.1555213603396808d-2
    call gen_oh(6, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.5033564271075117d+0
    b = 0.8999205842074875d-1
    v = 0.1802239128008525d-2
    call gen_oh(6, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    n = n - 1
  end subroutine LD0590

!> generate 974 grid points in Lebedev scheme
  pure subroutine LD0974(x, y, z, w, n)
    implicit none
    double precision, intent(out) :: x(:), y(:), z(:), w(:)
    integer, intent(out) :: n
    double precision :: a, b, v

    n = 1
    v = 0.1438294190527431d-3
    call gen_oh(1, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    v = 0.1125772288287004d-2
    call gen_oh(3, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.4292963545341347d-1
    v = 0.4948029341949241d-3
    call gen_oh(4, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.1051426854086404d+0
    v = 0.7357990109125470d-3
    call gen_oh(4, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.1750024867623087d+0
    v = 0.8889132771304384d-3
    call gen_oh(4, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.2477653379650257d+0
    v = 0.9888347838921435d-3
    call gen_oh(4, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.3206567123955957d+0
    v = 0.1053299681709471d-2
    call gen_oh(4, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.3916520749849983d+0
    v = 0.1092778807014578d-2
    call gen_oh(4, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.4590825874187624d+0
    v = 0.1114389394063227d-2
    call gen_oh(4, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.5214563888415861d+0
    v = 0.1123724788051555d-2
    call gen_oh(4, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.6253170244654199d+0
    v = 0.1125239325243814d-2
    call gen_oh(4, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.6637926744523170d+0
    v = 0.1126153271815905d-2
    call gen_oh(4, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.6910410398498301d+0
    v = 0.1130286931123841d-2
    call gen_oh(4, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.7052907007457760d+0
    v = 0.1134986534363955d-2
    call gen_oh(4, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.1236686762657990d+0
    v = 0.6823367927109931d-3
    call gen_oh(5, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.2940777114468387d+0
    v = 0.9454158160447096d-3
    call gen_oh(5, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.4697753849207649d+0
    v = 0.1074429975385679d-2
    call gen_oh(5, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.6334563241139567d+0
    v = 0.1129300086569132d-2
    call gen_oh(5, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.5974048614181342d-1
    b = 0.2029128752777523d+0
    v = 0.8436884500901954d-3
    call gen_oh(6, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.1375760408473636d+0
    b = 0.4602621942484054d+0
    v = 0.1075255720448885d-2
    call gen_oh(6, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.3391016526336286d+0
    b = 0.5030673999662036d+0
    v = 0.1108577236864462d-2
    call gen_oh(6, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.1271675191439820d+0
    b = 0.2817606422442134d+0
    v = 0.9566475323783357d-3
    call gen_oh(6, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.2693120740413512d+0
    b = 0.4331561291720157d+0
    v = 0.1080663250717391d-2
    call gen_oh(6, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.1419786452601918d+0
    b = 0.6256167358580814d+0
    v = 0.1126797131196295d-2
    call gen_oh(6, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.6709284600738255d-1
    b = 0.3798395216859157d+0
    v = 0.1022568715358061d-2
    call gen_oh(6, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.7057738183256172d-1
    b = 0.5517505421423520d+0
    v = 0.1108960267713108d-2
    call gen_oh(6, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.2783888477882155d+0
    b = 0.6029619156159187d+0
    v = 0.1122790653435766d-2
    call gen_oh(6, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.1979578938917407d+0
    b = 0.3589606329589096d+0
    v = 0.1032401847117460d-2
    call gen_oh(6, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.2087307061103274d+0
    b = 0.5348666438135476d+0
    v = 0.1107249382283854d-2
    call gen_oh(6, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    a = 0.4055122137872836d+0
    b = 0.5674997546074373d+0
    v = 0.1121780048519972d-2
    call gen_oh(6, n, x(n:), y(n:), z(n:), w(n:), a, b, v)

    n = n - 1
  end subroutine LD0974
end module Lebedev