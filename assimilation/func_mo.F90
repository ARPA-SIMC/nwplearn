MODULE func_mo
USE kinds
USE missing_values
!USE phys_const
IMPLICIT NONE

REAL(fp_d),PARAMETER :: pi = 3.14159265358979

INTERFACE random_gaussian
  MODULE PROCEDURE random_gaussian_0, random_gaussian_1, random_gaussian_2
END INTERFACE random_gaussian

PUBLIC func_autocorr, func_2dfield, &
 random_gaussian_f, random_gaussian, init_random_seed

CONTAINS

!> Compute an autocorrelation-like function of the argument.
!! Function supported at the moment (\a func argument) are
!! gauss (Gaussian), exp (exponential and delta (Dirac).
ELEMENTAL FUNCTION func_autocorr(x, func) RESULT(y)
REAL(fp_d),INTENT(in) :: x !< argument of the function
CHARACTER(len=*),INTENT(in) :: func !< type of function
REAL(fp_d) :: y

SELECT CASE(func)
CASE('gauss')
  y = EXP(-0.5_fp_d*x*x)
CASE('exp')
  y = EXP(-x)
CASE('delta')
  IF (x == 0.0_fp_d) THEN
    y = 1.0_fp_d
  ELSE
    y = 0.0_fp_d
  ENDIF
CASE default
  y = dmiss
END SELECT

END FUNCTION func_autocorr


!> Compute a function suitable for a 2-d scalar field.
!! Function supported at the moment (\a func argument) are
!! fourier (Fourier).
ELEMENTAL FUNCTION func_2dfield(x, y, func) RESULT(z)
REAL(fp_d),INTENT(in) :: x, y
CHARACTER(len=*),INTENT(in) :: func
REAL(fp_d) :: z

SELECT CASE(func)
CASE('fourier')
  z = SIN(pi*x)*SIN(pi*y)
CASE default
  z = 3.
END SELECT

END FUNCTION func_2dfield


!> Return a random number taken from a Gaussian distribution with zero average
!! and unit variance.
!! Scalar version.
RECURSIVE FUNCTION random_gaussian_f() RESULT(rg)
REAL(fp_d) :: rg

REAL(fp_d) :: r1, r2
REAL(fp_d),SAVE :: rgnext
LOGICAL,SAVE :: mustdo=.TRUE.

IF (mustdo) THEN
  CALL RANDOM_NUMBER(r1)
  CALL RANDOM_NUMBER(r2)
  IF (r1 <= 0.0_fp_d) THEN
    rg = random_gaussian_f()
    RETURN
  ENDIF

  rg = SQRT(-2.0_fp_d*LOG(r1))*COS(2.0_fp_d*pi*r2)
  rgnext = SQRT(-2.0_fp_d*LOG(r1))*SIN(2.0_fp_d*pi*r2)
  mustdo = .FALSE.
ELSE
  rg = rgnext
  mustdo = .TRUE.
ENDIF

END FUNCTION random_gaussian_f


!> Compute a random number taken from a Gaussian distribution with zero average
!! and unit variance.
!! Scalar version.
SUBROUTINE random_gaussian_0(rg)
REAL(fp_d) :: rg
rg = random_gaussian_f()
END SUBROUTINE random_gaussian_0


!> Compute a random number taken from a Gaussian distribution with zero average
!! and unit variance.
!! 1-d array version.
SUBROUTINE random_gaussian_1(rg)
REAL(fp_d) :: rg(:)
INTEGER :: i
DO i = 1, SIZE(rg)
  rg(i) = random_gaussian_f()
ENDDO
END SUBROUTINE random_gaussian_1


!> Compute a random number taken from a Gaussian distribution with zero average
!! and unit variance.
!! 2-d array version.
SUBROUTINE random_gaussian_2(rg)
REAL(fp_d) :: rg(:,:)
INTEGER :: i, j
DO j = 1, SIZE(rg, 2)
  DO i = 1, SIZE(rg, 1)
    rg(i,j) = random_gaussian_f()
  ENDDO
ENDDO
END SUBROUTINE random_gaussian_2


!> Initialize the Fortran random seed with a non-reproducible number
!! based on current system time.
!! It affects the intrinsic \a RANDOM_NUMBER subroutine as well as
!! the Gaussian random function/subroutines in this module.
SUBROUTINE init_random_seed()
INTEGER :: i, n, clock
INTEGER, DIMENSION(:), ALLOCATABLE :: seed

CALL RANDOM_SEED(size = n)
ALLOCATE(seed(n))

CALL SYSTEM_CLOCK(COUNT=clock)

seed = clock + 37 * (/ (i - 1, i = 1, n) /)
CALL RANDOM_SEED(PUT = seed)

DEALLOCATE(seed)
END SUBROUTINE init_random_seed

END MODULE func_mo
