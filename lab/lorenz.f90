MODULE lorenz_m
USE plplot
IMPLICIT NONE

TYPE lorenz_t
  PRIVATE
  REAL,ALLOCATABLE :: x(:), y(:), z(:)
  REAL :: t=0.
  REAL :: a=5., b=15., c=1.
  INTEGER :: n=0
  CONTAINS
  PROCEDURE :: advance => lorenz_advance
  PROCEDURE :: delete => lorenz_delete
!  PROCEDURE :: plot => lorenz_plot
END TYPE lorenz_t

INTERFACE lorenz_t
  MODULE PROCEDURE lorenz_new
END INTERFACE lorenz_t

PRIVATE
PUBLIC lorenz_t, lorenz_plot

CONTAINS

FUNCTION lorenz_new(x, y, z, n, a, b, c) RESULT(this)
REAL,INTENT(in) :: x, y, z
INTEGER,INTENT(in) :: n
REAL,INTENT(in),OPTIONAL :: a, b, c

TYPE(lorenz_t) :: this

ALLOCATE(this%x(0:n),this%y(0:n),this%z(0:n))
this%x(0) = x
this%y(0) = y
this%z(0) = z
IF (PRESENT(a)) this%a = a
IF (PRESENT(b)) this%b = b
IF (PRESENT(c)) this%c = c

END FUNCTION lorenz_new


SUBROUTINE lorenz_advance(this, dt)
CLASS(lorenz_t),INTENT(inout) :: this
REAL,INTENT(in) :: dt
REAL :: x, y, z

this%x(this%n+1) = this%x(this%n) + (-this%a*this%x(this%n)*dt) + &
 (this%a*this%y(this%n)*dt)
this%y(this%n+1) = this%y(this%n) + ( this%b*this%x(this%n)*dt) - &
 (       this%y(this%n)*dt) - (this%z(this%n)*this%x(this%n)*dt)
this%z(this%n+1) = this%z(this%n) + (-this%c*this%z(this%n)*dt) + &
 (this%x(this%n)*this%y(this%n)*dt)

this%n = this%n + 1
this%t = this%t + dt

END SUBROUTINE lorenz_advance


SUBROUTINE lorenz_delete(this)
CLASS(lorenz_t),INTENT(out) :: this

END SUBROUTINE lorenz_delete


SUBROUTINE lorenz_plot(this)
TYPE(lorenz_t) :: this(:)

INTEGER :: i

CALL pllab('x', 'z', 'Lorenz model')
DO i = 1, SIZE(this)
  CALL plcol0(MOD(i,14)+1)!14->15
  CALL plsym(REAL(this(i)%x(0:0),kind=plflt), REAL(this(i)%z(0:0),kind=plflt), &
   844)
  CALL plline(REAL(this(i)%x(:this(i)%n),kind=plflt), &
   REAL(this(i)%z(:this(i)%n),kind=plflt))
  CALL plsym(REAL(this(i)%x(this(i)%n:this(i)%n),kind=plflt), &
   REAL(this(i)%z(this(i)%n:this(i)%n),kind=plflt), &
   845)

ENDDO

END SUBROUTINE lorenz_plot

END MODULE lorenz_m


PROGRAM lorenz
USE lorenz_m
USE plplot
IMPLICIT NONE
INTEGER,PARAMETER :: nist=20, nstep=500
INTEGER :: n, m
REAL :: dt=0.02
TYPE(lorenz_t) :: lor_mod(nist)

!CALL plsdev('pngc')
!CALL plsfam(1, 1, 1)
!CALL plsfnam('lorenz-%n.png')
!CALL plsetopt('portrait', ' ')
CALL plsdev('pdfc')
CALL plsfnam('lorenz.pdf')
CALL plsetopt('portrait', ' ')
!  CALL plparseopts(PL_PARSE_FULL)
CALL plscol0 (0, 255, 255, 255)
CALL plscol0 (1, 0, 0, 0)
CALL plinit()


DO m = 1, nist
!    CALL PGSCI(MOD(m,15))
  lor_mod(m) = lorenz_t(REAL(m)-15., 1., REAL(m), n=nstep)
ENDDO
 
DO n = 1, nstep
  PRINT*,'step',n
  DO m = 1, nist
    CALL lor_mod(m)%advance(dt)
  ENDDO
  CALL plcol0(1)
  CALL plenv(REAL(-15.,kind=plflt), REAL(15.,kind=plflt), &
   REAL(0.,kind=plflt), REAL(30.,kind=plflt), 0, 0)
  CALL lorenz_plot(lor_mod)

ENDDO

DO m = 1, nist
  CALL lor_mod(m)%delete()
ENDDO

CALL plend()

END PROGRAM lorenz
