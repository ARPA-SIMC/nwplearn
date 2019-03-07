PROGRAM OSCILL
USE plplot
IMPLICIT NONE

INTEGER :: nn, jj, m, n, j, g1, g2, i
REAL(kind=plflt) :: dt, dx, om, xc, yc, pert(3)=(/1._plflt, .1_plflt, -1._plflt/)
COMPLEX(kind=plflt) :: u, ua, c
REAL(kind=plflt), ALLOCATABLE :: x(:), y(:)
!INTEGER :: PGOPEN, PGCURS
CHARACTER (LEN=4) :: cc
LOGICAL lismouse

PRINT*,'Inserisci il numero di passi temporali'
READ*,nn
ALLOCATE(x(0:nn),y(0:nn))

! Grafica
CALL plsdev('pdfc')
CALL plsfnam('oscill.pdf')
CALL plsetopt('portrait', ' ')
CALL plparseopts(PL_PARSE_FULL)
CALL plscol0 (0, 255, 255, 255)
CALL plscol0 (1, 0, 0, 0)
CALL plinit
!g1=PGOPEN('/XWIN')
!CALL PGASK(.FALSE.)
!CALL PGQINF('CURSOR',cc,i)
!lismouse = (cc == 'YES')
lismouse = .FALSE.

DO m=1,3
u=(1.,0.) ! Condizione iniziale
dt=1./nn
om=3.*8.*ATAN(1.) ! 3 oscillazioni (3*2*pi)
c=dt*CMPLX(0.,om) ! i*omega*dt

CALL plenv(0._plflt, 1._plflt, &
 -3._plflt*ABS(REAL(u,kind=plflt)),3._plflt*ABS(REAL(u,kind=plflt)),0,0)
CALL pllab('t', 'U(t)', 'Soluzione numerica dell''equazione dell''oscillatore')
DO n=0,nn
  x(n)=dt*n
ENDDO

CALL plwidth(4._plflt)
CALL plcol0(1)
CALL LEGENDA(u,0)
CALL ANALITICA(u,c,nn,y)
CALL plline(x,y)
CALL plwidth(1._plflt)

DO n=1,4
  u=(1.,0.) ! Condizione iniziale a n=0

  IF (n == 4 .AND. lismouse) THEN ! Leap frog, input della c.i. a n=1
!    i=PGCURS(xc,yc,cc)
!    g2=PGOPEN('/XWIN')
!    CALL PGASK(.FALSE.)
!    CALL PGSLCT(g2)
!    CALL PGENV(-2., 2.,-2., 2.,1,1)
!    CALL PGLAB('Parte reale', 'Parte immaginaria', &
!     'Scegli la condizione iniziale a n=1 per lo schema a 3 livelli')
!    ua=u*EXP(c)
!    CALL PGMOVE(-2.*REAL(ua),-2.*AIMAG(ua))
!    CALL PGDRAW(2.*REAL(ua),2.*AIMAG(ua)) 
!    CALL PGSFS(2)
!    CALL PGCIRC(0.,0.,1.)
!    CALL PGCIRC(0.,0.,2.)
!    CALL PGCIRC(REAL(ua),AIMAG(ua),.05)
!    i=PGCURS(xc,yc,cc)
!    IF (cc == 'D') THEN ! tasto centrale del mouse
!      ua=CMPLX(xc,yc)
!    ENDIF
!    CALL PGCLOS()
!    CALL PGSLCT(g1)
  ELSE
    ua=pert(m)*u*EXP(c) ! se non ho il mouse uso la c.i. analitica
  ENDIF

  CALL plcol0(1+n)
  CALL LEGENDA(u,n)
  CALL INTEGRA(u,ua,c,nn,n,y)
  CALL plline(x,y)
ENDDO
CALL plcol0(1)
ENDDO

!IF (lismouse) i=PGCURS(xc,yc,cc)
CALL plend()

END PROGRAM OSCILL


!*******************************************
SUBROUTINE LEGENDA(u,schema)
USE plplot
IMPLICIT NONE

COMPLEX(kind=plflt) :: u
INTEGER :: schema

REAL(kind=plflt) :: x(0:1),y(0:1)
INTEGER :: lw
CHARACTER (LEN=15) :: leg(0:4)=(/ &
 'Sol. analitica', &
 'Eulero        ', &
 'Backward      ', &
 'Trapezoidale  ', &
 'Leap frog     ' /)


y(0:1)=-2*REAL(u)+REAL(u)/5.*schema
x(0)=0.
x(1)=0.1
CALL plline(x,y)
!CALL PGQLW(lw)
!CALL plwid(2)
CALL plptex(x(1), y(1), 1._plflt, 0._plflt, 0._plflt, leg(schema))
!CALL PGSLW(lw)

END SUBROUTINE LEGENDA


!*******************************************
SUBROUTINE ANALITICA(u,c,nn,y)
USE plplot
IMPLICIT NONE

COMPLEX(kind=plflt) :: u, c
INTEGER :: nn
REAL(kind=plflt) :: y(0:nn)

COMPLEX(kind=plflt) :: ua
INTEGER :: n

y(0)=REAL(u, kind=plflt)
DO n=1,nn
  ua=u*EXP(c*n)
  y(n)=REAL(ua, kind=plflt)
ENDDO

END SUBROUTINE ANALITICA


!*******************************************
SUBROUTINE INTEGRA(u,u1,c,nn,schema,y)
USE plplot
IMPLICIT NONE

COMPLEX(kind=plflt) :: u, u1, c
INTEGER :: nn, schema
REAL(kind=plflt) :: y(0:nn)

COMPLEX(kind=plflt) :: un
INTEGER :: n
REAL(kind=plflt) :: a, b

IF (schema < 4) THEN ! a 2 livelli
  IF (schema == 1) THEN ! Eulero
    a=1.
  ELSE IF (schema == 2) THEN ! Implicito
    a=0.
  ELSE IF (schema == 3) THEN ! Trapezoidale
    a=0.5
  ENDIF
  b=1.-a

  y(0)=REAL(u)
  DO n=1,nn
    u=(u+c*(a*u))/(CMPLX(1.,0.)-c*b)
    y(n)=REAL(u)
  ENDDO

ELSE ! a 3 livelli, Leap frog
  y(0)=REAL(u)
  y(1)=REAL(u1)
  DO n=2,nn
    un=u+2*c*u1
    u=u1
    u1=un
    y(n)=REAL(un)
  ENDDO

ENDIF

END SUBROUTINE INTEGRA

