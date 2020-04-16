PROGRAM ADV1D
USE plplot
IMPLICIT NONE

INTEGER :: nn, jj, n, j, g1, g2, i, m
REAL(kind=plflt) :: dt, dx, xc, yc, cfl, c, limiti(2)
REAL(kind=plflt), ALLOCATABLE :: u0(:), u(:), x(:), ya(:), yn(:)
CHARACTER (LEN=4) :: cc
CHARACTER (LEN=17) :: schema(5)=(/'centrato ordine 2', 'upwind ordine 1  ', &
 'centrato ordine 4', 'upwind ordine 3  ', 'downwind ordine 1'/)
LOGICAL lismouse

PRINT*,'Inserisci il numero di passi spaziali e il numero di Courant'
PRINT*,'Esempio: 30,0.1'
READ*,jj,cfl
ALLOCATE(u0(-2:jj+2),u(-2:jj+2),x(0:jj),ya(0:jj),yn(0:jj))

dx=1./jj
DO j=0,jj
  x(j)=dx*j
ENDDO

c=1.
dt=cfl*dx/c
u0=0.
!u0(jj/4)=1.
DO j=1,jj/8
  u0(j)=EXP(-(((j-jj/16.)/(jj/16.))**2)) ! Condizione iniziale
ENDDO
CALL SETCC(u0,jj)
limiti=(/-3.*ABS(MAXVAL(u0)),3.*ABS(MAXVAL(u0))/)

! Grafica
CALL plsdev('pdfc')
CALL plsfnam('adv.pdf')
!CALL plsetopt('portrait', ' ')
!CALL plparseopts(PL_PARSE_FULL)
CALL plscol0 (0, 255, 255, 255)
CALL plscol0 (1, 0, 0, 0)
CALL plinit

DO m=1,SIZE(schema)
  u=u0
  DO n=1,NINT(2.0_plflt/dt*c) ! 2 cicli completi

    CALL ANALITICA(u0,c,dx,dt,nn,jj,n,ya)
    CALL INTEGRA(u,c,dx,dt,nn,jj,n,m,yn)
    IF (MOD(n,10) == 0) THEN ! graphics every 10 steps
      CALL plcol0(1)
      CALL plenv(0._plflt, 1._plflt, limiti(1), limiti(2),0,0)
      CALL pllab('x', 'U(x,t)', &
       'Soluzione numerica dell''equazione di avvezione, schema '//TRIM(schema(m)))

      CALL plcol0(1)
      CALL plline(x,ya)

      CALL plcol0(1+m)
      CALL plline(x,yn)
    ENDIF

    CALL SETCC(u,jj)

  ENDDO

ENDDO
CALL plend()

END PROGRAM ADV1D


!*******************************************
SUBROUTINE SETCC(u,jj)
USE plplot
IMPLICIT NONE

INTEGER :: jj
REAL(kind=plflt) :: u(-2:jj+2)

u(-1)=u(jj)
u(-2)=u(jj-1)
u(jj+1)=u(0)
u(jj+2)=u(1)
END

!*******************************************
SUBROUTINE ANALITICA(u0,c,dx,dt,nn,jj,n,y)
USE plplot
IMPLICIT NONE

INTEGER :: nn, jj, n
REAL(kind=plflt) :: u0(-2:jj+2), c, dx, dt, y(0:jj)

REAL(kind=plflt) :: pesi(2)
INTEGER :: j, diff, part

diff=-MOD(INT(c*dt*n*jj),jj)
pesi(1)=c*dt*n*jj-AINT(c*dt*n*jj)
pesi(2)=1.-pesi(1)

DO j=0,jj
  part=MOD(j+diff+jj,jj)
  y(j)=SUM(pesi*u0(part-1:part))
ENDDO

END SUBROUTINE ANALITICA


!*******************************************
SUBROUTINE INTEGRA(u,c,dx,dt,nn,jj,n,m,y)
USE plplot
IMPLICIT NONE

INTEGER :: nn, jj, n, m
REAL(kind=plflt) :: u(-2:jj+2), c, dx, dt, y(0:jj)

INTEGER :: j
REAL :: ut(0:jj)

IF (m == 1) THEN ! centrato ordine 2
  DO j=0,jj
    y(j)=u(j)-dt/(2.*dx)*(u(j+1)-u(j-1))
  ENDDO
ELSE IF (m == 2) THEN ! upwind ordine 1
  DO j=0,jj
    y(j)=u(j)-dt/dx*(u(j)-u(j-1))
  ENDDO
ELSE IF (m == 3) THEN ! centrato ordine 4
  DO j=0,jj
    y(j)=u(j)-dt/(12.*dx)*(8.*(u(j+1)-u(j-1))-(u(j+2)-u(j-2)))
  ENDDO
ELSE IF (m == 4) THEN ! upwind ordine 3
  DO j=0,jj
    y(j)=u(j)-dt/(6.*dx)*(2.*u(j+1)+3.*u(j)-6*u(j-1)+u(j-2))
  ENDDO
ELSE IF (m == 5) THEN ! downwind
  DO j=0,jj
    y(j)=u(j)-dt/dx*(u(j+1)-u(j))
  ENDDO
ENDIF
u(0:jj)=y

END SUBROUTINE INTEGRA

