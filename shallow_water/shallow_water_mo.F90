MODULE shallow_water_mo
  USE plplot
  IMPLICIT NONE

  INTEGER :: nx, ny, nt
  INTEGER :: nxu, nyv
  REAL :: g=9.8, h=1., dt, dx, dy
  REAL,ALLOCATABLE :: u(:,:), v(:,:), hp(:,:), ut(:,:), vt(:,:), hpt(:,:), &
       utemp(:,:), vtemp(:,:), hptemp(:,:)
  REAL(kind=plflt),ALLOCATABLE :: hx(:), hy(:)
  
  NAMELIST/model_config/nx, ny, nt, g, h, dt

  
CONTAINS

  SUBROUTINE init(filename)
    CHARACTER(len=*),INTENT(in) :: filename
    INTEGER :: i

    OPEN(10, file=filename, status='old')
    READ(10, nml=model_config)
    CLOSE(10)

    nxu = nx+1
    nyv = ny+1
    dx = 1./nx
    dy = 1./ny

    PRINT*,'nx, ny, nt, g, h, dt:',nx, ny, nt, g, h, dt
    PRINT*,'cfl:',SQRT(g*h)/(MIN(dx,dy)/dt)
    
    ALLOCATE(hp(nx,ny), u(nxu,ny), v(nx,nyv))
    ALLOCATE(hpt(nx,ny), ut(nxu,ny), vt(nx,nyv))
    ALLOCATE(hptemp(nx,ny), utemp(nxu,ny), vtemp(nx,nyv))
    ALLOCATE(hx(nx), hy(ny))
    
    DO i = 1, nx
       hx(i) = dx*(REAL(i)-0.5)
    ENDDO
    DO i = 1, ny
       hy(i) = dy*(REAL(i)-0.5)
    ENDDO
    
  END SUBROUTINE init


  SUBROUTINE rhs(u, v, hp, ut, vt, hpt)
    REAL,INTENT(in) :: u(:,:), v(:,:), hp(:,:) ! stato del modello
    REAL,INTENT(out) :: ut(:,:), vt(:,:), hpt(:,:) ! tendenze calcolate

  INTEGER :: i,j

    DO j = 1, ny
       DO i = 2, nx-1
          ut(i,j) = -g*(hp(i,j)-hp(i-1,j))/dx
       ENDDO
    ENDDO

    DO j = 2, ny-1
       DO i = 1, nx
          vt(i,j) = -g*(hp(i,j)-hp(i,j-1))/dy
       ENDDO
    ENDDO

    DO j = 1, ny
       DO i = 1, nx
          hpt(i,j) = -h*((u(i+1,j)-u(i,j))/dx + (v(i,j+1)-v(i,j))/dy)
       ENDDO
    ENDDO
  END SUBROUTINE rhs

  
  SUBROUTINE rhs_4_ord(u, v, hp, ut, vt, hpt)
  REAL,INTENT(in) :: u(:,:), v(:,:), hp(:,:) ! stato del modello
  REAL,INTENT(out) :: ut(:,:), vt(:,:), hpt(:,:) ! tendenze calcolate

  INTEGER :: i,j

  DO j = 1, ny
    ut(2,j) = -g*(hp(2,j)-hp(1,j))/dx
    DO i = 3, nx-2
      ut(i,j) = -g*(-hp(i+1,j)+8.*hp(i,j)-8.*hp(i-1,j)+hp(i-2,j))/(12.*dx)
    ENDDO
    ut(nx-1,j) = -g*(hp(nx-1,j)-hp(nx-2,j))/dx
  ENDDO

  vt(1:nx,2) = -g*(hp(1:nx,2)-hp(1:nx,1))/dy
  DO j = 3, ny-2
    DO i = 1, nx
      vt(i,j) = -g*(-hp(i,j+1)+8.*hp(i,j)-8.*hp(i,j-1)+hp(i,j-2))/(12.*dy)
    ENDDO
  ENDDO
  vt(1:nx,ny-1) = -g*(hp(1:nx,ny-1)-hp(1:nx,ny-2))/dy

  DO j = 1, ny, ny
    DO i = 1, nx, nx
      hpt(i,j) = -h*((u(i+1,j)-u(i,j))/dx + (v(i,j+1)-v(i,j))/dy)
    ENDDO
  ENDDO

  DO j = 2, ny-1
    DO i = 2, nx-1
      hpt(i,j) = -h*( &
       (-u(i+2,j)+8.*u(i+1,j)-8.*u(i,j)+u(i-1,j))/(12.*dx) + &
       (-v(i,j+2)+8.*v(i,j+1)-8.*v(i,j)+v(i,j-1))/(12.*dy))
    ENDDO
  ENDDO

  END SUBROUTINE rhs_4_ord


  SUBROUTINE eulero

    CALL rhs(u, v, hp, ut, vt, hpt)
    u(2:nx-1,:) = u(2:nx-1,:) + dt*ut(2:nx-1,:)
    v(:,2:ny-1) = v(:,2:ny-1) + dt*vt(:,2:ny-1)
    hp(:,:) = hp(:,:) + dt*hpt(:,:)

  END SUBROUTINE eulero


  SUBROUTINE eulerogen(u, v, hp, ut, vt, hpt, utemp, vtemp, hptemp, dtg)
    REAL,INTENT(in) :: u(:,:), v(:,:), hp(:,:) ! stato iniziale del modello
    REAL,INTENT(in) :: ut(:,:), vt(:,:), hpt(:,:) ! tendenze calcolate
    REAL,INTENT(out) :: utemp(:,:), vtemp(:,:), hptemp(:,:) ! stato nuovo del modello
    REAL,INTENT(in) :: dtg

    utemp(2:nx-1,:) = u(2:nx-1,:) + dtg*ut(2:nx-1,:)
    vtemp(:,2:ny-1) = v(:,2:ny-1) + dtg*vt(:,2:ny-1)
    hptemp(:,:) = hp(:,:) + dtg*hpt(:,:)

  END SUBROUTINE eulerogen


  SUBROUTINE runge_kutta_2

    CALL rhs_4_ord(u, v, hp, ut, vt, hpt)
    call eulerogen(u, v, hp, ut, vt, hpt, utemp, vtemp, hptemp, dt/2.)
    CALL rhs_4_ord(utemp, vtemp, hptemp, ut, vt, hpt)
    ! piu' speditiva ma u, v, hp assegnate a diverse variabili e'
    ! potenziale fonte di errori
    !    call eulerogen(u, v, hp, ut, vt, hpt, u, v, hp, dt)
    ! piu' corretto: metto il risultato in *temp e aggiorno
    call eulerogen(u, v, hp, ut, vt, hpt, utemp, vtemp, hptemp, dt)
    u(2:nx-1,:) = utemp(2:nx-1,:)
    v(:,2:ny-1) = vtemp(:,2:ny-1)
    hp(:,:) = hptemp(:,:)

  END SUBROUTINE runge_kutta_2


  SUBROUTINE boundary(ub, vb)
    REAL,INTENT(in) :: ub, vb

    u(1,:) = ub
    u(nx+1,:) = ub
    v(:,1) = vb
    v(:,ny+1) = vb

  END SUBROUTINE boundary

  
  SUBROUTINE plot_contour(title, levels)
    USE plplot
    character(len=*) :: title
    REAL(plflt) :: levels(:)
    
    CALL pl_setcontlabelparam(0.006_plflt, 0.6_plflt, 0.1_plflt, 1)
    CALL pl_setcontlabelformat(6,4)
    CALL plcol0(1)
    CALL plenv(0._plflt, 1._plflt, 0._plflt, 1._plflt, 0, 0)
    CALL pllab('x', 'y', TRIM(title))
    CALL plcol0(1)
    CALL plcont(REAL(hp, kind=plflt), 1, nx, 1, ny, levels, hx, hy)

  END SUBROUTINE plot_contour

    
END MODULE shallow_water_mo
