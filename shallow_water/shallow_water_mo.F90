MODULE shallow_water_mo
  IMPLICIT NONE

  INTEGER :: nx, ny, nt
  INTEGER :: nxu, nyv
  REAL :: g=9.8, h=1., dt, dx, dy
  REAL,ALLOCATABLE :: u(:,:), v(:,:), hp(:,:), ut(:,:), vt(:,:), hpt(:,:), &
       utemp(:,:), vtemp(:,:), hptemp(:,:)

  NAMELIST/model_config/nx, ny, nt, g, h, dt

  
CONTAINS

  SUBROUTINE init(filename)
    CHARACTER(len=*),INTENT(in) :: filename

    OPEN(10, file=filename, status='old')
    READ(10, nml=model_config)
    CLOSE(10)
    PRINT*,'nx, ny, nt, g, h, dt:',nx, ny, nt, g, h, dt

    nxu = nx+1
    nyv = ny+1
    dx = 1./nx
    dy = 1./ny
    ALLOCATE(hp(nx,ny), u(nxu,ny), v(nx,nyv))
    ALLOCATE(hpt(nx,ny), ut(nxu,ny), vt(nx,nyv))
    ALLOCATE(hptemp(nx,ny), utemp(nxu,ny), vtemp(nx,nyv))
    
  END SUBROUTINE init


  SUBROUTINE rhs
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

  
  SUBROUTINE eulero

    CALL rhs
    u(2:nx-1,:) = u(2:nx-1,:) + dt*ut(2:nx-1,:)
    v(:,2:ny-1) = v(:,2:ny-1) + dt*vt(:,2:ny-1)
    hp(:,:) = hp(:,:) + dt*hpt(:,:)

  END SUBROUTINE eulero


  SUBROUTINE eulerogen(dtg)
    REAL,INTENT(in) :: dtg

    CALL rhs
    utemp(2:nx-1,:) = u(2:nx-1,:) + dtg*ut(2:nx-1,:)
    vtemp(:,2:ny-1) = v(:,2:ny-1) + dtg*vt(:,2:ny-1)
    hptemp(:,:) = hp(:,:) + dtg*hpt(:,:)

  END SUBROUTINE eulero


  SUBROUTINE runge_kutta_2

    utemp = u
    
    call eulerogen(dt/2.)

    CALL rhs
    utemp(2:nx-1,:) = u(2:nx-1,:) + dt*ut(2:nx-1,:)
    vtemp(:,2:ny-1) = v(:,2:ny-1) + dt*vt(:,2:ny-1)
    hptemp(:,:) = hp(:,:) + dt*hpt(:,:)

  END SUBROUTINE eulero


  SUBROUTINE boundary(ub, vb)
    REAL,INTENT(in) :: ub, vb

    u(1,:) = ub
    u(nx+1,:) = ub
    v(:,1) = vb
    v(:,ny+1) = vb

  END SUBROUTINE boundary

END MODULE shallow_water_mo
