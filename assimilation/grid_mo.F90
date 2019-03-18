MODULE grid_mo
USE kinds
USE missing_values
USE coord_mo
IMPLICIT NONE

!> Class for defining a rectangular 2-d grid with points
!! equally spaced in each direction. The particular cases
!! of degenerated 1-d grids along x or y direction are supported.
TYPE,PUBLIC :: grid_t
  REAL(fp_d) :: err=dmiss !< error associated to the grid values
  INTEGER :: nx=0 !< number of points along x axis
  INTEGER :: ny=0 !< number of points along y axis
  TYPE(coord_t) :: llc=coord_miss !< coordinates of lower left corner
  TYPE(coord_t) :: urc=coord_miss !< coordinates of upper right corner
  TYPE(coord_t),ALLOCATABLE :: coord(:,:) !< coordinates of all grid points
  REAL(fp_d),ALLOCATABLE :: val(:,:) !< values at the grid points
  CONTAINS
  PROCEDURE :: alloc
  PROCEDURE :: compute_coord
  PROCEDURE :: interpol
  PROCEDURE :: getx
  PROCEDURE :: gety
  PROCEDURE :: getx2
  PROCEDURE :: gety2
  PROCEDURE :: export
  PROCEDURE :: import
  PROCEDURE :: plot_contour
END TYPE grid_t

!INTERFACE grid_t
!  MODULE PROCEDURE grid_t_new
!END INTERFACE grid_t

PRIVATE

CONTAINS


!> This subroutine allocates this%coord and stores in it
!! the coordinates of each point of the grid.
SUBROUTINE compute_coord(this)
CLASS(grid_t),INTENT(inout) :: this

INTEGER :: i

IF (ALLOCATED(this%coord)) RETURN

IF (this%nx > 0 .AND. this%ny > 0) THEN
  ALLOCATE(this%coord(this%nx,this%ny))

  this%coord(:,1)%x = get_coord(this%nx, this%llc%x, this%urc%x)
  DO i = 2, this%ny
    this%coord(:,i)%x = this%coord(:,1)%x
  ENDDO

  this%coord(1,:)%y = get_coord(this%ny, this%llc%y, this%urc%y)
  DO i = 2, this%nx
    this%coord(i,:)%y = this%coord(1,:)%y
  ENDDO

ENDIF

END SUBROUTINE compute_coord


!> This subroutine allocates the data array and the coordinate array
!! and computes the coordinates.
!! The members \a nx, \a ny, \a llc, \a urc must have
!! already been initialized.
!! The \a compute_coord method is used internally.
SUBROUTINE alloc(this)
CLASS(grid_t),INTENT(inout) :: this

CALL this%compute_coord()
IF (this%nx > 0 .AND. this%ny > 0) THEN
  ALLOCATE(this%val(this%nx,this%ny))
  this%val = dmiss
ENDIF

END SUBROUTINE alloc


!> Interpolate the grid value on a list of requested target points.
!! The target points are specified as an array of \a coord_t objects.
!! Bilinear interpolation is used for 2d grids and linear for 1d grids.
!! The result is an allocatable array, allocated to the size of \a coord.
FUNCTION interpol(this, coord)
CLASS(grid_t),INTENT(inout) :: this
TYPE(coord_t),INTENT(in) :: coord(:) !< target points
REAL(fp_d),ALLOCATABLE :: interpol(:)

INTEGER :: i, indx, indy

CALL this%compute_coord()

ALLOCATE(interpol(SIZE(coord)))

IF (this%nx > 1 .AND. this%ny > 1) THEN ! 2d

  DO i = 1, SIZE(coord)
    indx = FLOOR((coord(i)%x-this%llc%x)/ &
     ((this%urc%x-this%llc%x)/REAL(this%nx-1,kind=fp_d))) + 1
    indy = FLOOR((coord(i)%y-this%llc%y)/ &
     ((this%urc%y-this%llc%y)/REAL(this%ny-1,kind=fp_d))) + 1
    indx = MIN(MAX(1,indx),this%nx-1)
    indy = MIN(MAX(1,indy),this%ny-1)

    interpol(i) = hbilin(this%val(indx,indy), this%val(indx+1,indy), &
     this%val(indx+1,indy+1), this%val(indx,indy+1), &
     this%coord(indx,indy)%x, this%coord(indx,indy)%y, &
     this%coord(indx+1,indy+1)%x, this%coord(indx+1,indy+1)%y, &
     coord(i)%x, coord(i)%y)
  ENDDO

ELSE IF (this%nx > 1) THEN ! 1d

  DO i = 1, SIZE(coord)
    indx = FLOOR((coord(i)%x-this%llc%x)/ &
     ((this%urc%x-this%llc%x)/REAL(this%nx-1,kind=fp_d))) + 1
    indx = MIN(MAX(1,indx),this%nx-1)

    interpol(i) = hlin(this%val(indx,1), this%val(indx+1,1), &
     this%coord(indx,1)%x, this%coord(indx+1,1)%x, coord(i)%x)
  ENDDO

ELSE IF (this%ny > 1) THEN ! 1d

  DO i = 1, SIZE(coord)
    indy = FLOOR((coord(i)%y-this%llc%y)/ &
     ((this%urc%y-this%llc%y)/REAL(this%ny-1,kind=fp_d))) + 1
    indy = MIN(MAX(1,indy),this%ny-1)

    interpol(i) = hlin(this%val(1,indy), this%val(1,indy+1), &
     this%coord(1,indy)%y, this%coord(1,indy+1)%y, coord(i)%y)
  ENDDO

ENDIF

END FUNCTION interpol


!> Export grid definition, excluding grid values, to a namelist file.
SUBROUTINE export(this, filename)
CLASS(grid_t),INTENT(in) :: this
CHARACTER(len=*),INTENT(in) :: filename

REAL(fp_d) :: err
INTEGER :: nx, ny
REAL(fp_d) :: llcx,llcy,urcx,urcy
NAMELIST/grid_nml/err,nx,ny,llcx,llcy,urcx,urcy

err=this%err; nx=this%nx; ny=this%ny
llcx=this%llc%x; llcy=this%llc%y; urcx=this%urc%x; urcy=this%urc%y;

OPEN(200,file=filename)
WRITE(200,NML=grid_nml)
CLOSE(200)

END SUBROUTINE export


!> Import grid definition, excluding grid values, from a namelist file.
SUBROUTINE import(this, filename)
CLASS(grid_t),INTENT(out) :: this
CHARACTER(len=*),INTENT(in) :: filename

REAL(fp_d) :: err
INTEGER :: nx, ny
REAL(fp_d) :: llcx,llcy,urcx,urcy
NAMELIST/grid_nml/err,nx,ny,llcx,llcy,urcx,urcy

err=dmiss; nx=0; ny=0; llcx=dmiss; llcy=dmiss; urcx=dmiss; urcy=dmiss;

OPEN(200,file=filename,status='OLD')
READ(200,NML=grid_nml)
this%err=err; this%nx=nx; this%ny=ny
this%llc%x=llcx; this%llc%y=llcy; this%urc%x=urcx; this%urc%y=urcy;
! will work later
!this = grid_t(err, nx, ny, coord_t(llcx, llcy), cord_t(urcx, urcy))
CLOSE(200)

END SUBROUTINE import


FUNCTION getx(this) RESULT(x)
CLASS(grid_t),INTENT(in) :: this

REAL(fp_d),ALLOCATABLE :: x(:)

x(:) = get_coord(this%nx, this%llc%x, this%urc%x)

END FUNCTION getx


FUNCTION getx2(this) RESULT(x)
CLASS(grid_t),INTENT(in) :: this
REAL(fp_d),ALLOCATABLE :: x(:,:)

INTEGER :: i

ALLOCATE(x(this%nx,this%ny))
x(:,1) = get_coord(this%nx, this%llc%x, this%urc%x)
DO i = 2, this%ny
  x(:,i) = x(:,1)
ENDDO

END FUNCTION getx2


FUNCTION gety(this) RESULT(y)
CLASS(grid_t),INTENT(in) :: this
REAL(fp_d),ALLOCATABLE :: y(:)

y(:) = get_coord(this%ny, this%llc%y, this%urc%y)

END FUNCTION gety


FUNCTION gety2(this) RESULT(y)
CLASS(grid_t),INTENT(in) :: this
REAL(fp_d),ALLOCATABLE :: y(:,:)

INTEGER :: i

ALLOCATE(y(this%nx,this%ny))
y(1,:) = get_coord(this%ny, this%llc%y, this%urc%y)
DO i = 2, this%nx
  y(i,:) = y(1,:)
ENDDO

END FUNCTION gety2


! internal function for computing the coordinates along an axis
FUNCTION get_coord(n, xmin, xmax) RESULT(coord)
INTEGER,INTENT(in) :: n
REAL(fp_d),INTENT(in) :: xmin, xmax
REAL(fp_d) :: coord(n)

INTEGER :: i
REAL(fp_d) :: dcoord, frac

dcoord = xmax - xmin
coord(1) = xmin
coord(n) = xmax
DO i = 2, n-1
  frac = REAL(i-1, kind=fp_d)/REAL(n-1, kind=fp_d)
  coord(i) = xmin + frac*dcoord
ENDDO

END FUNCTION get_coord


ELEMENTAL FUNCTION hbilin (z1,z2,z3,z4,x1,y1,x3,y3,xp,yp) RESULT(zp)
REAL(fp_d),INTENT(in) :: z1,z2,z3,z4 ! Z values on the four points
REAL(fp_d),INTENT(in):: x1,y1 ! coordinate of the lower left point
REAL(fp_d),INTENT(in):: x3,y3 ! coordinate of the upper right point
REAL(fp_d),INTENT(in):: xp,yp ! coordinate of point where interpolate

REAL(fp_d) :: zp,p1,p2,z5,z6

p2 = (yp-y1)/(y3-y1)
p1 = (xp-x1)/(x3-x1)

z5 = (z4-z1)*p2+z1
z6 = (z3-z2)*p2+z2

zp = (z6-z5)*p1+z5
      
END FUNCTION hbilin


ELEMENTAL FUNCTION hlin (z1,z2,x1,x3,xp) RESULT(zp)
REAL(fp_d),INTENT(in) :: z1,z2 ! Z values on the two points
REAL(fp_d),INTENT(in):: x1 ! coordinate of the lower left point
REAL(fp_d),INTENT(in):: x3 ! coordinate of the upper right point
REAL(fp_d),INTENT(in):: xp ! coordinate of point where interpolate

REAL(fp_d) :: zp,p1

p1 = (xp-x1)/(x3-x1)

zp = (z2-z1)*p1+z1
      
END FUNCTION hlin


SUBROUTINE plot_contour(this, title, levels)
USE plplot
CLASS(grid_t),INTENT(in) :: this
CHARACTER(len=*),INTENT(in) :: title
REAL(kind=plflt),INTENT(in) :: levels(:)

REAL(kind=plflt) :: tr(6)

IF (.NOT.ALLOCATED(this%val)) RETURN
! Identity transform for plcont
tr(1) = 2._plflt/DBLE(SIZE(this%val,1)-1)
tr(2) = 0.0_plflt
tr(3) = -1.0_plflt
tr(4) = 0.0_plflt
tr(5) = 2._plflt/DBLE(SIZE(this%val,2)-1)
tr(6) = -1.0_plflt

CALL pl_setcontlabelparam(0.006_plflt, 0.6_plflt, 0.1_plflt, 1)
CALL pl_setcontlabelformat(6,4)
CALL plcol0(1)
CALL plenv(this%llc%x, this%urc%x, this%llc%y, this%urc%y, 0, 0)
CALL pllab('x', 'y', TRIM(title))
CALL plcol0(1)
!CALL plcont(this%val,1,SIZE(this%val,1),1,SIZE(this%val,2),levels,tr)
CALL plcont(this%val,1,SIZE(this%val,1),1,SIZE(this%val,2),levels,&
 this%coord%x,this%coord%y)

END SUBROUTINE plot_contour

END MODULE grid_mo
