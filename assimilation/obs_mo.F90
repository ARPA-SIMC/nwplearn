MODULE obs_mo
USE kinds
USE missing_values
USE coord_mo
USE func_mo
IMPLICIT NONE

!> Class for defining a single observation with associated
!! 2-d coordinates and error distribution variance.
TYPE,PUBLIC :: obs_t
  TYPE(coord_t) :: coord=coord_miss
  REAL(fp_d) :: val=dmiss
  REAL(fp_d) :: err=dmiss
END TYPE obs_t


!> Class for defining a set of observations with associated
!! default error distribution variance.
TYPE,PUBLIC :: obs_set_t
  REAL(fp_d) :: err=dmiss
  TYPE(obs_t),ALLOCATABLE :: obs(:)
  CONTAINS
  PROCEDURE :: export
  PROCEDURE :: import
  PROCEDURE :: compute_distance
  PROCEDURE :: perturb_gauss
END TYPE obs_set_t

PRIVATE

CONTAINS


!> Export an observation set to a namelist file.
SUBROUTINE export(this, filename)
CLASS(obs_set_t),INTENT(in) :: this
CHARACTER(len=*),INTENT(in) :: filename

INTEGER :: nobs
REAL(fp_d) :: default_err, x, y, val, err
NAMELIST/obs_set_nml/nobs,default_err
NAMELIST/obs_nml/x,y,val,err
INTEGER :: i

nobs = 0
IF (ALLOCATED(this%obs)) nobs = SIZE(this%obs)
default_err = this%err

OPEN(200,file=filename)
WRITE(200,NML=obs_set_nml)

DO i = 1, nobs
  x = this%obs(i)%coord%x
  y = this%obs(i)%coord%y
  val = this%obs(i)%val
  err = this%obs(i)%err
  WRITE(200,NML=obs_nml)
ENDDO

CLOSE(200)

END SUBROUTINE export


!> Import an observation set from a namelist file.
SUBROUTINE import(this, filename)
CLASS(obs_set_t),INTENT(out) :: this
CHARACTER(len=*),INTENT(in) :: filename

INTEGER :: nobs
REAL(fp_d) :: default_err, x, y, val, err
NAMELIST/obs_set_nml/nobs,default_err
NAMELIST/obs_nml/x,y,val,err
INTEGER :: i

nobs = 0
default_err = dmiss

OPEN(200,file=filename,status='OLD')
READ(200,NML=obs_set_nml)
this%err = default_err
ALLOCATE(this%obs(nobs))

DO i = 1, nobs
  err = dmiss
  READ(200,NML=obs_nml,END=100)
  IF (err == dmiss) err = default_err
  this%obs(i) = obs_t(coord_t(x, y), val, err)
ENDDO

100 CONTINUE
CLOSE(20)

END SUBROUTINE import


!> Compute a distance matrix for the observation set.
!! The result is allocated on the fly to a size of nobs x nobs.
FUNCTION compute_distance(this) RESULT(distance_mat)
CLASS(obs_set_t),INTENT(in) :: this
REAL(fp_d),ALLOCATABLE :: distance_mat(:,:)

INTEGER :: i,j,n

IF (.NOT.ALLOCATED(this%obs)) RETURN
n = SIZE(this%obs)
ALLOCATE(distance_mat(n,n))
distance_mat(1,1) = 0.0_fp_d
DO j = 2, n
  distance_mat(j,j) = 0.0_fp_d
  DO i = 1, j-1
    distance_mat(i,j) = this%obs(j)%coord%distance(this%obs(i)%coord)
    distance_mat(j,i) = distance_mat(i,j)
  ENDDO
ENDDO

END FUNCTION compute_distance


!> Perturb all the observations in the observation set
!! with a random Gaussian perturbation having zero average and
!! variance equal to the variance associated with each single
!! observation.
SUBROUTINE perturb_gauss(this)
CLASS(obs_set_t),INTENT(inout) :: this

REAL(fp_d) :: gauss_pert
REAL(fp_d) :: err
INTEGER :: i

IF (.NOT.ALLOCATED(this%obs)) RETURN
DO i = 1, SIZE(this%obs)
  CALL random_gaussian(gauss_pert)
  IF (this%obs(i)%err /= dmiss) THEN
    err = this%obs(i)%err
  ELSE
    err = this%err
  ENDIF
  this%obs(i)%val = this%obs(i)%val + gauss_pert*err
ENDDO

END SUBROUTINE perturb_gauss

END MODULE obs_mo
