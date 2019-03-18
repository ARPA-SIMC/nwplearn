MODULE coord_mo
USE kinds
USE missing_values
IMPLICIT NONE

TYPE,PUBLIC ::  coord_t
  REAL(fp_d) :: x=dmiss
  REAL(fp_d) :: y=dmiss
  CONTAINS
  PROCEDURE :: distance => coord_distance
!  PROCEDURE,NOPASS :: distance2 => coord_distance
END TYPE coord_t

TYPE(coord_t),PUBLIC,PARAMETER :: coord_miss=coord_t(dmiss,dmiss)

INTERFACE distance
  MODULE PROCEDURE coord_distance2
END INTERFACE

PUBLIC distance
PRIVATE

CONTAINS

ELEMENTAL FUNCTION coord_distance(this, other) RESULT(distance)
CLASS(coord_t),INTENT(in) :: this
CLASS(coord_t),INTENT(in) :: other
REAL(fp_d) :: distance

distance = SQRT((this%x-other%x)**2 + (this%y-other%y)**2)

END FUNCTION coord_distance


ELEMENTAL FUNCTION coord_distance2(this, other) RESULT(distance)
TYPE(coord_t),INTENT(in) :: this
TYPE(coord_t),INTENT(in) :: other
REAL(fp_d) :: distance

distance = SQRT((this%x-other%x)**2 + (this%y-other%y)**2)

END FUNCTION coord_distance2

END MODULE coord_mo
