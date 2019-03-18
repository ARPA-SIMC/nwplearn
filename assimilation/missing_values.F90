MODULE missing_values
USE kinds
IMPLICIT NONE

REAL, PARAMETER :: rmiss = HUGE(1.0) !< default single precision real
DOUBLE PRECISION, PARAMETER :: dmiss = HUGE(1.0D0) !< default double precision real
REAL(kind=fp_s), PARAMETER :: rsmiss = HUGE(1.0_fp_s) !< single precision IEEE real \a (kind=fp_s)
REAL(kind=fp_d), PARAMETER :: rdmiss = HUGE(1.0_fp_d) !< double precision IEEE real \a (kind=fp_d)
INTEGER, PARAMETER :: imiss = HUGE(0) !< default integer
INTEGER(kind=int_b), PARAMETER :: ibmiss = HUGE(0_int_b) !< 1-byte integer \a (kind=int_b)
INTEGER(kind=int_b), PARAMETER :: bmiss = ibmiss
INTEGER(kind=int_s), PARAMETER :: ismiss = HUGE(0_int_s) !< 2-byte integer \a (kind=int_s)
INTEGER(kind=int_l), PARAMETER :: ilmiss = HUGE(0_int_l) !< 4-byte integer \a (kind=int_l)
INTEGER(kind=int_ll), PARAMETER :: illmiss = HUGE(0_int_ll) !< 8-byte integer if supported \a (kind=int_ll)
CHARACTER(len=1), PARAMETER :: cmiss = char(0) !< character (any length)

END MODULE missing_values
