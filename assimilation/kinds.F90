MODULE kinds
IMPLICIT NONE

INTEGER, PARAMETER :: int_b    = SELECTED_INT_KIND(1) !< 1-byte integer (byte)
INTEGER, PARAMETER :: int_s    = SELECTED_INT_KIND(4) !< 2-byte integer (short)
INTEGER, PARAMETER :: int_l    = SELECTED_INT_KIND(8) !< 4-byte integer (long)
INTEGER, PARAMETER, PRIVATE :: &
 int_ll_t = SELECTED_INT_KIND(16)
!> 8-byte integer (long long) if supported, otherwise 4-byte integer
INTEGER, PARAMETER :: int_ll = &
 ( ( ( 1 + SIGN( 1, int_ll_t ) ) / 2 ) * int_ll_t ) + &
 ( ( ( 1 - SIGN( 1, int_ll_t ) ) / 2 ) * int_l    )

INTEGER, PARAMETER :: fp_s = SELECTED_REAL_KIND(6) !< single precision floating point (4 byte IEEE)
INTEGER, PARAMETER :: fp_d = SELECTED_REAL_KIND(15) !< double precision floating point (8 byte IEEE)
INTEGER, PARAMETER, PRIVATE :: fp_q_t = SELECTED_REAL_KIND(20)
!> quad precision floating point (16 byte IEEE) if supported, otherwise double precision floating point
INTEGER, PARAMETER :: fp_q = &
 ( ( ( 1 + SIGN( 1, fp_q_t ) ) / 2 ) * fp_q_t ) + &
 ( ( ( 1 - SIGN( 1, fp_q_t ) ) / 2 ) * fp_d )

!INTEGER, PARAMETER :: ptr_c = SIZEOF_PTR_C !< kind for an integer having the same size of a C pointer

END MODULE kinds
