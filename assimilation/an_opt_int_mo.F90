MODULE an_opt_int_mo
USE kinds
USE obs_mo
USE coord_mo
USE grid_mo
USE func_mo
IMPLICIT NONE


!> Class for performing an objective analysis using the optimal interpolation
!! technique.
!! The object must be filled with information about how to compute the
!! observation and background error covariance matrices.
TYPE,PUBLIC :: an_opt_int_t
  CHARACTER(len=20) :: o_mat_func = cmiss !< function for computing observation error covariance matrix O or R through \a func_mo module
  CHARACTER(len=20) :: b_mat_func = cmiss !< function for computing background error covariance matrix B through \a func_mo module
  REAL(fp_d) :: o_mat_fact=1.0_fp_d !< distance rescaling factor for argument of O matrix function
  REAL(fp_d) :: b_mat_fact=1.0_fp_d !< distance rescaling factor for argument of B matrix function
  CONTAINS
  PROCEDURE :: compute
END TYPE an_opt_int_t

PRIVATE

CONTAINS


!> Compute the optimal interpolation analysis given the covariance matrix
!! definitions stored in \a this, the observation set \a obs_set and the
!! background on grid \a bg.
!! The result is the analysed value on the same grid as the background
!! and it is stored in the \a an array allocated within the soubroutine.
!! If any error occurs, \a an will no be allocated.
SUBROUTINE compute(this, bg, obs_set, an)
CLASS(an_opt_int_t),INTENT(in) :: this
TYPE(grid_t),INTENT(inout) :: bg
TYPE(obs_set_t),INTENT(inout) :: obs_set
REAL(fp_d),INTENT(out),ALLOCATABLE :: an(:,:)

REAL(fp_d),ALLOCATABLE :: distoo(:,:), distob(:), bg_on_obs(:), &
 om(:,:), bm(:,:), brhs(:), a(:,:), work(:), ev(:)
!REAL(fp_d) :: distfact, freqx, freqy
INTEGER,ALLOCATABLE :: ipiv(:)
INTEGER :: i, j, n, ierr, lwork


IF (.NOT.ALLOCATED(obs_set%obs)) RETURN
n = SIZE(obs_set%obs)
! observation distance matrix
distoo = obs_set%compute_distance()

! observation error covariance matrix
om = obs_set%err*obs_set%err*func_autocorr(distoo/this%o_mat_fact, this%o_mat_func)
! background error covariance matrix
bm = bg%err*bg%err*func_autocorr(distoo/this%b_mat_fact, this%b_mat_func)
! cut small bm values
!WHERE (bm < obs_set%err*obs_set%err*0.01_fp_d)
!  bm = 0.0_fp_d
!END WHERE

! interpolate background on observation points
bg_on_obs = bg%interpol(obs_set%obs(:)%coord)
! allocate result
ALLOCATE(an(bg%nx,bg%ny))
! allocate lapack workspace
ALLOCATE(ipiv(n),work(4*n),ev(n))
lwork = SIZE(work)

a = om
CALL dsyev('N', 'U', n, a, n, ev, work, lwork, ierr)
!PRINT*,'eigenvalues',ierr,n,COUNT(ev > 0.0_fp_d),COUNT(ev < 0.0_fp_d)
!PRINT*,'eigenvalues',ev
DO j = 1, bg%ny
  DO i = 1, bg%nx
! a is overwritten at every call
    a = om + bm
! column vector with gridpoint-observations distance
    distob = distance(bg%coord(i,j), obs_set%obs(:)%coord)
! right hand side vector with gridpoint-observations background covariance
    brhs = bg%err*bg%err*func_autocorr(distob/this%b_mat_fact, this%b_mat_func)
! call lapack for inverting matrix
    CALL dposv('U', n, 1, a, n, brhs, n, ierr)
!    CALL dsysv('U', n, 1, a, n, ipiv, brhs, n, work, lwork, ierr)
    
! error condition
    IF (ierr /= 0) THEN
      PRINT*,'Error in lapack',ierr
      DEALLOCATE(an)
      RETURN
    ENDIF
! set analysis output
    an(i,j) = bg%val(i,j) + DOT_PRODUCT(brhs, obs_set%obs(:)%val-bg_on_obs(:))
  ENDDO
ENDDO

END SUBROUTINE compute


END MODULE an_opt_int_mo
