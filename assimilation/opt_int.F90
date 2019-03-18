PROGRAM do_opt_int
USE plplot
IMPLICIT NONE
CHARACTER(len=512) :: obsfile

CALL plsdev('pdfc')
CALL plsfnam('opt_int.pdf')
CALL plsetopt('portrait', ' ')
!CALL plparseopts(PL_PARSE_FULL)
CALL plscol0 (0, 255, 255, 255)
CALL plscol0 (1, 0, 0, 0)
CALL plinit()

CALL getarg(1,obsfile)
IF (obsfile == '') obsfile='obs_input.naml'
PRINT*,'Reading observations from ',TRIM(obsfile)

CALL opt_int(obsfile)

CALL plend()

END PROGRAM do_opt_int

SUBROUTINE opt_int(obsfile)
USE grid_mo
USE obs_mo
USE an_opt_int_mo
USE func_mo
USE coord_mo
USE plplot
IMPLICIT NONE
CHARACTER(len=*) :: obsfile

TYPE(obs_set_t) :: obs_set
TYPE(grid_t) :: bg
TYPE(an_opt_int_t) :: an_opt_int
REAL(fp_d),ALLOCATABLE :: an(:,:)
REAL(fp_d) :: clev(11)
INTEGER :: i

! levels for contouring
clev = (/( 9.0_plflt+0.2_plflt*(i-1),i=1,SIZE(clev))/)

! analytical functions to use for covariances
!an_opt_int%o_mat_func = 'delta'
an_opt_int%o_mat_func = 'gauss'
an_opt_int%o_mat_fact = 0.1_fp_d
an_opt_int%b_mat_func = 'gauss'
an_opt_int%b_mat_fact = 0.4_fp_d

! import background/analysis grid definition
CALL bg%import('grid.naml')
! allocate grid
CALL bg%alloc()
! set bg value
bg%val(:,:) = 10.0_fp_d
! compute grid coordinates
CALL bg%compute_coord()

! import observations
CALL obs_set%import(obsfile)

! compute analysis
CALL an_opt_int%compute(bg, obs_set, an)

IF (ALLOCATED(an)) THEN
! plot results
  bg%val = an
  CALL bg%plot_contour('Analysis', clev)
  CALL plstring(obs_set%obs(:)%coord%x, obs_set%obs(:)%coord%y, 'x')
ENDIF


END SUBROUTINE opt_int
