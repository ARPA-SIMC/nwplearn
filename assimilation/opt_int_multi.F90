PROGRAM do_opt_int_multi
USE plplot
IMPLICIT NONE

CALL plsdev('pdfc')
CALL plsfnam('opt_int_multi.pdf')
!CALL plsetopt('portrait', ' ')
!CALL plparseopts(PL_PARSE_FULL)
CALL plscol0 (0, 255, 255, 255)
CALL plscol0 (1, 0, 0, 0)
CALL plinit()

CALL opt_int_multi(5)
CALL opt_int_multi(10)
CALL opt_int_multi(20)
CALL opt_int_multi(50)
CALL opt_int_multi(100)
CALL opt_int_multi(200)
CALL opt_int_multi(500)
!CALL opt_int_multi(1000)

CALL plend()

END PROGRAM do_opt_int_multi

SUBROUTINE opt_int_multi(n)
USE grid_mo
USE obs_mo
USE an_opt_int_mo
USE func_mo
USE coord_mo
USE plplot
IMPLICIT NONE
INTEGER,INTENT(in) :: n !< number of observations

TYPE(obs_set_t) :: obs_set
TYPE(grid_t) :: bg
TYPE(an_opt_int_t) :: an_opt_int
REAL(fp_d),ALLOCATABLE :: an(:,:), truth(:,:)
REAL(fp_d) :: freqx, freqy
REAL(fp_d) :: clev(11)
INTEGER :: i

! levels for contouring
clev = (/( 9.0_plflt+0.2_plflt*(i-1),i=1,SIZE(clev))/)

CALL init_random_seed()
! "truth" function wavenumber in x and y
freqx = 4.0_fp_d
freqy = 4.0_fp_d
! observation error distribution variance
obs_set%err = 0.1_fp_d
! analytical functions to use for covariances
!an_opt_int%o_mat_func = 'delta'
an_opt_int%o_mat_func = 'gauss'
an_opt_int%o_mat_fact = 0.05_fp_d
an_opt_int%b_mat_func = 'gauss'
an_opt_int%b_mat_fact = 0.2_fp_d

CALL bg%import('grid.naml')
CALL bg%alloc()
bg%val(:,:) = 10.0_fp_d
CALL bg%compute_coord()

ALLOCATE(obs_set%obs(n))

! define random coordinates in unit square for observations
i = 1
DO WHILE(i <= n)
  CALL RANDOM_NUMBER(obs_set%obs(i)%coord%x)
  CALL RANDOM_NUMBER(obs_set%obs(i)%coord%y)
! ensure points are distant enough 
  IF (i > 1) THEN
    IF (ANY(distance(obs_set%obs(i)%coord, obs_set%obs(:i-1)%coord) < 0.02)) THEN
      CYCLE
    ENDIF
  ENDIF
  i = i + 1
ENDDO

! set observation to "truth" and perturb
obs_set%obs(:)%val = func_2dfield(freqx*obs_set%obs(:)%coord%x, &
 freqy*obs_set%obs(:)%coord%y, 'fourier') + 10.0_fp_d
CALL obs_set%perturb_gauss()

! compute analysis
CALL an_opt_int%compute(bg, obs_set, an)

IF (ALLOCATED(an)) THEN
! compute gridded "truth"
  truth = func_2dfield(freqx*bg%coord%x, freqy*bg%coord%y, 'fourier') + 10.0_fp_d
! compute error
  PRINT*,n,'Abs err:',SUM(ABS(an-truth))/SIZE(truth)

! plot results
  bg%val = an
  CALL bg%plot_contour('Analysis', clev)
  CALL plstring(obs_set%obs(:)%coord%x, obs_set%obs(:)%coord%y, 'x')
  bg%val = truth
  CALL bg%plot_contour('Truth', clev)
ENDIF


END SUBROUTINE opt_int_multi
