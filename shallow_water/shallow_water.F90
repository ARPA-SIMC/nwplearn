PROGRAM shallow_water
  USE shallow_water_mo
  USE plplot
  IMPLICIT NONE
  integer :: i, j
  REAL(kind=plflt) :: levels(21)


  CALL plsdev('pdfc')
  CALL plsfnam('shallow_water.pdf')
  CALL plscol0 (0, 255, 255, 255)
  CALL plscol0 (1, 0, 0, 0)
  CALL plinit()
  
  CALL init('sw.naml')

  u = 0.
  v = 0.
  hp = 0.
  ! condizione iniziale
  hp(NINT(0.25*nx):NINT(0.35*nx),NINT(0.25*ny):NINT(0.35*ny)) = 1.
  CALL boundary(0., 0.)
  ! imposto i livelli per il contour
  DO j = 1, size(levels)
     levels(j) = -1. + (j-1)*0.1_plflt
  ENDDO

  ! loop sul tempo
  DO i = 1, nt
     PRINT*,i,SUM(u)/((nx+1)*ny),SUM(v)/(nx*(ny+1)),SUM(hp)/(nx*ny)
     PRINT*,MINVAL(hp),MAXVAL(hp)
     CALL runge_kutta_2 ! integro nel tempo
     IF (mod(i,10) == 0) THEN
        call plot_contour('Runge Kutta 2', levels) ! disegno il grafico
     ENDIF
  ENDDO

  CALL plend()

END PROGRAM shallow_water
