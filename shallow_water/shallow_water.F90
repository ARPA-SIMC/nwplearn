PROGRAM shallow_water
  USE shallow_water_mo
  IMPLICIT NONE
  integer :: i
  
  CALL init('sw.naml')

  u = 1.
  v = 1.
  hp = 0.
  CALL boundary(0., 0.)
  
  ! loop sul tempo
  DO i = 1, nt
     PRINT*,i,SUM(u)/((nx+1)*ny),SUM(v)/(nx*(ny+1)),SUM(hp)/(nx*ny)
     CALL eulero
  ENDDO

END PROGRAM shallow_water
