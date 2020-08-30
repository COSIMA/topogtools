test_one.f90:484.6:

      rad2=maxval((cx-t_source(1))**2+(cy-t_source(2))**2+(cz-t_source(3))**2))
      1
Error: Unclassifiable statement at (1)
test_one.f90:452.14:

   do i=1,i_in
              1
Error: Symbol 'i_in' at (1) has no IMPLICIT type
test_one.f90:451.11:

do j=1,j_in
           1
Error: Symbol 'j_in' at (1) has no IMPLICIT type
test_one.f90:480.6:

      cx=cos(xc)*cos(yc)
      1
Error: Incompatible ranks 0 and 2 in assignment at (1)
test_one.f90:481.6:

      cy=sin(xc)*cos(yc)
      1
Error: Incompatible ranks 0 and 2 in assignment at (1)
test_one.f90:482.6:

      cz=sin(yc)
      1
Error: Incompatible ranks 0 and 2 in assignment at (1)
