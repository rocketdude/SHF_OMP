program eps
  !This program is to calculate the machine epsilon
  !Run the program and then inser the results to "eps" in main.f90
  integer*4 i 
  real*8 eps1,eps2 
  eps1 = 1 
  do while (.not. 1.+eps1 == 1.) 
     eps2 = eps1 
     eps1 = eps1/2 
  end do 
  print *,"The machine epsilon is",eps2 
  print *,"The 1 + machine epsilon is ",1+eps2 
end program eps 
