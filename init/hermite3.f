c           This function returns the 3rd order Hermite polynomial
c           with h(0)=1.; h'(0)=h'(1)=h(1)=0;

           real function hermite3(x)

              implicit none

              real, intent(in) :: x

              hermite3 = 1. - 3.*x*x + 2.*x*x*x

              return

           end function hermite3
