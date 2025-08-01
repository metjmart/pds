!     This version of makevtbar uses a modifed rankine vortex
!     to define the radial structure of the initial wind field

!     The result comes back in the vector vtbar
!     the r locations are assumed to be at 0,1*dr,2*dr,....(nr-1)*dr

      subroutine makevtbar(vtbar, nr, dr)

         implicit none

              real, dimension (nr), intent(out) :: vtbar

              real    :: r
              integer :: i,j

!              Local variables to define vortex:
              real :: Rmw
              real :: Vmax
              integer, intent(in) :: nr
              real, intent(in) :: dr
              real, parameter :: a = .4

              integer, parameter :: itMax = 20    ! Number of smoothing
                                                  ! interations

              Rmw = 90000.0
              Vmax = 20.0 * 1.09     ! Adjust upward because of smoothing
!             Vmax = 0.0 * 1.09     ! for TVPDS env sounding

         print *, 'In makevtbar.f, Vmax = ' , Vmax
         print *, 'In makevtbar.f, Rmw = ' , Rmw
         print *, 'In makevtbar.f, nr = ' , nr
         print *, 'In makevtbar.f, dr = ' , dr

         do j = 1, nr

            r = (j-1) * dr

            if (r <= Rmw) then
               vtbar(j) = Vmax * (r/Rmw)      ! Linear increase to Vmax
            else
               vtbar(j) = Vmax * (Rmw/r)**a   ! (1/r)**a falloff
            end if

         end do

!        I dont like cusps. Lets do a little (1:2:1) smoothing:
         do i = 1, itMax
            do j = 2, nr-1
               vtbar(j) = .25*(vtbar(j-1) + 2.*vtbar(j) + vtbar(j+1))
            end do
         end do

!        print *,vtbar

         return
      end subroutine makevtbar


