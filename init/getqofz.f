c        This function returns the moisture as a function
c        of height.  Uses linear interpolation.

c        Will work for altitudes outside of the data range, uses
c        linear extrapolation

c       NOTE: This here the Dunion "moist-tropical" sounding

      real function getqofz(z) result(q)
         implicit none

         integer, parameter :: levels = 16

         real, dimension(levels), parameter :: zData =
     &      (/
     &           10.,    125.,   810.,  1541.,
     &          3178.,  4438.,  5887.,  7596.,
     &          9690., 10949., 12417., 14202.,
     &         16589., 20727., 22139., 23971.
     &      /)

         real, dimension(levels), parameter :: qData =
     &      (/
     &         0.01862, 0.01847, 0.01526, 0.01196,
     &         0.00673, 0.00412, 0.00241, 0.00112,
     &         0.00033, 0.00004, 0.00001, 0.00001,
     &         0.     , 0.     , 0.     , 0.
     &      /)

         real, intent(in) :: z

         integer :: i

         real :: dz

         i = 1

10       i = i + 1

         if (zData(i) < z) go to 10

         dz = zData(i) - zData(i-1)

         q = (1. - (zData(i) - z)/dz)*qData(i) +
     &       (1. - (z - zData(i-1))/dz)*qData(i-1)

         return
      end function getqofz
