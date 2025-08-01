c     This function returns the temperature as a function
c     of height. Uses linear interpolation.

c     Will work for altitudes outside of the data range, uses
c     linear extrapolation

c     NOTE: This here is the Dunion "moist-tropical" sounding up to 50 mb

      real function gettofz(z) result(t)
         implicit none

         integer, parameter :: levels = 16
         real   , parameter :: tZero = 273.15

         real, dimension(levels), parameter :: zData =
     &      (/
     &           10.,    125.,   810.,  1541.,
     &          3178.,  4438.,  5887.,  7596.,
     &          9690., 10949., 12417., 14202.,
     &         16589., 20727., 22139., 23971.
     &      /)

         real, dimension(levels), parameter :: tCent =
     &      (/
     &          26.8,  26.5,  21.9,  17.6,
     &           8.9,   1.5,  -6.6, -17.1,
     &         -32.3, -42.4, -54.4, -67.2,
     &         -74.4, -63.0, -57.3, -54.0
     &      /)

         real, intent(in) :: z

         integer :: i

         real :: dz

         i = 1

10       i = i + 1

         if (zData(i) < z) go to 10

         dz = zData(i) - zData(i-1)

        t = (1. - (zdata(i) -  z) /dz)*tCent(i) +
     &      (1. - (z - zdata(i-1))/dz)*tCent(i-1) + tZero

        return
      end
