      !=================================================================
      ! Linearly interpolate to obtain a value for a dependent variable.
      !
      ! Input:
      !    nz     Number of levels in grid
      !    zg     Gridded coordinate values
      !    vg     Gridded ordinate   values
      !    z      Interpolate to this location
      !
      ! Output:
      !    v      Interpolated value of the dependent variable
      !=================================================================
      real function linterp (nz, zg, vg, z) result (v)
         implicit none

         integer,                intent(in) :: nz
         real   , dimension(nz), intent(in) :: zg
         real   , dimension(nz), intent(in) :: vg
         real   ,                intent(in) :: z

         real    :: slope
         integer :: k

         !------------------------------------
         ! Find interval in which z is located
         !------------------------------------
         k = 1
10       k = k + 1
         if (zg(k) < z) go to 10

         !---------------------------------
         ! Compute slope over grid interval
         !---------------------------------
         slope = (vg(k)-vg(k-1)) / (zg(k)-zg(k-1))
         !---------------------------------
         ! Extrapolate upward to find value
         !---------------------------------
         v = vg(k-1) + slope * (z - zg(k-1))
         return
      end function linterp
