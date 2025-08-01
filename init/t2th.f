      !==================================================
      ! Convert temperature to potential temperature.
      !
      ! Input:
      !    t     Temperature (degrees Kelvin)
      !    p     Pressure (Pascals)
      !
      ! Output:
      !    t2th  Potential temperature (degrees Kelvin)
      !==================================================
      real function t2th (t, p)
         implicit none
!        real, parameter  :: rcp = 287.04 / 1004.64
         real, parameter  :: rcp = 2. / 7.
         real, intent(in) :: t
         real, intent(in) :: p
         t2th = t * (100000. / p) ** rcp
      end function t2th
