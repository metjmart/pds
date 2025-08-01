      !==================================================
      ! Convert potential temperature to temperature.
      !
      ! Input:
      !    th    Potential temperature (degrees Kelvin)
      !    p     Pressure (Pascals)
      !
      ! Output:
      !    th2t  Temperature (degrees Kelvin)
      !==================================================
      real function th2t (th, p)
!        real, parameter  :: rcp = 287.04 / 1004.64
         real, parameter  :: rcp = 2. / 7.
         real, intent(in) :: th
         real, intent(in) :: p
         th2t = th * (.00001 * p) ** rcp
      end function th2t
