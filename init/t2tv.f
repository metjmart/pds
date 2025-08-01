      !==================================================
      ! Convert temperature to virtual temperature.
      !
      ! Input:
      !    t     Temperature (degrees Kelvin)
      !    q     Water vapor mixing ratio (kg H2O/kg air)
      !
      ! Output:
      !    t2tv  Virtual temperature (degrees Kelvin)
      !==================================================
      real function t2tv (t, q)
         implicit none
         real, parameter  :: epsil = .62198
         real, intent(in) :: t
         real, intent(in) :: q
         t2tv = t * (epsil+q) / (epsil*(1.+q))
!        t2tv = t * (1. + ((1.-epsil)/epsil)*q) ! approximation
!        t2tv = t * (1. + .607743*q)            ! approximation
      end function t2tv
