      !==================================================
      ! Calculate the Coriolis force at a given latitude.
      !
      !    Coriolis = 2 * omega * sin(latitude)
      !
      ! Input:
      !    lat          Latitude (degrees north)
      !
      ! Output:
      !    coriolis     Coriolis force (1/s)
      !==================================================
      real function coriolis (lat)
         real :: pi
         real :: deg2rad
         real :: omega
         real, intent(in) :: lat

         pi = 2. * asin(1.)
         deg2rad = pi / 180.
         omega = 2. * pi / 86164.098903691

         coriolis = 2. * omega * sin(deg2rad * lat)
      end function coriolis
