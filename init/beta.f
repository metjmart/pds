      !===================================================
      ! Calculate beta (the latitudinal variation of the
      ! Coriolis force) at a given latitude.
      !
      !    Beta = 2 * omega * cos(latitude) / (eathRadius)
      !
      ! Input:
      !    lat          Latitude (degrees north)
      !
      ! Output:
      !    beta         d(Coriolis)/dLat (1/s/m)
      !===================================================
      real function beta (lat)
         real :: pi
         real :: deg2rad
         real :: omega
         real, parameter  :: earthRad = 6378137.
         real, intent(in) :: lat

         pi = 2. * asin(1.)
         deg2rad = pi / 180.
         omega = 2. * pi / 86164.098903691

         beta = 2. * omega * cos(deg2rad * lat) / earthRad
      end function beta
