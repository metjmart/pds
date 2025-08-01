      function getq(t, p, rh) result(q)
         implicit none

         integer, parameter :: WRF    = 1
         integer, parameter :: BOLTON = 2
         integer, parameter :: method = WRF

         logical, parameter :: wrt_liquid = .false.

         real, parameter :: tZero  = 273.15   ! 0 deg C
         real, parameter :: tf     = tZero
         real, parameter :: t_ref  = 0.

         real, parameter :: esZero = 6.1121
         real, parameter :: eps    = .62198   ! Rd/Rv = ratio of dry/wet
                                              ! gas constants
         real, parameter :: mw_air = 28.966
         real, parameter :: mw_vap = 18.0152

         real, parameter :: a0     = 6.107799961
         real, parameter :: a1     = 4.436518521e-01
         real, parameter :: a2     = 1.428945805e-02
         real, parameter :: a3     = 2.650648471e-04
         real, parameter :: a4     = 3.031240396e-06
         real, parameter :: a5     = 2.034080948e-08
         real, parameter :: a6     = 6.136820929e-11

         real, parameter :: c1     = 9.09718
         real, parameter :: c2     = 3.56654
         real, parameter :: c3     = 0.876793
         real, parameter :: eis    = 6.1071

         real, parameter :: svp1   = esZero
         real, parameter :: svp2   = 17.67
         real, parameter :: svp3   = tZero - 243.5
         real, parameter :: svp4   = 23.33086
         real, parameter :: svp5   = 6111.72784
         real, parameter :: svp6   = .15215

         real, intent(in) :: t
         real, intent(in) :: p
         real, intent(in) :: rh
         real :: q

         real :: qs
         real :: es
         real :: ew
         real :: rhs
         real :: tk
         real :: t1
         real :: q1

         if (method == WRF) then
            ! Compute mixing ratio q (kg/kg) from pressure p (Pa),
            ! temperature t (K) and relative humidity rh (%).
            ! The reference temperature t_ref (C) is used to describe
            ! the temperature at which the liquid and ice phase change
            ! occurs.

            if (wrt_liquid) then
               es = svp1 * exp(svp2 * (t-tZero)/(t-svp3))
               qs = eps * (es / (.01*p-es))
               q  = max(.01*rh * qs, 0.)
            else
               t1 = t - tZero    ! T in deg C

               if (t1 < -200.) then
                  ! obviously dry
                  q = 0.
               else
                  ! first compute the ambient vapor pressure of water
                  if ( (t1 >= t_ref) .and. (t1 >= -47.) ) then
                     ! liquid phase eslo
                     ew = a0 + t1 * (a1 + t1 * (a2 + t1 * (a3 + t1 *
     &                   (a4 + t1 * (a5 + t1 * a6)))))
                  else if ( (t1 >= t_ref) .and. (t1 < -47.) ) then
                     ! liquid phase poor es
                     ew = svp1 * exp(svp2 * t1 / ( t1 + 243.5))
                  else
                     tk = t
                     rhs = -c1 * (tf / tk - 1.) - c2 * alog10(tf / tk) +
     &                      c3 * (1. - tk / tf) +      alog10(eis)
                     ew = 10. ** rhs
                  end if

                  ! Now that saturation vapor pressure has been obtained
                  ! compute local vapor pressure
                  ew = max(ew , 0.) * rh * 0.01

                  ! Compute the specific humidity using the partial
                  ! vapor pressures of water vapor (ew) and dry air
                  ! (p-ew).
                  ! The constants assume that the pressure is in hPa,
                  ! so we divide the pressures by 100.
                  q1 = mw_vap * ew
                  q1 = q1 / (q1 + mw_air * (p*.01 - ew))
                  q = max(q1 / (1. - q1), 0.)
               end if
            end if

         else if (method == BOLTON) then
            if (t > tZero) then
               es = svp1 * exp((svp2*(t-tZero)) / (t-svp3))
            else
               es = exp(svp4 - svp5/t + svp6*log(t))
            end if
            qs = eps*(es/(.01*p-es))
            q  = max(.01*rh * qs, 0.)
         end if

         return
      end function getq
