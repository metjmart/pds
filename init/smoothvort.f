c     This version of smoothVort generates a monopolar vorticity
c     distribution

      real function smoothVort(r)
         implicit none

         real, intent(in) :: r

         real, parameter :: zetaMax = .0015
         real, parameter :: width   = 45000.

         smoothVort = zetaMax * exp( - (r/width)**2 )

         return
      end function smoothVort
