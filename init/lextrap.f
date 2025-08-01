      !================================================================
      ! Linearly extrapolate to obtain a value for a dependent variable
      ! in a data-void region.
      !
      !    lextrap = Yhi  ---------------  Zhi
      !
      !              Ymid ---------------  Zmid             Ymid-Ylo
      !              Yavg  - - - - - - -   Zavg     slope = --------
      !              Ylo  ---------------  Zlo              Zmid-Zlo
      !
      ! Input:
      !    z-        Independent vertical coordinate values at low,
      !              middle, and high levels
      !    y-        Values of dependent variables at low and middle
      !              levels
      ! Output:
      !    Yhi       Value of the dependent variable at the high level
      !
      ! Revision history:
      !    Date        Programmer                 Description of change
      !    ----        ----------                 ---------------------
      !    12/2011     Craig Mattocks, UM-RSMAS   Wrote original code
      !================================================================
      real function lextrap (Zlo,Ylo, Zmid,Ymid, Zhi) result (Yhi)
         implicit none

         real, intent(in) :: Zlo, Zmid, Zhi
         real, intent(in) :: Ylo, Ymid

         real :: Zavg, Yavg
         real :: slope

         Zavg = .5 * (Zmid + Zlo)
         Yavg = .5 * (Ymid + Ylo)

         slope = (Ymid - Ylo) / (Zmid - Zlo)

         Yhi  = Yavg + slope * (Zhi - Zavg)
!        Yhi  = Ymid + slope * (Zhi - Zmid)

         return
      end function lextrap
