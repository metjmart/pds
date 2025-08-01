      subroutine evalgradErrYZ(errMax,pofyz,uofyz,tofyz,qofyz,
     &                         z,ny,nz,dy,fofy)

         implicit none

         real, dimension(ny,nz), intent(in) :: pofyz,tofyz,uofyz,qofyz
         real, dimension(nz)   , intent(in) :: z
         real, dimension(ny)   , intent(in) :: fofy
         real                  , intent(in) :: dy
         integer               , intent(in) :: ny, nz

         real, intent(out) :: errMax

         integer :: j, k

         real :: tempv, gradErr
         real, external :: t2tv

         real, parameter :: gasConst = 287.04

         errMax = 0.

         do j=2,ny-1
            do k=1,nz

               tempv = t2tv (tofyz(j,k), qofyz(j,k))

!              Compute error in terms of acceleration:

               gradErr = ((tempv*gasConst)/pofyz(j,k))*
     &                   (pofyz(j+1,k) - pofyz(j-1,k))/(2*dy) +
     &                   fofy(j)*uofyz(j,k)

               errMax = max(errMax, abs(gradErr))

            end do
         end do

         return
      end subroutine evalgradErrYZ
