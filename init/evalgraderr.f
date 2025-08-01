      subroutine evalgradErr(errMax, pofrz, vtofrz, tofrz,
     &                       qsound, z, nr, nz, dr, f)

         integer               , intent(in) :: nr, nz
         real, dimension(nr,nz), intent(in) :: pofrz, tofrz, vtofrz
         real, dimension(nz)   , intent(in) :: z, qsound
         real                  , intent(in) :: dr, f

         real, intent(out) :: errMax

         integer :: i, j

         real :: tempv, gradErr
         real, external :: t2tv

         real, parameter :: gasConst = 287.04

         errMax = 0.

         do j=1,nz
            do i=2,nr-1

               tempv = t2tv (tofrz(i,j), qsound(j))

c              Compute error in terms of acceleration:

               gradErr = ((tempv*gasConst)/pofrz(i,j))*
     &                   (pofrz(i+1,j) - pofrz(i-1,j))/(2*dr) -
     &           (f*vtofrz(i,j) + vtofrz(i,j)*vtofrz(i,j)/((i-1)*dr) )

               errMax = max(errMax, abs(gradErr))

            end do
         end do

         return
      end subroutine evalgradErr
