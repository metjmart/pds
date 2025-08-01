      subroutine inittofyz(tofyz, tSound, z, ny, nz)

         implicit none

         integer            , intent(in) :: ny, nz
         real, dimension(nz), intent(in) :: tSound
         real, dimension(nz), intent(in) :: z

         real, dimension(ny,nz), intent(out) :: tofyz

         integer :: j, k

         real :: temp

         do k = 1, nz
            temp = tSound(k)
            do j=1,ny
               tofyz(j,k) = temp
            end do
         end do

         return
      end subroutine inittofyz
