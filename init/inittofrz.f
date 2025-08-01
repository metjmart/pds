      subroutine inittofrz(tSound, tofrz, z, nr, nz)

         implicit none

         integer, intent(in) :: nr, nz
         real, dimension(nz), intent(in) :: z, tSound

         real, dimension(nr,nz), intent(out) :: tofrz

         real    :: temp
         integer :: i, k

         do k=1,nz
            temp = tSound(k)
            do i=1,nr
                tofrz(i,k) = temp
            end do
         end do

         return
      end subroutine inittofrz
