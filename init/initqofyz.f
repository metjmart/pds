      subroutine initqofyz (qofyz, qsound, z, ny, nz)

         integer            , intent(in) :: ny, nz
         real, dimension(nz), intent(in) :: qsound
         real, dimension(nz), intent(in) :: z

         real, dimension(ny,nz), intent(out) :: qofyz

         integer :: j, k

         do k=1,nz
            do j=1,ny
                qofyz(j,k) = qsound(k)
            end do
         end do

         return
      end subroutine initqofyz
