      subroutine updateqofyz(tofyz, qofyz, pofyz, rhsound, ny, nz)

         implicit none

         integer,                   intent(in) :: ny, nz
         real   , dimension(ny,nz), intent(in) :: tofyz, pofyz
         real   , dimension(nz)   , intent(in) :: rhsound

         real, dimension(ny,nz), intent(out) :: qofyz

         real, external :: getq

         integer :: j, k

         do k=1,nz
            do j=1,ny

               qofyz(j,k) = getq(tofyz(j,k), pofyz(j,k), rhsound(k))

!              if (j == 120) then
!                 print *,tofyz(j,k),pofyz(j,k),rhsound(k),qofyz(j,k)
!              end if

            end do
         end do

         return
      end subroutine updateqofyz
