      subroutine updatetofyz (tofyz, qofyz, pofyz, z, ny, nz)

         implicit none

         integer               , intent(in) :: ny,nz
         real, dimension(ny,nz), intent(in) :: pofyz, qofyz
         real, dimension(nz)   , intent(in) :: z

         real, dimension(ny,nz), intent(out) :: tofyz

         real :: dpdz
         real :: dz, p

         integer j, k

         real, parameter :: gasConst = 287.04
         real, parameter :: g = 9.80665

         real, external :: t2tv

!        This routine only works for regularly spaced z levels
         dz = z(2) - z(1)

         do k=1,nz
            do j=1,ny

               p = pofyz(j,k)

               if (k == 1) then

                  dpdz = (pofyz(j,3) - pofyz(j,1))/(2.*dz)
     &                - dz*(pofyz(j,3)+pofyz(j,1)-2.*pofyz(j,2))/(dz*dz)

               else if (k == nz) then

                  dpdz = (pofyz(j,nz) - pofyz(j,nz-2))/(2.*dz)
     &         + dz*(pofyz(j,nz)+pofyz(j,nz-2)-2.*pofyz(j,nz-1))/(dz*dz)

               else

                  dpdz = (pofyz(j,k+1) - pofyz(j,k-1))/(2.*dz)

               end if

               tofyz(j,k) = (-p*g) / (gasConst * t2tv(dpdz, qofyz(j,k)))

            end do
         end do

         return
      end subroutine updatetofyz
