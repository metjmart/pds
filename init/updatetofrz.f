      subroutine updatetofrz(tofrz, pofrz, qsound, z, nr, nz)

         implicit none

         integer               , intent(in) :: nr, nz
         real, dimension(nr,nz), intent(in) :: pofrz
         real, dimension(nz)   , intent(in) :: qsound, z

         real, dimension(nr,nz), intent(out) :: tofrz

         real :: dpdz
         real :: dz, p

         integer :: i,j

         real, parameter :: gasConst = 287.04
         real, parameter :: g = 9.80665

         real, external :: t2tv

!        This routine only works for regularly spaced z levels
         dz = z(2) - z(1)

         do j=1,nz
!           do i=1,nr-1  ! Don't touch outer column
            do i=1,nr    ! Yes, must change outer column too or you get
                         ! discontinuity

               p = pofrz(i,j)

               if (j == 1) then

                  dpdz = (pofrz(i,3) - pofrz(i,1))/(2.*dz)
     &                - dz*(pofrz(i,3)+pofrz(i,1)-2.*pofrz(i,2))/(dz*dz)

               else if (j == nz) then

                  dpdz = (pofrz(i,nz) - pofrz(i,nz-2))/(2.*dz)
     &         + dz*(pofrz(i,nz)+pofrz(i,nz-2)-2.*pofrz(i,nz-1))/(dz*dz)

               else

                  dpdz = (pofrz(i,j+1) - pofrz(i,j-1))/(2.*dz)

               end if

               tofrz(i,j) = (-p*g) / (gasConst * t2tv(dpdz, qsound(j)))

!              if (j == 10) then
!                 print *, i, tofrz(i,j),p ofrz(i,j)
!              end if

            end do
         end do

         return
      end subroutine updatetofrz
