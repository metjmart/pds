      subroutine makepofyz(pofyz,tofyz,qofyz,uofyz,phyd,
     &                     z,fofy,ny,nz,dy,jvort)

         implicit none

         integer               , intent(in) :: ny, nz, jvort

         real, dimension(ny,nz), intent(in) :: tofyz, uofyz, qofyz
         real, dimension(nz)   , intent(in) :: phyd, z
         real, dimension(ny)   , intent(in) :: fofy
         real                  , intent(in) :: dy

         real, dimension(ny,nz), intent(out) :: pofyz

         real :: aj, ajp1, ajm1

         integer j, k

         real, external :: t2tv

         real, parameter :: gasConst = 287.04

!        loop through the levels:

         do k=1,nz

            pofyz(jvort,k) = phyd(k)

!           Knowing tofyz, we can use a Crank-Nicholson formulation
!           to integrate p northwards at each level:

            do j=jvort,ny-1

               aj = (fofy(j)*uofyz(j,k)) /
     &              ( gasConst * t2tv(tofyz(j,k), qofyz(j,k)) )

               ajp1 = (fofy(j+1)*uofyz(j+1,k)) /
     &                ( gasConst * t2tv(tofyz(j+1,k), qofyz(j,k)) )

               pofyz(j+1,k)=pofyz(j,k)*(1. - .5*dy*aj)/(1. + .5*dy*ajp1)

            end do

!           Now integrate southward:

            do j=jvort,2,-1

               aj = (fofy(j)*uofyz(j,k)) /
     &              ( gasConst * t2tv(tofyz(j,k), qofyz(j,k)) )

               ajm1 = (fofy(j-1)*uofyz(j-1,k)) /
     &                ( gasConst * t2tv(tofyz(j-1,k), qofyz(j,k)) )

               pofyz(j-1,k)=pofyz(j,k)*(1. + .5*dy*aj)/(1. - .5*dy*ajm1)

            end do

         end do

         return
      end subroutine makepofyz
