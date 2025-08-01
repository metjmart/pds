!     This subroutine generates a matrix field for
!     the combined hydrostatic and barotropic pressure
!     as a function of r and z. Can be used as initializer
!     or as iterator.

      subroutine makepofrz(pofrz, tofrz, vtofrz, qsound, phyd,
     &                     z, nr, nz, dr, f)

         implicit none

         integer               , intent(in) :: nr, nz
         real, dimension(nr,nz), intent(in) :: tofrz, vtofrz
         real, dimension(nz)   , intent(in) :: phyd, z, qsound
         real                  , intent(in) :: dr, f

         real, dimension(nr,nz), intent(out) :: pofrz

         integer :: i, j

         real :: ri,rip1, ai,aip1

         real, external :: t2tv

         real, parameter :: gasConst = 287.04

!        print *, f, dr, gasConst

!        loop through the levels:

         do j=1,nz

            pofrz(nr,j) = phyd(j)

!           Knowing tofrz, we can use a Crank-Nicholson formulation
!           to integrate p inwards at each level:

            do i=nr-1, 1, -1

               ri   = (i-1) * dr
               rip1 = (i  ) * dr

               if (ri /= 0.) then
                  ai = (f*vtofrz(i,j) + (vtofrz(i,j)*vtofrz(i,j))/ri)/
     &                 ( gasConst* t2tv(tofrz(i,j), qsound(j)) )
               else
                  ai = 0.
               end if

               aip1 = (f*vtofrz(i+1,j)
     &              + (vtofrz(i+1,j)*vtofrz(i+1,j))/rip1)/
     &                ( gasConst * t2tv(tofrz(i+1,j), qsound(j)) )

               pofrz(i,j)=pofrz(i+1,j)*(1. - .5*dr*aip1)/(1. + .5*dr*ai)

!              if (j == 10) then
!                 print *,i,vtofrz(i,j),tofrz(i,j),pofrz(i,j)
!              end if

            end do

         end do

         return
      end subroutine makepofrz
