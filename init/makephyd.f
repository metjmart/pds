!     This subroutine computes the hydrostatic pressure,
!     based on the input temperature  and moisture profiles
!     as a function of z

!     Note: p0 is a reference presure at the SURFACE, regardless
!     of the location of the lowest level

      subroutine makephyd(tsound, qsound, phyd, z, nz, p0)

         implicit none

         integer            , intent(in ) :: nz
         real, dimension(nz), intent(in ) :: tsound, qsound, z
         real               , intent(in ) :: p0

         real, dimension(nz), intent(out) :: phyd

         integer :: j

         real :: dz, pj,pjm1, tj,tjm1

         real, parameter :: gasConst = 287.04
         real, parameter :: g = 9.80665

         real, external :: t2tv

         phyd(1) = p0

         do j=2,nz

            pjm1 = phyd(j-1)

            dz = z(j) - z(j-1)

            tjm1 = t2tv (tsound(j-1), qsound(j-1))
            tj   = t2tv (tsound(j  ), qsound(j  ))

c           Here we use a semi-implicit scheme to step forward:
            pj = pjm1*(1. - (g*dz)/(2.*gasConst*tjm1) ) /
     &                (1. + (g*dz)/(2.*gasConst*tj  ) )

            phyd(j) = pj

         end do

         return
      end
