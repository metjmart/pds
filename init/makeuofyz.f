
!     This subrouting generates idealized U(y,z) fields

      subroutine makeuofyz(uofyz, phyd, z, ny, nz, dy)

         implicit none

         integer,                intent(in) :: ny, nz
         real   ,                intent(in) :: dy
         real   , dimension(nz), intent(in) :: z, phyd

         real, dimension(ny,nz), intent(out) :: uofyz

         integer :: j, k
         real    :: yloc, zloc, p, x
         real    :: degdist

         real :: pi

         real, parameter :: metersPerDegree = 111317.1
         real, parameter :: d = 2.5

         real, parameter :: zcen  =  10000.
         real, parameter :: sy    = 500000.
         real, parameter :: sz    =   5000.

         real, parameter :: surf  = -5 
         real, parameter :: shear = 5
!         real, parameter :: shear = 15

         real :: ycen

         pi = 2. * asin(1.)
         ycen  = .5 * ny * dy

         do j=1,ny

            yloc = dy * j

            do k=1,nz

               zloc = z(k)

!              Use log-pressure-height formula for U for best comparison
               p = phyd(k)

!              Formula for cosine shear; must add v in sheargeneric.f
!              if you want idealized meridional wind
               x = 1.592 * (4.929 - log10(p))
              if (p > 85000.) then
                 uofyz(j,k) = surf
              else if (p < 20000.) then
                 uofyz(j,k) = surf + shear
              else
                 uofyz(j,k) = surf + .5 * shear * (1. - cos(pi*x))
              end if

!              Here is a Gaussian-shaped jet:
!              uofyz(j,k) = 20. * exp( -((yloc-ycen)/sy)**2.
!     &                                -((zloc-zcen)/sz)**2. )

!              Here is the Tuleya and Kurihara barotropic flow:
!              degdist = (yloc - ycen) / metersPerDegree
!              uofyz(j,k) = uofyz(j,k) - 2.5*tanh(degdist/d)

!              if (j == 120) then
!                 print *, 'uofyz: ', k, z(k), phyd(k), uofyz(j,k)
!              end if

            end do
         end do

         return
      end subroutine makeuofyz

