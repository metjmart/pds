!        This subroutine generates a balanced axisymmetric, baroclinic
!        vortex. First we define a reference state at the outer
!        boundary. Then we generate the desired velocity field
!        of the vortex. Then we integrate inwards from the outer
!        boundary to compute the pressure field. Then we update
!        the temperature field to enforce hydrostatic balance.
!        The we repeat these last two steps until the error in the
!        pressure gradient balance is really small.

      subroutine makeVortex(pofrz,vtofrz,tofrz,tsound,qsound,phyd,
     &                    vtbar,z,nr,nz,dr,f,p0)

         integer, intent(in) :: nz, nr
         real, dimension(nz) :: z, phyd,tsound,qsound
         real, dimension(nr) :: vtbar
         real, dimension(nr,nz) :: pofrz, vtofrz, tofrz

         real :: dr, p0, f

         integer :: i

         real :: errmax

         integer, parameter :: itMax = 5  ! Number of iterations

         call makevtofrz(vtofrz,vtbar,z,nr,nz,dr)
!        print *, 'makeVortex: done mvtofrz'

         call inittofrz(tsound,tofrz,z,nr,nz)
!        print *, 'makeVortex: done inittofrz'

         do i=1,itMax

            call makepofrz(pofrz,tofrz,vtofrz,qsound,phyd,z,nr,nz,dr,f)
!           print *, 'makeVortex: done makepofrz'

            call updatetofrz(tofrz,pofrz,qsound,z,nr,nz)
!           print *,'makeVortex: done update'

            call evalgraderr(errmax,pofrz,vtofrz,tofrz,
     &                       qsound,z,nr,nz,dr,f)
!           print *, 'makeVortex: done eval'

            write (*, "('makeVortex: on iteration ', i2,
     &                  ' gradient balance errmax = ', 1pg14.7)")
     &            i, errmax

         end do

         return
      end subroutine makeVortex
