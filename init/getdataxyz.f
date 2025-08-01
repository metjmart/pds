
!        Given a point x,y z, this returns the velocity, pressure,
!        and temp values in the vortex. It does linear interpolation
!        from data above and below, on the prescribed z levels of 
!        the axisymmetric fields


        subroutine getdataxyz(xloc,yloc,zloc,xvel,yvel,pres,temp,
     &                                pofrz,vtofrz,tofrz,z,nr,nz,dr)

        implicit none

        integer :: nz,nr
        real :: xvel,yvel,pres,xloc,yloc,zloc,temp,dr
        real, dimension(nr,nz) :: pofrz, vtofrz, tofrz
        real, dimension(nz) :: z

        integer :: zlevabove,zlevbelow
        real :: xvelabove,yvelabove,presabove,tempabove
        real :: xvelbelow,yvelbelow,presbelow,tempbelow
        real :: dz, zabove, zbelow

        integer :: i

!        First, figure out which zlevs are above and below the point:

        i = 1
        look_for_zlevs: do

          i=i+1
          if (z(i).ge.zloc) then
             zlevabove = i
             zlevbelow = i-1
             zabove = z(i)
             zbelow = z(i-1)
             dz = zabove-zbelow
             exit look_for_zlevs
          end if

          if (i.gt.nz) then
             print *,'failure in getdataxyz to find zlevs'
             exit look_for_zlevs
          end if

        end do look_for_zlevs

c        Now get data from above and below:

        call getdatazlev(xloc,yloc,zlevabove,xvelabove,yvelabove,
     &                                presabove,tempabove,
     &                                pofrz,vtofrz,tofrz,nr,nz,dr)

        call getdatazlev(xloc,yloc,zlevbelow,xvelbelow,yvelbelow,
     &                                presbelow,tempbelow,
     &                                pofrz,vtofrz,tofrz,nr,nz,dr)


c        Then interpolate to actual zloc:

        pres=(1.0 - ((zabove - zloc)/dz))*presabove 
     &                        + (1.0 - ((zloc - zbelow)/dz))*presbelow
        
        temp=(1.0 - ((zabove - zloc)/dz))*tempabove 
     &                        + (1.0 - ((zloc - zbelow)/dz))*tempbelow

        xvel=(1.0 - ((zabove - zloc)/dz))*xvelabove 
     &                        + (1.0 - ((zloc - zbelow)/dz))*xvelbelow

        yvel=(1.0 - ((zabove - zloc)/dz))*yvelabove 
     &                        + (1.0 - ((zloc - zbelow)/dz))*yvelbelow


        return
        end

