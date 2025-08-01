
!        Given a point x,y, at level number zlev, this puts the
!        appropiate values into pres, temp, xvel, and yvel
!        given the data in the arrays pforz, vtofrz, tofrz.
!         We assume the center of the vortex is at [x,y]=[0,0]


        subroutine getdatazlev(xloc,yloc,zlev,xvel,yvel,pres,temp,
     &                                pofrz,vtofrz,tofrz,nr,nz,dr)

        implicit none

        integer :: nz,nr,zlev
        real :: xvel,yvel,pres,xloc,yloc,temp, dr
        real, dimension(nr,nz) :: pofrz, vtofrz, tofrz

        real :: pi,pip1,vi,vip1,ri,rip1, ti, tip1
        real :: angle,v,r

        integer :: i

c        First, figure out v, p, and t by linear interpolation:

        r=sqrt(xloc*xloc + yloc*yloc)

        i = int(r/dr) + 1

        if (i.gt.nr-1) then
          i=nr-1
          r=nr*dr
        endif

        ri   = dr * (i-1)
        rip1 = dr * (i  )

        pi=pofrz(i,zlev)
        pip1=pofrz(i+1,zlev)

        vi=vtofrz(i,zlev)
        vip1=vtofrz(i+1,zlev)

        ti=tofrz(i,zlev)
        tip1=tofrz(i+1,zlev)

        pres=(1.0 - ((rip1 - r)/dr))*pip1 + (1.0 - ((r - ri)/dr))*pi

        v=(1.0 - ((rip1 - r)/dr))*vip1 + (1.0 - ((r - ri)/dr))*vi

        if (xloc.eq.0.0) then
           if (yloc.ge.0.0) then
               angle = 3.1415927/2.0
           else
              angle = -3.1415927/2.0
           end if
        else
          angle = atan(yloc/xloc)
        end if

        xvel=-v*sin(angle)*sign(1.0,xloc)
        yvel=v*cos(angle)*sign(1.0,xloc)

        temp=(1.0 - ((rip1 - r)/dr))*tip1 + (1.0 - ((r - ri)/dr))*ti


        return
        end


