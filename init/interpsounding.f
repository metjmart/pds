
!   This subroutine interpolates the sounding data to the 
!   defined z levels.

      subroutine interpsounding(z,tsound,qsound,usound,vsound, nz,
     &                          zdata,tdata,qdata,udata,vdata, dl)

      implicit none

      integer :: nz
      integer :: dl  ! datalength
      real, dimension(nz) :: z,tsound,qsound,usound,vsound
      real, dimension(dl) :: zdata,tdata,qdata,udata,vdata
      
      integer :: k

      real :: zloc

      real, external :: linterp
      real, external :: lextrap

      do k=1,nz

        zloc=z(k)

        if (zloc<zdata(1)) then  ! level is below sounding data range
          tsound(k)=lextrap(zdata(2),tdata(2),zdata(1),tdata(1),zloc)
          qsound(k)=lextrap(zdata(2),qdata(2),zdata(1),qdata(1),zloc)
          usound(k)=lextrap(zdata(2),udata(2),zdata(1),udata(1),zloc)
          vsound(k)=lextrap(zdata(2),vdata(2),zdata(1),vdata(1),zloc)
        
        else if (zloc > zdata(dl)) then  ! level is above
          tsound(k)=lextrap(zdata(dl-1),tdata(dl-1),zdata(dl),
     &                      tdata(dl),zloc)
          qsound(k)=lextrap(zdata(dl-1),qdata(dl-1),zdata(dl),
     &                      qdata(dl),zloc)
          usound(k)=lextrap(zdata(dl-1),udata(dl-1),zdata(dl),
     &                      udata(dl),zloc)
          vsound(k)=lextrap(zdata(dl-1),vdata(dl-1),zdata(dl),
     &                      vdata(dl),zloc)

        else
          tsound(k)=linterp(dl,zdata,tdata,zloc)
          qsound(k)=linterp(dl,zdata,qdata,zloc)
          usound(k)=linterp(dl,zdata,udata,zloc)
          vsound(k)=linterp(dl,zdata,vdata,zloc)

        end if

      end do

      return

      end subroutine interpsounding

