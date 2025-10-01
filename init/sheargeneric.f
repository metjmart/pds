      program sheargeneric

         implicit none

!        Geophysical constants
         real :: pi
         real :: deg2rad
         real, parameter :: earthRad = 6378137.
         real, parameter :: g = 9.80665
         real, parameter :: tZero = 273.15
         real, parameter :: metersPerDegree = 111317.1

         integer, parameter :: iswater = 16
         integer, parameter :: map_proj = 0
         character (len=4), parameter :: mminlu = 'USGS'

         integer :: i,j,k,jp1
         real :: truelat1,truelat2
         integer :: ix ! Dimensions of WRF input data
         integer :: jx ! grid
         integer :: kx ! nz

         real, allocatable :: vtofrz(:,:), pofrz(:,:), tofrz(:,:),
     &                        thetaofrz(:,:), rhoofrz(:,:),
     &                        pvtofrz(:,:), absvortr(:,:),
     &                        absvortz(:,:), nsq(:,:),
     &                        dvdr(:,:), dvdz(:,:), 
     &                        dthetadr(:,:), dthetadz(:,:)

         ! Sounding variables
         real, allocatable :: z(:)
         real, allocatable :: pSound(:)
         real, allocatable :: tSound(:)
         real, allocatable :: tdSound(:)
         real, allocatable :: rhSound(:)
         real, allocatable :: qSound(:)
         real, allocatable :: usound(:)
         real, allocatable :: vsound(:)
         real, allocatable :: vtbar(:)
         real, allocatable :: zdata(:)
         real, allocatable :: tdata(:)
         real, allocatable :: qdata(:)
         real, allocatable :: udata(:)
         real, allocatable :: vdata(:)
         real :: dz

         real :: pres,temp, xLoc,yLoc,zLoc, xvel,yvel

!        Fields for 2D balanced shear flow
         real, allocatable :: uofyz(:,:),vofyz(:,:),pofyz(:,:)
         real, allocatable :: tofyz(:,:),qofyz(:,:)
         real, allocatable :: thetaofyz(:,:),rhofyz(:,:)
         real, allocatable :: fofy(:)
         real :: latVal, bet

!        And here are some for making the WRF fields:
         real, allocatable :: u(:,:,:),v(:,:,:),q(:,:,:),
     &                        p(:,:,:),t(:,:,:),theta(:,:,:)
         real, allocatable :: zetaFull(:), zetaHalf(:)
         real, allocatable :: tmn(:,:),ht(:,:),tsk(:,:),lu_index(:,:),
     &                        xland(:,:),lat(:,:),lon(:,:),msft(:,:),
     &                        msfu(:,:),msfv(:,:),msftx(:,:),
     &                        msfty(:,:),msfux(:,:),msfuy(:,:),
     &                        msfvx(:,:),msfvy(:,:),
     &                        msfvx_inv(:,:),
     &                        sina(:,:),cosa(:,:),e(:,:),f(:,:),
     &                        pSurface(:,:)
!    &,ivgtyp(:,:),isltyp(:,:),

         real    :: dLat, dLon, eval, fval, dZeta
         real    :: tDelta, x, rhfac
         real :: r, rand

         real, external :: getrh
         real, external :: getq
         real, external :: coriolis
         real, external :: beta
         real, external :: t2th

!        Here define values for input variables
         integer :: nx = 361             ! size of WRF domain
         integer :: ny = 361             ! size of WRF domain
         real :: cen_lat = 20.           ! lat at center of domain
         real :: cen_lon = -60.          ! lon at center of domain
         real :: p0 = 101500.0           ! reference surface pressure
         integer :: iVort = 180          ! grid location of initial vortex
         integer :: jVort = 180          ! grid location of initial vortex
         real :: dx = 18000.             ! WRF domain grid spacing
         real :: dy = 18000.             ! WRF domain grid spacing
         real :: zTop = 26000.           ! top of this domain, should be a bit higher than WRF
         real :: dr = 4000.              ! radial spacing for axisymmetric vortex
         integer :: nr = 401             ! radial grid size for vortex
         integer :: nz = 89              ! vertical grid size for data set 
         integer :: ispds = 1            ! set to one for point-downscaling
         integer :: betaPlane = 0        ! set to one for beta plane (real-shear only)
         real :: noiseAmp = 0.5          ! some parameters for adding noise to initial winds
!        real :: noiseAmp = 0.0          ! for TVPDS env soundings
         real :: noiseRad = 100000.      ! central radius of annulus of noise
         real :: noiseWidth = 30000.     ! width of annulus of noise
         integer :: adjQ = 0             ! set to one to adjust RH to constant across domain
         real :: tskval = 301.15         ! domain-wide SST (can make spatially varying below)
         integer :: adjSST = 0           ! adjust SST to follow atmospheric variations
         integer :: adjTlev = 41         ! level at which for SST to follow variations
         integer :: datalength = 14      ! length of the input sounding data
         integer :: idealwinds = 1       ! for defining idealized wind profiles; over-writes data

         truelat1 = cen_lat  ! historical
         truelat2 = cen_lat  ! historical
         ix = nx  ! Dimensions of WRF input data
         jx = ny  ! 
         kx = nz  ! (Historical reasons for these variable names)

!        *** Allocate our arrays now that we have dimension variables
         allocate(vtofrz(nr,nz))
         allocate(pofrz(nr,nz))
         allocate(tofrz(nr,nz))
         allocate(thetaofrz(nr,nz))
         allocate(rhoofrz(nr,nz))
         allocate(pvtofrz(nr,nz))
         allocate(absvortr(nr,nz))
         allocate(absvortz(nr,nz))
         allocate(nsq(nr,nz))
         allocate(dvdr(nr,nz))
         allocate(dvdz(nr,nz))
         allocate(dthetadr(nr,nz))
         allocate(dthetadz(nr,nz))
         allocate(z(nz))
         allocate(pSound(nz))
         allocate(qSound(nz))
         allocate(tSound(nz))
         allocate(rhSound(nz))
         allocate(usound(nz))
         allocate(vsound(nz))

         allocate(zdata(datalength))
         allocate(tdata(datalength))
         allocate(qdata(datalength))
         allocate(udata(datalength))
         allocate(vdata(datalength))

         allocate(vtbar(nr))
         allocate(uofyz(ny,nz))
         allocate(vofyz(ny,nz))
         allocate(pofyz(ny,nz))
         allocate(tofyz(ny,nz))
         allocate(qofyz(ny,nz))
         allocate(thetaofyz(ny,nz))
         allocate(rhofyz(ny,nz))
         allocate(fofy(ny))
         allocate(u(ix,jx,kx))
         allocate(v(ix,jx,kx))
         allocate(q(ix,jx,kx))
         allocate(p(ix,jx,kx))
         allocate(t(ix,jx,kx))
         allocate(theta(ix,jx,kx))
         allocate(zetaFull(kx))
         allocate(zetaHalf(kx))
         allocate(tmn(ix,jx))
         allocate(ht(ix,jx))
         allocate(tsk(ix,jx))
         allocate(lu_index(ix,jx))
!        allocate(ivgtyp(ix,jx))
!        allocate(isltyp(ix,jx))
         allocate(xland(ix,jx))
         allocate(lat(ix,jx))
         allocate(lon(ix,jx))
         allocate(msft(ix,jx))
         allocate(msfu(ix,jx))
         allocate(msfv(ix,jx))
         allocate(msftx(ix,jx))
         allocate(msfty(ix,jx))
         allocate(msfux(ix,jx))
         allocate(msfuy(ix,jx))
         allocate(msfvx(ix,jx))
         allocate(msfvy(ix,jx))
         allocate(msfvx_inv(ix,jx))
         allocate(sina(ix,jx))
         allocate(cosa(ix,jx))
         allocate(e(ix,jx))
         allocate(f(ix,jx))
         allocate(pSurface(ix,jx))


      print *,              'cen_lat = ' , cen_lat
      print *,              'cen_lon = ' , cen_lon
      print *,              'iVort = ' , iVort
      print *,              'jVort = ' , jVort
      print *,              'nx = ' , nx
      print *,              'ny = ' , ny
      print *,              'dx = ' , dx
      print *,              'zTop = ' , zTop
      print *,              'dr = ' , dr
      print *,              'nr = ' , nr
      print *,              'nz = ' , nz
      print *,              'ispds = ' , ispds
      print *,              'betaPlane = ' , betaPlane
      print *,              'noiseAmp = ' , noiseAmp
      print *,              'noiseRad = ' , noiseRad
      print *,              'noiseWidth = ' , noiseWidth
      print *,              'adjQ = ' , adjQ
      print *,              'tskval = ' ,  tskval
      print *,              'adjSST = ' , adjSST
      print *,              'adjTlev = ' , adjTlev
      print *,              'idealwinds = ' , idealwinds
      print *,              'datalength = ' , datalength

!     Now read in sounding data:

        open(unit=53,file='sounding.input',status='unknown')
        do i = 1,datalength
          read(53,*) zdata(i),tdata(i),qdata(i),udata(i),vdata(i)
          print *, zdata(i),tdata(i),qdata(i),udata(i),vdata(i)
        end do

!        *** All data is now read in...
!**********************************************************************************************

         pi = 2. * asin(1.)
         deg2rad = pi / 180.

!        Axisymmetric vortex output:
         open (unit=20, file='z.out'           , blank='null')
         open (unit=60, file='vofz.out'        , blank='null')
         open (unit=21, file='pofrz.out'       , blank='null')
         open (unit=22, file='tofrz.out'       , blank='null')
         open (unit=23, file='vtofrz.out'       , blank='null')
         open (unit=24, file='thetaofrz.out'   , blank='null')
         open (unit=25, file='rhoofrz.out'     , blank='null')
         open (unit=26, file='pvtofrz.out'      , blank='null')

!        Balanced shear flow output:
         open (unit=27, file='uofyz.out'       , blank='null')
         open (unit=28, file='pofyz.out'       , blank='null')
         open (unit=29, file='tofyz.out'       , blank='null')
         open (unit=30, file='qofyz.out'       , blank='null')
         open (unit=31, file='thetaofyz.out'   , blank='null')
         open (unit=32, file='rhofyz.out'      , blank='null')

!        Files which will be WRF input
         open (unit=12, file='input_sounding'  , blank='null')
         open (unit=14, file='input_sounding_v', blank='null')
         open (unit=15, file='input_sounding_u', blank='null')
         open (unit=16, file='input_sounding_o', blank='null')

         dz = zTop / (nz-1)
         print *, 'axisymmetric vortex dr = ', dr, ' dz = ', dz

!        Define altitudes of data levels for the vortex and
!        the balanced shear
         do k=1,nz
            z(k) = (k-1) * dz
         end do

!        This value of f is for center of vortex and domain:
         fval = coriolis(cen_lat)
!        fval = 5.0000e-5
!        eval = beta(cen_lat) * earthRad
         eval = 0.
         bet = beta(cen_lat)
         print*, 'tskval, fval, eval, beta = ', tskval, fval, eval, bet

!        In the first part of this program, we use subroutines to
!        generate an axisymmetric, balanced, hurricane-like vortex
!        First thing is to get reference presssure profile at center of domain
!        T and q profiles no longer live in the sounding function

!        Here we interpolate the data from the sounding file to the data levels:

         call interpsounding(z,tsound,qSound,usound,vsound,nz,
     &                       zdata,tdata,qdata,udata,vdata,datalength)
         print *,'done interpsounding'

!        Now compute pressure profile at center of domain

         call makephyd(tsound,qSound,psound,z,nz,p0)
         print *,'done makephyd'

!        compute tDelta in case it is needed to adjust SST
         tDelta = tskval - tSound(adjTlev)
         print *, 'tDelta = ', tDelta

!        Compute RH profile of central sounding
         print *,'Recomputed sounding:'
         do k = nz, 1, -1
            rhSound(k) = getrh(tSound(k), pSound(k), qSound(k))
            write (*, "('z = ', f8.2, 3x, 'p = ', f9.2, 3x,
     &                  't = ', f6.2, 3x, 'q = ', f9.6, 3x,
     &               'rh = ', f6.2,3x,'u = ',f7.2,3x,'v = ',f7.2,3x)")
     &            z(k), pSound(k), tSound(k), qSound(k), rhSound(k),
     &                  usound(k), vsound(k)
         end do


!        **************************
!        Some code to modify q/rh profiles:
!        ***** comment out if you want the original RH profile or
!        make amplitude of bump be 0.0
!        do i=1,2   ! iterate twice

!           do k=1,nz
!              rhSound(k) = getrh(tSound(k), pSound(k), qSound(k))
!           end do

!           do k = 1,nz
!              if ((pSound(k) <= 85000.) .and. (pSound(k) >= 40000.)) then
!                 Gaussian-shaped bump of RH around 600 mb
!                 rhfac = 1. + .2*exp( -((pSound(k) - 60000.)/15000.)**2 )
!                 rhSound(k) = min(rhSound(k)*rhfac, 95.)  ! limit initial RH to 95%
!                 qSound (k) = getq(tSound(k), pSound(k), rhSound(k))
!              end if
!           end do

!           call makePhyd(tSound,qSound,pSound,z,nz,p0)

!        end do

!        Re-compute RH profile of central sounding
!        do k=1,nz
!           rhSound(k)=getrh(tSound(k),pSound(k),qSound(k))
!           print *,'z=',z(k),'p=',pSound(k),'rh=',rhSound(k)
!        end do
!        End modify RH code
!        **********************************

         call makeVortex(pofrz,vtofrz,tofrz,tSound,qSound,pSound,
     &                vtbar,z,nr,nz,dr,fval,p0)
         print *,'done makeVortex'

!        Compute density, potential temp
         do k=1,nz
            do i=1,nr
               rhoofrz  (i,k) = pofrz(i,k) / (287.04 * tofrz(i,k))
               thetaofrz(i,k) = t2th (tofrz(i,k), pofrz(i,k))
            end do
         end do
         print *,'done computing density and theta'

!        Compute derivative quantities. Only works for regularly spaced.
         do k=2,nz-1
            do i=2,nr-1
               dvdz(i,k) = (vtofrz(i,k+1) - vtofrz(i,k-1))/(2.*dz)
               dvdr(i,k) = (vtofrz(i+1,k) - vtofrz(i-1,k))/(2.*dr)
               absvortz(i,k) = fval +
     &                         (vtofrz(i+1,k) - vtofrz(i-1,k))/(2.*dr)
     &                                       + vtofrz(i,k)/(i*dr)
               absvortr(i,k) = -dvdz(i,k)
               dthetadr(i,k) = (thetaofrz(i+1,k)-thetaofrz(i-1,k)) /
     &                         (2.*dr)
               dthetadz(i,k) = (thetaofrz(i,k+1)-thetaofrz(i,k-1)) /
     &                         (2.*dz)
               nsq(i,k) = (g/thetaofrz(i,k)) * dthetadz(i,k)
            end do
         end do
         print *,'done computing derivatives'

!        Compute symmetric PV
         do k=1,nz
            do i=1,nr
               pvtofrz(i,k)=(absvortr(i,k)*dthetadr(i,k)
     &                    + absvortz(i,k)*dthetadz(i,k)) / rhoofrz(i,k)
               if (i ==  1) pvtofrz(i,k) = pvtofrz(2   ,k)
               if (i == nr) pvtofrz(i,k) = pvtofrz(nr-1,k)
               if (k ==  1) pvtofrz(i,k) = pvtofrz(i   ,2)
               if (k == nz) pvtofrz(i,k) = pvtofrz(i,nz-1)
            end do
         end do
         print *,'done computing symmetriic PV'

!        New code: Now, change axisymmetric thermodynamic vortex data
!        to a perturbation about the basic state.
         do k=1,nz
            do i=1,nr
               pofrz(i,k) = pofrz(i,k) - pofrz(nr,k)
               tofrz(i,k) = tofrz(i,k) - tofrz(nr,k)
               thetaofrz(i,k) = thetaofrz(i,k) - thetaofrz(nr,k)
            end do
         end do
         print *,'done computing perturbations'

!        Print these out to ascii files for examination:
         do k=1,nz
100         format(1000(1x,1pe12.5))
            write(20,100) (z(k))
            write(60,100) (vsound(k))
            write(21,100) (pofrz(i,k),i=1,nr)
            write(22,100) (tofrz(i,k),i=1,nr)
            write(23,100) (vtofrz(i,k),i=1,nr)
            write(24,100) (thetaofrz(i,k),i=1,nr)
            write(25,100) (rhoofrz(i,k),i=1,nr)
            write(26,100) (pvtofrz(i,k),i=1,nr)
         end do
         print *,'done writing diagnostic files'

!        go to 1000

!        Second part of program: Compute a balanced shear flow in (y,z)

!        If beta plane,
!        get Coriolis as function of y
         dLat = dy / metersPerDegree
         do j=1,ny
            latVal = cen_lat + dLat*(-(ny-1)/2. + j) ! unstaggered in y
            yLoc   = 0.      + dy  *(-(ny-1)/2. + j) ! unstaggered in y
            if (betaPlane == 1) then
!              fofy(j) = coriolis(latVal)   ! not really "beta plane"
               fofy(j) = fval + bet*yLoc    ! really beta plane
            else
               fofy(j) = fval
            end if
!           print *,j,'fofy=',fofy(j)
         end do
         print *,'done computing coriolis force'

!        Central sounding is same as for outer boundary of vortex.

!        Now make zonal wind field and balance it.

         call makeShear(uofyz,pofyz,tofyz,qofyz,
     &                  tSound,qSound,rhSound,pSound,z,
     &                  usound,idealwinds,
     &                  ny,nz,dy,fofy,jVort,ispds,adjQ)
         print *,'done making shear'

         do k=1,nz
            do j=1,ny
               thetaofyz(j,k) = t2th (tofyz(j,k), pofyz(j,k))
               rhofyz(j,k)=getrh(tofyz(j,k), pofyz(j,k), qofyz(j,k))
            end do
         end do

!        Print these out to ascii files for examination:
         do k=1,nz
            write(27,100) (uofyz(j,k),j=1,ny)
            write(28,100) (pofyz(j,k),j=1,ny)
            write(29,100) (tofyz(j,k),j=1,ny)
            write(30,100) (qofyz(j,k),j=1,ny)
            write(31,100) (thetaofyz(j,k),j=1,ny)
            write(32,100) (rhofyz(j,k),j=1,ny)
         end do
         print *,'done writing diagnostic files'

!        go to 1000

!        Third part of program is to generate 3D fields for WRF input

!        NOTE: the use here of half levels and full levels is legacy
!        but is also consistent with the WRF model, so we'll
!        keep for now

!        dZeta is vertical grid spacing on WRF data grid
         dZeta = zTop / (kx-1)
         print *, 'dZeta = ', dZeta

!        NOTE: for now we are assuming WRF data grid vertical levels
!        Are half levels for the balanced shear flow grid
         zetaFull( 1) = 0.
         zetaHalf(kx) = 0.
         do k=2,kx
            zetaFull(k) = zetaFull(1) + dZeta * (k-1)
!           print *, 'zfull=', k, zetaFull(k), zetaFull(k)-zetaFull(k-1)
         end do

!        calculate half levels based on full levels
         do k = 1,kx-1
            zetaHalf(k) = .5 * (zetaFull(k) + zetaFull(k+1))
         end do

!        This is code to compute full levels from half levels:
!        zetaFull(1) = 0
!        do k = 2,nz
!           zetaFull(k)=zetaFull(k-1)+2.*(zetaHalf(k-1)-zetaFull(k-1))
!        end do


!        Now, initialize surface data:

!        xland (real number) 1.0 for land, 2.0 for water (get it?)
!        tsk surface temperature (land or ocean)
!        tmn temperature at bottom of soil. This is not getting into
!           the model now but we may figure it out later.
!        ht  terrain height
!        lu_index Land use index. See USGS data in LANDUSE.TBL in run directory

         do j=1,jx
            do i=1,ix

!              Stuff that should always (?) be constant everywhere:
               sina(i,j) = 0.
               cosa(i,j) = 1.
               e(i,j) = eval
               f(i,j) = fval
               msft(i,j) = 1.
               msfu(i,j) = 1.
               msfv(i,j) = 1.
               msftx(i,j) = 1.
               msfux(i,j) = 1.
               msfvx(i,j) = 1.
               msfvx_inv(i,j) = 1.
               msfty(i,j) = 1.
               msfuy(i,j) = 1.
               msfvy(i,j) = 1.
               ht(i,j) = 0.         
               xland(i,j)  = 2
               lu_index(i,j) = 16
               if (adjSST == 1) then
                 tsk(i,j) = tofyz(j,adjTlev) + tDelta
                 tmn(i,j) = tofyz(j,adjTlev) + tDelta
                 if (i == 5) print *,'j = ', j, ' tsk = ', tsk(i,j)
               else
                 tsk(i,j) = tskval
                 tmn(i,j) = tskval
               end if
            
            end do
         end do

!        initialize the lat and lon grids on mass grid (1:nx-1,1:ny-1)
         dLat = dy / metersPerDegree
         dLon = dx / metersPerDegree !/ cos(deg2rad * cen_lat)

         print *,cen_lat,cen_lon,dLat,dLon

         do j= 1,jx-1
            do i = 1,ix-1

               lat(i,j) = cen_lat + dLat*( -jx/2. + j )
               lon(i,j) = cen_lon + dLon*( -ix/2. + i )

               yLoc = 0. + dy*( -jx/2. + j )

               if (betaPlane == 1) then
!                 f(i,j) = coriolis(lat(i,j))   ! not really beta plane
                  f(i,j) = fval + bet * yLoc    ! really beta plane
               end if

!              if (i == 10) then
!                 print *,'latlon:',i,j,lat(i,j),lon(i,j),f(i,j)
!              end if

            end do
         end do

!        go to 1000


!        Now, fill 3D arrays. Have to account for horizontal staggers.
!        AND vertical staggers

!        First we will loop over all points to fill in the balanced shear
!        flow; we will put the vortex in afterwards with its own ix,jx,kx loops

         do j=1,jx-1
            jp1 = min(j+1, jx)
            do i=1,ix

               pSurface(i,j) = .5 * (pofyz(j,1) + pofyz(jp1,1))

               do k=1,kx-1

!                 Thermo and zonal wind data are at cell centers of (y,z) grid
!                 REMEMBER for now we are assuming that WRF data grid vertical
!                 levelsi (zetaHalf) are at centers of "z" levels

                  p(i,j,k) = .25 * ( pofyz(j,k  ) + pofyz(jp1,k  )
     &                             + pofyz(j,k+1) + pofyz(jp1,k+1) )

                  t(i,j,k) = .25 * ( tofyz(j,k  ) + tofyz(jp1,k  )
     &                             + tofyz(j,k+1) + tofyz(jp1,k+1) )

                  q(i,j,k) = .25 * ( qofyz(j,k  ) + qofyz(jp1,k  )
     &                             + qofyz(j,k+1) + qofyz(jp1,k+1) )

                  u(i,j,k) = .25 * ( uofyz(j,k  ) + uofyz(jp1,k  )
     &                             + uofyz(j,k+1) + uofyz(jp1,k+1) )

!              Putting in meridional wind: It cannot vary spatially.
!              Either use vsound from file, or idealized or nothing.
!              Must be nothing for realshear/channel.

                  if (idealwinds.eq.1) then

!                   No v:
                    v(i,j,k) = 0.
                    !v(i,j,k) = 5.0    
!                   Or pressure-dependent v for point-downscaling:
!                   (this is rotating shear vector)
!                    x = 1.592*(4.929 - log10(p(i,j,k)))
!                    if (p(i,j,k) > 85000.) then
!                       v(i,j,k) = 0.
!                    else if (p(i,j,k) < 20000.) then
!                       v(i,j,k) = 0.
!                    else
!                       v(i,j,k) = 1. * sin(pi*x)
!                    end if

                  else       ! install vsound(k):

                    v(i,j,k) = vsound(k)

                  endif

                end do

              end do
           end do

!        Now, add the vortex wind fields and temperature
!        and pressure perturbations

         call RANDOM_SEED()

         do j=1,jx-1
            do i=1,ix

               ! Decide where x and y are, based on
               ! array type. Then, go up from bottom to top

               ! First, the mass grid

               xLoc = dx * ( -iVort + i )
               yLoc = dy * ( -jVort + j )

               do k=1,kx-1
                  zLoc = zetaHalf(k)
                  call getdataxyz(xLoc,yLoc,zLoc,xvel,yvel,pres,temp,
     &                            pofrz,vtofrz,tofrz,z,nr,nz,dr)

                  p(i,j,k) = p(i,j,k) + pres

                  t(i,j,k) = t(i,j,k) + temp

!                 WRF data set needs theta, not temp:
                  theta(i,j,k) = t2th (t(i,j,k), p(i,j,k))
               end do

!              print *,'t grid: ',i,j, xLoc,yLoc

               zLoc = 0.
               call getdataxyz(xLoc,yLoc,zLoc,xvel,yvel,pres,temp,
     &                         pofrz,vtofrz,tofrz,z,nr,nz,dr)
   
               pSurface(i,j) = pSurface(i,j) + pres
            end do
         end do

            do j=1,jx
            do i=1,ix

               ! Decide where x and y are, based on
               ! array type. Then, go up from bottom to top

               ! Now, U grid

               xLoc = dx * ( -iVort - .5 + i )
               yLoc = dy * ( -jVort + j )
               r = sqrt(xLoc*xLoc + yLoc*yLoc)

               do k=1,kx-1
                  zLoc = zetaHalf(k)

                  call getdataxyz(xLoc,yLoc,zLoc,xvel,yvel,pres,temp,
     &                            pofrz,vtofrz,tofrz,z,nr,nz,dr)

                  u(i,j,k) = u(i,j,k) + xvel
                  
                  call RANDOM_SEED()
                  call RANDOM_NUMBER(rand)
                  
                     u(i,j,k) = u(i,j,k) + (rand-.5) * noiseAmp *
     &                       exp( -((r - noiseRad)/noiseWidth)**2. )*
     &                       exp( - (2.*zLoc/10000)**3. )
               end do
!              print *,'u grid: ', i,j, xLoc,yLoc

               ! Last, V grid

               xLoc = dx * ( -iVort + i )
               yLoc = dy * ( -jVort - .5 + j )
               r = sqrt(xLoc*xLoc + yLoc*yLoc)

               do k=1,kx-1
                  zLoc = zetaHalf(k)

                  call getdataxyz(xLoc,yLoc,zLoc,xvel,yvel,pres,temp,
     &                            pofrz,vtofrz,tofrz,z,nr,nz,dr)

                   v(i,j,k) = v(i,j,k) + yvel

                  call RANDOM_SEED()
                  call RANDOM_NUMBER(rand)

                     v(i,j,k) = v(i,j,k) + (rand-.5)*noiseAmp*
     &                        exp(-((r - noiseRad)/noiseWidth)**2. )*
     &                        exp( - (2.*zLoc/10000)**3. )
               end do

!               print *,'v grid: ',i,j, xLoc,yLoc, v(i,j,1)

            end do
         end do

!        *************************************************************

!        write out data for WRF input_sounding
         print *,'Writing out data for WRF input_sounding...'

         do j=1,jx-1
            do i=1,ix-1
               write(12,32) i,j,pSurface(i,j)/100.,
     &                      t2th(tsk(i,j), pSurface(i,j)),
     &                      tsk(i,j), q(i,j,1)*1000.
 32            format(i4,1x,i4,3(f8.2,1x))
 33            format(3(f8.2,3x),2(f10.5,3x))
               do k=1,kx-1
                  write(12,33) zetaHalf(k),theta(i,j,k),q(i,j,k)*1000.,
     &                         u(i,j,k),v(i,j,k)
               end do
               write(12,33) -999., -999., -999., -999., -999.
            end do
         end do

         ! write v

         do j=1,jx
            do i=1,ix-1
               do k=1,kx-1
                  write(14,34) v(i,j,k)
 34              format(f10.5)
               end do
            end do
         end do

         ! write u

         do j=1,jx-1
            do i=1,ix
               do k=1,kx-1
                  write(15,34) u(i,j,k)
               end do
            end do
         end do

! write f e cosa sina msft msfu msfv

         write(16,36)cen_lat,cen_lon,truelat1,truelat2,map_proj,iswater
 36      format(f7.3,1x,f8.3,1x,2(f7.3,1x),i2,1x,i2)

         do j=1,jx
            do i=1,ix
               write(16,35) e(i,j),f(i,j),
     &                      msft(i,j),msfu(i,j),msfv(i,j),
     &                      msftx(i,j),msfux(i,j),msfvx(i,j),
     &                      msfty(i,j),msfuy(i,j),msfvy(i,j),
     &                      msfvx_inv(i,j),
     &                      sina(i,j),cosa(i,j), lat(i,j),lon(i,j),
     &                      xland(i,j),lu_index(i,j),ht(i,j),tmn(i,j),
     &                      tsk(i,j)

35             format(2(e14.7,1x),9(f7.3,1x),6(f9.3,1x))

               if ((i == 1) .and. (j == 1)) then
                  write(*,35) e(i,j),f(i,j),
     &                        msft(i,j),msfu(i,j),msfv(i,j),
     &                        msftx(i,j),msfux(i,j),msfvx(i,j),
     &                        msfty(i,j),msfuy(i,j),msfvy(i,j),
     &                        msfvx_inv(i,j),
     &                        sina(i,j),cosa(i,j), lat(i,j),lon(i,j),
     &                        xland(i,j),lu_index(i,j),ht(i,j),tmn(i,j),
     &                        tsk(i,j)

                  print *, e(i,j),f(i,j),msft(i,j),msfu(i,j),msfv(i,j),
     &                     msftx(i,j),msfux(i,j),msfvx(i,j),
     &                     msfty(i,j),msfuy(i,j),msfvy(i,j),
     &                     msfvx_inv(i,j),
     &                     sina(i,j),cosa(i,j), lat(i,j),lon(i,j),
     &                     xland(i,j),lu_index(i,j),ht(i,j),tmn(i,j),
     &                     tsk(i,j)

               end if
            end do
         end do

1000     continue

      end program sheargeneric

!        subroutine write_grads(var,ix,jx,kx,length)
!        integer:: length,ix,jx,kx,k
!        real :: var(ix,jx,kx)
!
!        do k=1,kx
!         length=length+1
!         write(33,rec=length)var(:,:,K)
!        end do
!        return
!        end


