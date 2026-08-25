
        subroutine makevtofrz(vtofrz,vtbar,z,nr,nz,dr)

        implicit none

        integer, intent(in) :: nr,nz
        real, dimension(nr,nz), intent(out) :: vtofrz
        real, dimension(nr), intent(out) :: vtbar
        real, dimension(nz), intent(in) :: z
        real, intent(in) :: dr

        integer :: i,j
        real :: ztop,zhere,fac,r,Lz
        real :: alpha,zmax,Lzup,Lzdown,rcut
        
        Lzdown = 2500.0     ! length scale of decay downward
        Lzup = 3500         ! length scale of decay upward
        alpha = 1.8         ! exponent for Gaussian-like decay
        zmax = 1500.0       ! level of max winds

        rcut = 600000.0     ! radius that winds are forced to zero


        call makevtbar(vtbar,nr,dr)

        do j=1,nz

          zhere=z(j)

          do i=1,nr
         
            r=i*dr

            if (zhere.gt.zmax) then
              Lz=Lzup
            else
              Lz=Lzdown
            endif

            fac = exp( -(abs(zhere-zmax)**alpha)/(alpha * Lz**alpha) )

            fac=fac*exp( -(r/rcut)**4.0 )

            vtofrz(i,j) = fac*vtbar(i)

          end do
        end do  

        return
        end          



