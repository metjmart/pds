!     This subroutine makes a balanced shear flow state.

      subroutine makeShear (uofyz,pofyz,tofyz,qofyz,
     &                      tSound,qSound,rhSound,phyd,z,
     &                      usound,idealwinds,
     &                      ny,nz,dy,fofy,jvort,ispds,adjq)

         implicit none

         integer :: nz,ny,jvort,ispds,adjq,idealwinds
         real, dimension(ny,nz) :: uofyz, pofyz, tofyz, qofyz
         real, dimension(nz) :: tSound,qSound,rhSound,usound,phyd,z
         real, dimension(ny) :: fofy

         real :: dy

         integer :: i,j,k

         real :: errmax

         ! phyd is already computed

         if (idealwinds.eq.1) then
           call makeuofyz(uofyz,phyd,z,ny,nz,dy)
           print *, 'makeShear: done makeuofyz'
         else 
           do j=1,ny
             do k=1,nz

               uofyz(j,k)=usound(k)

             enddo 
           enddo
         endif

         call inittofyz(tofyz,tSound,z,ny,nz)
         print *, 'makeShear: done inittofyz'

         call initqofyz(qofyz,qSound,z,ny,nz)
         print *, 'makeShear: done initqofyz'

         if (ispds /= 1) then

            do i=1,5

               call makepofyz(pofyz,tofyz,qofyz,uofyz,
     &                        phyd,z,fofy,ny,nz,dy,jvort)
               print *, 'makeShear: done makepofyz'

               call updatetofyz(tofyz,qofyz,pofyz,z,ny,nz)
               print *, 'makeShear: done updatetofyz'

               if (adjq == 1) then
                  call updateqofyz(tofyz,qofyz,pofyz,rhSound,ny,nz)
                  print *, 'makeShear: done updateqofyz'
               end if

               call evalgraderryz(errmax,pofyz,uofyz,tofyz,qofyz,
     &                            z,ny,nz,dy,fofy)
               print *, 'makeShear: done evalyz'

               write (*, "('makeShear: on iteration ', i2,
     &                     ' balance errmax = ', 1pg14.7)")
     &               i, errmax

            end do

         else

!           For fake shear, just fill with central pressure profile
            do k=1,nz
               do j=1,ny
                  pofyz(j,k) = phyd(k)
               end do
            end do

         end if

         return
      end subroutine makeShear

