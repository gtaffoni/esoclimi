************************************************************
*   this updates auxiliary seasonal vector and matrix                 *
************************************************************ 
      subroutine update_seasonal_qties_aux(is, veg, Vegmat)
     

      implicit none
      integer is, j
      include 'parEBM.h'
      real*8 Vegmat(Ns,N), veg(N)

      do j=1,N
         Vegmat(is,j)=veg(j)    ! update temperature matrix used in f_ice
      end do 
	            

      return
      end




************************************************************
*   this writes on file needed auxiliary seasonal data              *
************************************************************
      subroutine write_seasons_aux(i, is, nlastorbits, t2, veg, fo, phi, 
     >     phiDEG)
      
      implicit none
      integer i, is, nlastorbits
      include 'parEBM.h'
      real*8 veg(N), fo(N), t2
      real*8 phi(N) ! latitude in radians 
      real*8 phiDEG(N) ! latitude in degrees 

      integer j

      real*8 Aterm, f_ice

      
c     store season-latitude-vegetation of the last 5 orbits after convergence
      if( nlastorbits .ge.1 .and. nlastorbits .le. 5) then   
         do j=1,N
            write(24,30) (i+is/dfloat(Ns)),phiDEG(j),veg(j)*(1.-fo(j)), 
     >           veg(j) 
         enddo
      endif
      
c     store season-latitude-temperature of the last orbit of the simulations (5th orbit after convergence)
      if( nlastorbits .eq. 5 ) then 
         do j=1,N 
            write(28,30) (i+is/dfloat(Ns)),phiDEG(j),veg(j)*(1.-fo(j)),
     >           veg(j)
         enddo
      endif 
30      format(f14.8,1x,f14.8,1x,f14.8,1x,f14.8)
      
      
      return
      end
      


************************************************************
*   this calculates and updates mean annual quantities     *
*       for the current orbit                              *
************************************************************
      subroutine calculate_annual_means_aux(i, phi, fo,
     >     vegmat, annualglobalVegZonal, annualglobalVeg)


      implicit none
      include 'parEBM.h'
      include 'functions.h'

      
      integer i, is, j
      real*8 veg(N), vegmat(Ns,N), fo(N)
      real*8 phi(N)          ! latitude in degrees

      real*8 annualglobalVegZonal, annualglobalVeg
      real*8 GM,wt,sw,nw,sumlatcos, GMZ, sumcoslat


      integer k, ks, kl

c       sum of the weights of planet area
       sumcoslat=0.
       do j=1,N
          sumcoslat=sumcoslat+dcos(phi(j)) 
       end do 

c        area-weighted global mean for vegetation
       GM=0.
       GMZ=0.
       sw=0.
       nw=0.
       do j=1,N 
          do k=1,Ns
             wt=dcos(phi(j))
             GMZ=GMZ+vegmat(k,j)*(1.-fo(j))*wt
             GM=GM+vegmat(k,j)*(1.-fo(j))
             sw=sw+wt
             nw=nw+1.
          end do
       end do
       annualglobalVegZonal=GMZ/sw
       annualglobalVeg=GM/nw


      return
      end
          


************************************************************
*   this writes on file needed auxiliary seasonal data              *
************************************************************
      subroutine out_annual_veg_aux(i, nprompt, t2, annualglobalVeg, 
     >     annualglobalVegZonal)
      implicit none
      include 'parEBM.h'

      integer i, nprompt
      real*8 t2, annualglobalVeg, annualglobalVegZonal

      if( (i/nprompt)*nprompt .eq. i) then  
         write(*,7022) annualglobalVeg,annualglobalVegZonal
 7022     format('         ','<Veg>=',f9.4,' <Veg, Zonal>',f9.4)
      endif


      write(30, 7023) t2/Porb, annualglobalVeg, annualglobalVegZonal
 7023 format(f16.8, 1x, f16.8,1x,f16.8)

      return
      end

************************************************************
*   final output and closing of files                     *
************************************************************
      subroutine make_output_aux()
      implicit none
      include 'parEBM.h'

* no final output for the moment - think about it      
      close(24)
      close(28)
      close(30)

      
      return
      end
      
