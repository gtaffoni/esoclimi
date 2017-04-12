
c        MEAN ANNUAL (ORBITAL) ZONAL MEAN
      subroutine ZonalMean(matrix,ZM)
      include 'parEBM.h' 
      real*8 matrix(Ns,N),ZM(N) 
      integer j,i
      do j=1,N
         ZM(j)=0.
         do i=1,Ns
            ZM(j)=ZM(j)+matrix(i,j)/Ns
         end do
      end do
      return
      end
      
c        MEAN ANNUAL (ORBITAL) AREA-WEIGHTED GLOBAL MEAN
      subroutine GlobalMean(matrix,phi,GM)
      include 'parEBM.h' 
      real*8 matrix(Ns,N),phi(N)
      real*8 wt,sw,GM 
      integer j,i
      GM=0.
      sw=0.
      do j=1,N 
         do i=1,Ns
            wt=dcos(phi(j))
            GM=GM+matrix(i,j)*wt
            sw=sw+wt
         end do
      end do
      GM=GM/sw
      return
      end

c        MEAN ANNUAL (ORBITAL) AREA-WEIGHTED GLOBAL MEAN
      subroutine HemisphereMean(matrix,phi,GM)
      include 'parEBM.h' 
      real*8 matrix(Ns,N),phi(N)
      real*8 wt,sw,GM 
      integer j,i
      GM=0.
      sw=0. 
      do j=1,N  
        wt=dcos(phi(j)) 
        if(phi(j).gt.0.) then
          do i=1,Ns 
          GM=GM+matrix(i,j)*wt
          sw=sw+wt
          end do
        end if
        if(phi(j).eq.0.) then
          GM=GM+matrix(i,j)*0.5*wt
          sw=sw+0.5*wt
        end if  
      end do 
      GM=GM/sw
      return
      end


      real*8 function dlat(j)
      integer j
      include 'parEBM.h' 
      if(xmin < -1.5) then 
           fi1 = dsin(xmin+dx*(j-1)) ! integration in latitude
           fi2 = dsin(xmin+dx*j)
        else
           fi1 = xmin+dx*(j-1)       ! integration in sin(latitude)
           fi2 = xmin+dx*j
      endif 
      dlat=fi2-fi1  
      return
      end





		

***********************************************************************
*     Define zonal ocean fractions 
*     The routine also returns the mean global ocean fraction fom
***********************************************************************
      subroutine geography(x,fo,fom)

      include 'parEBM.h'  ! the parameter gg is read here

      real*8 xo(18),foo(18),fo1(18),fo2(18),fo3(18)
      real*8 fo(N), x(N),fom 
      integer i,j,k
      real*8 x1,x2,dx1,dx2,dx3
      
 
      data xo/-85.,-75.,-65.,-55.,-45.,-35.,
     >        -25.,-15.,-5.,5.,15.,25.,
     >         35.,45.,55.,65.,75.,85./

c     Present Earth geography from WK97
      data fo1/0.001, 0.246, 0.896, 0.992, 0.970, 0.888,
     >         0.769, 0.780, 0.764, 0.772, 0.736, 0.624,
     >         0.572, 0.475, 0.428, 0.294, 0.713, 0.934/

c     Equatorial geography from WK97  
      data fo2/1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
     >         1.000, 0.250, 0.000, 0.000, 0.250, 1.000,
     >         1.000, 1.000, 1.000, 1.000, 1.000, 1.000/

c     Polar geography from WK97
      data fo3/0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     >         0.367, 1.000, 1.000, 1.000, 1.000, 1.000,
     >         1.000, 1.000, 1.000, 1.000, 1.000, 1.000/

      do i=1,18
      if(gg.eq.1) foo(i)=fo1(i) ! present Earth
      if(gg.eq.2) foo(i)=fo2(i) ! equatorial geography
      if(gg.eq.3) foo(i)=fo3(i) ! polar geography
      end do


      do i=1, 18
         xo(i) = ( (xo(i)-5) * pi2/360.)
*         foo(i) = foo(i) * dx /
*     >        (dsin(xo(i)+ 5.*pi2/360) -  dsin(xo(i)- 5.*pi2/360))/2. !fattore d'area
      end do 
      
      
      if(xmin < -1.5) then    ! caso x=phi   
 
      do i=1,N

         xx = x(i)
         if(xx.lt.xo(1)) then
            fo(i) = foo(1)
         else if(xx.ge.xo(18)) then
            fo(i) = foo(18)
         else
            k=-1
            do j=1,18
               if(xo(j).gt. xx .and. k.eq.-1) then
                  k=j-1
               endif
            end do
            
            x1 = x(i)-dx/2
            x2 = x(i)+dx/2
            dx1=0.0
            dx2=0.0
            dx3=0.0
           
            if(x1.ge.xo(k) .and. x2.le.xo(k+1)) then
               fo(i) = foo(k)
            else if(x1.lt.xo(k) .and. x2.le.xo(k+1)) then
               dx1 =  (xo(k) - x1)/2
               dx2 =  (dx2 - xo(k))/2
               if(k.gt.1) then
                  fo(i) = (foo(k-1)*dx1 + foo(k)*dx2)/(dx1+dx2)
               else
                  fo(i) = foo(k)
               endif
            else if(x1.ge.xo(k) .and. x2.gt.xo(k+1)) then
               dx1 =  (xo(k+1) - x1)/2
               dx2 =  (x2 - xo(k+1))/2
               fo(i) = (foo(k)*dx1 + foo(k+1)*dx2)/(dx1+dx2)
            else
               dx1 = (xo(k) - x1)/2
               dx2 = (x2 - xo(k+1))/2
               dx3 = (xo(k+1) -xo(k))/2
               if(k.ge.1) then
                  fo(i) = (foo(k-1)*dx1 + foo(k+1)*dx3 + 
     >                 foo(k)*dx2)/(dx1+dx2+dx3)
               else
                  fo(i) = (foo(k)*dx1 + foo(k+1)*dx3 + 
     >                 foo(k)*dx2)/(dx1+dx2+dx3)
               endif
            endif
            
*            print 111,i,k,x1,xx,x2,xo(k),xo(k+1), dx1, dx2, dx3
* 111        format(2(1x,i4),8(1x,f12.6))
         endif

      end do 

      else   !  caso x=sin(phi) 
 
      do i=1,N

         xx = dasin(x(i))
         if(xx.lt.xo(1)) then
            fo(i) = foo(1)
         else if(xx.ge.xo(18)) then
            fo(i) = foo(18)
         else
            k=-1
            do j=1,18
               if(xo(j).gt. xx .and. k.eq.-1) then
                  k=j-1
               endif
            end do
            
            x1 = dasin(x(i)-dx/2)
            x2 = dasin(x(i)+dx/2)
            dx1=0.0
            dx2=0.0
            dx3=0.0
           
            if(x1.ge.xo(k) .and. x2.le.xo(k+1)) then
               fo(i) = foo(k)
            else if(x1.lt.xo(k) .and. x2.le.xo(k+1)) then
               dx1 =  (dsin(xo(k)) - dsin(x1))/2
               dx2 =  (dsin(x2) - dsin(xo(k)))/2
               if(k.gt.1) then
                  fo(i) = (foo(k-1)*dx1 + foo(k)*dx2)/(dx1+dx2)
               else
                  fo(i) = foo(k)
               endif
            else if(x1.ge.xo(k) .and. x2.gt.xo(k+1)) then
               dx1 =  (dsin(xo(k+1)) - dsin(x1))/2
               dx2 =  (dsin(x2) - dsin(xo(k+1)))/2
               fo(i) = (foo(k)*dx1 + foo(k+1)*dx2)/(dx1+dx2)
            else
               dx1 =  (dsin(xo(k)) - dsin(x1))/2
               dx2 =  (dsin(x2) - dsin(xo(k+1)))/2
               dx3 = (dsin(xo(k+1)) -dsin(xo(k)))/2
               if(k.ge.1) then
                  fo(i) = (foo(k-1)*dx1 + foo(k+1)*dx3 + 
     >                 foo(k)*dx2)/(dx1+dx2+dx3)
               else
                  fo(i) = (foo(k)*dx1 + foo(k+1)*dx3 + 
     >                 foo(k)*dx2)/(dx1+dx2+dx3)
               endif
            endif
            
*            print 111,i,k,x1,xx,x2,xo(k),xo(k+1), dx1, dx2, dx3
* 111        format(2(1x,i4),8(1x,f12.6))
         endif

      end do 

      end if ! fine caso x=sin(phi)


      fom=0.
      do i=1,N
               if(xmin < -1.5) then 
                  fi1 = dsin(xmin+dx*(i-1)) !integration in phi
                  fi2 = dsin(xmin+dx*i)
               else
                  fi1 = xmin+dx*(i-1) !integration in sin(phi)
                  fi2 = xmin+dx*i
               end if
      dfi=fi2-fi1
      fom=fom+ dfi*fo(i)/2.0
      enddo

      
      return
      end 




************************************************************
*  calculate Liquid Water Temperature Range
************************************************************
      subroutine LWTR(pressTOT,Tice,Tvap)
      real*8 pressTOT,Tice,Tvap,P  
      include 'parEBM.h'  
      integer i,NW1,NW2
      parameter(NW1=4)
      parameter(NW2=14)

c     data from CRC Handbook of Chemistry and Physics
      
      real*8 P1(NW1),T1(NW1),P2(NW2),T2(NW2)
      
      data P1/5.,10.,1000.,10000./  ! in mb       
      data T1/273.15,273.15,273.15,273.09/
      
      data P2/5.,10.,20.,30.,50.,100.,300.,1000.,2000.,
     > 3000.,5000.,7000.,8300.,10000./ ! in mb
      data T2/273.15,273.15,288.,299.55,306.03,318.97,342.26,
     > 372.78,393.39,  ! primo valore estrapolato
     > 407.,425.,438.,445.,454./ ! questi ultimi presi da engineeringtoolbox.com (non si trovano nel Handbook)
     
      P=pressTOT*1.e-2 ! convert Pa to mb 
      
      do i=1,NW1-1
      if(P.ge.P1(i).and.P.lt.P1(i+1)) then
      Tice=T1(i)+
     > (T1(i+1)-T1(i))*(P-P1(i))/(P1(i+1)-P1(i)) 
      end if
      end do 
      
      if(P.gt.P2(NW2)) then
      Tvap=T2(NW2)
      print *,'Extrapolated Tvap=',Tvap
      return
      end if      
      
      do i=1,NW2-1
      if(P.ge.P2(i).and.P.lt.P2(i+1)) then 
      Tvap=T2(i)+
     > (T2(i+1)-T2(i))*(P-P2(i))/(P2(i+1)-P2(i))
      end if
      end do  
      
      return
      end



****************************************************************
c          calculate saturated H2O pressure at temperature T (K)
****************************************************************
        real*8 function P_H2O(T)
	real*8 T 

	P_H2O=dexp(77.3450+0.0057*T-7235./T)/T**8.2
c       function found at
c       http://www.engineeringtoolbox.com/water-vapor-saturation-pressure-air-d_689.html
c       gives result in Pa, T must be in K
c       checked with CRC Handbook of Chemistry 
	
c	P_H2O=P_H2O/1.d5
c       conversion from Pa to Bar	
	
	return
	end



!**************************************************************
! Quadrilinear interpolator used for Top-of-Atmosphere albedo
!**************************************************************	  
      real*8 function ToaQuadLinInterp(T, p, Z, as)

      real*8 T, p, Z, as 
      integer iT, ip, iz, ias

      include 'parEBM.h' 

      real*8 num1111, num1112, num1121, num1211
      real*8 num2111, num2211, num2221, num2222
      real*8 num2212, num1222, num1212, num1221
      real*8 num2122, num2112, num2121, num1122, den

      common /interp_TOA/ Matrix_TOA,T_TOA,p_TOA,z_TOA,as_TOA 

      call hunt(T_TOA,NT_TOA,T,iT)
      call hunt(p_TOA,Np_TOA,p,ip)
      call hunt(z_TOA,Nz_TOA,Z,iz)
      call hunt(as_TOA,Nas_TOA,as,ias)

      x1=T_TOA(iT)
      x2=T_TOA(iT+1)
      y1=p_TOA(ip)
      y2=p_TOA(ip+1)
      z1=z_TOA(iz)
      z2=z_TOA(iz+1)
      w1=as_TOA(ias)
      w2=as_TOA(ias+1)

      dx1 = (x2-T)
      dy1 = (y2-p)
      dz1 = (z2-Z)
      dw1 = (w2-as)

      dx2 = (T-x1)
      dy2 = (p-y1)
      dz2 = (Z-z1)
      dw2 = (as-w1)
 
      den = (x2-x1)*(y2-y1)*(z2-z1)*(w2-w1)

      num1111 = Matrix_TOA(iT,ip,iz,ias) * dx1 * dy1 * dz1 * dw1
      num2222 = Matrix_TOA(iT+1,ip+1,iz+1,ias+1) * dx2 * dy2 * dz2 *dw2
      num2111 = Matrix_TOA(iT+1,ip,iz,ias) * dx2 * dy1 * dz1 * dw1
      num1211 = Matrix_TOA(iT,ip+1,iz,ias) * dx1 * dy2 * dz1 * dw1
      num1121 = Matrix_TOA(iT,ip,iz+1,ias) * dx1 * dy1 * dz2 * dw1
      num1112 = Matrix_TOA(iT,ip,iz,ias+1) * dx1 * dy1 * dz1 * dw2

      num2211 = Matrix_TOA(iT+1,ip+1,iz,ias) * dx2 * dy2 * dz1 * dw1
      num2221 = Matrix_TOA(iT+1,ip+1,iz+1,ias) * dx2 * dy2 * dz2 * dw1
     
      num1221 = Matrix_TOA(iT,ip+1,iz+1,ias) * dx1 * dy2 * dz2 * dw1
      num1222 = Matrix_TOA(iT,ip+1,iz+1,ias+1) * dx1 * dy2 * dz2 * dw2

      num1122 = Matrix_TOA(iT,ip,iz+1,ias+1) * dx1 * dy1 * dz2 * dw2

      num2122 = Matrix_TOA(iT+1,ip,iz+1,ias+1) * dx2 * dy1 * dz2 * dw2
      num2112 = Matrix_TOA(iT+1,ip,iz,ias+1) * dx2 * dy1 * dz1 * dw2
      num2212 = Matrix_TOA(iT+1,ip+1,iz,ias+1) * dx2 * dy2 * dz1 * dw2

      num1212 = Matrix_TOA(iT,ip+1,iz,ias+1) * dx1 * dy2 * dz1 * dw2
      num2121 = Matrix_TOA(iT+1,ip,iz+1,ias) * dx2 * dy1 * dz2 * dw1

      
      ToaQuadLinInterp = (num1111 + num2222 + num2111 + num1211 +
     >     num1121 + num1112 + num2211 + num2221 + num1221 + num1222 +
     >     num1122 + num2122 + num2112 + num2212 + num1212 + num2121) / 
     >     den

* per la formula si veda ad es 
* http://en.wikipedia.org/wiki/Bilinear_interpolation 

      return
      end

!**************************************************************
! Used by the quadrilinear interpolator in TOA albedo
!**************************************************************	   
      SUBROUTINE HUNT(XX,N,X,JLO)
      
      real*8 XX(N)
      !DIMENSION XX(N)
      real*8 X
      LOGICAL ASCND
      
      ASCND=XX(N).GT.XX(1)
      IF(JLO.LE.0.OR.JLO.GT.N)THEN
        JLO=0
        JHI=N+1
        GO TO 3
      ENDIF
      INC=1
      IF(X.GE.XX(JLO).EQV.ASCND)THEN
1       JHI=JLO+INC
        IF(JHI.GT.N)THEN
          JHI=N+1
        ELSE IF(X.GE.XX(JHI).EQV.ASCND)THEN
          JLO=JHI
          INC=INC+INC
          GO TO 1
        ENDIF
      ELSE
        JHI=JLO
2       JLO=JHI-INC
        IF(JLO.LT.1)THEN
          JLO=0
        ELSE IF(X.LT.XX(JLO).EQV.ASCND)THEN
          JHI=JLO
          INC=INC+INC
          GO TO 2
        ENDIF
      ENDIF
3     IF(JHI-JLO.EQ.1)RETURN
      JM=(JHI+JLO)/2
      IF(X.GT.XX(JM).EQV.ASCND)THEN
        JLO=JM
      ELSE
        JHI=JM
      ENDIF
      GO TO 3
      END

!***********************************************************


        real*8 function DiurnalMeanMu(time,x,HourD)
        implicit none
        include 'parEBM.h'
        real*8 time
	real*8 phi, x, xx
	real*8 nu, xl, Ls
	real*8 sind, delta, HourD, cosH, HH

        if(xmin<-1.5) then ! integration in latitude
           phi = x
           xx = dsin(x)
        else
           xx = x               !integration in sin(latitude)
           phi =  dasin(xx)
        endif 

        if(xx.lt.-1.0) xx=-1.0
        if(xx.gt. 1.0) xx= 1.0

        Ls=nu(time)-nu0 ! nu(0.)
        xl   =  dmod(Ls+Ls0,pi2)
        sind= -dsin(OBLIQUITY)*dcos(xl)
        delta=dasin(sind)

	cosH = -dtan(phi)*dtan(delta) !!il problema e' QUI!!!!! 
        if(cosH.gt.1.0d0) cosH =  1.0d0;
        if(cosH.lt.-1.0d0) cosH = -1.0d0;

	HH=dacos(cosH)
        
	if(dabs(HH).gt.0.) then
        DiurnalMeanMu=dsin(phi)*dsin(delta) +
     >          dcos(phi)*dcos(delta)*dsin(HH)/HH 
        else
	DiurnalMeanMu=0.
	end if 
        
	HourD=HH
        return
        end
     
