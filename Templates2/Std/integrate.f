


************************************************************************************* 
*         Calculate the derivates of the diffusion equation
*         This subroutine is called by the Runge-Kutta integrator
*************************************************************************************  

        subroutine derivs(time,fo,T,dT)
        implicit none
        include 'parEBM.h'

        real*8 time, fi1,fi2,dfi
        real*8 T(N),dT(N),flux(N)
        real*8 der(0:N),xk(0:N)
        real*8 x(N),out(N),alb(N)
        real*8 term1, term2 , tmed
        integer i

        real*8 Aterm, Iterm, Sterm, Dterm, Cterm
	real*8 fo(N), Ddry
	
	real*8 xxx
	 
	
	do i=1,N
	 x(i)= xmin+dx*(i-1)+dx/2.0    ! latitudine 
	 out(i)=  Iterm(i,T(i)) ! tabella OLR calcolata con i modelli ccm 
	end do 
	
	
	if(albedoFile.eq.'none') then  ! calculate albedo if albedoFile does not exist
	
        if(xmin<-1.5) then  ! latitude
           do i=1,N 
           alb(i)=Aterm(time,dsin(x(i)),T(i),i,fo)  
	   enddo
	else  ! sin(latitude)
	   do i=1,N 
           alb(i)=Aterm(time,x(i),T(i),i,fo)  
	   enddo
	endif
	
	else  ! read the albedo array if albedoFile exists 
	
	   open(unit=78,file=albedoFile,status='unknown')
              do i=1,N
                 read(78,*) xxx,alb(i) 
              enddo
           close(78)  
	
	end if  !  END ALBEDO

        der(0)=0.0 
        der(N)=0.0
        xk(0)=0.0
        xk(N)=0.0

        if(xmin<-1.5) then ! integration in latitude
           do i=1,N-1
              Tmed=(T(i+1)+T(i))/2
              der(i)=(T(i+1)-T(i))/dx/dcos(x(i)+dx/2)
          xk(i)=Dterm(time,x(i),Tmed,Ddry)*(1-dsin( x(i)+dx/2 )**2)
              flux(i)=(xk(i)*der(i)-xk(i-1)*der(i-1))/dx/
     >             dcos(x(i))
           end do
           flux(N)=-xk(N-1)*der(N-1)/dx/dcos(x(N))
        else
           do i=1,N-1
              Tmed=(T(i+1)+T(i))/2
              der(i)=(T(i+1)-T(i))/dx
          xk(i)=Dterm(time,x(i),Tmed,Ddry)*(1-(x(i)+dx/2)**2)
              flux(i)=(xk(i)*der(i)-xk(i-1)*der(i-1))/dx
           end do
           flux(N)=-xk(N-1)*der(N-1)/dx
        endif
       

        do i=1,N
           term2=Sterm(time,x(i))*(1-alb(i))  
           dT(i)=(flux(i)-out(i)+term2)/Cterm(T(i),time,i,fo)
        enddo

        if(idebug.eq.1) then 
           print *, ' '
           print *, 't=',time, N
           do i=1,N
              print *, i, x(i),T(i), dT(i), Sterm(time,x(i)), out(i)
              write(20,*), time, i, x(i),T(i), dT(i), 
     >                     Sterm(time,x(i)), out(i)
           enddo
        endif

        return
        end



		




*********************************************************
*       RUNGE-KUTTA INTEGRATOR and related routines
*********************************************************
      SUBROUTINE odeint(ystart,fo,nvar,x1,x2,eps,h1,hmin,nok,nbad)
      INTEGER nbad,nok,nvar,KMAXX,MAXSTP,NMAX
      REAL*8 eps,h1,hmin,x1,x2,ystart(nvar),TINY
*      EXTERNAL derivs,rkqs
      PARAMETER (MAXSTP=100000,NMAX=500,KMAXX=2000,TINY=1.e-30)
      INTEGER i,kmax,kount,nstp
      REAL*8 dxsav,h,hdid,hnext,x,xsav,dydx(NMAX),xp(KMAXX),y(NMAX),
     *yp(NMAX,KMAXX),yscal(NMAX)
      real*8 fo(nvar)
      COMMON /path/ kmax,kount,dxsav,xp,yp
      x=x1
      h=sign(h1,x2-x1)
      nok=0
      nbad=0
      kount=0
      do 11 i=1,nvar
        y(i)=ystart(i)
11    continue
      if (kmax.gt.0) xsav=x-2.*dxsav
      do 16 nstp=1,MAXSTP 
        call derivs(x,fo,y,dydx) 
        do 12 i=1,nvar
          yscal(i)=abs(y(i))+abs(h*dydx(i))+TINY
12      continue
        if(kmax.gt.0)then
          if(abs(x-xsav).gt.abs(dxsav)) then
            if(kount.lt.kmax-1)then
              kount=kount+1
              xp(kount)=x
              do 13 i=1,nvar
                yp(i,kount)=y(i)
13            continue
              xsav=x
            endif
          endif
        endif
        if((x+h-x2)*(x+h-x1).gt.0.) h=x2-x
        call rkqs(y,fo,dydx,nvar,x,h,eps,yscal,hdid,hnext)
        if(hdid.eq.h)then
          nok=nok+1
        else
          nbad=nbad+1
        endif
        if((x-x2)*(x2-x1).ge.0.)then
          do 14 i=1,nvar
            ystart(i)=y(i)
14        continue
          if(kmax.ne.0)then
            kount=kount+1
            xp(kount)=x
            do 15 i=1,nvar
              yp(i,kount)=y(i)
15          continue
          endif
          return
        endif
        if(abs(hnext).lt.hmin) then
        print *,'stepsize smaller than minimum in odeint',hnext 
        stop
        endif
        h=hnext
16    continue
      print *, 'too many steps in odeint'
      return
      END


      SUBROUTINE rkqs(y,fo,dydx,n,x,htry,eps,yscal,hdid,hnext)
      INTEGER n,NMAX
      REAL*8 eps,hdid,hnext,htry,x,dydx(n),y(n),yscal(n)
*      EXTERNAL derivs
      PARAMETER (NMAX=500)
CU    USES derivs,rkck
      INTEGER i
      REAL*8 errmax,h,xnew,yerr(NMAX),ytemp(NMAX),SAFETY,PGROW,PSHRNK,
     *ERRCON
      PARAMETER (SAFETY=0.9,PGROW=-.2,PSHRNK=-.25,ERRCON=1.89e-4)
      real*8 fo(n)
      h=htry
1     call rkck(y,fo,dydx,n,x,h,ytemp,yerr)
      errmax=0.
      do 11 i=1,n
        errmax=max(errmax,abs(yerr(i)/yscal(i)))
11    continue
      errmax=errmax/eps
      if(errmax.gt.1.)then
        h=SAFETY*h*(errmax**PSHRNK)
        if(h.lt.0.1*h)then
          h=.1*h
        endif
        xnew=x+h
        if(xnew.eq.x) print *, 'stepsize underflow in rkqs'
        goto 1
      else
        if(errmax.gt.ERRCON)then
          hnext=SAFETY*h*(errmax**PGROW)
        else
          hnext=5.*h
        endif
        hdid=h
        x=x+h
        do 12 i=1,n
          y(i)=ytemp(i)
12      continue
        return
      endif
      END


      SUBROUTINE rkck(y,fo,dydx,n,x,h,yout,yerr)
      INTEGER n,NMAX
      REAL*8 h,x,dydx(n),y(n),yerr(n),yout(n)
*      EXTERNAL derivs
      PARAMETER (NMAX=500)
CU    USES derivs
      INTEGER i
      REAL*8 ak2(NMAX),ak3(NMAX),ak4(NMAX),ak5(NMAX),ak6(NMAX),
     *ytemp(NMAX),A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,B51,B52,B53,
     *B54,B61,B62,B63,B64,B65,C1,C3,C4,C6,DC1,DC3,DC4,DC5,DC6
      PARAMETER (A2=.2,A3=.3,A4=.6,A5=1.,A6=.875,B21=.2,B31=3./40.,
     *B32=9./40.,B41=.3,B42=-.9,B43=1.2,B51=-11./54.,B52=2.5,
     *B53=-70./27.,B54=35./27.,B61=1631./55296.,B62=175./512.,
     *B63=575./13824.,B64=44275./110592.,B65=253./4096.,C1=37./378.,
     *C3=250./621.,C4=125./594.,C6=512./1771.,DC1=C1-2825./27648.,
     *DC3=C3-18575./48384.,DC4=C4-13525./55296.,DC5=-277./14336.,
     *DC6=C6-.25)
      real*8 fo(n)
      do 11 i=1,n
        ytemp(i)=y(i)+B21*h*dydx(i)
11    continue
      call derivs(x+A2*h,fo,ytemp,ak2)
      do 12 i=1,n
        ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))
12    continue
      call derivs(x+A3*h,fo,ytemp,ak3)
      do 13 i=1,n
        ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
13    continue
      call derivs(x+A4*h,fo,ytemp,ak4)
      do 14 i=1,n
        ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+B54*ak4(i))
14    continue
      call derivs(x+A5*h,fo,ytemp,ak5)
      do 15 i=1,n
        ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+B64*ak4(i)+
     *B65*ak5(i))
15    continue
      call derivs(x+A6*h,fo,ytemp,ak6)
      do 16 i=1,n
        yout(i)=y(i)+h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+C6*ak6(i))
16    continue
      do 17 i=1,n
        yerr(i)=h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)+DC6*
     *ak6(i))
17    continue
      return
      END


