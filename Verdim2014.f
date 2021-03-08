      program verdim 
      common e,fyd,fcd,np,xp,yp,xg,yg,lx,ly,area,jox,joy,joxy, 
     * na,max,may,nb,xb,yb,perc,nrc,il(5),gs,gc,as 
     * ,epsc2,epscu,a1,a2,const 
      real*8 epsc2(5),epscu(5),a1(5),a2(5),const,fck 
      real*8 nn,y,x,x2,x3,x4,xy,x2y,ddx 
      real*8 e,fyd(100),fcd(5),xp(100),yp(100),xg,yg,area,xb(100), 
     * yb(100),perc(100) 
      real*8 jox,joy,joxy,lx,ly,na,max,may 
      real*8 jx,jy,jxy,as,dx,dy,gc,gs,sox,soy,sx,sy,xmax,xmin,ymax,ymin 
      character arq*12 
      integer op 
   10 write(*,11) 
   11 format(////,'>>>>>> qual o nome do arquivo de dados ?') 
      read(*,'(a12)',err=10) arq 
      if(arq.eq.' ') stop 
      ir=1 
      iw=0 
      open(unit=1,file=arq,status='old',form='formatted', 
     *access='sequential') 
      read(ir,21) op 
      read(ir,20) gc,gs 
      read(ir,21) np 
      do 30 j1=1,np 
         read(ir,20) xp(j1),yp(j1) 
   30 continue 
      read(ir,21) nrc 
      do 40 j1=1,nrc 
         read(ir,24) fcd(j1),il(j1) 
         fcd(j1)=fcd(j1)/gc 
   40 continue 
      read(ir,22) nb,e,as 
      do 50 j1=1,nb 
         read(ir,23) xb(j1),yb(j1),fyd(j1),perc(j1) 
         fyd(j1)=fyd(j1)/gs 
   50 continue 
      read(ir,20) na,max,may 
      const=210000.d0/e 
      do 55 j1=1,nrc 
          fck = fcd(j1)*gc*const 
          if(fck.le.50) then 
              epsc2(j1)=0.002d0 
              epscu(j1)=0.0035d0 
              a1(j1)=1000.d0 
              a2(j1)=250000.d0 
          else 
              epsc2(j1)=0.002d0+0.000085d0*(fck-50.d0)**0.53 
              epscu(j1)=0.0026d0+0.035d0*((90-fck)/100)**4 
              nn=1.4d0+23.4d0*((90-fck)/100)**4 
              ddx=epsc2(j1)/1000.d0 
              x=0.d0 
              x2=0.d0 
              x3=0.d0 
              x4=0.d0 
              xy=0.d0 
              x2y=0.d0 
              do 51 j2=1,1000 
                  y=1.d0-(1.d0-x/epsc2(j1))**nn 
                  x2=x2+x*x 
                  x3=x3+x*x*x 
                  x4=x4+x*x*x*x 
                  xy=xy+x*y 
                  x2y=x2y+x*x*y 
                  x=x+ddx 
   51         continue               
              a1(j1)=(x2y-x4*xy/x3)/(x3-x2*x4/x3) 
              a2(j1)=-(x2y-x3*a1(j1))/x4 
          end if 
   55 continue              
      xmax=-1.e11 
      ymax=-1.e11 
      xmin=1e11 
      ymin=1e11 
      do 60 j1=1,np 
         if(xp(j1).gt.xmax) xmax=xp(j1) 
         if(xp(j1).lt.xmin) xmin=xp(j1) 
         if(yp(j1).gt.ymax) ymax=yp(j1) 
         if(yp(j1).lt.ymin) ymin=yp(j1) 
   60 continue 
      lx=xmax-xmin 
      ly=ymax-ymin 
      np1=np-1 
      area=0 
      sx=0 
      sy=0 
      jx=0 
      jy=0 
      jxy=0 
      do 70 j1=1,np1 
         dx=xp(j1+1)-xp(j1) 
         dy=yp(j1+1)-yp(j1) 
         area=area+(xp(j1)+dx/2.)*dy 
         sx=sx+(xp(j1)*(yp(j1)+dy/2.)+dx*(yp(j1)/2.+dy/3.))*dy 
         sy=sy+(xp(j1)*(xp(j1)+dx)+dx*dx/3.)*dy/2. 
         jx=jx+(xp(j1)*(yp(j1)*(dy+yp(j1))+dy*dy/3.)+dx*(yp(j1)* 
     .    (yp(j1)/2.+dy/1.5)+dy*dy/4.))*dy 
         jy=jy+(dx**3/4.+xp(j1)*(dx*dx+xp(j1)*(1.5*dx+xp(j1))))*dy/3. 
         jxy=jxy+(xp(j1)*(xp(j1)*(yp(j1)+dy/2.)+dx*(yp(j1)+dy/1.5))+ 
     .    dx*dx*(yp(j1)/3.+dy/4.))*dy/2. 
   70 continue 
      xg=sy/area 
      yg=sx/area 
      sox=sx-yg*area 
      soy=sy-xg*area 
      jox=jx-area*yg*yg 
      joy=jy-area*xg*xg 
      joxy=jxy-xg*yg*area 
      do 80 j1=1,np 
         xp(j1)=xp(j1)-xg 
         yp(j1)=yp(j1)-yg 
   80 continue 
      do 90 j1=1,nb 
         xb(j1)=xb(j1)-xg 
         yb(j1)=yb(j1)-yg 
   90 continue 
      call ajustl(op) 
      close(1) 
      go to 10 
   20 format(8f10.0) 
   21 format(8i10) 
   22 format(i10,7f10.0) 
   23 format(4f10.0) 
   24 format(f10.0,9i10) 
      end 
 
      subroutine ajustl(op) 
      integer op 
      real*8 nr,mrx,mry,lam,lammin,nrmin,mrxmin,mrymin 
      real*8 r(3,2),rt(3,3),dp(3) 
      common e,fyd,fcd,np,xp,yp,xg,yg,lx,ly,area,jox,joy,joxy, 
     *   na,max,may,nb,xb,yb,perc,nrc,il(5),gs,gc,as 
     *   ,epsc2,epscu,a1,a2 
      real*8 epsc2(5),epscu(5),a1(5),a2(5),const 
      real*8 e,fyd(100),fcd(5),xp(100),yp(100),xg,yg,area,xb(100), 
     *   yb(100),perc(100) 
      real*8 jox,joy,joxy,lx,ly,na,max,may 
      real*8 alfa0,alph,alpg,as1,as2,as3,as,b,bas,c,ca,ca0,car,x, 
     *   epsi,epss,fs,gc,gs,pjx,pjy,sa,sa0,ss,tol,tole,pi,pi2,graus 
      data pi,pi2,graus/3.1415926535897932385,1.5707963267948966192, 
     *   57.29577951308232/ 
      tolmin=1.d0 
      k=0 
      iw=0 
      tole=1.d-8 
      bas=max*max+may*may+na*na 
      lam=1.d0 
      alfa0=0 
      if(jox.eq.joy.and.abs(joxy).gt.1e-5) then 
         alfa0=pi2 
      else 
         if(jox.ne.joy) alfa0=atan(-2*joxy/(jox-joy))/2. 
      endif 
      ca0=cos(alfa0) 
      sa0=sin(alfa0) 
      ca0=ca0*ca0 
      sa0=sa0*sa0 
      ss=joxy*sin(2*alfa0) 
      pjx=jox*ca0+joy*sa0-ss 
      pjy=joy*ca0+jox*sa0+ss 
      alph=0 
      if(max.eq.0) then 
         alph=-dsign(pi2,may) 
      else 
         alph=atan(may*pjx/(max*pjy)) 
         if(max.gt.0) alph=alph+pi 
      endif 
      alfa0=alph+alfa0 
      uu=0.5d0 
  300 uu=(pi+uu)**5 
      uu=uu-int(uu) 
      k0=0 
      x=(lx+ly)*uu 
      if(op.eq.1) then 
         as1=abs(max)/(0.4*ly*fyd(1)) 
         as2=abs(may)/(0.4*lx*fyd(1)) 
         if(na.gt.0) then 
            as3=na/fyd(1) 
         else 
            as3=dmax1(0d0,(na-fcd(1)*area)/fyd(1)) 
         endif 
         as=as1+as2+as3 
      endif 
      alph=alfa0 
  890 call esfor(e,fyd,fcd,np,xp,yp,nb,xb,yb,perc,x,alph,as,b,c, 
     *epss,epsi,nr,mrx,mry,r,nrc,il,epsc2,epscu,a1,a2) 
      dp(1)=max-lam*mrx 
      dp(2)=may-lam*mry 
      dp(3)=na-lam*nr 
      tol=sqrt((dp(1)**2+dp(2)**2+dp(3)**2)/bas) 
      if(tol.le.tole) go to 900 
      k=k+1 
      k0=k0+1 
      if(k0.le.50) go to 301 
      if(tol.lt.tolmin) then 
        tolmin=tol 
        mrxmin=mrx 
        mrymin=mry 
        nrmin=nr 
        asmin=as 
        epssmin=epss 
        epsimin=epsi 
        alphmin=alph 
        lammin=lam 
      end if 
      go to 300 
  301 ca=cos(alph) 
      sa=sin(alph) 
      rt(1,1)=lam*(r(1,1)*ca-r(2,1)*sa) 
      rt(1,2)=lam*(-mry) 
      rt(2,1)=lam*(r(1,1)*sa+r(2,1)*ca) 
      rt(2,2)=lam*mrx 
      rt(3,1)=lam*r(3,1) 
      rt(3,2)=0 
      if(op.eq.2) then 
         rt(1,3)=mrx 
         rt(2,3)=mry 
         rt(3,3)=nr 
      else 
         rt(1,3)=r(1,2)*ca-r(2,2)*sa 
         rt(2,3)=r(1,2)*sa+r(2,2)*ca 
         rt(3,3)=r(3,2) 
      endif 
      call pivo(rt,dp,iver) 
      if(iver.eq.1) go to 300 
      x=x+dp(1) 
      if(op.eq.2) then 
         lam=lam+dp(3) 
      else 
         as1=as+dp(3) 
               as=as1 
      endif 
      alph=alph+dp(2) 
      if(abs(alph).gt.2*pi) alph=sign(mod(alph,2*pi),alph) 
      if(k.lt.10000) go to 890 
      write(*,*) ">>>> nao convergiu",tolmin 
      mrx=mrxmin 
      mry=mrymin 
      nr=nrmin 
      as=asmin 
      epss=epssmin 
      epsi=epsimin 
      alph=alphmin 
      lam=lammin 
  900 write(iw,500) 
      write(iw,500) 
      write(iw,501) 
      if(op.eq.2) then 
         write(iw,503) 
      else 
         write(iw,502) 
      endif 
      write(iw,501) 
      write(iw,500) 
      write(iw,501) 
      write(iw,504) 
      write(iw,505) max,may,na 
      write(iw,506) mrx,mry,nr 
      write(iw,501) 
      write(iw,500) 
      write(iw,501) 
      if(op.eq.2) then 
         fs=1./lam 
         write(iw,511) fs 
      else 
         write(iw,507) as 
      endif 
      write(iw,508) epss 
      write(iw,509) epsi 
      if(abs(alph).gt.2*pi) alph=sign(mod(alph,2*pi),alph) 
      alpg=graus*alph 
      write(iw,510) alpg 
      write(iw,501) 
	  write(iw,500) 
      write(iw,500) 
      write(*,512) 
      read(*,513) car 
      return 
  500 format(1x,78('*'))     
  501 format(1x,'**',74x,'**') 
  502 format(1x,'**     dimensionamento de secao de concreto armado a so 
     *licitacoes normais   **') 
  503 format(1x,'**      verificacao de secao de concreto armado a solic 
     *itacoes normais      **') 
  504 format(1x,'**                                           mx         
     *  my          n      **') 
  505 format(1x,'**    esforcos atuantes de calculo:   ',3e12.4,'  **') 
  506 format(1x,'**    esforcos resistentes de calculo:',3e12.4,'  **') 
  507 format(1x,'**',12x,'area total de armadura:               ', 
     *e15.4,9x,'**') 
  508 format(1x,'**',12x,'deformacao na fibra superior da secao:', 
     *e15.4,9x,'**') 
  509 format(1x,'**',12x,'deformacao na fibra inferior da secao:', 
     *e15.4,9x,'**') 
  510 format(1x,'**',12x,'inclinacao da linha neutra:           ', 
     *e15.4,9x,'**') 
  511 format(1x,'**',12x'reserva:                              ', 
     *e15.4,9x,'**') 
  512 format(/,25x,' tecle <enter> para continuar',/) 
  513 format(a1) 
      end 
 
      subroutine pivo(a,b,iver) 
      real*8 a(3,3),b(3) 
      integer ii(18),iver 
      real*8 aux1,aux2,dum1,dum2,dum3,dum4,dum5 
      data ii/1,2,3,1,3,2,2,1,3,2,3,1,3,1,2,3,2,1/ 
      iver=0 
      do 10 j1=1,6 
         jj=3*j1 
         i=ii(jj-2) 
         j=ii(jj-1) 
         k=ii(jj) 
         if(a(i,i).eq.0) go to 10 
         aux1=a(j,i)/a(i,i) 
         aux2=a(k,i)/a(i,i) 
         dum1=a(k,k)-a(i,k)*aux2 
         dum2=a(k,j)-a(i,j)*aux2 
         dum3=b(k)-b(i)*aux2 
         aux2=a(j,j)-a(i,j)*aux1 
         if(aux2.eq.0) go to 10 
         dum4=(b(j)-b(i)*aux1)/aux2 
         dum5=(a(j,k)-a(i,k)*aux1)/aux2 
         aux1=dum1-dum2*dum5 
         if(aux1.eq.0) go to 10 
         b(k)=(dum3-dum2*dum4)/aux1 
         b(j)=dum4-dum5*b(3) 
         b(i)=(b(i)-b(j)*a(i,j)-b(k)*a(i,k))/a(i,i) 
         go to 20 
   10 continue 
      iver=1 
   20 return 
      end 
  

      subroutine esfor(e,fyd,fc,np,xp,yp,nb,xb,yb,perc,x,alfa,as,b,c, 
     *epss,epsi,nrz,mrx,mry,r,nrc,il,epsc2,epscu,a1,a2) 
      real*8 nrzti,nrzt,mrks,mret,nrz,mrx,mry,ks1i,ks2i,ks1ii,ks2ii 
      real*8 xp(*),yp(*),xb(*),yb(*),perc(*),ksp(100),etp(100), 
     *     ksb(100),etb(100),r(3,2),fyd(*),fc(*),epsp(100) 
      real*8 epsc2(*),epscu(*),a1(*),a2(*) 
      real*8 x,dd,alfa,as,b,blx,c,ca,clx,dum1,dum2,e,eps0,epsb,eps1, 
     *epsi,epss,et,et01,et12,et1i,et1ii,et2i,et2ii,fcd,sa,sig,tetia, 
     *tetic,tetsa,tetsc,tksia,tksic,tkssa,tkssc 
      integer il(*) 
      ca=cos(alfa) 
      sa=sin(alfa) 
      tetsc=0 
      tetic=0 
      do 10 j1=1,np 
         ksp(j1)=xp(j1)*ca+yp(j1)*sa 
         etp(j1)=-xp(j1)*sa+yp(j1)*ca 
         if(etp(j1).gt.tetsc) then 
            tetsc=etp(j1) 
            tkssc=ksp(j1) 
         endif 
         if(etp(j1).lt.tetic) then 
            tetic=etp(j1) 
            tksic=ksp(j1) 
         endif 
   10 continue 
      tetsa=0 
      tetia=0 
      do 20 j1=1,nb 
         ksb(j1)=xb(j1)*ca+yb(j1)*sa 
         etb(j1)=-xb(j1)*sa+yb(j1)*ca 
         if(etb(j1).gt.tetsa) then 
            tetsa=etb(j1) 
            tkssa=ksb(j1) 
         endif 
         if(etb(j1).lt.tetia) then 
            tetia=etb(j1) 
            tksia=ksb(j1) 
         endif 
   20 continue 
      h=tetsc-tetic 
      dd=tetsc-tetia 
      x23=epscu(1)/(0.01d0+epscu(1))*dd 
      if(x.lt.x23) then 
         b=-0.01d0/(dd-x) 
         c=b*(x-tetsc) 
         blx=-0.01d0/(dd-x)**2 
         clx=(dd-tetsc)*blx 
         epsi=0.01d0 
         epss=-0.01d0*x/(dd-x) 
      else 
         if(x.lt.h) then 
            b=-epscu(1)/x 
            c=-epscu(1)-b*tetsc 
            blx=epscu(1)/x**2 
            clx=-tetsc*blx 
            epss=-epscu(1) 
            epsi=dmax1(epscu(1)*(dd-x)/x,0.d0) 
         else 
            if(x.gt.1.d150) then 
               b=0 
	           c=-epsc2(1) 
               blx=1.d-100 
               clx=1.d-100 
               epss=-epsc2(1) 
               epsi=-epsc2(1) 
               go to 29 
            end if 
            b=-epsc2(1)/(x-(epscu(1)-epsc2(1))/epscu(1)*h) 
            c=b*(x-tetsc) 
            blx=epsc2(1)/(x-(epscu(1)-epsc2(1))/epscu(1)*h)**2 
            clx=(tetsc-(epscu(1)-epsc2(1))/epscu(1)*h)*blx 
            epss=-epsc2(j1)*x/(x-(epscu(1)-epsc2(1))/epscu(1)*h) 
            epsi=-epsc2(j1)*(x-h)/(x-(epscu(1)-epsc2(1))/epscu(1)*h) 
         end if 
      end if 
   29 continue    
      do 30 j1=1,np 
         epsp(j1)=b*etp(j1)+c 
   30 continue 
      nrzt=0 
      mrks=0 
      mret=0 
      do 50 j1=1,3 
         do 50 j2=1,2 
            r(j1,j2)=0 
   50 continue 
      do 60 j1=1,nb 
         epsb=b*etb(j1)+c 
         call aco(e,epsb,fyd(j1),sig,et) 
         dum1=perc(j1)*as*et*(blx*etb(j1)+clx) 
         dum2=perc(j1)*sig 
         nrzti=as*dum2 
         nrzt=nrzt+nrzti 
         mrks=mrks+nrzti*etb(j1) 
         mret=mret-nrzti*ksb(j1) 
         r(1,1)=r(1,1)+dum1*etb(j1) 
         r(1,2)=r(1,2)+dum2*etb(j1) 
         r(2,1)=r(2,1)-dum1*ksb(j1) 
         r(2,2)=r(2,2)-dum2*ksb(j1) 
         r(3,1)=r(3,1)+dum1 
         r(3,2)=r(3,2)+dum2 
   60 continue 
      if(abs(epss-epsi).le.1e-10) then 
         if(epss.ge.0) go to 100 
         do 70 j1=1,nrc 
            fcd=fc(j1) 
            if(j1.eq.1) then 
               np1=1 
            else 
               np1=il(j1-1) 
            endif 
            np2=il(j1)-1 
            call centra(np1,np2,fcd,b,c,epss,ksp,etp,nrzt,mrks,mret, 
     *         blx,clx,r,epsc2(j1),a1(j1),a2(j1)) 
   70    continue 
      else 
         if(epss.ge.0.and.epsi.ge.0) go to 100 
         et01=-c/b 
         do 80 j1=1,nrc 
            et12=(-epsc2(j1)-c)/b 
            fcd=fc(j1) 
            if(j1.eq.1) then 
               np1=1 
            else 
               np1=il(j1-1) 
            endif 
            np2=il(j1)-1 
            do 80 j2=np1,np2 
               eps0=epsp(j2) 
               eps1=epsp(j2+1) 
               if(eps0.eq.eps1) go to 80 
               if(eps0.ge.0.and.eps1.ge.0) go to 80 
               call difer(j2,et01,et12,ksp,etp,eps0,eps1,ks1i,et1i,ks2i, 
     *            et2i,ks1ii,et1ii,ks2ii,et2ii,epsc2(j1)) 
               call regi(fcd,b,c,ks1i,et1i,ks2i,et2i,nrzt,mrks,mret, 
     *            blx,clx,r,a1(j1),a2(j1)) 
               call regii(fcd,ks1ii,et1ii,ks2ii,et2ii,nrzt,mrks,mret) 
   80    continue 
      endif 
  100 nrz=nrzt 
      mrx=mrks*ca-mret*sa 
      mry=mrks*sa+mret*ca 
      return 
      end 
 
 
      subroutine aco(e,epsb,fyd,sig,et) 
      real*8 epsb,fyd 
      real*8 a,b,c,dum1,e,eps1,eps2,et,sig 
      eps2=fyd/e 
      if(abs(epsb).le.eps2) then 
         sig=e*epsb 
         et=e 
      else 
         sig=sign(fyd,epsb) 
         et=0 
      endif 
      return 
      end 
 
 
      subroutine difer(i,et01,et12,ksp,etp,eps0,eps1,ks1i,et1i,ks2i, 
     *   et2i,ks1ii,et1ii,ks2ii,et2ii,epsc2) 
      integer t01,t12 
      real*8 etp(*),ksp(*),ks1i,ks2i,ks1ii,ks2ii,ks01,ks12 
      real*8 det,det01,det12,dksdet,dum1,dum2,eps0,eps1,et01,et12,et1i, 
     * et1ii,et2i,et2ii,epsc2 
      t01=0 
      t12=0 
      ks1i=0 
      et1i=0 
      ks2i=0 
      et2i=0 
      ks1ii=0 
      et1ii=0 
      ks2ii=0 
      et2ii=0 
      i2=i+1 
      det=etp(i2)-etp(i) 
      dksdet=(ksp(i2)-ksp(i))/det 
      dum1=et01-etp(i) 
      dum2=et12-etp(i) 
	  ks01=ksp(i)+dum1*dksdet 
      ks12=ksp(i)+dum2*dksdet 
      det01=dum1/det 
      det12=dum2/det 
      if(det01.gt.0.and.det01.lt.1) t01=1 
      if(det12.gt.0.and.det12.lt.1) t12=1 
      if(eps0.lt.eps1) then 
         t01=-t01 
         t12=-t12 
      endif 
      if(t01.eq.0.and.t12.eq.0) then 
         if(eps0.lt.0) then 
            if(eps0.gt.-epsc2) then 
               ks1i=ksp(i) 
               et1i=etp(i) 
               ks2i=ksp(i2) 
               et2i=etp(i2) 
            else 
               ks1ii=ksp(i) 
               et1ii=etp(i) 
               ks2ii=ksp(i2) 
               et2ii=etp(i2) 
            endif 
         endif 
      else 
         if(t01.eq.1) then 
            ks1i=ks01 
            et1i=et01 
            if(t12.eq.1) then 
               ks2i=ks12 
               et2i=et12 
               ks1ii=ks12 
               et1ii=et12 
               ks2ii=ksp(i2) 
               et2ii=etp(i2) 
            else 
               ks2i=ksp(i2) 
               et2i=etp(i2) 
            endif 
         else 
            if(t01.eq.-1) then 
               ks2i=ks01 
               et2i=et01 
               if(t12.eq.-1) then 
                  ks1i=ks12 
                  et1i=et12 
                  ks2ii=ks12 
                  et2ii=et12 
                  ks1ii=ksp(i) 
                  et1ii=etp(i) 
               else 
                  ks1i=ksp(i) 
                  et1i=etp(i) 
               endif 
            else 
               if(t12.eq.1) then 
                  ks1i=ksp(i) 
                  et1i=etp(i) 
                  ks2i=ks12 
                  et2i=et12 
                  ks1ii=ks12 
                  et1ii=et12 
                  ks2ii=ksp(i2) 
                  et2ii=etp(i2) 
               else 
                  ks1i=ks12 
                  et1i=et12 
                  ks2i=ksp(i2) 
                  et2i=etp(i2) 
                  ks1ii=ksp(i) 
                  et1ii=etp(i) 
                  ks2ii=ks12 
                  et2ii=et12 
               endif 
            endif 
         endif 
      endif 
      return 
      end 
 
 
      subroutine centra(np1,np2,fcd,b,c,epss,ksp,etp,nrzt,mrks,mret, 
     *   blx,clx,r,epsc2,a1,a2) 
      real*8 a1,a2 
      real*8 ksp(*),etp(*),r(3,2),nrzt,mrks,mret 
      real*8 b,blx,c,clx,epss,fcd,epsc2 
      if(epss.ge.0) return 
      if(epss.lt.-epsc2) go to 20 
      do 10 j1=np1,np2 
         j2=j1+1 
         call regi(fcd,b,c,ksp(j1),etp(j1),ksp(j2),etp(j2),nrzt,mrks, 
     *      mret,blx,clx,r,a1,a2) 
   10 continue 
      return 
   20 do 30 j1=np1,np2 
         j2=j1+1 
         call regii(fcd,ksp(j1),etp(j1),ksp(j2),etp(j2),nrzt,mrks,mret) 
   30 continue 
      return 
      end 
 
 
      subroutine regi(fcd,b,c,ks1,et1,ks2,et2,nrzt,mrks,mret,blx,clx,r 
     *,a1,a2) 
      real*8 a1,a2 
      real*8 r(3,2),ks1,ks2,nrzt,mrks,mret 
      real*8 b,blx,ble,br,c,clx,cle,d0,d1,d2,det,det1,det2,det3,dks, 
     * dks2,e0,e1,e2,et1,et2,fcd,g00,g01,g02,g03,g10,g11,g12 
      if(ks1.eq.0.and.et1.eq.0.and.ks2.eq.0.and.et2.eq.0) return 
      ble=2.d0*a2*b 
      cle=2.d0*a2*c+a1 
      d0=c*a1+a2*c*c 
      d1=b*cle 
      d2=a2*b*b 
      e0=cle*clx 
      e1=ble*clx+cle*blx 
      e2=ble*blx 
      br=.85*fcd 
      dks=ks2-ks1 
      det=et2-et1 
      det1=det/2. 
      det2=det*det 
      det3=det2*det 
      dks2=dks*dks 
      g00=(ks1+dks/2.)*det 
      g01=(ks1*(et1+det1)+dks*(et1/2.+det/3.))*det 
      g02=(ks1*(et1*(det+et1)+det2/3.)+dks*(et1*(et1/2+det/1.5)+ 
     *   det2/4.))*det 
      g03=(ks1*(et1*(det2+et1*(1.5*det+et1))+det3/4.)+dks*(et1* 
     *   (0.75*det2+et1*(det+et1/2.))+det3/5.))*det 
      g10=(ks1*(ks1+dks)+dks2/3)*det1 
      g11=(ks1*(ks1*(et1+det1)+dks*(et1+det/1.5))+dks2*(et1/3.+det/4.))* 
     *   det1 
      g12=(ks1*(ks1*(et1*(et1+det)+det2/3.)+dks*(et1*(et1+det/0.75)+ 
     *   det2/2.))+dks2*(et1*(et1/3.+det1)+det2/5.))*det1 
      nrzt=nrzt+br*(d0*g00+d1*g01+d2*g02) 
      mrks=mrks+br*(d0*g01+d1*g02+d2*g03) 
      mret=mret-br*(d0*g10+d1*g11+d2*g12) 
      r(1,1)=r(1,1)+br*(e0*g01+e1*g02+e2*g03) 
      r(2,1)=r(2,1)-br*(e0*g10+e1*g11+e2*g12) 
      r(3,1)=r(3,1)+br*(e0*g00+e1*g01+e2*g02) 
      return 
      end 
 
 
      subroutine regii(fcd,ks1,et1,ks2,et2,nrzt,mrks,mret) 
      real*8 ks1,ks2,nrzt,mrks,mret 
      real*8 det,dks,et1,et2,fc,fcd,g00,g01,g10 
      if(ks1.eq.0.and.et1.eq.0.and.ks2.eq.0.and.et2.eq.0) return 
      dks=ks2-ks1 
      det=et2-et1 
      g00=(ks1+dks/2.)*det 
      g01=(ks1*(et1+det/2.)+dks*(et1/2.+det/3.))*det 
      g10=(ks1*(ks1+dks)+dks*dks/3.)*det/2. 
      fc=0.85*fcd 
      nrzt=nrzt-fc*g00 
      mrks=mrks-fc*g01 
      mret=mret+fc*g10 
      return 
      end 
 
 
 
 
 