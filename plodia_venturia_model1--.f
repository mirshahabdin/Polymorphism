      program plodia_venturia_model1
C
C Model with larval host, adult host & adult parasitoid
C
      implicit double precision (a-h,o-z)
      character*2 realtostring
      parameter (neq=5,ngrid=3000,lwork=8*neq+21+ngrid,liwork=20)
      parameter (nrdens=neq,lrcont=5000*(5*nrdens+2),licont=nrdens+1)
      parameter (k=5,m=256,npts=(2*k+1)*m,nstore=5000)
      dimension y(neq),work(lwork),iwork(liwork),rpar(6),ipar(2)
      real p(m),w1(4*m),w2(m),q(409)
      common /corer/rcont(lrcont) 
      common /corei/icont(licont)
      common /h_aliases/iHL,iSHL,iSHLP,iHA
      common /p_aliases/iPA
      common /timelags/tauE,tauL,tauP,tauA
      common /mortality/DE,DL,DP,DA
      common /survival/sigmaE,sigmaL,sigmaP,sigmaA,prod
      common /comp/C
      common /paras/tauPL,tauPA,DPL,DPA,sigmaPL,sigmaPA
      common /pfparas/a,rm,rk
      common /pinocparams/rinocpval,rinocpbegin
      external fcn,solout
      if(npts.gt.nstore)stop
c==================
c variable aliases
c==================
      iHL=1
      iSHL=2
      iSHLP=3
      iHA=4
      iPA=5
c==================
c parameter values
c==================
c time lags
      tauE = 4.3D0
      tauL = 25.0d0
      tauP = 7.0D0
      tauA = 5.5D0
c
      tauPL = 20.0d0
      tauPA = 2.0d0
c mortality
      DE = 0.017D0
      DL = 0.0D0
      DP = 0.0D0
      DA = 0.1D0 !?
c
      DPL = 0.1d0
      DPA = 0.1d0
c survival
      sigmaE = exp(-DE*tauE)
      sigmaL = exp(-DL*tauL)
      sigmaP = exp(-DP*tauP)
      sigmaA = exp(-DA*tauA)
c
      sigmaPL = exp(-DPL*tauPL)
      sigmaPA = exp(-DPA*tauPA)
c fecundity
      prod = 21.0d0
c competition
      C = 2.5d-4 !1.0d-4 !4.01D-5
c parasitism
      rk = 1.0d0
      a = 0.01d0
      rm = 0.0d0
c
      rinocpval = 2.0d0
      rinocpbegin = 200.0d0
c
      niloop=1
      njloop=1
c
      open(10,file='results/fig_a_vs_rk_plodiaventuria.data')
      do 600 iloop=1,niloop
c         a=(0.1d0*iloop)/(1.0d0*niloop)
c         rpar(1)=a
         do 500 jloop=1,njloop
c            rk=(0.05d0*jloop)/(1.0d0*njloop)
c            rpar(2)=rk
            print *,'Starting iteration ',(iloop-1)*njloop+jloop
            open(2,file='plodia_venturia_graph'
     +           //realtostring((iloop-1)*njloop+jloop)//'.data')
c=========================
c dimension of the system
c=========================
            n=neq
c===========================================
c output routine is used during integration
c===========================================
            iout=1
c============================================
c initial values and endpoint of integration
c============================================
            t=0.0d0
            do 5 i=1,n
               if(i.eq.iSHL.or.i.eq.iSHLP)then
                  y(i) = 1.0d0
               else
                  y(i) = 0.0d0
               endif
 5          continue
            tend=1.0d0*nstore
c============================================
c required (relative and absolute) tolerance 
c============================================
            itol=0
            rtol=1.0D-8
            atol=rtol
c==========================================
c default values for numerical parameters
c==========================================
            idid=0
            do 10 i=1,20
               iwork(i)=0
               work(i)=0.0d0
 10         continue
            work(6)=1.0d0
            iwork(1)=1000000
            iwork(6)=ngrid
            do 12 i=1,ngrid-1
               work(20+i)=i*0.1d0
 12         continue
            work(20+ngrid)=300.0d0         
c===============================
c call of the subroutine retard 
c===============================
            ipar(1)=0
            ipar(2)=0
            call retard(n,fcn,t,y,tend,rtol,atol,itol,solout,iout,
     &           work,lwork,iwork,liwork,lrcont,licont,rpar,ipar,idid)
c            write (*,99) rtol,(iwork(j),j=17,20)
c 99         format('     tol=',D8.2,'   fcn=',I8,' step=',I5,
c     &           ' accpt=',I5,' rejct=',I5)
            if(idid.ne.1)then
c               if(idid.eq.2)then
c                  hostperiod1=512.0d0
c                  paraperiod1=512.0d0
c                  hostperiod2=512.0d0
c                  paraperiod2=512.0d0
c                  goto 300
c               else
                  hostperiod1=0.0d0
                  paraperiod1=0.0d0
                  hostperiod2=0.0d0
                  paraperiod2=0.0d0
                  goto 300
c               endif
            endif
c=======================
c calculate periodicity
c=======================
            do 200 i=1,2
               rewind 2
               open(1,file='store_model1_profile.data')
               do 101 j=1,nstore
                  read (2,'(3f50.8)') t,store1,store2
                  if(j.gt.(nstore-npts).and.i.eq.1)then
                     write (1,'(f50.8)') store1 
                  endif
                  if(j.gt.(nstore-npts).and.i.eq.2)then
                     write (1,'(f50.8)') store2 
                  endif
 101           continue
               rewind 1
               call spctrm(p,m,k,.true.,w1,w2)
               close(1)
               open(20,file='spectral_model1_profile.data')
               pmax1=0.0d0
               pmax2=0.0d0
               period1=2.0d0*m
               period2=2.0d0*m
               jnew=0
               do 150 j=2,m 
                  jnew=jnew+1
                  q(jnew)=p(j)+p(j+1)+p(j+2)+p(j+3)+p(j+4)
                  p(jnew)=p(j)
                  write(20,'(2f50.8)')(j-1)/(2.0d0*m),p(jnew)
                  if(p(jnew).gt.pmax1)then
                     pmax2=pmax1
                     period2=period1
                     period1=2.0d0*m/(j-1)
                     pmax1=p(jnew)
                  else
                     if(p(jnew).gt.pmax2)then
                        period2=2.0d0*m/(j-1)
                        pmax2=p(jnew)
                     endif
                  endif
 150           continue
               if(i.eq.1)then
                  hostperiod1=period1
                  hostperiod2=period2
               else
                  paraperiod1=period1
                  paraperiod2=period2
               endif
               close(20)
 200        continue
c================================
c print statistics for iteration
c================================
 300        write(*,'(10f25.8,2i3)')(rpar(i),i=1,6),hostperiod1,
     &           paraperiod1,hostperiod2,paraperiod2,ipar(1),ipar(2)
            close(2)
 500     continue
 600  continue
      close(10)
      stop
      end
      
      subroutine solout (nr,told,t,y,n,rpar,ipar,irtrn)
c====
c     PRINTS SOLUTION AT EQUIDISTANT OUTPUT-POINTS
c====
      implicit double precision (a-h,o-z)
      dimension y(n),rpar(6),ipar(2)
      common /h_aliases/iHL,iSHL,iSHLP,iHA
      common /p_aliases/iPA
      common /pinocparams/rinocpval,rinocpbegin
      common /intern/tout,icountHA,icountPA
      external phi
      if(y(iPA).lt.-1.0d0)then
         print *,'error - negative parasitoid #s'
         irtrn = -1
         goto 100
      endif
      if(y(iHA).lt.-1.0d0)then
         print *,'error - negative host #s'
         irtrn = -1
         goto 100
      endif
      if (nr.eq.1) then
         rpar(3)=1.0d12
         rpar(4)=1.0d12
         rpar(5)=0.0d0
         rpar(6)=0.0d0
         write (2,'(3f50.8)') t,y(iHA),y(iPA)
         tout=t+1.0d0
         icountHA=0
         icountPA=0
      else
 10      continue
         if (t.ge.tout) then
            ylagstore1=ylag(iHL,tout,phi,rpar,ipar)
            ylagstore2=ylag(iHA,tout,phi,rpar,ipar)
            ylagstore3=ylag(iPA,tout,phi,rpar,ipar)
            write (2,'(3f50.8)') tout,ylagstore2,ylagstore3 
            if(ylagstore1.lt.1.0d0)then
               icountHA=icountHA+1
               if(icountHA.gt.45)then 
                  ipar(1)=1
c                  y(iHA)=0.0d0
c                  print *,'host # < 1 for',icountHA,' days'
c                  irtrn=-1
               endif
            else
               icountHA=0
            endif
            if(ylagstore2.lt.1.0d0.and.tout.ge.rinocpbegin)then
               icountPA=icountPA+1
               if(icountPA.gt.35)then 
                  ipar(2)=1
c                  y(iPA)=0.0d0
c                  print *,'parasitoid # < 1 for',icountPA,' days'
c                  irtrn=-1
               endif
            else
               icountPA=0
            endif
            if(tout.ge.4700.0d0)then
               rpar(3)=min(rpar(3),ylagstore1)
               rpar(4)=min(rpar(4),ylagstore2)
               rpar(5)=max(rpar(5),ylagstore1)
               rpar(6)=max(rpar(6),ylagstore2)
            endif
            tout=tout+1.0d0
            goto 10
         endif
      endif
 100  continue 
      return
      end
      
      subroutine fcn(n,t,y,f,rpar,ipar)
      implicit double precision (a-h,o-z)
      dimension y(n),f(n)
      common /h_aliases/iHL,iSHL,iSHLP,iHA
      common /p_aliases/iPA
      common /timelags/tauE,tauL,tauP,tauA
      common /mortality/DE,DL,DP,DA
      common /survival/sigmaE,sigmaL,sigmaP,sigmaA,prod
      common /comp/C
      common /paras/tauPL,tauPA,DPL,DPA,sigmaPL,sigmaPA
      external phi
c=====
      RHL = prod*ylag(iHA,t-tauE,phi,rpar,ipar)*sigmaE
      rMHL = prod*ylag(iHA,t-tauE-tauL,phi,rpar,ipar)*sigmaE*
     +     sigmaL*y(iSHL)
      RHA = prod*ylag(iHA,t-tauE-tauL-tauP,phi,rpar,ipar)*
     +     sigmaE*sigmaL*sigmaP*
     +     ylag(iSHL,t-tauP,phi,rpar,ipar)+rinoch(t)
      rMHA = prod*ylag(iHA,t-tauE-tauL-tauP-tauA,phi,rpar,ipar)*
     +     sigmaE*sigmaL*sigmaP*sigmaA*
     +     ylag(iSHL,t-tauP-tauA,phi,rpar,ipar)+
     +     rinoch(t-tauA)*sigmaA
c=====
      RPA = rinocp(t)+pf(ylag(iPA,t-tauPL,phi,rpar,ipar))*
     +     ylag(iHL,t-tauPL,phi,rpar,ipar)*sigmaPL
      rMPA = rinocp(t-tauPA)*sigmaPA+
     +     pf(ylag(iPA,t-tauPL-tauPA,phi,rpar,ipar))*
     +     ylag(iHL,t-tauPL-tauPA,phi,rpar,ipar)*sigmaPL*sigmaPA
c=====
      f(iHL) = RHL-(pf(y(iPA))+C*y(iHL)+DL)*y(iHL)-rMHL
      f(iSHL) = (pf(ylag(iPA,t-tauL,phi,rpar,ipar))+
     +     C*ylag(iHL,t-tauL,phi,rpar,ipar)+
     +     -pf(y(iPA))-C*y(iHL))*y(iSHL)
      f(iSHLP) = (pf(ylag(iPA,t-tauL,phi,rpar,ipar))+
     +     -pf(y(iPA)))*y(iSHLP)
      f(iHA) = RHA-DA*y(iHA)-rMHA
c=====
      f(iPA) = RPA-DPA*y(iPA)-rMPA 
      return
      end

      double precision function phi(i,t,rpar,ipar)
      implicit double precision (a-h,o-z)
      common /h_aliases/iHL,iSHL,iSHLP,iHA
      common /p_aliases/iPA
      if(i.eq.iSHL.or.i.eq.iSHLP)then
         phi = 1.0d0
      else
         phi = 0.0d0
      endif
      return
      end

      double precision function rinoch(tdumh)
      implicit double precision (a-h,o-z)
      rinochval = 10.0d0
      rinochbegin = 0.0d0
      rinochend = 1.0d0
      if(tdumh.gt.rinochbegin.and.tdumh.lt.rinochend)then
         rinoch = rinochval
      else
         rinoch = 0.0d0
      endif
      return
      end

      double precision function rinocp(tdump)
      implicit double precision (a-h,o-z)
      common /pinocparams/rinocpval,rinocpbegin
      rinocpend = rinocpbegin+1.0d0
      if(tdump.gt.rinocpbegin.and.tdump.lt.rinocpend)then
         rinocp = rinocpval
      else
         rinocp = 0.0d0
      endif
      return
      end

      double precision function pf(padum)
      implicit double precision (a-h,o-z)
      common /pfparas/a,rm,rk
      arg = 1.0+a*(padum**(1.0-rm))/rk
      if (arg.lt.0.0d0) then 
c         print *,'pf error, -ve log arg'
      else
         pf = rk*log(arg)
      endif
      return
      end

      character*2 function realtostring(anumber)
      double precision anumber
      integer number1, i
      character*2 output
      i = 0
      output = ''
      number1 = int(anumber)
c      print *,anumber
c      print *,number1
      do 900 i = 1,2 
         output = CHAR(48 + mod(number1,10)) // output
         number1 = number1/10
 900  continue
      realtostring = output//'       '
      return 
      end
