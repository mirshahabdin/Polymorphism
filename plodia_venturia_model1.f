      program plodia_venturia_model1
C
C Model with larval host, adult host & adult parasitoid
C
      implicit double precision (a-h,o-z)
      character*2 realtostring
      parameter (nh=3,nhstrain=1,npstrain=1,neq=nh*nhstrain+npstrain)
      parameter (ngrid=3000,lwork=8*neq+21+ngrid,liwork=20)
      parameter (nrdens=neq,lrcont=5000*(5*nrdens+2),licont=nrdens+1)
      parameter (k=5,m=256,npts=(2*k+1)*m,nstore=5000)
      dimension y(neq),work(lwork),iwork(liwork),rpar(6),ipar(2)
      real p(m),w1(4*m),w2(m),q(409)
      common /corer/rcont(lrcont)
      common /corei/icont(licont)
      common /h_aliases/iHL,iSHL,iHA
      common /p_aliases/iPA
      common /timelags/tauE,tauL,tauP,tauA
      common /mortality/DE,DL,DP,DA
      common /survival/sigmaE,sigmaL,sigmaP,sigmaA,bHA(nhstrain),bHA0
      common /enc/eta(nhstrain,npstrain),Aenc
      common /mutant/eps
      common /cost/cH,cP
      common /comp/C
      common /paras/tauPL,tauPA,DPL(npstrain),DPL0,DPA,sigmaPL(npstrain),sigmaPA
      common /pfparas/a,rm,rk
      common /pinocparams/rinocpval,rinocpbegin
      common /strains/resi(nhstrain),viru(npstrain)
      external fcn,solout
      if(npts.gt.nstore)stop
c==================
c variable aliases
c==================
      iHL=1
      iSHL=2
      iHA=3
      iPA=4
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
      DA = 0.1D0
c
      DPL = 0.1d0
      DPL0 = 0.01d0
      DPA = 0.1d0
c mutants
      if (nhstrain.gt.1.and.nhstrain.gt.1) then
         eps = 1.d-5
      else
         eps = 0.d0
      endif
c cost
c      cH
c      cP
c survival
      sigmaE = exp(-DE*tauE)
      sigmaL = exp(-DL*tauL)
      sigmaP = exp(-DP*tauP)
      sigmaA = exp(-DA*tauA)
c
      sigmaPL = exp(-DPL*tauPL)
      sigmaPA = exp(-DPA*tauPA)
c fecundity
      bHA = 21.0d0
      bHA0 = 80.d0
c competition
      C = 2.0d-4 !1.0d-4 !4.01D-5
c parasitism
      rk = 0.01d0
      a = 0.01d0
      rm = 0.0d0
c encapsulation
      Aenc = 0.5d0
      eta = 0.d0
      do 800 ih=1,nhstrain
      do 801 ip=1,npstrain
         eta(ih,ip) = 1.d0/(1.d0+exp(-Aenc*(resi(ih)-viru(ip))))
801   continue
800   continue
      
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
            open(2,file='results/plodia_venturia_graph'
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
               if(i.eq.iSHL)then
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
               open(1,file='results/store_model1_profile.data')
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
               open(20,file='results/spectral_model1_profile.data')
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
      common /h_aliases/iHL,iSHL,iHA
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
         write (2,'(3f20.8)') t,y(iHA),y(iPA)
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
      parameter (nh=3,nhstrain=1,npstrain=1)
      common /h_aliases/iHL,iSHL,iHA
      common /p_aliases/iPA
      common /timelags/tauE,tauL,tauP,tauA
      common /mortality/DE,DL,DP,DA
      common /mutant/eps
      common /survival/sigmaE,sigmaL,sigmaP,sigmaA,bHA(nhstrain),bHA0
      common /comp/C
      common /paras/tauPL,tauPA,DPL(npstrain),DPL0,DPA,sigmaPL(npstrain),sigmaPA
      common /enc/eta(nhstrain,npstrain),Aenc
      external phi

c=====
      do 501 ih=1,nhstrain
         RHL = (1.d0-eps)*bHA(ih)*sigmaE*
     +     ylag((ih-1)*nh+iHA,t-tauE,phi,rpar,ipar)
         do 502 il=1,ih-1
            RHL = RHL+eps/(nhstrain-1)*bHA(il)*sigmaE*
     +     ylag((il-1)*nh+iHA,t-tauE,phi,rpar,ipar)
502   continue
         do 503 il=ih+1,nhstrain
            RHL = RHL+eps/(nhstrain-1)*bHA(il)*sigmaE*
     +     ylag((il-1)*nh+iHA,t-tauE,phi,rpar,ipar)
503   continue
         rMHL = bHA(ih)*sigmaE*sigmaL*y(iSHL)*
     +     ylag((ih-1)*nh+iHA,t-tauE-tauL,phi,rpar,ipar)
         DHL = (pf(y((ih-1)*nh+iPA))+C*y((ih-1)*nh+iHL)+DL)*
     +     y(ih*nh+iHL)
         RHA = bHA(ih)*sigmaE*sigmaL*sigmaP*
     +     ylag((ih-1)*nh+iHA,t-tauE-tauL-tauP,phi,rpar,ipar)*
     +     ylag((ih-1)*nh+iSHL,t-tauP,phi,rpar,ipar)+rinoch(t)
         rMHA = bHA(ih)*sigmaE*sigmaL*sigmaP*sigmaA*
     +     ylag((ih-1)*nh+iHA,t-tauE-tauL-tauP-tauA,phi,rpar,ipar)*
     +     ylag((ih-1)*nh+iSHL,t-tauP-tauA,phi,rpar,ipar)+
     +     rinoch(t-tauA)*sigmaA
         DHA = DA*y((ih-1)*nh+iHA)
         f((ih-1)*nh+iHL) = RHL-rMHL-DHL
         f((ih-1)*nh+iHA) = RHA-rMHA-DHA
         f((ih-1)*nh+iSHL) = 
     +     (pf(ylag((ih-1)*nh+iPA,t-tauL,phi,rpar,ipar))+
     +     C*ylag((ih-1)*nh+iHL,t-tauL,phi,rpar,ipar)-
     +     pf(y((ih-1)*nh+iPA))-C*y((ih-1)*nh+iHL))*y((ih-1)*nh+iSHL)

501   continue
c=====
      do 504 ip=1,npstrain
         iauxn = nh*nhstrain
         RPA = rinocp(t)
         RPA = 0.d0
         do 505 ih=1,nhstrain
            RPA = RPA+(1.d0-eps)*(1.d0-eta(ih,ip))*
     +     pf(ylag(iauxn+ip,t-tauPL,phi,rpar,ipar))*
     +     ylag((ih-1)*nh+iHL,t-tauPL,phi,rpar,ipar)*
     +     sigmaPL(ip)
         do 506 il=1,ip-1
            RPA = RPA+eps/(npstrain-1)*(1.d0-eta(ih,il))*
     +     pf(ylag(iauxn+il,t-tauPL,phi,rpar,ipar))*
     +     ylag((ih-1)*nh+iHL,t-tauPL,phi,rpar,ipar)*
     +     sigmaPL(il)
506   continue
         do 507 il=ip+1,npstrain
            RPA = RPA+eps/(npstrain-1)*(1.d0-eta(ih,il))*
     +     pf(ylag(iauxn+il,t-tauPL,phi,rpar,ipar))*
     +     ylag((ih-1)*nh+iHL,t-tauPL,phi,rpar,ipar)*
     +     sigmaPL(il)
507   continue
            rMPA = rinocp(t-tauPA)*sigmaPA
            rMPA = 0.d0
            rMPA = rMPA+(1.d0-eta(ih,ip))*
     +     pf(ylag(iauxn+ip,t-tauPL-tauPA,phi,rpar,ipar))*
     +     ylag((ih-1)*nh+iHL,t-tauPL-tauPA,phi,rpar,ipar)*
     +     sigmaPL(ip)*sigmaPA
505   continue
         DPA = DPA*y(iauxn+ip)
         f(iauxn+ip) = RPA-rMPA-DPA
         write(*,*) t
504   continue
c=====
      return
      end

      double precision function phi(i,t,rpar,ipar)
      implicit double precision (a-h,o-z)
      common /h_aliases/iHL,iSHL,iHA
      common /p_aliases/iPA
      if(i.eq.iSHL)then
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
