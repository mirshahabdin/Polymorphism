      subroutine spctrm(p,m,k,ovrlap,w1,w2)
      integer k,m
      real p(m),w1(4*m),w2(m)
      logical ovrlap ! true for overlapping segments, false otherwise
c
c Reads data from input unit 1 and returns as p(j) the data's power 
c (mean square amplitude) at frequency (j-1)/(2*m) cycles per gridpoint 
c for j=1,2,...,m based on (2*k+1)*m data points (if ovrlap is set .true.) 
c or 4*k*m data points (if ovrlap is set .false.).
c The number of segments of the data is 2*k in both cases: the routine calls 
c four1 k times, each call with 2 partitions each of 2*m real data points.
c w1(1:4*m) and w2(1:m) are user-supplied workspaces.
c
      integer j,j2,joff,joffn,kk,m4,m43,m44,mm,minit
      real den,sumw,w
      minit=m
c      window(j)=1.0d0-abs(((j-1)-facm)*facp) ! Bartlett window
c      window(j)=1.0d0 ! square window
c      window(j)=1.0d0-(((j-1)-facm)*facp)**2 ! Welch window
      mm=2*m
      m4=4*m
      m44=m4+4
      m43=m4+3
      den=0.0d0
c      facm=minit
c      facp=1.0d0/minit
      sumw=0.0d0
c
c accumulate squared sum of the weights
c
      do 10 j=1,mm
         sumw=sumw+(window(j,minit))**2
 10   continue
c
c initialise spectrum
c
      do 15 j=1,m
         p(j)=0.0d0
 15   continue
c
c initialise the "save" half-buffer
c
      if(ovrlap)then
         read (1,*) (w2(j),j=1,m)
      endif
c
      do 50 kk=1,k
         do 30 joff=-1,0,1
            if(ovrlap)then
               do 20 j=1,m
                  w1(joff+j+j)=w2(j)
 20            continue
               read (1,*) (w2(j),j=1,m)
               joffn=joff+mm
               do 25 j=1,m
                  w1(joffn+j+j)=w2(j)
 25            continue
            else
               read (1,*) (w1(j),j=joff+2,m4,2)
            endif
 30      continue
         do 40 j=1,mm
            j2=j+j
            w=window(j,minit)
            w1(j2)=w1(j2)*w
            w1(j2-1)=w1(j2-1)*w
 40      continue
         call four1(w1,mm,1)
         p(1)=p(1)+w1(1)**2+w1(2)**2
         do 45 j=2,m
            j2=j+j
            p(j)=p(j)+0.5d0*(w1(j2)**2+w1(j2-1)**2+
     +           w1(m44-j2)**2+w1(m43-j2)**2)
 45      continue
         den=den+sumw
 50   continue
      den=m4*den
c
c normalise the output
c
      do 60 j=1,m
         p(j)=p(j)/den
 60   continue
      return
      end

      real function window(j,minit)
      integer j,minit
      real facm,facp
      facm=minit
      facp=1.0d0/minit
c      window=1.0d0-abs(((j-1)-facm)*facp) ! Bartlett window
c      window=0.5d0*(1.0d0-cos(3.14159265359d0*(j-1)*facp)) ! Hann window
      window=1.0d0 ! square window
c      window=1.0d0-(((j-1)-facm)*facp)**2 ! Welch window
      return
      end

      subroutine four1(datam,nn,isign)
      integer isign,nn
      real datam(2*nn)
c
c Replaces datam(1:2*nn) by its discrete Fourier transform, if isign is 
c input as 1; or replaces datam(1:2*nn) by nn times its inverse discrete
c Fourier transform, if isign is input as -1. 
c datam is a real array of length 2*nn (nn complex numbers). 
c nn must be an integer power of 2 - need to check for this!
c
      integer i,istep,j,m,mmax,n
      real tempi,tempr
      double precision theta,wi,wpi,wpr,wr,wtemp
      n=2*nn
c
c Bit-reversal, exchanging the pairs of complex numbers
c
      j=1
      do 10 i=1,n,2
         if(j.gt.i)then
            tempr=datam(j)
            tempi=datam(j+1)
            datam(j)=datam(i)
            datam(j+1)=datam(i+1)
            datam(i)=tempr
            datam(i+1)=tempi
         endif
         m=n/2 
 1       if((m.ge.2).and.(j.gt.m))then
            j=j-m
            m=m/2
            goto 1
         endif
         j=j+m
 10   continue
c
c Danielson-Lanczos section
c Outer loop executed log_2 nn times
c
      mmax=2
 20   if(n.gt.mmax)then
         istep=2*mmax
         theta=6.28318530717959d0/(isign*mmax)
         wpr=-2.0d0*(sin(0.5d0*theta)**2)
         wpi=sin(theta)
         wr=1.0d0
         wi=0.0d0
         do 30 m=1,mmax,2
            do 25 i=m,n,istep
               j=i+mmax ! Danielson-Lanczos formula:
               tempr=sngl(wr)*datam(j)-sngl(wi)*datam(j+1)
               tempi=sngl(wr)*datam(j+1)+sngl(wi)*datam(j)
               datam(j)=datam(i)-tempr
               datam(j+1)=datam(i+1)-tempi
               datam(i)=datam(i)+tempr
               datam(i+1)=datam(i+1)+tempi   
 25         continue
            wtemp=wr ! trigonometric recurrence
            wr=wr*wpr-wi*wpi+wr
            wi=wi*wpr+wtemp*wpi+wi
 30      continue
         mmax=istep
         goto 20
      endif
      return
      end
