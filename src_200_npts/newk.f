********************************************************************
* Kerl.f             Compute the kernel library                    *
*                                                                  *
* Modified a little by YG + DG, from VAX version to SUN version    *
*                                                                  *
*           Input files:  sprm1bin                                 *
*                         libreria1.inn                            *
*           Output file:  libre010.b                               *
*                                                    July 1999     *
*                                                                  *
********************************************************************

c  computes kernels for moment tensor inversion using 
c  band passed very long period data.
c  version for azimutally invariant kernels

      program kerl

cfb
c     VAR DECLARATION
      
      implicit none
      integer jst(24),jco(24),jkl(24),npus,nmax,ir1,irr,nst,j,i
      real synt(200),sla,slo,elev,dts,tin,pmin,fc,dis,depp
      character*10 name
c     character chdep*3  bf chdep used only in interattive mode


c      dimension synt(200),jst(24),jco(24),jkl(24)
c      character*10 name
c      character chdep*3
cfb
      common/stat/sla(3),slo(3),elev
      common/lplp/dts,npus,tin,pmin,fc,nmax

      data jst/1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3/
      data jco/1,1,1,1,2,3,3,3,3,1,1,1,2,2,2,3,3,1,2,2,2,3,3,3/
      data jkl/1,2,3,5,4,1,2,3,5,2,3,4,2,3,4,5,6,6,2,3,6,2,3,6/
c---  
      fc=1.0e24
      npus=200

      pmin=30.

      tin=0.
      dts=10.
      elev=0.
c---  

c     write(6,'(/'' Enter the depth (in two digits) : '')')
c     write(6,'(/'' -- Such as:  08  10  21 --  '',$)')
c     read(5,'(a3)' ) chdep
c     name='libre0'//chdep(1:2)//'.b'
c     name='libre'//chdep//'.b'
      name='libre200.b'

      depp=200.
c      depth=33
c      name='libre010.b'

c---
      nmax=1000              !!! all modes
c     nmax=1                 !!! only fundamental modes
c---
      open(8,file=name,form='unformatted',access='direct'
     .      ,status='new',recl=4*(6+200))
c      open(2,file='libreria.log',status='new')
c---  
      open(7,file='sprm1bin',form='unformatted',access='direct'
     .      ,status='old',recl=4*(1024))
      call setprem(7,0)
      call findd(depp,ir1)
c---
      irr=0
      open(1,file='libreria1.inn',status='old')
1     read(1,'(i4,7f9.4)',end=99) nst,dis,(sla(j),slo(j),j=1,3)
c---
      do i=1,24
c      write(2,'(4i5))') nst,jst(i),jco(i),jkl(i)
      call mkkern(synt,ir1,jst(i),jco(i),jkl(i))
      irr=irr+1  
      write(8,rec=irr) nst,dis,i,jst(i),jco(i),jkl(i),synt
      enddo
c---
      go to 1
99    end

C--------------------------------------------------------------
C--------------------------------------------------------------

      subroutine mkkern(synt,ir1,icode,icomp,kl)

cfb
c     VAR DECLARATION

      implicit none
      integer npus,nmax,nord,lord,ir1,icode,icomp,kl,lmax,i,ll,ifgot,nn
      real sla,slo,elev,dts,tin,pmin,fc,amp,aker,ar1,ar2,r1,r2,synt(200)
      real slat,slon,sour(9),per,tpi,pi
      character*1 ityp

c      dimension synt(200),sour(9)
cfb
c      character*1 ityp
      common/stat/sla(3),slo(3),elev
      common/lplp/dts,npus,tin,pmin,fc,nmax
      common/premdata/
     .       amp(12),aker(10),ar1(4),ar2(4),r1,r2,nord,ityp,lord

      pi=acos(-1.)
      tpi=2.*pi
      lmax=1000
      slat=sla(icode)
      slon=slo(icode)
      do i=1,9
        sour(i)=0.0
      enddo
      sour(kl+3)=fc
c---
      do i=1,npus
        synt(i)=0.0
      enddo
c---  
      nord=0              !!  S mode selection and summation
      do 2 nn=1,nmax
c      write(5,*)'            S',nord
        lord=0
        do 1 ll=1,lmax
          if(nord.le.1.and.lord.eq.1) go to 1
          ityp='S'
          call getprem(ir1,ifgot)
          per=tpi/amp(1)
          if(ifgot.ne.1.or.per.lt.pmin) then
            if(lord.eq.0) go to 1
            if(lord.eq.1) go to 3
            go to 2
          endif
        call admode(synt,slat,slon,sour,icomp)
1     lord=lord+1
2     nord=nord+1
3     continue
c---  
      if(icomp.eq.1) go to 6   !!  T mode selection and summation
      nord=0
      do 5 nn=1,nmax
c      write(5,*)'            T',nord
        lord=1
        do 4 ll=1,lmax
          if(nord.eq.0.and.lord.eq.1) go to 4
          ityp='T'
          call getprem(ir1,ifgot)
          per=tpi/amp(1)
          if(ifgot.ne.1.or.per.lt.pmin) then
            if(lord.eq.1) go to 6
            go to 5
          endif
        call admode(synt,slat,slon,sour,icomp)
4     lord=lord+1
5     nord=nord+1
6     continue
c---  
      return
      end

C--------------------------------------------------------------------
C--------------------------------------------------------------------

      subroutine admode(synt,slat,slon,sour,icomp)

cfb   VAR DECLARATION
 
      implicit none
      integer npus,nmax,nord,lord,icomp,isig,n,i,im
      real sour(9),synt(200),deom(1000),dts,tin,pmin,fc,amp,aker
      real slat,slon,r1,r2,ar1,ar2
      complex vs(1000),vr(1000),apex,ci,ce,com,cin,omeg
      character*1 ityp

c      dimension sour(9),synt(200),deom(1000)
c      complex vs(1000),vr(1000),apex,ci,ce,com,cin,omeg
cfb
      common/lplp/dts,npus,tin,pmin,fc,nmax
      common/premdata/
     .       amp(12),aker(10),ar1(4),ar2(4),r1,r2,nord,ityp,lord


      call cpsplit(apex,vr,vs,slat,slon,sour,icomp)

      isig=1
      if(icomp.eq.2) isig=-1
      ci=(0.,1.)

      if(lord.le.10) then     !!  apply complete split formulation
        call cpdom(deom)
        n=2*lord+1
        do im=1,n
          omeg=cmplx(amp(1)+deom(im))
          ce=ci*omeg-.5*omeg*amp(2)
          com=cexp(ce*tin)*vs(im)*vr(im)
          cin=cexp(ce*dts)
          synt(1)=isig*real(com)+synt(1)
          do i=2,npus
            com=com*cin
            synt(i)=isig*real(com)+synt(i)
          enddo
        enddo
      return
      endif

      ce=amp(1)*(ci-.5*amp(2))    !!  apply asymptotic formulation
      com=cexp(ce*tin)*apex
      cin=cexp(ce*dts)
      synt(1)=isig*real(com)+synt(1)
      do i=2,npus
        com=com*cin
        synt(i)=isig*real(com)+synt(i)
      enddo

      return
      end

c--------------------------------------------------------------------

      subroutine cpsplit(apex,vr,vs,slat,slon,sour,icomp)

cfb   VAR DECLARATION
   
      implicit none
      integer no,it,lo,j,icomp
      real rec(5),x(2000),xp(2000),xcosec(2000),sour(9)
      real slat,slon,amp,aker,ar1,ar2,r1,r2
      complex ss,vr(1),vs(1),apex

c      dimension rec(5),x(2000),xp(2000),xcosec(2000),sour(9)
c      complex ss,vr(1),vs(1),apex
cfb

      common/premdata/amp(12),aker(10),ar1(4),ar2(4),r1,r2,no,it,lo

      apex=(0.,0.)
      rec(1)=slat
      rec(2)=slon
      rec(3)=0.
      rec(4)=0.
      rec(5)=0.
      rec(icomp+2)=1.
      call reprem(rec,vr,x,xp,xcosec)
      call soprem(sour,vs,x,xp,xcosec)
      if(lo.ne.0) then
        do j=1,lo+1
          vr(2*lo+2-j)=vr(lo+2-j)
          vs(2*lo+2-j)=vs(lo+2-j)
        enddo
        do j=1,lo
          if(mod(j,2).eq.0) ss=(1.,0.)
          if(mod(j,2).eq.1) ss=(-1.,0.)
          vr(lo+1-j)=ss*conjg(vr(lo+1+j))
          vs(lo+1-j)=ss*conjg(vs(lo+1+j))
        enddo
      endif
      do j=1,2*lo+1
        apex=apex+vr(j)*vs(j)
      enddo

      return
      end

c--------------------------------------------------------------------

      subroutine cpdom(deom)

cfb   VAR DECLARATION
     
      implicit none
      integer nord,lord,im,m
      real deom(1),amp,aker,ar1,ar2,r1,r2,alf,bet,rll
      character*1 ityp

c      dimension deom(1)
cfb

      common/premdata/
     .   amp(12),aker(10),ar1(4),ar2(4),r1,r2,nord,ityp,lord

      if(lord.eq.0) then
        deom(1)=0.0
        return
      endif
      alf=amp(1)*amp(6)
      bet=amp(1)*amp(12)
      rll=lord*(lord+1)
      do im=1,2*lord+1
        m=im-lord-1
        deom(im)=m*bet+alf*(1.-float(3*m**2)/rll)
      enddo
      return
      end

C--------------------------------------------------------------------
C--------------------------------------------------------------------

      SUBROUTINE SETPREM(LUMP,LUPER)

cfb   VAR DECLARATION
 
      implicit none
      integer mtot,lnts,lntt,indsfr,kntsfr,indtor,knttor
      integer nbatch,lump1,luper1,lump,luper
      real rmp
cfb

      COMMON/SETKER/MTOT,LNTS,LNTT,RMP(48),INDSFR(330),KNTSFR(330)
     1   ,INDTOR(300),KNTTOR(300),NBATCH,LUMP1,LUPER1

      LUPER1=LUPER
      LUMP1=LUMP

      CALL BFFI(LUMP1,MTOT,3+48,1)

      CALL BFFI(LUMP1,INDSFR,660,2)

      CALL BFFI(LUMP1,INDTOR,600,3)

      NBATCH=(MTOT+255)/256
      RETURN
      END


C--------------------------------------------------------------


      SUBROUTINE FINDD(DEPTH,IND)
      
cfb   VAR DECLARATION
 
      implicit none
      integer ind,ini
      real r,dum1,dum2,rs,h,hb,depth
cfb

      COMMON/SETKER/DUM1(3),R(48),DUM2(1263)
      INI=1
      RS=1.E0-(DEPTH/6371.0)
   10 INI=INI+1
      IF(INI.GT.48) PAUSE 'ERROR IN FINDD'
      H=RS-R(INI-1)
      HB=RS-R(INI)
      IF(H*HB.GT.0.E0) GO TO 10
      IND=INI-1
      RETURN
      END


C--------------------------------------------------------------


      SUBROUTINE GETPREM(IR1,IFGOT)
      
cfb   VAR DECLARATION
      
      implicit none
      integer mtot,lnts,lntt,indsfr,kntsfr,indtor,knttor
      integer nbatch,lump1,luper1,nord,lord,ir1,ifgot,ibtchl
      integer ibatch,ivar,irecmp,irc1,irc2,ind,i,mdnum
      real rmp,b,bker,eig1,eig2,amp
      real aker,ar1,ar2,r1,r2
      character iq*1   ! cfb
cfb
      
      COMMON/SETKER/MTOT,LNTS,LNTT,RMP(48),INDSFR(330),KNTSFR(330)
     1   ,INDTOR(300),KNTTOR(300),NBATCH,LUMP1,LUPER1
      COMMON/GETMOD/B(256,12),BKER(256,10),EIG1(256,4),EIG2(256,4)
      COMMON/PREMDATA/
     .   AMP(12),AKER(10),AR1(4),AR2(4),R1,R2,NORD,IQ,LORD
      DATA IBTCHL/0/
      SAVE IBTCHL

      IFGOT=0
      IF(IQ.NE.'S') GOTO 10
      IF(NORD.GE.KNTSFR(LORD+1)) RETURN
      MDNUM=INDSFR(LORD+1)+NORD
cfb      IQ=3
      IQ='3'
      GOTO 20

   10 IF(IQ.NE.'T') PAUSE ' ERROR 1 IN GETMOD '
      IF(NORD.GE.KNTTOR(LORD)) RETURN
      MDNUM=INDTOR(LORD)+NORD
cfb      IQ=2
      IQ='2'

   20 IBATCH=(MDNUM+255)/256
      IF(IBATCH.EQ.IBTCHL) GOTO 50

      IBTCHL=IBATCH
      DO 30 IVAR=1,3
      IRECMP=3+NBATCH*(IVAR-1)+IBATCH
   30 CALL BFFI(LUMP1,B(1,4*(IVAR-1)+1),4*256,IRECMP)

      IRC1=IRECMP+NBATCH*IR1
      IRC2=NBATCH+IRC1
      CALL BFFI(LUMP1,EIG1,4*256,IRC1)
      CALL BFFI(LUMP1,EIG2,4*256,IRC2)

   50 IND=MDNUM-256*(IBATCH-1)
      DO 60 I=1,12
   60 AMP(I)=B(IND,I)
      DO 80 I=1,4
      AR1(I)=EIG1(IND,I)
   80 AR2(I)=EIG2(IND,I)
   90 IFGOT=1
      R1=RMP(IR1)
      R2=RMP(IR1+1)
      RETURN
      END


C--------------------------------------------------------------


C     SOPREM RETURNS THE EXITATION OF NORMAL MODES IN EVEC
C     INPUT: SOUR(I),I=1,9
C            LATITUDE,LONGITUDE,DEPTH,XM(1-6)
C     OUTPUT:COMPLEX EVEC(LORD+1) WHERE LORD IS ANGULAR ORDER
C     WORKSPACES X,XP,XCOSEC HAVE TO BE DIMENSIONED (LORD+1)
C     EVEC(1) - M=0   EVEC(2) - M=1   EVEC(LORD+1) - M=LORD
C     THIS IS A SINGLE PRECISION VERSION OF EXTAN
C     EVEC(-M)=(-1)**M*CONJG(EVEC(M)

      SUBROUTINE SOPREM(SOUR,EVEC,X,XP,XCOSEC)

cfb   VAR DECLARATION

      implicit none
      integer no,l,lp1,i,j
      real sour(9),x(1),xp(1),xcosec(1),om,qu,vaa,gvel,ellip,dum1
      real dum2,u,up,v,vp,u2,u2p,v2,v2p,r1,r2,haa,pi,pib2,rs,rotsp
      real theta,phi,h,hn,rhn,hn2,hn3,a,b,vv,vvp,uu,uup,b0,b1,b2,es
      real st,cosec,cot,fl3,xm,y,yp,ct
      complex evec(3),fact,efac,e(6)
      character iq*1  ! cfb

c      COMPLEX EVEC(3),FACT,EFAC,E(6)
c      DIMENSION SOUR(9),X(1),XP(1),XCOSEC(1)
cfb

      COMMON/PREMDATA/OM,QU,VAA,HAA,GVEL,ELLIP,DUM1(5),ROTSP
     1 ,DUM2(10),U,UP,V,VP,U2,U2P,V2,V2P,R1,R2,NO,IQ,L
      DATA PI,PIB2/3.141592      ,1.570796      /

      RS=1.0-SOUR(3)/6371.0
      THETA=PIB2-SOUR(1)*PI/180.E0
      PHI=SOUR(2)*PI/180.E0
      H=RS-R1
      HN=R2-R1
      RHN=1.E0/HN
      HN2=RHN*RHN
      HN3=RHN*HN2
      A=HN3*(HN*(VP+V2P)+
     1         2.E0*(V-V2))
      B=HN2*(3.E0*(V2-V)-
     1         HN*(V2P+2.D0*VP))
      VV=V+H*(VP+H*(B+H*A))
      VVP=VP+H*(2.E0*B+3.E0*H*A)
cfb      IF(IQ.EQ.2) GOTO 20
      IF(IQ.EQ.'2') GOTO 20
      A=HN3*(HN*(UP+U2P)+
     1         2.E0*(U-U2))
      B=HN2*(3.E0*(U2-U)-
     1         HN*(U2P+2.E0*UP))
      UU=U+H*(UP+H*(B+H*A))
      UUP=UP+H*(2.E0*B+3.E0*H*A)

   20 B0=SQRT(FLOAT(2*L+1)/(4.E0*PI))
      B1=-.5D0*B0*SQRT(FLOAT(L*(L+1)))
      B2=0.E0
      IF(L.GT.0) B2=-.25*B1*SQRT(FLOAT((L-1)*(L+2)))
      LP1=L+1
      DO 1 I=1,LP1
    1 EVEC(I)=(0.,0.)
      ES=VVP+(UU-VV)/RS
      IF(THETA.NE.0.E0.AND.THETA.NE.PI) GOTO 100
      FACT=(1.,0.)
      IF(THETA.EQ.PI) FACT=(-1.,0.)
cfb       IF(IQ.EQ.2) GOTO 110
      IF(IQ.EQ.'2') GOTO 110
      EVEC(1)=SOUR(4)*B0*UUP+(SOUR(5)+SOUR(6))
     1     *(UU-.5E0*FLOAT(L*(L+1))*VV)*B0/RS
      IF(L.EQ.0) RETURN
      EVEC(2)=FACT*ES*B1*(SOUR(7)-FACT*CMPLX(0.E0,SOUR(8)))
      IF(L.EQ.1) RETURN
      EVEC(3)=(SOUR(5)-SOUR(6))*2.E0*B2*VV/RS+SOUR(9)*
     1     FACT*CMPLX(0.E0,-4.E0*B2*VV/RS)
      RETURN

  110 IF(L.EQ.0) RETURN
      EVEC(2)=-FACT*B1*ES*(FACT*CMPLX(0.E0,SOUR(7))-SOUR(8))
      IF(L.EQ.1) RETURN
      EVEC(3)=-FACT*B2*VV*2.E0/RS*CMPLX(0.E0,SOUR(5)-SOUR(6))
     1     -4.E0*B2*VV*SOUR(9)/RS
      RETURN
  100 CALL LEGRG(THETA,L,L,X,XP,XCOSEC)
      CT=COS(THETA)
      ST=SIN(THETA)
      COSEC=1.E0/ST
      COT=CT/ST
      FL3=L*(L+1)
      FACT=(1.,0.)
      EFAC=CEXP(CMPLX(0.E0,-PHI))
      DO 11 I=1,LP1
      XM=I-1
      Y=X(I)
      YP=XP(I)
cfb      IF(IQ.EQ.2) GOTO 21
      IF(IQ.EQ.'2') GOTO 21
      E(1)=UUP*Y
      E(2)=(UU*Y-VV*(COT*YP-((XM*COSEC)**2-FL3)*Y))/RS
      E(3)=(UU*Y+VV*(COT*YP-((XM*COSEC)**2)*Y))/RS
      E(4)=ES*YP
      E(5)=CMPLX(0.E0,-XM*ES*COSEC*Y)
      E(6)=CMPLX(0.E0,-2.E0*XM*VV*COSEC*(YP-COT*Y)/RS)
      GOTO 22
   21 E(1)=(0.,0.)
      E(2)=CMPLX(0.E0,-VV*XM*COSEC*(YP-COT*Y)/RS)
      E(3)=-E(2)
      E(4)=CMPLX(0.E0,-XM*COSEC*Y*ES)
      E(5)=-YP*ES
      E(6)=VV*(2.E0*COT*YP-(2.E0*((XM*COSEC)**2)-FL3)*Y)/RS
   22 DO 23 J=1,6
   23 EVEC(I)=EVEC(I)+SOUR(J+3)*E(J)
      EVEC(I)=EVEC(I)*FACT
      FACT=FACT*EFAC
   11 CONTINUE
      RETURN
      END


C--------------------------------------------------------------


C     REPREM CALCULATES THE RECEIVER VECTOR REVEC
C     INPUT:  RECE(5) , LAT,LONG,VERT,N-S,E-W  (????)
C     OUTPUT: REVEC(LORD+1)  REVEC(1) - M=0 ,ETC.
C     WORKSPACES HAVE TO BE ALLOCATED TO (LORD+1)
C     REVEC(-M)=(-1)**M*(REVEC(M)

      SUBROUTINE REPREM(REC,REVEC,X,XP,XCOSEC)
      
cfb   VAR DECLARATION
 
      implicit none
      integer dum1,dum2,no,l,lp1,i
      real rec(5),x(1),xcosec(1),xp(1),om,qu,vaa,haa,gvel,ellip
      real u2,u2p,v2,v2p,r1,r2,rotsp,pi,pib2,theta,phi,cosec,xm
      real u,up,v,vp   !! bf: var rotps never used
      complex revec(1),fact,efac
      character iq*1   ! cfb

c      COMPLEX REVEC(1),FACT,EFAC
c      DIMENSION REC(5),X(1),XCOSEC(1),XP(1)
cfb
      COMMON/PREMDATA/OM,QU,VAA,HAA,GVEL,ELLIP,DUM1(5),ROTSP
     1 ,DUM2(10),U,UP,V,VP,U2,U2P,V2,V2P,R1,R2,NO,IQ,L
      DATA PI,PIB2/3.1415926     ,1.5707963     /

      THETA=PIB2-REC(1)*PI/180.E0
      PHI=REC(2)*PI/180.E0
      COSEC=1.E0/SIN(THETA)
      FACT=(1.,0.)
      EFAC=CEXP(CMPLX(0.E0,PHI))
      CALL LEGRG(THETA,L,L,X,XP,XCOSEC)
      LP1=L+1
      DO 10 I=1,LP1
      XM=I-1
cfb      IF(IQ.EQ.2) GOTO 20

      IF(IQ.EQ.'2') GOTO 20
      REVEC(I)=(REC(3)*VAA*X(I)+REC(4)*HAA*XP(I)
     1       +REC(5)*CMPLX(0.E0,XM*COSEC*HAA)*X(I))*FACT
      GOTO 9
   20 REVEC(I)=(CMPLX(0.E0,XM*COSEC*X(I)*REC(4))
     1        -REC(5)*XP(I))*HAA*FACT
    9 FACT=FACT*EFAC
   10 CONTINUE
      RETURN
      END


C--------------------------------------------------------------


C     THIS IS A SINGLE PRECISION VERSION OF LEGROT
C     WHICH CALCULATES AND ROTATES SPHERICAL HARMONICS.

      SUBROUTINE LEGRG(THETA,L,MMAX,X,XP,XC)
      
cfb   VAR DECLARATION

      implicit none
      integer l,mmax,num,i,ind
      real d(3,1001),x(1),xp(1),xc(1),theta,pi,fac,fac1,rm

c     DIMENSION D(3,1001),X(1),XP(1),XC(1)
cfb

      DATA PI/3.1415926         /

      CALL ROTMG(1,L,THETA,D,3,1001)
      FAC=SQRT(FLOAT(2*L+1)/(4.E0*PI))
      FAC1=SQRT(FLOAT(L*(L+1)))
      NUM=MMAX+1
      DO 10 I=1,NUM
      RM=FLOAT(I-1)
      IND=L+I
      X(I)=FAC*D(2,IND)
      XP(I)=-FAC*(FAC1*D(3,IND)+RM*D(2,IND)/SIN(THETA))
   10 CONTINUE
      RETURN
      END


C--------------------------------------------------------------


C     THIS IS A SINGLE PRECISION VERSION OF ROTMX2
C     WHICH ROTATES A SECOND RANK TENSOR (I THINK).

      SUBROUTINE ROTMG(NMAX,L,THETA,D,ID1,ID2)

cfb   VAR DECLRATION
 
      implicit none
      integer n,nmax,id1,id2,isup,nm,nmp1,lp1,lm1,lp2,nrow,in1ct,in2ct
      integer nmaxp1,lmn,im1,m1,nm2,nit,i,mult,nm2p1,in1
      integer in2,im2ct,im2,l,m2,im1ct
      real theta,big,small,dlbig,dlsml,pi,th,shth,chth,sth,cth,dlogf
      real dlogs,rm1,temp,rmod,sgn,csum,fac,sfac,d(3,1001),t1

c     DIMENSION D(ID1,ID2)
cfb
      DATA BIG,SMALL,DLBIG,DLSML/1.E35,1.E-35,35.E0,-35.E0/
      DATA PI/3.1415926         /

cuk!!!!!!!!!      FLOAT(N)=N
      TH=THETA
      IF((TH.GT.PI).OR.(TH.LT.0.E0)) STOP 'ILLEGAL ARG IN ROTMX2'
      IF(L.NE.0) GOTO 350
      D(1,1)=1.E0
      RETURN
350   ISUP=1
      IF(TH.LE.PI/2.E0) GOTO 310
      TH=PI-TH
      ISUP=-1
310   NM=2*L+1
      NMP1=NM+1
      LP1=L+1
      LM1=L-1
      LP2=L+2
      NROW=2*NMAX+1
      NMAXP1=NMAX+1
      LMN=L-NMAX
      IF(TH.NE.0.E0) GOTO 320
      DO 330 IM1CT=1,NROW
      IM1=IM1CT+LMN
      DO 330 IM2=LP1,NM
      D(IM1CT,IM2)=0.E0
      IF(IM1.EQ.IM2) D(IM1CT,IM2)=1.E0
330   CONTINUE
      GOTO 400
320   CONTINUE

C                           ZERO L.H.S. OF MATRIX

      DO 340 IM1=1,NROW
      DO 340 IM2=1,LP1
340   D(IM1,IM2)=0.E0

C                           SET UP PARAMETERS

      SHTH=SIN(0.5E0*TH)
      CHTH=COS(0.5E0*TH)
      STH=2.E0*SHTH*CHTH
      CTH=2.E0*CHTH*CHTH-1.E0
      DLOGF=ALOG10(CHTH/SHTH)
      DLOGS=ALOG10(SHTH)

C                      ITERATE FROM LAST COLUMN USING 1. AS STARTING VALUE

      DO 10 IM1CT=1,NROW
      IM1=IM1CT+LMN
      M1=IM1-LP1
      RM1=M1
      NM2=MIN0(IM1-1,NM-IM1)
      D(IM1CT,NM)=1.E0
      IF(NM2.EQ.0) GOTO 10
      DO 20 NIT=1,NM2
      M2=L-NIT
      IM2=M2+LP1
      IF(M2.NE.LM1) GOTO 70
      T1=0.E0
      GOTO 30
70    T1=-SQRT(FLOAT((IM2+1)*(L-M2-1)))*D(IM1CT,IM2+2) 
30    D(IM1CT,IM2)=T1-(2.E0/STH)*(CTH*FLOAT(M2+1)-RM1)
     1    *D(IM1CT,IM2+1)
      D(IM1CT,IM2)=D(IM1CT,IM2)/SQRT(FLOAT(IM2*(L-M2)))
      TEMP=D(IM1CT,IM2)
      RMOD=ABS(TEMP)
      IF(RMOD.LT.BIG) GOTO 20
      IF(NIT.EQ.NM2) GOTO 20
      D(IM1CT,NIT+1)=DLBIG
      D(IM1CT,IM2)=D(IM1CT,IM2)/BIG
      D(IM1CT,IM2+1)=D(IM1CT,IM2+1)/BIG
20    CONTINUE
10    CONTINUE

C                              SET UP NORMALIZATION FOR RIGHTMOST COLUMN

      T1=FLOAT(2*L)*DLOGS
      IF(LMN.EQ.0) GOTO 720
      DO 710 I=1,LMN
      M1=I-L
      T1=DLOGF+0.5E0*ALOG10(FLOAT(LP1-M1)/FLOAT(L+M1))+T1
710   CONTINUE
720   D(1,1)=T1
      IF(NROW.EQ.1) GOTO 730
      DO 110 IM1CT=2,NROW
      M1=IM1CT-NMAXP1
110   D(IM1CT,1)=DLOGF+0.5E0*ALOG10(FLOAT(L-M1+1)/FLOAT(L+M1))
     1     +D(IM1CT-1,1)

730   SGN=-1.E0
      IF((LMN/2)*2.NE.LMN) SGN=1.E0

C                                       RENORMALIZE ROWS

      DO 120 IM1CT=1,NROW
      IM1=IM1CT+LMN
      SGN=-SGN
      CSUM=D(IM1CT,1)
      MULT=1
520   IF(ABS(CSUM).LT.DLBIG) GOTO 510
      MULT=MULT*2
      CSUM=0.5*CSUM
      GOTO 520
510   FAC=10.E0**CSUM
      SFAC=SMALL/FAC
      NM2=MIN0(IM1-1,NM-IM1)
      NM2P1=NM2+1
      DO 130 IM2=1,NM2P1
      IF((D(IM1CT,IM2+1).EQ.0.E0).OR.(IM2.GE.NM2)) GOTO 250
      CSUM=CSUM*FLOAT(MULT)+D(IM1CT,IM2+1)
      MULT=1
220   IF(ABS(CSUM).LT.DLBIG) GOTO 210
      MULT=MULT*2
      CSUM=0.5E0*CSUM
      GOTO 220
210   FAC=10.E0**CSUM
      SFAC=SMALL/FAC
250   IN2=NMP1-IM2
      DO 270 I=1,MULT
      TEMP=D(IM1CT,IN2)
      RMOD=ABS(TEMP)
      IF(RMOD.GT.SFAC) GOTO 260
      D(IM1CT,IN2)=0.E0
      GOTO 130
260   D(IM1CT,IN2)=D(IM1CT,IN2)*FAC
270   CONTINUE
      D(IM1CT,IN2)=SGN*D(IM1CT,IN2)
130   CONTINUE
120   CONTINUE

C                                        FILL REST OF MATRIX

400   IF(ISUP.GT.0) GOTO 410
      SGN=-1.E0
      IF((LMN/2)*2.NE.LMN) SGN=1.E0
      DO 420 IM1CT=1,NROW
      SGN=-SGN
      IM1=IM1CT+LMN
      NM2=MIN0(IM1,NMP1-IM1)
      DO 420 IN2=1,NM2
      IM2=NMP1-IN2
420   D(IM1CT,IN2)=SGN*D(IM1CT,IM2)
      DO 430 IM1CT=1,NROW
      IM1=IM1CT+LMN
      IN1=NMP1-IM1
      IN1CT=IN1-LMN
      SGN=-1.E0
      NM2=MIN0(IM1,IN1)
      DO 440 NIT=1,NM2
      SGN=-SGN
      IM2=1+NM2-NIT
      IN2=NMP1-IM2
      IM2CT=IM2-LMN
      IN2CT=IN2-LMN
      D(IN1CT,IN2)=SGN*D(IM1CT,IM2)
      IF(IN2CT.GT.NROW) GOTO 440
      D(IM2CT,IM1)=D(IN1CT,IN2)
      D(IN2CT,IN1)=D(IM1CT,IM2)
440   CONTINUE
430   CONTINUE
      RETURN
410   DO 450 IM1CT=1,NROW
      IM1=IM1CT+LMN
      IN1=NMP1-IM1
      IN1CT=IN1-LMN
      SGN=-1.E0
      NM2=MIN0(IM1,IN1)
      DO 460 NIT=1,NM2
      SGN=-SGN
      IM2=NM-NM2+NIT
      IN2=NMP1-IM2
      IM2CT=IM2-LMN
      IN2CT=IN2-LMN
      D(IN1CT,IN2)=SGN*D(IM1CT,IM2)
      IF(IM2CT.GT.NROW) GOTO 460
      D(IM2CT,IM1)=D(IN1CT,IN2)
      D(IN2CT,IN1)=D(IM1CT,IM2)
460   CONTINUE
450   CONTINUE
      RETURN
      END


C--------------------------------------------------------------

      subroutine bffi(iu,ibaf,nbyte,irec)
      
cfb   VAR DECLARATION

      implicit none
      
      integer ibaf(1),iu,j,irec,nbyte
c      dimension ibaf(1)
cfb

      read(iu,rec=irec) (ibaf(j),j=1,nbyte)


      return
      end

C------------------------  END ! ------------------------------
