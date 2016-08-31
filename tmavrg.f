c this routine calculate all time average of AVRG fields, and OP AS A ARG FILE
       INCLUDE  'PARAMS.CMD'
       INCLUDE  'CHPARM.CMD'
       INCLUDE  'CNTROL.CMD'
       INCLUDE  'FILCON.CMD'
       INCLUDE  'SEGTAB.CMD'
       INCLUDE  'NETDEP.CMD'
       INCLUDE  'BALANC.CMD'
       INCLUDE  'LOCPTR.CMD'
       INCLUDE  'MSCAL.CMD'
      integer,parameter::fmax=300
      integer,parameter::MXSPEC=40
      integer itmp(4),nti(2)
      character*4,allocatable:: SPNAME(:,:)
      CHARACTER*60 NAM0(fmax) ! input/output file names
      character*4 fname(10)
      character note(60)*4,names*10
      logical O3,GRD 
      integer,allocatable::ndate2(:),ndlast2(:)
      real,allocatable::ttime2(:),ttlast2(:)
      real,allocatable:: A1(:,:,:,:),tim(:,:)
      real,allocatable:: tm(:,:,:)
      integer,allocatable:: yyjul(:,:)
      real SCR(40000),r(4)
      NUAV=41
      narg=iARGc ()
      if(narg.lt.1.or.narg.gt.3)stop  
     + 'Time Meam avrg File Logical_O3,GRD'

      do i=1,narg
        call getarg(i,nam0(i))
      enddo
	if(narg.eq.1) then
	  O3=.false.
          GRD=.false.
	elseif(narg.eq.2)then
	  read(nam0(2),*)O3
          GRD=.false.
        else
	  read(nam0(3),*)GRD
	endif
      narg=2
      nam0(2)=trim(nam0(1))//'T'
      do i=1,narg
          open(i+10,file=trim(nam0(i)),
     +      form='unformatted',convert='BIG_ENDIAN',STATUS='unknown')
        print*,trim(nam0(i))
      enddo
      nfile=narg-1
      iout=narg+10

      allocate(ndate2(nfile))
      allocate(ndlast2(nfile))
      allocate(ttime2(nfile))
      allocate(ttlast2(nfile))
      DO IRD=1,nfile
        READ (IRD+10) fname, note, NOSEG, NOSPEC,
     +    NDATE2(ird), TTIME2(ird),
     $    NDLAST2(ird), TTLAST2(ird)
       enddo
      allocate(SPNAME(10,NOSPEC))
       ndate=minval(ndate2)
       ndlast=maxval(ndlast2)
       ttime=minval(ttime2)
       ttlast=maxval(ttlast2)
       numHr=(ndlast-ndate)*24+(ttlast-ttime)

       do ird=1,nfile
C
C--REGION DESCRIPTION HEADER
        READ (IRD+10) XUTM, YUTM, NZONE, XORG, YORG, DELTAX, DELTAY,
     $    NOXG, NOYG, NOZ, NVLOW, NVUP, DZSURF, DZMINL, DZMINU
!       print*, XUTM, YUTM, NZONE, XORG, YORG, DELTAX, DELTAY

!       if(ird.eq.1)then
!         write(*,*) XUTM, YUTM, NZONE, XORG, YORG, DELTAX, DELTAY,
!    $      NOXG, NOYG, NOZ, NVLOW, NVUP, DZSURF, DZMINL, DZMINU
!       endif
C
C--SEGMENT DESCRIPTION HEADER
        READ (IRD+10) (Itmp(j), J=1,4)
!       print*, (Itmp(j), J=1,4)
C
C--SPECIES DESCRIPTION HEADER
        nti(IRD)=0
        READ (IRD+10) ((SPNAME(I,J), I=1,10), J=1,NOSPEC)
!       print*, ((SPNAME(I,J), I=1,10), J=1,NOSPEC)
        do 
          READ (ird+10,END=30,err=30)i1,t1,i2,t2
!       print*,i1,t1,i2,t2
          DO  L=1,NOSPEC
            DO  K=1,NOZ
              READ (ird+10,END=30,err=30)
            enddo !k
          enddo !l
        nti(IRD)=nti(IRD)+1
        enddo !it
C
C         FIRST, WRITE TIME INTERVAL
C
30    rewind(IRD+10)
      do i=1,4!skip the header
              READ(10+ird)
      enddo
      ENDDO ! next IRD input file
      NXY=NOXG*NOYG
      NT=nti(1)!(NDLAST-NDATE)*24+(TTLAST-TTIME)
      allocate(A1(NXY,NOZ,NOSPEC,NT))
      allocate(tm(NXY,NOZ,NOSPEC))
      allocate(yyjul(2,NT),tim(2,NT))

      write(narg+10) fname, note, NOSEG, NOSPEC, NDATE, TTIME,
     $  NDATE, TTIME+1
      write(narg+10) XUTM, YUTM, NZONE, XORG, YORG, DELTAX, DELTAY,
     $  NOXG, NOYG, NOZ, NVLOW, NVUP, DZSURF, DZMINL, DZMINU
       write(narg+10)1,1,NOXG,NOYG
       write(narg+10)((SPNAME(I,J), I=1,10), J=1,NOSPEC)
	iO3=0
	do J=1,NOSPEC
	  if(SPNAME(1,J)(1:1).eq.'O'.and.SPNAME(2,J)(1:1).eq.'3')then
	    iO3=J
	    exit
	  endif
	enddo
	if(iO3.eq.0)O3=.false.
      do  ird=1,nfile
        icount=0
        do it=1,NT
          READ (ird+10,END=33,err=33) yyjul(1,it),tim(1,it), yyjul(2,it),tim(2,it)
          DO  L=1,NOSPEC
            DO  K=1,NOZ
              READ(10+ird,err=33)ISEG,(SPNAME(I,L),I=1,10),A1(:,k,l,it)
            enddo !k
          enddo !l
          icount=icount+1
        enddo !it
      enddo !ird
      goto 34
33    print*,'end of file'
34    if(NT.ne.icount) stop 'NT not right'
      A1=A1/icount
      DO  L=1,NOSPEC
        DO  K=1,NOZ
          DO  I=1,NXY
            tm(I,K,L)=sum(A1(I,K,L,:))
          enddo !j
        enddo !k
       enddo !l
       if(O3) then
       L=iO3
        DO  K=1,NOZ
          DO  I=1,NXY
            tm(I,K,L)=maxval(A1(I,K,L,:))*icount
          enddo !j
        enddo !k
        endif
      ird=narg+10
      write(ird) NDATE, TTIME, NDATE, TTIME+1
      DO  L=1,NOSPEC
       DO  K=1,NOZ
        write(ird)ISEG,(SPNAME(I,L),I=1,10),(tm(I,k,l),I=1,NXY)
       enddo !k
      enddo !l
      if(GRD)then
      data XLONC/120.99/PHIC/23.61/
      iutmzon=(180+XLONC)/6+1
      call utmgeo(0,iutmzon,rx4,ry4,XLONC,PHIC)
      Xcent=rx4*1000.           !center point UTM-m
      Ycent=ry4*1000.           !center point UTM-m
        r1=xorg+Xcent
        r2=xorg+real(noxg-1)*deltax+Xcent
        r3=yorg+Ycent
        r4=yorg+real(noyg-1)*deltay+Ycent
!       do i=1,4
!         r(i)=r(i)/1000.
!       enddo
!     call utmgeo(1,iutmzon,r(1),r(3),r1,r3)
!     call utmgeo(1,iutmzon,r(2),r(4),r2,r4)
      DO  L=1,NOSPEC
        do i=1,10
          names(i:i)=SPNAME(I,L)(1:1)
        enddo
        open(10,file=trim(nam0(1))//'_'//trim(names)//'.grd')
        do i=1,nxy
          SCR(i)=tm(i,1,L)
        enddo
        call surfer(SCR,noxg,noyg,10,r1,r2,r3,r4)
        close(10)
      enddo
      endif
      end
      subroutine surfer(SCR,ni,nj,n,xmin,xmax,ymin,ymax)
      real SCR(40000)
      real cmin,cmax

      CMIN=99999.
      CMAX=-99999.
      J1=1
      J2=ni*nj
      DO JJ=1,J2
          CMIN=AMIN1(CMIN,SCR(JJ))
          CMAX=AMAX1(CMAX,SCR(JJ))
      enddo

      WRITE(n,'(a4)')'DSAA'
      WRITE(n,*)ni,nj
      WRITE(n,*)xmin,xmax
      WRITE(n,*)ymin,ymax
      WRITE(n,*)CMIN,CMAX
      I1=NI/10+1
      I2=MOD(NI,10)
      IF(I2.EQ.0)I1=I1-1

      K1=J1
      DO J=1,NJ
        DO KK=1,I1
          K2=K1+9
          IF(KK.EQ.I1.AND.I2.NE.0) K2=K1+I2-1
          WRITE(N,'(1x,10(1PG13.5E3))')(SCR(I),I=K1,K2)
          K1=K2+1
        enddo
        WRITE(N,*)
      enddo
      close(n)
      RETURN
      END subroutine surfer

