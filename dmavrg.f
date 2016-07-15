c this routine calculate daily mean of average of fields, and OP AS A ARG FILE
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
      integer,allocatable::ndate2(:),ndlast2(:)
      real,allocatable::ttime2(:),ttlast2(:)
      real,allocatable:: A1(:,:,:,:,:),tm(:)
      real,allocatable:: D1(:,:,:,:)

      NUAV=41
      narg=iARGc ()
      if(narg.ne.1)stop  'input avrg_file '
      do i=1,narg
        call getarg(i,nam0(i))
      enddo
      narg=2
      nam0(2)=trim(nam0(1))//'D'
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
      allocate(SPNAME(10,NOSPEC+1))
       ndate=minval(ndate2)
       ndlast=maxval(ndlast2)
       ttime=minval(ttime2)
       ttlast=maxval(ttlast2)
       do ird=1,nfile
C
C--REGION DESCRIPTION HEADER
        READ (IRD+10) XUTM, YUTM, NZONE, XORG, YORG, DELTAX, DELTAY,
     $    NOXG, NOYG, NOZ, NVLOW, NVUP, DZSURF, DZMINL, DZMINU
!      print*, XUTM, YUTM, NZONE, XORG, YORG, DELTAX, DELTAY,
!    $      NOXG, NOYG, NOZ, NVLOW, NVUP, DZSURF, DZMINL, DZMINU

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
      READ (IRD+10) ((SPNAME(I,J), I=1,10), J=1,NOSPEC)
      iO3=0
      do J=1,NOSPEC
        if(SPNAME(1,J)(1:1).eq.'O'.and.SPNAME(2,J)(1:1).eq.'3') then
          iO3=J
          exit
        endif
      enddo
      if(iO3.eq.0) then
       print*, ((SPNAME(I,J), I=1,10), J=1,NOSPEC)
       stop
      endif
      nt=0
      minj=99999
      maxj=-99999
      do 
        READ (ird+10,END=30,err=30)i1,t1,i2,t2
        if(i1.eq.0)i1=i2
        minj=min(minj,i1)
        maxj=max(maxj,i1)
!       print*,i1,t1,i2,t2
        DO  L=1,NOSPEC
          DO  K=1,NOZ
            READ(10+ird,err=30)
          enddo !k
        enddo !l
        nt=nt+1
      enddo !it
C
C         FIRST, WRITE TIME INTERVAL
C
30    rewind(IRD+10)
      do i=1,4
        read(IRD+10)
      enddo
      ENDDO ! next IRD input file
      NXY=NOXG*NOYG
      allocate(tm(0:16))
      allocate(A1(NXY,NOZ,NOSPEC,minj:maxj,0:23))
      allocate(D1(NXY,NOZ,NOSPEC+1,minj:maxj))
      ird=1 
      do it=1,nt
        READ (ird+10,END=31,err=31)i1,t1,i2,t2
        if(i1.eq.0)i1=i2
        ih=int(t1)
        DO  L=1,NOSPEC
          DO  K=1,NOZ
            READ(10+ird,err=31)ISEG,(SPNAME(I,L),I=1,10),A1(:,k,l,i1,ih)
          ENDDO
        ENDDO
      ENDDO
      print*,'normal end'
      goto 32
31    print*,'wrong NT'  
32    continue
      DO  L=1,NOSPEC
        do i=1,10
          names(i:i)=SPNAME(I,L)(1:1)
        enddo
        names=trim(names)//'D' 
        do i=1,10
          SPNAME(I,L)(1:1)=names(i:i)
        enddo
      enddo
      L=NOSPEC+1
        names(1:10)='O3_DE'//'     ' 
        do i=1,10
          SPNAME(I,L)(1:1)=names(i:i)
        enddo
      write(iout) fname, note, NOSEG, NOSPEC+1, minj, 0.,
     $ maxj, 23.
      write(narg+10) XUTM, YUTM, NZONE, XORG, YORG, DELTAX, DELTAY,
     $  NOXG, NOYG, NOZ, NVLOW, NVUP, DZSURF, DZMINL, DZMINU
      write(iout)1,1,NOXG,NOYG
      write(iout)((SPNAME(I,J), I=1,10), J=1,NOSPEC+1)
      D1=0.
        DO  K=1,NOZ
          DO  I=1,NXY
            DO J=minj,maxj
      DO  L=1,NOSPEC
        if(L.ne.iO3) then
              if(sum(A1(I,K,L,J,:),mask=A1(I,K,L,J,:).gt.0).gt.0)then
                nh=0
                suma=0
                do ih=0,23
                  if(A1(I,K,L,J,ih).le.0)cycle
                  suma=suma+A1(I,K,L,J,ih)
                  nh=nh+1
                enddo
                if(nh.gt.0)D1(I,K,L,J)=suma/nh
              endif
        else
            D1(I,K,L,J)=maxval(A1(I,K,L,J,10:18))
            do ih=0,16
              tm(ih)=0
              do ii=0,7
                tm(ih)=tm(ih)+A1(I,K,L,J,ih+ii)/8
              enddo
            enddo
            D1(I,K,NOSPEC+1,J)=maxval(tm(10:16))
        endif
            ENDDO
          enddo !j
        enddo !k
      enddo !l
      do jjj=minj,maxj
        if(sum(D1(:,:,:,jjj)).eq.0)cycle
        write(iout) jjj, 0., jjj, 23
      DO  L=1,NOSPEC+1
       DO  K=1,NOZ
        write(iout)ISEG,(SPNAME(I,L),I=1,10),(D1(I,k,l,jjj),I=1,NXY)
       enddo !k
      enddo !l
      enddo !it
      end
