c this routine calculate maximun of average of fields, and OP AS A ARG FILE
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
      character*4 note(60)
      integer,allocatable::ndate2(:),ndlast2(:),jul(:)
      real,allocatable::ttime2(:),ttlast2(:),ttmax(:)
      real,allocatable:: A1(:,:,:,:),tm(:,:,:),mxt(:),tim(:)


      NUAV=41
      narg=iARGc ()
      if(narg.ne.1)stop  'input avrg_file '
      do i=1,narg
        call getarg(i,nam0(i))
      enddo
      narg=2
      nam0(2)=trim(nam0(1))//'M'
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
!       print*, ((SPNAME(I,J), I=1,10), J=1,NOSPEC)
        nt=0
      do 
        READ (ird+10,END=30,err=30)i1,t1,i2,t2
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
      allocate(tm(NXY,NOZ,NOSPEC))
      allocate(A1(NXY,NOZ,NOSPEC,NT))
      allocate(ttmax(NXY))
      allocate(mxt(NT),jul(NT),tim(NT))
      ird=1 
      do it=1,nt
        READ(ird+10,END=31,err=31)jul(it),tim(it),i2,t2
!     print*, i1,t1,i2,t2
        DO  L=1,NOSPEC
          DO  K=1,NOZ
            READ(10+ird,err=31)ISEG,(SPNAME(I,L),I=1,10),A1(:,k,l,it)
          ENDDO
        ENDDO
      ENDDO
      print*,'normal end'
      goto 32
31    print*,'wrong NT'  
32    continue
      i1=jul(1)
      t1=tim(1)
      open(1,file=trim(nam0(2))//'98.dat',status='unknown')
      write(iout) fname, note, NOSEG, NOSPEC, i1, t1,
     $ i2, t2
      write(iout) XUTM, YUTM, NZONE, XORG, YORG, DELTAX, DELTAY,
     $  NOXG, NOYG, NOZ, NVLOW, NVUP, DZSURF, DZMINL, DZMINU
      write(iout)1,1,NOXG,NOYG
      write(iout)((SPNAME(I,J), I=1,10), J=1,NOSPEC)
      mx98=max(1,int(0.02*nt)) !maximum 98%
!      DO  L=1,NOSPEC
!        DO  K=1,NOZ
!          DO  I=1,NXY
!            tm(I,K,L)=maxval(A1(I,K,L,:))
!          enddo !j
!        enddo !k
!      enddo !l
      print*, i1,t1,i2,t2
      write(iout) i1,t1,i2,t2
      DO  L=1,NOSPEC
      do I=1,10
        nam0(2)(i:i)=SPNAME(I,L)(1:1)
      enddo
       DO  K=1,NOZ
       do I=1,NXY
        do it=1,nt
          mxt(it)= jul(it)*100+tim(it)
        enddo
        CALL SORT98(nt, A1(I,K,L,:),mxt,v,tv)
        tm(I,k,L)=v
        ttmax(I)=tv
       enddo
!      CALL SORT98(nt*NXY,A1(:,K,L,:),cmax)
        write(iout)ISEG,(SPNAME(I,L),I=1,10),(tm(I,k,l),I=1,NXY)
       enddo !k
      write(*,'(I5,A10,2F10.3)')L,nam0(2)(1:10),maxval(tm(:,1,l))
      write(1,'(I5,A10,2F10.3)')L,nam0(2)(1:10),maxval(tm(:,1,l))
      call sortNXY(NXY,tm(:,1,L),ttmax)  
      do i=1,10
        ical=ttmax(i)/100
        call caldate(ical)
      write(1,*)ttmax(i),ical
      enddo
      enddo !L
      end
      SUBROUTINE SORT98(M,SC,TM,v,t)
      DIMENSION TM(M),SC(M),TT(M)
      data small_v/-9999/
      if(maxval(SC).eq.small_v) stop 'max=small_v'
      M98=amax1(1.,0.02*real(M))
      if(M98.gt.1)then
        do i=1,M98-1
          SC(maxloc(SC(1:M)))=small_v
        enddo
      endif
      v=maxval(SC(1:M))
      do i=1,M
        if(SC(i).eq.v) t=TM(i)
      enddo
      return
      end

      SUBROUTINE SORTNXY(M,SC,TM)
      DIMENSION TM(M),SC(M),t1(M),t2(M),t3(M)
      data small_v/-9999/
      if(maxval(SC).eq.small_v) stop 'max=small_v'
      t1=SC  
      t2=TM  
      do i=1,10
        v=maxval(t1(1:M))
        do k=1,M
          if(t1(k).eq.v)TM(i)=t2(k)
        enddo
        do j=i,M
          if(t2(j).eq.TM(i)) t1(j)=small_v
        enddo
      enddo
      return
      end
      

