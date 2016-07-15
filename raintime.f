c this routine calculate rain hr FLAG of CR fields, and OP AS A ARG FILE
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
      integer itmp(4),nti(2)
      character*4,allocatable:: SPNAME(:,:)
      CHARACTER*60 NAM0(fmax) ! input/output file names
      character*4 fname(10)
      character note(60)*4, names*10
      logical,allocatable::rain(:,:) !xy, and time 
      integer,allocatable::ndate2(:),ndlast2(:)
      real,allocatable::ttime2(:),ttlast2(:)
      real,allocatable:: A1(:,:,:,:),tim(:,:)
      real,allocatable:: dz1(:,:)
      real,allocatable:: tm(:,:,:)
      real,allocatable:: r(:,:)
      integer,allocatable:: yyjul(:,:)

      NUAV=41
      narg=iARGc ()
      if(narg.ne.1)stop  'RAIN_Time CR File'

      do i=1,narg
        call getarg(i,nam0(i))
      enddo
      narg=2
      nam0(2)=trim(nam0(1))//'F' !flags of rain
!     read the coord and level 1 length (m)
      nam0(3)=nam0(1)(1:len_trim(nam0(1))-2)//'3d'
          open(10,file=trim(nam0(3)),
     +      form='unformatted',convert='BIG_ENDIAN',STATUS='unknown')
      read(10)fname, note, NOSEG, NOSPEC
      allocate(SPNAME(10,NOSPEC))
      READ (10) XUTM, YUTM, NZONE, XORG, YORG, DELTAX, DELTAY,
     $    NOXG, NOYG, NOZ, NVLOW, NVUP, DZSURF, DZMINL, DZMINU
      READ (10) (Itmp(j), J=1,4)
      READ (10) ((SPNAME(I,J), I=1,10), J=1,NOSPEC)
      do I=1,10
        names(i:i)=SPNAME(I,1)
      enddo
      if(trim(names).ne.'ZGRID_M') stop 'wrong dz file'
      NXY=NOXG*NOYG
      allocate(dz1(NXY,NOZ))
      READ (10)i1,t1,i2,t2
      do k=1,NOZ
        READ(10)ISEG,(SPNAME(I,1),I=1,10),dz1(:,k)
      enddo
      close(10)
      deallocate(SPNAME)
!      print*,maxval(dz1),minval(dz1)

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
      allocate(rain(NXY,NT))
      allocate(r(NXY,NT))
      allocate(tm(NXY,NOZ,NOSPEC))
      allocate(yyjul(2,NT),tim(2,NT))
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
      rain=.false.
      r=0
	irain=0
	do J=1,NOSPEC
          do I=1,10
            names(I:I)=SPNAME(I,J)(1:1)
          enddo
          if(trim(names).eq.'RAINW_GpM3') then
	    irain=J
	    exit
	  endif
	enddo
      if(irain.eq.0)stop 'RAINW_GpM3 not found '
      do it=1,nt
      do i=1,NXY
        r(i,it)=sum(A1(i,:,irain,it)*dz1(i,:))/1000. !mm rain in it hour 
        if(r(i,it).ge.0.1) rain(i,it)=.true.
!       if(r(i,it).ge.0.1) r(i,it)=1.
      enddo  
      enddo  
      NOSPEC=2  
      NOZ=1  
      ird=narg+10
      names='RAIN_FALG'
      do I=1,10
        SPNAME(I,irain)(1:1)=names(I:I)
      enddo
      names='RAIN_mm'
      do I=1,10
        SPNAME(I,irain+1)(1:1)=names(I:I)
      enddo
      
      
      write(iout) fname, note, NOSEG, NOSPEC, NDATE, TTIME,
     $  NDLAST, TTLAST
      write(iout) XUTM, YUTM, NZONE, XORG, YORG, DELTAX, DELTAY,
     $  NOXG, NOYG,NOZ,NVLOW, NVUP, DZSURF, DZMINL, DZMINU
      write(iout)1,1,NOXG,NOYG
      write(iout)((SPNAME(I,J), I=1,10), J=irain,irain+1)
      do it=1,NT
        write(iout)yyjul(1,it),tim(1,it), yyjul(2,it),tim(2,it)
       write(iout)ISEG,(SPNAME(I,irain),I=1,10),rain(:,it)
       write(iout)ISEG,(SPNAME(I,irain+1),I=1,10),r(:,it)
      enddo !it
      stop
      end
