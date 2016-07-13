       INCLUDE  'PARAMS.CMD'
       INCLUDE  'CHPARM.CMD'
       INCLUDE  'CNTROL.CMD'
       INCLUDE  'FILCON.CMD'
       INCLUDE  'SEGTAB.CMD'
       INCLUDE  'NETDEP.CMD'
       INCLUDE  'BALANC.CMD'
       INCLUDE  'LOCPTR.CMD'
       INCLUDE  'MSCAL.CMD'
       character*4 NAMAV,MSPEC,MRUNID
      integer,parameter::MXSPEC=40
      character*4  SPNAME(10,MXSPEC)
      integer,parameter::fmax=300
      integer Itmp(4)
      CHARACTER*100 NAM0(fmax) ! input/output file names

      integer,allocatable::ndate2(:),ndlast2(:),yyjjjhhBE(:,:)
      real,allocatable::ttime2(:),ttlast2(:)
      real,allocatable:: A1(:,:,:,:,:),tim(:)
      integer,allocatable:: yyjul(:),n(:)

      NUAV=41
      narg=iARGc ()
c      if(narg.lt.3) stop 'cbin_avrg File1 File2 File3... CombineFile'
      if(narg.gt.fmax+1) stop 'File too many'

      do i=1,narg
        call getarg(i,nam0(i))
      enddo

      do i=1,narg
        if(i.eq.narg)then
          open(i+10,file=trim(nam0(i)),
     +      form='unformatted',convert='BIG_ENDIAN',STATUS='unknown')
        else
          open(i+10,file=trim(nam0(i)),
     +      form='unformatted',convert='BIG_ENDIAN',STATUS='old')
        endif
      enddo
      nfile=narg-1
      iout=narg+10

      allocate(ndate2(nfile))
      allocate(ndlast2(nfile))
      allocate(ttime2(nfile))
      allocate(ttlast2(nfile))
      allocate(yyjjjhhBE(2,nfile))

      DO IRD=1,nfile
        READ (IRD+10) NAMAV, MRUNID, NOSEG, NOSPEC,
     +    NDATE2(ird), TTIME2(ird),
     $    NDLAST2(ird), TTLAST2(ird)
       enddo
        print*,NDATE2,TTIME2
        print*,NDLAS
       NDATE=NDATE2(1)
       NDLAST=NDLAST2(nfile)
       TTIME=TTIME2(1)
       TTLAST=TTLAST2(nfile)
      write(narg+10) NAMAV, MRUNID, NOSEG, NOSPEC, NDATE, TTIME,
     $  NDLAST, TTLAST

       do ird=1,nfile
C
C--REGION DESCRIPTION HEADER
        READ (IRD+10) XUTM, YUTM, NZONE, XORG, YORG, DELTAX, DELTAY,
     $    NOXG, NOYG, NOZ, NVLOW, NVUP, DZSURF, DZMINL, DZMINU

C--SEGMENT DESCRIPTION HEADER
	NOSEG=1
        READ (IRD+10) ((Itmp(j), J=1,4), I=1, NOSEG)
C--SPECIES DESCRIPTION HEADER
        READ (IRD+10) ((SPNAME(I,J), I=1,10), J=1,NOSPEC)
      ENDDO ! next IRD input file
      NXY=NOXG*NOYG
      NT=(NDLAST-NDATE)*24+(TTLAST-TTIME)
      allocate(A1(narg,NXY,NOZ,NOSPEC,NT))
      allocate(yyjul(NT+1),tim(NT+1),n(nt))
C
C         FIRST, WRITE TIME INTERVAL
C
      write(narg+10) XUTM, YUTM, NZONE, XORG, YORG, DELTAX, DELTAY,
     $  NOXG, NOYG, NOZ, NVLOW, NVUP, DZSURF, DZMINL, DZMINU
      write(narg+10)1,1,NOXG,NOYG
      write(narg+10)((SPNAME(I,J), I=1,10), J=1,NOSPEC)

      ird=1
      do it=1,NT
        n(it)=ird
        READ (ird+10,END=33) yyjul(it),tim(it), NDATE, TAVE
        !print*,ird, yyjul(it),tim(it), NDATE, TAVE
        DO  L=1,NOSPEC
          DO  K=1,NOZ
            READ(10+ird)ISEG,(MRUNID(I),I=1,10),(A1(ird,I,K,L,it),I=1,NXY)
          enddo
        enddo
        cycle
33      continue
        close(ird)
        if(ird.eq.narg-1)exit
        ird=ird+1
        READ (ird+10,END=33) yyjul(it),tim(it), NDATE, TAVE
        !print*,ird, yyjul(it),tim(it), NDATE, TAVE
        DO  L=1,NOSPEC
          DO  K=1,NOZ
            READ(10+ird)ISEG,(MRUNID(I),I=1,10),(A1(ird,I,K,L,it),I=1,NXY)
          enddo
        enddo
      enddo
      yyjul(NT+1)=yyjul(NT)
      tim(NT+1)=tim(NT)+1
      if(tim(NT+1).eq.24) then
        yyjul(NT+1)=yyjul(NT+1)+1
        tim(NT+1)=0
      endif
      ird=narg+10
      do it=1,NT
        write(ird)yyjul(it),tim(it),yyjul(it+1),tim(it+1)
        DO  L=1,NOSPEC
          DO  K=1,NOZ
            write(ird)ISEG,(SPNAME(I,L),I=1,10),A1(n(it),:,k,l,it)
          enddo !k
        enddo !l
      enddo !it
      end
