! file1 over file2 with same time frame
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
      integer itmp(4)
      CHARACTER*60 NAM0(fmax),FNAME ! input/output file names

      integer,allocatable::ndate2(:),ndlast2(:)
      real,allocatable::ttime2(:),ttlast2(:),A1(:,:,:)

      NUAV=41
      narg=iARGc ()
c      if(narg.lt.3) stop 'cbin_avrg File1 File2 File3... CombineFile'
      if(narg.gt.fmax+1) stop 'File too many'

      do i=1,narg
        call getarg(i,nam0(i))
      enddo

      do i=1,narg
        do j=1,60
          if(nam0(i)(j:j).ne.' ')IEND=J
        enddo
        if(i.eq.narg)then
          open(i+10,file=trim(nam0(i)),
     +      form='unformatted',convert='BIG_ENDIAN',STATUS='unknown')
        else
          open(i+10,file=trim(nam0(i)),
     +      form='unformatted',convert='BIG_ENDIAN',STATUS='old')
        endif
        print*,nam0(i)
      enddo
      nfile=narg-1
      iout=narg+10

      allocate(ndate2(nfile))
      allocate(ndlast2(nfile))
      allocate(ttime2(nfile))
      allocate(ttlast2(nfile))

      DO IRD=1,nfile
        READ (IRD+10) NAMAV, MRUNID, NOSEG, NOSPEC,
     +    NDATE2(ird), TTIME2(ird),
     $    NDLAST2(ird), TTLAST2(ird)
       enddo

       do ird=1,nfile
C
C--REGION DESCRIPTION HEADER
        READ (IRD+10) XUTM, YUTM, NZONE, XORG, YORG, DELTAX, DELTAY,
     $    NOXG, NOYG, NOZ, NVLOW, NVUP, DZSURF, DZMINL, DZMINU
        NXY=NOXG*NOYG
      allocate(A1(NOSPEC,NOZ,NXY))
C
C--SEGMENT DESCRIPTION HEADER
        READ (IRD+10) (Itmp(j), J=1,4)
C
C--SPECIES DESCRIPTION HEADER
        READ (IRD+10) ((MSPEC(I,J), I=1,10), J=1,NOSPEC)
C
C         FIRST, WRITE TIME INTERVAL
C

      ENDDO ! next IRD input file
!check the dates and times
      DO IRD=1,nfile
        ifirst=1
        do while(.true.)
          READ (10+ird,end=22) NDATE1, TAVE1, NDATE, TAVE
       if(ifirst.eq.1) then 
          write(*,*)ird,   NDATE1, TAVE1, NDATE, TAVE
        ndate2(ird)=NDATE1
        ttime2(ird)=TAVE1
        ifirst=0
        endif
          DO  L=1,NOSPEC
            DO  K=1,NOZ
              READ(10+ird,end=21)
21            continue
            enddo !K
          enddo !L
      enddo
22     ndlast2(ird)=NDATE
       ttlast2(ird)=TAVE
          write(*,*)ird,   NDATE1, TAVE1, NDATE, TAVE
       rewind(ird+10)
      enddo
       ndate=minval(ndate2)
       ndlast=maxval(ndlast2)
       ttime=minval(ttime2)
       ttlast=maxval(ttlast2)
          write(*,*)0,   NDATE, TTIME, NDLAST, TTLAST
      beg=NDATE2(2)*100+TTIME2(2)

      DO IRD=1,nfile
        do ih=1,4
        READ (IRD+10) 
        enddo
      enddo
      write(narg+10) NAMAV, MRUNID, NOSEG, NOSPEC, NDATE, TTIME,
     $  NDLAST, TTLAST
      write(narg+10) XUTM, YUTM, NZONE, XORG, YORG, DELTAX, DELTAY,
     $  NOXG, NOYG, NOZ, NVLOW, NVUP, DZSURF, DZMINL, DZMINU
       write(narg+10)(Itmp(j) , J=1,4)
       write(narg+10)((MSPEC(I,J), I=1,10), J=1,NOSPEC)

      ird=1
      ird2=2
34      do while(.true.)
          READ (10+ird,end=33) NDATE1, TAVE1, NDATE, TAVE
          write(narg+10)NDATE1, TAVE1, NDATE, TAVE
          write(*,*)ird,NDATE1, TAVE1, NDATE, TAVE  !!
          datetime=NDATE1*100+TAVE1
          DO  L=1,NOSPEC
            DO  K=1,NOZ
         READ(10+ird,end=35)ISEG,(MSPEC(I,L),I=1,10),(A1(L,K,I),I=1,NXY)
35          write(narg+10)iseg,(mspec(i,l),i=1,10),(a1(L,K,i),i=1,nxy)
            enddo !K
          enddo !L
          if(ird.eq.1.and.datetime.ge.beg) then 
!skip the overlapped period
          READ (10+ird2) NDATE1, TAVE1, NDATE, TAVE
          write(*,*)0,   NDATE1, TAVE1, NDATE, TAVE
          DO  L=1,NOSPEC
            DO  K=1,NOZ
              READ(10+ird2)ISEG,(MSPEC(I,L),I=1,10),(A1(L,K,I),I=1,NXY)
            enddo !K
          enddo !L
          endif
        enddo !T
33      close(10+ird)
      if(ird.eq.1) then
        ird=2
        goto 34
      elseif(ird.eq.2) then
        stop 'normal stop'
      endif
      end
