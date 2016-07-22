c this routine calculate K-direction mean of AVRG fields, and OP AS A ARG FILE
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


      integer,allocatable::ndate2(:),ndlast2(:)
      real,allocatable::ttime2(:),ttlast2(:)
      real,allocatable:: A1(:,:,:,:),tim(:)
      integer,allocatable:: yyjul(:)

      NUAV=41
      narg=iARGc ()
      if(narg.ne.1) stop 'slim_avrg File1 (UAM avrg File)'

      do i=1,narg
        call getarg(i,nam0(i))
      enddo
      narg=2
      nam0(2)=trim(nam0(1))//'L'
      do i=1,narg
          open(i+10,file=trim(nam0(i)),
     +      form='unformatted',convert='BIG_ENDIAN',STATUS='unknown')
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
       ndate=minval(ndate2)
       ndlast=maxval(ndlast2)
       ttime=minval(ttime2)
       ttlast=maxval(ttlast2)
      allocate(SPNAME(10,NOSPEC))



       do ird=1,nfile
C
C--REGION DESCRIPTION HEADER
        READ (IRD+10) XUTM, YUTM, NZONE, XORG, YORG, DELTAX, DELTAY,
     $    NOXG, NOYG, NOZ, NVLOW, NVUP, DZSURF, DZMINL, DZMINU

!        if(ird.eq.1)then
!          write(*,*) XUTM, YUTM, NZONE, XORG, YORG, DELTAX, DELTAY,
!     $      NOXG, NOYG, NOZ, NVLOW, NVUP, DZSURF, DZMINL, DZMINU
!        endif
C
C--SEGMENT DESCRIPTION HEADER
        READ (IRD+10) ((Itmp(j), J=1,4), I=1, NOSEG)
c        if(ird.eq.1)then
c          write(iout)((Idum, J=1,4), I=1, NOSEG)
c        endif
C
C--SPECIES DESCRIPTION HEADER
        READ (IRD+10) ((SPNAME(I,J), I=1,10), J=1,NOSPEC)
c        if(ird.eq.1)then
c          write(iout) ((MSPEC(I,J), I=1,10), J=1,NOSPEC)
c        endif
      nti(IRD)=0
        do 
          READ (ird+10,END=30)i1,t1,i2,t2
!       print*,i1,t1,i2,t2
          DO  L=1,NOSPEC*NOZ
              READ(10+ird)
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
C
C         FIRST, WRITE TIME INTERVAL
C
      allocate(A1(NXY,NOZ,NOSPEC,NT))
      allocate(yyjul(NT+1),tim(NT+1))
      ird=1
        do it=1,NT
        READ (ird+10,END=33) yyjul(it),tim(it),yyjul(it+1),tim(it+1) 
          DO  L=1,NOSPEC
            DO  K=1,NOZ
          READ(10+ird)ISEG,(SPNAME(I,L),I=1,10),(A1(I,k,l,it),I=1,NXY)
            enddo !k
          enddo !l
        enddo !it
33    continue
      nt_actual=nt
      if(it.le.nt)then !if final, it must be nt+1
        nt_actual=it-2
        NDLAST=yyjul(nt_actual)
        TTLAST=tim(nt_actual)
        endif
 
      write(narg+10) fname, note, NOSEG, NOSPEC, NDATE, TTIME,
     $  NDLAST, TTLAST
      write(narg+10) XUTM, YUTM, NZONE, XORG, YORG, DELTAX, DELTAY,
     $  NOXG, NOYG, 1, NVLOW, NVUP, DZSURF, DZMINL, DZMINU
       write(narg+10)1,1,NOXG,NOYG
       write(narg+10)((SPNAME(I,J), I=1,10), J=1,NOSPEC)
	ird=narg+10
        do it=1,NT_actual
        write(ird)yyjul(it),tim(it),yyjul(it+1),tim(it+1)
          DO  L=1,NOSPEC
            DO  K=1,1
               write(ird)ISEG,(SPNAME(I,L),I=1,10),
     +         (sum(A1(I,:,l,it))/NOZ,I=1,NXY)
            enddo !k
          enddo !l
        enddo !it
      end
