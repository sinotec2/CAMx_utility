c davrg0 file1 file2 file1-file2 (no matter spec. name)
       INCLUDE  'PARAMS.CMD'
       INCLUDE  'CNTROL.CMD'
       INCLUDE  'FILCON.CMD'
       INCLUDE  'SEGTAB.CMD'
       INCLUDE  'NETDEP.CMD'
       INCLUDE  'BALANC.CMD'
       INCLUDE  'LOCPTR.CMD'
       INCLUDE  'MSCAL.CMD'

C     COMDECK.CHPARM
C-----------------------------------------------------------------------
C
C     ***  /CHPARM/ CONTAINS CHEMISTRY PARAMETERS
C


      COMMON /CHPARM/ NSMAX,
     $   LREAC(32), LSSIC(32), LSSBC(32),
     $     BDSL (32), BDSU (32), BDNL (32), BDNU (32), DVRES(32)
C
      COMMON /CHPARM/ NRMAX,
     $  RKN  (90), RK   (90), LPHOT(90), LTDEP(90),
     $  ACTEN(90), TREF (90), TRINV(90), LATDEP(40)
C
      COMMON /CHPARM/ NCMAX,
     $  MCOEF(10,10), CVAL(10)
C
      LOGICAL  LREAC, LSSIC, LSSBC
      LOGICAL LPHOT, LTDEP
C
      CHARACTER*60 NAM0(5),FNAME ! input/output file names
      integer,parameter::MXSPEC=40
      character*4,allocatable:: SPNAME(:,:)

      real,allocatable:: A1(:,:,:,:,:),tim(:,:)
      integer nt(2)
      integer,allocatable:: yyjul(:,:)

      NUAV=41
	narg=iARGc ()
	if(narg.ne.3) stop 'addAVRG File1 File2 (File1+File2)'
	do i=1,narg
		call getarg(i,nam0(i))
	enddo
	NAM0(4)='BIG_ENDIAN'

	do i=1,2
          open(i,file=trim(nam0(i)),
     +      form='unformatted',convert=NAM0(4),STATUS='old')
          print*,trim(nam0(i))
	enddo

          open(3,file=trim(nam0(3)),
     +      form='unformatted',convert=NAM0(4))
          print*,trim(nam0(3))

	DO IRD=1,2
      READ (IRD,END=998) NAMAV, MRUNID, NOSEG, NOSPEC, NDATE, TTIME,
     $  NDLAST, TTLAST
	IF(IRD.EQ.1)IBD1=NDATE*100+TTIME
	IF(IRD.EQ.1)IED1=NDLAST*100+TTLAST
	IF(IRD.EQ.2)IBD2=NDATE*100+TTIME
	IF(IRD.EQ.2)IED2=NDLAST*100+TTLAST
        enddo
      allocate(SPNAME(10,NOSPEC))
	DO IRD=1,2
C
C--REGION DESCRIPTION HEADER
      READ (IRD ) XUTM, YUTM, NZONE, XORG, YORG, DELTAX, DELTAY,
     $  NOXG, NOYG, NOZ, NVLOW, NVUP, DZSURF, DZMINL, DZMINU
C
C--SEGMENT DESCRIPTION HEADER
      READ (IRD ) (Idum    , J=1,4)
C
C--SPECIES DESCRIPTION HEADER
      READ (IRD ) ((SPNAME(I,J), I=1,10), J=1,NOSPEC)
C
C         FIRST, WRITE TIME INTERVAL
C
        nt(IRD)=0
        do 
          read(IRD,end=29)i1,t1,i2,t2
        print*,i1,t1
		DO  L=1,NOSPEC
		DO  K=1,NOZ
		 READ (ird,END=998,ERR=998)
		ENDDO !LOOP K
		ENDDO !LOOP SPECIES
        nt(IRD)=nt(IRD)+1
      enddo
29    rewind(IRD)
      ENDDO ! next IRD input file
      NXY=NOXG*NOYG
      if(nt(1).ne.nt(2)) print*,'nt different',nt(1),nt(2)
	IBD=MAX0(IBD1,IBD2)
	IED=MIN0(IED1,IED2)
	IF(IBD.GT.IED)	STOP 'IBD.GT.IED'
	NDATE =IBD/100
	TTIME =MOD(IBD,100)
	NDLAST=IED/100
	TTLAST=MOD(IED,100)
		write(3) NAMAV, MRUNID, NOSEG, NOSPEC, NDATE, TTIME,
     $  NDLAST, TTLAST
		write(3) XUTM, YUTM, NZONE, XORG, YORG, DELTAX, DELTAY,
     $  NOXG, NOYG, NOZ, NVLOW, NVUP, DZSURF, DZMINL, DZMINU
		write(3)1,1,NOXG,NOYG
		write(3)((SPNAME(I,J), I=1,10), J=1,NOSPEC)
      NTm=min(nt(1),nt(2))
      allocate(A1(3,NXY,NOZ,NOSPEC,NTm))
      allocate(yyjul(2,NTm+1),tim(2,NTm+1))
      do ird=1,2
        do i=1,4
   	READ (ird)
        enddo
101	READ (ird,END=998 ) NDATE1, TAVE1, NDATE, TAVE
	ID=NDATE1*100+TAVE1
	IF(ID.LT.IBD) THEN
		DO  L=1,NOSPEC
		DO  K=1,NOZ
		 READ (ird,END=998,ERR=997)
		ENDDO !LOOP K
		ENDDO !LOOP SPECIES
997		GOTO 101
	ENDIF
        backspace(ird)
	do it=1,NTm
	READ (ird,END=998 ) yyjul(1,it),tim(1,it),yyjul(2,it),tim(2,it)
!print*,ird,yyjul(it),tim(it), NDATE, TAVE
      DO  L=1,NOSPEC
      DO  K=1,NOZ
       READ(ird,END=994)ISEG,(SPNAME(I,L),I=1,10),A1(ird,:,k,l,it)
	ENDDO !LOOP K
	ENDDO !LOOP SPECIES
	enddo !it
        close(ird)
	enddo !ird
994    A1(3,:,:,:,:) = A1(1,:,:,:,:) + A1(2,:,:,:,:)
	ird=3
	do it=1,NTm
	write(ird)yyjul(1,it),tim(1,it),yyjul(2,it),tim(2,it)
!R	print*,ird,yyjul(it),tim(it),yyjul(it+1),tim(it+1)
      DO  L=1,NOSPEC
      DO  K=1,NOZ
       write(ird)ISEG,(SPNAME(I,L),I=1,10),A1(ird,:,k,l,it)
	ENDDO !LOOP K
	ENDDO !LOOP SPECIES
	enddo !it
        stop 'normal stop'
998     stop 'end of file'
        END
