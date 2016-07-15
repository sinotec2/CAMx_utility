c this routine filter the depn file and transfer to deposition results(S/N/T,
dd/wd/td(kg/ha/month, year), lc(10^(-6)eq/liter)

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
      integer itmp(4),nti(2),ISN(2,4)
      character*4,allocatable:: SPNAME(:,:)
      CHARACTER*60 NAM0(fmax) ! input/output file names
      character*4 fname(10)
      character*4 note(60),var(4),PSN(3)
      character*10 names
      real unit(5,5)

      integer,allocatable::ndate2(:),ndlast2(:)
      real,allocatable::ttime2(:),ttlast2(:)
      real,allocatable:: A1(:,:,:,:),B(:,:,:,:),tim(:)
      integer,allocatable:: yyjul(:)

      NUAV=41
      narg=iARGc ()
      if(narg.ne.1) stop 'depo depn-File resultD '

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
      data var/'DD','WD','LC','TD'/PSN/'PSO4','PNO3','SnN'/
      var(:)(1:3)='_'//var(:)(1:2)
      ISN=0
      do mm=1,2
        do J=1,NOSPEC
          do I=1,10
            names(I:I)=SPNAME(I,J)(1:1)
          enddo
          if(names(1:4).ne.PSN(mm))cycle
          do K=1,3
            if(trim(names).eq.PSN(mm)//trim(var(K)))then
              ISN(mm,K)=J
              exit
            endif
          enddo
        enddo !L
      enddo
!     do K=1,3
!     print*,(ISN(mm,K),(SPNAME(I,ISN(mm,K)),I=1,10),mm=1,2)
!     enddo
       
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
        SPNAME=''
        do K=1,4
        do mm=1,3
        J=(K-1)*3+mm
        names(1:7)=PSN(mm)//trim(var(K))
          do I=1,10
            SPNAME(I,J)(1:1)=names(I:I)
          enddo
        enddo
        enddo
        NOSPEC=J
!     print*, ((SPNAME(I,J), I=1,10), J=1,NOSPEC)
      allocate(B(NXY,NOZ,NOSPEC,NT))
        unit=0.001 !g/ha ->kg/ha
        unit(1,3)= 1000*1000/96*2 ! g/l ->micro eq/l
        unit(2,3)= 1000*1000/62*1 !
        do K=1,3
          do mm=1,2 ! S and N
            J=(K-1)*3+mm
            B(:,:,J,:)=A1(:,:,ISN(mm,K),:)*unit(mm,K)
          enddo
          mm=3!SnN
          J=(K-1)*3+mm
          Js=(K-1)*3+1
          Jn=(K-1)*3+2
          B(:,:,J,:)=B(:,:,Js,:)+B(:,:,Jn,:)
        enddo
        K=4
        do mm=1,3 ! S and N
          Jd=(1-1)*3+mm
          Jw=(2-1)*3+mm
          J=(K-1)*3+mm
          B(:,:,J,:)=B(:,:,Jd,:)+B(:,:,Jw,:)
        enddo

      write(narg+10) fname, note, NOSEG, NOSPEC, NDATE, TTIME,
     $  NDLAST, TTLAST
      write(narg+10) XUTM, YUTM, NZONE, XORG, YORG, DELTAX, DELTAY,
     $  NOXG, NOYG, NOZ, NVLOW, NVUP, DZSURF, DZMINL, DZMINU
       write(narg+10)1,1,NOXG,NOYG
       write(narg+10)((SPNAME(I,J), I=1,10), J=1,NOSPEC)
	ird=narg+10
        do it=1,NT_actual
        write(ird)yyjul(it),tim(it),yyjul(it+1),tim(it+1)
          DO  L=1,NOSPEC
            DO  K=1,NOZ
               write(ird)ISEG,(SPNAME(I,L),I=1,10),
     +         (B(I,K,l,it),I=1,NXY)
            enddo !k
          enddo !l
        enddo !it
      end
