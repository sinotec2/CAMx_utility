
      integer Itmp(4)
      character*4,allocatable:: SPNAME(:,:)
      integer,parameter::fmax=300
      CHARACTER*60 NAM0(fmax) ! input/output file names
      character*4 fname(10)
      character note(60)*4,names*10
      logical,allocatable::LAND(:)
      real,allocatable:: A1(:,:,:,:),tim(:,:)
      integer,allocatable:: yyjul(:,:)

      narg=iARGc ()
      if(narg.ne.1)stop'landmask avrg_File'
      do i=1,narg
        call getarg(i,nam0(i))
      enddo
      nam0(2)=trim(nam0(1))//'K'
      IRD=11
      IWT=12
      do i=1,2
          open(i+10,file=trim(nam0(i)),
     +      form='unformatted',convert='BIG_ENDIAN',STATUS='unknown')
        print*,trim(nam0(i))
      enddo
      READ (IRD) fname, note, NOSEG, NOSPEC,
     +    NDATE, TTIME, NDLAST, TTLAST
      READ (IRD) XUTM, YUTM, NZONE, XORG, YORG, DELTAX, DELTAY,
     $    NOXG, NOYG, NOZ, NVLOW, NVUP, DZSURF, DZMINL, DZMINU
      READ (IRD) (Itmp(j), J=1,4)
      allocate(SPNAME(10,NOSPEC))
      READ (IRD) ((SPNAME(I,J), I=1,10), J=1,NOSPEC)
      nt=0
        do 
          READ (IRD,END=30,err=30)i1,t1,i2,t2
          DO  L=1,NOSPEC
            DO  K=1,NOZ
              READ (IRD,END=30,err=30)
            enddo !k
          enddo !l
        nt=nt+1
        enddo !it
C
C         FIRST, WRITE TIME INTERVAL
C
30    rewind(IRD)
      do i=1,4!skip the header
        READ(IRD)
      enddo
      NXY=NOXG*NOYG
      allocate(LAND(NXY))
      call landmask(NOXG,NOYG,LAND,xorg,yorg,deltax)
      allocate(A1(NXY,NOZ,NOSPEC,NT))
      allocate(yyjul(2,NT),tim(2,NT))
      do it=1,NT
        READ (IRD,END=33,err=33) yyjul(1,it),tim(1,it), yyjul(2,it),tim(2,it)
        DO  L=1,NOSPEC
          DO  K=1,NOZ
            READ(IRD,err=33)ISEG,(SPNAME(I,L),I=1,10),A1(:,k,l,it)
          enddo !k
        enddo !l
      enddo !it
      goto 34
33    stop 'end of file or error'
34    print*,'normal reading'        
      do I=1,NXY
       if(.not.LAND(I))A1(I,:,:,:)=0.
      enddo

      WRITE(IWT) fname, note, NOSEG, NOSPEC,
     +    NDATE, TTIME, NDLAST, TTLAST
      WRITE(IWT) XUTM, YUTM, NZONE, XORG, YORG, DELTAX, DELTAY,
     $    NOXG, NOYG, NOZ, NVLOW, NVUP, DZSURF, DZMINL, DZMINU
      WRITE(IWT) (Itmp(j), J=1,4)
      WRITE(IWT) ((SPNAME(I,J), I=1,10), J=1,NOSPEC)
      do it=1,NT
        WRITE(IWT) yyjul(1,it),tim(1,it), yyjul(2,it),tim(2,it)
        print*, yyjul(1,it),tim(1,it), yyjul(2,it),tim(2,it)
        DO  L=1,NOSPEC
          DO  K=1,NOZ
           WRITE(IWT)ISEG,(SPNAME(I,L),I=1,10),A1(:,k,l,it)
          enddo !k
        enddo !l
      enddo !it
      stop
      end 

      subroutine landmask(M,N,LAND,xorg,yorg,delta)
      logical LAND(M,N)
      LAND=.false.
      data XLONC/120.99/PHIC/23.61/
      iutmzon=(180+XLONC)/6+1
      call utmgeo(0,iutmzon,rx4,ry4,XLONC,PHIC)
      Xcent=rx4           !center point UTM-km
      Ycent=ry4           !center point UTM-km

C--READ THE LAND USE DATA
        open(1,file='/home/kuang/bin/dict.txt',STATUS='OLD')
        read(1,*)
        do
        read(1,*,end=102)i,j,K
        ix=((i-Xcent)*1000-xorg)/delta+1
        iy=((j-Ycent)*1000-yorg)/delta+1
        if((ix-1)*(ix-M).gt.0)cycle
        if((iy-1)*(iy-N).gt.0)cycle
        LAND(ix,iy)=.true.
        enddo
102     close(1)
!       open(1,file='check.dat',status='unknown')
!       do j=N,1,-1
!       write(1,'(100I1)')(LAND(i,j),i=1,M)
!       enddo
!       close(1)
        return

      end

