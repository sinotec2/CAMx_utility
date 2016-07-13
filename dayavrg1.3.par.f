      program dayavrg
c-----v1.0-2013/11/23
c     **改寫自camx.com提供的bin2asc
c-----v1.3-2014/12/12
c     **搭配davrg1.1.f(PM2.5改為PM25,物種8種)

      integer,parameter::mxspec=200
      integer nx,ny,nz
      integer stdate,endate
      real,allocatable::conc(:,:,:),maxc(:,:),o38(:,:),pm(:,:,:)
      real,allocatable::conc2(:,:,:,:)
      integer,allocatable::yyjjjhh(:,:)
      real orgx,orgy,utmx,utmy,dx,dy
      real a(200,200),r(4)
      real SCR(40000)
      character*40 fn,rtfn
      character*10 ftype,spcname
      character*4 fname(10), spname(10,MXSPEC)
      character*4 note(60)
      character*5 oup(4)
      character*6 adate
      real XLONC,PHIC
      data XLONC/120.99/PHIC/23.61/
      iutmzon=(180+XLONC)/6+1
      call utmgeo(0,iutmzon,rx4,ry4,XLONC,PHIC)
      Xcent=rx4*1000.           !center point UTM-m
      Ycent=ry4*1000.           !center point UTM-m
      oup(1)='o3h  '
      oup(2)='o38h '
      oup(3)='pm10d'
      oup(4)='pm25d'

      narg=iARGc ()
      if(narg.ne.3)then
        if(narg.eq.1) then
          stdate=-1
          endate=999999
          call getarg(1,fn)
        else
          stop'dayavrg file [start_date(yymmdd) end_date(mmdd)]'
        endif
      else
        call getarg(1,fn)
        call getarg(2,adate)
        read(adate,'(3i2)')iy,im,id
        stdate=(iy*100+im)*100+id
        call juldate(stdate)
        call getarg(3,adate(1:4))
        read(adate(1:4),'(2i2)')im,id
        endate=(iy*100+im)*100+id
        call juldate(endate)
      endif

      do i=1,40
        if(fn(i:i).eq.'.')then
          rtfn(1:i-1)=fn(1:i-1)
          irt=i-1
          exit
        endif
      enddo

c-----讀入模擬值
      open (610,file=trim(fn),form='unformatted',status='old',
     +  convert='big_endian')
      read (610) fname,note,nseg,nspec,idate,begtim,jdate,endtim
      write(ftype,'(10a1)') (fname(i),i=1,10)
      read (610) rr,rr,iutm,xorg,yorg,dx,dy,noxg,noyg,
     &          nz,nzlo,nzup,rhts,rhtl,rhtu

      read (610) i1,j1,nx1,ny1
      read (610) ((spname(m,l),m=1,10),l=1,nspec)
      IO3=0
      IPc=0
      IPf=0
      do l=1,nspec
        do m=1,10
          spcname(m:m)=spname(m,l)(1:1)
        enddo
        if(spcname(1:10).eq.'O3        ')IO3=l
        if(spcname(1:10).eq.'PM10      ')IPc=l
        if(spcname(1:10).eq.'PM25      ')IPf=l
      enddo
      if(IO3*IPc*IPf.eq.0) stop 'spec not found !'
      nt=0
      do
        read(610,end=611)
        do kl=1,nz*nspec
          read(610)
        enddo        
        nt=nt+1
      enddo
611   rewind(610)
      do i=1,4 !read header
          read(610)
      enddo

      if(nspec .lt. 1) then
        write(*,*) 'no of species less than 1'
        stop
      endif

      nlay=nz
      nxy=noxg*noyg
      allocate(conc(nxy,nspec,nt))
      allocate(maxc(nxy,4))  !4個物種 o3max1h,o3max8h,pm10maxd,pm25maxd
      allocate(yyjjjhh(2,nt))  !2物種(pm10與pm25)

      idold=0
      do it=1,nt
        read(610) ndate1, tave1, ndate, tave
        yyjjjhh(1,it)=ndate1*100+tave1
        yyjjjhh(2,it)=ndate*100+tave
        do l=1,nspec
        read(610)ii,(spname(i,l),i=1,10),conc(:,l,it)
        if(nz.gt.1)then
          do k=2,nz
            read(610)
          enddo
        endif
        enddo
      enddo
      close(610)
      itb=1
      do it=1,nt
        jd=yyjjjhh(1,it)/100
        ih=mod(yyjjjhh(1,it),100)
        if(jd.ge.stdate.and.ih.lt.(24-16)) then
        itb=it
        exit
        endif
      enddo
      ite=nt
      do it=nt,1,-1
        jd=yyjjjhh(1,it)/100
        ih=mod(yyjjjhh(1,it),100)
        if(jd.le.endate.and.ih.gt.16) then
        ite=it
        exit
        endif
      enddo
      nday=(ite-itb)/24+1
      allocate(o38(nxy,nt))
      allocate(pm(nxy,nday,2))
      allocate(conc2(nxy,nspec,nday,0:23))
      id=1 
      do it=itb,ite+24,24
        if(it.gt.nt)exit
        do ih=0,23
          conc2(:,:,id,ih)=conc(:,:,it+ih)
        enddo
        id=id+1  
      enddo
      do i=1,nxy
      do id=1,nday
        pm(i,id,1)=sum(conc2(i,IPc,id,:))/24
        pm(i,id,2)=sum(conc2(i,IPf,id,:))/24
      enddo
      enddo
      nh=mod(yyjjjhh(1,ite),100)
      if(nh.lt.24.and.nh.gt.8)pm(:,nday,:)=pm(:,nday,:)*24/nh
      do i=1,nxy
        maxc(i,3)=maxval(pm(i,:,1))
        maxc(i,4)=maxval(pm(i,:,2))
      enddo
      do i=1,nxy
        maxc(i,1)=maxval(conc(i,IO3,:))
      enddo
      do it=4,nt-4
        do i=1,nxy
          sumo=0
          do ii=-3,4
            sumo=sumo+conc(i,IO3,it+ii)/8
          enddo
          o38(i,it)=sumo
        enddo
      enddo
      do i=1,nxy
        maxc(i,2)=maxval(o38(i,:))
      enddo
        r1=xorg+Xcent
        r2=xorg+real(noxg-1)*dx+Xcent
        r3=yorg+Ycent
        r4=yorg+real(noyg-1)*dy+Ycent
!       do i=1,4
!         r(i)=r(i)/1000.
!       enddo
!     call utmgeo(1,iutmzon,r(1),r(3),r1,r3)
!     call utmgeo(1,iutmzon,r(2),r(4),r2,r4)
c-----寫出.grd檔
      do ii=1,4
        open(10,file=trim(oup(ii))//'_'//rtfn(1:irt)//'.grd')
        do i=1,nxy  
            SCR(i)=maxc(i,ii)
        enddo
        call surfer(SCR,noxg,noyg,10,r1,r2,r3,r4)
        close(10)
      enddo
      write(*,*)'dayavrg finished.'
      end
c----------------------------------------------------
      subroutine surfer(SCR,ni,nj,n,xmin,xmax,ymin,ymax)  
!SCR:1 維陣列，ni(nj):i(j)網格數，n:寫出檔案之ifile，xmin:x座標的最小值
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



