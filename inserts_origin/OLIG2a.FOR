        program oligo
c-----------
c universal
c-----------
c DBfile,mism1(ins only),"C","R"
c------------------------------------------------
        character str*35000,seq(300000000),seq1(300000000)
        character s1*50,s2*2,olig*30,oligs*25,seqnam(80),seqnam1(80)
        CHARACTER bas(5),basi(5),mm,s3,s4,seqins(500),seqins1(500)
        integer pos(2,1000)
        DATA BAS/'G', 'A', 'C', 'T', 'N'/
        DATA BASI/'C', 'T', 'G', 'A', 'N'/

        open(unit=1,file='ins_ext.seq',status='old')
2721    read(1,272)seqins
        if(seqins(1).eq.'>')goto 2721
        ins1=0
        do i=1,500
        if(seqins(i).eq.' ')ins1=i+1
        if(seqins(i).eq.' ')ins1a=i
        if(seqins(i).eq.' ')goto 3331
        enddo
3331    ins2=0
        do ii=i+2,500
        if((seqins(ii).eq.' ').and.(ins2.eq.0))ins2a=ii-2
        if((seqins(ii).eq.' ').and.(ins2.eq.0))ins2=ii-1
        if(seqins(ii).eq.' ')goto 3332
        enddo

3332    iseqin=0
        do i=1,500
         if(seqins(i).ne.' ')then
        iseqin=iseqin+1
        seqins1(iseqin)=seqins(i)
         endif
        enddo
        close(unit=1)

        ipos=0

c-------------------------
        length=ins2-ins1+1
c-------------------------

        do i=1,length
        oligs(i:i)=seqins(ins1+i-1)
        enddo
        do i=1,5
        olig(i:i)=seqins(ins1+i-7)
        enddo
        do i=1,length
        olig(i+5:i+5)=seqins(ins1+i-1)
        enddo

        i=igetarg(1,s1)
        i=igetarg(2,s2)
        mmis1=0
        if(s2.eq.'1 ')mmis1=1
        if(s2.eq.'2 ')mmis1=2
        if(s2.eq.'3 ')mmis1=3
        if(s2.eq.'4 ')mmis1=4
        if(s2.eq.'5 ')mmis1=5
        if(s2.eq.'6 ')mmis1=6
        if(s2.eq.'7 ')mmis1=7
        if(s2.eq.'8 ')mmis1=8
        if(s2.eq.'9 ')mmis1=9
        if(s2.eq.'10')mmis1=10
        if(s2.eq.'11')mmis1=11
        if(s2.eq.'12')mmis1=12
        mmis2=100
        i=igetarg(3,s3)
        icomp=0
        if((s3.eq.'c').or.(s3.eq.'C'))icomp=1
        i=igetarg(4,s4)
        iran=0
        if((s4.eq.'r').or.(s4.eq.'R'))iran=1

      if(iran.eq.0)GOTO 288
c random sampling
      open(unit=1,file='ran.tmp',status='unknown')
      iii=0
      read(1,8366,err=8365)iii
8365  iii=iii+111
      call amrset(iii)
      close(unit=1)
      flen=length-1
      do i=1,length
      ii1=ran(i1,i2)*(flen)+1
      ii2=ran(i1,i2)*(flen)+1
      mm=oligs(ii1:ii1)
      oligs(ii1:ii1)=oligs(ii2:ii2)
      oligs(ii2:ii2)=mm
      enddo
      do i=1,length
      olig(i+5:i+5)=oligs(i:i)
      enddo
      open(unit=1,file='ran.tmp',status='unknown')
      iii=30000.*ran(i1,i2)
      write(1,8366)iii
8366  format(i5)
      close(unit=1)
c random sampling end
288   continue

        open(unit=2,file='out.txt',status='unknown')

        ichro=0
        iend=0
        iseq=0
c chromosome
        open(unit=1,file=s1,status='old',recl=35000)
        read(1,272)seqnam
272     format(500a1)
10      read(1,2,end=12)str
        if(str(1:1).eq.'>')goto 112
2       format(a35000)
        do i=1,35000
        if(str(i:i).eq.'a')str(i:i)='A'
        if(str(i:i).eq.'t')str(i:i)='T'
        if(str(i:i).eq.'g')str(i:i)='G'
        if(str(i:i).eq.'c')str(i:i)='C'
        if(str(i:i).ne.' ')then
        iseq=iseq+1
        seq(iseq)=str(i:i)
        endif
        enddo
        goto 10

12      close(unit=1)
        iend=1

112      do i=1,80
         seqnam1(i)=str(i:i)
         enddo
        ichro=ichro+1
      do i=1,iseq
      if(seq(i).eq.'a')seq(i)='A'
      if(seq(i).eq.'t')seq(i)='T'
      if(seq(i).eq.'g')seq(i)='G'
      if(seq(i).eq.'c')seq(i)='C'
      enddo

c complement
        do ii=1,iseq
        seq1(ii)=seq(ii)
        enddo
C---------------------------------
      j=iseq+1
      do 826 i=1,iseq/2
      j=j-1
      mm=seq1(i)
      seq1(i)=seq1(j)
826   seq1(j)=mm
      do 827 i=1,iseq
      do 828 j=1,5
828   if(seq1(i).eq.bas(j))goto 8271
      seq1(i)='N'
      goto 827
8271  seq1(i)=basi(j) 
827   continue
c---------------------------------

c      do i=1,75
c      if((seqnam(i).eq.'r').and.(seqnam(i+1).eq.'e').and.
c     *(seqnam(i+2).eq.' ').and.(seqnam(i+3).eq.'a').and.
c     *(seqnam(i+4).eq.'c').and.(seqnam(i+5).eq.'u'))goto 1010
c      enddo

c DIR
c i = start of 5upstr

        do 222 i=1,iseq
        mism=0
        do j=6,length+5
         if(seq(i+j-1).ne.oligs(j-5:j-5))mism=mism+1
         if(mism.gt.mmis1)goto 222
        enddo
        mism2=0
        do j=1,length+5
         if(seq(i+j-1).ne.olig(j:j))mism2=mism2+1
         if(mism2.gt.mmis2)goto 222
        enddo
c wider context +-10
        mmmis=0
        iii=0
        do ij=ins1a-10,ins2a+10
        iii=iii+1
        if(seqins1(ij).ne.seq(i+iii-11+5))mmmis=mmmis+1
c        print *,seqins1(ij),seq(i+iii-11+5)
        enddo
c        print *,ins1a,ins2a,mmmis
        if(mmmis.le.10)goto 222
c===================
        kkk1=i-140
        if(kkk1.lt.1)kkk1=1
        kkk2=i+151
        if(kkk2.gt.iseq)kkk2=iseq
         if(iran.eq.0)then
        print 3491,i,mism,' d ',(seqnam(iii),iii=1,60)
        print 221,(seq(j),j=kkk1,kkk2)
221     format(70a1)
        print *
         endif

        write(2,3491)i,mism,' d ',(seqnam(iii),iii=1,60)
        write(2,7771)(olig(iii:iii),iii=6,length+5)
        write(2,7771)(seq(iii),iii=i+5,i+length+4)
7771    format(80a1)
        write(2,*)
        ipos=ipos+1
        pos(1,ipos)=mism
        pos(2,ipos)=i

222     continue

c COMP
             if(icomp.eq.1)then
        do 333 i=1,iseq
        mism=0
        do j=6,length+5
         if(seq1(i+j-1).ne.oligs(j-5:j-5))mism=mism+1
         if(mism.gt.mmis1)goto 333
        enddo
        mism2=0
        do j=1,length+5
         if(seq1(i+j-1).ne.olig(j:j))mism2=mism2+1
         if(mism2.gt.mmis2)goto 333
        enddo
c wider context +-10
        mmmis=0
        iii=0
        do ij=ins1a-10,ins2a+10
        iii=iii+1
        if(seqins1(ij).ne.seq1(i+iii-11+5))mmmis=mmmis+1
c        print *,seqins1(ij),seq1(i+iii-11+5)
        enddo
c        print *,ins1a,ins2a,mmmis
        if(mmmis.le.10)goto 333
c===================
        kkk1=i-140
        if(kkk1.lt.1)kkk1=1
        kkk2=i+151
        if(kkk2.gt.iseq)kkk2=iseq
         if(iran.eq.0)then
        print 3491,iseq-i+1,mism,' c ',(seqnam(iii),iii=1,60)
3491    format('>',i7,i2,a3,60a1)
        print 221,(seq1(j),j=kkk1,kkk2)
        print *
         endif

        write(2,3491)iseq-i+1,mism,' c ',(seqnam(iii),iii=1,60)
        write(2,7771)(olig(iii:iii),iii=6,length+5)
        write(2,7771)(seq1(iii),iii=i+5,i+length+4)
        write(2,*)
        ipos=ipos+1
        pos(1,ipos)=mism
        pos(2,ipos)=iseq-i+1

333     continue
             endif

c find max weight

1010    iseq=0
        if(iend.eq.1)goto 999

         do i=1,80
         seqnam(i)=seqnam1(i)
         enddo

        goto 10

999     open(unit=3,file='weight.txt',status='unknown')
        call endfil(3)
           if(ipos.ne.0)then
        iii=1
         if(ipos.gt.1)then
        do i=2,ipos  
        if(pos(1,i).lt.pos(1,iii))iii=i
        enddo
         endif
        if(iran.eq.1)then
        write(3,655)pos(1,iii)
        else
        write(3,656)pos(1,iii)
        endif
           else
        iii=-1
        if(iran.eq.1)then
        write(3,655)iii
655     format('R',i4)
        else
        write(3,655)iii
656     format('O',i4)
        endif
           endif

        end

**********************************************************
      function ran(i1,i2)
c     Portable Random number generator by George Marsaglia
      double precision u(97), c, cd, cm
      integer i97, j97
      logical set
      common /raset1/ u, c, cd, cm, i97, j97, set
      double precision  uni

      if ( .not. set ) then
           write(6,'(a)') 'amrand not initd'
           stop
      endif

      uni = U(I97) - U(J97)
      if ( uni .lt. 0.d0 ) uni = uni + 1.d0
      U(I97) = uni
      I97 = I97 - 1
      if (I97 .eq. 0) I97 = 97
      J97 = J97 - 1
      if (J97 .eq. 0) J97 = 97
      C = C - CD
      if ( C .lt. 0.d0 ) C = C + CM
      uni = uni - C
      if ( uni .lt. 0.d0 ) uni = uni + 1.d0
      ran = uni

      end

      subroutine amrset(iseed)
      integer iseed
      double precision  u(97), c, cd, cm
      integer i97, j97
      logical set
      common /raset1/ u, c, cd, cm, i97, j97, set
      integer is1, is2
      integer IS1MAX
      integer IS2MAX
      integer i,j,k,l,m
      real s,t
      integer ii, jj
      data IS1MAX, IS2MAX /31328, 30081/
      data set /.false./

      is1 = max((iseed / IS2MAX)+1, 1)
      is1 = min(is1, IS1MAX)
      is2 = max(1, mod(iseed, IS2MAX)+1)
      is2 = min(is2, IS2MAX)
      i = mod(is1/177, 177) + 2
      j = mod(is1    , 177) + 2
      k = mod(is2/169, 178) + 1
      l = mod(is2    , 169)
      do 200 ii = 1, 97
           s = 0.0d0
           t = 0.5d0
           do 100 jj = 1, 24
                m = mod(mod(i*j, 179)*k, 179)
                i = j
                j = k
                k = m
                l = mod(53*l+1, 169)
                if (mod(l*m, 64) .ge. 32) s = s + t
                t = 0.5d0 * t
  100      continue
           u(ii) = s
  200 continue

      c  = 362436.d0   / 16777216.d0
      cd = 7654321.d0  / 16777216.d0
      cm = 16777213.d0 / 16777216.d0

      i97 = 97
      j97 = 33

      set = .true.
      end
**********************************************************
      subroutine endfil(n)
c выход на конец файла
1     read(n,*,end=2)
      goto 1
2     backspace n
      end

