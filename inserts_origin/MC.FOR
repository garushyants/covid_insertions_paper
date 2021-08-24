        program mc

        open(unit=3,file='weight.txt',status='old')
        fmon=0.
        all=0
        read(3,655)iiio
655     format(1x,i4)
777     read(3,655,end=888)iiir
        all=all+1
         if(iiir.ne.-1)then
        if(iiio.ge.iiir)fmon=fmon+1.
         endif
        goto 777
888     continue
        write(3,*)
        write(3,999)fmon/all
999     format('P = ',f5.2)
        end

**********************************************************
      subroutine endfil(n)
c выход на конец файла
1     read(n,*,end=2)
      goto 1
2     backspace n
      end

