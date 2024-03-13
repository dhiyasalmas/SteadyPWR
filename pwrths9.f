      program steadypwr
c  
c  program simulasi thermal hydraulik PWR
c  Model paling sederhana : analisis thermal hydraulik di jalur
c  pendingin
c
      common /g1/rpelet,rclin,rclout,rcool,qfuel,coflow,msr,msz
      common /g2/wzrfmt,wpufmt,porfp,baeffu
      common /g3/tavpel,tavcld,tavgap,tavcol
      common /sth1/xmass(1000),tempr(1000),press(1000),
     &       hgt(1000),presdp(1000),drhead(1000)
      common /sth2/tpinp,tpout,xdens(1000),xvel(1000),
     &       xxx(1000),zzz(1000),dlp(1000),pwrb(1000),
     &       dxx(1000),dzz(1000)
      common /sth3/maxb,gtot,psist,hcheam,graf
      common /sth4/maxc,maxm,maxp,tpc(50),denct(50),visct(50),conct(50)
     &           ,cpct(50),tpm(50),denmt(50),conmt(50),cpmt(50)
c reads input data from file name 'initherm.inp'
      call inp1
c initializes the variables and parameters of the system
      call init1
      call thermcoo1
      write(6,*)'selesai thermcool'
      call prdrp2
      write(6,*)'selesai prdrop'
      stop
      end
c
c
      subroutine thermcoo1
c
c  Menghitung temperatur di jalur pendingin
c
      common /sth1/xmass(1000),tempr(1000),press(1000),
     &       hgt(1000),presdp(1000),drhead(1000)
      common /sth2/tpinp,tpout,xdens(1000),xvel(1000),
     &       xxx(1000),zzz(1000),dlp(1000),pwrb(1000),
     &       dxx(1000),dzz(1000)
      common /sth3/maxb,gtot,psist,hcheam,graf
      common /sth4/maxc,maxm,maxp,tpc(50),denct(50),visct(50),conct(50)
     &           ,cpct(50),tpm(50),denmt(50),conmt(50),cpmt(50)
c
c    menghitung/update temperatur di seluruh jalur pendingin
c
      tempr(1)=tpinp
      do 50 i=1,maxb-2 
        tptmp=tempr(i)
        htp=entsbw(psist,tptmp)
        cpwtp=cpw(htp,tpinp)
        tempr(i+1)=tempr(i)+pwrb(i)/gtot/cpwtp
        write(6,*)'i+1:',i+1,tempr(i+1),tptmp,houttp,cpwtp
   50 continue
c
c  menghitung pressure drop grafitasi
c
      drhead1=0.0
      drhead2=0.0
      hcheam=10.0
      graf=10.0
c
c  bagian updtream
c
      do 70 i=1,14
        tptmp=tempr(i)
        htp=entsbw(psist,tptmp)
        dentp=denw(htp,tptmp)
        drhead1=drhead1+dentp*graf*dzz(i)
        if (i.eq.9) then
          drhead1=drhead1+dentp*graf*hcheam
        endif
   70 continue
c
c  bagian downstream
c
      do 90 i=15,27
        tptmp=tempr(i)
        htp=entsbw(psist,tptmp)
        dentp=denw(htp,tptmp)
        drhead2=drhead2+dentp*graf*dzz(i)
        if (i.eq.23) then
          drhead2=drhead2+dentp*graf*hcheam
        endif
   90 continue
      write(6,*)'driving head 1 dan 2 ',drhead1,drhead2
      return
      stop
      end 
c
c
c    Menghitung pressure drop akibat friksi
c
      subroutine prdrp2
      common /sth1/xmass(1000),tempr(1000),press(1000),
     &       hgt(1000),presdp(1000),drhead(1000)
      common /sth2/tpinp,tpout,xdens(1000),xvel(1000),
     &       xxx(1000),zzz(1000),dlp(1000),pwrb(1000),
     &       dxx(1000),dzz(1000)
      common /sth3/maxb,gtot,psist,hcheam,graf
      common /sth4/maxc,maxm,maxp,tpc(50),denct(50),visct(50),conct(50)
     &        ,cpct(50),tpm(50),denmt(50),conmt(50),cpmt(50)
      power=3.0E8
      vol=3.0
      area=2.0
      crat=0.35
      tptmp=tempr(4)
      write(6,*)'tptmp:',tptmp
      htp=entsbw(psist,tptmp)
      write(6,*)'htp:',htp
      cptp=cpw(htp,tptmp)
      write(6,*)'cptp:',cptp
      dentp=denw(htp,tptmp)
      write(6,*)'dentp:',dentp
      flowt=power/cptp/40.0/area/crat
      write(6,*)'flowt:',flowt
      drhead1=dentp*graf*dzz(4)
      p=0.0125
      d=0.010
      s=0.001
      write(6,*)'tptmp,htp,cptp,dentp,flowt,drhead1'
      write(6,*)tptmp,htp,cptp,dentp,flowt,drhead1
      write(6,*)'j,htp,cptp,de,denw1,visw1,re,fsm,pdrop'
      pdrop=0.0
      flowt1=flowt*crat*area
      dpipe=0.4
      apipe=3.14159*dpipe**2/4.0
c
c    sebelum teras reaktor
c
      write(6,*)'j,de,velo,re,fsm,dpdrop,pdrop'
      do 50 j=1,3
        tptmp=tempr(j)
        htp=entsbw(psist,tptmp)
        cptp=cpw(htp,tptmp)
        de=dpipe
        denw1=denw(htp,psist)
        visw1=visw(htp,psist)
        velo=flowt1/denw1/apipe
        re=denw1*velo*de/visw1
        fsm=fricpt(re)
        dpdrop=velo**2/(2.0*de)*fsm*denw1*(dxx(j)+dzz(j))
        pdrop=pdrop+dpdrop
        write(6,*)j,de,velo,re,fsm,dpdrop,pdrop
   50 continue
      rco=p/(3.14159)**0.5
      areaco=3.14159*(rco**2-d**2/4.0)
      perim=2*3.14159*(rco+d/2)
      de=4*areaco/perim
c
c    bagian teras reaktor
c
c     write(6,*)'pitch rco de:',p,rco,de
      do 100 j=4,7
        tptmp=tempr(j)
        htp=entsbw(psist,tptmp)
        cptp=cpw(htp,tptmp)
        denw1=denw(htp,psist)
        visw1=visw(htp,psist)
        velo=flowt/denw1
        re=denw1*velo*de/visw1
        fsm=fricpt(re)
c       pdrop=pdrop+flowt**2/(2.0*de)*fsm/denw1*dzz(j)
c       dpdrop=velo**2/(2.0*de)*fsm/denw1*dzz(j)
        dpdrop=velo**2/(2.0*de)*fsm*denw1*dzz(j)
        pdrop=pdrop+dpdrop
        write(6,*)j,de,velo,re,fsm,dpdrop,pdrop
  100 continue
c
c   antara teras dengan steam generator
c
      do 150 j=8,17
        tptmp=tempr(j)
        htp=entsbw(psist,tptmp)
        cptp=cpw(htp,tptmp)
        de=dpipe
        denw1=denw(htp,psist)
        visw1=visw(htp,psist)
        velo=flowt1/denw1/apipe
        re=denw1*velo*de/visw1
        fsm=fricpt(re)
        dpdrop=velo**2/(2.0*de)*fsm*denw1*(dxx(j)+dzz(j))
        pdrop=pdrop+dpdrop
        write(6,*)j,de,velo,re,fsm,dpdrop,pdrop
  150 continue
      rsgi=0.5
      rsgo=1.0
      areaco=3.14159*(rsgo**2-rsgi**2)
      perim=2*3.14159*(rsgo+rsgi)
      de=4*areaco/perim
      flowt2=flowt/3
c
c    antara steam generator dengan pompa
c
c     write(6,*)'pitch rco de:',p,rco,de
      do 170 j=18,21
        tptmp=tempr(j)
        htp=entsbw(psist,tptmp)
        cptp=cpw(htp,tptmp)
        denw1=denw(htp,psist)
        visw1=visw(htp,psist)
        velo=flowt2/denw1
        re=denw1*velo*de/visw1
        fsm=fricpt(re)
        dpdrop=velo**2/(2.0*de)*fsm*denw1*abs(dzz(j))
        pdrop=pdrop+dpdrop
        write(6,*)j,de,velo,re,fsm,dpdrop,pdrop
  170 continue
  110 return
      stop
      end 
c     
c     INITIALISASI
c
c     blok 1-3 pipa pendingin antara pompa utama dengan teras
c     blok 4-7 teras reaktor dibagi 4 blok
c     blok 8-17 pipa pendingin antara teras sampai SG
c     blok 18-21 SG dibagi 4 bagian
c     blok 22-27 pipa pendingin antara SG dengan pompa utama
c     blok 28 pompa utama
c
c
      subroutine init1
      common /sth1/xmass(1000),tempr(1000),press(1000),
     &       hgt(1000),presdp(1000),drhead(1000)
      common /sth2/tpinp,tpout,xdens(1000),xvel(1000),
     &       xxx(1000),zzz(1000),dlp(1000),pwrb(1000),
     &       dxx(1000),dzz(1000)
      common /sth3/maxb,gtot,psist,hcheam,graf
      common /sth4/maxc,maxm,maxp,tpc(50),denct(50),visct(50),conct(50)
     &           ,cpct(50),tpm(50),denmt(50),conmt(50),cpmt(50)
      power=3.0E8
      tpinp=300.0
      tpout=340.0
      hinlet=entsbw(psist,tpinp)
      houtl=entsbw(psist,tpout)
      cpwin=cpw(hinlet,tpinp)
      cpwout=cpw(houtl,tpout)
      gtot=power/(tpout-tpinp)/(cpwin+cpwout)*2.0
      write(6,*)'hinlet, houtlet: ',hinlet,houtl
      write(6,*)'cpinlet, cpoutlet: ',cpwin,cpwout
      maxb=28
      psist=2.0E7
      hsat=entw(psist)
      tpfu1=800.0
      write(6,*)'h sat :',hsat
c     cpfuel1=cpfni(tpfu1)
c     tpcld1=400.0
c     cpcld1=cpm(tpcld1)
c     write(6,*)'cpcld :',cpcld1
c     write(6,*)'cp fu nit :',cpfuel1
      xxx(1)=1.0
      xxx(2)=2.0
      xxx(3)=2.5
      xxx(4)=2.5
      xxx(5)=2.5
      xxx(6)=2.5
      xxx(7)=2.5
      xxx(8)=2.5
      xxx(9)=2.5
      xxx(10)=2.5
      xxx(11)=2.5
      xxx(12)=2.5
      xxx(13)=2.0
      xxx(14)=1.0
      xxx(15)=0
      xxx(16)=-1.0
      xxx(17)=-2.0
      xxx(18)=-2.5
      xxx(19)=-2.5
      xxx(20)=-2.5
      xxx(21)=-2.5
      xxx(22)=-2.5
      xxx(23)=-2.5
      xxx(24)=-2.5
      xxx(25)=-2.5
      xxx(26)=-2.0
      xxx(27)=-1.0
      xxx(28)=0.0
      dxx(1)=1.0
      dxx(2)=1.0
      dxx(3)=0.0
      dxx(4)=0.0
      dxx(5)=0.0
      dxx(6)=0.0
      dxx(7)=0.0
      dxx(8)=0.0
      dxx(9)=0.0
      dxx(10)=0.0
      dxx(11)=0.0
      dxx(12)=0.0
      dxx(13)=1.0
      dxx(14)=1.0
      dxx(15)=1.0
      dxx(16)=1.0
      dxx(17)=1.0
      dxx(18)=0.0 
      dxx(19)=0.0 
      dxx(20)=0.0 
      dxx(21)=0.0 
      dxx(22)=0.0 
      dxx(23)=0.0 
      dxx(24)=0.0 
      dxx(25)=0.0 
      dxx(26)=1.0
      dxx(27)=1.0
      dxx(28)=1.0
      pwrb(1)=0.0
      pwrb(2)=0.0
      pwrb(3)=0.0
      pwrb(4)=power/4
      pwrb(5)=power/4
      pwrb(6)=power/4
      pwrb(7)=power/4
      pwrb(8)=0.0
      pwrb(9)=0.0
      pwrb(10)=0.0
      pwrb(11)=0.0
      pwrb(12)=0.0
      pwrb(13)=0.0
      pwrb(14)=0.0
      pwrb(15)=0.0
      pwrb(16)=0.0
      pwrb(17)=0.0
      pwrb(18)=-power/4
      pwrb(19)=-power/4
      pwrb(20)=-power/4
      pwrb(21)=-power/4
      pwrb(22)=0.0
      pwrb(23)=0.0
      pwrb(24)=0.0
      pwrb(25)=0.0
      pwrb(26)=0.0
      pwrb(27)=0.0
      pwrb(28)=0.0
      zzz(1)=0.0
      zzz(2)=0.0
      zzz(3)=0.5
      zzz(4)=1.50
      zzz(5)=2.50
      zzz(6)=3.50
      zzz(7)=4.50
      zzz(8)=5.50
      zzz(9)=6.50
      zzz(10)=7.50
      zzz(11)=8.50
      zzz(12)=9.50
      zzz(13)=10.00
      zzz(14)=10.00
      zzz(15)=10.00
      zzz(16)=10.00
      zzz(17)=10.00
      zzz(18)=9.50
      zzz(19)=8.50
      zzz(20)=7.50
      zzz(21)=6.50
      zzz(22)=5.00
      zzz(23)=3.00
      zzz(24)=1.50
      zzz(25)=0.50
      zzz(26)=0.00
      zzz(27)=0.00
      zzz(28)=0.00
      dzz(1)=0.0
      dzz(2)=0.0
      dzz(3)=1.0
      dzz(4)=1.00
      dzz(5)=1.00
      dzz(6)=1.00
      dzz(7)=1.00
      dzz(8)=1.00
      dzz(9)=1.00
      dzz(10)=1.00
      dzz(11)=1.00
      dzz(12)=1.00
      dzz(13)=0.00
      dzz(14)=0.00
      dzz(15)=0.00
      dzz(16)=0.00
      dzz(17)=0.00
      dzz(18)=1.00
      dzz(19)=1.00
      dzz(20)=1.00
      dzz(21)=1.00
      dzz(22)=2.00
      dzz(23)=2.00
      dzz(24)=1.00
      dzz(25)=1.00
      dzz(26)=0.00
      dzz(27)=0.00
      dzz(28)=0.00
      return
      stop
      end
c
c  TETAPAN-TETAPAN UNTUK PWR
c
c
c     function conw(hw,pw)
c     xx=hw/5.815e+5
c     conw=0.57374+0.25361*xx-0.14547*xx**2+0.013875*xx**3
c     return
c     end
c
c     menghitung densitas air dengan masukan entalpi dan tekanan
c
      function denw(hw,pw)
      if (hw.le.6.513e5) then
       denw=(999.65+4.9737e-7*pw)+(-2.5847e-10+6.1767e-19*pw)*hw**2+
     &     (1.2696e-22-4.9223e-31*pw)*hw**4
      else
       denw=(1488.64+1.3389e-6*pw)+(1.4695e9+8.85736*pw)/(hw-3.20372e6-
     &      1.20483e-2*pw)
      endif
      return
      end
c
c
c
c     menghitung entalpi air dengan masukan tekanan pada kondisi
c     saturasi
c
      function entw(pw)
      entw=5.7474e+5+2.09206e-1*pw-2.8051e-8*pw**2+2.38098e-15*pw**3-
     &     1.0042e-22*pw**4+1.6587e-30*pw**5
      return
      end
c
c
c
c     menghitung entalpi air sub cooled dengan masukan tekanan  dan
c     temperatur
c
      function entsbw(pw,tw1)
      tw=tw1+273.15
      temst=temsw(pw)
      entst=entw(pw)
      cpst=cpw(entst,pw)
      ent1=entst+(tw-temst)*cpst
      tvtst=temw(ent1,pw)
      if (abs(tvtst-tw)/tw.lt.0.001) goto 50
      do 40 i=1,25
       ent1=ent1+(tw-temst)/(tvtst-temst)*(tw-tvtst)*cpst 
       tvtst=temw(ent1,pw)
       if (abs(tvtst-tw)/tw.lt.0.001) goto 50
   40 continue
      write(6,*)'entsbw do not convergent!'
      stop
   50 entsbw=ent1
      return
      end
c
c
c
c     menghitung temperatur air sub cooled dengan masukan tekanan  dan
c     entalpi
c
      function temw(hw,pw)
      temw=(2.7291e2-1.5954e-7*pw)+(2.3949e-4-5.1963e-13*pw)*hw+
     &      (5.9660e-12+1.2064e-18*pw)*hw**2-
     &      (1.3147e-17+5.6026e-25*pw)*hw**3
      return
      end
c
c
c
c     menghitung cp  air sub cooled dengan masukan tekanan  dan
c     entalpi
c
      function cpw(hw,pw)
      a0p=2.3949e-4-5.1963e-13*pw
      a1p=1.1932e-11+2.4127e-18*pw
      a2p=-3.9441e-17-1.6808e-24*pw
      cpw=1.0/(a0p+a1p*hw+a2p*hw**2)
      return
      end
c
c
c
c     menghitung temperatur saturasi  air dengan masukan tekanan 
c
      function temsw(pw)
      entw1=entw(pw)
      temsw=temw(entw1,pw)
      return
      end
c
c
c
c     menghitung viskositas air sub cooled dengan masukan tekanan  dan
c     entalpi
c
      function visw(hw,pw)
      if (hw.le.2.76e5) then
       xx=8.5813e-6*(hw-4.2659e4)
       ee=6.4845e-6*(hw-5.5359e4)
       visw=1.2995e-3-9.2640e-4*xx+3.8105e-4*xx**2-8.2194e-5*xx**3+
     &      7.0224e-6*xx**4-(-6.5959e-12+6.763e-12*ee-2.8883e-12*ee**2+
     &      4.4525e-13*ee**3)*(pw-6.8946e5)
       goto 50
      endif
      if (hw.lt.3.94e5) then
       e0h=1.4526e-3-6.9881e-9*hw+1.521e-14*hw**2-1.2303e-20*hw**3
       e1h=-3.8064e-11+3.9285e-16*hw-1.2586e-21*hw**2+1.2860e-27*hw**3 
       visw=e0h+e1h*(pw-6.8946e5)
      else
       zz=3.8921e-6*(hw-4.0147e5)
       visw=3.026e-4-1.8366e-4*zz+7.5671e-5*zz**2-1.6479e-5*zz**3+
     &      1.4165e-6*zz**4
      endif
   50 return
      end
c
c     subroutine untuk masukkan data2
c
      subroutine inp1
c     INCLUDE "combsg3.f"
      common /g1/rpelet,rclin,rclout,rcool,qfuel,coflow,msr,msz
      common /sth4/maxc,maxm,maxp,tpc(50),denct(50),visct(50),conct(50)
     &           ,cpct(50),tpm(50),denmt(50),conmt(50),cpmt(50)
      open(2,file='initherm.inp')
      read(2,*)
      read(2,*)maxc,maxm,maxp
      write(6,*)'maxc,maxm,maxp'
      write(6,*)maxc,maxm,maxp
      read(2,*)
      read(2,*)(tpc(i),i=1,maxc) 
      read(2,*)
      read(2,*)(denct(i),i=1,maxc) 
      read(2,*)
      read(2,*)(visct(i),i=1,maxc) 
      read(2,*)
      read(2,*)(conct(i),i=1,maxc) 
      read(2,*)
      read(2,*)(cpct(i),i=1,maxc) 
      read(2,*)
      read(2,*)(tpm(i),i=1,maxm) 
      read(2,*)
      read(2,*)(denmt(i),i=1,maxm) 
      read(2,*)
      read(2,*)(conmt(i),i=1,maxm) 
      read(2,*)
      read(2,*)(cpmt(i),i=1,maxm) 
      write(6,*)'(tpc(i),i=1,maxc)' 
      write(6,*)(tpc(i),i=1,maxc) 
      write(6,*)'(denct(i),i=1,maxc)' 
      write(6,*)(denct(i),i=1,maxc) 
      write(6,*)'(visct(i),i=1,maxc)' 
      write(6,*)(visct(i),i=1,maxc) 
      write(6,*)'(conct(i),i=1,maxc)' 
      write(6,*)(conct(i),i=1,maxc) 
      write(6,*)'(cpct(i),i=1,maxc)' 
      write(6,*)(cpct(i),i=1,maxc) 
      write(6,*)'(tpm(i),i=1,maxm)' 
      write(6,*)(tpm(i),i=1,maxm) 
      write(6,*)'(denmt(i),i=1,maxm)' 
      write(6,*)(denmt(i),i=1,maxm) 
      write(6,*)'(conmt(i),i=1,maxm)' 
      write(6,*)(conmt(i),i=1,maxm) 
      write(6,*)'(cpmt(i),i=1,maxm)' 
      write(6,*)(cpmt(i),i=1,maxm) 
      close(2)
      return
      end
c
c     faktor friksi laminer
c
      function fricpl(re)
      fricpl=16.0/re
      return
      end
c
c     faktor friksi turbulen
c
      function fricpt(re)
      fricpt=0.0791*re**(-0.25)
      return
      end
c
c
