c 
c 
c ******************************************************************
c 
c
      integer function  opfrac(iat,ion,t,ane,pf,fra)
c     =======================================
c
c     Partition functions and ionization fractions
c     computed from the Opacity Project tables.
c     The routine returns values only fro the species
c     coinsidered by the Opacity Project, i.e.,
c
c     H,He,C,N,O,Ne,Na,Mg,Al,Si,S,Ar,Ca,Cr.Mn,Fe, and Ni
c
c     Input:
c      IAT   - atomic number
c      ION   - ionic charge (=1 for neutrals, =2 for once ionized, etc)
c      T     - temperature
c      ANE   - electron density
c
c     Output:
c      U     - partition function
c      FRA   - fraction of the given ion
c
c
      implicit real*8 (a-h,o-z)
      parameter (mtemp =100,
     *           melec = 60,
     *           mion1 = 30,
     *           mion2 = 32,
     *           mdat  = 17,
     *           mstag = 258)
      data ierr / 0 /
      parameter (inp=71)
      dimension frac0(mion2),ioo(mion2),idat(mion1)
      dimension gg(mion1,mdat),g0(mion1),z0(mion2)
      dimension uu(mion1,mdat),u0(mion1)
      dimension indxat(mion1,mdat),indxa(mion1)
      dimension u6(6),u7(7),u8(8),u10(10),u11(11)
      dimension u12(12),u13(13),u14(14),u16(16),u18(18),u20(20)
      dimension u24(24),u25(25),u26(26),u28(28)
      common/pfoptb/pfop(mtemp,melec,mstag),pfophm(mtemp,melec),
     *             frac(mtemp,melec,mstag),
     *             frop(mtemp,melec,mstag),itemp(mtemp),
     *             ntt
      equivalence (u6(1),uu(1,3)),(u7(1),uu(1,4)),(u8(1),uu(1,5))
      equivalence (u10(1),uu(1,6)),(u11(1),uu(1,7)),(u12(1),uu(1,8))
      equivalence (u13(1),uu(1,9)),(u14(1),uu(1,10)),(u16(1),uu(1,11))
      equivalence (u18(1),uu(1,12)),(u20(1),uu(1,13)),(u24(1),uu(1,14))
      equivalence (u25(1),uu(1,15)),(u26(1),uu(1,16)),(u28(1),uu(1,17))
      data idat   / 1, 2, 0, 0, 0, 3, 4, 5, 0, 6,
     *              7, 8, 9,10, 0,11, 0,12, 0,13,
     *              0, 0, 0,14,15,16, 0,17, 0, 0/ 
      data indxat /1,2,28*0,
     *            3,4,5,27*0,
     *            6,7,8,9,10,11,12,23*0,
     *            13,14,15,16,17,18,19,20,22*0,
     *            21,22,23,24,25,26,27,28,29,21*0,
     *            30,31,32,33,34,35,36,37,38,39,40,19*0,
     *            41,42,43,44,45,46,47,48,49,50,51,52,18*0,
     *            53,54,55,56,57,58,59,60,61,62,63,64,65,17*0,
     *            66,67,68,69,70,71,72,73,74,75,76,77,78,79,16*0,
     *            80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,15*0,
     *            95,96,97,98,99,100,101,102,103,104,105,106,107,
     *            108,109,110,111,13*0,
     *            112,113,114,115,116,117,118,119,120,121,122,123,
     *            124,125,126,127,128,129,130,11*0,
     *            131,132,133,134,135,136,137,138,139,140,141,142,
     *            143,144,145,146,147,148,149,150,151,9*0,
     *            152,153,154,155,156,157,158,159,160,161,162,163,
     *            164,165,166,167,168,169,170,171,172,173,174,175,
     *            176,5*0,
     *            177,178,179,180,181,182,183,184,185,186,187,188,
     *            189,190,191,192,193,194,195,196,197,198,199,200,
     *            201,202,4*0,
     *            203,204,205,206,207,208,209,210,211,212,213,214,
     *            215,216,217,218,219,220,221,222,223,224,225,226,
     *            227,228,229,3*0,
     *            230,231,232,233,234,235,236,237,238,239,240,241,
     *            242,243,244,245,246,247,248,249,250,251,252,253,
     *            254,255,256,257,258,0/
      data gg/2.,29*0.,
     *        2.,1.,28*0.,
     *        2.,1.,2.,1.,6.,9.,24*0.,
     *        2.,1.,2.,1.,6.,9.,4.,23*0.,
     *        2.,1.,2.,1.,6.,9.,4.,9.,22*0.,
     *        2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,20*0.,
     *        2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,19*0.,
     *        2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,1.,18*0.,
     *        2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,1.,6.,17*0.,
     *        2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,1.,6.,9.,16*0.,
     *        2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,1.,6.,9.,4.,9.,14*0.,
     *        2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,1.,6.,9.,4.,9.,6.,1.,
     *        12*0.,2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,1.,6.,9.,4.,9.,
     *        6.,1.,2.,1.,10*0.,2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,1.,
     *        6.,9.,4.,9.,6.,1.,10.,21.,28.,25.,6.,7.,6*0.,
     *        2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,1.,6.,9.,4.,9.,
     *           6.,1.,10.,21.,28.,25.,6.,7.,6.,5*0.,
     *        2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,1.,6.,9.,4.,9.,
     *           6.,1.,10.,21.,28.,25.,6.,25.,30.,25.,4*0.,
     *        2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,1.,6.,9.,4.,9.,
     *           6.,1.,10.,21.,28.,25.,6.,25.,28.,21.,10.,21.,0.,0./
      data uu(1,1)/109.6787/
      data uu(1,2)/198.3108/
      data uu(2,2)/438.9089/
      data u6/90.82,196.665,386.241,520.178,3162.395,3952.061/
      data u7/117.225,238.751,382.704,624.866,789.537,4452.758,5380.089/
      data u8/109.837,283.24,443.086,624.384,918.657,1114.008,5963.135,
     *        7028.393/
      data u10/173.93,330.391,511.8,783.3,1018.,1273.8,1671.792,
     *         1928.462,9645.005,10986.876/
      data u11/41.449,381.395,577.8,797.8,1116.2,1388.5,1681.5,2130.8,
     *         2418.7,11817.061,13297.676/
      data u12/61.671,121.268,646.41,881.1,1139.4,1504.3,1814.3,2144.7,
     *         2645.2,2964.4,14210.261,15829.951/
      data u13/48.278,151.86,229.446,967.8,1239.8,1536.3,1947.3,2295.4,
     *         2663.4,3214.8,3565.6,16825.022,18584.138/
      data u14/65.748,131.838,270.139,364.093,1345.1,1653.9,1988.4,
     *         2445.3,2831.9,3237.8,3839.8,4222.4,19661.693,21560.63/
      data u16/83.558,188.2,280.9,381.541,586.2,710.184,2265.9,2647.4,
     *         3057.7,3606.1,4071.4,4554.3,5255.9,5703.6,26002.663,
     *         28182.535/
      data u18/127.11,222.848,328.6,482.4,605.1,734.04,1002.73,1157.08,
     *         3407.3,3860.9,4347.,4986.6,5533.8,6095.5,6894.2,7404.4,
     *         33237.173,35699.936/
      data u20/49.306,95.752,410.642,542.6,681.6,877.4,1026.,1187.6,
     *         1520.64,1704.047,4774.,5301.,5861.,6595.,7215.,7860.,
     *         8770.,9338.,41366.,44177.41/
      data u24/54.576,132.966,249.7,396.5,560.2,731.02,1291.9,1490.,
     *         1688.,1971.,2184.,2404.,2862.,3098.52,8151.,8850.,
     *         9560.,10480.,11260.,12070.,13180.,13882.,60344.,63675.9/
      data u25/59.959,126.145,271.55,413.,584.,771.1,961.44,1569.,
     *         1789.,2003.,2307.,2536.,2771.,3250.,3509.82,9152.,
     *         9872.,10620.,11590.,12410.,13260.,14420.,15162.,
     *         65660.,69137.4/
      data u26/63.737,130.563,247.22,442.,605.,799.,1008.,1218.38,
     *         1884.,2114.,2341.,2668.,2912.,3163.,3686.,3946.82,
     *         10180.,10985.,11850.,12708.,13620.,14510.,15797.,
     *         16500.,71203.,74829.6/
      data u28/61.6,146.542,283.8,443.,613.5,870.,1070.,1310.,1560.,
     *         1812.,2589.,2840.,3100.,3470.,3740.,4020.,4606.,
     *         4896.2,12430.,13290.,14160.,15280.,16220.,17190.,
     *         18510.,19351.,82984.,86909.4/
c
c     If the routine is called with IAT=0, initialization:
c     read the ionization fractions form the Opacity Project tables,
c     and evaluate partition functions (assuming that the partition
c     function for the highest ion is eqaul to the statistical weight of
c     the ground state.
c     The table of partition functions, in the OP temperatures and
c     electron densities, are stored in the array PFOP
c
C     Read is done in a loop over the OP species
c
      if(iat.gt.0) go to 50
      do 40 iatnum=1,28
      if(idat(iatnum).eq.0) go to 40
c
      g0(iatnum+1)=1.
      indxa(iatnum+1)=indxat(iatnum+1,idat(iatnum))
      do i=1,iatnum
        ig0=iatnum-i+1
        g0(ig0)=gg(i,idat(iatnum))
        indxa(i)=indxat(i,idat(iatnum))
        u0(i)=uu(i,idat(iatnum))*1000.
      enddo
c
c     initializion of partition functions by the statistical weights of
c     the ground state
c
      do i=1,iatnum+1
         indx=indxa(i)
         pf0=g0(i)
         do it=1,mtemp
            do ieind=1,melec
               pfop(it,ieind,indx)=pf0
            enddo
         enddo
      enddo
c
      if(iatnum.eq.1) open(inp,file='atomic/ioniz.dat',status='old')
      read(inp,*)
      read(inp,*) it0,it1,itstp
      ntt=(it1-it0)/itstp+1
c
      do 30 it=1,ntt
         read(inp,*) itt,ie0,ie1,iestp
         itemp(it)=itt
         net=(ie1-ie0)/iestp+1
         t=exp(2.3025851*0.025*itt)
         safac0=sqrt(t)*t/2.07d-16
         tkcm=0.69496*t
         do 20 ie=1,net
            read(inp,601) iee,ion0,ion1,
     *                    (ioo(i),frac0(i),i=ion0+2,min(ion1,ion0+3)+2)
            ane=exp(2.3025851*0.25*iee)
            safac=safac0/ane
            nio=ion1-ion0
            if(nio.ge.3) then
               nlin=nio/4
               do ilin=1,nlin
                  read(inp,602) (ioo(i),frac0(i),
     *                 i=ion0+4*ilin+2,min(ion1,ion0+4*ilin+3)+2)
               end do
            end if
            ieind=iee/2
            ion0=ion0+2
            ion1=ion1+2
            do 10 ionn=ion0,ion1
              iation=iatnum+2-ionn
              if(ionn.lt.iatnum+2) then
                 if(ionn.eq.ion0) then
                    z0(ionn)=g0(iation)
                  else
                    z0(ionn)=frac0(ionn)/frac0(ionn-1)*safac*z0(ionn-1)
                    z0(ionn)=z0(ionn)*exp(-u0(iation)/tkcm)
                 endif
                 pfop(it,ieind,indxa(iation))=z0(ionn)
                 frac(it,ieind,indxa(iation))=frac0(ionn)
              else
                 u0hm=6090.5
                 z0hm=frac0(ionn)/frac0(ionn-1)*safac
                 z0hm=z0hm*exp(-u0hm/tkcm)
                 pfophm(it,ieind)=z0hm
              end if
   10       continue
   20    continue
   30 continue
   40 continue
  601 format(3i4,2x,4(i4,1x,e9.3))
  602 format(14x,4(i4,1x,e9.3))
      return
C
C     ----------------------------------------------------------
C
C     If the routine is called with IAT>0, evaluate the partition
C     function of atom IAT, ion ION,
C     for temperature T and electron density ANE
C     Evaluation is done by interpolation from previously computed
C     Opacity Project values
C
   50 continue
c
      xt=log10(t)
      kt0=2*int(20.*xt)
      xne=log10(ane)
      kn0=int(2.*xne)
c
      iatnum=iat
      if(idat(iatnum).eq.0) then
         if(ierr.lt.10) write(6,600) iatnum
         ierr=ierr+1
         iatnum=-1
  600    format(' Opfrac Warning: no data for element no. ',i3)
         opfrac=-1
         return
      end if
      indx=indxat(ion,idat(iatnum))
      if(kt0.lt.itemp(1)) then
         kt1=1
         write(6,611) t
  611    format(' (FRACOP) Extrapol. in T (low)',f7.0)
         go to 120
      endif
      if(kt0.ge.itemp(ntt)) then
         kt1=ntt-1
         write(6,612) t
  612    format(' (FRACOP) Extrapol. in T (high)',f12.0)
         go to 120
      endif
      do 110 it=1,ntt
         if(kt0.eq.itemp(it)) then
            kt1=it
            go to 120
         endif
  110 continue
  120 continue
      if(kn0.lt.1) then
         kn1=1
         go to 130
      endif
      if(kn0.ge.60) then
         kn1=59
         write(6,614) xne
  614    format(' (FRACOP) Extrapol. in Ne (high)',f9.4)
         go to 130
      endif
      kn1=kn0
  130 continue
      xt1=0.025*itemp(kt1)
      dxt=0.05
      at1=(xt-xt1)/dxt
      xn1=0.5*kn1
      dxn=0.5
      an1=(xne-xn1)/dxn
      x11=pfop(kt1,kn1,indx)
      x21=pfop(kt1+1,kn1,indx)
      x12=pfop(kt1,kn1+1,indx)
      x22=pfop(kt1+1,kn1+1,indx)
      x1221=x11*x21*x12*x22
      if(x1221.eq.0.) then
         xx1=x11+at1*(x21-x11)
         xx2=x12+at1*(x22-x12)
         rrx=xx1+an1*(xx2-xx1)
       else
         x11=log10(x11)
         x21=log10(x21)
         x12=log10(x12)
         x22=log10(x22)
         xx1=x11+at1*(x21-x11)
         xx2=x12+at1*(x22-x12)
         rrx=xx1+an1*(xx2-xx1)
         rrx=exp(2.3025851*rrx)
      endif
      pf=rrx
c
c     the next section - the true ionization fractions 
c
      x11=frac(kt1,kn1,indx)
      x21=frac(kt1+1,kn1,indx)
      x12=frac(kt1,kn1+1,indx)
      x22=frac(kt1+1,kn1+1,indx)
      x1221=x11*x21*x12*x22
      if(x1221.eq.0.) then
         xx1=x11+at1*(x21-x11)
         xx2=x12+at1*(x22-x12)
         rrx=xx1+an1*(xx2-xx1)
       else
         x11=log10(x11)
         x21=log10(x21)
         x12=log10(x12)
         x22=log10(x22)
         xx1=x11+at1*(x21-x11)
         xx2=x12+at1*(x22-x12)
         rrx=xx1+an1*(xx2-xx1)
         rrx=exp(2.3025851*rrx)
      endif
      fra=rrx
c
      opfrac=0
      return
      end
