      subroutine EOS(temp0, den0, pres0, eg0, cs0, f)
      implicit none
      save
      include 'vector_eos.dek'
      double precision abar,zbar
      common /xxx/ abar, zbar
c
      double precision temp0, den0, pres0, eg0, cs0
      integer f
c
      double precision delta, diff
      integer i
c
      jlo_eos = 1
      jhi_eos = 1
      abar_row(1) = abar
      zbar_row(1) = zbar
c
c     initial test
      temp_row(1) = 1.0e4
      den_row(1)  = den0
      call helmeos
      diff = eg0 - etot_row(1)
c
      if ((diff < 0.0).or.(diff <= 0.1*eg0)) then
         temp0 = temp_row(1)
         pres0 = ptot_row(1)
         cs0   = cs_row(1)
         f = 0
         return
      endif
c
c      if (temp0 < 1.0e4) then
c         temp0 = 1.0e4
c      endif
c
      temp_row(1) = 5.0e7
      den_row(1)  = den0
 999  continue
      call helmeos
      diff = (eg0-etot_row(1))/eg0
c
      i = 0
      do while (dabs(diff) > 1.0e-5)
         delta = -(etot_row(1) - eg0)/det_row(1)
         temp_row(1) = temp_row(1) + delta
         if (temp_row(1) < 1.0e4) then 
            f = 1
            return
         endif
         call helmeos
         diff = (eg0-etot_row(1))/eg0
         i = i + 1 
         if (eosfail) then
c            write(*,*) "FATAL!!!!!!!! EOS"
            f = 1
            return
         endif
         if (i > 50) then
            f = 1
            return
         endif
      enddo
c
      temp0 = temp_row(1)
      pres0 = ptot_row(1)
      cs0   = cs_row(1)
      f = 0
c
      return
      end
c
      subroutine EOS_RETURN(temp0, den0, pres0, cs0, eg0)
      implicit none
      save
      include 'vector_eos.dek'
      double precision abar,zbar
      common /xxx/ abar, zbar
c
      double precision temp0, den0, pres0, cs0, eg0
c
      jlo_eos = 1
      jhi_eos = 1
c
      if (temp0 <= 0.0) then
         temp0 = 1.0e4
      endif
c
      abar_row(1) = abar
      zbar_row(1) = zbar
      temp_row(1) = temp0
      den_row(1)  = den0
c
      call helmeos
c
      pres0 = ptot_row(1)
      cs0   = cs_row(1)
      eg0   = etot_row(1)
c
      return
      end
c
      subroutine EOS_RECOVER(temp0, den0, pres0, eg0, cs0, f)
      implicit none
      save
      include 'vector_eos.dek'
      double precision abar,zbar
      common /xxx/ abar, zbar
c
      double precision temp0, den0, pres0, eg0, cs0
      integer f
c
      double precision delta, diff
      integer i
c
c
c..here the pipeline is only 1 element long
      jlo_eos = 1
      jhi_eos = 1
c
      abar_row(1) = abar
      zbar_row(1) = zbar
      temp_row(1) = 1.0e4
      den_row(1)  = den0
      call helmeos

      temp0 = temp_row(1)
      pres0 = ptot_row(1)
      cs0   = cs_row(1)
      if (eg0-etot_row(1) < 0.0) then
         eg0 = etot_row(1)
      endif
c
      return
c
      diff = (eg0-etot_row(1))/eg0
      if ((diff < 0.0).or.(dabs(diff) < 1.0e-3)) then
         temp0 = temp_row(1)
         pres0 = ptot_row(1)
         cs0   = cs_row(1)
      else
         write(*,*) 'test ', diff, eg0, etot_row(1), 0.1*eg0
         stop
      endif
c
      return
      end
c
c
      subroutine EOS9(temp0, den0, pres0, eg0, cs0, f)
      implicit none
      save
      include 'vector_eos.dek'
      double precision abar,zbar
      common /xxx/ abar, zbar
c
      double precision temp0, den0, pres0, eg0, cs0
      double precision temp_new, tol, t0
      integer f
c
      double precision delta, diff
      integer i
c
      tol = 1.0e-3
      t0  = 1.001e3
c
      jlo_eos = 1
      jhi_eos = 1
      abar_row(1) = abar
      zbar_row(1) = zbar
      temp_row(1) = temp0
      den_row(1)  = den0
c
      call helmeos
c
      temp_new = temp_row(1) - (etot_row(1) - eg0) / det_row(1)
      temp_new = min(10.0*temp_row(1), 
     &     max(0.1*temp_row(1), temp_new))
      diff = dabs((temp_new-temp_row(1))/temp_row(1))
c
      temp_row(1) = temp_new
      if (temp_new < t0) then
         temp_row(1) = t0
         diff = 0.1*tol
      endif
c
      i = 2
      do while ( (i <= 50) .and. (diff > tol) )
         call helmeos
         temp_new = temp_row(1) - (etot_row(1) - eg0) / det_row(1)
         temp_new = min(10.0*temp_row(1), 
     &        max(0.1*temp_row(1), temp_new))
         diff = dabs((temp_new-temp_row(1))/temp_row(1))

         temp_row(1) = temp_new
         if (temp_new < t0) then
            temp_row(1) = t0
            diff = 0.1*tol
         endif

         i = i + 1
      enddo

      if (i > 49) then
         write(*,*) "FATAL!!!!!!!! EOS"
         stop
      endif
c
      call helmeos
c
      temp0 = temp_row(1)
      pres0 = ptot_row(1)
      cs0   = cs_row(1)
      f = 0
c
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine EOSX(temp0, den0, pres0, eg0, cs0, flagcoulomb)
      implicit none
      double precision abar,zbar
      common /xxx/ abar, zbar
c
      double precision temp0, den0, pres0, eg0, cs0
      double precision temp_new, tol, t0
      double precision temp1, eg1, pres1, dedt, cs1
      logical          eosfail
      double precision flagcoulomb
c
      double precision delta, diff
      integer i
c
      tol = 1.0e-3
      t0  = 1.001e3
c
      temp1 = temp0
      call helmeos2(temp1, den0, abar, zbar, pres1, eg1,
     &     dedt, cs1, eosfail, flagcoulomb)
c
      temp_new = temp1 - (eg1 - eg0) / dedt
      temp_new = min(10.0*temp1, 
     &     max(0.1*temp1, temp_new))
      diff = dabs((temp_new-temp1)/temp1)
c
      temp1 = temp_new
c      if (temp_new < t0) then
c         temp1 = t0
c         diff = 0.1*tol
c      endif
c
      i = 2
      do while ( (i <= 50) .and. (diff > tol) )
         call helmeos2(temp1, den0, abar, zbar, pres1, eg1,
     &        dedt, cs1, eosfail, flagcoulomb)
         temp_new = temp1 - (eg1 - eg0) / dedt
         temp_new = min(10.0*temp1, 
     &        max(0.1*temp1, temp_new))
         diff = dabs((temp_new-temp1)/temp1)
c
         temp1 = temp_new
         if (temp_new < t0) then
            temp1 = t0
            diff = 0.1*tol
         endif

         i = i + 1
      enddo

      if (i > 49) then
         write(*,*) "FATAL!!!!!!!! EOS A"
         write(*,*) temp0, den0, eg0
         stop
      endif
c
      call helmeos2(temp1, den0, abar, zbar, pres1, eg1,
     &     dedt, cs1, eosfail, flagcoulomb)
c
      temp0 = temp1
      pres0 = pres1
      cs0   = cs1
c
      return
      end
c
      subroutine EOSX_return(temp0, den0, pres0, eg0, cs0,
     &      flagcoulomb)
      implicit none
      double precision abar,zbar
      common /xxx/ abar, zbar
c
      double precision temp0, den0, pres0, eg0, cs0
      double precision temp1, eg1, pres1, dedt, cs1
      logical          eosfail
      double precision flagcoulomb
c
      call helmeos2(temp0, den0, abar, zbar, pres1, eg1,
     &     dedt, cs1, eosfail, flagcoulomb)
c
      pres0 = pres1
      cs0   = cs1
      eg0   = eg1
c
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine helm_init
      implicit none
      save
      include 'vector_eos.dek'
c
c..
c..tests the helmholtz eos routine
c..
c..ionmax  = number of isotopes in the network
c..xmass   = mass fractions
c..ymass   = molar fractions
c..aion    = number of nucleons
c..zion    = number of protons

      integer          ionmax
c      parameter        (ionmax=3)
c      parameter        (ionmax=1)
      parameter        (ionmax=2)
      double precision xmass(ionmax),ymass(ionmax),
     1                 aion(ionmax),zion(ionmax),temp,den,abar,zbar
      common /xxx/ abar, zbar
c
C c..set the mass fractions, z's and a's of the composition
C c..hydrogen
C       xmass(1) = 0.75d0
C       aion(1)  = 1.0d0
C       zion(1)  = 1.0d0

C c..helium
C       xmass(2) = 0.23d0
C       aion(2)  = 4.0d0
C       zion(2)  = 2.0d0

C c..carbon 12
C       xmass(3) = 0.02d0
C       aion(3)  = 12.0d0
C       zion(3)  = 6.0d0

ccc pure Fe56 version
c      xmass(1) = 1.0d0
c      aion(1)  = 56.0d0
c      zion(1)  = 28.0d0


CCCC CO White Dwarf
C c..carbon 12
       xmass(1) = 0.5d0
       aion(1)  = 12.0d0
       zion(1)  = 6.0d0
C c..oxygen 16
       xmass(2) = 0.5d0
       aion(2)  = 16.0d0
       zion(2)  = 8.0d0

c..get abar, zbar and a few other composition variables
      call azbar(xmass,aion,zion,ionmax,
     1           ymass,abar,zbar)
c
      return
      end
ccccc
      subroutine azbar(xmass,aion,zion,ionmax,
     1                 ymass,abar,zbar)
      implicit none
      save

c..this routine calculates composition variables for an eos routine

c..input:
c..mass fractions     = xmass(1:ionmax)
c..number of nucleons = aion(1:ionmax)
c..charge of nucleus  = zion(1:ionmax)
c..number of isotopes = ionmax

c..output:
c..molar abundances        = ymass(1:ionmax), 
c..mean number of nucleons = abar
c..mean nucleon charge     = zbar


c..declare
      integer          i,ionmax
      double precision xmass(ionmax),aion(ionmax),zion(ionmax),
     1                 ymass(ionmax),abar,zbar,zbarxx,ytot1

      zbarxx  = 0.0d0
      ytot1   = 0.0d0
      do i=1,ionmax
       ymass(i) = xmass(i)/aion(i)
       ytot1    = ytot1 + ymass(i)
       zbarxx   = zbarxx + zion(i) * ymass(i)
      enddo
      abar   = 1.0d0/ytot1
      zbar   = zbarxx * abar
      return
      end


      subroutine pretty_eos_out(whose)
      implicit none
      save
      include 'vector_eos.dek'
c..
c..writes a pretty output for the eos tester
c..
c..declare
      integer     j
      character*7 whose

c..popular formats
01    format(1x,t2,a,t11,'total',t24,'ion',t34,'e- & e+',
     1       t46,'radiation',t58,'coulomb')
02    format(1x,t2,a,1p6e12.4)
03    format(1x,t2,a6,1pe12.4,t22,a6,1pe12.4,
     1         t42,a6,1pe12.4,t62,a6,1pe12.4)



      do j=jlo_eos,jhi_eos


c..the input 
      write(6,03) 'temp =',temp_row(1),'den  =',den_row(1),
     1            'abar =',abar_row(1),'zbar =',zbar_row(1)
      write(6,*) ' ' 


c..and the output
c..first the totals from each of the components
      write(6,01)  whose
      write(6,02) 'pres =',
     1            ptot_row(j),pion_row(j),pele_row(j),
     2            prad_row(j),pcou_row(j)
      write(6,02) 'ener =',
     1            etot_row(j),eion_row(j),eele_row(j),
     2            erad_row(j),ecou_row(j)
      write(6,02) 'entr =',
     1            stot_row(j),sion_row(j),sele_row(j),
     2            srad_row(j),scou_row(j)

c..derivatives of the totals with respect to the input variables
      write(6,*)  ' '
      write(6,03) 'dp/dd=',dpd_row(j),'dp/dt=',dpt_row(j),
     1            'dp/da=',dpa_row(j),'dp/dz=',dpz_row(j)
      write(6,03) 'de/dd=',ded_row(j),'de/dt=',det_row(j),
     1            'de/da=',dea_row(j),'de/dz=',dez_row(j)
      write(6,03) 'ds/dd=',dsd_row(j),'ds/dt=',dst_row(j),
     1            'ds/da=',dsa_row(j),'ds/dz=',dsz_row(j)


c..derivatives of the electron-positron compoenets with
c..respect to the input variables
      write(6,*) ' ' 
      write(6,03) 'dpepd=',dpepd_row(j),'dpept=',dpept_row(j),
     1            'dpepa=',dpepa_row(j),'dpepz=',dpepz_row(j)
      write(6,03) 'deepd=',deepd_row(j),'deept=',deept_row(j),
     1            'deepa=',deepa_row(j),'deepz=',deepz_row(j)
      write(6,03) 'dsepd=',dsepd_row(j),'dsept=',dsept_row(j),
     1            'dsepa=',dsepa_row(j),'dsepz=',dsepz_row(j)


c..the thermodynamic consistency relations, these should all be
c..at the floating poiint limit of zero
      write(6,*) ' ' 
      write(6,03) 'maxw1=',dse_row(j),'maxw2=',dpe_row(j),
     1            'maxw3=',dsp_row(j)


c..number density of electrons, poistrons, matter electrons, and ions
      write(6,03) 'xne  =',xne_row(j),'xnp  =',xnp_row(j),
     1            'xnem =',xnem_row(j),'xni  =',xni_row(j)


c..derivatibves of the electron number density with 
c..respect to the input variables
      write(6,03) 'dxned=',dxned_row(j),'dxnet=',dxnet_row(j),
     1            'dxnea=',dxnea_row(j),'dxnez=',dxnez_row(j)


c..electron chemical potential, positron chemical potential
c..and derivatives of electron chemical potential with respect
c..to the input variables
      write(6,03) 'eta  =',etaele_row(j),'etap =',etapos_row(j)
      write(6,03) 'detad=',detad_row(j),'detat=',detat_row(j),
     1            'detaa=',detaa_row(j),'detaz=',detaz_row(j)


c..specific heats, and ratio of electostatic to thermal energy
      write(6,03) 'cp   =',cp_row(j),'cv   =',cv_row(j),
     1            'plasg=',plasg_row(j)

c..the 3 gammas and the sound speed
      write(6,03) 'gam1 =',gam1_row(j),'gam2 =',gam2_row(j),
     1            'gam3 =',gam3_row(j),'csond=',cs_row(j)
      write(6,*) ' '

      enddo
      return
      end





      subroutine helmeos
      implicit none
      save
      include 'vector_eos.dek'


c..given a temperature temp [K], density den [g/cm**3], and a composition 
c..characterized by abar and zbar, this routine returns most of the other 
c..thermodynamic quantities. of prime interest is the pressure [erg/cm**3], 
c..specific thermal energy [erg/gr], the entropy [erg/g/K], along with 
c..their derivatives with respect to temperature, density, abar, and zbar.
c..other quantites such the normalized chemical potential eta (plus its
c..derivatives), number density of electrons and positron pair (along 
c..with their derivatives), adiabatic indices, specific heats, and 
c..relativistically correct sound speed are also returned.
c..
c..this routine assumes planckian photons, an ideal gas of ions,
c..and an electron-positron gas with an arbitrary degree of relativity
c..and degeneracy. interpolation in a table of the helmholtz free energy
c..is used to return the electron-positron thermodynamic quantities.
c..all other derivatives are analytic.
c..
c..references: cox & giuli chapter 24 ; timmes & swesty apj 1999


c..declare
      double precision pi,amu,kerg,clight,avo,qe,h,ssol,asol
      parameter       (pi      = 3.1415926535897932384d0,
     1                  amu    = 1.6605402d-24,
     2                  kerg   = 1.380658d-16,
     3                  clight = 2.99792458d10, 
     4                  avo    = 6.0221367d23,
     5                  qe     = 4.8032068d-10,  
     6                  h      = 6.6260755d-27,
     7                  ssol   = 5.67051d-5,
     8                  asol   = 4.0d0 * ssol / clight)

      integer          i,j
      double precision x,y,zz,zzi,deni,tempi,xni,dxnidd,dxnida,
     1                 dpepdt,dpepdd,deepdt,deepdd,dsepdd,dsepdt,
     2                 dpraddd,dpraddt,deraddd,deraddt,dpiondd,dpiondt,
     3                 deiondd,deiondt,dsraddd,dsraddt,dsiondd,dsiondt,
     4                 dse,dpe,dsp,kt,ktinv,prad,erad,srad,pion,eion,
     5                 sion,xnem,pele,eele,sele,pres,ener,entr,dpresdd,
     6                 dpresdt,denerdd,denerdt,dentrdd,dentrdt,cv,cp,
     7                 gam1,gam2,gam3,chit,chid,nabad,sound,etaele,
     8                 detadt,detadd,xnefer,dxnedt,dxnedd,s,
     9                 temp,den,abar,zbar,ytot1,ye,
     &                 sioncon,forth,forpi,kergavo,ikavo,asoli3,light2

      parameter        (sioncon = (2.0d0 * pi * amu * kerg)/(h*h),
     1                  forth   = 4.0d0/3.0d0,
     2                  forpi   = 4.0d0 * pi,
     3                  kergavo = kerg * avo, 
     4                  ikavo   = 1.0d0/kergavo,
     5                  asoli3  = asol/3.0d0,
     6                  light2  = clight * clight)

c..for the abar derivatives
      double precision dpradda,deradda,dsradda,
     1                 dpionda,deionda,dsionda,
     2                 dpepda,deepda,dsepda,
     3                 dpresda,denerda,dentrda,
     4                 detada,dxneda


c..for the zbar derivatives
      double precision dpraddz,deraddz,dsraddz,
     1                 dpiondz,deiondz,dsiondz,
     2                 dpepdz,deepdz,dsepdz,
     3                 dpresdz,denerdz,dentrdz,
     4                 detadz,dxnedz


c..for the tables, in general
      integer          imax,jmax
c      parameter        (imax = 211, jmax = 71)
      parameter        (imax = 271, jmax = 101)
      double precision d(imax),t(jmax)

c..for the helmholtz free energy tables
      double precision f(imax,jmax),fd(imax,jmax),
     1                 ft(imax,jmax),fdd(imax,jmax),ftt(imax,jmax),
     2                 fdt(imax,jmax),fddt(imax,jmax),fdtt(imax,jmax),
     3                 fddtt(imax,jmax)

c..for the pressure derivative with density ables
      double precision dpdf(imax,jmax),dpdfd(imax,jmax),
     1                 dpdft(imax,jmax),dpdfdd(imax,jmax),
     2                 dpdftt(imax,jmax),dpdfdt(imax,jmax)

c..for chemical potential tables
      double precision ef(imax,jmax),efd(imax,jmax),
     1                 eft(imax,jmax),efdd(imax,jmax),eftt(imax,jmax),
     2                 efdt(imax,jmax)

c..for the number density tables
      double precision xf(imax,jmax),xfd(imax,jmax),
     1                 xft(imax,jmax),xfdd(imax,jmax),xftt(imax,jmax),
     2                 xfdt(imax,jmax)

c..for the interpolations
      integer          iat,jat
      double precision tlo,thi,tstp,tstpi,dlo,dhi,dstp,dstpi,
     1                 tsav,dsav,free,df_d,df_t,df_dd,df_tt,df_dt
      double precision dth,dt2,dti,dt2i,dd,dd2,ddi,dd2i,xt,xd,mxt,mxd,
     1                 si0t,si1t,si2t,si0mt,si1mt,si2mt,
     2                 si0d,si1d,si2d,si0md,si1md,si2md,
     3                 dsi0t,dsi1t,dsi2t,dsi0mt,dsi1mt,dsi2mt,
     4                 dsi0d,dsi1d,dsi2d,dsi0md,dsi1md,dsi2md,
     5                 ddsi0t,ddsi1t,ddsi2t,ddsi0mt,ddsi1mt,ddsi2mt,
     6                 ddsi0d,ddsi1d,ddsi2d,ddsi0md,ddsi1md,ddsi2md,
     7                 z,psi0,dpsi0,ddpsi0,psi1,dpsi1,ddpsi1,psi2,
     8                 dpsi2,ddpsi2,din,h5,fi(36),
     9                 xpsi0,xdpsi0,xpsi1,xdpsi1,h3,
     1                 w0t,w1t,w2t,w0mt,w1mt,w2mt,
     2                 w0d,w1d,w2d,w0md,w1md,w2md

c..for storing the differences
      double precision dt_sav(jmax),dt2_sav(jmax),
     1                 dti_sav(jmax),dt2i_sav(jmax),
     2                 dd_sav(imax),dd2_sav(imax),
     3                 ddi_sav(imax),dd2i_sav(imax)


c..for the coulomb corrections
      double precision dsdd,dsda,lami,inv_lami,lamida,lamidd,
     1                 plasg,plasgdd,plasgdt,plasgda,plasgdz,
     1                 a1,b1,c1,d1,e1,a2,b2,c2,
     3                 ecoul,decouldd,decouldt,decoulda,decouldz,
     4                 pcoul,dpcouldd,dpcouldt,dpcoulda,dpcouldz,
     5                 scoul,dscouldd,dscouldt,dscoulda,dscouldz,
     6                 tmelt,tfermi,rhocond,z2,x1,x2,third,esqu
      parameter        (a1    = -0.898004d0, 
     1                  b1    =  0.96786d0, 
     2                  c1    =  0.220703d0, 
     3                  d1    = -0.86097d0,
     4                  e1    =  2.5269d0, 
     5                  a2    =  0.29561d0, 
     6                  b2    =  1.9885d0,    
     7                  c2    =  0.288675d0,
     8                  third = 1.0d0/3.0d0,
     9                  esqu  = qe * qe)


c..for initialization
      integer          ifirst
      data             ifirst/0/ 


c..quintic hermite polynomial statement functions
c..psi0 and its derivatives
      psi0(z)   = z**3 * ( z * (-6.0d0*z + 15.0d0) -10.0d0) + 1.0d0
      dpsi0(z)  = z**2 * ( z * (-30.0d0*z + 60.0d0) - 30.0d0)
      ddpsi0(z) = z* ( z*( -120.0d0*z + 180.0d0) -60.0d0)

c..psi1 and its derivatives
      psi1(z)   = z* ( z**2 * ( z * (-3.0d0*z + 8.0d0) - 6.0d0) + 1.0d0)
      dpsi1(z)  = z*z * ( z * (-15.0d0*z + 32.0d0) - 18.0d0) +1.0d0
      ddpsi1(z) = z * (z * (-60.0d0*z + 96.0d0) -36.0d0)

c..psi2  and its derivatives
      psi2(z)   = 0.5d0*z*z*( z* ( z * (-z + 3.0d0) - 3.0d0) + 1.0d0)
      dpsi2(z)  = 0.5d0*z*( z*(z*(-5.0d0*z + 12.0d0) - 9.0d0) + 2.0d0)
      ddpsi2(z) = 0.5d0*(z*( z * (-20.0d0*z + 36.0d0) - 18.0d0) + 2.0d0)

c..biquintic hermite polynomial statement function
      h5(i,j,w0t,w1t,w2t,w0mt,w1mt,w2mt,w0d,w1d,w2d,w0md,w1md,w2md)=
     1       fi(1)  *w0d*w0t   + fi(2)  *w0md*w0t
     2     + fi(3)  *w0d*w0mt  + fi(4)  *w0md*w0mt
     3     + fi(5)  *w0d*w1t   + fi(6)  *w0md*w1t
     4     + fi(7)  *w0d*w1mt  + fi(8)  *w0md*w1mt
     5     + fi(9)  *w0d*w2t   + fi(10) *w0md*w2t
     6     + fi(11) *w0d*w2mt  + fi(12) *w0md*w2mt
     7     + fi(13) *w1d*w0t   + fi(14) *w1md*w0t
     8     + fi(15) *w1d*w0mt  + fi(16) *w1md*w0mt
     9     + fi(17) *w2d*w0t   + fi(18) *w2md*w0t
     &     + fi(19) *w2d*w0mt  + fi(20) *w2md*w0mt
     1     + fi(21) *w1d*w1t   + fi(22) *w1md*w1t
     2     + fi(23) *w1d*w1mt  + fi(24) *w1md*w1mt
     3     + fi(25) *w2d*w1t   + fi(26) *w2md*w1t
     4     + fi(27) *w2d*w1mt  + fi(28) *w2md*w1mt
     5     + fi(29) *w1d*w2t   + fi(30) *w1md*w2t
     6     + fi(31) *w1d*w2mt  + fi(32) *w1md*w2mt
     7     + fi(33) *w2d*w2t   + fi(34) *w2md*w2t
     8     + fi(35) *w2d*w2mt  + fi(36) *w2md*w2mt



c..cubic hermite polynomial statement functions
c..psi0 & derivatives
      xpsi0(z)  = z * z * (2.0d0*z - 3.0d0) + 1.0
      xdpsi0(z) = z * (6.0d0*z - 6.0d0)

c..psi1 & derivatives
      xpsi1(z)  = z * ( z * (z - 2.0d0) + 1.0d0)
      xdpsi1(z) = z * (3.0d0*z - 4.0d0) + 1.0d0


c..bicubic hermite polynomial statement function
      h3(i,j,w0t,w1t,w0mt,w1mt,w0d,w1d,w0md,w1md) = 
     1       fi(1)  *w0d*w0t   +  fi(2)  *w0md*w0t 
     2     + fi(3)  *w0d*w0mt  +  fi(4)  *w0md*w0mt
     3     + fi(5)  *w0d*w1t   +  fi(6)  *w0md*w1t 
     4     + fi(7)  *w0d*w1mt  +  fi(8)  *w0md*w1mt
     5     + fi(9)  *w1d*w0t   +  fi(10) *w1md*w0t 
     6     + fi(11) *w1d*w0mt  +  fi(12) *w1md*w0mt
     7     + fi(13) *w1d*w1t   +  fi(14) *w1md*w1t 
     8     + fi(15) *w1d*w1mt  +  fi(16) *w1md*w1mt



c..popular format statements
01    format(1x,5(a,1pe11.3))
02    format(1x,a,1p4e16.8)
03    format(1x,4(a,1pe11.3))
04    format(1x,4(a,i4))


c..do this stuff once
      if (ifirst .eq. 0) then
       ifirst = 1

c..open the table
c      open(unit=2,file='helm_table_large.dat',status='old')
c      open(unit=2,file='heos/helm_table_large.dat',status='old')
      open(unit=2,file='data.code/helm_table_large.dat',status='old')

c..read the helmholtz free energy table
c       tlo   = 4.0d0
c       thi   = 11.0d0
c       tstp  = (thi - tlo)/float(jmax-1)
c       tstpi = 1.0d0/tstp
c       dlo   = -10.0d0
c       dhi   = 11.0d0
c
       tlo   = 3.0d0
       thi   = 13.0d0
       tstp  = (thi - tlo)/float(jmax-1)
       tstpi = 1.0d0/tstp
       dlo   = -12.0d0
       dhi   = 15.0d0
c
       dstp  = (dhi - dlo)/float(imax-1)
       dstpi = 1.0d0/dstp
       do j=1,jmax
        tsav = tlo + (j-1)*tstp
        t(j) = 10.0d0**(tsav)
        do i=1,imax
         dsav = dlo + (i-1)*dstp
         d(i) = 10.0d0**(dsav)
         read(2,*) f(i,j),fd(i,j),ft(i,j),fdd(i,j),ftt(i,j),fdt(i,j),
     1            fddt(i,j),fdtt(i,j),fddtt(i,j)
        enddo
       enddo

c..read the pressure derivative with density table
       do j=1,jmax
        do i=1,imax
         read(2,*) dpdf(i,j),dpdfd(i,j),dpdft(i,j),dpdfdt(i,j)
        enddo
       enddo

c..read the electron chemical potential table
       do j=1,jmax
        do i=1,imax
         read(2,*) ef(i,j),efd(i,j),eft(i,j),efdt(i,j)
        enddo
       enddo

c..read the number density table
       do j=1,jmax
        do i=1,imax
         read(2,*) xf(i,j),xfd(i,j),xft(i,j),xfdt(i,j)
        enddo
       enddo

c..construct the temperature and density deltas and their inverses 
       do j=1,jmax-1
        dth          = t(j+1) - t(j)
        dt2         = dth * dth
        dti         = 1.0d0/dth
        dt2i        = 1.0d0/dt2
        dt_sav(j)   = dth
        dt2_sav(j)  = dt2
        dti_sav(j)  = dti
        dt2i_sav(j) = dt2i
       end do
       do i=1,imax-1
        dd          = d(i+1) - d(i)
        dd2         = dd * dd
        ddi         = 1.0d0/dd
        dd2i        = 1.0d0/dd2
        dd_sav(i)   = dd
        dd2_sav(i)  = dd2
        ddi_sav(i)  = ddi
        dd2i_sav(i) = dd2i
       enddo

       close(unit=2)
       write(6,*)
       write(6,*) 'finished reading eos table'
       write(6,04) 'imax=',imax,' jmax=',jmax
       write(6,03) 'temp(1)   =',t(1),' temp(jmax)   =',t(jmax)
       write(6,03) 'ye*den(1) =',d(1),' ye*den(imax) =',d(imax)
       write(6,*)

      end if



c..start of vectorization loop, normal executaion starts here
      eosfail = .false.
      do j=jlo_eos,jhi_eos

       if (temp_row(j) .le. 0.0) stop 'temp less than 0 in helmeos'
       if (den_row(j)  .le. 0.0) stop 'den less than 0 in helmeos'

       temp  = temp_row(j)
       den   = den_row(j)
       abar  = abar_row(j)
       zbar  = zbar_row(j)
       ytot1 = 1.0d0/abar
       ye    = ytot1 * zbar

c..initialize
       deni    = 1.0d0/den
       tempi   = 1.0d0/temp 
       kt      = kerg * temp
       ktinv   = 1.0d0/kt


c..radiation section:
       prad    = asoli3 * temp * temp * temp * temp
       dpraddd = 0.0d0
       dpraddt = 4.0d0 * prad*tempi
       dpradda = 0.0d0
       dpraddz = 0.0d0

       erad    = 3.0d0 * prad*deni
       deraddd = -erad*deni
       deraddt = 3.0d0 * dpraddt*deni
       deradda = 0.0d0
       deraddz = 0.0d0


       srad    = (prad*deni + erad)*tempi
       dsraddd = (dpraddd*deni - prad*deni*deni + deraddd)*tempi
       dsraddt = (dpraddt*deni + deraddt - srad)*tempi
       dsradda = 0.0d0
       dsraddz = 0.0d0

c..ion section:
       xni     = avo * ytot1 * den
       dxnidd  = avo * ytot1
       dxnida  = -xni * ytot1

       pion    = xni * kt
       dpiondd = dxnidd * kt
       dpiondt = xni * kerg
       dpionda = dxnida * kt 
       dpiondz = 0.0d0

       eion    = 1.5d0 * pion*deni
       deiondd = (1.5d0 * dpiondd - eion)*deni
       deiondt = 1.5d0 * dpiondt*deni
       deionda = 1.5d0 * dpionda*deni
       deiondz = 0.0d0
    

       x       = abar*abar*sqrt(abar) * deni/avo
       s       = sioncon * temp
       z       = x * s * sqrt(s)
       y       = log(z)
       sion    = (pion*deni + eion)*tempi + kergavo * ytot1 * y
       dsiondd = (dpiondd*deni - pion*deni*deni + deiondd)*tempi
     1            - kergavo * deni * ytot1
       dsiondt = (dpiondt*deni + deiondt)*tempi - 
     1           (pion*deni + eion) * tempi*tempi 
     2           + 1.5d0 * kergavo * tempi*ytot1
       x       = avo*kerg/abar
       dsionda = (dpionda*deni + deionda)*tempi 
     1           + kergavo*ytot1*ytot1* (2.5d0 - y)
       dsiondz = 0.0d0


c..electron-positron section:
c..assume complete ionization 
       xnem    = xni * zbar

c..enter the table with ye*den
       din = ye*den

c..bomb proof the input
       if (temp .gt. t(jmax)) then
        write(6,01) 'temp=',temp,' t(jmax)=',t(jmax)
        write(6,*) 'temp too hot, off grid'       
        write(6,*) 'setting eosfail to true and returning'
        call flush(6)
        eosfail = .true.
        return
       end if
       if (temp .lt. t(1)) then
        write(6,01) 'temp=',temp,' t(1)=',t(1)
        write(6,*) 'temp too cold, off grid'
        write(6,*) 'setting eosfail to true and returning'
        call flush(6)
        eosfail = .true.
        return
       end if
       if (din  .gt. d(imax)) then
        write(6,01) 'den*ye=',din,' d(imax)=',d(imax)
        write(6,*) 'ye*den too big, off grid'
        write(6,*) 'setting eosfail to true and returning'
        call flush(6)
        eosfail = .true.
        return
       end if
       if (din  .lt. d(1)) then
        write(6,01) 'ye*den=',din,' d(1)=',d(1)
        write(6,*) 'ye*den too small, off grid'
        write(6,*) 'setting eosfail to true and returning'
        call flush(6)
        eosfail = .true.
        return
       end if

c..hash locate this temperature and density
       jat = int((log10(temp) - tlo)*tstpi) + 1
       jat = max(1,min(jat,jmax-1))
       iat = int((log10(din) - dlo)*dstpi) + 1
       iat = max(1,min(iat,imax-1))


c..access the table locations only once
       fi(1)  = f(iat,jat)
       fi(2)  = f(iat+1,jat)
       fi(3)  = f(iat,jat+1)
       fi(4)  = f(iat+1,jat+1)
       fi(5)  = ft(iat,jat)
       fi(6)  = ft(iat+1,jat)
       fi(7)  = ft(iat,jat+1)
       fi(8)  = ft(iat+1,jat+1)
       fi(9)  = ftt(iat,jat)
       fi(10) = ftt(iat+1,jat)
       fi(11) = ftt(iat,jat+1)
       fi(12) = ftt(iat+1,jat+1)
       fi(13) = fd(iat,jat)
       fi(14) = fd(iat+1,jat)
       fi(15) = fd(iat,jat+1)
       fi(16) = fd(iat+1,jat+1)
       fi(17) = fdd(iat,jat)
       fi(18) = fdd(iat+1,jat)
       fi(19) = fdd(iat,jat+1)
       fi(20) = fdd(iat+1,jat+1)
       fi(21) = fdt(iat,jat)
       fi(22) = fdt(iat+1,jat)
       fi(23) = fdt(iat,jat+1)
       fi(24) = fdt(iat+1,jat+1)
       fi(25) = fddt(iat,jat)
       fi(26) = fddt(iat+1,jat)
       fi(27) = fddt(iat,jat+1)
       fi(28) = fddt(iat+1,jat+1)
       fi(29) = fdtt(iat,jat)
       fi(30) = fdtt(iat+1,jat)
       fi(31) = fdtt(iat,jat+1)
       fi(32) = fdtt(iat+1,jat+1)
       fi(33) = fddtt(iat,jat)
       fi(34) = fddtt(iat+1,jat)
       fi(35) = fddtt(iat,jat+1)
       fi(36) = fddtt(iat+1,jat+1)


c..various differences
       xt  = max( (temp - t(jat))*dti_sav(jat), 0.0d0)
       xd  = max( (din - d(iat))*ddi_sav(iat), 0.0d0)
       mxt = 1.0d0 - xt
       mxd = 1.0d0 - xd

c..the six density and six temperature basis functions
       si0t =   psi0(xt)
       si1t =   psi1(xt)*dt_sav(jat)
       si2t =   psi2(xt)*dt2_sav(jat)

       si0mt =  psi0(mxt)
       si1mt = -psi1(mxt)*dt_sav(jat)
       si2mt =  psi2(mxt)*dt2_sav(jat)

       si0d =   psi0(xd)
       si1d =   psi1(xd)*dd_sav(iat)
       si2d =   psi2(xd)*dd2_sav(iat)

       si0md =  psi0(mxd)
       si1md = -psi1(mxd)*dd_sav(iat)
       si2md =  psi2(mxd)*dd2_sav(iat)

c..derivatives of the weight functions
       dsi0t =   dpsi0(xt)*dti_sav(jat)
       dsi1t =   dpsi1(xt)
       dsi2t =   dpsi2(xt)*dt_sav(jat)

       dsi0mt = -dpsi0(mxt)*dti_sav(jat)
       dsi1mt =  dpsi1(mxt)
       dsi2mt = -dpsi2(mxt)*dt_sav(jat)

       dsi0d =   dpsi0(xd)*ddi_sav(iat)
       dsi1d =   dpsi1(xd)
       dsi2d =   dpsi2(xd)*dd_sav(iat)

       dsi0md = -dpsi0(mxd)*ddi_sav(iat)
       dsi1md =  dpsi1(mxd)
       dsi2md = -dpsi2(mxd)*dd_sav(iat)

c..second derivatives of the weight functions
       ddsi0t =   ddpsi0(xt)*dt2i_sav(jat)
       ddsi1t =   ddpsi1(xt)*dti_sav(jat)
       ddsi2t =   ddpsi2(xt)

       ddsi0mt =  ddpsi0(mxt)*dt2i_sav(jat)
       ddsi1mt = -ddpsi1(mxt)*dti_sav(jat)
       ddsi2mt =  ddpsi2(mxt)

c       ddsi0d =   ddpsi0(xd)*dd2i_sav(iat)
c       ddsi1d =   ddpsi1(xd)*ddi_sav(iat)
c       ddsi2d =   ddpsi2(xd)

c       ddsi0md =  ddpsi0(mxd)*dd2i_sav(iat)
c       ddsi1md = -ddpsi1(mxd)*ddi_sav(iat)
c       ddsi2md =  ddpsi2(mxd)


c..the free energy
       free  = h5(iat,jat,
     1         si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt,
     2         si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

c..derivative with respect to density
       df_d  = h5(iat,jat,
     1         si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt,
     2         dsi0d,  dsi1d,  dsi2d,  dsi0md,  dsi1md,  dsi2md)

c..derivative with respect to temperature
       df_t = h5(iat,jat,
     1         dsi0t,  dsi1t,  dsi2t,  dsi0mt,  dsi1mt,  dsi2mt,
     2         si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

c..derivative with respect to density**2
c       df_dd = h5(iat,jat,
c     1         si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt,
c     2         ddsi0d, ddsi1d, ddsi2d, ddsi0md, ddsi1md, ddsi2md)

c..derivative with respect to temperature**2
       df_tt = h5(iat,jat,
     1       ddsi0t, ddsi1t, ddsi2t, ddsi0mt, ddsi1mt, ddsi2mt,
     2         si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

c..derivative with respect to temperature and density
       df_dt = h5(iat,jat,
     1         dsi0t,  dsi1t,  dsi2t,  dsi0mt,  dsi1mt,  dsi2mt,
     2         dsi0d,  dsi1d,  dsi2d,  dsi0md,  dsi1md,  dsi2md)



c..now get the pressure derivative with density, chemical potential, and 
c..electron positron number densities
c..get the interpolation weight functions
       si0t   =  xpsi0(xt)
       si1t   =  xpsi1(xt)*dt_sav(jat)

       si0mt  =  xpsi0(mxt)
       si1mt  =  -xpsi1(mxt)*dt_sav(jat)

       si0d   =  xpsi0(xd)
       si1d   =  xpsi1(xd)*dd_sav(iat)

       si0md  =  xpsi0(mxd)
       si1md  =  -xpsi1(mxd)*dd_sav(iat)


c..derivatives of weight functions
       dsi0t  = xdpsi0(xt)*dti_sav(jat)
       dsi1t  = xdpsi1(xt)

       dsi0mt = -xdpsi0(mxt)*dti_sav(jat)
       dsi1mt = xdpsi1(mxt)

       dsi0d  = xdpsi0(xd)*ddi_sav(iat)
       dsi1d  = xdpsi1(xd)

       dsi0md = -xdpsi0(mxd)*ddi_sav(iat)
       dsi1md = xdpsi1(mxd)


c..look in the pressure derivative only once
       fi(1)  = dpdf(iat,jat)
       fi(2)  = dpdf(iat+1,jat)
       fi(3)  = dpdf(iat,jat+1)
       fi(4)  = dpdf(iat+1,jat+1)
       fi(5)  = dpdft(iat,jat)
       fi(6)  = dpdft(iat+1,jat)
       fi(7)  = dpdft(iat,jat+1)
       fi(8)  = dpdft(iat+1,jat+1)
       fi(9)  = dpdfd(iat,jat)
       fi(10) = dpdfd(iat+1,jat)
       fi(11) = dpdfd(iat,jat+1)
       fi(12) = dpdfd(iat+1,jat+1)
       fi(13) = dpdfdt(iat,jat)
       fi(14) = dpdfdt(iat+1,jat)
       fi(15) = dpdfdt(iat,jat+1)
       fi(16) = dpdfdt(iat+1,jat+1)

c..pressure derivative with density
       dpepdd  = h3(iat,jat,
     1                 si0t,   si1t,   si0mt,   si1mt,
     2                 si0d,   si1d,   si0md,   si1md)
       dpepdd  = max(ye * dpepdd,0.0d0)



c..look in the electron chemical potential table only once
       fi(1)  = ef(iat,jat)
       fi(2)  = ef(iat+1,jat)
       fi(3)  = ef(iat,jat+1)
       fi(4)  = ef(iat+1,jat+1)
       fi(5)  = eft(iat,jat)
       fi(6)  = eft(iat+1,jat)
       fi(7)  = eft(iat,jat+1)
       fi(8)  = eft(iat+1,jat+1)
       fi(9)  = efd(iat,jat)
       fi(10) = efd(iat+1,jat)
       fi(11) = efd(iat,jat+1)
       fi(12) = efd(iat+1,jat+1)
       fi(13) = efdt(iat,jat)
       fi(14) = efdt(iat+1,jat)
       fi(15) = efdt(iat,jat+1)
       fi(16) = efdt(iat+1,jat+1)


c..electron chemical potential etaele
       etaele  = h3(iat,jat,
     1               si0t,   si1t,   si0mt,   si1mt,
     2               si0d,   si1d,   si0md,   si1md)

c..derivative with respect to density
       x       = h3(iat,jat,
     1               si0t,   si1t,   si0mt,   si1mt,
     2              dsi0d,  dsi1d,  dsi0md,  dsi1md)
       detadd  = ye * x

c..derivative with respect to temperature
       detadt  = h3(iat,jat,
     1              dsi0t,  dsi1t,  dsi0mt,  dsi1mt,
     2               si0d,   si1d,   si0md,   si1md)

c..derivative with respect to abar and zbar
      detada = -x * din * ytot1
      detadz =  x * den * ytot1



c..look in the number density table only once
       fi(1)  = xf(iat,jat)
       fi(2)  = xf(iat+1,jat)
       fi(3)  = xf(iat,jat+1)
       fi(4)  = xf(iat+1,jat+1)
       fi(5)  = xft(iat,jat)
       fi(6)  = xft(iat+1,jat)
       fi(7)  = xft(iat,jat+1)
       fi(8)  = xft(iat+1,jat+1)
       fi(9)  = xfd(iat,jat)
       fi(10) = xfd(iat+1,jat)
       fi(11) = xfd(iat,jat+1)
       fi(12) = xfd(iat+1,jat+1)
       fi(13) = xfdt(iat,jat)
       fi(14) = xfdt(iat+1,jat)
       fi(15) = xfdt(iat,jat+1)
       fi(16) = xfdt(iat+1,jat+1)

c..electron + positron number densities
      xnefer   = h3(iat,jat,
     1               si0t,   si1t,   si0mt,   si1mt,
     2               si0d,   si1d,   si0md,   si1md)

c..derivative with respect to density
      x        = h3(iat,jat,
     1               si0t,   si1t,   si0mt,   si1mt,
     2              dsi0d,  dsi1d,  dsi0md,  dsi1md)
      x = max(x,0.0d0)
      dxnedd   = ye * x

c..derivative with respect to temperature
      dxnedt   = h3(iat,jat,
     1              dsi0t,  dsi1t,  dsi0mt,  dsi1mt,
     2               si0d,   si1d,   si0md,   si1md)

c..derivative with respect to abar and zbar
      dxneda = -x * din * ytot1
      dxnedz =  x  * den * ytot1
       


c..the desired electron-positron thermodynamic quantities

c..dpepdd at high temperatures and low densities is below the
c..floating point limit of the subtraction of two large terms.
c..since dpresdd doesn't enter the maxwell relations at all, use the
c..bicubic interpolation done above instead of this one
       x       = din * din
       pele    = x * df_d
       dpepdt  = x * df_dt
c       dpepdd  = ye * (x * df_dd + 2.0d0 * din * df_d)
       s       = dpepdd/ye - 2.0d0 * din * df_d
       dpepda  = -ytot1 * (2.0d0 * pele + s * din)
       dpepdz  = den*ytot1*(2.0d0 * din * df_d  +  s)


       x       = ye * ye
       sele    = -df_t * ye
       dsepdt  = -df_tt * ye
       dsepdd  = -df_dt * x
       dsepda  = ytot1 * (ye * df_dt * din - sele)
       dsepdz  = -ytot1 * (ye * df_dt * den  + df_t)


       eele    = ye*free + temp * sele
       deepdt  = temp * dsepdt
       deepdd  = x * df_d + temp * dsepdd
       deepda  = -ye * ytot1 * (free +  df_d * din) + temp * dsepda
       deepdz  = ytot1* (free + ye * df_d * den) + temp * dsepdz




c..coulomb section:
c..initialize


        pcoul    = 0.0d0
        dpcouldd = 0.0d0
        dpcouldt = 0.0d0
        dpcoulda = 0.0d0
        dpcouldz = 0.0d0
        ecoul    = 0.0d0
        decouldd = 0.0d0
        decouldt = 0.0d0
        decoulda = 0.0d0
        decouldz = 0.0d0
        scoul    = 0.0d0
        dscouldd = 0.0d0
        dscouldt = 0.0d0
        dscoulda = 0.0d0
        dscouldz = 0.0d0


c..uniform background corrections only 
c..from yakovlev & shalybkov 1989 
c..lami is the average ion seperation
c..plasg is the plasma coupling parameter
        z        = forth * pi
        s        = z * xni
        dsdd     = z * dxnidd
        dsda     = z * dxnida

        lami     = 1.0d0/s**third
        inv_lami = 1.0d0/lami
        z        = -third * lami
        lamidd   = z * dsdd/s
        lamida   = z * dsda/s

        plasg    = zbar*zbar*esqu*ktinv*inv_lami
        z        = -plasg * inv_lami 
        plasgdd  = z * lamidd
        plasgda  = z * lamida
        plasgdt  = -plasg*ktinv * kerg
        plasgdz  = 2.0d0 * plasg/zbar


c..yakovlev & shalybkov 1989 equations 82, 85, 86, 87
         if (plasg .ge. 1.0) then
          x        = plasg**(0.25d0) 
          y        = avo * ytot1 * kerg 
          ecoul    = y * temp * (a1*plasg + b1*x + c1/x + d1)
          pcoul    = third * den * ecoul
          scoul    = -y * (3.0d0*b1*x - 5.0d0*c1/x
     1              + d1 * (log(plasg) - 1.0d0) - e1)

          y        = avo*ytot1*kt*(a1 + 0.25d0/plasg*(b1*x - c1/x))
          decouldd = y * plasgdd 
          decouldt = y * plasgdt + ecoul/temp
          decoulda = y * plasgda - ecoul/abar
          decouldz = y * plasgdz

          y        = third * den
          dpcouldd = third * ecoul + y*decouldd
          dpcouldt = y * decouldt
          dpcoulda = y * decoulda
          dpcouldz = y * decouldz


          y        = -avo*kerg/(abar*plasg)*(0.75d0*b1*x+1.25d0*c1/x+d1)
          dscouldd = y * plasgdd
          dscouldt = y * plasgdt
          dscoulda = y * plasgda - scoul/abar
          dscouldz = y * plasgdz


c..yakovlev & shalybkov 1989 equations 102, 103, 104
         else if (plasg .lt. 1.0) then
          x        = plasg*sqrt(plasg)
          y        = plasg**b2
          z        = c2 * x - third * a2 * y
          pcoul    = -pion * z
          ecoul    = 3.0d0 * pcoul/den
          scoul    = -avo/abar*kerg*(c2*x -a2*(b2-1.0d0)/b2*y)

          s        = 1.5d0*c2*x/plasg - third*a2*b2*y/plasg
          dpcouldd = -dpiondd*z - pion*s*plasgdd
          dpcouldt = -dpiondt*z - pion*s*plasgdt
          dpcoulda = -dpionda*z - pion*s*plasgda
          dpcouldz = -dpiondz*z - pion*s*plasgdz

          s        = 3.0d0/den
          decouldd = s * dpcouldd - ecoul/den
          decouldt = s * dpcouldt
          decoulda = s * dpcoulda
          decouldz = s * dpcouldz

          s        = -avo*kerg/(abar*plasg)*(1.5d0*c2*x-a2*(b2-1.0d0)*y)
          dscouldd = s * plasgdd
          dscouldt = s * plasgdt
          dscoulda = s * plasgda - scoul/abar
          dscouldz = s * plasgdz
         end if



c..bomb proof
        x   = prad + pion + pele + pcoul
        if (x .le. 0.0) then

c         write(6,*) 
c         write(6,*) 'coulomb corrections are causing a negative pressure'
c         write(6,*) 'setting all coulomb corrections to zero'
c         write(6,*) 

         pcoul    = 0.0d0
         dpcouldd = 0.0d0
         dpcouldt = 0.0d0
         dpcoulda = 0.0d0
         dpcouldz = 0.0d0
         ecoul    = 0.0d0
         decouldd = 0.0d0
         decouldt = 0.0d0
         decoulda = 0.0d0
         decouldz = 0.0d0
         scoul    = 0.0d0
         dscouldd = 0.0d0
         dscouldt = 0.0d0
         dscoulda = 0.0d0
         dscouldz = 0.0d0
        end if





c..sum all the components
       pres    = prad + pion + pele + pcoul
       ener    = erad + eion + eele + ecoul
       entr    = srad + sion + sele + scoul

       dpresdd = dpraddd + dpiondd + dpepdd + dpcouldd 
       dpresdt = dpraddt + dpiondt + dpepdt + dpcouldt
       dpresda = dpradda + dpionda + dpepda + dpcoulda
       dpresdz = dpraddz + dpiondz + dpepdz + dpcouldz

       denerdd = deraddd + deiondd + deepdd + decouldd
       denerdt = deraddt + deiondt + deepdt + decouldt
       denerda = deradda + deionda + deepda + decoulda
       denerdz = deraddz + deiondz + deepdz + decouldz

       dentrdd = dsraddd + dsiondd + dsepdd + dscouldd
       dentrdt = dsraddt + dsiondt + dsepdt + dscouldt
       dentrda = dsradda + dsionda + dsepda + dscoulda
       dentrdz = dsraddz + dsiondz + dsepdz + dscouldz

c..the temperature and density exponents (c&g 9.81 9.82) 
c..the specific heat at constant volume (c&g 9.92)
c..the third adiabatic exponent (c&g 9.93)
c..the first adiabatic exponent (c&g 9.97) 
c..the second adiabatic exponent (c&g 9.105)
c..the specific heat at constant pressure (c&g 9.98) 
c..and relativistic formula for the sound speed (c&g 14.29)
       zz    = pres*deni
       zzi   = den/pres
       chit  = temp/pres * dpresdt
       chid  = dpresdd*zzi
       cv    = denerdt
       x     = zz * chit/(temp * cv)
       gam3  = x + 1.0d0
       gam1  = chit*x + chid
       nabad = x/gam1
       gam2  = 1.0d0/(1.0d0 - nabad)
       cp    = cv * gam1/chid
       z     = 1.0d0 + (ener + light2)*zzi
       sound = clight * sqrt(gam1/z)


c..maxwell relations; each is zero if the consistency is perfect
C        x   = den * den
C        dse = temp*dentrdt/denerdt - 1.0d0
C        dpe = (denerdd*x + temp*dpresdt)/pres - 1.0d0
C        dsp = -dentrdd*x/dpresdt - 1.0d0


c..store this row
        ptot_row(j)   = pres
C         dpt_row(j)    = dpresdt
C         dpd_row(j)    = dpresdd
C         dpa_row(j)    = dpresda   
C         dpz_row(j)    = dpresdz

        etot_row(j)   = ener
        det_row(j)    = denerdt
C         ded_row(j)    = denerdd
C         dea_row(j)    = denerda   
C         dez_row(j)    = denerdz

C         stot_row(j)   = entr 
C         dst_row(j)    = dentrdt
C         dsd_row(j)    = dentrdd
C         dsa_row(j)    = dentrda        
C         dsz_row(j)    = dentrdz

C         prad_row(j)   = prad
C         erad_row(j)   = erad
C         srad_row(j)   = srad 

C         pion_row(j)   = pion
C         eion_row(j)   = eion
C         sion_row(j)   = sion 
C         xni_row(j)    = xni

C         pele_row(j)   = pele
C         ppos_row(j)   = 0.0d0
C         dpept_row(j)  = dpepdt
C         dpepd_row(j)  = dpepdd
C         dpepa_row(j)  = dpepda  
C         dpepz_row(j)  = dpepdz

C         eele_row(j)   = eele
C         epos_row(j)   = 0.0d0
C         deept_row(j)  = deepdt
C         deepd_row(j)  = deepdd
C         deepa_row(j)  = deepda   
C         deepz_row(j)  = deepdz

C         sele_row(j)   = sele 
C         spos_row(j)   = 0.0d0
C         dsept_row(j)  = dsepdt 
C         dsepd_row(j)  = dsepdd 
C         dsepa_row(j)  = dsepda        
C         dsepz_row(j)  = dsepdz

C         xnem_row(j)   = xnem
C         xne_row(j)    = xnefer
C         dxnet_row(j)  = dxnedt
C         dxned_row(j)  = dxnedd
C         dxnea_row(j)  = dxneda
C         dxnez_row(j)  = dxnedz
C         xnp_row(j)    = 0.0d0

C         etaele_row(j) = etaele
C         detat_row(j)  = detadt
C         detad_row(j)  = detadd
C         detaa_row(j)  = detada
C         detaz_row(j)  = detadz
C         etapos_row(j) = 0.0d0

C         pcou_row(j)   = pcoul
C         ecou_row(j)   = ecoul
C         scou_row(j)   = scoul 
C         plasg_row(j)  = plasg

C         dse_row(j)    = dse
C         dpe_row(j)    = dpe
C         dsp_row(j)    = dsp

C         cv_row(j)     = cv
C         cp_row(j)     = cp
C         gam1_row(j)   = gam1
C         gam2_row(j)   = gam2
C         gam3_row(j)   = gam3
        cs_row(j)     = sound

c..end of vectorization loop
      enddo
      return
      end

      subroutine helmeos2(temp, den, abar, zbar, pres, ener,
     &     denerdt, sound, eosfail, flagcoulomb)
      implicit none
CC      save
cc      include 'vector_eos.dek'
      logical          eosfail
      double precision flagcoulomb

c..given a temperature temp [K], density den [g/cm**3], and a composition 
c..characterized by abar and zbar, this routine returns most of the other 
c..thermodynamic quantities. of prime interest is the pressure [erg/cm**3], 
c..specific thermal energy [erg/gr], the entropy [erg/g/K], along with 
c..their derivatives with respect to temperature, density, abar, and zbar.
c..other quantites such the normalized chemical potential eta (plus its
c..derivatives), number density of electrons and positron pair (along 
c..with their derivatives), adiabatic indices, specific heats, and 
c..relativistically correct sound speed are also returned.
c..
c..this routine assumes planckian photons, an ideal gas of ions,
c..and an electron-positron gas with an arbitrary degree of relativity
c..and degeneracy. interpolation in a table of the helmholtz free energy
c..is used to return the electron-positron thermodynamic quantities.
c..all other derivatives are analytic.
c..
c..references: cox & giuli chapter 24 ; timmes & swesty apj 1999


c..declare
      double precision pi,amu,kerg,clight,avo,qe,h,ssol,asol
      parameter       (pi      = 3.1415926535897932384d0,
     1                  amu    = 1.6605402d-24,
     2                  kerg   = 1.380658d-16,
     3                  clight = 2.99792458d10, 
     4                  avo    = 6.0221367d23,
     5                  qe     = 4.8032068d-10,  
     6                  h      = 6.6260755d-27,
     7                  ssol   = 5.67051d-5,
     8                  asol   = 4.0d0 * ssol / clight)

      integer          i,j
      double precision x,y,zz,zzi,deni,tempi,xni,dxnidd,dxnida,
     1                 dpepdt,dpepdd,deepdt,deepdd,dsepdd,dsepdt,
     2                 dpraddd,dpraddt,deraddd,deraddt,dpiondd,dpiondt,
     3                 deiondd,deiondt,dsraddd,dsraddt,dsiondd,dsiondt,
     4                 dse,dpe,dsp,kt,ktinv,prad,erad,srad,pion,eion,
     5                 sion,xnem,pele,eele,sele,pres,ener,entr,dpresdd,
     6                 dpresdt,denerdd,denerdt,dentrdd,dentrdt,cv,cp,
     7                 gam1,gam2,gam3,chit,chid,nabad,sound,etaele,
     8                 detadt,detadd,xnefer,dxnedt,dxnedd,s,
     9                 temp,den,abar,zbar,ytot1,ye,
     &                 sioncon,forth,forpi,kergavo,ikavo,asoli3,light2

      parameter        (sioncon = (2.0d0 * pi * amu * kerg)/(h*h),
     1                  forth   = 4.0d0/3.0d0,
     2                  forpi   = 4.0d0 * pi,
     3                  kergavo = kerg * avo, 
     4                  ikavo   = 1.0d0/kergavo,
     5                  asoli3  = asol/3.0d0,
     6                  light2  = clight * clight)

c..for the abar derivatives
      double precision dpradda,deradda,dsradda,
     1                 dpionda,deionda,dsionda,
     2                 dpepda,deepda,dsepda,
     3                 dpresda,denerda,dentrda,
     4                 detada,dxneda


c..for the zbar derivatives
      double precision dpraddz,deraddz,dsraddz,
     1                 dpiondz,deiondz,dsiondz,
     2                 dpepdz,deepdz,dsepdz,
     3                 dpresdz,denerdz,dentrdz,
     4                 detadz,dxnedz


c..for the interpolations
      integer          iat,jat
      double precision tlo,thi,tstp,tstpi,dlo,dhi,dstp,dstpi,
     1                 tsav,dsav,free,df_d,df_t,df_dd,df_tt,df_dt
      double precision dth,dt2,dti,dt2i,dd,dd2,ddi,dd2i,xt,xd,mxt,mxd,
     1                 si0t,si1t,si2t,si0mt,si1mt,si2mt,
     2                 si0d,si1d,si2d,si0md,si1md,si2md,
     3                 dsi0t,dsi1t,dsi2t,dsi0mt,dsi1mt,dsi2mt,
     4                 dsi0d,dsi1d,dsi2d,dsi0md,dsi1md,dsi2md,
     5                 ddsi0t,ddsi1t,ddsi2t,ddsi0mt,ddsi1mt,ddsi2mt,
     6                 ddsi0d,ddsi1d,ddsi2d,ddsi0md,ddsi1md,ddsi2md,
     7                 z,psi0,dpsi0,ddpsi0,psi1,dpsi1,ddpsi1,psi2,
     8                 dpsi2,ddpsi2,din,h5,fi(36),
     9                 xpsi0,xdpsi0,xpsi1,xdpsi1,h3,
     1                 w0t,w1t,w2t,w0mt,w1mt,w2mt,
     2                 w0d,w1d,w2d,w0md,w1md,w2md

c..for the coulomb corrections
      double precision dsdd,dsda,lami,inv_lami,lamida,lamidd,
     1                 plasg,plasgdd,plasgdt,plasgda,plasgdz,
     1                 a1,b1,c1,d1,e1,a2,b2,c2,
     3                 ecoul,decouldd,decouldt,decoulda,decouldz,
     4                 pcoul,dpcouldd,dpcouldt,dpcoulda,dpcouldz,
     5                 scoul,dscouldd,dscouldt,dscoulda,dscouldz,
     6                 tmelt,tfermi,rhocond,z2,x1,x2,third,esqu
      parameter        (a1    = -0.898004d0, 
     1                  b1    =  0.96786d0, 
     2                  c1    =  0.220703d0, 
     3                  d1    = -0.86097d0,
     4                  e1    =  2.5269d0, 
     5                  a2    =  0.29561d0, 
     6                  b2    =  1.9885d0,    
     7                  c2    =  0.288675d0,
     8                  third = 1.0d0/3.0d0,
     9                  esqu  = qe * qe)

      save
c..for the tables, in general
      integer          imax,jmax
c      parameter        (imax = 211, jmax = 71)
      parameter        (imax = 271, jmax = 101)
      double precision d(imax),t(jmax)

c..for the helmholtz free energy tables
      double precision f(imax,jmax),fd(imax,jmax),
     1                 ft(imax,jmax),fdd(imax,jmax),ftt(imax,jmax),
     2                 fdt(imax,jmax),fddt(imax,jmax),fdtt(imax,jmax),
     3                 fddtt(imax,jmax)

c..for the pressure derivative with density ables
      double precision dpdf(imax,jmax),dpdfd(imax,jmax),
     1                 dpdft(imax,jmax),dpdfdd(imax,jmax),
     2                 dpdftt(imax,jmax),dpdfdt(imax,jmax)

c..for chemical potential tables
      double precision ef(imax,jmax),efd(imax,jmax),
     1                 eft(imax,jmax),efdd(imax,jmax),eftt(imax,jmax),
     2                 efdt(imax,jmax)

c..for the number density tables
      double precision xf(imax,jmax),xfd(imax,jmax),
     1                 xft(imax,jmax),xfdd(imax,jmax),xftt(imax,jmax),
     2                 xfdt(imax,jmax)

c..for storing the differences
      double precision dt_sav(jmax),dt2_sav(jmax),
     1                 dti_sav(jmax),dt2i_sav(jmax),
     2                 dd_sav(imax),dd2_sav(imax),
     3                 ddi_sav(imax),dd2i_sav(imax)

c..for initialization
      integer          ifirst
      data             ifirst/0/ 
CC      save ifirst

c..quintic hermite polynomial statement functions
c..psi0 and its derivatives
      psi0(z)   = z**3 * ( z * (-6.0d0*z + 15.0d0) -10.0d0) + 1.0d0
      dpsi0(z)  = z**2 * ( z * (-30.0d0*z + 60.0d0) - 30.0d0)
      ddpsi0(z) = z* ( z*( -120.0d0*z + 180.0d0) -60.0d0)

c..psi1 and its derivatives
      psi1(z)   = z* ( z**2 * ( z * (-3.0d0*z + 8.0d0) - 6.0d0) + 1.0d0)
      dpsi1(z)  = z*z * ( z * (-15.0d0*z + 32.0d0) - 18.0d0) +1.0d0
      ddpsi1(z) = z * (z * (-60.0d0*z + 96.0d0) -36.0d0)

c..psi2  and its derivatives
      psi2(z)   = 0.5d0*z*z*( z* ( z * (-z + 3.0d0) - 3.0d0) + 1.0d0)
      dpsi2(z)  = 0.5d0*z*( z*(z*(-5.0d0*z + 12.0d0) - 9.0d0) + 2.0d0)
      ddpsi2(z) = 0.5d0*(z*( z * (-20.0d0*z + 36.0d0) - 18.0d0) + 2.0d0)

c..biquintic hermite polynomial statement function
      h5(i,j,w0t,w1t,w2t,w0mt,w1mt,w2mt,w0d,w1d,w2d,w0md,w1md,w2md)=
     1       fi(1)  *w0d*w0t   + fi(2)  *w0md*w0t
     2     + fi(3)  *w0d*w0mt  + fi(4)  *w0md*w0mt
     3     + fi(5)  *w0d*w1t   + fi(6)  *w0md*w1t
     4     + fi(7)  *w0d*w1mt  + fi(8)  *w0md*w1mt
     5     + fi(9)  *w0d*w2t   + fi(10) *w0md*w2t
     6     + fi(11) *w0d*w2mt  + fi(12) *w0md*w2mt
     7     + fi(13) *w1d*w0t   + fi(14) *w1md*w0t
     8     + fi(15) *w1d*w0mt  + fi(16) *w1md*w0mt
     9     + fi(17) *w2d*w0t   + fi(18) *w2md*w0t
     &     + fi(19) *w2d*w0mt  + fi(20) *w2md*w0mt
     1     + fi(21) *w1d*w1t   + fi(22) *w1md*w1t
     2     + fi(23) *w1d*w1mt  + fi(24) *w1md*w1mt
     3     + fi(25) *w2d*w1t   + fi(26) *w2md*w1t
     4     + fi(27) *w2d*w1mt  + fi(28) *w2md*w1mt
     5     + fi(29) *w1d*w2t   + fi(30) *w1md*w2t
     6     + fi(31) *w1d*w2mt  + fi(32) *w1md*w2mt
     7     + fi(33) *w2d*w2t   + fi(34) *w2md*w2t
     8     + fi(35) *w2d*w2mt  + fi(36) *w2md*w2mt



c..cubic hermite polynomial statement functions
c..psi0 & derivatives
      xpsi0(z)  = z * z * (2.0d0*z - 3.0d0) + 1.0
      xdpsi0(z) = z * (6.0d0*z - 6.0d0)

c..psi1 & derivatives
      xpsi1(z)  = z * ( z * (z - 2.0d0) + 1.0d0)
      xdpsi1(z) = z * (3.0d0*z - 4.0d0) + 1.0d0


c..bicubic hermite polynomial statement function
      h3(i,j,w0t,w1t,w0mt,w1mt,w0d,w1d,w0md,w1md) = 
     1       fi(1)  *w0d*w0t   +  fi(2)  *w0md*w0t 
     2     + fi(3)  *w0d*w0mt  +  fi(4)  *w0md*w0mt
     3     + fi(5)  *w0d*w1t   +  fi(6)  *w0md*w1t 
     4     + fi(7)  *w0d*w1mt  +  fi(8)  *w0md*w1mt
     5     + fi(9)  *w1d*w0t   +  fi(10) *w1md*w0t 
     6     + fi(11) *w1d*w0mt  +  fi(12) *w1md*w0mt
     7     + fi(13) *w1d*w1t   +  fi(14) *w1md*w1t 
     8     + fi(15) *w1d*w1mt  +  fi(16) *w1md*w1mt



c..popular format statements
01    format(1x,5(a,1pe11.3))
02    format(1x,a,1p4e16.8)
03    format(1x,4(a,1pe11.3))
04    format(1x,4(a,i4))


c..do this stuff once
      if (ifirst .eq. 0) then
       ifirst = 1

c..open the table
c      open(unit=2,file='helm_table_large.dat',status='old')
c       open(unit=2,file='heos/helm_table_large.dat',status='old')
       open(unit=2,file='data.code/helm_table_large.dat',status='old')

c..read the helmholtz free energy table
c       tlo   = 4.0d0
c       thi   = 11.0d0
c       tstp  = (thi - tlo)/float(jmax-1)
c       tstpi = 1.0d0/tstp
c       dlo   = -10.0d0
c       dhi   = 11.0d0
c
       tlo   = 3.0d0
       thi   = 13.0d0
       tstp  = (thi - tlo)/float(jmax-1)
       tstpi = 1.0d0/tstp
       dlo   = -12.0d0
       dhi   = 15.0d0
c
       dstp  = (dhi - dlo)/float(imax-1)
       dstpi = 1.0d0/dstp
       do j=1,jmax
        tsav = tlo + (j-1)*tstp
        t(j) = 10.0d0**(tsav)
        do i=1,imax
         dsav = dlo + (i-1)*dstp
         d(i) = 10.0d0**(dsav)
         read(2,*) f(i,j),fd(i,j),ft(i,j),fdd(i,j),ftt(i,j),fdt(i,j),
     1            fddt(i,j),fdtt(i,j),fddtt(i,j)
        enddo
       enddo

c..read the pressure derivative with density table
       do j=1,jmax
        do i=1,imax
         read(2,*) dpdf(i,j),dpdfd(i,j),dpdft(i,j),dpdfdt(i,j)
        enddo
       enddo

c..read the electron chemical potential table
       do j=1,jmax
        do i=1,imax
         read(2,*) ef(i,j),efd(i,j),eft(i,j),efdt(i,j)
        enddo
       enddo

c..read the number density table
       do j=1,jmax
        do i=1,imax
         read(2,*) xf(i,j),xfd(i,j),xft(i,j),xfdt(i,j)
        enddo
       enddo

c..construct the temperature and density deltas and their inverses 
       do j=1,jmax-1
        dth          = t(j+1) - t(j)
        dt2         = dth * dth
        dti         = 1.0d0/dth
        dt2i        = 1.0d0/dt2
        dt_sav(j)   = dth
        dt2_sav(j)  = dt2
        dti_sav(j)  = dti
        dt2i_sav(j) = dt2i
       end do
       do i=1,imax-1
        dd          = d(i+1) - d(i)
        dd2         = dd * dd
        ddi         = 1.0d0/dd
        dd2i        = 1.0d0/dd2
        dd_sav(i)   = dd
        dd2_sav(i)  = dd2
        ddi_sav(i)  = ddi
        dd2i_sav(i) = dd2i
       enddo

       close(unit=2)
       write(6,*)
       write(6,*) 'finished reading eos table'
       write(6,04) 'imax=',imax,' jmax=',jmax
       write(6,03) 'temp(1)   =',t(1),' temp(jmax)   =',t(jmax)
       write(6,03) 'ye*den(1) =',d(1),' ye*den(imax) =',d(imax)
       write(6,*)

      end if



c..start of vectorization loop, normal executaion starts here
      eosfail = .false.
c      do j=jlo_eos,jhi_eos

       if (temp .le. 0.0) stop 'temp less than 0 in helmeos2'
       if (den  .le. 0.0) stop 'den less than 0 in helmeos2'

C       temp  = temp
C       den   = den
C       abar  = abar
C       zbar  = zbar
       ytot1 = 1.0d0/abar
       ye    = ytot1 * zbar

c..initialize
       deni    = 1.0d0/den
       tempi   = 1.0d0/temp 
       kt      = kerg * temp
       ktinv   = 1.0d0/kt


c..radiation section:
       prad    = asoli3 * temp * temp * temp * temp
       dpraddd = 0.0d0
       dpraddt = 4.0d0 * prad*tempi
       dpradda = 0.0d0
       dpraddz = 0.0d0

       erad    = 3.0d0 * prad*deni
       deraddd = -erad*deni
       deraddt = 3.0d0 * dpraddt*deni
       deradda = 0.0d0
       deraddz = 0.0d0


       srad    = (prad*deni + erad)*tempi
       dsraddd = (dpraddd*deni - prad*deni*deni + deraddd)*tempi
       dsraddt = (dpraddt*deni + deraddt - srad)*tempi
       dsradda = 0.0d0
       dsraddz = 0.0d0

c..ion section:
       xni     = avo * ytot1 * den
       dxnidd  = avo * ytot1
       dxnida  = -xni * ytot1

       pion    = xni * kt
       dpiondd = dxnidd * kt
       dpiondt = xni * kerg
       dpionda = dxnida * kt 
       dpiondz = 0.0d0

       eion    = 1.5d0 * pion*deni
       deiondd = (1.5d0 * dpiondd - eion)*deni
       deiondt = 1.5d0 * dpiondt*deni
       deionda = 1.5d0 * dpionda*deni
       deiondz = 0.0d0
    

       x       = abar*abar*sqrt(abar) * deni/avo
       s       = sioncon * temp
       z       = x * s * sqrt(s)
       y       = log(z)
       sion    = (pion*deni + eion)*tempi + kergavo * ytot1 * y
       dsiondd = (dpiondd*deni - pion*deni*deni + deiondd)*tempi
     1            - kergavo * deni * ytot1
       dsiondt = (dpiondt*deni + deiondt)*tempi - 
     1           (pion*deni + eion) * tempi*tempi 
     2           + 1.5d0 * kergavo * tempi*ytot1
       x       = avo*kerg/abar
       dsionda = (dpionda*deni + deionda)*tempi 
     1           + kergavo*ytot1*ytot1* (2.5d0 - y)
       dsiondz = 0.0d0


c..electron-positron section:
c..assume complete ionization 
       xnem    = xni * zbar

c..enter the table with ye*den
       din = ye*den

c..bomb proof the input
       if (temp .gt. t(jmax)) then
        write(6,01) 'temp=',temp,' t(jmax)=',t(jmax)
        write(6,*) 'temp too hot, off grid'       
        write(6,*) 'setting eosfail to true and returning'
        call flush(6)
        eosfail = .true.
        return
       end if
       if (temp .lt. t(1)) then
        write(6,01) 'temp=',temp,' t(1)=',t(1)
        write(6,*) 'temp too cold, off grid'
        write(6,*) 'setting eosfail to true and returning'
        call flush(6)
        eosfail = .true.
        return
       end if
       if (din  .gt. d(imax)) then
        write(6,01) 'den*ye=',din,' d(imax)=',d(imax)
        write(6,*) 'ye*den too big, off grid'
        write(6,*) 'setting eosfail to true and returning'
        call flush(6)
        eosfail = .true.
        return
       end if
       if (din  .lt. d(1)) then
        write(6,01) 'ye*den=',din,' d(1)=',d(1)
        write(6,*) 'ye*den too small, off grid'
        write(6,*) 'setting eosfail to true and returning'
        call flush(6)
        eosfail = .true.
        return
       end if

c..hash locate this temperature and density
       jat = int((log10(temp) - tlo)*tstpi) + 1
       jat = max(1,min(jat,jmax-1))
       iat = int((log10(din) - dlo)*dstpi) + 1
       iat = max(1,min(iat,imax-1))


c..access the table locations only once
       fi(1)  = f(iat,jat)
       fi(2)  = f(iat+1,jat)
       fi(3)  = f(iat,jat+1)
       fi(4)  = f(iat+1,jat+1)
       fi(5)  = ft(iat,jat)
       fi(6)  = ft(iat+1,jat)
       fi(7)  = ft(iat,jat+1)
       fi(8)  = ft(iat+1,jat+1)
       fi(9)  = ftt(iat,jat)
       fi(10) = ftt(iat+1,jat)
       fi(11) = ftt(iat,jat+1)
       fi(12) = ftt(iat+1,jat+1)
       fi(13) = fd(iat,jat)
       fi(14) = fd(iat+1,jat)
       fi(15) = fd(iat,jat+1)
       fi(16) = fd(iat+1,jat+1)
       fi(17) = fdd(iat,jat)
       fi(18) = fdd(iat+1,jat)
       fi(19) = fdd(iat,jat+1)
       fi(20) = fdd(iat+1,jat+1)
       fi(21) = fdt(iat,jat)
       fi(22) = fdt(iat+1,jat)
       fi(23) = fdt(iat,jat+1)
       fi(24) = fdt(iat+1,jat+1)
       fi(25) = fddt(iat,jat)
       fi(26) = fddt(iat+1,jat)
       fi(27) = fddt(iat,jat+1)
       fi(28) = fddt(iat+1,jat+1)
       fi(29) = fdtt(iat,jat)
       fi(30) = fdtt(iat+1,jat)
       fi(31) = fdtt(iat,jat+1)
       fi(32) = fdtt(iat+1,jat+1)
       fi(33) = fddtt(iat,jat)
       fi(34) = fddtt(iat+1,jat)
       fi(35) = fddtt(iat,jat+1)
       fi(36) = fddtt(iat+1,jat+1)


c..various differences
       xt  = max( (temp - t(jat))*dti_sav(jat), 0.0d0)
       xd  = max( (din - d(iat))*ddi_sav(iat), 0.0d0)
       mxt = 1.0d0 - xt
       mxd = 1.0d0 - xd

c..the six density and six temperature basis functions
       si0t =   psi0(xt)
       si1t =   psi1(xt)*dt_sav(jat)
       si2t =   psi2(xt)*dt2_sav(jat)

       si0mt =  psi0(mxt)
       si1mt = -psi1(mxt)*dt_sav(jat)
       si2mt =  psi2(mxt)*dt2_sav(jat)

       si0d =   psi0(xd)
       si1d =   psi1(xd)*dd_sav(iat)
       si2d =   psi2(xd)*dd2_sav(iat)

       si0md =  psi0(mxd)
       si1md = -psi1(mxd)*dd_sav(iat)
       si2md =  psi2(mxd)*dd2_sav(iat)

c..derivatives of the weight functions
       dsi0t =   dpsi0(xt)*dti_sav(jat)
       dsi1t =   dpsi1(xt)
       dsi2t =   dpsi2(xt)*dt_sav(jat)

       dsi0mt = -dpsi0(mxt)*dti_sav(jat)
       dsi1mt =  dpsi1(mxt)
       dsi2mt = -dpsi2(mxt)*dt_sav(jat)

       dsi0d =   dpsi0(xd)*ddi_sav(iat)
       dsi1d =   dpsi1(xd)
       dsi2d =   dpsi2(xd)*dd_sav(iat)

       dsi0md = -dpsi0(mxd)*ddi_sav(iat)
       dsi1md =  dpsi1(mxd)
       dsi2md = -dpsi2(mxd)*dd_sav(iat)

c..second derivatives of the weight functions
       ddsi0t =   ddpsi0(xt)*dt2i_sav(jat)
       ddsi1t =   ddpsi1(xt)*dti_sav(jat)
       ddsi2t =   ddpsi2(xt)

       ddsi0mt =  ddpsi0(mxt)*dt2i_sav(jat)
       ddsi1mt = -ddpsi1(mxt)*dti_sav(jat)
       ddsi2mt =  ddpsi2(mxt)

c       ddsi0d =   ddpsi0(xd)*dd2i_sav(iat)
c       ddsi1d =   ddpsi1(xd)*ddi_sav(iat)
c       ddsi2d =   ddpsi2(xd)

c       ddsi0md =  ddpsi0(mxd)*dd2i_sav(iat)
c       ddsi1md = -ddpsi1(mxd)*ddi_sav(iat)
c       ddsi2md =  ddpsi2(mxd)


c..the free energy
       free  = h5(iat,jat,
     1         si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt,
     2         si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

c..derivative with respect to density
       df_d  = h5(iat,jat,
     1         si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt,
     2         dsi0d,  dsi1d,  dsi2d,  dsi0md,  dsi1md,  dsi2md)

c..derivative with respect to temperature
       df_t = h5(iat,jat,
     1         dsi0t,  dsi1t,  dsi2t,  dsi0mt,  dsi1mt,  dsi2mt,
     2         si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

c..derivative with respect to density**2
c       df_dd = h5(iat,jat,
c     1         si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt,
c     2         ddsi0d, ddsi1d, ddsi2d, ddsi0md, ddsi1md, ddsi2md)

c..derivative with respect to temperature**2
       df_tt = h5(iat,jat,
     1       ddsi0t, ddsi1t, ddsi2t, ddsi0mt, ddsi1mt, ddsi2mt,
     2         si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

c..derivative with respect to temperature and density
       df_dt = h5(iat,jat,
     1         dsi0t,  dsi1t,  dsi2t,  dsi0mt,  dsi1mt,  dsi2mt,
     2         dsi0d,  dsi1d,  dsi2d,  dsi0md,  dsi1md,  dsi2md)



c..now get the pressure derivative with density, chemical potential, and 
c..electron positron number densities
c..get the interpolation weight functions
       si0t   =  xpsi0(xt)
       si1t   =  xpsi1(xt)*dt_sav(jat)

       si0mt  =  xpsi0(mxt)
       si1mt  =  -xpsi1(mxt)*dt_sav(jat)

       si0d   =  xpsi0(xd)
       si1d   =  xpsi1(xd)*dd_sav(iat)

       si0md  =  xpsi0(mxd)
       si1md  =  -xpsi1(mxd)*dd_sav(iat)


c..derivatives of weight functions
       dsi0t  = xdpsi0(xt)*dti_sav(jat)
       dsi1t  = xdpsi1(xt)

       dsi0mt = -xdpsi0(mxt)*dti_sav(jat)
       dsi1mt = xdpsi1(mxt)

       dsi0d  = xdpsi0(xd)*ddi_sav(iat)
       dsi1d  = xdpsi1(xd)

       dsi0md = -xdpsi0(mxd)*ddi_sav(iat)
       dsi1md = xdpsi1(mxd)


c..look in the pressure derivative only once
       fi(1)  = dpdf(iat,jat)
       fi(2)  = dpdf(iat+1,jat)
       fi(3)  = dpdf(iat,jat+1)
       fi(4)  = dpdf(iat+1,jat+1)
       fi(5)  = dpdft(iat,jat)
       fi(6)  = dpdft(iat+1,jat)
       fi(7)  = dpdft(iat,jat+1)
       fi(8)  = dpdft(iat+1,jat+1)
       fi(9)  = dpdfd(iat,jat)
       fi(10) = dpdfd(iat+1,jat)
       fi(11) = dpdfd(iat,jat+1)
       fi(12) = dpdfd(iat+1,jat+1)
       fi(13) = dpdfdt(iat,jat)
       fi(14) = dpdfdt(iat+1,jat)
       fi(15) = dpdfdt(iat,jat+1)
       fi(16) = dpdfdt(iat+1,jat+1)

c..pressure derivative with density
       dpepdd  = h3(iat,jat,
     1                 si0t,   si1t,   si0mt,   si1mt,
     2                 si0d,   si1d,   si0md,   si1md)
       dpepdd  = max(ye * dpepdd,0.0d0)



c..look in the electron chemical potential table only once
       fi(1)  = ef(iat,jat)
       fi(2)  = ef(iat+1,jat)
       fi(3)  = ef(iat,jat+1)
       fi(4)  = ef(iat+1,jat+1)
       fi(5)  = eft(iat,jat)
       fi(6)  = eft(iat+1,jat)
       fi(7)  = eft(iat,jat+1)
       fi(8)  = eft(iat+1,jat+1)
       fi(9)  = efd(iat,jat)
       fi(10) = efd(iat+1,jat)
       fi(11) = efd(iat,jat+1)
       fi(12) = efd(iat+1,jat+1)
       fi(13) = efdt(iat,jat)
       fi(14) = efdt(iat+1,jat)
       fi(15) = efdt(iat,jat+1)
       fi(16) = efdt(iat+1,jat+1)


c..electron chemical potential etaele
       etaele  = h3(iat,jat,
     1               si0t,   si1t,   si0mt,   si1mt,
     2               si0d,   si1d,   si0md,   si1md)

c..derivative with respect to density
       x       = h3(iat,jat,
     1               si0t,   si1t,   si0mt,   si1mt,
     2              dsi0d,  dsi1d,  dsi0md,  dsi1md)
       detadd  = ye * x

c..derivative with respect to temperature
       detadt  = h3(iat,jat,
     1              dsi0t,  dsi1t,  dsi0mt,  dsi1mt,
     2               si0d,   si1d,   si0md,   si1md)

c..derivative with respect to abar and zbar
      detada = -x * din * ytot1
      detadz =  x * den * ytot1



c..look in the number density table only once
       fi(1)  = xf(iat,jat)
       fi(2)  = xf(iat+1,jat)
       fi(3)  = xf(iat,jat+1)
       fi(4)  = xf(iat+1,jat+1)
       fi(5)  = xft(iat,jat)
       fi(6)  = xft(iat+1,jat)
       fi(7)  = xft(iat,jat+1)
       fi(8)  = xft(iat+1,jat+1)
       fi(9)  = xfd(iat,jat)
       fi(10) = xfd(iat+1,jat)
       fi(11) = xfd(iat,jat+1)
       fi(12) = xfd(iat+1,jat+1)
       fi(13) = xfdt(iat,jat)
       fi(14) = xfdt(iat+1,jat)
       fi(15) = xfdt(iat,jat+1)
       fi(16) = xfdt(iat+1,jat+1)

c..electron + positron number densities
      xnefer   = h3(iat,jat,
     1               si0t,   si1t,   si0mt,   si1mt,
     2               si0d,   si1d,   si0md,   si1md)

c..derivative with respect to density
      x        = h3(iat,jat,
     1               si0t,   si1t,   si0mt,   si1mt,
     2              dsi0d,  dsi1d,  dsi0md,  dsi1md)
      x = max(x,0.0d0)
      dxnedd   = ye * x

c..derivative with respect to temperature
      dxnedt   = h3(iat,jat,
     1              dsi0t,  dsi1t,  dsi0mt,  dsi1mt,
     2               si0d,   si1d,   si0md,   si1md)

c..derivative with respect to abar and zbar
      dxneda = -x * din * ytot1
      dxnedz =  x  * den * ytot1
       


c..the desired electron-positron thermodynamic quantities

c..dpepdd at high temperatures and low densities is below the
c..floating point limit of the subtraction of two large terms.
c..since dpresdd doesn't enter the maxwell relations at all, use the
c..bicubic interpolation done above instead of this one
       x       = din * din
       pele    = x * df_d
       dpepdt  = x * df_dt
c       dpepdd  = ye * (x * df_dd + 2.0d0 * din * df_d)
       s       = dpepdd/ye - 2.0d0 * din * df_d
       dpepda  = -ytot1 * (2.0d0 * pele + s * din)
       dpepdz  = den*ytot1*(2.0d0 * din * df_d  +  s)


       x       = ye * ye
       sele    = -df_t * ye
       dsepdt  = -df_tt * ye
       dsepdd  = -df_dt * x
       dsepda  = ytot1 * (ye * df_dt * din - sele)
       dsepdz  = -ytot1 * (ye * df_dt * den  + df_t)


       eele    = ye*free + temp * sele
       deepdt  = temp * dsepdt
       deepdd  = x * df_d + temp * dsepdd
       deepda  = -ye * ytot1 * (free +  df_d * din) + temp * dsepda
       deepdz  = ytot1* (free + ye * df_d * den) + temp * dsepdz




c..coulomb section:
c..initialize


        pcoul    = 0.0d0
        dpcouldd = 0.0d0
        dpcouldt = 0.0d0
        dpcoulda = 0.0d0
        dpcouldz = 0.0d0
        ecoul    = 0.0d0
        decouldd = 0.0d0
        decouldt = 0.0d0
        decoulda = 0.0d0
        decouldz = 0.0d0
        scoul    = 0.0d0
        dscouldd = 0.0d0
        dscouldt = 0.0d0
        dscoulda = 0.0d0
        dscouldz = 0.0d0


c..uniform background corrections only 
c..from yakovlev & shalybkov 1989 
c..lami is the average ion seperation
c..plasg is the plasma coupling parameter
        z        = forth * pi
        s        = z * xni
        dsdd     = z * dxnidd
        dsda     = z * dxnida

        lami     = 1.0d0/s**third
        inv_lami = 1.0d0/lami
        z        = -third * lami
        lamidd   = z * dsdd/s
        lamida   = z * dsda/s

        plasg    = zbar*zbar*esqu*ktinv*inv_lami
        z        = -plasg * inv_lami 
        plasgdd  = z * lamidd
        plasgda  = z * lamida
        plasgdt  = -plasg*ktinv * kerg
        plasgdz  = 2.0d0 * plasg/zbar


c..yakovlev & shalybkov 1989 equations 82, 85, 86, 87
         if (plasg .ge. 1.0) then
          x        = plasg**(0.25d0) 
          y        = avo * ytot1 * kerg 
          ecoul    = y * temp * (a1*plasg + b1*x + c1/x + d1)
          pcoul    = third * den * ecoul
          scoul    = -y * (3.0d0*b1*x - 5.0d0*c1/x
     1              + d1 * (log(plasg) - 1.0d0) - e1)

          y        = avo*ytot1*kt*(a1 + 0.25d0/plasg*(b1*x - c1/x))
          decouldd = y * plasgdd 
          decouldt = y * plasgdt + ecoul/temp
          decoulda = y * plasgda - ecoul/abar
          decouldz = y * plasgdz

          y        = third * den
          dpcouldd = third * ecoul + y*decouldd
          dpcouldt = y * decouldt
          dpcoulda = y * decoulda
          dpcouldz = y * decouldz


          y        = -avo*kerg/(abar*plasg)*(0.75d0*b1*x+1.25d0*c1/x+d1)
          dscouldd = y * plasgdd
          dscouldt = y * plasgdt
          dscoulda = y * plasgda - scoul/abar
          dscouldz = y * plasgdz


c..yakovlev & shalybkov 1989 equations 102, 103, 104
         else if (plasg .lt. 1.0) then
          x        = plasg*sqrt(plasg)
          y        = plasg**b2
          z        = c2 * x - third * a2 * y
          pcoul    = -pion * z
          ecoul    = 3.0d0 * pcoul/den
          scoul    = -avo/abar*kerg*(c2*x -a2*(b2-1.0d0)/b2*y)

          s        = 1.5d0*c2*x/plasg - third*a2*b2*y/plasg
          dpcouldd = -dpiondd*z - pion*s*plasgdd
          dpcouldt = -dpiondt*z - pion*s*plasgdt
          dpcoulda = -dpionda*z - pion*s*plasgda
          dpcouldz = -dpiondz*z - pion*s*plasgdz

          s        = 3.0d0/den
          decouldd = s * dpcouldd - ecoul/den
          decouldt = s * dpcouldt
          decoulda = s * dpcoulda
          decouldz = s * dpcouldz

          s        = -avo*kerg/(abar*plasg)*(1.5d0*c2*x-a2*(b2-1.0d0)*y)
          dscouldd = s * plasgdd
          dscouldt = s * plasgdt
          dscoulda = s * plasgda - scoul/abar
          dscouldz = s * plasgdz
         end if


c..off coulomb by tanikawa FROM
         if (flagcoulomb .eq. 0.) then
            pcoul    = 0.0d0
            dpcouldd = 0.0d0
            dpcouldt = 0.0d0
            dpcoulda = 0.0d0
            dpcouldz = 0.0d0
            ecoul    = 0.0d0
            decouldd = 0.0d0
            decouldt = 0.0d0
            decoulda = 0.0d0
            decouldz = 0.0d0
            scoul    = 0.0d0
            dscouldd = 0.0d0
            dscouldt = 0.0d0
            dscoulda = 0.0d0
            dscouldz = 0.0d0
         endif
c..off coulomb by tanikawa TO


c..bomb proof
        x   = prad + pion + pele + pcoul
        if (x .le. 0.0) then

c         write(6,*) 
c         write(6,*) 'coulomb corrections are causing a negative pressure'
c         write(6,*) 'setting all coulomb corrections to zero'
c         write(6,*) 

         pcoul    = 0.0d0
         dpcouldd = 0.0d0
         dpcouldt = 0.0d0
         dpcoulda = 0.0d0
         dpcouldz = 0.0d0
         ecoul    = 0.0d0
         decouldd = 0.0d0
         decouldt = 0.0d0
         decoulda = 0.0d0
         decouldz = 0.0d0
         scoul    = 0.0d0
         dscouldd = 0.0d0
         dscouldt = 0.0d0
         dscoulda = 0.0d0
         dscouldz = 0.0d0
        end if





c..sum all the components
       pres    = prad + pion + pele + pcoul
       ener    = erad + eion + eele + ecoul
       entr    = srad + sion + sele + scoul

       dpresdd = dpraddd + dpiondd + dpepdd + dpcouldd 
       dpresdt = dpraddt + dpiondt + dpepdt + dpcouldt
       dpresda = dpradda + dpionda + dpepda + dpcoulda
       dpresdz = dpraddz + dpiondz + dpepdz + dpcouldz

       denerdd = deraddd + deiondd + deepdd + decouldd
       denerdt = deraddt + deiondt + deepdt + decouldt
       denerda = deradda + deionda + deepda + decoulda
       denerdz = deraddz + deiondz + deepdz + decouldz

       dentrdd = dsraddd + dsiondd + dsepdd + dscouldd
       dentrdt = dsraddt + dsiondt + dsepdt + dscouldt
       dentrda = dsradda + dsionda + dsepda + dscoulda
       dentrdz = dsraddz + dsiondz + dsepdz + dscouldz

c..the temperature and density exponents (c&g 9.81 9.82) 
c..the specific heat at constant volume (c&g 9.92)
c..the third adiabatic exponent (c&g 9.93)
c..the first adiabatic exponent (c&g 9.97) 
c..the second adiabatic exponent (c&g 9.105)
c..the specific heat at constant pressure (c&g 9.98) 
c..and relativistic formula for the sound speed (c&g 14.29)
       zz    = pres*deni
       zzi   = den/pres
       chit  = temp/pres * dpresdt
       chid  = dpresdd*zzi
       cv    = denerdt
       x     = zz * chit/(temp * cv)
       gam3  = x + 1.0d0
       gam1  = chit*x + chid
       nabad = x/gam1
       gam2  = 1.0d0/(1.0d0 - nabad)
       cp    = cv * gam1/chid
       z     = 1.0d0 + (ener + light2)*zzi
       sound = clight * sqrt(gam1/z)


c..maxwell relations; each is zero if the consistency is perfect
C        x   = den * den
C        dse = temp*dentrdt/denerdt - 1.0d0
C        dpe = (denerdd*x + temp*dpresdt)/pres - 1.0d0
C        dsp = -dentrdd*x/dpresdt - 1.0d0


c..store this row
C        ptot_row(j)   = pres
C         dpt_row(j)    = dpresdt
C         dpd_row(j)    = dpresdd
C         dpa_row(j)    = dpresda   
C         dpz_row(j)    = dpresdz

C        etot_row(j)   = ener
C        det_row(j)    = denerdt
C         ded_row(j)    = denerdd
C         dea_row(j)    = denerda   
C         dez_row(j)    = denerdz

C         stot_row(j)   = entr 
C         dst_row(j)    = dentrdt
C         dsd_row(j)    = dentrdd
C         dsa_row(j)    = dentrda        
C         dsz_row(j)    = dentrdz

C         prad_row(j)   = prad
C         erad_row(j)   = erad
C         srad_row(j)   = srad 

C         pion_row(j)   = pion
C         eion_row(j)   = eion
C         sion_row(j)   = sion 
C         xni_row(j)    = xni

C         pele_row(j)   = pele
C         ppos_row(j)   = 0.0d0
C         dpept_row(j)  = dpepdt
C         dpepd_row(j)  = dpepdd
C         dpepa_row(j)  = dpepda  
C         dpepz_row(j)  = dpepdz

C         eele_row(j)   = eele
C         epos_row(j)   = 0.0d0
C         deept_row(j)  = deepdt
C         deepd_row(j)  = deepdd
C         deepa_row(j)  = deepda   
C         deepz_row(j)  = deepdz

C         sele_row(j)   = sele 
C         spos_row(j)   = 0.0d0
C         dsept_row(j)  = dsepdt 
C         dsepd_row(j)  = dsepdd 
C         dsepa_row(j)  = dsepda        
C         dsepz_row(j)  = dsepdz

C         xnem_row(j)   = xnem
C         xne_row(j)    = xnefer
C         dxnet_row(j)  = dxnedt
C         dxned_row(j)  = dxnedd
C         dxnea_row(j)  = dxneda
C         dxnez_row(j)  = dxnedz
C         xnp_row(j)    = 0.0d0

C         etaele_row(j) = etaele
C         detat_row(j)  = detadt
C         detad_row(j)  = detadd
C         detaa_row(j)  = detada
C         detaz_row(j)  = detadz
C         etapos_row(j) = 0.0d0

C         pcou_row(j)   = pcoul
C         ecou_row(j)   = ecoul
C         scou_row(j)   = scoul 
C         plasg_row(j)  = plasg

C         dse_row(j)    = dse
C         dpe_row(j)    = dpe
C         dsp_row(j)    = dsp

C         cv_row(j)     = cv
C         cp_row(j)     = cp
C         gam1_row(j)   = gam1
C         gam2_row(j)   = gam2
C         gam3_row(j)   = gam3
C        cs_row(j)     = sound

c..end of vectorization loop
C      enddo

      return
      end



