      subroutine burner(time_start,time_end,dtsav, &
                        tin,din,vin,zin,ein,xin, &
                        tout,dout,eout,xout, &
                        conserv,nok,nbad)
      include 'implno.dek'
      include 'const.dek'
      include 'network.dek'

! input:
! time_start starting time
! time_end   ending time
! dtsa v     if non-zero, timestep interval to save output
! tin        initial temperature
! din        initial density
! vin        initial velocity                                                                                           
! zin        initial position
! ein        initial internal energy
! xin        initial composition vector

! output:
! tout       final temperature
! dout       final density
! eout       final internal energy
! xout       final composition vector
! conserv    1 - sum of mass fractions
! nok        number of good time steps taken
! nbad       number of bad timesteps attempted


! declare the pass
      integer          nok,nbad
      double precision time_start,time_end,dtsav,tin,din,vin,zin,ein,xin(*), &
                       tout,dout,eout,xout(*),conserv

! local variables
      integer          i
      double precision abar,zbar,wbar,ye,xcess

! for the integration driver
      integer          kount
      double precision stptry,stpmin,ys2(abignet*nzmax), &
                       odescal,tol

! usually adequate
      parameter        (tol     = 1.0d-6, &
                        odescal = 1.0d-8)

! for very accurate integrations
!      parameter        (tol     = 1.0d-10, &
!                        odescal = 1.0d-12)


      external         aprox13,saprox13,baprox13,daprox13


!      external         forder_ma28
!      external         forder_umf
!      external         forder_y12m
!      external         forder_ludcmp
!      external         forder_leqs
!      external         forder_lapack
!      external         forder_gift
!      external         forder_biconj

!      external         rosen_ma28
!      external         rosen_umf
!      external         rosen_y12m
!      external         rosen_ludcmp
!      external         rosen_leqs
!      external         rosen_lapack
!      external         rosen_gift
!      external         rosen_biconj

      external         stifbs_ma28
!      external         stifbs_umf
!      external         stifbs_y12m
!      external         stifbs_ludcmp
!      external         stifbs_leqs
!      external         stifbs_lapack
!      external         stifbs_gift
!      external         stifbs_biconj



! load the mass fractions
      xmass(ionbeg:ionend) = xin(ionbeg:ionend)

! get abar, zbar and a few other composition variables
      call azbarapprox13(xmass,aion,zion,wion,ionmax, &
           ymass,abar,zbar,wbar,ye,xcess)



! stuff the initial conditions into ys2
       ys2(ionbeg:ionend) = ymass(ionbeg:ionend)
       ys2(iener) = ein
       ys2(itemp) = tin
       ys2(iden)  = din
       ys2(ivelx) = vin
       ys2(iposx) = zin


! set the step size and minimum step size to try.
      if (one_step) then
       stptry = time_end
       stpmin = time_end * 1.0d-20
      else
       stptry = max(time_start * 1.0d-10,1.0d-16)
       stpmin = stptry * 1.0d-12
      end if


! scale the energy, temperature and density equations
      call aprox13(stptry,ys2,zwork1)
      energy_scale = max(ein, clight*clight, abs(sdot-sneut))
      inv_energy_scale = 1.0d0/energy_scale
      ys2(iener) = ein/energy_scale

      temp_scale     = 1.0d10
      inv_temp_scale = 1.0d0/temp_scale
      ys2(itemp)     = tin * inv_temp_scale

      den_scale     = din * 10.0d0
      inv_den_scale = 1.0d0/den_scale
      ys2(iden)     = din * inv_den_scale


      if (expansion) temp_stop = temp_stop * inv_temp_scale



! integrate the aprox13 network
       call netint(time_start,stptry,stpmin,time_end,dtsav,ys2, &
                   tol,neqs,nok,nbad,kount,odescal, &
!                  aprox13,saprox13,baprox13,forder_ma28)
!                  aprox13,saprox13,baprox13,forder_umf)
!                  aprox13,saprox13,baprox13,forder_y12m)
!                  aprox13,daprox13,baprox13,forder_ludcmp)
!                  aprox13,daprox13,baprox13,forder_leqs)
!                  aprox13,daprox13,baprox13,forder_lapack)
!                  aprox13,daprox13,baprox13,forder_gift)
!                  aprox13,saprox13,baprox13,forder_biconj)
!                  aprox13,saprox13,baprox13,rosen_ma28)
!                  aprox13,saprox13,baprox13,rosen_umf)
!                  aprox13,saprox13,baprox13,rosen_y12m)
!                  aprox13,daprox13,baprox13,rosen_ludcmp)
!                  aprox13,daprox13,baprox13,rosen_leqs)
!                  aprox13,daprox13,baprox13,rosen_lapack)
!                  aprox13,daprox13,baprox13,rosen_gift)
!                  aprox13,saprox13,baprox13,rosen_biconj)
                  aprox13,saprox13,baprox13,stifbs_ma28)
!                  aprox13,saprox13,baprox13,stifbs_umf)
!                  aprox13,saprox13,baprox13,stifbs_y12m)
!                  aprox13,daprox13,baprox13,stifbs_ludcmp)
!                  aprox13,daprox13,baprox13,stifbs_leqs)
!                  aprox13,daprox13,baprox13,stifbs_lapack)
!                  aprox13,daprox13,baprox13,stifbs_gift)
!                  aprox13,saprox13,baprox13,stifbs_biconj)



! set the output composition
! positive definite molar abundances
      do i=1,ionmax
       ys2(i) = min(1.0d0,max(ys2(i),1.0d-30))
      enddo
      xout(1:ionmax) = ys2(1:ionmax) * aion(1:ionmax)
      tout = ys2(itemp)
      dout = ys2(iden)
!      eout = ys2(iener)
      eout = ys2(iener) * energy_scale

! set the mass non-conservation
      conserv = 0.0d0
      do i=ionbeg,ionend
       conserv = conserv + xout(i)
      enddo
      conserv = 1.0d0 - conserv

      return
      end



!---------------------------------------------------------------------







!---------------------------------------------------------------------
! this file contains aprox13 network

! routine aprox13 sets up the odes
! routine rhs evaluates the right hand sides
! routine daprox13 sets up the dense aprox13 jacobian
! routine baprox13 builds the nonzero locations for saprox13
! routine saprox13 sets up the sparse aprox13 jacobian
! routine aprox13rat generates the reaction rates for routine aprox13
! routine aprox13tab generates the raw rates using table interpolation
! routine screen_aprox13 applies screening corrections to the raw rates
! routine init_aprox13 initializes the aprox13 network




      subroutine aprox13(tt,y,dydt)
      include 'implno.dek'
      include 'const.dek'
      include 'burn_common.dek'
      include 'network.dek'
      include 'vector_eos.dek'

! this routine sets up the system of ode's for the aprox13 nuclear reactions.
! this is an alpha chain + heavy ion network with (a,p)(p,g) links
!
! isotopes: he4,  c12,  o16,  ne20, mg24, si28, s32,
!           ar36, ca40, ti44, cr48, fe52, ni56


! declare the pass
      real*8            ::  tt,y(*),dydt(*)


! local variables
      logical,parameter :: deriva = .false.
      integer           :: i
      real*8            :: enuc,taud,taut,abar,zbar,ye,velx,posx,snuda,snudz


! prepare
      call prepare_for_net_eval(tt,y,abar,zbar,ye)


! get the reaction rates
       if (use_tables .eq. 1) then
        call aprox13tab
       else
        call aprox13rat
       end if


! do the screening here because the corrections depend on the composition

       call screen_aprox13(y)


! get the right hand side of the odes

       call rhs(y,ratdum,dydt,deriva)


! instantaneous energy generation rate

       call ener_gener_rate(dydt,enuc)
       sdot = enuc


! get the neutrino losses
       call sneut5(btemp,bden,abar,zbar, &
                   sneut,dsneutdt,dsneutdd,snuda,snudz)



! if we are doing a pure network, we are done

      if (pure_network .eq. 1) return



! append an energy equation
      dydt(iener) = (enuc - sneut) * inv_energy_scale


! the type of temperature and density ode's depend
! on the burn mode:


! hydrostatic or single step cases
      if (hydrostatic  .or.  one_step  .or. trho_hist) then
       dydt(itemp) = 0.0d0
       dydt(iden)  = 0.0d0
       dydt(ivelx) = 0.0d0
       dydt(iposx) = 0.0d0



! adiabatic expansion or contraction
      else if (expansion) then
       taud = 446.0d0/sqrt(den0)
       taut = 3.0d0 * taud
       dydt(itemp) = -psi * y(itemp)/taut 
       dydt(iden)  = -psi * y(iden)/taud  
       dydt(ivelx) = 0.0d0
       dydt(iposx) = 0.0d0


! self heating
      else if (self_heat_const_den) then

! call an eos
       temp_row(1) = y(itemp)
       den_row(1)  = y(iden)
       abar_row(1) = abar
       zbar_row(1) = zbar
       jlo_eos = 1
       jhi_eos = 1

!       call helmeos
       call helmeosaprox13

! temperature equation that is self-consistent with an eos
       dydt(itemp) = dydt(iener)*energy_scale/cv_row(1) * inv_temp_scale
       dydt(iden)  = 0.0d0
       dydt(ivelx) = 0.0d0
       dydt(iposx) = 0.0d0



! self heating at constant pressure
      else if (self_heat_const_pres) then
       dydt(itemp) = dydt(iener)*energy_scale/cp_row(1) * inv_temp_scale
       dydt(iden)  = 0.0d0
!       dydt(iden) = -dpt_row(1)/dpd_row(1) * dydt(itemp) * inv_den_scale
       dydt(ivelx) = 0.0d0
       dydt(iposx) = 0.0d0


! end of burning modes
      end if


      return
      end subroutine aprox13





      subroutine rhs(y,rate,dydt,deriva)
      include 'implno.dek'
      include 'network.dek'

! evaluates the right hand side of the aprox13 odes

! declare the pass
      logical   ::  deriva
      real*8    ::  y(*),rate(*),dydt(*)


! local variables
      integer           :: i
      real*8, parameter :: sixth = 1.0d0/6.0d0

! quad precision dydt sums
      real*16           :: a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17
      real*16           :: qray(abignet)



! zero the odes
      dydt(1:neqs) = 0.0d0
      qray(1:neqs) = 0.0q0


! he4 reactions
! heavy ion reactions
      a1  = 0.5d0 * y(ic12) * y(ic12) * rate(ir1212) 
      a2  = 0.5d0 * y(ic12) * y(io16) * rate(ir1216) 
      a3  = 0.56d0 * 0.5d0 * y(io16) * y(io16) * rate(ir1616)

      qray(ihe4) =  qray(ihe4) + a1 + a2 + a3

! (a,g) and (g,a) reactions
      a1  = -0.5d0 * y(ihe4) * y(ihe4) * y(ihe4) * rate(ir3a) 
      a2  =  3.0d0 * y(ic12) * rate(irg3a) 
      a3  = -y(ihe4)  * y(ic12) * rate(ircag) 
      a4  =  y(io16)  * rate(iroga) 
      a5  = -y(ihe4)  * y(io16) * rate(iroag) 
      a6  =  y(ine20) * rate(irnega) 
      a7  = -y(ihe4)  * y(ine20) * rate(irneag) 
      a8  =  y(img24) * rate(irmgga) 
      a9  = -y(ihe4)  * y(img24)* rate(irmgag) 
      a10 =  y(isi28) * rate(irsiga) 
      a11 = -y(ihe4)  * y(isi28)*rate(irsiag) 
      a12 =  y(is32)  * rate(irsga)

      qray(ihe4) =  qray(ihe4) + a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10 + a11 + a12

      a1  = -y(ihe4)  * y(is32) * rate(irsag) 
      a2  =  y(iar36) * rate(irarga) 
      a3  = -y(ihe4)  * y(iar36)*rate(irarag) 
      a4  =  y(ica40) * rate(ircaga) 
      a5  = -y(ihe4)  * y(ica40)*rate(ircaag) 
      a6  =  y(iti44) * rate(irtiga) 
      a7  = -y(ihe4)  * y(iti44)*rate(irtiag) 
      a8  =  y(icr48) * rate(ircrga) 
      a9  = -y(ihe4)  * y(icr48)*rate(ircrag) 
      a10 =  y(ife52) * rate(irfega) 
      a11 = -y(ihe4)  * y(ife52) * rate(irfeag) 
      a12 =  y(ini56) * rate(irniga)

      qray(ihe4) =  qray(ihe4) + a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10 + a11 + a12

! (a,p)(p,g) and (g,p)(p,a) reactions

      if (.not.deriva) then

       a1  =  0.34d0*0.5d0*y(io16)*y(io16)*rate(irs1)*rate(ir1616) 
       a2  = -y(ihe4)  * y(img24) * rate(irmgap)*(1.0d0-rate(irr1))
       a3  =  y(isi28) * rate(irsigp) * rate(irr1) 
       a4  = -y(ihe4)  * y(isi28) * rate(irsiap)*(1.0d0-rate(irs1)) 
       a5  =  y(is32)  * rate(irsgp) * rate(irs1) 
       a6  = -y(ihe4)  * y(is32) * rate(irsap)*(1.0d0-rate(irt1)) 
       a7  =  y(iar36) * rate(irargp) * rate(irt1) 
       a8  = -y(ihe4)  * y(iar36) * rate(irarap)*(1.0d0-rate(iru1)) 
       a9  =  y(ica40) * rate(ircagp) * rate(iru1) 
       a10 = -y(ihe4)  * y(ica40) * rate(ircaap)*(1.0d0-rate(irv1)) 
       a11 =  y(iti44) * rate(irtigp) * rate(irv1)
       a12 = -y(ihe4)  * y(iti44) * rate(irtiap)*(1.0d0-rate(irw1)) 
       a13 =  y(icr48) * rate(ircrgp) * rate(irw1) 
       a14 = -y(ihe4)  * y(icr48) * rate(ircrap)*(1.0d0-rate(irx1)) 
       a15 =  y(ife52) * rate(irfegp) * rate(irx1) 
       a16 = -y(ihe4)  * y(ife52) * rate(irfeap)*(1.0d0-rate(iry1))
       a17 =  y(ini56) * rate(irnigp) * rate(iry1)

       qray(ihe4) =  qray(ihe4) + a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10 & 
                                + a11 + a12 + a13 + a14 + a15 + a16 + a17

      else
       a1  =  0.34d0*0.5d0*y(io16)*y(io16) * ratdum(irs1) * rate(ir1616) 
       a2  =  0.34d0*0.5d0*y(io16)*y(io16) * rate(irs1) * ratdum(ir1616)
       a3  = -y(ihe4)*y(img24) * rate(irmgap)*(1.0d0 - ratdum(irr1)) 
       a4  =  y(ihe4)*y(img24) * ratdum(irmgap)*rate(irr1)
       a5  =  y(isi28) * ratdum(irsigp) * rate(irr1) 
       a6  =  y(isi28) * rate(irsigp) * ratdum(irr1)
       a7  = -y(ihe4)*y(isi28) * rate(irsiap)*(1.0d0 - ratdum(irs1)) 
       a8  = -y(ihe4)*y(isi28) * ratdum(irsiap) * rate(irs1)
       a9  =  y(is32)  * ratdum(irsgp) * rate(irs1)
       a10 =  y(is32)  * rate(irsgp) * ratdum(irs1)

       qray(ihe4) =  qray(ihe4) + a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10

       a1  = -y(ihe4)*y(is32) * rate(irsap)*(1.0d0 - ratdum(irt1)) 
       a2  =  y(ihe4)*y(is32) * ratdum(irsap)*rate(irt1)
       a3  =  y(iar36) * ratdum(irargp) * rate(irt1) 
       a4  =  y(iar36) * rate(irargp) * ratdum(irt1)
       a5  = -y(ihe4)*y(iar36) * rate(irarap)*(1.0d0 - ratdum(iru1))
       a6  =  y(ihe4)*y(iar36) * ratdum(irarap)*rate(iru1)
       a7  =  y(ica40) * ratdum(ircagp) * rate(iru1)
       a8  =  y(ica40) * rate(ircagp) * ratdum(iru1)
       a9  = -y(ihe4)*y(ica40) * rate(ircaap)*(1.0d0-ratdum (irv1)) 
       a10 =  y(ihe4)*y(ica40) * ratdum(ircaap)*rate(irv1)
       a11 =  y(iti44) * ratdum(irtigp) * rate(irv1)
       a12 =  y(iti44) * rate(irtigp) * ratdum(irv1)

       qray(ihe4) =  qray(ihe4) + a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10 + a11 + a12

       a1  = -y(ihe4)*y(iti44) * rate(irtiap)*(1.0d0 - ratdum(irw1))
       a2  =  y(ihe4)*y(iti44) * ratdum(irtiap)*rate(irw1)
       a3  =  y(icr48) * ratdum(ircrgp) * rate(irw1) 
       a4  =  y(icr48) * rate(ircrgp) * ratdum(irw1)
       a5  = -y(ihe4)*y(icr48) * rate(ircrap)*(1.0d0 - ratdum(irx1))
       a6  =  y(ihe4)*y(icr48) * ratdum(ircrap)*rate(irx1)
       a7  =  y(ife52) * ratdum(irfegp) * rate(irx1) 
       a8  =  y(ife52) * rate(irfegp) * ratdum(irx1)
       a9  = -y(ihe4)*y(ife52) * rate(irfeap)*(1.0d0 - ratdum(iry1))
       a10 =  y(ihe4)*y(ife52) * ratdum(irfeap)*rate(iry1)
       a11 =  y(ini56) * ratdum(irnigp) * rate(iry1) 
       a12 =  y(ini56) * rate(irnigp) * ratdum(iry1)

       qray(ihe4) =  qray(ihe4) + a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10 + a11 + a12
      end if


! c12 reactions
      a1 = -y(ic12) * y(ic12) * rate(ir1212) 
      a2 = -y(ic12) * y(io16) * rate(ir1216) 
      a3 =  sixth * y(ihe4) * y(ihe4) * y(ihe4) * rate(ir3a) 
      a4 = -y(ic12) * rate(irg3a) 
      a5 = -y(ic12) * y(ihe4) * rate(ircag) 
      a6 =  y(io16) * rate(iroga) 

      qray(ic12) =  qray(ic12) + a1 + a2 + a3 + a4 + a5 + a6


! o16 reactions
      a1 = -y(ic12) * y(io16) * rate(ir1216) 
      a2 = -y(io16) * y(io16) * rate(ir1616) 
      a3 =  y(ic12) * y(ihe4) * rate(ircag) 
      a4 = -y(io16) * y(ihe4) * rate(iroag) 
      a5 = -y(io16) * rate(iroga) 
      a6 =  y(ine20) * rate(irnega) 

      qray(io16) =  qray(io16) + a1 + a2 + a3 + a4 + a5 + a6 


! ne20 reactions
      a1 =  0.5d0 * y(ic12) * y(ic12) * rate(ir1212) 
      a2 =  y(io16) * y(ihe4) * rate(iroag) 
      a3 = -y(ine20) * y(ihe4) * rate(irneag) 
      a4 = -y(ine20) * rate(irnega) 
      a5 =  y(img24) * rate(irmgga) 

      qray(ine20) =  qray(ine20) + a1 + a2 + a3 + a4 + a5


! mg24 reactions
      a1 =  0.5d0 * y(ic12) * y(io16) * rate(ir1216) 
      a2 =  y(ine20) * y(ihe4) * rate(irneag) 
      a3 = -y(img24) * y(ihe4) * rate(irmgag) 
      a4 = -y(img24) * rate(irmgga) 
      a5 =  y(isi28) * rate(irsiga)

      qray(img24) =  qray(img24) + a1 + a2 + a3 + a4 + a5 

      if (.not.deriva) then
       a1 = -y(img24) * y(ihe4) * rate(irmgap)*(1.0d0-rate(irr1)) 
       a2 =  y(isi28) * rate(irr1) * rate(irsigp)

       qray(img24) =  qray(img24) + a1 + a2

      else
       a1 = -y(img24)*y(ihe4) * rate(irmgap)*(1.0d0 - ratdum(irr1))
       a2 =  y(img24)*y(ihe4) * ratdum(irmgap)*rate(irr1)
       a3 =  y(isi28) * ratdum(irr1) * rate(irsigp) 
       a4 =  y(isi28) * rate(irr1) * ratdum(irsigp)

       qray(img24) =  qray(img24) + a1 + a2 + a3 + a4
      end if



! si28 reactions
      a1 =  0.5d0 * y(ic12) * y(io16) * rate(ir1216) 
      a2 =  0.56d0 * 0.5d0*y(io16) * y(io16) * rate(ir1616) 
      a3 =  y(img24) * y(ihe4) * rate(irmgag) 
      a4 = -y(isi28) * y(ihe4) * rate(irsiag) 
      a5 = -y(isi28) * rate(irsiga) 
      a6 =  y(is32)  * rate(irsga)

      qray(isi28) =  qray(isi28) + a1 + a2 + a3 + a4 + a5 + a6

      if (.not.deriva) then

       a1 =  0.34d0*0.5d0*y(io16)*y(io16)*rate(irs1)*rate(ir1616) 
       a2 =  y(img24) * y(ihe4) * rate(irmgap)*(1.0d0-rate(irr1)) 
       a3 = -y(isi28) * rate(irr1) * rate(irsigp) 
       a4 = -y(isi28) * y(ihe4) * rate(irsiap)*(1.0d0-rate(irs1)) 
       a5 =  y(is32)  * rate(irs1) * rate(irsgp)

       qray(isi28) =  qray(isi28) + a1 + a2 + a3 + a4 + a5

      else
       a1  =  0.34d0*0.5d0*y(io16)*y(io16) * ratdum(irs1)*rate(ir1616) 
       a2  =  0.34d0*0.5d0*y(io16)*y(io16) * rate(irs1)*ratdum(ir1616)
       a3  =  y(img24)*y(ihe4) * rate(irmgap)*(1.0d0 - ratdum(irr1))
       a4  = -y(img24)*y(ihe4) * ratdum(irmgap)*rate(irr1)
       a5  = -y(isi28) * ratdum(irr1) * rate(irsigp) 
       a6  = -y(isi28) * rate(irr1) * ratdum(irsigp)
       a7  = -y(isi28)*y(ihe4) * rate(irsiap)*(1.0d0 - ratdum(irs1))
       a8  =  y(isi28)*y(ihe4) * ratdum(irsiap)*rate(irs1)
       a9  = y(is32) * ratdum(irs1) * rate(irsgp)
       a10 = y(is32) * rate(irs1) * ratdum(irsgp)

       qray(isi28) =  qray(isi28) + a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10
      end if



! s32 reactions
      a1 =  0.1d0 * 0.5d0*y(io16) * y(io16) * rate(ir1616) 
      a2 =  y(isi28) * y(ihe4) * rate(irsiag) 
      a3 = -y(is32) * y(ihe4) * rate(irsag) 
      a4 = -y(is32) * rate(irsga) 
      a5 =  y(iar36) * rate(irarga)

      qray(is32) =  qray(is32) + a1 + a2 + a3 + a4 + a5


      if (.not.deriva) then
       a1 =  0.34d0*0.5d0*y(io16)*y(io16)* rate(ir1616)*(1.0d0-rate(irs1)) 
       a2 =  y(isi28) * y(ihe4) * rate(irsiap)*(1.0d0-rate(irs1)) 
       a3 = -y(is32) * rate(irs1) * rate(irsgp) 
       a4 = -y(is32) * y(ihe4) * rate(irsap)*(1.0d0-rate(irt1)) 
       a5 =  y(iar36) * rate(irt1) * rate(irargp)

       qray(is32) =  qray(is32) + a1 + a2 + a3 + a4 + a5

      else
       a1  =  0.34d0*0.5d0*y(io16)*y(io16) * rate(ir1616)*(1.0d0-ratdum(irs1))
       a2  = -0.34d0*0.5d0*y(io16)*y(io16) * ratdum(ir1616)*rate(irs1)
       a3  =  y(isi28)*y(ihe4) * rate(irsiap)*(1.0d0-ratdum(irs1))
       a4  = -y(isi28)*y(ihe4) * ratdum(irsiap)*rate(irs1) 
       a5  = -y(is32) * ratdum(irs1) * rate(irsgp) 
       a6  = -y(is32) * rate(irs1) * ratdum(irsgp)
       a7  = -y(is32)*y(ihe4) * rate(irsap)*(1.0d0-ratdum(irt1))
       a8  =  y(is32)*y(ihe4) * ratdum(irsap)*rate(irt1)
       a9  =  y(iar36) * ratdum(irt1) * rate(irargp)
       a10 =  y(iar36) * rate(irt1) * ratdum(irargp)

       qray(is32) =  qray(is32) + a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10
      end if


! ar36 reactions
      a1 =  y(is32)  * y(ihe4) * rate(irsag) 
      a2 = -y(iar36) * y(ihe4) * rate(irarag)
      a3 = -y(iar36) * rate(irarga) 
      a4 =  y(ica40) * rate(ircaga)

      qray(iar36) =  qray(iar36) + a1 + a2 + a3 + a4

      if (.not.deriva) then
       a1 = y(is32)  * y(ihe4) * rate(irsap)*(1.0d0-rate(irt1)) 
       a2 = -y(iar36) * rate(irt1) * rate(irargp) 
       a3 = -y(iar36) * y(ihe4) * rate(irarap)*(1.0d0-rate(iru1)) 
       a4 =  y(ica40) * rate(ircagp) * rate(iru1)

       qray(iar36) =  qray(iar36) + a1 + a2 + a3 + a4

      else
       a1 =  y(is32)*y(ihe4) * rate(irsap)*(1.0d0 - ratdum(irt1)) 
       a2 = -y(is32)*y(ihe4) * ratdum(irsap)*rate(irt1)
       a3 = -y(iar36) * ratdum(irt1) * rate(irargp) 
       a4 = -y(iar36) * rate(irt1) * ratdum(irargp)
       a5 = -y(iar36)*y(ihe4) * rate(irarap)*(1.0d0-ratdum(iru1))
       a6 =  y(iar36)*y(ihe4) * ratdum(irarap)*rate(iru1)
       a7 =  y(ica40) * ratdum(ircagp) * rate(iru1) 
       a8 =  y(ica40) * rate(ircagp) * ratdum(iru1)

       qray(iar36) =  qray(iar36) + a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 
      end if


! ca40 reactions
      a1 =  y(iar36) * y(ihe4) * rate(irarag)
      a2 = -y(ica40) * y(ihe4) * rate(ircaag)
      a3 = -y(ica40) * rate(ircaga)
      a4 =  y(iti44) * rate(irtiga)

      qray(ica40) =  qray(ica40) + a1 + a2 + a3 + a4

      if (.not.deriva) then
       a1 =  y(iar36) * y(ihe4) * rate(irarap)*(1.0d0-rate(iru1)) 
       a2 = -y(ica40) * rate(ircagp) * rate(iru1) 
       a3 = -y(ica40) * y(ihe4) * rate(ircaap)*(1.0d0-rate(irv1)) 
       a4 =  y(iti44) * rate(irtigp) * rate(irv1)

      qray(ica40) =  qray(ica40) + a1 + a2 + a3 + a4

      else
       a1 =  y(iar36)*y(ihe4) * rate(irarap)*(1.0d0-ratdum(iru1))
       a2 = -y(iar36)*y(ihe4) * ratdum(irarap)*rate(iru1)
       a3 = -y(ica40) * ratdum(ircagp) * rate(iru1) 
       a4 = -y(ica40) * rate(ircagp) * ratdum(iru1)
       a5 = -y(ica40)*y(ihe4) * rate(ircaap)*(1.0d0-ratdum(irv1)) 
       a6 =  y(ica40)*y(ihe4) * ratdum(ircaap)*rate(irv1)
       a7 =  y(iti44) * ratdum(irtigp) * rate(irv1) 
       a8 =  y(iti44) * rate(irtigp) * ratdum(irv1)

       qray(ica40) =  qray(ica40) + a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 
      end if


! ti44 reactions
      a1 =  y(ica40) * y(ihe4) * rate(ircaag) 
      a2 = -y(iti44) * y(ihe4) * rate(irtiag) 
      a3 = -y(iti44) * rate(irtiga) 
      a4 =  y(icr48) * rate(ircrga)

      qray(iti44) =  qray(iti44) + a1 + a2 + a3 + a4

      if (.not.deriva) then
       a1 =  y(ica40) * y(ihe4) * rate(ircaap)*(1.0d0-rate(irv1)) 
       a2 = -y(iti44) * rate(irv1) * rate(irtigp) 
       a3 = -y(iti44) * y(ihe4) * rate(irtiap)*(1.0d0-rate(irw1)) 
       a4 =  y(icr48) * rate(irw1) * rate(ircrgp)

      qray(iti44) =  qray(iti44) + a1 + a2 + a3 + a4

      else
       a1 =  y(ica40)*y(ihe4) * rate(ircaap)*(1.0d0-ratdum(irv1)) 
       a2 = -y(ica40)*y(ihe4) * ratdum(ircaap)*rate(irv1)
       a3 = -y(iti44) * ratdum(irv1) * rate(irtigp) 
       a4 = -y(iti44) * rate(irv1) * ratdum(irtigp)
       a5 = -y(iti44)*y(ihe4) * rate(irtiap)*(1.0d0-ratdum(irw1))
       a6 =  y(iti44)*y(ihe4) * ratdum(irtiap)*rate(irw1)
       a7 =  y(icr48) * ratdum(irw1) * rate(ircrgp)
       a8 =  y(icr48) * rate(irw1) * ratdum(ircrgp)

       qray(iti44) =  qray(iti44) + a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 
      end if


! cr48 reactions
      a1 =  y(iti44) * y(ihe4) * rate(irtiag) 
      a2 = -y(icr48) * y(ihe4) * rate(ircrag) 
      a3 = -y(icr48) * rate(ircrga) 
      a4 =  y(ife52) * rate(irfega)

      qray(icr48) =  qray(icr48) + a1 + a2 + a3 + a4

      if (.not.deriva) then
       a1 =  y(iti44) * y(ihe4) * rate(irtiap)*(1.0d0-rate(irw1)) 
       a2 = -y(icr48) * rate(irw1) * rate(ircrgp) 
       a3 = -y(icr48) * y(ihe4) * rate(ircrap)*(1.0d0-rate(irx1)) 
       a4 =  y(ife52) * rate(irx1) * rate(irfegp)

       qray(icr48) =  qray(icr48) + a1 + a2 + a3 + a4

      else
        a1 =  y(iti44)*y(ihe4) * rate(irtiap)*(1.0d0-ratdum(irw1)) 
        a2 = -y(iti44)*y(ihe4) * ratdum(irtiap)*rate(irw1)
        a3 = -y(icr48) * ratdum(irw1) * rate(ircrgp) 
        a4 = -y(icr48) * rate(irw1) * ratdum(ircrgp)
        a5 = -y(icr48)*y(ihe4) * rate(ircrap)*(1.0d0-ratdum(irx1))
        a6 =  y(icr48)*y(ihe4) * ratdum(ircrap)*rate(irx1)
        a7 =  y(ife52) * ratdum(irx1) * rate(irfegp) 
        a8 =  y(ife52) * rate(irx1) * ratdum(irfegp)

       qray(icr48) =  qray(icr48) + a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 
      end if


! fe52 reactions
      a1 =  y(icr48) * y(ihe4) * rate(ircrag) 
      a2 = -y(ife52) * y(ihe4) * rate(irfeag) 
      a3 = -y(ife52) * rate(irfega) 
      a4 =  y(ini56) * rate(irniga)

      qray(ife52) =  qray(ife52) + a1 + a2 + a3 + a4

      if (.not.deriva) then
       a1 =  y(icr48) * y(ihe4) * rate(ircrap)*(1.0d0-rate(irx1)) 
       a2 = -y(ife52) * rate(irx1) * rate(irfegp) 
       a3 = -y(ife52) * y(ihe4) * rate(irfeap)*(1.0d0-rate(iry1)) 
       a4 =  y(ini56) * rate(iry1) * rate(irnigp)

       qray(ife52) =  qray(ife52) + a1 + a2 + a3 + a4

      else
       a1 =  y(icr48)*y(ihe4) * rate(ircrap)*(1.0d0-ratdum(irx1))
       a2 = -y(icr48)*y(ihe4) * ratdum(ircrap)*rate(irx1)
       a3 = -y(ife52) * ratdum(irx1) * rate(irfegp) 
       a4 = -y(ife52) * rate(irx1) * ratdum(irfegp)
       a5 = -y(ife52)*y(ihe4) * rate(irfeap)*(1.0d0-ratdum(iry1))
       a6 =  y(ife52)*y(ihe4) * ratdum(irfeap)*rate(iry1)
       a7 =  y(ini56) * ratdum(iry1) * rate(irnigp)
       a8 =  y(ini56) * rate(iry1) * ratdum(irnigp)

       qray(ife52) =  qray(ife52) + a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 
      end if


! ni56 reactions
      a1 =  y(ife52) * y(ihe4) * rate(irfeag) 
      a2 = -y(ini56) * rate(irniga) 

      qray(ini56) =  qray(ini56) + a1 + a2 

      if (.not.deriva) then
       a1 =  y(ife52) * y(ihe4) * rate(irfeap)*(1.0d0-rate(iry1)) 
       a2 = -y(ini56) * rate(iry1) * rate(irnigp)

       qray(ini56) =  qray(ini56) + a1 + a2 

      else
       a1 =  y(ife52)*y(ihe4) * rate(irfeap)*(1.0d0-ratdum(iry1))
       a2 = -y(ife52)*y(ihe4) * ratdum(irfeap)*rate(iry1)
       a3 = -y(ini56) * ratdum(iry1) * rate(irnigp)
       a4 = -y(ini56) * rate(iry1) * ratdum(irnigp)

       qray(ini56) =  qray(ini56) + a1 + a2 + a3 + a4
      end if


! now set the double precision return argument dydt

      dydt(1:ionmax) = qray(1:ionmax)


      return
      end subroutine rhs







       subroutine daprox13(tt,y,dfdy,nlog,nphys)
      include 'implno.dek'
      include 'const.dek'
      include 'burn_common.dek'
      include 'network.dek'
      include 'vector_eos.dek'

! this routine sets up the dense aprox13 jacobian


! declare the pass
      integer   :: nlog,nphys
      real*8    :: tt,y(*),dfdy(nphys,nphys)

! local variables
      logical, parameter :: deriva = .true.
      integer            :: i,j
      real*8             :: abar,zbar,ye,taud,taut, b1,snuda,snudz,enuc,velx,posx,zz
      real*8, parameter  :: enuc_conv2 = -avo*clight*clight


! prepare
      call prepare_for_net_eval(tt,y,abar,zbar,ye)
       dfdy(1:nlog,1:nlog) = 0.0d0


! get the reaction rates
       if (use_tables .eq. 1) then
        call aprox13tab
       else
        call aprox13rat
       end if

! do the screening here because the corrections depend on the composition

       call screen_aprox13(y)


! do the jacobian elements for the isotopes
      
      call dfdy_isotopes_aprox13(tt,y,dfdy,nlog,nphys)



! if we are doing a pure network, we are done

      if (pure_network .eq. 1) return


! temperature dependence of the rate equations

      call rhs(y,dratdumdt,zwork1,deriva)
      dfdy(1:ionmax,itemp) = zwork1(1:ionmax)  * inv_temp_scale


! density dependence of the rate equations

      call rhs(y,dratdumdd,zwork1,deriva)
      dfdy(1:ionmax,iden) = zwork1(1:ionmax) * inv_den_scale


! append the energy generation rate jacobian elements
      do j=1,ionmax
       do i=1,ionmax
        dfdy(iener,j) = dfdy(iener,j) + dfdy(i,j)*mion(i)
       enddo
       dfdy(iener,j) = dfdy(iener,j) * enuc_conv2  * inv_energy_scale
       dfdy(iener,itemp) = dfdy(iener,itemp) + dfdy(j,itemp)*mion(j) 
       dfdy(iener,iden)  = dfdy(iener,iden) + dfdy(j,iden)*mion(j) 
      enddo
      dfdy(iener,itemp) = dfdy(iener,itemp) * enuc_conv2 * inv_energy_scale * inv_temp_scale
      dfdy(iener,iden)  = dfdy(iener,iden) * enuc_conv2  * inv_energy_scale * inv_den_scale
      dsdotdt = dfdy(iener,itemp) * (energy_scale * temp_scale)
      dsdotdd = dfdy(iener,iden) * (energy_scale * den_scale)


! account for the thermal neutrino losses
      call sneut5(btemp,bden,abar,zbar,sneut,dsneutdt,dsneutdd,snuda,snudz)

      do j=1,ionmax
       b1 = ((aion(j) - abar)*abar*snuda + (zion(j) - zbar)*abar*snudz) * inv_energy_scale
       dfdy(iener,j) = dfdy(iener,j) - b1
      enddo
      dfdy(iener,itemp) = dfdy(iener,itemp) - dsneutdt * inv_energy_scale * inv_temp_scale
      dfdy(iener,iden)  = dfdy(iener,iden)  - dsneutdd * inv_energy_scale * inv_den_scale




! for hydrostatic or one step or trho_hist burns
! all the temperature and density jacobian elements are zero,
! so there is nothing to do since dfdy is initialized to zero.


! adiabatic expansion
      if (expansion) then
       taud = 446.0d0/sqrt(den0)
       taut = 3.0d0 * taud
       dfdy(itemp,itemp) = -psi/taut 
       dfdy(iden,iden)   = -psi/taud 


! for self-heating
      else if (self_heat_const_den) then


! call an eos
       temp_row(1) = y(itemp)
       den_row(1)  = y(iden)
       abar_row(1) = abar
       zbar_row(1) = zbar
       jlo_eos = 1
       jhi_eos = 1

!       call helmeos
       call helmeosaprox13


! d(itemp)/d(yi)
      zz = energy_scale/cv_row(1) 
      dfdy(itemp,1:ionmax) = zz*dfdy(iener,1:ionmax)

! d(itemp)/d(temp)
       dfdy(itemp,itemp) = zz*dfdy(iener,itemp) 

! d(itemp)/d(den)
       dfdy(itemp,iden) = zz*dfdy(iener,iden) 



! for self-heating at constant pressure
      else if (self_heat_const_pres) then


! d(itemp)/d(yi)
      zz = energy_scale/cp_row(1)
      dfdy(itemp,1:ionmax) = zz*dfdy(iener,1:ionmax)

! d(itemp)/d(temp)
       dfdy(itemp,itemp) = zz*dfdy(iener,itemp) 

! d(itemp)/d(den)
       dfdy(itemp,iden) = zz*dfdy(iener,iden) 


! end of the burning mode ifs
      end if

      return
      end subroutine daprox13




      subroutine dfdy_isotopes_aprox13(tt,y,dfdy,nlog,nphys)
      include 'implno.dek'
      include 'const.dek'
      include 'burn_common.dek'
      include 'network.dek'
      include 'vector_eos.dek'

! this routine sets up the dense aprox21 jacobian for the isotopes
!
! declare the pass
      integer   :: nlog,nphys
      real*8    :: tt,y(*),dfdy(nphys,nphys)



! he4 jacobian elements
      dfdy(ihe4,ihe4)  = -1.5d0 * y(ihe4) * y(ihe4) * ratdum(ir3a) &
               - y(ic12)  * ratdum(ircag) &
               - y(io16)  * ratdum(iroag) &
               - y(ine20) * ratdum(irneag) &
               - y(img24) * ratdum(irmgag) &
               - y(isi28) * ratdum(irsiag) &
               - y(is32)  * ratdum(irsag) &
               - y(iar36) * ratdum(irarag) &
               - y(ica40) * ratdum(ircaag) &
               - y(iti44) * ratdum(irtiag) &
               - y(icr48) * ratdum(ircrag) &
               - y(ife52) * ratdum(irfeag)

      dfdy(ihe4,ihe4)  = dfdy(ihe4,ihe4) &
               - y(img24) * ratdum(irmgap) * (1.0d0-ratdum(irr1)) &
               - y(isi28) * ratdum(irsiap) * (1.0d0-ratdum(irs1)) &
               - y(is32) * ratdum(irsap)   * (1.0d0-ratdum(irt1)) &
               - y(iar36) * ratdum(irarap) * (1.0d0-ratdum(iru1)) &
               - y(ica40) * ratdum(ircaap) * (1.0d0-ratdum(irv1)) &
               - y(iti44) * ratdum(irtiap) * (1.0d0-ratdum(irw1)) &
               - y(icr48) * ratdum(ircrap) * (1.0d0-ratdum(irx1)) &
               - y(ife52) * ratdum(irfeap) * (1.0d0-ratdum(iry1))

      dfdy(ihe4,ic12)  = y(ic12) * ratdum(ir1212) &
                  + 0.5d0 * y(io16) * ratdum(ir1216) &
                  + 3.0d0 * ratdum(irg3a) &
                  - y(ihe4) * ratdum(ircag)

      dfdy(ihe4,io16)  = 0.5d0 * y(ic12) * ratdum(ir1216) &
                  + 1.12d0 * 0.5d0 * y(io16) * ratdum(ir1616) &
                  + 0.68d0*ratdum(irs1)*0.5d0*y(io16)*ratdum(ir1616) &
                  + ratdum(iroga) &
                  - y(ihe4) * ratdum(iroag)

      dfdy(ihe4,ine20) =  ratdum(irnega) &
                  - y(ihe4) * ratdum(irneag)

      dfdy(ihe4,img24) =   ratdum(irmgga) &
                  - y(ihe4) * ratdum(irmgag) &
                  - y(ihe4) * ratdum(irmgap) * (1.0d0-ratdum(irr1))

      dfdy(ihe4,isi28) =   ratdum(irsiga) &
                  - y(ihe4) * ratdum(irsiag) &
                  - y(ihe4) * ratdum(irsiap) * (1.0d0-ratdum(irs1)) &
                  + ratdum(irr1) * ratdum(irsigp)

      dfdy(ihe4,is32)  =   ratdum(irsga) &
                  - y(ihe4) * ratdum(irsag) &
                  - y(ihe4) * ratdum(irsap) * (1.0d0-ratdum(irt1)) &
                  + ratdum(irs1) * ratdum(irsgp)

      dfdy(ihe4,iar36) =   ratdum(irarga) &
                  - y(ihe4) * ratdum(irarag) &
                  - y(ihe4) * ratdum(irarap) * (1.0d0-ratdum(iru1)) &
                  + ratdum(irt1) * ratdum(irargp)

      dfdy(ihe4,ica40) =   ratdum(ircaga) &
                  - y(ihe4) * ratdum(ircaag) &
                  - y(ihe4) * ratdum(ircaap) * (1.0d0-ratdum(irv1)) &
                  + ratdum(iru1) * ratdum(ircagp)

      dfdy(ihe4,iti44) =   ratdum(irtiga) &
                  - y(ihe4) * ratdum(irtiag) &
                  - y(ihe4) * ratdum(irtiap) * (1.0d0-ratdum(irw1)) &
                  + ratdum(irv1) * ratdum(irtigp)

      dfdy(ihe4,icr48) =   ratdum(ircrga) &
                  - y(ihe4) * ratdum(ircrag) &
                  - y(ihe4) * ratdum(ircrap) * (1.0d0-ratdum(irx1)) &
                  + ratdum(irw1) * ratdum(ircrgp)

      dfdy(ihe4,ife52) =   ratdum(irfega) &
                  - y(ihe4) * ratdum(irfeag) &
                  - y(ihe4) * ratdum(irfeap) * (1.0d0-ratdum(iry1)) &
                  + ratdum(irx1) * ratdum(irfegp)

      dfdy(ihe4,ini56) =   ratdum(irniga) &
                  + ratdum(iry1) * ratdum(irnigp)


! c12 jacobian elements
      dfdy(ic12,ihe4) = 0.5d0*y(ihe4)*y(ihe4)*ratdum(ir3a) &
                - y(ic12) * ratdum(ircag)

      dfdy(ic12,ic12) = -2.0d0 * y(ic12) * ratdum(ir1212) &
                 - y(io16) * ratdum(ir1216) &
                 - ratdum(irg3a) &
                 - y(ihe4) * ratdum(ircag)

      dfdy(ic12,io16) = -y(ic12) * ratdum(ir1216) &
                + ratdum(iroga)



! o16 jacobian elements
      dfdy(io16,ihe4) = y(ic12)*ratdum(ircag) &
               - y(io16)*ratdum(iroag)

      dfdy(io16,ic12) = -y(io16)*ratdum(ir1216) &
                + y(ihe4)*ratdum(ircag)

      dfdy(io16,io16) = -y(ic12) * ratdum(ir1216) &
                 - 2.0d0 * y(io16) * ratdum(ir1616) &
                 - y(ihe4) * ratdum(iroag) &
                 - ratdum(iroga)

      dfdy(io16,ine20) = ratdum(irnega)



! ne20 jacobian elements
      dfdy(ine20,ihe4)  = y(io16) * ratdum(iroag) &
                 - y(ine20) * ratdum(irneag)

      dfdy(ine20,ic12)  = y(ic12) * ratdum(ir1212)

      dfdy(ine20,io16)  = y(ihe4) * ratdum(iroag)

      dfdy(ine20,ine20) = -y(ihe4) * ratdum(irneag) &
                  - ratdum(irnega)

      dfdy(ine20,img24) = ratdum(irmgga)



! mg24 jacobian elements
      dfdy(img24,ihe4)  = y(ine20) * ratdum(irneag) &
                  -y(img24) * ratdum(irmgag) &
                  -y(img24) * ratdum(irmgap) * (1.0d0-ratdum(irr1))

      dfdy(img24,ic12)  = 0.5d0 * y(io16) * ratdum(ir1216)

      dfdy(img24,io16)  = 0.5d0 * y(ic12) * ratdum(ir1216)

      dfdy(img24,ine20) = y(ihe4) * ratdum(irneag)

      dfdy(img24,img24) = -y(ihe4) * ratdum(irmgag) &
                  - ratdum(irmgga) &
                  - y(ihe4) * ratdum(irmgap) * (1.0d0-ratdum(irr1))

      dfdy(img24,isi28) = ratdum(irsiga) &
                 + ratdum(irr1) * ratdum(irsigp)



! si28 jacobian elements
      dfdy(isi28,ihe4)  = y(img24) * ratdum(irmgag) &
                 - y(isi28) * ratdum(irsiag) &
                 + y(img24) * ratdum(irmgap) * (1.0d0-ratdum(irr1)) &
                 - y(isi28) * ratdum(irsiap) * (1.0d0-ratdum(irs1))

      dfdy(isi28,ic12)  = 0.5d0 * y(io16) * ratdum(ir1216)

      dfdy(isi28,io16)  =   0.5d0 * y(ic12) * ratdum(ir1216) &
                   + 1.12d0 * 0.5d0*y(io16) * ratdum(ir1616) &
                   + 0.68d0*0.5d0*y(io16)*ratdum(irs1)*ratdum(ir1616)

      dfdy(isi28,img24) = y(ihe4) * ratdum(irmgag) &
                 + y(ihe4) * ratdum(irmgap) * (1.0d0-ratdum(irr1))

      dfdy(isi28,isi28) = -y(ihe4) * ratdum(irsiag) &
                  - ratdum(irsiga) &
                  - ratdum(irr1) * ratdum(irsigp) &
                  - y(ihe4) * ratdum(irsiap) * (1.0d0-ratdum(irs1))

      dfdy(isi28,is32)  = ratdum(irsga) &
                 + ratdum(irs1) * ratdum(irsgp)



! s32 jacobian elements
      dfdy(is32,ihe4)  = y(isi28) * ratdum(irsiag) &
                - y(is32) * ratdum(irsag) &
                + y(isi28) * ratdum(irsiap) * (1.0d0-ratdum(irs1)) &
                - y(is32) * ratdum(irsap) * (1.0d0-ratdum(irt1))

      dfdy(is32,io16)  = 0.68d0*0.5d0*y(io16)*ratdum(ir1616)*(1.0d0-ratdum(irs1)) &
                  + 0.2d0 * 0.5d0*y(io16) * ratdum(ir1616)

      dfdy(is32,isi28) = y(ihe4) * ratdum(irsiag) &
                 + y(ihe4) * ratdum(irsiap) * (1.0d0-ratdum(irs1))

      dfdy(is32,is32)  = -y(ihe4) * ratdum(irsag) &
                 - ratdum(irsga) &
                 - ratdum(irs1) * ratdum(irsgp) &
                 - y(ihe4) * ratdum(irsap) * (1.0d0-ratdum(irt1))

      dfdy(is32,iar36) = ratdum(irarga) &
                + ratdum(irt1) * ratdum(irargp)


! ar36 jacobian elements
      dfdy(iar36,ihe4)  = y(is32)  * ratdum(irsag) &
                - y(iar36) * ratdum(irarag) &
                + y(is32)  * ratdum(irsap) * (1.0d0-ratdum(irt1)) &
                - y(iar36) * ratdum(irarap) * (1.0d0-ratdum(iru1))

      dfdy(iar36,is32)  = y(ihe4) * ratdum(irsag) &
                  + y(ihe4) * ratdum(irsap) * (1.0d0-ratdum(irt1))

      dfdy(iar36,iar36) = -y(ihe4) * ratdum(irarag) &
                  - ratdum(irarga) &
                  - ratdum(irt1) * ratdum(irargp) &
                  - y(ihe4) * ratdum(irarap) * (1.0d0-ratdum(iru1))

      dfdy(iar36,ica40) = ratdum(ircaga) &
                        + ratdum(ircagp) * ratdum(iru1)



! ca40 jacobian elements
      dfdy(ica40,ihe4)   = y(iar36) * ratdum(irarag) &
                  - y(ica40) * ratdum(ircaag) &
                  + y(iar36) * ratdum(irarap)*(1.0d0-ratdum(iru1)) &
                  - y(ica40) * ratdum(ircaap)*(1.0d0-ratdum(irv1))

      dfdy(ica40,iar36)  = y(ihe4) * ratdum(irarag) &
                  + y(ihe4) * ratdum(irarap)*(1.0d0-ratdum(iru1))

      dfdy(ica40,ica40)  = -y(ihe4) * ratdum(ircaag) &
                   - ratdum(ircaga) &
                   - ratdum(ircagp) * ratdum(iru1) &
                   - y(ihe4) * ratdum(ircaap)*(1.0d0-ratdum(irv1))

      dfdy(ica40,iti44)  = ratdum(irtiga) &
                          + ratdum(irtigp) * ratdum(irv1)



! ti44 jacobian elements
      dfdy(iti44,ihe4)   = y(ica40) * ratdum(ircaag) &
                  - y(iti44) * ratdum(irtiag) &
                  + y(ica40) * ratdum(ircaap)*(1.0d0-ratdum(irv1)) &
                  - y(iti44) * ratdum(irtiap)*(1.0d0-ratdum(irw1))

      dfdy(iti44,ica40)  = y(ihe4) * ratdum(ircaag) &
                  + y(ihe4) * ratdum(ircaap)*(1.0d0-ratdum(irv1))

      dfdy(iti44,iti44)  = -y(ihe4) * ratdum(irtiag) &
                   - ratdum(irtiga) &
                   - ratdum(irv1) * ratdum(irtigp) &
                   - y(ihe4) * ratdum(irtiap)*(1.0d0-ratdum(irw1))

      dfdy(iti44,icr48)  = ratdum(ircrga) &
                  + ratdum(irw1) * ratdum(ircrgp)



! cr48 jacobian elements
      dfdy(icr48,ihe4)  = y(iti44) * ratdum(irtiag) &
                 - y(icr48) * ratdum(ircrag) &
                 + y(iti44) * ratdum(irtiap)*(1.0d0-ratdum(irw1)) &
                 - y(icr48) * ratdum(ircrap)*(1.0d0-ratdum(irx1))

      dfdy(icr48,iti44) = y(ihe4) * ratdum(irtiag) &
                 + y(ihe4) * ratdum(irtiap)*(1.0d0-ratdum(irw1))

      dfdy(icr48,icr48) = -y(ihe4) * ratdum(ircrag) &
                  - ratdum(ircrga) &
                  - ratdum(irw1) * ratdum(ircrgp) &
                  - y(ihe4) * ratdum(ircrap)*(1.0d0-ratdum(irx1))

      dfdy(icr48,ife52) = ratdum(irfega) &
                 + ratdum(irx1) * ratdum(irfegp)



! fe52 jacobian elements
      dfdy(ife52,ihe4)  = y(icr48) * ratdum(ircrag) &
                 - y(ife52) * ratdum(irfeag) &
                 + y(icr48) * ratdum(ircrap) * (1.0d0-ratdum(irx1)) &
                 - y(ife52) * ratdum(irfeap) * (1.0d0-ratdum(iry1))

      dfdy(ife52,icr48) = y(ihe4) * ratdum(ircrag) &
                 + y(ihe4) * ratdum(ircrap) * (1.0d0-ratdum(irx1))

      dfdy(ife52,ife52) = -y(ihe4) * ratdum(irfeag) &
                   - ratdum(irfega) &
                   - ratdum(irx1) * ratdum(irfegp) &
                   - y(ihe4) * ratdum(irfeap) * (1.0d0-ratdum(iry1))

      dfdy(ife52,ini56) = ratdum(irniga) &
                 + ratdum(iry1) * ratdum(irnigp)



! ni56 jacobian elements
      dfdy(ini56,ihe4)  = y(ife52) * ratdum(irfeag) &
                 + y(ife52) * ratdum(irfeap) * (1.0d0-ratdum(iry1))

      dfdy(ini56,ife52) = y(ihe4) * ratdum(irfeag) &
                 + y(ihe4) * ratdum(irfeap) * (1.0d0-ratdum(iry1))

      dfdy(ini56,ini56) = -ratdum(irniga) &
                  - ratdum(iry1) * ratdum(irnigp)


      return
      end subroutine dfdy_isotopes_aprox13






      subroutine baprox13(iloc,jloc,nzo,np)
      include 'implno.dek'
      include 'network.dek'

! this routine builds the nonzero matrix locations for saprox13
! input is the integer arrys iloc and jloc, both of dimension np, that
! on output contain nzo matrix element locations.
!
! declare the pass
      integer    :: np,iloc(np),jloc(np),nzo


! communicate with saprox13
      integer          neloc
      parameter        (neloc=172*nzmax)
      integer          eloc(neloc),nterms
      common /elca13/  eloc,nterms


! local variables
      integer   ::  i


! initialize
      nterms = 0
      nzo    = 0
      eloc(1:neloc) = 0
      call tree_init(neqs)



! tag the nonzero locations
! he4 jacobian elements
      call tree(ihe4,ihe4,eloc,neloc,nterms,nzo,iloc,jloc,np)
      call tree(ihe4,ic12,eloc,neloc,nterms,nzo,iloc,jloc,np)
      call tree(ihe4,io16,eloc,neloc,nterms,nzo,iloc,jloc,np)
      call tree(ihe4,ine20,eloc,neloc,nterms,nzo,iloc,jloc,np)
      call tree(ihe4,img24,eloc,neloc,nterms,nzo,iloc,jloc,np)
      call tree(ihe4,isi28,eloc,neloc,nterms,nzo,iloc,jloc,np)
      call tree(ihe4,is32,eloc,neloc,nterms,nzo,iloc,jloc,np)
      call tree(ihe4,iar36,eloc,neloc,nterms,nzo,iloc,jloc,np)
      call tree(ihe4,ica40,eloc,neloc,nterms,nzo,iloc,jloc,np)
      call tree(ihe4,iti44,eloc,neloc,nterms,nzo,iloc,jloc,np)
      call tree(ihe4,icr48,eloc,neloc,nterms,nzo,iloc,jloc,np)
      call tree(ihe4,ife52,eloc,neloc,nterms,nzo,iloc,jloc,np)
      call tree(ihe4,ini56,eloc,neloc,nterms,nzo,iloc,jloc,np)

! c12 jacobian elements
      call tree(ic12,ihe4,eloc,neloc,nterms,nzo,iloc,jloc,np)
      call tree(ic12,ic12,eloc,neloc,nterms,nzo,iloc,jloc,np)
      call tree(ic12,io16,eloc,neloc,nterms,nzo,iloc,jloc,np)

! o16 jacobian elements
      call tree(io16,ihe4,eloc,neloc,nterms,nzo,iloc,jloc,np)
      call tree(io16,ic12,eloc,neloc,nterms,nzo,iloc,jloc,np)
      call tree(io16,io16,eloc,neloc,nterms,nzo,iloc,jloc,np)
      call tree(io16,ine20,eloc,neloc,nterms,nzo,iloc,jloc,np)

! ne20 jacobian elements
      call tree(ine20,ihe4,eloc,neloc,nterms,nzo,iloc,jloc,np)
      call tree(ine20,ic12,eloc,neloc,nterms,nzo,iloc,jloc,np)
      call tree(ine20,io16,eloc,neloc,nterms,nzo,iloc,jloc,np)
      call tree(ine20,ine20,eloc,neloc,nterms,nzo,iloc,jloc,np)
      call tree(ine20,img24,eloc,neloc,nterms,nzo,iloc,jloc,np)

! mg24 jacobian elements
      call tree(img24,ihe4,eloc,neloc,nterms,nzo,iloc,jloc,np)
      call tree(img24,ic12,eloc,neloc,nterms,nzo,iloc,jloc,np)
      call tree(img24,io16,eloc,neloc,nterms,nzo,iloc,jloc,np)
      call tree(img24,ine20,eloc,neloc,nterms,nzo,iloc,jloc,np)
      call tree(img24,img24,eloc,neloc,nterms,nzo,iloc,jloc,np)
      call tree(img24,isi28,eloc,neloc,nterms,nzo,iloc,jloc,np)

! si28 jacobian elements
      call tree(isi28,ihe4,eloc,neloc,nterms,nzo,iloc,jloc,np)
      call tree(isi28,ic12,eloc,neloc,nterms,nzo,iloc,jloc,np)
      call tree(isi28,io16,eloc,neloc,nterms,nzo,iloc,jloc,np)
      call tree(isi28,img24,eloc,neloc,nterms,nzo,iloc,jloc,np)
      call tree(isi28,isi28,eloc,neloc,nterms,nzo,iloc,jloc,np)
      call tree(isi28,is32,eloc,neloc,nterms,nzo,iloc,jloc,np)

! s32 jacobian elements
      call tree(is32,ihe4,eloc,neloc,nterms,nzo,iloc,jloc,np)
      call tree(is32,io16,eloc,neloc,nterms,nzo,iloc,jloc,np)
      call tree(is32,isi28,eloc,neloc,nterms,nzo,iloc,jloc,np)
      call tree(is32,is32,eloc,neloc,nterms,nzo,iloc,jloc,np)
      call tree(is32,iar36,eloc,neloc,nterms,nzo,iloc,jloc,np)

! ar36 jacobian elements
      call tree(iar36,ihe4,eloc,neloc,nterms,nzo,iloc,jloc,np)
      call tree(iar36,is32,eloc,neloc,nterms,nzo,iloc,jloc,np)
      call tree(iar36,iar36,eloc,neloc,nterms,nzo,iloc,jloc,np)
      call tree(iar36,ica40,eloc,neloc,nterms,nzo,iloc,jloc,np)

! ca40 jacobian elements
      call tree(ica40,ihe4,eloc,neloc,nterms,nzo,iloc,jloc,np)
      call tree(ica40,iar36,eloc,neloc,nterms,nzo,iloc,jloc,np)
      call tree(ica40,ica40,eloc,neloc,nterms,nzo,iloc,jloc,np)
      call tree(ica40,iti44,eloc,neloc,nterms,nzo,iloc,jloc,np)

! ti44 jacobian elements
      call tree(iti44,ihe4,eloc,neloc,nterms,nzo,iloc,jloc,np)
      call tree(iti44,ica40,eloc,neloc,nterms,nzo,iloc,jloc,np)
      call tree(iti44,iti44,eloc,neloc,nterms,nzo,iloc,jloc,np)
      call tree(iti44,icr48,eloc,neloc,nterms,nzo,iloc,jloc,np)

! cr48 jacobian elements
      call tree(icr48,ihe4,eloc,neloc,nterms,nzo,iloc,jloc,np)
      call tree(icr48,iti44,eloc,neloc,nterms,nzo,iloc,jloc,np)
      call tree(icr48,icr48,eloc,neloc,nterms,nzo,iloc,jloc,np)
      call tree(icr48,ife52,eloc,neloc,nterms,nzo,iloc,jloc,np)

! fe52 jacobian elements
      call tree(ife52,ihe4,eloc,neloc,nterms,nzo,iloc,jloc,np)
      call tree(ife52,icr48,eloc,neloc,nterms,nzo,iloc,jloc,np)
      call tree(ife52,ife52,eloc,neloc,nterms,nzo,iloc,jloc,np)
      call tree(ife52,ini56,eloc,neloc,nterms,nzo,iloc,jloc,np)

! ni56 jacobian elements
      call tree(ini56,ihe4,eloc,neloc,nterms,nzo,iloc,jloc,np)
      call tree(ini56,ife52,eloc,neloc,nterms,nzo,iloc,jloc,np)
      call tree(ini56,ini56,eloc,neloc,nterms,nzo,iloc,jloc,np)


! if we are doing a pure network, we are done

      if (pure_network .eq. 1) return



! temperature contributions
      do i=ionbeg,ionend
       call tree(i,itemp,eloc,neloc,nterms,nzo,iloc,jloc,np)
      end do

! density contributions
      do i=ionbeg,ionend
       call tree(i,iden,eloc,neloc,nterms,nzo,iloc,jloc,np)
      end do


! energy equation jacobian elements
      do i=ionbeg,ionend
       call tree(iener,i,eloc,neloc,nterms,nzo,iloc,jloc,np)
      enddo
      call tree(iener,iener,eloc,neloc,nterms,nzo,iloc,jloc,np)
      call tree(iener,itemp,eloc,neloc,nterms,nzo,iloc,jloc,np)
      call tree(iener,iden,eloc,neloc,nterms,nzo,iloc,jloc,np)
      call tree(iener,ivelx,eloc,neloc,nterms,nzo,iloc,jloc,np)
      call tree(iener,iposx,eloc,neloc,nterms,nzo,iloc,jloc,np)


! neutrino losses
      do i=ionbeg,ionend
       call tree(iener,i,eloc,neloc,nterms,nzo,iloc,jloc,np)
      enddo
      call tree(iener,itemp,eloc,neloc,nterms,nzo,iloc,jloc,np)
      call tree(iener,iden,eloc,neloc,nterms,nzo,iloc,jloc,np)



! the jacobian elements depend on the burning mode

! hydrostatic or single step
      if (hydrostatic .or. one_step .or. trho_hist) then
       call tree(itemp,itemp,eloc,neloc,nterms,nzo,iloc,jloc,np)
       call tree(iden,iden,eloc,neloc,nterms,nzo,iloc,jloc,np)
       call tree(ivelx,ivelx,eloc,neloc,nterms,nzo,iloc,jloc,np)
       call tree(iposx,iposx,eloc,neloc,nterms,nzo,iloc,jloc,np)



! adiabatic expansion
      else if (expansion) then
       call tree(itemp,itemp,eloc,neloc,nterms,nzo,iloc,jloc,np)
       call tree(iden,iden,eloc,neloc,nterms,nzo,iloc,jloc,np)
       call tree(ivelx,ivelx,eloc,neloc,nterms,nzo,iloc,jloc,np)
       call tree(iposx,iposx,eloc,neloc,nterms,nzo,iloc,jloc,np)


! self heating
      else if (self_heat_const_den) then
       do i=ionbeg,ionend
        call tree(itemp,i,eloc,neloc,nterms,nzo,iloc,jloc,np)
       enddo
       call tree(itemp,itemp,eloc,neloc,nterms,nzo,iloc,jloc,np)
       call tree(itemp,iden,eloc,neloc,nterms,nzo,iloc,jloc,np)
       call tree(iden,iden,eloc,neloc,nterms,nzo,iloc,jloc,np)
       call tree(ivelx,ivelx,eloc,neloc,nterms,nzo,iloc,jloc,np)
       call tree(iposx,iposx,eloc,neloc,nterms,nzo,iloc,jloc,np)


! self heating at constant pressure
      else if (self_heat_const_pres) then
       do i=1,ionmax
        call tree(itemp,i,eloc,neloc,nterms,nzo,iloc,jloc,np)
       enddo
       call tree(itemp,itemp,eloc,neloc,nterms,nzo,iloc,jloc,np)
       call tree(itemp,iden,eloc,neloc,nterms,nzo,iloc,jloc,np)
       do i=1,ionmax
        call tree(iden,i,eloc,neloc,nterms,nzo,iloc,jloc,np)
       enddo
       call tree(iden,itemp,eloc,neloc,nterms,nzo,iloc,jloc,np)
       call tree(iden,iden,eloc,neloc,nterms,nzo,iloc,jloc,np)
       call tree(ivelx,ivelx,eloc,neloc,nterms,nzo,iloc,jloc,np)
       call tree(iposx,iposx,eloc,neloc,nterms,nzo,iloc,jloc,np)
      end if



! store the number of non-zero matrix elements in common

      non_zero_elements = nzo


! write a diagnostic
!      write(6,*) ' '
!      write(6,*) nzo,' matrix elements'
!      write(6,*) nterms,' jacobian contributions'
!      write(6,*) ' '

      return
      end subroutine baprox13





      subroutine add_one(i,j,a,dfdy,nt)
      include 'implno.dek'
      include 'network.dek'

! this routine adds one jacaobian element d(y_i)/d(y_j) = a to 
! the sparse matrix array dfdy 
!

! declare the pass
       integer :: i,j,nt
       real*8  :: a
       real*16 :: dfdy(*)


! local variables
       integer :: iat

! communicate with saprox13
      integer          neloc
      parameter        (neloc=172*nzmax)
      integer          eloc(neloc),nterms
      common /elca13/  eloc,nterms

! add d(i)/d(j)
       nt        = nt + 1
       iat       = eloc(nt)
       dfdy(iat) = dfdy(iat) + a
       return
       end subroutine add_one




      subroutine add_two(i,j,a,dfdy,xsum,nt)
      include 'implno.dek'
      include 'network.dek'

! this routine adds one jacaobian element d(y_i)/d(y_j) = a to 
! the sparse matrix array dfdy and accumulates the energy generation
! sum d(yi)/dt * mion(i) in the xsum array.
!

! declare the pass
       integer :: i,j,nt
       real*16 :: a,dfdy(*),xsum(*)


! local variables
       integer :: iat


! communicate with the saprox13
      integer          neloc
      parameter        (neloc=172*nzmax)
      integer          eloc(neloc),nterms
      common /elca13/  eloc,nterms


! add d(i)/d(j) and piece of energy generation sum
       nt        = nt + 1
       iat       = eloc(nt)
       dfdy(iat) = dfdy(iat) + a
       xsum(j)   = xsum(j) + a * mion(i)

       return
       end subroutine add_two







      subroutine saprox13(tt,y,dfdy,num_nonzero)
      include 'implno.dek'
      include 'const.dek'
      include 'burn_common.dek'
      include 'network.dek'
      include 'vector_eos.dek'
      include 'sparse_matrix.dek'

! this routine sets up the sparse aprox13 jacobian.
! input is tt (irrelevant here) and the abundances y(1).
! output is the jacobian dfdy(num_nonzero).

! declare the pass
      integer  :: num_nonzero
      real*8   :: tt,y(*),dfdy(*)


! local variables
      logical, parameter :: deriva = .true.
      integer            :: i,nt,iat
      real*8             :: zz,abar,zbar,z2bar,taud,taut,snuda,snudz,ye,enucbt,enucbd
      real*8             :: xsum(abignet),x1,x2,x3,x4,x5,x6
      real*16            :: b(30),a1,a2,a4

! conversion factors for nuclear energy generation rate
      real*8, parameter :: enuc_conv2 = -avo*clight*clight


! communicate with the jacobian builder
      integer          neloc
      parameter        (neloc=172*nzmax)
      integer          eloc(neloc),nterms
      common /elca13/  eloc,nterms



! prepare
      call prepare_for_net_eval(tt,y,abar,zbar,ye)
      nt = 0
      qdfdy(1:num_nonzero) = 0.0q0
      qsum(1:ionmax)       = 0.0q0

! get the reaction rates
       if (use_tables .eq. 1) then
        call aprox13tab
       else
        call aprox13rat
       end if



! do the screening here because the corrections depend on the composition

      call screen_aprox13(y)



! he4 jacobian elements
! d(he4)/d(he4)
      b(1)  = -1.5d0 * y(ihe4) * y(ihe4) * ratdum(ir3a) 
      b(2)  = -y(ic12)  * ratdum(ircag) 
      b(3)  = -y(io16)  * ratdum(iroag) 
      b(4)  = -y(ine20) * ratdum(irneag) 
      b(5)  = -y(img24) * ratdum(irmgag) 
      b(6)  = -y(isi28) * ratdum(irsiag) 
      b(7)  = -y(is32)  * ratdum(irsag) 
      b(8)  = -y(iar36) * ratdum(irarag) 
      b(9)  = -y(ica40) * ratdum(ircaag) 
      b(10) = -y(iti44) * ratdum(irtiag) 
      b(11) = -y(icr48) * ratdum(ircrag) 
      b(12) = -y(ife52) * ratdum(irfeag)
      b(13) = -y(img24) * ratdum(irmgap) * (1.0d0-ratdum(irr1)) 
      b(14) = -y(isi28) * ratdum(irsiap) * (1.0d0-ratdum(irs1)) 
      b(15) = -y(is32)  * ratdum(irsap)  * (1.0d0-ratdum(irt1)) 
      b(16) = -y(iar36) * ratdum(irarap) * (1.0d0-ratdum(iru1)) 
      b(17) = -y(ica40) * ratdum(ircaap) * (1.0d0-ratdum(irv1)) 
      b(18) = -y(iti44) * ratdum(irtiap) * (1.0d0-ratdum(irw1)) 
      b(19) = -y(icr48) * ratdum(ircrap) * (1.0d0-ratdum(irx1)) 
      b(20) = -y(ife52) * ratdum(irfeap) * (1.0d0-ratdum(iry1))
      a1   = sum(b(1:20))
      call add_two(ihe4,ihe4,a1,qdfdy,qsum,nt)


! d(he4)/d(c12)
      b(1) =  y(ic12) * ratdum(ir1212) 
      b(2) =  0.5d0 * y(io16) * ratdum(ir1216) 
      b(3) =  3.0d0 * ratdum(irg3a) 
      b(4) = -y(ihe4) * ratdum(ircag)
      a1   = sum(b(1:4))
      call add_two(ihe4,ic12,a1,qdfdy,qsum,nt)

! d(he4)/d(o16)
      b(1) =  0.5d0 * y(ic12) * ratdum(ir1216) 
      b(2) =  1.12d0 * 0.5d0*y(io16) * ratdum(ir1616) 
      b(3) =  0.68d0 * ratdum(irs1) * 0.5d0*y(io16) * ratdum(ir1616) 
      b(4) =  ratdum(iroga) 
      b(5) = -y(ihe4) * ratdum(iroag) 
      a1   = sum(b(1:5))
      call add_two(ihe4,io16,a1,qdfdy,qsum,nt)

! d(he4)/d(ne20)
      b(1) =  ratdum(irnega) 
      b(2) = -y(ihe4) * ratdum(irneag)
      a1   = sum(b(1:2))
      call add_two(ihe4,ine20,a1,qdfdy,qsum,nt)

! d(he4)/d(mg24)
      b(1) =  ratdum(irmgga) 
      b(2) = -y(ihe4) * ratdum(irmgag) 
      b(3) = -y(ihe4) * ratdum(irmgap) * (1.0d0-ratdum(irr1))
      a1   = sum(b(1:3))
      call add_two(ihe4,img24,a1,qdfdy,qsum,nt)

! d(he4)/d(si28)
      b(1) =  ratdum(irsiga) 
      b(2) = -y(ihe4) * ratdum(irsiag) 
      b(3) = -y(ihe4) * ratdum(irsiap) * (1.0d0-ratdum(irs1)) 
      b(4) =  ratdum(irr1) * ratdum(irsigp)
      a1   = sum(b(1:4))
      call add_two(ihe4,isi28,a1,qdfdy,qsum,nt)

! d(he4)/d(s32)
      b(1) =  ratdum(irsga) 
      b(2) = -y(ihe4) * ratdum(irsag) 
      b(3) = -y(ihe4) * ratdum(irsap) * (1.0d0-ratdum(irt1)) 
      b(4) =  ratdum(irs1) * ratdum(irsgp)
      a1   = sum(b(1:4))
      call add_two(ihe4,is32,a1,qdfdy,qsum,nt)

! d(he4)/d(ar36)
      b(1)  =  ratdum(irarga) 
      b(2)  = -y(ihe4) * ratdum(irarag) 
      b(3)  = -y(ihe4) * ratdum(irarap) * (1.0d0-ratdum(iru1)) 
      b(4)  =  ratdum(irt1) * ratdum(irargp)
      a1   = sum(b(1:4))
      call add_two(ihe4,iar36,a1,qdfdy,qsum,nt)

! d(he4)/d(ca40)
      b(1) =  ratdum(ircaga) 
      b(2) = -y(ihe4) * ratdum(ircaag) 
      b(3) = -y(ihe4) * ratdum(ircaap) * (1.0d0-ratdum(irv1)) 
      b(4) =  ratdum(iru1) * ratdum(ircagp)
      a1   = sum(b(1:4))
      call add_two(ihe4,ica40,a1,qdfdy,qsum,nt)

! d(he4)/d(ti44)
      b(1) =  ratdum(irtiga) 
      b(2) = -y(ihe4) * ratdum(irtiag) 
      b(3) = -y(ihe4) * ratdum(irtiap) * (1.0d0-ratdum(irw1)) 
      b(4) =  ratdum(irv1) * ratdum(irtigp)
      a1   = sum(b(1:4))
      call add_two(ihe4,iti44,a1,qdfdy,qsum,nt)

! d(he4)/d(cr48)
      b(1) =  ratdum(ircrga) 
      b(2) = -y(ihe4) * ratdum(ircrag) 
      b(3) = -y(ihe4) * ratdum(ircrap) * (1.0d0-ratdum(irx1)) 
      b(4) =  ratdum(irw1) * ratdum(ircrgp)
      a1   = sum(b(1:4))
      call add_two(ihe4,icr48,a1,qdfdy,qsum,nt)

! d(he4)/d(fe52)
      b(1) =  ratdum(irfega) 
      b(2) = -y(ihe4) * ratdum(irfeag) 
      b(3) = -y(ihe4) * ratdum(irfeap) * (1.0d0-ratdum(iry1)) 
      b(4) =  ratdum(irx1) * ratdum(irfegp)
      a1   = sum(b(1:4))
      call add_two(ihe4,ife52,a1,qdfdy,qsum,nt)

! d(he4)/d(ni56)
      b(1) = ratdum(irniga) 
      b(2) = ratdum(iry1) * ratdum(irnigp)
      a1   = sum(b(1:2))
      call add_two(ihe4,ini56,a1,qdfdy,qsum,nt)


! c12 jacobian elements
! d(c12)/d(he4)
      b(1) =  0.5d0 * y(ihe4) * y(ihe4) * ratdum(ir3a) 
      b(2) = -y(ic12) * ratdum(ircag)
      a1   = sum(b(1:2))
      call add_two(ic12,ihe4,a1,qdfdy,qsum,nt)

! d(c12)/d(c12)
      b(1) = -2.0d0 * y(ic12) * ratdum(ir1212) 
      b(2) = -y(io16) * ratdum(ir1216) 
      b(3) = -ratdum(irg3a) 
      b(4) = -y(ihe4) * ratdum(ircag) 
      a1   = sum(b(1:4))
      call add_two(ic12,ic12,a1,qdfdy,qsum,nt)

! d(c12)/d(o16)
      b(1) = -y(ic12) * ratdum(ir1216) 
      b(2) =  ratdum(iroga)
      a1   = sum(b(1:2))
      call add_two(ic12,io16,a1,qdfdy,qsum,nt)



! o16 jacobian elements
! d(o16)/d(he4)
      b(1) =  y(ic12)*ratdum(ircag) 
      b(2) = -y(io16)*ratdum(iroag)
      a1   = sum(b(1:2))
      call add_two(io16,ihe4,a1,qdfdy,qsum,nt)

! d(o16)/d(c12)
      b(1) = -y(io16)*ratdum(ir1216) 
      b(2) =  y(ihe4)*ratdum(ircag)
      a1   = sum(b(1:2))
      call add_two(io16,ic12,a1,qdfdy,qsum,nt)

! d(o16)/d(o16)
      b(1) = -y(ic12) * ratdum(ir1216) 
      b(2) = -2.0d0 * y(io16) * ratdum(ir1616) 
      b(3) = -y(ihe4) * ratdum(iroag) 
      b(4) = -ratdum(iroga) 
      a1   = sum(b(1:4))
      call add_two(io16,io16,a1,qdfdy,qsum,nt)

! d(o16)/d(ne20)
      a1  = ratdum(irnega)
      call add_two(io16,ine20,a1,qdfdy,qsum,nt)



! ne20 jacobian elements
! d(ne20)/d(he4)
      b(1) =  y(io16) * ratdum(iroag) 
      b(2) = -y(ine20) * ratdum(irneag) 
      a1   = sum(b(1:2))
      call add_two(ine20,ihe4,a1,qdfdy,qsum,nt)

! d(ne20)/d(c12)
      a1  = y(ic12) * ratdum(ir1212)
      call add_two(ine20,ic12,a1,qdfdy,qsum,nt)

! d(ne20)/d(o16)
      a1  = y(ihe4) * ratdum(iroag)
      call add_two(ine20,io16,a1,qdfdy,qsum,nt)

! d(ne20)/d(ne20)
      b(1) = -y(ihe4) * ratdum(irneag) 
      b(2) = -ratdum(irnega)
      a1   = sum(b(1:2))
      call add_two(ine20,ine20,a1,qdfdy,qsum,nt)

! d(ne20)/d(mg24)
      a1  = ratdum(irmgga)
      call add_two(ine20,img24,a1,qdfdy,qsum,nt)


! mg24 jacobian elements
! d(mg24)/d(he4)
      b(1) =  y(ine20) * ratdum(irneag) 
      b(2) = -y(img24) * ratdum(irmgag) 
      b(3) = -y(img24) * ratdum(irmgap) * (1.0d0-ratdum(irr1))
      a1   = sum(b(1:3))
      call add_two(img24,ihe4,a1,qdfdy,qsum,nt)

! d(mg24)/d(c12)
      a1  = 0.5d0 * y(io16) * ratdum(ir1216)
      call add_two(img24,ic12,a1,qdfdy,qsum,nt)

! d(mg24)/d(o16)
      a1  = 0.5d0 * y(ic12) * ratdum(ir1216)
      call add_two(img24,io16,a1,qdfdy,qsum,nt)

! d(mg24)/d(ne20)
      a1  = y(ihe4) * ratdum(irneag)
      call add_two(img24,ine20,a1,qdfdy,qsum,nt)

! d(mg24)/d(mg24)
      b(1) = -y(ihe4) * ratdum(irmgag) 
      b(2) = -ratdum(irmgga) 
      b(3) = -y(ihe4) * ratdum(irmgap) * (1.0d0-ratdum(irr1))
      a1   = sum(b(1:3))
      call add_two(img24,img24,a1,qdfdy,qsum,nt)

! d(mg24)/d(si28)
      b(1) = ratdum(irsiga) 
      b(2) = ratdum(irr1) * ratdum(irsigp)
      a1   = sum(b(1:2))
      call add_two(img24,isi28,a1,qdfdy,qsum,nt)


! si28 jacobian elements
! d(si28)/d(he4)
      b(1) =  y(img24) * ratdum(irmgag) 
      b(2) = -y(isi28) * ratdum(irsiag) 
      b(3) =  y(img24) * ratdum(irmgap) * (1.0d0-ratdum(irr1)) 
      b(4) = -y(isi28) * ratdum(irsiap) * (1.0d0-ratdum(irs1))
      a1   = sum(b(1:4))
      call add_two(isi28,ihe4,a1,qdfdy,qsum,nt)

! d(si28)/d(c12)
      a1  =  0.5d0 * y(io16) * ratdum(ir1216)
      call add_two(isi28,ic12,a1,qdfdy,qsum,nt)

! d(si28)/d(o16)
      b(1) = 0.5d0 * y(ic12) * ratdum(ir1216) 
      b(2) = 1.12d0 * 0.5d0*y(io16) * ratdum(ir1616) 
      b(3) = 0.68d0 * 0.5d0*y(io16) * ratdum(irs1) * ratdum(ir1616)
      a1   = sum(b(1:3))
      call add_two(isi28,io16,a1,qdfdy,qsum,nt)

! d(si28)/d(mg24)
      b(1) =  y(ihe4) * ratdum(irmgag) 
      b(2) =  y(ihe4) * ratdum(irmgap) * (1.0d0-ratdum(irr1))
      a1   = sum(b(1:2))
      call add_two(isi28,img24,a1,qdfdy,qsum,nt)

! d(si28)/d(si28)
      b(1) =  -y(ihe4) * ratdum(irsiag) 
      b(2) = -ratdum(irsiga) 
      b(3) = -ratdum(irr1) * ratdum(irsigp) 
      b(4) = -y(ihe4) * ratdum(irsiap) * (1.0d0-ratdum(irs1))
      a1   = sum(b(1:4))
      call add_two(isi28,isi28,a1,qdfdy,qsum,nt)

! d(si28)/d(s32)
      b(1) = ratdum(irsga) 
      b(2) = ratdum(irs1) * ratdum(irsgp)
      a1   = sum(b(1:2))
      call add_two(isi28,is32,a1,qdfdy,qsum,nt)


! s32 jacobian elements
! d(s32)/d(he4)
      b(1) =  y(isi28) * ratdum(irsiag) 
      b(2) = -y(is32) * ratdum(irsag) 
      b(3) =  y(isi28) * ratdum(irsiap) * (1.0d0-ratdum(irs1)) 
      b(4) = -y(is32) * ratdum(irsap) * (1.0d0-ratdum(irt1))
      a1   = sum(b(1:4))
      call add_two(is32,ihe4,a1,qdfdy,qsum,nt)

! d(s32)/d(o16)
      b(1) = 0.68d0*0.5d0*y(io16)*ratdum(ir1616)*(1.0d0-ratdum(irs1)) 
      b(2) = 0.2d0 * 0.5d0*y(io16) * ratdum(ir1616)
      a1   = sum(b(1:2))
      call add_two(is32,io16,a1,qdfdy,qsum,nt)

! d(s32)/d(si28)
      b(1)  =y(ihe4) * ratdum(irsiag) 
      b(2) = y(ihe4) * ratdum(irsiap) * (1.0d0-ratdum(irs1))
      a1   = sum(b(1:2))
      call add_two(is32,isi28,a1,qdfdy,qsum,nt)

! d(s32)/d(s32)
      b(1) = -y(ihe4) * ratdum(irsag) 
      b(2) = -ratdum(irsga) 
      b(3) = -ratdum(irs1) * ratdum(irsgp) 
      b(4) = -y(ihe4) * ratdum(irsap) * (1.0d0-ratdum(irt1))
      a1   = sum(b(1:4))
      call add_two(is32,is32,a1,qdfdy,qsum,nt)

! d(s32)/d(ar36)
      b(1) = ratdum(irarga) 
      b(2) = ratdum(irt1) * ratdum(irargp)
      a1   = sum(b(1:2))
      call add_two(is32,iar36,a1,qdfdy,qsum,nt)


! ar36 jacobian elements
! d(ar36)/d(he4)
      b(1) =  y(is32)  * ratdum(irsag) 
      b(2) = -y(iar36) * ratdum(irarag) 
      b(3) =  y(is32)  * ratdum(irsap) * (1.0d0-ratdum(irt1)) 
      b(4) = -y(iar36) * ratdum(irarap) * (1.0d0-ratdum(iru1))
      a1   = sum(b(1:4))
      call add_two(iar36,ihe4,a1,qdfdy,qsum,nt)

! d(ar36)/d(s32)
      b(1) = y(ihe4) * ratdum(irsag) 
      b(2) = y(ihe4) * ratdum(irsap) * (1.0d0-ratdum(irt1))
      a1   = sum(b(1:2))
      call add_two(iar36,is32,a1,qdfdy,qsum,nt)

! d(ar36)/d(ar36)
      b(1) = -y(ihe4) * ratdum(irarag) 
      b(2) = -ratdum(irarga) 
      b(3) = -ratdum(irt1) * ratdum(irargp) 
      b(4) = -y(ihe4) * ratdum(irarap) * (1.0d0-ratdum(iru1))
      a1   = sum(b(1:4))
      call add_two(iar36,iar36,a1,qdfdy,qsum,nt)

! d(ar36)/d(ca40)
      b(1) = ratdum(ircaga) 
      b(2) = ratdum(ircagp) * ratdum(iru1)
      a1   = sum(b(1:2))
      call add_two(iar36,ica40,a1,qdfdy,qsum,nt)


! ca40 jacobian elements
! d(ca40)/d(he4)
      b(1)  =  y(iar36) * ratdum(irarag) 
      b(2)  = -y(ica40) * ratdum(ircaag) 
      b(3)  =  y(iar36) * ratdum(irarap)*(1.0d0-ratdum(iru1)) 
      b(4)  = -y(ica40) * ratdum(ircaap)*(1.0d0-ratdum(irv1))
      a1   = sum(b(1:4))
      call add_two(ica40,ihe4,a1,qdfdy,qsum,nt)

! d(ca40)/d(ar36)
      b(1) =  y(ihe4) * ratdum(irarag) 
      b(2) =  y(ihe4) * ratdum(irarap)*(1.0d0-ratdum(iru1))
      a1   = sum(b(1:2))
      call add_two(ica40,iar36,a1,qdfdy,qsum,nt)

! d(ca40)/d(ca40)
      b(1) = -y(ihe4) * ratdum(ircaag) 
      b(2) = -ratdum(ircaga) 
      b(3) = -ratdum(ircagp) * ratdum(iru1) 
      b(4) = -y(ihe4) * ratdum(ircaap)*(1.0d0-ratdum(irv1))
      a1   = sum(b(1:4))
      call add_two(ica40,ica40,a1,qdfdy,qsum,nt)

! d(ca40)/d(ti44)
      b(1) = ratdum(irtiga) 
      b(2) = ratdum(irtigp) * ratdum(irv1)
      a1   = sum(b(1:2))
      call add_two(ica40,iti44,a1,qdfdy,qsum,nt)



! ti44 jacobian elements
! d(ti44)/d(he4)
      b(1) =  y(ica40) * ratdum(ircaag) 
      b(2) = -y(iti44) * ratdum(irtiag) 
      b(3) =  y(ica40) * ratdum(ircaap)*(1.0d0-ratdum(irv1)) 
      b(4) = -y(iti44) * ratdum(irtiap)*(1.0d0-ratdum(irw1))
      a1   = sum(b(1:4))
      call add_two(iti44,ihe4,a1,qdfdy,qsum,nt)

! d(ti44)/d(ca40)
      b(1) =  y(ihe4) * ratdum(ircaag) 
      b(2) =  y(ihe4) * ratdum(ircaap)*(1.0d0-ratdum(irv1))
      a1   = sum(b(1:2))
      call add_two(iti44,ica40,a1,qdfdy,qsum,nt)

! d(ti44)/d(ti44)
      b(1) = -y(ihe4) * ratdum(irtiag) 
      b(2) = -ratdum(irtiga) 
      b(3) = -ratdum(irv1) * ratdum(irtigp) 
      b(4) = -y(ihe4) * ratdum(irtiap)*(1.0d0-ratdum(irw1))
      a1   = sum(b(1:4))
      call add_two(iti44,iti44,a1,qdfdy,qsum,nt)

! d(ti44)/d(cr48)
      b(1) = ratdum(ircrga) 
      b(2) = ratdum(irw1) * ratdum(ircrgp)
      a1   = sum(b(1:2))
      call add_two(iti44,icr48,a1,qdfdy,qsum,nt)



! cr48 jacobian elements
! d(cr48)/d(he4)
      b(1) =  y(iti44) * ratdum(irtiag) 
      b(2) = -y(icr48) * ratdum(ircrag) 
      b(3) =  y(iti44) * ratdum(irtiap)*(1.0d0-ratdum(irw1)) 
      b(4) = -y(icr48) * ratdum(ircrap)*(1.0d0-ratdum(irx1))
      a1   = sum(b(1:4))
      call add_two(icr48,ihe4,a1,qdfdy,qsum,nt)

! d(cr48)/d(ti44)
      b(1) =  y(ihe4) * ratdum(irtiag) 
      b(2) =  y(ihe4) * ratdum(irtiap)*(1.0d0-ratdum(irw1))
      a1   = sum(b(1:2))
      call add_two(icr48,iti44,a1,qdfdy,qsum,nt)

! d(cr48)/d(cr48)
      b(1) = -y(ihe4) * ratdum(ircrag) 
      b(2) = -ratdum(ircrga) 
      b(3) = -ratdum(irw1) * ratdum(ircrgp) 
      b(4) = -y(ihe4) * ratdum(ircrap)*(1.0d0-ratdum(irx1))
      a1   = sum(b(1:4))
      call add_two(icr48,icr48,a1,qdfdy,qsum,nt)

! d(cr48)/d(fe52)
      b(1) = ratdum(irfega) 
      b(2) = ratdum(irx1) * ratdum(irfegp)
      a1   = sum(b(1:2))
      call add_two(icr48,ife52,a1,qdfdy,qsum,nt)



! fe52 jacobian elements
! d(fe52)/d(he4)
      b(1) =  y(icr48) * ratdum(ircrag) 
      b(2) = -y(ife52) * ratdum(irfeag) 
      b(3) =  y(icr48) * ratdum(ircrap) * (1.0d0-ratdum(irx1)) 
      b(4) = -y(ife52) * ratdum(irfeap) * (1.0d0-ratdum(iry1)) 
      a1   = sum(b(1:4))
      call add_two(ife52,ihe4,a1,qdfdy,qsum,nt)

! d(fe52)/d(cr48)
      b(1) = y(ihe4) * ratdum(ircrag) 
      b(2) = y(ihe4) * ratdum(ircrap) * (1.0d0-ratdum(irx1))
      a1   = sum(b(1:2))
      call add_two(ife52,icr48,a1,qdfdy,qsum,nt)

! d(fe52)/d(fe52)
      b(1) = -y(ihe4) * ratdum(irfeag) 
      b(2) = -ratdum(irfega) 
      b(3) = -ratdum(irx1) * ratdum(irfegp) 
      b(4) = -y(ihe4) * ratdum(irfeap) * (1.0d0-ratdum(iry1))
      a1   = sum(b(1:4))
      call add_two(ife52,ife52,a1,qdfdy,qsum,nt)

! d(fe52)/d(ni56)
      b(1) = ratdum(irniga) 
      b(2) = ratdum(iry1) * ratdum(irnigp)
      a1   = sum(b(1:2))
      call add_two(ife52,ini56,a1,qdfdy,qsum,nt)


! ni56 jacobian elements
! d(ni56)/d(he4)
      b(1) =  y(ife52) * ratdum(irfeag) 
      b(2) =  y(ife52) * ratdum(irfeap) * (1.0d0-ratdum(iry1))
      a1   = sum(b(1:2))
      call add_two(ini56,ihe4,a1,qdfdy,qsum,nt)

! d(ni56)/d(fe52)
      b(1) = y(ihe4) * ratdum(irfeag) 
      b(2) = y(ihe4) * ratdum(irfeap) * (1.0d0-ratdum(iry1)) 
      a1   = sum(b(1:2))
      call add_two(ini56,ife52,a1,qdfdy,qsum,nt)

! d(ni56)/d(ni56)
      b(1) = -ratdum(irniga) 
      b(2) = -ratdum(iry1) * ratdum(irnigp)
      a1   = sum(b(1:2))
      call add_two(ini56,ini56,a1,qdfdy,qsum,nt)


! if we are doing a pure network, we are done, head to the error check

      if (pure_network .eq. 1) goto 678


! transfer the quad energy generation sum to double

      xsum(1:ionmax) = qsum(1:ionmax)


! temperature derivatives of the rate equations
! d(yi)/dtemp

      call rhs(y,dratdumdt,zwork1,deriva)
      do i=1,ionmax
       call add_one(i,itemp,zwork1(i)*inv_temp_scale,qdfdy,nt)
      enddo

! density derivatives of the rate equations
! d(yi)/d(den)

      call rhs(y,dratdumdd,zwork2,deriva)
      do i=1,ionmax
       call add_one(i,iden,zwork2(i)*inv_den_scale,qdfdy,nt)
      enddo

! energy jacobian elements
! d(iener)/d(yi)
      zz = enuc_conv2 * inv_energy_scale
      do i=1,ionmax
       call add_one(iener,i, zz*xsum(i),qdfdy,nt)
      enddo

! d(iener)/d(iener)
      call add_one(iener,iener, 0.0d0,qdfdy,nt)


! d(iener)/d(temp)
      x1 = sum(zwork1(1:ionmax)*mion(1:ionmax)) * enuc_conv2 *inv_energy_scale * inv_temp_scale
      call add_one(iener,itemp, x1,qdfdy,nt)
      dsdotdt = qdfdy(iat)


! d(iener)/d(den)
      x1 = sum(zwork2(1:ionmax)*mion(1:ionmax)) * enuc_conv2 * inv_energy_scale * inv_den_scale
      call add_one(iener,iden, x1,qdfdy,nt)
      dsdotdd = qdfdy(iat)


! d(iener)/d(ivelx)
       call add_one(iener,ivelx, 0.0d0,qdfdy,nt)


! d(iener)/d(iposx)
       call add_one(iener,iposx, 0.0d0,qdfdy,nt)





! account for the neutrino losses
      call sneut5(btemp,bden,abar,zbar,sneut,dsneutdt,dsneutdd,snuda,snudz)

! d(ener)/d(yi)
      do i=1,ionmax
       x1  = -((aion(i) - abar)*abar*snuda + (zion(i) - zbar)*abar*snudz) * inv_energy_scale
       call add_one(iener,i, x1,qdfdy,nt)
      enddo

! d(iener)/d(temp)
      x1 = -dsneutdt * inv_energy_scale * inv_temp_scale
      call add_one(iener,itemp, x1,qdfdy,nt)

! d(iener)/d(den)
      x1 = -dsneutdd * inv_energy_scale * inv_den_scale
      call add_one(iener,iden, a1,qdfdy,nt)




! the temperature and density jacobian elements
! depend on the burning mode


! hydrostatic
      if (hydrostatic  .or.  one_step  .or. trho_hist) then

       call add_one(itemp,itemp, 0.0d0,qdfdy,nt)
       call add_one(iden,iden,   0.0d0,qdfdy,nt)
       call add_one(ivelx,ivelx, 0.0d0,qdfdy,nt)
       call add_one(iposx,iposx, 0.0d0,qdfdy,nt)

! adiabatic expansion
      else if (expansion) then
       taud = 446.0d0/sqrt(den0)
       taut = 3.0d0 * taud
       x1  = -psi/taut 
       x2  = -psi/taud 

       call add_one(itemp,itemp, x1,qdfdy,nt)
       call add_one(iden,iden,   x2,qdfdy,nt)
       call add_one(ivelx,ivelx, 0.0d0,qdfdy,nt)
       call add_one(iposx,iposx, 0.0d0,qdfdy,nt)

! self heating
      else if (self_heat_const_den) then

! call an eos
       temp_row(1) = y(itemp)
       den_row(1)  = y(iden)
       abar_row(1) = abar
       zbar_row(1) = zbar
       jlo_eos = 1
       jhi_eos = 1

!       call helmeos
       call helmeosaprox13

! d(itemp)/d(yi)
       x1 = 1.0d0/cv_row(1) 
       x2  = enuc_conv2 * inv_temp_scale * x1
       x3 = sum(zwork1(1:ionmax)*mion(1:ionmax)) * enuc_conv2 * inv_temp_scale
       x4 = (x3 - dsneutdt) * x1
       x5 = sum(zwork2(1:ionmax)*mion(1:ionmax))*enuc_conv2 * inv_temp_scale
       x6 = (x5 - dsneutdd) * x1

       do i=1,ionmax
        call add_one(itemp,i, x2*xsum(i),qdfdy,nt)
       enddo
       call add_one(itemp,itemp, x4,qdfdy,nt)
       call add_one(itemp,iden,  x6,qdfdy,nt)
       call add_one(iden,iden,   0.0d0,qdfdy,nt)
       call add_one(ivelx,ivelx, 0.0d0,qdfdy,nt)
       call add_one(iposx,iposx, 0.0d0,qdfdy,nt)

! self heating at constant pressure
      else if (self_heat_const_pres) then

! temperature
       x1 = 1.0d0/cp_row(1)
       x2  = enuc_conv2 * inv_temp_scale * x1
       x3 = sum(zwork1(1:ionmax)*mion(1:ionmax)) * enuc_conv2 * inv_temp_scale
       x4 = (x3 - dsneutdt) * x1
       x5 = sum(zwork2(1:ionmax)*mion(1:ionmax))*enuc_conv2  * inv_temp_scale
       x6 = (x5 - dsneutdd) * x1

       do i=1,ionmax
        call add_one(itemp,i, x2*xsum(i),qdfdy,nt)
       enddo
       call add_one(itemp,itemp, x4,qdfdy,nt)
       call add_one(itemp,iden,  x6,qdfdy,nt)


! density
       x1 = 0.0d0
!       x1 = 1.0d0/cp_row(1)
       x2 = -dpt_row(1)/dpd_row(1) * x1  * enuc_conv2 * inv_den_scale
       x3 = sum(zwork1(1:ionmax)*mion(1:ionmax)) * enuc_conv2 
       x4 = (x3 - dsneutdt) * x1
       x5 = sum(zwork2(1:ionmax)*mion(1:ionmax)) * enuc_conv2 
       x6 = (x5 - dsneutdd) * x1

       do i=1,ionmax
        call add_one(iden,i, x2*xsum(i),qdfdy,nt)
       enddo
       call add_one(iden,itemp, x4,qdfdy,nt)
       call add_one(iden,iden,  x6,qdfdy,nt)

       call add_one(ivelx,ivelx, 0.0d0,qdfdy,nt)
       call add_one(iposx,iposx, 0.0d0,qdfdy,nt)

! end of burning mode ifs
      end if

678   continue


! set the double precision version of the return argument

      dfdy(1:nt) = qdfdy(1:nt)

! bullet check the counting
      if (nt .ne. nterms) then
       write(6,*) 'nt =',nt,'  nterms =',nterms
       write(6,*) 'error in routine saprox19: nt .ne. nterms'
       stop 'error in routine saprox19'
      end if


      return
      end subroutine saprox13





      subroutine aprox13rat
      include 'implno.dek'
      include 'burn_common.dek'
      include 'network.dek'

! this routine generates unscreened
! nuclear reaction rates for the aprox13 network.

! declare
      integer          i
      double precision rrate,drratedt,drratedd


! zero the rates
      ratraw(1:nrat) = 0.0d0
      dratrawdt(1:nrat) = 0.0d0
      dratrawdd(1:nrat) = 0.0d0


      if (btemp .lt. 1.0e6) return


! get the temperature factors
      call tfactors(btemp)


! triple alpha to c12
      call rate_tripalf(btemp,bden, &
           ratraw(ir3a),dratrawdt(ir3a),dratrawdd(ir3a), &
           ratraw(irg3a),dratrawdt(irg3a),dratrawdd(irg3a))

! c12(a,g)o16
      call rate_c12ag(btemp,bden, &
           ratraw(ircag),dratrawdt(ircag),dratrawdd(ircag), &
           ratraw(iroga),dratrawdt(iroga),dratrawdd(iroga))

! c12 + c12
      call rate_c12c12(btemp,bden, &
           ratraw(ir1212),dratrawdt(ir1212),dratrawdd(ir1212), &
           rrate,drratedt,drratedd)

! c12 + o16
      call rate_c12o16(btemp,bden, &
           ratraw(ir1216),dratrawdt(ir1216),dratrawdd(ir1216), &
           rrate,drratedt,drratedd)

! o16 + o16
      call rate_o16o16(btemp,bden, &
           ratraw(ir1616),dratrawdt(ir1616),dratrawdd(ir1616), &
           rrate,drratedt,drratedd)

! o16(a,g)ne20
      call rate_o16ag(btemp,bden, &
           ratraw(iroag),dratrawdt(iroag),dratrawdd(iroag), &
           ratraw(irnega),dratrawdt(irnega),dratrawdd(irnega))

! ne20(a,g)mg24
      call rate_ne20ag(btemp,bden, &
           ratraw(irneag),dratrawdt(irneag),dratrawdd(irneag), &
           ratraw(irmgga),dratrawdt(irmgga),dratrawdd(irmgga))

! mg24(a,g)si28
      call rate_mg24ag(btemp,bden, &
           ratraw(irmgag),dratrawdt(irmgag),dratrawdd(irmgag), &
           ratraw(irsiga),dratrawdt(irsiga),dratrawdd(irsiga))

! mg24(a,p)al27
      call rate_mg24ap(btemp,bden, &
           ratraw(irmgap),dratrawdt(irmgap),dratrawdd(irmgap), &
           ratraw(iralpa),dratrawdt(iralpa),dratrawdd(iralpa))

! al27(p,g)si28
      call rate_al27pg(btemp,bden, &
           ratraw(iralpg),dratrawdt(iralpg),dratrawdd(iralpg), &
           ratraw(irsigp),dratrawdt(irsigp),dratrawdd(irsigp))

! si28(a,g)s32
      call rate_si28ag(btemp,bden, &
           ratraw(irsiag),dratrawdt(irsiag),dratrawdd(irsiag), &
           ratraw(irsga),dratrawdt(irsga),dratrawdd(irsga))

! si28(a,p)p31
      call rate_si28ap(btemp,bden, &
           ratraw(irsiap),dratrawdt(irsiap),dratrawdd(irsiap), &
           ratraw(irppa),dratrawdt(irppa),dratrawdd(irppa))

! p31(p,g)s32
      call rate_p31pg(btemp,bden, &
           ratraw(irppg),dratrawdt(irppg),dratrawdd(irppg), &
           ratraw(irsgp),dratrawdt(irsgp),dratrawdd(irsgp))

! s32(a,g)ar36
      call rate_s32ag(btemp,bden, &
           ratraw(irsag),dratrawdt(irsag),dratrawdd(irsag), &
           ratraw(irarga),dratrawdt(irarga),dratrawdd(irarga))

! s32(a,p)cl35
      call rate_s32ap(btemp,bden, &
           ratraw(irsap),dratrawdt(irsap),dratrawdd(irsap), &
           ratraw(irclpa),dratrawdt(irclpa),dratrawdd(irclpa))

! cl35(p,g)ar36
      call rate_cl35pg(btemp,bden, &
           ratraw(irclpg),dratrawdt(irclpg),dratrawdd(irclpg), &
           ratraw(irargp),dratrawdt(irargp),dratrawdd(irargp))

! ar36(a,g)ca40
      call rate_ar36ag(btemp,bden, &
           ratraw(irarag),dratrawdt(irarag),dratrawdd(irarag), &
           ratraw(ircaga),dratrawdt(ircaga),dratrawdd(ircaga))

! ar36(a,p)k39
      call rate_ar36ap(btemp,bden, &
           ratraw(irarap),dratrawdt(irarap),dratrawdd(irarap), &
           ratraw(irkpa),dratrawdt(irkpa),dratrawdd(irkpa))

! k39(p,g)ca40
      call rate_k39pg(btemp,bden, &
           ratraw(irkpg),dratrawdt(irkpg),dratrawdd(irkpg), &
           ratraw(ircagp),dratrawdt(ircagp),dratrawdd(ircagp))

! ca40(a,g)ti44
      call rate_ca40ag(btemp,bden, &
           ratraw(ircaag),dratrawdt(ircaag),dratrawdd(ircaag), &
           ratraw(irtiga),dratrawdt(irtiga),dratrawdd(irtiga))

! ca40(a,p)sc43
      call rate_ca40ap(btemp,bden, &
           ratraw(ircaap),dratrawdt(ircaap),dratrawdd(ircaap), &
           ratraw(irscpa),dratrawdt(irscpa),dratrawdd(irscpa))

! sc43(p,g)ti44
      call rate_sc43pg(btemp,bden, &
           ratraw(irscpg),dratrawdt(irscpg),dratrawdd(irscpg), &
           ratraw(irtigp),dratrawdt(irtigp),dratrawdd(irtigp))

! ti44(a,g)cr48
      call rate_ti44ag(btemp,bden, &
           ratraw(irtiag),dratrawdt(irtiag),dratrawdd(irtiag), &
           ratraw(ircrga),dratrawdt(ircrga),dratrawdd(ircrga))

! ti44(a,p)v47
      call rate_ti44ap(btemp,bden, &
           ratraw(irtiap),dratrawdt(irtiap),dratrawdd(irtiap), &
           ratraw(irvpa),dratrawdt(irvpa),dratrawdd(irvpa))

! v47(p,g)cr48
      call rate_v47pg(btemp,bden, &
           ratraw(irvpg),dratrawdt(irvpg),dratrawdd(irvpg), &
           ratraw(ircrgp),dratrawdt(ircrgp),dratrawdd(ircrgp))

! cr48(a,g)fe52
      call rate_cr48ag(btemp,bden, &
           ratraw(ircrag),dratrawdt(ircrag),dratrawdd(ircrag), &
           ratraw(irfega),dratrawdt(irfega),dratrawdd(irfega))

! cr48(a,p)mn51
      call rate_cr48ap(btemp,bden, &
           ratraw(ircrap),dratrawdt(ircrap),dratrawdd(ircrap), &
           ratraw(irmnpa),dratrawdt(irmnpa),dratrawdd(irmnpa))

! mn51(p,g)fe52
      call rate_mn51pg(btemp,bden, &
           ratraw(irmnpg),dratrawdt(irmnpg),dratrawdd(irmnpg), &
           ratraw(irfegp),dratrawdt(irfegp),dratrawdd(irfegp))

! fe52(a,g)ni56
      call rate_fe52ag(btemp,bden, &
           ratraw(irfeag),dratrawdt(irfeag),dratrawdd(irfeag), &
           ratraw(irniga),dratrawdt(irniga),dratrawdd(irniga))

! fe52(a,p)co55
      call rate_fe52ap(btemp,bden, &
           ratraw(irfeap),dratrawdt(irfeap),dratrawdd(irfeap), &
           ratraw(ircopa),dratrawdt(ircopa),dratrawdd(ircopa))

! co55(p,g)ni56
      call rate_co55pg(btemp,bden, &
           ratraw(ircopg),dratrawdt(ircopg),dratrawdd(ircopg), &
           ratraw(irnigp),dratrawdt(irnigp),dratrawdd(irnigp))


! write out the rates
!      write(6,133) btemp,bden
! 133  format(1x,1pe12.4)
!      do i=1,nrat
!       write(6,134) ratnam(i),ratraw(i)
! 134   format(1x,a,'  ',1pe14.6,1pe11.3,1pe14.6)
!      enddo
!      read(5,*)



! for a strict alpha chain with only (a,g) and (g,a) reactions
! shut down the (a,p) (p,g) and (g,p) (p,a) rates

! irr1
!           ratraw(irmgag) = 0.0d0
!           ratraw(irsiga) = 0.0d0
!      ratraw(irmgap)    = 0.0d0  ; dratrawdt(irmgap) = 0.0d0  ; dratrawdd(irmgap) = 0.0d0
!      ratraw(iralpa)    = 0.0d0  ; dratrawdt(iralpa) = 0.0d0  ; dratrawdd(iralpa) = 0.0d0
!      ratraw(iralpg)    = 0.0d0  ; dratrawdt(iralpg) = 0.0d0  ; dratrawdd(iralpg) = 0.0d0
!      ratraw(irsigp)    = 0.0d0  ; dratrawdt(irsigp) = 0.0d0  ; dratrawdd(irsigp) = 0.0d0

! irs1
!           ratraw(ir1616) = 0.0d0
!           ratraw(irsiag) = 0.0d0
!           ratraw(irsga) = 0.0d0
!      ratraw(irsiap) = 0.0d0  ; dratrawdt(irsiap) = 0.0d0  ; dratrawdd(irsiap) = 0.0d0
!      ratraw(irppa)  = 0.0d0  ; dratrawdt(irppa)  = 0.0d0  ; dratrawdd(irppa)  = 0.0d0
!      ratraw(irppg)  = 0.0d0  ; dratrawdt(irppg)  = 0.0d0  ; dratrawdd(irppg)  = 0.0d0
!      ratraw(irsgp)  = 0.0d0  ; dratrawdt(irsgp)  = 0.0d0  ; dratrawdd(irsgp)  = 0.0d0

! irt1 
!           ratraw(irsag) = 0.0d0
!           ratraw(irarga) = 0.0d0
!      ratraw(irsap)  = 0.0d0  ;  dratrawdt(irsap)  = 0.0d0  ; dratrawdd(irsap)  = 0.0d0
!      ratraw(irclpa) = 0.0d0  ;  dratrawdt(irclpa) = 0.0d0  ; dratrawdd(irclpa) = 0.0d0
!      ratraw(irclpg) = 0.0d0  ;  dratrawdt(irclpg) = 0.0d0  ; dratrawdd(irclpg) = 0.0d0
!      ratraw(irargp) = 0.0d0  ;  dratrawdt(irargp) = 0.0d0  ; dratrawdd(irargp) = 0.0d0

! iru1 
!           ratraw(irarag) = 0.0d0
!           ratraw(ircaga) = 0.0d0
!      ratraw(irarap) = 0.0d0  ; dratrawdt(irarap) = 0.0d0  ; dratrawdd(irarap) = 0.0d0
!      ratraw(irkpa)  = 0.0d0  ;  dratrawdt(irkpa)  = 0.0d0 ; dratrawdd(irkpa)  = 0.0d0
!      ratraw(irkpg)  = 0.0d0  ; dratrawdt(irkpg)  = 0.0d0  ; dratrawdd(irkpg)  = 0.0d0
!      ratraw(ircagp) = 0.0d0  ; dratrawdt(ircagp) = 0.0d0  ; dratrawdd(ircagp) = 0.0d0

! irv1 
!           ratraw(ircaag) = 0.0d0
!           ratraw(irtiga) = 0.0d0
!      ratraw(ircaap) = 0.0d0  ; dratrawdt(ircaap) = 0.0d0  ; dratrawdd(ircaap) = 0.0d0
!      ratraw(irscpa) = 0.0d0  ; dratrawdt(irscpa) = 0.0d0  ; dratrawdd(irscpa) = 0.0d0
!      ratraw(irscpg) = 0.0d0  ; dratrawdt(irscpg) = 0.0d0  ; dratrawdd(irscpg) = 0.0d0
!      ratraw(irtigp) = 0.0d0  ; dratrawdt(irtigp) = 0.0d0  ; dratrawdd(irtigp) = 0.0d0

! irw1 
!           ratraw(irtiag) = 0.0d0
!           ratraw(ircrga) = 0.0d0
!      ratraw(irtiap) = 0.0d0  ; dratrawdt(irtiap) = 0.0d0  ; dratrawdd(irtiap) = 0.0d0
!      ratraw(irvpa)  = 0.0d0  ; dratrawdt(irvpa)  = 0.0d0  ; dratrawdd(irvpa)  = 0.0d0
!      ratraw(irvpg)  = 0.0d0  ; dratrawdt(irvpg)  = 0.0d0  ; dratrawdd(irvpg)  = 0.0d0
!      ratraw(ircrgp) = 0.0d0  ; dratrawdt(ircrgp) = 0.0d0  ; dratrawdd(ircrgp) = 0.0d0


! irx1 
!           ratraw(ircrag) = 0.0d0
!           ratraw(irfega) = 0.0d0
!      ratraw(ircrap) = 0.0d0  ; dratrawdt(ircrap) = 0.0d0  ; dratrawdd(ircrap) = 0.0d0
!      ratraw(irmnpa) = 0.0d0  ; dratrawdt(irmnpa) = 0.0d0  ; dratrawdd(irmnpa) = 0.0d0
!      ratraw(irmnpg) = 0.0d0  ; dratrawdt(irmnpg) = 0.0d0  ; dratrawdd(irmnpg) = 0.0d0
!      ratraw(irfegp) = 0.0d0  ; dratrawdt(irfegp) = 0.0d0  ; dratrawdd(irfegp) = 0.0d0


! iry1
!           ratraw(irfeag) = 0.0d0
!           ratraw(irniga) = 0.0d0
!      ratraw(irfeap) = 0.0d0  ; dratrawdt(irfeap) = 0.0d0  ; dratrawdd(irfeap) = 0.0d0
!      ratraw(ircopa) = 0.0d0  ; dratrawdt(ircopa) = 0.0d0  ; dratrawdd(ircopa) = 0.0d0
!      ratraw(ircopg) = 0.0d0  ; dratrawdt(ircopg) = 0.0d0  ; dratrawdd(ircopg) = 0.0d0
!      ratraw(irnigp) = 0.0d0  ; dratrawdt(irnigp) = 0.0d0  ; dratrawdd(irnigp) = 0.0d0



      return
      end subroutine aprox13rat





      subroutine aprox13tab
      include 'implno.dek'
      include 'burn_common.dek'
      include 'network.dek'

! uses tables instead of analytical expressions to evaluate the
! raw reaction rates. a cubic polynomial is hardwired for speed.

      integer          i,j,imax,iat,mp,ifirst
      parameter        (mp = 4)
      double precision tlo,thi,tstp,bden_sav,btemp_sav, &
                       x,x1,x2,x3,x4,a,b,c,d,e,f,g,h,p,q, &
                       alfa,beta,gama,delt

      double precision psi0,dpsi0,psi1,dpsi1,h3, &
                       z,zm1, &
                       t1,t2,t3,t4, &
                       ff1,ff2,ff3,ff4, &
                       g1,g2,g3,g4, &
                       xxbtemp

      integer          iswitch

      data             ifirst/0/


! cubic hermite basis functions and their derivatives as statement functions
      psi0(z)  = z * z * (2.0d0*z - 3.0d0) + 1.0
      dpsi0(z) = z * (6.0d0*z - 6.0d0)

      psi1(z)  = z * ( z * (z - 2.0d0) + 1.0d0)
      dpsi1(z) = z * (3.0d0*z - 4.0d0) + 1.0d0


! bicubic hermite polynomial as a statement function
      h3(ff1,t1,ff2,t2,ff3,t3,ff4,t4) = ff1*t1 + ff2*t2 + ff3*t3 + ff4*t4


! pick the interpolation method
      iswitch = 1


! make the table
      if (ifirst .eq. 0) then
       ifirst = 1


! set the log temperature loop limits
       thi  = 10.0d0
       tlo  = 6.0d0
       imax = int(thi-tlo)*per_decade + 1
       if (imax .gt. nrattab) stop 'imax too small in aprox13tab'
       tstp = (thi - tlo)/dfloat(imax-1)


! save the input
       btemp_sav = btemp
       bden_sav  = bden


! form the table
       bden = 1.0d0
       do i=1,imax
        btemp = tlo + dfloat(i-1)*tstp
        btemp = 10.0d0**(btemp)
        call aprox13rat
        ttab(i) = btemp
        do j=1,nrat
         rattab(j,i)    = ratraw(j)
         drattabdt(j,i) = dratrawdt(j)
         drattabdd(j,i) = dratrawdd(j)
        enddo
       enddo

! restore the input
       bden  = bden_sav
       btemp = btemp_sav
      end if


! normal execution starts here
! set the density dependence vector
      dtab(ircag)  = bden
      dtab(iroga)  = 1.0d0
      dtab(ir3a)   = bden*bden
      dtab(irg3a)  = 1.0d0
      dtab(ir1212) = bden
      dtab(ir1216) = bden
      dtab(ir1616) = bden
      dtab(iroag)  = bden
      dtab(irnega) = 1.0d0
      dtab(irneag) = bden
      dtab(irmgga) = 1.0d0
      dtab(irmgag) = bden
      dtab(irsiga) = 1.0d0
      dtab(irmgap) = bden
      dtab(iralpa) = bden
      dtab(iralpg) = bden
      dtab(irsigp) = 1.0d0
      dtab(irsiag) = bden
      dtab(irsga)  = 1.0d0
      dtab(irppa)  = bden
      dtab(irsiap) = bden
      dtab(irppg)  = bden
      dtab(irsgp)  = 1.0d0
      dtab(irsag)  = bden
      dtab(irarga) = 1.0d0
      dtab(irsap)  = bden
      dtab(irclpa) = bden
      dtab(irclpg) = bden
      dtab(irargp) = 1.0d0
      dtab(irarag) = bden
      dtab(ircaga) = 1.0d0
      dtab(irarap) = bden
      dtab(irkpa)  = bden
      dtab(irkpg)  = bden
      dtab(ircagp) = 1.0d0
      dtab(ircaag) = bden
      dtab(irtiga) = 1.0d0
      dtab(ircaap) = bden
      dtab(irscpa) = bden
      dtab(irscpg) = bden
      dtab(irtigp) = 1.0d0
      dtab(irtiag) = bden
      dtab(ircrga) = 1.0d0
      dtab(irtiap) = bden
      dtab(irvpa)  = bden
      dtab(irvpg)  = bden
      dtab(ircrgp) = 1.0d0
      dtab(ircrag) = bden
      dtab(irfega) = 1.0d0
      dtab(ircrap) = bden
      dtab(irmnpa) = bden
      dtab(irmnpg) = bden
      dtab(irfegp) = 1.0d0
      dtab(irfeag) = bden
      dtab(irniga) = 1.0d0
      dtab(irfeap) = bden
      dtab(ircopa) = bden
      dtab(ircopg) = bden
      dtab(irnigp) = 1.0d0



! cubic hermite interpolation
      if (iswitch .eq. 1) then

      xxbtemp = log10(btemp)

! hash locate the temperature
      iat = int((xxbtemp - tlo)/tstp) + 1
      iat = max(1,min(iat,imax - 1))

! blending functions 
       z   = max( (xxbtemp - log10(ttab(iat)))/tstp, 0.0d0)
       zm1 = 1.0d0 - z

! basis functions
       ff1  =  psi0(z)
       ff2  =  psi0(zm1)
       ff3  =  psi1(z)*tstp
       ff4  =  -psi1(zm1)*tstp

! basis function derivatives
       g1  =  dpsi0(z)/tstp
       g2  =  -dpsi0(zm1)/tstp
       g3  =  dpsi1(z)
       g4  =  dpsi1(zm1)


! crank off the raw reaction rates
      do j=1,nrat
       ratraw(j)    = h3(rattab(j,iat),ff1, rattab(j,iat+1),ff2, &
                         drattabdt(j,iat),ff3, drattabdt(j,iat+1),ff4) * dtab(j)

       dratrawdt(j) = h3(rattab(j,iat),g1, rattab(j,iat+1),g2, &
                         drattabdt(j,iat),g3, drattabdt(j,iat+1),g4) * dtab(j)

       dratrawdd(j) = h3(rattab(j,iat),g1, rattab(j,iat+1),g2, &
                         drattabdd(j,iat),g3, drattabdd(j,iat+1),g4)
      end do


! cubic lagrange interpolation
      else if (iswitch .eq. 2) then


! hash locate 
      iat = int((log10(btemp) - tlo)/tstp) + 1
      iat = max(1,min(iat - mp/2 + 1,imax - mp + 1))


! setup the lagrange interpolation coefficients for a cubic
      x  = btemp
      x1 = ttab(iat)
      x2 = ttab(iat+1)
      x3 = ttab(iat+2)
      x4 = ttab(iat+3)
      a  = x - x1
      b  = x - x2
      c  = x - x3
      d  = x - x4
      e  = x1 - x2
      f  = x1 - x3
      g  = x1 - x4
      h  = x2 - x3
      p  = x2 - x4
      q  = x3 - x4
      alfa =  b*c*d/(e*f*g)
      beta = -a*c*d/(e*h*p)
      gama =  a*b*d/(f*h*q)
      delt = -a*b*c/(g*p*q)

! crank off the raw reaction rates
      do j=1,nrat
       ratraw(j) = (alfa*rattab(j,iat) &
                  + beta*rattab(j,iat+1) &
                  + gama*rattab(j,iat+2) &
                  + delt*rattab(j,iat+3) &
                    ) * dtab(j)

       dratrawdt(j) = (alfa*drattabdt(j,iat) &
                     + beta*drattabdt(j,iat+1) &
                     + gama*drattabdt(j,iat+2) &
                     + delt*drattabdt(j,iat+3) &
                       ) * dtab(j)

       dratrawdd(j) =  alfa*drattabdd(j,iat) &
                     + beta*drattabdd(j,iat+1) &
                     + gama*drattabdd(j,iat+2) &
                     + delt*drattabdd(j,iat+3)

      enddo

! end of interpolation cases
      end if

! hand finish the three body reactions
      dratrawdd(ir3a) = bden * dratrawdd(ir3a)

      return
      end subroutine aprox13tab






      subroutine screen_aprox13(y)
      include 'implno.dek'
      include 'burn_common.dek'
      include 'network.dek'

! this routine computes the screening factors
! and applies them to the raw reaction rates,
! producing the final reaction rates used by the
! right hand sides and jacobian matrix elements

! this routine assumes screen_on = 1 or = 0 has been set at a higher
! level presumably in the top level driver


! declare the pass
      real*8  :: y(*)

! declare
      integer :: i,jscr,init=1
      real*8  :: sc1a,sc1adt,sc1add,sc2a,sc2adt,sc2add, &
                 sc3a,sc3adt,sc3add,abar,zbar,ye,z2bar, &
                 denom,denomdt,denomdd,xx,zz


! initialize
       ratdum(1:nrat) = ratraw(1:nrat)  ; dratdumdt(1:nrat) = dratrawdt(1:nrat)  ; dratdumdd(1:nrat) = dratrawdd(1:nrat)
       scfac(1:nrat)  = 1.0d0  ; dscfacdt(1:nrat)  = 0.0d0  ; dscfacdd(1:nrat)  = 0.0d0


! if screening is on
      if (screen_on .ne. 0) then


! with the passed composition compute abar, zbar and z**2 bar

       call azbar_simple(y,aion,zion,ionbeg,ionmax,abar,zbar,ye)
       z2bar  = abar * sum( zion(1:ionmax) * zion(1:ionmax) * y(1:ionmax) )


! first the always fun triple alpha and its inverse
      jscr = 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(ihe4),aion(ihe4),zion(ihe4),aion(ihe4), &
                   jscr,init,sc1a,sc1adt,sc1add)

      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(ihe4),aion(ihe4),4.0d0,8.0d0, &
                   jscr,init,sc2a,sc2adt,sc2add)

      sc3a   = sc1a * sc2a
      sc3adt = sc1adt*sc2a + sc1a*sc2adt
      sc3add = sc1add*sc2a + sc1a*sc2add

      ratdum(ir3a)    = ratraw(ir3a) * sc3a
      dratdumdt(ir3a) = dratrawdt(ir3a)*sc3a + ratraw(ir3a)*sc3adt
      dratdumdd(ir3a) = dratrawdd(ir3a)*sc3a + ratraw(ir3a)*sc3add

      scfac(ir3a) = sc3a  ; dscfacdt(ir3a)  = sc3adt  ; dscfacdd(ir3a)  = sc3add

      ratdum(irg3a)    = ratraw(irg3a) * sc3a
      dratdumdt(irg3a) = dratrawdt(irg3a)*sc3a + ratraw(irg3a)*sc3adt
      dratdumdd(irg3a) = dratrawdd(irg3a)*sc3a + ratraw(irg3a)*sc3add

      scfac(irg3a)  = sc3a  ; dscfacdt(irg3a)  = sc3adt  ; dscfacdd(irg3a)  = sc3add


! c12 to o16
! c12(a,g)o16
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(ic12),aion(ic12),zion(ihe4),aion(ihe4), &
                   jscr,init,sc1a,sc1adt,sc1add)

      ratdum(ircag)     = ratraw(ircag) * sc1a
      dratdumdt(ircag)  = dratrawdt(ircag)*sc1a + ratraw(ircag)*sc1adt
      dratdumdd(ircag)  = dratrawdd(ircag)*sc1a + ratraw(ircag)*sc1add

      scfac(ircag)  = sc1a  ; dscfacdt(ircag)   = sc1adt  ; dscfacdt(ircag)   = sc1add

      ratdum(iroga)     = ratraw(iroga) * sc1a
      dratdumdt(iroga)  = dratrawdt(iroga)*sc1a + ratraw(iroga)*sc1adt
      dratdumdd(iroga)  = dratrawdd(iroga)*sc1a + ratraw(iroga)*sc1add

      scfac(iroga)  = sc1a  ; dscfacdt(iroga)   = sc1adt  ; dscfacdt(iroga)   = sc1add


! c12 + c12
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(ic12),aion(ic12),zion(ic12),aion(ic12), &
                   jscr,init,sc1a,sc1adt,sc1add)


      ratdum(ir1212)    = ratraw(ir1212) * sc1a
      dratdumdt(ir1212) = dratrawdt(ir1212)*sc1a + ratraw(ir1212)*sc1adt
      dratdumdd(ir1212) = dratrawdd(ir1212)*sc1a + ratraw(ir1212)*sc1add

      scfac(ir1212)     = sc1a  ; dscfacdt(ir1212)  = sc1adt  ; dscfacdd(ir1212)  = sc1add



! c12 + o16
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(ic12),aion(ic12),zion(io16),aion(io16), &
                   jscr,init,sc1a,sc1adt,sc1add)

      ratdum(ir1216)    = ratraw(ir1216) * sc1a
      dratdumdt(ir1216) = dratrawdt(ir1216)*sc1a + ratraw(ir1216)*sc1adt
      dratdumdd(ir1216) = dratrawdd(ir1216)*sc1a + ratraw(ir1216)*sc1add

      scfac(ir1216)     = sc1a  ; dscfacdt(ir1216)  = sc1adt  ; dscfacdd(ir1216)  = sc1add


! o16 + o16
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(io16),aion(io16),zion(io16),aion(io16), &
                   jscr,init,sc1a,sc1adt,sc1add)


      ratdum(ir1616)    = ratraw(ir1616) * sc1a
      dratdumdt(ir1616) = dratrawdt(ir1616)*sc1a + ratraw(ir1616)*sc1adt
      dratdumdd(ir1616) = dratrawdd(ir1616)*sc1a + ratraw(ir1616)*sc1add

      scfac(ir1616)     = sc1a  ; dscfacdt(ir1616)  = sc1adt  ; dscfacdd(ir1616)  = sc1add



! o16 to ne20
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(io16),aion(io16),zion(ihe4),aion(ihe4), &
                   jscr,init,sc1a,sc1adt,sc1add)


! o16(a,g)ne20
      ratdum(iroag)    = ratraw(iroag) * sc1a
      dratdumdt(iroag) = dratrawdt(iroag)*sc1a + ratraw(iroag)*sc1adt
      dratdumdd(iroag) = dratrawdd(iroag)*sc1a + ratraw(iroag)*sc1add

      scfac(iroag)  = sc1a  ; dscfacdt(iroag)  = sc1adt  ; dscfacdd(iroag)  = sc1add

      ratdum(irnega)    = ratraw(irnega) * sc1a
      dratdumdt(irnega) = dratrawdt(irnega)*sc1a + ratraw(irnega)*sc1adt
      dratdumdd(irnega) = dratrawdd(irnega)*sc1a + ratraw(irnega)*sc1add

      scfac(irnega)  = sc1a  ; dscfacdt(irnega)  = sc1adt  ; dscfacdd(irnega)  = sc1add


! ne20 to mg24
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(ine20),aion(ine20),zion(ihe4),aion(ihe4), &
                   jscr,init,sc1a,sc1adt,sc1add)


! ne20(a,g)mg24
      ratdum(irneag)    = ratraw(irneag) * sc1a
      dratdumdt(irneag) = dratrawdt(irneag)*sc1a + ratraw(irneag)*sc1adt
      dratdumdd(irneag) = dratrawdd(irneag)*sc1a + ratraw(irneag)*sc1add

      scfac(irneag) = sc1a  ; dscfacdt(irneag)  = sc1adt  ; dscfacdd(irneag)  = sc1add

      ratdum(irmgga)    = ratraw(irmgga) * sc1a
      dratdumdt(irmgga) = dratrawdt(irmgga)*sc1a + ratraw(irmgga)*sc1adt
      dratdumdd(irmgga) = dratrawdd(irmgga)*sc1a + ratraw(irmgga)*sc1add

      scfac(irmgga) = sc1a  ; dscfacdt(irmgga)  = sc1adt  ; dscfacdd(irmgga)  = sc1add



! mg24 to si28
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(img24),aion(img24),zion(ihe4),aion(ihe4), &
                   jscr,init,sc1a,sc1adt,sc1add)


! mg24(a,g)si28
      ratdum(irmgag)    = ratraw(irmgag) * sc1a
      dratdumdt(irmgag) = dratrawdt(irmgag)*sc1a + ratraw(irmgag)*sc1adt
      dratdumdd(irmgag) = dratrawdd(irmgag)*sc1a + ratraw(irmgag)*sc1add

      scfac(irmgag) = sc1a ; dscfacdt(irmgag)  = sc1adt  ; dscfacdd(irmgag)  = sc1add

      ratdum(irsiga)    = ratraw(irsiga) * sc1a
      dratdumdt(irsiga) = dratrawdt(irsiga)*sc1a + ratraw(irsiga)*sc1adt
      dratdumdd(irsiga) = dratrawdd(irsiga)*sc1a + ratraw(irsiga)*sc1add

      scfac(irsiga) = sc1a ; dscfacdt(irsiga)  = sc1adt  ; dscfacdd(irsiga)  = sc1add


! mg24(a,p)al27
      ratdum(irmgap)    = ratraw(irmgap) * sc1a
      dratdumdt(irmgap) = dratrawdt(irmgap)*sc1a + ratraw(irmgap)*sc1adt
      dratdumdd(irmgap) = dratrawdd(irmgap)*sc1a + ratraw(irmgap)*sc1add

      scfac(irmgap)     = sc1a  ; dscfacdt(irmgap)  = sc1adt  ; dscfacdd(irmgap)  = sc1add

      ratdum(iralpa)    = ratraw(iralpa) * sc1a
      dratdumdt(iralpa) = dratrawdt(iralpa)*sc1a + ratraw(iralpa)*sc1adt
      dratdumdd(iralpa) = dratrawdd(iralpa)*sc1a + ratraw(iralpa)*sc1add

      scfac(iralpa)     = sc1a  ; dscfacdt(iralpa)  = sc1adt  ; dscfacdd(iralpa)  = sc1add


      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   13.0d0,27.0d0,1.0d0,1.0d0, &
                   jscr,init,sc1a,sc1adt,sc1add)

! al27(p,g)si28
      ratdum(iralpg)    = ratraw(iralpg) * sc1a
      dratdumdt(iralpg) = dratrawdt(iralpg)*sc1a + ratraw(iralpg)*sc1adt
      dratdumdd(iralpg) = dratrawdd(iralpg)*sc1a + ratraw(iralpg)*sc1add

      scfac(iralpg)     = sc1a  ; dscfacdt(iralpg)  = sc1adt  ; dscfacdd(iralpg)  = sc1add

      ratdum(irsigp)    = ratraw(irsigp) * sc1a
      dratdumdt(irsigp) = dratrawdt(irsigp)*sc1a + ratraw(irsigp)*sc1adt
      dratdumdd(irsigp) = dratrawdd(irsigp)*sc1a + ratraw(irsigp)*sc1add

      scfac(irsigp)     = sc1a  ; dscfacdt(irsigp)  = sc1adt  ; dscfacdd(irsigp)  = sc1add


! si28 to s32
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(isi28),aion(isi28),zion(ihe4),aion(ihe4), &
                   jscr,init,sc1a,sc1adt,sc1add)


! si28(a,g)s32
      ratdum(irsiag)    = ratraw(irsiag) * sc1a
      dratdumdt(irsiag) = dratrawdt(irsiag)*sc1a + ratraw(irsiag)*sc1adt
      dratdumdd(irsiag) = dratrawdd(irsiag)*sc1a + ratraw(irsiag)*sc1add

      scfac(irsiag)     = sc1a  ; dscfacdt(irsiag)  = sc1adt  ; dscfacdd(irsiag)  = sc1add

      ratdum(irsga)    = ratraw(irsga) * sc1a
      dratdumdt(irsga) = dratrawdt(irsga)*sc1a + ratraw(irsga)*sc1adt
      dratdumdd(irsga) = dratrawdd(irsga)*sc1a + ratraw(irsga)*sc1add

      scfac(irsga)     = sc1a  ; dscfacdt(irsga)  = sc1adt  ; dscfacdd(irsga)  = sc1add


! si28(a,p)p31
      ratdum(irsiap)    = ratraw(irsiap) * sc1a
      dratdumdt(irsiap) = dratrawdt(irsiap)*sc1a + ratraw(irsiap)*sc1adt
      dratdumdd(irsiap) = dratrawdd(irsiap)*sc1a + ratraw(irsiap)*sc1add

      scfac(irsiap)     = sc1a  ; dscfacdt(irsiap)  = sc1adt  ; dscfacdd(irsiap)  = sc1add

      ratdum(irppa)     = ratraw(irppa) * sc1a
      dratdumdt(irppa)  = dratrawdt(irppa)*sc1a  + ratraw(irppa)*sc1adt
      dratdumdd(irppa)  = dratrawdd(irppa)*sc1a  + ratraw(irppa)*sc1add

      scfac(irppa)      = sc1a  ; dscfacdt(irppa)   = sc1adt  ; dscfacdd(irppa)   = sc1add


      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   15.0d0,31.0d0,1.0d0,1.0d0, &
                   jscr,init,sc1a,sc1adt,sc1add)

! p31(p,g)s32
      ratdum(irppg)     = ratraw(irppg) * sc1a
      dratdumdt(irppg)  = dratrawdt(irppg)*sc1a + ratraw(irppg)*sc1adt
      dratdumdd(irppg)  = dratrawdd(irppg)*sc1a + ratraw(irppg)*sc1add

      scfac(irppg)      = sc1a  ; dscfacdt(irppg)   = sc1adt  ; dscfacdd(irppg)   = sc1add

      ratdum(irsgp)     = ratraw(irsgp) * sc1a
      dratdumdt(irsgp)  = dratrawdt(irsgp)*sc1a + ratraw(irsgp)*sc1adt
      dratdumdd(irsgp)  = dratrawdd(irsgp)*sc1a + ratraw(irsgp)*sc1add

      scfac(irsgp)      = sc1a  ; dscfacdt(irsgp)   = sc1adt  ; dscfacdd(irsgp)   = sc1add



! s32 to ar36
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(is32),aion(is32),zion(ihe4),aion(ihe4), &
                   jscr,init,sc1a,sc1adt,sc1add)


! s32(a,g)ar36
      ratdum(irsag)     = ratraw(irsag) * sc1a
      dratdumdt(irsag)  = dratrawdt(irsag)*sc1a + ratraw(irsag)*sc1adt
      dratdumdd(irsag)  = dratrawdd(irsag)*sc1a + ratraw(irsag)*sc1add

      scfac(irsag)      = sc1a ; dscfacdt(irsag)   = sc1adt ; dscfacdd(irsag)   = sc1add

      ratdum(irarga)     = ratraw(irarga) * sc1a
      dratdumdt(irarga)  = dratrawdt(irarga)*sc1a + ratraw(irarga)*sc1adt
      dratdumdd(irarga)  = dratrawdd(irarga)*sc1a + ratraw(irarga)*sc1add

      scfac(irarga)      = sc1a ; dscfacdt(irarga)   = sc1adt ; dscfacdd(irarga)   = sc1add

! s32(a,p)cl35
      ratdum(irsap)     = ratraw(irsap) * sc1a
      dratdumdt(irsap)  = dratrawdt(irsap)*sc1a + ratraw(irsap)*sc1adt
      dratdumdd(irsap)  = dratrawdd(irsap)*sc1a + ratraw(irsap)*sc1add

      scfac(irsap)      = sc1a  ; dscfacdt(irsap)   = sc1adt  ; dscfacdd(irsap)   = sc1add

      ratdum(irclpa)    = ratraw(irclpa) * sc1a
      dratdumdt(irclpa) = dratrawdt(irclpa)*sc1a + ratraw(irclpa)*sc1adt
      dratdumdd(irclpa) = dratrawdd(irclpa)*sc1a + ratraw(irclpa)*sc1add

      scfac(irclpa)     = sc1a  ; dscfacdt(irclpa)  = sc1adt  ; dscfacdt(irclpa)  = sc1add


      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   17.0d0,35.0d0,1.0d0,1.0d0, &
                   jscr,init,sc1a,sc1adt,sc1add)

! cl35(p,g)ar36
      ratdum(irclpg)    = ratraw(irclpg) * sc1a
      dratdumdt(irclpg) = dratrawdt(irclpg)*sc1a + ratraw(irclpg)*sc1adt
      dratdumdd(irclpg) = dratrawdd(irclpg)*sc1a + ratraw(irclpg)*sc1add

      scfac(irclpg)     = sc1a  ; dscfacdt(irclpg)  = sc1adt  ; dscfacdd(irclpg)  = sc1add

      ratdum(irargp)    = ratraw(irargp) * sc1a
      dratdumdt(irargp) = dratrawdt(irargp)*sc1a + ratraw(irargp)*sc1adt
      dratdumdd(irargp) = dratrawdd(irargp)*sc1a + ratraw(irargp)*sc1add

      scfac(irargp)     = sc1a  ; dscfacdt(irargp)  = sc1adt  ; dscfacdd(irargp)  = sc1add



! ar36 to ca40
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(iar36),aion(iar36),zion(ihe4),aion(ihe4), &
                   jscr,init,sc1a,sc1adt,sc1add)


! ar36(a,g)ca40
      ratdum(irarag)    = ratraw(irarag) * sc1a
      dratdumdt(irarag) = dratrawdt(irarag)*sc1a + ratraw(irarag)*sc1adt
      dratdumdd(irarag) = dratrawdd(irarag)*sc1a + ratraw(irarag)*sc1add

      scfac(irarag)     = sc1a  ; dscfacdt(irarag)  = sc1adt  ; dscfacdd(irarag)  = sc1add

      ratdum(ircaga)    = ratraw(ircaga) * sc1a
      dratdumdt(ircaga) = dratrawdt(ircaga)*sc1a + ratraw(ircaga)*sc1adt
      dratdumdd(ircaga) = dratrawdd(ircaga)*sc1a + ratraw(ircaga)*sc1add

      scfac(ircaga)     = sc1a  ; dscfacdt(ircaga)  = sc1adt  ; dscfacdd(ircaga)  = sc1add


! ar36(a,p)k39
      ratdum(irarap)    = ratraw(irarap) * sc1a
      dratdumdt(irarap) = dratrawdt(irarap)*sc1a + ratraw(irarap)*sc1adt
      dratdumdd(irarap) = dratrawdd(irarap)*sc1a + ratraw(irarap)*sc1add

      scfac(irarap)     = sc1a  ; dscfacdt(irarap)  = sc1adt  ; dscfacdd(irarap)  = sc1add

      ratdum(irkpa)     = ratraw(irkpa) * sc1a
      dratdumdt(irkpa)  = dratrawdt(irkpa)*sc1a  + ratraw(irkpa)*sc1adt
      dratdumdd(irkpa)  = dratrawdd(irkpa)*sc1a  + ratraw(irkpa)*sc1add

      scfac(irkpa)      = sc1a  ; dscfacdt(irkpa)   = sc1adt  ; dscfacdd(irkpa)   = sc1add


      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   19.0d0,39.0d0,1.0d0,1.0d0, &
                   jscr,init,sc1a,sc1adt,sc1add)

! k39(p,g)ca40
      ratdum(irkpg)     = ratraw(irkpg) * sc1a
      dratdumdt(irkpg)  = dratrawdt(irkpg)*sc1a  + ratraw(irkpg)*sc1adt
      dratdumdd(irkpg)  = dratrawdd(irkpg)*sc1a  + ratraw(irkpg)*sc1add

      scfac(irkpg)      = sc1a  ; dscfacdt(irkpg)   = sc1adt  ; dscfacdd(irkpg)   = sc1add

      ratdum(ircagp)     = ratraw(ircagp) * sc1a
      dratdumdt(ircagp)  = dratrawdt(ircagp)*sc1a  + ratraw(ircagp)*sc1adt
      dratdumdd(ircagp)  = dratrawdd(ircagp)*sc1a  + ratraw(ircagp)*sc1add

      scfac(ircagp)      = sc1a  ; dscfacdt(ircagp)   = sc1adt  ; dscfacdd(ircagp)   = sc1add



! ca40 to ti44
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(ica40),aion(ica40),zion(ihe4),aion(ihe4), &
                   jscr,init,sc1a,sc1adt,sc1add)


! ca40(a,g)ti44
      ratdum(ircaag)    = ratraw(ircaag) * sc1a
      dratdumdt(ircaag) = dratrawdt(ircaag)*sc1a + ratraw(ircaag)*sc1adt
      dratdumdd(ircaag) = dratrawdd(ircaag)*sc1a + ratraw(ircaag)*sc1add

      scfac(ircaag)     = sc1a  ; dscfacdt(ircaag)  = sc1adt  ; dscfacdd(ircaag)  = sc1add

      ratdum(irtiga)    = ratraw(irtiga) * sc1a
      dratdumdt(irtiga) = dratrawdt(irtiga)*sc1a + ratraw(irtiga)*sc1adt
      dratdumdd(irtiga) = dratrawdd(irtiga)*sc1a + ratraw(irtiga)*sc1add

      scfac(irtiga)     = sc1a  ; dscfacdt(irtiga)  = sc1adt  ; dscfacdd(irtiga)  = sc1add


! ca40(a,p)sc43
      ratdum(ircaap)    = ratraw(ircaap) * sc1a
      dratdumdt(ircaap) = dratrawdt(ircaap)*sc1a + ratraw(ircaap)*sc1adt
      dratdumdd(ircaap) = dratrawdd(ircaap)*sc1a + ratraw(ircaap)*sc1add

      scfac(ircaap)     = sc1a  ; dscfacdt(ircaap)  = sc1adt  ; dscfacdd(ircaap)  = sc1add

      ratdum(irscpa)    = ratraw(irscpa) * sc1a
      dratdumdt(irscpa) = dratrawdt(irscpa)*sc1a + ratraw(irscpa)*sc1adt
      dratdumdd(irscpa) = dratrawdd(irscpa)*sc1a + ratraw(irscpa)*sc1add

      scfac(irscpa)     = sc1a  ; dscfacdt(irscpa)  = sc1adt  ; dscfacdd(irscpa)  = sc1add


      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   21.0d0,43.0d0,1.0d0,1.0d0, &
                   jscr,init,sc1a,sc1adt,sc1add)

! sc43(p,g)ti44
      ratdum(irscpg)    = ratraw(irscpg) * sc1a
      dratdumdt(irscpg) = dratrawdt(irscpg)*sc1a + ratraw(irscpg)*sc1adt
      dratdumdd(irscpg) = dratrawdd(irscpg)*sc1a + ratraw(irscpg)*sc1add

      scfac(irscpg)     = sc1a  ; dscfacdt(irscpg)  = sc1adt  ; dscfacdd(irscpg)  = sc1add

      ratdum(irtigp)    = ratraw(irtigp) * sc1a
      dratdumdt(irtigp) = dratrawdt(irtigp)*sc1a + ratraw(irtigp)*sc1adt
      dratdumdd(irtigp) = dratrawdd(irtigp)*sc1a + ratraw(irtigp)*sc1add

      scfac(irtigp)     = sc1a  ; dscfacdt(irtigp)  = sc1adt  ; dscfacdd(irtigp)  = sc1add



! ti44 to cr48
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(iti44),aion(iti44),zion(ihe4),aion(ihe4), &
                   jscr,init,sc1a,sc1adt,sc1add)


! ti44(a,g)cr48
      ratdum(irtiag)    = ratraw(irtiag) * sc1a
      dratdumdt(irtiag) = dratrawdt(irtiag)*sc1a + ratraw(irtiag)*sc1adt
      dratdumdd(irtiag) = dratrawdd(irtiag)*sc1a + ratraw(irtiag)*sc1add

      scfac(irtiag)     = sc1a  ; dscfacdt(irtiag)  = sc1adt  ; dscfacdd(irtiag)  = sc1add

      ratdum(ircrga)    = ratraw(ircrga) * sc1a
      dratdumdt(ircrga) = dratrawdt(ircrga)*sc1a + ratraw(ircrga)*sc1adt
      dratdumdd(ircrga) = dratrawdd(ircrga)*sc1a + ratraw(ircrga)*sc1add

      scfac(ircrga)     = sc1a  ; dscfacdt(ircrga)  = sc1adt  ; dscfacdd(ircrga)  = sc1add

! ti44(a,p)v47
      ratdum(irtiap)    = ratraw(irtiap) * sc1a
      dratdumdt(irtiap) = dratrawdt(irtiap)*sc1a + ratraw(irtiap)*sc1adt
      dratdumdd(irtiap) = dratrawdd(irtiap)*sc1a + ratraw(irtiap)*sc1add

      scfac(irtiap)  = sc1a  ; dscfacdt(irtiap)  = sc1adt  ; dscfacdd(irtiap)  = sc1add

      ratdum(irvpa)     = ratraw(irvpa) * sc1a
      dratdumdt(irvpa)  = dratrawdt(irvpa)*sc1a  + ratraw(irvpa)*sc1adt
      dratdumdd(irvpa)  = dratrawdd(irvpa)*sc1a  + ratraw(irvpa)*sc1add

      scfac(irvpa)      = sc1a  ; dscfacdt(irvpa)   = sc1adt  ; dscfacdd(irvpa)   = sc1add


      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   23.0d0,47.0d0,1.0d0,1.0d0, &
                   jscr,init,sc1a,sc1adt,sc1add)

! v47(p,g)cr48
      ratdum(irvpg)     = ratraw(irvpg) * sc1a
      dratdumdt(irvpg)  = dratrawdt(irvpg)*sc1a  + ratraw(irvpg)*sc1adt
      dratdumdd(irvpg)  = dratrawdd(irvpg)*sc1a  + ratraw(irvpg)*sc1add

      scfac(irvpg)      = sc1a  ; dscfacdt(irvpg)   = sc1adt  ; dscfacdd(irvpg)   = sc1add

      ratdum(ircrgp)     = ratraw(ircrgp) * sc1a
      dratdumdt(ircrgp)  = dratrawdt(ircrgp)*sc1a  + ratraw(ircrgp)*sc1adt
      dratdumdd(ircrgp)  = dratrawdd(ircrgp)*sc1a  + ratraw(ircrgp)*sc1add

      scfac(ircrgp)      = sc1a  ; dscfacdt(ircrgp)   = sc1adt  ; dscfacdd(ircrgp)   = sc1add



! cr48 to fe52
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(icr48),aion(icr48),zion(ihe4),aion(ihe4), &
                   jscr,init,sc1a,sc1adt,sc1add)


! cr48(a,g)fe52
      ratdum(ircrag)    = ratraw(ircrag) * sc1a
      dratdumdt(ircrag) = dratrawdt(ircrag)*sc1a + ratraw(ircrag)*sc1adt
      dratdumdd(ircrag) = dratrawdd(ircrag)*sc1a + ratraw(ircrag)*sc1add

      scfac(ircrag)     = sc1a  ; dscfacdt(ircrag)  = sc1adt  ; dscfacdd(ircrag)  = sc1add

      ratdum(irfega)    = ratraw(irfega) * sc1a
      dratdumdt(irfega) = dratrawdt(irfega)*sc1a + ratraw(irfega)*sc1adt
      dratdumdd(irfega) = dratrawdd(irfega)*sc1a + ratraw(irfega)*sc1add

      scfac(irfega)     = sc1a  ; dscfacdt(irfega)  = sc1adt  ; dscfacdd(irfega)  = sc1add


! cr48(a,p)mn51
      ratdum(ircrap)    = ratraw(ircrap) * sc1a
      dratdumdt(ircrap) = dratrawdt(ircrap)*sc1a + ratraw(ircrap)*sc1adt
      dratdumdd(ircrap) = dratrawdd(ircrap)*sc1a + ratraw(ircrap)*sc1add

      scfac(ircrap)     = sc1a  ; dscfacdt(ircrap)  = sc1adt  ; dscfacdd(ircrap)  = sc1add

      ratdum(irmnpa)    = ratraw(irmnpa) * sc1a
      dratdumdt(irmnpa) = dratrawdt(irmnpa)*sc1a + ratraw(irmnpa)*sc1adt
      dratdumdd(irmnpa) = dratrawdd(irmnpa)*sc1a + ratraw(irmnpa)*sc1add

      scfac(irmnpa)     = sc1a  ; dscfacdt(irmnpa)  = sc1adt  ; dscfacdd(irmnpa)  = sc1add


      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   25.0d0,51.0d0,1.0d0,1.0d0, &
                   jscr,init,sc1a,sc1adt,sc1add)

! mn51(p,g)fe52
      ratdum(irmnpg)    = ratraw(irmnpg) * sc1a
      dratdumdt(irmnpg) = dratrawdt(irmnpg)*sc1a + ratraw(irmnpg)*sc1adt
      dratdumdd(irmnpg) = dratrawdd(irmnpg)*sc1a + ratraw(irmnpg)*sc1add

      scfac(irmnpg)     = sc1a  ; dscfacdt(irmnpg)  = sc1adt  ; dscfacdd(irmnpg)  = sc1add

      ratdum(irfegp)    = ratraw(irfegp) * sc1a
      dratdumdt(irfegp) = dratrawdt(irfegp)*sc1a + ratraw(irfegp)*sc1adt
      dratdumdd(irfegp) = dratrawdd(irfegp)*sc1a + ratraw(irfegp)*sc1add

      scfac(irfegp)     = sc1a  ; dscfacdt(irfegp)  = sc1adt  ; dscfacdd(irfegp)  = sc1add



! fe52 to ni56
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(ife52),aion(ife52),zion(ihe4),aion(ihe4), &
                   jscr,init,sc1a,sc1adt,sc1add)


! fe52(a,g)ni56
      ratdum(irfeag)    = ratraw(irfeag) * sc1a
      dratdumdt(irfeag) = dratrawdt(irfeag)*sc1a + ratraw(irfeag)*sc1adt
      dratdumdd(irfeag) = dratrawdd(irfeag)*sc1a + ratraw(irfeag)*sc1add

      scfac(irfeag)     = sc1a  ; dscfacdt(irfeag)  = sc1adt  ; dscfacdd(irfeag)  = sc1add

      ratdum(irniga)    = ratraw(irniga) * sc1a
      dratdumdt(irniga) = dratrawdt(irniga)*sc1a + ratraw(irniga)*sc1adt
      dratdumdd(irniga) = dratrawdd(irniga)*sc1a + ratraw(irniga)*sc1add

      scfac(irniga)     = sc1a  ; dscfacdt(irniga)  = sc1adt  ; dscfacdd(irniga)  = sc1add


! fe52(a,p)co55
      ratdum(irfeap) = ratraw(irfeap) * sc1a
      dratdumdt(irfeap) = dratrawdt(irfeap)*sc1a + ratraw(irfeap)*sc1adt
      dratdumdd(irfeap) = dratrawdd(irfeap)*sc1a + ratraw(irfeap)*sc1add

      scfac(irfeap)     = sc1a  ; dscfacdt(irfeap)  = sc1adt  ; dscfacdd(irfeap)  = sc1add

      ratdum(ircopa)    = ratraw(ircopa) * sc1a
      dratdumdt(ircopa) = dratrawdt(ircopa)*sc1a + ratraw(ircopa)*sc1adt
      dratdumdd(ircopa) = dratrawdd(ircopa)*sc1a + ratraw(ircopa)*sc1add

      scfac(ircopa)     = sc1a  ; dscfacdt(ircopa)  = sc1adt  ; dscfacdd(ircopa)  = sc1add


      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   27.0d0,55.0d0,1.0d0,1.0d0, &
                   jscr,init,sc1a,sc1adt,sc1add)


! co55(p,g)ni56
      ratdum(ircopg)    = ratraw(ircopg) * sc1a
      dratdumdt(ircopg) = dratrawdt(ircopg)*sc1a + ratraw(ircopg)*sc1adt
      dratdumdd(ircopg) = dratrawdd(ircopg)*sc1a + ratraw(ircopg)*sc1add

      scfac(ircopg)     = sc1a  ; dscfacdt(ircopg)  = sc1adt  ; dscfacdd(ircopg)  = sc1add

      ratdum(irnigp)    = ratraw(irnigp) * sc1a
      dratdumdt(irnigp) = dratrawdt(irnigp)*sc1a + ratraw(irnigp)*sc1adt
      dratdumdd(irnigp) = dratrawdd(irnigp)*sc1a + ratraw(irnigp)*sc1add

      scfac(irnigp)     = sc1a  ; dscfacdt(irnigp)  = sc1adt  ; dscfacdd(irniga)  = sc1add



! reset the screen initialization flag
      init = 0

! end of screening block if
      end if


! now form those lovely dummy proton link rates

! mg24(a,p)27al(p,g)28si
      ratdum(irr1)     = 0.0d0  ; dratdumdt(irr1)  = 0.0d0  ; dratdumdd(irr1)  = 0.0d0
      denom    = ratdum(iralpa) + ratdum(iralpg)
      denomdt  = dratdumdt(iralpa) + dratdumdt(iralpg)
      denomdd  = dratdumdd(iralpa) + dratdumdd(iralpg)
      if (denom .gt. 1.0d-30) then
       zz = 1.0d0/denom
       ratdum(irr1)    = ratdum(iralpa)*zz
       dratdumdt(irr1) = (dratdumdt(iralpa) - ratdum(irr1)*denomdt)*zz
       dratdumdd(irr1) = (dratdumdd(iralpa) - ratdum(irr1)*denomdd)*zz
      end if

! si28(a,p)p31(p,g)s32
      ratdum(irs1)     = 0.0d0  ; dratdumdt(irs1)  = 0.0d0  ; dratdumdd(irs1)  = 0.0d0
      denom    = ratdum(irppa) + ratdum(irppg)
      denomdt  = dratdumdt(irppa) + dratdumdt(irppg)
      denomdd  = dratdumdd(irppa) + dratdumdd(irppg)
      if (denom .gt. 1.0d-30) then
       zz = 1.0d0/denom
       ratdum(irs1)    = ratdum(irppa)*zz
       dratdumdt(irs1) = (dratdumdt(irppa) - ratdum(irs1)*denomdt)*zz
       dratdumdd(irs1) = (dratdumdd(irppa) - ratdum(irs1)*denomdd)*zz
      end if

! s32(a,p)cl35(p,g)ar36
      ratdum(irt1)     = 0.0d0  ; dratdumdt(irt1)  = 0.0d0  ; dratdumdd(irt1)  = 0.0d0
      denom    = ratdum(irclpa) + ratdum(irclpg)
      denomdt  = dratdumdt(irclpa) + dratdumdt(irclpg)
      denomdd  = dratdumdd(irclpa) + dratdumdd(irclpg)
      if (denom .gt. 1.0d-30) then
       zz = 1.0d0/denom
       ratdum(irt1)    = ratdum(irclpa)*zz
       dratdumdt(irt1) = (dratdumdt(irclpa) - ratdum(irt1)*denomdt)*zz
       dratdumdd(irt1) = (dratdumdd(irclpa) - ratdum(irt1)*denomdd)*zz
      end if

! ar36(a,p)k39(p,g)ca40
      ratdum(iru1)     = 0.0d0  ; dratdumdt(iru1)  = 0.0d0  ; dratdumdd(iru1)  = 0.0d0
      denom    = ratdum(irkpa) + ratdum(irkpg)
      denomdt  = dratdumdt(irkpa) + dratdumdt(irkpg)
      denomdd  = dratdumdd(irkpa) + dratdumdd(irkpg)
      if (denom .gt. 1.0d-30) then
       zz   = 1.0d0/denom
       ratdum(iru1)   = ratdum(irkpa)*zz
       dratdumdt(iru1) = (dratdumdt(irkpa) - ratdum(iru1)*denomdt)*zz
       dratdumdd(iru1) = (dratdumdd(irkpa) - ratdum(iru1)*denomdd)*zz
      end if

! ca40(a,p)sc43(p,g)ti44
      ratdum(irv1)     = 0.0d0  ; dratdumdt(irv1)  = 0.0d0  ; dratdumdd(irv1)  = 0.0d0
      denom    = ratdum(irscpa) + ratdum(irscpg)
      denomdt  = dratdumdt(irscpa) + dratdumdt(irscpg)
      denomdd  = dratdumdd(irscpa) + dratdumdd(irscpg)
      if (denom .gt. 1.0d-30) then
       zz  = 1.0d0/denom
       ratdum(irv1)    = ratdum(irscpa)*zz
       dratdumdt(irv1) = (dratdumdt(irscpa) - ratdum(irv1)*denomdt)*zz
       dratdumdd(irv1) = (dratdumdd(irscpa) - ratdum(irv1)*denomdd)*zz
      end if

! ti44(a,p)v47(p,g)cr48
      ratdum(irw1)    = 0.0d0  ; dratdumdt(irw1) = 0.0d0  ; dratdumdd(irw1) = 0.0d0
      denom    = ratdum(irvpa) + ratdum(irvpg)
      denomdt  = dratdumdt(irvpa) + dratdumdt(irvpg)
      denomdd  = dratdumdd(irvpa) + dratdumdd(irvpg)
      if (denom .gt. 1.0d-30) then
       zz = 1.0d0/denom
       ratdum(irw1)    = ratdum(irvpa)*zz
       dratdumdt(irw1) = (dratdumdt(irvpa) - ratdum(irw1)*denomdt)*zz
       dratdumdd(irw1) = (dratdumdd(irvpa) - ratdum(irw1)*denomdd)*zz
      end if

! cr48(a,p)mn51(p,g)fe52
      ratdum(irx1)    = 0.0d0  ; dratdumdt(irx1) = 0.0d0  ; dratdumdd(irx1) = 0.0d0
      denom    = ratdum(irmnpa) + ratdum(irmnpg)
      denomdt  = dratdumdt(irmnpa) + dratdumdt(irmnpg)
      denomdd  = dratdumdd(irmnpa) + dratdumdd(irmnpg)
      if (denom .gt. 1.0d-30) then
       zz = 1.0d0/denom
       ratdum(irx1)    = ratdum(irmnpa)*zz
       dratdumdt(irx1) = (dratdumdt(irmnpa) - ratdum(irx1)*denomdt)*zz
       dratdumdd(irx1) = (dratdumdd(irmnpa) - ratdum(irx1)*denomdd)*zz
      endif

! fe52(a,p)co55(p,g)ni56
      ratdum(iry1)    = 0.0d0  ; dratdumdt(iry1) = 0.0d0  ; dratdumdd(iry1) = 0.0d0
      denom    = ratdum(ircopa) + ratdum(ircopg)
      denomdt  = dratdumdt(ircopa) + dratdumdt(ircopg)
      denomdd  = dratdumdd(ircopa) + dratdumdd(ircopg)
      if (denom .gt. 1.0d-30) then
       zz = 1.0d0/denom
       ratdum(iry1)    = ratdum(ircopa)*zz
       dratdumdt(iry1) = (dratdumdt(ircopa) - ratdum(iry1)*denomdt)*zz
       dratdumdd(iry1) = (dratdumdd(ircopa) - ratdum(iry1)*denomdd)*zz
      end if



! debugs
!      do i=1,nrat
!       if (ratdum(i) .lt. 0.0) then
!        write(6,110) i,ratnam(i),ratraw(i),scfac(i),ratdum(i)
! 110    format(1x,i4,' ',a,' ',1p3e12.4)
!        stop 'negative rate'
!       end if
!      enddo

!      write(6,*) btemp,bden
!      write(6,*) abar,zbar
!      do i=1,nrat
!       write(6,111) i,ratnam(i),ratraw(i),scfac(i),ratdum(i)
!       write(6,111) i,ratnam(i),ratdum(i),dratdumdt(i),dratdumdd(i)
! 111   format(1x,i4,' ',a,' ',1p3e14.6)
!      enddo
!      read(5,*)

      return
      end subroutine screen_aprox13







      subroutine init_aprox13
      include 'implno.dek'
      include 'const.dek'
      include 'network.dek'

! this routine initializes stuff for the aprox13 network
!
! declare
      integer           :: i
      real*8, parameter :: mev2erg = ev2erg*1.0d6, mev2gr  = mev2erg/clight**2


! zero all the isotope and rate pointers
      call init_isotope_rate_pointers


! set the size of the network and the number of rates
      idnet   = idaprox13

      ionmax  = 13
      iener   = ionmax + 1
      itemp   = ionmax + 2
      iden    = ionmax + 3
      ivelx   = ionmax + 4
      iposx   = ionmax + 5
      neqs    = iposx


      nrat    = 67
      netname = 'aprox13'


! set the id numbers of the elements

      ionbeg = 1
      ionend = ionmax

      ihe4  = ionbeg
      ic12  = ionbeg + 1
      io16  = ionbeg + 2
      ine20 = ionbeg + 3
      img24 = ionbeg + 4
      isi28 = ionbeg + 5
      is32  = ionbeg + 6
      iar36 = ionbeg + 7
      ica40 = ionbeg + 8
      iti44 = ionbeg + 9
      icr48 = ionbeg + 10
      ife52 = ionbeg + 11
      ini56 = ionbeg + 12


! set the names of the elements
      ionam(ihe4)  = 'he4  '
      ionam(ic12)  = 'c12  '
      ionam(io16)  = 'o16  '
      ionam(ine20) = 'ne20 '
      ionam(img24) = 'mg24 '
      ionam(isi28) = 'si28 '
      ionam(is32)  = 's32  '
      ionam(iar36) = 'ar36 '
      ionam(ica40) = 'ca40 '
      ionam(iti44) = 'ti44 '
      ionam(icr48) = 'cr48 '
      ionam(ife52) = 'fe52 '
      ionam(ini56) = 'ni56 '

      ionam(iener) = 'ener '
      ionam(itemp) = 'temp '
      ionam(iden)  = 'den  '
      ionam(ivelx) = 'velx '
      ionam(iposx) = 'posx '


! set the number of nucleons in the element
      aion(ihe4)  = 4.0d0
      aion(ic12)  = 12.0d0
      aion(io16)  = 16.0d0
      aion(ine20) = 20.0d0
      aion(img24) = 24.0d0
      aion(isi28) = 28.0d0
      aion(is32)  = 32.0d0
      aion(iar36) = 36.0d0
      aion(ica40) = 40.0d0
      aion(iti44) = 44.0d0
      aion(icr48) = 48.0d0
      aion(ife52) = 52.0d0
      aion(ini56) = 56.0d0

! set the number of protons in the element
      zion(ihe4)  = 2.0d0
      zion(ic12)  = 6.0d0
      zion(io16)  = 8.0d0
      zion(ine20) = 10.0d0
      zion(img24) = 12.0d0
      zion(isi28) = 14.0d0
      zion(is32)  = 16.0d0
      zion(iar36) = 18.0d0
      zion(ica40) = 20.0d0
      zion(iti44) = 22.0d0
      zion(icr48) = 24.0d0
      zion(ife52) = 26.0d0
      zion(ini56) = 28.0d0


! set the binding energy of the element
      bion(ihe4)  =  28.29603d0
      bion(ic12)  =  92.16294d0
      bion(io16)  = 127.62093d0
      bion(ine20) = 160.64788d0
      bion(img24) = 198.25790d0
      bion(isi28) = 236.53790d0
      bion(is32)  = 271.78250d0
      bion(iar36) = 306.72020d0
      bion(ica40) = 342.05680d0
      bion(iti44) = 375.47720d0
      bion(icr48) = 411.46900d0
      bion(ife52) = 447.70800d0
      bion(ini56) = 484.00300d0


! set the number of neutrons and mass
       nion(ionbeg:ionend) = aion(ionbeg:ionend) - zion(ionbeg:ionend)

! mass of each isotope
       mion(ionbeg:ionend) = nion(ionbeg:ionend)*mn + zion(ionbeg:ionend)*mp - bion(ionbeg:ionend)*mev2gr

! molar mass
       wion(ionbeg:ionend) = avo * mion(ionbeg:ionend)

! a common approximation
      wion(ionbeg:ionend) = aion(ionbeg:ionend)

! set the partition functions - statistical weights, ground-state only here
      wpart(ionbeg:ionend) = 1.0d0


! set the id numbers of the reaction rates
      ir3a   = 1
      irg3a  = 2
      ircag  = 3
      iroga  = 4
      ir1212 = 5
      ir1216 = 6
      ir1616 = 7
      iroag  = 8
      irnega = 9
      irneag = 10
      irmgga = 11
      irmgag = 12
      irsiga = 13
      irmgap = 14
      iralpa = 15
      iralpg = 16
      irsigp = 17
      irsiag = 18
      irsga  = 19
      irsiap = 20
      irppa  = 21
      irppg  = 22
      irsgp  = 23
      irsag  = 24
      irarga = 25
      irsap  = 26
      irclpa = 27
      irclpg = 28
      irargp = 29
      irarag = 30
      ircaga = 31
      irarap = 32
      irkpa  = 33
      irkpg  = 34
      ircagp = 35
      ircaag = 36
      irtiga = 37
      ircaap = 38
      irscpa = 39
      irscpg = 40
      irtigp = 41
      irtiag = 42
      ircrga = 43
      irtiap = 44
      irvpa  = 45
      irvpg  = 46
      ircrgp = 47
      ircrag = 48
      irfega = 49
      ircrap = 50
      irmnpa = 51
      irmnpg = 52
      irfegp = 53
      irfeag = 54
      irniga = 55
      irfeap = 56
      ircopa = 57
      ircopg = 58
      irnigp = 59

      irr1   = 60
      irs1   = 61
      irt1   = 62
      iru1   = 63
      irv1   = 64
      irw1   = 65
      irx1   = 66
      iry1   = 67


! set the names of the reaction rates
      ratnam(ir3a)   = 'r3a  '
      ratnam(irg3a)  = 'rg3a '
      ratnam(ircag)  = 'rcag '
      ratnam(ir1212) = 'r1212'
      ratnam(ir1216) = 'r1216'
      ratnam(ir1616) = 'r1616'
      ratnam(iroga)  = 'roga '
      ratnam(iroag)  = 'roag '
      ratnam(irnega) = 'rnega'
      ratnam(irneag) = 'rneag'
      ratnam(irmgga) = 'rmgga'
      ratnam(irmgag) = 'rmgag'
      ratnam(irsiga) = 'rsiga'
      ratnam(irmgap) = 'rmgap'
      ratnam(iralpa) = 'ralpa'
      ratnam(iralpg) = 'ralpg'
      ratnam(irsigp) = 'rsigp'
      ratnam(irsiag) = 'rsiag'
      ratnam(irsga)  = 'rsga '
      ratnam(irsiap) = 'rsiap'
      ratnam(irppa)  = 'rppa '
      ratnam(irppg)  = 'rppg '
      ratnam(irsgp)  = 'rsgp '
      ratnam(irsag)  = 'rsag '
      ratnam(irarga) = 'rarga'
      ratnam(irsap)  = 'rsap '
      ratnam(irclpa) = 'rclpa'
      ratnam(irclpg) = 'rclpg'
      ratnam(irargp) = 'rargp'
      ratnam(irarag) = 'rarag'
      ratnam(ircaga) = 'rcaga'
      ratnam(irarap) = 'rarap'
      ratnam(irkpa)  = 'rkpa '
      ratnam(irkpg)  = 'rkpg '
      ratnam(ircagp) = 'rcagp'
      ratnam(ircaag) = 'rcaag'
      ratnam(irtiga) = 'rtiga'
      ratnam(ircaap) = 'rcaap'
      ratnam(irscpa) = 'rscpa'
      ratnam(irscpg) = 'rscpg'
      ratnam(irtigp) = 'rtigp'
      ratnam(irtiag) = 'rtiag'
      ratnam(ircrga) = 'rcrga'
      ratnam(irtiap) = 'rtiap'
      ratnam(irvpa)  = 'rvpa '
      ratnam(irvpg)  = 'rvpg '
      ratnam(ircrgp) = 'rcrgp'
      ratnam(ircrag) = 'rcrag'
      ratnam(irfega) = 'rfega'
      ratnam(ircrap) = 'rcrap'
      ratnam(irmnpa) = 'rmnpa'
      ratnam(irmnpg) = 'rmnpg'
      ratnam(irfegp) = 'rfegp'
      ratnam(irfeag) = 'rfeag'
      ratnam(irniga) = 'rniga'
      ratnam(irfeap) = 'rfeap'
      ratnam(ircopa) = 'rcopa'
      ratnam(ircopg) = 'rcopg'
      ratnam(irnigp) = 'rnigp'

      ratnam(irr1)   = 'irr1 '
      ratnam(irs1)   = 'irs1 '
      ratnam(irt1)   = 'irt1 '
      ratnam(iru1)   = 'iru1 '
      ratnam(irv1)   = 'irv1 '
      ratnam(irw1)   = 'irw1 '
      ratnam(irx1)   = 'irx1 '
      ratnam(iry1)   = 'iry1 '


      return
      end subroutine init_aprox13
!---------------------------------------------------------------------







!---------------------------------------------------------------------
      subroutine netint(start,stptry,stpmin,stopp,stpsav,bc, &
                        eps,ylogi,nok,nbad,kount,odescal, &
                        derivs,jakob,bjakob,steper)
      include 'implno.dek'
      include 'burn_common.dek'
      include 'network.dek'



! ode integrator for stiff odes
! tuned for nnuclear reacton networks

! input:
! start    = beginning integration point
! stptry   = suggested first step size
! stpmin   = minimum allowable step size
! stopp    = ending integration point
! stpsav   = fixed time steps to save output ; if zero integration decides step
! bc       = initial conditions, array of physical dimension yphys
! eps      = desired fraction error during the integration
! odescal  = error scaling factor
! derivs   = name of the routine that contains the odes
! jakob    = name of the routine that contains the jacobian of the odes
! bjakob   = name of the routine that sets the pointers of the sparse jacobian
! steper   = name of the routine that will take a single step

! output:
! nok      = number of succesful steps taken
! nbad     = number of bad steps taken, bad but retried and then succesful
! kount    = total number of steps taken



! declare the pass
      external       derivs,jakob,bjakob,steper
      integer    ::  ylogi,nok,nbad,kount
      real*8     ::  start,stptry,stpmin,stopp,stpsav,bc(ylogi),eps,odescal


! local variables
      integer            :: i,nstp,idt
      integer, parameter :: nmax = abignet*nzmax, stpmax=200000
      real*8             :: yscal(nmax),y(nmax),dydx(nmax), &
                            tsav,cons,x,h,hdid,hnext,xx
      real*8, parameter  :: tiny=1.0d-15


! for timesteps to make pleasing plots 
      real*8             :: ylast(nmax), ratio,xfloor,ychangemax,ynew,yold,yy,dtx


! initialize
      if (ylogi  .gt. nmax) stop 'ylogi > nmax in routine netint'
      x      = start
      h      = sign(stptry,stopp-start)
      nok    = 0
      nbad   = 0
      kount  = 0
      idt    = 0


! store the first step
      y(1:ylogi) = bc(1:ylogi)
      tsav = x - 2.0d0 * stpsav


! take at most stpmax steps
      do nstp=1,stpmax

!      write(6,123) nstp,h
!123   format(1x,i6,1p3e24.16)


! store the old abundances
!       ylast(1:ionmax) = y(1:ionmax)


! positive definite abundance fractions
       y(1:ionmax) = min(1.0d0, max(y(1:ionmax),1.0d-30))


! renorm the abundances
       if (renorm_abundances .eq. 1) then
        zwork1(1:ionmax) = y(1:ionmax)*aion(1:ionmax)
        call indexx(ionmax,zwork1,izwork1)      
        xx = 0.0d0
        do i=1,ionmax
         xx = xx + zwork1(izwork1(i))
        end do
        xx = 1.0d0/xx
        zwork2(1:ionmax) = xx * zwork1(1:ionmax)
        y(1:ionmax) = min(1.0d0,max( zwork2(1:ionmax)/aion(1:ionmax),1.0d-30))
       end if



! get the right hand sides and form scaling vector
       call derivs(x,y,dydx)
       yscal(1:ylogi) = max(odescal,abs(y(1:ylogi)))


! store intermediate results
       if ( (abs(stpsav) - abs(x-tsav)) .le. tiny) then
        kount  = kount+1


! detailed file and print
        if (iprint_files .eq. 1) call net_output(kount,x,y,derivs)
        if (iprint_screen .eq. 1) then 
         call net_tty(kount,x,y,idt,hdid)
         tsav = x
        end if
       end if


! if the step can overshoot the stop point cut it
       if ((x+h-stopp)*(x+h-start) .gt. 0.0d0) h = stopp - x
       if ((stpsav .ne. 0.0d0) .and. (h .gt. (tsav - x + stpsav))) h = tsav + stpsav - x


! do an integration step

        call steper(y,dydx,ylogi,x,h,eps,yscal,hdid,hnext, &
                    derivs,jakob,bjakob,nstp,idt)


! tally if we did or did not do the step
       if (hdid.eq.h) then
        nok = nok+1
       else
        nbad = nbad+1
       end if



! this is the normal exit point
       if ( (nstp .eq. stpmax)                             .or. &
            (x-stopp)*(stopp-start).ge. 0.0d0              .or. &
            (psi .ge.  1.0 .and. y(itemp) .lt. temp_stop)  .or. &
            (psi .le. -1.0 .and. y(itemp) .gt. temp_stop)  .or. &
            (xmass_stop .gt. 0.0d0 .and. y(id_stop)*aion(id_stop) .lt. xmass_stop)   ) then

!        write(6,*) 'bailing'
!        write(6,*) (nstp .eq. stpmax)
!        write(6,*) nstp,stpmax
!        write(6,*)
!        write(6,*) (x-stopp)*(stopp-start).ge. 0.0d0
!        write(6,*) x,start,stopp
!        write(6,*)
!        write(6,*) (psi .ge.  1.0 .and. y(itemp) .lt. temp_stop) 
!        write(6,*) psi,y(itemp),temp_stop
!        write(6,*)
!        write(6,*) (psi .le. -1.0 .and. y(itemp) .gt. temp_stop) 
!        write(6,*) psi,y(itemp),temp_stop
!        write(6,*) 
!        write(6,*)  (y(id_stop)*aion(id_stop) .lt. xmass_stop)         
!        write(6,*) id_stop,y(id_stop),aion(id_stop),xmass_stop


! renorm the abundances
!        zwork1(1:ionmax) = y(1:ionmax)*aion(1:ionmax)
!        xx = 1.0d0/sum(zwork1(1:ionmax))
!        zwork2(1:ionmax) = xx * zwork1(1:ionmax)
!        y(1:ionmax) = min(1.0d0,max( zwork2(1:ionmax)/aion(1:ionmax),1.0d-30))


! save the final step
        bc(1:ylogi) = y(1:ylogi)
        kount = kount+1
        if (iprint_files .eq. 1) call net_output(kount,x,y,derivs)
        if (iprint_screen .eq. 1) call net_tty(kount,x,y,idt,hdid)
        return
       end if



! set the step size for the next iteration; stay above stpmin
!       dtx = 1.0e30
!        if  (kount .ne. 1) then
!        ratio      = 1.0d30
!        xfloor     = 1.0e-10
!        ychangemax = 0.05d0
!        do i=1,ionmax
!         ynew = max(y(i),1.0d-20)
!         yold = ylast(i)
!         yy   = abs(ynew - yold)
!         if (yy*ratio .gt. yold .and. yold .ge. xfloor) ratio=yold/yy
!        enddo
!        dtx = min(ratio*h*ychangemax,2.0d0*h)
!        end if
!
!       h = min(hnext,dtx)

! limit timestep changes to a factor of two
!       h = min(hnext,2.0d0*h)


! normal timestep choice
       h = hnext


! die 
       if (abs(h).lt.stpmin) then

        write(6,210) 'nstp=',nstp
        write(6,210) 'nok nbad',nok,nbad
        write(6,220) 'attempted time step =',stptry
        write(6,220) 'current time step =',h
        write(6,220) 'temperature y(itemp) =',y(itemp)
        write(6,220) 'temperature btemp    =',btemp
        write(6,220) 'density y(iden) =',y(iden)
        write(6,220) 'density bden    =',bden
        write(6,220) 'input composition:'
        write(6,230) (bc(i), i=1,ylogi)
        write(6,220) 'current composition:'
        write(6,230) (y(i), i=1,ylogi)
        stop 'h < stpmin in netint'

210     format(1x,a,4i6)
220     format(1x,a,1p3e24.16)
230     format(1x,1p3e24.16)

       end if


! back for another iteration or death
      enddo
      stop 'more than stpmax steps required in netint'
      end

!---------------------------------------------------------------------







!---------------------------------------------------------------------
      subroutine net_input(time_start,time_end,tin,din,vin,zin,ein,xin)
      include 'implno.dek'
      include 'const.dek'
      include 'vector_eos.dek'
      include 'burn_common.dek'
      include 'network.dek'


! declare the pass
      double precision time_start,time_end,tin,din,vin,zin,ein,xin(*)


! local variables
      character*80  :: string,word
      logical       :: err
      integer       :: i,j,k,ibtype,ictype,igues,kkase,ians,getnam
      double precision xneut,xh1,xhe4,xc12,xc13,xn14,xo16,xne20,xne22,xmg24,&
                       xsi28,xs32,xfe52,xfe54,xfe56,xni56,zye,xx,abar,zbar, &
                       wbar,xcess,ye,ye_orig,xmup,xmun,qdum,a,z,xelem, &
                       andgrev,value2

! bigbang specifics
      real*8            :: fac


! popular format statements
01    format(1x,a,a,a)
02    format(1x,a,'=',1pe10.3,' ',a,'=',1pe10.3,' ', &
                a,'=',1pe10.3,' ',a,'=',1pe10.3,' ', &
                a,'=',1pe10.3)
03    format(a)
04    format(1x,a,'=',i2,' ',a,'=',i2,' ', &
                a,'=',i2,' ',a,'=',i2,' ', &
                a,'=',i2)



! initialize the common block variables
      call net_initialize


! inititailize local variables
      ibtype     = 0
      ictype     = 0
      time_start = 0.0d0
      time_end   = 0.0d0
      bpres      = 0.0d0
      tin        = 0.0d0
      din        = 0.0d0
      vin        = 0.0d0
      zin        = 0.0d0
      zye        = 0.0d0
      xin(1:ionmax) = 1.0d-30


!---------------------------------------------------------------------------


! get the burn type
10    continue
      if (idnet .eq. idbigbang) then
       ibtype = 7
      else if (idnet .eq. idtorch) then

       write(6,01) 'give burning mode, <cr> for hydrostatic:'
       write(6,01) '     0 = stop'
       write(6,01) '     1 = onestep'
       write(6,01) '     2 = hydrostatic - constant temperature & density'
       write(6,01) '     3 = adiabatic expansion'
       write(6,01) '     4 = self-heat at constant density'
       write(6,01) '     5 = self heat at constant pressure'
       write(6,01) '     6 = self-heat pressure-temp trajectory'
       write(6,01) '     7 = big bang'
       write(6,01) '     8 = detonation (disabled)'
       write(6,01) '     9 = temp-den trajectory file'

       string = ' ' 
       read(5,03) string
       if (string(1:1) .eq. ' ') then
        ibtype = 2
       else
        ibtype = value2(string,err)
       if (err) goto 10
       if (ibtype .lt. 0 .or. ibtype .gt. 9) goto 10
       end if 

      else

       write(6,01) 'give burning mode, <cr> for hydrostatic:'
       write(6,01) '     0 = stop'
       write(6,01) '     1 = onestep'
       write(6,01) '     2 = hydrostatic - constant temperature & density'
       write(6,01) '     3 = adiabatic expansion'
       write(6,01) '     4 = self-heat at constant density'
       write(6,01) '     5 = self heat at constant pressure'
       write(6,01) '     9 = temp-den trajectory file'

       string = ' ' 
       read(5,03) string
       if (string(1:1) .eq. ' ') then
        ibtype = 2
       else
        ibtype = value2(string,err)
       if (err) goto 10
       if (ibtype .lt. 0 .or. ibtype .gt. 6) goto 10
       end if 
      end if


! set the burn type logical
      if (ibtype .eq. 0) then
       stop 'normal termination'
      else if (ibtype .eq. 1) then
       one_step = .true.
      else if (ibtype .eq. 2) then
       hydrostatic = .true.
      else if (ibtype .eq. 3) then
       expansion = .true.
      else if (ibtype .eq. 4) then
       self_heat_const_den = .true.
      else if (ibtype .eq. 5) then
       self_heat_const_pres = .true.
      else if (ibtype .eq. 6) then
       pt_hist = .true.
      else if (ibtype .eq. 7) then
       bbang = .true.
      else if (ibtype .eq. 8) then
       write(6,*)
       write(6,*) 'detonation mode disabled'
       write(6,*)
       goto 10
      else if (ibtype .eq. 9) then
       trho_hist = .true.
      else
       goto 10
      end if

!---------------------------------------------------------------------------



! general options
11    continue

      if (.not.bbang) then

       if (idnet .eq. idtorch) then

12      write(6,*)
        write(6,04) 'set general switches:'
        write(6,04) 'screening on',screen_on
        write(6,04) 'use rate tables',use_tables
        write(6,04) 'use weak reactions',weak_on
        write(6,04) 'use lmp weak rates',lmp_on
        write(6,04) 'pure network',pure_network
        write(6,04) 'print output files',iprint_files
        write(6,04) 'print to screen',iprint_screen
        write(6,*)
        write(6,01) '<cr> if these are ok, otherwise give 7 integer vector =>'

        string = ' ' 
        read(5,03) string
        if (string(1:1) .ne. ' ') then
         read(string,*) screen_on, use_tables, weak_on, lmp_on, pure_network, iprint_files, iprint_screen
         goto 12
        end if

       else

13      write(6,*)
        write(6,04) 'set general switches on=1 or off=0:'
        write(6,04) 'screening on',screen_on
        write(6,04) 'use rate tables',use_tables
        write(6,04) 'pure network',pure_network
        write(6,04) 'print output files',iprint_files
        write(6,04) 'print to screen',iprint_screen
        write(6,*)
        write(6,01) '<cr> if these are ok, otherwise enter 0 =>'

        string = ' ' 
        read(5,03) string
        if (string(1:1) .ne. ' ') then
         read(string,*) screen_on, use_tables, pure_network, iprint_files, iprint_screen
         goto 13
        end if
       end if

       if (screen_on .lt. 0 .or. screen_on .gt. 1) goto 11
       if (use_tables .lt. 0 .or. use_tables .gt. 1) goto 11
       if (weak_on .lt. 0 .or. weak_on .gt. 1) goto 11
       if (pure_network .lt. 0 .or. pure_network .gt. 1) goto 11
       if (iprint_files .lt. 0 .or. iprint_files .gt. 1) goto 11
       if (iprint_screen .lt. 0 .or. iprint_screen .gt. 1) goto 11
      end if

!---------------------------------------------------------------------------


! get the bigbang parameters; set default to wmap 2008 (5 year) values
      if (bbang) then

       eta1    = 6.23d-10
       xnnu    = 3.0d0
       hubble  = 70.5d0
       cmbtemp = 2.725d0

16     write(6,*)
       write(6,02) 'bigbang parameters - wmap 2005 values:'
       write(6,*)
       write(6,02) 'eta',eta1
       write(6,02) 'number of neutrino families',xnnu
       write(6,02) 'hubble constant',hubble
       write(6,02) 'present cmb temperature',cmbtemp
       write(6,*)
       write(6,01) '<cr> if ok, enter 0 to enter values, enter 1 to exit program =>'

       string = ' '
       read(5,03) string
       if (string(1:1) .ne. ' ') then
        i = value2(string,err)
        if (err) goto 16
        if (i .eq. 0) then
         write(6,01) 'give eta, xnu, hubble, and cmbtemp  =>'
         read(5,*) eta1, xnnu, hubble, cmbtemp
         goto 16
        else
         stop 'normal termination'
        end if
       end if
      end if


!---------------------------------------------------------------------------

! get an alternative the stopping condition; when the
! mass fraction of a given isotope falls below a given level


!      write(6,*)
!      string = ' ' 
!      write(6,*) 'stop when an isotope falls below a given abundance?', &
!                  ' <cr> = no, any other input = yes :'
!      read(5,03) string
!      if (string(1:1) .eq. ' ') then
       name_stop  = 'he4 '
       xmass_stop = -1.0d30
       id_stop    = ihe4
!      else 
!15     write(6,*) 'give the name of the isotope and the mass fraction'
!       write(6,*) 'for example: c12 0.50'
!       read(5,03) string
!       j = 1
!       i = getnam(string,word,j)
!       name_stop = word(1:5)
!       i = getnam(string,word,j)
!       xmass_stop = value2(word,err)
!       if (err) goto 15
!       write(6,*) name_stop,xmass_stop

! check that the name_stop isotope is in the network
!       do i=1,ionmax
!        if (ionam(i) .eq. name_stop) then
!         id_stop = i
!         goto 16
!        end if
!       enddo
!       write(6,*)
!       write(6,*) 'name_stop>',name_stop,'< not in network'
!       write(6,*)
!       goto 15
!       stop ' bad name for stopping isotope'
!      end if
! 16   continue


!---------------------------------------------------------------------------


! get the initial thermodynamics
      write(6,*)
      if (self_heat_const_pres) then
       write(6,01) 'give the ending time, temperature, pressure =>'
       read(5,*)  time_end,tin,bpres

! remove the initial temperature as a user-input; too many people enter inappropriate values
      else if (bbang) then
       write(6,01) 'give the ending time in seconds (suggest 1e6) =>'
       read(5,*)  time_end
       tin= 1.0d11

! thermodynamic profile being given
      else if (trho_hist .or. pt_hist) then
       write(6,*) 'give the trajectory file =>'
       read(5,03) trho_file

       call update2(time_start,tin,din)
       time_start    = zwork1(1)
       time_end      = zwork1(2)
       zye           = zwork1(3)

      else 
       write(6,01) 'give the ending time, temperature, density =>'
       read(5,*)  time_end,tin,din
      end if

! limit the temperature since the rates are invalid much above t9=100
       tin = min(1.0d11,tin)



! set variables based on the burn type

! adiabatic expansion
! psi =  1 is an adiabatic expansion, -1 in an adiabatic implosion

      if (expansion) then
       psi       = 1.0d0
!       psi       = -1.0d0
       den0      = din
       temp0     = tin
       temp_stop = 1.0d7
!       temp_stop = 1.0d10
       if ( (psi .ge. 1.0  .and. temp_stop .ge. tin)  .or. &
            (psi .le. -1.0 .and. temp_stop .le. tin)) &
          stop 'bad adiabatic temp_stop in routine burner'
      end if

!---------------------------------------------------------------------------



! set the initial composition
      if (bbang) then

! set the initial n and p abundances; equation 3 of wagoner et al 1967
       fac = exp((mn - mp)*clight**2/(kerg*tin))
       xneut = 1.0d0/(1.0d0 + fac)
       xh1   = 1.0d0 - xneut

! set the density from the temperature and eta1
       f1  = 30.0d0 * zeta3/pi**4 * asol/(kerg*avo)
       din = f1 * eta1 * tin**3

      else
 20    write(6,*)
       write(6,01) 'give initial composition, <cr> for specify initial composition :'
       write(6,01) '     ictype = 0 = read from file'
       write(6,01) '              1 = solar abundances'
       write(6,01) '              2 = nse'
       write(6,01) '              3 = specify initial composition at terminal'

       string = ' ' 
       read(5,03) string
       if (string(1:1) .eq. ' ') then
        ictype = 3
       else
        ictype = value2(string,err)
        if (err) goto 20
        if (ictype .lt. 0 .or. ictype .gt. 3) goto 20
       end if 
      end if



      if (ictype .eq. 0) then
       call read_abundances_from_file(xin)
      
      else if (ictype .eq. 1) then
       do i=1,ionmax
        xin(i) = andgrev(ionam(i),z,a,xelem)
       enddo
       if (iprot .ne. 0) xin(iprot) = andgrev('h1   ',z,a,xelem)

      else if (ictype .eq. 2) then
       if (zye .eq. 0.0) zye   = 0.5d0
       igues = 1
       call nse(tin,din,zye,igues,1,1,xin,xmun,xmup,0)

      else if (ictype .eq. 3) then
        write(6,01) &
        'n h1 he4 c12 c13 n14 o16 ne20 ne22 mg24 si28 s32 fe52 fe54 fe56 ni56 =>'
        read(5,*) xneut,xh1,xhe4,xc12,xc13,xn14,xo16,xne20,xne22,xmg24,xsi28, &
                  xs32,xfe52,xfe54,xfe56,xni56
      end if


! set the composition variables
      if (ictype .eq. 3 .or. bbang) then
       if (ineut .ne. 0) xin(ineut) = xneut
       if (ih1   .ne. 0) xin(ih1)   = xh1
       if (iprot .ne. 0) xin(iprot) = xh1
       if (ih1 .ne. 0 .and. iprot .ne. 0) xin(iprot) = 0.0d0
       if (ihe4  .ne. 0) xin(ihe4)  = xhe4
       if (ic12  .ne. 0) xin(ic12)  = xc12
       if (ic13  .ne. 0) xin(ic13)  = xc13
       if (in14  .ne. 0) xin(in14)  = xn14
       if (io16  .ne. 0) xin(io16)  = xo16
       if (ine20 .ne. 0) xin(ine20) = xne20
       if (ine22 .ne. 0) xin(ine22) = xne22
       if (img24 .ne. 0) xin(img24) = xmg24
       if (isi28 .ne. 0) xin(isi28) = xsi28
       if (is32  .ne. 0) xin(is32)  = xs32
       if (ife52 .ne. 0) xin(ife52) = xfe52
       if (ife54 .ne. 0) xin(ife54) = xfe54
       if (ife56 .ne. 0) xin(ife56) = xfe56
       if (ini56 .ne. 0) xin(ini56) = xni56
      end if


! write out the input composition so far
!      write(6,02) (ionam(i),xin(i), i=1,ionmax)
!      read(5,*)


! normalize the composition
      xin(1:ionmax) = min(1.0d0,max(xin(1:ionmax),1.0d-30))
      xx = 1.0d0/sum(xin(1:ionmax))
      xin(1:ionmax) = min(1.0d0,max(xin(1:ionmax) * xx,1.0d-30))

!      write(6,*) 'post norm', ionmax,xin(ih1)
!      write(6,02) (ionam(i),xin(i), i=1,ionmax)
!      read(5,*)


!---------------------------------------------------------------------------


! get the thermodynamic state of the final mixture

        call azbarapprox13(xin,aion,zion,wion,ionmax, &
             zwork1,abar,zbar,wbar,ye_orig,xcess)


! get the thermodynamic state
      temp_row(1) = tin
      den_row(1)  = din
      ptot_row(1) = bpres
      abar_row(1) = abar
      zbar_row(1) = zbar
      jlo_eos = 1
      jhi_eos = 1

      if (self_heat_const_pres .or. pt_hist) then
       den_row(1)  = bpres * abar/(avo * kerg * tin)
       call invert_helm_pt
       din = den_row(1)
      else
!       call helmeos
         call helmeosaprox13
       bpres = ptot_row(1)
      endif
      ein   = etot_row(1)



!---------------------------------------------------------------------------

! get the output root file name
      write(6,*)  ' '
      write(6,01) 'give output root name, <cr> for default "foo_"=>'
      read(5,03) hfile
      if (hfile(1:2) .eq. '  ')  hfile = 'foo_'





!---------------------------------------------------------------------------


! write out the final input
        write(6,*)
        write(6,02) 'time_start',time_start,'time_end',time_end
        write(6,02) 'tin',tin,'din',din,'bpres',bpres,'ein',ein

! largest mass fractions
        call indexx(ionmax,xin,izwork1)
        j = min(20,ionmax)
        k = max(ionmax-19,1)
        write(6,*) j,' largest mass fractions'
        do i=ionmax,k,-1
         if (xin(izwork1(i)) .gt. 1.0e-12) &
            write(6,02) ionam(izwork1(i)),xin(izwork1(i))
        end do

! nonconservation, abar, zbar of the mixture
        xx = sum(xin(1:ionmax)) 
        write(6,02) '1-sum',1.0d0 - xx
        write(6,02) 'abar',abar,'zbar',zbar,'ye',zbar/abar
        write(6,*)


!---------------------------------------------------------------------------


! there is probably a better place for this
! if requested, adjust the number of equations being solved

      if (pure_network .eq. 1) then
       neqs  = ionmax
       btemp = tin
       bden  = din
      end if


! say how many isotopes and rates 
       write(6,*)
       write(6,222) 'ionbeg=',ionbeg,' ionend=',ionend,' ionmax=',ionmax,' neqs=',neqs,' nrates=',nrat 
222    format(1x,5(a,i5))
       write(6,*)

      return
      end
!---------------------------------------------------------------------







!---------------------------------------------------------------------
! this routine contains auxillary network routine

! routines for a tree construction to mark nonzero matrix locations
! routine screen6 computes screening factors
! routine screen5 computes screening factors
! routine snupp computes neutrino loss rates for the pp chain
! routine snucno computes neutrino loss rates for the cno cycles
! routine sneut5 computes neutrino loss rates
! routine ifermi12 does an inverse fermi integral of order 1/2
! routine zfermim12 does an inverse fermi integral of order -1/2

! routine ecapnuc02 computes electron capture rates
! routine ecapnuc computes electron capture rates
! routine mazurek computes ni56 electron capture rates
! routine time_scales computes various timescales
! routine ener_gener_rate computes the instantaneous energy generation rate

! routine read_partition_functions does what it says on the rauscher & thielemann data
! routine get_partition_functions interpolates the partition functions

! routine init_isotope_rate_pointers initializes integer name & reaction rate pointers for the networks
! routine prepare_for_net_eval does some prep work before evaluating the right-hand-side and jacobian evaluations




! interfaces to a balanced tree sort
      subroutine tree_init(n)
      implicit none
      common/locatdat/nmax
      integer nmax,n
      nmax=n+1
      call avlinit(30*nmax+2048)
      return
      end


      subroutine tree(i,j,eloc,neloc,nterm,nzo,iloc,jloc,np)
      implicit none
      common/locatdat/nmax
      integer nmax
      integer i,j,neloc,eloc(neloc),nterm,nzo,np,iloc(np),jloc(np)
      integer iat,nzo_old

      nzo_old = nzo
      nterm   = nterm + 1
      if (nterm .gt. neloc) then
       write(6,10) 'nterm=',nterm
       write(6,10) 'neloc=',neloc
 10    format(1x,a,' ',i6)
       stop 'nterm > neloc in routine tree'
      end if

      call avlinsert(i*nmax+j,iat,nzo)
      eloc(nterm) = iat

      if (nzo .gt. np) stop 'nzo > np in routine tree'
      if (nzo .ne. nzo_old) then
       eloc(nterm) = nzo
       iloc(nzo) = i
       jloc(nzo) = j
      end if
      return
      end



      subroutine tree_out(irow,icol,nzo,np)
      implicit none
      common/locatdat/nmax
      integer nmax,np,nzo,i,irow(np),icol(np)
      call avlgetlist(np,icol,nzo)
      call avlfree()
      do i=1,nzo
         irow(i)=icol(i)/nmax
         icol(i)=icol(i)-irow(i)*nmax
      enddo
      return
      end





!    AVL sort
!   
!    In 1960 two Russian mathematicians, Georgii Maksimovich
!    Adel'son-Vel'skii and Evgenii Mikhailovich Landis developed a
!    technique for keeping a binary search tree balanced as items are
!    inserted into it.  called AVL trees.
!   
!    efficiently sort integers in N log N operations
!   
!    implemetation taken from
!    http://www.moorpark.cc.ca.us/~csalazar/cs20/nonlin.txt (buggy)
!    see also
!    http://swww.ee.uwa.edu.au/~plsd210/ds/AVL.html
!    http://www.purists.org/georg/avltree/ (my favorite site)
!   
!    implemented by Alexander Heger, 20010129
!    avldelete   by Alexander Heger, 20010205

!=======================================================================
!=======================================================================

      MODULE avl_data
      implicit none
!      integer, parameter :: maxavldata = 65536
      integer :: maxavldata
      integer, parameter :: maxavlindex = 4
      integer, parameter :: NULL = 0
      integer, parameter :: l_LEFT = 1
      integer, parameter :: l_RIGHT = l_LEFT+1 ! do not change
      integer, parameter :: l_BAL = 3
      integer, parameter :: l_KEY = 4
      integer, parameter :: i_ROOT = 1
      integer, parameter :: i_NODEOFFSET = 1
      integer, parameter :: l_ROOT = l_RIGHT
      integer, parameter :: l_RIGHTHEAVY = 1 ! do not change
      integer, parameter :: l_BALANCED = 0   ! do not change
      integer, parameter :: l_LEFTHEAVY = -l_RIGHTHEAVY
      integer, parameter :: l_UNBALANCED = l_BALANCED+1
      integer, parameter :: l_GARBAGE = l_LEFT
!    tree data
      integer :: maxel
      integer :: garbage
!      integer, dimension(maxavlindex,maxavldata) :: ichild
      integer, allocatable, dimension(:,:) :: ichild
      SAVE
      END MODULE avl_data

!=======================================================================

      MODULE avl_stack
      implicit none
      integer, parameter :: maxdepth = 48
      integer, parameter :: i_STACKBASE = 1
      integer, dimension(maxdepth) :: istack,lrstack
      integer :: ipstack
      END MODULE avl_stack

!=======================================================================

      subroutine avlinit(nmax)
      USE avl_data
      implicit none
      integer, intent(IN) :: nmax

      SELECT CASE (nmax)
      CASE (1:)
         maxavldata=nmax+1
      CASE (0)
         maxavldata=1024
      END SELECT
      IF (nmax >= 0) THEN
         CALL avlfree()
         ALLOCATE(ichild(maxavlindex,maxavldata))
      ENDIF

!    initialize root pointer and zero number of elements
      ichild(l_ROOT,i_ROOT)=NULL
      maxel=0
!    initialize garbage list
      garbage=NULL

      end

!=======================================================================

      subroutine avlfree()
      USE avl_data
      implicit none
      IF (ALLOCATED(ichild)) DEALLOCATE(ichild)
      end

!=======================================================================

      subroutine avlgetlist(nmax,list,n)
      USE avl_data
      USE avl_stack
      implicit none

!    some constants
      integer, intent(IN) :: nmax
      integer, dimension(nmax), intent(OUT) :: list
      integer, intent(OUT) :: n

!    running variables
      integer :: i, lr, ii

      n=0

      i=ichild(l_ROOT,i_ROOT)
      IF (i == NULL) RETURN

      IF (nmax < maxel) THEN
         WRITE(*,"(' [AVL LIST] ERROR: list too small for data.')")
         n=-1
         RETURN
      ENDIF

!    recursively traverse tree to get sorted list of key values
      ipstack=i_STACKBASE-1
      lr=l_LEFT

      DO
         IF (lr <= l_LEFT) THEN
!    add left branch
            ii=ichild(l_LEFT,i)
            IF (ii /= NULL) THEN
               ipstack=ipstack+1
               istack(ipstack)=i
               lrstack(ipstack)=l_RIGHT
               i=ii
               lr=l_LEFT
               CYCLE
            ENDIF
         ENDIF
         IF (lr <= l_RIGHT) THEN
!    add node
            n=n+1
            list(n)=ichild(l_KEY,i)
!    add right branch
            ii=ichild(l_RIGHT,i)
            IF (ii /= NULL) THEN
               ipstack=ipstack+1
               istack(ipstack)=i
               lrstack(ipstack)=l_RIGHT+1
               i=ii
               lr=l_LEFT
               CYCLE
            ENDIF
         ENDIF

         IF (ipstack < i_STACKBASE) EXIT
         i=istack(ipstack)
         lr=lrstack(ipstack)
         ipstack=ipstack-1
      ENDDO

      IF (n /= maxel) THEN
         WRITE(*,"(' [AVL LIST] ERROR in AVL data.')")
         n=-1
         RETURN
      ENDIF

      END

!=======================================================================

      subroutine avltree()
      USE avl_data
      USE avl_stack
      implicit none

!    some constants
      character*(*), parameter :: form = "(I5)"
      integer, parameter :: nwidth = 5
      integer, parameter :: nmax = 1024

!    running variables
      integer, dimension(nmax) :: level, index
      character*(nwidth*nmax),dimension(5,maxdepth+1) :: line
      character*(nwidth) :: item
      integer :: i, lr, ii
      integer :: n, maxlevel
      integer :: l1,l2

      IF (nmax < maxel) THEN
         WRITE(*,"(' [AVL TREE] ERROR: too much data.')")
         RETURN
      ENDIF

      n=0
      maxlevel=0

      i=ichild(l_ROOT,i_ROOT)
      IF (i == NULL) RETURN

!    recursively traverse tree to get sorted list of key values
      ipstack=i_STACKBASE-1
      lr=l_LEFT

      DO
         IF (lr <= l_LEFT) THEN
!    add left branch
            ii=ichild(l_LEFT,i)
            IF (ii /= NULL) THEN
               ipstack=ipstack+1
               istack(ipstack)=i
               lrstack(ipstack)=l_RIGHT
               i=ii
               lr=l_LEFT
               CYCLE
            ENDIF
         ENDIF
         IF (lr <= l_RIGHT) THEN
!    add node
            n=n+1
            index(n)=i
            level(n)=ipstack+1
            maxlevel=MAX(maxlevel,ipstack+1)
!    add right branch
            ii=ichild(l_RIGHT,i)
            IF (ii /= NULL) THEN
               ipstack=ipstack+1
               istack(ipstack)=i
               lrstack(ipstack)=l_RIGHT+1
               i=ii
               lr=l_LEFT
               CYCLE
            ENDIF
         ENDIF

         IF (ipstack < i_STACKBASE) EXIT
         i=istack(ipstack)
         lr=lrstack(ipstack)
         ipstack=ipstack-1
      ENDDO

      line=" "
      DO i=1,n
         l1=1+nwidth*(i-1)
         l2=nwidth*i
         write(item,form) index(i)
         line(1,level(i))(l1:l2)=item
         write(item,form) ichild(l_KEY,index(i))
         line(2,level(i))(l1:l2)=item
         write(item,form) ichild(l_BAL,index(i))
         line(3,level(i))(l1:l2)=item
         write(item,form) ichild(l_LEFT,index(i))
         line(4,level(i))(l1:l2)=item
         write(item,form) ichild(l_RIGHT,index(i))
         line(5,level(i))(l1:l2)=item
         line(1,maxlevel+1)(l1:l2)='------------'
      ENDDO

      WRITE(*,"(A)") line(1,maxlevel+1)(1:nwidth*n)//"|"
      DO ii=1,maxlevel
         DO i=1,5
            WRITE(*,"(A)") line(i,ii)(1:nwidth*n)//"|"
         ENDDO
      ENDDO
      WRITE(*,"(A)") line(1,maxlevel+1)(1:nwidth*n)//"|"

      END

!=======================================================================

      subroutine avlincrease
      USE avl_data
      implicit none

      integer, allocatable, dimension(:,:) :: ichild_temp
!      integer, dimension(:,:) :: ichild_temp
      ALLOCATE(ichild_temp(maxavlindex,maxavldata))
      ichild_temp(:,:)=ichild(:,:)
      DEALLOCATE(ichild)
      maxavldata=maxavldata*2
      ALLOCATE(ichild(maxavlindex,maxavldata))
      ichild(:,:)=ichild_temp(:,:)
      DEALLOCATE(ichild_temp)
      WRITE (*,"(' [AVL INCREASE] INFO: now ',I8,' elements.')") &
           maxavldata
      end

!=======================================================================

      subroutine avlinsert(key,ijk,nzo)
      USE avl_data
      USE avl_stack

      implicit none

      integer, intent(IN) :: key
      integer, intent(OUT):: ijk,nzo

      integer :: i, ii, lr, irevbal, ic, ip, lrs, lri

      ijk = 1

      i=i_root
      lr=l_ROOT
      ii=ichild(lr,i)
      ipstack=i_STACKBASE
      istack(ipstack)=i
      lrstack(ipstack)=lr

!    find location and insert
      DO WHILE (ii /= NULL)
         i=ii

         ijk = i
!         write(6,*) ijk

         SELECT CASE (key - ichild(l_KEY,i))
         CASE (0)
!            write(6,*) 'same as ',ijk-1
            ijk = ijk - 1
            RETURN              ! element already present: done.
         CASE (:-1)
            lr=l_LEFT
         CASE DEFAULT
            lr=l_RIGHT
         END SELECT
         ipstack=ipstack+1
         istack(ipstack)=i
         lrstack(ipstack)=lr
         ii=ichild(lr,i)

      ENDDO

!    initialize new element
      maxel=maxel+1

      nzo = maxel

      IF (garbage /= NULL) THEN
         ii=garbage
         garbage=ichild(l_GARBAGE,garbage)
      ELSE
         IF (maxel == maxavldata-1) CALL avlincrease
         ii=maxel+i_NODEOFFSET
      ENDIF
      ichild(lr,i)=ii
      ichild(l_KEY,ii)=key
      ichild(l_BAL,ii)=l_BALANCED
      ichild(l_LEFT:l_RIGHT,ii)=NULL

!    balance tree
      irevbal=l_UNBALANCED
      DO WHILE ((ipstack > i_STACKBASE) .AND. (irevbal /= l_BALANCED))
         i=istack(ipstack)
         lr=lrstack(ipstack)
         ipstack=ipstack-1
         lri=(l_LEFT+l_RIGHT)-lr ! pointer to the opposite direction
                                ! sign for balance determination
         lrs=2*lr-(l_LEFT+l_RIGHT)
         SELECT CASE (ichild(l_BAL,i)*lrs)
         CASE (l_LEFTHEAVY)
            ichild(l_BAL,i)=l_BALANCED
            irevbal=l_BALANCED
         CASE (l_BALANCED)
            ichild(l_BAL,i)=lrs
         CASE DEFAULT
!    update tree
            ic=ichild(lr,i)
            IF (ichild(l_BAL,ic) == lrs) THEN
!    single rotation
               ichild(l_BAL,i)=l_BALANCED
               ichild(l_BAL,ic)=l_BALANCED
               ichild(lr,i)=ichild(lri,ic)
               ichild(lri,ic)=i
               ichild(lrstack(ipstack),istack(ipstack))=ic
            ELSE IF (ichild(l_BAL,ic) == -lrs) THEN
!    double rotation
               ip=ichild(lri,ic)
               SELECT CASE (ichild(l_BAL,ip)*lrs)
               CASE (l_LEFTHEAVY)
                  ichild(l_BAL,i)=l_BALANCED
                  ichild(l_BAL,ic)=lrs
               CASE (l_BALANCED)
                  ichild(l_BAL,i)=l_BALANCED
                  ichild(l_BAL,ic)=l_BALANCED
               CASE DEFAULT
                  ichild(l_BAL,i)=-lrs
                  ichild(l_BAL,ic)=l_BALANCED
               END SELECT
               ichild(l_BAL,ip)=l_BALANCED
               ichild(lri,ic)=ichild(lr,ip)
               ichild(lr,ip)=ic
               ichild(lr,i)=ichild(lri,ip)
               ichild(lri,ip)=i
               ichild(lrstack(ipstack),istack(ipstack))=ip
            ENDIF
            irevbal=l_BALANCED
         END SELECT
      ENDDO

      END
!=======================================================================

      subroutine avldelete(key)
      USE avl_data
      USE avl_stack

      implicit none

      integer, intent(IN) :: key

      integer :: i, ii, lr, irevbal, ic, ip, lrs, lri, i0, ipstack0

      lr=l_ROOT
      ipstack=i_STACKBASE
      istack(ipstack)=i_ROOT
      lrstack(ipstack)=l_ROOT
      i=ichild(l_ROOT,i_ROOT)

!    find location to delete
      DO
         IF (i == NULL) RETURN ! element not present
         SELECT CASE (key - ichild(l_KEY,i))
         CASE (0)
            EXIT              ! element found
         CASE (:-1)
            lr=l_LEFT
         CASE DEFAULT
            lr=l_RIGHT
         END SELECT
         ipstack=ipstack+1
         istack(ipstack)=i
         lrstack(ipstack)=lr
         i=ichild(lr,i)
      ENDDO
      i0=i
!    find closest element to replace it
!    decide whether to take left or right branch
      SELECT CASE (ichild(l_BAL,i))
      CASE (l_LEFTHEAVY)
         lri=l_LEFT
      CASE (l_RIGHTHEAVY)
         lri=l_RIGHT
      CASE DEFAULT
         lri=(l_LEFT+l_RIGHT)-lr
      END SELECT

!    now search for it
      ii=ichild(lri,i)
      IF (ii /= NULL) THEN
         ipstack0=ipstack
!    go one step in lrx direction
         ipstack=ipstack+1
         istack(ipstack)=i
         lrstack(ipstack)=lri
!    now seach element most in opposite direction
         lr=(l_LEFT+l_RIGHT)-lri
         i=ii
         ii=ichild(lr,i)
         DO WHILE (ii /= NULL)
            ipstack=ipstack+1
            istack(ipstack)=i
            lrstack(ipstack)=lr
            i=ii
            ii=ichild(lr,i)
         ENDDO
!    found element to swap
!    do swap
         ic=ichild(lri,i)
         ichild(lri,i)=ichild(lri,i0)
         ichild(lr,i)=ichild(lr,i0)
         ichild(l_BAL,i)=ichild(l_BAL,i0)
         ichild(lrstack(ipstack0),istack(ipstack0))=i
!    CORRECT STACK
         istack(ipstack0+1)=i
         ichild(lrstack(ipstack),istack(ipstack))=ic
      ELSE
!    element last of chain
!    just remove it
         ic=NULL
      ENDIF
!    move rest of branch one level up
      ichild(lrstack(ipstack),istack(ipstack))=ic

!    (i): balance=balance - lrs
!    start regular re-balancing loop
      irevbal=l_UNBALANCED

      DO WHILE ((ipstack > i_STACKBASE) .AND. (irevbal /= l_BALANCED))
         i=istack(ipstack)
         lr=lrstack(ipstack)
         ipstack=ipstack-1
!   
         lri=(l_LEFT+l_RIGHT)-lr
         lrs=2*lr-(l_LEFT+l_RIGHT)

         SELECT CASE (ichild(l_BAL,i)*lrs)
         CASE (l_RIGHTHEAVY)
            ichild(l_BAL,i)=l_BALANCED
         CASE (l_BALANCED)
            ichild(l_BAL,i)=-lrs
            irevbal=l_BALANCED
         CASE DEFAULT
!    update tree
            ic=ichild(lri,i)
            IF (ichild(l_BAL,ic) == lrs) THEN
!    double rotation
               ip=ichild(lr,ic)
               SELECT CASE (ichild(l_BAL,ip)*lrs)
               CASE (l_RIGHTHEAVY)
                  ichild(l_BAL,i)=l_BALANCED
                  ichild(l_BAL,ic)=-lrs
               CASE (l_LEFTHEAVY)
                  ichild(l_BAL,i)=lrs
                  ichild(l_BAL,ic)=l_BALANCED
               CASE DEFAULT
                  ichild(l_BAL,i)=l_BALANCED
                  ichild(l_BAL,ic)=l_BALANCED
               END SELECT
               ichild(l_BAL,ip)=l_BALANCED
               ichild(lri,i)=ichild(lr,ip)
               ichild(lr,ip)=i
               ichild(lr,ic)=ichild(lri,ip)
               ichild(lri,ip)=ic
               ichild(lrstack(ipstack),istack(ipstack))=ip
            ELSE
!    single rotation
               IF (ichild(l_BAL,ic) == l_BALANCED) THEN
                  ichild(l_BAL,i)=-lrs
                  ichild(l_BAL,ic)=lrs
                  irevbal=l_BALANCED
               ELSE
                  ichild(l_BAL,i)=l_BALANCED
                  ichild(l_BAL,ic)=l_BALANCED
               ENDIF
               ichild(lri,i)=ichild(lr,ic)
               ichild(lr,ic)=i
               ichild(lrstack(ipstack),istack(ipstack))=ic
            ENDIF
         END SELECT
      ENDDO

!    free element
      ichild(l_GARBAGE,i0)=garbage
      garbage=i0
      maxel=maxel-1

      END










      subroutine screen6(jscreen, &
                         temp,den,zbar,abar,z2bar, &
                         z1,a1,z2,a2, &
                         scor,scordt,scordd)
      include 'implno.dek'
      include 'const.dek'
      include 'network.dek'

! this subroutine calculates screening factors and their derivatives
! for nuclear reaction rates in the weak, intermediate and strong regimes.
! based on graboske, dewit, grossman and cooper apj 181 457 1973 for
! weak screening. based on alastuey and jancovici apj 226 1034 1978,
! with plasma parameters from itoh et al apj 234 1079 1979, for strong
! screening.

! vector version

! input:
! jscreen = length of vector
! temp    = temperature
! den     = density
! zbar    = mean charge per nucleus
! abar    = mean number of nucleons per nucleus
! z2bar   = mean square charge per nucleus
! z1 a1   = charge and number in the entrance channel
! z2 a2   = charge and number in the exit channel

! output:
! scor    = screening correction
! scordt  = derivative of screening correction with temperature
! scordd  = derivative of screening correction with density


! declare the pass
      integer          jscreen
      double precision temp,den,zbar,abar,z2bar, &
                       z1(jscreen),a1(jscreen),z2(jscreen),a2(jscreen), &
                       scor(jscreen),scordt(jscreen),scordd(jscreen)


! local variables
      integer          i
      double precision bb,cc,dccdt,dccdd, &
                       pp,dppdt,dppdd,qq,dqqdt,dqqdd,rr,drrdt,drrdd, &
                       ss,dssdt,dssdd,tt,dttdt,dttdd,uu,duudt,duudd, &
                       vv,dvvdt,dvvdd,a3,da3,tempi,dtempi,deni, &
                       qlam0z,qlam0zdt,qlam0zdd, &
                       h12w,dh12wdt,dh12wdd,h12,dh12dt,dh12dd, &
                       h12x,dh12xdt,dh12xdd,alfa,beta, &
                       taufac,taufacdt,gamp,gampdt,gampdd, &
                       gamef,gamefdt,gamefdd, &
                       tau12,tau12dt,alph12,alph12dt,alph12dd, &
                       xlgfac,dxlgfacdt,dxlgfacdd, &
                       gamp14,gamp14dt,gamp14dd, &
                       xni,dxnidd,ytot


! screening variables
! zs13    = (z1+z2)**(1./3.)
! zhat    = combination of z1 and z2 raised to the 5/3 power
! zhat2   = combination of z1 and z2 raised to the 5/12 power
! lzav    = log of effective charge
! aznut   = combination of a1,z1,a2,z2 raised to 1/3 power


      integer          nscreen_max
      parameter        (nscreen_max = 2*abignet + 40)


      double precision zs13(nscreen_max),zhat(nscreen_max), &
                       zhat2(nscreen_max),lzav(nscreen_max), &
                       aznut(nscreen_max),zs13inv(nscreen_max), &
                       fac1(nscreen_max),fac2(nscreen_max), &
                       h12_vec(nscreen_max), &
                       dh12dt_vec(nscreen_max), &
                       dh12dd_vec(nscreen_max)


! parameter fact is the cube root of 2
      double precision  x13,x14,x53,x532,x512,fact,co2,gamefx,gamefs, &
                        blend_frac
      parameter        (x13        = 1.0d0/3.0d0, &
                        x14        = 1.0d0/4.0d0, &
                        x53        = 5.0d0/3.0d0, &
                        x532       = 5.0d0/32.0d0, &
                        x512       = 5.0d0/12.0d0, &
                        fact       = 1.25992104989487d0, &
                        co2        = x13 * 4.248719d3, &
                        gamefx     = 0.3d0, &
                        gamefs     = 0.8d0, &
                        blend_frac = 0.05d0)


! compute and store the more expensive screening factors
      if (init_screen .eq. 1) then
       init_screen = 0
       write(6,*) 'begin store of screening factors'

       if (jscreen .gt. nscreen_max) &
       stop 'jscreen > nscreen_max in screen6'

       do i=1,jscreen
        zs13(i)    = (z1(i) + z2(i))**x13
        zs13inv(i) = 1.0d0/zs13(i)
        zhat(i)    = (z1(i) + z2(i))**x53  - z1(i)**x53 - z2(i)**x53
        zhat2(i)   = (z1(i) + z2(i))**x512 - z1(i)**x512 -z2(i)**x512
        lzav(i)    = x53 * log(z1(i)*z2(i)/(z1(i) + z2(i)))
        aznut(i)   = (z1(i)**2*z2(i)**2*a1(i)*a2(i)/(a1(i)+a2(i)))**x13
        fac1(i)    = 0.896434d0 * zhat(i)
        fac2(i)    = 3.44740d0  * zhat2(i)
       enddo
       write(6,*) 'done store of screening factors'
      endif


! calculate average plasma
      ytot     = 1.0d0/abar
      rr       = den * ytot
      tempi   = 1.0d0/temp
      dtempi  = -tempi*tempi
      deni    = 1.0d0/den

      pp       = sqrt(rr*tempi*(z2bar + zbar))
      qq       = 0.5d0/pp *(z2bar + zbar)
      dppdt    = qq*rr*dtempi
      dppdd    = qq*ytot*tempi

      qlam0z   = 1.88d8 * tempi * pp
      qlam0zdt = 1.88d8 * (dtempi*pp + tempi*dppdt)
      qlam0zdd = 1.88d8 * tempi * dppdd

      taufac   = co2 * tempi**x13
      taufacdt = -x13*taufac*tempi

      qq      = rr * zbar
      xni     = qq**x13
      dxnidd  = x13 * xni * deni

      gamp    = 2.27493d5 * tempi * xni
      gampdt  = 2.27493d5 * dtempi * xni
      gampdd  = 2.27493d5 * tempi * dxnidd



! calculate individual screening factors, start the pipeline
      do i=1,jscreen

       bb       = z1(i) * z2(i)
       qq       = fact * bb * zs13inv(i)
       gamef    = qq * gamp
       gamefdt  = qq * gampdt
       gamefdd  = qq * gampdd

       tau12    = taufac * aznut(i)
       tau12dt  = taufacdt * aznut(i)

       qq       = 1.0d0/tau12
       alph12   = gamef * qq
       alph12dt = (gamefdt - alph12*tau12dt) * qq
       alph12dd = gamefdd * qq



! limit alph12 to 1.6 to prevent unphysical behavior.
! this should really be replaced by a pycnonuclear reaction rate formula
       if (alph12 .gt. 1.6) then
        alph12   = 1.6d0
        alph12dt = 0.0d0
        alph12dd = 0.0d0

        gamef    = 1.6d0 * tau12
        gamefdt  = 1.6d0 * tau12dt
        gamefdd  = 0.0d0

        qq       = zs13(i)/(fact * bb)
        gamp     = gamef * qq
        gampdt   = gamefdt * qq
        gampdd   = 0.0d0
       end if



! weak screening regime
       h12w    = bb * qlam0z
       dh12wdt = bb * qlam0zdt
       dh12wdd = bb * qlam0zdd

       h12     = h12w
       dh12dt  = dh12wdt
       dh12dd  = dh12wdd



! intermediate and strong sceening regime
       if (gamef .gt. gamefx) then

        qq       = sqrt(gamp)
        gamp14   = sqrt(qq)
        rr       = 1.0d0/gamp
        qq       = 0.25d0*gamp14*rr
        gamp14dt = qq * gampdt
        gamp14dd = qq * gampdd

        cc       =   gamp * fac1(i) - gamp14 * fac2(i) &
                   - 0.5551d0   * (log(gamp) + lzav(i)) &
                   - 2.996d0

        dccdt    =   gampdt * fac1(i)  - gamp14dt * fac2(i) &
                   - 0.5551d0*rr*gampdt

        dccdd    =   gampdd * fac1(i) - gamp14dd * fac2(i) &
                   - 0.5551d0*rr*gampdd

        qq     = alph12 * alph12
        a3     = qq * alph12
        da3    = 3.0d0 * qq

        qq     = 0.014d0 + 0.0128d0*alph12
        dqqdt  = 0.0128d0*alph12dt
        dqqdd  = 0.0128d0*alph12dd

        rr     = x532 - alph12*qq
        drrdt  = -(alph12dt*qq + alph12*dqqdt)
        drrdd  = -(alph12dd*qq + alph12*dqqdd)

        ss     = tau12*rr
        dssdt  = tau12dt*rr + tau12*drrdt
        dssdd  = tau12*drrdd

        tt     =  -0.0098d0 + 0.0048d0*alph12
        dttdt  = 0.0048d0*alph12dt
        dttdd  = 0.0048d0*alph12dd

        uu     =  0.0055d0 + alph12*tt
        duudt  = alph12dt*tt + alph12*dttdt
        duudd  = alph12dd*tt + alph12*dttdd

        vv   = gamef * alph12 * uu
        dvvdt= gamefdt*alph12*uu + gamef*(alph12dt*uu + alph12*duudt)
        dvvdd= gamefdd*alph12*uu + gamef*(alph12dd*uu + alph12*duudd)

        h12     = cc - a3 * (ss + vv)
        rr      = da3 * (ss + vv)
        dh12dt  = dccdt - rr*alph12dt - a3*(dssdt + dvvdt)
        dh12dd  = dccdd - rr*alph12dd - a3*(dssdd + dvvdd)

        rr     =  1.0d0 - 0.0562d0*a3
        ss     =  -0.0562d0*da3
        drrdt  = ss*alph12dt
        drrdd  = ss*alph12dd

        if (rr .ge. 0.77d0) then
         xlgfac    = rr
         dxlgfacdt = drrdt
         dxlgfacdd = drrdd
        else
         xlgfac    = 0.77d0
         dxlgfacdt = 0.0d0
         dxlgfacdd = 0.0d0
        end if


        h12    = log(xlgfac) + h12
        rr     = 1.0d0/xlgfac
        dh12dt = rr*dxlgfacdt + dh12dt
        dh12dd = rr*dxlgfacdd + dh12dd


        if (gamef .le. gamefs) then
         rr     =  2.0d0*(gamefs - gamef)
         drrdt  = -2.0d0*gamefdt
         drrdd  = -2.0d0*gamefdd

         ss     = 2.0d0*(gamef-gamefx)
         dssdt  = 2.0d0*gamefdt
         dssdd  = 2.0d0*gamefdd


! store current values for possible blending
         h12x    = h12
         dh12xdt = dh12dt
         dh12xdd = dh12dd

         vv     = h12
         h12    = h12w*rr + vv*ss
         dh12dt = dh12wdt*rr + h12w*drrdt + dh12dt*ss + vv*dssdt
         dh12dd = dh12wdd*rr + h12w*drrdd + dh12dd*ss + vv*dssdd


! blend the transition region - from bill paxton
        if (gamefs - gamef .lt. blend_frac*(gamefs - gamefx)) then
          alfa   = (gamefs - gamef) / (blend_frac*(gamefs - gamefx))
          alfa   = 0.5d0 * (1d0 - cos(pi*alfa))
          beta   = 1.0d0 - alfa
          h12    = alfa * h12 + beta * h12x
          dh12dt = alfa * dh12dt + beta * dh12xdt
          dh12dd = alfa * dh12dd + beta * dh12xdd
         end if
        end if

! end of intermediate and strong screening if
       end if


! store what we got
       h12           = max(min(h12,300.0d0),0.0d0)
       if (h12 .eq. 300.0d0) then
        dh12dt = 0.0d0
        dh12dt = 0.0d0
       end if

       scor(i)   = exp(h12)
       scordt(i) = scor(i) * dh12dt
       scordd(i) = scor(i) * dh12dd


!       h12_vec(i)    = h12
!       dh12dt_vec(i) = dh12dt
!       dh12dt_vec(i) = dh12dd

!       scor(i)   = exp(h12_vec(i))
!       scordt(i) = scor(i) * dh12dt_vec(i)
!       scordd(i) = scor(i) * dh12dd_vec(i)


! end of individual screening pipeline
      end do



! crank the exponential
!       scor(1:jscreen)   = exp(h12_vec(1:jscreen))
!       scordt(1:jscreen) = scor(1:jscreen) * dh12dt_vec(1:jscreen)
!       scordd(1:jscreen) = scor(1:jscreen) * dh12dd_vec(1:jscreen)

      return
      end







      subroutine screen5(temp,den,zbar,abar,z2bar, &
                         z1,a1,z2,a2,jscreen,init, &
                         scor,scordt,scordd)
      include 'implno.dek'
      include 'const.dek'
      include 'network.dek'

! this subroutine calculates screening factors and their derivatives
! for nuclear reaction rates in the weak, intermediate and strong regimes.
! based on graboske, dewit, grossman and cooper apj 181 457 1973 for
! weak screening. based on alastuey and jancovici apj 226 1034 1978,
! with plasma parameters from itoh et al apj 234 1079 1979, for strong
! screening.

! input:
! temp    = temperature
! den     = density
! zbar    = mean charge per nucleus
! abar    = mean number of nucleons per nucleus
! z2bar   = mean square charge per nucleus
! z1 a1   = charge and number in the entrance channel
! z2 a2   = charge and number in the exit channel
! jscreen = counter of which reaction is being calculated
! init    = flag to compute the more expensive functions just once

! output:
! scor    = screening correction
! scordt  = derivative of screening correction with temperature
! scordd  = derivative of screening correction with density


! declare the pass
      integer          jscreen,init
      double precision temp,den,zbar,abar,z2bar,z1,a1,z2,a2, &
                       scor,scordt,scordd


! local variables
      double precision aa,daadt,daadd,bb,cc,dccdt,dccdd, &
                       pp,dppdt,dppdd,qq,dqqdt,dqqdd,rr,drrdt,drrdd, &
                       ss,dssdt,dssdd,tt,dttdt,dttdd,uu,duudt,duudd, &
                       vv,dvvdt,dvvdd,a3,da3,tempi,dtempi,deni, &
                       qlam0z,qlam0zdt,qlam0zdd, &
                       h12w,dh12wdt,dh12wdd,h12,dh12dt,dh12dd, &
                       h12x,dh12xdt,dh12xdd,alfa,beta, &
                       taufac,taufacdt,gamp,gampdt,gampdd, &
                       gamef,gamefdt,gamefdd, &
                       tau12,tau12dt,alph12,alph12dt,alph12dd, &
                       xlgfac,dxlgfacdt,dxlgfacdd, &
                       gamp14,gamp14dt,gamp14dd, &
                       xni,dxnidd,ytot, &
                       temp_old,den_old,zbar_old,abar_old


! screening variables
! zs13    = (z1+z2)**(1./3.)
! zhat    = combination of z1 and z2 raised to the 5/3 power
! zhat2   = combination of z1 and z2 raised to the 5/12 power
! lzav    = log of effective charge
! aznut   = combination of a1,z1,a2,z2 raised to 1/3 power


      integer          nscreen_max
      parameter        (nscreen_max = 2*abignet + 40)

      double precision zs13(nscreen_max),zhat(nscreen_max), &
                       zhat2(nscreen_max),lzav(nscreen_max), &
                       aznut(nscreen_max),zs13inv(nscreen_max)


! parameter fact is the cube root of 2
      double precision  x13,x14,x53,x532,x512,fact,co2,gamefx,gamefs, &
                        blend_frac
      parameter        (x13    = 1.0d0/3.0d0, &
                        x14    = 1.0d0/4.0d0, &
                        x53    = 5.0d0/3.0d0, &
                        x532   = 5.0d0/32.0d0, &
                        x512   = 5.0d0/12.0d0, &
                        fact   = 1.25992104989487d0, &
                        co2    = x13 * 4.248719d3, &
                        gamefx = 0.3d0, &
                        gamefs = 0.8d0, &
                        blend_frac = 0.05d0)


      data     temp_old/-1.0d0/, den_old/-1.0d0/, &
               zbar_old/-1.0d0/, abar_old/-1.0d0/




! compute and store the more expensive screening factors
      if (init .eq. 1) then
       if (jscreen .gt. nscreen_max) &
       stop 'jscreen > nscreen_max in screen5'
       zs13(jscreen)    = (z1 + z2)**x13
       zs13inv(jscreen) = 1.0d0/zs13(jscreen)
       zhat(jscreen)    = (z1 + z2)**x53  - z1**x53 - z2**x53
       zhat2(jscreen)   = (z1 + z2)**x512 - z1**x512 -z2**x512
       lzav(jscreen)    = x53 * log(z1*z2/(z1 + z2))
       aznut(jscreen)   = (z1**2 * z2**2 * a1*a2 / (a1 + a2))**x13
      endif


! calculate average plasma, if need be
      if (temp_old .ne. temp .or. &
          den_old  .ne. den  .or. &
          zbar_old  .ne. zbar  .or. &
          abar_old  .ne. abar ) then

       temp_old = temp
       den_old  = den
       zbar_old  = zbar
       abar_old  = abar

       ytot     = 1.0d0/abar
       rr       = den * ytot
       tempi   = 1.0d0/temp
       dtempi  = -tempi*tempi
       deni    = 1.0d0/den

       pp       = sqrt(rr*tempi*(z2bar + zbar))
       qq       = 0.5d0/pp *(z2bar + zbar)
       dppdt    = qq*rr*dtempi
       dppdd    = qq*ytot*tempi

       qlam0z   = 1.88d8 * tempi * pp
       qlam0zdt = 1.88d8 * (dtempi*pp + tempi*dppdt)
       qlam0zdd = 1.88d8 * tempi * dppdd

       taufac   = co2 * tempi**x13
       taufacdt = -x13*taufac*tempi

       qq      = rr*zbar
       xni     = qq**x13
       dxnidd  = x13 * xni * deni

       aa     = 2.27493d5 * tempi * xni
       daadt  = 2.27493d5 * dtempi * xni
       daadd  = 2.27493d5 * tempi * dxnidd
      end if


! calculate individual screening factors
      bb       = z1 * z2
      gamp     = aa
      gampdt   = daadt
      gampdd   = daadd

      qq       = fact * bb * zs13inv(jscreen)
      gamef    = qq * gamp
      gamefdt  = qq * gampdt
      gamefdd  = qq * gampdd

      tau12    = taufac * aznut(jscreen)
      tau12dt  = taufacdt * aznut(jscreen)

      qq       = 1.0d0/tau12
      alph12   = gamef * qq
      alph12dt = (gamefdt - alph12*tau12dt) * qq
      alph12dd = gamefdd * qq



! limit alph12 to 1.6 to prevent unphysical behavior.
! this should really be replaced by a pycnonuclear reaction rate formula
      if (alph12 .gt. 1.6) then
       alph12   = 1.6d0
       alph12dt = 0.0d0
       alph12dd = 0.0d0

       gamef    = 1.6d0 * tau12
       gamefdt  = 1.6d0 * tau12dt
       gamefdd  = 0.0d0

       qq       = zs13(jscreen)/(fact * bb)
       gamp     = gamef * qq
       gampdt   = gamefdt * qq
       gampdd   = 0.0d0
      end if



! weak screening regime
      h12w    = bb * qlam0z
      dh12wdt = bb * qlam0zdt
      dh12wdd = bb * qlam0zdd

      h12     = h12w
      dh12dt  = dh12wdt
      dh12dd  = dh12wdd



! intermediate and strong sceening regime
      if (gamef .gt. gamefx) then

       gamp14   = gamp**x14
       rr       = 1.0d0/gamp
       qq       = 0.25d0*gamp14*rr
       gamp14dt = qq * gampdt
       gamp14dd = qq * gampdd

       cc       =   0.896434d0 * gamp * zhat(jscreen) &
                  - 3.44740d0  * gamp14 * zhat2(jscreen) &
                  - 0.5551d0   * (log(gamp) + lzav(jscreen)) &
                  - 2.996d0

       dccdt    =   0.896434d0 * gampdt * zhat(jscreen) &
                  - 3.44740d0  * gamp14dt * zhat2(jscreen) &
                  - 0.5551d0*rr*gampdt

       dccdd    =   0.896434d0 * gampdd * zhat(jscreen) &
                  - 3.44740d0  * gamp14dd * zhat2(jscreen) &
                  - 0.5551d0*rr*gampdd

       a3     = alph12 * alph12 * alph12
       da3    = 3.0d0 * alph12 * alph12

       qq     = 0.014d0 + 0.0128d0*alph12
       dqqdt  = 0.0128d0*alph12dt
       dqqdd  = 0.0128d0*alph12dd

       rr     = x532 - alph12*qq
       drrdt  = -(alph12dt*qq + alph12*dqqdt)
       drrdd  = -(alph12dd*qq + alph12*dqqdd)

       ss     = tau12*rr
       dssdt  = tau12dt*rr + tau12*drrdt
       dssdd  = tau12*drrdd

       tt     =  -0.0098d0 + 0.0048d0*alph12
       dttdt  = 0.0048d0*alph12dt
       dttdd  = 0.0048d0*alph12dd

       uu     =  0.0055d0 + alph12*tt
       duudt  = alph12dt*tt + alph12*dttdt
       duudd  = alph12dd*tt + alph12*dttdd

       vv   = gamef * alph12 * uu
       dvvdt= gamefdt*alph12*uu + gamef*alph12dt*uu + gamef*alph12*duudt
       dvvdd= gamefdd*alph12*uu + gamef*alph12dd*uu + gamef*alph12*duudd

       h12     = cc - a3 * (ss + vv)
       rr      = da3 * (ss + vv)
       dh12dt  = dccdt - rr*alph12dt - a3*(dssdt + dvvdt)
       dh12dd  = dccdd - rr*alph12dd - a3*(dssdd + dvvdd)

       rr     =  1.0d0 - 0.0562d0*a3
       ss     =  -0.0562d0*da3
       drrdt  = ss*alph12dt
       drrdd  = ss*alph12dd

       if (rr .ge. 0.77d0) then
        xlgfac    = rr
        dxlgfacdt = drrdt
        dxlgfacdd = drrdd
       else
        xlgfac    = 0.77d0
        dxlgfacdt = 0.0d0
        dxlgfacdd = 0.0d0
       end if


       h12    = log(xlgfac) + h12
       rr     = 1.0d0/xlgfac
       dh12dt = rr*dxlgfacdt + dh12dt
       dh12dd = rr*dxlgfacdd + dh12dd


       if (gamef .le. gamefs) then
        rr     =  2.0d0*(gamefs - gamef)
        drrdt  = -2.0d0*gamefdt
        drrdd  = -2.0d0*gamefdd

        ss     = 2.0d0*(gamef - gamefx)
        dssdt  = 2.0d0*gamefdt
        dssdd  = 2.0d0*gamefdd


! store current values for possible blending
        h12x    = h12
        dh12xdt = dh12dt
        dh12xdd = dh12dd

        vv     = h12
        h12    = h12w*rr + vv*ss
        dh12dt = dh12wdt*rr + h12w*drrdt + dh12dt*ss + vv*dssdt
        dh12dd = dh12wdd*rr + h12w*drrdd + dh12dd*ss + vv*dssdd

! blend the transition region - from bill paxton
       if (gamefs - gamef .lt. blend_frac*(gamefs - gamefx)) then
         alfa   = (gamefs - gamef) / (blend_frac*(gamefs - gamefx))
         alfa   = 0.5d0 * (1d0 - cos(pi*alfa))
         beta   = 1.0d0 - alfa
         h12    = alfa * h12 + beta * h12x
         dh12dt = alfa * dh12dt + beta * dh12xdt
         dh12dd = alfa * dh12dd + beta * dh12xdd
        end if
       end if


! end of intermediate and strong screening if
      end if


! machine limit the output
      h12    = max(min(h12,300.0d0),0.0d0)
      scor   = exp(h12)
      if (h12 .eq. 300.0d0) then
       scordt = 0.0d0
       scordd = 0.0d0
      else
       scordt = scor * dh12dt
       scordd = scor * dh12dd
      end if

!      write(6,111) 'weak =',h12w,' total =',h12,
!     1             ' 1-ratio =',1.0d0-h12w/h12,' correction',scor
! 111  format(1x,4(a,1pe13.6))
!      read(5,*)

      return
      end









      double precision function snupp(yp,ratepp,ybe7,ratebeec, &
                                      yb8,rateb8epnu)
      include 'implno.dek'
      include 'const.dek'

! computes approximate neutrino losses from pp chain reactions
! see page 142 of astro 289j notes for these loss formulas

! input:
! yp         = proton molar abbundance
! ratepp     = pp reaction rate
! ybe7       = be7 molar abundance
! ratebeec   = be7 electron capture reaction rate
! yb8        = b8 molar abundance
! rateb8epnu = b8 decay reaction rate


! declare the pass
      double precision yp,ratepp,ybe7,ratebeec,yb8,rateb8epnu


! local variables
      double precision pp1nu,pp2nu,pp3nu,conv
      parameter        (conv = ev2erg*1.0d6*avo)


! nu losses from p(p,e-nu)h2
      pp1nu  = yp*yp*ratepp * 0.5d0 * 0.263d0


! nu losses from be7(n=>p)li7
      pp2nu  = ybe7 * ratebeec * 0.81d0


! nu losses from b8(p=>n)be8=>2a
      pp3nu  = yb8 * rateb8epnu * 7.73d0

! sum the pp-chain neutrino losses and convert to erg/g/s
      snupp  = (pp1nu + pp2nu + pp3nu) * conv

      return
      end



      double precision function snucno(yn13,bc13,bn13,yo14,bn14,bo14, &
                                       yo15,bn15,bo15,yf17,bo17,bf17, &
                                       yf18,bo18,bf18)
      include 'implno.dek'
      include 'const.dek'

! computes approximate neutrino losses from cno cycle  reactions
! see page 142 of astro 289j notes for these loss formulas

! input:
! yn13 = n13 molar abundance
! bc13 = c13 binding energy in mev
! bn13 = n13 binding energy in mev
! yo14 = o14 molar abundance
! bn14 = n14 binding energy in mev
! bo14 = o14 binding energy in mev
! yo15 = o15 molar abundance
! bn15 = n15 binding energy in mev
! bo15 = o15 binding energy in mev
! yf17 = f17 molar abundance
! bo17 = o17 binding energy in mev
! bf17 = f17 binding energy in mev
! yf18 = f18 molar abundance
! bo18 = o18 binding energy in mev
! bf18 = f18 binding energy in mev


! declare the pass
      double precision yn13,bc13,bn13,yo14,bn14,bo14, &
                       yo15,bn15,bo15,yf17,bo17,bf17, &
                       yf18,bo18,bf18

! local variables
      double precision sum,sum2,enu13n,enu14o,enu15o,enu17f,enu18f, &
                       conv,lntwo,tm1,tm2,tm3,tm4,tm5
      parameter        (conv  = ev2erg*1.0d6*avo, &
                        lntwo = 0.693147181d0, &
                        tm1   = lntwo/597.9d0, &
                        tm2   = lntwo/70.606d0, &
                        tm3   = lntwo/124.0, &
                        tm4   = lntwo/64.49, &
                        tm5   = lntwo/6586.2)


! 13n(e+nu)13c
      sum    = bc13 - bn13 - 0.782d0 - 1.022d0
      sum    = 1.0d0 + sum/0.511d0
      sum2   = sum*sum
      enu13n = 0.5d0 * sum * 0.511d0 * (1.0d0 - 1.0d0/sum2) &
               * (1.0d0 - 1.0d0/(4.0d0*sum) - 1.0d0/(9.0d0*sum2))
      enu13n = yn13 * enu13n * tm1



! hot cno cycle 14o(e+nu)14n
      sum    = bn14 - bo14 - 0.782d0 - 1.022d0
      sum    = 1.0d0 + sum/0.511d0
      sum2   = sum*sum
      enu14o = 0.5d0 * sum * 0.511d0 * (1.0d0 - 1.0d0/sum2) &
               * (1.0d0 - 1.0d0/(4.0d0*sum) - 1.0d0/(9.0d0*sum2))
      enu14o = yo14 * enu14o * tm2


! 15o(e+nu)15n
      sum    = bn15 - bo15 - 0.782d0 - 1.022d0
      sum    = 1.0d0 + sum/0.511d0
      sum2   = sum*sum
      enu15o = 0.5d0 * sum * 0.511d0 * (1.0d0 - 1.0d0/sum2) &
               * (1.0d0 - 1.0d0/(4.0d0*sum) - 1.0d0/(9.0d0*sum2))
      enu15o = yo15 * enu15o * tm3


! 17f(e+nu)17o
      sum    = bo17 - bf17 - 0.782d0 - 1.022d0
      sum    = 1.0d0 + sum/0.511d0
      sum2   = sum*sum
      enu17f = 0.5d0 * sum * 0.511d0 * (1.0d0 - 1.0d0/sum2) &
               * (1.0d0 - 1.0d0/(4.0d0*sum) - 1.0d0/(9.0d0*sum2))
      enu17f  = yf17 * enu17f * tm4


! 18f(e+nu)18o
      sum    = bo18 - bf18 - 0.782d0 - 1.022d0
      sum    = 1.0d0 + sum/0.511d0
      sum2   = sum*sum
      enu18f = 0.5d0 * sum * 0.511d0 * (1.0d0 - 1.0d0/sum2) &
               * (1.0d0 - 1.0d0/(4.0d0*sum) - 1.0d0/(9.0d0*sum2))
      enu18f = yf18 * enu18f * tm5


! sum the cno cycle losses and convert to erg/g/s
      snucno = (enu13n + enu14o + enu15o + enu17f + enu18f) * conv

      return
      end








      subroutine sneut5(temp,den,abar,zbar, &
                        snu,dsnudt,dsnudd,dsnuda,dsnudz)
      include 'implno.dek'
      include 'const.dek'

! this routine computes neutrino losses from the analytic fits of
! itoh et al. apjs 102, 411, 1996, and also returns their derivatives.

! input:
! temp = temperature
! den  = density
! abar = mean atomic weight
! zbar = mean charge

! output:
! snu    = total neutrino loss rate in erg/g/sec
! dsnudt = derivative of snu with temperature
! dsnudd = derivative of snu with density
! dsnuda = derivative of snu with abar
! dsnudz = derivative of snu with zbar


! declare the pass
      double precision temp,den,abar,zbar, &
                       snu,dsnudt,dsnudd,dsnuda,dsnudz

! local variables
      double precision spair,spairdt,spairdd,spairda,spairdz, &
                       splas,splasdt,splasdd,splasda,splasdz, &
                       sphot,sphotdt,sphotdd,sphotda,sphotdz, &
                       sbrem,sbremdt,sbremdd,sbremda,sbremdz, &
                       sreco,srecodt,srecodd,srecoda,srecodz

      double precision t9,xl,xldt,xlp5,xl2,xl3,xl4,xl5,xl6,xl7,xl8,xl9, &
                       xlmp5,xlm1,xlm2,xlm3,xlm4,xlnt,cc,den6,tfermi, &
                       a0,a1,a2,a3,b1,b2,c00,c01,c02,c03,c04,c05,c06, &
                       c10,c11,c12,c13,c14,c15,c16,c20,c21,c22,c23,c24, &
                       c25,c26,dd00,dd01,dd02,dd03,dd04,dd05,dd11,dd12, &
                       dd13,dd14,dd15,dd21,dd22,dd23,dd24,dd25,b,c,d,f0, &
                       f1,deni,tempi,abari,zbari,f2,f3,z,xmue,ye, &
                       dum,dumdt,dumdd,dumda,dumdz, &
                       gum,gumdt,gumdd,gumda,gumdz


! pair production
      double precision rm,rmdd,rmda,rmdz,rmi,gl,gldt, &
                       zeta,zetadt,zetadd,zetada,zetadz,zeta2,zeta3, &
                       xnum,xnumdt,xnumdd,xnumda,xnumdz, &
                       xden,xdendt,xdendd,xdenda,xdendz, &
                       fpair,fpairdt,fpairdd,fpairda,fpairdz, &
                       qpair,qpairdt,qpairdd,qpairda,qpairdz

! plasma
      double precision gl2,gl2dt,gl2dd,gl2da,gl2dz,gl12,gl32,gl72,gl6, &
                       ft,ftdt,ftdd,ftda,ftdz,fl,fldt,fldd,flda,fldz, &
                       fxy,fxydt,fxydd,fxyda,fxydz

! photo
      double precision tau,taudt,cos1,cos2,cos3,cos4,cos5,sin1,sin2, &
                       sin3,sin4,sin5,last,xast, &
                       fphot,fphotdt,fphotdd,fphotda,fphotdz, &
                       qphot,qphotdt,qphotdd,qphotda,qphotdz

! brem
      double precision t8,t812,t832,t82,t83,t85,t86,t8m1,t8m2,t8m3,t8m5, &
                       t8m6, &
                       eta,etadt,etadd,etada,etadz,etam1,etam2,etam3, &
                       fbrem,fbremdt,fbremdd,fbremda,fbremdz, &
                       gbrem,gbremdt,gbremdd,gbremda,gbremdz, &
                       u,gm1,gm2,gm13,gm23,gm43,gm53,v,w,fb,gt,gb, &
                       fliq,fliqdt,fliqdd,fliqda,fliqdz, &
                       gliq,gliqdt,gliqdd,gliqda,gliqdz

! recomb
      double precision ifermi12,zfermim12,nu,nudt,nudd,nuda,nudz, &
                       nu2,nu3,bigj,bigjdt,bigjdd,bigjda,bigjdz



! numerical constants
      double precision fac1,fac2,fac3,oneth,twoth,con1,sixth,iln10
      parameter        (fac1   = 5.0d0 * pi / 3.0d0, &
                        fac2   = 10.0d0 * pi, &
                        fac3   = pi / 5.0d0, &
                        oneth  = 1.0d0/3.0d0, &
                        twoth  = 2.0d0/3.0d0, &
                        con1   = 1.0d0/5.9302d0, &
                        sixth  = 1.0d0/6.0d0, &
                        iln10  = 4.342944819032518d-1)


! theta is sin**2(theta_weinberg) = 0.2319 plus/minus 0.00005 (1996)
! xnufam is the number of neutrino flavors = 3.02 plus/minus 0.005 (1998)
! change theta and xnufam if need be, and the changes will automatically
! propagate through the routine. cv and ca are the vector and axial currents.

      double precision theta,xnufam,cv,ca,cvp,cap,tfac1,tfac2,tfac3, &
                       tfac4,tfac5,tfac6
      parameter        (theta  = 0.2319d0, &
                        xnufam = 3.0d0, &
                        cv     = 0.5d0 + 2.0d0 * theta, &
                        cvp    = 1.0d0 - cv, &
                        ca     = 0.5d0, &
                        cap    = 1.0d0 - ca, &
                        tfac1  = cv*cv + ca*ca + &
                                 (xnufam-1.0d0) * (cvp*cvp+cap*cap), &
                        tfac2  = cv*cv - ca*ca + &
                                 (xnufam-1.0d0) * (cvp*cvp - cap*cap), &
                        tfac3  = tfac2/tfac1, &
                        tfac4  = 0.5d0 * tfac1, &
                        tfac5  = 0.5d0 * tfac2, &
                        tfac6  = cv*cv + 1.5d0*ca*ca + (xnufam - 1.0d0)* &
                                 (cvp*cvp + 1.5d0*cap*cap))



! initialize
      spair   = 0.0d0
      spairdt = 0.0d0
      spairdd = 0.0d0
      spairda = 0.0d0
      spairdz = 0.0d0

      splas   = 0.0d0
      splasdt = 0.0d0
      splasdd = 0.0d0
      splasda = 0.0d0
      splasdz = 0.0d0

      sphot   = 0.0d0
      sphotdt = 0.0d0
      sphotdd = 0.0d0
      sphotda = 0.0d0
      sphotdz = 0.0d0

      sbrem   = 0.0d0
      sbremdt = 0.0d0
      sbremdd = 0.0d0
      sbremda = 0.0d0
      sbremdz = 0.0d0

      sreco   = 0.0d0
      srecodt = 0.0d0
      srecodd = 0.0d0
      srecoda = 0.0d0
      srecodz = 0.0d0

      snu     = 0.0d0
      dsnudt  = 0.0d0
      dsnudd  = 0.0d0
      dsnuda  = 0.0d0
      dsnudz  = 0.0d0

      if (temp .lt. 1.0e7) return


! to avoid lots of divisions
      deni  = 1.0d0/den
      tempi = 1.0d0/temp
      abari = 1.0d0/abar
      zbari = 1.0d0/zbar


! some composition variables
      ye    = zbar*abari
      xmue  = abar*zbari




! some frequent factors
      t9     = temp * 1.0d-9
      xl     = t9 * con1
      xldt   = 1.0d-9 * con1
      xlp5   = sqrt(xl)
      xl2    = xl*xl
      xl3    = xl2*xl
      xl4    = xl3*xl
      xl5    = xl4*xl
      xl6    = xl5*xl
      xl7    = xl6*xl
      xl8    = xl7*xl
      xl9    = xl8*xl
      xlmp5  = 1.0d0/xlp5
      xlm1   = 1.0d0/xl
      xlm2   = xlm1*xlm1
      xlm3   = xlm1*xlm2
      xlm4   = xlm1*xlm3

      rm     = den*ye
      rmdd   = ye
      rmda   = -rm*abari
      rmdz   = den*abari
      rmi    = 1.0d0/rm

      a0     = rm * 1.0d-9
      a1     = a0**oneth
      zeta   = a1 * xlm1
      zetadt = -a1 * xlm2 * xldt
      a2     = oneth * a1*rmi * xlm1
      zetadd = a2 * rmdd
      zetada = a2 * rmda
      zetadz = a2 * rmdz

      zeta2 = zeta * zeta
      zeta3 = zeta2 * zeta




! pair neutrino section
! for reactions like e+ + e- => nu_e + nubar_e

! equation 2.8
      gl   = 1.0d0 - 13.04d0*xl2 +133.5d0*xl4 +1534.0d0*xl6 +918.6d0*xl8
      gldt = xldt*(-26.08d0*xl +534.0d0*xl3 +9204.0d0*xl5 +7348.8d0*xl7)

! equation 2.7

      a1     = 6.002d19 + 2.084d20*zeta + 1.872d21*zeta2
      a2     = 2.084d20 + 2.0d0*1.872d21*zeta

      if (t9 .lt. 10.0) then
       b1     = exp(-5.5924d0*zeta)
       b2     = -b1*5.5924d0
      else
       b1     = exp(-4.9924d0*zeta)
       b2     = -b1*4.9924d0
      end if

      xnum   = a1 * b1
      c      = a2*b1 + a1*b2
      xnumdt = c*zetadt
      xnumdd = c*zetadd
      xnumda = c*zetada
      xnumdz = c*zetadz

      if (t9 .lt. 10.0) then
       a1   = 9.383d-1*xlm1 - 4.141d-1*xlm2 + 5.829d-2*xlm3
       a2   = -9.383d-1*xlm2 + 2.0d0*4.141d-1*xlm3 - 3.0d0*5.829d-2*xlm4
      else
       a1   = 1.2383d0*xlm1 - 8.141d-1*xlm2
       a2   = -1.2383d0*xlm2 + 2.0d0*8.141d-1*xlm3
      end if

      b1   = 3.0d0*zeta2

      xden   = zeta3 + a1
      xdendt = b1*zetadt + a2*xldt
      xdendd = b1*zetadd
      xdenda = b1*zetada
      xdendz = b1*zetadz

      a1      = 1.0d0/xden
      fpair   = xnum*a1
      fpairdt = (xnumdt - fpair*xdendt)*a1
      fpairdd = (xnumdd - fpair*xdendd)*a1
      fpairda = (xnumda - fpair*xdenda)*a1
      fpairdz = (xnumdz - fpair*xdendz)*a1


! equation 2.6
      a1     = 10.7480d0*xl2 + 0.3967d0*xlp5 + 1.005d0
      a2     = xldt*(2.0d0*10.7480d0*xl + 0.5d0*0.3967d0*xlmp5)
      xnum   = 1.0d0/a1
      xnumdt = -xnum*xnum*a2

      a1     = 7.692d7*xl3 + 9.715d6*xlp5
      a2     = xldt*(3.0d0*7.692d7*xl2 + 0.5d0*9.715d6*xlmp5)

      c      = 1.0d0/a1
      b1     = 1.0d0 + rm*c

      xden   = b1**(-0.3d0)

      d      = -0.3d0*xden/b1
      xdendt = -d*rm*c*c*a2
      xdendd = d*rmdd*c
      xdenda = d*rmda*c
      xdendz = d*rmdz*c

      qpair   = xnum*xden
      qpairdt = xnumdt*xden + xnum*xdendt
      qpairdd = xnum*xdendd
      qpairda = xnum*xdenda
      qpairdz = xnum*xdendz



! equation 2.5
      a1    = exp(-2.0d0*xlm1)
      a2    = a1*2.0d0*xlm2*xldt

      spair   = a1*fpair
      spairdt = a2*fpair + a1*fpairdt
      spairdd = a1*fpairdd
      spairda = a1*fpairda
      spairdz = a1*fpairdz

      a1      = spair
      spair   = gl*a1
      spairdt = gl*spairdt + gldt*a1
      spairdd = gl*spairdd
      spairda = gl*spairda
      spairdz = gl*spairdz

      a1      = tfac4*(1.0d0 + tfac3 * qpair)
      a2      = tfac4*tfac3

      a3      = spair
      spair   = a1*a3
      spairdt = a1*spairdt + a2*qpairdt*a3
      spairdd = a1*spairdd + a2*qpairdd*a3
      spairda = a1*spairda + a2*qpairda*a3
      spairdz = a1*spairdz + a2*qpairdz*a3




! plasma neutrino section
! for collective reactions like gamma_plasmon => nu_e + nubar_e
! equation 4.6

      a1   = 1.019d-6*rm
      a2   = a1**twoth
      a3   = twoth*a2/a1

      b1   =  sqrt(1.0d0 + a2)
      b2   = 1.0d0/b1

      c00  = 1.0d0/(temp*temp*b1)

      gl2   = 1.1095d11 * rm * c00

      gl2dt = -2.0d0*gl2*tempi
      d     = rm*c00*b2*0.5d0*b2*a3*1.019d-6
      gl2dd = 1.1095d11 * (rmdd*c00  - d*rmdd)
      gl2da = 1.1095d11 * (rmda*c00  - d*rmda)
      gl2dz = 1.1095d11 * (rmdz*c00  - d*rmdz)


      gl    = sqrt(gl2)
      gl12  = sqrt(gl)
      gl32  = gl * gl12
      gl72  = gl2 * gl32
      gl6   = gl2 * gl2 * gl2


! equation 4.7
      ft   = 2.4d0 + 0.6d0*gl12 + 0.51d0*gl + 1.25d0*gl32
      gum  = 1.0d0/gl2
      a1   =(0.25d0*0.6d0*gl12 +0.5d0*0.51d0*gl +0.75d0*1.25d0*gl32)*gum
      ftdt = a1*gl2dt
      ftdd = a1*gl2dd
      ftda = a1*gl2da
      ftdz = a1*gl2dz


! equation 4.8
      a1   = 8.6d0*gl2 + 1.35d0*gl72
      a2   = 8.6d0 + 1.75d0*1.35d0*gl72*gum

      b1   = 225.0d0 - 17.0d0*gl + gl2
      b2   = -0.5d0*17.0d0*gl*gum + 1.0d0

      c    = 1.0d0/b1
      fl   = a1*c

      d    = (a2 - fl*b2)*c
      fldt = d*gl2dt
      fldd = d*gl2dd
      flda = d*gl2da
      fldz = d*gl2dz


! equation 4.9 and 4.10
      cc   = log10(2.0d0*rm)
      xlnt = log10(temp)

      xnum   = sixth * (17.5d0 + cc - 3.0d0*xlnt)
      xnumdt = -iln10*0.5d0*tempi
      a2     = iln10*sixth*rmi
      xnumdd = a2*rmdd
      xnumda = a2*rmda
      xnumdz = a2*rmdz

      xden   = sixth * (-24.5d0 + cc + 3.0d0*xlnt)
      xdendt = iln10*0.5d0*tempi
      xdendd = a2*rmdd
      xdenda = a2*rmda
      xdendz = a2*rmdz


! equation 4.11
      if (abs(xnum) .gt. 0.7d0  .or.  xden .lt. 0.0d0) then
       fxy   = 1.0d0
       fxydt = 0.0d0
       fxydd = 0.0d0
       fxydz = 0.0d0
       fxyda = 0.0d0

      else

       a1  = 0.39d0 - 1.25d0*xnum - 0.35d0*sin(4.5d0*xnum)
       a2  = -1.25d0 - 4.5d0*0.35d0*cos(4.5d0*xnum)

       b1  = 0.3d0 * exp(-1.0d0*(4.5d0*xnum + 0.9d0)**2)
       b2  = -b1*2.0d0*(4.5d0*xnum + 0.9d0)*4.5d0

       c   = min(0.0d0, xden - 1.6d0 + 1.25d0*xnum)
       if (c .eq. 0.0) then
        dumdt = 0.0d0
        dumdd = 0.0d0
        dumda = 0.0d0
        dumdz = 0.0d0
       else
        dumdt = xdendt + 1.25d0*xnumdt
        dumdd = xdendd + 1.25d0*xnumdd
        dumda = xdenda + 1.25d0*xnumda
        dumdz = xdendz + 1.25d0*xnumdz
       end if

       d   = 0.57d0 - 0.25d0*xnum
       a3  = c/d
       c00 = exp(-1.0d0*a3**2)

       f1  = -c00*2.0d0*a3/d
       c01 = f1*(dumdt + a3*0.25d0*xnumdt)
       c02 = f1*(dumdd + a3*0.25d0*xnumdd)
       c03 = f1*(dumda + a3*0.25d0*xnumda)
       c04 = f1*(dumdz + a3*0.25d0*xnumdz)

       fxy   = 1.05d0 + (a1 - b1)*c00
       fxydt = (a2*xnumdt -  b2*xnumdt)*c00 + (a1-b1)*c01
       fxydd = (a2*xnumdd -  b2*xnumdd)*c00 + (a1-b1)*c02
       fxyda = (a2*xnumda -  b2*xnumda)*c00 + (a1-b1)*c03
       fxydz = (a2*xnumdz -  b2*xnumdz)*c00 + (a1-b1)*c04

      end if



! equation 4.1 and 4.5
      splas   = (ft + fl) * fxy
      splasdt = (ftdt + fldt)*fxy + (ft+fl)*fxydt
      splasdd = (ftdd + fldd)*fxy + (ft+fl)*fxydd
      splasda = (ftda + flda)*fxy + (ft+fl)*fxyda
      splasdz = (ftdz + fldz)*fxy + (ft+fl)*fxydz

      a2      = exp(-gl)
      a3      = -0.5d0*a2*gl*gum

      a1      = splas
      splas   = a2*a1
      splasdt = a2*splasdt + a3*gl2dt*a1
      splasdd = a2*splasdd + a3*gl2dd*a1
      splasda = a2*splasda + a3*gl2da*a1
      splasdz = a2*splasdz + a3*gl2dz*a1

      a2      = gl6
      a3      = 3.0d0*gl6*gum

      a1      = splas
      splas   = a2*a1
      splasdt = a2*splasdt + a3*gl2dt*a1
      splasdd = a2*splasdd + a3*gl2dd*a1
      splasda = a2*splasda + a3*gl2da*a1
      splasdz = a2*splasdz + a3*gl2dz*a1


      a2      = 0.93153d0 * 3.0d21 * xl9
      a3      = 0.93153d0 * 3.0d21 * 9.0d0*xl8*xldt

      a1      = splas
      splas   = a2*a1
      splasdt = a2*splasdt + a3*a1
      splasdd = a2*splasdd
      splasda = a2*splasda
      splasdz = a2*splasdz




! photoneutrino process section
! for reactions like e- + gamma => e- + nu_e + nubar_e
!                    e+ + gamma => e+ + nu_e + nubar_e
! equation 3.8 for tau, equation 3.6 for cc,
! and table 2 written out for speed
      if (temp .ge. 1.0d7  .and. temp .lt. 1.0d8) then
       tau  =  log10(temp * 1.0d-7)
       cc   =  0.5654d0 + tau
       c00  =  1.008d11
       c01  =  0.0d0
       c02  =  0.0d0
       c03  =  0.0d0
       c04  =  0.0d0
       c05  =  0.0d0
       c06  =  0.0d0
       c10  =  8.156d10
       c11  =  9.728d8
       c12  = -3.806d9
       c13  = -4.384d9
       c14  = -5.774d9
       c15  = -5.249d9
       c16  = -5.153d9
       c20  =  1.067d11
       c21  = -9.782d9
       c22  = -7.193d9
       c23  = -6.936d9
       c24  = -6.893d9
       c25  = -7.041d9
       c26  = -7.193d9
       dd01 =  0.0d0
       dd02 =  0.0d0
       dd03 =  0.0d0
       dd04 =  0.0d0
       dd05 =  0.0d0
       dd11 = -1.879d10
       dd12 = -9.667d9
       dd13 = -5.602d9
       dd14 = -3.370d9
       dd15 = -1.825d9
       dd21 = -2.919d10
       dd22 = -1.185d10
       dd23 = -7.270d9
       dd24 = -4.222d9
       dd25 = -1.560d9

      else if (temp .ge. 1.0d8  .and. temp .lt. 1.0d9) then
       tau   =  log10(temp * 1.0d-8)
       cc   =  1.5654d0
       c00  =  9.889d10
       c01  = -4.524d8
       c02  = -6.088d6
       c03  =  4.269d7
       c04  =  5.172d7
       c05  =  4.910d7
       c06  =  4.388d7
       c10  =  1.813d11
       c11  = -7.556d9
       c12  = -3.304d9
       c13  = -1.031d9
       c14  = -1.764d9
       c15  = -1.851d9
       c16  = -1.928d9
       c20  =  9.750d10
       c21  =  3.484d10
       c22  =  5.199d9
       c23  = -1.695d9
       c24  = -2.865d9
       c25  = -3.395d9
       c26  = -3.418d9
       dd01 = -1.135d8
       dd02 =  1.256d8
       dd03 =  5.149d7
       dd04 =  3.436d7
       dd05 =  1.005d7
       dd11 =  1.652d9
       dd12 = -3.119d9
       dd13 = -1.839d9
       dd14 = -1.458d9
       dd15 = -8.956d8
       dd21 = -1.549d10
       dd22 = -9.338d9
       dd23 = -5.899d9
       dd24 = -3.035d9
       dd25 = -1.598d9

      else if (temp .ge. 1.0d9) then
       tau  =  log10(t9)
       cc   =  1.5654d0
       c00  =  9.581d10
       c01  =  4.107d8
       c02  =  2.305d8
       c03  =  2.236d8
       c04  =  1.580d8
       c05  =  2.165d8
       c06  =  1.721d8
       c10  =  1.459d12
       c11  =  1.314d11
       c12  = -1.169d11
       c13  = -1.765d11
       c14  = -1.867d11
       c15  = -1.983d11
       c16  = -1.896d11
       c20  =  2.424d11
       c21  = -3.669d9
       c22  = -8.691d9
       c23  = -7.967d9
       c24  = -7.932d9
       c25  = -7.987d9
       c26  = -8.333d9
       dd01 =  4.724d8
       dd02 =  2.976d8
       dd03 =  2.242d8
       dd04 =  7.937d7
       dd05 =  4.859d7
       dd11 = -7.094d11
       dd12 = -3.697d11
       dd13 = -2.189d11
       dd14 = -1.273d11
       dd15 = -5.705d10
       dd21 = -2.254d10
       dd22 = -1.551d10
       dd23 = -7.793d9
       dd24 = -4.489d9
       dd25 = -2.185d9
      end if

      taudt = iln10*tempi


! equation 3.7, compute the expensive trig functions only one time
      cos1 = cos(fac1*tau)
      cos2 = cos(fac1*2.0d0*tau)
      cos3 = cos(fac1*3.0d0*tau)
      cos4 = cos(fac1*4.0d0*tau)
      cos5 = cos(fac1*5.0d0*tau)
      last = cos(fac2*tau)

      sin1 = sin(fac1*tau)
      sin2 = sin(fac1*2.0d0*tau)
      sin3 = sin(fac1*3.0d0*tau)
      sin4 = sin(fac1*4.0d0*tau)
      sin5 = sin(fac1*5.0d0*tau)
      xast = sin(fac2*tau)

      a0 = 0.5d0*c00 &
           + c01*cos1 + dd01*sin1 + c02*cos2 + dd02*sin2 &
           + c03*cos3 + dd03*sin3 + c04*cos4 + dd04*sin4 &
           + c05*cos5 + dd05*sin5 + 0.5d0*c06*last

      f0 =  taudt*fac1*(-c01*sin1 + dd01*cos1 - c02*sin2*2.0d0 &
           + dd02*cos2*2.0d0 - c03*sin3*3.0d0 + dd03*cos3*3.0d0 &
           - c04*sin4*4.0d0 + dd04*cos4*4.0d0 &
           - c05*sin5*5.0d0 + dd05*cos5*5.0d0) &
           - 0.5d0*c06*xast*fac2*taudt

      a1 = 0.5d0*c10 &
           + c11*cos1 + dd11*sin1 + c12*cos2 + dd12*sin2 &
           + c13*cos3 + dd13*sin3 + c14*cos4 + dd14*sin4 &
           + c15*cos5 + dd15*sin5 + 0.5d0*c16*last

      f1 = taudt*fac1*(-c11*sin1 + dd11*cos1 - c12*sin2*2.0d0 &
           + dd12*cos2*2.0d0 - c13*sin3*3.0d0 + dd13*cos3*3.0d0 &
           - c14*sin4*4.0d0 + dd14*cos4*4.0d0 - c15*sin5*5.0d0 &
           + dd15*cos5*5.0d0) - 0.5d0*c16*xast*fac2*taudt

      a2 = 0.5d0*c20 &
           + c21*cos1 + dd21*sin1 + c22*cos2 + dd22*sin2 &
           + c23*cos3 + dd23*sin3 + c24*cos4 + dd24*sin4 &
           + c25*cos5 + dd25*sin5 + 0.5d0*c26*last

      f2 = taudt*fac1*(-c21*sin1 + dd21*cos1 - c22*sin2*2.0d0 &
           + dd22*cos2*2.0d0 - c23*sin3*3.0d0 + dd23*cos3*3.0d0 &
           - c24*sin4*4.0d0 + dd24*cos4*4.0d0 - c25*sin5*5.0d0 &
           + dd25*cos5*5.0d0) - 0.5d0*c26*xast*fac2*taudt

! equation 3.4
      dum   = a0 + a1*zeta + a2*zeta2
      dumdt = f0 + f1*zeta + a1*zetadt + f2*zeta2 + 2.0d0*a2*zeta*zetadt
      dumdd = a1*zetadd + 2.0d0*a2*zeta*zetadd
      dumda = a1*zetada + 2.0d0*a2*zeta*zetada
      dumdz = a1*zetadz + 2.0d0*a2*zeta*zetadz

      z      = exp(-cc*zeta)

      xnum   = dum*z
      xnumdt = dumdt*z - dum*z*cc*zetadt
      xnumdd = dumdd*z - dum*z*cc*zetadd
      xnumda = dumda*z - dum*z*cc*zetada
      xnumdz = dumdz*z - dum*z*cc*zetadz

      xden   = zeta3 + 6.290d-3*xlm1 + 7.483d-3*xlm2 + 3.061d-4*xlm3

      dum    = 3.0d0*zeta2
      xdendt = dum*zetadt - xldt*(6.290d-3*xlm2 &
               + 2.0d0*7.483d-3*xlm3 + 3.0d0*3.061d-4*xlm4)
      xdendd = dum*zetadd
      xdenda = dum*zetada
      xdendz = dum*zetadz

      dum      = 1.0d0/xden
      fphot   = xnum*dum
      fphotdt = (xnumdt - fphot*xdendt)*dum
      fphotdd = (xnumdd - fphot*xdendd)*dum
      fphotda = (xnumda - fphot*xdenda)*dum
      fphotdz = (xnumdz - fphot*xdendz)*dum


! equation 3.3
      a0     = 1.0d0 + 2.045d0 * xl
      xnum   = 0.666d0*a0**(-2.066d0)
      xnumdt = -2.066d0*xnum/a0 * 2.045d0*xldt

      dum    = 1.875d8*xl + 1.653d8*xl2 + 8.449d8*xl3 - 1.604d8*xl4
      dumdt  = xldt*(1.875d8 + 2.0d0*1.653d8*xl + 3.0d0*8.449d8*xl2 &
               - 4.0d0*1.604d8*xl3)

      z      = 1.0d0/dum
      xden   = 1.0d0 + rm*z
      xdendt =  -rm*z*z*dumdt
      xdendd =  rmdd*z
      xdenda =  rmda*z
      xdendz =  rmdz*z

      z      = 1.0d0/xden
      qphot = xnum*z
      qphotdt = (xnumdt - qphot*xdendt)*z
      dum      = -qphot*z
      qphotdd = dum*xdendd
      qphotda = dum*xdenda
      qphotdz = dum*xdendz

! equation 3.2
      sphot   = xl5 * fphot
      sphotdt = 5.0d0*xl4*xldt*fphot + xl5*fphotdt
      sphotdd = xl5*fphotdd
      sphotda = xl5*fphotda
      sphotdz = xl5*fphotdz

      a1      = sphot
      sphot   = rm*a1
      sphotdt = rm*sphotdt
      sphotdd = rm*sphotdd + rmdd*a1
      sphotda = rm*sphotda + rmda*a1
      sphotdz = rm*sphotdz + rmdz*a1

      a1      = tfac4*(1.0d0 - tfac3 * qphot)
      a2      = -tfac4*tfac3

      a3      = sphot
      sphot   = a1*a3
      sphotdt = a1*sphotdt + a2*qphotdt*a3
      sphotdd = a1*sphotdd + a2*qphotdd*a3
      sphotda = a1*sphotda + a2*qphotda*a3
      sphotdz = a1*sphotdz + a2*qphotdz*a3

      if (sphot .le. 0.0) then
       sphot   = 0.0d0
       sphotdt = 0.0d0
       sphotdd = 0.0d0
       sphotda = 0.0d0
       sphotdz = 0.0d0
      end if





! bremsstrahlung neutrino section
! for reactions like e- + (z,a) => e- + (z,a) + nu + nubar
!                    n  + n     => n + n + nu + nubar
!                    n  + p     => n + p + nu + nubar
! equation 4.3

      den6   = den * 1.0d-6
      t8     = temp * 1.0d-8
      t812   = sqrt(t8)
      t832   = t8 * t812
      t82    = t8*t8
      t83    = t82*t8
      t85    = t82*t83
      t86    = t85*t8
      t8m1   = 1.0d0/t8
      t8m2   = t8m1*t8m1
      t8m3   = t8m2*t8m1
      t8m5   = t8m3*t8m2
      t8m6   = t8m5*t8m1


      tfermi = 5.9302d9*(sqrt(1.0d0+1.018d0*(den6*ye)**twoth)-1.0d0)

! "weak" degenerate electrons only
      if (temp .gt. 0.3d0 * tfermi) then

! equation 5.3
       dum   = 7.05d6 * t832 + 5.12d4 * t83
       dumdt = (1.5d0*7.05d6*t812 + 3.0d0*5.12d4*t82)*1.0d-8

       z     = 1.0d0/dum
       eta   = rm*z
       etadt = -rm*z*z*dumdt
       etadd = rmdd*z
       etada = rmda*z
       etadz = rmdz*z

       etam1 = 1.0d0/eta
       etam2 = etam1 * etam1
       etam3 = etam2 * etam1


! equation 5.2
       a0    = 23.5d0 + 6.83d4*t8m2 + 7.81d8*t8m5
       f0    = (-2.0d0*6.83d4*t8m3 - 5.0d0*7.81d8*t8m6)*1.0d-8
       xnum  = 1.0d0/a0

       dum   = 1.0d0 + 1.47d0*etam1 + 3.29d-2*etam2
       z     = -1.47d0*etam2 - 2.0d0*3.29d-2*etam3
       dumdt = z*etadt
       dumdd = z*etadd
       dumda = z*etada
       dumdz = z*etadz

       c00   = 1.26d0 * (1.0d0+etam1)
       z     = -1.26d0*etam2
       c01   = z*etadt
       c02   = z*etadd
       c03   = z*etada
       c04   = z*etadz

       z      = 1.0d0/dum
       xden   = c00*z
       xdendt = (c01 - xden*dumdt)*z
       xdendd = (c02 - xden*dumdd)*z
       xdenda = (c03 - xden*dumda)*z
       xdendz = (c04 - xden*dumdz)*z

       fbrem   = xnum + xden
       fbremdt = -xnum*xnum*f0 + xdendt
       fbremdd = xdendd
       fbremda = xdenda
       fbremdz = xdendz


! equation 5.9
       a0    = 230.0d0 + 6.7d5*t8m2 + 7.66d9*t8m5
       f0    = (-2.0d0*6.7d5*t8m3 - 5.0d0*7.66d9*t8m6)*1.0d-8

       z     = 1.0d0 + rm*1.0d-9
       dum   = a0*z
       dumdt = f0*z
       z     = a0*1.0d-9
       dumdd = z*rmdd
       dumda = z*rmda
       dumdz = z*rmdz

       xnum   = 1.0d0/dum
       z      = -xnum*xnum
       xnumdt = z*dumdt
       xnumdd = z*dumdd
       xnumda = z*dumda
       xnumdz = z*dumdz

       c00   = 7.75d5*t832 + 247.0d0*t8**(3.85d0)
       dd00  = (1.5d0*7.75d5*t812 + 3.85d0*247.0d0*t8**(2.85d0))*1.0d-8

       c01   = 4.07d0 + 0.0240d0 * t8**(1.4d0)
       dd01  = 1.4d0*0.0240d0*t8**(0.4d0)*1.0d-8

       c02   = 4.59d-5 * t8**(-0.110d0)
       dd02  = -0.11d0*4.59d-5 * t8**(-1.11d0)*1.0d-8

       z     = den**(0.656d0)
       dum   = c00*rmi  + c01  + c02*z
       dumdt = dd00*rmi + dd01 + dd02*z
       z     = -c00*rmi*rmi
       dumdd = z*rmdd + 0.656d0*c02*den**(-0.454d0)
       dumda = z*rmda
       dumdz = z*rmdz

       xden  = 1.0d0/dum
       z      = -xden*xden
       xdendt = z*dumdt
       xdendd = z*dumdd
       xdenda = z*dumda
       xdendz = z*dumdz

       gbrem   = xnum + xden
       gbremdt = xnumdt + xdendt
       gbremdd = xnumdd + xdendd
       gbremda = xnumda + xdenda
       gbremdz = xnumdz + xdendz


! equation 5.1
       dum    = 0.5738d0*zbar*ye*t86*den
       dumdt  = 0.5738d0*zbar*ye*6.0d0*t85*den*1.0d-8
       dumdd  = 0.5738d0*zbar*ye*t86
       dumda  = -dum*abari
       dumdz  = 0.5738d0*2.0d0*ye*t86*den

       z       = tfac4*fbrem - tfac5*gbrem
       sbrem   = dum * z
       sbremdt = dumdt*z + dum*(tfac4*fbremdt - tfac5*gbremdt)
       sbremdd = dumdd*z + dum*(tfac4*fbremdd - tfac5*gbremdd)
       sbremda = dumda*z + dum*(tfac4*fbremda - tfac5*gbremda)
       sbremdz = dumdz*z + dum*(tfac4*fbremdz - tfac5*gbremdz)




! liquid metal with c12 parameters (not too different for other elements)
! equation 5.18 and 5.16

      else
       u     = fac3 * (log10(den) - 3.0d0)
       a0    = iln10*fac3*deni

! compute the expensive trig functions of equation 5.21 only once
       cos1 = cos(u)
       cos2 = cos(2.0d0*u)
       cos3 = cos(3.0d0*u)
       cos4 = cos(4.0d0*u)
       cos5 = cos(5.0d0*u)

       sin1 = sin(u)
       sin2 = sin(2.0d0*u)
       sin3 = sin(3.0d0*u)
       sin4 = sin(4.0d0*u)
       sin5 = sin(5.0d0*u)

! equation 5.21
       fb =  0.5d0 * 0.17946d0  + 0.00945d0*u + 0.34529d0 &
             - 0.05821d0*cos1 - 0.04969d0*sin1 &
             - 0.01089d0*cos2 - 0.01584d0*sin2 &
             - 0.01147d0*cos3 - 0.00504d0*sin3 &
             - 0.00656d0*cos4 - 0.00281d0*sin4 &
             - 0.00519d0*cos5

       c00 =  a0*(0.00945d0 &
             + 0.05821d0*sin1       - 0.04969d0*cos1 &
             + 0.01089d0*sin2*2.0d0 - 0.01584d0*cos2*2.0d0 &
             + 0.01147d0*sin3*3.0d0 - 0.00504d0*cos3*3.0d0 &
             + 0.00656d0*sin4*4.0d0 - 0.00281d0*cos4*4.0d0 &
             + 0.00519d0*sin5*5.0d0)


! equation 5.22
       ft =  0.5d0 * 0.06781d0 - 0.02342d0*u + 0.24819d0 &
             - 0.00944d0*cos1 - 0.02213d0*sin1 &
             - 0.01289d0*cos2 - 0.01136d0*sin2 &
             - 0.00589d0*cos3 - 0.00467d0*sin3 &
             - 0.00404d0*cos4 - 0.00131d0*sin4 &
             - 0.00330d0*cos5

       c01 = a0*(-0.02342d0 &
             + 0.00944d0*sin1       - 0.02213d0*cos1 &
             + 0.01289d0*sin2*2.0d0 - 0.01136d0*cos2*2.0d0 &
             + 0.00589d0*sin3*3.0d0 - 0.00467d0*cos3*3.0d0 &
             + 0.00404d0*sin4*4.0d0 - 0.00131d0*cos4*4.0d0 &
             + 0.00330d0*sin5*5.0d0)


! equation 5.23
       gb =  0.5d0 * 0.00766d0 - 0.01259d0*u + 0.07917d0 &
             - 0.00710d0*cos1 + 0.02300d0*sin1 &
             - 0.00028d0*cos2 - 0.01078d0*sin2 &
             + 0.00232d0*cos3 + 0.00118d0*sin3 &
             + 0.00044d0*cos4 - 0.00089d0*sin4 &
             + 0.00158d0*cos5

       c02 = a0*(-0.01259d0 &
             + 0.00710d0*sin1       + 0.02300d0*cos1 &
             + 0.00028d0*sin2*2.0d0 - 0.01078d0*cos2*2.0d0 &
             - 0.00232d0*sin3*3.0d0 + 0.00118d0*cos3*3.0d0 &
             - 0.00044d0*sin4*4.0d0 - 0.00089d0*cos4*4.0d0 &
             - 0.00158d0*sin5*5.0d0)


! equation 5.24
       gt =  -0.5d0 * 0.00769d0  - 0.00829d0*u + 0.05211d0 &
             + 0.00356d0*cos1 + 0.01052d0*sin1 &
             - 0.00184d0*cos2 - 0.00354d0*sin2 &
             + 0.00146d0*cos3 - 0.00014d0*sin3 &
             + 0.00031d0*cos4 - 0.00018d0*sin4 &
             + 0.00069d0*cos5

       c03 = a0*(-0.00829d0 &
             - 0.00356d0*sin1       + 0.01052d0*cos1 &
             + 0.00184d0*sin2*2.0d0 - 0.00354d0*cos2*2.0d0 &
             - 0.00146d0*sin3*3.0d0 - 0.00014d0*cos3*3.0d0 &
             - 0.00031d0*sin4*4.0d0 - 0.00018d0*cos4*4.0d0 &
             - 0.00069d0*sin5*5.0d0)


       dum   = 2.275d-1 * zbar * zbar*t8m1 * (den6*abari)**oneth
       dumdt = -dum*tempi
       dumdd = oneth*dum*deni
       dumda = -oneth*dum*abari
       dumdz = 2.0d0*dum*zbari

       gm1   = 1.0d0/dum
       gm2   = gm1*gm1
       gm13  = gm1**oneth
       gm23  = gm13 * gm13
       gm43  = gm13*gm1
       gm53  = gm23*gm1


! equation 5.25 and 5.26
       v  = -0.05483d0 - 0.01946d0*gm13 + 1.86310d0*gm23 - 0.78873d0*gm1
       a0 = oneth*0.01946d0*gm43 - twoth*1.86310d0*gm53 + 0.78873d0*gm2

       w  = -0.06711d0 + 0.06859d0*gm13 + 1.74360d0*gm23 - 0.74498d0*gm1
       a1 = -oneth*0.06859d0*gm43 - twoth*1.74360d0*gm53 + 0.74498d0*gm2


! equation 5.19 and 5.20
       fliq   = v*fb + (1.0d0 - v)*ft
       fliqdt = a0*dumdt*(fb - ft)
       fliqdd = a0*dumdd*(fb - ft) + v*c00 + (1.0d0 - v)*c01
       fliqda = a0*dumda*(fb - ft)
       fliqdz = a0*dumdz*(fb - ft)

       gliq   = w*gb + (1.0d0 - w)*gt
       gliqdt = a1*dumdt*(gb - gt)
       gliqdd = a1*dumdd*(gb - gt) + w*c02 + (1.0d0 - w)*c03
       gliqda = a1*dumda*(gb - gt)
       gliqdz = a1*dumdz*(gb - gt)


! equation 5.17
       dum    = 0.5738d0*zbar*ye*t86*den
       dumdt  = 0.5738d0*zbar*ye*6.0d0*t85*den*1.0d-8
       dumdd  = 0.5738d0*zbar*ye*t86
       dumda  = -dum*abari
       dumdz  = 0.5738d0*2.0d0*ye*t86*den

       z       = tfac4*fliq - tfac5*gliq
       sbrem   = dum * z
       sbremdt = dumdt*z + dum*(tfac4*fliqdt - tfac5*gliqdt)
       sbremdd = dumdd*z + dum*(tfac4*fliqdd - tfac5*gliqdd)
       sbremda = dumda*z + dum*(tfac4*fliqda - tfac5*gliqda)
       sbremdz = dumdz*z + dum*(tfac4*fliqdz - tfac5*gliqdz)

      end if




! recombination neutrino section
! for reactions like e- (continuum) => e- (bound) + nu_e + nubar_e
! equation 6.11 solved for nu
      xnum   = 1.10520d8 * den * ye /(temp*sqrt(temp))
      xnumdt = -1.50d0*xnum*tempi
      xnumdd = xnum*deni
      xnumda = -xnum*abari
      xnumdz = xnum*zbari

! the chemical potential
      nu   = ifermi12(xnum)

! a0 is d(nu)/d(xnum)
      a0 = 1.0d0/(0.5d0*zfermim12(nu))
      nudt = a0*xnumdt
      nudd = a0*xnumdd
      nuda = a0*xnumda
      nudz = a0*xnumdz

      nu2  = nu * nu
      nu3  = nu2 * nu

! table 12
      if (nu .ge. -20.0  .and. nu .lt. 0.0) then
       a1 = 1.51d-2
       a2 = 2.42d-1
       a3 = 1.21d0
       b  = 3.71d-2
       c  = 9.06e-1
       d  = 9.28d-1
       f1 = 0.0d0
       f2 = 0.0d0
       f3 = 0.0d0
      else if (nu .ge. 0.0  .and. nu .le. 10.0) then
       a1 = 1.23d-2
       a2 = 2.66d-1
       a3 = 1.30d0
       b  = 1.17d-1
       c  = 8.97e-1
       d  = 1.77d-1
       f1 = -1.20d-2
       f2 = 2.29d-2
       f3 = -1.04d-3
      end if


! equation 6.7, 6.13 and 6.14
      if (nu .ge. -20.0  .and.  nu .le. 10.0) then

       zeta   = 1.579d5*zbar*zbar*tempi
       zetadt = -zeta*tempi
       zetadd = 0.0d0
       zetada = 0.0d0
       zetadz = 2.0d0*zeta*zbari

       c00    = 1.0d0/(1.0d0 + f1*nu + f2*nu2 + f3*nu3)
       c01    = f1 + f2*2.0d0*nu + f3*3.0d0*nu2
       dum    = zeta*c00
       dumdt  = zetadt*c00 + zeta*c01*nudt
       dumdd  = zeta*c01*nudd
       dumda  = zeta*c01*nuda
       dumdz  = zetadz*c00 + zeta*c01*nudz


       z      = 1.0d0/dum
       dd00   = dum**(-2.25)
       dd01   = dum**(-4.55)
       c00    = a1*z + a2*dd00 + a3*dd01
       c01    = -(a1*z + 2.25*a2*dd00 + 4.55*a3*dd01)*z


       z      = exp(c*nu)
       dd00   = b*z*(1.0d0 + d*dum)
       gum    = 1.0d0 + dd00
       gumdt  = dd00*c*nudt + b*z*d*dumdt
       gumdd  = dd00*c*nudd + b*z*d*dumdd
       gumda  = dd00*c*nuda + b*z*d*dumda
       gumdz  = dd00*c*nudz + b*z*d*dumdz


       z   = exp(nu)
       a1  = 1.0d0/gum

       bigj   = c00 * z * a1
       bigjdt = c01*dumdt*z*a1 + c00*z*nudt*a1 - c00*z*a1*a1 * gumdt
       bigjdd = c01*dumdd*z*a1 + c00*z*nudd*a1 - c00*z*a1*a1 * gumdd
       bigjda = c01*dumda*z*a1 + c00*z*nuda*a1 - c00*z*a1*a1 * gumda
       bigjdz = c01*dumdz*z*a1 + c00*z*nudz*a1 - c00*z*a1*a1 * gumdz


! equation 6.5
       z     = exp(zeta + nu)
       dum   = 1.0d0 + z
       a1    = 1.0d0/dum
       a2    = 1.0d0/bigj

       sreco   = tfac6 * 2.649d-18 * ye * zbar**13 * den * bigj*a1
       srecodt = sreco*(bigjdt*a2 - z*(zetadt + nudt)*a1)
       srecodd = sreco*(1.0d0*deni + bigjdd*a2 - z*(zetadd + nudd)*a1)
       srecoda = sreco*(-1.0d0*abari + bigjda*a2 - z*(zetada+nuda)*a1)
       srecodz = sreco*(14.0d0*zbari + bigjdz*a2 - z*(zetadz+nudz)*a1)

      end if


! convert from erg/cm^3/s to erg/g/s
! comment these out to duplicate the itoh et al plots

      spair   = spair*deni
      spairdt = spairdt*deni
      spairdd = spairdd*deni - spair*deni
      spairda = spairda*deni
      spairdz = spairdz*deni

      splas   = splas*deni
      splasdt = splasdt*deni
      splasdd = splasdd*deni - splas*deni
      splasda = splasda*deni
      splasdz = splasdz*deni

      sphot   = sphot*deni
      sphotdt = sphotdt*deni
      sphotdd = sphotdd*deni - sphot*deni
      sphotda = sphotda*deni
      sphotdz = sphotdz*deni

      sbrem   = sbrem*deni
      sbremdt = sbremdt*deni
      sbremdd = sbremdd*deni - sbrem*deni
      sbremda = sbremda*deni
      sbremdz = sbremdz*deni

      sreco   = sreco*deni
      srecodt = srecodt*deni
      srecodd = srecodd*deni - sreco*deni
      srecoda = srecoda*deni
      srecodz = srecodz*deni


! the total neutrino loss rate
      snu    =  splas + spair + sphot + sbrem + sreco
      dsnudt =  splasdt + spairdt + sphotdt + sbremdt + srecodt
      dsnudd =  splasdd + spairdd + sphotdd + sbremdd + srecodd
      dsnuda =  splasda + spairda + sphotda + sbremda + srecoda
      dsnudz =  splasdz + spairdz + sphotdz + sbremdz + srecodz

      return
      end






      double precision function ifermi12(f)
      include 'implno.dek'

! this routine applies a rational function expansion to get the inverse
! fermi-dirac integral of order 1/2 when it is equal to f.
! maximum error is 4.19d-9.   reference: antia apjs 84,101 1993

! declare
      integer          i,m1,k1,m2,k2
      double precision f,an,a1(12),b1(12),a2(12),b2(12),rn,den,ff


! load the coefficients of the expansion
      data  an,m1,k1,m2,k2 /0.5d0, 4, 3, 6, 5/
      data  (a1(i),i=1,5)/ 1.999266880833d4,   5.702479099336d3, &
           6.610132843877d2,   3.818838129486d1, &
           1.0d0/
      data  (b1(i),i=1,4)/ 1.771804140488d4,  -2.014785161019d3, &
           9.130355392717d1,  -1.670718177489d0/
      data  (a2(i),i=1,7)/-1.277060388085d-2,  7.187946804945d-2, &
                          -4.262314235106d-1,  4.997559426872d-1, &
                          -1.285579118012d0,  -3.930805454272d-1, &
           1.0d0/
      data  (b2(i),i=1,6)/-9.745794806288d-3,  5.485432756838d-2, &
                          -3.299466243260d-1,  4.077841975923d-1, &
                          -1.145531476975d0,  -6.067091689181d-2/


      if (f .lt. 4.0d0) then
       rn  = f + a1(m1)
       do i=m1-1,1,-1
        rn  = rn*f + a1(i)
       enddo
       den = b1(k1+1)
       do i=k1,1,-1
        den = den*f + b1(i)
       enddo
       ifermi12 = log(f * rn/den)

      else
       ff = 1.0d0/f**(1.0d0/(1.0d0 + an))
       rn = ff + a2(m2)
       do i=m2-1,1,-1
        rn = rn*ff + a2(i)
       enddo
       den = b2(k2+1)
       do i=k2,1,-1
        den = den*ff + b2(i)
       enddo
       ifermi12 = rn/(den*ff)
      end if
      return
      end






      double precision function zfermim12(x)
      include 'implno.dek'

! this routine applies a rational function expansion to get the fermi-dirac
! integral of order -1/2 evaluated at x. maximum error is 1.23d-12.
! reference: antia apjs 84,101 1993

! declare
      integer          i,m1,k1,m2,k2
      double precision x,an,a1(12),b1(12),a2(12),b2(12),rn,den,xx

! load the coefficients of the expansion
      data  an,m1,k1,m2,k2 /-0.5d0, 7, 7, 11, 11/
      data  (a1(i),i=1,8)/ 1.71446374704454d7,    3.88148302324068d7, &
                           3.16743385304962d7,    1.14587609192151d7, &
                           1.83696370756153d6,    1.14980998186874d5, &
                           1.98276889924768d3,    1.0d0/
      data  (b1(i),i=1,8)/ 9.67282587452899d6,    2.87386436731785d7, &
                           3.26070130734158d7,    1.77657027846367d7, &
                           4.81648022267831d6,    6.13709569333207d5, &
                           3.13595854332114d4,    4.35061725080755d2/
      data (a2(i),i=1,12)/-4.46620341924942d-15, -1.58654991146236d-12, &
                          -4.44467627042232d-10, -6.84738791621745d-8, &
                          -6.64932238528105d-6,  -3.69976170193942d-4, &
                          -1.12295393687006d-2,  -1.60926102124442d-1, &
                          -8.52408612877447d-1,  -7.45519953763928d-1, &
                           2.98435207466372d0,    1.0d0/
      data (b2(i),i=1,12)/-2.23310170962369d-15, -7.94193282071464d-13, &
                          -2.22564376956228d-10, -3.43299431079845d-8, &
                          -3.33919612678907d-6,  -1.86432212187088d-4, &
                          -5.69764436880529d-3,  -8.34904593067194d-2, &
                          -4.78770844009440d-1,  -4.99759250374148d-1, &
                           1.86795964993052d0,    4.16485970495288d-1/


      if (x .lt. 2.0d0) then
       xx = exp(x)
       rn = xx + a1(m1)
       do i=m1-1,1,-1
        rn = rn*xx + a1(i)
       enddo
       den = b1(k1+1)
       do i=k1,1,-1
        den = den*xx + b1(i)
       enddo
       zfermim12 = xx * rn/den
!
      else
       xx = 1.0d0/(x*x)
       rn = xx + a2(m2)
       do i=m2-1,1,-1
        rn = rn*xx + a2(i)
       enddo
       den = b2(k2+1)
       do i=k2,1,-1
        den = den*xx + b2(i)
       enddo
       zfermim12 = sqrt(x)*rn/den
      end if
      return
      end








      subroutine mazurek(btemp,bden,y56,ye,rn56ec,sn56ec)
      include 'implno.dek'

! this routine evaluates mazurel's 1973 fits for the ni56 electron
! capture rate rn56ec and neutrino loss rate sn56ec

! input:
! y56 = nickel56 molar abundance
! ye  = electron to baryon number, zbar/abar

! output:
! rn56ec = ni56 electron capture rate
! sn56ec = ni56 neutrino loss rate

! declare
      integer          ifirst,jp,kp,jr,jd,ii,ik,ij,j,k
      double precision btemp,bden,y56,ye,rn56ec,sn56ec, &
                       rnt(2),rne(2,7),datn(2,6,7), &
                       tv(7),rv(6),rfdm(4),rfd0(4),rfd1(4),rfd2(4), &
                       tfdm(5),tfd0(5),tfd1(5),tfd2(5), &
                       t9,r,rfm,rf0,rf1,rf2,dfacm,dfac0,dfac1,dfac2, &
                       tfm,tf0,tf1,tf2,tfacm,tfac0,tfac1,tfac2

! initialize
      data  rv /6.0, 7.0, 8.0, 9.0, 10.0, 11.0/
      data  tv /2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0/
      data (((datn(ii,ik,ij),ik=1,6),ij=1,7),ii=1,1) / &
          -3.98, -2.84, -1.41,  0.20,  1.89,  3.63, &
          -3.45, -2.62, -1.32,  0.22,  1.89,  3.63, &
          -2.68, -2.30, -1.19,  0.27,  1.91,  3.62, &
          -2.04, -1.87, -1.01,  0.34,  1.94,  3.62, &
          -1.50, -1.41, -0.80,  0.45,  1.99,  3.60, &
          -1.00, -0.95, -0.54,  0.60,  2.06,  3.58, &
          -0.52, -0.49, -0.21,  0.79,  2.15,  3.55 /
      data (((datn(ii,ik,ij),ik=1,6),ij=1,7),ii=2,2) / &
          -3.68, -2.45, -0.80,  1.12,  3.13,  5.19, &
          -2.91, -2.05, -0.64,  1.16,  3.14,  5.18, &
          -1.95, -1.57, -0.40,  1.24,  3.16,  5.18, &
          -1.16, -0.99, -0.11,  1.37,  3.20,  5.18, &
          -0.48, -0.40,  0.22,  1.54,  3.28,  5.16, &
           0.14,  0.19,  0.61,  1.78,  3.38,  5.14, &
           0.75,  0.78,  1.06,  2.07,  3.51,  5.11 /
      data  ifirst /0/

! first time; calculate the cubic interp parameters for ni56 electron capture
      if (ifirst .eq. 0) then
       ifirst = 1
       do k=2,4
        rfdm(k)=1./((rv(k-1)-rv(k))*(rv(k-1)-rv(k+1))*(rv(k-1)-rv(k+2)))
        rfd0(k)=1./((rv(k)-rv(k-1))*(rv(k)-rv(k+1))*(rv(k)-rv(k+2)))
        rfd1(k)=1./((rv(k+1)-rv(k-1))*(rv(k+1)-rv(k))*(rv(k+1)-rv(k+2)))
        rfd2(k)=1./((rv(k+2)-rv(k-1))*(rv(k+2)-rv(k))*(rv(k+2)-rv(k+1)))
       enddo
       do j=2,5
        tfdm(j)=1./((tv(j-1)-tv(j))*(tv(j-1)-tv(j+1))*(tv(j-1)-tv(j+2)))
        tfd0(j)=1./((tv(j)-tv(j-1))*(tv(j)-tv(j+1))*(tv(j)-tv(j+2)))
        tfd1(j)=1./((tv(j+1)-tv(j-1))*(tv(j+1)-tv(j))*(tv(j+1)-tv(j+2)))
        tfd2(j)=1./((tv(j+2)-tv(j-1))*(tv(j+2)-tv(j))*(tv(j+2)-tv(j+1)))
       enddo
      end if

! calculate ni56 electron capture and neutrino loss rates
      rn56ec = 0.0
      sn56ec = 0.0
      if ( (btemp .lt. 2.0e9) .or. (bden*ye .lt. 1.0e6)) return
      t9    = max(btemp,1.4d10) * 1.0d-9
      r     = max(6.0d0,min(11.0d0,log10(bden*ye)))
      jp    = min(max(2,int(0.5d0*t9)),5)
      kp    = min(max(2,int(r)-5),4)
      rfm   = r - rv(kp-1)
      rf0   = r - rv(kp)
      rf1   = r - rv(kp+1)
      rf2   = r - rv(kp+2)
      dfacm = rf0*rf1*rf2*rfdm(kp)
      dfac0 = rfm*rf1*rf2*rfd0(kp)
      dfac1 = rfm*rf0*rf2*rfd1(kp)
      dfac2 = rfm*rf0*rf1*rfd2(kp)
      tfm   = t9 - tv(jp-1)
      tf0   = t9 - tv(jp)
      tf1   = t9 - tv(jp+1)
      tf2   = t9 - tv(jp+2)
      tfacm = tf0*tf1*tf2*tfdm(jp)
      tfac0 = tfm*tf1*tf2*tfd0(jp)
      tfac1 = tfm*tf0*tf2*tfd1(jp)
      tfac2 = tfm*tf0*tf1*tfd2(jp)

! evaluate the spline fits
      do jr = 1,2
       do jd = jp-1,jp+2
        rne(jr,jd) =   dfacm*datn(jr,kp-1,jd) + dfac0*datn(jr,kp,jd) &
                     + dfac1*datn(jr,kp+1,jd) + dfac2*datn(jr,kp+2,jd)
       enddo
       rnt(jr) =  tfacm*rne(jr,jp-1) + tfac0*rne(jr,jp) &
                + tfac1*rne(jr,jp+1) + tfac2*rne(jr,jp+2)
      enddo

! set the output
      rn56ec = 10.0d0**rnt(1)
      sn56ec = 6.022548d+23 * 8.18683d-7 * y56 * 10.0d0**rnt(2)
      return
      end






      subroutine ecapnuc02(temp,den,abar,zbar,rpen,rnep,spenc,snepc)
      include 'implno.dek'
      include 'vector_eos.dek'

! given the temperature temp, desnity den, mean atomic number abar and mean
! atomic charge zbar, this routine returns 
! electron capture on protons rpen (captures/sec/proton),
! positron capture on neutrons rnep (captures/sec/neutron),
! and their associated neutrino energy loss rates
! spenc (erg/sec/proton) and snepc (erg/sec/neutron)

! declare the pass
      double precision temp,den,abar,zbar,rpen,rnep,spenc,snepc


! local variables
      integer          iflag
      double precision t9,t5,qn,etaef,etael,zetan,eta,etael2, &
                       etael3,etael4,f1l,f2l,f3l,f4l,f5l,f1g, &
                       f2g,f3g,f4g,f5g,exmeta,eta2,eta3,eta4, &
                       fac0,fac1,fac2,fac3,rie1,rie2,facv0,facv1, &
                       facv2,facv3,facv4,rjv1,rjv2,spen,snep, &
                       pi2,exeta,zetan2,f0,etael5, &
                       qn1,ft,twoln,cmk5,cmk6,bk,pi,qn2,c2me, &
                       xmp,xmn,qndeca,tmean,etakep
      parameter        (qn1    = -2.0716446d-06, &
                        ft     = 1083.9269d0, &
                        twoln  = 0.6931472d0, &
                        cmk5   = 1.3635675d-49, &
                        cmk6   = 2.2993864d-59, &
                        bk     = 1.38062e-16, &
                        pi     = 3.1415927d0, &
                        pi2    = pi * pi, &
                        qn2    = 2.0716446d-06, &
                        c2me   = 8.1872665d-07, &
                        xmp    = 1.6726485d-24, &
                        xmn    = 1.6749543d-24, &
                        qndeca = 1.2533036d-06, &
                        tmean  = 886.7d0)
!                       tmean  = 935.14d0)



! tmean and qndeca are the mean lifetime and decay energy of the neutron
! xmp,xnp are masses of the p and n in grams.
! c2me is the constant used to convert the neutrino energy
! loss rate from mec2/s (as in the paper) to ergs/particle/sec.

! initialize
      rpen  = 0.0d0
      rnep  = 0.0d0
      spen  = 0.0d0
      snep  = 0.0d0
      t9    = temp * 1.0d-9
      iflag = 0
      qn    = qn1


! call an eos to get the chemical potential
       temp_row(1) = temp
       den_row(1)  = den
       abar_row(1) = abar
       zbar_row(1) = zbar
       jlo_eos = 1
       jhi_eos = 1
!       call helmeos
       call helmeosaprox13
       etakep = etaele_row(1)


! chemical potential including the electron rest mass
      etaef = etakep + c2me/bk/temp


! iflag=1 is for electrons,  iflag=2 is for positrons
502   iflag = iflag + 1
      if (iflag.eq.1) etael = qn2/bk/temp
      if (iflag.eq.2) etael = c2me/bk/temp
      if (iflag.eq.2) etaef = -etaef

      t5    = temp*temp*temp*temp*temp
      zetan = qn/bk/temp
      eta   = etaef - etael

! protect from overflowing with large eta values
      if (eta .le. 6.8e+02) then
       exeta = exp(eta)
      else
       exeta = 0.0d0
      end if
      etael2 = etael*etael
      etael3 = etael2*etael
      etael4 = etael3*etael
      etael5 = etael4*etael
      zetan2 = zetan*zetan
      if (eta .le. 6.8e+02) then
       f0 = log(1.0d0 + exeta)
      else
       f0 = eta
      end if

! if eta le. 0., the following fermi integrals apply
      f1l = exeta
      f2l = 2.0d0   * f1l
      f3l = 6.0d0   * f1l
      f4l = 24.0d0  * f1l
      f5l = 120.0d0 * f1l

! if eta gt. 0., the following fermi integrals apply:
      f1g = 0.0d0
      f2g = 0.0d0
      f3g = 0.0d0
      f4g = 0.0d0
      f5g = 0.0d0
      if (eta .gt. 0.0) then
       exmeta = dexp(-eta)
       eta2   = eta*eta
       eta3   = eta2*eta
       eta4   = eta3*eta
       f1g = 0.5d0*eta2 + 2.0d0 - exmeta
       f2g = eta3/3.0d0 + 4.0d0*eta + 2.0d0*exmeta
       f3g = 0.25d0*eta4 + 0.5d0*pi2*eta2 + 12.0d0 - 6.0d0*exmeta
       f4g = 0.2d0*eta4*eta + 2.0d0*pi2/3.0d0*eta3 + 48.0d0*eta &
             + 24.0d0*exmeta
       f5g = eta4*eta2/6.0d0 + 5.0d0/6.0d0*pi2*eta4 &
             + 7.0d0/6.0d0*pi2*eta2  + 240.0d0 -120.d0*exmeta
       end if

! factors which are multiplied by the fermi integrals
      fac3 = 2.0d0*zetan + 4.0d0*etael
      fac2 = 6.0d0*etael2 + 6.0d0*etael*zetan + zetan2
      fac1 = 4.0d0*etael3 + 6.0d0*etael2*zetan + 2.0d0*etael*zetan2
      fac0 = etael4 + 2.0d0*zetan*etael3 + etael2*zetan2

! electron capture rates onto protons with no blocking
      rie1 = f4l + fac3*f3l + fac2*f2l + fac1*f1l + fac0*f0
      rie2 = f4g + fac3*f3g + fac2*f2g + fac1*f1g + fac0*f0

! neutrino emission rate for electron capture:
      facv4 = 5.0d0*etael + 3.0d0*zetan
      facv3 = 10.0d0*etael2 + 12.0d0*etael*zetan + 3.0d0*zetan2
      facv2 = 10.0d0*etael3 + 18.0d0*etael2*zetan &
              + 9.0d0*etael*zetan2 + zetan2*zetan
      facv1 = 5.0d0*etael4 + 12.0d0*etael3*zetan &
              + 9.0d0*etael2*zetan2 + 2.0d0*etael*zetan2*zetan
      facv0 = etael5 + 3.0d0*etael4*zetan &
              + 3.0d0*etael3*zetan2 + etael2*zetan2*zetan
      rjv1  = f5l + facv4*f4l + facv3*f3l &
              + facv2*f2l + facv1*f1l + facv0*f0
      rjv2  = f5g + facv4*f4g + facv3*f3g &
              + facv2*f2g + facv1*f1g + facv0*f0

! for electrons capture onto protons
      if (iflag.eq.2) go to 503
      if (eta.gt.0.) go to 505
      rpen  = twoln*cmk5*t5*rie1/ft
      spen  = twoln*cmk6*t5*temp*rjv1/ft
      spenc = twoln*cmk6*t5*temp*rjv1/ft*c2me
      go to 504
505   rpen = twoln*cmk5*t5*rie2/ft
      spen = twoln*cmk6*t5*temp*rjv2/ft
      spenc = twoln*cmk6*t5*temp*rjv2/ft*c2me
504   continue
      qn = qn2
      go to 502

! for positrons capture onto neutrons
503   if (eta.gt.0.) go to 507
      rnep  = twoln*cmk5*t5*rie1/ft
      snep  = twoln*cmk6*t5*temp*rjv1/ft
      snepc = twoln*cmk6*t5*temp*rjv1/ft*c2me
!      if (rho.lt.1.0e+06) snep=snep+qndeca*xn(9)/xmn/tmean
      go to 506
507   rnep  = twoln*cmk5*t5*rie2/ft
      snep  = twoln*cmk6*t5*temp*rjv2/ft
      snepc = twoln*cmk6*t5*temp*rjv2/ft*c2me
!      if (rho.lt.1.0e+06) snep=snep+qndeca*xn(9)/xmn/tmean
506   continue
      return
      end







      subroutine ecapnuc(etakep,temp,rpen,rnep,spenc,snepc)
      include 'implno.dek'

! given the electron degeneracy parameter etakep (chemical potential
! without the electron's rest mass divided by kt) and the temperature
! temp, this routine calculates rates for
! electron capture on protons rpen (captures/sec/proton),
! positron capture on neutrons rnep (captures/sec/neutron),
! and their associated neutrino energy loss rates
! spenc (erg/sec/proton) and snepc (erg/sec/neutron)

! declare the pass
      double precision etakep,temp,rpen,rnep,spenc,snepc


! local variables
      integer          iflag
      double precision t9,t5,qn,etaef,etael,zetan,eta,etael2, &
                       etael3,etael4,f1l,f2l,f3l,f4l,f5l,f1g, &
                       f2g,f3g,f4g,f5g,exmeta,eta2,eta3,eta4, &
                       fac0,fac1,fac2,fac3,rie1,rie2,facv0,facv1, &
                       facv2,facv3,facv4,rjv1,rjv2,spen,snep, &
                       pi2,exeta,zetan2,f0,etael5,bktinv, &
                       qn1,ftinv,twoln,cmk5,cmk6,bk,pi,qn2,c2me, &
                       xmp,xmn,qndeca,tmean
      parameter        (qn1    = -2.0716446d-06, &
                        ftinv  = 1.0d0/1083.9269d0, &
                        twoln  = 0.6931472d0, &
                        cmk5   = 1.3635675d-49, &
                        cmk6   = 2.2993864d-59, &
                        bk     = 1.38062e-16, &
                        pi     = 3.1415927d0, &
                        pi2    = pi * pi, &
                        qn2    = 2.0716446d-06, &
                        c2me   = 8.1872665d-07, &
                        xmp    = 1.6726485d-24, &
                        xmn    = 1.6749543d-24, &
                        qndeca = 1.2533036d-06, &
                        tmean  = 886.7d0)
!     3                  tmean  = 935.14d0)

      double precision third,sixth
      parameter        (third = 1.0d0/3.0d0, &
                        sixth = 1.0d0/6.0d0)



! tmean and qndeca are the mean lifetime and decay energy of the neutron
! xmp,xnp are masses of the p and n in grams.
! c2me is the constant used to convert the neutrino energy
! loss rate from mec2/s (as in the paper) to ergs/particle/sec.

! initialize
      rpen   = 0.0d0
      rnep   = 0.0d0
      spen   = 0.0d0
      snep   = 0.0d0
      t9     = temp * 1.0d-9
      bktinv = 1.0d0/(bk *temp)
      iflag  = 0
      qn     = qn1


! chemical potential including the electron rest mass
      etaef = etakep + c2me*bktinv


! iflag=1 is for electrons,  iflag=2 is for positrons
502   iflag = iflag + 1
      if (iflag.eq.1) etael = qn2*bktinv
      if (iflag.eq.2) then
       etael = c2me*bktinv
       etaef = -etaef
      endif

      t5    = temp*temp*temp*temp*temp
      zetan = qn*bktinv
      eta   = etaef - etael

! protect from overflowing with large eta values
      if (eta .le. 6.8e+02) then
       exeta = exp(eta)
      else
       exeta = 0.0d0
      end if
      etael2 = etael*etael
      etael3 = etael2*etael
      etael4 = etael3*etael
      etael5 = etael4*etael
      zetan2 = zetan*zetan
      if (eta .le. 6.8e+02) then
       f0 = log(1.0d0 + exeta)
      else
       f0 = eta
      end if

! if eta le. 0., the following fermi integrals apply
      f1l = exeta
      f2l = 2.0d0   * f1l
      f3l = 6.0d0   * f1l
      f4l = 24.0d0  * f1l
      f5l = 120.0d0 * f1l

! if eta gt. 0., the following fermi integrals apply:
      f1g = 0.0d0
      f2g = 0.0d0
      f3g = 0.0d0
      f4g = 0.0d0
      f5g = 0.0d0
      if (eta .gt. 0.0) then
       exmeta = dexp(-eta)
       eta2   = eta*eta
       eta3   = eta2*eta
       eta4   = eta3*eta
       f1g = 0.5d0*eta2 + 2.0d0 - exmeta
       f2g = eta3*third + 4.0d0*eta + 2.0d0*exmeta
       f3g = 0.25d0*eta4 + 0.5d0*pi2*eta2 + 12.0d0 - 6.0d0*exmeta
       f4g = 0.2d0*eta4*eta + 2.0d0*pi2*third*eta3 + 48.0d0*eta &
             + 24.0d0*exmeta
       f5g = eta4*eta2*sixth + 5.0d0*sixth*pi2*eta4 &
             + 7.0d0*sixth*pi2*eta2  + 240.0d0 -120.d0*exmeta
       end if

! factors which are multiplied by the fermi integrals
      fac3 = 2.0d0*zetan + 4.0d0*etael
      fac2 = 6.0d0*etael2 + 6.0d0*etael*zetan + zetan2
      fac1 = 4.0d0*etael3 + 6.0d0*etael2*zetan + 2.0d0*etael*zetan2
      fac0 = etael4 + 2.0d0*zetan*etael3 + etael2*zetan2

! electron capture rates onto protons with no blocking
      rie1 = f4l + fac3*f3l + fac2*f2l + fac1*f1l + fac0*f0
      rie2 = f4g + fac3*f3g + fac2*f2g + fac1*f1g + fac0*f0

! neutrino emission rate for electron capture:
      facv4 = 5.0d0*etael + 3.0d0*zetan
      facv3 = 10.0d0*etael2 + 12.0d0*etael*zetan + 3.0d0*zetan2
      facv2 = 10.0d0*etael3 + 18.0d0*etael2*zetan &
              + 9.0d0*etael*zetan2 + zetan2*zetan
      facv1 = 5.0d0*etael4 + 12.0d0*etael3*zetan &
              + 9.0d0*etael2*zetan2 + 2.0d0*etael*zetan2*zetan
      facv0 = etael5 + 3.0d0*etael4*zetan &
              + 3.0d0*etael3*zetan2 + etael2*zetan2*zetan
      rjv1  = f5l + facv4*f4l + facv3*f3l &
              + facv2*f2l + facv1*f1l + facv0*f0
      rjv2  = f5g + facv4*f4g + facv3*f3g &
              + facv2*f2g + facv1*f1g + facv0*f0

! for electrons capture onto protons
      if (iflag.eq.2) go to 503
      if (eta.gt.0.) go to 505
      rpen  = twoln*cmk5*t5*rie1*ftinv
      spen  = twoln*cmk6*t5*temp*rjv1*ftinv
      spenc = twoln*cmk6*t5*temp*rjv1*ftinv*c2me
      go to 504
505   rpen = twoln*cmk5*t5*rie2*ftinv
      spen = twoln*cmk6*t5*temp*rjv2*ftinv
      spenc = twoln*cmk6*t5*temp*rjv2*ftinv*c2me
504   continue
      qn = qn2
      go to 502

! for positrons capture onto neutrons
503   if (eta.gt.0.) go to 507
      rnep  = twoln*cmk5*t5*rie1*ftinv
      snep  = twoln*cmk6*t5*temp*rjv1*ftinv
      snepc = twoln*cmk6*t5*temp*rjv1*ftinv*c2me
      go to 506
507   rnep  = twoln*cmk5*t5*rie2*ftinv
      snep  = twoln*cmk6*t5*temp*rjv2*ftinv
      snepc = twoln*cmk6*t5*temp*rjv2*ftinv*c2me
506   continue
      return
      end






      subroutine time_scales(tt,dd,taud,tau_nse,tau_qse)
      include 'implno.dek'

! input:
! tt = temperature
! dd = desnity
!
! output:
! taud    = e-folding timescale for density in an adiabatic expansion
! tau_nse = timescale to reach nse
! tau_qse = timescale to reach qse


! declare the pass
      double precision tt,dd,taud,tau_nse,tau_qse

! local variables
      double precision t9,tmin
      parameter        (tmin = 2.5d9)

! go
      if (tt .gt. tmin) then
       t9 = tt * 1.0d-9
       tau_nse = dd**(0.2d0) * exp(179.7d0/t9 - 40.5d0)
       tau_qse = exp(149.7d0/t9 - 39.15d0)
      else
       tau_nse = 1.0d20
       tau_qse = 1.0d20
      end if

      taud = 446.0d0/sqrt(dd)

      return
      end



      subroutine ener_gener_rate(dydt,enuc)
      include 'implno.dek'
      include 'const.dek'
      include 'network.dek'

! computes the instantaneous energy generation rate

! declare the pass
       real*8  :: dydt(*),enuc

! local variables
       integer :: i
       real*16 :: a1,a2,a3,xx

! conversion factors for the nuclear energy generation rate
! detlap is the mass excess of the proton in MeV
! detlan is the mass excess of the neutron in MeV

      real*8, parameter :: enuc_conv  = ev2erg*1.0d6*avo, &
                           enuc_conv2 = -avo*clight*clight, &
                           deltap     = 7.288969d0, &
                           deltan     = 8.071323d0


! this form misses n <-> p differences
!      enuc = 0.0d0
!      do i=1,ionmax
!       enuc = enuc + dydt(i) * bion(i)
!       enddo
!      enuc = enuc * enuc_conv


! this atomic mass excess form gets the n <-> p and e- and e+ parts 
!      enuc = 0.0d0
!      do i=1,ionmax
!       enuc = enuc + dydt(i) * (bion(i) - zion(i)*deltap - nion(i)*deltan)
!      enddo
!      enuc = enuc * enuc_conv


! this form is closest to e = m c**2 and also gets all the pieces
! minimize round-off by using quad and promotion rules

!      enuc = enuc_conv2 * sum(dydt(1:ionmax)*mion(1:ionmax))

      xx = 0.0q0
      do i=1,ionmax
       a1 = dydt(i) 
       a2 = mion(i)
       xx = xx + a1*a2
      enddo
      enuc = enuc_conv2 * xx



      return
      end






      subroutine read_partition_functions
      include 'implno.dek'
      include 'const.dek'
      include 'network.dek'
!
! this routine reads and initializes the partition function data
!
! declare

      character*132    :: string
      integer          :: i,j,izpf,iapf,inpf
      double precision :: xx


! initialize
      pf(1:pf_tmax,1:pf_zmax,1:pf_nmax) = -101.0d0 
      pf_gspin(1:pf_zmax,1:pf_nmax)     = -101.0d0

! open the partition function data file and read the header
      open(unit=22, file='partition_function_frdm.dat', status='old')
      do i=1,4
       read(22,11) string
11     format(a)
      enddo

! read the file
      do i=1,40000
       read(22,11,end=30) string
       if (string(1:12) .eq. 'END OF TABLE') goto 30
       read(22,*) izpf,iapf,xx
       inpf = iapf - izpf
       if (izpf .lt. 1) stop 'izpf < 1 during read of partition functions'
       if (inpf .lt. 1) stop 'inpf < 1 during read of partition functions'
       if (izpf .gt. pf_zmax) stop 'izpf > pf_izmax during read of partition functions'
       if (inpf .gt. pf_nmax) stop 'inpf > pf_inmax during read of partition functions'
       pf_gspin(izpf,inpf) = xx
       read(22,*) (pf(j,izpf,inpf), j=1,24)
      end do

! close up
30    close(unit=22)


! set the temperatures on which the partition functions are calculated
      pf_t9 = (/0.01d0, 0.15d0,  0.2d0, 0.3d0, 0.4d0, 0.5d0, 0.6d0, 0.7d0, &
                 0.8d0,  0.9d0,  1.0d0, 1.5d0, 2.0d0, 2.5d0, 3.0d0, 3.5d0, & 
                 4.0d0,  4.5d0,  5.0d0, 6.0d0, 7.0d0, 8.0d0, 9.0d0, 10.0d0/)

      return
      end




      subroutine get_partition_functions(t9)
      include 'implno.dek'
      include 'burn_common.dek'
      include 'vector_eos.dek'
      include 'network.dek'


! this routine generates partition function values

! declare the pass
     double precision :: t9


! declare locals
! for the partition function interpolation
      integer            :: i,jat
      integer, parameter :: norder=2
      double precision   :: x,x1,x2,a,b,e,alfa,beta,stat_weight


! waprt(i) is the temperature dependent partition functions
! zwork1(i) is the ratio of wpart(i) to the ground state partition function
! zwork2(i) is the temperature derivative of zwork1
! initialize


      zwork1(ionbeg:ionend) = 1.0d0
      zwork2(ionbeg:ionend) = 0.0d0
      wpart(ionbeg:ionend)  = 1.0d0

      if (ineut .ne. 0) then
       zwork1(ineut) = 1.0d0
       zwork2(ineut) = 0.0d0
       wpart(ineut)  = 2.0d0
      end if

      if (iprot .ne. 0) then
       zwork1(iprot) = 1.0d0
       zwork2(iprot) = 0.0d0
       wpart(iprot)  = 2.0d0
      end if

      if (ih1 .ne. 0) then
       zwork1(ih1) = 1.0d0
       zwork2(ih1) = 0.0d0
       wpart(ih1)  = 2.0d0
      end if

      if (ihe4 .ne. 0) then
       zwork1(ihe4) = 1.0d0
       zwork2(ihe4) = 0.0d0
       wpart(ihe4)  = 1.0d0
      end if



! for every isotope not named n, p or alpha

      do i=ionbeg,ionend

! bullet check
       if (int(zion(i)) .lt. 1) stop 'izpf < 1 during read of partition functions'
       if (int(nion(i)) .lt. 1) stop 'izpf < 1 during read of partition functions'
       if (int(zion(i)) .gt. pf_zmax) stop 'izpf > pf_izmax during read of partition functions'
       if (int(nion(i)) .gt. pf_nmax) stop 'inpf > pf_inmax during read of partition functions'


! if isotope is in the network 
       if (pf_gspin(int(zion(i)),int(nion(i))) .lt. 0.0d0) cycle

!       write(6,*) int(zion(i)),int(nion(i)),pf_gspin(int(zion(i)),int(nion(i)))


        stat_weight = 2.0d0 * pf_gspin(int(zion(i)),int(nion(i))) + 1.0d0

! clamp at both ends of the temperature grid
       if (t9 .le. pf_t9(1)) then
        zwork1(i) = pf(1,int(zion(i)),int(nion(i)))
        zwork2(i) = 0.0d0
        wpart(i)  = zwork1(i) * stat_weight

       else if (t9 .ge. pf_t9(pf_tmax)) then
        zwork1(i) = pf(pf_tmax,int(zion(i)),int(nion(i)))
        zwork2(i) = 0.0d0
        wpart(i)  = zwork1(i) * stat_weight

! linear larange interpolation

       else
        call locate(pf_t9,pf_tmax,t9,jat)
        jat = max(1,min(jat - norder/2 + 1,pf_tmax - norder + 1))

!        write(6,*) jat,t9
!        write(6,*) pf_t9(jat),pf_t9(jat+1)
!        write(6,*) pf(jat,int(zion(i)),int(nion(i))),pf(jat+1,int(zion(i)),int(nion(i)))

        x  = t9
        x1 = pf_t9(jat)
        x2 = pf_t9(jat+1)
        a  = x - x1
        b  = x - x2
        e  = x1 - x2
        alfa =  b/e
        beta = -a/e
        zwork1(i) = alfa * pf(jat,int(zion(i)),int(nion(i))) + beta * pf(jat+1,int(zion(i)),int(nion(i)))
        zwork2(i) = 0.0d0
        wpart(i)  = zwork1(i) * stat_weight
       end if

!       write(6,156) i, ionam(i), zwork1(i), wpart(i)
!156    format(1x,i4,' ',a,' ',1p2e12.4)
!       read(5,*)

      enddo
      return
      end







      subroutine init_isotope_rate_pointers
      include 'implno.dek'
      include 'network.dek'

! initializes all the integer pointers for the isotopes and rates


! zero all the isotope pointers
      ih1=0   ; iprot=0 ; ineut=0 ; ihe4=0  ; ih2=0   ; ih3=0   ; ihe3=0  ; ili6=0  ; ili7=0  ; ili8=0
      ibe7=0  ; ibe9=0  ; ib8=0   ; ib9=0   ; ib10=0  ; ib11=0  ; ic11=0  ; ic12=0  ; ic13=0  ; ic14=0
      in12=0  ; in13=0  ; in14=0  ; in15=0  ; io14=0  ; io15=0  ; io16=0  ; io17=0  ; io18=0  ; if17=0
      if18=0  ; if19=0  ; ine18=0 ; ine19=0 ; ine20=0 ; ine21=0 ; ine22=0 ; ine23=0 ; ina20=0 ; ina21=0
      ina22=0 ; ina23=0 ; img22=0 ; img23=0 ; img24=0 ; img25=0 ; img26=0 ; ial25=0 ; ial26=0 ; ial27=0
      ial28=0 ; isi26=0 ; isi27=0 ; isi28=0 ; isi29=0 ; isi30=0 ; ip28=0  ; ip29=0  ; ip30=0  ; ip31=0
      is30=0  ; is31=0  ; is32=0  ; is33=0  ; is34=0  ; icl35=0 ; icl37=0 ; iar36=0 ; iar37=0 ; iar38=0
      ik39=0  ; ica40=0 ; isc43=0 ; isc45=0 ; iti44=0 ; iti46=0 ; iti48=0 ; iti50=0 ; iv46=0  ; iv47=0
      iv48=0  ; iv51=0  ; icr47=0 ; icr48=0 ; icr49=0 ; icr50=0 ; icr51=0 ; icr52=0 ; icr53=0 ; icr54=0
      icr56=0 ; imn50=0 ; imn51=0 ; imn52=0 ; imn55=0 ; ife52=0 ; ife54=0 ; ife55=0 ; ife56=0 ; ife57=0
      ife58=0 ; ife59=0 ; ico54=0 ; ico55=0 ; ico56=0 ; ico59=0 ; ini56=0 ; ini57=0 ; ini58=0 ; ini59=0
      ini64=0 ; ini66=0 ; icu63=0 ; izn60=0 ; izn64=0


! zero all the rate pointers

      ir3a=0   ;  irg3a=0  ;  ircag=0  ;  ir1212=0 ;  ir1216=0 ;  iroga=0 
      iroag=0  ;  irnega=0 ;  irneag=0 ;  irmgga=0 ;  irmgag=0 ;  irsiga=0   
      irmgap=0 ;  iralpa=0 ;  iralpg=0 ;  irsigp=0 ;  irsiag=0 ;  irsga=0    
      irsiap=0 ;  irppa=0  ;  irppg=0  ;  irsgp=0  ;  irsag=0  ;  irarga=0   
      irsap=0  ;  irclpa=0 ;  irclpg=0 ;  irargp=0 ;  irarag=0 ;  ircaga=0   
      irarap=0 ;  irkpa=0  ;  irkpg=0  ;  ircagp=0 ;  ircaag=0 ;  irtiga=0   
      ircaap=0 ;  irscpa=0 ;  irscpg=0 ;  irtigp=0 ;  irtiag=0 ;  ircrga=0   
      irtiap=0 ;  irvpa=0  ;  irvpg=0  ;  ircrgp=0 ;  ircrag=0 ;  irfega=0   
      ircrap=0 ;  irmnpa=0 ;  irmnpg=0 ;  irfegp=0 ;  irfeag=0 ;  irniga=0   
      ir1616=0 ;  iraan=0

      irpp=0   ;  ir33=0   ;  ircpg=0  ;  irnpg=0  ;  iropg=0   ;  irnag=0 
      irfeap=0 ;  ircopa=0 ;  ircopg=0 ;  irnigp=0 ;  irfepg=0  ;  ircogp=0 
      ir52ng=0 ;  ir53gn=0 ;  ir53ng=0 ;  ir54gn=0 ;  irheng=0  ;  irhegn=0   
      irhng=0  ;  irhegp=0 ;  irpen=0  ;  irnep=0  ;  irn56ec=0 ;  ifa=0      
      ifg=0    

      ir54ng=0      ; ir55gn=0       ; ir55ng=0      ; ir56gn=0    
      irfe54ap=0    ; irco57pa=0     ; irfe56pg=0    ; irco57gp=0 
      irfe56_aux1=0 ;  irfe56_aux2=0 ; irfe56_aux3=0 ; irfe56_aux4=0

      irbepg=0 ;  irb8ep=0 ;  irbeec=0 ;  irpep=0 ;  irhep=0

      irc12pg=0  ;  irn13gp=0  ;  irn13enu=0  ;  irc13pg=0   ; irn14gp=0 
      irn14pg=0  ;  iro15gp=0  ;  iro15enu=0  ;  irn15pa=0   ; irc12ap=0    
      irn15pg=0  ;  iro16gp=0  ;  iro16pg=0   ;  irf17gp=0   ; irf17enu=0   
      iro17pa=0  ;  irn14ap=0  ;  iro17pg=0   ;  irf18gp=0   ; irf18enu=0   
      iro18pa=0  ;  irn15ap=0  ;  iro18pg=0   ;  irf19gp=0   ; irf19pa=0    
      iro16ap=0  ;  irn13pg=0  ;  iro14gp=0   ;  iro14enu=0  ; iro14ap=0    
      irf17pa=0  ;  irf17pg=0  ;  irne18gp=0  ;  irne18enu=0 ; irf18pa=0    
      iro15ap=0  ;  irne18ap=0 ;  iro15ag=0   ;  irne19ga=0  ; irne19pg=0   
      irna20gp=0 ;  irsi26ap=0 ;  irne19enu=0 ;  irne20pg=0  ; irna21gp=0

      irweakn=0 ;     irweakp=0 ;     irpng=0 ;       irdgn=0 ;      irdpn=0 
      ir2pnp=0 ;      irdng=0 ;       irtgn=0 ;       irddp=0 ;      irtpd=0       
      irdpg=0 ;       irhe3gp=0 ;     irddn=0 ;       irhe3nd=0 ;    irtpn=0       
      irhe3np=0 ;     irddg=0 ;       irhe4gd=0 ;     irtpg=0 ;      irhe4gp=0 
      irtdn=0 ;       irhe4nd=0 ;     irtt2n=0 ;      irhe42nt=0 ;   irhe3ng=0     
      irhe4gn=0 ;     irhe3dp=0 ;     irhe4pd=0 ;     irhe3td=0 ;    irhe4dt=0     
      irhe3tnp=0 ;    irhe4pnt=0 ;    ir33inv=0 ;     irhe4npg=0 ;   irli6gpn=0    
      irdag=0 ;       irli6ga=0 ;     irtan=0 ;       irli6na=0 ;    irli6pa=0     
      irhe3ap=0 ;     irtag=0 ;       irli7ga=0 ;     irli7pa=0 ;    ir2he4ga=0 
      irli7dn=0 ;     ir2and=0 ;      irli7t2n=0 ;    ir2a2nt=0 ;    irli7he3np=0 
      ir2apnhe3=0 ;   irhe3ag=0 ;     irbe7ga=0 ;     irli6pg=0 ;    irbe7gp=0   
      irli7pn=0 ;     irbe7np=0 ;     irbe7dp=0 ;     ir2apd=0 ;     irbe7tnp=0    
      ir2apnt=0 ;     irbe7he32p=0 ;  ir2a2phe3=0 ;   irbe9pa=0 ;    irli6ap=0     
      irbe7pg=0 ;     irb8gp=0 ;      irli6ag=0 ;     irb10ga=0 ;    irli7an=0     
      irb10na=0 ;     irbe9pg=0 ;     irb10gp=0 ;     irb10pa=0 ;    irbe7ap=0     
      irli7ag=0 ;     irb11ga=0 ;     irb11pa=0 ;     ir2he4ap=0 ;   irbe7ag=0     
      irc11ga=0 ;     irb11pn=0 ;     irc11np=0 ;     irb8ap=0 ;     irc11pa=0      
      irb10pg=0 ;     irc11gp=0 ;     irc11na=0 ;     ir3an=0 ;      irbe9an=0     
      irc12na=0 ;     irb11pg=0 ;     irc12gp=0 ;     irb11ap=0 ;    irc14pa=0     
      irli6ng=0 ;     irli7gn=0


      irh3em=0 ;    irc11em=0 ;   iraanbe7=0 ;  irdap=0 ;     irtap=0 ;   irb10ng=0 
      irbe9dn=0 ;   irb10dp=0 ;   irb8npa=0 ;   irbe9pda=0 ;  ir2adp=0 ;  irb11p2a=0 
      irc11n2a=0


      irh3pg=0 ;    irli7pag=0 ;  irgaan=0 ;    irbe9pd=0 ;   ir3ap=0 ;    ir1212p=0 ;   irna23pc=0 
      ir1212n=0 ;   irmg23nc=0 ;  ir1212a=0 ;   irne20ac=0 ;  ir1616n=0 ;  irs31no=0 ;   ir1616p=0    
      irp31po=0 ;   ir1616a=0 ;   irsi28ao=0 ;  ir1616d=0 ;   irp30do=0 ;  ir1216n=0 ;   irsi27no=0   
      ir1216p=0 ;   iral27po=0 ;  ir1216a=0 ;   irmg24ao=0 ;  irli7ng=0 ;  irli8gn=0 ;   irbe9pn=0    
      irb9np=0 ;    irb10an=0 ;   irn13na=0 ;   irb11an=0 ;   irn14na=0 ;  irn14pa=0 ;   irc11ap=0    
      irc11pg=0 ;   irn12gp=0 ;   irc12an=0 ;   iro15na=0 ;   irc13an=0 ;  iro16na=0 ;   irc13pn=0    
      irn13np=0 ;   irc14pg=0 ;   irn15gp=0 ;   irc14pn=0 ;   irn14np=0 ;  irc14ag=0 ;   iro18ga=0    
      iro16pa=0 ;   irn13ap=0 ;   irn14pn=0 ;   iro14np=0 ;   irn14an=0 ;  irf17na=0 ;   irn14ag=0    
      irf18ga=0 ;   irn15pn=0 ;   iro15np=0 ;   irn15an=0 ;   irf18na=0 ;  irn15ag=0 ;   irf19ga=0    
      iro14ag=0 ;   irne18ga=0 ;  iro17ag=0 ;   irne21ga=0 ;  iro17an=0 ;  irne20na=0 ;  iro18ag=0    
      irne22ga=0 ;  iro18an=0 ;   irne21na=0

      irne20pa=0 ;  irf17ap=0 ;   irf18pg=0 ;   irne19gp=0 ;  irf19pg=0 ;   irne20gp=0 ;  irf19pn=0    
      irne19np=0 ;  irf19ap=0 ;   irne22pa=0 ;  irna22na=0 ;  irf19an=0 ;   irna23pa=0 ;  irne20ap=0    
      irne20ng=0 ;  irne21gn=0 ;  irne21pg=0 ;  irna22gp=0 ;  irne21ag=0 ;  irmg25ga=0 ;  irne21an=0   
      irmg24na=0 ;  irne22pg=0 ;  irna23gp=0 ;  irne22ag=0 ;  irmg26ga=0 ;  irna22np=0 ;  irne22pn=0   
      irne22an=0 ;  irmg25na=0 ;  irna21pg=0 ;  irmg22gp=0 ;  irmg24pa=0 ;  irna21ap=0 ;  irna22pg=0   
      irmg23gp=0 ;  irna23pg=0 ;  irmg24gp=0 ;  irna23pn=0 ;  irmg23np=0 ;  irmg24pg=0 ;  iral25gp=0   
      iral27pa=0 ;  irmg24ap=0 ;  irmg25pg=0 ;  iral26gp=0 ;  irmg25ap=0 ;  iral28pa=0 ;  irmg25ag=0   
      irsi29ga=0 ;  irmg25an=0 ;  irsi28na=0 ;  irmg26pg=0 ;  iral27gp=0 ;  irmg26ag=0 ;  irsi30ga=0   
      irmg26an=0 ;  irsi29na=0 ;  iral25pg=0 ;  irsi26gp=0 ;  iral26pg=0 ;  irsi27gp=0 ;  iral27pg=0   
      irsi28gp=0 ;  iral27an=0 ;  irp30na=0 ;   irsi27pg=0 ;  irp28gp=0 ;   irsi28pg=0 ;  irp29gp=0    
      irsi29pg=0 ;  irp30gp=0 ;   irsi30pg=0 ;  irp31gp=0 ;   irr1=0 ;      irs1=0 ;      irt1=0       
      iru1=0 ;      irv1=0 ;      irw1=0 ;      irx1=0 ;      iry1=0 ;      ir1f54=0 ;    ir2f54=0     
      ir3f54=0 ;    ir4f54=0 ;    ir5f54=0 ;    ir6f54=0 ;    ir7f54=0 ;    ir8f54=0 ;    iralf1=0     
      iralf2=0 ;    irsi2ni=0 ;   irni2si=0 ;   iralam1=0 ;   irdelta1=0 ;  iralam2=0 ;   irdelta2=0 
      irdebug=0

      return
      end








      subroutine read_abundances_from_file(xin)
      include 'implno.dek'
      include 'network.dek'

! this routine reads a flat file of isotope names and mass fractions 

! declare the pass
      real*8            :: xin(*)


! for reading the abundance file
      character*80      :: string,word,atoken,aval
      logical           :: err
      integer           :: i,j,k,ip,getnam
      real*8            :: x,value2


! popular format statements
01    format(1x,a,a,a)
02    format(a)


! get file name
10    write(6,*) 'give network file name =>'
      read(5,2) string
      open(unit=2, file=string(1:len_trim(string)), status='old', err=10)


! for every line in the file
      do i=1,10000
       read(2,2,end=40) string

! replace any tab with a space
       j = len_trim(string)
       do k=1,j
        if (iachar(string(k:k)) .eq. 9) string(k:k) = ' '
       enddo

! strip leading blanks and any comments, which are preceeded by a !
       word = adjustl(string)
       j = len_trim(word)
       do k=1,j
        if (iachar(word(k:k)) .eq. 33) then
         do ip=k,j
          word(ip:ip) = ' '
         end do
        end if
       end do

! get another line if this one is now empty 
       j = len_trim(word)
       if (j .eq. 0) cycle

!       write(6,*) trim(word)
!       read(5,*)

! chop the string to get the name and value
       ip = 1
       j = getnam(word,atoken,ip)
       j = getnam(word,aval,ip)
       x = value2(aval,err)
       if (err) then
        write(6,*)
        write(6,*) 'could not convert ',trim(aval),' to a number, skipping'
        write(6,*)
        goto 30
       end if

!       write(6,*) trim(atoken(1:5)),x
!       read(5,*)

! see if the isotope is in the network
       do j=1,ionmax
        if (ionam(j)(1:5) .eq. atoken(1:5)) then
         xin(j) = x
!         write(6,*) 'found it'
         goto 30
        end if
       enddo

! isotope not in network, check special cases
       if (atoken(1:5) .eq. 'n   ') then
        xin(ineut) = x
        goto 30
       else if (atoken(1:5) .eq. 'h1   ') then
        if (ih1 .ne. 0) xin(ih1) = x
        goto 30
       else if (atoken(1:5) .eq. 'p    ') then
        xin(iprot) = x
        goto 30
       else if (atoken(1:5) .eq. 'd    ') then
        if (ih2 .ne. 0) xin(ih2) = x
        goto 30
       else if (atoken(1:5) .eq. 't    ') then
        if (ih3 .ne. 0) xin(ih3) = x
        goto 30
       end if

! for sure not in the network list
       write(6,*)
       write(6,20) atoken(1:5),x
20     format(1x,'did not find isotope in network: ', a,' ',1pe13.5,' ,skipping it.')
       write(6,*)


! back for another line in the file 
30     continue
      enddo


! close file and return
40    close(unit=2)
      return 
      end










      subroutine prepare_for_net_eval(tt,y,abar,zbar,ye)
      include 'implno.dek'
      include 'const.dek'
      include 'burn_common.dek'
      include 'network.dek'
      include 'vector_eos.dek'

! does some prep work before evaluating the right-hand-sides and jacobian evaluations
! for a network integration 

! declare the pass
      real*8  :: tt,y(*),abar,zbar,ye


! positive definite molar fractions
      y(1:ionmax) = min(1.0d0,max(y(1:ionmax),1.0d-30))


! generate abar and zbar for this composition
      call azbar_simple(y,aion,zion,ionbeg,ionmax,abar,zbar,ye)


! for evolution with energy, temp, den equations
      if (pure_network .eq. 0) then

! temperature and density from a history file
       if (trho_hist) then
        call update2(tt,y(itemp),y(iden))
        y(itemp) = y(itemp) * inv_temp_scale
        y(iden)  = y(iden)  * inv_den_scale
       end if

   
! for pressure evolutions
       if (self_heat_const_pres .or. pt_hist) then
        jlo_eos = 1
        jhi_eos = 1
        ptot_row(1) = bpres
        temp_row(1) = y(itemp)*temp_scale
        abar_row(1) = abar
        zbar_row(1) = zbar
        den_row(1)  = y(iden)*den_scale
        call invert_helm_pt_quiet
        y(iden) = den_row(1) * inv_den_scale
       end if

! bounded temperatures and densities
       y(itemp) = min(1.0d11,max(y(itemp)*temp_scale,1.0d4)) * inv_temp_scale
       y(iden)  = min(1.0d11,max(y(iden)*den_scale,1.0d-10)) * inv_den_scale

! set the common block temperature and density
       btemp = y(itemp)*temp_scale
       bden  = y(iden)*den_scale


! or for pure network
      else if (pure_network .eq. 1) then

! get the new temperature and density if need be from a history file
       if (trho_hist) call update2(tt,btemp,bden)

      end if
      return
      end subroutine prepare_for_net_eval

!---------------------------------------------------------------------







!---------------------------------------------------------------------
! function wien1
! function dwien1dx
! function wien2
! function dwien2dx
! function func1
! function dfunc1dx
! function func2
! function dfunc2dx
! function bb_qromb
! function bb_trapzd
! routine bb_polint does polynomial interpolation





      double precision function wien1(x)
      include 'implno.dek'
      include 'const.dek'

! this is the function given in
! weinberg's "gravitation and cosmology" page 537, equation 15.6.40

! declare the pass
      double precision x

! communicate xcom via common block
      double precision xcom
      common /tes1/    xcom

! local variables
      external         func1
      double precision func1,f1,con


! the integration limits ylo and yhi, along with the integration
! tolerance tol, will give at least 9 significant figures of precision

      double precision ylo,yhi,tol
      parameter        (ylo = 1.0d-6, &
                        yhi = 50.0d0, &
                        tol = 1.0d-10, &
                        con = 45.0d0/(2.0d0*pi*pi*pi*pi))


! for quadrature
      integer          nquad,ifirst
      parameter        (nquad = 100)
      double precision xquad(nquad),wquad(nquad)
      data             ifirst/0/



! initialization of the quadrature abcissas and weights
      if (ifirst .eq. 0) then
       ifirst = 1
       call bb_gauleg(ylo,yhi,xquad,wquad,nquad)
      end if


! don't do any integration if x is large enough
      if (x .gt. 50.0) then
       wien1 = 1.0d0

! do the integration
      else
       xcom = x
!       call bb_qromb(func1,ylo,yhi,tol,f1)
       call bb_qgaus(func1,xquad,wquad,nquad,f1)
       wien1 = 1.0d0 + con * f1
      end if

      return
      end





      double precision function dwien1dx(x)
      include 'implno.dek'
      include 'const.dek'

! this is the derivative with respect to x of the function given in
! weinberg's "gravitation and cosmology" page 537, equation 15.6.40

! declare the pass
      double precision x

! communicate xcom via common block
      double precision xcom
      common /tes1/    xcom

! local variables
      external         dfunc1dx
      double precision dfunc1dx,df1,con


! the integration limits ylo and yhi, along with the integration
! tolerance tol, will give at least 9 significant figures of precision

      double precision ylo,yhi,tol
      parameter        (ylo = 1.0d-6, &
                        yhi = 50.0d0, &
                        tol = 1.0d-10, &
                        con = 45.0d0/(2.0d0*pi*pi*pi*pi))


! for quadrature
      integer          nquad,ifirst
      parameter        (nquad = 100)
      double precision xquad(nquad),wquad(nquad)
      data             ifirst/0/


! initialization of the quadrature abcissas and weights
      if (ifirst .eq. 0) then
       ifirst = 1
       call bb_gauleg(ylo,yhi,xquad,wquad,nquad)
      end if


! don't do any integration if x is large enough
      if (x .gt. 50.0) then
       dwien1dx = 0.0d0

! do the integration
      else
       xcom = x
!       call bb_qromb(dfunc1dx,ylo,yhi,tol,df1)
       call bb_qgaus(dfunc1dx,xquad,wquad,nquad,df1)
       dwien1dx = con * df1
      end if

      return
      end







      double precision function wien2(x)
      include 'implno.dek'
      include 'const.dek'

! this is the function given in
! weinberg's "gravitation and cosmology" page 537, equation 15.6.40

! declare the pass
      double precision x


! communicate xcom via common block
      double precision xcom
      common /tes1/    xcom


! communicate the number of neutrino families
! using 2 families of neutrinos duplicates the time-temperature
! table in wienberg's "gravitation and cosmology", page 540, table 15.4

      double precision xnnu
      common /nufam/   xnnu


! local variables
      external         func2
      double precision func2,f2,wien1


! the integration limits ylo and yhi, along with the integration
! tolerance tol, will give at least 9 significant figures of precision

      double precision ylo,yhi,tol,con1,con2,con3,fthirds
      parameter        (ylo     = 1.0d-6, &
                        yhi     = 50.0d0, &
                        tol     = 1.0d-10, &
                        con2    = 4.0d0/11.0d0, &
                        con3    = 30.0d0/(pi*pi*pi*pi), &
                        fthirds = 4.0d0/3.0d0)


! for quadrature
      integer          nquad,ifirst
      parameter        (nquad = 100)
      double precision xquad(nquad),wquad(nquad)
      data             ifirst/0/


! initialization of the quadrature abcissas and weights
      if (ifirst .eq. 0) then
       ifirst = 1
       call bb_gauleg(ylo,yhi,xquad,wquad,nquad)

! a constant that depends on the number of neutrino families
       con1 = xnnu * 7.0d0/8.0d0
      end if



! don't do any integration if x is large enough
      if (x .gt. 50.0) then
       wien2 = 1.0d0 + con1*con2**fthirds

! do the integration
      else
       xcom = x
!       call bb_qromb(func2,ylo,yhi,tol,f2)
       call bb_qgaus(func2,xquad,wquad,nquad,f2)
       wien2 = 1.0d0 + con1 * (con2 * wien1(x))**fthirds + con3 * f2
      end if

      return
      end






      double precision function dwien2dx(x)
      include 'implno.dek'
      include 'const.dek'

! this is the derivative with respect to x of the function given in
! weinberg's "gravitation and cosmology" page 537, equation 15.6.40

! declare the pass
      double precision x



! communicate xcom via common block
      double precision xcom
      common /tes1/    xcom


! xnnu communicates the number of neutrino families
! using 2 families of neutrinos duplicates the time-temperature
! table in wienberg's "gravitation and cosmology", page 540, table 15.4

      double precision xnnu
      common /nufam/   xnnu


! local variables
      external         dfunc2dx
      double precision dfunc2dx,w1,dw1,f2,df2,wien1,dwien1dx


! the integration limits ylo and yhi, along with the integration
! tolerance tol, will give at least 9 significant figures of precision

      double precision ylo,yhi,tol,con1,con2,con3,fthirds,third
      parameter        (ylo     = 1.0d-6, &
                        yhi     = 50.0d0, &
                        tol     = 1.0d-10, &
                        con2    = 4.0d0/11.0d0, &
                        con3    = 30.0d0/(pi*pi*pi*pi), &
                        fthirds = 4.0d0/3.0d0, &
                        third   = 1.0d0/3.0d0)


! for quadrature
      integer          nquad,ifirst
      parameter        (nquad = 100)
      double precision xquad(nquad),wquad(nquad)
      data             ifirst/0/


! initialization of the quadrature abcissas and weights
      if (ifirst .eq. 0) then
       ifirst = 1
       call bb_gauleg(ylo,yhi,xquad,wquad,nquad)

! a constant that depends on the number of neutrino families
       con1 = xnnu * 7.0d0/8.0d0
      end if



! don't do any integration if x is large enough
      if (x .gt. 50.0) then
       dwien2dx = 0.0d0

! do the integration
      else
       xcom = x
!       call bb_qromb(dfunc2dx,ylo,yhi,tol,df2)
       call bb_qgaus(dfunc2dx,xquad,wquad,nquad,df2)
       w1   = wien1(x)
       dw1  = dwien1dx(x)
!       w2   = 1.0d0 + con1*(con2*w1)**fthirds + con3*f2
       dwien2dx = fthirds*con1*(con2*w1)**third * con2*dw1 + con3*df2
      end if

      return
      end






      double precision function func1(y)
      include 'implno.dek'
      include 'const.dek'

! this is the integrand of the function given in
! weinberg's "gravitation and cosmology" page 537, equation 15.6.40

! declare the pass
      double precision y

! local variables
      double precision y2,x2,aa,aalim
      parameter        (aalim = 200.0d0)


! communicate xcom via common block
      double precision xcom
      common /tes1/    xcom

      func1 = 0.0d0
      if (aa .le. aalim) then
       y2    = y * y
       x2    = xcom * xcom
       aa    = sqrt(y2 + x2)
       func1 = (aa + y2/(3.0d0*aa)) * y2 / (exp(aa) + 1.0d0)
      end if
      return
      end






      double precision function dfunc1dx(y)
      include 'implno.dek'
      include 'const.dek'

! this is the derivative with respect to x of the integrand in the function
! given by weinberg's "gravitation and cosmology" page 537, equation 15.6.40

! declare the pass
      double precision y

! local variables
      double precision y2,x2,aa,daa,zz,denom,ddenom,f1,aalim
      parameter        (aalim = 200.0d0)


! communicate xcom via common block
      double precision xcom
      common /tes1/    xcom

      dfunc1dx = 0.0d0
      if (aa .le. aalim) then
       y2       = y * y
       x2       = xcom * xcom
       aa       = sqrt(y2 + x2)
       daa      = xcom/aa
       zz       = exp(aa)
       denom    = zz + 1.0d0
       ddenom   = zz * daa
       f1       = (aa + y2/(3.0d0*aa)) * y2 / denom
       dfunc1dx = (1.0d0 - y2/(3.0d0*aa**2)) * daa * y2 / denom &
                  - f1/denom * ddenom
      end if
      return
      end





      double precision function func2(y)
      include 'implno.dek'
      include 'const.dek'

! this is the integrand of the function given in
! weinberg's "gravitation and cosmology" page 539, equation 15.6.48

! declare the pass
      double precision y

! local variables
      double precision y2,x2,aa,aalim
      parameter        (aalim = 200.0d0)


! communicate xcom via common block
      double precision xcom
      common /tes1/    xcom

      func2 = 0.0d0
      if (aa .le. aalim) then
       y2    = y * y
       x2    = xcom * xcom
       aa    = sqrt(y2 + x2)
       func2 = aa * y2 / (exp(aa) + 1.0d0)
      end if
      return
      end






      double precision function dfunc2dx(y)
      include 'implno.dek'
      include 'const.dek'

! this is the derivative with respect to x of the integrand
! of the function given in weinberg's "gravitation and cosmology"
! page 539, equation 15.6.48

! declare the pass
      double precision y

! local variables
      double precision y2,x2,aa,daa,zz,denom,ddenom,f2,aalim
      parameter        (aalim = 200.0d0)


! communicate xcom via common block
      double precision xcom
      common /tes1/    xcom

      dfunc2dx = 0.0d0
      if (aa .le. aalim) then
       y2       = y * y
       x2       = xcom * xcom
       aa       = sqrt(y2 + x2)
       daa      = xcom/aa
       zz       = exp(aa)
       denom    = zz + 1.0d0
       ddenom   = zz * daa
       f2       = aa * y2 / denom
       dfunc2dx = (daa*y2 - f2*ddenom)/denom
      end if
      return
      end












      subroutine bb_qromb(func,a,b,eps,ss)
      include 'implno.dek'

! returns as ss the integral of the function func from a to b with fractional
! accuracy eps. integration by romberg's method of order 2k where e.g k=2 is
! simpson's rule.
!
! jmax limits the total number of steps; k is the
! the number of points used in the extrapolation; arrays s and h store the
! trapazoidal approximations and their relative step sizes.

! declare
      external          func
      integer           j,jmax,jmaxp,k,km
      parameter         (jmax=20, jmaxp=jmax+1, k=5, km=k-1)
      double precision  a,b,ss,s(jmaxp),h(jmaxp),eps,dss,func

      h(1) = 1.0d0
      do j=1,jmax
       call bb_trapzd(func,a,b,s(j),j)
       if (j .ge. k) then
        call bb_polint(h(j-km),s(j-km),k,0.0d0,ss,dss)
        if (abs(dss) .le. eps*abs(ss)) return
       end if
       s(j+1) = s(j)
       h(j+1) = 0.25d0 * h(j)
      enddo

!      write(6,*) ' after ',jmax,' iterations '
!      write(6,*) ' of trying to integrate between ',a,' and ',b
!      write(6,*) ' and fractional accuracy ',eps
!      write(6,*) ' the integral is ',ss
!      write(6,*) ' and error estimate ',dss
!      write(6,*) ' so that abs(dss) ',abs(dss),
!     1           ' > eps*abs(ss)',eps*abs(ss)
      ss = 0.0d0
!      stop       'too many steps in qromb'
      return
      end





      subroutine bb_trapzd(func,a,b,s,n)
      include 'implno.dek'

! this routine computes the n'th stage of refinement of an extended
! trapazoidal rule. func is input as the name of a function to be
! integrated between limits a and b. when n=1 the routine returns as s
! the crudest estimate of the integral of func(x)dx from a to b.
! subsequent calls with n=2,3... will improve the accuracy of s by adding
! 2**(n-2) additional interior points. s should not be modified between
! sequential calls.
!
! this routine is the workhorse of all the following closed formula
! integration routines.
!
! local it  is the number of points to be added on the next call
! local del is the step size.
!
! declare
      external          func
      integer           n,it,j
      double precision  func,a,b,s,del,x,sum,tnm

! go
      if (n.eq.1) then
       s  = 0.5d0 * (b-a) * ( func(a) + func(b) )
      else
       it  = 2**(n-2)
       tnm = it
       del = (b-a)/tnm
       x   = a + (0.5d0 *del)
       sum = 0.0d0
       do j=1,it
        sum = sum + func(x)
        x   = x + del
       enddo
       s  = 0.5d0 * (s + (b-a)*sum/tnm)
      end if
      return
      end







      subroutine bb_polint(xa,ya,n,x,y,dy)
      include 'implno.dek'

! given arrays xa and ya of length n and a value x, this routine returns a
! value y and an error estimate dy. if p(x) is the polynomial of degree n-1
! such that ya = p(xa) then the returned value is y = p(x)

! declare
      integer          n,nmax,ns,i,m
      parameter        (nmax=10)
      double precision xa(n),ya(n),x,y,dy,c(nmax),d(nmax),dif,dift, &
                       ho,hp,w,den


! find the index ns of the closest table entry; initialize the c and d tables
      ns  = 1
      dif = abs(x - xa(1))
      do i=1,n
       dift = abs(x - xa(i))
       if (dift .lt. dif) then
        ns  = i
        dif = dift
       end if
       c(i)  = ya(i)
       d(i)  = ya(i)
      enddo

! first guess for y
      y = ya(ns)

! for each column of the table, loop over the c's and d's and update them
      ns = ns - 1
      do m=1,n-1
       do i=1,n-m
        ho   = xa(i) - x
        hp   = xa(i+m) - x
        w    = c(i+1) - d(i)
        den  = ho - hp
        if (den .eq. 0.0) stop ' 2 xa entries are the same in polint'
        den  = w/den
        d(i) = hp * den
        c(i) = ho * den
       enddo

! after each column is completed, decide which correction c or d, to add
! to the accumulating value of y, that is, which path to take in the table
! by forking up or down. ns is updated as we go to keep track of where we
! are. the last dy added is the error indicator.
       if (2*ns .lt. n-m) then
        dy = c(ns+1)
       else
        dy = d(ns)
        ns = ns - 1
       end if
       y = y + dy
      enddo
      return
      end





      subroutine bb_qgaus(func,x,w,n,ss)
      include 'implno.dek'

! returns as ss the quadrature summation of the function func as determined
! by the abcissas x and weights w.

! declare
      external          func
      integer           i,n
      double precision  func,ss,x(n),w(n)

      ss = 0.0d0
      do i=1,n
       ss = ss + w(i)*func(x(i))
      enddo
      return
      end





      subroutine bb_gauleg(x1,x2,x,w,n)
      include 'implno.dek'
      include 'const.dek'

! given the lower and upper limits of integration x1 and x2, and given n,
! this routine returns arrays x and w of length n, containing the
! abscissas and weights of the gauss-legendre n-point quadrature formula

! declare
      integer            i,m,j,n
      double precision   x1,x2,x(n),w(n),eps,xm,xl,p1,p2,p3,pp,z,z1
      parameter          (eps=1.0e-14)


! roots are symmetric in the interval so we only have to find half of them
      m = (n+1)/2
      xm = 0.5d0 * (x2 + x1)
      xl = 0.5d0 * (x2 - x1)

! loop over the desired roots and make a slick guess at each one
      do i=1,m
       z = cos(3.141592653589d0 * (i-0.25d0)/(n + 0.5d0))

! newton do while loop
1      continue
       p1 = 1.0d0
       p2 = 0.0d0

! loop the recurrence relation to get the legendre polynomial at z
       do j=1,n
        p3 = p2
        p2 = p1
        p1 = ((2.0d0 * j - 1.0d0) * z * p2 - (j - 1.0d0)*p3)/j
       enddo

! p1 is now the desired legendre polynomial. pp is the derivative.
       pp = n * (z*p1 - p2)/(z*z - 1.0d0)
       z1 = z
       z  = z1 - p1/pp
       if (abs(z-z1) .gt. eps) goto  1

! scale to the users interval
       x(i)     = xm - xl*z
       x(n+1-i) = xm + xl * z
       w(i)     = 2.0d0 * xl/((1.0d0 - z*z)*pp*pp)
       w(n+1-i) = w(i)
      enddo
      return
      end


!---------------------------------------------------------------------







!---------------------------------------------------------------------
      subroutine net_decay_abund(xout)
      include 'implno.dek'
      include 'burn_common.dek'
      include 'network.dek'

! decays the composition

! declare the pass
      double precision xout(*)

! local variables
      character*80     decayed
      integer          i,lenstr

      integer          nsol
      parameter        (nsol = 286)
      character*5      namsol(nsol)
      integer          nzsol(nsol),nasol(nsol)
      double precision xstable(nsol),ag(nsol)



! popular format statements
 01   format(a,'decayed.dat')
 02   format(1x,i4,i4,1p2e12.4,a6)



! for the file name and open it
      write(decayed,01) hfile(1:lenstr(hfile,80))
      call sqeeze(decayed)
      open(unit=51,file=decayed,status='unknown')


! convert to integers
      do i=1,ionmax
       izwork1(i) = int(zion(i))
       izwork2(i) = int(aion(i))
      enddo

! do the work
      call decay_andgrev(ionmax,izwork1,izwork2,xout, &
                         nsol,namsol,nzsol,nasol,xstable,ag)

!      call decay_lodders(ionmax,izwork1,izwork2,xout,
!     1                   nsol,namsol,nzsol,nasol,xstable,ag)


! write it out
      write(51,02) (nzsol(i), nasol(i), xstable(i), &
                     ag(i), namsol(i), i=1,nsol)


! close up shop
      close(unit=51)
      return
      end


!---------------------------------------------------------------------







!---------------------------------------------------------------------
! this file contains abundance routines
! function andgrev returns the solar abundance of an isotope or element
! routine decay_andgrev reduces an abundance vector to the stable isotopes


      double precision function andgrev(nam,z,a,xelem)
      include 'implno.dek'

! anders and grevesse 1989 solar abundances from h1 to u238
! or the lodders 2003 solar abundances from h1 to u238

! input:
! name of the isotope nam

! output
! mass fraction andgrev
! charge z
! number of nucleons a
! elemental mass fraction associated with this isotope xelem

! declare the pass
      character*(*)    nam
      double precision z,a,xelem

! for the solar abundance data
      integer          solsiz
      parameter        (solsiz = 286)
      character*5      namsol(solsiz)
      integer          izsol(solsiz),iasol(solsiz),jcode(solsiz)
      double precision sol(solsiz)

! local variables
      integer          i,j,ifirst,jbeg,jend
      double precision sum,zsol,yesol

      data ifirst/0/


! bring in the solar abundance data

! names of the stable isotopes
      data (namsol(j), j=1,120) / &
       'h1   ','h2   ','he3  ','he4  ','li6  ','li7  ','be9  ','b10  ', &
       'b11  ','c12  ','c13  ','n14  ','n15  ','o16  ','o17  ','o18  ', &
       'f19  ','ne20 ','ne21 ','ne22 ','na23 ','mg24 ','mg25 ','mg26 ', &
       'al27 ','si28 ','si29 ','si30 ','p31  ','s32  ','s33  ','s34  ', &
       's36  ','cl35 ','cl37 ','ar36 ','ar38 ','ar40 ','k39  ','k40  ', &
       'k41  ','ca40 ','ca42 ','ca43 ','ca44 ','ca46 ','ca48 ','sc45 ', &
       'ti46 ','ti47 ','ti48 ','ti49 ','ti50 ','v50  ','v51  ','cr50 ', &
       'cr52 ','cr53 ','cr54 ','mn55 ','fe54 ','fe56 ','fe57 ','fe58 ', &
       'co59 ','ni58 ','ni60 ','ni61 ','ni62 ','ni64 ','cu63 ','cu65 ', &
       'zn64 ','zn66 ','zn67 ','zn68 ','zn70 ','ga69 ','ga71 ','ge70 ', &
       'ge72 ','ge73 ','ge74 ','ge76 ','as75 ','se74 ','se76 ','se77 ', &
       'se78 ','se80 ','se82 ','br79 ','br81 ','kr78 ','kr80 ','kr82 ', &
       'kr83 ','kr84 ','kr86 ','rb85 ','rb87 ','sr84 ','sr86 ','sr87 ', &
       'sr88 ','y89  ','zr90 ','zr91 ','zr92 ','zr94 ','zr96 ','nb93 ', &
       'mo92 ','mo94 ','mo95 ','mo96 ','mo97 ','mo98 ','mo100','ru96 '/

      data (namsol(j), j=121,240) / &
       'ru98 ','ru99 ','ru100','ru101','ru102','ru104','rh103','pd102', &
       'pd104','pd105','pd106','pd108','pd110','ag107','ag109','cd106', &
       'cd108','cd110','cd111','cd112','cd113','cd114','cd116','in113', &
       'in115','sn112','sn114','sn115','sn116','sn117','sn118','sn119', &
       'sn120','sn122','sn124','sb121','sb123','te120','te122','te123', &
       'te124','te125','te126','te128','te130','i127 ','xe124','xe126', &
       'xe128','xe129','xe130','xe131','xe132','xe134','xe136','cs133', &
       'ba130','ba132','ba134','ba135','ba136','ba137','ba138','la138', &
       'la139','ce136','ce138','ce140','ce142','pr141','nd142','nd143', &
       'nd144','nd145','nd146','nd148','nd150','sm144','sm147','sm148', &
       'sm149','sm150','sm152','sm154','eu151','eu153','gd152','gd154', &
       'gd155','gd156','gd157','gd158','gd160','tb159','dy156','dy158', &
       'dy160','dy161','dy162','dy163','dy164','ho165','er162','er164', &
       'er166','er167','er168','er170','tm169','yb168','yb170','yb171', &
       'yb172','yb173','yb174','yb176','lu175','lu176','hf174','hf176'/

      data (namsol(j), j=241,286) / &
       'hf177','hf178','hf179','hf180','ta180','ta181','w180 ','w182 ', &
       'w183 ','w184 ','w186 ','re185','re187','os184','os186','os187', &
       'os188','os189','os190','os192','ir191','ir193','pt190','pt192', &
       'pt194','pt195','pt196','pt198','au197','hg196','hg198','hg199', &
       'hg200','hg201','hg202','hg204','tl203','tl205','pb204','pb206', &
       'pb207','pb208','bi209','th232','u235 ','u238'/


! anders & grevesse 1989 solar mass fractions
        data (sol(i),i=1,45)/ &
           7.0573E-01, 4.8010E-05, 2.9291E-05, 2.7521E-01, 6.4957E-10, &
           9.3490E-09, 1.6619E-10, 1.0674E-09, 4.7301E-09, 3.0324E-03, &
           3.6501E-05, 1.1049E-03, 4.3634E-06, 9.5918E-03, 3.8873E-06, &
           2.1673E-05, 4.0515E-07, 1.6189E-03, 4.1274E-06, 1.3022E-04, &
           3.3394E-05, 5.1480E-04, 6.7664E-05, 7.7605E-05, 5.8052E-05, &
           6.5301E-04, 3.4257E-05, 2.3524E-05, 8.1551E-06, 3.9581E-04, &
           3.2221E-06, 1.8663E-05, 9.3793E-08, 2.5320E-06, 8.5449E-07, &
           7.7402E-05, 1.5379E-05, 2.6307E-08, 3.4725E-06, 4.4519E-10, &
           2.6342E-07, 5.9898E-05, 4.1964E-07, 8.9734E-07, 1.4135E-06/

        data (sol(i),i=46,90)/ &
             2.7926E-09, 1.3841E-07, 3.8929E-08, 2.2340E-07, 2.0805E-07, &
             2.1491E-06, 1.6361E-07, 1.6442E-07, 9.2579E-10, 3.7669E-07, &
             7.4240E-07, 1.4863E-05, 1.7160E-06, 4.3573E-07, 1.3286E-05, &
             7.1301E-05, 1.1686E-03, 2.8548E-05, 3.6971E-06, 3.3579E-06, &
             4.9441E-05, 1.9578E-05, 8.5944E-07, 2.7759E-06, 7.2687E-07, &
             5.7528E-07, 2.6471E-07, 9.9237E-07, 5.8765E-07, 8.7619E-08, &
             4.0593E-07, 1.3811E-08, 3.9619E-08, 2.7119E-08, 4.3204E-08, &
             5.9372E-08, 1.7136E-08, 8.1237E-08, 1.7840E-08, 1.2445E-08, &
             1.0295E-09, 1.0766E-08, 9.1542E-09, 2.9003E-08, 6.2529E-08/

        data (sol(i),i=91,135)/ &
             1.1823E-08, 1.1950E-08, 1.2006E-08, 3.0187E-10, 2.0216E-09, &
             1.0682E-08, 1.0833E-08, 5.4607E-08, 1.7055E-08, 1.1008E-08, &
             4.3353E-09, 2.8047E-10, 5.0468E-09, 3.6091E-09, 4.3183E-08, &
             1.0446E-08, 1.3363E-08, 2.9463E-09, 4.5612E-09, 4.7079E-09, &
             7.7706E-10, 1.6420E-09, 8.7966E-10, 5.6114E-10, 9.7562E-10, &
             1.0320E-09, 5.9868E-10, 1.5245E-09, 6.2225E-10, 2.5012E-10, &
             8.6761E-11, 5.9099E-10, 5.9190E-10, 8.0731E-10, 1.5171E-09, &
             9.1547E-10, 8.9625E-10, 3.6637E-11, 4.0775E-10, 8.2335E-10, &
             1.0189E-09, 1.0053E-09, 4.5354E-10, 6.8205E-10, 6.4517E-10/

        data (sol(i),i=136,180)/ &
             5.3893E-11, 3.9065E-11, 5.5927E-10, 5.7839E-10, 1.0992E-09, &
             5.6309E-10, 1.3351E-09, 3.5504E-10, 2.2581E-11, 5.1197E-10, &
             1.0539E-10, 7.1802E-11, 3.9852E-11, 1.6285E-09, 8.6713E-10, &
             2.7609E-09, 9.8731E-10, 3.7639E-09, 5.4622E-10, 6.9318E-10, &
             5.4174E-10, 4.1069E-10, 1.3052E-11, 3.8266E-10, 1.3316E-10, &
             7.1827E-10, 1.0814E-09, 3.1553E-09, 4.9538E-09, 5.3600E-09, &
             2.8912E-09, 1.7910E-11, 1.6223E-11, 3.3349E-10, 4.1767E-09, &
             6.7411E-10, 3.3799E-09, 4.1403E-09, 1.5558E-09, 1.2832E-09, &
             1.2515E-09, 1.5652E-11, 1.5125E-11, 3.6946E-10, 1.0108E-09/

        data (sol(i),i=181,225)/ &
             1.2144E-09, 1.7466E-09, 1.1240E-08, 1.3858E-12, 1.5681E-09, &
             7.4306E-12, 9.9136E-12, 3.5767E-09, 4.5258E-10, 5.9562E-10, &
             8.0817E-10, 3.6533E-10, 7.1757E-10, 2.5198E-10, 5.2441E-10, &
             1.7857E-10, 1.7719E-10, 2.9140E-11, 1.4390E-10, 1.0931E-10, &
             1.3417E-10, 7.2470E-11, 2.6491E-10, 2.2827E-10, 1.7761E-10, &
             1.9660E-10, 2.5376E-12, 2.8008E-11, 1.9133E-10, 2.6675E-10, &
             2.0492E-10, 3.2772E-10, 2.9180E-10, 2.8274E-10, 8.6812E-13, &
             1.4787E-12, 3.7315E-11, 3.0340E-10, 4.1387E-10, 4.0489E-10, &
             4.6047E-10, 3.7104E-10, 1.4342E-12, 1.6759E-11, 3.5397E-10/

        data (sol(i),i=226,270)/ &
             2.4332E-10, 2.8557E-10, 1.6082E-10, 1.6159E-10, 1.3599E-12, &
             3.2509E-11, 1.5312E-10, 2.3624E-10, 1.7504E-10, 3.4682E-10, &
             1.4023E-10, 1.5803E-10, 4.2293E-12, 1.0783E-12, 3.4992E-11, &
             1.2581E-10, 1.8550E-10, 9.3272E-11, 2.4131E-10, 1.1292E-14, &
             9.4772E-11, 7.8768E-13, 1.6113E-10, 8.7950E-11, 1.8989E-10, &
             1.7878E-10, 9.0315E-11, 1.5326E-10, 5.6782E-13, 5.0342E-11, &
             5.1086E-11, 4.2704E-10, 5.2110E-10, 8.5547E-10, 1.3453E-09, &
             1.1933E-09, 2.0211E-09, 8.1702E-13, 5.0994E-11, 2.1641E-09, &
             2.2344E-09, 1.6757E-09, 4.8231E-10, 9.3184E-10, 2.3797E-12/

        data (sol(i),i=271,286)/ &
             1.7079E-10, 2.8843E-10, 3.9764E-10, 2.2828E-10, 5.1607E-10, &
             1.2023E-10, 2.7882E-10, 6.7411E-10, 3.1529E-10, 3.1369E-09, &
             3.4034E-09, 9.6809E-09, 7.6127E-10, 1.9659E-10, 3.8519E-13, &
             5.3760E-11/


! charge of the stable isotopes

        data (izsol(i),i=1,117)/ &
         1,   1,   2,   2,   3,   3,   4,   5,   5,   6,   6,   7,   7, &
         8,   8,   8,   9,  10,  10,  10,  11,  12,  12,  12,  13,  14, &
        14,  14,  15,  16,  16,  16,  16,  17,  17,  18,  18,  18,  19, &
        19,  19,  20,  20,  20,  20,  20,  20,  21,  22,  22,  22,  22, &
        22,  23,  23,  24,  24,  24,  24,  25,  26,  26,  26,  26,  27, &
        28,  28,  28,  28,  28,  29,  29,  30,  30,  30,  30,  30,  31, &
        31,  32,  32,  32,  32,  32,  33,  34,  34,  34,  34,  34,  34, &
        35,  35,  36,  36,  36,  36,  36,  36,  37,  37,  38,  38,  38, &
        38,  39,  40,  40,  40,  40,  40,  41,  42,  42,  42,  42,  42/

        data (izsol(i),i=118,234)/ &
        42,  42,  44,  44,  44,  44,  44,  44,  44,  45,  46,  46,  46, &
        46,  46,  46,  47,  47,  48,  48,  48,  48,  48,  48,  48,  48, &
        49,  49,  50,  50,  50,  50,  50,  50,  50,  50,  50,  50,  51, &
        51,  52,  52,  52,  52,  52,  52,  52,  52,  53,  54,  54,  54, &
        54,  54,  54,  54,  54,  54,  55,  56,  56,  56,  56,  56,  56, &
        56,  57,  57,  58,  58,  58,  58,  59,  60,  60,  60,  60,  60, &
        60,  60,  62,  62,  62,  62,  62,  62,  62,  63,  63,  64,  64, &
        64,  64,  64,  64,  64,  65,  66,  66,  66,  66,  66,  66,  66, &
        67,  68,  68,  68,  68,  68,  68,  69,  70,  70,  70,  70,  70/

        data (izsol(i),i=235,286)/ &
        70,  70,  71,  71,  72,  72,  72,  72,  72,  72,  73,  73,  74, &
        74,  74,  74,  74,  75,  75,  76,  76,  76,  76,  76,  76,  76, &
        77,  77,  78,  78,  78,  78,  78,  78,  79,  80,  80,  80,  80, &
        80,  80,  80,  81,  81,  82,  82,  82,  82,  83,  90,  92,  92/


! number of nucleons (protons and neutrons) in the stable isotopes

        data (iasol(i),i=1,117)/ &
         1,   2,   3,   4,   6,   7,   9,  10,  11,  12,  13,  14,  15, &
        16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28, &
        29,  30,  31,  32,  33,  34,  36,  35,  37,  36,  38,  40,  39, &
        40,  41,  40,  42,  43,  44,  46,  48,  45,  46,  47,  48,  49, &
        50,  50,  51,  50,  52,  53,  54,  55,  54,  56,  57,  58,  59, &
        58,  60,  61,  62,  64,  63,  65,  64,  66,  67,  68,  70,  69, &
        71,  70,  72,  73,  74,  76,  75,  74,  76,  77,  78,  80,  82, &
        79,  81,  78,  80,  82,  83,  84,  86,  85,  87,  84,  86,  87, &
        88,  89,  90,  91,  92,  94,  96,  93,  92,  94,  95,  96,  97/

        data (iasol(i),i=118,234)/ &
        98, 100,  96,  98,  99, 100, 101, 102, 104, 103, 102, 104, 105, &
       106, 108, 110, 107, 109, 106, 108, 110, 111, 112, 113, 114, 116, &
       113, 115, 112, 114, 115, 116, 117, 118, 119, 120, 122, 124, 121, &
       123, 120, 122, 123, 124, 125, 126, 128, 130, 127, 124, 126, 128, &
       129, 130, 131, 132, 134, 136, 133, 130, 132, 134, 135, 136, 137, &
       138, 138, 139, 136, 138, 140, 142, 141, 142, 143, 144, 145, 146, &
       148, 150, 144, 147, 148, 149, 150, 152, 154, 151, 153, 152, 154, &
       155, 156, 157, 158, 160, 159, 156, 158, 160, 161, 162, 163, 164, &
       165, 162, 164, 166, 167, 168, 170, 169, 168, 170, 171, 172, 173/

        data (iasol(i),i=235,286)/ &
       174, 176, 175, 176, 174, 176, 177, 178, 179, 180, 180, 181, 180, &
       182, 183, 184, 186, 185, 187, 184, 186, 187, 188, 189, 190, 192, &
       191, 193, 190, 192, 194, 195, 196, 198, 197, 196, 198, 199, 200, &
       201, 202, 204, 203, 205, 204, 206, 207, 208, 209, 232, 235, 238/



! jcode tells the type progenitors each stable species can have.
! jcode = 0 if the species is the only stable one of that a
!       = 1 if the species can have proton-rich progenitors
!       = 2 if the species can have neutron-rich progenitors
!       = 3 if the species can only be made as itself (eg k40)

        data (jcode(i),i=1,117)/ &
         0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, &
         0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, &
         0,   0,   0,   0,   0,   0,   2,   0,   0,   1,   0,   2,   0, &
         3,   0,   1,   0,   0,   0,   2,   2,   0,   1,   0,   1,   0, &
         2,   3,   0,   1,   0,   0,   2,   0,   1,   0,   0,   2,   0, &
         1,   0,   0,   0,   2,   0,   0,   1,   0,   0,   0,   2,   0, &
         0,   1,   0,   0,   2,   2,   0,   1,   1,   0,   2,   2,   2, &
         0,   0,   1,   1,   1,   0,   2,   2,   0,   2,   1,   1,   1, &
         0,   0,   0,   0,   2,   2,   2,   0,   1,   1,   0,   3,   0/

        data (jcode(i),i=118,234)/ &
         2,   2,   1,   1,   0,   1,   0,   2,   2,   0,   1,   1,   0, &
         2,   2,   2,   0,   0,   1,   1,   1,   0,   2,   2,   2,   2, &
         1,   2,   1,   1,   1,   1,   0,   0,   0,   2,   2,   2,   0, &
         2,   1,   1,   1,   3,   0,   2,   2,   2,   0,   1,   1,   1, &
         0,   3,   0,   2,   2,   2,   0,   1,   1,   1,   0,   3,   0, &
         2,   3,   0,   1,   1,   0,   2,   0,   1,   0,   2,   0,   0, &
         2,   2,   1,   0,   1,   0,   1,   2,   2,   0,   0,   1,   1, &
         0,   2,   0,   2,   2,   0,   1,   1,   1,   0,   2,   0,   2, &
         0,   1,   1,   0,   0,   2,   2,   0,   1,   1,   0,   0,   0/

        data (jcode(i),i=235,286)/ &
         2,   2,   0,   3,   1,   1,   0,   0,   0,   2,   3,   0,   1, &
         0,   0,   2,   2,   0,   2,   1,   1,   1,   0,   0,   2,   2, &
         0,   0,   1,   1,   0,   0,   2,   2,   0,   1,   1,   0,   0, &
         0,   0,   2,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0/



!      include 'solar_data_lodders_2003.f'



! sum; stuff residual into hydrogen
      if (ifirst .eq. 0) then
       ifirst = 1
       sum = 0.0d0
       do j=1,solsiz
        sum = sum + sol(j)
       enddo
       sum = 1.0d0 - sum
       sol(1) = sol(1) + sum

       sum  = 0.0d0
       do j=1,solsiz
        if (izsol(j) .ge. 3) then
         sum  = sum + sol(j)
        endif
       enddo
       zsol = sum

       sum = 0.0d0
       do j=1,solsiz
        if (izsol(j) .ge. 3) then
         sum = sum + float(izsol(j))/float(iasol(j))*sol(j)
        endif
       enddo
       yesol = sum
      end if

! straight sweep
      andgrev  = 0.0d0
      z        = 0.0d0
      a        = 0.0d0

      if (len(nam) .lt. 5) stop 'nam < 5 characters in routine andgrev'
      do i=1,solsiz
       if ( namsol(i)(1:5) .eq. nam(1:5) ) then

! load the mass fraction, charge, and number of nucleons
        andgrev = sol(i)
        z       = float(izsol(i))
        a       = float(iasol(i))

! load the elmental mass fraction associated with this isotope
        xelem = 0.0d0
        jbeg = max(1,i-12)
        jend = min(i+12,solsiz)
        do j=jbeg,jend
         if (izsol(j) .eq. z) xelem = xelem + sol(j)
        enddo

! bail
        return
       end if
      enddo

!      write(6,*) 'warning: no such entry ',nam(1:5)

      return
      end







      subroutine decay_andgrev(nin,nz,na,xin, &
                               nout,namout,nzout,naout,xout,ag)
      include 'implno.dek'

! converts a radioactive nucleosynthesis down to their final
! stable mass fractions

! declare the pass
      integer          nin,nout
      character*5      namout(nout)
      integer          nz(nin),na(nin),nzout(nout),naout(nout)
      double precision xin(nin),xout(nout),ag(nout)


! for the solar abundance data
      integer          solsiz
      parameter        (solsiz = 286)
      character*5      namsol(solsiz)
      integer          izsol(solsiz),iasol(solsiz),jcode(solsiz)
      double precision sol(solsiz)


! local variables
      integer          i,j
      double precision termx,sum,xx


! bring in the solar abundance data

! names of the stable isotopes
      data (namsol(j), j=1,120) / &
       'h1   ','h2   ','he3  ','he4  ','li6  ','li7  ','be9  ','b10  ', &
       'b11  ','c12  ','c13  ','n14  ','n15  ','o16  ','o17  ','o18  ', &
       'f19  ','ne20 ','ne21 ','ne22 ','na23 ','mg24 ','mg25 ','mg26 ', &
       'al27 ','si28 ','si29 ','si30 ','p31  ','s32  ','s33  ','s34  ', &
       's36  ','cl35 ','cl37 ','ar36 ','ar38 ','ar40 ','k39  ','k40  ', &
       'k41  ','ca40 ','ca42 ','ca43 ','ca44 ','ca46 ','ca48 ','sc45 ', &
       'ti46 ','ti47 ','ti48 ','ti49 ','ti50 ','v50  ','v51  ','cr50 ', &
       'cr52 ','cr53 ','cr54 ','mn55 ','fe54 ','fe56 ','fe57 ','fe58 ', &
       'co59 ','ni58 ','ni60 ','ni61 ','ni62 ','ni64 ','cu63 ','cu65 ', &
       'zn64 ','zn66 ','zn67 ','zn68 ','zn70 ','ga69 ','ga71 ','ge70 ', &
       'ge72 ','ge73 ','ge74 ','ge76 ','as75 ','se74 ','se76 ','se77 ', &
       'se78 ','se80 ','se82 ','br79 ','br81 ','kr78 ','kr80 ','kr82 ', &
       'kr83 ','kr84 ','kr86 ','rb85 ','rb87 ','sr84 ','sr86 ','sr87 ', &
       'sr88 ','y89  ','zr90 ','zr91 ','zr92 ','zr94 ','zr96 ','nb93 ', &
       'mo92 ','mo94 ','mo95 ','mo96 ','mo97 ','mo98 ','mo100','ru96 '/

      data (namsol(j), j=121,240) / &
       'ru98 ','ru99 ','ru100','ru101','ru102','ru104','rh103','pd102', &
       'pd104','pd105','pd106','pd108','pd110','ag107','ag109','cd106', &
       'cd108','cd110','cd111','cd112','cd113','cd114','cd116','in113', &
       'in115','sn112','sn114','sn115','sn116','sn117','sn118','sn119', &
       'sn120','sn122','sn124','sb121','sb123','te120','te122','te123', &
       'te124','te125','te126','te128','te130','i127 ','xe124','xe126', &
       'xe128','xe129','xe130','xe131','xe132','xe134','xe136','cs133', &
       'ba130','ba132','ba134','ba135','ba136','ba137','ba138','la138', &
       'la139','ce136','ce138','ce140','ce142','pr141','nd142','nd143', &
       'nd144','nd145','nd146','nd148','nd150','sm144','sm147','sm148', &
       'sm149','sm150','sm152','sm154','eu151','eu153','gd152','gd154', &
       'gd155','gd156','gd157','gd158','gd160','tb159','dy156','dy158', &
       'dy160','dy161','dy162','dy163','dy164','ho165','er162','er164', &
       'er166','er167','er168','er170','tm169','yb168','yb170','yb171', &
       'yb172','yb173','yb174','yb176','lu175','lu176','hf174','hf176'/

      data (namsol(j), j=241,286) / &
       'hf177','hf178','hf179','hf180','ta180','ta181','w180 ','w182 ', &
       'w183 ','w184 ','w186 ','re185','re187','os184','os186','os187', &
       'os188','os189','os190','os192','ir191','ir193','pt190','pt192', &
       'pt194','pt195','pt196','pt198','au197','hg196','hg198','hg199', &
       'hg200','hg201','hg202','hg204','tl203','tl205','pb204','pb206', &
       'pb207','pb208','bi209','th232','u235 ','u238'/


! anders & grevesse 1989 solar mass fractions
        data (sol(i),i=1,45)/ &
           7.0573E-01, 4.8010E-05, 2.9291E-05, 2.7521E-01, 6.4957E-10, &
           9.3490E-09, 1.6619E-10, 1.0674E-09, 4.7301E-09, 3.0324E-03, &
           3.6501E-05, 1.1049E-03, 4.3634E-06, 9.5918E-03, 3.8873E-06, &
           2.1673E-05, 4.0515E-07, 1.6189E-03, 4.1274E-06, 1.3022E-04, &
           3.3394E-05, 5.1480E-04, 6.7664E-05, 7.7605E-05, 5.8052E-05, &
           6.5301E-04, 3.4257E-05, 2.3524E-05, 8.1551E-06, 3.9581E-04, &
           3.2221E-06, 1.8663E-05, 9.3793E-08, 2.5320E-06, 8.5449E-07, &
           7.7402E-05, 1.5379E-05, 2.6307E-08, 3.4725E-06, 4.4519E-10, &
           2.6342E-07, 5.9898E-05, 4.1964E-07, 8.9734E-07, 1.4135E-06/

        data (sol(i),i=46,90)/ &
             2.7926E-09, 1.3841E-07, 3.8929E-08, 2.2340E-07, 2.0805E-07, &
             2.1491E-06, 1.6361E-07, 1.6442E-07, 9.2579E-10, 3.7669E-07, &
             7.4240E-07, 1.4863E-05, 1.7160E-06, 4.3573E-07, 1.3286E-05, &
             7.1301E-05, 1.1686E-03, 2.8548E-05, 3.6971E-06, 3.3579E-06, &
             4.9441E-05, 1.9578E-05, 8.5944E-07, 2.7759E-06, 7.2687E-07, &
             5.7528E-07, 2.6471E-07, 9.9237E-07, 5.8765E-07, 8.7619E-08, &
             4.0593E-07, 1.3811E-08, 3.9619E-08, 2.7119E-08, 4.3204E-08, &
             5.9372E-08, 1.7136E-08, 8.1237E-08, 1.7840E-08, 1.2445E-08, &
             1.0295E-09, 1.0766E-08, 9.1542E-09, 2.9003E-08, 6.2529E-08/

        data (sol(i),i=91,135)/ &
             1.1823E-08, 1.1950E-08, 1.2006E-08, 3.0187E-10, 2.0216E-09, &
             1.0682E-08, 1.0833E-08, 5.4607E-08, 1.7055E-08, 1.1008E-08, &
             4.3353E-09, 2.8047E-10, 5.0468E-09, 3.6091E-09, 4.3183E-08, &
             1.0446E-08, 1.3363E-08, 2.9463E-09, 4.5612E-09, 4.7079E-09, &
             7.7706E-10, 1.6420E-09, 8.7966E-10, 5.6114E-10, 9.7562E-10, &
             1.0320E-09, 5.9868E-10, 1.5245E-09, 6.2225E-10, 2.5012E-10, &
             8.6761E-11, 5.9099E-10, 5.9190E-10, 8.0731E-10, 1.5171E-09, &
             9.1547E-10, 8.9625E-10, 3.6637E-11, 4.0775E-10, 8.2335E-10, &
             1.0189E-09, 1.0053E-09, 4.5354E-10, 6.8205E-10, 6.4517E-10/

        data (sol(i),i=136,180)/ &
             5.3893E-11, 3.9065E-11, 5.5927E-10, 5.7839E-10, 1.0992E-09, &
             5.6309E-10, 1.3351E-09, 3.5504E-10, 2.2581E-11, 5.1197E-10, &
             1.0539E-10, 7.1802E-11, 3.9852E-11, 1.6285E-09, 8.6713E-10, &
             2.7609E-09, 9.8731E-10, 3.7639E-09, 5.4622E-10, 6.9318E-10, &
             5.4174E-10, 4.1069E-10, 1.3052E-11, 3.8266E-10, 1.3316E-10, &
             7.1827E-10, 1.0814E-09, 3.1553E-09, 4.9538E-09, 5.3600E-09, &
             2.8912E-09, 1.7910E-11, 1.6223E-11, 3.3349E-10, 4.1767E-09, &
             6.7411E-10, 3.3799E-09, 4.1403E-09, 1.5558E-09, 1.2832E-09, &
             1.2515E-09, 1.5652E-11, 1.5125E-11, 3.6946E-10, 1.0108E-09/

        data (sol(i),i=181,225)/ &
             1.2144E-09, 1.7466E-09, 1.1240E-08, 1.3858E-12, 1.5681E-09, &
             7.4306E-12, 9.9136E-12, 3.5767E-09, 4.5258E-10, 5.9562E-10, &
             8.0817E-10, 3.6533E-10, 7.1757E-10, 2.5198E-10, 5.2441E-10, &
             1.7857E-10, 1.7719E-10, 2.9140E-11, 1.4390E-10, 1.0931E-10, &
             1.3417E-10, 7.2470E-11, 2.6491E-10, 2.2827E-10, 1.7761E-10, &
             1.9660E-10, 2.5376E-12, 2.8008E-11, 1.9133E-10, 2.6675E-10, &
             2.0492E-10, 3.2772E-10, 2.9180E-10, 2.8274E-10, 8.6812E-13, &
             1.4787E-12, 3.7315E-11, 3.0340E-10, 4.1387E-10, 4.0489E-10, &
             4.6047E-10, 3.7104E-10, 1.4342E-12, 1.6759E-11, 3.5397E-10/

        data (sol(i),i=226,270)/ &
             2.4332E-10, 2.8557E-10, 1.6082E-10, 1.6159E-10, 1.3599E-12, &
             3.2509E-11, 1.5312E-10, 2.3624E-10, 1.7504E-10, 3.4682E-10, &
             1.4023E-10, 1.5803E-10, 4.2293E-12, 1.0783E-12, 3.4992E-11, &
             1.2581E-10, 1.8550E-10, 9.3272E-11, 2.4131E-10, 1.1292E-14, &
             9.4772E-11, 7.8768E-13, 1.6113E-10, 8.7950E-11, 1.8989E-10, &
             1.7878E-10, 9.0315E-11, 1.5326E-10, 5.6782E-13, 5.0342E-11, &
             5.1086E-11, 4.2704E-10, 5.2110E-10, 8.5547E-10, 1.3453E-09, &
             1.1933E-09, 2.0211E-09, 8.1702E-13, 5.0994E-11, 2.1641E-09, &
             2.2344E-09, 1.6757E-09, 4.8231E-10, 9.3184E-10, 2.3797E-12/

        data (sol(i),i=271,286)/ &
             1.7079E-10, 2.8843E-10, 3.9764E-10, 2.2828E-10, 5.1607E-10, &
             1.2023E-10, 2.7882E-10, 6.7411E-10, 3.1529E-10, 3.1369E-09, &
             3.4034E-09, 9.6809E-09, 7.6127E-10, 1.9659E-10, 3.8519E-13, &
             5.3760E-11/


! charge of the stable isotopes

        data (izsol(i),i=1,117)/ &
         1,   1,   2,   2,   3,   3,   4,   5,   5,   6,   6,   7,   7, &
         8,   8,   8,   9,  10,  10,  10,  11,  12,  12,  12,  13,  14, &
        14,  14,  15,  16,  16,  16,  16,  17,  17,  18,  18,  18,  19, &
        19,  19,  20,  20,  20,  20,  20,  20,  21,  22,  22,  22,  22, &
        22,  23,  23,  24,  24,  24,  24,  25,  26,  26,  26,  26,  27, &
        28,  28,  28,  28,  28,  29,  29,  30,  30,  30,  30,  30,  31, &
        31,  32,  32,  32,  32,  32,  33,  34,  34,  34,  34,  34,  34, &
        35,  35,  36,  36,  36,  36,  36,  36,  37,  37,  38,  38,  38, &
        38,  39,  40,  40,  40,  40,  40,  41,  42,  42,  42,  42,  42/

        data (izsol(i),i=118,234)/ &
        42,  42,  44,  44,  44,  44,  44,  44,  44,  45,  46,  46,  46, &
        46,  46,  46,  47,  47,  48,  48,  48,  48,  48,  48,  48,  48, &
        49,  49,  50,  50,  50,  50,  50,  50,  50,  50,  50,  50,  51, &
        51,  52,  52,  52,  52,  52,  52,  52,  52,  53,  54,  54,  54, &
        54,  54,  54,  54,  54,  54,  55,  56,  56,  56,  56,  56,  56, &
        56,  57,  57,  58,  58,  58,  58,  59,  60,  60,  60,  60,  60, &
        60,  60,  62,  62,  62,  62,  62,  62,  62,  63,  63,  64,  64, &
        64,  64,  64,  64,  64,  65,  66,  66,  66,  66,  66,  66,  66, &
        67,  68,  68,  68,  68,  68,  68,  69,  70,  70,  70,  70,  70/

        data (izsol(i),i=235,286)/ &
        70,  70,  71,  71,  72,  72,  72,  72,  72,  72,  73,  73,  74, &
        74,  74,  74,  74,  75,  75,  76,  76,  76,  76,  76,  76,  76, &
        77,  77,  78,  78,  78,  78,  78,  78,  79,  80,  80,  80,  80, &
        80,  80,  80,  81,  81,  82,  82,  82,  82,  83,  90,  92,  92/


! number of nucleons (protons and neutrons) in the stable isotopes

        data (iasol(i),i=1,117)/ &
         1,   2,   3,   4,   6,   7,   9,  10,  11,  12,  13,  14,  15, &
        16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28, &
        29,  30,  31,  32,  33,  34,  36,  35,  37,  36,  38,  40,  39, &
        40,  41,  40,  42,  43,  44,  46,  48,  45,  46,  47,  48,  49, &
        50,  50,  51,  50,  52,  53,  54,  55,  54,  56,  57,  58,  59, &
        58,  60,  61,  62,  64,  63,  65,  64,  66,  67,  68,  70,  69, &
        71,  70,  72,  73,  74,  76,  75,  74,  76,  77,  78,  80,  82, &
        79,  81,  78,  80,  82,  83,  84,  86,  85,  87,  84,  86,  87, &
        88,  89,  90,  91,  92,  94,  96,  93,  92,  94,  95,  96,  97/

        data (iasol(i),i=118,234)/ &
        98, 100,  96,  98,  99, 100, 101, 102, 104, 103, 102, 104, 105, &
       106, 108, 110, 107, 109, 106, 108, 110, 111, 112, 113, 114, 116, &
       113, 115, 112, 114, 115, 116, 117, 118, 119, 120, 122, 124, 121, &
       123, 120, 122, 123, 124, 125, 126, 128, 130, 127, 124, 126, 128, &
       129, 130, 131, 132, 134, 136, 133, 130, 132, 134, 135, 136, 137, &
       138, 138, 139, 136, 138, 140, 142, 141, 142, 143, 144, 145, 146, &
       148, 150, 144, 147, 148, 149, 150, 152, 154, 151, 153, 152, 154, &
       155, 156, 157, 158, 160, 159, 156, 158, 160, 161, 162, 163, 164, &
       165, 162, 164, 166, 167, 168, 170, 169, 168, 170, 171, 172, 173/

        data (iasol(i),i=235,286)/ &
       174, 176, 175, 176, 174, 176, 177, 178, 179, 180, 180, 181, 180, &
       182, 183, 184, 186, 185, 187, 184, 186, 187, 188, 189, 190, 192, &
       191, 193, 190, 192, 194, 195, 196, 198, 197, 196, 198, 199, 200, &
       201, 202, 204, 203, 205, 204, 206, 207, 208, 209, 232, 235, 238/



! jcode tells the type progenitors each stable species can have.
! jcode = 0 if the species is the only stable one of that a
!       = 1 if the species can have proton-rich progenitors
!       = 2 if the species can have neutron-rich progenitors
!       = 3 if the species can only be made as itself (eg k40)

        data (jcode(i),i=1,117)/ &
         0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, &
         0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, &
         0,   0,   0,   0,   0,   0,   2,   0,   0,   1,   0,   2,   0, &
         3,   0,   1,   0,   0,   0,   2,   2,   0,   1,   0,   1,   0, &
         2,   3,   0,   1,   0,   0,   2,   0,   1,   0,   0,   2,   0, &
         1,   0,   0,   0,   2,   0,   0,   1,   0,   0,   0,   2,   0, &
         0,   1,   0,   0,   2,   2,   0,   1,   1,   0,   2,   2,   2, &
         0,   0,   1,   1,   1,   0,   2,   2,   0,   2,   1,   1,   1, &
         0,   0,   0,   0,   2,   2,   2,   0,   1,   1,   0,   3,   0/

        data (jcode(i),i=118,234)/ &
         2,   2,   1,   1,   0,   1,   0,   2,   2,   0,   1,   1,   0, &
         2,   2,   2,   0,   0,   1,   1,   1,   0,   2,   2,   2,   2, &
         1,   2,   1,   1,   1,   1,   0,   0,   0,   2,   2,   2,   0, &
         2,   1,   1,   1,   3,   0,   2,   2,   2,   0,   1,   1,   1, &
         0,   3,   0,   2,   2,   2,   0,   1,   1,   1,   0,   3,   0, &
         2,   3,   0,   1,   1,   0,   2,   0,   1,   0,   2,   0,   0, &
         2,   2,   1,   0,   1,   0,   1,   2,   2,   0,   0,   1,   1, &
         0,   2,   0,   2,   2,   0,   1,   1,   1,   0,   2,   0,   2, &
         0,   1,   1,   0,   0,   2,   2,   0,   1,   1,   0,   0,   0/

        data (jcode(i),i=235,286)/ &
         2,   2,   0,   3,   1,   1,   0,   0,   0,   2,   3,   0,   1, &
         0,   0,   2,   2,   0,   2,   1,   1,   1,   0,   0,   2,   2, &
         0,   0,   1,   1,   0,   0,   2,   2,   0,   1,   1,   0,   0, &
         0,   0,   2,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0/



!      include 'solar_data_lodders_2003.f'


! initialize
      if (nout .lt. solsiz) stop 'not < solsiz in routine decay_andgrev'
      do i=1,solsiz
       xout(i)   = 0.0d0
       namout(i) = namsol(i)
       nzout(i)  = izsol(i)
       naout(i)  = iasol(i)
      enddo


! start the conversion
      do 400 i=1,nin

! for every isotope in the solar list
       do 390  j=1,solsiz
        if (na(i) .ne. iasol(j)) goto 390
        if (jcode(j) .eq. 0) goto 350
        if (nz(i).ge.izsol(j) .and. jcode(j).eq.1) goto 350
        if (nz(i).le.izsol(j) .and. jcode(j).eq.2) goto 350
        if (nz(i).eq.izsol(j) .and. jcode(j).eq.3) goto 350
        goto 390
350     termx   = xin(i)
        xout(j) = xout(j) + termx

! record the isotope that makes the largest contribution
! to this stable isotope
!        if (termx .le. prod(j)) goto 389
!        prod(j)  = termx
!        nzprod(j) = nz(i)
!        naprod(j) = na(i)
!389     continue


        goto 400
390    continue

400   continue


! scaled to solar
      do i=1,solsiz
       ag(i)   = xout(i)/sol(i)
      enddo

      return
      end

!---------------------------------------------------------------------







!---------------------------------------------------------------------
      subroutine net_output(kount,x,y,derivs)
      include 'implno.dek'
      include 'const.dek'
      include 'vector_eos.dek'
      include 'burn_common.dek'
      include 'network.dek'

! writes the output

! declare the pass
      external         derivs
      integer          kount
      double precision x,y(*)


! local variables
      character*8      atim
      character*9      adat
      character*80     string
      integer          ii,k,kk,j,lop,ilop,jrem,kb,ke,nn,lenstr,ncolumns
      parameter        (ncolumns = 12)
      double precision xcons,ycons,yex,xdum(abignet), &
                       abar,zbar,wbar,ye,xcess,zero,tdum,ddum,pdum, &
                       ener,denerdt,xx,posx,velx
      parameter        (zero = 0.0d0)


! for nse
      integer          igues
      double precision xmun,xmup,t9,tau_nse,tau_qse,taud


! popular format statements
01    format(1x,'*',t13,a,t33,a,t47,a,t61,a,t75,a,t89,a, &
                    t103,a,t117,a,t131,a,t145,a,t159,a, & 
                    t173,a,t187,a,t201,a,t215,a)
03    format(a50,i4.4,a2,i8,a)
04    format(1x,i6,1pe20.12,1p15e14.6)
05    format(1x,i6,1pe20.12,1p12e14.6)
07    format(1x,'* ',a,5(a,1pe11.3))


!      write(6,*) kount,neqs,nzone

! initialize the files with their headers
      if (kount .eq. 1) then



! for every spatial zone
       do k=1,max(1,nzone)
        kk = neqs*(k-1)

! logical unit 22 records the energetics
!  AT comments out FROM.
!        write(string,03) hfile,0,'_z',k,'.dat'
!        call sqeeze(string)
!        call today(adat,atim)
!        open (unit=22, file=string, status='unknown')


! logical unit 23 records the thermodynamics
!        write(string,03) hfile,1,'_z',k,'.dat'
!        call sqeeze(string)
!        open (unit=23, file=string, status='unknown')


!         write(22,01) adat,atim
!         write(23,01) adat,atim

!        if (one_step) then
!         write(22,07) 'one_step:','  btemp=',btemp,' bden=',bden
!         write(23,07) 'one_step:','  btemp=',btemp,' bden=',bden

!        else if (hydrostatic) then
!         write(22,07) 'hydrostatic:','  btemp=',btemp,' bden=',bden
!         write(23,07) 'hydrostatic:','  btemp=',btemp,' bden=',bden

!        else if (expansion) then
!         write(22,07) 'expansion:','  temp0=',temp0,' den0=',den0, &
!                     ' temp_stop=',temp_stop * temp_scale
!         write(23,07) 'expansion:','  temp0=',temp0,' den0=',den0, &
!                     ' temp_stop=',temp_stop * temp_scale

!        else if (self_heat_const_den) then
!         write(22,07) 'self_heat:','   temp0=',y(itemp+kk), &
!                                  '   den0=',y(iden+kk)
!         write(23,07) 'self_heat:','   temp0=',y(itemp+kk), &
!                                  '   den0=',y(iden+kk)

!        else if (self_heat_const_pres) then
!         write(22,07) 'self_heat:','   temp0=',y(itemp+kk), &
!                                  '   den0=',y(iden+kk)
!         write(23,07) 'self_heat:','   temp0=',y(itemp+kk), &
!                                  '   den0=',y(iden+kk)

!        else if (bbang) then
!         write(22,07) 'big bang:','   eta=',eta1, &
!                                  '   Nnu=',xnnu, &
!                                  '   H0=',hubble, &
!                                  '   Tcmb=',cmbtemp
!         write(23,07) 'big_bang:','   eta=',eta1, &
!                                  '   Nnu=',xnnu, &
!                                  '   H0=',hubble, &
!                                  '   Tcmb=',cmbtemp


!        else if (trho_hist) then
!         call update2(zero,tdum,ddum)
!         write(22,07) 'trho_hist:','  mass interior =',mint, &
!                                   '  shell mass =',mshell
!         write(23,07) 'trho_hist:',' mass interior =',mint, &
!                                   '  shell mass =',mshell

!        end if


!        write(22,01) 'time','temp','den','ener','sdot','sneut', &
!                     's-snu','ye','1-sum'

!        write(23,01) 'time','pos','vel','temp','den','pres','ener', &
!                     'entr','cs'


! close up the files
!        close(unit=22)
!        close(unit=23)
!  AT comments out TO.


! end of spatial loop
       enddo



! done writing thermodynamic headers



! for every spatial zone
       do k=1,max(1,nzone)
        kk = neqs*(k-1)


! write out the isotopic mass fractions in blocks of ncolumns
! lop is how many groups of ncolumns exist; jrem is the remainder
        lop  = ionmax/ncolumns
        jrem  = ionmax - ncolumns*lop
        do ilop = 1,lop+1
         kb = 1 + ncolumns*(ilop-1)
         ke = ncolumns + ncolumns*(ilop-1)
         if (ilop .eq. lop+1  .and. jrem .eq. 0) goto 50
         if (ilop .eq. lop+1) ke = ionmax


! logical unit 34 records the abundance evolution
! open the output file
!  AT comments out FROM.
!         write(string,03) hfile,ilop+1,'_z',k,'.dat'
!         call sqeeze(string)
!         open (unit=34, file=string, status='unknown')

!         write(34,01) adat,atim


!        if (one_step) then
!         write(34,07) 'one_step:','  btemp=',btemp,' bden=',bden

!        else if (hydrostatic) then
!         write(34,07) 'hydrostatic:','  btemp=',btemp,' bden=',bden

!        else if (expansion) then
!         write(34,07) 'expansion:','  temp0=',temp0,' den0=',den0, &
!                                  ' temp_stop=',temp_stop

!        else if (self_heat_const_den) then
!         write(34,07) 'self_heat:','   temp0=',y(itemp+kk), &
!                                  '   den0=',y(iden+kk)

!        else if (bbang) then
!         write(34,07) 'big bang:','   eta=',eta1, &
!                                  '   Nnu=',xnnu, &
!                                  '   H0=',hubble, &
!                                  '   Tcmb=',cmbtemp

!        else if (trho_hist) then
!         call update2(zero,tdum,ddum)
!         write(34,07) 'trho_hist:','  mass interior =',mint, &
!                                   '  shell mass =',mshell

!        end if

!        write(34,01) 'time',(ionam(nn), nn=kb,ke)

!        close(unit=34)
!  AT comments out TO.
 50     continue
       enddo

! end of the spatial loop
      enddo


!       write(6,*) 'wrote mass fraction headers'

! end of the file initialization
      end if

!      write(6,*) 'done with initialization'







! normal execution starts here

! for any time point

! for every spatial zone
      do k=1,max(1,nzone)
       kk = neqs*(k-1)

! open the files in append mode 

! energetics file
!  AT comments out FROM.
!       write(string,03) hfile,0,'_z',k,'.dat'
!       call sqeeze(string)
!       open (unit=22, file=string, status='old', position='append')


! thermodynamics file
!       write(string,03) hfile,1,'_z',k,'.dat'
!       call sqeeze(string)
!       open (unit=23, file=string, status='old', position='append')
!  AT comments out TO.


! form the mass fractions
       do j=1,ionmax
        xdum(j) = min(1.0d0,max(y(j+kk)*aion(j),1.0d-30))
       enddo


! mass conservation
       call indexx(ionmax,xdum,izwork2)      
       xx = 0.0d0
       do ii=1,ionmax
        xx = xx + xdum(izwork2(ii))
       end do
       xcons = 1.0d0 - xx


! get ye using normalized mass fractions
       xx = 1.0d0/xx
       xdum(1:ionmax) = min(1.0d0,max(xx*xdum(1:ionmax),1.0d-30))


! get abar, zbar and a few other composition variables
       call azbarapprox13(xdum(ionbeg),aion(ionbeg),zion(ionbeg),wion(ionbeg),ionmax, &
            zwork1(ionbeg),abar,zbar,wbar,yex,xcess)


! call an eos
       if (pure_network .eq. 0) then
        temp_row(1) = y(itemp+kk)*temp_scale
        den_row(1)  = y(iden+kk)*den_scale
       else
        temp_row(1) = btemp
        den_row(1)  = bden
       end if
       if (trho_hist) call update2(x,temp_row(1),den_row(1))
       abar_row(1) = abar
       zbar_row(1) = zbar
       jlo_eos = 1
       jhi_eos = 1
!       call helmeos
       call helmeosaprox13


! get the right hand sides, exact energy generation rate and so on

        call derivs(x,y,zwork1)
        if (pure_network .eq. 0) then
         ener = y(iener + kk) * energy_scale
         denerdt = zwork1(iener + kk) * energy_scale
        else
         ener = etot_row(1)
         denerdt = sdot - sneut
        end if



! figure some time scales
       call time_scales(temp_row(1),den_row(1),taud,tau_nse,tau_qse)


! compute the chemical potentials
!       do j=1,ionmax
!        chem_pot(j) = abar*((zion(j) - zbar)*deionz_row(1) &
!                          + (aion(j) - abar)*deiona_row(1))
!       end do
!       chem_sum = sum(chem_pot(1:ionmax)*zwork1(1:ionmax))



! and write what we found
!  AT comments out FROM.
!       write(22,05) kount,x,temp_row(1),den_row(1), &
!                    ener, sdot,sneut, denerdt, & 
!                    yex,xcons
!  AT comments out TO.


       if (iposx .eq. 0) then 
        posx = 0.0d0
       else
        posx = y(iposx+kk)
       end if
       if (ivelx .eq. 0) then 
        velx = 0.0d0
       else
        velx = y(ivelx+kk)
       end if

!  AT comments out FROM.
!       write(23,05) kount,x,posx,velx, &
!                   temp_row(1),den_row(1),ptot_row(1), &
!                   ener,stot_row(1),cs_row(1)

! close up the files
!       close(unit=22)
!       close(unit=23)
!  AT comments out TO.

! end of spatial loop
      end do

!      write(6,*) 'done with thermo file'




! for every spatial zone
      do k=1,max(1,nzone)
       kk = neqs*(k-1)

! write out the isotopic mass fractions in blocks of ncolumns
! lop is how many groups of ncolumns exist; jrem is the remainder
       lop  = ionmax/ncolumns
       jrem  = ionmax - ncolumns*lop
       do ilop = 1,lop+1
        kb = 1 + ncolumns*(ilop-1)
        ke = ncolumns + ncolumns*(ilop-1)
        if (ilop .eq. lop+1  .and. jrem .eq. 0) goto 70
        if (ilop .eq. lop+1) ke = ionmax


! open the output file in position mode
! abundance evolution file
!  AT comments out FROM.
!        write(string,03) hfile,ilop+1,'_z',k,'.dat'
!        call sqeeze(string)
!        open (unit=34, file=string, status='old', position='append')

!        write(34,04) kount,x,(y(nn+kk)*aion(nn), nn=kb,ke)
!!        write(34,04) kount,x,(y(nn+kk), nn=kb,ke)

!        close(unit=34)
!  AT comments out TO.
70      continue
       enddo


! end of spatial zone loop
      enddo

!      write(6,*) 'done with mass fractions file'



! all done
      return
      end




      subroutine net_tty(kount,x,y,idt,hdid)
      include 'implno.dek'
      include 'burn_common.dek'
      include 'network.dek'

! writes to the screen

! declare the pass
      integer :: kount,idt
      real*8  :: x,y(*),hdid

! declare local variables
      character*5 :: cdtname
      integer     :: ii
      real*8      :: xx,cons,abar,zbar,wbar,ye,xc,ttz,ddz


! form the mass fractions and nonconservation (cascading the sum)
      zwork2(1:ionmax) = y(1:ionmax) * aion(1:ionmax)
      call indexx(ionmax,zwork2,izwork2)      
      xx = 0.0d0
      do ii=1,ionmax
       xx = xx + zwork2(izwork2(ii))
      end do
!      xx   = sum(zwork2(1:ionmax))
      cons = 1.0d0 - xx

      call azbarapprox13(zwork2,aion,zion,wion,ionmax, &
           zwork1,abar,zbar,wbar,ye,xc)

      ttz = -1.0d0
      ddz = -1.0d0
      if (pure_network .eq. 0) then
       ttz = y(itemp)*temp_scale
       ddz = y(iden)*den_scale
      else
       ttz = btemp
       ddz = bden
      end if
      if (trho_hist) call update2(x,ttz,ddz)

      cdtname = 'time '
      if (idt .ne. 0) cdtname = ionam(idt)

      call indexx(ionmax,zwork2,izwork1)

      write(6,100) kount,' time',x, &
                   ' dt(',cdtname,')',hdid, &
                   (ionam(izwork1(ii)),zwork2(izwork1(ii)),ii=ionmax,ionmax-2,-1), &
                   ' temp',ttz,' den',ddz,' enuc',sdot-sneut, &
                   ' ye',ye,' sum',cons

100   format(1x,i6,' ',a,1pe11.4,a,a,a,1pe11.4,3(a6,1pe10.3),2(a,1pe9.2),1(a,1pe10.2),1(a,1pe9.2),1(a,1pe10.2))

     
      return
      end
!---------------------------------------------------------------------







!---------------------------------------------------------------------
! this file contains routines to write out the abundaces
! subroutine net_final_abund writes out the full isotopic abundance set
! subroutine net_final_element_abund writes out the elemental anbundances



      subroutine net_final_abund(xout)
      include 'implno.dek'
      include 'burn_common.dek'
      include 'network.dek'

! writes out the final composition

! declare the pass
      double precision xout(*)

! local variables
      character*80     final
      integer          i,lenstr

! popular format statements
 01   format(a,'final.dat')
 02   format(1x,i4,i4,1pe14.6,a6)
 07   format(1x,a,1pe14.6,a,1pe14.6)



! for the file name and open it
      write(final,01) hfile(1:lenstr(hfile,80))
      call sqeeze(final)
      open(unit=51,file=final,status='unknown')
      write(51,07) 'mass interior =',mint,'  shell mass =',mshell


! convert to integers
      izwork1(1:ionmax) = int(zion(1:ionmax))
      izwork2(1:ionmax) = int(aion(1:ionmax))

! write it out
      write(51,02) (izwork1(i), izwork2(i), xout(i), ionam(i), i=1,ionmax)


! close up shop
      close(unit=51)
      return
      end






      subroutine net_final_element_abund(xout)
      include 'implno.dek'
      include 'burn_common.dek'
      include 'network.dek'

! writes out the final composition

! declare the pass
      double precision xout(*)

! local variables
      character*80     final
      integer          i,j,lenstr
      double precision sum_x


! popular format statements
 01   format(a,'element_final.dat')
 02   format(1x,i4,1pe14.6,' ',a)
 07   format(1x,a,1pe14.6,a,1pe14.6)


! for the file name and open it
      write(final,01) hfile(1:lenstr(hfile,80))
      call sqeeze(final)
      open(unit=51,file=final,status='unknown')
      write(51,07) 'mass interior =',mint,'  shell mass =',mshell


! convert to integers
      izwork1(1:ionmax) = int(zion(1:ionmax))
      zwork1(1:ionmax)  = 0.0d0

! normal
! not sure why this doesn't work - indirect addressing?
!      zwork1(izwork1(1:ionmax)) = zwork1(izwork1(1:ionmax)) + xout(1:ionmax)

      do i=1,ionmax
       if (izwork1(i) .eq. 0) cycle
       zwork1(izwork1(i)) = zwork1(izwork1(i)) + xout(i)
      enddo 

      write(51,02) (i,zwork1(i),zsymb(i), i=1,zmax)




! close up shop
      close(unit=51)
      return
      end


!---------------------------------------------------------------------







!---------------------------------------------------------------------
      subroutine update2(tt,temp,den)
      include 'implno.dek'
      include 'network.dek'

! this routine evaluates the temperature temp and density den of
! as a function of time tt.

! declare the pass
      double precision tt,den,temp



! local variables, norder sets the order of the
! interpolation (2 points = linear, 3 = quadratic ...)

      integer          i,j,k,ntime,ntmax,jat,norder
      parameter        (ntmax=2000, norder=2)
      double precision ztime(ntmax),zden(ntmax),ztemp(ntmax),dy

      character*80     string,word
      integer          ipos,getnam,kiso
      double precision x,value

      integer          init
      data             init/0/


! stuff to do once
      if (init .eq. 0) then
       init = 1
       open (unit=17, file=trho_file, status='old')
       do i=1,ntmax
        read(17,*,end=10) ztime(i),ztemp(i),zden(i)
!        ztemp(i) = 1.0d9 * ztemp(i)
        ntime = ntime + 1
       enddo
       stop 'more than ntmax points in update2'
 10    close (unit=17)



! reset the zero point
!       zwork1(1) = ztime(1)
!       do i=1,ntime
!        ztime(i) = ztime(i) - zwork1(1)
!       end do

! store the beginning and end points in the coomon block zum array
       zwork1(1) = ztime(1)
       zwork1(2) = ztime(ntime)

      end if




! locate and interpolate to get the temperature and density
      if (tt .lt. ztime(1)) then
       temp = ztemp(1)
       den  = zden(1)
      else if (tt .gt. ztime(ntime)) then
       temp = ztemp(ntime)
       den  = zden(ntime)
      else
       call locate(ztime,ntime,tt,jat)
       jat = max(1,min(jat - norder/2 + 1,ntime - norder + 1))
       call up_polint(ztime(jat),ztemp(jat),norder,tt,temp,dy)
       call up_polint(ztime(jat),zden(jat),norder,tt,den,dy)
      end if

! bound the temperature, since a lot of rates go nuts
! above t9=10

      temp = max(1.0d7,min(temp,1.0d10))

      return
      end




      subroutine up_polint(xa,ya,n,x,y,dy)
      implicit none
      save


! given arrays xa and ya of length n and a value x, this routine returns a
! value y and an error estimate dy. if p(x) is the polynomial of degree n-1
! such that ya = p(xa) ya then the returned value is y = p(x)

! declare
      integer          n,nmax,ns,i,m
      parameter        (nmax=10)
      double precision xa(n),ya(n),x,y,dy,c(nmax),d(nmax),dif,dift, &
                       ho,hp,w,den

! find the index ns of the closest table entry; initialize the c and d tables
      ns  = 1
      dif = abs(x - xa(1))
      do i=1,n
       dift = abs(x - xa(i))
       if (dift .lt. dif) then
        ns  = i
        dif = dift
       end if
       c(i)  = ya(i)
       d(i)  = ya(i)
      enddo

! first guess for y
      y = ya(ns)

! for each column of the table, loop over the c's and d's and update them
      ns = ns - 1
      do m=1,n-1
       do i=1,n-m
        ho   = xa(i) - x
        hp   = xa(i+m) - x
        w    = c(i+1) - d(i)
        den  = ho - hp
        if (den .eq. 0.0) stop ' 2 xa entries are the same in polint'
        den  = w/den
        d(i) = hp * den
        c(i) = ho * den
       enddo

! after each column is completed, decide which correction c or d, to add
! to the accumulating value of y, that is, which path to take in the table
! by forking up or down. ns is updated as we go to keep track of where we
! are. the last dy added is the error indicator.
       if (2*ns .lt. n-m) then
        dy = c(ns+1)
       else
        dy = d(ns)
        ns = ns - 1
       end if
       y = y + dy
      enddo
      return
      end







      double precision function up_zbrent(func,x1,x2,tol,niter)
      implicit none
      save

! using brent's method this routine finds the root of a function func
! between the limits x1 and x2. the root is when accuracy is less than tol.
!
! note: eps the the machine floating point precision

! declare
      external          func
      integer           niter,itmax,iter
      parameter         (itmax = 100)
      double precision  func,x1,x2,tol,a,b,c,d,e,fa, &
                        fb,fc,xm,tol1,p,q,r,s,eps
      parameter         (eps=3.0d-15)

! initialize
      niter = 0
      a     = x1
      b     = x2
      fa    = func(a)
      fb    = func(b)
      if ( (fa .gt. 0.0  .and. fb .gt. 0.0)  .or. &
           (fa .lt. 0.0  .and. fb .lt. 0.0)       ) then
       write(6,100) x1,fa,x2,fb
 100   format(1x,' x1=',1pe11.3,' f(x1)=',1pe11.3,/, &
            1x,' x2=',1pe11.3,' f(x2)=',1pe11.3)
       stop 'root not bracketed in routine up_zbrent'
      end if
      c  = b
      fc = fb

! rename a,b,c and adjusting bound interval d
      do iter =1,itmax
       niter = niter + 1
       if ( (fb .gt. 0.0  .and. fc .gt. 0.0)  .or. &
            (fb .lt. 0.0  .and. fc .lt. 0.0)      ) then
        c  = a
        fc = fa
        d  = b-a
        e  = d
       end if
       if (abs(fc) .lt. abs(fb)) then
        a  = b
        b  = c
        c  = a
        fa = fb
        fb = fc
        fc = fa
       end if
       tol1 = 2.0d0 * eps * abs(b) + 0.5d0 * tol
       xm   = 0.5d0 * (c-b)

! convergence check
       if (abs(xm) .le. tol1 .or. fb .eq. 0.0) then
        up_zbrent = b
        return
       end if

! attempt quadratic interpolation
       if (abs(e) .ge. tol1 .and. abs(fa) .gt. abs(fb)) then
        s = fb/fa
        if (a .eq. c) then
         p = 2.0d0 * xm * s
         q = 1.0d0 - s
        else
         q = fa/fc
         r = fb/fc
         p = s * (2.0d0 * xm * q *(q-r) - (b-a)*(r - 1.0d0))
         q = (q - 1.0d0) * (r - 1.0d0) * (s - 1.0d0)
        end if

! check if in bounds
        if (p .gt. 0.0) q = -q
        p = abs(p)

! accept interpolation
        if (2.0d0*p .lt. min(3.0d0*xm*q - abs(tol1*q),abs(e*q))) then
         e = d
         d = p/q

! or bisect
        else
         d = xm
         e = d
        end if

! bounds decreasing to slowly use bisection
       else
        d = xm
        e = d
       end if

! move best guess to a
       a  = b
       fa = fb
       if (abs(d) .gt. tol1) then
        b = b + d
       else
        b = b + sign(tol1,xm)
       end if
       fb = func(b)
      enddo
      stop 'too many iterations in routine up_zbrent'
      end



!---------------------------------------------------------------------







!---------------------------------------------------------------------
! reaction rate library

! torch rates
! li7(t,n)   a(an,g)    be9(p,d)    be9(p,n)    b10(a,n)   b11(a,n)
! n14(p,a)   c11(p,g)   c12(a,n)    c13(a,n)    c13(p,n)   c14(a,g)
! c14(p,n)   c14(p,g)   o16(p,a)    n14(p,n)    n14(a,n)   n15(p,n)
! n15(a,n)   n15(a,g)   o14(a,g)    o17(a,g)    o17(a,n)   o18(a,g)
! o18(a,n)   ne20(p,a)  f18(p,g)    f19(p,g)    f19(p,n)   f19(a,p)
! na22(n,a)  ne20(p,g)  na23(p,a)   ne20(n,g)   ne21(p,g)  ne21(a,g)
! ne22(p,g)  ne22(a,g)  na22(n,p)   ne22(a,n)   na21(p,g)  mg24(p,a)
! ne21(a,n)  na22(p,g)  na23(p,g)   na23(p,n)   mg24(p,g)  al27(p,a)
! mg25(p,g)  mg25(a,p)  mg25(a,g)   mg25(a,n)   mg26(p,g)  mg26(a,g)
! mg26(a,n)  al25(p,g)  al26(p,g)   al27(a,n)   si27(p,g)  si28(p,g)
! si29(p,g)  si30(p,g)

! bigbang rates:
! n(e-nu)p   p(e-,nu)n  d(p,n)      d(n,g)      d(d,p)     d(d,n)
! t(p,n)     d(d,g)     t(p,g)      t(d,n)      t(t,2n)    he3(d,p)
! he3(t,d)   he3(t,np)  he4(np,g)   he4(d,g)    he4(t,n)   li6(p,he3)
! li6(n,g)   li7(d,n)   lit(t,2n)   li7(he3,np) li6(p,g)   li7(p,n)
! be7(d,p)   be7(t,np)  be7(3he,2p) li6(a,g)    li7(a,n)   be9(p,g)
! b10(p,a)   li7(a,g)   b11(p,a)    be7(a,g)    b11(p,n)   b8(a,p)
! b10(p,g)   c11(n,a)   be9(a,n)    b11(p,g)    b11(a,p)

! pp123 rates:
! p(p,e+nu)  p(n,g)     d(p,g)      he3(n,g)    he3+he3    he3(a,g)
! be7(e-,nu) be7(p,g)   li7(p,g)    li7(p,a)    b8(e+,nu)

! cno rates:
! c12(p,g)   n13(e-nu)  c13(p,g)    n14(p,g)    o15(e-nu)  n14(a,g)
! n15(p,g)   n15(p,a)   o16(p,g)    o17(p,a)    o17(p,g)   o18(p,a)
! o18(p,g)   f17(e-nu)  f18(e-nu)   f19(p,a)

! hot cno rates
! n13(p,g)   o14(e-nu)  o14(a,p)    o15(a,g)    f17(p,g)   ne18(e-nu)
! f18(p,a)   ne18(a,p)  ne19(p,g)   ne19(e-nu)  si26(a,p)

! alfa chain rates:
! a(aa,g)    c12(a,g)   c12+c12     c12+o16     o16+o16    o16(a,g)
! ne20(a,g)  ne20(a,g)  mg24(a,g)   mg24(a,p)   al27(p,g)  si28(a,g)
! si28(a,p)  p31(p,g)   s32(a,g)    s32(a,p)    cl35(p,g)  ar36(a,g)
! ar36(a,p)  k39(p,g)   ca40(a,g)   ca40(a,p)   sc43(p,g)  ti44(a,g)
! ti44(a,p)  v47(p,g)   cr48(a,g)   cr(a,p)     mn51(p,g)  fe52(a,g)
! fe52(a,p)  co55(p,g)

! photodisintegration rates:
! fe52(n,g) fe53(n,g)  fe54(p,g)    fe54(n,g)   fe54(a,p)   fe55(n,g)  






      subroutine tfactors(temp)
      include 'implno.dek'
      include 'tfactors.dek'

! sets various popular temperature factors into common block
! this routine must be called before any of the rates are called

! declare the pass
      double precision temp

! all these are in common block

      t9    = temp * 1.0d-9
      t92   = t9*t9
      t93   = t9*t92
      t94   = t9*t93
      t95   = t9*t94
      t96   = t9*t95

      t912  = sqrt(t9)
      t932  = t9*t912
      t952  = t9*t932
      t972  = t9*t952

      t913  = t9**oneth
      t923  = t913*t913
      t943  = t9*t913
      t953  = t9*t923
      t973  = t953*t923
      t9113 = t973*t943

      t914  = t9**(0.25d0)
      t934  = t914*t914*t914
      t954  = t9*t914
      t974  = t9*t934

      t915  = t9**onefif
      t935  = t915*t915*t915
      t945  = t915 * t935
      t965  = t9 * t915

      t916  = t9**onesix
      t976  = t9 * t916
      t9i76 = 1.0d0/t976

      t917  = t9**onesev
      t927  = t917*t917
      t947  = t927*t927

      t918  = sqrt(t914)
      t938  = t918*t918*t918
      t958  = t938*t918*t918

      t9i   = 1.0d0/t9
      t9i2  = t9i*t9i
      t9i3  = t9i2*t9i

      t9i12 = 1.0d0/t912
      t9i32 = t9i*t9i12
      t9i52 = t9i*t9i32
      t9i72 = t9i*t9i52

      t9i13 = 1.0d0/t913
      t9i23 = t9i13*t9i13
      t9i43 = t9i*t9i13
      t9i53 = t9i*t9i23

      t9i14 = 1.0d0/t914
      t9i34 = t9i14*t9i14*t9i14
      t9i54 = t9i*t9i14

      t9i15 = 1.0d0/t915
      t9i35 = t9i15*t9i15*t9i15
      t9i45 = t9i15 * t9i35
      t9i65 = t9i*t9i15

      t9i17 = 1.0d0/t917
      t9i27 = t9i17*t9i17
      t9i47 = t9i27*t9i27

      t9i18 = 1.0d0/t918
      t9i38 = t9i18*t9i18*t9i18
      t9i58 = t9i38*t9i18*t9i18

      return
      end





      subroutine rate_aan(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,bb,dbb,cc,dcc,dd,ddd

! he4(an,g)be9
      aa  = 1.0d0 + 0.344*t9
      bb  = t92 * aa
      dbb = 2.0d0 * t9 * aa + t92*0.344

      cc  = 1.0d0/bb
      dcc = -cc*cc*dbb

      dd  = 2.59e-6 * exp(-1.062*t9i)
      ddd = dd*1.062*t9i2

      term    = cc * dd
      dtermdt = dcc*dd + cc*ddd

! rates

      fr    = den * den * term
      dfrdt = den * den * dtermdt * 1.0d-9
      dfrdd = 2.0d0 * den * term

      rev      = 5.84e19 * t93 * exp(-18.260*t9i)
      drevdt   = rev*(3.0d0*t9i + 18.260*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_be9pd(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,q1
      parameter        (q1 = 1.0d0/0.2704d0)


! be9(p,d)be8 =>2a
      aa  = 2.11e+11 * t9i23 * exp(-10.359*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*10.359*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0  + 0.04*t913 + 1.09*t923 + 0.307*t9 &
            + 3.21*t943 + 2.30*t953
      dbb = oneth*0.04*t9i23 + twoth*1.09*t9i13 + 0.307 &
            + fourth*3.21*t913 + fiveth*2.30*t923

      cc  = 5.79e+08 * t9i * exp(-3.046*t9i)
      dcc = cc*(-t9i + 3.046*t9i2)

      dd  = 8.50e+08 * t9i34 * exp(-5.800*t9i)
      ddd = dd*(-0.75d0*t9i + 5.800*t9i2)

      term    = aa*bb + cc + dd
      dtermdt = daa*bb + aa*dbb + dcc + ddd

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 8.07e-11 * t9i32 *exp(-7.555*t9i)
      drevdt   = rev*(-1.5d0*t9i + 7.555*t9i2)

      rr    = den * den * rev * term
      drrdt = den * den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 2.0d0 * den * rev * term

      return
      end






      subroutine rate_be9pn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,bb,dbb,cc,dcc,zz,dzz


! be9(p,n)b9
      aa  = 5.58e7*(1.0d0 + 0.042*t912 + 0.985*t9)
      daa = 5.58e7*(0.5d0*0.042*t9i12 + 0.985)

      zz  = exp(-21.473*t9i)
      dzz = zz*21.473*t9i2

      bb  = aa * zz
      dbb = daa*zz + aa*dzz

      cc  = 1.02e+09 * t9i32 * exp(-26.725*t9i)
      dcc = cc*(-1.5d0*t9i + 26.725*t9i2)

      term    = bb + cc
      dtermdt = dbb + dcc

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term


      bb  = 0.998 * aa
      dbb = 0.998 * daa

      cc  = 0.998 * 1.02e+09 * t9i32 * exp(-5.252*t9i)
      dcc = cc*(-1.5d0*t9i + 5.252*t9i2)

      term    = bb + cc
      dtermdt = dbb + dcc

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end





      subroutine rate_b10an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,q1
      parameter        (q1 = 1.0d0/91.948921d0)


! b10(a,n)n13
      term    = 1.20e+13 * t9i23 * exp(-27.989*t9i13 - t92*q1)
      dtermdt = -twoth*term*t9i &
                + term*(oneth*27.989*t9i43 - 2.0d0*t9*q1)

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 9.34 * exp(-12.287*t9i)
      drevdt   = rev*12.287*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev*term

      return
      end





      subroutine rate_b11an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,q1
      parameter        (q1 = 1.0d0/0.0196d0)


! b11(a,n)n14
      aa  = 6.97e+12 * t9i23 * exp(-28.234*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*28.234*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.015*t913 + 8.115*t923 + 0.838*t9 &
            + 39.804*t943 + 10.456*t953
      dbb = oneth*0.015*t9i23 + twoth*8.115*t9i13 + 0.838 &
            + fourth*39.804*t913 + fiveth*10.456*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 1.79 * t9i32 * exp(-2.827*t9i)
      ddd  = dd*(-1.5d0*t9i + 2.827*t9i2)

      ee   = 1.71e+03 * t9i32 * exp(-5.178*t9i)
      dee  = ee*(-1.5d0*t9i + 5.178*t9i2)

      ff   = 4.49e+06 * t935 * exp(-8.596*t9i)
      dff  = ff*(0.6d0*t9i + 8.596*t9i2)

      term    = cc + dd + ee + ff
      dtermdt = dcc + ddd + dee + dff

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.67 * exp(-1.835*t9i)
      drevdt   = rev*1.835*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev*term

      return
      end





      subroutine rate_n14pa(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,bb,dbb,cc,dcc,dd,ddd, &
                       t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56,zz


! n14(p,a)b11
       aa     = 1.0d0 + 0.0478*t9
       bb     = aa**twoth
       dbb    = twoth*bb/aa*0.0478

       zz     = 1.0d0/bb
       cc     = aa + 7.56e-03*t953*zz
       dcc    = 0.0478 + (fiveth*7.56e-3*t923 - 7.56e-3*t953*zz*dbb)*zz

       zz     = 1.0d0/cc
       t9a    = t9*zz
       dt9a   = (1.0d0 - t9a*dcc)*zz

       zz      = dt9a/t9a
       t9a13   = t9a**oneth
       dt9a13  = oneth*t9a13*zz

       t9a56   = t9a**fivsix
       dt9a56  = fivsix*t9a56*zz

       dd      = 2.63e+16 * t9a56 * t9i32 * exp(-31.883/t9a13)
       ddd     = dd*(dt9a56/t9a56 - 1.5d0*t9i &
                 + 31.883/t9a13**2 * dt9a13)

       term    = dd * exp(-33.915*t9i)
       dtermdt = term*(ddd/dd + 33.915*t9i2)

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 0.272 * dd
      drevdt   = 0.272 * ddd

      rr    = den * rev
      drrdt = den * drevdt * 1.0d-9
      drrdd = rev

      return
      end





      subroutine rate_c11pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,q1
      parameter        (q1 = 1.0d0/2.647129d0)


! c11(p,g)n12
      aa  = 4.24e+04 * t9i23 * exp(-13.658*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*13.658*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0  + 0.031*t913 + 3.11*t923 + 0.665*t9 &
            + 4.61*t943 + 2.50*t953
      dbb = oneth*0.031*t9i23 + twoth*3.11*t9i13 + 0.665 &
            + fourth*4.61*t913 + fiveth*2.50*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 8.84e+03 * t9i32 * exp(-7.021*t9i)
      ddd  = dd*(-1.5d0*t9i + 7.021*t9i2)

      term    = cc + dd
      dtermdt = dcc + ddd

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 2.33e+10 * t932 * exp(-6.975*t9i)
      drevdt   = rev*(1.5d0*t9i + 6.975*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_c12an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,bb,dbb


! c12(a,n)o15
      aa  = 2.48e7 * (1.0d0 + 0.188*t912 + 0.015*t9)
      daa = 2.48e7 * (0.5d0*0.188*t9i12 + 0.015)

      bb  = exp(-98.661*t9i)
      dbb = bb*98.661*t9i2

      term    = aa * bb
      dtermdt = daa*bb + aa*dbb


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rr    = den * 1.41 * aa
      drrdt = den * 1.41 * daa  * 1.0d-9
      drrdd = 1.41 * aa

      return
      end





      subroutine rate_c13an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,gg,dgg,q1
      parameter        (q1 = 1.0d0/1.648656d0)


! c13(a,n)o16
      aa  = 6.77e+15 * t9i23 * exp(-32.329*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*32.329*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.013*t913 + 2.04*t923 + 0.184*t9
      dbb = oneth*0.013*t9i23 + twoth*2.04*t9i13 + 0.184

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 3.82e+05 * t9i32 * exp(-9.373*t9i)
      ddd  = dd*(-1.5d0*t9i + 9.373*t9i2)

      ee   = 1.41e+06 * t9i32 * exp(-11.873*t9i)
      dee  = ee*(-1.5d0*t9i + 11.873*t9i2)

      ff   = 2.0e+09 * t9i32 * exp(-20.409*t9i)
      dff  = ff*(-1.5d0*t9i + 20.409*t9i2)

      gg   = 2.92e+09 * t9i32 * exp(-29.283*t9i)
      dgg  = gg*(-1.5d0*t9i + 29.283*t9i2)

      term    = cc + dd + ee + ff + gg
      dtermdt = dcc + ddd + dee + dff + dgg

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev    = 5.79e+00 * exp(-25.711*t9i)
      drevdt = rev*25.711*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_c13pn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,bb,dbb


! c13(p,n)n13
      aa  = 1.88e+08*(1.0d0 - 0.167*t912 + 0.037*t9)
      daa = 1.88e+08*(0.037 - 0.5d0*0.167*t9i12)

      bb  = exp(-34.846*t9i)
      dbb = bb*34.846*t9i2

      term    = aa * bb
      dtermdt = daa*bb + aa*dbb

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rr    = den * 0.998 * aa
      drrdt = den * 0.998 * daa * 1.0d-9
      drrdd = 0.998 * aa

      return
      end






      subroutine rate_c14ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,q1
      parameter        (q1 = 1.0d0/7.086244d0)


! c14(a,g)o18
      aa  = 1.528e+09 * t9i23 * exp(-32.513*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*32.513*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.0128*t913 - 0.869*t923 - 0.0779*t9 &
            + 0.321*t943 + 0.0732*t953
      dbb = oneth*0.0128*t9i23 - twoth*0.869*t9i13 - 0.0779 &
            + fourth*0.321*t913 + fiveth*0.0732*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 3.375e+08 * t9i2 * exp(-32.513*t9i13)
      ddd  = dd*(-2.0d0*t9i + oneth*32.513*t9i43)

      ee   = 9.29e-08 * t9i32 * exp(-2.048*t9i)
      dee  = ee*(-1.5d0*t9i + 2.048*t9i2)

      ff   = 2.77e+03 * t9i45 * exp(-9.876*t9i)
      dff  = ff*(-0.8d0*t9i + 9.876*t9i2)

      term    = cc + dd + ee + ff
      dtermdt = dcc + ddd + dee + dff

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 5.42e+10 * t932 * exp(-72.262*t9i)
      drevdt   = rev*(1.5d0*t9i + 72.262*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_c14pn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,bb,dbb,cc,dcc,zz,dzz


! c14(p,n)n14
      aa  = 7.19e+05*(1.0d0 + 0.361*t912 + 0.502*t9)
      daa = 7.19e+05*(0.5d0*0.361*t9i12 + 0.502)

      zz  = exp(-7.263*t9i)
      dzz = zz*7.263*t9i2

      bb  = aa * zz
      dbb = daa*zz + aa*dzz

      cc  = 3.34e+08 * t9i12 * exp(-12.246*t9i)
      dcc = cc*(-0.5d0*t9i + 12.246*t9i2)

      term    = bb + cc
      dtermdt = dbb + dcc

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      cc  = 3.34e+08 * t9i12 * exp(-4.983*t9i)
      dcc = cc*(-0.5d0*t9i + 4.983*t9i2)

      rr    = den * 0.333 * (aa + cc)
      drrdt = den * 0.333 * (daa + dcc) * 1.0d-9
      drrdd = 0.333 * (aa + cc)

      return
      end






      subroutine rate_c14pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,q1
      parameter        (q1 = 1.0d0/32.729841d0)


! c14(p,g)n14
      aa  = 6.80e+06 * t9i23 * exp(-13.741*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*13.741*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.03*t913 + 0.503*t923 + 0.107*t9 &
            + 0.213*t943 + 0.115*t953
      dbb = oneth*0.03*t9i23 + twoth*0.503*t9i13 + 0.107 &
            + fourth*0.213*t913 + fiveth*0.115*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 5.36e+03 * t9i32 * exp(-3.811*t9i)
      ddd  = dd*(-1.5d0*t9i + 3.811*t9i2)

      ee   = 9.82e+04 * t9i13 * exp(-4.739*t9i)
      dee  = ee*(-oneth*t9i + 4.739*t9i2)

      term    = cc + dd + ee
      dtermdt = dcc + ddd + dee

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 9.00e+09 * t932 * exp(-118.452*t9i)
      drevdt   = rev*(1.5d0*t9i + 118.452*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_o16pa(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,bb,dbb,cc,dcc,dd,ddd, &
                       t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56,zz


! o16(p,a)n13
       aa     = 1.0d0 + 0.0776*t9
       bb     = aa**twoth
       dbb    = twoth*bb/aa*0.0776

       zz     = 1.0d0/bb
       cc     = aa + 0.0264*t953*zz
       dcc    = 0.0776 + (fiveth*0.0264*t923 - 0.0264*t953*zz*dbb)*zz

       zz     = 1.0d0/cc
       t9a    = t9*zz
       dt9a   = (1.0d0 - t9a*dcc)*zz

       zz      = dt9a/t9a
       t9a13   = t9a**oneth
       dt9a13  = oneth*t9a13*zz

       t9a56   = t9a**fivsix
       dt9a56  = fivsix * t9a56*zz

       dd      = 1.88e+18 * t9a56 * t9i32 * exp(-35.829/t9a13)
       ddd     = dd*(dt9a56/t9a56 - 1.5d0*t9i &
                 + 35.829/t9a13**2 * dt9a13)

       term    = dd * exp(-60.561*t9i)
       dtermdt = term*(ddd/dd + 60.561*t9i2)

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 0.172 * dd
      drevdt   = 0.172 * ddd

      rr    = den * rev
      drrdt = den * drevdt * 1.0d-9
      drrdd = rev

      return
      end





      subroutine rate_n14pn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,bb,dbb


! n14(p,n)o14
      aa  = 6.74e+07 * (1.0d0 + 0.658*t912 + 0.379*t9)
      daa = 6.74e+07 * (0.5d0*0.658*t9i12 + 0.379)

      bb  = exp(-68.762*t9i)
      dbb = bb*68.762*t9i2

      term    = aa * bb
      dtermdt = daa*bb + aa*dbb

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rr    = den * 2.99 * aa
      drrdt = den * 2.99 * daa * 1.0d-9
      drrdd = 2.99 * aa

      return
      end





      subroutine rate_n14an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,bb,dbb,cc,dcc,zz,dzz,q1
      parameter        (q1 = 1.0d0/7.828804d0)


! n14(a,n)f17
      aa  = 5.24e9*(1.0d0 - 1.15*t912 + 0.365*t9)
      daa = 5.24e9*(0.365 - 0.5d0*1.15*t9i12)

      zz  = exp(-t92*q1)
      dzz = -zz*2.0d0*t9*q1

      bb  = aa * zz
      dbb = daa*zz + aa*dzz

      cc   = 3.28e10 * t9i32 * exp(-1.5766e1*t9i)
      dcc  = cc*(-1.5d0*t9i + 1.5766e1*t9i2)

      term     = (bb + cc) * exp(-54.942*t9i)
      dtermdt  = term*((dbb+dcc)/(bb+cc) + 54.942*t9i2)

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      term     = 1.48 * (bb + cc)
      dtermdt  = 1.48 * (dbb + dcc)

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end




      subroutine rate_n15pn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,t9a,aa,daa,bb,dbb


! n15(p,n)o15
      t9a = min(t9,10.0d0)
      aa  = 3.51e+08 * (1.0d0 + 0.452*t912 - 0.191*t9a)
      if (t9a .eq. 10.0) then
       daa = 3.51e+08 * 0.5d0*0.452*t9i12
      else
       daa = 3.51e+08 * (0.5d0*0.452*t9i12 - 0.191)
      end if

      bb  = exp(-41.032*t9i)
      dbb = bb*41.032*t9i2

      term    = aa * bb
      dtermdt = daa*bb + aa*dbb

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      term    = 0.998 * aa
      dtermdt = 0.998 * daa

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end




      subroutine rate_n15an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,bb,dbb


! n15(a,n)f18
      aa  = 3.14e8 * (1.0d0 - 0.641*t912 + 0.108*t9)
      daa = 3.14e8 * (0.108 - 0.5d0*0.641*t9i12)

      bb  = exp(-74.479*t9i)
      dbb = bb*74.479*t9i2

      term     = aa * bb
      dtermdt  = daa*bb + aa*dbb

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      term     = 2.0d0 * aa
      dtermdt  = 2.0d0 * daa

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end





      subroutine rate_n15ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,q1
      parameter        (q1 = 1.0d0/0.379456d0)


! n15(a,g)f19
      aa  = 2.54e+10 * t9i23 * exp(-36.211*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*36.211*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.012*t913 + 1.69*t923 + 0.136*t9 &
            + 1.91*t943 + 0.391*t953
      dbb = oneth*0.012*t9i23 + twoth*1.69*t9i13 + 0.136 &
            + fourth*1.91*t913 + fiveth*0.391*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 9.83e-03 * t9i32 * exp(-4.232*t9i)
      ddd  = dd*(-1.5d0*t9i + 4.232*t9i2)

      ee   = 1.52e+03 * t9 * exp(-9.747*t9i)
      dee  = ee*(t9i + 9.747*t9i2)

      term    = cc + dd + ee
      dtermdt = dcc + ddd + dee

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev    = 5.54e+10 * t932 * exp(-46.578*t9i)
      drevdt = rev*(1.5d0*t9i + 46.578*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_o14ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,q1
      parameter        (q1 = 1.0d0/0.514089d0)


! o14(a,g)ne18
      aa  = 9.47e+08 * t9i23 * exp(-39.388*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*39.388*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.011*t913 + 1.974*t923 + 0.146*t9 &
            + 3.036*t943 + 0.572*t953
      dbb = oneth*0.011*t9i23 + twoth*1.974*t9i13 + 0.146 &
            + fourth*3.036*t913 + fiveth*0.572*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 1.16e-01 * t9i32 * exp(-11.733*t9i)
      ddd  = dd*(-1.5d0*t9i + 11.733*t9i2)

      ee   = 3.39e+01 * t9i32 * exp(-22.609*t9i)
      dee  = ee*(-1.5d0*t9i + 22.609*t9i2)

      ff   = 9.10e-03 * t95 * exp(-12.159*t9i)
      dff  = ff*(5.0d0*t9i + 12.159*t9i2)

      term    = cc + dd + ee + ff
      dtermdt = dcc + ddd + dee + dff

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev    = 5.42e+10 * t932 * exp(-59.328*t9i)
      drevdt = rev*(1.5d0*t9i + 59.328*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_o17ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56, &
                       ft9a,dft9a,fpt9a,dfpt9a,gt9x,dgt9x,zz


! o17(a,g)ne21
       aa    = 1.0d0 + 0.1646*t9
       zz    = 1.0d0/aa
       t9a   = t9*zz
       dt9a  = (1.0d0 - t9a*0.1646)*zz

       zz     = dt9a/t9a
       t9a13  = t9a**oneth
       dt9a13 = oneth*t9a13*zz

       t9a56  = t9a**fivsix
       dt9a56 = fivsix * t9a56*zz

       aa     = 0.786/t9a
       daa    = -aa*zz
       bb     = aa**3.51
       dbb    = 3.51*bb/aa * daa
       ft9a   = exp(-bb)
       dft9a  = -ft9a*dbb

       aa     = t9a/1.084
       bb     = aa**1.69
       dbb    = 1.69*bb/aa * dt9a/1.084
       fpt9a  = exp(-bb)
       dfpt9a = -fpt9a*dbb

       aa     = oneth*exp(-10.106*t9i)
       daa    = aa*10.106*t9i2
       gt9x   = 1.0d0 + aa
       dgt9x  = daa

       zz     = 1.0d0/gt9x
       aa     = 1.73e17 * fpt9a*zz
       daa    = (1.73e17*dfpt9a - aa*dgt9x)*zz

       bb     = 3.50e15 * ft9a*zz
       dbb    = (3.50e15*dft9a - bb*dgt9x)*zz

       term    = (aa+bb) * t9a56 * t9i32 * exp(-39.914/t9a13)
       dtermdt = term*((daa+dbb)/(aa+bb) + dt9a56/t9a56 &
                 - 1.5d0*t9i + 39.914/t9a13**2 * dt9a13)

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 8.63e+10 * t932 * exp(-85.305*t9i)
      drevdt   = rev*(1.5d0*t9i + 85.305*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_o17an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,bb,dbb,cc,dcc,dd, &
                       t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56,gt9x,dgt9x,zz


! o17(a,n)ne20
       aa     = 1.0d0 + 0.0268*t9
       bb     = aa**twoth
       dbb    = twoth*bb/aa*0.0268

       zz     = 1.0d0/bb
       cc     = aa + 0.0232*t953*zz
       dcc    = 0.0268 + (fiveth*0.0232*t923 - 0.0232*t953*zz*dbb)*zz

       zz     = 1.0d0/cc
       t9a    = t9*zz
       dt9a   = (1.0d0 - t9a*dcc)*zz

       zz      = dt9a/t9a
       t9a13   = t9a**oneth
       dt9a13  = oneth*t9a13*zz

       t9a56   = t9a**fivsix
       dt9a56  = fivsix * t9a56*zz

       dd     = oneth*exp(-10.106*t9i)
       gt9x   = 1.0d0 + dd
       dgt9x  = dd*10.106*t9i2

       term      = 1.03e+18/gt9x * t9a56 * t9i32 * exp(-39.914/t9a13)
       dtermdt   = term*(-dgt9x/gt9x + dt9a56/t9a56 &
                   - 1.5d0*t9i + 39.914/t9a13**2 * dt9a13)

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.86e+01 * exp(-6.852*t9i)
      drevdt   = rev*6.852*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_o18ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,gg,dgg,hh,dhh,theta,q1
      parameter        (theta = 0.1d0, &
                        q1    = 1.0d0/0.117649d0)


! o18(a,g)ne22
! giessen et al 1994 nuc phys a 567, 146 for t9 less than 0.3
! cf88 otherwise

      if (t9.lt.0.3) then
       aa   = 1.066d-41 * t9i32 * exp(-5.507d-01*t9i)
       daa  = aa*(-1.5d0*t9i + 5.507d-1*t9i2)

       bb   = 1.852d-13 * t9i32 * exp(-2.070*t9i)
       dbb  = bb*(-1.5d0*t9i + 2.070*t9i2)

       cc   = 1.431d-02 * t9i32 * exp(-4.462*t9i)
       dcc  = cc*(-1.5d0*t9i + 4.462*t9i2)

       dd   = 2.055d-04 * t9i32 * exp(-5.374*t9i)
       ddd  = dd*(-1.5d0*t9i + 5.374*t9i2)

       ee   = 5.332d+00 * t9i32 * exp(-6.285*t9i)
       dee  = ee*(-1.5d0*t9i + 6.285*t9i2)

       ff   = 1.457d+00 * t9i32 * exp(-7.121*t9i)
       dff  = ff*(-1.5d0*t9i + 7.121*t9i2)

       gg   = 3.121d-02 * t9i32 * exp(-7.292*t9i)
       dgg  = gg*(-1.5d0*t9i + 7.292*t9i2)

       hh   = 6.23d+03 * t9 * exp(-16.987*t9i)
       dhh  = hh*(t9i + 16.987*t9i2)

       term    = aa + bb + cc + dd + ee + ff + gg + hh
       dtermdt = daa + dbb + dcc + ddd + dee + dff + dgg + dhh


      else
       aa  = 1.82d+12 * t9i23 * exp(-40.057*t9i13 - t92*q1)
       daa = aa*(-twoth*t9i + oneth*40.057*t9i43 - 2.0d0*t9*q1)

       bb  = 1.0d0 + 0.01*t913 + 0.988*t923 + 0.072*t9 &
             + 3.17*t943 + 0.586*t953
       dbb = oneth*0.01*t9i23 + twoth*0.988*t9i13 + 0.072 &
            + fourth*3.17*t913 + fiveth*0.586*t923

       cc   = aa * bb
       dcc  = daa*bb + aa*dbb

       dd   = 7.54 * t9i32 * exp(-6.228*t9i)
       ddd  = dd*(-1.5d0*t9i + 6.228*t9i2)

       ee   = 34.8 * t9i32 * exp(-7.301*t9i)
       dee  = ee*(-1.5d0*t9i + 7.301*t9i2)

       ff   = 6.23d+03 * t9 * exp(-16.987*t9i)
       dff  = ff*(t9i + 16.987*t9i2)

       gg   = theta * 1.0d-11 * t9i32 * exp(-1.994*t9i)
       dgg  = gg*(-1.5d0*t9i + 1.994*t9i2)

       term    = cc + dd + ee + ff + gg
       dtermdt = dcc + ddd + dee + dff + dgg
      end if

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 5.85d+10 * t932 * exp(-112.208*t9i)
      drevdt   = rev*(1.5d0*t9i + 112.208*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end







      subroutine rate_o18an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,bb,dbb,cc,dcc,dd, &
                       ee,dee,ff,dff,gg,dgg,hh,dhh,ft9a,dft9a,gt9,dgt9, &
                       t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56,gt9i,zz


! o18(a,n)ne21
      aa     = 1.0d0 + 0.0483*t9
      bb     = aa**twoth
      dbb    = twoth*bb/aa*0.0483

      zz     = 1.0d0/bb
      cc     = aa + 0.00569*t953*zz
      dcc    = 0.0483 + (fiveth*0.00569*t923 - 0.00569*t953*zz*dbb)*zz

      zz     = 1.0d0/cc
      t9a    = t9*zz
      dt9a   = (1.0d0 - t9a*dcc)*zz

      zz      = dt9a/t9a
      t9a13   = t9a**oneth
      dt9a13  = oneth*t9a13*zz

      t9a56   = t9a**fivsix
      dt9a56  = fivsix * t9a56*zz

      dd     = 5.0d0 * exp(-23.002*t9i)
      gt9    = 1.0d0 + dd
      gt9i   = 1.0d0/gt9
      dgt9   = dd*23.002*t9i2

      ee     = 0.431/t9a
      dee    = -ee*zz
      ff     = ee**3.89
      dff    = 3.89*ff/ee*dee
      ft9a   = exp(-ff)
      dft9a  = -ft9a*dff

      gg     = 7.22e+17 * ft9a*gt9i * t9a56 * t9i32 * exp(-40.056/t9a13)
      dgg    = gg*(-dff - gt9i*dgt9 + dt9a56/t9a56 - 1.5d0*t9i &
               + 40.056/t9a13**2 *dt9a13)

      hh     = 150.31 / gt9 * exp(-8.045*t9i)
      dhh    = hh*(-gt9i*dgt9 + 8.045*t9i2)

      term     = gg + hh
      dtermdt  = dgg + dhh

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

! must protect the 8.045*t9i from overflow, so write it this way

      gg     = 7.22e+17*gt9i * t9a56 * t9i32 &
               * exp(-ff - 40.056/t9a13 + 8.045*t9i)
      dgg    = gg*(-gt9i*dgt9 + dt9a56/t9a56 - 1.5d0*t9i &
                - dff + 40.056/t9a13**2*dt9a13 - 8.045*t9i2)

      hh     = 150.31 * gt9i
      dhh    = -hh*gt9i*dgt9

      term    = 0.784 * (gg + hh)
      dtermdt = 0.784 * (dgg + dhh)

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end




      subroutine rate_ne20pa(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,bb,dbb,cc,dcc,dd,ddd, &
                       ee,dee,ff,dff,t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56, &
                       zz,t9b


! ne20(p,a)f17
       aa     = 1.0d0 + 0.0612*t9
       bb     = aa**twoth
       dbb    = twoth*bb/aa*0.0612

       zz     = 1.0d0/bb
       cc     = aa + 0.013*t953*zz
       dcc    = 0.0612 + (fiveth*0.013*t923 - 0.013*t953*zz*dbb)*zz

       zz     = 1.0d0/cc
       t9a    = t9*zz
       dt9a   = (1.0d0 - t9a*dcc)*zz

       zz      = dt9a/t9a
       t9a13   = t9a**oneth
       dt9a13  = oneth*t9a13*zz

       t9a56   = t9a**fivsix
       dt9a56  = fivsix * t9a56*zz

       t9b     = min(t9,10.0d0)
       dd      = 5.31 + 0.544*t9b - 0.0523*t9b*t9b
       ddd     = 0.544 - 2.0d0*0.0523*t9b
       if (t9b .eq. 10.0) ddd     = 0.0d0

       ee      = 3.25e19 * dd * t9a56 * t9i32 * exp(-43.176/t9a13)
       dee     = ee*(ddd/dd + dt9a56/t9a56 - 1.5d0*t9i &
                 +   43.176/t9a13**2 * dt9a13)

       ff      = exp(-47.969*t9i)
       dff     = ff*47.969*t9i2

       term    = ee * ff
       dtermdt = dee*ff + ee*dff


! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 0.0537 * ee
      drevdt   = 0.0537 * dee

      rr    = den * rev
      drrdt = den * drevdt * 1.0d-9
      drrdd = rev

      return
      end





      subroutine rate_f18pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,gg,dgg


! f18(p,g)ne19
! wiescher and kettner, apj 263, 891 1982
      aa  = 1.658e+7 * t9i23 * exp(-18.06*t9i13)
      daa = aa*(-twoth*t9i + oneth*18.06*t9i43)

      bb  = 4.604 + 0.106*t913 + 0.053*t923 + 0.009*t9 &
            - 0.036*t943 - 0.015*t953
      dbb = oneth*0.106*t9i23 + twoth*0.053*t9i13 + 0.009 &
            - fourth*0.036*t913 - fiveth*0.015*t923

! for temps greater than about t9 = 20, bb goes negative
      if (bb .le. 0.0) then
       bb = 0.0d0
       dbb = 0.0d0
      end if

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 4.55e-14 * t9i32* exp(-0.302*t9i)
      ddd  = dd*(-1.5d0*t9i + 0.302*t9i2)

      ee   = 327.0 * t9i32 * exp(-3.84*t9i)
      dee  = ee*(-1.5d0*t9i + 3.84*t9i2)

      ff   = 1.32e+04 * t9i32 * exp(-5.22*t9i)
      dff  = ff*(-1.5d0*t9i + 5.22*t9i2)

      gg   = 93.0 * t9i32 * exp(-4.29*t9i)
      dgg  = gg*(-1.5d0*t9i + 4.29*t9i2)

      term    = cc + dd + ee + ff + gg
      dtermdt = dcc + ddd + dee + dff + dgg

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 2.73e+10 * t932 * exp(-74.396*t9i)
      drevdt   = rev*(1.5d0*t9i + 74.396*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_f19pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,gg,dgg,hh,dhh,zz,q1
      parameter        (q1 = 1.0d0/0.173056d0)


! f19(p,g)ne20
      aa  = 6.04e+07 * t9i23 * exp(-18.113*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*18.113*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.023*t913 + 2.06*t923 + 0.332*t9 &
            + 3.16*t943 + 1.30*t953
      dbb = oneth*0.023*t9i23 + twoth*2.06*t9i13 + 0.332 &
            + fourth*3.16*t913 + fiveth*1.30*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 6.32e+02 * t9i32 * exp(-3.752*t9i)
      ddd  = dd*(-1.5d0*t9i + 3.752*t9i2)

      ee   = 7.56e+04 * t9i27 * exp(-5.722*t9i)
      dee  = ee*(-twosev*t9i + 5.722*t9i2)

      ff   = 7.0*exp(-16.44*t9i)
      dff  = ff*16.44*t9i2

      gg   = 4.0 * exp(-2.09*t9i)
      dgg  = gg*2.09*t9i2

      hh   = 1.0d0 + ff + gg
      dhh  = dff + dgg

      zz      = 1.0d0/hh
      term    = (cc + dd + ee)*zz
      dtermdt = (dcc + ddd + dee - term*dhh)*zz

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.7e+10 * t932 * exp(-149.093*t9i)
      drevdt   = rev*(1.5d0*t9i + 149.093*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_f19pn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,bb,dbb


! f19(p,n)ne19
      aa  = 1.27e+08 * (1.0d0 - 0.147*t912 + 0.069*t9)
      daa = 1.27e+08 * (0.069 - 0.5d0*0.147*t9i12)

      bb  = exp(-46.659*t9i)
      dbb = bb*46.659*t9i2

      term    = aa * bb
      dtermdt = daa*bb + aa*dbb

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      term    = 0.998 * aa
      dtermdt = 0.998 * daa

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end





      subroutine rate_f19ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,q1
      parameter        (q1 = 1.0d0/0.405769)


! f19(a,p)ne22
      aa  = 4.50e+18 * t9i23 * exp(-43.467*t9i13  - t92*q1)
      daa = -twoth*aa*t9i + aa*(oneth*43.467*t9i43 - 2.0d0*t9*q1)

      bb   = 7.98e+04 * t932 * exp(-12.760*t9i)
      dbb  = 1.5d0*bb*t9i + bb*12.760*t9i2

      term    = aa + bb
      dtermdt = daa + dbb

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 6.36 * exp(-19.439*t9i)
      drevdt   = rev*19.439*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_na22na(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa, &
                       t9b,t9b2,t9b3


! na22(n,a)f19
      t9b  = min(t9,10.0d0)
      t9b2 = t9b*t9b
      t9b3 = t9b2*t9b
      aa  = 1.0d0 + 0.8955*t9b - 0.05645*t9b2 + 7.302e-04*t9b3
      daa = 0.8955 - 2.0d0*0.05645*t9b + 3.0d0*7.302e-4*t9b2
      if (t9b .eq. 10.0) daa = 0.0d0

      term     = 1.21e6 * exp(aa)
      dtermdt  = term*daa

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.10 * exp(-22.620*t9i)
      drevdt   = rev*22.620*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_ne20pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ff,gg,dgg,zz


! ne20(p,g)na21
      aa  = 9.55e+06 * exp(-19.447*t9i13)
      daa = aa*oneth*19.447*t9i43

      bb  = 1.0d0 + 0.0127*t9i23
      dbb = -twoth*0.0127*t9i53

      cc  = t92 * bb * bb
      dcc = 2.0d0*cc*t9i + 2.0d0*t92*bb*dbb

      zz  = 1.0d0/cc
      dd  = aa*zz
      ddd = (daa - dd*dcc)*zz

      aa  = 2.05e+08 * t9i23 * exp(-19.447*t9i13)
      daa = aa*(-twoth*t9i + oneth*19.447*t9i43)

      bb  = sqrt (t9/0.21)
      dbb = 0.5d0/(bb * 0.21)

      cc  = 2.67 * exp(-bb)
      dcc = -cc*dbb

      ff  = 1.0d0 + cc

      gg  = aa*ff
      dgg = daa*ff + aa*dcc


      aa  = 18.0 * t9i32 * exp(-4.242*t9i)
      daa = aa*(-1.5d0*t9i + 4.242*t9i2)

      bb  = 10.2 * t9i32 * exp(-4.607*t9i)
      dbb = bb*(-1.5d0*t9i + 4.607*t9i2)

      cc  = 3.6e+04 * t9i14 * exp(-11.249*t9i)
      dcc = cc*(-0.25d0*t9i + 11.249*t9i2)

      term    = dd + gg + aa + bb + cc
      dtermdt = ddd + dgg + daa + dbb + dcc


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 4.63e+09 * t932 * exp(-28.216*t9i)
      drevdt   = rev*(1.5d0*t9i + 28.216*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_na23pa(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,gg,dgg,hh,dhh,theta,q1,q2
      parameter        (theta = 0.1d0, &
                        q1    = 1.0d0/0.0169d0, &
                        q2    = 1.0d0/0.017161d0)


! na23(p,a)ne20
! el eid & champagne 1995
      if (t9 .le. 2.0) then
       aa  = 1.26d+10 * t9i23 * exp(-20.758*t9i13 - t92*q1)
       daa = -twoth*aa*t9i + aa*(oneth*20.758*t9i43 - 2.0d0*t9*q1)

       bb  = 1.0d0  + 0.02*t913 - 13.8*t923 - 1.93*t9 &
             + 234.0*t943 + 83.6*t953
       dbb = oneth*0.02*t9i23 - twoth*13.8*t9i13 - 1.93 &
             + fourth*234.0*t913 + fiveth*83.6*t923

       cc   = aa * bb
       dcc  = daa*bb + aa*dbb

       dd   = 4.38 * t9i32 * exp(-1.979*t9i)
       ddd  = -1.5d0*dd*t9i + dd*1.979*t9i2

       ee   = 6.50d+06 * (t9**(-1.366d0)) * exp(-6.490*t9i)
       dee  = -1.366d0*ee*t9i + ee*6.490*t9i2

       ff   = 1.19d+08 * (t9**(-1.055d0)) * exp(-11.411*t9i)
       dff  = -1.055d0*ff*t9i + ff*11.411*t9i2

       gg   = theta * 9.91d-14 * t9i32 * exp(-0.418*t9i)
       dgg  = -1.5d0*gg*t9i + gg*0.418*t9i2

       term    = cc + dd + ee + ff + gg
       dtermdt = dcc + ddd + dee + dff + dgg



! cf88 + one term from gorres, wiesher & rolfs 1989, apj 343, 365
      else
       aa  = 8.56d+09 * t9i23 * exp(-20.766*t9i13 - t92*q2)
       daa = -twoth*aa*t9i + aa*(oneth*20.766*t9i43 - 2.0d0*t9*q2)

       bb  = 1.0d0  + 0.02*t913 + 8.21*t923 + 1.15*t9 &
             + 44.36*t943 + 15.84*t953
       dbb = oneth*0.02*t9i23 + twoth*8.21*t9i13 + 1.15 &
             + fourth*44.36*t913 + fiveth*15.84*t923

       cc   = aa * bb
       dcc  = daa*bb + aa*dbb

       dd   = 4.02 * t9i32 * exp(-1.99*t9i)
       ddd  = -1.5d0*dd*t9i + dd*1.99*t9i2

       ee   = 1.18d+04 * t9i54 * exp(-3.148*t9i)
       dee  = -1.25d0*ee*t9i + ee*3.148*t9i2

       ff   = 8.59d+05 * t943 * exp(-4.375*t9i)
       dff  = fourth*ff*t9i + ff*4.375*t9i2

       gg   = theta * 3.06d-12 * t9i32 * exp(-0.447*t9i)
       dgg  = -1.5d0*gg*t9i + gg*0.447*t9i2

       hh   = theta * 0.820 * t9i32 * exp(-1.601*t9i)
       dhh  = -1.5d0*hh*t9i + hh*1.601*t9i2

       term    = cc + dd + ee + ff + gg + hh
       dtermdt = dcc + ddd + dee + dff + dgg + dhh
      end if


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.25 * exp(-27.606*t9i)
      drevdt   = rev*27.606*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end






      subroutine rate_ne20ng(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc


! ne20(n,g)ne21
! wm88 Apj 239, 943; fit over range of experimental data, constant otherwise

      if (t9 .lt. 5.8025d-2) then
       term    = 5.449d+03
       dtermdt = 0.0d0

      else if (t9 .gt. 1.1605) then
       term    = 6.977d+04
       dtermdt = 0.0d0

      else if (t9 .ge. 5.8025d-2 .and. t9 .le. 2.9012d-1) then
       term    = 4.7219d+3 + 2.5248d+4*t9 - 2.7448d+5*t92 &
                 + 9.2848d+5*t93
       dtermdt = 2.5248d+4 - 2.0d0*2.7448d+5*t9 &
                 + 3.0d0*9.2848d+5*t92

      else

       aa  = 1.802d+04 * (t9/0.348)**4.43
       daa = 4.43 * aa * t9i

       bb  = -5.931 * (t9-0.348) + 1.268 * (t9-0.348)**2
       dbb = -5.931 + 2.0d0*1.268*(t9 - 0.348)

       cc  = exp(bb)
       dcc = cc*dbb

       term    = aa * cc
       dtermdt = daa*cc + aa*dcc

      end if

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      =  4.650d+09 * t932 * exp(-78.46*t9i)
      drevdt   = rev*(1.5d0*t9i + 78.46*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_ne21pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,gg,dgg,hh,dhh,xx,dxx,zz, &
                       theta,q1
      parameter        (theta = 0.1d0, &
                        q1    = 1.0d0/0.003364d0)


! ne21(p,g)na22

! el eid & champagne 1995

      if (t9.le.2.0) then
       aa  = 3.4d+08 * t9i23 * exp(-19.41*t9i13)
       daa = aa*(-twoth*t9i + oneth*19.41*t9i43)

       bb  = (16.7*t9 - 1.0)**2
       dbb = 2.0d0*(16.7*t9 - 1.0)*16.7

       cc  = 0.56 * exp(-bb)
       dcc = -cc*dbb

       dd  = 1.0d0 + cc
       ddd = dcc

       ee   = aa * dd
       dee  = daa*dd + aa*ddd

       ff   = 6.12 * t9i32 * exp(-1.403*t9i)
       dff  = ff*(-1.5d0*t9i + 1.403*t9i2)

       gg   = 1.35d+04 * t9i32 * exp(-3.008*t9i)
       dgg  = gg*(-1.5d0*t9i + 3.008*t9i2)

       aa   = t9**0.67
       daa  = 0.67*aa*t9i
       zz   = 1.0d0/aa

       hh   = 3.12d+06 * t9**(-0.72) * exp(-8.268*zz)
       dhh  = hh*(-0.72d0*t9i + 8.268*zz*zz*daa)

       xx   = theta * 1.1d-03 * t9i32 * exp(-1.114*t9i)
       dxx  = xx*(-1.5d0*t9i + 1.114*t9i2)

       term    = ee + ff + gg + hh + xx
       dtermdt = dee + dff + dgg + dhh + dxx


! cf88
      else

       aa  = theta * 2.95d+08 * t9i23 * exp(-19.462*t9i13 -t92*q1)
       daa = aa*(-twoth*t9i + oneth*19.462*t9i43 - 2.0d0*t9*q1)

       bb  = 1.0d0 + 0.021*t913 + 13.29*t923 + 1.99*t9 &
             + 124.1*t943 + 47.29*t953
       dbb = oneth*0.021*t9i23 + twoth*13.29*t9i13 + 1.99 &
             + fourth*124.1*t913 + fiveth*47.29*t923

       cc   = aa * bb
       dcc  = daa*bb + aa*dbb

       dd   = theta * 7.80d-01 * t9i32 * exp(-1.085*t9i)
       ddd  = dd*(-1.5d0*t9i + 1.085*t9i2)

       ee   = 4.37d+08 * t9i23 * exp(-19.462*t9i13)
       dee  = ee*(-twoth*t9i + oneth*19.462*t9i43)

       ff   = 5.85 * t9i32 * exp(-1.399*t9i)
       dff  = ff*(-1.5d0*t9i + 1.399*t9i2)

       gg   = 1.29d+04 * t9i32 * exp(-3.009*t9i)
       dgg  = gg*(-1.5d0*t9i + 3.009*t9i2)

       hh   = 3.15d+05 * t9i35 * exp(-5.763*t9i)
       dhh  = hh*(-0.6d0*t9i + 5.763*t9i2)

       term    = cc + dd + ee + ff + gg + hh
       dtermdt = dcc + ddd + dee + dff + dgg + dhh
      end if


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.06d+10 * t932 * exp(-78.194*t9i)
      drevdt   = rev*(1.5d0*t9i + 78.194*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_ne21ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,gg,dgg,hh,dhh,zz, &
                       t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56


! ne21(a,g)mg25
       aa    = 1.0d0 + 0.0537*t9
       zz    = 1.0d0/aa

       t9a   = t9*zz
       dt9a  = (1.0d0 - t9a*0.0537)*zz

       zz     = dt9a/t9a
       t9a13  = t9a**oneth
       dt9a13 = oneth*t9a13*zz

       t9a56  = t9a**fivsix
       dt9a56 = fivsix * t9a56*zz

       aa     = 8.72e-03*t9 - 6.87e-04*t92 + 2.15e-05*t93
       daa    = 8.72e-3 - 2.0d0*6.87e-4*t9 + 3.0d0*2.15e-5*t92

       bb     = 1.52e-04 * exp(-46.90*t9i13*aa)
       dbb    = bb*46.90*(oneth*t9i43*aa - t9i13*daa)

       cc     = 1.5*exp(-4.068*t9i)
       dcc    =  cc*4.068*t9i2

       gg     = 2.0 * exp(-20.258*t9i)
       dgg    = gg*20.258*t9i2

       hh     = 1.0d0 + cc + gg
       dhh    = dcc + dgg

       zz     = 1.0d0/hh
       dd     = bb*zz
       ddd    = (dbb - dd*dhh)*zz

       aa     = 4.94e+19 * t9a56 * t9i32 * exp(-46.89/t9a13)
       daa    = aa*(dt9a56/t9a56 - 1.5d0*t9i &
                    + 46.89/t9a13**2 * dt9a13)

       bb     =  2.66e+07 * t9i32 * exp(-22.049*t9i)
       dbb    = bb*(-1.5d0*t9i + 22.049*t9i2)

       cc     = aa + bb
       dcc    = daa + dbb

       term    = dd * cc
       dtermdt = ddd*cc + dd*dcc


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev    = 4.06e+10 * t932 * exp(-114.676*t9i)
      drevdt = rev*(1.5d0*t9i + 114.676*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end







      subroutine rate_ne21an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,zz, &
                       t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56


! ne21(a,n)mg24
       aa    = 1.0d0 + 0.0537*t9
       zz    = 1.0d0/aa

       t9a   = t9*zz
       dt9a  = (1.0d0 - t9a*0.0537)*zz

       zz     = dt9a/t9a
       t9a13  = t9a**oneth
       dt9a13 = oneth*t9a13*zz

       t9a56  = t9a**fivsix
       dt9a56 = fivsix * t9a56*zz

       aa     = 4.94e+19 * t9a56 * t9i32 * exp(-46.89/t9a13)
       daa    = aa*(dt9a56/t9a56 - 1.5d0*t9i &
                    + 46.89/t9a13**2 * dt9a13)

       bb     =  2.66e+07 * t9i32 * exp(-22.049*t9i)
       dbb    = bb*(-1.5d0*t9i + 22.049*t9i2)

       cc     = 2.0d0*exp(-20.258*t9i)
       dcc    = cc*20.258*t9i2

       dd     = 1.5*exp(-4.068*t9i)
       ddd    = dd*4.068*t9i2

       ee     = 1.0d0 + cc + dd
       dee    = dcc + ddd

       zz      = 1.0d0/ee
       term    = (aa + bb)*zz
       dtermdt = (daa + dbb - term*dee)*zz


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 12.9 * exp(-29.606*t9i)
      drevdt   = rev*29.606*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end






      subroutine rate_ne22pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,gg,dgg,theta
      parameter        (theta = 0.1d0)


! ne22(p,g)na23

! el eid & champagne 1995

      if (t9.le.2.0) then
       aa  = 1.05d+09 * t9i23 * exp(-19.431*t9i13)
       daa = aa*(-twoth*t9i + oneth*19.431*t9i43)

       bb  = 1.24d-09 * t9i32 * exp(-0.414*t9i)
       dbb = bb*(-1.5d0*t9i + 0.414*t9i2)

       cc  = 2.90d-02 * t9i32 * exp(-1.752*t9i)
       dcc = cc*(-1.5d0*t9i + 1.752*t9i2)

       dd  = 9.30d+04 * t9**(-1.174) * exp(-5.100*t9i)
       ddd = dd*(-1.174*t9i + 5.100*t9i2)

       ee   = 5.71d+05 * t9**(0.249) * exp(-7.117*t9i)
       dee  = ee*(0.249*t9i + 7.117*t9i2)

       ff   = theta * 3.25d-04 * t9i32 * exp(-0.789*t9i)
       dff  = ff*(-1.5d0*t9i + 0.789*t9i2)

       gg   = theta * 0.10 * t9i32 * exp(-1.161*t9i)
       dgg  = gg*(-1.5d0*t9i + 1.161*t9i2)

       term    = aa + bb + cc + dd + ee + ff + gg
       dtermdt = daa + dbb + dcc + ddd + dee + dff + dgg


! cf88
      else

       aa  = 1.15d+09 * t9i23 * exp(-19.475*t9i13)
       daa = aa*(-twoth*t9i + oneth*19.475*t9i43)

       bb  = 9.77d-12 * t9i32 * exp(-0.348*t9i)
       dbb = bb*(-1.5d0*t9i + 0.348*t9i2)

       cc   = 8.96d+03 * t9i32 * exp(-4.84*t9i)
       dcc  = cc*(-1.5d0*t9i + 4.84*t9i2)

       dd   = 6.52d+04 * t9i32 * exp(-5.319*t9i)
       ddd  = dd*(-1.5d0*t9i + 5.319*t9i2)

       ee   = 7.97d+05 * t9i12 * exp(-7.418*t9i)
       dee  = ee*(-0.5d0*t9i + 7.418*t9i2)

       ff   = theta * 1.63d-01 * t9i32 * exp(-1.775*t9i)
       dff  = ff*(-1.5d0*t9i + 1.775*t9i2)

       term    = aa + bb + cc + dd + ee + ff
       dtermdt = daa + dbb + dcc + ddd + dee + dff

      end if


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 4.67d+09 * t932 * exp(-102.048*t9i)
      drevdt   = rev*(1.5d0*t9i + 102.048*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_ne22ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,res1,dres1, &
                       ft9a,dft9a,fpt9a,dfpt9a,gt9x,dgt9x, &
                       t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56, &
                       rdmass,res2,zz
      parameter        (rdmass = 22.0d0*4.0d0/26.0d0, &
                        res2   = -11.604d0 * 22.0d0/26.0d0)


! ne22(a,g)mg26
! kappeler 1994 apj 437, 396

      if (t9 .lt. 1.25) then

       res1 = 1.54d-01*(t9*rdmass)**(-1.5)
       dres1 = -1.5d0 * res1 * t9i

       aa    = 1.7d-36 * res1 * exp(res2*t9i*0.097)
       daa   = aa/res1*dres1 - aa*res2*0.097*t9i2

       bb    = 1.5d-7 * res1 * exp(res2*t9i*0.400)
       dbb   = bb/res1*dres1 - bb * res2 * 0.400 * t9i2

       cc    = 0.5 * res1 * 3.7d-2 * exp(res2*t9i*0.633)
       dcc   = cc/res1*dres1 - cc*res2*0.633*t9i2

       dd    = res1 * 3.6d+1 * exp(res2*t9i*0.828)
       ddd   = dd/res1*dres1 - dd*res2*0.828*t9i2

       term    = aa + bb + cc + dd
       dtermdt = daa + dbb + dcc + ddd


! cf88
      else
       aa    = 1.0d0 + 0.0548*t9
       zz    = 1.0d0/aa

       t9a   = t9*zz
       dt9a  = (1.0d0 - t9a*0.0548)*zz

       zz     = dt9a/t9a
       t9a13  = t9a**oneth
       dt9a13 = oneth*t9a13*zz

       t9a56  = t9a**fivsix
       dt9a56 = fivsix * t9a56*zz

       aa     = 0.197/t9a
       daa    = -aa*zz
       bb     = aa**4.82
       dbb    = 4.82*bb/aa * daa
       ft9a   = exp(-bb)
       dft9a  = -ft9a*dbb

       aa     = t9a/0.249
       bb     = aa**2.31
       dbb    = 2.31*bb/aa * dt9a/0.249
       fpt9a  = exp(-bb)
       dfpt9a = -fpt9a*dbb

       aa     = 5.0d0*exp(-14.791*t9i)
       daa    = aa*14.791*t9i2
       gt9x   = 1.0d0 + aa
       dgt9x  = daa

       zz     = 1.0d0/gt9x
       aa     = 4.16e19 * fpt9a*zz
       daa    = (4.16e19*dfpt9a - aa*dgt9x)*zz

       bb     = 2.08e16 * ft9a*zz
       dbb    = (2.08e16*dft9a - bb*dgt9x)*zz

       term    = (aa+bb) * t9a56 * t9i32 * exp(-47.004/t9a13)
       dtermdt = term*((daa+dbb)/(aa+bb) + dt9a56/t9a56 &
                       - 1.5d0*t9i + 47.004/t9a13**2 * dt9a13)
      end if


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev    = 6.15d+10 * t932 * exp(-123.151*t9i)
      drevdt = rev*(1.5d0*t9i + 123.151*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_na22np(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa


! na22(n,p)ne22
      aa  = 1.0d0 - 3.037e-02*t9 + 8.380e-03*t92 - 7.101e-04*t93
      daa =  -3.037e-02 + 2.0d0*8.380e-03*t9 - 3.0d0*7.101e-04*t92

      term    = 1.24e+08 * exp(aa)
      dtermdt = term*daa

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev     = 7.01*exp(-42.059*t9i)
      drevdt  = rev*42.059*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_ne22an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,gg,dgg,ft9a,dft9a,gt9x,dgt9x, &
                       t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56,res1,res2, &
                       zz
      parameter        (res1 = 2.4731857150793075d-2, &
                        res2 = -9.8187694549560547d0)

! note: res1=1.54d-1*(88./26.)**(-1.5)   res2=-11.604*(22./26.)

! ne22(a,n)mg25
! kappeler 1994 apj 437, 396 ; wiescher suggest only 828 kev, ignore 633 kev

      if (t9 .lt. 0.6) then
       term = res1*1.64d+02 * t9i32 * exp(t9i*0.828*res2)
       dtermdt = -1.5d0*term*t9i - term*res2*0.828*t9i2

! cf88
      else
       aa    = 1.0d0 + 0.0548*t9
       zz    = 1.0d0/aa

       t9a   = t9*zz
       dt9a  = (1.0d0 - t9a*0.0548)*zz

       zz     = dt9a/t9a
       t9a13  = t9a**oneth
       dt9a13 = oneth*t9a13*zz

       t9a56  = t9a**fivsix
       dt9a56 = fivsix * t9a56*zz

       aa     = 0.197/t9a
       daa    = -aa*zz
       bb     = aa**4.82
       dbb    = 4.82*bb/aa * daa
       ft9a   = exp(-bb)
       dft9a  = -ft9a*dbb

       gg     = bb
       dgg    = dbb

       aa     = 5.0d0*exp(-14.791*t9i)
       daa    = aa*14.791*t9i2
       gt9x   = 1.0d0 + aa
       dgt9x  = daa

       zz     = 1.0d0/gt9x
       aa     = ft9a*zz
       daa    = (dft9a - aa*dgt9x)*zz

       bb     = 4.16e+19 * t9a56 * t9i32 * exp(-47.004/t9a13)
       dbb    = bb*(dt9a56/t9a56 - 1.5d0*t9i &
                    + 47.004/t9a13**2 * dt9a13)

       cc     = aa*bb
       dcc    = daa*bb + aa*dbb

       dd      = 1.44e-04*zz * exp(-5.577*t9i)
       ddd     = -dd*zz*dgt9x + dd*5.577*t9i2

       term    = cc + dd
       dtermdt = dcc + ddd
      end if


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev    = 7.833d-5
      drevdt = 0.0d0
      if (t9 .gt. 0.008) then
       rev    = 0.544 * exp(5.577*t9i)
       drevdt = -rev*5.577*t9i2
      end if

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end






      subroutine rate_na21pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,q1
      parameter        (q1 = 1.0d0/0.133956d0)


! na21(p,g)mg22
      aa  = 1.41e+05 * t9i23 * exp(-20.739*t9i13 -  t92*q1)
      daa = aa*(-twoth*t9i + oneth*20.739*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.020*t913 + 4.741*t923 + 0.667*t9 &
            + 16.380*t943 + 5.858*t953
      dbb = oneth*0.020*t9i23 + twoth*4.741*t9i13 + 0.667 &
            + fourth*16.380*t913 + fiveth*5.858*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 6.72e+02 * t9i34 * exp(-2.436*t9i)
      ddd  = dd*(-0.75d0*t9i + 2.436*t9i2)

      term    = cc + dd
      dtermdt = dcc + ddd

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 7.44e+10 * t932 * exp(-63.790*t9i)
      drevdt   = rev*(1.5d0*t9i + 63.790*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_mg24pa(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,gg,t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56,zz


! mg24(p,a)na21
       aa     = 1.0d0 + 0.127*t9
       zz     = 1.0d0/aa

       t9a    = t9*zz
       dt9a   = (1.0d0 - t9a*0.127)*zz

       zz      = dt9a/t9a
       t9a13   = t9a**oneth
       dt9a13  = oneth*t9a13*zz

       t9a56   = t9a**fivsix
       dt9a56  = fivsix * t9a56*zz

       gg      = min(t9,12.0d0)
       aa      = 4.43 + 3.31*gg - 0.229*gg*gg
       daa     = 3.31 - 2.0d0*0.229*gg
       if (gg .eq. 12.0) daa = 0.0d0

       bb      = 1.81e21 * t9a56 * t9i32 * exp(-49.967/t9a13)
       dbb     = bb*(dt9a56/t9a56 - 1.5d0*t9i &
                     + 49.967/t9a13**2 * dt9a13)

       cc      = aa*bb
       dcc     = daa*bb + aa*dbb

       dd      = exp(-79.843*t9i)
       ddd     = dd*79.843*t9i2

       term    = cc * dd
       dtermdt = dcc*dd + cc*ddd

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 0.0771 * cc
      drevdt   = 0.0771 * dcc

      rr    = den * rev
      drrdt = den * drevdt * 1.0d-9
      drrdd = rev

      return
      end






      subroutine rate_na22pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb


! na22(p,g)mg23
      aa  = 9.63e-05 * t932 * exp(-0.517*t9i)
      daa = aa*(1.5d0*t9i + 0.517*t9i2)

      bb  = 2.51e+04 * t9 * exp(-2.013*t9i)
      dbb = bb*(t9i + 2.013*t9i2)

      term    = aa + bb
      dtermdt = daa + dbb

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.27e+10 * t932 * exp(-87.933*t9i)
      drevdt   = rev*(1.5d0*t9i + 87.933*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_na23pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,gg,dgg,hh,hhi,xx,dxx, &
                       theta,q1
      parameter        (theta = 0.1d0, &
                        q1    = 1.0d0/0.088209d0)


! na23(p,g)mg24

! el eid & champagne 1995
      if (t9 .le. 2.0) then
       aa  = 2.47d+09 * t9i23 * exp(-20.758*t9i13)
       daa = aa*(-twoth*t9i + oneth*20.758*t9i43)

       bb  = 9.19d+01 * t9i32 * exp(-2.789*t9i)
       dbb = bb*(-1.5d0*t9i + 2.789*t9i2)

       cc  = 1.72d+04 * t9i32 * exp(-3.433*t9i)
       dcc = cc*(-1.5d0*t9i + 3.433*t9i2)

       dd  = 3.44d+04 * t9**0.323 * exp(-5.219*t9i)
       ddd = dd*(0.323*t9i + 5.219*t9i2)

       ee   = theta * 2.34d-04 * t9i32 * exp(-1.590*t9i)
       dee  = ee*(-1.5d0*t9i + 1.590*t9i2)

       term    = aa + bb + cc + dd + ee
       dtermdt = daa + dbb + dcc + ddd + dee


! cf88 + gorres, wiesher & rolfs 1989, apj 343, 365
      else

       aa  = 2.93d+08 * t9i23 * exp(-20.766*t9i13 - t92*q1)
       daa = aa*(-twoth*t9i + oneth*20.766*t9i43 - 2.0d0*t9*q1)

       bb  = 1.0d0 + 0.02*t913 + 1.61*t923 + 0.226*t9 &
             + 4.94*t943 + 1.76*t953
       dbb = oneth*0.02*t9i23 + twoth*1.61*t9i13 + 0.226 &
            + fourth*4.94*t913 + fiveth*1.76*t923

       xx  = aa * bb
       dxx = daa*bb + aa*dbb

       cc   = 9.34d+01 * t9i32 * exp(-2.789*t9i)
       dcc  = cc*(-1.5d0*t9i + 2.789*t9i2)

       dd   = 1.89d+04 * t9i32 * exp(-3.434*t9i)
       ddd  = dd*(-1.5d0*t9i + 3.434*t9i2)

       ee   = 5.1d+04 * t915 * exp(-5.51*t9i)
       dee  = ee*(0.2d0*t9i + 5.51*t9i2)

       ff   = theta * 0.820 * t9i32 * exp(-1.601*t9i)
       dff  = ff*(-1.5d0*t9i + 1.601*t9i2)

       gg   = 1.5 * exp(-5.105*t9i)
       dgg  = gg*5.105*t9i2

       hh   = 1.0d0 + gg
       hhi  = 1.0d0/hh

       term    = (xx + cc + dd + ee + ff) * hhi
       dtermdt = (dxx + dcc + ddd + dee + dff - term*dgg)*hhi

      end if


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 7.49d+10 * t932 * exp(-135.665*t9i)
      drevdt   = rev*(1.5d0*t9i + 135.665*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_na23pn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,bb,dbb,cc,dcc, &
                       t9a,dt9a,t9a32,dt9a32,zz


! na23(p,n)mg24
      aa  = 1.0d0 + 0.141*t9
      zz  = 1.0d0/aa

      t9a = t9*zz
      dt9a = (1.0d0 - t9a*0.141)*zz

      aa    = sqrt(t9a)
      t9a32 = t9a * aa
      dt9a32 = 1.5d0 * aa * dt9a

      bb   = 9.29d8 * (1.0d0 - 0.881d0 * t9a32 * t9i32)
      dbb  = -9.29d8 * 0.881d0 * t9i32*(dt9a32 - 1.5d0*t9a32*t9i)

      cc   = exp(-56.173*t9i)
      dcc  = cc*56.173*t9i2

      term    = bb * cc
      dtermdt = dbb*cc + bb*dcc

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      term    = 0.998 * bb
      dtermdt = 0.998 * dbb

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end






      subroutine rate_mg24pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,gg,ggi


! mg24(p,g)al25
      aa  = 5.60e+08 * t9i23 * exp(-22.019*t9i13)
      daa = aa*(-twoth*t9i + oneth*22.019*t9i43)

      bb  = 1.0d0 + 0.019*t913 - 0.173*t923 - 0.023*t9
      dbb = oneth*0.019*t9i23 - twoth*0.173*t9i13 - 0.023

! stop negative rates above t9 = 10
      if (bb .le. 0.0) then
       bb  = 0.0d0
       dbb = 0.0d0
      end if

      cc  = aa * bb
      dcc = daa*bb + aa*dbb

      dd   = 1.48e+03 * t9i32 * exp(-2.484*t9i)
      ddd  = dd*(-1.5d0*t9i + 2.484*t9i2)

      ee   = 4.00e+03 * exp(-4.180*t9i)
      dee  = ee*4.180*t9i2

      ff   = 5.0 * exp(-15.882*t9i)
      dff  = ff*15.882*t9i2

      gg   = 1.0d0 + ff
      ggi  = 1.0d0/gg

      term    = (cc + dd + ee) * ggi
      dtermdt = (dcc + ddd + dee - term*dff)*ggi


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.13e+09 * t932 * exp(-26.358*t9i)
      drevdt   = rev*(1.5d0*t9i + 26.358*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_al27pa(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,gg,dgg,theta
      parameter        (theta = 0.1d0)


! al27(p,a)mg24
! champagne 1996
      aa  = 4.71d+05 * t9i23 * exp(-23.25*t9i13 - 3.57*t92)
      daa = -twoth*aa*t9i + aa*(oneth*23.25*t9i43 - 2.0d0*3.57*t9)

      bb  = 1.0d0 + 0.018*t913 - 7.29*t923 - 0.914*t9 &
            + 77.2*t943 + 24.6*t953
      dbb = oneth*0.018*t9i23 - twoth*7.29*t9i13 - 0.914 &
            + fourth*77.2*t913 + fiveth*24.6*t923

      cc  = aa * bb
      dcc = daa*bb + aa*dbb

      dd  = 2.23d+04 * (t9**3.989) * exp(-2.148 * t9**(-1.293))
      ddd = 3.989*dd*t9i + 1.293*dd*2.148*t9**(-2.293)

      ee   = 0.17 * 1.29d-09 * t9i32 * exp(-0.836*t9i)
      dee  = -1.5d0*ee*t9i + ee*0.836*t9i2

      ff   = theta * 2.73d-03 * t9i32 * exp(-2.269*t9i)
      dff  = -1.5d0*ff*t9i + ff*2.269*t9i2

      gg   = theta * 2.60d-02 * t9i32 * exp(-2.492*t9i)
      dgg  = -1.5d0*gg*t9i + gg*2.492*t9i2

      term    = cc + dd + ee + ff + gg
      dtermdt = dcc + ddd + dee + dff + dgg

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.81*exp(-18.572*t9i)
      drevdt   = rev*18.572*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_mg25pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,gg,dgg,q1,q2
      parameter        (q1 = 1.0d0/0.0036d0, &
                        q2 = 1.0d0/169.0d0)


! mg25(p,g)al26
      aa  = 3.57e+09 * t9i23 * exp(-22.031*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*22.031*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.019*t913 + 7.669*t923 + 1.015*t9 &
            + 167.4*t943 + 56.35*t953
      dbb = oneth*0.019*t9i23 + twoth*7.669*t9i13 + 1.015 &
            + fourth*167.4*t913 + fiveth*56.35*t923

      cc  = aa * bb
      dcc = daa*bb + aa*dbb

      dd  = 3.07e-13 * t9i32 * exp(-0.435*t9i)
      ddd = dd*(-1.5d0*t9i + 0.435*t9i2)

      ee   = 1.94e-07 * t9i32 * exp(-0.673*t9i)
      dee  = ee*(-1.5d0*t9i + 0.673*t9i2)

      ff   = 3.15e-05 * t9**(-3.40)* exp(-1.342*t9i - t92*q2)
      dff  = ff*(-3.40d0*t9i + 1.342*t9i2 - 2.0d0*t9*q2)

      gg   = 1.77e+04 * t958 * exp(-3.049*t9i - t92*q2)
      dgg  = gg*(0.625*t9i + 3.049*t9i2 - 2.0d0*t9*q2)

      term    = cc + dd + ee + ff + gg
      dtermdt = dcc + ddd + dee + dff + dgg

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.03e+10 * t932 * exp(-73.183*t9i)
      drevdt   = rev*(1.5d0*t9i + 73.183*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_mg25ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc


! mg25(a,p)si28
      aa  = -23.271*t9i13 + 6.46*t9 - 2.39*t92 + 0.506*t93 &
            - 6.04e-2*t94 + 3.75e-3*t95 - 9.38e-5*t96

      daa = oneth*23.271*t9i43 + 6.46 - 2.0d0*2.39*t9 + 3.0d0*0.506*t92 &
            - 4.0d0*6.04e-2*t93 + 5.0d0*3.75e-3*t94 - 6.0d0*9.38e-5*t95

      bb  = 3.23e8 * t9i23 * exp(aa)
      dbb  = -twoth*bb*t9i + bb*daa

! dbb/bb
      cc   = -twoth*t9i + daa

      term    = bb * exp(-13.995*t9i)
      dtermdt = term*cc + term*13.995*t9i2

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 2.86 * bb
      drevdt   = 2.86 * dbb

      rr    = den * rev
      drrdt = den * drevdt * 1.0d-9
      drrdd = rev

      return
      end





      subroutine rate_mg25ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,gt9x,dgt9x,t9a,dt9a,t9a13,dt9a13, &
                       t9a56,dt9a56,zz


! mg25(a,g)si29
      aa    = 1.0d0 + 0.0630*t9
      zz    = 1.0d0/aa

      t9a   = t9*zz
      dt9a  = (1.0d0 - t9a*0.0630)*zz

      zz     = dt9a/t9a
      t9a13  = t9a**oneth
      dt9a13 = oneth*t9a13*zz

      t9a56  = t9a**fivsix
      dt9a56 = fivsix * t9a56*zz

      aa     = oneth*10.0d0*exp(-13.180*t9i)
      daa    = aa*13.180*t9i2
      gt9x   = 1.0d0 + aa
      dgt9x  = daa

      bb     = 1.0d0/gt9x
      dbb    = -bb*bb*dgt9x

      cc     = 3.59e+20 * bb * t9a56 * t9i32 * exp(-53.41/t9a13)
      dcc    = cc*(dbb*gt9x + dt9a56/t9a56 - 1.5d0*t9i &
                   + 53.41/t9a13**2 * dt9a13)

      dd     = 0.0156*t9 - 1.79e-03*t92 + 9.08e-05*t93
      ddd    = 0.0156 - 2.0d0*1.79e-03*t9 + 3.0d0*9.08e-05*t92

      ee     = 5.87e-04*exp(-53.42*t9i13*dd)
      dee    = ee*53.42*(oneth*t9i43*dd - t9i13*ddd)

      term    = cc * ee
      dtermdt = dcc*ee + cc*dee

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev    = 1.90e+11 * t932 * exp(-129.128*t9i)
      drevdt = rev*(1.5d0*t9i + 129.128*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_mg25an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,zz, &
                       gt9x,dgt9x,t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56


! mg25(a,n)si28
      aa    = 1.0d0 + 0.0630*t9
      zz    = 1.0d0/aa

      t9a   = t9*zz
      dt9a  = (1.0d0 - t9a*0.0630)*zz

      zz     = dt9a/t9a
      t9a13  = t9a**oneth
      dt9a13 = oneth*t9a13*zz

      t9a56  = t9a**fivsix
      dt9a56 = fivsix * t9a56*zz

      aa     = oneth*10.0d0*exp(-13.180*t9i)
      daa    = aa*13.180*t9i2
      gt9x   = 1.0d0 + aa
      dgt9x  = daa

      bb     = 1.0d0/gt9x
      dbb    = -bb*bb*dgt9x

      term    = 3.59e+20 * bb * t9a56 * t9i32 * exp(-53.41/t9a13)
      dtermdt = term*(dbb*gt9x + dt9a56/t9a56 - 1.5d0*t9i &
                + 53.41/t9a13**2 * dt9a13)

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 20.0*exp(-30.792*t9i)
      drevdt   = rev*30.792*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end




      subroutine rate_mg26pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,gg,dgg,hh,dhh,theta
      parameter        (theta = 0.1d0)


! mg26(p,g)al27
! champagne 1996
      aa  = 8.54d-12 * t9i32 * exp(-0.605*t9i)
      daa = aa*(-1.5d0*t9i + 0.605*t9i2)

      bb  = 2.75d-06 * t9i32 * exp(-1.219*t9i)
      dbb = bb*(-1.5d0*t9i + 1.219*t9i2)

      cc  = 1.30d-02 * t9i32 * exp(-1.728*t9i)
      dcc = cc*(-1.5d0*t9i + 1.728*t9i2)

      dd  = 8.06d+00 * t9i32 * exp(-2.537*t9i)
      ddd = dd*(-1.5d0*t9i + 2.537*t9i2)

      ee   = 1.45d+03 * t9i32 * exp(-3.266*t9i)
      dee  = ee*(-1.5d0*t9i + 3.266*t9i2)

      ff   = 4.03d+04 * t9i32 * exp(-3.784*t9i)
      dff  = ff*(-1.5d0*t9i + 3.784*t9i2)

      gg   = 8.82d+04 * t9**(-0.21) * exp(-4.194*t9i)
      dgg  = gg*(-0.21*t9i + 4.194*t9i2)

      hh   = theta * 1.93d-05 * t9i32 * exp(-1.044*t9i)
      dhh  = hh*(-1.5d0*t9i + 1.044*t9i2)

      term    = aa + bb + cc + dd + ee + ff + gg + hh
      dtermdt = daa + dbb + dcc + ddd + dee + dff + dgg + dhh


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.14d+09 * t932 * exp(-95.99*t9i)
      drevdt   = rev*(1.5d0*t9i + 95.99*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_mg26ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,gt9x,dgt9x,t9a,dt9a,t9a13,dt9a13, &
                       t9a56,dt9a56,zz


! mg26(a,g)si30
      aa    = 1.0d0 + 0.0628*t9
      zz    = 1.0d0/aa

      t9a   = t9*zz
      dt9a  = (1.0d0 - t9a*0.0628)*zz

      zz     = dt9a/t9a
      t9a13  = t9a**oneth
      dt9a13 = oneth*t9a13*zz

      t9a56  = t9a**fivsix
      dt9a56 = fivsix * t9a56*zz

      aa     = 5.0d0*exp(-20.990*t9i)
      daa    = aa*20.990*t9i2
      gt9x   = 1.0d0 + aa
      dgt9x  = daa

      bb     = 1.0d0/gt9x
      dbb    = -bb*bb*dgt9x

      cc     = 2.93e+20 * bb * t9a56 * t9i32 * exp(-53.505/t9a13)
      dcc    = cc*(dbb*gt9x + dt9a56/t9a56 - 1.5d0*t9i &
                   + 53.505/t9a13**2 * dt9a13)

      dd     = 0.0751*t9 - 0.0105*t92 + 5.57e-04*t93
      ddd    = 0.0751 - 2.0d0*0.0105*t9 + 3.0d0*5.57e-04*t92

      ee     = 4.55e-2 * exp(-53.51*t9i13*dd)
      dee    = ee*53.51*(oneth*t9i43*dd - t9i13*ddd)

      term    = cc * ee
      dtermdt = dcc*ee + cc*dee

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev    = 6.38e+10 * t932 * exp(-123.52*t9i)
      drevdt = rev*(1.5d0*t9i + 123.52*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_mg26an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,zz, &
                       gt9x,dgt9x,t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56


! mg26(a,n)si29
      aa    = 1.0d0 + 0.0628*t9
      zz    = 1.0d0/aa

      t9a   = t9*zz
      dt9a  = (1.0d0 - t9a*0.0628)*zz

      zz     = dt9a/t9a
      t9a13  = t9a**oneth
      dt9a13 = oneth*t9a13*zz

      t9a56  = t9a**fivsix
      dt9a56 = fivsix * t9a56*zz

      aa     = 5.0d0*exp(-20.990*t9i)
      daa    = aa*20.990*t9i2
      gt9x   = 1.0d0 + aa
      dgt9x  = daa

      bb     = 1.0d0/gt9x
      dbb    = -bb*bb*dgt9x

      term   = 2.93e+20 * bb * t9a56 * t9i32 * exp(-53.505/t9a13)
      dtermdt= term*(dbb*gt9x + dt9a56/t9a56 - 1.5d0*t9i &
               + 53.505/t9a13**2 * dt9a13)

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.68*exp(-0.401*t9i)
      drevdt   = rev*0.401*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_al25pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee

! al25(p,g)si26
! coc et al 1995 a&a 299, 479 , case b

      aa  = 8.98d+1 * t9i32 * exp(-4.874*t9i)
      daa = aa*(-1.5d0*t9i + 4.874*t9i2)

      bb  = 1.568d+3 * t9i32 * exp(-9.632*t9i)
      dbb = bb*(-1.5d0*t9i + 9.632*t9i2)

      cc  = 2.42d+8 * t9i23 * exp(-23.18*t9i13)
      dcc = cc*(-twoth*t9i + oneth*23.18*t9i43)

      dd  = 4.10d-02 * t9i32 * exp(-1.741*t9i)
      ddd = dd*(-1.5d0*t9i + 1.741*t9i2)

      ee  = 2.193d+3 * t9i32 * exp(-4.642*t9i)
      dee = ee*(-1.5d0*t9i + 4.642*t9i2)

      term    = aa + bb + cc + dd + ee
      dtermdt = daa + dbb + dcc + ddd + dee


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.117d+11 * t932 * exp(-64.048*t9i)
      drevdt   = rev*(1.5d0*t9i + 64.048*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end







      subroutine rate_al26pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,theta
      parameter        (theta = 0.1d0)


! al26(p,g)si27
! coc et al 1995 a&a 299, 479

      aa  = 1.53d+9 * t9**(-1.75) * exp(-23.19*t9i13)
      daa = aa*(-1.75*t9i + oneth*23.19*t9i43)

      bb  = theta*8.7d-7 * t9i32 * exp(-0.7845*t9i)
      dbb = bb*(-1.5d0*t9i + 0.7845*t9i2)

      cc  = theta*1.00d-3 * t9i32 * exp(-1.075*t9i)
      dcc = cc*(-1.5d0*t9i + 1.075*t9i2)

      dd  = 9.00d+00 * t9i32 * exp(-2.186*t9i)
      ddd = dd*(-1.5d0*t9i + 2.186*t9i2)

      ee  = 5.05d+02 * t9i32 * exp(-3.209*t9i)
      dee = ee*(-1.5d0*t9i + 3.209*t9i2)

      ff  = 9.45d+03 * t9i * exp(-4.008*t9i)
      dff = ff*(-t9i + 4.008*t9i2)

      term    = aa + bb + cc + dd + ee + ff
      dtermdt = daa + dbb + dcc + ddd + dee + dff


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.46d+10 * t932 * exp(-86.621*t9i)
      drevdt   = rev*(1.5d0*t9i + 86.621*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_al27an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb


! al27(a,n)p30

      aa    = 8.2e+04*exp(-30.588*t9i)
      daa   = aa*30.588*t9i2

      bb    = 5.21e+05 * t974 * exp(-33.554*t9i)
      dbb   = 1.75d0*bb*t9i + bb*33.554*t9i2

      term    = aa + bb
      dtermdt = daa + dbb

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term


      aa  = 5.21e+05 * t974 * exp(-2.966*t9i)
      daa = aa*(1.75d0*t9i + 2.966*t9i2)

      rev      = 6.75d0 * (8.20e4 + aa)
      drevdt   = 6.75d0 * daa

      rr    = den * rev
      drrdt = den * drevdt * 1.0d-9
      drrdd = rev

      return
      end





      subroutine rate_si27pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd


! si27(p,g)p28

      aa  = 1.64e+09 * t9i23 * exp(-24.439*t9i13)
      daa = aa*(-twoth*t9i + oneth*24.439*t9i43)

      bb  = 2.00e-08 * t9i32 * exp(-0.928*t9i)
      dbb = bb*(-1.5d0*t9i + 0.928*t9i2)

      cc  = 1.95e-02 * t9i32 * exp(-1.857*t9i)
      dcc = cc*(-1.5d0*t9i + 1.857*t9i2)

      dd  = 3.70e+02 * t9i47 * exp(-3.817*t9i)
      ddd = dd*(-foursev*t9i + 3.817*t9i2)

      term    = aa + bb + cc + dd
      dtermdt = daa + dbb + dcc + ddd


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.62e+10 * t932 * exp(-23.960*t9i)
      drevdt   = rev*(1.5d0*t9i + 23.960*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_si28pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,xx,dxx,q1
      parameter        (q1 = 1.0d0/8.4681d0)



! si28(p,g)p29

! champagne et al 96

      if (t9.le.5.0) then

       aa  = 8.44d+08 * t9i23 * exp(-24.389*t9i13 - t92*q1)
       daa = aa*(-twoth*t9i + oneth*24.389*t9i43 - 2.0d0*t9*q1)

       bb  = 1.0d0 + 0.17*t913 + 0.113*t923 + 0.0135*t9 &
             + 0.194*t943 + 0.0591*t953
       dbb = oneth*0.17*t9i23 + twoth*0.113*t9i13 + 0.0135 &
            + fourth*0.194*t913 + fiveth*0.0591*t923

       xx  = aa * bb
       dxx = daa*bb + aa*dbb

       cc   = 2.92d+02 * t9i32 * exp(-4.157*t9i)
       dcc  = cc*(-1.5d0*t9i + 4.157*t9i2)

       dd   = 4.30d+05 * t9i32 * exp(-18.51*t9i)
       ddd  = dd*(-1.5d0*t9i + 18.51*t9i2)

       ee   = 6.05d+03 * t9i32 * exp(-18.17*t9i)
       dee  = ee*(-1.5d0*t9i + 18.17*t9i2)

       term    = xx + cc + dd + ee
       dtermdt = dxx + dcc + ddd + dee


! cf88
      else

       aa  = 1.64d+08 * t9i23 * exp(-24.449*t9i13 - t92*q1)
       daa = aa*(-twoth*t9i + oneth*24.449*t9i43 - 2.0d0*t9*q1)

       bb  = 1.0d0 + 0.017*t913 - 4.11*t923 - 0.491*t9 &
             + 5.22*t943 + 1.58*t953
       dbb = oneth*0.017*t9i23 - twoth*4.11*t9i13 - 0.491 &
            + fourth*5.22*t913 + fiveth*1.58*t923

       xx  = aa * bb
       dxx = daa*bb + aa*dbb

       cc   = 3.52d+02 * t9i32 * exp(-4.152*t9i)
       dcc  = cc*(-1.5d0*t9i + 4.152*t9i2)

       dd   = 6.3d+05 * t9i32 * exp(-18.505*t9i)
       ddd  = dd*(-1.5d0*t9i + 18.505*t9i2)

       ee   = 1.69d+03 * exp(-14.518*t9i)
       dee  = ee*14.518*t9i2

       term    = xx + cc + dd + ee
       dtermdt = dxx + dcc + ddd + dee

      end if


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 9.46d+09 * t932 * exp(-31.879*t9i)
      drevdt   = rev*(1.5d0*t9i + 31.879*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_si29pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,xx,dxx,q1
      parameter        (q1 = 1.0d0/0.065536d0)



! si29(p,g)p30

      aa  = 3.26e+09 * t9i23 * exp(-24.459*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*24.459*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.017*t913 + 4.27*t923 + 0.509*t9 &
            + 15.40*t943 + 4.67*t953
      dbb = oneth*0.017*t9i23 + twoth*4.27*t9i13 + 0.509 &
           + fourth*15.40*t913 + fiveth*4.67*t923

      xx  = aa * bb
      dxx = daa*bb + aa*dbb

      cc   = 2.98e+03 * t9i32 * exp(-3.667*t9i)
      dcc  = cc*(-1.5d0*t9i + 3.667*t9i2)

      dd   = 3.94e+04 * t9i32 * exp(-4.665*t9i)
      ddd  = dd*(-1.5d0*t9i + 4.665*t9i2)

      ee   = 2.08e+04 * t912 * exp(-8.657*t9i)
      dee  = ee*(0.5d0*t9i + 8.657*t9i2)

      term    = xx + cc + dd + ee
      dtermdt = dxx + dcc + ddd + dee

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.26e+10 * t932 * exp(-65.002*t9i)
      drevdt   = rev*(1.5d0*t9i + 65.002*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_si30pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,xx,dxx,q1
      parameter        (q1 = 1.0d0/0.4489d0)



! si30(p,g)p31

      aa  = 4.25e8 * t9i23 * exp(-24.468*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*24.468*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.017*t913 + 0.150*t923 + 0.018*t9 &
            + 5.53*t943 + 1.68*t953
      dbb = oneth*0.017*t9i23 + twoth*0.150*t9i13 + 0.018 &
           + fourth*5.53*t913 + fiveth*1.68*t923

      xx  = aa * bb
      dxx = daa*bb + aa*dbb

      cc   = 1.86e4 * t9i32 * exp(-5.601*t9i)
      dcc  = cc*(-1.5d0*t9i + 5.601*t9i2)

      dd   = 3.15e5 * t9i32 * exp(-6.961*t9i)
      ddd  = dd*(-1.5d0*t9i + 6.961*t9i2)

      ee   = 2.75e5 * t9i12 * exp(-10.062*t9i)
      dee  = ee*(-0.5d0*t9i + 10.062*t9i2)

      term    = xx + cc + dd + ee
      dtermdt = dxx + dcc + ddd + dee

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 9.50e9 * t932 * exp(-84.673*t9i)
      drevdt   = rev*(1.5d0*t9i + 84.673*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end








      subroutine rate_weaknp(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision aa,daa,bb,dbb,cc,dcc, &
                       zm1,zm2,zm3,zm4,zm5, &
                       c1,c2
      parameter        (c1 = 1.0d0/5.93d0, &
                        c2 = 0.98d0/885.7d0)


! free decay of neutrons, n(e-nu)p and p(e-,nu)n
! fit formula from schramm and wagoner annual review 1977
! currently accepted best value for the neutron lifetime,
! 886.7 (+/- 1.9) seconds. P.R. Huffman et al., Nature, 6 January 2000.
! world average 885.7 +/- 0.8.

      zm1   = t9 * c1
      zm2   = zm1*zm1
      zm3   = zm1*zm2
      zm4   = zm1*zm3
      zm5   = zm1*zm4

      aa   = 27.512*zm5 + 36.492*zm4 + 11.108*zm3 &
             - 6.382*zm2 + 0.565*zm1 + 1.0d0
      daa  = (5.0d0*27.512*zm4 + 4.0d0*36.492*zm3 + 3.0d0*11.108*zm2 &
             - 2.0d0*6.382*zm1 + 0.565)*c1


! n=>p
      fr    = c2 * aa
      dfrdt = c2 * daa * 1.0d-9
      dfrdd = 0.0d0


      aa  = 27.617*zm5 + 34.181*zm4 + 18.059*zm3 &
            - 16.229*zm2 + 5.252*zm1
      daa = (5.0d0*27.617*zm4 + 4.0d0*34.181*zm3 + 3.0d0*18.059*zm2 &
            - 2.0d0*16.229*zm1 + 5.252)*c1

      bb = exp(-2.531d0/zm1)
      dbb = bb*2.531d0/zm2*c1

      cc  = aa*bb
      dcc = daa*bb + aa*dbb

! p=>n
      rr    = c2 * cc
      drrdt = c2 * dcc * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_dpn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc


! d(p,n)2p
       aa   = 3.35e7 * exp(-3.720*t9i13)
       daa  = aa*oneth*3.720d0*t9i43

       bb   = 1.0d0 + 0.784*t913 + 0.346*t923 + 0.690*t9
       dbb  = oneth*0.784*t9i23 + twoth*0.346*t9i13 + 0.690

       term    = aa * bb
       dtermdt = daa * bb + aa * dbb

! rate
      cc = exp(-25.815*t9i)
      dcc = cc*25.815*t9i2

      fr    = den * cc * term
      dfrdt = den * (dcc*term + cc*dtermdt) * 1.0d-9
      dfrdd = cc * term

      rev      =  4.24e-10 * t9i32
      drevdt   = -1.5d0*rev*t9i

      rr    = den**2 * rev * term
      drrdt = den**2 * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 2.0d0*den * rev * term

      return
      end




      subroutine rate_dng(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,c1
      parameter        (c1 = 66.2d0*18.9d0)


! d(n,g)t
      term    = 66.2 * (1.0d0  + 18.9*t9)
      dtermdt = c1

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      =  1.63e+10 * t9i32 * exp(-72.62*t9i)
      drevdt   = rev*(-1.5d0*t9i + 72.62*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_ddp(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb



! d(d,p)t
      aa  = 4.13e8 * t9i23 * exp(-4.258*t9i13)
      daa = -twoth*aa*t9i + oneth*aa*4.258*t9i43

      bb  = 1.0d0 + 0.098*t913 + 4.39e-2*t923 + 3.01e-2*t9 &
            + 0.543*t943 + 0.946*t953
      dbb = oneth*0.098*t9i23 + twoth*4.39e-2*t9i13 + 3.01e-2 &
            + fourth*0.543*t913 + fiveth*0.946*t923

      term    = aa * bb
      dtermdt = daa*bb + aa*dbb

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.73 * exp(-46.798*t9i)
      drevdt   = rev*46.798*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_ddn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb



! d(d,n)he3
      aa  = 3.88e8 * t9i23 * exp(-4.258*t9i13)
      daa = -twoth*aa*t9i + oneth*aa*4.258*t9i43

      bb  = 1.0d0 + 0.098*t913 + 0.418*t923 + 0.287*t9 &
            + 0.638*t943 + 1.112*t953
      dbb = oneth*0.098*t9i23 + twoth*0.418*t9i13 + 0.287 &
            + fourth*0.638*t913 + fiveth*1.112*t923

      term    = aa * bb
      dtermdt = daa*bb + aa*dbb

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.730 * exp(-37.935*t9i)
      drevdt   = rev*37.935*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end




      subroutine rate_tpn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa


! t(p,n)he3
      term    = 7.07e8 * (1.0d0 - 0.15*t912 + 0.098*t9)
      dtermdt = 7.07e8 * (-0.5d0*0.15*t9i12 + 0.098)

      aa  = exp(-8.863*t9i)
      daa = aa*8.863*t9i2

! rate
      fr    = den * aa * term
      dfrdt = den * (daa*term + aa*dtermdt) * 1.0d-9
      dfrdd = aa * term

      rev      = 0.998
      drevdt   = 0.0d0

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end




      subroutine rate_ddg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb


! d(d,g)he4
      aa  = 4.84e+01 * t9i23 * exp(-4.258*t9i13)
      daa = aa*(-twoth*t9i + oneth*4.258*t9i43)

      bb  = 1.0d0 + 0.098*t913 - 0.203*t923 - 0.139*t9 &
            + 0.106*t943 + 0.185*t953
      dbb = oneth*0.098*t9i23 - twoth*0.203*t9i13 - 0.139 &
            + fourth*0.106*t913 + fiveth*0.185*t923

      term    = aa * bb
      dtermdt = daa*bb + aa*dbb

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 4.53e+10 * t932 * exp(-276.729*t9i)
      drevdt   = rev*(1.5d0*t9i + 276.729*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_tpg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb


! t(p,g)he4
      aa  = 2.20e+04 * t9i23 * exp(-3.869*t9i13)
      daa = aa*(-twoth*t9i + oneth*3.869*t9i43)

      bb  = 1. + 0.108*t913 + 1.68*t923 + 1.26*t9 &
            + 0.551*t943 + 1.06*t953
      dbb = oneth*0.108*t9i23 + twoth*1.68*t9i13 + 1.26 &
            + fourth*0.551*t913 + fiveth*1.06*t923

      term    = aa * bb
      dtermdt = daa*bb + aa*dbb

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 2.61e+10 * t932 * exp(-229.932*t9i)
      drevdt   = rev*(1.5d0*t9i + 229.932*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_tdn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc,q1
      parameter        (q1 = 1.0d0/0.0144d0)


! t(d,n)he4 ; the "dt" reaction
      aa  = 8.09e+10 * t9i23 * exp(-4.524*t9i13 - t92*q1)
      daa = -twoth*aa*t9i + aa*(oneth*4.524*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.092*t913 + 1.80*t923 + 1.16*t9 &
            + 10.52*t943 + 17.24*t953
      dbb = oneth*0.092*t9i23 + twoth*1.80*t9i13 + 1.16 &
            + fourth*10.52*t913 + fiveth*17.24*t923

      cc  = 8.73e+08 * t9i23 * exp(-0.523*t9i)
      dcc = -twoth*cc*t9i + cc*0.523*t9i2

      term    = aa * bb + cc
      dtermdt = daa*bb + aa*dbb + dcc

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 5.54*exp(-204.117*t9i)
      drevdt   = rev*204.117*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_tt2n(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb


! t(t,2n)he4
      aa  = 1.67e+09 * t9i23 * exp(-4.872*t9i13)
      daa = aa*(-twoth*t9i + oneth*4.872*t9i43)

      bb  = 1.0d0 + 0.086*t913 - 0.455*t923 - 0.272*t9 &
            + 0.148*t943 + 0.225*t953
      dbb = oneth*0.086*t9i23 - twoth*0.455*t9i13 - 0.272 &
            + fourth*0.148*t913 + fiveth*0.225*t923

      term    = aa * bb
      dtermdt = daa*bb + aa*dbb

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.38e-10 * t9i32 * exp(-131.504*t9i)
      drevdt   = rev*(-1.5d0*t9i + 131.504*t9i2)

      rr    = den * den * rev * term
      drrdt = den * den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 2.0d0 * den * rev * term

      return
      end





      subroutine rate_he3dp(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc,q1
      parameter        (q1 = 1.0d0/0.099225d0)


! he3(d,p)he4
      aa  = 5.86e+10 * t9i23 * exp(-7.181*t9i13 - t92*q1)
      daa = -twoth*aa*t9i + aa*(oneth*7.181*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.058*t913 + 0.142*t923 + 0.0578*t9 &
            + 2.25*t943 + 2.32*t953
      dbb = oneth*0.058*t9i23 + twoth*0.142*t9i13 + 0.0578 &
            + fourth*2.25*t913 + fiveth*2.32*t923

      cc  = 4.36e+08 * t9i12 * exp(-1.72*t9i)
      dcc = -0.5d0*cc*t9i + cc*1.72*t9i2

      term    = aa * bb + cc
      dtermdt = daa*bb + aa*dbb + dcc

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 5.55*exp(-212.980*t9i)
      drevdt   = rev*212.980*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end




      subroutine rate_he3td(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,t9a,dt9a, &
                       t9a13,dt9a13,t9a56,dt9a56,zz


! he3(t,d)he4
      aa       = 1.0d0 + 0.128*t9
      zz       = 1.0d0/aa

      t9a      = t9*zz
      dt9a     = (1.0d0 - t9a*0.128)*zz

      zz      = dt9a/t9a
      t9a13   = t9a**oneth
      dt9a13  = oneth*t9a13*zz

      t9a56    = t9a**fivsix
      dt9a56   = fivsix*t9a56*zz

      term     = 5.46e+09 * t9a56 * t9i32 * exp(-7.733/t9a13)
      dtermdt = term*(dt9a56/t9a56 - 1.5d0*t9i &
                      + 7.733/t9a13**2 * dt9a13)

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.60*exp(-166.182*t9i)
      drevdt   = rev*166.182*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end




      subroutine rate_he3tnp(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,t9a,dt9a, &
                       t9a13,dt9a13,t9a56,dt9a56,zz


! he3(t,np)he4
      aa       = 1.0d0 + 0.115d0 * t9
      zz       = 1.0d0/aa

      t9a      = t9*zz
      dt9a     = (1.0d0 - t9a * 0.115d0)*zz

      zz      = dt9a/t9a
      t9a13   = t9a**oneth
      dt9a13  = oneth*t9a13*zz

      t9a56    = t9a**fivsix
      dt9a56   = fivsix*t9a56*zz

      term     = 7.71d9 * t9a56 * t9i32 * exp(-7.733d0/t9a13)
      dtermdt = term*(dt9a56/t9a56 - 1.5d0*t9i &
                      + 7.733d0/t9a13**2 * dt9a13)

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.39d-10 * t9i32 * exp(-140.367*t9i)
      drevdt   = rev*(-1.5d0 * t9i + 140.367d0 * t9i2)

      rr    = den * den * rev * term
      drrdt = den * den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 2.0d0 * den * rev * term

      return
      end




      subroutine rate_he4npg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb


! he4(np,g)li6
      aa  = 4.62e-6 * t9i2 * exp(-19.353*t9i)
      daa = aa*(-2.0d0*t9i + 19.353*t9i2)

      bb  = 1.0d0 + 0.075*t9
      dbb = 0.075

      term    = aa * bb
      dtermdt = daa*bb + aa*dbb

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 7.22e19 * t93 * exp(-42.933*t9i)
      drevdt   = rev*(3.0d0*t9i + 42.933*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_he4dg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc


! he4(d,g)li6
      aa  = 3.01e1 * t9i23 * exp(-7.423*t9i13)
      daa = aa*(-twoth*t9i + oneth*7.423*t9i43)

      bb  = 1.0d0 + 0.056*t913 - 4.85*t923 + 8.85*t9 &
            - 0.585*t943 - 0.584*t953
      dbb = oneth*0.056*t9i23 - twoth*4.85*t9i13 + 8.850 &
            - fourth*0.585*t913 - fiveth*0.584*t923

! rate goes negative for t9 greater than about 15, so try this
      if (bb .le. 0.0) then
       bb = 0.0d0
       dbb = 0.0d0
      end if

      cc =  8.55e1 * t9i32 * exp(-8.228*t9i)
      dcc = cc*(-1.5d0*t9i + 8.228*t9i2)

      term    = aa * bb + cc
      dtermdt = daa*bb + aa*dbb + dcc


! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.53e10 * t932 * exp(-17.1180*t9i)
      drevdt   = rev*(1.5d0*t9i + 17.1180*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_he4tn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,t9a,dt9a,t9a32,dt9a32,zz


! he4(t,n)li6
      aa     = 1.0d0 + 49.180*t9
      zz     = 1.0d0/aa

      t9a   = t9*zz
      dt9a   = (1.0d0 - t9a*49.180)*zz

      t9a32  = t9a * sqrt(t9a)
      dt9a32 = 1.5d0*t9a32/t9a * dt9a

      aa     = 1.80e8 * exp(-55.4940*t9i)
      daa    = aa*55.4940*t9i2

      bb     = 1.0d0 - 0.2610 * t9a32 * t9i32
      dbb    = -0.2610*(-1.5d0*t9a32*t9i52 + dt9a32*t9i32)

      cc     = aa*bb
      dcc    = daa*bb + aa*dbb

      dd     = 2.72e9 * t9i32 * exp(-57.8840*t9i)
      ddd    = dd*(-1.5d0*t9i + 57.8840*t9i2)

      term    = cc + dd
      dtermdt = dcc + ddd

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 2.72e9 * t9i32 * exp(-2.39*t9i)
      drevdt   = rev*(-1.5d0*t9i + 2.39*t9i2)

      term = 0.935*(1.80e8*bb + rev)
      dtermdt = 0.935*(1.80e8*dbb + drevdt)

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end




      subroutine rate_li6phe3(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,q1
      parameter        (q1 = 1.0d0/30.25d0)


! li6(p,he3)he4
      aa  = 3.73e10 * t9i23 * exp(-8.413*t9i13 - t92*q1)
      daa = -twoth*aa*t9i + aa*(oneth*8.413*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.050*t913 - 0.061*t923 - 0.0210*t9 &
            + 0.0060*t943 + 0.0050*t953
      dbb = oneth*0.050*t9i23 - twoth*0.061*t9i13 - 0.0210 &
            + fourth*0.0060*t913 + fiveth*0.0050*t923

      cc  = 1.33e10 * t9i32 * exp(-17.7630*t9i)
      dcc = -1.5d0*cc*t9i + cc*17.7630*t9i2

      dd  =  1.29e9 * t9i * exp(-21.82*t9i)
      ddd = -dd*t9i + dd*21.82*t9i2

      term    = aa * bb + cc + dd
      dtermdt = daa*bb + aa*dbb + dcc + ddd

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.070d0 * exp(-46.6310*t9i)
      drevdt   = rev*46.6310*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_li6ng(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt


! li6(n,g)li7
! malaney-fowler 1989

      term    = 5.10e3
      dtermdt = 0.0d0

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.19e10 * t932 * exp(-84.17*t9i)
      drevdt   = rev*(1.5d0*t9i + 84.17*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_he4tg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb


! he4(t,g)li7
      aa  = 8.67e5 * t9i23 * exp(-8.08*t9i13)
      daa = aa*(-twoth*t9i + oneth*8.08*t9i43)

      bb  = 1.0d0 + 0.052*t913 - 0.448*t923 - 0.165*t9 &
            + 0.144*t943 + 0.134*t953
      dbb = oneth*0.052*t9i23 - twoth*0.448*t9i13 - 0.165 &
            + fourth*0.144*t913 + fiveth*0.134*t923

      term    = aa * bb
      dtermdt = daa*bb + aa*dbb

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.11e10 * t932 * exp(-28.64*t9i)
      drevdt   = rev*(1.5d0*t9i + 28.64*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end



      subroutine rate_li7dn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt


! li7(d,n)2a
      term    = 2.92e11 * t9i23 * exp(-10.259*t9i13)
      dtermdt = term*(-twoth*t9i + oneth*10.259*t9i43)

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 9.95e-10 * t9i32 * exp(-175.476*t9i)
      drevdt   = rev*(-1.5d0*t9i + 175.476*t9i2)

      rr    = den * den * rev * term
      drrdt = den * den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 2.0d0 * den * rev * term

      return
      end




      subroutine rate_li7tn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt


! li7(t,n)be9
! malaney and fowler (apjl, 345, l5, 1989)
      term    = 1.46d+11 * t9i23 * exp(-11.333*t9i13)
      dtermdt = -twoth*term*t9i + oneth*term*11.333*t9i43

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 0.0d0
      drevdt   = 0.0d0

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_li7t2n(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt


! li7(t,2n)2a
      term    = 8.81e11 * t9i23 * exp(-11.333*t9i13)
      dtermdt = term*(-twoth*t9i + oneth*11.333*t9i43)

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.22e-19 * t9i3 * exp(-102.864*t9i)
      drevdt   = rev*(-3.0d0*t9i + 102.864*t9i2)

      rr    = den**3 * rev * term
      drrdt = den**3 * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 3.0d0 * den**2 * rev * term

      return
      end




      subroutine rate_li7he3np(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt


! li7(he3,np)2a
      term    = 1.11e13 * t9i23 * exp(-17.989*t9i13)
      dtermdt = term*(-twoth*t9i + oneth*17.989*t9i43)

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 6.09e-20 * t9i3 * exp(-111.727*t9i)
      drevdt   = rev*(-3.0d0*t9i + 111.727*t9i2)

      rr    = den**3 * rev * term
      drrdt = den**3 * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 3.0d0 * den**2 * rev * term

      return
      end




      subroutine rate_li6pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,bb,dbb,cc,dcc, &
                       t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56,zz


! li6(p,g)be7
      if (t9 .gt. 10.0) then
       t9a  = 1.0d0
       dt9a = 0.0d0
      else
       aa   = 1.0d0 - 0.0969*t9

       bb   = aa**(-twoth)
       dbb  = twoth*bb/aa*0.0969

       cc   = aa + 0.0284*t953*bb
       dcc  = -0.0969 + 0.0284*(fiveth*t923*bb + t953*dbb)

       zz   = 1.0d0/cc
       t9a  = t9*zz
       dt9a = (1.0d0 - t9a*dcc)*zz
      end if

      zz     = dt9a/t9a
      t9a13  = t9a**oneth
      dt9a13 = oneth*t9a13*zz

      t9a56  = t9a**fivsix
      dt9a56 = fivsix*t9a56*zz

      term    = 6.69e+05 * t9a56 * t9i32 * exp(-8.413/t9a13)
      dtermdt = term*(dt9a56/t9a56 - 1.5d0*t9i &
                      + 8.413/t9a13**2 * dt9a13)

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.19e+10 * t932 * exp(-65.054*t9i)
      drevdt   = rev*(1.5d0*t9i + 65.054*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0


      return
      end




      subroutine rate_li7pn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,bb,dbb


! li7(p,n)be7
      aa  = 5.15e+09 * exp(-1.167*t913 - 19.081*t9i)
      daa = aa*(-oneth*1.167*t9i23 + 19.081*t9i2)

      bb  = 7.84e+09 * t9i32 * exp(-22.832*t9i)
      dbb = -1.5d0*bb*t9i + bb*22.832*t9i2

      term    = aa + bb
      dtermdt = daa + dbb

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      aa  = 5.15e+09 * exp(-1.167*t913)
      daa = -aa*oneth*1.167*t9i23

      bb  = 0.998 * 7.84e+09 * t9i32 * exp(-3.751*t9i)
      dbb = -1.5d0*bb*t9i + bb*3.751*t9i2

      term    = aa + bb
      dtermdt = daa + dbb
      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end




      subroutine rate_li7ng(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa


! li7(n,g)li8
! apj 372, 1
      aa      = 4.26d+03 * t9i32 * exp(-2.576*t9i)
      daa     = aa*(-1.5d0*t9i + 2.576*t9i2)

      term    = 3.144d+03 + aa
      dtermdt = daa

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev    = 1.2923d+10 * t932 * exp(-2.359d+01*t9i)
      drevdt = rev*(1.5d0*t9i + 2.359d+01*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_be7dp(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt


! be7(d,p)2a
      term    = 1.07e12 * t9i23 * exp(-12.428*t9i13)
      dtermdt = term*(-twoth*t9i + oneth*12.428*t9i43)

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 9.97e-10 * t9i32 * exp(-194.557*t9i)
      drevdt   = rev*(-1.5d0*t9i + 194.557*t9i2)

      rr    = den * den * rev * term
      drrdt = den * den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 2.0d0 * den * rev * term

      return
      end




      subroutine rate_be7tnp(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt


! be7(t,np)2a
      term    = 2.91e12 * t9i23 * exp(-13.729*t9i13)
      dtermdt = term*(-twoth*t9i + oneth*13.729*t9i43)

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 6.09e-20 * t9i3 * exp(-121.944*t9i)
      drevdt   = rev*(-3.0d0*t9i + 121.944*t9i2)

      rr    = den**3 * rev * term
      drrdt = den**3 * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 3.0d0 * den**2 * rev * term

      return
      end



      subroutine rate_be7he32p(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt


! be7(he3,2p)2a
      term    = 6.11e13 * t9i23 * exp(-21.793*t9i13)
      dtermdt = term*(-twoth*t9i + oneth*21.793*t9i43)

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.22e-19 * t9i3 * exp(-130.807*t9i)
      drevdt   = rev*(-3.0d0*t9i + 130.807*t9i2)

      rr    = den**3 * rev * term
      drrdt = den**3 * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 3.0d0 * den**2 * rev * term

      return
      end




      subroutine rate_be9pa(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,q1
      parameter        (q1 = 1.0d0/0.2704d0)


! be9(p,a)li6
      aa  = 2.11e11 * t9i23 * exp(-10.359*t9i13 - t92*q1)
      daa = -twoth*aa*t9i + aa*(oneth*10.359*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0 + 0.04*t913 + 1.09*t923 + 0.307*t9 &
            + 3.21*t943 + 2.30*t953
      dbb  = oneth*0.04*t9i23 + twoth*1.09*t9i13 + 0.307 &
             + fourth*3.21*t913 + fiveth*2.30*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 4.51e8 * t9i * exp(-3.046*t9i)
      ddd  = -dd*t9i + dd*3.046*t9i2

      ee   = 6.70e8 * t9i34 * exp(-5.160*t9i)
      dee  = -0.75d0*ee*t9i + ee*5.160*t9i2

      term    = cc + dd + ee
      dtermdt = dcc + ddd + dee

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 6.18e-1 * exp(-24.674*t9i)
      drevdt   = rev*24.674*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end




      subroutine rate_li6ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,q1
      parameter        (q1 = 1.0d0/1.758276d0)


! li6(a,g)b10
      aa  = 4.06e6 * t9i23 * exp(-18.790*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*18.790*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0 + 0.022*t913 + 1.54*t923 + 0.239*t9 &
             +  2.20*t943 + 0.869*t953
      dbb  = oneth*0.022*t9i23 + twoth*1.54*t9i13 + 0.239 &
             + fourth*2.20*t913 + fiveth*0.869*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 1910.0 * t9i32 * exp(-3.484*t9i)
      ddd  = dd*(-1.5d0*t9i + 3.484*t9i2)

      ee   = 1.01e4 * t9i * exp(-7.269*t9i)
      dee  = ee*(-t9i + 7.269*t9i2)

      term    = cc + dd + ee
      dtermdt = dcc + ddd + dee

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.58e10 * t932 * exp(-51.753*t9i)
      drevdt   = rev*(1.5d0*t9i + 51.753*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_li7an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt

! li7(a,n)b10
      term    = 3.84e8 * exp(-32.382*t9i)
      dtermdt = term*32.382*t9i2

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.84e8*1.32
      drevdt   = 0.0d0

      rr    = den * rev
      drrdt = 0.0d0
      drrdd = rev

      return
      end





      subroutine rate_be9pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,q1
      parameter        (q1 = 1.0d0/0.715716d0)


! be9(p,g)b10
      aa  = 1.33e7 * t9i23 * exp(-10.359*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*10.359*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0 + 0.04*t913 + 1.52*t923 + 0.428*t9 &
             + 2.15*t943 + 1.54*t953
      dbb  = oneth*0.04*t9i23 + twoth*1.52*t9i13 + 0.428 &
             + fourth*2.15*t913 + fiveth*1.54*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 9.64e4 * t9i32 * exp(-3.445*t9i)
      ddd  = dd*(-1.5d0*t9i + 3.445*t9i2)

      ee   = 2.72e6 * t9i32 * exp(-10.62*t9i)
      dee  = ee*(-1.5d0*t9i + 10.62*t9i2)

      term    = cc + dd + ee
      dtermdt = dcc + ddd + dee

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 9.73e9 * t932 * exp(-76.427*t9i)
      drevdt   = rev*(1.5d0*t9i + 76.427*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_b10pa(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,q1
      parameter        (q1 = 1.0d0/19.377604d0)

! b10(p,a)li7
      aa  = 1.26e11 * t9i23 * exp(-12.062*t9i13 - t92*q1)
      daa = -twoth*aa*t9i + aa*(oneth*12.062*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0 + 0.035*t913 - 0.498*t923 - 0.121*t9 &
             + 0.3*t943 + 0.184*t953
      dbb  = oneth*0.035*t9i23 - twoth*0.498*t9i13 - 0.121 &
             + fourth*0.3*t913 + fiveth*0.184*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 2.59e9 * t9i * exp(-12.260*t9i)
      ddd  = -dd*t9i + dd*12.260*t9i2

      term    = cc + dd
      dtermdt = dcc + ddd

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 7.54e-01 * exp(-13.301*t9i)
      drevdt   = rev*13.301*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end




      subroutine rate_li7ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,q1
      parameter        (q1 = 1.0d0/17.598025d0)


! li7(a,g)b11
      aa  = 3.55e7 * t9i23 * exp(-19.161*t9i13 -t92*q1)
      daa = aa*(-twoth*t9i + oneth*19.161*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0 + 0.022*t913 + 0.775*t923 + 0.118*t9 &
             + 0.884*t943 + 0.342*t953
      dbb  = oneth*0.022*t9i23 + twoth*0.775*t9i13 + 0.118 &
             + fourth*0.884*t913 + fiveth*0.342*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 3.33e2 * t9i32 * exp(-2.977*t9i)
      ddd  = dd*(-1.5d0*t9i + 2.977*t9i2)

      ee   = 4.10e4 * t9i * exp(-6.227*t9i)
      dee  = ee*(-t9i + 6.227*t9i2)

      term    = cc + dd + ee
      dtermdt = dcc + ddd + dee

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 4.02e10 * t932 * exp(-100.538*t9i)
      drevdt   = rev*(1.5d0*t9i + 100.538*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_b11pa(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,ff,dff,q1
      parameter        (q1 = 1.0d0/2.702736d0)

! b11(p,a)be8=>2a
      aa  = 2.20e12 * t9i23 * exp(-12.095*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*12.095*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0  + 0.034*t913 + 0.14*t923 + 0.034*t9 &
             + 0.19*t943 + 0.116*t953
      dbb  = oneth*0.034*t9i23 + twoth*0.14*t9i13 + 0.034 &
             + fourth*0.19*t913 + fiveth*0.116*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb


      dd   = 4.03e6 * t9i32 * exp(-1.734*t9i)
      ddd  = dd*(-1.5d0*t9i + 1.734*t9i2)

      ee   = 6.73e9 * t9i32 * exp(-6.262*t9i)
      dee  = ee*(-1.5d0*t9i + 6.262*t9i2)

      ff   = 3.88e9*t9i * exp(-14.154*t9i)
      dff  = ff*(-t9i + 14.154*t9i2)

      term    = cc + dd + ee + ff
      dtermdt = dcc + ddd +dee + dff

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.5e-10* t9i32 *exp(-100.753*t9i)
      drevdt   = rev*(-1.5d0*t9i + 100.753*t9i2)

      rr    = den * den * rev * term 
      drrdt = den * den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 2.0d0 * den * rev * term 

      return
      end




      subroutine rate_be7ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,q1
      parameter        (q1 = 1.0d0/22.743361d0)

! be7(a,g)c11
      aa  = 8.45e+07 * t9i23 * exp(-23.212*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*23.212*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0 + 0.018*t913 + 0.488*t923 + 0.061*t9 &
             + 0.296*t943 + 0.095*t953
      dbb  = oneth*0.018*t9i23 + twoth*0.488*t9i13 + 0.061 &
             + fourth*0.296*t913 + fiveth*0.095*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 1.25e+04 * t9i32 * exp(-6.510*t9i)
      ddd  = dd*(-1.5d0*t9i + 6.510*t9i2)

      ee   = 1.29e+05 * t9i54 * exp(-10.039*t9i)
      dee  = ee*(-1.25d0*t9i + 10.039*t9i2)

      term    = cc + dd + ee
      dtermdt = dcc + ddd + dee

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 4.02e+10 * t932 * exp(-87.539*t9i)
      drevdt   = rev*(1.5d0*t9i + 87.539*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end



      subroutine rate_b11pn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,bb,dbb


! b11(p,n)c11
      aa  = 1.69e8*(1.0d0 - 0.048*t912 + 0.010*t9)
      daa = 1.69e8*(-0.5d0*0.048*t9i12 + 0.010)

      bb  = exp(-32.080*t9i)
      dbb = bb*32.080*t9i2

      term    = aa*bb
      dtermdt = daa*bb + aa*dbb

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      term     = 0.998*aa
      dtermdt  = 0.998*daa

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end





      subroutine rate_b8ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,ff,dff,gg,dgg,hh,dhh, &
                       pp,dpp,qq,dqq,drr


! b8(a,p)c11
      aa  = 1.67e-09 * t9i32 * exp(-1.079*t9i)
      daa = -1.5d0*aa*t9i + aa*1.079*t9i2

      bb  = 9.55e-08 * t9i32 * exp(-1.311*t9i)
      dbb = -1.5d0*bb*t9i + bb*1.311*t9i2

      cc  = 1.98e-01 * t9i32 * exp(-2.704*t9i)
      dcc = -1.5d0*cc*t9i + cc*2.704*t9i2

      dd  = 1.34e+00 * t9i32 * exp(-4.282*t9i)
      ddd = -1.5d0*dd*t9i + dd*4.282*t9i2

      ee  = 3.22e+04 * t9i32 * exp(-6.650*t9i)
      dee = -1.5d0*ee*t9i + ee*6.650*t9i2

      ff  = 2.33e+05 * t9i32 * exp(-8.123*t9i)
      dff = -1.5d0*ff*t9i + ff*8.123*t9i2

      gg  = 2.55e+06 * t9i32 * exp(-11.99*t9i)
      dgg = -1.5d0*gg*t9i + gg*11.99*t9i2

      hh  = 9.90e+06 * t9i32 * exp(-13.50*t9i)
      dhh = -1.5d0*hh*t9i + hh*13.50*t9i2

      pp  = 1.41e+06 * t9i32 * exp(-16.51*t9i)
      dpp = -1.5d0*pp*t9i + pp*16.51*t9i2

      qq  = 1.99e+07 * t9i32 * exp(-18.31*t9i)
      dqq = -1.5d0*qq*t9i + qq*18.31*t9i2

      rr  = 6.01e+07 * t9i32 * exp(-20.63*t9i)
      drr = -1.5d0*rr*t9i + rr*20.63*t9i2

      term    = aa + bb + cc + dd + ee + ff + gg + hh + pp + qq + rr
      dtermdt = daa + dbb + dcc + ddd + dee + dff + dgg + dhh &
                + dpp + dqq + drr

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.101 * exp(-85.95*t9i)
      drevdt   = rev*85.95*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev*term

      return
      end





      subroutine rate_b10pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,q1
      parameter        (q1 = 1.0d0/19.377604d0)


! b10(p,g)c11
      aa  = 4.61e+05 * t9i23 * exp(-12.062*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*12.062*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0 + 0.035*t913 + 0.426*t923 + 0.103*t9 &
             + 0.281*t943 + 0.173*t953
      dbb  = oneth*0.035*t9i23 + twoth*0.426*t9i13 + 0.103 &
             + fourth*0.281*t913 + fiveth*0.173*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 1.93e+05 * t9i32 * exp(-12.041*t9i)
      ddd  = dd*(-1.5d0*t9i + 12.041*t9i2)

      ee   = 1.14e+04 * t9i32 * exp(-16.164*t9i)
      dee  = ee*(-1.5d0*t9i + 16.164*t9i2)

      term    = cc + dd + ee
      dtermdt = dcc + ddd + dee

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.03e+10 * t932 * exp(-100.840*t9i)
      drevdt   = rev*(1.5d0*t9i + 100.840*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_c11na(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd


! c11(n,a)be8=>2a   get a real rate in here

!      fr    = den * 7.0e4
!      dfrdt = 0.0d0
!      dfrdd = 7.0e4

      fr    = 0.0d0
      dfrdt = 0.0d0
      dfrdd = 0.0d0

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end





      subroutine rate_be9an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,gg,dgg,hh,dhh


! be9(a,n)c12
! Wrean 94 Phys Rev C (1994) vol 49, #2, 1205
      aa  = 6.476d+13 * t9i23 * exp(-23.8702*t9i13)
      daa = -twoth*aa*t9i + oneth*aa*23.8702*t9i43

      bb  = (1.0d0 - 0.3270*t913)
      dbb = -oneth*0.3270*t9i23

! rate goes negative for t9 greater than about 15, so try this
      if (bb .le. 0.0) then
       bb  = 0.0d0
       dbb = 0.0d0
      end if

      cc  = aa*bb
      dcc = daa*bb + aa*dbb

      dd  = 6.044d-3*t9i32*exp(-1.401*t9i)
      ddd = -1.5d0*dd*t9i + dd*1.401*t9i2

      ee  = 7.268*t9i32*exp(-2.063*t9i)
      dee = -1.5d0*ee*t9i + ee*2.063*t9i2

      ff  = 3.256d+4*t9i32*exp(-3.873*t9i)
      dff = -1.5d0*ff*t9i + ff*3.873*t9i2

      gg  = 1.946d+5*t9i32*exp(-4.966*t9i)
      dgg = -1.5d0*gg*t9i + gg*4.966*t9i2

      hh  = 1.838e9*t9i32*exp(-15.39*t9i)
      dhh = -1.5d0*hh*t9i + hh*15.39*t9i2

      term    = cc + dd + ee + ff + gg + hh
      dtermdt = dcc + ddd + dee + dff + dgg + dhh


! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 10.3d0 * exp(-66.160*t9i)
      drevdt   = rev*66.160d0*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term


      return
      end




      subroutine rate_b11pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,q1
      parameter        (q1 = 1.0d0/0.57121d0)


! b11(p,g)c12
      aa  = 4.62e+07 * t9i23 * exp(-12.095*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*12.095*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0 + 0.035*t913 + 3.0*t923 + 0.723*t9 &
             + 9.91*t943 + 6.07*t953
      dbb  = oneth*0.035*t9i23 + twoth*3.0*t9i13 + 0.723 &
             + fourth*9.91*t913 + fiveth*6.07*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 7.89e+03 * t9i32 * exp(-1.733*t9i)
      ddd  = dd*(-1.5d0*t9i + 1.733*t9i2)

      ee   = 9.68e+04 * t9i15 * exp(-5.617*t9i)
      dee  = ee*(-0.2d0*t9i + 5.617*t9i2)

      term    = cc + dd + ee
      dtermdt = dcc + ddd + dee

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 7.01e+10 * t932 * exp(-185.173*t9i)
      drevdt   = rev*(1.5d0*t9i + 185.173*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_b11ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,ff,dff,q1
      parameter        (q1 = 1.0d0/0.120409d0)


! b11(a,p)c14
      aa  = 5.37e+11 * t9i23 * exp(-28.234*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*28.234*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0 + 0.015*t913 + 5.575*t923 + 0.576*t9 &
             + 15.888*t943 + 4.174*t953
      dbb  = oneth*0.015*t9i23 + twoth*5.575*t9i13 + 0.576 &
             + fourth*15.888*t913 + fiveth*4.174*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 5.44e-03 * t9i32 * exp(-2.827*t9i)
      ddd  = dd*(-1.5d0*t9i + 2.827*t9i2)

      ee   = 3.36e+02 * t9i32 * exp(-5.178*t9i)
      dee  = ee*(-1.5d0*t9i + 5.178*t9i2)

      ff   = 5.32e+06 * t9i38 * exp(-11.617*t9i)
      dff  = ff*(-0.375*t9i + 11.617*t9i2)

      term    = cc + dd + ee + ff
      dtermdt = dcc + ddd + dee + dff

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.10e+01 * exp(-9.098*t9i)
      drevdt   = rev*9.098*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end






      subroutine rate_pp(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,bb,dbb


! p(p,e+nu)d
      if (t9 .le. 3.0) then
       aa   = 4.01d-15 * t9i23 * exp(-3.380d0*t9i13)
       daa  = aa*(-twoth*t9i + oneth*3.380d0*t9i43)

       bb   = 1.0d0 + 0.123d0*t913 + 1.09d0*t923 + 0.938d0*t9
       dbb  = oneth*0.123d0*t9i23 + twoth*1.09d0*t9i13 + 0.938d0

       term    = aa * bb
       dtermdt = daa * bb + aa * dbb

      else
       term    = 1.1581136d-15
       dtermdt = 0.0d0
      end if

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end






      subroutine rate_pep(temp,den,ye,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,ye,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,bb,dbb


! p(e-p,nu)d
      if (t9 .le. 3.0) then

       aa   = 1.36e-20 * t9i76 * exp(-3.380*t9i13)
       daa  = aa*(-sevsix*t9i + oneth*3.380d0*t9i43)

       bb   = (1.0d0 - 0.729d0*t913 + 9.82d0*t923)
       dbb  = -oneth*0.729d0*t9i23 + twoth*9.82d0*t9i13

       term    = aa * bb
       dtermdt = daa * bb + aa * dbb

      else
       term    = 7.3824387e-21
       dtermdt = 0.0d0
      end if

! rate
      fr    = ye * den * den * term
      dfrdt = ye * den * den * dtermdt * 1.0d-9
      dfrdd = ye * 2.0d0 * den * term

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end





      subroutine rate_hep(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,bb,dbb


! he3(p,e+nu)he4
      if (t9 .le. 3.0) then

       aa   = 8.78e-13 * t9i23 * exp(-6.141d0*t9i13)
       daa  = aa*(-twoth*t9i + oneth*6.141d0*t9i43)

       term    = aa
       dtermdt = daa

      else
       term    = 5.9733434e-15
       dtermdt = 0.0d0
      end if

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end







      subroutine rate_png(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa


! p(n,g)d
! smith,kawano,malany 1992

       aa      = 1.0d0 - 0.8504*t912 + 0.4895*t9 &
                 - 0.09623*t932 + 8.471e-3*t92 &
                 - 2.80e-4*t952

       daa     =  -0.5d0*0.8504*t9i12 + 0.4895 &
                 - 1.5d0*0.09623*t912 + 2.0d0*8.471e-3*t9 &
                 - 2.5d0*2.80e-4*t932

       term    = 4.742e4 * aa
       dtermdt = 4.742e4 * daa


! wagoner,schramm 1977
!      aa      = 1.0d0 - 0.86*t912 + 0.429*t9
!      daa     =  -0.5d0*0.86*t9i12 + 0.429

!      term    = 4.4d4 * aa
!      dtermdt = 4.4d4 * daa



! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 4.71d+09 * t932 * exp(-25.82*t9i)
      drevdt   = rev*(1.5d0*t9i + 25.82*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_dpg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb


! d(p,g)he3
      aa      = 2.24d+03 * t9i23 * exp(-3.720*t9i13)
      daa     = aa*(-twoth*t9i + oneth*3.720*t9i43)

      bb      = 1.0d0 + 0.112*t913 + 3.38*t923 + 2.65*t9
      dbb     = oneth*0.112*t9i23 + twoth*3.38*t9i13 + 2.65

      term    = aa * bb
      dtermdt = daa * bb + aa * dbb


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.63d+10 * t932 * exp(-63.750*t9i)
      drevdt   = rev*(1.5d0*t9i + 63.750*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_he3ng(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt


! he3(n,g)he4
      term    = 6.62 * (1.0d0 + 905.0*t9)
      dtermdt = 5.9911d3

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 2.61d+10 * t932 * exp(-238.81*t9i)
      drevdt   = rev*(1.5d0*t9i + 238.81*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_he3he3(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb


! he3(he3,2p)he4
      aa   = 6.04d+10 * t9i23 * exp(-12.276*t9i13)
      daa  = aa*(-twoth*t9i + oneth*12.276*t9i43)

      bb   = 1.0d0 + 0.034*t913 - 0.522*t923 - 0.124*t9 &
             + 0.353*t943 + 0.213*t953
      dbb  = oneth*0.034*t9i23 - twoth*0.522*t9i13 - 0.124 &
             + fourth*0.353*t913 + fiveth*0.213*t923

      term    = aa * bb
      dtermdt = daa*bb + aa*dbb

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.39e-10 * t9i32 * exp(-149.230*t9i)
      drevdt   = rev*(-1.5d0*t9i + 149.230*t9i2)

      rr    = den * den * rev * term
      drrdt = den * den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 2.0d0 * den * rev * term

      return
      end





      subroutine rate_he3he4(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,t9a,dt9a, &
                       t9a13,dt9a13,t9a56,dt9a56,zz


! he3(he4,g)be7
      aa      = 1.0d0 + 0.0495*t9
      daa     = 0.0495

      zz      = 1.0d0/aa
      t9a     = t9*zz
      dt9a    = (1.0d0 - t9a*daa)*zz

      zz      = dt9a/t9a
      t9a13   = t9a**oneth
      dt9a13  = oneth*t9a13*zz

      t9a56   = t9a**fivsix
      dt9a56  = fivsix*t9a56*zz

      term    = 5.61d+6 * t9a56 * t9i32 * exp(-12.826/t9a13)
      dtermdt = term*(dt9a56/t9a56 - 1.5d0*t9i &
                      + 12.826/t9a13**2 * dt9a13)

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.11e+10 * t932 * exp(-18.423*t9i)
      drevdt   = rev*(1.5d0*t9i + 18.423*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_be7em(temp,den,ye,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,ye,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,bb,dbb

! be7(e-,nu+g)li7
      if (t9 .le. 3.0) then
       aa  = 0.0027 * t9i * exp(2.515e-3*t9i)
       daa = -aa*t9i - aa*2.515e-3*t9i2

       bb  = 1.0d0 - 0.537*t913 + 3.86*t923 + aa
       dbb = -oneth*0.537*t9i23 + twoth*3.86*t9i13 + daa

       term    = 1.34e-10 * t9i12 * bb
       dtermdt = -0.5d0*term*t9i + 1.34e-10*t9i12*dbb

      else
       term    = 0.0d0
       dtermdt = 0.0d0
      endif

! rates
      fr    = ye * den * term
      dfrdt = ye * den * dtermdt * 1.0d-9
      dfrdd = ye * term

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end




      subroutine rate_be7pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb


! be7(p,g)b8
      aa      = 3.11e+05 * t9i23 * exp(-10.262*t9i13)
      daa     = aa*(-twoth*t9i + oneth*10.262*t9i43)

      bb      = 2.53e+03 * t9i32 * exp(-7.306*t9i)
      dbb     = bb*(-1.5d0*t9i + 7.306*t9i2)

      term    = aa + bb
      dtermdt = daa + dbb


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.30e+10 * t932 * exp(-1.595*t9i)
      drevdt   = rev*(1.5d0*t9i + 1.595*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_li7pag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,bb,dbb,cc,dcc,dd,ddd, &
                       t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56,zz, &
                       term1,dterm1,term2,dterm2,rev1,drev1dt,rev2,drev2dt,q1
      parameter        (q1 = 1.0d0/2.876416d0)



! 7li(p,g)8be=>2a
      aa   = 1.56d5 * t9i23 * exp(-8.472d0 * t9i13 - t92*q1)
      daa  = aa*(-twoth*t9i + oneth*8.472*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0 + 0.049d0*t913 + 2.498d0*t923 + 0.86d0*t9 &
             + 3.518d0*t943 + 3.08*t953
      dbb  = oneth*0.049*t9i23 + twoth*2.498*t9i13 + 0.86 &
             + fourth*3.518*t913 + fiveth*3.08*t923

      cc   = aa*bb
      dcc  = daa*bb + aa*dbb

      dd   =  1.55d6 * t9i32 * exp(-4.478d0 * t9i)
      ddd  = dd*(-1.5d0*t9i + 4.478*t9i2)

      term1  = cc + dd
      dterm1 = dcc + ddd

      rev1    = 6.55e+10 * t932 * exp(-200.225*t9i)
      drev1dt = rev1*(1.5d0*t9i + 200.225*t9i2)


! 7li(p,a)a
      aa     = 1.0d0 + 0.759*t9

      zz     = 1.0d0/aa
      t9a    = t9*zz
      dt9a   = (1.0d0 - t9a*0.759)*zz

      zz     = dt9a/t9a
      t9a13  = t9a**oneth
      dt9a13 = oneth*t9a13*zz

      t9a56  = t9a**fivsix
      dt9a56 = fivsix*t9a56*zz

      aa     = 1.096e+09 * t9i23 * exp(-8.472*t9i13)
      daa    = aa*(-twoth*t9i + oneth*8.472*t9i43)

      bb     = -4.830e+08 * t9a56 * t9i32 * exp(-8.472/t9a13)
      dbb    = bb*(dt9a56/t9a56 - 1.5d0*t9i + 8.472/t9a13**2*dt9a13)

      cc     = 1.06e+10 * t9i32 * exp(-30.442*t9i)
      dcc    = cc*(-1.5d0*t9i + 30.442*t9i2)

      term2   = aa + bb + cc
      dterm2  = daa + dbb + dcc

      rev2    = 4.69 * exp(-201.291*t9i)
      drev2dt = aa*201.291*t9i2




! sum the two forward rates (per f35 cf88)

      term = term1 + term2
      dtermdt = dterm1 + dterm2


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rr    = den * rev2 * term
      drrdt = den * (drev2dt*term + rev2*dtermdt) * 1.0d-9
      drrdd = rev2 * term

      return
      end




      subroutine rate_b8ep(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision lntwo,halflife,con
      parameter        (lntwo    = 0.6931471805599453d0, &
                        halflife = 0.77d0, &
                        con      = lntwo/halflife)


! b8(e+,nu)be8 => 2a

      fr    = con
      dfrdt = 0.0d0
      dfrdd = 0.0d0

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end






      subroutine rate_c12pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,q1
      parameter        (q1 = 1.0d0/2.25d0)


! c12(p,g)13n
      aa   = 2.04e+07 * t9i23 * exp(-13.69*t9i13 - t92*q1)
      daa  = aa*(-twoth*t9i + oneth*13.69*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0 + 0.03*t913 + 1.19*t923 + 0.254*t9 &
             + 2.06*t943 + 1.12*t953
      dbb  = oneth*0.03*t9i23 + twoth*1.19*t9i13 + 0.254 &
             + fourth*2.06*t913 + fiveth*1.12*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 1.08e+05 * t9i32 * exp(-4.925*t9i)
      ddd  = dd*(-1.5d0*t9i + 4.925*t9i2)

      ee   = 2.15e+05 * t9i32 * exp(-18.179*t9i)
      dee  = ee*(-1.5d0*t9i + 18.179*t9i2)

      term    = cc + dd + ee
      dtermdt = dcc + ddd + dee

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 8.84e+09 * t932 * exp(-22.553*t9i)
      drevdt   = rev*(1.5d0*t9i + 22.553*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_n13em(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision lntwo,halflife,con
      parameter        (lntwo    = 0.6931471805599453d0, &
                        halflife = 597.9d0, &
                        con      = lntwo/halflife)

! n13(e-nu)c13
      fr    = con
      dfrdt = 0.0d0
      dfrdd = 0.0d0

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end




      subroutine rate_c13pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,q1
      parameter        (q1 = 1.0d0/4.0d0)


! c13(p,g)13n
      aa   = 8.01e+07 * t9i23 * exp(-13.717*t9i13 - t92*q1)
      daa  = aa*(-twoth*t9i + oneth*13.717*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0 + 0.030*t913 + 0.958*t923 + 0.204*t9 &
             + 1.39*t943 + 0.753*t953
      dbb  = oneth*0.030*t9i23 + twoth*0.958*t9i13 + 0.204 &
             + fourth*1.39*t913 + fiveth*0.753*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 1.21e+06 * t9i65 * exp(-5.701*t9i)
      ddd  = dd*(-sixfif*t9i + 5.701*t9i2)

      term    = cc + dd
      dtermdt = dcc + ddd


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.19e+10 * t932 * exp(-87.621*t9i)
      drevdt   = rev*(1.5d0*t9i + 87.621*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_n14pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,q1
      parameter        (q1 = 1.0d0/10.850436d0)


! n14(p,g)o15
      aa  = 4.90e+07 * t9i23 * exp(-15.228*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*15.228*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0 + 0.027*t913 - 0.778*t923 - 0.149*t9 &
             + 0.261*t943 + 0.127*t953
      dbb  = oneth*0.027*t9i23 - twoth*0.778*t9i13 - 0.149 &
             + fourth*0.261*t913 + fiveth*0.127*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 2.37e+03 * t9i32 * exp(-3.011*t9i)
      ddd  = dd*(-1.5d0*t9i + 3.011*t9i2)

      ee   = 2.19e+04 * exp(-12.530*t9i)
      dee  = ee*12.530*t9i2

      term    = cc + dd + ee
      dtermdt = dcc + ddd + dee

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev    = 2.70e+10 * t932 * exp(-84.678*t9i)
      drevdt = rev*(1.5d0*t9i + 84.678*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end



      subroutine rate_o15em(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision lntwo,halflife,con
      parameter        (lntwo    = 0.6931471805599453d0, &
                        halflife = 122.24d0, &
                        con      = lntwo/halflife)

! o15(e-nu)n15
      fr    = con
      dfrdt = 0.0d0
      dfrdd = 0.0d0

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end






      subroutine rate_n14ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,ff,dff,q1
      parameter        (q1 = 1.0d0/0.776161d0)


! n14(a,g)f18
      aa  = 7.78d+09 * t9i23 * exp(-36.031*t9i13- t92*q1)
      daa = aa*(-twoth*t9i + oneth*36.031*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0 + 0.012*t913 + 1.45*t923 + 0.117*t9 &
             + 1.97*t943 + 0.406*t953
      dbb  = oneth*0.012*t9i23 + twoth*1.45*t9i13 + 0.117 &
             + fourth*1.97*t913 + fiveth*0.406*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 2.36d-10 * t9i32 * exp(-2.798*t9i)
      ddd  = dd*(-1.5d0*t9i + 2.798*t9i2)

      ee   = 2.03 * t9i32 * exp(-5.054*t9i)
      dee  = ee*(-1.5d0*t9i + 5.054*t9i2)

      ff   = 1.15d+04 * t9i23 * exp(-12.310*t9i)
      dff  = ff*(-twoth*t9i + 12.310*t9i2)

      term    = cc + dd + ee + ff
      dtermdt = dcc + ddd + dee + dff

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 5.42e+10 * t932 * exp(-51.236*t9i)
      drevdt   = rev*(1.5d0*t9i + 51.236*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_n15pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,ff,dff,q1
      parameter        (q1 = 1.0d0/0.2025d0)


! n15(p,g)o16
      aa  = 9.78e+08 * t9i23 * exp(-15.251*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*15.251*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0  + 0.027*t913 + 0.219*t923 + 0.042*t9 &
             + 6.83*t943 + 3.32*t953
      dbb  = oneth*0.027*t9i23 + twoth*0.219*t9i13 + 0.042 &
             + fourth*6.83*t913 + fiveth*3.32*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 1.11e+04*t9i32*exp(-3.328*t9i)
      ddd  = dd*(-1.5d0*t9i + 3.328*t9i2)

      ee   = 1.49e+04*t9i32*exp(-4.665*t9i)
      dee  = ee*(-1.5d0*t9i + 4.665*t9i2)

      ff   = 3.8e+06*t9i32*exp(-11.048*t9i)
      dff  = ff*(-1.5d0*t9i + 11.048*t9i2)

      term    = cc + dd + ee + ff
      dtermdt = dcc + ddd + dee + dff

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.62e+10 * t932 * exp(-140.734*t9i)
      drevdt   = rev*(1.5d0*t9i + 140.734*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_n15pa(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,ff,dff,gg,dgg, &
                       theta,q1
      parameter        (theta = 0.1d0, &
                        q1    = 1.0d0/0.272484d0)


! n15(p,a)c12
      aa  = 1.08d+12*t9i23*exp(-15.251*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*15.251*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0 + 0.027*t913 + 2.62*t923 + 0.501*t9 &
             + 5.36*t943 + 2.60*t953
      dbb  = oneth*0.027*t9i23 + twoth*2.62*t9i13 + 0.501 &
             + fourth*5.36*t913 + fiveth*2.60*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 1.19d+08 * t9i32 * exp(-3.676*t9i)
      ddd  = dd*(-1.5d0*t9i + 3.676*t9i2)

      ee   = 5.41d+08 * t9i12 * exp(-8.926*t9i)
      dee  = ee*(-0.5d0*t9i + 8.926*t9i2)

      ff   = theta * 4.72d+08 * t9i32 * exp(-7.721*t9i)
      dff  = ff*(-1.5d0*t9i + 7.721*t9i2)

      gg   = theta * 2.20d+09 * t9i32 * exp(-11.418*t9i)
      dgg  = gg*(-1.5d0*t9i + 11.418*t9i2)

      term    = cc + dd + ee + ff + gg
      dtermdt = dcc + ddd + dee + dff + dgg

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 7.06d-01*exp(-57.625*t9i)
      drevdt   = rev*57.625*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end




      subroutine rate_o16pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,zz


! o16(p,g)f17
      aa  = exp(-0.728*t923)
      daa = -twoth*aa*0.728*t9i13

      bb  = 1.0d0 + 2.13 * (1.0d0 - aa)
      dbb = -2.13*daa

      cc  = t923 * bb
      dcc = twoth*cc*t9i + t923*dbb

      dd   = exp(-16.692*t9i13)
      ddd  = oneth*dd*16.692*t9i43

      zz   = 1.0d0/cc
      ee   = dd*zz
      dee  = (ddd - ee*dcc)*zz

      term    = 1.50d+08 * ee
      dtermdt = 1.50d+08 * dee


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.03e+09*t932*exp(-6.968*t9i)
      drevdt   = rev*(1.5d0*t9i + 6.968*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_o17pa(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,res1,dres1,res2,dres2,res3,dres3, &
                       res4,dres4,res5,dres5,res6,dres6,zz, &
                       theta,q1,q2
      parameter        (theta = 0.1d0, &
                        q1    = 1.0d0/0.319225d0, &
                        q2    = 1.0d0/0.0016d0)



! o17(p,a)n14
! rate from jeff blackmons thesis, includes terms from fowler 75,
! landre 1990 (a&a 240, 85), and new results
! use rev factor from cf88 rate

      aa  = 1.53d+07 * t9i23 * exp(-16.712*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*16.712*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0 + 0.025*t913 + 5.39*t923 + 0.940*t9 &
             + 13.5*t943 + 5.98*t953
      dbb  = oneth*0.025*t9i23 + twoth*5.39*t9i13 + 0.940 &
             + fourth*13.5*t913 + fiveth*5.98*t923

      res1  = aa * bb
      dres1 = daa*bb + aa*dbb

      res2  = 2.92d+06 * t9 * exp(-4.247*t9i)
      dres2 = res2*(t9i + 4.247*t9i2)


      aa    = 0.479 * t923 + 0.00312
      daa   = twoth*0.479*t9i13

      bb    = aa*aa
      dbb   = 2.0d0 * aa * daa

      cc    =  1.78d+05 * t9i23 * exp(-16.669*t9i13)
      dcc   = cc*(-twoth*t9i + oneth*16.669*t9i43)

      zz    = 1.0d0/bb
      res3  = cc*zz
      dres3 = (dcc - res3*dbb)*zz

      res4  = 8.68d+10 * t9 * exp(-16.667*t9i13 - t92*q2)
      dres4 = res4*(t9i + oneth*16.667*t9i43 - 2.0d0*t9*q2)

      res5  = 9.22d-04 * t9i32 * exp(-0.767*t9i)
      dres5 = res5*(-1.5d0*t9i + 0.767*t9i2)

      res6  = theta * 98.0 * t9i32 * exp(-2.077*t9i)
      dres6 = res6*(-1.5d0*t9i + 2.077*t9i2)

      term    = res1 + res2 + res3 + res4 + res5 + res6
      dtermdt = dres1 + dres2 + dres3 + dres4 + dres5 + dres6

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 0.676 * exp(-13.825*t9i)
      drevdt   = rev*13.825*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end






      subroutine rate_o17pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,ff,dff,gg,dgg, &
                       t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56, &
                       zz,theta
      parameter        (theta = 0.1d0)


! o17(p,g)18f
! from landre et al 1990 a&a 240, 85
      aa     = 1.0d0 + 2.69*t9
      zz     = 1.0d0/aa

      t9a    = t9*zz
      dt9a   = (1.0d0 - t9a*2.69)*zz

      zz     = dt9a/t9a
      t9a13  = t9a**oneth
      dt9a13 = oneth*t9a13*zz

      t9a56  = t9a**fivsix
      dt9a56 = fivsix*t9a56*zz

      aa  = 7.97d+07 * t9a56 * t9i32 * exp(-16.712/t9a13)
      daa = aa*(dt9a56/t9a56 - 1.5d0*t9i + 16.712/t9a13**2*dt9a13)

      bb  = 1.0d0  + 0.025*t913 - 0.051*t923 - 8.82d-3*t9
      dbb = oneth*0.025*t9i23 - twoth*0.051*t9i13 - 8.82d-3
      if (bb .le. 0.0) then
       bb  = 0.0d0
       dbb = 0.0d0
      end if

      cc  = 1.51d+08 * t9i23 * exp(-16.712*t9i13)
      dcc = cc*(-twoth*t9i + oneth*16.712*t9i43)

      dd  = bb*cc
      ddd = dbb*cc + bb*dcc

      ee  = 1.56d+5 * t9i * exp(-6.272*t9i)
      dee = ee*(-t9i + 6.272*t9i2)

      ff  = 2.0d0 * theta * 3.16d-05 * t9i32 * exp(-0.767*t9i)
      dff = ff*(-1.5d0*t9i + 0.767*t9i2)

      gg  = theta * 98.0 * t9i32 * exp(-2.077*t9i)
      dgg = gg*(-1.5d0*t9i + 2.077*t9i2)

      term    = aa + dd + ee + ff + gg
      dtermdt = daa + ddd + dee + dff + dgg

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.66d+10 * t932 * exp(-65.061*t9i)
      drevdt   = rev*(1.5d0*t9i + 65.061*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_o18pa(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,ff,dff,gg,dgg,q1
      parameter        (q1 = 1.0d0/1.852321d0)


! o18(p,a)n15
      aa  = 3.63e+11 * t9i23 * exp(-16.729*t9i13 - t92*q1)
      daa = -twoth*aa*t9i + aa*(oneth*16.729*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.025*t913 + 1.88*t923 + 0.327*t9 &
            + 4.66*t943 + 2.06*t953
      dbb = oneth*0.025*t9i23 + twoth*1.88*t9i13 + 0.327 &
            + fourth*4.66*t913 + fiveth*2.06*t923

      cc  = aa * bb
      dcc = daa*bb + aa*dbb

      dd  = 9.90e-14 * t9i32 * exp(-0.231*t9i)
      ddd = -1.5d0*dd*t9i + dd*0.231*t9i2

      ee  = 2.66e+04 * t9i32 * exp(-1.670*t9i)
      dee = -1.5d0*ee*t9i + ee*1.670*t9i2

      ff  = 2.41e+09 * t9i32 * exp(-7.638*t9i)
      dff = -1.5d0*ff*t9i + ff*7.638*t9i2

      gg  = 1.46e+09 * t9i * exp(-8.310*t9i)
      dgg = -gg*t9i + gg*8.310*t9i2

      term    = cc + dd + ee + ff + gg
      dtermdt = dcc + ddd + dee + dff + dgg


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.66e-01 * exp(-46.191*t9i)
      drevdt   = rev*46.191*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end




      subroutine rate_o18pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,ff,dff,q1
      parameter        (q1 = 1.0d0/0.019321d0)


! o18(p,g)19f
      aa  = 3.45e+08 * t9i23 * exp(-16.729*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*16.729*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.025*t913 + 2.26*t923 + 0.394*t9 &
            + 30.56*t943 + 13.55*t953
      dbb = oneth*0.025*t9i23 + twoth*2.26*t9i13 + 0.394 &
            + fourth*30.56*t913 + fiveth*13.55*t923

      cc  = aa*bb
      dcc = daa*bb + aa*dbb

      dd  = 1.25e-15 * t9i32 * exp(-0.231*t9i)
      ddd = dd*(-1.5d0*t9i + 0.231*t9i2)

      ee  = 1.64e+02 * t9i32 * exp(-1.670*t9i)
      dee = ee*(-1.5d0*t9i + 1.670*t9i2)

      ff  = 1.28e+04 * t912 * exp(-5.098*t9i)
      dff = ff*(0.5d0*t9i + 5.098*t9i2)

      term    = cc + dd + ee + ff
      dtermdt = dcc + ddd + dee + dff


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 9.20e+09 * t932 * exp(-92.769*t9i)
      drevdt   = rev*(1.5d0*t9i + 92.769*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end







      subroutine rate_f17em(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision lntwo,halflife,con
      parameter        (lntwo    = 0.6931471805599453d0, &
                        halflife = 64.49d0, &
                        con      = lntwo/halflife)

! f17(e-nu)o17
      fr    = con
      dfrdt = 0.0d0
      dfrdd = 0.0d0

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end






      subroutine rate_f18em(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision lntwo,halflife,con
      parameter        (lntwo    = 0.6931471805599453d0, &
                        halflife = 6586.2d0, &
                        con      = lntwo/halflife)

! f18(e-nu)o18
      fr    = con
      dfrdt = 0.0d0
      dfrdd = 0.0d0

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end





      subroutine rate_f19pa(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,ff,dff,gg,dgg,hh,dhh,q1
      parameter        (q1 = 1.0d0/0.714025d0)


! f19(p,a)o16
      aa  = 3.55d+11 * t9i23 * exp(-18.113*t9i13 - t92*q1)
      daa = -twoth*aa*t9i + aa*(oneth*18.113*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.023*t913 + 1.96*t923 + 0.316*t9 &
            + 2.86*t943 + 1.17*t953
      dbb = oneth*0.023*t9i23 + twoth*1.96*t9i13 + 0.316 &
            + fourth*2.86*t913 + fiveth*1.17*t923

      cc  = aa * bb
      dcc = daa*bb + aa*dbb

      dd  = 3.67d+06 * t9i32 * exp(-3.752*t9i)
      ddd = -1.5d0*dd*t9i + dd*3.752*t9i2

      ee  = 3.07d+08 * exp(-6.019*t9i)
      dee = ee*6.019*t9i2

      ff  = 4.0*exp(-2.090*t9i)
      dff = ff*2.090*t9i2

      gg  = 7.0*exp(-16.440*t9i)
      dgg = gg*16.440*t9i2

      hh  = 1.0d0 + ff + gg
      dhh = dff + dgg

      term    = (cc + dd + ee)/hh
      dtermdt = ((dcc + ddd + dee) - term*dhh)/hh


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 6.54e-01 * exp(-94.159*t9i)
      drevdt   = rev*94.159*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_n13pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,q1
      parameter        (q1 = 1.0d0/0.69288976d0)


! n13(p,g)o14
! Keiner et al 1993 Nucl Phys A552, 66
      aa  = -1.727d+7 * t9i23 * exp(-15.168*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*15.168*t9i43 -2.0d0*t9*q1)

      bb  = 1.0d0 + 0.027*t913 - 17.54*t923 - 3.373*t9 &
            + 0.0176*t943 + 0.766d-2*t953
      dbb = oneth*0.027*t9i23 - twoth*17.54*t9i13 - 3.373 &
            + fourth*0.0176*t913 + fiveth*0.766d-2*t923

      cc  = aa*bb
      dcc = daa*bb + aa*dbb

      dd  = 3.1d+05 * t9i32 * exp(-6.348*t9i)
      ddd = dd*(-1.5d0*t9i + 6.348*t9i2)

      term    = cc + dd
      dtermdt = dcc + ddd

! goes negative below about t7=1.5
! note cf88 rate stays positive
      if (term .lt. 0.0) then
       term    = 0.0d0
       dtermdt = 0.0d0
      end if


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.57d+10*t932*exp(-53.706*t9i)
      drevdt   = rev*(1.5d0*t9i + 53.706*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_o14em(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision lntwo,halflife,con
      parameter        (lntwo    = 0.6931471805599453d0, &
                        halflife = 70.606d0, &
                        con      = lntwo/halflife)

! o14(e-nu)n14
      fr    = con
      dfrdt = 0.0d0
      dfrdd = 0.0d0

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end






      subroutine rate_o14ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,ff,dff,q1
      parameter        (q1 = 1.0d0/0.514089d0)


! o14(a,p)f17
      aa  = 1.68e+13 * t9i23 * exp(-39.388*t9i13- t92*q1)
      daa = -twoth*aa*t9i + aa*(oneth*39.388*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.011*t913 + 13.117*t923 + 0.971*t9 &
            + 85.295*t943 + 16.061*t953
      dbb = oneth*0.011*t9i23 + twoth*13.117*t9i13 + 0.971 &
            + fourth*85.295*t913 + fiveth*16.061*t923

      cc  = aa * bb
      dcc = daa*bb + aa*dbb

      dd  = 3.31e+04 * t9i32 * exp(-11.733*t9i)
      ddd = -1.5d0*dd*t9i + dd*11.733*t9i2

      ee  = 1.79e+07 * t9i32 * exp(-22.609*t9i)
      dee = -1.5d0*ee*t9i + ee*22.609*t9i2

      ff  = 9.00e+03 * t9113 * exp(-12.517*t9i)
      dff = elvnth*ff*t9i + ff*12.517*t9i2

      term    = cc + dd + ee + ff
      dtermdt = dcc + ddd + dee + dff


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 4.93e-01*exp(-13.820*t9i)
      drevdt   = rev*13.820*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end




      subroutine rate_o15ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,ff,dff,gg,dgg,hh,dhh, &
                       q1,q2,q3
      parameter        (q1 = 1.0d0/9.0d0, &
                        q2 = 1.0d0/3.751969d0, &
                        q3 = 1.0d0/64.0d0)


! o15(a,g)ne19

      aa  = 3.57d+11 * t9i23 * exp(-39.584d+0*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*39.584d0*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.011*t913 - 0.273*t923 - 0.020*t9
      dbb = oneth*0.011*t9i23 - twoth*0.273*t9i13 - 0.020

      cc  = aa*bb
      dcc = daa*bb + aa*dbb

      dd  = 5.10d+10 * t9i23 * exp(-39.584d+0*t9i13 - t92*q2)
      ddd = dd*(-twoth*t9i + oneth*39.584*t9i43 - 2.0d0*t9*q2)

      ee  = 1.0d0 + 0.011*t913 + 1.59*t923 + 0.117*t9 &
            + 1.81*t943 + 0.338*t953
      dee = oneth*0.011*t9i23 + twoth*1.59*t9i13 + 0.117 &
            + fourth*1.81*t913 + fiveth*0.338*t923

      ff  = dd*ee
      dff = ddd*ee + dd*dee

      gg  = 3.95d-1 * t9i32 * exp(-5.849*t9i)
      dgg = gg*(-1.5d0*t9i + 5.849*t9i2)

      hh  = 1.90d+1 * t9**2.85 * exp(-7.356*t9i - t92*q3)
      dhh = hh*(2.85*t9i + 7.356*t9i2 - 2.0d0*t9*q3)


      term    = cc + ff + gg + hh
      dtermdt = dcc + dff + dgg + dhh


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 5.54e+10 * t932 * exp(-40.957*t9i)
      drevdt   = rev*(1.5d0*t9i + 40.957*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_f17pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee


! f17(p,g)ne18
! wiescher and kettner, ap. j., 263, 891 (1982)

      aa  = 1.66e+07 * t9i23 * exp(-18.03*t9i13)
      daa = aa*(-twoth*t9i + oneth*18.03*t9i43)

      bb  = 2.194 + 0.050*t913 - 0.376*t923 - 0.061*t9 &
            + 0.026*t943 + 0.011*t953
      dbb = oneth*0.050*t9i23 - twoth*0.376*t9i13 - 0.061 &
            + fourth*0.026*t913 + fiveth*0.011*t923

      cc  = aa*bb
      dcc = daa*bb + aa*dbb

      dd  = 839.0 * t9i32 * exp(-6.93*t9i)
      ddd = dd*(-1.5d0*t9i + 6.93*t9i2)

      ee  = 33.56 * t9i32 * exp(-7.75*t9i)
      dee = ee*(-1.5d0*t9i + 7.75*t9i2)

      term    = cc + dd + ee
      dtermdt = dcc + ddd + dee

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.087e+11 * t932 * exp(-45.501*t9i)
      drevdt   = rev*(1.5d0*t9i + 45.501*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_ne18em(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision lntwo,halflife,con
      parameter        (lntwo    = 0.6931471805599453d0, &
                        halflife = 1.672d0, &
                        con      = lntwo/halflife)

! ne18(e-nu)f18
      fr    = con
      dfrdt = 0.0d0
      dfrdd = 0.0d0

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end





      subroutine rate_f18pa(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,ff,dff


! f18(p,a)o15
! wiescher and kettner, ap. j., 263, 891 (1982)

      aa  = 1.66e-10 * t9i32 * exp(-0.302*t9i)
      daa = aa*(-1.5d0*t9i + 0.302*t9i2)

      bb  = 1.56e+05 * t9i32 * exp(-3.84*t9i)
      dbb = bb*(-1.5d0*t9i + 3.84*t9i2)

      cc  = 1.36e+06 * t9i32 * exp(-5.22*t9i)
      dcc = cc*(-1.5d0*t9i + 5.22*t9i2)

      dd  = 8.1e-05 * t9i32 * exp(-1.05*t9i)
      ddd = dd*(-1.5d0*t9i + 1.05*t9i2)

      ee  = 8.9e-04 * t9i32 * exp(-1.51*t9i)
      dee = ee*(-1.5d0*t9i + 1.51*t9i2)

      ff  = 3.0e+05 * t9i32 * exp(-4.29*t9i)
      dff = ff*(-1.5d0*t9i + 4.29*t9i2)

      term    = aa + bb + cc + dd + ee + ff
      dtermdt = daa + dbb + dcc + ddd + dee + dff


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 4.93e-01 * exp(-33.433*t9i)
      drevdt   = rev*33.433*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_ne18ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,zz

      double precision z1,a1,ztot,ared,r,c1,c2,c3,c4
      parameter        (z1   = 10.0d0, &
                        a1   = 18.0d0, &
                        ztot = 2.0d0 * z1, &
                        ared = 4.0d0*a1/(4.0d0 + a1), &
                        r    = 5.1566081196876965d0, &
                        c1   = 4.9080044545315392d10, &
                        c2   = 4.9592784569936502d-2, &
                        c3   = 1.9288564401521285d1, &
                        c4   = 4.6477847042196437d1)

! note:
!      r    = 1.09 * a1**oneth + 2.3
!      c1   = 7.833e9 * 0.31 * ztot**fourth/(ared**fivsix)
!      c2   = 0.08617 * 0.1215 * sqrt(ared*r**3/ztot)
!      c3   = 2.0d0 * 0.52495 * sqrt(ared*r*ztot)
!      c4   = 4.2487 * (ztot**2*ared)**oneth


! ne18ap(a,p)na21
! was a call to aprate

      aa  = 1.0d0 + c2*t9
      zz  = c2/aa

      bb  = aa**fivsix
      dbb = fivsix*bb*zz

      cc  = t923 * bb
      dcc = twoth*cc*t9i + t923 * dbb

      dd = aa**oneth
      ddd = oneth*dd*zz

      ee  = t9i13 * dd
      dee = -oneth*ee*t9i + t9i13 * ddd

      zz      = 1.0d0/cc
      term    = c1*zz * exp(c3 - c4*ee)
      dtermdt = -term*(zz*dcc + c4*dee)

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev    = 0.0d0
      drevdt = 0.0d0

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end




      subroutine rate_ne19pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,ff,dff,gg,dgg,q1
      parameter        (q1 = 1.0d0/1.304164d0)


! ne19(p,g)na20

      aa  = 1.71d+6 * t9i23 * exp(-19.431d0*t9i13)
      daa = aa*(-twoth*t9i + oneth*19.431*t9i43)

      bb  = 1.0d0 + 0.021*t913 + 0.130*t923 + 1.95d-2*t9 &
            + 3.86d-2*t943 + 1.47d-02*t953
      dbb = oneth*0.021*t9i23 + twoth*0.130*t9i13 + 1.95d-2 &
            + fourth*3.86d-2*t913 + fiveth*1.47d-2*t923

      cc  = aa*bb
      dcc = daa*bb + aa*dbb


      dd  = 1.89d+5 * t9i23 * exp(-19.431d0*t9i13 - t92*q1)
      ddd = dd*(-twoth*t9i + oneth*19.431*t9i43 - 2.0d0*t9*q1)

      ee  = 1.0d0 + 0.021*t913 + 2.13*t923 + 0.320*t9 &
            + 2.80*t943 + 1.07*t953
      dee = oneth*0.021*t9i23 + twoth*2.13*t9i13 + 0.320 &
            + fourth*2.80*t913 + fiveth*1.07*t923

      ff  = dd*ee
      dff = ddd*ee + dd*dee

      gg  = 8.45d+3 * t9i54 * exp(-7.64d0*t9i)
      dgg = gg*(-fivfour*t9i + 7.64d0*t9i2)


      term    = cc + ff + gg
      dtermdt = dcc + dff + dgg


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 7.39e+09 * t932 * exp(-25.519*t9i)
      drevdt   = rev*(1.5d0*t9i + 25.519*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_ne19em(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision lntwo,halflife,con
      parameter        (lntwo    = 0.6931471805599453d0, &
                        halflife = 17.3982d0, &
                        con      = lntwo/halflife)

! ne18(e-nu)f18
      fr    = con
      dfrdt = 0.0d0
      dfrdd = 0.0d0

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end





      subroutine rate_si26ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,zz

      double precision z1,a1,ztot,ared,r,c1,c2,c3,c4
      parameter        (z1   = 14.0d0, &
                        a1   = 26.0d0, &
                        ztot = 2.0d0 * z1, &
                        ared = 4.0d0*a1/(4.0d0 + a1), &
                        r    = 5.5291207145640335d0, &
                        c1   = 7.3266779970543091d10, &
                        c2   = 4.7895369289991982d-02, &
                        c3   = 2.4322657793918662d1, &
                        c4   = 5.9292366232997814d1)

! note:
!      r    = 1.09 * a1**oneth + 2.3
!      c1   = 7.833e9 * 0.31 * ztot**fourth/(ared**fivsix)
!      c2   = 0.08617 * 0.1215 * sqrt(ared*r**3/ztot)
!      c3   = 2.0d0 * 0.52495 * sqrt(ared*r*ztot)
!      c4   = 4.2487 * (ztot**2*ared)**oneth



! si26ap(a,p)p29
! was a call to aprate

      aa  = 1.0d0 + c2*t9
      zz  = c2/aa

      bb  = aa**fivsix
      dbb = fivsix*bb*zz

      cc  = t923 * bb
      dcc = twoth*cc*t9i + t923 * dbb

      dd = aa**oneth
      ddd = oneth*dd*zz

      ee  = t9i13 * dd
      dee = -oneth*ee*t9i + t9i13 * ddd

      zz      = 1.0d0/cc
      term    = c1*zz * exp(c3 - c4*ee)
      dtermdt = -term*(zz*dcc + c4*dee)

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 0.0d0
      drevdt   = 0.0d0

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end





      subroutine rate_c12ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,gg,dgg,hh,dhh,f1,df1,f2,df2, &
                       zz,q1
      parameter        (q1 = 1.0d0/12.222016d0)


! c12(a,g)o16
      aa   = 1.0d0 + 0.0489d0*t9i23
      daa  = -twoth*0.0489d0*t9i53

      bb   = t92*aa*aa
      dbb  = 2.0d0*(bb*t9i + t92*aa*daa)

      cc   = exp(-32.120d0*t9i13 - t92*q1)
      dcc  = cc * (oneth*32.120d0*t9i43 - 2.0d0*t9*q1)

      dd   = 1.0d0 + 0.2654d0*t9i23
      ddd  = -twoth*0.2654d0*t9i53

      ee   = t92*dd*dd
      dee  = 2.0d0*(ee*t9i + t92*dd*ddd)

      ff   = exp(-32.120d0*t9i13)
      dff  = ff * oneth*32.120d0*t9i43

      gg   = 1.25d3 * t9i32 * exp(-27.499*t9i)
      dgg  = gg*(-1.5d0*t9i + 27.499*t9i2)

      hh   = 1.43d-2 * t95 * exp(-15.541*t9i)
      dhh  = hh*(5.0d0*t9i + 15.541*t9i2)

      zz   = 1.0d0/bb
      f1   = cc*zz
      df1  = (dcc - f1*dbb)*zz

      zz   = 1.0d0/ee
      f2   = ff*zz
      df2  = (dff - f2*dee)*zz

      term    = 1.04d8*f1  + 1.76d8*f2 + gg + hh
      dtermdt = 1.04d8*df1 + 1.76d8*df2 + dgg + dhh


! 1.7 times cf88 value
      term     = 1.7d0 * term
      dtermdt  = 1.7d0 * dtermdt

      fr    = term * den
      dfrdt = dtermdt * den * 1.0d-9
      dfrdd = term

      rev    = 5.13d10 * t932 * exp(-83.111*t9i)
      drevdt = rev*(1.5d0*t9i + 83.111*t9i2)

      rr     = rev * term
      drrdt  = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd  = 0.0d0

      return
      end






      subroutine rate_tripalf(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,r2abe,dr2abedt,rbeac, &
                       drbeacdt,aa,daa,bb,dbb,cc,dcc,dd,ddd,ee,dee, &
                       ff,dff,xx,dxx,yy,dyy,zz,dzz,uu,vv,f1,df1,rc28, &
                       q1,q2
      parameter        (rc28   = 0.1d0, &
                        q1     = 1.0d0/0.009604d0, &
                        q2     = 1.0d0/0.055225d0)



! triple alfa to c12
! this is a(a,g)be8
      aa    = 7.40d+05 * t9i32 * exp(-1.0663*t9i)
      daa   = aa*(-1.5d0*t9i  + 1.0663*t9i2)

      bb    = 4.164d+09 * t9i23 * exp(-13.49*t9i13 - t92*q1)
      dbb   = bb*(-twoth*t9i + oneth*13.49*t9i43 - 2.0d0*t9*q1)

      cc    = 1.0d0 + 0.031*t913 + 8.009*t923 + 1.732*t9 &
              + 49.883*t943 + 27.426*t953
      dcc   = oneth*0.031*t9i23 + twoth*8.009*t9i13 + 1.732 &
              + fourth*49.883*t913 + fiveth*27.426*t923

      r2abe    = aa + bb * cc
      dr2abedt = daa + dbb*cc + bb*dcc


! this is be8(a,g)c12
      dd    = 130.0d0 * t9i32 * exp(-3.3364*t9i)
      ddd   = dd*(-1.5d0*t9i + 3.3364*t9i2)

      ee    = 2.510d+07 * t9i23 * exp(-23.57*t9i13 - t92*q2)
      dee   = ee*(-twoth*t9i + oneth*23.57*t9i43 - 2.0d0*t9*q2)

      ff    = 1.0d0 + 0.018*t913 + 5.249*t923 + 0.650*t9 + &
              19.176*t943 + 6.034*t953
      dff   = oneth*0.018*t9i23 + twoth*5.249*t9i13 + 0.650 &
              + fourth*19.176*t913 + fiveth*6.034*t923

      rbeac    = dd + ee * ff
      drbeacdt = ddd + dee * ff + ee * dff


! a factor
      xx    = rc28 * 1.35d-07 * t9i32 * exp(-24.811*t9i)
      dxx   = xx*(-1.5d0*t9i + 24.811*t9i2)


! high temperature rate
      if (t9.gt.0.08) then
       term    = 2.90d-16 * r2abe * rbeac + xx
       dtermdt =   2.90d-16 * dr2abedt * rbeac &
                 + 2.90d-16 * r2abe * drbeacdt &
                 + dxx


! low temperature rate
      else
       uu   = 0.8d0*exp(-(0.025*t9i)**3.263)
       yy   = 0.2d0 + uu
!       yy   = 0.01 + 0.2d0 + uu
       dyy  = uu * 3.263*(0.025*t9i)**2.263 * (0.025*t9i2)
       vv   = 4.0d0*exp(-(t9/0.025)**9.227)
       zz   = 1.0d0 + vv
       dzz  = vv * 9.227*(t9/0.025)**8.227 * 40.0d0
       aa   = 1.0d0/zz
       f1   = 0.01d0 + yy * aa
!       f1   = yy * aa
       df1  = (dyy - f1*dzz)*aa
       term = 2.90d-16 * r2abe * rbeac * f1 +  xx
       dtermdt =   2.90d-16 * dr2abedt * rbeac * f1 &
                 + 2.90d-16 * r2abe * drbeacdt * f1 &
                 + 2.90d-16 * r2abe * rbeac * df1 &
                 + dxx
      end if


! rates
!      term    = 1.2d0 * term
!      dtermdt = 1.2d0 * term

      fr    = term * den * den
      dfrdt = dtermdt * den * den * 1.0d-9
      dfrdd = 2.0d0 * term * den

      rev    = 2.00d+20*t93*exp(-84.424*t9i)
      drevdt = rev*(3.0d0*t9i + 84.424*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end







      subroutine rate_c12c12(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56, &
                       aa,zz


! c12 + c12 reaction
      aa      = 1.0d0 + 0.0396*t9
      zz      = 1.0d0/aa

      t9a     = t9*zz
      dt9a    = (1.0d0 -  t9a*0.0396)*zz

      zz      = dt9a/t9a
      t9a13   = t9a**oneth
      dt9a13  = oneth*t9a13*zz

      t9a56   = t9a**fivsix
      dt9a56  = fivsix*t9a56*zz

      term    = 4.27d+26 * t9a56 * t9i32 * &
                exp(-84.165/t9a13 - 2.12d-03*t93)
      dtermdt = term*(dt9a56/t9a56 - 1.5d0*t9i &
                      + 84.165/t9a13**2*dt9a13 - 6.36d-3*t92)

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end






      subroutine rate_c12c12npa(temp,den, &
                      fr1,dfr1dt,dfr1dd,rr1,drr1dt,drr1dd, &
                      fr2,dfr2dt,dfr2dd,rr2,drr2dt,drr2dd, &
                      fr3,dfr3dt,dfr3dd,rr3,drr3dt,drr3dd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den, &
                      fr1,dfr1dt,dfr1dd,rr1,drr1dt,drr1dd, &
                      fr2,dfr2dt,dfr2dd,rr2,drr2dt,drr2dd, &
                      fr3,dfr3dt,dfr3dd,rr3,drr3dt,drr3dd

! locals
      double precision term,dtermdt,rev,drevdt,t9a,dt9a,t9a13,dt9a13, &
                       t9a56,dt9a56,aa,daa,bb,dbb,cc,dcc,dd,ddd,zz, &
                       b24n,db24n,b24p,db24p,b24a,db24a


! c12(c12,n)mg23
! c12(c12,p)na23
! c12(c12,a)ne20


      aa      = 1.0d0 + 0.0396*t9
      zz      = 1.0d0/aa

      t9a     = t9 * zz
      dt9a    = (1.0d0 -  t9a*0.0396) * zz

      zz      = dt9a/t9a
      t9a13   = t9a**oneth
      dt9a13  = oneth*t9a13*zz

      t9a56   = t9a**fivsix
      dt9a56  = fivsix*t9a56*zz

      aa = 4.27d+26 * t9a56 * t9i32 * &
           exp(-84.165/t9a13 - 2.12d-03*t93)

      daa = aa * (dt9a56/t9a56 - 1.5d0*t9i &
                + 84.165/t9a13**2*dt9a13 - 6.36d-3*t92)


! neutron branching from dayras switkowski and woosley 1976
      if (t9 .ge. 1.5) then

       bb    =  0.055 * exp(0.976 - 0.789*t9)
       dbb   = -bb*0.789

       b24n  = 0.055  - bb
       db24n = -dbb

      else

       bb    = 1.0d0 + 0.0789*t9 + 7.74*t92
       dbb   = 0.0789 + 2.0d0*7.74*t9

       cc    = 0.766*t9i3
       dcc   = -3.0d0*cc*t9i

       dd    = bb * cc
       ddd   = dbb*cc + bb*dcc

       b24n  = 0.859*exp(-dd)
       db24n = -b24n*ddd
      end if


! proton branching ratio
      if (t9.gt.3.) then
        b24p  = oneth*(1.0d0 - b24n)
        db24p = -oneth*db24n

        b24a  = 2.0d0 * b24p
        db24a = 2.0d0 * db24p

       else
        b24p  = 0.5d0*(1.0d0 - b24n)
        db24p = -0.5d0*db24n

        b24a  = b24p
        db24a = db24p

       end if

! rates

! c12(c12,n)mg23
      term    = aa * b24n
      dtermdt = daa*b24n + aa*db24n
      fr1     = den * term
      dfr1dt  = den * dtermdt * 1.0d-9
      dfr1dd  = term

      rev    = 0.0d0
      drevdt = 0.0d0
      if (t9 .gt. 0.1) then
       rev    = 3.93 * exp(30.16100515d0*t9i)
       drevdt = -rev*30.16100515d0*t9i2
      end if
      rr1    = den * rev * term
      drr1dt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drr1dd = rev*term


! c12(c12,p)na23
      term    = aa * b24p
      dtermdt = daa*b24p + aa*db24p
      fr2     = den * term
      dfr2dt  = den * dtermdt * 1.0d-9
      dfr2dd  = term

      rev    = 0.0d0
      drevdt = 0.0d0
      if (t9 .gt. 0.1) then
       rev    = 3.93 * exp(-25.98325915d0*t9i)
       drevdt = rev*25.98325915d0*t9i2
      end if
      rr2    = den * rev * term
      drr2dt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drr2dd = rev*term


! c12(c12,a)ne20
      term    = aa * b24a
      dtermdt = daa*b24a + aa*db24a
      fr3     = den * term
      dfr3dt  = den * dtermdt * 1.0d-9
      dfr3dd  = term

      rev    = 0.0d0
      drevdt = 0.0d0
      if (t9 .gt. 0.1) then
       rev    = 2.42 * exp(-53.576110995d0*t9i)
       drevdt = rev*53.576110995d0*t9i2
      end if
      rr3    = den * rev * term
      drr3dt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drr3dd = rev*term

      return
      end






      subroutine rate_c12o16(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,t9a,dt9a,t9a13,dt9a13,t9a23,dt9a23, &
                       t9a56,dt9a56,aa,daa,bb,dbb,cc,dcc,zz


! c12 + o16 reaction; see cf88 references 47-4

      if (t9.ge.0.5) then
       aa     = 1.0d0 + 0.055*t9
       zz     = 1.0d0/aa

       t9a    = t9*zz
       dt9a   = (1.0d0 - t9a*0.055)*zz

       zz     = dt9a/t9a
       t9a13  = t9a**oneth
       dt9a13 = oneth*t9a13*zz

       t9a23  = t9a13*t9a13
       dt9a23 = 2.0d0 * t9a13 * dt9a13

       t9a56  = t9a**fivsix
       dt9a56 = fivsix*t9a56*zz

       aa      = exp(-0.18*t9a*t9a)
       daa     = -aa * 0.36 * t9a * dt9a

       bb      = 1.06d-03*exp(2.562*t9a23)
       dbb     = bb * 2.562 * dt9a23

       cc      = aa + bb
       dcc     = daa + dbb

       zz      = 1.0d0/cc
       term    = 1.72d+31 * t9a56 * t9i32 * exp(-106.594/t9a13) * zz
       dtermdt = term*(dt9a56/t9a56 - 1.5d0*t9i &
                       + 106.594/t9a23*dt9a13 - zz*dcc)

      else
!       term    = 2.6288035d-29
       term    = 0.0d0
       dtermdt = 0.0d0
      endif


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end




      subroutine rate_c12o16npa(temp,den, &
                      fr1,dfr1dt,dfr1dd,rr1,drr1dt,drr1dd, &
                      fr2,dfr2dt,dfr2dd,rr2,drr2dt,drr2dd, &
                      fr3,dfr3dt,dfr3dd,rr3,drr3dt,drr3dd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                      fr1,dfr1dt,dfr1dd,rr1,drr1dt,drr1dd, &
                      fr2,dfr2dt,dfr2dd,rr2,drr2dt,drr2dd, &
                      fr3,dfr3dt,dfr3dd,rr3,drr3dt,drr3dd

! locals
      double precision term,dtermdt,rev,drevdt,t9a,dt9a,t9a13,dt9a13, &
                       t9a23,dt9a23,t9a56,dt9a56,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,b27n,b27p,b24a,zz


! c12(o16,n)si27
! c12(o16,p)al27
! c12(o16,a)mg24


      if (t9.ge.0.5) then
       aa     = 1.0d0 + 0.055*t9
       zz     = 1.0d0/aa

       t9a    = t9*zz
       dt9a   = (1.0d0 - t9a*0.055)*zz

       zz     = dt9a/t9a
       t9a13  = t9a**oneth
       dt9a13 = oneth*t9a13*zz

       t9a23  = t9a13*t9a13
       dt9a23 = 2.0d0 * t9a13 * dt9a13

       t9a56  = t9a**fivsix
       dt9a56 = fivsix*t9a56*zz

       aa     = exp(-0.18*t9a*t9a)
       daa    = -aa * 0.36 * t9a * dt9a

       bb     = 1.06d-03*exp(2.562*t9a23)
       dbb    = bb * 2.562 * dt9a23

       cc     = aa + bb
       dcc    = daa + dbb

       zz     = 1.0d0/cc
       dd     = 1.72d+31 * t9a56 * t9i32 * exp(-106.594/t9a13) *zz
       ddd    = dd*(dt9a56/t9a56 - 1.5d0*t9i &
                 + 106.594/t9a23 * dt9a13 - dcc*zz)

      else
!       dd     = 2.6288035d-29
       dd     = 0.0d0
       ddd    = 0.0d0
      endif


! branching ratios from pwnsz data
        b27n = 0.1d0
        b27p = 0.5d0
        b24a = 0.4d0


! rates

! c12(o16,n)si27
      term    = dd * b27n
      dtermdt = ddd * b27n
      fr1     = den * term
      dfr1dt  = den * dtermdt * 1.0d-9
      dfr1dd  = term

      rev    = 0.0d0
      drevdt = 0.0d0
      if (t9 .gt. 0.1) then
       rev    = 1.58d0 * exp(4.8972467d0*t9i)
       drevdt = -rev*4.8972467d0*t9i2
      end if
      rr1    = den * rev * term
      drr1dt = den*(drevdt*term + rev*dtermdt)*1.0d-9
      drr1dd = rev*term


! c12(o16,p)al27
      term    = dd * b27p
      dtermdt = ddd * b27p
      fr2     = den * term
      dfr2dt  = den * dtermdt * 1.0d-9
      dfr2dd  = term

      rev    = 0.0d0
      drevdt = 0.0d0
      if (t9 .gt. 0.1) then
       rev    = 1.58d0 * exp(-59.9970745d0*t9i)
       drevdt = rev*59.9970745d0*t9i2
      end if
      rr2    = den * rev * term
      drr2dt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drr2dd = rev*term


! c12(o16,a)mg24
      term    = dd * b24a
      dtermdt = ddd * b24a
      fr3     = den * term
      dfr3dt  = den * dtermdt * 1.0d-9
      dfr3dd  = term

      rev    = 0.0d0
      drevdt = 0.0d0
      if (t9 .gt. 0.1) then
       rev    = 2.83d0 * exp(-78.5648345d0*t9i)
       drevdt = rev*78.5648345d0*t9i2
      end if
      rr3    = den * rev * term
      drr3dt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drr3dd = rev*term


      return
      end








      subroutine rate_o16o16(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt


! o16 + o16
      term  = 7.10d36 * t9i23 * &
              exp(-135.93 * t9i13 - 0.629*t923 &
                   - 0.445*t943 + 0.0103*t9*t9)

      dtermdt = -twoth*term*t9i &
                + term * (oneth*135.93*t9i43 - twoth*0.629*t9i13 &
                          - fourth*0.445*t913 + 0.0206*t9)


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end






      subroutine rate_o16o16npad(temp,den, &
                      fr1,dfr1dt,dfr1dd,rr1,drr1dt,drr1dd, &
                      fr2,dfr2dt,dfr2dd,rr2,drr2dt,drr2dd, &
                      fr3,dfr3dt,dfr3dd,rr3,drr3dt,drr3dd, &
                      fr4,dfr4dt,dfr4dd,rr4,drr4dt,drr4dd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                      fr1,dfr1dt,dfr1dd,rr1,drr1dt,drr1dd, &
                      fr2,dfr2dt,dfr2dd,rr2,drr2dt,drr2dd, &
                      fr3,dfr3dt,dfr3dd,rr3,drr3dt,drr3dd, &
                      fr4,dfr4dt,dfr4dd,rr4,drr4dt,drr4dd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa, &
                       b32n,db32n,b32p,db32p,b32a,db32a,b32d,db32d, &
                       ezro,dezro,dlt,ddlt,xxt,dxxt,thrs,dthrs


! o16(o16,n)s31
! o16(o16,p)p31
! o16(o16,a)si28
! o16(o16,d)p30

      aa  = 7.10d36 * t9i23 * exp(-135.93 * t9i13 - 0.629*t923 &
                   - 0.445*t943 + 0.0103*t9*t9)

      daa = -twoth*aa*t9i + aa * (oneth*135.93*t9i43 - twoth*0.629*t9i13 &
                          - fourth*0.445*t913 + 0.0206*t9)


! branching ratios highly uncertain;  guessed using fcz 1975
! deuteron channel is endoergic. apply error function cut-off.
       ezro = 3.9*t923
       dezro = twoth*ezro*t9i

       dlt  = 1.34*t9**fivsix
       ddlt = fivsix*dlt*t9i

       xxt  = 2.0d0*(2.406 - ezro)/dlt
       dxxt = -(2.0d0*dezro + xxt*ddlt)/dlt

       call fowthrsh(xxt,thrs,dthrs)

       b32d  = 0.05d0*thrs
       db32d = 0.05d0*dthrs*dxxt

! shut down deutrium channel
!       b32d  = 0.0d0
!       db32d = 0.0d0

       b32n  = 0.1d0
       db32n = 0.0d0

       b32a  = 0.25d0
       db32a = 0.0d0

       b32p  = 1.0d0 - b32d - b32a - b32n
       db32p = -db32d



! rates

! o16(o16,n)s31
      term    = aa * b32n
      dtermdt = daa*b32n + aa*db32n
      fr1     = den * term
      dfr1dt  = den * dtermdt * 1.0d-9
      dfr1dd  = term

      rev    = 0.0d0
      drevdt = 0.0d0
      if (t9 .gt. 0.1) then
       rev    = 5.92 * exp(-16.8038228d0*t9i)
       drevdt = rev*16.8038228d0*t9i2
      end if
      rr1    = den * rev * term
      drr1dt = den*(drevdt*term + rev*dtermdt) * 1.0d-9
      drr1dd = rev*term

! o16(o16,p)p31
      term    = aa * b32p
      dtermdt = daa*b32p + aa*db32p
      fr2     = den * term
      dfr2dt  = den * dtermdt * 1.0d-9
      dfr2dd  = term

      rev    = 0.0d0
      drevdt = 0.0d0
      if (t9 .gt. 0.1) then
       rev    = 5.92*exp(-89.0788286d0*t9i)
       drevdt = rev*89.0788286d0*t9i2
      end if
      rr2    = den * rev * term
      drr2dt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drr2dd = rev*term


! o16(o16,a)si28
      term    = aa * b32a
      dtermdt = daa*b32a + aa*db32a
      fr3     = den * term
      dfr3dt  = den * dtermdt * 1.0d-9
      dfr3dd  = term

      rev    = 0.0d0
      drevdt = 0.0d0
      if (t9 .gt. 0.1) then
       rev    = 3.46*exp(-111.3137212d0*t9i)
       drevdt = rev*111.3137212d0*t9i2
      end if
      rr3    = den * rev * term
      drr3dt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drr3dd = rev*term


! o16(o16,d)p30
      term    = aa * b32d
      dtermdt = daa*b32d + aa*db32d
      fr4     = den * term
      dfr4dt  = den * dtermdt * 1.0d-9
      dfr4dd  = term

      rev    = 0.0d0
      drevdt = 0.0d0
      if (t9 .gt. 0.1) then
       rev    = 0.984*exp(27.9908982d0*t9i)
       drevdt = -rev*27.9908982d0*t9i2
      end if
      rr4    = den * rev * term
      drr4dt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drr4dd = rev*term

      return
      end






      subroutine fowthrsh(x,thrs,dthrs)
      include 'implno.dek'

! fowler threshold fudge function.
! err func rational (abramowitz p.299)7.1.25 and its derivative

! declare
      double precision x,thrs,dthrs, &
                       ag,z,z2,t,t2,t3,tt,er,aa,daa,dt,dtt,der

      ag   = sign(1.0d0,x)
      z    = abs(x)
      z2   = z*z
      aa   = 1.0d0 + 0.47047d0*z

      t    = 1.0d0/aa
      dt   = -t*t*0.47047*ag

      t2   = t*t
      t3   = t2*t

      tt   = 0.3480242d0*t - 0.0958798d0*t2 + 0.7478556d0*t3
      dtt  = dt * (0.3480242d0 - 2.0d0*0.0958798d0*t &
                   + 3.0d0*0.7478556d0*t2)

      thrs  = 0.5d0
      dthrs = -5.6452433937900004d-1

      if (z .ne. 0) then
       aa   = exp(-z2)
       daa  = -2.0d0*aa*z*ag

       er   = 1.0d0 - tt * aa
       der  = -dtt*aa - tt*daa

       thrs  = 0.5d0 * (1.0d0 - ag*er)
       dthrs = -0.5d0*ag*der
      end if
      return
      end





      subroutine rate_o16ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,term1,dterm1,aa,daa,bb,dbb, &
                       cc,dcc,term2,dterm2,rev,drevdt,q1
      parameter        (q1 = 1.0d0/2.515396d0)



! o16(a,g)ne20
      term1   = 9.37d9 * t9i23 * exp(-39.757*t9i13 - t92*q1)
      dterm1  = term1*(-twoth*t9i + oneth*39.757*t9i43 - 2.0d0*t9*q1)

      aa      = 62.1 * t9i32 * exp(-10.297*t9i)
      daa     = aa*(-1.5d0*t9i + 10.297*t9i2)

      bb      = 538.0d0 * t9i32 * exp(-12.226*t9i)
      dbb     = bb*(-1.5d0*t9i + 12.226*t9i2)

      cc      = 13.0d0 * t92 * exp(-20.093*t9i)
      dcc     = cc*(2.0d0*t9i + 20.093*t9i2)

      term2   = aa + bb + cc
      dterm2  = daa + dbb + dcc

      term    = term1 + term2
      dtermdt = dterm1 + dterm2


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 5.65d+10*t932*exp(-54.937*t9i)
      drevdt   = rev*(1.5d0*t9i + 54.937*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_ne20ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,term1,dterm1,aa,daa,bb,dbb, &
                       term2,dterm2,term3,dterm3,rev,drevdt,zz,rc102,q1
      parameter        (rc102 = 0.1d0, &
                        q1    = 1.0d0/4.923961d0)



! ne20(a,g)mg24

      aa   = 4.11d+11 * t9i23 * exp(-46.766*t9i13 - t92*q1)
      daa  = aa*(-twoth*t9i + oneth*46.766*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0 + 0.009*t913 + 0.882*t923 + 0.055*t9 &
             + 0.749*t943 + 0.119*t953
      dbb  = oneth*0.009*t9i23 + twoth*0.882*t9i13 + 0.055 &
             + fourth*0.749*t913 + fiveth*0.119*t923

      term1  = aa * bb
      dterm1 = daa * bb + aa * dbb


      aa   = 5.27d+03 * t9i32 * exp(-15.869*t9i)
      daa  = aa*(-1.5d0*t9i + 15.869*t9i2)

      bb   = 6.51d+03 * t912 * exp(-16.223*t9i)
      dbb  = bb*(0.5d0*t9i + 16.223*t9i2)

      term2  = aa + bb
      dterm2 = daa + dbb


      aa   = 42.1 * t9i32 * exp(-9.115*t9i)
      daa  = aa*(-1.5d0*t9i + 9.115*t9i2)

      bb   =  32.0 * t9i23 * exp(-9.383*t9i)
      dbb  = bb*(-twoth*t9i + 9.383*t9i2)

      term3  = rc102 * (aa + bb)
      dterm3 = rc102 * (daa + dbb)


      aa  = 5.0d0*exp(-18.960*t9i)
      daa = aa*18.960*t9i2

      bb  = 1.0d0 + aa
      dbb = daa

      zz      = 1.0d0/bb
      term    = (term1 + term2 + term3)*zz
      dtermdt = ((dterm1 + dterm2 + dterm3) - term*dbb)*zz


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 6.01d+10 * t932 * exp(-108.059*t9i)
      drevdt   = rev*(1.5d0*t9i + 108.059*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_mg24ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,bb,dbb,cc,dcc,dd,ddd,ee,dee, &
                       ff,dff,gg,dgg,hh,hhi,rev,drevdt,rc121
      parameter        (rc121 = 0.1d0)



! 24mg(a,g)28si

      aa    = 4.78d+01 * t9i32 * exp(-13.506*t9i)
      daa   = aa*(-1.5d0*t9i + 13.506*t9i2)

      bb    =  2.38d+03 * t9i32 * exp(-15.218*t9i)
      dbb   = bb*(-1.5d0*t9i + 15.218*t9i2)

      cc    = 2.47d+02 * t932 * exp(-15.147*t9i)
      dcc   = cc*(1.5d0*t9i + 15.147*t9i2)

      dd    = rc121 * 1.72d-09 * t9i32 * exp(-5.028*t9i)
      ddd   = dd*(-1.5d0*t9i + 5.028*t9i2)

      ee    = rc121* 1.25d-03 * t9i32 * exp(-7.929*t9i)
      dee   = ee*(-1.5d0*t9i + 7.929*t9i2)

      ff    = rc121 * 2.43d+01 * t9i * exp(-11.523*t9i)
      dff   = ff*(-t9i + 11.523*t9i2)

      gg    = 5.0d0*exp(-15.882*t9i)
      dgg   = gg*15.882*t9i2

      hh    = 1.0d0 + gg
      hhi   = 1.0d0/hh

      term    = (aa + bb + cc + dd + ee + ff) * hhi
      dtermdt = (daa + dbb + dcc + ddd + dee + dff - term*dgg) * hhi


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 6.27d+10 * t932 * exp(-115.862*t9i)
      drevdt   = rev*(1.5d0*t9i + 115.862*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_mg24ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,bb,dbb,cc,dcc,dd,ddd,ee,dee, &
                       ff,dff,gg,dgg,term1,dterm1,term2,dterm2, &
                       rev,drevdt,rc148,q1
      parameter        (rc148 = 0.1d0, &
                        q1    = 1.0d0/0.024649d0)



! 24mg(a,p)al27
      aa     = 1.10d+08 * t9i23 * exp(-23.261*t9i13 - t92*q1)
      daa    = -twoth*aa*t9i + aa*(23.261*t9i43 - 2.0d0*t9*q1)

      bb     =  1.0d0 + 0.018*t913 + 12.85*t923 + 1.61*t9 &
               + 89.87*t943 + 28.66*t953
      dbb    = oneth*0.018*t9i23 + twoth*12.85*t9i13 + 1.61 &
                + fourth*89.87*t913 + fiveth*28.66*t923

      term1  = aa * bb
      dterm1 = daa * bb + aa * dbb

      aa     = 129.0d0 * t9i32 * exp(-2.517*t9i)
      daa    = -1.5d0*aa*t9i + aa*2.517*t9i2

      bb     = 5660.0d0 * t972 * exp(-3.421*t9i)
      dbb    = 3.5d0*bb*t9i +  bb*3.421*t9i2

      cc     = rc148 * 3.89d-08 * t9i32 * exp(-0.853*t9i)
      dcc    = -1.5d0*cc*t9i + cc*0.853*t9i2

      dd     = rc148 * 8.18d-09 * t9i32 * exp(-1.001*t9i)
      ddd    = -1.5d0*dd*t9i + dd*1.001*t9i2

      term2  = aa + bb + cc + dd
      dterm2 = daa + dbb + dcc + ddd

      ee     = oneth*exp(-9.792*t9i)
      dee    = ee*9.792*t9i2

      ff     =  twoth * exp(-11.773*t9i)
      dff    = ff*11.773*t9i2

      gg     = 1.0d0 + ee + ff
      dgg    = dee + dff

      term    = (term1 + term2)/gg
      dtermdt = ((dterm1 + dterm2) - term*dgg)/gg


! the rates
      rev      = 1.81 * exp(-18.572*t9i)
      drevdt   = rev*18.572*t9i2

      fr    = den * rev * term
      dfrdt = den * (drevdt * term + rev * dtermdt) * 1.0d-9
      dfrdd = rev * term

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end






      subroutine rate_al27pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,gg,dgg


! al27(p,g)si28
! champagne 1996

      aa  = 1.32d+09 * t9i23 * exp(-23.26*t9i13)
      daa = aa*(-twoth*t9i + oneth*23.26*t9i43)

      bb  = 3.22d-10 * t9i32 * exp(-0.836*t9i)*0.17
      dbb = bb*(-1.5d0*t9i + 0.836*t9i2)

      cc  = 1.74d+00 * t9i32 * exp(-2.269*t9i)
      dcc = cc*(-1.5d0*t9i + 2.269*t9i2)

      dd  = 9.92d+00 * t9i32 * exp(-2.492*t9i)
      ddd = dd*(-1.5d0*t9i + 2.492*t9i2)

      ee  = 4.29d+01 * t9i32 * exp(-3.273*t9i)
      dee = ee*(-1.5d0*t9i + 3.273*t9i2)

      ff  = 1.34d+02 * t9i32 * exp(-3.654*t9i)
      dff = ff*(-1.5d0*t9i + 3.654*t9i2)

      gg  = 1.77d+04 * (t9**0.53) * exp(-4.588*t9i)
      dgg = gg*(0.53*t9i + 4.588*t9i2)

      term    = aa + bb + cc + dd + ee + ff + gg
      dtermdt = daa + dbb + dcc + ddd + dee + dff + dgg


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.13d+11 * t932 * exp(-134.434*t9i)
      drevdt   = rev*(1.5d0*t9i + 134.434*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end







      subroutine rate_al27pg_old(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,bb,dbb,cc,dcc,dd,ddd,ee,dee, &
                       ff,dff,gg,dgg,hh,dhh,xx,dxx,yy,dyy,zz,dzz,pp, &
                       rev,drevdt,rc147,q1
      parameter        (rc147 = 0.1d0, &
                        q1    = 1.0d0/0.024025d0)



! 27al(p,g)si28  cf88

      aa  = 1.67d+08 * t9i23 * exp(-23.261*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*23.261*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.018*t913 + 5.81*t923 + 0.728*t9 &
            + 27.31*t943 + 8.71*t953
      dbb = oneth*0.018*t9i23 + twoth*5.81*t9i13 + 0.728 &
            + fourth*27.31*t913 + fiveth*8.71*t923

      cc  = aa*bb
      dcc = daa*bb + aa*dbb

      dd  = 2.20d+00 * t9i32 * exp(-2.269*t9i)
      ddd = dd*(-1.5d0*t9i + 2.269*t9i2)

      ee  = 1.22d+01 * t9i32 * exp(-2.491*t9i)
      dee = ee*(-1.5d0*t9i + 2.491*t9i2)

      ff  =  1.50d+04 * t9 * exp(-4.112*t9i)
      dff = ff*(t9i + 4.112*t9i2)

      gg  = rc147 * 6.50d-10 * t9i32 * exp(-0.853*t9i)
      dgg = gg*(-1.5d0*t9i + 0.853*t9i2)

      hh  = rc147 * 1.63d-10 * t9i32 * exp(-1.001*t9i)
      dhh = hh*(-1.5d0*t9i + 1.001*t9i2)

      xx     = oneth*exp(-9.792*t9i)
      dxx    = xx*9.792*t9i2

      yy     =  twoth * exp(-11.773*t9i)
      dyy    = yy*11.773*t9i2

      zz     = 1.0d0 + xx + yy
      dzz    = dxx + dyy

      pp      = 1.0d0/zz
      term    = (cc + dd + ee + ff + gg + hh)*pp
      dtermdt = ((dcc + ddd + dee + dff + dgg + dhh) - term*dzz)*pp


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.13d+11*t932*exp(-134.434*t9i)
      drevdt   = rev*(1.5d0*t9i + 134.434*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_si28ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! si28(a,g)s32
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 6.340d-2*z + 2.541d-3*z2 - 2.900d-4*z3
      if (z .eq. 10.0) then
       daa = 0
      else
       daa   = 6.340d-2 + 2.0d0*2.541d-3*t9 - 3.0d0*2.900d-4*t92
      end if

      term    = 4.82d+22 * t9i23 * exp(-61.015 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 61.015*t9i13*(oneth*t9i*aa - daa))

! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 6.461d+10 * t932 * exp(-80.643*t9i)
      drevdt   = rev*(1.5d0*t9i + 80.643*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_si28ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! si28(a,p)p31
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 2.798d-3*z + 2.763d-3*z2 - 2.341d-4*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 2.798d-3 + 2.0d0*2.763d-3*t9 - 3.0d0*2.341d-4*t92
      end if

      term    = 4.16d+13 * t9i23 * exp(-25.631 * t9i13 * aa)
      dtermdt = -twoth*term*t9i + term*25.631*t9i13*(oneth*t9i*aa - daa)


! the rates
      rev      = 0.5825d0 * exp(-22.224*t9i)
      drevdt   = rev*22.224*t9i2

      fr    = den * rev * term
      dfrdt = den * (drevdt * term + rev * dtermdt) * 1.0d-9
      dfrdd = rev * term

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end





      subroutine rate_p31pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! p31(p,g)s32
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 1.928d-1*z - 1.540d-2*z2 + 6.444d-4*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 1.928d-1 - 2.0d0*1.540d-2*t9 + 3.0d0*6.444d-4*t92
      end if

      term    = 1.08d+16 * t9i23 * exp(-27.042 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 27.042*t9i13*(oneth*t9i*aa - daa))


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.764d+10 * t932 * exp(-102.865*t9i)
      drevdt   = rev*(1.5d0*t9i + 102.865*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_s32ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! s32(a,g)ar36
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 4.913d-2*z + 4.637d-3*z2 - 4.067d-4*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 4.913d-2 + 2.0d0*4.637d-3*t9 - 3.0d0*4.067d-4*t92
      end if

      term    = 1.16d+24 * t9i23 * exp(-66.690 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 66.690*t9i13*(oneth*t9i*aa - daa))


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 6.616d+10 * t932 * exp(-77.080*t9i)
      drevdt   = rev*(1.5d0*t9i + 77.080*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_s32ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! s32(a,p)cl35
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 1.041d-1*z - 1.368d-2*z2 + 6.969d-4*z3
      if (z .eq. 10) then
       daa = 0.0d0
      else
       daa   = 1.041d-1 - 2.0d0*1.368d-2*t9 + 3.0d0*6.969d-4*t92
      end if

      term    = 1.27d+16 * t9i23 * exp(-31.044 * t9i13 * aa)
      dtermdt = -twoth*term*t9i + term*31.044*t9i13*(oneth*t9i*aa - daa)


! the rates
      rev      = 1.144 * exp(-21.643*t9i)
      drevdt   = rev*21.643*t9i2

      fr    = den * rev * term
      dfrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      dfrdd = rev * term

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end





      subroutine rate_cl35pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt


! cl35(p,g)ar36
      aa    = 1.0d0 + 1.761d-1*t9 - 1.322d-2*t92 + 5.245d-4*t93
      daa   = 1.761d-1 - 2.0d0*1.322d-2*t9 + 3.0d0*5.245d-4*t92


      term    =  4.48d+16 * t9i23 * exp(-29.483 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 29.483*t9i13*(oneth*t9i*aa - daa))


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 7.568d+10*t932*exp(-98.722*t9i)
      drevdt   = rev*(1.5d0*t9i + 98.722*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_ar36ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! ar36(a,g)ca40
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 1.458d-1*z - 1.069d-2*z2 + 3.790d-4*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 1.458d-1 - 2.0d0*1.069d-2*t9 + 3.0d0*3.790d-4*t92
      end if

      term    = 2.81d+30 * t9i23 * exp(-78.271 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 78.271*t9i13*(oneth*t9i*aa - daa))


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 6.740d+10 * t932 * exp(-81.711*t9i)
      drevdt   = rev*(1.5d0*t9i + 81.711*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_ar36ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! ar36(a,p)k39
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 4.826d-3*z - 5.534d-3*z2 + 4.021d-4*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 4.826d-3 - 2.0d0*5.534d-3*t9 + 3.0d0*4.021d-4*t92
      end if

      term    = 2.76d+13 * t9i23 * exp(-34.922 * t9i13 * aa)
      dtermdt = -twoth*term*t9i + term*34.922*t9i13*(oneth*t9i*aa - daa)


! the rates
      rev      = 1.128*exp(-14.959*t9i)
      drevdt   = rev*14.959*t9i2

      fr    = den * rev * term
      dfrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      dfrdd = rev * term

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end





      subroutine rate_k39pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! k39(p,g)ca40
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 1.622d-1*z - 1.119d-2*z2 + 3.910d-4*z3
      if (z .eq. 10) then
       daa = 0.0d0
      else
       daa   = 1.622d-1 - 2.0d0*1.119d-2*t9 + 3.0d0*3.910d-4*t92
      end if

      term    = 4.09d+16 * t9i23 * exp(-31.727 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 31.727*t9i13*(oneth*t9i*aa - daa))


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 7.600d+10 * t932 * exp(-96.657*t9i)
      drevdt   = rev*(1.5d0*t9i + 96.657*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_ca40ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! ca40(a,g)ti44
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 1.650d-2*z + 5.973d-3*z2 - 3.889d-04*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 1.650d-2 + 2.0d0*5.973d-3*z - 3.0d0*3.889d-4*z2
      end if

      term    = 4.66d+24 * t9i23 * exp(-76.435 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 76.435*t9i13*(oneth*t9i*aa - daa))


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 6.843d+10 * t932 * exp(-59.510*t9i)
      drevdt   = rev*(1.5d0*t9i + 59.510*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_ca40ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! ca40(a,p)sc43
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 - 1.206d-2*z + 7.753d-3*z2 - 5.071d-4*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = -1.206d-2 + 2.0d0*7.753d-3*t9 - 3.0d0*5.071d-4*t92
      end if

      term    = 4.54d+14 * t9i23 * exp(-32.177 * t9i13 * aa)
      dtermdt = -twoth*term*t9i + term*32.177*t9i13*(oneth*t9i*aa - daa)


! the rates
      rev      = 2.229 * exp(-40.966*t9i)
      drevdt   = rev*40.966*t9i2

      fr    = den * rev * term
      dfrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      dfrdd = rev * term

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end





      subroutine rate_sc43pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! sc43(p,g)ca40
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 1.023d-1*z - 2.242d-3*z2 - 5.463d-5*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 1.023d-1 - 2.0d0*2.242d-3*t9 - 3.0d0*5.463d-5*t92
      end if

      term    = 3.85d+16 * t9i23 * exp(-33.234 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 33.234*t9i13*(oneth*t9i*aa - daa))


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.525d+11 * t932 * exp(-100.475*t9i)
      drevdt   = rev*(1.5d0*t9i + 100.475*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_ti44ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! ti44(a,g)cr48
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 1.066d-1*z - 1.102d-2*z2 + 5.324d-4*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 1.066d-1 - 2.0d0*1.102d-2*t9 + 3.0d0*5.324d-4*t92
      end if

      term    = 1.37d+26 * t9i23 * exp(-81.227 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 81.227*t9i13*(oneth*t9i*aa - daa))


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 6.928d+10*t932*exp(-89.289*t9i)
      drevdt   = rev*(1.5d0*t9i + 89.289*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_ti44ap(temp,den, &
                            fr,dfrdt,dfrdd, &
                            rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                       fr,dfrdt,dfrdd, &
                       rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! ti44(a,p)v47
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 2.655d-2*z - 3.947d-3*z2 + 2.522d-4*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 2.655d-2 - 2.0d0*3.947d-3*t9 + 3.0d0*2.522d-4*t92
      end if

      term    = 6.54d+20 * t9i23 * exp(-66.678 * t9i13 * aa)
      dtermdt = -twoth*term*t9i + term*66.678*t9i13*(oneth*t9i*aa - daa)


! the rates
      rev      = 1.104 * exp(-4.723*t9i)
      drevdt   = rev*4.723*t9i2

      fr    = den * rev * term
      dfrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      dfrdd = rev * term

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end





      subroutine rate_v47pg(temp,den, &
                            fr,dfrdt,dfrdd, &
                            rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                       fr,dfrdt,dfrdd, &
                       rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! v47(p,g)cr48
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 9.979d-2*z - 2.269d-3*z2 - 6.662d-5*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 9.979d-2 - 2.0d0*2.269d-3*t9 - 3.0d0*6.662d-5*t92
      end if

      term    = 2.05d+17 * t9i23 * exp(-35.568 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 35.568*t9i13*(oneth*t9i*aa - daa))


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 7.649d+10*t932*exp(-93.999*t9i)
      drevdt   = rev*(1.5d0*t9i + 93.999*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_cr48ag(temp,den, &
                            fr,dfrdt,dfrdd, &
                            rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                       fr,dfrdt,dfrdd, &
                       rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! cr48(a,g)fe52
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 6.325d-2*z - 5.671d-3*z2 + 2.848d-4*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 6.325d-2 - 2.0d0*5.671d-3*t9 + 3.0d0*2.848d-4*t92
      end if

      term    = 1.04d+23 * t9i23 * exp(-81.420 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 81.420*t9i13*(oneth*t9i*aa - daa))


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 7.001d+10 * t932 * exp(-92.177*t9i)
      drevdt   = rev*(1.5d0*t9i + 92.177*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_cr48ap(temp,den, &
                            fr,dfrdt,dfrdd, &
                            rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                       fr,dfrdt,dfrdd, &
                       rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! cr48(a,p)mn51
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 1.384d-2*z + 1.081d-3*z2 - 5.933d-5*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 1.384d-2 + 2.0d0*1.081d-3*t9 - 3.0d0*5.933d-5*t92
      end if

      term    = 1.83d+26 * t9i23 * exp(-86.741 * t9i13 * aa)
      dtermdt = -twoth*term*t9i + term*86.741*t9i13*(oneth*t9i*aa - daa)


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 0.6087*exp(-6.510*t9i)
      drevdt   = rev*6.510*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_mn51pg(temp,den, &
                            fr,dfrdt,dfrdd, &
                            rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                       fr,dfrdt,dfrdd, &
                       rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! mn51(p,g)fe52
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 8.922d-2*z - 1.256d-3*z2 - 9.453d-5*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 8.922d-2 - 2.0d0*1.256d-3*t9 - 3.0d0*9.453d-5*t92
      end if

      term    = 3.77d+17 * t9i23 * exp(-37.516 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 37.516*t9i13*(oneth*t9i*aa - daa))


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.150d+11*t932*exp(-85.667*t9i)
      drevdt   = rev*(1.5d0*t9i + 85.667*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_fe52ag(temp,den, &
                            fr,dfrdt,dfrdd, &
                            rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                       fr,dfrdt,dfrdd, &
                       rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! fe52(a,g)ni56
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 7.846d-2*z - 7.430d-3*z2 + 3.723d-4*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 7.846d-2 - 2.0d0*7.430d-3*t9 + 3.0d0*3.723d-4*t92
      end if

      term    = 1.05d+27 * t9i23 * exp(-91.674 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 91.674*t9i13*(oneth*t9i*aa - daa))


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 7.064d+10*t932*exp(-92.850*t9i)
      drevdt   = rev*(1.5d0*t9i + 92.850*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_fe52ap(temp,den, &
                            fr,dfrdt,dfrdd, &
                            rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                       fr,dfrdt,dfrdd, &
                       rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! fe52(a,p)co55
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 1.367d-2*z + 7.428d-4*z2 - 3.050d-5*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 1.367d-2 + 2.0d0*7.428d-4*t9 - 3.0d0*3.050d-5*t92
      end if

      term    = 1.30d+27 * t9i23 * exp(-91.674 * t9i13 * aa)
      dtermdt = -twoth*term*t9i + term*91.674*t9i13*(oneth*t9i*aa - daa)


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 0.4597*exp(-9.470*t9i)
      drevdt   = rev*9.470*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_co55pg(temp,den, &
                            fr,dfrdt,dfrdd, &
                            rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                       fr,dfrdt,dfrdd, &
                       rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! co55(p,g)ni56
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 9.894d-2*z - 3.131d-3*z2 - 2.160d-5*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 9.894d-2 - 2.0d0*3.131d-3*t9 - 3.0d0*2.160d-5*t92
      end if

      term    = 1.21d+18 * t9i23 * exp(-39.604 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 39.604*t9i13*(oneth*t9i*aa - daa))


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.537d+11*t932*exp(-83.382*t9i)
      drevdt   = rev*(1.5d0*t9i + 83.382*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_fe52ng(temp,den, &
                            fr,dfrdt,dfrdd, &
                            rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                       fr,dfrdt,dfrdd, &
                       rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,tq2


! fe52(n,g)fe53
      tq2     = t9 - 0.348d0
      term    = 9.604d+05 * exp(-0.0626*tq2)
      dtermdt = -term*0.0626

! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 2.43d+09 * t932 * exp(-123.951*t9i)
      drevdt   = rev*(1.5d0*t9i + 123.951*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_fe53ng(temp,den, &
                            fr,dfrdt,dfrdd, &
                            rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                       fr,dfrdt,dfrdd, &
                       rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,tq1,tq10,dtq10,tq2


! fe53(n,g)fe54
      tq1   = t9/0.348
      tq10  = tq1**0.10
      dtq10 = 0.1d0*tq10/(0.348*tq1)
      tq2   = t9 - 0.348d0

      term    = 1.817d+06 * tq10 * exp(-0.06319*tq2)
      dtermdt = term/tq10*dtq10 - term*0.06319

! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.56d+11 * t932 * exp(-155.284*t9i)
      drevdt   = rev*(1.5d0*t9i + 155.284*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end



      subroutine rate_fe54ng(temp,den, &
                            fr,dfrdt,dfrdd, &
                            rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                       fr,dfrdt,dfrdd, &
                       rr,drrdt,drrdd

! locals
      double precision aa,daa,bb,dbb,term,dtermdt


! fe54(n,g)fe55
      aa   =  2.307390d+01 - 7.931795d-02 * t9i + 7.535681d+00 * t9i13 &
            - 1.595025d+01 * t913 + 1.377715d+00 * t9 - 1.291479d-01 * t953 &
            + 6.707473d+00 * log(t9)

      daa  =  7.931795d-02 * t9i2 - oneth * 7.535681d+00 * t9i43 &
            - oneth * 1.595025d+01 *t9i23 + 1.377715d+00 - fiveth * 1.291479d-01 *t923 &
            + 6.707473d+00 * t9i

      if (aa .lt. 200.0) then
       term    = exp(aa)
       dtermdt = term*daa*1.0d-9
      else
       term    = exp(200.0d0)
       dtermdt = 0.0d0
      end if

      bb  = 4.800293d+09 * t932 * exp(-1.078986d+02 * t9i)
      dbb = bb*(1.5d0*t9i + 1.078986d+02 * t9i2)

! reverse rate
      rr    = term*bb
      drrdt = dtermdt*bb + term*dbb*1.0d-9
      drrdd = 0.0d0

! forward rate
      dfrdd = term
      fr    = term*den
      dfrdt = dtermdt*den

      return
      end



      subroutine rate_fe54pg(temp,den, &
                            fr,dfrdt,dfrdd, &
                            rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                       fr,dfrdt,dfrdd, &
                       rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,z,z2,z3


! fe54(p,g)co55
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 9.593d-2*z - 3.445d-3*z2 + 8.594d-5*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 9.593d-2 - 2.0d0*3.445d-3*t9 + 3.0d0*8.594d-5*t92
      end if

      term    = 4.51d+17 * t9i23 * exp(-38.483 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 38.483*t9i13*(oneth*t9i*aa - daa))


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 2.400d+09 * t932 * exp(-58.605*t9i)
      drevdt   = rev*(1.5d0*t9i + 58.605*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end



      subroutine rate_fe54ap(temp,den, &
                            fr,dfrdt,dfrdd, &
                            rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                       fr,dfrdt,dfrdd, &
                       rr,drrdt,drrdd

! locals
      double precision aa,daa,bb,dbb,term,dtermdt


! fe54(a,p)co57
      aa   =  3.97474900d+01 - 6.06543100d+00 * t9i + 1.63239600d+02 * t9i13 &
            - 2.20457700d+02 * t913 + 8.63980400d+00 * t9 - 3.45841300d-01 * t953 &
            + 1.31464200d+02 * log(t9)

      daa  =  6.06543100d+00 * t9i2 - oneth * 1.63239600d+02 * t9i43 &
            - oneth * 2.20457700d+02 * t9i23 + 8.63980400d+00 - fiveth * 3.45841300d-01 * t923 &
            + 1.31464200d+02  * t9i

      if (aa .lt. 200.0) then
       term    = exp(aa)
       dtermdt = term*daa*1.0d-9
      else
       term    = exp(200.0d0)
       dtermdt = 0.0d0
      end if

      bb  = 2.16896000d+00  * exp(-2.05631700d+01 * t9i)
      dbb = bb * 2.05631700d+01 * t9i2

! reverse rate
      drrdd = term
      rr    = term*den
      drrdt = dtermdt*den

! forward rate
      fr    = rr*bb
      dfrdt = drrdt*bb + rr*dbb*1.0d-9
      dfrdd = drrdd*bb

      return
      end







      subroutine rate_fe55ng(temp,den, &
                            fr,dfrdt,dfrdd, &
                            rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                       fr,dfrdt,dfrdd, &
                       rr,drrdt,drrdd

! locals
      double precision aa,daa,bb,dbb,term,dtermdt


! fe55(n,g)fe56
      aa   =  1.954115d+01 - 6.834029d-02 * t9i + 5.379859d+00 * t9i13 &
            - 8.758150d+00 * t913 + 5.285107d-01 * t9 - 4.973739d-02  * t953 &
            + 4.065564d+00  * log(t9)

      daa  =  6.834029d-02 * t9i2 - oneth * 5.379859d+00 * t9i43 &
            - oneth * 8.758150d+00 * t9i23 + 5.285107d-01 - fiveth * 4.973739d-02  *t923 &
            + 4.065564d+00  * t9i

      if (aa .lt. 200.0) then
       term    = exp(aa)
       dtermdt = term*daa*1.0d-9
      else
       term    = exp(200.0d0)
       dtermdt = 0.0d0
      end if

      bb  = 7.684279d+10  * t932 * exp(-1.299472d+02  * t9i)
      dbb = bb*(1.5d0*t9i + 1.299472d+02 * t9i2)

! reverse rate
      rr    = term*bb
      drrdt = dtermdt*bb + term*dbb*1.0d-9
      drrdd = 0.0d0

! forward rate
      dfrdd = term
      fr    = term*den
      dfrdt = dtermdt*den

      return
      end




      subroutine rate_fe56pg(temp,den, &
                            fr,dfrdt,dfrdd, &
                            rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                       fr,dfrdt,dfrdd, &
                       rr,drrdt,drrdd

! locals
      double precision aa,daa,bb,dbb,term,dtermdt


! fe56(p,g)co57
   
      aa   =  1.755960d+02 - 7.018872d+00 * t9i + 2.800131d+02 * t9i13 &
            - 4.749343d+02 * t913 + 2.683860d+01 * t9 - 1.542324d+00  * t953 &
            + 2.315911d+02  * log(t9)

      daa  =  7.018872d+00 * t9i2 - oneth * 2.800131d+02 * t9i43 &
            - oneth * 4.749343d+02 * t9i23 + 2.683860d+01 - fiveth * 1.542324d+00  *t923 &
            + 2.315911d+02  * t9i

      if (aa .lt. 200.0) then
       term    = exp(aa)
       dtermdt = term*daa*1.0d-9
      else
       term    = exp(200.0d0)
       dtermdt = 0.0d0
      end if

      bb  = 2.402486d+09 * t932 * exp(-6.995192d+01 * t9i)
      dbb = bb*(1.5d0*t9i + 6.995192d+01 * t9i2)


! reverse rate
      rr    = term*bb
      drrdt = dtermdt*bb + term*dbb*1.0d-9
      drrdd = 0.0d0

! forward rate
      dfrdd = term
      fr    = term*den
      dfrdt = dtermdt*den

      return
      end
!---------------------------------------------------------------------







!---------------------------------------------------------------------
      subroutine net_initialize
      include 'implno.dek'
      include 'network.dek'

! initializes quantities

! local variables
      integer   i


! general options
      screen_on      = 1
      use_tables     = 1
      weak_on        = 1
      lmp_on         = 1
      pure_network   = 0

! rate controls
      write_rates    = 0
      shut_down_ng   = 0
      shut_down_pn   = 0
      shut_down_pg   = 0
      shut_down_ap   = 0
      shut_down_an   = 0
      shut_down_ag   = 0
      shut_down_specials = 0
      renorm_abundances = 1


! printing information
      iprint_files  = 1
!      iprint_screen = 1
!      iprint_files  = 0
      iprint_screen = 0


! inititailize the burn type logicals
      one_step             = .false.
      hydrostatic          = .false.
      expansion            = .false.
      self_heat_const_den  = .false.
      self_heat_const_pres = .false.
      pt_hist              = .false.
      bbang                = .false.
      detonation           = .false.
      trho_hist            = .false.


! adiabatic expanion off
      psi       = 0.0d0
      temp_stop = 1.0d30


! scaling values for the energy, temperature, and density evolution equations
      energy_scale     = 1.0d0
      inv_energy_scale = 1.0d0
      temp_scale       = 1.0d0
      inv_temp_scale   = 1.0d0
      den_scale        = 1.0d0
      inv_den_scale    = 1.0d0

      one_step = .true.

      return
      end
!---------------------------------------------------------------------







!---------------------------------------------------------------------
      subroutine net_summary(tstep,tin,din,ein,tout,dout,eout,conserv,nbad,nok,xout)
      include 'implno.dek'
      include 'timers.dek'
      include 'vector_eos.dek'
      include 'burn_common.dek'
      include 'network.dek'


! writes out a summary of the network run

! declare the pass
      integer :: nbad,nok
      real*8  :: tstep,tin,din,ein,tout,dout,eout,conserv,xout(*)


! local variables
      character*80      :: summary
      integer           :: i,j,k,lenstr,ioff
      real*8            :: ye
      real*8, parameter :: abundance_print_threshold = 1.0d30


! popular format statements
 01   format(a,'summary.dat')
 02   format(1x,a,'=',1pe10.3,' ',a,'=',1pe10.3,' ', &
                a,'=',1pe10.3,' ',a,'=',1pe10.3,' ', &
                a,'=',1pe10.3)
 03   format(1x,a,1pe20.12)
 04   format(1x,a,':',/, &
             1x,3(a,1pe20.12),/, &
             1x,3(a,1pe20.12),/, &
             1x,1(a,1pe11.3),2(a,i5))
 08   format(1x,a,1pe10.3,a)
 09   format(1x,a,i2,a)



! construct the file name and open it
       write(summary,01) hfile(1:lenstr(hfile,80))
       call sqeeze(summary)
       open(unit=41,file=summary,status='unknown')


       write(6,*) ' '
       write(6,04) trim(netname), &
                   ' tin =',tin,' din =',din,' ein =',ein, &
                   ' tout=',tout,' dout=',dout,' eout=',eout, &
                   ' sum =',conserv,' nbad=',nbad,' nok=',nok
       write(6,*) ' '

       write(41,*) ' '
       write(41,04) trim(netname), &
                   ' tin =',tin,' din =',din,' ein =',ein, &
                   ' tout=',tout,' dout=',dout,' eout=',eout, &
                   ' sum =',conserv,' nbad=',nbad,' nok=',nok
       write(41,*) ' '



! write out the biggest mass fractions
       xout(ionbeg:ionend) = max(xout(ionbeg:ionend),1.0d-30)
       call indexx(ionmax,xout(ionbeg),izwork1(ionbeg))
       ioff = ionbeg - 1
       
       if (abundance_print_threshold .le. 1.0  .and. abundance_print_threshold .gt. 0.0) then
        do i=ionmax,1,-1
         if (xout(izwork1(i)+ioff) .lt. abundance_print_threshold) then
          k = i + 1
          write(6,08)  'mass fractions larger than ',abundance_print_threshold
          write(41,08) 'mass fractions larger than ',abundance_print_threshold
          goto 20
         end if
        end do
       else
        j = min(20,ionmax)
        k = max(ionmax-19,ionbeg)
        write(6,09)  'top ',j,' mass fractions:'
        write(41,09) 'top ',j,' mass fractions:'
       end if

 20   continue


      ye  = sum(zion(1:ionmax)*xout(1:ionmax)/aion(1:ionmax))

       write(6,02) (ionam(izwork1(i)+ioff),xout(izwork1(i)+ioff), i=ionmax,k,-1)
       if (iprot .ne. 0 .and. ineut .ne. 0) then
        write(6,02) ionam(iprot),xout(iprot), &
                    ionam(ineut),xout(ineut), &
                    ionam(ihe4),xout(ihe4), &
                     'ye   ',ye
       end if
       write(6,*) ' '

       write(41,02) (ionam(izwork1(i)+ioff),xout(izwork1(i)+ioff), i=ionmax,k,-1)
       if (iprot .ne. 0 .and. ineut .ne. 0) then
        write(41,02) ionam(iprot),xout(iprot), &
                     ionam(ineut),xout(ineut), &
                     ionam(ihe4),xout(ihe4), &
                     'ye   ',ye
       end if
       write(41,*) ' '



! end the clock
      call zsecond(timtot)
      timtot = timtot - timzer
      call timlap(timtot,hours,minuts,secs,msecs)
      write(6,100) hours,minuts,secs,msecs
      write(41,100) hours,minuts,secs,msecs
 100  format(1x,'cpu time : ',i2.2,' hrs  ',i2.2,' min  ', &
                              i2.2,' sec  ',i6,' usec',/,/)


! close up shop
      close(unit=41)
      return
      end
!---------------------------------------------------------------------







!---------------------------------------------------------------------
! this file contains routines that sort, search and select parts of arrays:
!
! index and rank makers:
! routine indexx constructs a sort index for a real*8 array
! routine indexx constructs a sort index for a real*16 array
! routine locate find the array position corresponding to a given value


      subroutine indexx(n,arr,indx)
      include 'implno.dek'
!
! indexes an array arr(1:n). that is it outputs the array indx(1:n) such
! that arr(indx(j)) is in ascending order for j=1...n. the input quantities
! are not changed.
!
! declare
      integer          n,indx(n),m,nstack
      parameter        (m=7, nstack = 50)
      integer          i,indxt,ir,itemp,j,jstack,k,l,istack(nstack)
      double precision arr(n),a
!
! initialize
      do 11 j=1,n
       indx(j) = j
11    continue
      jstack = 0
      l      = 1
      ir     = n
!
! insertion sort when subbarray small enough
1     if (ir - l .lt. m) then
       do 13 j=l+1,ir
        indxt = indx(j)
        a     = arr(indxt)
        do 12 i=j-1,l,-1
         if (arr(indx(i)) .le. a) go to 2
         indx(i+1) = indx(i)
12      continue
        i = l - 1
2       indx(i+1) = indxt
13     continue
!
! pop stack and begin a new round of partitioning
       if (jstack .eq. 0) return
       ir     = istack(jstack)
       l      = istack(jstack-1)
       jstack = jstack - 2
!
! choose median of left, center and right elements as partitioning element
! also rearrange so that a(l+1) < a(l) < a(ir)
      else
       k         = (l + ir)/2
       itemp     = indx(k)
       indx(k)   = indx(l+1)
       indx(l+1) = itemp

       if (arr(indx(l)) .gt. arr(indx(ir))) then
        itemp    = indx(l)
        indx(l)  = indx(ir)
        indx(ir) = itemp
       end if


       if(arr(indx(l+1)).gt.arr(indx(ir)))then
        itemp=indx(l+1)
        indx(l+1)=indx(ir)
        indx(ir)=itemp
       endif
       if(arr(indx(l)).gt.arr(indx(l+1)))then
        itemp=indx(l)
        indx(l)=indx(l+1)
        indx(l+1)=itemp
       endif

!
! initialize pointers for partitioning
       i     = l + 1
       j     = ir
       indxt = indx(l+1)
       a     = arr(indxt)
3      continue
       i = i + 1
       if (arr(indx(i)) .lt. a) go to 3
4      continue
       j = j - 1
       if (arr(indx(j)) .gt. a) go to 4
       if (j .lt. i) go to 5
       itemp   = indx(i)
       indx(i) = indx(j)
       indx(j) = itemp
       go to 3
!
5      indx(l+1) = indx(j)
       indx(j)   = indxt
       jstack    = jstack + 2
!
! push pointers to larger subarray on stack
       if (jstack .gt. nstack) stop 'jstack > nstack in routine indexx'
       if (ir - i + 1  .ge.  j - l) then
        istack(jstack)   = ir
        istack(jstack-1) = i
        ir               = j - 1
       else
        istack(jstack)   = j-1
        istack(jstack-1) = l
        l                = i
       end if
      end if
      go to 1
      end






      subroutine indexx16(n,arr,indx)
      include 'implno.dek'

! indexes a real*16 array arr(1:n). outputs the array indx(1:n) such
! that arr(indx(j)) is in ascending order for j=1...n.
! input quantities are not changed.
!
! declare
      integer      n,indx(n),m,nstack
      parameter    (m=7, nstack = 50)
      integer      i,indxt,ir,itemp,j,jstack,k,l,istack(nstack)
      real*16      arr(n),a


! initialize
      do j=1,n
       indx(j) = j
      enddo
      jstack = 0
      l      = 1
      ir     = n

! insertion sort when subbarray small enough
1     if (ir - l .lt. m) then
       do j=l+1,ir
        indxt = indx(j)
        a     = arr(indxt)
        do i=j-1,l,-1
         if (arr(indx(i)) .le. a) goto 2
         indx(i+1) = indx(i)
        enddo
        i = l - 1
2       indx(i+1) = indxt
       enddo

! pop stack and begin a new round of partitioning
       if (jstack .eq. 0) return
       ir     = istack(jstack)
       l      = istack(jstack-1)
       jstack = jstack - 2

! choose median of left, center and right elements as partitioning element
! also rearrange so that a(l+1) < a(l) < a(ir)
      else
       k         = (l + ir)/2
       itemp     = indx(k)
       indx(k)   = indx(l+1)
       indx(l+1) = itemp

       if (arr(indx(l)) .gt. arr(indx(ir))) then
        itemp    = indx(l)
        indx(l)  = indx(ir)
        indx(ir) = itemp
       end if


       if(arr(indx(l+1)).gt.arr(indx(ir)))then
        itemp=indx(l+1)
        indx(l+1)=indx(ir)
        indx(ir)=itemp
       endif
       if(arr(indx(l)).gt.arr(indx(l+1)))then
        itemp=indx(l)
        indx(l)=indx(l+1)
        indx(l+1)=itemp
       endif


! initialize pointers for partitioning
       i     = l + 1
       j     = ir
       indxt = indx(l+1)
       a     = arr(indxt)
3      continue
       i = i + 1
       if (arr(indx(i)) .lt. a) goto 3
4      continue
       j = j - 1
       if (arr(indx(j)) .gt. a) goto 4
       if (j .lt. i) goto 5
       itemp   = indx(i)
       indx(i) = indx(j)
       indx(j) = itemp
       goto 3

5      indx(l+1) = indx(j)
       indx(j)   = indxt
       jstack    = jstack + 2

! push pointers to larger subarray on stack
       if (jstack .gt. nstack) stop 'jstack > nstack in routine indexx'
       if (ir - i + 1  .ge.  j - l) then
        istack(jstack)   = ir
        istack(jstack-1) = i
        ir               = j - 1
       else
        istack(jstack)   = j-1
        istack(jstack-1) = l
        l                = i
       end if
      end if
      goto 1
      end






      subroutine locate(xx,n,x,j)
      implicit none
      save


! given an array xx of length n, and a value of x, this routine returns
! a value j such that x is between xx(j) and xx(j+1). the array xx must be
! monotonic. j=0 or j=n indicates that x is out of range. bisection is used
! to find the entry

! declare
      integer           n,j,jl,ju,jm
      double precision  xx(n),x

! initialize
      jl = 0
      ju = n+1

! compute a midpoint, and replace either the upper or lower limit
 10   if (ju-jl .gt. 1) then
       jm = (ju+jl)/2
       if ( (xx(n) .ge. xx(1)) .eqv. (x .ge. xx(jm)) ) then
        jl = jm
       else
        ju = jm
       end if
       goto 10
      end if
      if (x .eq. xx(1))then
        j = 1
      else if(x .eq. xx(n))then
        j = n - 1
      else
        j = jl
      end if
      return
      end
!---------------------------------------------------------------------







!---------------------------------------------------------------------
      subroutine cjsolve(kkase,xmass_up,temp_up,den_up,mach, &
                         qx,xmass_det,ener_up,pres_up,cs_up, &
                         vel_det,vel_x,temp_x,den_x,ener_x,pres_x,cs_x)

      include 'implno.dek'
      include 'const.dek'
      include 'vector_eos.dek'
      include 'network.dek'

! solves the hugoniot and rayleigh relations for a detonation or a shock.
! an nse distribution is assumed for the detonated material.

! input:
! kkase     = 1 for a chapman-jouget detonation
!           = 2 for a strong point driven detonation
!           = 3 for a weak point driven detonation
!           = 4 for a shock wave
! xmass_up  = upstream composition vector
! temp_up   = temperature of upstream material
! den_up    = density of upstream material
! mach      = mach number of shock or detonation


! output:
! qx        = energy/gram from burning
! xmass_det = burned composition
! ener_up   = energy of upstream material
! pres_up   = pressure of upstream material
! cs_up     = sound speed of upstream material
! vel_det   = speed of detonation or shock front = upstream speed
! vel_x     = speed behind detonation or shock = downstream speed
! temp_x    = temperature of downstream material
! den_x     = density of downstream material
! ener_x    = internal energy of downstream material
! pres_x    = pressure of downstream material
! cs_x      = sound speed of downstream material


! declare the pass
      integer          kkase
      double precision xmass_up(1),temp_up,den_up,mach,qx,xmass_det(1), &
                       ener_up,pres_up,cs_up,vel_det,vel_x,temp_x,den_x, &
                       ener_x,pres_x,cs_x


! common block communication with the routine cjfunc
      integer          kase
      double precision xmup(abignet),den1,temp1,p1,u1,v1,ye1,vel1,cs1, &
                       mach1,qburn,den2,temp2,p2,u2,v2,vel2,cs2
      common /cjstate/ xmup,den1,temp1,p1,u1,v1,ye1,vel1,cs1, &
                       mach1,qburn,den2,temp2,p2,u2,v2,vel2,cs2, &
                       kase


! locals
      external         cjfunc
      character*8      type
      logical          check
      integer          i,n,ntrial,ntaken,nfev,nstrong,nsmax,nweak,nwmax
      parameter        (ntrial = 60, n=2, nsmax = 10, nwmax=10)
      double precision cv1,g1,g1p1,g1m1,msq,v1mv2,x(n),f(n),dum, &
                       xl,xx,gsq,den2_sav,abar,zbar,wbar,ye,xcess, &
                       tolx,tolf,conv
      parameter        (tolf = 1.0d-6, tolx = 1.0d-6)
      parameter        (conv = ev2erg * 1.0d6 * avo)



! check the input
      if (kkase .lt. 1  .or. kkase .gt. 4) then
       write(6,*) 'kkase =',kkase
       stop 'kkase in cjsolve is invalid'
      end if


! transfer the input to common
      kase    = kkase
      temp1   = temp_up
      den1    = den_up
      mach1   = mach
      do i=1,ionmax
       xmup(i) = xmass_up(i)
      enddo

! load the eos conditions
      call azbarapprox13(xmass_up(ionbeg),aion(ionbeg),zion(ionbeg),wion(ionbeg),ionmax, &
           ymass(ionbeg),abar,zbar,wbar,ye,xcess)

      temp_row(1) = temp_up
      den_row(1)  = den_up
      abar_row(1) = abar
      zbar_row(1) = zbar
      jlo_eos = 1
      jhi_eos = 1


! call an eos
!      call helmeos
      call helmeosaprox13


! set upstream thermodynamic conditions
      p1      = ptot_row(1)
      u1      = etot_row(1)
      cs1     = cs_row(1)
      cv1     = cv_row(1)
      g1      = gam1_row(1)
      v1      = 1.0d0/den1
      ye1     = zbar/abar
      pres_up = p1
      ener_up = u1
      cs_up   = cs1



! now we start making initial guesses for downstream temperature and density

! for the detonation cases here is a guess for the energy generated
! a pure si28 burned composition seems robust

      if (kase .eq. 1  .or. kase .eq. 2  .or. kase .eq. 3) then
       qburn  = (1.0d0 - xmass_up(isi28))/aion(isi28)*bion(isi28)*conv
      end if



! an initial guess for the density and temperature of the burned material
! from landau & lifshitz fluid dynamics, 129.15
! with gamma2 about gamma1 and cv2 about 1/8 cv1
! this is nearly exact for the ions

      if (kase .eq. 1  .or. kase .eq. 2  .or. kase .eq. 3) then
       g1p1  = g1 + 1.0d0
       den2  = den1 * g1p1/g1
       temp2 = 2.0d0 * g1 * qburn / (8.0d0 * cv1 *g1p1)


! limit temp2 so nse material is not photodisintegrated back to helium
       temp2 = min(5.0d9, max(temp1,temp2))


! modify density guess for strong or weak detonations
! the weak detonation appears to be a much stronger attractor

        if (kase .eq. 2) then
         den2 = 1.2d0 * den2
        else if (kase .eq. 3) then
         den2 = 0.9 * den2
        end if


! an initial guess for the density and temperature behind a shock wave
! from landau & lifshitz fluid dynamics, eq 85.7 to 85.10
! its exact for the ions

      else
       g1p1  = g1 + 1.0d0
       g1m1  = g1 - 1.0d0
       msq   = mach1*mach1
       den2  = den1  * g1p1*msq/(g1m1*msq + 2.0d0)
       temp2 = temp1 * (2.0d0*g1*msq - g1m1) * &
                       (g1m1*msq + 2.0d0) / (g1p1*g1p1*msq)
      end if


! loop to here with a new den2 if the strong or weak point is not proper
       nstrong  = 0
       nweak    = 0
       den2_sav = den2
100    continue


! done with the initial guesses section




! root find of the rayleigh line and the hugoniot
      x(1)   = den2
      x(2)   = temp2

      call xnewt_cj(ntrial,x,n,tolx,tolf,ntaken,check,nfev,cjfunc)


      if (check .or. ntaken .eq. ntrial) then
       write(6,*)
       write(6,*) 'check convergence of xnewt_cj root find'
       write(6,*)
       write(6,*) 'iterations taken =',ntaken
       write(6,*) 'function evals =',nfev
       write(6,111) 'roots =',x(1),x(2)
 111   format(1x,a,1p2e14.6)
      end if


! with the converged values, get the return arguments
      call cjfunc(dum,x,f)



! set the return arguments
      do i=1,ionmax
       xmass_det(i) = xmass(i)
      enddo
      vel1    = den2/den1 * vel2
      vel_det = vel1
      vel_x   = vel2
      temp_x  = temp2
      den_x   = den2
      pres_x  = p2
      ener_x  = u2
      cs_x    = cs2
      qx      = qburn



! the strong point solution needs to be checked before returning
      if (kase .eq. 2  .and. nstrong .gt. nsmax) then
       write(6,*) ' '
       write(6,*) 'warning: did not find strong point solution'
       write(6,*) ' '
       return
      end if

      if (kase .eq. 2) then
       if (vel_x .ge. cs_x) then
        nstrong  = nstrong + 1
        den2     = 1.1d0 * den2_sav
        den2_sav = den2
        goto 100
       end if
      end if


! the weak point solution needs to be checked before returning
      if (kase .eq. 3  .and. nweak .gt. nwmax) then
       write(6,*) ' '
       write(6,*) 'warning: did not find weak point solution'
       write(6,*) ' '
       return
      end if

      if (kase .eq. 3) then
       if (vel_x .le. cs_x) then
        nweak    = nweak + 1
        den2     = 0.9d0 * den2_sav
        den2_sav = den2
        goto 100
       end if
      end if

! normal bail point
      return
      end







      subroutine cjfunc(x,y,f)
      include 'implno.dek'
      include 'const.dek'
      include 'vector_eos.dek'
      include 'network.dek'


! this routine returns the functions to do the root find on
! input is x (not relevant here) and y, a vector of the unknowns,
! y(1) is the density and y(2) is the temperature.
! output is f, a vector of residuals to be minimized. f(1) is the
! rayleigh line and f(2) is the hugoniot.

! declare the pass
      double precision x,y(*),f(*)


! common block communication with the routine cjfunc
      integer          kase
      double precision xmup(abignet),den1,temp1,p1,u1,v1,ye1,vel1,cs1, &
                       mach1,qburn,den2,temp2,p2,u2,v2,vel2,cs2
      common /cjstate/ xmup,den1,temp1,p1,u1,v1,ye1,vel1,cs1, &
                       mach1,qburn,den2,temp2,p2,u2,v2,vel2,cs2, &
                       kase


! locals
      integer          i,newguess,iprint
      double precision xmnse(abignet),xmunn,xmupp,abar,zbar,wbar, &
                       ye,xcess,conv
      parameter        (conv = ev2erg * 1.0d6 * avo)



! map the input vector, bail if its hosed
      den2  = y(1)
      temp2 = y(2)

      if (den2 .lt. 0.0 .or. temp2 .le. 0.0) then
       f(1) = 1.0d30
       f(2) = 1.0d30
       return
      end if


! load the nse composition
      if (kase .eq. 1  .or. kase .eq. 2  .or. kase .eq. 3) then
       iprint   = 0
       newguess = 1
       call nse(temp2,den2,ye1,newguess,1,1,xmnse,xmunn,xmupp,iprint)
      end if


! get the energy generated
      qburn = 0.0d0
      if (kase .eq. 1  .or. kase .eq. 2  .or. kase .eq. 3) then
       do i=1,ionmax
        qburn  = qburn + (xmnse(i) - xmup(i))/aion(i) *bion(i)
       end do
       qburn = qburn * conv
      end if



! get the eos
! for the detonation cases use the nse composition
! for the shock case, the shocked composition is the upstream composition

      if (kase .eq. 1  .or. kase .eq. 2  .or. kase .eq. 3) then
       do i=1,ionmax
        xmass(i) = xmnse(i)
       enddo
      else
       do i=1,ionmax
        xmass(i) = xmup(i)
       enddo
      end if
      call azbarapprox13(xmass(ionbeg),aion(ionbeg),zion(ionbeg),wion(ionbeg),ionmax, &
           ymass(ionbeg),abar,zbar,wbar,ye,xcess)

      temp_row(1) = temp2
      den_row(1)  = den2
      abar_row(1) = abar
      zbar_row(1) = zbar
      jlo_eos = 1
      jhi_eos = 1

!      call helmeos
      call helmeosaprox13


! set the downstream thermodynamic variables
      v2    = 1.0d0/den_row(1)
      p2    = ptot_row(1)
      u2    = etot_row(1)
      cs2   = cs_row(1)



! for a chapman-jouget detonation vel2 is the burned sound speed
! for a strong or weak detonation or a shock wave the upstream
! mach number mach1 is specified

      if (kase .eq. 1) then
        vel2  = cs2
      else
       vel2   = mach1 * cs1 * den1/den2
      endif


! the rayleigh line, glassman, page 227, eq. 6, fickett & davis page 17
      f(1) = (p2 - p1) - (vel2*vel2*den2*den2) * (v1 - v2)


! the specific internal energy hugoniot, glassman, page 229, eq. 11
      f(2) = u1 + qburn - u2 + 0.5d0 *( (p1+p2) * (v1-v2) )


! scale the functions for better behavior of the root finder
      f(1)  = f(1) * v1/p1
      f(2)  = f(2) / (p1 * v1)


      return
      end






      subroutine xnewt_cj(ntrial,x,n,tolx,tolf,ntaken,check,nfev,func)
      include 'implno.dek'

! given an initial guess x(1:n) for the root of n equations, this routine
! finds the root by a globally convergent newtons method. the vector of
! functions to be zeroed, called fvec(1:n) in the routine below, is
! returned by the user supplied routine func. the output quantity check
! is false on nomal return and true if the routine has converged to a
! local minimum of the function xfminx_cj. if so, try restarting from a
! different initial guess.
!
! np is the maximum number of equations n
! ntrial is the maximum number of iterations to try
! ntaken is the number of iterations done
! tolf sets the convergence on function values
! tolmin sets the criterion for deciding wheather spurious convergence to
!        a false minimum of xfminx_cj has occured
! tolx is the convergence criteria on deltax
! stpmx is the scaled maximum step length allowed in the line searches
! nfev is the number of function evaluations


! declare the pass
      external         func
      logical          check
      integer          ntrial,n,ntaken,nfev
      double precision x(n),tolf,tolx


! locals
      integer          np
      parameter        (np=4)
      integer          nn,i,its,j,indx(np)
      double precision fvec(np),tolmin,stpmx,d,den,f, &
                       fold,stpmax,sum,temp,test,fjac(np,np),g(np), &
                       p(np),xold(np),xfminx_cj,dum
      parameter        (tolmin = 1.0d-12, &
                        stpmx = 2.0d0)


! common block communicates values from routine xfminx_cj
      common /newtcj/  fvec,nn


! initialize
      if (n .gt. np) stop 'n > np in routine xnewt'
      nn   = n
      f    = xfminx_cj(x,func)
      nfev = 1

!  test for the initial guess being a root, using a more stringent tolf
      test = 0.0d0
      do i=1,n
       if (abs(fvec(i)) .gt. test) test = abs(fvec(i))
      enddo
      if (test .lt. 0.01*tolf) then
       check = .false.
       return
      end if


! get stpmax for the line search
      sum = 0.0d0
      do i=1,n
       sum = sum + x(i)*x(i)
      enddo
      stpmax = stpmx * max(sqrt(sum),dfloat(n))


! start of iteration loop; get the jacobian
      do its = 1, ntrial
       ntaken = its

! second order accurate jacobian
       call jac_cj(dum,x,fjac,n,n,np,np,func)
       nfev = nfev + 2*n + 1

! compute grad f for the line searches
       do i=1,n
        sum = 0.0d0
        do j=1,n
         sum = sum + fjac(j,i)*fvec(j)
        enddo
        g(i) = sum
       enddo


! store x, and f and form right hand sides
       do i=1,n
        xold(i) = x(i)
       enddo
       fold = f
       do i=1,n
        p(i) = -fvec(i)
       enddo


! solve the linear systems
       call ludcmp(fjac,n,np,indx,d)
       call lubksb(fjac,n,np,indx,p)


! line search returns new x and f
! it also gets fvec at the new x when it calls xfminx_cj
       call lnsrch_cj(n,xold,fold,g,p,x,f,stpmax,check,nfev,func)


! test for convergence on function value
       test = 0.0d0
       do i=1,n
        if (abs(fvec(i)) .gt. test) test = abs(fvec(i))
       enddo
       if (test .lt. tolf) then
        check = .false.
        return
       end if

! check for zero gradiant, i.e. spurious convergence
       if (check) then
        test = 0.0d0
        den  = max(f, 0.5d0 * n)
        do i=1,n
         temp = abs(g(i)) * max(abs(x(i)),1.0d0)/den
         if (temp .gt. test) test = temp
        enddo
        if (test .lt. tolmin) then
         check = .true.
        else
         check = .false.
        end if
        return
       end if

! test for convergence on deltax
       test = 0.0d0
       do i=1,n
        temp = (abs(x(i)-xold(i)))/max(abs(x(i)),1.0d0)
        if (temp .gt. test) test = temp
       enddo
       if (test .lt. tolx) return

! back for another iteration
      enddo
      check = .true.
      return
      end





      subroutine lnsrch_cj(n,xold,fold,g,p,x,f,stpmax,check,nfev,func)
      include 'implno.dek'

! given an n dimensional point xold(1:n), the value of the function fold
! and the gradient g(1:n) at the point, and a direction p(1:n), this routine
! finds a new point x(1:n) along the direction of p from xold where the
! function xfminx_cj has decreased "sufficiently". the new function value is
! returned in f. stpmax is an input quanity that limits the length of the
! steps so that the function is not evaluated in regions where it is
! undefined or subject to overflow. p is usually the newton direction. the
! output quantity check is false on normal exit, and true when x is too
! close to xold. in a minimization routine, this usually signals
! convergence and can be ignored. however, in a root finding routine, the
! calling routine should check wheather the convergence is spurious.

! declare the pass
      external         func
      logical          check
      integer          n,nfev
      double precision f,fold,stpmax,g(n),p(n),x(n),xold(n)


! locals
      integer          i
      double precision xfminx_cj,alf,tolx,a,alam,alam2,alamin,b, &
                       disc,f2,rhs1,rhs2,slope,sum,temp,test,tmplam
      parameter        (alf=1.0d-4, tolx=3.0d-13)

! alf ensures sufficient decrease in the function value
! tolx is the convergence criterion on deltax


! initialize and scale if the attempted step is too big
      check = .false.
      sum   = 0.0d0
      do i=1,n
       sum = sum + p(i)*p(i)
      enddo
      sum = sqrt(sum)
      if (sum .gt. stpmax) then
       do i=1,n
        p(i) = p(i) * stpmax/sum
       enddo
      end if
      slope = 0.0d0
      do i=1,n
       slope = slope + g(i)*p(i)
      enddo
      if (slope .ge. 0.0) stop 'roundoff problem in lnsrch_cj'


! compute lambda_min
      test = 0.0d0
      do i=1,n
       temp = abs(p(i))/max(abs(xold(i)),1.0d0)
       if (temp .gt. test) test = temp
      enddo
      alamin = tolx/test


! always try a full newton step, start of iteration loop
      alam = 1.0d0
1     continue
      do i=1,n
       x(i) = xold(i) + alam*p(i)
      enddo
      f    = xfminx_cj(x,func)
      nfev = nfev + 1

! convergence on deltax, for root finding, the calling routine
! should verify the convergence
      if (alam .lt. alamin) then
       do i=1,n
        x(i) = xold(i)
       enddo
       check = .true.
       return

! sufficient function decrease
      else if (f .le. fold + alf*alam*slope) then
       return

! backtrack
      else
       if (alam .eq. 1.0) then
        tmplam = -slope / (2.0d0 * (f-fold-slope))
       else
        rhs1 = f  - fold - alam*slope
        rhs2 = f2 - fold - alam2*slope
        a    = (rhs1/alam**2 - rhs2/alam2**2)/(alam-alam2)
        b    = (-alam2*rhs1/alam**2 + alam*rhs2/alam2**2) / (alam-alam2)
        if (a .eq. 0.0) then
         tmplam = -slope/(2.0d0 * b)
        else
         disc = b*b - 3.0d0 * a * slope
         if (disc .lt. 0.0) then
          tmplam = 0.5d0 * alam
         else if (b .le. 0.0) then
          tmplam = (-b + sqrt(disc)) / (3.0d0 * a)
         else
          tmplam = -slope/(b + sqrt(disc))
         end if
        end if
        if (tmplam .gt. 0.5d0*alam) tmplam = 0.5d0*alam
       end if
      end if

! store for the next trip through
      alam2 = alam
      f2    = f
      alam  = max(tmplam, 0.1d0*alam)
      go to 1
      end





      double precision function xfminx_cj(x,func)
      include 'implno.dek'

! returns f = 0.5 f dot f at x. func is a user supplied routine of the
! functions to be root found.
!
! common block communicates values back to routine xnewt

! declare
      external         func
      integer          nn,np,i
      parameter        (np = 4)
      double precision x(1),fvec(np),sum,dum
      common /newtcj/   fvec,nn

      call func(dum,x,fvec)
      sum = 0.0d0
      do i=1,nn
       sum = sum + fvec(i)*fvec(i)
      enddo
      xfminx_cj = 0.5d0 * sum
      return
      end






      subroutine jac_cj(x,y,dfdy,mcol,nrow,mmax,nmax,derivs)
      include 'implno.dek'

! this routine computes a second order accurate jacobian matrix
! of the function contained in the routine derivs.
!
! input is the point x and the the vector y at which to compute the
! jacobian dfdy.
!
! y has logical dimension nrow and physical dimension nmax,
! dfdy has logical dimension (mcol,nrow) and physical dimension (mmax,nmax)
!
! uses 2*nrow + 1 function evaluations


! declare the pass
      external         derivs
      integer          mcol,nrow,mmax,nmax
      double precision x,y(nmax),dfdy(mmax,nmax)


! locals
      integer          i,j,imax
      parameter        (imax = 4)
      double precision fminus(imax),fplus(imax),rel,ax,temp,h,hinv
      parameter        (rel = 3.162278d-8, ax = 1.0e-16)

! check
       if (nrow .gt. imax) stop 'nrow > imax in jac_cj'


! for each row, get the right stepsize
      do j=1,nrow
       temp = y(j)
       h    = rel * max(abs(y(j)),ax)
       y(j) = temp + h
       h    = y(j) - temp
       call derivs(x,y,fplus)
       y(j) = temp

       temp = y(j)
       y(j) = temp - h
       h    = temp - y(j)
       call derivs(x,y,fminus)
       y(j) = temp

! compute the jth row of the jacobian
        hinv = 1.0d0/(2.0d0 * h)
        do i=1,mcol
         dfdy(i,j) = (fplus(i) - fminus(i)) * hinv
        enddo
       enddo

! restore the original state
      call derivs(x,y,fplus)
      return
      end
!---------------------------------------------------------------------







!---------------------------------------------------------------------
! this file contains dense and special linear equation a*x=b solvers:
!
! lu decomposition:
! routine ludcmp does a pivoting lower-upper decomposition
! routine lubksb does the backsubstitution from ludcmp
! routine luinv inverts a matrix using ludcmp and lubksb
! routine ludet gets the determinant of a matrix using ludcmp



      subroutine ludcmp(a,n,np,indx,d)
      include 'implno.dek'
!
! given th matrix a(n,n), with physical dimsnsions a(np,ap) this routine
! replaces a by the lu decompostion of a row-wise permutation of itself.
! input are a,n,np. output is a, indx which records the row
! permutations effected by the partial pivoting, and d which is 1 if
! the number of interchanges is even, -1 if odd.
! use routine lubksb to solve a system of linear equations.
!
! nmax is the largest expected value of n
!
! declare
      integer          n,np,indx(np),nmax,i,j,k,imax
      parameter        (nmax=500)
      double precision a(np,np),d,tiny,vv(nmax),aamax,sum,dum
      parameter        (tiny=1.0d-20)

! bullet check
      if (np .gt. nmax) then
       write(6,*) 'np=',np,' nmax=',nmax
       stop 'np > nmax in routine ludcmp'
      end if

! vv stores the implicit scaling of each row
! loop over the rows to get the scaling information
      d = 1.0d0
      do i=1,n
       aamax = 0.0d0
       do j=1,n
        if (abs(a(i,j)) .gt. aamax) aamax = abs(a(i,j))
       enddo
       if (aamax .eq. 0.0) stop 'singular matrix in ludcmp'
       vv(i) = 1.0d0/aamax
      enddo

! for each column apply crouts method; see equation 2.3.12
      do j=1,n
       do i=1,j-1
        sum = a(i,j)
        do k=1,i-1
         sum = sum - a(i,k)*a(k,j)
        enddo
        a(i,j) = sum
       enddo

! find the largest pivot element
       aamax = 0.0d0
       do i=j,n
        sum=a(i,j)
        do k=1,j-1
         sum = sum - a(i,k)*a(k,j)
        enddo
        a(i,j) = sum
        dum = vv(i)*abs(sum)
        if (dum .ge. aamax) then
         imax  = i
         aamax = dum
        end if
       enddo

! if we need to interchange rows
       if (j .ne. imax) then
        do k=1,n
         dum       = a(imax,k)
         a(imax,k) = a(j,k)
         a(j,k)    = dum
        enddo
        d          = -d
        vv(imax)   = vv(j)
       end if

! divide by the pivot element
       indx(j) = imax
       if (a(j,j) .eq. 0.0) a(j,j) = tiny
       if (j .ne. n) then
        dum = 1.0d0/a(j,j)
        do i=j+1,n
         a(i,j) = a(i,j)*dum
        enddo
       end if

! and go back for another column of crouts method
      enddo
      return
      end




      subroutine lubksb(a,n,np,indx,b)
      include 'implno.dek'
!
! solves a set of n linear equations ax=b. a is input in its lu decomposition
! form, determined by the routine above ludcmp. indx is input as the
! permutation vector also returned by ludcmp. b is input as the right hand
! side vector and returns with the solution vector x.
! a,n ans np are not modified by this routine and thus can be left in place
! for successive calls (i.e matrix inversion)
!
!
! declare
      integer           n,np,indx(np),i,ii,j,ll
      double precision  a(np,np),b(np),sum

! when ii is > 0, ii becomes the index of the first nonzero element of b
! this is forward substitution of equation 2.3.6, and unscamble in place
      ii = 0
      do i=1,n
       ll = indx(i)
       sum = b(ll)
       b(ll) = b(i)
       if (ii .ne. 0) then
        do j=ii,i-1
         sum = sum - a(i,j) * b(j)
        enddo

! nonzero element was found, so dos the sums in the loop above
       else if (sum .ne. 0.0) then
        ii  = i
       end if
       b(i) = sum
      enddo

! back substitution equation 2.3.7
      do i = n,1,-1
       sum = b(i)
       if (i .lt. n) then
        do j=i+1,n
         sum = sum - a(i,j) * b(j)
        enddo
       end if
       b(i) = sum/a(i,i)
      enddo
      return
      end




      subroutine luinv(a,n,np,indx,y)
      include 'implno.dek'
!
! this routine takes as input the n by n matrix a, of physical dimension
! np by np and on output fills y with the inverse of a
!
! declare
      integer           n,np,i,j,indx(np)
      double precision  a(np,np),y(np,np),d

! set y to the identity matrix
      do j=1,n
       do i=1,n
        y(i,j) = 0.0d0
       enddo
      enddo
      do i=1,n
       y(i,i) = 1.0d0
      enddo

! decomp and backsubstitute each column
      call ludcmp(a,n,np,indx,d)
      do j=1,n
       call lubksb(a,n,np,indx,y(1,j))
      enddo
      return
      end






      double precision function ludet(a,n,np,indx)
      include 'implno.dek'
!
! this function takes as input the n by n matrix a, of physical dimension
! np by np and on output returns the determinate
! be carefull of raspy overflows
!
! declare
      integer           n,np,j,indx(np)
      double precision  a(np,np),d

! decomp
      call ludcmp(a,n,np,indx,d)
      do j=1,n
       d = d * a(j,j)
      enddo
      ludet = d
      return
      end
!---------------------------------------------------------------------







!---------------------------------------------------------------------

      subroutine net_pzextr_mesa(iest,xest,yest,yz,dy,nv)
      include 'implno.dek'

! use polynomial extrapolation to evaluate nv functions at x=0 by fitting
! a polynomial to a sequence of estimates with progressively smaller values
! x=xest, and corresponding function vectors yest(1:nv). the call is number
! iest in the sequence of calls. extrapolated function values are output as
! yz(1:nv), and their estimated error is output as dy(1:nv)

! declare the pass
     integer :: iest,nv
     real*8  :: xest,dy(nv),yest(nv),yz(nv)


! locals; qcol and x must be "saved" between successive calls
! nmax is the maximum number of isotopes

      integer            :: j,k1
      integer, parameter :: nmax=3200, imax=13
      real*8             :: delta,f1,f2,q,d(nmax),qcol(nmax,imax),x(imax)


! sanity checks

      if (iest .gt. imax) stop 'iest > imax in net_pzextr'
      if (nv .gt. nmax) stop 'nv > nmax in net_pzextr'



            x(iest) = xest
            do j=1,nv
               dy(j) = yest(j)
               yz(j) = yest(j)
            end do

            ! store first estimate in first column
            if (iest .eq. 1) then
               do j=1,nv
                  qcol(j,1) = yest(j)
               end do
            else
               do j=1,nv
                  d(j) = yest(j)
               end do
               do k1=1,iest-1
                  delta = 1.0d0/(x(iest-k1) - xest)
                  f1    = xest * delta
                  f2    = x(iest-k1) * delta
                  ! propagate tableu 1 diagonal more
                  do j=1,nv
                     q          = qcol(j,k1)
                     qcol(j,k1) = dy(j)
                     delta      = d(j) - q
                     dy(j)      = f1*delta
                     d(j)       = f2*delta
                     yz(j)      = yz(j) + dy(j)
                  end do
               end do
               do j=1,nv
                  qcol(j,iest) = dy(j)
               end do
            end if

      return
      end




      subroutine net_pzextr(iest,xest,yest,yz,dy,nv)
      include 'implno.dek'

! use polynomial extrapolation to evaluate nv functions at x=0 by fitting
! a polynomial to a sequence of estimates with progressively smaller values
! x=xest, and corresponding function vectors yest(1:nv). the call is number
! iest in the sequence of calls. extrapolated function values are output as
! yz(1:nv), and their estimated error is output as dy(1:nv)

! declare the pass
     integer :: iest,nv
     real*8  :: xest,dy(nv),yest(nv),yz(nv)


! locals; qcol and x must be "saved" between successive calls
! nmax is the maximum number of isotopes

      integer            :: j,k1
      integer, parameter :: nmax=3200, imax=13
      real*8             :: delta,f1,f2,q,d(nmax),qcol(nmax,imax),x(imax)


! sanity checks

      if (iest .gt. imax) stop 'iest > imax in net_pzextr'
      if (nv .gt. nmax) stop 'nv > nmax in net_pzextr'


! save current independent variables
      x(iest) = xest
      do j=1,nv
       dy(j) = yest(j)
       yz(j) = yest(j)
      enddo

! store first estimate in first column
      if (iest .eq. 1) then
       do j=1,nv
        qcol(j,1) = yest(j)
       enddo
      else
       do j=1,nv
        d(j) = yest(j)
       enddo
       do k1=1,iest-1
        delta = 1.0d0/(x(iest-k1) - xest)
        f1    = xest * delta
        f2    = x(iest-k1) * delta

! propagate tableu 1 diagonal more
        do j=1,nv
         q          = qcol(j,k1)
         qcol(j,k1) = dy(j)
         delta      = d(j) - q
         dy(j)      = f1*delta
         d(j)       = f2*delta
         yz(j)      = yz(j) + dy(j)
        enddo
       enddo
       do j=1,nv
        qcol(j,iest) = dy(j)
       enddo
      end if
      return
      end
!---------------------------------------------------------------------







!---------------------------------------------------------------------
      subroutine stifbs_ma28(y,dydx,nv,x,htry,eps,yscal,hdid,hnext, &
                             derivs,jakob,bjakob,nstp,ierr)
      include 'implno.dek'
      include 'sparse_matrix.dek'

! for sparse analytic jacobians, ma28 linear algebra
!
! semi-implicit extrapolation step for stiff ode's with monitoring
! of local truncation error to adjust stepsize. inputs are the dependent
! variable vector y(nv) and its derivative dydx(nv) at the starting of the
! independent variable x. also input are the stepsize to be attempted htry,
! the required accuracy eps, and the vector yscal against which the error is
! scaled. on output, y and x are replaced by their new values, hdid is the
! stepsize actually accomplished, and hnext is the estimated next stepsize.
! dervs is a user supplied function that computes the right hand side of
! the equations.
!
! 1/scalmx is the maximum increase in the step size allowed
!
! declare
      external         derivs,jakob,bjakob
      logical          first,reduct
      integer          nv,nmax,kmaxx,imax,ierr
      parameter        (nmax  = iodemax, &
                        kmaxx = 7, &
                        imax  = kmaxx+1)
      integer          i,iq,k,kk,km,kmax,kopt,nvold,nseq(imax)
      double precision y(nv),dydx(nv),x,htry,eps,yscal(nv),hdid,hnext, &
                       eps1,epsold,errmax,fact,h,red,scale,work,wrkmin, &
                       xest,xnew,a(imax),alf(kmaxx,kmaxx),err(kmaxx), &
                       yerr(nmax),ysav(nmax),yseq(nmax),safe1,safe2, &
                       redmax,redmin,tiny,scalmx,dum
      parameter        (safe1 = 0.25d0, safe2 = 0.7d0, redmax=1.0d-5, &
                        redmin = 0.7d0, tiny = 1.0d-30, &
                        scalmx = 0.1d0)
!     2                  scalmx = 0.5d0)
!     2                  scalmx = 0.3d0)
!     2                  scalmx = 0.2d0)



! for jacobian pictures
      character*20     string
      integer          nstp,ifirst,j
      double precision ans13(13,13),anydt(13),sum


! for the ma28 package
      integer          n5,n8,initmat
      parameter        (n5 = 5*iodemax, n8=8*iodemax)
      integer          ikeep(n5),iw(n8),flag
      double precision w(iodemax),u
      common /ma2c3/   w,u,iw,ikeep,flag



! assume that the independent variable is not explicit in the odes
      data             first/.true./, epsold/-1.0d0/, nvold/-1/
      data             nseq /2, 6, 10, 14, 22, 34, 50, 70/
      data             ifirst/0/, initmat/1/

! bullet check
      if (nv .gt. iodemax) stop 'nv bigger than iodemax in stifbs_ma28'


! initialize
!      if (initmat .eq. 1) then
!       initmat = 0
!       write(6,*) 'initializing pivots'
!      end if



! normal execution starts here
! a new tolerance or a new number of variables, so reinitialize
      if (eps .ne. epsold  .or.  nv .ne. nvold) then

!       write(6,*) 'initializing pivots in stifbs_ma28'

! initialize pivots
       nzo     = 0
       do i=1,naij
        iloc(i) = 0
        jloc(i) = 0
       end do

! get the sparse pattern
       call bjakob(iloc,jloc,nzo,naij)


! copy the location
       do i=1,nzo
        ivect(i) = iloc(i)
        jvect(i) = jloc(i)
       enddo


! force the diagonal to be the pivot elements
       do i=1,nzo
        amat(i) = 1.0d-10
        if (ivect(i) .eq. jvect(i)) amat(i) = 1.0d0
       enddo
       u  = 0.1d0
       call ma28ad(nv,nzo,amat,naij,iloc,naij,jloc,u,ikeep,iw,w,flag)
       if (flag .lt. 0 .and. flag .ne. -14 ) then
        write(6,*) 'error in ma28ad flag',flag
        stop 'error in ma28ad in stifbs_ma28'
       end if


       hnext = -1.0e29
       xnew  = -1.0e29
       eps1  = safe1 * eps

! compute the work coefficients a_k
       a(1)  = nseq(1) + 1
       do k=1,kmaxx
        a(k+1) = a(k) + nseq(k+1)
       enddo

! compute alf(k,q)
       do iq=2,kmaxx
        do k=1,iq-1
         alf(k,iq) = eps1**((a(k+1) - a(iq+1)) / &
                     ((a(iq+1) - a(1) + 1.0d0) * (2*k + 1)))
        enddo
       enddo
       epsold = eps
       nvold  = nv

! add cost of jacobians to work coefficients
       a(1)   = nv + a(1)
       do k=1,kmaxx
        a(k+1) = a(k) + nseq(k+1)
       enddo


! determine optimal row number for convergence
       do kopt=2,kmaxx-1
        if (a(kopt+1) .gt. a(kopt)*alf(kopt-1,kopt)) go to 01
       enddo
01     kmax = kopt
      end if




! save the starting values
      h = htry
      do i=1,nv
       ysav(i)  = y(i)
      enddo

! get the sparse jacobian in sparse_dfdy
      call jakob(x,y,sparse_dfdy,nzo)


! here are the jacobian values
!      if (x .gt. 5.0d10  .and. ifirst .eq. 0) then
!       ifirst = 1
!       do i=1,nzo
!        ans13(ivect(i),jvect(i)) = sparse_dfdy(i)
!       enddo
!       write(6,*)
!       write(6,102) ((ans13(i,j), i=1,nv), j=1,nv)
! 102   format(1x,1p13e11.3)
!       write(6,*)

!       read(5,*)
!      end if



! here is a picture of the jacobian
!      if (x .gt. 1.0d-8 .and. ifirst .eq. 0) then
!       ifirst = 1
!       write(6,*) 'making jacobian picture data at time=',x
!       write(string,111) 'jac_',nstp,'.dat'
! 111   format(a,i4.4,a)
!       call sqeeze(string)
!       open(unit=47,file=string,status='unknown')
!       do i=1,nzo
!        write(47,112) i,ivect(i),jvect(i),sparse_dfdy(i)
! 112    format(1x,3i5,1pe15.6)
!       enddo
!       close(unit=47)
!      end if





! a new stepsize or a new integration, re-establish the order window
      if (h .ne. hnext  .or.  x .ne. xnew) then
       first = .true.
       kopt = kmax
      end if
      reduct = .false.

! evaluate the sequence of semi implicit midpoint rules
! this loop is run through a minimum of two times
02    do 18 k=1,kmax

! see how many orders are needed
!       write(6,*) k,h

       xnew = x + h
       if (xnew .eq. x) stop 'step too small in routine stifbs_ma28'

       call simpr_ma28(ysav,dydx,nv,x,h,nseq(k),yseq,derivs,nstp)
       xest = (h/nseq(k))**2
       call net_pzextr(k,xest,yseq,y,yerr,nv)


! compute normalized error estimate and record which component is setting the scale
       if (k .ne. 1) then
        errmax = tiny
        ierr   = 0
        do i=1,nv
!        errmax = max(errmax,abs(yerr(i)/yscal(i)))
         dum = abs(yerr(i)/yscal(i))
         if (dum .ge. errmax) then
          errmax = dum
          ierr = i
         end if
        enddo

        errmax   = errmax/eps
        km = k - 1
        err(km) = (errmax/safe1)**(1.0d0/(2*km+1))
       end if

! in order window
       if (k .ne. 1  .and. (k .ge. kopt-1  .or. first)) then

! converged
        if (errmax .lt. 1.0) goto 04

!        write(6,114) ierr,yerr(ierr),yscal(ierr),errmax*eps,errmax
! 114    format(1x,i4,1p6e12.4)


! possible step size reductions
        if (k .eq. kmax  .or.  k .eq. kopt + 1) then
         red = safe2/err(km)
         go to 03
        else if (k .eq. kopt) then
         if (alf(kopt-1,kopt) .lt. err(km)) then
          red = 1.0d0/err(km)
          go to 03
         end if
        else if (kopt .eq. kmax) then
         if (alf(km,kmax-1) .lt. err(km)) then
          red = alf(km,kmax-1) * safe2/err(km)
          go to 03
         end if
        else if (alf(km,kopt) .lt. err(km)) then
         red = alf(km,kopt-1)/err(km)
         go to 03
        end if
       end if

18    continue

! reduce stepsize by at least redmin and at most redmax
03    red    = min(red,redmin)
      red    = max(red,redmax)
      h      = h * red
      ierr   = 0
      reduct = .true.
      go to 2

! successful step; get optimal row for convergence and corresponding stepsize
04    x = xnew
      hdid = h
      first = .false.
      wrkmin = 1.0e35
      do kk=1,km
       fact = max(err(kk),scalmx)
       work = fact * a(kk+1)
       if (work .lt. wrkmin) then
        scale  = fact
        wrkmin = work
        kopt   = kk + 1
       end if
      enddo

! check for possible order increase, but not if stepsize was just reduced
      hnext = h/scale
      if (kopt .ge. k  .and.  kopt .ne. kmax  .and.  .not.reduct) then
       fact = max(scale/alf(kopt-1,kopt),scalmx)
       if (a(kopt+1)*fact .le. wrkmin) then
        hnext = h/fact
        kopt = kopt + 1
       end if
      end if
      return
      end






      subroutine simpr_ma28(y,dydx,n,xs,htot,nstep,yout,derivs,nstp)
      include 'implno.dek'
      include 'sparse_matrix.dek'
!
! an implicit midpoint stepper, for ma28 sparse linear algebra.
!
! declare
      external         derivs
      integer          nmax,n,nstep,nmaxx
      parameter        (nmaxx=iodemax)
      integer          i,nn
      double precision y(n),dydx(n),xs,htot, &
                       yout(n),h,x,del(nmaxx),ytemp(nmaxx)

! for the ma28 package
      integer          n5,n8
      parameter        (n5 = 5*iodemax, n8=8*iodemax)
      integer          ikeep(n5),iw(n8),flag
      double precision w(iodemax),u
      common /ma2c3/   w,u,iw,ikeep,flag

! for jacobian pictures
      character*20     string
      integer          nstp



! stepsize this trip, and make the a matrix
      h = htot/nstep
      do i=1,nzo
       amat(i) = -h * sparse_dfdy(i)
       if (ivect(i) .eq. jvect(i)) amat(i) = 1.0d0 + amat(i)
      enddo

! here is a picture of the matrix being decomposed
!      write(string,111) 'mat_',nstp,'.dat'
! 111  format(a,i4.4,a)
!      call sqeeze(string)
!      open(unit=47,file=string,status='unknown')
!      do i=1,nzo
!       write(47,112) i,ivect(i),jvect(i),amat(i)
! 112   format(1x,3i5,1pe15.6)
!      enddo
!      close(unit=47)



! symbolic decomp, full partial pivot
!       do i=1,nzo
!        iloc(i) = ivect(i)
!        jloc(i) = jvect(i)
!       enddo
!      u  = 1.0d0
!      call ma28ad(n,nzo,amat,naij,iloc,naij,jloc,u,ikeep,iw,w,flag)
!      if (flag .lt. 0 .and. flag .ne. -14 ) then
!       write(6,*) 'error in ma28ad flag',flag
!       stop 'error in ma28ad in stifbs_ma28'
!      end if
!      do i=1,nzo
!       amat(i) = -h * sparse_dfdy(i)
!       if (ivect(i) .eq. jvect(i)) amat(i) = 1.0d0 + amat(i)
!      enddo



! numeric decomp
      call ma28bd(n,nzo,amat,naij,ivect,jvect,jloc,ikeep,iw,w,flag)

      if (flag .lt. 0 .and. flag .ne. -14 ) then
       write(6,*) 'flag',flag
       stop 'error in ma28bd in simpr_ma28'
      end if


! use yout as temporary storage; the first step
      do i=1,n
       yout(i) = h * dydx(i)
      enddo

      call ma28cd(n,amat,naij,jloc,ikeep,yout,w,1)

      do i=1,n
       del(i)   = yout(i)
       ytemp(i) = y(i) + del(i)
      enddo
      x = xs + h
      call derivs(x,ytemp,yout)


! use yout as temporary storage; general step
      do nn=2,nstep
       do i=1,n
        yout(i) = h*yout(i) - del(i)
       enddo

       call ma28cd(n,amat,naij,jloc,ikeep,yout,w,1)

       do i=1,n
        del(i)   = del(i) + 2.0d0 * yout(i)
        ytemp(i) = ytemp(i) + del(i)
       enddo
       x = x + h
       call derivs(x,ytemp,yout)
      enddo


! take the last step
      do i=1,n
       yout(i) = h * yout(i) - del(i)
      enddo

      call ma28cd(n,amat,naij,jloc,ikeep,yout,w,1)

      do i=1,n
       yout(i) = ytemp(i) + yout(i)
      enddo
      return
      end
!---------------------------------------------------------------------







!---------------------------------------------------------------------
! this file contains the harwell ma28 sparse matrix routines
!
! easy to use front end routines:
! routine ma28ad does the symbolic and numeric lu decomp
! routine ma28bd does the numeric lu decomp of ma28ad
! routine ma28cd solves the system of equations directly
! routine ma28id solves the system iteratively
! routine ma28dd does some pointer work
!
! these are the hardball routines:
! routine ma30ad does core symbolic and numeric lu decomp
! routine ma30bd does the numeric decomp the sparse pattern
! routine ma30cd solves the linear system
! routine ma30dd does garbage collection
!
! support hardball routines:
! routine ma28int1 does some common block initialization
! roytine ma28int2 does some common block initialization
! routine ma28int3 does some common block initialization
! routine mc20ad sort a matrix
! routine mc23ad does the block triangularization pointers
! routine mc22ad reorders the off diagonal blocks based on the pivot info
! routine mc21a front end of mc21b
! routine mc21b pernutes the rows to get a zero free diagonal
! routine mc13d front end for mc13e
! routine mc13e permutes a lower triangular block
! routine mc24ad gets a growth rate of fillin
!
! initialization routines (was block data routines)
! routine ma28int1 initializes the ma28 routine flags
! routine ma28int2 initializes the ma28 routine flags
! routine ma28int3 initializes the ma28 routine flags
!
! never called:
! routine mc20bd
!
!
!
!
!
      subroutine ma28ad(n,nz,a,licn,irn,lirn,icn,u,ikeep,iw,w,iflag)
      include 'implno.dek'
!
! this subroutine performs the lu factorization of a.
!
! input:
! n     order of matrix ... not altered by subroutine
! nz    number of non-zeros in input matrix ... not altered by subroutine
! a     is a real array  length licn.  holds non-zeros of matrix on entry
!       and non-zeros of factors on exit.  reordered by mc20a/ad and
!       mc23a/ad and altered by ma30a/ad
! licn  integer  length of arrays a and icn ... not altered by subroutine
! irn   integer array of length lirn.  holds row indices on input.
!       used as workspace by ma30a/ad to hold column orientation of matrix
! lirn  integer  length of array irn ... not altered by the subroutine
! icn   integer array of length licn.  holds column indices on entry
!       and column indices of decomposed matrix on exit. reordered by
!       mc20a/ad and mc23a/ad and altered by ma30a/ad.
! u     real variable  set by user to control bias towards numeric or
!       sparsity pivoting.  u=1.0 gives partial pivoting while u=0. does
!       not check multipliers at all.  values of u greater than one are
!       treated as one while negative values are treated as zero.  not
!       altered by subroutine.
! ikeep integer array of length 5*n  used as workspace by ma28a/ad
!       it is not required to be set on entry and, on exit, it contains
!       information about the decomposition. it should be preserved between
!       this call and subsequent calls to ma28b/bd or ma28c/cd.
!       ikeep(i,1),i=1,n  holds the total length of the part of row i
!       in the diagonal block.
!       row ikeep(i,2),i=1,n  of the input matrix is the ith row in
!       pivot order.
!       column ikeep(i,3),i=1,n  of the input matrix is the ith column
!       in pivot order.
!       ikeep(i,4),i=1,n  holds the length of the part of row i in
!       the l part of the l/u decomposition.
!       ikeep(i,5),i=1,n  holds the length of the part of row i in the
!       off-diagonal blocks.  if there is only one diagonal block,
!       ikeep(1,5) will be set to -1.
! iw    integer array of length 8*n.  if the option nsrch .le. n is
!       used, then the length of array iw can be reduced to 7*n.
! w     real array  length n.  used by mc24a/ad both as workspace and to
!       return growth estimate in w(1).  the use of this array by ma28a/ad
!       is thus optional depending on common block logical variable grow.
! iflag integer variable  used as error flag by routine.  a positive
!       or zero value on exit indicates success.  possible negative
!       values are -1 through -14.
!
! declare
      integer          n,nz,licn,lirn,iflag,irn(lirn),icn(licn), &
                       ikeep(n,5),iw(n,8),i,j1,j2,jj,j,length, &
                       newpos,move,newj1,jay,knum,ii,i1,iend
      double precision a(licn),u,w(n)
!
! common and private variables. common block ma28f/fd is used merely
! to communicate with common block ma30f/fd  so that the user
! need not declare this common block in his main program.
!
! the common block variables are:
! lp,mp    default value 6 (line printer).  unit number
!          for error messages and duplicate element warning resp.
! nlp,mlp  unit number for messages from ma30a/ad and
!          mc23a/ad resp.  set by ma28a/ad to value of lp.
! lblock   logical  default value true.  if true mc23a/ad is used
!          to first permute the matrix to block lower triangular form.
! grow     logical  default value true.  if true then an estimate
!          of the increase in size of matrix elements during l/u
!          decomposition is given by mc24a/ad.
! eps,rmin,resid  variables not referenced by ma28a/ad.
! irncp,icncp  set to number of compresses on arrays irn and icn/a
! minirn,minicn  minimum length of arrays irn and icn/a; for success on
!                future runs.
! irank  integer   estimated rank of matrix.
! mirncp,micncp,mirank,mirn,micn communicate between ma30f/fd and ma28f/fd
!                                values of above named variables with
!                                somewhat similar names.
! abort1,abort2  logical variables with default value true.  if false
!                then decomposition will be performed even if the matrix is
!                structurally or numerically singular respectively.
! aborta,abortb  logical variables used to communicate values of
!                 abort1 and abort2 to ma30a/ad.
! abort  logical  used to communicate value of abort1 to mc23a/ad.
! abort3  logical variable not referenced by ma28a/ad.
! idisp  integer array  length 2.  used to communicate information
!        on decomposition between this call to ma28a/ad and subsequent
!        calls to ma28b/bd and ma28c/cd.  on exit, idisp(1) and
!        idisp(2) indicate position in arrays a and icn of the
!        first and last elements in the l/u decomposition of the
!        diagonal blocks, respectively.
! numnz  integer  structural rank of matrix.
! num    integer  number of diagonal blocks.
! large  integer  size of largest diagonal block.
!
!
      logical          grow,lblock,abort,abort1,abort2,abort3,aborta, &
                       abortb,lbig,lbig1
      integer          idisp(2),lp,mp,irncp,icncp,minirn,minicn, &
                       irank,ndrop,maxit,noiter,nsrch,istart, &
                       ndrop1,nsrch1,nlp,mirncp,micncp,mirank, &
                       mirn,micn,mlp,numnz,num,large,lpiv(10), &
                       lnpiv(10),mapiv,manpiv,iavpiv,ianpiv,kountl, &
                       ifirst
      double precision tol,themax,big,dxmax,errmax,dres,cgce, &
                       tol1,big1,upriv,rmin,eps,resid,zero
!
      common /ma28ed/  lp,mp,lblock,grow
      common /ma28fd/  eps,rmin,resid,irncp,icncp,minirn,minicn, &
                       irank,abort1,abort2
      common /ma28gd/  idisp
      common /ma28hd/  tol,themax,big,dxmax,errmax,dres,cgce, &
                       ndrop,maxit,noiter,nsrch,istart,lbig
      common /ma30id/  tol1,big1,ndrop1,nsrch1,lbig1
      common /ma30ed/  nlp,aborta,abortb,abort3
      common /ma30fd/  mirncp,micncp,mirank,mirn,micn
      common /mc23bd/  mlp,numnz,num,large,abort
      common /lpivot/  lpiv,lnpiv,mapiv,manpiv,iavpiv,ianpiv,kountl
      data zero        /0.0d0/, ifirst/0/
!
! format statements
99990 format(1x,'error return from ma28a/ad because')
99991 format(1x,'error return from ma30a/ad')
99992 format(1x,'error return from mc23a/ad')
99993 format(1x,'duplicate element in position',i8,' ',i8, &
                'with value ',1pe22.14)
99994 format (1x,i6,'element with value',1pe22.14,'is out of range',/, &
              1x,'with indices',i8,' ',i8)
99995 format(1x,'error return from ma28a/ad; indices out of range')
99996 format(1x,'lirn too small = ',i10)
99997 format(1x,'licn too small = ',i10)
99998 format(1x,'nz non positive = ',i10)
99999 format(1x,'n out of range = ',i10)
!
!
! initialization and transfer of information between common blocks
      if (ifirst .eq. 0) then
       ifirst = 1
       call ma28int1
       call ma28int2
       call ma28int3
      end if
      iflag  = 0
      aborta = abort1
      abortb = abort2
      abort  = abort1
      mlp    = lp
      nlp    = lp
      tol1   = tol
      lbig1  = lbig
      nsrch1 = nsrch
!
! upriv private copy of u is used in case it is outside
      upriv = u
!
! simple data check on input variables and array dimensions.
      if (n .gt. 0) go to 10
      iflag = -8
      if (lp .ne. 0) write (lp,99999) n
      go to 210
10    if (nz .gt. 0) go to 20
      iflag = -9
      if (lp .ne. 0) write (lp,99998) nz
      go to 210
20    if (licn .ge. nz) go to 30
      iflag = -10
      if (lp .ne. 0) write (lp,99997) licn
      go to 210
30    if (lirn .ge. nz) go to 40
      iflag = -11
      if (lp .ne. 0) write (lp,99996) lirn
      go to 210
!
! data check to see if all indices lie between 1 and n.
40    do 50 i=1,nz
       if (irn(i) .gt. 0 .and. irn(i) .le. n .and. icn(i) .gt. 0 .and. &
           icn(i) .le. n) go to 50
       if (iflag .eq. 0 .and. lp .ne. 0) write (lp,99995)
       iflag = -12
       if (lp .ne. 0) write (lp,99994) i,a(i),irn(i),icn(i)
50    continue
      if (iflag .lt. 0) go to 220
!
! sort matrix into row order.
      call mc20ad(n,nz,a,icn,iw,irn,0)
!
! part of ikeep is used here as a work-array.  ikeep(i,2) is the last row to
! have a non-zero in column i.  ikeep(i,3) is the off-set of column i from
! the start of the row.
      do 60 i=1,n
       ikeep(i,2) = 0
       ikeep(i,1) = 0
60    continue
!
! check for duplicate elements .. summing any such entries and printing a
! warning message on unit mp. move is equal to the number of duplicate
! elements found; largest element in the matrix is themax; j1 is position in
! arrays of first non-zero in row.
      move   = 0
      themax = zero
      j1     = iw(1,1)
      do 130 i=1,n
       iend = nz + 1
       if (i .ne. n) iend = iw(i+1,1)
       length = iend - j1
       if (length .eq. 0) go to 130
       j2 = iend - 1
       newj1 = j1 - move
       do 120 jj=j1,j2
        j = icn(jj)
        themax = max(themax,abs(a(jj)))
        if (ikeep(j,2) .eq. i) go to 110
!
! first time column has ocurred in current row.
        ikeep(j,2) = i
        ikeep(j,3) = jj - move - newj1
        if (move .eq. 0) go to 120
!
! shift necessary because of previous duplicate element.
        newpos = jj - move
        a(newpos) = a(jj)
        icn(newpos) = icn(jj)
        go to 120
!
! duplicate element.
110     move = move + 1
        length = length - 1
        jay = ikeep(j,3) + newj1
        if (mp .ne. 0) write (mp,99993) i,j,a(jj)
        a(jay) = a(jay) + a(jj)
        themax = max(themax,abs(a(jay)))
120    continue
       ikeep(i,1) = length
       j1 = iend
130    continue
!
! knum is actual number of non-zeros in matrix with any multiple entries
! counted only once
      knum = nz - move
      if (.not.lblock) go to 140
!
! perform block triangularisation.
      call mc23ad(n,icn,a,licn,ikeep,idisp,ikeep(1,2), &
                  ikeep(1,3),ikeep(1,5),iw(1,3),iw)
      if (idisp(1) .gt. 0) go to 170
      iflag = -7
      if (idisp(1) .eq. -1) iflag = -1
      if (lp .ne. 0) write (lp,99992)
      go to 210
!
! block triangularization not requested. move structure to end of data arrays
! in preparation for ma30a/ad; set lenoff(1) to -1 and set permutation arrays.
140   do 150 i=1,knum
       ii = knum - i + 1
       newpos = licn - i + 1
       icn(newpos) = icn(ii)
       a(newpos) = a(ii)
150   continue
      idisp(1) = 1
      idisp(2) = licn - knum + 1
      do 160 i=1,n
       ikeep(i,2) = i
       ikeep(i,3) = i
160   continue
      ikeep(1,5) = -1
170   if (lbig) big1 = themax
      if (nsrch .le. n) go to 180
!
! perform l/u decomosition on diagonal blocks.
      call ma30ad(n,icn,a,licn,ikeep,ikeep(1,4),idisp, &
                 ikeep(1,2),ikeep(1,3),irn,lirn,iw(1,2),iw(1,3),iw(1,4), &
                 iw(1,5),iw(1,6),iw(1,7),iw(1,8),iw,upriv,iflag)
      go to 190
!
! this call if used if nsrch has been set less than or equal n; in this case,
! two integer work arrays of length can be saved.
180    call ma30ad(n,icn,a,licn,ikeep,ikeep(1,4),idisp, &
                 ikeep(1,2),ikeep(1,3),irn,lirn,iw(1,2),iw(1,3),iw(1,4), &
                 iw(1,5),iw,iw,iw(1,6),iw,upriv,iflag)
!
! transfer common block information.
190   minirn = max0(mirn,nz)
      minicn = max0(micn,nz)
      irncp = mirncp
      icncp = micncp
      irank = mirank
      ndrop = ndrop1
      if (lbig) big = big1
      if (iflag .ge. 0) go to 200
      if (lp .ne. 0) write (lp,99991)
      go to 210
!
! reorder off-diagonal blocks according to pivot permutation.
200   i1 = idisp(1) - 1
      if (i1 .ne. 0) call mc22ad(n,icn,a,i1,ikeep(1,5),ikeep(1,2), &
                               ikeep(1,3),iw,irn)
      i1 = idisp(1)
      iend = licn - i1 + 1
!
! optionally calculate element growth estimate.
      if (grow) call mc24ad(n,icn,a(i1),iend,ikeep,ikeep(1,4),w)
!
! increment growth estimate by original maximum element.
      if (grow) w(1) = w(1) + themax
      if (grow .and. n .gt. 1) w(2) = themax
!
! set flag if the only error is due to duplicate elements.
      if (iflag .ge. 0 .and. move .ne. 0) iflag = -14
      go to 220
210   if (lp .ne. 0) write (lp,99990)
220   return
      end
!
!
!
!
!
      subroutine ma28bd(n,nz,a,licn,ivect,jvect,icn,ikeep,iw,w,iflag)
      include 'implno.dek'
!
! this subroutine factorizes a matrix with the same pattern as that
! previously factorized by ma28a/ad.
!
! input is :
! n      order of matrix  not altered by subroutine.
! nz     number of non-zeros in input matrix  not altered by subroutine.
! a      array  length licn.  holds non-zeros of matrix on entry and
!        non-zeros of factors on exit.  reordered by ma28d/dd and altered by
!        subroutine ma30b/bd.
! licn   integer  length of arrays a and icn.  not altered by subroutine.
! ivect,jvect  integer arrays of length nz.  hold row and column
!        indices of non-zeros respectively.  not altered by subroutine.
! icn    integer array of length licn.  same array as output from ma28a/ad.
!        unchanged by ma28b/bd.
! ikeep  integer array of length 5*n.  same array as output from
!        ma28a/ad.  unchanged by ma28b/bd.
! iw     integer array  length 5*n.  used as workspace by ma28d/dd and
!        ma30b/bd.
! w      array  length n.  used as workspace by ma28d/dd,ma30b/bd and
!        (optionally) mc24a/ad.
! iflag  error flag with positive or zero value indicating success.
!
!
! declare
      integer          n,nz,licn,iw(n,5),iflag,ikeep(n,5),ivect(nz), &
                       jvect(nz),icn(licn),i1,iend,idup
      double precision a(licn),w(n)
!
! private and common variables: unless otherwise stated common block
! variables are as in ma28a/ad. those variables referenced by ma28b/bd are
! mentioned below.
!
! lp,mp  used as in ma28a/ad as unit number for error and
!        warning messages, respectively.
! nlp    variable used to give value of lp to ma30e/ed.
! eps    real/double precision  ma30b/bd will output a positive value
!        for iflag if any modulus of the ratio of pivot element to the
!        largest element in its row (u part only) is less than eps (unless
!        eps is greater than 1.0 when no action takes place).
! rmin   variable equal to the value of this minimum ratio in cases where
!        eps is less than or equal to 1.0.
! meps,mrmin variables used by the subroutine to communicate between common
!         blocks ma28f/fd and ma30g/gd.
!
! declare
      logical          grow,lblock,aborta,abortb,abort1,abort2,abort3, &
                       lbig,lbig1
      integer          idisp(2),mp,lp,irncp,icncp,minirn,minicn,irank, &
                       ndrop,maxit,noiter,nsrch,istart,nlp,ndrop1,nsrch1
      double precision eps,meps,rmin,mrmin,resid,tol,themax,big,dxmax, &
                       errmax,dres,cgce,tol1,big1
      common /ma28ed/  mp,lp,lblock,grow
      common /ma28fd/  eps,rmin,resid,irncp,icncp,minirn,minicn,irank, &
                       abort1,abort2
      common /ma28gd/  idisp
      common /ma28hd/  tol,themax,big,dxmax,errmax,dres,cgce,ndrop, &
                       maxit,noiter,nsrch,istart,lbig
      common /ma30ed/  nlp,aborta,abortb,abort3
      common /ma30gd/  meps,mrmin
      common /ma30id/  tol1,big1,ndrop1,nsrch1,lbig1
!
! formats
99994 format(1x,'error return from ma28b/bd because')
99995 format(1x,'error return from ma30b/bd')
99996 format(1x,'licn too small = ',i10)
99997 format(1x,'nz non positive = ',i10)
99998 format(1x,'n out of range = ',i10)
99999 format(1x,'error return from ma28b/bd with iflag=',i4,/, &
             1x,i7,' entries dropped from structure by ma28a/ad')
!
!
! check to see if elements were dropped in previous ma28a/ad call.
      if (ndrop .eq. 0) go to 10
      iflag = -15
      write (6,99999) iflag,ndrop
      go to 70
10    iflag = 0
      meps  = eps
      nlp   = lp
!
! simple data check on variables.
      if (n .gt. 0) go to 20
      iflag = -11
      if (lp .ne. 0) write (lp,99998) n
      go to 60
20    if (nz .gt. 0) go to 30
      iflag = -10
      if (lp .ne. 0) write (lp,99997) nz
      go to 60
30    if (licn .ge. nz) go to 40
      iflag = -9
      if (lp .ne. 0) write (lp,99996) licn
      go to 60
!
!
40     call ma28dd(n,a,licn,ivect,jvect,nz,icn,ikeep,ikeep(1,4), &
                   ikeep(1,5),ikeep(1,2),ikeep(1,3),iw(1,3),iw, &
                   w(1),iflag)
!
! themax is largest element in matrix
      themax = w(1)
      if (lbig) big1 = themax
!
! idup equals one if there were duplicate elements, zero otherwise.
      idup = 0
      if (iflag .eq. (n+1)) idup = 1
      if (iflag .lt. 0) go to 60
!
! perform row-gauss elimination on the structure received from ma28d/dd
      call ma30bd(n,icn,a,licn,ikeep,ikeep(1,4),idisp, &
                  ikeep(1,2),ikeep(1,3),w,iw,iflag)
!
! transfer common block information.
      if (lbig) big1 = big
      rmin = mrmin
      if (iflag .ge. 0) go to 50
      iflag = -2
      if (lp .ne. 0) write (lp,99995)
      go to 60
!
! optionally calculate the growth parameter.
50    i1   = idisp(1)
      iend = licn - i1 + 1
      if (grow) call mc24ad(n,icn,a(i1),iend,ikeep,ikeep(1,4),w)
!
! increment estimate by largest element in input matrix.
      if (grow) w(1) = w(1) + themax
      if (grow .and. n .gt. 1) w(2) = themax
!
! set flag if the only error is due to duplicate elements.
      if (idup .eq. 1 .and. iflag .ge. 0) iflag = -14
      go to 70
60    if (lp .ne. 0) write (lp,99994)
70    return
      end
!
!
!
!
!
      subroutine ma28cd(n,a,licn,icn,ikeep,rhs,w,mtype)
      include 'implno.dek'
!
! uses the factors from ma28a/ad or ma28b/bd to solve a system of equations
!
! input:
! n     order of matrix  not altered by subroutine.
! a     array  length licn.  the same array as most recent call to ma28a/ad
!       or ma28b/bd.
! licn  length of arrays a and icn.  not altered by subroutine.
! icn   integer array of length licn.  same array as output from ma28a/ad.
!       unchanged by ma28c/cd.
! ikeep integer array of length 5*n.  same array as output from ma28a/ad.
!       unchanged by ma28c/cd.
! rhs   array  length n.  on entry, it holds the right hand side.
!       on exit, the solution vector.
! w     array  length n. used as workspace by ma30c/cd.
! mtype integer  used to tell ma30c/cd to solve the direct equation
!       (mtype=1) or its transpose (mtype .ne. 1).
!
! resid  variable returns maximum residual of equations where pivot was zero.
! mresid variable used by ma28c/cd to communicate with ma28f/fd and ma30h/hd.
! idisp  integer array ; the same as that used by ma28a/ad. un changed.
!
! declare
      logical          abort1,abort2
      integer          n,licn,idisp(2),icn(licn),ikeep(n,5), &
                       irncp,icncp,minirn,minicn,irank,mtype
      double precision a(licn),rhs(n),w(n),resid,mresid,eps,rmin
      common /ma28fd/  eps,rmin,resid,irncp,icncp,minirn,minicn, &
                       irank,abort1,abort2
      common /ma28gd/  idisp
      common /ma30hd/  mresid
!
! this call performs the solution of the set of equations.
      call ma30cd(n,icn,a,licn,ikeep,ikeep(1,4),ikeep(1,5), &
                  idisp,ikeep(1,2),ikeep(1,3),rhs,w,mtype)
!
! transfer common block information.
      resid = mresid
      return
      end
!
!
!
!
!
      subroutine ma28id(n,nz,aorg,irnorg,icnorg,licn,a,icn, &
                        ikeep,rhs,x,r,w,mtype,prec,iflag)
      include 'implno.dek'
!
! this subroutine uses the factors from an earlier call to ma28a/ad
! or ma28b/bd to solve the system of equations with iterative refinement.
!
! parameters are:
!
! n    order of the matrix. it is not altered by the subroutine.
! nz   number of entries in the original matrix.  not altered by subroutine.
!      for this entry the original matrix must have been saved in
!      aorg,irnorg,icnorg where entry aorg(k) is in row irnorg(k) and
!      column icnorg(k), k=1,...nz.  information about the factors of a
!      is communicated to this subroutine via the parameters licn, a, icn
!      and ikeep where:
! aorg   array of length nz.  not altered by ma28i/id.
! irnorg array of length nz.  not altered by ma28i/id.
! icnorg array of length nz.  not altered by ma28i/id.
! licn   equal to the length of arrays a and icn. not altered
! a    array of length licn. it must be unchanged since the last call
!      to ma28a/ad or ma28b/bd. it is not altered by the subroutine.
! icn, ikeep are the arrays (of lengths licn and 5*n, respectively) of
!      the same names as in the previous all to ma28a/ad. they should be
!      unchanged since this earlier call. not altered.
!
! other parameters are as follows:
! rhs array of length n. the user must set rhs(i) to contain the
!     value of the i th component of the right hand side. not altered.
!
! x   array of length n. if an initial guess of the solution is
!     given (istart equal to 1), then the user must set x(i) to contain
!     the value of the i th component of the estimated solution.  on
!     exit, x(i) contains the i th component of the solution vector.
! r   array of length n. it need not be set on entry.  on exit, r(i)
!     contains the i th component of an estimate of the error if maxit
!     is greater than 0.
! w is an array of length n. it is used as workspace by ma28i/id.
! mtype must be set to determine whether ma28i/id will solve a*x=rhs
!      (mtype equal to 1) or at*x=rhs (mtype ne 1, zero say). not altered.
! prec should be set by the user to the relative accuracy required. the
!      iterative refinement will terminate if the magnitude of the
!      largest component of the estimated error relative to the largest
!      component in the solution is less than prec. not altered.
! iflag is a diagnostic flag which will be set to zero on successful
!       exit from ma28i/id, otherwise it will have a non-zero value. the
!       non-zero value iflag can have on exit from ma28i/id are ...
!       -16    indicating that more than maxit iteartions are required.
!       -17    indicating that more convergence was too slow.
!
! declare
      integer          n,nz,licn,mtype,iflag,icnorg(nz),irnorg(nz), &
                       ikeep(n,5),icn(licn),i,iterat,nrow,ncol
      double precision a(licn),aorg(nz),rhs(n),r(n),x(n),w(n),prec, &
                       d,dd,conver,zero
!
! common block communication
      logical          lblock,grow,lbig
      integer          lp,mp,ndrop,maxit,noiter,nsrch,istart
      double precision tol,themax,big,dxmax,errmax,dres,cgce
      common /ma28ed/  lp,mp,lblock,grow
      common /ma28hd/  tol,themax,big,dxmax,errmax,dres,cgce, &
                       ndrop,maxit,noiter,nsrch,istart,lbig
      data             zero /0.0d0/
!
!
! formats
99998 format(1x,'error return from ma28i with iflag = ', i3,/, &
             1x,'convergence rate of',1pe9.2,'too slow',/, &
             1x,'maximum acceptable rate set to ',1pe9.2)
99999 format(1x,'error return from ma28i/id with iflag = ',i3,/, &
             1x,'more than',i5,'iterations required')
!
!
! initialization of noiter,errmax and iflag.
      noiter = 0
      errmax = zero
      iflag   = 0
!
! jump if a starting vector has been supplied by the user.
      if (istart .eq. 1) go to 20
!
! make a copy of the right-hand side vector.
      do 10 i=1,n
       x(i) = rhs(i)
10    continue
!
! find the first solution.
      call ma28cd(n,a,licn,icn,ikeep,x,w,mtype)
!
! stop the computations if   maxit=0.
20    if (maxit .eq. 0) go to 160
!
! calculate the max-norm of the first solution.
      dd = 0.0d0
      do 30 i=1,n
       dd = max(dd,abs(x(i)))
30    continue
      dxmax = dd
!
! begin the iterative process.
      do 120 iterat=1,maxit
       d = dd
!
! calculate the residual vector.
       do 40 i=1,n
        r(i) = rhs(i)
40     continue
       if (mtype .eq. 1) go to 60
       do 50 i=1,nz
        nrow = irnorg(i)
        ncol = icnorg(i)
        r(ncol) = r(ncol) - aorg(i)*x(nrow)
50     continue
       go to 80
!
! mtype=1.
60     do 70 i=1,nz
        nrow = irnorg(i)
        ncol = icnorg(i)
        r(nrow) = r(nrow) - aorg(i)*x(ncol)
70     continue
80     dres = 0.0d0
!
! find the max-norm of the residual vector.
       do 90 i=1,n
        dres = max(dres,abs(r(i)))
90     continue
!
! stop the calculations if the max-norm of the residual vector is zero.
       if (dres .eq. 0.0) go to 150
!
! calculate the correction vector.
       noiter = noiter + 1
       call ma28cd(n,a,licn,icn,ikeep,r,w,mtype)
!
! find the max-norm of the correction vector.
       dd = 0.0d0
       do 100 i=1,n
        dd = max(dd,abs(r(i)))
100    continue
!
! check the convergence.
       if (dd .gt. d*cgce .and. iterat .ge. 2) go to 130
       if (dxmax*10.0d0 + dd .eq. dxmax*10.0d0) go to 140
!
! attempt to improve the solution.
       dxmax = 0.0d0
       do 110 i=1,n
        x(i) = x(i) + r(i)
        dxmax = max(dxmax,abs(x(i)))
110    continue
!
! check the stopping criterion; end of iteration loop
       if (dd .lt. prec*dxmax) go to 140
120   continue
!
! more than maxit iterations required.
      iflag = -16
      write (lp,99999) iflag,maxit
      go to 140
!
! convergence rate unacceptably slow.
130   iflag = -17
      conver = dd/d
      write (lp,99998) iflag,conver,cgce
!
! the iterative process is terminated.
140   errmax = dd
150   continue
160   return
      end
!
!
!
!
!
      subroutine ma28dd(n,a,licn,ivect,jvect,nz,icn,lenr,lenrl, &
                        lenoff,ip,iq,iw1,iw,w1,iflag)
      include 'implno.dek'
!
! this subroutine need never be called by the user directly.
! it sorts the user's matrix into the structure of the decomposed
! form and checks for the presence of duplicate entries or
! non-zeros lying outside the sparsity pattern of the decomposition
! it also calculates the largest element in the input matrix.
!
! declare
      logical          lblock,grow,blockl
      integer          n,licn,nz,iw(n,2),idisp(2),icn(licn),ivect(nz), &
                       jvect(nz),ip(n),iq(n),lenr(n),iw1(n,3),lenrl(n), &
                       lenoff(n),iflag,lp,mp,i,ii,jj,inew,jnew,iblock, &
                       iold,jold,j1,j2,idisp2,idummy,jdummy,midpt,jcomp
      double precision a(licn),zero,w1,aa
!
      common /ma28ed/  lp,mp,lblock,grow
      common /ma28gd/  idisp
      data             zero /0.0d0/
!
! formats
99997 format(1x,'element',i6,' ',i6,' was not in l/u pattern')
99998 format(1x,'non-zero',i7,' ',i6,' in zero off-diagonal block')
99999 format(1x,'element ',i6,' with value ',1pe22.14,/, &
             1x,'has indices',i8,' ',i8,'out of range')
!
! iw1(i,3)  is set to the block in which row i lies and the inverse
! permutations to ip and iq are set in iw1(.,1) and iw1(.,2) resp.
! pointers to beginning of the part of row i in diagonal and off-diagonal
! blocks are set in iw(i,2) and iw(i,1) resp.
      blockl  = lenoff(1) .ge. 0
      iblock  = 1
      iw(1,1) = 1
      iw(1,2) = idisp(1)
      do 10 i=1,n
       iw1(i,3) = iblock
       if (ip(i) .lt. 0) iblock = iblock + 1
       ii = iabs(ip(i)+0)
       iw1(ii,1) = i
       jj = iq(i)
       jj = iabs(jj)
       iw1(jj,2) = i
       if (i .eq. 1) go to 10
       if (blockl) iw(i,1) = iw(i-1,1) + lenoff(i-1)
       iw(i,2) = iw(i-1,2) + lenr(i-1)
10    continue
!
! place each non-zero in turn into its correct location in the a/icn array.
      idisp2 = idisp(2)
      do 170 i=1,nz
       if (i .gt. idisp2) go to 20
       if (icn(i) .lt. 0) go to 170
20     iold = ivect(i)
       jold = jvect(i)
       aa = a(i)
!
! dummy loop for following a chain of interchanges. executed nz times.
       do 140 idummy=1,nz
        if (iold .le. n .and. iold .gt. 0 .and. jold .le. n .and. &
            jold .gt. 0) go to 30
        if (lp .ne. 0) write (lp,99999) i, a(i), iold, jold
        iflag = -12
        go to 180
30      inew = iw1(iold,1)
        jnew = iw1(jold,2)
!
! are we in a valid block and is it diagonal or off-diagonal?
        if (iw1(inew,3)-iw1(jnew,3)) 40, 60, 50
40      iflag = -13
        if (lp .ne. 0) write (lp,99998) iold, jold
        go to 180
50      j1 = iw(inew,1)
        j2 = j1 + lenoff(inew) - 1
        go to 110
!
! element is in diagonal block.
60      j1 = iw(inew,2)
        if (inew .gt. jnew) go to 70
        j2 = j1 + lenr(inew) - 1
        j1 = j1 + lenrl(inew)
        go to 110
70      j2 = j1 + lenrl(inew)
!
! binary search of ordered list  .. element in l part of row.
        do 100 jdummy=1,n
         midpt = (j1+j2)/2
         jcomp = iabs(icn(midpt)+0)
         if (jnew-jcomp) 80, 130, 90
80       j2 = midpt
         go to 100
90       j1 = midpt
100     continue
        iflag = -13
        if (lp .ne. 0) write (lp,99997) iold, jold
        go to 180
!
! linear search ... element in l part of row or off-diagonal blocks.
110     do 120 midpt=j1,j2
         if (iabs(icn(midpt)+0) .eq. jnew) go to 130
120     continue
        iflag = -13
        if (lp .ne. 0) write (lp,99997) iold, jold
        go to 180
!
! equivalent element of icn is in position midpt.
130     if (icn(midpt) .lt. 0) go to 160
        if (midpt .gt. nz .or. midpt .le. i) go to 150
        w1 = a(midpt)
        a(midpt) = aa
        aa = w1
        iold = ivect(midpt)
        jold = jvect(midpt)
        icn(midpt) = -icn(midpt)
140    continue
!
150    a(midpt) = aa
       icn(midpt) = -icn(midpt)
       go to 170
160    a(midpt) = a(midpt) + aa
!
! set flag for duplicate elements; end of big loop
       iflag = n + 1
170   continue
!
! reset icn array  and zero elements in l/u but not in a. get max a element
180   w1 = zero
      do 200 i=1,idisp2
       if (icn(i) .lt. 0) go to 190
       a(i) = zero
       go to 200
190    icn(i) = -icn(i)
       w1 = max(w1,abs(a(i)))
200   continue
      return
      end
!
!
!
!
!
      subroutine ma30ad(nn,icn,a,licn,lenr,lenrl,idisp,ip,iq, &
                        irn,lirn,lenc,ifirst,lastr,nextr,lastc, &
                        nextc,iptr,ipc,u,iflag)
      include 'implno.dek'
!
!
! if the user requires a more convenient data interface then the ma28
! package should be used.  the ma28 subroutines call the ma30 routines after
! checking the user's input data and optionally using mc23a/ad to permute the
! matrix to block triangular form.
!
! this package of subroutines (ma30a/ad, ma30b/bd, ma30c/cd and ma30d/dd)
! performs operations pertinent to the solution of a general sparse n by n
! system of linear equations (i.e. solve ax=b). structually singular matrices
! are permitted including those with row or columns consisting entirely of
! zeros (i.e. including rectangular matrices).  it is assumed that the
! non-zeros of the matrix a do not differ widely in size. if necessary a
! prior call of the scaling subroutine mc19a/ad may be made.
!
! a discussion of the design of these subroutines is given by duff and reid
! (acm trans math software 5 pp 18-35,1979 (css 48)) while fuller details of
! the implementation are given in duff (harwell report aere-r 8730,1977).
! the additional pivoting option in ma30a/ad and the use of drop tolerances
! (see common block ma30i/id) were added to the package after joint work with
! duff,reid,schaumburg,wasniewski and zlatev, harwell report css 135, 1983.
!
! ma30a/ad performs the lu decomposition of the diagonal blocks of the
! permutation paq of a sparse matrix a, where input permutations p1 and q1
! are used to define the diagonal blocks.  there may be non-zeros in the
! off-diagonal blocks but they are unaffected by ma30a/ad. p and p1 differ
! only within blocks as do q and q1. the permutations p1 and q1 may be found
! by calling mc23a/ad or the matrix may be treated as a single block by
! using p1=q1=i. the matrix non-zeros should be held compactly by rows,
! although it should be noted that the user can supply the matrix by columns
! to get the lu decomposition of a transpose.
!
! this description of the following parameters should also be consulted for
! further information on most of the parameters of ma30b/bd and ma30c/cd:
!
! n    is an integer variable which must be set by the user to the order
!      of the matrix.  it is not altered by ma30a/ad.
!
! icn  is an integer array of length licn. positions idisp(2) to
!      licn must be set by the user to contain the column indices of
!      the non-zeros in the diagonal blocks of p1*a*q1. those belonging
!      to a single row must be contiguous but the ordering of column
!      indices with each row is unimportant. the non-zeros of row i
!      precede those of row i+1,i=1,...,n-1 and no wasted space is
!      allowed between the rows.  on output the column indices of the
!      lu decomposition of paq are held in positions idisp(1) to
!      idisp(2), the rows are in pivotal order, and the column indices
!      of the l part of each row are in pivotal order and precede those
!      of u. again there is no wasted space either within a row or
!      between the rows. icn(1) to icn(idisp(1)-1), are neither
!      required nor altered. if mc23a/ad been called, these will hold
!      information about the off-diagonal blocks.
!
! a    is a real/double precision array of length licn whose entries
!      idisp(2) to licn must be set by the user to the  values of the
!      non-zero entries of the matrix in the order indicated by  icn.
!      on output a will hold the lu factors of the matrix where again
!      the position in the matrix is determined by the corresponding
!      values in icn. a(1) to a(idisp(1)-1) are neither required nor altered.
!
! licn is an integer variable which must be set by the user to the
!      length of arrays icn and a. it must be big enough for a and icn
!      to hold all the non-zeros of l and u and leave some "elbow
!      room".  it is possible to calculate a minimum value for licn by
!      a preliminary run of ma30a/ad. the adequacy of the elbow room
!      can be judged by the size of the common block variable icncp. it
!      is not altered by ma30a/ad.
!
! lenr is an integer array of length n.  on input, lenr(i) should
!      equal the number of non-zeros in row i, i=1,...,n of the
!      diagonal blocks of p1*a*q1. on output, lenr(i) will equal the
!      total number of non-zeros in row i of l and row i of u.
!
! lenrl is an integer array of length n. on output from ma30a/ad,
!       lenrl(i) will hold the number of non-zeros in row i of l.
!
! idisp is an integer array of length 2. the user should set idisp(1)
!       to be the first available position in a/icn for the lu
!       decomposition while idisp(2) is set to the position in a/icn of
!       the first non-zero in the diagonal blocks of p1*a*q1. on output,
!       idisp(1) will be unaltered while idisp(2) will be set to the
!       position in a/icn of the last non-zero of the lu decomposition.
!
! ip    is an integer array of length n which holds a permutation of
!       the integers 1 to n.  on input to ma30a/ad, the absolute value of
!       ip(i) must be set to the row of a which is row i of p1*a*q1. a
!       negative value for ip(i) indicates that row i is at the end of a
!       diagonal block.  on output from ma30a/ad, ip(i) indicates the row
!       of a which is the i th row in paq. ip(i) will still be negative
!       for the last row of each block (except the last).
!
! iq    is an integer array of length n which again holds a
!       permutation of the integers 1 to n.  on input to ma30a/ad, iq(j)
!       must be set to the column of a which is column j of p1*a*q1. on
!       output from ma30a/ad, the absolute value of iq(j) indicates the
!       column of a which is the j th in paq.  for rows, i say, in which
!       structural or numerical singularity is detected iq(i) is negated.
!
! irn  is an integer array of length lirn used as workspace by ma30a/ad.
!
! lirn is an integer variable. it should be greater than the
!      largest number of non-zeros in a diagonal block of p1*a*q1 but
!      need not be as large as licn. it is the length of array irn and
!      should be large enough to hold the active part of any block,
!      plus some "elbow room", the  a posteriori  adequacy of which can
!      be estimated by examining the size of common block variable irncp.
!
! lenc,ifirst,lastr,nextr,lastc,nextc
!      are all integer arrays of length n which are used as workspace by
!      ma30a/ad.  if nsrch is set to a value less than or equal to n, then
!      arrays lastc and nextc are not referenced by ma30a/ad and so can be
!      dummied in the call to ma30a/ad.
!
! iptr,ipc are integer arrays of length n; used as workspace by ma30a/ad.
!
! u    is a real/double precision variable which should be set by the
!      user to a value between 0. and 1.0. if less than zero it is
!      reset to zero and if its value is 1.0 or greater it is reset to
!      0.9999 (0.999999999 in d version).  it determines the balance
!      between pivoting for sparsity and for stability, values near
!      zero emphasizing sparsity and values near one emphasizing
!      stability. we recommend u=0.1 as a posible first trial value.
!      the stability can be judged by a later call to mc24a/ad or by
!      setting lbig to .true.
!
! iflag is an integer variable. it will have a non-negative value if
!       ma30a/ad is successful. negative values indicate error
!       conditions while positive values indicate that the matrix has
!       been successfully decomposed but is singular. for each non-zero
!       value, an appropriate message is output on unit lp.  possible
!       non-zero values for iflag are
!  -1   the matrix is structually singular with rank given by irank in
!       common block ma30f/fd.
!  +1   if, however, the user wants the lu decomposition of a
!       structurally singular matrix and sets the common block variable
!       abort1 to .false., then, in the event of singularity and a
!       successful decomposition, iflag is returned with the value +1
!       and no message is output.
!  -2   the matrix is numerically singular (it may also be structually
!       singular) with estimated rank given by irank in common block ma30f/fd.
!  +2   the  user can choose to continue the decomposition even when a
!       zero pivot is encountered by setting common block variable
!       abort2 to .false.  if a singularity is encountered, iflag will
!       then return with a value of +2, and no message is output if the
!       decomposition has been completed successfully.
!  -3   lirn has not been large enough to continue with the
!       decomposition.  if the stage was zero then common block variable
!       minirn gives the length sufficient to start the decomposition on
!       this block.  for a successful decomposition on this block the user
!       should make lirn slightly (say about n/2) greater than this value.
!  -4   licn not large enough to continue with the decomposition.
!  -5   the decomposition has been completed but some of the lu factors
!       have been discarded to create enough room in a/icn to continue
!       the decomposition. the variable minicn in common block ma30f/fd
!       then gives the size that licn should be to enable the
!       factorization to be successful.  if the user sets common block
!       variable abort3 to .true., then the subroutine will exit
!       immediately instead of destroying any factors and continuing.
!  -6   both licn and lirn are too small. termination has been caused by
!       lack of space in irn (see error iflag= -3), but already some of
!       the lu factors in a/icn have been lost (see error iflag= -5).
!       minicn gives the minimum amount of space required in a/icn for
!       decomposition up to this point.
!
!
! declare
      logical          abort1,abort2,abort3,lbig
      integer          nn,licn,lirn,iptr(nn),pivot,pivend,dispc, &
                       oldpiv,oldend,pivrow,rowi,ipc(nn),idisp(2), &
                       colupd,icn(licn),lenr(nn),lenrl(nn),ip(nn), &
                       iq(nn),lenc(nn),irn(lirn),ifirst(nn),lastr(nn), &
                       nextr(nn),lastc(nn),nextc(nn),lpiv(10),lnpiv(10), &
                       msrch,nsrch,ndrop,kk,mapiv,manpiv,iavpiv,ianpiv, &
                       kountl,minirn,minicn,morei,irank,irncp,icncp, &
                       iflag,ibeg,iactiv,nzrow,num,nnm1,i,ilast,nblock, &
                       istart,irows,n,ising,lp,itop,ii,j1,jj,j,indrow, &
                       j2,ipos,nzcol,nzmin,nz,isw,isw1,jcost, &
                       isrch,ll,ijfir,idummy,kcost,ijpos,ipiv,jpiv,i1, &
                       i2,jpos,k,lc,nc,lr,nr,lenpiv,ijp1,nz2,l,lenpp, &
                       nzpc,iii,idrop,iend,iop,jnew,ifill,jdiff,jnpos, &
                       jmore,jend,jroom,jbeg,jzero,idispc,kdrop, &
                       ifir,jval,jzer,jcount,jdummy,jold
      double precision a(licn),u,au,umax,amax,zero,pivrat,pivr, &
                       tol,big,anew,aanew,scale
!
      common /ma30ed/  lp,abort1,abort2,abort3
      common /ma30fd/  irncp,icncp,irank,minirn,minicn
      common /ma30id/  tol,big,ndrop,nsrch,lbig
      common /lpivot/  lpiv,lnpiv,mapiv,manpiv,iavpiv,ianpiv,kountl
!
      data             umax/0.999999999d0/
      data             zero /0.0d0/
!
!
! formats
99992 format(1x,'to continue set lirn to at least',i8)
99993 format(1x,'at stage',i5,'in block',i5,'with first row',i5, &
                'and last row',i5)
99994 format(1x,'error return from ma30a/ad lirn and licn too small')
99995 format(1x,'error return from ma30a/ad ; lirn not big enough')
99996 format(1x,'error return from ma30a/ad ; licn not big enough')
99997 format(1x,'lu decomposition destroyed to create more space')
99998 format(1x,'error return from ma30a/ad;', &
                'matrix is numerically singular')
99999 format(1x,'error return from ma30a/ad;', &
                'matrix is structurally singular')
!
!
!
! initialize
      msrch = nsrch
      ndrop = 0
      do 1272 kk=1,10
       lnpiv(kk) = 0
       lpiv(kk)  = 0
1272  continue
      mapiv  = 0
      manpiv = 0
      iavpiv = 0
      ianpiv = 0
      kountl = 0
      minirn = 0
      minicn = idisp(1) - 1
      morei  = 0
      irank  = nn
      irncp  = 0
      icncp  = 0
      iflag  = 0
      u      = min(u,umax)
      u      = max(u,zero)
!
! ibeg is the position of the next pivot row after elimination step using it.
! iactiv is the position of the first entry in the active part of a/icn.
! nzrow is current number of non-zeros in active and unprocessed part of row
! file icn.
      ibeg   = idisp(1)
      iactiv = idisp(2)
      nzrow  = licn - iactiv + 1
      minicn = nzrow + minicn
!
! count the number of diagonal blocks and set up pointers to the beginnings of
! the rows. num is the number of diagonal blocks.
      num = 1
      iptr(1) = iactiv
      if (nn .eq. 1) go to 20
      nnm1 = nn - 1
      do 10 i=1,nnm1
       if (ip(i) .lt. 0) num = num + 1
       iptr(i+1) = iptr(i) + lenr(i)
10    continue
!
! ilast is the last row in the previous block.
20    ilast = 0
!
! lu decomposition of block nblock starts
! each pass on this loop performs lu decomp on one of the diagonal blocks.
      do 1000 nblock=1,num
       istart = ilast + 1
       do 30 irows=istart,nn
        if (ip(irows) .lt. 0) go to 40
30     continue
       irows = nn
40     ilast = irows
!
! n is the number of rows in the current block.
! istart is the index of the first row in the current block.
! ilast is the index of the last row in the current block.
! iactiv is the position of the first entry in the block.
! itop is the position of the last entry in the block.
       n = ilast - istart + 1
       if (n .ne. 1) go to 90
!
! code for dealing with 1x1 block.
       lenrl(ilast) = 0
       ising = istart
       if (lenr(ilast) .ne. 0) go to 50
!
! block is structurally singular.
       irank = irank - 1
       ising = -ising
       if (iflag .ne. 2 .and. iflag .ne. -5) iflag = 1
       if (.not.abort1) go to 80
       idisp(2) = iactiv
       iflag = -1
       if (lp .ne. 0) write (lp,99999)
       go to 1120
!
50     scale = abs(a(iactiv))
       if (scale .eq. zero) go to 60
       if (lbig) big = max(big,scale)
       go to 70
60     ising = -ising
       irank = irank - 1
       iptr(ilast) = 0
       if (iflag .ne. -5) iflag = 2
       if (.not.abort2) go to 70
       idisp(2) = iactiv
       iflag    = -2
       if (lp .ne. 0) write (lp,99998)
       go to 1120
70     a(ibeg)       = a(iactiv)
       icn(ibeg)     = icn(iactiv)
       iactiv        = iactiv + 1
       iptr(istart)  = 0
       ibeg          = ibeg + 1
       nzrow         = nzrow - 1
80     lastr(istart) = istart
       ipc(istart)   = -ising
       go to 1000
!
! non-trivial block.
90     itop = licn
       if (ilast .ne. nn) itop = iptr(ilast+1) - 1
!
! set up column oriented storage.
       do 100 i=istart,ilast
        lenrl(i) = 0
        lenc(i)  = 0
100    continue
       if (itop-iactiv .lt. lirn) go to 110
       minirn = itop - iactiv + 1
       pivot  = istart - 1
       go to 1100
!
! calculate column counts.
110    do 120 ii=iactiv,itop
        i       = icn(ii)
        lenc(i) = lenc(i) + 1
120    continue
!
! set up column pointers so that ipc(j) points to position after end of
! column j in column file.
       ipc(ilast) = lirn + 1
       j1         = istart + 1
       do 130 jj=j1,ilast
        j      = ilast - jj + j1 - 1
        ipc(j) = ipc(j+1) - lenc(j+1)
130    continue
       do 150 indrow=istart,ilast
        j1 = iptr(indrow)
        j2 = j1 + lenr(indrow) - 1
        if (j1 .gt. j2) go to 150
        do 140 jj=j1,j2
         j         = icn(jj)
         ipos      = ipc(j) - 1
         irn(ipos) = indrow
         ipc(j)    = ipos
140     continue
150    continue
!
! dispc is the lowest indexed active location in the column file.
       dispc  = ipc(istart)
       nzcol  = lirn - dispc + 1
       minirn = max0(nzcol,minirn)
       nzmin  = 1
!
! initialize array ifirst.  ifirst(i) = +/- k indicates that row/col k has i
! non-zeros.  if ifirst(i) = 0,there is no row or column with i non zeros.
       do 160 i=1,n
        ifirst(i) = 0
160    continue
!
! compute ordering of row and column counts. first run through columns (from
! column n to column 1).
       do 180 jj=istart,ilast
        j  = ilast - jj + istart
        nz = lenc(j)
        if (nz .ne. 0) go to 170
        ipc(j) = 0
        go to 180
170     if (nsrch .le. nn) go to 180
        isw        = ifirst(nz)
        ifirst(nz) = -j
        lastc(j)   = 0
        nextc(j)   = -isw
        isw1       = iabs(isw)
        if (isw .ne. 0) lastc(isw1) = j
180    continue
!
! now run through rows (again from n to 1).
       do 210 ii=istart,ilast
        i  = ilast - ii + istart
        nz = lenr(i)
        if (nz .ne. 0) go to 190
        iptr(i)  = 0
        lastr(i) = 0
        go to 210
190     isw        = ifirst(nz)
        ifirst(nz) = i
        if (isw .gt. 0) go to 200
        nextr(i) = 0
        lastr(i) = isw
        go to 210
200     nextr(i)   = isw
        lastr(i)   = lastr(isw)
        lastr(isw) = i
210    continue
!
!
! start of main elimination loop
!
! first find the pivot using markowitz criterion with stability control.
! jcost is the markowitz cost of the best pivot so far,.. this pivot is in
! row ipiv and column jpiv.
!
       do 980 pivot=istart,ilast
        nz2   = nzmin
        jcost = n*n
!
! examine rows/columns in order of ascending count.
        do 340 l=1,2
         pivrat = zero
         isrch  = 1
         ll     = l
!
! a pass with l equal to 2 is only performed in the case of singularity.
         do 330 nz=nz2,n
          if (jcost .le. (nz-1)**2) go to 420
          ijfir = ifirst(nz)
          if (ijfir) 230, 220, 240
220       if (ll .eq. 1) nzmin = nz + 1
          go to 330
230       ll    = 2
          ijfir = -ijfir
          go to 290
240       ll = 2
!
!  scan rows with nz non-zeros.
          do 270 idummy=1,n
           if (jcost .le. (nz-1)**2) go to 420
           if (isrch .gt. msrch) go to 420
           if (ijfir .eq. 0) go to 280
!
! row ijfir is now examined.
           i     = ijfir
           ijfir = nextr(i)
!
! first calculate multiplier threshold level.
           amax = zero
           j1   = iptr(i) + lenrl(i)
           j2   = iptr(i) + lenr(i) - 1
           do 250 jj=j1,j2
            amax = max(amax,abs(a(jj)))
250        continue
           au    = amax*u
           isrch = isrch + 1
!
! scan row for possible pivots
           do 260 jj=j1,j2
            if (abs(a(jj)) .le. au .and. l .eq. 1) go to 260
            j     = icn(jj)
            kcost = (nz-1)*(lenc(j)-1)
            if (kcost .gt. jcost) go to 260
            pivr = zero
            if (amax .ne. zero) pivr = abs(a(jj))/amax
            if (kcost .eq. jcost .and. (pivr .le. pivrat .or. &
                nsrch .gt. nn+1)) go to 260
!
! best pivot so far is found.
            jcost = kcost
            ijpos = jj
            ipiv  = i
            jpiv  = j
            if (msrch .gt. nn+1 .and. jcost .le. (nz-1)**2) go to 420
            pivrat = pivr
260        continue
270       continue
!
! columns with nz non-zeros now examined.
280       ijfir = ifirst(nz)
          ijfir = -lastr(ijfir)
290       if (jcost .le. nz*(nz-1)) go to 420
          if (msrch .le. nn) go to 330
          do 320 idummy=1,n
           if (ijfir .eq. 0) go to 330
           j     = ijfir
           ijfir = nextc(ijfir)
           i1    = ipc(j)
           i2    = i1 + nz - 1
!
! scan column j
           do 310 ii=i1,i2
            i     = irn(ii)
            kcost = (nz-1)*(lenr(i)-lenrl(i)-1)
            if (kcost .ge. jcost) go to 310
!
! pivot has best markowitz count so far ... now check its suitability on
! numeric grounds by examining the other non-zeros in its row.
            j1 = iptr(i) + lenrl(i)
            j2 = iptr(i) + lenr(i) - 1
!
! we need a stability check on singleton columns because of possible problems
! with underdetermined systems.
            amax = zero
            do 300 jj=j1,j2
             amax = max(amax,abs(a(jj)))
             if (icn(jj) .eq. j) jpos = jj
300         continue
            if (abs(a(jpos)) .le. amax*u .and. l .eq. 1) go to 310
            jcost = kcost
            ipiv  = i
            jpiv  = j
            ijpos = jpos
            if (amax .ne. zero) pivrat = abs(a(jpos))/amax
            if (jcost .le. nz*(nz-1)) go to 420
310        continue
320       continue
330      continue
!
! in the event of singularity; must make sure all rows and columns are tested.
! matrix is numerically or structurally singular; it will be diagnosed later.
         msrch = n
         irank = irank - 1
340     continue
!
! assign rest of rows and columns to ordering array. matrix is singular.
        if (iflag .ne. 2 .and. iflag .ne. -5) iflag = 1
        irank = irank - ilast + pivot + 1
        if (.not.abort1) go to 350
        idisp(2) = iactiv
        iflag = -1
        if (lp .ne. 0) write (lp,99999)
        go to 1120
350     k = pivot - 1
        do 390 i=istart,ilast
         if (lastr(i) .ne. 0) go to 390
         k        = k + 1
         lastr(i) = k
         if (lenrl(i) .eq. 0) go to 380
         minicn = max0(minicn,nzrow+ibeg-1+morei+lenrl(i))
         if (iactiv-ibeg .ge. lenrl(i)) go to 360
         call ma30dd(a,icn,iptr(istart),n,iactiv,itop,.true.)
!
! check now to see if ma30d/dd has created enough available space.
         if (iactiv-ibeg .ge. lenrl(i)) go to 360
!
! create more space by destroying previously created lu factors.
         morei = morei + ibeg - idisp(1)
         ibeg = idisp(1)
         if (lp .ne. 0) write (lp,99997)
         iflag = -5
         if (abort3) go to 1090
360      j1 = iptr(i)
         j2 = j1 + lenrl(i) - 1
         iptr(i) = 0
         do 370 jj=j1,j2
          a(ibeg)   = a(jj)
          icn(ibeg) = icn(jj)
          icn(jj)   = 0
          ibeg      = ibeg + 1
370      continue
         nzrow = nzrow - lenrl(i)
380      if (k .eq. ilast) go to 400
390     continue
400     k = pivot - 1
        do 410 i=istart,ilast
         if (ipc(i) .ne. 0) go to 410
         k      = k + 1
         ipc(i) = k
         if (k .eq. ilast) go to 990
410     continue
!
! the pivot has now been found in position (ipiv,jpiv) in location ijpos in
! row file. update column and row ordering arrays to correspond with removal
! of the active part of the matrix.
420     ising = pivot
        if (a(ijpos) .ne. zero) go to 430
!
! numerical singularity is recorded here.
        ising = -ising
        if (iflag .ne. -5) iflag = 2
        if (.not.abort2) go to 430
        idisp(2) = iactiv
        iflag = -2
        if (lp .ne. 0) write (lp,99998)
        go to 1120
430     oldpiv = iptr(ipiv) + lenrl(ipiv)
        oldend = iptr(ipiv) + lenr(ipiv) - 1
!
! changes to column ordering.
        if (nsrch .le. nn) go to 460
        colupd = nn + 1
        lenpp  = oldend-oldpiv+1
        if (lenpp .lt. 4) lpiv(1) = lpiv(1) + 1
        if (lenpp.ge.4 .and. lenpp.le.6) lpiv(2) = lpiv(2) + 1
        if (lenpp.ge.7 .and. lenpp.le.10) lpiv(3) = lpiv(3) + 1
        if (lenpp.ge.11 .and. lenpp.le.15) lpiv(4) = lpiv(4) + 1
        if (lenpp.ge.16 .and. lenpp.le.20) lpiv(5) = lpiv(5) + 1
        if (lenpp.ge.21 .and. lenpp.le.30) lpiv(6) = lpiv(6) + 1
        if (lenpp.ge.31 .and. lenpp.le.50) lpiv(7) = lpiv(7) + 1
        if (lenpp.ge.51 .and. lenpp.le.70) lpiv(8) = lpiv(8) + 1
        if (lenpp.ge.71 .and. lenpp.le.100) lpiv(9) = lpiv(9) + 1
        if (lenpp.ge.101) lpiv(10) = lpiv(10) + 1
        mapiv  = max0(mapiv,lenpp)
        iavpiv = iavpiv + lenpp
        do 450 jj=oldpiv,oldend
         j        = icn(jj)
         lc       = lastc(j)
         nc       = nextc(j)
         nextc(j) = -colupd
         if (jj .ne. ijpos) colupd = j
         if (nc .ne. 0) lastc(nc) = lc
         if (lc .eq. 0) go to 440
         nextc(lc) = nc
         go to 450
440      nz  = lenc(j)
         isw = ifirst(nz)
         if (isw .gt. 0) lastr(isw) = -nc
         if (isw .lt. 0) ifirst(nz) = -nc
450     continue
!
! changes to row ordering.
460     i1 = ipc(jpiv)
        i2 = i1 + lenc(jpiv) - 1
        do 480 ii=i1,i2
         i  = irn(ii)
         lr = lastr(i)
         nr = nextr(i)
         if (nr .ne. 0) lastr(nr) = lr
         if (lr .le. 0) go to 470
         nextr(lr) = nr
         go to 480
470      nz = lenr(i) - lenrl(i)
         if (nr .ne. 0) ifirst(nz) = nr
         if (nr .eq. 0) ifirst(nz) = lr
480     continue
!
! move pivot to position lenrl+1 in pivot row and move pivot row to the
! beginning of the available storage. the l part and the pivot in the old
! copy of the pivot row is nullified while, in the strictly upper triangular
! part, the column indices, j say, are overwritten by the corresponding
! entry of iq (iq(j)) and iq(j) is set to the negative of the displacement of
! the column index from the pivot entry.
        if (oldpiv .eq. ijpos) go to 490
        au          = a(oldpiv)
        a(oldpiv)   = a(ijpos)
        a(ijpos)    = au
        icn(ijpos)  = icn(oldpiv)
        icn(oldpiv) = jpiv
!
! check if there is space available in a/icn to hold new copy of pivot row.
490     minicn = max0(minicn,nzrow+ibeg-1+morei+lenr(ipiv))
        if (iactiv-ibeg .ge. lenr(ipiv)) go to 500
        call ma30dd(a,icn,iptr(istart),n,iactiv,itop,.true.)
        oldpiv = iptr(ipiv) + lenrl(ipiv)
        oldend = iptr(ipiv) + lenr(ipiv) - 1
!
! check now to see if ma30d/dd has created enough available space.
        if (iactiv-ibeg .ge. lenr(ipiv)) go to 500
!
! create more space by destroying previously created lu factors.
        morei = morei + ibeg - idisp(1)
        ibeg = idisp(1)
        if (lp .ne. 0) write (lp,99997)
        iflag = -5
        if (abort3) go to 1090
        if (iactiv-ibeg .ge. lenr(ipiv)) go to 500
!
! there is still not enough room in a/icn.
        iflag = -4
        go to 1090
!
! copy pivot row and set up iq array.
500     ijpos = 0
        j1    = iptr(ipiv)
        do 530 jj=j1,oldend
         a(ibeg)   = a(jj)
         icn(ibeg) = icn(jj)
         if (ijpos .ne. 0) go to 510
         if (icn(jj) .eq. jpiv) ijpos = ibeg
         icn(jj) = 0
         go to 520
510      k       = ibeg - ijpos
         j       = icn(jj)
         icn(jj) = iq(j)
         iq(j)   = -k
520      ibeg    = ibeg + 1
530     continue
!
        ijp1       = ijpos + 1
        pivend     = ibeg - 1
        lenpiv     = pivend - ijpos
        nzrow      = nzrow - lenrl(ipiv) - 1
        iptr(ipiv) = oldpiv + 1
        if (lenpiv .eq. 0) iptr(ipiv) = 0
!
! remove pivot row (including pivot) from column oriented file.
        do 560 jj=ijpos,pivend
         j       = icn(jj)
         i1      = ipc(j)
         lenc(j) = lenc(j) - 1
!
! i2 is last position in new column.
         i2 = ipc(j) + lenc(j) - 1
         if (i2 .lt. i1) go to 550
         do 540 ii=i1,i2
          if (irn(ii) .ne. ipiv) go to 540
          irn(ii) = irn(i2+1)
          go to 550
540      continue
550      irn(i2+1) = 0
560     continue
        nzcol = nzcol - lenpiv - 1
!
! go down the pivot column and for each row with a non-zero add the
! appropriate multiple of the pivot row to it. we loop on the number of
! non-zeros in the pivot column since ma30d/dd may change its actual position.
        nzpc = lenc(jpiv)
        if (nzpc .eq. 0) go to 900
        do 840 iii=1,nzpc
         ii = ipc(jpiv) + iii - 1
         i  = irn(ii)
!
! search row i for non-zero to be eliminated, calculate multiplier, and place
! it in position lenrl+1 in its row. idrop is the number of non-zero entries
! dropped from row i because these fall beneath tolerance level.
         idrop = 0
         j1    = iptr(i) + lenrl(i)
         iend  = iptr(i) + lenr(i) - 1
         do 570 jj=j1,iend
          if (icn(jj) .ne. jpiv) go to 570
!
! if pivot is zero, rest of column is and so multiplier is zero.
          au = zero
          if (a(ijpos) .ne. zero) au = -a(jj)/a(ijpos)
          if (lbig) big = max(big,abs(au))
          a(jj)    = a(j1)
          a(j1)    = au
          icn(jj)  = icn(j1)
          icn(j1)  = jpiv
          lenrl(i) = lenrl(i) + 1
          go to 580
570      continue
!
! jump if pivot row is a singleton.
580      if (lenpiv .eq. 0) go to 840
!
! now perform necessary operations on rest of non-pivot row i.
         rowi = j1 + 1
         iop  = 0
!
! jump if all the pivot row causes fill-in.
         if (rowi .gt. iend) go to 650
!
! perform operations on current non-zeros in row i. innermost loop.
         lenpp = iend-rowi+1
         if (lenpp .lt. 4) lnpiv(1) = lnpiv(1) + 1
         if (lenpp.ge.4 .and. lenpp.le.6) lnpiv(2) = lnpiv(2) + 1
         if (lenpp.ge.7 .and. lenpp.le.10) lnpiv(3) = lnpiv(3) + 1
         if (lenpp.ge.11 .and. lenpp.le.15) lnpiv(4) = lnpiv(4) + 1
         if (lenpp.ge.16 .and. lenpp.le.20) lnpiv(5) = lnpiv(5) + 1
         if (lenpp.ge.21 .and. lenpp.le.30) lnpiv(6) = lnpiv(6) + 1
         if (lenpp.ge.31 .and. lenpp.le.50) lnpiv(7) = lnpiv(7) + 1
         if (lenpp.ge.51 .and. lenpp.le.70) lnpiv(8) = lnpiv(8) + 1
         if (lenpp.ge.71 .and. lenpp.le.100) lnpiv(9) = lnpiv(9) + 1
         if (lenpp.ge.101) lnpiv(10) = lnpiv(10) + 1
         manpiv = max0(manpiv,lenpp)
         ianpiv = ianpiv + lenpp
         kountl = kountl + 1
         do 590 jj=rowi,iend
          j = icn(jj)
          if (iq(j) .gt. 0) go to 590
          iop    = iop + 1
          pivrow = ijpos - iq(j)
          a(jj)  = a(jj) + au*a(pivrow)
          if (lbig) big = max(abs(a(jj)),big)
          icn(pivrow) = -icn(pivrow)
          if (abs(a(jj)) .lt. tol) idrop = idrop + 1
590      continue
!
! jump if no non-zeros in non-pivot row have been removed because these are
! beneath the drop-tolerance  tol.
         if (idrop .eq. 0) go to 650
!
! run through non-pivot row compressing row so that only non-zeros greater
! than tol are stored. all non-zeros less than tol are also removed from the
! column structure.
         jnew = rowi
         do 630 jj=rowi,iend
          if (abs(a(jj)) .lt. tol) go to 600
          a(jnew)   = a(jj)
          icn(jnew) = icn(jj)
          jnew      = jnew + 1
          go to 630
!
! remove non-zero entry from column structure.
600       j = icn(jj)
          i1 = ipc(j)
          i2 = i1 + lenc(j) - 1
          do 610 ii=i1,i2
           if (irn(ii) .eq. i) go to 620
610       continue
620       irn(ii) = irn(i2)
          irn(i2) = 0
          lenc(j) = lenc(j) - 1
          if (nsrch .le. nn) go to 630
!
! remove column from column chain and place in update chain.
          if (nextc(j) .lt. 0) go to 630
!
! jump if column already in update chain.
          lc       = lastc(j)
          nc       = nextc(j)
          nextc(j) = -colupd
          colupd   = j
          if (nc .ne. 0) lastc(nc) = lc
          if (lc .eq. 0) go to 622
          nextc(lc) = nc
          go to 630
622       nz = lenc(j) + 1
          isw = ifirst(nz)
          if (isw .gt. 0) lastr(isw) = -nc
          if (isw .lt. 0) ifirst(nz) = -nc
630      continue
         do 640 jj=jnew,iend
          icn(jj) = 0
640      continue
!
! the value of idrop might be different from that calculated earlier because,
! we may have dropped some non-zeros which were not modified by the pivot row.
         idrop   = iend + 1 - jnew
         iend    = jnew - 1
         lenr(i) = lenr(i) - idrop
         nzrow   = nzrow - idrop
         nzcol   = nzcol - idrop
         ndrop   = ndrop + idrop
650      ifill   = lenpiv - iop
!
! jump is if there is no fill-in.
         if (ifill .eq. 0) go to 750
!
! now for the fill-in.
         minicn = max0(minicn,morei+ibeg-1+nzrow+ifill+lenr(i))
!
! see if there is room for fill-in. get maximum space for row i in situ.
         do 660 jdiff=1,ifill
          jnpos = iend + jdiff
          if (jnpos .gt. licn) go to 670
          if (icn(jnpos) .ne. 0) go to 670
660      continue
!
! there is room for all the fill-in after the end of the row so it can be
! left in situ. next available space for fill-in.
         iend = iend + 1
         go to 750
!
! jmore spaces for fill-in are required in front of row.
670      jmore = ifill - jdiff + 1
         i1    = iptr(i)
!
! look in front of the row to see if there is space for rest of the fill-in.
         do 680 jdiff=1,jmore
          jnpos = i1 - jdiff
          if (jnpos .lt. iactiv) go to 690
          if (icn(jnpos) .ne. 0) go to 700
680      continue
690      jnpos = i1 - jmore
         go to 710
!
! whole row must be moved to the beginning of available storage.
700      jnpos = iactiv - lenr(i) - ifill
!
! jump if there is space immediately available for the shifted row.
710      if (jnpos .ge. ibeg) go to 730
         call ma30dd(a,icn,iptr(istart),n,iactiv,itop,.true.)
         i1    = iptr(i)
         iend  = i1 + lenr(i) - 1
         jnpos = iactiv - lenr(i) - ifill
         if (jnpos .ge. ibeg) go to 730
!
! no space available; try to create some by trashing previous lu decomposition.
         morei = morei + ibeg - idisp(1) - lenpiv - 1
         if (lp .ne. 0) write (lp,99997)
         iflag = -5
         if (abort3) go to 1090
!
! keep record of current pivot row.
         ibeg      = idisp(1)
         icn(ibeg) = jpiv
         a(ibeg)   = a(ijpos)
         ijpos     = ibeg
         do 720 jj=ijp1,pivend
          ibeg      = ibeg + 1
          a(ibeg)   = a(jj)
          icn(ibeg) = icn(jj)
  720    continue
         ijp1   = ijpos + 1
         pivend = ibeg
         ibeg   = ibeg + 1
         if (jnpos .ge. ibeg) go to 730
!
! this still does not give enough room.
         iflag = -4
         go to 1090
730      iactiv = min0(iactiv,jnpos)
!
! move non-pivot row i.
         iptr(i) = jnpos
         do 740 jj=i1,iend
          a(jnpos)   = a(jj)
          icn(jnpos) = icn(jj)
          jnpos      = jnpos + 1
          icn(jj)    = 0
740      continue
!
! first new available space.
         iend  = jnpos
750      nzrow = nzrow + ifill
!
! innermost fill-in loop which also resets icn.
         idrop = 0
         do 830 jj=ijp1,pivend
          j = icn(jj)
          if (j .lt. 0) go to 820
          anew  = au*a(jj)
          aanew = abs(anew)
          if (aanew .ge. tol) go to 760
          idrop  = idrop + 1
          ndrop  = ndrop + 1
          nzrow  = nzrow - 1
          minicn = minicn - 1
          ifill  = ifill - 1
          go to 830
760       if (lbig) big = max(aanew,big)
          a(iend)   = anew
          icn(iend) = j
          iend      = iend + 1
!
! put new entry in column file.
          minirn = max0(minirn,nzcol+lenc(j)+1)
          jend   = ipc(j) + lenc(j)
          jroom  = nzpc - iii + 1 + lenc(j)
          if (jend .gt. lirn) go to 770
          if (irn(jend) .eq. 0) go to 810
770       if (jroom .lt. dispc) go to 780
!
! compress column file to obtain space for new copy of column.
          call ma30dd(a,irn,ipc(istart),n,dispc,lirn,.false.)
          if (jroom .lt. dispc) go to 780
          jroom = dispc - 1
          if (jroom .ge. lenc(j)+1) go to 780
!
! column file is not large enough.
          go to 1100
!
! copy column to beginning of file.
780       jbeg   = ipc(j)
          jend   = ipc(j) + lenc(j) - 1
          jzero  = dispc - 1
          dispc  = dispc - jroom
          idispc = dispc
          do 790 ii=jbeg,jend
           irn(idispc) = irn(ii)
           irn(ii) = 0
           idispc  = idispc + 1
790       continue
          ipc(j) = dispc
          jend   = idispc
          do 800 ii=jend,jzero
           irn(ii) = 0
800       continue
810       irn(jend) = i
          nzcol     = nzcol + 1
          lenc(j)   = lenc(j) + 1
!
! end of adjustment to column file.
          go to 830
!
820       icn(jj) = -j
830      continue
         if (idrop .eq. 0) go to 834
         do 832 kdrop=1,idrop
          icn(iend) = 0
          iend = iend + 1
832      continue
834      lenr(i) = lenr(i) + ifill
!
! end of scan of pivot column.
840     continue
!
!
! remove pivot column from column oriented storage; update row ordering arrays.
        i1 = ipc(jpiv)
        i2 = ipc(jpiv) + lenc(jpiv) - 1
        nzcol = nzcol - lenc(jpiv)
        do 890 ii=i1,i2
         i       = irn(ii)
         irn(ii) = 0
         nz      = lenr(i) - lenrl(i)
         if (nz .ne. 0) go to 850
         lastr(i) = 0
         go to 890
850      ifir       = ifirst(nz)
         ifirst(nz) = i
         if (ifir) 860, 880, 870
860      lastr(i) = ifir
         nextr(i) = 0
         go to 890
870      lastr(i)    = lastr(ifir)
         nextr(i)    = ifir
         lastr(ifir) = i
         go to 890
880      lastr(i) = 0
         nextr(i) = 0
         nzmin    = min0(nzmin,nz)
890     continue
!
! restore iq and nullify u part of old pivot row. record the column
! permutation in lastc(jpiv) and the row permutation in lastr(ipiv).
900     ipc(jpiv)   = -ising
        lastr(ipiv) = pivot
        if (lenpiv .eq. 0) go to 980
        nzrow      = nzrow - lenpiv
        jval       = ijp1
        jzer       = iptr(ipiv)
        iptr(ipiv) = 0
        do 910 jcount=1,lenpiv
         j         = icn(jval)
         iq(j)     = icn(jzer)
         icn(jzer) = 0
         jval      = jval + 1
         jzer      = jzer + 1
910     continue
!
! adjust column ordering arrays.
        if (nsrch .gt. nn) go to 920
        do 916 jj=ijp1,pivend
         j  = icn(jj)
         nz = lenc(j)
         if (nz .ne. 0) go to 914
         ipc(j) = 0
         go to 916
914      nzmin = min0(nzmin,nz)
916     continue
        go to 980
920     jj = colupd
        do 970 jdummy=1,nn
         j = jj
         if (j .eq. nn+1) go to 980
         jj = -nextc(j)
         nz = lenc(j)
         if (nz .ne. 0) go to 924
         ipc(j) = 0
         go to 970
924      ifir     = ifirst(nz)
         lastc(j) = 0
         if (ifir) 930, 940, 950
930      ifirst(nz)  = -j
         ifir        = -ifir
         lastc(ifir) = j
         nextc(j)    = ifir
         go to 970
940      ifirst(nz) = -j
         nextc(j)   = 0
         go to 960
950      lc          = -lastr(ifir)
         lastr(ifir) = -j
         nextc(j)    = lc
         if (lc .ne. 0) lastc(lc) = j
960      nzmin = min0(nzmin,nz)
970     continue
980    continue
!
! that was the end of main elimination loop
!
!
! reset iactiv to point to the beginning of the next block.
990    if (ilast .ne. nn) iactiv = iptr(ilast+1)
1000  continue
!
! that was the end of deomposition of block
!
!
! record singularity (if any) in iq array.
      if (irank .eq. nn) go to 1020
      do 1010 i=1,nn
       if (ipc(i) .lt. 0) go to 1010
       ising     = ipc(i)
       iq(ising) = -iq(ising)
       ipc(i)    = -ising
1010  continue
!
!
! run through lu decomposition changing column indices to that of new order
! and permuting lenr and lenrl arrays according to pivot permutations.
1020  istart = idisp(1)
      iend = ibeg - 1
      if (iend .lt. istart) go to 1040
      do 1030 jj=istart,iend
       jold    = icn(jj)
       icn(jj) = -ipc(jold)
1030  continue
1040  do 1050 ii=1,nn
       i        = lastr(ii)
       nextr(i) = lenr(ii)
       iptr(i)  = lenrl(ii)
1050  continue
      do 1060 i=1,nn
       lenrl(i) = iptr(i)
       lenr(i)  = nextr(i)
1060  continue
!
! update permutation arrays ip and iq.
      do 1070 ii=1,nn
       i        = lastr(ii)
       j        = -ipc(ii)
       nextr(i) = iabs(ip(ii)+0)
       iptr(j)  = iabs(iq(ii)+0)
1070  continue
      do 1080 i=1,nn
       if (ip(i) .lt. 0) nextr(i) = -nextr(i)
       ip(i) = nextr(i)
       if (iq(i) .lt. 0) iptr(i) = -iptr(i)
       iq(i) = iptr(i)
1080  continue
      ip(nn)   = iabs(ip(nn)+0)
      idisp(2) = iend
      go to 1120
!
!
! error returns
1090  idisp(2) = iactiv
      if (lp .eq. 0) go to 1120
      write (lp,99996)
      go to 1110
1100  if (iflag .eq. -5) iflag = -6
      if (iflag .ne. -6) iflag = -3
      idisp(2) = iactiv
      if (lp .eq. 0) go to 1120
      if (iflag .eq. -3) write (lp,99995)
      if (iflag .eq. -6) write (lp,99994)
1110  pivot = pivot - istart + 1
      write (lp,99993) pivot,nblock,istart,ilast
      if (pivot .eq. 0) write (lp,99992) minirn
1120  return
      end
!
!
!
!
!
      subroutine ma30bd(n,icn,a,licn,lenr,lenrl,idisp,ip,iq,w,iw,iflag)
      include 'implno.dek'
!
! ma30b/bd performs the lu decomposition of the diagonal blocks of a new
! matrix paq of the same sparsity pattern, using information from a previous
! call to ma30a/ad. the entries of the input matrix  must already be in their
! final positions in the lu decomposition structure.  this routine executes
! about five times faster than ma30a/ad.
!
! parameters (see also ma30ad):
! n   is an integer variable set to the order of the matrix.
!
! icn is an integer array of length licn. it should be unchanged
!     since the last call to ma30a/ad. it is not altered by ma30b/bd.
!
! a   is a real/double precision array of length licn the user must set
!     entries idisp(1) to idisp(2) to contain the entries in the
!     diagonal blocks of the matrix paq whose column numbers are held
!     in icn, using corresponding positions. note that some zeros may
!     need to be held explicitly. on output entries idisp(1) to
!     idisp(2) of array a contain the lu decomposition of the diagonal
!     blocks of paq. entries a(1) to a(idisp(1)-1) are neither
!     required nor altered by ma30b/bd.
!
! licn is an integer variable which must be set by the user to the
!      length of arrays a and icn. it is not altered by ma30b/bd.
!
! lenr,lenrl are integer arrays of length n. they should be
!     unchanged since the last call to ma30a/ad. not altered by ma30b/bd.
!
! idisp is an integer array of length 2. it should be unchanged since
!       the last call to ma30a/ad. it is not altered by ma30b/bd.
!
! ip,iq are integer arrays of length n. they should be unchanged
!       since the last call to ma30a/ad. not altered by ma30b/bd.
!
! w   is a array of length n which is used as workspace by ma30b/bd.
!
! iw  is an integer array of length n which is used as workspace by ma30b/bd.
!
! iflag  is an integer variable. on output from ma30b/bd, iflag has
!        the value zero if the factorization was successful, has the
!        value i if pivot i was very small and has the value -i if an
!        unexpected singularity was detected at stage i of the decomposition.
!
! declare
      logical          abort1,abort2,abort3,stab,lbig
      integer          n,licn,iflag,iw(n),idisp(2),pivpos,icn(licn), &
                       lenr(n),lenrl(n),ip(n),iq(n),ndrop,nsrch,ising, &
                       i,istart,ifin,ilend,j,ipivj,jfin,jay,jayjay, &
                       jj,lp
      double precision a(licn),w(n),au,eps,rowmax,zero,one,rmin,tol,big
!
      common /ma30ed/  lp,abort1,abort2,abort3
      common /ma30id/  tol,big,ndrop,nsrch,lbig
      common /ma30gd/  eps,rmin
      data             zero /0.0d0/, one /1.0d0/
!
! formats
99999 format(1x,'error return from ma30b/bd singularity in row',i8)
!
! initialize
      stab  = eps .le. one
      rmin  = eps
      ising = 0
      iflag = 0
      do 10 i=1,n
       w(i) = zero
10    continue
!
! set up pointers to the beginning of the rows.
      iw(1) = idisp(1)
      if (n .eq. 1) go to 25
      do 20 i=2,n
       iw(i) = iw(i-1) + lenr(i-1)
20    continue
!
! start  of main loop
! at step i, row i of a is transformed to row i of l/u by adding appropriate
! multiples of rows 1 to i-1. using row-gauss elimination.
! istart is beginning of row i of a and row i of l.
! ifin is end of row i of a and row i of u.
! ilend is end of row i of l.
25    do 160 i=1,n
       istart = iw(i)
       ifin   = istart + lenr(i) - 1
       ilend  = istart + lenrl(i) - 1
       if (istart .gt. ilend) go to 90
!
! load row i of a into vector w.
       do 30 jj=istart,ifin
        j    = icn(jj)
        w(j) = a(jj)
30     continue
!
! add multiples of appropriate rows of  i to i-1  to row i.
! ipivj is position of pivot in row j.
       do 70 jj=istart,ilend
        j     = icn(jj)
        ipivj = iw(j) + lenrl(j)
        au    = -w(j)/a(ipivj)
        if (lbig) big = max(abs(au),big)
        w(j) = au
!
! au * row j (u part) is added to row i.
        ipivj = ipivj + 1
        jfin  = iw(j) + lenr(j) - 1
        if (ipivj .gt. jfin) go to 70
!
! innermost loop.
        if (lbig) go to 50
        do 40 jayjay=ipivj,jfin
         jay    = icn(jayjay)
         w(jay) = w(jay) + au*a(jayjay)
40      continue
        go to 70
50      do 60 jayjay=ipivj,jfin
         jay    = icn(jayjay)
         w(jay) = w(jay) + au*a(jayjay)
         big    = max(abs(w(jay)),big)
60      continue
70     continue
!
! reload w back into a (now l/u)
       do 80 jj=istart,ifin
        j     = icn(jj)
        a(jj) = w(j)
        w(j)  = zero
80     continue
!
! now perform the stability checks.
90     pivpos = ilend + 1
       if (iq(i) .gt. 0) go to 140
!
! matrix had singularity at this point in ma30a/ad.
! is it the first such pivot in current block ?
       if (ising .eq. 0) ising = i
!
! does current matrix have a singularity in the same place ?
       if (pivpos .gt. ifin) go to 100
       if (a(pivpos) .ne. zero) go to 170
!
! it does .. so set ising if it is not the end of the current block
! check to see that appropriate part of l/u is zero or null.
100    if (istart .gt. ifin) go to 120
       do 110 jj=istart,ifin
        if (icn(jj) .lt. ising) go to 110
        if (a(jj) .ne. zero) go to 170
110    continue
120    if (pivpos .le. ifin) a(pivpos) = one
       if (ip(i) .gt. 0 .and. i .ne. n) go to 160
!
! end of current block ... reset zero pivots and ising.
       do 130 j=ising,i
        if ((lenr(j)-lenrl(j)) .eq. 0) go to 130
        jj    = iw(j) + lenrl(j)
        a(jj) = zero
130    continue
       ising = 0
       go to 160
!
! matrix had non-zero pivot in ma30a/ad at this stage.
140    if (pivpos .gt. ifin) go to 170
       if (a(pivpos) .eq. zero) go to 170
       if (.not.stab) go to 160
       rowmax = zero
       do 150 jj=pivpos,ifin
        rowmax = max(rowmax,abs(a(jj)))
150    continue
       if (abs(a(pivpos))/rowmax .ge. rmin) go to 160
       iflag = i
       rmin  = abs(a(pivpos))/rowmax
160   continue
      go to 180
!
! error return
170   if (lp .ne. 0) write (lp,99999) i
      iflag = -i
180   return
      end
!
!
!
!
!
      subroutine ma30cd(n,icn,a,licn,lenr,lenrl,lenoff,idisp,ip, &
                        iq,x,w,mtype)
      include 'implno.dek'
!
!
! ma30c/cd uses the factors produced by ma30a/ad or ma30b/bd to solve
! ax=b or a transpose x=b when the matrix p1*a*q1 (paq) is block lower
! triangular (including the case of only one diagonal block).
!
! parameters:
! n  is an integer variable set to the order of the matrix. it is not
!    altered by the subroutine.
!
! icn is an integer array of length licn. entries idisp(1) to
!     idisp(2) should be unchanged since the last call to ma30a/ad. if
!     the matrix has more than one diagonal block, then column indices
!     corresponding to non-zeros in sub-diagonal blocks of paq must
!     appear in positions 1 to idisp(1)-1. for the same row those
!     entries must be contiguous, with those in row i preceding those
!     in row i+1 (i=1,...,n-1) and no wasted space between rows.
!     entries may be in any order within each row. not altered by ma30c/cd.
!
! a  is a real/double precision array of length licn.  entries
!    idisp(1) to idisp(2) should be unchanged since the last call to
!    ma30a/ad or ma30b/bd.  if the matrix has more than one diagonal
!    block, then the values of the non-zeros in sub-diagonal blocks
!    must be in positions 1 to idisp(1)-1 in the order given by icn.
!    it is not altered by ma30c/cd.
!
! licn  is an integer variable set to the size of arrays icn and a.
!       it is not altered by ma30c/cd.
!
! lenr,lenrl are integer arrays of length n which should be
!      unchanged since the last call to ma30a/ad. not altered by ma30c/cd.
!
! lenoff  is an integer array of length n. if the matrix paq (or
!         p1*a*q1) has more than one diagonal block, then lenoff(i),
!         i=1,...,n should be set to the number of non-zeros in row i of
!         the matrix paq which are in sub-diagonal blocks.  if there is
!         only one diagonal block then lenoff(1) may be set to -1, in
!         which case the other entries of lenoff are never accessed. it is
!         not altered by ma30c/cd.
!
! idisp  is an integer array of length 2 which should be unchanged
!        since the last call to ma30a/ad. it is not altered by ma30c/cd.
!
!! ip,iq are integer arrays of length n which should be unchanged
!        since the last call to ma30a/ad. they are not altered by ma30c/cd.
!
! x   is a real/double precision array of length n. it must be set by
!     the user to the values of the right hand side vector b for the
!     equations being solved.  on exit from ma30c/cd it will be equal
!     to the solution x required.
!
! w  is a real/double precision array of length n which is used as
!    workspace by ma30c/cd.
!
! mtype is an integer variable which must be set by the user. if
!      mtype=1, then the solution to the system ax=b is returned; any
!      other value for mtype will return the solution to the system a
!      transpose x=b. it is not altered by ma30c/cd.
!
! declare
      logical          neg,nobloc
      integer          n,licn,idisp(2),icn(licn),lenr(n),lenrl(n), &
                       lenoff(n),ip(n),iq(n),mtype,ii,i,lt,ifirst, &
                       iblock,ltend,jj,j,iend,j1,ib,iii,j2,jpiv, &
                       jpivp1,ilast,iblend,numblk,k,j3,iback,lj2,lj1
      double precision a(licn),x(n),w(n),wii,wi,resid,zero
      common /ma30hd/  resid
      data zero       /0.0d0/
!
! final value of resid is the max residual for inconsistent set of equations.
      resid = zero
!
! nobloc is .true. if subroutine block has been used previously and is .false.
! otherwise.  the value .false. means that lenoff will not be subsequently
! accessed.
      nobloc = lenoff(1) .lt. 0
      if (mtype .ne. 1) go to 140
!
! now solve   a * x = b. neg is used to indicate when the last row in a block
! has been reached.  it is then set to true whereafter backsubstitution is
! performed on the block.
      neg = .false.
!
! ip(n) is negated so that the last row of the last block can be recognised.
! it is reset to its positive value on exit.
      ip(n) = -ip(n)
!
! preorder vector ... w(i) = x(ip(i))
      do 10 ii=1,n
       i     = ip(ii)
       i     = iabs(i)
       w(ii) = x(i)
10    continue
!
! lt is the position of first non-zero in current row of off-diagonal blocks.
! ifirst holds the index of the first row in the current block.
! iblock holds the position of the first non-zero in the current row
! of the lu decomposition of the diagonal blocks.
      lt     = 1
      ifirst = 1
      iblock = idisp(1)
!
! if i is not the last row of a block, then a pass through this loop adds the
! inner product of row i of the off-diagonal blocks and w to w and performs
! forward elimination using row i of the lu decomposition.   if i is the last
! row of a block then, after performing these aforementioned operations,
! backsubstitution is performed using the rows of the block.
      do 120 i=1,n
       wi = w(i)
       if (nobloc) go to 30
       if (lenoff(i) .eq. 0) go to 30
!
! operations using lower triangular blocks.
! ltend is the end of row i in the off-diagonal blocks.
       ltend = lt + lenoff(i) - 1
       do 20 jj=lt,ltend
        j  = icn(jj)
        wi = wi - a(jj)*w(j)
20     continue
!
! lt is set the beginning of the next off-diagonal row.
! set neg to .true. if we are on the last row of the block.
       lt = ltend + 1
30     if (ip(i) .lt. 0) neg = .true.
       if (lenrl(i) .eq. 0) go to 50
!
! forward elimination phase.
! iend is the end of the l part of row i in the lu decomposition.
       iend = iblock + lenrl(i) - 1
       do 40 jj=iblock,iend
        j  = icn(jj)
        wi = wi + a(jj)*w(j)
40     continue
!
! iblock is adjusted to point to the start of the next row.
50     iblock = iblock + lenr(i)
       w(i)   = wi
       if (.not.neg) go to 120
!
! back substitution phase.
! j1 is position in a/icn after end of block beginning in row ifirst
! and ending in row i.
       j1 = iblock
!
! are there any singularities in this block?  if not, continue
       ib = i
       if (iq(i) .gt. 0) go to 70
       do 60 iii=ifirst,i
        ib = i - iii + ifirst
        if (iq(ib) .gt. 0) go to 70
        j1    = j1 - lenr(ib)
        resid = max(resid,abs(w(ib)))
        w(ib) = zero
60     continue
!
! entire block is singular.
       go to 110
!
!
! each pass through this loop performs the back-substitution
! operations for a single row, starting at the end of the block and
! working through it in reverse order.
! j2 is end of row ii. j1 is beginning of row ii. jpiv is the position of the
! pivot in row ii. jump out if row ii of u has no non-zeros.
70     do 100 iii=ifirst,ib
        ii     = ib - iii + ifirst
        j2     = j1 - 1
        j1     = j1 - lenr(ii)
        jpiv   = j1 + lenrl(ii)
        jpivp1 = jpiv + 1
        if (j2 .lt. jpivp1) go to 90
        wii = w(ii)
        do 80 jj=jpivp1,j2
         j   = icn(jj)
         wii = wii - a(jj)*w(j)
80      continue
        w(ii) = wii
90      w(ii) = w(ii)/a(jpiv)
100    continue
110    ifirst = i + 1
       neg    = .false.
120   continue
!
! reorder solution vector ... x(i) = w(iqinverse(i))
      do 130 ii=1,n
       i    = iq(ii)
       i    = iabs(i)
       x(i) = w(ii)
130   continue
      ip(n) = -ip(n)
      go to 320
!
!
! now solve  atranspose * x = b. preorder vector ... w(i)=x(iq(i))
140   do 150 ii=1,n
       i     = iq(ii)
       i     = iabs(i)
       w(ii) = x(i)
150   continue
!
! lj1 points to the beginning the current row in the off-diagonal blocks.
! iblock is initialized to point to beginning of block after the last one
! ilast is the last row in the current block.
! iblend points to the position after the last non-zero in the current block.
      lj1    = idisp(1)
      iblock = idisp(2) + 1
      ilast  = n
      iblend = iblock
!
! each pass through this loop operates with one diagonal block and
! the off-diagonal part of the matrix corresponding to the rows
! of this block.  the blocks are taken in reverse order and the
! number of times the loop is entered is min(n,no. blocks+1).
      do 290 numblk=1,n
       if (ilast .eq. 0) go to 300
       iblock = iblock - lenr(ilast)
!
! this loop finds the index of the first row in the current block. it is
! first and iblock is set to the position of the beginning of this first row.
       do 160 k=1,n
        ii = ilast - k
        if (ii .eq. 0) go to 170
        if (ip(ii) .lt. 0) go to 170
        iblock = iblock - lenr(ii)
160    continue
170    ifirst = ii + 1
!
! j1 points to the position of the beginning of row i (lt part) or pivot
       j1 = iblock
!
! forward elimination. each pass through this loop performs the operations
! for one row of the block.  if the corresponding entry of w is zero then the
! operations can be avoided.
       do 210 i=ifirst,ilast
        if (w(i) .eq. zero) go to 200
!
!  jump if row i singular.
        if (iq(i) .lt. 0) go to 220
!
! j2 first points to the pivot in row i and then is made to point to the
! first non-zero in the u transpose part of the row.
! j3 points to the end of row i.
        j2 = j1 + lenrl(i)
        wi = w(i)/a(j2)
        if (lenr(i)-lenrl(i) .eq. 1) go to 190
        j2 = j2 + 1
        j3 = j1 + lenr(i) - 1
        do 180 jj=j2,j3
         j    = icn(jj)
         w(j) = w(j) - a(jj)*wi
180     continue
190     w(i) = wi
200     j1   = j1 + lenr(i)
210    continue
       go to 240
!
! deals with rest of block which is singular.
220    do 230 ii=i,ilast
        resid = max(resid,abs(w(ii)))
        w(ii) = zero
230    continue
!
! back substitution. this loop does the back substitution on the rows of the
! block in the reverse order doing it simultaneously on the l transpose part
! of the diagonal blocks and the off-diagonal blocks.
! j1 points to the beginning of row i.
! j2 points to the end of the l transpose part of row i.
240    j1 = iblend
       do 280 iback=ifirst,ilast
        i  = ilast - iback + ifirst
        j1 = j1 - lenr(i)
        if (lenrl(i) .eq. 0) go to 260
        j2 = j1 + lenrl(i) - 1
        do 250 jj=j1,j2
         j    = icn(jj)
         w(j) = w(j) + a(jj)*w(i)
250     continue
260     if (nobloc) go to 280
!
! operations using lower triangular blocks.
! lj2 points to the end of row i of the off-diagonal blocks.
! lj1 points to the beginning of row i of the off-diagonal blocks.
        if (lenoff(i) .eq. 0) go to 280
        lj2 = lj1 - 1
        lj1 = lj1 - lenoff(i)
        do 270 jj=lj1,lj2
         j = icn(jj)
         w(j) = w(j) - a(jj)*w(i)
270     continue
280    continue
       iblend = j1
       ilast = ifirst - 1
290   continue
!
! reorder solution vector ... x(i)=w(ipinverse(i))
300   do 310 ii=1,n
       i    = ip(ii)
       i    = iabs(i)
       x(i) = w(ii)
310   continue
320   return
      end
!
!
!
!
!
!
      subroutine ma30dd(a,icn,iptr,n,iactiv,itop,reals)
      include 'implno.dek'
!
! this subroutine performs garbage collection operations on the arrays a,
! icn and irn. iactiv is the first position in arrays a/icn from which the
! compress starts.  on exit, iactiv equals the position of the first entry
! in the compressed part of a/icn
!
      logical          reals
      integer          n,itop,iptr(n),icn(itop), &
                       irncp,icncp,irank,minirn,minicn,j,k,kn,kl, &
                       jpos,iactiv
      double precision a(itop)
      common /ma30fd/  irncp,icncp,irank,minirn,minicn
!
!
      if (reals) icncp = icncp + 1
      if (.not.reals) irncp = irncp + 1
!
! set the first non-zero entry in each row to the negative of the
! row/col number and hold this row/col index in the row/col
! pointer.  this is so that the beginning of each row/col can
! be recognized in the subsequent scan.
      do 10 j=1,n
       k = iptr(j)
       if (k .lt. iactiv) go to 10
       iptr(j) = icn(k)
       icn(k) = -j
10    continue
      kn = itop + 1
      kl = itop - iactiv + 1
!
! go through arrays in reverse order compressing to the back so
! that there are no zeros held in positions iactiv to itop in icn.
! reset first entry of each row/col and pointer array iptr.
      do 30 k=1,kl
       jpos = itop - k + 1
       if (icn(jpos) .eq. 0) go to 30
       kn = kn - 1
       if (reals) a(kn) = a(jpos)
       if (icn(jpos) .ge. 0) go to 20
!
! first non-zero of row/col has been located
       j         = -icn(jpos)
       icn(jpos) = iptr(j)
       iptr(j)   = kn
20     icn(kn)   = icn(jpos)
30    continue
      iactiv = kn
      return
      end
!
!
!
!
!
      subroutine ma28int1
      include 'implno.dek'
!
!
! lp,mp are used by the subroutine as the unit numbers for its warning
!       and diagnostic messages. default value for both is 6 (for line
!       printer output). the user can either reset them to a different
!       stream number or suppress the output by setting them to zero.
!       while lp directs the output of error diagnostics from the
!       principal subroutines and internally called subroutines, mp
!       controls only the output of a message which warns the user that he
!       has input two or more non-zeros a(i), . . ,a(k) with the same row
!       and column indices.  the action taken in this case is to proceed
!       using a numerical value of a(i)+...+a(k). in the absence of other
!       errors, iflag will equal -14 on exit.
! lblock is a logical variable which controls an option of first
!        preordering the matrix to block lower triangular form (using
!        harwell subroutine mc23a). the preordering is performed if lblock
!        is equal to its default value of .true. if lblock is set to
!        .false. , the option is not invoked and the space allocated to
!        ikeep can be reduced to 4*n+1.
! grow is a logical variable. if it is left at its default value of
!      .true. , then on return from ma28a/ad or ma28b/bd, w(1) will give
!      an estimate (an upper bound) of the increase in size of elements
!      encountered during the decomposition. if the matrix is well
!      scaled, then a high value for w(1), relative to the largest entry
!      in the input matrix, indicates that the lu decomposition may be
!      inaccurate and the user should be wary of his results and perhaps
!      increase u for subsequent runs.  we would like to emphasise that
!      this value only relates to the accuracy of our lu decomposition
!      and gives no indication as to the singularity of the matrix or the
!      accuracy of the solution.  this upper bound can be a significant
!      overestimate particularly if the matrix is badly scaled. if an
!      accurate value for the growth is required, lbig (q.v.) should be
!      set to .true.
! eps,rmin are real variables. if, on entry to ma28b/bd, eps is less
!      than one, then rmin will give the smallest ratio of the pivot to
!      the largest element in the corresponding row of the upper
!      triangular factor thus monitoring the stability of successive
!      factorizations. if rmin becomes very large and w(1) from
!      ma28b/bd is also very large, it may be advisable to perform a
!      new decomposition using ma28a/ad.
! resid is a real variable which on exit from ma28c/cd gives the value
!       of the maximum residual over all the equations unsatisfied because
!       of dependency (zero pivots).
! irncp,icncp are integer variables which monitor the adequacy of "elbow
!      room" in irn and a/icn respectively. if either is quite large (say
!      greater than n/10), it will probably pay to increase the size of
!      the corresponding array for subsequent runs. if either is very low
!      or zero then one can perhaps save storage by reducing the size of
!      the corresponding array.
! minirn,minicn are integer variables which, in the event of a
!      successful return (iflag ge 0 or iflag=-14) give the minimum size
!      of irn and a/icn respectively which would enable a successful run
!      on an identical matrix. on an exit with iflag equal to -5, minicn
!      gives the minimum value of icn for success on subsequent runs on
!      an identical matrix. in the event of failure with iflag= -6, -4,
!      -3, -2, or -1, then minicn and minirn give the minimum value of
!      licn and lirn respectively which would be required for a
!      successful decomposition up to the point at which the failure occurred.
! irank is an integer variable which gives an upper bound on the rank of
!       the matrix.
! abort1 is a logical variable with default value .true.  if abort1 is
!        set to .false.  then ma28a/ad will decompose structurally singular
!        matrices (including rectangular ones).
! abort2 is a logical variable with default value .true.  if abort2 is
!        set to .false. then ma28a/ad will decompose numerically singular
!        matrices.
! idisp is an integer array of length 2. on output from ma28a/ad, the
!       indices of the diagonal blocks of the factors lie in positions
!       idisp(1) to idisp(2) of a/icn. this array must be preserved
!       between a call to ma28a/ad and subsequent calls to ma28b/bd,
!       ma28c/cd or ma28i/id.
! tol is a real variable.  if it is set to a positive value, then any
!     non-zero whose modulus is less than tol will be dropped from the
!     factorization.  the factorization will then require less storage
!     but will be inaccurate.  after a run of ma28a/ad with tol positive
!     it is not possible to use ma28b/bd and the user is recommended to
!     use ma28i/id to obtain the solution.  the default value for tol is 0.0.
! themax is a real variable.  on exit from ma28a/ad, it will hold the
!        largest entry of the original matrix.
! big is a real variable. if lbig has been set to .true., big will hold
!     the largest entry encountered during the factorization by ma28a/ad
!      or ma28b/bd.
! dxmax is a real variable. on exit from ma28i/id, dxmax will be set to
!       the largest component of the solution.
! errmax is a real variable.  on exit from ma28i/id, if maxit is
!        positive, errmax will be set to the largest component in the
!        estimate of the error.
! dres is a real variable.  on exit from ma28i/id, if maxit is positive,
!      dres will be set to the largest component of the residual.
! cgce is a real variable. it is used by ma28i/id to check the
!      convergence rate.  if the ratio of successive corrections is
!      not less than cgce then we terminate since the convergence
!      rate is adjudged too slow.
! ndrop is an integer variable. if tol has been set positive, on exit
!      from ma28a/ad, ndrop will hold the number of entries dropped from
!      the data structure.
! maxit is an integer variable. it is the maximum number of iterations
!      performed by ma28i/id. it has a default value of 16.
! noiter is an integer variable. it is set by ma28i/id to the number of
!      iterative refinement iterations actually used.
! nsrch is an integer variable. if nsrch is set to a value less than n,
!      then a different pivot option will be employed by ma28a/ad.  this
!      may result in different fill-in and execution time for ma28a/ad.
!      if nsrch is less than or equal to n, the workspace array iw can be
!      reduced in length.  the default value for nsrch is 32768.
! istart is an integer variable. if istart is set to a value other than
!      zero, then the user must supply an estimate of the solution to
!      ma28i/id.  the default value for istart is zero.
! lbig is a logical variable. if lbig is set to .true., the value of the
!     largest element encountered in the factorization by ma28a/ad or
!     ma28b/bd is returned in big.  setting lbig to .true.  will
!     increase the time for ma28a/ad marginally and that for ma28b/bd
!     by about 20%.  the default value for lbig is .false.
!
! declare
      logical          lblock,grow,abort1,abort2,lbig
      integer          lp,mp,irncp,icncp,minirn,minicn,irank, &
                       ndrop,maxit,noiter,nsrch,istart
      double precision eps,rmin,resid,tol,themax,big,dxmax, &
                       errmax,dres,cgce
      common /ma28ed/  lp,mp,lblock,grow
      common /ma28fd/  eps,rmin,resid,irncp,icncp,minirn,minicn, &
                       irank,abort1,abort2
      common /ma28hd/  tol,themax,big,dxmax,errmax,dres,cgce, &
                       ndrop,maxit,noiter,nsrch,istart,lbig
!
      eps    = 1.0d-4
      tol    = 0.0d0
      cgce   = 0.5d0
      maxit  = 16
      lp     = 6
      mp     = 6
!      nsrch  = 1
      nsrch  = 32768
      istart = 0
!      lblock = .true.
      lblock = .false.
      grow   = .false.
      lbig   = .false.
      abort1 = .true.
      abort2 = .true.
      return
      end
!
!
!
!
!
      subroutine ma28int2
      include 'implno.dek'
!
!
! common block ma30e/ed holds control parameters
!     common /ma30ed/ lp, abort1, abort2, abort3
! the integer lp is the unit number to which the error messages are
! sent. lp has a default value of 6.  this default value can be
! reset by the user, if desired.  a value of 0 suppresses all
! messages.
! the logical variables abort1,abort2,abort3 are used to control the
! conditions under which the subroutine will terminate.
! if abort1 is .true. then the subroutine will exit  immediately on
! detecting structural singularity.
! if abort2 is .true. then the subroutine will exit immediately on
! detecting numerical singularity.
! if abort3 is .true. then the subroutine will exit immediately when
! the available space in a/icn is filled up by the previously decomposed,
! active, and undecomposed parts of the matrix.
!
! the default values for abort1,abort2,abort3 are set to .true.,.true.
! and .false. respectively.
!
!
! the variables in the common block ma30f/fd are used to provide the
! user with information on the decomposition.
! common /ma30fd/ irncp, icncp, irank, minirn, minicn
!
! irncp and icncp are integer variables used to monitor the adequacy
! of the allocated space in arrays irn and a/icn respectively, by
! taking account of the number of data management compresses
! required on these arrays. if irncp or icncp is fairly large (say
! greater than n/10), it may be advantageous to increase the size
! of the corresponding array(s).  irncp and icncp are initialized
! to zero on entry to ma30a/ad and are incremented each time the
! compressing routine ma30d/dd is entered.
!
! icncp is the number of compresses on a/icn.
! irncp is the number of compresses on irn.
!
! irank is an integer variable which gives an estimate (actually an
! upper bound) of the rank of the matrix. on an exit with iflag
! equal to 0, this will be equal to n.
!
! minirn is an integer variable which, after a successful call to
! ma30a/ad, indicates the minimum length to which irn can be
! reduced while still permitting a successful decomposition of the
! same matrix. if, however, the user were to decrease the length
! of irn to that size, the number of compresses (irncp) may be
! very high and quite costly. if lirn is not large enough to begin
! the decomposition on a diagonal block, minirn will be equal to
! the value required to continue the decomposition and iflag will
! be set to -3 or -6. a value of lirn slightly greater than this
! (say about n/2) will usually provide enough space to complete
! the decomposition on that block. in the event of any other
! failure minirn gives the minimum size of irn required for a
! successful decomposition up to that point.
!
! minicn is an integer variable which after a successful call to
! ma30a/ad, indicates the minimum size of licn required to enable
! a successful decomposition. in the event of failure with iflag=
! -5, minicn will, if abort3 is left set to .false., indicate the
! minimum length that would be sufficient to prevent this error in
! a subsequent run on an identical matrix. again the user may
! prefer to use a value of icn slightly greater than minicn for
! subsequent runs to avoid too many conpresses (icncp). in the
! event of failure with iflag equal to any negative value except
! -4, minicn will give the minimum length to which licn could be
! reduced to enable a successful decomposition to the point at
! which failure occurred.  notice that, on a successful entry
! idisp(2) gives the amount of space in a/icn required for the
! decomposition while minicn will usually be slightly greater
! because of the need for "elbow room".  if the user is very
! unsure how large to make licn, the variable minicn can be used
! to provide that information. a preliminary run should be
! performed with abort3 left set to .false. and licn about 3/2
! times as big as the number of non-zeros in the original matrix.
! unless the initial problem is very sparse (when the run will be
! successful) or fills in extremely badly (giving an error return
! with iflag equal to -4), an error return with iflag equal to -5
! should result and minicn will give the amount of space required
! for a successful decomposition.
!
!
! common block ma30g/gd is used by the ma30b/bd entry only.
!    common /ma30gd/ eps, rmin
! eps is a real/double precision variable. it is used to test for
! small pivots. its default value is 1.0e-4 (1.0d-4 in d version).
! if the user sets eps to any value greater than 1.0, then no
! check is made on the size of the pivots. although the absence of
! such a check would fail to warn the user of bad instability, its
! absence will enable ma30b/bd to run slightly faster. an  a
! posteriori  check on the stability of the factorization can be
! obtained from mc24a/ad.
!
! rmin is a real/double precision variable which gives the user some
! information about the stability of the decomposition.  at each
! stage of the lu decomposition the magnitude of the pivot apiv
! is compared with the largest off-diagonal entry currently in its
! row (row of u), rowmax say. if the ratio min (apiv/rowmax)
! where the minimum is taken over all the rows, is less than eps
! then rmin is set to this minimum value and iflag is returned
! with the value +i where i is the row in which this minimum
! occurs.  if the user sets eps greater than one, then this test
! is not performed. in this case, and when there are no small
! pivots rmin will be set equal to eps.
!
!
! common block ma30h/hd is used by ma30c/cd only.
!    common /ma30hd/ resid
! resid is a real/double precision variable. in the case of singular
! or rectangular matrices its final value will be equal to the
! maximum residual for the unsatisfied equations; otherwise its
! value will be set to zero.
!
!
! common  block ma30i/id controls the use of drop tolerances, the
! modified pivot option and the the calculation of the largest
! entry in the factorization process. this common block was added
! to the ma30 package in february, 1983.
!    common /ma30id/ tol, big, ndrop, nsrch, lbig
!
! tol is a real/double precision variable.  if it is set to a positive
! value, then ma30a/ad will drop from the factors any non-zero
! whose modulus is less than tol.  the factorization will then
! require less storage but will be inaccurate.  after a run of
! ma30a/ad where entries have been dropped, ma30b/bd  should not
! be called.  the default value for tol is 0.0.
!
! big is a real/double precision variable.  if lbig has been set to
! .true., big will be set to the largest entry encountered during
! the factorization.
! ndrop is an integer variable. if tol has been set positive, on exit
! from ma30a/ad, ndrop will hold the number of entries dropped
! from the data structure.
!
! nsrch is an integer variable. if nsrch is set to a value less than
! or equal to n, then a different pivot option will be employed by
! ma30a/ad.  this may result in different fill-in and execution
! time for ma30a/ad. if nsrch is less than or equal to n, the
! workspace arrays lastc and nextc are not referenced by ma30a/ad.
! the default value for nsrch is 32768.
! lbig is a logical variable. if lbig is set to .true., the value of
! the largest entry encountered in the factorization by ma30a/ad
! is returned in big.  setting lbig to .true.  will marginally
! increase the factorization time for ma30a/ad and will increase
! that for ma30b/bd by about 20%.  the default value for lbig is
! .false.
!
! declare
      logical          abort1,abort2,abort3,lbig
      integer          lp,ndrop,nsrch
      double precision eps,rmin,tol,big
!
      common /ma30ed/ lp,abort1,abort2,abort3
      common /ma30gd/ eps,rmin
      common /ma30id/ tol,big,ndrop,nsrch,lbig
!
      eps    = 1.0d-4
      tol    = 0.0d0
      big    = 0.0d0
      lp     = 6
      nsrch  = 32768
      lbig   = .false.
      abort1 = .true.
      abort2 = .true.
      abort3 = .false.
      return
      end
!
!
!
!
!
      subroutine ma28int3
      logical         abort
      integer         lp,numnz,num,large
      common /mc23bd/ lp,numnz,num,large,abort
      lp       = 6
      abort    = .false.
      return
      end
!
!
!
!
!
      subroutine mc20ad(nc,maxa,a,inum,jptr,jnum,jdisp)
      include 'implno.dek'
!
! sorts a matrix into row order
!
      integer          nc,maxa,inum(maxa),jnum(maxa),jptr(nc),jdisp, &
                       null,j,k,kr,ice,jce,ja,jb,i,loc,icep,jcep
      double precision a(maxa),ace,acep
!
! go
      null = -jdisp
      do 60 j=1,nc
       jptr(j) = 0
60    continue
!
! count the number of elements in each column.
      do 120 k=1,maxa
       j       = jnum(k) + jdisp
       jptr(j) = jptr(j) + 1
120   continue
!
! set the jptr array
      k = 1
      do 150 j=1,nc
       kr      = k + jptr(j)
       jptr(j) = k
       k       = kr
150   continue
!
! reorder the elements into column order; an in-place sort of order maxa.
!  jce is the current entry.
      do 230 i=1,maxa
       jce = jnum(i) + jdisp
       if (jce .eq. 0) go to 230
       ace = a(i)
       ice = inum(i)
!
! clear the location vacated.
       jnum(i) = null
!
! chain from current entry to store items.
       do 200 j=1,maxa
        loc        = jptr(jce)
        jptr(jce)  = jptr(jce) + 1
        acep       = a(loc)
        icep       = inum(loc)
        jcep       = jnum(loc)
        a(loc)     = ace
        inum(loc)  = ice
        jnum(loc)  = null
        if (jcep .eq. null) go to 230
        ace = acep
        ice = icep
        jce = jcep + jdisp
200    continue
230   continue
!
! reset jptr vector.
      ja = 1
      do 250 j=1,nc
       jb      = jptr(j)
       jptr(j) = ja
       ja      = jb
250   continue
      return
      end
!
!
!
!
!
      subroutine mc23ad(n,icn,a,licn,lenr,idisp,ip,iq,lenoff,iw,iw1)
      include 'implno.dek'
!
! performs the block triangularization
!
! declare
      logical          abort
      integer          n,licn,idisp(2),iw1(n,2),icn(licn),lenr(n), &
                       ip(n),iq(n),lenoff(n),iw(n,5),lp,numnz,num, &
                       large,i,ii,ibeg,iend,i1,i2,k,iblock,jnpos, &
                       ilend,inew,irowe,irowb,leni,nz,j,jj,iold,jold, &
                       jnew
      double precision a(licn)
      common /mc23bd/  lp,numnz,num,large,abort
!
! formats
180   format(1x,'matrix is structurally singular, rank = ',i6)
200   format(1x,'licn not big enough increase by ',i6)
220   format(1x,'error return from mc23ad because')
!
! set pointers iw(*,1) to beginning of the rows and set lenoff equal to lenr.
      iw1(1,1)  = 1
      lenoff(1) = lenr(1)
      if (n .eq. 1) go to 20
      do 10 i=2,n
       lenoff(i) = lenr(i)
       iw1(i,1)  = iw1(i-1,1) + lenr(i-1)
10    continue
!
! idisp(1) points to the first position in a/icn after the off-diagonal blocks
! and untreated rows.
20    idisp(1) = iw1(n,1) + lenr(n)
!
! find row permutation ip to make diagonal zero-free.
      call mc21a(n,icn,licn,iw1,lenr,ip,numnz,iw)
!
! possible error return for structurally singular matrices.
      if (numnz .ne. n  .and.  abort) go to 170
!
! iw1(*,2) and lenr are permutations of iw1(*,1) and lenr/lenoff suitable for
! entry to mc13d since matrix with these row pointer and length arrays has
! maximum number of non-zeros on the diagonal.
      do 30 ii=1,n
       i         = ip(ii)
       iw1(ii,2) = iw1(i,1)
       lenr(ii)  = lenoff(i)
30    continue
!
! find symmetric permutation iq to block lower triangular form.
      call mc13d(n,icn,licn,iw1(1,2),lenr,iq,iw(1,4),num,iw)
      if (num .ne. 1) go to 60
!
! action taken if matrix is irreducible. whole matrix is just moved to the
! end of the storage.
      do 40 i=1,n
       lenr(i) = lenoff(i)
       ip(i)   = i
       iq(i)   = i
40    continue
      lenoff(1) = -1
!
! idisp(1) is the first position after the last element in the off-diagonal
! blocks and untreated rows.
      nz       = idisp(1)-1
      idisp(1) = 1
!
! idisp(2) is position in a/icn of the first element in the diagonal blocks.
      idisp(2) = licn - nz + 1
      large    = n
      if (nz .eq. licn) go to 230
      do 50 k=1,nz
       j       = nz - k + 1
       jj      = licn - k + 1
       a(jj)   = a(j)
       icn(jj) = icn(j)
50    continue
      go to 230
!
! data structure reordered. form composite row permutation:ip(i) = ip(iq(i)).
60    do 70 ii=1,n
       i        = iq(ii)
       iw(ii,1) = ip(i)
70    continue
      do 80 i=1,n
       ip(i) = iw(i,1)
80    continue
!
! run through blocks in reverse order separating diagonal blocks which are
! moved to the end of the storage.  elements in off-diagonal blocks are left
! in place unless a compress is necessary.
! ibeg indicates the lowest value of j for which icn(j) has been
!      set to zero when element in position j was moved to the
!      diagonal block part of storage.
! iend is position of first element of those treated rows which are in
!      diagonal blocks.
! large is the dimension of the largest block encountered so far.
! num is the number of diagonal blocks.
! i1 is first row (in permuted form) of block iblock.
! i2 is last row (in permuted form) of block iblock.
      ibeg  = licn + 1
      iend  = licn + 1
      large = 0
      do 150 k=1,num
       iblock = num - k + 1
       i1 = iw(iblock,4)
       i2 = n
       if (k .ne. 1) i2 = iw(iblock+1,4) - 1
       large = max0(large,i2-i1+1)
!
! go through the rows of block iblock in the reverse order.
       do 140 ii=i1,i2
        inew = i2 - ii + i1
!
! we now deal with row inew in permuted form (row iold in original matrix).
        iold = ip(inew)
!
! if there is space to move up diagonal block portion of row go to 110
        if (iend-idisp(1) .ge. lenoff(iold)) go to 110
!
! in-line compress.; moves separated off-diagonal elements and untreated rows
! to front of storage.
        jnpos = ibeg
        ilend = idisp(1)-1
        if (ilend .lt. ibeg) go to 190
        do 90 j=ibeg,ilend
         if (icn(j) .eq. 0) go to 90
         icn(jnpos) = icn(j)
         a(jnpos)   = a(j)
         jnpos      = jnpos + 1
90      continue
        idisp(1) = jnpos
        if (iend-jnpos .lt. lenoff(iold)) go to 190
        ibeg = licn + 1
!
! reset pointers to the beginning of the rows.
        do 100 i=2,n
         iw1(i,1) = iw1(i-1,1) + lenoff(i-1)
100     continue
!
! row iold is now split into diag. and off-diag. parts.
110     irowb = iw1(iold,1)
        leni  = 0
        irowe = irowb+lenoff(iold)-1
!
! backward scan of whole of row iold (in original matrix).
        if (irowe .lt. irowb) go to 130
        do 120 jj=irowb,irowe
         j    = irowe - jj + irowb
         jold = icn(j)
!
! iw(.,2) holds the inverse permutation to iq.; it was set to this in mc13d.
         jnew = iw(jold,2)
!
! if (jnew.lt.i1) then element is in off-diagonal block and so is left in situ.
         if (jnew .lt. i1) go to 120
!
! element is in diagonal block and is moved to the end of the storage.
         iend      = iend-1
         a(iend)   = a(j)
         icn(iend) = jnew
         ibeg      = min0(ibeg,j)
         icn(j)    = 0
         leni      = leni + 1
120     continue
        lenoff(iold) = lenoff(iold) - leni
130     lenr(inew)   = leni
140    continue
       ip(i2)        = -ip(i2)
150   continue
!
! resets ip(n) to positive value.
! idisp(2) is position of first element in diagonal blocks.
      ip(n)    = -ip(n)
      idisp(2) = iend
!
! this compress used to move all off-diagonal elements to the front of storage.
      if (ibeg .gt. licn) go to 230
      jnpos = ibeg
      ilend = idisp(1) - 1
      do 160 j=ibeg,ilend
       if (icn(j) .eq. 0) go to 160
       icn(jnpos) = icn(j)
       a(jnpos)   = a(j)
       jnpos      = jnpos + 1
160   continue
!
! idisp(1) is first position after last element of off-diagonal blocks.
      idisp(1) = jnpos
      go to 230
!
! error return
170   if (lp .ne. 0) write(lp,180) numnz
      idisp(1) = -1
      go to 210
190   if (lp .ne. 0) write(lp,200) n
      idisp(1) = -2
210   if (lp .ne. 0) write(lp,220)
230   return
      end
!
!
!
!
!
      subroutine mc22ad(n,icn,a,nz,lenrow,ip,iq,iw,iw1)
      include 'implno.dek'
!
! reorders the off diagonal blocks based on the pivot information
!
! declare
      integer          n,nz,iw(n,2),icn(nz),lenrow(n),ip(n),iq(n), &
                       iw1(nz),i,jj,iold,j2,length,j,ipos,jval, &
                       ichain,newpos,jnum
      double precision a(nz),aval
!
! go
      if (nz .le. 0) go to 1000
      if (n  .le. 0) go to 1000
!
! set start of row i in iw(i,1) and lenrow(i) in iw(i,2)
      iw(1,1) = 1
      iw(1,2) = lenrow(1)
      do 10 i=2,n
       iw(i,1) = iw(i-1,1) + lenrow(i-1)
       iw(i,2) = lenrow(i)
10    continue
!
! permute lenrow according to ip.  set off-sets for new position of row iold
! in iw(iold,1) and put old row indices in iw1 in positions corresponding to
! the new position of this row in a/icn.
      jj = 1
      do 20 i=1,n
       iold      = ip(i)
       iold      = iabs(iold)
       length    = iw(iold,2)
       lenrow(i) = length
       if (length .eq. 0) go to 20
       iw(iold,1) = iw(iold,1) - jj
       j2 = jj + length - 1
       do 15 j=jj,j2
        iw1(j) = iold
15     continue
       jj = j2 + 1
20    continue
!
! set inverse permutation to iq in iw(.,2).
      do 30 i=1,n
       iold       = iq(i)
       iold       = iabs(iold)
       iw(iold,2) = i
30    continue
!
! permute a and icn in place, changing to new column numbers.
! main loop; each pass through this loop places a closed chain of column
! indices in their new (and final) positions ... this is recorded by
! setting the iw1 entry to zero so that any which are subsequently
! encountered during this major scan can be bypassed.
      do 200 i=1,nz
       iold = iw1(i)
       if (iold .eq. 0) go to 200
       ipos = i
       jval = icn(i)
!
! if row iold is in same positions after permutation go to 150.
       if (iw(iold,1) .eq. 0) go to 150
       aval = a(i)
!
! chain loop; each pass through this loop places one (permuted) column index
! in its final position  .. viz. ipos.
! newpos is the original position in a/icn of the element to be placed
! in position ipos.  it is also the position of the next element in the chain.
       do 100 ichain=1,nz
        newpos = ipos + iw(iold,1)
        if (newpos .eq. i) go to 130
        a(ipos)   = a(newpos)
        jnum      = icn(newpos)
        icn(ipos) = iw(jnum,2)
        ipos      = newpos
        iold      = iw1(ipos)
        iw1(ipos) = 0
100    continue
130    a(ipos)   = aval
150    icn(ipos) = iw(jval,2)
200   continue
1000  return
      end
!
!
!
!
!
      subroutine mc21a(n,icn,licn,ip,lenr,iperm,numnz,iw)
      include 'implno.dek'
      integer n,licn,ip(n),icn(licn),lenr(n),iperm(n),iw(n,4),numnz
      call mc21b(n,icn,licn,ip,lenr,iperm,numnz,iw(1,1),iw(1,2),iw(1,3), &
                 iw(1,4))
      return
      end
!
!
!
!
!
      subroutine mc21b(n,icn,licn,ip,lenr,iperm,numnz,pr,arp,cv,out)
      include 'implno.dek'
!
! does a row permutation to make the diagonal zero free
!
! pr(i) is the previous row to i in the depth first search.
!      it is used as a work array in the sorting algorithm.
!      elements (iperm(i),i) i=1, ... n  are non-zero at the end of the
!      algorithm unless n assignments have not been made.  in which case
! (iperm(i),i) will be zero for n-numnz entries.
! cv(i)  is the most recent row extension at which column i was visited.
! arp(i) is one less than the number of non-zeros in row i
!        which have not been scanned when looking for a cheap assignment.
! out(i) is one less than the number of non-zeros in row i
!        which have not been scanned during one pass through the main loop.
!
! declare
      integer n,licn,ip(n),icn(licn),lenr(n),iperm(n),pr(n),cv(n), &
              arp(n),out(n),i,jord,j,in1,in2,k,ii,ioutk,j1,kk,numnz
!
! initialization of arrays.
      do 10 i=1,n
       arp(i)   = lenr(i)-1
       cv(i)    = 0
       iperm(i) = 0
10    continue
      numnz=0
!
! main loop. each pass round this loop either results in a new assignment
! or gives a row with no assignment.
      do 130 jord=1,n
       j     = jord
       pr(j) = -1
       do 100 k=1,jord
!
! look for a cheap assignment
        in1 = arp(j)
        if (in1 .lt. 0) go to 60
        in2 = ip(j) + lenr(j)-1
        in1 = in2 - in1
        do 50 ii=in1,in2
         i = icn(ii)
         if (iperm(i) .eq. 0) go to 110
50      continue
!
! no cheap assignment in row.
! begin looking for assignment chain starting with row j.
        arp(j) = -1
60      out(j) = lenr(j)-1
!
! c inner loop.  extends chain by one or backtracks.
        do 90 kk=1,jord
         in1 = out(j)
         if (in1 .lt. 0) go to 80
         in2 = ip(j)+lenr(j)-1
         in1 = in2 - in1
!
! forward scan.
         do 70 ii=in1,in2
          i = icn(ii)
          if (cv(i) .eq. jord) go to 70
!
! column i has not yet been accessed during this pass.
          j1      = j
          j       = iperm(i)
          cv(i)   = jord
          pr(j)   = j1
          out(j1) = in2-ii-1
          go to 100
70       continue
!
! backtracking step.
80       j = pr(j)
         if (j .eq. -1) go to 130
90      continue
100    continue
!
! new assignment is made.
110    iperm(i) = j
       arp(j)   = in2 - ii - 1
       numnz    = numnz + 1
       do 120 k=1,jord
        j = pr(j)
        if (j .eq. -1) go to 130
        ii       = ip(j) + lenr(j) - out(j) - 2
        i        = icn(ii)
        iperm(i) = j
120    continue
130   continue
!
! if matrix is structurally singular, we now complete the permutation iperm.
      if (numnz .eq. n) return
      do 140 i=1,n
       arp(i) = 0
140   continue
      k = 0
      do 160 i=1,n
       if (iperm(i) .ne. 0) go to 150
       k      = k + 1
       out(k) = i
       go to 160
150    j      = iperm(i)
       arp(j) = i
160   continue
      k = 0
      do 170 i=1,n
       if (arp(i) .ne. 0) go to 170
       k            = k+1
       ioutk        = out(k)
       iperm(ioutk) = i
170   continue
      return
      end
!
!
!
!
!
      subroutine mc13d(n,icn,licn,ip,lenr,ior,ib,num,iw)
      include 'implno.dek'
      integer n,licn,ip(n),icn(licn),lenr(n),ior(n),ib(n),iw(n,3),num
      call mc13e(n,icn,licn,ip,lenr,ior,ib,num,iw(1,1),iw(1,2),iw(1,3))
      return
      end
!
!
!
!
!
      subroutine mc13e(n,icn,licn,ip,lenr,arp,ib,num,lowl,numb,prev)
      include 'implno.dek'
!
!  arp(i) is one less than the number of unsearched edges leaving
!         node i.  at the end of the algorithm it is set to a
!         permutation which puts the matrix in block lower
!         triangular form.
! ib(i)   is the position in the ordering of the start of the ith
!         block.  ib(n+1-i) holds the node number of the ith node
!         on the stack.
! lowl(i) is the smallest stack position of any node to which a path
!         from node i has been found.  it is set to n+1 when node i
!         is removed from the stack.
! numb(i) is the position of node i in the stack if it is on
!         it, is the permuted order of node i for those nodes
!         whose final position has been found and is otherwise zero.
! prev(i) is the node at the end of the path when node i was
!         placed on the stack.
!
! declare
      integer n,licn,stp,dummy,ip(n),icn(licn),lenr(n),arp(n),ib(n), &
              lowl(n),numb(n),prev(n),icnt,num,nnm1,j,iv,ist,i1,i2, &
              ii,iw,ist1,lcnt,i,isn,k
!
!
! icnt is number of nodes whose positions in final ordering have been found.
! num is the number of blocks that have been found.
      icnt = 0
      num  = 0
      nnm1 = n + n-1
!
! initialization of arrays.
      do 20 j=1,n
       numb(j) = 0
       arp(j)  = lenr(j)-1
20    continue
!
! look for a starting node
! ist is the number of nodes on the stack ... it is the stack pointer.
      do 120 isn=1,n
       if (numb(isn) .ne. 0) go to 120
       iv  = isn
       ist = 1
!
! put node iv at beginning of stack.
       lowl(iv) = 1
       numb(iv) = 1
       ib(n)    = iv
!
! the body of this loop puts a new node on the stack or backtracks.
       do 110 dummy=1,nnm1
        i1 = arp(iv)
!
! have all edges leaving node iv been searched.
        if (i1 .lt. 0) go to 60
        i2 = ip(iv) + lenr(iv) - 1
        i1 = i2 - i1
!
! look at edges leaving node iv until one enters a new node or all edges are
! exhausted.
        do 50 ii=i1,i2
         iw = icn(ii)
         if (numb(iw) .eq. 0) go to 100
         lowl(iv) = min0(lowl(iv),lowl(iw))
50      continue
!
! there are no more edges leaving node iv.
        arp(iv) = -1
!
! is node iv the root of a block.
60      if (lowl(iv) .lt. numb(iv)) go to 90
!
! order nodes in a block.
        num  = num + 1
        ist1 = n + 1 - ist
        lcnt = icnt + 1
!
! peel block off the top of the stack starting at the top and working down to
! the root of the block.
        do 70 stp=ist1,n
         iw       = ib(stp)
         lowl(iw) = n + 1
         icnt     = icnt + 1
         numb(iw) = icnt
         if (iw .eq. iv) go to 80
70      continue
80      ist     = n - stp
        ib(num) = lcnt
!
! are there any nodes left on the stack.
        if (ist .ne. 0) go to 90
!
! have all the nodes been ordered.
        if (icnt .lt. n) go to 120
        go to 130
!
! backtrack to previous node on path.
90      iw = iv
        iv = prev(iv)
!
! update value of lowl(iv) if necessary.
        lowl(iv) = min0(lowl(iv),lowl(iw))
        go to 110
!
! put new node on the stack.
100     arp(iv)  = i2 - ii - 1
        prev(iw) = iv
        iv       = iw
        ist      = ist+1
        lowl(iv) = ist
        numb(iv) = ist
        k        = n+1-ist
        ib(k)    = iv
110    continue
120   continue
!
! put permutation in the required form.
130   do 140 i=1,n
       ii      = numb(i)
       arp(ii) = i
140   continue
      return
      end
!
!
!
!
!
      subroutine mc24ad(n,icn,a,licn,lenr,lenrl,w)
      include 'implno.dek'
!
! computes the gwoth rate of fill in
!
      integer          n,licn,icn(licn),lenr(n),lenrl(n),i,j0,j2,j1,jj,j
      double precision a(licn),w(n),amaxl,wrowl,amaxu,zero
      data             zero/0.0d0/
!
! initialize
      amaxl = zero
      do 10 i=1,n
       w(i) = zero
10    continue
      j0 = 1
      do 100 i=1,n
       if (lenr(i) .eq. 0) go to 100
       j2=j0+lenr(i)-1
       if (lenrl(i) .eq. 0) go to 50
!
! calculation of 1-norm of l.
       j1 = j0 + lenrl(i) - 1
       wrowl=zero
       do 30 jj=j0,j1
        wrowl = wrowl + abs(a(jj))
30     continue
!
! amaxl is the maximum norm of columns of l so far found.
       amaxl = max(amaxl,wrowl)
       j0    = j1 + 1
!
! calculation of norms of columns of u (max-norms).
50     j0 = j0 + 1
       if (j0 .gt. j2) go to 90
       do 80 jj=j0,j2
        j    = icn(jj)
        w(j) = max(abs(a(jj)),w(j))
80     continue
90     j0 = j2 + 1
100   continue
!
! amaxu is set to maximum max-norm of columns of u.
      amaxu = zero
      do 200 i=1,n
       amaxu = max(amaxu,w(i))
200   continue
!
! grofac is max u max-norm times max l 1-norm.
      w(1) = amaxl*amaxu
      return
      end
!
!
!
!
!
      subroutine mc20bd(nc,maxa,a,inum,jptr)
      include 'implno.dek'
!
! never called
!
!
! declare
      integer          nc,maxa,inum(maxa),jptr(nc),kmax,jj,j,klo,kor, &
                       kdummy,ice,k,ik
      double precision a(maxa),ace
!
! go
      kmax=maxa
      do 35 jj=1,nc
       j   = nc + 1 - jj
       klo = jptr(j)+1
       if (klo .gt. kmax) go to 30
       kor=kmax
!
! items kor, kor+1, .... ,kmax are in order
       do 25 kdummy=klo,kmax
        ace = a(kor-1)
        ice = inum(kor-1)
        do 10 k=kor,kmax
         ik = inum(k)
         if (iabs(ice) .le. iabs(ik)) go to 20
         inum(k-1) = ik
         a(k-1)    = a(k)
10      continue
        k         = kmax+1
20      inum(k-1) = ice
        a(k-1)    = ace
        kor       = kor-1
25     continue
30     kmax = klo - 2
35    continue
      return
      end
!---------------------------------------------------------------------







!---------------------------------------------------------------------
      subroutine read_helm_table
      include 'implno.dek'
      include 'helm_table_storage.dek'

! this routine reads the helmholtz eos file, and
! must be called once before the helmeos routine is invoked.

! declare local variables
      integer          i,j
      double precision tsav,dsav,dth,dt2,dti,dt2i,dt3i, &
                       dd,dd2,ddi,dd2i,dd3i


! open the file (use softlinks to input the desired table)

!       open(unit=19,file='helm_table.dat',status='old')
!       open(unit=19,file='aprox13/helm_table.dat',status='old')
       open(unit=19,file='data.code/helm_table.dat',status='old')


! for standard table limits
       tlo   = 3.0d0
       thi   = 13.0d0
       tstp  = (thi - tlo)/float(jmax-1)
       tstpi = 1.0d0/tstp
       dlo   = -12.0d0
       dhi   = 15.0d0
       dstp  = (dhi - dlo)/float(imax-1)
       dstpi = 1.0d0/dstp

! read the helmholtz free energy and its derivatives
       do j=1,jmax
        tsav = tlo + (j-1)*tstp
        t(j) = 10.0d0**(tsav)
        do i=1,imax
         dsav = dlo + (i-1)*dstp
         d(i) = 10.0d0**(dsav)
         read(19,*) f(i,j),fd(i,j),ft(i,j),fdd(i,j),ftt(i,j),fdt(i,j), &
                  fddt(i,j),fdtt(i,j),fddtt(i,j)
        enddo
       enddo


! read the pressure derivative with density table
       do j=1,jmax
        do i=1,imax
         read(19,*) dpdf(i,j),dpdfd(i,j),dpdft(i,j),dpdfdt(i,j)
        enddo
       enddo

! read the electron chemical potential table
       do j=1,jmax
        do i=1,imax
         read(19,*) ef(i,j),efd(i,j),eft(i,j),efdt(i,j)
        enddo
       enddo

! read the number density table
       do j=1,jmax
        do i=1,imax
         read(19,*) xf(i,j),xfd(i,j),xft(i,j),xfdt(i,j)
        enddo
       enddo

! close the file and write a summary message
      close(unit=19)


! construct the temperature and density deltas and their inverses
       do j=1,jmax-1
        dth          = t(j+1) - t(j)
        dt2         = dth * dth
        dti         = 1.0d0/dth
        dt2i        = 1.0d0/dt2
        dt3i        = dt2i*dti
        dt_sav(j)   = dth
        dt2_sav(j)  = dt2
        dti_sav(j)  = dti
        dt2i_sav(j) = dt2i
        dt3i_sav(j) = dt3i
       end do
       do i=1,imax-1
        dd          = d(i+1) - d(i)
        dd2         = dd * dd
        ddi         = 1.0d0/dd
        dd2i        = 1.0d0/dd2
        dd3i        = dd2i*ddi
        dd_sav(i)   = dd
        dd2_sav(i)  = dd2
        ddi_sav(i)  = ddi
        dd2i_sav(i) = dd2i
        dd3i_sav(i) = dd3i
       enddo



!      write(6,*)
!      write(6,*) 'finished reading eos table'
!      write(6,04) 'imax=',imax,' jmax=',jmax
!04    format(1x,4(a,i4))
!      write(6,03) 'temp(1)   =',t(1),' temp(jmax)   =',t(jmax)
!      write(6,03) 'ye*den(1) =',d(1),' ye*den(imax) =',d(imax)
!03    format(1x,4(a,1pe11.3))
!      write(6,*)

      return
      end







!      subroutine helmeos
      subroutine helmeosaprox13
      include 'implno.dek'
      include 'const.dek'
      include 'vector_eos.dek'
      include 'helm_table_storage.dek'


! given a temperature temp [K], density den [g/cm**3], and a composition
! characterized by abar and zbar, this routine returns most of the other
! thermodynamic quantities. of prime interest is the pressure [erg/cm**3],
! specific thermal energy [erg/gr], the entropy [erg/g/K], along with
! their derivatives with respect to temperature, density, abar, and zbar.
! other quantites such the normalized chemical potential eta (plus its
! derivatives), number density of electrons and positron pair (along
! with their derivatives), adiabatic indices, specific heats, and
! relativistically correct sound speed are also returned.
!
! this routine assumes planckian photons, an ideal gas of ions,
! and an electron-positron gas with an arbitrary degree of relativity
! and degeneracy. interpolation in a table of the helmholtz free energy
! is used to return the electron-positron thermodynamic quantities.
! all other derivatives are analytic.
!
! references: cox & giuli chapter 24 ; timmes & swesty apj 1999


! declare
      integer          i,j
      double precision temp,den,abar,zbar,ytot1,ye, &
                       x,y,zz,zzi,deni,tempi,xni,dxnidd,dxnida, &
                       dpepdt,dpepdd,deepdt,deepdd,dsepdd,dsepdt, &
                       dpraddd,dpraddt,deraddd,deraddt,dpiondd,dpiondt, &
                       deiondd,deiondt,dsraddd,dsraddt,dsiondd,dsiondt, &
                       dse,dpe,dsp,kt,ktinv,prad,erad,srad,pion,eion, &
                       sion,xnem,pele,eele,sele,pres,ener,entr,dpresdd, &
                       dpresdt,denerdd,denerdt,dentrdd,dentrdt,cv,cp, &
                       gam1,gam2,gam3,chit,chid,nabad,sound,etaele, &
                       detadt,detadd,xnefer,dxnedt,dxnedd,s

      double precision pgas,dpgasdd,dpgasdt,dpgasda,dpgasdz, &
                       egas,degasdd,degasdt,degasda,degasdz, &
                       sgas,dsgasdd,dsgasdt,dsgasda,dsgasdz, &
                       cv_gas,cp_gas,gam1_gas,gam2_gas,gam3_gas, &
                       chit_gas,chid_gas,nabad_gas,sound_gas


      double precision sioncon,forth,forpi,kergavo,ikavo,asoli3,light2
      parameter        (sioncon = (2.0d0 * pi * amu * kerg)/(h*h), &
                        forth   = 4.0d0/3.0d0, &
                        forpi   = 4.0d0 * pi, &
                        kergavo = kerg * avo, &
                        ikavo   = 1.0d0/kergavo, &
                        asoli3  = asol/3.0d0, &
                        light2  = clight * clight)

! for the abar derivatives
      double precision dpradda,deradda,dsradda, &
                       dpionda,deionda,dsionda, &
                       dpepda,deepda,dsepda, &
                       dpresda,denerda,dentrda, &
                       detada,dxneda

! for the zbar derivatives
      double precision dpraddz,deraddz,dsraddz, &
                       dpiondz,deiondz,dsiondz, &
                       dpepdz,deepdz,dsepdz, &
                       dpresdz,denerdz,dentrdz, &
                       detadz,dxnedz

! for the interpolations
      integer          iat,jat
      double precision free,df_d,df_t,df_dd,df_tt,df_dt
      double precision xt,xd,mxt,mxd, &
                       si0t,si1t,si2t,si0mt,si1mt,si2mt, &
                       si0d,si1d,si2d,si0md,si1md,si2md, &
                       dsi0t,dsi1t,dsi2t,dsi0mt,dsi1mt,dsi2mt, &
                       dsi0d,dsi1d,dsi2d,dsi0md,dsi1md,dsi2md, &
                       ddsi0t,ddsi1t,ddsi2t,ddsi0mt,ddsi1mt,ddsi2mt, &
                       ddsi0d,ddsi1d,ddsi2d,ddsi0md,ddsi1md,ddsi2md, &
                       z,psi0,dpsi0,ddpsi0,psi1,dpsi1,ddpsi1,psi2, &
                       dpsi2,ddpsi2,din,h5,fi(36), &
                       xpsi0,xdpsi0,xpsi1,xdpsi1,h3, &
                       w0t,w1t,w2t,w0mt,w1mt,w2mt, &
                       w0d,w1d,w2d,w0md,w1md,w2md


! for the uniform background coulomb correction
      double precision dsdd,dsda,lami,inv_lami,lamida,lamidd, &
                       plasg,plasgdd,plasgdt,plasgda,plasgdz, &
                       ecoul,decouldd,decouldt,decoulda,decouldz, &
                       pcoul,dpcouldd,dpcouldt,dpcoulda,dpcouldz, &
                       scoul,dscouldd,dscouldt,dscoulda,dscouldz, &
                       a1,b1,c1,d1,e1,a2,b2,c2,third,esqu
      parameter        (a1    = -0.898004d0, &
                        b1    =  0.96786d0, &
                        c1    =  0.220703d0, &
                        d1    = -0.86097d0, &
                        e1    =  2.5269d0, &
                        a2    =  0.29561d0, &
                        b2    =  1.9885d0, &
                        c2    =  0.288675d0, &
                        third =  1.0d0/3.0d0, &
                        esqu  =  qe * qe)


! quintic hermite polynomial statement functions
! psi0 and its derivatives
      psi0(z)   = z**3 * ( z * (-6.0d0*z + 15.0d0) -10.0d0) + 1.0d0
      dpsi0(z)  = z**2 * ( z * (-30.0d0*z + 60.0d0) - 30.0d0)
      ddpsi0(z) = z* ( z*( -120.0d0*z + 180.0d0) -60.0d0)


! psi1 and its derivatives
      psi1(z)   = z* ( z**2 * ( z * (-3.0d0*z + 8.0d0) - 6.0d0) + 1.0d0)
      dpsi1(z)  = z*z * ( z * (-15.0d0*z + 32.0d0) - 18.0d0) +1.0d0
      ddpsi1(z) = z * (z * (-60.0d0*z + 96.0d0) -36.0d0)


! psi2  and its derivatives
      psi2(z)   = 0.5d0*z*z*( z* ( z * (-z + 3.0d0) - 3.0d0) + 1.0d0)
      dpsi2(z)  = 0.5d0*z*( z*(z*(-5.0d0*z + 12.0d0) - 9.0d0) + 2.0d0)
      ddpsi2(z) = 0.5d0*(z*( z * (-20.0d0*z + 36.0d0) - 18.0d0) + 2.0d0)


! biquintic hermite polynomial statement function
      h5(i,j,w0t,w1t,w2t,w0mt,w1mt,w2mt,w0d,w1d,w2d,w0md,w1md,w2md)= &
             fi(1)  *w0d*w0t   + fi(2)  *w0md*w0t &
           + fi(3)  *w0d*w0mt  + fi(4)  *w0md*w0mt &
           + fi(5)  *w0d*w1t   + fi(6)  *w0md*w1t &
           + fi(7)  *w0d*w1mt  + fi(8)  *w0md*w1mt &
           + fi(9)  *w0d*w2t   + fi(10) *w0md*w2t &
           + fi(11) *w0d*w2mt  + fi(12) *w0md*w2mt &
           + fi(13) *w1d*w0t   + fi(14) *w1md*w0t &
           + fi(15) *w1d*w0mt  + fi(16) *w1md*w0mt &
           + fi(17) *w2d*w0t   + fi(18) *w2md*w0t &
           + fi(19) *w2d*w0mt  + fi(20) *w2md*w0mt &
           + fi(21) *w1d*w1t   + fi(22) *w1md*w1t &
           + fi(23) *w1d*w1mt  + fi(24) *w1md*w1mt &
           + fi(25) *w2d*w1t   + fi(26) *w2md*w1t &
           + fi(27) *w2d*w1mt  + fi(28) *w2md*w1mt &
           + fi(29) *w1d*w2t   + fi(30) *w1md*w2t &
           + fi(31) *w1d*w2mt  + fi(32) *w1md*w2mt &
           + fi(33) *w2d*w2t   + fi(34) *w2md*w2t &
           + fi(35) *w2d*w2mt  + fi(36) *w2md*w2mt



! cubic hermite polynomial statement functions
! psi0 & derivatives
      xpsi0(z)  = z * z * (2.0d0*z - 3.0d0) + 1.0
      xdpsi0(z) = z * (6.0d0*z - 6.0d0)


! psi1 & derivatives
      xpsi1(z)  = z * ( z * (z - 2.0d0) + 1.0d0)
      xdpsi1(z) = z * (3.0d0*z - 4.0d0) + 1.0d0


! bicubic hermite polynomial statement function
      h3(i,j,w0t,w1t,w0mt,w1mt,w0d,w1d,w0md,w1md) = &
             fi(1)  *w0d*w0t   +  fi(2)  *w0md*w0t &
           + fi(3)  *w0d*w0mt  +  fi(4)  *w0md*w0mt &
           + fi(5)  *w0d*w1t   +  fi(6)  *w0md*w1t &
           + fi(7)  *w0d*w1mt  +  fi(8)  *w0md*w1mt &
           + fi(9)  *w1d*w0t   +  fi(10) *w1md*w0t &
           + fi(11) *w1d*w0mt  +  fi(12) *w1md*w0mt &
           + fi(13) *w1d*w1t   +  fi(14) *w1md*w1t &
           + fi(15) *w1d*w1mt  +  fi(16) *w1md*w1mt



! popular format statements
01    format(1x,5(a,1pe11.3))
02    format(1x,a,1p4e16.8)
03    format(1x,4(a,1pe11.3))
04    format(1x,4(a,i4))



! start of pipeline loop, normal execution starts here
      eosfail = .false.
      do j=jlo_eos,jhi_eos

!       if (temp_row(j) .le. 0.0) stop 'temp less than 0 in helmeos'
!       if (den_row(j)  .le. 0.0) stop 'den less than 0 in helmeos'

       temp  = temp_row(j)
       den   = den_row(j)
       abar  = abar_row(j)
       zbar  = zbar_row(j)
       ytot1 = 1.0d0/abar
       ye    = max(1.0d-16,ytot1 * zbar)



! initialize
       deni    = 1.0d0/den
       tempi   = 1.0d0/temp
       kt      = kerg * temp
       ktinv   = 1.0d0/kt


! radiation section:
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


! ion section:
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


! sackur-tetrode equation for the ion entropy of
! a single ideal gas characterized by abar
        x       = abar*abar*sqrt(abar) * deni/avo
        s       = sioncon * temp
        z       = x * s * sqrt(s)
        y       = log(z)

!        y       = 1.0d0/(abar*kt)
!        yy      = y * sqrt(y)
!        z       = xni * sifac * yy
!        etaion  = log(z)


        sion    = (pion*deni + eion)*tempi + kergavo * ytot1 * y
        dsiondd = (dpiondd*deni - pion*deni*deni + deiondd)*tempi &
                   - kergavo * deni * ytot1
        dsiondt = (dpiondt*deni + deiondt)*tempi - &
                  (pion*deni + eion) * tempi*tempi &
                  + 1.5d0 * kergavo * tempi*ytot1
        x       = avo*kerg/abar
        dsionda = (dpionda*deni + deionda)*tempi &
                  + kergavo*ytot1*ytot1* (2.5d0 - y)
        dsiondz = 0.0d0



! electron-positron section:


! assume complete ionization
        xnem    = xni * zbar


! enter the table with ye*den
        din = ye*den


! bomb proof the input
        if (temp .gt. t(jmax)) then
         write(6,01) 'temp=',temp,' t(jmax)=',t(jmax)
         write(6,*) 'temp too hot, off grid'
         write(6,*) 'setting eosfail to true and returning'
         eosfail = .true.
         return
        end if
        if (temp .lt. t(1)) then
         write(6,01) 'temp=',temp,' t(1)=',t(1)
         write(6,*) 'temp too cold, off grid'
         write(6,*) 'setting eosfail to true and returning'
         eosfail = .true.
         return
        end if
        if (din  .gt. d(imax)) then
         write(6,01) 'den*ye=',din,' d(imax)=',d(imax)
         write(6,*) 'ye*den too big, off grid'
         write(6,*) 'setting eosfail to true and returning'
         eosfail = .true.
         return
        end if
        if (din  .lt. d(1)) then
         write(6,01) 'ye*den=',din,' d(1)=',d(1)
         write(6,*) 'ye*den too small, off grid'
         write(6,*) 'setting eosfail to true and returning'
         eosfail = .true.
         return
        end if

! hash locate this temperature and density
        jat = int((log10(temp) - tlo)*tstpi) + 1
        jat = max(1,min(jat,jmax-1))
        iat = int((log10(din) - dlo)*dstpi) + 1
        iat = max(1,min(iat,imax-1))


! access the table locations only once
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


! various differences
        xt  = max( (temp - t(jat))*dti_sav(jat), 0.0d0)
        xd  = max( (din - d(iat))*ddi_sav(iat), 0.0d0)
        mxt = 1.0d0 - xt
        mxd = 1.0d0 - xd

! the six density and six temperature basis functions
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

! derivatives of the weight functions
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

! second derivatives of the weight functions
        ddsi0t =   ddpsi0(xt)*dt2i_sav(jat)
        ddsi1t =   ddpsi1(xt)*dti_sav(jat)
        ddsi2t =   ddpsi2(xt)

        ddsi0mt =  ddpsi0(mxt)*dt2i_sav(jat)
        ddsi1mt = -ddpsi1(mxt)*dti_sav(jat)
        ddsi2mt =  ddpsi2(mxt)

!        ddsi0d =   ddpsi0(xd)*dd2i_sav(iat)
!        ddsi1d =   ddpsi1(xd)*ddi_sav(iat)
!        ddsi2d =   ddpsi2(xd)

!        ddsi0md =  ddpsi0(mxd)*dd2i_sav(iat)
!        ddsi1md = -ddpsi1(mxd)*ddi_sav(iat)
!        ddsi2md =  ddpsi2(mxd)


! the free energy
        free  = h5(iat,jat, &
                si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt, &
                si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

! derivative with respect to density
        df_d  = h5(iat,jat, &
                si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt, &
                dsi0d,  dsi1d,  dsi2d,  dsi0md,  dsi1md,  dsi2md)


! derivative with respect to temperature
        df_t = h5(iat,jat, &
                dsi0t,  dsi1t,  dsi2t,  dsi0mt,  dsi1mt,  dsi2mt, &
                si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

! derivative with respect to density**2
!        df_dd = h5(iat,jat,
!     1          si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt,
!     2          ddsi0d, ddsi1d, ddsi2d, ddsi0md, ddsi1md, ddsi2md)

! derivative with respect to temperature**2
        df_tt = h5(iat,jat, &
              ddsi0t, ddsi1t, ddsi2t, ddsi0mt, ddsi1mt, ddsi2mt, &
                si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

! derivative with respect to temperature and density
        df_dt = h5(iat,jat, &
                dsi0t,  dsi1t,  dsi2t,  dsi0mt,  dsi1mt,  dsi2mt, &
                dsi0d,  dsi1d,  dsi2d,  dsi0md,  dsi1md,  dsi2md)



! now get the pressure derivative with density, chemical potential, and
! electron positron number densities
! get the interpolation weight functions
        si0t   =  xpsi0(xt)
        si1t   =  xpsi1(xt)*dt_sav(jat)

        si0mt  =  xpsi0(mxt)
        si1mt  =  -xpsi1(mxt)*dt_sav(jat)

        si0d   =  xpsi0(xd)
        si1d   =  xpsi1(xd)*dd_sav(iat)

        si0md  =  xpsi0(mxd)
        si1md  =  -xpsi1(mxd)*dd_sav(iat)


! derivatives of weight functions
        dsi0t  = xdpsi0(xt)*dti_sav(jat)
        dsi1t  = xdpsi1(xt)

        dsi0mt = -xdpsi0(mxt)*dti_sav(jat)
        dsi1mt = xdpsi1(mxt)

        dsi0d  = xdpsi0(xd)*ddi_sav(iat)
        dsi1d  = xdpsi1(xd)

        dsi0md = -xdpsi0(mxd)*ddi_sav(iat)
        dsi1md = xdpsi1(mxd)


! look in the pressure derivative only once
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

! pressure derivative with density
        dpepdd  = h3(iat,jat, &
                       si0t,   si1t,   si0mt,   si1mt, &
                       si0d,   si1d,   si0md,   si1md)
        dpepdd  = max(ye * dpepdd,1.0d-30)



! look in the electron chemical potential table only once
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


! electron chemical potential etaele
        etaele  = h3(iat,jat, &
                     si0t,   si1t,   si0mt,   si1mt, &
                     si0d,   si1d,   si0md,   si1md)


! derivative with respect to density
        x       = h3(iat,jat, &
                     si0t,   si1t,   si0mt,   si1mt, &
                    dsi0d,  dsi1d,  dsi0md,  dsi1md)
        detadd  = ye * x

! derivative with respect to temperature
        detadt  = h3(iat,jat, &
                    dsi0t,  dsi1t,  dsi0mt,  dsi1mt, &
                     si0d,   si1d,   si0md,   si1md)

! derivative with respect to abar and zbar
       detada = -x * din * ytot1
       detadz =  x * den * ytot1



! look in the number density table only once
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

! electron + positron number densities
       xnefer   = h3(iat,jat, &
                     si0t,   si1t,   si0mt,   si1mt, &
                     si0d,   si1d,   si0md,   si1md)

! derivative with respect to density
       x        = h3(iat,jat, &
                     si0t,   si1t,   si0mt,   si1mt, &
                    dsi0d,  dsi1d,  dsi0md,  dsi1md)
       x = max(x,1.0d-30)
       dxnedd   = ye * x

! derivative with respect to temperature
       dxnedt   = h3(iat,jat, &
                    dsi0t,  dsi1t,  dsi0mt,  dsi1mt, &
                     si0d,   si1d,   si0md,   si1md)

! derivative with respect to abar and zbar
       dxneda = -x * din * ytot1
       dxnedz =  x  * den * ytot1


! the desired electron-positron thermodynamic quantities

! dpepdd at high temperatures and low densities is below the
! floating point limit of the subtraction of two large terms.
! since dpresdd doesn't enter the maxwell relations at all, use the
! bicubic interpolation done above instead of the formally correct expression
        x       = din * din
        pele    = x * df_d
        dpepdt  = x * df_dt
!        dpepdd  = ye * (x * df_dd + 2.0d0 * din * df_d)
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




! coulomb section:

! uniform background corrections only
! from yakovlev & shalybkov 1989
! lami is the average ion seperation
! plasg is the plasma coupling parameter

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


! yakovlev & shalybkov 1989 equations 82, 85, 86, 87
        if (plasg .ge. 1.0) then
         x        = plasg**(0.25d0)
         y        = avo * ytot1 * kerg
         ecoul    = y * temp * (a1*plasg + b1*x + c1/x + d1)
         pcoul    = third * den * ecoul
         scoul    = -y * (3.0d0*b1*x - 5.0d0*c1/x &
                    + d1 * (log(plasg) - 1.0d0) - e1)

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


! yakovlev & shalybkov 1989 equations 102, 103, 104
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


! bomb proof
        x   = prad + pion + pele + pcoul
        y   = erad + eion + eele + ecoul
        z   = srad + sion + sele + scoul

!        write(6,*) x,y,z
        if (x .le. 0.0 .or. y .le. 0.0 .or. z .le. 0.0) then
!        if (x .le. 0.0 .or. y .le. 0.0) then
!        if (x .le. 0.0) then

!         write(6,*)
!         write(6,*) 'coulomb corrections are causing a negative pressure'
!         write(6,*) 'setting all coulomb corrections to zero'
!         write(6,*)

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


! sum all the gas components
       pgas    = pion + pele + pcoul
       egas    = eion + eele + ecoul
       sgas    = sion + sele + scoul

       dpgasdd = dpiondd + dpepdd + dpcouldd
       dpgasdt = dpiondt + dpepdt + dpcouldt
       dpgasda = dpionda + dpepda + dpcoulda
       dpgasdz = dpiondz + dpepdz + dpcouldz

       degasdd = deiondd + deepdd + decouldd
       degasdt = deiondt + deepdt + decouldt
       degasda = deionda + deepda + decoulda
       degasdz = deiondz + deepdz + decouldz

       dsgasdd = dsiondd + dsepdd + dscouldd
       dsgasdt = dsiondt + dsepdt + dscouldt
       dsgasda = dsionda + dsepda + dscoulda
       dsgasdz = dsiondz + dsepdz + dscouldz




! add in radiation to get the total
       pres    = prad + pgas
       ener    = erad + egas
       entr    = srad + sgas

       dpresdd = dpraddd + dpgasdd
       dpresdt = dpraddt + dpgasdt
       dpresda = dpradda + dpgasda
       dpresdz = dpraddz + dpgasdz

       denerdd = deraddd + degasdd
       denerdt = deraddt + degasdt
       denerda = deradda + degasda
       denerdz = deraddz + degasdz

       dentrdd = dsraddd + dsgasdd
       dentrdt = dsraddt + dsgasdt
       dentrda = dsradda + dsgasda
       dentrdz = dsraddz + dsgasdz


! for the gas
! the temperature and density exponents (c&g 9.81 9.82)
! the specific heat at constant volume (c&g 9.92)
! the third adiabatic exponent (c&g 9.93)
! the first adiabatic exponent (c&g 9.97)
! the second adiabatic exponent (c&g 9.105)
! the specific heat at constant pressure (c&g 9.98)
! and relativistic formula for the sound speed (c&g 14.29)

       zz        = pgas*deni
       zzi       = den/pgas
       chit_gas  = temp/pgas * dpgasdt
       chid_gas  = dpgasdd*zzi
       cv_gas    = degasdt
       x         = zz * chit_gas/(temp * cv_gas)
       gam3_gas  = x + 1.0d0
       gam1_gas  = chit_gas*x + chid_gas
       nabad_gas = x/gam1_gas
       gam2_gas  = 1.0d0/(1.0d0 - nabad_gas)
       cp_gas    = cv_gas * gam1_gas/chid_gas
       z         = 1.0d0 + (egas + light2)*zzi
       sound_gas = clight * sqrt(gam1_gas/z)



! for the totals
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



! maxwell relations; each is zero if the consistency is perfect
       x   = den * den

       dse = temp*dentrdt/denerdt - 1.0d0

       dpe = (denerdd*x + temp*dpresdt)/pres - 1.0d0

       dsp = -dentrdd*x/dpresdt - 1.0d0


! store this row
        ptot_row(j)   = pres
        dpt_row(j)    = dpresdt
        dpd_row(j)    = dpresdd
        dpa_row(j)    = dpresda
        dpz_row(j)    = dpresdz

        etot_row(j)   = ener
        det_row(j)    = denerdt
        ded_row(j)    = denerdd
        dea_row(j)    = denerda
        dez_row(j)    = denerdz

        stot_row(j)   = entr
        dst_row(j)    = dentrdt
        dsd_row(j)    = dentrdd
        dsa_row(j)    = dentrda
        dsz_row(j)    = dentrdz


        pgas_row(j)   = pgas
        dpgast_row(j) = dpgasdt
        dpgasd_row(j) = dpgasdd
        dpgasa_row(j) = dpgasda
        dpgasz_row(j) = dpgasdz

        egas_row(j)   = egas
        degast_row(j) = degasdt
        degasd_row(j) = degasdd
        degasa_row(j) = degasda
        degasz_row(j) = degasdz

        sgas_row(j)   = sgas
        dsgast_row(j) = dsgasdt
        dsgasd_row(j) = dsgasdd
        dsgasa_row(j) = dsgasda
        dsgasz_row(j) = dsgasdz


        prad_row(j)   = prad
        dpradt_row(j) = dpraddt
        dpradd_row(j) = dpraddd
        dprada_row(j) = dpradda
        dpradz_row(j) = dpraddz

        erad_row(j)   = erad
        deradt_row(j) = deraddt
        deradd_row(j) = deraddd
        derada_row(j) = deradda
        deradz_row(j) = deraddz

        srad_row(j)   = srad
        dsradt_row(j) = dsraddt
        dsradd_row(j) = dsraddd
        dsrada_row(j) = dsradda
        dsradz_row(j) = dsraddz


        pion_row(j)   = pion
        dpiont_row(j) = dpiondt
        dpiond_row(j) = dpiondd
        dpiona_row(j) = dpionda
        dpionz_row(j) = dpiondz

        eion_row(j)   = eion
        deiont_row(j) = deiondt
        deiond_row(j) = deiondd
        deiona_row(j) = deionda
        deionz_row(j) = deiondz

        sion_row(j)   = sion
        dsiont_row(j) = dsiondt
        dsiond_row(j) = dsiondd
        dsiona_row(j) = dsionda
        dsionz_row(j) = dsiondz

        xni_row(j)    = xni

        pele_row(j)   = pele
        ppos_row(j)   = 0.0d0
        dpept_row(j)  = dpepdt
        dpepd_row(j)  = dpepdd
        dpepa_row(j)  = dpepda
        dpepz_row(j)  = dpepdz

        eele_row(j)   = eele
        epos_row(j)   = 0.0d0
        deept_row(j)  = deepdt
        deepd_row(j)  = deepdd
        deepa_row(j)  = deepda
        deepz_row(j)  = deepdz

        sele_row(j)   = sele
        spos_row(j)   = 0.0d0
        dsept_row(j)  = dsepdt
        dsepd_row(j)  = dsepdd
        dsepa_row(j)  = dsepda
        dsepz_row(j)  = dsepdz

        xnem_row(j)   = xnem
        xne_row(j)    = xnefer
        dxnet_row(j)  = dxnedt
        dxned_row(j)  = dxnedd
        dxnea_row(j)  = dxneda
        dxnez_row(j)  = dxnedz
        xnp_row(j)    = 0.0d0
        zeff_row(j)   = zbar

        etaele_row(j) = etaele
        detat_row(j)  = detadt
        detad_row(j)  = detadd
        detaa_row(j)  = detada
        detaz_row(j)  = detadz
        etapos_row(j) = 0.0d0

        pcou_row(j)   = pcoul
        dpcout_row(j) = dpcouldt
        dpcoud_row(j) = dpcouldd
        dpcoua_row(j) = dpcoulda
        dpcouz_row(j) = dpcouldz

        ecou_row(j)   = ecoul
        decout_row(j) = decouldt
        decoud_row(j) = decouldd
        decoua_row(j) = decoulda
        decouz_row(j) = decouldz

        scou_row(j)   = scoul
        dscout_row(j) = dscouldt
        dscoud_row(j) = dscouldd
        dscoua_row(j) = dscoulda
        dscouz_row(j) = dscouldz

        plasg_row(j)  = plasg

        dse_row(j)    = dse
        dpe_row(j)    = dpe
        dsp_row(j)    = dsp

        cv_gas_row(j)    = cv_gas
        cp_gas_row(j)    = cp_gas
        gam1_gas_row(j)  = gam1_gas
        gam2_gas_row(j)  = gam2_gas
        gam3_gas_row(j)  = gam3_gas
        nabad_gas_row(j) = nabad_gas
        cs_gas_row(j)    = sound_gas

        cv_row(j)     = cv
        cp_row(j)     = cp
        gam1_row(j)   = gam1
        gam2_row(j)   = gam2
        gam3_row(j)   = gam3
        nabad_row(j)  = nabad
        cs_row(j)     = sound

! end of pipeline loop
      enddo
      return
      end
!---------------------------------------------------------------------







!---------------------------------------------------------------------
! this file contains routines to invert the helmholtz eos
!
! routine invert_helm_pt is used when the pressure and temperature are given
! routine invert_helm_pt_quiet is as above, but supresses all error messages
! routine invert_helm_pd is used when the pressure and density are given
! routine invert_helm_et is used when the energy and temperature are given
! routine invert_helm_ed is used when the energy and density are given
! routine invert_helm_st is used when the entropy and temperature are given
! routine invert_helm_st_quiet is as above, but supresses all error messages



      subroutine invert_helm_pt
      include 'implno.dek'
      include 'const.dek'
      include 'vector_eos.dek'


! given the pressure, temperature, and composition
! find everything else

! it is assumed that ptot_row(j), temp_row(j), abar_row(j),
! zbar_row(j), and the pipe limits (jlo_eos:jhi_eos), have
! been set before calling this routine.

! on input den_row(j) conatins a guess for the density,
! on output den_row(j) contains the converged density.

! To get the greatest speed advantage, the eos should be fed a
! large pipe of data to work on.


! local variables
      integer          i,j,jlo_save,jhi_save
      double precision den,f,df,dennew,eostol,fpmin
      parameter        (eostol = 1.0d-8, &
                        fpmin  = 1.0d-14)


! initialize
      jlo_save = jlo_eos
      jhi_save = jhi_eos
      do j=jlo_eos, jhi_eos
       eoswrk01(j) = 0.0d0
       eoswrk02(j) = 0.0d0
       eoswrk03(j) = ptot_row(j)
       eoswrk04(j) = den_row(j)
      end do


! do the first newton loop with all elements in the pipe
!      call helmeos
      call helmeosaprox13

      do j = jlo_eos, jhi_eos

       f     = ptot_row(j)/eoswrk03(j) - 1.0d0
       df    = dpd_row(j)/eoswrk03(j)
       eoswrk02(j) = f/df

! limit excursions to factor of two changes
       den    = den_row(j)
       dennew = min(max(0.5d0*den,den - eoswrk02(j)),2.0d0*den)

! compute the error
       eoswrk01(j)  = abs((dennew - den)/den)

! store the new density, keep it within the table limits
       den_row(j)  = min(1.0d14,max(dennew,1.0d-11))
      enddo



! now loop over each element of the pipe individually
      do j = jlo_save, jhi_save

       do i=2,40

        if (eoswrk01(j) .lt. eostol .or. &
            abs(eoswrk02(j)) .le. fpmin) goto 20

        jlo_eos = j
        jhi_eos = j

!        call helmeos
        call helmeosaprox13

        f     = ptot_row(j)/eoswrk03(j) - 1.0d0
        df    = dpd_row(j)/eoswrk03(j)
        eoswrk02(j) = f/df

! limit excursions to factor of two changes
        den    = den_row(j)
        dennew = min(max(0.5d0*den,den - eoswrk02(j)),2.0d0*den)

! compute the error
        eoswrk01(j)  = abs((dennew - den)/den)

! store the new density, keep it within the table limits
        den_row(j)  = min(1.0d14,max(dennew,1.0d-11))

! end of netwon loop
       end do


! we did not converge if we land here
      write(6,*)
      write(6,*) 'newton-raphson failed in routine invert_helm_pt'
      write(6,*) 'pipeline element',j
      write(6,01) 'pwant  =',eoswrk03(j),' temp =',temp_row(j)
 01   format(1x,5(a,1pe16.8))
      write(6,01) 'error =',eoswrk01(j), &
                  '  eostol=',eostol,'  fpmin =',fpmin
      write(6,01) 'den   =',den_row(j),'  denold=',eoswrk04(j)
      write(6,01) 'f/df  =',eoswrk02(j),' f   =',f,    ' df    =',df
      write(6,*)
      stop 'could not find a density in routine invert_helm_pt'



! land here if newton loop converged, back for another pipe element
 20    continue
      end do



! call eos one more time with the converged value of the density

      jlo_eos = jlo_save
      jhi_eos = jhi_save

!      call helmeos
      call helmeosaprox13

      return
      end







      subroutine invert_helm_pt_quiet
      include 'implno.dek'
      include 'const.dek'
      include 'vector_eos.dek'


! given the pressure, temperature, and composition
! find everything else

! it is assumed that ptot_row(j), temp_row(j), abar_row(j),
! zbar_row(j), and the pipe limits (jlo_eos:jhi_eos), have
! been set before calling this routine.

! on input den_row(j) conatins a guess for the density,
! on output den_row(j) contains the converged density.

! To get the greatest speed advantage, the eos should be fed a
! large pipe of data to work on.

! this version is quiet on all errors


! local variables
      integer          i,j,jlo_save,jhi_save
      double precision den,f,df,dennew,eostol,fpmin
      parameter        (eostol = 1.0d-8, &
                        fpmin  = 1.0d-14)


! initialize
      jlo_save = jlo_eos
      jhi_save = jhi_eos
      do j=jlo_eos, jhi_eos
       eoswrk01(j) = 0.0d0
       eoswrk02(j) = 0.0d0
       eoswrk03(j) = ptot_row(j)
       eoswrk04(j) = den_row(j)
      end do


! do the first newton loop with all elements in the pipe
!      call helmeos
      call helmeosaprox13

      do j = jlo_eos, jhi_eos

       f     = ptot_row(j)/eoswrk03(j) - 1.0d0
       df    = dpd_row(j)/eoswrk03(j)
       eoswrk02(j) = f/df

! limit excursions to factor of two changes
       den    = den_row(j)
       dennew = min(max(0.5d0*den,den - eoswrk02(j)),2.0d0*den)

! compute the error
       eoswrk01(j)  = abs((dennew - den)/den)

! store the new density, keep it within the table limits
       den_row(j)  = min(1.0d14,max(dennew,1.0d-11))

      enddo



! now loop over each element of the pipe individually
      do j = jlo_save, jhi_save

       do i=2,40

        if (eoswrk01(j) .lt. eostol .or. &
            abs(eoswrk02(j)) .le. fpmin) goto 20

        jlo_eos = j
        jhi_eos = j

!        call helmeos
        call helmeosaprox13

        f     = ptot_row(j)/eoswrk03(j) - 1.0d0
        df    = dpd_row(j)/eoswrk03(j)
        eoswrk02(j) = f/df

! limit excursions to factor of two changes
        den    = den_row(j)
        dennew = min(max(0.5d0*den,den - eoswrk02(j)),2.0d0*den)

! compute the error
        eoswrk01(j)  = abs((dennew - den)/den)

! store the new density, keep it within the table limits
        den_row(j)  = min(1.0d14,max(dennew,1.0d-11))

! end of netwon loop
       end do


! we did not converge if we land here
!      write(6,*)
!      write(6,*) 'newton-raphson failed in routine invert_helm_pt'
!      write(6,*) 'pipeline element',j
!      write(6,01) 'pres  =',eoswrk03(j)
! 01   format(1x,5(a,1pe16.8))
!      write(6,01) 'error =',eoswrk01(j),
!     1            '  eostol=',eostol,'  fpmin =',fpmin
!      write(6,01) 'den   =',den_row(j),'  denold=',eoswrk04(j)
!      write(6,01) 'f/df  =',eoswrk02(j),' f   =',f,    ' df    =',df
!      write(6,*)
!      stop 'could not find a density in routine invert_helm_pt'



! land here if newton loop converged, back for another pipe element
 20    continue
      end do



! call eos one more time with the converged value of the density

      jlo_eos = jlo_save
      jhi_eos = jhi_save

!      call helmeos
      call helmeosaprox13

      return
      end







      subroutine invert_helm_pd
      include 'implno.dek'
      include 'const.dek'
      include 'vector_eos.dek'


! given the pressure, density, and composition
! find everything else

! it is assumed that ptot_row(j), den_row(j), abar_row(j),
! zbar_row(j), and the pipe limits (jlo_eos:jhi_eos), have
! been set before calling this routine.

! on input temp_row(j) conatins a guess for the temperature,
! on output temp_row(j) contains the converged temperature.

! To get the greatest speed advantage, the eos should be fed a
! large pipe of data to work on.


! local variables
      integer          i,j,jlo_save,jhi_save
      double precision tmpold,tmp,f,df,tmpnew,eostol,fpmin
      parameter        (eostol = 1.0d-8, &
                        fpmin  = 1.0d-14)


! initialize
      jlo_save = jlo_eos
      jhi_save = jhi_eos
      do j=jlo_eos, jhi_eos
       eoswrk01(j) = 0.0d0
       eoswrk02(j) = 0.0d0
       eoswrk03(j) = ptot_row(j)
       eoswrk04(j) = temp_row(j)
      end do


! do the first newton loop with all elements in the pipe
!      call helmeos
      call helmeosaprox13

      do j = jlo_eos, jhi_eos

       f     = ptot_row(j)/eoswrk03(j) - 1.0d0
       df    = dpt_row(j)/eoswrk03(j)
       eoswrk02(j) = f/df

! limit excursions to factor of two changes
       tmp    = temp_row(j)
       tmpnew = min(max(0.5d0*tmp,tmp - eoswrk02(j)),2.0d0*tmp)

! compute the error
       eoswrk01(j)  = abs((tmpnew - tmp)/tmp)

! store the new density, keep it within the table limits
       temp_row(j)  = min(1.0d14,max(tmpnew,1.0d-11))

      enddo



! now loop over each element of the pipe individually
      do j = jlo_save, jhi_save

       do i=2,40

        if (eoswrk01(j) .lt. eostol .or. &
            abs(eoswrk02(j)) .le. fpmin) goto 20

        jlo_eos = j
        jhi_eos = j

!        call helmeos
        call helmeosaprox13

        f     = ptot_row(j)/eoswrk03(j) - 1.0d0
        df    = dpt_row(j)/eoswrk03(j)
        eoswrk02(j) = f/df

! limit excursions to factor of two changes
        tmp    = temp_row(j)
        tmpnew = min(max(0.5d0*tmp,tmp - eoswrk02(j)),2.0d0*tmp)

! compute the error
        eoswrk01(j)  = abs((tmpnew - tmp)/tmp)

! store the new density, keep it within the table limits
        temp_row(j)  = min(1.0d13,max(tmpnew,1.0d3))

! end of netwon loop
       end do


! we did not converge if we land here
      write(6,*)
      write(6,*) 'newton-raphson failed in routine invert_helm_pd'
      write(6,*) 'pipeline element',j
      write(6,01) 'pwant  =',eoswrk03(j)
 01   format(1x,5(a,1pe16.8))
      write(6,01) 'error =',eoswrk01(j), &
                  '  eostol=',eostol,'  fpmin =',fpmin
      write(6,01) 'tmp   =',temp_row(j),'  tmpold=',eoswrk04(j)
      write(6,01) 'f/df  =',eoswrk02(j),' f   =',f,    ' df    =',df
      write(6,*)
      stop 'could not find a temperature in routine invert_helm_pd'



! land here if newton loop converged, back for another pipe element
 20    continue
      end do



! call eos one more time with the converged value of the density

      jlo_eos = jlo_save
      jhi_eos = jhi_save

!      call helmeos
      call helmeosaprox13

      return
      end






      subroutine invert_helm_pd_quiet
      include 'implno.dek'
      include 'const.dek'
      include 'vector_eos.dek'


! given the pressure, density, and composition
! find everything else

! it is assumed that ptot_row(j), den_row(j), abar_row(j),
! zbar_row(j), and the pipe limits (jlo_eos:jhi_eos), have
! been set before calling this routine.

! on input temp_row(j) conatins a guess for the temperature,
! on output temp_row(j) contains the converged temperature.

! To get the greatest speed advantage, the eos should be fed a
! large pipe of data to work on.


! local variables
      integer          i,j,jlo_save,jhi_save
      double precision tmpold,tmp,f,df,tmpnew,eostol,fpmin
      parameter        (eostol = 1.0d-8, &
                        fpmin  = 1.0d-14)


! initialize
      jlo_save = jlo_eos
      jhi_save = jhi_eos
      do j=jlo_eos, jhi_eos
       eoswrk01(j) = 0.0d0
       eoswrk02(j) = 0.0d0
       eoswrk03(j) = ptot_row(j)
       eoswrk04(j) = temp_row(j)
      end do


! do the first newton loop with all elements in the pipe
!      call helmeos
      call helmeosaprox13

      do j = jlo_eos, jhi_eos

       f     = ptot_row(j)/eoswrk03(j) - 1.0d0
       df    = dpt_row(j)/eoswrk03(j)
       eoswrk02(j) = f/df

! limit excursions to factor of two changes
       tmp    = temp_row(j)
       tmpnew = min(max(0.5d0*tmp,tmp - eoswrk02(j)),2.0d0*tmp)

! compute the error
       eoswrk01(j)  = abs((tmpnew - tmp)/tmp)

! store the new density, keep it within the table limits
       temp_row(j)  = min(1.0d14,max(tmpnew,1.0d-11))

      enddo



! now loop over each element of the pipe individually
      do j = jlo_save, jhi_save

       do i=2,40

        if (eoswrk01(j) .lt. eostol .or. &
            abs(eoswrk02(j)) .le. fpmin) goto 20

        jlo_eos = j
        jhi_eos = j

!        call helmeos
        call helmeosaprox13

        f     = ptot_row(j)/eoswrk03(j) - 1.0d0
        df    = dpt_row(j)/eoswrk03(j)
        eoswrk02(j) = f/df

! limit excursions to factor of two changes
        tmp    = temp_row(j)
        tmpnew = min(max(0.5d0*tmp,tmp - eoswrk02(j)),2.0d0*tmp)

! compute the error
        eoswrk01(j)  = abs((tmpnew - tmp)/tmp)

! store the new density, keep it within the table limits
        temp_row(j)  = min(1.0d13,max(tmpnew,1.0d3))

! end of netwon loop
       end do


! we did not converge if we land here
!      write(6,*)
!      write(6,*) 'newton-raphson failed in routine invert_helm_pd'
!      write(6,*) 'pipeline element',j
!      write(6,01) 'pwant  =',eoswrk03(j)
! 01   format(1x,5(a,1pe16.8))
!      write(6,01) 'error =',eoswrk01(j),
!     1            '  eostol=',eostol,'  fpmin =',fpmin
!      write(6,01) 'tmp   =',temp_row(j),'  tmpold=',eoswrk04(j)
!      write(6,01) 'f/df  =',eoswrk02(j),' f   =',f,    ' df    =',df
!      write(6,*)
!      stop 'could not find a temperature in routine invert_helm_pd'



! land here if newton loop converged, back for another pipe element
 20    continue
      end do



! call eos one more time with the converged value of the density

      jlo_eos = jlo_save
      jhi_eos = jhi_save

!      call helmeos
      call helmeosaprox13

      return
      end






      subroutine invert_helm_et
      include 'implno.dek'
      include 'const.dek'
      include 'vector_eos.dek'


! given the specific internal energy, temperature, and composition,
! find everything else

! it is assumed that etot_row(j), temp_row(j), abar_row(j),
! zbar_row(j), and the pipe limits (jlo_eos:jhi_eos), have
! been set before calling this routine.

! on input den_row(j) conatins a guess for the density,
! on output den_row(j) contains the converged density.

! To get the greatest speed advantage, the eos should be fed a
! large pipe of data to work on.


! local variables
      integer          i,j,jlo_save,jhi_save
      double precision den,f,df,dennew,eostol,fpmin
      parameter        (eostol = 1.0d-8, &
                        fpmin  = 1.0d-14)


! initialize
      jlo_save = jlo_eos
      jhi_save = jhi_eos
      do j=jlo_eos, jhi_eos
       eoswrk01(j) = 0.0d0
       eoswrk02(j) = 0.0d0
       eoswrk03(j) = etot_row(j)
       eoswrk04(j) = den_row(j)
      end do


! do the first newton loop with all elements in the pipe
!      call helmeos
      call helmeosaprox13

      do j = jlo_eos, jhi_eos

       f     = etot_row(j)/eoswrk03(j) - 1.0d0
       df    = ded_row(j)/eoswrk03(j)
       eoswrk02(j) = f/df

! limit excursions to factor of two changes
       den    = den_row(j)
       dennew = min(max(0.5d0*den,den - eoswrk02(j)),2.0d0*den)

! compute the error
       eoswrk01(j)  = abs((dennew - den)/den)

! store the new density, keep it within the table limits
       den_row(j)  = min(1.0d14,max(dennew,1.0d-11))

      enddo



! now loop over each element of the pipe individually
      do j = jlo_save, jhi_save

       do i=2,40

        if (eoswrk01(j) .lt. eostol .or. &
            abs(eoswrk02(j)) .le. fpmin) goto 20

        jlo_eos = j
        jhi_eos = j

!        call helmeos
        call helmeosaprox13

        f     = etot_row(j)/eoswrk03(j) - 1.0d0
        df    = ded_row(j)/eoswrk03(j)
        eoswrk02(j) = f/df

! limit excursions to factor of two changes
        den    = den_row(j)
        dennew = min(max(0.5d0*den,den - eoswrk02(j)),2.0d0*den)

! compute the error
        eoswrk01(j)  = abs((dennew - den)/den)

! store the new density, keep it within the table limits
        den_row(j)  = min(1.0d14,max(dennew,1.0d-11))

! end of netwon loop
       end do


! we did not converge if we land here
      write(6,*)
      write(6,*) 'newton-raphson failed in routine invert_helm_et'
      write(6,*) 'pipeline element',j
      write(6,01) 'ewant  =',eoswrk03(j)
 01   format(1x,5(a,1pe16.8))
      write(6,01) 'error =',eoswrk01(j), &
                  '  eostol=',eostol,'  fpmin =',fpmin
      write(6,01) 'den   =',den_row(j),'  denold=',eoswrk04(j)
      write(6,01) 'f/df  =',eoswrk02(j),' f   =',f,    ' df    =',df
      write(6,*)
      stop 'could not find a density in routine invert_helm_et'



! land here if newton loop converged, back for another pipe element
 20    continue
      end do



! call eos one more time with the converged value of the density

      jlo_eos = jlo_save
      jhi_eos = jhi_save

!      call helmeos
      call helmeosaprox13

      return
      end





      subroutine invert_helm_ed
      include 'implno.dek'
      include 'const.dek'
      include 'vector_eos.dek'


! given the specific internal energy density, density, and composition
! find everything else

! it is assumed that etot_row(j), den_row(j), abar_row(j),
! zbar_row(j), and the pipe limits (jlo_eos:jhi_eos), have
! been set before calling this routine.

! on input temp_row(j) conatins a guess for the temperature,
! on output temp_row(j) contains the converged temperature.

! To get the greatest speed advantage, the eos should be fed a
! large pipe of data to work on.


! local variables
      integer          i,j,jlo_save,jhi_save
      double precision tmpold,tmp,f,df,tmpnew,eostol,fpmin
      parameter        (eostol = 1.0d-8, &
                        fpmin  = 1.0d-14)



! initialize
      jlo_save = jlo_eos
      jhi_save = jhi_eos
      do j=jlo_eos, jhi_eos
       eoswrk01(j) = 0.0d0
       eoswrk02(j) = 0.0d0
       eoswrk03(j) = etot_row(j)
       eoswrk04(j) = temp_row(j)
      end do


! do the first newton loop with all elements in the pipe

!      call helmeos
      call helmeosaprox13

      do j = jlo_eos, jhi_eos

       f     = etot_row(j)/eoswrk03(j) - 1.0d0
       df    = det_row(j)/eoswrk03(j)

       eoswrk02(j) = f/df

! limit excursions to factor of two changes
       tmp    = temp_row(j)
       tmpnew = min(max(0.5d0*tmp,tmp - eoswrk02(j)),2.0d0*tmp)

! compute the error
       eoswrk01(j)  = abs((tmpnew - tmp)/tmp)

! store the new temperature, keep it within the table limits
       temp_row(j)  = min(1.0d14,max(tmpnew,1.0d-11))


      enddo



! now loop over each element of the pipe individually
      do j = jlo_save, jhi_save

       do i=2,40

        if (eoswrk01(j) .lt. eostol .or. &
            abs(eoswrk02(j)) .le. fpmin) goto 20

        jlo_eos = j
        jhi_eos = j

!        call helmeos
        call helmeosaprox13

        f     = etot_row(j)/eoswrk03(j) - 1.0d0
        df    = det_row(j)/eoswrk03(j)

        eoswrk02(j) = f/df

! limit excursions to factor of two changes
        tmp    = temp_row(j)
        tmpnew = min(max(0.5d0*tmp,tmp - eoswrk02(j)),2.0d0*tmp)

! compute the error
        eoswrk01(j)  = abs((tmpnew - tmp)/tmp)

! store the new density, keep it within the table limits
        temp_row(j)  = min(1.0d13,max(tmpnew,1.0d3))


! end of netwon loop
       end do


! we did not converge if we land here
      write(6,*)
      write(6,*) 'newton-raphson failed in routine invert_helm_ed'
      write(6,*) 'pipeline element',j
      write(6,01) 'ewant  =',eoswrk03(j)
 01   format(1x,5(a,1pe16.8))
      write(6,01) 'error =',eoswrk01(j), &
                  '  eostol=',eostol,'  fpmin =',fpmin
      write(6,01) 'tmp   =',temp_row(j),'  tmpold=',eoswrk04(j)
      write(6,01) 'f/df  =',eoswrk02(j),' f   =',f,    ' df    =',df
      write(6,*)
      stop 'could not find a temperature in routine invert_helm_ed'



! land here if newton loop converged, back for another pipe element
 20    continue
      end do



! call eos one more time with the converged value of the density

      jlo_eos = jlo_save
      jhi_eos = jhi_save

!      call helmeos
      call helmeosaprox13

      return
      end





      subroutine invert_helm_st
      include 'implno.dek'
      include 'const.dek'
      include 'vector_eos.dek'


! given the entropy, temperature, and composition
! find everything else

! it is assumed that stot_row(j), temp_row(j), abar_row(j),
! zbar_row(j), and the pipe limits (jlo_eos:jhi_eos), have
! been set before calling this routine.

! on input den_row(j) conatins a guess for the density,
! on output den_row(j) contains the converged density.

! To get the greatest speed advantage, the eos should be fed a
! large pipe of data to work on.

! this version is quiet on all errors


! local variables
      integer          i,j,jlo_save,jhi_save
      double precision den,f,df,dennew,eostol,fpmin
      parameter        (eostol = 1.0d-8, &
                        fpmin  = 1.0d-14)


! initialize
      jlo_save = jlo_eos
      jhi_save = jhi_eos
      do j=jlo_eos, jhi_eos
       eoswrk01(j) = 0.0d0
       eoswrk02(j) = 0.0d0
       eoswrk03(j) = stot_row(j)
       eoswrk04(j) = den_row(j)
      end do


! do the first newton loop with all elements in the pipe
!      call helmeos
      call helmeosaprox13

      do j = jlo_eos, jhi_eos

       f     = stot_row(j)/eoswrk03(j) - 1.0d0
       df    = dsd_row(j)/eoswrk03(j)
       eoswrk02(j) = f/df

! limit excursions to factor of two changes
       den    = den_row(j)
       dennew = min(max(0.5d0*den,den - eoswrk02(j)),2.0d0*den)

! compute the error
       eoswrk01(j)  = abs((dennew - den)/den)

! store the new density, keep it within the table limits
       den_row(j)  = min(1.0d14,max(dennew,1.0d-11))

      enddo



! now loop over each element of the pipe individually
      do j = jlo_save, jhi_save

       do i=2,40

        if (eoswrk01(j) .lt. eostol .or. &
            abs(eoswrk02(j)) .le. fpmin) goto 20

        jlo_eos = j
        jhi_eos = j

!        call helmeos
        call helmeosaprox13

        f     = stot_row(j)/eoswrk03(j) - 1.0d0
        df    = dsd_row(j)/eoswrk03(j)
        eoswrk02(j) = f/df

! limit excursions to factor of two changes
        den    = den_row(j)
        dennew = min(max(0.5d0*den,den - eoswrk02(j)),2.0d0*den)

! compute the error
        eoswrk01(j)  = abs((dennew - den)/den)

! store the new density, keep it within the table limits
        den_row(j)  = min(1.0d14,max(dennew,1.0d-11))

! end of netwon loop
       end do


! we did not converge if we land here
      write(6,*)
      write(6,*) 'newton-raphson failed in routine invert_helm_st'
      write(6,*) 'pipeline element',j
      write(6,01) 'entr  =',eoswrk03(j)
 01   format(1x,5(a,1pe16.8))
      write(6,01) 'error =',eoswrk01(j), &
                  '  eostol=',eostol,'  fpmin =',fpmin
      write(6,01) 'den   =',den_row(j),'  denold=',eoswrk04(j)
      write(6,01) 'f/df  =',eoswrk02(j),' f   =',f,    ' df    =',df
      write(6,*)
      stop 'could not find a density in routine invert_helm_st'



! land here if newton loop converged, back for another pipe element
 20    continue
      end do



! call eos one more time with the converged value of the density

      jlo_eos = jlo_save
      jhi_eos = jhi_save

!      call helmeos
      call helmeosaprox13

      return
      end



      subroutine invert_helm_st_quiet
      include 'implno.dek'
      include 'const.dek'
      include 'vector_eos.dek'


! given the entropy, temperature, and composition
! find everything else

! it is assumed that stot_row(j), temp_row(j), abar_row(j),
! zbar_row(j), and the pipe limits (jlo_eos:jhi_eos), have
! been set before calling this routine.

! on input den_row(j) conatins a guess for the density,
! on output den_row(j) contains the converged density.

! To get the greatest speed advantage, the eos should be fed a
! large pipe of data to work on.

! this version is quiet on all errors


! local variables
      integer          i,j,jlo_save,jhi_save
      double precision den,f,df,dennew,eostol,fpmin
      parameter        (eostol = 1.0d-8, &
                        fpmin  = 1.0d-14)


! initialize
      jlo_save = jlo_eos
      jhi_save = jhi_eos
      do j=jlo_eos, jhi_eos
       eoswrk01(j) = 0.0d0
       eoswrk02(j) = 0.0d0
       eoswrk03(j) = stot_row(j)
       eoswrk04(j) = den_row(j)
      end do


! do the first newton loop with all elements in the pipe
!      call helmeos
      call helmeosaprox13

      do j = jlo_eos, jhi_eos

       f     = stot_row(j)/eoswrk03(j) - 1.0d0
       df    = dsd_row(j)/eoswrk03(j)
       eoswrk02(j) = f/df

! limit excursions to factor of two changes
       den    = den_row(j)
       dennew = min(max(0.5d0*den,den - eoswrk02(j)),2.0d0*den)

! compute the error
       eoswrk01(j)  = abs((dennew - den)/den)

! store the new density, keep it within the table limits
       den_row(j)  = min(1.0d14,max(dennew,1.0d-11))

      enddo



! now loop over each element of the pipe individually
      do j = jlo_save, jhi_save

       do i=2,40

        if (eoswrk01(j) .lt. eostol .or. &
            abs(eoswrk02(j)) .le. fpmin) goto 20

        jlo_eos = j
        jhi_eos = j

!        call helmeos
        call helmeosaprox13

        f     = stot_row(j)/eoswrk03(j) - 1.0d0
        df    = dsd_row(j)/eoswrk03(j)
        eoswrk02(j) = f/df

! limit excursions to factor of two changes
        den    = den_row(j)
        dennew = min(max(0.5d0*den,den - eoswrk02(j)),2.0d0*den)

! compute the error
        eoswrk01(j)  = abs((dennew - den)/den)

! store the new density, keep it within the table limits
        den_row(j)  = min(1.0d14,max(dennew,1.0d-11))

! end of netwon loop
       end do


! we did not converge if we land here
!      write(6,*)
!      write(6,*) 'newton-raphson failed in routine invert_helm_pt'
!      write(6,*) 'pipeline element',j
!      write(6,01) 'pres  =',eoswrk03(j)
! 01   format(1x,5(a,1pe16.8))
!      write(6,01) 'error =',eoswrk01(j),
!     1            '  eostol=',eostol,'  fpmin =',fpmin
!      write(6,01) 'den   =',den_row(j),'  denold=',eoswrk04(j)
!      write(6,01) 'f/df  =',eoswrk02(j),' f   =',f,    ' df    =',df
!      write(6,*)
!      stop 'could not find a density in routine invert_helm_pt'



! land here if newton loop converged, back for another pipe element
 20    continue
      end do



! call eos one more time with the converged value of the density

      jlo_eos = jlo_save
      jhi_eos = jhi_save

!      call helmeos
      call helmeosaprox13

      return
      end

!---------------------------------------------------------------------







!---------------------------------------------------------------------
!
! some system and glue utility routines
!
! routine today gets the date and clock time
! routine zsecond get the elasped cpu time
! routine lenstr finds the non-blank length of a string
! routine sqeeze compresses a string
! routine timlap converts total seconds into hours, minutes, seconds

! routine casedn converts a string to lower case
! routine uttoday gets the ut time out of a machine

! routine getnam is a basic string parser
! function value converts charcter strings to real numbers, stop on error
! function value2 converts charcter strings to real numbers, flag on error






      subroutine today(adat,atim)
      implicit none

! forms date and time strings

! declare the pass
      character*8  atim
      character*9  adat


! local variables
      character*3  amon(12)
      character*8  date
      character*10 time
      character*5  zone
      integer      idat(3),itim(3),values(8)

      data amon/ 'jan' , 'feb' , 'mar' , 'apr' , 'may' , 'jun' , &
                 'jul' , 'aug' , 'sep' , 'oct' , 'nov' , 'dec' /


! format statements for the time and date
113   format(i2.2,':',i2.2,':',i2.2)
114   format(i2.2,a3,i4.4)


! initialize
      adat=' '
      atim=' '


! f77 way, keep capitals for maximum portability
!      call ITIME(itim)
!      call IDATE(idat)
!      write(atim,113) itim
!      write(adat,114) idat(1),amon(idat(2)),idat(3)


! f90 way
!      call date_and_time(date,time,zone,values)
!      write(atim,113) values(5),values(6),values(7)
!      write(adat,114) values(3),amon(values(2)),values(1)


      return
      end






      subroutine zsecond(time)
!
! this routine gets the elapsed time of a job from the machine
!
! declare
!      external         ETIME
      double precision time
!      real             ETIME,tarray(2),ttt


! initialize
      time = 0.0d0

! f77 way, keep capitals for maximum portability
!      time = ETIME(tarray)

! f90 doesn't have a way to get the cputime, only wall clock time. grrr.

! f95 intrinsic
!      call cpu_time(ttt)
!      time = ttt


      return
      end





      integer function lenstr(string,istrln)
      include 'implno.dek'
!
! lenstr returns the non blank length length of the string.
!
! declare
      integer       istrln,i
      character*(*) string


      lenstr=0
      do i=istrln,1,-1
       if (string(i:i).ne. ' ') then
        if (ichar(string(i:i)).ne. 0 )then
         lenstr=i
         goto  20
        end if
       end if
      enddo
20    return
      end




      subroutine sqeeze(line)
      include 'implno.dek'
!
! this routine takes line and removes all blanks, such as
! those from writing to string with fortran format statements
!
! declare
      character*(*)  line
      character*1    achar
      integer        l,n,k,lend,lsiz,lenstr


! find the end of the line
      lsiz = len(line)
      lend = lenstr(line,lsiz)
      n    = 0
      l    = 0

! do the compression in place
10    continue
      l = l + 1
      achar = line(l:l)
      if (achar .eq. ' ') goto 10
      n = n + 1
      line(n:n) = achar
      if (l .lt. lend) goto 10

! blank the rest of the line
      do k=n+1,lsiz
       line(k:k) = ' '
      enddo
      return
      end





      subroutine timlap(tlap,hours,minut,sec,msec)
      include 'implno.dek'
!
! this routines converts seconds to hours, minutes, seconds and microseconds
!
! declare
      integer          hours,minut,sec,msec
      double precision tlap,x
      msec  = 0
      sec   = 0
      minut = 0
      hours = 0
      sec   = int(tlap)
      msec  = 1.0d6 * (tlap-sec)
      if (sec .ge. 60) then
       x  = dble(sec)/60.0d0
       minut = int(x)
      end if
      sec = sec - minut*60
      if (minut .ge. 60) then
       x = dble(minut)/60.0d0
       hours = int(x)
      end if
      minut = minut - hours*60
      return
      end






      subroutine casedn(string)
      implicit none
      save

! this routine converts an ascii string to all lower case.

! declare
      character*(*)  string
      integer        i,x,biga,bigz,change
      parameter     (biga = 65, bigz = 90, change = 32)

      do i=1,len(string)
       x = ichar(string(i:i))
       if (x .ge. biga  .and.  x .le. bigz ) then
        x           = x + change
        string(i:i) = char(x)
       end if
      enddo
      return
      end




      subroutine uttoday(iyear,imonth,iday,ihour,iminute,isecond)
      implicit none
      save

! this routine gets the UT date and time out of a machine.

! output:
! iyear   = integer year
! imonth  = integar month, between 1-12
! iday    = integer day of month, between 1-31
! ihour   = hours past midnight, between 0-23
! iminute = minutes past the hour, 0-59
! isecond = seconds past the minute, 0-59

! declare
      integer iyear,imonth,iday,ihour,iminute,isecond

! local variables
!      external      time
!      integer       time
!c      external      time_
!c      integer       time_
!      integer*4     stime,tarray(9)


      isecond = 0
      iminute = 0
      ihour   = 0
      iday    = 0
      imonth  = 0
      iyear   = 0


! system and compiler dependent
!      stime = TIME()
!      call GMTIME(stime,tarray)

! xlf mac
!      stime = time_()
!      call gmtime_(stime,tarray)

! intel mac
!      stime = time()
!      call gmtime(stime,tarray)

!      isecond = tarray(1)
!      iminute = tarray(2)
!      ihour   = tarray(3)
!      iday    = tarray(4)
!      imonth  = tarray(5) + 1
!      iyear   = tarray(6) + 1900

      return
      end






      integer function getnam(string,word,ipos)
      include 'implno.dek'

! this routine is a basic string parser, only spaces and equal signs are
! used as delimiters.

! input is the string to chop and where to start the parse from ipos.

! output is the parsed word and the updated position on the string ipos.
! the length of the word is returned as getnam. if no word is found,
! meaning that there are no more words in the string, then getnam is
! returned as zero and the word is filled with spaces.

! declare
      character*(*) string,word
      integer       kpos,ipos,k,lend,nbegin


! get the length of the input line, blank word and save ipos
      lend = len(string)
      word = ' '
      kpos = ipos


! find begining of the word; if nothing found return getnam as zero
      do k = kpos,lend
       nbegin = k
       if (string(k:k).ne.' ' .and. string(k:k).ne.'=') goto  25
      enddo
      getnam = 0
      return

! find end of the word
 25   continue
      do k = nbegin,lend
       ipos = k
       if( string(k:k).eq.' ' .or. string(k:k).eq.'=') goto  35
      enddo

! build the word, set getnam and return
 35   continue
      word   = string(nbegin:ipos-1)
      getnam = ipos - nbegin
      return
      end




      double precision function value(string)
      include 'implno.dek'

! this routine takes a character string and converts it to a real number.
! on error during the conversion, a fortran stop is issued

! declare
      logical          pflag
      character*(*)    string
      character*1      plus,minus,decmal,blank,se,sd,se1,sd1
      integer          noblnk,long,ipoint,power,psign,iten,j,z,i
      double precision x,sign,factor,rten,temp
      parameter        (plus = '+'  , minus = '-' , decmal = '.'   , &
                        blank = ' ' , se = 'e'    , sd = 'd'       , &
                        se1 = 'E'   , sd1 = 'D'   , rten =  10.0, &
                        iten = 10                                   )

! initialize
      x      =  0.0d0
      sign   =  1.0d0
      factor =  rten
      pflag  =  .false.
      noblnk =  0
      power  =  0
      psign  =  1
      long   =  len(string)


! remove any leading blanks and get the sign of the number
      do z = 1,7
       noblnk = noblnk + 1
       if ( string(noblnk:noblnk) .eq. blank) then
        if (noblnk .gt. 6 ) goto  30
       else
        if (string(noblnk:noblnk) .eq. plus) then
         noblnk = noblnk + 1
        else if (string(noblnk:noblnk) .eq. minus) then
         noblnk = noblnk + 1
         sign =  -1.0d0
        end if
        goto 10
       end if
      enddo


! main number conversion loop
 10   continue
      do i = noblnk,long
       ipoint = i + 1


! if a blank character then we are done
       if ( string(i:i) .eq. blank ) then
        x     = x * sign
        value = x
        return


! if an exponent character, process the whole exponent, and return
       else if (string(i:i).eq.se  .or. string(i:i).eq.sd .or. &
                string(i:i).eq.se1 .or. string(i:i).eq.sd1   ) then
        if (x .eq. 0.0 .and. ipoint.eq.2)     x = 1.0d0
        if (sign .eq. -1.0 .and. ipoint.eq.3) x = 1.0d0
        if (string(ipoint:ipoint) .eq. plus) ipoint = ipoint + 1
        if (string(ipoint:ipoint) .eq. minus) then
         ipoint = ipoint + 1
         psign = -1
        end if
        do z = ipoint,long
         if (string(z:z) .eq. blank)  then
          x = sign * x * rten**(power*psign)
          value = x
          return
         else
          j = ichar(string(z:z)) - 48
          if ( (j.lt.0) .or. (j.gt.9) ) goto 30
          power= (power * iten)  + j
         end if
        enddo


! if an ascii number character, process ie
       else if (string(i:i) .ne. decmal) then
        j = ichar(string(i:i)) - 48
        if ( (j.lt.0) .or. (j.gt.9) ) goto 30
        if (.not.(pflag) ) then
         x = (x*rten) + j
        else
         temp   = j
         x      = x + (temp/factor)
         factor = factor * rten
         goto 20
        end if

! must be a decimal point if none of the above
! check that there are not two decimal points!
       else
        if (pflag) goto 30
        pflag = .true.
       end if
 20    continue
      end do

! if we got through the do loop ok, then we must be done
      x     = x * sign
      value = x
      return


! error processing the number
 30   write(6,40) long,string(1:long)
 40   format(' error converting the ',i4,' characters ',/, &
             ' >',a,'< ',/, &
             ' into a real number in function value')
      stop ' error in routine value'
      end






      double precision function value2(string,err)
      include 'implno.dek'
!
! this routine takes a character string and converts it to a real number.
! on error during the conversion, the err flag is set to true
!
! declare
      logical          pflag,err
      character*(*)    string
      character*1      plus,minus,decmal,blank,se,sd,se1,sd1
      integer          noblnk,long,ipoint,power,psign,iten,j,z,i
      double precision x,sign,factor,rten,temp
      parameter        (plus = '+'  , minus = '-' , decmal = '.'   , &
                        blank = ' ' , se = 'e'    , sd = 'd'       , &
                        se1 = 'E'   , sd1 = 'D'   , rten =  10.0, &
                        iten = 10                                   )

! initialize
      err    = .false.
      x      =  0.0d0
      sign   =  1.0d0
      factor =  rten
      pflag  =  .false.
      noblnk =  0
      power  =  0
      psign  =  1
      long   =  len(string)

! remove any leading blanks and get the sign of the number
      do z = 1,7
       noblnk = noblnk + 1
       if ( string(noblnk:noblnk) .eq. blank) then
        if (noblnk .gt. 6 ) goto  30
       else
        if (string(noblnk:noblnk) .eq. plus) then
         noblnk = noblnk + 1
        else if (string(noblnk:noblnk) .eq. minus) then
         noblnk = noblnk + 1
         sign =  -1.0d0
        end if
        goto 10
       end if
      enddo

! main number conversion loop
10    do i = noblnk,long
       ipoint = i + 1

! if a blank character then we are done
       if ( string(i:i) .eq. blank ) then
        x     = x * sign
        value2 = x
        return

! if an exponent character, process the whole exponent, and return
       else if (string(i:i).eq.se  .or. string(i:i).eq.sd .or. &
                string(i:i).eq.se1 .or. string(i:i).eq.sd1   ) then
        if (x .eq. 0.0 .and. ipoint.eq.2)     x = 1.0d0
        if (sign .eq. -1.0 .and. ipoint.eq.3) x = 1.0d0
        if (string(ipoint:ipoint) .eq. plus) ipoint = ipoint + 1
        if (string(ipoint:ipoint) .eq. minus) then
         ipoint = ipoint + 1
         psign = -1
        end if
        do z = ipoint,long
         if (string(z:z) .eq. blank)  then
          x = sign * x * rten**(power*psign)
          value2 = x
          return
         else
          j = ichar(string(z:z)) - 48
          if ( (j.lt.0) .or. (j.gt.9) ) goto 30
          power= (power * iten)  + j
         end if
        enddo

! if an ascii number character, process ie
       else if (string(i:i) .ne. decmal) then
        j = ichar(string(i:i)) - 48
        if ( (j.lt.0) .or. (j.gt.9) ) goto 30
        if (.not.(pflag) ) then
         x = (x*rten) + j
        else
         temp   = j
         x      = x + (temp/factor)
         factor = factor * rten
         goto 20
        end if

! must be a decimal point if none of the above
! check that there are not two decimal points!
       else
        if (pflag) goto 30
        pflag = .true.
       end if
20    continue
      end do

! if we got through the do loop ok, then we must be done
      x      = x * sign
      value2 = x
      return

! error processing the number
30    err = .true.
      return
      end
!---------------------------------------------------------------------







!---------------------------------------------------------------------
      subroutine azbarapprox13(xmass,aion,zion,wion,ionmax, &
           ymass,abar,zbar,wbar,ye,nxcess)
      include 'implno.dek'

! this routine calculates composition variables

! input:
! mass fractions               = xmass(1:ionmax)  dimensionless
! number of nucleons           = aion(1:ionmax)   dimensionless
! charge of nucleus            = zion(1:ionmax)   dimensionless
! atomic weight or molar mass  = wion(1:ionmax)    g/mole
! number of isotopes           = ionmax
!
! output:
! molar abundances        = ymass(1:ionmax)   mole/g
! mean number of nucleons = abar              dimensionless
! mean nucleon charge     = zbar              dimensionless
! mean weight             = wbar              g/mole
! electron fraction       = ye                mole/g
! neutron excess          = xcess


! declare the pass
      integer          ionmax
      double precision xmass(ionmax),aion(ionmax),zion(ionmax), &
                       wion(ionmax),ymass(ionmax),abar,zbar,wbar, &
                       ye,nxcess


! local variables
      double precision sum1


! molar abundances
      ymass(1:ionmax) = xmass(1:ionmax)/wion(1:ionmax)

! mean molar mass
      wbar  = 1.0d0/sum(ymass(1:ionmax))

! mean number of nucleons
      sum1  = sum(aion(1:ionmax)*ymass(1:ionmax))
      abar  = wbar * sum1

! mean charge
      ye  = sum(zion(1:ionmax)*ymass(1:ionmax))
      zbar  = wbar * ye

! neutron excess
      nxcess = sum1 - 2.0d0 * ye

      return
      end




      subroutine azbar_simple(ymass,aion,zion,ionbeg,ionend, &
                              abar,zbar,ye)
      include 'implno.dek'

! this routine calculates composition variables

! input:
! molar abundances             = ymass(1:ionmax)   mole/g
! number of nucleons           = aion(1:ionmax)   dimensionless
! charge of nucleus            = zion(1:ionmax)   dimensionless
! start of isotopes            = ionbeg
! end of isotopes              = ionend
!
! output:
! mean number of nucleons = abar              dimensionless
! mean nucleon charge     = zbar              dimensionless
! electron fraction       = ye                mole/g


! declare the pass
      integer :: ionbeg,ionend
      real*8  :: ymass(ionbeg:ionend),aion(ionbeg:ionend),zion(ionbeg:ionend),abar,zbar,ye


! local variables
      double precision sum1

! mean number of nucleons, numerator assumes sum xmass = 1.0d0
      abar  = 1.0d0/sum(ymass(ionbeg:ionend))

! mean electron fraction and mean charge
      ye  = sum(zion(ionbeg:ionend)*ymass(ionbeg:ionend))
      zbar  = abar * ye

      return
      end



!---------------------------------------------------------------------







!---------------------------------------------------------------------
      subroutine nse(tt,dd,yye,newguess,ipartition,icoulomb, &
                     xmass_out,xmun,xmup,iprint)
      include 'implno.dek'
      include 'const.dek'
      include 'network.dek'


! this routine puts a chosen reaction network into its nse distribution.

! input:
! tt  = temperature
! dd  = density
! yye = electron mole number
! newguess   = 0 = use the old initial guess
!            = 1 = a new initial guess is made
! ipartition = 0 = use temperature independent partiction functions
!            = 1 = use temperature dependent partition functions
! icoulomb   = 0 = do not use coulomb corrections
!            = 1 = use coulomb corrections
! iprint = print flag


! output:
! xmass_out  = output nse mass fractions
! xmun        = chemical potential of neutrons
! xmup        = chemical potential of protons



! declare the pass
      integer :: newguess,ipartition,icoulomb,iprint
      real*8  :: tt,dd,yye,xmass_out(1),xmun,xmup


! communicate
      real*8 ::       temp,den,ye_want,beta,mu_c_p,yei(abignet),mu_c(abignet)
      common /nsec1/  temp,den,ye_want,beta,mu_c_p,yei,mu_c


! locals
      external           :: nsefunc
      logical            :: check
      integer            :: nfev,ntaken,i
      integer, parameter :: ntrial = 200, n = 2
      real*8             :: x(n),resid(n),amass,fac1,fac2,dum,xne,ge,a3,gi,sqrtgi
      real*8, parameter  :: tolf    = 1.0d-10,              tolx  = 1.0d-12, &
                            twopi   = 2.0d0*pi,             esqu  = qe*qe, &
                            forthpi = 4.0d0 * pi/3.0d0,     third = 1.0d0/3.0d0, &
                            fivth   = 5.0d0/3.0d0,          a1    = -0.9052d0, &
                            a2      = 0.6322d0,             a2inv = 1.0d0/a2, &
                            rt3     = 1.7320508075688772d0, half_rt3 = 0.5d0 * rt3

! for the initial guess
      real*8, parameter :: hinv = 1.0d0/h, mev2erg = ev2erg*1.0d6, &
                           mev2gr  = mev2erg/clight**2,  a56 = 56.0d0, &
                           z56 = 28.0d0, n56 = 28.0d0,  b56  = 4.8398d2


! initialize a common quantity
      yei(1:ionmax) = zion(1:ionmax)/aion(1:ionmax)


! fill the common block
      temp    = tt
      den     = dd
      ye_want = yye
      beta    = 1.0d0/(kerg * temp)


! get the partition functions
      wpart(ionbeg:ionmax) = 1.0d0
      if (ipartition .eq. 1) call get_partition_functions(tt*1.0d-9)


! set the coulomb corrections
! calder et al, apj 656 313 2007, eq a1

      mu_c_p = 0.0d0
      mu_c(ionbeg:ionmax) = 0.0d0
      if (icoulomb .eq. 1) then
       xne    = ye_want * avo * den
       ge     = esqu * beta * (forthpi * xne)**third
       a3     = -0.5d0*sqrt(3.0d0) - a1 / sqrt(a2)

       do i=1,ionmax
        gi      = zion(i)**(fivth)  * ge
        sqrtgi  = sqrt(gi)
        mu_c(i) = a1*(sqrt(gi*(a2+gi)) &
                   - a2*log(sqrt(gi*a2inv) + sqrt(1.0d0 + gi*a2inv)) &
                   + 2.0d0*a3*(sqrtgi - atan(sqrtgi)))
       enddo
       mu_c_p = mu_c(iprot)
      end if



! here is an initial guess for the neutron and proton chemical potentials,
! (x1) and (x2) respectively. obtained by setting xmass(ini56) = 1,
! setting mup = mun, and inverting the saha equation.
! here we use pure ni56 as it emprically appears to be a
! robust guess for all temp, rho, ye combinations.

      if (newguess .eq. 1) then
       newguess = 0
       amass  = n56*mn + z56*mp - b56*mev2gr
       fac1   = a56/(avo * den)
       fac2   = twopi/(beta*h) * amass*hinv
       fac2   = fac2 * sqrt(fac2)
       x(1)   = -(log(fac1*fac2)/beta + b56*ev2erg*1.0d6)/a56
       x(2)   = x(1)
      end if



! root find on mass and charge conservation for
! the chemical potentials of protons and neutrons

      call xnewt_nse(ntrial,x,n,tolx,tolf,ntaken,check,nfev,nsefunc)


! be sure we converged
      if (check .or. ntaken .eq. ntrial) then
       write(6,*)
       write(6,*) 'check convergence of root finder'
       write(6,*)
      end if


! some optional diagnostics
      if (iprint .eq. 1) then
       write(6,*)
       write(6,110) 'iterations taken             =',ntaken
       write(6,110) 'function evals               =',nfev
       write(6,111) 'roots                        =',x(1),x(2)
       call nsefunc(dum,x,resid)
       write(6,111) 'mass conservation   residual =',resid(1)
       write(6,111) 'charge conservation residual =',resid(2)

 110   format(1x,a,i4)
 111   format(1x,a,1p2e14.6)
      end if


! fill the output array using the converged values
      call nsefunc(dum,x,resid)


! bound the output nse abundances
      xmass_out(1:ionmax) = min(1.0d0,max(xmass_nse(1:ionmax),1.0d-30))
      xmun = x(1)
      xmup = x(2)

      return
      end







      subroutine nsefunc(x,y,f)
      include 'implno.dek'
      include 'const.dek'
      include 'network.dek'


! this routine returns the root find functions.
! input is the point x and y a vector of the unknowns.
! output is the vector of root find functions f, which should be the
! zero vector upon convergence.

! y(1) is input as the neutron chemical potential
! y(2) is input as the proton chemical potential


! declare the pass
      real*8 :: x,y(*),f(*)


! locals
      integer           :: i,indx(abignet),j,ifirst
      real*8            ::  ye,mu,fac1,fac2,fac3,xmsum,sum2,deninv,ww
      real*8, parameter :: hinv    = 1.0d0/h, twopih  = 2.0d0 * pi/h, &
                           mev2erg = ev2erg*1.0d6, mev2gr  = mev2erg/clight**2

! communicate
      real*8 ::       temp,den,ye_want,beta,mu_c_p,yei(abignet),mu_c(abignet)
      common /nsec1/  temp,den,ye_want,beta,mu_c_p,yei,mu_c


! chemical potential and mass fraction of each isotope
! hartmann et al, apj 297 837 1985, eq 2
! calder et al, apj 656 313 2007, eq a1

      deninv = 1.0d0/den
      ww     = twopih/beta

      do i=1,ionmax
       mu           = nion(i)*y(1) + zion(i)*y(2)
       fac1         = mion(i) * deninv * wpart(i)
       fac2         = ww * mion(i) * hinv
       fac2         = fac2*sqrt(fac2)
       fac3         = exp( beta * (mu + bion(i)*mev2erg) - mu_c(i) + zion(i)*mu_c_p)
       xmass_nse(i) = fac1 * fac2 * fac3
      enddo


! sum the mass fractions - mass conservation
! minimize roundoff by adding in ascending order

      call indexx(ionmax,xmass_nse,izwork1) 
      xmsum = 0.0d0
      do i=1,ionmax
       xmsum = xmsum + xmass_nse(izwork1(i))
      end do
!      xmsum = sum(xmass_nse(1:ionmax))


! form ye - charge conservation
! this formulation assumes mion(i) = aion(i)*amu

      zwork1(1:ionmax) = yei(1:ionmax) * xmass_nse(1:ionmax)
      call indexx(ionmax,zwork1,izwork1) 
      ye = 0.0d0
      do i=1,ionmax
       ye = ye + zwork1(izwork1(i))
      end do
!      ye = sum(yei(1:ionmax) * xmass_nse(1:ionmax))



! this formulation does not assume mion(i) = aion(i)*amu
! seitenzahl et al 2007

!      sum2 = sum(amu/mion(1:ionmax) * ((ye_want - 1.0d0)*zion(1:ionmax) &
!                 + ye_want * nion(1:ionmax)) * xmass_nse(1:ionmax) )


! mass and charge conservation are the requirements
      f(1) = xmsum - 1.0d0
      f(2) = ye - ye_want

!      f(2) = sum2

      return
      end





      subroutine nsejac(x,y,f,dfdy,n,np)
      include 'implno.dek'
      include 'const.dek'
      include 'network.dek'

! this routine returns the functions and the jacobian to do the root find on
! input is x, and y(n) a vector of the unknowns. output is f(n)
! and its jacobian dfdy(np,np).

! y(1) is the neutron chemical potential
! y(2) is the proton chemical potential


! declare the pass
      integer :: n,np
      real*8  :: x,y(n),f(n),dfdy(np,np)

! locals
      integer :: indx(abignet),i,j
      real*8  :: xmbn(abignet),xmbp(abignet),fac6(abignet), &
                 mu,mubn,mubp,fac1,fac2,fac3,fac4,fac5, &
                 xmsum,xmsumbn,xmsumbp,ye,yebn,yebp,deninv,ww,sum2,sum2bn,sum2bp
      real*8, parameter :: hinv  = 1.0d0/h, twopih  = 2.0d0 * pi/h, &
                           mev2erg = ev2erg*1.0d6, mev2gr  = mev2erg/clight**2

! communicate
      real*8 ::      temp,den,ye_want,beta,mu_c_p,yei(abignet),mu_c(abignet)
      common /nsec1/ temp,den,ye_want,beta,mu_c_p,yei,mu_c



! chemical potential and mass fraction of each isotope
! hartmann et al, apj 297 837 1985, eq 2
! calder et al, apj 656 313 2007, eq a1

      deninv = 1.0d0/den
      ww     = twopih/beta

! loop over isotopes
      do i=1,ionmax
       mu       = nion(i) * y(1) + zion(i) * y(2)
       mubn     = nion(i)
       mubp     = zion(i)

       fac1     = mion(i) * deninv * wpart(i)
       fac2     = ww * mion(i) * hinv
       fac2     = fac2 * sqrt(fac2)
       fac3     = exp( beta * (mu + bion(i) * ev2erg * 1.0d6) &
                           - mu_c(i) + zion(i)*mu_c_p)
       fac4     = fac1 * fac2 * fac3

       xmass_nse(i) = fac4
       xmbn(i)      = fac4 * beta * mubn
       xmbp(i)      = fac4 * beta * mubp

      enddo


! sum the mass fractions - mass conservation
      xmsum   = sum(xmass_nse(1:ionmax))
      xmsumbn = sum(xmbn(1:ionmax))
      xmsumbp = sum(xmbp(1:ionmax))


! form ye - charge conservation 
! this formulation assumes mion(i) = aion(i)*amu
      ye   = sum(yei(1:ionmax)*xmass_nse(1:ionmax))
      yebn = sum(yei(1:ionmax)*xmbn(1:ionmax))
      yebp = sum(yei(1:ionmax)*xmbp(1:ionmax))


! this formulation does not assume mion(i) = aion(i)*amu
! seitenzahl et al 2007

!     fac6(1:ionmax) = amu/mion(1:ionmax) * ((ye_want - 1.0d0)*zion(1:ionmax) + ye_want*nion(1:ionmax))
!     sum2   = sum(fac6(1:ionmax) * xmass_nse(1:ionmax)) 
!     sum2bn = sum(fac6(1:ionmax) * xmbn(1:ionmax)) 
!     sum2bp = sum(fac6(1:ionmax) * xmbp(1:ionmax)) 


! mass and charge conservation are the requirements
      f(1) = xmsum - 1.0d0
      f(2) = ye - ye_want

!      f(2) = sum2


! jacobian
      dfdy(1,1) = xmsumbn ; dfdy(1,2) = xmsumbp
      dfdy(2,1) = yebn    ; dfdy(2,2) = yebp

!      dfdy(2,1) = sum2bn  ; dfdy(2,2) = sum2bp

      return
      end





      subroutine xnewt_nse(ntrial,x,n,tolx,tolf,ntaken,check,nfev,func)
      include 'implno.dek'


! given an initial guess x(1:n) for the root of n equations, this routine
! finds the root by a globally convergent newtons method. the vector of
! functions to be zeroed, called fvec(1:n) in the routine below, is
! returned by the user supplied routine func. the output quantity check
! is false on nomal return and true if the routine has converged to a
! local minimum of the function xfminx_nse. if so, try restarting from a
! different initial guess.

! np is the maximum number of equations n
! ntrial is the maximum number of iterations to try
! ntaken is the number of iterations done
! tolf sets the convergence on function values
! tolmin sets the criterion for deciding wheather spurious convergence to
!        a false minimum of xfminx_nse has occured
! tolx is the convergence criteria on deltax
! stpmx is the scaled maximum step length allowed in the line searches
! nfev is the number of function evaluations


! declare the pass
      external :: func
      logical  :: check
      integer  :: ntrial,n,ntaken,nfev
      real*8   :: x(n),tolx,tolf


! common block communicates values from routine xfminx_nse
      integer, parameter :: np = 4
      integer  ::  nn
      real*8   :: fvec(np)
      common /newtnse/ fvec,nn

! locals
      integer :: i,its,j,indx(np)
      real*8  :: d,den,f,fold,stpmax,xsum,temp,test, &
                 fjac(np,np),g(np),p(np),xold(np),xfminx_nse,dum
      real*8, parameter :: tolmin = 1.0d-12, stpmx = 2.0d0



! initialize
      if (n .gt. np) stop 'n > np in routine xnewt'
      nn   = n ; f = xfminx_nse(x,func) ; nfev = 1 ; ntaken = 0

!  test for the initial guess being a root, using a more stringent tolf
      test = maxval(abs(fvec(1:n)))
      if (test .lt. 0.01*tolf) then
       check = .false.
       return
      end if


! get stpmax for the line search
      xsum = sum(x(1:n)*x(1:n))
      stpmax = stpmx * max(sqrt(xsum),dfloat(n))


! start of iteration loop; get the jacobian
      do its = 1, ntrial
       ntaken = its


! second order accurate numerical jacobian
!       call jac_nse(dum,x,fjac,n,n,np,np,func)
!       nfev = nfev + 2*n + 1


! analytic jacobian
       call nsejac(dum,x,fvec,fjac,n,np)
       nfev = nfev + 1


! compute grad f for the line searches
       do i=1,n
        g(i) = sum(fjac(1:n,i)*fvec(1:n))
       enddo


! store x, and f and form right hand sides
       xold(1:n) = x(1:n)  ; fold = f  ;  p(1:n) = -fvec(1:n)


! solve the linear systems

!       write(6,*)
!       write(6,112) x(1),x(2)
!       write(6,112) fvec(1),fvec(2)
!       write(6,112) fjac(1,1),fjac(1,2),fjac(2,1),fjac(2,2)
! 112   format(1x,1p2e14.6)


       call ludcmp(fjac,n,np,indx,d)
       call lubksb(fjac,n,np,indx,p)


! line search returns new x and f
! it also gets fvec at the new x when it calls xfminx_nse

       call lnsrch_nse(n,xold,fold,g,p,x,f,stpmax,check,nfev,func)

!       write(6,112) x(1),x(2)
!       write(6,112) f
!       write(6,*)

! test for convergence on function value
       test = maxval(abs(fvec(1:n)))
       if (test .lt. tolf) then
        check = .false.
        return
       end if

! check for zero gradiant, i.e. spurious convergence
       if (check) then
        den  = max(f, 0.5d0 * n)
        test = maxval( abs(g(1:n) * max(abs(x(:1:n)),1.0d0)/den) )
        if (test .lt. tolmin) then
         check = .true.
        else
         check = .false.
        end if
        return
       end if

! test for convergence on deltax
       test = maxval((abs(x(1:n)-xold(1:n)))/max(abs(x(1:n)),1.0d0))
       if (test .lt. tolx) return

!       write(6,*) its,test

! back for another iteration
      enddo
      check = .true.
      return
      end




      subroutine lnsrch_nse(n,xold,fold,g,p,x,f,stpmax,check,nfev,func)
      include 'implno.dek'

! given an n dimensional point xold(1:n), the value of the function fold
! and the gradient g(1:n) at the point, and a direction p(1:n), this routine
! finds a new point x(1:n) along the direction of p from xold where the
! function xfminx_nse has decreased "sufficiently". the new function value is
! returned in f. stpmax is an input quanity that limits the length of the
! steps so that the function is not evaluated in regions where it is
! undefined or subject to overflow. p is usually the newton direction. the
! output quantity check is false on normal exit, and true when x is too
! close to xold. in a minimization routine, this usually signals
! convergence and can be ignored. however, in a root finding routine, the
! calling routine should check wheather the convergence is spurious.


! declare the pass
      external :: func
      logical  :: check
      integer  :: n,nfev
      real*8   :: f,fold,stpmax,g(n),p(n),x(n),xold(n)


! locals
      integer :: i
      real*8  :: xfminx_nse,a,alam,alam2,alamin,b,disc,f2,rhs1, &
                 rhs2,slope,xsum,temp,test,tmplam

! alf ensures decreases in the function value, tolx is the convergence criterion on deltax
      real*8, parameter :: alf  = 1.0d-4, tolx = 3.0d-12, alam_start = 1.0d0


! initialize and scale if the attempted step is too big
      check = .false.
      xsum = sqrt(sum(p(1:n)*p(1:n)))

      if (xsum .gt. stpmax) then
       xsum = 1.0d0/xsum
       p(1:n) = p(1:n) * stpmax*xsum
      end if

      slope = sum(g(1:n)*p(1:n))
      if (slope .ge. 0.0) stop 'roundoff problem in lnsrch_nse'


! compute lambda_min
      test = maxval(abs(p(1:n))/max(abs(xold(1:n)),1.0d0))
      alamin = tolx/test


! always try a full newton step, start of iteration loop
      alam = alam_start
1     continue
      do i=1,n
       x(i) = xold(i) + alam*p(i)


! for the nse problem make sure the neutron and proton
! chemical potentials are less than or equal to zero
! hmm for low ye (0.2) and high density (1e14), the neutron
! chemical potential does indeed go positive. let's try
! letting it go wherever it wants

!       if (x(i) .gt. 0.0) x(i) = xold(i) + 1.0d-4*alam*p(i)
       if (x(i) .gt. 0.0) x(i) = 1.0d-10
      enddo

      f    = xfminx_nse(x,func)
      nfev = nfev + 1


! convergence on deltax, for root finding, the calling routine
! should verify the convergence
      if (alam .lt. alamin) then
       x(1:n) = xold(1:n)
       check = .true.
       return

! sufficient function decrease
      else if (f .le. fold + alf*alam*slope) then
       return

! backtrack
      else
       if (alam .eq. alam_start) then
        tmplam = -slope / (2.0d0 * (f-fold-slope))
       else
        rhs1 = f  - fold - alam*slope
        rhs2 = f2 - fold - alam2*slope
        a    = (rhs1/alam**2 - rhs2/alam2**2)/(alam-alam2)
        b    = (-alam2*rhs1/alam**2 + alam*rhs2/alam2**2) / (alam-alam2)
        if (a .eq. 0.0) then
         tmplam = -slope/(2.0d0 * b)
        else
         disc = b*b - 3.0d0 * a * slope
         if (disc .lt. 0.0) then
          tmplam = 0.5d0 * alam
         else if (b .le. 0.0) then
          tmplam = (-b + sqrt(disc)) / (3.0d0 * a)
         else
          tmplam = -slope/(b + sqrt(disc))
         end if
        end if
        if (tmplam .gt. 0.5d0*alam) tmplam = 0.5d0*alam
       end if
      end if

! store for the next trip through
      alam2 = alam
      f2    = f
      alam  = max(tmplam, 0.1d0*alam)
      goto 1
      end



      real*8 function xfminx_nse(x,func)
      include 'implno.dek'


! returns f = 0.5 f dot f at x. func is a user supplied routine of the
! functions to be root found.

! declare the pass
      external :: func
      real*8   :: x(*)

! locals
      integer :: i
      real*8  :: dum

! common block communicates values back to routine xnewt
      integer :: nn
      integer, parameter :: np = 4
      real*8 :: fvec(np)
      common /newtnse/ fvec,nn

      call func(dum,x,fvec)
      xfminx_nse = 0.5d0 * sum(fvec(1:nn)*fvec(1:nn))
      return
      end



      subroutine jac_nse(x,y,dfdy,mcol,nrow,mmax,nmax,derivs)
      include 'implno.dek'

! this routine computes a second order accurate jacobian matrix
! of the function contained in the routine derivs.
!
! input is the point x and the the vector y(nrow) at which to compute the
! jacobian dfdy(nrow,mcol).
!
! uses 2*nrow + 1 function evaluations


! declare the pass
      external :: derivs
      integer  :: mcol,nrow,mmax,nmax
      real*8   :: x,y(nmax),dfdy(mmax,nmax)


! locals
      integer :: i,j
      integer, parameter :: imax = 4
      real*8 :: fminus(imax),fplus(imax),temp,h,hinv
      real*8, parameter :: rel = 3.162278d-8, ax  = 1.0d-16


! check
       if (nrow .gt. imax) stop 'nrow > imax in jacobian2'


! for each row, get the right stepsize
      do j=1,nrow
       temp = y(j)
       h    = rel * max(abs(y(j)),ax)
       y(j) = temp + h
       h    = y(j) - temp
       call derivs(x,y,fplus)
       y(j) = temp

       temp = y(j)
       y(j) = temp - h
       h    = temp - y(j)
       call derivs(x,y,fminus)
       y(j) = temp

! compute the jth row of the jacobian
        hinv = 1.0d0/(2.0d0 * h)
        do i=1,mcol
         dfdy(i,j) = (fplus(i) - fminus(i)) * hinv
        enddo
       enddo


! restore the original state
      call derivs(x,y,fplus)
      return
      end
