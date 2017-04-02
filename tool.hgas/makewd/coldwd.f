      program coldwd
      include 'implno.dek'
      include 'const.dek'

c..
c..this program generats cold white dwarf models
c..

c..declare
      character*80     string
      integer          jj,nstep,ionmax
      parameter        (ionmax=2)
      double precision mass,xn,kpoly,alfa,surf,m53,m43,
     1                 temp,den,abar,zbar,dlo,dhi,dstp,presc,dpresdd,
     2                 ener,denerdd,denold,
     3                 xmass(ionmax),ymass(ionmax),
     4                 aion(ionmax),zion(ionmax), t_dyn, den_mean

      double precision con1,third
      parameter        (con1  = 4.0d0/3.0d0 * pi,
     1                  third = 1.0d0/3.0d0)


c..declare the integration technology
      external          derv,rkqc
      integer           xdim,ydim,kount,nok,nbad,iprint,i
      parameter         (xdim=600, ydim=2)
      double precision  xrk(xdim),yrk(ydim,xdim),bc(ydim),dydx,stptry,
     1                  stpmin,tol,odescal,start,sstop,dxsav, gamma


c..stop integrating when the density or pressure drop below these numbers
      double precision stopden,stop_pres
      common  /sbyiol/ stopden,stop_pres



c..communicate if general relativistic corrections are included,
c..and when to start a using a new central density

      double precision denc,ye
      logical          genrel,new_model
      common /grcm1/   denc,ye,genrel,new_model  




c..formats
 01   format(1x,1p8e12.4)
 02   format(1x,1p5e12.4,i5)
 03   format(a) 
 04   format(1x,a,t10,a,t22,a,t34,a,t46,a,t58,a,t70,a)
 05   format(1x,t5,a,t18,a,t30,a,t42,a,t54,a,t66,a,t78,a)




c..set a flag if general relativistic corrections are to be used
c      genrel = .true.
      genrel = .false.



c..set the composition, half carbon half oxygen
      xmass(1) = 0.5d0
      aion(1)  = 12.0
      zion(1)  = 6.0

      xmass(2) = 0.5d0
      aion(2)  = 16.0
      zion(2)  = 8.0



c..get abar and zbar      
      call azbar(xmass,aion,zion,ionmax,
     1           ymass,abar,zbar)

      ye = zbar/abar



c..get the central density and find the central pressure
c      write(6,*) 'give the central density =>'
c      read(5,*) denc
c      denc = 4.15e7 (1.15Mo)
      open(unit=917, file="input.dat", status="old")
c      denc = 3.17e7
      read(917,*) denc
      close(917)

      call fergas(denc,ye,presc,dpresdd,ener,denerdd)


     
c..get the output file name
C      write(6,*) 'give output file name =>'
C      read(5,03) string

      open(unit=3, file="poly.dat", status='unknown')



c..for comparing the central density to total mass 
c..nstep evenly spaced points in log space

      dlo   = 4.0d0
      dhi   = 11.0d0
      nstep = 200 
      dstp = (dhi - dlo)/(nstep -1)
      gamma = 1.4d0


c.sweep over the central density

!      write(3,05) 'density','massg',
!     1             'n=3/2 mass','n=3 mass','radius'



c      do jj = 1,nstep
c       denc = dlo + float(jj-1) * dstp
c       denc = 10.0d0**denc


c..cold ideal fermi gas
c..get the central pressure
       call fergas(denc,ye,presc,dpresdd,ener,denerdd)



c..integrate in radius (in cm) from start to sstop
      start   = 1.0e4
      stptry  = start
      stpmin  = 1.0d0
      sstop   = 30.0 * rsol


c..set the initial conditions: 
c..bc(1) is the gravitational mass 
c..bc(2) is the central pressure
c..these are expansions about zero radius

      bc(1)   = con1 * start**3 * denc
      bc(2)   = presc - 0.5d0 * con1 * g * start**2 * denc**2


c..control parameters
      new_model = .true.
      stopden   = max(1.0d0, 1.0e-6 * denc)


c..cold ideal fermi gas
c..get the stopping pressure
      call fergas(stopden,ye,stop_pres,dpresdd,ener,denerdd)


c..more control parameters
      tol     = 1.0e-6
      dxsav   = 0.0d0
      odescal = 1.0d0
      iprint  = 0



c..and integrate the odes in routine derv

      call fermint(start,stptry,stpmin,sstop,bc,
     1             tol,dxsav,xdim,
     2             xrk,yrk,xdim,ydim,xdim,ydim,
     3             nok,nbad,kount,odescal,iprint,
     4             derv,rkqc)
c..for dynamical timescale t_dyn = 1.0/sqrt(G*den_mean)
c..      double precision den_mean, t_dyn

c..write out the profile
c      write(3,04) '  i','radius','mass_grav','pressure',
c     1            'density'
      denold = denc
      write(3,*) kount
      write(3,01) gamma
      write(3,01) xrk(kount)
      do i=1,kount

c..cold ideal fermi gas
       den = denold
       call invert_fergas(den,ye,yrk(2,i),dpresdd,ener,denerdd)
       denold = den
       
       write(3,01) xrk(i),den,yrk(4,i),1.0e7,ener,yrk(1,i)/msol
      enddo
      write(3,01) yrk(1,kount)/msol
c..output for terminal confirmation (Sato corrections)
      write(6,*) "Mass of WD is..."
      write(6,01) yrk(1,kount)/msol
c..Dynamical timescale
      den_mean = (3.0d0*yrk(1,kount))/(4.0d0*pi*xrk(kount)**3)
      t_dyn = sqrt(g*den_mean)
      t_dyn = 1.0/t_dyn;
c      write(3,01) t_dyn




c..write a summary and the n=3/2 and n=3 polytrope answers

c..n=3/2 polytrope
c      xn    = 1.5d0
c      kpoly = 1.004e13 * ye**(5.0d0*third)
c      surf  = 2.71406d0
c      alfa  = sqrt((xn+1.0d0)*kpoly/(4.0d0*pi*g)*denc**(1.0d0/xn-1.0d0))
c      m53   = 4.0d0 * pi * alfa**3 * denc * surf/msol

c..n=3 polytrope
c      xn    = 3.0d0
c      kpoly = 1.2435e15 * ye**(4.0d0*third)
c      surf  = 2.01824d0
c      alfa  = sqrt((xn+1.0d0)*kpoly/(4.0d0*pi*g)*denc**(1.0d0/xn-1.0d0))
c      m43   = 4.0d0 * pi * alfa**3 * denc * surf/msol


c      write(3,02) xrk(kount), den, ener, temp

c      write(6,02) denc,yrk(1,kount)/msol,m53,m43,xrk(kount)/rsol,kount


c..back for another central density
c      enddo

      close(unit=3)
      stop 'normal termination'
      end






      subroutine derv(x,y,dydx)
      include 'implno.dek'
      include 'const.dek'

c..this routine sets up the continuity and hydrostatic equilibrium ode's. 
c..x is the radial coordinate, y(1) is the gravitational mass, 
c..y(2) is the pressure

c..declare the pass
      double precision x,y(1),dydx(1)


c..communicate if general relativistic corrections are included,
c..and when to start a using a new central density

      double precision denc,ye
      logical          genrel,new_model
      common /grcm1/   denc,ye,genrel,new_model  



c..local variables
      integer          i,niter
      double precision massg,den,pres,dpresdd,ener,denerdd,cor,
     1                 denold,con1,c2
      parameter        (con1   = 4.0d0 * pi,
     1                  c2     = clight*clight)


c..map the input vector
      massg = y(1)
      pres  = y(2)


c..set the density to start with for new models
      if (new_model) then
       denold    = denc
       new_model = .false.
      endif 


c..cold ideal fermi gas
      den = denold
      call invert_fergas(den,ye,pres,dpresdd,ener,denerdd) 
      denold = den 



c..here is d(massg)/dr
      dydx(1) = con1 * x**2 * den  


c..here is d(press)/dr
      if (genrel) then
       cor = (1.0d0 + pres/(den*c2)) * 
     1       (1.0d0 + (con1*pres*x**3)/(massg*c2)) /
     2       (1.0d0 - (2.0d0*g*massg)/(x*c2))
      else
       cor = 1.0d0
      end if

      dydx(2) = -g * massg/x**2 * den * cor 

      return
      end





      subroutine fergas(den,ye,pres,dpdd,ener,dedd)
      include 'implno.dek'
      include 'const.dek'

c..the eos for a cold fermi gas 
c..see, for example, cox & guili equation 24.157 or chandra chapter 10, eq. 19-20 

c..input is the density den and the mean charge to mean weight ratio ye.
c..output is the pressure (in erg/cm**3), the pressure derivative with density,
c..the energy (in erg/g), and the energy derivative with density. 

c..declare the pass
      double precision den,ye,pres,dpdd,ener,dedd

c..local variables
      double precision x,dxdd,dpdx,fac1,fac2,fac3,dfac1,dfac2,dfac3,
     1                 gac1,gac2,dgac1,dgac2,dedx,deni,x2,x3,
     2                 third,twothird,fourthird,eightthird,
     3                 lambda,lam3,xcon,pcon

      parameter        (third      = 1.0d0/3.0d0,
     1                  twothird   = 2.0d0 * third,
     2                  fourthird  = 4.0d0 * third,
     3                  eightthird = 8.0d0 * third)

      parameter        (lambda = hbar/(me*clight),
     1                  lam3   = lambda*lambda*lambda,
     2                  xcon   = 3.0d0 * pi*pi * lam3 * avo, 
     3                  pcon   = me*clight*clight/(lam3*8.0d0* pi*pi) )


      if (den .le. 0.0) stop 'bad pass in routine fergas'

      deni = 1.0d0/den 

      x    = (xcon*den*ye)**third
      dxdd = third * x * deni

      x2   = x * x
      x3   = x2 * x


c..the pressure in erg/cm**3
c..note: fac3 is a way of writing arc-sinh(x)
      fac1 = sqrt(1.0d0 + x2)
      fac2 = twothird * x2  - 1.0d0
      fac3 = log(x + fac1) 
      pres = pcon * (x*fac1*fac2 + fac3)


c..pressure derivative with density
      dfac1 = x/fac1
      dfac2 = fourthird*x
      dfac3 = (1.0d0 + dfac1)/(x + fac1)
      dpdx  = pcon * (fac1*fac2 + x*dfac1*fac2 + x*fac1*dfac2 + dfac3)
      dpdd  = dpdx * dxdd


c..the internal energy in erg/cm**3
      gac1  = eightthird * x3
      gac2  = fac1 - 1.0d0
      ener  = pcon*gac1*gac2 - pres

c..energy derivative with density
      dgac1 = 8.0d0 * x2
      dgac2 = dfac1
      dedx  = pcon*(dgac1*gac2 + gac1*dgac2) - dpdx
      dedd  = dedx * dxdd

c..convert the energies into erg/g
      ener = ener * deni
      dedd = (dedd - ener) * deni

      return
      end








      subroutine invert_fergas(den,ye,pres,dpdd,ener,dedd)
      include 'implno.dek'

c..given the pressure, ye, and a guess for the density,
c..find the density and dpdd


c..declare the pass
      double precision den,ye,pres,dpdd,ener,dedd


c..local variables
      integer          i,niter        
      double precision pres1,f,df,ratio,dennew,z,denold,
     1                 eostol,fpmin
      parameter        (eostol = 1.0d-6,
     1                  fpmin  = 1.0d-14)


c..save the initial guess
      denold = den

c..newton loop
      do i=1,100
       call fergas(den,ye,pres1,dpdd,ener,dedd)
       f  = pres1/pres - 1.0d0
       df = dpdd/pres
       if (df .eq. 0.0) goto 11
       ratio  = f/df
       dennew = den - ratio
       z      = abs((dennew - den)/den)
       den    = dennew
       niter  = i
       if (z .lt. eostol .or. abs(ratio) .le. fpmin) goto 20
      enddo


 11   write(6,*) 
      write(6,*) 'newton-raphson failed in routine derv'
      write(6,01) 'pres  =',pres
 01   format(1x,5(a,1pe16.8))
      write(6,01) 'z     =',z,'  eostol=',eostol,'  fpmin =',fpmin
      write(6,01) 'den   =',den,'  denold=',denold 
      write(6,01) 'f/df  =',ratio,' f   =',f,    ' df    =',df
      write(6,*) 
      stop 'could not find a density in routine invert_fergas'
 20   continue


c..call it one more time with the converged value of the density
      call fergas(den,ye,pres1,dpdd,ener,dedd)

      return
      end





      subroutine azbar(xmass,aion,zion,ionmax,
     1                 ymass,abar,zbar)
      include 'implno.dek'

c..this routine calculates composition variables for an eos routine

c..input:
c..mass fractions     = xmass(1:ionmax)
c..number of nucleons = aion(1:ionmax)
c..charge of nucleus  = zion(1:ionmax)
c..number of isotopes = ionmax
c..
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








      subroutine fermint(start,stptry,stpmin,stopp,bc,
     1                   eps,dxsav,kmax, 
     2                   xrk,yrk,xphys,yphys,xlogi,ylogi,
     3                   nok,nbad,kount,odescal,iprint,
     4                   derivs,steper)  
      include 'implno.dek'



c..generic ode driver, tuned a bit for cold white dwarfs


c..declare the pass
      external         derivs,steper
      integer          xphys,yphys,xlogi,ylogi,nok,nbad,
     1                 kmax,kount,iprint
      double precision start,stptry,stpmin,stopp,bc(yphys),eps,
     2                 dxsav,xrk(xphys),yrk(yphys,xphys),odescal 


c..local variables
      integer          nmax,stpmax,i,j,nstp
      parameter        (nmax = 20, stpmax=10000)  
      double precision yscal(nmax),y(nmax),dydx(nmax),  
     1                 x,xsav,h,hdid,hnext,zero,one,tiny
      parameter        (zero=0.0d0, one=1.0d0, tiny=1.0d-15)



c..stop integrating when the density or pressure drop below these numbers
      double precision stopden,stop_pres
      common  /sbyiol/ stopden,stop_pres



c..here are the format statements for printouts as we integrate
 100     format(1x,i4,1p12e10.2)



c..initialize   
      if (ylogi .gt. yphys) stop 'ylogi > yphys in routine fermint'
      if (yphys .gt. nmax)  stop 'yphys > nmax in routine fermint'
      x     = start   
      h     = sign(stptry,stopp-start) 
      nok   = 0 
      nbad  = 0
      kount = 0   


c..store the first step 
      do i=1,ylogi
       y(i) = bc(i)  
      enddo
      xsav = x - 2.0d0 * dxsav



c..take at most stpmax steps
      do nstp=1,stpmax
       call derivs(x,y,dydx)


c..scaling vector used to monitor accuracy  
       do i=1,ylogi

c..constant fractional accuracy        
c        yscal(i) = abs(y(i)) + tiny

c..const. frac. cept near zero  
c        yscal(i)=abs(y(i)) + abs(h*dydx(i)) + tiny  

c..step size dependent accuracy   
c..         yscal(i) = abs(odescal * h * dydx(i)) + tiny  

c..for stiffs
        yscal(i) = max(odescal,y(i))

c..strait scaling (decrease to get more time steps)
c        yscal(i) = odescal
       enddo



c..store intermediate results   
       if (kmax .gt. 0) then
        if ( (abs(dxsav) - abs(x-xsav)) .le. tiny) then 
         if ( kount .lt. (kmax-1) ) then  
          kount         = kount+1  
          xrk(kount)    = x   
          do i=1,ylogi 
           yrk(i,kount) = y(i)
          enddo
          if (iprint .eq. 1) then
           write(6,100) kount,xrk(kount),(yrk(j,kount), j=1,ylogi)
          end if
          xsav=x 
         end if
        end if  
       end if



c..if the step can overshoot the stop point or the dxsav increment then cut it
       if ((x+h-stopp)*(x+h-start) .gt. zero) h = stopp - x  
       if (dxsav.ne.zero .and. h.gt.(xsav-x+dxsav)) h = xsav + dxsav-x


c..do an integration step
       call steper(y,dydx,ylogi,x,h,eps,yscal,hdid,hnext,derivs)   
       if (hdid.eq.h) then
        nok = nok+1   
       else 
        nbad = nbad+1 
       end if




c..this is the normal exit point, save the final step   
       if (nstp .eq. stpmax .or. (x-stopp)*(stopp-start) .ge. zero
     1     .or. y(2) .le. stop_pres) then
        do i=1,ylogi  
         bc(i) = y(i) 
        enddo
        if (kmax.ne.0) then   
         kount         = kount+1  
         xrk(kount)    = x   
         do i=1,ylogi 
          yrk(i,kount) = y(i) 
         enddo
         if (iprint .eq. 1) then
           write(6,100) kount,xrk(kount),(yrk(j,kount), j=1,ylogi)
         end if
        end if
        return  
       end if


c..set the step size for the next iteration; stay above stpmin
       h=hnext
       if (abs(hnext).lt.stpmin) return
c       if (abs(hnext).lt.stpmin) stop 'hnext < stpmin in fermint'


c..back for another iteration or death
      enddo
      write(6,*) '> than stpmax steps required in fermint' 
      return
      end










c..explicit: 
c..routine rk4 is a fourth order runge-kutta stepper 
c..routine rkqc is a step doubling rk4 integrator; plug into odeint 




       subroutine rk4(y,dydx,n,x,h,yout,derivs)  
       include 'implno.dek' 
c..  
c..given values for the variables y(1:n) and their derivatives dydx(1:n) known 
c..at x, use the fourth order runge-kutta method to advance the solution over 
c..an interval h and return the incremented variables in yout(1:n) (which need 
c..not be a distinct array from y). one supplies the routine derivs which  
c..evaluates the right hand side of the ode's.    
c..  
c..declare   
       external          derivs  
       integer           n,nmax,i    
       parameter         (nmax = 2000)   
       double precision  x,h,y(n),dydx(n),yout(n), 
     1                   yt(nmax),dyt(nmax),dym(nmax), 
     2                   hh,h6,xh    

c..initialize the step sizes and weightings  
       hh = h*0.5d0 
       h6 = h/6.0d0  
       xh = x + hh   

c..the first step    
       do i=1,n   
        yt(i) = y(i) + hh*dydx(i)   
       enddo

c..the second step   
       call derivs(xh,yt,dyt)    
       do i=1,n   
        yt(i) = y(i) + hh*dyt(i)    
       enddo

c..the third step    
       call derivs(xh,yt,dym)    
       do i=1,n   
        yt(i)  = y(i) + h*dym(i) 
        dym(i) = dyt(i) + dym(i)    
       enddo

c..the fourth step and accumulate the increments with the proper weights 
       call derivs(x+h,yt,dyt)   
       do i=1,n   
        yout(i) = y(i) + h6*(dydx(i) +dyt(i) + 2.0d0*dym(i))  
       enddo
       return    
       end   





      subroutine rkqc(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs)   
      include 'implno.dek' 

c..fifth order, step doubling, runge-kutta ode integrator with monitering of  
c..local truncation errors. input are the vector y of length n, which has a  
c..known the derivative dydx at the point x, the step size to be attempted  
c..htry, the required accuracy eps, and the vector yscal against which the  
c..error is to be scaled. on output, y and x are replaced by their new values,
c..hdid is the step size that was actually accomplished, and hnext is the  
c..estimated next step size. derivs is a user supplied routine that computes  
c..the right hand side of the first order system of ode's. plug into odeint. 

c..declare   
      external         derivs    
      integer          n,nmax,i  
      parameter        (nmax = 2000)   
      double precision x,htry,eps,hdid,hnext,y(n),dydx(n),yscal(n),  
     1                 ytemp(nmax),ysav(nmax),dysav(nmax),fcor,safety, 
     2                 errcon,pgrow,pshrnk,xsav,h,hh,errmax  
      parameter        (fcor=1.0d0/15.0d0, pgrow = -0.2d0,  
     1                  pshrnk = -0.25d0,  safety=0.9d0,  errcon=6.0e-4) 

c..note errcon = (4/safety)**(1/pgrow)   
c..nmax is the maximum number of differential equations 


c..save the initial values  
      h      = htry 
      xsav   =  x 
      do i=1,n    
       ysav(i)  = y(i) 
       dysav(i) = dydx(i) 
      enddo

c..take two half steps   
1     hh = 0.5d0*h   
      call rk4(ysav,dysav,n,xsav,hh,ytemp,derivs)    
      x  = xsav + hh  
      call derivs(x,ytemp,dydx)  
      call rk4(ytemp,dydx,n,x,hh,y,derivs)   
      x  = xsav + h   
      if (x .eq. xsav) stop 'stepsize not significant in rkqc'  

c..now take the large step   
      call rk4(ysav,dysav,n,xsav,h,ytemp,derivs) 

c..ytemp is the error estimate   
      errmax = 0.0 
      do i=1,n    
       ytemp(i) = y(i) - ytemp(i)   
       errmax   = max(errmax,abs(ytemp(i)/yscal(i)))    
      enddo
      errmax     = errmax/eps  

c..truncation error too big, reduce the step size and try again  
      if (errmax .gt. 1.0) then 
       h = safety * h * (errmax**pshrnk)  
       go to  1   

c..truncation within limits, compute the size of the next step   
      else   
       hdid = h   
       if (errmax.gt.errcon) then 
        hnext = safety * h * (errmax**pgrow) 
       else 
        hnext = 4.0d0 * h 
       end if    
      end if  

c..mop up the fifth order truncation error   
      do i=1,n    
       y(i) = y(i) + ytemp(i)*fcor  
      enddo
      return 
      end 







