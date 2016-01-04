!      subroutine nse(tt,dd,yye,newguess,xmass_out,iprint)
      subroutine nse0(tt,dd,yye,newguess,xmass_out,iprint)
      include 'implno_nse.dek'
      include 'const_nse.dek'
      include 'network_nse.dek'


c..given the temperature tt, density dd, and electron mole number ye
c..this routine puts a chosen reaction network into its nse distribution.

c..input:
c..tt = temperature
c..dd = density
c..ye = electron mol number
c..newguess = 1 = a new initial guess is made  
c..         = 0 = values from the previous call are used as the initial guess  
c..iprint = print flag   

c..output:
c..xmass_out  = nse mass fraction 


c..declare the pass
      integer          newguess,iprint
      double precision tt,dd,yye,xmass_out(1)


c..communicate
      double precision temp,den,ye_want,beta
c      common /nsec1/   temp,den,ye_want,beta
      common /nsec1_nse/   temp,den,ye_want,beta


c..locals
!      external         nsefunc
      external         nsefunc0
      logical          check
      integer          ntrial,nfev,ntaken,n,i
      parameter        (ntrial = 200, n = 2)
      double precision x(n),amass,fac1,fac2,tolf,tolx,twopi,
     1                 dum,resid(n)
      parameter        (tolf = 1.0d-8, tolx = 1.0d-14, twopi=2.0d0*pi)



c..fill the common block
      temp    = tt
      den     = dd
      ye_want = yye
      beta    = 1.0d0/(kerg * temp)



c..set the partition functions
c..these are generally temperature dependent 
c..here i'll just take the ground state

      do i=1,ionmax
       wpart(i) = 1.0d0
      enddo
      wpart(ineut) = 2.0d0
      wpart(iprot) = 2.0d0




c..here is an initial guess for the neutron and proton chemical potentials,
c..(x1) and (x2) respectively. obtained by setting xmass(ini56) = 1, 
c..setting mup = mun, and inverting the saha equation. 
c..all nse networks should at least have ni56 and this appears to be a 
c..decent guess for all temp, rho, ye combinations.
      
      if (newguess .eq. 1) then
       newguess = 0
       i      = ini56
       amass  = aion(i) * amu
       fac1   = aion(i)/(avo * den) * wpart(i)
       fac2   = (twopi/(beta*h) * amass/h )**1.5d0 
       x(1)   = -(log(fac1*fac2)/beta + bion(i)*ev2erg*1.0d6)/aion(i)
       x(2)   = x(1) 
      end if



c..root find on mass and charge conservation for
c..the chemical potentials of protons and neutrons

!      call xnewt_nse(ntrial,x,n,tolx,tolf,ntaken,check,nfev,nsefunc)
!      call xnewt_nse(ntrial,x,n,tolx,tolf,ntaken,check,nfev,nsefunc0)
      call xnewt_nse0(ntrial,x,n,tolx,tolf,ntaken,check,nfev,nsefunc0)


c..be sure we converged
      if (check .or. ntaken .eq. ntrial) then
       write(6,*)
       write(6,*) 'check convergence of root finder'
       write(6,*)
      end if


c..some optional diagnostics
      if (iprint .eq. 1) then
       write(6,*) 
       write(6,110) 'iterations taken             =',ntaken
       write(6,110) 'function evals               =',nfev
       write(6,111) 'roots                        =',x(1),x(2)
!       call nsefunc(dum,x,resid)
       call nsefunc0(dum,x,resid)
       write(6,111) 'mass conservation   residual =',resid(1)
       write(6,111) 'charge conservation residual =',resid(2)

 110   format(1x,a,i4)
 111   format(1x,a,1p2e14.6)
      end if 


      if (check .or. ntaken .eq. ntrial) stop



c..fill the output array using the converged values
!      call nsefunc(dum,x,resid)
      call nsefunc0(dum,x,resid)
      do i=1,ionmax
       xmass_out(i) = xmass(i)
      enddo

      return
      end







!      subroutine nsefunc(x,y,f)
      subroutine nsefunc0(x,y,f)
      include 'implno_nse.dek'
      include 'const_nse.dek'
      include 'network_nse.dek'


c..this routine returns the root find functions.
c..input is the point x and y a vector of the unknowns. 
c..output is the vector of root find functions f, which should be the
c..zero vector upon convergence.

c..y(1) is input as the neutron chemical potential
c..y(2) is input as the proton chemical potential


c..declare the pass
      double precision x,y(*),f(*)


c..locals
      integer          i,indx(abignet),j
      double precision ye,mu,amass,fac1,fac2,fac3,sum,twopi
      parameter        (twopi = 2.0d0 * pi)


c..communicate
      double precision temp,den,ye_want,beta
c      common /nsec1/   temp,den,ye_want,beta
      common /nsec1_nse/   temp,den,ye_want,beta



c..chemical potential and mass fraction of each isotope
c..hartmann et al, 1985 apj 297 837, eq 2
c..take the mass of each isotope to be amu * aion, otherwise a mass formula

      do i=1,ionmax
       mu       = (aion(i) - zion(i))*y(1) + zion(i)*y(2)
       amass    = aion(i) * amu
       fac1     = aion(i)/(avo * den) * wpart(i)
       fac2     = ( twopi/(beta*h) * amass/h )**1.5d0
       fac3     = exp( beta * (mu + bion(i)*ev2erg*1.0d6) )
       xmass(i) = fac1 * fac2 * fac3
      enddo


c..sum the mass fractions in ascending order to minimize roundoff
!      call indexx(ionmax,xmass,indx)
      call indexx0(ionmax,xmass,indx)
      sum   = 0.0d0
      do i=1,ionmax
       sum   = sum + xmass(indx(i))
      enddo

c..sum the mass fractions to form ye
      ye = 0.0d0
      do i=1,ionmax
       j = indx(i)
       ye = ye + zion(j)/aion(j) * xmass(j)
      enddo

c..mass and charge conservation are the requirements
      f(1) = sum - 1.0d0
      f(2) = ye - ye_want

      return
      end








!      subroutine nsejac(x,y,f,dfdy,n,np)
      subroutine nsejac0(x,y,f,dfdy,n,np)
      include 'implno_nse.dek'
      include 'const_nse.dek'
      include 'network_nse.dek'

c..this routine returns the functions and the jacobian to do the root find on
c..input is x, and y(n) a vector of the unknowns. output is f(n) 
c..and its jacobian dfdy(np,np).

c..y(1) is the neutron chemical potential
c..y(2) is the proton chemical potential


c..declare the pass
      integer          n,np
      double precision x,y(*),f(*),dfdy(np,np)


c..locals
      integer          indx(ionmax),i,j
      double precision mu,mubn,mubp,amass,fac1,fac2,fac3,fac4,fac5,
     2                 xmbn(ionmax),xmbp(ionmax),sum,sumbn,sumbp,
     3                 ye,yebn,yebp,twopi
      parameter        (twopi = 2.0d0 * pi)


c..communicate
      double precision temp,den,ye_want,beta
c      common /nsec1/   temp,den,ye_want,beta
      common /nsec1_nse/   temp,den,ye_want,beta


c..chemical potential and mass fraction of each isotope
c..hartmann et al, 1985 apj 297 837, eq 2
c..take the mass of each isotope to be amu * aion, otherwise a mass formula

      do i=1,ionmax
       mu       = (aion(i) - zion(i)) * y(1) + zion(i) * y(2)
       mubn     = aion(i) - zion(i)
       mubp     = zion(i)

       amass    = aion(i) * amu

       fac1     = aion(i)/(avo * den) * wpart(i)
       fac2     = ( twopi/(beta*h) * amass/h )**1.5d0
       fac3     = exp( beta * (mu + bion(i) * ev2erg * 1.0d6) )
       fac4     = fac1 * fac2 * fac3 

       xmass(i) = fac4
       xmbn(i)  = fac4 * beta * mubn   
       xmbp(i)  = fac4 * beta * mubp   

      enddo


c..sum the mass fractions in ascending order to minimize roundoff
!      call indexx(ionmax,xmass,indx)
      call indexx0(ionmax,xmass,indx)
      sum   = 0.0d0
      sumbn = 0.0d0
      sumbp = 0.0d0
      do i=1,ionmax
       j     = indx(i)
       sum   = sum   + xmass(j)
       sumbn = sumbn + xmbn(j)
       sumbp = sumbp + xmbp(j)
      enddo


c..sum the mass fractions to form ye
      ye   = 0.0d0
      yebn = 0.0d0
      yebp = 0.0d0
      do i=1,ionmax
       j    = indx(i)
       fac5 = zion(j)/aion(j)
       ye   = ye   + fac5 * xmass(j)
       yebn = yebn + fac5 * xmbn(j)
       yebp = yebp + fac5 * xmbp(j)
      enddo


c..mass and charge conservation are the requirements
      f(1) = sum - 1.0d0
      f(2) = ye - ye_want


c..jacobian 
      dfdy(1,1) = sumbn
      dfdy(1,2) = sumbp
      dfdy(2,1) = yebn
      dfdy(2,2) = yebp

      return
      end







!      subroutine xnewt_nse(ntrial,x,n,tolx,tolf,ntaken,check,nfev,func) 
      subroutine xnewt_nse0(ntrial,x,n,tolx,tolf,ntaken,check,nfev,func) 
      include 'implno_nse.dek' 


c..given an initial guess x(1:n) for the root of n equations, this routine
c..finds the root by a globally convergent newtons method. the vector of 
c..functions to be zeroed, called fvec(1:n) in the routine below, is 
c..returned by the user supplied routine func. the output quantity check 
c..is false on nomal return and true if the routine has converged to a 
c..local minimum of the function xfminx_nse. if so, try restarting from a 
c..different initial guess. 

c..np is the maximum number of equations n
c..ntrial is the maximum number of iterations to try
c..ntaken is the number of iterations done
c..tolf sets the convergence on function values 
c..tolmin sets the criterion for deciding wheather spurious convergence to 
c..       a false minimum of xfminx_nse has occured
c..tolx is the convergence criteria on deltax 
c..stpmx is the scaled maximum step length allowed in the line searches 
c..nfev is the number of function evaluations


c..declare the pass
      external         func 
      logical          check 
      integer          ntrial,n,ntaken,nfev 
      double precision x(n),tolx,tolf


c..common block communicates values from routine xfminx_nse
      integer          nn,np 
      parameter        (np = 4) 
      double precision fvec(np)
c      common /newtnse/ fvec,nn 
      common /newtnse_nse/ fvec,nn 

c..locals
      integer          i,its,j,indx(np)
      double precision tolmin,stpmx,d,den,f,fold,stpmax,sum,temp,test,
!     1                 fjac(np,np),g(np),p(np),xold(np),xfminx_nse,dum 
     1                 fjac(np,np),g(np),p(np),xold(np),xfminx_nse0,dum 

      parameter        (tolmin = 1.0d-12, 
     1                  stpmx = 2.0d0) 



c..initialize 
      if (n .gt. np) stop 'n > np in routine xnewt' 
      nn     = n 
!      f      = xfminx_nse(x,func)
      f      = xfminx_nse0(x,func)
      nfev   = 1
      ntaken = 0


c.. test for the initial guess being a root, using a more stringent tolf
      test = 0.0d0 
      do i=1,n 
       if (abs(fvec(i)) .gt. test) test = abs(fvec(i)) 
      enddo
      if (test .lt. 0.01*tolf) then 
       check = .false. 
       return 
      end if 


c..get stpmax for the line search
      sum = 0.0d0 
      do i=1,n 
       sum = sum + x(i)*x(i)
      enddo
      stpmax = stpmx * max(sqrt(sum),dfloat(n)) 


c..start of iteration loop; get the jacobian 
      do its = 1, ntrial
       ntaken = its


c..second order accurate numerical jacobian
c       call jac_nse(dum,x,fjac,n,n,np,np,func) 
c       nfev = nfev + 2*n + 1


c..analytic jacobian
!       call nsejac(dum,x,fvec,fjac,n,np) 
       call nsejac0(dum,x,fvec,fjac,n,np) 
       nfev = nfev + 1


c..compute grad f for the line searches 
       do i=1,n 
        sum = 0.0d0 
        do j=1,n 
         sum = sum + fjac(j,i)*fvec(j) 
        enddo
        g(i) = sum 
       enddo


c..store x, and f and form right hand sides 
       do i=1,n 
        xold(i) = x(i) 
       enddo
       fold = f 
       do i=1,n 
        p(i) = -fvec(i) 
       enddo


c..solve the linear systems 
!       call ludcmp(fjac,n,np,indx,d) 
       call ludcmp0(fjac,n,np,indx,d) 
!       call lubksb(fjac,n,np,indx,p) 
       call lubksb0(fjac,n,np,indx,p) 


c..line search returns new x and f 
c..it also gets fvec at the new x when it calls xfminx_nse 
!       call lnsrch_nse(n,xold,fold,g,p,x,f,stpmax,check,nfev,func) 
       call lnsrch_nse0(n,xold,fold,g,p,x,f,stpmax,check,nfev,func) 


c..test for convergence on function value 
       test = 0.0d0 
       do i=1,n 
        if (abs(fvec(i)) .gt. test) test = abs(fvec(i)) 
       enddo
c       write(6,*) its,test,tolf
       if (test .lt. tolf) then 
        check = .false. 
c        write(6,*) 'returning on converged function value',test,tolf
        return 
       end if 

c..check for zero gradiant, i.e. spurious convergence 
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

c..test for convergence on deltax 
       test = 0.0d0 
       do i=1,n 
        temp = (abs(x(i)-xold(i)))/max(abs(x(i)),1.0d0) 
        if (temp .gt. test) test = temp 
       enddo
       if (test .lt. tolx) return 

c       write(6,*) its,test,tolx

c..back for another iteration
      enddo
      check = .true.
      return
      end 





!      subroutine lnsrch_nse(n,xold,fold,g,p,x,f,stpmax,check,nfev,func) 
      subroutine lnsrch_nse0(n,xold,fold,g,p,x,f,stpmax,check,nfev,func) 
      include 'implno_nse.dek' 

c..given an n dimensional point xold(1:n), the value of the function fold
c..and the gradient g(1:n) at the point, and a direction p(1:n), this routine 
c..finds a new point x(1:n) along the direction of p from xold where the 
c..function xfminx_nse has decreased "sufficiently". the new function value is 
c..returned in f. stpmax is an input quanity that limits the length of the 
c..steps so that the function is not evaluated in regions where it is 
c..undefined or subject to overflow. p is usually the newton direction. the 
c..output quantity check is false on normal exit, and true when x is too 
c..close to xold. in a minimization routine, this usually signals  
c..convergence and can be ignored. however, in a root finding routine, the  
c..calling routine should check wheather the convergence is spurious. 


c..declare the pass
      external         func 
      logical          check 
      integer          n,nfev 
      double precision f,fold,stpmax,g(n),p(n),x(n),xold(n)
       

c..locals
      integer          i
!      double precision xfminx_nse,a,alam,alam2,alamin,b,disc,f2,rhs1,
      double precision xfminx_nse0,a,alam,alam2,alamin,b,disc,f2,rhs1,
     1                 rhs2,slope,sum,temp,test,tmplam,
     2                 alf,tolx 
      parameter        (alf  = 1.0d-4, 
     1                  tolx = 3.0d-15) 


c..alf ensures sufficient decrease in the function value, tolx is the  
c..convergence criterion on deltax 


c..initialize and scale if the attempted step is too big
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
      if (slope .ge. 0.0) stop 'roundoff problem in lnsrch_nse'


c..compute lambda_min 
      test = 0.0d0 
      do i=1,n 
       temp = abs(p(i))/max(abs(xold(i)),1.0d0) 
       if (temp .gt. test) test = temp 
      enddo
      alamin = tolx/test 


c..always try a full newton step, start of iteration loop 
      alam = 1.0d0 
1     continue 
      do i=1,n 
       x(i) = xold(i) + alam*p(i) 
      enddo


!      f    = xfminx_nse(x,func)
      f    = xfminx_nse0(x,func)
      nfev = nfev + 1 



c..convergence on deltax, for root finding, the calling routine
c..should verify the convergence
      if (alam .lt. alamin) then 
       do i=1,n 
        x(i) = xold(i) 
       enddo
       check = .true. 
       return 

c..sufficient function decrease 
      else if (f .le. fold + alf*alam*slope) then 
       return 

c..backtrack 
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

c..store for the next trip through
      alam2 = alam 
      f2    = f 
      alam  = max(tmplam, 0.1d0*alam) 
      goto 1 
      end 






!      double precision function xfminx_nse(x,func) 
      double precision function xfminx_nse0(x,func) 
      include 'implno_nse.dek' 


c..returns f = 0.5 f dot f at x. func is a user supplied routine of the  
c..functions to be root found. 

c..declare the pass
      external         func 
      double precision x(1)


c..locals
      integer          i 
      double precision sum,dum 


c..common block communicates values back to routine xnewt
      integer          nn,np 
      parameter        (np = 4) 
      double precision fvec(np)
c      common /newtnse/ fvec,nn 
      common /newtnse_nse/ fvec,nn 



      call func(dum,x,fvec) 

      sum = 0.0d0 
      do i=1,nn 
       sum = sum + fvec(i)*fvec(i) 
      enddo

!      xfminx_nse = 0.5d0 * sum 
      xfminx_nse0 = 0.5d0 * sum 
      return 
      end 






!      subroutine jac_nse(x,y,dfdy,mcol,nrow,mmax,nmax,derivs)
      subroutine jac_nse0(x,y,dfdy,mcol,nrow,mmax,nmax,derivs)
      include 'implno_nse.dek'

c..this routine computes a second order accurate jacobian matrix 
c..of the function contained in the routine derivs.
c..
c..input is the point x and the the vector y(nrow) at which to compute the 
c..jacobian dfdy(mcol,nrow). 
c..
c..uses 2*nrow + 1 function evaluations 


c..declare the pass
      external         derivs
      integer          mcol,nrow,mmax,nmax
      double precision x,y(nmax),dfdy(mmax,nmax)
       

c..locals
      integer          i,j,imax
      parameter        (imax = 4)
      double precision fminus(imax),fplus(imax),rel,ax,temp,h,hinv
      parameter        (rel = 3.162278d-8, 
     1                  ax  = 1.0d-16)


c..check
       if (nrow .gt. imax) stop 'nrow > imax in jacobian2'


c..for each row, get the right stepsize
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

c..compute the jth row of the jacobian
        hinv = 1.0d0/(2.0d0 * h)
        do i=1,mcol
         dfdy(i,j) = (fplus(i) - fminus(i)) * hinv
        enddo
       enddo

c..restore the original state
      call derivs(x,y,fplus)
      return
      end








c..lu decomposition: 
c..routine ludcmp does a pivoting lower-upper decomposition  
c..routine lubksb does the backsubstitution from ludcmp 



!      subroutine ludcmp(a,n,np,indx,d) 
      subroutine ludcmp0(a,n,np,indx,d) 
      implicit none
      save

c..given the matrix a(n,n), with physical dimensions a(np,ap) this routine
c..replaces a by the lu decompostion of a row-wise permutation of itself. 
c..input are a,n,np. output is a, indx which records the row 
c..permutations effected by the partial pivoting, and d which is 1 if 
c..the number of interchanges is even, -1 if odd. 
c..use routine lubksb to solve a system of linear equations.
c.. 
c..nmax is the largest expected value of n 

c..declare 
      integer          n,np,indx(np),nmax,i,j,k,imax 
      parameter        (nmax=500) 
      double precision a(np,np),d,tiny,vv(nmax),aamax,sum,dum 
      parameter        (tiny=1.0d-20) 


c..vv stores the implicit scaling of each row 
c..loop over the rows to get the scaling information 
      d = 1.0d0 
      do i=1,n 
       aamax = 0.0d0 
       do j=1,n 
        if (abs(a(i,j)) .gt. aamax) aamax = abs(a(i,j)) 
       enddo
       if (aamax .eq. 0.0) stop 'singular matrix in ludcmp' 
       vv(i) = 1.0d0/aamax 
      enddo

c..for each column apply crouts method; see equation 2.3.12 
      do j=1,n 
       do i=1,j-1 
        sum = a(i,j) 
        do k=1,i-1 
         sum = sum - a(i,k)*a(k,j) 
        enddo
        a(i,j) = sum 
       enddo

c..find the largest pivot element 
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

c..if we need to interchange rows 
       if (j .ne. imax) then 
        do k=1,n 
         dum       = a(imax,k) 
         a(imax,k) = a(j,k) 
         a(j,k)    = dum 
        enddo
        d          = -d 
        vv(imax)   = vv(j) 
       end if 

c..divide by the pivot element 
       indx(j) = imax 
       if (a(j,j) .eq. 0.0) a(j,j) = tiny 
       if (j .ne. n) then 
        dum = 1.0d0/a(j,j) 
        do i=j+1,n 
         a(i,j) = a(i,j)*dum 
        enddo
       end if 

c..and go back for another column of crouts method
      enddo
      return 
      end 





!      subroutine lubksb(a,n,np,indx,b) 
      subroutine lubksb0(a,n,np,indx,b) 
      implicit none
      save

c..solves a set of n linear equations ax=b. a is input in its lu decomposition 
c..form, determined by the routine above ludcmp. indx is input as the 
c..permutation vector also returned by ludcmp. b is input as the right hand 
c..side vector and returns with the solution vector x. 
c..a,n ans np are not modified by this routine and thus can be left in place 
c..for successive calls (i.e matrix inversion) 

c..declare 
      integer           n,np,indx(np),i,ii,j,ll 
      double precision  a(np,np),b(np),sum 

c..when ii is > 0, ii becomes the index of the first nonzero element of b 
c..this is forward substitution of equation 2.3.6, and unscamble in place
      ii = 0 
      do i=1,n 
       ll = indx(i) 
       sum = b(ll) 
       b(ll) = b(i) 
       if (ii .ne. 0) then 
        do j=ii,i-1 
         sum = sum - a(i,j) * b(j) 
        enddo

c..nonzero element was found, so dos the sums in the loop above 
       else if (sum .ne. 0.0) then 
        ii  = i 
       end if 
       b(i) = sum 
      enddo

c..back substitution equation 2.3.7 
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







!      subroutine indexx(n,arr,indx) 
      subroutine indexx0(n,arr,indx) 
      include 'implno_nse.dek' 
c.. 
c..indexes an array arr(1:n). that is it outputs the array indx(1:n) such 
c..that arr(indx(j)) is in ascending order for j=1...n. the input quantities 
c..are not changed. 
c.. 
c..declare 
      integer          n,indx(n),m,nstack 
      parameter        (m=7, nstack = 50) 
      integer          i,indxt,ir,itemp,j,jstack,k,l,istack(nstack) 
      double precision arr(n),a 
c.. 
c..initialize 
      do 11 j=1,n 
       indx(j) = j 
11    continue 
      jstack = 0 
      l      = 1 
      ir     = n 
c.. 
c..insertion sort when subbarray small enough 
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
c.. 
c..pop stack and begin a new round of partitioning 
       if (jstack .eq. 0) return 
       ir     = istack(jstack) 
       l      = istack(jstack-1) 
       jstack = jstack - 2 
c.. 
c..choose median of left, center and right elements as partitioning element 
c..also rearrange so that a(l+1) < a(l) < a(ir) 
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
 
c.. 
c..initialize pointers for partitioning 
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
c.. 
5      indx(l+1) = indx(j) 
       indx(j)   = indxt 
       jstack    = jstack + 2 
c.. 
c..push pointers to larger subarray on stack 
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







      subroutine init_network0
      include 'implno_nse.dek'
      include 'network_nse.dek'

c..this routine initializes stuff for a network

c..declare
      integer          i


c..for easy zeroing of the isotope pointers
      integer          isotp(nisotp) 
      equivalence      (isotp(1),ih1)


c..zero all the isotope pointers
      do i=1,nisotp
       isotp(i)   = 0
      enddo


c..set the size of the network and the number of rates
ccc      ionmax  = 47
      ionmax  = 13


c..set the id numbers of the elements
ccc      ihe3  = 1
ccc      ihe4  = 2
ccc      ic12  = 3
ccc      ic13  = 4
ccc      in13  = 5
ccc      in14  = 6
ccc      io14  = 7
ccc      io15  = 8
ccc      io16  = 9
ccc      io17  = 10
ccc      io18  = 11
ccc      if17  = 12
ccc      if18  = 13
ccc      if19  = 14
ccc      ine18 = 15
ccc      ine19 = 16
ccc      ine20 = 17
ccc      img22 = 18
ccc      img24 = 19
ccc      ial27 = 20
ccc      isi28 = 21
ccc      ip31  = 22
ccc      is30  = 23
ccc      is32  = 24
ccc      icl35 = 25
ccc      iar36 = 26
ccc      ik39  = 27
ccc      ica40 = 28
ccc      iti44 = 29
ccc      icr48 = 30
ccc      icr49 = 31
ccc      icr50 = 32
ccc      icr51 = 33
ccc      icr52 = 34
ccc      icr53 = 35
ccc      icr54 = 36
ccc      ife52 = 37
ccc      ife54 = 38
ccc      ife55 = 39
ccc      ife56 = 40
ccc      ife57 = 41
ccc      ife58 = 42
ccc      ico55 = 43
ccc      ini56 = 44
ccc      ini58 = 45  
ccc      ineut = 46
ccc      iprot = 47
      ihe4  = 1
      ic12  = 2
      io16  = 3
      ine20 = 4
      img24 = 5
      isi28 = 6
      is32  = 7
      iar36 = 8
      ica40 = 9
      iti44 = 10
      icr48 = 11
      ife52 = 12
      ini56 = 13


c..set the names of the elements
ccc      ionam(ihe3)  = 'he3 '
ccc      ionam(ihe4)  = 'he4 '
ccc      ionam(ic12)  = 'c12 '
ccc      ionam(ic13)  = 'c13 '
ccc      ionam(in13)  = 'n13 '
ccc      ionam(in14)  = 'n14 '
ccc      ionam(io14)  = 'o14 '
ccc      ionam(io15)  = 'o15 '
ccc      ionam(io16)  = 'o16 '
ccc      ionam(io17)  = 'o17 '
ccc      ionam(io18)  = 'o18 '
ccc      ionam(if17)  = 'f17 '
ccc      ionam(if18)  = 'f18 '
ccc      ionam(if19)  = 'f19 '
ccc      ionam(ine18) = 'ne18'
ccc      ionam(ine19) = 'ne19'
ccc      ionam(ine20) = 'ne20'
ccc      ionam(img22) = 'mg22'
ccc      ionam(img24) = 'mg24'
ccc      ionam(ial27) = 'al27'
ccc      ionam(isi28) = 'si28'
ccc      ionam(ip31)  = 'p31 '
ccc      ionam(is30)  = 's30 '
ccc      ionam(is32)  = 's32 '
ccc      ionam(icl35) = 'cl35'
ccc      ionam(iar36) = 'ar36'
ccc      ionam(ik39)  = 'k39 '
ccc      ionam(ica40) = 'ca40'
ccc      ionam(iti44) = 'ti44'
ccc      ionam(icr48) = 'cr48'
ccc      ionam(icr49) = 'cr49'
ccc      ionam(icr50) = 'cr50'
ccc      ionam(icr51) = 'cr51'
ccc      ionam(icr52) = 'cr52'
ccc      ionam(icr53) = 'cr53'
ccc      ionam(icr54) = 'cr54'
ccc      ionam(ife52) = 'fe52'
ccc      ionam(ife54) = 'fe54'
ccc      ionam(ife55) = 'fe55'
ccc      ionam(ife56) = 'fe56'
ccc      ionam(ife57) = 'fe57'
ccc      ionam(ife58) = 'fe58'
ccc      ionam(ico55) = 'co55'
ccc      ionam(ini56) = 'ni56'
ccc      ionam(ini58) = 'ni58'
ccc      ionam(ineut) = 'neut'
ccc      ionam(iprot) = 'prot'
      ionam(ihe4)  = 'he4 '
      ionam(ic12)  = 'c12 '
      ionam(io16)  = 'o16 '
      ionam(ine20) = 'ne20'
      ionam(img24) = 'mg24'
      ionam(isi28) = 'si28'
      ionam(is32)  = 's32 '
      ionam(iar36) = 'ar36'
      ionam(ica40) = 'ca40'
      ionam(iti44) = 'ti44'
      ionam(icr48) = 'cr48'
      ionam(ife52) = 'fe52'
      ionam(ini56) = 'ni56'



c..set the number of nucleons in the element
ccc      aion(ihe3)  = 3.0d0
ccc      aion(ihe4)  = 4.0d0
ccc      aion(ic12)  = 12.0d0
ccc      aion(ic13)  = 13.0d0
ccc      aion(in13)  = 13.0d0
ccc      aion(in14)  = 14.0d0
ccc      aion(io14)  = 14.0d0
ccc      aion(io15)  = 15.0d0
ccc      aion(io16)  = 16.0d0
ccc      aion(io17)  = 17.0d0
ccc      aion(io18)  = 18.0d0
ccc      aion(if17)  = 17.0d0
ccc      aion(if18)  = 18.0d0
ccc      aion(if19)  = 19.0d0
ccc      aion(ine18) = 18.0d0
ccc      aion(ine19) = 19.0d0
ccc      aion(ine20) = 20.0d0
ccc      aion(img22) = 22.0d0
ccc      aion(img24) = 24.0d0
ccc      aion(ial27) = 27.0d0
ccc      aion(isi28) = 28.0d0
ccc      aion(ip31)  = 31.0d0
ccc      aion(is30)  = 30.0d0
ccc      aion(is32)  = 32.0d0
ccc      aion(icl35) = 35.0d0
ccc      aion(iar36) = 36.0d0
ccc      aion(ik39)  = 39.0d0
ccc      aion(ica40) = 40.0d0
ccc      aion(iti44) = 44.0d0
ccc      aion(icr48) = 48.0d0
ccc      aion(icr49) = 49.0d0
ccc      aion(icr50) = 50.0d0
ccc      aion(icr51) = 51.0d0
ccc      aion(icr52) = 52.0d0
ccc      aion(icr53) = 53.0d0
ccc      aion(icr54) = 54.0d0
ccc      aion(ife52) = 52.0d0  
ccc      aion(ife54) = 54.0d0  
ccc      aion(ife55) = 55.0d0  
ccc      aion(ife56) = 56.0d0  
ccc      aion(ife57) = 57.0d0  
ccc      aion(ife58) = 58.0d0  
ccc      aion(ico55) = 55.0d0
ccc      aion(ini56) = 56.0d0  
ccc      aion(ini58) = 58.0d0
ccc      aion(ineut) = 1.0d0
ccc      aion(iprot) = 1.0d0
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



c..set the number of protons in the element
ccc      zion(ihe3)  = 2.0d0
ccc      zion(ihe4)  = 2.0d0
ccc      zion(ic12)  = 6.0d0
ccc      zion(ic13)  = 6.0d0
ccc      zion(in13)  = 7.0d0
ccc      zion(in14)  = 7.0d0
ccc      zion(io14)  = 8.0d0
ccc      zion(io15)  = 8.0d0
ccc      zion(io16)  = 8.0d0
ccc      zion(io17)  = 8.0d0
ccc      zion(io18)  = 8.0d0
ccc      zion(if17)  = 9.0d0
ccc      zion(if18)  = 9.0d0
ccc      zion(if19)  = 9.0d0
ccc      zion(ine18) = 10.0d0
ccc      zion(ine19) = 10.0d0
ccc      zion(ine20) = 10.0d0
ccc      zion(img22) = 12.0d0
ccc      zion(img24) = 12.0d0
ccc      zion(ial27) = 13.0d0
ccc      zion(isi28) = 14.0d0
ccc      zion(ip31)  = 15.0d0
ccc      zion(is30)  = 16.0d0
ccc      zion(is32)  = 16.0d0
ccc      zion(icl35) = 17.0d0
ccc      zion(iar36) = 18.0d0
ccc      zion(ik39)  = 19.0d0
ccc      zion(ica40) = 20.0d0
ccc      zion(iti44) = 22.0d0
ccc      zion(icr48) = 24.0d0
ccc      zion(icr49) = 24.0d0
ccc      zion(icr50) = 24.0d0
ccc      zion(icr51) = 24.0d0
ccc      zion(icr52) = 24.0d0
ccc      zion(icr53) = 24.0d0
ccc      zion(icr54) = 24.0d0
ccc      zion(ife52) = 26.0d0  
ccc      zion(ife54) = 26.0d0  
ccc      zion(ife55) = 26.0d0  
ccc      zion(ife56) = 26.0d0  
ccc      zion(ife57) = 26.0d0  
ccc      zion(ife58) = 26.0d0  
ccc      zion(ico55) = 27.0d0
ccc      zion(ini56) = 28.0d0  
ccc      zion(ini58) = 28.0d0  
ccc      zion(ineut) = 0.0d0
ccc      zion(iprot) = 1.0d0
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


c..set the number of neutrons
       do i=1,ionmax
        nion(i) = aion(i) - zion(i)
       enddo



c..set the binding energy of the element
ccc      bion(ihe3)  = 7.71819d0
ccc      bion(ihe4)  = 28.29603d0 
ccc      bion(ic12)  = 92.16294d0
ccc      bion(ic13)  = 97.1060d0
ccc      bion(in13)  = 94.1030d0
ccc      bion(in14)  = 104.65998d0
ccc      bion(io14)  = 98.7310d0
ccc      bion(io15)  = 111.9530d0
ccc      bion(io16)  = 127.62093d0
ccc      bion(io17)  = 131.7600d0
ccc      bion(io18)  = 139.8040d0 
ccc      bion(if17)  = 128.2170d0
ccc      bion(if18)  = 137.3670d0
ccc      bion(if19)  = 147.7980d0
ccc      bion(ine18) = 132.1390d0
ccc      bion(ine19) = 143.7780d0
ccc      bion(ine20) = 160.64788d0 
ccc      bion(img22) = 168.5750d0
ccc      bion(img24) = 198.2579d0
ccc      bion(ial27) = 224.9480d0
ccc      bion(isi28) = 236.5379d0
ccc      bion(ip31)  = 262.9120d0
ccc      bion(is30)  = 243.6810d0
ccc      bion(is32)  = 271.7825d0
ccc      bion(icl35) = 298.2050d0
ccc      bion(iar36) = 306.7202d0
ccc      bion(ik39)  = 333.7180d0
ccc      bion(ica40) = 342.0568d0
ccc      bion(iti44) = 375.4772d0
ccc      bion(icr48) = 411.469d0
ccc      bion(icr49) = 422.0370d0
ccc      bion(icr50) = 435.0370d0
ccc      bion(icr51) = 444.2980d0
ccc      bion(icr52) = 456.3370d0
ccc      bion(icr53) = 464.2760d0
ccc      bion(icr54) = 473.9950d0
ccc      bion(ife52) = 447.708d0
ccc      bion(ife54) = 471.7696d0
ccc      bion(ife55) = 481.0480d0
ccc      bion(ife56) = 492.2450d0
ccc      bion(ife57) = 499.8910d0
ccc      bion(ife58) = 509.9350d0
ccc      bion(ico55) = 476.8150d0
ccc      bion(ini56) = 484.003d0 
ccc      bion(ini58) = 506.4450d0
ccc      bion(ineut) = 0.0d0
ccc      bion(iprot) = 0.0d0
      bion(ihe4)  = 28.29603d0 
      bion(ic12)  = 92.16294d0
      bion(io16)  = 127.62093d0
      bion(ine20) = 160.64788d0 
      bion(img24) = 198.2579d0
      bion(isi28) = 236.5379d0
      bion(is32)  = 271.7825d0
      bion(iar36) = 306.7202d0
      bion(ica40) = 342.0568d0
      bion(iti44) = 375.4772d0
      bion(icr48) = 411.469d0
      bion(ife52) = 447.708d0
      bion(ini56) = 484.003d0 

      return
      end


