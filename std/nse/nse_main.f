      subroutine setup_nse()
      include 'implno.dek'
      
      call init_network
      
      end subroutine setup_nse
      
      subroutine solve_nse(tin,din,xout)
      include 'implno.dek'
      double precision tin,din,xout(47)
      double precision tt,dd,yy
      integer newguess,iprint
      
      tt       = tin
      dd       = din
      yy       = 0.5
      newguess = 1
      iprint   = 0
      call nse(tt,dd,yy,newguess,xout,iprint)
  
      return
      end subroutine solve_nse

