subroutine setup_iso7()
  include 'implno.dek'
  include 'network.dek'

  integer          nok,nbad
  double precision time_start,time_end,conserv,tin,din,ein,vin,zin,xin(7), &
       tout,dout,eout,xout(7)

  call init_iso7
  call read_helm_table
  call net_input(time_start,time_end,tin,din,vin,zin,ein,xin)
  hydrostatic = .true.
!  one_step = .true.

!  write(*,*) "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx"
end subroutine setup_iso7

subroutine solve_iso7(tstep, tin, din, xin, eout, mode)
  include 'implno.dek'
  include 'network.dek'

  double precision tstart,tstep,conserv,tin,din,ein,vin,zin,xin(18)
  double precision tout,dout,eout,xout(18),edum,dtsav
  integer          nok,nbad,mode
  double precision abar,zbar,wbar,xcess,ye_orig

  tstart = 0.0d0
  dtsav  = 0.0d0
  vin    = 0.0d0
  zin    = 0.0d0
  ein    = 0.0d0

  call burner(tstart,tstep,dtsav, &
       tin,din,vin,zin,ein,xin, &
       tout,dout,eout,xout, &
       conserv,nok,nbad)

  xin(1:ionmax) = xout(1:ionmax)

  return
end subroutine solve_iso7
