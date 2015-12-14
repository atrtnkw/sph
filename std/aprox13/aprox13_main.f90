subroutine setup_aprox13()
  include 'implno.dek'
  include 'network.dek'

  integer          nok,nbad
  double precision time_start,time_end,conserv,tin,din,ein,vin,zin,xin(13), &
       tout,dout,eout,xout(13),dtsav

  call init_aprox13
  call read_helm_table
  call net_input0(time_start,time_end,tin,din,vin,zin,ein,xin)
  hydrostatic = .true.
!  one_step = .true.

!  write(*,*) "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx"
  write(*,*) "Start nuclear reaction."
end subroutine setup_aprox13

subroutine solve_aprox13(tstep, tin, din, xin, deout, mode)
  include 'implno.dek'
  include 'network.dek'

  double precision tstart,tstep,conserv,tin,din,ein,vin,zin,xin(18)
  double precision tout,dout,eout,xout(18),edum,dtsav,deout
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
!  deout         = sdot - sneut
  deout         = eout

  return
end subroutine solve_aprox13
