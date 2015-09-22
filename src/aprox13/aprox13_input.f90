!---------------------------------------------------------------------
      subroutine net_input0(time_start,time_end,tin,din,vin,zin,ein,xin)
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

! limit the temperature since the rates are invalid much above t9=100
       tin = min(1.0d11,tin)

!---------------------------------------------------------------------------

! there is probably a better place for this
! if requested, adjust the number of equations being solved

      if (pure_network .eq. 1) then
       neqs  = ionmax
       btemp = tin
       bden  = din
      end if

      return
      end
!---------------------------------------------------------------------
