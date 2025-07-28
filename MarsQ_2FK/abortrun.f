      subroutine abortrun(routine,ident,errmsg,vnam0,ier0,vnam1,ier1
     &                                        ,iw,kunit)
c
c ---------------------------------------------------------
c
c  Print warning or abort job if error detected
c
c ---------------------------------------------------------
c
      character*16  codenam
      character*(*) routine
      character*(*) vnam0, vnam1
      character*(*) errmsg
c
      character*16 filnam
      logical      filxist, filopen
c
c
c
c 1.0 Define the error label for termination
c
c 1.1 Set file names
c
      codenam  = 'fftdriver'
      filnam   = 'Error_Messages'
c
c
c 1.1 open the Message file
c
      inquire(file   = filnam, opened = filopen
     &                       , exist  = filxist
     &                       , iostat = ier2)
      if(.not. filxist) filopen  = .false.
      if(      filxist) then
         if(filopen) then
            close(unit = kunit)
            filopen  = .false.
         endif
      endif
c
      open(unit =  kunit, file   =  filnam,     form     = 'formatted'
     &                  , status = 'unknown',   position = 'append')
c
c
c
c 2.0 Write out error label and number
c
c 2.1 Serious warning
c
      if    (iw .lt. -2) then
         write(kunit,1000) iw,routine,ident,errmsg,vnam0,ier0,vnam1,ier1
         close(unit = kunit)
         return
c
c
c 2.2 Serious warning
c
      elseif(iw .eq. -2) then
         write(kunit,1100) routine,ident,errmsg,vnam0,ier0,vnam1,ier1
         close(unit = kunit)
         return
c
c
c 2.3 Important warning
c
      elseif(iw .eq. -1) then
         write(kunit,1200) routine,ident,errmsg,vnam0,ier0,vnam1,ier1
         close(unit = kunit)
         return
c
c
c 2.4 Warning
c
      elseif(iw .eq.  0) then
         write(kunit,1300) routine,ident,errmsg,vnam0,ier0,vnam1,ier1
         close(unit = kunit)
         return
c
c
c 2.5 Fatal error
c
      elseif(iw .eq. +1) then
         write(kunit,2000) codenam
         write(kunit,2010) routine,ident,errmsg,vnam0,ier0,vnam1,ier1
         close(unit = kunit)
         stop 'Fatal error'
c
c
c 2.6 Completed calculation
c
      elseif(iw .gt. +1) then
         write(kunit,2100) codenam
         write(kunit,2110) errmsg,vnam0,ier0,vnam1,ier1
         close(unit = kunit)
         stop 'Calculation complete'
      endif
c
c
c
c 3.0 Return and end
c
      return
c
 1000 format(/,2x,'***** Serious Warning: status number'
     &        ,2x,i2,1x,'*****'
     &      ,/,8x,'Non fatal error in subroutine',2x,a8,1x,'at location'
     &        ,1x,i4,/
     &        ,8x,a36,' :',4x,a8,' = ',i8,4x,a8,' = ',i8,/)
 1100 format(/,2x,'***** Serious Warning:  '
     &        ,1x,'Non fatal error in subroutine',2x,a8,1x,'at location'
     &        ,1x,i4,1x,'*****',/
     &        ,8x,a36,' :',4x,a8,' = ',i8,4x,a8,' = ',i8,/)
 1200 format(/,8x,'Important Warning:'
     &        ,1x,'Non fatal error in subroutine',2x,a8,1x,'at location'
     &        ,1x,i4,/
     &        ,8x,a36,' :',4x,a8,' = ',i8,4x,a8,' = ',i8,/)
 1300 format(/,8x,'Warning:'
     &        ,1x,'Non fatal error in subroutine',2x,a8,1x,'at location'
     &        ,1x,i4,/
     &        ,8x,a36,' :',4x,a8,' = ',i8,4x,a8,' = ',i8)
 2000 format(/,2x,'***** Aborting  run in',1x,a16,1x,'*****')
 2010 format(/,8x,'Fatal error in subroutine    ',2x,a8,1x,'at location'
     &        ,1x,i4,/
     &        ,1x,'Error: ',a36,' :',4x,a8,' = ',i8,4x,a8,' = ',i8,/)
 2020 format(    'fatal')
 2100 format(/,2x,'***** Completed run in',1x,a16,1x,'*****')
 2110 format(/,8x,a36,' :',4x,a8,' = ',i8,4x,a8,' = ',i8,/)
      end
