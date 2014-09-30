 subroutine siesta_tdhyb
 use m_variables, only : tdhyb, tdproptyp, speconly
 implicit none
!
 call timer( 'tdhyb', 1 )
!
 if(tdhyb) then
   call siesta_dminit()
   call siesta_tdinit()
   if(tdproptyp.eq.1) then
     call siesta_evolut()
   else
     call siesta_evolut_wfk()
   endif
 endif
!
 if(tdhyb.or.speconly) call siesta_spectra()
!
 call timer( 'tdhyb', 2 )
!
 end subroutine siesta_tdhyb
