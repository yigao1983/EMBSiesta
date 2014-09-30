 subroutine siesta_tdinit
 use m_variables, only : chkspows
 implicit none
!
 call timer( 'tdinit', 1 )
!
 call overlap_powers()
!
 if(chkspows) call check_spowers()
!
 call state_stabiliz()
!
 call timer( 'tdinit', 2 )
!
 end subroutine siesta_tdinit
