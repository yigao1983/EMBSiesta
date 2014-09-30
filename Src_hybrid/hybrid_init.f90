 subroutine hybrid_init
 use sys,            only : die
 use parallel,       only : Node
 use m_variables,    only : pseudoatm, pseudonfr, stagenfr
 use siesta_options, only : maxsav
 use nearfieldsubs,  only : ConfigMater
 implicit none
!
 if(Node.eq.0) write(6,'(/a)') 'hybrid_init: Reading nearfield and TDDFT input...'
! 
 call get_nfvars()
 call get_tdvars()
!
 if(pseudoatm .and. pseudonfr) &
  call die('ERROR: Atoms and nearfield cannot be both pseudo!')
!
 if(maxsav.lt.1) maxsav = 1
 stagenfr = 'scf'
!
 call ConfigMater()
!
 end subroutine hybrid_init
