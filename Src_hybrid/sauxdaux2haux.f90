 subroutine sauxdaux2haux(iscf, newS, Saux, Daux, Haux)
 use precision,           only : dp
 use m_spin,              only : nspin
 use m_variables,         only : nuotot, nuoloc
 use m_setup_hamiltonian, only : setup_hamiltonian
 implicit none
 integer,  intent(in)  :: iscf
 logical,  intent(in)  :: newS
 real(dp), intent(in)  :: Saux(nuotot,nuoloc)
 real(dp), intent(in)  :: Daux(nuotot,nuoloc,nspin)
 real(dp), intent(out) :: Haux(nuotot,nuoloc,nspin)
! Local variables
 logical  :: first, last
!
 if(newS) call saux2ssparse(Saux)
!
 call daux2dsparse(Daux)
!
 first = .false.
  last = .false.
 call setup_hamiltonian( first, last, iscf )
!
 call hsparse2haux(Haux)
!
 end subroutine sauxdaux2haux
