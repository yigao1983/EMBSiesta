 subroutine renormdaux(Daux)
 use precision,   only : dp
 use sys,         only : die
 use alloc,       only : re_alloc, de_alloc
 use parallel,    only : Node, Nodes
 use atomlist,    only : qtot
 use parallelsubs,only : LocalToGlobalOrb
 use m_spin,      only : nspin
 use m_variables, only : absocctol
 use m_variables, only : nonscfdm, nonscfwf
 use m_variables, only : natot, nuotot, nuoloc, nuoocc, nuolococc
 use m_variables, only : homolvl, homochg, lumolvl, lumochg
 use m_variables, only : occ, wfk
 use m_variables, only : nuowfkadd, nuooffset, addchg, wfkadd
 use m_variables, only : Saux, Saux_red
 implicit none
 real(dp), intent(out) :: Daux(nuotot,nuoloc,nspin)
! Local variables
 integer  :: ispin, ie, io, jo, ko, lo, iio
 integer  :: ioadd(nspin)
 real(dp) :: norm, qnew
 real(dp) :: popu(natot)
!
 call timer( 'renormdaux', 1 )
!
 if(nonscfdm) then
   do ispin = 1,nspin
      do ie = 1,nuotot
         if(ie.eq.homolvl(ispin)) occ(ie,ispin) = occ(ie,ispin) - homochg(ispin)
         if(ie.eq.lumolvl(ispin)) occ(ie,ispin) = occ(ie,ispin) + lumochg(ispin)
      enddo
   enddo
 endif
!
 if(.not.associated(Saux_red)) call re_alloc(Saux_red, 1, nuotot, 1, nuotot, name='Saux_red', routine='renormdaux')
 call collectdata_dp(nuotot, nuoloc, nspin, Saux, Saux_red)
!
 if(nonscfwf) then
   if(.not.associated(wfkadd)) then
     call die('added wavefunction not associated')
   elseif(.not.associated(nuowfkadd) .or. maxval(nuowfkadd).le.0) then
     call die('0-orbital added wavefunction')
   endif
   do ispin = 1,nspin
      ie = 1
      do while(ie.le.nuotot .and. dabs(occ(ie,ispin)).gt.1.0e-9_dp)
         ie = ie + 1
      enddo
      ioadd(ispin) = min(ie,nuotot)
      do io = 1,nuoloc
         call LocalToGlobalOrb(io, Node, Nodes, iio)
         if(iio.eq.ioadd(ispin)) then
           wfk(:,io,ispin) = 0.0_dp
           do jo = 1,nuowfkadd(ispin)
              wfk(jo+nuooffset(ispin),io,ispin) = wfkadd(jo,ispin)
           enddo
!          norm = 0.0_dp
!          do ko = 1,nuotot
!             do lo = 1,nuotot
!                norm = norm + wfk(lo,io,ispin)*Saux_red(lo,ko)*wfk(ko,io,ispin)
!             enddo
!          enddo
           norm = sum(wfk(:,io,ispin)*matmul(Saux_red,wfk(:,io,ispin)))
           if(dabs(norm).lt.1.e-9_dp) call die('0 wavefunction added')
           wfk(:,io,ispin) = 1.0_dp/dsqrt(norm) * wfk(:,io,ispin)
         endif
      enddo
      occ(ioadd(ispin),ispin) = occ(ioadd(ispin),ispin) + addchg(ispin)
   enddo
 endif
! Gamma point eigen-energies and occupancies
 do ispin = 1,nspin
    do io = 1,nuotot
       if(occ(io,ispin).gt.absocctol) nuoocc = io
    enddo
 enddo
!
 nuolococc = 0
 do io = 1,nuoloc
    call LocalToGlobalOrb(io, Node, Nodes, iio)
    if(iio.le.nuoocc) nuolococc = io
 enddo
!
 call psi2daux(wfk, Daux)
 call daux2dsparse(Daux)
!
 call mullikpop(nuotot, nuoloc, nspin, Daux, natot, popu)
!
 qnew = sum(popu)
 if(abs(qnew-qtot) .gt. 1.0e-9_dp) then
   qtot = qnew
   if(Node.eq.0) write(6,'(/a,e20.7)') 'Total number of electrons changed to: ', qnew
 endif
!
 call de_alloc(Saux_red, name='Saux_red')
!
 call timer( 'renormdaux', 2 )
!
 end subroutine renormdaux
