 subroutine siesta_dminit
 use precision,       only : dp
 use sys,             only : die
 use siesta_geom,     only : na_u
 use atomlist,        only : no_u, no_l, iaorb, lasto
 use alloc,           only : re_alloc, de_alloc
 use files,           only : slabel
 use m_variables,     only : absocctol
 use m_variables,     only : natot, nuotot, nuoloc, nuoocc, nuolococc
 use m_variables,     only : eig, occ
 use m_variables,     only : Saux, Haux, Daux, wfk
 use m_variables,     only : wfk_red
 use m_variables,     only : nuoind, spnind
 use parallel,        only : Node, Nodes
 use parallelsubs,    only : LocalToGlobalOrb
 use sparse_matrices, only : maxnh, numh, listh, listhptr
 use sparse_matrices, only : H, S, Dscf
 use densematrix,     only : psi
 use m_eo,            only : eo, qo
 use m_gamma,         only : gamma
 use m_spin,          only : nspin
 implicit none
! Local variables
 integer  :: ispin, ia, io, jo, j, ind, iio, iu, ipsi
 external :: collectdata_dp
!
 call timer( "dminit", 1 )
! Gamma point only
 if(.not.gamma) &
  call die('ERROR: TD dynamics has not been implemented for multiple k-points yet!')
!
  natot = na_u
 nuotot = no_u
 nuoloc = no_l
!
 if(associated(eig)) call de_alloc(eig,name="eig")
 call re_alloc(eig,1,nuotot,1,nspin,name="eig", routine="siesta_dminit")
 if(associated(occ)) call de_alloc(occ,name="occ")
 call re_alloc(occ,1,nuotot,1,nspin,name="occ", routine="siesta_dminit")
!
 if(associated(Saux)) call de_alloc(Saux,name="Saux")
 call re_alloc(Saux,1,nuotot,1,nuoloc,name="Saux",routine="siesta_dminit")
!
 if(associated(Haux)) call de_alloc(Haux,name="Haux")
 call re_alloc(Haux,1,nuotot,1,nuoloc,1,nspin,name="Haux",routine="siesta_dminit")
 if(associated(Daux)) call de_alloc(Daux,name="Daux")
 call re_alloc(Daux,1,nuotot,1,nuoloc,1,nspin,name="Daux",routine="siesta_dminit")
!
 if(associated(wfk))  call de_alloc(wfk,name="wfk")
 call re_alloc(wfk,1,nuotot,1,nuoloc,1,nspin,name="wfk",routine="siesta_dminit")
! Gamma point eigen-energies and occupancies
 do ispin = 1,nspin
    do io = 1,nuotot
       eig(io,ispin) = eo(io,ispin,1)
       occ(io,ispin) = qo(io,ispin,1)
       if(occ(io,ispin).gt.absocctol) nuoocc = io
    enddo
 enddo
!
 nuolococc = 0
 do io = 1,nuoloc
    call LocalToGlobalOrb(io, Node, Nodes, iio)
    if(iio.le.nuoocc) nuolococc = io
 enddo
! Restore dense from sparse
 do ispin = 1,nspin
    do io = 1,nuoloc
       do jo = 1,nuotot
          if(ispin.eq.1) &
          Saux(jo,io)       = 0.0_dp
          Haux(jo,io,ispin) = 0.0_dp
          Daux(jo,io,ispin) = 0.0_dp
       enddo
    enddo
    do io = 1,nuoloc
       do j = 1,numh(io)
          ind = listhptr(io) + j
           jo = listh(ind)
          if(ispin.eq.1) &
          Saux(jo,io)       = Saux(jo,io)       + S(ind)
          Haux(jo,io,ispin) = Haux(jo,io,ispin) + H(ind,ispin)
          Daux(jo,io,ispin) = Daux(jo,io,ispin) + Dscf(ind,ispin)
       enddo
    enddo
 enddo
 ipsi = 1
 do ispin = 1,nspin
    do io = 1,nuoloc
       do jo = 1,nuotot
          wfk(jo,io,ispin) = psi(ipsi)
          ipsi = ipsi + 1
       enddo
    enddo
 enddo
!
 if(associated(wfk_red)) call de_alloc(wfk_red,name="wfk_red")
 call re_alloc(wfk_red,1,nuotot,1,nuotot,1,nspin,name="wfk_red",routine="siesta_dminit")
!
 call collectdata_dp(nuotot, nuoloc, nspin, wfk, wfk_red)
!
 do ispin = 1,nspin
    do io = 1,nuotot
       if(ispin.eq.spnind .and. io.eq.nuoind) then
         if(Node.eq.0) then
           open(unit=50,file=trim(slabel)//'.WFKOUT_ascii')
           do jo = 1,nuotot
              write(50,'(e20.7)') wfk_red(jo,io,ispin)
           enddo
           close(50)
         endif
       endif
    enddo
 enddo
!
 call timer( "dminit", 2 )
!
 end subroutine siesta_dminit
