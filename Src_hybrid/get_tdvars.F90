 subroutine get_tdvars
 use parallel,    only : Node
 use files,       only : label_length, slabel
 use alloc,       only : re_alloc
 use m_spin,      only : nspin
 use m_variables, only : simplempi, fieldtype, tdproptyp
 use m_variables, only : maxnstb, maxntim, abstolstb
 use m_variables, only : tadiabatic, eext, estrength, dtim, sigma
 use m_variables, only : tdhyb, speconly, pseudoatm
 use m_variables, only : nonscfdm, nonscfwf
 use m_variables, only : chkspows, chknormdm, chkstabl
 use m_variables, only : homolvl, lumolvl, homochg, lumochg
 use m_variables, only : nuowfkadd, nuooffset, wfkadd, addchg
 use m_variables, only : writetdpop
 use m_variables, only : nuoind, spnind
#ifdef MPI
 use mpi_siesta
#endif 
 implicit none
 integer :: ispin
#ifdef MPI
 integer :: MPIerror
#endif
 external :: redtdv
!
 if(Node.eq.0) call redtdv()
#ifdef MPI
 call MPI_Bcast(   slabel,label_length,MPI_character,       0,MPI_Comm_World,MPIerror)
 call MPI_Bcast(simplempi,           1,MPI_logical,         0,MPI_Comm_World,MPIerror)
 call MPI_Bcast(  maxnstb,           1,MPI_integer,         0,MPI_Comm_World,MPIerror)
 call MPI_Bcast(  maxntim,           1,MPI_integer,         0,MPI_Comm_World,MPIerror)
 call MPI_Bcast(   nuoind,           1,MPI_integer,         0,MPI_Comm_World,MPIerror)
 call mpi_bcast(   spnind,           1,MPI_integer,         0,mpi_comm_world,mpierror)
 call mpi_bcast(fieldtype,           1,MPI_integer,         0,mpi_comm_world,mpierror)
 call mpi_bcast(tdproptyp,           1,MPI_integer,         0,mpi_comm_world,mpierror)
 call MPI_Bcast(abstolstb,           1,MPI_double_precision,0,MPI_Comm_World,MPIerror)
 call MPI_Bcast(     eext,           3,MPI_double_precision,0,MPI_Comm_World,MPIerror)
 call MPI_Bcast(estrength,           1,MPI_double_precision,0,MPI_Comm_World,MPIerror)
 call MPI_Bcast(tadiabatic,          1,MPI_double_precision,0,MPI_Comm_World,MPIerror)
 call MPI_Bcast(     dtim,           1,MPI_double_precision,0,MPI_Comm_World,MPIerror)
 call MPI_Bcast(    sigma,           1,MPI_double_precision,0,MPI_Comm_World,MPIerror)
 call MPI_Bcast(    tdhyb,           1,MPI_logical,         0,MPI_Comm_World,MPIerror)
 call MPI_Bcast( speconly,           1,MPI_logical,         0,MPI_Comm_World,MPIerror)
 call MPI_Bcast(pseudoatm,           1,MPI_logical,         0,MPI_Comm_World,MPIerror)
 call MPI_Bcast( nonscfdm,           1,MPI_logical,         0,MPI_Comm_World,MPIerror)
 call MPI_Bcast( nonscfwf,           1,MPI_logical,         0,MPI_Comm_World,MPIerror)
 call MPI_Bcast( chkspows,           1,MPI_logical,         0,MPI_Comm_World,MPIerror)
 call MPI_Bcast(chknormdm,           1,MPI_logical,         0,MPI_Comm_World,MPIerror)
 call MPI_Bcast( chkstabl,           1,MPI_logical,         0,MPI_Comm_World,MPIerror)
 call MPI_Bcast(writetdpop,          1,mpi_logical,         0,mpi_comm_world,mpierror)
 if(Node.ne.0 .and. nonscfdm) then
   if(.not.associated(homolvl))   &
    call re_alloc(homolvl,   1, nspin, name='homolvl',   routine='get_tdvars')
   if(.not.associated(lumolvl))   &
    call re_alloc(lumolvl,   1, nspin, name='lumolvl',   routine='get_tdvars')
   if(.not.associated(homochg))   &
    call re_alloc(homochg,   1, nspin, name='homochg',   routine='get_tdvars')
   if(.not.associated(lumochg))   &
    call re_alloc(lumochg,   1, nspin, name='lumochg',   routine='get_tdvars')
 endif
 if(nonscfdm) then
   call MPI_Bcast(homolvl(1), nspin, MPI_integer,         0,MPI_Comm_World,MPIerror)
   call MPI_Bcast(homochg(1), nspin, MPI_double_precision,0,MPI_Comm_World,MPIerror)
   call MPI_Bcast(lumolvl(1), nspin, MPI_integer,         0,MPI_Comm_World,MPIerror)
   call MPI_Bcast(lumochg(1), nspin, MPI_double_precision,0,MPI_Comm_World,MPIerror)
 endif
 if(Node.ne.0 .and. nonscfwf) then
   if(.not.associated(nuowfkadd)) &
    call re_alloc(nuowfkadd, 1, nspin, name='nuowfkadd', routine='get_tdvars')
   if(.not.associated(nuooffset)) &
    call re_alloc(nuooffset, 1, nspin, name='nuooffset', routine='get_tdvars')
   if(.not.associated(addchg))    &
    call re_alloc(addchg,    1, nspin, name='chgadd',    routine='get_tdvars')
 endif
 if(nonscfwf) then
   call MPI_Bcast(nuowfkadd(1), nspin, MPI_integer,         0,MPI_Comm_World,MPIerror)
   call MPI_Bcast(nuooffset(1), nspin, MPI_integer,         0,MPI_Comm_World,MPIerror)
   call MPI_Bcast(   addchg(1), nspin, MPI_double_precision,0,MPI_Comm_World,MPIerror)
 endif
 if(Node.ne.0 .and. nonscfwf) then
   if(.not.associated(wfkadd) .and. maxval(nuowfkadd).gt.0) &
    call re_alloc(wfkadd, 1, maxval(nuowfkadd), 1, nspin, name='wfkadd', routine='get_tdvars')
 endif
 if(nonscfwf) then
   if(maxval(nuowfkadd).gt.0) &
   call MPI_Bcast(wfkadd(1,1),maxval(nuowfkadd)*nspin,MPI_double_precision,0,MPI_Comm_World,MPIerror)
 endif
#endif
!
 end subroutine get_tdvars
