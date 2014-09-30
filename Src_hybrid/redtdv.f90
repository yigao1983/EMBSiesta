 subroutine redtdv
 use fdf
 use precision,   only : dp
 use sys,         only : die
 use parse,       only : parsed_line, digest, match
 use parse,       only : names, values
 use files,       only : slabel
 use alloc,       only : re_alloc
 use m_spin,      only : nspin
 use m_variables, only : slabel_default, simplempi_default, fieldtype_default
 use m_variables, only : maxnstb_default, maxntim_default, maxntim_default, tdproptyp_default
 use m_variables, only : abstolstb_default, tdhyb_default, dtim_default, sigma_default
 use m_variables, only : speconly_default, pseudoatm_default
 use m_variables, only : tadiabatic_default, eext_default
 use m_variables, only : nonscfdm_default, nonscfwf_default
 use m_variables, only : chkspows_default, chknormdm_default, chkstabl_default
 use m_variables, only : homolvl_default, homochg_default
 use m_variables, only : lumolvl_default, lumochg_default
 use m_variables, only : nuowfkadd_default, nuooffset_default, addchg_default
 use m_variables, only : writetdpop_default
 use m_variables, only : nuoind_default, spnind_default, addchg_default
 use m_variables, only : simplempi, fieldtype
 use m_variables, only : maxnstb, maxntim, maxntim, tdproptyp
 use m_variables, only : abstolstb, tdhyb, dtim, sigma
 use m_variables, only : speconly, pseudoatm, tadiabatic, eext, estrength
 use m_variables, only : nonscfdm, nonscfwf
 use m_variables, only : chkspows, chknormdm, chkstabl
 use m_variables, only : homolvl, homochg, lumolvl, lumochg
 use m_variables, only : nuowfkadd, nuooffset, wfkadd, addchg
 use m_variables, only : writetdpop
 use m_variables, only : nuoind, spnind
 implicit none
! Local variable
 integer  :: iu, ispin, io, ix
 real(dp) :: cfactor
 logical  :: found
 type(block), pointer        :: bp => null()
 type(parsed_line), pointer  :: p  => null()
 character(len=132)          :: line, eunits
!
    slabel = fdf_string( 'SystemLabel',            slabel_default)
!
 simplempi = fdf_boolean('TDSimpleMPI',         simplempi_default)
!
   maxnstb = fdf_integer('MaxStabilizIterations', maxnstb_default)
   maxntim = fdf_integer('MaxEvolutionSteps',     maxntim_default)
    nuoind = fdf_integer('OutputOrbitalIndex',     nuoind_default)
    spnind = fdf_integer('OutputSpinIndex',        spnind_default)
 fieldtype = fdf_integer('TDElectricFieldType', fieldtype_default)
 tdproptyp = fdf_integer('TDPropagationType',   tdproptyp_default)
 abstolstb = fdf_double( 'Stabiliz.Tolerance',  abstolstb_default)
     tdhyb = fdf_boolean('TDHybridDynamics',        tdhyb_default)
tadiabatic = fdf_double( 'AdiabaticTime',      tadiabatic_default)
      dtim = fdf_double( 'Evolution.TimeStep',       dtim_default)
     sigma = fdf_double( 'Evolution.DecayRate',     sigma_default)
  speconly = fdf_boolean('SpectraOnly',          speconly_default)
 pseudoatm = fdf_boolean('PseudoAtoms',         pseudoatm_default)
!
  nonscfdm = fdf_boolean('NonSCFDM',             nonscfdm_default)
  nonscfwf = fdf_boolean('NonSCFWF',             nonscfwf_default)
  chkspows = fdf_boolean('CheckOverlapPowers',   chkspows_default)
 chknormdm = fdf_boolean('CheckNormDM',         chknormdm_default)
  chkstabl = fdf_boolean('CheckStablization',    chkstabl_default)
 writetdpop= fdf_boolean('WriteTDPopulations', writetdpop_default)
!
 if(tdproptyp.ne.1) simplempi = .true.
!
!if(fdf_block('ElectricField', iu)) then
!  read(iu,*) eext
!else
!  eext = eext_default
!endif
!estrength = dsqrt(sum(eext(:)*eext(:)))
!
 found = fdf_block('TDExternalElectricField',bp)
 if (found) then
   loop: DO
     if (.not. fdf_bline(bp,line)) exit loop
     p => digest(line)
     if (.not. match(p,"vvvn")) &
       call die("Wrong format in TDElectricField block")
     eunits = names(p,1)
     cfactor = fdf_convfac(eunits,'Ry/Bohr/e')
     do ix = 1,3
        eext(ix) = values(p,ix) * cfactor
     enddo
   enddo loop
 else
   eext(:) = eext_default
 endif ! (found)
 estrength = dsqrt(sum(eext(:)*eext(:)))
!
 if(nonscfdm) then
   if(.not.associated(homolvl))   &
    call re_alloc(homolvl, 1, nspin, name='homolvl',   routine='redtdv')
   if(.not.associated(lumolvl))   &
    call re_alloc(lumolvl, 1, nspin, name='lumolvl',   routine='redtdv')
   if(.not.associated(homochg))   &
    call re_alloc(homochg, 1, nspin, name='homochg',   routine='redtdv')
   if(.not.associated(lumochg))   &
    call re_alloc(lumochg, 1, nspin, name='lumochg',   routine='redtdv')
 endif
 if(nonscfwf) then
   if(.not.associated(nuowfkadd)) &
    call re_alloc(nuowfkadd, 1, nspin, name='nuowfkadd', routine='redtdv')
   if(.not.associated(nuooffset)) &
    call re_alloc(nuooffset, 1, nspin, name='nuooffset', routine='redtdv')
   if(.not.associated(addchg))    &
    call re_alloc(addchg,    1, nspin, name='addchg',    routine='redtdv')
 endif
!
 if(nonscfdm) then
   if(fdf_block('HomoRenormalization', iu)) then
     do ispin = 1,nspin
        read(iu,*) homolvl(ispin), homochg(ispin)
     enddo
   else
     homolvl = homolvl_default
     homochg = homochg_default
   endif
!
   if(fdf_block('LumoRenormalization', iu)) then
     do ispin = 1,nspin
        read(iu,*) lumolvl(ispin), lumochg(ispin)
     enddo
   else
     lumolvl = lumolvl_default
     lumochg = lumochg_default
   endif
 endif
!
 if(nonscfwf) then
   if(fdf_block('AddWavefunction', iu)) then
     do ispin = 1,nspin
        read(iu,*) nuowfkadd(ispin), nuooffset(ispin), addchg(ispin)
     enddo
     if(.not.associated(wfkadd)) &
      call re_alloc(wfkadd, 1, maxval(nuowfkadd), 1, nspin, name='wfkadd', routine='redtdv')
     open(unit=50,file=trim(slabel)//'.WFKADD_ascii')
     do ispin = 1,nspin
        if(nuowfkadd(ispin).gt.0) then
          do io = 1,nuowfkadd(ispin)
             read(50,*) wfkadd(io,ispin)
          enddo
        endif
     enddo
     close(50)
   else
     nuowfkadd = nuowfkadd_default
     nuooffset = nuooffset_default
        addchg =    addchg_default
   endif
 endif
!
 end subroutine redtdv
