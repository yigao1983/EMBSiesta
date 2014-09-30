      Program HYBSIESTA

      use m_siesta_init
      use m_siesta_analysis
      use m_siesta_move
      use m_siesta_end
      use m_siesta_forces
      USE m_steps, only: inicoor, fincoor
      ! Yi: Determine if calculating spectrum only
      use m_variables, only: speconly

      implicit none

      integer :: istep
      logical :: relaxd

!------------------------------------------------------------------------- BEGIN
      call siesta_init()
      call hybrid_init()

      ! Begin of coordinate relaxation iteration
      relaxd = .false.
      istep  = inicoor
      DO WHILE ((istep.le.fincoor) .AND. (.not. relaxd))
        ! Yi: If not only spectra, do GS-SCF
        if (.not. speconly) call siesta_forces( istep )
        ! Yi: For a specific configuration, TD dynamics
        call siesta_tdhyb()
        ! Yi: If not only spectra, do relaxation
        if (.not. speconly) then

           call siesta_move( istep, relaxd )
           if (.not. relaxd) then
             istep = istep + 1
           endif
        ! Yi: Otherwise relaxed to quit loop
        else

           relaxd = .true.

        endif

      ENDDO
      ! End of coordinate-relaxation loop 

      call siesta_analysis( relaxd )

      call siesta_end()
!--------------------------------------------------------------------------- END
      END program hybsiesta
