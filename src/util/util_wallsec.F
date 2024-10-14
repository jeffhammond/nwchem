      double precision function util_wallsec_timer()
      implicit none
#ifdef NO_MPIF
#include "tcgmsg.fh"
      util_wallsec_timer = TCGTIME()
#else
      double precision MPI_Wtime
      external MPI_Wtime
      util_wallsec_timer = MPI_Wtime()
#endif
      end function

      double precision function util_wallsec()
      implicit none
!
! routine to give wall clock seconds since start of execution
! delta times are handled by all applications
!
      logical first_call
      double precision first_time,cur_time

      save first_time
      save first_call
      data first_call /.true./

! why is this here? it is overwritten later...
      util_wallsec = util_wallsec_timer()

      if (first_call) then
         first_time = util_wallsec_timer()
         first_call = .false.
      endif

      cur_time = util_wallsec_timer()
      util_wallsec = cur_time - first_time

      end
