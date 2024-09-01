      subroutine int_1e_time(ibas, jbas, integ_type)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "bas.fh"
#include "util.fh"
!     
!     arguments
!     
      integer ibas, jbas        ! [input] bra and ket basis sets
      character*(*) integ_type  ! [input] Name of integrals to compute
!     
!     local variables
!     
      integer nshell_i, nshell_j
      integer ishell, jshell, mem1, max1e
      integer ijshell, ilo, ihi, jlo, jhi
      integer l_buf, l_scr
      integer k_buf, k_scr
      integer type
      double precision :: t0, t1
      if (inp_compare(.false., integ_type, 'kinetic')) then
         type = 1
      else if (inp_compare(.false., integ_type, 'potential')) then
         type = 2
      else if (inp_compare(.false., integ_type, 'overlap')) then
         type = 3
      else if (inp_compare(.false., integ_type, 'pvp')) then
         type = 4
      else if (inp_compare(.false., integ_type, 'so'))then 
         type = 5 
      else if (inp_compare(.false., integ_type, 'cos_chg_pot'))then 
         type = 6 
      else if (inp_compare(.false., integ_type, 'bq_pot'))then 
         type = 7
      else if (inp_compare(.false., integ_type, 'potential_beta')) then
         type = 8
      else
         write(6,*) ' integ_type = ', integ_type
         call errquit('unknown integ_type', 0, INT_ERR)
      end if
      if (.not. bas_numcont(ibas, nshell_i))
     $   call errquit('bas_numcont failed for ibas', ibas,BASIS_ERR)
      if (.not. bas_numcont(jbas, nshell_j))
     $   call errquit('bas_numcont failed for jbas', jbas,BASIS_ERR)
      call int_mem_1e(max1e, mem1)
      if (type .eq. 5) then  
         max1e = 3*max1e 
         mem1 = 3*mem1 
      endif 
      if (.not. MA_push_get(MT_DBL,max1e,'buf',l_buf,k_buf))
     $   call errquit('MA failed', max1e, MA_ERR)
      if (.not. MA_push_get(MT_DBL, mem1,'scr',l_scr,k_scr))
     $   call errquit('MA failed', mem1, MA_ERR)

      t0 = util_wallsec()

      do jshell = 1, nshell_j
         do ishell = 1, nshell_i
            if (type .eq. 1) then
               call int_1eke (jbas, jshell, ibas, ishell, 
     $              mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf))
            else if (type .eq. 2) then
               call int_1epe (jbas, jshell, ibas, ishell, 
     $              mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf))
            else if (type .eq. 3) then
               call int_1eov (jbas, jshell, ibas, ishell, 
     $              mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf))
            else if (type .eq. 4) then
               call int_1epvpe (jbas, jshell, ibas, ishell, 
     $              mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf))
            else if (type .eq. 5) then 
               call intso_1e (jbas, jshell, ibas, ishell,
     $              mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf))
            else if (type .eq. 6) then 
               call int_1epot (jbas, jshell, ibas, ishell,
     $              mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf))
            else if (type .eq. 7) then 
               call int_1epot1 (jbas, jshell, ibas, ishell,
     $              mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf))
            else
               call errquit('invalid type?',type,CALC_ERR)
            end if
         end do
      end do

      t1 = util_wallsec()
      write(6,'(a20,f30.15)') 'int_1e_time: ', t1-t0

      if (.not. MA_pop_stack(l_scr)) 
     $   call errquit('pop failed', 0, GA_ERR)
      if (.not. MA_pop_stack(l_buf))
     $   call errquit('pop failed', 0, GA_ERR)
      end
