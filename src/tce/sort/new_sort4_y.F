      subroutine tce_sort_4_y(u,s,a,b,c,d,i,j,k,l,f,async)
      implicit none
      integer, intent(in) :: a,b,c,d,async
      integer, intent(in) :: i,j,k,l
      double precision, intent(out) :: s(a*b*c*d)
      double precision, intent(in) :: u(a*b*c*d)
      double precision, intent(in) :: f
      integer :: version
      version = 1000*i+100*j+10*k+l
      print*,'BEFORE',phase
      if (version.eq.1234) then
        call transpose_1234_y(u,s,a,b,c,d,f,async)
      elseif (version.eq.1243) then
        call transpose_1243_y(u,s,a,b,c,d,f,async)
      elseif (version.eq.1324) then
        call transpose_1324_y(u,s,a,b,c,d,f,async)
      elseif (version.eq.1342) then
        call transpose_1342_y(u,s,a,b,c,d,f,async)
      elseif (version.eq.1423) then
        call transpose_1423_y(u,s,a,b,c,d,f,async)
      elseif (version.eq.1432) then
        call transpose_1432_y(u,s,a,b,c,d,f,async)
      elseif (version.eq.2134) then
        call transpose_2134_y(u,s,a,b,c,d,f,async)
      elseif (version.eq.2143) then
        call transpose_2143_y(u,s,a,b,c,d,f,async)
      elseif (version.eq.2314) then
        call transpose_2314_y(u,s,a,b,c,d,f,async)
      elseif (version.eq.2341) then
        call transpose_2341_y(u,s,a,b,c,d,f,async)
      elseif (version.eq.2413) then
        call transpose_2413_y(u,s,a,b,c,d,f,async)
      elseif (version.eq.2431) then
        call transpose_2431_y(u,s,a,b,c,d,f,async)
      elseif (version.eq.3124) then
        call transpose_3124_y(u,s,a,b,c,d,f,async)
      elseif (version.eq.3142) then
        call transpose_3142_y(u,s,a,b,c,d,f,async)
      elseif (version.eq.3214) then
        call transpose_3214_y(u,s,a,b,c,d,f,async)
      elseif (version.eq.3241) then
        call transpose_3241_y(u,s,a,b,c,d,f,async)
      elseif (version.eq.3412) then
        call transpose_3412_y(u,s,a,b,c,d,f,async)
      elseif (version.eq.3421) then
        call transpose_3421_y(u,s,a,b,c,d,f,async)
      elseif (version.eq.4123) then
        call transpose_4123_y(u,s,a,b,c,d,f,async)
      elseif (version.eq.4132) then
        call transpose_4132_y(u,s,a,b,c,d,f,async)
      elseif (version.eq.4213) then
        call transpose_4213_y(u,s,a,b,c,d,f,async)
      elseif (version.eq.4231) then
        call transpose_4231_y(u,s,a,b,c,d,f,async)
      elseif (version.eq.4312) then
        call transpose_4312_y(u,s,a,b,c,d,f,async)
      elseif (version.eq.4321) then
        call transpose_4321_y(u,s,a,b,c,d,f,async)
      else
        print*,'something is wrong...'
      endif
      print*,'AFTER',phase
      end
        subroutine transpose_1234_y(u,s,d1,d2,d3,d4,f,async)
        implicit none
        integer, intent(in) :: d1,d2,d3,d4,async
        integer :: old_offset,new_offset
        integer :: j1,j2,j3,j4
        double precision, intent(out) :: s(d1*d2*d3*d4)
        double precision, intent(in) :: u(d1*d2*d3*d4)
        double precision, intent(in) :: f
        !!$acc parallel loop async(async)
        !do concurrent(j2 = 1:d2, j1 = 1:d1, j3 = 1:d3, j4 = 1:d4)
        !    old_offset = j4+d4*(j3-1+d3*(j2-1+d2*(j1-1)))
        !    new_offset = j4+d4*(j3-1+d3*(j2-1+d2*(j1-1)))
        !    s(new_offset) = u(old_offset) * f
        !enddo
        s = u * f
        print*,'transpose_1234_y should never be called'
        stop 1234
        end
        subroutine transpose_1243_y(u,s,d1,d2,d3,d4,f,async)
        implicit none
        integer, intent(in) :: d1,d2,d3,d4,async
        integer :: old_offset,new_offset
        integer :: j1,j2,j3,j4
        double precision, intent(out) :: s(d1*d2*d3*d4)
        double precision, intent(in) :: u(d1*d2*d3*d4)
        double precision, intent(in) :: f
        !$acc parallel loop async(async)
        do concurrent(j1 = 1:d1, j2 = 1:d2, j3 = 1:d3, j4 = 1:d4)
            !s(j3+d3*(j4-1+d4*(j2-1+d2*(j1-1)))) = f*u(j4+d4*(j3-1+d3*(j2-1+d2*(j1-1))))
            old_offset = j4+d4*(j3-1+d3*(j2-1+d2*(j1-1)))
            new_offset = j3+d3*(j4-1+d4*(j2-1+d2*(j1-1)))
            s(new_offset) = u(old_offset) * f
        enddo
        end
        subroutine transpose_1324_y(u,s,d1,d2,d3,d4,f,async)
        implicit none
        integer, intent(in) :: d1,d2,d3,d4,async
        integer :: old_offset,new_offset
        integer :: j1,j2,j3,j4
        double precision, intent(out) :: s(d1*d2*d3*d4)
        double precision, intent(in) :: u(d1*d2*d3*d4)
        double precision, intent(in) :: f
        !$acc parallel loop async(async)
        do concurrent(j1 = 1:d1, j3 = 1:d3, j2 = 1:d2, j4 = 1:d4)
            old_offset = j4+d4*(j3-1+d3*(j2-1+d2*(j1-1)))
            new_offset = j4+d4*(j2-1+d2*(j3-1+d3*(j1-1)))
            s(new_offset) = u(old_offset) * f
        enddo
        end
        subroutine transpose_1342_y(u,s,d1,d2,d3,d4,f,async)
        implicit none
        integer, intent(in) :: d1,d2,d3,d4,async
        integer :: old_offset,new_offset
        integer :: j1,j2,j3,j4
        double precision, intent(out) :: s(d1*d2*d3*d4)
        double precision, intent(in) :: u(d1*d2*d3*d4)
        double precision, intent(in) :: f
        !$acc parallel loop async(async)
        do concurrent(j1 = 1:d1, j3 = 1:d3, j2 = 1:d2, j4 = 1:d4)
            old_offset = j4+d4*(j3-1+d3*(j2-1+d2*(j1-1)))
            new_offset = j2+d2*(j4-1+d4*(j3-1+d3*(j1-1)))
            s(new_offset) = u(old_offset) * f
        enddo
        end
        subroutine transpose_1423_y(u,s,d1,d2,d3,d4,f,async)
        implicit none
        integer, intent(in) :: d1,d2,d3,d4,async
        integer :: old_offset,new_offset
        integer :: j1,j2,j3,j4
        double precision, intent(out) :: s(d1*d2*d3*d4)
        double precision, intent(in) :: u(d1*d2*d3*d4)
        double precision, intent(in) :: f
        !$acc parallel loop async(async)
        do concurrent(j1 = 1:d1, j2 = 1:d2, j4 = 1:d4, j3 = 1:d3)
            old_offset = j4+d4*(j3-1+d3*(j2-1+d2*(j1-1)))
            new_offset = j3+d3*(j2-1+d2*(j4-1+d4*(j1-1)))
            s(new_offset) = u(old_offset) * f
        enddo
        end
        subroutine transpose_1432_y(u,s,d1,d2,d3,d4,f,async)
        implicit none
        integer, intent(in) :: d1,d2,d3,d4,async
        integer :: old_offset,new_offset
        integer :: j1,j2,j3,j4
        double precision, intent(out) :: s(d1*d2*d3*d4)
        double precision, intent(in) :: u(d1*d2*d3*d4)
        double precision, intent(in) :: f
        !$acc parallel loop async(async)
        do concurrent(j1 = 1:d1, j3 = 1:d3, j4 = 1:d4, j2 = 1:d2)
            old_offset = j4+d4*(j3-1+d3*(j2-1+d2*(j1-1)))
            new_offset = j2+d2*(j3-1+d3*(j4-1+d4*(j1-1)))
            s(new_offset) = u(old_offset) * f
        enddo
        end
        subroutine transpose_2134_y(u,s,d1,d2,d3,d4,f,async)
        implicit none
        integer, intent(in) :: d1,d2,d3,d4,async
        integer :: old_offset,new_offset
        integer :: j1,j2,j3,j4
        double precision, intent(out) :: s(d1*d2*d3*d4)
        double precision, intent(in) :: u(d1*d2*d3*d4)
        double precision, intent(in) :: f
        !$acc parallel loop async(async)
        do concurrent(j1 = 1:d1, j2 = 1:d2, j3 = 1:d3, j4 = 1:d4)
            old_offset = j4+d4*(j3-1+d3*(j2-1+d2*(j1-1)))
            new_offset = j4+d4*(j3-1+d3*(j1-1+d1*(j2-1)))
            s(new_offset) = u(old_offset) * f
        enddo
        end
        subroutine transpose_2143_y(u,s,d1,d2,d3,d4,f,async)
        implicit none
        integer, intent(in) :: d1,d2,d3,d4,async
        integer :: old_offset,new_offset
        integer :: j1,j2,j3,j4
        double precision, intent(out) :: s(d1*d2*d3*d4)
        double precision, intent(in) :: u(d1*d2*d3*d4)
        double precision, intent(in) :: f
        !$acc parallel loop async(async)
        do concurrent(j1 = 1:d1, j2 = 1:d2, j3 = 1:d3, j4 = 1:d4)
            old_offset = j4+d4*(j3-1+d3*(j2-1+d2*(j1-1)))
            new_offset = j3+d3*(j4-1+d4*(j1-1+d1*(j2-1)))
            s(new_offset) = u(old_offset) * f
        enddo
        end
        subroutine transpose_2314_y(u,s,d1,d2,d3,d4,f,async)
        implicit none
        integer, intent(in) :: d1,d2,d3,d4,async
        integer :: old_offset,new_offset
        integer :: j1,j2,j3,j4
        double precision, intent(out) :: s(d1*d2*d3*d4)
        double precision, intent(in) :: u(d1*d2*d3*d4)
        double precision, intent(in) :: f
        !$acc parallel loop async(async)
        do concurrent(j2 = 1:d2, j1 = 1:d1, j3 = 1:d3, j4 = 1:d4)
            old_offset = j4+d4*(j3-1+d3*(j2-1+d2*(j1-1)))
            new_offset = j4+d4*(j1-1+d1*(j3-1+d3*(j2-1)))
            s(new_offset) = u(old_offset) * f
        enddo
        end
        subroutine transpose_2341_y(u,s,d1,d2,d3,d4,f,async)
        implicit none
        integer, intent(in) :: d1,d2,d3,d4,async
        integer :: old_offset,new_offset
        integer :: j1,j2,j3,j4
        double precision, intent(out) :: s(d1*d2*d3*d4)
        double precision, intent(in) :: u(d1*d2*d3*d4)
        double precision, intent(in) :: f
        !$acc parallel loop async(async)
        do concurrent(j2 = 1:d2, j3 = 1:d3, j1 = 1:d1, j4 = 1:d4)
            old_offset = j4+d4*(j3-1+d3*(j2-1+d2*(j1-1)))
            new_offset = j1+d1*(j4-1+d4*(j3-1+d3*(j2-1)))
            s(new_offset) = u(old_offset) * f
        enddo
        end
        subroutine transpose_2413_y(u,s,d1,d2,d3,d4,f,async)
        implicit none
        integer, intent(in) :: d1,d2,d3,d4,async
        integer :: old_offset,new_offset
        integer :: j1,j2,j3,j4
        double precision, intent(out) :: s(d1*d2*d3*d4)
        double precision, intent(in) :: u(d1*d2*d3*d4)
        double precision, intent(in) :: f
        !$acc parallel loop async(async)
        do concurrent(j2 = 1:d2, j1 = 1:d1, j4 = 1:d4, j3 = 1:d3)
            old_offset = j4+d4*(j3-1+d3*(j2-1+d2*(j1-1)))
            new_offset = j3+d3*(j1-1+d1*(j4-1+d4*(j2-1)))
            s(new_offset) = u(old_offset) * f
        enddo
        end
        subroutine transpose_2431_y(u,s,d1,d2,d3,d4,f,async)
        implicit none
        integer, intent(in) :: d1,d2,d3,d4,async
        integer :: old_offset,new_offset
        integer :: j1,j2,j3,j4
        double precision, intent(out) :: s(d1*d2*d3*d4)
        double precision, intent(in) :: u(d1*d2*d3*d4)
        double precision, intent(in) :: f
        !$acc parallel loop async(async)
        do concurrent(j2 = 1:d2, j3 = 1:d3, j1 = 1:d1, j4 = 1:d4)
            old_offset = j4+d4*(j3-1+d3*(j2-1+d2*(j1-1)))
            new_offset = j1+d1*(j3-1+d3*(j4-1+d4*(j2-1)))
            s(new_offset) = u(old_offset) * f
        enddo
        end
        subroutine transpose_3124_y(u,s,d1,d2,d3,d4,f,async)
        implicit none
        integer, intent(in) :: d1,d2,d3,d4,async
        integer :: old_offset,new_offset
        integer :: j1,j2,j3,j4
        double precision, intent(out) :: s(d1*d2*d3*d4)
        double precision, intent(in) :: u(d1*d2*d3*d4)
        double precision, intent(in) :: f
        !$acc parallel loop async(async)
        do concurrent(j2 = 1:d2, j3 = 1:d3, j4 = 1:d4, j1 = 1:d1)
            old_offset = j4+d4*(j3-1+d3*(j2-1+d2*(j1-1)))
            new_offset = j4+d4*(j2-1+d2*(j1-1+d1*(j3-1)))
            s(new_offset) = u(old_offset) * f
        enddo
        end
        subroutine transpose_3142_y(u,s,d1,d2,d3,d4,f,async)
        implicit none
        integer, intent(in) :: d1,d2,d3,d4,async
        integer :: old_offset,new_offset
        integer :: j1,j2,j3,j4
        double precision, intent(out) :: s(d1*d2*d3*d4)
        double precision, intent(in) :: u(d1*d2*d3*d4)
        double precision, intent(in) :: f
        !$acc parallel loop async(async)
        do concurrent(j2 = 1:d2, j3 = 1:d3, j4 = 1:d4, j1 = 1:d1)
            old_offset = j4+d4*(j3-1+d3*(j2-1+d2*(j1-1)))
            new_offset = j2+d2*(j4-1+d4*(j1-1+d1*(j3-1)))
            s(new_offset) = u(old_offset) * f
        enddo
        end
        subroutine transpose_3214_y(u,s,d1,d2,d3,d4,f,async)
        implicit none
        integer, intent(in) :: d1,d2,d3,d4,async
        integer :: old_offset,new_offset
        integer :: j1,j2,j3,j4
        double precision, intent(out) :: s(d1*d2*d3*d4)
        double precision, intent(in) :: u(d1*d2*d3*d4)
        double precision, intent(in) :: f
        !$acc parallel loop async(async)
        do concurrent(j2 = 1:d2, j3 = 1:d3, j4 = 1:d4, j1 = 1:d1)
            old_offset = j4+d4*(j3-1+d3*(j2-1+d2*(j1-1)))
            new_offset = j4+d4*(j1-1+d1*(j2-1+d2*(j3-1)))
            s(new_offset) = u(old_offset) * f
        enddo
        end
        subroutine transpose_3241_y(u,s,d1,d2,d3,d4,f,async)
        implicit none
        integer, intent(in) :: d1,d2,d3,d4,async
        integer :: old_offset,new_offset
        integer :: j1,j2,j3,j4
        double precision, intent(out) :: s(d1*d2*d3*d4)
        double precision, intent(in) :: u(d1*d2*d3*d4)
        double precision, intent(in) :: f
        !$acc parallel loop async(async)
        do concurrent(j2 = 1:d2, j3 = 1:d3, j4 = 1:d4, j1 = 1:d1)
            old_offset = j4+d4*(j3-1+d3*(j2-1+d2*(j1-1)))
            new_offset = j1+d1*(j4-1+d4*(j2-1+d2*(j3-1)))
            s(new_offset) = u(old_offset) * f
        enddo
        end
        subroutine transpose_3412_y(u,s,d1,d2,d3,d4,f,async)
        implicit none
        integer, intent(in) :: d1,d2,d3,d4,async
        integer :: old_offset,new_offset
        integer :: j1,j2,j3,j4
        double precision, intent(out) :: s(d1*d2*d3*d4)
        double precision, intent(in) :: u(d1*d2*d3*d4)
        double precision, intent(in) :: f
        !$acc parallel loop async(async)
        do concurrent(j1 = 1:d1, j3 = 1:d3, j2 = 1:d2, j4 = 1:d4)
            old_offset = j4+d4*(j3-1+d3*(j2-1+d2*(j1-1)))
            new_offset = j2+d2*(j1-1+d1*(j4-1+d4*(j3-1)))
            s(new_offset) = u(old_offset) * f
        enddo
        end
        subroutine transpose_3421_y(u,s,d1,d2,d3,d4,f,async)
        implicit none
        integer, intent(in) :: d1,d2,d3,d4,async
        integer :: old_offset,new_offset
        integer :: j1,j2,j3,j4
        double precision, intent(out) :: s(d1*d2*d3*d4)
        double precision, intent(in) :: u(d1*d2*d3*d4)
        double precision, intent(in) :: f
        !$acc parallel loop async(async)
        do concurrent(j2 = 1:d2, j3 = 1:d3, j1 = 1:d1, j4 = 1:d4)
            old_offset = j4+d4*(j3-1+d3*(j2-1+d2*(j1-1)))
            new_offset = j1+d1*(j2-1+d2*(j4-1+d4*(j3-1)))
            s(new_offset) = u(old_offset) * f
        enddo
        end
        subroutine transpose_4123_y(u,s,d1,d2,d3,d4,f,async)
        implicit none
        integer, intent(in) :: d1,d2,d3,d4,async
        integer :: old_offset,new_offset
        integer :: j1,j2,j3,j4
        double precision, intent(out) :: s(d1*d2*d3*d4)
        double precision, intent(in) :: u(d1*d2*d3*d4)
        double precision, intent(in) :: f
        !$acc parallel loop async(async)
        do concurrent(j2 = 1:d2, j3 = 1:d3, j4 = 1:d4, j1 = 1:d1)
            old_offset = j4+d4*(j3-1+d3*(j2-1+d2*(j1-1)))
            new_offset = j3+d3*(j2-1+d2*(j1-1+d1*(j4-1)))
            s(new_offset) = u(old_offset) * f
        enddo
        end
        subroutine transpose_4132_y(u,s,d1,d2,d3,d4,f,async)
        implicit none
        integer, intent(in) :: d1,d2,d3,d4,async
        integer :: old_offset,new_offset
        integer :: j1,j2,j3,j4
        double precision, intent(out) :: s(d1*d2*d3*d4)
        double precision, intent(in) :: u(d1*d2*d3*d4)
        double precision, intent(in) :: f
        !$acc parallel loop async(async)
        do concurrent(j2 = 1:d2, j3 = 1:d3, j4 = 1:d4, j1 = 1:d1)
            old_offset = j4+d4*(j3-1+d3*(j2-1+d2*(j1-1)))
            new_offset = j2+d2*(j3-1+d3*(j1-1+d1*(j4-1)))
            s(new_offset) = u(old_offset) * f
        enddo
        end
        subroutine transpose_4213_y(u,s,d1,d2,d3,d4,f,async)
        implicit none
        integer, intent(in) :: d1,d2,d3,d4,async
        integer :: old_offset,new_offset
        integer :: j1,j2,j3,j4
        double precision, intent(out) :: s(d1*d2*d3*d4)
        double precision, intent(in) :: u(d1*d2*d3*d4)
        double precision, intent(in) :: f
        !$acc parallel loop async(async)
        do concurrent(j2 = 1:d2, j3 = 1:d3, j4 = 1:d4, j1 = 1:d1)
            old_offset = j4+d4*(j3-1+d3*(j2-1+d2*(j1-1)))
            new_offset = j3+d3*(j1-1+d1*(j2-1+d2*(j4-1)))
            s(new_offset) = u(old_offset) * f
        enddo
        end
        subroutine transpose_4231_y(u,s,d1,d2,d3,d4,f,async)
        implicit none
        integer, intent(in) :: d1,d2,d3,d4,async
        integer :: old_offset,new_offset
        integer :: j1,j2,j3,j4
        double precision, intent(out) :: s(d1*d2*d3*d4)
        double precision, intent(in) :: u(d1*d2*d3*d4)
        double precision, intent(in) :: f
        !$acc parallel loop async(async)
        do concurrent(j2 = 1:d2, j3 = 1:d3, j4 = 1:d4, j1 = 1:d1)
            old_offset = j4+d4*(j3-1+d3*(j2-1+d2*(j1-1)))
            new_offset = j1+d1*(j3-1+d3*(j2-1+d2*(j4-1)))
            s(new_offset) = u(old_offset) * f
        enddo
        end
        subroutine transpose_4312_y(u,s,d1,d2,d3,d4,f,async)
        implicit none
        integer, intent(in) :: d1,d2,d3,d4,async
        integer :: old_offset,new_offset
        integer :: j1,j2,j3,j4
        double precision, intent(out) :: s(d1*d2*d3*d4)
        double precision, intent(in) :: u(d1*d2*d3*d4)
        double precision, intent(in) :: f
        !$acc parallel loop async(async)
        do concurrent(j3 = 1:d3, j1 = 1:d1, j4 = 1:d4, j2 = 1:d2)
            old_offset = j4+d4*(j3-1+d3*(j2-1+d2*(j1-1)))
            new_offset = j2+d2*(j1-1+d1*(j3-1+d3*(j4-1)))
            s(new_offset) = u(old_offset) * f
        enddo
        end
        subroutine transpose_4321_y(u,s,d1,d2,d3,d4,f,async)
        implicit none
        integer, intent(in) :: d1,d2,d3,d4,async
        integer :: old_offset,new_offset
        integer :: j1,j2,j3,j4
        double precision, intent(out) :: s(d1*d2*d3*d4)
        double precision, intent(in) :: u(d1*d2*d3*d4)
        double precision, intent(in) :: f
        !$acc parallel loop async(async)
        do concurrent(j2 = 1:d2, j3 = 1:d3, j1 = 1:d1, j4 = 1:d4)
            old_offset = j4+d4*(j3-1+d3*(j2-1+d2*(j1-1)))
            new_offset = j1+d1*(j2-1+d2*(j3-1+d3*(j4-1)))
            s(new_offset) = u(old_offset) * f
        enddo
        end
