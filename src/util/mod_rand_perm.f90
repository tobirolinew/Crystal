!***********************************************************************
      module mod_rand_perm
!***********************************************************************
      use mod_swap, only : swap
      implicit none
      contains
      !*****************************************************************
         function init_perm(n) result(p)
      !*****************************************************************
         implicit none
         integer, allocatable :: &
            p(:)
         integer :: &
            n,i
         !***
         p=[(i,i=1,n)]
         !***
         end
      !*****************************************************************
         subroutine rand_perm(p)
      !*****************************************************************
         implicit none
         integer :: &
            p(:)
         real :: &
            rnd
         integer :: &
            i,j
         !***
         do i=size(p),2,-1
            call random_number(rnd)
            j=1+int(rnd*i)
            call swap(p(i),p(j))
         end do
         !***
         end
      !***
      end
