!***********************************************************************
      module mod_sort
!***********************************************************************
      use iso_fortran_env, only: &
         r32=>real32,&
         r64=>real64,&
         i32=>int32,&
         i64=>int64
      implicit none
      contains
      !*****************************************************************
         function get_num_sort_perm(x,asc,g,ndigs,nlow) result(perm)
      !*****************************************************************
         use mod_catch, only : catch_error
         implicit none
         type lview_
            integer(i64), allocatable :: &
               g(:)
            integer :: &
               nlow
         end type
         type(lview_) :: lview
         class(*) :: &
            x(:)
         logical :: &
            asc
         class(*), optional :: &
            g(:)
         integer, optional :: &
            ndigs,&
            nlow
         integer, allocatable :: &
            perm(:)
         real(r64), allocatable :: &
            r(:)
         integer :: &
            n,i,j,pi,pj
         !***
         call init()
         do i=1,lview%nlow
            do j=i+1,n
               pi=perm(i)
               pj=perm(j)
               associate(gpi => lview%g(pi),&
                         gpj => lview%g(pj),&
                         rpi => r(pi),&
                         rpj => r(pj))
                  if(gpi.ne.gpj) cycle
                  if(asc.and.rpi.le.rpj) cycle
                  if(.not.asc.and.rpi.ge.rpj) cycle
                  perm(i)=pj
                  perm(j)=pi
               end associate
            end do
         end do
         !***
         contains
         !**************************************************************
            subroutine init()
         !**************************************************************
            implicit none
            real(r64) :: &
               fac
            !***
            !n
            n=size(x)
            !nlow
            lview%nlow=n
            if(present(nlow)) lview%nlow=nlow
            !g
            lview%g=spread(1_i64,1,n)
            if(present(g)) then
               select type(g)
               type is(integer(i32))
                  lview%g=g
               type is(integer(i64))
                  lview%g=g
               class default
                  call catch_error(&
                           err=.true.,&
                           msg='unknown g type.',&
                           proc='get_num_sort_perm')
               end select
            end if
            !r
            select type(x)
            type is(integer(i32))
               r=x
            type is(integer(i64))
               r=x
            type is(real(r32))
               r=x
            type is(real(r64))
               r=x
            class default
               call catch_error(&
                        err=.true.,&
                        msg='unknown array type.',&
                        proc='get_num_sort_perm')
            end select
            if(present(ndigs)) then
               fac=10**ndigs
               r=nint(fac*r)/fac
            endif
            !perm
            perm=[(i,i=1,n)]
            !***
            end
         !***
         end
      !*****************************************************************
         function &
               get_num_lex_sort_perm(x,asc,g,ndigs,nlow) result(perm)
      !*****************************************************************
         use mod_catch, only : catch_error
         implicit none
         type lview_
            integer(i64), allocatable :: &
               g(:)
            integer :: &
               nlow
         end type
         type(lview_) :: lview
         class(*) :: &
            x(..)
         logical :: &
            asc
         class(*), optional :: &
            g(:)
         integer, optional :: &
            ndigs,&
            nlow
         integer, allocatable :: &
            perm(:)
         real(r64), allocatable :: &
            rflat(:,:)
         integer :: &
            n,i,j,pi,pj
         !***
         call init()
         do i=1,lview%nlow
            do j=i+1,n
               pi=perm(i)
               pj=perm(j)
               associate(gpi => lview%g(pi),&
                         gpj => lview%g(pj))
                  if(gpi.ne.gpj) cycle
                  if(lex_less(pi,pj).eqv.asc) cycle
                  perm(i)=pj
                  perm(j)=pi
               end associate
            end do
         end do
         !***
         contains
         !**************************************************************
            subroutine init()
         !**************************************************************
            implicit none
            real(r64) :: &
               fac
            integer :: &
               d,m,p
            !***
            !size
            associate(s=>shape(x))
               d=size(s)
               n=s(d)
               m=product(s(:d-1))
               p=n*m
            end associate
            !rflat
            select rank(x)
            rank(1)
#              include "dep/sort_select_inc.f90"
            rank(2)
#              include "dep/sort_select_inc.f90"
            rank(3)
#              include "dep/sort_select_inc.f90"
            rank(4)
#              include "dep/sort_select_inc.f90"
            rank(5)
#              include "dep/sort_select_inc.f90"
            rank(6)
#              include "dep/sort_select_inc.f90"
            rank(7)
#              include "dep/sort_select_inc.f90"
            rank(8)
#              include "dep/sort_select_inc.f90"
            rank(9)
#              include "dep/sort_select_inc.f90"
            rank(10)
#              include "dep/sort_select_inc.f90"
            rank(11)
#              include "dep/sort_select_inc.f90"
            rank(12)
#              include "dep/sort_select_inc.f90"
            rank(13)
#              include "dep/sort_select_inc.f90"
            rank(14)
#              include "dep/sort_select_inc.f90"
            rank(15)
#              include "dep/sort_select_inc.f90"
            rank default
               call catch_error(&
                        err=.true.,&
                        msg='unsupported rank.',&
                        proc='get_num_lex_sort_perm')
            end select
            if(present(ndigs)) then
               fac=10**ndigs
               rflat=nint(fac*rflat)/fac
            endif
            !nlow
            lview%nlow=n
            if(present(nlow)) lview%nlow=nlow
            !g
            lview%g=spread(1_i64,1,n)
            if(present(g)) then
               select type(g)
               type is(integer(i32))
                  lview%g=g
               type is(integer(i64))
                  lview%g=g
               class default
                  call catch_error(&
                           err=.true.,&
                           msg='unknown g type.',&
                           proc='get_num_lex_sort_perm')
               end select
            end if
            !perm
            perm=[(i,i=1,n)]
            !***
            end
         !**************************************************************
            function lex_less(pi,pj) result(less)
         !**************************************************************
            implicit none
            integer :: &
               pi,pj,k
            logical :: &
               less
            !***
            less=.false.
            k=findloc(rflat(:,pi).ne.rflat(:,pj),value=.true.,dim=1)
            if(k.eq.0) return
            less=rflat(k,pi).lt.rflat(k,pj)
            !***
            end
         !***
         end
      !***
      end