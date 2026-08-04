!***********************************************************************
      subroutine update_wcoors(lpo_ptr,wdrive)
!***********************************************************************
      use mod_catch, only : catch_error
      use mod_lattice, only : lpo
      implicit none
      type(lpo), pointer :: &
         lpo_ptr
      character(*) :: &
         wdrive
      !***
      select case(wdrive)
         case('normal')
            call update_wcoors_normal()
         case default
            call catch_error(&
                     err=.true.,&
                     msg='Unknown "wdrive".',&
                     proc='update_wcoors')
      end select
      !***
      contains
      !*****************************************************************
         subroutine update_wcoors_normal()
      !*****************************************************************
         use mod_global, only : &
            ncoors, start_wcoors, step_wcoors, fix_wcoors
         implicit none
         type(lpo), pointer :: &
            p
         integer :: &
            i
         !***
         lpo_ptr%wcoors=spread(0.d0,1,ncoors)
         p=>lpo_ptr%parent
         associate(w  => lpo_ptr%wcoors,&
                   L  => lpo_ptr%lcoors,&
                   w0 => start_wcoors,&
                   s  => step_wcoors)
            do i=1,ncoors
               if(fix_wcoors(i)) then
                  if(associated(p)) then
                     w(i)=p%wcoors(i)
                  else
                     w(i)=w0(i)
                  endif
               else
                  w(i)=w0(i)+L(i)*s(i)
               endif
            end do
         end associate
         !***
         end
      !***
      end
