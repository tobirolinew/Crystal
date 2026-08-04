!***********************************************************************
      subroutine update_lcoors(lpo_ptr,ldrive)
!***********************************************************************
      use mod_catch, only : catch_error
      use mod_lattice, only : lpo
      implicit none
      type(lpo), pointer :: &
         lpo_ptr
      character(*) :: &
         ldrive
      !***
      select case(ldrive)
         case('normal')
            call update_lcoors_normal()
         case default
            call catch_error(&
                     err=.true.,&
                     msg='Unknown "ldrive".',&
                     proc='update_lcoors')
      end select
      !***
      contains
      !*****************************************************************
         subroutine update_lcoors_normal()
      !*****************************************************************
         use mod_global, only : &
               ncoors, start_wcoors, low_wcoors, &
               up_wcoors, step_wcoors, fix_wcoors
         implicit none
         integer :: &
            lcLi,&
            lcUi,&
            lci0,&
            i
         !***
         do i=1,ncoors
            associate(ini_wci => start_wcoors(i),&
                      wci     => lpo_ptr%wcoors(i),&
                      Li      => low_wcoors(i),&
                      Ui      => up_wcoors(i),&
                      si      => step_wcoors(i),&
                      fixi    => fix_wcoors(i),&
                      lci     => lpo_ptr%lcoors(i))
               if(fixi) then
                  lci=0
               else
                  lcLi=ceiling((Li-ini_wci)/si)
                  lcUi=floor((Ui-ini_wci)/si)
                  lci0=nint((wci-ini_wci)/si)
                  lci=max(lcLi,min(lci0,lcUi))
               endif
            end associate
         end do
         !***
         end
      !***
      end
