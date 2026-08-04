!***********************************************************************
      subroutine customize_rigidMC(queue)
!***********************************************************************
      use mod_catch, only : catch_error
      use mod_lattice, only : prior_queue
      implicit none
      type(prior_queue) :: &
         queue
      !***
      select case(custom_drive)
         case('normal')
            call customize_normal()
         case default
            call catch_error(&
                     err=.true.,&
                     msg='Unknown custom_drive "'//custom_drive//'."',&
                     proc='customize_rigidMC')
      end select
      !***
      contains
      !*****************************************************************
         subroutine customize_normal()
      !*****************************************************************
         use mod_f2b, only : pot_f2b_homo
         use mod_lattice, only : lpo
         implicit none
         type(lpo), pointer :: &
            lpo_ptr
         type(geo), pointer :: &
            geo_ptr
         logical :: &
            abnorm
         integer :: &
            i
         !***
         !$omp parallel do default(private) &
         !$omp shared(queue,control,layout,refframe,landscape,sentinel)
         do i=1,queue%length
            lpo_ptr=>queue%at(i)
            call update_wcoors(lpo_ptr,wdrive='normal')
            allocate(geo_ptr)
            call geo_ptr%init()
            call exchange_wcoors(lpo_ptr,geo_ptr,'lpo->geo')
            call geo_ptr%renormalize()
            call geo_ptr%freeze()
            call exchange_wcoors(lpo_ptr,geo_ptr,'geo->lpo')
            call update_lcoors(lpo_ptr,ldrive='normal')
            call update_wcoors(lpo_ptr,wdrive='normal')
            call exchange_wcoors(lpo_ptr,geo_ptr,'lpo->geo')
            call geo_ptr%rcoors2carts()
            associate(rad_relax => control%rad_relax,&
                      parent    => lpo_ptr%parent,&
                      lcoors    => lpo_ptr%lcoors)
               if(rad_relax) then
                  abnorm=.false.
               else
                  abnorm=geo_ptr%get_abnorm_flag().ne.0
               endif
               if(abnorm) then
                  deallocate(geo_ptr)
                  lpo_ptr%geo_uptr=>null()
               else
                  lpo_ptr%geo_uptr=>geo_ptr
               endif
            end associate
         end do
         !$omp end parallel do
         !***
         end
      !***
      end
