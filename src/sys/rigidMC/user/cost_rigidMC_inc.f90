!***********************************************************************
      subroutine cost_rigidMC(queue)
!***********************************************************************
      use mod_catch, only : catch_error
      use mod_lattice, only : prior_queue
      implicit none
      type(prior_queue) :: &
         queue
      select case(cost_drive)
         case('f2b_homo')
            call cost_f2b_homo()
         case default
            call catch_error(&
                     err=.true.,&
                     msg='Unknown cost_drive "'//cost_drive//'".',&
                     proc='cost_rigidMC')
      end select
      contains
      !*****************************************************************
         subroutine cost_f2b_homo()
      !*****************************************************************
         use mod_lattice, only : lpo
         use mod_f2b, only : pot_f2b_homo
         implicit none
         type(lpo), pointer :: &
            lpo_ptr
         type(geo), pointer :: &
            geo_ptr
         real(r64) :: &
            ener
         integer :: &
            i
         !$omp parallel do default(private) &
         !$omp shared(queue,layout,control)
         do i=1,queue%length
            lpo_ptr=>queue%at(i)
            geo_ptr=>resolve_geo_uptr(lpo_ptr%geo_uptr)
            if(control%rad_relax) then
               call geo_ptr%optimize(rad_relax=.true.,&
                                     dump=.false.)
               call exchange_wcoors(lpo_ptr,geo_ptr,'geo->lpo')
               lpo_ptr%cost=geo_ptr%ener
            else
               ener=pot_f2b_homo(geo_ptr%carts,layout%abounds)
               geo_ptr%ener=ener
               lpo_ptr%cost=ener
            endif
            call geo_ptr%reparametrize('assert_clean')
            deallocate(geo_ptr%carts)
         end do
         !$omp end parallel do
         end
      end
