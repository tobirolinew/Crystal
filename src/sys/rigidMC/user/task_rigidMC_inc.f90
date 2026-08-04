!***********************************************************************
      subroutine task_rigidMC(queue)
!***********************************************************************
      use mod_catch, only : catch_error
      use mod_lattice, only : prior_queue
      implicit none
      type(prior_queue) :: &
         queue
      write(*,'(/6x,a/)') 'Calling task_rigidMC() ...'
      select case(task_drive)
         case('msopt')
            call task_msopt()
         case default
            call catch_error(&
                     err=.true.,&
                     msg='Unknown task_drive "'//task_drive//'".',&
                     proc='task_rigidMC')
      end select
      contains
      !*****************************************************************
         subroutine task_msopt()
      !*****************************************************************
         use mod_lattice, only : lpo
         implicit none
         type(prior_queue) :: &
            tmp_queue,&
            targ_queue
         type(lpo), pointer :: &
            lpo_ptr
         type(geo), pointer :: &
            geo_ptr
         integer :: &
            ndone,&
            next_pcent,&
            i,&
            ntargets
         !***
         associate(max_targets => control%max_targets,&
                   max_ener    => control%max_ener,&
                   length      => queue%length)
            if(length.eq.0) then
               write(*,'(12x,a)') &
                       'Remark: no points to optimize!'
               return
            endif
            do i=1,length
               lpo_ptr=>queue%at(i)
               call tmp_queue%push(lpo_ptr,-lpo_ptr%cost)
            end do
            do i=1,length
               if(targ_queue%length.ge.max_targets) exit
               if(tmp_queue%weight_top().lt.-max_ener) exit
               lpo_ptr=>tmp_queue%pop()
               call targ_queue%push(lpo_ptr,-lpo_ptr%cost)
            end do
            ntargets=targ_queue%length
            if(ntargets.eq.0) then
               write(*,'(12x,a)') &
                       'Remark: no points below max_ener to optimize!'
               return
            endif
            write(*,'(9x,a,i0)') &
                    'Number of lattice points to be optimized: ',&
                    ntargets
            write(*,'(/9x,a/)') 'Optimization progress:'
         end associate
         ndone=0
         next_pcent=10
         open(1,file='done.log')
         !$omp parallel do default(private) &
         !$omp shared(targ_queue,ntargets,ndone,control) &
         !$omp shared(next_pcent) schedule(dynamic,1)
         do i=1,ntargets
            lpo_ptr=>targ_queue%at(i)
            geo_ptr=>resolve_geo_uptr(lpo_ptr%geo_uptr)
            call geo_ptr%reparametrize('assert_clean')
            call geo_ptr%optimize()
            !$omp critical
            ndone=ndone+1
            write(1,*) ndone
            associate(pcent=>100*ndone/ntargets)
               if(pcent.ge.next_pcent) then
                  write(*,'(11x,i3,a)') next_pcent,'% done'
                  next_pcent=next_pcent+10
               endif
            end associate
            !$omp end critical
         end do
         !$omp end parallel do
         call tmp_queue%delete()
         do i=1,ntargets
            lpo_ptr=>targ_queue%heap(i)%obj
            geo_ptr=>resolve_geo_uptr(lpo_ptr%geo_uptr)
            call tmp_queue%push(lpo_ptr,-geo_ptr%rgmax)
         end do
         do i=1,ntargets
            targ_queue%heap(i)%obj=>tmp_queue%pop()
         end do
         call tmp_queue%delete()
         write(*,'(/9x,a/)') 'Registration progress:'
         ndone=0
         next_pcent=10
         !$omp parallel do default(private) &
         !$omp shared(targ_queue,ntargets,control,landscape) &
         !$omp shared(ndone,next_pcent) schedule(dynamic,1)
         do i=1,ntargets-1
            lpo_ptr=>targ_queue%heap(i)%obj
            geo_ptr=>resolve_geo_uptr(lpo_ptr%geo_uptr)
            call register(geo_ptr,last=.false.)
            lpo_ptr%geo_uptr=>null()
            !$omp critical
            ndone=ndone+1
            associate(pcent=>100*ndone/ntargets)
               if(pcent.ge.next_pcent) then
                  write(*,'(11x,i3,a)') next_pcent,'% done'
                  next_pcent=next_pcent+10
               endif
            end associate
            !$omp end critical
         end do
         !$omp end parallel do
         lpo_ptr=>targ_queue%heap(ntargets)%obj
         geo_ptr=>resolve_geo_uptr(lpo_ptr%geo_uptr)
         call register(geo_ptr,last=.true.)
         lpo_ptr%geo_uptr=>null()
         call targ_queue%delete()
         close(1)
         !***
         end
      end
