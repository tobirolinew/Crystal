!***********************************************************************
      module mod_crystal
!***********************************************************************
      use iso_fortran_env, only : &
         i64=>int64, r64=>real64
      use mod_lattice, only : &
         contrac, prior_queue, lock_set_, cost_set_
      use mod_global, only : &
         cut
      use omp_lib
      implicit none
      public :: &
         crystal,&
         cut,&
         branch_queue,&
         head_queue
      type(prior_queue) :: &
         seed_queue,&
         branch_queue,&
         head_queue,&
         lock_queue,&
         drop_queue
      type(lock_set_) :: &
         lock_set
      type(cost_set_) :: &
         cost_set
      integer :: &
         nthreads,&
         iter,&
         seed_size
      real(r64) :: &
         min_cost,&
         seed_min_cost
      real(r64) :: &
         propag_time=0.d0,&
         score_time=0.d0,&
         recast_time=0.d0,&
         postproc_time=0.d0
      contains
      !*****************************************************************
         subroutine crystal()
      !*****************************************************************
         implicit none
         !***
         call report('start')
         call configure()
         call report('header')
         iter=0
         do
            call propagate()
            call report('iter')
            if(branch_queue%is_empty()) exit
            call score()
            call recast()
            iter=iter+1
         end do
         call lock_set%reset()
         call report('postprocess')
         call postprocess()
         call report('final')
         !***
         end
      !*****************************************************************
         subroutine configure()
      !*****************************************************************
         use mod_sys, only : init_sys, load, globalize
         use mod_global, only : validate_globals
         implicit none
         !***
         call prepare()
         call init_sys()
         call load()
         call globalize()
         call validate_globals()
         !***
         contains
         !**************************************************************
            subroutine prepare()
         !**************************************************************
            use mod_catch, only : catch_error
            use mod_proj_plot, only : plot
            use mod_sys, only : sys_name
            use mod_block_input, only : &
                  get_clarg_filename, read_input, parse_keyval
            implicit none
            character(:), allocatable :: &
               inp_file
            !***
            inp_file=get_clarg_filename(iarg=1)
            call read_input(inp_file)
            call parse_keyval(inp_file,'Generic','cut',rval=cut)
            call parse_keyval(inp_file,'Generic','contrac',rval=contrac)
            call parse_keyval(&
                     inp_file,'Generic','sys_name',val=sys_name)
            call parse_keyval(&
                     inp_file,'Generic','nthreads',ival=nthreads)
            call parse_keyval(inp_file,'Generic','plot',ival=plot)
            call omp_set_num_threads(nthreads)
            call report('prepare')
            !***
            end
         !***
         end
      !*****************************************************************
         subroutine propagate()
      !*****************************************************************
         use mod_sys, only : customize
         implicit none
         !***
         associate(propag_time0=>omp_get_wtime())
            call collect()
            call customize(branch_queue)
            call prune()
            propag_time=propag_time-propag_time0+omp_get_wtime()
         end associate
         !***
         contains
         !**************************************************************
            subroutine collect()
         !**************************************************************
            use mod_global, only : ncoors
            use mod_lattice, only : lpo, instantiate
            implicit none
            type(lpo), pointer :: &
               lpo_ptr,&
               lpo_kid_ptr
            integer &
               i
            !***
            if(lock_set%n_keys_stored.eq.0) then
               min_cost=huge(0.d0)
               associate(lcoors=>spread(0,1,ncoors))
                  lpo_ptr=>instantiate(lcoors)
               end associate
               call branch_queue%push(lpo_ptr,weight=0.d0)
            else
               seed_min_cost=huge(0.d0)
               seed_size=seed_queue%length
               do while(.not.seed_queue%is_empty())
                  lpo_ptr=>seed_queue%pop()
                  seed_min_cost=min(seed_min_cost,lpo_ptr%cost)
                  associate(lcoors=>lpo_ptr%lcoors)
                     do i=-ncoors,ncoors
                        if(i.eq.0) cycle
                        lpo_kid_ptr=>instantiate(lcoors,i,lpo_ptr)
                        if(.not.associated(lpo_kid_ptr)) cycle
                        call branch_queue%push(lpo_kid_ptr,weight=0.d0)
                     end do
                  end associate
               end do
            endif
            !***
            end
         !**************************************************************
            subroutine prune()
         !**************************************************************
            use mod_catch, only : catch_error
            use mod_lattice, only : lpo
            implicit none
            logical :: &
               remove
            type(lpo), pointer :: &
               lpo_ptr
            integer(i64) :: &
               idx,stat
            integer :: &
               i,none
            !***
            i=1
            do while(i.le.branch_queue%length)
               lpo_ptr=>branch_queue%at(i)
               remove=.true.
               if(associated(lpo_ptr%geo_uptr)) then
                  idx=lock_set%get_index(lpo_ptr)
                  remove=idx.ne.-1
               endif
               if(remove) then
                  call branch_queue%erase(i)
                  deallocate(lpo_ptr%geo_uptr,stat=none)
                  deallocate(lpo_ptr)
               else
                  i=i+1
                  call lock_set%store_key(lpo_ptr,stat)
                  call catch_error(&
                           err=stat.eq.-1,&
                           msg='failed to extend "lock_set".',&
                           proc='pr')
                  call lock_queue%push(lpo_ptr,0.d0)
               endif
            end do
            !***
            end            
         end
      !*****************************************************************
         subroutine score()
      !*****************************************************************
         use mod_sys, only : cost
         implicit none
         !***
         associate(score_time0=>omp_get_wtime())
            call cost(branch_queue)
            score_time=score_time-score_time0+omp_get_wtime()
         end associate
         !***
         end
      !*****************************************************************
         subroutine recast()
      !*****************************************************************
         use mod_catch, only: catch_error
         use mod_lattice, only : lpo
         implicit none
         type(lpo), pointer :: &
            lpo_ptr
         integer(i64) :: &
            stat
         !***
         associate(recast_time0=>omp_get_wtime())
            do while(.not.branch_queue%is_empty())
               lpo_ptr=>branch_queue%pop()
               if(lpo_ptr%cost.gt.max(cut,min_cost)) cycle
               if(contrac.ne.0.d0) then
                  call cost_set%store_key(lpo_ptr%cost,stat,&
                           existing_key_is_error=.true.)
                  call catch_error(&
                           err=stat.eq.-1,&
                           msg='failed to extend "cost_set".',&
                           proc='recast')
                  if(stat.eq.-2) cycle
               endif
               call seed_queue%push(lpo_ptr,0.d0)
               if(lpo_ptr%cost.lt.cut) &
                  call head_queue%push(lpo_ptr,-lpo_ptr%cost)
            end do
            min_cost=-head_queue%weight_top()
            recast_time=recast_time-recast_time0+omp_get_wtime()
         end associate
         !***
         end
      !*****************************************************************
         subroutine postprocess()
      !*****************************************************************
         use mod_sys, only : task
         use mod_proj_plot, only : create_proj_plots
         implicit none
         !***
         associate(postproc_time0=>omp_get_wtime())
            call create_proj_plots(head_queue)
            call task(head_queue)
            postproc_time=postproc_time-postproc_time0+omp_get_wtime()
         end associate
         !***
         end
      !*****************************************************************
         subroutine report(mode)
      !*****************************************************************
         use mod_sys, only : sys_name
         use mod_catch, only : catch_error
         implicit none
         character(*) :: &
            mode
         select case(mode)
            case('start')
               write(*,'(/3x,a)') &
            "##########################################################"
                write(*,'(3x,a)') &
            "#                                                        #"
                write(*,'(3x,a)') &
            "#   $$$$$$                          $$            $$     #"
                write(*,'(3x,a)') &
            "#  $$   $$                          $$            $$     #"
                write(*,'(3x,a)') &
            "#  $$      $$$$$$ $$    $$ $$$$$$ $$$$$$$ $$$$$$$ $$     #"
                write(*,'(3x,a)') &
            "#  $$      $$  $$ $$    $$ $$       $$         $$ $$     #"
                write(*,'(3x,a)') &
            "#  $$      $$     $$    $$ $$$$$$   $$    $$$$$$$ $$     #"
                write(*,'(3x,a)') &
            "#  $$   $$ $$     $$    $$     $$   $$    $$   $$ $$     #"
                write(*,'(3x,a)') &
            "#   $$$$$$ $$      $$$$$$$ $$$$$$   $$$$$ $$$$$$$ $$$$$  #"
                write(*,'(3x,a)') &
            "#                       $$                               #"
                write(*,'(3x,a)') &
            "#                 $$    $$                               #"
                write(*,'(3x,a)') &
            "#                 $$$$$$$$                               #"
                write(*,'(3x,a)') &
            "#                                                        #"
                write(*,'(3x,a/)') &
            "##########################################################"
            case('prepare')
               write(*,'(3x,a/)') 'Calling configure::prepare() ...'
               write(*,'(6x,a,es0.16)') &
                       'cut          : ',cut
               write(*,'(6x,a,es0.16)') &
                       'contrac      : ',contrac
               write(*,'(6x,a)') &
                       'sys_name     : '//sys_name
               write(*,'(6x,a,i0)') &
                       'nthreads     : ',nthreads
            case('header')
               write(*,'(3x,a/)') 'Executing the lattice traversal ...'
               write(*,'(5a11,a12,4x,a)') &
                       '#', 'seed', 'branch', 'head', &
                       'lock', 'min_cost', 'min_cost(seed)'
            case('iter')
               if(iter.eq.0.or.seed_size.eq.0) return
               write(*,'(5i11$)') &
                  iter,&
                  seed_size,&
                  branch_queue%length,&
                  head_queue%length,&
                  lock_set%n_keys_stored
               if(dabs(min_cost).lt.1.d4) then
                  write(*,'(f12.2$)') min_cost
               else
                  write(*,'(es12.2$)') min_cost
               endif
               if(dabs(seed_min_cost).lt.1.d4) then
                  write(*,'(f12.2)') seed_min_cost
               else
                  write(*,'(es12.2)') seed_min_cost
               endif
            case('postprocess')
               write(*,'(/3x,a)') 'Calling postprocess() ...'
            case('final')
               associate(sum_time=>propag_time+score_time+&
                                   recast_time+postproc_time)
                  write(*,'(/3x,a/)') 'Timing results [s]:'
                  write(*,'(6x,a,f0.1)') 'propapate   : ',propag_time
                  write(*,'(6x,a,f0.1)') 'score       : ',score_time
                  write(*,'(6x,a,f0.1)') 'recast      : ',recast_time
                  write(*,'(6x,a,f0.1)') 'postprocess : ',postproc_time
                  write(*,'(6x,a,f0.1)') 'total       : ',sum_time
                  write(*,*)
               end associate
            case default
               call catch_error(&
                        err=.true.,&
                        msg='unknown mode "'//trim(mode)//'".',&
                        proc='report')
         end select
         end
      end
