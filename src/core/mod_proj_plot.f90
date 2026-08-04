!***********************************************************************
      module mod_proj_plot
!***********************************************************************
      use iso_fortran_env, only : r64=>real64
      use mod_lattice, only : prior_queue
      use mod_string, only : expand_tilde
      implicit none
      private
      public :: create_proj_plots, &
                  plot
      type proj_plot
         integer, allocatable :: coor_inds(:)
         type(prior_queue) :: proj_queue
         contains
            procedure :: &
               set_proj_queue,&
               write_input_file
      end type
      type(proj_plot), allocatable :: &
         proj_plots(:)
      integer :: &
         plot
      character(:), allocatable :: &
         crys_dir,&
         cost_unit_id,&
         inp_dir,&
         cmd
      contains
#        define ARR_TYPE type(proj_plot)
#        define APPEND append_proj_plot
#        include "../util/dep/generic_append_inc.f90"
      !*****************************************************************
         subroutine create_proj_plots(queue)
      !*****************************************************************
         implicit none
         type(prior_queue) :: &
            queue
         integer :: &
            i
         !***
         if(plot.eq.0) return
         write(*,'(/6x,a/)') &
            'Creating input files for ProjectionPlotter ...'
         call init_proj_plots()
         do i=1,size(proj_plots)
            associate(pplot=>proj_plots(i))
               call pplot%set_proj_queue(queue)
               call pplot%write_input_file()
            end associate
         end do
         write(*,'(/6x,a/)') &
            'Calling ProjectionPlotter.py ...'
         call execute_command_line(cmd)
         !***
         end
      !*****************************************************************
         subroutine init_proj_plots()
      !*****************************************************************
         use mod_xinfo, only : get_xdir
         use mod_catch, only : catch_error
         use mod_global, only : &
            ncoors, fix_wcoors, cost_unit
         implicit none
         character(:), allocatable :: &
            py,&
            home
         type(proj_plot), pointer :: &
            proj_plot_ptr
         character(:), allocatable :: &
            out_dir,&
            proj_plotter
         logical :: &
            ex
         integer :: &
            i,j
         !***
         do i=1,ncoors
            if(fix_wcoors(i)) cycle
            do j=1,i
               if(fix_wcoors(j)) cycle
               allocate(proj_plot_ptr)
               if(i.eq.j) then
                  proj_plot_ptr%coor_inds=[i]
               else
                  proj_plot_ptr%coor_inds=[i,j]
               endif
               call append_proj_plot(&
                        arr=proj_plots,&
                        elem=proj_plot_ptr,&
                        last=proj_plot_ptr)
            end do
         end do
         cost_unit_id=cost_unit
         if(allocated(cost_unit)) then
            select case(cost_unit)
               case('cm-1')
                  cost_unit_id='1'
               case('eV')
                  cost_unit_id='2'
               case('Eh')
                  cost_unit_id='3'
               case('kcal*mol-1')
                  cost_unit_id='4'
               case('kJ*mol-1')
                  cost_unit_id='5'
            end select
         endif
         associate(base_dir=>'plot/')
            inp_dir=base_dir//'inp/'
            out_dir=base_dir//'out/'
            call execute_command_line('rm -rf '//base_dir)
            call execute_command_line('mkdir '//base_dir)
            call execute_command_line('mkdir '//inp_dir)
            call execute_command_line('mkdir '//out_dir)
            crys_dir=get_xdir()
            proj_plotter=crys_dir//'/src/core/ProjectionPlotter.py'
            inquire(file=expand_tilde(proj_plotter),exist=ex)
            call catch_error(err=.not.ex,&
                              msg='ProjectionPlotter.py not found.',&
                              proc='init_proj_plots')
            call get_environment_variable('HOME',home)
            py='/gpfs1/home/t/r/troy5/.conda/envs/crystal/bin/python'
            inquire(file=expand_tilde(py),exist=ex)
            cmd=py//' '//proj_plotter//' '//inp_dir//' '//out_dir
            call catch_error(err=.not.ex,&
                              msg='Python is not functional.',&
                              proc='init_proj_plots')
         end associate
         !***
         end
      !*****************************************************************
         subroutine set_proj_queue(this,queue)
      !*****************************************************************
         use mod_lattice, only : prior_queue, lpo
         implicit none
         class(proj_plot) :: &
            this
         type(prior_queue) :: &
            queue
         type(lpo), pointer :: &
            lpo_ptr,&
            prlpo_ptr
         integer :: &
            i,j
         !***
         associate(coor_inds  => this%coor_inds,&
                   proj_queue => this%proj_queue)
            do i=1,queue%length
               lpo_ptr=>queue%at(i)
               associate(key     => lpo_ptr%lcoors(coor_inds),&
                         prqlen  => proj_queue%length,&
                         cost    => lpo_ptr%cost)
                  do j=1,prqlen
                     prlpo_ptr=>proj_queue%at(j)
                     associate(prlpo_key  => prlpo_ptr%lcoors(coor_inds))
                     associate(prlpo_cost => prlpo_ptr%cost)
                        if(all(prlpo_key.eq.key)) then
                           if(cost.lt.prlpo_cost) &
                              call proj_queue%replace_obj_at(j,lpo_ptr)
                           exit
                        endif
                     end associate
                     end associate
                  end do
                  if(j.eq.prqlen+1) &
                     call proj_queue%push(lpo_ptr,0.d0)
               end associate
            end do
         end associate
         !***
         end
      !*****************************************************************
         subroutine write_input_file(this)
      !*****************************************************************
         use mod_lattice, only : prior_queue, lpo
         use mod_global, only : lab_wcoors, unit_wcoors, cost_label
         implicit none
         class(proj_plot) :: &
            this
         type(lpo), pointer :: &
            proj_lpo_ptr
         character(:), allocatable :: &
            inp_file
         integer :: &
            i
         !***
         associate(coor_inds  => this%coor_inds,&
                   proj_queue => this%proj_queue)
         associate(ncinds     => size(coor_inds))
            inp_file=inp_dir//'plot_'
            do i=1,ncinds
               associate(idx=>coor_inds(i))
                  inp_file=&
                     inp_file//&
                     trim(lab_wcoors(idx))//'_'//&
                     trim(unit_wcoors(idx))//'_'
               end associate
            end do
            inp_file=inp_file//cost_label//'_'//cost_unit_id//'.txt'
            write(*,'(9x,a)') inp_file
            open(1,file=expand_tilde(inp_file))
            do i=1,proj_queue%length
               proj_lpo_ptr=>proj_queue%at(i)
               associate(cost   => proj_lpo_ptr%cost,&
                         wcoors => proj_lpo_ptr%wcoors(coor_inds))
                  write(1,'(3es26.16)') cost,wcoors
               end associate
            end do
            close(1)
         end associate
         end associate
         !***
         end
       end
