!***********************************************************************
      module mod_sys
!***********************************************************************
      use iso_fortran_env
      use mod_lattice, only : prior_queue
      implicit none
      abstract interface
         !*** pload
         subroutine pload()
         end
         !*** pglobalize
         subroutine pglobalize()
         end
         !*** pcustomize
         subroutine pcustomize(queue)
         use mod_lattice, only : prior_queue
         type(prior_queue) :: queue
         end
         !*** pcost
         subroutine pcost(queue)
         use mod_lattice, only : prior_queue
         type(prior_queue) :: queue
         end
         !*** ptask
         subroutine ptask(queue)
         use mod_lattice, only : prior_queue
         type(prior_queue) :: queue
         end
      end interface
      procedure(pload), pointer :: &
         load
      procedure(pglobalize), pointer :: &
         globalize
      procedure(pcustomize), pointer :: &
         customize
      procedure(pcost), pointer :: &
         cost
      procedure(ptask), pointer :: &
         task
      character(:), allocatable :: &
         sys_name
      contains
      !*****************************************************************
         subroutine init_sys()
      !*****************************************************************
         use mod_catch, only : catch_error
         use mod_user_rigidMC, only : &
            load_rigidMC,&
            globalize_rigidMC,&
            customize_rigidMC,&
            cost_rigidMC,&
            task_rigidMC
         implicit none
         select case(sys_name)
            case('rigidMC')
               load=>load_rigidMC
               globalize=>globalize_rigidMC
               customize=>customize_rigidMC
               cost=>cost_rigidMC
               task=>task_rigidMC
            case default
               call catch_error(&
                        err=.true.,&
                        msg='Unknown "sys_name".',&
                        proc='init_sys')
         end select
         end
      end
