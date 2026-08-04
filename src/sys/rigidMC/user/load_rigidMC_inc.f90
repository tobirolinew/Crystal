!***********************************************************************
      subroutine load_rigidMC()
!***********************************************************************
      use mod_block_input, only : &
         buffers,&
         schema,&
         store,&
         extend_schema,&
         print_schema,&
         print_store,&
         schema2store,&
         get_key_list
      implicit none
      type(schema) :: &
         sc
      type(store) :: &
         st
      !***
      call report('header')
      call report('start')
      call report('schema')
      call init_schema()
      call report('store')
      call schema2store(sc,st)
      call report('show')
      call report('user')
      call init_user_drives()
      call report('adjust')
      call adjust_keys()
      call report('workspace')
      call init_workspace()
      !***
      contains
      !*****************************************************************
         subroutine init_schema()
      !*****************************************************************
         implicit none
         integer :: &
            i
         !***
         sc%filename=buffers(1)%filename
         !geometry
         call get_key_list(sc%filename,'Geometry','monomer',sc%keys)
         do i=1,size(sc%keys)
            call extend_schema(sc,section='Geometry',key='',&
                     valtype='string',unit_flag=.false.)
         end do
         call extend_schema(sc,section='Geometry',key='cluster',&
                            valtype='string',unit_flag=.false.)
         !customize
         call extend_schema(sc,section='Customize',key='custom_drive',&
                  valtype='string',unit_flag=.false.)
         call extend_schema(sc,section='Customize',key='ang_step',&
                  valtype='real',unit_flag=.true.)
         call extend_schema(sc,section='Customize',key='rad_step',&
                  valtype='string',unit_flag=.true.)
         !cost
         call extend_schema(sc,section='Cost',key='cost_label',&
                  valtype='string',unit_flag=.false.)
         call extend_schema(sc,section='Cost',key='cost_unit',&
                  valtype='string',unit_flag=.false.)
         call extend_schema(sc,section='Cost',key='cost_spec',&
                  valtype='string',unit_flag=.false.)
         call extend_schema(sc,section='Cost',key='cost_drive',&
                  valtype='string',unit_flag=.false.)
         !task
         call extend_schema(sc,section='Task',key='task_drive',&
                  valtype='string',unit_flag=.false.)
         call extend_schema(sc,section='Task',key='max_order',&
                  valtype='string',unit_flag=.false.)
         call extend_schema(sc,section='Task',key='stac_point_file',&
                  valtype='string',unit_flag=.false.)
         call extend_schema(sc,section='Task',key='max_targets',&
                  valtype='string',unit_flag=.false.)
         call extend_schema(sc,section='Task',key='max_ener',&
                  valtype='string',unit_flag=.false.)
         call extend_schema(sc,section='Task',key='skip_flat',&
                  valtype='integer',unit_flag=.false.)
         !***
         end
      !*****************************************************************
         subroutine init_user_drives()
      !*****************************************************************
         use mod_catch, only : catch_error
         implicit none
         integer :: &
            i
         !***
         do i=1,size(st%keys)
            associate(key => st%keys(i),&
                      val => st%vals(i))
               select case(trim(key))
                  case('custom_drive')
                     custom_drive=trim(val)
                  case('cost_drive')
                     cost_drive=trim(val)
                  case('task_drive')
                     task_drive=trim(val)
               end select
            end associate
         end do
         !***
         call catch_error(&
            err=.not.allocated(custom_drive),&
            msg='missing required key "custom_drive"',&
            proc='load_rigidMC::init_user_drives')
         call catch_error(&
            err=.not.allocated(cost_drive),&
            msg='missing required key "cost_drive"',&
            proc='load_rigidMC::init_user_drives')
         call catch_error(&
            err=.not.allocated(task_drive),&
            msg='missing required key "task_drive"',&
            proc='load_rigidMC::init_user_drives')
         !***
         end
      !*****************************************************************
         subroutine adjust_keys()
      !*****************************************************************
         use mod_block_input, only : store
         implicit none
         integer :: &
            i
         !*** 
         if(.not.allocated(st%keys)) return
         !***
         do i=1,size(st%keys)
            associate(key => st%keys(i),&
                      val => st%vals(i))
               select case(trim(key))
                  case('cost_label')
                     key='pot_label'
                  case('cost_unit')
                     key='pot_unit'
                  case('cost_spec')
                     key='pot_spec'
                  case('cost_drive')
                     key='pot_drive'
               end select
            end associate
         end do
         !*** 
         end
      !*****************************************************************
         subroutine init_workspace()
      !*****************************************************************
         implicit none
         call report('control')
         call control%init(st)
         call report('layout')
         call layout%init()
         call report('refframe')
         call refframe%init()
         call report('sentinel')
         call sentinel%init()
         call report('landscape')
         call landscape%init()
         end
      !*****************************************************************
         subroutine report(mode)
      !*****************************************************************
         use mod_catch, only : catch_error
         implicit none
         character(*) :: &
            mode
         !***
         select case(mode)
            case('header')
               write(*,'(/3x,a/)') &
                       'Units of quantities used throughout this file:'
               write(*,'(6x,a)') &
                       'cut      : cm-1'
               write(*,'(6x,a)') &
                       'contrac  : cm-1'
               write(*,'(6x,a)') &
                       'objfunc  : cm-1'
               write(*,'(6x,a)') &
                       'min_cost : cm-1'
               write(*,'(6x,a)') &
                       'energy   : cm-1'
               write(*,'(6x,a)') &
                       'XYZ      : A'
               write(*,'(6x,a)') &
                       'gradnorm : '//&
                       'mixed (ObjFunc in cm-1, coords in A or deg)'
               write(*,'(6x,a)') &
                       'eigval   : '//&
                       'mixed (ObjFunc in cm-1, coords in A or deg)'
            case('start')
               write(*,'(/3x,a/)') &
                       'Calling load_rigicMC() ...'
            case('schema')
               write(*,'(6x,a/)') &
                       'Calling init_schema() ...'
            case('store')
               write(*,'(6x,a)') &
                       'Calling schema2store() ...'
            case('show')
               call print_store(st,nspaces=6)
            case('user')
               write(*,'(/6x,a)') &
                       'Calling init_user_drives() ...'
            case('adjust')
               write(*,'(/6x,a)') &
                       'Calling adjust_keys() ...'
            case('workspace')
               write(*,'(/6x,a)') &
                       'Calling init_workspace() ...'
            case('control')
               write(*,'(/9x,a)') &
                       'Calling control%init() ...'
            case('layout')
               write(*,'(/9x,a)') &
                       'Calling layout%init() ...'
            case('refframe')
               write(*,'(/9x,a)') &
                       'Calling refframe%init() ...'
            case('sentinel')
               write(*,'(/9x,a)') &
                        'Calling sentinel%init() ...'
            case('landscape')
               write(*,'(/9x,a)') &
                       'Calling landscape%init() ...'
            case default
               call catch_error(&
                        err=.true.,&
                        msg='unknown report mode "'//mode//'".',&
                        proc='report')
         end select
         !***
         end
      !***
      end
