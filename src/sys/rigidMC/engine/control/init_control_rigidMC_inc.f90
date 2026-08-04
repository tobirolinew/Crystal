!***********************************************************************
      subroutine init_control(this,st)
!***********************************************************************
      use mod_block_input, only : store
      use mod_catch, only : catch_error
      use mod_f2b, only : read_f2b_pars
      implicit none
      class(control_) :: &
            this
      type(store) :: &
         st
      integer :: &
         i,ios
      !***
      call catch_local_errors('arrays')
      associate(t=>this)
         t%ref_xyz_files=[character(0)::]
         t%start_xyz_file=''
         do i=1,size(st%keys)
            associate(key  => st%keys(i),&
                      val  => st%vals(i),&
                      ival => st%ivals(i),&
                      rval => st%rvals(i),&
                      unit => st%units(i))
               if(index(key,'monomer').eq.1) then
                  call append_string(arr=t%ref_xyz_files,elem=val)
               else
                  select case(trim(key))
                     case('cluster')
                        t%start_xyz_file=trim(val)
                     case('pot_label')
                        t%pot_label=trim(val)
                     case('pot_unit')
                        t%pot_unit=trim(val)
                     case('pot_spec')
                        t%pot_spec=trim(val)
                     case('pot_drive')
                        t%pot_drive=trim(val)
                     case('ang_step')
                        t%ang_step=rval
                        t%ang_unit=trim(unit)
                     case('rad_step')
                        if(trim(val).eq.'relax') then
                           t%rad_relax=.true.
                           t%rad_step=0.d0
                        else
                           t%rad_relax=.false.
                           read(val,*,iostat=ios) t%rad_step
                           call catch_local_errors('norelax')
                        endif
                        t%rad_unit=trim(unit)
                     case('skip_flat')
                        call catch_local_errors('skip_flat')
                        t%skip_flat=ival.eq.1
                     case('stac_point_file')
                        t%stac_point_file=trim(val)
                     case('max_targets')
                        t%max_targets=huge(0)
                        if(trim(val).ne.'none') then
                           read(val,*,iostat=ios) t%max_targets
                           call catch_local_errors('max_targets')
                        endif
                     case('max_order')
                        t%max_order=huge(0)
                        if(trim(val).ne.'none') then
                           read(val,*,iostat=ios) t%max_order
                           call catch_local_errors('max_order')
                        endif
                     case('max_ener')
                        t%max_ener=huge(0.d0)
                        if(trim(val).ne.'none') then
                           read(val,*,iostat=ios) t%max_ener
                           call catch_local_errors('max_ener')
                        endif
                  end select
               endif
            end associate
         end do
         select case(t%pot_drive)
            case('f2b_homo')
               call read_f2b_pars(t%pot_spec)
         end select
         call catch_local_errors('final')
      end associate
      !***
      contains
      !*****************************************************************
         subroutine catch_local_errors(opt)
      !*****************************************************************
         use mod_string, only : to_string
         implicit none
         character(*) :: &
            opt
         !***
         associate(t=>this)
            select case(opt)
               case('arrays')
                  call catch_error(&
                           err=.not.allocated(st%keys),&
                           msg='array "st%keys" is not allocated.',&
                           proc='control%init')
                  call catch_error(&
                           err=.not.allocated(st%vals),&
                           msg='array "st%vals" is not allocated.',&
                           proc='control%init')
                  call catch_error(&
                           err=.not.allocated(st%ivals),&
                           msg='array "st%ivals" is not allocated.',&
                           proc='control%init')
                  call catch_error(&
                           err=.not.allocated(st%rvals),&
                           msg='array "st%rvals" is not allocated.',&
                           proc='control%init')
               case('norelax')
                  call catch_error(&
                           err=ios.ne.0,&
                           msg='unable to extract a real value from '//&
                               '"rad_step".',&
                           proc='control%init')
               case('max_order')
                  call catch_error(&
                           err=ios.ne.0.or.t%max_order.lt.0,&
                           msg='invalid value detected for '//&
                               '"max_order".',&
                           proc='control%init',&
                           comment='Only nonnegative integers '//&
                                   'are allowed.')
               case('skip_flat')
                  call catch_error(&
                           err=all(st%ivals(i).ne.[0,1]),&
                           msg='unable to extract a binary (0/1) '//&
                               'integer value from "skip_flat".',&
                           proc='control%init')
               case('max_targets')
                  call catch_error(&
                           err=t%max_targets.lt.0,&
                           msg='unable to extract a nonnegative '//&
                               'integer value from "max_targets".',&
                           proc='control%init')
               case('max_ener')
                  call catch_error(&
                           err=ios.ne.0,&
                           msg='unable to extract a real value from '//&
                               '"max_ener".',&
                           proc='control%init')
               case('final')
                  call catch_error(&
                           err=.not.allocated(t%ref_xyz_files),&
                           msg='no monomer filenames found in input.',&
                           proc='control%init')
                  call catch_error(&
                           err=size(t%ref_xyz_files).lt.2,&
                           msg='insufficent of monomer filenames '//&
                               '(nmons = '//to_string(&
                               size(t%ref_xyz_files))//')'//&
                               'found in input.',&
                           proc='control%init')
                  call catch_error(&
                           err=.not.allocated(t%start_xyz_file),&
                           msg='missing required key "cluster".',&
                           proc='control%init')
                  call catch_error(&
                           err=.not.allocated(t%pot_label),&
                           msg='missing required key "pot_label".',&
                           proc='control%init')
                  call catch_error(&
                           err=.not.allocated(t%pot_unit),&
                           msg='missing required key "pot_unit"',&
                           proc='control%init')
                  call catch_error(&
                           err=.not.allocated(t%pot_spec),&
                           msg='missing required key "pot_spec"',&
                           proc='control%init')
                  call catch_error(&
                           err=.not.allocated(t%stac_point_file),&
                           msg='missing required key "stac_point_file"',&
                           proc='control%init')
                  call catch_error(&
                           err=t%ang_step.le.0.d0,&
                           msg='invalid or missing "ang_step"',&
                           proc='control%init')
                  call catch_error(&
                           err=.not.t%rad_relax.and.t%rad_step.le.0.d0,&
                           msg='invalid or missing "rad_step"',&
                           proc='control%init')
                  call catch_error(&
                           err=.not.allocated(t%ang_unit),&
                           msg='missing required key "ang_unit"',&
                           proc='control%init')
                  call catch_error(&
                           err=.not.allocated(t%rad_unit),&
                           msg='missing required key "rad_unit"',&
                           proc='control%init')
                  select case(t%ang_unit)
                     case('deg')
                     case default
                        call catch_error(&
                                 err=.true.,&
                                 msg='unknown "ang_unit"',&
                                 proc='control%init')
                  end select
                  select case(t%rad_unit)
                     case('A')
                     case default
                        call catch_error(&
                                 err=.true.,&
                                 msg='unknown "rad_unit"',&
                                 proc='control%init')
                  end select
                  select case(t%pot_unit)
                     case('cm-1')
                     case default
                        call catch_error(&
                                 err=.true.,&
                                 msg='unknown "pot_unit"',&
                                 proc='control%init')
                  end select
               case default
                  call catch_error(&
                           err=.true.,&
                           msg='unknown catch option',&
                           proc='control%init')
            end select
      end associate
         !***
         end
      !*****************************************************************
      !  subroutine append_string()
      !*****************************************************************
#        define STR_ARR_TYPE character(:)
#        define APPEND append_string
#        include "../../../../util/dep/generic_append_inc.f90"
         !end
      end
