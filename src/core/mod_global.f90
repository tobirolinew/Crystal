!***********************************************************************
      module mod_global
!***********************************************************************
      use iso_fortran_env, only : r64=>real64
      integer :: &
         ncoors=0
      real(r64) :: &
         cut=huge(0.d0)
      logical, allocatable :: &
         fix_wcoors(:)
      real(r64), allocatable :: &
         low_wcoors(:),&
         up_wcoors(:),&
         start_wcoors(:),&
         step_wcoors(:)
      character(:), allocatable :: &
         cost_label,&
         cost_unit,&
         lab_wcoors(:),&
         unit_wcoors(:)
      contains
      !*****************************************************************
         subroutine validate_globals()
      !*****************************************************************
         use mod_catch, only : catch_error
         implicit none
         call catch_error(err=ncoors.le.0,&
                  msg='"ncoors" must be greater than zero.',&
                  proc='validate_globals')
         call catch_error(err=.not.allocated(fix_wcoors),&
                  msg='Array "fix_wcoors" is not allocated.',&
                  proc='validate_globals')
         call catch_error(err=.not.allocated(low_wcoors),&
                  msg='Array "low_wcoors" is not allocated.',&
                  proc='validate_globals')
         call catch_error(err=.not.allocated(up_wcoors),&
                  msg='Array "up_wcoors" is not allocated.',&
                  proc='validate_globals')
         call catch_error(err=.not.allocated(start_wcoors),&
                  msg='Array "start_wcoors" is not allocated.',&
                  proc='validate_globals')
         call catch_error(err=.not.allocated(step_wcoors),&
                  msg='Array "step_wcoors" is not allocated.',&
                  proc='validate_globals')
         call catch_error(err=.not.allocated(lab_wcoors),&
                  msg='Array "lab_wcoors" is not allocated.',&
                  proc='validate_globals')
         call catch_error(err=.not.allocated(unit_wcoors),&
                  msg='Array "unit_wcoors" is not allocated.',&
                  proc='validate_globals')
         call catch_error(err=.not.allocated(cost_label),&
                  msg='Character variable "cost_label" '//&
                       'is not allocated.',&
                  proc='validate_globals')
         call catch_error(err=.not.allocated(cost_unit),&
                  msg='Character variable "cost_unit" '//&
                      'is not allocated.',&
                  proc='validate_globals')
         end
      end
