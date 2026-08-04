!***********************************************************************
      program main
!***********************************************************************
      use mod_stack_unlim, only : stack_unlim
      use mod_crystal
      implicit none
      call stack_unlim()
      call crystal()
      end
