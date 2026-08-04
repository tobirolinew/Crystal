!***********************************************************************
      module mod_lpo
!***********************************************************************
      use iso_fortran_env, only : r64=>real64
      implicit none
      type :: lpo
         real(r64) :: &
            cost
         integer, allocatable :: &
            lcoors(:)
         real(r64), allocatable :: &
            wcoors(:)
         type(lpo), pointer :: &
            parent
         class(*), pointer :: &
            geo_uptr
      end type
      end
