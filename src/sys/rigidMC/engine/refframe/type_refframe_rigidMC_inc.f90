!**********************************************************************
      type refframe_
!**********************************************************************
         character(2), allocatable :: &
            symbs(:)
         real(r64), allocatable :: &
            carts(:,:)
         integer, allocatable :: &
            anums(:),&
            nedofs(:)
         type(pgroup), allocatable :: &
            pgroups(:)
         integer :: &
            max_nedof
         contains
            procedure :: init=>init_refframe
      end type
      !end
