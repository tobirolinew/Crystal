!**********************************************************************
      type layout_ 
!**********************************************************************
      integer :: &
         nmons,&
         natoms
      integer, allocatable :: &
         mon_natoms(:),&
         abounds(:,:),&
         mon_bins(:),&
         mon_inds(:)
      contains
         procedure :: init=>init_layout
      end type
      !end
