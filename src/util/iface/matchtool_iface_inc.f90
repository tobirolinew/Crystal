      abstract interface
         subroutine matchtool_iface(links,matches,weights,cost)
         use iso_fortran_env, only: r64=>real64
         implicit none
         integer :: &
            links(:,:)
         real(r64), optional :: &
            weights(:),&
            cost
         integer, allocatable :: &
            matches(:)
         end
      end interface
