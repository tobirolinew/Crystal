   !********************************************************************
      function pot_f2b_homo(carts,abounds) result(ener)
   !********************************************************************
      implicit none
      integer :: &
         abounds(:,:)
      real(r64) :: &
         carts(:,:),&
         ener
      complex(r64), allocatable :: &
         ccarts(:,:)
      complex(r64) :: &
         cener
      !***
      ccarts=carts
      cener=cpot_f2b_homo(ccarts,abounds)
      ener=dreal(cener)
      !***
      end
