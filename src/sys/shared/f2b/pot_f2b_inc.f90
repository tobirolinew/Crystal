   !********************************************************************
      function pot_f2b(carts) result(ener)
   !********************************************************************
      implicit none
      real(r64) :: &
         carts(:,:),&
         ener
      complex(r64), allocatable :: &
         ccarts(:,:)
      complex(r64) :: &
         cener
      !***
      ccarts=carts
      cener=cpot_f2b(ccarts)
      ener=dreal(cener)
      !***
      end
