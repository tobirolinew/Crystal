   !********************************************************************
      subroutine grad_f2b_homo(carts,jac,abounds,fixed,g)
   !********************************************************************
      implicit none
      integer :: &
         abounds(:,:)
      real(r64) :: &
         carts(:,:),&
         jac(:,:,:)
      logical :: &
         fixed(:)
      real(r64), allocatable :: &
         g(:)
      complex(r64), allocatable :: &
         cartsc(:,:)
      integer :: &
         k
      !***
      associate(nv=>size(fixed))
         g=spread(0.d0,1,nv)
         do k=1,nv
            if(fixed(k)) cycle
            cartsc=dcmplx(carts,cfdstep*jac(:,:,k))
            g(k)=aimag(cpot_f2b_homo(cartsc,abounds))/cfdstep
         end do
      end associate
      !***
      end
