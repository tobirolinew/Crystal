   !********************************************************************
      function cpot_f2b_homo(carts,abounds) result(ener)
   !********************************************************************
      implicit none
      integer :: &
         abounds(:,:)
      complex(r64) :: &
         carts(:,:),&
         ener,&
         dim_ener
      complex(r64), allocatable :: &
         dim_carts(:,:)
      integer :: &
         m,n,Lm,Um,Ln,Un
      !***
      ener=0.d0
      associate(nmons=>size(abounds,1))
         do m=1,nmons
            Lm=abounds(m,1)
            Um=abounds(m,2)
            do n=m+1,nmons
               Ln=abounds(n,1)
               Un=abounds(n,2)
               dim_carts=reshape([carts(:,Lm:Um),&
                                  carts(:,Ln:Un)],[3,Un-Ln+Um-Lm+2])
               dim_ener=cpot_f2b(dim_carts)
               ener=ener+dim_ener
            end do
         end do
      end associate
      !***
      end
