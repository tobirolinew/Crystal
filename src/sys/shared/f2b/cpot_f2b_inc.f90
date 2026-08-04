   !********************************************************************
      function cpot_f2b(carts) result(ener)
   !********************************************************************
      implicit none
      complex(r64) :: &
         carts(:,:),&
         ener,&
         site_carts(3,nbof_sites)
      real(r64) :: &
         sing_tol=1.d-6
      !***
      site_carts=0.d0
      site_carts(:,:nbof_atoms)=carts(:,:nbof_atoms)
      if(nbof_sites.gt.nbof_atoms) &
         call add_OA_sites(site_carts,sing_tol,init=.false.)
      ener=kcalmolcm*f2b_func(site_carts)
      !***
      end
