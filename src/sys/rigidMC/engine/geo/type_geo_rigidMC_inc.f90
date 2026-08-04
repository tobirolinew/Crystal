!**********************************************************************
      type geo
!**********************************************************************
      integer, allocatable :: &
         anchors(:)
      real(r64), allocatable :: &
         carts(:,:),&
         rcoors(:),&
         axes(:,:,:),&
         grad(:),&
         sort_dists(:),&
         sens_hess_eigvals(:),&
         hess_eigvals(:),&
         hess_eigvects(:,:)
      real(r64) :: &
         ener=huge(0.d0),&
         rgmax=huge(0.d0)
      logical :: &
         conv=.false.
      contains
         procedure :: &
            !init
            init=>init_geo,&
            !manip
            rcoors2carts,&
            carts2rcoors,&
            reorient,&
            freeze,&
            renormalize,&
            align,&
            !reparam
            reparametrize,&
            adapt_anchors,&
            adapt_axes,&
            !sync
            get_anchors,&
            get_collin_flag,&
            get_axes,&
            get_fix_flag,&
            get_fix_rcoors,&
            get_abnorm_flag,&
            get_flat_flag,&
            get_saddle_order,&
            get_pgroup_name,&
            get_grmsdLB,&
            get_sort_dists,&
            get_jac_carts,&
            get_gmat,&
            !calc
            calc_ener,&
            calc_grad,&
            calc_hess_eigmodes,&
            calc_sens_hess_eigvals,&
            !screen
            optimize,&
            read_carts,&
            generate,&
            print_summary
      end type
!     end
