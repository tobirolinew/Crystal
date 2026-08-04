!***********************************************************************
      function grad(rcoors,fixed) result(g)
!***********************************************************************
      use mod_f2b, only : grad_f2b_homo
      use mod_catch, only : catch_error
      implicit none
      type(geo) :: &
         snapshot
      real(r64) :: &
         rcoors(:)
      logical :: &
         fixed(:)
      real(r64), allocatable :: &
         g(:),&
         jac(:,:,:)
      !***
      associate(pot_drive => control%pot_drive,&
                nrcoors   => landscape%nrcoors,&
                abounds   => layout%abounds)
         snapshot%rcoors=this%rcoors
         snapshot%carts=this%carts
         this%rcoors=rcoors(:nrcoors)
         call this%rcoors2carts()
         jac=this%get_jac_carts(fixed)
         select case(pot_drive)
            case('f2b_homo')
               call grad_f2b_homo(&
                        this%carts,jac,layout%abounds,fixed,g)
            case default
               call catch_error(&
                       err=.true.,&
                       msg='Unknown pot_drive "'//'"'//&
                           pot_drive//'"',&
                       proc='grad')
         end select
         this%rcoors=snapshot%rcoors
         this%carts=snapshot%carts
      end associate
      !***
      end
