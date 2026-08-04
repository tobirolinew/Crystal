!***********************************************************************
      function pot(rcoors) result(ener)
!***********************************************************************
      use mod_f2b, only : pot_f2b_homo
      use mod_catch, only : catch_error
      implicit none
     !type(geo) :: &
     !   this
      type(geo) :: &
         snapshot
      real(r64) :: &
         rcoors(:)
      real(r64) :: &
         ener
      !***
      associate(pot_drive => control%pot_drive,&
                nrcoors   => landscape%nrcoors)
         snapshot%rcoors=this%rcoors
         snapshot%carts=this%carts
         this%rcoors=rcoors(:nrcoors)
         call this%rcoors2carts()
         select case(pot_drive)
            case('f2b_homo')
               ener=pot_f2b_homo(this%carts,layout%abounds)
            case default
               call catch_error(&
                       err=.true.,&
                       msg='Unknown pot_drive "'//'"'//&
                           pot_drive//'"',&
                       proc='pot')
         end select
         this%rcoors=snapshot%rcoors
         this%carts=snapshot%carts
      end associate
      !***
      end
