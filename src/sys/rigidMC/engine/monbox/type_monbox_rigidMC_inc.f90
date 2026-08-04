!**********************************************************************
      type monbox
!**********************************************************************
      integer :: &
         m                = 0,&
         ibase            = 0
      type(geo), pointer :: &
         targ_geo         => null()
      integer, pointer :: &
         first            => null(),&
         last             => null(),&
         natoms           => null(),&
         bin              => null(),&
         nedof            => null()
      character(:), pointer :: &
         ref_xyz_file     => null()
      character(2), pointer :: &
         symbs(:)         => null()
      integer, pointer :: &
         anums(:)         => null()
      real(r64), pointer :: &
         ref_carts(:,:)   => null(),&
         carts(:,:)       => null(),&
         rcoors(:)        => null(),&
         low_rcoors(:)    => null(),&
         up_rcoors(:)     => null(),&
         step_rcoors(:)   => null(),&
         gtol_rcoors(:)   => null(),&
         alpha            => null(),&
         beta             => null(),&
         gamma            => null(),&
         phi              => null(),&
         theta            => null(),&
         rho              => null()
      character(:), pointer :: &
         lab_rcoors(:)    => null(),&
         unit_rcoors(:)   => null()
      type(pgroup), pointer :: &
         pgroup           => null()
      type(coorbox) :: &
          calpha,&
          cbeta,&
          cgamma,&
          crho,&
          ctheta,&
          cphi
      contains
         procedure :: &
            init=>init_monbox
      end type
      !end
