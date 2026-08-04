!**********************************************************************
      type coorbox
!**********************************************************************
      integer :: &
         idx    =  0
      real(r64), pointer :: &
         val    => null()
      real(r64), pointer :: &
         low    => null()
      real(r64), pointer :: &
         up     => null()
      real(r64), pointer :: &
         step   => null()
      real(r64), pointer :: &
         gtol   => null()
      character(:), pointer :: &
         lab    => null()
      character(:), pointer :: &
         unit   => null()
      contains
         procedure :: &
            init=>init_coorbox
      end type
      !end
