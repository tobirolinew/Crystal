!**********************************************************************
      type landscape_
!**********************************************************************
      integer :: &
         nrcoors
      real(r64) :: &
         rho_ubound
      real(r64), allocatable :: &
         low_rcoors(:),&
         up_rcoors(:),&
         step_rcoors(:),&
         gtol_rcoors(:),&
         dvtol_rcoors(:),&
         scal_rcoors(:)
      character(:), allocatable :: &
         lab_rcoors(:),&
         unit_rcoors(:)
      type(geo), pointer :: &
         start_geo
      type(geo_prior_queue) :: &
         registry
      contains
         procedure :: init=>init_landscape
      end type
      !end
