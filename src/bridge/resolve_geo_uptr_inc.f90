!***********************************************************************
      function resolve_geo_uptr(geo_uptr) result(geo_ptr)
!***********************************************************************
      use mod_catch, only : catch_error
      implicit none
      class(*), pointer :: &
         geo_uptr
      type(GEO_TYPE), pointer :: &
         geo_ptr
      !***
      select type(geo_uptr)
         class is(GEO_TYPE)
            geo_ptr=>geo_uptr
         class default
            call catch_error(&
                     err=.true.,&
                     msg='Unknown GEO_TYPE.',&
                     proc='resolve_geo_uptr')
      end select
      !***
      end
