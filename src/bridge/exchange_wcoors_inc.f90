!***********************************************************************
      subroutine exchange_wcoors(lpo_ptr,geo_ptr,mode)
!***********************************************************************
      use mod_lattice, only : lpo
      use mod_catch, only : catch_error
      implicit none
      type(lpo), pointer :: &
           lpo_ptr
      type(GEO_TYPE), pointer :: &
           geo_ptr
      character(*) :: &
         mode
      !***
      select case(mode)
         case('lpo->'//GEO_TYPE_STR)
            geo_ptr%WCOORS=lpo_ptr%wcoors
         case(GEO_TYPE_STR//'->lpo')
            lpo_ptr%wcoors=geo_ptr%WCOORS
         case default
            call catch_error(&
                     err=.true.,&
                     msg='Unknown exchange mode "'//trim(mode)//'".',&
                     proc='exchange_wcoors')
      end select
      !***
      end
#     undef GEO_TYPE
#     undef WCOORS
