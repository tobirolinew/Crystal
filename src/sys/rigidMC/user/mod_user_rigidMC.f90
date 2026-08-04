!***********************************************************************
      module mod_user_rigidMC
!***********************************************************************
      use iso_fortran_env, only : r64=>real64
      use mod_engine_rigidMC
      implicit none
      private
      public :: &
         load_rigidMC,&
         globalize_rigidMC,&
         customize_rigidMC,&
         cost_rigidMC,&
         task_rigidMC
      character(:), public, allocatable :: &
         custom_drive,&
         cost_drive,&
         task_drive
      contains
#        include "load_rigidMC_inc.f90"
#        include "globalize_rigidMC_inc.f90"
#        include "customize_rigidMC_inc.f90"
#        include "cost_rigidMC_inc.f90"
#        include "task_rigidMC_inc.f90"
      !****************************************************************
      !  subroutine update_wcoors()
      !****************************************************************
#        include "../../../bridge/update_wcoors_inc.f90"
         !end
      !****************************************************************
      !  subroutine update_lcoors()
      !****************************************************************
#        include "../../../bridge/update_lcoors_inc.f90"
         !end
      !****************************************************************
      !  subroutine exchange_wcoors()
      !****************************************************************
#        define GEO_TYPE geo
#        define GEO_TYPE_STR 'geo'
#        define WCOORS rcoors
#        include "../../../bridge/exchange_wcoors_inc.f90"
        !end
      !****************************************************************
      !  subroutine resolve_geo_uptr()
      !****************************************************************
#        define GEO_TYPE geo
#        include "../../../bridge/resolve_geo_uptr_inc.f90"
        !end
      end
