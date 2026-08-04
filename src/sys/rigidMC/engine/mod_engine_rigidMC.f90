!***********************************************************************
      module mod_engine_rigidMC
!***********************************************************************
      use iso_fortran_env, only : r64=>real64
      use mod_point_group, only : pgop_, pgroup
      implicit none

      real(r64), parameter :: &
         deg2rad=dacos(-1.d0)/180.d0
      !***
#     include "control/type_control_rigidMC_inc.f90"
#     include "layout/type_layout_rigidMC_inc.f90"
#     include "refframe/type_refframe_rigidMC_inc.f90"
#     include "sentinel/type_sentinel_rigidMC_inc.f90"
#     define PQ_TYPE geo_prior_queue
#     define PQ_OBJ_TYPE type(geo), pointer
#     define PQ_WEIGHT_TYPE real(r64)
#     include "../../../util/dep/prior_queue_inc.f90"
#     include "landscape/type_landscape_rigidMC_inc.f90"
#     include "geo/type_geo_rigidMC_inc.f90"
#     include "coorbox/type_coorbox_rigidMC_inc.f90"
#     include "monbox/type_monbox_rigidMC_inc.f90"
      !***
      type(control_), target :: &
         control
      type(layout_), target :: &
         layout
      type(refframe_), target :: &
         refframe
      type(landscape_), target :: &
         landscape
      type(sentinel_), target :: &
         sentinel
      !***
      contains
      !*****************************************************************
#        define PQ_TYPE geo_prior_queue
      !*****************************************************************
#        define PQ_OBJ_TYPE type(geo), pointer
#        define PQ_WEIGHT_TYPE real(r64)
#        define PQ_IMPL
#        include "../../../util/dep/prior_queue_inc.f90"
         !end
#        include "control/init_control_rigidMC_inc.f90"
#        include "layout/init_layout_rigidMC_inc.f90"
#        include "refframe/init_refframe_rigidMC_inc.f90"
#        include "sentinel/init_sentinel_rigidMC_inc.f90"
#        include "landscape/init_landscape_rigidMC_inc.f90"
#        include "coorbox/init_coorbox_rigidMC_inc.f90"
#        include "monbox/init_monbox_rigidMC_inc.f90"
#        include "geo/init_geo_rigidMC_inc.f90"
#        include "geo/manip_geo_rigidMC_inc.f90"
#        include "geo/reparam_geo_rigidMC_inc.f90"
#        include "geo/sync_geo_rigidMC_inc.f90"
#        include "geo/calc_geo_rigidMC_inc.f90"
#        include "geo/assess_geo_rigidMC_inc.f90"
#        include "geo/screen_geo_rigidMC_inc.f90"
#        include "geo/stream_geo_rigidMC_inc.f90"
      !***
      end
