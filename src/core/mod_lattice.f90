#     include "dep/mod_lpo_inc.f90"
#     include "dep/mod_lpo_lock_set_inc.f90"
#     include "dep/mod_lpo_cost_set_inc.f90"
!***********************************************************************
      module mod_lattice
!***********************************************************************
      use iso_fortran_env, only : r64=>real64
      use mod_lpo, only : lpo
      use mod_lpo_lock_set, only : lock_set_=>ffh_t
      use mod_lpo_cost_set, only : contrac, cost_set_=>ffh_t
      !***
#     define PQ_TYPE prior_queue
#     define PQ_OBJ_TYPE type(lpo), pointer
#     define PQ_WEIGHT_TYPE real(r64)
#     include "../util/dep/prior_queue_inc.f90"
      !***
      contains
      !*****************************************************************
#     define PQ_TYPE prior_queue
      !*****************************************************************
#     define PQ_OBJ_TYPE type(lpo), pointer
#     define PQ_WEIGHT_TYPE real(r64)
#     define PQ_IMPL
#     include "../util/dep/prior_queue_inc.f90"
         !end
      !*****************************************************************
         function instantiate(lcoors,sidx,parent) result(lpo_ptr)
      !*****************************************************************
         use mod_global, only : fix_wcoors
         implicit none
         integer :: &
            lcoors(:)
         integer, optional :: &
            sidx
         type(lpo), optional, pointer :: &
            parent
         type(lpo), pointer :: &
            lpo_ptr
         integer :: &
            i
         !***
         lpo_ptr=>null()
         i=0
         if(present(sidx)) then
            i=iabs(sidx)
            if(fix_wcoors(i)) return
         endif
         allocate(lpo_ptr)
         lpo_ptr%lcoors=lcoors
         if(i.ne.0) &
            lpo_ptr%lcoors(i)=lcoors(i)+sign(1,sidx)
         lpo_ptr%geo_uptr=>null()
         if(present(parent)) then
            lpo_ptr%parent=>parent
         else
            lpo_ptr%parent=>null()
         endif
         !***
         end
      end
