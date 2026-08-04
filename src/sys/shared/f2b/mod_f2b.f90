!***********************************************************************
      module mod_f2b
!***********************************************************************
      use mod_euclid, only : cross_product
      use iso_fortran_env, only : r64=>real64, int64
      implicit none
      private
      public :: read_f2b_pars, pot_f2b, pot_f2b_homo, grad_f2b_homo
      logical flex_asym
      integer &
         nbof_mons,&
         nbof_atoms,&
         nbof_sites,&
         nbof_OA_sites
      real(r64),parameter :: &
         ehkcalmol=627.509578263d0,&
         kcalmolcm=349.755029988d0,&
         ehcm=219474.631363d0,&
         bohrA=0.529177210903d0
      character, allocatable :: &
         symbs(:)*2
      real(r64), allocatable :: &
         iso_carts(:,:),&
         w_OA(:,:),&
         B_intra(:),&
         S_intra(:,:)
      integer, allocatable :: &
         lst_OA(:,:),&
         Lambda_intra(:,:),&
         nbof_mon_atoms(:),&
         nbof_mon_OA_sites(:)
#        include "type_f2b_inc.f90"
      type(site), allocatable :: &
         sites(:)
      type(site_pair), allocatable :: &
         site_pairs(:,:)
      real(r64), parameter :: &
         cfdstep=1.d-30
      contains
#        include "pot_f2b_homo_inc.f90"
#        include "pot_f2b_inc.f90"
#        include "cpot_f2b_homo_inc.f90"
#        include "cpot_f2b_inc.f90"
#        include "f2b_func_inc.f90"
#        include "scsqrt_inc.f90"
#        include "add_OA_sites_inc.f90"
#        include "read_f2b_pars_inc.f90"
#        include "grad_f2b_homo_inc.f90"
      !***
      end
