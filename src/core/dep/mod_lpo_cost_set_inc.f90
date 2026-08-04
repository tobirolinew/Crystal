!***********************************************************************
      module mod_lpo_cost_set
!***********************************************************************
      use iso_fortran_env, only : int32, int64, r64=>real64
      implicit none
      real(r64), public :: &
         contrac
#     define FFH_KEY_TYPE real(r64)
#     define FFH_CUSTOM_KEYS_EQUAL
#     define FFH_CUSTOM_CONVERT_KEY
#     define FFH_ENABLE_INT64
#     include "../../vendor/ffhash/ffhash_inc.f90"
      !*****************************************************************
         pure function convert_key(key) result(buf)
      !*****************************************************************
         integer, parameter :: &
            nbytes=storage_size(1_int64)/8
         real(r64), intent(in) :: &
            key
         integer(int64) :: &
            ikey
         character(len=nbytes) :: &
            buf
         !***
         ikey=nint(key/contrac,int64)
         buf=transfer(ikey,buf)
         !***
         end
      !*****************************************************************
         pure logical function keys_equal(a,b)
      !*****************************************************************
         real(r64), intent(in) :: &
            a,b
         integer(int64) :: &
            ia,ib
         !***
         ia=nint(a/contrac,int64)
         ib=nint(b/contrac,int64)
         keys_equal=ia.eq.ib
         !***
         end
      end
