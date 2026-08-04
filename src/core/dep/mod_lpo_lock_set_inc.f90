!***********************************************************************
      module mod_lpo_lock_set
!***********************************************************************
      use iso_fortran_env, only : int32, int64
      use mod_lpo, only : lpo
      implicit none
#     define FFH_KEY_TYPE type(lpo)
#     define FFH_CUSTOM_KEYS_EQUAL
#     define FFH_CUSTOM_CONVERT_KEY
#     define FFH_ENABLE_INT64
#     include "../../vendor/ffhash/ffhash_inc.f90"
      !*****************************************************************
         pure logical function keys_equal(a,b)
      !*****************************************************************
         type(lpo), intent(in) :: &
            a,b
         !***
         keys_equal=size(a%lcoors).eq.size(b%lcoors)
         if(keys_equal) keys_equal=all(a%lcoors.eq.b%lcoors)
         !***
         end
      !*****************************************************************
         pure function convert_key(key) result(buf)
      !*****************************************************************
         type(lpo), intent(in) :: &
            key
         character(len=4*size(key%lcoors)) :: &
            buf
         !***
         buf=transfer(key%lcoors,buf)
         !***
         end
      end
