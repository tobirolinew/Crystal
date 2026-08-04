!***********************************************************************
      module mod_block_input
!***********************************************************************
      use iso_fortran_env, only : r64=>real64
      use mod_catch, only : catch_error
      implicit none
      private
      public :: &
         buffer,&
         buffers,&
         get_buf_ptr,&
         get_clarg_filename,&
         read_input,&
         parse_keyval,&
         count_key,&
         get_key_list,&
         schema,&
         store,&
         extend_schema,&
         print_schema,&
         print_store,&
         schema2store
      type buffer
         character(:), allocatable :: &
            filename,&
            rows(:)
      end type
      type store
         character(:), allocatable :: &
            keys(:),&
            vals(:),&
            units(:)
         integer, allocatable :: &
            ivals(:)
         real(r64), allocatable :: &
            rvals(:)
      end type
      type schema
         character(:), allocatable :: &
            filename,&
            sections(:),&
            keys(:),&
            valtypes(:)
         logical, allocatable :: &
            unit_flags(:)
      end type
      type(buffer), target, allocatable :: &
         buffers(:)
      contains
#        include "dep/gather_block_input_inc.f90"
#        include "dep/popul_block_input_inc.f90"
      !*****************************************************************
      !  subroutine append_buffer()
      !*****************************************************************
#        define ARR_TYPE type(buffer)
#        define APPEND append_buffer
#        include "dep/generic_append_inc.f90"
         !end
      !*****************************************************************
      !  subroutine append_string()
      !*****************************************************************
#        define STR_ARR_TYPE
#        define APPEND append_string
#        include "dep/generic_append_inc.f90"
         !end
      !*****************************************************************
      !  subroutine append_r64()
      !*****************************************************************
#        define ARR_TYPE real(r64)
#        define APPEND append_r64
#        include "dep/generic_append_inc.f90"
         !end
      !*****************************************************************
      !  subroutine append_i32()
      !*****************************************************************
#        define ARR_TYPE integer
#        define APPEND append_i32
#        include "dep/generic_append_inc.f90"
         !end
      !*****************************************************************
      !  subroutine append_logical()
      !*****************************************************************
#        define ARR_TYPE logical
#        define APPEND append_logical
#        include "dep/generic_append_inc.f90"
         !end
      end
