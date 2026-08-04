!**********************************************************************
      subroutine APPEND(arr,elem,last)
!**********************************************************************
      implicit none
#ifdef STR_ARR_TYPE
      character(:), allocatable,target :: arr(:)
      character(*) :: elem
      character(:), allocatable :: tmp(:)
      character(:), pointer, optional :: last
      integer :: mlen
      !***
      if(.not.allocated(arr)) arr=[character(0)::]
      mlen=max(len(arr),len(elem))
      tmp=[character(mlen)::arr,elem]
      call move_alloc(tmp,arr)
      !***
#else
      !***
      ARR_TYPE,allocatable,target :: arr(:)
      ARR_TYPE :: elem
      ARR_TYPE,pointer,optional :: last
      ARR_TYPE,allocatable :: tmp(:)
      integer :: n,none
      allocate(arr(0),stat=none)
      n=size(arr)
      allocate(tmp(n+1))
      tmp(:n)=arr
      tmp(n+1)=elem
      deallocate(arr)
      call move_alloc(tmp,arr)
      !***
#endif
      !***
      if(present(last)) then
         associate(n=>size(arr))
            last=>arr(n)
         end associate
      endif
      !***
      end
#     undef STR_ARR_TYPE
#     undef ARR_TYPE
#     undef APPEND
