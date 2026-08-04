#ifdef __GFORTRAN__
#  define PQ_EXPAND(x) x
#  define PQ_PASTE(a,b) a/**/b
#  define PQ_CONCAT(a,b) PQ_PASTE(PQ_EXPAND(a),b)
#else
#  define PQ_PASTE(a,b) a##b
#  define PQ_CONCAT(a,b) PQ_PASTE(a,b)
#endif
#ifndef PQ_IMPL
!***********************************************************************
!     prior_queue decelarations
!***********************************************************************
      ! weighted object slot
      type :: PQ_CONCAT(PQ_TYPE,_wobj_)
         PQ_OBJ_TYPE :: &
            obj
         PQ_WEIGHT_TYPE :: &
            weight
      end type
      ! binary max-heap of weighted object pointers
      type, public :: PQ_TYPE
         type(PQ_CONCAT(PQ_TYPE,_wobj_)), allocatable :: &
            heap(:)
         integer :: &
            length=0
         contains
            procedure, public :: &
               is_empty       => PQ_CONCAT(PQ_TYPE,_is_empty),&
               at             => PQ_CONCAT(PQ_TYPE,_at),&
               weight_at      => PQ_CONCAT(PQ_TYPE,_weight_at),&
               replace_obj_at => PQ_CONCAT(PQ_TYPE,_replace_obj_at),&
               top            => PQ_CONCAT(PQ_TYPE,_top),&
               weight_top     => PQ_CONCAT(PQ_TYPE,_weight_top),&
               push           => PQ_CONCAT(PQ_TYPE,_push),&
               pop            => PQ_CONCAT(PQ_TYPE,_pop),&
               erase          => PQ_CONCAT(PQ_TYPE,_erase),&
               delete         => PQ_CONCAT(PQ_TYPE,_delete)
            procedure, private :: &
               sift           => PQ_CONCAT(PQ_TYPE,_sift),&
               is_in_range    => PQ_CONCAT(PQ_TYPE,_is_in_range)
      end type
      !end
#endif
#ifdef PQ_IMPL
!***********************************************************************
      logical function PQ_CONCAT(PQ_TYPE,_is_empty)(this) result(res)
!***********************************************************************
      ! True when the queue holds no object.
      implicit none
      class(PQ_TYPE) :: &
            this
      !***
      res=this%length.eq.0
      !***
      end
!***********************************************************************
      function PQ_CONCAT(PQ_TYPE,_at)(this,i) result(obj)
!***********************************************************************
      ! Returns a pointer to the i-th object of the heap (1-based).
      implicit none
      class(PQ_TYPE) :: &
         this
      integer :: &
         i
      PQ_OBJ_TYPE :: &
         obj
      !***
      obj=>this%heap(i)%obj
      !***
      end
!***********************************************************************
      function PQ_CONCAT(PQ_TYPE,_weight_at)(this,i) result(weight)
!***********************************************************************
      ! Returns the weight of the i-th slot, -huge out of range.
      implicit none
      class(PQ_TYPE) :: &
         this
      integer :: &
         i
      PQ_WEIGHT_TYPE :: &
         weight
      !***
      if(this%is_in_range(i)) then
         weight=this%heap(i)%weight
      else
         weight=-huge(0.d0)
      endif
      !***
      end
!***********************************************************************
      subroutine PQ_CONCAT(PQ_TYPE,_replace_obj_at)(this,i,obj)
!***********************************************************************
      ! Retargets the i-th slot without touching its weight.
      use mod_catch, only : catch_error
      implicit none
      class(PQ_TYPE) :: &
         this
      integer :: &
         i
      PQ_OBJ_TYPE :: &
         obj
      !***
      call catch_error(err=.not.this%is_in_range(i),&
                           msg='Invalid index.',&
                           proc='prior_queue%replace_obj_at')
      this%heap(i)%obj=>obj
      !***
      end
!***********************************************************************
      function PQ_CONCAT(PQ_TYPE,_top)(this) result(obj)
!***********************************************************************
      ! Returns a pointer to the heaviest object.
      implicit none
      class(PQ_TYPE) :: &
         this
      PQ_OBJ_TYPE :: &
         obj
      !***
      obj=>this%at(1)
      !***
      end
!***********************************************************************
      function PQ_CONCAT(PQ_TYPE,_weight_top)(this) result(weight)
!***********************************************************************
      ! Returns the largest weight in the queue.
      implicit none
      class(PQ_TYPE) :: &
         this
      PQ_WEIGHT_TYPE :: &
         weight
      !***
      weight=this%weight_at(1)
      !***
      end
!***********************************************************************
      subroutine PQ_CONCAT(PQ_TYPE,_push)(this,obj,weight)
!***********************************************************************
      ! Inserts obj with the given weight. Doubles capacity when full.
      implicit none
      integer, parameter :: &
         min_capac=1000
      class(PQ_TYPE) :: &
         this
      PQ_OBJ_TYPE :: &
         obj
      PQ_WEIGHT_TYPE :: &
         weight
      integer :: &
         i
      type(PQ_CONCAT(PQ_TYPE,_wobj_)), allocatable :: &
         new_heap(:)
      !***
      associate(length=>this%length)
         if(.not.allocated(this%heap)) then
            allocate(this%heap(min_capac))
         else
            associate(heap_size=>size(this%heap))
               if(length.eq.heap_size) then
                  associate(new_capac=>max(min_capac,2*heap_size))
                     allocate(new_heap(new_capac))
                     new_heap(:length)=this%heap
                     call move_alloc(new_heap,this%heap)
                  end associate
               endif
            end associate
         endif
         i=length+1
         length=i
         associate(heap_entry=>this%heap(i))
            heap_entry%obj=>obj
            heap_entry%weight=weight
         end associate
         call this%sift(i,'up')
      end associate
      !***
      end
!***********************************************************************
      function PQ_CONCAT(PQ_TYPE,_pop)(this) result(obj)
!***********************************************************************
      ! Removes and returns the heaviest object.
      use mod_catch, only : catch_error
      implicit none
      class(PQ_TYPE) :: &
         this
      PQ_OBJ_TYPE :: &
         obj
      !***
      call catch_error(&
               err=this%is_empty(),&
               msg='empty queue.',&
               proc='prior_queue%pop')
      obj=>this%top()
      call this%erase()
      !***
      end
!***********************************************************************
      subroutine PQ_CONCAT(PQ_TYPE,_erase)(this,i)
!***********************************************************************
      ! Removes the i-th object (the top by default). Does not
      ! deallocate the object itself.
      use mod_catch, only : catch_error
      implicit none
      class(PQ_TYPE) :: &
         this
      integer, optional :: &
         i
      integer :: &
         iloc
      !***
      associate(length => this%length,&
                heap   => this%heap)
         if(present(i)) then
            iloc=i
         else
            iloc=1
         endif
         call catch_error(&
                  err=.not.this%is_in_range(iloc),&
                  msg='invalid index.',&
                  proc='prior_queue%erase')
         heap(iloc)=heap(length)
         length=length-1
         call this%sift(iloc,'up')
         call this%sift(iloc,'down')
         if(present(i)) i=iloc
      end associate
      !***
      end
!***********************************************************************
      subroutine PQ_CONCAT(PQ_TYPE,_delete)(this)
!***********************************************************************
      ! Deallocates the heap array. Does not deallocate the objects.
      implicit none
      class(PQ_TYPE) :: &
            this
      !***
      deallocate(this%heap)
      this%length=0
      !***
      end
!***********************************************************************
      subroutine PQ_CONCAT(PQ_TYPE,_sift)(this,i,dir)
!***********************************************************************
      ! Restores the heap property from slot i, upwards or downwards.
      implicit none
      class(PQ_TYPE) :: &
         this
      character(*) :: &
         dir
      type(PQ_CONCAT(PQ_TYPE,_wobj_)) :: &
         last
      integer :: &
         inext,ileft,iright,&
         i
      !***
      if(.not.this%is_in_range(i)) return
      associate(heap=>this%heap)
         last=heap(i)
         select case(dir)
         case('up')
            do while(i.ne.1)
               inext=i/2
               associate(next=>heap(inext))
                  if(next%weight.ge.last%weight) exit
                  heap(i)=next
                  i=inext
               end associate
            end do
         case('down')
            associate(half_length=>this%length/2)
               do while(i.le.half_length)
                  ileft=2*i
                  iright=ileft+1
                  associate(wleft  => this%weight_at(ileft),&
                            wright => this%weight_at(iright))
                     if(wright.gt.wleft) then
                        inext=iright
                     else
                        inext=ileft
                     endif
                     associate(next=>heap(inext))
                        if(last%weight.ge.next%weight) exit
                        heap(i)=next
                        i=inext
                     end associate
                  end associate
               end do
            end associate
         end select
         heap(i)=last
      end associate
      !***
      end
!***********************************************************************
      logical function PQ_CONCAT(PQ_TYPE,_is_in_range)(this,i) result(res)
!***********************************************************************
      ! True when i addresses an occupied slot.
      implicit none
      class(PQ_TYPE) :: &
         this
      integer :: &
         i
      !***
      res=i.ge.1.and.i.le.this%length
      !***
      end
#endif

#undef PQ_TYPE
#undef PQ_OBJ_TYPE
#undef PQ_WEIGHT_TYPE
#undef PQ_IMPL
#undef PQ_EXPAND
#undef PQ_PASTE
#undef PQ_CONCAT
