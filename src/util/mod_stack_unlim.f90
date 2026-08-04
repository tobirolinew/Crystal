!***********************************************************************
      module mod_stack_unlim
!***********************************************************************
      implicit none
      private
      public :: &
         stack_unlim
      contains
      !*****************************************************************
         subroutine stack_unlim()
      !*****************************************************************
         use iso_c_binding
         implicit none
         interface
            function c_getrlimit(res,rlim) bind(c,name='getrlimit')
            import :: c_int,c_ptr
            integer(c_int), value :: res
            type(c_ptr), value :: rlim
            integer(c_int) :: c_getrlimit
            end
            function c_setrlimit(res,rlim) bind(c,name='setrlimit')
            import :: c_int,c_ptr
            integer(c_int), value :: res
            type(c_ptr), value :: rlim
            integer(c_int) :: c_setrlimit
            end
            function c_getenv(name) bind(c,name='getenv')
            import :: c_ptr
            type(c_ptr), value :: name
            type(c_ptr) :: c_getenv
            end
            function c_setenv(name,val,over) bind(c,name='setenv')
            import :: c_ptr,c_int
            type(c_ptr), value :: name,val
            integer(c_int), value :: over
            integer(c_int) :: c_setenv
            end
            function c_execv(path,argv) bind(c,name='execv')
            import :: c_ptr,c_int
            type(c_ptr), value :: path
            type(c_ptr) :: argv(*)
            integer(c_int) :: c_execv
            end
            function c_readlink(path,buf,sz) bind(c,name='readlink')
            import :: c_ptr,c_long
            type(c_ptr), value :: path,buf
            integer(c_long), value :: sz
            integer(c_long) :: c_readlink
            end
         end interface
         integer(c_int), parameter :: &
            RLIMIT_STACK=3
         integer(c_long), parameter :: &
            RLIM_INFINITY=-1_c_long,&
            MAINSTK=256_c_long*1024_c_long**3
         integer(c_long), target :: &
            rlim(2)
         integer(c_long) :: &
            nlink
         integer(c_int) :: &
            ierr,i,l,nargs
         character(kind=c_char,len=8192), &
            allocatable, target :: &
            av(:)
         type(c_ptr), allocatable :: &
            argv(:)
         character(kind=c_char,len=16), target :: &
            flag,one,ompsz,ompval,proc
         character(kind=c_char,len=4096), target :: &
            exe
         !***
         flag='CRYSTAL_STACKED'//c_null_char
         if(.not.c_associated(c_getenv(c_loc(flag)))) then
            one='1'//c_null_char
            ierr=c_setenv(c_loc(flag),c_loc(one),1_c_int)
            ompsz='OMP_STACKSIZE'//c_null_char
            ompval='128M'//c_null_char
            ierr=c_setenv(c_loc(ompsz),c_loc(ompval),0_c_int)
            ierr=c_getrlimit(RLIMIT_STACK,c_loc(rlim))
            rlim(1)=MAINSTK
            if(rlim(2).ne.RLIM_INFINITY.and.rlim(2).lt.MAINSTK) &
               rlim(1)=rlim(2)
            ierr=c_setrlimit(RLIMIT_STACK,c_loc(rlim))
            nargs=command_argument_count()
            allocate(av(0:nargs),argv(0:nargs+1))
            do i=0,nargs
               call get_command_argument(i,av(i),l)
               av(i)(l+1:l+1)=c_null_char
               argv(i)=c_loc(av(i))
            end do
            argv(nargs+1)=c_null_ptr
            proc='/proc/self/exe'//c_null_char
            nlink=c_readlink(c_loc(proc),c_loc(exe),4095_c_long)
            if(nlink.gt.0_c_long.and.nlink.lt.4096_c_long) then
               exe(nlink+1:nlink+1)=c_null_char
            else
               exe='/proc/self/exe'//c_null_char
            endif
            ierr=c_execv(c_loc(exe),argv)
         endif
         !***
         end
      !***
      end
