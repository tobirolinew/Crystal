!***********************************************************************
      module mod_xinfo
!***********************************************************************
      use iso_c_binding
      implicit none
      interface
      !*****************************************************************
         function readlink(path,buf,bufsiz) bind(C,name="readlink")
      !*****************************************************************
         import :: &
            c_char,&
            c_size_t,&
            c_long
         character(c_char), intent(in) :: &
            path(*)
         character(c_char), intent(out) :: &
            buf(*)
         integer(c_size_t), value :: &
            bufsiz
         integer(c_long) :: &
            readlink
         end
      !***
      end interface
      contains
      !*****************************************************************
         function get_xpath() result(path)
      !*****************************************************************
         implicit none
         character(len=:), allocatable :: &
            path
         character(c_char) :: &
            cbuf(4096)
         integer(c_long) :: &
            n
         !***
         n=readlink("/proc/self/exe"//c_null_char,&
                    cbuf,&
                    size(cbuf,kind=c_size_t))
         if(n.le.0) then
            path=""
            return
         endif
         path=transfer(cbuf(1:n),repeat(" ",int(n)))
         !***
         end
      !*****************************************************************
         function get_xdir() result(dir)
      !*****************************************************************
         implicit none
         character(len=:), allocatable :: &
            dir
         integer :: &
            pos
         !***
         dir=get_xpath()
         pos=scan(dir,"/",back=.true.)
         if(pos.gt.1) then
            dir=dir(:pos-1)
         else
            dir="/"
         endif
         !***
         end
      !***
      end
