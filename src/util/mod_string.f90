!**********************************************************************
      module mod_string
!**********************************************************************
      use iso_fortran_env, only: &
         r32=>real32,&
         r64=>real64,&
         i32=>int32,&
         i64=>int64
      implicit none
      contains
      !****************************************************************
         function expand_tilde(path) result(res)
      !****************************************************************
         implicit none
         character(*) :: &
            path
         character(:), allocatable :: &
            res,&
            home
         integer :: &
            n
         !***
         res=trim(path)
         if(len(res).eq.0) return
         if(res(1:1).ne.'~') return
         call get_environment_variable('HOME',length=n)
         if(n.eq.0) return
         allocate(character(n) :: home)
         call get_environment_variable('HOME',home)
         res=home//res(2:)
         !***
         end
      !****************************************************************
         function to_string(x,nddigs,exp_form) result(str)
      !****************************************************************
         use mod_catch, only : catch_error
         class(*) :: &
            x
         integer, optional :: &
            nddigs
         logical, optional :: &
            exp_form
         type lview_
            integer :: nddigs=0
            logical :: exp_form=.false.
         end type
         type(lview_) :: &
            lview
         character(:), allocatable :: &
            str,&
            frm,&
            spref
         real(r64) :: &
            r
         integer :: &
            wstr,&
            wfrm,&
            ln,&
            k
         !***
         if(present(nddigs)) &
            lview%nddigs=nddigs
         if(present(exp_form)) &
            lview%exp_form=exp_form
         select type(x)
         type is(logical)
            str=transfer(merge('T','F',x),'x')
            return
         type is(character(*))
            str=x
            return
         type is(real(r32))
            r=x
         type is(real(r64))
            r=x
         type is(integer(i32))
            r=x
         type is(integer(i64))
            r=x
         class default
            call catch_error(err=.true.,&
                              msg='Unsupported type.',&
                              proc='tostring')
         end select
         associate(ndd   => lview%nddigs,&
                   eform => lview%exp_form,&
                   rred => merge(1.d0,dabs(r),r.eq.0.d0))
            spref=transfer(merge('es','f ',eform),'xx')
            k=ceiling(dabs(dlog10(rred)))
            wstr=ndd+2*k+6
            wfrm=wstr+ndd+7
            allocate(character(len=wfrm)::frm)
            write(frm,'("(",a,i0,".",i0,")")') &
               trim(spref),wstr,ndd
            allocate(character(len=wstr)::str)
            write(str,frm) r
            str=trim(adjustl(str))
            ln=len(str)
            if(str(ln:ln).eq.'.') str=str(:ln-1)
            if(str.eq.'-0') str='0'
         end associate
         end
      end
