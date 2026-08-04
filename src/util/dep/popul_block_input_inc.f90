!***********************************************************************
      subroutine extend_schema(sc,section,key,valtype,unit_flag)
!***********************************************************************
      implicit none
      type(schema) :: &
         sc
      character(*) :: &
         section,&
         key,&
         valtype
      logical :: &
         unit_flag
      !***
      call append_string(arr=sc%sections,elem=section)
      if(len(key).gt.0) call append_string(arr=sc%keys,elem=key)
      call append_string(arr=sc%valtypes,elem=valtype)
      call append_logical(arr=sc%unit_flags,elem=unit_flag)
      !***
      end
!***********************************************************************
      subroutine schema2store(sc,st)
!***********************************************************************
      use mod_catch, only : catch_error
      implicit none
      type(schema) :: &
         sc
      type(store) :: &
         st
      character(:), allocatable :: &
         section,&
         key,&
         valtype,&
         sval,&
         unit
      real(r64) :: &
         rval
      logical :: &
         unit_flag
      integer :: &
         i,&
         ival
      !***
      call catch_local_errors('precond')
      do i=1,size(sc%keys)
         section=trim(sc%sections(i))
         key=trim(sc%keys(i))
         valtype=trim(sc%valtypes(i))
         unit_flag=sc%unit_flags(i)
         call append_string(arr=st%keys,elem=key)
         select case(valtype)
            case('string')
               if(unit_flag) then
                  call parse_keyval(sc%filename,section,key,&
                                    val=sval,unit=unit)
                  call append_string(arr=st%units,elem=unit)
               else
                  call parse_keyval(sc%filename,section,key,val=sval)
                  call append_string(arr=st%units,elem='')
               endif
               call append_string(arr=st%vals,elem=sval)
               call append_i32(arr=st%ivals,elem=huge(0))
               call append_r64(arr=st%rvals,elem=huge(0.d0))
            case('integer')
               if(unit_flag) then
                  call parse_keyval(sc%filename,section,key,&
                                    ival=ival,unit=unit)
                  call append_string(arr=st%units,elem=unit)
               else
                  call parse_keyval(sc%filename,section,key,ival=ival)
                  call append_string(arr=st%units,elem='')
               endif
               call append_i32(arr=st%ivals,elem=ival)
               call append_r64(arr=st%rvals,elem=huge(0.d0))
               call append_string(arr=st%vals,elem='')
            case('real')
               if(unit_flag) then
                  call parse_keyval(sc%filename,section,key,&
                                    rval=rval,unit=unit)
                  call append_string(arr=st%units,elem=unit)
               else
                  call parse_keyval(sc%filename,section,key,&
                                    rval=rval)
                  call append_string(arr=st%units,elem='')
               endif
               call append_r64(arr=st%rvals,elem=rval)
               call append_string(arr=st%vals,elem='')
               call append_i32(arr=st%ivals,elem=huge(0))
            case default
               call catch_local_errors('valtype')
         end select
      end do
      !***
      contains
      !*****************************************************************
         subroutine catch_local_errors(opt)
      !*****************************************************************
         use mod_catch, only : catch_error
         implicit none
         character(*) :: &
            opt
         !***
         select case(opt)
            case('precond')
               call catch_error(&
                  err=.not.allocated(sc%filename),&
                  msg='schema%filename is not allocated.',&
                  proc='schema2store')
               call catch_error(&
                  err=.not.allocated(sc%sections),&
                  msg='schema%sections is not allocated.',&
                  proc='schema2store')
               call catch_error(&
                  err=.not.allocated(sc%keys),&
                  msg='schema%keys is not allocated.',&
                  proc='schema2store')
               call catch_error(&
                  err=.not.allocated(sc%valtypes),&
                  msg='schema%valtypes is not allocated.',&
                  proc='schema2store')
               call catch_error(&
                  err=.not.allocated(sc%unit_flags),&
                  msg='schema%unit_flags is not allocated.',&
                  proc='schema2store')
               call catch_error(&
                  err=size(sc%sections).ne.size(sc%keys).or.&
                      size(sc%valtypes).ne.size(sc%keys).or.&
                      size(sc%unit_flags).ne.size(sc%keys),&
                  msg='inconsistent schema array sizes.',&
                  proc='schema2store')
               call catch_error(&
                  err=allocated(st%keys),&
                  msg='store%keys must be unallocated on entry.',&
                  proc='schema2store')
               call catch_error(&
                  err=allocated(st%vals),&
                  msg='store%vals must be unallocated on entry.',&
                  proc='schema2store')
               call catch_error(&
                  err=allocated(st%ivals),&
                  msg='store%ivals must be unallocated on entry.',&
                  proc='schema2store')
               call catch_error(&
                  err=allocated(st%rvals),&
                  msg='store%rvals must be unallocated on entry.',&
                  proc='schema2store')
               call catch_error(&
                  err=allocated(st%units),&
                  msg='store%unit must be unallocated on entry.',&
                  proc='schema2store')
            case('valtype')
               call catch_error(&
                  err=.true.,&
                  msg='unknown value type "'//trim(valtype)//'".',&
                  proc='schema2store')
            case default
               call catch_error(&
                  err=.true.,&
                  msg='unknown catch option"'//trim(opt)//'".',&
                  proc='schema2store')
         end select
         !***
         end
      end
!***********************************************************************
      subroutine print_schema(sc,nspaces)
!***********************************************************************
      implicit none
      type(schema) :: &
         sc
      character(:), allocatable :: &
         p
      integer :: &
         i,&
         nspaces
      !***
      p=repeat(' ',nspaces)
      write(*,'(/a)') p//'*** SCHEMA:'
      write(*,'(a,4x,a)') p,'---'
      write(*,'(a,4x,a)') p,'filename : '//trim(sc%filename)
      write(*,'(a,4x,a)') p,'---'
      do i=1,size(sc%keys)
         associate(sect      => sc%sections(i),&
                   key       => sc%keys(i),&
                   valtype   => sc%valtypes(i),&
                   unit_flag => sc%unit_flags(i))
            write(*,'(a,4x,a)') p,'section  : '//trim(sect)
            write(*,'(a,4x,a)') p,'key      : '//trim(key)
            write(*,'(a,4x,a)') p,'type     : '//trim(valtype)
            write(*,'(a,4x,a)') p,'unit     : '//&
                                   merge('yes','no ',unit_flag)
            write(*,'(a,4x,a)') p,'---'
         end associate
      end do
      !***
      end
!***********************************************************************
      subroutine print_store(st,nspaces)
!***********************************************************************
      implicit none
      type(store) :: &
         st
      character(:), allocatable :: &
         p
      integer :: &
         i,&
         nspaces
      !***
      p=repeat(' ',nspaces)
      write(*,'(/a)') p//'*** STORE:'
      write(*,'(a,4x,a)') p,'---'
      do i=1,size(st%keys)
         associate(key  => st%keys(i),&
                   sval => st%vals(i),&
                   ival => st%ivals(i),&
                   rval => st%rvals(i),&
                   unit => st%units(i))
            write(*,'(a,4x,a)') p,'key  : '//trim(key)
            if(len_trim(sval).ne.0) then
               write(*,'(a,4x,a)') p,'val  : '//trim(sval)
            elseif(ival.ne.huge(0)) then
               write(*,'(a,4x,a,i0)') p,'val  : ',ival
            else
               write(*,'(a,4x,a,es0.16)') p,'val  : ',rval
               if(len_trim(unit).ne.0) &
                  write(*,'(a,4x,a)') p,'unit : '//trim(unit)
            endif
            write(*,'(a,4x,a)') p,'---'
         end associate
      end do
      !***
      end
