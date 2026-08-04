!***********************************************************************
      function get_clarg_filename(iarg) result(filename)
!***********************************************************************
      implicit none
      character(:), allocatable :: &
         filename
      integer :: &
         iarg,&
         m
      !***
      call get_command_argument(iarg,length=m)
      allocate(character(m)::filename)
      call get_command_argument(iarg,filename)
      !***
      end
!***********************************************************************
      subroutine read_input(filename)
!***********************************************************************
      use mod_string, only : expand_tilde
      implicit none
      logical :: &
         eof
      character(*) :: &
         filename
      character(:), allocatable :: &
         row
      type(buffer), pointer :: &
         buf_ptr
      !***
      allocate(buf_ptr)
      call append_buffer(arr=buffers,elem=buf_ptr,last=buf_ptr)
      buf_ptr%filename=filename
      open(1,file=expand_tilde(filename))
      do
         call get_next_row(funit=1,row=row,eof=eof)
         if(eof) exit
         call append_string(arr=buf_ptr%rows,elem=row)
      end do
      close(1)
      !***
      contains
      !*****************************************************************
         subroutine get_next_row(funit,row,eof)
      !*****************************************************************
         use iso_fortran_env, only : iostat_eor, iostat_end
         implicit none
         integer :: &
            funit
         character(:), allocatable :: &
            row
         logical, optional :: &
            eof
         integer :: &
            ios
         character :: &
            ch
         !***
         row=''
         if(present(eof)) eof=.false.
         do
            read(funit,'(a)',advance='no',iostat=ios) ch
            select case(ios)
               case(0)
                  row=row//ch
               case(iostat_eor)
                  exit
               case(iostat_end)
                  if(present(eof)) eof=.true.
                  exit
               case default
                  call catch_error(&
                     err=.true.,&
                     msg='read failed.',&
                     proc='get_next_row')
            end select
         end do
         !***
         end
      end
!***********************************************************************
      function get_buf_ptr(filename) result(buf_ptr)
!***********************************************************************
      implicit none
      character(*) :: &
         filename
      type(buffer), pointer :: &
         buf_ptr
      integer :: &
         i,&
         n
      !***
      buf_ptr=>null()
      n=0
      do i=1,size(buffers)
         associate(act_buf=>buffers(i))
            if(act_buf%filename.eq.filename) then
               buf_ptr=>act_buf
               n=n+1
            endif
         end associate
      end do
      call catch_error(&
              err=.not.associated(buf_ptr),&
              msg='unknown buffer filename.',&
              proc='get_buf_ptr')
      call catch_error(&
              err=n.ne.1,&
              msg='duplicated buffer filename.',&
              proc='get_buf_ptr')
      !***
      end
!***********************************************************************
      function count_key(filename,section,key) result(cnt)
!***********************************************************************
      implicit none
      character(*) :: &
         filename,&
         section,&
         key
      character(:), allocatable :: &
         row
      type(buffer), pointer :: &
         buf_ptr
      integer &
         i,&
         j,&
         ind,&
         cnt
      !***
      buf_ptr=>get_buf_ptr(filename)
      associate(rows=>buf_ptr%rows)
         do i=1,size(rows)
            row=rows(i)
            ind=index(row,section)
            if(ind.ne.0) exit
         end do
         cnt=0
         do j=i+1,size(rows)
            row=rows(j)
            if(len_trim(row).eq.0) cycle
            if(index(adjustl(row),'#').eq.1) exit
            ind=index(row,key)
            if(ind.ne.0) cnt=cnt+1
         end do
      end associate
      !***
      end
!***********************************************************************
      subroutine &
              parse_keyval(filename,section,key,val,ival,rval,unit)
!***********************************************************************
      implicit none
      character(*) :: &
         filename,&
         section,&
         key
      character(:), allocatable, optional :: &
         val,&
         unit
      integer, optional :: &
         ival
      real(r64), optional :: &
         rval
      type(buffer), pointer :: &
         buf_ptr
      character(:), allocatable :: &
         row,&
         tval,&
         cand_unit,&
         fpart
      logical :: &
         found
      integer :: &
         i,j,&
         ind,&
         istart,&
         iend,&
         ios,&
         vpres,&
         ipres,&
         rpres
      !***
      buf_ptr=>get_buf_ptr(filename)
      associate(rows=>buf_ptr%rows)
         do i=1,size(rows)
            row=trim(rows(i))
            ind=index(row,section)
            if(ind.ne.0) exit
         end do
         found=.false.
         do j=i+1,size(rows)
            row=trim(rows(j))
            if(len_trim(row).eq.0) cycle
            if(index(adjustl(row),'#').eq.1) exit
            ind=index(row,'=')
            if(ind.eq.0) cycle
            fpart=row(:ind-1)
            if(present(unit)) then
               istart=index(fpart,'[')+1
               iend=index(fpart,']')-1
               cand_unit=fpart(istart:iend)
               fpart=fpart(:istart-2)
            endif
            fpart=adjustl(fpart)
            if(fpart.ne.key) cycle
            if(present(unit)) unit=cand_unit
            tval=row(ind+1:)
            vpres=merge(1,0,present(val))
            ipres=merge(1,0,present(ival))
            rpres=merge(1,0,present(rval))
            call catch_error(&
                     err=vpres+ipres+rpres.ne.1,&
                     msg='exactly one optional argument '//&
                         'is required.',&
                     proc='parse_keyval')
            if(vpres.eq.1) then
               val=trim(adjustl(tval))
               ios=0
            elseif(ipres.eq.1) then
               read(tval,*,iostat=ios) ival
            elseif(rpres.eq.1) then
               read(tval,*,iostat=ios) rval
            endif
            call catch_error(&
                    err=ios.ne.0,&
                    msg='value conversion failed for key '//&
                        '"'//trim(key)//'".',&
                    proc='parse_keyval')
            found=.true.
         end do
         call catch_error(&
                 err=.not.found,&
                 msg='key "'//key//'" not found in section "'//&
                     section//'".',&
                 proc='parse_keyval')
      end associate
      !***
      end
!***********************************************************************
      subroutine get_key_list(filename,section,base_key,keys)
!***********************************************************************
      use mod_string, only : to_string
      implicit none
      character(*) :: &
         filename,&
         section,&
         base_key
      character(:), allocatable :: &
         keys(:)
      integer :: &
         i
      character(:), allocatable :: &
         key,&
         sect
      !***
      allocate(sect,source=section) !NEEDED to circumvent a GNU bug
      associate(n=>count_key(filename,sect,base_key))
         do i=1,n
            key=base_key//to_string(i)
            call append_string(arr=keys,elem=key)
         end do
      end associate
      !***
      end
