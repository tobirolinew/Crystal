#ifdef DYNAMIC
#  define READ_XYZ_NAME read_xyz_dynamic
#  define EXTRA allocatable
#  define ROW :
#  define COL :
#else
#  define READ_XYZ_NAME read_xyz
#  define EXTRA target
#  define ROW *
#  define COL 3
#endif
!***********************************************************************
      subroutine &
           READ_XYZ_NAME(xyz_file,anums,symbs,carts,title,natoms,matoms)
!***********************************************************************
      use mod_catch, only : catch_error
      use mod_pertab
      use mod_string, only : expand_tilde
      implicit none
      integer, optional :: &
         natoms,&
         matoms
      character(:), allocatable, optional :: &
         title
      integer :: &
         i,j,&
         n,&
         ios
      character(*) :: &
         xyz_file
      character(2), EXTRA :: &
         symbs(ROW)
      integer, EXTRA :: &
         anums(ROW)
      real(r64), EXTRA :: &
         carts(COL,ROW)
      character :: &
         ch
      !***
      open(1,file=expand_tilde(xyz_file),status='old',iostat=ios)
      call catch_error(err=ios.ne.0,&
               msg='cannot open "xyz_file".',&
               proc='read_xyz')
      read(1,*,iostat=ios) n
      if(present(natoms)) natoms=n
      if(present(matoms)) &
         call catch_error(err=n.gt.matoms,&
                  msg='number of atoms is larger than the '//&
                      'allocated array size.',&
                  proc='read_xyz')
      call catch_error(err=ios.ne.0,&
               msg='failed to read number of atoms.',&
               proc='read_xyz')
      call catch_error(err=n.le.0,&
               msg='number of atoms must be positive.',&
               proc='read_xyz')
      call catch_error(err=n.gt.huge(0),&
               msg='unreasonably large number of atoms.',&
               proc='read_xyz')
      if(present(title)) then
         title= ''
         do
            read(1,'(a)',advance='no',iostat=ios) ch
            if(ios.lt.0.or.ch.eq.new_line('a')) exit
            title=title//ch
         end do
      else
         read(1,*,iostat=ios)
      endif
# ifdef DYNAMIC 
      symbs=spread('  ',1,n)
      anums=spread(0,1,n)
      carts=spread([0.d0,0.d0,0.d0],2,n)
# endif
      do i=1,n
         associate(symb => symbs(i),&
                   anum => anums(i))
            read(1,*,iostat=ios) symb,carts(:,i)
            call catch_error(err=ios.ne.0,&
                     msg='failed to parse atomic symbol '//&
                         'and coordinates.',&
                     proc='read_xyz')
            do j=0,nelems
               if(symb.eq.elem_symbs(j)) exit
            end do
            call catch_error(err=j.eq.nelems+1,&
                     msg='unknown atomic symbol "'//symb//'"',&
                     proc='read_xyz')
            anum=j
         end associate
      end do
      close(1)
      do i=1,n
         do j=1,i-1
            associate(d=>norm2(carts(:,i)-carts(:,j)))
               call catch_error(err=d.le.epsilon(0.d0),&
                        msg='overlaping atoms detected.',&
                        proc='read_xyz')
            end associate
         end do
      end do
      !***
      end
#undef DYNAMIC
#undef READ_XYZ_NAME
#undef EXTRA
#undef ROW
#undef COL
