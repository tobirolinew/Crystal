!***********************************************************************
      subroutine read_carts(this,xyz_file)
!***********************************************************************
      use mod_euclid, only : read_xyz_dynamic
      implicit none
      character(*) :: &
         xyz_file
      class(geo) :: &
         this
      integer, allocatable :: &
         anums(:)
      character(2), allocatable :: &
         symbs(:)
      !***
      call catch_local_errors('file')
      call read_xyz_dynamic(xyz_file,anums,symbs,this%carts)
      call catch_local_errors('composition')
      !***
      contains
      !*****************************************************************
         subroutine catch_local_errors(opt)
      !*****************************************************************
         use mod_catch, only : catch_error
         use mod_string, only : expand_tilde
      implicit none
         character(*) :: &
            opt
         logical :: &
            ex
         !***
         select case(opt)
            case('file')
               inquire(file=expand_tilde(xyz_file),exist=ex)
               call catch_error(&
                        err=.not.ex,&
                        msg='file "'//trim(xyz_file)//&
                            '" does not exist.',&
                        proc='geo%read_carts')
            case('composition')
               call catch_error(&
                        err=size(symbs).ne.size(refframe%symbs),&
                        msg='the atom count of "start_geo" is '//&
                            'inconsistent with that of the '//&
                            'reference frame.',&
                        proc='geo%read_carts')
               call catch_error(&
                        err=any(symbs.ne.refframe%symbs),&
                        msg='the atomic symbols of "start_geo" '//&
                            'are inconsistent with those of the '//&
                            'reference frame.',&
                        proc='geo%read_carts')
            case default
               call catch_error(&
                        err=.true.,&
                        msg='unknown catch option "'//opt//'".',&
                        proc='geo%read_carts')
         end select
         !***
         end
      !*** 
      end
!***********************************************************************
      subroutine print_summary(this,npad)
!***********************************************************************
      use mod_string, only : to_string
      implicit none
      class(geo), target :: &
         this
      integer :: &
         npad
      logical, allocatable :: &
         fix_rcoors(:)
      character(:), allocatable :: &
         clead
      !***
      allocate(character(npad)::clead)
      clead(:npad)=''
      fix_rcoors=this%get_fix_rcoors('reparam')
      call print_energy()
      call print_pgname()
      call print_saddle_order()
      call print_xyz()
      call print_diag_hess()
      !***
      contains
      !*****************************************************************
         subroutine print_energy()
      !*****************************************************************
         implicit none
         !***
         if(this%ener.ne.huge(0.d0)) &
            write(*,'(/a,f0.3)') clead//'Energy: ',this%ener
         !***
         end
      !*****************************************************************
         subroutine print_pgname()
      !*****************************************************************
         implicit none
         character(:), allocatable :: &
            pgname
         !***
         pgname=this%get_pgroup_name()
         write(*,'(/a)') clead//'Point group: '//pgname
         !***
         end
      !*****************************************************************
         subroutine print_saddle_order()
      !*****************************************************************
         implicit none
         logical :: &
            flat_flag
         integer :: &
            order
         !***
         associate(t => this)
            if(.not.allocated(t%hess_eigvals)) return
            order=t%get_saddle_order()
            write(*,'(/a$)') &
                    clead//'Saddle order: '//to_string(order)
            if(allocated(t%sens_hess_eigvals)) then
               flat_flag=t%get_flat_flag()
               if(flat_flag) write(*,'(a$)') ' (FLAT)'
            endif
            write(*,*)
         end associate
         !***
         end
      !*****************************************************************
         subroutine print_xyz()
      !*****************************************************************
         implicit none
         integer :: &
            k
         !***
         write(*,'(/a/)') clead//'XYZ coordinates:'
         do k=1,layout%natoms
            write(*,'(18x,a,3f16.10)') &
               refframe%symbs(k),this%carts(:,k)
         end do
         !***
         end
      !*****************************************************************
         subroutine print_diag_hess()
      !*****************************************************************
         implicit none
         character(:), allocatable :: &
            rowlab
         integer :: &
            ndashes,&
            lsize
         !***
         associate(t        => this,&
                   nrcoors  => landscape%nrcoors,&
                   labs     => landscape%lab_rcoors)
            if(.not.allocated(t%hess_eigvals)) return
            !***
            joint: block
               lsize=max(len(labs),10)
               rowlab=repeat(' ',lsize)
            end block joint
            !***
            header: block
               if(allocated(t%sens_hess_eigvals)) then
                  write(*,'(/a)') &
                        clead//'Hessian eigenvalues & eigenvectors '//&
                               'with sensitivities '
                  write(*,'(a/)') &
                        clead//'(fixed coordinates with stars):'
               else
                  write(*,'(/a/)') &
                     clead//'Hessian eigenvalues & eigenvectors '//&
                            '(fixed coordinates with stars):'
               endif
            end block header
            !***
            rows: block
               !***
               integer :: k,kstart,kend
               do kstart=1,nrcoors,5
                  !***
                  init: block
                     kend=min(kstart+4,nrcoors)
                     ndashes=lsize+12*(kend-kstart)+12
                  end block init
                  !***
                  eigvals: block
                     write(rowlab,'(a)') 'eigval:'
                     write(*,'(3x,a,*(f12.4))') &
                        clead//rowlab,t%hess_eigvals(kstart:kend)
                     if(allocated(t%sens_hess_eigvals)) then
                        write(rowlab,'(a)') 'sens:'
                        write(*,'(3x,a,*(f12.4))') &
                           clead//rowlab,&
                           t%sens_hess_eigvals(kstart:kend)
                     endif
                     write(*,'(3x,a)') clead//repeat('-',ndashes)
                  end block eigvals
                  !***
                  coords: block
                     character(12) :: cmode
                     write(rowlab,'(a)') 'coord'
                     write(*,'(3x,a$)') clead//rowlab
                     do k=kstart,kend
                        write(cmode,'(a,i0)') 'eigvec #',k
                        cmode=adjustr(cmode)
                        write(*,'(a$)') cmode
                     end do
                     write(*,*)
                     write(*,'(3x,a)') clead//repeat('-',ndashes)
                  end block coords
                  !***
                  eigvects: block
                     integer :: l
                     do l=1,nrcoors
                        write(rowlab,'(a$)') &
                           trim(labs(l))//&
                           merge('*',' ',fix_rcoors(l))
                        write(*,'(3x,a,*(f12.4))') &
                           clead//rowlab,(t%hess_eigvects(l,k),&
                           k=kstart,kend)
                     end do
                     write(*,*)
                  end block eigvects
                  !***
               end do
               !***
            end block rows
            !***
         end associate
         !***
         end
      !***
      end
