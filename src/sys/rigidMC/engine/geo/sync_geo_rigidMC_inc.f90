!***********************************************************************
      function get_collin_flag(this) result(collin_flag)
!***********************************************************************
      use mod_euclid, only : analyze_shape
      implicit none
      class(geo) :: &
         this
      logical :: &
         collin_flag
      integer :: &
         rank
      !***
      associate(t          => this,&
                collin_tol => sentinel%collin_tol)
         call analyze_shape(t%carts,rank,sing_tol=collin_tol)
         collin_flag=rank.eq.1
      end associate
      !***
      end
!***********************************************************************
      function get_anchors(this) result(anchors)
!***********************************************************************
      implicit none
      class(geo) :: &
         this
      real(r64) :: &
         anchors(3)
      integer :: &
         m
      !***
      if(allocated(this%anchors)) then
         anchors=this%anchors
      else
         associate(nmons     => layout%nmons,&
                   nedofs    => refframe%nedofs,&
                   inds      => [(m,m=1,layout%nmons)],&
                   a1        => anchors(1),&
                   a2        => anchors(2),&
                   a3        => anchors(3))
            a3=maxloc(nedofs,dim=1)
            a2=maxloc(nedofs,mask=inds.ne.a3,dim=1)
            a1=findloc(inds.ne.a2.and.inds.ne.a3,.true.,dim=1)
         end associate
      endif
      !***
      end
!***********************************************************************
      function get_axes(this,axclass,m) result(axes)
!***********************************************************************
      use mod_catch, only : catch_error
      implicit none
      class(geo) :: &
         this
      character(*) :: &
         axclass
      real(r64), allocatable :: &
         axes(:,:)
      integer, optional :: &
         m
      !***
      associate(t => this,&
                i => [1.d0,0.d0,0.d0],&
                j => [0.d0,1.d0,0.d0],&
                k => [0.d0,0.d0,1.d0])
         select case(axclass)
            case('standard')
               axes=reshape([k,j,k,k,j,k],[3,6])
            case('alt:u')
               axes=reshape([j,i,k,k,j,k],[3,6])
            case('alt:v')
               axes=reshape([k,j,k,j,i,k],[3,6])
            case('alt:uv')
               axes=reshape([j,i,k,j,i,k],[3,6])
            case('dynamic')
               call catch_error(&
                        err=.not.present(m),&
                        msg='monomer index m is required '//&
                            'for axis class "dynamic".',&
                        proc='geo%get_axes')
               if(allocated(t%axes)) then
                  axes=t%axes(:,:,m)
               else
                  axes=reshape([k,j,k,k,j,k],[3,6])
               endif
            case default
               call catch_error(&
                        err=.true.,&
                        msg='unknown axis class "'//axclass//'".',&
                        proc='geo%get_axes')
         end select
      end associate
      !***
      end
!***********************************************************************
      recursive function get_fix_rcoors(this,mode) result(fix_rcoors)
!***********************************************************************
      use mod_catch, only: catch_error
      implicit none
      class(geo) :: &
         this
      character(*) :: &
         mode
      logical :: &
         fix_flag
      logical, allocatable :: &
         fix_rcoors(:)
      integer :: &
         m
      !***
      associate(t         => this,&
                nmons     => layout%nmons,&
                nrcoors   => landscape%nrcoors,&
                rad_relax => control%rad_relax)
         select case(mode)
            case('normal')
               do m=1,nmons
                  fix_flag=t%get_fix_flag('alpha',m)
                  call append_logical(fix_rcoors,fix_flag)
                  fix_flag=t%get_fix_flag('beta',m)
                  call append_logical(fix_rcoors,fix_flag)
                  fix_flag=t%get_fix_flag('gamma',m)
                  call append_logical(fix_rcoors,fix_flag)
                  fix_flag=t%get_fix_flag('phi',m)
                  call append_logical(fix_rcoors,fix_flag)
                  fix_flag=t%get_fix_flag('theta',m)
                  call append_logical(fix_rcoors,fix_flag)
                  fix_flag=t%get_fix_flag('rho',m)
                  call append_logical(fix_rcoors,fix_flag)
               end do
            case('rad_relax')
               fix_rcoors=spread(.true.,1,nrcoors)
               do m=1,nmons
                  associate(rho_idx=>6*m)
                     fix_rcoors(rho_idx)=t%get_fix_flag('rho',m)
                  end associate
               end do
            case('crystal')
               call t%reparametrize('assert_clean')
               fix_rcoors=t%get_fix_rcoors('normal')
               if(rad_relax) &
                  fix_rcoors=fix_rcoors.or.&
                        .not.t%get_fix_rcoors('rad_relax')
            case('reparam')
               call t%reparametrize('run')
               fix_rcoors=t%get_fix_rcoors('normal')
               call t%reparametrize('clean')
            case default
               call catch_error(&
                        err=.true.,&
                        msg='unknown mode "'//mode//'".',&
                        proc='geo%fix_rcoors')
         end select
      end associate
      !***
      contains
      !*****************************************************************
      !  subroutine append_logical()
      !*****************************************************************
#        define ARR_TYPE logical
#        define APPEND append_logical
#        include "../../../../util/dep/generic_append_inc.f90"
      end
!***********************************************************************
      function get_fix_flag(this,rclab,m) result(fix_flag)
!***********************************************************************
      use mod_catch, only : catch_error
      implicit none
      class(geo) :: &
         this
      character(*) :: &
         rclab
      logical :: &
         fix_flag,&
         collin
      integer, allocatable :: &
         anchors(:)
      integer :: &
         m,&
         ianch
      !***
      associate(t         => this,&
                nedof     => refframe%nedofs(m),&
                max_nedof => refframe%max_nedof)
         collin=t%get_collin_flag()
         anchors=t%get_anchors()
         ianch=findloc(anchors,m,dim=1)
         select case(rclab)
            case('alpha','beta')
               select case(ianch)
                  case(3)
                     fix_flag=.true.
                  case default
                     fix_flag=nedof.eq.3
               end select
            case('gamma')
               select case(ianch)
                  case(3)
                     fix_flag=.true.
                  case default
                     fix_flag=nedof.le.5
               end select
            case('phi')
               select case(ianch)
                  case(1)
                     if(collin) then
                        fix_flag=.false.
                     else
                        fix_flag=max_nedof.eq.3
                     endif
                  case(2)
                     if(collin) then
                        fix_flag=anchors(1).gt.0.or.nedof.eq.3
                     else
                        fix_flag=max_nedof.le.5
                     endif
                  case(3)
                     fix_flag=.true.
                  case default
                     fix_flag=.false.
               end select
            case('theta')
               select case(ianch)
                  case(3)
                     fix_flag=.true.
                  case(2)
                     fix_flag=max_nedof.eq.3
                  case default
                     fix_flag=.false.
               end select
            case('rho')
               select case(ianch)
                  case(3)
                     fix_flag=.true.
                  case default
                     fix_flag=.false.
               end select
            case default
               call catch_error(&
                        err=.true.,&
                        msg='unknown coordinate label "'//&
                             rclab//'".',&
                        proc='geo%fix_flag')
         end select
      end associate
      !***
      end
