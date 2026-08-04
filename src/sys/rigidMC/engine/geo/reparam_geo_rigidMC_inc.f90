!***********************************************************************
      subroutine reparametrize(this,mode)
!***********************************************************************
      use mod_catch, only : catch_error
      implicit none
      class(geo), target :: &
         this
      character(*) :: &
         mode
      !***
      associate(t=>this)
         select case(mode)
            case('run')
               call t%adapt_anchors()
               call t%reorient()
               call t%carts2rcoors()
               call t%adapt_axes()
               call t%freeze()
            case('assert_clean')
               call catch_local_errors('alloc_axes')
               call catch_local_errors('alloc_anchors')
            case('clean')
               call catch_local_errors('dealloc_axes')
               call catch_local_errors('dealloc_anchors')
               deallocate(t%axes,t%anchors)
               call t%reorient()
               call t%carts2rcoors()
               call t%freeze()
            case default
               call catch_local_errors('unknown')
         end select
      end associate
      !***
      contains
      !*****************************************************************
         subroutine catch_local_errors(opt)
      !*****************************************************************
         implicit none
         character(*) :: &
            opt
         !***
         associate(t=>this)
            select case(opt)
               case('alloc_axes')
                  call catch_error(&
                           err=allocated(t%axes),&
                           msg='array "axes" must be unallocated',&
                           proc='geo%reparametrize')
               case('alloc_anchors')
                  call catch_error(&
                           err=allocated(t%anchors),&
                           msg='array "anchors" must be unallocated',&
                           proc='geo%reparametrize')
               case('dealloc_axes')
                  call catch_error(&
                           err=.not.allocated(t%axes),&
                           msg='attempt to deallocate unallocated '//&
                               'array "axes".',&
                           proc='geo%reparametrize')
               case('dealloc_anchors')
                  call catch_error(&
                           err=.not.allocated(t%axes),&
                           msg='attempt to deallocate unallocated '//&
                               'array "anchors".',&
                           proc='geo%reparametrize')
               case('unknown')
                  call catch_error(&
                           err=.true.,&
                           msg='unknown mode "'//mode//'".',&
                           proc='geo%reparametrize')
               case default
                  call catch_error(&
                           err=.true.,&
                           msg='unknown catch option "'//opt//'".',&
                           proc='geo%reparametrize')
            end select
         end associate
         !***
         end
      end
!***********************************************************************
      subroutine adapt_axes(this)
!***********************************************************************
      use mod_euclid, only : gen_rigid_trfo
      implicit none
      class(geo), target :: &
         this
      character(:), allocatable :: &
         axclass
      real(r64), pointer :: &
         axes(:,:)
      type(monbox) :: &
         mb
      integer :: &
         m
      !***
      associate(t        => this,&
                nmons    => layout%nmons,&
                sing_tol => sentinel%sing_tol,&
                rig_tol  => sentinel%rig_tol)
         t%axes=spread(t%get_axes('standard'),3,nmons)
         do m=1,nmons
            axes=>t%axes(:,:,m)
            call mb%init(m,t)
            associate(alt_rot  => .not.t%get_fix_flag('alpha',m),&
                      alt_tran => .not.t%get_fix_flag('phi',m))
               axclass='alt:'
               if(alt_rot) axclass=axclass//'u'
               if(alt_tran) axclass=axclass//'v'
               if(axclass.ne.'alt:') &
                  axes=t%get_axes(axclass)
               call gen_rigid_trfo('carts2rcoors',axes,mb%rcoors,&
                        mb%ref_carts,mb%carts,sing_tol,rig_tol)
               if(alt_rot) &
                  call update(ax1=1,ang1=mb%alpha,&
                              ax2=2,ang2=mb%beta)
               if(alt_tran) &
                  call update(ax1=4,ang1=mb%phi,&
                              ax2=5,ang2=mb%theta)
            end associate
         end do
      end associate
      !***
      contains
      !*****************************************************************
         subroutine update(ax1,ang1,ax2,ang2)
      !*****************************************************************
         use mod_swap, only : swap
         use mod_euclid, only : get_rot_mat, std_axis_angle
         integer :: &
            ax1,ax2
         real(r64) :: &
            ang1,ang2
         !***
         associate(i => [1.d0,0.d0,0.d0],&
                   j => [0.d0,1.d0,0.d0])
            call swap(ang1,ang2)
            axes(:,ax1)=matmul(get_rot_mat(j,ang2),i)
            call std_axis_angle(axes(:,ax1),ang1)
            axes(:,ax2)=j
         end associate
         !***
         end
      !***
      end
!***********************************************************************
      subroutine adapt_anchors(this)
!***********************************************************************
      implicit none
      class(geo), target :: &
         this
      real(r64), target, allocatable :: &
         rot_list(:,:,:),&
         tran_list(:,:)
      !***
      associate(t         => this,&
                max_nedof => refframe%max_nedof)
         if(max_nedof.eq.6) then
            t%anchors=t%get_anchors()
         else
            call set_rigid_arrays()
            call find_best_anchors()
         endif
      end associate
      !***
      contains
      !*****************************************************************
         subroutine find_best_anchors()
      !*****************************************************************
         implicit none
         integer, allocatable :: &
            anchors(:)
         real(r64), allocatable :: & 
            cand_smdets(:)
         integer, allocatable :: &
            cand_anchors(:,:)
         real(r64) :: &
            max_smdet
         integer :: &
            ncands,i,ibest,&
            a1,a2,a3,a1ini,a1c
         !***
         associate(t          => this,&
                   nmons      => layout%nmons,&
                   max_nedof  => refframe%max_nedof,&
                   nedofs     => refframe%nedofs,&
                   rsmdet_tol => sentinel%rsmdet_tol,&
                   mcands     => layout%nmons**3+1)
            cand_smdets=spread(0.d0,1,mcands)
            cand_anchors=spread([0,0,0],2,mcands)
            ncands=0
            do a3=nmons,1,-1
               do a2=nmons,1,-1
                  if(a2.eq.a3) cycle
                  if(nedofs(a2).gt.nedofs(a3)) cycle
                  a1ini=merge(0,1,max_nedof.eq.5)
                  do a1=nmons,a1ini,-1
                     if(any(a1.eq.[a2,a3])) cycle
                     a1c=max(a1,1)
                     if(nedofs(a1c).gt.nedofs(a2)) cycle
                     ncands=ncands+1
                     anchors=[a1,a2,a3]
                     cand_anchors(:,ncands)=anchors
                     cand_smdets(ncands)=get_smdet(anchors)
                     if(max_nedof.eq.5) exit
                  end do
               end do
            end do
            ibest=0
            if(ncands.gt.0) then
               cand_smdets=cand_smdets(:ncands)
               cand_anchors=cand_anchors(:,:ncands)
               max_smdet=maxval(cand_smdets)
               do i=1,ncands
                  if(cand_smdets(i).gt.rsmdet_tol*max_smdet) then
                     ibest=i
                     exit
                  endif
               end do
            endif
            if(ibest.eq.0) then
               t%anchors=t%get_anchors()
            else
               t%anchors=cand_anchors(:,ibest)
            endif
         end associate
         !***
         end
      !*****************************************************************
         function get_smdet(anchors) result(smdet)
      !*****************************************************************
         use mod_euclid, only : gen_spherical_decomp, cross_product,&
                                normalize
         implicit none
         real(r64), allocatable :: &
            axes(:,:)
         integer :: &
            anchors(3)
         real(r64) :: &
            smdet,&
            tran(3),&
            phi,theta,rho
         real(r64), pointer :: &
            rot(:,:),&
            tran3(:),&
            old_tran(:)
         !***
         associate(max_nedof => refframe%max_nedof,&
                   anc1      => anchors(1),&
                   anc2      => anchors(2),&
                   anc3      => anchors(3),&
                   sing_tol  => sentinel%sing_tol,&
                   deg2rad   => dacos(-1.d0)/180.d0)
            tran3=>tran_list(:,anc3)
            select case(max_nedof)
               case(5)
                  rot=>rot_list(:,:,anc3)
                  old_tran=>tran_list(:,anc2)
                  tran=matmul(transpose(rot),old_tran-tran3)
                  axes=this%get_axes('standard')
                  call gen_spherical_decomp(tran,axes(:,:3),&
                           phi,theta,rho)
                  smdet=rho**4*dsin(deg2rad*theta)**2
               case(3)
                  associate(v => tran_list(:,anc1)-tran3,&
                            u => normalize(tran_list(:,anc2)-tran3,&
                                    sing_tol))
                     smdet=dot_product(v,v)*&
                           dot_product(cross_product(v,u),&
                                       cross_product(v,u))
                  end associate
            end select
         end associate
         end
      !*****************************************************************
         subroutine set_rigid_arrays()
      !*****************************************************************
         use mod_euclid, only : gen_rigid_kernel, gen_rigid_trfo
         implicit none
         real(r64), allocatable :: &
            axes(:,:)
         type(monbox) :: &
            mb
         integer :: &
            m
         !***
         associate(nmons=>layout%nmons)
            allocate(rot_list(3,3,nmons))
            allocate(tran_list(3,nmons))
            axes=this%get_axes('standard')
            do m=1,nmons
               call mb%init(m,this)
               call gen_rigid_trfo(&
                        'carts2rcoors',&
                        axes,mb%rcoors,&
                        mb%ref_carts,&
                        mb%carts,&
                        sentinel%sing_tol,&
                        sentinel%rig_tol)
               call gen_rigid_kernel(&
                        axes=axes,&
                        rcoors=mb%rcoors,&
                        rot=rot_list(:,:,m),&
                        tran=tran_list(:,m))
            end do
         end associate
         !***
         end
      !***
      end
