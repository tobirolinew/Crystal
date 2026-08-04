!***********************************************************************
      function get_abnorm_flag(this) result(abnorm_flag)
!***********************************************************************
      use mod_graph, only : is_connected
      use mod_pertab, only : vdw_rads, cov_rads
      implicit none
      class(geo) :: &
         this
      integer :: &
         abnorm_flag
      type(monbox) :: &
         mb1,mb2
      integer, allocatable :: &
         C(:,:)
      real(r64), allocatable :: &
         D(:,:),&
         Scov(:,:),&
         Svdw(:,:)
      integer :: &
         m1,m2
      !*** 
      associate(t       => this,&
                nmons   => layout%nmons,&
                lam_cov => sentinel%lam_cov,&
                lam_vdw => sentinel%lam_vdw)
         allocate(C(nmons,nmons),source=0)
         abnorm_flag=-1
         do m1=1,nmons
            call mb1%init(m1,t)
            associate(carts1 => mb1%carts,&
                      cov1   => cov_rads(mb1%anums),&
                      vdw1   => vdw_rads(mb1%anums),&
                      n1     => mb1%natoms)
               do m2=1,m1-1
                  call mb2%init(m2,t)
                  associate(carts2 => mb2%carts,&
                            cov2   => cov_rads(mb2%anums),&
                            vdw2   => vdw_rads(mb2%anums),&
                            n2     => mb2%natoms)
                     D=norm2(spread(carts1,3,n2)-&
                             spread(carts2,2,n1),dim=1)
                     Scov=spread(cov1,2,n2)+spread(cov2,1,n1)
                     Svdw=spread(vdw1,2,n2)+spread(vdw2,1,n1)
                     if(any(D.lt.lam_cov*Scov)) return
                     if(any(D.le.lam_vdw*Svdw)) then
                        C(m1,m2)=1
                        C(m2,m1)=1
                     endif
                  end associate
               end do
            end associate
         end do
         abnorm_flag=merge(0,1,is_connected(C))
      end associate
      !*** 
      end
!***********************************************************************
      function get_pgroup_name(this) result(pgname)
!***********************************************************************
      use mod_point_group, only : detect_point_group
      implicit none
      class(geo) :: &
         this
      real(r64), allocatable :: &
         tmp_carts(:,:)
      character(:), allocatable :: &
         pgname
      call catch_local_errors()
      !***
      associate(t        => this,&
                anums    => refframe%anums,&
                pgrp_tol => sentinel%pgrp_tol)
         tmp_carts=t%carts
         call detect_point_group(&
                  anums,tmp_carts,pgname,pgrp_tol=pgrp_tol)
      end associate
      !***
      contains
      !*****************************************************************
         subroutine catch_local_errors()
      !*****************************************************************
      use mod_catch, only: catch_error
         implicit none
         !***
         associate(t      => this,&
                   natoms => layout%natoms)
            call catch_error(&
                     err=.not.allocated(t%carts),&
                     msg='array "carts" is unallocated.',&
                     proc='geo%get_pgroup_name')
            call catch_error(&
                     err=any(shape(t%carts).ne.[3,natoms]),&
                     msg='array "carts" has incorrect shape.',&
                     proc='geo%get_pgroup_name')
         end associate
         !***
         end
      !***
      end
!***********************************************************************
      function get_flat_flag(this) result(flat_flag)
!***********************************************************************
      implicit none
      class(geo) :: &
         this
      logical :: &
         flat_flag
      integer :: &
         i
      !***
      call catch_local_errors()
      associate(t        => this,&
                nrcoors  => landscape%nrcoors,&
                flat_fac => sentinel%flat_fac)
         flat_flag=.false.
         do i=1,nrcoors
            associate(abs_eigval  => dabs(t%hess_eigvals(i)),&
                      sens_eigval => t%sens_hess_eigvals(i))
               if(abs_eigval.ne.0.d0.and.&
                  abs_eigval.lt.flat_fac*sens_eigval) then
                  flat_flag=.true.
                  exit
               endif
            end associate
         end do
      end associate
      !***
      contains
      !*****************************************************************
         subroutine catch_local_errors()
      !*****************************************************************
         use mod_catch, only: catch_error
         implicit none
         !***
         associate(t        => this,&
                   nrcoors  => landscape%nrcoors)
            call catch_error(&
                     err=.not.allocated(t%hess_eigvals),&
                     msg='array "hess_eigvals" is unallocated.',&
                     proc='geo%get_flat_flag')
            call catch_error(&
                     err=.not.allocated(t%sens_hess_eigvals),&
                     msg='array "sens_hess_eigvals" is unallocated.',&
                     proc='geo%get_flat_flag')
            call catch_error(&
                     err=any(nrcoors.ne.[size(t%hess_eigvals),&
                                         size(t%sens_hess_eigvals)]),&
                     msg='size mismatch among arrays "hess_eigvals" '//&
                         'and "sens_hess_eigvals"',&
                     proc='geo%get_flat_flag')
         end associate
         !***
         end
      !***
      end
!***********************************************************************
      function get_saddle_order(this) result(order)
!***********************************************************************
      implicit none
      class(geo) :: &
         this
      integer :: &
         order
      !***
      call catch_local_errors()
      associate(t        => this,&
                flat_fac => sentinel%flat_fac)
         order=count(t%hess_eigvals.lt.-flat_fac*t%sens_hess_eigvals)
      end associate 
      !***
      contains
      !*****************************************************************
         subroutine catch_local_errors()
      !*****************************************************************
         use mod_catch, only: catch_error
         implicit none
         !***
         associate(t       => this,&
                   nrcoors => landscape%nrcoors)
            call catch_error(&
                     err=.not.allocated(t%hess_eigvals),&
                     msg='array "hess_eigvals" is not allocated.',&
                     proc='geo%get_saddle_order')
            call catch_error(&
                     err=size(t%hess_eigvals).ne.nrcoors,&
                     msg='array "hess_eigvals" has incorrect size.',&
                     proc='geo%get_saddle_order')
            call catch_error(&
                     err=.not.allocated(t%hess_eigvects),&
                     msg='array "hess_eigvects" is not allocated.',&
                     proc='geo%get_saddle_order')
            call catch_error(&
                     err=any(nrcoors.ne.shape(t%hess_eigvects)),&
                     msg='array "hess_eigvects" has incorrect size.',&
                     proc='geo%get_saddle_order')
            call catch_error(&
                     err=.not.allocated(t%sens_hess_eigvals),&
                     msg='array "sens_hess_eigvals" is not allocated.',&
                     proc='geo%get_saddle_order')
            call catch_error(&
                     err=size(t%sens_hess_eigvals).ne.nrcoors,&
                     msg='array "sens_hess_eigvals" has '//&
                         'incorrect size.',&
                     proc='geo%get_saddle_order')
         end associate
         !***
         end
      !***
      end
!***********************************************************************
      function get_grmsdLB(this,other) result(grmsdLB)
!***********************************************************************
      implicit none
      class(geo) :: &
         this,&
         other
      real(r64) :: &
         grmsdLB
      integer :: &
         npairs
      !***
      if(.not.allocated(this%sort_dists)) then
         write(*,*) 'WTF???????????????????????????????'
         stop
      endif
      if(.not.allocated(other%sort_dists)) then
         write(*,*) 'WTF???????????????????????????????'
         stop
      endif
      associate(sdistsA => this%sort_dists,&
                sdistsB => other%sort_dists,&
                natoms  => layout%natoms)
         npairs=natoms*(natoms-1)/2
         grmsdLB=norm2(sdistsA-sdistsB)/dsqrt(dble(npairs))
      end associate
      !***
      end
!***********************************************************************
      function get_sort_dists(this) result(sort_dists)
!***********************************************************************
      use mod_sort, only: get_num_lex_sort_perm
      implicit none
      class(geo) :: &
         this
      real(r64), allocatable :: &
         sort_dists(:),&
         tmp(:,:)
      integer, allocatable :: &
         perm(:)
      integer :: &
         npairs,&
         max_mon_bin,&
         max_anum,&
         i,&
         j,&
         k
      !***
      associate(t        => this,&
                anums    => refframe%anums,&
                mon_bins => layout%mon_bins,&
                mon_inds => layout%mon_inds,&
                natoms   => layout%natoms)
         npairs=natoms*(natoms-1)/2
         max_mon_bin=maxval(mon_bins)
         max_anum=maxval(anums)
         tmp=spread([0.d0,0.d0,0.d0],2,npairs)
         k=0
         do i=1,natoms-1
            associate(mi => mon_bins(mon_inds(i)),&
                      ai => anums(i))
               do j=i+1,natoms
                  k=k+1
                  associate(mj => mon_bins(mon_inds(j)),&
                            aj => anums(j))
                     tmp(1,k)=max_mon_bin*(min(mi,mj)-1)+max(mi,mj)
                     tmp(2,k)=max_anum*(min(ai,aj)-1)+max(ai,aj)
                     tmp(3,k)=norm2(t%carts(:,i)-t%carts(:,j))
                  end associate
               end do
            end associate
         end do
         perm=get_num_lex_sort_perm(tmp,asc=.true.)
         sort_dists=tmp(3,perm)
      end associate
      !***
      end
!***********************************************************************
      function get_jac_carts(this,fixed) result(jac)
!***********************************************************************
      use mod_euclid, only : gen_rigid_jac
      implicit none
      class(geo) :: &
         this
      logical :: &
         fixed(:)
      real(r64), allocatable :: &
         jac(:,:,:)
      type(monbox) :: &
         mb
      real(r64), allocatable :: &
         axes(:,:)
      integer :: &
         k,n,m
      !***
      associate(nrcoors => landscape%nrcoors,&
                natoms  => layout%natoms)
         jac=spread(spread(spread(0.d0,1,3),2,natoms),3,nrcoors)
         do k=1,nrcoors
            if(fixed(k)) cycle
            m=(k-1)/6+1
            n=k-6*(m-1)
            axes=this%get_axes('dynamic',m)
            call mb%init(m,this)
            call gen_rigid_jac(axes,mb%rcoors,mb%ref_carts,n,&
                               jac(:,mb%first:mb%last,k))
         end do
      end associate
      !***
      end
!***********************************************************************
      function get_gmat(this,fixed) result(gmat)
!***********************************************************************
      use mod_catch, only : catch_error
      implicit none
      class(geo) :: &
         this
      logical :: &
         fixed(:)
      real(r64), allocatable :: &
         gmat(:,:),&
         jac(:,:,:)
      real(r64) :: &
         tr
      integer :: &
         k,l,nfree
      !***
      call this%rcoors2carts()
      jac=this%get_jac_carts(fixed)
      associate(nrcoors=>landscape%nrcoors)
         gmat=spread(spread(0.d0,1,nrcoors),2,nrcoors)
         do k=1,nrcoors
            if(fixed(k)) cycle
            do l=k,nrcoors
               if(fixed(l)) cycle
               gmat(k,l)=sum(jac(:,:,k)*jac(:,:,l))
               gmat(l,k)=gmat(k,l)
            end do
         end do
         nfree=count(.not.fixed)
         call catch_error(&
                  err=nfree.eq.0,&
                  msg='there are no free coordinates.',&
                  proc='geo%get_gmat')
         tr=sum([(gmat(k,k),k=1,nrcoors)],mask=.not.fixed)
         gmat=gmat*dble(nfree)/tr
      end associate
      !***
      end
