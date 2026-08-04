!***********************************************************************
      subroutine rcoors2carts(this)
!***********************************************************************
      use mod_euclid, only : gen_rigid_trfo
      implicit none
      class(geo) :: &
         this
      type(monbox) :: &
         mb
      real(r64), allocatable :: &
         axes(:,:)
      integer :: &
         m
      !***
      this%carts=spread([0.d0,0.d0,0.d0],2,layout%natoms)
      do m=1,layout%nmons
         axes=this%get_axes('dynamic',m)
         call mb%init(m,this)
         call gen_rigid_trfo('rcoors2carts',axes,mb%rcoors,&
                  mb%ref_carts,mb%carts)
      end do
      !***
      end
!***********************************************************************
      subroutine carts2rcoors(this)
!***********************************************************************
      use mod_euclid, only : gen_rigid_trfo
      implicit none
      class(geo) :: &
         this
      type(monbox) :: &
         mb
      real(r64), allocatable :: &
         axes(:,:)
      integer :: &
         m
      !***
      associate(t        => this,&
                nrcoors  => landscape%nrcoors,&
                nmons    => layout%nmons,&
                sing_tol => sentinel%sing_tol,&
                rig_tol  => sentinel%rig_tol)
         t%rcoors=spread(0.d0,1,nrcoors)
         do m=1,nmons
            axes=t%get_axes('dynamic',m)
            call mb%init(m,this)
            call gen_rigid_trfo('carts2rcoors',axes,mb%rcoors,&
                     mb%ref_carts,mb%carts,sing_tol,rig_tol)
         end do
      end associate
      !***
      end
!***********************************************************************
      subroutine reorient(this)
!***********************************************************************
      use mod_euclid, only : gen_rigid_trfo
      implicit none
      class(geo) :: &
         this
      integer, allocatable :: &
         anchors(:)
      real(r64), allocatable :: &
         axes(:,:)
      logical :: &
         collin
      !***
      associate(t=>this)
         anchors=t%get_anchors()
         axes=t%get_axes('standard')
         collin=t%get_collin_flag()
      end associate
      associate(nmons     => layout%nmons,&
                max_nedof => refframe%max_nedof,&
                alast     => anchors(3),&
                a2last    => anchors(2),&
                a3last    => anchors(1))
         call remove(alast,'phi')
         call remove(alast,'theta')
         call remove(alast,'rho')
         call remove(alast,'alpha')
         call remove(alast,'beta')
         call remove(alast,'gamma')
         select case(max_nedof)
            case(5)
               if(.not.collin) call remove(a2last,'phi')
            case(3)
               call remove(a2last,'phi')
               call remove(a2last,'theta')
               if(a3last.ne.0.and..not.collin) &
                  call remove(a3last,'phi')
         end select
      end associate
      if(collin) call align_collin()
      contains
      !*****************************************************************
         subroutine align_collin()
      !*****************************************************************
         use mod_pertab, only : iso_masses
         use mod_euclid, only : analyze_shape
         implicit none
         type(monbox) :: &
            mb
         real(r64) :: &
            eigvects(3,3),&
            com(3)
         integer :: &
            m
         !***
         associate(t          => this,&
                   nmons      => layout%nmons,&
                   collin_tol => sentinel%collin_tol)
            call analyze_shape(t%carts,eigvects=eigvects,&
                     sing_tol=collin_tol)
            t%carts=matmul(transpose(eigvects),t%carts)
            do m=1,nmons
               call mb%init(m,t)
               associate(carts  => mb%carts,&
                         masses => iso_masses(mb%anums),&
                         natoms => mb%natoms)
                  com=matmul(carts,masses)/sum(masses)
                  carts=carts-spread(com,2,natoms)
                  call analyze_shape(carts,eigvects=eigvects,&
                           sing_tol=collin_tol)
                  carts=matmul(transpose(eigvects),carts)
                  carts=carts+spread([0.d0,0.d0,&
                           dsign(norm2(com),com(3))],2,natoms)
               end associate
            end do
         end associate
         !***
         end
      !*****************************************************************
         subroutine remove(m,rcname)
      !*****************************************************************
         use mod_euclid, only : gen_rigid_trfo
         implicit none
         character(*) :: &
            rcname
         type(monbox) :: &
            mb
         real(r64) :: &
            inv_rcoors(6)
         real(r64), allocatable :: &
            inv_ref_carts(:,:)
         integer :: &
            m
         !***
         if(m.le.0) return
         call mb%init(m,this)
         call gen_rigid_trfo('carts2rcoors',axes,rcoors=mb%rcoors,&
                  ref_carts=mb%ref_carts,carts=mb%carts,&
                  tol=sentinel%sing_tol,rig_tol=sentinel%rig_tol)
         select case(rcname)
         case('alpha')
            inv_rcoors=[0.d0,0.d0,-mb%alpha,0.d0,0.d0,0.d0]
         case('beta')
            inv_rcoors=[0.d0,-mb%beta,0.d0,0.d0,0.d0,0.d0]
         case('gamma')
            inv_rcoors=[0.d0,0.d0,-mb%alpha-mb%gamma,0.d0,0.d0,0.d0]
         case('phi')
            inv_rcoors=[0.d0,0.d0,-mb%phi,0.d0,0.d0,0.d0]
         case('theta')
            inv_rcoors=[0.d0,-mb%theta,0.d0,0.d0,0.d0,0.d0]
         case('rho')
            inv_rcoors=[0.d0,0.d0,0.d0,0.d0,0.d0,-mb%rho]
         case default
            call catch_local_errors('rcname',rcname)
         end select
         inv_ref_carts=this%carts
         call gen_rigid_trfo('rcoors2carts',axes,rcoors=inv_rcoors,&
                  ref_carts=inv_ref_carts,carts=this%carts)
         end
      !*****************************************************************
         subroutine catch_local_errors(opt,rcname)
      !*****************************************************************
         use mod_catch, only : catch_error
         implicit none
         character(*) :: &
            opt,rcname
         select case(opt)
            case('rcname')
               call catch_error(&
                        err=.true.,&
                        msg='incorrect rcname value "'//rcname//'".',&
                        proc='geo%reorient')
            case('rhofix')
            case default
               call catch_error(&
                        err=.true.,&
                        msg='unknown option "'//opt//'".',&
                        proc='geo%reorient')
         end select
         end
      !***
      end
!***********************************************************************
      subroutine freeze(this)
!***********************************************************************
      use mod_catch, only : catch_error
      implicit none
      class(geo) :: &
         this
      integer, allocatable :: &
         anchors(:)
      logical, allocatable :: &
         fix_rcoors(:)
      integer :: &
         m
      !***
      anchors=this%get_anchors()
      fix_rcoors=this%get_fix_rcoors('normal')
      do m=1,layout%nmons
         if(m.eq.anchors(3)) cycle
         call catch_error(&
                  err=fix_rcoors(6*m),&
                  msg='only the last radial coordinate '//&
                      'can be fixed in module "rigidMC".',&
                  proc='geo%freeze')
      end do
      where(fix_rcoors)
         this%rcoors=0.d0
      end where
      !***
      end
!***********************************************************************
      subroutine renormalize(this)
!***********************************************************************
      use mod_catch, only : catch_error
      use mod_euclid, only : wrap_angle
      implicit none
      class(geo) :: &
         this
      type(monbox) :: &
         mb
      integer :: &
         m
      !***
      call this%reparametrize('assert_clean')
      do m=1,layout%nmons
         call mb%init(m,this)
         associate(alpha => mb%alpha,&
                   beta  => mb%beta,&
                   gamma => mb%gamma,&
                   phi   => mb%phi,&
                   theta => mb%theta,&
                   rho   => mb%rho)
            call wrap_angle(beta)
            if(beta.lt.0.d0) then
               beta=-beta
               alpha=alpha+180.d0
               gamma=gamma+180.d0
            endif
            if(rho.lt.0.d0) then
               rho=-rho
               theta=180.d0-theta
               phi=phi+180.d0
            endif
            call wrap_angle(theta)
            if(theta.lt.0.d0) then
               theta=-theta
               phi=phi+180.d0
            endif
            call wrap_angle(alpha)
            call wrap_angle(gamma)
            call wrap_angle(phi)
         end associate
      end do
      !***
      end
!***********************************************************************
      subroutine align(this,other,update,grmsd,grmsd_tol)
!***********************************************************************
      use mod_graph, only: bipartite_matching
      use mod_rbalign, only: config_rbalign, best_align
      implicit none
      type :: tier
         integer, allocatable :: &
            links(:,:),&
            matches(:)
         real(r64), allocatable :: &
            weights(:)
         real(r64) :: &
            cost
      end type
      class(geo), target :: &
         this,&
         other
      logical :: &
         update
      real(r64), pointer :: &
         cartsB(:,:)
      real(r64), optional :: &
         grmsd_tol
      real(r64) :: &
         grmsd
      !***
      call config_rbalign(matchtool=bipartite_matching)
      if(present(grmsd_tol)) &
         call config_rbalign(grmsd_tol=grmsd_tol)
      associate(cartsA => other%carts,&
                anums  => refframe%anums,&
                natoms => layout%natoms)
         if(update) then
            cartsB=>this%carts
         else
            allocate(cartsB(3,natoms))
            cartsB=this%carts
         endif
         call best_align(anums,cartsA,anums,cartsB,grmsd)
         deallocate(cartsB)
      end associate
      !***
      end
