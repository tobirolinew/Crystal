!***********************************************************************
      subroutine optimize(this,rad_relax,polish,dump,npad)
!***********************************************************************
      use mod_quasi_newton, only : &
         bfgs, newton_powell
      implicit none
      type lview_
         logical :: &
            rad_relax=.false.,&
            dump=.false.,&
            polish=.false.
         integer :: &
            npad=0
      end type
      type(lview_) :: &
         lview
      class(geo) :: &
         this
      logical, optional :: &
         dump,&
         rad_relax,&
         polish
      integer, optional :: &
         npad
      real(r64), allocatable :: &
         prev_axes(:,:,:),&
         hess0(:,:)
      logical, allocatable :: &
         fix_rcoors(:)
      character(:), allocatable :: &
         clead
      logical :: &
         conv_flag,&
         stable
      real(r64) :: &
         recond,&
         dev,&
         rgmax,&
         best_rgmax
      integer :: &
         nsharps,&
         i,k,m
      !***
      if(present(rad_relax)) lview%rad_relax=rad_relax
      if(present(polish)) lview%polish=polish
      if(present(dump)) lview%dump=dump
      if(present(npad)) lview%npad=npad
      associate(t            => this,&
                max_order    => control%max_order,&
                mreconds     => sentinel%mreconds,&
                recond_tol   => sentinel%recond_tol,&
                mbfgs_relax  => sentinel%mbfgs_relax,&
                mrest_relax  => sentinel%mrest_relax,&
                msharps      => sentinel%msharps,&
                mbfgs        => sentinel%mbfgs,&
                mnp          => sentinel%mnp,&
                gtol_rcoors  => landscape%gtol_rcoors,&
                scal_rcoors  => landscape%scal_rcoors,&
                dvtol_rcoors => landscape%dvtol_rcoors,&
                nmons        => layout%nmons)
         call report('header')
         if(lview%rad_relax) then
            call t%reparametrize('assert_clean')
            fix_rcoors=t%get_fix_rcoors('rad_relax')
            do i=1,mrest_relax
               call bfgs(func=pot,grad=wrap_grad,v=t%rcoors,&
                         f=t%ener,conv_flag=conv_flag,&
                         fixed=fix_rcoors,&
                         gtols=gtol_rcoors,&
                         dvtols=dvtol_rcoors,&
                         mbfgs=mbfgs_relax*(nmons-1),&
                         hess0=hess(t%rcoors,fix_rcoors,'forward'),&
                         scals=scal_rcoors,&
                         dump=lview%dump,npad=lview%npad)
               call t%rcoors2carts()
               if(t%get_abnorm_flag().ne.0) then
                  t%ener=huge(0.d0)
                  exit
               endif
               if(conv_flag) exit
            end do
         else
            prev_axes=spread(t%get_axes('standard'),3,nmons)
            conv_flag=.false.
            nsharps=0
            best_rgmax=huge(0.d0)
            do i=1,mreconds
               call t%rcoors2carts()
               call t%reparametrize('run')
               !***
               fix_rcoors=t%get_fix_rcoors('normal')
               associate(naxes => size(t%axes,2))
                  recond=0.d0
                  do m=1,nmons
                     do k=1,naxes
                        associate(axes  => t%axes(:,k,m),&
                                  paxes => prev_axes(:,k,m))
                           dev=min(norm2(axes-paxes),&
                                   norm2(axes+paxes))
                           recond=max(recond,dev)
                        end associate
                     end do
                  end do
                  stable=i.gt.1.and.recond.lt.recond_tol
               end associate
               if(stable.and.conv_flag) exit
               call report('recond')
               select case(max_order)
                  case(0)
                     hess0=t%get_gmat(fix_rcoors)
                     call bfgs(func=pot,grad=wrap_grad,&
                               v=t%rcoors,f=t%ener,&
                               conv_flag=conv_flag,&
                               fixed=fix_rcoors,&
                               gtols=gtol_rcoors,&
                               dvtols=dvtol_rcoors,&
                               hess0=hess0,&
                               mbfgs=mbfgs,&
                               scals=scal_rcoors,&
                               dump=lview%dump,&
                               rgmax=rgmax,npad=lview%npad)
                  case default
                     hess0=hess(t%rcoors,fix_rcoors,'forward')
                     call newton_powell(&
                               func=pot,grad=wrap_grad,&
                               v=t%rcoors,f=t%ener,&
                               conv_flag=conv_flag,&
                               fixed=fix_rcoors,&
                               gtols=gtol_rcoors,&
                               dvtols=dvtol_rcoors,&
                               mnp=mnp,&
                               hess0=hess0,&
                               stepcorr=dispcorr,&
                               rgmax=rgmax,&
                               dump=lview%dump,npad=lview%npad)
                     fix_rcoors=t%get_fix_rcoors('rad_relax')
                     hess0=t%get_gmat(fix_rcoors)
                     call bfgs(func=pot,grad=wrap_grad,&
                               v=t%rcoors,f=t%ener,&
                               fixed=fix_rcoors,&
                               gtols=gtol_rcoors,&
                               hess0=hess0,&
                               mbfgs=mbfgs_relax*(nmons-1),&
                               scals=scal_rcoors,&
                               dump=lview%dump,npad=lview%npad)
               end select
               call t%rcoors2carts()
               if(t%get_abnorm_flag().ne.0) then
                  t%ener=huge(0.d0)
                  conv_flag=.false.
                  exit
               endif
               if(rgmax.lt.0.9d0*best_rgmax) then
                  nsharps=0
               else
                  nsharps=nsharps+1
               endif
               best_rgmax=rgmax
               if(stable.and.nsharps.ge.msharps) exit
               prev_axes=t%axes
            end do
            if(lview%polish) then
               fix_rcoors=t%get_fix_rcoors('normal')
               call report('polish')
               hess0=hess(t%rcoors,fix_rcoors,'forward')
               call newton_powell(&
                      func=pot,grad=wrap_grad,&
                      v=t%rcoors,f=t%ener,&
                      conv_flag=conv_flag,&
                      fixed=fix_rcoors,&
                      gtols=gtol_rcoors,&
                      dvtols=dvtol_rcoors,&
                      mnp=mnp,&
                      hess0=hess0,&
                      stepcorr=dispcorr,&
                      rgmax=rgmax,&
                      dump=lview%dump,npad=lview%npad)
               call t%rcoors2carts()
            endif
            call t%reparametrize('clean')
            call t%calc_grad()
            t%rgmax=maxval(dabs(t%grad)/gtol_rcoors)
            write(2,'(i5,3x,l1)') i,conv_flag
         endif
         t%conv=conv_flag
      end associate
      !***
      contains
      !*****************************************************************
      !  function pot()
      !*****************************************************************
#        include "pot_geo_rigidMC_inc.f90"
        !end
      !*****************************************************************
      !  function grad()
      !*****************************************************************
#        include "grad_geo_rigidMC_inc.f90"
        !end
      !*****************************************************************
      !  function hess()
      !*****************************************************************
#        include "hess_geo_rigidMC_inc.f90"
        !end
      !*****************************************************************
      !  subroutine dispcorr()/disptry()
      !*****************************************************************
#        include "disp_geo_rigidMC_inc.f90"
        !end
      !*****************************************************************
         function wrap_grad(func,v,fixed) result(g)
      !*****************************************************************
         implicit none
         include 'iface/func_iface_inc.f90'
         real(r64) :: &
            v(:)
         logical, optional :: &
            fixed(:)
         real(r64), allocatable :: &
            g(:)
         !***
         g=grad(v,fixed)
         !***
         end
      !*****************************************************************
         subroutine report(mode)
      !*****************************************************************
         use mod_catch, only : catch_error
         implicit none
         character(*) mode
         !***
         if(.not.lview%dump) return
         select case(mode)
            case('header')
               clead=repeat(' ',lview%npad)
               if(lview%rad_relax) then
                  write(*,'(/a)') clead//&
                          'Performing a radial optimization ...'
               else
                  write(*,'(/a)') clead//&
                          'Performing a reconditioned optimization ...'
               endif
            case('recond')
               write(*,'(/3x,a,i0,a)') clead//'Reconditioning #',i,':'
            case default
               call catch_error(&
                        err=.true.,&
                        msg='unknown report mode "'//mode//'".',&
                        proc='geo%optimize')
            case('polish')
               write(*,'(/3x,a)') &
                  clead//'Polishing the BFGS results:'
         end select
         !***
         end
      !***
      end
!***********************************************************************
      subroutine register(this,last,is_unique,is_adequate,dump)
!***********************************************************************
      implicit none
      abstract interface
         subroutine sub()
         end
         logical function fun()
         end
      end interface
      procedure(fun), &
         optional, pointer :: &
            is_unique,&
            is_adequate
      procedure(sub), &
         optional, pointer :: &
            dump
      type(geo), pointer :: &
         this
      logical :: &
         last,&
         add
      integer :: &
         idup
      !***
      idup=0
      this%sort_dists=this%get_sort_dists()
      !$omp critical
      if(present(is_unique)) then
         add=is_unique()
      else
         add=is_unique_default()
      endif
      !$omp end critical
      if(add) then
         if(present(is_adequate)) then
            add=is_adequate()
         else
            add=is_adequate_default()
         endif
      endif
      if(add) then
         !$omp critical
         if(present(is_unique)) then
            add=is_unique()
         else
            add=is_unique_default()
         endif
         if(add) call append_registry()
         !$omp end critical
      endif
      if(.not.add) deallocate(this)
      if(last) then
         if(present(dump)) then
            call dump()
         else
            call dump_default()
         endif
      endif
      !***
      contains
      !*****************************************************************
         subroutine append_registry()
      !*****************************************************************
         implicit none
         type(geo), pointer :: &
            geo_ptr
         !***
         associate(registry=>landscape%registry)
            if(idup.gt.0) then
               geo_ptr=>registry%at(idup)
               call registry%erase(idup)
               deallocate(geo_ptr)
            endif
            call registry%push(this,-this%ener)
         end associate
         !***
         end
      !*****************************************************************
         function is_unique_default() result(unique)
      !*****************************************************************
         implicit none
         logical :: &
            unique
         type(geo), pointer :: &
            other
         integer :: &
            j
         !***
         unique=.true.
         idup=0
         associate(same_ener_rtol => sentinel%same_ener_rtol,&
                   grmsdLB_tol    => sentinel%grmsdLB_tol,&
                   registry       => landscape%registry)
            do j=1,registry%length
               other=>registry%at(j)
               if(dabs(this%ener-other%ener).ge.&
                  same_ener_rtol*dabs(this%ener)) cycle
               if(this%get_grmsdLB(other).ge.grmsdLB_tol) cycle
               unique=this%rgmax.lt.other%rgmax
               if(unique) idup=j
               return
            end do
         end associate
         !***
         end
      !*****************************************************************
         function is_adequate_default() result(adequate)
      !*****************************************************************
         implicit none
         logical :: &
            adequate
         !***
         associate(t         => this,&
                   max_order => control%max_order,&
                   skip_flat => control%skip_flat)
            adequate=.false.
            if(t%get_abnorm_flag().ne.0) return
            if(.not.t%conv) return
            if(max_order.eq.0) then
               call t%optimize(polish=.true.)
               if(t%get_abnorm_flag().ne.0) return
               call t%calc_grad()
               t%sort_dists=t%get_sort_dists()
               if(.not.t%conv) return
            endif
            call t%calc_hess_eigmodes()
            call t%calc_sens_hess_eigvals()
            if(skip_flat.and.t%get_flat_flag()) return
            if(t%get_saddle_order().gt.max_order) return
            adequate=.true.
         end associate
         !***
         end
      !*****************************************************************
         subroutine dump_default()
      !*****************************************************************
         use mod_string, only : to_string
         use mod_euclid, only : write_xyz
         implicit none
         type(geo_prior_queue) :: &
            registry0
         type(geo), pointer :: &
            geo_ptr
         character(:), allocatable :: &
            title,&
            pgname,&
            ftag
         integer :: &
            i,&
            order,&
            length
         !***
         associate(registry        => landscape%registry,&
                   stac_point_file => control%stac_point_file,&
                   symbs           => refframe%symbs)
            length=registry%length
            write(*,*)
            write(*,'(9x,a,i0)') &
               'Registered stationary points (SP): ',length
            write(*,*)
            registry0=registry
            do i=1,length
               geo_ptr=>registry0%pop()
               write(*,'(12x,a,i0,a)') 'SP #',i,':'
               call geo_ptr%print_summary(npad=15)
            end do
            write(*,'(9x,a)') &
                    'Printing the Cartesians into "'//&
                    trim(stac_point_file)//'" ...'
            registry0=registry
            do i=1,length
               geo_ptr=>registry0%pop()
               pgname=geo_ptr%get_pgroup_name()
               ftag=''
               if(geo_ptr%get_flat_flag()) ftag='FLAT'
               order=geo_ptr%get_saddle_order()
               title=to_string(geo_ptr%ener,3)//'  '//&
                     pgname//'  order = '//&
                     to_string(order)//' '//ftag
               call write_xyz(stac_point_file,title,&
                              symbs,geo_ptr%carts,append=i.gt.1)
            end do
         end associate
         !***
         end
      !***
      end
!***********************************************************************
      subroutine generate(this,rad_scan,delta)
!***********************************************************************
      implicit none
      type lview_
         logical :: &
            rad_scan=.true.
         real(r64) :: &
            delta=0.1d0
      end type
      type(lview_) :: &
         lview
      class(geo) :: &
         this
      logical, optional :: &
         rad_scan 
      real(r64), optional :: &
         delta
      !***
      if(present(rad_scan)) &
         lview%rad_scan=rad_scan
      if(present(delta)) &
         lview%delta=delta
      do
         call randomize()
         if(lview%rad_scan) call radial_scan()
         if(this%get_abnorm_flag().eq.0) exit
      end do
      !***
      contains
      !*****************************************************************
         subroutine randomize()
      !*****************************************************************
         implicit none
         type(monbox) :: &
            mb
         real(r64) :: &
            rnd
         integer :: &
            m
         !***
         associate(nmons      => layout%nmons,&
                   rho_ubound => landscape%rho_ubound)
            do m=1,nmons
               call mb%init(m,this)
               associate(alpha  => mb%alpha,&
                         beta   => mb%beta,&
                         gamma  => mb%gamma,&
                         phi    => mb%phi,&
                         theta  => mb%theta,&
                         rho    => mb%rho,&
                         calpha => mb%calpha,&
                         cbeta  => mb%cbeta,&
                         cgamma => mb%cgamma,&
                         cphi   => mb%cphi,&
                         ctheta => mb%ctheta,&
                         crho   => mb%crho)
                  !alpha
                  call random_number(rnd)
                  alpha=(calpha%up-calpha%low)*rnd+calpha%low
                  !beta
                  call random_number(rnd)
                  beta=(cbeta%up-cbeta%low)*rnd+cbeta%low
                  !gamma
                  call random_number(rnd)
                  gamma=(cgamma%up-cgamma%low)*rnd+cgamma%low
                  !phi
                  call random_number(rnd)
                  phi=(cphi%up-cphi%low)*rnd+cphi%low
                  !theta
                  call random_number(rnd)
                  theta=(ctheta%up-ctheta%low)*rnd+ctheta%low
                  !rho
                  if(lview%rad_scan) then
                     rho=(nmons-m)*rho_ubound
                  else
                     call random_number(rnd)
                     rho=(crho%up-crho%low)*rnd+crho%low
                  endif
               end associate
            end do
         end associate
         call this%rcoors2carts()
         !***
         end
      !*****************************************************************
         subroutine radial_scan()
      !*****************************************************************
         use mod_rand_perm, only : init_perm, rand_perm
         implicit none
         type(monbox) :: &
            mb
         integer, allocatable :: &
            perm(:)
         real(r64) :: &
            min_ener,&
            rho_min
         integer :: &
            m
         !***
         associate(t          => this,&
                   n          => layout%nmons-1,&
                   rho_ubound => landscape%rho_ubound)
            perm=init_perm(n)
            call rand_perm(perm)
            do m=1,n     
               call mb%init(m=perm(m),targ_geo=this)
               associate(rho=>mb%rho)
                  min_ener=huge(0.d0)
                  rho=0.d0
                  do
                     rho=rho+lview%delta
                     if(rho.gt.rho_ubound) exit
                     call t%rcoors2carts()
                     if(t%get_abnorm_flag().eq.-1) cycle
                     call t%calc_ener()
                     if(min_ener.gt.t%ener) then
                        min_ener=t%ener
                        rho_min=rho
                     endif
                  end do
                  rho=rho_min
                  call t%rcoors2carts()
               end associate
            end do
         end associate
         !***
         end
      !***
      end
