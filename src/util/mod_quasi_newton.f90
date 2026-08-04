!***********************************************************************
      module mod_quasi_newton
!***********************************************************************
      !  DESCRIPTION:
      !*****************************************************************
      !
      !  Local optimization framework providing Newton-Powell and BFGS
      !  quasi-Newton methods for stationary-point searches.
      !
      !*****************************************************************
      use iso_fortran_env, only : r64=>real64
      implicit none
      public :: bfgs, newton_powell
      private
      contains
      !*****************************************************************
         subroutine newton_powell(func,grad,v,f,conv_flag,fixed,gtols,&
                                  dvtols,mnp,dump,npad,hess0,niter,&
                                  rgmax,rdvmax,stepcorr,updtol,&
                                  mbtrust,rtrust,mintrust,maxtrust)
      !*****************************************************************
      !  DESCRIPTION:
      !*****************************************************************
      !
      !  This subroutine performs Newton-Powell optimization using
      !  exact and/or Powell-updated (diagonalized) Hessians.
      !  The step is the plain Newton step, regularized by updtol.
      !  Its correction is left to the optional stepcorr routine,
      !  which also receives the Hessian eigensystem.
      !
      !*****************************************************************
         use mod_catch, only : catch_error
         implicit none
         include "iface/func_iface_inc.f90"
         include "iface/grad_iface_inc.f90"
         include "iface/stepcorr_iface_inc.f90"
         include "iface/lapack_iface_inc.f90"
         optional :: &
            stepcorr
         !local view of optionals
         type lview_
            integer :: &
               mnp=100,&
               npad=0,&
               mbtrust=40
            real(r64) :: &
               rgmax=0.d0,&
               rdvmax=0.d0,&
               updtol=1.d-10,&
               rtrust=0.75d0,&
               mintrust=0.1d0,&
               maxtrust=0.5d0
            logical, allocatable :: &
               fixed(:)
            logical :: &
               dump=.false.
         end type
         type(lview_) :: &
            lview
         real(r64) :: &
            v(:),&
            f,&
            gtols(:)
         logical :: &
            conv_flag,&
            dump
         real(r64), optional :: &
            dvtols(:),&
            hess0(:,:),&
            rgmax,&
            rdvmax,&
            updtol,&
            rtrust,&
            mintrust,&
            maxtrust
         logical, optional :: &
            fixed(:)
         integer, optional :: &
            mnp,&
            npad,&
            niter,&
            mbtrust
         real(r64), allocatable :: &
            eigvals(:),&
            eigvects(:,:),&
            dv(:),&
            g(:),&
            dg(:)
         character(:), allocatable :: &
            clead
         real(r64) :: &
            emax
         integer :: &
            i,j
         !***
         associate(nv=>size(v))
            !set options
            if(present(mnp)) lview%mnp=mnp
            if(present(npad)) lview%npad=npad
            if(present(updtol)) lview%updtol=updtol
            if(present(mbtrust)) lview%mbtrust=mbtrust
            if(present(rtrust)) lview%rtrust=rtrust
            if(present(mintrust)) lview%mintrust=mintrust
            if(present(maxtrust)) lview%maxtrust=maxtrust

            !set fixed coordinates
            if(present(fixed)) then
               lview%fixed=fixed
            else
               lview%fixed=spread(.false.,1,nv)
            endif

            !check convergence tolerances
            call catch_error(&
                     err=nv.ne.size(gtols),&
                     msg='Tolerance size mismatch.',&
                     proc='newton_powell')
            if(present(dvtols)) &
               call catch_error(&
                        err=nv.ne.size(dvtols),&
                        msg='Tolerance size mismatch.',&
                        proc='newton_powell')

            !initialize gradient
            g=spread(0.d0,1,nv)
            lview%rdvmax=0.d0

            !output header
            clead=repeat(' ',lview%npad)
            if(dump) then
               if(present(dvtols)) then
                  write(*,'(/a,a11,3a13)') &
                     clead,'#','objfunc','rgmax','rdvmax'
               else
                  write(*,'(/a,a11,2a13)') &
                     clead,'#','objfunc','rgmax'
               endif
            endif

            !Newton-Powell iteration
            do i=1,lview%mnp
               if(present(niter)) niter=i
               !function value
               f=func(v)

               !gradient update
               dg=g
               g=grad(func,v,lview%fixed)
               dg=g-dg
               where(lview%fixed) g=0.d0
               lview%rgmax=maxval(dabs(g)/gtols)
               if(present(rgmax)) rgmax=lview%rgmax

               !if no displacement tolerance is present,
               !a converged gradient needs no Hessian
               if(.not.present(dvtols)) then
                  conv_flag=lview%rgmax.lt.1.d0
                  if(conv_flag) return
               endif

               !Hessian update
               call update_diag_hess()
               emax=maxval(dabs(eigvals))

               !project the gradient
               dv=matmul(transpose(eigvects),g)

               do j=1,nv
                  associate(dvj => dv(j),&
                            ej  => eigvals(j))
                     if(dabs(ej).lt.lview%updtol*emax) then
                        !no step for flat directions
                        dvj=0.d0
                     else
                        !Newton step
                        dvj=dvj/ej
                     endif
                  end associate
               end do

               !Newton step (back in the input coordinates)
               dv=-matmul(eigvects,dv)

               !respect fixed coordinates
               where(lview%fixed) dv=0.d0

               !external correction of the step
               if(present(stepcorr)) &
                  call stepcorr(v,dv,eigvals,eigvects)

               !convergence measure of the step
               lview%rdvmax=0.d0
               if(present(dvtols)) &
                  lview%rdvmax=maxval(dabs(dv)/dvtols)
               if(present(rdvmax)) rdvmax=lview%rdvmax

               !print diagnostics
               if(dump) then
                  if(present(dvtols)) then
                     write(*,'(a,i11,f13.2,2es13.3)') &
                           clead,i,f,lview%rgmax,lview%rdvmax
                  else
                     write(*,'(a,i11,f13.2,es13.3)') &
                           clead,i,f,lview%rgmax
                  endif
               endif

               !convergence test
               conv_flag=max(lview%rgmax,lview%rdvmax).lt.1.d0
               if(conv_flag) return

               !no step is possible any more
               if(all(dv.eq.0.d0)) return

               !scale step to ensure trust radius
               call trust_scale()

               !update variables
               v=v+dv
            end do
         end associate
         !*** 
         contains 
         !**************************************************************
            subroutine update_diag_hess()
         !**************************************************************
         !  DESCRIPTION:
         !**************************************************************
         !
         !  This subroutine performs the Powell Symmetric Broyden (PSB) 
         !  update of the Hessian matrix and returns its eigenvalues
         !  and eigenvectors.
         !
         !**************************************************************
            implicit none
            real(r64), allocatable :: &
               hess(:,:),&
               u(:),&
               work(:)
            real(r64) :: &
               dvn2,&
               udu
            integer :: &
               i,info
            !***
            associate(nv=>size(v))
               if(allocated(eigvals)) then
                  dvn2=dot_product(dv,dv)
                  if(dvn2.lt.lview%updtol) return
                  hess=matmul(eigvects*spread(eigvals,1,nv),&
                              transpose(eigvects))
                  u=dg-matmul(hess,dv)
                  udu=dot_product(u,dv)
                  hess=hess+&
                     (spread(u,1,nv)*spread(dv,2,nv)+&
                     spread(dv,1,nv)*spread(u,2,nv))/dvn2-&
                     udu*spread(dv,1,nv)*spread(dv,2,nv)/dvn2**2
               elseif(present(hess0)) then
                  call catch_error(&
                           err=any(shape(hess0).ne.[nv,nv]),&
                           msg='hess0 has incorrect shape.',&
                           proc='newton_powell')
                  hess=hess0
               else
                  hess=spread(spread(0.d0,1,nv),2,nv)
                  forall(i=1:nv) hess(i,i)=1.d0
               endif
               eigvals=spread(0.d0,1,nv)
               eigvects=hess
               work=spread(0.d0,1,3*nv-1)
               call dsyev('V','U',nv,eigvects,nv,eigvals,work,&
                          size(work),info)
               call catch_error(&
                        err=info.ne.0,&
                        msg='dsyev crashed.',&
                        proc='newton_powell')
            end associate
       
            end

         !**************************************************************
            subroutine trust_scale()
         !**************************************************************
         !  DESCRIPTION:
         !**************************************************************
         !
         !  This subroutine shortens the step until the actual change of
         !  the objective function agrees with the one predicted by the
         !  local quadratic model.
         !
         !**************************************************************
            implicit none
            real(r64), allocatable :: &
               w(:)
            real(r64) :: &
               t1,t2,s,&
               df,dfpred,&
               curv,sbase
            integer :: &
               k
            !***
            t1=dot_product(g,dv)
            w=matmul(transpose(eigvects),dv)
            t2=dot_product(w,eigvals*w)
            s=1.d0
            do k=1,lview%mbtrust
               df=func(v+s*dv)-f
               dfpred=s*t1+0.5d0*s**2*t2
               if(df*dfpred.ge.lview%rtrust*dfpred**2) exit
               curv=2.d0*(df-s*t1)/s**2
               sbase=0.d0
               if(curv.ne.0.d0) sbase=-t1/curv
               s=min(max(sbase,lview%mintrust*s),lview%maxtrust*s)
            end do
            dv=s*dv
            !***
            end
         !***
         end

      !*****************************************************************
         subroutine bfgs(func,grad,v,f,fixed,mbfgs,gtols,dvtols,hess0,&
                         dump,npad,scals,rgmax,rdvmax,conv_flag,&
                         updtol)
      !*****************************************************************
      !  DESCRIPTION:
      !*****************************************************************
      !
      !  This subroutine performs quasi-Newton minimization of a
      !  multivariate function using the Broyden-Fletcher-Goldfarb-
      !  Shanno (BFGS) update scheme, applied to the diagonalized
      !  Hessian.
      !  A coordinate is converged when its gradient falls below gtols
      !  and the displacement that gradient implies falls below dvtols.
      !
      !*****************************************************************
         use mod_catch, only : catch_error
         implicit none
         include "iface/func_iface_inc.f90"
         include "iface/grad_iface_inc.f90"
         include "iface/lapack_iface_inc.f90"
         type lview_
            integer :: &
               mbfgs=100,&
               npad=0
            logical :: &
               conv_flag=.false.
            real(r64) :: &
               rgmax=0.d0,&
               rdvmax=0.d0,&
               updtol=1.d-10
            real(r64), allocatable :: &
               scals(:)
            logical, allocatable :: &
               fixed(:)
         end type
         type(lview_) :: &
            lview
         real(r64) :: &
            v(:),&
            f,&
            gtols(:),&
            fdd,&
            emax,&
            resdv
         logical :: &
            dump
         real(r64), optional :: &
            dvtols(:),&
            scals(:),&
            hess0(:,:),&
            rgmax,&
            rdvmax,&
            updtol
         logical, optional :: &
            fixed(:),&
            conv_flag
         integer, optional :: &
            mbfgs,&
            npad
         real(r64), allocatable :: &
            eigvals(:),&
            eigvects(:,:),&
            dv(:),&
            g(:),&
            dg(:)
         character(:), allocatable :: &
            clead
         integer :: &
            i,j
         !***
         associate(nv=>size(v))
            !read options
            if(present(mbfgs)) lview%mbfgs=mbfgs
            if(present(npad)) lview%npad=npad
            if(present(updtol)) lview%updtol=updtol

            !set scales
            if(present(scals)) then
               lview%scals=scals
            else
               lview%scals=spread(1.d0,1,nv)
            endif

            !set fixed coordinates
            if(present(fixed)) then
               lview%fixed=fixed
            else
               lview%fixed=spread(.false.,1,nv)
            endif


            !check convergence tolerances
            call catch_error(&
                     err=nv.ne.size(gtols),&
                     msg='Tolerance size mismatch.',&
                     proc='bfgs')
            if(present(dvtols)) &
               call catch_error(&
                        err=nv.ne.size(dvtols),&
                        msg='Tolerance size mismatch.',&
                        proc='bfgs')

            !check scale factors
            call catch_error(&
                     err=any(dabs(lview%scals).lt.epsilon(0.d0)),&
                     msg='Zero scale factor passed.',&
                     proc='bfgs')

            !initialize the Hessian eigensystem
            call update_diag_hess()

            !initial state
            f=func(v)
            g=grad(func,v,lview%fixed)

            !output header
            clead=repeat(' ',lview%npad)
            if(dump) then
               if(present(dvtols)) then
                  write(*,'(/a,a11,3a13)') &
                     clead,'#','objfunc','rgmax','rdvmax'
               else
                  write(*,'(/a,a11,2a13)') &
                     clead,'#','objfunc','rgmax'
               endif
            endif

            !BFGS iteration
            do i=1,lview%mbfgs
               emax=maxval(dabs(eigvals))

               !project gradient
               dv=matmul(transpose(eigvects),g)
               do j=1,nv
                  associate(dvj => dv(j),&
                            ej  => eigvals(j))
                     if(ej.lt.lview%updtol*emax) then
                        !no step for flat directions
                        dvj=0.d0
                     else
                        !Newton step
                        dvj=dvj/ej
                     endif
                  end associate
               end do
               dv=-matmul(eigvects,dv)

               !respect fixed coordinates
               where(lview%fixed) dv=0.d0

               !step rescaling
               resdv=1.d0/max(maxval(dabs(dv/lview%scals)),1.d0)
               dv=resdv*dv

               !directional derivative
               fdd=dot_product(g,dv)

               !line search
               call blsrch(func,v,dv,f,fdd)

               !gradient update
               dg=g
               g=grad(func,v,lview%fixed)
               dg=g-dg
               where(lview%fixed) g=0.d0

               !convergence measures
               lview%rgmax=maxval(dabs(g)/gtols)
               lview%rdvmax=0.d0
               if(present(dvtols)) &
                  lview%rdvmax=maxval(dabs(dv)/dvtols)
               if(present(rgmax)) rgmax=lview%rgmax
               if(present(rdvmax)) rdvmax=lview%rdvmax

               !print diagnostics
               if(dump) then
                  if(present(dvtols)) then
                     write(*,'(a,i11,f13.2,2es13.3)') &
                        clead,i,f,lview%rgmax,lview%rdvmax
                  else
                     write(*,'(a,i11,f13.2,es13.3)') &
                        clead,i,f,lview%rgmax
                  endif
               endif

               !convergence test
               lview%conv_flag=&
                  max(lview%rgmax,lview%rdvmax).lt.1.d0
               if(present(conv_flag)) conv_flag=lview%conv_flag
               if(lview%conv_flag) return

               !no step is possible any more
               if(all(dv.eq.0.d0)) return

               !BFGS update of the diagonalized Hessian
               call update_diag_hess()
            end do

            !final status
            if(present(conv_flag)) conv_flag=.false.
         end associate
         !***
         contains
         !**************************************************************
            subroutine update_diag_hess()
         !**************************************************************
         !  DESCRIPTION:
         !**************************************************************
         !
         !  This subroutine performs the BFGS update of the Hessian
         !  matrix and returns its eigenvalues and eigenvectors.
         !
         !**************************************************************
            implicit none
            real(r64), allocatable :: &
               hess(:,:),&
               hdv(:),&
               work(:)
            real(r64) :: &
               sigma,&
               dhd
            integer :: &
               i,info
            !***
            associate(nv  => size(v),&
                      tol => lview%updtol)
               if(allocated(eigvals)) then
                  hess=matmul(eigvects*spread(eigvals,1,nv),&
                              transpose(eigvects))
                  hdv=matmul(hess,dv)
                  dhd=dot_product(dv,hdv)
                  if(dhd.le.tol*norm2(dv)*norm2(hdv)) return
                  sigma=dot_product(dg,dv)
                  if(sigma.le.tol*norm2(dg)*norm2(dv)) return
                  hess=hess+&
                     spread(dg,1,nv)*spread(dg,2,nv)/sigma-&
                     spread(hdv,1,nv)*spread(hdv,2,nv)/dhd
               elseif(present(hess0)) then
                  call catch_error(&
                           err=any(shape(hess0).ne.[nv,nv]),&
                           msg='hess0 has incorrect shape.',&
                           proc='bfgs')
                  hess=hess0
               else
                  hess=spread(spread(0.d0,1,nv),2,nv)
                  forall(i=1:nv) hess(i,i)=1.d0
               endif
               eigvals=spread(0.d0,1,nv)
               eigvects=hess
               work=spread(0.d0,1,3*nv-1)
               call dsyev('V','U',nv,eigvects,nv,eigvals,work,&
                          size(work),info)
               eigvals=dabs(eigvals)
               call catch_error(&
                        err=info.ne.0,&
                        msg='dsyev crashed.',&
                        proc='bfgs')
            end associate
            !***
            end
         !***
         end
      !*****************************************************************
         subroutine blsrch(func,v,dv,f,fdd,mblsrch,agtol,minlam,&
                           minfac,maxfac)
      !*****************************************************************
      !  DESCRIPTION:
      !*****************************************************************
      !
      !  Backtracking line search using quadratic/cubic interpolation
      !  (see also Numerical Recipes in Fortran 77, pp. 377-378).
      !
      !*****************************************************************
         implicit none
         include "iface/func_iface_inc.f90"
         type lview_
            integer :: &
               mblsrch=100
            real(r64) :: &
               agtol=1.d-4,&
               minlam=1.d-6,&
               minfac=0.1d0,&
               maxfac=0.5d0
         end type
         type(lview_) :: &
            lview
         real(r64), optional :: &
            agtol,&
            minlam,&
            minfac,&
            maxfac
         integer, optional :: &
            mblsrch
         real(r64) :: &
            v(:),&
            dv(:),&
            f,&
            fdd,&
            lam,&
            lam_prev,&
            f0,&
            dlam,&
            q,&
            q_prev,&
            a,&
            b,&
            disc
         real(r64), allocatable :: &
            v0(:)
         integer :: &
            i
         !*** 
         !set options
         if(present(mblsrch)) lview%mblsrch=mblsrch
         if(present(agtol)) lview%agtol=agtol
         if(present(minlam)) lview%minlam=minlam
         if(present(minfac)) lview%minfac=minfac
         if(present(maxfac)) lview%maxfac=maxfac
         block
            !initial state
            v0=v; f0=f; lam=1.d0

            !backtracking loop
            do i=1,lview%mblsrch
               v=v0
               !no step is possible any more
               if(lam.lt.lview%minlam) then
                  f=f0
                  exit
               endif
               v=v+lam*dv
               f=func(v)
               if(f.le.f0+lview%agtol*lam*fdd) exit
               q=(f-f0-lam*fdd)/lam**2
               if(lam.eq.1.d0) then
                  a=0.d0
                  b=q
               else
                  dlam=lam-lam_prev
                  a=(q-q_prev)/dlam
                  b=(lam*q_prev-lam_prev*q)/dlam
               endif
               lam_prev=lam
               q_prev=q
               if(a.eq.0.d0) then
                  lam=-fdd/(2.d0*b)
               else
                  disc=b**2-3.d0*a*fdd
                  lam=(-b+dsqrt(disc))/(3.d0*a)
               endif
               lam=max(min(lam,lview%maxfac*lam_prev),&
                       lview%minfac*lam_prev)
            end do

            !final step
            dv=v-v0
         end block
         !*** 
         end
      !***
      end
