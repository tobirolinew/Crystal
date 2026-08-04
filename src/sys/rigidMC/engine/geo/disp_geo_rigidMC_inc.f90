!***********************************************************************
      subroutine dispcorr(v,dv,eigvals,eigvects)
!***********************************************************************
      implicit none
      type(geo) :: &
         snapshot
      logical :: &
         bad,&
         skip
      real(r64) :: &
         v(:),&
         dv(:),&
         eigvals(:),&
         eigvects(:,:),&
         eta,&
         low_eta,&
         up_eta
      real(r64), allocatable :: &
         v0(:),&
         p(:)
      integer :: &
         k,twice
      !***
      associate(t      => this,&
                scals  => landscape%scal_rcoors,&
                mdcbis => sentinel%mdcbis)
         snapshot%rcoors=t%rcoors
         snapshot%carts=t%carts
         p=matmul(transpose(eigvects),dv)
         up_eta=norm2(eigvals*p)/minval(scals)
         v0=v
         eta=0.d0
         low_eta=0.d0
         skip=.false.
         do twice=1,2
            do k=1,mdcbis
               call disptry(v0,p,eigvals,eigvects,eta,dv,bad)
               if(bad) then
                  low_eta=eta
               else
                  up_eta=eta
               endif
               if(twice.eq.1) then
                  if(.not.bad) exit
                  eta=max(2.d0*eta,up_eta)
               else
                  eta=0.5d0*(low_eta+up_eta)
               endif
            end do
            if(twice.eq.1) then
               skip=bad
               if(skip.or.up_eta.le.0.d0) exit
               eta=0.5d0*(low_eta+up_eta)
            endif
         end do
         if(skip) then
            dv=0.d0
         elseif(up_eta.gt.0.d0) then
            call disptry(v0,p,eigvals,eigvects,up_eta,dv,bad)
         endif
         t%rcoors=snapshot%rcoors
         t%carts=snapshot%carts
      end associate
      !***
      end
!***********************************************************************
      subroutine disptry(v,p,eigvals,eigvects,eta,dv,bad)
!***********************************************************************
      implicit none
      real(r64) :: &
         v(:),&
         p(:),&
         eigvals(:),&
         eigvects(:,:),&
         eta,&
         dv(:),&
         dvmax
      logical :: &
         bad
      !***
      associate(t           => this,&
                scals       => landscape%scal_rcoors,&
                abs_eigvals => dabs(eigvals))
         dv=matmul(eigvects,p*abs_eigvals/(abs_eigvals+eta))
         where(fix_rcoors) dv=0.d0
         dvmax=maxval(dabs(dv)/scals)
         bad=dvmax.gt.1.d0
         if(bad) return
         t%rcoors=v+dv
         call t%rcoors2carts()
         bad=t%get_abnorm_flag().ne.0
      end associate
      !***
      end
