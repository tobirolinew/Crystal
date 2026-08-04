!***********************************************************************
      subroutine calc_ener(this)
!***********************************************************************
      implicit none
      class(geo) :: &
         this
      !***
      this%ener=pot(this%rcoors)
      !***
      contains
      !*****************************************************************
      !  function pot()
      !*****************************************************************
#        include "pot_geo_rigidMC_inc.f90"
        !end
      !***
      end
!***********************************************************************
      subroutine calc_grad(this)
!***********************************************************************
      implicit none
      class(geo) :: &
         this
      logical, allocatable :: &
         fix_rcoors(:)
      !***
      associate(t => this)
         call t%reparametrize('run')
         fix_rcoors=t%get_fix_rcoors('normal')
         t%grad=grad(t%rcoors,fix_rcoors)
         call t%reparametrize('clean')
      end associate
      !***
      contains
      !*****************************************************************
      !  function grad()
      !*****************************************************************
#        include "grad_geo_rigidMC_inc.f90"
        !end
      !***
      end
!***********************************************************************
      subroutine calc_hess_eigmodes(this)
!***********************************************************************
      use mod_sort, only: get_num_sort_perm
      use mod_catch, only : catch_error
      implicit none
      include 'iface/lapack_iface_inc.f90'
      logical, allocatable :: &
         fix_rcoors(:)
      class(geo) :: &
         this
      real(r64), allocatable :: &
         work(:)
      integer, allocatable :: &
         perm(:)
      integer :: &
         nfix,info
      !***
      associate(t=>this)
         call t%reparametrize('run')
         fix_rcoors=t%get_fix_rcoors('normal')
         t%hess_eigvects=hess(t%rcoors,fix_rcoors)
         associate(nv => size(t%rcoors))
            t%hess_eigvals=spread(0.d0,1,nv)
            work=spread(0.d0,1,3*nv-1)
            call dsyev('V','U',nv,t%hess_eigvects,nv,t%hess_eigvals,&
                       work,size(work),info)
         end associate
         call catch_error(&
                  err=info.ne.0,&
                  msg='DSYEV failed.',&
                  proc='geo%calc_hess_eigmodes')
         call t%reparametrize('clean')
         perm=get_num_sort_perm(dabs(t%hess_eigvals),asc=.true.)
         nfix=count(t%get_fix_rcoors('normal'))
         t%hess_eigvals(perm(:nfix))=0.d0
      end associate
      !***
      contains
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
      !***
      end
!***********************************************************************
      subroutine calc_sens_hess_eigvals(this)
!***********************************************************************
      implicit none
      class(geo) :: &
         this
      logical, allocatable :: &
         fix_rcoors(:)
      real(r64), allocatable :: &
         gpert(:),&
         dgrad(:),&
         res_dgrad(:),&
         proj_dgrad(:),&
         abs_devs(:,:)
      integer :: &
         i,ipert,n
      !***
      call catch_local_errors()
      associate(t          => this,&
                nrcoors    => landscape%nrcoors,&
                f2doffsets => sentinel%f2doffsets,&
                sing_tol   => sentinel%sing_tol)
         init: block
            call t%reparametrize('run')
            fix_rcoors=t%get_fix_rcoors('normal')
         end block init
         calc_abs_devs: block
          n=count(.not.fix_rcoors)
          allocate(abs_devs(n,nrcoors))
          ipert=0
          do i=1,nrcoors
             if(fix_rcoors(i)) cycle
             ipert=ipert+1
             associate(rcoor    => t%rcoors(i),&
                       eps      => f2doffsets(i),&
                       eigvals  => t%hess_eigvals,&
                       eigvects => t%hess_eigvects,&
                       row      => t%hess_eigvects(i,:))
               rcoor=rcoor+eps
               gpert=grad(t%rcoors,fix_rcoors)
               dgrad=gpert-t%grad
               rcoor=rcoor-eps
               res_dgrad=dgrad-eps*matmul(eigvects,eigvals*row)
               proj_dgrad=matmul(res_dgrad,eigvects)
               abs_devs(ipert,:)=&
                  dabs((2.d0*row*proj_dgrad-res_dgrad(i)*row**2)/eps)
             end associate
          end do
        end block calc_abs_devs         
         finalize: block
            t%sens_hess_eigvals=norm2(abs_devs,dim=1)/dsqrt(dble(n))
            call t%rcoors2carts()
            call t%reparametrize('clean')
         end block finalize
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
                     msg='array "hess_eigvals" is not allocated.',&
                     proc='geo%calc_sens_hess_eigvals')
            call catch_error(&
                     err=size(t%hess_eigvals).ne.nrcoors,&
                     msg='array "hess_eigvals" has incorrect size.',&
                     proc='geo%calc_sens_hess_eigvals')
            call catch_error(&
                     err=.not.allocated(t%hess_eigvects),&
                     msg='array "hess_eigvects" is not allocated.',&
                     proc='geo%calc_sens_hess_eigvals')
            call catch_error(&
                     err=any(nrcoors.ne.shape(t%hess_eigvects)),&
                     msg='array "hess_eigvects" has incorrect size.',&
                     proc='geo%calc_sens_hess_eigvals')
         end associate
         !***
         end
      !*****************************************************************
      !  function grad()
      !*****************************************************************
#        include "grad_geo_rigidMC_inc.f90"
        !end
      !***
      end
