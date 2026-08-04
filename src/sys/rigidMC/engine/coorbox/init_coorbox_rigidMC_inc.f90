!***********************************************************************
      subroutine init_coorbox(this,idx,targ_geo)
!***********************************************************************
      use mod_catch, only : catch_error
      implicit none
      class(coorbox) :: &
         this
      integer :: &
         idx
      type(geo), target :: &
         targ_geo
      !***
      associate(t   => this,&
                lan => landscape)
         !check
         call catch_error(&
                  err=idx.lt.1.or.idx.gt.lan%nrcoors,&
                  msg='Index "idx" is out of range',&
                  proc='coorbox%init')
         !index
         t%idx=idx
         !binds
         t%val    => targ_geo%rcoors(idx)
         t%low    => lan%low_rcoors(idx)
         t%up     => lan%up_rcoors(idx)
         t%step   => lan%step_rcoors(idx)
         t%gtol   => lan%gtol_rcoors(idx)
         t%lab    => lan%lab_rcoors(idx)
         t%unit   => lan%unit_rcoors(idx)
      end associate
      !***
      end
