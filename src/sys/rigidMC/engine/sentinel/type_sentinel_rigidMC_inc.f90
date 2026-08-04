!**********************************************************************
      type sentinel_
!**********************************************************************
      integer :: &
         mreconds=200,&
         mdcbis=40,&
         msharps=5,&
         mbfgs_relax=10,&
         mrest_relax=3,&
         mbfgs=10,&
         mnp=10
      real(r64) :: &
         ang_f2doff=1.d-3,&
         rad_f2doff=1.d-5,&
         rad_gtol=0.05d0,&
         ang_gtol=1.d-3,&
         rad_dvtol=1.d-3,&
         ang_dvtol=0.01d0,&
         ang_scal=30.d0,&
         rad_scal=0.5d0,&
         sing_tol=1.d-12,&
         rig_tol=1.d-5,&
         flat_fac=2.d0,&
         collin_tol=1.d-3,&
         pgrp_tol=5.d-2,&
         rsmdet_tol=0.9d0,&
         same_ener_rtol=1.d-3,&
         grmsdLB_tol=0.025d0,&
         lam_cov=1.5d0,&
         lam_vdw=1.5d0,&
         recond_tol=0.05d0
      real(r64), allocatable :: &
         f2doffsets(:)
      contains
         procedure :: &
            init=>init_sentinel
      end type
      !end
