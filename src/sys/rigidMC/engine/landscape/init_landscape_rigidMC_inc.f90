!***********************************************************************
      subroutine init_landscape(this)
!***********************************************************************
      use mod_pertab, only : vdw_rads
      use mod_string, only : to_string
      implicit none
      class(landscape_) :: &
         this
      type(monbox) :: &
         mb
      character(:), allocatable :: &
         vmi
      integer :: &
         m,i
      !***
      associate(t         => this,&
                nmons     => layout%nmons,&
                ang_step  => control%ang_step,&
                rad_step  => control%rad_step,& 
                ang_unit  => control%ang_unit,&
                rad_unit  => control%rad_unit,&
                nedofs    => refframe%nedofs,&
                ang_gtol  => sentinel%ang_gtol,&
                rad_gtol  => sentinel%rad_gtol,&
                ang_dvtol => sentinel%ang_dvtol,&
                rad_dvtol => sentinel%rad_dvtol,&
                ang_scal  => sentinel%ang_scal,&
                rad_scal  => sentinel%rad_scal,&
                lam       => sentinel%lam_vdw,&
                xyz_file  => control%start_xyz_file)
         !nrcoors
         t%nrcoors=6*nmons
         !rho_ubound
         t%rho_ubound=0.d0
         do m=1,nmons
            call mb%init(m)
            t%rho_ubound=t%rho_ubound+&
                  !monomer extent:
                     2.d0*maxval(&
                        !distances from COM:
                           norm2(mb%ref_carts,dim=1)+&
                        !effective vdW radii
                           lam*vdw_rads(mb%anums))
         end do
         !low_rcoors 
         associate(v=>[-180.d0,&  !alpha
                          0.d0,&  !beta
                       -180.d0,&  !gamma
                       -180.d0,&  !phi
                          0.d0,&  !theta
                          0.d0])  !rho
            t%low_rcoors=[(v,m=1,nmons)]
         end associate
         !up_rcoors 
         associate(v=>[180.d0,&   !alpha
                       180.d0,&   !beta
                       180.d0,&   !gamma
                       180.d0,&   !phi
                       180.d0,&   !theta
                 t%rho_ubound])   !rho
            t%up_rcoors=[(v,m=1,nmons)]
         end associate
         !step_rcoors
          associate(v=>[ang_step,& !alpha
                        ang_step,& !beta
                        ang_step,& !gamma
                        ang_step,& !phi
                        ang_step,& !theta
                        rad_step]) !rho
            t%step_rcoors=[(v,m=1,nmons)]
         end associate
         !gtol_rcoors
         associate(v=>[ang_gtol,& !alpha
                       ang_gtol,& !beta
                       ang_gtol,& !gamma
                       ang_gtol,& !phi
                       ang_gtol,& !theta
                       rad_gtol]) !rho
            t%gtol_rcoors=[(v,m=1,nmons)]
         end associate
         !dvtol_rcoors
         associate(v=>[ang_dvtol,& !alpha
                       ang_dvtol,& !beta
                       ang_dvtol,& !gamma
                       ang_dvtol,& !phi
                       ang_dvtol,& !theta
                       rad_dvtol]) !rho
            t%dvtol_rcoors=[(v,m=1,nmons)]
         end associate
         !scal_rcoors
         associate(v=>[ang_scal,& !alpha
                       ang_scal,& !beta
                       ang_scal,& !gamma
                       ang_scal,& !phi
                       ang_scal,& !theta
                       rad_scal]) !rho
            t%scal_rcoors=[(v,m=1,nmons)]
         end associate
         !lab_rcoors
         t%lab_rcoors=[character(0)::]
         associate(v=>['alpha',&
                       'beta ',&
                       'gamma',&
                       'phi  ',&
                       'theta',&
                       'rho  '])
            do m=1,nmons
               do i=1,6
                  vmi=trim(v(i))//to_string(m)
                  call append_string(arr=t%lab_rcoors,elem=vmi)
               end do
            end do
         end associate
         !unit_rcoors
         t%unit_rcoors=[character(0)::]
         do m=1,nmons
            !angular:
            do i=1,5
               call append_string(arr=t%unit_rcoors,elem=ang_unit)
            end do
            !radial:
            call append_string(arr=t%unit_rcoors,elem=rad_unit)
         end do
         !start_geo
         allocate(t%start_geo)
         call t%start_geo%init(&
                  xyz_file,opt=1,hess=.true.,&
                  dump=.true.,npad=12)
         call t%start_geo%reparametrize('assert_clean')
      end associate
      !***
      contains
      !*****************************************************************
      !  subroutine append_string()
      !*****************************************************************
#        define STR_ARR_TYPE character(:)
#        define APPEND append_string
#        include "../../../../util/dep/generic_append_inc.f90"
         !end
      end
