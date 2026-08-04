!**********************************************************************
      type control_
!**********************************************************************
      integer :: &
         max_targets,&
         max_order
      real(r64) :: &
         ang_step,&
         rad_step,&
         max_ener
      character(:), allocatable :: &
         ref_xyz_files(:),&
         start_xyz_file,&
         pot_label,&
         pot_unit,&
         pot_drive,&
         pot_spec,&
         ang_unit,&
         rad_unit,&
         stac_point_file
      logical :: &
         rad_relax,&
         skip_flat
      contains
         procedure :: init=>init_control
      end type
      !end
