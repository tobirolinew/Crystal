!***********************************************************************
      subroutine globalize_rigidMC()
!***********************************************************************
      use mod_global, only : &
         ncoors,&
         fix_wcoors,&
         start_wcoors,&
         low_wcoors,&
         up_wcoors,&
         step_wcoors,&
         cost_label,&
         cost_unit,&
         lab_wcoors,&
         unit_wcoors
      implicit none
      !***
      write(*,'(3x,a/)') 'Calling globalize_rigidMC() ...'
      associate(con       => control,&
                lan       => landscape,&
                start_geo => landscape%start_geo)
          ncoors=lan%nrcoors
          cost_label=con%pot_label
          cost_unit=con%pot_unit
          start_wcoors=start_geo%rcoors
          fix_wcoors=start_geo%get_fix_rcoors('crystal')
          low_wcoors=lan%low_rcoors
          up_wcoors=lan%up_rcoors
          step_wcoors=lan%step_rcoors
          lab_wcoors=lan%lab_rcoors
          unit_wcoors=lan%unit_rcoors
      end associate
      !***
      end
