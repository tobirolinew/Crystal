!***********************************************************************
      subroutine init_monbox(this,m,targ_geo)
!***********************************************************************
      use mod_catch, only : catch_error
      implicit none
      class(monbox) :: &
         this
      character(:), pointer :: &
         tmp(:)
      type(geo), target, optional :: &
         targ_geo
      integer :: &
         m
      !***
      associate(t   => this,&
                con => control,&
                lay => layout,&
                ref => refframe,&
                lan => landscape)
         !check
         call catch_error(&
                  err=m.lt.1.or.m.gt.lay%nmons,&
                  msg='Index "m" is out of range',&
                  proc='monbox%init')
         !indices
         t%m=m
         t%ibase=6*(m-1)
         !scalar parameters
         t%ref_xyz_file => con%ref_xyz_files(m)
         t%first        => lay%abounds(m,1)
         t%last         => lay%abounds(m,2)
         t%natoms       => lay%mon_natoms(m)
         t%bin          => lay%mon_bins(m)
         t%nedof        => ref%nedofs(m)
         t%pgroup       => ref%pgroups(m)
         !atomic data   
         t%symbs        => ref%symbs(t%first:t%last)
         t%anums        => ref%anums(t%first:t%last)
         t%ref_carts    => ref%carts(:,t%first:t%last)
         !parameters related to targ_geo (optional)
         if(present(targ_geo)) then
            !conditional allocations
            !coorboxes
            associate(i=>t%ibase)
               call t%calpha%init(i+1,targ_geo)
               call t%cbeta%init(i+2,targ_geo)
               call t%cgamma%init(i+3,targ_geo)
               call t%cphi%init(i+4,targ_geo)         
               call t%ctheta%init(i+5,targ_geo)
               call t%crho%init(i+6,targ_geo)
            end associate
            !targ_geo
            t%targ_geo => targ_geo
            t%carts    => targ_geo%carts(:,t%first:t%last)
            !rigid coordinates
            t%alpha    => t%calpha%val
            t%beta     => t%cbeta%val
            t%gamma    => t%cgamma%val
            t%phi      => t%cphi%val
            t%theta    => t%ctheta%val
            t%rho      => t%crho%val
            !direct rigid data arrays
            associate(rcfirst => t%ibase+1,&
                      rclast  => t%ibase+6)
               t%rcoors       => targ_geo%rcoors(rcfirst:rclast)
               t%low_rcoors   => lan%low_rcoors(rcfirst:rclast)
               t%up_rcoors    => lan%up_rcoors(rcfirst:rclast)
               t%step_rcoors  => lan%step_rcoors(rcfirst:rclast)
               t%gtol_rcoors  => lan%step_rcoors(rcfirst:rclast)
               associate(lab_len  => len(lan%lab_rcoors))
                  allocate(character(lab_len)::t%lab_rcoors(0))
                  tmp             => t%lab_rcoors
                  t%lab_rcoors    => lan%lab_rcoors(rcfirst:rclast)
                  deallocate(tmp)
               end associate
               associate(unit_len => len(lan%unit_rcoors))
                  allocate(character(unit_len)::t%unit_rcoors(0))
                  tmp             => t%unit_rcoors
                  t%unit_rcoors   => lan%unit_rcoors(rcfirst:rclast)
                  deallocate(tmp)
               end associate
            end associate
         endif
      end associate
      !***
      end
