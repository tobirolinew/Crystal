!***********************************************************************
      subroutine init_refframe(this)
!***********************************************************************
      use mod_euclid, only : read_xyz
      use mod_point_group, only : pgroup, get_pgroup
      implicit none
      class(refframe_) :: &
         this
      type(monbox) :: &
         mb
      integer :: &
         m
      !***
      associate(t       => this,&
                files   => control%ref_xyz_files,&
                nmons   => layout%nmons,&
                natoms  => layout%natoms)
         t%symbs=spread('  ',1,natoms)
         t%anums=spread(0,1,natoms)
         t%carts=spread([0.d0,0.d0,0.d0],2,natoms)
         t%nedofs=spread(0,1,nmons)
         t%pgroups=[(pgroup(),m=1,nmons)]
         do m=1,nmons
            associate(pg    => t%pgroups(m),&
                      nedof => t%nedofs(m))
               call mb%init(m)
               !geometry
               call read_xyz(files(m),mb%anums,mb%symbs,mb%ref_carts)
               !point group
               pg=get_pgroup(mb%anums,mb%ref_carts)
               !number of external DOFs
               select case(pg%pgname)
                  case('Kh')
                     nedof=3
                  case('Dih','Civ')
                     nedof=5
                  case default
                     nedof=6
               end select
            end associate
         end do
         t%max_nedof=maxval(t%nedofs)
      end associate
      !***
      end
