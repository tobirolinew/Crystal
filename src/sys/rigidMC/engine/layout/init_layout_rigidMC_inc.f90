!***********************************************************************
      subroutine init_layout(this)
!***********************************************************************
      use mod_catch, only : catch_error
      use mod_string, only : expand_tilde
      implicit none
      class(layout_) :: &
         this
      integer :: &
         m,&
         n,&
         first,&
         last,&
         ios
      !***
      associate(t=>this)
         !nmons
         t%nmons=size(control%ref_xyz_files)
         !mon_natoms
         t%mon_natoms=spread(0,1,t%nmons)
         associate(files=>control%ref_xyz_files)
            do m=1,t%nmons
               open(1,file=expand_tilde(files(m)),&
                    status='old',iostat=ios)
               read(1,*,iostat=ios) t%mon_natoms(m)
               close(1)
               call catch_error(&
                        err=ios.ne.0,&
                        msg='failed to read atom count from '//&
                            'file "'//files(m)//'."',&
                        proc='layout%init')
            end do
         end associate
         !mon_bins
         t%mon_bins=spread(0,1,t%nmons)
         associate(files=>control%ref_xyz_files)
            do m=1,t%nmons
               n=findloc(files(:m-1),files(m),dim=1)
               if(n.eq.0) then
                  t%mon_bins(m)=maxval(t%mon_bins(:m-1))+1
               else
                  t%mon_bins(m)=t%mon_bins(n)
               endif
            end do
         end associate
         !natoms
         t%natoms=sum(t%mon_natoms)
         !abounds
         t%abounds=spread(spread(0,1,t%nmons),2,2)
         last=0
         do m=1,t%nmons
            first=last+1
            last=first+t%mon_natoms(m)-1
            t%abounds(m,:)=[first,last]
         end do
         !mon_inds
         t%mon_inds=spread(0,1,t%natoms)
         do m=1,t%nmons
            first=t%abounds(m,1)
            last=t%abounds(m,2)
            t%mon_inds(first:last)=m
         end do
      end associate
      !***
      end
