!***********************************************************************
      subroutine init_sentinel(this)
!***********************************************************************
      implicit none
      class(sentinel_) :: &
         this
      integer &
         m
      !***
      associate(nmons      => layout%nmons,&
                ang_f2doff => sentinel%ang_f2doff,&
                rad_f2doff => sentinel%rad_f2doff)
         !f2d_offsets 
         associate(v=>[ang_f2doff,& !alpha
                       ang_f2doff,& !beta
                       ang_f2doff,& !gamma
                       ang_f2doff,& !phi
                       ang_f2doff,& !theta
                       rad_f2doff]) !rho
            this%f2doffsets=[(v,m=1,nmons)]
         end associate
      end associate
      !***
      end
