!***********************************************************************
      function hess(rcoors,fixed,mode) result(h)
!***********************************************************************
      use mod_catch, only : catch_error
      implicit none
      real(r64) :: &
         rcoors(:)
      logical :: &
         fixed(:)
      character(*), optional :: &
         mode
      character(:), allocatable :: &
         lmode
      real(r64), allocatable :: &
         h(:,:),&
         prcoors(:),&
         g0(:),&
         gp(:),&
         gm(:)
      integer :: &
         k
      !***
      associate(pot_drive  => control%pot_drive,&
                nrcoors    => landscape%nrcoors,&
                f2doffsets => sentinel%f2doffsets)
         lmode='central'
         if(present(mode)) lmode=mode
         h=spread(spread(0.d0,1,nrcoors),2,nrcoors)
         select case(pot_drive)
            case('f2b_homo')
               if(lmode.eq.'forward') g0=grad(rcoors,fixed)
               do k=1,nrcoors
                  if(fixed(k)) cycle
                  associate(off=>f2doffsets(k))
                     prcoors=rcoors
                     prcoors(k)=rcoors(k)+off
                     gp=grad(prcoors,fixed)
                     select case(lmode)
                        case('central')
                           prcoors(k)=rcoors(k)-off
                           gm=grad(prcoors,fixed)
                           h(:,k)=(gp-gm)/(2.d0*off)
                        case('forward')
                           h(:,k)=(gp-g0)/off
                        case default
                           call catch_error(&
                                   err=.true.,&
                                   msg='Unknown mode "'//lmode//'".',&
                                   proc='hess')
                     end select
                  end associate
               end do
            case default
               call catch_error(&
                       err=.true.,&
                       msg='Unknown pot_drive "'//'"'//&
                            pot_drive//'"',&
                       proc='hess')
         end select
         h=0.5d0*(h+transpose(h))
      end associate
      !***
      end
