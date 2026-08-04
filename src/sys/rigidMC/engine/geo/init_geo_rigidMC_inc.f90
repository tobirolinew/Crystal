!***********************************************************************
      subroutine init_geo(this,xyz_file,opt,hess,dump,npad)
!***********************************************************************
      use mod_string, only : to_string
      use mod_catch, only : catch_error
      use mod_f2b
      implicit none
      type lview_
         integer :: opt=0
         logical :: hess=.false.
         logical :: dump=.false.
         integer :: npad=9
      end type
      type(lview_) :: &
         lview
      class(geo) :: &
         this
      character(*), optional :: &
         xyz_file
      integer, optional :: &
         opt
      logical, optional :: &
         hess,&
         dump
      integer, optional :: &
         npad
      character(:), allocatable :: &
         clead
      !***
      if(present(opt))  lview%opt=opt
      if(present(hess)) lview%hess=hess
      if(present(dump)) lview%dump=dump
      if(present(npad)) lview%npad=npad
      associate(t         => this,&
                max_order => control%max_order,&
                natoms    => layout%natoms,&
                nrcoors   => landscape%nrcoors)
         clead=repeat(' ',lview%npad)
         t%carts=spread([0.d0,0.d0,0.d0],2,natoms)
         t%rcoors=spread(0.d0,1,nrcoors)
         if(present(xyz_file)) then
            if(xyz_file.ne.'none') then
               call t%read_carts(xyz_file)
            else
               call t%generate()
            endif
            call t%reorient()
            call t%carts2rcoors()
            call t%freeze()
         endif
         call t%calc_ener()
         if(lview%dump) then
            write(*,'(/a)') clead//&
                    'Geometry before optimization:'
            call t%print_summary(npad=lview%npad+3)
         endif
         if(lview%hess) then 
            call t%calc_grad()
            call t%calc_hess_eigmodes()
            call t%calc_sens_hess_eigvals()
         endif
         if(lview%dump) & 
            call t%print_summary(npad=lview%npad+3)
         select case(lview%opt)
            case(0)
               !no optimization
            case(-1,1)
               call t%optimize(&
                        rad_relax=lview%opt.eq.-1,&
                        dump=lview%dump,&
                        npad=npad,&
                        polish=.true.)
            case default
               call catch_error(&
                        err=.true.,&
                        msg='incorrect "opt" value'//&
                             to_string(opt)//'".',&
                        proc='geo%init',&
                        comment='only -1, 0, or 1 is allowed.')
         end select
         if(lview%hess) then 
            call t%calc_grad()
            call t%calc_hess_eigmodes()
            call t%calc_sens_hess_eigvals()
         endif
         if(lview%dump) then
            write(*,'(/a)') clead//&
                    'Geometry after optimization:'
            call t%print_summary(npad=lview%npad+3)
         endif
      end associate
      !***
      end
