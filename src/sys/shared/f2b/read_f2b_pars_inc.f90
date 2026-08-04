   !********************************************************************
      subroutine read_f2b_pars(f2b_file_name)
   !********************************************************************
      use mod_string, only : expand_tilde
      implicit none
      character(*) :: &
         f2b_file_name
      character(1000) :: &
        row
      real(r64) :: &
         sing_tol=1.d-6
      real(r64) :: &
         val
      complex(r64), allocatable :: &
         iso_cartsc(:,:)
      integer :: &
         length,&
         i,j,k,l,m
      !***
      !Open the F2B file
      open(1,file=expand_tilde(f2b_file_name))
      !Read the number of effective sites and monomers
      call read_colon_int(funit=1,val=nbof_sites,nskip=4)
      call read_colon_int(funit=1,val=nbof_mons,nskip=1)
      call skip_rows(funit=1,nskip=2)
      !Save the number of monomer atoms
      allocate(nbof_mon_atoms(nbof_mons))
      do m=1,nbof_mons
         call read_colon_int(funit=1,val=nbof_mon_atoms(m),nskip=1)
      end do
      !Save the off-atomic sites
      read(1,'(/a)') row
      if(index(row,'Number of off-atomic sites').ne.0) then
         allocate(nbof_mon_OA_sites(nbof_mons))
         do m=1,nbof_mons
            call read_colon_int(funit=1,&
                     val=nbof_mon_OA_sites(m),nskip=1)
         end do
         call skip_rows(funit=1,nskip=2)
      endif
      call skip_rows(funit=1,nskip=2)
      !Read the Cartesian coordinates of the
      !reference structure
      allocate(symbs(nbof_sites))
      allocate(iso_carts(3,nbof_sites))
      do i=1,nbof_sites
         read(1,'(a)') row
         read(row,*) symbs(i),iso_carts(:,i)
      end do
      !Save the number of atoms and off-atomic sites
      nbof_OA_sites=count(symbs(:)(1:1).eq.'X')
      nbof_atoms=nbof_sites-nbof_OA_sites
      read(1,'(/a)') row
      !Read the connectivity list of off-atomic sites
      if(nbof_OA_sites.ne.0) then
         if(index(row,"Off-atomic connectivity list").eq.0) then
            write(*,'(/2x,a/)') &
                    'Error: missing off-atomic connectivity list.'
            stop
         endif
         block
            allocate(lst_OA(3,nbof_OA_sites))
            allocate(w_OA(3,nbof_OA_sites))
            call skip_rows(funit=1,nskip=2)
            do i=1,nbof_OA_sites
               read(1,*) lst_OA(:,i)
            end do
            !Calculate the internal coordinates of the
            !off-atomic sites
            iso_cartsc=iso_carts
            call add_OA_sites(iso_cartsc,sing_tol,init=.true.)
            call skip_rows(funit=1,nskip=2)
         end block
      endif
      !Read the polynomial parameters of term "U_intra"
      call skip_rows(funit=1,nskip=3)
      call init_poly(Lambda_intra,S_intra,B_intra,&
                     read_length=.true.,degree=2)
      call skip_rows(funit=1,nskip=5)
      !Read the parameters of the charge polynomials
      allocate(sites(nbof_sites))
      do i=1,nbof_sites
         associate(s=>sites(i))
            call skip_rows(funit=1,nskip=2)
            call init_poly(s%lambda_elst,s%s_elst,s%b_elst,&
                           read_length=.true.,degree=1)
         end associate
      end do
      read(1,'(a)') row
      backspace(1)
      if(len_trim(row).eq.0) call skip_rows(funit=1,nskip=1)
      !Read the delta1 damping factors
      call read_colon_int(funit=1,val=length,nskip=2)
      call skip_rows(funit=1,nskip=4)
      allocate(site_pairs(nbof_sites,nbof_sites))
      do l=1,length
         call get_triplet(1,i,j,val)
         site_pairs(i,j)%delta1=val
      end do
      read(1,'(/a)') row
      flex_asym=index(row,'flexible').ne.0
      !Read the A12 coefficients
      call read_colon_int(funit=1,val=length,nskip=3)
      if(length.gt.0) call skip_rows(funit=1,nskip=4)
      do l=1,length
         call get_triplet(1,i,j,val)
         site_pairs(i,j)%A12=val
      end do
      !Read the delta(4+2k) (k=3,2,1,0) damping factors
      do k=3,0,-1
         call read_colon_int(funit=1,val=length,nskip=3)
         if(length.gt.0) call skip_rows(funit=1,nskip=4)
         do l=1,length
            call get_triplet(1,i,j,val)
            associate(sp=>site_pairs(i,j))
               if(k.eq.3) then
                  sp%delta10=val
               elseif(k.eq.2) then
                  sp%delta8=val
               elseif(k.eq.1) then
                  sp%delta6=val
               elseif(k.eq.0) then
                  sp%delta4=val
               endif
            end associate
         end do
      end do
      !Read the C(4+2k) (k=3,2,1,0) coefficients
      call skip_rows(funit=1,nskip=1)
      if(flex_asym) then
         do k=3,0,-1
            call skip_rows(funit=1,nskip=2)
            do i=1,nbof_sites
               call skip_rows(funit=1,nskip=2)
               associate(s=>sites(i))
                  if(k.eq.3) then
                     call init_poly(&
                              s%lambda_disp10,s%s_disp10,&
                              s%b_disp10,read_length=.true.,degree=1)
                  elseif(k.eq.2) then
                     call init_poly(&
                              s%lambda_disp8,s%s_disp8,s%b_disp8,&
                              read_length=.true.,degree=1)
                  elseif(k.eq.1) then
                     call init_poly(&
                              s%lambda_disp6,s%s_disp6,s%b_disp6,&
                              read_length=.true.,degree=1)
                  elseif(k.eq.0) then
                     call init_poly(&
                              s%lambda_disp4,s%s_disp4,s%b_disp4,&
                              read_length=.true.,degree=1)
                  endif
               end associate
            end do
         end do
      else
         do k=3,0,-1
            call read_colon_int(funit=1,val=length,nskip=2)
            if(length.gt.0) call skip_rows(funit=1,nskip=4)
            do l=1,length
               call get_triplet(1,i,j,val)
               associate(sp=>site_pairs(i,j))
                  if(k.eq.3) then
                     sp%C10=val
                  elseif(k.eq.2) then
                     sp%C8=val
                  elseif(k.eq.1) then
                     sp%C6=val
                  elseif(k.eq.0) then
                     sp%C4=val
                  endif
               end associate
            end do
            call skip_rows(funit=1,nskip=1)
         end do
      endif
      !Read the beta(i,j) parameters
      call read_colon_int(funit=1,val=length,nskip=4)
      call skip_rows(funit=1,nskip=4)
      do l=1,length
         call get_triplet(1,i,j,val)
         site_pairs(i,j)%beta=val
      end do
      !Read the gamma(i,j) parameters
      call read_colon_int(funit=1,val=length,nskip=3)
      if(length.gt.0) call skip_rows(funit=1,nskip=4)
      do l=1,length
         call get_triplet(1,i,j,val)
         site_pairs(i,j)%gamma=val
      end do
      call skip_rows(funit=1,nskip=3)
      !Read the preexponential polynomials
      do
         read(1,'(a)') row
         associate(neqs=>count([(row(k:k).eq.'=',&
                               k=1,len_trim(row))]))
            if(neqs.ne.2) exit
         end associate
         associate(ieq1 => index(row,'='),&
                   icom => index(row,','),&
                   irp  => index(row,')'))
         associate(ieq2 => ieq1+index(row(ieq1+1:),'='))
            read(row(ieq1+1:icom-1),*) i
            read(row(ieq2+1:irp-1),*) j
         end associate
         end associate
         call skip_rows(funit=1,nskip=1)
         call order(i,j)
         associate(sp=>site_pairs(i,j))
            call init_poly(sp%lambda_exp,sp%s_exp,sp%b_exp,&
                           read_length=.true.,degree=1)
         end associate
         call skip_rows(funit=1,nskip=1)
      end do
      do i=1,nbof_sites-1
         do j=i+1,nbof_sites
            associate(sp=>site_pairs(i,j))
               if(.not.allocated(sp%lambda_exp)) &
                  call init_poly(&
                     sp%lambda_exp,sp%s_exp,sp%b_exp,&
                     read_length=.false.,degree=1)
            end associate
         end do
      end do
      !Read the alpha parameters
      call read_colon_int(funit=1,val=length,nskip=3)
      call skip_rows(funit=1,nskip=3)
      if(length.gt.0) call skip_rows(funit=1,nskip=1)
      do l=1,length
         read(1,*) i,sites(i)%alpha
      end do
      if(length.gt.0) call skip_rows(funit=1,nskip=3)
      !Read the delta3 parameters
      call read_colon_int(funit=1,val=length,nskip=0)
      if(length.gt.0) call skip_rows(funit=1,nskip=4)
      do l=1,length
         call get_triplet(1,i,j,val)
         site_pairs(i,j)%delta3=val
      end do
      do i=1,nbof_sites-1
         do j=i+1,nbof_sites
            site_pairs(j,i)=site_pairs(i,j)
         end do
      end do
      close(1)
      !***
      contains
      !*****************************************************************
         subroutine read_colon_int(funit,val,nskip)
      !*****************************************************************
         implicit none
         character(1000) :: &
            row
         integer :: &
            funit,&
            nskip,&
            val
         !***
         call skip_rows(funit,nskip)
         read(funit,'(a)') row
         associate(start=>index(row,':')+1)
            read(row(start:),*) val
         end associate
         end
      !*****************************************************************
         subroutine skip_rows(funit,nskip)
      !*****************************************************************
         implicit none
         character(1000) row
         integer :: &
            i,&
            funit,&
            nskip
         !***
         do i=1,nskip
            read(funit,'(a)') row
         end do
         !***
         end
      !*****************************************************************
         subroutine get_triplet(funit,i,j,val)
      !*****************************************************************
         implicit none
         character(1000) :: &
            str
         real(r64) :: &
            val
         integer :: &
            funit,&
            i,j
         !***
         read(funit,*) i,j,str
         call order(i,j)
         if(trim(str).eq.'d_huge') then
            val=huge(0.d0)
         else
            read(str,*) val
         endif
         !***
         end
      !*****************************************************************
         subroutine order(i,j)
      !*****************************************************************
         implicit none
         integer :: &
            i,j,&
            ii,jj
         !***
         ii=min(i,j)
         jj=max(i,j)
         i=ii
         j=jj
         !***
         end
      !*****************************************************************
         subroutine init_poly(lambda,s,b,read_length,degree)
      !*****************************************************************
         implicit none
         integer, allocatable :: &
            lambda(:,:)
         logical :: &
            read_length
         real(r64), allocatable :: &
            s(:,:),&
            b(:)
         integer :: &
            j,p,&
            nskip0,&
            degree,&
            length
         !***
         nskip0=0
         if(read_length) then
            call read_colon_int(funit=1,val=length,nskip=0)
            if(degree.eq.1) nskip0=1
         else
            length=0
         endif
         associate(eff_length=>max(length,1))
            allocate(lambda(3*degree,eff_length))
            allocate(s(degree,eff_length),source=0.d0)
            allocate(b(eff_length),source=0.d0)
            if(length.eq.0) then
               call skip_rows(funit=1,nskip=nskip0)
               lambda(:,1)=[([1,2,0],p=1,degree)]
            else
               call skip_rows(funit=1,nskip=4)
               do j=1,length
                  read(1,*) lambda(:,j),s(:,j),b(j)
               end do
            endif
         end associate
         !***
         end
      !***
      end
