!***********************************************************************
      module mod_point_group
!***********************************************************************
      use iso_fortran_env, only : r64=>real64
      use mod_euclid, only : get_rot_mat, get_refl_mat, cross_product
      implicit none
      private
      public :: pgop_, pgroup, detect_point_group, get_pgroup
      !***
      type pgop_
        real(r64), allocatable :: &
           v(:)
        integer, allocatable :: &
           pinv(:)
        integer :: grade=0
        integer :: degree=1
      end type
      type pgroup
         type(pgop_), allocatable :: &
            pgops(:)
         character(:), allocatable :: &
            pgname
      end type
      contains
      !*****************************************************************
         function get_pgroup(anums,carts,pgrp_tol) result(pg)
      !*****************************************************************
         implicit none
         integer :: &
            anums(:)
         type(pgroup) :: &
            pg
         real(r64) :: &
            carts(:,:)
         real(r64), optional :: &
            pgrp_tol
         !***
         call detect_point_group(anums,carts,pg%pgname,pg%pgops,pgrp_tol=pgrp_tol)
         !***
         end
      !*****************************************************************
         subroutine detect_point_group(anums,carts,pgname,pgops,&
                                       eqinds,pgrp_tol,sdir_tol,&
                                       sing_tol,mtries,dump)
      !*****************************************************************
         use mod_catch, only : catch_error
         implicit none
         type lview_
            type(pgop_), allocatable :: &
               pgops(:)
            integer, allocatable :: &
               eqinds(:)
            real(r64) :: &
               pgrp_tol=1.d-4,&
               sdir_tol=1.d-6,&
               sing_tol=1.d-5
            integer :: &
               mtries=20
            logical :: &
               dump=.false.
         end type
         type(lview_) :: &
            lview
         integer :: &
            anums(:)
         real(r64) :: &
            carts(:,:)
         type(pgop_), allocatable, optional :: &
            pgops(:)
         integer, allocatable, optional :: &
            eqinds(:)
         real(r64), optional :: &
            pgrp_tol,&
            sdir_tol,&
            sing_tol
         integer, optional :: &
            mtries
         logical, optional :: &
            dump
         character(:), allocatable :: &
            gshape
         character(:), allocatable :: &
            pgname
         real(r64) :: &
            tol
         logical :: &
            closed
         integer :: &
            itry, none
         !***
         if(present(pgrp_tol)) lview%pgrp_tol=pgrp_tol
         if(present(sdir_tol)) lview%sdir_tol=sdir_tol
         if(present(sing_tol)) lview%sing_tol=sing_tol
         if(present(mtries)) lview%mtries=mtries
         if(present(dump)) lview%dump=dump
         call catch_error(&
                  err=size(anums).eq.0.or.&
                      size(anums).ne.size(carts,2),&
                  msg='inconsistent sizes.',&
                  proc='detect_point_group')
         associate(L=>lview)
            call realign_structure(anums,carts,gshape,L%sing_tol)
            tol=L%pgrp_tol
            do itry=1,L%mtries
               deallocate(L%pgops,stat=none)
               call find_equiv_atoms(anums,carts,L%eqinds,tol)
               call find_pgops(L%eqinds,carts,gshape,L%pgops,&
                               tol,L%sdir_tol)
               call find_point_group(L%pgops,pgname,closed=closed)
               if(closed) exit
               tol=0.5d0*tol
            end do
            if(L%dump) call print_pgops(anums,L%eqinds,carts,L%pgops)
            if(itry.gt.1) pgname=pgname//'?'
            if(present(eqinds)) eqinds=L%eqinds
            if(present(pgops)) pgops=L%pgops
         end associate
         end
      !*****************************************************************
         subroutine realign_structure(anums,carts,gshape,sing_tol)
      !*****************************************************************
         use mod_catch, only : catch_error
         use mod_pertab, only : iso_masses
         use mod_euclid, only : analyze_shape
         implicit none
         integer :: &
            anums(:)
         character(:), allocatable :: &
            gshape
         real(r64) :: &
            carts(:,:), sing_tol,&
            eigvects(3,3),&
            COM(3)
         real(r64), allocatable :: &
            trf_carts(:,:)
         integer :: &
            rank
         !***
         associate(natoms => size(anums),&
                   M      => iso_masses(anums))
            COM=matmul(carts,M)/sum(M)
            carts=carts-spread(COM,2,natoms)
            call analyze_shape(carts,rank,eigvects,sing_tol=sing_tol)
            if(natoms.eq.1) then
               gshape='spherical'
            else
               call catch_error(&
                        err=rank.eq.0,&
                        msg='invalid geometry passed.',&
                        proc='realign_structure')
               if(rank.eq.1) then
                  gshape='linear'
                  trf_carts=matmul(transpose(eigvects),carts)
                  if(norm2(carts-trf_carts).le.&
                     norm2(carts+trf_carts)) then
                     carts=trf_carts
                  else
                     carts=-trf_carts
                  endif
               else
                  gshape='nonlinear'
               endif
            endif
         end associate
         !***
         end
      !*****************************************************************
         subroutine find_equiv_atoms(anums,carts,eqinds,pgrp_tol)
      !*****************************************************************
         real(r64) :: &
            carts(:,:),&
            pgrp_tol
         integer ::&
            anums(:),&
            ind,&
            i,j,k,m
         integer, allocatable :: &
            eqinds(:)
         real(r64), allocatable :: &
            dists(:,:)
         logical, allocatable :: &
            free(:)
         !***
         associate(natoms=>size(anums))
            eqinds=spread(0,1,natoms)
            free=spread(.true.,1,natoms)
            dists=norm2(spread(carts,3,natoms)-&
                        spread(carts,2,natoms),1)
            ind=0
            do i=1,natoms
               associate(Zi => anums(i),&
                         Di => dists(:,i),&
                         Ei => eqinds(i))
                  if(Ei.ne.0) cycle
                  ind=ind+1
                  Ei=ind
                  do j=i+1,natoms
                     associate(Zj => anums(j),&
                               Dj => dists(:,j),&
                               Ej => eqinds(j))
                        if(Zi.ne.Zj) cycle
                        free=.true.
                        do k=1,natoms
                           associate(Zk  => anums(k),&
                                     Dik => Di(k))
                              m=minloc(dabs(Dik-Dj),dim=1,&
                                       mask=anums.eq.Zk.and.free)
                              if(dabs(Dik-Dj(m)).gt.pgrp_tol) exit
                              free(m)=.false.
                           end associate
                        end do
                        if(k.eq.natoms+1) Ej=ind
                     end associate
                  end do
               end associate
            end do
         end associate
         !***
         end
      !*****************************************************************
         subroutine find_pgops(eqinds,carts,gshape,pgops,&
                               pgrp_tol,sdir_tol)
      !*****************************************************************
         implicit none
         integer :: &
            eqinds(:)
         real(r64) :: &
            carts(:,:),&
            pgrp_tol,&
            sdir_tol
         character(*) :: &
            gshape
         type(pgop_), allocatable :: &
            pgops(:)
         type sdir
            real(r64) :: v(3)
         end type
         type(sdir), allocatable, target :: &
            sdirs(:)
         real(r64) :: &
               eY(3)=[0.d0,1.d0,0.d0],&
               eZ(3)=[0.d0,0.d0,1.d0]
         real(r64) :: &
            Pij(3), Pmij(3), Pjk(3),&
            cPiPj(3), cPijPjk(3)
         integer, parameter :: &
            inf=huge(0)
         logical :: &
            newrot,&
            skip_cross
         integer :: &
            natoms,&
            max_order,&
            pass,&
            i,j,k,m
         !***
         natoms=size(eqinds)
         call extend_pgops(grade=1,v=eZ)
         call extend_pgops(grade=-2,v=eZ)
         if(gshape.eq.'spherical'.or.gshape.eq.'linear') then
            call extend_pgops(grade=-inf,v=eZ)
            call extend_pgops(grade=inf,v=eZ)
            call extend_pgops(grade=-1,v=eZ)
            if(gshape.eq.'spherical') then
               call extend_pgops(grade=inf,v=eY)
            else
               call extend_pgops(grade=2,v=eY)
            endif
            call extend_pgops(grade=-1,v=eY)
         else
            max_order=maxval(&
               [(count(eqinds.eq.eqinds(i)),i=1,natoms)])
            skip_cross=.false.
            do pass=1,3
               do i=1,natoms
                  associate(Ei => eqinds(i),&
                            Pi => carts(:,i))
                     do j=1,i-1
                        associate(Ej => eqinds(j),&
                                  Pj => carts(:,j))
                           if(pass.eq.3) then
                              if(skip_cross) cycle
                              cPiPj=cross_product(Pi,Pj)
                              if(norm2(cPiPj).gt.pgrp_tol) then
                                 call extend_sdirs(v=cPiPj)
                                 skip_cross=.true.
                              endif
                              cycle
                           endif
                           if(Ei.ne.Ej) cycle
                           Pij=Pj-Pi
                           if(pass.eq.1) then
                              call extend_sdirs(v=Pij)
                              cycle
                           endif
                           Pmij=Pj+Pi
                           call extend_sdirs(v=Pmij)
                           do k=1,j-1
                              associate(Ek => eqinds(k),&
                                        Pk => carts(:,k))
                                 if(Ek.ne.Ei) cycle
                                 Pjk=Pk-Pj
                                 cPijPjk=cross_product(Pij,Pjk)
                                 call extend_sdirs(v=cPijPjk)
                              end associate
                           end do
                        end associate
                     end do
                  end associate
               end do
            end do
            do i=1,size(sdirs)
               associate(sdv=>sdirs(i)%v)
                  call extend_pgops(grade=-1,v=sdv)
                  do m=max_order,2,-1
                     call extend_pgops(grade=m,v=sdv,new=newrot)
                     if(newrot) exit
                  end do
                  do m=max_order,3,-1
                     call extend_pgops(grade=-m,v=sdv,new=newrot)
                     if(newrot) exit
                  end do
               end associate
            end do
         endif
         !***
         contains
         !**************************************************************
            subroutine extend_sdirs(v)
         !**************************************************************
            implicit none
            real(r64) :: &
               v(:)
            type(sdir), pointer :: &
               sdir_ptr
            real(r64) :: &
               dplus, dminus
            integer :: &
               i, none
            !***
            associate(vnorm=>norm2(v))
               if(vnorm.lt.pgrp_tol) return
               v=v/vnorm
            end associate
            allocate(sdirs(0),stat=none)
            associate(nsdirs=>size(sdirs))
               do i=1,nsdirs
                  sdir_ptr=>sdirs(i)
                  associate(vi=>sdir_ptr%v)
                     dplus=norm2(v-vi)
                     dminus=norm2(v+vi)
                     if(min(dplus,dminus).lt.sdir_tol) exit
                  end associate
               end do
               if(i.eq.nsdirs+1) then
                  allocate(sdir_ptr)
                  sdir_ptr%v=v
                  if(allocated(sdirs)) then
                     sdirs=[sdirs,sdir_ptr]
                  else
                     sdirs=[sdir_ptr]
                  endif
               endif
            end associate
            !***
            end
         !**************************************************************
            subroutine extend_pgops(grade,v,new) !come back here
         !**************************************************************
            integer :: &
               grade
            logical, optional :: &
               new
            type(pgop_), pointer :: &
               pgop_ptr
            real(r64) :: &
               v(:), trfo(3,3), refl(3,3),&
               phi, det
            real(r64), allocatable :: &
               trfo_carts(:,:)
            integer, allocatable :: &
               perm(:),&
               pinv(:)
            integer :: &
               degree, degree_max, red_grade,&
               red_degree, none, inv, i
            !***
            if(present(new)) new=.false.
            allocate(pgops(0),stat=none)
            associate(improp => grade.lt.0,&
                      order  => iabs(grade),&
                      inf    => huge(0))
               if(any(order.eq.[1,inf])) then
                  degree_max=1
               else
                  degree_max=order-1
               endif
               do degree=1,degree_max
                  associate(npgops => size(pgops),&
                            div    => gcd(grade,degree))
                     red_grade=grade/div
                     red_degree=degree/div
                     if(grade.ne.-2.and.red_grade.eq.-2) cycle
                     phi=360.d0*red_degree/iabs(red_grade)
                     trfo=get_rot_mat(v,phi)
                     if(improp) then
                        refl=get_refl_mat(v)
                        trfo=matmul(trfo,refl)
                     endif
                     trfo_carts=matmul(trfo,carts)
                     perm=get_sym_perm(trfo_carts)
                     if(size(perm).eq.0) then
                        if(degree.eq.1) then
                           return
                        else
                           cycle
                        endif
                     endif
                     if(present(new)) new=.true.
                     det=dot_product(trfo(:,1),&
                             cross_product(trfo(:,2),trfo(:,3)))
                     inv=merge(1,0,det.lt.0.d0)
                     pinv=[perm,inv]
                     if(gshape.eq.'nonlinear') then
                        do i=1,npgops
                           associate(pgop=>pgops(i))
                              if(size(pgop%pinv).ne.size(pinv)) cycle
                              if(all(pgop%pinv.eq.pinv)) exit
                           end associate
                        end do
                        if(i.le.npgops) cycle
                     endif
                     allocate(pgop_ptr)
                     pgop_ptr%v=v
                     pgop_ptr%pinv=pinv
                     pgop_ptr%grade=red_grade
                     pgop_ptr%degree=red_degree
                     if(allocated(pgops)) then
                        pgops=[pgops,pgop_ptr]
                     else
                        pgops=[pgop_ptr]
                     endif
                  end associate
               end do
            end associate
            !***
            end
         !**************************************************************
            function gcd(m,k) result(d)
         !**************************************************************
            implicit none
            integer :: &
               m, k, d, g
            !***
            d=m
            g=k
            do while(g.ne.0)
               associate(t=>mod(d,g))
                  d=g
                  g=t
               end associate
            end do
            d=iabs(d)
            !***
            end
         !**************************************************************
            function get_sym_perm(trfo_carts) result(perm)
         !**************************************************************
            real(r64) :: &
               trfo_carts(:,:)
            integer, allocatable :: &
               perm(:)
            logical, allocatable :: &
               free(:)
            integer &
               i,j
            !***
            perm=spread(0,1,natoms)
            free=spread(.true.,1,natoms)
            do i=1,natoms
               associate(Pi => carts(:,i),&
                         Ei => eqinds(i))
                  do j=1,natoms
                     associate(Qj   => trfo_carts(:,j),&
                               Ej   => eqinds(j),&
                               Fj   => free(j))
                        if(Ei.ne.Ej.or..not.Fj) cycle
                        if(norm2(Pi-Qj).lt.pgrp_tol) exit
                     end associate
                  end do
                  if(j.eq.natoms+1) then
                     perm=[integer::]
                     return
                  endif
                  free(j)=.false.
                  perm(j)=i
               end associate
            end do
            !***
            end
         !***
         end
      !*****************************************************************
         subroutine find_point_group(pgops,pgname,closed)
      !*****************************************************************
         use mod_string, only : to_string
         implicit none
         type(pgop_) :: pgops(:)
         logical :: &
            closed
         character(:), allocatable :: &
            pgname, smax_order, smax_rot_order
         integer :: &
            pgsize, max_order, max_rot_order,&
            nInv, nSigma, nInf, nC2, nC3, nC4, nC5
         !***
         associate(grades  => pgops%grade,&
                   degrees => pgops%degree,&
                   inf     => huge(0))
            max_order=maxval(iabs(grades))
            max_rot_order=maxval(grades)
            nInv=count(grades.eq.-2)
            nSigma=count(grades.eq.-1)
            nInf=count(grades.eq.inf)
            nC2=count(grades.eq.2.and.degrees.eq.1)
            nC3=count(grades.eq.3.and.degrees.eq.1)
            nC4=count(grades.eq.4.and.degrees.eq.1)
            nC5=count(grades.eq.5.and.degrees.eq.1)
            !convert "max_order" to string
            smax_order=''
            if(max_order.ne.inf) &
               smax_order=to_string(max_order)
            !convert "max_rot_order" to string
            smax_rot_order=''
            if(max_rot_order.ne.inf) &
               smax_rot_order=to_string(max_rot_order)
            !K group
            if(nInf.gt.1) then
               pgname='Kh'
               pgsize=7
            !linear (Civ/Dih) groups
            elseif(nInf.eq.1) then
               !Dih
               if(nInv.eq.1) then
                  pgname='Dih'
                  pgsize=7
               !Civ (L would be better)
               else
                  pgname='Civ'
                  pgsize=3
               endif
            !I groups
            elseif(nC5.gt.2) then
               !Ih
               if(nInv.eq.1) then
                  pgname='Ih'
                  pgsize=120
               !I
               else
                  pgname='I'
                  pgsize=60
               endif
            !O groups
            elseif(nC4.gt.2) then
               !Oh
               if(nInv.eq.1) then
                  pgname='Oh'
                  pgsize=48
               !O
               else
                  pgname='O'
                  pgsize=24
               endif
            !T groups
            elseif(nC3.gt.2) then
               !Th
               if(nInv.eq.1) then
                  pgname='Th'
                  pgsize=24
               !Td
               elseif(nSigma.gt.0) then
                  pgname='Td'
                  pgsize=24
               !T
               else
                  pgname='T'
                  pgsize=12
               endif
            !D groups
            elseif(nC2.gt.1) then
               !Dn
               if(nSigma.eq.0) then
                  pgname='D'//smax_rot_order
                  pgsize=2*max_rot_order
               !Dnd
               elseif(nSigma.eq.max_rot_order) then
                  pgname='D'//smax_rot_order//'d'
                  pgsize=4*max_rot_order
               !Dnh
               else!if(nSigma.eq.max_order+1) then
                  pgname='D'//smax_order//'h'
                  pgsize=4*max_rot_order
               endif
            !S groups (including Ci = S2)
            elseif(max_order.gt.max_rot_order) then
               !Sn
               pgname='S'//smax_order
               pgsize=max_order
            !C groups
            else
               !Cn
               if(nSigma.eq.0) then
                  pgname='C'//smax_order
                  pgsize=max_order
               !Cs/Cnh
               elseif(nSigma.eq.1) then
                  !Cs
                  if(max_rot_order.eq.1) then
                     pgname='Cs'
                     pgsize=2
                  !Cnh
                  else
                     pgname='C'//smax_order//'h'
                     pgsize=2*max_order
                  endif
               !Cnv
               else!if(nSigma.eq.max_order) then
                  pgname='C'//smax_order//'v'
                  pgsize=2*max_order
               endif
            endif
         end associate
         if(pgname.eq.'S2') pgname='Ci'
         !***
         closed=pgsize.eq.size(pgops)
         end
      !*****************************************************************
         subroutine print_pgops(anums,eqinds,carts,pgops)
      !*****************************************************************
         implicit none
         real(r64) :: &
            carts(:,:)
         integer :: &
            anums(:),&
            eqinds(:)
         type(pgop_), target :: &
            pgops(:)
         integer :: &
            natoms
         type(pgop_), pointer :: &
            pgo_ptr
         character(5) :: &
            label
         integer :: &
            i, order, quasi_grade
         !***
         natoms=size(anums)
         write(*,'(/a/)') 'Geometry:'
         do i=1,natoms
            write(*,'(i5,3x,3f12.6,i5)') &
               anums(i),carts(:,i),eqinds(i)
         end do
         write(*,'(/a/)') 'Point-group operations:'
         do i=1,size(pgops)
            pgo_ptr=>pgops(i)
            associate(v      => pgo_ptr%v,&
                      inf    => huge(0),&
                      grade  => pgo_ptr%grade,&
                      degree => pgo_ptr%degree)
               order=iabs(grade)
               quasi_grade=merge(0,grade,order.eq.inf)
               if(grade.lt.-1) then
                  if(grade.eq.-inf) then
                     label='inf'
                  else
                     write(label,'(i0)') order
                  endif
                  label='S'//label
               elseif(grade.eq.-2) then
                  label='i'
               elseif(grade.eq.-1) then
                  label='sigma'
               elseif(grade.eq.1) then
                  label='E'
               elseif(grade.gt.1) then
                  if(grade.eq.inf) then
                     label='inf'
                  else
                     write(label,'(i0)') order
                  endif
                  label='C'//label
               endif
               write(*,'(i5,3x,a,2(3x,i3)$)') &
                        i,label,quasi_grade,degree
               write(*,'(3x,3f12.6)') v
            end associate
         end do
         write(*,*)
         !***
         end
      end
