!***********************************************************************
      module mod_rbalign
!***********************************************************************
      use iso_fortran_env, only: r64=>real64
      use mod_graph, only: bipartite_matching
      implicit none
      include "iface/matchtool_iface_inc.f90"
      private
      public :: config_rbalign, best_align, get_grmsd
      type rbaconfig_
         real(r64) :: &
            min_bias_tol=1.d-6,&
            rep_tol=1.d-6,&
            grmsd_tol=huge(0.d0)
         integer :: &
            mreps=10,&
            mcmcalls=huge(0)
         logical :: &
            biased=.true.,&
            mirror=.true.
         procedure(matchtool_iface), pointer, nopass :: &
            matchtool=>bipartite_matching
      end type
      type(rbaconfig_), save :: &
         rbaconfig
      contains
      !*****************************************************************
         subroutine &
               config_rbalign(grmsd_tol,min_bias_tol,rep_tol,&
                              mreps,mcmcalls,biased,mirror,matchtool)
      !*****************************************************************
         use mod_catch, only: catch_error
         implicit none
         real(r64), optional :: &
            grmsd_tol,&
            min_bias_tol,&
            rep_tol
         integer, optional :: &
            mreps,&
            mcmcalls
         logical, optional :: &
            biased,&
            mirror
         procedure(matchtool_iface), optional :: &
            matchtool
         !***
         if(present(grmsd_tol)) then
            call catch_error(&
                     err=grmsd_tol.lt.epsilon(0.d0).or.&
                         grmsd_tol.gt.huge(0.d0),&
                     msg='incorrect grmsd_tol value.',&
                     proc='config_rbalign')
            rbaconfig%grmsd_tol=grmsd_tol
         endif
         if(present(min_bias_tol)) then
            call catch_error(&
                     err=min_bias_tol.lt.epsilon(0.d0).or.&
                         min_bias_tol.gt.huge(0.d0),&
                     msg='incorrect min_bias_tol value.',&
                     proc='config_rbalign')
            rbaconfig%min_bias_tol=min_bias_tol
         endif
         if(present(rep_tol)) then
            call catch_error(&
                     err=rep_tol.lt.epsilon(0.d0).or.&
                         rep_tol.gt.huge(0.d0),&
                     msg='incorrect rep_tol value.',&
                     proc='config_rbalign')
            rbaconfig%rep_tol=rep_tol
         endif
         if(present(mreps)) then
            call catch_error(&
                     err=mreps.le.0,&
                     msg='mreps must be positive.',&
                     proc='config_rbalign')
            rbaconfig%mreps=mreps
         endif
         if(present(mcmcalls)) then
            call catch_error(&
                     err=mcmcalls.le.0,&
                     msg='mcmcalls must be positive.',&
                     proc='config_rbalign')
            rbaconfig%mcmcalls=mcmcalls
         endif
         if(present(biased)) &
            rbaconfig%biased=biased
         if(present(mirror)) &
            rbaconfig%mirror=mirror
         if(present(matchtool)) &
            rbaconfig%matchtool=>matchtool
         !***
         end
      !*****************************************************************
         function get_grmsd(eqindsA,cartsA,eqindsB,cartsB) result(grmsd)
      !*****************************************************************
         implicit none
         integer :: &
            eqindsA(:),&
            eqindsB(:)
         real(r64) :: &
            cartsA(:,:),&
            cartsB(:,:),&
            grmsd
         real(r64), allocatable :: &
            tmp_cartsB(:,:)
         !***
         tmp_cartsB=cartsB
         call best_align(eqindsA,cartsA,eqindsB,tmp_cartsB,grmsd)
         !***
         end
      !*****************************************************************
         subroutine best_align(eqindsA,cartsA,eqindsB,cartsB,grmsd)
      !*****************************************************************
         use mod_catch, only: catch_error
         use mod_euclid, only: kabsch, get_uni_rot_mat
         implicit none
         integer :: &
            eqindsA(:),&
            eqindsB(:),&
            npoints,&
            nreps,&
            ncmcalls,i
         real(r64) :: &
            cartsA(:,:),&
            cartsB(:,:),&
            rot(3,3),&
            tran(3),&
            grmsd,&
            rmsd,&
            lbgrmsd,&
            ran
         integer, allocatable :: &
            links(:,:),&
            bip_conn(:,:),&
            matches(:)
         real(r64), allocatable :: &
            weights(:),&
            gcartsB(:,:)
         !***
         npoints=size(cartsA,2)
         call catch_error(&
                  err=size(cartsB,2).ne.npoints,&
                  msg='size incompatibility detected.',&
                  proc='best_align')
         call reseed()
         bip_conn=&
            get_bip_conn(eqindsA,cartsA,eqindsB,cartsB)
         lbgrmsd=grmsd_lbound(eqindsA,cartsA,eqindsB,cartsB)
         if(lbgrmsd.gt.rbaconfig%grmsd_tol) then
            grmsd=lbgrmsd
            return
         endif
         if(rbaconfig%min_bias_tol.eq.-huge(0.d0)) then
            grmsd=huge(0.d0)
            return
         endif
         !***
         associate(mreps     => rbaconfig%mreps,&
                   grmsd_tol => rbaconfig%grmsd_tol,&
                   rep_tol   => rbaconfig%rep_tol,&
                   mcmcalls  => rbaconfig%mcmcalls,&
                   mirror    => rbaconfig%mirror)
            grmsd=huge(0.d0)
            ncmcalls=0
            outer: do
               rot=get_uni_rot_mat()
               cartsB=matmul(rot,cartsB)
               if(mirror) then
                  call random_number(ran)
                  if(ran.lt.0.5d0) cartsB=-cartsB
               endif
               do
                  if(ncmcalls.eq.mcmcalls) exit outer
                  ncmcalls=ncmcalls+1
                  call sparsify_exgraph(bip_conn,links,&
                                        cartsA,cartsB,weights)
                  call rbaconfig%matchtool(links,matches,weights)
                  call catch_error(&
                           err=.not.allocated(matches),&
                           msg='array matches is not allocated.',&
                           proc='best_align')
                  call catch_error(&
                           err=size(matches).ne.npoints.or.&
                               any(matches.eq.0),&
                           msg='unsolvable matching problem.',&
                           proc='best_align')
                  bip_conn=bip_conn(:,matches)
                  cartsB=cartsB(:,matches)
                  call kabsch(cartsA,cartsB,rot,tran)
                  cartsB=matmul(rot,cartsB)+spread(tran,2,npoints)
                  rmsd=norm2(cartsA-cartsB)/dsqrt(dble(npoints))
                  if(all(matches.eq.[(i,i=1,npoints)])) exit
               end do
               if(grmsd.gt.rmsd+rep_tol) then
                  grmsd=rmsd
                  gcartsB=cartsB
                  nreps=1
               elseif(dabs(grmsd-rmsd).lt.rep_tol) then
                  nreps=nreps+1
               endif
               if(nreps.ge.mreps.or.&
                  grmsd_tol.ne.huge(0.d0).and.&
                  grmsd.lt.grmsd_tol) exit
            end do outer
            cartsB=gcartsB
         end associate
         !***
         end
      !*****************************************************************
         subroutine reseed()
      !*****************************************************************
         implicit none
         integer :: &
            n
         integer, allocatable :: &
            seed(:)
         !***
         call random_seed(size=n)
         allocate(seed(n),source=1)
         call random_seed(put=seed)
         !***
         end
      !*****************************************************************
         function grmsd_lbound(eqindsA,cartsA,eqindsB,cartsB) result(lb)
      !*****************************************************************
         implicit none
         include 'iface/lapack_iface_inc.f90'
         real(r64) :: &
            cartsA(:,:),&
            cartsB(:,:),&
            lb,&
            gramA(3,3),&
            gramB(3,3),&
            centA(3),&
            centB(3),&
            eigvalsA(3),&
            eigvalsB(3),&
            work(15)
         integer :: &
            eqindsA(:),&
            eqindsB(:),&
            npoints,&
            info,&
            i
         !***
         npoints=size(cartsA,2)
         centA=sum(cartsA,2)/npoints
         centB=sum(cartsB,2)/npoints
         lb=0.d0
         associate(sA    => cartsA-spread(centA,2,npoints),&
                   sB    => cartsB-spread(centB,2,npoints),&
                   first => minval(eqindsA),&
                   last  => maxval(eqindsA))
            do i=first,last
               associate(mA => merge(1,0,eqindsA.eq.i),&
                         mB => merge(1,0,eqindsB.eq.i))
                  if(sum(mA).eq.0) cycle
                  gramA=matmul(sA*spread(mA,1,3),transpose(sA))
                  gramB=matmul(sB*spread(mB,1,3),transpose(sB))
               end associate
               call dsyev('N','U',3,gramA,3,eigvalsA,work,20,info)
               call dsyev('N','U',3,gramB,3,eigvalsB,work,20,info)
               eigvalsA=max(0.d0,eigvalsA)
               eigvalsB=max(0.d0,eigvalsB)
               lb=lb+norm2(dsqrt(eigvalsA)-dsqrt(eigvalsB))**2
            end do
         end associate
         lb=dsqrt(lb/npoints)
         !***
         end
      !*****************************************************************
         subroutine &
               sparsify_exgraph(bip_conn,links,cartsA,cartsB,weights)
      !*****************************************************************
         implicit none
         real(r64), optional :: &
            cartsA(:,:),&
            cartsB(:,:)
         real(r64), optional, allocatable :: &
            weights(:)
         integer, allocatable :: &
            links(:,:),&
            bip_conn(:,:)
         integer :: &
            npoints,&
            medges,&
            nedges,&
            i,&
            j
         !***
         npoints=size(bip_conn,1)
         medges=npoints**2
         links=spread([0,0],2,medges)
         if(present(weights)) &
            weights=spread(0.d0,1,medges)
         nedges=0
         do i=1,npoints
            do j=1,npoints
               if(bip_conn(i,j).eq.0) cycle
               nedges=nedges+1
               links(:,nedges)=[i,j]
               if(present(weights)) &
                  weights(nedges)=norm2(cartsA(:,i)-cartsB(:,j))**2
            end do
         end do
         links=links(:,:nedges)
         if(present(weights)) weights=weights(:nedges)
         !***
         end
      !*****************************************************************
         function get_bip_conn(eqindsA,cartsA,eqindsB,cartsB) &
                               result(bip_conn)
      !*****************************************************************
         use mod_catch, only: catch_error
         use mod_sort, only: get_num_sort_perm
         implicit none
         real(r64) :: &
            cartsA(:,:),&
            cartsB(:,:),&
            bias_tol,&
            bias_fac
         integer :: &
            eqindsA(:),&
            eqindsB(:),&
            npoints,&
            i,j
         integer, allocatable :: &
            links(:,:),&
            bip_conn(:,:),&
            matches(:),&
            p(:)
         real(r64), allocatable :: &
            fpA(:,:),&
            fpB(:,:)
         !***
         npoints=size(cartsA,2)
         associate(eA => spread(eqindsA,2,npoints),&
                   eB => spread(eqindsB,1,npoints))
            bip_conn=merge(1,0,eA.eq.eB)
            call sparsify_exgraph(bip_conn,links)
            call rbaconfig%matchtool(links,matches)
            call catch_error(&
                     err=any(matches.eq.0),&
                     msg='unsolvable matching problem.',&
                     proc='best_align')
            eqindsB=eqindsB(matches)
            cartsB=cartsB(:,matches)
         end associate
         associate(cA2 => spread(cartsA,2,npoints),&
                   cA3 => spread(cartsA,3,npoints),&
                   cB2 => spread(cartsB,2,npoints),&
                   cB3 => spread(cartsB,3,npoints))
            fpA=norm2(cA2-cA3,1)
            fpB=norm2(cB2-cB3,1)
            do i=1,npoints
               p=get_num_sort_perm(fpA(:,i),.true.,eqindsA)
               fpA(:,i)=fpA(p,i)
               p=get_num_sort_perm(fpB(:,i),.true.,eqindsB)
               fpB(:,i)=fpB(p,i)
            end do
         end associate
         if(rbaconfig%biased) then
           !rbaconfig%min_bias_tol=&
           !   2.d0*dsqrt(dble(npoints))*rbaconfig%grmsd_tol
            bias_tol=rbaconfig%min_bias_tol
            bias_fac=1.5d0
            do
               bip_conn=0
               do i=1,npoints
                  associate(eAi  => eqindsA(i),&
                            fpAi => fpA(:,i))
                     do j=1,npoints
                        associate(eBj  => eqindsB(j),&
                                  fpBj => fpB(:,j))
                           if(eAi.ne.eBj) cycle
                           if(any(dabs(fpAi-fpBj).gt.bias_tol)) cycle
                           bip_conn(i,j)=1
                        end associate
                     end do
                  end associate
               end do
               call sparsify_exgraph(bip_conn,links)
               call rbaconfig%matchtool(links,matches)
               if(size(matches).eq.npoints.and.&
                  all(matches.ne.0)) exit
              !rbaconfig%min_bias_tol=-huge(0.d0)
              !exit
               bias_tol=bias_fac*bias_tol
            end do
         endif
         !***
         end
      !***
      end
