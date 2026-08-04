!***********************************************************************
      module mod_euclid
!***********************************************************************
      use iso_fortran_env, only : r64=>real64
      implicit none
      include 'iface/lapack_iface_inc.f90'
      private
      public :: &
         read_xyz, read_xyz_dynamic, write_xyz, cross_product,&
         normalize, orthog_proj, bond_angle, torsion_angle,&
         analyze_shape, get_rot_mat, get_refl_mat, std_axis_angle,&
         get_uni_rot_mat, wrap_angle, kabsch, gen_euler_decomp,&
         gen_spherical_decomp, gen_rigid_kernel, gen_rigid_trfo,&
         get_drot_mat, gen_rigid_jac
      interface cross_product
         module procedure cross_product_r
         module procedure cross_product_c
      end interface
      contains
      !*****************************************************************
      !  subroutine read_xyz(&
      !                 xyz_file,anums,symbs,carts,title,natoms,matoms)
      !*****************************************************************
#        include "dep/read_xyz_inc.f90"
        !end
      !*****************************************************************
      !  subroutine read_xyz_dynamic(&
      !                  xyz_file,anums,symbs,carts,title,natoms,matoms)
      !*****************************************************************
#        define DYNAMIC
#        include "dep/read_xyz_inc.f90"
        !end
      !*****************************************************************
         subroutine write_xyz(xyz_file,title,symbs,carts,append)
      !*****************************************************************
         use mod_string, only : expand_tilde
         implicit none
         logical :: &
            append
         character(*) :: &
            xyz_file,&
            title
         character(2) :: &
            symbs(:)
         real(r64) :: &
            carts(:,:)
         integer :: &
            i
         !***
         if(append) then
            open(1,file=expand_tilde(xyz_file),position='append')
         else
            open(1,file=expand_tilde(xyz_file))
         endif
         if(append) write(1,'(a)')
         associate(natoms=>size(symbs))
            write(1,'(i0)') natoms 
            write(1,'(a)') title
            do i=1,natoms
               write(1,'(a,3f15.6)') symbs(i),carts(:,i)
            end do
         end associate
         close(1)
         !***
         end
      !*****************************************************************
         function cross_product_r(x,y) result(res)
      !*****************************************************************
         implicit none
         real(r64) &
            x(3),&
            y(3),&
            res(3)
         !***
         res(1)=x(2)*y(3)-x(3)*y(2)
         res(2)=x(3)*y(1)-x(1)*y(3)
         res(3)=x(1)*y(2)-x(2)*y(1)
         !***
         end
      !*****************************************************************
         function cross_product_c(x,y) result(res)
      !*****************************************************************
         implicit none
         complex(r64) :: &
            x(3),&
            y(3),&
            res(3)
         !***
         res(1)=x(2)*y(3)-x(3)*y(2)
         res(2)=x(3)*y(1)-x(1)*y(3)
         res(3)=x(1)*y(2)-x(2)*y(1)
         !***
         end
      !*****************************************************************
         function normalize(v,tol) result(nv)
      !*****************************************************************
         use mod_catch, only : catch_error
         implicit none
         type lview_
            real(r64) :: tol=1.d-6
         end type
         type(lview_) :: &
            lview
         real(r64), optional :: &
            tol
         real(r64) :: &
            v(:)
         real(r64), allocatable :: &
            nv(:)
         real(r64) :: &
            vnorm
         !***
         if(present(tol)) lview%tol=tol
         vnorm=norm2(v)
         call catch_error(&
                  err=vnorm.lt.lview%tol,&
                  msg='Zero-length vector detected.',&
                  proc='normalize')
         nv=v/vnorm
         !***
         end
      !*****************************************************************
         function orthog_proj(u,v,tol) result(oproj)
      !*****************************************************************
         implicit none
         type lview_
            real(r64) :: tol=1.d-10
         end type
         type(lview_) :: &
            lview
         real(r64), optional :: &
            tol
         real(r64) :: &
            u(:),&
            v(:)
         real(r64), allocatable :: &
            oproj(:)
         !***
         associate(d=>normalize(u,lview%tol))
            if(present(tol)) lview%tol=tol
            oproj=v-dot_product(v,d)*d
         end associate
         !***
         end
      !*****************************************************************
         function bond_angle(v1,v2,v3,tol) result(bondang)
      !*****************************************************************
         implicit none
         type lview_
            real(r64) :: tol=1.d-6
         end type
         type(lview_) :: &
            lview
         real(r64), optional :: &
            tol
         real(r64) :: &
            v1(3),&
            v2(3),&
            v3(3)
         real(r64) :: &
            bondang,&
            b
         !***
         if(present(tol)) lview%tol=tol
         associate(rad2deg => 180.d0/dacos(-1.d0),&
                   r       => normalize(v1-v2,lview%tol),&
                   q       => normalize(v3-v2,lview%tol))
            b=dot_product(r,q)
            b=max(-1.d0,min(b,1.d0))
            bondang=rad2deg*dacos(b)
         end associate
         !***
         end
      !*****************************************************************
         function torsion_angle(v1,v2,v3,v4,tol) result(torang)
      !*****************************************************************
         implicit none
         type lview_
            real(r64) :: tol=1.d-6
         end type
         type(lview_) :: &
            lview
         real(r64), optional :: &
            tol
         real(r64) :: &
            v1(3),&
            v2(3),&
            v3(3),&
            v4(3),&
            t(3),&
            pt(3),&
            s(3),&
            u(3)
         real(r64) :: &
            torang,&
            a,c
         !***
         if(present(tol)) lview%tol=tol
         associate(rad2deg => 180.d0/dacos(-1.d0),&
                   r       => v1-v2,&
                   q       => v3-v2,&
                   p       => v4-v3)
            t=normalize(q,lview%tol)
            pt=orthog_proj(t,p,lview%tol)
            s=normalize(pt,lview%tol)
            u=cross_product(s,t)
            a=dot_product(r,s)
            c=dot_product(r,u)
            torang=rad2deg*datan2(c,a)
         end associate
         !***
         end
      !*****************************************************************
         subroutine analyze_shape(carts,rank,eigvects,avg_shift,&
                                  sing_tol)
      !*****************************************************************
         use mod_catch, only : catch_error
         implicit none
         include 'iface/lapack_iface_inc.f90'
         type lview_
            real(r64) :: sing_tol=1.d-5
            logical :: avg_shift=.false.
         end type
         type(lview_) :: &
            lview
         real(r64), optional :: &
            eigvects(3,3),&
            sing_tol
         integer, optional :: &
            rank
         logical, optional :: &
            avg_shift
         real(r64) :: &
            carts(:,:),&
            eigvals(3),&
            cov(3,3),&
            work(8),&
            avg(3)
         real(r64), allocatable :: &
            rel_carts(:,:)
         integer :: &
            info,&
            n,m
         !***
         if(present(sing_tol)) lview%sing_tol=sing_tol
         if(present(avg_shift)) &
            lview%avg_shift=avg_shift
         m=size(carts,1)
         n=size(carts,2)
         call catch_local_errors('shape')
         rel_carts=carts
         if(lview%avg_shift) then
            avg=sum(carts,dim=2)/n
            rel_carts=rel_carts-spread(avg,2,n)
         endif
         cov=matmul(rel_carts,transpose(rel_carts))/n
         call dsyev('V','U',3,cov,3,eigvals,work,8,info)
         call catch_local_errors('dsyev')
         if(present(eigvects)) eigvects=cov
         if(present(rank)) rank=count(eigvals.gt.lview%sing_tol)
         !***
         contains
         !**************************************************************
            subroutine catch_local_errors(opt)
         !**************************************************************
            use mod_catch, only : catch_error
            use mod_string, only : to_string
            implicit none
            character(*) :: &
               opt
            !***
            select case(opt)
               case('shape')
                  if(m.ne.3.or.n.le.0) &
                     call catch_error(&
                              err=.true.,&
                              msg='invalid shape ('//&
                                  to_string(m)//','//&
                                  to_string(n)//')'//&
                                  'provided for array "a".',&
                              proc='get_shape_rank',&
                              comment='valid shapes are (3,n>0).')
               case('dsyev')
                  call catch_error(&
                           err=info.ne.0,&
                           msg='dsyev crashed.',&
                           proc='get_shape_rank')
               case default
                  call catch_error(&
                           err=.true.,&
                           msg='unknown catch option.',&
                           proc='get_shape_rank')
            end select
            end
         !***
         end
      !*****************************************************************
         function get_rot_mat(axis,ang,tol) result(rot)
      !*****************************************************************
         implicit none
         type lview_
            real(r64) :: tol=1.d-12
         end type
         type(lview_) :: &
            lview
         real(r64), optional :: &
            tol
         real(r64) :: &
            axis(3),&
            rot(3,3),&
            K(3,3)
         real(r64) :: &
            ang,&
            rang
         integer :: &
            i
         !***
         if(present(tol)) lview%tol=tol
         associate(deg2rad => dacos(-1.d0)/180.d0,&
                   d       => normalize(axis,lview%tol))
            rang=deg2rad*ang
            K(1,:)=(/0.d0,-d(3),d(2)/)
            K(2,:)=(/d(3),0.d0,-d(1)/)
            K(3,:)=(/-d(2),d(1),0.d0/)
            rot=dsin(rang)*K+(1.d0-dcos(rang))*matmul(K,K)
            do i=1,3
               rot(i,i)=rot(i,i)+1.d0
            end do
         end associate
         !***
         end
      !*****************************************************************
         function get_refl_mat(vec,tol) result(refl)
      !*****************************************************************
         implicit none
         type lview_
            real(r64) :: tol=1.d-12
         end type
         type(lview_) :: &
            lview
         real(r64), optional :: &
            tol
         real(r64) :: &
            vec(3),&
            refl(3,3),&
            dresh(3,1)
         integer :: &
            i
         !***
         if(present(tol)) lview%tol=tol
         associate(d=>normalize(vec,lview%tol))
            dresh=reshape(d,[3,1])
            refl=-2.d0*matmul(dresh,transpose(dresh))
            do i=1,3
               refl(i,i)=refl(i,i)+1.d0
            end do
         end associate
         !***
         end
      !*****************************************************************
         subroutine std_axis_angle(axis,angle,tol)
      !*****************************************************************
         implicit none
         type lview_
            real(r64) :: tol=1.d-6
         end type
         type(lview_) :: &
            lview
         real(r64), optional :: &
            tol
         real(r64) :: &
            axis(3),&
            angle
         integer &
            idx
         !***
         if(present(tol)) &
            lview%tol=tol
         idx=findloc(dabs(axis).gt.lview%tol,value=.true.,dim=1)
         if(idx.eq.0) return
         if(axis(idx).lt.0.d0) then
            axis=-axis
            angle=-angle
         endif
         !***
         end         
      !*****************************************************************
         function get_uni_rot_mat() result(rot)
      !*****************************************************************
         implicit none
         real(r64), parameter :: &
            pi=dacos(-1.d0),&
            rad2deg=180.d0/pi
         real(r64) :: &
            rands(3),&
            rot(3,3),&
            axis(3),&
            ang
         !***
         call random_number(rands)
         associate(z   => 2.d0*rands(1)-1.d0,&
                   phi => 2.d0*pi*rands(2),&
                   w   => 1.d0-2.d0*dsqrt(rands(3)))
            axis(1:2)=dsqrt(1.d0-z*z)*[dcos(phi),dsin(phi)]
            axis(3)=z
            ang=rad2deg*dacos(w)
            rot=get_rot_mat(axis,ang)
         end associate
         !***
         end         
      !*****************************************************************
         subroutine wrap_angle(ang)
      !*****************************************************************
         implicit none
         real(r64) :: &
            ang
         ang=modulo(ang+180.d0,360.d0)-180.d0
         end
      !*****************************************************************
         subroutine kabsch(a,b,r,t)
      !*****************************************************************
         implicit none
         real(r64), optional :: &
            t(3)
         real(r64) :: &
            a(:,:),&
            b(:,:),&
            ta(3),&
            tb(3),&
            h(3,3),&
            u(3,3),&
            vT(3,3),&
            p(3),&
            s(3),&
            r(3,3),&
            work(15)
         real(r64) :: & 
            f
         integer :: &
            n,m,info
         !***
         call catch_local_errors('mismatch')
         m=size(a,1)
         n=size(a,2)
         call catch_local_errors('shape')
         if(present(t)) then
            ta=sum(a,dim=2)/n
            tb=sum(b,dim=2)/n
         else
            ta=0.d0
            tb=0.d0
         endif
         h=matmul(a-spread(ta,2,n),transpose(b-spread(tb,2,n)))
         call dgesvd('A','A',3,3,h,3,s,u,3,vT,3,work,15,info)
         call catch_local_errors('dgesvd')
         r=matmul(u,vT)
         p=cross_product(r(:,1),r(:,2))
         f=dot_product(r(:,3),p)-1.d0
         call dger(3,3,f,u(1,3),1,vT(3,1),3,r,3)
         if(present(t)) t=ta-matmul(r,tb)
         !***
         contains
         !**************************************************************
            subroutine catch_local_errors(opt)
         !**************************************************************
            use mod_catch, only : catch_error
            use mod_string, only : to_string
            implicit none
            character(*) :: &
               opt
            !***
            select case(opt)
               case('mismatch')
                  call catch_error(&
                           err=any(shape(a).ne.shape(b)),&
                           msg='shape mismatch between arrays '//&
                               '"a" and "b".',&
                           proc='kabsch')
               case('shape')
                  if(m.ne.3.or.n.le.0) &
                     call catch_error(&
                              err=.true.,&
                              msg='invalid shape ('//&
                                  to_string(m)//','//&
                                  to_string(n)//')'//&
                                  'provided for array "a".',&
                              proc='kabsch',&
                              comment='valid shapes are (3,n>0).')
               case('dgesvd')
                  call catch_error(&
                           err=info.ne.0,&
                           msg='dgesvd crashed.',&
                           proc='kabsch')
               case default
                  call catch_error(&
                           err=.true.,&
                           msg='unknown catch option.',&
                           proc='kabsch')
            end select
            !***
            end
         !***
         end
      !*****************************************************************
         subroutine &
                  gen_euler_decomp(R,axes,omega1,omega2,omega3,sing,tol)
      !*****************************************************************
         implicit none
         type lview_
            real(r64) :: tol=1.d-10
            logical :: sing
         end type
         type(lview_) :: &
            lview
         real(r64), optional :: &
            tol
         logical, optional :: &
            sing
         real(r64) :: &
            R(3,3),&
            axes(3,3),&
            Rd3(3),&
            R2(3,3),&
            Rtd2(3),&
            Rtd1(3),&
            R2d3(3),&
            R2td1(3)
         real(r64) :: &
            sigma,&
            omega1,&
            omega2,&
            omega3
         !***
         if(present(tol)) lview%tol=tol
         associate(d1 => normalize(axes(:,1),lview%tol),&
                   d2 => normalize(axes(:,2),lview%tol),&
                   d3 => normalize(axes(:,3),lview%tol))
            Rd3=matmul(R,d3)
            sigma=dot_product(d1,Rd3)
            omega2=rotang(d2,d3,d1,lview%tol,Rd3)
            lview%sing=dabs(dabs(sigma)-1.d0).lt.lview%tol
            if(lview%sing) then
               Rtd2=matmul(transpose(R),d2)
               omega1=dsign(1.d0,sigma)*rotang(d3,Rtd2,d2,lview%tol)
               omega3=0.d0
            else
               R2=get_rot_mat(d2,omega2)
               Rtd1=matmul(transpose(R),d1)
               R2d3=matmul(R2,d3)
               R2td1=matmul(transpose(R2),d1)
               omega1=rotang(d1,R2d3,Rd3,lview%tol)
               omega3=rotang(d3,Rtd1,R2td1,lview%tol)
            endif
            if(present(sing)) sing=lview%sing
         end associate
         !***
         end
      !*****************************************************************
         subroutine gen_spherical_decomp(t,axes,omega1,omega2,rho,tol)
      !*****************************************************************
         implicit none
         type lview_
            real(r64) :: tol=1.d-6
            logical :: sing
         end type
         type(lview_) :: &
            lview
         real(r64), optional :: &
            tol
         real(r64) :: &
            t(3),&
            axes(3,3),&
            R2(3,3),&
            Rd3(3),R2d3(3),&
            rho,&
            omega1,&
            omega2
         !***
         if(present(tol)) lview%tol=tol
         rho=norm2(t)
         if(rho.lt.lview%tol) then
            omega1=0.d0
            omega2=0.d0
         else
            associate(d1 => normalize(axes(:,1),lview%tol),&
                      d2 => normalize(axes(:,2),lview%tol),&
                      d3 => normalize(axes(:,3),lview%tol))
               Rd3=normalize(t,lview%tol)
               omega2=rotang(d2,d3,d1,lview%tol,Rd3)
               R2=get_rot_mat(d2,omega2)
               R2d3=matmul(R2,d3)
               omega1=rotang(d1,R2d3,Rd3,lview%tol)
            end associate
         endif
         !***
         end
      !*****************************************************************
         function rotang(d,w1,w2,tol,w3) result(omega)
      !*****************************************************************
         use mod_catch, only : catch_error
         implicit none
         real(r64), optional :: &
            w3(3),&
            tol
         type lview_
            real(r64) :: tol=1.d-6
         end type
         type(lview_) :: &
            lview
         real(r64) :: &
            d(3),&
            w1(3),w2(3),&
            omega,&
            a,b,c,&
            omega0,&
            disc,s,&
            omega_corr
         !***
         if(present(tol)) lview%tol=tol
         associate(rad2deg => 180.d0/dacos(-1.d0),&
                   dcw1    => cross_product(d,w1),&
                   dcw2    => cross_product(d,w2))
            a=dot_product(dcw1,dcw2)
            b=dot_product(dcw1,w2)
            omega0=rad2deg*datan2(b,a)
            omega=omega0
            if(.not.present(w3)) return
            c=a+dot_product(w3-w1,w2)
            disc=a**2+b**2-c**2
            call catch_error(&
                     err=disc.lt.-lview%tol,&
                     msg='nonreal rotation angle found.',&
                     proc='rotang')
            s=dsqrt(dabs(disc))
            omega_corr=rad2deg*datan2(s,c)
            omega=omega0+omega_corr !simply choose +omega_corr
            omega=modulo(omega+180.d0,360.d0)-180.d0
         end associate
         !***
         end
      !*****************************************************************
         subroutine gen_rigid_kernel(axes,rcoors,rot,tran,crot)
      !*****************************************************************
         implicit none
         real(r64) :: &
            axes(3,6),&
            rcoors(6),&
            rot(3,3),&
            tran(3)
         real(r64), optional :: &
            crot(3,3)
         real(r64) :: &
            rots(3,3,5)
         integer :: &
            i
         !***
         do i=1,5
            rots(:,:,i)=get_rot_mat(axes(:,i),rcoors(i))
         end do
         associate(rot1  => rots(:,:,1),&
                   rot2  => rots(:,:,2),&
                   rot3  => rots(:,:,3),&
                   rot4  => rots(:,:,4),&
                   rot5  => rots(:,:,5),&
                   axis6 => axes(:,6),&
                   rho   => rcoors(6))
            rot=matmul(rot1,matmul(rot2,rot3))
            if(present(crot)) then
               crot=matmul(rot4,rot5)
               tran=rho*matmul(crot,axis6)
            else
               tran=rho*matmul(matmul(rot4,rot5),axis6)
            endif
         end associate
         end
      !*****************************************************************
         function get_drot_mat(axis,ang) result(drot)
      !*****************************************************************
         implicit none
         real(r64) :: &
            axis(3),&
            ang,&
            drot(3,3),&
            K(3,3),&
            rang
         !***
         associate(deg2rad => dacos(-1.d0)/180.d0,&
                   d       => normalize(axis,1.d-12))
            rang=deg2rad*ang
            K(1,:)=[0.d0,-d(3),d(2)]
            K(2,:)=[d(3),0.d0,-d(1)]
            K(3,:)=[-d(2),d(1),0.d0]
            drot=deg2rad*(dcos(rang)*K+dsin(rang)*matmul(K,K))
         end associate
         !***
         end
      !*****************************************************************
         subroutine gen_rigid_jac(axes,rcoors,ref_carts,k,jac)
      !*****************************************************************
         use mod_catch, only : catch_error
         implicit none
         real(r64) :: &
            axes(3,6),&
            rcoors(6),&
            ref_carts(:,:),&
            jac(:,:)
         integer :: &
            k
         real(r64) :: &
            rots(3,3,5),&
            mat(3,3),&
            tvec(3)
         integer :: &
            i,n
         !***
         n=size(ref_carts,2)
         do i=1,5
            rots(:,:,i)=get_rot_mat(axes(:,i),rcoors(i))
         end do
         associate(rot1 => rots(:,:,1),&
                   rot2 => rots(:,:,2),&
                   rot3 => rots(:,:,3),&
                   rot4 => rots(:,:,4),&
                   rot5 => rots(:,:,5),&
                   rho  => rcoors(6),&
                   ax6  => axes(:,6))
            select case(k)
               case(1)
                  mat=matmul(get_drot_mat(axes(:,1),rcoors(1)),&
                             matmul(rot2,rot3))
                  jac=matmul(mat,ref_carts)
               case(2)
                  mat=matmul(rot1,&
                             matmul(get_drot_mat(axes(:,2),rcoors(2)),&
                                    rot3))
                  jac=matmul(mat,ref_carts)
               case(3)
                  mat=matmul(rot1,&
                             matmul(rot2,&
                                get_drot_mat(axes(:,3),rcoors(3))))
                  jac=matmul(mat,ref_carts)
               case(4)
                  tvec=rho*matmul(&
                        matmul(get_drot_mat(axes(:,4),rcoors(4)),rot5),&
                        ax6)
                  jac=spread(tvec,2,n)
               case(5)
                  tvec=rho*matmul(&
                        matmul(rot4,get_drot_mat(axes(:,5),rcoors(5))),&
                        ax6)
                  jac=spread(tvec,2,n)
               case(6)
                  tvec=matmul(matmul(rot4,rot5),ax6)
                  jac=spread(tvec,2,n)
               case default
                  call catch_error(err=.true.,&
                          msg='invalid local rcoor index.',&
                          proc='gen_rigid_jac')
            end select
         end associate
         !***
         end
      !*****************************************************************
         subroutine gen_rigid_trfo(&
                            mode,axes,rcoors,ref_carts,carts,tol,rig_tol)
      !*****************************************************************
         use mod_catch, only : catch_error
         implicit none
         type lview_
            real(r64) :: &
               tol=1.d-12
            real(r64) :: &
               rig_tol=1.d-4
         end type
         type(lview_) :: &
            lview
         real(r64), optional :: &
            tol,&
            rig_tol
         character(*) :: &
            mode
         real(r64) :: &
            axes(3,6),&
            rcoors(6),&
            ref_carts(:,:),&
            carts(:,:)
         real(r64) :: &
            rot(3,3),&
            tran(3)
         integer :: &
            n,m
         !***
         m=size(ref_carts,1)
         n=size(ref_carts,2)
         call catch_local_errors('arrays')
         if(present(tol)) lview%tol=tol
         if(present(rig_tol)) lview%rig_tol=rig_tol
         select case(mode)
            case('rcoors2carts')
               call gen_rigid_kernel(axes,rcoors,rot,tran)
               carts=matmul(rot,ref_carts)+spread(tran,2,n)
            case('carts2rcoors')
               call kabsch(carts,ref_carts,rot,tran)
               call catch_local_errors('grmsd')
               associate(rot_axes  => axes(:,:3),&
                         tran_axes => axes(:,4:))
                  call gen_euler_decomp(&
                           rot,rot_axes,omega1=rcoors(1),&
                           omega2=rcoors(2),omega3=rcoors(3),&
                           tol=lview%tol)
                  call gen_spherical_decomp(&
                           tran,tran_axes,omega1=rcoors(4),&
                           omega2=rcoors(5),rho=rcoors(6),&
                           tol=lview%tol)
               end associate
            case default
               call catch_local_errors('unknown')
         end select
         !***
         contains
         !**************************************************************
            subroutine catch_local_errors(opt)
         !**************************************************************
            use mod_string, only : to_string
            implicit none
            character(*) :: &
               opt
            real(r64), allocatable :: &
               fit_carts(:,:)
            real(r64) :: &
               grmsd
            !***
            select case(opt)
               case('arrays')
                  if(m.ne.3.or.n.le.0) &
                     call catch_error(&
                              err=.true.,&
                              msg='invalid shape ('//&
                                  to_string(m)//','//&
                                  to_string(n)//') provided for '//&
                                  'array "ref_carts".',&
                              proc='rigid_trfo')
                  if(mode.eq.'carts2rcoors') &
                     call catch_error(&
                              err=any(shape(ref_carts).ne.&
                                      shape(carts)),&
                              msg='shape mismatch between arrays '//&
                                  '"ref_carts" and "carts".',&
                              proc='rigid_trfo')
               case('grmsd')
                  fit_carts=matmul(rot,ref_carts)+spread(tran,2,n)
                  grmsd=norm2(carts-fit_carts)/dsqrt(dble(n))
                  if(grmsd.gt.lview%rig_tol) &
                     call catch_error(&
                              err=grmsd.gt.lview%rig_tol,&
                              msg='structural divergence exceeds '//&
                                  'the threshold (grmsd = '//&
                                  to_string(grmsd,2,&
                                            exp_form=.true.)//').',&
                              proc='rigid_trfo')
               case('unknown')
                  call catch_error(&
                           err=.true.,&
                           msg='unknown mode "'//mode//'".',&
                           proc='rigid_trfo')
               case('default')
                  call catch_error(&
                           err=.true.,&
                           msg='unknown option "'//opt//'".',&
                           proc='rigid_trfo')
            end select
            !***
            end
         !***
         end
      end
