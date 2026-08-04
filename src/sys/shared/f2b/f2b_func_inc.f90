   !********************************************************************
      function f2b_func(carts) result(U)
   !********************************************************************
      implicit none
      complex(r64) :: &
         U,&
         U_intra,&
         U_els,&
         U_disp,&
         U_exp,&
         U_ind
      real(r64) :: &
         damp_tol=1.d-15,&
         std_tol=1.d-12
      integer :: &
         maxiter=30
      complex(r64) :: &
         carts(:,:),&
         r(nbof_sites,nbof_sites),&
         q(nbof_sites),&
         C4(nbof_sites,nbof_sites),&
         C6(nbof_sites,nbof_sites),&
         C8(nbof_sites,nbof_sites),&
         C10(nbof_sites,nbof_sites),&
         eps(3,nbof_sites),&
         eta(3,nbof_sites),&
         omega(3,3,nbof_sites,nbof_sites),&
         mu0(3,nbof_sites),&
         mu(3,nbof_sites)
      integer :: &
         i,j,k,l,m
      !***
      !Calculate the entries of the distance matrix
      r=0.d0
      do i=1,nbof_sites
         do j=i+1,nbof_sites
            r(i,j)=scsqrt(sum((carts(:,i)-carts(:,j))**2))/bohrA
            r(j,i)=r(i,j)
         end do
      end do
      !Set term U_intra
      U_intra=poly2(Lambda_intra,S_intra,B_intra)
     !write(*,'(a,f20.10)') 'U_intra :',U_intra*kcalmolcm
      !Calculate the partial charges
      do i=1,nbof_sites
         associate(s=>sites(i))
            q(i)=poly1(s%lambda_elst,s%s_elst,s%b_elst)
         end associate
      end do
      !Set term U_els
      U_els=0.d0
      do i=1,nbof_sites
         associate(qi=>q(i))
            do j=i+1,nbof_sites
               associate(qj  => q(j),&
                         sp  => site_pairs(i,j),&
                         rij => r(i,j))
               associate(f1  => tang_toennies(1,sp%delta1,rij))
                  U_els=U_els+qi*qj*f1/rij
               end associate
               end associate
            end do
         end associate
      end do
      U_els=ehkcalmol*U_els
     !write(*,'(a,f20.10)') 'U_els   :',U_els*kcalmolcm
      !Calculate the disp coefficients
      if(flex_asym) then
         do i=1,nbof_sites
            associate(s     => sites(i),&
                      C4ii  => C4(i,i),&
                      C6ii  => C6(i,i),&
                      C8ii  => C8(i,i),&
                      C10ii => C10(i,i))
               C4ii=poly1(s%lambda_disp4,s%s_disp4,s%b_disp4)
               C6ii=poly1(s%lambda_disp6,s%s_disp6,s%b_disp6)
               C8ii=poly1(s%lambda_disp8,s%s_disp8,s%b_disp8)
               C10ii=poly1(s%lambda_disp10,s%s_disp10,s%b_disp10)
               do j=1,i-1
                  associate(C4jj  => C4(j,j),&
                            C6jj  => C6(j,j),&
                            C8jj  => C8(j,j),&
                            C10jj => C10(j,j),&
                            C4ij  => C4(i,j),&
                            C6ij  => C6(i,j),&
                            C8ij  => C8(i,j),&
                            C10ij => C10(i,j))
                     C4ij=scsqrt(C4ii*C4jj)
                     C6ij=scsqrt(C6ii*C6jj)
                     C8ij=scsqrt(C8ii*C8jj)
                     C10ij=scsqrt(C10ii*C10jj)
                     C4(j,i)=C4ij
                     C6(j,i)=C6ij
                     C8(j,i)=C8ij
                     C10(j,i)=C10ij
                  end associate
               end do
            end associate
         end do
      else
         c4=site_pairs%c4
         c6=site_pairs%c6
         c8=site_pairs%c8
         c10=site_pairs%c10
      endif
      !Set term U_disp
      U_disp=0.d0
      do i=1,nbof_sites
         do j=i+1,nbof_sites
            associate(sp    => site_pairs(i,j),&
                      rij   => r(i,j),&
                      C4ij  => C4(i,j),&
                      C6ij  => C6(i,j),&
                      C8ij  => C8(i,j),&
                      C10ij => C10(i,j))
            associate(f4    => tang_toennies(4,sp%delta4,rij),&
                      f6    => tang_toennies(6,sp%delta6,rij),&
                      f8    => tang_toennies(8,sp%delta8,rij),&
                      f10   => tang_toennies(10,sp%delta10,rij))
            associate(rec2  => rij**(-2))
            associate(rec4  => rec2**2)
            associate(rec6  => rec2*rec4)
            associate(rec8  => rec2*rec6)
            associate(rec10 => rec2*rec8)
            associate(rec12 => rec2*rec10)
               U_disp=U_disp-C4ij*f4*rec4
               U_disp=U_disp-C6ij*f6*rec6
               U_disp=U_disp-C8ij*f8*rec8
               U_disp=U_disp-C10ij*f10*rec10
               U_disp=U_disp+sp%A12*rec12
            end associate
            end associate
            end associate
            end associate
            end associate
            end associate
            end associate
            end associate
         end do
      end do
     !write(*,'(a,f20.10)') 'U_disp  :',U_disp*kcalmolcm
      !Set term U_exp
      U_exp=0.d0
      do i=1,nbof_sites
         do j=i+1,nbof_sites
            associate(sp  => site_pairs(i,j),&
                      rij => r(i,j))
            associate(h   => poly1(sp%lambda_exp,&
                             sp%s_exp,sp%b_exp))
               U_exp=U_exp+h*exp(sp%gamma-sp%beta*rij)
            end associate
            end associate
         end do
      end do
     !write(*,'(a,f20.10)') 'U_exp   :',U_exp*kcalmolcm
      !Set term U_ind
      omega=0.d0
      eps=0.d0
      do i=1,nbof_sites
         associate(alpha => sites(i)%alpha,&
                   pi    => carts(:,i),&
                   epsi  => eps(:,i))
            do j=1,nbof_sites
               associate(sp  => site_pairs(i,j),&
                         pj  => carts(:,j),&
                         rij => r(i,j))
               associate(f3  => tang_toennies(3,sp%delta3,rij))
                  if(abs(dreal(f3)).lt.damp_tol) cycle
                  associate(q0   => sites(j)%b_elst(1),&
                            f1   => tang_toennies(1,sp%delta1,rij),&
                            v    => (pi-pj)/bohrA,&
                            rec2 => rij**(-2))
                  associate(rec3 => rec2/rij)
                  associate(rec5 => rec2*rec3)
                        epsi=epsi+q0*f1*rec3*v
                        do k=1,3
                           associate(vk=>v(k))
                              do l=1,3
                                 associate(o => omega(k,l,i,j),&
                                           vl => v(l))
                                    o=3.d0*f3*rec5*vk*vl
                                    if(k.eq.l) o=o-f3*rec3
                                    o=alpha*o
                                 end associate
                              end do
                           end associate
                        end do
                  end associate
                  end associate
                  end associate
               end associate
               end associate
            end do
            eta(:,i)=alpha*epsi
         end associate
      end do
      mu0=0.d0
      do m=1,maxiter
         mu=eta
         do i=1,nbof_sites
            associate(mui=>mu(:,i))
               do j=1,nbof_sites
                  mui=mui+matmul(omega(:,:,i,j),mu0(:,j))
               end do
            end associate
         end do
         associate(std=>scsqrt(sum((mu-mu0)**2))/dsqrt(dble(nbof_sites)))
            if(dreal(std).lt.std_tol) exit
         end associate
         mu0=mu
      end do
      U_ind=-0.5d0*ehkcalmol*sum(eps*mu)
     !write(*,'(a,f20.10)') 'U_ind   :',U_ind*kcalmolcm
      U=U_intra+U_els+U_disp+U_exp+U_ind
      !***
      contains
      !*****************************************************************
         function poly2(Lambda,S,B) result(p2val)
      !*****************************************************************
         implicit none
         integer :: &
            Lambda(:,:)
         complex(r64) :: &
            p2val
         real(r64) :: &
            S(:,:),&
            B(:)
         integer :: &
            m
         !***
         p2val=0.d0
         associate(L_B=>size(B))
            do m=1,L_B
               associate(Lam1 => Lambda(1,m),&
                         Lam2 => Lambda(2,m),&
                         Lam3 => Lambda(3,m),&
                         Lam4 => Lambda(4,m),&
                         Lam5 => Lambda(5,m),&
                         Lam6 => Lambda(6,m),&
                         Bm   => B(m),&
                         S1m  => S(1,m),&
                         S2m  => S(2,m))
               associate(rL12 => r(Lam1,Lam2),&
                         rL45 => r(Lam4,Lam5))
                  p2val=p2val+Bm*(rL12-S1m)**Lam3*(rL45-S2m)**Lam6
               end associate
               end associate
            end do
         end associate
         !***
         end
      !*****************************************************************
         function poly1(lambda,s,b) result(p1val)
      !*****************************************************************
         implicit none
         integer :: &
            lambda(:,:)
         complex(r64) :: &
            p1val
         real(r64) :: &
            s(:,:),&
            b(:)
         integer :: &
            m
         !***
         p1val=0.d0
         associate(L_b=>size(b))
            do m=1,L_b
               associate(lam1 => Lambda(1,m),&
                         lam2 => Lambda(2,m),&
                         lam3 => Lambda(3,m),&
                         bm   => b(m),&
                         s1m  => s(1,m))
               associate(rL12 => r(Lam1,Lam2))
                  p1val=p1val+bm*(rL12-s1m)**lam3
               end associate
               end associate
            end do
         end associate
         !***
         end
      !*****************************************************************
         function tang_toennies(n,d,x) result(ttval)
      !*****************************************************************
         implicit none
         real(r64) :: &
            d
         complex(r64) :: &
            x,&
            ttval,&
            s,&
            t
         integer :: &
            n,&
            k
         !***
         if(d.eq.huge(0.d0)) then
            ttval=1.d0
         else
            associate(dx=>d*x)
               t=1.d0
               do k=1,n
                  t=t*dx/dble(k)
               end do
               s=0.d0
               k=n
               do
                  k=k+1
                  t=t*dx/dble(k)
                  s=s+t
                  if(dreal(t).le.std_tol*dreal(s)) exit
               end do
               ttval=exp(-dx)*s
            end associate
         endif
         !***
         end
      !***
      end
