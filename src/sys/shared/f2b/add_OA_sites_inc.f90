   !********************************************************************
      subroutine add_OA_sites(carts,tol,init)
   !********************************************************************
      implicit none
      logical :: &
         init
      complex(r64) :: &
         carts(:,:)
      real(r64) :: &
         tol
      integer :: &
         i,m
      !***
      do i=nbof_atoms+1,nbof_sites
         m=i-nbof_atoms
         associate(j => lst_OA(1,m),&
                   k => lst_OA(2,m),&
                   l => lst_OA(3,m))
         associate(w => w_OA(:,m),&
                   a => carts(:,i),&
                   b => carts(:,j),&
                   c => carts(:,k))
         associate(s => a-b,&
                   t => a-c,&
                   u => c-b)
            if(l.eq.0) then
               if(init) then
                  associate(snorm => scsqrt(sum(s**2)),&
                            tnorm => scsqrt(sum(t**2)))
                  associate(denom => snorm + tnorm)
                     if(dreal(denom).lt.tol) then
                        write(*,'(2x,a/)') &
                           'Error: coincident atomic '//&
                           'positions detected in '//&
                           'subroutine '//&
                           '"add_OA_sites".'
                        stop
                     endif
                     w(1)=dreal(tnorm/denom)
                     w(2:3)=0.d0
                  end associate
                  end associate
               else
                  a=b+w(1)*u
               endif
            else
               associate(d    => carts(:,l))
               associate(v    => d-b)
               associate(Cuv  => cross_product(u,v))
               associate(norm => scsqrt(sum(Cuv**2)))
                  if(dreal(norm).lt.tol) then
                     write(*,'(2x,a/)') &
                        'Error: collinear bonds '//&
                        'detected in subroutine '//&
                        '"add_OA_sites".'
                     stop
                  endif
                  associate(nCuv => Cuv/norm)
                     if(init) then
                        associate(Csv => cross_product(s,v),&
                                  Cus => cross_product(u,s))
                           w(1)=dreal(sum(nCuv*Csv)/norm)
                           w(2)=dreal(sum(nCuv*Cus)/norm)
                           w(3)=dreal(sum(nCuv*s))
                        end associate
                     else
                        a=b+w(1)*u+w(2)*v+w(3)*nCuv
                     endif
                  end associate
               end associate
               end associate
               end associate
               end associate
            endif
         end associate
         end associate
         end associate
      end do
      !***
      end
