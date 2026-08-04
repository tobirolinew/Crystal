!***********************************************************************
      function scsqrt(z)
!***********************************************************************
      implicit none
      complex(r64) :: &
         z,&
         scsqrt
      real(r64), parameter :: &
         split=134217729.d0
      real(r64) :: &
         x,y,ax,t,e,c,ahi,alo,bhi,blo,p,pe,rr,h,pc,pce,tt,hh
      !***
      x=dreal(z)
      y=dimag(z)
      ax=dabs(x)
      if(y.eq.0.d0) then
         if(x.ge.0.d0) then
            scsqrt=dcmplx(dsqrt(ax),0.d0)
         else
            scsqrt=dcmplx(0.d0,dsqrt(ax))
         endif
         return
      endif
      if(x.eq.0.d0) then
         t=dsqrt(0.5d0*dabs(y))
         scsqrt=dcmplx(t,dsign(t,y))
         return
      endif
      if(dabs(y).lt.1.d-9*ax) then
         t=dsqrt(ax)
         c=split*t
         ahi=c-(c-t)
         alo=t-ahi
         p=t*t
         pe=((ahi*ahi-p)+2.d0*ahi*alo)+alo*alo
         rr=(ax-p)-pe
         e=0.5d0*rr/t
         tt=t+e
         h=0.5d0*y/t
         c=split*h
         ahi=c-(c-h)
         alo=h-ahi
         c=split*t
         bhi=c-(c-t)
         blo=t-bhi
         pc=h*t
         pce=((ahi*bhi-pc)+ahi*blo+alo*bhi)+alo*blo
         hh=h+(((0.5d0*y-pc)-pce)-h*e)/tt
      else
         t=dsqrt(0.5d0*(dsqrt(x*x+y*y)+ax))
         tt=t
         hh=0.5d0*y/t
      endif
      if(x.ge.0.d0) then
         scsqrt=dcmplx(tt,hh)
      else
         scsqrt=dcmplx(dabs(hh),dsign(tt,y))
      endif
      !***
      end
