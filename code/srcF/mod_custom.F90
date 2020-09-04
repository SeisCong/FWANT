module custom_mod

! This module is used to custom the 3D geology model
!
! Author: Wei ZHANG     Email: zhangwei.zw@gmail.com
! Copyright (C) 2006 Wei ZHANG

!*****************************************************************************
!
! $Date: 2009-01-16 11:21:43 -0500 (Fri, 16 Jan 2009) $
! $Revision: 54 $
! $LastChangedBy: zhangw $
!
!*****************************************************************************

use constants_mod
implicit none
!private

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------
subroutine custom_media(x,y,z,vp,vs,dens,Qatt,Qf0,model)
real(SP),intent(in) :: x,y,z
real(SP),dimension(2),intent(out) :: Vp,Vs,dens,Qatt
real,intent(out) :: Qf0
character (len=*),intent(in) :: model
select case (trim(model))
case ('ak135')
   call custom_media_ak135(x,y,z,Vp,Vs,dens,Qatt,Qf0)
case ('ESG2006')
   call custom_media_ESG2006(x,y,z,Vp,Vs,dens,Qatt,Qf0)
case default
   print *, "there is no custom_media_"//trim(model)
   stop 1
end select
end subroutine custom_media

subroutine custom_media_ESG2006(x,y,z,vp,vs,dens,Qatt,Qf0)
real,intent(in) :: x,y,z
real,dimension(2),intent(out) :: Vp,Vs,dens,Qatt
real,intent(out) :: Qf0

!bechmark2 of ESG2006
real x0,y0,z0
real D

x0=x;y0=y;z0=z
z0=z0
D=max(-z0,0.0)

vs  = 300.0 + 19.0 * sqrt(D)
vp  = 1450.0 + 1.2 * D 
dens = 2140.0 + 0.125 * D 
!D: depth in m
!VS: S velocity in m/s
!VP: P velocity in m/s
!RHO: MASS DENSITY in kg/m^3


!ATTENUATION:
!Qmu=50
!Qkappa=infinite
!or
Qatt  = 50
!QP  = QS * (3/4) / ( VS^2/VP^2 ) 
!where Qmu=shear quality factor, Qkappa=bulk quality factor
end subroutine custom_media_ESG2006

subroutine custom_media_ak135(x,y,z,vp,vs,dens,Qatt,Qf0)
real,intent(in) :: x,y,z
real,dimension(2),intent(out) :: Vp,Vs,dens,Qatt
real,intent(out) :: Qf0

real,dimension(2,10) :: alpha,beta,rho,Q
real,dimension(10) :: H
integer,dimension(1) :: p
real L1,L2

integer n,k1,k2

Qf0=1.0

alpha(:,1 )=(/ 5800.0, 5800.0 /)
alpha(:,2 )=(/ 5800.0, 6500.0 /)
alpha(:,3 )=(/ 6500.0, 8040.0 /)
alpha(:,4 )=(/ 8040.0, 8040.0 /)
alpha(:,5 )=(/ 8050.0, 8050.0 /)
alpha(:,6 )=(/ 8175.0, 8175.0 /)
alpha(:,7 )=(/ 8300.0, 8300.0 /)
alpha(:,8 )=(/ 8482.5, 8482.5 /)
alpha(:,9 )=(/ 8665.0, 8665.0 /)
alpha(:,10)=(/ 8847.5, 8847.5 /)

beta(:,1 )=(/ 3460.0, 3460.0 /)
beta(:,2 )=(/ 3460.0, 3850.0 /)
beta(:,3 )=(/ 3850.0, 4480.0 /)
beta(:,4 )=(/ 4490.0, 4490.0 /)
beta(:,5 )=(/ 4500.0, 4500.0 /)
beta(:,6 )=(/ 4509.0, 4509.0 /)
beta(:,7 )=(/ 4518.0, 4523.0 /)
beta(:,8 )=(/ 4609.0, 4609.0 /)
beta(:,9 )=(/ 4696.0, 4696.0 /)
beta(:,10)=(/ 4783.0, 4783.0 /) 

rho(:,1 )=(/ 2720.0, 2720.0 /)
rho(:,2 )=(/ 2720.0, 2920.0 /)
rho(:,3 )=(/ 2920.0, 3319.8 /)
rho(:,4 )=(/ 3345.5, 3345.5 /)
rho(:,5 )=(/ 3371.3, 3371.3 /)
rho(:,6 )=(/ 3398.5, 3398.5 /)
rho(:,7 )=(/ 3425.8, 3425.8 /)
rho(:,8 )=(/ 3456.1, 3456.1 /)
rho(:,9 )=(/ 3486.4, 3486.4 /)
rho(:,10)=(/ 3516.7, 3516.7 /) 

Q(:,1 )=(/ 2000.0, 2000.0 /)
Q(:,2 )=(/ 2000.0, 2000.0 /)
Q(:,3 )=(/ 2000.0, 2000.0 /)
Q(:,4 )=(/ 2000.0, 2000.0 /)
Q(:,5 )=(/ 2000.0, 2000.0 /)
Q(:,6 )=(/ 2000.0, 2000.0 /)
Q(:,7 )=(/ 2000.0, 2000.0 /)
Q(:,8 )=(/ 2000.0, 2000.0 /)
Q(:,9 )=(/ 2000.0, 2000.0 /)
Q(:,10)=(/ 2000.0, 2000.0 /) 

H =(/ 0.0, -20.0, -35.0, -77.5, -120.0, -165.0, -210.0, -260.0, -310.0, -360.0 /)*1000

if (z>=H(1)) then
   Vp=alpha(2,1)
   Vs=beta (2,1)
   dens=rho(2,1)
   Qatt=Q  (2,1)
elseif (z<=H(10)) then
   Vp=alpha(1,10)
   Vs=beta (1,10)
   dens=rho(1,10)
   Qatt=Q  (1,10)
else
   p=minloc(abs(z-H))
   n=p(1)
   if (abs(z-H(n))<=SEIS_EQUAL) then
      Vp(1)=alpha(1,n); Vp(2)=alpha(2,n) 
      Vs(1)=beta (1,n); Vs(2)=beta (2,n) 
      dens(1)=rho(1,n); dens(2)=rho(2,n) 
      Qatt(1)=Q  (1,n); Qatt(2)=Q  (2,n) 
   else
      if (z>H(n)) then
         k1=n-1; k2=n
      else
         k1=n; k2=n+1
      end if
      L1=(z-H(k2))/(H(k1)-H(k2))
      L2=(z-H(k1))/(H(k2)-H(k1))
      Vp(:)   =alpha(2,k1)*L1+alpha(1,k2)*L2
      Vs(:)   =beta (2,k1)*L1+beta (1,k2)*L2
      dens(:) =rho  (2,k1)*L1+rho  (1,k2)*L2
      Qatt(:) =Q    (2,k1)*L1+Q    (1,k2)*L2
   end if
end if

call custom_media_perturb(x,y,z,Vp,Vs,dens,Qatt)

end subroutine custom_media_ak135

subroutine custom_media_perturb(x,y,z,Vp,Vs,dens,Qatt)
real,intent(in) :: x,y,z
real,dimension(2),intent(out) :: Vp,Vs,dens,Qatt

real x0,y0,z0,dx,dy,dz
real fct

fct=0.05
x0=400e3; y0=200e3; z0=-750;
dx=10e3; dy=10e3; dz=760;

if ( (x-x0)**2/dx**2 + (y-y0)**2/dy**2 <=1.0 .and. (z-z0)**2<=dz**2 ) then
   Vs=Vs*(1.0+fct)
end if
end subroutine custom_media_perturb

end module custom_mod

! vim:ft=fortran:ts=4:sw=4:nu:et:ai:
