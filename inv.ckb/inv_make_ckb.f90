program inv_make_ckb
! modified from synthetic_ckb.f in Ray-based Finite-frequency Tomo Package.

implicit none
real,parameter :: PI=3.1415296
integer :: mx,my,mz,mxc,myc,mzc
integer :: nx1,ny1,nz1,nx2,ny2,nz2
integer :: ixb,ixe,jyb,jye,kzb,kze
integer :: ci1,ci2,cj1,cj2,ck1,ck2
integer :: wid1,wid2,wid3,hlf1,hlf2,hlf3,len1,len2,len3
integer :: edgx1,edgx2,edgy1,edgy2,edgz1,edgz2
integer :: i,j,k,i1,j1,k1,fid
integer :: itype
real :: dvmax,dvmin1,dvmin2,dvmin3
real :: tmpx,tmpy,tmpz
real,allocatable,dimension(:,:,:) :: tmpckb

write(*,*) '  enter type of checkerboard pattern'
write(*,*) '   1) linear square'
write(*,*) '   2) linear cylinder'
write(*,*) '   3) linear spher'
write(*,*) '   4) cos^2 cylinder'
write(*,*) '   5) cos^2 spher'
read(*,*)    itype

write(*,*)'  enter dimension size of the inversion blocks'
read(*,*)    mx,my,mz

write(*,*)'    enter center (mxc,myc,mzc) of block, which is negative'
read(*,*) mxc,myc,mzc

write(*,*)'   maximum perturbation'
read(*,*) dvmax

write(*,*)'   x: size, edge of checkerborad and dvmin'
read(*,*) wid1, edgx1,edgx2, dvmin1

write(*,*)'   y: size, edge of checkerborad and dvmin'
read(*,*) wid2, edgy1,edgy2,dvmin2

write(*,*)'   z: size, edge of checkerborad and dvmin'
read(*,*) wid3, edgz1,edgz2,dvmin3

write(*,*)'    index of first block of ckb'
read(*,*) nx1, ny1, nz1
write(*,*)'   index of last block of ckb'
read(*,*) nx2, ny2, nz2

allocate(tmpckb(mx,my,mz)); tmpckb=0.0

! -- parameters -- 
hlf1=int((wid1+1)/2); hlf2=int((wid2+1)/2); hlf3=int((wid3+1)/2)
len1=wid1+edgx1+edgx2;len2=wid2+edgy1+edgy2;len3=wid3+edgz1+edgz2;

ixb = mxc-hlf1+1; ixe = ixb+wid1-1
jyb = myc-hlf2+1; jye = jyb+wid2-1
kzb = mzc-hlf3+1; kze = kzb+wid3-1

ck1 = kzb-edgz1; ck2 = kze+edgz2
cj1 = jyb-edgy1; cj2 = jye+edgy2
ci1 = ixb-edgx1; ci2 = ixe+edgx2

! center checkerboard
select case (itype)
case (1)
    do k=kzb,kze
       k1=min(k-kzb+1,kze-k+1)
       tmpz=( (hlf3-k1)*dvmin3+k1*dvmax ) /hlf3
       do j=jyb,jye
          j1=min(j-jyb+1,jye-j+1)
          tmpy=( (hlf2-j1)*dvmin2+j1*dvmax ) /hlf2
          do i=ixb,ixe
             i1=min(i-ixb+1,ixe-i+1)
             tmpx=( (hlf1-i1)*dvmin1+i1*dvmax ) /hlf1
             tmpckb(i,j,k)=min(min(tmpx,tmpy),tmpz)
          end do
       end do
    end do
case (2)
    do k=kzb,kze
       k1=min(k-kzb+1,kze-k+1)
       tmpz=( (hlf3-k1)*dvmin3+k1*dvmax ) /hlf3
       do j=jyb,jye
          j1=hlf2-min(j-jyb+1,jye-j+1)
          do i=ixb,ixe
             i1=hlf1-min(i-ixb+1,ixe-i+1)
             tmpy=sqrt(i1**2.0+j1**2.0)
             !if (tmpy>hlf1) then
             !   tmpx=0.0
             !else
                tmpx=( (hlf1-tmpy)*dvmax+tmpy*dvmin1 ) /hlf1
             !end if
             tmpckb(i,j,k)=min(tmpx,tmpz)
          end do
       end do
    end do
case (3)
    do k=kzb,kze
       k1=hlf3-min(k-kzb+1,kze-k+1)
       do j=jyb,jye
          j1=hlf2-min(j-jyb+1,jye-j+1)
          do i=ixb,ixe
             i1=hlf1-min(i-ixb+1,ixe-i+1)
             tmpy=sqrt(i1**2.0+j1**2.0+k1**2.0)
             tmpx=( (hlf1-i1)*dvmax+i1*dvmin1 ) /hlf1
             tmpckb(i,j,k)=tmpx
          end do
       end do
    end do
case (4)
    do k=kzb,kze
       k1=min(k-kzb+1,kze-k+1)
       tmpz=( (hlf3-k1)*dvmin3+k1*dvmax ) /hlf3
       do j=jyb,jye
          j1=hlf2-min(j-jyb+1,jye-j+1)
          tmpy=cos(PI*j1/real(hlf2)/2.0)**2
          do i=ixb,ixe
             i1=hlf1-min(i-ixb+1,ixe-i+1)
             tmpx=cos(PI*i1/real(hlf1)/2.0)**2
             tmpckb(i,j,k)=min(tmpx*tmpy*dvmax,tmpz)
          end do
       end do
    end do
case (5)
    do k=kzb,kze
       k1=hlf3-min(k-kzb+1,kze-k+1)
       tmpz=cos(PI*k1/real(hlf3)/2.0)**2
       do j=jyb,jye
          j1=hlf2-min(j-jyb+1,jye-j+1)
          tmpy=cos(PI*j1/real(hlf2)/2.0)**2
          do i=ixb,ixe
             i1=hlf1-min(i-ixb+1,ixe-i+1)
             tmpx=cos(PI*i1/real(hlf1)/2.0)**2
             tmpckb(i,j,k)=tmpx*tmpy*tmpz*dvmax
          end do
       end do
    end do
case default
    call error_except("itype error")
end select

! extend along x
do k=kzb,kze
   do j=jyb,jye
      do i=ci1-1,nx1,-1
         tmpckb(i,j,k)=-tmpckb(i+len1,j,k)
      end do
      do i=ci2+1,nx2
         tmpckb(i,j,k)=-tmpckb(i-len1,j,k)
      end do
   end do
end do
    
! extend along y
do k=kzb,kze
   do i=nx1,nx2
      do j=cj1-1,ny1,-1
         tmpckb(i,j,k)=-tmpckb(i,j+len2,k)
      end do
      do j=cj2+1,ny2
         tmpckb(i,j,k)=-tmpckb(i,j-len2,k)
      end do
   end do
end do

! extend along z
do j=ny1,ny2
   do i=nx1,nx2
      do k=ck1-1,nz1,-1
         tmpckb(i,j,k)=-tmpckb(i,j,k+len3)
      end do
      do k=ck2+1,nz2
         tmpckb(i,j,k)=-tmpckb(i,j,k-len3)
      end do
   end do
end do

! output
fid=5001
open(fid,file="block_ckb_info.dat")
write(fid,*)'  enter dimension size of the inversion blocks'
write(fid,*)    mx,my,mz

write(fid,*)'    enter center (mxc,myc,mzc) of block, which is negative'
write(fid,*) mxc,myc,mzc

write(fid,*)'   maximum perturbation'
write(fid,*) dvmax

write(*,*)'   x: size, edge of checkerborad and dvmin'
write(*,*) wid1, edgx1,edgx2, dvmin1

write(*,*)'   y: size, edge of checkerborad and dvmin'
write(*,*) wid2, edgy1,edgy2,dvmin2

write(*,*)'   z: size, edge of checkerborad and dvmin'
write(*,*) wid3, edgz1,edgz2,dvmin3

write(*,*)'   length of checkerborad'
write(*,*) len1,len2,len3

write(fid,*)'    index of first block of ckb'
write(fid,*) nx1, ny1, nz1
write(fid,*)'   index of last block of ckb'
write(fid,*) nx2, ny2, nz2
close(fid)

open(fid,file="block_ckb.dat")
do k=1,mz
do j=1,my
do i=1,mx
   write(fid,*) tmpckb(i,j,k)
end do
end do
end do
close(fid)

deallocate(tmpckb)

!----------------------------------------------------------
contains
!----------------------------------------------------------

subroutine error_except(msg)
  character (len=*),intent(in) :: msg
  print *, trim(msg)
  stop 1
end subroutine error_except

end program inv_make_ckb
