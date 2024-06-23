!pvch.f90
!pvwin用可視化，分子位置ファイル，色ファイルの出力
!必要ファイル：position.dat，syuuki.dat
!pvwin用のデータの形式については，pvwinのＨＰにある
!使い方を参照ください．

double precision :: posx,posy,posz,nkoss
real  ::  pomx,pomy,pomz
!dimension :: mcolor(nkosu)

moltype = 1
!manual change
ndat= 200
ntime0 = 0
ndt = 1
ns = 5
no = 14
ng = 2
ng2 = 10

open(1,file='posit.dat')
open(2,file='pos.dat')
open(3,file='mask.dat')
open(4,file='syuuki.dat')

read(4,*)nkoss
read(4,*)vlx2
read(4,*)vly2
read(4,*)vlz2
write(2,'(3I7)') moltype,nmol,ndat
write(2,'(3F15.5)') vlx2,vly2,vlz2
write(2,'(2I7)') ntime0, ndt


nmol = nkoss
do k22 = 1,ndat
  do j = 1, nkoss
      read(1,'(I6,3D15.7)')num1, ponx, pony, ponz
          pomx = ponx
          pomy = pony
          pomz = ponz
      write(2,'(3E15.7)') pomx, pomy, pomz
  end do
end do
do i = 1,ndat
  do i75 = 1,nkoss
     write(3,'(I7)')ns
  end do
end do