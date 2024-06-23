!-----------------------------------------------------------------------
! example.f90
! ナノ粒子の初期状態を分子動力学法で作成する．
! 必要ファイル: random0.dat, random1.dat
!-----------------------------------------------------------------------

! パラメータ（もともとparate.datにあったもの）
module parameters
    implicit none
    integer, parameter :: typmol = 3    ! 分子の種類数
    !integer, parameter :: comp = typmol*(typmol-1)/2    ! 相互作用の組み合わせ数 = typmolC2
    integer, parameter :: nummol(typmol) = (/150, 180, 150/)     ! 各分子の数
    !integer, parameter :: nkoss = 150 ! 全分子数
    ! 周期長さの単位がよくわからない → [Å]?    3.92*3.92*8.70[nm^3]にしたい
    double precision, parameter :: stdist(typmol) = (/3.92d0, 5.0d0, 3.92d0/) ! 格子定数
    double precision, parameter :: xsyul0 = stdist(1) * 5 ! x方向の周期境界長さ(無次元)
    double precision, parameter :: ysyul0 = stdist(1) * 5 ! y方向の周期境界長さ(無次元）
    double precision, parameter :: zsyul0 = 60.000d0 ! z方向の周期境界長さ(無次元）
    ! アルゴンの分子の質量39.95[kg/mol]
    double precision, parameter :: bunsi3(typmol) = (/195.084d-3, 39.950d-3, 195.084d-3/) ! 分子の質量  kg/mol
    ! 分子のLennard-Jonesパラメータσ(無次元) Pt-Pt σ=2.54[Å], Ar-Ar σ=3.40[Å], Ar-Pt σ=2.97[Å]
    double precision, parameter :: sig(typmol, typmol)  = reshape([2.54d0, 2.97d0, 2.54d0, &
                                                                   2.97d0, 3.40d0, 2.97d0, &
                                                                   2.54d0, 2.97d0, 2.54d0], [3, 3])
    ! 分子のLennard-Jonesパラメータε(無次元) Pt-Pt ε=109.2[10^(-21)J], Ar-Ar ε=1.67[10^(-21)J], Ar-Pt ε=13.50[10^(-21)J]
    double precision, parameter :: eps(typmol, typmol)  = reshape([109.2d-5, 13.50d-5, 109.2d-5, &
                                                                   13.50d-5, 1.670d-5, 13.50d-5, &
                                                                   109.2d-5, 13.50d-5, 109.2d-5], [3, 3])
    ! 接触角
    double precision, parameter :: cang(typmol, typmol) = reshape([1.00d0, 0.05d0, 1.00d0, &
                                                                   0.05d0, 1.00d0, 0.05d0, &
                                                                   1.00d0, 0.05d0, 1.00d0], [3, 3])
    double precision, parameter :: atemp1 = 150d0 ! 系内（目標）温度  K
    double precision, parameter :: cutoff33 = 3.300d0 ! カットオフ長さ/σ
    integer, parameter :: maxstp = 5000 ! 最大ステップ数
    double precision, parameter :: avoga = 6.022d+23 ! アボガドロ数
    double precision, parameter :: boltz = 1.3806662d-23 ! ボルツマン定数
    double precision, parameter :: pi = 3.141592654d0 ! 円周率
    double precision, parameter :: dt = 5.00d0 ! 無次元時間ステップ(無次元，有次元の場合fs)
end module parameters

! 変数
module variable
    use parameters
    implicit none
    double precision :: zmass(typmol) = 0.0d0 ! 分子の質量
    double precision :: cforce(typmol, typmol) = 0.0d0 ! 力の計算の係数
    ! double precision :: sig, eps ! 分子のLennard-Jonesパラメータ，σ，ε
    double precision, dimension(typmol, typmol) :: xcutof, ycutof, zcutof ! ポテンシャルのカットオフ長さx,y,z方向
    integer :: nowstp ! 現在のステップ数
    double precision :: xsyul, ysyul, zsyul ! ポテンシャルのカットオフ長さx,y,z方向，x,y,z方向の周期長さ
end module variable

module molecules_module
    use parameters
    implicit none

    !構造体の定義
    type :: mol_type
        double precision :: posx, posy, posz
        double precision :: velx, vely, velz
        double precision :: forx, fory, forz  
        double precision :: poten, ukine 
    end type mol_type

    type mol_array
        type(mol_type), allocatable :: array(:)
    end type mol_array

    type(mol_array), dimension(typmol) :: mol  ! 上Pt, 中Ar, 下Pt

end module molecules_module

module pvch
    use parameters
    implicit none
    integer, parameter :: tlnkoss = nummol(1) + nummol(2) + nummol(3)
    integer, parameter :: moltype = 1
    integer, parameter :: ndat = 200
    integer, parameter :: ntime0 = 0
    integer, parameter :: ndt = 1
end module pvch

program main
    use variable, only: nowstp
    use parameters
    use molecules_module
    use pvch
    implicit none
    integer :: i, j

    ! 配列初期化
    allocate(mol(1)%array(nummol(1)))
    allocate(mol(2)%array(nummol(2)))
    allocate(mol(3)%array(nummol(3)))

    ! 読み込み用乱数ファイル
        open(1,file='random1.dat', status='old')
        open(2,file='random0.dat', status='old')
    ! 各分子の位置データの出力
        open(10,file='posit_PtUp.dat')
        open(11,file='posit_Ar.dat')
        open(12,file='posit_PtDw.dat')
    ! 可視化用のpvch.fを移植 
        open(15,file='pos.dat')
    ! 各分子の速度データの出力
        open(20,file='veloc_PtUp.dat')
        open(21,file='veloc_Ar.dat')
        open(22,file='veloc_PtDw.dat')
    ! 系のエネルギーデータの出力
        open(30,file='energy_PtUp.dat')
        open(31,file='energy_Ar.dat')
        open(32,file='energy_PtDw.dat')
    ! 系の温度データの出力
        open(40,file='tempe.dat')
    ! 系の周期長さの出力
        open(50,file='syuuki.dat')



    ! 各分子の最終位置データの出力
        open(80,file='finpos.dat')
    !　分子の色
        open(90,file='mask.dat')

    write(15,'(3I7)') moltype, tlnkoss, ndat
    do i = 1,ndat
        do j = 1, nummol(1)
            write(90,'(I7)') 7      ! 黄色
        end do
        do j = 1, nummol(2)
            write(90,'(I7)') 14     ! 赤色
        end do
        do j = 1, nummol(3)
            write(90,'(I7)') 7      ! 黄色
        end do
    end do
    

    nowstp = 0

    call seting ! 各分子の初期位置，初期速度などの設定
    !!call cortra ! 系内の全分子の並進速度の補正
    !!call jyusin ! 系内の全分子の重心の補正
    !call scale ! 系内の全分子の温度の補正
    !call record ! データの出力１
    !call record2 ! データの出力２

    do i = 1, maxstp
        nowstp = i
        ! ステップ数が500の倍数のとき
        if (mod(nowstp,500) == 0) then
            write(6,*) nowstp
        endif
        ! ステップ数が100の倍数のとき
        if (nowstp >= maxstp/5 .and. mod(nowstp,100) == 0) then
        	!!call cortra	! 系内の全分子の並進速度の補正
            !!call jyusin	! 系内の全分子の重心の補正
        	!call scale ! 系内の全分子の温度の補正
        endif

        call calcu ! 各分子に働く力，速度，位置の分子動力学計算
        call bound ! 境界条件の付与
        
        ! ステップ数が100の倍数+1のとき
        if(mod(nowstp, 100) == 1) then
          call record ! データの出力１
          call record2 ! データの出力２
        endif
    end do
    call record3 ! データの出力３

contains
    subroutine seting ! 各分子の初期位置，初期速度などの設定
        use variable
        use parameters
        use molecules_module
        implicit none
        integer :: num, i, j, k
        double precision :: x = 0.0000d0, y = 0.0000d0, z = 0.0000d0
        double precision :: ofstx1, ofsty1, ofstz1
        double precision :: ran, alpha, beta, cr
        double precision :: vx, vy, vz

        do i = 1, typmol
            zmass(i) = 1.00d+26*bunsi3(i)/avoga
            do j = 1, typmol
                cforce(i, j) = cang(i, j)*24.00d0*eps(i, j)/sig(i, j)
            end do 
        end do
        
        xsyul = xsyul0
        ysyul = ysyul0
        zsyul = zsyul0
        write(50,*)xsyul
        write(50,*)ysyul
        write(50,*)zsyul
        write(15,*)xsyul0, ysyul0, zsyul0
        write(15,*)ntime0, ndt

        num = 0

        do i = 1, typmol
            do j = 1, typmol
                xcutof(i, j) = xsyul - cutoff33*sig(i, j)
                ycutof(i, j) = ysyul - cutoff33*sig(i, j)
                zcutof(i, j) = zsyul - cutoff33*sig(i, j)
            end do
        end do

        !上段のPt配置
        ofstx1 = stdist(1)*0.25d0
        ofsty1 = stdist(1)*0.25d0
        ofstz1 = zsyul0 - stdist(1)
        do k = 1,3
            z = ofstz1 + dble(k-1)*stdist(1)*0.5d0
            do i = 1,10
                x = ofstx1 + dble(i-1)*stdist(1)*0.5d0
                do j = 1,5
                    if(mod(k,2) == 0) then  !z偶数
                        if(mod(i,2) == 0) then
                            y = ofsty1 + dble(j-1)*stdist(1)   !x偶数
                        else
                            y = ofsty1 + dble(j-1)*stdist(1) + stdist(1)*0.5d0    !x奇数
                        endif
                    else    !z奇数
                        if(mod(i,2) == 0) then
                            y = ofsty1 + dble(j-1)*stdist(1) + stdist(1)*0.5d0    !x偶数
                        else
                            y = ofsty1 + dble(j-1)*stdist(1)   !x奇数
                        endif
                    endif
                    num = num + 1
                    mol(1)%array(num)%posx = x
                    mol(1)%array(num)%posy = y
                    mol(1)%array(num)%posz = z
                end do
            end do
        end do

        !中段のAr配置
        num = 0
        ofstx1 = ysyul0*0.50d0 - stdist(2)*1.25d0
        ofsty1 = ysyul0*0.50d0 - stdist(2)*1.25d0
        ofstz1 = zsyul0*0.50d0 - stdist(2)*2.0d0
        do k = 1,10
            z = ofstz1 + dble(k-1)*stdist(2)*0.5d0
            do i = 1,6
                x = ofstx1 + dble(i-1)*stdist(2)*0.5d0
                do j = 1,3
                    if(mod(k,2) == 0) then  !z偶数
                        if(mod(i,2) == 0) then
                            y = ofsty1 + dble(j-1)*stdist(2)   !x偶数
                        else
                            y = ofsty1 + dble(j-1)*stdist(2) + stdist(2)*0.5d0    !x奇数
                        endif
                    else    !z奇数
                        if(mod(i,2) == 0) then
                            y = ofsty1 + dble(j-1)*stdist(2) + stdist(2)*0.5d0    !x偶数
                        else
                            y = ofsty1 + dble(j-1)*stdist(2)   !x奇数
                        endif
                    endif
                    num = num + 1
                    mol(2)%array(num)%posx = x
                    mol(2)%array(num)%posy = y
                    mol(2)%array(num)%posz = z
                end do
            end do
        end do

        !下段のPt配置
        num = 0
        ofstx1 = stdist(3)*0.25d0
        ofsty1 = stdist(3)*0.25d0
        ofstz1 = 0.0000d0
        do k = 1,3
            z = ofstz1 + dble(k-1)*stdist(3)*0.5d0
            do i = 1,10
                x = ofstx1 + dble(i-1)*stdist(3)*0.5d0
                do j = 1,5
                    if(mod(k,2) == 0) then  !z偶数
                        if(mod(i,2) == 0) then
                            y = ofsty1 + dble(j-1)*stdist(3)   !x偶数
                        else
                            y = ofsty1 + dble(j-1)*stdist(3) + stdist(3)*0.5d0    !x奇数
                        endif
                    else    !z奇数
                        if(mod(i,2) == 0) then
                            y = ofsty1 + dble(j-1)*stdist(3) + stdist(3)*0.5d0    !x偶数
                        else
                            y = ofsty1 + dble(j-1)*stdist(3)   !x奇数
                        endif
                    endif
                    num = num + 1
                    mol(3)%array(num)%posx = x
                    mol(3)%array(num)%posy = y
                    mol(3)%array(num)%posz = z
                end do
            end do
        end do

        cr = 1.00d-6
        do i = 1, nummol(2)
            read(1,*)ran
            alpha = pi*ran
            read(2,*)ran
            beta = 2.000d0*pi*ran
            vx = dsin(alpha)*dcos(beta)*cr
            vy = dsin(alpha)*dsin(beta)*cr
            vz = dcos(alpha)*cr
            mol(2)%array(i)%velx = vx
            mol(2)%array(i)%vely = vy
            mol(2)%array(i)%velz = vz
        end do
    end subroutine seting

    subroutine scale ! 系内の全分子の温度の補正
        use variable, !only: velx, vely, velz, zmass
        use parameters
        use molecules_module
        implicit none
        double precision :: temptp, vel2, aimtem, aimnot, baiss
        integer :: i, j

        do j = 1, typmol
            temptp = 0.0d0
            do i = 1, nummol(j)
                vel2 = mol(j)%array(i)%velx*mol(j)%array(i)%velx &
                     + mol(j)%array(i)%vely*mol(j)%array(i)%vely &
                     + mol(j)%array(i)%velz*mol(j)%array(i)%velz
                temptp = temptp + vel2
            end do
            temptp = temptp / nummol(j) / 1.000d+16
            aimtem = atemp1
            aimnot = 3.0d0 * boltz * aimtem / zmass(j)
            baiss = dsqrt(aimnot / temptp)

            ! 速度ベクトルのスケーリング
            do i = 1, nummol(j)
                mol(j)%array(i)%velx = mol(j)%array(i)%velx * baiss
                mol(j)%array(i)%vely = mol(j)%array(i)%vely * baiss
                mol(j)%array(i)%velz = mol(j)%array(i)%velz * baiss
            end do
        end do
    end subroutine scale

    subroutine calcu ! 各分子に働く力，速度，位置の分子動力学計算
        use variable
        use parameters
        use molecules_module
        implicit none
        integer :: i, i1, j1
        double precision :: divx, divy, divz, dist
        double precision :: dit2, dit4, dit6, dit8, dit12, dit14
        double precision :: ppp, force, forcex, forcey, forcez
        double precision :: vxene, vyene, vzene, vene

        do j = 1, typmol
            do i = 1, nummol(j)
                mol(j)%array(i)%forx  = 0.0000d0
                mol(j)%array(i)%fory  = 0.0000d0
                mol(j)%array(i)%forz  = 0.0000d0
                mol(j)%array(i)%poten = 0.0000d0
                mol(j)%array(i)%ukine = 0.0000d0
            end do
        end do

        ! 分子間の相互作用力 → ポテンシャルエネルギー
        ! 同じ分子同士の影響
        do i = 1, typmol
            do j = j, typmol
                do i1 = 1, nummol(i)
                    do j1 = 1, nummol(j)
                        if(i == j .and. i1 == j1) then
                            cycle
                        else
                            divx = mol(i)%array(i1)%posx - mol(j)%array(j1)%posx
                            divy = mol(i)%array(i1)%posy - mol(j)%array(j1)%posy
                            divz = mol(i)%array(i1)%posz - mol(j)%array(j1)%posz
            
                            if (divx < -xcutof(i, j)) then
                                divx = divx + xsyul
                            else if(divx > xcutof(i, j)) then
                                divx = divx - xsyul
                            endif
            
                            if (divy < -ycutof(i, j)) then
                                divy = divy + ysyul
                            else if(divy > ycutof(i, j)) then
                                divy = divy - ysyul
                            endif
            
                            if (divz < -zcutof(i, j)) then
                                divz = divz + zsyul
                            else if(divz > zcutof(i, j)) then
                                divz = divz - zsyul
                            endif
            
                            divx = divx / sig(i, j)
                            divy = divy / sig(i, j)
                            divz = divz / sig(i, j)

                            if (divx < -cutoff33 .or. divx > cutoff33) then
                                cycle
                            endif
            
                            if (divy < -cutoff33 .or. divy > cutoff33) then
                                cycle
                            endif
            
                            if (divz < -cutoff33 .or. divz > cutoff33) then
                                cycle
                            endif

                            dit2 = divx*divx + divy*divy + divz*divz
                            dist = dsqrt(dit2)
            
                            if(dist > cutoff33) then
                                cycle
                            endif
            
                            dit4   = dit2*dit2
                            dit6   = dit4*dit2
                            dit8   = dit4*dit4
                            dit12  = dit6*dit6
                            dit14  = dit8*dit6
                            ppp    = 4.00d0*eps(i, j)*(1.00d0/dit12-1.00d0/dit6)
                            force  = cforce(i, j)*(-2.00d0/dit14+1.00d0/dit8)
                            forcex = -force*divx/zmass(i)
                            mol(i)%array(i1)%forx = mol(i)%array(i1)%forx + forcex
                            forcex = -force*divx/zmass(j)
                            mol(j)%array(j1)%forx = mol(j)%array(j1)%forx - forcex
                            forcey = -force*divy/zmass(i)
                            mol(i)%array(i1)%fory = mol(i)%array(i1)%fory + forcey
                            forcey = -force*divy/zmass(j)
                            mol(j)%array(j1)%fory = mol(j)%array(j1)%fory - forcey
                            forcez = -force*divz/zmass(i)
                            mol(i)%array(i1)%forz = mol(i)%array(i1)%forz + forcez
                            forcez = -force*divz/zmass(j)
                            mol(j)%array(j1)%forz = mol(j)%array(j1)%forz - forcez
                            mol(i)%array(i1)%poten = mol(i)%array(i1)%poten + ppp*0.500d0
                            mol(j)%array(j1)%poten = mol(j)%array(j1)%poten + ppp*0.500d0
                        end if
                    end do
                end do
            end do
        end do

        ! ! 異なる分子同士の影響
        ! do i = 1, typmol
        !     if (i /= 2) then     ! Ar-Ptの処理    配列は1が上Pt, 2がAr, 3がしたPt
        !         do i1 = 1, nummol(i)       ! Pt
        !             do i2 = 1, nummol(2)    ! Ar
        !                 divx = mol(i)%array(i1)%posx - mol(2)%array(i2)%posx
        !                 divy = mol(i)%array(i1)%posy - mol(2)%array(i2)%posy
        !                 divz = mol(i)%array(i1)%posz - mol(2)%array(i2)%posz

        !                 if (divx < -xcutof) then             ! cutof, syulは1がPt-Pt, 2がAr-Ar, 3がPt-Ar
        !                     divx = divx + xsyul
        !                 else if(divx > xcutof) then
        !                     divx = divx - xsyul
        !                 endif

        !                 if (divy < -ycutof) then
        !                     divy = divy + ysyul
        !                 else if(divy > ycutof) then
        !                     divy = divy - ysyul
        !                 endif

        !                 if (divz < -zcutof) then
        !                     divz = divz + zsyul
        !                 else if(divz > zcutof) then
        !                     divz = divz - zsyul
        !                 endif

        !                 divx = divx / sig33(4)
        !                 divy = divy / sig33(4)
        !                 divz = divz / sig33(4)

        !                 if (divx < -cutoff33 .or. divx > cutoff33) then
        !                     cycle
        !                 endif

        !                 if (divy < -cutoff33 .or. divy > cutoff33) then
        !                     cycle
        !                 endif

        !                 if (divz < -cutoff33 .or. divz > cutoff33) then
        !                     cycle
        !                 endif

        !                 dit2 = divx*divx + divy*divy + divz*divz
        !                 dist = dsqrt(dit2)

        !                 if(dist > cutoff33) then
        !                     cycle
        !                 endif

        !                 dit4   = dit2*dit2
        !                 dit6   = dit4*dit2
        !                 dit8   = dit4*dit4
        !                 dit12  = dit6*dit6
        !                 dit14  = dit8*dit6
        !                 ppp    = 4.00d0*eps33(4)*(1.00d0/dit12-1.00d0/dit6)
        !                 force  = cforce(4)*(-2.00d0/dit14+1.00d0/dit8)   ! dit13, dit7じゃない？
        !                 forcex = -force*divx/zmass(i)
        !                 mol(i)%array(i1)%forx = mol(i)%array(i1)%forx + forcex
        !                 forcex = -force*divx/zmass(2)
        !                 mol(2)%array(i2)%forx = mol(2)%array(i2)%forx - forcex
        !                 forcey = -force*divy/zmass(i)
        !                 mol(i)%array(i1)%fory = mol(i)%array(i1)%fory + forcey
        !                 forcey = -force*divy/zmass(2)
        !                 mol(2)%array(i2)%fory = mol(2)%array(i2)%fory - forcey
        !                 forcez = -force*divz/zmass(i)
        !                 mol(i)%array(i1)%forz = mol(i)%array(i1)%forz + forcez
        !                 forcez = -force*divz/zmass(2)
        !                 mol(2)%array(i2)%forz = mol(2)%array(i2)%forz - forcez
        !                 mol(i)%array(i1)%poten = mol(i)%array(i1)%poten + ppp*0.500d0
        !                 mol(2)%array(i2)%poten = mol(2)%array(i2)%poten + ppp*0.500d0
        !             end do
        !         end do
        !     end if    
        ! end do
        
        do j = 1, typmol
            do i = 1, nummol(j)
                vxene = mol(j)%array(i)%velx + mol(j)%array(i)%forx*0.500d0*dt
                vyene = mol(j)%array(i)%vely + mol(j)%array(i)%fory*0.500d0*dt
                vzene = mol(j)%array(i)%velz + mol(j)%array(i)%forz*0.500d0*dt
                vene = vxene*vxene + vyene*vyene + vzene*vzene
                mol(j)%array(i)%ukine = 0.500d0*zmass(j)*vene
            end do

            do i = 1, nummol(j)
                mol(j)%array(i)%velx = mol(j)%array(i)%velx + mol(j)%array(i)%forx*dt
                mol(j)%array(i)%vely = mol(j)%array(i)%vely + mol(j)%array(i)%fory*dt
                mol(j)%array(i)%velz = mol(j)%array(i)%velz + mol(j)%array(i)%forz*dt
                mol(j)%array(i)%posx = mol(j)%array(i)%posx + mol(j)%array(i)%velx*dt
                mol(j)%array(i)%posy = mol(j)%array(i)%posy + mol(j)%array(i)%vely*dt
                mol(j)%array(i)%posz = mol(j)%array(i)%posz + mol(j)%array(i)%velz*dt
            end do
        end do

        do i = 1, 150
            mol(1)%array(i)%velx = 0.0000d0
            mol(1)%array(i)%vely = 0.0000d0
            mol(1)%array(i)%velz = 0.0000d0
            mol(1)%array(i)%forx = 0.0000d0
            mol(1)%array(i)%fory = 0.0000d0
            mol(1)%array(i)%forz = 0.0000d0
        end do

        do i = 1, 150
            mol(3)%array(i)%velx = 0.0000d0
            mol(3)%array(i)%vely = 0.0000d0
            mol(3)%array(i)%velz = 0.0000d0
            mol(3)%array(i)%forx = 0.0000d0
            mol(3)%array(i)%fory = 0.0000d0
            mol(3)%array(i)%forz = 0.0000d0
        end do

    end subroutine calcu

    subroutine bound ! 境界条件の付与
        use variable, only: xsyul, ysyul, zsyul
        use parameters
        use molecules_module
        implicit none
        integer :: i, j
        do j = 1, typmol
            do i = 1, nummol(j)
                if(mol(j)%array(i)%posx < 0.00d0) then
                    mol(j)%array(i)%posx = mol(j)%array(i)%posx + xsyul
                else if(mol(j)%array(i)%posx > xsyul) then
                    mol(j)%array(i)%posx = mol(j)%array(i)%posx - xsyul
                endif

                if(mol(j)%array(i)%posy < 0.00d0) then
                    mol(j)%array(i)%posy = mol(j)%array(i)%posy + ysyul
                else if(mol(j)%array(i)%posy > ysyul) then
                    mol(j)%array(i)%posy = mol(j)%array(i)%posy - ysyul
                endif

                ! z方向は周期境界の補正を行わない
                ! if(mol(j)%array(i)%posz < 0.00d0) then
                !     mol(j)%array(i)%posz = mol(j)%array(i)%posz + zsyul
                ! else if(mol(j)%array(i)%posz > zsyul) then
                !     mol(j)%array(i)%posz = mol(j)%array(i)%posz - zsyul
                ! endif
            end do
        end do
    end subroutine bound

    subroutine record ! データの出力１
        use variable !only: posx, posy, posz, velx, vely, velz
        use parameters
        use molecules_module
        implicit none
        integer :: i, j

        do i = 1, nummol(1)
            write(10, '(I6, 3E15.7)') i, mol(1)%array(i)%posx, mol(1)%array(i)%posy, mol(1)%array(i)%posz
            write(20, '(I6, 3E15.7)') i, mol(1)%array(i)%velx, mol(1)%array(i)%vely, mol(1)%array(i)%velz
        end do
        do i = 1, nummol(2)
            write(11, '(I6, 3E15.7)') i, mol(2)%array(i)%posx, mol(2)%array(i)%posy, mol(2)%array(i)%posz
            write(21, '(I6, 3E15.7)') i, mol(2)%array(i)%velx, mol(2)%array(i)%vely, mol(2)%array(i)%velz
        end do
        do i = 1, nummol(3)
            write(12, '(I6, 3E15.7)') i, mol(3)%array(i)%posx, mol(3)%array(i)%posy, mol(3)%array(i)%posz
            write(22, '(I6, 3E15.7)') i, mol(3)%array(i)%velx, mol(3)%array(i)%vely, mol(3)%array(i)%velz
        end do

        do j = 1, typmol
            do i = 1, nummol(j)
                write(15, '(3E15.7)') mol(j)%array(i)%posx, mol(j)%array(i)%posy, mol(j)%array(i)%posz
            end do
        end do
    end subroutine record

    subroutine record2 ! データの出力２
        use variable !only: poten, ukine
        use parameters
        use molecules_module
        implicit none
        double precision, dimension(typmol) :: totene, totpot, totkin, temp
        integer :: i, j

        ! エネルギーの合計計算
        do j = 1, typmol
            do i = 1, nummol(j)
                totpot = totpot + mol(j)%array(i)%poten
                totkin = totkin + mol(j)%array(i)%ukine
            end do

            ! エネルギーを大きな数で割る処理（正規化や単位変換のため）
            totpot(j) = totpot(j) / 1.00d16
            totkin(j) = totkin(j) / 1.00d16
            totene(j) = totpot(j) + totkin(j)

            ! 温度計算
            temp(j) = 2.0d0 * totkin(j) / (3.0d0 * dble(nummol(j)) * boltz)
        end do

        write(30, '(I6, 4E15.7)') (nowstp+99)*int(dt), totene(1), totpot(1), totkin(1)
        write(31, '(I6, 4E15.7)') (nowstp+99)*int(dt), totene(2), totpot(2), totkin(2)
        write(32, '(I6, 4E15.7)') (nowstp+99)*int(dt), totene(3), totpot(3), totkin(3)   
        write(40, '(4E15.7)') temp(1), temp(2), temp(3)
    end subroutine record2

    subroutine record3
        use variable !only: posx, posy, posz, velx, vely, velz
        use parameters
        use molecules_module
        implicit none
        integer :: i, j
        do j = 1, typmol
            do i = 1, nummol(j)
                write(50, '(I6, 6E15.7)') & 
                i, mol(j)%array(i)%posx, mol(j)%array(i)%posy, mol(j)%array(i)%posz, &
                   mol(j)%array(i)%velx, mol(j)%array(i)%vely, mol(j)%array(i)%velz
            end do
        end do
    end subroutine record3
end program main