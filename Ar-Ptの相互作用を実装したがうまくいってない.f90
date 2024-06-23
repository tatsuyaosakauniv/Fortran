!-----------------------------------------------------------------------
! example.f90
! ナノ粒子の初期状態を分子動力学法で作成する．
! 必要ファイル: random0.dat, random1.dat
!-----------------------------------------------------------------------

! パラメータ（もともとparate.datにあったもの）
module parameters
    implicit none
    integer, parameter :: typmol = 3    !分子の種類数
    integer, parameter :: comp = 3!typmol*(typmol-1)/2    ! 相互作用の組み合わせ数 = typmolC2
    integer, parameter :: nkoss = 50 ! 1種類の分子数
    ! 周期長さの単位がよくわからない → [Å]?    3.92*3.92*8.70[nm^3]にしたい
    double precision, parameter :: xsyul0 = 39.200D0 ! x方向の周期境界長さ(無次元)
    double precision, parameter :: ysyul0 = 39.200D0 ! y方向の周期境界長さ(無次元）
    double precision, parameter :: zsyul0 = 87.000D0 ! z方向の周期境界長さ(無次元）
    ! アルゴンの分子の質量39.95[kg/mol]
    double precision, parameter :: bunsi3(typmol) = (/39.950D-3, 195.084D-3, 39.950D-3/) ! 分子の質量  kg/mol
    ! Ar-Ar σ=3.40[Å], Pt-Pt σ=2.54[Å], Ar-Pt σ=2.97[Å]
    double precision, parameter :: sig33(comp) = (/3.40D0, 2.54D0, 2.97D0/) ! 分子のLennard-Jonesパラメータσ(無次元)
    ! Ar-Ar ε=1.67[10^(-21)J], Pt-Pt ε=109.2[10^(-21)J], Ar-Pt ε=13.50[10^(-21)J]
    double precision, parameter :: eps33(comp) = (/1.67D-5, 109.2D-5, 13.50D-5/) ! 分子のLennard-Jonesパラメータε(無次元)
    double precision, parameter :: atemp1 = 150D0 ! 系内（目標）温度  K
    double precision, parameter :: cutoff33 = 3.300D0 ! カットオフ長さ/σ
    integer, parameter :: maxstep = 20000 ! 最大ステップ数
    double precision, parameter :: avoga = 6.022D+23 ! アボガドロ数
    double precision, parameter :: boltz = 1.3806662D-23 ! ボルツマン定数
    double precision, parameter :: pi = 3.141592654D0 ! 円周率
    double precision, parameter :: dt = 5.00D0 ! 無次元時間ステップ(無次元，有次元の場合fs)
end module parameters

! 変数
module variable
    use parameters, !only: typmol, comp
    implicit none
    double precision, dimension(typmol) :: zmass, cforce! = (/0.0D0, 0.0D0/) ! 分子の質量，力の計算の係数
    double precision :: sig, eps ! 分子のLennard-Jonesパラメータ，σ，ε
    double precision, dimension(typmol) :: xcutof, ycutof, zcutof ! ポテンシャルのカットオフ長さx,y,z方向
    integer :: nowstp ! 現在のステップ数
    double precision, dimension(typmol) :: xsyul, ysyul, zsyul ! ポテンシャルのカットオフ長さx,y,z方向，x,y,z方向の周期長さ
end module variable

module molecules_module
    use parameters, !only: nkoss, typmol
    use variable
    implicit none

    !構造体の定義
    type :: mol_type
        double precision :: posx, posy, posz
        double precision :: velx, vely, velz
        double precision :: forx  = 0.0000D0, fory  = 0.0000D0, forz  = 0.0000D0
        double precision :: poten = 0.0000D0, ukine = 0.0000D0
    end type mol_type

    type(mol_type), dimension(nkoss, typmol) :: mol  !上Pt, 中Ar, 下Pt

end module molecules_module

module pvch
    use parameters, only: nkoss, typmol
    implicit none
    integer, parameter :: tlnkoss = nkoss * typmol
    integer, parameter :: moltype = 1
    integer, parameter :: ndat = 200
    integer, parameter :: ntime0 = 0
    integer, parameter :: ndt = 1
end module pvch

program main
    use parameters
    use variable!, only: nowstp
    use molecules_module
    use pvch
    implicit none

    integer :: i, j

    ! 読み込み用乱数ファイル
        open(1,file='random1.dat', status='old')
        open(2,file='random0.dat', status='old')
        open(3,file='posit.dat')
    ! 各分子の速度データの出力
        open(4,file='veloc.dat')
    ! 系のエネルギーデータの出力
        open(7,file='energy.dat')
    ! 系の温度データの出力
    !   open(8,file='tempe.dat')
    ! 系の周期長さの出力
        open(9,file='syuuki.dat')
    ! 各分子の最終位置データの出力
       open(10,file='finpos.dat')
    !   write(9,*)nkoss
    ! 可視化用のpvch.fを移植
       open(11,file='pos.dat')
    !　分子の色
       open(12,file='mask.dat')

    write(11,'(3I7)') moltype, tlnkoss, ndat
    do i = 1,ndat
        do j = 1, nkoss
            write(12,'(I7)') 1
        end do
        do j = 1, nkoss
            write(12,'(I7)') 14
        end do
        do j = 1, nkoss
            write(12,'(I7)') 1
        end do
    end do

    nowstp = 0

    call seting ! 各分子の初期位置，初期速度などの設定
    !call cortra
    !call scale2
    !call jyusin
    call record ! データの出力１
    call record2 ! データの出力２

    do i = 1, maxstep
        nowstp = i
        ! ステップ数が500の倍数のとき
        if (mod(nowstp,500) == 0) then
            write(6,*) nowstp
        endif
        ! ステップ数が100の倍数のとき
        if (mod(nowstp,100) == 0) then
            !call cortra     ! 系内の全分子の並進速度の補正
            !call scale2     ! 系内の全分子の温度の補正
            !call jyusin     ! 系内の全分子の重心の補正
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
        use parameters
        use variable
        use molecules_module
        use pvch
        implicit none
        integer :: num, i, j, k, ios
        double precision :: x = 0.0000D0, y = 0.0000D0, z = 0.0000D0
        double precision :: ofstx1 = 0.0000D0, ofsty1 = 0.0000D0, ofstz1
        double precision :: stdist, ran, alpha, beta, cr
        double precision :: vx, vy, vz

        write(9,*)xsyul0
        write(9,*)ysyul0
        write(9,*)zsyul0
        write(11,*)xsyul0, ysyul0, zsyul0
        write(11,*)ntime0, ndt

        do i = 1, typmol
            xsyul(i) = xsyul0
            ysyul(i) = ysyul0
            zsyul(i) = zsyul0
            zmass(i) = 1.00D+26*bunsi3(i)/avoga
        end do

        do i = 1, comp
            sig   = sig33(i)
            eps   = eps33(i)
            cforce(i) = 24.00D0*eps/sig     ! ポテンシャルを偏微分して得られる相互作用力の係数
            xcutof(i) = xsyul(i) - cutoff33*sig33(i)
            ycutof(i) = ysyul(i) - cutoff33*sig33(i)
            zcutof(i) = zsyul(i) - cutoff33*sig33(i)
        end do

        !上段のPt配置
        num = 0
        ofstz1 = 76.0000D0
        stdist = 4.00D0     !原子間距離
        do k = 1,2
            z = ofstz1 + dble(k-1)*stdist*0.5D0
            do i = 1,5
                x = ofstx1 + dble(i-1)*stdist*0.5D0
                do j = 1,5
                    if(mod(k,2) == 0) then  !z偶数
                        if(mod(i,2) == 0) then
                            y = ofsty1 + dble(j-1)*stdist   !x偶数
                        else
                            y = ofsty1 + dble(j-1)*stdist + stdist*0.5D0    !x奇数
                        endif
                    else    !z奇数
                        if(mod(i,2) == 0) then
                            y = ofsty1 + dble(j-1)*stdist + stdist*0.5D0    !x偶数
                        else
                            y = ofsty1 + dble(j-1)*stdist   !x奇数
                        endif
                    endif
                    num = num + 1
                    mol(num, 1)%posx = x
                    mol(num, 1)%posy = y
                    mol(num, 1)%posz = z
                end do
            end do
        end do

        !中段のAr配置
        num = 0
        ofstz1 = 38.0000D0
        stdist = 4.00D0     !原子間距離
        do k = 1,2
            z = ofstz1 + dble(k-1)*stdist*0.5D0
            do i = 1,5
                x = ofstx1 + dble(i-1)*stdist*0.5D0
                do j = 1,5
                    if(mod(k,2) == 0) then  !z偶数
                        if(mod(i,2) == 0) then
                            y = ofsty1 + dble(j-1)*stdist   !x偶数
                        else
                            y = ofsty1 + dble(j-1)*stdist + stdist*0.5D0    !x奇数
                        endif
                    else    !z奇数
                        if(mod(i,2) == 0) then
                            y = ofsty1 + dble(j-1)*stdist + stdist*0.5D0    !x偶数
                        else
                            y = ofsty1 + dble(j-1)*stdist   !x奇数
                        endif
                    endif
                    num = num + 1
                    mol(num, 2)%posx = x
                    mol(num, 2)%posy = y
                    mol(num, 2)%posz = z
                end do
            end do
        end do

        !下段のPt配置
        num = 0
        ofstz1 = 0.0000D0
        do k = 1,2
            z = ofstz1 + dble(k-1)*stdist*0.5D0
            do i = 1,5
                x = ofstx1 + dble(i-1)*stdist*0.5D0
                do j = 1,5
                    if(mod(k,2) == 0) then  !z偶数
                        if(mod(i,2) == 0) then
                            y = ofsty1 + dble(j-1)*stdist   !x偶数
                        else
                            y = ofsty1 + dble(j-1)*stdist + stdist*0.5D0    !x奇数
                        endif
                    else    !z奇数
                        if(mod(i,2) == 0) then
                            y = ofsty1 + dble(j-1)*stdist + stdist*0.5D0    !x偶数
                        else
                            y = ofsty1 + dble(j-1)*stdist   !x奇数
                        endif
                    endif
                    num = num + 1
                    mol(num, 3)%posx = x
                    mol(num, 3)%posy = y
                    mol(num, 3)%posz = z
                end do
            end do
        end do

        ! Reset the file position again to read velocities
    !close(1)

    ! Read and set molecule velocities
    cr = 1.00D-6
    do i = 1, typmol
        do j = 1, nkoss
            if(i == 2) then
            read(1, *) ran
            !if (ios /= 0) then
            !    print *, "Error reading random1.dat"
            !    stop
            !end if
            alpha = pi * ran
            read(2, *) ran
            !if (ios /= 0) then
            !    print *, "Error reading random0.dat"
            !    stop
            !end if
            beta = 2.000D0 * pi * ran
            vx = dsin(alpha) * dcos(beta) * cr
            vy = dsin(alpha) * dsin(beta) * cr
            vz = dcos(alpha) * cr
            mol(j, i)%velx = vx
            mol(j, i)%vely = vy
            mol(j, i)%velz = vz
            else if (i /= 2) then
                mol(j, i)%velx = 0.0000D0
                mol(j, i)%vely = 0.0000D0
                mol(j, i)%velz = 0.0000D0
            end if
        end do
    end do

    ! Close the file after reading velocities
    close(1)
    close(2)
    end subroutine seting

    subroutine cortra ! 系内の全分子の並進速度の補正
        use parameters
        use variable, !only: velx, vely, velz
        use molecules_module
        implicit none
        double precision, dimension(typmol) :: trvx = 0.0D0, trvy = 0.0D0, trvz = 0.0D0
        integer :: i, j

        ! 速度ベクトルの成分の平均を計算
        do i = 1, typmol
            do j = 1, nkoss
                trvx(i) = trvx(i) + mol(j, i)%velx
                trvy(i) = trvy(i) + mol(j, i)%vely
                trvz(i) = trvz(i) + mol(j, i)%velz
            end do
        end do

        do i = 1, typmol
            trvx(i) = trvx(i) / nkoss
            trvy(i) = trvy(i) / nkoss
            trvz(i) = trvz(i) / nkoss
        end do

        ! 速度ベクトルから平均を引いて中心補正
        do i = 1, typmol
            do j = 1, nkoss
                mol(j, i)%velx = mol(j, i)%velx - trvx(i)
                mol(j, i)%vely = mol(j, i)%vely - trvy(i)
                mol(j, i)%velz = mol(j, i)%velz - trvz(i)
            end do
        end do
    end subroutine cortra

    subroutine jyusin ! 系内の全分子の重心の補正
        use parameters
        use variable, !only: posx, posy, posz
        use molecules_module
        implicit none
        double precision :: cmsx, cmsy, cmsz
        double precision, dimension(typmol) :: tcmsx = 0.0000D0, tcmsy = 0.0000D0, tcmsz = 0.0000D0
        integer :: i, j

        cmsx = xsyul0 * 0.500D0
        cmsy = ysyul0 * 0.500D0
        cmsz = zsyul0 * 0.500D0

        !座標の平均を計算
        do i = 1, typmol
            do j = 1, nkoss
                tcmsx(i) = tcmsx(i) + mol(j, i)%posx
                tcmsy(i) = tcmsy(i) + mol(j, i)%posy
                tcmsz(i) = tcmsz(i) + mol(j, i)%posz
            end do
        end do

        do i = 1, typmol
            tcmsx(i) = cmsx - tcmsx(i)/dble(nkoss)
            tcmsy(i) = cmsy - tcmsy(i)/dble(nkoss)
            tcmsz(i) = cmsz - tcmsz(i)/dble(nkoss)
        end do
        
        do i = 1, typmol
            do j = 1, nkoss
                mol(j, i)%posx = mol(j, i)%posx + tcmsx(i)
                mol(j, i)%posy = mol(j, i)%posy + tcmsy(i)
                mol(j, i)%posz = mol(j, i)%posz + tcmsz(i)
            end do
        end do
    end subroutine jyusin

    subroutine scale2 ! 系内の全分子の温度の補正
        use parameters
        use variable, !only: velx, vely, velz, zmass
        use molecules_module
        implicit none
        double precision :: temptp, vel2, aimtem, aimnot, baiss
        integer :: i, j

        temptp = 0.0d0
        do i = 1, typmol
            do j = 1, nkoss
                vel2 = mol(j, i)%velx*mol(j, i)%velx + mol(j, i)%vely*mol(j, i)%vely &
                                                     + mol(j, i)%velz*mol(j, i)%velz
                temptp = temptp + vel2
            end do
            temptp = temptp / nkoss / 1.000d+16
            aimtem = atemp1
            aimnot = 3.0d0 * boltz * aimtem / zmass(i)
            baiss = dsqrt(aimnot / temptp)
        end do

        ! 速度ベクトルのスケーリング
        do i = 1, typmol
            do j = 1, nkoss
                mol(j, i)%velx = mol(j, i)%velx * baiss
                mol(j, i)%vely = mol(j, i)%vely * baiss
                mol(j, i)%velz = mol(j, i)%velz * baiss
            end do
        end do
    end subroutine scale2

    subroutine calcu ! 各分子に働く力，速度，位置の分子動力学計算
        use parameters
        use variable
        use molecules_module
        implicit none
        integer :: i, j, i1, i2
        double precision :: divx, divy, divz, dist
        double precision :: dit2, dit4, dit6, dit8, dit12, dit14
        double precision :: ppp, force, forcex, forcey, forcez
        double precision :: vxene, vyene, vzene, vene

        ! 分子間の相互作用力 → ポテンシャルエネルギー
        ! 同じ分子同士の影響
        ! 上段と下段のPt
        do i = 1, typmol
            if(i /= 2) then
                do i1 = 1, nkoss
                    do i2 = i1+1, nkoss
                        divx = mol(i1, i)%posx - mol(i2, i)%posx
                        divy = mol(i1, i)%posy - mol(i2, i)%posy
                        divz = mol(i1, i)%posz - mol(i2, i)%posz
        
                        if (divx < -xcutof(1)) then
                            divx = divx + xsyul(i)
                        else if(divx > xcutof(1)) then
                            divx = divx - xsyul(i)
                        endif
        
                        if (divy < -ycutof(1)) then
                            divy = divy + ysyul(i)
                        else if(divy > ycutof(1)) then
                            divy = divy - ysyul(i)
                        endif
        
                        if (divz < -zcutof(1)) then
                            divz = divz + zsyul(i)
                        else if(divz > zcutof(1)) then
                            divz = divz - zsyul(i)
                        endif
        
                        divx = divx / sig33(1)
        
                        ! if (divx < -cutoff33 .or. divx > cutoff33) then
                        !     cycle
                        ! endif
        
                        ! divy = divy/sig
                        ! if (divy < -cutoff33 .or. divy > cutoff33) then
                        !     cycle
                        ! endif
        
                        ! divz = divz/sig
                        ! if (divz < -cutoff33 .or. divz > cutoff33) then
                        !     cycle
                        ! endif

                        if (divx > cutoff33) then
                            cycle
                        endif
        
                        if (divx < -cutoff33) then
                            cycle
                        endif
        
                        divy = divy/sig33(1)
                        if (divy > cutoff33) then
                            cycle
                        endif
        
                        if (divy < -cutoff33) then
                            cycle
                        endif
        
                        divz = divz/sig33(1)
                        if (divz > cutoff33) then
                            cycle
                        endif
        
                        if (divz < -cutoff33) then
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
                        ppp    = 4.00D0*eps33(1)*(1.00D0/dit12-1.00D0/dit6)
                        force  = cforce(1)*(-2.00D0/dit14+1.00D0/dit8)   ! dit13, dit7じゃない？
                        forcex = -force*divx/zmass(i)
                        forcey = -force*divy/zmass(i)
                        forcez = -force*divz/zmass(i)
                        mol(i1, i)%forx = mol(i1, i)%forx + forcex
                        mol(i2, i)%forx = mol(i2, i)%forx - forcex
                        mol(i1, i)%fory = mol(i1, i)%fory + forcey
                        mol(i2, i)%fory = mol(i2, i)%fory - forcey
                        mol(i1, i)%forz = mol(i1, i)%forz + forcez
                        mol(i2, i)%forz = mol(i2, i)%forz - forcez
                        mol(i1, i)%poten = mol(i1, i)%poten + ppp*0.500D0
                        mol(i2, i)%poten = mol(i2, i)%poten + ppp*0.500D0
                    end do
                end do
            end if
        end do

        !Ar
        do i1 = 1, nkoss
            do i2 = i1+1, nkoss
                divx = mol(i1, 2)%posx - mol(i2, 2)%posx
                divy = mol(i1, 2)%posy - mol(i2, 2)%posy
                divz = mol(i1, 2)%posz - mol(i2, 2)%posz

                if (divx < -xcutof(2)) then
                    divx = divx + xsyul(2)
                else if(divx > xcutof(2)) then
                    divx = divx - xsyul(2)
                endif

                if (divy < -ycutof(2)) then
                    divy = divy + ysyul(2)
                else if(divy > ycutof(2)) then
                    divy = divy - ysyul(2)
                endif

                if (divz < -zcutof(2)) then
                    divz = divz + zsyul(2)
                else if(divz > zcutof(2)) then
                    divz = divz - zsyul(2)
                endif

                divx = divx / sig33(2)

                ! if (divx < -cutoff33 .or. divx > cutoff33) then
                !     cycle
                ! endif

                ! divy = divy/sig
                ! if (divy < -cutoff33 .or. divy > cutoff33) then
                !     cycle
                ! endif

                ! divz = divz/sig
                ! if (divz < -cutoff33 .or. divz > cutoff33) then
                !     cycle
                ! endif

                if (divx > cutoff33) then
                    cycle
                endif

                if (divx < -cutoff33) then
                    cycle
                endif

                divy = divy/sig33(2)
                if (divy > cutoff33) then
                    cycle
                endif

                if (divy < -cutoff33) then
                    cycle
                endif

                divz = divz/sig33(2)
                if (divz > cutoff33) then
                    cycle
                endif

                if (divz < -cutoff33) then
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
                ppp    = 4.00D0*eps33(2)*(1.00D0/dit12-1.00D0/dit6)
                force  = cforce(2)*(-2.00D0/dit14+1.00D0/dit8)   ! dit13, dit7じゃない？
                forcex = -force*divx/zmass(2)
                forcey = -force*divy/zmass(2)
                forcez = -force*divz/zmass(2)
                mol(i1, 2)%forx = mol(i1, 2)%forx + forcex
                mol(i2, 2)%forx = mol(i2, 2)%forx - forcex
                mol(i1, 2)%fory = mol(i1, 2)%fory + forcey
                mol(i2, 2)%fory = mol(i2, 2)%fory - forcey
                mol(i1, 2)%forz = mol(i1, 2)%forz + forcez
                mol(i2, 2)%forz = mol(i2, 2)%forz - forcez
                mol(i1, 2)%poten = mol(i1, 2)%poten + ppp*0.500D0
                mol(i2, 2)%poten = mol(i2, 2)%poten + ppp*0.500D0
            end do
        end do

        ! ! 異なる分子同士の影響
        ! do i = 1, typmol
        !     if (i /= 2) then     ! Ar-Ptの処理    配列は1が上Pt, 2がAr, 3がしたPt
        !         do i1 = 1, nkoss
        !             do i2 = 1, nkoss
        !                 divx = mol(i1, i)%posx - mol(i2, 2)%posx
        !                 divy = mol(i1, i)%posy - mol(i2, 2)%posy
        !                 divz = mol(i1, i)%posz - mol(i2, 2)%posz

        !                 if (divx < -xcutof(3)) then             ! cutof, syulは1がPt-Pt, 2がAr-Ar, 3がPt-Ar
        !                     divx = divx + xsyul(3)
        !                 else if(divx > xcutof(3)) then
        !                     divx = divx - xsyul(3)
        !                 endif

        !                 if (divy < -ycutof(3)) then
        !                     divy = divy + ysyul(3)
        !                 else if(divy > ycutof(3)) then
        !                     divy = divy - ysyul(3)
        !                 endif

        !                 if (divz < -zcutof(3)) then
        !                     divz = divz + zsyul(3)
        !                 else if(divz > zcutof(3)) then
        !                     divz = divz - zsyul(3)
        !                 endif

        !                 divx = divx / sig33(3)

        !                 if (divx < -cutoff33 .or. divx > cutoff33) then
        !                     cycle
        !                 endif

        !                 divy = divy / sig33(3)
        !                 if (divy < -cutoff33 .or. divy > cutoff33) then
        !                     cycle
        !                 endif

        !                 divz = divz / sig33(3)
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
        !                 ppp    = 4.00D0*eps33(3)*(1.00D0/dit12-1.00D0/dit6)
        !                 force  = cforce(3)*(-2.00D0/dit14+1.00D0/dit8)   ! dit13, dit7じゃない？
        !                 forcex = -force*divx/zmass(i)
        !                 forcey = -force*divy/zmass(i)
        !                 forcez = -force*divz/zmass(i)
        !                 mol(i1, i)%forx = mol(i1, i)%forx + forcex
        !                 mol(i2, 2)%forx = mol(i2, 2)%forx - forcex
        !                 mol(i1, i)%fory = mol(i1, i)%fory + forcey
        !                 mol(i2, 2)%fory = mol(i2, 2)%fory - forcey
        !                 mol(i1, i)%forz = mol(i1, i)%forz + forcez
        !                 mol(i2, 2)%forz = mol(i2, 2)%forz - forcez
        !                 mol(i1, i)%poten = mol(i1, i)%poten + ppp*0.500D0
        !                 mol(i2, 2)%poten = mol(i2, 2)%poten + ppp*0.500D0
        !             end do
        !         end do
        !     end if    
        ! end do

        ! 数値積分(蛙跳び法)
        do i = 1, typmol
            do j = 1, nkoss
                vxene = mol(j, i)%velx + mol(j, i)%forx*0.500D0*dt
                vyene = mol(j, i)%vely + mol(j, i)%fory*0.500D0*dt
                vzene = mol(j, i)%velz + mol(j, i)%forz*0.500D0*dt
                vene = vxene*vxene + vyene*vyene + vzene*vzene
                mol(j, i)%ukine = 0.500D0*zmass(i)*vene
            end do
        end do

        do i = 1, typmol
            do j = 1, nkoss
                mol(j, i)%velx = mol(j, i)%velx + mol(j, i)%forx*dt;
                mol(j, i)%vely = mol(j, i)%vely + mol(j, i)%fory*dt;
                mol(j, i)%velz = mol(j, i)%velz + mol(i, i)%forz*dt;
                mol(j, i)%posx = mol(j, i)%posx + mol(j, i)%velx*dt;
                mol(j, i)%posy = mol(j, i)%posy + mol(j, i)%vely*dt;
                mol(j, i)%posz = mol(j, i)%posz + mol(j, i)%velz*dt;
            end do
        end do

    end subroutine calcu

    subroutine bound ! 境界条件の付与
        use parameters
        use variable, only: xsyul, ysyul, zsyul
        use molecules_module
        implicit none
        integer :: i, j
        do i = 1, typmol
            do j = 1, nkoss
                if(mol(j, i)%posx < 0.00D0) then
                    mol(j, i)%posx = mol(j, i)%posx + xsyul(i)
                else if(mol(j, i)%posx > xsyul(i)) then
                    mol(j, i)%posx = mol(j, i)%posx - xsyul(i)
                endif

                if(mol(j, i)%posy < 0.00D0) then
                    mol(j, i)%posy = mol(j, i)%posy + ysyul(i)
                else if(mol(j, i)%posy > ysyul(i)) then
                    mol(j, i)%posy = mol(j, i)%posy - ysyul(i)
                endif

                if(mol(j, i)%posz < 0.00D0) then
                    mol(j, i)%posz = mol(j, i)%posz + zsyul(i)
                else if(mol(j, i)%posz > zsyul(i)) then
                    mol(j, i)%posz = mol(j, i)%posz - zsyul(i)
                endif
            end do
        end do
    end subroutine bound

    subroutine record ! データの出力１
        use parameters
        use variable !only: posx, posy, posz, velx, vely, velz
        use molecules_module
        implicit none
        integer :: i, j
        do i = 1, typmol
            do j = 1, nkoss
                write(3, '(I6, 3D15.7)') j, mol(j, i)%posx, mol(j, i)%posy, mol(j, i)%posz
                write(11,'(3E15.7)') mol(j, i)%posx, mol(j, i)%posy, mol(j, i)%posz
                write(4, '(I6, 3D15.7)') j, mol(j, i)%velx, mol(j, i)%vely, mol(j, i)%velz
            end do
        end do
    end subroutine record

    subroutine record2 ! データの出力２
        use parameters
        use variable !only: poten, ukine
        use molecules_module
        implicit none
        double precision, dimension(typmol) :: totene = 0.00D0, totpot = 0.00D0, totkin = 0.00D0, temp
        integer :: i, j

        ! エネルギーの合計計算
        do i = 1, typmol
            do j = 1, nkoss
                totpot(i) = totpot(i) + mol(j, i)%poten
                totkin(i) = totkin(i) + mol(j, i)%ukine
            end do
        end do

        ! エネルギーを大きな数で割る処理（正規化や単位変換のため）
        do i = 1, typmol
            totpot(i) = totpot(i) / 1.00d16
            totkin(i) = totkin(i) / 1.00d16
            totene(i) = totpot(i) + totkin(i)
        end do

        ! 温度計算
        temp = 2.0d0 * totkin / (3.0d0 * dble(nkoss) * boltz)

        do i = 1, typmol
            write(7, '(4D15.7)') totene(i), totpot(i), totkin(i), temp(i)
        end do
    end subroutine record2

    subroutine record3
        use parameters
        use variable !only: posx, posy, posz, velx, vely, velz
        use molecules_module
        implicit none
        integer :: i, j
        do i = 1, typmol
            do j = 1, nkoss
                write(9, '(I6, 6D15.7)') j, mol(j, i)%posx, mol(j, i)%posy, mol(j, i)%posz, &
                                            mol(j, i)%velx, mol(j, i)%vely, mol(j, i)%velz
            end do
        end do
        
    end subroutine record3
end program main