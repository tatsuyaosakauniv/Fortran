!-----------------------------------------------------------------------
! example.f90
! ナノ粒子の初期状態を分子動力学法で作成する．
! 必要ファイル: random0.dat, random1.dat
!-----------------------------------------------------------------------

! パラメータ（もともとparate.datにあったもの）
module parameters
    implicit none
    integer, parameter :: nkoss = 1200 ! 全分子数
    ! 周期長さの単位がよくわからない → [Å]?    3.92*3.92*8.70[nm^3]にしたい
    double precision, parameter :: xsyul0 = 39.200D0 ! x方向の周期境界長さ(無次元)
    double precision, parameter :: ysyul0 = 39.200D0 ! y方向の周期境界長さ(無次元）
    double precision, parameter :: zsyul0 = 87.000D0 ! z方向の周期境界長さ(無次元）
    ! アルゴンの分子の質量39.95[kg/mol]
    double precision, parameter :: bunsi3(2) = (/39.950D-3, 195.084D-3/) ! 分子の質量  kg/mol
    ! Ar-Ar σ=3.40[Å], Pt-Pt σ=2.54[Å], Ar-Pt σ=2.97[Å]
    double precision, parameter :: sig33(3) = (/3.40D0, 2.54D0, 2.97D0/) ! 分子のLennard-Jonesパラメータσ(無次元)
    ! Ar-Ar ε=1.67[10^(-21)J], Pt-Pt ε=109.2[10^(-21)J], Ar-Pt ε=13.50[10^(-21)J]
    double precision, parameter :: eps33(3) = (/1.67D-21, 1.092D-19, 1.350D-20/) ! 分子のLennard-Jonesパラメータε(無次元)
    double precision, parameter :: atemp1 = 150D0 ! 系内（目標）温度  K
    double precision, parameter :: cutoff33 = 3.300D0 ! カットオフ長さ/σ
    integer, parameter :: maxstep = 500 ! 最大ステップ数
    double precision, parameter :: avoga = 6.022D+23 ! アボガドロ数
    double precision, parameter :: boltz = 1.3806662D-23 ! ボルツマン定数
    double precision, parameter :: pi = 3.141592654D0 ! 円周率
    double precision, parameter :: dt = 5.00D0 ! 無次元時間ステップ(無次元，有次元の場合fs)
end module parameters

! 変数
module variable
    use parameters
    implicit none
    !double precision, dimension(nkoss) :: posx, posy, posz, ptempx, ptempy, ptempz ! 各分子の座標(x,y,z方向）
    !double precision, dimension(nkoss) :: velx, vely, velz, vtempx, vtempy, vtempz ! 各分子の速度（x,y,z方向）
    !double precision, dimension(nkoss) :: forx = 0.0000D0, fory = 0.0000D0, forz = 0.0000D0 ! 各分子に働く力（x,y,z方向）
    !double precision, dimension(nkoss) :: poten = 0.0000D0, ukine = 0.0000D0 ! 各分子のポテンシャルエネルギー，運動エネルギー
    double precision :: zmass, cforce ! 分子の質量，力の計算の係数
    double precision :: sig, eps ! 分子のLennard-Jonesパラメータ，σ，ε
    double precision :: xcutof, ycutof, zcutof ! ポテンシャルのカットオフ長さx,y,z方向
    integer :: nowstp ! 現在のステップ数
    double precision :: xsyul, ysyul, zsyul ! ポテンシャルのカットオフ長さx,y,z方向，x,y,z方向の周期長さ
end module variable

module molecules_module
    implicit none

    !構造体の定義
    type :: mol_type
        double precision :: posx, posy, posz
        double precision :: velx, vely, velz
        double precision :: forx  = 0.0000D0, fory  = 0.0000D0, forz  = 0.0000D0
        double precision :: poten = 0.0000D0, ukine = 0.0000D0
    end type mol_type

    type(mol_type), dimension(1200 , 3) :: mol  !上Pt, 中Ar, 下Pt

end module molecules_module

program main
    use variable, only: nowstp
    use parameters
    use molecules_module
    implicit none

    integer :: i

    ! 読み込み用乱数ファイル
        open(1,file='random1.dat')
        open(2,file='random0.dat')
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


    nowstp = 0

    call seting ! 各分子の初期位置，初期速度などの設定
    !call cortra ! 系内の全分子の並進速度の補正
    !call jyusin ! 系内の全分子の重心の補正
    !call scale ! 系内の全分子の温度の補正
    call record ! データの出力１
    !call record2 ! データの出力２

    do i = 1, maxstep
        nowstp = i
        ! ステップ数が500の倍数のとき
        if (mod(nowstp,500) == 0) then
            write(6,*) nowstp
        endif
        ! ステップ数が100の倍数のとき
        !if (mod(nowstp,100) == 0) then
        !	call cortra	! 系内の全分子の並進速度の補正
        !	call scale ! 系内の全分子の温度の補正
        !	call jyusin	! 系内の全分子の重心の補正
        !endif

        !call calcu ! 各分子に働く力，速度，位置の分子動力学計算
        !call bound ! 境界条件の付与
        
        ! ステップ数が100の倍数+1のとき
        !if(mod(nowstp, 100) == 1) then
         !  call record ! データの出力１
         !  call record2 ! データの出力２
        !endif
    end do

    call record3 ! データの出力３

contains
    subroutine seting ! 各分子の初期位置，初期速度などの設定
        use variable
        use parameters
        implicit none
        integer :: num, i, j, k
        double precision :: x = 0.0000D0, y = 0.0000D0, z = 0.0000D0
        double precision :: ofstx1 = 0.0000D0, ofsty1 = 0.0000D0, ofstz1 = 80.0000D0
        double precision :: stdist, ran, alpha, beta, cr
        double precision :: vx, vy, vz

        zmass = 1.00D+26*bunsi3/avoga
        sig   = sig33(1)
        eps   = eps33(1)
        cforce = 24.00D0*eps/sig
        xsyul = xsyul0
        ysyul = ysyul0
        zsyul = zsyul0
        write(9,*)xsyul
        write(9,*)ysyul
        write(9,*)zsyul

        num = 0

        xcutof = xsyul - cutoff33*sig33(1)
        ycutof = ysyul - cutoff33*sig33(1)
        zcutof = zsyul - cutoff33*sig33(1)

        stdist = 4.00D0     !原子間距離
        do k = 1,6
            z = ofstz1 + dble(k-1)*stdist*0.5D0
            do i = 1,20
                x = ofstx1 + dble(i-1)*stdist*0.5D0
                do j = 1,10
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
                    !write(10,*)posz(num)
                end do
            end do
        end do

        num = 0
        cr = 1.00D-6
        do i = 1, nkoss
            read(1,*)ran
            alpha = pi*ran
            read(2,*)ran
            beta = 2.000D0*pi*ran
            vx = 0.000D0 !dsin(alpha)*dcos(beta)*cr
            vy = 0.000D0 !dsin(alpha)*dsin(beta)*cr
            vz = 0.000D0 !dcos(alpha)*cr
            num = num + 1
            velx(num) = vx
            vely(num) = vy
            velz(num) = vz
        end do
    end subroutine seting

    subroutine cortra ! 系内の全分子の並進速度の補正
        use variable, only: velx, vely, velz
        use parameters
        implicit none
        double precision :: trvx = 0.0D0, trvy = 0.0D0, trvz = 0.0D0
        integer :: i, j

        ! 速度ベクトルの成分の平均を計算
        do i = 1, nkoss
            trvx = trvx + velx(i)
            trvy = trvy + vely(i)
            trvz = trvz + velz(i)
        end do

        trvx = trvx / nkoss
        trvy = trvy / nkoss
        trvz = trvz / nkoss

        ! 速度ベクトルから平均を引いて中心補正
        do j = 1, nkoss
            velx(j) = velx(j) - trvx
            vely(j) = vely(j) - trvy
            velz(j) = velz(j) - trvz
        end do
    end subroutine cortra

    subroutine jyusin ! 系内の全分子の重心の補正
        use variable, only: posx, posy, posz
        use parameters
        implicit none
        double precision :: cmsx, cmsy, cmsz
        double precision :: tcmsx = 0.0000D0, tcmsy = 0.0000D0, tcmsz = 0.0000D0
        integer :: i, j

        cmsx = xsyul0 * 0.500D0
        cmsy = ysyul0 * 0.500D0
        cmsz = zsyul0 * 0.500D0

        do i = 1, nkoss
            tcmsx = tcmsx + posx(i)
            tcmsy = tcmsy + posy(i)
            tcmsz = tcmsz + posz(i)
        end do

        tcmsx = cmsx - tcmsx/dble(nkoss)
        tcmsy = cmsy - tcmsy/dble(nkoss)
        tcmsz = cmsz - tcmsz/dble(nkoss)
        
        do j = 1, nkoss
            posx(j) = posx(j) + tcmsx
            posy(j) = posy(j) + tcmsy
            posz(j) = posz(j) + tcmsz
        end do
    end subroutine jyusin

    subroutine scale ! 系内の全分子の温度の補正
        use variable, only: velx, vely, velz, zmass
        use parameters
        implicit none
        double precision :: temptp, vel2, aimtem, aimnot, baiss
        integer :: i

        temptp = 0.0d0
        do i = 1, nkoss
            vel2 = velx(i)*velx(i) + vely(i)*vely(i) + velz(i)*velz(i)
            temptp = temptp + vel2
        end do

        temptp = temptp / nkoss / 1.000d+16
        aimtem = atemp1
        aimnot = 3.0d0 * boltz * aimtem / zmass
        baiss = dsqrt(aimnot / temptp)

        ! 速度ベクトルのスケーリング
        do i = 1, nkoss
            velx(i) = velx(i) * baiss
            vely(i) = vely(i) * baiss
            velz(i) = velz(i) * baiss
        end do
    end subroutine scale

    subroutine calcu ! 各分子に働く力，速度，位置の分子動力学計算
        use variable
        use parameters
        implicit none
        integer :: i, i1, i2
        double precision :: divx, divy, divz, dist
        double precision :: dit2, dit4, dit6, dit8, dit12, dit14
        double precision :: ppp, force, forcex, forcey, forcez
        double precision :: vxene, vyene, vzene, vene

        do i1 = 1, nkoss
            do i2 = i1+1, nkoss
                divx = posx(i1)-posx(i2)
                divy = posy(i1)-posy(i2)
                divz = posz(i1)-posz(i2)

                if (divx < -xcutof) then
                    divx = divx + xsyul
                else if(divx > xcutof) then
                    divx = divx - xsyul
                endif

                if (divy < -ycutof) then
                    divy = divy + ysyul
                else if(divy > ycutof) then
                    divy = divy - ysyul
                endif

                if (divz < -zcutof) then
                    divz = divz + zsyul
                else if(divz > zcutof) then
                    divz = divz - zsyul
                endif

                divx = divx / sig

                if (divx < -cutoff33 .or. divx > cutoff33) then
                    cycle
                endif

                divy = divy/sig
                if (divy < -cutoff33 .or. divy > cutoff33) then
                    cycle
                endif

                divz = divz/sig
                if (divz < -cutoff33 .or. divz > cutoff33) then
                    cycle
                endif

                dit2 = divx*divx + divy*divy + divz*divz
                dist = dsqrt(dit2)

                if(dist > cutoff33) then
                    cycle
                endif

                dit4 = dit2*dit2
                dit6 = dit4*dit2
                dit8 = dit4*dit4
                dit12 = dit6*dit6
                dit14 = dit8*dit6
                ppp = 4.00D0*eps*(1.00D0/dit12-1.00D0/dit6)
                force = cforce*(-2.00D0/dit14+1.00D0/dit8)
                forcex = -force*divx/zmass
                forcey = -force*divy/zmass
                forcez = -force*divz/zmass
                forx(i1) = forx(i1) + forcex
                forx(i2) = forx(i2) - forcex
                fory(i1) = fory(i1) + forcey
                fory(i2) = fory(i2) - forcey
                forz(i1) = forz(i1) + forcez
                forz(i2) = forz(i2) - forcez
                poten(i1) = poten(i1) + ppp*0.500D0
                poten(i2) = poten(i2) + ppp*0.500D0
            end do
        end do

        do i = 1, nkoss
            vxene = velx(i) + forx(i)*0.500D0*dt
            vyene = vely(i) + fory(i)*0.500D0*dt
            vzene = velz(i) + forz(i)*0.500D0*dt
            vene = vxene*vxene + vyene*vyene + vzene*vzene
            ukine(i) = 0.500D0*zmass*vene
        end do

        do i = 1, nkoss
            posx(i) = posx(i) + velx(i)*dt;
            posy(i) = posy(i) + vely(i)*dt;
            posz(i) = posz(i) + velz(i)*dt;
            velx(i) = velx(i) + forx(i)*dt;
            vely(i) = vely(i) + fory(i)*dt;
            velz(i) = velz(i) + forz(i)*dt;
        end do

    end subroutine calcu

    subroutine bound ! 境界条件の付与
        use variable, only: posx, posy, posz, xsyul, ysyul, zsyul
        use parameters
        implicit none
        integer :: i
        do i = 1, nkoss
            if(posx(i) < 0.00D0) then
                posx(i) = posx(i) + xsyul
            else if(posx(i) > xsyul) then
                posx(i) = posx(i) - xsyul
            endif

            if(posy(i) < 0.00D0) then
                posy(i) = posy(i) + ysyul
            else if(posy(i) > ysyul) then
                posy(i) = posy(i) - ysyul
            endif

            if(posz(i) < 0.00D0) then
                posz(i) = posz(i) + zsyul
            else if(posz(i) > zsyul) then
                posz(i) = posz(i) - zsyul
            endif
        end do
    end subroutine bound

    subroutine record ! データの出力１
        use variable, only: posx, posy, posz, velx, vely, velz
        use parameters
        implicit none
        integer :: i
        do i = 1, nkoss
            write(3, '(I6, 3D15.7)') i, posx(i), posy(i), posz(i) !- 40.00D0
            if(i == nkoss) then
                write(6,*) nkoss
                write(6,*) 'ok'
            end if
            write(4, '(I6, 3D15.7)') i, velx(i), vely(i), velz(i)
        end do
    end subroutine record

    subroutine record2 ! データの出力２
        use variable, only: poten, ukine
        use parameters
        implicit none
        double precision :: totene = 0.00D0, totpot = 0.00D0, totkin = 0.00D0, temp
        integer :: i

        ! エネルギーの合計計算
        do i = 1, nkoss
            totpot = totpot + poten(i)
            totkin = totkin + ukine(i)
        end do

        ! エネルギーを大きな数で割る処理（正規化や単位変換のため）
        totpot = totpot / 1.00d16
        totkin = totkin / 1.00d16
        totene = totpot + totkin

        ! 温度計算
        temp = 2.0d0 * totkin / (3.0d0 * dble(nkoss) * boltz)

        write(7, '(4D15.7)') totene, totpot, totkin, temp
    end subroutine record2

    subroutine record3
        use variable, only: posx, posy, posz, velx, vely, velz
        use parameters
        implicit none
        integer :: i
        do i = 1, nkoss
            write(9, '(I6, 6D15.7)') i, posx(i), posy(i), posz(i), velx(i), vely(i), velz(i)
        end do
    end subroutine record3
end program main