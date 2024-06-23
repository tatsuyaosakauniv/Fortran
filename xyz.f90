!-----------------------------------------------------------------------
! example.f90
! ナノ粒子の初期状態を分子動力学法で作成する．
! 必要ファイル: random0.dat, random1.dat
!-----------------------------------------------------------------------

! パラメータ（もともとparate.datにあったもの）
module parameters
    implicit none
    integer, parameter :: stpMax = 20000 ! 最大ステップ数
    double precision, parameter :: dt = 1.00d0 ! 無次元時間ステップ(無次元，有次元の場合fs)
    double precision, parameter :: tempAr = 100d0 ! 系内（目標）温度  K
    double precision, parameter :: angCon = 0.05d0 ! 接触角

    integer, parameter :: TYPMOL = 3 ! 分子の種類数
    integer, parameter :: COMP = 3 ! 相互作用の組み合わせ数 = TYPMOLC2
    integer, parameter :: numx(TYPMOL) = [10, 6, 10]
    integer, parameter :: numy(TYPMOL) = [ 5, 3,  5]
    integer, parameter :: numz(TYPMOL) = [ 4, 10, 4]
    integer, parameter :: nummol(TYPMOL) = [numx(1)*numy(1)*numz(1), numx(2)*numy(2)*numz(2), numx(3)*numy(3)*numz(3)] ! 各分子の数
    

    double precision, parameter :: STDIST(TYPMOL) = [3.92d0, 6.0d0, 3.92d0] ! 格子定数(無次元)[Å]
    double precision, parameter :: xsyul0 = STDIST(1) * numy(1) ! x方向の周期境界長さ(無次元)
    double precision, parameter :: ysyul0 = STDIST(1) * numy(1) ! y方向の周期境界長さ(無次元）
    double precision, parameter :: zsyul0 = 60.000d0 ! z方向の周期境界長さ(無次元）
    
    double precision, parameter :: CUTOFF = 3.300d0 ! カットオフ長さ/σ
    double precision, parameter :: AVOGA = 6.022d+23 ! アボガドロ数
    double precision, parameter :: BOLTZ = 1.3806662d-23 ! ボルツマン定数 [J/K]
    double precision, parameter :: PI = 3.141592654d0 ! 円周率

    ! 分子の質量
    !double precision, parameter :: bunsi(TYPMOL) = [195.084d-3, 39.950d-3, 195.084d-3] ! 分子の質量  kg/mol   
    double precision, parameter :: MASS(TYPMOL) = [32.395d0, 6.6340d0, 32.395d0] ! 分子の質量（無次元） * 10d-26 [kg/個]
    ! Lennard-Jonesパラメータ
    double precision, parameter :: SIG(COMP+1) = [2.475d0, 3.4000d0, 2.4750d0, 2.9375d0]  ! σ(無次元)
    double precision, parameter :: EPS(COMP+1) = [83.31d-5, 1.666d-5, 83.31d-5, 11.78d-5] ! ε(無次元)

    ! Langevin法用
    double precision, parameter :: tempLanPt(TYPMOL) = [150d0, 0d0, 50d0] ! Langevin法を用いるPtの温度  真ん中は使わない
    double precision, parameter :: DIRAC = 1.054571817d-34 ! ディラック定数 [J･s]
    double precision, parameter :: DEBTMP = 240d0 ! Debye温度 [K]
    double precision, parameter :: OMEGA = BOLTZ * DEBTMP / DIRAC * 1.000d-11 ! Debye定数(無次元)
    double precision, parameter :: DAMP =  MASS(1) * PI * OMEGA / 6.000d0 ! ダンパーの減衰係数(無次元)

end module parameters

! 変数
module variable
    use parameters
    implicit none
    double precision :: forCoef(4) = [0.0d0, 0.0d0, 0.0d0, 0.0d0] ! 力の計算の係数
    double precision :: xcutof, ycutof, zcutof ! ポテンシャルのカットオフ長さx,y,z方向
    integer :: stpNow ! 現在のステップ数
    double precision :: xsyul, ysyul, zsyul ! ポテンシャルのカットオフ長さx,y,z方向，x,y,z方向の周期長さ
    double precision :: stddev ! 標準偏差
    double precision, allocatable :: rforce(:,:)   ! ランダム力用
    double precision, allocatable :: dforce(:,:)   ! ランダム力用
    double precision :: rnd_2
    logical :: isOdd = .true. ! 乱数のsinとcosを交互に出すためのフラグ

    contains

    function Random() result(rnd_) ! 乱数生成用の関数　　呼び出す時に何回も計算しなくていいように下の関数と分けた
        use parameters
        implicit none
        double precision :: rnd_, u1, u2
        
        ! 呼び出される回数が偶数か奇数かによってsinかcosかを使い分ける
        if(isOdd) then
            call random_number(u1)
            call random_number(u2)
            rnd_2 = sqrt(-2.0d0 * log(u1)) * sin(2.0d0 * PI * u2)
            rnd_  = sqrt(-2.0d0 * log(u1)) * cos(2.0d0 * PI * u2)
            isOdd = .false.
        else
           rnd_ = rnd_2
           isOdd = .true.
        end if

    end function random

    function getStddev(T_) result(stddev_) ! 設定温度を引数とし、ランダム力計算用の標準偏差を出力する関数
        use parameters
        implicit none
        double precision, intent(in) :: T_
        double precision :: stddev_

        stddev_ = dsqrt(2.000d0 * BOLTZ * T_ * DAMP / dt)  ! 無次元

    end function getStddev
end module variable

module molecules_struct
    use parameters
    implicit none

    !構造体の定義
    type :: mol_info
        double precision :: posx, posy, posz
        double precision :: velx, vely, velz
        double precision :: accx, accy, accz  
        double precision :: poten, kinet 
    end type mol_info

    type mol_typ
        type(mol_info), allocatable :: mol(:)
    end type mol_typ

    type(mol_typ), dimension(TYPMOL) :: typ  ! 上Pt, 中Ar, 下Pt

end module molecules_struct

module forPVwin
    use parameters
    implicit none
    integer, parameter :: tlnkoss = nummol(1) + nummol(2) + nummol(3)
    integer, parameter :: moltype = 1
    integer, parameter :: ndat = int(stpMax/100)
    integer, parameter :: ntime0 = 0
    integer, parameter :: ndt = 1
end module forPVwin

program main
    use variable, only: stpNow
    use parameters
    use molecules_struct
    use forPVwin
    implicit none
    integer :: i, j

    ! 配列初期化
    allocate(typ(1)%mol(nummol(1)))
    allocate(typ(2)%mol(nummol(2)))
    allocate(typ(3)%mol(nummol(3)))

    ! 読み込み用乱数ファイル
        open(1,file='random0.dat', status='old')
        open(2,file='random1.dat', status='old')
        !open(6,*) は使用できない
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
        open(35,file='energy_all.dat')
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
        do j = 1, int(nummol(1)/numz(1))
            write(90,'(I7)') 15      ! 白色
        end do
        do j = int(nummol(1)/numz(1)) + 1, nummol(1)
            write(90,'(I7)') 14      ! 赤色
        end do
        do j = 1, nummol(2)
            write(90,'(I7)') 7       ! 黄色
        end do
        do j = 1, int(nummol(3)/numz(3))
            write(90,'(I7)') 15      ! 白色
        end do
        do j = int(nummol(3)/numz(3)) + 1, nummol(3)
            write(90,'(I7)') 0       ! 青色
        end do
    end do
    
    stpNow = 0

    call seting ! 各分子の初期位置，初期速度などの設定
    !call cortra ! 系内の全分子の並進速度の補正
    !call jyusin ! 系内の全分子の重心の補正
    !call scale ! 系内の全分子の温度の補正
    !call record_pos_vel! データの出力１
    !call record_energy_temp ! データの出力２

    do i = 1, stpMax
        stpNow = i
        ! ステップ数が500の倍数のとき
        if (mod(stpNow,500) == 0) then
            write(6,'(3X, I7, 1X, A3, I7)') stpNow, ' / ', stpMax
        endif
        ! ステップ数が100の倍数のとき
        if (stpNow <= int(stpMax*0.3) .and. mod(stpNow,100) == 0) then
        	!call cortra	! 系内の全分子の並進速度の補正
            !call jyusin	! 系内の全分子の重心の補正
        	call scale ! 系内の全分子の温度の補正
        endif

        call calcu ! 各分子に働く力，速度，位置の分子動力学計算
        call bound ! 境界条件の付与
        
        ! ステップ数が100の倍数+1のとき
        if(mod(stpNow, 100) == 1) then
          call record_pos_vel! データの出力１
          call record_energy_temp ! データの出力２
        endif
    end do

    call record_finpos_vel ! データの出力３

    contains
    subroutine seting ! 各分子の初期位置，初期速度などの設定
        use variable
        use parameters
        use molecules_struct
        implicit none
        integer :: num, i, j, k
        double precision :: x = 0.0000d0, y = 0.0000d0, z = 0.0000d0
        double precision :: ofstx, ofsty, ofstz
        double precision :: ran, alpha, beta, cr
        double precision :: vx, vy, vz

        do i = 1, COMP
            forCoef(i) = 24.00d0*EPS(i)/SIG(i)
        end do
        forCoef(4) = angCon*24.00d0*EPS(4)/SIG(4)
        
        xsyul = xsyul0
        ysyul = ysyul0
        zsyul = zsyul0
        write(50,*)xsyul
        write(50,*)ysyul
        write(50,*)zsyul
        write(15,*)xsyul0, ysyul0, zsyul0
        write(15,*)ntime0, ndt

        num = 0

        do i = 1, TYPMOL
            xcutof = xsyul - CUTOFF*SIG(i)
            ycutof = ysyul - CUTOFF*SIG(i)
            zcutof = zsyul - CUTOFF*SIG(i)
        end do

        !上段のPt配置
        ofstx = STDIST(1)*0.25d0
        ofsty = STDIST(1)*0.25d0
        ofstz = zsyul0 - STDIST(1)*0.25d0
        do k = 1,numz(1)
            z = ofstz - dble(k-1)*STDIST(1)*0.5d0
            do i = 1,numx(1)
                x = ofstx + dble(i-1)*STDIST(1)*0.5d0
                do j = 1,numy(1)
                    if(mod(k,2) == 0) then  !z偶数
                        if(mod(i,2) == 0) then
                            y = ofsty + dble(j-1)*STDIST(1)   !x偶数
                        else
                            y = ofsty + dble(j-1)*STDIST(1) + STDIST(1)*0.5d0    !x奇数
                        endif
                    else    !z奇数
                        if(mod(i,2) == 0) then
                            y = ofsty + dble(j-1)*STDIST(1) + STDIST(1)*0.5d0    !x偶数
                        else
                            y = ofsty + dble(j-1)*STDIST(1)   !x奇数
                        endif
                    endif
                    num = num + 1
                    typ(1)%mol(num)%posx = x
                    typ(1)%mol(num)%posy = y
                    typ(1)%mol(num)%posz = z
                end do
            end do
        end do

        !中段のAr配置
        num = 0
        ofstx = xsyul0*0.50d0 - STDIST(2)*(0.25d0*(numx(2)  -1))
        ofsty = ysyul0*0.50d0 - STDIST(2)*(0.25d0*(numy(2)*2-1))
        ofstz = zsyul0*0.50d0 - STDIST(2)*(0.25d0*(numz(2)  -1))
        do k = 1,numz(2)
            z = ofstz + dble(k-1)*STDIST(2)*0.5d0
            do i = 1,numx(2)
                x = ofstx + dble(i-1)*STDIST(2)*0.5d0
                do j = 1,numy(2)
                    if(mod(k,2) == 0) then  !z偶数
                        if(mod(i,2) == 0) then
                            y = ofsty + dble(j-1)*STDIST(2)   !x偶数
                        else
                            y = ofsty + dble(j-1)*STDIST(2) + STDIST(2)*0.5d0    !x奇数
                        endif
                    else    !z奇数
                        if(mod(i,2) == 0) then
                            y = ofsty + dble(j-1)*STDIST(2) + STDIST(2)*0.5d0    !x偶数
                        else
                            y = ofsty + dble(j-1)*STDIST(2)   !x奇数
                        endif
                    endif
                    num = num + 1
                    typ(2)%mol(num)%posx = x
                    typ(2)%mol(num)%posy = y
                    typ(2)%mol(num)%posz = z
                end do
            end do
        end do

        !下段のPt配置
        num = 0
        ofstx = STDIST(3)*0.25d0
        ofsty = STDIST(3)*0.25d0
        ofstz = STDIST(3)*0.25d0
        do k = 1,numz(3)
            z = ofstz + dble(k-1)*STDIST(3)*0.5d0
            do i = 1,numx(3)
                x = ofstx + dble(i-1)*STDIST(3)*0.5d0
                do j = 1,numy(3)
                    if(mod(k,2) == 0) then  !z偶数
                        if(mod(i,2) == 0) then
                            y = ofsty + dble(j-1)*STDIST(3)   !x偶数
                        else
                            y = ofsty + dble(j-1)*STDIST(3) + STDIST(3)*0.5d0    !x奇数
                        endif
                    else    !z奇数
                        if(mod(i,2) == 0) then
                            y = ofsty + dble(j-1)*STDIST(3) + STDIST(3)*0.5d0    !x偶数
                        else
                            y = ofsty + dble(j-1)*STDIST(3)   !x奇数
                        endif
                    endif
                    num = num + 1
                    typ(3)%mol(num)%posx = x
                    typ(3)%mol(num)%posy = y
                    typ(3)%mol(num)%posz = z
                end do
            end do
        end do

        cr = 1.00d-6
        do j = 1, TYPMOL
            do i = 1, nummol(j)
                !read(1,*)ran
                call random_number(ran)
                alpha = PI*ran
                !read(2,*)ran
                call random_number(ran)
                beta = 2.000d0*PI*ran
                vx = dsin(alpha)*dcos(beta)*cr
                vy = dsin(alpha)*dsin(beta)*cr
                vz = dcos(alpha)*cr
                typ(j)%mol(i)%velx = vx
                typ(j)%mol(i)%vely = vy
                typ(j)%mol(i)%velz = vz
            end do
        end do

        do j = 1, TYPMOL
            do i = 1, int(nummol(j)/numz(j))
                if(j == 2) then
                    cycle
                else
                    typ(j)%mol(i)%velx = 0.000d0
                    typ(j)%mol(i)%vely = 0.000d0
                    typ(j)%mol(i)%velz = 0.000d0
                end if
            end do
        end do
    end subroutine seting

    subroutine cortra ! 系内の全分子の並進速度の補正
        use variable, !only: velx, vely, velz
        use parameters
        implicit none
        double precision :: trvx , trvy, trvz
        integer :: i, j

        ! 速度ベクトルの成分の平均を計算
        do j = 1, TYPMOL
            trvx = 0.0d0
            trvy = 0.0d0
            trvz = 0.0d0
            do i = 1, nummol(j)
                trvx = trvx + typ(j)%mol(i)%velx
                trvy = trvy + typ(j)%mol(i)%vely
                trvz = trvz + typ(j)%mol(i)%velz
            end do

            trvx = trvx / nummol(j)
            trvy = trvy / nummol(j)
            trvz = trvz / nummol(j)

            ! 速度ベクトルから平均を引いて中心補正
            do i = 1, nummol(j)
                typ(j)%mol(i)%velx = typ(j)%mol(i)%velx - trvx
                typ(j)%mol(i)%vely = typ(j)%mol(i)%vely - trvy
                typ(j)%mol(i)%velz = typ(j)%mol(i)%velz - trvz
            end do
        end do
    end subroutine cortra

    subroutine jyusin ! 系内の全分子の重心の補正
        use variable, !only: posx, posy, posz
        use parameters
        use molecules_struct
        implicit none
        double precision :: cmsx, cmsy, cmsz
        double precision :: tcmsx, tcmsy, tcmsz
        integer :: i, j

        cmsx = xsyul0 * 0.500d0
        cmsy = ysyul0 * 0.500d0
        cmsz = zsyul0 * 0.500d0

        do j = 1, TYPMOL
            tcmsx = 0.0000d0
            tcmsy = 0.0000d0
            tcmsz = 0.0000d0
            do i = 1, nummol(j)
                tcmsx = tcmsx + typ(j)%mol(i)%posx
                tcmsy = tcmsy + typ(j)%mol(i)%posy
                tcmsz = tcmsz + typ(j)%mol(i)%posz
            end do

            tcmsx = cmsx - tcmsx/dble(nummol(j))
            tcmsy = cmsy - tcmsy/dble(nummol(j))
            tcmsz = cmsz - tcmsz/dble(nummol(j))

            do i = 1, nummol(j)
                typ(j)%mol(i)%posx = typ(j)%mol(i)%posx + tcmsx
                typ(j)%mol(i)%posy = typ(j)%mol(i)%posy + tcmsy
                typ(j)%mol(i)%posz = typ(j)%mol(i)%posz + tcmsz
            end do
        end do
    end subroutine jyusin

    subroutine scale ! 系内の全分子の温度の補正
        use variable, !only: velx, vely, velz, MASS
        use parameters
        use molecules_struct
        implicit none
        double precision :: temptp, vel2, aimtem, aimnot, baiss
        integer :: i
        integer :: j = 2 ! Arのみ

        temptp = 0.000d0
        do i = 1, nummol(j)
            vel2 = typ(j)%mol(i)%velx**2 + typ(j)%mol(i)%vely**2 + typ(j)%mol(i)%velz**2
            temptp = temptp + vel2
        end do
        temptp = temptp / nummol(j) * 1.000d-16
        aimtem = tempAr
        aimnot = 3.000d0 * BOLTZ * aimtem / MASS(j)
        baiss = dsqrt(aimnot / temptp)

        ! 速度ベクトルのスケーリング
        ! Arのみ
        do i = 1, nummol(j)
            typ(j)%mol(i)%velx = typ(j)%mol(i)%velx * baiss
            typ(j)%mol(i)%vely = typ(j)%mol(i)%vely * baiss
            typ(j)%mol(i)%velz = typ(j)%mol(i)%velz * baiss
        end do
    end subroutine scale

    subroutine calcu ! 各分子に働く力，速度，位置の分子動力学計算
        use variable
        use parameters
        use molecules_struct
        implicit none
        integer :: i, j, k, i1, i2
        double precision :: divx, divy, divz, dist
        double precision :: dit2, dit4, dit6, dit8, dit12, dit14
        double precision :: ppp, force, accelx, accely, accelz
        double precision :: vxene, vyene, vzene, vene
        double precision :: rnd

        do j = 1, TYPMOL
            do i = 1, nummol(j)
                typ(j)%mol(i)%accx  = 0.0000d0
                typ(j)%mol(i)%accy  = 0.0000d0
                typ(j)%mol(i)%accz  = 0.0000d0
                typ(j)%mol(i)%poten = 0.0000d0
                typ(j)%mol(i)%kinet = 0.0000d0
            end do
        end do

        do j = 1, TYPMOL
            if(j == 2) then
                cycle
            end if

            !-----------------------ランダム力を入れるとエネルギーが保存しない---------------------
            !ランダム力生成
            allocate(rforce(nummol(j), 3))
            allocate(dforce(nummol(j), 3))
            do k = 1, 3     !xyz
                do i = int(nummol(j)/numz(j)) + 1, int(nummol(j) *2/numz(j))   !Phantom層のみ
                    rnd = Random()
                    rforce(i,k) = rnd * getStddev(tempLanPt(j))
                end do
            end do

            do i = int(nummol(j)/numz(j)) + 1, int(nummol(j) *2/numz(j))
                dforce(i,1) = - DAMP*typ(j)%mol(i)%velx * 1.000d-10
                dforce(i,2) = - DAMP*typ(j)%mol(i)%vely * 1.000d-10
                dforce(i,3) = - DAMP*typ(j)%mol(i)%velz * 1.000d-10
            end do

            ! ランダム力とダンパー力を追加
            do i = int(nummol(j)/numz(j)) + 1, int(nummol(j) *2/numz(j))
                typ(j)%mol(i)%accx = (rforce(i,1) + dforce(i,1)) / MASS(j) * 1.000d+6
                typ(j)%mol(i)%accy = (rforce(i,2) + dforce(i,2)) / MASS(j) * 1.000d+6
                typ(j)%mol(i)%accz = (rforce(i,3) + dforce(i,3)) / MASS(j) * 1.000d+6
            end do

            deallocate(rforce)
            deallocate(dforce)
        end do

        ! 分子間の相互作用力 → ポテンシャルエネルギー
        ! 同じ分子同士の影響
        do i = 1, TYPMOL
            do i1 = 1, nummol(i)
                do i2 = i1+1, nummol(i)
                    divx = typ(i)%mol(i1)%posx - typ(i)%mol(i2)%posx
                    divy = typ(i)%mol(i1)%posy - typ(i)%mol(i2)%posy
                    divz = typ(i)%mol(i1)%posz - typ(i)%mol(i2)%posz
    
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
    
                    divx = divx / SIG(i)
                    divy = divy / SIG(i)
                    divz = divz / SIG(i)

                    if (divx < -CUTOFF .or. divx > CUTOFF) then
                        cycle
                    endif
    
                    if (divy < -CUTOFF .or. divy > CUTOFF) then
                        cycle
                    endif
    
                    if (divz < -CUTOFF .or. divz > CUTOFF) then
                        cycle
                    endif

                    dit2 = divx*divx + divy*divy + divz*divz
                    dist = dsqrt(dit2)
    
                    if(dist > CUTOFF) then
                        cycle
                    endif
    
                    dit4   = dit2*dit2
                    dit6   = dit4*dit2
                    dit8   = dit4*dit4
                    dit12  = dit6*dit6
                    dit14  = dit8*dit6
                    ppp    = 4.00d0*EPS(i)*(1.00d0/dit12-1.00d0/dit6)
                    force  = forCoef(i)*(-2.00d0/dit14+1.00d0/dit8)
                    accelx = -force*divx/MASS(i)
                    accely = -force*divy/MASS(i)
                    accelz = -force*divz/MASS(i)
                    typ(i)%mol(i1)%accx = typ(i)%mol(i1)%accx + accelx
                    typ(i)%mol(i2)%accx = typ(i)%mol(i2)%accx - accelx
                    typ(i)%mol(i1)%accy = typ(i)%mol(i1)%accy + accely
                    typ(i)%mol(i2)%accy = typ(i)%mol(i2)%accy - accely
                    typ(i)%mol(i1)%accz = typ(i)%mol(i1)%accz + accelz
                    typ(i)%mol(i2)%accz = typ(i)%mol(i2)%accz - accelz
                    typ(i)%mol(i1)%poten = typ(i)%mol(i1)%poten + ppp*0.500d0
                    typ(i)%mol(i2)%poten = typ(i)%mol(i2)%poten + ppp*0.500d0
                end do
            end do
        end do

        ! 異なる分子同士の影響  ! Ar-Ptの処理    配列は1が上Pt, 2がAr, 3がしたPt
        do i = 1, TYPMOL
            if(i == 2) then
                cycle
            end if
            do i1 = 1, nummol(i)       ! Pt
                do i2 = 1, nummol(2)    ! Ar
                    divx = typ(i)%mol(i1)%posx - typ(2)%mol(i2)%posx
                    divy = typ(i)%mol(i1)%posy - typ(2)%mol(i2)%posy
                    divz = typ(i)%mol(i1)%posz - typ(2)%mol(i2)%posz

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

                    divx = divx / SIG(4)
                    divy = divy / SIG(4)
                    divz = divz / SIG(4)

                    if (divx < -CUTOFF .or. divx > CUTOFF) then
                        cycle
                    endif

                    if (divy < -CUTOFF .or. divy > CUTOFF) then
                        cycle
                    endif

                    if (divz < -CUTOFF .or. divz > CUTOFF) then
                        cycle
                    endif

                    dit2 = divx*divx + divy*divy + divz*divz
                    dist = dsqrt(dit2)

                    if(dist > CUTOFF) then
                        cycle
                    endif

                    dit4   = dit2*dit2
                    dit6   = dit4*dit2
                    dit8   = dit4*dit4
                    dit12  = dit6*dit6
                    dit14  = dit8*dit6
                    ppp    = angCon*4.00d0*EPS(4)*(1.00d0/dit12-1.00d0/dit6)    ! 異分子間ではangCon(接触角)を忘れずに
                    force  = forCoef(4)*(-2.00d0/dit14+1.00d0/dit8)
                    accelx = -force*divx/MASS(i)
                    typ(i)%mol(i1)%accx = typ(i)%mol(i1)%accx + accelx
                    accelx = -force*divx/MASS(2)
                    typ(2)%mol(i2)%accx = typ(2)%mol(i2)%accx - accelx
                    accely = -force*divy/MASS(i)
                    typ(i)%mol(i1)%accy = typ(i)%mol(i1)%accy + accely
                    accely = -force*divy/MASS(2)
                    typ(2)%mol(i2)%accy = typ(2)%mol(i2)%accy - accely
                    accelz = -force*divz/MASS(i)
                    typ(i)%mol(i1)%accz = typ(i)%mol(i1)%accz + accelz
                    accelz = -force*divz/MASS(2)
                    typ(2)%mol(i2)%accz = typ(2)%mol(i2)%accz - accelz
                    typ(i)%mol(i1)%poten = typ(i)%mol(i1)%poten + ppp*0.500d0
                    typ(2)%mol(i2)%poten = typ(2)%mol(i2)%poten + ppp*0.500d0
                end do
            end do
        end do

        ! Arの計算      
        do i = 1, nummol(2)
            vxene = typ(2)%mol(i)%velx + typ(2)%mol(i)%accx*0.500d0*dt
            vyene = typ(2)%mol(i)%vely + typ(2)%mol(i)%accy*0.500d0*dt
            vzene = typ(2)%mol(i)%velz + typ(2)%mol(i)%accz*0.500d0*dt
            vene = vxene*vxene + vyene*vyene + vzene*vzene
            typ(2)%mol(i)%kinet = 0.500d0*MASS(2)*vene
        end do

        do i = 1, nummol(2)
            typ(2)%mol(i)%velx = typ(2)%mol(i)%velx + typ(2)%mol(i)%accx*dt
            typ(2)%mol(i)%vely = typ(2)%mol(i)%vely + typ(2)%mol(i)%accy*dt
            typ(2)%mol(i)%velz = typ(2)%mol(i)%velz + typ(2)%mol(i)%accz*dt
            typ(2)%mol(i)%posx = typ(2)%mol(i)%posx + typ(2)%mol(i)%velx*dt
            typ(2)%mol(i)%posy = typ(2)%mol(i)%posy + typ(2)%mol(i)%vely*dt
            typ(2)%mol(i)%posz = typ(2)%mol(i)%posz + typ(2)%mol(i)%velz*dt
        end do

        ! Ptの計算
        do j = 1, TYPMOL
            if(j == 2) then     ! Arの場合を除外
                cycle
            end if

            ! 運動エネルギー計算
            do i = 1, nummol(j)
                vxene = typ(j)%mol(i)%velx + typ(j)%mol(i)%accx*0.500d0*dt
                vyene = typ(j)%mol(i)%vely + typ(j)%mol(i)%accy*0.500d0*dt
                vzene = typ(j)%mol(i)%velz + typ(j)%mol(i)%accz*0.500d0*dt
                vene = vxene**2 + vyene**2 + vzene**2
                typ(j)%mol(i)%kinet = 0.500d0*MASS(j)*vene
            end do

            ! 固定層
            do i = 1, int(nummol(j)/numz(j))
                typ(j)%mol(i)%velx = 0.0000d0
                typ(j)%mol(i)%vely = 0.0000d0
                typ(j)%mol(i)%velz = 0.0000d0
            end do

            ! その他の層
            do i = int(nummol(j)/numz(j)) + 1, int(nummol(j))
                typ(j)%mol(i)%velx = typ(j)%mol(i)%velx + dt * typ(j)%mol(i)%accx
                typ(j)%mol(i)%vely = typ(j)%mol(i)%vely + dt * typ(j)%mol(i)%accy
                typ(j)%mol(i)%velz = typ(j)%mol(i)%velz + dt * typ(j)%mol(i)%accz
                typ(j)%mol(i)%posx = typ(j)%mol(i)%posx + dt * typ(j)%mol(i)%velx
                typ(j)%mol(i)%posy = typ(j)%mol(i)%posy + dt * typ(j)%mol(i)%vely
                typ(j)%mol(i)%posz = typ(j)%mol(i)%posz + dt * typ(j)%mol(i)%velz
            end do
        end do
    end subroutine calcu

    subroutine bound ! 境界条件の付与
        use variable, only: xsyul, ysyul, zsyul
        use parameters
        use molecules_struct
        implicit none
        integer :: i, j
        do j = 1, TYPMOL
            do i = 1, nummol(j)
                if(typ(j)%mol(i)%posx < 0.00d0) then
                    typ(j)%mol(i)%posx = typ(j)%mol(i)%posx + xsyul
                else if(typ(j)%mol(i)%posx > xsyul) then
                    typ(j)%mol(i)%posx = typ(j)%mol(i)%posx - xsyul
                endif

                if(typ(j)%mol(i)%posy < 0.00d0) then
                    typ(j)%mol(i)%posy = typ(j)%mol(i)%posy + ysyul
                else if(typ(j)%mol(i)%posy > ysyul) then
                    typ(j)%mol(i)%posy = typ(j)%mol(i)%posy - ysyul
                endif

                ! z方向は周期境界の補正を行わない
                ! if(typ(j)%mol(i)%posz < 0.00d0) then
                !     typ(j)%mol(i)%posz = typ(j)%mol(i)%posz + zsyul
                ! else if(typ(j)%mol(i)%posz > zsyul) then
                !     typ(j)%mol(i)%posz = typ(j)%mol(i)%posz - zsyul
                ! endif
            end do
        end do
    end subroutine bound

    subroutine record_pos_vel! データの出力１
        use variable
        use parameters
        use molecules_struct
        implicit none
        integer :: i, j

        do i = 1, nummol(1)
            ! posit_PtUp.dat
            write(10, '(I6, 3E15.7)') i, typ(1)%mol(i)%posx, typ(1)%mol(i)%posy, typ(1)%mol(i)%posz
            ! veloc_PtUp.dat
            write(20, '(I6, 3E15.7)') i, typ(1)%mol(i)%velx, typ(1)%mol(i)%vely, typ(1)%mol(i)%velz
        end do
        do i = 1, nummol(2)
            ! posit_Ar.dat
            write(11, '(I6, 3E15.7)') i, typ(2)%mol(i)%posx, typ(2)%mol(i)%posy, typ(2)%mol(i)%posz
            ! veloc_Ar.dat
            write(21, '(I6, 3E15.7)') i, typ(2)%mol(i)%velx, typ(2)%mol(i)%vely, typ(2)%mol(i)%velz
        end do
        do i = 1, nummol(3)
            ! posit_PtDw.dat
            write(12, '(I6, 3E15.7)') i, typ(3)%mol(i)%posx, typ(3)%mol(i)%posy, typ(3)%mol(i)%posz
            ! veloc_PtDw.dat
            write(22, '(I6, 3E15.7)') i, typ(3)%mol(i)%velx, typ(3)%mol(i)%vely, typ(3)%mol(i)%velz
        end do

        ! 可視化用
        do j = 1, TYPMOL
            do i = 1, nummol(j)
                ! pos.dat
                write(15, '(3E15.7)') typ(j)%mol(i)%posx, typ(j)%mol(i)%posy, typ(j)%mol(i)%posz
            end do
        end do
    end subroutine record_pos_vel

    subroutine record_energy_temp ! データの出力２
        use variable
        use parameters
        use molecules_struct
        implicit none
        double precision, dimension(TYPMOL) :: totEne, totPot, totKin, temp
        double precision :: allEne, allPot, allKin
        integer :: i, j


        allEne = 0.000d0
        allPot = 0.000d0
        allKin = 0.000d0
        ! エネルギーの合計計算
        do j = 1, TYPMOL
            totEne(j) = 0.000d0
            totPot(j) = 0.000d0
            totKin(j) = 0.000d0
            do i = 1, nummol(j)
                totPot(j) = totPot(j) + typ(j)%mol(i)%poten
            end do
            if(j == 2) then
                do i = 1, nummol(j)     ! Ar
                    totKin(j) = totKin(j) + typ(j)%mol(i)%kinet
                end do
            else
                do i = int(nummol(j)/numz(j)) + 1, nummol(j)    ! Ptの場合は固定層を除く
                    totKin(j) = totKin(j) + typ(j)%mol(i)%kinet
                end do
            end if

            ! エネルギーを大きな数で割る処理（正規化や単位変換のため）
            totPot(j) = totPot(j) * 1.000d-16
            totKin(j) = totKin(j) * 1.000d-16
            totEne(j) = totPot(j) + totKin(j)

            if(j == 2) then
                ! Arの温度計算
                temp(j) = 2.0d0 * totKin(j) / (3.0d0 * dble(nummol(j)) * BOLTZ)
            else
                ! Ptの温度計算
                temp(j) = 2.0d0 * totKin(j) / (3.0d0 * dble(nummol(j) *2/numz(j)) * BOLTZ)
            end if

            allEne = allEne + totEne(j)
            allPot = allPot + totPot(j)
            allKin = allKin + totKin(j)
        end do

        write(30, '(I6, 4E15.7)') (stpNow+99)*int(dt), totEne(1), totPot(1), totKin(1)  ! energy_PtUp.dat
        write(31, '(I6, 4E15.7)') (stpNow+99)*int(dt), totEne(2), totPot(2), totKin(2)  ! energy_Ar.dat
        write(32, '(I6, 4E15.7)') (stpNow+99)*int(dt), totEne(3), totPot(3), totKin(3)  ! energy_PtDw.dat
        write(35, '(I6, 4E15.7)') (stpNow+99)*int(dt), allEne, allPot, allKin           ! energy_all.dat
        write(40, '(I6, 4E15.7)') (stpNow+99)*int(dt), temp(1), temp(2), temp(3)        ! tempe.dat
    end subroutine record_energy_temp

    subroutine record_finpos_vel
        use variable
        use parameters
        use molecules_struct
        implicit none
        integer :: i, j
        do j = 1, TYPMOL
            do i = 1, nummol(j)
                ! syuuki.dat
                write(50, '(I6, 6E15.7)') & 
                i, typ(j)%mol(i)%posx, typ(j)%mol(i)%posy, typ(j)%mol(i)%posz, &
                   typ(j)%mol(i)%velx, typ(j)%mol(i)%vely, typ(j)%mol(i)%velz
            end do
        end do
    end subroutine record_finpos_vel
end program main