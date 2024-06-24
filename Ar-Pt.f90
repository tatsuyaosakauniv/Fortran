!-----------------------------------------------------------------------
! example.f90
! ナノ粒子の初期状態を分子動力学法で作成する．
! 必要ファイル: random0.dat, random1.dat
!-----------------------------------------------------------------------

! パラメータ（もともとparate.datにあったもの）
module parameters
    implicit none
!--------------- よく変更する---------------------!
    ! 計算時間
    double precision, parameter :: timeScaling = 0.005d0 ! スケーリング時間[ns]
    double precision, parameter :: timeRelax   = 0.005d0 ! 緩和計算時間[ns]
    double precision, parameter :: timeMeasure = 0.010d0 ! データ計測時間[ns]
    double precision, parameter :: timeSum = timeScaling + timeRelax + timeMeasure ! 全計算時間
    
    double precision, parameter :: dt = 1.00d0 ! 無次元時間ステップ(無次元，有次元の場合fs)
    double precision, parameter :: tempAr = 175d0 ! 系内（目標）温度  K
    double precision, parameter :: angCon = 0.05d0 ! 接触角

    ! ステップ
    integer, parameter :: stpScaling = int(timeScaling/dt*1.0d+6) ! スケーリングステップ
    integer, parameter :: stpRelax   = int(timeRelax  /dt*1.0d+6) ! 緩和計算ステップ
    integer, parameter :: stpMeasure = int(timeMeasure/dt*1.0d+6) ! データ計測ステップ
    integer, parameter :: stpMax = int(timeSum/dt*1.0d+6) ! 最大ステップ数

    ! 分子の数，配置
    integer, parameter :: TYPMOL = 3 ! 分子の種類数
    integer, parameter :: COMP = 3 ! 相互作用の組み合わせ数 = TYPMOLC2
    integer, parameter :: numx(TYPMOL) = [10, 6, 10]
    integer, parameter :: numy(TYPMOL) = [ 5, 3,  5]
    integer, parameter :: numz(TYPMOL) = [ 4, 15, 4]
    integer, parameter :: nummol(TYPMOL) = [numx(1)*numy(1)*numz(1), numx(2)*numy(2)*numz(2), numx(3)*numy(3)*numz(3)] ! 各分子の数
    
    ! 境界条件
    double precision, parameter :: STDIST(TYPMOL) = [3.92d0, 6.0d0, 3.92d0] ! 格子定数(無次元)[Å]
    double precision, parameter :: xsyul0 = STDIST(1) * numy(1) ! x方向の周期境界長さ(無次元)
    double precision, parameter :: ysyul0 = STDIST(1) * numy(1) ! y方向の周期境界長さ(無次元）
    double precision, parameter :: zsyul0 = (STDIST(1)*numz(1)+STDIST(2)*numz(2)+STDIST(3)*numz(3))*0.5d0 ! z方向の周期境界長さ(無次元）
    double precision, parameter :: syul0(3) = [xsyul0, ysyul0, zsyul0]
    
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

    ! Langevin法
    double precision, parameter :: tempLanPt(TYPMOL) = [200d0, 0d0, 150d0] ! Langevin法を用いるPtの温度  真ん中は使わない
    double precision, parameter :: DIRAC = 1.054571817d-34 ! ディラック定数 [J･s]
    double precision, parameter :: DEBTMP = 240d0 ! Debye温度 [K]
    double precision, parameter :: OMEGA = 3.14212728482d+13 ! BOLTZ * DEBTMP / DIRAC * 1.000d-11 ! Debye定数 (有次元)
    double precision, parameter :: DAMP =  5.32967075080d-12 ! MASS(1) * PI * OMEGA / 6.000d0 ! ダンパーの減衰係数 (有次元)
    ! 熱流束
    double precision, parameter :: areaPt = STDIST(1)*STDIST(1)*int(numx(1)*0.5)*numy(1)
    !double precision, 

end module parameters

! 変数
module variable
    use parameters
    implicit none
    integer :: stpNow ! 現在のステップ数

    ! 分子間力
    double precision :: forCoef(4) = [0.0d0, 0.0d0, 0.0d0, 0.0d0] ! 力の計算の係数

    ! カットオフ
    double precision :: cutof(3) ! ポテンシャルのカットオフ長さx,y,z方向
    double precision :: syul(3) ! ポテンシャルのカットオフ長さx,y,z方向，x,y,z方向の周期長さ

    ! Langevin法
    double precision :: rndForce(nummol(1), 3, TYPMOL)   ! ランダム力用
    double precision :: dmpForce(nummol(1), 3, TYPMOL)   ! ダンパー力用
    double precision :: stddev ! 標準偏差
    double precision :: rnd_2
    logical :: isOdd = .true. ! 乱数のsinとcosを交互に出すためのフラグ

    ! 熱流束
    double precision :: interForce(nummol(2), 3) ! 熱流束を計算するための相互作用力
    double precision :: heatPhantom(TYPMOL) ! Phantom層からの熱輸送量
    double precision :: heatSl_Lq(TYPMOL) ! 固液界面での熱輸送量
    double precision :: fluxPt ! 熱流束
    double precision :: tempLayer(numz(1),TYPMOL) ! Arは使わない
    ! double precision :: tempLayerLw(numz(3))

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

        stddev_ = dsqrt(2.000d0 * BOLTZ * T_ * DAMP / dt * 1.000d+33)  ! 無次元 stddev -> 10^9

    end function getStddev
end module variable

module molecules_struct
    use parameters
    implicit none

    !構造体の定義
    type :: mol_info
        double precision :: pos(3)
        double precision :: vel(3)
        double precision :: acc(3)
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
    use parameters
    use variable, only: stpNow
    use molecules_struct
    use forPVwin
    implicit none
    integer :: i, j

    ! 配列初期化
    allocate(typ(1)%mol(nummol(1)))
    allocate(typ(2)%mol(nummol(2)))
    allocate(typ(3)%mol(nummol(3)))

    ! 読み込み用乱数ファイル -> call random_numberを使えば良い？
        !open(1,file='random0.dat', status='old')
        !open(2,file='random1.dat', status='old')

        !open(6,*) は使用できない
    ! aaa
        open(10,file='出力ファイル/posit_PtUp.dat')
        open(11,file='出力ファイル/posit_Ar.dat')
        open(12,file='出力ファイル/posit_PtLw.dat')
    ! 可視化用のpvch.fを移植 
        open(15,file='出力ファイル/pos.dat')
    ! 各分子の速度データの出力
        open(20,file='出力ファイル/veloc_PtUp.dat')
        open(21,file='出力ファイル/veloc_Ar.dat')
        open(22,file='出力ファイル/veloc_PtLw.dat')
    ! 系のエネルギーデータの出力
        open(30,file='出力ファイル/energy_PtUp.dat')
        open(31,file='出力ファイル/energy_Ar.dat')
        open(32,file='出力ファイル/energy_PtLw.dat')
        open(35,file='出力ファイル/energy_all.dat')
    ! 系の温度データの出力
        open(40,file='出力ファイル/tempe.dat')
        open(41,file='出力ファイル/tempe_PtUp_Layer.dat')
        open(42,file='出力ファイル/tempe_PtLw_Layer.dat')
    ! 系の周期長さの出力
        open(50,file='出力ファイル/syuuki.dat')
    ! 熱流束のデータ
        open(60,file='出力ファイル/heatflux.dat')
        
        open(70,file='出力ファイル/force.dat')

    ! 各分子の最終位置データの出力
        open(80,file='出力ファイル/finpos.dat')
    !　分子の色
        open(90,file='出力ファイル/mask.dat')

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
    
    !write(6,*) OMEGA, DAMP
    ! ターミナルに表示
    write(6,*) ''
    write(6, '(A17, F7.4, A2)') 'Scaling Time: ', timeScaling, 'ns'
    write(6, '(A17, F7.4, A2)') 'Relaxation Time: ', timeRelax, 'ns'
    write(6, '(A17, F7.4, A2)') 'Measure Time: ', timeMeasure, 'ns'
    write(6,*) ''
    write(6,'(A16)') '- Scaling Step -'
    write(6,*) ''
    
    
    stpNow = 0

    call seting ! 各分子の初期位置，初期速度などの設定

    do i = 1, stpMax
        stpNow = i

        ! ターミナルに表示
        if(stpNow == stpScaling) then
            write(6,*) ''
            write(6,'(A19)') '- Relaxation Step -'
            write(6,*) ''
        end if
        if(stpNow == stpScaling+stpRelax) then
            write(6,*) ''
            write(6,'(A16)') '- Measure Step -'
            write(6,*) ''
        end if

        ! ステップ数が500の倍数のとき
        if (mod(stpNow,500) == 0) then
            write(6,'(3X, I7, 1X, A3, I7)') stpNow, ' / ', stpMax
        endif

        ! スケーリング
        if (stpNow <= stpScaling .and. mod(stpNow,100) == 0) then
        	call scaling ! 系内の全分子の温度の補正
        endif

        call calcu ! 各分子に働く力，速度，位置の分子動力学計算
        call bound ! 境界条件の付与

        if (stpNow >= stpMeasure) then
            call calc_heatFlux
            if(mod(stpNow, 100) == 1) then
                call record_heatflux ! 熱流束を記録
            end if
        end if
        
        ! ステップ数が100の倍数+1のとき
        if(mod(stpNow, 100) == 1) then
          call record_pos_vel ! 位置，速度を記録
          call record_energy_temp ! エネルギー，温度を記録
        end if
    end do

    call record_finpos_vel ! 最終状態の分子の位置と速度を記録

    contains
    subroutine seting ! 各分子の初期位置，初期速度などの設定
        use parameters
        use variable
        use molecules_struct
        implicit none
        integer :: num, i, j, k
        double precision :: coord(3) = 0.0000D0 ! xyz座標
        double precision :: ofst(3)
        double precision :: ran, alpha, beta, cr
        double precision :: v(3)

        do i = 1, COMP
            forCoef(i) = 24.00d0*EPS(i)/SIG(i)
        end do
            forCoef(4) = angCon*24.00d0*EPS(4)/SIG(4)
        
        syul(:) = syul0(:)
        write(50,*)syul(1)
        write(50,*)syul(2)
        write(50,*)syul(3)
        write(15,*)syul0(1), syul0(2), syul0(3)
        write(15,*)ntime0, ndt

        do i = 1, TYPMOL
            cutof(:) = syul(:) - CUTOFF*SIG(i)
        end do

        num = 0

        !上段のPt配置
        ofst(1) = STDIST(1)*0.25d0
        ofst(2) = STDIST(1)*0.25d0
        ofst(3) = zsyul0 - STDIST(1)*0.25d0
        do k = 1,numz(1)
            coord(3) = ofst(3) - dble(k-1)*STDIST(1)*0.5d0
            do i = 1,numx(1)
                coord(1) = ofst(1) + dble(i-1)*STDIST(1)*0.5d0
                do j = 1,numy(1)
                    if(mod(k,2) == 0) then  !z偶数
                        if(mod(i,2) == 0) then
                            coord(2) = ofst(2) + dble(j-1)*STDIST(1)   !x偶数
                        else
                            coord(2) = ofst(2) + dble(j-1)*STDIST(1) + STDIST(1)*0.5d0    !x奇数
                        endif
                    else    !z奇数
                        if(mod(i,2) == 0) then
                            coord(2) = ofst(2) + dble(j-1)*STDIST(1) + STDIST(1)*0.5d0    !x偶数
                        else
                            coord(2) = ofst(2) + dble(j-1)*STDIST(1)   !x奇数
                        endif
                    endif
                    num = num + 1
                    typ(1)%mol(num)%pos(:) = coord(:)
                end do
            end do
        end do

        !中段のAr配置
        num = 0
        ofst(1) = syul0(1)*0.50d0 - STDIST(2)*(0.25d0*(numx(2)  -1))
        ofst(2) = syul0(2)*0.50d0 - STDIST(2)*(0.25d0*(numy(2)*2-1))
        ofst(3) = syul0(3)*0.50d0 - STDIST(2)*(0.25d0*(numz(2)  -1))
        do k = 1,numz(2)
            coord(3) = ofst(3) + dble(k-1)*STDIST(2)*0.5d0
            do i = 1,numx(2)
                coord(1) = ofst(1) + dble(i-1)*STDIST(2)*0.5d0
                do j = 1,numy(2)
                    if(mod(k,2) == 0) then  !z偶数
                        if(mod(i,2) == 0) then
                            coord(2) = ofst(2) + dble(j-1)*STDIST(2)   !x偶数
                        else
                            coord(2) = ofst(2) + dble(j-1)*STDIST(2) + STDIST(2)*0.5d0    !x奇数
                        endif
                    else    !z奇数
                        if(mod(i,2) == 0) then
                            coord(2) = ofst(2) + dble(j-1)*STDIST(2) + STDIST(2)*0.5d0    !x偶数
                        else
                            coord(2) = ofst(2) + dble(j-1)*STDIST(2)   !x奇数
                        endif
                    endif
                    num = num + 1
                    typ(2)%mol(num)%pos(:) = coord(:)
                end do
            end do
        end do

        !下段のPt配置
        num = 0
        ofst(1) = STDIST(3)*0.25d0
        ofst(2) = STDIST(3)*0.25d0
        ofst(3) = STDIST(3)*0.25d0
        do k = 1,numz(3)
            coord(3) = ofst(3) + dble(k-1)*STDIST(3)*0.5d0
            do i = 1,numx(3)
                coord(1) = ofst(1) + dble(i-1)*STDIST(3)*0.5d0
                do j = 1,numy(3)
                    if(mod(k,2) == 0) then  !z偶数
                        if(mod(i,2) == 0) then
                            coord(2) = ofst(2) + dble(j-1)*STDIST(3)   !x偶数
                        else
                            coord(2) = ofst(2) + dble(j-1)*STDIST(3) + STDIST(3)*0.5d0    !x奇数
                        endif
                    else    !z奇数
                        if(mod(i,2) == 0) then
                            coord(2) = ofst(2) + dble(j-1)*STDIST(3) + STDIST(3)*0.5d0    !x偶数
                        else
                            coord(2) = ofst(2) + dble(j-1)*STDIST(3)   !x奇数
                        endif
                    endif
                    num = num + 1
                    typ(3)%mol(num)%pos(:) = coord(:)
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
                v(1) = dsin(alpha)*dcos(beta)*cr
                v(2) = dsin(alpha)*dsin(beta)*cr
                v(3) = dcos(alpha)*cr
                typ(j)%mol(i)%vel(:) = v(:)
            end do
        end do

        do j = 1, TYPMOL
            do i = 1, int(nummol(j)/numz(j))
                if(j == 2) then
                    cycle
                else
                    typ(j)%mol(i)%vel(:) = 0.000d0
                end if
            end do
        end do
    end subroutine seting

    subroutine cortra ! 系内の全分子の並進速度の補正
        use parameters
        use variable
        implicit none
        double precision :: trv(3)
        integer :: i, j

        ! 速度ベクトルの成分の平均を計算
        do j = 1, TYPMOL
            trv(:) = 0.0d0
            do i = 1, nummol(j)
                trv(:) = trv(:) + typ(j)%mol(i)%vel(:)
            end do

            trv(:) = trv(:) / nummol(j)

            ! 速度ベクトルから平均を引いて中心補正
            do i = 1, nummol(j)
                typ(j)%mol(i)%vel(:) = typ(j)%mol(i)%vel(:) - trv(:)
            end do
        end do
    end subroutine cortra

    subroutine jyusin ! 系内の全分子の重心の補正
        use parameters
        use variable
        use molecules_struct
        implicit none
        double precision :: cms(3)
        double precision :: tcms(3)
        integer :: i, j

        cms(:) = syul0(:) * 0.500d0

        do j = 1, TYPMOL
            tcms(:) = 0.0000d0
            do i = 1, nummol(j)
                tcms(:) = tcms(:) + typ(j)%mol(i)%pos(:)
            end do

            tcms(:) = cms(:) - tcms(:)/dble(nummol(j))

            do i = 1, nummol(j)
                typ(j)%mol(i)%pos(:) = typ(j)%mol(i)%pos(:) + tcms(:)
            end do
        end do
    end subroutine jyusin

    subroutine scaling ! 系内の全分子の温度の補正
        use parameters
        use variable
        use molecules_struct
        implicit none
        double precision :: temptp, vel2, aimtem, aimnot, baiss
        integer :: i
        integer :: j = 2 ! Arのみ

        temptp = 0.000d0
        do i = 1, nummol(j)
            vel2 = typ(j)%mol(i)%vel(1)**2 + typ(j)%mol(i)%vel(2)**2 + typ(j)%mol(i)%vel(3)**2
            temptp = temptp + vel2
        end do
        temptp = temptp / nummol(j) * 1.000d-16
        aimtem = tempAr
        aimnot = 3.000d0 * BOLTZ * aimtem / MASS(j)
        baiss = dsqrt(aimnot / temptp)

        ! 速度ベクトルのスケーリング
        ! Arのみ
        do i = 1, nummol(j)
            typ(j)%mol(i)%vel(:) = typ(j)%mol(i)%vel(:) * baiss
        end do
    end subroutine scaling

    subroutine calcu ! 各分子に働く力，速度，位置の分子動力学計算
        use parameters
        use variable
        use molecules_struct
        implicit none
        integer :: i, j, k, i1, i2
        double precision :: div(3), dist
        double precision :: dit2, dit4, dit6, dit8, dit12, dit14
        double precision :: ppp, force, forVec(3)
        double precision :: vene(3), sumvene
        double precision :: rnd

        do j = 1, TYPMOL
            do i = 1, nummol(j)
                typ(j)%mol(i)%acc(:) = 0.0000d0
                typ(j)%mol(i)%poten  = 0.0000d0
                typ(j)%mol(i)%kinet  = 0.0000d0
            end do
        end do
        interForce(:,:) = 0.000d0

        ! PtのPhantom層はダンパー力とランダム力を付与
        do j = 1, TYPMOL
            if(j == 2) then
                cycle
            end if

            do i = int(nummol(j)/numz(j)) + 1, 2*int(nummol(j)/numz(j)) ! Phantom層のみ
                do k = 1, 3
                    rnd = Random() 
                    ! ランダム力
                    rndForce(i,k,j) = rnd * getStddev(tempLanPt(j)) * 1.000d-9 ! 標準偏差の有次元化
                    ! ダンパー力
                    dmpForce(i,k,j) = - DAMP * typ(j)%mol(i)%vel(k) * 1.000d+5 ! 速度の有次元化
                end do
            end do

            ! ランダム力とダンパー力を追加
            do i = int(nummol(j)/numz(j)) + 1, 2*int(nummol(j)/numz(j))         ! 加速度の無次元化 10^-20
                typ(j)%mol(i)%acc(:) = (rndForce(i,:,j)*1.0d+9 + dmpForce(i,:,j)*1.0d+9) / MASS(j)*1.000d-3
            end do
        end do

        ! 分子間の相互作用力 → ポテンシャルエネルギー
        ! 同じ分子同士の影響
        do j = 1, TYPMOL
            do i1 = 1, nummol(j)
                do i2 = i1+1, nummol(j)
                    div(:) = typ(j)%mol(i1)%pos(:) - typ(j)%mol(i2)%pos(:)
    
                    ! カットオフ
                    do k = 1, 3
                        if (div(k) < -cutof(k)) then
                            div(k) = div(k) + syul(k)
                        else if(div(k) > cutof(k)) then
                            div(k) = div(k) - syul(k)
                        endif
        
                        div(k) = div(k) / SIG(j)

                        if (abs(div(k)) > CUTOFF) then
                            cycle
                        endif
                    end do

                    dit2 = div(1)**2 + div(2)**2 + div(3)**2
                    dist = dsqrt(dit2)
    
                    if(dist > CUTOFF) then
                        cycle
                    endif
    
                    dit4   = dit2*dit2
                    dit6   = dit4*dit2
                    dit8   = dit4*dit4
                    dit12  = dit6*dit6
                    dit14  = dit8*dit6
                    ppp    = 4.00d0*EPS(j)*(1.00d0/dit12-1.00d0/dit6)
                    typ(j)%mol(i1)%poten = typ(j)%mol(i1)%poten + ppp*0.500d0
                    typ(j)%mol(i2)%poten = typ(j)%mol(i2)%poten + ppp*0.500d0

                    force  = forCoef(j)*(-2.00d0/dit14+1.00d0/dit8)
                    forVec(:) = -force*div(:)
                    typ(j)%mol(i1)%acc(:) = typ(j)%mol(i1)%acc(:) + forVec(:)/MASS(j)
                    typ(j)%mol(i2)%acc(:) = typ(j)%mol(i2)%acc(:) - forVec(:)/MASS(j)
                end do
            end do
        end do

        ! 異なる分子同士の影響  ! Ar-Ptの処理    配列は1が上Pt, 2がAr, 3がしたPt
        do j = 1, TYPMOL
            if(j == 2) then
                cycle
            end if
            do i1 = 1, nummol(j)       ! Pt
                do i2 = 1, nummol(2)    ! Ar
                    div(:) = typ(j)%mol(i1)%pos(:) - typ(2)%mol(i2)%pos(:)

                    do k = 1, 3
                        if (div(k) < -cutof(k)) then
                            div(k) = div(k) + syul(k)
                        else if(div(k) > cutof(k)) then
                            div(k) = div(k) - syul(k)
                        endif
        
                        div(k) = div(k) / SIG(j)

                        if (abs(div(k)) > CUTOFF) then
                            cycle
                        endif
                    end do

                    dit2 = div(1)**2 + div(2)**2 + div(3)**2
                    dist = dsqrt(dit2)

                    if(dist > CUTOFF) then
                        cycle
                    endif

                    dit4   = dit2*dit2
                    dit6   = dit4*dit2
                    dit8   = dit4*dit4
                    dit12  = dit6*dit6
                    dit14  = dit8*dit6
                    ppp    = angCon*4.00d0*EPS(4)*(1.00d0/dit12-1.00d0/dit6) ! 異分子間ではangCon(接触角)を忘れずに
                    typ(j)%mol(i1)%poten = typ(j)%mol(i1)%poten + ppp*0.500d0
                    typ(2)%mol(i2)%poten = typ(2)%mol(i2)%poten + ppp*0.500d0

                    force  = forCoef(4)*(-2.00d0/dit14+1.00d0/dit8)
                    do k = 1, 3
                        forVec(k) = -force*div(k)
                        interForce(i1,k) = interForce(i1,k) + forVec(k) ! 無次元なことに注意
                    end do
                    typ(j)%mol(i1)%acc(:) = typ(j)%mol(i1)%acc(:) + forVec(:)/MASS(j)
                    typ(2)%mol(i2)%acc(:) = typ(2)%mol(i2)%acc(:) - forVec(:)/MASS(2)
                end do
            end do
        end do

        ! Arの計算  
        ! 運動エネルギー計算    
        do i = 1, nummol(2)
            vene(:) = typ(2)%mol(i)%vel(:) + typ(2)%mol(i)%acc(:)*0.500d0*dt
            sumvene = vene(1)**2 + vene(2)**2 + vene(3)**2
            typ(2)%mol(i)%kinet = 0.500d0*MASS(2)*sumvene
        end do

        ! 数値積分 (蛙跳び法)
        do i = 1, nummol(2)
            typ(2)%mol(i)%vel(:) = typ(2)%mol(i)%vel(:) + typ(2)%mol(i)%acc(:)*dt
            typ(2)%mol(i)%pos(:) = typ(2)%mol(i)%pos(:) + typ(2)%mol(i)%vel(:)*dt
        end do

        ! Ptの計算
        do j = 1, TYPMOL
            if(j == 2) then ! Arの場合を除外
                cycle
            end if

            ! 運動エネルギー計算
            do i = 1, nummol(j)
                vene(:) = typ(j)%mol(i)%vel(:) + typ(j)%mol(i)%acc(:)*0.500d0*dt
                sumvene = vene(1)**2 + vene(2)**2 + vene(3)**2
                typ(j)%mol(i)%kinet = 0.500d0*MASS(j)*sumvene
            end do

            ! 固定層
            do i = 1, int(nummol(j)/numz(j))
                typ(j)%mol(i)%vel(:) = 0.0000d0
            end do

            ! その他の層
            do i = int(nummol(j)/numz(j)) + 1, int(nummol(j))
                typ(j)%mol(i)%vel(:) = typ(j)%mol(i)%vel(:) + dt * typ(j)%mol(i)%acc(:)
                typ(j)%mol(i)%pos(:) = typ(j)%mol(i)%pos(:) + dt * typ(j)%mol(i)%vel(:)
            end do
        end do
    end subroutine calcu

    subroutine bound ! 境界条件の付与
        use parameters
        use variable
        use molecules_struct
        implicit none
        integer :: i, j
        do j = 1, TYPMOL
            do i = 1, nummol(j)
                if(typ(j)%mol(i)%pos(1) < 0.00d0) then
                    typ(j)%mol(i)%pos(1) = typ(j)%mol(i)%pos(1) + syul(1)
                else if(typ(j)%mol(i)%pos(1) > syul(1)) then
                    typ(j)%mol(i)%pos(1) = typ(j)%mol(i)%pos(1) - syul(1)
                endif

                if(typ(j)%mol(i)%pos(2) < 0.00d0) then
                    typ(j)%mol(i)%pos(2) = typ(j)%mol(i)%pos(2) + syul(2)
                else if(typ(j)%mol(i)%pos(2) > syul(2)) then
                    typ(j)%mol(i)%pos(2) = typ(j)%mol(i)%pos(2) - syul(2)
                endif

                ! z方向は周期境界の補正を行わない
                ! if(typ(j)%mol(i)%pos(3) < 0.00d0) then
                !     typ(j)%mol(i)%pos(3) = typ(j)%mol(i)%pos(3) + syul(3)
                ! else if(typ(j)%mol(i)%pos(3) > syul(3)) then
                !     typ(j)%mol(i)%pos(3) = typ(j)%mol(i)%pos(3) - syul(3)
                ! endif
            end do
        end do
    end subroutine bound

    subroutine calc_heatFlux
        use parameters
        use variable
        use molecules_struct
        implicit none
        integer :: i, j, k

        heatPhantom(:) = 0.000d0
        heatSl_Lq = 0.000d0

        do j = 1, TYPMOL
            if (j == 2) then
                cycle
            else
                do k = 1, 3
                    do i = int(nummol(j)/numz(j)) + 1, 2*int(nummol(j)/numz(j)) ! Phantom層              ! 速さの有次元化 10^5
                        heatPhantom(j) = heatPhantom(j) + (rndForce(i,k,j) + dmpForce(i,k,j))*typ(j)%mol(i)%vel(k)*1.000d+5
                    end do
                end do
                        heatPhantom(j) = heatPhantom(j) / (areaPt * 1.000d-20) ! 面積の有次元化 10^-20
            end if
        end do

        do j = 1, TYPMOL
            if (j == 2) then
                cycle
            else
                do k = 1, 3
                    !do i = (numz(j)-1)*int(nummol(j)/numz(j)) + 1, nummol(j) ! 固液界面層
                    do i = 1, nummol(2)
                        heatSl_Lq(j) = heatSl_Lq(j) + interForce(i,k)*1.000d-10*typ(j)%mol(i)%vel(k)*1.000d+5
                    end do
                end do
                        heatSl_Lq(j) = heatSl_Lq(j) / (areaPt * 1.000d-20)
            end if
        end do
    
    end subroutine calc_heatFlux

    subroutine record_pos_vel ! 位置，速度を記録
        use parameters
        use variable
        use molecules_struct
        implicit none
        integer :: i, j

        do i = 1, nummol(1)
            ! posit_PtUp.dat
            write(10, '(I6, 3E15.7)') i, typ(1)%mol(i)%pos(1), typ(1)%mol(i)%pos(2), typ(1)%mol(i)%pos(3)
            ! veloc_PtUp.dat
            write(20, '(I6, 3E15.7)') i, typ(1)%mol(i)%vel(1), typ(1)%mol(i)%vel(2), typ(1)%mol(i)%vel(3)
        end do
        do i = 1, nummol(2)
            ! posit_Ar.dat
            write(11, '(I6, 3E15.7)') i, typ(2)%mol(i)%pos(1), typ(2)%mol(i)%pos(2), typ(2)%mol(i)%pos(3)
            ! veloc_Ar.dat
            write(21, '(I6, 3E15.7)') i, typ(2)%mol(i)%vel(1), typ(2)%mol(i)%vel(2), typ(2)%mol(i)%vel(3)
        end do
        do i = 1, nummol(3)
            ! posit_PtLw.dat
            write(12, '(I6, 3E15.7)') i, typ(3)%mol(i)%pos(1), typ(3)%mol(i)%pos(2), typ(3)%mol(i)%pos(3)
            ! veloc_PtLw.dat
            write(22, '(I6, 3E15.7)') i, typ(3)%mol(i)%vel(1), typ(3)%mol(i)%vel(2), typ(3)%mol(i)%vel(3)
        end do

        do i = int(nummol(1)/numz(1)) + 1, int(nummol(1) *2/numz(1))
            write(70, '(9E15.7)') rndForce(i,1,1), rndForce(i,2,1), rndForce(i,3,1), dmpForce(i,1,1), dmpForce(i,2,1), dmpForce(i,3,1), interForce(i,1), interForce(i,2), interForce(i,3)
        end do

        ! 可視化用
        do j = 1, TYPMOL
            do i = 1, nummol(j)
                ! pos.dat
                write(15, '(3E15.7)') typ(j)%mol(i)%pos(1), typ(j)%mol(i)%pos(2), typ(j)%mol(i)%pos(3)
            end do
        end do
    end subroutine record_pos_vel

    subroutine record_energy_temp ! エネルギー，温度を記録
        use parameters
        use variable
        use molecules_struct
        implicit none
        double precision, dimension(TYPMOL) :: totEne, totPot, totKin, temp
        double precision :: allEne, allPot, allKin, kintmp
        integer :: i, j, k

        allEne = 0.000d0
        allPot = 0.000d0
        allKin = 0.000d0
        totEne(:) = 0.000d0
        totPot(:) = 0.000d0
        totKin(:) = 0.000d0
        temp(:) = 0.000d0

        do j = 1, TYPMOL
            ! ポテンシャル
            do i = 1, nummol(j)
                totPot(j) = totPot(j) + typ(j)%mol(i)%poten
            end do
            totPot(j) = totPot(j) * 1.000d-16

            ! 運動エネルギー
            if (j == 2) then
                ! Ar
                do i = 1, nummol(j)
                    totKin(j) = totKin(j) + typ(j)%mol(i)%kinet
                end do
                totKin(j) = totKin(j) * 1.000d-16
                temp(j) = 2.0d0 * totKin(j) / (3.0d0 * dble(nummol(j)) * BOLTZ)
            else
                ! Pt
                do k = 1, numz(j) ! Ptの層の数
                    kintmp = 0.000d0
                    do i = (k-1)*int(nummol(j)/numz(j)) + 1, k*int(nummol(j)/numz(j))
                        kintmp = kintmp + typ(j)%mol(i)%kinet
                    end do

                    kintmp = kintmp * 1.000d-16
                    tempLayer(k, j) = 2.0d0 * kintmp / (3.0d0 * dble(nummol(j)/numz(j)) * BOLTZ)

                    if (k /= 1) then ! 固定層は除く
                        temp(j) = temp(j) + tempLayer(k, j)
                    end if
                end do

                temp(j) = temp(j) / dble(numz(j)-1)

            end if

            totEne(j) = totPot(j) + totKin(j)
            allEne = allEne + totEne(j)
            allPot = allPot + totPot(j)
            allKin = allKin + totKin(j)

        end do

        write(30, '(I6, 4E15.7)') (stpNow+99)*int(dt), totEne(1), totPot(1), totKin(1)  ! energy_PtUp.dat
        write(31, '(I6, 4E15.7)') (stpNow+99)*int(dt), totEne(2), totPot(2), totKin(2)  ! energy_Ar.dat
        write(32, '(I6, 4E15.7)') (stpNow+99)*int(dt), totEne(3), totPot(3), totKin(3)  ! energy_PtLw.dat
        write(35, '(I6, 4E15.7)') (stpNow+99)*int(dt), allEne, allPot, allKin           ! energy_all.dat
        write(40, '(I6, 4E15.7)') (stpNow+99)*int(dt), temp(1), temp(2), temp(3)        ! tempe.dat

        !!!!!!!!!Pt層を増やすとき必ず変更すること!!!!!!!!!
        write(41, '(I6, 4E15.7)') (stpNow+99)*int(dt), tempLayer(1,1), tempLayer(2,1), tempLayer(3,1), tempLayer(4,1)
        write(42, '(I6, 4E15.7)') (stpNow+99)*int(dt), tempLayer(1,3), tempLayer(2,3), tempLayer(3,3), tempLayer(4,3)
    end subroutine record_energy_temp

    subroutine record_heatflux ! 熱流束を記録
        use parameters
        use variable
        use molecules_struct
        implicit none

        write(60,'(I6, 4E15.7)') (stpNow+99)*int(dt), heatPhantom(1), heatPhantom(3), heatSl_Lq(1), heatSl_Lq(3)
    end subroutine record_heatflux

    subroutine record_finpos_vel ! 最終状態の分子の位置と速度を記録
        use parameters
        use variable
        use molecules_struct
        implicit none
        integer :: i, j
        do j = 1, TYPMOL
            do i = 1, nummol(j)
                ! syuuki.dat
                write(50, '(I6, 6E15.7)') & 
                i, typ(j)%mol(i)%pos(1), typ(j)%mol(i)%pos(2), typ(j)%mol(i)%pos(3), &
                   typ(j)%mol(i)%vel(1), typ(j)%mol(i)%vel(2), typ(j)%mol(i)%vel(3)
            end do
        end do
    end subroutine record_finpos_vel
end program main