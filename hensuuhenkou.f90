! 変数名のルール
! num → 数字, 数 (number)
! mol → 分子 (molecule)
! temp → 温度 (temepature)
! time → 時間 (time)
! step → ステップ (step)
! top → 上壁面 (top layer)
! btm → 下壁面 (bottom layer)
! coef → 係数 (coefficient)
! heat → 熱量 (heat)
! flux → 熱流束 (heat flux)
! pot → ポテンシャルエネルギー (potential energy)
! kin → 運動エネルギー (kinetic energy)
! len → 長さ (length)


module enum
    integer, parameter :: X = 1, Y = 2, Z = 3, XYZ = 3 ! xyz方向
    integer, parameter :: Pt = 2, Ar = 3, PtAr = 3 ! 分子の種類
    integer, parameter :: TOP = 1, BTM = 2 ! Ptの上下壁面
end module enum

module parameters
    use enum
    implicit none
!--------------- よく変更する---------------------!
    ! 計算時間
    double precision, parameter :: time_scaling = 0.001d0 ! スケーリング時間[ns]
    double precision, parameter :: time_relax   = 0.001d0 ! 緩和計算時間[ns]
    double precision, parameter :: time_measure = 0.001d0 ! データ計測時間[ns]
    double precision, parameter :: time_all = time_scaling + time_relax + time_measure ! 全計算時間
    
    double precision, parameter :: dt = 1.00d0 ! 無次元時間ステップ(無次元，有次元の場合fs)
    double precision, parameter :: tau = 1.00d0 ! 測定間隔[fs]

    double precision, parameter :: temp_Ar = 110d0 ! 系内（目標）温度  K
    double precision, parameter :: INTER_STG = 0.10d0 ! 相互作用強さ (interaction strength)

    character(len=20) :: dir_name = 'hensuu'  ! ファイル名

    ! ステップ
    integer, parameter :: step_scaling = int(time_scaling/dt*1.0d+6) ! スケーリングステップ
    integer, parameter :: step_relax   = int(time_relax  /dt*1.0d+6) ! 緩和計算ステップ
    integer, parameter :: step_measure = int(time_measure/dt*1.0d+6) ! データ計測ステップ
    integer, parameter :: step_all = step_scaling + step_relax + step_measure ! 最大ステップ数

    ! 分子の数，配置
    integer, parameter :: NUM_X(PtAr) = [20,20,12] !
    integer, parameter :: NUM_Y(PtAr) = [10,10, 6] !
    integer, parameter :: NUM_Z(PtAr) = [ 6, 6,25] !
    integer, parameter :: NUM_MOL(PtAr) = [NUM_X(TOP)*NUM_Y(TOP)*NUM_Z(TOP), NUM_X(BTM)*NUM_Y(BTM)*NUM_Z(BTM), NUM_X(Ar)*NUM_Y(Ar)*NUM_Z(Ar)] ! 各分子の数
    integer, parameter :: NUM_MOL_ALL = NUM_MOL(TOP) + NUM_MOL(BTM) + NUM_MOL(Ar)
    
    ! 境界条件
    double precision, parameter :: STD_DIST(PtAr) = [3.92d0, 3.92d0, 5.0d0] ! 格子定数(無次元)[Å] (stndard distance)
    double precision, parameter :: THICK_Pt_TOP = STD_DIST(TOP)*(NUM_Z(TOP)*0.5d0+0.25d0) !
    double precision, parameter :: THICK_Pt_BTM = STD_DIST(BTM)*(NUM_Z(BTM)*0.5d0+0.25d0) !
    double precision, parameter :: THICK_Ar = STD_DIST(Ar)*NUM_Z(Ar)*0.5d0 !
    double precision, parameter :: THICK(PtAr) = [THICK_Pt_TOP, THICK_Pt_BTM, THICK_Ar] !
    double precision, parameter :: BND_LEN_X0 = STD_DIST(Pt) * NUM_Y(Pt) ! x方向の周期境界長さ(無次元) (boundary length)
    double precision, parameter :: BND_LEN_Y0 = STD_DIST(Pt) * NUM_Y(Pt) ! y方向の周期境界長さ(無次元）
    double precision, parameter :: BND_LEN_Z0 = 87.0d0 ! THICK(X)+THICK(Y)+THICK(Z) ! z方向の周期境界長さ(無次元）
    double precision, parameter :: BND_LEN0(XYZ) = [BND_LEN_X0, BND_LEN_Y0, BND_LEN_Z0]
    
    double precision, parameter :: CUTOFF = 3.300d0 ! カットオフ長さ/σ
    double precision, parameter :: AVOGA = 6.022d+23 ! アボガドロ数
    double precision, parameter :: BOLTZ = 1.3806662d-23 ! ボルツマン定数 [J/K]
    double precision, parameter :: PI = 3.141592654d0 ! 円周率

    ! 分子の質量
    !double precision, parameter :: bunsi(PtAr) = [195.084d-3, 195.084d-3, 39.950d-3] ! 分子の質量  kg/mol   
    double precision, parameter :: MASS(PtAr) = [32.395d0, 32.395d0, 6.6340d0] ! 分子の質量（無次元） * 10d-26 [kg/個]
    ! Lennard-Jonesパラメータ
    double precision, parameter :: SIG(4) = [2.540d0, 2.540d0, 3.400d0, 2.970d0]  ! σ(無次元) *1.0d-10
    double precision, parameter :: EPS(4) = [109.2d-5, 109.2d-5, 1.666d-5, 13.49d-5] ! ε(無次元) *1.0d-16

    ! 粒子登録法
    integer, parameter :: step_update = 40 ! 更新ステップ

    ! Langevin法
    double precision, parameter :: temp_Langevin(Pt) = [120d0, 100d0] ! Langevin法を用いるPtの温度  真ん中は使わない
    double precision, parameter :: DIRAC = 1.054571817d-34 ! ディラック定数 [J･s]
    double precision, parameter :: DEBTMP = 240d0 ! Debye温度 [K]
    double precision, parameter :: OMEGA = 3.14212728482d+13 ! BOLTZ * DEBTMP / DIRAC * 1.000d-11 ! Debye定数 (有次元)
    double precision, parameter :: DAMP =  5.32967075080d-12 ! MASS(1) * PI * OMEGA / 6.000d0 ! ダンパーの減衰係数 (有次元)
    ! 熱流束
    double precision, parameter :: AREA_Pt = STD_DIST(Pt)*STD_DIST(Pt)*int(NUM_X(Pt)*0.5)*NUM_Y(Pt) !
    
    integer, parameter :: NUM_DIV_Ar = 15 ! Arの温度分布の分割数

end module parameters

! 変数
module variable
    use enum
    use parameters
    implicit none
    integer :: step_now ! 現在のステップ数

    ! 分子間力
    double precision :: coef_force(4) = [0.0d0, 0.0d0, 0.0d0, 0.0d0] ! 力の計算の係数

    ! カットオフ
    double precision :: cutof(XYZ) ! ポテンシャルのカットオフ長さx,y,z方向
    double precision :: bnd_len(XYZ) ! ポテンシャルのカットオフ長さx,y,z方向，x,y,z方向の周期長さ

    ! Langevin法
    double precision :: force_rand(NUM_MOL(Pt),Pt,XYZ)   ! ランダム力用
    double precision :: force_damp(NUM_MOL(Pt),Pt,XYZ)   ! ダンパー力用
    double precision :: stddev ! 標準偏差
    double precision :: rnd_2
    logical :: isOdd = .true. ! 乱数のsinとcosを交互に出すためのフラグ

    ! 熱流束
    double precision :: force_inter(NUM_MOL(Pt),Pt,XYZ) ! 熱流束を計算するための相互作用力 z方向のみを使う
    double precision :: heat_phantom(Pt)! = 0.000d0 ! Phantom層からの熱輸送量
    double precision :: heat_interface(Pt)! = 0.000d0 ! 固液界面での熱輸送量
    double precision :: pressure(Pt) ! 圧力

    double precision :: temp_layer_Pt(NUM_Z(Pt),Pt) = 0.000d0 ! Ptの層ごとの温度
    double precision :: temp_layer_Ar(NUM_DIV_Ar) = 0.000d0 ! z方向に分割した領域内のArの温度
    double precision :: zdiv = THICK(Ar)/NUM_DIV_Ar ! Arの温度分布の分割距離

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

        stddev_ = dsqrt(2.000d0 * DAMP * BOLTZ * T_ / dt * 1.000d+33)  ! 無次元 stddev -> 有次元では10^9

    end function getStddev
end module variable

module molecules_struct
    use enum
    use parameters
    implicit none

    !構造体の定義
    type :: mol_info
        double precision :: pos(XYZ)
        double precision :: vel(XYZ), vmid(XYZ) ! vmidは運動エネルギー計算で必要になる v(t)
        double precision :: acc(XYZ)
        double precision :: pot, kin
    end type mol_info

    type mol_typ
        type(mol_info), allocatable :: mol(:)
    end type mol_typ

    type(mol_typ), dimension(PtAr) :: typ  ! 上Pt, 中Ar, 下Pt

end module molecules_struct

! module forPVwin
!     use enum
!     use parameters
!     implicit none
!     integer, parameter :: NUM_MOL_ALL = NUM_MOL(1) + NUM_MOL(2) + NUM_MOL(3)
!     integer, parameter :: moltype = 1
!     integer, parameter :: ndat = int(step_all/100)
!     integer, parameter :: ntime0 = 0
!     integer, parameter :: ndt = 1
! end module forPVwin

program main
    use enum
    use parameters
    use variable
    use molecules_struct
    ! use forPVwin
    implicit none
    double precision :: time
    integer :: i, j

    ! 日付と時刻
    ! character(len=8) :: date
    character(len=100) :: filepath

    filepath = '/home/kawaguchi/' // trim(dir_name)

    ! 配列初期化
    allocate(typ(TOP)%mol(NUM_MOL(TOP)))
    allocate(typ(BTM)%mol(NUM_MOL(BTM)))
    allocate(typ(Ar)%mol(NUM_MOL(Ar)))

        !open(6,*) は使用できない
    ! 各分子の位置データの出力
        open(10,file=trim(filepath) // '/posit_Pt_top.dat')
        open(11,file=trim(filepath) // '/posit_Pt_bottom.dat')
        open(12,file=trim(filepath) // '/posit_Ar.dat')
    ! 可視化用のpvch.fを移植 
        ! open(15,file=trim(filepath) // '/pos.dat')
        open(16,file=trim(filepath) // '/ovito.dat')
    ! 各分子の速度データの出力
        open(20,file=trim(filepath) // '/veloc_Pt_top.dat')
        open(21,file=trim(filepath) // '/veloc_Pt_bottom.dat')
        open(22,file=trim(filepath) // '/veloc_Ar.dat')
    ! 系のエネルギーデータの出力
        open(30,file=trim(filepath) // '/energy_Pt_top.dat')
        open(31,file=trim(filepath) // '/energy_Pt_bottom.dat')
        open(32,file=trim(filepath) // '/energy_Ar.dat')
        open(35,file=trim(filepath) // '/energy_all.dat')
    ! 系の温度データの出力
        open(40,file=trim(filepath) // '/tempe.dat')
        open(41,file=trim(filepath) // '/tempe_Pt_top_Layer.dat')
        open(42,file=trim(filepath) // '/tempe_Pt_bottom_Layer.dat')
        open(43,file=trim(filepath) // '/tempe_Ar_Layer.dat')
        
        open(45,file=trim(filepath) // '/tempe_Layer.dat')
    ! 系の周期長さの出力
        ! open(50,file=trim(filepath) // '/syuuki.dat')
    ! 熱流束のデータ
        open(60,file=trim(filepath) // '/heatflux.dat')
        open(61,file=trim(filepath) // '/pressure.dat')
        
        open(70,file=trim(filepath) // '/force_phantom.dat')
        open(71,file=trim(filepath) // '/force_interface.dat')

    ! 各分子の最終位置データの出力
        open(80,file=trim(filepath) // '/finpos.dat')
    !　分子の色
        ! open(90,file=trim(filepath) // '/mask.dat')

    ! write(15,'(3I7)') moltype, NUM_MOL_ALL, ndat
    ! do i = 1,ndat
    !     do j = 1, int(NUM_MOL(1)/NUM_Z(1))
    !         write(90,'(I7)') 15      ! 白色
    !     end do
    !     do j = int(NUM_MOL(1)/NUM_Z(1)) + 1, NUM_MOL(1)
    !         write(90,'(I7)') 14      ! 赤色
    !     end do
    !     do j = 1, NUM_MOL(2)
    !         write(90,'(I7)') 7       ! 黄色
    !     end do
    !     do j = 1, int(NUM_MOL(3)/NUM_Z(3))
    !         write(90,'(I7)') 15      ! 白色
    !     end do
    !     do j = int(NUM_MOL(3)/NUM_Z(3)) + 1, NUM_MOL(3)
    !         write(90,'(I7)') 0       ! 青色
    !     end do
    ! end do
    
    ! ターミナルに表示
        write(6,*) ''
        write(6, '(A18, F7.4, A, I8, A)') 'Scaling Time : ', time_scaling, 'ns ', step_scaling, 'step'
        write(6, '(A18, F7.4, A, I8, A)') 'Relaxation Time : ', time_relax, 'ns ', step_relax, 'step'
        write(6, '(A18, F7.4, A, I8, A)') 'Measure Time : ', time_measure, 'ns ', step_measure, 'step'
        write(6,*) ''
        write(6, '(A8, I5)') 'Ar : ', NUM_MOL(Ar) !
        write(6, '(A8, I5, A2)') 'Pt : ', NUM_MOL(TOP), '*2' !
        write(6, '(A8, I5)') 'total : ', NUM_MOL_ALL
        write(6, '(A23, F4.2)') 'Coefficient Strength : ', INTER_STG
        write(6,*) ''
        write(6,*) '----------------------- Scaling Step -----------------------'
        write(6,'(I8, A4)') step_scaling, 'step'
        write(6,*) '' 
    
    step_now = 0

    call seting ! 各分子の初期位置，初期速度などの設定
    call ovito

    do i = 1, step_all
        step_now = i

        ! ターミナルに表示
        if(step_now == step_scaling+1) then
            write(6,*) ''
            write(6,*) '---------------------- Relaxation Step ----------------------'
            write(6,'(I8, A4)') step_relax, 'step'
            write(6,*) ''
        end if
        if(step_now == step_scaling + step_relax+1) then
            write(6,*) ''
            write(6,*) '----------------------- Measure Step ------------------------'
            write(6,'(I8, A4)') step_measure, 'step'
            write(6,*) ''
        end if

        call cpu_time(time)

        ! ステップ数が500の倍数のとき
        if (mod(step_now,500) == 0) then
            write(6,'(3X, I7, 1X, A3, I9, A1, F15.7, A3)') step_now, ' / ', step_all, '(', time, ' s)'
        endif

        ! スケーリング
        if (step_now <= step_scaling .and. mod(step_now,100) == 0) then
            call scaling ! 系内の全分子の温度の補正
        endif

        call calcu ! 各分子に働く力，速度，位置の分子動力学計算
        call bound ! 境界条件の付与

        ! ステップ数が100の倍数+1のとき
        if(mod(step_now, 100) == 1) then
            call record_pos_vel ! 位置，速度を記録
            call record_energy_temp ! エネルギー，温度を記録

            if(step_now <= step_scaling + step_relax) then
                call ovito
            end if

            if(step_now > step_scaling + step_relax) then
                call record_pressure_heatflux ! 熱流束を記録
            end if
        end if
    end do

    call calcu
    call record_energy_temp
    call record_pressure_heatflux
    call record_final ! 最終状態を記録

    contains
    subroutine seting ! 各分子の初期位置，初期速度などの設定
        use enum
        use parameters
        use variable
        use molecules_struct
        implicit none
        integer :: num, i, j, k
        double precision :: coord(XYZ) = 0.0000d0 ! xyz座標 (coordinate)
        double precision :: offset(XYZ)
        double precision :: ran, alpha, beta, cr
        double precision :: v(XYZ)

        do i = 1, PtAr
            coef_force(i) = 24.00d0*EPS(i)/SIG(i)  ! 無次元なことに注意 
        end do
            coef_force(4) = INTER_STG*24.00d0*EPS(4)/SIG(4)
        
        bnd_len(:) = BND_LEN0(:)
        ! write(50,*)bnd_len(1)
        ! write(50,*)bnd_len(2)
        ! write(50,*)bnd_len(3)
        ! write(15,*)BND_LEN0(1), BND_LEN0(2), BND_LEN0(3)
        ! write(15,*)ntime0, ndt

        do i = 1, PtAr
            cutof(:) = bnd_len(:) - CUTOFF*SIG(i)
        end do

        num = 0

        ! 上段のPt配置
        offset(X) = STD_DIST(TOP)*0.25d0
        offset(Y) = STD_DIST(TOP)*0.25d0
        offset(Z) = BND_LEN_Z0 - STD_DIST(TOP)*0.25d0
        do k = 1,NUM_Z(TOP)
            coord(Z) = offset(Z) - dble(k-1)*STD_DIST(TOP)*0.5d0
            do i = 1,NUM_X(TOP)
                coord(X) = offset(X) + dble(i-1)*STD_DIST(TOP)*0.5d0
                do j = 1,NUM_Y(TOP)
                    if(mod(k,2) == 0) then  !z偶数
                        if(mod(i,2) == 0) then
                            coord(Y) = offset(Y) + dble(j-1)*STD_DIST(TOP)   !x偶数
                        else
                            coord(Y) = offset(Y) + dble(j-1)*STD_DIST(TOP) + STD_DIST(TOP)*0.5d0    !x奇数
                        endif
                    else    !z奇数
                        if(mod(i,2) == 0) then
                            coord(Y) = offset(Y) + dble(j-1)*STD_DIST(TOP) + STD_DIST(TOP)*0.5d0    !x偶数
                        else
                            coord(Y) = offset(Y) + dble(j-1)*STD_DIST(TOP)   !x奇数
                        endif
                    endif
                    num = num + 1
                    typ(TOP)%mol(num)%pos(:) = coord(:)
                end do
            end do
        end do

        ! 下段のPt配置
        num = 0
        offset(X) = STD_DIST(BTM)*0.25d0
        offset(Y) = STD_DIST(BTM)*0.25d0
        offset(Z) = STD_DIST(BTM)*0.25d0
        do k = 1,NUM_Z(BTM)
            coord(Z) = offset(Z) + dble(k-1)*STD_DIST(BTM)*0.5d0
            do i = 1,NUM_X(BTM)
                coord(X) = offset(X) + dble(i-1)*STD_DIST(BTM)*0.5d0
                do j = 1,NUM_Y(BTM)
                    if(mod(k,2) == 0) then  !z偶数
                        if(mod(i,2) == 0) then
                            coord(Y) = offset(Y) + dble(j-1)*STD_DIST(BTM)   !x偶数
                        else
                            coord(Y) = offset(Y) + dble(j-1)*STD_DIST(BTM) + STD_DIST(BTM)*0.5d0    !x奇数
                        endif
                    else    !z奇数
                        if(mod(i,2) == 0) then
                            coord(Y) = offset(Y) + dble(j-1)*STD_DIST(BTM) + STD_DIST(BTM)*0.5d0    !x偶数
                        else
                            coord(Y) = offset(Y) + dble(j-1)*STD_DIST(BTM)   !x奇数
                        endif
                    endif
                    num = num + 1
                    typ(BTM)%mol(num)%pos(:) = coord(:)
                end do
            end do
        end do

        ! 中段のAr配置
        num = 0
        offset(X) = BND_LEN0(X)*0.50d0 - STD_DIST(Ar)*(0.25d0*(NUM_X(Ar)  -1))
        offset(Y) = BND_LEN0(Y)*0.50d0 - STD_DIST(Ar)*(0.25d0*(NUM_Y(Ar)*2-1))
        offset(Z) = BND_LEN0(Z)*0.50d0 - STD_DIST(Ar)*(0.25d0*(NUM_Z(Ar)  -1))
        do k = 1,NUM_Z(3)
            coord(Z) = offset(Z) + dble(k-1)*STD_DIST(Ar)*0.5d0
            do i = 1,NUM_X(3)
                coord(X) = offset(X) + dble(i-1)*STD_DIST(Ar)*0.5d0
                do j = 1,NUM_Y(3)
                    if(mod(k,2) == 0) then  !z偶数
                        if(mod(i,2) == 0) then
                            coord(Y) = offset(Y) + dble(j-1)*STD_DIST(Ar)   !x偶数
                        else
                            coord(Y) = offset(Y) + dble(j-1)*STD_DIST(Ar) + STD_DIST(Ar)*0.5d0    !x奇数
                        endif
                    else    !z奇数
                        if(mod(i,2) == 0) then
                            coord(Y) = offset(Y) + dble(j-1)*STD_DIST(Ar) + STD_DIST(Ar)*0.5d0    !x偶数
                        else
                            coord(Y) = offset(Y) + dble(j-1)*STD_DIST(Ar)   !x奇数
                        endif
                    endif
                    num = num + 1
                    typ(Ar)%mol(num)%pos(:) = coord(:)
                end do
            end do
        end do

        cr = 1.00d-6
        do j = 1, PtAr
            do i = 1, NUM_MOL(j)
                call random_number(ran)
                alpha = PI*ran
                call random_number(ran)
                beta = 2.000d0*PI*ran
                v(1) = dsin(alpha)*dcos(beta)*cr
                v(2) = dsin(alpha)*dsin(beta)*cr
                v(3) = dcos(alpha)*cr
                typ(j)%mol(i)%vel(:) = v(:)
            end do
        end do

        do j = 1, Pt
            do i = 1, int(NUM_MOL(j)/NUM_Z(j))        
                typ(j)%mol(i)%vel(:) = 0.000d0
            end do
        end do
    end subroutine seting

    subroutine scaling ! 系内の全分子の温度の補正
        use enum
        use parameters
        use variable
        use molecules_struct
        implicit none
        double precision :: temptp, vel2, aimtem, aimnot, baiss
        integer :: i

        temptp = 0.000d0
        do i = 1, NUM_MOL(Ar)
            vel2 = typ(Ar)%mol(i)%vel(X)**2 + typ(Ar)%mol(i)%vel(Y)**2 + typ(Ar)%mol(i)%vel(Z)**2
            temptp = temptp + vel2
        end do
        temptp = temptp / NUM_MOL(Ar) * 1.000d-16
        aimtem = temp_Ar
        aimnot = 3.000d0 * BOLTZ * aimtem / MASS(Ar)
        baiss = dsqrt(aimnot / temptp)

        ! 速度ベクトルのスケーリング
        ! Arのみ
        do i = 1, NUM_MOL(Ar)
            typ(Ar)%mol(i)%vel(:) = typ(Ar)%mol(i)%vel(:) * baiss
        end do
    end subroutine scaling

    subroutine ovito
        use enum
        use parameters
        use variable
        use molecules_struct
        implicit none
        integer :: i, j
        double precision :: t
        
        t = dble(step_now)*dt*1.0d-6 * 100 ! 100ステップおきに計測　単位は [ns]

        write(16, *) NUM_MOL_ALL
        write(16, '(A8, F4.1, A13, F4.1, A13, F4.1, A39, F6.3)') 'Lattice=', BND_LEN_X0, '0.0 0.0 0.0 ', &
                BND_LEN_Y0, '0.0 0.0 0.0 ', BND_LEN_Z0, '"Properties="species:S1:pos:R:3" Time=', t
        do j = 1, Pt
            do i = 1, NUM_MOL(j)
                write(16, '(A2, 3E15.7)') 'Pt', typ(j)%mol(i)%pos(X), typ(j)%mol(i)%pos(Y), typ(j)%mol(i)%pos(Z)
            end do
        end do
            do i = 1, NUM_MOL(Ar)
                write(16, '(A2, 3E15.7)') 'Ar', typ(Ar)%mol(i)%pos(X), typ(Ar)%mol(i)%pos(Y), typ(Ar)%mol(i)%pos(Z)
            end do
    end subroutine ovito

    subroutine calcu ! 各分子に働く力，速度，位置の分子動力学計算
        use enum
        use parameters
        use variable
        use molecules_struct
        use omp_lib
        implicit none
        integer :: i, j, k, i1, i2
        ! 分子間相互作用
        double precision :: div(3), dist
        double precision :: dit2, dit4, dit6, dit8, dit12, dit14
        double precision :: ppp, force, forVec(3)
        double precision :: vel, sumvene
        ! 粒子登録法
        double precision :: vmax(5000, PtAr), rbook(5000, PtAr)
        logical :: judge_same(5000, 5000, PtAr), judge_diff(NUM_MOL(1), NUM_MOL(2), Pt)
        integer, parameter :: SAFE = 2 ! 安全率
        ! Langevin法
        double precision :: rnd
    
        do j = 1, PtAr
            do i = 1, NUM_MOL(j)
                typ(j)%mol(i)%acc(:) = 0.0000d0
                typ(j)%mol(i)%pot  = 0.0000d0
                typ(j)%mol(i)%kin  = 0.0000d0
            end do
        end do
        force_inter(:,:,:) = 0.000d0
    
        ! 粒子登録法のために出す速度
        do j = 1, PtAr
            do i = 1, NUM_MOL(j)
                vel = dsqrt(typ(j)%mol(i)%vel(X)**2 + typ(j)%mol(i)%vel(Y)**2 + typ(j)%mol(i)%vel(Z)**2)
                if(vmax(i,j) < vel) then
                    vmax(i,j) = vel
                end if
            end do
        end do

        ! 粒子登録法
        ! 最大速度を更新
        if(mod(step_now,step_update) == 1) then
            do j = 1, PtAr
                do i = 1, NUM_MOL(j)
                    rbook(i,j) = CUTOFF*SIG(j) + SAFE*step_update*vmax(i,j)*dt
                    vmax(i,j) = 0.000d0
                end do

                do i1 = 1, NUM_MOL(j)
                    do i2 = i1+1, NUM_MOL(j)
                        div(:) = typ(j)%mol(i1)%pos(:) - typ(j)%mol(i2)%pos(:)
 
                        do k = 1, 2
                            if(div(k) < -bnd_len(k)/2) then
                                div(k) = div(k) + bnd_len(k)
                            else if(div(k) > bnd_len(k)/2) then
                                div(k) = div(k) - bnd_len(k)
                            end if
                        end do

                        dit2 = div(1)**2 + div(2)**2 + div(3)**2
                        dist = dsqrt(dit2)

                        if(rbook(i1,j) > dist .or. rbook(i2,j) > dist) then
                            judge_same(i1,i2,j) = .true.
                        else
                            judge_same(i1,i2,j) = .false.
                        end if
                    end do
                end do
            end do
        end if

        if(mod(step_now,step_update) == 1) then
            do j = 1, Pt
                do i1 = 1, NUM_MOL(j) ! Pt
                    do i2 = 1, NUM_MOL(Ar) ! Ar
                        div(:) = typ(j)%mol(i1)%pos(:) - typ(Ar)%mol(i2)%pos(:)
 
                        do k = 1, 2
                            if(div(k) < -bnd_len(k)/2) then
                                div(k) = div(k) + bnd_len(k)
                            else if(div(k) > bnd_len(k)/2) then
                                div(k) = div(k) - bnd_len(k)
                            end if
                        end do

                        dit2 = div(X)**2 + div(Y)**2 + div(Z)**2
                        dist = dsqrt(dit2)

                        if(rbook(i2,Ar) > dist) then     ! Arの条件のみでよい？
                            judge_diff(i1,i2,j) = .true.
                        else
                            judge_diff(i1,i2,j) = .false.
                        end if
                    end do
                end do
            end do
        end if

        ! 分子間の相互作用力 → ポテンシャルエネルギー
        ! 同じ分子同士の影響
        do j = 1, PtAr
            do i1 = 1, NUM_MOL(j)
                do i2 = i1+1, NUM_MOL(j)
                    if(.not. judge_same(i1,i2,j)) then
                        cycle
                    else
                        div(:) = typ(j)%mol(i1)%pos(:) - typ(j)%mol(i2)%pos(:)
                        ! カットオフ
                        do k = 1, XYZ
                            if (div(k) < -cutof(k)) then
                                div(k) = div(k) + bnd_len(k)
                            else if(div(k) > cutof(k)) then
                                div(k) = div(k) - bnd_len(k)
                            endif
            
                            div(k) = div(k) / SIG(j)
        
                            if (abs(div(k)) > CUTOFF) then
                                cycle
                            endif
                        end do
        
                        dit2 = div(X)**2 + div(Y)**2 + div(Z)**2
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
                        typ(j)%mol(i1)%pot = typ(j)%mol(i1)%pot + ppp*0.500d0
                        typ(j)%mol(i2)%pot = typ(j)%mol(i2)%pot + ppp*0.500d0
        
                        force  = coef_force(j)*(-2.00d0/dit14+1.00d0/dit8)
                        forVec(:) = -force*div(:)
                        typ(j)%mol(i1)%acc(:) = typ(j)%mol(i1)%acc(:) + forVec(:)/MASS(j)
                        typ(j)%mol(i2)%acc(:) = typ(j)%mol(i2)%acc(:) - forVec(:)/MASS(j)
                    end if
                end do
            end do
        end do

        ! 異なる分子同士の影響  ! Ar-Ptの処理    配列は1が上Pt, 2がAr, 3がしたPt
        do j = 1, Pt
            do i1 = 1, NUM_MOL(j)       ! Pt
                do i2 = 1, NUM_MOL(Ar)    ! Ar
                    if(.not. judge_diff(i1,i2,j)) then
                        cycle
                    else
                        div(:) = typ(j)%mol(i1)%pos(:) - typ(Ar)%mol(i2)%pos(:)
            
                        do k = 1, XYZ
                            if (div(k) < -cutof(k)) then
                                div(k) = div(k) + bnd_len(k)
                            else if(div(k) > cutof(k)) then
                                div(k) = div(k) - bnd_len(k)
                            endif
            
                            div(k) = div(k) / SIG(4)
            
                            if (abs(div(k)) > CUTOFF) then
                                cycle
                            endif
                        end do
            
                        dit2 = div(X)**2 + div(Y)**2 + div(Z)**2
                        dist = dsqrt(dit2)
            
                        if(dist > CUTOFF) then
                            cycle
                        endif
            
                        dit4   = dit2*dit2
                        dit6   = dit4*dit2
                        dit8   = dit4*dit4
                        dit12  = dit6*dit6
                        dit14  = dit8*dit6
                        ppp    = INTER_STG*4.00d0*EPS(4)*(1.00d0/dit12-1.00d0/dit6) ! 異分子間ではangCon(接触角)を忘れずに
                        typ(j)%mol(i1)%pot = typ(j)%mol(i1)%pot + ppp*0.500d0
                        typ(Ar)%mol(i2)%pot = typ(Ar)%mol(i2)%pot + ppp*0.500d0
            
                        force  = coef_force(4)*(-2.00d0/dit14+1.00d0/dit8)
                        forVec(:) = -force*div(:)
                        force_inter(i1,j,:) = force_inter(i1,j,:) - forVec(:) ! 無次元なことに注意　符号が逆な気がする
            
                        typ(j)%mol(i1)%acc(:) = typ(j)%mol(i1)%acc(:) + forVec(:)/MASS(j)
                        typ(Ar)%mol(i2)%acc(:) = typ(Ar)%mol(i2)%acc(:) - forVec(:)/MASS(Ar)
                    end if
                end do
            end do
        end do
    
        ! PtのPhantom層はダンパー力とランダム力を付与
        do j = 1, Pt
            do i = int(NUM_MOL(j)/NUM_Z(j)) + 1, 2*int(NUM_MOL(j)/NUM_Z(j)) ! Phantom層のみ
                do k = 1, XYZ
                    rnd = Random() 
                    ! ランダム力
                    force_rand(i,j,k) = rnd * getStddev(temp_Langevin(j)) * 1.000d-9 ! 標準偏差の有次元化
                    ! ダンパー力
                    force_damp(i,j,k) = -DAMP * typ(j)%mol(i)%vel(k) * 1.000d+5 ! 速度の有次元化
                end do
            end do
    
            ! ランダム力とダンパー力を追加
            do i = int(NUM_MOL(j)/NUM_Z(j)) + 1, 2*int(NUM_MOL(j)/NUM_Z(j))         ! 加速度の無次元化 10^-20
                typ(j)%mol(i)%acc(:) = typ(j)%mol(i)%acc(:) + (force_rand(i,j,:)*1.0d+9 + force_damp(i,j,:)*1.0d+9) / MASS(j)*1.000d-3 ! -9+26-20 = -3
            end do
        end do

        ! 運動エネルギー計算
        do j = 1, PtAr
            do i = 1, NUM_MOL(j)
                typ(j)%mol(i)%vmid(:) = typ(j)%mol(i)%vel(:) + typ(j)%mol(i)%acc(:)*0.500d0*dt ! vel(t) = vel(t-dt/2) + acc(t)*dt/2
                sumvene = typ(j)%mol(i)%vmid(X)**2 + typ(j)%mol(i)%vmid(Y)**2 + typ(j)%mol(i)%vmid(Z)**2
                typ(j)%mol(i)%kin = 0.500d0*MASS(j)*sumvene
            end do
        end do
    
        ! パラメータモジュールで初期化するとうまくいかなかったのでここで初期化
        if(step_now == step_scaling + step_relax+1) then 
            heat_phantom(:) = 0.000d0
            heat_interface(:) = 0.000d0
        end if

        ! 熱流束（未完成）を積算
        if(step_now > step_scaling + step_relax) then
            do j = 1, Pt
                do i = int(NUM_MOL(j)/NUM_Z(j)) + 1, 2*int(NUM_MOL(j)/NUM_Z(j)) ! Phantom層  
                    heat_phantom(j) = heat_phantom(j) + &
                        &  (force_rand(i,j,X) + force_damp(i,j,X)) * typ(j)%mol(i)%vmid(X) * 1.000d+5 + &
                        &  (force_rand(i,j,Y) + force_damp(i,j,Y)) * typ(j)%mol(i)%vmid(Y) * 1.000d+5 + &
                        &  (force_rand(i,j,Z) + force_damp(i,j,Z)) * typ(j)%mol(i)%vmid(Z) * 1.000d+5 ! 速さの有次元化 10^5
                end do
                
                do i = 1, NUM_MOL(j) ! Pt分子全体
                    heat_interface(j) = heat_interface(j) + &
                        &  force_inter(i,j,X) * 1.000d-6 * typ(j)%mol(i)%vmid(X) * 1.000d+5 + &
                        &  force_inter(i,j,Y) * 1.000d-6 * typ(j)%mol(i)%vmid(Y) * 1.000d+5 + &
                        &  force_inter(i,j,Z) * 1.000d-6 * typ(j)%mol(i)%vmid(Z) * 1.000d+5
                end do   
            end do
        end if

        ! 数値積分 (蛙跳び法)
        ! Ptの計算
        do j = 1, Pt
            ! 固定層
            do i = 1, int(NUM_MOL(j)/NUM_Z(j))
                typ(j)%mol(i)%vel(:) = 0.0000d0
            end do

            ! その他の層
            do i = int(NUM_MOL(j)/NUM_Z(j)) + 1, int(NUM_MOL(j))
                typ(j)%mol(i)%vel(:) = typ(j)%mol(i)%vel(:) + typ(j)%mol(i)%acc(:) * dt
                typ(j)%mol(i)%pos(:) = typ(j)%mol(i)%pos(:) + typ(j)%mol(i)%vel(:) * dt
            end do
        end do

        ! Arの計算
            do i = 1, NUM_MOL(Ar)
                typ(Ar)%mol(i)%vel(:) = typ(Ar)%mol(i)%vel(:) + typ(Ar)%mol(i)%acc(:) * dt   ! vel(t+dt/2) = vel(t-dt/2) + acc(t)*dt
                typ(Ar)%mol(i)%pos(:) = typ(Ar)%mol(i)%pos(:) + typ(Ar)%mol(i)%vel(:) * dt   ! pos(t+dt)   = pos(t)      + vel(t+dt/2)*dt
            end do
    end subroutine calcu

    subroutine bound ! 境界条件の付与
        use enum
        use parameters
        use variable
        use molecules_struct
        implicit none
        integer :: i, j
        do j = 1, PtAr
            do i = 1, NUM_MOL(j)
                if(typ(j)%mol(i)%pos(X) < 0.00d0) then
                    typ(j)%mol(i)%pos(X) = typ(j)%mol(i)%pos(X) + bnd_len(X)
                else if(typ(j)%mol(i)%pos(X) > bnd_len(X)) then
                    typ(j)%mol(i)%pos(X) = typ(j)%mol(i)%pos(X) - bnd_len(X)
                endif

                if(typ(j)%mol(i)%pos(Y) < 0.00d0) then
                    typ(j)%mol(i)%pos(Y) = typ(j)%mol(i)%pos(Y) + bnd_len(Y)
                else if(typ(j)%mol(i)%pos(Y) > bnd_len(Y)) then
                    typ(j)%mol(i)%pos(Y) = typ(j)%mol(i)%pos(Y) - bnd_len(Y)
                endif
            end do
        end do
    end subroutine bound

    subroutine record_pos_vel ! 位置，速度を記録
        use enum
        use parameters
        use variable
        use molecules_struct
        implicit none
        integer :: i, j
    
        do i = 1, NUM_MOL(TOP)
            ! posit_Pt_top.dat
            write(10, '(I6, 3E15.7)') i, typ(TOP)%mol(i)%pos(X), typ(TOP)%mol(i)%pos(Y), typ(TOP)%mol(i)%pos(Z)
            ! veloc_Pt_top.dat
            write(20, '(I6, 3E15.7)') i, typ(TOP)%mol(i)%vel(X), typ(TOP)%mol(i)%vel(Y), typ(TOP)%mol(i)%vel(Z)
        end do
        do i = 1, NUM_MOL(BTM)
            ! posit_Pt_bottom.dat
            write(11, '(I6, 3E15.7)') i, typ(BTM)%mol(i)%pos(X), typ(BTM)%mol(i)%pos(Y), typ(BTM)%mol(i)%pos(Z)
            ! veloc_Pt_bottom.dat
            write(21, '(I6, 3E15.7)') i, typ(BTM)%mol(i)%vel(X), typ(BTM)%mol(i)%vel(Y), typ(BTM)%mol(i)%vel(Z)
        end do
        do i = 1, NUM_MOL(Ar)
            ! posit_Ar.dat
            write(12, '(I6, 3E15.7)') i, typ(Ar)%mol(i)%pos(X), typ(Ar)%mol(i)%pos(Y), typ(Ar)%mol(i)%pos(Z)
            ! veloc_Ar.dat
            write(22, '(I6, 3E15.7)') i, typ(Ar)%mol(i)%vel(X), typ(Ar)%mol(i)%vel(Y), typ(Ar)%mol(i)%vel(Z)
        end do

        do i = int(NUM_MOL(TOP)/NUM_Z(TOP)) + 1, 2*int(NUM_MOL(TOP)/NUM_Z(TOP)) ! Phantom層
            write(70, '(I6, 6E15.7)') i, force_rand(i,1,X), force_rand(i,1,Y), force_rand(i,1,Z), &
                                         force_damp(i,1,X), force_damp(i,1,Y), force_damp(i,1,Z)
        end do
    
        do i = 1, NUM_MOL(TOP)!(NUM_Z(TOP)-1)*int(NUM_MOL(TOP)/NUM_Z(TOP)) + 1, NUM_MOL(TOP) ! 固液界面層
            write(71, '(I6, 3E15.7)') i, force_inter(i,1,X), force_inter(i,1,Y), force_inter(i,1,Z) ! 上Pt, Ar, 下Pt
        end do
    
        ! ! 可視化用
        ! do j = 1, PtAr
        !     do i = 1, NUM_MOL(j)
        !         ! pos.dat
        !         write(15, '(3E15.7)') typ(j)%mol(i)%pos(X), typ(j)%mol(i)%pos(Y), typ(j)%mol(i)%pos(Z)
        !     end do
        ! end do
    end subroutine record_pos_vel
    
    subroutine record_energy_temp ! エネルギー，温度を記録
        use enum
        use parameters
        use variable
        use molecules_struct
        implicit none
        double precision, dimension(PtAr) :: totEne, totPot, totKin, temp
        double precision :: allEne, allPot, allKin
        double precision :: kinPtTmp, kinArTmp(NUM_DIV_Ar)
        integer :: cnt(NUM_DIV_Ar)
        double precision :: temp_layer_Pt_(NUM_Z(TOP),PtAr)
        double precision :: temp_layer_Ar_(NUM_DIV_Ar)
        integer :: i, j, k
        integer :: step
        step = step_now-step_relax-step_scaling
    
        allEne = 0.000d0
        allPot = 0.000d0
        allKin = 0.000d0
        totEne(:) = 0.000d0
        totPot(:) = 0.000d0
        totKin(:) = 0.000d0
        temp(:) = 0.000d0
        temp_layer_Pt_(:,:) = 0.000d0
        temp_layer_Ar_(:) = 0.000d0
    
        do j = 1, PtAr
            ! ポテンシャル
            do i = 1, NUM_MOL(j)
                totPot(j) = totPot(j) + typ(j)%mol(i)%pot
            end do
            totPot(j) = totPot(j) * 1.000d-16
        end do
    
        ! 運動エネルギー
        ! Pt
        do j = 1, Pt
            do k = 2, NUM_Z(j) ! Ptの層の数
                kinPtTmp = 0.000d0
                do i = (k-1)*int(NUM_MOL(j)/NUM_Z(j)) + 1, k*int(NUM_MOL(j)/NUM_Z(j))
                    kinPtTmp = kinPtTmp + typ(j)%mol(i)%kin
                end do

                kinPtTmp = kinPtTmp * 1.000d-16
                totKin(j) = totKin(j) + kinPtTmp
                temp_layer_Pt_(k, j) = 2.0d0 * kinPtTmp / (3.0d0 * dble(NUM_MOL(j)/NUM_Z(j)) * BOLTZ)
                temp(j) = temp(j) + temp_layer_Pt_(k, j)
            end do

            temp(j) = temp(j) / dble(NUM_Z(j)-1)
        end do

        ! Ar
        kinArTmp(:) = 0.000d0
        cnt(:) = 0
        do i = 1, NUM_MOL(Ar)
            do k = 1, NUM_DIV_Ar
                if ( (k-1)*zdiv <= (typ(Ar)%mol(i)%pos(Z) - THICK(BTM)) .and. (typ(Ar)%mol(i)%pos(Z) - THICK(BTM)) < k*zdiv) then
                    kinArTmp(k) = kinArTmp(k) + typ(Ar)%mol(i)%kin
                    cnt(k) = cnt(k) + 1
                    cycle
                end if
            end do
        end do

        do k = 1, NUM_DIV_Ar
            kinArTmp(k) = kinArTmp(k) * 1.000d-16
            totKin(Ar) = totKin(Ar) + kinArtmp(k)
            temp_layer_Ar_(k) = 2.0d0 * kinArTmp(k) / (3.0d0 * dble(cnt(k)) * BOLTZ)
            temp(Ar) = temp(Ar) + temp_layer_Ar_(k)
        end do

        temp(Ar) = temp(Ar) / dble(NUM_DIV_Ar)

        do i = 1, PtAr
            totEne(j) = totPot(j) + totKin(j)
            allEne = allEne + totEne(j)
            allPot = allPot + totPot(j)
            allKin = allKin + totKin(j)
        end do
    
        ! Pt
        do j = 1, Pt
            do i = 2, NUM_Z(j)
                temp_layer_Pt(i,j) = temp_layer_Pt(i,j) + temp_layer_Pt_(i,j)
            end do
        end do
    
        ! Ar
            do i = 1, NUM_DIV_Ar
                temp_layer_Ar(i) = temp_layer_Ar(i) + temp_layer_Ar_(i)
            end do

        write(30, '(5E15.7)') step_now*int(dt)*1.0d-6, totEne(TOP), totPot(TOP), totKin(TOP) ! energy_Pt_top.dat
        write(31, '(5E15.7)') step_now*int(dt)*1.0d-6, totEne(BTM), totPot(BTM), totKin(BTM) ! energy_Pt_bottom.dat
        write(32, '(5E15.7)') step_now*int(dt)*1.0d-6, totEne(Ar),  totPot(Ar),  totKin(Ar)    ! energy_Ar.dat
        write(35, '(5E15.7)') step_now*int(dt)*1.0d-6, allEne, allPot, allKin           ! energy_all.dat
        write(40, '(5E15.7)') step_now*int(dt)*1.0d-6, temp(TOP), temp(Ar), temp(BTM)        ! tempe.dat
    
        ! ! !!!!!!!!!Pt層を増やすとき必ず変更すること!!!!!!!!!
        ! write(41, '(5E15.7)') step*int(dt)*1.0d-6, temp_layer_Pt_(1,1), temp_layer_Pt_(2,1), temp_layer_Pt_(3,1), temp_layer_Pt_(4,1)
        ! write(42, '(5E15.7)') step*int(dt)*1.0d-6, temp_layer_Pt_(1,3), temp_layer_Pt_(2,3), temp_layer_Pt_(3,3), temp_layer_Pt_(4,3)
        ! write(43, '(16E15.7)') step*int(dt)*1.0d-6, temp_layer_Ar_(1),  temp_layer_Ar_(2),  temp_layer_Ar_(3),  temp_layer_Ar_(4),  &
        !                                           & temp_layer_Ar_(5),  temp_layer_Ar_(6),  temp_layer_Ar_(7),  temp_layer_Ar_(8),  &
        !                                           & temp_layer_Ar_(9),  temp_layer_Ar_(10), temp_layer_Ar_(11), temp_layer_Ar_(12), &
        !                                           & temp_layer_Ar_(13), temp_layer_Ar_(14), temp_layer_Ar_(15)
    end subroutine record_energy_temp
    
    subroutine record_pressure_heatflux ! 熱流束を記録
        use enum
        use parameters
        use variable
        use molecules_struct
        implicit none
        integer :: i, j
        double precision :: flux_phantom(Pt) ! Phantom層からの熱輸送量
        double precision :: flux_phantom(Pt) ! 固液界面での熱輸送量
        integer :: step
        step = step_now-step_relax-step_scaling
    
        flux_phantom(:) = 0.000d0
        flux_interface(:) = 0.000d0
        pressure(:) = 0.000d0
    
        do j = 1, Pt
            flux_phantom(j) = heat_phantom(j) / (AREA_Pt * 1.000d-20) * 1.000d-15 ! 速さの有次元化 10^5
            flux_interface(j) = heat_interface(j) / (AREA_Pt * 1.000d-20) * 1.000d-15 ! 面積の有次元化 10^-20
            
            do i = 1, NUM_MOL(j)
                if(j == 1) then
                    pressure(TOP) = pressure(TOP) + force_inter(i,TOP,Z)*1.000d-6 / (AREA_Pt * 1.000d-20) *1.000d-6 ! [MPa]　圧力のオーダーはあってそう
                else
                    pressure(BTM) = pressure(BTM) - force_inter(i,BTM,Z)*1.000d-6 / (AREA_Pt * 1.000d-20) *1.000d-6
                end if
            end do
        end do
    
        write(60,'(5E15.7)') step*int(dt)*1.0d-6, flux_phantom(TOP), flux_phantom(BTM), flux_interface(TOP), flux_interface(BTM)
        write(61,'(3E15.7)') step*int(dt)*1.0d-6, pressure(TOP), pressure(BTM)
    end subroutine record_pressure_heatflux
    
    subroutine record_final ! 最終状態を記録
        use enum
        use parameters
        use variable
        use molecules_struct
        implicit none
        integer :: i, j
        double precision :: disz
    
        ! Pt
        do j = 1, Pt
            do i = 2, NUM_Z(j)
                temp_layer_Pt(i,j) = temp_layer_Pt(i,j) / dble(step_measure/(tau/dt))
            end do
        end do
        ! Ar
        do i = 1, NUM_DIV_Ar
            temp_layer_Ar(i) = temp_layer_Ar(i) / dble(step_measure/(tau/dt))
        end do

        ! do j = 1, PtAr
        !     do i = 1, NUM_MOL(j)
        !         ! syuuki.dat
        !         write(50, '(I6, 6E15.7)') & 
        !         i, typ(j)%mol(i)%pos(X), typ(j)%mol(i)%pos(Y), typ(j)%mol(i)%pos(Z), &
        !            typ(j)%mol(i)%vel(X), typ(j)%mol(i)%vel(Y), typ(j)%mol(i)%vel(Z)
        !     end do
        ! end do
    
        disz = STD_DIST(BTM)*0.25d0
        do i = 2, NUM_Z(BTM)   ! 固定層は除外
            disz = disz + STD_DIST(BTM)*0.5d0
            write(45, '(2E15.7)') disz, temp_layer_Pt(i,BTM)
        end do
    
        disz = disz - zdiv * 0.5d0
        do i = 1, NUM_DIV_Ar
            disz = disz + zdiv
            write(45, '(2E15.7)') disz, temp_layer_Ar(i)
        end do
    
        disz = BND_LEN_Z0 - STD_DIST(TOP)*0.25d0
        do i = 2, NUM_Z(TOP)
            disz = disz - STD_DIST(TOP)*0.5d0
            write(45, '(2E15.7)') disz, temp_layer_Pt(i,TOP)
        end do
    
    end subroutine record_final
end program main