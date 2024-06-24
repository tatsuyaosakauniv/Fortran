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