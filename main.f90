program main
    use parameters
    use variable, !only: stpNow
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
    ! 各分子の位置データの出力
        open(10,file='/Users/tatsuya/fortran/output/posit_PtUp.dat')
        open(11,file='/Users/tatsuya/fortran/output/posit_Ar.dat')
        open(12,file='/Users/tatsuya/fortran/output/posit_PtLw.dat')
    ! 可視化用のpvch.fを移植 
        open(15,file='/Users/tatsuya/fortran/exe/pos.dat')
    ! 各分子の速度データの出力
        open(20,file='/Users/tatsuya/fortran/output/veloc_PtUp.dat')
        open(21,file='/Users/tatsuya/fortran/output/veloc_Ar.dat')
        open(22,file='/Users/tatsuya/fortran/output/veloc_PtLw.dat')
    ! 系のエネルギーデータの出力
        open(30,file='/Users/tatsuya/fortran/output/energy_PtUp.dat')
        open(31,file='/Users/tatsuya/fortran/output/energy_Ar.dat')
        open(32,file='/Users/tatsuya/fortran/output/energy_PtLw.dat')
        open(35,file='/Users/tatsuya/fortran/output/energy_all.dat')
    ! 系の温度データの出力
        open(40,file='/Users/tatsuya/fortran/output/tempe.dat')
        open(41,file='/Users/tatsuya/fortran/output/tempe_PtUp_Layer.dat')
        open(42,file='/Users/tatsuya/fortran/output/tempe_Ar_Layer.dat')
        open(43,file='/Users/tatsuya/fortran/output/tempe_PtLw_Layer.dat')
        
        open(45,file='/Users/tatsuya/fortran/output/tempe_Layer.dat')
    ! 系の周期長さの出力
        open(50,file='/Users/tatsuya/fortran/output/syuuki.dat')
    ! 熱流束のデータ
        open(60,file='/Users/tatsuya/fortran/output/heatflux.dat')
        open(61,file='/Users/tatsuya/fortran/output/pressure.dat')
        
        open(70,file='/Users/tatsuya/fortran/output/force_phantom.dat')
        open(71,file='/Users/tatsuya/fortran/output/force_interface.dat')

    ! 各分子の最終位置データの出力
        open(80,file='/Users/tatsuya/fortran/output/finpos.dat')
    !　分子の色
        open(90,file='/Users/tatsuya/fortran/exe/mask.dat')

    write(15,'(3I7)') moltype, tlnkoss, ndat
    do i = 1,ndat
        do j = 1, int(nummol(1)/numz(1))
            write(90,'(I7)') 15      ! 白色
        end do
        do j = int(nummol(1)/numz(1)) + 1, nummol(1)
            write(90,'(I7)') 14      ! 赤色
        end do
        do j = 1, nummol(2)
            write(90,'(I7)') 1       ! 黄色
        end do
        do j = 1, int(nummol(3)/numz(3))
            write(90,'(I7)') 15      ! 白色
        end do
        do j = int(nummol(3)/numz(3)) + 1, nummol(3)
            write(90,'(I7)') 0       ! 青色
        end do
    end do
    
    ! ターミナルに表示
    write(6,*) ''
    write(6, '(A17, F7.4, A2)') 'Scaling Time: ', timeScaling, 'ns'
    write(6, '(A17, F7.4, A2)') 'Relaxation Time: ', timeRelax, 'ns'
    write(6, '(A17, F7.4, A2)') 'Measure Time: ', timeMeasure, 'ns'
    write(6,*) ''
    write(6,'(A16)') '- Scaling Step -'
    write(6,*) '' 
    
    write(6,*) stpScaling, stpRelax, stpMax
    write(6,*) xsyul0, ysyul0, zsyul0
    write(6,*) zdiv

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
        if(stpNow == stpRelax) then
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

        if(stpNow >= stpRelax .and. mod(stpNow, int(tau/dt)) == 0) then
            call record_energy_temp ! エネルギー，温度を記録
            call record_pressure_heatflux ! 熱流束を記録
        end if
        
        ! ステップ数が100の倍数+1のとき
        if(mod(stpNow, 100) == 1) then
            call record_pos_vel ! 位置，速度を記録
        end if
    end do

    call record_final ! 最終状態を記録

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
            forCoef(i) = 24.00d0*EPS(i)/SIG(i)  ! 無次元なことに注意 *1.0d-6
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

        do j = 1, TYPMOL
            if(j == 2) then
                cycle
            end if

            do i = 1, int(nummol(j)/numz(j))        
                typ(j)%mol(i)%vel(:) = 0.000d0
            end do
        end do
    end subroutine seting

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
end program main