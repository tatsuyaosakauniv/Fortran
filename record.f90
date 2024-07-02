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

    do i = int(nummol(1)/numz(1)) + 1, 2*int(nummol(1)/numz(1)) ! Phantom層
        write(70, '(I6, 6E15.7)') i, rndForce(i,1,1), rndForce(i,2,1), rndForce(i,3,1), dmpForce(i,1,1), dmpForce(i,2,1), dmpForce(i,3,1)
    end do

    do i = 1, nummol(1)!(numz(1)-1)*int(nummol(1)/numz(1)) + 1, nummol(1) ! 固液界面層
        write(71, '(I6, 3E15.7)') i, interForce(i,1), interForce(i,2), interForce(i,3) ! 上Pt, Ar, 下Pt
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
    double precision :: allEne, allPot, allKin
    double precision :: kinPtTmp, kinArTmp(numDivAr)
    integer :: cnt(numDivAr)
    double precision :: tempLayerPt_(numz(1),TYPMOL)
    double precision :: tempLayerAr_(numDivAr)
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
            ! do i = 1, nummol(j)
            !     totKin(j) = totKin(j) + typ(j)%mol(i)%kinet
            ! end do
            ! totKin(j) = totKin(j) * 1.000d-16
            ! temp(j) = 2.0d0 * totKin(j) / (3.0d0 * dble(nummol(j)) * BOLTZ)

            kinArTmp(:) = 0.000d0
            cnt(:) = 0
            do i = 1, nummol(j)
                do k = 1, numDivAr
                    if ( (k-1)*zdiv <= (typ(j)%mol(i)%pos(3) - thick(3)) .and. (typ(j)%mol(i)%pos(3) - thick(3)) < k*zdiv) then
                        kinArTmp(k) = kinArTmp(k) + typ(j)%mol(i)%kinet
                        cnt(k) = cnt(k) + 1
                        cycle
                    end if
                end do
            end do

            do k = 1, numDivAr
                kinArTmp(k) = kinArTmp(k) * 1.000d-16
                totKin(j) = totKin(j) + kinArtmp(k)
                tempLayerAr_(k) = 2.0d0 * kinArTmp(k) / (3.0d0 * dble(cnt(k)) * BOLTZ)
                temp(j) = temp(j) + tempLayerAr_(k)
            end do

            temp(j) = temp(j) / dble(numDivAr)
        else
            ! Pt
            do k = 2, numz(j) ! Ptの層の数
                kinPtTmp = 0.000d0
                do i = (k-1)*int(nummol(j)/numz(j)) + 1, k*int(nummol(j)/numz(j))
                    kinPtTmp = kinPtTmp + typ(j)%mol(i)%kinet
                end do

                kinPtTmp = kinPtTmp * 1.000d-16
                totKin(j) = totKin(j) + kinPtTmp
                tempLayerPt_(k, j) = 2.0d0 * kinPtTmp / (3.0d0 * dble(nummol(j)/numz(j)) * BOLTZ)
                temp(j) = temp(j) + tempLayerPt_(k, j)
            end do

            temp(j) = temp(j) / dble(numz(j)-1)
        end if

        totEne(j) = totPot(j) + totKin(j)
        allEne = allEne + totEne(j)
        allPot = allPot + totPot(j)
        allKin = allKin + totKin(j)
    end do

    do j = 1, TYPMOL
        if(j == 2) then
            do i = 1, numDivAr
                tempLayerAr(i) = tempLayerAr(i) + tempLayerAr_(i)
            end do
        else
            do i = 2, numz(j)
                tempLayerPt(i,j) = tempLayerPt(i,j) + tempLayerPt_(i,j)
            end do
        end if
    end do

    write(30, '(I6, 4E15.7)') (stpNow+99)*int(dt), totEne(1), totPot(1), totKin(1)  ! energy_PtUp.dat
    write(31, '(I6, 4E15.7)') (stpNow+99)*int(dt), totEne(2), totPot(2), totKin(2)  ! energy_Ar.dat
    write(32, '(I6, 4E15.7)') (stpNow+99)*int(dt), totEne(3), totPot(3), totKin(3)  ! energy_PtLw.dat
    write(35, '(I6, 4E15.7)') (stpNow+99)*int(dt), allEne, allPot, allKin           ! energy_all.dat
    write(40, '(I6, 4E15.7)') (stpNow+99)*int(dt), temp(1), temp(2), temp(3)        ! tempe.dat

    !!!!!!!!!Pt層を増やすとき必ず変更すること!!!!!!!!!
    write(41, '(I6, 4E15.7)') (stpNow+99)*int(dt), tempLayerPt_(1,1), tempLayerPt_(2,1), tempLayerPt_(3,1), tempLayerPt_(4,1)
    write(42, '(I6, 15E15.7)') (stpNow+99)*int(dt), tempLayerAr_(1), tempLayerAr_(2), tempLayerAr_(3), tempLayerAr_(4), tempLayerAr_(5), tempLayerAr_(6), tempLayerAr_(7), tempLayerAr_(8), tempLayerAr_(9), tempLayerAr_(10), tempLayerAr_(11), tempLayerAr_(12), tempLayerAr_(13), tempLayerAr_(14), tempLayerAr_(15)
    write(43, '(I6, 4E15.7)') (stpNow+99)*int(dt), tempLayerPt_(1,3), tempLayerPt_(2,3), tempLayerPt_(3,3), tempLayerPt_(4,3)
end subroutine record_energy_temp

subroutine record_pressure_heatflux ! 熱流束を記録
    use parameters
    use variable
    use molecules_struct
    implicit none
    integer :: i, j
    integer :: stp
    stp = stpNow-stpRelax

    pressure(:) = 0.000d0

    do j = 1, TYPMOL
        if (j == 2) then
            cycle
        end if

        do i = int(nummol(j)/numz(j)) + 1, 2*int(nummol(j)/numz(j)) ! Phantom層
            
            heatPhantom(j) = heatPhantom(j) + (rndForce(i,3,j) + dmpForce(i,3,j))*typ(j)%mol(i)%vtmp(3) * 1.000d+5 / (areaPt * 1.000d-20) * tau * 1.000d-15 ! 速さの有次元化 10^5
        end do

        ! do i = (numz(j)-1)*int(nummol(j)/numz(j)) + 1, nummol(j) ! 固液界面層
        !     heatInterface(j) = heatInterface(j) + interForce(i,j)*1.000d-6 * typ(j)%mol(i)%vtmp(3) * 1.000d+5 / (areaPt * 1.000d-20) * tau * 1.000d-15 ! 面積の有次元化 10^-20
        ! end do

        do i = 1, nummol(j)
            heatInterface(j) = heatInterface(j) + interForce(i,j)*1.000d-6 * typ(j)%mol(i)%vtmp(3) * 1.000d+5 / (areaPt * 1.000d-20) * tau * 1.000d-15 ! 面積の有次元化 10^-20
        end do
        
        do i = 1, nummol(j)
            if(j == 1) then
                pressure(1) = pressure(1) + interForce(i,1)*1.000d-6 / (areaPt * 1.000d-20) *1.000d-6 ! [MPa]　圧力のオーダーはあってそう
            else
                pressure(3) = pressure(3) - interForce(i,3)*1.000d-6 / (areaPt * 1.000d-20) *1.000d-6
            end if
        end do
    end do

    write(60,'(I6, 6E15.7)') (stp)*int(dt), heatPhantom(1), heatPhantom(3), heatInterface(1), heatInterface(3)
    write(61,'(I6, 3E15.7)') (stp)*int(dt), pressure(1), pressure(2), pressure(3)
end subroutine record_pressure_heatflux

subroutine record_final ! 最終状態を記録
    use parameters
    use variable
    use molecules_struct
    implicit none
    integer :: i, j
    double precision :: disz

    do j = 1, TYPMOL
        if(j == 2) then
            do i = 1, numDivAr
                tempLayerAr(i) = tempLayerAr(i) / dble((stpMax)/100)
            end do
        else
            do i = 2, numz(j)
                tempLayerPt(i,j) = tempLayerPt(i,j) / dble((stpMax)/100)
            end do
        end if
    end do
    
    do j = 1, TYPMOL
        do i = 1, nummol(j)
            ! syuuki.dat
            write(50, '(I6, 6E15.7)') & 
            i, typ(j)%mol(i)%pos(1), typ(j)%mol(i)%pos(2), typ(j)%mol(i)%pos(3), &
               typ(j)%mol(i)%vel(1), typ(j)%mol(i)%vel(2), typ(j)%mol(i)%vel(3)
        end do
    end do

    disz = STDIST(3)*0.25d0
    do i = 1, numz(3)
        disz = disz + STDIST(3)*0.5d0
        write(45, '(2E15.7)') disz, tempLayerPt(i,3)
    end do

    disz = disz - zdiv * 0.5d0
    do i = 1, numDivAr
        disz = disz + zdiv
        write(45, '(2E15.7)') disz, tempLayerAr(i)
    end do

    disz = zsyul0 - STDIST(1)*0.25d0
    do i = 1, numz(1)
        disz = disz - STDIST(1)*0.5d0
        write(45, '(2E15.7)') disz, tempLayerPt(i,1)
    end do


end subroutine record_final