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