subroutine cortra ! 系内の全分子の並進速度の補正
    use variable, !only: velx, vely, velz
    use parameters
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
    use variable
    use parameters
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