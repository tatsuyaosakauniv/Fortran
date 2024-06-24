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