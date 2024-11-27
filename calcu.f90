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
    interForce(:,:,:) = 0.000d0

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
    
                    div(k) = div(k) / SIG(4)

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
                forVec(:) = -force*div(:)
                interForce(i1,:,j) = interForce(i1,:,j) - forVec(:) ! 無次元なことに注意　符号が逆な気がする

                typ(j)%mol(i1)%acc(:) = typ(j)%mol(i1)%acc(:) + forVec(:)/MASS(j)
                typ(2)%mol(i2)%acc(:) = typ(2)%mol(i2)%acc(:) - forVec(:)/MASS(2)
            end do
        end do
    end do

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
                dmpForce(i,k,j) = -DAMP * typ(j)%mol(i)%vtmp(k) * 1.000d+5 ! 速度の有次元化
            end do
        end do

        ! ランダム力とダンパー力を追加
        do i = int(nummol(j)/numz(j)) + 1, 2*int(nummol(j)/numz(j))         ! 加速度の無次元化 10^-20
            typ(j)%mol(i)%acc(:) = typ(j)%mol(i)%acc(:) + (rndForce(i,:,j)*1.0d+9 + dmpForce(i,:,j)*1.0d+9) / MASS(j)*1.000d-3 ! -9+26-20 = -3
        end do
    end do

    ! 運動エネルギー計算
    do j = 1, TYPMOL
        do i = 1, nummol(j)
            typ(j)%mol(i)%vtmp(:) = typ(j)%mol(i)%vel(:) + typ(j)%mol(i)%acc(:)*0.500d0*dt ! vel(t) = vel(t-dt/2) + acc(t)*dt/2
            sumvene = typ(j)%mol(i)%vtmp(1)**2 + typ(j)%mol(i)%vtmp(2)**2 + typ(j)%mol(i)%vtmp(3)**2
            typ(j)%mol(i)%kinet = 0.500d0*MASS(j)*sumvene
        end do
    end do

    ! 熱流束（未完成）を積算
    if(stpNow > stpRelax) then
        do j = 1, TYPMOL
            if(stpNow == stpRelax+1) then   ! パラメータモジュールで初期化するとうまくいかなかったのでここで初期化
                heatPhantom(:) = 0.000d0
                heatInterface(:) = 0.000d0
            end if

            if (j == 2) then
                cycle
            end if
    
            do i = int(nummol(j)/numz(j)) + 1, 2*int(nummol(j)/numz(j)) ! Phantom層  
                heatPhantom(j) = heatPhantom(j) + ((rndForce(i,1,j) + dmpForce(i,1,j))*typ(j)%mol(i)%vtmp(1) +  (rndForce(i,2,j) + dmpForce(i,2,j))*typ(j)%mol(i)%vtmp(2) + (rndForce(i,3,j) + dmpForce(i,3,j))*typ(j)%mol(i)%vtmp(3)) * 1.000d+5 ! 速さの有次元化 10^5
            end do
    
            do i = 1, nummol(j) ! Pt分子全体
                heatInterface(j) = heatInterface(j) + (interForce(i,1,j)*1.000d-6 * typ(j)%mol(i)%vtmp(1)*1.000d+5 + interForce(i,2,j)*1.000d-6 * typ(j)%mol(i)%vtmp(2)*1.000d+5 + interForce(i,3,j)*1.000d-6 * typ(j)%mol(i)%vtmp(3)*1.000d+5)
            end do
        end do
    end if

    ! 数値積分 (蛙跳び法)
    do j = 1, TYPMOL
        ! Arの計算
        if(j == 2) then
            do i = 1, nummol(2)
                typ(2)%mol(i)%vel(:) = typ(2)%mol(i)%vel(:) + typ(2)%mol(i)%acc(:) * dt   ! vel(t+dt/2) = vel(t-dt/2) + acc(t)*dt
                typ(2)%mol(i)%pos(:) = typ(2)%mol(i)%pos(:) + typ(2)%mol(i)%vel(:) * dt   ! pos(t+dt)   = pos(t)      + vel(t+dt/2)*dt
            end do
        else
        ! Ptの計算
            ! 固定層
            do i = 1, int(nummol(j)/numz(j))
                typ(j)%mol(i)%vel(:) = 0.0000d0
            end do

            ! その他の層
            do i = int(nummol(j)/numz(j)) + 1, int(nummol(j))
                typ(j)%mol(i)%vel(:) = typ(j)%mol(i)%vel(:) + typ(j)%mol(i)%acc(:) * dt
                typ(j)%mol(i)%pos(:) = typ(j)%mol(i)%pos(:) + typ(j)%mol(i)%vel(:) * dt
            end do
        end if
    end do
end subroutine calcu