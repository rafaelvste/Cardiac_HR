    ! Functions and subroutines:
    module funcs_subrs
        implicit none
        
    contains
    
        !HR model
        subroutine hr_model(b, X_f, Y_f, Z_f, D_X, D_Y, D_Z)
            implicit none
            
            ! Declaration of variables:
            real(kind=8), intent(in) :: b, X_f, Y_f, Z_f
            real(kind=8), intent(out) :: D_X, D_Y, D_Z
            real(kind=8) :: a, c, d, r, s, Xr, Ie
            
            ! HR model parameters:
            a = 1.0d0
            c = 1.0d0
            d = 5.0d0
            r = 0.01d0
            s = 4.0d0
            Xr = -1.6d0
            Ie = 2.5d0
            
            ! HR equations:
            D_X = Y_f - a*(X_f**3.0d0) + b*(X_f**2.0d0) - Z_f + Ie
            D_Y = c - d*(X_f**2.0d0) - Y_f
            D_Z = r*(s*(X_f - Xr) - Z_f)
            
        end subroutine hr_model
        
        !4th order Runge-Kutta:
        subroutine rk4(h, b, X_o, Xf)
            implicit none
            
            !Declaration of variables:
            real(kind=8), intent(in) :: h, b
            real(kind=8), dimension(3), intent(in) :: X_o
            real(kind=8), dimension(3), intent(inout) :: Xf
            integer :: i, j
            real(kind=8), dimension(3) :: DX
            real(kind=8), dimension(4,3) :: k

            ! Iterate RK slopes:
            do i = 1, 4
                ! Calculate the HR equations (x = X(1), y = X(2), z = X(3)):
                call hr_model(b, Xf(1), Xf(2), Xf(3), DX(1), DX(2), DX(3))
                ! Calculate the RK4 equations for x, y and z:
                k(i,:) = DX         
                if ((i == 1) .or. (i == 2)) then
                    Xf = X_o + (h/2.0d0)*k(i,:)
                else if (i == 3) then
                    Xf = X_o + h*k(i,:)
                else if (i == 4) then
                    Xf = X_o + (h/6.0d0)*(k(1,:) + 2.0d0*(k(2,:) + k(3,:)) + k(4,:))
                end if
            end do

        end subroutine rk4
        
        ! Iterate time and measure the APDs and ISIs for a value of b
        subroutine bifurc_APD_ISI(b)
            implicit none
            
            ! Declaration of variables:
            real(kind=8), intent(in) :: b
            integer :: t1, t1_max
            real(kind=8) :: rk_step, t_max, t_trans, t, t_apd, t_isi, APD, ISI
            real(kind=8), dimension(3) :: Xo, X
            
            ! Discretize time:
            t_max = 2000.0d0 !Max time
            t_trans = 1000.0d0 !Transient
            rk_step = 0.01d0 !RK time step in the model time units
            t1_max = nint(t_max/rk_step) !Max number of time steps
            t = 0.0d0
            
            t_apd = 0.0d0 !Store the start of the APD interval
            t_isi = 0.0d0 !Store the start of the ISI interval
            
            ! Iterate time:
            do t1 = 1, t1_max
                t = t + rk_step
                ! Initial conditions (x = X(1), y = X(2), z = X(3)):
                if (t1 == 1) then 
                    Xo(1) = 0.0d0
                    Xo(2) = 0.0d0
                    Xo(3) = 0.0d0
                    X = 0.0d0
                ! Update variables:
                else if (t1 > 1) then 
                    Xo = X
                end if
                ! 4th order Runge-Kutta with HR model:
                call rk4(rk_step, b, Xo, X)
                !After transient:
                if (t >= t_trans) then
                    !Calculate the APD:
                    if ((X(1)*Xo(1) < 0.0d0) .and. (X(1) > Xo(1))) then
                        t_apd = t
                    else if ((X(1)*Xo(1) < 0.0d0) .and. (X(1) < Xo(1))) then
                        if (t_apd > 0.0d0) then
                            APD = t - t_apd
                            t_apd = 0.0d0
                            !Write to file:
                            write(11,*) b, APD
                        end if
                    end if
                    !Calculate the ISI:
                    if ((X(1)*Xo(1) < 0.0d0) .and. (X(1) < Xo(1))) then
                        if (t_isi == 0.0d0) then
                            t_isi = t
                        else if (t_isi > 0.0d0) then
                            ISI = t - t_isi
                            t_isi = t
                            !Write to file:
                            write(12,*) b, ISI
                        end if
                    end if
                end if
            end do
            
        end subroutine bifurc_APD_ISI
        
    end module funcs_subrs
    
    !Main program:
    program HR_cell_bifurc_APD_ISI
        use funcs_subrs
        implicit none

        !Calculates the action potential duration (APD) and interspike interval (ISI)
        !bifurcation diagrams of the HR model. The HR model is solved with a 4th order
        !Runge-Kutta method with a fixed time step. Outputs the bifurcation diagrams
        !by varying parameter b.
        
        ! Declaration of variables:
        integer :: n_b, n_b_max
        real(kind=8) :: b_0, b_f, db, b_p
        
        ! Start and end of the parameter b range:
        b_0 = 0.0d0 
        b_f = 4.0d0

        ! Interval between points:
        db = 0.001d0

        ! Number of points in the b axis of the diagram:
        n_b_max = nint((b_f-b_0)/db)
        
        ! Open output files:
        open(unit=11, file="bifurc_apd.dat", status="replace")
        open(unit=12, file="bifurc_isi.dat", status="replace")
        
        ! Iterate paramater b:
        b_p = b_0 - db
        do n_b = 1, n_b_max
            b_p = b_p + db
            ! Determine maximum and minimum ISI in the time series:
            call bifurc_APD_ISI(b_p)
        end do
        
        ! CLose output files:
        close(unit=11)
        close(unit=12)
        
        ! Gnuplot script
        ! Plot the APD bifucation diagram against parameter b:
        open(unit=13, file="bifurc_apd.plt", status="replace")
            write(13, '(a)') "set terminal wxt"
            write(13, '(a)') "set grid"
            write(13, '(a)') "unset key"           
            write(13, '(a)') "set tics font ', 15'"
            write(13, '(a)') "set pointsize 0.4"   
            write(13, '(a)') "set xlabel 'b' font ', 20'"
            write(13, '(a)') "set ylabel 'APD' font ', 20'"
            write(13, '(a)') "plot 'bifurc_apd.dat' pointtype 21 lt 8"
        close(unit=13)
        call system("gnuplot -p bifurc_apd.plt")
        
        ! Plot the ISI bifucation diagram against parameter b:
        open(unit=14, file="bifurc_isi.plt", status="replace")
            write(14, '(a)') "set terminal wxt"
            write(14, '(a)') "set grid"
            write(14, '(a)') "unset key"           
            write(14, '(a)') "set tics font ', 15'"
            write(14, '(a)') "set pointsize 0.4"       
            write(14, '(a)') "set xlabel 'b' font ', 20'"
            write(14, '(a)') "set ylabel 'ISI' font ', 20'"
            write(14, '(a)') "plot 'bifurc_isi.dat' pointtype 21 lt 8"
        close(unit=14)      
        call system("gnuplot -p bifurc_isi.plt")
        
    end program HR_cell_bifurc_APD_ISI