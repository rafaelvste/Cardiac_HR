    ! Functions and subroutines:
    module funcs_subrs
        implicit none
        
    contains
    
        !HR model
        subroutine hr_model(Ie, b, X_f, Y_f, Z_f, D_X, D_Y, D_Z)
            implicit none
            
            ! Declaration of variables:
            real(kind=8), intent(in) :: Ie, b, X_f, Y_f, Z_f
            real(kind=8), intent(out) :: D_X, D_Y, D_Z
            real(kind=8) :: a, c, d, r, s, Xr
            
            ! HR model parameters:
            a = 1.0d0
            c = 1.0d0
            d = 5.0d0
            r = 0.01d0
            s = 4.0d0
            Xr = -1.6d0
            
            ! HR equations:
            D_X = Y_f - a*(X_f**3.0d0) + b*(X_f**2.0d0) - Z_f + Ie
            D_Y = c - d*(X_f**2.0d0) - Y_f
            D_Z = r*(s*(X_f - Xr) - Z_f)
            
        end subroutine hr_model
        
        !4th order Runge-Kutta:
        subroutine rk4(h, Ie, b, X_o, Xf)
            implicit none
            
            !Declaration of variables:
            real(kind=8), intent(in) :: h, Ie, b
            real(kind=8), dimension(3), intent(in) :: X_o
            real(kind=8), dimension(3), intent(inout) :: Xf
            integer :: i, j
            real(kind=8), dimension(3) :: DX
            real(kind=8), dimension(4,3) :: k

            ! Iterate RK slopes:
            do i = 1, 4
                ! Calculate the HR equations (x = X(1), y = X(2), z = X(3)):
                call hr_model(Ie, b, Xf(1), Xf(2), Xf(3), DX(1), DX(2), DX(3))
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
        
        ! Iterate time and determine maximum and minimum ISI for a point (Ie,b)
        subroutine measure_ISI(Ie, b, max_ISI, min_ISI, n_ISIs)
            implicit none
            
            ! Declaration of variables:
            real(kind=8), intent(in) :: Ie, b
            integer, intent(out) :: n_ISIs
            real(kind=8), intent(out) :: max_ISI, min_ISI
            integer :: t1, t1_max, t1_trans, n_ISIs_proxy
            real(kind=8) :: rk_step, t_max, t_trans, t, t_isi, ISI
            real(kind=8), allocatable, dimension(:) :: V_ISI
            real(kind=8), dimension(3) :: Xo, X
            
            ! Discretize time:
            t_max = 10000.0d0 !Max time
            t_trans = 5000.0d0 !Transient
            rk_step = 0.01d0 !RK time step in the model time units
            t1_max = nint(t_max/rk_step) !Max number of time steps
            t1_trans = nint(t_trans/rk_step) !Transient number of time steps
            t = 0.0d0
            
            ! Vector to store the series of ISIs found in the time series:
            allocate(V_ISI(t1_max-t1_trans)) !Number of ISIs always <= to the time steps considered
            V_ISI = -1.0d0 
            
            n_ISIs = 0 !Count the number of ISIs in the time series
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
                call rk4(rk_step, Ie, b, Xo, X)
                !After transient:
                if (t >= t_trans) then
                    !Calculate the ISI:
                    if ((X(1)*Xo(1) < 0.0d0) .and. (X(1) < Xo(1))) then
                        if (t_isi == 0.0d0) then
                            t_isi = t
                        else if (t_isi > 0.0d0) then
                            n_ISIs = n_ISIs + 1
                            ISI = t - t_isi
                            t_isi = t
                            V_ISI(n_ISIs) = ISI
                        end if
                    end if
                end if
            end do

            ! Get maximum and minimum ISI in the series of ISIs:
            if (n_ISIs /= 0) then !Only when ISIs were found
                max_ISI = maxval(V_ISI)
                do n_ISIs_proxy = 1, (t1_max-t1_trans)
                    if (V_ISI(n_ISIs_proxy) == -1.0d0) then
                        V_ISI(n_ISIs_proxy) = max_ISI
                    end if
                end do
                min_ISI = minval(V_ISI)
            else
                max_ISI = 0.0d0
                min_ISI = 0.0d0
            end if
            
        end subroutine measure_ISI
        
    end module funcs_subrs
    
    !Main program:
    program HR_cell_phase_diagram_ISI
        use funcs_subrs
        implicit none

        !Calculates the phase diagram of the HR model. The HR model is solved with a
        !4th order Runge-Kutta method with a fixed time step. The interspike interval
        !(ISI) is measured in the time series for each point (Ie,b). Maximum and minimum
        !ISIs are obtained. The maximum ISI for each (Ie,b) is plotted associated
        !with a color scale
        
        ! Declaration of variables:
        integer :: n_Ie, n_Ie_max, n_b, n_b_max, number_ISIs
        real(kind=8) :: Ie_p, Ie_0, Ie_f, dIe, b_p, b_0, b_f, db, ISI_max, ISI_min
        
        ! Start and end of the parameter Ie range:
        Ie_0 = 1.0d0
        Ie_f = 6.0d0
        
        ! Start and end of the parameter b range:
        b_0 = 0.0d0 
        b_f = 5.0d0
        
        ! Interval between points:
        dIe = 0.1d0
        db = 0.1d0

        ! Number of points in the Ie axis of the diagram:
        n_Ie_max = nint((Ie_f-Ie_0)/dIe) + 1
        
        ! Number of points in the T axis of the diagram:
        n_b_max = nint((b_f-b_0)/db) + 1
        
        ! Open output files:
        open(unit=11, file="diag_ISI_max.dat", status="replace")
        open(unit=12, file="diag_ISI_min.dat", status="replace")
        
        !Iterate parameter Ie:
        Ie_p = Ie_0 - dIe
        do n_Ie = 1, n_Ie_max
            Ie_p = Ie_p + dIe
            ! Iterate paramater b:
            b_p = b_0 - db
            do n_b = 1, n_b_max
                b_p = b_p + db
                ! Determine maximum and minimum ISI in the time series:
                call measure_ISI(Ie_p, b_p, ISI_max, ISI_min, number_ISIs)
                ! Write to files (only when ISIs were found in the time series):
                if (number_ISIs /= 0) then
                    write(11,*) b_p, Ie_p, ISI_max
                    write(12,*) b_p, Ie_p, ISI_min
                end if
            end do
        end do
        
        ! CLose output files:
        close(unit=11)
        close(unit=12)
        
        ! Gnuplot scripts
        ! Plot the phase diagram with the maximum ISI in a png file:
        open(unit=13, file="diag_ISI_max.plt", status="replace")
            write(13, '(a)') "set terminal png enhanced size 1920,1080"
            write(13, '(a)') "set size square"
            write(13, '(a)') "unset key"
            write(13, '(a)') "set tics font ',20'"
            write(13, '(a)') "set cbrange[0:700]"
            write(13, '(a)') "set yrange[1.0:5.7]"
            write(13, '(a)') "set xrange[0.0:5.0]"
            write(13, '(a)') "set xlabel 'b' font ',30' offset '3,0,0'"
            write(13, '(a)') "set ylabel 'I_{e}' font ',30'"
            write(13, '(a)') "set cblabel 'ISI_{max}' font ',30'"
            write(13, '(a)') "set output 'diag_ISI_max.png'"
            write(13, '(a)') "plot 'diag_ISI_max.dat' using 1:2:3 with points pointtype 7 pointsize 4 palette"
        close(unit=13)

        call system("gnuplot -p diag_ISI_max.plt")
        
    end program HR_cell_phase_diagram_ISI