    ! Functions and subroutines:
    module funcs_subrs
        implicit none
        
    contains
    
        !HR model
        subroutine hr_model(X_f, Y_f, Z_f, I_t, D_X, D_Y, D_Z)
            implicit none
            
            ! Declaration of variables:
            real(kind=8), intent(in) :: X_f, Y_f, Z_f, I_t
            real(kind=8), intent(out) :: D_X, D_Y, D_Z  
            real(kind=8) :: a, b, c, d, r, s, Xr, Ie
            
            ! HR model parameters:
            a = 1.0d0
            b = 1.0d0
            c = 1.0d0
            d = 5.0d0
            r = 0.05d0
            s = 4.0d0
            Xr = -1.6d0
            Ie = 5.2885d0
            
            ! HR equations:
            D_X = Y_f - a*(X_f**3.0d0) + b*(X_f**2.0d0) - Z_f + Ie + I_t
            D_Y = c - d*(X_f**2.0d0) - Y_f
            D_Z = r*(s*(X_f - Xr) - Z_f)
            
        end subroutine hr_model
        
        !4th order Runge-Kutta:
        subroutine rk4(h, X_o, I_t, Xf)
            implicit none
            
            ! Declaration of variables:
            real(kind=8), intent(in) :: h
            real(kind=8), dimension(3), intent(in) :: X_o
            real(kind=8), dimension(3), intent(inout) :: Xf
            integer :: i, j
            real(kind=8) :: I_t
            real(kind=8), dimension(3) :: DX
            real(kind=8), dimension(4,3) :: k

            ! Iterate RK slopes:
            do i = 1, 4
                ! Calculate the HR equations (x = X(1), y = X(2), z = X(3)):
                call hr_model(Xf(1), Xf(2), Xf(3), I_t, DX(1), DX(2), DX(3))
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
        
    end module funcs_subrs
    
    !Main program:
    program HR_RK4
        use funcs_subrs
        implicit none
        
        !Simulates a single Hindmarsh-Rose cell using the 4th order Runge-Kutta method
        !with a fixed time step. An external periodic current pulse of adjustable duration
        !is also implemented. The membrane potential is used to approximate the
        !electrocardiogram (ECG). Outputs the membrane potential x(t), external current I(t),
        !phase space and ECG(t). 
        
        ! Declaration of variables:
        integer :: t1, t1_max, St, dur, dt, t_trans, D_t1
        real(kind=8) :: rk_step, t_max, t, It, I_stim, S_c, S_e, dist, ECG
        real(kind=8), allocatable, dimension(:) :: Vm
        real(kind=8), dimension(3) :: Xo, X
        
        ! Discretize time:
        t_max = 5000.0d0 !Max time
        rk_step = 0.01d0 !RK time step in the model time units
        t1_max = nint(t_max/rk_step) !Max number of time steps
        t = 0.0d0

        ! Stimulus parameters:
        I_stim = 0.0d0 !Intensity of the pulse
        dur = 5000 !Duration of the pulse (in discrete time steps, i.e., dur = [dur in the model time units] / h)
        dt = 300000 !Interval between pulses (dt = [period in the model time units] / h)
        t_trans = dt !Transient before stimulation (in discrete time steps)
        
        ! ECG parameters:
        S_c = 0.0d0 !Position of the cell
        S_e = 1.0d0 !Position of the electrode
        dist = sqrt((S_c-S_e)**2) !Distance between cell and electrode
        allocate(Vm(t1_max))
        
        ! Open output files:
        open(unit=11, file="potential.dat", status="replace")
        open(unit=12, file="current.dat", status="replace")     
        open(unit=13, file="phase_space.dat", status="replace")
        
        ! Iterate time:
        do t1 = 1, t1_max
            t = t + rk_step
            ! Initial conditions (x = X(1), y = X(2), z = X(3)):
            if (t1 == 1) then
                St = 0 !Count the time steps for stimulation
                Xo(1) = 0.0d0
                Xo(2) = 0.0d0
                Xo(3) = 0.0d0
                X = 0.0d0
            ! Update variables:
            else if (t1 > 1) then 
                Xo = X          
            end if
            ! Periodic stimulation:
            if (t1 >= t_trans) then
                St = St + 1
                It = 0.0d0
                if (St <= dur) then
                    It = I_stim
                end if
                if (St == dt) then
                    St = 0
                end if
            end if
            ! 4th order Runge-Kutta with HR model:
            call rk4(rk_step, Xo, It, X)
            ! Write to file the 1st variable x(t):
            write(11,*) t, X(1)
            write(12,*) t, It
            ! Save the membrane potential to calculate the ECG/LFP:
            Vm(t1) = X(1)           
            ! Write to file the phase space (x,y,z):
            if (t >= 1000.0d0) then !Ignore transient
                write(13,*) X(1), X(2), X(3)
            end if
        end do
        
        ! CLose output files:
        close(unit=11)
        close(unit=12)
        close(unit=13)

        !Calculate the ECG:
        open(unit=14, file="ECG.dat", status="replace")
        t = 0.0d0
        D_t1 = 1 !D_t1 = D_t/h for Vm(t+D_t)-Vm(t), where D_t1=1 => D_t=h (minimum)
        do t1 = 1, t1_max-D_t1
           t = t + rk_step
           ECG = Vm(t1+D_t1)-Vm(t1)
           write(14,*) t, ECG
        end do
        close(unit=14)

        ! Gnuplot scripts
        ! Plot the potential x(t):
        open(unit=15, file="potential.plt", status="replace")
            write(15, '(a)') "set term wxt"
            write(15, '(a)') "set tics font ', 12'"
            write(15, '(a)') "set xlabel 't' font ', 16'"
            write(15, '(a)') "set ylabel 'x' font ', 16'"
            write(15, '(a)') "set grid"
            write(15, '(a)') "plot 'potential.dat' w lines lt 6 lw 1.5"
    !       write(15, '(a)') "replot 'potential.dat'" !Use to plot the discretized result
        close(unit=15)
        call system("gnuplot -p potential.plt")
        
        ! Plot the phase space (x,y,z):
        open(unit=16, file="phase_space.plt", status="replace")
            write(16, '(a)') "set term wxt"     
            write(16, '(a)') "set xlabel 'x' font ', 16'"
            write(16, '(a)') "set ylabel 'y' font ', 16'"
            write(16, '(a)') "set zlabel 'z' font ', 16'"
            write(16, '(a)') "set grid"         
            write(16, '(a)') "splot 'phase_space.dat' w lines lt 6 lw 1.5"
        close(unit=16)      
        call system("gnuplot -p phase_space.plt")

        ! Plot the electrocardiogram ECG(t):
        open(unit=17, file="ECG.plt", status="replace")
            write(17, '(a)') "set term wxt"
            write(17, '(a)') "set tics font ', 12'"
            write(17, '(a)') "set xlabel 't' font ', 16'"
            write(17, '(a)') "set ylabel 'ECG' font ', 16'"
            write(17, '(a)') "set grid"
            write(17, '(a)') "plot 'ECG.dat' w lines lt 6 lw 1.5"
        close(unit=17)
        call system("gnuplot -p ECG.plt")
        
    end program HR_RK4