    ! Functions and subroutines:
    module funcs_subrs
        implicit none
        
    contains
 
        !KTz model and Jacobian
        subroutine ktz_tanh_jacob(Xr, T, X_o, Y_o, Z_o, time, trans, X_f, Y_f, Z_f, Jacob)
            implicit none
        
			! Declaration of variables:
            real(kind=8), intent(in) :: Xr, T, X_o, Y_o, Z_o
			integer, intent(in) :: time, trans
            real(kind=8), intent(out) :: X_f, Y_f, Z_f
			real(kind=8), dimension(3,3), intent(out) :: Jacob
			real(kind=8) :: arg, k, del, lamb, Ie, a11, a12
    
			! KTz model parameters:
			k = 0.6d0
			del = 0.001d0
			lamb = 0.001d0
			Ie = 0.0d0
	
			! KTz equations:
            arg = (X_o - k * Y_o + Z_o + Ie) / T
            X_f = tanh(arg) !Hyperbolic KTz
            Y_f = X_o
            Z_f = (1.0d0 - del) * Z_o - lamb * (X_o - Xr)

			! Calculate the Jacobian:
			Jacob = 0.0d0			
		    if (time >= trans) then !Only after transient
				! Elements of the Jacobian matrix:
		        a11 = (1.0d0 - X_f**2.0d0) / T
		        a12 = -k * a11
				! Jacobian matrix:
				Jacob(1,1) = a11
				Jacob(1,2) = a12
				Jacob(1,3) = a11
				Jacob(2,1) = 1.0d0
				Jacob(3,1) = -lamb
				Jacob(3,3) = 1.0d0 - del
		    end	if
			
        end subroutine ktz_tanh_jacob
		
		!LU decomposition (similar to SciPy and Matlab lu() function)
		subroutine lu(M_input, O_1, O_2)
			implicit none
			
			real(kind=8), dimension(3,3), intent(in) :: M_input
			real(kind=8), dimension(3,3), intent(out) :: O_1, O_2
			integer :: u, info
			real(kind=8), dimension(3) :: row
			real(kind=8), dimension(3,3) :: Low, Per
			integer, dimension(3) :: ipiv
			external :: DGEMM, DGETRF2 !External subroutines for matrix product and LU decomposition
			
			! LAPACK subroutine that receives a matrix in Low, decomposes it into lower
			! and upper triangular (O_2) and returns the results in Low, ignoring the 1s
			! on the diagonal of the lower. ipiv is a vector with the indices of the
			! permuted rows:
			Low = M_input
			call DGETRF2(3,3,Low,3,ipiv,info)

			!if (info /= 0) stop "lufact: info /= 0"
			
			! Removes the upper triangular (O_2) matrix from Low and restores the 1s
			! to the diagonal in Low. It also creates an identity matrix in Per
			! for the permutations:
			O_2 = 0.0d0
			Per = 0.0d0
			do u = 1, 3
				O_2(u,u:3) = Low(u,u:3)
				Low(u,u:3) = 0.0d0
				Low(u,u) = 1.0d0
				Per(u,u) = 1.0d0
			end do
			
			! Use the ipiv vector to rewrite the identity as the permutation matrix (Per):
			do u = 1, 3
				row = Per(u,:)
				Per(u,:) = Per(ipiv(u),:)
				Per(ipiv(u),:) = row
			end do
			
			! Calculate the remaining (permuted lower) matrix O_1 from the LU decomposition as
			! O_1 = Per_transpose*Lower_triangular
			call DGEMM('T','N',3,3,3,1.0d0,Per,3,Low,3,0.0d0,O_1,3)
			
		end subroutine lu
            
		!Iterate time and determine the Lyapunov spectrum for a point (Xr,T)
		subroutine lyapunov(Xr, T, exp_lyap)
			implicit none
			
			! Declaration of variables:
			real(kind=8), intent(in) :: Xr, T
			real(kind=8), dimension(3), intent(out) :: exp_lyap
			integer :: t1, t_max, t_trans, u
			real(kind=8) :: Xo, Yo, Zo, X, Y, Z
			real(kind=8), dimension(3) :: diag_O2, sum_log_diag
			real(kind=8), dimension(3,3) :: M, M_dgemm, O1, O2
			external :: DGEMM !External subroutine for matrix poduct
			
			! Additional parameters:
			t_max = 100000 !Total time steps
			t_trans = 50000 !Transient
			
			sum_log_diag = 0.0d0
			
			! Iterate time:
			do t1 = 1, t_max
				! Initial conditions:
				if (t1 == 1) then
					Xo = 0.0d0
					Yo = 0.0d0
					Zo = 0.0d0
				! Update variables:
				else if (t1 > 1) then
					Xo = X
					Yo = Y
					Zo = Z
				end if
				! Calculate the KTz equations and the Jacobian:
				call ktz_tanh_jacob(Xr, T, Xo, Yo, Zo, t1, t_trans, X, Y, Z, M);
				! Calculate the Lyapunov exponents using the Eckmann-Ruelle method:
				if (t1 >= t_trans) then
					! Multiply the Jacobian by the remaining matrix from the LU
					! decomposition:
					M_dgemm = M
					if (t1 > t_trans) then
						! DGEMM for matrix product (faster than matmul):
						! M = M*O1
						call DGEMM('N','N',3,3,3,1.0d0,M_dgemm,3,O1,3,0.0d0,M,3)
					end if
					! LU decomposition (retuns the upper triangular matrix O2 and
					! the remaining matrix O1):
					call lu(M, O1, O2)
					! Get the diagonal elements of the upper triangular matrix:
					do u = 1, 3
						diag_O2(u) = O2(u,u)
					end do
					! Store in a sum the log of each diagonal element:
					sum_log_diag = sum_log_diag + log(dabs(diag_O2))
				end if
			end do
			
			! Divide the sum of the diagonal elements for the the total time steps
			! (+ 1 because it includes the transient):
			exp_lyap = sum_log_diag / (t_max - t_trans + 1)*1.0d0			  
			
		end subroutine lyapunov
		
    end module funcs_subrs
    
    ! Main program:
    program ktz_tanh_lyap_eck_diag_m
        use funcs_subrs
        implicit none
        
		!Calculates the Lyapunov exponent over the phase diagram of the KTz model.
		!The Lyapunov exponent is calculate using the Eckmann-Ruelle method for
		!each point (Xr,T), where the Jacobian matrix is triangularized by LU  
		!decompostion at each iteration using Lapack/MKL subroutines.
		!Outputs the main exponent for each (Xr,T).

        ! Declaration of variables:
        integer :: n_T, n_T_max, n_Xr, n_Xr_max
        real(kind=8) :: T_p, T_0, T_f, Xr_p, Xr_0, Xr_f, d_T, dXr
		real(kind=8), dimension(3) :: lyap_exp
		
		!Start and end of the parameter T range:
		T_0 = 1.0e-8
		T_f = 0.65d0
		
		!Start and end of the parameter Xr range:
		Xr_0 = -0.35d0
		Xr_f = -0.35d0
		
		!Interval between points:
		d_T = 0.001d0
		dXr = 0.001d0
		  
		! Number of points in the T axis of the diagram:
		n_T_max = nint((T_f-T_0)/d_T) + 1

		! Number of points in the Xr axis of the diagram:
		n_Xr_max = nint((Xr_f-Xr_0)/dXr) + 1

		! Open output file:
		open(unit=10, file="exp_lyap.dat", status="replace")

		! Iterate paramater Xr:
		Xr_p = Xr_0 - dXr
		do n_Xr = 1, n_Xr_max
			Xr_p = Xr_p + dXr
			! Iterate parameter T:
			T_p = T_0 - d_T	
			do n_T = 1, n_T_max
				T_p = T_p + d_T
				! Obtain the Lyapunov spectrum for (Xr,T):
				call lyapunov(Xr_p, T_p, lyap_exp)
				! Write to file the main Lyapunov exponent:
				!write(10,*) Xr_p, T_p, lyap_exp(1)
				! Spectrum vs T:
				write(10,*) T_p, lyap_exp(:)
			end do
		end do

        ! Close output files:
        close(unit=10)
		
    end program ktz_tanh_lyap_eck_diag_m