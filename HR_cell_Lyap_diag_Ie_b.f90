MODULE FUNCS_SUBRS
    IMPLICIT NONE

CONTAINS

    ! 4th order Runge-Kutta:
    SUBROUTINE RK4(N, X, Y, H, FCN, B, Ie)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(KIND=8), INTENT(IN) :: B, Ie      
        REAL(KIND=8), INTENT(INOUT) :: X, Y(N), H
        INTERFACE
            SUBROUTINE FCN(N, X, Y, YPRIME, B, Ie)
                INTEGER, INTENT(IN) :: N
                REAL(KIND=8), INTENT(IN) :: X, Y(N), B, Ie
                REAL(KIND=8), INTENT(OUT) :: YPRIME(N)
            END SUBROUTINE FCN
        END INTERFACE

        REAL(KIND=8) :: K1(N), K2(N), K3(N), K4(N), YTEMP(N)

        K1 = 0.0D0
        K2 = 0.0D0
        K3 = 0.0D0
        K4 = 0.0D0
        YTEMP = 0.0D0

        ! Calculate K1
        CALL FCN(N, X, Y, K1, B, Ie)

        ! Calculate K2
        YTEMP = Y + 0.5D0*H*K1
        CALL FCN(N, X + 0.5D0*H, YTEMP, K2, B, Ie)

        ! Calculate K3
        YTEMP = Y + 0.5D0*H*K2
        CALL FCN(N, X + 0.5D0*H, YTEMP, K3, B, Ie)

        ! Calculate K4
        YTEMP = Y + H*K3
        CALL FCN(N, X + H, YTEMP, K4, B, Ie)

        ! Update Y and X
        Y = Y + (H/6.0D0)*(K1 + 2.0D0*K2 + 2.0D0*K3 + K4)
        X = X + H
        
    END SUBROUTINE RK4

    ! Hindmarsh-Rose model and linearized equations
    SUBROUTINE FCN(N, X, Y, YPRIME, B, Ie)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(KIND=8), INTENT(IN) :: X, Y(N), B, Ie
        REAL(KIND=8), INTENT(OUT) :: YPRIME(N)
        INTEGER :: I
        REAL(KIND=8) :: A, C, D, R, S, Xr

        YPRIME = 0.0D0

        ! Hindmarsh-Rose parameters
        A = 1.0D0
        C = 1.0D0
        D = 5.0D0
        R = 0.05D0
        S = 4.0D0
        Xr = -1.6D0

        ! Hindmarsh-Rose equations of motion
        YPRIME(1) = Y(2) - A*Y(1)**3 + B*Y(1)**2 - Y(3) + Ie
        YPRIME(2) = C - D*Y(1)**2 - Y(2)
        YPRIME(3) = R*(S*(Y(1) - Xr) - Y(3))

        ! 3 COPIES OF LINEARIZED EQUATIONS OF MOTION
        !(product of the Jacobin matrix with a vector del = (del_X, del_Y, del_Z))        
        DO I = 0, 2
            YPRIME(4+I) = Y(7+I) - A*3.0D0*Y(1)**2 * Y(4+I) + 2.0D0*B*Y(1)*Y(4+I) - Y(10+I)
            YPRIME(7+I) = -D*2.0D0*Y(1)*Y(4+I) - Y(7+I)
            YPRIME(10+I) = R*S*Y(4+I) - R*Y(10+I)
        END DO
        
    END SUBROUTINE FCN

END MODULE FUNCS_SUBRS

! Main program
PROGRAM ODE
    USE FUNCS_SUBRS
    IMPLICIT NONE
    
    ! Calculates the Lyapunov spectrum for the Hindmarsh-Rose model using the Wolf method. Outputs the three exponents for different external current Ie and parameter b. The code was adapted from Fortran 77 directly from the Appendix A of the 1985 paper.

    ! Declaration of variables
    INTEGER, PARAMETER :: N = 3
    INTEGER, PARAMETER :: NN = 12
    INTEGER :: NSTEP, I, J, K, L, TRANS, NB, NB_MAX, NIe, NIe_MAX
    REAL(KIND=8) :: STPSZE, X, XEND, B, B0, BF, DB, Ie, Ie0, IeF, DIe, tic, toc
    REAL(KIND=8), DIMENSION(NN) :: Y
    REAL(KIND=8), DIMENSION(N) :: ZNORM, GSC, CUM

    CALL CPU_TIME(tic)

    ! Range of parameter b
    B0 = 1.3D0
    BF = 2.0D0
    DB = 0.0004D0
    NB_MAX = NINT((BF - B0) / DB)

    ! Range of parameter Ie
    Ie0 = 5.22D0
    IeF = 5.28D0
    DIe = 0.00004D0
    NIe_MAX = NINT((IeF - Ie0) / DIe)

    NSTEP = 1000000
    STPSZE = 0.01D0
    TRANS = 800000
    X = 0.0D0

    ! Open output files
    OPEN(UNIT = 10, FILE = 'lyapunov_spectrum_vs_b_Ie.dat', STATUS = 'REPLACE')

    ! Loop over parameter Ie
    Ie = Ie0 - DIe
    DO NIe = 1, NIe_MAX
        Ie = Ie + DIe

        ! Loop over parameter b
        B = B0 - DB
        DO NB = 1, NB_MAX
            B = B + DB

            ! Initialize the nonlinear system
            Y(1) = 0.0D0
            Y(2) = 0.0D0
            Y(3) = 0.0D0

            ! Initialize the linear system (orthonormal frame)
            DO I = N+1, NN
                Y(I) = 0.0D0
            END DO
            DO I = 1, N
                Y((N+1)*I) = 1.0D0
            END DO
            CUM = 0.0D0

            ! Main integration loop
            DO I = 1, NSTEP
                XEND = STPSZE * REAL(I, KIND=8)

                ! 4th Runge-Kutta integration
                CALL RK4(NN, X, Y, STPSZE, FCN, B, Ie)

                ! Gram-Schmidt orthonormalization

                ! Normalize the first vector
                ZNORM(1) = 0.0D0
                DO J = 1, N
                    ZNORM(1) = ZNORM(1) + Y(N*J+1)**2
                END DO
                ZNORM(1) = SQRT(ZNORM(1))
                DO J = 1, N
                    Y(N*J+1) = Y(N*J+1) / ZNORM(1)
                END DO

                ! Generate the new orthonormal set of vectors
                DO J = 2, N

                    ! Gram-Schmidt coefficients
                    DO K = 1, J - 1
                        GSC(K) = 0.0D0
                        DO L = 1, N
                            GSC(K) = GSC(K) + Y(N*L+J) * Y(N*L+K)
                        END DO
                    END DO

                    ! Construct the new vector
                    DO K = 1, N
                        DO L = 1, (J-1)
                            Y(N*K+J) = Y(N*K+J) - GSC(L) * Y(N*K+L)
                        END DO
                    END DO

                    ! Normalize
                    ZNORM(J) = 0.0D0
                    DO K = 1, N
                        ZNORM(J) = ZNORM(J) + Y(N*K+J)**2
                    END DO
                    ZNORM(J) = SQRT(ZNORM(J))

                    DO K = 1, N
                        Y(N*K+J) = Y(N*K+J) / ZNORM(J)
                    END DO

                END DO

                ! After transient
                IF (I >= TRANS) THEN

                    ! Accumulate the logarithms of the norms
                    DO K = 1, N
                        CUM(K) = CUM(K) + LOG(ZNORM(K))
                    END DO

                END IF

            END DO

            ! Calculate and write Lyapunov exponents for Ie and b
            WRITE(10, *) B, Ie, (CUM(J) / (NSTEP*STPSZE), J = 1, N)

        END DO

    END DO

    CLOSE(10)

    CALL CPU_TIME(toc)

    PRINT *, "Execution Time: ", toc - tic

END PROGRAM ODE
