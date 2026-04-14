! @@@@@

MODULE MODAL_DG_1D_GLOBALS

  IMPLICIT NONE

  ! COUNTERS
  INTEGER                                                                          :: ORD                                    ! INDEX OF THE DEGREE FOR LEGENDRE POLYNOMIALS
  INTEGER                                                                          :: I                                      ! SPATIAL (X-AXIS) INDEX
  INTEGER                                                                          :: N                                      ! TIME INDEX
  INTEGER                                                                          :: I0                                     ! GRID REFINEMENT INDEX

  ! REFINEMENT CONTROL
  INTEGER, PARAMETER                                                               :: NLEV       = 6                         ! NUMBER OF GRID LEVELS

  ! ORDER OF THE DISCONTINUOUS GALERKIN METHOD AND POLYNOMIALS
  INTEGER, PARAMETER                                                               :: DG_ORD     = 3                         ! ORDER OF DG SOLVER
  INTEGER, PARAMETER                                                               :: POLY_ORD   = DG_ORD - 1                ! ORDER OF ADAPTIVE POLYNOMIALS
  INTEGER, PARAMETER                                                               :: COMMON_ORD = 4                         ! ORDER OF COMMON POLYNOMIALS
  INTEGER, PARAMETER                                                               :: NQ         = 6                         ! NUMBER OF QUADRATURE POINTS
  INTEGER, PARAMETER                                                               :: LAMBDA_SW  = 1                         ! 0 : EXACT / 1 : COMPUTED FOR SIN( PI * X )
  INTEGER, PARAMETER                                                               :: WENO_SW    = 1                         ! 0 : WITHOUT WENO  / 1 : WITH WENO  

  ! USEFUL CONSTANTS AND FUNCTIONS
  REAL( 8 ), PARAMETER                                                             :: PI  = 4.0D0 * DATAN( 1.0D0 )
  REAL( 8 ), PARAMETER                                                             :: CFL = 0.4D0 / DBLE( 2 * POLY_ORD + 1 ) ! CFL NUMBER

  ! SPATIAL VARIABLES
  INTEGER                                                                          :: IMAX                                   ! ACTIVE NUMBER OF SPATIAL CELLS
  INTEGER,   PARAMETER                                                             :: I_BC     = 1                           ! GHOST BOUNDARY CELLS
  REAL( 8 ), PARAMETER                                                             :: X_LEFT   = - 1.0D0                     ! LEFT BOUNDARY ON X-AXIS
  REAL( 8 ), PARAMETER                                                             :: X_RIGHT  = + 1.0D0                     ! RIGHT BOUNDARY ON X-AXIS
  REAL( 8 ), PARAMETER                                                             :: X_LENGTH = DABS( X_RIGHT - X_LEFT )    ! LENGTH OF DOMAIN
  INTEGER,   PARAMETER                                                             :: IMAX_MAX = 40 * ( 2 ** NLEV )


  REAL( 8 ), DIMENSION( - I_BC + 1 : IMAX_MAX + I_BC )                             :: X                                      ! CELL CENTERS
  REAL( 8 ), DIMENSION( - I_BC + 1 : IMAX_MAX + I_BC )                             :: DX                                     ! CELL VOLUMES
  REAL( 8 )                                                                        :: MAXDX                                  ! MAXIMAL CELL VOLUME

  ! TIME VARIABLES
  INTEGER,   PARAMETER                                                             :: RK_ORD    = 4                          ! ORDER OF RUNGE-KUTTA SOLVER
  INTEGER,   PARAMETER                                                             :: NMAX      = 99999                      ! NUMBER OF THE MARCHING STEP
  REAL( 8 ), PARAMETER                                                             :: TIME_OUT  = 2.0D0                      ! FINAL TIME OF THE SOLUTION
  REAL( 8 )                                                                        :: DT                                     ! TIME STEP SIZE
  REAL( 8 )                                                                        :: FLAG_0, FLAG_1                         ! TIME FLAGS

  ! CPU TIME CHECKERS
  INTEGER                                                                          :: RATE, TIME_0, TIME_1

  ! VARIABLES FOR QUADRATURES
  REAL( 8 ), DIMENSION( 1 : NQ, 1 : 2 )                                            :: Q

  ! MASS COEFFICIENT MATRICES
  REAL( 8 ), DIMENSION( 0 : COMMON_ORD, 0 : COMMON_ORD )                           :: M
  REAL( 8 ), DIMENSION( 0 : POLY_ORD, 0 : POLY_ORD, - I_BC + 1 : IMAX_MAX + I_BC ) :: M_STAGE

  ! CELL INTERFACE CONTRIBUTIONS
  REAL( 8 ), DIMENSION( 0 : POLY_ORD, - I_BC + 1 : IMAX_MAX + I_BC )               :: L_B, R_B

  ! CELL-WISE SHAPE PARAMETERS
  REAL( 8 ), DIMENSION( - I_BC + 1 : IMAX_MAX + I_BC )                             :: LAMBDA

  ! DEGREES OF FREEDOM IN THE MOMENTS BY THE L2 PROJECTION
  REAL( 8 ), DIMENSION( 0 : COMMON_ORD, - I_BC + 1 : IMAX_MAX + I_BC )             :: OLD_DEG
  REAL( 8 ), DIMENSION( 0 : COMMON_ORD, 1 : IMAX_MAX )                             :: NEW_DEG

  ! CONVEX SUMMATION OF DEGREES OF FREEDOM AND LEGENDRE POLYNOMIALS
  REAL( 8 )                                                                        :: TEMP
  REAL( 8 ), DIMENSION( 1 : IMAX_MAX )                                             :: U, NEW_U

  ! ERROR STORAGE FOR CONVERGENCE STUDY
  REAL( 8 ), DIMENSION( 0 : NLEV )                                                 :: L2_ERR, LINF_ERR

  ! PLOTTING VARIABLES
  INTEGER, PARAMETER                                                               :: WIDTH = 36
  INTEGER                                                                          :: PAD
  INTEGER                                                                          :: LEN_ROW
  CHARACTER( LEN = 80 )                                                            :: ROW

END MODULE MODAL_DG_1D_GLOBALS

! @@@@@

PROGRAM MODAL_DG_1D

  USE MODAL_DG_1D_GLOBALS

  IMPLICIT NONE

  PRINT *, "================================================================================================="
  PRINT *, "                            SCALAR CONSERVATION LAWS SOLVER. (ADVECTION)                         "
  PRINT *, "================================================================================================="
  PRINT *, "                                    CALCULATIONS HAVE STARTED.                                   "
  PRINT *, "================================================================================================="
  PRINT *, " "
  
  CALL GET_QUADRATURE( Q )
  
  PRINT *, "------------------------------------"
  PRINT *, "CFL NUMBER:", REAL( CFL )
  PRINT *, "LAMBDA  SWITCH:", LAMBDA_SW
  PRINT *, "LIMITER SWITCH:", WENO_SW
  PRINT *, "------------------------------------"
  PRINT *, " "

  CALL SYSTEM_CLOCK( COUNT_RATE = RATE )
  CALL SYSTEM_CLOCK( TIME_0 )

  DO I0 = 0, NLEV

      PRINT *, "------------------------------------"
      PRINT *, "GRID LEVEL: ", I0
      IMAX =   ( 2 ** ( I0  ) ) * 4 
      PRINT *, "NUMBER OF CELLS: ", IMAX

      FLAG_0 = 0.0D0
      FLAG_1 = FLAG_0

      CALL GET_SPATIAL( I_BC, IMAX, X_LEFT, X_LENGTH, X, DX, MAXDX )
      CALL GET_MASS   ( COMMON_ORD, Q, M )
      CALL GET_INITIAL( PI, COMMON_ORD, Q, M, I_BC, IMAX, X, DX, OLD_DEG )

      N = 1

      DO WHILE ( ( FLAG_0 < TIME_OUT ) .AND. ( N <= NMAX ) )

          FLAG_1 = FLAG_0

          ! COMPUTE THE TIME STEP SIZE BY USING THE FIRST-DERIVATIVE OF FLUX
          CALL GET_TIMESTEP( I_BC, IMAX, DX, FLAG_1, TIME_OUT, CFL, DT )

          ! UPDATE THE VARIABLE BY USING STAGE-WISE ADAPTIVE RK SOLVER
          CALL GET_RK( RK_ORD, COMMON_ORD, POLY_ORD, Q, M, I_BC, IMAX, X, DX, DT, OLD_DEG, NEW_DEG )

          ! UPDATE THE REAL-TIME
          FLAG_1 = FLAG_1 + DT
          FLAG_0 = FLAG_1

          ! RECEIVE THE UPDATED DATA
          DO I = 1, IMAX

              DO ORD = 0, COMMON_ORD

                  OLD_DEG( ORD, I ) = NEW_DEG( ORD, I )

              END DO

          END DO

          N = N + 1

      END DO

      DO I = 1, IMAX

          TEMP = 0.0D0

          DO ORD = 0, COMMON_ORD

              TEMP = TEMP + OLD_DEG( ORD, I ) * POLY_C( ORD, 0.0D0 )

          END DO

          NEW_U( I ) = TEMP

      END DO

      ! SET UP THE EXACT SOLUTION AT TIME_OUT FOR THIS GRID
      DO I = 1, IMAX

          U( I ) = EXT_U( X( I ), PI, TIME_OUT )

      END DO

      ! COMPUTE VARIOUS ERRORS AND STORE THEM FOR CONVERGENCE STUDY
      PRINT *, "--------- ERRORS (FINAL) -----------"

      CALL GET_L2    ( I_BC, IMAX, DX, U, NEW_U, L2_ERR  ( I0 ) )
      CALL GET_LINFTY(       IMAX,     U, NEW_U, LINF_ERR( I0 ) )

      PRINT '(A18, ES12.4E3)', " L_2 ERROR:     ", L2_ERR( I0 )
      PRINT '(A18, ES12.4E3)', " L_INFTY ERROR: ", LINF_ERR( I0 )
      PRINT *, "------------------------------------"
      PRINT *, " "

  END DO

  CALL SYSTEM_CLOCK( TIME_1 )

  PRINT *, "---------- L^INFTY ORDERS ----------"
  PRINT '(A, 1X, F12.3)', "1ST STAGE ORDER: ", REAL( LOG2( LINF_ERR( 0 ) / LINF_ERR( 1 ) ) )
  PRINT '(A, 1X, F12.3)', "2ND STAGE ORDER: ", REAL( LOG2( LINF_ERR( 1 ) / LINF_ERR( 2 ) ) )
  PRINT '(A, 1X, F12.3)', "3RD STAGE ORDER: ", REAL( LOG2( LINF_ERR( 2 ) / LINF_ERR( 3 ) ) )
  PRINT '(A, 1X, F12.3)', "4TH STAGE ORDER: ", REAL( LOG2( LINF_ERR( 3 ) / LINF_ERR( 4 ) ) )
  PRINT '(A, 1X, F12.3)', "5TH STAGE ORDER: ", REAL( LOG2( LINF_ERR( 4 ) / LINF_ERR( 5 ) ) )
  PRINT '(A, 1X, F12.3)', "6TH STAGE ORDER: ", REAL( LOG2( LINF_ERR( 5 ) / LINF_ERR( 6 ) ) )
  PRINT *, "------------------------------------"

  PRINT *, "-------------- L^2 ORDERS ----------"
  PRINT '(A, 1X, F12.3)', "1ST STAGE ORDER: ", REAL( LOG2( L2_ERR( 0 ) / L2_ERR( 1 ) ) )
  PRINT '(A, 1X, F12.3)', "2ND STAGE ORDER: ", REAL( LOG2( L2_ERR( 1 ) / L2_ERR( 2 ) ) )
  PRINT '(A, 1X, F12.3)', "3RD STAGE ORDER: ", REAL( LOG2( L2_ERR( 2 ) / L2_ERR( 3 ) ) )
  PRINT '(A, 1X, F12.3)', "4TH STAGE ORDER: ", REAL( LOG2( L2_ERR( 3 ) / L2_ERR( 4 ) ) )
  PRINT '(A, 1X, F12.3)', "5TH STAGE ORDER: ", REAL( LOG2( L2_ERR( 4 ) / L2_ERR( 5 ) ) )
  PRINT '(A, 1X, F12.3)', "6TH STAGE ORDER: ", REAL( LOG2( L2_ERR( 5 ) / L2_ERR( 6 ) ) )
  PRINT *, "------------------------------------"

  PRINT *, "------------------------------------"
  PRINT '(A, F12.3, 1X, A)', "COMPUTATIONAL TIME : ", REAL( TIME_1 - TIME_0 ) / REAL( RATE ), "SEC"
  PRINT *, "------------------------------------"
  PRINT *, " "

  PRINT *, "================================================================================================="
  PRINT *, "                                   CALCULATIONS HAVE FINISHED.                                   "
  PRINT *, "================================================================================================="

  STOP

CONTAINS

  ! @@@@@

  REAL( 8 ) FUNCTION POLY_C( ORD, X )

    IMPLICIT NONE

    INTEGER,   INTENT( IN ) :: ORD
    REAL( 8 ), INTENT( IN ) :: X
    !----------------------------------------- CALCULATIONS HAVE STARTED ----------------------------------------!

    ! SET UP COMMON ORTHONORMAL BASES
    IF ( ORD .EQ. 0 ) THEN

        POLY_C = 1.0D0

    ELSEIF ( ORD .EQ. 1 ) THEN

        POLY_C = X

    ELSEIF ( ORD .EQ. 2 ) THEN

        POLY_C = X ** 2 - 1.0D0 / 12.0D0

    ELSEIF ( ORD .EQ. 3 ) THEN

        POLY_C = X ** 3 - ( 3.0D0 / 20.0D0 ) * X

    ELSEIF ( ORD .EQ. 4 ) THEN

        POLY_C = X ** 4 - ( 3.0D0 / 14.0D0 ) * X ** 2 + 3.0D0 / 560.0D0

    END IF

    RETURN
    !----------------------------------------- CALCULATIONS HAVE FINISHED ---------------------------------------!
  END FUNCTION POLY_C

  REAL( 8 ) FUNCTION D_POLY_C( ORD, X )

    IMPLICIT NONE

    INTEGER,   INTENT( IN ) :: ORD
    REAL( 8 ), INTENT( IN ) :: X
    !----------------------------------------- CALCULATIONS HAVE STARTED ----------------------------------------!

    ! SET UP DERIVATIVES OF COMMON ORTHONORMAL BASES
    IF ( ORD .EQ. 0 ) THEN

        D_POLY_C = 0.0D0

    ELSEIF ( ORD .EQ. 1 ) THEN

        D_POLY_C = 1.0D0

    ELSEIF ( ORD .EQ. 2 ) THEN

        D_POLY_C = 2.0D0 * X

    ELSEIF ( ORD .EQ. 3 ) THEN

        D_POLY_C = 3.0D0 * X ** 2 - 3.0D0 / 20.0D0

    ELSEIF ( ORD .EQ. 4 ) THEN

        D_POLY_C = 4.0D0 * X ** 3 - ( 3.0D0 / 7.0D0 ) * X

    END IF

    RETURN
    !----------------------------------------- CALCULATIONS HAVE FINISHED ---------------------------------------!
  END FUNCTION D_POLY_C

  REAL( 8 ) FUNCTION D2_POLY_C( ORD, X )

    IMPLICIT NONE

    INTEGER,   INTENT( IN ) :: ORD
    REAL( 8 ), INTENT( IN ) :: X
    !----------------------------------------- CALCULATIONS HAVE STARTED ----------------------------------------!

    ! SET UP SECOND DERIVATIVES OF COMMON ORTHONORMAL BASES
    IF ( ORD .EQ. 0 ) THEN

        D2_POLY_C = 0.0D0

    ELSEIF ( ORD .EQ. 1 ) THEN

        D2_POLY_C = 0.0D0

    ELSEIF ( ORD .EQ. 2 ) THEN

        D2_POLY_C = 2.0D0

    ELSEIF ( ORD .EQ. 3 ) THEN

        D2_POLY_C = 6.0D0 * X

    ELSEIF ( ORD .EQ. 4 ) THEN

        D2_POLY_C = 12.0D0 * X ** 2 - 3.0D0 / 7.0D0

    END IF

    RETURN
    !----------------------------------------- CALCULATIONS HAVE FINISHED ---------------------------------------!
  END FUNCTION D2_POLY_C

  REAL( 8 ) FUNCTION D3_POLY_C( ORD, X )

    IMPLICIT NONE

    INTEGER,   INTENT( IN ) :: ORD
    REAL( 8 ), INTENT( IN ) :: X
    !----------------------------------------- CALCULATIONS HAVE STARTED ----------------------------------------!

    ! SET UP THIRD DERIVATIVES OF COMMON ORTHONORMAL BASES
    IF ( ORD .EQ. 0 ) THEN

        D3_POLY_C = 0.0D0

    ELSEIF ( ORD .EQ. 1 ) THEN

        D3_POLY_C = 0.0D0

    ELSEIF ( ORD .EQ. 2 ) THEN

        D3_POLY_C = 0.0D0

    ELSEIF ( ORD .EQ. 3 ) THEN

        D3_POLY_C = 6.0D0

    ELSEIF ( ORD .EQ. 4 ) THEN

        D3_POLY_C = 24.0D0 * X

    END IF

    RETURN
    !----------------------------------------- CALCULATIONS HAVE FINISHED ---------------------------------------!
  END FUNCTION D3_POLY_C

  ! @@@@@

  ! @@@@@

  REAL( 8 ) FUNCTION POLY( ORD, X, LAMBDA )

    IMPLICIT NONE

    INTEGER,   INTENT( IN ) :: ORD
    REAL( 8 ), INTENT( IN ) :: X, LAMBDA
    !----------------------------------------- CALCULATIONS HAVE STARTED ----------------------------------------!

    ! SET UP STAGE-WISE ADAPTIVE BASES
    IF ( ORD .EQ. 0 ) THEN

        POLY = 1.0D0

    ELSEIF ( ORD .EQ. 1 ) THEN

        POLY = X + ( 2.0D0 / 3.0D0 ) * ( LAMBDA ) * ( X ** 3 )

    ELSEIF ( ORD .EQ. 2 ) THEN

        POLY = ( X ** 2 )                                  &
             + ( 1.0D0 / 3.0D0 ) * ( LAMBDA ) * ( X ** 4 ) &
             - ( 1.0D0 / 12.0D0 )                          &
             - ( LAMBDA ) / 240.0D0

    END IF

    RETURN
    !----------------------------------------- CALCULATIONS HAVE FINISHED ---------------------------------------!
  END FUNCTION POLY

  REAL( 8 ) FUNCTION D_POLY( ORD, X, LAMBDA )

    IMPLICIT NONE

    INTEGER,   INTENT( IN ) :: ORD
    REAL( 8 ), INTENT( IN ) :: X, LAMBDA
    !----------------------------------------- CALCULATIONS HAVE STARTED ----------------------------------------!

    ! SET UP DERIVATIVES OF STAGE-WISE ADAPTIVE BASES
    IF ( ORD .EQ. 0 ) THEN

        D_POLY = 0.0D0

    ELSEIF ( ORD .EQ. 1 ) THEN

        D_POLY = 1.0D0 + 2.0D0 * ( LAMBDA ) * ( X ** 2 )

    ELSEIF ( ORD .EQ. 2 ) THEN

        D_POLY = 2.0D0 * X + ( 4.0D0 / 3.0D0 ) * ( LAMBDA ) * ( X ** 3 )

    END IF

    RETURN
    !----------------------------------------- CALCULATIONS HAVE FINISHED ---------------------------------------!
  END FUNCTION D_POLY

  REAL( 8 ) FUNCTION D2_POLY( ORD, X, LAMBDA )

    IMPLICIT NONE

    INTEGER,   INTENT( IN ) :: ORD
    REAL( 8 ), INTENT( IN ) :: X, LAMBDA
    !----------------------------------------- CALCULATIONS HAVE STARTED ----------------------------------------!

    ! SET UP SECOND DERIVATIVES OF STAGE-WISE ADAPTIVE BASES
    IF ( ORD .EQ. 0 ) THEN

        D2_POLY = 0.0D0

    ELSEIF ( ORD .EQ. 1 ) THEN

        D2_POLY = 4.0D0 * ( LAMBDA ) * X

    ELSEIF ( ORD .EQ. 2 ) THEN

        D2_POLY = 2.0D0 + 4.0D0 * ( LAMBDA ) * ( X ** 2 )

    END IF

    RETURN
    !----------------------------------------- CALCULATIONS HAVE FINISHED ---------------------------------------!
  END FUNCTION D2_POLY

  ! @@@@@

  ! @@@@@

  REAL( 8 ) FUNCTION FLUX( U )

    IMPLICIT NONE

    REAL( 8 ), INTENT( IN ) :: U
    !----------------------------------------- CALCULATIONS HAVE STARTED ----------------------------------------!

    ! SET UP THE FLUX
    FLUX = U

    RETURN
    !----------------------------------------- CALCULATIONS HAVE FINISHED ---------------------------------------!
  END FUNCTION FLUX

  SUBROUTINE GET_QUADRATURE( Q )

    IMPLICIT NONE

    REAL( 8 ), INTENT( OUT ) :: Q( 1 : NQ, 1 : 2 )
    !----------------------------------------- CALCULATIONS HAVE STARTED ----------------------------------------!

    ! ! 6TH-POINT GAUSS-LEGENDRE QUADRATURE ON [ -1.0D0 / 2.0D0, 1.0D0 / 2.0D0 ]
    Q( 1, 1 ) = + 0.4662347571015760D0
    Q( 2, 1 ) = + 0.3306046932331323D0
    Q( 3, 1 ) = + 0.1193095930415985D0
    Q( 4, 1 ) = - Q( 3, 1 )
    Q( 5, 1 ) = - Q( 2, 1 )
    Q( 6, 1 ) = - Q( 1, 1 )

    Q( 1, 2 ) = 0.0856622461895852D0
    Q( 2, 2 ) = 0.1803807865240693D0
    Q( 3, 2 ) = 0.2339569672863455D0
    Q( 4, 2 ) = Q( 3, 2 )
    Q( 5, 2 ) = Q( 2, 2 )
    Q( 6, 2 ) = Q( 1, 2 )

    RETURN
    !----------------------------------------- CALCULATIONS HAVE FINISHED ---------------------------------------!
  END SUBROUTINE GET_QUADRATURE
 
  SUBROUTINE GET_MASS( MAX_ORD, Q, M )

    IMPLICIT NONE

    INTEGER,   INTENT( IN )  :: MAX_ORD
    REAL( 8 ), INTENT( IN )  :: Q( 1 : NQ, 1 : 2 )

    REAL( 8 ), INTENT( OUT ) :: M( 0 : MAX_ORD, 0 : MAX_ORD )

    ! COUNTERS
    INTEGER                  :: ORD_0, ORD_1, QUAD
    !----------------------------------------- CALCULATIONS HAVE STARTED ----------------------------------------!

    ! COMPUTE THE INVERSE OF THE COMMON MASS MATRIX
    DO ORD_0 = 0, MAX_ORD

        DO ORD_1 = 0, MAX_ORD

            M( ORD_0, ORD_1 ) = 0.0D0

            DO QUAD = 1, NQ

                M( ORD_0, ORD_1 ) = M( ORD_0, ORD_1 ) &
                                  + POLY_C( ORD_0, Q( QUAD, 1 ) ) &
                                  * POLY_C( ORD_1, Q( QUAD, 1 ) ) * Q( QUAD, 2 )

            END DO

            IF ( M( ORD_0, ORD_1 ) .GE. 1.0D-12 ) THEN

                M( ORD_0, ORD_1 ) = 1.0D0 / M( ORD_0, ORD_1 )

            ELSE

                M( ORD_0, ORD_1 ) = 0.0D0

            END IF

        END DO

    END DO

    RETURN
    !----------------------------------------- CALCULATIONS HAVE FINISHED ---------------------------------------!
  END SUBROUTINE GET_MASS

  SUBROUTINE GET_SPATIAL( I_BC, IMAX, X_LEFT, X_LENGTH, X, DX, MAXDX )

    IMPLICIT NONE

    INTEGER,   INTENT( IN )  :: I_BC, IMAX
    REAL( 8 ), INTENT( IN )  :: X_LEFT, X_LENGTH

    REAL( 8 ), INTENT( OUT ) :: X    ( - I_BC + 1 : IMAX + I_BC )
    REAL( 8 ), INTENT( OUT ) :: DX   ( - I_BC + 1 : IMAX + I_BC )
    REAL( 8 ), INTENT( OUT ) :: MAXDX

    ! COUNTER
    INTEGER                  :: I

    ! CALCULATION VARIABLE
    REAL( 8 ), DIMENSION( 0 : IMAX ) :: TEMP_X
    !----------------------------------------- CALCULATIONS HAVE STARTED ----------------------------------------!

    ! GENERATE PHYSICAL CELL CENTER LOCATIONS AND VOLUMES OF CELLS
    DO I = 0, IMAX

        TEMP_X( I ) = X_LEFT + DBLE( I ) * DBLE( X_LENGTH / IMAX )

    END DO

    DO I = 1, IMAX

        X ( I ) = ( TEMP_X( I ) + TEMP_X( I - 1 ) ) / 2.0D0
        DX( I ) = DABS( TEMP_X( I ) - TEMP_X( I - 1 ) )

    END DO

    MAXDX = MAXVAL( DX( 1 : IMAX ) )

    ! SET UP GHOST LOCATIONS AND VOLUMES OF CELLS
    DO I = 1, I_BC

        X ( - I + 1 ) = X( IMAX - I + 1 ) - X_LENGTH
        X ( IMAX + I ) = X( I )           + X_LENGTH

        DX( - I + 1 ) = DX( IMAX - I + 1 )
        DX( IMAX + I ) = DX( I )

    END DO

    RETURN
    !----------------------------------------- CALCULATIONS HAVE FINISHED ---------------------------------------!
  END SUBROUTINE GET_SPATIAL

  SUBROUTINE GET_BOUNDARY( MAX_ORD, I_BC, IMAX, DEG )

    IMPLICIT NONE

    INTEGER, INTENT( IN )      :: MAX_ORD
    INTEGER, INTENT( IN )      :: I_BC, IMAX
    REAL( 8 ), INTENT( INOUT ) :: DEG( 0 : MAX_ORD, - I_BC + 1 : IMAX + I_BC )

    ! COUNTERS
    INTEGER                    :: ORD
    INTEGER                    :: I
    !----------------------------------------- CALCULATIONS HAVE STARTED ----------------------------------------!

    ! APPLY PERIODIC CONDITIONS TO GET GHOST CELLS
    DO I = 1, I_BC

        DO ORD = 0, MAX_ORD

            DEG( ORD, - I + 1 ) = DEG( ORD, IMAX - I + 1 )
            DEG( ORD, IMAX + I ) = DEG( ORD, I )

        END DO

    END DO

    RETURN
    !----------------------------------------- CALCULATIONS HAVE FINISHED ---------------------------------------!
  END SUBROUTINE GET_BOUNDARY

  ! U_0가 사실은 DOF 자리에 있음. 
  SUBROUTINE GET_INITIAL( PI, MAX_ORD, Q, M, I_BC, IMAX, X, DX, U_0 )

    IMPLICIT NONE

    REAL( 8 ), INTENT( IN )  :: PI
    INTEGER,   INTENT( IN )  :: MAX_ORD
    INTEGER,   INTENT( IN )  :: I_BC, IMAX
    REAL( 8 ), INTENT( IN )  :: Q ( 1 : NQ, 1 : 2 )
    REAL( 8 ), INTENT( IN )  :: M ( 0 : MAX_ORD, 0 : MAX_ORD )
    REAL( 8 ), INTENT( IN )  :: X ( - I_BC + 1 : IMAX + I_BC )
    REAL( 8 ), INTENT( IN )  :: DX( - I_BC + 1 : IMAX + I_BC )

    REAL( 8 ), INTENT( OUT ) :: U_0( 0 : MAX_ORD, - I_BC + 1 : IMAX + I_BC )

    ! COUNTERS
    INTEGER                  :: ORD
    INTEGER                  :: I
    INTEGER                  :: QUAD
    !----------------------------------------- CALCULATIONS HAVE STARTED ----------------------------------------!

    ! COMPUTE L2-PROJECTION OF DATA INTO THE COMMON FINITE ELEMENT SPACE
    DO I = 1, IMAX

        DO ORD = 0, MAX_ORD

            U_0( ORD, I ) = 0.0D0

            DO QUAD = 1, NQ

                U_0( ORD, I ) = U_0( ORD, I ) &
                              + INI_U( X( I ) + Q( QUAD, 1 ) * DX( I ), PI ) &
                              * POLY_C( ORD, Q( QUAD, 1 ) ) * Q( QUAD, 2 )

            END DO

            U_0( ORD, I ) = M( ORD, ORD ) * U_0( ORD, I )

        END DO

    END DO

    CALL GET_BOUNDARY( MAX_ORD, I_BC, IMAX, U_0 )

    RETURN
    !----------------------------------------- CALCULATIONS HAVE FINISHED ---------------------------------------!
  END SUBROUTINE GET_INITIAL

  ! 람다를 설정하고, 그 각 셀의 람다에 따른 mass matrix를 만드는 서브루틴 
  ! 경계에서의 값도 업데이트 함. 
  SUBROUTINE GET_STAGE_DATA( POLY_ORD, COMMON_ORD, Q, I_BC, IMAX, X, DX, COM_DEG, LAMBDA, M_STAGE, L_B, R_B )

    IMPLICIT NONE

    INTEGER,   INTENT( IN )    :: POLY_ORD, COMMON_ORD
    INTEGER,   INTENT( IN )    :: I_BC, IMAX
    REAL( 8 ), INTENT( IN )    :: Q      ( 1 : NQ, 1 : 2 )
    REAL( 8 ), INTENT( IN )    :: X      ( - I_BC + 1 : IMAX + I_BC )
    REAL( 8 ), INTENT( IN )    :: DX     ( - I_BC + 1 : IMAX + I_BC )
    REAL( 8 ), INTENT( INOUT ) :: COM_DEG( 0 : COMMON_ORD, - I_BC + 1 : IMAX + I_BC )

    REAL( 8 ), INTENT( OUT )   :: LAMBDA ( - I_BC + 1 : IMAX + I_BC )
    REAL( 8 ), INTENT( OUT )   :: M_STAGE( 0 : POLY_ORD, 0 : POLY_ORD, - I_BC + 1 : IMAX + I_BC )
    REAL( 8 ), INTENT( OUT )   :: L_B    ( 0 : POLY_ORD, - I_BC + 1 : IMAX + I_BC )
    REAL( 8 ), INTENT( OUT )   :: R_B    ( 0 : POLY_ORD, - I_BC + 1 : IMAX + I_BC )

    ! COUNTERS
    INTEGER                    :: I, ORD, ORD_0, ORD_1, QUAD
    !----------------------------------------- CALCULATIONS HAVE STARTED ----------------------------------------!

    IF ( LAMBDA_SW .EQ. 0 ) THEN

        CALL GET_EXACT_LAMBDA( I_BC, IMAX, DX, LAMBDA )

    ELSEIF ( LAMBDA_SW .EQ. 1 ) THEN

        CALL GET_LAMBDA( COMMON_ORD, I_BC, IMAX, X, DX, COM_DEG, LAMBDA )

    END IF

    DO I = 1, IMAX

        DO ORD_0 = 0, POLY_ORD

            DO ORD_1 = 0, POLY_ORD

                M_STAGE( ORD_0, ORD_1, I ) = 0.0D0

                DO QUAD = 1, NQ

                    M_STAGE( ORD_0, ORD_1, I ) = M_STAGE( ORD_0, ORD_1, I ) &
                                                + POLY( ORD_0, Q( QUAD, 1 ), LAMBDA( I ) ) &
                                                * POLY( ORD_1, Q( QUAD, 1 ), LAMBDA( I ) ) * Q( QUAD, 2 )

                END DO

                IF ( M_STAGE( ORD_0, ORD_1, I ) .GE. 1.0D-12 ) THEN

                    M_STAGE( ORD_0, ORD_1, I ) = 1.0D0 / M_STAGE( ORD_0, ORD_1, I )

                ELSE

                    M_STAGE( ORD_0, ORD_1, I ) = 0.0D0

                END IF

            END DO

        END DO

        DO ORD = 0, POLY_ORD

            L_B( ORD, I ) = POLY( ORD, - 1.0D0 / 2.0D0, LAMBDA( I ) )
            R_B( ORD, I ) = POLY( ORD, + 1.0D0 / 2.0D0, LAMBDA( I ) )

        END DO

    END DO

    DO I = 1, I_BC

        L_B( :, - I + 1 )  = L_B( :, IMAX - I + 1 )
        R_B( :, - I + 1 )  = R_B( :, IMAX - I + 1 )
        L_B( :, IMAX + I ) = L_B( :, I )
        R_B( :, IMAX + I ) = R_B( :, I )

        M_STAGE( :, :, - I + 1 )  = M_STAGE( :, :, IMAX - I + 1 )
        M_STAGE( :, :, IMAX + I ) = M_STAGE( :, :, I )

    END DO

    RETURN
    !----------------------------------------- CALCULATIONS HAVE FINISHED ---------------------------------------!
  END SUBROUTINE GET_STAGE_DATA

  ! space를 변화시키는 서브루틴
  ! 다항식 기저 --> 비다항식 기저(계수를 변화시키는 것) 
  SUBROUTINE PROJECT_COMMON_TO_ADAPTIVE( COMMON_ORD, POLY_ORD, Q, M_STAGE, I_BC, IMAX, COM_DEG, LAMBDA, ADP_DEG )

    IMPLICIT NONE

    INTEGER,   INTENT( IN )  :: COMMON_ORD, POLY_ORD
    INTEGER,   INTENT( IN )  :: I_BC, IMAX
    REAL( 8 ), INTENT( IN )  :: Q      ( 1 : NQ, 1 : 2 )
    REAL( 8 ), INTENT( IN )  :: M_STAGE( 0 : POLY_ORD, 0 : POLY_ORD, - I_BC + 1 : IMAX + I_BC )
    REAL( 8 ), INTENT( IN )  :: COM_DEG( 0 : COMMON_ORD, - I_BC + 1 : IMAX + I_BC )
    REAL( 8 ), INTENT( IN )  :: LAMBDA ( - I_BC + 1 : IMAX + I_BC )

    REAL( 8 ), INTENT( OUT ) :: ADP_DEG( 0 : POLY_ORD, - I_BC + 1 : IMAX + I_BC )

    ! COUNTERS
    INTEGER                  :: I, ORD, K, QUAD

    ! CALCULATION VARIABLES
    REAL( 8 )                :: U_VAL
    !----------------------------------------- CALCULATIONS HAVE STARTED ----------------------------------------!

    DO I = 1, IMAX

        DO ORD = 0, POLY_ORD

            ADP_DEG( ORD, I ) = 0.0D0

            DO QUAD = 1, NQ

                U_VAL = 0.0D0

                DO K = 0, COMMON_ORD

                    U_VAL = U_VAL + COM_DEG( K, I ) * POLY_C( K, Q( QUAD, 1 ) )

                END DO

                ADP_DEG( ORD, I ) = ADP_DEG( ORD, I ) &
                                  + U_VAL * POLY( ORD, Q( QUAD, 1 ), LAMBDA( I ) ) * Q( QUAD, 2 )

            END DO

            ADP_DEG( ORD, I ) = M_STAGE( ORD, ORD, I ) * ADP_DEG( ORD, I )

        END DO

    END DO

    CALL GET_BOUNDARY( POLY_ORD, I_BC, IMAX, ADP_DEG )

    RETURN
    !----------------------------------------- CALCULATIONS HAVE FINISHED ---------------------------------------!
  END SUBROUTINE PROJECT_COMMON_TO_ADAPTIVE

  !비디항식 기저--> 다항식 기저
  SUBROUTINE PROJECT_ADAPTIVE_TO_COMMON( COMMON_ORD, POLY_ORD, Q, M, I_BC, IMAX, ADP_DEG, LAMBDA, COM_DEG )

    IMPLICIT NONE

    INTEGER,   INTENT( IN )  :: COMMON_ORD, POLY_ORD
    INTEGER,   INTENT( IN )  :: I_BC, IMAX
    REAL( 8 ), INTENT( IN )  :: Q      ( 1 : NQ, 1 : 2 )
    REAL( 8 ), INTENT( IN )  :: M      ( 0 : COMMON_ORD, 0 : COMMON_ORD )
    REAL( 8 ), INTENT( IN )  :: ADP_DEG( 0 : POLY_ORD, - I_BC + 1 : IMAX + I_BC )
    REAL( 8 ), INTENT( IN )  :: LAMBDA ( - I_BC + 1 : IMAX + I_BC )

    REAL( 8 ), INTENT( OUT ) :: COM_DEG( 0 : COMMON_ORD, 1 : IMAX )

    ! COUNTERS
    INTEGER                  :: I, ORD, K, QUAD

    ! CALCULATION VARIABLES
    REAL( 8 )                :: U_VAL
    !----------------------------------------- CALCULATIONS HAVE STARTED ----------------------------------------!

    DO I = 1, IMAX

        DO K = 0, COMMON_ORD

            COM_DEG( K, I ) = 0.0D0

            DO QUAD = 1, NQ

                U_VAL = 0.0D0

                DO ORD = 0, POLY_ORD

                    U_VAL = U_VAL + ADP_DEG( ORD, I ) * POLY( ORD, Q( QUAD, 1 ), LAMBDA( I ) )

                END DO

                COM_DEG( K, I ) = COM_DEG( K, I ) &
                                + U_VAL * POLY_C( K, Q( QUAD, 1 ) ) * Q( QUAD, 2 )

            END DO

            COM_DEG( K, I ) = M( K, K ) * COM_DEG( K, I )

        END DO

    END DO

    RETURN
    !----------------------------------------- CALCULATIONS HAVE FINISHED ---------------------------------------!
  END SUBROUTINE PROJECT_ADAPTIVE_TO_COMMON

  ! 람다 만들고, 그에 따른 경계값 생성 + 질량 행열 생성
  ! 원래 가지고 있던 르잔드르 기저 공간에서 비 다항식 기저 공간으로 사영 
  ! 비 다항식 공간에서 리미터 적용 이후 flux 계산하기 
  ! 비 다항식 공간에서 다항식 공간으로 사영 
  SUBROUTINE GET_STAGE_RHS( COMMON_ORD, POLY_ORD, Q, M, I_BC, IMAX, X, DX, COM_DEG, RHS_DEG )

    IMPLICIT NONE

    INTEGER,   INTENT( IN )    :: COMMON_ORD, POLY_ORD
    INTEGER,   INTENT( IN )    :: I_BC, IMAX
    REAL( 8 ), INTENT( IN )    :: Q      ( 1 : NQ, 1 : 2 )
    REAL( 8 ), INTENT( IN )    :: M      ( 0 : COMMON_ORD, 0 : COMMON_ORD )
    REAL( 8 ), INTENT( IN )    :: X      ( - I_BC + 1 : IMAX + I_BC )
    REAL( 8 ), INTENT( IN )    :: DX     ( - I_BC + 1 : IMAX + I_BC )
    REAL( 8 ), INTENT( INOUT ) :: COM_DEG( 0 : COMMON_ORD, - I_BC + 1 : IMAX + I_BC )

    REAL( 8 ), INTENT( OUT )   :: RHS_DEG( 0 : COMMON_ORD, 1 : IMAX )

    ! CALCULATION VARIABLES
    REAL( 8 ), DIMENSION( 0 : POLY_ORD, - I_BC + 1 : IMAX + I_BC ) :: ADP_DEG
    REAL( 8 ), DIMENSION( 0 : POLY_ORD, 1 : IMAX )                 :: RHS_ADP
    !----------------------------------------- CALCULATIONS HAVE STARTED ----------------------------------------!

    ! DETERMINE THE CELL-WISE STAGE DATA
    CALL GET_STAGE_DATA( POLY_ORD, COMMON_ORD, Q, I_BC, IMAX, X, DX, COM_DEG, LAMBDA, M_STAGE, L_B, R_B )

    ! REPROJECT THE COMMON STAGE STATE ONTO THE ADAPTIVE STAGE BASIS
    CALL PROJECT_COMMON_TO_ADAPTIVE( COMMON_ORD, POLY_ORD, Q, M_STAGE, I_BC, IMAX, COM_DEG, LAMBDA, ADP_DEG )

    ! 여기서 RHS_ADP 이거 뭐지? 
    ! COMPUTE THE ADAPTIVE RESIDUAL
    CALL GET_RESIDUAL( POLY_ORD, Q, M_STAGE, L_B, R_B, I_BC, IMAX, DX, LAMBDA, ADP_DEG, RHS_ADP )
    
    ! 다시 돌아가야하는 이유가 있을까? 
    ! REPROJECT THE ADAPTIVE RESIDUAL BACK TO THE COMMON BASIS
    CALL PROJECT_ADAPTIVE_TO_COMMON( COMMON_ORD, POLY_ORD, Q, M, I_BC, IMAX, &
    ADP_DEG = RESHAPE_ADP( RHS_ADP, I_BC, IMAX, POLY_ORD ), LAMBDA = LAMBDA, COM_DEG = RHS_DEG )

    RETURN
    !----------------------------------------- CALCULATIONS HAVE FINISHED ---------------------------------------!
  END SUBROUTINE GET_STAGE_RHS

  ! 람다, 질량행열, 경계값을 업데이트하고 그 업데이트 한 것을 다항식 기저에 반영한다. 
  SUBROUTINE GET_ACCEPTED_STATE( COMMON_ORD, POLY_ORD, Q, M, I_BC, IMAX, X, DX, COM_DEG )

    IMPLICIT NONE

    INTEGER,   INTENT( IN )    :: COMMON_ORD, POLY_ORD
    INTEGER,   INTENT( IN )    :: I_BC, IMAX
    REAL( 8 ), INTENT( IN )    :: Q      ( 1 : NQ, 1 : 2 )
    REAL( 8 ), INTENT( IN )    :: M      ( 0 : COMMON_ORD, 0 : COMMON_ORD )
    REAL( 8 ), INTENT( IN )    :: X      ( - I_BC + 1 : IMAX + I_BC )
    REAL( 8 ), INTENT( IN )    :: DX     ( - I_BC + 1 : IMAX + I_BC )

    REAL( 8 ), INTENT( INOUT ) :: COM_DEG( 0 : COMMON_ORD, - I_BC + 1 : IMAX + I_BC )

    ! COUNTERS
    INTEGER                                                        :: ORD
    INTEGER                                                        :: I

    ! CALCULATION VARIABLES
    REAL( 8 ), DIMENSION( 0 : POLY_ORD, - I_BC + 1 : IMAX + I_BC ) :: ADP_DEG
    REAL( 8 ), DIMENSION( 0 : COMMON_ORD, 1 : IMAX )               :: TEMP_COM
    !----------------------------------------- CALCULATIONS HAVE STARTED ----------------------------------------!

    ! DETERMINE THE CELL-WISE ACCEPTED BASIS
    CALL GET_STAGE_DATA( POLY_ORD, COMMON_ORD, Q, I_BC, IMAX, X, DX, COM_DEG, LAMBDA, M_STAGE, L_B, R_B )

    ! 이건 왜 하는 걸까? --> 추측하건데 마지막에 람다가 업데이트 되었을 때 그 때 결론을 얻기 위해서가 아닐까? 
    ! REPROJECT THE ACCEPTED STATE ONTO THE ADAPTIVE SPACE
    CALL PROJECT_COMMON_TO_ADAPTIVE( COMMON_ORD, POLY_ORD, Q, M_STAGE, I_BC, IMAX, COM_DEG, LAMBDA, ADP_DEG )

    ! STORE THE ACCEPTED ADAPTIVE STATE BACK IN THE COMMON COORDINATES
    CALL PROJECT_ADAPTIVE_TO_COMMON( COMMON_ORD, POLY_ORD, Q, M, I_BC, IMAX, ADP_DEG, LAMBDA, TEMP_COM )
    ! 이거 왜 위에다가 com_deg넣는 거 안하고 이렇게 따로 하는 걸까? 
    DO I = 1, IMAX

        DO ORD = 0, COMMON_ORD

            COM_DEG( ORD, I ) = TEMP_COM( ORD, I )

        END DO

    END DO

    CALL GET_BOUNDARY( COMMON_ORD, I_BC, IMAX, COM_DEG )

    RETURN
    !----------------------------------------- CALCULATIONS HAVE FINISHED ---------------------------------------!
  END SUBROUTINE GET_ACCEPTED_STATE

  ! get_ boundary 랑 다른 게 없는 듯..?  
  FUNCTION RESHAPE_ADP( RHS_ADP, I_BC, IMAX, POLY_ORD ) RESULT( ADP_OUT )

    IMPLICIT NONE

    INTEGER,   INTENT( IN ) :: I_BC, IMAX, POLY_ORD
    REAL( 8 ), INTENT( IN ) :: RHS_ADP( 0 : POLY_ORD, 1 : IMAX )

    REAL( 8 )               :: ADP_OUT( 0 : POLY_ORD, - I_BC + 1 : IMAX + I_BC )

    INTEGER                 :: ORD, I
    !----------------------------------------- CALCULATIONS HAVE STARTED ----------------------------------------!

    ADP_OUT( :, : ) = 0.0D0

    DO I = 1, IMAX

        DO ORD = 0, POLY_ORD

            ADP_OUT( ORD, I ) = RHS_ADP( ORD, I )

        END DO

    END DO

    DO I = 1, I_BC

        ADP_OUT( :,  - I + 1 ) = ADP_OUT( :, IMAX - I + 1 )
        ADP_OUT( :, IMAX + I ) = ADP_OUT( :, I )

    END DO

    RETURN
    !----------------------------------------- CALCULATIONS HAVE FINISHED ---------------------------------------!
  END FUNCTION RESHAPE_ADP

  ! rk 하는데 람다 업데이트를 매번 step마다 하고, 
  ! 마지막에도 dof가 바뀌었을 것이니, 한 번 더 해서 새로운 dof 를 만드는 코드(시간에 따라)
  SUBROUTINE GET_RK( RK_ORD, COMMON_ORD, POLY_ORD, Q, M, I_BC, IMAX, X, DX, DT, OLD_DEG, NEW_DEG )

    IMPLICIT NONE

    INTEGER,   INTENT( IN )    :: RK_ORD
    INTEGER,   INTENT( IN )    :: COMMON_ORD
    INTEGER,   INTENT( IN )    :: POLY_ORD
    INTEGER,   INTENT( IN )    :: I_BC, IMAX
    REAL( 8 ), INTENT( IN )    :: Q      ( 1 : NQ, 1 : 2 )
    REAL( 8 ), INTENT( IN )    :: M      ( 0 : COMMON_ORD, 0 : COMMON_ORD )
    REAL( 8 ), INTENT( IN )    :: X      ( - I_BC + 1 : IMAX + I_BC )
    REAL( 8 ), INTENT( IN )    :: DX     ( - I_BC + 1 : IMAX + I_BC )
    REAL( 8 ), INTENT( IN )    :: DT

    REAL( 8 ), INTENT( INOUT ) :: OLD_DEG( 0 : COMMON_ORD, - I_BC + 1 : IMAX + I_BC )
    REAL( 8 ), INTENT( OUT )   :: NEW_DEG( 0 : COMMON_ORD, 1 : IMAX )

    ! COUNTERS
    INTEGER                                                                 :: ORD
    INTEGER                                                                 :: I

    ! CALCULATION VARIABLES
    REAL( 8 ), DIMENSION( 0 : COMMON_ORD, - I_BC + 1 : IMAX + I_BC, 1 : 3 ) :: TEMP_DEG
    REAL( 8 ), DIMENSION( 0 : COMMON_ORD, 1 : IMAX )                        :: F_FLUX
    REAL( 8 ), DIMENSION( 0 : COMMON_ORD, - I_BC + 1 : IMAX + I_BC )        :: ACC_DEG
    !----------------------------------------- CALCULATIONS HAVE STARTED ----------------------------------------!

    IF ( RK_ORD .EQ. 3 ) THEN

        ! 3RD-ORDER TVD RK SOLVER
        !==========================================================================================================!
        ! STAGE 1
        !==========================================================================================================!

        CALL GET_STAGE_RHS( COMMON_ORD, POLY_ORD, Q, M, I_BC, IMAX, X, DX, OLD_DEG, F_FLUX )

        DO I = 1, IMAX

            DO ORD = 0, COMMON_ORD

              TEMP_DEG( ORD, I, 1 ) = + 1.0D0 * OLD_DEG( ORD, I ) &
                                      + 1.0D0 * DT * F_FLUX( ORD, I )

            END DO

        END DO

        CALL GET_BOUNDARY( COMMON_ORD, I_BC, IMAX, TEMP_DEG( :, :, 1 ) )

        !==========================================================================================================!
        ! STAGE 2
        !==========================================================================================================!

        CALL GET_STAGE_RHS( COMMON_ORD, POLY_ORD, Q, M, I_BC, IMAX, X, DX, TEMP_DEG( :, :, 1 ), F_FLUX )

        DO I = 1, IMAX

            DO ORD = 0, COMMON_ORD

              TEMP_DEG( ORD, I, 2 ) = + 3.0D0 * OLD_DEG( ORD, I ) &
                                      + 1.0D0 * TEMP_DEG( ORD, I, 1 ) &
                                      + 1.0D0 * DT * F_FLUX( ORD, I )

              TEMP_DEG( ORD, I, 2 ) = TEMP_DEG( ORD, I, 2 ) / 4.0D0

            END DO

        END DO

        CALL GET_BOUNDARY( COMMON_ORD, I_BC, IMAX, TEMP_DEG( :, :, 2 ) )

        !==========================================================================================================!
        ! STAGE 3
        !==========================================================================================================!

        CALL GET_STAGE_RHS( COMMON_ORD, POLY_ORD, Q, M, I_BC, IMAX, X, DX, TEMP_DEG( :, :, 2 ), F_FLUX )

        DO I = 1, IMAX

            DO ORD = 0, COMMON_ORD

                NEW_DEG( ORD, I ) = + 1.0D0 * OLD_DEG( ORD, I ) &
                                    + 2.0D0 * TEMP_DEG( ORD, I, 2 ) &
                                    + 2.0D0 * DT * F_FLUX( ORD, I )

                NEW_DEG( ORD, I ) = NEW_DEG( ORD, I ) / 3.0D0

            END DO

        END DO

    ELSEIF ( RK_ORD .EQ. 4 ) THEN

        ! 4TH-ORDER NON-TVD RK SOLVER
        !==========================================================================================================!
        ! STAGE 1
        !==========================================================================================================!

        CALL GET_STAGE_RHS( COMMON_ORD, POLY_ORD, Q, M, I_BC, IMAX, X, DX, OLD_DEG, F_FLUX )

        DO I = 1, IMAX

            DO ORD = 0, COMMON_ORD

                TEMP_DEG( ORD, I, 1 ) = + 1.0D0 * OLD_DEG( ORD, I ) &
                                        + 0.5D0 * DT * F_FLUX( ORD, I )

            END DO

        END DO

        CALL GET_BOUNDARY( COMMON_ORD, I_BC, IMAX, TEMP_DEG( :, :, 1 ) )

        !==========================================================================================================!
        ! STAGE 2
        !==========================================================================================================!

        CALL GET_STAGE_RHS( COMMON_ORD, POLY_ORD, Q, M, I_BC, IMAX, X, DX, TEMP_DEG( :, :, 1 ), F_FLUX )

        DO I = 1, IMAX

            DO ORD = 0, COMMON_ORD

                TEMP_DEG( ORD, I, 2 ) = + 1.0D0 * OLD_DEG( ORD, I ) &
                                        + ( DT / 2.0D0 ) * F_FLUX( ORD, I )

            END DO

        END DO

        CALL GET_BOUNDARY( COMMON_ORD, I_BC, IMAX, TEMP_DEG( :, :, 2 ) )

        !==========================================================================================================!
        ! STAGE 3
        !==========================================================================================================!

        CALL GET_STAGE_RHS( COMMON_ORD, POLY_ORD, Q, M, I_BC, IMAX, X, DX, TEMP_DEG( :, :, 2 ), F_FLUX )

        DO I = 1, IMAX

            DO ORD = 0, COMMON_ORD

                TEMP_DEG( ORD, I, 3 ) = + 1.0D0 * OLD_DEG( ORD, I ) &
                                        + 1.0D0 * DT * F_FLUX( ORD, I )

            END DO

        END DO

        CALL GET_BOUNDARY( COMMON_ORD, I_BC, IMAX, TEMP_DEG( :, :, 3 ) )

        !==========================================================================================================!
        ! STAGE 4
        !==========================================================================================================!

        CALL GET_STAGE_RHS( COMMON_ORD, POLY_ORD, Q, M, I_BC, IMAX, X, DX, TEMP_DEG( :, :, 3 ), F_FLUX )

        DO I = 1, IMAX

            DO ORD = 0, COMMON_ORD

                NEW_DEG( ORD, I ) = - 1.0D0 * OLD_DEG( ORD, I ) &
                                    + 1.0D0 * TEMP_DEG( ORD, I, 1 ) &
                                    + 2.0D0 * TEMP_DEG( ORD, I, 2 ) &
                                    + 1.0D0 * TEMP_DEG( ORD, I, 3 ) &
                                    + ( DT / 2.0D0 ) * F_FLUX( ORD, I )

                NEW_DEG( ORD, I ) = NEW_DEG( ORD, I ) / 3.0D0

            END DO

        END DO

    END IF

    DO I = 1, IMAX

        DO ORD = 0, COMMON_ORD

            ACC_DEG( ORD, I ) = NEW_DEG( ORD, I )

        END DO

    END DO

    CALL GET_BOUNDARY( COMMON_ORD, I_BC, IMAX, ACC_DEG )
    CALL GET_ACCEPTED_STATE( COMMON_ORD, POLY_ORD, Q, M, I_BC, IMAX, X, DX, ACC_DEG )

    DO I = 1, IMAX

        DO ORD = 0, COMMON_ORD

            NEW_DEG( ORD, I ) = ACC_DEG( ORD, I )

        END DO

    END DO

    RETURN
    !----------------------------------------- CALCULATIONS HAVE FINISHED ---------------------------------------!
  END SUBROUTINE GET_RK
  
  ! @@@@@

  ! @@@@@

  REAL( 8 ) FUNCTION LOG2( X )

    IMPLICIT NONE

    REAL( 8 ), INTENT( IN ) :: X
    !----------------------------------------- CALCULATIONS HAVE STARTED ----------------------------------------!

    LOG2 = LOG( X ) / LOG( 2.0D0 )

    RETURN
    !----------------------------------------- CALCULATIONS HAVE FINISHED ---------------------------------------!
  END FUNCTION LOG2

  SUBROUTINE GET_L2( I_BC, IMAX, DX, U_EXACT, U_NUM, ERR_L2 )

    IMPLICIT NONE

    INTEGER,                                          INTENT( IN )  :: I_BC, IMAX
    REAL( 8 ), DIMENSION( - I_BC + 1 : IMAX + I_BC ), INTENT( IN )  :: DX
    REAL( 8 ), DIMENSION( 1 : IMAX ),                 INTENT( IN )  :: U_EXACT, U_NUM

    REAL( 8 ),                                        INTENT( OUT ) :: ERR_L2

    INTEGER                                                         :: I, ORD, QUAD
    REAL( 8 )                                                       :: SUM
    REAL( 8 )                                                       :: XQ
    REAL( 8 )                                                       :: UH
    REAL( 8 )                                                       :: UE

    REAL( 8 ), DIMENSION( 1 : 6, 1 : 2 )                            :: Q_ERR
    !----------------------------------------- CALCULATIONS HAVE STARTED ----------------------------------------!

    ! 6-POINT GAUSS-LEGENDRE QUADRATURE ON [ -1.0D0 / 2.0D0, 1.0D0 / 2.0D0 ]
    Q_ERR( 1, 1 ) = + 0.4662347571015760D0
    Q_ERR( 2, 1 ) = + 0.3306046932331323D0
    Q_ERR( 3, 1 ) = + 0.1193095930415985D0
    Q_ERR( 4, 1 ) = - Q_ERR( 3, 1 )
    Q_ERR( 5, 1 ) = - Q_ERR( 2, 1 )
    Q_ERR( 6, 1 ) = - Q_ERR( 1, 1 )

    Q_ERR( 1, 2 ) = 0.0856622461895852D0
    Q_ERR( 2, 2 ) = 0.1803807865240693D0
    Q_ERR( 3, 2 ) = 0.2339569672863455D0
    Q_ERR( 4, 2 ) = Q_ERR( 3, 2 )
    Q_ERR( 5, 2 ) = Q_ERR( 2, 2 )
    Q_ERR( 6, 2 ) = Q_ERR( 1, 2 )

    SUM = 0.0D0

    DO I = 1, IMAX

        DO QUAD = 1, 6

            XQ = X( I ) + Q_ERR( QUAD, 1 ) * DX( I )

            UH = 0.0D0
            DO ORD = 0, COMMON_ORD
                UH = UH + OLD_DEG( ORD, I ) * POLY_C( ORD, Q_ERR( QUAD, 1 ) )
            END DO

            UE = EXT_U( XQ, PI, TIME_OUT )

            SUM = SUM + ( UE - UH ) ** 2 * Q_ERR( QUAD, 2 ) * DX( I )

        END DO

    END DO

    ERR_L2 = DSQRT( SUM )

    RETURN
    !----------------------------------------- CALCULATIONS HAVE FINISHED ---------------------------------------!
  END SUBROUTINE GET_L2

  SUBROUTINE GET_LINFTY( IMAX, U_EXACT, U_NUM, ERR_LINF )

    IMPLICIT NONE

    INTEGER,  INTENT( IN )                          :: IMAX
    REAL( 8 ), DIMENSION( 1 : IMAX ), INTENT( IN )  :: U_EXACT, U_NUM

    REAL( 8 ),                        INTENT( OUT ) :: ERR_LINF

    INTEGER                                         :: I
    REAL( 8 )                                       :: DIFF
    !----------------------------------------- CALCULATIONS HAVE STARTED ----------------------------------------!

    ERR_LINF = 0.0D0

    DO I = 1, IMAX

        DIFF = DABS( U_EXACT( I ) - U_NUM( I ) )

        IF ( DIFF .GT. ERR_LINF ) THEN

            ERR_LINF = DIFF

        END IF

    END DO

    RETURN
    !----------------------------------------- CALCULATIONS HAVE FINISHED ---------------------------------------!
  END SUBROUTINE GET_LINFTY

  ! @@@@@

  ! @@@@@

  REAL( 8 ) FUNCTION INI_U( X, PI )

    IMPLICIT NONE

    REAL( 8 ), INTENT( IN ) :: X, PI
    !----------------------------------------- CALCULATIONS HAVE STARTED ----------------------------------------!

    ! SET UP THE INITIAL FUNCTION
    INI_U = DSIN( PI * X ) ** 1
    ! INI_U = DSIN( PI * X ) ** 2
    ! INI_U = 1.0D0 + 0.5D0 * DSIN( 4.0D0 * PI * X ) ** 1
    ! INI_U = 1.0D0 + 0.5D0 * DSIN( 4.0D0 * PI * X ) ** 2
    ! INI_U = ( DSIN( PI * X ) + DCOS( PI * X ) ) ** 1
    ! INI_U = ( DSIN( PI * X ) + DCOS( PI * X ) ) ** 2

    RETURN
    !----------------------------------------- CALCULATIONS HAVE FINISHED ---------------------------------------!
  END FUNCTION INI_U

  REAL( 8 ) FUNCTION EXT_U( X, PI, TIME_OUT )

    IMPLICIT NONE

    REAL( 8 ), INTENT( IN ) :: X, PI
    REAL( 8 ), INTENT( IN ) :: TIME_OUT
    !----------------------------------------- CALCULATIONS HAVE STARTED ----------------------------------------!

    ! SET UP THE EXACT SOLUTION
    EXT_U = DSIN( PI * ( X - TIME_OUT ) ) ** 1
    ! EXT_U = DSIN( PI * ( X - TIME_OUT ) ) ** 2
    ! EXT_U = 1.0D0 + 0.5D0 * DSIN( 4.0D0 * PI * ( X - TIME_OUT ) ) ** 1
    ! EXT_U = 1.0D0 + 0.5D0 * DSIN( 4.0D0 * PI * ( X - TIME_OUT ) ) ** 2
    ! EXT_U = ( DSIN( PI * ( X - TIME_OUT ) ) + DCOS( PI * ( X - TIME_OUT ) ) ) ** 1
    ! EXT_U = ( DSIN( PI * ( X - TIME_OUT ) ) + DCOS( PI * ( X - TIME_OUT ) ) ) ** 2

    RETURN
    !----------------------------------------- CALCULATIONS HAVE FINISHED ---------------------------------------!
  END FUNCTION EXT_U

  ! @@@@@

  ! @@@@@

  SUBROUTINE GET_TIMESTEP( I_BC, IMAX, DX, FLAG_1, TIME_OUT, CFL, DT )

    IMPLICIT NONE

    INTEGER,   INTENT( IN )  :: I_BC, IMAX
    REAL( 8 ), INTENT( IN )  :: DX( - I_BC + 1 : IMAX + I_BC )
    REAL( 8 ), INTENT( IN )  :: FLAG_1, TIME_OUT, CFL

    REAL( 8 ), INTENT( OUT ) :: DT
    !----------------------------------------- CALCULATIONS HAVE STARTED ----------------------------------------!

    ! DETERMINE THE TIME STEP SIZE
    DT = CFL * MAXVAL( DX ) 
    
    ! RESTRICTION
    IF ( ( FLAG_1 + DT ) .GE. TIME_OUT ) THEN

        DT = TIME_OUT - FLAG_1

    END IF

    RETURN
    !----------------------------------------- CALCULATIONS HAVE FINISHED ---------------------------------------!
  END SUBROUTINE GET_TIMESTEP
  
  SUBROUTINE GET_EXACT_LAMBDA( I_BC, IMAX, DX, LAMBDA )

    IMPLICIT NONE   

    INTEGER,   INTENT( IN )  :: I_BC, IMAX
    REAL( 8 ), INTENT( IN )  :: DX    ( - I_BC + 1 : IMAX + I_BC )
    REAL( 8 ), INTENT( OUT ) :: LAMBDA( - I_BC + 1 : IMAX + I_BC )

    ! COUNTERS
    INTEGER                  :: I
    !----------------------------------------- CALCULATIONS HAVE STARTED ----------------------------------------!

    DO I = 1, IMAX

        LAMBDA( I ) = - ( PI * DX( I ) ) ** 2 / 4.0D0

    END DO

    DO I = 1, I_BC

        LAMBDA( - I + 1 )  = LAMBDA( IMAX - I + 1 )
        LAMBDA( IMAX + I ) = LAMBDA( I )

    END DO

    RETURN
    !----------------------------------------- CALCULATIONS HAVE FINISHED ---------------------------------------!
  END SUBROUTINE GET_EXACT_LAMBDA

  ! @@@@@

  ! @@@@@
  
  SUBROUTINE GET_LAMBDA( MAX_ORD, I_BC, IMAX, X, DX, DEG, LAMBDA )

    IMPLICIT NONE

    INTEGER,   INTENT( IN )    :: MAX_ORD
    INTEGER,   INTENT( IN )    :: I_BC, IMAX
    REAL( 8 ), INTENT( IN )    :: X     ( - I_BC + 1 : IMAX + I_BC )
    REAL( 8 ), INTENT( IN )    :: DX    ( - I_BC + 1 : IMAX + I_BC )
    REAL( 8 ), INTENT( INOUT ) :: DEG   ( 0 : MAX_ORD, - I_BC + 1 : IMAX + I_BC )
    REAL( 8 ), INTENT( OUT )   :: LAMBDA( - I_BC + 1 : IMAX + I_BC )

    ! COUNTERS
    INTEGER                    :: ORD
    INTEGER                    :: I
    INTEGER                    :: LEFT_I, RIGHT_I

    ! CALCULATION VARIABLES
    REAL( 8 )                  :: U_1
    REAL( 8 )                  :: U_3
    !----------------------------------------- CALCULATIONS HAVE STARTED ----------------------------------------!

    CALL GET_BOUNDARY( MAX_ORD, I_BC, IMAX, DEG )

    ! COMPUTE CELL-WISE SHAPE PARAMETERS
    DO I = 1, IMAX

        U_1 = 0.0D0
        U_3 = 0.0D0

        DO ORD = 0, MAX_ORD

            U_1 = U_1 + DEG( ORD, I ) * D_POLY_C ( ORD, 0.0D0 ) 
            U_3 = U_3 + DEG( ORD, I ) * D3_POLY_C( ORD, 0.0D0 ) 

        END DO

        LAMBDA( I ) = U_3 / ( 4.0D0 * U_1 )

        ! IF( DABS( U_1 ) .LE. 1D-12 ) THEN
        
        !    WRITE( *, * ) I, REAL( U_1 ), REAL( DX( I ) ** 1 )
        !    WRITE( *, * ) I, REAL( U_3 ), REAL( DX( I ) ** 3 )
           
        !    STOP

        ! END IF

        write(*, *) i,"U_1 = ", real(U_1),real(MAXVAL(DX))
        write(*, *) i,"U_3 = ",  real(U_3),real(MAXVAL(DX)**3)
        ! write(*, *) i,real(LAMBDA(I)),real(MAXVAL(DX)**2/4 )

    END DO

    stop

    DO I = 1, I_BC

        LAMBDA( - I + 1 )  = LAMBDA( IMAX - I + 1 )
        LAMBDA( IMAX + I ) = LAMBDA( I )

    END DO

    ! LAMBDA( : ) = 0.0D0

    RETURN
    !----------------------------------------- CALCULATIONS HAVE FINISHED ---------------------------------------!
  END SUBROUTINE GET_LAMBDA

  SUBROUTINE GET_RESIDUAL( POLY_ORD, Q, M_STAGE, L_B, R_B, I_BC, IMAX, DX, LAMBDA, DEG, F_FLUX )

    IMPLICIT NONE

    INTEGER,   INTENT( IN )    :: POLY_ORD
    INTEGER,   INTENT( IN )    :: I_BC, IMAX
    REAL( 8 ), INTENT( IN )    :: Q      ( 1 : NQ, 1 : 2 )
    REAL( 8 ), INTENT( IN )    :: M_STAGE( 0 : POLY_ORD, 0 : POLY_ORD, - I_BC + 1 : IMAX_MAX + I_BC )
    REAL( 8 ), INTENT( IN )    :: L_B    ( 0 : POLY_ORD, - I_BC + 1 : IMAX_MAX + I_BC )
    REAL( 8 ), INTENT( IN )    :: R_B    ( 0 : POLY_ORD, - I_BC + 1 : IMAX_MAX + I_BC )
    REAL( 8 ), INTENT( IN )    :: DX     ( - I_BC + 1 : IMAX + I_BC )
    REAL( 8 ), INTENT( IN )    :: LAMBDA ( - I_BC + 1 : IMAX + I_BC )

    REAL( 8 ), INTENT( INOUT ) :: DEG    ( 0 : POLY_ORD, - I_BC + 1 : IMAX + I_BC )
    REAL( 8 ), INTENT( OUT )   :: F_FLUX ( 0 : POLY_ORD, 1 : IMAX )

    ! COUNTERS
    INTEGER                                                        :: ORD
    INTEGER                                                        :: I
    INTEGER                                                        :: QUAD

    ! VOLUME INTEGRATION VARIABLES
    REAL( 8 ), DIMENSION( 1 : NQ )                                 :: QUAD_SUM
    REAL( 8 ), DIMENSION( 1 : NQ, 1 : IMAX )                       :: QUAD_U

    ! VOLUME INTEGRATION VARIABLE
    REAL( 8 ), DIMENSION( 0 : POLY_ORD, 1 : IMAX )                 :: F_FLUX_V

    ! CELL INTERFACE VALUES
    REAL( 8 ), DIMENSION( 0 : 1 )                                  :: TEMP
    REAL( 8 ), DIMENSION( 0 : IMAX + 1 )                           :: U_M
    REAL( 8 ), DIMENSION( - 1 : IMAX )                             :: U_P

    ! FLUX SPLITTING VARIABLE
    REAL( 8 ), DIMENSION( 0 : IMAX )                               :: HAT_F

    ! BOUNDARY INTEGRATION VARIABLE
    REAL( 8 ), DIMENSION( 0 : POLY_ORD, 1 : IMAX )                 :: F_FLUX_B

    ! SINGULARITY DETECTOR
    REAL( 8 ), DIMENSION( 1 : IMAX )                               :: THETA
    REAL( 8 ), DIMENSION( 0 : POLY_ORD, - I_BC + 1 : IMAX + I_BC ) :: DEG_REC
    !----------------------------------------- CALCULATIONS HAVE STARTED ----------------------------------------!

    DEG_REC( :, : ) = DEG( :, : )

    ! APPLY LIMITING ALGORITHM CELL BY CELL
    IF (WENO_SW .EQ. 1) THEN

        CALL GET_DETECTOR( POLY_ORD, I_BC, IMAX, DEG, THETA )

        DO I = 1, IMAX

            IF ( THETA( I ) .GT. 1.0D0 ) THEN

                CALL GET_WENO( POLY_ORD, Q, M_STAGE, I_BC, IMAX, I, DEG, DEG_REC )

            END IF

        END DO
    
    END IF 

    CALL GET_BOUNDARY( POLY_ORD, I_BC, IMAX, DEG_REC )

    !============================================================================================================!
    ! VOLUME INTEGRAL PART
    !============================================================================================================!

    ! COMPUTE THE VOLUME INTEGRAL BY USING QUADRATURES
    DO I = 1, IMAX

        QUAD_SUM( : ) = 0.0D0

        DO ORD = 0, POLY_ORD

            DO QUAD = 1, NQ

                QUAD_SUM( QUAD ) = QUAD_SUM( QUAD ) + DEG_REC( ORD, I ) * POLY( ORD, Q( QUAD, 1 ), LAMBDA( I ) )

            END DO

        END DO

        DO QUAD = 1, NQ

            QUAD_U( QUAD, I ) = QUAD_SUM( QUAD )

        END DO

    END DO

    DO I = 1, IMAX

        DO ORD = 0, POLY_ORD

            F_FLUX_V( ORD, I ) = 0.0D0

            DO QUAD = 1, NQ

                F_FLUX_V( ORD, I ) = F_FLUX_V( ORD, I ) &
                                   + FLUX( QUAD_U( QUAD, I ) ) * D_POLY( ORD, Q( QUAD, 1 ), LAMBDA( I ) ) * Q( QUAD, 2 )

            END DO

        END DO

    END DO

    !============================================================================================================!
    ! BOUNDARY INTEGRAL PART
    !============================================================================================================!

    ! FOR GIVEN DEGREES OF FREEDOM, COMPUTE CONVEX SUMMATIONS AT CELL INTERFACES
    DO I = 0, IMAX + 1

        TEMP( : ) = 0.0D0

        DO ORD = 0, POLY_ORD

            ! LEFT BOUNDARY
            TEMP( 0 ) = TEMP( 0 ) + DEG_REC( ORD, I ) * L_B( ORD, I )

            ! RIGHT BOUNDARY
            TEMP( 1 ) = TEMP( 1 ) + DEG_REC( ORD, I ) * R_B( ORD, I )

        END DO

        U_P( I - 1 ) = TEMP( 0 )
        U_M( I + 0 ) = TEMP( 1 )

    END DO

    ! TAKE THE APPROPRIATE NUMERICAL TRACE
    ! AS THE AVERAGED MONOTONE FLUX, APPLY LAX-FRIEDRICHS (LF) FLUX
    DO I = 0, IMAX

        HAT_F( I ) = FLUX( U_M( I ) ) + FLUX( U_P( I ) ) + 1.0D0 * U_M( I ) - 1.0D0 * U_P( I )
        HAT_F( I ) = HAT_F( I ) / 2.0D0

    END DO

    ! COMPUTE THE BOUNDARY INTEGRAL BY USING THE FUNDAMENTAL THEOREM OF CALCULUS
    DO I = 1, IMAX

        DO ORD = 0, POLY_ORD

            F_FLUX_B( ORD, I ) = HAT_F( I - 1 ) * L_B( ORD, I ) - HAT_F( I + 0 ) * R_B( ORD, I )

        END DO

    END DO

    !============================================================================================================!
    ! GATHERING ALL INTEGRALS
    !============================================================================================================!

    ! COMPUTE LOCAL RESIDUALS BY USING CALCULATED INTEGRAL VALUES
    DO I = 1, IMAX

        DO ORD = 0, POLY_ORD

            F_FLUX( ORD, I ) = M_STAGE( ORD, ORD, I ) * ( F_FLUX_B( ORD, I ) + F_FLUX_V( ORD, I ) ) / DX( I )

        END DO

    END DO

    RETURN
    !----------------------------------------- CALCULATIONS HAVE FINISHED ---------------------------------------!
  END SUBROUTINE GET_RESIDUAL

  SUBROUTINE GET_WENO( POLY_ORD, Q, M_STAGE, I_BC, IMAX, TARGET_I, DEG_0, DEG )

    IMPLICIT NONE
  
    INTEGER,   INTENT( IN )                   :: POLY_ORD
    INTEGER,   INTENT( IN )                   :: I_BC, IMAX
    INTEGER,   INTENT( IN )                   :: TARGET_I
    REAL( 8 ), INTENT( IN )                   :: Q( 1 : NQ, 1 : 2 )
    REAL( 8 ), INTENT( IN )                   :: M_STAGE( 0 : POLY_ORD, 0 : POLY_ORD, - I_BC + 1 : IMAX + I_BC )
    REAL( 8 ), INTENT( IN )                   :: DEG_0( 0 : POLY_ORD, - I_BC + 1 : IMAX + I_BC )

    REAL( 8 ), INTENT( INOUT )                :: DEG( 0 : POLY_ORD, - I_BC + 1 : IMAX + I_BC )

    ! COUNTERS
    INTEGER                                   :: I
    INTEGER                                   :: J
    INTEGER                                   :: ORD
    INTEGER                                   :: QUAD
    INTEGER                                   :: L

    ! LINEAR COEFFICIENTS FOR THE WENO TYPE RECONSTRUCTION
    REAL( 8 ), DIMENSION( - 1 : 1 )           :: LINEAR_COEFF

    ! CANDIDATE POLYNOMIALS
    REAL( 8 ), DIMENSION( - 1 : 1, 1 : NQ )   :: U_CAND
    REAL( 8 ), DIMENSION( - 1 : 1, 1 : NQ )   :: D1_CAND
    REAL( 8 ), DIMENSION( - 1 : 1, 1 : NQ )   :: D2_CAND

    ! CELL-AVERAGE APPROXIMATIONS (ON TARGET CELL)
    REAL( 8 ), DIMENSION( - 1 : 1 )           :: AVG

    ! FINAL RECONSTRUCTED POLYNOMIAL VALUES AT QUADRATURE NODES
    REAL( 8 ), DIMENSION( 1 : NQ )            :: NEW_U

    ! PARAMETERS FOR SMOOTHNESS INDICATORS
    REAL( 8 )                                 :: EPS
    REAL( 8 ), DIMENSION( - 1 : 1 )           :: BETA
    REAL( 8 )                                 :: TAU
    REAL( 8 ), DIMENSION( - 1 : 1 )           :: ALPHA
    REAL( 8 )                                 :: ALPHA_SUM
    REAL( 8 ), DIMENSION( - 1 : 1 )           :: OMEGA

    ! TEMPORARY SCALARS
    REAL( 8 )                                 :: XI_J
    REAL( 8 )                                 :: XQ
    REAL( 8 )                                 :: VAL
    REAL( 8 )                                 :: D1VAL
    REAL( 8 )                                 :: D2VAL
    !----------------------------------------- CALCULATIONS HAVE STARTED ----------------------------------------!
    
    EPS = 1D-40

    ! SET UP LINEAR COEFFICIENTS
    LINEAR_COEFF( - 1 ) = 1.0D-3
    LINEAR_COEFF( + 1 ) = LINEAR_COEFF( - 1 )
    LINEAR_COEFF( + 0 ) = 1.0D0 - LINEAR_COEFF( - 1 ) - LINEAR_COEFF( + 1 )

    I = TARGET_I

    ! BUILD CANDIDATE POLYNOMIALS ON TARGET CELL I
    DO L = - 1, 1

        J = I + L

        DO QUAD = 1, NQ
                
            ! PHYSICAL LOCATION OF QUADRATURE POINT
            XQ = X( I ) + Q( QUAD, 1 ) * DX( I )

            ! MAP THIS PHYSICAL POINT INTO CELL J'S LOCAL COORDINATE
            XI_J = ( XQ - X( J ) ) / DX( J )

            ! EVALUATE DG POLYNOMIAL AND ITS DERIVATIVES OF CELL J AT XQ
            VAL   = 0.0D0
            D1VAL = 0.0D0
            D2VAL = 0.0D0

            DO ORD = 0, POLY_ORD

                VAL   = VAL   + DEG_0( ORD, J ) * POLY   ( ORD, XI_J, LAMBDA( J ) )
                D1VAL = D1VAL + DEG_0( ORD, J ) * D_POLY ( ORD, XI_J, LAMBDA( J ) )
                D2VAL = D2VAL + DEG_0( ORD, J ) * D2_POLY( ORD, XI_J, LAMBDA( J ) )

            END DO

            U_CAND ( L, QUAD ) = VAL
            D1_CAND( L, QUAD ) = D1VAL
            D2_CAND( L, QUAD ) = D2VAL

        END DO

        ! APPROXIMATE CELL AVERAGE ON TARGET CELL I FOR CANDIDATE L
        AVG( L ) = 0.0D0

        DO QUAD = 1, NQ

            AVG( L ) = AVG( L ) + U_CAND( L, QUAD ) * Q( QUAD, 2 )

        END DO

    END DO

    ! MODIFY CANDIDATES TO MATCH TARGET CELL AVERAGE
    DO QUAD = 1, NQ

        U_CAND( - 1, QUAD ) = U_CAND( - 1, QUAD ) - AVG( - 1 ) + AVG( 0 )
        U_CAND( + 0, QUAD ) = U_CAND( + 0, QUAD )
        U_CAND( + 1, QUAD ) = U_CAND( + 1, QUAD ) - AVG( + 1 ) + AVG( 0 )

    END DO

    ! COMPUTE LOCAL SMOOTHNESS INDICATORS
    DO L = - 1, 1

        BETA( L ) = 0.0D0

        DO QUAD = 1, NQ

            BETA( L ) = BETA( L ) &
                      + ( D1_CAND( L, QUAD ) ** 2 ) * Q( QUAD, 2 ) &
                      + ( D2_CAND( L, QUAD ) ** 2 ) * Q( QUAD, 2 )

        END DO

    END DO

    TAU = DABS( BETA( + 1 ) - BETA( - 1 ) )

    ! COMPUTE NONLINEAR WEIGHTS
    ALPHA( : ) = 0.0D0

    DO L = - 1, 1
        
        ! WENO - JS
        ! ALPHA( L ) = LINEAR_COEFF( L ) / ( ( BETA( L ) + EPS ) ** 2 )

        ! WENO - Z
        ALPHA( L ) = LINEAR_COEFF( L ) * ( 1.0D0 + ( TAU /  ( BETA( L ) + EPS ) ) ** 2 )

    END DO

    ALPHA_SUM = ALPHA( - 1 ) + ALPHA( 0 ) + ALPHA( + 1 )

    OMEGA( : ) = 0.0D0

    DO L = - 1, 1

        OMEGA( L ) = ALPHA( L ) / ALPHA_SUM
            
    END DO

    ! BUILD RECONSTRUCTED POLYNOMIAL VALUES NEW_U
    DO QUAD = 1, NQ

        NEW_U( QUAD ) = 0.0D0

        DO L = - 1, 1

            NEW_U( QUAD ) = NEW_U( QUAD ) + OMEGA( L ) * U_CAND( L, QUAD )

        END DO

    END DO

    ! EVALUATE DEGREES OF FREEDOM (THE MOMENTS) BY USING THE MODIFIED DG POLYNOMIALS
    DO ORD = 1, POLY_ORD

        VAL = 0.0D0

        DO QUAD = 1, NQ

            VAL = VAL + NEW_U( QUAD ) * POLY( ORD, Q( QUAD, 1 ), LAMBDA( I ) ) * Q( QUAD, 2 )

        END DO

        DEG( ORD, I ) = M_STAGE( ORD, ORD, I ) * VAL

    END DO
  
    RETURN
    !----------------------------------------- CALCULATIONS HAVE FINISHED ---------------------------------------!
  END SUBROUTINE GET_WENO

  SUBROUTINE GET_DETECTOR( POLY_ORD, I_BC, IMAX, DEG, THETA )

    IMPLICIT NONE

    INTEGER,   INTENT( IN )                            :: POLY_ORD
    INTEGER,   INTENT( IN )                            :: I_BC
    INTEGER,   INTENT( IN )                            :: IMAX
    REAL( 8 ), INTENT( INOUT )                         :: DEG  ( 0 : POLY_ORD, - I_BC + 1 : IMAX + I_BC )
    
    REAL( 8 ), INTENT( OUT )                           :: THETA( 1 : IMAX )

    ! COUNTERS
    INTEGER                                            :: I
    INTEGER                                            :: ORD

    ! PARAMETERS FOR INDICATOR
    INTEGER,  PARAMETER                                :: P   = 2
    REAL( 8 ), PARAMETER                               :: XI  = 1.0D0

    ! FIRST AND SECOND DERIVATIVES AT CELL CENTERS
    REAL( 8 ), DIMENSION( - I_BC + 1 : IMAX + I_BC )  :: C_DX1_U
    REAL( 8 ), DIMENSION( - I_BC + 1 : IMAX + I_BC )  :: C_DX2_U

    ! SMOOTHNESS INDICATORS 
    REAL( 8 ), DIMENSION( - 2 : 2 )                   :: IS
    !----------------------------------------- CALCULATIONS HAVE STARTED ----------------------------------------!
    
    CALL GET_BOUNDARY( POLY_ORD, I_BC, IMAX, DEG )
    
    ! RECONSTRUCT FIRST AND SECOND DERIVATIVES AT CELL CENTER ( XI = 0.0D0 )
    DO I = - I_BC + 1, IMAX + I_BC

        C_DX1_U( I ) = 0.0D0
        C_DX2_U( I ) = 0.0D0

        DO ORD = 0, POLY_ORD

            C_DX1_U( I ) = C_DX1_U( I ) &
                         + DEG( ORD, I ) * D_POLY ( ORD, 0.0D0, LAMBDA( I ) )

            C_DX2_U( I ) = C_DX2_U( I ) &
                         + DEG( ORD, I ) * D2_POLY( ORD, 0.0D0, LAMBDA( I ) )

        END DO

    END DO

    DO I = 1, IMAX

        ! LOCAL DETECTORS
        IS( - 2 ) = DABS( ( ( C_DX1_U( I - 1 ) / DX( I - 1 ) ** 1 ) / 1.0D0 ) * DX( I - 1 ) ** 0 ) ** P &
                  + XI &
                  * DABS( ( ( C_DX2_U( I - 1 ) / DX( I - 1 ) ** 2 ) / 2.0D0 ) * DX( I - 1 ) ** 0 ) ** P

        IS( + 2 ) = DABS( ( ( C_DX1_U( I + 1 ) / DX( I + 1 ) ** 1 ) / 1.0D0 ) * DX( I + 1 ) ** 0 ) ** P &
                  + XI &
                  * DABS( ( ( C_DX2_U( I + 1 ) / DX( I + 1 ) ** 2 ) / 2.0D0 ) * DX( I + 1 ) ** 0 ) ** P

        ! GLOBAL DETECTORS
        IS( - 1 ) = DABS( ( ( C_DX1_U( I - 1 ) + C_DX1_U( I + 0 ) ) / DX( I - 1 ) ** 1 ) / 2.0D0 * DX( I - 1 ) ** 0 ) ** P &
                  + XI &
                  * DABS( ( ( C_DX2_U( I - 1 ) + C_DX2_U( I + 0 ) ) / DX( I - 1 ) ** 2 ) / 4.0D0 * DX( I - 1 ) ** 0 ) ** P
      
        IS( + 1 ) = DABS( ( ( C_DX1_U( I + 0 ) + C_DX1_U( I + 1 ) ) / DX( I + 1 ) ** 1 ) / 2.0D0 * DX( I + 1 ) ** 0 ) ** P &
                  + XI &
                  * DABS( ( ( C_DX2_U( I + 0 ) + C_DX2_U( I + 1 ) ) / DX( I + 1 ) ** 2 ) / 4.0D0 * DX( I + 1 ) ** 0 ) ** P

        ! SET UP THE SINGULARITY DETECTOR SCHEME
        THETA( I ) = ( 1.0D0 / 2.0D0 ) * ( MAX( IS( - 1 ), IS( + 1 ) ) / ( MIN( IS( - 2 ), IS( + 2 ) ) + 1.0D-40 ) )
    
    END DO

    RETURN
    !----------------------------------------- CALCULATIONS HAVE FINISHED ---------------------------------------!
  END SUBROUTINE GET_DETECTOR

  ! @@@@@

END PROGRAM MODAL_DG_1D