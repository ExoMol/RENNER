      SUBROUTINE TR15SS ( XR )
      IMPLICIT REAL*8 (A-H,O-Z)
C            :
C DATE       : 15.03.1994
C AUTHOR     : PER JENSEN
C UPDATES    :
C LANGUAGE   : FORTRAN
C PURPOSE    : SPECIAL ADAPTION OF THE TRRFEQ ROUTINE. MADE NECESSARY
C            : BY THE INTRODUCTION OF THE GAMMA ANGLE IN RENNER
      REAL*8 XR(9),XDUMM(15)
C
      XDUMM( 1) = XR( 1)
      XDUMM( 2) = XR( 2)
      XDUMM( 3) = XR( 2)
      XDUMM( 4) = XR( 3)
      XDUMM( 5) = XR( 3)
      XDUMM( 6) = XR( 4)
      XDUMM( 7) = XR( 5)
      XDUMM( 8) = XR( 5)
      XDUMM( 9) = XR( 6)
      XDUMM(10) = XR( 6)
      XDUMM(11) = XR( 7)
      XDUMM(12) = XR( 7)
      XDUMM(13) = XR( 8)
      XDUMM(14) = XR( 8)
      XDUMM(15) = XR( 9)
C
      CALL TRRFEQ ( XDUMM )
C
      XR( 1) = XDUMM( 1)
      XR( 2) = XDUMM( 2)
      XR( 3) = XDUMM( 4)
      XR( 4) = XDUMM( 6)
      XR( 5) = XDUMM( 7)
      XR( 6) = XDUMM( 9)
      XR( 7) = XDUMM(11)
      XR( 8) = XDUMM(13)
      XR( 9) = XDUMM(15)
      RETURN
      END
C
C
      SUBROUTINE TR15SA ( XR )
      IMPLICIT REAL*8 (A-H,O-Z)
C            :
C DATE       : 15.03.1994
C AUTHOR     : PER JENSEN
C UPDATES    :
C LANGUAGE   : FORTRAN
C PURPOSE    : SPECIAL ADAPTION OF THE TRRFEQ ROUTINE. MADE NECESSARY
C            : BY THE INTRODUCTION OF THE GAMMA ANGLE IN RENNER
      REAL*8 XR(8),XDUMM(15)
C
      XDUMM( 1) =  0.0D+00
      XDUMM( 2) =  XR( 1)
      XDUMM( 3) = -XR( 1)
      XDUMM( 4) =  XR( 2)
      XDUMM( 5) = -XR( 2)
      XDUMM( 6) =  XR( 3)
      XDUMM( 7) =  XR( 4)
      XDUMM( 8) = -XR( 4)
      XDUMM( 9) =  XR( 5)
      XDUMM(10) = -XR( 5)
      XDUMM(11) =  XR( 6)
      XDUMM(12) = -XR( 6)
      XDUMM(13) =  XR( 7)
      XDUMM(14) = -XR( 7)
      XDUMM(15) =  XR( 8)
C
      CALL TRRFEQ ( XDUMM )
C
      XR( 1)  =  XDUMM( 2)
      XR( 2)  =  XDUMM( 4)
      XR( 3)  =  XDUMM( 6)
      XR( 4)  =  XDUMM( 7)
      XR( 5)  =  XDUMM( 9)
      XR( 6)  =  XDUMM(11)
      XR( 7)  =  XDUMM(13)
      XR( 8)  =  XDUMM(15)
      RETURN
      END
C
C
      SUBROUTINE TR15A  ( XR1, XR2 )
      IMPLICIT REAL*8 (A-H,O-Z)
C            :
C DATE       : 15.03.1994
C AUTHOR     : PER JENSEN
C UPDATES    :
C LANGUAGE   : FORTRAN
C PURPOSE    : SPECIAL ADAPTION OF THE TRRFEQ ROUTINE. MADE NECESSARY
C            : BY THE INTRODUCTION OF THE GAMMA ANGLE IN RENNER
      REAL*8 XR1(9), XR2(6) , XDUMM(15)
C
      XDUMM( 1)  = XR1( 1)
      XDUMM( 2)  = XR1( 2)
      XDUMM( 3)  = XR2( 1)
      XDUMM( 4)  = XR1( 3)
      XDUMM( 5)  = XR2( 2)
      XDUMM( 6)  = XR1( 4)
      XDUMM( 7)  = XR1( 5)
      XDUMM( 8)  = XR2( 3)
      XDUMM( 9)  = XR1( 6)
      XDUMM(10)  = XR2( 4)
      XDUMM(11)  = XR1( 7)
      XDUMM(12)  = XR2( 5)
      XDUMM(13)  = XR1( 8)
      XDUMM(14)  = XR2( 6)
      XDUMM(15)  = XR1( 9)
C
      CALL TRRFEQ ( XDUMM )
C
      XR1( 1) =   XDUMM( 1)
      XR1( 2) =   XDUMM( 2)
      XR2( 1) =   XDUMM( 3)
      XR1( 3) =   XDUMM( 4)
      XR2( 2) =   XDUMM( 5)
      XR1( 4) =   XDUMM( 6)
      XR1( 5) =   XDUMM( 7)
      XR2( 3) =   XDUMM( 8)
      XR1( 6) =   XDUMM( 9)
      XR2( 4) =   XDUMM(10)
      XR1( 7) =   XDUMM(11)
      XR2( 5) =   XDUMM(12)
      XR1( 8) =   XDUMM(13)
      XR2( 6) =   XDUMM(14)
      XR1( 9) =   XDUMM(15)
      RETURN
      END
C
C
      SUBROUTINE TR10   ( XR )
      IMPLICIT REAL*8 (A-H,O-Z)
C            :
C DATE       : 15.03.1994
C AUTHOR     : PER JENSEN
C UPDATES    :
C LANGUAGE   : FORTRAN
C PURPOSE    : SPECIAL ADAPTION OF THE TRRFEQ ROUTINE. MADE NECESSARY
C            : BY THE INTRODUCTION OF THE GAMMA ANGLE IN RENNER
      REAL*8 XR(9),XDUMM(15)
C
      XDUMM( 1)   = 0.0D+00
      XDUMM( 2)   = XR( 1)
      XDUMM( 3)   = XR( 2)
      XDUMM( 4)   = XR( 3)
      XDUMM( 5)   = XR( 4)
      XDUMM( 6)   = XR( 5)
      XDUMM( 7)   = XR( 6)
      XDUMM( 8)   = XR( 7)
      XDUMM( 9)   = XR( 8)
      XDUMM(10)   = XR( 9)
      XDUMM(11)   = 0.0D+00
      XDUMM(12)   = 0.0D+00
      XDUMM(13)   = 0.0D+00
      XDUMM(14)   = 0.0D+00
      XDUMM(15)   = 0.0D+00
C
      CALL TRRFEQ ( XDUMM )
C
      XR( 1)  =   XDUMM( 2) 
      XR( 2)  =   XDUMM( 3) 
      XR( 3)  =   XDUMM( 4) 
      XR( 4)  =   XDUMM( 5) 
      XR( 5)  =   XDUMM( 6) 
      XR( 6)  =   XDUMM( 7) 
      XR( 7)  =   XDUMM( 8) 
      XR( 8)  =   XDUMM( 9) 
      XR( 9)  =   XDUMM(10) 
      RETURN
      END
C
C
      SUBROUTINE TR06SS ( XR )
      IMPLICIT REAL*8 (A-H,O-Z)
C            :
C DATE       : 15.03.1994
C AUTHOR     : PER JENSEN
C UPDATES    :
C LANGUAGE   : FORTRAN
C PURPOSE    : SPECIAL ADAPTION OF THE TRRFEQ ROUTINE. MADE NECESSARY
C            : BY THE INTRODUCTION OF THE GAMMA ANGLE IN RENNER
      REAL*8 XR(4),XDUMM(15)
C
      XDUMM( 1)  = XR( 1)
      XDUMM( 2)  = XR( 2)
      XDUMM( 3)  = XR( 2)
      XDUMM( 4)  = XR( 3)
      XDUMM( 5)  = XR( 3)
      XDUMM( 6)  = XR( 4)
      XDUMM( 7)  = 0.0D+00
      XDUMM( 8)  = 0.0D+00
      XDUMM( 9)  = 0.0D+00
      XDUMM(10)  = 0.0D+00
      XDUMM(11)  = 0.0D+00
      XDUMM(12)  = 0.0D+00
      XDUMM(13)  = 0.0D+00
      XDUMM(14)  = 0.0D+00
      XDUMM(15)  = 0.0D+00
C
      CALL TRRFEQ ( XDUMM )
C
      XR( 1) =   XDUMM( 1)  
      XR( 2) =   XDUMM( 3)  
      XR( 3) =   XDUMM( 5)  
      XR( 4) =   XDUMM( 6)  
      RETURN
      END
C
C
      SUBROUTINE TR06A  ( XR1, XR2 )
      IMPLICIT REAL*8 (A-H,O-Z)
C            :
C DATE       : 15.03.1994
C AUTHOR     : PER JENSEN
C UPDATES    :
C LANGUAGE   : FORTRAN
C PURPOSE    : SPECIAL ADAPTION OF THE TRRFEQ ROUTINE. MADE NECESSARY
C            : BY THE INTRODUCTION OF THE GAMMA ANGLE IN RENNER
      REAL*8 XR1(4), XR2(2), XDUMM(15)
C
      XDUMM( 1) = XR1( 1)
      XDUMM( 2) = XR1( 2)
      XDUMM( 3) = XR2( 1)
      XDUMM( 4) = XR1( 3)
      XDUMM( 5) = XR2( 2)
      XDUMM( 6) = XR1( 4)
      XDUMM( 7) = 0.0D+00
      XDUMM( 8) = 0.0D+00
      XDUMM( 9) = 0.0D+00
      XDUMM(10) = 0.0D+00
      XDUMM(11) = 0.0D+00
      XDUMM(12) = 0.0D+00
      XDUMM(13) = 0.0D+00
      XDUMM(14) = 0.0D+00
      XDUMM(15) = 0.0D+00
C
      CALL TRRFEQ ( XDUMM )
C
      XR1( 1) =   XDUMM( 1)
      XR1( 2) =   XDUMM( 2)
      XR2( 1) =   XDUMM( 3)
      XR1( 3) =   XDUMM( 4)
      XR2( 2) =   XDUMM( 5)
      XR1( 4) =   XDUMM( 6)
      RETURN
      END
C
C
      SUBROUTINE TR15A2  ( X0, XR1, XR2 )
      IMPLICIT REAL*8 (A-H,O-Z)
C            :
C DATE       : 15.03.1994
C AUTHOR     : PER JENSEN
C UPDATES    :
C LANGUAGE   : FORTRAN
C PURPOSE    : SPECIAL ADAPTION OF THE TRRFEQ ROUTINE. MADE NECESSARY
C            : BY THE INTRODUCTION OF THE GAMMA ANGLE IN RENNER
      REAL*8 XR1(8), XR2(6),XDUMM(15)
C
      XDUMM( 1)  = X0
      XDUMM( 2)  = XR1( 1)
      XDUMM( 3)  = XR2( 1)
      XDUMM( 4)  = XR1( 2)
      XDUMM( 5)  = XR2( 2)
      XDUMM( 6)  = XR1( 3)
      XDUMM( 7)  = XR1( 4)
      XDUMM( 8)  = XR2( 3)
      XDUMM( 9)  = XR1( 5)
      XDUMM(10)  = XR2( 4)
      XDUMM(11)  = XR1( 6)
      XDUMM(12)  = XR2( 5)
      XDUMM(13)  = XR1( 7)
      XDUMM(14)  = XR2( 6)
      XDUMM(15)  = XR1( 8)
C
      CALL TRRFEQ ( XDUMM )
C
      X0      =   XDUMM( 1)  
      XR1( 1) =   XDUMM( 2)  
      XR2( 1) =   XDUMM( 3)  
      XR1( 2) =   XDUMM( 4)  
      XR2( 2) =   XDUMM( 5)  
      XR1( 3) =   XDUMM( 6)  
      XR1( 4) =   XDUMM( 7)  
      XR2( 3) =   XDUMM( 8)  
      XR1( 5) =   XDUMM( 9)  
      XR2( 4) =   XDUMM(10)  
      XR1( 6) =   XDUMM(11)  
      XR2( 5) =   XDUMM(12)  
      XR1( 7) =   XDUMM(13)  
      XR2( 6) =   XDUMM(14)  
      XR1( 8) =   XDUMM(15)  
      RETURN
      END
