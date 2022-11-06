C**********************************
      SUBROUTINE LESS(I,J,L,NFIL53)
C**********************************
      INTEGER I,J
      INTEGER NFIL53
      LOGICAL L
C
      REAL*8 XLAB11,XLAB21,E11,DE1,LIST1,INT1
      REAL*8 XLAB12,XLAB22,E12,DE2,LIST2,INT2
C
      READ(NFIL53,REC=I)
     1    XLAB11,XLAB21,E11,DE1,LIST1,INT1  
      READ(NFIL53,REC=J)
     1    XLAB12,XLAB22,E12,DE2,LIST2,INT2
C
      IF (DE1.LT.DE2) THEN
        L=.TRUE.
      ELSE
        L=.FALSE.
      ENDIF
C
      RETURN 
      END
C********************************
      SUBROUTINE SWAP(I,J,NFIL53)
C********************************
      INTEGER I,J
      INTEGER NFIL53
C
      REAL*8 XLAB11,XLAB21,E11,DE1,LIST1,INT1
      REAL*8 XLAB12,XLAB22,E12,DE2,LIST2,INT2
C     
      READ(NFIL53,REC=I)
     1    XLAB11,XLAB21,E11,DE1,LIST1,INT1
      READ(NFIL53,REC=J)
     1    XLAB12,XLAB22,E12,DE2,LIST2,INT2
C
      WRITE(NFIL53,REC=J)
     1    XLAB11,XLAB21,E11,DE1,LIST1,INT1
      WRITE(NFIL53,REC=I)
     1    XLAB12,XLAB22,E12,DE2,LIST2,INT2         
C 
      RETURN
      END
C**********************************
      SUBROUTINE SORT(NRAUS,NFIL53)
C**********************************
C
C     Sortieren nach QuickSort
C
      INTEGER MAXSTK
      PARAMETER (MAXSTK=100)
C
      INTEGER STACK(MAXSTK,2)
      INTEGER STK
      INTEGER NRAUS
      INTEGER NFIL53
      INTEGER L,R,M
      INTEGER A,E
      LOGICAL RES
C
      A=1
      E=NRAUS
      STK=0
C
10    IF ((E-A).LT.1) GOTO 20
      L=A
      R=E
      M= (L+R)/2
C
30    CALL LESS(L,M,RES,NFIL53)
      IF (RES) THEN
        L=L+1
        GOTO 30
      ENDIF
40    CALL LESS(M,R,RES,NFIL53)
      IF (RES) THEN
        R=R-1
        GOTO 40
      ENDIF
      IF (R.EQ.M) THEN
        M=L
      ELSE
        IF (L.EQ.M) M=R
      ENDIF
      CALL SWAP(R,L,NFIL53)
      IF (L.LT.M) L=L+1
      IF (M.LT.R) R=R-1
      IF ((L.NE.M).OR.(M.NE.R)) GOTO 30
C
      STK=STK+1
      IF (STK.GT.MAXSTK) THEN
        STK=MAXSTK
        WRITE(*,1000)
      ENDIF
      STACK(STK,1)=M+1
      STACK(STK,2)=E
C     A=1
      E=M-1            
      GOTO 10
C
20    IF (STK.NE.0) THEN
        A=STACK(STK,1)
        E=STACK(STK,2)
        STK=STK-1
        GOTO 10
      ENDIF
C 
1000  FORMAT(1H0,'Stack Overflow while sorting the',
     1       ' transitions')
      RETURN
      END      
