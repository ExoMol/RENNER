c****************************************************
      SUBROUTINE DIAROT(HMAT,
     1                  EIGVEC,EIGVAL,
     2                  WSPACE,IWORK,IFAIL,
     3                  EIGMAX,
     4                  IDIMTO,IROMAX,MAXPSY,INDSYM,
     5                  JX2VAL)
C****************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      include 'rensys.h'
      INTEGER PRTINT
C
      INTEGER IDIMTO
      INTEGER IROMAX
      INTEGER MAXPSY
      INTEGER INDSYM
C
      REAL*8 HMAT(IDIMTO,IDIMTO)
      REAL*8 EIGVEC(IDIMTO,IDIMTO),EIGVAL(IDIMTO)
      REAL*8 WSPACE(8*IDIMTO)
      INTEGER IWORK(5*IDIMTO)
      INTEGER IFAIL(IDIMTO)
C
      REAL*8 EIGMAX(2*MAXPSY,IROMAX,2)
C
      CHARACTER*1 JOBZ,RANGE,UPLO
      REAL*8 VL,VU
      INTEGER IL,IU
      REAL*8 ABSTOL
      INTEGER M
      INTEGER INFO
C
C     Evtl. Ausgabe der Matrix
C
      IF (ITEST .GT. 0) THEN
            WRITE (NFIL6,9000)
            DO 210 II=1,IDIMTO
            WRITE (NFIL6,9010)
            WRITE (NFIL6,9020) (HMAT(JJ,II),JJ=1,IDIMTO)
210         CONTINUE
      ENDIF
C
C     Diagonalisierung
C
      IF (IPRINT .LE. 0) THEN
        JOBZ='N'
        RANGE='A'
        UPLO='U'
        VL=0.0
        VU=0.0
        IL=1
        IU=1
        ABSTOL=0
        M=0
        CALL DSYEVX(JOBZ,RANGE,UPLO,
     1              IDIMTO,HMAT,IDIMTO,
     2              VL,VU,IL,IU,
     3              ABSTOL,
     4              M,EIGVAL,EIGVEC,IDIMTO,
     5              WSPACE,8*IDIMTO,IWORK,
     6              IFAIL,INFO)
        IF (INFO .NE. 0) THEN
          WRITE (NFIL6,6100)
          STOP
        ENDIF
      ELSE
        JOBZ='V'
        RANGE='A'
        UPLO='U'
        VL=0.0
        VU=0.0
        IL=1
        IU=1
        ABSTOL=0
        M=0
        CALL DSYEVX(JOBZ,RANGE,UPLO,
     1              IDIMTO,HMAT,IDIMTO,
     2              VL,VU,IL,IU,
     3              ABSTOL,
     4              M,EIGVAL,EIGVEC,IDIMTO,
     5              WSPACE,8*IDIMTO,IWORK,
     6              IFAIL,INFO)
        IF (INFO .NE. 0) THEN
          WRITE (NFIL6,6100)
          STOP
        ENDIF
      ENDIF
C
C     Eigenwerte speichern, maximalen Koeffizienten ermitteln
C
      DO 10 II=1,IDIMTO
        EIGMAX(INDSYM,II,1)=EIGVAL(II)
        EIGMAX(INDSYM,II,2)=DFLOAT(IDAMAX(IDIMTO,EIGVEC(1,II),1))
10    CONTINUE
C
C
C
      RETURN
6100  FORMAT(1H0,'  RENTEL.MRO.ERR  MATRIX DIAGONALIZATION  ',
     1             'FAILED IN DIAROT'/)
9000  FORMAT('0','  SYMMETRIZED MATRIX BLOCK FROM DIAROT',//)
9010  FORMAT('0')
9020  FORMAT(1H ,5F12.3)
      END
