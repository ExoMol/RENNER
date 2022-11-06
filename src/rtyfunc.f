      SUBROUTINE YP(IOSC,IW1,IW2,OVER,H)
C     --------------------(P1*Y1+Y1*P1)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 RMK(2),AAS(2)
      INTEGER NOBAS,MBASP1,LENIW
      INTEGER IOSC
      INTEGER LR,LC,IR,IC,JR,JC
      REAL*8 H(NOBAS,NOBAS),OVER(2,MBASP1,MBASP1),
     1      RN,RJ,RP1,RP2,RP3,AL,RP,RLN,RLS
      INTEGER IW1(LENIW),IW2(LENIW)
      COMMON/MORSE/RMK,AAS
      COMMON/MODIM/NOBAS,MBASP1,LENIW
C
      DO 1 LR=1,NOBAS
      IR=IW1(LR)
      JR=IW2(LR)
      DO 1 LC=1,NOBAS
      IC=IW1(LC)
      JC=IW2(LC)
      H(LR,LC)=0.0E+00
      IF((IR.EQ.IC).OR.(JR.NE.JC))GOTO 1
      RJ=ABS(DFLOAT(IR-IC))
      RN=RLN(IR,IC)
      RP1=1.0D0-RJ+RJ*(2.0D0*RN+RJ+1.0D0)/(RMK(IOSC)*2.0D0)
      RP1=RP1*OVER(IOSC,IR+1,IC+1)*RLS(IR,IC)
      RP=RP1
      H(LR,LC)=RP
1     CONTINUE
      RETURN
      END
C
C
      SUBROUTINE Y2P(IOSC,IW1,IW2,OVER,H)
C     --------------------(Y1*Y1*P1+P1*Y1*Y1)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 RMK(2),AAS(2)
      INTEGER NOBAS,MBASP1,LENIW
      INTEGER IOSC
      INTEGER LR,LC,IR,IC,JR,JC
      REAL*8 H(NOBAS,NOBAS),OVER(2,MBASP1,MBASP1),
     1      RN,RJ,RP1,RP2,RP3,AL,RP,RLN,RLS
      INTEGER IW1(LENIW),IW2(LENIW)
      COMMON/MORSE/RMK,AAS
      COMMON/MODIM/NOBAS,MBASP1,LENIW
C
      DO 1 LR=1,NOBAS
      IR=IW1(LR)
      JR=IW2(LR)
      DO 1 LC=1,NOBAS
      IC=IW1(LC)
      JC=IW2(LC)
      H(LR,LC)=0.0E+00
      IF((IR.EQ.IC).OR.(JR.NE.JC))GOTO 1
      RJ=ABS(DFLOAT(IR-IC))
      RN=RLN(IR,IC)
      AL=(RMK(IOSC)-RN)*2.0D0-1.0D0
      RP1=1.0D0-RJ*(AL-RJ)/RMK(IOSC)
      RP=RJ**4+AL*AL*(RJ+1.0D0)*RJ
      RP3=RJ*(AL-RJ)*(2.0D0*RN+1.0D0)-AL*RJ*RJ*(2.0D0*RJ+1.0D0)
      RP1=(RP1+(RP+RP3)/(8.0D0*RMK(IOSC)*RMK(IOSC)))*
     1      OVER(IOSC,IR+1,IC+1)*RLS(IR,IC)
      RP=RP1
      H(LR,LC)=RP
1     CONTINUE
      RETURN
      END
C
C
      SUBROUTINE P1Y2(IOSC1,IOSC2,IW1,IW2,OVER,H)
C     --------------------(P1*Y2)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 RMK(2),AAS(2)
      INTEGER NOBAS,MBASP1,LENIW
      INTEGER IOSC1,IOSC2
      INTEGER LR,LC,IR,IC,JR,JC
      REAL*8 H(NOBAS,NOBAS),OVER(2,MBASP1,MBASP1),
     1      RN,RJ,RP1,RP2,RP,RLN,RLS
      INTEGER IW1(LENIW),IW2(LENIW)
      COMMON/MORSE/RMK,AAS
      COMMON/MODIM/NOBAS,MBASP1,LENIW
C
      DO 1 LR=1,NOBAS
      IR=IW1(LR)
      JR=IW2(LR)
      DO 1 LC=1,NOBAS
      IC=IW1(LC)
      JC=IW2(LC)
      H(LR,LC)=0.0E+00
      IF(IR.EQ.IC)GOTO 1
      RP1=0.5E0*OVER(IOSC1,IR+1,IC+1)*RLS(IR,IC)
      IF(JR.NE.JC)GOTO 2
      RP2=(DFLOAT(JR)+0.5E0)/RMK(IOSC2)
      GOTO 3
2     RP2=-OVER(IOSC2,JR+1,JC+1)/(2.0D0*AAS(IOSC2)*RMK(IOSC2))
3     RP=RP1*RP2
      H(LR,LC)=RP
1     CONTINUE
      RETURN
      END
C
C
      SUBROUTINE P1Y2Y2(IOSC1,IOSC2,IW1,IW2,OVER,H)
C     --------------------(P1*Y2*Y2)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 RMK(2),AAS(2)
      INTEGER NOBAS,MBASP1,LENIW
      INTEGER IOSC1,IOSC2
      INTEGER LR,LC,IR,IC,JR,JC
      REAL*8 H(NOBAS,NOBAS),OVER(2,MBASP1,MBASP1),
     1      RN,RJ,RP1,RP2,RP,RLN,RLS
      INTEGER IW1(LENIW),IW2(LENIW)
      COMMON/MORSE/RMK,AAS
      COMMON/MODIM/NOBAS,MBASP1,LENIW
C
      DO 1 LR=1,NOBAS
      IR=IW1(LR)
      JR=IW2(LR)
      DO 1 LC=1,NOBAS
      IC=IW1(LC)
      JC=IW2(LC)
      H(LR,LC)=0.0E+00
      IF(IR.EQ.IC)GOTO 1
      RP1=0.5E0*OVER(IOSC1,IR+1,IC+1)*RLS(IR,IC)
      IF(JR.NE.JC)GOTO 2
      RP2=(DFLOAT(JR)+0.5E0)/RMK(IOSC2)
      GOTO 3
2     RN=RLN(JR,JC)
      RJ=ABS(DFLOAT(JR-JC))
      RP2=(RJ-1.E0)/(2.0D0*AAS(IOSC2)*RMK(IOSC2))
      RP2=RP2-RJ*(2.0D0*RN+RJ+1.0D0)/(4.0D0*RMK(IOSC2)*RMK(IOSC2))
      RP2=RP2*OVER(IOSC2,JR+1,JC+1)
3     RP=RP1*RP2
      H(LR,LC)=RP
1     CONTINUE
      RETURN
      END
C
C
      SUBROUTINE Y3P(IOSC,IW1,IW2,OVER,H)
C     --------------------(Y1*Y1*Y1*P1+P1*Y1*Y1*Y1)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 RMK(2),AAS(2)
      INTEGER NOBAS,MBASP1,LENIW
      INTEGER IOSC
      INTEGER LR,LC,IR,IC,JR,JC
      REAL*8 H(NOBAS,NOBAS),OVER(2,MBASP1,MBASP1),
     1      RQ,RN,RJ,RP1,RP2,RP3,AL,RP,RLN,RLS
      INTEGER IW1(LENIW),IW2(LENIW)
      COMMON/MORSE/RMK,AAS
      COMMON/MODIM/NOBAS,MBASP1,LENIW
C
      DO 1 LR=1,NOBAS
      IR=IW1(LR)
      JR=IW2(LR)
      DO 1 LC=1,NOBAS
      IC=IW1(LC)
      JC=IW2(LC)
      H(LR,LC)=0.0E+00
      IF((IR.EQ.IC).OR.(JR.NE.JC))GOTO 1
      RJ=ABS(DFLOAT(IR-IC))
      RN=RLN(IR,IC)
      AL=2.0D0*(RMK(IOSC)-RN) - 1.0D0
      RQ=RN+AL-RJ
      RP1=1.0D0-3.0D0*RJ*(AL-RJ)/(2.0D0*RMK(IOSC))
      RP1=RP1
     1     +3.0D0*(RJ**4+AL**2*RJ*(RJ+1.0D0)
     1     -AL*RJ**2*(2.0D0*RJ+1.0D0)
     1     +RJ*(AL-RJ)*(2.0D0*RN+1.0D0))/(8.0D0*RMK(IOSC)**2)
      RP2=(AL+3.0D0)*(0.5E0*RJ*(RJ-1.0D0)*(RQ+2.0D0)*(RQ+1.0D0)
     1     -RJ**2*RN*(RQ+1.0D0)+0.5E0*RJ*(RJ+1.0D0)*RN*(RN-1.0D0)
     1     +RN*RQ+2.0D0*RJ*(RQ+1.0D0)*(RN+RQ+2.0D0)
     1     -2.0D0*RN*RJ*(RN+RQ)+(RN+1.0D0)*(RQ+1.0D0)+(RN+RQ+1.0D0)**2)
      RP2=RP2
     1     +2.0D0*(0.5E0*RJ*(RJ-1.0D0)*(RQ+2.0D0)*(RQ+1.0D0)*RN
     1     -RJ**2*(RQ+1.0D0)*RN*(RN-1.0D0)+0.5E0*RJ*(RJ+1.0D0)*
     1     RN*(RN-1.0D0)*(RN-2.0D0)
     1     +RJ*(RQ+1.0D0)*RN*(2.0D0*RN+3.0D0*RQ+4.0D0)
     1     -RJ*RN*(RN-1.0D0)*(2.0D0*RN+3.0D0*RQ-1.0D0)
     1     +RN*(RN+1.0D0)*(RQ+1.0D0)
     1     +RN*(RN+2.0D0*RQ+1.0D0)*(RN+RQ+1.0D0)
     1     +RQ*(2.0D0*RN+RQ+1.0D0))
      RP2=RP2-
     1     (RJ*(RJ-1.0D0)*(RJ-2.0D0)*(RQ+3.0D0)*(RQ+2.0D0)*(RQ+1.0D0)
     1 /6.0D0
     1     -0.5E0*RJ**2*(RJ-1.0D0)*(RQ+2.0D0)*(RQ+1.0D0)*RN
     1     +0.5E0*RJ**2*(RJ+1.0D0)*(RQ+1.0D0)*RN*(RN-1.0D0)
     1     -RJ*(RJ+1.0D0)*(RJ+2.0D0)*RN*(RN-1.0D0)*(RN-2.0D0)/6.0D0
     1     +1.5E0*RJ*(RJ-1.0D0)*(RQ+2.0D0)*(RQ+1.0D0)*(RN+RQ+3.0D0))
      RP2=RP2-
     1     (-3.0D0*RJ**2*(RQ+1.0D0)*RN*(RN+RQ+1.0D0)
     1     +1.5E0*RJ*(RJ+1.0D0)*(RN-1.0D0)*RN*(RN+RQ-1.0D0)
     1     +RJ*(RQ+1.0D0)*((RN+2.0D0)*(RQ+2.0D0)
     1     +2.0D0*(RN+RQ+2.0D0)*(RN+RQ+3.0D0)+(RN+1.0D0)*(RQ+1.0D0)
     1     +(RN+RQ+1.0D0)**2+RN*RQ)
     1     -RJ*RN*((RN+1.0D0)*(RQ+1.0D0)+2.0D0*(RN+RQ)*(RN+RQ+1.0D0)
     1     +RN*RQ+(RN+RQ-1.0D0)**2+(RN-1.0D0)*(RQ-1.0D0))
     1     +2.0D0*(RN+1.0D0)*(RQ+1.0D0)*(RN+RQ+2.0D0)
     1     +(RN+RQ+1.0D0)*((RN+1.0D0)*(RQ+1.0D0)
     1     +(RN+RQ+1.0D0)**2+RN*RQ)+2.0D0*RN*RQ*(RN+RQ))
      RP1=RP1+RP2/(8.0D0*RMK(IOSC)**3)
      RP2=RP1*OVER(IOSC,IR+1,IC+1)*RLS(IR,IC)
      RP=RP2
      H(LR,LC)=RP
1     CONTINUE
      RETURN
      END
C
C
      SUBROUTINE YPY(IOSC1,IOSC2,IW1,IW2,OVER,H)
C     -----------------Y2*(P1*Y1+Y1*P1)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 RMK(2),AAS(2)
      INTEGER NOBAS,MBASP1,LENIW
      INTEGER IOSC1,IOSC2
      INTEGER LR,LC,IR,IC,JR,JC
      REAL*8 H(NOBAS,NOBAS),OVER(2,MBASP1,MBASP1),
     1      RN,RJ,RP1,RP2,RP3,AL,RP,RLN,RLS
      INTEGER IW1(LENIW),IW2(LENIW)
      COMMON/MORSE/RMK,AAS
      COMMON/MODIM/NOBAS,MBASP1,LENIW
C
      DO 1 LR=1,NOBAS
      IR=IW1(LR)
      JR=IW2(LR)
      DO 1 LC=1,NOBAS
      IC=IW1(LC)
      JC=IW2(LC)
      H(LR,LC)=0.0E+00
      IF(IR.EQ.IC)GOTO 1
      RJ=ABS(DFLOAT(IR-IC))
      RN=RLN(IR,IC)
      RP1=1.0D0-RJ+RJ*(2.0D0*RN+RJ+1.0D0)/(RMK(IOSC1)*2.0D0)
      RP1=RP1*OVER(IOSC1,IR+1,IC+1)*RLS(IR,IC)
      IF(JR.NE.JC)GOTO 2
      RP2=(DFLOAT(JR)+0.5E0)/RMK(IOSC2)
      GOTO 3
2     RP2=-OVER(IOSC2,JR+1,JC+1)/(2.0D0*AAS(IOSC2)*RMK(IOSC2))
3     RP=RP1*RP2
      H(LR,LC)=RP
1     CONTINUE
      RETURN
      END
C
C
      SUBROUTINE PY3(IOSC1,IOSC2,IW1,IW2,OVER,H)
C     --------------------(P1*Y2*Y2*Y2)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 RMK(2),AAS(2)
      INTEGER NOBAS,MBASP1,LENIW
      INTEGER IOSC1,IOSC2
      INTEGER LR,LC,IR,IC,JR,JC
      REAL*8 H(NOBAS,NOBAS),OVER(2,MBASP1,MBASP1),
     1      RN,RJ,RP1,RP2,RP,RLN,RLS
      INTEGER IW1(LENIW),IW2(LENIW)
      COMMON/MORSE/RMK,AAS
      COMMON/MODIM/NOBAS,MBASP1,LENIW
C
      DO 1 LR=1,NOBAS
      IR=IW1(LR)
      JR=IW2(LR)
      DO 1 LC=1,NOBAS
      IC=IW1(LC)
      JC=IW2(LC)
      H(LR,LC)=0.0E+00
      IF(IR.EQ.IC)GOTO 1
      RP1=0.5E0*OVER(IOSC1,IR+1,IC+1)*RLS(IR,IC)
      IF(JR.NE.JC)GOTO 2
      RN=DFLOAT(JR)
      RP2=1.5E0*(RN+0.5E0)*(RN+0.5E0)/(RMK(IOSC2)*RMK(IOSC2))
      RP2=RP2-1.0D0/(8.0D0*RMK(IOSC2)*RMK(IOSC2))
      RP2=RP2-RN*(2.0D0*RN*RN+3.0D0*RN+1.0D0)/(4.0D0*RMK(IOSC2)**3)
      GOTO 3
2     RJ=ABS(DFLOAT(JR-JC))
      RN=RLN(JR,JC)
      RP2=(RJ-1.0D0)*(2.0D0-RJ)/(4.0D0*AAS(IOSC2)*RMK(IOSC2))
      RP=(RJ-2.0D0)*(2.0D0*RJ+1.0D0)*(2.0D0*RN+RJ+1.0D0)
      RP2=(RP2+RP/(8.0D0*AAS(IOSC2)*RMK(IOSC2)**2))
     1                           *OVER(IOSC2,JR+1,JC+1)
3     RP=RP1*RP2
      H(LR,LC)=RP
1     CONTINUE
      RETURN
      END
C
C
      SUBROUTINE Y2PY(IOSC1,IOSC2,IW1,IW2,OVER,H)
C     --------------------(Y1*Y1*P1+P1*Y1*Y1)*Y2
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 RMK(2),AAS(2)
      INTEGER NOBAS,MBASP1,LENIW
      INTEGER IOSC1,IOSC2
      INTEGER LR,LC,IR,IC,JR,JC
      REAL*8 H(NOBAS,NOBAS),OVER(2,MBASP1,MBASP1),
     1      RN,RJ,RP1,RP2,RP3,AL,RP,RLN,RLS
      INTEGER IW1(LENIW),IW2(LENIW)
      COMMON/MORSE/RMK,AAS
      COMMON/MODIM/NOBAS,MBASP1,LENIW
C
      DO 1 LR=1,NOBAS
      IR=IW1(LR)
      JR=IW2(LR)
      DO 1 LC=1,NOBAS
      IC=IW1(LC)
      JC=IW2(LC)
      H(LR,LC)=0.0E+00
      IF(IR.EQ.IC)GOTO 1
      RJ=ABS(DFLOAT(IR-IC))
      RN=RLN(IR,IC)
      AL=(RMK(IOSC1)-RN)*2.0D0-1.0D0
      RP1=1.0D0-RJ*(AL-RJ)/RMK(IOSC1)
      RP=RJ**4+AL*AL*(RJ+1.0D0)*RJ
      RP3=RJ*(AL-RJ)*(2.0D0*RN+1.0D0)-AL*RJ*RJ*(2.0D0*RJ+1.0D0)
      RP1=(RP1+(RP+RP3)/(8.0D0*RMK(IOSC1)*RMK(IOSC1)))*
     1      OVER(IOSC1,IR+1,IC+1)*RLS(IR,IC)
      IF(JR.NE.JC)GOTO 2
      RP2=(DFLOAT(JR)+0.5E0)/RMK(IOSC2)
      GOTO 3
2     RP2=-OVER(IOSC2,JR+1,JC+1)/(2.0D0*AAS(IOSC2)*RMK(IOSC2))
3     RP=RP1*RP2
      H(LR,LC)=RP
1     CONTINUE
      RETURN
      END
C
C
      SUBROUTINE YPY2(IOSC1,IOSC2,IW1,IW2,OVER,H)
C     --------------------(P1*Y1+Y1*P1)*Y2*Y2
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 RMK(2),AAS(2)
      INTEGER NOBAS,MBASP1,LENIW
      INTEGER IOSC1,IOSC2
      INTEGER LR,LC,IR,IC,JR,JC
      REAL*8 H(NOBAS,NOBAS),OVER(2,MBASP1,MBASP1),
     1      RN,RJ,RP1,RP2,RP3,AL,RP,RLN,RLS
      INTEGER IW1(LENIW),IW2(LENIW)
      COMMON/MORSE/RMK,AAS
      COMMON/MODIM/NOBAS,MBASP1,LENIW
C
      DO 1 LR=1,NOBAS
      IR=IW1(LR)
      JR=IW2(LR)
      DO 1 LC=1,NOBAS
      IC=IW1(LC)
      JC=IW2(LC)
      H(LR,LC)=0.0E+00
      IF(IR.EQ.IC)GOTO 1
      RJ=ABS(DFLOAT(IR-IC))
      RN=RLN(IR,IC)
      RP1=1.0D0-RJ+RJ*(2.0D0*RN+RJ+1.0D0)/(RMK(IOSC1)*2.0D0)
      RP1=RP1*OVER(IOSC1,IR+1,IC+1)*RLS(IR,IC)
      IF(JR.NE.JC)GOTO 2
      RN=DFLOAT(JR)
      RP2=(RN+0.5E0)/RMK(IOSC2)
      GOTO 3
2     RN=RLN(JR,JC)
      RJ=ABS(DFLOAT(JR-JC))
      RP2=(RJ-1.0D0)/(2.0D0*AAS(IOSC2)*RMK(IOSC2))
      RP2=RP2-RJ*(2.0D0*RN+RJ+1.0D0)/(4.0D0*AAS(IOSC2)*RMK(IOSC2)**2)
      RP2=RP2*OVER(IOSC2,JR+1,JC+1)
3     RP=RP1*RP2
      H(LR,LC)=RP
1     CONTINUE
      RETURN
      END
C
C
      SUBROUTINE YONLY(IOSC,IW1,IW2,OVER,H)
C------------------------------(Y)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER IOSC
      REAL*8 RMK(2),AAS(2)
      INTEGER NOBAS,MBASP1,LENIW
      INTEGER LR,LC,IR,IC,JR,JC
      REAL*8 H(NOBAS,NOBAS),OVER(2,MBASP1,MBASP1),
     1      RN,RJ,RP1,RP2,RP3,AL,RP,RLN,RLS
      INTEGER IW1(LENIW),IW2(LENIW)
      COMMON/MORSE/RMK,AAS
      COMMON/MODIM/NOBAS,MBASP1,LENIW
C
      DO 1 LR=1,NOBAS
      IR=IW1(LR)
      JR=IW2(LR)
      DO 1 LC=1,NOBAS
      IC=IW1(LC)
      JC=IW2(LC)
      H(LR,LC)=0.0E+00
      IF(JR.NE.JC)GOTO 1
      IF(IR.NE.IC)GOTO 2
      RP1=(DFLOAT(IR)+0.5E0)/RMK(IOSC)
      GOTO 3
2     RP1=-OVER(IOSC,IR+1,IC+1)/(2.0D0*AAS(IOSC)*RMK(IOSC))
3     RP=RP1
      H(LR,LC)=RP
1     CONTINUE
      RETURN
      END
C
C
      SUBROUTINE YY(IOSC,IW1,IW2,OVER,H)
C     --------------------(Y*Y)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 RMK(2),AAS(2)
      INTEGER NOBAS,MBASP1,LENIW
      INTEGER IOSC
      INTEGER LR,LC,IR,IC,JR,JC
      REAL*8 H(NOBAS,NOBAS),OVER(2,MBASP1,MBASP1),
     1      RN,RJ,RP1,RP2,RP3,AL,RP,RLN,RLS
      INTEGER IW1(LENIW),IW2(LENIW)
      COMMON/MORSE/RMK,AAS
      COMMON/MODIM/NOBAS,MBASP1,LENIW
C
      DO 1 LR=1,NOBAS
      IR=IW1(LR)
      JR=IW2(LR)
      DO 1 LC=1,NOBAS
      IC=IW1(LC)
      JC=IW2(LC)
      H(LR,LC)=0.0E+00
      IF(JR.NE.JC)GOTO 1
      IF(IR.NE.IC)GOTO 2
      RN=DFLOAT(IR)
      RP1=(RN+0.5E0)/RMK(IOSC)
      GOTO 3
2     RN=RLN(IR,IC)
      RJ=ABS(DFLOAT(IR-IC))
      RP1=(RJ-1.0D0)/(2.0D0*AAS(IOSC)*RMK(IOSC))
      RP1=RP1-RJ*(2.0D0*RN+RJ+1.0D0)/(4.0D0*AAS(IOSC)
     1                            *RMK(IOSC)*RMK(IOSC))
      RP1=RP1*OVER(IOSC,IR+1,IC+1)
3     RP=RP1
      H(LR,LC)=RP
1     CONTINUE
      RETURN
      END
C
C
      SUBROUTINE Y1Y2(IOSC1,IOSC2,IW1,IW2,OVER,H)
C     --------------------(Y1*Y2)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 RMK(2),AAS(2)
      INTEGER NOBAS,MBASP1,LENIW
      INTEGER IOSC1,IOSC2
      INTEGER LR,LC,IR,IC,JR,JC
      REAL*8 H(NOBAS,NOBAS),OVER(2,MBASP1,MBASP1),
     1      RN,RJ,RP1,RP2,RP3,AL,RP,RLN,RLS
      INTEGER IW1(LENIW),IW2(LENIW)
      COMMON/MORSE/RMK,AAS
      COMMON/MODIM/NOBAS,MBASP1,LENIW
C
      DO 1 LR=1,NOBAS
      IR=IW1(LR)
      JR=IW2(LR)
      DO 1 LC=1,NOBAS
      IC=IW1(LC)
      JC=IW2(LC)
      H(LR,LC)=0.0E+00
      IF(IR.NE.IC)GOTO 2
      RP1=(DFLOAT(IR)+0.5E0)/RMK(IOSC1)
      GOTO 3
2     RP1=-OVER(IOSC1,IR+1,IC+1)/(2.0D0*AAS(IOSC1)*RMK(IOSC1))
3     IF(JR.NE.JC)GOTO 4
      RP2=(DFLOAT(JR)+0.5E0)/RMK(IOSC2)
      GOTO 5
4     RP2=-OVER(IOSC2,JR+1,JC+1)/(2.0D0*AAS(IOSC2)*RMK(IOSC2))
5     RP=RP1*RP2
      H(LR,LC)=RP
1     CONTINUE
      RETURN
      END
C
C
      SUBROUTINE YYY(IOSC,IW1,IW2,OVER,H)
C     --------------------(Y1*Y1*Y1)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 RMK(2),AAS(2)
      INTEGER NOBAS,MBASP1,LENIW
      INTEGER IOSC
      INTEGER LR,LC,IR,IC,JR,JC
      REAL*8 H(NOBAS,NOBAS),OVER(2,MBASP1,MBASP1),
     1      RN,RJ,RP1,RP2,RP3,AL,RP,RLN,RLS
      INTEGER IW1(LENIW),IW2(LENIW)
      COMMON/MORSE/RMK,AAS
      COMMON/MODIM/NOBAS,MBASP1,LENIW
C
      DO 1 LR=1,NOBAS
      IR=IW1(LR)
      JR=IW2(LR)
      DO 1 LC=1,NOBAS
      IC=IW1(LC)
      JC=IW2(LC)
      H(LR,LC)=0.0E+00
      IF(JR.NE.JC)GOTO 1
      IF(IR.NE.IC)GOTO 2
      RN=DFLOAT(IR)
      RP1=1.5E0*(RN+0.5E0)*(RN+0.5E0)/(RMK(IOSC)*RMK(IOSC))
      RP1=RP1-1.0D0/(8.0D0*RMK(IOSC)*RMK(IOSC))
      RP1=RP1-RN*(2.0D0*RN*RN+3.0D0*RN+1.0D0)/(4.0D0*RMK(IOSC)**3)
      GOTO 3
2     RJ=ABS(DFLOAT(IR-IC))
      RN=RLN(IR,IC)
      RP1=(RJ-1.0D0)*(2.0D0-RJ)/(4.0D0*AAS(IOSC)*RMK(IOSC))
      RP=(RJ-2.0D0)*(2.0D0*RJ+1.0D0)*(2.0D0*RN+RJ+1.0D0)
      RP1=(RP1+RP/(8.0D0*AAS(IOSC)*RMK(IOSC)**2))
     1                            *OVER(IOSC,IR+1,IC+1)
3     RP=RP1
      H(LR,LC)=RP
1     CONTINUE
      RETURN
      END
C
C
      SUBROUTINE Y1Y1Y2(IOSC1,IOSC2,IW1,IW2,OVER,H)
C     --------------------(Y1*Y1*Y2)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 RMK(2),AAS(2)
      INTEGER NOBAS,MBASP1,LENIW
      INTEGER IOSC1,IOSC2
      INTEGER LR,LC,IR,IC,JR,JC
      REAL*8 H(NOBAS,NOBAS),OVER(2,MBASP1,MBASP1),
     1      RN,RJ,RP1,RP2,RP3,AL,RP,RLN,RLS
      INTEGER IW1(LENIW),IW2(LENIW)
      COMMON/MORSE/RMK,AAS
      COMMON/MODIM/NOBAS,MBASP1,LENIW
C
      DO 1 LR=1,NOBAS
      IR=IW1(LR)
      JR=IW2(LR)
      DO 1 LC=1,NOBAS
      IC=IW1(LC)
      JC=IW2(LC)
      H(LR,LC)=0.0E+00
      IF(IR.NE.IC)GOTO 2
      RP1=(DFLOAT(IR)+0.5E0)/RMK(IOSC1)
      GOTO 3
2     RJ=ABS(DFLOAT(IR-IC))
      RN=RLN(IR,IC)
      RP1=(RJ-1.0D0)/(2.0D0*AAS(IOSC1)*RMK(IOSC1))
      RP=RJ*(2.0D0*RN+RJ+1.0D0)/(4.0D0*AAS(IOSC1)*RMK(IOSC1)
     1                                         *RMK(IOSC1))
      RP1=(RP1-RP)*OVER(IOSC1,IR+1,IC+1)
3     IF(JR.NE.JC)GOTO 4
      RP2=(DFLOAT(JR)+0.5E0)/RMK(IOSC2)
      GOTO 5
4     RP2=-OVER(IOSC2,JR+1,JC+1)/(2.0D0*AAS(IOSC2)*RMK(IOSC2))
5     RP=RP1*RP2
      H(LR,LC)=RP
1     CONTINUE
      RETURN
      END
C
C
      SUBROUTINE YYYY(IOSC,IW1,IW2,OVER,H)
C     --------------------(Y1*Y1*Y1*Y1)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 RMK(2),AAS(2)
      INTEGER NOBAS,MBASP1,LENIW
      INTEGER IOSC
      INTEGER LR,LC,IR,IC,JR,JC
      REAL*8 H(NOBAS,NOBAS),OVER(2,MBASP1,MBASP1),
     1      RN,RJ,RP1,RP2,RP3,AL,RP,RLN,RLS
      INTEGER IW1(LENIW),IW2(LENIW)
      COMMON/MORSE/RMK,AAS
      COMMON/MODIM/NOBAS,MBASP1,LENIW
C
      DO 1 LR=1,NOBAS
      IR=IW1(LR)
      JR=IW2(LR)
      DO 1 LC=1,NOBAS
      IC=IW1(LC)
      JC=IW2(LC)
      H(LR,LC)=0.0E+00
      IF(JR.NE.JC)GOTO 1
      IF(IR.NE.IC)GOTO 2
      RN=DFLOAT(IR)
      RP=(1.5E0*(RN+0.5E0)**2+3.0D0/8.0D0)/(RMK(IOSC)*RMK(IOSC))
      RP1=(RN*(2.0D0*RN*RN+3.0D0*RN+3.0D0)+1.0D0)/(4.0D0*RMK(IOSC)**3)
      RP1=RP-RP1
      GOTO 3
2     RN=RLN(IR,IC)
      RJ=ABS(DFLOAT(IR-IC))
      RP=(RJ-1.0D0)*(RJ-2.0D0)*(RJ-3.0D0)/(12.0D0*AAS(IOSC)
     1                                        *RMK(IOSC))
      RP1=(RJ-1.0D0)*(RJ*RJ-3.0D0*RJ-2.0D0)*(2.0D0*RN+RJ+1.0D0)
      RP1=(RP-RP1/(8.0D0*AAS(IOSC)*RMK(IOSC)*RMK(IOSC)))
     1                               *OVER(IOSC,IR+1,IC+1)
3     RP=RP1
      H(LR,LC)=RP
1     CONTINUE
      RETURN
      END
C
C
      SUBROUTINE Y12Y22(IOSC1,IOSC2,IW1,IW2,OVER,H)
C     --------------------(Y1*Y1*Y2*Y2)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 RMK(2),AAS(2)
      INTEGER NOBAS,MBASP1,LENIW
      INTEGER IOSC1,IOSC2
      INTEGER LR,LC,IR,IC,JR,JC
      REAL*8 H(NOBAS,NOBAS),OVER(2,MBASP1,MBASP1),
     1      RN,RJ,RP1,RP2,RP3,AL,RP,RLN,RLS
      INTEGER IW1(LENIW),IW2(LENIW)
      COMMON/MORSE/RMK,AAS
      COMMON/MODIM/NOBAS,MBASP1,LENIW
C
      DO 1 LR=1,NOBAS
      IR=IW1(LR)
      JR=IW2(LR)
      DO 1 LC=1,NOBAS
      IC=IW1(LC)
      JC=IW2(LC)
      H(LR,LC)=0.0E+00
      IF(IR.NE.IC)GOTO 2
      RP1=(DFLOAT(IR)+0.5E0)/RMK(IOSC1)
      GOTO 3
2     RJ=ABS(DFLOAT(IR-IC))
      RN=RLN(IR,IC)
      RP=(RJ-1.0D0)/(2.0D0*AAS(IOSC1)*RMK(IOSC1))
      RP1=RJ*(2.0D0*RN+RJ+1.0D0)/(4.0D0*AAS(IOSC1)*RMK(IOSC1)
     1                                         *RMK(IOSC1))
      RP1=(RP-RP1)*OVER(IOSC1,IR+1,IC+1)
3     IF(JR.NE.JC)GOTO 4
      RP2=(DFLOAT(JR)+0.5E0)/RMK(IOSC2)
      GOTO 5
4     RN=RLN(JR,JC)
      RJ=ABS(DFLOAT(JR-JC))
      RP=(RJ-1.0D0)/(2.0D0*AAS(IOSC2)*RMK(IOSC2))
      RP2=RJ*(2.0D0*RN+RJ+1.0D0)/(4.0D0*AAS(IOSC2)*RMK(IOSC2)
     1                                         *RMK(IOSC2))
      RP2=(RP-RP2)*OVER(IOSC2,JR+1,JC+1)
5     RP=RP1*RP2
      H(LR,LC)=RP
1     CONTINUE
      RETURN
      END
C
C
      SUBROUTINE Y13Y2(IOSC1,IOSC2,IW1,IW2,OVER,H)
C---------------(Y1*Y1*Y1*Y2)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 RMK(2),AAS(2)
      INTEGER NOBAS,MBASP1,LENIW
      INTEGER IOSC1,IOSC2
      INTEGER LR,LC,IR,IC,JR,JC
      REAL*8 H(NOBAS,NOBAS),OVER(2,MBASP1,MBASP1),
     1      RN,RJ,RP1,RP2,RP3,AL,RP,RLN,RLS
      INTEGER IW1(LENIW),IW2(LENIW)
      COMMON/MORSE/RMK,AAS
      COMMON/MODIM/NOBAS,MBASP1,LENIW
C
      DO 1 LR=1,NOBAS
      IR=IW1(LR)
      JR=IW2(LR)
      DO 1 LC=1,NOBAS
      IC=IW1(LC)
      JC=IW2(LC)
      H(LR,LC)=0.0E+00
      IF(IR.NE.IC)GOTO 2
      RN=DFLOAT(IR)
      RP1=1.5E0*(RN+0.5E0)*(RN+0.5E0)/(RMK(IOSC1)*RMK(IOSC1))
      RP1=RP1-1.0D0/(8.0D0*RMK(IOSC1)*RMK(IOSC1))
      RP1=RP1-RN*(2.0D0*RN*RN+3.0D0*RN+1.0D0)/(4.0D0*RMK(IOSC1)**3)
      GOTO 3
2     RJ=ABS(DFLOAT(IR-IC))
      RN=RLN(IR,IC)
      RP1=(RJ-1.0D0)*(2.0D0-RJ)/(4.0D0*AAS(IOSC1)*RMK(IOSC1))
      RP=(RJ-2.0D0)*(2.0D0*RJ+1.0D0)*(2.0D0*RN+RJ+1.0D0)
      RP1=(RP1+RP/(8.0D0*AAS(IOSC1)*RMK(IOSC1)*RMK(IOSC1)))*
     1      OVER(IOSC1,IR+1,IC+1)
3     IF(JR.NE.JC)GOTO 4
      RP2=(DFLOAT(JR)+0.5E0)/RMK(IOSC2)
      GOTO 5
4     RP2=-OVER(IOSC2,JR+1,JC+1)/(2.0D0*AAS(IOSC2)*RMK(IOSC2))
5     RP=RP1*RP2
      H(LR,LC)=RP
1     CONTINUE
      RETURN
      END
C
C
      SUBROUTINE PP(IOSC,IW1,IW2,OVER,H)
C     --------------------(P*P)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 RMK(2),AAS(2)
      INTEGER NOBAS,MBASP1,LENIW
      INTEGER IOSC
      INTEGER LR,LC,IR,IC,JR,JC
      REAL*8 H(NOBAS,NOBAS),OVER(2,MBASP1,MBASP1),
     1      RN,RJ,RP1,RP2,RP3,AL,RP,RLN,RLS
      INTEGER IW1(LENIW),IW2(LENIW)
      COMMON/MORSE/RMK,AAS
      COMMON/MODIM/NOBAS,MBASP1,LENIW
C
      DO 1 LR=1,NOBAS
      IR=IW1(LR)
      JR=IW2(LR)
      DO 1 LC=1,NOBAS
      IC=IW1(LC)
      JC=IW2(LC)
      H(LR,LC)=0.0E+00
      IF(JR.NE.JC)GOTO 1
      IF(IC.NE.IR)GOTO 2
      RN=DFLOAT(IC)
      RP=-AAS(IOSC)*AAS(IOSC)*(RN+0.5E0)*(2.0D0*(RMK(IOSC)-RN)
     1                                          -1.0D0)/2.0D0
      GOTO 3
2     RJ=ABS(DFLOAT(IC-IR))
      RN=RLN(IR,IC)
      RP=0.5E0*AAS(IOSC)*OVER(IOSC,IR+1,IC+1)*(RMK(IOSC)*(RJ-1.0D0)
     1     -0.5E0*RJ*(2.0D0*RN+RJ+1.0D0))
3     H(LR,LC)=RP
1     CONTINUE
      RETURN
      END
C
C
      SUBROUTINE P1P2(IOSC1,IOSC2,IW1,IW2,OVER,H)
C     --------------------(P1*P2)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 RMK(2),AAS(2)
      INTEGER NOBAS,MBASP1,LENIW
      INTEGER IOSC1,IOSC2
      INTEGER LR,LC,IR,IC,JR,JC
      REAL*8 H(NOBAS,NOBAS),OVER(2,MBASP1,MBASP1),
     1      RN,RJ,RP1,RP2,RP3,AL,RP,RLN,RLS
      INTEGER IW1(LENIW),IW2(LENIW)
      COMMON/MORSE/RMK,AAS
      COMMON/MODIM/NOBAS,MBASP1,LENIW
C
      DO 1 LR=1,NOBAS
      IR=IW1(LR)
      JR=IW2(LR)
      DO 1 LC=1,NOBAS
      IC=IW1(LC)
      JC=IW2(LC)
      H(LR,LC)=0.0E+00
      IF((IR.EQ.IC).OR.(JR.EQ.JC))GOTO 1
      RP1=0.5E0*OVER(IOSC1,IR+1,IC+1)*RLS(IR,IC)
      RP2=0.5E0*OVER(IOSC2,JR+1,JC+1)*RLS(JR,JC)
      RP=RP1*RP2
      H(LR,LC)=RP
1     CONTINUE
      RETURN
      END
C
      SUBROUTINE PYYPP(IOSC1,IOSC2,IW1,IW2,OVER,H)
C     --------------------(P1*Y1+Y1*P1)*P2)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 RMK(2),AAS(2)
      INTEGER NOBAS,MBASP1,LENIW
      INTEGER IOSC1,IOSC2
      INTEGER LR,LC,IR,IC,JR,JC
      REAL*8 H(NOBAS,NOBAS),OVER(2,MBASP1,MBASP1),
     1      RN,RJ,RP1,RP2,RP3,AL,RP,RLN,RLS
      INTEGER IW1(LENIW),IW2(LENIW)
      COMMON/MORSE/RMK,AAS
      COMMON/MODIM/NOBAS,MBASP1,LENIW
C
      DO 1 LR=1,NOBAS
      IR=IW1(LR)
      JR=IW2(LR)
      DO 1 LC=1,NOBAS
      IC=IW1(LC)
      JC=IW2(LC)
      H(LR,LC)=0.0E+00
      IF((IR.EQ.IC).OR.(JR.EQ.JC))GOTO 1
      RJ=ABS(DFLOAT(IR-IC))
      RN=RLN(IR,IC)
      RP1=1.0D0-RJ+RJ*(2.0D0*RN+RJ+1.0D0)/(RMK(IOSC1)*2.0D0)
      RP1=RP1*OVER(IOSC1,IR+1,IC+1)*RLS(IR,IC)
      RP2=0.5E0*OVER(IOSC2,JR+1,JC+1)*RLS(JR,JC)
      RP=RP1*RP2
      H(LR,LC)=RP
1     CONTINUE
      RETURN
      END
C
C
      SUBROUTINE YYPP(IOSC1,IOSC2,IW1,IW2,OVER,H)
C     --------------------((Y1*Y1*P1+P1*Y1*Y1)*P2)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 RMK(2),AAS(2)
      INTEGER NOBAS,MBASP1,LENIW
      INTEGER IOSC1,IOSC2
      INTEGER LR,LC,IR,IC,JR,JC
      REAL*8 H(NOBAS,NOBAS),OVER(2,MBASP1,MBASP1),
     1      RN,RJ,RP1,RP2,RP3,AL,RP,RLN,RLS
      INTEGER IW1(LENIW),IW2(LENIW)
      COMMON/MORSE/RMK,AAS
      COMMON/MODIM/NOBAS,MBASP1,LENIW
C
      DO 1 LR=1,NOBAS
      IR=IW1(LR)
      JR=IW2(LR)
      DO 1 LC=1,NOBAS
      IC=IW1(LC)
      JC=IW2(LC)
      H(LR,LC)=0.0E+00
      IF((IR.EQ.IC).OR.(JR.EQ.JC))GOTO 1
      RJ=ABS(DFLOAT(IR-IC))
      RN=RLN(IR,IC)
      AL=(RMK(IOSC1)-RN)*2.0D0-1.0D0
      RP1=1.0D0-RJ*(AL-RJ)/RMK(IOSC1)
      RP=RJ**4+AL*AL*(RJ+1.0D0)*RJ
      RP3=RJ*(AL-RJ)*(2.0D0*RN+1.0D0)-AL*RJ*RJ*(2.0D0*RJ+1.0D0)
      RP1=(RP1+(RP+RP3)/(8.0D0*RMK(IOSC1)*RMK(IOSC1)))*
     1      OVER(IOSC1,IR+1,IC+1)*RLS(IR,IC)
      RP2=0.5E0*OVER(IOSC2,JR+1,JC+1)*RLS(JR,JC)
      RP=RP1*RP2
      H(LR,LC)=RP
1     CONTINUE
      RETURN
      END
C
C
      SUBROUTINE YPPY(IOSC1,IOSC2,IW1,IW2,OVER,H)
C     -----------------------(Y1*P1+P1*Y1)*(Y2*P2+P2*Y2)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 RMK(2),AAS(2)
      INTEGER NOBAS,MBASP1,LENIW
      INTEGER IOSC1,IOSC2
      INTEGER LR,LC,IR,IC,JR,JC
      REAL*8 H(NOBAS,NOBAS),OVER(2,MBASP1,MBASP1),
     1      RN,RJ,RP1,RP2,RP3,AL,RP,RLN,RLS
      INTEGER IW1(LENIW),IW2(LENIW)
      COMMON/MORSE/RMK,AAS
      COMMON/MODIM/NOBAS,MBASP1,LENIW
C
      DO 1 LR=1,NOBAS
      IR=IW1(LR)
      JR=IW2(LR)
      DO 1 LC=1,NOBAS
      IC=IW1(LC)
      JC=IW2(LC)
      H(LR,LC)=0.0E+00
      IF((IR.EQ.IC).OR.(JR.EQ.JC))GOTO 1
      RJ=ABS(DFLOAT(IR-IC))
      RN=RLN(IR,IC)
      RP1=1.0D0-RJ+RJ*(2.0D0*RN+RJ+1.0D0)/(2.0D0*RMK(IOSC1))
      RP1=RP1*OVER(IOSC1,IR+1,IC+1)*RLS(IR,IC)
C
      RJ=ABS(DFLOAT(JR-JC))
      RN=RLN(JR,JC)
      RP2=1.0D0-RJ+RJ*(2.0D0*RN+RJ+1.0D0)/(2.0D0*RMK(IOSC2))
      RP2=RP2*OVER(IOSC2,JR+1,JC+1)*RLS(JR,JC)
      RP=RP1*RP2
      H(LR,LC)=RP
1     CONTINUE
      RETURN
      END
C
C
      DOUBLE PRECISION FUNCTION RLN(IR,IC)
      IMPLICIT REAL*8 (A-H,O-Z)
      IM=IR-IC
      IF(IM.LT.0)GOTO 1
      RLN=DFLOAT(IC)
      RETURN
1     RLN=DFLOAT(IR)
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION RLS(IR,IC)
      IMPLICIT REAL*8 (A-H,O-Z)
      IM=(IR-IC)/IABS(IR-IC)
      RLS=DFLOAT(IM)
      RETURN
      END
