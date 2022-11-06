      SUBROUTINE ECALC ( MWORK, NINPRE, NRECUN, WORK, IWRK,
     1                   IWRKST , NREC2 , SKPSYM , MXJVP1 , ISO ,
     1                   IBOOK , ENRGMX, COENMX )
C*********************************************************************
C
C       PROGRAM FOR CALCULATING THE ROTATION-VIBRATION ENERGIES OF
C       A TRIATOMIC MOLECULE EXHIBITING THE RENNER-TELLER EFFECT
C       AFTER THE MORBID SCHEME OF PER JENSEN
C
C*********************************************************************
C
C   THE FOLLOWING PARAMETERS ARE INITIALLY SET:
C
C       MWORK      :  LENGTH OF WORK ARRAY IN UNITS OF REALS
C       NINPRE     :  NUMBER OF INTEGERS OCCUPYING THE SAME STORAGE
C                  :  AS ONE REAL*8
C       NRECUN     :  NUMBER OF UNITS FOR RECORD LENGTH (FOR DIRECT
C                  :  ACCESS FILES) NEEDED TO OBTAIN A RECORD LENGTH
C                  :  OF ONE REAL*8
C
C*********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      INTEGER IBOOK ( MXJVP1 , 2 , 2 , 12 )
C
      REAL*8 M1,M2,M3,M,U1,U3,U13,V
C
      REAL*8 WORK(MWORK)
C
      REAL*8 AMASS(3,10)
C
      REAL*8 RHOE,FA1,FA2,FA3,FA4,FA5,FA6,FA7,FA8,
     1       AA1,AA3,
     2       F1A1,F2A1,F3A1,F4A1,F1A3,F2A3,F3A3,F4A3,
     3       F11,F1A11,F2A11,F3A11,F33,F1A33,F2A33,F3A33,
     4       F13,F1A13,F2A13,F3A13,
     5       F111,F1A111,F2A111,F333,F1A333,F2A333,
     6       F113,F1A113,F2A113,F133,F1A133,F2A133,
     7       F1111,FA1111,F3333,FA3333,F1113,FA1113,
     8       F1333,FA1333,F1133,FA1133,
     8       RE12 , RE32 , RHOREF , VMIN
C
      REAL*8 ETRIAL , RHOMAX , PNM1 , HBASE , HSTEP , EGUESS ,
     1      PREC
C
      REAL*8 THRSH1 , THRSH2 , THRSH3 , THRSH4 , THRSH5 ,
     1      THRSH6 , THRSH7 , THRSH8 , THRSH9 , THRSHX ,
     2      VELLGT , PLANCK , AVOGNO , DEGRAD , RADDEG ,
     3      PI
C
      REAL*8   B11,B13,B111,B133,B113,
     2      B1111,B1333,B1113,B1133,
     3      B11111,B13333,B11113,B11333,B11133,
     4      B31,B33,B311,B333,B313,
     5      B3111,B3333,B3113,B3133,
     6      B31111,B33333,B31113,B31333,B31133
C
      REAL*8   CR1,CR3,CR11,CR33,CR13,
     2      CR111,CR333,CR113,CR133,
     3      CR1111,CR3333,CR1113,CR1333,CR1133
C
      INTEGER NFIL1 , NFIL2 , NFIL3 , NFIL4 , NFIL5 , NFIL6 ,
     5       NFIL7 , NFIL8 , NFIL9 , NFIL10 ,
     6       NFIL11 , NFIL12 , NFIL13 , NFIL14 , NFIL15 ,
     7       NFIL16 , NFIL17 , NFIL18 , NFIL19 , NFIL20 ,
     6       ITEST  , IPRINT , NSTNR , NSTNIN , IREST ,
     7       IISOT , IQUAS , ISYMS , NISOT , NQUAS ,
     8       NUMQUA , NOPTIT , NOPTIM , IOBSER , NOBSER,
     9       ISOMAX , NATTS  , V0TYPE , IVAR(128) , PARMAX ,
     1       NUMPAR , PRTINT
C
      INTEGER V1 ,V2, V3 , V2MXP1, V2P1,
     1       NSTINT , NSERIN , NSERP , NSERQ , KQUA , NTEST ,
     2       NSEPP2 , NSEQP1 , MBASIS ,
     3       MDIM , NFSYM0, NFASY0, NFSYMJ, NFASYJ,
     4       KSTYPA(2) , LSTYPA(2) , JMAX , V2MAX , JMAXP1 ,
     5       V2MP11 , V2MP12
C
      INTEGER IQUANT(9,10)
C
      LOGICAL SYMM,POTSYM,ENRGOK
C
      REAL*8 RMK(2),AAS(2)
      REAL*8 ENRGMX ( ISOMAX )
      REAL*8 COENMX ( ISOMAX )
      INTEGER NOBAS,MBASP1,LENIW,JDIMST(4),LINDEX(4)
      LOGICAL SKPSYM ( 2 , 2 , MXJVP1 , ISOMAX )
      CHARACTER*130 ELINE
      CHARACTER*1 CABU
      CHARACTER*4 SYMSYM(4),ASYSYM(4)
      INTEGER IWRK(NINPRE*MWORK),IWRKST,IWRKEN,INTWRK,
     1       NRECS,NRSYM0,NRASY0,NRSYMJ,NRASYJ,
     2       ISSYM0,ISASY0,ISSYMJ,ISASYJ,IENRGY,
     3       IF1ST,IF2ST,IV0ST,IRRST,IGFST,IPHIST,IDERST,
     4       IJACST,IRHSST,IXTXST,IXTYST,ICSSST,IASSST,IIA1ST,
     5       IIA2ST,IRTIRR,LENREC,       IVALST,IPHIL,IPHIR,
     6       IDERR,IDERL,ISTOST,NPHI,LENPHI,LENVAL,IIW1ST,IIW3ST,
     7       IOVEST,IVMAST,IHMAST,NOVA,NOVB,IAMAST,IEVMST,
     8       IUMAST,IEWKST,I,NFCT,LENIAR,LENINT,IRTIST,
     9       ICOMST,IEVIST,IDOMST,IDIM,NORECS,
     1       ITMAST,IESTST,ISMAST,IEWRST,INDBST,NAMDIM,
     1       IHTRST,IWMAST,V2VALU,LZERSY(2),LZERAS(2),LONESY(2),
     2       LONEAS(2),LTWOSY(2),LTWOAS(2),LISTRT(2),LESTRT(2),
     3       LEBEND(2),LTMAST(2),LSPOSY(2),LSPOAS(2)
      REAL*8 RECS
C
C
      INTEGER IRO(4)
      INTEGER DIMR,DIMC
      INTEGER NABSYM,NABSR,NABSC
      CHARACTER*1 PSEUR,PSEUC
      LOGICAL AABB
      COMMON /AMPDI/ IDPMA,IDPTR
      INTEGER IDPMA,IDPTR
      COMMON /DIP/ IACT
      INTEGER IEIGMX,IEIGTP
      INTEGER INTFLG
      LOGICAL INTFST
      LOGICAL DPFL
      INTEGER IEGVEC(2,2,2,2)
      INTEGER IEIGST,IEIGAZ
      INTEGER NRAUS
C
      INTEGER ITRNSY
C
C
      include 'lzcomp.h'
      include 'isotop.h'
      include 'molcul.h'
      include 'value.h'
      include 'rentel.h'
      include 'integ.h'
      include 'dimen.h'
      include 'lsfit.h'
      include 'rensys.h'
      include 'spnorb.h'
      include 'bcoeff.h'
      include 'crcoef.h'
      include 'morse.h'
      include 'modim.h'
      include 'nummor.h'
C
      include 'intens.h'
C
C
      DATA SYMSYM/'A1  ','B2  ','B1  ','A2  '/
      DATA ASYSYM/'A''  ','A"  ','    ','    '/
      CALL CLOCKV ( OLDVEC , OLDTIM , 1 , 2 )
C
C
      ITRNSY=0
C
C
C*********************************************************************
C
C       SEQUENCE FOR CALCULATING THE RHO DEPENDENT FUNCTIONS
C       NECESSARY FOR THE RENNER CALCULATION.
C
C*********************************************************************
C
C       PREPARE TWO DIAGONAL MATRIX BLOCKS, ONE FOR EACH RENNER-
C       TELLER SURFACE
C
C*********************************************************************
C
C
C*********************************************************************
C
C       OPEN DIRECT ACCESS FILES FOR STORING BENDING WAVEFUNCTIONS
C
C*********************************************************************
C
      NUNIT=NFIL7
      OPEN(UNIT=NFIL7,STATUS='SCRATCH',ACCESS='DIRECT',
     1     RECL=2*NSTINT*NRECUN,ERR=1009)
      NUNIT=NFIL28
      OPEN(UNIT=NFIL28,STATUS='SCRATCH',ACCESS='DIRECT',
     1     RECL=2*NSTINT*NRECUN,ERR=1009)
C
C*********************************************************************
C
C       OPEN DIRECT ACCESS FILES FOR STORING MORSE OSCILLATOR
C       FUNCTIONS
C
C*********************************************************************
C
      NUNIT=NFIL37
      OPEN(UNIT=NFIL37,STATUS='SCRATCH',ACCESS='DIRECT',
     1     RECL=2*NSTSTR*NRECUN,ERR=1009)
C
      NUNIT=NFIL38
      OPEN(UNIT=NFIL38,STATUS='SCRATCH',ACCESS='DIRECT',
     1     RECL=2*NSTSTR*NRECUN,ERR=1009)
C
      NUNIT=NFIL39
      OPEN(UNIT=NFIL39,STATUS='SCRATCH',ACCESS='DIRECT',
     1     RECL=2*NSTSTR*NRECUN,ERR=1009)
C
      NUNIT=NFIL40
      OPEN(UNIT=NFIL40,STATUS='SCRATCH',ACCESS='DIRECT',
     1     RECL=2*NSTSTR*NRECUN,ERR=1009)
C
      REWIND NFIL33
C
C***********************************************************************
C
C       FIRST CALCULATE RHO DEPENDENT FUNCTIONS FOR BOTH SURFACES
C
C
C THE RHO DEPENDENT FUNCTIONS ARE STORED ON DISK AS FOLLOWS
C
C NFIL1: FUNCTIONS NECESSARY TO DO A J=0 CALCULATION FOR AN ABA
C        MOLECULE.
C NFIL2: ADDITIONAL FUNCTIONS NECESSARY TO DO A J=0 CALCULATION
C        FOR AN ABC MOLECULE.
C NFIL3: FUNCTIONS NECESSARY TO DO A J>0 CALCULATION FOR AN ABA
C        MOLECULE.
C NFIL4: ADDITIONAL FUNCTIONS NECESSARY TO DO A J>0 CALCULATION
C        FOR AN ABC MOLECULE.
C
C*********************************************************************
C
      NFSYM0=66
      KFSYM0=NFSYM0*NSTINT
      NFASY0=47
      KFASY0=NFASY0*NSTINT
      NFSYMJ=80
      KFSYMJ=NFSYMJ*NSTINT
      NFASYJ=58
      KFASYJ=NFASYJ*NSTINT
      NFINTA=15
      KFINTA=NFINTA*NSTINT
      INTSYM=58
      INTASY=49
C
C*********************************************************************
C
      IV0ST=IWRKST
      ISSYM0=IV0ST+NSTINT
      ISASY0=ISSYM0+KFSYM0
      ISSYMJ=ISASY0+KFASY0
      ISASYJ=ISSYMJ+KFSYMJ
      ISINTA=ISASYJ+KFASYJ
      IDGAST=ISINTA+KFINTA
      ICOSST=IDGAST+(JMAXP1*NSTINT)
      ISINST=ICOSST+(JMAXP1*NSTINT)
      IMUZST=ISINST+(JMAXP1*NSTINT)
      IWRKEN=IMUZST+(JMAXP1*NSTINT)
      IF (IWRKEN .GT. MWORK) GOTO 1001
C
C*********************************************************************
C
C RHOFCT CALCULATES THE APPROPRIATE RHO DEPENDENT FUNCTIONS AND
C STORES THEM ON UNITS NFIL1, NFIL2, NFIL3, AND NFIL4.
C
C*********************************************************************
C
      CALL RHOFCT ( WORK(IV0ST)  , WORK(ISSYM0) , WORK(ISASY0) ,
     1              WORK(ISSYMJ) , WORK(ISASYJ) , WORK(ISINTA) ,
     1              WORK(IDGAST) ,
     1              WORK(ICOSST) , WORK(ISINST) , WORK(IMUZST) ,
     1              KFSYM0 , KFASY0 , KFSYMJ , KFASYJ , KFINTA )
C
      CALL PRTIME( 'RHOFCT' , OLDTIM , OLDVEC )
C
C***********************************************************************
C
      REWIND NFIL9
      REWIND NFIL10
      REWIND NFIL29
C
      DO 9876 ISURF=1,2
C
      IF (IPRINT .NE. 0) WRITE (NFIL6,6650) ISURF
C
      NSURF=ISURF
C
      CALL CSHEET ( ISURF )
      CALL PRTIME( 'CSHEET' , OLDTIM , OLDVEC )
C
      IF (ISURF .EQ. 1) THEN
           NFILBN=NFIL7
           NFILBI=NFIL9
           NFILST=NFIL15
C
           NFILD0=NFIL17
           NFILD1=NFIL20
           REWIND NFIL17
           REWIND NFIL20
C
           V2MAX=IQUANT(2,ISO)
           MBASIS=IQUANT(3,ISO)
           LSTYPA(1)=IQUANT(4,ISO)
           LSTYPA(2)=IQUANT(5,ISO)
      ELSE
           NFILBN=NFIL28
           NFILBI=NFIL9
           NFILST=NFIL27
C
           NFILD0=NFIL19
           NFILD1=NFIL21
           REWIND NFIL19
           REWIND NFIL21
C
           V2MAX=IQUANT(6,ISO)
           MBASIS=IQUANT(7,ISO)
           LSTYPA(1)=IQUANT(8,ISO)
           LSTYPA(2)=IQUANT(9,ISO)
      ENDIF
      V2MXP1=V2MAX+1
      JMAXP1=JMAX+1
C
C
      CALL PRTIME( 'ZZZPOT' , OLDTIM , OLDVEC )
      CALL OPTPOT
      CALL PRTIME( 'OPTPOT' , OLDTIM , OLDVEC )
C
      CALL STROPT
      CALL PRTIME( 'STROPT' , OLDTIM , OLDVEC )
C
C*********************************************************************
C
C THE NF VARIABLES CONTAIN THE NUMBER OF FUNCTIONS WHOSE VALUES
C ARE TO BE STORED IN NFIL1, NFIL2, NFIL3, AND NFIL4, RESPECTIVELY.
C
C WE NOW DIVIDE UP THE WORK ARRAY TO HOLD THE VARIOUS VARIABLES.
C WE INITIALLY RESERVE SPACE FOR THE PURE BENDING ENERGIES
C WHICH WE SHALL CALCULATE IN THE SUBROUTINE WAVFUN
C
C*********************************************************************
C
C WE NOW RESERVE SPACE IN THE WORK ARRAY FOR THE PURE BENDING
C ENERGIES AND FOR THE FUNCTIONS NECESSARY FOR THE NUMEROV-COOLEY
C NUMERICAL INTEGRATION. THESE FUNCTIONS ARE: F1, F2, V0+UBEND,
C IRR0, AND G. FOR DEFINITIONS SEE P. JENSEN, COMP. PHYS. REP. 1,
C 1-55 (1983).
C
C*********************************************************************
C
      IENRGY=IWRKST
      IV0ST=IENRGY+V2MXP1*JMAXP1
      IF1ST=IV0ST+NSTINT
      IF2ST=IF1ST+NSTINT
      IRRST=IF2ST+NSTINT
      IGFST=IRRST+NSTINT
      IWRKEN=IGFST+NSTINT
      IF (IWRKEN .GE. MWORK) GOTO 1002
C
C*********************************************************************
C
C NUMFCT CALCULATES THE VALUES OF THESE FUNCTIONS
C
C*********************************************************************
C
      CALL NUMFCT ( WORK(IV0ST) , WORK(IF1ST) , WORK(IF2ST) ,
     1             WORK(IRRST) , WORK(IGFST) )
      CALL PRTIME( 'NUMFCT' , OLDTIM , OLDVEC )
C
C*********************************************************************
C
C
C*********************************************************************
C
C APART FROM THE SPACE ALREADY RESERVED IN THE WORK ARRAY
C WE NOW RESERVE SPACE FOR ONE BENDING WAVEFUNCTION (THE ONE
C THAT IS CURRENTLY BEING CALCULATED), FOR ITS DERIVATIVE
C AND FOR VARIOUS ARRAYS USED IN CARRYING OUT THE SERIES
C SOLUTION AROUND RHO=0.
C
C*********************************************************************
C
      IRTIST=IWRKEN
      IVPVST=IRTIST+NSTINT
      IVMVST=IVPVST+NSTINT
      IPHIST=IVMVST+NSTINT
      IDERST=IPHIST+NSTINT
      NSEPP2=NSERP+2
      NSEQP1=NSERQ+1
      IJACST=IDERST+NSTINT
      IRHSST=IJACST+NSERIN*NSEPP2
      IXTXST=IRHSST+NSERIN
      IXTYST=IXTXST+NSEPP2*NSEPP2
      ICSSST=IXTYST+NSEPP2
      IASSST=ICSSST+NSEPP2
      IIA1ST=IASSST+NSEQP1
      LIA1ST=ILCONV(IIA1ST,NINPRE)
      LENIAR=NSEPP2
      IIA2ST=IIA1ST+LENIAR
      LIA2ST=ILCONV(IIA2ST,NINPRE)
      ICOSGM=IIA2ST+LENIAR
      ISINGM=ICOSGM+(NSTINT*JMAXP1)
      IWRKEN=ISINGM+(NSTINT*JMAXP1)
      IF (IWRKEN .GT. MWORK) GOTO 1003
C
C*********************************************************************
C
C WAVFUN CALCULATES THE BENDING FUNCTIONS AND STORES THEM ON
C UNIT NFILBN.
C
C*********************************************************************
C
      CALL WAVFUN ( WORK(IENRGY) , WORK(IV0ST)  , WORK(IF1ST)  ,
     1              WORK(IF2ST)  , WORK(IRRST)  , WORK(IGFST)  ,
     2              WORK(IRTIST) , WORK(IVPVST) , WORK(IVMVST) ,
     2              WORK(IPHIST) , WORK(IDERST) ,
     3              WORK(IJACST) , WORK(IRHSST) , WORK(IXTXST) ,
     4              WORK(IXTYST) , WORK(ICSSST) , WORK(IASSST) ,
     5              IWRK(LIA1ST) , IWRK(LIA2ST) ,
     6              WORK(ICOSGM) , WORK(ISINGM) ,
     6              LENIAR       , NFILBN       )
      CALL PRTIME( 'WAVFUN' , OLDTIM , OLDVEC )
C
      LEBEND(ISURF)=IV0ST-IENRGY
      WRITE (NFIL33) (WORK(II), II=IENRGY,IV0ST-1)
C
C*********************************************************************
C
C WE NOW SET UP FOR CALCULATING THE INTEGRALS OVER THE BENDING
C COORDINATE NECESSARY FOR DOING A J=0 CALCULATION.
C WE SET UP SPACE FOR
C                       - THE INTEGRALS.
C                       - AS MANY VALUES OF THE RHO DEPENDENT
C                         THAT ONE NEEDS AT ANY GIVEN TIME.
C                       - A LEFT AND A RIGHT WAVEFUNCTION AND
C                         THEIR DERIVATIVES.
C
C*********************************************************************
C
      NZERSY=76
      NZERAS=53
      NONESY=26
      NONEAS=22
      NTWOSY=18
      NTWOAS=12
      NSPOSY=9
      NSPOAS=6
C
      NDIPSY=9
      NDIPAS=6
C
      MBASP1=MBASIS+1
      NOBAS=MBASP1*(MBASP1+1)/2
      IESTST=IPHIST
      IROTSY=IESTST+NOBAS
      IROTAS=IROTSY+KFSYMJ
      IVALSY=IROTAS+KFASYJ
      IVALAS=IVALSY+KFSYM0
      ICOSGM=IVALAS+KFASY0
      ISINGM=ICOSGM+(NSTINT*JMAXP1)
      IDGMDR=ISINGM+(NSTINT*JMAXP1)
      IAMUZZ=IDGMDR+(NSTINT*JMAXP1)
      IAMURR=IAMUZZ+NSTINT
      IPHIL=IAMURR+NSTINT
      IDERL=IPHIL+NSTINT
      IMUZZ=IDERL+NSTINT
      IPHIR=IMUZZ+NSTINT
      IDERR=IPHIR+NSTINT
      ILZS1=IDERR+NSTINT
      ILZS2=ILZS1+NSTINT
      ILZS3=ILZS2+NSTINT
      IACCUM=ILZS3+NSTINT
      IACCSO=IACCUM+NZERSY*JMAXP1*V2MXP1*V2MXP1
      IF (NZERAS .GT. NZERSY)
     1    IACCSO=IACCUM+NZERAS*JMAXP1*V2MXP1*V2MXP1
      LENACC=IACCSO-IACCUM
      IWRKEN=IACCSO+NSPOSY*JMAXP1*V2MXP1*V2MXP1
      IF (NSPOAS .GT. NSPOSY)
     1    IWRKEN=IACCSO+NSPOAS*JMAXP1*V2MXP1*V2MXP1
      LENASO=IWRKEN-IACCSO
C
C     Platz fuer Dipolmomente
C
      IDIPAC=IWRKEN
      IF (NDIPSY.GT.NDIPAS) THEN
         IWRKEN=IDIPAC+ NDIPSY*JMAXP1*V2MXP1*V2MXP1
      ELSE
         IWRKEN=IDIPAC+ NDIPAS*JMAXP1*V2MXP1*V2MXP1
      ENDIF
      LENDIP=IWRKEN-IDIPAC
      IF (IWRKEN .GT. MWORK) GOTO 1004
C
C*********************************************************************
C
C FIRST WE DO THE DELTA K=0 BENDING MATRIX ELEMENTS
C
C*********************************************************************
C
      CALL INTZER ( WORK(IGFST) , WORK(IROTSY) , WORK(IROTAS) ,
     1             WORK(IVALSY) , WORK(IVALAS) , WORK(ICOSGM) ,
     1             WORK(ISINGM) , WORK(IDGMDR) , WORK(IAMUZZ) ,
     1   WORK(IAMURR) , WORK(IPHIL)  , WORK(IDERL)  , WORK(IMUZZ),
     1             WORK(IPHIR)  , WORK(IDERR)  , WORK(ILZS1)  ,
     1             WORK(ILZS2)  , WORK(ILZS3)  ,
     1             WORK(IACCUM) , WORK(IACCSO) , WORK(IDIPAC),
     4             NFILBI , NFILBN , NFILD0,
     1             KFSYM0 , KFASY0 , KFSYMJ , KFASYJ ,
     1             LENACC , LENASO , LENDIP,
     1             LZERSY(ISURF) , LZERAS(ISURF),
     1             NZERSY , NZERAS , 
     1             NDIPSY,NDIPAS,
     1             LSPOSY(ISURF), LSPOAS(ISURF),
     1             NSPOSY , NSPOAS )
      CALL PRTIME( 'INTZER' , OLDTIM , OLDVEC )
C
C*********************************************************************
C
C THEN WE DO THE DELTA K=2 BENDING MATRIX ELEMENTS
C
C*********************************************************************
C
      ICOSGM=IVALSY
      ISINGM=ICOSGM+(NSTINT*JMAXP1)
      IPHIL=ISINGM+(NSTINT*JMAXP1)
      IDERL=IPHIL+NSTINT
      IPHIR=IDERL+NSTINT
      IDERR=IPHIR+NSTINT
      IACCUM=IDERR+NSTINT
      IWRKEN=IACCUM+2*NTWOSY*JMAXP1*V2MXP1*V2MXP1
      IF (NTWOAS .GT. NTWOSY)
     1    IWRKEN=IACCUM+2*NTWOAS*JMAXP1*V2MXP1*V2MXP1
      LENACC=IWRKEN-IACCUM
      IF (IWRKEN .GT. MWORK) GOTO 1016
C
      CALL INTTWO ( WORK(IROTSY) , WORK(IROTAS) , WORK(ICOSGM) ,
     1              WORK(ISINGM) ,
     1              WORK(IPHIL)  , WORK(IDERL)  , WORK(IPHIR)  ,
     2              WORK(IDERR)  , WORK(IACCUM) ,
     4             NFILBI ,          NFILBN ,
     1             KFSYMJ , KFASYJ ,
     1             LENACC , LTWOSY(ISURF) , LTWOAS(ISURF),
     1             NTWOSY , NTWOAS )
      CALL PRTIME( 'INTTWO' , OLDTIM , OLDVEC )
C
C*********************************************************************
C
C FINALLY WE DO THE DELTA K=1 BENDING MATRIX ELEMENTS
C
C*********************************************************************
C
      IDGKDR=IACCUM
      ILZS3=IDGKDR+JMAXP1*NSTINT
      IACCUM=ILZS3+NSTINT
      IWRKEN=IACCUM+2*NONESY*JMAXP1*V2MXP1*V2MXP1
      IF (NONEAS .GT. NONESY)
     1    IWRKEN=IACCUM+2*NONEAS*JMAXP1*V2MXP1*V2MXP1
      LENACC=IWRKEN-IACCUM
C
C     Platz fuer Dipolmomente 
C
      IDIPAC=IWRKEN
      IF (NDIPSY.GT.NDIPAS) THEN
         IWRKEN=IDIPAC + 2*NDIPSY*JMAXP1*V2MXP1*V2MXP1
      ELSE
         IWRKEN=IDIPAC + 2*NDIPAS*JMAXP1*V2MXP1*V2MXP1
      ENDIF
      LENDIP=IWRKEN-IDIPAC
      IF (IWRKEN .GT. MWORK) GOTO 1015
C
      CALL INTONE ( WORK(IROTSY) , WORK(IROTAS) , WORK(ICOSGM) ,
     1              WORK(ISINGM) , WORK(IDGKDR) ,
     1              WORK(IPHIL)  , WORK(IDERL)  , WORK(IPHIR)  ,
     2              WORK(IDERR)  , WORK(ILZS3)  , WORK(IACCUM) ,
     4             NFILBI ,          NFILBN ,
     1             KFSYMJ , KFASYJ ,
     1             LENACC , LONESY(ISURF) , LONEAS(ISURF),
     1             NONESY , NONEAS,
     2             WORK(IDIPAC),LENDIP,
     3             NDIPSY,NDIPAS,
     4             NFILD1 )
      CALL PRTIME( 'INTONE' , OLDTIM , OLDVEC )
C
C*********************************************************************
C
C SET UP PARAMETERS FOR THE CALCULATION OF MORSE OSCILLATOR MATRIX
C ELEMENTS
C
C*********************************************************************
C
      LENIAR=NOBAS
      MDIM=V2MXP1*(LSTYPA(1)+LSTYPA(2))
      IIW1ST=IROTSY
      LIW1ST=ILCONV(IIW1ST,NINPRE)
      IIW3ST=IIW1ST+LENIAR
      LIW3ST=ILCONV(IIW3ST,NINPRE)
      IOVEST=IIW3ST+LENIAR
      LENIW=LENIAR
      IWRKEN=IOVEST+2*MBASP1*MBASP1
      IF (IWRKEN .GE. MWORK) GOTO 1005
C
C*********************************************************************
C
      CALL MOPARM ( IWRK(LIW1ST) , IWRK(LIW3ST) , WORK(IOVEST) ,
     1              RK11 , RK31 , RK12 , RK32 ,
     1              KSTY11 , KSTY21 , KSTY12 , KSTY22 )
      CALL PRTIME( 'MOPARM' , OLDTIM , OLDVEC )
C
C*********************************************************************
C
      NOPERS=38
      ITMAST=IWRKEN
      ISMAST=ITMAST+NOBAS*NOBAS
      IVMAST=ISMAST+NOBAS*NOBAS*NOPERS
      IUMAST=IVMAST+NOBAS*NOBAS
      IEWRST=IUMAST+NOBAS*2
      IHTRST=IEWRST+NOBAS
      IEIGVS=IHTRST+NOBAS
      IWSPST=IEIGVS+NOBAS*NOBAS
      IIWRKS=IWSPST+8*NOBAS
      LIWRKS=ILCONV(IIWRKS,NINPRE)
      IFAILS=IIWRKS+5*NOBAS
      LFAILS=ILCONV(IFAILS,NINPRE)
      ICOMST=IFAILS+NOBAS
      LCOMST=ILCONV(ICOMST,NINPRE)
      INDBST=ICOMST+LENIAR
      LNDBST=ILCONV(INDBST,NINPRE)
      IDOMST=INDBST+LENIAR
      LDOMST=ILCONV(IDOMST,NINPRE)
      IWRKEN=IDOMST+LENIAR
      IF (IWRKEN .GE. MWORK) GOTO 1006
C
C*********************************************************************
C
      CALL STRTCH ( IWRK(LIW1ST) , IWRK(LIW3ST) , WORK(IOVEST) ,
     1             WORK(ITMAST) , WORK(IESTST) , WORK(ISMAST) ,
     1             WORK(IVMAST) , WORK(IUMAST) , WORK(IEWRST) ,
     1             WORK(IHTRST) , WORK(IEIGVS) , WORK(IWSPST) ,
     1             IWRK(LIWRKS) , IWRK(LFAILS) , IWRK(LCOMST) , 
     1             IWRK(LNDBST) , IWRK(LDOMST) )
      CALL PRTIME( 'STRTCH' , OLDTIM , OLDVEC )
C
      LESTRT(ISURF)=IROTSY-IESTST
      WRITE (NFIL33) (WORK(II), II=IESTST,IROTSY-1)
C
C*********************************************************************
C
      LTMAST(ISURF)=ISMAST-ITMAST
      WRITE (NFIL30) (WORK(III), III=ITMAST,ISMAST-1)
      IWMAST=IUMAST
      IWRKEN=IWMAST+NOBAS*NOBAS
      IF (IWRKEN .GE. MWORK) GOTO 1007
C
C*********************************************************************
C
      CALL GENSTR( IWRK(LIW1ST) , IWRK(LIW3ST) , WORK(IOVEST) ,
     1             WORK(ITMAST) , WORK(IESTST) , WORK(ISMAST) ,
     1             WORK(IVMAST) , WORK(IWMAST) , NOPERS )
      CALL PRTIME( 'GENSTR' , OLDTIM , OLDVEC )
C
C*********************************************************************
C
C
      LISTRT(ISURF)=IVMAST-ISMAST
      WRITE (NFIL10) (WORK(III), III=ISMAST,IVMAST-1)
C
C*********************************************************************
C
9876  CONTINUE
C
C*********************************************************************
C
C     CALCULATION OF DIAGONAL BLOCKS ENDS HERE; WE CONTINUE
C     WITH THE BLOCK COUPLING THE TWO ELECTRONIC STATES
C
C*********************************************************************
C
      IF (IPRINT .NE. 0) WRITE (NFIL6,6651)
C
C*********************************************************************
C
C  CALCULATE BENDING INTEGRALS NECESSARY FOR RENNER-TELLER
C  INTERACTION
C
C*********************************************************************
C
      NINSYM=45
      NINASY=33
C
      NFAB0=NFIL22
      NFAB1=NFIL23
      REWIND NFIL22
      REWIND NFIL23 
C
      IROTAS=IROTSY+KFSYMJ
      IVALSY=IROTAS+KFASYJ
      IVALAS=IVALSY+KFSYM0
      ICOSGM=IVALAS+KFASY0
      ISINGM=ICOSGM+(NSTINT*JMAXP1)
      IDGMDR=ISINGM+(NSTINT*JMAXP1)
      IPHIL=IDGMDR+(NSTINT*JMAXP1)
      IDERL=IPHIL+NSTINT
      IPHIR=IDERL+NSTINT
      IDERR=IPHIR+NSTINT
      ILZS1=IDERR+NSTINT
      ILZS2=ILZS1+NSTINT
      ILZS3=ILZS2+NSTINT
      IACCUM=ILZS3+NSTINT
      V2MP11=IQUANT(2,ISO)+1
      V2MP12=IQUANT(6,ISO)+1
      IACCSO=IACCUM+NINSYM*JMAXP1*V2MP11*V2MP12
      IF (NINASY .GT. NINSYM)
     1    IACCSO=IACCUM+NINASY*JMAXP1*V2MP11*V2MP12
      LENACC=IACCSO-IACCUM
      IWRKEN=IACCSO+NSPOSY*JMAXP1*V2MP11*V2MP12
      IF (NSPOAS .GT. NSPOSY)
     1    IWRKEN=IACCSO+NSPOAS*JMAXP1*V2MP11*V2MP12
      LENASO=IWRKEN-IACCSO
      IF (IWRKEN .GE. MWORK) GOTO 1012
C
C     Platz fuer Dipolmomente
C
      IDIPAC = IWRKEN
      IF (NDIPSY.GT.NDIPAS) THEN
        IWRKEN = IDIPAC + JMAXP1*V2MP11*V2MP12*NDIPSY
      ELSE
        IWRKEN = IDIPAC + JMAXP1*V2MP11*V2MP12*NDIPAS
      ENDIF
      LENDIP=IWRKEN-IDIPAC 
      IF (IWRKEN .GE. MWORK) GOTO 1012    
C
C*********************************************************************
C
C FIRST WE DO THE DELTA K=0 BENDING MATRIX ELEMENTS
C
C*********************************************************************
C
      CALL INTPM0 ( WORK(IROTSY) , WORK(IROTAS) ,
     1             WORK(IVALSY) , WORK(IVALAS) , WORK(ICOSGM) ,
     1             WORK(ISINGM) , WORK(IDGMDR) ,
     1             WORK(IPHIL)  , WORK(IDERL)  , WORK(IPHIR)  ,
     2             WORK(IDERR)  , WORK(ILZS1)  , WORK(ILZS2)  ,
     1             WORK(ILZS3)  , WORK(IACCUM) , WORK(IACCSO) ,
     4             NFILBI , V2MP11 , V2MP12 ,
     1             KFSYM0 , KFASY0 , KFSYMJ , KFASYJ ,
     1             LENACC , IZERSY , IZERAS,
     1             NINSYM , NINASY , LENASO, ISPOSY , ISPOAS ,
     1             NSPOSY , NSPOAS,
     2             WORK(IDIPAC),LENDIP,
     3             NFAB0,
     4             NDIPSY,NDIPAS )
      CALL PRTIME( 'INTPM0' , OLDTIM , OLDVEC )
C
C*********************************************************************
C
C THEN WE DO THE DELTA K=2 BENDING MATRIX ELEMENTS
C
C*********************************************************************
C
      ICOSGM=IVALSY
      ISINGM=ICOSGM+(NSTINT*JMAXP1)
      IPHIL=ISINGM+(NSTINT*JMAXP1)
      IDERL=IPHIL+NSTINT
      IPHIR=IDERL+NSTINT
      IDERR=IPHIR+NSTINT
      IACCUM=IDERR+NSTINT
      IWRKEN=IACCUM+2*NTWOSY*JMAXP1*V2MP11*V2MP12
      IF (NTWOAS .GT. NTWOSY)
     1    IWRKEN=IACCUM+2*NTWOAS*JMAXP1*V2MP11*V2MP12
      LENACC=IWRKEN-IACCUM
      IF (IWRKEN .GT. MWORK) GOTO 1016
C
      CALL INTPM2 ( WORK(IROTSY) , WORK(IROTAS) , WORK(ICOSGM) ,
     1              WORK(ISINGM) ,
     1              WORK(IPHIL)  , WORK(IDERL)  , WORK(IPHIR)  ,
     2              WORK(IDERR)  , WORK(IACCUM) ,
     4              NFILBI , V2MP11 , V2MP12 ,
     1             KFSYMJ , KFASYJ ,
     1             LENACC , ITWOSY , ITWOAS,
     1             NTWOSY , NTWOAS )
      CALL PRTIME( 'INTPM2' , OLDTIM , OLDVEC )
C
C*********************************************************************
C
C FINALLY WE DO THE DELTA K=1 BENDING MATRIX ELEMENTS
C
C*********************************************************************
C
      IDGKDR=IACCUM
      ILZS3=IDGKDR+JMAXP1*NSTINT
      IACCUM=ILZS3+NSTINT
      IWRKEN=IACCUM+2*NONESY*JMAXP1*V2MP11*V2MP12
      IF (NONEAS .GT. NONESY)
     1    IWRKEN=IACCUM+2*NONEAS*JMAXP1*V2MP11*V2MP12
      LENACC=IWRKEN-IACCUM
C
C     Platz fuer Dipolmomente
C
      IDIPAC=IWRKEN
      IF (NDIPSY.GT.NDIPAS) THEN 
         IWRKEN=IDIPAC + 2*NDIPSY*JMAXP1*V2MP11*V2MP12
      ELSE
         IWRKEN=IDIPAC + 2*NDIPAS*JMAXP1*V2MP11*V2MP12
      ENDIF
      LENDIP=IWRKEN-IDIPAC
      IF (IWRKEN .GT. MWORK) GOTO 1015
C
      CALL INTPM1 ( WORK(IROTSY) , WORK(IROTAS) , WORK(ICOSGM) ,
     1              WORK(ISINGM) , WORK(IDGKDR) ,
     1              WORK(IPHIL)  , WORK(IDERL)  , WORK(IPHIR)  ,
     2              WORK(IDERR)  , WORK(ILZS3)  , WORK(IACCUM) ,
     4              NFILBI , V2MP11 , V2MP12 ,
     1             KFSYMJ , KFASYJ ,
     1             LENACC , IONESY , IONEAS,
     1             NONESY , NONEAS,
     1             WORK(IDIPAC),LENDIP,
     2             NFAB1,
     3             NDIPSY,NDIPAS )
      CALL PRTIME( 'INTPM1' , OLDTIM , OLDVEC )
C
C*********************************************************************
C
C  CALCULATE TABLE VALUES OF MORSE OSCILLATOR EIGENFUNCTIONS
C
C*********************************************************************
C
      IPHIST=IROTSY
      IDERST=IPHIST+NSTSTR
      IWRKEN=IDERST+NSTSTR
      IF (IWRKEN .GE. MWORK) GOTO 1013
C
      MBASI1=IQUANT(3,ISO)
      MBASI2=IQUANT(7,ISO)
      CALL GENMOS ( WORK(IPHIST) , WORK(IDERST),
     1              RK11 , RK31 , RK12 , RK32 ,
     1              MBASI1 , MBASI2 )
      CALL PRTIME( 'GENMOS' , OLDTIM , OLDVEC )
C
C*********************************************************************
C
C  CALCULATE STRETCHING MATRIX ELEMENTS FOR RENNER-TELLER
C  INTERACTION
C
C*********************************************************************
C
      NINTOP=33
      MBALP1=MBASI1+1
      MBARP1=MBASI2+1
      NOVIB2=(MBASI1+1)*(MBASI2+1)
      NOBAL=(MBASI1+1)*(MBASI1+2)/2
      NOBAR=(MBASI2+1)*(MBASI2+2)/2
      NOVIBS = MAX(MBALP1,MBARP1)
      NOBAS  = MAX(NOBAR,NOBAL)
      IWFRS1=IROTSY
      IWFLS1=IWFRS1+NSTSTR
      IWFRS3=IWFLS1+NSTSTR
      IWFLS3=IWFRS3+NSTSTR
      IDFRS1=IWFLS3+NSTSTR
      IDFLS1=IDFRS1+NSTSTR
      IDFRS3=IDFLS1+NSTSTR
      IDFLS3=IDFRS3+NSTSTR
      IY1P1S=IDFLS3+NSTSTR
      IY1P2S=IY1P1S+NSTSTR
      IY1P3S=IY1P2S+NSTSTR
      IY1P4S=IY1P3S+NSTSTR
      IY3P1S=IY1P4S+NSTSTR
      IY3P2S=IY3P1S+NSTSTR
      IY3P3S=IY3P2S+NSTSTR
      IY3P4S=IY3P3S+NSTSTR
      IY1P0M=IY3P4S+NSTSTR
      IY1P1M=IY1P0M+NOVIB2
      IY1P2M=IY1P1M+NOVIB2
      IY1P3M=IY1P2M+NOVIB2
      IY1P4M=IY1P3M+NOVIB2
      IY3P0M=IY1P4M+NOVIB2
      IY3P1M=IY3P0M+NOVIB2
      IY3P2M=IY3P1M+NOVIB2
      IY3P3M=IY3P2M+NOVIB2
      IY3P4M=IY3P3M+NOVIB2
      IPY1P0=IY3P4M+NOVIB2
      IPY1P1=IPY1P0+NOVIB2
      IPY1P2=IPY1P1+NOVIB2
      IPY1P3=IPY1P2+NOVIB2
      IPY3P0=IPY1P3+NOVIB2
      IPY3P1=IPY3P0+NOVIB2
      IPY3P2=IPY3P1+NOVIB2
      IPY3P3=IPY3P2+NOVIB2
      IVMAST=IPY3P3+NOVIB2
      ITMASL=IVMAST+NOBAL*NOBAR
      ITMASR=ITMASL+NOBAL*NOBAL
      ISMAST=ITMASR+NOBAR*NOBAR
      IWMAST=ISMAST+NINTOP*NOBAL*NOBAR
      ISTIW1=IWMAST+NOBAR*NOBAL
      LSTIW1=ILCONV(ISTIW1,NINPRE)
      ISTIW3=ISTIW1+NOBAS
      LSTIW3=ILCONV(ISTIW3,NINPRE)
      IWRKEN=ISTIW3+NOBAS
      IF (IWRKEN .GE. MWORK) GOTO 1014
      REWIND NFIL30
      READ (NFIL30) (WORK(III), III=ITMASL,ITMASL+LTMAST(1)-1)
      READ (NFIL30) (WORK(III), III=ITMASR,ITMASR+LTMAST(2)-1)
C
      CALL GENYPN ( WORK(IWFRS1) , WORK(IWFLS1) , WORK(IWFRS3) ,
     .              WORK(IWFLS3) , WORK(IDFRS1) , WORK(IDFLS1) ,
     .              WORK(IDFRS3) , WORK(IDFLS3) , WORK(IY1P1S) ,
     .              WORK(IY1P2S) , WORK(IY1P3S) , WORK(IY1P4S) ,
     .              WORK(IY3P1S) , WORK(IY3P2S) , WORK(IY3P3S) ,
     .              WORK(IY3P4S) , WORK(IY1P0M) , WORK(IY1P1M) ,
     .              WORK(IY1P2M) , WORK(IY1P3M) , WORK(IY1P4M) ,
     .              WORK(IY3P0M) , WORK(IY3P1M) , WORK(IY3P2M) ,
     .              WORK(IY3P3M) , WORK(IY3P4M) , WORK(IPY1P0) ,
     .              WORK(IPY1P1) , WORK(IPY1P2) , WORK(IPY1P3) ,
     .              WORK(IPY3P0) , WORK(IPY3P1) , WORK(IPY3P2) ,
     .              WORK(IPY3P3) , WORK(IVMAST) ,
     .              WORK(ITMASL) , WORK(ITMASR) , WORK(ISMAST) ,
     .              WORK(IWMAST) , IWRK(LSTIW1) , IWRK(LSTIW3) ,
     .              NINTOP , NOBAL , NOBAR , MBASI1 , MBASI2 ,
     .              MBALP1 , MBARP1 , NOVIBS )
      CALL PRTIME( 'GENYPN' , OLDTIM , OLDVEC )
C
      LUINTS=IWMAST-ISMAST
      WRITE (NFIL10) (WORK(III), III=ISMAST,IWMAST-1)
C
C*********************************************************************
C
C  START LOOP OVER J = 0 THROUGH JMAX
C
C*********************************************************************
C
      IF (SYMM) THEN
         LZERAS(1)=LZERSY(1)
         LONEAS(1)=LONESY(1)
         LTWOAS(1)=LTWOSY(1)
         LZERAS(2)=LZERSY(2)
         LONEAS(2)=LONESY(2)
         LTWOAS(2)=LTWOSY(2)
      ENDIF
C
C*********************************************************************
C
C  DETERMINE RECORD LENGTH NECESSARY TO SAVE K-BLOCKS
C
C*********************************************************************
C
C
      IF (.NOT. SYMM) THEN
            NBLCK1=(IQUANT(2,ISO)+1)*IQUANT(4,ISO)
            NBLCK2=(IQUANT(6,ISO)+1)*IQUANT(8,ISO)
      ELSE
            NBLCK1=MAX0(IQUANT(4,ISO),IQUANT(5,ISO))
     1            *(IQUANT(2,ISO)+1)
            NBLCK2=MAX0(IQUANT(8,ISO),IQUANT(9,ISO))
     1            *(IQUANT(6,ISO)+1)
      ENDIF
C
C     LREC34: maximale Groesse eines a- oder b-Blockes
C
      LREC34=MAX0(NBLCK1,NBLCK2)
C
C     LREC35: maximale Blockgroesse fuer beide pseudoelektronischen
C     Funktionen (= eines zusammengesetzten Blockes) .
C     Datei fuer Eigenvektoren der Kontraktion wird geoeffnet.
C
      LREC35=NBLCK1+NBLCK2
      NUNIT=NFIL35
      OPEN(UNIT=NFIL35,STATUS='SCRATCH',ACCESS='DIRECT',
     1     RECL=LREC35*NRECUN,ERR=1009)
      MXFCCO=0
      NRFCCO=1
C
C     MAXPSY: Zahl der irreduziblen Darstellungen gleicher Paritaet
C
      MAXPSY=2
      IF (.NOT.SYMM) MAXPSY=1
C
C     Beschreiber der Eigenfunktionen
C
      IASSGN=IWRKST
C
C     Zahl der uebernommenen Eigenfunktionen
C
      INFCTS=IASSGN+JMAXP1*MAXPSY*2*LREC35*5
      KNFCTS=2*INFCTS-1
      IF (INFCTS.GE.MWORK) GOTO 1017
C
C     Positionen der Eigenvektoren
C
      IPSFCT=INFCTS+JMAXP1*4
C
      IONERC=IPSFCT+JMAXP1*4
      IF (IONERC.GE.MWORK) GOTO 1017
C
      IHMAST=IONERC+LREC34*LREC34
      IF (IHMAST .GE. MWORK) GOTO 1017
      ISMAST=IHMAST+LREC35*LREC35
      IF (ISMAST .GE. MWORK) GOTO 1017
C
C     Berechnung der kontrahierten Funktionen
C
C     Beginn der Schleifen
C
      DO 7010 KP1=1,JMAXP1
      DO 7008 ITAUP1=1,2
      DO 7009 IFLAG=1,MAXPSY
C
C     Matrix loeschen
C
      DO 7001 II=0,(LREC35*LREC35-1)
         WORK(IHMAST+II)=0.00
7001  CONTINUE
C
      JP1=KP1
      J=JP1-1
      ITAU=ITAUP1-1
C
      IF (MOD(J,2) .EQ. 0) THEN
            IROTE=J/2+1
            IROTO=J/2
      ELSE
            IROTE=(J+1)/2
            IROTO=IROTE
      ENDIF
      IF (MOD(J+ITAU,2) .NE. 0) THEN
               IROTE1=IROTE-1
               IROTE2=IROTE
      ELSE
               IROTE1=IROTE
               IROTE2=IROTE-1
      ENDIF
      IROTO1=IROTO
      IROTO2=IROTO
C
      INDSYM=2*(IFLAG-1)+ITAUP1
C
      IF (.NOT. SYMM) THEN
            NBLCK1=(IROTE1+IROTO1)*V2MP11*IQUANT(4,ISO)
            NBLCK2=(IROTE2+IROTO2)*V2MP12*IQUANT(8,ISO)
	    NDIMST=LSTYPA(1)
            NDIMRO=LSTYPA(1)
            NDIMCO=LSTYPA(1)
            ISTOFF=0
            JSTOFF=0
      ELSE
	    IF (IFLAG .EQ. 1) THEN
               NBLCK1=(IROTE1*IQUANT(4,ISO)+IROTO1*IQUANT(5,ISO))
     1               *V2MP11
               NBLCK2=(IROTE2*IQUANT(8,ISO)+IROTO2*IQUANT(9,ISO))
     1               *V2MP12
	    ELSE
               NBLCK1=(IROTE1*IQUANT(5,ISO)+IROTO1*IQUANT(4,ISO))
     1               *V2MP11
	       NBLCK2=(IROTE2*IQUANT(9,ISO)+IROTO2*IQUANT(8,ISO))
     1               *V2MP12
	    ENDIF
      ENDIF
      IRODIM=NBLCK1+NBLCK2
C
      IBOOK(JP1,ITAUP1,IFLAG,1)=IRODIM
      IBOOK(JP1,ITAUP1,IFLAG,2)=NBLCK1
      IBOOK(JP1,ITAUP1,IFLAG,3)=IROTE1
      IBOOK(JP1,ITAUP1,IFLAG,4)=IROTO1
      IBOOK(JP1,ITAUP1,IFLAG,5)=IROTE2
      IBOOK(JP1,ITAUP1,IFLAG,6)=IROTO2
C
C     Berechnung der aktuellen Blockgroessen (g.o.)
C
      NABLK=IQUANT(2,ISO)+1
      NBBLK=IQUANT(6,ISO)+1
      IF (SYMM) THEN
        IF (MOD(MOD(KP1-1,2)+(IFLAG-1),2).EQ.0) THEN
          NABLK=NABLK*IQUANT(4,ISO)
          NBBLK=NBBLK*IQUANT(8,ISO)
        ELSE
          NABLK=NABLK*IQUANT(5,ISO)
          NBBLK=NBBLK*IQUANT(9,ISO)
        ENDIF
      ELSE
        NABLK=NABLK*IQUANT(4,ISO)
        NBBLK=NBBLK*IQUANT(8,ISO)
      ENDIF
      IF (KP1.EQ.1) THEN
        IF (MOD(J+ITAU,2).EQ.0) THEN 
          NBBLK=0
        ELSE 
          NABLK=0
        ENDIF
      ENDIF
C*********************************************************************
C
C  CONSTRUCT DIAGONAL BLOCK FOR SURFACES 1 AND 2
C
C*********************************************************************
C
C
      REWIND NFIL9
      REWIND NFIL10
      REWIND NFIL33
C
      DO 324 ISURF=1,2
C
      NSURF=ISURF
      IF (ISURF.EQ.1) THEN
        NBLK=NABLK
      ELSE
        NBLK=NBBLK
      ENDIF
      IF (ISURF .EQ. 1) THEN
         V2MAX=V2MP11-1
         V2MXP1=V2MAX+1
         LSTYPA(1)=IQUANT(4,ISO)
         LSTYPA(2)=IQUANT(5,ISO)
         KSTYPA(1)=KSTY11
         KSTYPA(2)=KSTY21
         NOBAS=NOBAL
         ITAU=ITAUP1-1
      ELSE
         V2MAX=V2MP12-1
         V2MXP1=V2MAX+1
         LSTYPA(1)=IQUANT(8,ISO)
         LSTYPA(2)=IQUANT(9,ISO)
         KSTYPA(1)=KSTY12
         KSTYPA(2)=KSTY22
         NOBAS=NOBAR
         ITAU=2-ITAUP1
      ENDIF
C
      ISYMS0=ISMAST+LISTRT(ISURF)
      IF (ISYMS0 .GE. MWORK) GOTO 1010
      READ (NFIL10) (WORK(III), III=ISMAST,ISYMS0-1)
      IASYS0=ISYMS0+LZERSY(ISURF)
      IF (IASYS0 .GE. MWORK) GOTO 1010
      READ (NFIL9) (WORK(III), III=ISYMS0,IASYS0-1)
      IF (.NOT.SYMM) THEN
          ISYMS2=IASYS0+LZERAS(ISURF)
          IF (ISYMS2 .GE. MWORK) GOTO 1010
          READ (NFIL9) (WORK(III), III=IASYS0,ISYMS2-1)
      ELSE
          ISYMS2=IASYS0
      ENDIF
      IASYS2=ISYMS2+LTWOSY(ISURF)
      IF (IASYS2 .GE. MWORK) GOTO 1010
      READ (NFIL9) (WORK(III), III=ISYMS2,IASYS2-1)
      IF (.NOT.SYMM) THEN
          IENRGS=IASYS2+LTWOAS(ISURF)
          IF (IENRGS .GE. MWORK) GOTO 1010
          READ (NFIL9) (WORK(III), III=IASYS2,IENRGS-1)
      ELSE
          IENRGS=IASYS2
      ENDIF
      IESTST=IENRGS+LEBEND(ISURF)
      IF (IESTST .GE. MWORK) GOTO 1010
      READ (NFIL33) (WORK(III), III=IENRGS,IESTST-1)
      IWRKEN=IESTST+LESTRT(ISURF)
      IF (IWRKEN .GE. MWORK) GOTO 1010
      READ (NFIL33) (WORK(III), III=IESTST,IWRKEN-1)
C
C*********************************************************************
C
C  ADJUST ENERGY ZERO FOR SURFACE 2 RELATIVE TO SURFACE 1
C
C*********************************************************************
C
      IF (ISURF .EQ. 2) THEN
      DO 372 III=IESTST,IWRKEN-1
372   WORK(III)=WORK(III)+VMINS2-VMINS1
      ENDIF
C
C*********************************************************************
C
C  SET UP THE EVEN K CORNER OF THE HAMILTONIAN MATRIX BLOCK
C
C*********************************************************************
C
      IF (ISURF .EQ. 1) THEN
         IROTE=IROTE1
         IROTO=IROTO1
         NOFSET=0
      ELSE
         IROTE=IROTE2
         IROTO=IROTO2
         NOFSET=NABLK
      ENDIF
      IEO=0
      IF (IFLAG .EQ. 1) THEN
         NDIMST=LSTYPA(1)
         ISTOFF=0
      ELSE
         NDIMST=LSTYPA(2)
         ISTOFF=KSTYPA(1)
      ENDIF
      iF (ISURF .EQ. 1) THEN
         IBOOK(JP1,ITAUP1,IFLAG,7)=NDIMST
         IBOOK(JP1,ITAUP1,IFLAG,11)=V2MXP1
      ELSE
         IBOOK(JP1,ITAUP1,IFLAG,9)=NDIMST
         IBOOK(JP1,ITAUP1,IFLAG,12)=V2MXP1
      ENDIF
      IF (NBLK .NE. 0) THEN
      CALL DELK0 ( WORK(IENRGS) , WORK(IESTST) , WORK(ISYMS0) ,
     1             WORK(IASYS0) , WORK(ISMAST) , WORK(IHMAST) ,
     2             NOFSET , NDIMST , JP1,KP1    , IRODIM , IEO    ,
     3             ITAU   , LZERSY(ISURF) , LZERAS(ISURF) , ISTOFF ,
     1             NOBAS  , NOPERS , NZERSY , NZERAS , LREC35 )
      ENDIF
C*********************************************************************
C
C  SET UP THE ODD K CORNER OF THE HAMILTONIAN MATRIX BLOCK
C
C*********************************************************************
321   CONTINUE
      IEO=1
      IF (ISURF .EQ. 1) THEN
         NOFSET=0
      ELSE
         NOFSET=NABLK
      ENDIF
      IF (IFLAG .EQ. 1 .AND. SYMM) THEN
            NDIMST=LSTYPA(2)
            ISTOFF=KSTYPA(1)
      ELSE
            NDIMST=LSTYPA(1)
            ISTOFF=0
      ENDIF
      IF (ISURF .EQ. 1) THEN
         IBOOK(JP1,ITAUP1,IFLAG,8)=NDIMST
      ELSE
         IBOOK(JP1,ITAUP1,IFLAG,10)=NDIMST
      ENDIF
      IF (NBLK .NE. 0) THEN
      CALL DELK0 ( WORK(IENRGS) , WORK(IESTST) , WORK(ISYMS0) ,
     1             WORK(IASYS0) , WORK(ISMAST) , WORK(IHMAST) ,
     2             NOFSET , NDIMST , JP1,KP1    , IRODIM , IEO    ,
     3             ITAU   , LZERSY(ISURF) , LZERAS(ISURF) , ISTOFF ,
     1             NOBAS  , NOPERS , NZERSY , NZERAS , LREC35)
      ENDIF
340   CONTINUE
      ISYMS1=ISYMS0
      IASYS1=ISYMS1+LONESY(ISURF)
      IF (IASYS1 .GE. MWORK) GOTO 1010
      READ (NFIL9) (WORK(III), III=ISYMS1,IASYS1-1)
      IF (.NOT.SYMM) THEN
          IWRKEN=IASYS1+LONEAS(ISURF)
          IF (IWRKEN .GE. MWORK) GOTO 1010
          READ (NFIL9) (WORK(III), III=IASYS1,IWRKEN-1)
      ELSE
          IWRKEN=IASYS1
      ENDIF
      IF (IWRKEN .GE. MWORK) GOTO 1010
324   CONTINUE
C*********************************************************************
C
C  CALCULATE THE RENNER-TELLER INTERACTION TERMS
C
C*********************************************************************
C
C  K EVEN -- K EVEN BLOCK
C
      ITAU = ITAUP1 - 1
C
      ISMART=ISMAST
      ISYFC0=ISMART+LUINTS
      IF (ISYFC0 .GE. MWORK) GOTO 1010
      READ (NFIL10) (WORK(III), III=ISMART,ISYFC0-1)
      IASFC0=ISYFC0+IZERSY
      IF (ITRNSY .GE. MWORK) GOTO 1010
      READ (NFIL9) (WORK(III), III=ISYFC0,IASFC0-1)
C
      IF (.NOT.SYMM) THEN
	    ISYFC2=IASFC0+IZERAS
            IF (ISYFC2.GE. MWORK) GOTO 1010
            READ (NFIL9) (WORK(III), III=IASFC0,ISYFC2-1)
      ELSE
	    ISYFC2=IASFC0
      ENDIF
C
      IASFC2=ISYFC2+ITWOSY
      IF (IASFC2 .GE. MWORK) GOTO 1010
      READ (NFIL9) (WORK(III), III=ISYFC2,IASFC2-1)
C
      IF (.NOT.SYMM) THEN
	    IWRKEN=IASFC2+ITWOAS
            IF (IWRKEN .GE. MWORK) GOTO 1010
            READ (NFIL9) (WORK(III), III=IASFC2,IWRKEN-1)
      ELSE
	    IWRKEN=IASFC2
            IF (IWRKEN .GE. MWORK) GOTO 1010
      ENDIF
C
      IF (NABLK*NBBLK.GT.0) THEN
        NOFSRO=0
        NOFSCO=NABLK
        IF (IFLAG .EQ. 1) THEN
          NDIMRO=IQUANT(4,ISO)
          NDIMCO=IQUANT(8,ISO)
          ISTOFF=0
          JSTOFF=0
        ELSE
          NDIMRO=IQUANT(5,ISO)
          NDIMCO=IQUANT(9,ISO)
          ISTOFF=KSTY11
          JSTOFF=KSTY12
        ENDIF
        IEO=0
C
      CALL RNPM0 ( WORK(ISYFC0) , WORK(IASFC0) ,
     1             WORK(ISMART) ,
     1             WORK(IHMAST) , NOFSRO , NOFSCO ,
     1             NDIMRO , NDIMCO , V2MP11 , V2MP12 , JP1 ,KP1,
     2             IRODIM , IEO    , ITAU   , IZERSY , IZERAS ,
     3             ISTOFF , JSTOFF , NOBAL  , NOBAR  , NINTOP ,
     4             NINSYM,NINASY,LREC35)
C
C  K ODD -- K ODD BLOCK
C
      IF (IFLAG .EQ. 1) THEN
      NOFSRO=0
      NOFSCO=NABLK
      IF (.NOT.SYMM) THEN
           NDIMRO=IQUANT(4,ISO)
           NDIMCO=IQUANT(8,ISO)
           ISTOFF=0
           JSTOFF=0
      ELSE
           NDIMRO=IQUANT(5,ISO)
           NDIMCO=IQUANT(9,ISO)
           ISTOFF=KSTY11
           JSTOFF=KSTY12
      ENDIF
      ELSE
           NOFSRO=0
	   NOFSCO=NABLK
	   NDIMRO=IQUANT(4,ISO)
	   NDIMCO=IQUANT(8,ISO)
	   ISTOFF=0
	   JSTOFF=0
      ENDIF
      IEO=1
C
      CALL RNPM0 ( WORK(ISYFC0) , WORK(IASFC0) ,
     1             WORK(ISMART) ,
     1             WORK(IHMAST) , NOFSRO , NOFSCO ,
     1             NDIMRO , NDIMCO , V2MP11 , V2MP12 , JP1 ,KP1,
     2             IRODIM , IEO    , ITAU   , IZERSY , IZERAS ,
     3             ISTOFF , JSTOFF , NOBAL  , NOBAR  , NINTOP ,
     4             NINSYM , NINASY, LREC35 )
      ENDIF
C
C     Groesse der Matrix ermitteln
C
      IF (MOD(KP1,2).EQ.1) THEN
        IABLK=IBOOK(KP1,ITAUP1,IFLAG,7)*IBOOK(KP1,ITAUP1,IFLAG,11)
        IBBLK=IBOOK(KP1,ITAUP1,IFLAG,9)*IBOOK(KP1,ITAUP1,IFLAG,12)
      ELSE
        IABLK=IBOOK(KP1,ITAUP1,IFLAG,8)*IBOOK(KP1,ITAUP1,IFLAG,11)
        IBBLK=IBOOK(KP1,ITAUP1,IFLAG,10)*IBOOK(KP1,ITAUP1,IFLAG,12)
      ENDIF
      IF (KP1.EQ.1) THEN
        IF (ITAUP1.EQ.1) THEN
          NABLK=IABLK
          NBBLK=0
        ELSE
          NABLK=0
          NBBLK=IBBLK
        ENDIF
      ELSE
        NABLK=IABLK
        NBBLK=IBBLK
      ENDIF
      IDTOT=NABLK+NBBLK
C     Matrizen eingrenzen und Platz fuer Workspace schaffen
      IWRKEN=ISMAST+LREC35*LREC35
      IF (IWRKEN.GT.MWORK) GOTO 1017
      IWORK=IWRKEN+8*IDTOT
      IF (IWORK.GT.MWORK) GOTO 1017
      IWRKEI=IWORK+5*IDTOT
      IF (IWRKEI.GT.MWORK) GOTO 1017
      IWRKFA=IWRKEI+IDTOT
      IF (IWRKFA.GT.MWORK) GOTO 1017
      IWRKPA=IWRKFA+IDTOT
      IF (IWRKPA.GT.MWORK) GOTO 1017
      IWRKPB=IWRKPA+IDTOT
      IF (IWRKPB.GT.MWORK) GOTO 1017
      IWRKES=IWRKPB+IDTOT
      IF (IWRKES.GT.MWORK) GOTO 1017
      IWRKAB=IWRKES+2*IDTOT
      IF (IWRKAB.GT.MWORK) GOTO 1017
      IWRKDA=IWRKAB+2*IDTOT
      IF (IWRKDA.GT.MWORK) GOTO 1017
      IWRKUP=IWRKDA+6*IDTOT
      IF (IWRKUP.GT.MWORK) GOTO 1017
C
C     Kontrahierung aufrufen
C
      Call CONTRA(WORK(IHMAST),
     1            WORK(IWRKEI),WORK(IWRKFA),LREC34,
     1            LREC35,IDTOT,NABLK,NBBLK,
     1            WORK(IWRKEN),WORK(IWORK),
     2            WORK(IASSGN),MAXPSY,
     2            WORK(INFCTS),
     2            IFLAG,ITAUP1,KP1,
     3            COENMX,ISOMAX,ISO,
     5            WORK(IWRKPA),WORK(IWRKPB),WORK(IWRKES),WORK(IWRKAB),
     6            WORK(IWRKDA),
     6            IBOOK,MXJVP1,NMINP1,
     7            WORK(IPSFCT),MXFCCO,NRFCCO,
     8            SYMM)
C
7009  CONTINUE
7008  CONTINUE
7010  CONTINUE
C
C     BERECHNUNG DER DIPOLMOMENTMATRIZEN (g.o.)
C
      IF (IACT.GE.1) THEN
      DPFL=.TRUE.
      IWRKEN=ISMAST
      REWIND NFIL10
      ISTR11=IWRKEN
      ISTR22=ISTR11+LISTRT(1)
      ISTR12=ISTR22+LISTRT(2)
      IWRKEN=ISTR12+LUINTS
      IF (IWRKEN.GT.MWORK) GOTO 1017
C
      CALL READFI(NFIL10,WORK(ISTR11),LISTRT(1))
      CALL READFI(NFIL10,WORK(ISTR22),LISTRT(2))
      CALL READFI(NFIL10,WORK(ISTR12),LUINTS   )
C
      REWIND NFIL17
      REWIND NFIL20
      REWIND NFIL19 
      REWIND NFIL21
      REWIND NFAB0
      REWIND NFAB1      
C
      IB110S=IWRKEN
      IB110A=IB110S+1*NDIPSY*JMAXP1*V2MP11*V2MP11
      IB220S=IB110A+1*NDIPAS*JMAXP1*V2MP11*V2MP11
      IB220A=IB220S+1*NDIPSY*JMAXP1*V2MP12*V2MP12
      IB111S=IB220A+1*NDIPAS*JMAXP1*V2MP12*V2MP12
      IB111A=IB111S+2*NDIPSY*JMAXP1*V2MP11*V2MP11
      IB221S=IB111A+2*NDIPAS*JMAXP1*V2MP11*V2MP11
      IB221A=IB221S+2*NDIPSY*JMAXP1*V2MP12*V2MP12
      IB120S=IB221A+2*NDIPAS*JMAXP1*V2MP12*V2MP12
      IB120A=IB120S+1*NDIPSY*JMAXP1*V2MP11*V2MP12
      IB121S=IB120A+1*NDIPAS*JMAXP1*V2MP11*V2MP12
      IB121A=IB121S+2*NDIPSY*JMAXP1*V2MP11*V2MP12
      IWRKEN=IB121A+2*NDIPAS*JMAXP1*V2MP11*V2MP12
      IF (IWRKEN.GT.MWORK) GOTO 1017
C
      CALL READFI(NFIL17,WORK(IB110S),1*NDIPSY*JMAXP1*V2MP11*V2MP11)
      IF (.NOT.SYMM)
     1   CALL READFI(NFIL17,WORK(IB110A),
     2               1*NDIPAS*JMAXP1*V2MP11*V2MP11)
      CALL READFI(NFIL20,WORK(IB111S),2*NDIPSY*(JMAXP1-1)*V2MP11*V2MP11)
      IF (.NOT.SYMM)
     1   CALL READFI(NFIL20,WORK(IB111A),
     2               2*NDIPAS*(JMAXP1-1)*V2MP11*V2MP11)
      CALL READFI(NFIL19,WORK(IB220S),1*NDIPSY*JMAXP1*V2MP12*V2MP12)
      IF (.NOT.SYMM)
     1   CALL READFI(NFIL19,WORK(IB220A),
     2               1*NDIPAS*JMAXP1*V2MP12*V2MP12)
      CALL READFI(NFIL21,WORK(IB221S),2*NDIPSY*(JMAXP1-1)*V2MP12*V2MP12)
      IF (.NOT.SYMM) 
     1   CALL READFI(NFIL21,WORK(IB221A),
     2               2*NDIPAS*(JMAXP1-1)*V2MP12*V2MP12)
      CALL READFI(NFAB0 ,WORK(IB120S),1*NDIPSY*JMAXP1*V2MP11*V2MP12)
      IF (.NOT.SYMM)
     1   CALL READFI(NFAB0 ,WORK(IB120A),
     2               1*NDIPAS*JMAXP1*V2MP11*V2MP12)
      CALL READFI(NFAB1 ,WORK(IB121S),2*NDIPSY*(JMAXP1-1)*V2MP11*V2MP12)
      IF (.NOT.SYMM)  
     1   CALL READFI(NFAB1 ,WORK(IB121A),
     2               2*NDIPAS*(JMAXP1-1)*V2MP11*V2MP12)
C
      IWRKU=IWRKEN
C
      NUNIT=NFIL41
      OPEN(UNIT=NFIL41,STATUS='SCRATCH',ACCESS='DIRECT',
     1     RECL=MXFCCO*MXFCCO*NRECUN,ERR=1009)
      NUNIT=NFIL42
      OPEN(UNIT=NFIL42,STATUS='SCRATCH',ACCESS='DIRECT',
     1     RECL=MXFCCO*MXFCCO*NRECUN,ERR=1009)
      NUNIT=NFIL43
      OPEN(UNIT=NFIL43,STATUS='SCRATCH',ACCESS='DIRECT',
     1     RECL=MXFCCO*MXFCCO*NRECUN,ERR=1009)
C
      ITAUR=0
      ITAUC=1
      DO 9112 IFLAGR=1,MAXPSY
      IF (.NOT.SYMM) THEN
         IFLAGC=IFLAGR
      ELSE
         IFLAGC=3-IFLAGR  
      ENDIF 
      INDSYR=2*(IFLAGR-1)+ITAUR+1
      INDSYC=2*(IFLAGC-1)+ITAUC+1
      DO 9113 KDEL=-1,1,1
      DO 9114 K1P1=1,JMAXP1
      K1VAL=K1P1-1
      K2P1=K1P1+KDEL
      K2VAL=K2P1-1
      IF ((K2P1.LT.1).OR.(K2P1.GT.JMAXP1)) GOTO 9114
C
      IF (K1P1.EQ.1) THEN
         NPARRM=1
      ELSE
         NPARRM=0
      ENDIF
      IF (K2P1.EQ.1) THEN
         NPARCM=1
      ELSE
         NPARCM=0
      ENDIF
      DO 9115 NPARR=0,NPARRM
      DO 9116 NPARC=0,NPARCM
C 
      NABLKR=V2MP11
      NABLKC=V2MP11
      NBBLKR=V2MP12
      NBBLKC=V2MP12
      IF (SYMM) THEN
        IF (((MOD(K1P1,2).EQ.1).AND.(IFLAGR.EQ.1)).OR.
     1     ((MOD(K1P1,2).EQ.0).AND.(IFLAGR.EQ.2))) THEN
          NABLKR=NABLKR*IQUANT(4,ISO)
          NBBLKR=NBBLKR*IQUANT(8,ISO)
        ELSE
          NABLKR=NABLKR*IQUANT(5,ISO)
          NBBLKR=NBBLKR*IQUANT(9,ISO)
        ENDIF
      ELSE
          NABLKR=NABLKR*IQUANT(4,ISO)
          NBBLKR=NBBLKR*IQUANT(8,ISO)
      ENDIF
      IF (SYMM) THEN
        IF (((MOD(K2P1,2).EQ.1).AND.(IFLAGC.EQ.1)).OR.
     1     ((MOD(K2P1,2).EQ.0).AND.(IFLAGC.EQ.2))) THEN
          NABLKC=NABLKC*IQUANT(4,ISO)
          NBBLKC=NBBLKC*IQUANT(8,ISO)
        ELSE
          NABLKC=NABLKC*IQUANT(5,ISO)
          NBBLKC=NBBLKC*IQUANT(9,ISO)
        ENDIF
      ELSE
          NABLKC=NABLKC*IQUANT(4,ISO)
          NBBLKC=NBBLKC*IQUANT(8,ISO)
      ENDIF
      IF (K1P1.EQ.1) THEN
        IF (MOD(NPARR+ITAUR,2).EQ.0) THEN 
          NBBLKR=0
        ELSE
          NABLKR=0
        ENDIF
      ENDIF
      IF (K2P1.EQ.1) THEN
        IF (MOD(NPARC+ITAUC,2).EQ.0) THEN 
          NBBLKC=0
        ELSE
          NABLKC=0
        ENDIF
      ENDIF
      DIMR=NABLKR+NBBLKR
      DIMC=NABLKC+NBBLKC
      IF (NABLKR.GT.NBBLKR) THEN
        NBLKR=NABLKR
      ELSE
        NBLKR=NBBLKR
      ENDIF
      IF (NABLKC.GT.NBLKR) NBLKR=NABLKC
      IF (NABLKC.GT.NBBLKC) THEN
        NBLKC=NABLKC
      ELSE
        NBLKC=NBBLKC
      ENDIF
      IF (NBBLKR.GT.NBLKC) NBLKC=NBBLKR
C
      IWRKEN=IWRKU
      IBLTOT=IWRKEN
      IBLKAB=IBLTOT+LREC35*LREC35
      IWRKBT=IBLKAB+NBLKC*NBLKR
      IWRKWM=IWRKBT+LREC35*MXFCCO
      IWRKEN=IWRKWM+LREC35*MXFCCO
      IF (IWRKEN.GT.MWORK) GOTO 1017    
C
      IF (((MOD(K1P1,2).EQ.1).AND.(IFLAGR.EQ.1)).OR.
     1    ((MOD(K1P1,2).EQ.0).AND.(IFLAGR.EQ.2))) THEN
          NABSYM=1
      ELSE
          NABSYM=2
      ENDIF
      IF (.NOT.SYMM) NABSYM=1
      NABSR=NABSYM
      IF (KDEL.EQ.0) THEN
        NABSC=3-NABSR
      ELSE
        NABSC=NABSR
      ENDIF
      IF (.NOT.SYMM) NABSC=1
C 
      IF ((KDEL.EQ.(-1)).OR.(KDEL.EQ.(+1))) THEN
C
         IF ((NABLKR*NABLKC).GT.0) THEN
           PSEUR='a'
           PSEUC='a'
           INTSYM=2*NDIPSY*(JMAXP1-1)*V2MP11*V2MP11
           INTASY=2*NDIPAS*(JMAXP1-1)*V2MP11*V2MP11
           IF (NABSYM.EQ.1) THEN
             ISTOFF=0
             JSTOFF=0
           ELSE
             ISTOFF=KSTY11
             JSTOFF=KSTY11
           ENDIF
           IF (NABSYM.EQ.1) THEN
              NDIMRO=IQUANT(4,ISO)
              NDIMCO=IQUANT(4,ISO)
           ELSE
              NDIMRO=IQUANT(5,ISO)
              NDIMCO=IQUANT(5,ISO)
           ENDIF
           AABB=.TRUE.
           CALL MUEDK1(WORK(IB111S),INTSYM,
     1                 WORK(IB111A),INTASY,
     1                 WORK(ISTR11),NOBAL,NOBAL,NOPERS,
     1                 ISTOFF,JSTOFF,
     1                 WORK(IBLKAB),NABLKR,NABLKC,
     1                 NDIMRO,NDIMCO,
     1                 V2MP11,V2MP11,
     1                 JMAXP1,K1P1,KDEL,
     1                 NDIPSY,NDIPAS,
     1                 SYMM,
     2                 AABB,
     1                 NABSR,NABSC,
     1                 PSEUR,PSEUC)
           CALL BLKCOP(WORK(IBLTOT),LREC35,LREC35,
     1                 WORK(IBLKAB),NABLKR,NABLKC,
     1                 0,0,.FALSE.)
         ENDIF
C
         IF ((NBBLKR*NBBLKC).GT.0) THEN
           PSEUR='b'
           PSEUC='b'
           INTSYM=2*NDIPSY*(JMAXP1-1)*V2MP12*V2MP12
           INTASY=2*NDIPAS*(JMAXP1-1)*V2MP12*V2MP12
           IF (NABSYM.EQ.1) THEN
              ISTOFF=0
              JSTOFF=0
              NDIMRO=IQUANT(8,ISO)
              NDIMCO=IQUANT(8,ISO)
           ELSE
              ISTOFF=KSTY12
              JSTOFF=KSTY12
              NDIMRO=IQUANT(9,ISO)
              NDIMCO=IQUANT(9,ISO)
           ENDIF
           AABB=.TRUE.
           CALL MUEDK1(WORK(IB221S),INTSYM,
     1                 WORK(IB221A),INTASY,
     1                 WORK(ISTR22),NOBAR,NOBAR,NOPERS,
     1                 ISTOFF,JSTOFF,
     1                 WORK(IBLKAB),NBBLKR,NBBLKC,
     1                 NDIMRO,NDIMCO,
     1                 V2MP12,V2MP12,
     1                 JMAXP1,K1P1,KDEL,
     1                 NDIPSY,NDIPAS,
     1                 SYMM,
     2                 AABB,
     3                 NABSR,NABSC,
     1                 PSEUR,PSEUC)
           CALL BLKCOP(WORK(IBLTOT),LREC35,LREC35,
     1                 WORK(IBLKAB),NBBLKR,NBBLKC,
     1                 NABLKR,NABLKC,.FALSE.)
         ENDIF
C
         IF ((NABLKR*NBBLKC).GT.0) THEN
           PSEUR='a'
           PSEUC='b' 
           INTSYM=2*NDIPSY*(JMAXP1-1)*V2MP11*V2MP12
           INTASY=2*NDIPAS*(JMAXP1-1)*V2MP11*V2MP12
           IF (NABSYM.EQ.1) THEN
             ISTOFF=0
             JSTOFF=0 
             NDIMRO=IQUANT(4,ISO)
             NDIMCO=IQUANT(8,ISO)      
           ELSE
             ISTOFF=KSTY11
             JSTOFF=KSTY12
             NDIMRO=IQUANT(5,ISO)
             NDIMCO=IQUANT(9,ISO)
           ENDIF
           AABB=.FALSE.
           CALL MUEDK1(WORK(IB121S),INTSYM,
     1                 WORK(IB121A),INTASY,
     1                 WORK(ISTR12),NOBAL,NOBAR,NINTOP,
     1                 ISTOFF,JSTOFF,
     1                 WORK(IBLKAB),NABLKR,NBBLKC,
     1                 NDIMRO,NDIMCO,
     1                 V2MP11,V2MP12,
     1                 JMAXP1,K1P1,KDEL,
     1                 NDIPSY,NDIPAS,
     1                 SYMM,
     2                 AABB,
     1                 NABSR,NABSC,
     1                 PSEUR,PSEUC)
           CALL BLKCOP(WORK(IBLTOT),LREC35,LREC35,
     1                 WORK(IBLKAB),NABLKR,NBBLKC,
     1                 0,NABLKC,.FALSE.)
         ENDIF
C
         IF ((NBBLKR*NABLKC).GT.0) THEN
           PSEUR='b'
           PSEUC='a'
           INTSYM=2*NDIPSY*(JMAXP1-1)*V2MP11*V2MP12
           INTASY=2*NDIPAS*(JMAXP1-1)*V2MP11*V2MP12
           IF (NABSYM.EQ.1) THEN
              ISTOFF=0
              JSTOFF=0
              NDIMRO=IQUANT(4,ISO)
              NDIMCO=IQUANT(8,ISO)             
           ELSE
              ISTOFF=KSTY11
              JSTOFF=KSTY12
              NDIMRO=IQUANT(5,ISO)
              NDIMCO=IQUANT(9,ISO)
           ENDIF
           AABB=.FALSE. 
           CALL MUEDK1(WORK(IB121S),INTSYM,
     1                 WORK(IB121A),INTASY,
     1                 WORK(ISTR12),NOBAL,NOBAR,NINTOP,
     1                 ISTOFF,JSTOFF,
     1                 WORK(IBLKAB),NABLKC,NBBLKR,
     1                 NDIMRO,NDIMCO,
     1                 V2MP11,V2MP12,
     1                 JMAXP1,K2P1,(-1)*KDEL,
     1                 NDIPSY,NDIPAS,
     1                 SYMM,
     1                 AABB,
     2                 NABSR,NABSC,
     1                 PSEUR,PSEUC)
           CALL BLKCOP(WORK(IBLTOT),LREC35,LREC35,
     1                 WORK(IBLKAB),NABLKC,NBBLKR,
     1                 NABLKR,0,.TRUE.)
         ENDIF
C
      ENDIF
C
      IF (KDEL.EQ.0) THEN
C
         IF ((NABLKR*NABLKC).GT.0) THEN
           PSEUR='a'
           PSEUC='a'
           INTSYM=1*NDIPSY*JMAXP1*V2MP11*V2MP11
           INTASY=1*NDIPAS*JMAXP1*V2MP11*V2MP11
           IF (.NOT.SYMM) THEN
              ISTOFF=0
              JSTOFF=0
              NDIMRO=IQUANT(4,ISO)
              NDIMCO=IQUANT(4,ISO)
           ELSE
              IF (NABSYM.EQ.1) THEN
                ISTOFF=0
                JSTOFF=KSTY11
                NDIMRO=IQUANT(4,ISO)
                NDIMCO=IQUANT(5,ISO)             
              ELSE
                ISTOFF=KSTY11
                JSTOFF=0
                NDIMRO=IQUANT(5,ISO)
                NDIMCO=IQUANT(4,ISO)
              ENDIF
           ENDIF
           AABB=.TRUE.
           CALL MUEDK0(WORK(IB110S),INTSYM,
     1                 WORK(IB110A),INTASY,
     1                 WORK(ISTR11),NOBAL,NOBAL,NOPERS,
     1                 ISTOFF,JSTOFF,
     1                 WORK(IBLKAB),NABLKR,NABLKC,
     1                 NDIMRO,NDIMCO,
     1                 V2MP11,V2MP11,
     1                 JMAXP1,K1P1,
     1                 NDIPSY,NDIPAS,
     1                 SYMM,
     2                 AABB,
     1                 NABSR,NABSC,
     1                 PSEUR,PSEUC)
           CALL BLKCOP(WORK(IBLTOT),LREC35,LREC35,
     1                 WORK(IBLKAB),NABLKR,NABLKC,
     1                 0,0,.FALSE.) 
         ENDIF
C
         IF ((NBBLKR*NBBLKC).GT.0) THEN
           PSEUR='b'
           PSEUC='b'
           INTSYM=1*NDIPSY*JMAXP1*V2MP12*V2MP12
           INTASY=1*NDIPAS*JMAXP1*V2MP12*V2MP12
           IF (.NOT.SYMM) THEN
              ISTOFF=0
              JSTOFF=0
              NDIMRO=IQUANT(8,ISO)
              NDIMCO=IQUANT(8,ISO)
           ELSE
              IF (NABSYM.EQ.1) THEN
                ISTOFF=0
                JSTOFF=KSTY12
                NDIMRO=IQUANT(8,ISO)
                NDIMCO=IQUANT(9,ISO)
              ELSE
                ISTOFF=KSTY12
                JSTOFF=0
                NDIMRO=IQUANT(9,ISO)
                NDIMCO=IQUANT(8,ISO)
              ENDIF
           ENDIF
           AABB=.TRUE.
           CALL MUEDK0(WORK(IB220S),INTSYM,
     1                 WORK(IB220A),INTASY,
     1                 WORK(ISTR22),NOBAR,NOBAR,NOPERS,
     1                 ISTOFF,JSTOFF,
     1                 WORK(IBLKAB),NBBLKR,NBBLKC,
     1                 NDIMRO,NDIMCO,
     1                 V2MP12,V2MP12,
     1                 JMAXP1,K1P1,
     1                 NDIPSY,NDIPAS,
     1                 SYMM,
     1                 AABB,
     2                 NABSR,NABSC,
     1                 PSEUR,PSEUC)
           CALL BLKCOP(WORK(IBLTOT),LREC35,LREC35,
     1                 WORK(IBLKAB),NBBLKR,NBBLKC,
     1                 NABLKR,NABLKC,.FALSE.)
         ENDIF
C
         IF ((NABLKR*NBBLKC).GT.0) THEN
           PSEUR='a'
           PSEUC='b'
           INTSYM=1*NDIPSY*JMAXP1*V2MP11*V2MP12
           INTASY=1*NDIPAS*JMAXP1*V2MP11*V2MP12
           IF (.NOT.SYMM) THEN
             ISTOFF=0
             JSTOFF=0
             NDIMRO=IQUANT(4,ISO)
             NDIMCO=IQUANT(8,ISO)
           ELSE
             IF (NABSYM.EQ.1) THEN
               ISTOFF=0
               JSTOFF=KSTY12
               NDIMRO=IQUANT(4,ISO)
               NDIMCO=IQUANT(9,ISO)
             ELSE
               ISTOFF=KSTY11
               JSTOFF=0
               NDIMRO=IQUANT(5,ISO)
               NDIMCO=IQUANT(8,ISO)
             ENDIF
           ENDIF
           AABB=.FALSE.
           CALL MUEDK0(WORK(IB120S),INTSYM,
     1                 WORK(IB120A),INTASY,
     1                 WORK(ISTR12),NOBAL,NOBAR,NINTOP,
     1                 ISTOFF,JSTOFF,
     1                 WORK(IBLKAB),NABLKR,NBBLKC,
     1                 NDIMRO,NDIMCO,
     1                 V2MP11,V2MP12,
     1                 JMAXP1,K1P1,
     1                 NDIPSY,NDIPAS,
     1                 SYMM,
     2                 AABB,
     1                 NABSR,NABSC,
     1                 PSEUR,PSEUC)
           CALL BLKCOP(WORK(IBLTOT),LREC35,LREC35,
     1                 WORK(IBLKAB),NABLKR,NBBLKC,
     1                 0,NABLKC,.FALSE.)
         ENDIF
C
         IF ((NBBLKR*NABLKC).GT.0) THEN
           PSEUR='b'
           PSEUC='a'
           INTSYM=1*NDIPSY*JMAXP1*V2MP11*V2MP12
           INTASY=1*NDIPAS*JMAXP1*V2MP11*V2MP12
           IF (.NOT.SYMM) THEN
              ISTOFF=0
              JSTOFF=0
              NDIMRO=IQUANT(4,ISO)
              NDIMCO=IQUANT(8,ISO)
           ELSE
              IF (NABSYM.EQ.1) THEN
                 ISTOFF=KSTY11
                 JSTOFF=0
                 NDIMRO=IQUANT(5,ISO)
                 NDIMCO=IQUANT(8,ISO)               
              ELSE
                 ISTOFF=0
                 JSTOFF=KSTY12
                 NDIMRO=IQUANT(4,ISO)
                 NDIMCO=IQUANT(9,ISO)
              ENDIF
           ENDIF
           AABB=.FALSE.
           CALL MUEDK0(WORK(IB120S),INTSYM,
     1                 WORK(IB120A),INTASY,
     1                 WORK(ISTR12),NOBAL,NOBAR,NINTOP,
     1                 ISTOFF,JSTOFF,
     1                 WORK(IBLKAB),NABLKC,NBBLKR,
     1                 NDIMRO,NDIMCO,
     1                 V2MP11,V2MP12,
     1                 JMAXP1,K1P1,
     1                 NDIPSY,NDIPAS,
     1                 SYMM,
     1                 AABB,
     2                 NABSR,NABSC,
     1                 PSEUR,PSEUC)
           CALL BLKCOP(WORK(IBLTOT),LREC35,LREC35,
     1                 WORK(IBLKAB),NABLKC,NBBLKR,
     1                 NABLKR,0,.TRUE.)
         ENDIF
C
      ENDIF
C
      CALL BSSTRN(NPARR+1,K1P1,NPARC+1,K2P1,
     1            INDSYR,INDSYC,
     2            WORK(IPSFCT),WORK(INFCTS),
     2            JMAXP1,
     2            NROW,NCOL,
     2            DIMR,DIMC,
     3            WORK(IBLTOT),WORK(IWRKBT),WORK(IWRKWM),
     4            LREC35,MXFCCO,
     5            NFIL35)
C
      CALL RCNRDP(IFLAGR,
     1            NPARR,NPARC,
     2            K1P1,K2P1,KDEL,
     3            NRCNR)
C
      IF (IDPTR.EQ.1) THEN
         CALL DPTRA1(WORK(IBLTOT),LREC35,
     1               NROW,NCOL,
     2               NFIL6,
     3               INDSYR,INDSYC,
     4               K1P1,KDEL,
     5               NPARR,NPARC)
      ENDIF
C
      IF (IDPTR.EQ.2) THEN
         CALL DPTRA2(WORK(IBLTOT),LREC35,
     1               WORK(IASSGN),JMAXP1,MAXPSY,
     1               NROW,NCOL,
     2               NFIL6,
     3               INDSYR,INDSYC,
     3               NABSR,NABSC,
     4               K1P1,KDEL,
     5               NPARR,NPARC,
     6               DPFL)
      ENDIF
C
      IF (KDEL.EQ.(-1)) NFDP=NFIL41
      IF (KDEL.EQ.0   ) NFDP=NFIL42
      IF (KDEL.EQ.1   ) NFDP=NFIL43 
C
      CALL WRTDIP(WORK(IBLTOT),LREC35,
     1            MXFCCO,NFDP,NRCNR) 
C
9116  CONTINUE
9115  CONTINUE
9114  CONTINUE
9113  CONTINUE    
9112  CONTINUE
      IF (IDPTR.EQ.2) THEN
C        lediglich Dipolmomentmatrixelemente 
C        sollten berechnet werden
         GOTO 2323 
      ENDIF
      ENDIF
C 
c     NUNIT=NFIL34
c     OPEN(UNIT=NFIL34,STATUS='SCRATCH',ACCESS='DIRECT',
c    1     RECL=MXFCCO*MXFCCO*NRECUN,ERR=1009)
C
C     Verschiedene Berechnungen
C
      N2XS=NINT(XMULTI)-1
      ESS=DFLOAT(N2XS)/2.0D+00
      FACTS=SQRT((2.0D+00*ESS+1.0D+00)*ESS*(ESS+1.0D+00))
      J2RMAX=IABS(2*(JMAXP1-1)-N2XS)
      J2RMIN=MOD(N2XS,2)
      NMINMX=IABS(J2RMAX-N2XS)/2
      NMAXMX=    (J2RMAX+N2XS)/2
C
C     Maximale Laenge IEIGMX eines Eigenvektors ermitteln
C
      IEIGMX=0
      DO 8121 JX2VAL=J2RMIN,J2RMAX,2
         NMIN=IABS(JX2VAL-N2XS)/2
         NMAX=    (JX2VAL+N2XS)/2
         DO 8122 INDSYM=1,2*MAXPSY
            IEIGTP=0 
            CALL MATPOS(NMIN,NMAX,NMAX+1,WORK(INFCTS),
     1                  JMAXP1,INDSYM,IEIGTP)
            IF (IEIGTP.GT.IEIGMX) IEIGMX=IEIGTP
8122     CONTINUE       
8121  CONTINUE  
C
C     Dateien fuer Eigenvektoren fuer die 
C     Intensitaetsberechnungen oeffnen
C
      NUNIT=NFIL51
      MENREC=(IEIGMX+2)*NRECUN
      OPEN(UNIT=NFIL51,STATUS='SCRATCH',ACCESS='DIRECT',
     1     RECL=MENREC,ERR=1009)
      NUNIT=NFIL52
      OPEN(UNIT=NFIL52,STATUS='SCRATCH',ACCESS='DIRECT',
     1     RECL=MENREC,ERR=1009)
      INTFLG=1
      NRAUS=0 
      PARFUN=0   
C
C     Schleife ueber J
C
C
C  WE DO NOT LOOP OVER J; WE LOOP OVER 2*J (JX2) WHICH IS ALWAYS
C  INTEGRAL
C  THE HIGHEST POSSIBLE VALUE OF J IS ABS(NMAX-S). IN THE PROGRAM
C  NMAX IS UNFORTUNATELY CALLED JMAX. THE REAL JMAX*2 WE CALL J2RMAX
C
C
C**********************************************************
      DO 8120 JX2VAL=J2RMIN,J2RMAX,2
C**********************************************************
      XJVAL=DFLOAT(JX2VAL)/2.0D+00
C
      IEIGST=0
      IF (JX2VAL.EQ.J2RMIN) THEN
        INTFST=.TRUE.
      ELSE
        INTFST=.FALSE.
      ENDIF
      IF (INTFLG.EQ.1) THEN
        INTFLG=0
      ELSE
        INTFLG=1
      ENDIF
C
C     DETERMINE RANGE OF N VALUES NEEDED FOR THIS J VALUE
C
      NMIN=IABS(JX2VAL-N2XS)/2
      NMAX=    (JX2VAL+N2XS)/2
      NMINP1=NMIN+1
      NMAXP1=NMAX+1
C
      NUNIT=NFIL34
      OPEN(UNIT=NFIL34,STATUS='SCRATCH',ACCESS='DIRECT',
     1     RECL=MXFCCO*MXFCCO*NRECUN,ERR=1009)
C
C*********************************************************************
C
      DO 200 JP1=NMINP1,NMAXP1
C
      J=JP1-1
C
      DO 180 ITAUP1=1,2
C
      ITAU=ITAUP1-1
C
      IF (MOD(J,2) .EQ. 0) THEN
            IROTE=J/2+1
            IROTO=J/2
      ELSE
            IROTE=(J+1)/2
            IROTO=IROTE
      ENDIF
      IF (MOD(J+ITAU,2) .NE. 0) THEN
               IROTE1=IROTE-1
               IROTE2=IROTE
      ELSE
               IROTE1=IROTE
               IROTE2=IROTE-1
      ENDIF
      IROTO1=IROTO
      IROTO2=IROTO
C
C*********************************************************************
C
      DO 362 IFLAG=1,MAXPSY
C
C*********************************************************************
C
      IDX=JINDEX ( JX2VAL )
      IF (SKPSYM(ITAUP1,IFLAG,IDX,ISO)) GOTO 362
C
      INDSYM=2*(IFLAG-1)+ITAUP1
C
      DO 411 K1P1=1,JP1
      KP1=K1P1
C
      DO 420 KDEL=0,2,1
      K2P1=K1P1+KDEL
      IF (K2P1.GT.JP1) GOTO 420
C
C     MATRIZEN loeschen
C
      DO 421 II=0,(LREC35*LREC35-1)
         WORK(IHMAST+II)=0
421   CONTINUE
      DO 422 II=0,(LREC34*LREC34-1)
         WORK(IONERC+II)=0
422   CONTINUE
C
      IF (.NOT. SYMM) THEN
            NBLCK1=(IROTE1+IROTO1)*V2MP11*IQUANT(4,ISO)
            NBLCK2=(IROTE2+IROTO2)*V2MP12*IQUANT(8,ISO)
	    NDIMST=LSTYPA(1)
            NDIMRO=LSTYPA(1)
            NDIMCO=LSTYPA(1)
            ISTOFF=0
            JSTOFF=0
      ELSE
	    IF (IFLAG .EQ. 1) THEN
               NBLCK1=(IROTE1*IQUANT(4,ISO)+IROTO1*IQUANT(5,ISO))
     1               *V2MP11
               NBLCK2=(IROTE2*IQUANT(8,ISO)+IROTO2*IQUANT(9,ISO))
     1               *V2MP12
	    ELSE
               NBLCK1=(IROTE1*IQUANT(5,ISO)+IROTO1*IQUANT(4,ISO))
     1               *V2MP11
	       NBLCK2=(IROTE2*IQUANT(9,ISO)+IROTO2*IQUANT(8,ISO))
     1               *V2MP12
	    ENDIF
      ENDIF
      IRODIM=NBLCK1+NBLCK2
C
      IBOOK(JP1,ITAUP1,IFLAG,1)=IRODIM
      IBOOK(JP1,ITAUP1,IFLAG,2)=NBLCK1
      IBOOK(JP1,ITAUP1,IFLAG,3)=IROTE1
      IBOOK(JP1,ITAUP1,IFLAG,4)=IROTO1
      IBOOK(JP1,ITAUP1,IFLAG,5)=IROTE2
      IBOOK(JP1,ITAUP1,IFLAG,6)=IROTO2
C
C     BERECHNUNG DER BLOCKGROESSEN (g.o)
C
      ITAU=ITAUP1-1
      NABLKR=IQUANT(2,ISO)+1
      NBBLKR=IQUANT(6,ISO)+1
      IF (SYMM) THEN
        IF (MOD(MOD(K1P1-1,2)+(IFLAG-1),2).EQ.0) THEN
          NABLKR=NABLKR*IQUANT(4,ISO)
          NBBLKR=NBBLKR*IQUANT(8,ISO)
        ELSE
          NABLKR=NABLKR*IQUANT(5,ISO)
          NBBLKR=NBBLKR*IQUANT(9,ISO)
        ENDIF
      ELSE
        NABLKR=NABLKR*IQUANT(4,ISO)
        NBBLKR=NBBLKR*IQUANT(8,ISO)
      ENDIF
      IF (K1P1.EQ.1) THEN
        IF (MOD(J+ITAU,2).EQ.0) THEN
          NBBLKR=0
        ELSE
          NABLKR=0
        ENDIF
      ENDIF
      NABLKC=IQUANT(2,ISO)+1
      NBBLKC=IQUANT(6,ISO)+1
      IF (SYMM) THEN
        IF (MOD(MOD(K2P1-1,2)+(IFLAG-1),2).EQ.0) THEN
          NABLKC=NABLKC*IQUANT(4,ISO)
          NBBLKC=NBBLKC*IQUANT(8,ISO)
        ELSE
          NABLKC=NABLKC*IQUANT(5,ISO)
          NBBLKC=NBBLKC*IQUANT(9,ISO)
        ENDIF
      ELSE
        NABLKC=NABLKC*IQUANT(4,ISO)
        NBBLKC=NBBLKC*IQUANT(8,ISO)
      ENDIF
      IF (K2P1.EQ.1) THEN
        IF (MOD(J+ITAU,2).EQ.0) THEN
          NBBLKC=0
        ELSE
          NABLKC=0
        ENDIF
      ENDIF

C*********************************************************************
C
C  CONSTRUCT DIAGONAL BLOCK FOR SURFACES 1 AND 2
C
C*********************************************************************
C
C
      REWIND NFIL9
      REWIND NFIL10
      REWIND NFIL33
C
      DO 323 ISURF=1,2
C
      NSURF=ISURF
      IF (NSURF.EQ.1) THEN
        NBLKR=NABLKR
        NBLKC=NABLKC
      ELSE
        NBLKR=NBBLKR
        NBLKC=NBBLKC
      ENDIF
      IF (ISURF .EQ. 1) THEN
         V2MAX=V2MP11-1
         V2MXP1=V2MAX+1
         LSTYPA(1)=IQUANT(4,ISO)
         LSTYPA(2)=IQUANT(5,ISO)
         KSTYPA(1)=KSTY11
         KSTYPA(2)=KSTY21
         NOBAS=NOBAL
         ITAU=ITAUP1-1
      ELSE
         V2MAX=V2MP12-1
         V2MXP1=V2MAX+1
         LSTYPA(1)=IQUANT(8,ISO)
         LSTYPA(2)=IQUANT(9,ISO)
         KSTYPA(1)=KSTY12
         KSTYPA(2)=KSTY22
         NOBAS=NOBAR
         ITAU=2-ITAUP1
      ENDIF
C
      ISYMS0=ISMAST+LISTRT(ISURF)
      IF (ISYMS0 .GE. MWORK) GOTO 1010
      READ (NFIL10) (WORK(III), III=ISMAST,ISYMS0-1)
      IASYS0=ISYMS0+LZERSY(ISURF)
      IF (IASYS0 .GE. MWORK) GOTO 1010
      READ (NFIL9) (WORK(III), III=ISYMS0,IASYS0-1)
      IF (.NOT.SYMM) THEN
          ISYMS2=IASYS0+LZERAS(ISURF)
          IF (ISYMS2 .GE. MWORK) GOTO 1010
          READ (NFIL9) (WORK(III), III=IASYS0,ISYMS2-1)
      ELSE
          ISYMS2=IASYS0
      ENDIF
      IASYS2=ISYMS2+LTWOSY(ISURF)
      IF (IASYS2 .GE. MWORK) GOTO 1010
      READ (NFIL9) (WORK(III), III=ISYMS2,IASYS2-1)
      IF (.NOT.SYMM) THEN
          IENRGS=IASYS2+LTWOAS(ISURF)
          IF (IENRGS .GE. MWORK) GOTO 1010
          READ (NFIL9) (WORK(III), III=IASYS2,IENRGS-1)
      ELSE
          IENRGS=IASYS2
      ENDIF
      IESTST=IENRGS+LEBEND(ISURF)
      IF (IESTST .GE. MWORK) GOTO 1010
      READ (NFIL33) (WORK(III), III=IENRGS,IESTST-1)
      IWRKEN=IESTST+LESTRT(ISURF)
      IF (IWRKEN .GE. MWORK) GOTO 1010
      READ (NFIL33) (WORK(III), III=IESTST,IWRKEN-1)
C
C*********************************************************************
C
C  ADJUST ENERGY ZERO FOR SURFACE 2 RELATIVE TO SURFACE 1
C
C*********************************************************************
C
      IF (ISURF .EQ. 2) THEN
      DO 373 III=IESTST,IWRKEN-1
373   WORK(III)=WORK(III)+VMINS2-VMINS1
      ENDIF
C
C*********************************************************************
C
C  SET UP THE EVEN K CORNER OF THE HAMILTONIAN MATRIX BLOCK
C
C*********************************************************************
C
      IF (ISURF .EQ. 1) THEN
         IROTE=IROTE1
         IROTO=IROTO1
         NOFSET=0
      ELSE
         IROTE=IROTE2
         IROTO=IROTO2
         NOFSET=NABLKR
      ENDIF
      IEO=0
      IF (IFLAG .EQ. 1) THEN
         NDIMST=LSTYPA(1)
         ISTOFF=0
      ELSE
         NDIMST=LSTYPA(2)
         ISTOFF=KSTYPA(1)
      ENDIF
      IF (ISURF .EQ. 1) THEN
         IBOOK(JP1,ITAUP1,IFLAG,7)=NDIMST
         IBOOK(JP1,ITAUP1,IFLAG,11)=V2MXP1
      ELSE
         IBOOK(JP1,ITAUP1,IFLAG,9)=NDIMST
         IBOOK(JP1,ITAUP1,IFLAG,12)=V2MXP1
      ENDIF
      IF (((NBLKR*NBLKC).NE. 0).AND.(KDEL.EQ.0)) THEN
      CALL DELK0 ( WORK(IENRGS) , WORK(IESTST) , WORK(ISYMS0) ,
     1             WORK(IASYS0) , WORK(ISMAST) , WORK(IHMAST) ,
     2             NOFSET , NDIMST , JP1,KP1    , IRODIM , IEO    ,
     3             ITAU   , LZERSY(ISURF) , LZERAS(ISURF) , ISTOFF ,
     1             NOBAS  , NOPERS , NZERSY , NZERAS , LREC35)
      ENDIF
      IF (ISURF .EQ. 1) THEN
         NOFSRO=0
         NOFSCO=0
      ELSE
         NOFSRO=NABLKR
         NOFSCO=NABLKC
      ENDIF
      IF (IFLAG .EQ. 1) THEN
         NDIMRO=LSTYPA(1)
         NDIMCO=LSTYPA(1)
         ISTOFF=0
         JSTOFF=0
      ELSE
         NDIMRO=LSTYPA(2)
         NDIMCO=LSTYPA(2)
         ISTOFF=KSTYPA(1)
         JSTOFF=KSTYPA(1)
      ENDIF

      IF (((NBLKR*NBLKC).NE. 0).AND.(KDEL.EQ.2)) THEN
      CALL DELK2 (                               WORK(ISYMS2) ,
     1             WORK(IASYS2) , WORK(ISMAST) , WORK(IHMAST) ,
     2             NOFSRO , NOFSCO , NDIMRO , NDIMCO , JP1    ,
     3             KP1,KDEL,
     3             IRODIM , IEO    , ITAU   , LTWOSY(ISURF) ,
     4             LTWOAS(ISURF) , ISTOFF , JSTOFF ,
     1             NOBAS  , NOPERS , NTWOSY , NTWOAS , LREC35)
      ENDIF
C
C*********************************************************************
C
C  SET UP THE ODD K CORNER OF THE HAMILTONIAN MATRIX BLOCK
C
C*********************************************************************
C
320   CONTINUE
      IEO=1
      IF (ISURF .EQ. 1) THEN
         NOFSET=0
      ELSE
         NOFSET=NABLKR
      ENDIF
      IF (IFLAG .EQ. 1 .AND. SYMM) THEN
            NDIMST=LSTYPA(2)
            ISTOFF=KSTYPA(1)
      ELSE
            NDIMST=LSTYPA(1)
            ISTOFF=0
      ENDIF
      IF (ISURF .EQ. 1) THEN
         IBOOK(JP1,ITAUP1,IFLAG,8)=NDIMST
      ELSE
         IBOOK(JP1,ITAUP1,IFLAG,10)=NDIMST
      ENDIF
      IF (((NBLKR*NBLKC) .NE. 0).AND.(KDEL.EQ.0)) THEN
      CALL DELK0 ( WORK(IENRGS) , WORK(IESTST) , WORK(ISYMS0) ,
     1             WORK(IASYS0) , WORK(ISMAST) , WORK(IHMAST) ,
     2             NOFSET , NDIMST , JP1,KP1, IRODIM , IEO    ,
     3             ITAU   , LZERSY(ISURF) , LZERAS(ISURF) , ISTOFF ,
     1             NOBAS  , NOPERS , NZERSY , NZERAS , LREC35)
      ENDIF
      IF (ISURF.EQ.1) THEN
        NOFSRO=0
        NOFSCO=0
      ELSE
        NOFSRO=NABLKR
        NOFSCO=NABLKC
      ENDIF
      IF (IFLAG .EQ. 1 .AND. SYMM) THEN
            NDIMRO=LSTYPA(2)
            NDIMCO=LSTYPA(2)
            ISTOFF=KSTYPA(1)
            JSTOFF=KSTYPA(1)
      ELSE
	   NDIMRO=LSTYPA(1)
           NDIMCO=LSTYPA(1)
	   ISTOFF=0
           JSTOFF=0
      ENDIF
      IF (((NBLKR*NBLKC) .NE. 0).AND.(KDEL.EQ.2)) THEN
      CALL DELK2 (                               WORK(ISYMS2) ,
     1             WORK(IASYS2) , WORK(ISMAST) , WORK(IHMAST) ,
     2             NOFSRO , NOFSCO , NDIMRO , NDIMCO , JP1    ,
     3             KP1,KDEL,
     3             IRODIM , IEO    , ITAU   , LTWOSY(ISURF) ,
     4             LTWOAS(ISURF) , ISTOFF , JSTOFF ,
     1             NOBAS  , NOPERS , NTWOSY , NTWOAS , LREC35)
      ENDIF
C
C*********************************************************************
C
C  SET UP DELTA K = 1 MATRIX ELEMENTS CONNECTING THE TWO CORNERS
C
C*********************************************************************
C
341   CONTINUE
      ISYMS1=ISYMS0
      IASYS1=ISYMS1+LONESY(ISURF)
      IF (IASYS1 .GE. MWORK) GOTO 1010
      READ (NFIL9) (WORK(III), III=ISYMS1,IASYS1-1)
      IF (.NOT.SYMM) THEN
          IWRKEN=IASYS1+LONEAS(ISURF)
          IF (IWRKEN .GE. MWORK) GOTO 1010
          READ (NFIL9) (WORK(III), III=IASYS1,IWRKEN-1)
      ELSE
          IWRKEN=IASYS1
      ENDIF
      IF (IWRKEN .GE. MWORK) GOTO 1010
      IF (ISURF .EQ. 1) THEN
         NOFSRO=0
         NOFSCO=0
      ELSE
         NOFSRO=NABLKR
         NOFSCO=NABLKC
      ENDIF
C *****  DO DELTA K = 1 MATRIX ELEMENTS FOR LEFT K EVEN
      IEO=0
      IF (SYMM) THEN
          IF (IFLAG .EQ. 2) THEN
	       NDIMRO=LSTYPA(2)
               NDIMCO=LSTYPA(1)
	       ISTOFF=KSTYPA(1)
	       JSTOFF=0
           ELSE
               NDIMRO=LSTYPA(1)
               NDIMCO=LSTYPA(2)
               ISTOFF=0
               JSTOFF=KSTYPA(1)
           ENDIF
      ELSE
               NDIMRO=LSTYPA(1)
               NDIMCO=LSTYPA(1)
               ISTOFF=0
               JSTOFF=0
      ENDIF
      IF (((NBLKR*NBLKC).NE.0).AND.(KDEL.EQ.1)) THEN
      CALL DELK1 (                               WORK(ISYMS1) ,
     1             WORK(IASYS1) , WORK(ISMAST) , WORK(IHMAST) ,
     2             NOFSRO , NOFSCO , NDIMRO , NDIMCO , JP1    ,
     3             KP1,KDEL,
     3             IRODIM , IEO    , ITAU   , LONESY(ISURF) ,
     4             LONEAS(ISURF) , ISTOFF , JSTOFF ,
     1             NOBAS  , NOPERS , NONESY , NONEAS , LREC35)
      ENDIF
C *****  DO DELTA K = 1 MATRIX ELEMENTS FOR LEFT K ODD
358   IEO=1
      IF (SYMM) THEN
          IF (IFLAG .EQ. 1) THEN
	       NDIMRO=LSTYPA(2)
               NDIMCO=LSTYPA(1)
	       ISTOFF=KSTYPA(1)
	       JSTOFF=0
           ELSE
               NDIMRO=LSTYPA(1)
               NDIMCO=LSTYPA(2)
               ISTOFF=0
               JSTOFF=KSTYPA(1)
           eNDIF
      ENDIF
      IF (((NBLKR*NBLKC).NE.0).AND.(KDEL.EQ.1)) THEN
      CALL DELK1 (                               WORK(ISYMS1) ,
     1             WORK(IASYS1) , WORK(ISMAST) , WORK(IHMAST) ,
     2             NOFSRO , NOFSCO , NDIMRO , NDIMCO , JP1    ,
     3             KP1,KDEL,
     3             IRODIM , IEO    , ITAU   , LONESY(ISURF) ,
     4             LONEAS(ISURF) , ISTOFF , JSTOFF ,
     1             NOBAS  , NOPERS , NONESY , NONEAS , LREC35 )
      ENDIF
360   CONTINUE
323   CONTINUE
C
C*********************************************************************
C
C  CALCULATE THE RENNER-TELLER INTERACTION TERMS
C
C*********************************************************************
C
C  K EVEN -- K EVEN BLOCK
C
      ITAU = ITAUP1 - 1
C
      ISMART=ISMAST
      ISYFC0=ISMART+LUINTS
      IF (ISYFC0 .GE. MWORK) GOTO 1010
      READ (NFIL10) (WORK(III), III=ISMART,ISYFC0-1)
      IASFC0=ISYFC0+IZERSY
      IF (ITRNSY .GE. MWORK) GOTO 1010
      READ (NFIL9) (WORK(III), III=ISYFC0,IASFC0-1)
C
      IF (.NOT.SYMM) THEN
	    ISYFC2=IASFC0+IZERAS
            IF (ISYFC2.GE. MWORK) GOTO 1010
            READ (NFIL9) (WORK(III), III=IASFC0,ISYFC2-1)
      ELSE
	    ISYFC2=IASFC0
      ENDIF
C
      IASFC2=ISYFC2+ITWOSY
      IF (IASFC2 .GE. MWORK) GOTO 1010
      READ (NFIL9) (WORK(III), III=ISYFC2,IASFC2-1)
C
      IF (.NOT.SYMM) THEN
	    IWRKEN=IASFC2+ITWOAS
            IF (IWRKEN .GE. MWORK) GOTO 1010
            READ (NFIL9) (WORK(III), III=IASFC2,IWRKEN-1)
      ELSE
	    IWRKEN=IASFC2
            IF (IWRKEN .GE. MWORK) GOTO 1010
      ENDIF
C
      IF (((NABLKR*NBBLKC) .GT. 0).OR.((NABLKC*NBBLKR).GT.0)) THEN
        NOFSRO=0
        NOFSCO=NABLKC
        IF (IFLAG .EQ. 1) THEN
          NDIMRO=IQUANT(4,ISO)
          NDIMCO=IQUANT(8,ISO)
          ISTOFF=0
          JSTOFF=0
        ELSE
          NDIMRO=IQUANT(5,ISO)
          NDIMCO=IQUANT(9,ISO)
          ISTOFF=KSTY11
          JSTOFF=KSTY12
        ENDIF
        IEO=0
C
      IF (KDEL.EQ.0) THEN
      CALL RNPM0 ( WORK(ISYFC0) , WORK(IASFC0) ,
     1             WORK(ISMART) ,
     1             WORK(IHMAST) , NOFSRO , NOFSCO ,
     1             NDIMRO , NDIMCO , V2MP11 , V2MP12 , JP1 ,KP1,
     2             IRODIM , IEO    , ITAU   , IZERSY , IZERAS ,
     3             ISTOFF , JSTOFF , NOBAL  , NOBAR  , NINTOP ,
     4             NINSYM,NINASY,LREC35)
      ENDIF
      IF (KDEL.EQ.2) THEN
      IF ((NABLKR*NBBLKC).GT.0) THEN
      CALL RNPM2 ( WORK(ISYFC2) , WORK(IASFC2) , WORK(ISMART) ,
     1             WORK(IHMAST) , NOFSRO , NOFSCO ,
     1             NDIMRO , NDIMCO , V2MP11 , V2MP12 , JP1 ,
     3             KP1,KDEL,
     2             IRODIM , IEO    , ITAU   , ITWOSY , ITWOAS ,
     3             ISTOFF , JSTOFF , NOBAL  , NOBAR  , NINTOP ,
     4             NTWOSY , NTWOAS , LREC35)
      ENDIF
      IF ((NBBLKR*NABLKC).GT.0) THEN
        NOFSRO=0
        NOFSCO=0
      CALL RNPM2 ( WORK(ISYFC2) , WORK(IASFC2) , WORK(ISMART) ,
     1             WORK(IONERC) , NOFSRO , NOFSCO ,
     1             NDIMRO , NDIMCO , V2MP11 , V2MP12 , JP1 ,
     3             K2P1,(-1)*KDEL,
     2             IRODIM , IEO    , ITAU   , ITWOSY , ITWOAS ,
     3             ISTOFF , JSTOFF , NOBAL  , NOBAR  , NINTOP ,
     4             NTWOSY , NTWOAS , LREC34)
      ENDIF
      ENDIF
C
C
C  K ODD -- K ODD BLOCK
C
      IF (IFLAG .EQ. 1) THEN
      NOFSRO=0
      NOFSCO=NABLKC
      IF (.NOT.SYMM) THEN
           NDIMRO=IQUANT(4,ISO)
           NDIMCO=IQUANT(8,ISO)
           ISTOFF=0
           JSTOFF=0
      ELSE
           NDIMRO=IQUANT(5,ISO)
           NDIMCO=IQUANT(9,ISO)
           ISTOFF=KSTY11
           JSTOFF=KSTY12
      ENDIF
      ELSE
           NOFSRO=0
	   NOFSCO=NABLKC
	   NDIMRO=IQUANT(4,ISO)
	   NDIMCO=IQUANT(8,ISO)
	   ISTOFF=0
	   JSTOFF=0
      ENDIF
      IEO=1
C
      IF (KDEL.EQ.0) THEN
      CALL RNPM0 ( WORK(ISYFC0) , WORK(IASFC0) ,
     1             WORK(ISMART) ,
     1             WORK(IHMAST) , NOFSRO , NOFSCO ,
     1             NDIMRO , NDIMCO , V2MP11 , V2MP12 , JP1 ,
     2             KP1,
     2             IRODIM , IEO    , ITAU   , IZERSY , IZERAS ,
     3             ISTOFF , JSTOFF , NOBAL  , NOBAR  , NINTOP ,
     4             NINSYM , NINASY, LREC35)
      ENDIF
      IF (KDEL.EQ.2) THEN
      IF ((NABLKR*NBBLKC).GT.0) THEN
      CALL RNPM2 ( WORK(ISYFC2) , WORK(IASFC2) , WORK(ISMART) ,
     1             WORK(IHMAST) , NOFSRO , NOFSCO ,
     1             NDIMRO , NDIMCO , V2MP11 , V2MP12 , JP1 ,
     3             KP1,KDEL,
     2             IRODIM , IEO    , ITAU   , ITWOSY , ITWOAS ,
     3             ISTOFF , JSTOFF , NOBAL  , NOBAR  , NINTOP ,
     4             NTWOSY , NTWOAS , LREC35)
      ENDIF
      IF ((NBBLKR*NABLKC).GT.0) THEN
        NOFSRO=0
        NOFSCO=0
      CALL RNPM2 ( WORK(ISYFC2) , WORK(IASFC2) , WORK(ISMART) ,
     1             WORK(IONERC) , NOFSRO , NOFSCO ,
     1             NDIMRO , NDIMCO , V2MP11 , V2MP12 , JP1 ,
     3             K2P1,(-1)*KDEL,
     2             IRODIM , IEO    , ITAU   , ITWOSY , ITWOAS ,
     3             ISTOFF , JSTOFF , NOBAL  , NOBAR  , NINTOP ,
     4             NTWOSY , NTWOAS , LREC34)
      ENDIF
      ENDIF
C
      ISYFC1=ISYFC0
      IASFC1=ISYFC1+IONESY
      IF (IASFC1 .GE. MWORK) GOTO 1010
      READ (NFIL9) (WORK(III), III=ISYFC1,IASFC1-1)
C
      IF (.NOT.SYMM) THEN
	    IWRKEN=IASFC1+IONEAS
            IF (IWRKEN .GE. MWORK) GOTO 1010
            READ (NFIL9) (WORK(III), III=IASFC1,IWRKEN-1)
      ELSE
	    IWRKEN=IASFC1
            IF (IWRKEN .GE. MWORK) GOTO 1010
      ENDIF
C
C
C  K ODD -- K EVEN BLOCK
C
      IF (IFLAG .EQ. 1) THEN
      NOFSRO=0
      NOFSCO=NABLKC
      IF (.NOT.SYMM) THEN
           NDIMRO=IQUANT(4,ISO)
           NDIMCO=IQUANT(8,ISO)
           ISTOFF=0
           JSTOFF=0
      ELSE
           NDIMRO=IQUANT(5,ISO)
           NDIMCO=IQUANT(8,ISO)
           ISTOFF=KSTY11
           JSTOFF=0
      ENDIF
      ELSE
           NOFSRO=0
           NOFSCO=NABLKC
	   NDIMRO=IQUANT(4,ISO)
	   NDIMCO=IQUANT(9,ISO)
	   ISTOFF=0
	   JSTOFF=KSTY12
      ENDIF
      IEO=1
C
      IF (KDEL.EQ. 1) THEN
      IF ((NABLKR*NBBLKC).GT.0) THEN
      CALL RNPM1 ( WORK(ISYFC1) , WORK(IASFC1) , WORK(ISMART) ,
     1             WORK(IHMAST) , NOFSRO , NOFSCO , NDIMRO , NDIMCO ,
     1             V2MP11 , V2MP12 , JP1,KP1,KDEL    , IRODIM ,
     2             IEO    , ITAU   , IONESY , IONEAS ,
     1             ISTOFF , JSTOFF , NOBAL , NOBAR , NINTOP ,
     1             NONESY , NONEAS , LREC35)
      ENDIF
      IF ((NBBLKR*NABLKC).GT.0) THEN
      IF (IFLAG .EQ. 1) THEN
      NOFSRO=0
      NOFSCO=NABLKC
      IF (.NOT.SYMM) THEN
           NDIMRO=IQUANT(4,ISO)
           NDIMCO=IQUANT(8,ISO)
           ISTOFF=0
           JSTOFF=0
      ELSE
           NDIMRO=IQUANT(4,ISO)
           NDIMCO=IQUANT(9,ISO)
           ISTOFF=0
           JSTOFF=KSTY12
      ENDIF
      ELSE
	    NOFSRO=0
            NOFSCO=NABLKC
	    NDIMRO=IQUANT(5,ISO)
	    NDIMCO=IQUANT(8,ISO)
	    ISTOFF=KSTY11
	    JSTOFF=0
      ENDIF
      IEO=0
        NOFSRO=0
        NOFSCO=0
      CALL RNPM1 ( WORK(ISYFC1) , WORK(IASFC1) , WORK(ISMART) ,
     1             WORK(IONERC) , NOFSRO , NOFSCO , NDIMRO , NDIMCO ,
     1             V2MP11 , V2MP12 , JP1,
     3             K2P1,(-1)*KDEL    , IRODIM ,
     2             IEO    , ITAU   , IONESY , IONEAS ,
     1             ISTOFF , JSTOFF , NOBAL , NOBAR , NINTOP ,
     1             NONESY , NONEAS , LREC34)

      ENDIF
      ENDIF
C
C
C  K EVEN -- K ODD BLOCK
C
      IF (IFLAG .EQ. 1) THEN
      NOFSRO=0
      NOFSCO=NABLKC
      IF (.NOT.SYMM) THEN
           NDIMRO=IQUANT(4,ISO)
           NDIMCO=IQUANT(8,ISO)
           ISTOFF=0
           JSTOFF=0
      ELSE
           NDIMRO=IQUANT(4,ISO)
           NDIMCO=IQUANT(9,ISO)
           ISTOFF=0
           JSTOFF=KSTY12
      ENDIF
      ELSE
	    NOFSRO=0
            NOFSCO=NABLKC
	    NDIMRO=IQUANT(5,ISO)
	    NDIMCO=IQUANT(8,ISO)
	    ISTOFF=KSTY11
	    JSTOFF=0
      ENDIF
      IEO=0
C
      IF (KDEL.EQ.1) THEN
      IF ((NABLKR*NBBLKC).GT.0) THEN
      CALL RNPM1 ( WORK(ISYFC1) , WORK(IASFC1) , WORK(ISMART) ,
     1             WORK(IHMAST) , NOFSRO , NOFSCO , NDIMRO , NDIMCO ,
     1             V2MP11 , V2MP12 , JP1,KP1,KDEL, IRODIM ,
     2             IEO    , ITAU   , IONESY , IONEAS ,
     1             ISTOFF , JSTOFF , NOBAL , NOBAR , NINTOP ,
     1             NONESY , NONEAS , LREC35)
      ENDIF
      IF ((NBBLKR*NABLKC).GT.0) THEN
      IF (IFLAG .EQ. 1) THEN
      NOFSRO=0
      NOFSCO=NABLKC
      IF (.NOT.SYMM) THEN
           NDIMRO=IQUANT(4,ISO)
           NDIMCO=IQUANT(8,ISO)
           ISTOFF=0
           JSTOFF=0
      ELSE
           NDIMRO=IQUANT(5,ISO)
           NDIMCO=IQUANT(8,ISO)
           ISTOFF=KSTY11
           JSTOFF=0
      ENDIF
      ELSE
           NOFSRO=0
           NOFSCO=NABLKC
	   NDIMRO=IQUANT(4,ISO)
	   NDIMCO=IQUANT(9,ISO)
	   ISTOFF=0
	   JSTOFF=KSTY12
      ENDIF
      IEO=1
        NOFSRO=0
        NOFSCO=0
      CALL RNPM1 ( WORK(ISYFC1) , WORK(IASFC1) , WORK(ISMART) ,
     1             WORK(IONERC) , NOFSRO , NOFSCO , NDIMRO , NDIMCO ,
     1             V2MP11 , V2MP12 , JP1,
     3             K2P1,(-1)*KDEL, IRODIM ,
     2             IEO    , ITAU   , IONESY , IONEAS ,
     1             ISTOFF , JSTOFF , NOBAL , NOBAR , NINTOP ,
     1             NONESY , NONEAS , LREC34)
      ENDIF
      ENDIF
C
      ENDIF
C
C     ba-BLOCK anfuegen
C
      CALL BABLK(NABLKR,NABLKC,
     1           NBBLKR,NBBLKC,
     2           KDEL,
     3           WORK(IHMAST),
     4           WORK(IONERC),
     5           LREC34,
     6           LREC35)
C
C     AUSGABE
C
C    1           'JP1=',JP1,'K1P1=',K1P1,'KDEL=',KDEL
C     CALL AUS(WORK(IHMAST),LREC35,NFIL6)
C
C     TRANSFORMIERN
C
      DIMR=NABLKR+NBBLKR
      DIMC=NABLKC+NBBLKC
C
      IWRKT = ISMAST
      IWRKW  = IWRKT  + LREC35*MXFCCO
      IWRKEN = IWRKW  + LREC35*MXFCCO
      IF (IWRKEN.GT.MWORK) GOTO 1010
C
      CALL BSSTRN(JP1,K1P1,JP1,K2P1,
     1       INDSYM,INDSYM,
     2       WORK(IPSFCT),WORK(INFCTS),
     3       JMAXP1,
     4       NROW,NCOL,
     4       DIMR,DIMC,
     5       WORK(IHMAST),WORK(IWRKT),WORK(IWRKW),
     6       LREC35,MXFCCO,
     7       NFIL35)
C
C     AUSGABE
C
c    1           'JP1=',JP1,'K1P1=',K1P1,'KDEL=',KDEL,
C    2           'NROW=',NROW,'NCOL=',NCOL
c     CALL AUS(WORK(IHMAST),LREC35,NFIL6)
C
C
C     SPEICHERN
C
      CALL RECPOS(JMAXP1,JP1,KP1,INDSYM,KDEL,NMINP1,IREC)
      CALL SCHREI(IREC,
     1            WORK(IHMAST),
     2            LREC35,MXFCCO,
     3            NFIL34)
C
C     ENDE DER SCHLEIFEN
C
420   CONTINUE
411   CONTINUE
362   CONTINUE
C
180   CONTINUE
200   CALL PRTIMN( J , OLDTIM , OLDVEC )
C
C
C
C  WE START BY READING IN THE BEND AND STRETCH INTEGRALS
C
      REWIND NFIL10
      REWIND NFIL29
C
      ISMAPP=ISMAST
      ISMSPP=ISMAPP+LISTRT(1)
      IF (ISMSPP .GE. MWORK) GOTO 1010
      READ (NFIL10) (WORK(III), III=ISMAPP,ISMSPP-1)
      IASMPP=ISMSPP+LSPOSY(1)
      IF (IASMPP .GE. MWORK) GOTO 1010
      READ (NFIL29) (WORK(III), III=ISMSPP,IASMPP-1)
      IF (.NOT.SYMM) THEN
          IWRKEN=IASMPP+LSPOAS(1)
          IF (IWRKEN .GE. MWORK) GOTO 1010
          READ (NFIL29) (WORK(III), III=IASMPP,IWRKEN-1)
      ELSE
          IWRKEN=IASMPP
      ENDIF
C
      ISMAMM=IWRKEN
      ISMSMM=ISMAMM+LISTRT(2)
      IF (ISMSMM .GE. MWORK) GOTO 1010
      READ (NFIL10) (WORK(III), III=ISMAMM,ISMSMM-1)
      IASMMM=ISMSMM+LSPOSY(2)
      IF (IASMMM .GE. MWORK) GOTO 1010
      READ (NFIL29) (WORK(III), III=ISMSMM,IASMMM-1)
      IF (.NOT.SYMM) THEN
          IWRKEN=IASMMM+LSPOAS(2)
          IF (IWRKEN .GE. MWORK) GOTO 1010
          READ (NFIL29) (WORK(III), III=IASMMM,IWRKEN-1)
      ELSE
          IWRKEN=IASMMM
      ENDIF
C
      ISMAMP=IWRKEN
      ISMSMP=ISMAMP+LUINTS
      IF (ISMSMP .GE. MWORK) GOTO 1010
      READ (NFIL10) (WORK(III), III=ISMAMP,ISMSMP-1)
      IASMMP=ISMSMP+ISPOSY
      IF (IASMMP .GE. MWORK) GOTO 1010
      READ (NFIL29) (WORK(III), III=ISMSMP,IASMMP-1)
C
      IF (.NOT.SYMM) THEN
	    IWRKEN=IASMMP+ISPOAS
            IF (IWRKEN .GE. MWORK) GOTO 1010
            READ (NFIL29) (WORK(III), III=IASMMP,IWRKEN-1)
      ELSE
	    IWRKEN=IASMMP
      ENDIF
C
C     Transformationsmatrix
C
      IWRKT = IWRKEN
      IF (IWRKT.GT.MWORK) GOTO 1010
C
C     K-Matrix ( a und b )
C
      IWRKK = IWRKT + LREC35 * MXFCCO
      IF (IWRKK.GT.MWORK) GOTO 1010
C
C     Hilfsmatrix fuer die Transformation
C
      IWRKW = IWRKK + LREC35 * LREC35
      IF (IWRKW.GT.MWORK) GOTO 1010
C
C     Einzelner K-Block ( a oder b )
C
      IWRKKB = IWRKW + LREC35 * MXFCCO
C
C     Maximalgroesse der Gesamtmatrix ermitteln
C
      DO 8201 II=1,4
        IRO(II)=0
8201  CONTINUE
      DO 8200 ITAUP1=1,2
      DO 8200 IPSYM=1,MAXPSY
        INDSYM=2*(IPSYM-1)+ITAUP1
        CALL MATPOS(NMIN,NMAX,NMAX+1,WORK(INFCTS),
     1              JMAXP1,INDSYM,IRO(INDSYM))
8200  CONTINUE
      IROMAX=MAX0(IRO(1),IRO(2),IRO(3),IRO(4))
C
C     Platz fuer Eigenwerte und andere Daten schaffen
C
      IWRKED = IWRKKB + LREC35 * LREC35
      IF (IWRKED.GT.MWORK) GOTO 1010
C
C     SYMMETRY -- LOOPS OVER ITAUP1 AND IPSYM
C
      DO 8110 ITAUP1=1,2
C
      DO 8100 IPSYM=1,MAXPSY
      IDX=JINDEX ( JX2VAL )
      IF (SKPSYM(ITAUP1,IPSYM,IDX,ISO)) GOTO 8100
C
C
      INDSYM=2*(IPSYM-1)+ITAUP1
C
C     Groesse IRODIM der Gesamtmatrix ermitteln
C
      CALL MATPOS(NMIN,NMAX,NMAX+1,WORK(INFCTS),JMAXP1,INDSYM,IRODIM)
C
C     Platz fuer diese Matrix schaffen
C
      IGMAST = IWRKED + 2*MAXPSY*IROMAX*2
      IF ((IGMAST+IRODIM*IRODIM).GT.MWORK) GOTO 1010
      DO 8130 I=IGMAST,IGMAST+IRODIM*IRODIM
         WORK(I)=0.0D+00
8130  CONTINUE
C
C     Work-arrays fuer die Diagonalisierung
C
      IWRKAR=IGMAST + IRODIM * IRODIM
      IWORK= IWRKAR + 8*IRODIM
      IWRKFA=IWORK  + 5*IRODIM
      IF (IWORK.GT.MWORK) GOTO 1010
C
C     Platz fuer Eigenwerte und Eigenvektoren
C
      IWEVAL=IWRKFA + IRODIM
      IWEVEC=IWEVAL + IRODIM
      IINTST=IWEVEC+IRODIM*IRODIM
      IF (IINTST.GT.MWORK) GOTO 1010
C
C     Dimensionen eines K-Blocks berechnen
C
      NMINP1=NMIN+1
      IABLKE=IBOOK(NMINP1,ITAUP1,IPSYM, 7)
     1      *IBOOK(NMINP1,ITAUP1,IPSYM,11)
      IABLKO=IBOOK(NMINP1,ITAUP1,IPSYM, 8)
     1      *IBOOK(NMINP1,ITAUP1,IPSYM,11)
      IBBLKE=IBOOK(NMINP1,ITAUP1,IPSYM, 9)
     1      *IBOOK(NMINP1,ITAUP1,IPSYM,12)
      IBBLKO=IBOOK(NMINP1,ITAUP1,IPSYM,10)
     1      *IBOOK(NMINP1,ITAUP1,IPSYM,12)
      KDIMMX=MAX0(IABLKE+IBBLKE,IABLKO+IBBLKO)
C
C     Matrix aufstellen
C
      CALL MAKMAT(WORK(IWRKT),
     1       WORK(IGMAST),WORK(IWRKK),WORK(IWRKW),WORK(IWRKKB),
     2       IRODIM,LREC34,LREC35,NMIN,NMAX,
     3       IABLKE,IBBLKE,IABLKO,IBBLKO,
     4       ITAUP1,INDSYM,WORK(INFCTS),
     5       WORK(ISMAPP),WORK(ISMSPP),WORK(IASMPP),
     5       WORK(ISMAMM),WORK(ISMSMM),WORK(IASMMM),
     5       WORK(ISMAMP),WORK(ISMSMP),WORK(IASMMP),
     6       LSPOSY(1),LSPOAS(1),
     7       LSPOSY(2),LSPOAS(2),
     8       ISPOSY,ISPOAS,
     9       NOBAL,NOBAR,
     1       NOPERS,NINTOP,
     2       FACTS,JX2VAL,N2XS,ESS,XJVAL,
     3       IBOOK,MXJVP1,IPSYM,
     4       KSTY11,KSTY12,
     5       NSPOSY,NSPOAS,SYMM,
     6       WORK(IPSFCT),MXFCCO,NRFCCO)
C
C     Matrix diagonalisieren und Resultate speichern
C
      CALL DIAROT(WORK(IGMAST),
     1            WORK(IWEVEC),WORK(IWEVAL),
     2            WORK(IWRKAR),WORK(IWORK),
     3            WORK(IWRKFA),
     4            WORK(IWRKED),
     5            IRODIM,IROMAX,MAXPSY,INDSYM,
     6            JX2VAL)
C     
C     Kleinsten Eigenwert in E00 speichern
C     
      IF ((JX2VAL.EQ.J2RMIN).AND.(INDSYM.EQ.1)) THEN
        E00=WORK(IWRKED) 
        IF (IPRINT.GT.0) WRITE (NFIL6,6600) E00
      ENDIF
C
C     Berechnung der Zustandsumme PARFUN
C
      IF (IACT.GE.1) THEN
      DO 9001 IENZZ=1,IRODIM
        CALL LABEL(INDSYM,IENZZ,NMIN,
     1             WORK(IASSGN),IWRK(KNFCTS),WORK(IWRKED),
     2             MAXPSY,JMAXP1,LREC35,IROMAX,
     2             NFIX,KVALU,V2VALU,NVALU,ISURF,CABU,ENZZ)
        ENZZ=WORK(IWEVAL+IENZZ-1)
        ENZZ=ENZZ-E00
        ENZZ=(JX2VAL+1)*DEXP(ENZZ*(-1.0D+00)*
     1       PLANCK*VELLGT/BOLTZ/TEMPRA)
        IF ((INDSYM.EQ.1).OR.(INDSYM.EQ.4)) THEN
          ENZZ=ENZZ*IGNS(1,ISO)
        ELSE
          ENZZ=ENZZ*IGNS(2,ISO)
        ENDIF
        IF (ISURF.EQ.ISURF) THEN 
           PARFUN=PARFUN+ENZZ
        ENDIF
9001  CONTINUE
      ENDIF 
C
C     Eigenvektoren fuer Intensitaetsberechnungen speichern
C
      IF (IACT.GE.1) THEN
        IEIGAZ=0
        IEGVEC(INTFLG+1,IPSYM,ITAUP1,1)=IEIGST
1115    IEIGAZ=IEIGAZ+1
        IF (IEIGAZ.GT.IRO(INDSYM)) THEN
          IEIGAZ=IEIGAZ-1
          GOTO 1114 
        ENDIF
        CALL LABEL(INDSYM,IEIGAZ,NMIN,
     1             WORK(IASSGN),IWRK(KNFCTS),WORK(IWRKED),
     2             MAXPSY,JMAXP1,LREC35,IROMAX,
     2             NFIX,KVALU,V2VALU,NVALU,ISURF,CABU,EOUT)
        EOUT=EOUT-E00
        IF (EOUT.GT.ENRGMX(ISO)) THEN 
          IEIGAZ=IEIGAZ-1
          GOTO 1114
        ENDIF
        CALL PACK(JX2VAL,NFIX,KVALU,V2VALU,
     1            NVALU,ISURF,CABU,INDSYM,XLABEL)
        IF (INTFLG.EQ.0) THEN
           NFEIG=NFIL51
        ELSE
           NFEIG=NFIL52
        ENDIF
        WRITE (NFEIG,REC=IEIGST+IEIGAZ)
     1        EOUT,XLABEL,(WORK(IWEVEC+(IEIGAZ-1)*IRODIM+II-1),
     2        II=1,IRODIM)
        GOTO 1115            
1114    IEIGST=IEIGST+IEIGAZ
        IEGVEC(INTFLG+1,IPSYM,ITAUP1,2)=IEIGAZ
      ENDIF
C
C     Store energies on NFIL16 for LSQ fit
C
      KKMAX=MIN0(IRODIM,NREC2)
      DO 1112 KK=1,KKMAX
1112  WORK(KK)=WORK(IWEVAL+KK-1)
      DO 1113 KK=KKMAX+1,NREC2
1113  WORK(KK)=0.0D+00
      IRLSQ=LSQREC(JX2VAL,ITAUP1,IPSYM,MAXPSY)
      WRITE(NFIL16,REC=IRLSQ) (WORK(KK),KK=1,NREC2)
C
8100  CONTINUE
8110  CONTINUE
C
C     Linienstaerken berechnen
C
      IF ((IACT.GE.1).AND.(JX2VAL.NE.0)) THEN
        IINTE1=IINTST
        IINTE2=IINTE1+IEIGMX
        IINTHP=IINTE2+IEIGMX
        IINTDP=IINTHP+MXFCCO
        IF ((IINTDP+MXFCCO*MXFCCO).GT.MWORK) GOTO 1010
        CALL CINTNS(WORK(IINTE1),WORK(IINTE2),
     1              WORK(IINTHP),WORK(IINTDP),
     1              IEGVEC,
     1              NRAUS,IWRK(KNFCTS),
     1              JMAXP1,IEIGMX,MXFCCO,
     2              MAXPSY,NRECUN,
     3              JX2VAL,N2XS,
     4              INTFLG,INTFST,
     5              NFIL6,NFIL53,
     6              NFIL51,NFIL52,
     7              NFIL41,NFIL42,NFIL43)
      ENDIF

C
C     Resultate ausgeben
C
      IF (IPRINT .LE. 0) GOTO 8119
C
      JDIMTO=IRO(1)+IRO(2)+IRO(3)+IRO(4)
C
      IF (MOD(JX2VAL,2) .EQ. 0) THEN
	 JACTU=JX2VAL/2
         WRITE (NFIL6,6502) JACTU
      ELSE
         WRITE (NFIL6,6501) JX2VAL
      ENDIF
C
C
      DO 182 II=1,4
182   LINDEX(II)=0
C
      II=0
      DO 184 JJ=1,4
        IF (IRO(JJ) .NE. 0) THEN
          II=II+1
          LINDEX(II)=JJ
        ENDIF
184   CONTINUE
      NCOL=II
C
      DO 290 II=1,130
290   ELINE(II:II)=' '
      DO 300 II=1,9
300   ELINE(II:II)='-'
      DO 620 JJ=1,NCOL
      IPLLO=10+(JJ-1)*29
      IPLHI=IPLLO+28
      DO 310 II=IPLLO,IPLHI
310   ELINE(II:II)='-'
620   CONTINUE
      WRITE (NFIL6,6300) ELINE
      DO 330 II=1,130
330   ELINE(II:II)=' '
      ELINE(1:1)=':'
      ELINE(9:9)=':'
      DO 640 II=1,NCOL
      JJ=38+(II-1)*29
640   ELINE(JJ:JJ)=':'
      WRITE (NFIL6,6300) ELINE
      DO 350 II=1,NCOL
      JJ=22+(II-1)*29
      KINDEX=LINDEX(II)
      IF (SYMM) THEN
            WRITE (ELINE(JJ:JJ+3),6400) SYMSYM(KINDEX)
      ELSE
            WRITE (ELINE(JJ:JJ+3),6400) ASYSYM(KINDEX)
      ENDIF
350   CONTINUE
      WRITE (NFIL6,6300) ELINE
      DO 660 II=1,130
660   ELINE(II:II)=' '
      ELINE(1:1)=':'
      ELINE(9:9)=':'
      DO 680 II=1,NCOL
      JJ=38+(II-1)*29
680   ELINE(JJ:JJ)=':'
      WRITE (NFIL6,6300) ELINE
      DO 661 II=1,130
661   ELINE(II:II)=' '
      ELINE(1:1)=':'
      ELINE(9:9)=':'
      DO 681 II=1,NCOL
      JJ=38+(II-1)*29
      WRITE (ELINE(JJ-28:JJ-1),6111)
681   ELINE(JJ:JJ)=':'
      WRITE (NFIL6,6300) ELINE
      DO 390 II=1,9
390   ELINE(II:II)='-'
      DO 410 JJ=1,NCOL
      IPLLO=10+(JJ-1)*29
      IPLHI=IPLLO+28
      DO 700 II=IPLLO,IPLHI
700   ELINE(II:II)='-'
410   CONTINUE
      WRITE (NFIL6,6300) ELINE
C
      IDIMMX=IROMAX
C
      DO 412 II=1,130
412   ELINE(II:II)=' '
      DO 190 II=1,IDIMMX
        ELINE(1:1)=':'
        ELINE(9:9)=':'
        KINDEX=II
        WRITE (ELINE(3:6),6000) KINDEX
C
        KOFFS=0
        ENRGOK=.FALSE.
        DO 188 JJ=1,NCOL
          KINDEX=LINDEX(JJ)
          IF (SYMM) THEN
            IF (KINDEX.EQ.2) THEN
              KINDEX=3
            ELSE
              IF (KINDEX.EQ.3) KINDEX=2
            ENDIF
          ENDIF
          IPLLO=10+(JJ-1)*29
          IPLHI=IPLLO+26
          WRITE (ELINE(IPLLO:IPLHI),6200)
          IF (IRO(KINDEX) .GE. II) THEN
               CALL LABEL(KINDEX,II,NMIN,
     1                    WORK(IASSGN),IWRK(KNFCTS),WORK(IWRKED),
     2                    MAXPSY,JMAXP1,LREC35,IROMAX,
     3                    NFIX,KVALU,V2VALU,NVALU,ISURF,CABU,EOUT)
               EOUT=EOUT-E00
               WRITE (ELINE(IPLLO:IPLHI),6100) NFIX,KVALU,V2VALU,
     1                      NVALU,ISURF,CABU,EOUT
               IF (EOUT.LE.ENRGMX(ISO)) ENRGOK=.TRUE.
          ENDIF
          ELINE(IPLHI+2:IPLHI+2)=':'
188     KOFFS=KOFFS+JDIMST(KINDEX)
        IF (ENRGOK) THEN
            WRITE (NFIL6,6300) ELINE
        ELSE
            GOTO 192
        ENDIF
190   CONTINUE
192   DO 720 II=1,9
720   ELINE(II:II)='-'
      DO 740 JJ=1,NCOL
        IPLLO=10+(JJ-1)*29
        IPLHI=IPLLO+28
        DO 430 II=IPLLO,IPLHI
430     ELINE(II:II)='-'
740   CONTINUE
      WRITE (NFIL6,6300) ELINE
8119  CLOSE(NFIL34)
8120  CALL PRTIMJ( JX2VAL , OLDTIM , OLDVEC )
10    CONTINUE
C
C     Uebergaenge sortieren                      
C
      IF (IACT.GE.1) THEN 
        CALL SORT(NRAUS,NFIL53)
      ENDIF
C
C     Intensitaeten berechnen
C
      IF (IACT.GE.1) THEN  
        CALL CALINT(NRAUS,PARFUN,ISO,IACT)
      ENDIF
C
C     Uebergaenge ausgeben
C
      IF (IACT.GE.1) THEN 
        CALL OUTINT(NRAUS,NFIL53,NFIL6,SYMM)
      ENDIF 
C
C*********************************************************************
C
2323  CONTINUE       
C
C     Dateien mit Eigenvektoren schliessen
C
c     IF (IACT.GE.1) THEN 
        CLOSE(NFIL51)
        CLOSE(NFIL52)
c     ENDIF
C
C     Dateien mit Dipolmomentmatrixelementen schliessen
C
      IF (IACT.GE.1) THEN 
        CLOSE(NFIL41)
        CLOSE(NFIL42)
        CLOSE(NFIL43)
      ENDIF
C
C     Datei mit Uebergaengen schliessen
C
      IF (IACT.GE.1) THEN
        CLOSE(NFIL53)
      ENDIF
C
      CLOSE(NFIL35)
      CLOSE(NFIL7)
      CLOSE(NFIL28)
      CLOSE(NFIL37)
      CLOSE(NFIL38)
      CLOSE(NFIL39)
      CLOSE(NFIL40)
C
      RETURN
6000  FORMAT(I4)
6100  FORMAT(4I3,1X,I1,A1,1X,F11.5)
6111  FORMAT('  N KA V2 NS          E     ')
6200  FORMAT(27X)
6300  FORMAT(1H ,1X,A130)
6400  FORMAT(A4)
6500  FORMAT(1H0,12('*'),'  J = ',I3,' ENERGIES  ',12('*'),/)
6501  FORMAT(1H0,12(' '),'*******************************',/,
     1       1H ,12(' '),'****                       ****',/,
     1       1H ,12(' '),'****      J   = ',I3,'/2      ****',/,
     1       1H ,12(' '),'****                       ****',/,
     1       1H ,12(' '),'*******************************')
6502  FORMAT(1H0,12(' '),'*******************************',/,
     1       1H ,12(' '),'****                       ****',/,
     1       1H ,12(' '),'****      J   = ',I3,'        ****',/,
     1       1H ,12(' '),'****                       ****',/,
     1       1H ,12(' '),'*******************************')
6600  FORMAT(1H0,12('*'),'  ZERO POINT ENERGY  ',12('*'),//,
     1      12X,'E0 = ',F12.5//)
6650  FORMAT(1H1,11X,'********************************************',/,
     1           12X,'********************************************',/,
     1           12X,'****                                    ****',/,
     1           12X,'****      PARAMETERS FOR SURFACE ',I1,
     1               '      ****',/,
     1           12X,'****                                    ****',/,
     1           12X,'********************************************',/,
     1           12X,'********************************************')
6651  FORMAT(1H1,11X,'********************************************',/,
     1           12X,'********************************************',/,
     1           12X,'****                                    ****',/,
     1           12X,'****       INTERACTION PARAMETERS       ****',/,
     1           12X,'****                                    ****',/,
     1           12X,'********************************************',/,
     1           12X,'********************************************')
1001  WRITE (NFIL6,2001)
2001  FORMAT(1H0,'RENNER.CNT.ERR  INSUFFICIENT WORK SPACE FOR',
     1          ' INITIAL FUNCTION CALCULATION')
      STOP
1002  WRITE (NFIL6,2002)
2002  FORMAT(1H0,'RENNER.CNT.ERR  INSUFFICIENT WORK SPACE FOR',
     1          ' CALCULATION OF F1, F2 ETC.')
      STOP
1003  WRITE (NFIL6,2003)
2003  FORMAT(1H0,'RENNER.CNT.ERR  INSUFFICIENT WORK SPACE FOR',
     1          ' NUMEROV-COOLEY INTEGRATION')
      STOP
1004  WRITE (NFIL6,2004)
2004  FORMAT(1H0,'RENNER.CNT.ERR  INSUFFICIENT WORK SPACE FOR',
     1          ' BENDING INTEGRAL CALCULATION (DELTA K=0)')
      STOP
1005  WRITE (NFIL6,2005)
2005  FORMAT(1H0,'RENNER.CNT.ERR  INSUFFICIENT WORK SPACE FOR',
     1          ' SETTING UP MORSE OSCILLATOR PARAMETERS')
      STOP
1006  WRITE (NFIL6,2006)
2006  FORMAT(1H0,'RENNER.CNT.ERR  INSUFFICIENT WORK SPACE FOR',
     1          ' DIAGONALIZING THE STRETCHING MATRIX')
      STOP
1007  WRITE (NFIL6,2007)
2007  FORMAT(1H0,'RENNER.CNT.ERR  INSUFFICIENT WORK SPACE FOR',
     1          ' CALCULATING STRETCHING MATRIX ELEMENTS')
      STOP
1009  WRITE (NFIL6,2009) NUNIT
2009  FORMAT(1H0,'RENNER.CNT.ERR  SCRATCH FILE COULD NOT BE O',
     1          'PENED, UNIT NUMBER IS ',I2)
      STOP
1010  WRITE (NFIL6,2010) J
2010  FORMAT(1H0,'RENNER.CNT.ERR  INSUFFICIENT WORK SPACE FOR',
     1          ' ROTATIONAL CALCULATION AT N = ',I3)
      STOP
1012  WRITE (NFIL6,2012)
2012  FORMAT(1H0,'RENNER.CNT.ERR  INSUFFICIENT WORK SPACE FOR',
     1          ' RENNER-TELLER BENDING INTEGRAL CALCULATION')
      STOP
1013  WRITE (NFIL6,2013)
2013  FORMAT(1H0,'RENNER.CNT.ERR  INSUFFICIENT WORK SPACE FOR',
     1          ' CALCULATION OF MORSE OSCILLATOR FUNCTIONS')
      STOP
1014  WRITE (NFIL6,2014)
2014  FORMAT(1H0,'RENNER.CNT.ERR  INSUFFICIENT WORK SPACE FOR',
     1       ' CALCULATION OF R-T STRETCHING MATRIX ELEMENTS')
1015  WRITE (NFIL6,2015)
2015  FORMAT(1H0,'RENNER.CNT.ERR  INSUFFICIENT WORK SPACE FOR',
     1          ' BENDING INTEGRAL CALCULATION (DELTA K=1)')
      STOP
1016  WRITE (NFIL6,2016)
2016  FORMAT(1H0,'RENNER.CNT.ERR  INSUFFICIENT WORK SPACE FOR',
     1          ' BENDING INTEGRAL CALCULATION (DELTA K=2)')
      STOP
1017  WRITE (NFIL6,2017)
2017  FORMAT(1H0,'RENNER.CNT.ERR  INSUFFICIENT WORK SPACE FOR',
     1          ' STORING INTERMEDIATE ENERGIES')
      STOP
      END
C
C
      SUBROUTINE WRMATX ( AWRK , HMAT , IRODIM , LREC34 ,
     1                   J , ITAUP1 , IPSYM )
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 AWRK ( LREC34 ) , HMAT ( IRODIM , IRODIM )
      INTEGER PRTINT
C
      include 'rensys.h'
C
      DO 10 II=1,LREC34
10    AWRK(II)=0.0D+00
      IREC0=(4*J+2*(ITAUP1-1)+IPSYM-1)*LREC34
      DO 30 II=1,IRODIM
      DO 20 KK=1,IRODIM
20    AWRK(KK)=HMAT(KK,II)
      IREC=IREC0+II
30    WRITE (NFIL34,REC=IREC) (AWRK(KK),  KK=1,LREC34)
      RETURN
      END
C
C
      SUBROUTINE RDMATX ( AWRK   , HMAT   , IRODIM , LREC34 ,
     1                    NOFSET , IROBLK , NVAL , ITAUP1 , IPSYM )
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 AWRK ( LREC34 ) , HMAT ( IRODIM , IRODIM )
      INTEGER PRTINT
C
      include 'rensys.h'
C
      IREC0=(4*NVAL+2*(ITAUP1-1)+IPSYM-1)*LREC34
      DO 30 II=1,IROBLK
      IREC=IREC0+II
      READ (NFIL34,REC=IREC) (AWRK(KK),  KK=1,LREC34)
      NOFFSC=NOFSET+II
      DO 20 KK=1,IROBLK
20    HMAT(NOFSET+KK,NOFFSC)=AWRK(KK)
30    CONTINUE
      RETURN
      END
C
C
      INTEGER FUNCTION ILCONV (I,N)
      ILCONV=N*(I-1)+1
      RETURN
      END
C
C
      SUBROUTINE PRTIME( ROUTIN , OLDTIM , OLDVEC )
      IMPLICIT REAL*8 (A-H,O-Z)
      DOUBLE PRECISION OLDTIM,NEWTIM,DELTIM,TIM1S,TIM2S
      CHARACTER*6 ROUTIN
      INTEGER PRTINT
C
      include 'rensys.h'
C
      IF (IPRINT .EQ. 0) RETURN
      CALL CLOCKV ( VECTIM , NEWTIM , 1 , 2 )
      DELTIM = NEWTIM - OLDTIM
      DELVEC = VECTIM - OLDVEC
      TIM1S=DELTIM/1.0D+03
      TIM2S=NEWTIM/1.0D+03
      TIM3S=DELVEC/1.0D+03
      TIM4S=VECTIM/1.0D+03
      WRITE (NFIL6,2007) ROUTIN, TIM1S, TIM3S, TIM2S, TIM4S
2007  FORMAT(1H0,'RENNER.TIM.INF  FOR ROUTINE ',A6,' : '/
     1   1H ,'                CPU TIME           ',F20.8,' SECONDS'/
     1   1H ,'                VECTOR TIME        ',F20.8,' SECONDS'/
     1   1H ,'                TOTAL CPU TIME     ',F20.8,' SECONDS'/
     1   1H ,'                TOTAL VECTOR TIME  ',F20.8,' SECONDS')
      OLDTIM=NEWTIM
      OLDVEC=VECTIM
      RETURN
      END
C
C
      SUBROUTINE PRTIMJ( JX2VAL , OLDTIM , OLDVEC )
      IMPLICIT REAL*8 (A-H,O-Z)
      DOUBLE PRECISION OLDTIM,NEWTIM,DELTIM,TIM1S,TIM2S
      INTEGER PRTINT
C
      include 'rensys.h'
C
      IF (IPRINT .EQ. 0) RETURN
      CALL CLOCKV ( VECTIM , NEWTIM , 1 , 2 )
      DELTIM = NEWTIM - OLDTIM
      DELVEC = VECTIM - OLDVEC
      TIM1S=DELTIM/1.0D+03
      TIM2S=NEWTIM/1.0D+03
      TIM3S=DELVEC/1.0D+03
      TIM4S=VECTIM/1.0D+03
C
      IF (MOD(JX2VAL,2) .EQ. 0) THEN
	 JACTU=JX2VAL/2
         WRITE (NFIL6,2007) JACTU, TIM1S, TIM3S, TIM2S, TIM4S
      ELSE
         WRITE (NFIL6,2008) JX2VAL, TIM1S, TIM3S, TIM2S, TIM4S
      ENDIF
C
2007  FORMAT(1H0,'RENNER.TIM.INF  FOR CALCULATING J = ',I3,' : '/
     1   1H ,'                CPU TIME           ',F20.8,' SECONDS'/
     1   1H ,'                VECTOR TIME        ',F20.8,' SECONDS'/
     1   1H ,'                TOTAL CPU TIME     ',F20.8,' SECONDS'/
     1   1H ,'                TOTAL VECTOR TIME  ',F20.8,' SECONDS')
2008  FORMAT(1H0,'RENNER.TIM.INF  FOR CALCULATING J = ',I3,'/2 : '/
     1   1H ,'                CPU TIME           ',F20.8,' SECONDS'/
     1   1H ,'                VECTOR TIME        ',F20.8,' SECONDS'/
     1   1H ,'                TOTAL CPU TIME     ',F20.8,' SECONDS'/
     1   1H ,'                TOTAL VECTOR TIME  ',F20.8,' SECONDS')
      OLDTIM=NEWTIM
      OLDVEC=VECTIM
      RETURN
      END
C
C
      SUBROUTINE PRTIMN( J , OLDTIM , OLDVEC )
      IMPLICIT REAL*8 (A-H,O-Z)
      DOUBLE PRECISION OLDTIM,NEWTIM,DELTIM,TIM1S,TIM2S
      INTEGER PRTINT
C
      include 'rensys.h'
C
      IF (IPRINT .EQ. 0) RETURN
      CALL CLOCKV ( VECTIM , NEWTIM , 1 , 2 )
      DELTIM = NEWTIM - OLDTIM
      DELVEC = VECTIM - OLDVEC
      TIM1S=DELTIM/1.0D+03
      TIM2S=NEWTIM/1.0D+03
      TIM3S=DELVEC/1.0D+03
      TIM4S=VECTIM/1.0D+03
      WRITE (NFIL6,2007) J, TIM1S, TIM3S, TIM2S, TIM4S
2007  FORMAT(1H0,'RENNER.TIM.INF  FOR CALCULATING N = ',I3,
     1   ' MATRICES : '/
     1   1H ,'                CPU TIME           ',F20.8,' SECONDS'/
     1   1H ,'                VECTOR TIME        ',F20.8,' SECONDS'/
     1   1H ,'                TOTAL CPU TIME     ',F20.8,' SECONDS'/
     1   1H ,'                TOTAL VECTOR TIME  ',F20.8,' SECONDS')
      OLDTIM=NEWTIM
      OLDVEC=VECTIM
      RETURN
      END
C
C
      SUBROUTINE PRJACO ( IPRM , PXNAME )
      IMPLICIT REAL*8 (A-H,O-Z)
      DOUBLE PRECISION NEWTIM,TIM1S
      CHARACTER*6 ROUTIN,PXNAME
      INTEGER PRTINT
C
      include 'rensys.h'
C
      CALL CLOCKV ( VECTIM , NEWTIM , 1 , 2 )
      TIM2S=NEWTIM/1.0D+03
      TIM4S=VECTIM/1.0D+03
      WRITE (NFIL6,2007) IPRM, PXNAME, TIM2S, TIM4S
2395  FORMAT('  JACOBIAN CALCULATION COMPLETED FOR PARAMETER NO.',I2,
     1       ', ',A8)
2007  FORMAT(1H0,'RENNER.CNT.INF  JACOBIAN CALCULATION COMPLETED FOR'/
     1       1H ,'                PARAMETER NO.',I3,', ',A8/
     1   1H ,'                TOTAL CPU TIME     ',F20.8,' SECONDS'/
     1   1H ,'                TOTAL VECTOR TIME  ',F20.8,' SECONDS')
      RETURN
      END
C
C
      SUBROUTINE PRTITO
      IMPLICIT REAL*8 (A-H,O-Z)
      DOUBLE PRECISION NEWTIM,TIM1S
      CHARACTER*6 ROUTIN
      INTEGER PRTINT
C
      include 'rensys.h'
C
      CALL CLOCKV ( VECTIM , NEWTIM , 1 , 2 )
      TIM2S=NEWTIM/1.0D+03
      TIM4S=VECTIM/1.0D+03
      WRITE (NFIL6,2007) TIM2S, TIM4S
2007  FORMAT(1H0,'RENNER.CNT.INF  PROGRAM TERMINATED NORMALLY ',/,
     1       1H ,'                ENTIRE CALCULATION USED : '/
     1   1H ,'                TOTAL CPU TIME     ',F20.8,' SECONDS'/
     1   1H ,'                TOTAL VECTOR TIME  ',F20.8,' SECONDS')
      RETURN
      END
C***************************************************************
      SUBROUTINE LABEL(INDSYM,IEIGNR,NMIN,
     1                 ASSGN,NFCTS,EIGMAX,
     2                 MAXPSY,JMAXP1,LREC35,IROMAX,
     3                 NFIX,KVALU,V2VALU,NVALU,ISURF,CABU,EOUT)
C***************************************************************
      REAL*8 ASSGN(JMAXP1,2*MAXPSY,LREC35,5)
      INTEGER NFCTS(JMAXP1,4)
      REAL*8  EIGMAX(2*MAXPSY,IROMAX,2)
C
      INTEGER INDSYM,IEIGNR,NMIN
C
      INTEGER MAXPSY,JMAXP1,LREC35,IROMAX
C
      INTEGER NFIX,KVALU,V2VALU,NVALU,ISURF
      INTEGER ICABU
      CHARACTER*1 CABU
      REAL*8  EOUT
C
      INTEGER MAXCOE
      INTEGER AKTPOS,N,K
      INTEGER AKTALT,NALT,KALT
C
      EOUT  =      EIGMAX(INDSYM,IEIGNR,1)
      MAXCOE=IDINT(EIGMAX(INDSYM,IEIGNR,2))
      AKTPOS=1
      N=NMIN
      K=0
C
10    NALT=N
      KALT=K
      AKTALT=AKTPOS
C
      IF ((K.EQ.0).AND.(MOD(N,2).EQ.1)) THEN
        IF (MOD(INDSYM,2).EQ.1) THEN
          INDSYH=INDSYM+1
        ELSE
          INDSYH=INDSYM-1
        ENDIF
      ELSE
        INDSYH=INDSYM
      ENDIF
      AKTPOS=AKTPOS+NFCTS(K+1,INDSYH)
      K=K+1
      IF (K.GT.N) THEN
        K=0
        N=N+1
      ENDIF
      IF (AKTPOS.LE.MAXCOE) GOTO 10
C
      MAXCOE=MAXCOE-AKTALT+1
C
      NFIX=NALT
      KVALU=KALT
      IF ((KVALU.EQ.0).AND.(MOD(NFIX,2).EQ.1)) THEN
        IF (MOD(INDSYM,2).EQ.1) THEN
          INDSYH=INDSYM+1
        ELSE
          INDSYH=INDSYM-1
        ENDIF
      ELSE
        INDSYH=INDSYM
      ENDIF
      KVP1=KVALU+1
      ISURF=NINT(ASSGN(KVP1,INDSYH,MAXCOE,1))
      NVALU=NINT(ASSGN(KVP1,INDSYH,MAXCOE,2))
      V2VALU=NINT(ASSGN(KVP1,INDSYH,MAXCOE,3))
      ICABU=NINT(ASSGN(KVP1,INDSYH,MAXCOE,4))
      IF (ICABU.EQ.1) THEN
        CABU='A'
      ELSE
        CABU='B'
      ENDIF
C
      RETURN
      END
C
      SUBROUTINE BABLK(NABLKR,NABLKC,
     1                 NBBLKR,NBBLKC,
     2                 KDEL,
     3                 HMAT,
     4                 KMAT,
     5                 LREC34,LREC35)
C
      REAL*8 KMAT(LREC34,LREC34)
      REAL*8 HMAT(LREC35,LREC35)
C
      IF (KDEL.EQ.0) THEN
        DO 10 II=NABLKC+1,NABLKC+NBBLKC
        DO 10 JJ=1,NABLKR
          HMAT(II,JJ)=HMAT(JJ,II)
10      CONTINUE
      ELSE
        DO 20 II=1,NBBLKR
        DO 20 JJ=1,NABLKC
          HMAT(NABLKR+II,JJ)=KMAT(JJ,II)
20      CONTINUE
      ENDIF
C
      RETURN
      END
C
      SUBROUTINE SCHREI(IREC,
     1                  HMAT,
     2                  LREC35,MXFCCO,
     3                  NFIL34)
C
      REAL*8 HMAT(LREC35,LREC35)
C
      WRITE (NFIL34,REC=IREC)
     1      ((HMAT(JJ,II),JJ=1,MXFCCO),II=1,MXFCCO)
C
      RETURN
      END
C
C
      INTEGER FUNCTION LSQREC ( JX2 , ITAUP1 , IPSYM , MAXPSY )
      IRST=MOD(JX2,2)
      IXJ=(JX2-IRST)/2
      LSQREC=2*MAXPSY*IXJ+MAXPSY*(ITAUP1-1)+IPSYM
      RETURN
      END
C
C
      INTEGER FUNCTION JINDEX ( JX2 )
      IRST=MOD(JX2,2)
      JINDEX=(JX2-IRST)/2+1
      RETURN
      END
C
      SUBROUTINE READFI(NFILNR,ARR,LENGTH)
C
      INTEGER NFILNR
      INTEGER LENGTH
      REAL*8  ARR(LENGTH)
C
      READ (NFILNR) (ARR(II),II=1,LENGTH)
C
      RETURN
      END
C
      SUBROUTINE BLKCOP(TOTM,TOTR,TOTC,
     1                     M,   R,   C,
     2                  ROFF,COFF,
     3                  TRANS)
C
      INTEGER TOTR,TOTC,
     1           R,   C,
     1        ROFF,COFF
      REAL*8 TOTM(TOTR,TOTC),
     1          M(   R,   C)
      LOGICAL TRANS
C
      DO 10 II=1,R
      DO 20 JJ=1,C
         IF (TRANS) THEN
            TOTM(ROFF+JJ,COFF+II)=M(II,JJ)
         ELSE
            TOTM(ROFF+II,COFF+JJ)=M(II,JJ)
         ENDIF
20    CONTINUE
10    CONTINUE
C     
      RETURN 
      END
C
      SUBROUTINE WRTDIP(FELD,LREC35,MXFCCO,NFDP,NRCDP)
C
      INTEGER LREC35 
      INTEGER MXFCCO,NFDP,NRCDP
      REAL*8  FELD(LREC35,LREC35)
C
      WRITE (NFDP,REC=NRCDP)
     1      ((FELD(II,JJ),II=1,MXFCCO),JJ=1,MXFCCO)
C
      RETURN
      END
C
      SUBROUTINE RCNRDP(IFLAGR,
     1                  NPARR,NPARC,
     1                  K1P1,K2P1,KDEL,
     1                  NRCNR)
C
      INTEGER IFLAGR
      INTEGER NPARR,NPARC
      INTEGER K1P1,K2P1,KDEL 
      INTEGER NRCNR
C
      IF (KDEL.EQ.0) THEN 
        IF ((K1P1-1).EQ.0) THEN
          NRCNR=(NPARR*2+NPARC)
     1         +4*(IFLAGR-1)+1
        ELSE
          NRCNR=(K1P1-1)*2
     1         +(IFLAGR-1)
     1         +9
        ENDIF 
      ENDIF
C
      IF (KDEL.EQ.1) THEN 
        IF ((K1P1-1).EQ.0) THEN
          NRCNR=NPARR
     1         +2*(IFLAGR-1)+1
        ELSE
          NRCNR=(K1P1-1)*2
     1         +(IFLAGR-1)
     1         +5
        ENDIF
      ENDIF
C
      IF (KDEL.EQ.(-1)) THEN 
        IF ((K2P1-1).EQ.0) THEN
          NRCNR=NPARC
     1         +2*(IFLAGR-1)+1
        ELSE
          NRCNR=(K2P1-1)*2
     1         +(IFLAGR-1)
     1         +5
        ENDIF
      ENDIF
C
      RETURN
      END
C*************************************
      SUBROUTINE DPTRA1(DIPMAT,LREC35,
     1                 NROW,NCOL,
     2                 NFIL6,
     3                 INDSYR,INDSYC,
     4                 K1P1,KDEL,
     5                 NPARR,NPARC)
C*************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      REAL*8 DIPMAT(LREC35,LREC35)
      INTEGER LREC35
      INTEGER NROW,NCOL
      INTEGER NFIL6
      INTEGER INDSYR,INDSYC
      INTEGER K1P1,KDEL
      INTEGER NPARR,NPARC
C
      INTEGER NRC
      INTEGER K2P1
      CHARACTER*2 SYR,SYC
      INTEGER NR,NCU,NCL
C
      NRC=10
      WRITE (NFIL6,*)
     1      '***** Transformed Dipole Moment Matrix ****'
      IF (INDSYR.EQ.1) SYR='A1'
      IF (INDSYR.EQ.3) SYR='B2'
      IF (INDSYR.EQ.2) SYR='B1'
      IF (INDSYR.EQ.4) SYR='A2'
      IF (INDSYC.EQ.1) SYC='A1'
      IF (INDSYC.EQ.3) SYC='B2'
      IF (INDSYC.EQ.2) SYC='B1'
      IF (INDSYC.EQ.4) SYC='A2'
      K2P1=K1P1+KDEL
      WRITE (NFIL6,1000) SYR,SYC,K1P1-1,K2P1-1 
      IF (K1P1.EQ.1) THEN
        IF (NPARR.EQ.0) THEN
           WRITE(NFIL6,*) '   N_ROW even   '
        ELSE
           WRITE(NFIL6,*) '   N_ROW odd    '
        ENDIF
      ENDIF
      IF (K2P1.EQ.1) THEN
        IF (NPARC.EQ.0) THEN
           WRITE(NFIL6,*) '   N_COL even   '
        ELSE
           WRITE(NFIL6,*) '   N_COL odd    '
        ENDIF
      ENDIF
      NCL=1 
2000  NCU=NCL+NRC-1
      IF (NCU.GT.NCOL) NCU=NCOL
      IF (NCL.NE.1) WRITE(NFIL6,*)
      WRITE(NFIL6,1010) (J,J=NCL,NCU)
      DO 2020 NR=1,NROW
         WRITE(NFIL6,1020) NR,(DIPMAT(NR,J),J=NCL,NCU)
2020  CONTINUE     
      NCL=NCL+NRC
      IF (NCU.NE.NCOL) GOTO 2000    
C
1000  FORMAT(' SYM_ROW=',A2,' SYM_COL=',A2,' K_ROW=',I4,' K_COL=',I4,/) 
1010  FORMAT('     ',10I10)
1020  FORMAT(I5,10F10.4)
C
      RETURN
      END

C*************************************
      SUBROUTINE DPTRA2(DIPMAT,LREC35,
     1                 ASSGN,JMAXP1,MAXPSY,
     1                 NROW,NCOL,
     2                 NFIL6,
     3                 INDSYR,INDSYC,
     3                 NABSR,NABSC,
     4                 K1P1,KDEL,
     5                 NPARR,NPARC,
     6                 DPFL)
C*************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      REAL*8 ASSGN(JMAXP1,2*MAXPSY,LREC35,5)
      REAL*8 DIPMAT(LREC35,LREC35)
      INTEGER LREC35
      INTEGER NROW,NCOL
      INTEGER NFIL6
      INTEGER INDSYR,INDSYC
      INTEGER NABSR,NABSC
      INTEGER K1P1,KDEL
      INTEGER NPARR,NPARC
      LOGICAL DPFL
C
      INTEGER K2P1
      CHARACTER*2 SYR,SYC
      CHARACTER*2 STR,STC
      INTEGER FILSUR,FILNS,FILST
      INTEGER IVEC
      REAL*8  EDIFF,EL,E00,DP
C
      FILSUR=1
      FILNS=1
      FILV2=0
      FILST=1
      E00=ASSGN(1,1,1,5)
C
      IF (INDSYR.EQ.1) SYR='A1'
      IF (INDSYR.EQ.3) SYR='B2'
      IF (INDSYR.EQ.2) SYR='B1'
      IF (INDSYR.EQ.4) SYR='A2'
C
      IF (INDSYC.EQ.1) SYC='A1'
      IF (INDSYC.EQ.3) SYC='B2'
      IF (INDSYC.EQ.2) SYC='B1'
      IF (INDSYC.EQ.4) SYC='A2'
C
      IF (NABSR.EQ.1) THEN
         STR='A1'
      ELSE
         STR='B2'
      ENDIF
      IF (NABSC.EQ.1) THEN
         STC='A1'
      ELSE
         STC='B2'
      ENDIF
      IF (MAXPSY.EQ.1) THEN
         STC='A'''
         STR='A'''
      ENDIF
C
      K2P1=K1P1+KDEL
C
      IF ((K1P1.EQ.1).AND.(NPARR.EQ.1)) THEN
        IF (MOD(INDSYR,2).EQ.1) THEN
          INDSHR=INDSYR+1
        ELSE
          INDSHR=INDSYR-1
        ENDIF
      ELSE
        INDSHR=INDSYR
      ENDIF
      IF ((K2P1.EQ.1).AND.(NPARC.EQ.1)) THEN
        IF (MOD(INDSYC,2).EQ.1) THEN
          INDSHC=INDSYC+1
        ELSE
          INDSHC=INDSYC-1
        ENDIF
      ELSE
        INDSHC=INDSYC
      ENDIF
C
      IF (DPFL) THEN
      WRITE(NFIL6,1000)
      DPFL=.FALSE.
      ENDIF
C
      IF ((NINT(ASSGN(K1P1,INDSHR,1,1)).EQ.FILSUR).AND.
     1    (NINT(ASSGN(K1P1,INDSHR,1,2)).EQ.FILNS ).AND.
     1    (NINT(ASSGN(K1P1,INDSHR,1,3)).EQ.FILV2 ).AND.
     1    (NINT(ASSGN(K1P1,INDSHR,1,4)).EQ.FILST )) THEN
      EL=ASSGN(K1P1,INDSHR,1,5)-E00
      DO 100 IVEC=1,NCOL,1
      DP=DIPMAT(1,IVEC)
      DP=DP*DP
      EDIFF=ASSGN(K2P1,INDSHC,IVEC,5)-E00-EL
      IF ((EDIFF.GT.0.0D+00).AND.(DP.GT.1.0D-5)) THEN
      WRITE (NFIL6,1010)
     1      K1P1-1,K2P1-1,
     2      NINT(ASSGN(K2P1,INDSHC,IVEC,1)),
     2      NINT(ASSGN(K2P1,INDSHC,IVEC,2)),
     2      NINT(ASSGN(K2P1,INDSHC,IVEC,3)),
     2      STC,
     2      EL,EDIFF,DP
      ENDIF
100   CONTINUE
      WRITE(NFIL6,*)
      ENDIF
C
1000  FORMAT(10X,'   KL','   KU',' SURU','  NSU','  V2U','  STU',
     1       '            EL','            DE','            DP')
1010  FORMAT(10X,5I5,2X,A2,1X,2F14.4,F14.8)
C
      RETURN
      END
C*****************************************
      SUBROUTINE AUS(FELD,LREC35,NFIL6)
C*****************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      REAL*8 FELD(LREC35,LREC35)
      INTEGER LREC35
      INTEGER NFIL6
C
      INTEGER NRC
      INTEGER NR,NCU,NCL
      INTEGER NROW,NCOL
C
      NROW=LREC35
      NCOL=LREC35
      NRC=10
      NCL=1 
2000  NCU=NCL+NRC-1
      IF (NCU.GT.NCOL) NCU=NCOL
      IF (NCL.NE.1) WRITE(NFIL6,*)
      WRITE(NFIL6,1010) (J,J=NCL,NCU)
      DO 2020 NR=1,NROW
         WRITE(NFIL6,1020) NR,(FELD(NR,J),J=NCL,NCU)
2020  CONTINUE     
      NCL=NCL+NRC
      IF (NCU.NE.NCOL) GOTO 2000    
C
1010  FORMAT('     ',10I15)
1020  FORMAT(I5,10F15.4)
C
      RETURN
      END
