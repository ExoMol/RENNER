      DOUBLE PRECISION FUNCTION TRIALE(K,J,ENERGY)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 M1,M2,M3,M,U1,U3,U13,V
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
      REAL*8 B11,B13,B111,B133,B113,ETEST(0:44,0:1),
     2      B1111,B1333,B1113,B1133,
     3      B11111,B13333,B11113,B11333,B11133,
     4      B31,B33,B311,B333,B313,
     5      B3111,B3333,B3113,B3133,
     6      B31111,B33333,B31113,B31333,B31133
C
      REAL*8 CR1,CR3,CR11,CR33,CR13,
     2      CR111,CR333,CR113,CR133,
     3      CR1111,CR3333,CR1113,CR1333,CR1133
C
      INTEGER NFIL1 , NFIL2 , NFIL3 , NFIL4 , NFIL5 ,
     5       NFIL6 , NFIL7 , NFIL8 , NFIL9 , NFIL10 ,
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
     4       KSTYPA(2) , LSTYPA(2) , JMAX , V2MAX , JMAXP1
C
      INTEGER IQUANT(9,10)
C
      LOGICAL SYMM
C
      REAL*8 RMK(2),AAS(2)
      INTEGER NOBAS,MBASP1,LENIW
C
C     FUNCTION GUESSES A NEW STARTING VALUE FOR EGUESS FOR USE IN THE
C     NUMEROV - COOLEY INTEGRATION PROCEDURE.
C     HOW IT WORKS I DON'T KNOW.  IT IS TAKEN FROM SYMMETRIC NRB
C     OF AITKEN HOY.
C
      REAL*8 ENERGY(V2MXP1,JMAXP1)
      INTEGER K,J
C
      include 'isotop.h'
      include 'value.h'
      include 'lzcomp.h'
      COMMON /MOLCUL/ RHOE,FA1,FA2,FA3,FA4,FA5,FA6,FA7,FA8,
     1               AA1,AA3,
     1               F11,F33,F13,F111,F333,F113,F133,
     1               F1111,F3333,F1113,F1333,F1133,
     2               Z1,F1A1,F2A1,F3A1,F4A1,F3,F1A3,F2A3,F3A3,F4A3,
     3               F1A11,F2A11,F3A11,F1A33,F2A33,F3A33,
     4               F1A13,F2A13,F3A13,
     5               F1A111,F2A111,F1A333,F2A333,
     6               F1A113,F2A113,F1A133,F2A133,
     7               FA1111,FA3333,FA1113,
     8               FA1333,FA1133, R12RF1, R32RF1, R12RF2, R32RF2,
     8               RE12 , RE32 , M1 , M2 , M3 , M ,
     9               U1 , U3 , U13 , V ,
     1               SYMM
      include 'rentel.h'
      include 'integ.h'
      include 'dimen.h'
      include 'lsfit.h'
      include 'rensys.h'
      include 'bcoeff.h'
      include 'crcoef.h'
      include 'morse.h'
      include 'modim.h'
C
       IF (NSURF .EQ. 1) THEN
       ETEST (   0 , 0    ) =      536.90519   
       ETEST (   0 , 1    ) =      594.85546  
       ETEST (   1 , 0    ) =     1550.53839 
       ETEST (   1 , 1    ) =     1445.69779  
       ETEST (   2 , 0    ) =     2681.51800  
       ETEST (   2 , 1    ) =     2316.23781  
       ETEST (   3 , 0    ) =     3964.43838  
       ETEST (   3 , 1    ) =     3470.78008  
       ETEST (   4 , 0    ) =     5377.30395  
       ETEST (   4 , 1    ) =     4808.56077  
       ETEST (   5 , 0    ) =     6896.05736  
       ETEST (   5 , 1    ) =     6267.17682  
       ETEST (   6 , 0    ) =     8503.46297  
       ETEST (   6 , 1    ) =     7822.95632  
       ETEST (   7 , 0    ) =    10187.57178  
       ETEST (   7 , 1    ) =     9462.19993  
       ETEST (   8 , 0    ) =    11940.06386  
       ETEST (   8 , 1    ) =    11175.11497  
       ETEST (   9 , 0    ) =    13755.11403  
       ETEST (   9 , 1    ) =    12954.48465  
       ETEST (  10 , 0    ) =    15628.60265  
       ETEST (  10 , 1    ) =    14795.08072 
       ETEST (  11 , 0    ) =    17557.54906 
       ETEST (  11 , 1    ) =    16693.14891 
       ETEST (  12 , 0    ) =    19539.70376 
       ETEST (  12 , 1    ) =    18645.94553 
       ETEST (  13 , 0    ) =    21573.26147 
       ETEST (  13 , 1    ) =    20651.36752 
       ETEST (  14 , 0    ) =    23656.66639 
       ETEST (  14 , 1    ) =    22707.68620 
       ETEST (  15 , 0    ) =    25788.48592 
       ETEST (  15 , 1    ) =    24813.37077 
       ETEST (  16 , 0    ) =    27967.33247 
       ETEST (  16 , 1    ) =    26966.97889 
       ETEST (  17 , 0    ) =    30191.81775 
       ETEST (  17 , 1    ) =    29167.09287 
       ETEST (  18 , 0    ) =    32460.52767 
       ETEST (  18 , 1    ) =    31412.28449 
       ETEST (  19 , 0    ) =    34772.00974 
       ETEST (  19 , 1    ) =    33701.09679 
       ETEST (  20 , 0    ) =    37124.76756 
       ETEST (  20 , 1    ) =    36032.03533 
       ETEST (  21 , 0    ) =    39517.25889 
       ETEST (  21 , 1    ) =    38403.56432 
       ETEST (  22 , 0    ) =    41947.89526 
       ETEST (  22 , 1    ) =    40814.10503 
       ETEST (  23 , 0    ) =    44415.04183 
       ETEST (  23 , 1    ) =    43262.03496 
       ETEST (  24 , 0    ) =    46917.01671 
       ETEST (  24 , 1    ) =    45745.68691 
       ETEST (  25 , 0    ) =    49452.08940 
       ETEST (  25 , 1    ) =    48263.34736 
       ETEST (  26 , 0    ) =    52018.47798 
       ETEST (  26 , 1    ) =    50813.25402 
       ETEST (  27 , 0    ) =    54614.34498 
       ETEST (  27 , 1    ) =    53393.59203 
       ETEST (  28 , 0    ) =    57237.79157 
       ETEST (  28 , 1    ) =    56002.48890 
       ETEST (  29 , 0    ) =    59886.85023 
       ETEST (  29 , 1    ) =    58638.00774 
       ETEST (  30 , 0    ) =    62559.47578 
       ETEST (  30 , 1    ) =    61298.13896 
       ETEST (  31 , 0    ) =    65253.53697 
       ETEST (  31 , 1    ) =    63980.79105 
       ETEST (  32 , 0    ) =    67966.81746 
c      ETEST (  32 , 1    ) =    66683.78503 
       ETEST (  32 , 1    ) =    66681.78503 
       ETEST (  33 , 0    ) =    70697.06082 
       ETEST (  33 , 1    ) =    69404.87099 
       ETEST (  34 , 0    ) =    73442.16946 
       ETEST (  34 , 1    ) =    72141.83168 
       ETEST (  35 , 0    ) =    76200.83773 
       ETEST (  35 , 1    ) =    74892.85818 
       ETEST (  36 , 0    ) =    78974.10997 
       ETEST (  36 , 1    ) =    77657.59747 
       ETEST (  37 , 0    ) =    81768.12925 
       ETEST (  37 , 1    ) =    80439.35493 
       ETEST (  38 , 0    ) =    84596.78676 
       ETEST (  38 , 1    ) =    83248.10468 
       ETEST (  39 , 0    ) =    87480.89172 
       ETEST (  39 , 1    ) =    86101.77827 
       ETEST (  40 , 0    ) =    90442.22845 
       ETEST (  40 , 1    ) =    89022.61679 
       ETEST (  41 , 0    ) =    93496.77017 
       ETEST (  41 , 1    ) =    92029.98105 
       ETEST (  42 , 0    ) =    96652.37619 
       ETEST (  42 , 1    ) =    95135.54511 
       ETEST (  43 , 0    ) =    99910.77953 
       ETEST (  43 , 1    ) =    98343.59077 
       ETEST (  44 , 0    ) =   103270.55920
       ETEST (  44 , 1    ) =   101653.88205 
       else
       ETEST (   0 , 0    ) =     2428.35818   
       ETEST (   0 , 1    ) =     3323.12082   
       ETEST (   1 , 0    ) =     4866.64269   
       ETEST (   1 , 1    ) =     5855.18207   
       ETEST (   2 , 0    ) =     7314.28845   
       ETEST (   2 , 1    ) =     8350.58590   
       ETEST (   3 , 0    ) =     9772.51991   
       ETEST (   3 , 1    ) =    10838.66357   
       ETEST (   4 , 0    ) =    12241.44286   
       ETEST (   4 , 1    ) =    13329.29591   
       ETEST (   5 , 0    ) =    14720.76426   
       ETEST (   5 , 1    ) =    15825.96392   
       ETEST (   6 , 0    ) =    17210.19120   
       ETEST (   6 , 1    ) =    18330.05657   
       ETEST (   7 , 0    ) =    19709.62354   
       ETEST (   7 , 1    ) =    20842.39129   
       ETEST (   8 , 0    ) =    22219.21786   
       ETEST (   8 , 1    ) =    23363.71905   
       ETEST (   9 , 0    ) =    24739.37343   
       ETEST (   9 , 1    ) =    25894.85494   
       ETEST (  10 , 0    ) =    27270.67464   
       ETEST (  10 , 1    ) =    28436.67311   
       ETEST (  11 , 0    ) =    29813.81383   
       ETEST (  11 , 1    ) =    30990.05646   
       ETEST (  12 , 0    ) =    32369.51119   
       ETEST (  12 , 1    ) =    33555.83689   
       ETEST (  13 , 0    ) =    34938.44248   
       ETEST (  13 , 1    ) =    36134.74119   
       ETEST (  14 , 0    ) =    37521.18067   
       ETEST (  14 , 1    ) =    38727.34874  
       ETEST (  15 , 0    ) =    40118.15331   
       ETEST (  15 , 1    ) =    41334.06244   
       ETEST (  16 , 0    ) =    42729.61519   
       ETEST (  16 , 1    ) =    43955.09218   
       ETEST (  17 , 0    ) =    45355.63390   
       ETEST (  17 , 1    ) =    46590.44847   
       ETEST (  18 , 0    ) =    47996.08525   
       ETEST (  18 , 1    ) =    49239.94371   
       ETEST (  19 , 0    ) =    50650.65567   
       ETEST (  19 , 1    ) =    51903.19837   
       ETEST (  20 , 0    ) =    53318.84884   
       ETEST (  20 , 1    ) =    54579.65000   
       ETEST (  21 , 0    ) =    55999.99439   
       ETEST (  21 , 1    ) =    57268.56322   
       ETEST (  22 , 0    ) =    58693.25714   
       ETEST (  22 , 1    ) =    59969.03929   
       ETEST (  23 , 0    ) =    61397.64550   
       ETEST (  23 , 1    ) =    62680.02450   
       ETEST (  24 , 0    ) =    64112.01843   
       ETEST (  24 , 1    ) =    65400.31648   
       ETEST (  25 , 0    ) =    66835.09021   
       ETEST (  25 , 1    ) =    68128.56818   
       ETEST (  26 , 0    ) =    69565.43263   
       ETEST (  26 , 1    ) =    70863.28889   
       ETEST (  27 , 0    ) =    72301.47423   
       ETEST (  27 , 1    ) =    73602.84220   
       ETEST (  28 , 0    ) =    75041.49614   
       ETEST (  28 , 1    ) =    76345.44021   
       ETEST (  29 , 0    ) =    77783.62417   
       ETEST (  29 , 1    ) =    79089.13366   
       ETEST (  30 , 0    ) =    80525.81623   
       ETEST (  30 , 1    ) =    81831.79708   
       ETEST (  31 , 0    ) =    83265.84474   
       ETEST (  31 , 1    ) =    84571.10861   
       ETEST (  32 , 0    ) =    86001.27383   
       ETEST (  32 , 1    ) =    87304.52661   
       ETEST (  33 , 0    ) =    88729.44001   
       ETEST (  33 , 1    ) =    90029.28456   
       ETEST (  34 , 0    ) =    91447.49015   
       ETEST (  34 , 1    ) =    92742.52014   
       ETEST (  35 , 0    ) =    94152.72613   
       ETEST (  35 , 1    ) =    95442.00371   
       ETEST (  36 , 0    ) =    96844.09523   
       ETEST (  36 , 1    ) =    98128.74648   
       ETEST (  37 , 0    ) =    99526.50771   
       ETEST (  37 , 1    ) =   100813.05299   
       ETEST (  38 , 0    ) =   102218.13667   
       ETEST (  38 , 1    ) =   103521.21595   
       ETEST (  39 , 0    ) =   104953.27298   
       ETEST (  39 , 1    ) =   106292.04344   
       ETEST (  40 , 0    ) =   107770.83622   
       ETEST (  40 , 1    ) =   109160.21061   
       ETEST (  41 , 0    ) =   110697.36069   
       ETEST (  41 , 1    ) =   112143.81186   
       ETEST (  42 , 0    ) =   113742.43265   
       ETEST (  42 , 1    ) =   115246.42263   
       ETEST (  43 , 0    ) =   116905.15455   
       ETEST (  43 , 1    ) =   118464.72177   
       ETEST (  44 , 0    ) =   120180.87728   
       ETEST (  44 , 1    ) =   121793.61236   
      endif
      if (k .le. 2) then
          triale = etest ( j-1 , k-1 )
      else
          triale = etest ( j-1 ,   1 )
      endif
      RETURN
C
      END
