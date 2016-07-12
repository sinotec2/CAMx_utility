C     COMDECK.CHPARM
C-----------------------------------------------------------------------
C
C     ***  /CHPARM/ CONTAINS CHEMISTRY PARAMETERS
C
      COMMON /CHPARM/ NSMAX,
     $  MSPEC(10,32), LREAC(32), LSSIC(32), LSSBC(32),
     $     BDSL (32), BDSU (32), BDNL (32), BDNU (32), DVRES(32)
C
      COMMON /CHPARM/ NRMAX,
     $  RKN  (90), RK   (90), LPHOT(90), LTDEP(90),
     $  ACTEN(90), TREF (90), TRINV(90), LATDEP(40)
C
      COMMON /CHPARM/ NCMAX,
     $  MCOEF(10,10), CVAL(10)
C
      LOGICAL  LREAC, LSSIC, LSSBC
      LOGICAL LPHOT, LTDEP
C
