MODULE physvars
!
! physical parameters 
!
      REAL,    PARAMETER :: BE    = 13.6056923      ! [eV] binding energy H 
      REAL,    PARAMETER :: BEKEV = BE*1.E-3        ! [keV] 
      REAL,    PARAMETER :: CC    = 2.9979E10       ! [cm/s] speed of light 
      REAL,    PARAMETER :: CC2   = CC*CC           ! [cm^2/s^2] CC squared 
      REAL,    PARAMETER :: CC3   = CC*CC*CC        ! [cm^2/s^2] CC squared 

      REAL,    PARAMETER :: AMUKEV= 0.931494012E6   ! [keV] 1 AMU 
      REAL,    PARAMETER :: AMUG  =1.660538783E-24  ! [g] AMU
      REAL,    PARAMETER :: MEKEV = 510.998902      ! [keV] electron mass 
      REAL,    PARAMETER :: MPKEV = 0.938271998E6   ! [keV] proton mass 
      REAL,    PARAMETER :: MDKEV =2.0141018*AMUKEV ! [keV] deuterium mas
      REAL,    PARAMETER :: MTKEV =3.0160492*AMUKEV ! [keV] triton masss
      REAL,    PARAMETER :: MPG   =1.0072765*AMUG   ! [g] deuterium mass
      REAL,    PARAMETER :: MDG   =2.0141018*AMUG   ! [g] deuterium mass
      REAL,    PARAMETER :: MTG   =3.0160492*AMUG   ! [g] triton mass

      REAL,    PARAMETER :: HBARC = 197.33E-10      ! [keV-cm ] hbar*c 197.5863737 eV nm ?? 
      REAL,    PARAMETER :: HBAR  = 6.5821E-19      ! [keV-s] hbar 
      REAL,    PARAMETER :: HBARC2= HBARC*HBARC     ! [keV^2-cm^2]
      REAL,    PARAMETER :: A0CM  = 5.29E-9         ! [cm] Bohr radius 
      REAL,    PARAMETER :: ALPHA = 1./137.035999679! 1/137.036 
! 
!     k_b = DEBYE*SQRT(zb*zb*nb*betab)           ! [1/cm]
      REAL,    PARAMETER :: DEBYE =4.2539E-5     ! Debye wave number
      REAL,    PARAMETER :: DEBYE2=DEBYE*DEBYE
!
!     om_b=(1.32155E+3)*SQRT(zb*zb*nb*AMUKEV/mb) ! [1/s]
      REAL,    PARAMETER :: OMEGI =1.32155E+3    ! plasma frequency
      REAL,    PARAMETER :: OMEGI2=OMEGI*OMEGI
!
!     ge=(6.1260E-15)*SQRT(ne)*betae**1.5        ! dimensionless
      REAL,    PARAMETER :: GECOEFF=6.1260E-15   ! plasma coupling

!     LOGICAL, PARAMETER :: SMALL_E_NLO=.TRUE.   ! T= LO + NLO
    LOGICAL, PARAMETER :: SMALL_E_NLO=.FALSE.    ! F= LO only 

!
! conversion factors
!
!      REAL,    PARAMETER :: KTOEV =8.61772E-5    ! K -> eV : [eV/K]
!      REAL,    PARAMETER :: CMTOA0=1.8867925E8   ! cm -> a0: [a0/cm]
!      REAL,    PARAMETER :: A0TOCM=1/CMTOA0      ! a0 -> cm: [cm/a0]
!      REAL,    PARAMETER :: CONVFACT=1.E-7/A0CM  ! keV/a0 -> keV/micron
!      REAL,    PARAMETER :: CONVFACT=CMTOA0*1.E-7! keV/a0 -> keV/micron
END MODULE physvars
