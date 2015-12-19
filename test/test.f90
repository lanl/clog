! Test suite program. 
!
! 1. To run test suite and compare with test.dat: write_dat=.FALSE.
! 2. To rewrite test.dat file set               : write_dat=.TRUE.
!
!
! This is the driver for the BPS coefficient subroutine.
!

    PROGRAM main
      USE bpsvars
      USE physvars
      USE allocatablevars 
      IMPLICIT NONE
      REAL    :: mp, ep, zp
      INTEGER :: nni 
      LOGICAL :: write_dat

      OPEN(UNIT=1,FILE='test.dat')
!     write_dat=.TRUE.  ! create test.dat
      write_dat=.FALSE. ! testing mode: read from test.dat

      CALL define_plasma(ep,mp,zp,nni)
      CALL data_write_test(write_dat,nni,ep,zp,mp,betab,zb,mb,nb)

      CLOSE(1)
    END PROGRAM main
!
!============================================
!
    SUBROUTINE define_plasma(ep, mp, zp, nni) !, te, ti, ne, nni)
    USE allocatablevars 
    USE physvars
      IMPLICIT NONE
      REAL,    INTENT(OUT) :: ep   ! projectile kinetic energy [keV]
      REAL,    INTENT(OUT) :: mp   ! projectile mass [keV]
      REAL,    INTENT(OUT) :: zp   ! projectile charge Z
      INTEGER, INTENT(OUT) :: nni  ! number of ion species
      REAL           :: te   ! electron temperature [keV]
      REAL           :: ti   ! ion temperature [keV]
      REAL           :: ne   ! [cm^-3]
!
!     REAL,    DIMENSION(:), ALLOCATABLE  :: betab  ! [1/keV]
!     REAL,    DIMENSION(:), ALLOCATABLE  :: mb     ! [keV]
!     REAL,    DIMENSION(:), ALLOCATABLE  :: nb     ! [cm^-3]
!     REAL,    DIMENSION(:), ALLOCATABLE  :: zb     ! [e]
!     REAL,    DIMENSION(:), ALLOCATABLE  :: gb, etab, mpb, mbpb

      nni=2  ! number of ion species

      ALLOCATE(betab(1:nni+1),zb(1:nni+1),mb(nni+1),nb(1:nni+1))   ! allocatablevars
      ALLOCATE(gb(1:nni+1),etab(1:nni+1),mpb(nni+1),mbpb(1:nni+1)) ! allocatablevars

!
! projectile parameters
!
      ep=3.54E3    ! Projectile energy  [keV]
      mp=4*MPKEV   ! Projectile mass    [keV]
      zp=2.        ! Projectile charge
      ne=1.E25     ! Electron number density    [cm^-3]
!
! plasma parameters
!
      te=10.       ! Electron temperature       [keV]
      ti=10.       ! Ion temperature            [keV]
      zb(1)=-1.    ! Species charges
      zb(2)=+1.    ! 
      zb(3)=+1.    !
      mb(1)=MEKEV  ! Species masses [keV]
      mb(2)=2*MPKEV!
      mb(3)=3*MPKEV!
!
! Construct density and temperature arrays
!
      nb(1)=1.                          ! ONLY FOR EQUIMOLAR DT
      nb(2:nni+1)=1./(zb(2:nni+1)*nni)  ! charge neutrality
      nb=nb*ne                          ! number density array [cm^-3]
      betab(1)=1./te                    ! inverse temp array   [keV^-1]
      betab(2:nni+1)=1./ti              !

    END SUBROUTINE define_plasma

!!=============
!
! DT plasma with alpha particle projectile
!
!      nni=2
!      ALLOCATE(mb(1:nni+1),nb(1:nni+1),zb(1:nni+1),betab(1:nni+1),gb(1:nni+1))
!      te=10.       ! Electron temperature       [keV]
!      ti=10.       ! Ion temperature            [keV]
!      ep=3.54E3    ! Projectile energy  [keV]
!      mp=4*MPKEV   ! Projectile mass    [keV]
!      zp=2.        ! Projectile charge
!      ne=1.E25     ! Electron number density    [cm^-3]
!
!      te=10.       ! Electron temperature       [keV]
!      ti=10.       ! Ion temperature            [keV]
!      zb(1)=-1.    ! Species charges
!      zb(2)=+1.    ! 
!      zb(3)=+1.    !
!      mb(1)=MEKEV  ! Species masses [keV]
!      mb(2)=2*MPKEV!
!      mb(3)=3*MPKEV!
!
! Construct density and temperature arrays
!
!      nb(1)=1.                          ! 
!      nb(2:nni+1)=1./(zb(2:nni+1)*nni)  ! charge neutrality
!      nb=nb*ne                          ! number density array [cm^-3]
!      betab(1)=1./te                    ! inverse temp array   [keV^-1]
!      betab(2:nni+1)=1./ti              !
!
!!=============

!
!============================================
!
    SUBROUTINE data_write_test(write_dat, nni, ep, zp, mp, betab, zb, mb, nb)
      IMPLICIT NONE
      LOGICAL,                     INTENT(IN)  :: write_dat
      INTEGER,                     INTENT(IN)  :: nni   !  number of ions
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: betab !  temp         [1/keV]
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: mb    !  mass array    [keV]
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: nb    !  density array [1/cc]
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: zb    !  charge array
                                                        !
                                                        ! Projectile  
      REAL,                        INTENT(IN)  :: ep    !  projectile energy [keV]
      REAL,                        INTENT(IN)  :: mp    !  projectile mass   [keV]
      REAL,                        INTENT(IN)  :: zp    !  projectile charge
!
! local variables
!
      IF (write_dat) THEN
         WRITE(6,*) 'generating test.dat file ...'
         WRITE(6,*) ''
         CALL write_plasma(nni,ep,zp,mp,betab,zb,mb,nb) 
         CALL write_param(nni,betab,nb)
         CALL write_dedx_bps(nni,ep,zp,mp,betab,zb,mb,nb)
         CALL write_acoeff_bps_mass(nni,ep,zp,mp,betab,zb,mb,nb)
         CALL write_bps_rate(nni,betab,zb,mb,nb)
     ELSE
         CALL test_plasma(nni,ep,zp,mp,betab,zb,mb,nb)         ! Did the plasma change?
         CALL test_param(nni,betab,nb)                         ! Did g, eta, etc. change?
         CALL test_dedx_bps(nni,ep,zp,mp,betab,zb,mb,nb)       ! Did stopping power change?
         CALL test_acoeff_bps_mass(nni,ep,zp,mp,betab,zb,mb,nb)! Did A-coeff change?
         CALL test_bps_rate(nni,betab,zb,mb,nb)                ! Did rates change?
!        CALL test_eifrac_int_1t(nni,ep,zp,mp,betab,zb,mb,nb,e0,ec)! Are energy fractions the same?
     ENDIF

    END SUBROUTINE data_write_test
!
! plasma definition
!
    SUBROUTINE write_plasma(nni, ep, zp, mp, betab, zb, mb, nb) 
      IMPLICIT NONE
                                                       ! Plasma:
      INTEGER,                     INTENT(IN)  :: nni  !  number of ions
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: betab!  temp array    [1/keV]
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: mb   !  mass array    [keV]
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: nb   !  density array [1/cc]
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: zb   !  charge array
                                                       !
                                                       ! Projectile  
      REAL,                        INTENT(IN)  :: ep   !  projectile energy [keV]
      REAL,                        INTENT(IN)  :: mp   !  projectile mass   [keV]
      REAL,                        INTENT(IN)  :: zp   !  projectile charge
      REAL :: te, ti

      te=1./betab(1)
      ti=1./betab(2)
      WRITE(6,'(16A)') 'plasma: DT alpha'
      WRITE(6,*) nni
      WRITE(6,'(2D25.14)')  te, ti
      WRITE(6,'(10D25.14)') nb
      WRITE(6,'(10D25.14)') mb
      WRITE(6,'(10D25.14)') zb
      WRITE(6,'(2D25.14)')  zp, mp
      WRITE(6,'(1D25.14)')  ep
      WRITE(6,*)''

      WRITE(1,'(16A)') 'plasma: DT alpha'
      WRITE(1,*) nni
      WRITE(1,'(2D25.14)')  te, ti
      WRITE(1,'(10D25.14)') nb
      WRITE(1,'(10D25.14)') mb
      WRITE(1,'(10D25.14)') zb
      WRITE(1,'(2D25.14)')  zp, mp
      WRITE(1,'(1D25.14)')  ep
    END SUBROUTINE write_plasma

    SUBROUTINE test_plasma(nni, ep, zp, mp, betab, zb, mb, nb) 
      IMPLICIT NONE
                                                       ! Plasma:
      INTEGER,                     INTENT(IN)  :: nni  !  number of ions
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: betab!  temp array    [1/keV]
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: mb   !  mass array    [keV]
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: nb   !  density array [1/cc]
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: zb   !  charge array
                                                       !
                                                       ! Projectile  
      REAL,                        INTENT(IN)  :: ep   !  projectile energy [keV]
      REAL,                        INTENT(IN)  :: mp   !  projectile mass   [keV]
      REAL,                        INTENT(IN)  :: zp   !  projectile charge

      REAL,    DIMENSION(:), ALLOCATABLE :: nb_dat, mb_dat, zb_dat
      CHARACTER*50     :: fname
      LOGICAL          :: pass=.TRUE. ! true or false
      REAL,  PARAMETER :: TOLERANCE=1.E-3
      real    :: err
      REAL    :: te, ti
      REAL    :: te_dat, ti_dat, mp_dat, ep_dat, zp_dat
      INTEGER :: nni_dat 

      WRITE(6,*) 'calling test_plasma ...'
      te=1./betab(1)
      ti=1./betab(2)

      READ(1,'(16A)') fname
      READ(1,*) nni_dat
      IF (nni .NE. nni_dat) THEN
        PRINT *,'ERROR: number of plasma species has changed'
        PRINT *,'ERROR: stopping now'
        STOP
      ENDIF

      ALLOCATE(mb_dat(1:nni_dat+1),nb_dat(1:nni_dat+1),zb_dat(1:nni_dat+1))

      READ(1,'(2D25.14)')  te_dat, ti_dat
      READ(1,'(10D25.14)') nb_dat
      READ(1,'(10D25.14)') mb_dat
      READ(1,'(10D25.14)') zb_dat
      READ(1,'(2D25.14)')  zp_dat, mp_dat
      READ(1,'(1D25.14)')  ep_dat

      err=100*ABS(te_dat-te)/(ABS(te_dat)+ABS(te))
      IF (err > TOLERANCE) THEN
        PRINT *, 'ERROR: Te'
        PRINT *, te
        PRINT *, te_dat
        PRINT *, err
        pass=.FALSE.
      ENDIF
      err=100*ABS(ti_dat-ti)/(ABS(ti_dat)+ABS(ti))
      IF (err > TOLERANCE) THEN
        PRINT *, 'ERROR: Ti'
        PRINT *, ti
        PRINT *, ti_dat
        PRINT *, err
        pass=.FALSE.
      ENDIF
      err=100*SUM(ABS(nb_dat-nb))/(SUM(ABS(nb_dat))+SUM(ABS(nb)))
      IF (err > TOLERANCE) THEN
        PRINT *, 'ERROR: nb'
        PRINT *, nb
        PRINT *, nb_dat
        PRINT *, err
        pass=.FALSE.
        err=100*SUM(ABS(mb_dat-mb))/(SUM(ABS(mb_dat))+SUM(ABS(mb)))
        PRINT *, err
        pass=.FALSE.
      ENDIF
      IF (err > TOLERANCE) THEN
        PRINT *, 'ERROR: mb'
        PRINT *, mb
        PRINT *, mb_dat
        PRINT *, err
        pass=.FALSE.
      ENDIF
      err=100*SUM(ABS(zb_dat-zb))/(SUM(ABS(zb_dat))+SUM(ABS(zb)))
      IF (err > TOLERANCE) THEN
        PRINT *, 'ERROR:zb'
        PRINT *, zb
        PRINT *, zb_dat
        PRINT *, err
        pass=.FALSE.
      ENDIF
      err=100*ABS(zp_dat-zp)/(ABS(zp_dat)+ABS(zp))
      IF (err > TOLERANCE) THEN
        PRINT *, 'ERROR: zp'
        PRINT *, zp
        PRINT *, zp_dat
        PRINT *, err
        pass=.FALSE.
      ENDIF
      err=100*ABS(mp_dat-mp)/(ABS(mp_dat)+ABS(mp))
      IF (err > TOLERANCE) THEN
        PRINT *, 'ERROR:mp'
        PRINT *, mp
        PRINT *, mp_dat
        PRINT *, err
        pass=.FALSE.
      ENDIF
      err=100*ABS(ep_dat-ep)/(ABS(ep_dat)+ABS(ep))
      IF (err > TOLERANCE) THEN
        PRINT *, 'ERROR:ep'
        PRINT *, ep
        PRINT *, ep_dat
        PRINT *, err
        pass=.FALSE.
      ENDIF
      IF (pass) WRITE(6,*) ' passed.'
      DEALLOCATE(mb_dat,nb_dat,zb_dat)
    END SUBROUTINE test_plasma
!
!============================================
!
    SUBROUTINE write_param(nni, betab, nb)
      IMPLICIT NONE
      INTEGER,                     INTENT(IN)  :: nni   !  number of ions
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: betab !  temp [1/keV]
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: nb    !  density array [1/cc]
      REAL,    DIMENSION(:), ALLOCATABLE :: gb
      REAL :: ge, gi, etae, ze

      ALLOCATE(gb(1:nni+1))      
      CALL param(nni, betab, nb, gb, ge, gi, etae, ze)
      WRITE(6,'(23A)') 'subroutine check: param'
      WRITE(6,'(10D25.14)') gb
      WRITE(6,'(2D25.14)')  ge, gi
      WRITE(6,'(2D25.14)')  etae, ze
      WRITE(1,'(23A)') 'subroutine check: param'
      WRITE(1,'(10D25.14)') gb
      WRITE(1,'(2D25.14)')  ge, gi
      WRITE(1,'(2D25.14)')  etae, ze
      DEALLOCATE(gb)
    END SUBROUTINE write_param

    SUBROUTINE test_param(nni, betab, nb)
      IMPLICIT NONE
      INTEGER,                     INTENT(IN)  :: nni   !  number of ions
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: betab !  temp [1/keV]
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: nb    !  density array [1/cc]

      REAL,  DIMENSION(:), ALLOCATABLE :: gb, gb_dat
      CHARACTER*50     :: fname
      LOGICAL          :: pass=.TRUE. ! true or false
      REAL,  PARAMETER :: TOLERANCE=1.E-3
      REAL :: err
      REAL :: ge, gi, etae, ze, ge_dat, gi_dat, etae_dat, ze_dat

      WRITE(6,*) 'calling: test_param ...'
      ALLOCATE(gb(1:nni+1),gb_dat(1:nni+1))

      CALL param(nni, betab, nb, gb, ge, gi, etae, ze)

      READ(1,'(23A)') fname
      READ(1,'(10D25.14)') gb_dat
      READ(1,'(2D25.14)')  ge_dat, gi_dat
      READ(1,'(2D25.14)')  etae_dat, ze_dat

      err=100*SUM(ABS(gb_dat-gb))/(SUM(ABS(gb_dat))+SUM(ABS(gb)))
      IF (err > TOLERANCE) THEN
        PRINT *, 'ERROR: gb'
        PRINT *, gb
        PRINT *, gb_dat
        PRINT *, err
        pass=.FALSE.
      ENDIF
      err=100*ABS(ge_dat-ge)
      IF (err > TOLERANCE) THEN
        PRINT *, 'ERROR: ge'
        PRINT *, ge
        PRINT *, ge_dat
        PRINT *, err
        pass=.FALSE.
      ENDIF
      err=100*ABS(gi_dat-gi)
      IF (err > TOLERANCE) THEN
        PRINT *, 'ERROR: gi'
        PRINT *, gi
        PRINT *, gi_dat
        PRINT *, err
        pass=.FALSE.
      ENDIF
      err=100*ABS(etae_dat-etae)
      IF (err > TOLERANCE) THEN
        PRINT *, 'ERROR: etae'
        PRINT *, etae
        PRINT *, etae_dat
        PRINT *, err
        pass=.FALSE.
      ENDIF
      err=100*ABS(ze_dat-ze)
      IF (err > TOLERANCE) THEN
        PRINT *, 'ERROR:ze'
        PRINT *, ze
        PRINT *, ze_dat
        PRINT *, err
        pass=.FALSE.
      ENDIF
      IF (pass) WRITE(6,*) ' passed.'
      DEALLOCATE(gb_dat)
    END SUBROUTINE test_param
!
!============================================
!

    SUBROUTINE write_dedx_bps(nni,ep,zp,mp,betab,zb,mb,nb)
      IMPLICIT NONE
                                                       ! Plasma:
      INTEGER,                     INTENT(IN)  :: nni  !  number of ions
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: betab!  temp array    [1/keV]
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: mb   !  mass array    [keV]
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: nb   !  density array [1/cc]
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: zb   !  charge array
                                                       !
                                                       ! Projectile  
      REAL,                        INTENT(IN)  :: ep   !  projectile energy [keV]
      REAL,                        INTENT(IN)  :: mp   !  projectile mass   [keV]
      REAL,                        INTENT(IN)  :: zp   !  projectile charge
      REAL :: dedx_tot,  dedx_i,  dedx_e, dedxc_tot, dedxc_i, dedxc_e
      REAL :: dedxq_tot, dedxq_i, dedxq_e

      CALL dedx_bps(nni,ep,zp,mp,betab,zb,mb,nb,dedx_tot,dedx_i,dedx_e,&
         dedxc_tot,dedxc_i,dedxc_e,dedxq_tot,dedxq_i,dedxq_e)
      WRITE(6,'(40A)') 'subroutine check: dedx_bps [MeV/micron]'
      WRITE(6,'(3D25.14)') dedx_tot,  dedx_i,  dedx_e
      WRITE(6,'(3D25.14)') dedxc_tot, dedxc_i, dedxc_e
      WRITE(6,'(3D25.14)') dedxq_tot, dedxq_i, dedxq_e
      WRITE(1,'(40A)') 'subroutine check: dedx_bps [MeV/micron]'
      WRITE(1,'(3D25.14)') dedx_tot,  dedx_i,  dedx_e
      WRITE(1,'(3D25.14)') dedxc_tot, dedxc_i, dedxc_e
      WRITE(1,'(3D25.14)') dedxq_tot, dedxq_i, dedxq_e
    END SUBROUTINE write_dedx_bps

    SUBROUTINE test_dedx_bps(nni,ep,zp,mp,betab,zb,mb,nb)
      IMPLICIT NONE
                                                       ! Plasma:
      INTEGER,                     INTENT(IN)  :: nni  !  number of ions
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: betab!  temp array    [1/keV]
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: mb   !  mass array    [keV]
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: nb   !  density array [1/cc]
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: zb   !  charge array
                                                       !
                                                       ! Projectile  
      REAL,                        INTENT(IN)  :: ep   !  projectile energy [keV]
      REAL,                        INTENT(IN)  :: mp   !  projectile mass   [keV]
      REAL,                        INTENT(IN)  :: zp   !  projectile charge

      CHARACTER*50     :: fname                                                               !
      LOGICAL          :: pass=.TRUE. ! true or false
      REAL,  PARAMETER :: TOLERANCE=1.E-3
      REAL :: err
      REAL :: dedx_tot,      dedx_i,     dedx_e
      REAL :: dedx_tot_dat,  dedx_i_dat, dedx_e_dat
      REAL :: dedxc_tot, dedxc_i, dedxc_tot_dat, dedxc_i_dat
      REAL :: dedxq_tot, dedxq_i, dedxq_tot_dat, dedxq_i_dat
      REAL :: dedxc_e, dedxq_e, dedxc_e_dat, dedxq_e_dat

      WRITE(6,*) 'calling: test_dedx_bps ...'

      CALL dedx_bps(nni,ep,zp,mp,betab,zb,mb,nb,dedx_tot,        &
              dedx_i,dedx_e,dedxc_tot,dedxc_i,dedxc_e,dedxq_tot, &
              dedxq_i,dedxq_e) !stopping power [keV/micron]

      READ(1,'(23A)')     fname
      READ(1,'(3D25.14)') dedx_tot_dat,  dedx_i_dat,  dedx_e_dat
      READ(1,'(3D25.14)') dedxc_tot_dat, dedxc_i_dat, dedxc_e_dat
      READ(1,'(3D25.14)') dedxq_tot_dat, dedxq_i_dat, dedxq_e_dat

      err=100*ABS(dedx_tot_dat-dedx_tot)/(ABS(dedx_tot_dat)+ABS(dedx_tot))
      IF (err > TOLERANCE) THEN
        PRINT *, 'ERROR: dedx_tot'
        PRINT *, dedx_tot
        PRINT *, dedx_tot_dat
        PRINT *, err
        pass=.FALSE.
      ENDIF
      err=100*ABS(dedx_i_dat-dedx_i)/(ABS(dedx_i_dat)+ABS(dedx_i))
      IF (err > TOLERANCE) THEN
        PRINT *, 'ERROR: dedx_i'
        PRINT *, dedx_i
        PRINT *, dedx_i_dat
        PRINT *, err
        pass=.FALSE.
      ENDIF
      err=100*ABS(dedx_e_dat-dedx_e)/(ABS(dedx_e_dat)+ABS(dedx_e))
      IF (err > TOLERANCE) THEN
        PRINT *, 'ERROR: dedx_e'
        PRINT *, dedx_e
        PRINT *, dedx_e_dat
        PRINT *, err
        pass=.FALSE.
      ENDIF
      err=100*ABS(dedxc_tot_dat-dedxc_tot)/(ABS(dedxc_tot_dat)+ABS(dedxc_tot))
      IF (err > TOLERANCE) THEN
        PRINT *, 'ERROR: dedxc_tot'
        PRINT *, dedxc_tot
        PRINT *, dedxc_tot_dat
        PRINT *, err
        pass=.FALSE.
      ENDIF
      err=100*ABS(dedxc_i_dat-dedxc_i)/(ABS(dedxc_i_dat)+ABS(dedxc_i))
      IF (err > TOLERANCE) THEN
        PRINT *, 'ERROR: dedxc_i'
        PRINT *, dedxc_i
        PRINT *, dedxc_i_dat
        PRINT *, err
        pass=.FALSE.
      ENDIF
      err=100*ABS(dedxc_e_dat-dedxc_e)/(ABS(dedxc_e_dat)+ABS(dedxc_e)) ! new
      IF (err > TOLERANCE) THEN
        PRINT *, 'ERROR: dedxc_e'
        PRINT *, dedxc_e
        PRINT *, dedxc_e_dat
        PRINT *, err
        pass=.FALSE.
      ENDIF
      err=100*ABS(dedxq_tot_dat-dedxq_tot)/(ABS(dedxq_tot_dat)+ABS(dedxq_tot))
      IF (err > TOLERANCE) THEN
        PRINT *, 'ERROR: dedxq_tot'
        PRINT *, dedxq_tot
        PRINT *, dedxq_tot_dat
        PRINT *, err
        pass=.FALSE.
      ENDIF
      err=100*ABS(dedxq_i_dat-dedxq_i)/(ABS(dedxq_i_dat)+ABS(dedxq_i))
      IF (err > TOLERANCE) THEN
        PRINT *, 'ERROR: dedxq_i'
        PRINT *, dedxq_i
        PRINT *, dedxq_i_dat
        PRINT *, err
        pass=.FALSE.
      ENDIF
      err=100*ABS(dedxq_e_dat-dedxq_e)/(ABS(dedxq_e_dat)+ABS(dedxq_e)) ! new
      IF (err > TOLERANCE) THEN
        PRINT *, 'ERROR: dedxq_e'
        PRINT *, dedxq_e
        PRINT *, dedxq_e_dat
        PRINT *, err
        pass=.FALSE.
      ENDIF
      IF (pass) WRITE(6,*) ' passed.'
    END SUBROUTINE test_dedx_bps
!
!============================================
!
    SUBROUTINE write_acoeff_bps_mass(nni,ep,zp,mp,betab,zb,mb,nb)
                                                         ! Plasma:
      INTEGER,                     INTENT(IN)  :: nni    !  number of ions
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: betab  !  temp array    [1/keV]
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: mb     !  mass array    [keV]
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: nb     !  density array [1/cc]
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: zb     !  charge array
                                                         !
                                                         ! Projectile  
      REAL,                        INTENT(IN)  :: ep     !  projectile energy [keV]
      REAL,                        INTENT(IN)  :: mp     !  projectile mass   [keV]
      REAL,                        INTENT(IN)  :: zp     !  projectile charge
      REAL :: a_tot, a_i, a_e, ac_tot, ac_i, ac_e, aq_tot, aq_i, aq_e
      REAL :: ac_s_e, ac_r_e, ac_s_i, ac_r_i

      CALL bps_acoeff_ei_mass(nni, ep, zp, mp, betab, zb, mb, nb,    &
         a_tot, a_i, a_e, ac_tot, ac_i, ac_e, aq_tot, aq_i, aq_e,&
         ac_s_e, ac_r_e, ac_s_i, ac_r_i) 
      WRITE(6,'(54A)') 'subroutine check: bps_acoeff_ei_mass [A]=[MeV/micron]'
      WRITE(6,'(3D25.14)') a_tot, a_i, a_e
      WRITE(6,'(3D25.14)') ac_tot, ac_i, ac_e
      WRITE(6,'(3D25.14)') aq_tot, aq_i, ac_e
      WRITE(1,'(54A)') 'subroutine check: bps_acoeff_ei_mass [A]=[MeV/micron]'
      WRITE(1,'(3D25.14)') a_tot, a_i, a_e
      WRITE(1,'(3D25.14)') ac_tot, ac_i, ac_e
      WRITE(1,'(3D25.14)') aq_tot, aq_i, aq_e
    END SUBROUTINE write_acoeff_bps_mass

    SUBROUTINE test_acoeff_bps_mass(nni,ep,zp,mp,betab,zb,mb,nb)
                                                         ! Plasma:
      INTEGER,                     INTENT(IN)  :: nni    !  number of ions
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: betab  !  temp array    [1/keV]
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: mb     !  mass array    [keV]
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: nb     !  density array [1/cc]
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: zb     !  charge array
                                                         !
                                                         ! Projectile  
      REAL,                        INTENT(IN)  :: ep     !  projectile energy [keV]
      REAL,                        INTENT(IN)  :: mp     !  projectile mass   [keV]
      REAL,                        INTENT(IN)  :: zp     !  projectile charge

                                                         !
      CHARACTER*50     :: fname                          !A-coeffs [MeV/micron]
      LOGICAL          :: pass=.TRUE. ! true or false
      REAL,  PARAMETER :: TOLERANCE=1.E-3
      REAL :: err
      REAL :: a_tot, a_tot_dat          ! electron + ion
      REAL :: a_i, a_i_dat              ! ion contribution
      REAL :: a_e, a_e_dat              ! electorn contribution
      REAL :: ac_tot, ac_tot_dat        ! classical
      REAL :: ac_i, ac_i_dat, ac_e_dat  ! classical
      REAL :: aq_tot, aq_tot_dat        ! quantum
      REAL :: aq_i, aq_i_dat, aq_e_dat  ! quantum
      REAL :: ac_s_e, ac_s_i            ! singular [classical=sing+reg]
      REAL :: ac_r_e, ac_r_i            ! regular  [classical=sing+reg]


      WRITE(6,*) 'calling: test_coeff_bps_mass ...'

      CALL bps_acoeff_ei_mass(nni, ep, zp, mp, betab, zb, mb, nb,    &
         a_tot, a_i, a_e, ac_tot, ac_i, ac_e, aq_tot, aq_i, aq_e,&
         ac_s_e, ac_r_e, ac_s_i, ac_r_i)

      READ(1,'(23A)') fname
      READ(1,'(3D25.14)') a_tot_dat, a_i_dat, a_e_dat
      READ(1,'(3D25.14)') ac_tot_dat, ac_i_dat, ac_e_dat
      READ(1,'(3D25.14)') aq_tot_dat, aq_i_dat, aq_e_dat

      err=100*ABS(a_tot_dat-a_tot)/(ABS(a_tot_dat)+ABS(a_tot))
      IF (err > TOLERANCE) THEN
        PRINT *, 'ERROR: a_tot'
        PRINT *, a_tot
        PRINT *, a_tot_dat
        PRINT *, err
        pass=.FALSE.
      ENDIF
      err=100*ABS(a_i_dat-a_i)/(ABS(a_i_dat)+ABS(a_i))
      IF (err > TOLERANCE) THEN
        PRINT *, 'ERROR: a_i'
        PRINT *, a_i
        PRINT *, a_i_dat
        PRINT *, err
        pass=.FALSE.
      ENDIF
      err=100*ABS(a_e_dat-a_e)/(ABS(a_e_dat)+ABS(a_e))
      IF (err > TOLERANCE) THEN
        PRINT *, 'ERROR: a_e'
        PRINT *, a_e
        PRINT *, a_e_dat
        PRINT *, err
        pass=.FALSE.
      ENDIF
      err=100*ABS(ac_tot_dat-ac_tot)/(ABS(ac_tot_dat)+ABS(ac_tot))
      IF (err > TOLERANCE) THEN
        PRINT *, 'ERROR: ac_tot'
        PRINT *, ac_tot
        PRINT *, ac_tot_dat
        PRINT *, err
        pass=.FALSE.
      ENDIF
      err=100*ABS(ac_i_dat-ac_i)/(ABS(ac_i_dat)+ABS(ac_i))
      IF (err > TOLERANCE) THEN
        PRINT *, 'ERROR: ac_i'
        PRINT *, ac_i
        PRINT *, ac_i_dat
        PRINT *, err
        pass=.FALSE.
      ENDIF
      err=100*ABS(ac_e_dat-ac_e)/(ABS(ac_e_dat)+ABS(ac_e))
      IF (err > TOLERANCE) THEN
        PRINT *, 'ERROR: ac_e'
        PRINT *, ac_e
        PRINT *, ac_e_dat
        PRINT *, err
        pass=.FALSE.
      ENDIF
      err=100*ABS(aq_tot_dat-aq_tot)/(ABS(aq_tot_dat)+ABS(aq_tot))
      IF (err > TOLERANCE) THEN
        PRINT *, 'ERROR: aq_tot'
        PRINT *, aq_tot
        PRINT *, aq_tot_dat
        PRINT *, err
        pass=.FALSE.
      ENDIF
      err=100*ABS(aq_i_dat-aq_i)/(ABS(aq_i_dat)+ABS(aq_i))
      IF (err > TOLERANCE) THEN
        PRINT *, 'ERROR: aq_i'
        PRINT *, aq_i
        PRINT *, aq_i_dat
        PRINT *, err
        pass=.FALSE.
      ENDIF
      err=100*ABS(aq_e_dat-aq_e)/(ABS(aq_e_dat)+ABS(aq_e))
      IF (err > TOLERANCE) THEN
        PRINT *, 'ERROR: aq_e'
        PRINT *, aq_e
        PRINT *, aq_e_dat
        PRINT *, err
        pass=.FALSE.
      ENDIF
      IF (pass) WRITE(6,*) ' passed.'
    END SUBROUTINE test_acoeff_bps_mass
!
!============================================  ***===***
!
    SUBROUTINE write_bps_rate(nni,betab,zb,mb,nb)
      IMPLICIT NONE
      INTEGER,                     INTENT(IN)  :: nni    !  Number of ions
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: betab  !  Temperature array [1/keV]
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: zb     !  Charge array 
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: mb     !  Mass array [keV]
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: nb     !  Number density array [cm^-3]

      REAL     :: ln_bps_mass, cei_mass, delta_mass, cei_sing_mass, cei_reg_mass, cei_qm_mass
      REAL     :: ln_bps_born, cei_born

      REAL     :: cei_tot, cei_e, ceic_tot, ceic_i, ceic_e
      REAL     :: ceiq_tot, ceiq_e, ceic_s_e
      REAL     :: ceic_r_e
      REAL,    DIMENSION(1:nni+1)  :: ceib
      INTEGER  :: ib


      CALL bps_rate_cei_born(nni, betab, zb, mb, nb, ln_bps_born, cei_born)

      WRITE(6,'(56A)') 'subroutine check: bps_rate_cab_mass, etc.  [C]=[1/cm^3 s]'
      WRITE(1,'(56A)') 'subroutine check: bps_rate_cab_mass, etc.  [C]=[1/cm^3 s]'

      WRITE(6,'(2D25.14)') ln_bps_born, cei_born
      WRITE(1,'(2D25.14)') ln_bps_born, cei_born

      CALL bps_rate_cei_mass(nni, betab, zb, mb, nb, ln_bps_mass,            & 
      delta_mass, cei_tot, cei_mass, cei_e, ceic_tot, ceic_i, ceic_e, ceiq_tot, &
      cei_qm_mass, ceiq_e, cei_sing_mass, ceic_s_e, cei_reg_mass , ceic_r_e, ceib)

      WRITE(6,'(2D25.14)') ln_bps_mass, delta_mass
      WRITE(6,'(3D25.14)') cei_tot, cei_mass, cei_e
      WRITE(6,'(3D25.14)') ceiq_tot, cei_qm_mass, ceiq_e
      WRITE(1,'(2D25.14)') ln_bps_mass, delta_mass
      WRITE(1,'(3D25.14)') cei_tot, cei_mass, cei_e
      WRITE(1,'(3D25.14)') ceiq_tot, cei_qm_mass, ceiq_e
      DO ib=1,nni+1
         WRITE(6,'(3D25.14)') ceib(ib)
         WRITE(1,'(3D25.14)') ceib(ib)
      ENDDO

!cei_sing_mass, ceic_s_e, cei_reg_mass , ceic_r_e, ceib



!      CALL bps_rate_cab_matrix(nni, betab, zb, mb, nb,    &
!        cab, cab_sing, cab_reg, cab_qm, c_tot, c_i, c_e, cc_tot, &
!        cc_i, cc_e, cq_tot, cq_i, cq_e, cc_s_i, cc_s_e, cc_r_i, cc_r_e)



    END SUBROUTINE write_bps_rate

    SUBROUTINE test_bps_rate(nni,betab,zb,mb,nb)
      IMPLICIT NONE
      INTEGER,                     INTENT(IN)  :: nni    !  Number of ions
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: betab  !  Temperature array [1/keV]
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: zb     !  Charge array 
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: mb     !  Mass array [keV]
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: nb     !  Number density array [cm^-3]

      CHARACTER*50     :: fname                          
      LOGICAL          :: pass=.TRUE. ! true or false
      REAL,  PARAMETER :: TOLERANCE=1.E-3
      REAL :: err
      REAL :: ln_bps_born, ln_bps_born_dat  ! BPS Coulomb log in extreme quantum limit
      REAL :: ln_bps_mass, ln_bps_mass_dat  ! BPS Coulomb log in arbitrary regime
      REAL :: cei_born, cei_born_dat        ! BPS rate in extreme quantum limit
      REAL :: cei_mass, cei_mass_dat        ! BPS rate in arbitrary regime
      REAL :: delta_mass, delta_mass_dat    !
      REAL :: cei_tot, cei_tot_dat          !
      REAL :: cei_e, cei_e_dat              !
      REAL :: ceiq_tot, ceiq_tot_dat        !
      REAL :: cei_qm_mass, cei_qm_mass_dat  !
      REAL :: ceiq_e, ceiq_e_dat            !
      REAL, DIMENSION(1:nni+1)  :: ceib, ceib_dat
      INTEGER :: ib

      REAL :: cei_reg_mass    ! BPS rate in arbitrary regime
      REAL :: cei_sing_mass   ! BPS rate in arbitrary regime
      REAL :: ceic_tot, ceic_i, ceic_e
      REAL :: ceic_s_e
      REAL :: ceic_r_e



      WRITE(6,*) 'calling: test_bps_rate ...'
      READ(1,'(64A)') fname
!
! fetching data from file and code for comparison ...
!
      CALL bps_rate_cei_born(nni, betab, zb, mb, nb, ln_bps_born, cei_born)
      READ(1,'(2D25.14)') ln_bps_born_dat, cei_born_dat

      err=100*ABS(ln_bps_born_dat-ln_bps_born)/(ABS(ln_bps_born_dat)+ABS(ln_bps_born))
      IF (err > TOLERANCE) THEN
      PRINT *, 'ERROR: ln_bps_mass'
      PRINT *, ln_bps_born
      PRINT *, ln_bps_born_dat
      PRINT *, err
      pass=.FALSE.
      ENDIF
      err=100*ABS(cei_born_dat-cei_born)/(ABS(cei_born_dat)+ABS(cei_born))
      IF (err > TOLERANCE) THEN
      PRINT *, 'ERROR: ln_bps_mass'
      PRINT *, cei_born
      PRINT *, cei_born_dat
      PRINT *, err
      pass=.FALSE.
      ENDIF
!
!
      CALL bps_rate_cei_mass(nni, betab, zb, mb, nb, ln_bps_mass,            & 
      delta_mass, cei_tot, cei_mass, cei_e, ceic_tot, ceic_i, ceic_e, ceiq_tot, &
      cei_qm_mass, ceiq_e, cei_sing_mass, ceic_s_e, cei_reg_mass , ceic_r_e, ceib)
      READ(1,'(2D25.14)') ln_bps_mass_dat, delta_mass_dat
      READ(1,'(3D25.14)') cei_tot_dat, cei_mass_dat, cei_e_dat
      READ(1,'(3D25.14)') ceiq_tot_dat, cei_qm_mass_dat, ceiq_e_dat
      DO ib=1,nni+1
         READ(1,'(3D25.14)') ceib_dat(ib)
      ENDDO

      err=100*ABS(ln_bps_mass_dat-ln_bps_mass)/(ABS(ln_bps_mass_dat)+ABS(ln_bps_mass))
      IF (err > TOLERANCE) THEN
      PRINT *, 'ERROR: ln_bps_mass'
      PRINT *, ln_bps_mass
      PRINT *, ln_bps_mass_dat
      PRINT *, err
      pass=.FALSE.
      ENDIF
      err=100*ABS(delta_mass_dat-delta_mass)/(ABS(delta_mass_dat)+ABS(delta_mass))
      IF (err > TOLERANCE) THEN
      PRINT *, 'ERROR: delta_mass'
      PRINT *, delta_mass
      PRINT *, delta_mass_dat
      PRINT *, err
      pass=.FALSE.
      ENDIF
      err=100*ABS(cei_tot_dat-cei_tot)/(ABS(cei_tot_dat)+ABS(cei_tot))
      IF (err > TOLERANCE) THEN
      PRINT *, 'ERROR: cei_tot'
      PRINT *, cei_tot
      PRINT *, cei_tot_dat
      PRINT *, err
      pass=.FALSE.
      ENDIF
      err=100*ABS(cei_mass_dat-cei_mass)/(ABS(cei_mass_dat)+ABS(cei_mass))
      IF (err > TOLERANCE) THEN
      PRINT *, 'ERROR: cei_mass'
      PRINT *, cei_mass
      PRINT *, cei_mass_dat
      PRINT *, err
      pass=.FALSE.
      ENDIF
      err=100*ABS(cei_e_dat-cei_e)/(ABS(cei_e_dat)+ABS(cei_e))
      IF (err > TOLERANCE) THEN
      PRINT *, 'ERROR: cei_e'
      PRINT *, cei_e
      PRINT *, cei_e_dat
      PRINT *, err
      pass=.FALSE.
      ENDIF
      err=100*ABS(ceiq_tot_dat-ceiq_tot)/(ABS(ceiq_tot_dat)+ABS(ceiq_tot))
      IF (err > TOLERANCE) THEN
      PRINT *, 'ERROR: ceiq_tot'
      PRINT *, ceiq_tot
      PRINT *, ceiq_tot_dat
      PRINT *, err
      pass=.FALSE.
      ENDIF
      err=100*ABS(cei_qm_mass_dat-cei_qm_mass)/(ABS(cei_qm_mass_dat)+ABS(cei_qm_mass))
      IF (err > TOLERANCE) THEN
      PRINT *, 'ERROR: cei_qm_mass'
      PRINT *, cei_qm_mass
      PRINT *, cei_qm_mass_dat
      PRINT *, err
      pass=.FALSE.
      ENDIF
      err=100*ABS(ceiq_e_dat-ceiq_e)/(ABS(ceiq_e_dat)+ABS(ceiq_e))
      IF (err > TOLERANCE) THEN
      PRINT *, 'ERROR: ceiq_e'
      PRINT *, ceiq_e
      PRINT *, ceiq_e_dat
      PRINT *, err
      pass=.FALSE.
      ENDIF
      DO ib=1, nni+1
        err=100*ABS(ceib_dat(ib)-ceib(ib))/(ABS(ceib_dat(ib))+ABS(ceib(ib)))
        IF (err > TOLERANCE) THEN
        PRINT *, 'ERROR: ceiq_e'
        PRINT *, ceib(ib)
        PRINT *, ceib_dat(ib)
        PRINT *, err
        pass=.FALSE.
        ENDIF
      ENDDO

!!      CALL bps_rate_cab_matrix(nni, betab, zb, mb, nb,    &
!!        cab, cab_sing, cab_reg, cab_qm, c_tot, c_i, c_e, cc_tot, &
!!        cc_i, cc_e, cq_tot, cq_i, cq_e, cc_s_i, cc_s_e, cc_r_i, cc_r_e)


      IF (pass) WRITE(6,*) ' passed.'


    END SUBROUTINE test_bps_rate
!
!============================================
!

