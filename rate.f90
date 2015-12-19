!
! bps_rate_cab_mass: returns the rate C_ab, the corresponding Coulomb
! logarithm, and the function Delta (C_reg = -1/2 + Delta) for a
! general electron mass. While Delta -> 0 as me -> 0, this limit is 
! not assumed in this subroutine.
!
!
! ROUTINES: bps_rate_cab_mass(nni, betab, zb, mb, nb, delta ...)
!           bps_rate_cei_mass(nni, betab, zb, mb, nb, cei, delta ...)
! 
! Assume a plasma composed of several species b, each separately in 
! thermal equilibrium with themselves but not necessarily with each 
! other[1]. This routine returns several useful components of the 
! corresponding rate C-coefficients.  introduced in Note [2] below (BPS).
! 
! The rate of change of energy density is given by 
!
!  dE_{ab}
!  -------- = -C_{ab} (T_a - T_b) 
!    dt
!
! with C_{ab} = ( C_{ab}_reg + C_{ab}_sing ) + C_{ab}_qm
!
!
! UNITS: C_{ab} has units of [1/cm^2*s] 
!
! INPUT: nni, betab, zb, mb, nb
!
! plasma input quantities:
! nni  : Number of total plasma species = number ion species + 1
! zb   : Charges of the plasma species. By convention zp(1) is the 
!      : electron plasma component. [dimensionless, Array]
! betab: Inverse temperatures of the plasma components. For an
!        electron-ion plasma, set betab(1)=1/T_e and all other
!        values of the array to 1/T_I [keV^-1].
! mb   : Masses of the plasma species [keV]. 
! nb   : Number densities of the plasma species [cm^-3]. 
!
!
! OUTPUT: ln_bps, delta, cab(1:nni+1,1:nni+1), cab_sing(:,:), cab_reg(:,:), cab_sqm(:,:)
!                   |
!                   -----> 0  as  me --> 0
!
! Each plasma component b makes a linear contribution A_b to the total 
! A-coefficient, i.e. A = sum_b A_b [5]. Each A_b in turn can be be 
! decomposed into a classical-quantum or electron-ion contributions.
!
! Coulomb logarithm     : ln_bps
! General rate coeff    : cei = \sum_i c_{ei} [cm^-3 s^-1]
! classical total       : delta --> 0  as  me --> 0
! Singular contribution : cei_sing
! Regular  contribution : cei_reg
! Quantum  contribution : cei_qm
!
! NOTES:
! [1] The temperatures T_b may therefore all differ. By convention
!     I take b=1 for the electron component of the plasma. A very
!     useful and interesting parameter regime is the one in which 
!     the ions have a common temperature T_I and the electron have
!     a temperature T_e, usually with T_e =/= T_I.  See also USAGE
!     and note [3] below.
!
! [2] L. Brown, D. Preston, and R. Singleton~Jr., 
!     "Charged Particle Motion in a Highly Ionized Plasma",
!     Physics Reports, 410 (2005) 237
!     [arXiv:physics/0501084]
!
!
      SUBROUTINE bps_rate_cab_mass(nni, ia, ib, betab, zb, mb, nb, &
        c_ab, c_ab_sing, c_ab_reg, c_ab_qm)
      USE mathvars
      USE physvars
        IMPLICIT NONE
        INTEGER,                             INTENT(IN)  :: nni     !  Number of ions
        INTEGER,                             INTENT(IN)  :: ia      !  Species number
        INTEGER,                             INTENT(IN)  :: ib      !  Species number
        REAL,    DIMENSION(1:nni+1),         INTENT(IN)  :: betab   !  Temperature array    [1/keV]
        REAL,    DIMENSION(1:nni+1),         INTENT(IN)  :: zb      !  Charge array 
        REAL,    DIMENSION(1:nni+1),         INTENT(IN)  :: mb      !  Mass array [keV]
        REAL,    DIMENSION(1:nni+1),         INTENT(IN)  :: nb      !  Number density array [cm^-3]
        REAL,                                INTENT(OUT) :: c_ab
        REAL,                                INTENT(OUT) :: c_ab_sing
        REAL,                                INTENT(OUT) :: c_ab_reg
        REAL,                                INTENT(OUT) :: c_ab_qm
        REAL,    DIMENSION(1:nni+1)  :: kb2, ob2
        REAL    ::  c_s, c_r, c_q
        REAL    :: kia2, kib2, bmia, bmib, nab, omi2, nab_reg
        INTEGER :: ibmax

        ibmax=nni+1
!
! initialize components of rate-coefficients
!
!
! construct plasma quantities: kb^2
!
        ob2=8*PI*A0CM*BEKEV*zb*zb*nb*CC2/mb
        kb2=8*PI*A0CM*BEKEV*zb*zb*nb*betab
        omi2=SUM(ob2(2:ibmax))

        kia2  =kb2(ia)
        bmia  =betab(ia)*mb(ia)
        kib2=kb2(ib)
        bmib=betab(ib)*mb(ib)
        nab=kia2*kib2*CC*SQRT(bmia*bmib)/(bmia + bmib)**1.5   !
        nab=nab/TWOPI**1.5                                    ! normalization
!
! C_{ab}-classical-singular 
!
        CALL cab_sing_mass(nni,ia,ib,betab,zb,mb,nb,c_s)
        c_ab_sing=nab*c_s 
!
! C_{ab}-classical-regular 
!
        CALL cab_reg_mass(nni,ia,ib,betab,zb,mb,nb,c_r)
        nab_reg=kb2(ia)*omi2/TWOPI
        nab_reg=nab_reg*SQRT(betab(ia)*mb(ia)/TWOPI)/CC
        c_ab_reg=nab_reg*c_r  
!
! C_{ab}-quantum
!
        CALL cab_qm_mass(nni,ia,ib,betab,zb,mb,c_q)
        c_ab_qm=nab*c_q
!
! C_{ab}-total
!
        c_ab=c_ab_sing + c_ab_reg + c_ab_qm
      END SUBROUTINE bps_rate_cab_mass
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Assembles the matrix C_{ab} of the rate coefficients
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!

      SUBROUTINE bps_rate_cab_matrix(nni, betab, zb, mb, nb,    &
        c_ab, c_ab_sing, c_ab_reg, c_ab_qm, c_tot, c_i, c_e, cc_tot, &
        cc_i, cc_e, cq_tot, cq_i, cq_e, cc_s_i, cc_s_e, cc_r_i, cc_r_e)
      USE physvars
      USE mathvars    
        IMPLICIT NONE                                             ! Plasma:
        INTEGER,                            INTENT(IN)  :: nni    !  number of ions
        REAL,    DIMENSION(1:nni+1),        INTENT(IN)  :: betab  !  temp array [1/keV]
        REAL,    DIMENSION(1:nni+1),        INTENT(IN)  :: zb     !  charge array
        REAL,    DIMENSION(1:nni+1),        INTENT(IN)  :: mb     !  mass array [keV]
        REAL,    DIMENSION(1:nni+1),        INTENT(IN)  :: nb     !  density [1/cc]
                                                                  !
                                                                  ! A-coeffs [MeV/micron]
        REAL,    DIMENSION(1:nni+1,1:nni+1),INTENT(OUT) :: c_ab
        REAL,    DIMENSION(1:nni+1,1:nni+1),INTENT(OUT) :: c_ab_sing
        REAL,    DIMENSION(1:nni+1,1:nni+1),INTENT(OUT) :: c_ab_reg
        REAL,    DIMENSION(1:nni+1,1:nni+1),INTENT(OUT) :: c_ab_qm
        REAL,    DIMENSION(1:nni+1),        INTENT(OUT) :: c_tot
        REAL,    DIMENSION(1:nni+1),        INTENT(OUT) :: c_i
        REAL,    DIMENSION(1:nni+1),        INTENT(OUT) :: c_e
        REAL,    DIMENSION(1:nni+1),        INTENT(OUT) :: cc_tot
        REAL,    DIMENSION(1:nni+1),        INTENT(OUT) :: cc_i
        REAL,    DIMENSION(1:nni+1),        INTENT(OUT) :: cc_e
        REAL,    DIMENSION(1:nni+1),        INTENT(OUT) :: cq_tot
        REAL,    DIMENSION(1:nni+1),        INTENT(OUT) :: cq_i
        REAL,    DIMENSION(1:nni+1),        INTENT(OUT) :: cq_e
        REAL,    DIMENSION(1:nni+1),        INTENT(OUT) :: cc_s_i
        REAL,    DIMENSION(1:nni+1),        INTENT(OUT) :: cc_s_e
        REAL,    DIMENSION(1:nni+1),        INTENT(OUT) :: cc_r_i
        REAL,    DIMENSION(1:nni+1),        INTENT(OUT) :: cc_r_e

        REAL    :: cab, cab_sing, cab_reg, cab_qm
        REAL    :: mp, zp
        INTEGER :: ia, ib

        c_i   = 0
        cc_s_i= 0
        cc_r_i= 0
        cc_i  = 0
        cq_i  = 0
        DO ia=1,nni+1
          mp=mb(ia)
          zp=zb(ia)
          DO ib=1,nni+1
            CALL bps_rate_cab_mass(nni, ia, ib, betab, zb, mb, nb, &
              cab, cab_sing, cab_reg, cab_qm)
            c_ab(ia,ib)     =cab
            c_ab_sing(ia,ib)=cab_sing
            c_ab_reg(ia,ib) =cab_reg
            c_ab_qm(ia,ib)  =cab_qm 
            IF (ib == 1) THEN
               c_e(ia)   = cab 
               cc_s_e(ia)= cab_sing
               cc_r_e(ia)= cab_reg
               cc_e(ia)  = cab_sing + cab_reg
               cq_e(ia)  = cab_qm
            ELSE
               c_i(ia)   = c_i(ia)    + cab 
               cc_s_i(ia)= cc_s_i(ia) + cab_sing
               cc_r_i(ia)= cc_r_i(ia) + cab_reg
               cc_i(ia)  = cc_i(ia)   + cab_sing + cab_reg
               cq_i(ia)  = cq_i(ia)   + cab_qm
            ENDIF
          ENDDO
          c_tot(ia) = c_e(ia)  + c_i(ia)
          cc_tot(ia)= cc_e(ia) + cc_i(ia)
          cq_tot(ia)= cq_e(ia) + cq_i(ia)
        ENDDO
      END SUBROUTINE bps_rate_cab_matrix
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Returns C_{e I} = \sum_i C_{e i} for backward compatibility
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
      SUBROUTINE bps_rate_cei_mass(nni, betab, zb, mb, nb, ln_bps,      & 
      delta, cei_tot, cei_i, cei_e, ceic_tot, ceic_i, ceic_e, ceiq_tot, &
      ceiq_i, ceiq_e, ceic_s_i, ceic_s_e, ceic_r_i , ceic_r_e, ceib)
      USE physvars
      USE mathvars    
        IMPLICIT NONE  
        INTEGER,                     INTENT(IN)  :: nni     !  Number of ions
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: betab   !  Temperature array    [1/keV]
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: zb      !  Charge array 
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: mb      !  Mass array [keV]
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: nb      !  Number density array [cm^-3]
        REAL,                        INTENT(OUT) :: ln_bps  !  Coulomb logarithm
        REAL,                        INTENT(OUT) :: delta   !  C_reg =-1/2+delta
        REAL,                        INTENT(OUT) :: cei_tot ! [cm^-3 s^-1]
        REAL,                        INTENT(OUT) :: cei_i   !
        REAL,                        INTENT(OUT) :: cei_e   !
        REAL,                        INTENT(OUT) :: ceic_tot!
        REAL,                        INTENT(OUT) :: ceic_i  !
        REAL,                        INTENT(OUT) :: ceic_e  !
        REAL,                        INTENT(OUT) :: ceiq_tot!
        REAL,                        INTENT(OUT) :: ceiq_i  !
        REAL,                        INTENT(OUT) :: ceiq_e  !
        REAL,                        INTENT(OUT) :: ceic_s_i!
        REAL,                        INTENT(OUT) :: ceic_s_e!
        REAL,                        INTENT(OUT) :: ceic_r_i!
        REAL,                        INTENT(OUT) :: ceic_r_e!
        REAL,    DIMENSION(1:nni+1), INTENT(OUT) :: ceib    !
        REAL,    DIMENSION(1:nni+1)  :: ob2
        REAL     :: omi2, cn, kb2e
        REAL     :: cab, c_ab_sing, c_ab_reg, c_ab_qm
        INTEGER  :: ia, ib, nnb
!
! initialize components of A-coefficients
!
        delta=0
        cei_tot =0  ! electron + ion
        cei_i   =0  ! ion contribution
        cei_e   =0  ! electron contribution
        ceic_tot=0  ! classical total
        ceic_e  =0  ! classical electron
        ceic_i  =0  ! classical ion
        ceiq_tot=0  ! quantum total
        ceiq_e  =0  ! quantum electron
        ceiq_i  =0  ! quantum ion
        ceic_s_i=0 
        ceic_s_e=0 
        ceic_r_i=0
        ceic_r_e=0
        ceib=0

        kb2e=8*PI*A0CM*BEKEV*zb(1)*zb(1)*nb(1)*betab(1)
        ob2=8*PI*A0CM*BEKEV*zb*zb*nb*CC2/mb
        omi2=SUM(ob2(2:nni+1))
        cn=kb2e*omi2/TWOPI
        cn=cn*SQRT(betab(1)*mb(1)/TWOPI)/CC

        ia=1
        NNB = nni+1    ! number of ions + electrons
        DO ib=1,nni+1  ! loop over electrons to calculate cei_e
           CALL bps_rate_cab_mass(nni, ia, ib, betab, zb, mb, nb, &
           cab, c_ab_sing, c_ab_reg, c_ab_qm)
           IF (ib .NE. 1) delta=delta+c_ab_reg/cn
           ceib(ib)=c_ab_sing + c_ab_reg + c_ab_qm
           CALL x_collect(ib, NNB, c_ab_sing, c_ab_reg, c_ab_qm,      &
           cei_tot, cei_i, cei_e, ceic_tot, ceic_i, ceic_e, ceiq_tot, &
           ceiq_i, ceiq_e, ceic_s_i, ceic_s_e, ceic_r_i, ceic_r_e)
        ENDDO
        delta=0.5+delta
        ln_bps=cei_i/cn
      END SUBROUTINE bps_rate_cei_mass
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! singular contribution c_s for non-zero electron mass
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Returns the rate coefficient C_{ia ib, sing}
!
      SUBROUTINE cab_sing_mass(nni, ia, ib, betab, zb, mb, nb, cab_sing)
      USE mathvars
      USE physvars
        IMPLICIT NONE
        INTEGER,                     INTENT(IN)  :: nni  !Number of ion species
        INTEGER,                     INTENT(IN)  :: ia
        INTEGER,                     INTENT(IN)  :: ib
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: zb   !Charge array
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: betab!Temperature array [1/keV]
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: mb   !Mass array [keV]
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: nb   !density array [1/cm^3]
        REAL,                        INTENT(OUT) :: cab_sing

        REAL    :: betae, ne, ge, mabc2, zia, zib, vab2
!
! note: om_b=(1.32155E+3)*SQRT(zb*zb*nb*AMUKEV/mb) ! Plasma frequency [1/s]
!       ge=(6.1260E-15)*SQRT(ne)/te**1.5

        betae=betab(1)
        ne=nb(1)
        ge=GECOEFF*SQRT(ne)*betae**1.5
        zia=zb(ia)
        zib=zb(ib)
        mabc2=mb(ia)*mb(ib)/(mb(ia) + mb(ib)) ! [keV]
        vab2=1./(betab(ia)*mb(ia)) + 1./(betab(ib)*mb(ib)) ! [dimensionless]
        cab_sing=-LOG(0.25*ge*ABS(zia*zib)/(mabc2*betae*vab2)) - 2*GAMMA 
      END SUBROUTINE cab_sing_mass
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! regular contribution for non-zero electron mass: normalization cei_reg(me=0)=-1/2
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
      FUNCTION dcab_reg(x, nni, ia, ib, betab, zb, mb, nb)
      USE physvars
      USE mathvars
      IMPLICIT NONE
      REAL,                        INTENT(IN)  :: x
      INTEGER,                     INTENT(IN)  :: nni    !  Number of ion species
      INTEGER,                     INTENT(IN)  :: ia     !  Species type
      INTEGER,                     INTENT(IN)  :: ib     !  Species type
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: zb     !  Charge array
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: betab  !  Temperature array [1/keV]
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: mb     !  Mass array [keV]
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: nb     !  density array [1/cm^3]
      REAL                                     :: dcab_reg
      REAL,    DIMENSION(1:nni+1) :: betamb, kbar2b
      REAL                        :: fr, fi, fabs, farg
      REAL                        :: mu, mu2, x2, rx, hx, kaic, rxic
!     REAL                        :: abca2, kaia
      INTEGER                     :: ic
!
! construct parameters
!
      mu=SUM(SQRT(betab(2:nni+1)*mb(2:nni+1)))  ! inverse thermal velocity
      mu2=mu*mu                                 ! mu^2
      betamb=SQRT(betab*mb)/mu                  !
      kbar2b=zb*zb*betab*nb/(betab(1)*nb(1))    ! kbar=k_b/k_e
!
! construct H(x)*x from  F_re, F_im, |F|, arg(F)
!
      CALL frfi(x,nni,kbar2b,betamb,fr,fi,fabs,farg)! F(x) 
      hx=2*(fi*LOG(fabs) + fr*farg)*x           ! H(x)*x
!
! construct spectral weight ratio R_ab(x)

!*!
!      rxic=0
!      DO ic=1,nni+1
!         kaic=kbar2b(ic)*betamb(ic)
!         IF (ic == ib) kaia=kaic*EXP(-betamb(ib)*betamb(ib)*x*x)
!         abca2 =betamb(ic)**2 - betamb(ia)**2
!         rxic = rxic + kaic*EXP(-abca2*x*x)
!      ENDDO
!      rx=kaia/rxic
!*!
      x2=x*x
      rxic=0
      kaic=EXP(-(betamb(ia)*betamb(ia)+betamb(ib)*betamb(ib))*x2)
      DO ic=1,nni+1
         rxic=rxic + kbar2b(ic)*betamb(ic)*EXP(-betamb(ic)*betamb(ic)*x2)
      ENDDO
      rx=kaic*kbar2b(ib)*betamb(ib)/rxic
!
! construct un-normalized rate integrand: normalization in cab_reg_mass
!
      dcab_reg=rx*hx
!
      END FUNCTION dcab_reg 
      SUBROUTINE cab_reg_mass(nni, ia, ib, betab, zb, mb, nb, c_r)
      USE mathvars
        IMPLICIT NONE
        INTEGER,                     INTENT(IN)  :: nni    !  Number of ion species
        INTEGER,                     INTENT(IN)  :: ia     !
        INTEGER,                     INTENT(IN)  :: ib     !
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: zb     !  Charge array
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: betab  !  Temperature array [1/keV]
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: mb     !  Mass array [keV]
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: nb     !  density array [1/cm^3]
        REAL,                        INTENT(OUT) :: c_r

        REAL,    PARAMETER :: UPM=0.7745966692E0
        REAL,    PARAMETER :: W13=0.5555555556E0, W2=0.8888888889E0

        REAL    :: y, dcab_reg
        REAL    :: xmin, xmax, x, dx, xm, xc
        REAL    :: cn, ne, betae
        INTEGER :: nmax, ic

        REAL,    DIMENSION(1:nni+1) :: ab
        REAL    :: mu, mu2
!
! integration cutoff determined by thermal velocity of ions

        mu=SUM(SQRT(betab(2:nni+1)*mb(2:nni+1)))  ! inverse thermal velocity
        mu2=mu*mu                                 ! mu^2
        ab=SQRT(betab*mb)/mu                      !
!
        xc=5 ! ti=1, te=0.01, 0.1, 1, 10, 100
        xc=1./MIN(ab(ia),ab(ib))
        xmin=0.
        xmax=5*xc 
        nmax=1000
        dx=(xmax-xmin)/nmax 
        x=xmin-dx
        c_r=0
        DO ic=1,nmax,2
!     
           x=x+2.E0*dx
           y=dcab_reg(x,nni,ia,ib,betab,zb,mb,nb)
           c_r=c_r + W2*y
!
           xm=x-dx*UPM
           y=dcab_reg(xm,nni,ia,ib,betab,zb,mb,nb)
           c_r=c_r + W13*y
!
           xm=x+dx*UPM
           y=dcab_reg(xm,nni,ia,ib,betab,zb,mb,nb)
           c_r=c_r + W13*y
        ENDDO
        c_r=c_r*dx
!
! construct normalization coefficient
!
        ne=nb(1)                                  !
        betae=betab(1)                            !
        cn=0
        DO ic=2,nni+1                                            !
           cn = cn + zb(ic)*zb(ic)*nb(ic)/(mb(ic)*ne*betae)      ! omega_I^2/kappa_e^2
        ENDDO                                                    ! 
        cn=-1./(PI*mu2*cn)                                       ! ke^2/(PI*mu^2 omI^2)
        c_r=c_r*cn 
      END SUBROUTINE cab_reg_mass

!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! quantum contribution for non-zero electron mass
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
      SUBROUTINE cab_qm_mass(nni, ia, ib, betab, zb, mb, qm)
      USE physvars
      USE mathvars
        IMPLICIT NONE
        INTEGER,                     INTENT(IN)  :: nni    !  Number of ion species
        INTEGER,                     INTENT(IN)  :: ia     !  
        INTEGER,                     INTENT(IN)  :: ib     !  
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: zb     !  Charge array
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: betab  !  Temperature array [1/keV]
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: mb     !  Mass array [keV]
        REAL,                        INTENT(OUT) :: qm

        REAL,    PARAMETER :: UPM=0.7745966692E0
        REAL,    PARAMETER :: W13=0.5555555556E0, W2=0.8888888889E0
        INTEGER, PARAMETER :: NMAX=1000
        REAL    :: xmin, xmax, dx, x, xm, y, dcei_qm, eta_ab, vab
        INTEGER :: ix

        vab=SQRT(1./(betab(ia)*mb(ia)) + 1./(betab(ib)*mb(ib)))
        eta_ab=ABS(zb(ia)*zb(ib))*2*BEKEV*A0CM/HBARC/vab

        xmin=0.
        xmax=15
        dx=(xmax-xmin)/NMAX
        x=xmin-dx
        qm=0
        DO ix=1,NMAX,2
!     
           x=x+2*dx
           y=dcei_qm(x, eta_ab)
           qm=qm + W2*y
!
           xm=x-dx*UPM
           y=dcei_qm(xm,eta_ab)
           qm=qm + W13*y
!
           xm=x+dx*UPM
           y=dcei_qm(xm,eta_ab)
           qm=qm + W13*y
        ENDDO
        qm=-0.5*qm*dx
      END SUBROUTINE cab_qm_mass

!
      FUNCTION dcei_qm(x, eta)
      IMPLICIT NONE
      REAL  :: x       ! xi
      REAL  :: eta     ! quantum parameter
      REAL  :: dcei_qm ! quantum integrand
      REAL  :: xh, repsilog
      xh=eta/SQRT(x)
      dcei_qm=EXP(-x/2)*repsilog(xh) 
      END FUNCTION dcei_qm
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Rate coefficient in extreme quantum limit with light electrons: 
! eta --> 0 and me -->0
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! bps_rate_cei_born: returns the rate C_eI for me --> 0 in the extreme 
! quantum limit (the Born approximation to the two-body scattering).
!
! ROUTINE: bps_rate_cei_born(nni, betab, zb, mb, nb, ln_bps_born, cei_born)
! 
! The rate of change of energy density is given by 
!
!  dE_{eI}
!  -------- = -C_{eI} (T_e - T_I) 
!    dt
!                k_e^2    [ betae me]^1/2
! with C_{eI} = ------- * [ --------]      * Omega_I^2 * ln_BPS
!                 2pi     [    2 pi ]
!
!                                            Omega_I^2 = \sum_i omega_i^2
!
!                 1  [   [      8 T_e^2     ]               ]
!       ln_BPS = --- [ ln[ ---------------- ]  - gamma - 1  ]
!                 2  [   [ hbar^2 omega_e^2 ]               ]
!
!
! UNITS: C_{eI} has units of [1/cm^2*s] 
!
! INPUT: nni, betab, zb, mb, nb
!
! plasma input quantities:
! nni  : Number of total plasma species = number ion species + 1
! zb   : Charges of the plasma species. By convention zp(1) is the 
!      : electron plasma component. [dimensionless, Array]
! betab: Inverse temperatures of the plasma components. For an
!        electron-ion plasma, set betab(1)=1/T_e and all other
!        values of the array to 1/T_I [keV^-1].
! mb   : Masses of the plasma species [keV]. 
! nb   : Number densities of the plasma species [cm^-3]. 
!
! OUTPUT: ln_bps_born, cei_born
!
      SUBROUTINE bps_rate_cei_born(nni, betab, zb, mb, nb, ln_bps_born, cei_born)
      USE bpsvars
      USE mathvars
      USE physvars
        IMPLICIT NONE
        INTEGER,                     INTENT(IN)  :: nni     !  Number of ions
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: betab   !  Temperature array    [1/keV]
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: zb      !  Charge array 
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: mb      !  Mass array [keV]
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: nb      !  Number density array [cm^-3]

        REAL                       , INTENT(OUT) :: ln_bps_born! BPS Coulomb log 
        REAL                       , INTENT(OUT) :: cei_born   ! equilibration rate

        REAL,    PARAMETER :: UPM=0.7745966692E0
        REAL,    PARAMETER :: W13=0.5555555556E0, W2=0.8888888889E0

        REAL    :: mi, ni, zi, omi2, ome2, cn
        REAL    :: me, ne, betae, te2, ke2
        INTEGER :: ib
!
! construct plasma quantities: k_e^2, omega_I^2
!
        ne=nb(1)
        me=mb(1)
        betae=betab(1)
        te2 =1/betae**2
        ke2 =8*PI*A0CM*BEKEV*ne*betae
        ome2=8*PI*A0CM*BEKEV*ne*CC2/me
        omi2=0
        DO ib=2,nni+1
           ni=nb(ib)
           mi=mb(ib)
           zi=zb(ib)
           omi2=omi2 + 8*PI*A0CM*BEKEV*ni*zi*zi*CC2/mi
        ENDDO
!
!                           om_I^2*k_e^2    [ beta_e m_e*c^2 ]^1/2    1
! construct prefactor: cn = ------------  * |--------------- |     * ---
!                               2 pi        [   2 pi         ]        c
!
        cn=omi2*ke2/TWOPI
        cn=cn*SQRT(betae*me/TWOPI)/CC
!
!                                          1   [   [    8 T_e^2    ]             ]
! construct BPS Coulomb log: ln_bps_born = --- [ ln[ ------------- ] - GAMMA -1  ]
!                                          2   [   [  hbar^2 ome^2 ]             ]

        ln_bps_born=0.5*(LOG(8*te2*CC2/(HBARC2*ome2)) -GAMMA - 1)
!
! construct rate: cei_born
!
        cei_born=cn*ln_bps_born
      END SUBROUTINE bps_rate_cei_born
!
!=============== extraneous subroutines [we may need them] ====================
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Old code segments for bps_cei_* [sum over i]
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Returns the rate coefficient C_{e ib, sing}
!
      SUBROUTINE bps_cei_sing(nni, ib, betab, zb, mb, nb, cei_sing)
      USE mathvars
      USE physvars 
       IMPLICIT NONE
        INTEGER,                     INTENT(IN):: nni  !Number of ion species
        INTEGER,                     INTENT(IN):: ib
        REAL,    DIMENSION(1:nni+1), INTENT(IN):: zb   !Charge array
        REAL,    DIMENSION(1:nni+1), INTENT(IN):: betab!Temperature array [1/keV]
        REAL,    DIMENSION(1:nni+1), INTENT(IN):: mb   !Mass array [keV]
        REAL,    DIMENSION(1:nni+1), INTENT(IN):: nb   !density array [1/cm^3]
        REAL,                        INTENT(OUT) :: cei_sing
        REAL    :: betae, me, ne, zi, ni, mi, betai
        REAL    :: mui, deli, ge, oib2, sing
!
! note: om_b=(1.32155E+3)*SQRT(zb*zb*nb*AMUKEV/mb) ! Plasma frequency [1/s]
!       ge=(6.1260E-15)*SQRT(ne)/te**1.5

        betae=betab(1)
        me=mb(1)
        ne=nb(1)

        cei_sing=0
        ge=GECOEFF*SQRT(ne)*betae**1.5
        zi=zb(ib)
        ni=nb(ib)
        mi=mb(ib)
        betai=betab(ib)

        mui =me/mi                  ! mu_i
        deli=mui*betae/betai        ! del_i
        mui =1 + mui                ! 1 + mu_i
        deli=1 + deli               ! 1 + delta_i

        oib2=8*PI*A0CM*CC2*ni*zi*zi*BEKEV/mi
        sing=-LOG(0.25*ge*ABS(zi)*mui/deli) - 2*GAMMA 
        sing=oib2*sing/deli**1.5
        cei_sing = cei_sing + sing
      END SUBROUTINE bps_cei_sing
!
!Returns the rate coefficient C_{e ib, reg}
!
      SUBROUTINE bps_cei_reg(nni, ib, betab, zb, mb, nb, c_r)
        IMPLICIT NONE
        INTEGER,                     INTENT(IN)  :: nni    !  Number of ion species
        INTEGER,                     INTENT(IN)  :: ib     !
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: zb     !  Charge array
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: betab  !  Temperature array [1/keV]
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: mb     !  Mass array [keV]
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: nb     !  density array [1/cm^3]
        REAL,                        INTENT(OUT) :: c_r

        REAL,    PARAMETER :: UPM=0.7745966692E0
        REAL,    PARAMETER :: W13=0.5555555556E0, W2=0.8888888889E0

        REAL    :: y, dcei_reg
        REAL    :: xmin, xmax, x, dx, xm, xc
        INTEGER :: nmax, ic
!
! integration cutoff determined by thermal velocity of ions
!
        xc=2. ! ti=1, te=0.01, 0.1, 1, 10, 100

        xmin=0.
        xmax=5*xc ! automate this choice later.
        nmax=1000
        dx=(xmax-xmin)/nmax 
        x=xmin-dx
        c_r=0
        DO ic=1,nmax,2
!     
           x=x+2.E0*dx
           y=dcei_reg(x,nni,ib,betab,zb,mb,nb)
           c_r=c_r + W2*y
!
           xm=x-dx*UPM
           y=dcei_reg(xm,nni,ib,betab,zb,mb,nb)
           c_r=c_r + W13*y
!
           xm=x+dx*UPM
           y=dcei_reg(xm,nni,ib,betab,zb,mb,nb)
           c_r=c_r + W13*y
        ENDDO
        c_r=c_r*dx
      END SUBROUTINE bps_cei_reg

      FUNCTION dcei_reg(x, nni, ib, betab, zb, mb, nb)
      USE physvars
      USE mathvars
      IMPLICIT NONE
      REAL,                        INTENT(IN)  :: x
      INTEGER,                     INTENT(IN)  :: nni    !  Number of ion species
      INTEGER,                     INTENT(IN)  :: ib     !  Species type
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: zb     !  Charge array
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: betab  !  Temperature array [1/keV]
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: mb     !  Mass array [keV]
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: nb     !  density array [1/cm^3]
      REAL                                     :: dcei_reg
      REAL,    DIMENSION(1:nni+1) :: ab, kbar2b
      REAL                        :: fr, fi, fabs, farg
      REAL                        :: mu, mu2, rx, hx, ae, ne, betae
      REAL                        :: kaic, kaib, rxic, abe2, cn
      INTEGER                     :: ic
!
! construct parameters
!
      ne=nb(1)                                  !
      betae=betab(1)                            !
      kbar2b=zb*zb*betab*nb/(betae*ne)          ! kbar=k_b/k_e
      mu=SUM(SQRT(betab(2:nni+1)*mb(2:nni+1)))  ! inverse thermal velocity
      mu2=mu*mu                                 ! mu^2
      ab=SQRT(betab*mb)/mu                      ! a_b=beta_b m_b/sum_i beta_i m_i
      ae=ab(1)                                  !
!
! construct H(x)*x from  F_re, F_im, |F|, arg(F)
!
      CALL frfi(x,nni,kbar2b,ab,fr,fi,fabs,farg)! F(x) 
      hx=2*(fi*LOG(fabs) + fr*farg)*x ! H(x)*x
!
! construct spectral weight ratio R_b(x)
! R_b(x)=exp(- ab^2 x^2) rho_b(x)/rho_tot(x)
!
      rxic=0
      DO ic=1,nni+1
         kaic=kbar2b(ic)*ab(ic)
         IF (ic == ib) kaib=kaic*EXP(-ab(ib)*ab(ib)*x*x)
         abe2 =ab(ic)**2 - ae**2
         rxic = rxic + kaic*EXP(-abe2*x*x)
      ENDDO
      rx=kaib/rxic
!
! construct un-normalized rate integrand
!
      dcei_reg=rx*hx
!
! construct normalization coefficient
!
      cn=0
      DO ic=2,nni+1                                            !
         cn = cn + zb(ic)*zb(ic)*nb(ic)/(mb(ic)*ne*betae)      ! omega_I^2/kappa_e^2
      ENDDO                                                    ! 
      cn=-1./(PI*mu2*cn)                                       ! ke^2/(PI*mu^2 omI^2)
!
      dcei_reg=dcei_reg*cn 
      END FUNCTION dcei_reg


      SUBROUTINE bps_cei_qm(nni, ib, betab, zb, mb, nb, qm)
      USE physvars
      USE mathvars
      IMPLICIT NONE
      INTEGER,                     INTENT(IN)  :: nni    !  Number of ion species
      INTEGER,                     INTENT(IN)  :: ib     !  
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: zb     !  Charge array
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: betab  !  Temperature array [1/keV]
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: mb     !  Mass array [keV]
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: nb     !  density array [1/cm^3]
      REAL,                        INTENT(OUT) :: qm

      REAL,    PARAMETER :: UPM=0.7745966692E0
      REAL,    PARAMETER :: W13=0.5555555556E0, W2=0.8888888889E0
      INTEGER, PARAMETER :: NMAX=1000

      REAL    :: betae, betai, me, mi, ne, ni, zi
      REAL    :: deli, deli3, mui, oi2, eta_ei
      REAL    :: xmin, xmax, dx, x, xm, y, dcei_qm
      REAL    :: vei, ve2, vi2
      INTEGER :: ix
      qm=0
      betae=betab(1)
      me=mb(1)
      ne=nb(1)
      ve2=1./(me*betae)
         zi    =zb(ib)
         ni    =nb(ib)
         mi    =mb(ib)
         betai =betab(ib)
         vi2   =1./(mi*betai)
         vei   =CC*SQRT(ve2 + vi2)
         eta_ei=2*BEKEV*A0CM*CC*ABS(zi)/(HBARC*vei)
!        oi2   =OMEGI2*zi*zi*ni*AMUKEV/mi 
         oi2   =8*PI*A0CM*BEKEV*ABS(zi)*ni*CC2/mi

         mui  =me/mi                  ! mu_i
         deli =mui*betae/betai        ! del_i
         mui  =1 + mui                ! 1 + mu_i
         deli =1 + deli               ! 1 + delta_i
         deli3=deli**1.5

         xmin=0.
         xmax=15
         dx=(xmax-xmin)/NMAX
         x=xmin-dx
         qm=0
         DO ix=1,NMAX,2
!     
            x=x+2*dx
            y=dcei_qm(x, eta_ei)
            qm=qm + W2*y
!
            xm=x-dx*UPM
            y=dcei_qm(xm,eta_ei)
            qm=qm + W13*y
!
            xm=x+dx*UPM
            y=dcei_qm(xm,eta_ei)
            qm=qm + W13*y
         ENDDO
         qm=-0.5*qm*dx   
         qm=oi2*qm/deli3
      END SUBROUTINE bps_cei_qm

!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! bps_rate_el: Calculates the *classical* electron-ion temperature equilibration 
! rate for the physical mass of the electron (i.e. heavy ions and a light electron). 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
      SUBROUTINE bps_rate_el(nni, betab, zb, mb, nb, clog, cei)
      USE bpsvars
      USE mathvars
      USE physvars
      IMPLICIT NONE
      INTEGER,                     INTENT(IN)  :: nni    !  Number of ions
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: betab  !  Temperature array    [1/keV]
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: zb     !  Charge array 
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: mb     !  Mass array [keV]
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: nb     !  Number density array [cm^-3]

      REAL                       , INTENT(OUT) :: clog   ! BPS Coulomb log 
      REAL                       , INTENT(OUT) :: cei    ! equilibration rate

!     REAL,    PARAMETER :: DEBYE=4.2539E-5      ! Used for Debye wave number
!     REAL,    PARAMETER :: OMEGI=1.32155E+3     ! Used For plasma frequency 
!     REAL,    PARAMETER :: AMUKEV=0.931494012E6 ! 1 AMU [keV]


      REAL    :: ge, omi2, oi2, ke2, cn
      INTEGER :: ib
!
! construct om_I^2
!
      omi2=0
      DO ib=2,nni+1 ! Plasma frequency om_I^2 [1/s^2]
         omi2=omi2 + OMEGI*OMEGI*zb(ib)*zb(ib)*nb(ib)*AMUKEV/mb(ib) 
      ENDDO
!
! calculate cei
!
!     clog=LOG(4./(zb(2)*ge))-2*GAMMA-0.5

      ge=(6.1260E-15)*SQRT(nb(1))*betab(1)**1.5
      cei=0
      DO ib=2,nni+1 ! Plasma frequency om_i^2 [1/s^2]
         oi2=OMEGI*OMEGI*zb(ib)*zb(ib)*nb(ib)*AMUKEV/mb(ib) 
         cei=cei - oi2*LOG(0.25*zb(ib)*ge)
      ENDDO
      clog=cei/omi2 - (2*GAMMA + 0.5)
      cei=cei - omi2*(2*GAMMA + 0.5)
!
! construct prefactor
!
      ke2=DEBYE*DEBYE*zb(1)*zb(1)*nb(1)*betab(1) ! Electron Debye wave number [1/cm^2]
      cn=ke2/(2*PI)
      cn=cn*SQRT(betab(1)*mb(1)/(2*PI))/CC
      cei=cei*cn
      END SUBROUTINE bps_rate_el

!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! single ion species for C_{eI, reg}. Used for debugging full bps_cei_reg_mass()
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
      SUBROUTINE bps_cei_reg_ep(nni, zb, betab, mb, y_ep) ! returns \Delta for Cei_R
      IMPLICIT NONE
      INTEGER,                     INTENT(IN)  :: nni    !  Number of ion species
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: zb     !  Charge array
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: betab  !  Temperature array [1/keV]
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: mb     !  Mass array [keV]
      REAL,                        INTENT(OUT) :: y_ep
      REAL,    PARAMETER :: DEBYE=4.2539E-5   ! Used for Debye wave number
      REAL,    PARAMETER :: OMEGI=1.32155E+3  ! Used For plasma frequency 
      REAL,    DIMENSION(:), ALLOCATABLE :: ob2, kb2

      REAL,    PARAMETER :: AMUKEV=0.931494012E6 ! 1 AMU [keV]
      REAL,    PARAMETER :: PI   =3.141592654    ! pi
      REAL,    PARAMETER :: CC=2.998E10          ! speed of light [cm/s]

      REAL,    PARAMETER :: UPM=0.7745966692E0
      REAL,    PARAMETER :: W13=0.5555555556E0, W2=0.8888888889E0

      REAL    :: xmin, xmax, x, dx, xm
      INTEGER :: nmax, ix
      REAL    :: hi_ep, te, ti, alpha, delta

      ALLOCATE(ob2(1:nni+1),kb2(1:nni+1))  
      xmin=0.
      xmax=10.
      nmax=1000
      dx=(xmax-xmin)/nmax 
      te=1./betab(1)
      ti=1./betab(2)
      alpha=zb(2)**2*te/ti
      delta=SQRT(betab(1)*mb(1))/SQRT(betab(2)*mb(2))
      y_ep=0
      x=xmin-dx
      DO ix=1,nmax,2
!     
         x=x+2.E0*dx
         CALL di2_ep(x,alpha,delta,hi_ep)
         y_ep=y_ep+W2*hi_ep
!
         xm=x-dx*UPM
         CALL di2_ep(xm,alpha,delta,hi_ep)
         y_ep=y_ep+W13*hi_ep
!
         xm=x+dx*UPM
         CALL di2_ep(xm,alpha,delta,hi_ep)
         y_ep=y_ep+W13*hi_ep
      ENDDO
      y_ep=y_ep*dx
      END SUBROUTINE bps_cei_reg_ep

      SUBROUTINE di2_ep(x, alpha, delta, dint2)
      IMPLICIT NONE
      REAL,    INTENT(IN) :: x, alpha, delta
      REAL,    INTENT(OUT):: dint2
      REAL                :: delta1, x2
      REAL                :: fr_ep, fi_ep, fabs_ep, farg_ep
      REAL,    PARAMETER  :: PI   =3.141592654    ! pi
      CALL frfi_ep(x,alpha,delta,fr_ep,fi_ep,fabs_ep,farg_ep)      
      dint2=2*(fi_ep*LOG(fabs_ep) + fr_ep*farg_ep)*x
      x2=x*x
      dint2=dint2*EXP(-x2)
      delta1=1-delta*delta
      dint2=dint2/(delta + alpha*EXP(-delta1*x2))
      dint2=dint2*(2./PI)
      END SUBROUTINE di2_ep

      SUBROUTINE frfi_ep(x, alpha, delta, fr, fi, fabs, farg)
      IMPLICIT NONE
      REAL,    INTENT(IN) :: x, alpha, delta
      REAL                :: fr, fi, fabs, farg
      REAL                :: daw, y
      REAL,    PARAMETER  :: SQPI =1.772453851    ! sqrt(pi)
      REAL,    PARAMETER  :: PI   =3.141592654    ! pi
      y=x*delta
      fr=1 - 2*y*daw(y)            ! electrons
      fr=fr + alpha*(1-2*x*daw(x)) ! plus ions 
      fi=y*Exp(-y*y)               ! electrons
      fi=fi + alpha*x*EXP(-x*x)    ! plus ions
      fi=fi*SQRT(PI)
      fabs=SQRT(fr*fr + fi*fi)
      farg=ATAN2(fi,fr)
      END SUBROUTINE frfi_ep
