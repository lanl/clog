! 
! ROUTINE: bps_acoeff_ab_mass(nni, ep, mp, zp, ia, ib, betab, zb, mb, nb, &
!            a_ab, a_ab_sing, a_ab_reg, a_ab_qm)
! 
! Assume a plasma composed of several species b, each separately in 
! thermal equilibrium with themselves but not necessarily with each 
! other[1]. This routine returns several useful components of the 
! corresponding A-coefficients introduced in Note [2] below (BPS).
! 
! UNITS: A_{pb} has units of [MeV/micron] (subject to change in updates)
! 
! THE PHYSICS:
! The various subsystems b will exchange coulomb energy and they will
! eventually equilibrate to a common temperature. The A-coefficients 
! introduced in Ref. [2] encode this coulomb energy exchange, exactly to 
! leading and next-to-leading orders in the plasma coupling constant g.
! See Refs. [3,4,5] for more details. For a weakly coupled plasma (g << 1), 
! the BPS calculation is essentially exact, and the error is O(g). Physical 
! properties of interest, such as the stopping power dE/dx and the temperature 
! equilibration rate between plasma species, can be obtained directly from 
! the A-coefficients. 
!
! USAGE:
! Since electrons are thousands of times lighter than ions, one of the most
! physically accessible regime is the in which the electrons have a temperature 
! T_e and the ions have a (possibly different) common temperature T_I. This 
! is why the output is organized into electron contributions and total ion 
! contributions (sum over all ions).
!
! INPUT: nni, ep, zp, mp, ia, ib, betab, zb, mb, nb
!
! Describe the incident projectile and the background plasma. 
!
! projectile input quantities:
! ep : classical kinetic energy of the projectile [keV]
! zp : charge of the projectile in units of Z_p [dimensionless]
! mp : mass of the projectile [keV], i.e. mp = mp[grams]*c^2
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
! ia   : First plasma species [usually the projectile]
! ib   : Second plasma species
!
! OUTPUT: a_ab, a_ab_sing, a_ab_reg, a_ab_qm
!
! Each plasma component b makes a linear contribution A_b to the total 
! A-coefficient, i.e. A = sum_b A_b [5]. Each A_b in turn can be be 
! decomposed into a classical-quantum or electron-ion contributions.
!
! classical electron  : ac_e
! classical ion       : ac_i [sum over all ions]
! classical total     : ac_tot = ac_e + ac_i 
! quantum   electron  : aq_e
! quantum   ion       : aq_i [sum over all ions]
! quantum   total     : a_tot = aq_e + aq_i
! total     electric  : a_e = ac_e + aq_e
! total     ion       : a_i = ac_i + aq_i
! total               ! a_tot = a_e + a_i
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
! [3] The code employs rationalized cgs units in which the dimensionless
!     plasma coupling parameter is defined by g = e^2 kappa/(4Pi*T); in 
!     these units the Debye wavenumber is determined by kappa^2 = e^2 n/T 
!     and the plasma frequency by omega^2 = e^2 n/m. A weakly coupled 
!     plasma is one for which g << 1, i.e. a plasma with thermal kinetic 
!     energy (of order the temperature T) dominates the coulomb potential 
!     energy (for particles separated by a Debye length). In the more 
!     common non-rationalized cgs units, we define g = e^2 kappa/T, with 
!     kappa^2 = 4 Pi e^2 n/T and omega^2 = 4 Pi e^2 n/m. 
!
! [4] For coulomb energy exchange processes, the leading and next-to-
!     leading order terms in the plasma coupling g are proportional to 
!     g^2*ln(g) and g^2, respectively. That is to say, for a property 
!     denoted by F, one can expand F in powers of g in the form:
!
!       F(g,eta) = A(eta)*g^2*ln(g) + B(eta)*g^2 + O(g^3) 
!                 = A(eta)*g^2*[ ln(C(eta)*g)+O(g) ],
!
!     where eta is the dimensionless quantum parameter (eta <<1 means
!     extreme classical scattering while eta >>1 means the extreme quantum 
!     limit). The relative error of BPS is therefore O(g). At the center of 
!     the sun g=0.4, and so the error of Ref. [1] is only of order 4% in 
!     this case. For the processes of charged stopping power and electron-
!     ion temperature equilibration, Ref. [1] calculates the corresponding
!     functions A(eta) and B(eta) exactly, including all orders in the 
!     two-body quantum scattering parameter eta = e^2/4Pi*hbar*v_thermal 
!     (this means that BPS gives the correct interpolation between the 
!     classical and quantum regimes, exact to leading and next-to-leading 
!     order). The O(g^3) terms physically correspond to 3-body correlations 
!     within the plasma, and for a sufficiently weak plasma these are 
!     negligible. For strongly coupled plasmas (g >> 1), all terms in a 
!     g-expansion are important and the BPS calculation is not applicable.
!
! [5] It makes sense to talk about separate linear contribution A_b
!     contributing *from* a given plasma component b only for a weakly 
!     coupled plasma. More exactly, A=sum_b A_b holds true only up to 
!     leading and next-to-leading order in the plasma coupling g. This is 
!     the order to which Ref. [2] calculates all quantities, and therefore 
!     BPS works to a consistent order in g.
!
!
! vp = projectile speed [cm/s]
!
!            [ 2 Ep ]
!    = c sqrt[ ---- ]  where mp and Ep are in keV and
!            [ mpc2 ]  c is the speed of light in cm/s
!
!
!        ep^2 kb^2
! c1  =  --------- = 2 zp**2 * Be * kb^2 * a0  [keV/cm]
!           4 Pi                               [MeV/micron]
!
!       [ betab mbc2 ]^1/2  
! c2 =  |----------- |      vp   [dimensionless]
!       [   2 Pi     ]      
!        
! where
!          e^2 
! Be  =  --------- = 13.606E-3   [keV]
!        8 Pi a0                 Bohr radius: a0=5.29E9 cm
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! main driver for A-coefficient for general quantum and electron-mass regimes
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
      SUBROUTINE bps_acoeff_ab_mass(nni, ep, mp, zp, ia, ib, betab, zb, mb, nb, &
            a_ab, a_ab_sing, a_ab_reg, a_ab_qm)
      USE physvars
      USE mathvars    
        IMPLICIT NONE                                             ! Plasma:
        INTEGER,                            INTENT(IN)  :: nni    !  number of ions
        REAL,                               INTENT(IN)  :: ep     !  energy input [keV]
        REAL,                               INTENT(IN)  :: mp     !  mass [keV]
        REAL,                               INTENT(IN)  :: zp     !  charge
        INTEGER,                            INTENT(IN)  :: ia     !  
        INTEGER,                            INTENT(IN)  :: ib     !  
        REAL,    DIMENSION(1:nni+1),        INTENT(IN)  :: betab  !  temp array [1/keV]
        REAL,    DIMENSION(1:nni+1),        INTENT(IN)  :: mb     !  mass array [keV]
        REAL,    DIMENSION(1:nni+1),        INTENT(IN)  :: nb     !  density [1/cc]
        REAL,    DIMENSION(1:nni+1),        INTENT(IN)  :: zb     !  charge array
                                                                  !
                                                                  ! A-coeffs [MeV/micron]
        REAL,                               INTENT(OUT) :: a_ab
        REAL,                               INTENT(OUT) :: a_ab_sing
        REAL,                               INTENT(OUT) :: a_ab_reg
        REAL,                               INTENT(OUT) :: a_ab_qm

        REAL,    DIMENSION(1:nni+1)  :: mpb, mbpb, kb2, ab
        REAL                         :: vp, zp2, k, k2, kd, kd2, a, b, eta
        REAL                         :: ac_r, ac_s, aq, c1, c2

        REAL, PARAMETER              :: EPS_SMALL_E=2.E-4
        REAL, PARAMETER              :: EPS_SMALL_E_SING=2.E-4
        REAL, PARAMETER              :: EPS_SMALL_E_REG=2.E-4
!
! initialize components of A-coefficients
!
        kb2=8*PI*A0CM*BEKEV*zb*zb*nb*betab
        kd2 = SUM(kb2)                ! [1/cm^2]
        kd  = SQRT(kd2)               ! [1/cm]
        k2  = kb2(1)                  ! [1/cm^2]
        k   = SQRT(k2)                ! [1/cm]   k = k_e
!
! Loop over charged plasma species
!
        mpb = mp*mb/(mp+mb)            ! [keV]
        mbpb= mb/mpb                   ! [dimensionless]
        vp =CC*SQRT(2*ep/mp)           ! [cm/s]
        zp2=zp**2                      ! [dimensionless]
                                       ! ab=(1/2) betab(ib)*mbc2(ib)*vp2/CC2
        ab  =0.5*betab*mb*vp*vp/CC2    ! [dimensionless] 
        IF (zb(ib) .NE. 0.) THEN
        a  =ab(ib)
        b  =-Log(2*betab(ib)*BEKEV*ABS(zp*zb(ib))*k*A0CM*mbpb(ib) )-2*GAMMA+2
        eta=ABS(zp*zb(ib))*2.1870E8/vp ! defined with projectile velocity vp
        c1=2*zp2*BEKEV*kb2(ib)*A0CM    ! [keV/cm] c1 = e_p^2 kappa_b^2/(4 Pi)
        c1=c1*1.E-7                    ! [MeV/micron]  
        c2=SQRT(a/PI)                  ! [dimensionless] 
                                       ! c2=SQRT(betab(ib)*mb(ib)/TWOPI)*vp/CC 
!
! A_{ab}-classical-singular 
!
        CALL a_sing_mass(a,b,ac_s) 
        a_ab_sing=c1*c2*ac_s
!
! A_{ab}-classical-regular 
!
        CALL a_reg_mass(nni,ia,ib,vp,k2,kb2,betab,mb,ac_r)
        a_ab_reg=c1*ac_r
!
! A_{ab}-quantum
!
        CALL a_quantum_mass(ia,ib,a,eta,aq) ! eta = dimensionless quantum param.
        a_ab_qm=c1*c2*aq
!
! A_{ab}-total
!
        a_ab=a_ab_sing + a_ab_reg + a_ab_qm
        ENDIF
      END SUBROUTINE bps_acoeff_ab_mass
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Assembles the matrix A_{ab} of the A-coefficients.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
      SUBROUTINE bps_acoeff_ab_matrix(nni, ep, betab, zb, mb, nb,    &
        a_ab, a_ab_sing, a_ab_reg, a_ab_qm, a_tot, a_i, a_e, ac_tot, &
        ac_i, ac_e, aq_tot, aq_i, aq_e, ac_s_i, ac_s_e, ac_r_i, ac_r_e)
      USE physvars
      USE mathvars    
        IMPLICIT NONE                                             ! Plasma:
        INTEGER,                            INTENT(IN)  :: nni    !  number of ions
        REAL,                               INTENT(IN)  :: ep     !  energy
        REAL,    DIMENSION(1:nni+1),        INTENT(IN)  :: betab  !  temp array [1/keV]
        REAL,    DIMENSION(1:nni+1),        INTENT(IN)  :: zb     !  charge array
        REAL,    DIMENSION(1:nni+1),        INTENT(IN)  :: mb     !  mass array [keV]
        REAL,    DIMENSION(1:nni+1),        INTENT(IN)  :: nb     !  density [1/cc]
                                                                  !
                                                                  ! A-coeffs [MeV/micron]
        REAL,    DIMENSION(1:nni+1,1:nni+1),INTENT(OUT) :: a_ab
        REAL,    DIMENSION(1:nni+1,1:nni+1),INTENT(OUT) :: a_ab_sing
        REAL,    DIMENSION(1:nni+1,1:nni+1),INTENT(OUT) :: a_ab_reg
        REAL,    DIMENSION(1:nni+1,1:nni+1),INTENT(OUT) :: a_ab_qm
        REAL,    DIMENSION(1:nni+1),        INTENT(OUT) :: a_tot
        REAL,    DIMENSION(1:nni+1),        INTENT(OUT) :: a_i
        REAL,    DIMENSION(1:nni+1),        INTENT(OUT) :: a_e
        REAL,    DIMENSION(1:nni+1),        INTENT(OUT) :: ac_tot
        REAL,    DIMENSION(1:nni+1),        INTENT(OUT) :: ac_i
        REAL,    DIMENSION(1:nni+1),        INTENT(OUT) :: ac_e
        REAL,    DIMENSION(1:nni+1),        INTENT(OUT) :: aq_tot
        REAL,    DIMENSION(1:nni+1),        INTENT(OUT) :: aq_i
        REAL,    DIMENSION(1:nni+1),        INTENT(OUT) :: aq_e
        REAL,    DIMENSION(1:nni+1),        INTENT(OUT) :: ac_s_i
        REAL,    DIMENSION(1:nni+1),        INTENT(OUT) :: ac_s_e
        REAL,    DIMENSION(1:nni+1),        INTENT(OUT) :: ac_r_i
        REAL,    DIMENSION(1:nni+1),        INTENT(OUT) :: ac_r_e

        REAL    :: aab, aab_sing, aab_reg, aab_qm
        REAL    :: mp, zp
        INTEGER :: ia, ib

        a_i   = 0
        ac_s_i= 0
        ac_r_i= 0
        ac_i  = 0
        aq_i  = 0
        DO ia=1,nni+1
          mp=mb(ia)
          zp=zb(ia)
          DO ib=1,nni+1
            CALL bps_acoeff_ab_mass(nni, ep, mp, zp, ia, ib, betab, zb, mb, nb, &
            aab, aab_sing, aab_reg, aab_qm)
            a_ab(ia,ib)     =aab
            a_ab_sing(ia,ib)=aab_sing
            a_ab_reg(ia,ib) =aab_reg
            a_ab_qm(ia,ib)  =aab_qm 
            IF (ib == 1) THEN
               a_e(ia)   = aab 
               ac_s_e(ia)= aab_sing
               ac_r_e(ia)= aab_reg
               ac_e(ia)  = aab_sing + aab_reg
               aq_e(ia)  = aab_qm
            ELSE
               a_i(ia)   = a_i(ia)    + aab 
               ac_s_i(ia)= ac_s_i(ia) + aab_sing
               ac_r_i(ia)= ac_r_i(ia) + aab_reg
               ac_i(ia)  = ac_i(ia)   + aab_sing + aab_reg
               aq_i(ia)  = aq_i(ia)   + aab_qm
            ENDIF
          ENDDO
          a_tot(ia) = a_e(ia)  + a_i(ia)
          ac_tot(ia)= ac_e(ia) + ac_i(ia)
          aq_tot(ia)= aq_e(ia) + aq_i(ia)
        ENDDO
      END SUBROUTINE bps_acoeff_ab_matrix
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Returns A_{p I} = \sum_i A_{p i} for backward compatibility
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
      SUBROUTINE bps_acoeff_ei_mass(nni, ep, zp, mp, betab, zb, mb, nb, &
            a_tot, a_i, a_e, ac_tot, ac_i, ac_e, aq_tot, aq_i, aq_e,&
            ac_s_i, ac_s_e, ac_r_i, ac_r_e)
      USE physvars
      USE mathvars    
      USE controlvars  
        IMPLICIT NONE                                      ! Plasma:
        INTEGER,                     INTENT(IN)  :: nni    !  number of ions
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: betab  !  temp array [1/keV]
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: mb     !  mass array [keV]
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: nb     !  density [1/cc]
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: zb     !  charge array
                                                           !
                                                           ! Projectile  
        REAL,                        INTENT(IN)  :: ep     !  projectile energy [keV]
        REAL,                        INTENT(IN)  :: mp     !  projectile mass   [keV]
        REAL,                        INTENT(IN)  :: zp     !  projectile charge
                                                           !
                                                           ! A-coeffs [MeV/micron]
        REAL,                        INTENT(OUT) :: a_tot  !  electron + ion
        REAL,                        INTENT(OUT) :: a_i    !  ion contribution
        REAL,                        INTENT(OUT) :: a_e    !  electron contribution
        REAL,                        INTENT(OUT) :: ac_tot !  classical
        REAL,                        INTENT(OUT) :: ac_i   !  classical
        REAL,                        INTENT(OUT) :: ac_e   !  classical
        REAL,                        INTENT(OUT) :: aq_tot !  quantum
        REAL,                        INTENT(OUT) :: aq_i   !  quantum
        REAL,                        INTENT(OUT) :: aq_e   !  quantum
        REAL,                        INTENT(OUT) :: ac_s_i
        REAL,                        INTENT(OUT) :: ac_s_e 
        REAL,                        INTENT(OUT) :: ac_r_i
        REAL,                        INTENT(OUT) :: ac_r_e

        REAL     :: adum, ac_s, ac_r, aq
        INTEGER  :: ia, ib, nnb
!
! initialize components of A-coefficients
!
        a_tot =0  ! electron + ion
        a_i   =0  ! ion contribution
        a_e   =0  ! electron contribution
        ac_tot=0  ! classical total
        ac_e  =0  ! classical electron
        ac_i  =0  ! classical ion
        aq_tot=0  ! quantum total
        aq_e  =0  ! quantum electron
        aq_i  =0  ! quantum ion
        ac_s_i=0 
        ac_s_e=0 
        ac_r_i=0
        ac_r_e=0

        NNB = nni+1                 ! number of ions + electrons
        ia=1
        DO ib=1,nni+1
        IF (zb(ib) .NE. 0.) THEN
            CALL bps_acoeff_ab_mass(nni, ep, mp, zp, ia, ib, betab, zb, mb, nb, &
            adum, ac_s, ac_r, aq)
            CALL x_collect(ib, NNB, ac_s, ac_r, aq,       &
            a_tot, a_i, a_e, ac_tot, ac_i, ac_e, aq_tot,  &
            aq_i, aq_e, ac_s_i, ac_s_e, ac_r_i, ac_r_e)
        ENDIF
        ENDDO
      END SUBROUTINE bps_acoeff_ei_mass
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! singular contribution for non-zero electron mass
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
      SUBROUTINE a_sing_mass(a, b, ac_s)
        REAL,    INTENT(IN)  :: a
        REAL,    INTENT(IN)  :: b
        REAL,    INTENT(OUT) :: ac_s
        REAL                 :: u0, u1, du, u, um
        INTEGER, PARAMETER :: NS=1000 ! integration regions: must be even
        REAL,    PARAMETER :: UPM=0.7745966692E0 ! parameters for Gaussian Quad
        REAL,    PARAMETER :: W13=0.5555555556E0, W2=0.8888888889E0
           ac_s=0
           u0=0
           u1=1
           du=(u1-u0)/NS
           u=u0-du
           DO iu=1,NS,2 ! Gaussian quadrature
              u=u+2.E0*du
              ac_s=ac_s+W2*dab_sing(u,a,b)
              um=u-du*UPM
              ac_s=ac_s+W13*dab_sing(um,a,b)
              um=u+du*UPM
              ac_s=ac_s+W13*dab_sing(um,a,b)
           ENDDO
           ac_s=ac_s*du
      END SUBROUTINE a_sing_mass
!
      FUNCTION dab_sing(u, a, b)
        IMPLICIT NONE
        REAL,        INTENT(IN)  :: u        ! [dimensionless]
        REAL,        INTENT(IN)  :: a        ! [dimensionless] 
                                             ! a=(1/2)*beta*mpc2*vp^2/C^2
        REAL,        INTENT(IN)  :: b        ! [dimensionless]
        REAL                     :: dab_sing ! [dimensionless]
        dab_sing=SQRT(u)*EXP(-a*u)*(-LOG(u/(1-u)) + b)
      END FUNCTION dab_sing
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! regular contribution for non-zero electron mass
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
      FUNCTION dab_reg(u, vp, ia, ib, nni, k2, kb2, betab, mb)
      USE mathvars
      USE physvars
        IMPLICIT NONE
        REAL,                        INTENT(IN)  :: u      ! [dimensionless]
        REAL,                        INTENT(IN)  :: vp     ! Projectile velocity [cm/s]
        INTEGER,                     INTENT(IN)  :: ia     ! Species number
        INTEGER,                     INTENT(IN)  :: ib     ! Species number
        INTEGER,                     INTENT(IN)  :: nni    ! Number of ion species
        REAL,                        INTENT(IN)  :: k2     ! Wave-number squared [1/cm^2]
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: kb2    ! Debye wavenumber squared [1/cm^2]
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: betab  ! Temperature array [1/keV]
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: mb     ! Mass array [keV]
        REAL                                     :: dab_reg! [dimensionless]
        REAL,    DIMENSION(1:nni+1) :: kbar2b, ab, ab2
        REAL                        :: fr, fi, fabs, farg, h, r_ib
        REAL                        :: kcb, bm_ic, bm_ib, a_ic, a_ib, ex, au
        INTEGER                     :: ic
        ab=SQRT(0.5*betab*mb)*vp/CC
        ab2=ab*ab
        kbar2b=kb2/k2
        CALL frfi(u,nni,kbar2b,ab,fr,fi,fabs,farg)
        h=2*(fr*farg + fi*LOG(fabs))*u
!*!
!        h=2*(fr*farg + fi*LOG(fabs))
!*!
!
! construct spectral weight ratio Rb=rho_b/rho_tot
!
!*!
!         r_ib=0
!         DO ic=1,nni+1
!           r_ib=r_ib + kbar2b(ib)*(ab(ic)/ab(ib))*EXP((ab2(ib)-ab2(ic))*u*u)
!         ENDDO
!         r_ib=1./r_ib
!*!
        r_ib=0
        bm_ib=betab(ib)*mb(ib)
        a_ib =ab(ib)*ab(ib)
        DO ic=1,nni+1
           kcb=kb2(ic)/k2
           bm_ic=betab(ic)*mb(ic)
           a_ic =ab(ic)*ab(ic)
           IF (ic == ib) THEN
              ex=1.
           ELSE
              au=(a_ic-a_ib)*u
              ex=EXP(-au)
           ENDIF
           r_ib=r_ib + kcb*SQRT(bm_ic/bm_ib)*ex
        ENDDO      
        r_ib=1./r_ib
!*!
!r_ib=1.
!*!
!
        dab_reg=-r_ib*h/TWOPI
      END FUNCTION dab_reg

      SUBROUTINE a_reg_mass(nni, ia, ib, vp, k2, kb2, betab, mb, ac_r)
      USE physvars
        IMPLICIT NONE
        INTEGER,                     INTENT(IN)  :: nni 
        INTEGER,                     INTENT(IN)  :: ia
        INTEGER,                     INTENT(IN)  :: ib
        REAL,                        INTENT(IN)  :: vp
        REAL,                        INTENT(IN)  :: k2
        REAL,                        INTENT(IN)  :: kb2
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: betab
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: mb
        REAL,                        INTENT(OUT) :: ac_r
        REAL,    DIMENSION(1:nni+1)  :: ab
!       INTEGER, PARAMETER :: NR=10 ! integration regions: must be even
        INTEGER, PARAMETER :: NR=100 ! integration regions: must be even
        REAL,    PARAMETER :: UPM=0.7745966692E0 ! parameters for Gaussian Quad
        REAL,    PARAMETER :: W13=0.5555555556E0, W2=0.8888888889E0
        REAL               :: u0, u1, du, u, um, dab_reg
        INTEGER            :: iu
        ab=SQRT(0.5*betab*mb)*vp/CC
        ac_r=0
        u0=0.0
        u1=1.
!       u1=MIN(1.,5/(ab(ib)**2)) ! support can lie << 1
        du=(u1-u0)/NR
        u=u0-du
        DO iu=1,NR,2 ! Gaussian quadrature
           u=u+2.*du
           ac_r=ac_r+W2*dab_reg(u,vp,ia,ib,nni,k2,kb2,betab,mb)
           um=u-du*UPM
           ac_r=ac_r+W13*dab_reg(um,vp,ia,ib,nni,k2,kb2,betab,mb)
           um=u+du*UPM
           ac_r=ac_r+W13*dab_reg(um,vp,ia,ib,nni,k2,kb2,betab,mb)
        ENDDO
        ac_r=ac_r*du
!*!
!        um=1
!        ac_r=dab_reg(um,vp,ia,ib,nni,k2,kb2,betab,mb)
!*!
      END SUBROUTINE a_reg_mass
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! quantum contribution for non-zero electron mass
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
      FUNCTION daq(u, a, eta)
      USE physvars
        IMPLICIT NONE
        REAL,                        INTENT(IN)  :: u          ! [dimensionless]
        REAL,                        INTENT(IN)  :: a          ! [dimensionless]
        REAL,                        INTENT(IN)  :: eta        ! [dimensionless]
        REAL                                     :: daq  ! [dimensionless]
        REAL            :: repsi, au, eu, au2, ap, am, psilog, ch, sh, csh
        eu=eta/u 
        psilog=repsi(eu) - LOG(eu)
        au =2*a*u
        au2=a*u*u
        ap = au-au2-a
        am =-au-au2-a
        ch =0.5*(EXP(ap)+EXP(am))
        sh =0.5*(EXP(ap)-EXP(am))
        csh=2*(ch - sh/au)/au
        daq=-psilog*csh
      END FUNCTION daq

      SUBROUTINE a_quantum_mass(ia, ib, a, eta, aq)
        IMPLICIT NONE
        INTEGER, INTENT(IN)  :: ia    ! species index
        INTEGER, INTENT(IN)  :: ib    ! species index
        REAL,    INTENT(IN)  :: a     ! [dimensionless] (1/2) betab mb vp^2
        REAL,    INTENT(IN)  :: eta   ! [dimensionless] ep eb/4pi hbar vp
        REAL,    INTENT(OUT) :: aq 
        REAL               :: u0, u1, du, u, um
        INTEGER, PARAMETER :: NQ=1000            ! integration regions quantum : must be even
        REAL,    PARAMETER :: UPM=0.7745966692E0 ! parameters for Gaussian Quad
        REAL,    PARAMETER :: W13=0.5555555556E0, W2=0.8888888889E0
        REAL    :: daq
        INTEGER :: iu
        aq=0
        u0=0.
        aq=0
        IF (ib == ia) THEN
           u0=0
           u1=4./SQRT(a)
        ELSE
           u0=1-10./SQRT(a)
           u0=MAX(0.,u0)  
           u1=1+10./SQRT(a)
        ENDIF
        du=(u1-u0)/NQ
        u=u0-du
        DO iu=1,NQ,2 ! Gaussian quadrature
           u=u+2.E0*du
           aq=aq+W2*daq(u,a,eta)
           um=u-du*UPM
           aq=aq+W13*daq(um,a,eta)
           um=u+du*UPM
           aq=aq+W13*daq(um,a,eta)
        ENDDO
        aq=aq*du
      END SUBROUTINE a_quantum_mass

!
!=============== low- and high-energy asymptotic limits ====================
!

! 
! ROUTINE: SUBROUTINE coeff_bps_small_E(nni, ep, zp, mp, betab, zb, mb, nb, &
!   a_tot_lim, a_i_lim, a_e_lim, ac_tot_lim, ac_i_lim, ac_e_lim, aq_tot_lim,&
!   aq_i_lim, aq_e_lim)
!
! The asymptotic low energy regime E_p << T: this routine returns several 
! useful components of the corresponding A-coefficients in the low energy
! regime.
!
! UNITS: A_b has units of [MeV/micron] (subject to change in updates)
!
! The incident projectile and the background plasma. 
!
! projectile input quantities:
! ep : classical kinetic energy of the projectile [keV]
! zp : charge of the projectile in units of Z_p [dimensionless]
! mp : mass of the projectile [keV], i.e. mp = mp[grams]*c^2
!
! plasma input quantities:
! nni  : Number of total plasma species = number ion species + 1
! zb   : Charges of the plasma species. By convention zp(1) is the 
!      : electron plasma component. [dimensionless, Array]
! betab: Inverse temperatures of the plasma components. For an
!        electron-ion plasma, set betab(1)=1/T_e and all other
!        values of the array to 1/T_I.
! mb   : Masses of the plasma species [keV]. 
! nb   : Number densities of the plasma species [cm^-3]. 
!
! OUTPUT: a_tot_lim, a_i_lim, a_e_lim, ac_tot_lim, ac_i+lim, ac_e_lim, 
!         aq_tot_lim, aq_i_lim, aq_e_lim
!
! classical electron  : ac_e_lim
! classical ion       : ac_i_lim [sum over all ions]
! classical total     : ac_tot_lim = ac_e_lim + ac_i_lim 
! quantum   electron  : aq_e_lim
! quantum   ion       : aq_i_lim [sum over all ions]
! quantum   total     : a_tot_lim = aq_e_lim + aq_i_lim
! total     electric  : a_e_lim = ac_e_lim + aq_e_lim
! total     ion       : a_i_lim = ac_i_lim + aq_i_lim
! total               ! a_tot_lim = a_e_lim + a_i_lim
!
      SUBROUTINE coeff_bps_small_E(nni, ep, zp, mp, betab, zb, mb, nb,   &
            a_tot_lim, a_i_lim, a_e_lim, ac_tot_lim, ac_i_lim, ac_e_lim, &
            aq_tot_lim, aq_i_lim, aq_e_lim, ac_s_i_lim, ac_s_e_lim,      &
            ac_r_i_lim, ac_r_e_lim)
      USE physvars
      USE mathvars      
        IMPLICIT NONE
        INTEGER,                     INTENT(IN)  :: nni       !  number of ions
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: betab     !  temp array [1/keV]
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: mb        !  mass array [keV]
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: nb        !  density [1/cc]
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: zb        !  charge array
        REAL,                        INTENT(IN)  :: ep        !  projectile energy [keV]
        REAL,                        INTENT(IN)  :: mp        !  projectile mass   [keV]
        REAL,                        INTENT(IN)  :: zp        !  projectile charge
                                                              ! A-coeffs [MeV/micron]
        REAL,                        INTENT(OUT) :: a_tot_lim !  electron + ion
        REAL,                        INTENT(OUT) :: a_i_lim   !  ion contribution
        REAL,                        INTENT(OUT) :: a_e_lim   !  electron contribution
        REAL,                        INTENT(OUT) :: ac_tot_lim!  classical
        REAL,                        INTENT(OUT) :: ac_i_lim  !  classical
        REAL,                        INTENT(OUT) :: ac_e_lim  !  classical
        REAL,                        INTENT(OUT) :: aq_tot_lim!  quantum
        REAL,                        INTENT(OUT) :: aq_i_lim  !  quantum
        REAL,                        INTENT(OUT) :: aq_e_lim  !  quantum
        REAL,                        INTENT(OUT) :: ac_s_i_lim!  singular
        REAL,                        INTENT(OUT) :: ac_s_e_lim!  singular
        REAL,                        INTENT(OUT) :: ac_r_i_lim!  regular
        REAL,                        INTENT(OUT) :: ac_r_e_lim!  regular

        REAL    :: ac_r_lim, ac_s_lim, aq_lim
        INTEGER :: ib, nnb
        REAL,    DIMENSION(1:nni+1)  :: kb2  ! [1/cm^2]
        REAL,    DIMENSION(1:nni+1)  :: ab   ! [dimensionless]
        REAL    :: vp                        ! [cm/s]
        REAL    :: vp2                       ! [cm^2/s^2]
        REAL    :: kd2                       ! [1/cm^2]
        REAL    :: a, zp2                    ! [dimensionless]
        REAL    :: c1                        ! [keV/cm]
        REAL    :: c2                        ! [dimensionless]

        nnb =nni+1
        vp  = CC*SQRT(2*ep/mp)      ! [cm/s]
        vp2 = vp*vp                 ! [cm^2/s^2]
!       kb2 = DEBYE2*zb*zb*nb*betab ! [1/cm^2]
        kb2 =8*PI*A0CM*BEKEV*zb*zb*nb*betab
        kd2 = SUM(kb2)              ! [1/cm^2]
        zp2 = zp**2                 ! [dimensionless]
        ab  =0.5*betab*mb*vp2/CC2   ! [dimensionless]
                                    ! ab=(1/2) betab(ib)*mbc2(ib)*vp2/CC2 
!
! initialize A-coefficients
!
        a_tot_lim =0  ! electron + ion
        a_i_lim   =0  ! ion contribution
        a_e_lim   =0  ! electron contribution
        ac_tot_lim=0  ! classical total
        ac_e_lim  =0  ! classical electron
        ac_i_lim  =0  ! classical ion
        aq_tot_lim=0  ! quantum total
        aq_e_lim  =0  ! quantum electron
        aq_i_lim  =0  ! quantum ion
        ac_s_i_lim=0 
        ac_s_e_lim=0 
        ac_r_i_lim=0
        ac_r_e_lim=0
!
        DO ib=1,nni+1
        IF (zb(ib) .NE. 0.) THEN

           a=ab(ib)                    ! [dimensionless] 
           c1=2*zp2*BEKEV*kb2(ib)*A0CM ! [keV/cm]
           c1=c1*1.E-7                 ! [MeV/micron]
           c2=SQRT(a/PI)               ! [dimensionless] 
!
! singular: asymptotic low energy form
!
           CALL a_sing_ib_small_E(nni,ib,ep,zp,mp,betab,zb,mb,nb,ac_s_lim)
           ac_s_lim=c1*c2*ac_s_lim
!
! regular: asymptotic low energy form
!
           CALL a_reg_ib_small_E(nni,ib,ep,zp,mp,betab,zb,mb,nb,ac_r_lim)
           ac_r_lim=c1*c2*ac_r_lim
!
! quantum: asymptotic low energy form
!
          CALL aq_ib_small_E(nni, ib, ep, zp, mp, betab, zb, mb, nb, aq_lim)
          aq_lim=c1*c2*aq_lim
!
! collect components
!
           CALL x_collect(ib,nnb,ac_s_lim,ac_r_lim,aq_lim,a_tot_lim,a_i_lim,   &
             a_e_lim,ac_tot_lim,ac_i_lim,ac_e_lim,aq_tot_lim,aq_i_lim,aq_e_lim,&
             ac_s_i_lim, ac_s_e_lim, ac_r_i_lim, ac_r_e_lim)
        ENDIF
        ENDDO
      END SUBROUTINE coeff_bps_small_E


      SUBROUTINE a_sing_ib_small_E(nni, ib, ep, zp, mp, betab, zb, mb, nb, ac_s_lim)
      USE physvars
      USE mathvars      
        IMPLICIT NONE
        INTEGER,                     INTENT(IN)  :: nni      !  number of ions
        INTEGER,                     INTENT(IN)  :: ib       !  plasma species number
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: betab    !  temp array [1/keV]
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: mb       !  mass array [keV]
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: nb       !  density [1/cc]
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: zb       !  charge array
        REAL,                        INTENT(IN)  :: ep       !  projectile energy [keV]
        REAL,                        INTENT(IN)  :: mp       !  projectile mass   [keV]
        REAL,                        INTENT(IN)  :: zp       !  projectile charge
        REAL,                        INTENT(OUT) :: ac_s_lim !  A-coeffs [MeV/micron]

        REAL,    DIMENSION(1:nni+1)  :: kb2  ! [1/cm^2]
        REAL,    DIMENSION(1:nni+1)  :: ab   ! [dimensionless]
        REAL,    DIMENSION(1:nni+1)  :: mpb  ! [keV]
        REAL,    DIMENSION(1:nni+1)  :: mbpb ! [dimensionless]
        REAL    :: vp                        ! [cm/s]
        REAL    :: vp2                       ! [cm^2/s^2]
        REAL    :: kd2, k2                   ! [1/cm^2]
        REAL    :: kd, k                     ! [1/cm]
        REAL    :: a, zp2                    ! [dimensionless]

        vp  = CC*SQRT(2*ep/mp)      ! [cm/s]
        vp2 = vp*vp                 ! [cm^2/s^2]
        kb2 =8*PI*A0CM*BEKEV*zb*zb*nb*betab
        kd2 = SUM(kb2)              ! [1/cm^2]
        kd  = SQRT(kd2)             ! [1/cm]
        k2  = kb2(1)                ! [1/cm^2]  k2=k_e^2
        k   = SQRT(k2)              ! [1/cm]    k =k_e
        zp2 = zp**2                 ! [dimensionless]
                                    ! ab=(1/2) betab(ib)*mbc2(ib)*vp2/CC2
        ab  =0.5*betab*mb*vp2/CC2   ! [dimensionless] 
        mpb = mp*mb/(mp+mb)         ! [keV]
        mbpb= mb/mpb                ! [dimensionless]

        IF (SMALL_E_NLO) THEN
          a=ab(ib)             
        ELSE
          a=0
        ENDIF
!
! singular: asymptotic low energy form
!
        ac_s_lim=-(2./3. - 2*a/5.)*(LOG(betab(ib)*BEKEV*ABS(zp*zb(ib))* &
        0.5*k*A0CM*mbpb(ib) ) + 2*GAMMA) + (4./15.)*a
      END SUBROUTINE a_sing_ib_small_E


      SUBROUTINE a_reg_ib_small_E_old(nni, ib, ep, zp, mp, betab, zb, mb, nb, ac_r_lim)
      USE physvars
      USE mathvars      
        IMPLICIT NONE
        INTEGER,                     INTENT(IN)  :: nni      !  number of ions
        INTEGER,                     INTENT(IN)  :: ib       !  plasma species number
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: betab    !  temp array [1/keV]
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: mb       !  mass array [keV]
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: nb       !  density [1/cc]
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: zb       !  charge array
        REAL,                        INTENT(IN)  :: ep       !  projectile energy [keV]
        REAL,                        INTENT(IN)  :: mp       !  projectile mass   [keV]
        REAL,                        INTENT(IN)  :: zp       !  projectile charge
        REAL,                        INTENT(OUT) :: ac_r_lim !  A-coeffs [MeV/micron]

        REAL,    DIMENSION(1:nni+1)  :: kb2  ! [1/cm^2]
        REAL,    DIMENSION(1:nni+1)  :: ab   ! [dimensionless]
        REAL,    DIMENSION(1:nni+1)  :: ab2  ! [dimensionless]
        REAL,    DIMENSION(1:nni+1)  :: mpb  ! [keV]
        REAL,    DIMENSION(1:nni+1)  :: mbpb ! [dimensionless]
        REAL    :: vp                        ! [cm/s]
        REAL    :: vp2                       ! [cm^2/s^2]
        REAL    :: kd2, k2                   ! [1/cm^2]
        REAL    :: kd, k                     ! [1/cm]
        REAL    :: a, zp2                    ! [dimensionless]
        REAL    :: ar1, ar2                  ! [dimensionless]

        vp  = CC*SQRT(2*ep/mp)      ! [cm/s]
        vp2 = vp*vp                 ! [cm^2/s^2]
        kb2 =8*PI*A0CM*BEKEV*zb*zb*nb*betab
        kd2 = SUM(kb2)              ! [1/cm^2]
        kd  = SQRT(kd2)             ! [1/cm]
        k2  = kb2(1)                ! [1/cm^2]  k2=k_e^2
        k   = SQRT(k2)              ! [1/cm]    k =k_e
        zp2 = zp**2                 ! [dimensionless]
                                    ! ab=(1/2) betab(ib)*mbc2(ib)*vp2/CC2
        ab  =0.5*betab*mb*vp2/CC2   ! [dimensionless] 
        mpb = mp*mb/(mp+mb)         ! [keV]
        mbpb= mb/mpb                ! [dimensionless]
        a=ab(ib)                    ! [dimensionless] 
!
! regular: asymptotic low energy form
!
        ab2 = SQRT(ab)              ! [dimensionless]
        ar1 =-2*SUM(kb2*ab)/k2/5.   ! [dimensionless] coeff for A_reg with E<<T
        ar2 = SUM(kb2*ab2)/k2       ! [dimensionless] coeff for A_reg with E<<T
        ar2 = ar2*ar2*PI/30.        !
        ac_r_lim=-((THIRD - 0.2*a)*(LOG(kd2/k2)+1) + ar1 + ar2 )
      END SUBROUTINE a_reg_ib_small_E_old



      SUBROUTINE a_reg_ib_small_E(nni, ib, ep, zp, mp, betab, zb, mb, nb, ac_r_lim)
      USE physvars
      USE mathvars      
        IMPLICIT NONE
        INTEGER,                     INTENT(IN)  :: nni      !  number of ions
        INTEGER,                     INTENT(IN)  :: ib       !  plasma species number
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: betab    !  temp array [1/keV]
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: mb       !  mass array [keV]
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: nb       !  density [1/cc]
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: zb       !  charge array
        REAL,                        INTENT(IN)  :: ep       !  projectile energy [keV]
        REAL,                        INTENT(IN)  :: mp       !  projectile mass   [keV]
        REAL,                        INTENT(IN)  :: zp       !  projectile charge
        REAL,                        INTENT(OUT) :: ac_r_lim !  A-coeffs [MeV/micron]

        REAL,    DIMENSION(1:nni+1)  :: kb2  ! [1/cm^2]
        REAL,    DIMENSION(1:nni+1)  :: ab   ! [dimensionless]
        REAL,    DIMENSION(1:nni+1)  :: mpb  ! [keV]
        REAL,    DIMENSION(1:nni+1)  :: mbpb ! [dimensionless]
        REAL    :: vp                        ! [cm/s]
        REAL    :: vp2                       ! [cm^2/s^2]
        REAL    :: kd2, k2                   ! [1/cm^2]
        REAL    :: kd, k                     ! [1/cm]
        REAL    :: a, zp2                    ! [dimensionless]
        REAL    :: ar1, ar2, logone ! [dimensionless]

        vp  = CC*SQRT(2*ep/mp)      ! [cm/s]
        vp2 = vp*vp                 ! [cm^2/s^2]
        kb2 =8*PI*A0CM*BEKEV*zb*zb*nb*betab
        kd2 = SUM(kb2)              ! [1/cm^2]
        kd  = SQRT(kd2)             ! [1/cm]
        k2  = kb2(1)                ! [1/cm^2]  k2=k_e^2
        k   = SQRT(k2)              ! [1/cm]    k =k_e
        zp2 = zp**2                 ! [dimensionless]
        logone=LOG(kd2/k2)+1        ! [dimensionless]
                                    ! ab=(1/2) betab(ib)*mbc2(ib)*vp2/CC2
        ab  =0.5*betab*mb*vp2/CC2   ! [dimensionless] 
        mpb = mp*mb/(mp+mb)         ! [keV]
        mbpb= mb/mpb                ! [dimensionless]

        IF (SMALL_E_NLO) THEN       ! NLO in v_p correction
          a=ab(ib)                  ! [dimensionless] 
          ar1=2*SUM(kb2*ab)/kd2/5.
          ar2=SUM(kb2*SQRT(2*ab))/kd2
          ar2=-PI*ar2*ar2/60.
        ELSE
          a=0    ! set NLO to zero
          ar1=0  !
          ar2=0  !
        ENDIF
        ac_r_lim=(-THIRD + a/10.)*(1 + LOG(kd2/k2)) + ar1 + ar2 ! LO + NLO
      END SUBROUTINE a_reg_ib_small_E

      SUBROUTINE aq_ib_small_E(nni, ib, ep, zp, mp, betab, zb, mb, nb, aq_lim)
      USE physvars
      USE mathvars      
        IMPLICIT NONE
        INTEGER,                     INTENT(IN)  :: nni      !  number of ions
        INTEGER,                     INTENT(IN)  :: ib       !  plasma species number
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: betab    !  temp array [1/keV]
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: mb       !  mass array [keV]
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: nb       !  density [1/cc]
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: zb       !  charge array
        REAL,                        INTENT(IN)  :: ep       !  projectile energy [keV]
        REAL,                        INTENT(IN)  :: mp       !  projectile mass   [keV]
        REAL,                        INTENT(IN)  :: zp       !  projectile charge
        REAL,                        INTENT(OUT) :: aq_lim   !  A-coeffs [MeV/micron]

        REAL,    DIMENSION(1:nni+1)  :: kb2  ! [1/cm^2]
        REAL,    DIMENSION(1:nni+1)  :: ab   ! [dimensionless]
        REAL,    DIMENSION(1:nni+1)  :: mpb  ! [keV]
        REAL,    DIMENSION(1:nni+1)  :: mbpb ! [dimensionless]
        REAL    :: vp                        ! [cm/s]
        REAL    :: vp2                       ! [cm^2/s^2]
        REAL    :: kd2, k2                   ! [1/cm^2]
        REAL    :: kd, k                     ! [1/cm]
        REAL    :: a, zp2                    ! [dimensionless]
        REAL    :: etbar, etbar2             ! [dimensionless]

        vp  = CC*SQRT(2*ep/mp)      ! [cm/s]
        vp2 = vp*vp                 ! [cm^2/s^2]
        kb2 =8*PI*A0CM*BEKEV*zb*zb*nb*betab
        kd2 = SUM(kb2)              ! [1/cm^2]
        kd  = SQRT(kd2)             ! [1/cm]
        k2  = kb2(1)                ! [1/cm^2]  k2=k_e^2
        k   = SQRT(k2)              ! [1/cm]    k =k_e
        zp2 = zp**2                 ! [dimensionless]
                                    ! ab=(1/2) betab(ib)*mbc2(ib)*vp2/CC2
        ab  =0.5*betab*mb*vp2/CC2   ! [dimensionless] 
        mpb = mp*mb/(mp+mb)         ! [keV]
        mbpb= mb/mpb                ! [dimensionless]
!
        a=ab(ib)                    ! [dimensionless] 
!
! quantum: asymptotic low energy form: etbar defined with thermal velocity
! [$\bar eta_{pb} = e_p e_b/4\pi \hbar \bar v_b $  with $\bar v_b^2 = 3 T_b/m_b$]
!
           etbar=4.2115E-3*ABS(zp*zb(ib))*SQRT(betab(ib)*mb(ib))
           etbar2=etbar*etbar
           IF (ib==1) THEN
              aq_lim=LOG(1.5*etbar2)/3. + GAMMA  ! electrons only eta_pe << 1
           ELSE
              aq_lim=-1./(27*etbar2)
           ENDIF
      END SUBROUTINE aq_ib_small_E

! 
! ROUTINE: SUBROUTINE coeff_bps_high_E(nni, ep, zp, mp, betab, zb, mb, nb, &
!   a_tot_lim, a_i_lim, a_e_lim, ac_tot_lim, ac_i_lim, ac_e_lim, aq_tot_lim,&
!   aq_i_lim, aq_e_lim)
!
! Returns high energy asymptotic regimes. For the ions this means E_p >> T.
! For electrons there are two regimes: 
! (i)  extreme high energy E_p >> (m_I/m_e)*T [coeff_bps_very_high_E]
! (ii) intermediate high energy T << E_p << (m_I/m_e)*T [this routine]
!
! UNITS: A_b has units of [MeV/micron] (subject to change in updates)
!
! The incident projectile and the background plasma. 
!
! projectile input quantities:
! ep : classical kinetic energy of the projectile [keV]
! zp : charge of the projectile in units of Z_p [dimensionless]
! mp : mass of the projectile [keV], i.e. mp = mp[grams]*c^2
!
! plasma input quantities:
! nni  : Number of total plasma species = number ion species + 1
! zb   : Charges of the plasma species. By convention zp(1) is the 
!      : electron plasma component. [dimensionless, Array]
! betab: Inverse temperatures of the plasma components. For an
!        electron-ion plasma, set betab(1)=1/T_e and all other
!        values of the array to 1/T_I.
! mb   : Masses of the plasma species [keV]. 
! nb   : Number densities of the plasma species [cm^-3]. 
!
! OUTPUT: a_tot_lim, a_i_lim, a_e_lim, ac_tot_lim, ac_i+lim, ac_e_lim, 
!         aq_tot_lim, aq_i_lim, aq_e_lim
!
! classical electron  : ac_e_lim
! classical ion       : ac_i_lim [sum over all ions]
! classical total     : ac_tot_lim = ac_e_lim + ac_i_lim 
! quantum   electron  : aq_e_lim
! quantum   ion       : aq_i_lim [sum over all ions]
! quantum   total     : a_tot_lim = aq_e_lim + aq_i_lim
! total     electric  : a_e_lim = ac_e_lim + aq_e_lim
! total     ion       : a_i_lim = ac_i_lim + aq_i_lim
! total               ! a_tot_lim = a_e_lim + a_i_lim
!
      SUBROUTINE coeff_bps_high_E(nni, ep, zp, mp, betab, zb, mb, nb,    &
            a_tot_lim, a_i_lim, a_e_lim, ac_tot_lim, ac_i_lim, ac_e_lim, &
            aq_tot_lim, aq_i_lim, aq_e_lim, ac_s_i_lim, ac_s_e_lim,      &
            ac_r_i_lim, ac_r_e_lim)
      USE physvars
      USE mathvars      
        IMPLICIT NONE
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
                                                               ! A-coeffs [MeV/micron]
        REAL,                        INTENT(OUT) :: a_tot_lim  !  electron + ion
        REAL,                        INTENT(OUT) :: a_i_lim    !  ion contribution
        REAL,                        INTENT(OUT) :: a_e_lim    !  electron contribution
        REAL,                        INTENT(OUT) :: ac_tot_lim !  classical
        REAL,                        INTENT(OUT) :: ac_i_lim   !  classical
        REAL,                        INTENT(OUT) :: ac_e_lim   !  classical
        REAL,                        INTENT(OUT) :: aq_tot_lim !  quantum
        REAL,                        INTENT(OUT) :: aq_i_lim   !  quantum
        REAL,                        INTENT(OUT) :: aq_e_lim   !  quantum
        REAL,                        INTENT(OUT) :: ac_s_i_lim
        REAL,                        INTENT(OUT) :: ac_s_e_lim
        REAL,                        INTENT(OUT) :: ac_r_i_lim
        REAL,                        INTENT(OUT) :: ac_r_e_lim

        REAL    :: ac_r_lim, ac_s_lim, aq_lim
        INTEGER :: ib, nnb
        REAL    :: vp                        ! [cm/s]
        REAL    :: vp2                       ! [cm^2/s^2]
        REAL    :: a, zp2                    ! [dimensionless]
        REAL    :: c1, cs                    ! [keV/cm]
        REAL    :: c2                        ! [dimensionless]
        REAL    :: k                         ! [1/cm]
        REAL    :: k2, ke2                   ! [1/cm^2]
        REAL    :: eta                       ! [dimensionless]
        REAL    :: te, mec2, mpec22          ! [keV] 
        REAL,    DIMENSION(1:nni+1)  :: ab   ! [dimensionless]
        REAL,    DIMENSION(1:nni+1)  :: mpb  ! [keV]
        REAL,    DIMENSION(1:nni+1)  :: kb2  ! [1/cm^2
!
! omi2, and omi needed only for 
! E >> (mI/me)*T: extreme high energy limit for electrons
!
!       REAL    :: omi2, omi                 ! [1/s^2, 1/s]]

        nnb = nni+1
        vp  = CC*SQRT(2*ep/mp)      ! [cm/s]
        vp2 = vp*vp                 ! [cm^2/s^2]
        kb2 =8*PI*A0CM*BEKEV*zb*zb*nb*betab
        k2  = kb2(1)                ! [1/cm^2]  k2=k_e^2
        k   = SQRT(k2)              ! [1/cm]    k =k_e
        zp2 = zp**2                 ! [dimensionless]
                                    ! ab=(1/2) betab(ib)*mb(ib)*vp2/CC2
        ab  =0.5*betab*mb*vp2/CC2   ! [dimensionless] 
        mpb = mp*mb/(mp+mb)         ! [keV]
!
! initialize A-coefficients
!
        a_tot_lim = 0  ! electron + ion
        a_i_lim   = 0  ! ion contribution
        a_e_lim   = 0  ! electron contribution
        ac_tot_lim= 0  ! classical total
        ac_e_lim  = 0  ! classical electron
        ac_i_lim  = 0  ! classical ion
        aq_tot_lim= 0  ! quantum total
        aq_e_lim  = 0  ! quantum electron
        aq_i_lim  = 0  ! quantum ion
        ac_s_i_lim=0
        ac_s_e_lim=0 
        ac_r_i_lim=0 
        ac_r_e_lim=0
!
        DO ib=1,nni+1
        IF (zb(ib) .NE. 0.) THEN

           a=ab(ib)                    ! [dimensionless] 
           c1=2*zp2*BEKEV*kb2(ib)*A0CM ! [keV/cm]
           c1=c1*1.E-7                 ! [MeV/micron]
           c2=SQRT(a/PI)               ! [dimensionless] 
!
! singular: asymptotic high energy form [need electrons]
!
! At this point I only have an asymptotic form for the total 
! electron contribution (classical + quantum = sing + reg + quantum).
! I'll write this expression to ac_s_lim for now. This will give
! aq_lim=0 and ac_lim=a_e.  FIX LATER. 
!
           IF (ib == 1) THEN              
!
! T << E << (mI/me)*T: intermediate high energy limit for electrons
!
              te=1./betab(ib)
              mpec22=mpb(ib)**2
              mec2  =mb(ib)
              ke2   =kb2(ib)
              ac_s_lim=THIRD*(LOG(8*te*mpec22/(mec2*HBARC**2*ke2)) - GAMMA -1)
              ac_s_lim=c1*c2*ac_s_lim              ! [MeV/micron]
           ELSE
              cs=0.5*c1/ab(ib)
              ac_s_lim=(LOG(ABS(zp*zb(ib))*BEKEV*& ! [dimensionless]
                k*A0CM*CC2/(mpb(ib)*vp2)) + GAMMA) !
              ac_s_lim=-cs*ac_s_lim                ! [MeV/micron]
           ENDIF
!
! regular: asymptotic high energy form [need electrons]
!
           IF (ib == 1) THEN
              ac_r_lim=0
           ELSE
               ac_r_lim=-0.25*c1/ab(ib)                       ! [MeV/micron]
           ENDIF
!
! quantum: asymptotic high energy form [need electrons]
!
           IF (ib == 1) THEN
              aq_lim=0 
           ELSE
              eta =ABS(zp*zb(ib))*2.1870E8/vp  ! [dimensionless] quantum parameter
              aq_lim=LOG(eta) + GAMMA
              aq_lim=aq_lim*c1/a/2
           ENDIF
           CALL x_collect(ib,nnb,ac_s_lim,ac_r_lim,aq_lim,a_tot_lim,a_i_lim,    &
             a_e_lim,ac_tot_lim,ac_i_lim,ac_e_lim,aq_tot_lim,aq_i_lim,aq_e_lim, &
             ac_s_i_lim, ac_s_e_lim, ac_r_i_lim, ac_r_e_lim)
        ENDIF
        ENDDO
!*!
!a_i_lim=c1*c2*LOG(SUM(kb2)/k2)/3     ! ***===*** 
!*!
      END SUBROUTINE coeff_bps_high_E
!
!======================================

! Extreme high energy limit for electrons: E >> (mI/me)*T
! The ions have the same limit in coeff_bps_high_E and coeff_bps_very_high_E.
!
! This limit is not really self-consistent since relativistic effects become
! important.
!
!
!======================================
      SUBROUTINE coeff_bps_very_high_E(nni, ep, zp, mp, betab, zb, mb, nb, &
            a_tot_lim, a_i_lim, a_e_lim, ac_tot_lim, ac_i_lim, ac_e_lim,   &
            aq_tot_lim, aq_i_lim, aq_e_lim, ac_s_i_lim, ac_s_e_lim,        &
            ac_r_i_lim, ac_r_e_lim)
      USE physvars
      USE mathvars      
        IMPLICIT NONE
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
                                                               ! A-coeffs [MeV/micron]
        REAL,                        INTENT(OUT) :: a_tot_lim  !  electron + ion
        REAL,                        INTENT(OUT) :: a_i_lim    !  ion contribution
        REAL,                        INTENT(OUT) :: a_e_lim    !  electron contribution
        REAL,                        INTENT(OUT) :: ac_tot_lim !  classical
        REAL,                        INTENT(OUT) :: ac_i_lim   !  classical
        REAL,                        INTENT(OUT) :: ac_e_lim   !  classical
        REAL,                        INTENT(OUT) :: aq_tot_lim !  quantum
        REAL,                        INTENT(OUT) :: aq_i_lim   !  quantum
        REAL,                        INTENT(OUT) :: aq_e_lim   !  quantum
        REAL,                        INTENT(OUT) :: ac_s_i_lim
        REAL,                        INTENT(OUT) :: ac_s_e_lim
        REAL,                        INTENT(OUT) :: ac_r_i_lim
        REAL,                        INTENT(OUT) :: ac_r_e_lim

        REAL    :: omi2, omi
        REAL    :: ac_r_lim, ac_s_lim, aq_lim
        INTEGER :: ib, nnb
        REAL    :: vp                        ! [cm/s]
        REAL    :: vp2                       ! [cm^2/s^2]
        REAL    :: a, zp2                    ! [dimensionless]
        REAL    :: c1, cs                    ! [keV/cm]
        REAL    :: c2                        ! [dimensionless]
        REAL    :: k                         ! [1/cm]
        REAL    :: k2                        ! [1/cm^2]
        REAL    :: eta                       ! [dimensionless]
        REAL,    DIMENSION(1:nni+1)  :: ab   ! [dimensionless]
        REAL,    DIMENSION(1:nni+1)  :: mpb  ! [keV]
        REAL,    DIMENSION(1:nni+1)  :: kb2  ! [1/cm^2
!
! omi2, and omi needed only for 
! E >> (mI/me)*T: extreme high energy limit for electrons
!
!       REAL    :: omi2, omi                 ! [1/s^2, 1/s]]

        nnb = nni+1
        vp  = CC*SQRT(2*ep/mp)      ! [cm/s]
        vp2 = vp*vp                 ! [cm^2/s^2]
        kb2 =8*PI*A0CM*BEKEV*zb*zb*nb*betab
        k2  = kb2(1)                ! [1/cm^2]  k2=k_e^2
        k   = SQRT(k2)              ! [1/cm]    k =k_e
        zp2 = zp**2                 ! [dimensionless]
                                    ! ab=(1/2) betab(ib)*mb(ib)*vp2/CC2
        ab  =0.5*betab*mb*vp2/CC2   ! [dimensionless] 
        mpb = mp*mb/(mp+mb)         ! [keV]
!
! initialize A-coefficients
!
        a_tot_lim = 0  ! electron + ion
        a_i_lim   = 0  ! ion contribution
        a_e_lim   = 0  ! electron contribution
        ac_tot_lim= 0  ! classical total
        ac_e_lim  = 0  ! classical electron
        ac_i_lim  = 0  ! classical ion
        aq_tot_lim= 0  ! quantum total
        aq_e_lim  = 0  ! quantum electron
        aq_i_lim  = 0  ! quantum ion
        ac_s_i_lim=0 
        ac_s_e_lim=0 
        ac_r_i_lim=0 
        ac_r_e_lim=0
!
        DO ib=1,nni+1
        IF (zb(ib) .NE. 0.) THEN

           a=ab(ib)                    ! [dimensionless] 
           c1=2*zp2*BEKEV*kb2(ib)*A0CM ! [keV/cm]
           c1=c1*1.E-7                 ! [MeV/micron]
           c2=SQRT(a/PI)               ! [dimensionless] 
!
! singular: asymptotic high energy form [need electrons]
!
! At this point I only have an asymptotic form for the total 
! electron contribution (classical + quantum = sing + reg + quantum).
! I'll write this expression to ac_s_lim for now. This will give
! aq_lim=0 and ac_lim=a_e.  FIX LATER. 
!
           IF (ib == 1) THEN              
!
! E >> (mI/me)*T: extreme high energy limit for electrons
!
              omi2=OMEGI2*zb(ib)*zb(ib)*nb(ib)*AMUKEV/mb(ib) ! [1/s^2]
              omi =SQRT(omi2)
              ac_s_lim=0.5*c1*LOG(2*mpb(ib)*vp2/(HBARC*CC*omi))/ab(ib)
           ELSE
              cs=0.5*c1/ab(ib)
              ac_s_lim=(LOG(ABS(zp*zb(ib))*BEKEV*& ! [dimensionless]
                k*A0CM*CC2/(mpb(ib)*vp2)) + GAMMA) !
              ac_s_lim=-cs*ac_s_lim                ! [MeV/micron]
           ENDIF
!
! regular: asymptotic high energy form [need electrons]
!
           IF (ib == 1) THEN
              ac_r_lim=0
           ELSE
               ac_r_lim=-0.25*c1/ab(ib)                       ! [MeV/micron]
           ENDIF
!
! quantum: asymptotic high energy form [need electrons]
!
           IF (ib == 1) THEN
              aq_lim=0 
           ELSE
              eta =ABS(zp*zb(ib))*2.1870E8/vp  ! [dimensionless] quantum parameter
              aq_lim=LOG(eta) + GAMMA
              aq_lim=aq_lim*c1/a/2
           ENDIF
           CALL x_collect(ib,nnb,ac_s_lim,ac_r_lim,aq_lim,a_tot_lim,a_i_lim,   &
             a_e_lim,ac_tot_lim,ac_i_lim,ac_e_lim,aq_tot_lim,aq_i_lim,aq_e_lim,&
             ac_s_i_lim, ac_s_e_lim, ac_r_i_lim, ac_r_e_lim)
        ENDIF
        ENDDO
      END SUBROUTINE coeff_bps_very_high_E


