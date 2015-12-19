!
! Stopping poweer models from the literature.
!
!

! Fraley et al.
! Phys of Fluids 17 (1974) 474
!
! Stopping power for an alpha in DT is given on p. 475. It depends only upon
! the alpha particle density and the electron temperature.
!
      SUBROUTINE dedx_fraley(nni, ep, mb, nb, betab, dedx_tot, dedx_i, dedx_e)
      USE physvars
      USE mathvars      
        IMPLICIT NONE                                         ! Plasma:
        INTEGER,                     INTENT(IN)  :: nni       !  number of ions
        REAL,                        INTENT(IN)  :: ep        !  projectile energy [keV]
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: mb        !  mass array [keV]
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: nb        !  number density array [cm^-3]
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: betab     !  temp array [1/keV]


                                                              !
                                                              ! dE/dx [MeV/micron]
        REAL,                        INTENT(OUT) :: dedx_tot  !  electron + ion
        REAL,                        INTENT(OUT) :: dedx_i    !  ion contribution
        REAL,                        INTENT(OUT) :: dedx_e    !  electron contribution

        REAL,    PARAMETER           :: units=1.E-7        ! 10^-7 => MeV/micron; 10^-4 keV/micron
        REAL,    PARAMETER           :: MKEV2G=1.78273E-30 ! convert mass in keV/c^2 to g
        REAL,    PARAMETER           :: E0=3500.           ! threshold alpha energy in keV
        REAL,    PARAMETER           :: RHO0=0.213         ! g/cm^3 density of solid DT
        REAL     :: te, u, rho
    
!       REAL,    DIMENSION(1:nni+1)  :: kkb2, aab, ab2, c1b, c2b, c3b,c4b, eetb
!
! initialize components of dE/dx
!

       te =1./betab(1) ! electron temp in keV
       u  =ep/E0       ! normalized alpha particle energy
       rho=SUM(mb(1:nni+1)*nb(1:nni+1))*MKEV2G  ! total mass density of DT plasma in g/cm^3

       dedx_e=23.2*(rho/RHO0)*(u**0.5/te**1.5)*(1 + 0.17*LOG(te*SQRT(RHO0/rho)))*E0*units
       dedx_i=0.047*(rho/RHO0)*(1/u)*(1 + 0.075*LOG(SQRT((te*RHO0)/(rho*u))))*E0*units

       dedx_tot = dedx_e + dedx_i

      END SUBROUTINE dedx_fraley

! simple model for the transition to strong coupling. Note: the coefficients d_i 
! and d_e should have an unknown energy dependence in the strong regime. 
!
!
      SUBROUTINE dedx_bps_strong1(nni, ep, zp, mp, betab, zb, mb, nb, d_i, d_e,  &
            dedx_tot, dedx_i, dedx_e, dedxc_tot, dedxc_i, dedxc_e, & 
            dedxq_tot, dedxq_i, dedxq_e)
      USE physvars
      USE mathvars      
        IMPLICIT NONE                                      ! Plasma:
        INTEGER,                     INTENT(IN)  :: nni    !  number of ions
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: betab  !  temp array [1/keV]
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: mb     !  mass array [keV]
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: nb     !  density [1/cc]
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: zb     !  charge array
        REAL,                        INTENT(IN)  :: d_i    !  strong coupling ion parameter
        REAL,                        INTENT(IN)  :: d_e    !  strong coupling electron parameter
                                                           !
                                                           ! Projectile  
        REAL,                        INTENT(IN)  :: ep     !  projectile energy [keV]
        REAL,                        INTENT(IN)  :: mp     !  projectile mass   [keV]
        REAL,                        INTENT(IN)  :: zp     !  projectile charge
                                                           !
                                                              ! dE/dx [MeV/micron]
        REAL,                        INTENT(OUT) :: dedx_tot  !  electron + ion
        REAL,                        INTENT(OUT) :: dedx_i    !  ion contribution
        REAL,                        INTENT(OUT) :: dedx_e    !  electron contribution
        REAL,                        INTENT(OUT) :: dedxc_tot !  classical
        REAL,                        INTENT(OUT) :: dedxc_i   !  classical
        REAL,                        INTENT(OUT) :: dedxc_e   !  classical
        REAL,                        INTENT(OUT) :: dedxq_tot !  quantum
        REAL,                        INTENT(OUT) :: dedxq_i   !  quantum
        REAL,                        INTENT(OUT) :: dedxq_e   !  quantum

        REAL                        :: st_i, st_e
        REAL,    DIMENSION(1:nni+1) :: gb   !  Plasma coupling array
        REAL                        :: ge   !  Electron coupling
        REAL                        :: gi   !  Average ion coupling
        REAL                        :: etae !  Quantum parameter for electron
        REAL                        :: ze   !  Fugacity


! Call weakly coupled BPS stopping power
!
        CALL dedx_bps(nni, ep, zp, mp, betab, zb, mb, nb,   &
             dedx_tot, dedx_i, dedx_e, dedxc_tot, dedxc_i,  & 
             dedxc_e, dedxq_tot, dedxq_i, dedxq_e) ! [MeV/micron]

! Calculate plasma coupling constant g
!
        CALL param(nni, betab, nb, gb, ge, gi, etae, ze)
        st_i = 1. + ge*d_i
        st_e = 1. + ge*d_e

        dedx_i    = dedx_i*st_i
        dedx_e    = dedx_e*st_e
        dedx_tot  = dedx_i + dedx_e

        dedxc_i   = dedxc_i*st_i
        dedxc_e   = dedxc_e*st_e
        dedxc_tot = dedxc_i + dedxc_e

        dedxq_i   = dedxq_i*st_i
        dedxq_e   = dedxq_e*st_e
        dedxq_tot = dedxq_i + dedxq_e

      END SUBROUTINE dedx_bps_strong1
