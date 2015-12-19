!
! Not finished...
!
! This program calculates the rate C_{eI} in three ways:
!  1. A-coeff folded into distribution function
!  2. Sum rule or Born approximation
!  3. Arbitrary mass ratio
!
      PROGRAM main
      USE allocatablevars
      USE bpsvars
      USE physvars
      USE mathvars
        IMPLICIT NONE
        REAL     :: te, ti, ne
        INTEGER  :: nni

        REAL     :: ln_bps_mass, cei_mass, cei_reg_mass
        REAL     :: cei_sing_mass, cei_qm_mass, delta_mass

        REAL     :: cei_born, ln_bps_born, ln_bps, delta

        REAL     :: cei_dist_ei, cei_sing_dist_ei, cei_reg_dist_ei, cei_qm_dist_ei
        REAL     :: cei_dist_ei1, cei_dist_ei2, cei_dist_ei3, cei_dist_ei4

        REAL     :: cei_tot, cei_i, cei_e, ceic_tot, ceic_i, ceic_e, ceiq_tot
        REAL     :: ceiq_i, ceiq_e, ceic_s_i, ceic_s_e, ceic_r_i , ceic_r_e


        REAL     :: x0, x1, x2, x3, x4, dx, x
        INTEGER  :: nmax, ix, ib

        REAL     :: ze, yei, norme, me, betae, em
        REAL     :: a_tot, a_i, a_e, ac_tot, ac_i, ac_e, aq_tot, aq_i, aq_e
        REAL     :: ac_s_i, ac_s_e, ac_r_i, ac_r_e

        REAL,    DIMENSION(:), ALLOCATABLE :: ceib
!*!
!       REAL     ::  a_tot_lim, a_i_lim, a_e_lim, ac_tot_lim, ac_i_lim, ac_e_lim
!       REAL     ::  aq_tot_lim, aq_i_lim, aq_e_lim, ep
!*!
!
! create test plasma
!
        CALL define_plasma(te,ti,ne,nni)     ! eH
!       CALL define_plasma_He(te,ti,ne,nni)  ! He
        ALLOCATE(ceib(1:nni+1))


        CALL write_output(te, ti, ne, nni, betab, zb, mb, nb,        &
          ln_bps_mass, ln_bps_born, delta_mass, cei_sing_mass,       &
          cei_reg_mass, cei_qm_mass, cei_mass, cei_born, cei_dist_ei)


!
! rate C_{eI} using sum-rule and general expression
!
        CALL bps_rate_cei_born(nni,betab,zb,mb,nb,ln_bps_born,cei_born)   ! born
        CALL bps_rate_cei_mass(nni, betab, zb, mb, nb, ln_bps,      & 
          delta, cei_tot, cei_i, cei_e, ceic_tot, ceic_i, ceic_e, ceiq_tot, &
          ceiq_i, ceiq_e, ceic_s_i, ceic_s_e, ceic_r_i , ceic_r_e, ceib)
!
! rate C_{eI} using distribution and A-coefficient integral.
!

        ib=1
        ze=zb(ib)
        me=mb(ib)
        ne=nb(ib)
        betae=betab(ib)
        norme=2*ne*CC*SQRT(betae/me)

!       x0=0
!       x1=30.
!       x1=0.01
!       nmax=1000
!       cei_dist_ei1=0
!       CALL bps_rate_dist(nni, betab, zb, mb, nb, cei_dist_ei, cei_sing_dist_ei,        &
!         cei_reg_dist_ei, cei_qm_dist_ei, x0, x1, nmax)
!       PRINT *, cei_born, cei_born/norme
!       PRINT *, cei_mass, cei_mass/norme
!       PRINT *, cei_dist_ei, cei_dist_ei/norme, (cei_dist_ei-cei_mass)/cei_mass
!       PRINT *, ""
!STOP

        x0=0
        x1=0.2
        x2=1.
        x3=3
        x4=30.
        nmax=1000
        cei_dist_ei1=0
        CALL bps_rate_dist(nni, betab, zb, mb, nb, cei_dist_ei1, cei_sing_dist_ei,        &
          cei_reg_dist_ei, cei_qm_dist_ei, x0, x1, nmax)
        cei_dist_ei2=0
        CALL bps_rate_dist(nni, betab, zb, mb, nb, cei_dist_ei1, cei_sing_dist_ei,        &
          cei_reg_dist_ei, cei_qm_dist_ei, x1, x2, nmax)
        cei_dist_ei3=0
        CALL bps_rate_dist(nni, betab, zb, mb, nb, cei_dist_ei3, cei_sing_dist_ei,        &
          cei_reg_dist_ei, cei_qm_dist_ei, x2, x3, nmax)
        cei_dist_ei4=0
        CALL bps_rate_dist(nni, betab, zb, mb, nb, cei_dist_ei4, cei_sing_dist_ei,        &
          cei_reg_dist_ei, cei_qm_dist_ei, x3, x4, nmax)
        cei_dist_ei=cei_dist_ei1 + cei_dist_ei2 + cei_dist_ei3 + cei_dist_ei4
       PRINT *, cei_born, cei_born/norme
       PRINT *, cei_mass, cei_mass/norme
       PRINT *, cei_dist_ei, cei_dist_ei/norme, (cei_dist_ei-cei_mass)/cei_mass
       PRINT *, ""
!STOP

!STOP

        ib=1
        ze=zb(ib)
        me=mb(ib)
        ne=nb(ib)
        betae=betab(ib)
        norme=2*ne*CC*SQRT(betae/me)

        x0=1.E-10
        x1=0.001
        nmax=500
        dx=(x1-x0)/nmax
        DO ix=0,nmax
           x=x0 + dx*ix
           em=x/betae
           CALL bps_acoeff_ei_mass(nni,em,ze,me,betab,zb,mb,nb,   &
              a_tot,a_i,a_e,ac_tot,ac_i,ac_e,aq_tot,aq_i,aq_e,&
              ac_s_i, ac_s_e, ac_r_i, ac_r_e)
           yei=x*EXP(-x)*a_i*SQRT(2./PI)*1.E7
           WRITE(6,*) ix, x, yei, yei*norme
           WRITE(1,*) ix, x, yei, yei*norme
        ENDDO

        x0=x1 + 0.00001
        x1=1.0
        nmax=2000
        dx=(x1-x0)/nmax
        DO ix=0,nmax
           x=x0 + dx*ix
           em=x/betae
           CALL bps_acoeff_ei_mass(nni,em,ze,me,betab,zb,mb,nb,   &
              a_tot,a_i,a_e,ac_tot,ac_i,ac_e,aq_tot,aq_i,aq_e,&
              ac_s_i, ac_s_e, ac_r_i, ac_r_e)
           yei=x*EXP(-x)*a_i*SQRT(2./PI)*1.E7
           WRITE(6,*) ix+nmax, x, yei, yei*norme
           WRITE(1,*) ix+nmax, x, yei, yei*norme
        ENDDO

        x0=x1 + 0.00001
        x1=30
        nmax=500
        dx=(x1-x0)/nmax
        DO ix=0,nmax
           x=x0 + dx*ix
           em=x/betae
           CALL bps_acoeff_ei_mass(nni,em,ze,me,betab,zb,mb,nb,   &
              a_tot,a_i,a_e,ac_tot,ac_i,ac_e,aq_tot,aq_i,aq_e,&
              ac_s_i, ac_s_e, ac_r_i, ac_r_e)
           yei=x*EXP(-x)*a_i*SQRT(2./PI)*1.E7
           WRITE(6,*) ix+2*nmax, x, yei, yei*norme
           WRITE(1,*) ix+2*nmax, x, yei, yei*norme
        ENDDO
!
! write output header again
!
        CALL write_output(te, ti, ne, nni, betab, zb, mb, nb,        &
          ln_bps_mass, ln_bps_born, delta_mass, cei_sing_mass,       &
          cei_reg_mass, cei_qm_mass, cei_mass, cei_born, cei_dist_ei)

! PRINT *, cei_mass/norme

        CALL close_output
        CALL close_plasma
      END PROGRAM main
!
!======================================================


!
! Returns the plasma species arrays: betab, mb, nb, zb
! Allocates other plasma arrays
!
      SUBROUTINE define_plasma(te, ti, ne, nni)
      USE allocatablevars 
      USE physvars
      USE controlvars
        IMPLICIT NONE
        REAL                                :: te     ! [keV]
        REAL                                :: ti     ! [keV]
        REAL                                :: ne     ! [cm^-3]
        INTEGER                             :: nni    ! number of ion species
!
!       REAL,    DIMENSION(:), ALLOCATABLE  :: betab  ! [1/keV]
!       REAL,    DIMENSION(:), ALLOCATABLE  :: mb     ! [keV]
!       REAL,    DIMENSION(:), ALLOCATABLE  :: nb     ! [cm^-3]
!       REAL,    DIMENSION(:), ALLOCATABLE  :: zb     ! [e]
!       REAL,    DIMENSION(:), ALLOCATABLE  :: gb, etab, mpb, mbpb

!
! open output files
!

        te=100.   ! keV
        ti=100.   ! keV
        ne=1.E25 ! cm^-3

        OPEN(UNIT=1, FILE="main.symmetry.heavyMe.010kev.dat") !
        OPEN(UNIT=1, FILE="main.symmetry.heavyMe.100kev.dat") !

!       OPEN(UNIT=1, FILE="main.symmetry.100kev.dat") !

        asymptotic_forms=.FALSE.
        asymptotic_forms=.TRUE.
        IF (asymptotic_forms) THEN
           small_E_form=.TRUE.
           large_E_form=.TRUE.
        ELSE
           small_E_form=.FALSE.
           large_E_form=.FALSE.
        ENDIF

        nni=1  ! number of ion species
        ALLOCATE(betab(1:nni+1),zb(1:nni+1),mb(nni+1),nb(1:nni+1))   ! allocatablevars
        ALLOCATE(gb(1:nni+1),etab(1:nni+1),mpb(nni+1),mbpb(1:nni+1)) ! allocatablevars

        zb(1)=-1.    ! Species charges
        zb(2)=+1.    ! 
        mb(1)=MEKEV  ! Species masses [keV]
        mb(2)=MPKEV  !
!
! Construct density and temperature arrays
!
        nb(1)=1.                          ! ONLY FOR EQUIMOLAR
        nb(2:nni+1)=1./(zb(2:nni+1)*nni)  ! charge neutrality
        nb=nb*ne                          ! number density array [cm^-3]
        betab(1)=1./te                    ! inverse temp array   [keV^-1]
        betab(2:nni+1)=1./ti              !
      END SUBROUTINE define_plasma

      SUBROUTINE define_plasma_He(te, ti, ne, nni)
      USE allocatablevars 
      USE physvars
      USE controlvars
        IMPLICIT NONE
        REAL                                :: te     ! [keV]
        REAL                                :: ti     ! [keV]
        REAL                                :: ne     ! [cm^-3]
        INTEGER                             :: nni    ! number of ion species
!
!       REAL,    DIMENSION(:), ALLOCATABLE  :: betab  ! [1/keV]
!       REAL,    DIMENSION(:), ALLOCATABLE  :: mb     ! [keV]
!       REAL,    DIMENSION(:), ALLOCATABLE  :: nb     ! [cm^-3]
!       REAL,    DIMENSION(:), ALLOCATABLE  :: zb     ! [e]
!       REAL,    DIMENSION(:), ALLOCATABLE  :: gb, etab, mpb, mbpb

!
! open output files
!

        te=100.   ! keV
        ti=100.   ! keV
        ne=1.E25 ! cm^-3

        OPEN(UNIT=1, FILE="main.symmetry.heavyMe.010kev.dat") !
        OPEN(UNIT=1, FILE="main.symmetry.heavyMe.100kev.dat") !

!       OPEN(UNIT=1, FILE="main.symmetry.100kev.dat") !

        asymptotic_forms=.FALSE.
        asymptotic_forms=.TRUE.
        IF (asymptotic_forms) THEN
           small_E_form=.TRUE.
           large_E_form=.TRUE.
        ELSE
           small_E_form=.FALSE.
           large_E_form=.FALSE.
        ENDIF

        nni=1  ! number of ion species
        ALLOCATE(betab(1:nni+1),zb(1:nni+1),mb(nni+1),nb(1:nni+1))   ! allocatablevars
        ALLOCATE(gb(1:nni+1),etab(1:nni+1),mpb(nni+1),mbpb(1:nni+1)) ! allocatablevars

        zb(2)=1.    ! Species charges
        zb(1)=1.    ! 
        mb(2)=MEKEV  ! Species masses [keV]
        mb(1)=MPKEV  !
!
! Construct density and temperature arrays
!
        nb(1)=1.                          ! ONLY FOR EQUIMOLAR
        nb(2:nni+1)=1./(zb(2:nni+1)*nni)  ! charge neutrality
        nb=nb*ne                          ! number density array [cm^-3]
        betab(1)=1./te                    ! inverse temp array   [keV^-1]
        betab(2:nni+1)=1./ti              !
      END SUBROUTINE define_plasma_He


      SUBROUTINE close_plasma
      USE allocatablevars
        DEALLOCATE(betab,zb,mb,nb)   ! allocatablevars
        DEALLOCATE(gb,etab,mpb)      ! allocatablevars
      END SUBROUTINE close_plasma

      SUBROUTINE write_output(te, ti, ne, nni, betab, zb, mb, nb,   &
      ln_bps_mass, ln_bps_born, delta_mass, cei_sing_mass,        &
      cei_reg_mass, cei_qm_mass, cei_mass, cei_born, cei_dist_ei)
      USE physvars
        IMPLICIT NONE
        REAL,                        INTENT(IN) :: te          ! [keV]
        REAL,                        INTENT(IN) :: ti          ! [keV]
        REAL,                        INTENT(IN) :: ne          ! [cm^-3]
        INTEGER,                     INTENT(IN) :: nni         ! number of ion species
        REAL,    DIMENSION(1:nni+1), INTENT(IN) :: betab       ! [1/keV]
        REAL,    DIMENSION(1:nni+1), INTENT(IN) :: zb          ! [e]
        REAL,    DIMENSION(1:nni+1), INTENT(IN) :: mb          ! [keV]
        REAL,    DIMENSION(1:nni+1), INTENT(IN) :: nb          ! [cm^-3]
        REAL,                        INTENT(IN) :: ln_bps_mass !
        REAL,                        INTENT(IN) :: ln_bps_born !
        REAL,                        INTENT(IN) :: delta_mass  !
        REAL,                        INTENT(IN) :: cei_sing_mass !
        REAL,                        INTENT(IN) :: cei_reg_mass  !
        REAL,                        INTENT(IN) :: cei_qm_mass   !
        REAL,                        INTENT(IN) :: cei_mass      !
        REAL,                        INTENT(IN) :: cei_born      !
        REAL,                        INTENT(IN) :: cei_dist_ei   !
        REAL  :: etae, ge, gi, ze
        REAL,    DIMENSION(1:nni+1) :: gb
        REAL     :: c_ei_born
        REAL     :: c_ei_tot, c_ei_i,c_ei_e, c_eic_tot, c_eic_i, c_eic_e
        REAL     :: c_eiq_tot, c_eiq_i, c_eiq_e, c_eic_s_i, c_eic_s_e
        REAL     :: c_eic_r_i ,c_eic_r_e
        REAL,    DIMENSION(1:nni+1)  :: c_eib
        REAL,    DIMENSION(1:nni+1,1:nni+1)  :: c_ab, c_ab_sing, c_ab_reg, c_ab_qm

!     
! write header
!     
        WRITE(6,'(A)') '#'
        WRITE(6,'(A)') '#'
        WRITE(6,'(A, 3X,A17)') '#','Plasma Parameters'
        WRITE(6,'(A, 4X,A28, X,D12.4, X,D12.4)')  '#','electron and ion temp [keV]:', te, ti
        WRITE(6,'(A, 4X,A28, X,D12.4, X,5D12.4)') '#','number density nb   [cm^-3]:', nb
        WRITE(6,'(A, 4X,A28, X,D12.4, X,5D12.4)') '#','mass array mb         [amu]:', mb/AMUKEV
        WRITE(6,'(A, 4X,A28, X,D12.4, X,5D12.4)') '#','mass array mb         [keV]:', mb
        WRITE(6,'(A, 4X,A28, X,D12.4, X,5D12.4)') '#','charge array zb            :', zb
        WRITE(1,'(A)') '#'  
        WRITE(1,'(A)') '#'
        WRITE(1,'(A, 3X,A17)') '#','Plasma Parameters'
        WRITE(1,'(A, 4X,A28, X,D12.4, X,D12.4)')  '#','electron and ion temp [keV]:', te, ti
        WRITE(1,'(A, 4X,A28, X,D12.4, X,5D12.4)') '#','number density nb   [cm^-3]:', nb
        WRITE(1,'(A, 4X,A28, X,D12.4, X,5D12.4)') '#','mass array mb         [amu]:', mb/AMUKEV
        WRITE(1,'(A, 4X,A28, X,D12.4, X,5D12.4)') '#','mass array mb         [keV]:', mb
        WRITE(1,'(A, 4X,A28, X,D12.4, X,5D12.4)') '#','charge array zb            :', zb
!
! Print plasma parameters
!
        CALL param(nni, betab, nb, gb, ge, gi, etae, ze)
        WRITE(6,'(A1, 2X,A9, 4X,A7, 4X,A7, 8X,A2, 7X,A2, 10X,A4, 9X,A9)') &
          '#','ne[cm^-3]','Te[keV]', 'Ti[keV]','ge','gi','etae','ze/2**1.5'
        WRITE(6,'(A1,11D12.4)') '#',ne, te, ti, ge, gi, etae, ze/2**1.5
        WRITE(1,'(A1, 2X,A9, 4X,A7, 4X,A7, 8X,A2, 7X,A2, 10X,A4, 9X,A9)') &
          '#','ne[cm^-3]','Te[keV]', 'Ti[keV]','ge','gi','etae','ze/2**1.5'
        WRITE(1,'(A1,11D12.4)') '#',ne, te, ti, ge, gi, etae, ze/2**1.5
!
! Print rates
!
        CALL bps_rate_cei_born(nni,betab,zb,mb,nb,ln_bps_born,c_ei_born)
        CALL bps_rate_cei_mass(nni, betab, zb, mb, nb, ln_bps_mass,   & 
          delta_mass, c_ei_tot, c_ei_i,c_ei_e,c_eic_tot,c_eic_i,c_eic_e,    &
          c_eiq_tot,c_eiq_i,c_eiq_e, c_eic_s_i,c_eic_s_e,c_eic_r_i ,c_eic_r_e, c_eib) 
        CALL bps_rate_cab_matrix(nni, betab, zb, mb, nb, &
            c_ab, c_ab_sing, c_ab_reg, c_ab_qm)

        WRITE(6,*) "-----------------------------------"
        WRITE(6,'(A14,E13.5,2X,A10)') "cei_born     :",c_ei_born, "[1/cm^3 s]"
        WRITE(6,'(A14,E13.5,2X,A10)') "cei_mass     :",c_ei_i,    "[1/cm^3 s]"
        WRITE(6,'(A14,E13.5,2X,A10)') "cab_mass_sum :",SUM(c_ab(1,2:nni+1)), "[1/cm^3 s]"
        WRITE(6,'(A14,E13.5)')        "% diff       :",100*ABS((c_ei_i-c_ei_born)/c_ei_i)
        WRITE(6,*) ""
        WRITE(6,'(A14,E13.5)')        "ln_bps_born  :",ln_bps_born
        WRITE(6,'(A14,E13.5)')        "ln_bps_mass  :",ln_bps_mass
        WRITE(6,'(A14,E13.5)')        "% diff       :",100*ABS((ln_bps_mass-ln_bps_born)/ln_bps_mass)
        WRITE(6,*) ""
        WRITE(6,'(A14,E13.5,2X,A10)') "cei_reg_mass :",c_eic_r_i,   "[1/cm^3 s]"
        WRITE(6,'(A14,E13.5,2X,A10)') "cei_sing_mass:",c_eic_s_i,   "[1/cm^3 s]"
        WRITE(6,'(A14,E13.5,2X,A10)') "cei_qm_mass  :",c_eiq_i,     "[1/cm^3 s]"
        WRITE(6,'(A14,E13.5,2X,A10)') "delta_mass   :",delta_mass
        WRITE(6,*) "-----------------------------------"
        WRITE(6,'(A14,E13.5,2X,A10)') "cab_reg_mass :",SUM(c_ab_reg(1,2:nni+1)),   "[1/cm^3 s]"
        WRITE(6,'(A14,E13.5,2X,A10)') "cab_sing_mass:",SUM(c_ab_sing(1,2:nni+1)),  "[1/cm^3 s]"
        WRITE(6,'(A14,E13.5,2X,A10)') "cab_qm_mass  :",SUM(c_ab_qm(1,2:nni+1)),    "[1/cm^3 s]"
        WRITE(6,*) "-----------------------------------"

        WRITE(1,*) "-----------------------------------"
        WRITE(1,'(A14,E13.5,2X,A10)') "cei_born     :",c_ei_born, "[1/cm^3 s]"
        WRITE(1,'(A14,E13.5,2X,A10)') "cei_mass     :",c_ei_i,    "[1/cm^3 s]"
        WRITE(1,'(A14,E13.5,2X,A10)') "cab_mass_sum :",SUM(c_ab(1,2:nni+1)), "[1/cm^3 s]"
        WRITE(1,'(A14,E13.5)')        "% diff       :",100*ABS((c_ei_i-c_ei_born)/c_ei_i)
        WRITE(1,*) ""
        WRITE(1,'(A14,E13.5)')        "ln_bps_born  :",ln_bps_born
        WRITE(1,'(A14,E13.5)')        "ln_bps_mass  :",ln_bps_mass
        WRITE(1,'(A14,E13.5)')        "% diff       :",100*ABS((ln_bps_mass-ln_bps_born)/ln_bps_mass)
        WRITE(1,*) ""
        WRITE(1,'(A14,E13.5,2X,A10)') "cei_reg_mass :",c_eic_r_i,   "[1/cm^3 s]"
        WRITE(1,'(A14,E13.5,2X,A10)') "cei_sing_mass:",c_eic_s_i,   "[1/cm^3 s]"
        WRITE(1,'(A14,E13.5,2X,A10)') "cei_qm_mass  :",c_eiq_i,     "[1/cm^3 s]"
        WRITE(1,'(A14,E13.5,2X,A10)') "delta_mass   :",delta_mass
        WRITE(1,*) "-----------------------------------"
        WRITE(1,'(A14,E13.5,2X,A10)') "cab_reg_mass :",SUM(c_ab_reg(1,2:nni+1)),   "[1/cm^3 s]"
        WRITE(1,'(A14,E13.5,2X,A10)') "cab_sing_mass:",SUM(c_ab_sing(1,2:nni+1)),  "[1/cm^3 s]"
        WRITE(1,'(A14,E13.5,2X,A10)') "cab_qm_mass  :",SUM(c_ab_qm(1,2:nni+1)),    "[1/cm^3 s]"
        WRITE(1,*) "-----------------------------------"


!
! quantities plotted
!
        WRITE(1,'(A77)') "# ix  x  x*EXP(-x)*a_i*SQRT(2./PI)*1.E7  x*EXP(-x)*a_i*norme*SQRT(2./PI)*1.E7"

!       WRITE(6,'(A5, 2X,A7, 5X,A5, 7X,A3, 9X,A3, 9X,A6, 6X,A4, 8X,A4, 8X,A6, 6X,A4, 8X,A4)') &
!         "#  iu", "ep[keV]", "a_tot", "a_i", "a_e", "ac_tot", "ac_i_lim", "ac_e", &
!         "aq_tot", "aq_i", "aq_e"
!       WRITE(1,'(A)') "#"
!       WRITE(1,'(A5, 2X,A7, 5X,A5, 7X,A3, 9X,A3, 9X,A6, 6X,A4, 8X,A4, 8X,A6, 6X,A4, 8X,A4)') &
!         "#  iu", "ep[keV]", "a_tot", "a_i", "a_e", "ac_tot", "ac_i", "ac_e", "aq_tot", "aq_i","aq_e"
      END SUBROUTINE write_output
!
      SUBROUTINE close_output
      IMPLICIT NONE
        CLOSE(1)
      END SUBROUTINE close_output


!
! This routine integrates the A-coefficient weighted by the velocity times
! the distribution function to find the C-coefficient. This method does not
! look symmetric in C_{ab} for a <--> b, although it *must* be symmetric.
!
      SUBROUTINE bps_rate_dist(nni, betab, zb, mb, nb, cei_dist_ei,     &
        cei_sing_dist_ei, cei_reg_dist_ei, cei_qm_dist_ei, xmin, xmax, nmax)
      USE bpsvars
      USE mathvars
      USE physvars
        IMPLICIT NONE
        INTEGER,                     INTENT(IN)  :: nni             !  Number of ions
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: betab           !  Temperature array    [1/keV]
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: zb              !  Charge array 
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: mb              !  Mass array [keV]
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: nb              !  Number density array [cm^-3]
        REAL,                        INTENT(OUT) :: cei_dist_ei     !  Equilibration rate [cm^-3 s^-1]
        REAL,                        INTENT(OUT) :: cei_sing_dist_ei!  Normalized singular
        REAL,                        INTENT(OUT) :: cei_reg_dist_ei !  Normalized regular [cm^-3 s^-1]
        REAL,                        INTENT(OUT) :: cei_qm_dist_ei  !  Normalized quantum
        REAL,                        INTENT(IN)  :: xmin, xmax      !  Integration limits
        INTEGER,                     INTENT(IN)  :: nmax            !  Number of integration regions [even]

        REAL,    PARAMETER :: UPM=0.7745966692E0
        REAL,    PARAMETER :: W13=0.5555555556E0, W2=0.8888888889E0
        REAL               :: a_tot, a_i, a_e, ac_tot, ac_i, ac_e, aq_tot, aq_i, aq_e
        REAL               :: ac_s_i, ac_s_e, ac_r_i, ac_r_e
        REAL               :: x, dx, xm, e, em
        INTEGER            :: ix
        REAL               :: ze, ne, me, betae, norme, yei
        INTEGER            :: ib

        cei_dist_ei     =0.
        cei_sing_dist_ei=0.
        cei_reg_dist_ei =0.
        cei_qm_dist_ei  =0.

        ib=1
        ze=zb(ib)
        me=mb(ib)
        ne=nb(ib)
        betae=betab(ib)
        norme=2*ne*CC*SQRT(betae/me)

        dx=(xmax-xmin)/NMAX
        x=xmin-dx
        DO ix=1,NMAX,2 ! Gaussian quadrature
           x=x+2.*dx
           e=x/betae
           CALL bps_acoeff_ei_mass(nni,e,ze,me,betab,zb,mb,nb,     &
              a_tot,a_i,a_e,ac_tot,ac_i,ac_e,aq_tot,aq_i,aq_e, &
              ac_s_i, ac_s_e, ac_r_i, ac_r_e) ! MeV/micron = 10^7 kev/cm
           yei=x*EXP(-x)*a_i
           cei_dist_ei=cei_dist_ei+W2*yei
!
           xm=x-dx*UPM
           em=xm/betae
           CALL bps_acoeff_ei_mass(nni,em,ze,me,betab,zb,mb,nb,    &
              a_tot,a_i,a_e,ac_tot,ac_i,ac_e,aq_tot,aq_i,aq_e, &
              ac_s_i, ac_s_e, ac_r_i, ac_r_e)
           yei=xm*EXP(-xm)*a_i
           cei_dist_ei=cei_dist_ei+W13*yei
!
           xm=x+dx*UPM
           em=xm/betae
           CALL bps_acoeff_ei_mass(nni,em,ze,me,betab,zb,mb,nb,    &
              a_tot,a_i,a_e,ac_tot,ac_i,ac_e,aq_tot,aq_i,aq_e, &
              ac_s_i, ac_s_e, ac_r_i, ac_r_e)
           yei=xm*EXP(-x)*a_i
           cei_dist_ei=cei_dist_ei+W13*yei
        ENDDO
        cei_dist_ei=cei_dist_ei*1.E7 ! convert A from MeV/micron to keV/cm.
        cei_dist_ei=cei_dist_ei*dx*SQRT(2./PI)
        cei_dist_ei=cei_dist_ei*norme
      END SUBROUTINE bps_rate_dist
