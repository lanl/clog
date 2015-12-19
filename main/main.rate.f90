!
! Rate equations when D and T are at different temperatures.
!
      PROGRAM main
      USE allocatablevars
      USE bpsvars
      USE physvars
      USE mathvars
        IMPLICIT NONE
        REAL     :: te, td, tt, ti, ne
        INTEGER  :: nni

        INTEGER, PARAMETER :: NTSMAX=5000
        REAL    :: y_t(0:NTSMAX,1:2), t_t(0:NTSMAX), y(1:2)
        REAL    :: t, tmin, tmax, dt
        INTEGER :: it, nts, nit
!
! DT plasma with alpha particle projectile
!
       CALL define_plasma(te,td,tt,ne,nni) ! e-DT
!
! Plasma parameters
!
        ti=0.5*(td+tt) ! use average ion temp
        CALL write_output(te,ti,ne,nni,betab,zb,mb,nb)
!
! parameter setup: nts=2000 nit=100
!
        nts=100
        nit=100

        nts=100
        nit=10
        tmin=0
        tmax=1.E-15
        tmax=1.2
        dt=(tmax-tmin)/(nts*nit)
!
! initial conditions
!
        t=tmin
        y(1)=td
        y(2)=tt

        it=0
        t_t(it)=t
        y_t(it,1)=y(1)
        y_t(it,2)=y(2)
!
! evolution
!      

        WRITE(6,*) it, t, y(1), y(2)
        WRITE(1,*) it, t, y(1), y(2)
        DO it=1,nts
           CALL time_step_1(y,t,dt,nit,nni,betab,zb,mb,nb)
           t_t(it)=t
           y_t(it,1)=y(1)
           y_t(it,2)=y(2)
        WRITE(6,*) it, t, y(1), y(2)
        WRITE(1,*) it, t, y(1), y(2)
        ENDDO

        CALL close_output
        CALL close_plasma
      END PROGRAM main
!
! Returns the plasma species arrays: betab, mb, nb, zb
! Allocates other plasma arrays
!
      SUBROUTINE define_plasma(te, td, tt, ne, nni)
      USE allocatablevars 
      USE physvars
        IMPLICIT NONE
        REAL                                :: te     ! [keV]
        REAL                                :: td     ! [keV]
        REAL                                :: tt     ! [keV]
        REAL                                :: ne     ! [cm^-3]
        INTEGER                             :: nni    ! number of ion species

!       REAL,    DIMENSION(:), ALLOCATABLE  :: betab  ! [1/keV]
!       REAL,    DIMENSION(:), ALLOCATABLE  :: mb     ! [keV]
!       REAL,    DIMENSION(:), ALLOCATABLE  :: nb     ! [cm^-3]
!       REAL,    DIMENSION(:), ALLOCATABLE  :: zb     ! [e]
!       REAL,    DIMENSION(:), ALLOCATABLE  :: gb, etab, mpb, mbpb

        te=10.   ! [keV] electron temp
        td=35.   ! [keV] deuteron temp
        tt=30.   ! [keV] triton   temp
        ne=1.E25 ! [cm^-3]

        nni=2  ! number of ion species
        ALLOCATE(betab(1:nni+1),zb(1:nni+1),mb(nni+1),nb(1:nni+1))   ! allocatablevars
        ALLOCATE(gb(1:nni+1),etab(1:nni+1),mpb(nni+1),mbpb(1:nni+1)) ! allocatablevars

        zb(1)=-1.      ! Species charges
        zb(2)=+1.      ! 
        zb(3)=+1.      ! 
        mb(1)=MEKEV    ! Species masses [keV]
        mb(2)=2*MPKEV  !
        mb(3)=3*MPKEV  !

        betab(1)=1./te ! inverse temp array   [keV^-1]
        betab(2)=1./td ! 
        betab(3)=1./tt ! 

!
! Construct density and temperature arrays
!
        nb(1)=1.                          ! ONLY FOR EQUIMOLAR DT
        nb(2:nni+1)=1./(zb(2:nni+1)*nni)  ! charge neutrality
        nb=nb*ne                          ! number density array [cm^-3]
      END SUBROUTINE define_plasma


      SUBROUTINE close_plasma
      USE allocatablevars
        DEALLOCATE(betab,zb,mb,nb)   ! allocatablevars
        DEALLOCATE(gb,etab,mpb)      ! allocatablevars
!       DEALLOCATE(mb,nb,zb,betab,gb)
!       DEALLOCATE(ab,etab)
!       DEALLOCATE(mpb,mbpb)
!       DEALLOCATE(kb2)
      END SUBROUTINE close_plasma

      SUBROUTINE write_output(te,ti,ne,nni,betab,zb,mb,nb)
      USE physvars
        IMPLICIT NONE
        REAL,                        INTENT(IN) :: te     ! [keV]
        REAL,                        INTENT(IN) :: ti     ! [keV]
        REAL,                        INTENT(IN) :: ne     ! [cm^-3]
        INTEGER,                     INTENT(IN) :: nni    ! number of ion species
        REAL,    DIMENSION(1:nni+1), INTENT(IN) :: betab  ! [1/keV]
        REAL,    DIMENSION(1:nni+1), INTENT(IN) :: zb     ! [e]
        REAL,    DIMENSION(1:nni+1), INTENT(IN) :: mb     ! [keV]
        REAL,    DIMENSION(1:nni+1), INTENT(IN) :: nb     ! [cm^-3]
        REAL,    DIMENSION(1:nni+1) :: gb
        REAL     :: etae, ge, gi, ze
        REAL     :: ln_bps_mass, delta_mass
        REAL     :: c_ei_born, ln_bps_born
        REAL     :: c_ei_tot, c_ei_i,c_ei_e, c_eic_tot, c_eic_i, c_eic_e
        REAL     :: c_eiq_tot, c_eiq_i, c_eiq_e, c_eic_s_i, c_eic_s_e
        REAL     :: c_eic_r_i ,c_eic_r_e
        REAL,    DIMENSION(1:nni+1)  :: c_eib
        REAL,    DIMENSION(1:nni+1,1:nni+1)  :: c_ab, c_ab_sing, c_ab_reg, c_ab_qm
!
! open output files
!
        OPEN(UNIT=1, FILE="main.rate.dat")           ! 
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
!      WRITE(6,'(A5, 2X,A7, 5X,A5, 7X,A3, 9X,A3, 9X,A6, 6X,A4, 8X,A4, 8X,A6, 6X,A4, 8X,A4)') &
!        "#  iu", "ep[keV]", "a_tot", "a_i", "a_e", "ac_tot", "ac_i_lim", "ac_e", &
!        "aq_tot", "aq_i", "aq_e"
!      WRITE(1,'(A)') "#"
!      WRITE(1,'(A5, 2X,A7, 5X,A5, 7X,A3, 9X,A3, 9X,A6, 6X,A4, 8X,A4, 8X,A6, 6X,A4, 8X,A4)') &
!        "#  iu", "ep[keV]", "a_tot", "a_i", "a_e", "ac_tot", "ac_i", "ac_e", "aq_tot", "aq_i","aq_e"
      END SUBROUTINE write_output
!
      SUBROUTINE close_output
      IMPLICIT NONE
      CLOSE(1)
      END SUBROUTINE close_output


      SUBROUTINE time_step_1(y,t,dt,nit,nni,betab,zb,mb,nb)
      IMPLICIT NONE
      REAL,    DIMENSION(1:2),     INTENT(INOUT) :: y
      REAL,                        INTENT(INOUT) :: t
      REAL,                        INTENT(INOUT) :: dt
      INTEGER,                     INTENT(IN)    :: nit
      INTEGER,                     INTENT(IN)    :: nni
      REAL,    DIMENSION(1:nni+1), INTENT(IN)    :: betab  !  temp array [1/keV]
      REAL,    DIMENSION(1:nni+1), INTENT(IN)    :: mb     !  mass array [keV]
      REAL,    DIMENSION(1:nni+1), INTENT(IN)    :: nb     !  density [1/cc]
      REAL,    DIMENSION(1:nni+1), INTENT(IN)    :: zb     !  charge array


      REAL,    DIMENSION(1:2)   :: yy, th, f
      REAL,    DIMENSION(1:2)   :: k1, k2, k3, k4
      REAL,    PARAMETER        :: d6=0.1666666666666666666
      INTEGER :: it

      th=t+0.5*dt
      DO it=1,nit
!
! first intermediate step
!
         CALL force_bps_1(y,t,f,nni,betab,zb,mb,nb)
         k1=dt*f
         yy =y+0.5*k1
!
! second intermediate step
!
         CALL force_bps_1(yy,th,f,nni,betab,zb,mb,nb)
         k2=dt*f
         yy =y+0.5*k2
!
! third intermediate step
!
         CALL force_bps_1(yy,th,f,nni,betab,zb,mb,nb)
         k3=dt*f
         yy =y+k3
!
! fourth intermediate step
!
         t=t+dt
         CALL force_bps_1(yy,t,f,nni,betab,zb,mb,nb)
         k4=dt*f
!
         y = y + d6*(k1 + 2*k2 + 2*k3 + k4)
      END DO
      END SUBROUTINE time_step_1

!
! Returns the "force" f(y,t,) for y'=f(y,v,t)
!
      SUBROUTINE force_bps_1(y,t,f,nni,betab,zb,mb,nb)
      IMPLICIT NONE
        REAL,    DIMENSION(1:2),     INTENT(IN)  :: y
        REAL,                        INTENT(IN)  :: t
        REAL,    DIMENSION(1:2),     INTENT(OUT) :: f
        INTEGER,                     INTENT(IN)  :: nni
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: betab  !  temp array [1/keV]
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: mb     !  mass array [keV]
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: nb     !  density [1/cc]
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: zb     !  charge array
        REAL :: c_bps, cv
        REAL,   DIMENSION(1:nni+1,1:nni+1) :: cab
        REAL,   DIMENSION(1:nni+1,1:nni+1) :: cab_sing
        REAL,   DIMENSION(1:nni+1,1:nni+1) :: cab_reg
        REAL,   DIMENSION(1:nni+1,1:nni+1) :: cab_qm


        cv   =1.5*nb(2)
        CALL bps_rate_cab_matrix(nni, betab, zb, mb, nb, &
          cab, cab_sing, cab_reg, cab_qm)
        c_bps=cab(2,3)/cv  ! D-T coupling rate
        c_bps=c_bps*1.E-12 ! convert to ps

        f(1)=-c_bps*(y(1)-y(2))
        f(2)=-c_bps*(y(2)-y(1))
      END SUBROUTINE force_bps_1
