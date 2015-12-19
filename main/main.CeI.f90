!
! Rate equations when D and T are at different temperatures.
!
      PROGRAM main
      USE physvars
      USE mathvars
      USE allocatablevars
        IMPLICIT NONE
        REAL     :: te, td, tt, ti, ne      
        INTEGER  :: nni

        REAL,   DIMENSION(1:3,1:3) :: cab, a_ab
        REAL,   DIMENSION(1:3,1:3) :: cab_sing, a_ab_sing
        REAL,   DIMENSION(1:3,1:3) :: cab_reg, a_ab_reg
        REAL,   DIMENSION(1:3,1:3) :: cab_qm, a_ab_qm

        REAL    :: cint1, cint2
        INTEGER :: ia, ib

        REAL    :: a_tot, a_i, a_e, ac_tot, ac_i, ac_e, aq_tot, aq_i
        REAL    ::  aq_e, ac_s_i, ac_s_e, ac_r_i, ac_r_e, ep, zp, mp
!
! DT plasma with alpha particle projectile
!
       CALL define_plasma(te,td,tt,ne,nni) ! e-DT

        ep=1.
        mp=4*MPKEV
        zp=2.
!
! Plasma parameters
!
        ti=0.5*(td+tt) ! use average ion temp
        CALL write_output(te,ti,ne,nni,betab,zb,mb,nb)
!*!
! ======================================================
!*!
        CALL bps_rate_cab_matrix(nni, betab, zb, mb, nb, &
          cab, cab_sing, cab_reg, cab_qm)

        PRINT *, "cab:"
        PRINT *, cab(1,1:nni+1)
        PRINT *, cab(2,1:nni+1)
        PRINT *, cab(3,1:nni+1)
        PRINT *, ""
!*!
! ======================================================
!*!
        ep=1./betab(1)
        CALL bps_acoeff_ab_matrix(nni,ep,betab,zb,mb,nb,a_ab,a_ab_sing,a_ab_reg,a_ab_qm)
        PRINT *, "A_ab:"
        PRINT *, a_ab(1,1:nni+1)
        PRINT *, a_ab(2,1:nni+1)
        PRINT *, a_ab(3,1:nni+1)
        PRINT *, ""
!*!
! ======================================================
!*!
        ia=1
        ib=2
        zp=zb(ia)
        mp=mb(ia)
        CALL bps_acoeff_ab_matrix(nni,ep,betab,zb,mb,nb,a_ab,a_ab_sing,a_ab_reg,a_ab_qm)
        CALL bps_acoeff_ei_mass(nni, ep, zp, mp, betab, zb, mb, nb, &
            a_tot, a_i, a_e, ac_tot, ac_i, ac_e, aq_tot, aq_i, aq_e, &
            ac_s_i, ac_s_e, ac_r_i, ac_r_e)
        PRINT *, ia, ib, a_tot, a_e, a_i
        PRINT *, ia, ib, a_ab(1,1)+a_ab(1,2)+a_ab(1,3), a_ab(1,1), a_ab(1,2)+a_ab(1,3)
        PRINT *, ""
!*!
! ======================================================
!*!
        ia=2
        ib=3
        PRINT *, ""
        CALL c_from_a_integral(nni,ia,ib,betab,zb,mb,nb,cint1)
        PRINT *, 0, 0, cab(ia,ib)
        PRINT *, ia, ib, cint1
        ia=3
        ib=2
        CALL c_from_a_integral(nni,ia,ib,betab,zb,mb,nb,cint2)
        PRINT *, ia, ib, cint2
        PRINT *, ia, ib, 100*ABS(cint1-cint2)/cint1
        PRINT *, ia, ib, 100*ABS(cab(ia,ib)-cint2)/cint1
        PRINT *, ""

        ia=1
        ib=2
        CALL c_from_a_integral(nni,ia,ib,betab,zb,mb,nb,cint1)
        PRINT *, 0, 0, cab(ia,ib) 
        PRINT *, ia, ib, cint1
        ia=2
        ib=1
        CALL c_from_a_integral(nni,ia,ib,betab,zb,mb,nb,cint2)
        PRINT *, ia, ib, cint2
        PRINT *, ia, ib, 100*ABS(cint1-cint2)/cint1
        PRINT *, ia, ib, 100*ABS(cab(ia,ib)-cint2)/cint1         

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
        td=10.   ! [keV] deuteron temp
        tt=10.   ! [keV] triton   temp
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

!        mb(1)=MEKEV    ! Species masses [keV]
!        mb(2)=2*MEKEV  !
!        mb(3)=3*MEKEV  !


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
        OPEN(UNIT=1, FILE="main.CeI.dat")
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

      SUBROUTINE dedx_fr_acoeff_bps(nni,ep,zp,mp,betab,zb,mb,nb,  &
            dedx_a_tot, dedx_a_i, dedx_a_e, dedxc_a_tot, dedxc_a_i, dedxc_a_e, & 
            dedxq_a_tot, dedxq_a_i, dedxq_a_e, dedxc_a_s_i, dedxc_a_s_e,       &
            dedxc_a_r_i, dedxc_a_r_e)


      USE physvars
      USE mathvars      
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
                                                              ! dE/dx [MeV/micron]
        REAL,                        INTENT(OUT) :: dedx_a_tot  !  electron + ion
        REAL,                        INTENT(OUT) :: dedx_a_i    !  ion contribution
        REAL,                        INTENT(OUT) :: dedx_a_e    !  electron contribution
        REAL,                        INTENT(OUT) :: dedxc_a_tot !  classical
        REAL,                        INTENT(OUT) :: dedxc_a_i   !  classical
        REAL,                        INTENT(OUT) :: dedxc_a_e   !  classical
        REAL,                        INTENT(OUT) :: dedxq_a_tot !  quantum
        REAL,                        INTENT(OUT) :: dedxq_a_i   !  quantum
        REAL,                        INTENT(OUT) :: dedxq_a_e   !  quantum
        REAL,                        INTENT(OUT) :: dedxc_a_s_i
        REAL,                        INTENT(OUT) :: dedxc_a_s_e
        REAL,                        INTENT(OUT) :: dedxc_a_r_i
        REAL,                        INTENT(OUT) :: dedxc_a_r_e

        REAL :: a_tot_p, a_i_p, a_e_p, ac_tot_p, ac_i_p, ac_e_p, aq_tot_p, aq_i_p, aq_e_p
        REAL :: ac_s_i_p, ac_s_e_p, ac_r_i_p, ac_r_e_p
        REAL :: a_tot_m, a_i_m, a_e_m, ac_tot_m, ac_i_m, ac_e_m, aq_tot_m, aq_i_m, aq_e_m
        REAL :: ac_s_i_m, ac_s_e_m, ac_r_i_m, ac_r_e_m
        REAL :: a_tot,   a_i,   a_e,   ac_tot,   ac_i,   ac_e,   aq_tot,   aq_i,   aq_e
        REAL :: ac_s_i,ac_s_e, ac_r_i,ac_r_e
        REAL :: epm, epp, dep, dep2, te, ti

        te  =1./betab(1)
        ti  =1./betab(2)
        dep =ep*1.E-2
        dep2=2*dep

        dedx_a_tot  = 0 ! electron + ion
        dedx_a_i    = 0 ! ion contribution
        dedx_a_e    = 0 ! electron contribution
        dedxc_a_tot = 0 ! classical
        dedxc_a_i   = 0 ! classical
        dedxc_a_e   = 0 ! classical
        dedxq_a_tot = 0 ! quantum
        dedxq_a_i   = 0 ! quantum
        dedxq_a_e   = 0 ! quantum
        dedxc_a_s_i = 0 
        dedxc_a_s_e = 0
        dedxc_a_r_i = 0
        dedxc_a_r_e = 0

        epp=ep+dep
        CALL bps_acoeff_ei_mass(nni,epp,zp,mp,betab,zb,mb,nb,a_tot_p,a_i_p, &
          a_e_p,ac_tot_p,ac_i_p,ac_e_p,aq_tot_p,aq_i_p,aq_e_p,ac_s_i_p,ac_s_e_p,&
          ac_r_i_p,ac_r_e_p)
!*!
!PRINT *, "**1:",ac_i_p, ac_e_p
!STOP
!*!

        epm=ep-dep
        CALL bps_acoeff_ei_mass(nni,epm,zp,mp,betab,zb,mb,nb,a_tot_m,a_i_m, &
          a_e_m,ac_tot_m,ac_i_m,ac_e_m,aq_tot_m,aq_i_m,aq_e_m,ac_s_i_m,ac_s_e_m,&
          ac_r_i_m,ac_r_e_m)

        CALL bps_acoeff_ei_mass(nni,ep,zp,mp,betab,zb,mb,nb,a_tot,a_i, &
          a_e,ac_tot,ac_i,ac_e,aq_tot,aq_i,aq_e,ac_s_i,ac_s_e,&
          ac_r_i,ac_r_e)

!       dedx_tot = (1 - te/ep)*a_tot  - te*(a_tot_p -a_tot_m)/(2*dep)
        dedx_a_i   = (1 - te/ep)*a_i    - ti*(a_i_p   -a_i_m)/dep2
        dedx_a_e   = (1 - te/ep)*a_e    - te*(a_e_p   -a_e_m)/dep2
        dedx_a_tot = dedx_a_e + dedx_a_i

!       dedxc_tot= (1 - te/ep)*ac_tot - te*(ac_tot_p-ac_tot_m)/(2*dep)
        dedxc_a_i  = (1 - te/ep)*ac_i   - ti*(ac_i_p  -ac_i_m)/dep2
        dedxc_a_e  = (1 - te/ep)*ac_e   - te*(ac_e_p  -ac_e_m)/dep2
        dedxc_a_tot= dedxc_a_e + dedxc_a_i

!       dedxq_tot= (1 - te/ep)*aq_tot - te*(aq_tot_p-aq_tot_m)/(2*dep)
        dedxq_a_i  = (1 - ti/ep)*aq_i   - ti*(aq_i_p  -aq_i_m)/dep2
        dedxq_a_e  = (1 - te/ep)*aq_e   - te*(aq_e_p  -aq_e_m)/dep2
        dedxq_a_tot= dedxq_a_e + dedxq_a_i
!
        dedxc_a_s_i=(1 - te/ep)*ac_s_i   - te*(ac_s_i_p  -ac_s_i_m)/dep2
        dedxc_a_s_e=(1 - te/ep)*ac_s_e   - te*(ac_s_e_p  -ac_s_e_m)/dep2
        dedxc_a_r_i=(1 - te/ep)*ac_r_i   - te*(ac_r_i_p  -ac_r_i_m)/dep2
        dedxc_a_r_e=(1 - te/ep)*ac_r_e   - te*(ac_r_e_p  -ac_r_e_m)/dep2


      END SUBROUTINE dedx_fr_acoeff_bps



!*!
! ====================================================================
!*!
      SUBROUTINE c_from_a_integral(nni, ia, ib, betab, zb, mb, nb, cint)
      USE bpsvars
      USE physvars
      USE mathvars
        IMPLICIT NONE
        INTEGER,                     INTENT(IN) :: nni    ! number of ion species
        INTEGER,                     INTENT(IN) :: ia     !  
        INTEGER,                     INTENT(IN) :: ib     !  
        REAL,    DIMENSION(1:nni+1), INTENT(IN) :: betab  ! [1/keV]
        REAL,    DIMENSION(1:nni+1), INTENT(IN) :: zb     ! [e]
        REAL,    DIMENSION(1:nni+1), INTENT(IN) :: mb     ! [keV]
        REAL,    DIMENSION(1:nni+1), INTENT(IN) :: nb     ! [cm^-3]
        REAL,                        INTENT(OUT):: cint   ! [cm^-3 s^-1]
        REAL  :: a_ab_sing
        REAL  :: a_ab_reg, a_ab_qm

        REAL,    PARAMETER :: UPM=0.7745966692E0
        REAL,    PARAMETER :: W13=0.5555555556E0, W2=0.8888888889E0
        REAL    ::xmin, xmax, dx, x, xm, y
        REAL    :: a_ab, ep, mp, zp, np, tp, cnvt
        INTEGER :: ic, nmax
        cint=0

        np=nb(ia)
        zp=zb(ia)
        mp=mb(ia)
        tp=1./betab(ia)

        xmin=0.
        xmax=10. ! automate this choice later.
        nmax=2000
        dx=(xmax-xmin)/nmax 
        x=xmin-dx
        cint=0
        DO ic=1,nmax,2
!     
           x=x+2.E0*dx
           ep=x*tp  ! following subroutine needs to return: A for general ia and ib
           CALL bps_acoeff_ab_mass(nni,ep,mp,zp,ia,ib,betab,zb,mb,nb,a_ab,a_ab_sing,a_ab_reg,a_ab_qm)
           y=a_ab*x*EXP(-x)
           cint=cint + W2*y
!
           xm=x-dx*UPM
           ep=xm*tp
           CALL bps_acoeff_ab_mass(nni,ep,mp,zp,ia,ib,betab,zb,mb,nb,a_ab,a_ab_sing,a_ab_reg,a_ab_qm)
           y=a_ab*xm*EXP(-xm)
           cint=cint + W13*y
!
           xm=x+dx*UPM
           ep=xm*tp
           CALL bps_acoeff_ab_mass(nni,ep,mp,zp,ia,ib,betab,zb,mb,nb,a_ab,a_ab_sing,a_ab_reg,a_ab_qm)
           y=a_ab*xm*EXP(-xm)
           cint=cint + W13*y
        ENDDO
        cint=cint*dx
        cint=SQRT(2./PI)*cint

        cint=cint*2*np*CC/SQRT(mp*tp)
        cnvt=1.E7         ! MeV/micron to keV/cm
        cint=cint*cnvt    ! [cm^-3 s^-1]

      END SUBROUTINE c_from_a_integral

!======================================
!
      SUBROUTINE coeff_bps2(nni, ep, zp, mp, betab, zb, mb, nb, &
            a_tot, a_i, a_e, ac_tot, ac_i, ac_e, aq_tot, aq_i, aq_e)
      USE physvars
      USE mathvars      
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

        REAL,    DIMENSION(1:nni+1)  :: mpb, mbpb, kb2, ab, ab2
        REAL                         :: vp, vp2, zp2, k, k2, kd, kd2, a, b, eta
        REAL                         :: ac_r, ac_s, aq
        REAL                         :: c1, c2
        INTEGER                      :: ib, nnb
!
! initialize components of A-coefficients
!
        a_tot = 0  ! electron + ion
        a_i   = 0  ! ion contribution
        a_e   = 0  ! electron contribution
        ac_tot= 0  ! classical total
        ac_e  = 0  ! classical electron
        ac_i  = 0  ! classical ion
        aq_tot= 0  ! quantum total
        aq_e  = 0  ! quantum electron
        aq_i  = 0  ! quantum ion
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
!       [ betab mbc2 ]^1/2  vp
! c2 =  |----------- |      --  [dimensionless]
!       [   2 Pi     ]      c2
!        
! where
!          e^2 
! Be  =  --------- = 13.606E-3   [keV]
!        8 Pi a0                 Bohr radius: a0=5.29E9 cm
!
        NNB = nni+1                 ! number of ions + electrons
        vp  = CC*SQRT(2*ep/mp)      ! [cm/s]
        vp2 = vp*vp                 ! [cm^2/s^2]
        kb2 = DEBYE2*zb*zb*nb*betab ! [1/cm^2]
        kd2 = SUM(kb2)              ! [1/cm^2]
        kd  = SQRT(kd2)             ! [1/cm]
        k2  = kb2(1)                ! [1/cm^2]
        k   = SQRT(k2)              ! [1/cm]   k = k_e
        zp2 = zp**2                 ! [dimensionless]
                                    ! ab=(1/2) betab(ib)*mbc2(ib)*vp2/CC2
        ab  =0.5*betab*mb*vp2/CC2   ! [dimensionless] 
        ab2 = SQRT(ab)              ! [dimensionless]
        mpb = mp*mb/(mp+mb)         ! [keV]
        mbpb= mb/mpb                ! [dimensionless]
!
! Loop over charged plasma species
!
        DO ib=1,nni+1
        IF (zb(ib) .NE. 0.) THEN
           a  =ab(ib)
           b  =-Log( 2*betab(ib)*BEKEV* ABS(zp*zb(ib))*k*A0CM*mbpb(ib) ) - & 
             2*GAMMA + 2
           eta=ABS(zp*zb(ib))*2.1870E8/vp ! defined with projectile velocity vp

           c1=2*zp2*BEKEV*kb2(ib)*A0CM ! [keV/cm]
           c1=c1*1.E-7               ! [MeV/micron]  
           c2=SQRT(a/PI)             ! [dimensionless] 
                                     ! c2=SQRT(betab(ib)*mb(ib)/TWOPI)*vp/CC 
!
! A-classical-singular 
!
           CALL a_sing2(a,b,ac_s) 
           ac_s=c1*c2*ac_s
!
! A-classical-regular 
!
           CALL a_reg2(ib,nni,vp,k2,kb2,betab,mb,ac_r)
           ac_r=c1*ac_r
!
! A-quantum
!
           CALL a_quantum2(ib,a,eta,aq) ! eta = dimensionless quantum param.
           aq=c1*c2*aq
!
! construct electron and ion components
!
           CALL a_collect(ib,NNB,ac_s,ac_r,aq,a_tot,a_e,a_i,&
             ac_tot,ac_e,ac_i,aq_tot,aq_e,aq_i)
        ENDIF
        ENDDO
      END SUBROUTINE coeff_bps2

!======================================
! Singular contribution: short distance
!
      SUBROUTINE a_sing2(a, b, ac_s)
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
              ac_s=ac_s+W2*dab_sing2(u,a,b)
              um=u-du*UPM
              ac_s=ac_s+W13*dab_sing2(um,a,b)
              um=u+du*UPM
              ac_s=ac_s+W13*dab_sing2(um,a,b)
           ENDDO
           ac_s=ac_s*du
      END SUBROUTINE a_sing2
!
      FUNCTION dab_sing2(u, a, b)
        IMPLICIT NONE
        REAL,        INTENT(IN)  :: u        ! [dimensionless]
        REAL,        INTENT(IN)  :: a        ! [dimensionless] 
                                             ! a=(1/2)*beta*mpc2*vp^2/C^2
        REAL,        INTENT(IN)  :: b        ! [dimensionless]
        REAL                     :: dab_sing2 ! [dimensionless]
        dab_sing2=SQRT(u)*EXP(-a*u)*(-LOG(u/(1-u)) + b)
      END FUNCTION dab_sing2
!
!======================================


!======================================
! Regular contribution: long distance
!
      SUBROUTINE a_reg2(ib, nni, vp, k2, kb2, betab, mb, ac_r)
        INTEGER,                     INTENT(IN)  :: ib
        INTEGER,                     INTENT(IN)  :: nni 
        REAL,                        INTENT(IN)  :: k2
        REAL,                        INTENT(IN)  :: kb2
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: betab
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: mb
        REAL,                        INTENT(OUT) :: ac_r
        REAL               :: u0, u1, du, u, um
        INTEGER, PARAMETER :: NR=10 ! integration regions: must be even
        REAL,    PARAMETER :: UPM=0.7745966692E0 ! parameters for Gaussian Quad
        REAL,    PARAMETER :: W13=0.5555555556E0, W2=0.8888888889E0
        ac_r=0
        u0=0.
        u1=1.
        du=(u1-u0)/NR
        u=u0-du
        DO iu=1,NR,2 ! Gaussian quadrature
           u=u+2.E0*du
           ac_r=ac_r+W2*dab_reg2(u,vp,ib,nni,k2,kb2,betab,mb)
           um=u-du*UPM
           ac_r=ac_r+W13*dab_reg2(um,vp,ib,nni,k2,kb2,betab,mb)
           um=u+du*UPM
           ac_r=ac_r+W13*dab_reg2(um,vp,ib,nni,k2,kb2,betab,mb)
        ENDDO
        ac_r=ac_r*du
      END SUBROUTINE a_reg2
!
      FUNCTION dab_reg2(u, vp, ib, nni, k2, kb2, betab, mb)
      USE mathvars
      USE physvars
        IMPLICIT NONE
        REAL,                        INTENT(IN)  :: u       !  [dimensionless]
        REAL,                        INTENT(IN)  :: vp      !  Projectile velocity [cm/s]
        INTEGER,                     INTENT(IN)  :: ib      !  Species number
        INTEGER,                     INTENT(IN)  :: nni     !  Number of ion species
        REAL,                        INTENT(IN)  :: k2      !  Wave-number squared [1/cm^2]
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: kb2     !  Debye wavenumber squared [1/cm^2]
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: betab   !  Temperature array [1/keV]
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: mb      !  Mass array [keV]
        REAL                                     :: dab_reg2 !  [dimensionless]
        REAL,    DIMENSION(1:nni+1) :: alfb, ab
        REAL                        :: fr, fi, fabs, farg, h
        REAL                        :: kcb, r_ib, bm_ic, bm_ib, a_ic, a_ib, ex, au
        INTEGER                     :: ic
        ab=SQRT(0.5*betab*mb)*vp/CC
        alfb=kb2/k2
        CALL frfi(u,nni,alfb,ab,fr,fi,fabs,farg)
        h=2*(fr*farg + fi*LOG(fabs))
!
! construct spectral weight ratio Rb=rho_b/rho_tot
!
        r_ib=0
        bm_ib=betab(ib)*mb(ib)
        a_ib =ab(ib)*ab(ib)
        DO ic=1,nni+1
           kcb=kb2(ic)/k2
           bm_ic=betab(ic)*mb(ic)
           a_ic =ab(ic)*ab(ic)
           IF (ic == ib) THEN
              ex=1
           ELSE
              au=(a_ic-a_ib)*u
              ex=EXP(-au)
           ENDIF
           r_ib=r_ib + kcb*SQRT(bm_ic/bm_ib)*ex
        ENDDO      
        r_ib=1/r_ib
        dab_reg2=-u*r_ib*h/TWOPI
      END FUNCTION dab_reg2
!
!======================================

!======================================
!
      SUBROUTINE a_quantum2(ib, a, eta, aq)
        IMPLICIT NONE
        INTEGER, INTENT(IN)  :: ib    ! species index
        REAL,    INTENT(IN)  :: a     ! [dimensionless] (1/2) betab mb vp^2
        REAL,    INTENT(IN)  :: eta   ! [dimensionless] ep eb/4pi hbar vp
        REAL,    INTENT(OUT) :: aq 
        REAL               :: u0, u1, du, u, um
        INTEGER, PARAMETER :: NQ=1000            ! integration regions quantum : must be even
        REAL,    PARAMETER :: UPM=0.7745966692E0 ! parameters for Gaussian Quad
        REAL,    PARAMETER :: W13=0.5555555556E0, W2=0.8888888889E0
        REAL    :: daq2
        INTEGER :: iu
        aq=0
        u0=0.
        aq=0
        IF (ib == 1) THEN
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
           aq=aq+W2*daq2(u,a,eta)
           um=u-du*UPM
           aq=aq+W13*daq2(um,a,eta)
           um=u+du*UPM
           aq=aq+W13*daq2(um,a,eta)
        ENDDO
        aq=aq*du
      END SUBROUTINE a_quantum2
!
      FUNCTION daq2(u, a, eta)
      USE physvars
        IMPLICIT NONE
        REAL,                        INTENT(IN)  :: u          ! [dimensionless]
        REAL,                        INTENT(IN)  :: a          ! [dimensionless]
        REAL,                        INTENT(IN)  :: eta        ! [dimensionless]
        REAL                                     :: daq2  ! [dimensionless]
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
        daq2=-psilog*csh
      END FUNCTION daq2
!======================================
