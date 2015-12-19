!
! A-coeff in the format:  n  Ep[keV]  A_tot[Mev/micron]  A_e  A_I
! 
!
!
    PROGRAM main
    USE allocatablevars
    USE bpsvars
      USE physvars
      USE mathvars
        IMPLICIT NONE
        REAL     :: te, ti, ne, ep, mp, zp

        REAL     :: a_tot, a_e, a_i, ac_tot, ac_e, ac_i, aq_tot, aq_e, aq_i
        REAL     :: ac_s_i, ac_s_e, ac_r_i, ac_r_e
        REAL     :: ep0, ep1, ep2, ep3, dep1, dep2, dep3
        INTEGER  :: iu, nep1, nep2, nep3


        REAL     :: a_tot_lim, a_e_lim, a_i_lim, ac_tot_lim, ac_e_lim, ac_i_lim
        REAL     :: aq_tot_lim, aq_e_lim, aq_i_lim
        REAL     :: ac_s_i_lim, ac_s_e_lim, ac_r_i_lim, ac_r_e_lim

        REAL     :: nni
!
! Plasma parameters
!
        te=10.       ! Electron temperature       [keV]
        ti=10.       ! Ion temperature            [keV]
        ne=1.E25     ! Electron number density    [cm^-3]
!
! Projectile parameters (alpha particle)
!
        ep=1.        ! Projectile energy  [keV]
        mp=4*MPKEV   ! Projectile mass    [keV]
        zp=2.        ! Projectile charge  [e]
!
! DT plasma with alpha particle projectile
!
        CALL define_plasma(te,ti,ne,nni)
!
! plot the regular and singular contributions
!
        CALL write_output(ep,mp,zp,te,ti,ne,nni,betab,zb,mb,nb)

        ep0=1.E-7   ! [keV]
        ep1=    50. ! [keV]
        ep2=  3500. ! [keV] approx intermediate high energy limit
        ep3=500000. ! [keV] extreme high energy limit

        nep1=1000
        nep2=1000
        nep3=1000
!*!
nep1=10
nep2=10
nep3=10
!*!
        dep1=(ep1-ep0)/nep1
        dep2=(ep2-ep1)/nep2
        dep3=(ep3-ep2)/nep3
!
! first region [ep0,ep1]
!
        DO iu=0,nep1
           ep=ep0 + dep1*iu
           CALL bps_acoeff_ei_mass(nni, ep, zp, mp, betab, zb, mb, nb, &
               a_tot, a_i, a_e, ac_tot, ac_i, ac_e, aq_tot, aq_i, aq_e, &
               ac_s_i, ac_s_e, ac_r_i, ac_r_e)
           WRITE(6, '(I5,10E12.4)') iu, ep, a_tot, a_e, a_i, ac_tot, ac_e, ac_i, aq_tot, aq_e, aq_i
           WRITE(1, '(I5,10E12.4)') iu, ep, a_tot, a_e, a_i, ac_tot, ac_e, ac_i, aq_tot, aq_e, aq_i

           ! asymptotic small-E form
           CALL coeff_bps_small_E(nni, ep, zp, mp, betab, zb, mb, nb,     &
             a_tot_lim, a_i_lim, a_e_lim, ac_tot_lim, ac_i_lim, ac_e_lim, &
             aq_tot_lim, aq_i_lim, aq_e_lim,ac_s_i_lim, ac_s_e_lim,      &
             ac_r_i_lim, ac_r_e_lim)
           WRITE(2, '(I5,10E12.4)') iu, ep, a_tot_lim, a_e_lim, a_i_lim, ac_tot_lim, ac_e_lim, &
             ac_i_lim, aq_tot_lim, aq_e_lim, aq_i_lim

           ! asymptotic large-E form
           CALL coeff_bps_high_E(nni, ep, zp, mp, betab, zb, mb, nb,     &
            a_tot_lim, a_i_lim, a_e_lim, ac_tot_lim, ac_i_lim, ac_e_lim, &
            aq_tot_lim, aq_i_lim, aq_e_lim)
           WRITE(3, '(I5,10E12.4)') iu, ep, a_tot_lim, a_e_lim, a_i_lim, ac_tot_lim, ac_e_lim, &
              ac_i_lim, aq_tot_lim, aq_e_lim, aq_i_lim
        ENDDO
!
! second region [ep1,ep2]
!
        DO iu=1,nep2
           ep=ep1 + dep2*iu
           CALL bps_acoeff_ei_mass(nni, ep, zp, mp, betab, zb, mb, nb, &
               a_tot, a_i, a_e, ac_tot, ac_i, ac_e, aq_tot, aq_i, aq_e,&
               ac_s_i, ac_s_e, ac_r_i, ac_r_e)
           WRITE(6, '(I5,10E12.4)') iu+nep1, ep, a_tot, a_e, a_i, ac_tot, ac_e, ac_i, aq_tot, aq_e, aq_i
           WRITE(1, '(I5,10E12.4)') iu+nep1, ep, a_tot, a_e, a_i, ac_tot, ac_e, ac_i, aq_tot, aq_e, aq_i

           ! asymptotic small-E form
           CALL coeff_bps_small_E(nni, ep, zp, mp, betab, zb, mb, nb,       &
             a_tot_lim, a_i_lim, a_e_lim, ac_tot_lim, ac_i_lim, ac_e_lim,   &
             aq_tot_lim, aq_i_lim, aq_e_lim,ac_s_i_lim, ac_s_e_lim,      &
            ac_r_i_lim, ac_r_e_lim)
           WRITE(2, '(I5,10E12.4)') iu+nep1, ep, a_tot_lim, a_e_lim, a_i_lim, ac_tot_lim, ac_e_lim, &
              ac_i_lim, aq_tot_lim, aq_e_lim, aq_i_lim

           ! asymptotic large-E form
           CALL coeff_bps_high_E(nni, ep, zp, mp, betab, zb, mb, nb,     &
            a_tot_lim, a_i_lim, a_e_lim, ac_tot_lim, ac_i_lim, ac_e_lim, &
            aq_tot_lim, aq_i_lim, aq_e_lim)
           WRITE(3, '(I5,10E12.4)') iu, ep, a_tot_lim, a_e_lim, a_i_lim, ac_tot_lim, ac_e_lim, &
              ac_i_lim, aq_tot_lim, aq_e_lim, aq_i_lim
        ENDDO
!
! third region [ep2,ep3]
!
        DO iu=1,nep3
           ep=ep2 + dep3*iu
           CALL bps_acoeff_ei_mass(nni, ep, zp, mp, betab, zb, mb, nb, &
               a_tot, a_i, a_e, ac_tot, ac_i, ac_e, aq_tot, aq_i, aq_e, &
               ac_s_i, ac_s_e, ac_r_i, ac_r_e)
           WRITE(6, '(I5,10E12.4)') iu+nep1+nep2, ep, a_tot, a_e, a_i, ac_tot, ac_e, ac_i, &
              aq_tot, aq_e, aq_i
           WRITE(1, '(I5,10E12.4)') iu+nep1+nep2, ep, a_tot, a_e, a_i, ac_tot, ac_e, ac_i, &
              aq_tot, aq_e, aq_i

           ! asymptotic small-E form
           CALL coeff_bps_small_E(nni, ep, zp, mp, betab, zb, mb, nb,       &
             a_tot_lim, a_i_lim, a_e_lim, ac_tot_lim, ac_i_lim, ac_e_lim,   &
             aq_tot_lim, aq_i_lim, aq_e_lim,ac_s_i_lim, ac_s_e_lim,      &
            ac_r_i_lim, ac_r_e_lim)
           WRITE(2, '(I5,10E12.4)') iu+nep1+nep2, ep, a_tot_lim, a_e_lim, a_i_lim, &
              ac_tot_lim, ac_e_lim, ac_i_lim, aq_tot_lim, aq_e_lim, aq_i_lim

           ! asymptotic large-E form
           CALL coeff_bps_high_E(nni, ep, zp, mp, betab, zb, mb, nb,     &
            a_tot_lim, a_i_lim, a_e_lim, ac_tot_lim, ac_i_lim, ac_e_lim, &
            aq_tot_lim, aq_i_lim, aq_e_lim)
           WRITE(3, '(I5,10E12.4)') iu, ep, a_tot_lim, a_e_lim, a_i_lim, ac_tot_lim, ac_e_lim, &
              ac_i_lim, aq_tot_lim, aq_e_lim, aq_i_lim
        ENDDO


! 
! list format: { {x,y,z}, {x,y,x}, .... , {x,y,z} }
!
!           IF (iu == 0) THEN
!              WRITE(1, '(A2,I5,A1,E12.4,A1,E12.4,A1,E12.4,A1,E12.4,A2)') & 
!                '{{', iu, ',', ep, ',', a_tot, ',', a_e, ',', a_i, '},'
!           ELSEif (iu == nep) THEN
!              WRITE(1, '(A1,I5,A1,E12.4,A1,E12.4,A1,E12.4,A1,E12.4,A2)') & 
!                '{', iu, ',', ep, ',', a_tot, ',', a_e, ',', a_i, '}}'
!           ELSE
!              WRITE(1, '(A1,I5,A1,E12.4,A1,E12.4,A1,E12.4,A1,E12.4,A2)') & 
!                '{', iu, ',', ep, ',', a_tot, ',', a_e, ',', a_i, '},'
!           ENDIF


        CALL close_output
        CALL close_plasma
      END PROGRAM main

!
! Returns the plasma species arrays: betab, mb, nb, zb
! Allocates other plasma arrays
!
    SUBROUTINE define_plasma(te, ti, ne, nni)
    USE allocatablevars 
    USE physvars
      IMPLICIT NONE
      REAL                                :: te     ! [keV]
      REAL                                :: ti     ! [keV]
      REAL                                :: ne     ! [cm^-3]
      INTEGER                             :: nni    ! number of ion species
!
!     REAL,    DIMENSION(:), ALLOCATABLE  :: betab  ! [1/keV]
!     REAL,    DIMENSION(:), ALLOCATABLE  :: mb     ! [keV]
!     REAL,    DIMENSION(:), ALLOCATABLE  :: nb     ! [cm^-3]
!     REAL,    DIMENSION(:), ALLOCATABLE  :: zb     ! [e]
!     REAL,    DIMENSION(:), ALLOCATABLE  :: gb, etab, mpb, mbpb

      nni=2  ! number of ion species

      ALLOCATE(betab(1:nni+1),zb(1:nni+1),mb(nni+1),nb(1:nni+1))   ! allocatablevars
      ALLOCATE(gb(1:nni+1),etab(1:nni+1),mpb(nni+1),mbpb(1:nni+1)) ! allocatablevars

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

    SUBROUTINE close_plasma
    USE allocatablevars
      DEALLOCATE(betab,zb,mb,nb)   ! allocatablevars
      DEALLOCATE(gb,etab,mpb)      ! allocatablevars
!     DEALLOCATE(mb,nb,zb,betab,gb)
!     DEALLOCATE(ab,etab)
!     DEALLOCATE(mpb,mbpb)
!     DEALLOCATE(kb2)
    END SUBROUTINE close_plasma

    SUBROUTINE write_output(ep, mp, zp, te, ti, ne, nni, betab, zb, mb, nb)
    USE physvars
      IMPLICIT NONE
      REAL,                        INTENT(IN) :: ep     ! [keV]
      REAL,                        INTENT(IN) :: mp     ! [keV]
      REAL,                        INTENT(IN) :: zp     ! [e]
      REAL,                        INTENT(IN) :: te     ! [keV]
      REAL,                        INTENT(IN) :: ti     ! [keV]
      REAL,                        INTENT(IN) :: ne     ! [cm^-3]
      INTEGER,                     INTENT(IN) :: nni    ! number of ion species
      REAL,    DIMENSION(1:nni+1), INTENT(IN) :: betab  ! [1/keV]
      REAL,    DIMENSION(1:nni+1), INTENT(IN) :: zb     ! [e]
      REAL,    DIMENSION(1:nni+1), INTENT(IN) :: mb     ! [keV]
      REAL,    DIMENSION(1:nni+1), INTENT(IN) :: nb     ! [cm^-3]


      REAL  :: etae, ge, gi, ze
      REAL,    DIMENSION(1:nni+1) :: gb 

      REAL  :: vp
      REAL,    DIMENSION(1:nni+1) :: ab, etab
!
! open output files
!
      OPEN(UNIT=1, FILE="gr001.dat")        ! 
      OPEN(UNIT=2, FILE="gr001.smallE.dat") ! 
      OPEN(UNIT=3, FILE="gr001.highE.dat")  ! 

!
! write header
!
      WRITE(6,'(A)') '#'
      WRITE(6,'(A, 3X,A21)') '#','Projectile Parameters'
      WRITE(6,'(A, 4X,A13, X,D12.4)') '#','charge      :', zp
      WRITE(6,'(A, 4X,A13, X,D12.4)') '#','mass   [amu]:', mp/AMUKEV
      WRITE(6,'(A, 4X,A13, X,D12.4)') '#','mass   [keV]:', mp
      WRITE(6,'(A, 4X,A13, X,D12.4)') '#','energy [keV]:', ep
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
      WRITE(1,'(A)') '#'
      WRITE(1,'(A)') '#'
      WRITE(1,'(A, 3X,A21)') '#','Projectile Parameters'
      WRITE(1,'(A, 4X,A13, X,D12.4)') '#','charge      :', zp
      WRITE(1,'(A, 4X,A13, X,D12.4)') '#','mass   [amu]:', mp/AMUKEV
      WRITE(1,'(A, 4X,A13, X,D12.4)') '#','mass   [keV]:', mp
      WRITE(1,'(A, 4X,A13, X,D12.4)') '#','energy [keV]:', ep

      WRITE(2,'(A)') '#'  
      WRITE(2,'(A)') '#'
      WRITE(2,'(A, 3X,A17)') '#','Plasma Parameters'
      WRITE(2,'(A, 4X,A28, X,D12.4, X,D12.4)')  '#','electron and ion temp [keV]:', te, ti
      WRITE(2,'(A, 4X,A28, X,D12.4, X,5D12.4)') '#','number density nb   [cm^-3]:', nb
      WRITE(2,'(A, 4X,A28, X,D12.4, X,5D12.4)') '#','mass array mb         [amu]:', mb/AMUKEV
      WRITE(2,'(A, 4X,A28, X,D12.4, X,5D12.4)') '#','mass array mb         [keV]:', mb
      WRITE(2,'(A, 4X,A28, X,D12.4, X,5D12.4)') '#','charge array zb            :', zb
      WRITE(2,'(A)') '#'
      WRITE(2,'(A)') '#'
      WRITE(2,'(A, 3X,A21)') '#','Projectile Parameters'
      WRITE(2,'(A, 4X,A13, X,D12.4)') '#','charge      :', zp
      WRITE(2,'(A, 4X,A13, X,D12.4)') '#','mass   [amu]:', mp/AMUKEV
      WRITE(2,'(A, 4X,A13, X,D12.4)') '#','mass   [keV]:', mp
      WRITE(2,'(A, 4X,A13, X,D12.4)') '#','energy [keV]:', ep

      WRITE(3,'(A)') '#'  
      WRITE(3,'(A)') '#'
      WRITE(3,'(A, 3X,A17)') '#','Plasma Parameters'
      WRITE(3,'(A, 4X,A28, X,D12.4, X,D12.4)')  '#','electron and ion temp [keV]:', te, ti
      WRITE(3,'(A, 4X,A28, X,D12.4, X,5D12.4)') '#','number density nb   [cm^-3]:', nb
      WRITE(3,'(A, 4X,A28, X,D12.4, X,5D12.4)') '#','mass array mb         [amu]:', mb/AMUKEV
      WRITE(3,'(A, 4X,A28, X,D12.4, X,5D12.4)') '#','mass array mb         [keV]:', mb
      WRITE(3,'(A, 4X,A28, X,D12.4, X,5D12.4)') '#','charge array zb            :', zb
      WRITE(3,'(A)') '#'
      WRITE(3,'(A)') '#'
      WRITE(3,'(A, 3X,A21)') '#','Projectile Parameters'
      WRITE(3,'(A, 4X,A13, X,D12.4)') '#','charge      :', zp
      WRITE(3,'(A, 4X,A13, X,D12.4)') '#','mass   [amu]:', mp/AMUKEV
      WRITE(3,'(A, 4X,A13, X,D12.4)') '#','mass   [keV]:', mp
      WRITE(3,'(A, 4X,A13, X,D12.4)') '#','energy [keV]:', ep

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
      WRITE(2,'(A1, 2X,A9, 4X,A7, 4X,A7, 8X,A2, 7X,A2, 10X,A4, 9X,A9)') &
        '#','ne[cm^-3]','Te[keV]', 'Ti[keV]','ge','gi','etae','ze/2**1.5'
      WRITE(2,'(A1,11D12.4)') '#',ne, te, ti, ge, gi, etae, ze/2**1.5
      WRITE(3,'(A1, 2X,A9, 4X,A7, 4X,A7, 8X,A2, 7X,A2, 10X,A4, 9X,A9)') &
        '#','ne[cm^-3]','Te[keV]', 'Ti[keV]','ge','gi','etae','ze/2**1.5'
      WRITE(3,'(A1,11D12.4)') '#',ne, te, ti, ge, gi, etae, ze/2**1.5


      vp  = CC*SQRT(2*ep/mp)       ! [cm/s]
      ab  =0.5*betab*mb*vp*vp/CC2  ! [dimensionless] (1/2) betab(ib)*mb(ib)*vp2/CC2
      etab=ABS(zp*zb)*2.1870E8/vp  ! [dimensionless]

      WRITE(6,'(A7, 3D12.4)') '# ab  =', ab
      WRITE(6,'(A7, 3D12.4)') '# etab=', etab
      WRITE(1,'(A7, 3D12.4)') '# ab  =', ab
      WRITE(1,'(A7, 3D12.4)') '# etab=', etab
      WRITE(2,'(A7, 3D12.4)') '# ab  =', ab
      WRITE(2,'(A7, 3D12.4)') '# etab=', etab
      WRITE(3,'(A7, 3D12.4)') '# ab  =', ab
      WRITE(3,'(A7, 3D12.4)') '# etab=', etab

      WRITE(6,'(A)') "#"
      WRITE(6,'(A5, 2X,A7, 5X,A5, 7X,A3, 9X,A3, 9X,A6, 6X,A4, 8X,A4, 8X,A6, 6X,A4, 8X,A4)') &
        "#  iu", "ep[keV]", "a_tot", "a_e", "a_i", "ac_tot", "ac_e_lim", "ac_i", &
        "aq_tot", "aq_e", "aq_i"
      WRITE(1,'(A)') "#"
      WRITE(1,'(A5, 2X,A7, 5X,A5, 7X,A3, 9X,A3, 9X,A6, 6X,A4, 8X,A4, 8X,A6, 6X,A4, 8X,A4)') &
        "#  iu", "ep[keV]", "a_tot", "a_e", "a_i", "ac_tot", "ac_e", "ac_i", "aq_tot", "aq_e","aq_i"
      WRITE(2,'(A)') "#"
      WRITE(2,'(A5, 2X,A7, 5X,A9, 3X,A7, 5X,A7, 5X,A10, 2X,A8, 4X,A8, 4X,A10, 2X,A8, 4X,A8)') &
        "#  iu", "ep[keV]", "a_tot_lim", "a_e_lim", "a_i_lim", "ac_tot_lim", "ac_e_lim", "ac_i_lim", &
        "aq_tot_lim", "aq_e_lim", "aq_i_lim"
      WRITE(3,'(A)') "#"
      WRITE(3,'(A5, 2X,A7, 5X,A9, 3X,A7, 5X,A7, 5X,A10, 2X,A8, 4X,A8, 4X,A10, 2X,A8, 4X,A8)') &
        "#  iu", "ep[keV]", "a_tot_lim", "a_e_lim", "a_i_lim", "ac_tot_lim", "ac_e_lim", "ac_i_lim", &
        "aq_tot_lim", "aq_e_lim", "aq_i_lim"

    END SUBROUTINE write_output

    SUBROUTINE close_output
    IMPLICIT NONE
      CLOSE(1)
      CLOSE(3)
    END SUBROUTINE close_output
