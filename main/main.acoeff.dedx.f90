!
! This is a driver to check the analytic evalulation dE_b/dx
! against the one obtained by differentiating A_b.
!
! dE_b        [            1                d                           ]
! ----(vp) =  [ 1  -  ------------- Sum_l -------  {\hat vp}^l A_b(vp)  ]
!  dx         [        beta_b mp vp       d vp^l                        ]
!
!
!             [          2 T_b   ]                T_b     dA_b
!          =  [ 1  -  ---------- ] * A_b(vp) -  ------- * ----(vp)
!             [         mp vp^2  ]               mp vp     dvp
!
!             [       T_b  ]                    dA_b
!          =  [ 1  -  ---- ] * A_b(Ep) -  T_b * ----(Ep)
!             [        Ep  ]                     dEp
!
!
! DT plasma with alpha particle projectile
!
    PROGRAM main
    USE physvars
    USE mathvars
    USE allocatablevars
        IMPLICIT NONE
        REAL     :: ep, zp, mp
        REAL     :: te, ti, ne
        INTEGER  :: nni

        REAL     :: dedx_tot, dedx_i, dedx_e, dedxc_tot, dedxc_i, dedxc_e
        REAL     :: dedxq_tot, dedxq_i, dedxq_e

        REAL     :: dedx_a_tot, dedx_a_i, dedx_a_e, dedxc_a_tot, dedxc_a_i 
        REAL     :: dedxc_a_e, dedxq_a_tot, dedxq_a_i, dedxq_a_e
        REAL     :: dedxc_a_s_i, dedxc_a_s_e, dedxc_a_r_i, dedxc_a_r_e

        REAL     :: diff_dedx_tot, diff_dedx_i, diff_dedx_e, diff_dedxc_tot
        REAL     :: diff_dedxc_i, diff_dedxc_e, diff_dedxq_tot, diff_dedxq_i  
        REAL     :: diff_dedxq_e 


        REAL     :: dep, ep0, ep1
        INTEGER  :: iu, nep

!
! DT plasma with alpha particle projectile
!
        CALL define_plasma(ep,zp,mp,te,ti,ne,nni)
!
! Plasma parameters
!
        ep0=1.E-7   ! [keV]
        ep1=3500.   ! [keV]
        nep=1000
        dep=(ep1-ep0)/nep
        CALL write_output(nep+1,ep,zp,mp,te,ti,ne,nni,betab,zb,mb,nb)



      iu=0
      WRITE(6,'(A)') "#"
      WRITE(6,'(A5, 2X,A7, 5X,A8, 4X,A6, 6X,A6, 6X,A9, 3X,A7, 4X,A8, 5X,A9, 3X,A7, 5X,A7)')    &
        "#  iu", "ep[keV]", "dedx_tot", "dedx_e", "dedx_i", "dedxc_tot", "dedxc_e", "dedxc_i", &
        "dedxq_tot", "dedxq_e", "dedxq_i"

     CALL dedx_bps(nni, ep, zp, mp, betab, zb, mb, nb,          &
        dedx_tot, dedx_i, dedx_e, dedxc_tot, dedxc_i, dedxc_e, & 
        dedxq_tot, dedxq_i, dedxq_e)
      WRITE(6,'(I5,10E12.4)') iu,ep,dedx_tot,dedx_e,dedx_i,dedxc_tot,dedxc_e,dedxc_i,&
        dedxq_tot,dedxq_e,dedxq_i 

      CALL dedx_fr_acoeff_bps(nni,ep,zp,mp,betab,zb,mb,nb,  &
        dedx_a_tot, dedx_a_i, dedx_a_e, dedxc_a_tot, dedxc_a_i, dedxc_a_e, & 
        dedxq_a_tot, dedxq_a_i, dedxq_a_e, dedxc_a_s_i, dedxc_a_s_e,       &
        dedxc_a_r_i, dedxc_a_r_e)
      WRITE(6,'(I5,10E12.4)') iu,ep,dedx_a_tot,dedx_a_e,dedx_a_i,dedxc_a_tot,dedxc_a_e,dedxc_a_i,&
        dedxq_a_tot,dedxq_a_e,dedxq_a_i 

      diff_dedx_tot =100*ABS((dedx_tot-dedx_a_tot)/dedx_a_tot)
      diff_dedx_e   =100*ABS((dedx_e-dedx_a_e)/dedx_a_e)
      diff_dedx_i   =100*ABS((dedx_i-dedx_a_i)/dedx_a_i)
      diff_dedxc_tot=100*ABS((dedxc_tot-dedxc_a_tot)/dedxc_a_tot)
      diff_dedxc_e  =100*ABS((dedxc_e-dedxc_a_e)/dedxc_a_e)
      diff_dedxc_i  =100*ABS((dedxc_i-dedxc_a_i)/dedxc_a_i)
      diff_dedxq_tot=100*ABS((dedxq_tot-dedxq_a_tot)/dedxq_a_tot)
      diff_dedxq_e  =100*ABS((dedxq_e-dedxq_a_e)/dedxq_a_e)
      diff_dedxq_i  =100*ABS((dedxq_i-dedxq_a_i)/dedxq_a_i)
      WRITE(6,'(I5,10E12.4)') iu,ep,diff_dedx_tot,diff_dedx_e,&
        diff_dedx_i,diff_dedxc_tot,diff_dedxc_e,diff_dedxc_i,&
        diff_dedxq_tot,diff_dedxq_e,dedxq_i 
STOP
!
! percent difference: 100*ABS[(new - old)/old ]
!

        diff_dedx_tot =ABS(dedx_tot-dedx_a_tot)/dedx_tot
        diff_dedx_i   =ABS(dedx_i-dedx_a_i)/dedx_i
        diff_dedx_e   =ABS(dedx_e-dedx_a_e)/dedx_e
        diff_dedxc_tot=ABS(dedxc_tot-dedxc_a_tot)/dedxc_tot
        diff_dedxc_i  =ABS(dedxc_i-dedxc_a_i)/dedxc_i
        diff_dedxc_e  =ABS(dedxc_e-dedxc_a_e)/dedxc_e
        diff_dedxq_tot=ABS(dedxq_tot-dedxq_a_tot)/dedxq_tot
        diff_dedxq_i  =ABS(dedxq_i-dedxq_a_i)/dedxq_i
        diff_dedxq_e  =ABS(dedxq_e-dedxq_a_e)/dedxq_e
!
!
      WRITE(6, '(I5,10E12.4)') iu,ep,diff_dedx_tot,diff_dedx_i,diff_dedx_e,  &
        diff_dedxc_tot,diff_dedxc_i,diff_dedxc_e,diff_dedxq_tot,diff_dedxq_i,&
        diff_dedxq_e 
        DO iu=0,nep
           ep=ep0 + dep*iu
           CALL dedx_bps(nni,ep,zp,mp,betab,zb,mb,nb,  &
               dedx_tot, dedx_i, dedx_e, dedxc_tot, dedxc_i, dedxc_e, & 
               dedxq_tot, dedxq_i, dedxq_e)
!
           CALL dedx_fr_acoeff_bps(nni,ep,zp,mp,betab,zb,mb,nb,  &
               dedx_a_tot, dedx_a_i, dedx_a_e, dedxc_a_tot, dedxc_a_i, dedxc_a_e, & 
               dedxq_a_tot, dedxq_a_i, dedxq_a_e)
!
           WRITE(6, '(I5,10E12.4)') iu, ep, dedx_tot, dedx_i, dedx_e, dedxc_tot, &
             dedxc_i, dedxc_e, dedxq_tot, dedxq_i, dedxq_e
           WRITE(6, '(I5,10E12.4)') iu, ep, dedx_a_tot, dedx_a_i, dedx_a_e, dedxc_a_tot, & 
               dedxc_a_i, dedxc_a_e, dedxq_a_tot, dedxq_a_i, dedxq_a_e
           WRITE(1, '(I5,10E12.4)') iu, ep, dedx_tot, dedx_i, dedx_e, dedxc_tot, &
             dedxc_i, dedxc_e, dedxq_tot, dedxq_i, dedxq_e
           WRITE(2, '(I5,10E12.4)') iu, ep, dedx_a_tot, dedx_a_i, dedx_a_e, dedxc_a_tot, & 
               dedxc_a_i, dedxc_a_e, dedxq_a_tot, dedxq_a_i, dedxq_a_e


           diff_dedx_tot =ABS(dedx_tot-dedx_a_tot)/dedx_tot
           WRITE(3, '(I5,10E12.4)') iu, ep, dedx_tot, dedx_a_tot, diff_dedx_tot
          


        ENDDO     

        CALL close_output
        CALL close_plasma
      END PROGRAM main
!
! Returns the plasma species arrays: betab, mb, nb, zb
! Allocates other plasma arrays
!
    SUBROUTINE define_plasma(ep, zp, mp, te, ti, ne, nni)
    USE allocatablevars 
    USE physvars
      IMPLICIT NONE
      REAL                                :: ep     ! projectile energy [keV]
      REAL                                :: zp     ! projectile charge
      REAL                                :: mp     ! projectile mass [keV]
      REAL                                :: te     ! plasma electtron temp [keV]
      REAL                                :: ti     ! plasma ion temp [keV]
      REAL                                :: ne     ! plasma electron no. density [cm^-3]
      INTEGER                             :: nni    ! number of ion species

!     REAL,    DIMENSION(:), ALLOCATABLE  :: betab  ! [1/keV]
!     REAL,    DIMENSION(:), ALLOCATABLE  :: mb     ! [keV]
!     REAL,    DIMENSION(:), ALLOCATABLE  :: nb     ! [cm^-3]
!     REAL,    DIMENSION(:), ALLOCATABLE  :: zb     ! [e]
!     REAL,    DIMENSION(:), ALLOCATABLE  :: gb, etab, mpb, mbpb

!
! projectile
!
      ep=1.      ! keV
      mp=4*MPKEV ! keV
      zp=2.      !
!
! plasma
!
      te=10.   ! keV
      ti=10.   ! keV
      ne=1.E25 ! cm^-3

      nni=2  ! number of ion species
      ALLOCATE(betab(1:nni+1),zb(1:nni+1),mb(nni+1),nb(1:nni+1))   ! allocatablevars
      ALLOCATE(gb(1:nni+1),etab(1:nni+1),mpb(nni+1),mbpb(1:nni+1)) ! allocatablevars

      zb(1)=-1.     ! Species charges
      zb(2)=+1.     ! 
      zb(3)=+1.     ! 
      mb(1)=MEKEV   ! Species masses [keV]
      mb(2)=2*MPKEV !
      mb(3)=3*MPKEV !
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

    SUBROUTINE write_output(nep, ep, zp, mp, te, ti, ne, nni, betab, zb, mb, nb)
    USE physvars
      IMPLICIT NONE
      INTEGER,                     INTENT(IN) :: nep    ! 
      REAL,                        INTENT(IN) :: ep     ! [keV]
      REAL,                        INTENT(IN) :: zp     ! 
      REAL,                        INTENT(IN) :: mp     ! [keV]
      REAL,                        INTENT(IN) :: te     ! [keV]
      REAL,                        INTENT(IN) :: ti     ! [keV]
      REAL,                        INTENT(IN) :: ne     ! [cm^-3]
      INTEGER,                     INTENT(IN) :: nni    ! number of ion species
      REAL,    DIMENSION(1:nni+1), INTENT(IN) :: betab  ! [1/keV]
      REAL,    DIMENSION(1:nni+1), INTENT(IN) :: zb     ! [e]
      REAL,    DIMENSION(1:nni+1), INTENT(IN) :: mb     ! [keV]
      REAL,    DIMENSION(1:nni+1), INTENT(IN) :: nb     ! [cm^-3]

      REAL  :: etae, ge, gi, ze, vp
      REAL,    DIMENSION(1:nni+1) :: gb, ab, etab
!
! open output files
!
      OPEN(UNIT=1, FILE="main.dedx.dat")           ! 
      OPEN(UNIT=2, FILE="main.dedx.fr.acoeff.dat") ! 
      OPEN(UNIT=3, FILE="main.dedx.diff.dat") ! 
!
! write header
!
!     
! write header
!     
      WRITE(6,'(A)') '#'  
      WRITE(6,'(A)') '#'
      WRITE(6,'(A, 3X,A21)') '#','Projectile Parameters'
      WRITE(6,'(A, 4X,A13, X,D12.4)') '#','charge      :', zp
      WRITE(6,'(A, 4X,A13, X,D12.4)') '#','mass   [amu]:', mp/AMUKEV
      WRITE(6,'(A, 4X,A13, X,D12.4)') '#','mass   [keV]:', mp
      WRITE(6,'(A, 4X,A13, X,D12.4)') '#','energy [keV]:', ep
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

      vp  = CC*SQRT(2*ep/mp)       ! [cm/s]
      ab  =0.5*betab*mb*vp*vp/CC2  ! [dimensionless] (1/2) betab(ib)*mb(ib)*vp2/CC2
      etab=ABS(zp*zb)*2.1870E8/vp  ! [dimensionless]

      WRITE(6,'(A7, 3D12.4)') '# ab  =', ab
      WRITE(6,'(A7, 3D12.4)') '# etab=', etab
      WRITE(1,'(A7, 3D12.4)') '# ab  =', ab
      WRITE(1,'(A7, 3D12.4)') '# etab=', etab
      WRITE(2,'(A7, 3D12.4)') '# ab  =', ab
      WRITE(2,'(A7, 3D12.4)') '# etab=', etab
!      
      WRITE(6,'(D15.5)') ne
      WRITE(6,'(D15.5)') te
      WRITE(6,'(D15.5)') ti
      WRITE(6,'(I10)')   nep
      WRITE(1,'(D15.5)') ne
      WRITE(1,'(D15.5)') te
      WRITE(1,'(D15.5)') ti
      WRITE(1,'(I10)')   nep
      WRITE(2,'(D15.5)') ne
      WRITE(2,'(D15.5)') te
      WRITE(2,'(D15.5)') ti
      WRITE(2,'(I10)')   nep

      WRITE(6,'(A5, 2X,A7, 5X,A5, 7X,A3, 9X,A3, 9X,A6, 6X,A4, 8X,A4, 8X,A6, 6X,A4, 8X,A4)') &
        "#  iu", "ep[keV]", "a_tot", "a_i", "a_e", "ac_tot", "ac_i_lim", "ac_e", &
        "aq_tot", "aq_i", "aq_e"
      WRITE(1,'(A)') "#"
      WRITE(1,'(A5, 2X,A7, 5X,A5, 7X,A3, 9X,A3, 9X,A6, 6X,A4, 8X,A4, 8X,A6, 6X,A4, 8X,A4)') &
        "#  iu", "ep[keV]", "a_tot", "a_i", "a_e", "ac_tot", "ac_i", "ac_e", "aq_tot", "aq_i","aq_e"
      WRITE(2,'(A)') "#"
      WRITE(2,'(A5, 2X,A7, 5X,A9, 3X,A7, 5X,A7, 5X,A10, 2X,A8, 4X,A8, 4X,A10, 2X,A8, 4X,A8)') &
        "#  iu", "ep[keV]", "a_tot_lim", "a_i_lim", "a_e_lim", "ac_tot_lim", "ac_i_lim", "ac_e_lim", &
        "aq_tot_lim", "aq_i_lim", "aq_e_lim"
    END SUBROUTINE write_output

    SUBROUTINE close_output
    IMPLICIT NONE
      CLOSE(1)
      CLOSE(2)
    END SUBROUTINE close_output

!
! calculates dedx from the A-coefficients.
!

      SUBROUTINE dedx_fr_acoeff_bps(nni,ep,zp,mp,betab,zb,mb,nb,  &
            dedx_a_tot, dedx_a_i, dedx_a_e, dedxc_a_tot, dedxc_a_i, dedxc_a_e, & 
            dedxq_a_tot, dedxq_a_i, dedxq_a_e)

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

        REAL :: a_tot_p, a_i_p, a_e_p, ac_tot_p, ac_i_p, ac_e_p, aq_tot_p, aq_i_p, aq_e_p
        REAL :: a_tot_m, a_i_m, a_e_m, ac_tot_m, ac_i_m, ac_e_m, aq_tot_m, aq_i_m, aq_e_m
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

        epp=ep+dep
        CALL bps_acoeff_ei_mass(nni,epp,zp,mp,betab,zb,mb,nb,a_tot_p,a_i_p, &
          a_e_p,ac_tot_p,ac_i_p,ac_e_p,aq_tot_p,aq_i_p,aq_e_p,ac_s_i,ac_s_e,&
          ac_r_i,ac_r_e)

        epm=ep-dep
        CALL bps_acoeff_ei_mass(nni,epm,zp,mp,betab,zb,mb,nb,a_tot_m,a_i_m, &
          a_e_m,ac_tot_m,ac_i_m,ac_e_m,aq_tot_m,aq_i_m,aq_e_m)

        CALL bps_acoeff_ei_mass(nni,ep, zp,mp,betab,zb,mb,nb,a_tot,a_i,& 
          a_e,ac_tot,ac_i,ac_e,aq_tot,aq_i,aq_e)

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
      END SUBROUTINE dedx_fr_acoeff_bps
