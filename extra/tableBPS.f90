
    MODULE bpsvars
!
! misc parameters
!
      INTEGER, PARAMETER :: R=3        ! thermal velocity parameter
      INTEGER, PARAMETER :: I=1        ! plasma species index
      REAL               :: K          ! arbitrary wave number units a0^-1 

! plasma parameters: values set in dedx_bps
!
      REAL,    DIMENSION(:), ALLOCATABLE :: kb2, ab, bb, cb, eb, fb, rb, ggb
      REAL,    DIMENSION(:), ALLOCATABLE :: ab2, etb, rmb0, rrb0, mb0, mpb0
      LOGICAL, DIMENSION(:), ALLOCATABLE :: lzb
      REAL    :: cp1, cp2, cp3, vth, vthc, mp0, kd
      INTEGER :: NNB  ! number of plasma species = ni+1
    END MODULE bpsvars

      MODULE allocatablevars
      INTEGER :: ni  ! number of ion species
!      REAL,    DIMENSION(:), ALLOCATABLE :: betab, zb, mb, nb
!      REAL,    DIMENSION(:), ALLOCATABLE :: gb, etab, mpb, mbpb
      REAL,    DIMENSION(1:3) :: betab, zb, mb, nb
      REAL,    DIMENSION(1:3) :: gb, etab, mpb, mbpb
      END MODULE allocatablevars

    MODULE controlvars
      LOGICAL :: asymptotic_forms=.TRUE.
      LOGICAL :: small_E_form=.TRUE.
      LOGICAL :: large_E_form=.TRUE.
    END MODULE controlvars

      MODULE mathvars
!
! mathematical constants to 20 significant figures
!
      REAL,    PARAMETER :: PI   =3.1415926535897932385 ! pi
      REAL,    PARAMETER :: TWOPI=6.2831853071795864769 ! 2pi
      REAL,    PARAMETER :: FPI3 =2.3561944901923449288 ! 4*pi/3
      REAL,    PARAMETER :: SQPI =1.7724538509055160273 ! sqrt(pi)
      REAL,    PARAMETER :: TSPI =1.1283791670955125739 ! 2/sqrt(pi)
      REAL,    PARAMETER :: GAMMA=0.57721566490153286061! Euler Gamma
      REAL,    PARAMETER :: ZETA3=1.2020569031595942854 ! zeta(3)
      REAL,    PARAMETER :: EXP2E=3.1722189581254505277 ! exp(2*GAMMA)
      REAL,    PARAMETER :: LOG2 =0.69314718055994530942! ln(2)
      REAL,    PARAMETER :: LOG4 =1.3862943611198906188 ! ln(4)
      REAL,    PARAMETER :: LOG8 =2.0794415416798359283 ! ln(8)
      REAL,    PARAMETER :: LOG16=2.7725887222397812377 ! ln(16)
      REAL,    PARAMETER :: THIRD=1./3.                 ! 1/3=0.333333....

      END MODULE mathvars
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

!===================================================================
      program tableBPS
      implicit none
      real*8  :: dedze, dedzi, de, E0, f, DT, ne
      real*8  :: rho(0:72), Te(0:54), EeV(0:100)
      integer :: ntabrho, ntabte, ntabene, ir, it, ie, id, iu
      real*8, dimension(0:72,0:54,0:100) :: tabi, tabe

      ntabrho = 63
      do ir = 0, ntabrho
        id = int( float(ir) / 9. )
        iu = ir - 9 * id
        rho(ir) = 0.0001 * (iu+1)*10.**id	! DT densities [g/cc]: 0.0001 - 1000
      enddo

      ntabte  = 45
      do it = 0, ntabte
        id = int( float(it) / 9. )
        iu = it - 9 * id
        Te(it) = (iu+1)*10.**id			! Temperatures Te = Ti [eV]: 1 - 100 keV
      enddo

      ntabene = 100
      E0 = 17.6e6 / 5.				! eV
      f = 10.**(-2./(ntabene-1.))
      de  = E0 * ( f - 1. ) / ( f**ntabene - 1. )
      EeV(0) = E0
      do ie = 1, ntabene - 1
        EeV(ie) = EeV(ie-1) - de * f**(ie-1)	! alpha energies [eV]: 3.52e6 - 0
      enddo
      EeV(ntabene) = 0.

!-----

      tabe = 0.
      tabi = 0.

      DT = 0.5 * ( 2.013553212745 + 3.0160492 ) * 1.66053904e-24	! DT mass [g]
      do ir = 0, ntabrho
        do it = 0, ntabte
          do ie = 0, ntabene - 1

            ne = rho(ir) / DT

            call dedx(ne,Te(it),EeV(ie), dedze,dedzi ) ! for DT plasma A=2.5, Z=1, Ti=Te

            tabe(ir,it,ie) = dedze	! MeV/um
            tabi(ir,it,ie) = dedzi

          enddo
          write(6,*) ir,'/',ntabrho, it,'/',ntabte
        enddo
      enddo

!-----
! the stopping power is forced to be positive when dedx(ie) < 0; dedx is lineraized
! between dedx=0 for EeV=0 and the latest posivive value dedx(ie-1) for EeV(ie-1)

      do ir = 0, ntabrho
        do it = 0, ntabte
          do ie = 0, ntabene

            if ( tabe(ir,it,ie) < 0. ) then
              tabe(ir,it,ie) = EeV(ie) * tabe(ir,it,ie-1) / EeV(ie-1)
            endif

            if ( tabi(ir,it,ie) < 0. ) then
              tabi(ir,it,ie) = EeV(ie) * tabi(ir,it,ie-1) / EeV(ie-1)
            endif

          enddo
        enddo
      enddo

!-----

      open(1,err=1,file='tableBPS',status='old',form='formatted')
      close(1,status='delete')
1     open(1,file='tableBPS',status='new',form='formatted')
      write(1,*) ntabrho, ntabte, ntabene
      write(1,*) ( rho(ir), ir = 0, ntabrho)
      write(1,*) (  Te(it), it = 0, ntabte )
      write(1,*) ( EeV(ie), ie = 0, ntabene)
      do ir = 0, ntabrho
        do it = 0, ntabte
          do ie = 0, ntabene
            write(1,*) - tabe(ir,it,ie)*1.e10, - tabi(ir,it,ie)*1.e10		! [eV/cm]
          enddo
        enddo
      enddo

      close(1)

      end program tableBPS

!-----------------------------------------------------------------------------------

      Subroutine dedx(nel,TeV,EeV, dedx_e,dedx_i )	! copied from main.dedx.f90
!
! This program finds the range of a charged partcile in a plasma, where the
! charged particle is created with threshold energy Ep. 
!
! The subroutine rk4(y,t,dt,nit) advances y at time t by nit iterations, each
! of which has a constant time-step dt. 
!
! nts = number of time steps
! nit = number of iterations per time step
!
      USE allocatablevars
      USE bpsvars
      USE physvars
      USE mathvars
        IMPLICIT NONE
        REAL    :: dedx_tot,  dedx_i,  dedx_e
        REAL    :: dedxc_tot, dedxc_i, dedxc_e
        REAL    :: dedxq_tot, dedxq_i, dedxq_e
        REAL    :: dz, path, TeV, EeV, nel

        REAL    :: te, ti, ne, ep, mp, zp, de, epp
        INTEGER :: nni

        te = TeV / 1000.	! Electron temperature [keV]
        ep = EeV / 1000.	! Projectile energy    [keV]
	ne = nel
        ti=te
        mp=4*MPKEV  ! Projectile mass    [keV]
        zp=2.       ! Projectile charge  [e]

!
! DT plasma with alpha particle projectile
!
        CALL define_plasma(te,ti,ne,nni)

           IF (ep <= 0. ) ep=1.E-5

           CALL dedx_bps(nni, ep, zp, mp, betab, zb, mb, nb,   &
             dedx_tot, dedx_i, dedx_e, dedxc_tot, dedxc_i, & 
             dedxc_e, dedxq_tot, dedxq_i, dedxq_e) ! [MeV/micron] with epp/1000. in MeV

        END subroutine dedx

!
! SUBROUTINE define_plasma:
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
!     REAL,    DIMENSION(:), ALLOCATABLE  :: betab  ! [1/keV]
!     REAL,    DIMENSION(:), ALLOCATABLE  :: mb     ! [keV]
!     REAL,    DIMENSION(:), ALLOCATABLE  :: nb     ! [cm^-3]
!     REAL,    DIMENSION(:), ALLOCATABLE  :: zb     ! [e]
!     REAL,    DIMENSION(:), ALLOCATABLE  :: gb, etab, mpb, mbpb

      nni=2  ! number of ion species

!mt      ALLOCATE(betab(1:nni+1),zb(1:nni+1),mb(nni+1),nb(1:nni+1))   ! allocatablevars
!mt      ALLOCATE(gb(1:nni+1),etab(1:nni+1),mpb(nni+1),mbpb(1:nni+1)) ! allocatablevars

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

      vp  = CC*SQRT(2*ep/mp)       ! [cm/s]
      ab  =0.5*betab*mb*vp*vp/CC2  ! [dimensionless] (1/2) betab(ib)*mb(ib)*vp2/CC2
      etab=ABS(zp*zb)*2.1870E8/vp  ! [dimensionless]

      WRITE(6,'(A7, 3D12.4)') '# ab  =', ab
      WRITE(6,'(A7, 3D12.4)') '# etab=', etab
      WRITE(1,'(A7, 3D12.4)') '# ab  =', ab
      WRITE(1,'(A7, 3D12.4)') '# etab=', etab

!      WRITE(6,'(A)') "#"
!      WRITE(6,'(A5, 2X,A7, 5X,A5, 7X,A3, 9X,A3, 9X,A6, 6X,A4, 8X,A4, 8X,A6, 6X,A4, 8X,A4)') &
!        "#  iu", "ep[keV]", "a_tot", "a_e", "a_i", "ac_tot", "ac_e_lim", "ac_i", &
!        "aq_tot", "aq_e", "aq_i"
!      WRITE(1,'(A)') "#"
!      WRITE(1,'(A5, 2X,A7, 5X,A5, 7X,A3, 9X,A3, 9X,A6, 6X,A4, 8X,A4, 8X,A6, 6X,A4, 8X,A4)') &
!        "#  iu", "ep[keV]", "a_tot", "a_e", "a_i", "ac_tot", "ac_e", "ac_i", "aq_tot", "aq_e","aq_i"

    END SUBROUTINE write_output
      FUNCTION ferf(x)
      IMPLICIT NONE
      REAL*8 x, ferf
      ferf=ERF(x)
      RETURN
      END
! 
! 1. Routines common to dedx_bps and acoeff_bps.
! 2. Routines useful for general purposes.
!
!======================================
! acoeff and dedx_bps
!
      SUBROUTINE x_collect(ib, ibmax, xc_s, xc_r, xq, &
        x_tot, x_i, x_e, xc_tot, xc_i, xc_e, xq_tot,  &
        xq_i, xq_e, xc_s_i, xc_s_e, xc_r_i, xc_r_e)
        IMPLICIT NONE
        INTEGER, INTENT(IN)    :: ib     ! species index
        INTEGER, INTENT(IN)    :: ibmax  ! species index maximum = NNB+1
        REAL,    INTENT(IN)    :: xc_s   ! singular contribution
        REAL,    INTENT(IN)    :: xc_r   ! regular  contribution
        REAL,    INTENT(IN)    :: xq     ! quantum  contribution
        REAL,    INTENT(INOUT) :: x_tot  ! running total over ions
        REAL,    INTENT(INOUT) :: x_i    ! running total over ions
        REAL,    INTENT(INOUT) :: x_e    ! electron component
        REAL,    INTENT(INOUT) :: xc_tot ! running total over ions
        REAL,    INTENT(INOUT) :: xc_i   ! running total over ions
        REAL,    INTENT(INOUT) :: xc_e   ! electron component
        REAL,    INTENT(INOUT) :: xq_tot ! running total over ions
        REAL,    INTENT(INOUT) :: xq_i   ! running total over ions
        REAL,    INTENT(INOUT) :: xq_e   ! electron component
        REAL,    INTENT(INOUT) :: xc_s_i ! singular ion
        REAL,    INTENT(INOUT) :: xc_s_e ! singular electron
        REAL,    INTENT(INOUT) :: xc_r_i ! regular  ion
        REAL,    INTENT(INOUT) :: xc_r_e ! regular  electron
        REAL :: xc_sr
        xc_sr=xc_s + xc_r
        IF (ib==1) THEN
           xc_e=xc_sr
           xq_e=xq
           x_e =xc_e + xq_e

           xc_s_e=xc_s
           xc_r_e=xc_r
        ELSE
           xc_i=xc_i + xc_sr
           xq_i=xq_i + xq
           x_i =x_i  + xc_sr + xq

           xc_s_i=xc_s_i + xc_s
           xc_r_i=xc_r_i + xc_r
        ENDIF
        IF (ib==ibmax) THEN
           xc_tot = xc_e + xc_i
           xq_tot = xq_e + xq_i
           x_tot  = xc_tot  + xq_tot
        ENDIF
      END SUBROUTINE x_collect

!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! to do: compare this with fri() in dedx.f90
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
      SUBROUTINE frfi(x,nni,kb2,ab,fr,fi,fabs,farg) 
      USE mathvars
        IMPLICIT NONE
        REAL,                        INTENT(IN)  :: x      ! x^2=(1/2)\beta_b m_b v_p^2
        INTEGER,                     INTENT(IN)  :: nni    ! Number of ion species
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: kb2    ! alpha(b)
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: ab     ! a(b)
        REAL,                        INTENT(OUT) :: fr     ! Dimensionless: multiply fr and
        REAL,                        INTENT(OUT) :: fi     ! fi by ke^2 for physical units
        REAL,                        INTENT(OUT) :: fabs   !
        REAL,                        INTENT(OUT) :: farg   !

        REAL, PARAMETER     :: XMIN=1.E-4
        REAL                :: xib, daw
        INTEGER             :: ib
        fr=0
        fi=0
        DO ib=1,nni+1
          xib=ab(ib)*x
          fr=fr + kb2(ib)*(1-2*xib*daw(xib))
          fi=fi + kb2(ib)*xib*EXP(-xib*xib)
        ENDDO
        fi=fi*SQRT(PI)
        fabs=SQRT(fr*fr + fi*fi)
        farg=ATAN2(fi,fr) ! branch cut: x=0 y<0
      END SUBROUTINE frfi

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
!

      SUBROUTINE a_collect(ib, ibmax, ac_s, ac_r, aq, a_tot, a_i, a_e, &
        ac_tot, ac_i, ac_e, aq_tot, aq_i, aq_e, ac_s_i, ac_s_e, ac_r_i, ac_r_e)
        IMPLICIT NONE
        INTEGER, INTENT(IN)    :: ib     ! species index
        INTEGER, INTENT(IN)    :: ibmax  ! species index maximum = NNB+1
        REAL,    INTENT(IN)    :: ac_s   ! singular contribution
        REAL,    INTENT(IN)    :: ac_r   ! regular  contribution
        REAL,    INTENT(IN)    :: aq     ! quantum  contribution
                                         !
        REAL,    INTENT(INOUT) :: a_tot  ! running total over ions
        REAL,    INTENT(INOUT) :: a_i    ! running total over ions
        REAL,    INTENT(INOUT) :: a_e    ! electron component
        REAL,    INTENT(INOUT) :: ac_tot ! running total over ions
        REAL,    INTENT(INOUT) :: ac_i   ! running total over ions
        REAL,    INTENT(INOUT) :: ac_e   ! electron component
        REAL,    INTENT(INOUT) :: aq_tot ! running total over ions
        REAL,    INTENT(INOUT) :: aq_i   ! running total over ions
        REAL,    INTENT(INOUT) :: aq_e   ! electron component
        REAL,    INTENT(INOUT) :: ac_s_i
        REAL,    INTENT(INOUT) :: ac_s_e
        REAL,    INTENT(INOUT) :: ac_r_i
        REAL,    INTENT(INOUT) :: ac_r_e
        REAL :: ac_sr
        ac_sr=ac_s + ac_r
        IF (ib==1) THEN
           ac_e=ac_sr
           aq_e=aq
           a_e =ac_e + aq_e

           ac_s_e=ac_s
           ac_r_e=ac_r
        ELSE
           ac_i=ac_i + ac_sr
           aq_i=aq_i + aq
           a_i =a_i  + ac_sr + aq

           ac_s_i=ac_s_i + ac_s
           ac_r_i=ac_r_i + ac_r
        ENDIF
        IF (ib==ibmax) THEN
           ac_tot = ac_e + ac_i
           aq_tot = aq_e + aq_i
           a_tot  = ac_tot  + aq_tot
        ENDIF
      END SUBROUTINE a_collect
!
! This is a fit to repsi(x)-ln(x) = Re Psi(1 + I*x) - ln(x), 
! where Psi is the PolyGamma function. The accuracy is 0.1%.
!
      FUNCTION repsilog(x)  
      USE bpsvars
      USE mathvars
        IMPLICIT NONE
        REAL, INTENT(IN) :: x
        REAL, PARAMETER :: XMIN=0.16E0, XMAX=1.5E0
        REAL, PARAMETER :: TZETA3=2.404113806319188 ! 2*ZETA(3)
        REAL, PARAMETER :: A=0.1E0, B=1.33333E0, C=1.125E0
        REAL :: repsilog
        REAL :: x2, x4
        IF (x .LE. XMIN) THEN
           x2=x**2
           repsilog=-GAMMA + ZETA3*x2 - LOG(x)
        ELSEIF (x .GE. XMAX) THEN
           x2=x**2
           x4=x2*x2
           repsilog=1./(12.*x2)+1./(120.*x4)
        ELSE
           x2=x*x
           repsilog=0.5E0*LOG(1 + (EXP2E*x2*x2 + TZETA3*x2)/(1+x2))
           repsilog=repsilog/(1 - A*EXP(-B*x - C/x)) - GAMMA - LOG(x)
        ENDIF
      END FUNCTION repsilog
!
! This is a fit to repsi(x) = Re Psi(1 + I*x), where
! Psi is the PolyGamma function. The accuracy is 0.1%.
!
      FUNCTION repsi(x)  
      USE bpsvars
      USE mathvars
        IMPLICIT NONE
        REAL, INTENT(IN) :: x
        REAL, PARAMETER :: XMIN=0.16E0, XMAX=1.5E0
        REAL, PARAMETER :: TZETA3=2.404113806319188 ! 2*ZETA(3)
        REAL, PARAMETER :: A=0.1E0, B=1.33333E0, C=1.125E0
        REAL :: repsi
        REAL :: x2, x4
        IF (x .LE. XMIN) THEN
           x2=x**2
           repsi=-GAMMA + ZETA3*x2
        ELSEIF (x .GE. XMAX) THEN
           x2=x**2
           x4=x2*x2
           repsi=LOG(x)+1./(12.*x2)+1./(120.*x4)
        ELSE
           x2=x*x
           repsi=0.5E0*LOG(1 + (EXP2E*x2*x2 + TZETA3*x2)/(1+x2))
           repsi=repsi/(1 - A*EXP(-B*x - C/x)) - GAMMA
        ENDIF
      END FUNCTION repsi
!
!              d
! repsi1(x) = --- Re[ Psi(1 + I*x) = -Im Psi'(1 + I*x). 
!             dx
!
      FUNCTION repsi1(x)  
        IMPLICIT NONE
        REAL, INTENT(IN) :: x
        REAL :: repsi1

        REAL, PARAMETER :: XMIN=0.14E0, X1=0.7E0 ,XMAX=1.9E0
        REAL, PARAMETER :: ZETA32=2.404113806319188 ! 2*ZETA(3)
        REAL, PARAMETER :: ZETA54=4.147711020573480 ! 4*ZETA(5)
        REAL, PARAMETER :: a0= 0.004211521868683916
        REAL, PARAMETER :: a1= 2.314767988469241000
        REAL, PARAMETER :: a2= 0.761843932767193200
        REAL, PARAMETER :: a3=-7.498711815965575000
        REAL, PARAMETER :: a4= 7.940030433629257000
        REAL, PARAMETER :: a5=-2.749533936429732000
        REAL, PARAMETER :: b0=-0.253862873373708200
        REAL, PARAMETER :: b1= 4.600929855835432000
        REAL, PARAMETER :: b2=-6.761540444078382000
        REAL, PARAMETER :: b3= 4.467238548899841000
        REAL, PARAMETER :: b4=-1.444390097613873500
        REAL, PARAMETER :: b5= 0.185954029179227070
        REAL :: xi
        IF ( x .LE. XMIN) THEN             ! x < xmin=0.14 
           repsi1=ZETA32*x - ZETA54*x*x*x  ! accurate to 0.1%
        ELSEIF (x .LE. x1) THEN
           repsi1=a5                       ! xmin < x < x1=0.7
           repsi1=a4 + repsi1*x            ! accurate to 0.002%
           repsi1=a3 + repsi1*x            ! a0 + a1*x + a2*x^2 +
           repsi1=a2 + repsi1*x            ! a3* x^3 + a4*x^4 + a5*x^5
           repsi1=a1 + repsi1*x
           repsi1=a0 + repsi1*x
        ELSEIF (x .LE. xmax) THEN          ! x1 < x < xmax=1.9
           repsi1=b5                       ! accurate to 0.1%
           repsi1=b4+repsi1*x              ! b0 + b1*x + b2*x^2 + 
           repsi1=b3+repsi1*x              ! b3*x^3 + b4*x^4 +
           repsi1=b2+repsi1*x              ! b5*x^5
           repsi1=b1+repsi1*x              
           repsi1=b0+repsi1*x            
        ELSE
           xi=1/x                          ! x > xmax=1.9
           repsi1=-1.E0/30.E0              ! accurate to 0.08%
           repsi1=repsi1*xi                ! 1/x - 1/6x^3 - 1/30x^5
           repsi1=-1.E0/6.D0 + repsi1*xi   ! 
           repsi1=repsi1*xi                !
           repsi1=1.E0 + repsi1*xi         !
           repsi1=repsi1*xi                !
        ENDIF
      END FUNCTION repsi1
!
!              d^2
! repsi2(x) = ---- Re[ Psi(1 + I*x) = -Re Psi''(1 + I*x). 
!             dx^2
!
      FUNCTION repsi2(x) 
        IMPLICIT NONE
        REAL, INTENT(IN) :: x
        REAL :: repsi2
        REAL, PARAMETER :: XMIN=0.18E0, X1=1.2E0 ,XMAX=2.5E0
        REAL, PARAMETER :: ZETA32=2.4041138063191880 ! 2*ZETA(3)
        REAL, PARAMETER :: ZETA512=12.44313306172044 ! 12*ZETA(5)
        REAL, PARAMETER :: ZETA730=30.250478321457685! 30*ZETA(7)
        REAL, PARAMETER :: a0= 2.42013533575662130
        REAL, PARAMETER :: a1=-0.41115258967949336
        REAL, PARAMETER :: a2=-8.09116694062588400
        REAL, PARAMETER :: a3=-24.9364824558827640
        REAL, PARAMETER :: a4=114.8109056152471800
        REAL, PARAMETER :: a5=-170.854545232781960
        REAL, PARAMETER :: a6=128.8402466765824700
        REAL, PARAMETER :: a7=-50.2459090010302060 
        REAL, PARAMETER :: a8= 8.09941032385266400
        REAL, PARAMETER :: b0= 4.98436272402513600
        REAL, PARAMETER :: b1=-16.6546466531059530
        REAL, PARAMETER :: b2= 20.6941300362041100
        REAL, PARAMETER :: b3=-13.3726837850936920 
        REAL, PARAMETER :: b4= 4.83094787289278800 
        REAL, PARAMETER :: b5=-0.92976482601030100
        REAL, PARAMETER :: b6= 0.07456475055097825
        REAL :: xi, xx
        IF ( x .LE. XMIN) THEN
           xx=x*x
           repsi2=ZETA32 - ZETA512*xx + ZETA730*xx*xx ! x < xmin=0.18
        ELSEIF (x .LE. x1) THEN                       ! accurate to 0.1%
           repsi2=a8                                  !
           repsi2=a7 + repsi2*x                       ! xmin < x < x1=1.2
           repsi2=a6 + repsi2*x                       ! accurate to 0.01% 
           repsi2=a5 + repsi2*x                       ! a0 + a1*x + a2*x^2 +
           repsi2=a4 + repsi2*x                       ! a3*x^3 + a4*x^4 +
           repsi2=a3 + repsi2*x                       ! a5*x^5 + a6*x^6 +
           repsi2=a2 + repsi2*x                       ! a7*x&7 + a8*x^8 
           repsi2=a1 + repsi2*x
           repsi2=a0 + repsi2*x 
        ELSEIF (x .LE. xmax) THEN                     ! x1 < x < xmax=2.5
           repsi2=b6                                  ! accurate to 0.2%
           repsi2=b5+repsi2*x                         ! b0 + b1*x + b2*x^2 +
           repsi2=b4+repsi2*x                         ! b3*x^3 + b4*x^4 +
           repsi2=b3+repsi2*x                         ! b5*x^5 + b6*x^6
           repsi2=b2+repsi2*x          
           repsi2=b1+repsi2*x          
           repsi2=b0+repsi2*x          
        ELSE
           xi=1/x                                     ! x > xmax=2.5
           xi=xi*xi                                   ! accurate to 0.07%
           repsi2= 1.E0/6.E0                          ! -1/x^2 + 1/2x^4 + 
           repsi2= 0.5E0 + repsi2*xi                  ! 1/6x^6
           repsi2=-1. + repsi2*xi                     !
           repsi2=repsi2*xi
        ENDIF
      END FUNCTION repsi2
!
!         /1
!         |              ln[1-u]
! j1(x) = | du e^{-x u} ---------
!         |              Sqrt[u]
!         /0
!
! Approximates the j1 integral by a rational function
!
      FUNCTION j1(x)
      USE bpsvars   
        IMPLICIT NONE
        REAL, INTENT(IN) :: x
        REAL             :: j1
        REAL,    PARAMETER :: SQPI =1.772453851    ! sqrt(pi)
        REAL,    PARAMETER :: LOG8 =2.079441542    ! ln(8)
        REAL,    PARAMETER :: LOG16=2.772588722    ! ln(16)
        INTEGER, PARAMETER         :: NNJ1=3, NMJ1=2*NNJ1-1 ! NMJ1=5
        REAL,    DIMENSION(0:NMJ1) :: J1B=(/ &
        -926.65E0  ,         & !b0
        787.016E0  ,         & !b1
        -329.764E0 ,         & !b2
        39.7406E0  ,         & !b3
        -0.173896E0,         & !b4
        -1.66913E0 /)          !b5
        REAL,    DIMENSION(0:NMJ1) :: J1A=(/ &
        -926.65E0,           & !a0=b0
        787.165E0  ,         & !a1
        -213.584E0 ,         & !a2
        -1.04219E0 ,         & !a3
        33.594E0   ,         & !a4
        -11.8391E0/)           !a5
        REAL, PARAMETER :: J1MM= (4./9.)*(4.-LOG8)    ! 0.8535815 
        REAL, PARAMETER :: J1BB=-(4.-LOG16)           !-1.2274113
        REAL, PARAMETER :: J1AA=-SQPI/2.              !-0.8862270
        REAL, PARAMETER :: J1CC= 0.1E0
        REAL, PARAMETER :: J1EE= 0.2E0
        REAL, PARAMETER :: J1GG=-3.*SQPI/8.
        REAL,    PARAMETER   :: XMIN=0.1D0, XMAX=20.D0
        REAL    :: x2, x4
        REAL    :: xx, ra, rc
        REAL    :: y, y3, y5
        INTEGER :: n
!
! analytic asymptotic forms
!
        IF (x .LE. XMIN) THEN
           j1=J1MM*x+J1BB
        ELSEIF (x .GT. XMAX) THEN
           y=SQRT(x)           ! x^1/2
           y3=x*y              ! x^3/2
           y5=y3*x             ! x^5/2
           j1=J1AA/y3 +J1GG/y5
        ELSE
!
! numerical asymptotic form
!
           x2=x*x
           x4=x**3.5
           j1=(J1MM*x+J1BB)/(J1CC*x4+1) + J1EE*J1AA*x2/(J1EE*x4+1)
!
! spline correction
! 
           ra=0.E0
           rc=0.E0
           xx=1.E0 
           DO n=0,NMJ1
              ra=ra+J1A(n)*xx
              rc=rc+J1B(n)*xx
              xx=x*xx
           ENDDO
           ra=ra+xx
           rc=rc+xx
           j1=j1*rc/ra
      ENDIF
 END FUNCTION j1
!
!         /1
!         |              
! j2(x) = | du e^{-x u} ln[1-u]*Sqrt[u]
!         |              
!         /0
!
! Approximates the j2 integral by a rational function
!
 FUNCTION j2(x)
   USE bpsvars
   IMPLICIT NONE
   REAL, INTENT(IN) :: x
   REAL             :: j2
   REAL,    PARAMETER :: SQPI =1.772453851    ! sqrt(pi)
   REAL,    PARAMETER :: GAMMA=0.577215665    ! Euler Gamma
   REAL,    PARAMETER :: LOG2 =0.6931471806   ! ln(2)
   REAL,    PARAMETER :: LOG8 =2.079441542    ! ln(8)
   INTEGER, PARAMETER         :: NNJ2=3, NMJ2=2*NNJ2-1 ! NMJ2=5
   REAL,    DIMENSION(0:NMJ2) :: J2B=(/ &
   87.1714E0 ,         & !b0
   -277.584E0,         & !b1
   329.082E0 ,         & !b2
   -180.982E0,         & !b3
   56.7202E0 ,         & !b4
   -8.60238E0/)          !b5
   REAL,    DIMENSION(0:NMJ2) :: J2A=(/ & 
   87.1714E0 ,         & !a0=b0
   -277.693E0,         & !a1
   329.801E0 ,         & !a2
   -184.219E0,         & !a3
   59.9325E0 ,         & !a4
   -10.1138E0/)          !a5
   REAL, PARAMETER ::  J2MM= 4.*(23.-15.*LOG2)/75. ! 0.6721489
   REAL, PARAMETER ::  J2BB=-4.*(4.-LOG8)/9.       !-0.8535815
   REAL, PARAMETER ::  J2AA=-3.*SQPI/4.            !-1.3293405 
   REAL, PARAMETER ::  J2CC= 0.5E0
   REAL, PARAMETER ::  J2EE= 0.2E0
   REAL, PARAMETER ::  J2GG=-15.*SQPI/16.
   REAL,    PARAMETER   :: XMIN=0.1, XMAX=30.D0
   REAL    :: x2, x4
   REAL    :: y, y5, y7
   REAL    :: xx, ra, rc
   INTEGER :: n
!
! analytic asymptotic forms
!
  IF (x .LE. XMIN) THEN
     j2=J2MM*x+J2BB
  ELSEIF (x .GT. XMAX) THEN
     y=SQRT(x)          ! x^1/2
     y5=x*x*y           ! x^5/2
     y7=y5*x            ! x^7/2
     j2=J2AA/y5 +J2GG/y7
  ELSE
!
! numerical asymptotic form
!
     x2=x*x
     x4=x**4.5
     j2=(J2MM*x+J2BB)/(J2CC*x4+1) + &
        J2EE*J2AA*x2/(J2EE*x4+1)
!
! spline correction
! 
     ra=0.E0
     rc=0.E0
     xx=1.E0
     DO n=0,NMJ2
        ra=ra+J2A(n)*xx
        rc=rc+J2B(n)*xx
        xx=x*xx
     ENDDO
     ra=ra+xx
     rc=rc+xx
     j2=j2*rc/ra
  ENDIF
 END FUNCTION j2
!
!         /1
!         |               ln[u]
! j3(x) = | du e^{-x u} --------
!         |              Sqrt[u]
!         /0
!
! Approximates the j3 integral by a rational function
!
 FUNCTION j3(x)
   USE bpsvars
   IMPLICIT NONE
   REAL, INTENT(IN) :: x
   REAL             :: j3
   REAL,    PARAMETER :: SQPI =1.772453851    ! sqrt(pi)
   REAL,    PARAMETER :: GAMMA=0.577215665    ! Euler Gamma
   REAL,    PARAMETER :: LOG4 =1.386294361    ! ln(4)
   REAL, PARAMETER ::  J3MM = 4./9.                ! 0.444444 
   REAL, PARAMETER ::  J3BB =-4.0E0                !-4.0       
   REAL, PARAMETER ::  J3AA1=-SQPI                 !-1.7724539 
   REAL, PARAMETER ::  J3AA2=-SQPI*(GAMMA+LOG4)    !-3.4802318 
   REAL, PARAMETER ::  J3CC = 0.1E0 
   REAL, PARAMETER ::  J3EE = 0.2E0 
   REAL, PARAMETER :: XMIN=0.4, XMAX=2.3D0
   REAL    :: y
!
! analytic asymptotic forms
!
  IF (x .LE. XMIN) THEN
     j3=J3MM*x+J3BB
  ELSEIF (x .GE. XMAX) THEN
     y=SQRT(x)
     j3=(J3AA1*LOG(x) + J3AA2)/y
  ELSE
     j3=-4.00005 - 0.0380953*SQRT(x) + 0.599082*x - 0.206644*x**1.5 +  0.0227552*x**2
  ENDIF
END FUNCTION j3
!
!         /1
!         |               
! j4(x) = | du e^{-x u} ln[u]*Sqrt[u]
!         |              
!         /0
!
! Approximates the j3 integral by a rational function
!
 FUNCTION j4(x)
   USE bpsvars
   IMPLICIT NONE
   REAL, INTENT(IN) :: x
   REAL             :: j4
   INTEGER, PARAMETER         :: NNJ4=3, NMJ4=2*NNJ4-1 ! NMJ4=5
   REAL,    PARAMETER :: SQPI =1.772453851    ! sqrt(pi)
   REAL,    PARAMETER :: GAMMA=0.577215665    ! Euler Gamma
   REAL,    PARAMETER :: LOG4 =1.386294361    ! ln(4)
   REAL, PARAMETER :: J4FF =-2./49.                    !-0.0408163E0
   REAL, PARAMETER :: J4MM = 4./25.                    ! 0.16 
   REAL, PARAMETER :: J4BB =-4./9.                     !-0.4444444E0
   REAL, PARAMETER :: J4AA1=-SQPI/2.                   !-0.8862269E0  
   REAL, PARAMETER :: J4AA2=SQPI*(2.-GAMMA-LOG4)/2.    ! 3.23383974E-02 
   REAL, PARAMETER :: J4CC = 0.1E0
   REAL, PARAMETER :: J4EE = 0.1E0 
   REAL,    PARAMETER   :: XMIN=0.18D0, XMAX=4.7D0
   REAL    :: y, y3
!
! analytic asymptotic forms
!
  IF (x .LE. XMIN) THEN
     j4=J4MM*x+J4BB
  ELSEIF (x .GE. XMAX) THEN
     y=SQRT(x)                   ! x^1/2
     y3=x*y                      ! x^3/2
     j4=(J4AA1*LOG(x) + J4AA2)/y3
  ELSE 
     j4=-0.446059 - 0.0023836*SQRT(x) + 0.207461*x - &
        0.0883521*x**1.5 +  0.011114*x**2
  ENDIF
 END FUNCTION j4
!
! The Dawson function takes the form
!
!
!          /x
!          |               
! daw(x) = | dy e^{y^2 -x^2}
!          |              
!          /0
!      
!        =(sqrt(pi)/2)*exp(-x^2)*erfi(x)
!
! daw() approximates Dawson's integral by rational
! functions with coefficients. At large and small
! arguments we use analytic approximations.
!
!
! For small x < XMIN we use the asymptotic form
!
!                 2x^3     4x^5
! daw(x) =  x  +  ----  + ----- + O(x^7) 
!                  3        15
!
! and for large x > XMAX we use 
!
!            1     1       3
! daw(x) =  --- + ---- + ----- + O(x^-7) 
!           2x    4x^3    8x^5
!
! The error is 0.03% for XMIN=0.4 and 0.01% XMAX=5.0. For 
! intermediate values, we approximate daw(x) as a rational 
! function of the form
!
!              x      x^6+b5*x^5+b4*x^4+b3*x^3+b2*x^2+b1*x+b0
! daw(x) = --------- ----------------------------------------
!          2 x^2 + 1  x^6+a5*x^5+a4*x^4+a3*x^3+a2*x^2+a1*x+b0
!
!
! With the values of bn and an chosen below, the error is 0.03%.
!
 FUNCTION daw(x)
!
! As x->0 and x->infty, the Dawson function takes the asymptotic
! forms daw(x)~x and daw(x)~1/(2x), respectively. The first 
! rational function "R(x)=x/(2x**2+1)" reproduces this behavior; 
! the 6-th order polynomial-ratio Q6(x) asymptotes to one at both 
! ends (to preserve the asymptotic form of the previous function), 
! with the coefficients a5, b5, ... b0 being chosen to provide 
! agreement with the exact Dawson integral at the values:
!
! x0=0.92413 daw(x0)=0.541044
! x1=0.2     daw(x1)=0.194751
! x2=0.5     daw(x2)=0.424436
! x3=0.7     daw(x3)=0.510504
! x4=1.2     daw(x4)=0.507273
! x5=1.4     daw(x5)=0.456507
! x6=1.6     daw(x6)=0.399940
! x7=2.0     daw(x7)=0.301340
! x8=3.0     daw(x8)=0.178271
! x9=4.0     daw(x9)=0.129348
! x10=8.0    daw(x10)=0.0630002
!
! See daw.nb for details.
!
   USE bpsvars
   IMPLICIT NONE
   REAL, INTENT(IN) :: x
   REAL  daw
   INTEGER, PARAMETER          :: NNDAW=3, NMDAW=2*NNDAW-1 ! NMDAW=5
   REAL,    PARAMETER, DIMENSION(0:NMDAW) :: DWB=(/ &
   5.73593880244318E0 ,  & !b0
   -6.73666007137766E0,  & !b1
   1.99794422787154E1 ,  & !b2
   -1.85506350260761E1,  & !b3
   1.22651360905700E1 ,  & !b4
   -4.67285812684807E0/)   !b5
   REAL,    PARAMETER, DIMENSION(0:NMDAW) :: DWA=(/ &
   5.73593880244318E0 ,  & !a0=b0
   -6.82372048950896E0,  & !a1
   1.33804115903096E1 ,  & !a2
   -1.42130723670491E1,  & !a3
   1.11714434417979E1 ,  & !a4
   -4.66303387468937E0/)   !a5
   REAL,    PARAMETER  :: XMIN=0.4D0, XMAX=5.D0
   REAL    :: x3, x5, xx, ra, rc
   INTEGER :: n

   IF (x .LE. XMIN) THEN
      x3=x*x*x
      x5=x3*x*x
      daw=x - 2.D0*x3/3.D0 + 4.D0*x5/15.D0
   ELSEIF (x .GE. XMAX) THEN
      x3=x*x*x
      x5=x3*x*x
      daw=1.D0/(2.D0*x)+1.D0/(4.D0*x3)+3.D0/(8.D0*x5)
   ELSE
      ra=0.E0
      rc=0.E0
      xx=1.E0
      DO n=0,NMDAW
         ra=ra+DWA(n)*xx
         rc=rc+DWB(n)*xx
         xx=x*xx
      ENDDO
      ra=ra+xx
      rc=rc+xx
      daw=x/(1.E0+2.E0*x*x)
      daw=daw*rc/ra
   ENDIF
 END FUNCTION daw
!
!======================================


!======================================
! general purpose routines
!
!
      SUBROUTINE param(nni, betab, nb, gb, ge, gi, etae, ze)
      USE bpsvars
      IMPLICIT NONE
      INTEGER,                     INTENT(IN)  :: nni  !  Number of ions
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: nb   !  Number density array [cm^-3]
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: betab!  Temperature array    [1/keV]
      REAL,    DIMENSION(1:nni+1), INTENT(OUT) :: gb   !  Plasma coupling array
      REAL                       , INTENT(OUT) :: ge   !  Electron coupling
      REAL                       , INTENT(OUT) :: gi   !  Average ion coupling
      REAL                       , INTENT(OUT) :: etae !  Quantum parameter for electron
      REAL                       , INTENT(OUT) :: ze   !  Fugacity
      REAL,    PARAMETER :: GECOEFF=6.1260E-15         !  g-coefficient
      REAL,    PARAMETER :: ETACOEF=0.164961           !  eta-coeff
      REAL,    PARAMETER :: ZECOEFF=5.2359E-27         !  ze-coefficient
      REAL,    PARAMETER :: MEOMP  =1836.15267247      !  me/mp ratio
      REAL,    PARAMETER :: HBARCMEVFM=197.3           !  hbar-c [MeV-fm]
      REAL,    PARAMETER :: HBARCKEVCM=1.973E-8        !  hbar-c [keV-cm]
      REAL    :: te, ti, ne
      INTEGER :: ib

      te=1./betab(1)
      ti=1./betab(2)
      ne=nb(1)
!
! Plasma Coupling
!
!     kb2 =8*PI*A0CM*BEKEV*betab*zb*zb*nb
!     gb=2*BEKEV*betab*A0CM*SQRT(kb2)
      gb=GECOEFF*SQRT(nb)*betab**1.5
      ge=gb(1)
      gi=0
      DO ib=2,nni+1
        gi=gi+gb(ib)**2
      ENDDO
      gi=SQRT(gi/nni)
!
! Quantum coefficient
!
      etae=ETACOEF/SQRT(te)
!
! degeneracy coefficient (electron fugacity)
!      
      ze=ZECOEFF*ne/te**1.5

      END SUBROUTINE param


!
! n=n(g,T) for g-contours in T-n plane
!
      SUBROUTINE n_tg(t, g, n)
      IMPLICIT NONE        
        REAL, INTENT(IN)  :: t  ! keV
        REAL, INTENT(IN)  :: g  ! 
        REAL, INTENT(OUT) :: n  ! 1/cc
        REAL, PARAMETER   :: COEFF=6.1216E-15
        n=((g/COEFF)**2)*(t**3) ! 1/cc
      END SUBROUTINE n_tg

!
! ne=ne(ze,T) for ze-contours in T-n plane
!
      SUBROUTINE n_tze(t, ze, n)
      IMPLICIT NONE        
        REAL, INTENT(IN)  :: t  ! keV
        REAL, INTENT(IN)  :: ze ! 
        REAL, INTENT(OUT) :: n  ! 1/cc
        REAL, PARAMETER   :: COEFF=5.2359E-27
        n=ze*(t**1.5)/COEFF ! 1/cc
      END SUBROUTINE n_tze

!
! holding patterns ....
!

      SUBROUTINE frfi_p(x, alpha, fr, fi, fabs, farg)
      IMPLICIT NONE
      REAL,    INTENT(IN) :: x, alpha
      REAL                :: fr, fi, fabs, farg
      REAL                :: daw
      REAL,    PARAMETER  :: SQPI =1.772453851    ! sqrt(pi)
      REAL,    PARAMETER  :: PI   =3.141592654    ! pi
      fr=1 + alpha*(1-2*x*daw(x))
      fi=alpha*SQPI*x*EXP(-x*x)
      fabs=SQRT(fr*fr + fi*fi)
      farg=ATAN2(fi,fr)
      END SUBROUTINE frfi_p

      SUBROUTINE frfi_one(x, fr, fi, fabs, farg)
      IMPLICIT NONE
      REAL,    INTENT(IN) :: x
      REAL                :: fr, fi, fabs, farg
      REAL                :: daw
      REAL,    PARAMETER  :: SQPI =1.772453851    ! sqrt(pi)
      REAL,    PARAMETER  :: PI   =3.141592654    ! pi
      fr=1-2*x*daw(x)
      fi=SQPI*x*EXP(-x*x)
      fabs=SQRT(fr*fr + fi*fi)
      farg=ATAN2(fi,fr)
      END SUBROUTINE frfi_one
!
!======================================
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


! 
! Robert Singleton
! - Santa Fe, Winter 2005 [find correct date]
! - Santa Fe, March 2009    [start  rewrite]
! - Santa Fe, November 2009 [finish rewrite]
! 
! ROUTINE: dedx_bps(nni, ep, zp, mp, betab, zb, mb, nb,             &
!            dedx_tot, dedx_i, dedx_e, dedxc_tot, dedxc_i, dedxc_e, & 
!            dedxq_tot, dedxq_i, dedxq_e)
! 
! Assume a plasma composed of several species b, each separately in 
! thermal equilibrium with themselves but not necessarily with each 
! other. This routine returns several useful components of the BPS
! stopping power dE/dx. 
! 
! UNITS: dE_b/dx has units of [MeV/micron] (subject to change in updates)
! 
! THE PHYSICS:
! See notes in acoeff.f90.
!
! INPUT: nni, betab, mb, nb, zb, ep, mp, zp
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
!
! OUTPUT: dedx_tot, dedx_i, dedx_e, dedxc_tot, dedxc_i, dedxc_e
!         dedxq_tot, dedxq_i, dedxq_e
!
! Each plasma component b makes a linear contribution dE_b/dx to the 
! total stopping power, i.e. dE/dx = sum_b dE_b/dx [this is true only
! to leading and next-to-leading order in the plasma coupling g]. 
! Each dE_b/dx in turn can be be decomposed into a classical-quantum 
! or electron-ion contributions.
!
! classical electron  : dedxc_e
! classical ion       : dedxc_i [sum over all ions]
! classical total     : dedxc_tot = dedxc_e + dedxc_i 
! quantum   electron  : dedxq_e
! quantum   ion       : dedxq_i [sum over all ions]
! quantum   total     : dedx_tot = dedxq_e + dedxq_i
! total     electric  : dedx_e = dedxc_e + dedxq_e
! total     ion       : dedx_i = dedxc_i + dedxq_i
! total               ! dedx_tot = dedx_e + dedx_i
!
      SUBROUTINE dedx_bps(nni, ep, zp, mp, betab, zb, mb, nb,  &
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

        REAL,    PARAMETER           :: units=1.E-7 ! 10^-7 => MeV/micron; 10^-4 keV/micron
        REAL,    DIMENSION(1:nni+1)  :: kkb2, aab, ab2, c1b, c2b, c3b,c4b, eetb
        REAL,    DIMENSION(1:nni+1)  :: mpb, mbpb, mtotb, mtotpb, mratpb, bpb, cpb
        REAL     :: vp, vp2, zp2, kk, k2, kd, kd2, c1, c2, c3, c4, eta, a, b, c, a2
        REAL     :: dedxc_s, dedxc_r, dedxc_r1, dedxc_r2, dedxq, mt, mr, xx
        REAL     :: dedxc_s_e, dedxc_r_e, dedxc_s_i,dedxc_r_i
!       REAL     :: ln_dedx_bps
        INTEGER  :: ib, NNB
!
! initialize components of dE/dx
!
        dedx_tot = 0  ! electron + ion
        dedx_i   = 0  ! ion contribution
        dedx_e   = 0  ! electron contribution
        dedxc_tot= 0  ! classical total
        dedxc_e  = 0  ! classical electron
        dedxc_i  = 0  ! classical ion
        dedxq_tot= 0  ! quantum total
        dedxq_e  = 0  ! quantum electron
        dedxq_i  = 0  ! quantum ion
        dedxc_s_i=0 
        dedxc_s_e=0
        dedxc_r_i=0
        dedxc_r_e=0
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
!       [ betab mbc2 ]^1/2   vp
! c2 =  |----------- |      ---  [dimensionless]
!       [   2 Pi     ]      c^2
!        
! where
!          e^2 
! Be  =  --------- = 13.606E-3   [keV]
!        8 Pi a0                 Bohr radius: a0=5.29E9 cm
!
        NNB  = nni+1                 ! number of ions + electrons
        vp   = CC*SQRT(2*ep/mp)      ! [cm/s]
        vp2  = vp*vp                 ! [cm^2/s^2]
!       kkb2 = DEBYE2*zb*zb*nb*betab ! [1/cm^2]
        kkb2 = 8*PI*BEKEV*betab*zb*zb*nb*A0CM ! inv Debye length squared [a0^-2]
        kd2  = SUM(kkb2)             ! [1/cm^2]
        kd   = SQRT(kd2)             ! [1/cm]
        k2   = kkb2(1)               ! [1/cm^2]
        kk   = SQRT(k2)              ! [1/cm]   kk = k_e
!*!
        kk=kd
!*!
        zp2  = zp**2                 ! [dimensionless]
        mpb  = mp*mb/(mp+mb)         ! [keV]
        mbpb = mb/mpb                ! [dimensionless]
        mtotb= mp + mb               ! [keV]
        mtotpb=mtotb/mp              ! [dimensionless]
        mratpb=mb/mp                 ! [dimensionless]
        aab    =0.5*betab*mb*vp2/CC2 ! A_b=(1/2)betab*mbc2*vp2/CC2
        ab2   = SQRT(aab)            ! sqrt(ab)
        bpb   =betab*mtotb*vp2/CC2   ! B_b
        cpb   =2 - 2*GAMMA - &       ! C_b
              LOG(2*BEKEV*ABS(zp*zb)*A0CM*kk*betab*mtotpb)
        eetb  =2*BEKEV*ABS(zp*zb)*A0CM/HBARC/(vp/CC) ! quantum parameter

        c1b=2*zp2*BEKEV*kkb2*A0CM     ! [keV/cm]
        c1b=c1b*1.E-7                 ! [MeV/micron]  
        c2b=SQRT(aab/PI)              ! [dimensionless] 
                                      ! c2b=SQRT(betab*mb/TWOPI)*vp/CC 
        c3b=1/(2*aab)                 ! c3b=1/betab*mb*vp^2
        c4b=1/(betab*mp*vp2/CC2)      ! c4b
!
! Loop over charged plasma species
!
        DO ib=1,nni+1
        IF (zb(ib) .NE. 0.) THEN
           c1 =c1b(ib)   ! [MeV/micron]
           c2 =c2b(ib)   ! [dimensionless]
           c3 =c3b(ib)   ! [dimensionless]
           c4 =c4b(ib)   ! [dimensionless]

           a2 =ab2(ib)   ! [dimensionless] ab2=sqrt(aab)
           b  =bpb(ib)
           c  =cpb(ib)

           a  =aab(ib )   ! [dimensionless] 
           eta=eetb(ib)   ! [dimensionless] quantum parameter
           mt =mtotpb(ib) ! M_{pb}/m_p=(m_p + m_b)/m_p
           mr =mratpb(ib) ! m_b/m_p
!
! A-classical-singular 
!
           CALL dedx_sing(a2,b,c,dedxc_s) 
           dedxc_s=c1*c2*c4*dedxc_s
!
! A-classical-regular: 
!
! NOTE: Eventually the regular contribution needs to be re-written.
! dedx_reg1 and only interfaces to the old dedx_reg code.
! 
! 
!
           CALL dedx_reg1(nni,ib,aab,kkb2,kk,dedxc_r1) 
           xx=-BEKEV*zp2/(PI*betab(ib)*mp*(vp2/CC2)*A0CM)*1.E-7
           dedxc_r1=dedxc_r1*xx ! [MeV/micron] xx=cp3/vth2/rb(ib)
!             
           CALL dedx_reg2(nni,ib,aab,kkb2,kk,dedxc_r2) 
           xx=(BEKEV*zp2/TWOPI/A0CM)*1.E-7
           dedxc_r2=dedxc_r2*xx
!*!
!PRINT *, "xxRef2:", ib, dedxc_r2
!*! 
!
           dedxc_r=dedxc_r1 + dedxc_r2
!
! A-quantum
!
           CALL dedx_quantum(a,eta,mt,mr,dedxq)
           dedxq=-c1*c2*c3*dedxq  ! [MeV/micron]
!
! coulomb logarithm
!
!*!
!          ln_dedx_bps=dedxc_s + dedxc_r + dedxq
!          ln_dedx_bps=LOG(ln_dedx_bps/(c1*c3))
!          PRINT *, ib, ln_dedx_bps
!*!
           CALL x_collect(ib,NNB,dedxc_s,dedxc_r,dedxq,dedx_tot,dedx_i,dedx_e,&
             dedxc_tot,dedxc_i,dedxc_e,dedxq_tot,dedxq_i,dedxq_e, dedxc_s_i,  & 
             dedxc_s_e, dedxc_r_i, dedxc_r_e)

        ENDIF
        ENDDO
      END SUBROUTINE dedx_bps
!
!
! quantum correction:
!
! This subroutine calculates the quantum correction dedxq.
!
!           /Infinity
!           |
! dedxq  =  | du  d_dedxq(u)
!           |
!           /0
!
! a = ab(ib)      = (1/2) betab*mbc2*(vp/c)^2
! e = eetb(ib)     = ep*ez/4pi hbar vp
! mt= M_pb(ib)/mp = mtotpb(ib)
! mr= mb(ib)/mp   = mratpb(ib)
!
      SUBROUTINE dedx_quantum(a, e, mt, mr, dedxq)
      USE bpsvars
      IMPLICIT NONE
        REAL, INTENT(IN) :: a
        REAL, INTENT(IN) :: e
        REAL, INTENT(IN) :: mt
        REAL, INTENT(IN) :: mr
        REAL, INTENT(OUT):: dedxq
!  
        REAL,    PARAMETER :: UPM=0.7745966692E0
        REAL,    PARAMETER :: W13=0.5555555556E0, W2=0.8888888889E0
        INTEGER, PARAMETER :: NG=10000  ! must be even
        REAL,    PARAMETER :: NN=30.E0
        REAL               :: u0, u1, du
        REAL               :: u, um, d_dedxq
        INTEGER            :: iu
        u0=1.E0 - NN/SQRT(a) ! u0 >= 0 
        u0=MAX(0.,u0)        ! 
        u1=1.E0 + NN/SQRT(a)
        du=(u1-u0)/NG
        dedxq=0.E0
        u=u0-du
        DO iu=1,NG,2
!     
           u=u+2.E0*du
           dedxq=dedxq+W2*d_dedxq(a,e,mt,mr,u)
!
           um=u-du*UPM
           dedxq=dedxq+W13*d_dedxq(a,e,mt,mr,um)
!
           um=u+du*UPM
           dedxq=dedxq+W13*d_dedxq(a,e,mt,mr,um)
        ENDDO
        dedxq=dedxq*du
      END SUBROUTINE dedx_quantum
!
      FUNCTION d_dedxq(a, e, mt, mr, u)
      IMPLICIT NONE
        REAL, INTENT(IN) :: a, e, mt, mr, u
        REAL, PARAMETER  :: AXMAX=0.05/10
        REAL             :: d_dedxq
        REAL             :: repsi, sh, ch
        REAL             :: ep, em, um1, up1, eu, a2u
        eu =ABS(e/u)
        a2u=2*a*u
        um1=u-1
        up1=u+1
        ep=EXP(-a*up1*up1)
        em=EXP(-a*um1*um1)
        sh=0.5E0*(em-ep)
!       IF (a .GT. AXMAX) THEN   ! sh and ch are 
          ch=0.5E0*(em+ep)       ! not sinh and cosh
          ch=(ch - sh/a2u)/u
!       ELSE                     ! small argument
!         a2=a*a                 ! limit: gives wrong result for electrons: fix later
!         a3=a2*a
!         a4=a2*a2
!         ch=4*a2*u/3 + (8*a4/15. - 4*a3/3.)*u*u*u
!         ch=ch*EXP(-a)
!      ENDIF
       d_dedxq=2*(repsi(eu) - LOG(eu))
       d_dedxq=d_dedxq*(mt*ch - mr*sh)
      END FUNCTION d_dedxq

!
! classical: singular short distance contribution
!
!
      SUBROUTINE dedx_sing(a, b, c, dedxc_s)  ! I_1(a,b,c)
      USE mathvars
      IMPLICIT NONE
        REAL, INTENT(IN)  :: a ! SQRT(A_b)
        REAL, INTENT(IN)  :: b ! 
        REAL, INTENT(IN)  :: c ! 
        REAL, INTENT(OUT) :: dedxc_s
        REAL    :: bc, a2, a3, erfa, expa, ferf
        REAL    :: j1, j2, j3, j4
        bc=b*c
        a2=a*a  ! A_b
        a3=a2*a
        erfa=SQPI*ferf(a)  ! see ferf.f
        expa=EXP(-a2)
        dedxc_s=erfa*((2-c)/a + bc/(2*a3))-bc*expa/a2
        dedxc_s=dedxc_s + j3(a2) - j1(a2) + b*(j2(a2) -j4(a2))
      END SUBROUTINE dedx_sing
!
      SUBROUTINE dedx_sing2(a, b, c, dedxc_s) ! I_1(a,b,c) via Gaussian quadrature.
      USE bpsvars
      IMPLICIT NONE
        REAL, INTENT(IN)  :: a, b, c
        REAL, INTENT(OUT) :: dedxc_s
  
        REAL,    PARAMETER :: UPM=0.7745966692E0
        REAL,    PARAMETER :: W13=0.5555555556E0, W2=0.8888888889E0
        INTEGER, PARAMETER :: NG=10000  ! must be even
        REAL               :: u0, u1, du
        REAL               :: u, um, d_dedxc_s2
        INTEGER            :: iu
        u0=0
        u1=1
        du=(u1-u0)/NG
        dedxc_s=0.E0
        u=u0-du
        DO iu=1,NG,2
!     
           u=u+2.E0*du
           dedxc_s=dedxc_s+W2*d_dedxc_s2(a,b,c,u)
!
           um=u-du*UPM
           dedxc_s=dedxc_s+W13*d_dedxc_s2(a,b,c,um)
!
           um=u+du*UPM
           dedxc_s=dedxc_s+W13*d_dedxc_s2(a,b,c,um)
        ENDDO
        dedxc_s=dedxc_s*du
      END SUBROUTINE dedx_sing2
!
      FUNCTION d_dedxc_s2(a, b, c, u)
      IMPLICIT NONE
        REAL, INTENT(IN) :: a, b, c, u
        REAL             :: d_dedxc_s2
        d_dedxc_s2= (c - LOG(u/(1-u)))*(b*SQRT(u) - 1/SQRT(u))
        d_dedxc_s2=d_dedxc_s2 + 2/SQRT(u)
        d_dedxc_s2=d_dedxc_s2*Exp(-a*a*u)
      END FUNCTION d_dedxc_s2
!
! regular contribution
!
!
! classical: regular long distance contribution
!
!
      SUBROUTINE dedx_reg1(nni, ib, aab, kkb2, kk, dedxc_r)
      USE physvars
      IMPLICIT NONE
        INTEGER,                     INTENT(IN)  :: nni
        INTEGER,                     INTENT(IN)  :: ib
        REAL,   DIMENSION(1:nni+1),  INTENT(IN)  :: aab  ! don't use bpsvars
        REAL,   DIMENSION(1:nni+1),  INTENT(IN)  :: kkb2 ! don't use bpsvars
        REAL,                        INTENT(IN)  :: kk  ! don't use bpsvars
        REAL,                        INTENT(OUT) :: dedxc_r
        REAL,   DIMENSION(1:nni+1)  :: abv  ! SQRT(abv)
        REAL,   DIMENSION(1:nni+1)  :: ka0b2
        REAL                        :: ka0
        REAL                        :: hii
        abv=SQRT(aab)
        ka0=kk*A0CM
        ka0b2=kkb2*A0CM**2
        CALL hi(nni,ib,ka0,ka0b2,abv,hii)
        dedxc_r=hii
      END SUBROUTINE dedx_reg1
!
!                 /1
!                 |
! dedx_reg2 = 2 * | du u*H(u*abv)
!                 |
!                 /0
!
      SUBROUTINE dedx_reg2(nni, ib, aab, kkb2, kk, dedxc_r)
      USE physvars
      IMPLICIT NONE
        INTEGER,                     INTENT(IN)  :: nni
        INTEGER,                     INTENT(IN)  :: ib
        REAL,   DIMENSION(1:nni+1),  INTENT(IN)  :: aab  ! don't use bpsvars
        REAL,   DIMENSION(1:nni+1),  INTENT(IN)  :: kkb2 ! don't use bpsvars
        REAL,                        INTENT(IN)  :: kk  ! don't use bpsvars
        REAL,                        INTENT(OUT) :: dedxc_r
        REAL, PARAMETER :: UPM=0.7745966692E0
        REAL, PARAMETER :: W13=0.5555555556E0, W2=0.8888888889E0
   
        INTEGER, PARAMETER :: NG=2000  ! NG must be even
        REAL,    PARAMETER :: U0=0.E0, U1=1.E0, DU=1.E0/NG


        REAL,   DIMENSION(1:nni+1)  :: abv  ! SQRT(abv)
        REAL,   DIMENSION(1:nni+1)  :: ka0b2, uu
        REAL                        :: ka0
        REAL                        :: hii, u, um
        INTEGER :: iu
        abv=SQRT(aab)
        ka0=kk*A0CM
        ka0b2=kkb2*A0CM**2    !/CMTOA0/CMTOA0

        dedxc_r=0
        u=U0-DU
        DO iu=1,NG,2
!
           u=u+2.E0*DU
           uu=u*abv
           CALL hi(nni,ib,ka0,ka0b2,uu,hii)
           dedxc_r=dedxc_r+W2*u*hii
!
           um=u-DU*UPM
           uu=um*abv
           CALL hi(nni,ib,ka0,ka0b2,uu,hii)
           dedxc_r=dedxc_r+W13*um*hii
!
           um=u+DU*UPM
           uu=um*abv
           CALL hi(nni,ib,ka0,ka0b2,uu,hii)
           dedxc_r=dedxc_r+W13*um*hii
        ENDDO
        dedxc_r=2*dedxc_r*DU
      END SUBROUTINE dedx_reg2
!
! The ib-th component of Hb is hi=Sqrt[Pi]*kb^2*xb*Exp[-xb^2]*H/Fi 
! with Fi=Sqrt[Pi]Sum_c kc^2*xc*Exp[-xc^2] where xb=abv(ib) and 
! kb^2=kkb2(ib)
!
      SUBROUTINE hi(nni, ib, kk, kkb2, abv, hii)
        IMPLICIT NONE
        INTEGER,                     INTENT(IN)  :: nni
        INTEGER,                     INTENT(IN)  :: ib
        REAL,                        INTENT(IN)  :: kk
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: kkb2
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: abv
        REAL,                        INTENT(OUT) :: hii
        REAL, PARAMETER :: EPS=1.E-6, XB2MAX=10.D0
        REAL    :: fr, fi, fabs, farg, h, ebc
        REAL    :: xb, xb2, xc, xc2, xbc
        INTEGER :: ic
        CALL fri(nni,abv,kkb2,fr,fi,fabs,farg)
        h=-2*(fi*LOG(fabs/kk**2) + fr*farg)
        xb =abv(ib)
        xb2=xb*xb
!
! Express Fi in terms of its sum and write 
!
!   hi = kb^2*xb*Exp[-xb^2]/Sum_c kc^2*xc*Exp[-xc^2]
!
! then calculate hi^-1 first, excluding kb^2*xb since
! some components of kb^2*xb might be zero.
        hii=0                                       
        DO ic=1,nni+1
           xc =abv(ic)                             
           xc2=xc*xc                              
           xbc=xb2-xc2                            
           IF (xbc .LT. XB2MAX) THEN                
              ebc=EXP(xb2-xc2)                    
           ELSE                                   
              hii=0   ! ebc is large and dominates the sum, 
              RETURN  ! and its inverse is almost zero
           ENDIF
           hii=hii + kkb2(ic)*xc*ebc                 
        ENDDO
        hii=kkb2(ib)*xb/hii  ! invert and include kb^2*xb
        hii=h*hii
      END SUBROUTINE hi
!
! Returns the dielectric function F in terms of the 
! real part, the imaginary part, the absolute value, 
! and the argument: fr, fi, fabs, farg
!
      SUBROUTINE fri(nni, xb, kb2, fr, fi, fabs, farg)
      USE mathvars
      IMPLICIT NONE
        INTEGER,                 INTENT(IN)  :: nni
        REAL, DIMENSION(1:nni+1),INTENT(IN)  :: xb, kb2
        REAL,                    INTENT(OUT) :: fr,  fi, fabs, farg
        REAL    :: x, daw, d
        INTEGER :: ib
        fr=0
        fi=0
        DO ib=1,nni+1
           x=xb(ib)
           d=daw(x)
           fr=fr + kb2(ib)*(1-2*x*d)
           fi=fi + kb2(ib)*x*EXP(-x*x)
        ENDDO
        fi=fi*SQPI
        fabs=SQRT(fr*fr + fi*fi)
        farg=ATAN2(fi,fr)
      END SUBROUTINE fri
!
! ============================================
!


!
!             /Infinity
!             |
! kqm1(x)  =  | dy   y exp(-y^2) [ln(x/y) - repsi(x/y)]
!             |
!            /0
!
!
 FUNCTION kqm1(x)  
   IMPLICIT NONE
   REAL, INTENT(IN) :: x
   REAL :: kqm1

   REAL, PARAMETER :: XMIN=0.15E0, XMAX=3.2E0 

   REAL, PARAMETER :: a0= 0.4329117486761496454549449429 ! 3*GAMMA/4
   REAL, PARAMETER :: a1= 1.2020569031595942854250431561 ! ZETA(3)
   REAL, PARAMETER :: a2= 0.1487967944177345026410993331 ! ZETA'(3)+GAMMA*
   REAL, PARAMETER :: b2=-0.0416666666666666666666666667 !-1/24  ZETA(3)/2
   REAL, PARAMETER :: b4=-0.0083333333333333333333333333 !-1/120
   REAL, PARAMETER :: b6=-0.0119047619047619047619047619 !-1/84
   REAL, PARAMETER :: c0= 0.25109815055856566000
   REAL, PARAMETER :: c1=-0.02989283169766254700
   REAL, PARAMETER :: c2= 0.03339880139150325000
   REAL, PARAMETER :: c3=-0.00799128953390392700
   REAL, PARAMETER :: c4= 0.00070251863606810650
   REAL, PARAMETER :: d0=-0.18373957827052560000
   REAL, PARAMETER :: d1=-0.33121125339783110000
   REAL, PARAMETER :: d2= 0.04022076263527408400
   REAL, PARAMETER :: d3=-0.00331897950305779480
   REAL, PARAMETER :: d4= 0.00012313574797356784
   REAL :: x2, lx, xi
   IF ( x .LE. XMIN) THEN                        ! x < xmin=0.15: to 0.06%
      x2=x*x                                     ! ln(x)/2 + 3*GAMMA/4 +
      lx=LOG(x)                                  ! ZETA(3)*X^2*ln(x) +
      kqm1=0.5E0*lx + a0 + a1*x2*lx + a2*x2      ! [ZETA'(3) + GAMMA*
                                                 ! ZETA(3)/2]*x^2
                                                 !
   ELSEIF ( x .GE. XMAX ) then                   ! x > xmax=3.2: to 0.12%
      xi=1/x                                     ! -1/24*x^2 - 1/120*x^4 -
      x2=xi*xi                                   ! 1/84*x^6
      kqm1=b6                                    !
      kqm1=b4 + kqm1*x2                          !
      kqm1=b2 + kqm1*x2                          !
      kqm1=kqm1*x2                               !
   ELSE
      xi=1/x                                     ! xmin < x , xmax
      lx=LOG(x)                                  ! fit accurate to 0.2%
      kqm1=c4                                    ! c0 + c1*x + c2*x^2 +
      kqm1=c3+kqm1*x                             ! c3*x^3 + c4*x^4 +
      kqm1=c2+kqm1*x                             ! d0*ln(x) +
      kqm1=c1+kqm1*x                             ! d1/x + d2/x^2 +
      kqm1=c0+kqm1*x + d0*lx                     ! d3/x^3 + d4/x^4
      lx=d4                                      !
      lx=d3+lx*xi                                !
      lx=d2+lx*xi                                !
      lx=d1+lx*xi                                !
      lx=lx*xi                                   !
      kqm1=kqm1+lx
   ENDIF
 END FUNCTION kqm1
!
!
!             /Infinity
!             |
! kqm3(x)  =  | dy   y^3 exp(-y^2) [ln(x/y) - repsi(x/y)]
!             |
!            /0
!
!
 FUNCTION kqm3(x)  
   IMPLICIT NONE
   REAL, INTENT(IN) :: x
   REAL :: kqm3

   REAL, PARAMETER :: XMIN=0.15E0, XMAX=2.5E0 
   REAL, PARAMETER :: a0= 0.1829117486761496454549449429 ! 3*GAMMA/4 - 1/4
   REAL, PARAMETER :: a2=-0.6010284515797971427073328102 !-ZETA(3)/2
   REAL, PARAMETER :: b2=-0.0833333333333333333333333333 !-1/12
   REAL, PARAMETER :: b4=-0.025                          !-1/40
   REAL, PARAMETER :: b6=-0.046875                       !-3/64
   REAL, PARAMETER :: c0= 0.691191700599840900000
   REAL, PARAMETER :: c1=-1.094849572900974000000
   REAL, PARAMETER :: c2= 0.318264617154601400000
   REAL, PARAMETER :: c3=-0.060275957444801354000
   REAL, PARAMETER :: c4= 0.005112428730167831000
   REAL, PARAMETER :: d0= 0.835543536679762600000
   REAL, PARAMETER :: d1= 0.047821976622976340000
   REAL, PARAMETER :: d2= 0.000053594881446931025
   REAL, PARAMETER :: d3=-0.000268040997573199600
   REAL, PARAMETER :: d4= 0.000015765134162582942
   REAL :: x2, lx, xi
   IF ( x .LE. XMIN) THEN                        ! x < xmin=0.15: to 0.1%
      x2=x*x                                     ! ln(x)/2 + 3*GAMMA/4 -1/4
      lx=LOG(x)                                  ! -[ZETA(3)/2]*x^2
      kqm3=0.5E0*lx + a0 + a2*x2                 ! 
                                                 ! 
                                                 !
   ELSEIF ( x .GE. XMAX ) then                   ! x > xmax=2.5: to 0.25%
      xi=1/x                                     ! -1/12*x^2 - 1/40*x^4 -
      x2=xi*xi                                   ! 3/64*x^6
      kqm3=b6                                    !
      kqm3=b4 + kqm3*x2                          !
      kqm3=b2 + kqm3*x2                          !
      kqm3=kqm3*x2                               !
   ELSE
      xi=1/x                                     ! xmin < x < xmax
      lx=LOG(x)                                  ! fit accurate to 0.04%
      kqm3=c4                                    ! c0 + c1*x + c2*x^2 +
      kqm3=c3+kqm3*x                             ! c3*x^3 + c4*x^4 +
      kqm3=c2+kqm3*x                             ! d0*ln(x) +
      kqm3=c1+kqm3*x                             ! d1/x + d2/x^2 +
      kqm3=c0+kqm3*x + d0*lx                     ! d3/x^3 + d4/x^4
      lx=d4                                      !
      lx=d3+lx*xi                                !
      lx=d2+lx*xi                                !
      lx=d1+lx*xi                                !
      lx=lx*xi                                   !
      kqm3=kqm3+lx
   ENDIF
 END FUNCTION kqm3
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

