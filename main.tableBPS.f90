      PROGRAM range
      USE allocatablevars
      USE bpsvars
      USE physvars
      USE mathvars
        IMPLICIT NONE
        INTEGER NTSMAX
        PARAMETER (NTSMAX=5000)

        REAL    :: dedx_tot,  dedx_i,  dedx_e
        REAL    :: dedxc_tot, dedxc_i, dedxc_e
        REAL    :: dedxq_tot, dedxq_i, dedxq_e

        REAL    :: te, ti, ne, mp, zp
        REAL    :: epp
        INTEGER :: nni
        REAL    :: ddee, E0, f, mDT
        REAL    :: rho(0:72), TTee(0:54), EeV(0:100)
        INTEGER :: ntabrho, ntabte, ntabene, ir, it, ie, id, iu
        REAL, DIMENSION(0:72,0:54,0:100) :: tabi, tabe
!
! Range parameters
!
!        tmax=1.0E-11 ! [s]
!        nts=100
!        nit=10
!        nn =300
!        dt=tmax/(nts*nit)
!        te=3.        ! Electron temperature       [keV]
!        ti=3.        ! Ion temperature            [keV]
!        ne=1.E25     ! Electron number density    [cm^-3]
!        ep=3540.     ! Projectile energy  [keV]
         mp=4*MPKEV    ! Projectile mass    [keV]
         zp=2.         ! Projectile charge  [e]


        ntabrho = 72
        do ir = 0, ntabrho
           id = int( float(ir) / 9. )
           iu = ir - 9 * id
           rho(ir) = 0.0001 * (iu+1)*10.**id ! DT densities [g/cc]: 0.0001 - 10000
        enddo

        ntabte  = 54
        do it = 0, ntabte
           id = int( float(it) / 9. )
           iu = it - 9 * id
           TTee(it) = (iu+1)*10.**id ! Temperatures Te = Ti [eV]: 1 - 1000 keV
        enddo

        ntabene = 100
        E0 = 17.6e6 / 5.
        ! eV
        f = 10.**(-2./(ntabene-1.))
        ddee  = E0 * ( f - 1. ) / ( f**ntabene - 1. )
        EeV(0) = E0
        do ie = 1, ntabene - 1
           EeV(ie) = EeV(ie-1) - ddee * f**(ie-1) ! alpha energies [eV]: 3.52e6 - 0
        enddo
        EeV(ntabene) = 0.

        tabe = 0.
        tabi = 0.

        mDT = 0.5 * ( 2.013553212745 + 3.0160492 ) * 1.66053904e-24 ! DT mass [g]
        do ir = 0, ntabrho
           do it = 0, ntabte
              do ie = 0, ntabene - 1

                 ne = rho(ir) / mDT
                 te=TTee(it) / 1.e3      ! Electron temperature  [keV]
                 ti=te
!
! DT plasma with alpha particle projectile
!
                 CALL define_plasma(te,ti,ne,nni)
                 
                 epp=EeV(ie)/1.e3 ! keV

                 IF (epp .EQ. 0) epp=1.E-5
                 
                 write(6,*) ne, te, epp

                 CALL dedx_bps(nni, epp, zp, mp, betab, zb, mb, nb,   &
                      dedx_tot, dedx_i, dedx_e, dedxc_tot, dedxc_i, & 
                      dedxc_e, dedxq_tot, dedxq_i, dedxq_e) ! [MeV/micron]

                 tabe(ir,it,ie) = dedx_e ! MeV/um
                 tabi(ir,it,ie) = dedx_i

                 DEALLOCATE(betab,zb,mb,nb)
                 DEALLOCATE(gb,etab,mpb,mbpb)

                 write(6,*) ne, te, epp, dedx_e, dedx_i

              enddo
              write(6,*) ir,'/',ntabrho, it,'/',ntabte
           enddo
        enddo

!       open(1,err=1,file='tableBPS',status='old',form='formatted')
!       close(1,status='delete')
        open(1,file='tableBPS',status='new',form='formatted')
        write(1,*) ntabrho, ntabte, ntabene
        write(1,*) ( rho(ir) , ir = 0, ntabrho)
        write(1,*) ( TTee(it), it = 0, ntabte )
        write(1,*) ( EeV(ie) , ie = 0, ntabene)
        do ir = 0, ntabrho
           do it = 0, ntabte
              do ie = 0, ntabene
                 write(1,*) - tabe(ir,it,ie)*1.e10, - tabi(ir,it,ie)*1.e10 ! [eV/cm]
              enddo
           enddo
        enddo

        close(1)

      END PROGRAM range
      
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

