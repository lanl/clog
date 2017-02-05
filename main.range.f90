      PROGRAM range
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
        INTEGER NTSMAX
        PARAMETER (NTSMAX=5000)

        REAL    :: y, v, t, dt, tmax
        REAL    :: xt(0:NTSMAX), vt(0:NTSMAX), tt(0:NTSMAX)
        REAL    :: dedx_tot,  dedx_i,  dedx_e
        REAL    :: dedxc_tot, dedxc_i, dedxc_e
        REAL    :: dedxq_tot, dedxq_i, dedxq_e
        INTEGER :: j, nts, nit

        REAL    :: te, ti, ne, ep, mp, zp
        REAL    :: de, epp
        INTEGER :: nni, nn
!
! Range parameters
!
                     ! 3 keV: range.002.out
        tmax=1.0E-11 ! [s]
        nts=100
        nit=10
        nn =300
        dt=tmax/(nts*nit)
        te=3.        ! [keV] Electron temperature       
        ti=3.        ! [keV] Ion temperature            
        ne=1.E25     ! [cm^-3] Electron number density    
        ep=3540.     ! [keV] Projectile energy  
        mp=4*MPKEV   ! [keV] Projectile mass    
        zp=2.        ! [e] Projectile charge  
        OPEN(1, FILE='range.001.out')
        OPEN(2, FILE='dedx.001.out')

        CALL define_plasma_DT(te,ti,ne,nni)

        CALL write_header(ep,mp,zp,te,ti,ne,nni,betab,zb,mb,nb)
!
! stopping power
!
        WRITE(2,*) ne
        WRITE(2,*) te
        WRITE(2,*) ti
        WRITE(2,*) nn
        WRITE(2,*) ep, mp, zp
        de=ep/nn
        epp=0
        DO j=0,nn
           epp=j*de
           IF (epp .EQ. 0) epp=1.E-5
           CALL dedx_bps(nni, epp, zp, mp, betab, zb, mb, nb,   &
             dedx_tot, dedx_i, dedx_e, dedxc_tot, dedxc_i, & 
             dedxc_e, dedxq_tot, dedxq_i, dedxq_e) ! [MeV/micron]
           WRITE (6,'(I6,E17.8,6E22.13)') j, epp, dedx_tot, dedx_i, dedx_e
           WRITE (2,'(I6,E17.8,6E22.13)') j, epp, dedx_tot, dedx_i, dedx_e
        END DO
        CLOSE (2)

! range
!
! initial conditions
!
        t=0                ! [s]
        y=0.               ! [cm]
        v=CC*SQRT(2*ep/mp) ! [cm/s]
        tt(0)=t
        xt(0)=y
        vt(0)=v
!
! evolution
!
        j=0
        ep=0.5*mp*(vt(j)/CC)**2
        CALL dedx_bps(nni, ep, zp, mp, betab, zb, mb, nb,   &
            dedx_tot, dedx_i, dedx_e, dedxc_tot, dedxc_i, & 
            dedxc_e, dedxq_tot, dedxq_i, dedxq_e) ! [MeV/micron]

        WRITE(1,*) ne
        WRITE(1,*) te
        WRITE(1,*) ti
        WRITE(1,*) tmax, nts, nit
        WRITE(1,*) ep, mp, zp
        WRITE (6,'(I6,E17.8,6D22.13)') j, tt(j), xt(j), vt(j), ep, dedx_tot, dedx_i, dedx_e
        WRITE (1,'(I6,E17.8,6D22.13)') j, tt(j), xt(j), vt(j), ep, dedx_tot, dedx_i, dedx_e
        DO j=1,nts
           CALL rk4(y,v,t,dt,nit,nni,zp,mp,dedx_tot,dedx_i,dedx_e)
           tt(j)=t
           xt(j)=y
           vt(j)=v
           ep   =0.5*mp*(vt(j)/CC)**2
           WRITE (6,'(I6,E17.8,6D22.13)') j, tt(j), xt(j), vt(j), ep, dedx_tot, dedx_i, dedx_e
           WRITE (1,'(I6,E17.8,6D22.13)') j, tt(j), xt(j), vt(j), ep, dedx_tot, dedx_i, dedx_e
        END DO
        CLOSE (1)
        END PROGRAM range

!
! SUBROUTINE rk4:
! Given the equation of motion
!
!   y'' = f(y,v,t)  ,
!
! this subroutine advances (y,v) at time t by nit iterations of
! length dt. The variable (y,v) and the time t are returned with
! their updated values. 
!
      SUBROUTINE rk4(x, v, t, dt, nit, nni, zp, mp, dedx_tot, dedx_i, dedx_e)
      USE allocatablevars ! betab(1:nni+1), zb(1:nni+1), mb(1:nni+1), nb(1:nni+1)
      USE physvars
      IMPLICIT NONE
      REAL,    INTENT(INOUT) :: x, v, t
      REAL,    INTENT(INOUT) :: dt
      INTEGER, INTENT(IN)    :: nit
      INTEGER, INTENT(IN)    :: nni
      REAL,    INTENT(IN)    :: zp, mp
      REAL,    INTENT(OUT)   :: dedx_tot, dedx_i, dedx_e ! [MeV/micron]
      REAL    :: xx, vv, th, f
      REAL    :: xk1, xk2, xk3, xk4
      REAL    :: vk1, vk2, vk3, vk4
      REAL    :: dedxc_tot, dedxc_i, dedxc_e, dedxq_tot, dedxq_i, dedxq_e
      REAL    :: ep
      INTEGER :: it
      REAL, PARAMETER :: d6=0.1666666666666666666

      dedx_tot=0
      dedx_i  =0
      dedx_e  =0
      th=t+0.5*dt
      DO it=1,nit
!
! first intermediate step
!
         CALL force(v,f,nni,zp,mp)
         xk1=dt*v
         vk1=dt*f
         xx =x+0.5*xk1
         vv =v+0.5*vk1
!
! second intermediate step
!
         CALL force(vv,f,nni,zp,mp)
         xk2=dt*vv
         vk2=dt*f
         xx =x+0.5*xk2
         vv =v+0.5*vk2
!
! third intermediate step
!
         CALL force(vv,f,nni,zp,mp)
         xk3=dt*vv
         vk3=dt*f
         xx =x+xk3
         vv =v+vk3
!
! fourth intermediate step
!
         t=t+dt
         CALL force(v,f,nni,zp,mp)
         xk4=dt*v
         vk4=dt*f
!
         x = x + d6*(xk1 + 2*xk2 + 2*xk3 + xk4)
         v = v + d6*(vk1 + 2*vk2 + 2*vk3 + vk4)
      END DO
      ep = 0.5*mp*(v/CC)**2 ! [keV] projectile energy
      CALL dedx_bps(nni, ep, zp, mp, betab, zb, mb, nb,   &
        dedx_tot, dedx_i, dedx_e, dedxc_tot, dedxc_i,     & 
        dedxc_e, dedxq_tot, dedxq_i, dedxq_e) ! [MeV/micron]
      END SUBROUTINE rk4
!
! SUBROUTINE force:
! Returns the force f(x,v,t,) for x''=f(x,v,t)
!
      SUBROUTINE force(v, f, nni, zp, mp)
      USE allocatablevars ! betab(1:nni+1), zb(1:nni+1), mb(1:nni+1), nb(1:nni+1)
      USE physvars
      IMPLICIT NONE
      REAL,    INTENT(IN) :: v
      REAL,    INTENT(OUT):: f
      INTEGER, INTENT(IN) :: nni
      REAL,    INTENT(IN) :: zp, mp
      REAL    :: ep
      REAL    :: dedx_tot, dedx_i, dedx_e, dedxc_tot, dedxc_i, dedxc_e 
      REAL    :: dedxq_tot, dedxq_i, dedxq_e
      REAL    :: KEV_TO_ERG=1.6022E-9
      REAL    :: f1, m1

      ep = 0.5*mp*(v/CC)**2  ! [keV] projectile energy
      CALL dedx_bps(nni, ep, zp, mp, betab, zb, mb, nb,   &
            dedx_tot, dedx_i, dedx_e, dedxc_tot, dedxc_i, & 
            dedxc_e, dedxq_tot, dedxq_i, dedxq_e) ! [MeV/micron]

      f   = -dedx_tot             ! [MeV/mircon] is a force]
      f   = f * 1.E7              ! [keV/cm]
      f1  = f * KEV_TO_ERG        ! [erg]
      m1  = mp * KEV_TO_ERG / CC2 ! [g]
      f   = f1 / m1               ! [dyne]
      END SUBROUTINE force

!
! SUBROUTINE define_plasma:
! Returns the plasma species arrays: betab, mb, nb, zb
! Allocates other plasma arrays
!
      SUBROUTINE define_plasma_DT(te, ti, ne, nni)
        
    USE allocatablevars
    USE physvars
      IMPLICIT NONE
      REAL                                :: te     ! [keV]
      REAL                                :: ti     ! [keV]
      REAL                                :: ne     ! [cm^-3]
      INTEGER                             :: nni    ! number of ion species
!     allocatablevars
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
    END SUBROUTINE define_plasma_DT


    SUBROUTINE write_header(ep, mp, zp, te, ti, ne, nni, betab, zb, mb, nb)
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

    END SUBROUTINE write_header
