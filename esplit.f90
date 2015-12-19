! Robert Singleton
! Feb 2008
! This package contains routines for calculating various
! integrated quantities involving the A-coeff, such as the
! energy splitting fraction and the function S(E).
!
! * 20-Oct-2009
!   Not finished. Needs non-equal temperature case. Compare with 
!   Deans Mathematica notebook physics/working/15-Esplit/
!   07-figs_from_dean/G1_Te300_TI300_E26.nb

!==================================================================
! dE/dx calculated numerically from A
!
!  dE           T              d A(E)              1
!  ---  = [1 -  - ] A(E)  - T  ------      with E= - m v^2
!  dx           E                dE                2
!
    SUBROUTINE dedx_a(nni, ep, zp, mp, betab, zb, mb, nb, dedx_tot, &
           dedx_i, dedx_e, dedx_ctot, dedx_ci, dedx_qtot, dedx_qi)
      USE bpsvars
      IMPLICIT NONE
                                                        ! Plasma:
      INTEGER,                     INTENT(IN)  :: nni      !  number of ions
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: betab    !  temp array    [1/keV]
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: mb       !  mass array    [keV]
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: nb       !  density array [1/cc]
      REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: zb       !  charge array
                                                           ! Projectile:
      REAL,                        INTENT(IN)  :: ep       !  projectile engy [keV]
      REAL,                        INTENT(IN)  :: mp       !  projectile mass [keV]
      REAL,                        INTENT(IN)  :: zp       !  projectile charge
                                                           ! Stopping Power
      REAL,                        INTENT(OUT) :: dedx_tot ! [keV/micron]
      REAL,                        INTENT(OUT) :: dedx_i   !
      REAL,                        INTENT(OUT) :: dedx_e   !
      REAL,                        INTENT(OUT) :: dedx_ctot!
      REAL,                        INTENT(OUT) :: dedx_ci  !
      REAL,                        INTENT(OUT) :: dedx_qtot!
      REAL,                        INTENT(OUT) :: dedx_qi  !

      REAL  :: a_tot, a_sumi, a_e, ac_tot, ac_sumi, aq_tot, aq_sumi
      REAL  :: a_tot2, a_sumi2, a_e2, ac_tot2, ac_sumi2, aq_tot2, aq_sumi2
      REAL  :: a_tot1, a_sumi1, a_e1, ac_tot1, ac_sumi1, aq_tot1, aq_sumi1
      REAL  :: dade_tot, dade_i, dade_e, dade_ctot, dade_ci, dade_qtot, dade_qi
      REAL  :: dem, ep2, ep1, te, ti

      te=1./betab(1)  ! [keV]
      ti=1./betab(2)  !
      dem=ep/100.     ! [keV]
      ep2=ep+dem
      CALL coeff_bps(nni,ep2,zp,mp,betab,zb,mb,nb, &
           a_tot2,a_sumi2,a_e2,ac_tot2,ac_sumi2,aq_tot2,aq_sumi2)

      ep1=ep-dem
      CALL coeff_bps(nni,ep1,zp,mp,betab,zb,mb,nb, &
           a_tot1,a_sumi1,a_e1,ac_tot1,ac_sumi1,aq_tot1,aq_sumi1)

      dade_tot   =(a_tot2-a_tot1)/(2*dem)            ! derivative of A at ep
      dade_i     =(a_sumi2-a_sumi1)/(2*dem)
      dade_e     =(a_e2-a_e1)/(2*dem)
      dade_ctot  =(ac_tot2-ac_tot1)/(2*dem)
      dade_ci    =(ac_sumi2-ac_sumi1)/(2*dem)
      dade_qtot  =(aq_tot2-aq_tot1)/(2*dem)
      dade_qi    =(aq_sumi2-aq_sumi1)/(2*dem)

      CALL coeff_bps(nni,ep,zp,mp,betab,zb,mb,nb, &
           a_tot,a_sumi,a_e,ac_tot,ac_sumi,aq_tot,aq_sumi)

      dedx_tot   =(1.-te/ep)*a_tot - te*dade_tot   ! stopping power at ep
      dedx_i     =(1.-te/ep)*a_sumi - te*dade_i
      dedx_e     =(1.-te/ep)*a_e - te*dade_e
      dedx_ctot  =(1.-te/ep)*ac_tot - te*dade_ctot
      dedx_ci    =(1.-te/ep)*ac_sumi - te*dade_ci
      dedx_qtot  =(1.-te/ep)*aq_tot - te*dade_qtot
      dedx_qi    =(1.-te/ep)*aq_sumi - te*dade_qi

    END SUBROUTINE dedx_a


!==================================================================
! AA1(E)= A_e(E)A_i(E)/[Te Ae(E) + Ti Ai(E)]
! 
!
SUBROUTINE aa1_int(nni, ep, zp, mp, betab, zb, mb, nb, aa1)
  USE bpsvars
  IMPLICIT NONE
                                                      ! Plasma:
   INTEGER,                     INTENT(IN)  :: nni    !  number of ions
   REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: betab  !  temp array    [1/keV]
   REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: mb     !  mass array    [keV]
   REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: nb     !  density array [1/cc]
   REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: zb     !  charge array
                                                      ! Projectile:
   REAL,                        INTENT(IN)  :: ep     !  projectile energy [keV]
   REAL,                        INTENT(IN)  :: mp     !  projectile mass   [keV]
   REAL,                        INTENT(IN)  :: zp     !  projectile charge
                                                      ! Integrand of S(E):
   REAL,                        INTENT(OUT) :: aa1    !  [1/keV]

   REAL  :: a_tot, a_i, a_e, ac_tot, ac_i, aq_tot, aq_i
   REAL  :: te, ti

   CALL coeff_bps(nni,ep,zp,mp,betab,zb,mb,nb,a_tot,a_i,a_e,&
     ac_tot,ac_i,aq_tot,aq_i)   ! A-coeff keV/micron

   te=1./betab(1)               ! elec temp [keV]
   ti=1./betab(2)               ! ion  temp [keV]
   aa1=a_e*a_i/(te*a_e + ti*a_i)! [A/T]=[keV/micron-keV]=[1/micron]

END SUBROUTINE aa1_int

!==================================================================
!
!        /E
!        |         A_e(E') + A_I(E')
! S(E) = | dE' -------------------------   by Gaussian quadrature.
!        |      te*A_e(E') + ti*A_I(E')
!       /0
!
 SUBROUTINE s_int(nni, zp, mp, betab, zb, mb, nb, ng, e1, e2, s)
  USE bpsvars
  IMPLICIT NONE
                                                     ! Plasma:
   INTEGER,                     INTENT(IN) :: nni    !  number of ions
   REAL,    DIMENSION(1:nni+1), INTENT(IN) :: betab  !  temp array    [1/keV]
   REAL,    DIMENSION(1:nni+1), INTENT(IN) :: mb     !  mass array    [keV]
   REAL,    DIMENSION(1:nni+1), INTENT(IN) :: nb     !  density array [1/cc]
   REAL,    DIMENSION(1:nni+1), INTENT(IN) :: zb     !  charge array
                                                     ! Projectile:
   REAL,                        INTENT(IN) :: mp     !  projectile mass   [keV]
   REAL,                        INTENT(IN) :: zp     !  projectile charge
                                                     !Integration 
   REAL,                        INTENT(IN) :: e1     ! lower limit [keV]
   REAL,                        INTENT(IN) :: e2     ! upper limit [keV]
   INTEGER,                     INTENT(IN) :: ng     ! number of sub-intervals
   REAL,                        INTENT(OUT):: s      ! [1]

   REAL,    PARAMETER :: XPM=0.7745966692E0
   REAL,    PARAMETER :: W13=0.5555555556E0, W2=0.8888888889E0
   REAL               :: ds, em, ep, de
   INTEGER            :: iu
   
   IF (MOD(ng,2) .NE. 0) THEN
      PRINT *, 'ERROR in s_int: the number of sub-intervals in s_int must be even.'
      STOP
   ENDIF      

   IF (e2 < 1.E-5) THEN
     s=0
   ELSE
     s=0
     ds=0
     de=(e2-e1)/ng
     ep=e1-de
     DO iu=1,ng,2
!
        ep=ep + 2.*de
        CALL ds_int(nni,ep,zp,mp,betab,zb,mb,nb,ds)
        s=s + W2*ds
!
        em=ep - de*XPM
        CALL ds_int(nni,em,zp,mp,betab,zb,mb,nb,ds)
        s=s + W13*ds
!
        em=ep + de*XPM
        CALL ds_int(nni,em,zp,mp,betab,zb,mb,nb,ds)
        s=s + W13*ds
     ENDDO
     s=s*de
   ENDIF
 END SUBROUTINE s_int
!
! integrand dS=A(E')/[Te Ae(E') + Ti Ai(E')]
!
SUBROUTINE ds_int(nni, ep, zp, mp, betab, zb, mb, nb, ds)
  USE bpsvars
  IMPLICIT NONE
                                                      ! Plasma:
   INTEGER,                     INTENT(IN)  :: nni    !  number of ions
   REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: betab  !  temp array    [1/keV]
   REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: mb     !  mass array    [keV]
   REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: nb     !  density array [1/cc]
   REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: zb     !  charge array
                                                      ! Projectile:
   REAL,                        INTENT(IN)  :: ep     !  projectile energy [keV]
   REAL,                        INTENT(IN)  :: mp     !  projectile mass   [keV]
   REAL,                        INTENT(IN)  :: zp     !  projectile charge
                                                      ! Integrand of S(E):
   REAL,                        INTENT(OUT) :: ds     !  [1/keV]


   REAL  :: a_tot, a_i, a_e, ac_tot, ac_i, aq_tot, aq_i
   REAL  :: te, ti

   CALL coeff_bps(nni,ep,zp,mp,betab,zb,mb,nb,a_tot,a_i,a_e,&
     ac_tot,ac_i,aq_tot,aq_i) ! A-coeff keV/micron

   te=1./betab(1)             ! elec temp [keV]
   ti=1./betab(2)             ! ion  temp [keV]
   ds=a_tot/(te*a_e + ti*a_i) ! [1/keV]
END SUBROUTINE ds_int

!==================================================================
!
!
!            /infty                            
!            |                                 
! I(Te,Ti) = | dE'  Sqrt(E') exp[-S(E';te,ti)] 
!            |                                 
!           /0                                 

 SUBROUTINE i_int(nni, zp, mp, betab, zb, mb, nb, ng, e1, e2, int, ionet)
  USE bpsvars
  USE mathvars
  IMPLICIT NONE
                                                     ! Plasma:
   INTEGER,                     INTENT(IN) :: nni    !  number of ions
   REAL,    DIMENSION(1:nni+1), INTENT(IN) :: betab  !  temp array    [1/keV]
   REAL,    DIMENSION(1:nni+1), INTENT(IN) :: mb     !  mass array    [keV]
   REAL,    DIMENSION(1:nni+1), INTENT(IN) :: nb     !  density array [1/cc]
   REAL,    DIMENSION(1:nni+1), INTENT(IN) :: zb     !  charge array
                                                     ! Projectile:
   REAL,                        INTENT(IN) :: mp     !  projectile mass   [keV]
   REAL,                        INTENT(IN) :: zp     !  projectile charge
                                                     !Integration 
   REAL,                        INTENT(IN) :: e1     ! lower limit [keV]
   REAL,                        INTENT(IN) :: e2     ! upper limit [keV]
   INTEGER,                     INTENT(IN) :: ng     ! number of sub-intervals
   REAL,                        INTENT(OUT):: int    ! [fm^3]
   REAL,                        INTENT(OUT):: ionet  ! [fm^3]

   REAL,    PARAMETER :: XPM=0.7745966692E0
   REAL,    PARAMETER :: W13=0.5555555556E0, W2=0.8888888889E0
   REAL               :: dint, ym, yp, dy, y2, y1, ep, em
   REAL               :: e0, s, ti
   INTEGER            :: iu, ngs

      IF (MOD(NG,2) .NE. 0) THEN
         PRINT *, 'ERROR in n_int: the number of sub-intervals in s_int must be even.'
         STOP
      ENDIF      
      int=0
      y2=SQRT(e2)
      y1=SQRT(e1)
      dy=(y2-y1)/NG
      yp=y1-dy
      e0=0.

    ti=1./betab(2)
    ionet=0.5*SQRT(PI)*ti**1.5

!     /infty                               /infty
!     |                                    |
! I = | dE'  Sqrt(E') exp[-S(E';te,ti)] = 2| y dy * y Exp[-S(y)]
!     |                                    |
!    /0                                   /0

   DO iu=1,NG,2
      yp=yp + 2.*dy
      ngs=12

      ep=yp*yp
      CALL s_int(nni,zp,mp,betab,zb,mb,nb,ngs,e0,ep,s) 
!s=yp*yp/ti
      dint = 2*yp*yp*EXP(-s)
      int=int + W2*dint

      ym=yp - dy*XPM
      em=ym*ym
      CALL s_int(nni,zp,mp,betab,zb,mb,nb,ngs,e0,em,s) 
!s=ym*ym/ti
      dint = 2*ym*ym*EXP(-s)
      int=int + W13*dint

      ym=yp + dy*XPM
      em=ym*ym
      CALL s_int(nni,zp,mp,betab,zb,mb,nb,ngs,e0,em,s) 
!s=ym*ym/ti
      dint = 2*ym*ym*EXP(-s)
      int=int + W13*dint
   ENDDO
   int=int*dy
 END SUBROUTINE i_int

!==================================================================
!
! N(Te,Ti)
!
 SUBROUTINE n_int(nni, zp, mp, betab, zb, mb, nb, ng, e1, e2, n, naprx)
  USE bpsvars
  USE mathvars
  IMPLICIT NONE
                                                     ! Plasma:
   INTEGER,                     INTENT(IN) :: nni    !  number of ions
   REAL,    DIMENSION(1:nni+1), INTENT(IN) :: betab  !  temp array    [1/keV]
   REAL,    DIMENSION(1:nni+1), INTENT(IN) :: mb     !  mass array    [keV]
   REAL,    DIMENSION(1:nni+1), INTENT(IN) :: nb     !  density array [1/cc]
   REAL,    DIMENSION(1:nni+1), INTENT(IN) :: zb     !  charge array
                                                     ! Projectile:
   REAL,                        INTENT(IN) :: mp     !  projectile mass   [keV]
   REAL,                        INTENT(IN) :: zp     !  projectile charge
                                                     !Integration 
   REAL,                        INTENT(IN) :: e1     ! lower limit [keV]
   REAL,                        INTENT(IN) :: e2     ! upper limit [keV]
   INTEGER,                     INTENT(IN) :: ng     ! number of sub-intervals
   REAL,                        INTENT(OUT):: n      ! [fm^3]
   REAL,                        INTENT(OUT):: naprx  ! [fm^3]

   REAL,    PARAMETER :: NCOEFF=1.0725E-10
   REAL,    PARAMETER :: XPM=0.7745966692E0
   REAL,    PARAMETER :: W13=0.5555555556E0, W2=0.8888888889E0
   REAL               :: dn, ym, yp, dy, y2, y1, ep, em
   REAL               :: e0, s, ti, mpp
   INTEGER            :: iu, ngs

   IF (MOD(NG,2) .NE. 0) THEN
      PRINT *, 'ERROR in n_int: the number of sub-intervals in s_int must be even.'
      STOP
   ENDIF      
   n=0
   dn=0
   y2=SQRT(e2)
   y1=SQRT(e1)
   dy=(y2-y1)/NG
   yp=y1-dy
   e0=0.
!
! First do the integral 
!
!
!     /infty                               /infty
!     |                                    |
! I = | dE'  Sqrt(E') exp[-S(E';te,ti)] = 2| dy y^2 Exp[-S(y)]
!     |                                    |
!    /0                                   /0


   DO iu=1,NG,2
      yp=yp + 2.*dy
      ngs=12

      ep=yp*yp
      CALL s_int(nni,zp,mp,betab,zb,mb,nb,ngs,e0,ep,s) 
      dn = 2*yp*yp*EXP(-s)
      n=n + W2*dn

      ym=yp - dy*XPM
      em=ym*ym
      CALL s_int(nni,zp,mp,betab,zb,mb,nb,ngs,e0,em,s) 
      dn = 2*ym*ym*EXP(-s)
      n=n + W13*dn

      ym=yp + dy*XPM
      em=ym*ym
      CALL s_int(nni,zp,mp,betab,zb,mb,nb,ngs,e0,em,s) 
      dn = 2*ym*ym*EXP(-s)
      n=n + W13*dn
   ENDDO
   n=n*dy
   mpp=mp**1.5
   n=NCOEFF/mpp/n
!
   ti=1./betab(2)
   naprx=0.5*SQRT(PI)*ti**1.5
   naprx=NCOEFF/mpp/naprx  ! or  naprx=1.2102E-10/mpp/ti**1.5
 END SUBROUTINE n_int



!==================================================================
!
! H(E;Te,Ti)
!
 SUBROUTINE h_int(nni, zp, mp, betab, zb, mb, nb, ng, e1, ep, h, honet)
  USE bpsvars
  IMPLICIT NONE
                                                     ! Plasma:
   INTEGER,                     INTENT(IN) :: nni    !  number of ions
   REAL,    DIMENSION(1:nni+1), INTENT(IN) :: betab  !  temp array    [1/keV]
   REAL,    DIMENSION(1:nni+1), INTENT(IN) :: mb     !  mass array    [keV]
   REAL,    DIMENSION(1:nni+1), INTENT(IN) :: nb     !  density array [1/cc]
   REAL,    DIMENSION(1:nni+1), INTENT(IN) :: zb     !  charge array
                                                     ! Projectile:
   REAL,                        INTENT(IN) :: ep     !  projectile energy [keV]
   REAL,                        INTENT(IN) :: mp     !  projectile mass   [keV]
   REAL,                        INTENT(IN) :: zp     !  projectile charge
                                                     !Integration 
   REAL,                        INTENT(IN) :: e1     ! lower limit [keV]
   INTEGER,                     INTENT(IN) :: ng     ! number of sub-intervals
   REAL,                        INTENT(OUT):: h      ! [fm^3]
   REAL,                        INTENT(OUT):: honet  ! [fm^3]

   REAL,    PARAMETER :: XPM=0.7745966692E0
   REAL,    PARAMETER :: W13=0.5555555556E0, W2=0.8888888889E0
   REAL               :: dh, ym, yp, dy, y2, y1, em1, ep1, emax
   REAL               :: e0, s, es1, int, ionet
   INTEGER            :: iu, ngs
!
! Return asymptotic value of 1 when ep > the I-integral's large upper limit.
!
   emax=8./betab(2)     ! 8*ti [keV]
   IF (ep >= emax) THEN
     h=1.
     honet=1. ! not correct in this case, but who cares
     RETURN
   ENDIF

!
! Perform integrals
!
   IF (MOD(NG,2) .NE. 0) THEN
      PRINT *, 'ERROR in n_int: the number of sub-intervals in s_int must be even.'
      STOP
   ENDIF      
   h=0
   y2=SQRT(ep)
   y1=SQRT(e1)
   dy=(y2-y1)/NG
   yp=y1-dy

!
! First calculate I=I(Te,Ti) for normalization
!
    ngs=12
    es1=0
    CALL i_int(nni, zp, mp, betab, zb, mb, nb, ngs, es1, emax, int, ionet)

!                 /E                                       /y
!              1  |                                    1   |
!H(E,Te,Ti) = --- | dE'  Sqrt(E') exp[-S(E';te,ti)] = --- 2| dy y^2 Exp[-S(y)]
!              I  |                                    I   |
!                /0                                       /0

   e0=0.
   DO iu=1,NG,2
      yp=yp + 2.*dy
      ngs=12

      ep1=yp*yp
      CALL s_int(nni,zp,mp,betab,zb,mb,nb,ngs,e0,ep1,s) 
!s=yp*yp/ti
      dh = 2*yp*yp*EXP(-s)
      h=h + W2*dh

      ym=yp - dy*XPM
      em1=ym*ym
      CALL s_int(nni,zp,mp,betab,zb,mb,nb,ngs,e0,em1,s) 
!s=ym*ym/ti
      dh = 2*ym*ym*EXP(-s)
      h=h + W13*dh

      ym=yp + dy*XPM
      em1=ym*ym
      CALL s_int(nni,zp,mp,betab,zb,mb,nb,ngs,e0,em1,s) 
!s=ym*ym/ti
      dh = 2*ym*ym*EXP(-s)
      h=h + W13*dh
   ENDDO
   h=h*dy
   h=h/int
   honet=h/ionet
 END SUBROUTINE h_int


!==================================================================
!
! I1(Te,Ti) cf. I1approx(Te,Ti)
!
 SUBROUTINE i1_int(nni, zp, mp, betab, zb, mb, nb, ng, e1, ep, i1, i1onet)
  USE bpsvars
  IMPLICIT NONE
                                                     ! Plasma:
   INTEGER,                     INTENT(IN) :: nni    !  number of ions
   REAL,    DIMENSION(1:nni+1), INTENT(IN) :: betab  !  temp array    [1/keV]
   REAL,    DIMENSION(1:nni+1), INTENT(IN) :: mb     !  mass array    [keV]
   REAL,    DIMENSION(1:nni+1), INTENT(IN) :: nb     !  density array [1/cc]
   REAL,    DIMENSION(1:nni+1), INTENT(IN) :: zb     !  charge array
                                                     ! Projectile:
   REAL,                        INTENT(IN) :: ep     !  projectile energy [keV]
   REAL,                        INTENT(IN) :: mp     !  projectile mass   [keV]
   REAL,                        INTENT(IN) :: zp     !  projectile charge
                                                     !Integration 
   REAL,                        INTENT(IN) :: e1     ! lower limit [keV]
   INTEGER,                     INTENT(IN) :: ng     ! number of sub-intervals
   REAL,                        INTENT(OUT):: i1     ! [dimensionless]
   REAL,                        INTENT(OUT):: i1onet ! [dimensionless]

   REAL,    PARAMETER :: COEFF=1.2092  ! 2 Pi/[3 Sqrt(3)]
   REAL,    PARAMETER :: XPM=0.7745966692E0
   REAL,    PARAMETER :: W13=0.5555555556E0, W2=0.8888888889E0
   REAL               :: ym, yp, dy, y2, y1, em1, ep1, di1
   REAL               :: e0, ec, te, ti
   INTEGER            :: iu, ngs
   

!
! Perform integrals
!
   IF (MOD(ng,2) .NE. 0) THEN
      PRINT *, 'ERROR in n_int: the number of sub-intervals in s_int must be even.'
      STOP
   ENDIF      
   i1=0
   y2=SQRT(ep)
   y1=SQRT(e1)
   dy=(y2-y1)/ng
   yp=y1-dy

!            /E0                             /y0
!            |      Ai(E) Ae(E)   H(E)       |          Ai(y) Ae(y)   H(y)       
! I1Te,Ti) = | dE   -----------  -----  =  2 | y dy *   -----------  -----
!            |      < T A(E) >    A(E)       |          < T A(y) >    A(y)     
!           /0                               /0

   e0=0.
   DO iu=1,ng,2
      yp=yp + 2.*dy
      ngs=6

      ep1=yp*yp
      CALL di1_int(nni, ep1, zp, mp, betab, zb, mb, nb, di1)
      i1=i1 + W2*2*yp*di1

      ym=yp - dy*XPM
      em1=ym*ym
      CALL di1_int(nni, em1, zp, mp, betab, zb, mb, nb, di1)
      i1=i1 + W13*2*ym*di1

      ym=yp + dy*XPM
      em1=ym*ym
      CALL di1_int(nni, em1, zp, mp, betab, zb, mb, nb, di1)
      i1=i1 + W13*2*ym*di1
   ENDDO
   i1=i1*dy
!
! I1oneT (analytic approx)
!
  te=1./betab(1)
  ti=1./betab(2)
  ec=500.        ! Ec=500 keV for Te=Ti=10 keV
  IF (te .eq. ti) THEN
    i1onet=(2*ec/te)*COEFF*0.6666
  ELSE
    i1onet=(2*ec/te)*COEFF*te**0.3333*(ti**0.6666-te**0.6666)/(ti-te)
  ENDIF
  i1onet=i1onet - (2*ec/te)*SQRT(ec/3500)  ! large-E correction

! Now multiply by E_c. The velocity v_c is given by (B62)
!    
!
!   ne=nb(1)
!   ge=(6.1260E-15)*SQRT(ne)/te**1.5
!   eta=(2.7212E-2)/SQRT(te)

 END SUBROUTINE i1_int


!==================================================================
!
!                                  Ai(E) Ae(E)   H(E)
! The integrand of dI1(E;Te,Ti) =  -----------  -----
!                                  < T A(E) >    A(E)
!
 SUBROUTINE di1_int(nni, ep, zp, mp, betab, zb, mb, nb, di1)
  USE bpsvars
  IMPLICIT NONE
                                                     ! Plasma:
   INTEGER,                     INTENT(IN) :: nni    !  number of ions
   REAL,    DIMENSION(1:nni+1), INTENT(IN) :: betab  !  temp array    [1/keV]
   REAL,    DIMENSION(1:nni+1), INTENT(IN) :: mb     !  mass array    [keV]
   REAL,    DIMENSION(1:nni+1), INTENT(IN) :: nb     !  density array [1/cc]
   REAL,    DIMENSION(1:nni+1), INTENT(IN) :: zb     !  charge array
                                                     ! Projectile:
   REAL,                        INTENT(IN) :: ep     !  projectile energy [keV]
   REAL,                        INTENT(IN) :: mp     !  projectile mass   [keV]
   REAL,                        INTENT(IN) :: zp     !  projectile charge
                                                     !Integrand of I1
   REAL,                        INTENT(OUT):: di1    ! [dimensionless]

   REAL               :: e1, h, honet, te, ti, tae
   REAL               :: a_tot, a_sumi, a_e, ac_tot, ac_sumi, aq_tot, aq_sumi
   INTEGER            :: ngs

!
! Calculate H=H(E,Te,Ti) for E=ep
!
    ngs=6
    e1=0
    h=1.
    CALL h_int(nni, zp, mp, betab, zb, mb, nb, ngs, e1, ep, h, honet)
!
! calculate corresponding Ai and Ae coefficients
!
    te=1./betab(1) ! electron temp [keV]
    ti=1./betab(2) ! ion temp in   [keV]
    CALL coeff_bps(nni, ep, zp, mp, betab, zb, mb, nb,      &
     a_tot, a_sumi, a_e, ac_tot, ac_sumi, aq_tot, aq_sumi)
    tae = te*a_e + ti*a_sumi 
    di1=a_e*a_sumi*h/tae/a_tot
 END SUBROUTINE di1_int


!==================================================================
!
!                                        Ai(E) Ae(E)
! The integrand of rate2.35 = dI2 = E * -----------  Exp[-S(E;Te,Ti)]
!                                        < T A(E) > 
!
 SUBROUTINE di2_int(nni, ep, zp, mp, betab, zb, mb, nb, di2)
  USE bpsvars
  IMPLICIT NONE
                                                     ! Plasma:
   INTEGER,                     INTENT(IN) :: nni    !  number of ions
   REAL,    DIMENSION(1:nni+1), INTENT(IN) :: betab  !  temp array    [1/keV]
   REAL,    DIMENSION(1:nni+1), INTENT(IN) :: mb     !  mass array    [keV]
   REAL,    DIMENSION(1:nni+1), INTENT(IN) :: nb     !  density array [1/cc]
   REAL,    DIMENSION(1:nni+1), INTENT(IN) :: zb     !  charge array
                                                     ! Projectile:
   REAL,                        INTENT(IN) :: ep     !  projectile energy [keV]
   REAL,                        INTENT(IN) :: mp     !  projectile mass   [keV]
   REAL,                        INTENT(IN) :: zp     !  projectile charge
                                                     !Integrand of dI2
   REAL,                        INTENT(OUT):: di2    ! [keV/micron]

   REAL               :: e2, e1, s, te, ti, tae
   REAL               :: a_tot, a_i, a_e, ac_tot, ac_i, aq_tot, aq_i
   INTEGER            :: ngs

!
! Calculate S=S(E,Te,Ti) for E=ep
!
    ngs=12
    e1=0
    e2=ep
    s=0.
    CALL s_int(nni, zp, mp, betab, zb, mb, nb, ngs, e1, e2, s) 
!
! Calculate corresponding Ai and Ae coefficients
!
    te=1./betab(1) ! electron temp [keV]
    ti=1./betab(2) ! ion temp in   [keV]
    CALL coeff_bps(nni, ep, zp, mp, betab, zb, mb, nb,      &
     a_tot, a_i, a_e, ac_tot, ac_i, aq_tot, aq_i)
    tae = te*a_e + ti*a_i 
    di2=ep*a_i*a_e*EXP(-s)/tae
 END SUBROUTINE di2_int

!==================================================================
!
!             /infty                               /infty
!             |       Ai(E) Ae(E)                  |            Ai(y) Ae(y)  
! I2(Te,Ti) = | dE E  -----------  exp[-S(E)] =  2 | y dy E(y) -----------  exp[-S(y)]
!             |       < T A(E) >                   |            < T A(y) >   
!            /0                                   /0

      SUBROUTINE i2_int(nni, zp, mp, betab, zb, mb, nb, ng, e1, e2, i2)
      USE bpsvars
      IMPLICIT NONE
                                                        ! Plasma:
      INTEGER,                     INTENT(IN) :: nni    !  number of ions
      REAL,    DIMENSION(1:nni+1), INTENT(IN) :: betab  !  temp array    [1/keV]
      REAL,    DIMENSION(1:nni+1), INTENT(IN) :: mb     !  mass array    [keV]
      REAL,    DIMENSION(1:nni+1), INTENT(IN) :: nb     !  density array [1/cc]
      REAL,    DIMENSION(1:nni+1), INTENT(IN) :: zb     !  charge array
                                                        ! Projectile:
      REAL,                        INTENT(IN) :: mp     !  projectile mass   [keV]
      REAL,                        INTENT(IN) :: zp     !  projectile charge
                                                        !Integration 
      REAL,                        INTENT(IN) :: e1     ! lower limit [keV]
      REAL,                        INTENT(IN) :: e2     ! upper limit [keV]
      INTEGER,                     INTENT(IN) :: ng     ! number of sub-intervals
      REAL,                        INTENT(OUT):: i2     ! [keV^2/micron]

      REAL,    PARAMETER :: XPM=0.7745966692E0
      REAL,    PARAMETER :: W13=0.5555555556E0, W2=0.8888888889E0
      REAL               :: ym, yp, dy, y2, y1, em1, ep1, di2
      REAL               :: e0
      INTEGER            :: iu, ngs
   
      IF (MOD(ng,2) .NE. 0) THEN
         PRINT *, 'ERROR in n_int: the number of sub-intervals in s_int must be even.'
         STOP
      ENDIF      
      i2=0
      y2=SQRT(e2)
      y1=SQRT(e1)
      dy=(y2-y1)/ng
      yp=y1-dy

      e0=0.
      DO iu=1,ng,2
        yp=yp + 2.*dy
        ngs=6

        ep1=yp*yp
        CALL di2_int(nni, ep1, zp, mp, betab, zb, mb, nb, di2) 
        i2=i2 + W2*2*yp*di2

        ym=yp - dy*XPM
        em1=ym*ym
        CALL di2_int(nni, em1, zp, mp, betab, zb, mb, nb, di2)
        i2=i2 + W13*2*ym*di2

        ym=yp + dy*XPM
        em1=ym*ym
        CALL di2_int(nni, em1, zp, mp, betab, zb, mb, nb, di2)
        i2=i2 + W13*2*ym*di2
      ENDDO
      i2=i2*dy ! i2 [keV^2/micron]

      END SUBROUTINE i2_int

      SUBROUTINE rate_i2_int(nni, zp, mp, betab, zb, mb, nb, ng, e1, e2, rate_i2)
      USE bpsvars
      USE physvars
      IMPLICIT NONE
                                                        ! Plasma:
      INTEGER,                     INTENT(IN) :: nni    !  number of ions
      REAL,    DIMENSION(1:nni+1), INTENT(IN) :: betab  !  temp array    [1/keV]
      REAL,    DIMENSION(1:nni+1), INTENT(IN) :: mb     !  mass array    [keV]
      REAL,    DIMENSION(1:nni+1), INTENT(IN) :: nb     !  density array [1/cc]
      REAL,    DIMENSION(1:nni+1), INTENT(IN) :: zb     !  charge array
                                                        ! Projectile:
      REAL,                        INTENT(IN) :: mp     !  projectile mass   [keV]
      REAL,                        INTENT(IN) :: zp     !  projectile charge
                                                        !Integration 
      REAL,                        INTENT(IN) :: e1     ! lower limit [keV]
      REAL,                        INTENT(IN) :: e2     ! upper limit [keV]
      INTEGER,                     INTENT(IN) :: ng     ! number of sub-intervals
      REAL,                        INTENT(OUT):: rate_i2! [1/s]
      REAL,    PARAMETER :: CFACT=2.8929E9

      REAL               :: i2, a
      REAL               :: int, ionet

      rate_i2=0.
      CALL i2_int(nni, zp, mp, betab, zb, mb, nb, ng, e1, e2, i2)
      CALL i_int(nni, zp, mp, betab, zb, mb, nb, ng, e1, e2, int, ionet)
      a=mp/AMUKEV
      a=SQRT(a)
      rate_i2=i2*CFACT/int/a ! s^-1
      END SUBROUTINE rate_i2_int


!==================================================================
! 
! fractional electron-ion energy splitting
! 
!
    SUBROUTINE eifrac_int_1t(nni, ep, e0, ec, zp, mp, betab, zb, mb, nb, fe, fi)
      USE bpsvars
      IMPLICIT NONE
                                                        ! Plasma:
      INTEGER,                     INTENT(IN) :: nni    !  number of ions
      REAL,    DIMENSION(1:nni+1), INTENT(IN) :: betab  !  temp array    [1/keV]
      REAL,    DIMENSION(1:nni+1), INTENT(IN) :: mb     !  mass array    [keV]
      REAL,    DIMENSION(1:nni+1), INTENT(IN) :: nb     !  density array [1/cc]
      REAL,    DIMENSION(1:nni+1), INTENT(IN) :: zb     !  charge array
                                                        ! Projectile:
      REAL,                        INTENT(IN) :: ep     !  initial projectile energy [keV]
      REAL,                        INTENT(IN) :: e0     !  threshold energy  [keV];usually ep=e0 
      REAL,                        INTENT(IN) :: ec     !  cross over energy [keV];A_e=A_I
      REAL,                        INTENT(IN) :: mp     !  projectile mass   [keV]
      REAL,                        INTENT(IN) :: zp     !  projectile charge
                                                        !
      REAL,                        INTENT(OUT):: fe     ! 
      REAL,                        INTENT(OUT):: fi     ! 

      REAL               :: te, ti, f 
      REAL               :: dfe, dfi, de, e1, e2, em, epp
      INTEGER            :: ng, iu

      REAL,    PARAMETER :: XPM=0.7745966692E0
      REAL,    PARAMETER :: W13=0.5555555556E0, W2=0.8888888889E0
      te=1./betab(1) ! electron temp [keV]
      ti=1./betab(2) ! ion temp in   [keV]
!*!
      f=ec ! place holder
!*!
      fe=0
      fi=0
      e1=0.
      e2=ep
      ng=10
      de=(e2-e1)/ng
      epp=e1-de
      DO iu=1,ng,2
        epp=epp + 2.*de

        em=epp
        CALL dfei_int_1t(nni,em,zp,mp,betab,zb,mb,nb,dfe,dfi) 
        fe=fe + W2*dfe
        fi=fi + W2*dfi

        em=epp - de*XPM
        CALL dfei_int_1t(nni,em,zp,mp,betab,zb,mb,nb,dfe,dfi) 
        fe=fe + W13*dfe
        fi=fi + W13*dfi

        em=epp + de*XPM
        CALL dfei_int_1t(nni,em,zp,mp,betab,zb,mb,nb,dfe,dfi) 
        fe=fe + W13*dfe
        fi=fi + W13*dfi
      ENDDO
      fe=fe*de/e0
      fi=fi*de/e0
    END SUBROUTINE eifrac_int_1t
    SUBROUTINE dfei_int_1t(nni,ep,zp,mp,betab,zb,mb,nb,dfe,dfi) 
      USE bpsvars
      USE mathvars
      IMPLICIT NONE
                                                        ! Plasma:
      INTEGER,                     INTENT(IN) :: nni    !  number of ions
      REAL,    DIMENSION(1:nni+1), INTENT(IN) :: betab  !  temp array    [1/keV]
      REAL,    DIMENSION(1:nni+1), INTENT(IN) :: mb     !  mass array    [keV]
      REAL,    DIMENSION(1:nni+1), INTENT(IN) :: nb     !  density array [1/cc]
      REAL,    DIMENSION(1:nni+1), INTENT(IN) :: zb     !  charge array
                                                        ! Projectile:
      REAL,                        INTENT(IN) :: ep     !  projectile energy [keV]
      REAL,                        INTENT(IN) :: mp     !  projectile mass   [keV]
      REAL,                        INTENT(IN) :: zp     !  projectile charge
                                                        !
      REAL,                        INTENT(OUT):: dfe    ! 
      REAL,                        INTENT(OUT):: dfi    ! 

      REAL               :: te, ti, a_tot, a_i, a_e, ac_tot, ac_i, aq_tot, aq_i
      REAL               :: tae, x1, x2, f, ferf

      te=1./betab(1) ! electron temp [keV]
      ti=1./betab(2) ! ion temp in   [keV]
      dfe=0.
      dfi=0.

      CALL coeff_bps(nni,ep,zp,mp,betab,zb,mb,nb, &
           a_tot,a_i,a_e,ac_tot,ac_i,aq_tot,aq_i)
      tae=te*a_e + ti*a_i
      dfe=te*a_e/tae
      dfi=ti*a_i/tae

      dfe=te*a_e/tae
      dfi=ti*a_i/tae

      x1=ep/ti
      x2=SQRT(x1)
      f=ferf(x2) - 2*x2*EXP(-x1)/SQPI
      dfe=dfe*f
      dfi=dfi*f
  END SUBROUTINE dfei_int_1t

!
! generalization to multiple temperatures
!
 SUBROUTINE eifrac_int(nni, ep, e0, ec, zp, mp, betab, zb, mb, nb, fe, fi)
     USE bpsvars
     IMPLICIT NONE
                                                        ! Plasma:
      INTEGER,                     INTENT(IN) :: nni    !  number of ions
      REAL,    DIMENSION(1:nni+1), INTENT(IN) :: betab  !  temp array    [1/keV]
      REAL,    DIMENSION(1:nni+1), INTENT(IN) :: mb     !  mass array    [keV]
      REAL,    DIMENSION(1:nni+1), INTENT(IN) :: nb     !  density array [1/cc]
      REAL,    DIMENSION(1:nni+1), INTENT(IN) :: zb     !  charge array
                                                        ! Projectile:
      REAL,                        INTENT(IN) :: ep     !  initial projectile energy [keV]
      REAL,                        INTENT(IN) :: e0     !  threshold energy  [keV];usually ep=e0 
      REAL,                        INTENT(IN) :: ec     !  cross over energy [keV];A_e=A_I
      REAL,                        INTENT(IN) :: mp     !  projectile mass   [keV]
      REAL,                        INTENT(IN) :: zp     !  projectile charge
                                                        !
      REAL,                        INTENT(OUT):: fe     ! 
      REAL,                        INTENT(OUT):: fi     ! 

      REAL               :: te, ti, f
      REAL               :: dfe, dfi, de, e1, e2, em, epp
      INTEGER            :: ng, iu

      REAL,    PARAMETER :: XPM=0.7745966692E0
      REAL,    PARAMETER :: W13=0.5555555556E0, W2=0.8888888889E0
      te=1./betab(1) ! electron temp [keV]
      ti=1./betab(2) ! ion temp in   [keV]

      fe=0
      fi=0
      e1=0.
      e2=ep
      ng=10
      de=(e2-e1)/ng
      epp=e1-de
      DO iu=1,ng,2
        epp=epp + 2.*de

        em=epp
        CALL dfei_int(nni,em,zp,mp,betab,zb,mb,nb,dfe,dfi) 
        fe=fe + W2*dfe
        fi=fi + W2*dfi

        em=epp - de*XPM
        CALL dfei_int(nni,em,zp,mp,betab,zb,mb,nb,dfe,dfi) 
        fe=fe + W13*dfe
        fi=fi + W13*dfi

        em=epp + de*XPM
        CALL dfei_int(nni,em,zp,mp,betab,zb,mb,nb,dfe,dfi) 
        fe=fe + W13*dfe
        fi=fi + W13*dfi
      ENDDO
      fe=fe*de/e0
      fi=fi*de/e0
!*!
      f=ec ! just a place holder
!*!
!      IF (te .NE. ti) THEN
!        e1=0.
!        e2=e0
!        ng=6
!        CALL i1_int(nni, zp, mp, betab, zb, mb, nb, ng, e1, e2, i1, i1onet)
!        f = i1 - (ti/ec)**1.5/(1. + (ti/ec)**1.5) + 0.5 + te/ec/2.
!        fe= fe + (ti-te)*f/e0
!        fi= fi + (te-ti)*f/e0
!      ENDIF
    END SUBROUTINE eifrac_int
    SUBROUTINE dfei_int(nni,ep,zp,mp,betab,zb,mb,nb,dfe,dfi) 
      USE bpsvars
      USE mathvars
      IMPLICIT NONE
                                                        ! Plasma:
      INTEGER,                     INTENT(IN) :: nni    !  number of ions
      REAL,    DIMENSION(1:nni+1), INTENT(IN) :: betab  !  temp array    [1/keV]
      REAL,    DIMENSION(1:nni+1), INTENT(IN) :: mb     !  mass array    [keV]
      REAL,    DIMENSION(1:nni+1), INTENT(IN) :: nb     !  density array [1/cc]
      REAL,    DIMENSION(1:nni+1), INTENT(IN) :: zb     !  charge array
                                                        ! Projectile:
      REAL,                        INTENT(IN) :: ep     !  projectile energy [keV]
      REAL,                        INTENT(IN) :: mp     !  projectile mass   [keV]
      REAL,                        INTENT(IN) :: zp     !  projectile charge
                                                        !
      REAL,                        INTENT(OUT):: dfe    ! 
      REAL,                        INTENT(OUT):: dfi    ! 

      REAL               :: te, ti, a_tot, a_i, a_e, ac_tot, ac_i, aq_tot, aq_i
      REAL               :: tae, x1, x2, f, ferf

      te=1./betab(1) ! electron temp [keV]
      ti=1./betab(2) ! ion temp in   [keV]
      dfe=0.
      dfi=0.

      CALL coeff_bps(nni,ep,zp,mp,betab,zb,mb,nb, &
           a_tot,a_i,a_e,ac_tot,ac_i,aq_tot,aq_i)
      tae=te*a_e + ti*a_i
      dfe=te*a_e/tae
      dfi=ti*a_i/tae

      x1=ep/ti
      x2=SQRT(x1)
      f=ferf(x2) - 2*x2*EXP(-x1)/SQPI
      dfe=dfe*f
      dfi=dfi*f
  END SUBROUTINE dfei_int

!
! The name "coeff_bps" for backward compatibility. Will remove later.
!
      SUBROUTINE coeff_bps(nni, ep, zp, mp, betab, zb, mb, nb, &
            a_tot, a_i, a_e, ac_tot, ac_i, ac_e, aq_tot, aq_i, aq_e)
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

        REAL,    DIMENSION(1:nni+1)  :: mpb, mbpb, kb2, ab, ab2
        REAL                         :: vp, vp2, zp2, k, k2, kd, kd2, a, b, eta
        REAL                         :: ac_r, ac_s, aq
        REAL                         :: c1, c2
        INTEGER                      :: ib, nnb

        REAL, PARAMETER              :: EPS_SMALL_E=2.E-4
        REAL, PARAMETER              :: EPS_SMALL_E_SING=2.E-4
        REAL, PARAMETER              :: EPS_SMALL_E_REG=2.E-4
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
!       kb2 = DEBYE2*zb*zb*nb*betab ! [1/cm^2]
        kb2 =8*PI*A0CM*BEKEV*zb*zb*nb*betab
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
           c1=2*zp2*BEKEV*kb2(ib)*A0CM    ! [keV/cm] c1 = e_p^2 kappa_b^2/(4 Pi)
           c1=c1*1.E-7                    ! [MeV/micron]  
           c2=SQRT(a/PI)                  ! [dimensionless] 
                                          ! c2=SQRT(betab(ib)*mb(ib)/TWOPI)*vp/CC 
!
! A-classical-singular 
!
             CALL a_sing_mass(a,b,ac_s) 
             ac_s=c1*c2*ac_s
!
! A-classical-regular 
!
           CALL a_reg(ib,nni,vp,k2,kb2,betab,mb,ac_r)
           ac_r=c1*ac_r
!
! A-quantum
!
           CALL a_quantum_mass(ib,a,eta,aq) ! eta = dimensionless quantum param.
           aq=c1*c2*aq
!
! construct electron and ion components: see common.f90 for a_collect.
!
           CALL x_collect(ib,NNB,ac_s,ac_r,aq,a_tot,a_i,a_e,&
             ac_tot,ac_i,ac_e,aq_tot,aq_i,aq_e)
        ENDIF
        ENDDO
      END SUBROUTINE coeff_bps

      SUBROUTINE a_reg(ib, nni, vp, k2, kb2, betab, mb, ac_r)
        IMPLICIT NONE
        INTEGER,                     INTENT(IN)  :: ib
        INTEGER,                     INTENT(IN)  :: nni 
        REAL,                        INTENT(IN)  :: vp
        REAL,                        INTENT(IN)  :: k2
        REAL,                        INTENT(IN)  :: kb2
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: betab
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: mb
        REAL,                        INTENT(OUT) :: ac_r
        INTEGER, PARAMETER :: NR=10 ! integration regions: must be even
        REAL,    PARAMETER :: UPM=0.7745966692E0 ! parameters for Gaussian Quad
        REAL,    PARAMETER :: W13=0.5555555556E0, W2=0.8888888889E0
        REAL               :: u0, u1, du, u, um, dei_reg
        INTEGER            :: iu
        ac_r=0
        u0=0.
        u1=1.
        du=(u1-u0)/NR
        u=u0-du
        DO iu=1,NR,2 ! Gaussian quadrature
           u=u+2.E0*du
           ac_r=ac_r+W2*dei_reg(u,vp,ib,nni,k2,kb2,betab,mb)
           um=u-du*UPM
           ac_r=ac_r+W13*dei_reg(um,vp,ib,nni,k2,kb2,betab,mb)
           um=u+du*UPM
           ac_r=ac_r+W13*dei_reg(um,vp,ib,nni,k2,kb2,betab,mb)
        ENDDO
        ac_r=ac_r*du
      END SUBROUTINE a_reg
