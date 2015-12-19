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
