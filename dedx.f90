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
