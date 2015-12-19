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
      REAL,    DIMENSION(:), ALLOCATABLE :: betab, zb, mb, nb
      REAL,    DIMENSION(:), ALLOCATABLE :: gb, etab, mpb, mbpb
    END MODULE allocatablevars

    MODULE controlvars
      LOGICAL :: asymptotic_forms=.TRUE.
      LOGICAL :: small_E_form=.TRUE.
      LOGICAL :: large_E_form=.TRUE.
    END MODULE controlvars


!
! D-T plasma:
!
!      PARAMETER (ZB(1)=-1)        ! charge of 1-st component
!      PARAMETER (ZB(2)=+1)        ! charge of 2-nd component
!      PARAMETER (ZB(3)=+1)        ! charge of 3-rd component
!      PARAMETER (MB(1)=MEKEV)     ! mass of   1-st component
!      PARAMETER (MB(2)=2*MPKEV)   ! mass of   2-nd component
!      PARAMETER (MB(3)=3*MPKEV)   ! mass of   3-rd component
!      PARAMETER (NB(1)=NE)        ! density of 1-st comp
!      PARAMETER (NB(2)=0.5*NE)    ! density of 2-nd comp
!      PARAMETER (NB(3)=0.5*NE)    ! density of 3-rd comp


! projectile - D plasma (NNB=2 or nni=1)
! 
!  nni=1         
!  ALLOCATE(mb(1:nni+1),nb(1:nni+1),zb(1:nni+1))
!
!  ee=10D0      ! projectile energy (MeV)
!  zp=2         ! projectile charge
!  mp=4*MPKEV   ! projectile mass
!  te=10.       ! temperature of electrons (keV)
!  ti=te        ! temperature of ions      (keV)
!  ne=1E26      ! electron nmr density  (cm^-3)
!  zb(1)=-1     ! species charges
!  zb(2)=+1     !
!  mb(1)=MEKEV  ! species masses
!  mb(2)=2*MPKEV!
!  nb(1)=1.     ! number density with charge neutrality
!  nb(2:nni+1)=1./(zb(2:nni+1)*nni)   


! projectile - H plasma (NNB=2 or nni=1)
! 
!  nni=1         
!  ALLOCATE(mb(1:nni+1),nb(1:nni+1),zb(1:nni+1))
!
!  ee=1E0       ! projectile energy (MeV)
!  zp=1         ! projectile charge
!  mp=MPKEV     ! projectile mass
!  te=1.        ! temperature of electrons (keV)
!  ti=te        ! temperature of ions      (keV)
!  ne=5E25      ! electron numr density (cm^-3)
!  zb(1)=-1     ! species charges
!  zb(2)=+1     !
!  mb(1)=MEKEV  ! species masses
!  mb(2)=MPKEV  !
!  nb(1)=1.     ! number density with charge neutrality
!  nb(2:nni+1)=1./(zb(2:nni+1)*nni)   

! projectile - H plasma with no electrons (NNB=0)
! 
!  nni=0         
!  ALLOCATE(mb(1:nni+1),nb(1:nni+1),zb(1:nni+1))
!
!  ee=1E0       ! projectile energy (MeV)
!  zp=1         ! projectile charge
!  mp=MPKEV     ! projectile mass
!  te=10.       ! temperature of electrons (keV)
!  ti=te        ! temperature of ions      (keV)
!  ne=5.E25     ! electron numr density (cm^-3)
!  zb(1)=1      ! species charges
!  mb(1)=MPKEV  ! species masses
!  nb(1)=1.     ! number density 

! projectile - Be plasma (NNB=2 or nni=1)
! 
!  nni=1         
!  ALLOCATE(mb(1:nni+1),nb(1:nni+1),zb(1:nni+1))
! 
!  ee=1.0E0     ! projectile energy (MeV)
!  zp=1         ! projectile charge
!  mp=MPKEV     ! projectile mass 
!  te=1.        ! temperature of electrons (keV)
!  ti=te        ! temperature of of ions   (keV)
!  ne=1.E24     ! electron numr density (cm^-3)
!  zb(1)=-1     ! species charges
!  zb(2)=+4     !
!  mb(1)=MEKEV  ! species mass
!  mb(2)=9*MPKEV!
!  nb(1)=1.     ! number density with charge neutrality
!  nb(2:nni+1)=1./(zb(2:nni+1)*nni)   

