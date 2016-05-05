SUBROUTINE CKINIT (LENIWK, LENRWK, LENCWK, LINC, LOUT, ICKWRK,&
     RCKWRK, CCKWRK)
!
!  START PROLOGUE
!
!  SUBROUTINE CKINIT (LENIWK, LENRWK, LENCWK, LINC, LOUT, ICKWRK,
!                     RCKWRK, CCKWRK)*
!     Reads the binary file and creates the internal work arrays
!     ICKWRK, CCKWRK, and RCKWORK.  CKINIT must be called before any
!     other CHEMKIN subroutine is called.  The work arrays must then
!     be made available as input to the other CHEMKIN subroutines.
!
!  INPUT
!     LENIWK - Length of the integer work array, ICKWRK.
!                   Data type - integer scalar
!     LENCWK - Length of the character work array, CCKWRK.
!              The minimum length of CCKWRK(*) is MM + KK.
!                   Data type - integer scalar
!     LENRWK - Length of the real work array, WORK.
!                   Data type - integer scalar
!     LINC  -  Logical file number for the binary file.
!                   Data type - integer scalar
!     LOUT  -  Output file for printed diagnostics.
!                   Data type - integer scalar
!
!  OUTPUT
!     ICKWRK - Array of integer workspace.
!                   Data type - integer array
!                   Dimension ICKWRK(*) at least LENIWK.
!     RCKWRK - Array of real work space.
!                   Data type - real array
!                   Dimension RCKWRK(*) at least LENRWK.
!     CCKWRK - Array of character work space.
!                   Data type - CHARACTER*16 array
!                   Dimension CCKWRK(*) at least LENCWK.
!
!  END PROLOGUE
!
  IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
!
  DIMENSION ICKWRK(*), RCKWRK(*)
  CHARACTER CCKWRK(*)*(*), VERS*16, PREC*16
  LOGICAL IOK, ROK, COK, KERR
  COMMON /CKCONS/ PREC, VERS, KERR, LENI, LENR, LENC
!
  INCLUDE 'ckstrt.fh'
!
  COMMON /MACH/ SMALL,BIG,EXPARG
  COMMON /CMIN/ CKMIN
!
!     Data about the machine dependent constants is carried in
!
!     COMMON/MACH/SMALL,BIG,EXPARG
!
  DATA RU,RUC,PA /8.314D7, 1.987, 1.01325D6/
!
!      THIS STATEMENT WILL NOT COMPILE, MACHIND-DEPENDENT CONSTANTS
!*****exponent range > +/-30
!      SMALL = 1.0D-30
!      BIG   = 1.0E+30
!*****END exponent range > +/-30
!*****exponent range > +/-300
  SMALL = 10.0D0**(-14)
  BIG   = 10.0D0**(+60)
!*****END exponent range > +/-300
  EXPARG = LOG(BIG)
!
  !if(myid == 0) then
  !   WRITE (LOUT,15)
  !end if
15 FORMAT (/1X,' CKLIB:  Chemical Kinetics Library',&
       /1X,'         CHEMKIN-II Version 4.9, April 1994',&
       /1X,'         REAL*8')
!
  CALL CKLEN (LINC, LOUT, LI, LR, LC)
!
  IOK = (LENIWK .GE. LI)
  ROK = (LENRWK .GE. LR)
  COK = (LENCWK .GE. LC)
  IF (.NOT.IOK .OR. .NOT.ROK .OR. .NOT.COK) THEN
     IF (.NOT. IOK) WRITE (LOUT, 300) LI
     IF (.NOT. ROK) WRITE (LOUT, 350) LR
     IF (.NOT. COK) WRITE (LOUT, 375) LC
     STOP
  ENDIF
!
  REWIND LINC
  READ (LINC, ERR=110) VERS, PREC, KERR
  READ (LINC, ERR=110) LENI, LENR, LENC, MM, KK, II,&
       MAXSP, MAXTB, MAXTP, NTHCF, NIPAR, NITAR,&
       NIFAR, NRV, NFL, NTB, NLT, NRL, NW, NCHRG,&
       NSTO, NOR, MAXORD, CKMN
!
  IF (LEN(CCKWRK(1)) .LT. 16) THEN
     WRITE (LOUT,475)
     STOP
  ENDIF
!
  NMM = MM
  NKK = KK
  NII = II
  MXSP = MAXSP
  MXTB = MAXTB
  MXTP = MAXTP
  MXOR = MAXORD
  NCP  = NTHCF
  NCP1 = NTHCF+1
  NCP2 = NTHCF+2
  NCP2T = NCP2*(MAXTP-1)
  NPAR = NIPAR
  NLAR = NITAR
  NFAR = NIFAR
  NTHB = NTB
  NLAN = NLT
  NFAL = NFL
  NREV = NRV
  NRLT = NRL
  NWL  = NW
  NRNU= NSTO
  NORD = NOR
  MXORD= MAXORD
  CKMIN= CKMN
!
!             APPORTION work arrays
!
!             SET  ICKWRK(*)=1  TO FLAG THAT CKINIT HAS BEEN CALLED
!
  ICKWRK(1) = 1
!
!             STARTING LOCATIONS OF INTEGER SPACE
!
  !! elemental composition of species
  IcNC = 2
  !! species phase array
  IcPH = IcNC + KK*MM
  !! species charge array
  IcCH = IcPH + KK
  !! # of temperatures for fit
  IcNT = IcCH + KK
  !! stoichiometric coefficients
  IcNU = IcNT + KK
  !! species numbers for the coefficients
  IcNK = IcNU + MAXSP*II
  !! # of non-zero coefficients  (<0=reversible, >0=irreversible)
  IcNS = IcNK + MAXSP*II
  !! # of reactants
  IcNR = IcNS + II
  !! Landau-Teller reaction numbers
  IcLT = IcNR + II
  !! Reverse Landau-Teller reactions
  IcRL = IcLT + NLAN
  !! Fall-off reaction numbers
  IcFL = IcRL + NRLT
  !! Fall-off option numbers
  IcFO = IcFL + NFAL
  !! Fall-off enhanced species
  IcKF = IcFO + NFAL
  !! Third-body reaction numbers
  IcTB = IcKF + NFAL
  !! number of 3rd bodies for above
  IcKN = IcTB + NTHB
  !! array of species #'s for above
  IcKT = IcKN + NTHB
  !! Reverse parameter reaction numbers
  IcRV = IcKT + MAXTB*NTHB
  !! Radiation wavelength reactions
  IcWL = IcRV + NREV
  !! Real stoichometry reactions
  IcRNU= IcWL + NWL
  !! Change of order reactions
  IcORD= IcRNU + NRNU
  !! Species for which there is a change of order
  IcKOR= IcORD + NORD
!
  ITOT = IcKOR + NORD*MXORD - 1
!
!             STARTING LOCATIONS OF CHARACTER SPACE
!
  !! start of element names
  IcMM = 1
  !! start of species names
  IcKK = IcMM + MM
  ITOC = IcKK + KK - 1
!
!             STARTING LOCATIONS OF REAL SPACE
!
  !! atomic weights
  NcAW = 1
  !! molecular weights
  NcWT = NcAW + MM
  !! temperature fit array for species
  NcTT = NcWT + KK
  !! thermodynamic coefficients
  NcAA = NcTT + MAXTP*KK
  !! Arrhenius coefficients (3)
  NcCO = NcAA + (MAXTP-1)*NCP2*KK
  !! Reverse coefficients
  NcRV = NcCO + (NPAR+1)*II
  !! Landau-Teller #'s for NLT reactions
  NcLT = NcRV + (NPAR+1)*NREV
  !! Reverse Landau-Teller #'s
  NcRL = NcLT + NLAR*NLAN
  !! Fall-off parameters for NFL reactions
  NcFL = NcRL + NLAR*NRLT
  !! 3rd body coef'nts for NTHB reactions
  NcKT = NcFL + NFAR*NFAL
  !! wavelength
  NcWL = NcKT + MAXTB*NTHB
  !! real stoichometric coefficients
  NcRNU= NcWL + NWL
  !! change of order for species/reactions
  NcKOR= NcRNU + NRNU*MXSP
  !! universal gas constant
  NcRU = NcKOR + NORD*MXORD
  !! universal gas constant in units
  NcRC = NcRU + 1
  !! pressure of one atmosphere
  NcPA = NcRC + 1
  !! intermediate temperaturD-dependent forward rates
  NcKF = NcPA + 1
  !! intermediate temperaturD-dependent reverse rates
  NcKR = NcKF + II
  !! internal work space of length kk
  NcK1 = NcKR + II
  !!          'ditto'
  NcK2 = NcK1 + KK
  !!          'ditto'
  NcK3 = NcK2 + KK
  !!          'ditto'
  NcK4 = NcK3 + KK
  NcI1 = NcK4 + KK
  NcI2 = NcI1 + II
  NcI3 = NcI2 + II
  NcI4 = NcI3 + II
  NTOT = NcI4 + II - 1
!
!        SET UNIVERSAL CONSTANTS IN CGS UNITS
!
  RCKWRK(NcRU) = RU
  RCKWRK(NcRC) = RUC
  RCKWRK(NcPA) = PA
!
  !!element names, !atomic weights
  READ (LINC,err=111) (CCKWRK(IcMM+M-1), RCKWRK(NcAW+M-1), M=1,MM)
!
  !!species names, !composition, !phase, !charge, !molec weight,
  !!# of fit temps, !array of temps, !fit coeff'nts
  READ (LINC,err=222) (CCKWRK(IcKK+K-1),&
       (ICKWRK(IcNC+(K-1)*MM + M-1),M=1,MM),&
       ICKWRK(IcPH+K-1),&
       ICKWRK(IcCH+K-1),&
       RCKWRK(NcWT+K-1),&
       ICKWRK(IcNT+K-1),&
       (RCKWRK(NcTT+(K-1)*MAXTP + L-1),L=1,MAXTP),&
       ((RCKWRK(NcAA+(L-1)*NCP2+(K-1)*NCP2T+N-1),&
       N=1,NCP2), L=1,(MAXTP-1)),    K = 1,KK)
!
  IF (II .EQ. 0) RETURN
!
  !!# spec,reactants, !Arr. coefficients, !stoic coef, !species numbers
  READ (LINC,end=100,err=333)&
       (ICKWRK(IcNS+I-1), ICKWRK(IcNR+I-1),&
       (RCKWRK(NcCO+(I-1)*(NPAR+1)+N-1), N=1,NPAR),&
       (ICKWRK(IcNU+(I-1)*MAXSP+N-1),&
       ICKWRK(IcNK+(I-1)*MAXSP+N-1), N=1,MAXSP),&
       I = 1,II)
!
!     PERTURBATION FACTOR
!
  DO  I = 1, II
     RCKWRK(NcCO + (I-1)*(NPAR+1) + NPAR) = 1.0
  end DO
!
  IF (NREV .GT. 0) READ (LINC,err=444)&
       (ICKWRK(IcRV+N-1), (RCKWRK(NcRV+(N-1)*(NPAR+1)+L-1),&
       L=1,NPAR), N = 1,NREV)
!
  IF (NFAL .GT. 0) READ (LINC,err=555)&
       (ICKWRK(IcFL+N-1), ICKWRK(IcFO+N-1), ICKWRK(IcKF+N-1),&
       (RCKWRK(NcFL+(N-1)*NFAR+L-1),L=1,NFAR),N=1,NFAL)
!
  IF (NTHB .GT. 0) READ (LINC,err=666)&
       (ICKWRK(IcTB+N-1), ICKWRK(IcKN+N-1),&
       (ICKWRK(IcKT+(N-1)*MAXTB+L-1),&
       RCKWRK(NcKT+(N-1)*MAXTB+L-1),L=1,MAXTB),N=1,NTHB)
!
  IF (NLAN .GT. 0) READ (LINC,err=777)&
       (ICKWRK(IcLT+N-1), (RCKWRK(NcLT+(N-1)*NLAR+L-1),L=1,NLAR),&
       N=1,NLAN)
!
  IF (NRLT .GT. 0) READ (LINC,err=888)&
       (ICKWRK(IcRL+N-1), (RCKWRK(NcRL+(N-1)*NLAR+L-1),L=1,NLAR),&
       N=1,NRLT)
!
  IF (NWL .GT. 0) READ (LINC,err=999)&
       (ICKWRK(IcWL+N-1), RCKWRK(NcWL+N-1), N=1,NWL)
!
  IF (NRNU .GT. 0) READ (LINC,err=1111)&
       (ICKWRK(IcRNU+N-1), (RCKWRK(NcRNU+(N-1)*MAXSP+L-1),L=1,MAXSP),&
       N=1,NRNU)
!
  IF (NORD .GT. 0) READ (LINC,err=2222)&
       (ICKWRK(IcORD+N-1), (ICKWRK(IcKOR+(N-1)*MXORD+L-1),&
       RCKWRK(NcKOR+(N-1)*MXORD+L-1),&
       L = 1, MXORD), N=1, NORD)
!
100 CONTINUE
  RETURN
!
110 WRITE (LOUT,*) ' Error reading binary file...'
  STOP
111 WRITE (LOUT,*) ' Error reading element data...'
  STOP
222 WRITE (LOUT,*) ' Error reading species data...'
  STOP
333 WRITE (LOUT,*) ' Error reading reaction data...'
  STOP
444 WRITE (LOUT,*) ' Error reading reverse Arrhenius parameters...'
  STOP
555 WRITE (LOUT,*) ' Error reading Fall-off data...'
  STOP
666 WRITE (LOUT,*) ' Error reading third-body data...'
  STOP
777 WRITE (LOUT,*) ' Error reading Landau-Teller data...'
  STOP
888 WRITE (LOUT,*) ' Error reading reverse Landau-Teller data...'
  STOP
999 WRITE (LOUT,*) ' Error reading Wavelength data...'
  STOP
1111 WRITE (LOUT,*) ' Error reading real stoichometric data...'
  STOP
2222 WRITE (LOUT,*) ' Error reading order data...'
  STOP
!
300 FORMAT (10X,'ICKWRK MUST BE DIMENSIONED AT LEAST ',I5)
350 FORMAT (10X,'RCKWRK MUST BE DIMENSIONED AT LEAST ',I5)
375 FORMAT (10X,'CCKWRK MUST BE DIMENSIONED AT LEAST ',I5)
475 FORMAT (10X,'CHARACTER LENGTH OF CCKWRK MUST BE AT LEAST 16 ')
END SUBROUTINE CKINIT

!
!----------------------------------------------------------------------C
!
SUBROUTINE CKLEN (LINC, LOUT, LI, LR, LC)
!
!  START PROLOGUE
!
!  SUBROUTINE CKLEN (LINC, LOUT, LENI, LENR, LENC)
!     Returns the lengths required for the work arrays.
!
!  INPUT
!
!     LINC  -  Logical file number for the binary file.
!                   Data type - integer scalar
!     LOUT  -  Output file for printed diagnostics.
!                   Data type - integer scalar
!
!  OUTPUT
!     LENI  -  Minimum length required for the integer work array.
!                   Data type - integer scalar
!     LENR  -  Minimum length required for the real work array.
!                   Data type - integer scalar
!     LENC  -  Minimum length required for the character work array.
!
!  END PROLOGUE
!
  IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
!
  PARAMETER (NLIST = 3)
  LOGICAL KERR, VOK, POK
  CHARACTER LIST(NLIST)*16, PREC*16, VERS*16
  COMMON /CKCONS/ PREC, VERS, KERR, LENI, LENR, LENC
!      DATA LIST/'1.9','2.0','2.1','2.2','2.3','2.4','2.5','2.6',
!     1          '2.7','2.8','2.9','3.0','3.1','3.2','3.3'/
  DATA LIST /'3.4','3.5','3.6'/
!
  VERS = ' '
  PREC = ' '
  LENI = 0
  LENR = 0
  LENC = 0
!
  KERR = .FALSE.
  REWIND LINC
  READ (LINC, ERR=999) VERS, PREC, KERR
!
  VOK = .FALSE.
  DO N = 1, NLIST
     IF (VERS .EQ. LIST(N)) VOK = .TRUE.
  enddo
!
  POK = .FALSE.
  IF (INDEX(PREC, 'DOUB') .GT. 0) POK = .TRUE.
!
  IF (KERR .OR. (.NOT.POK) .OR. (.NOT.VOK)) THEN
     IF (KERR) THEN
        WRITE (LOUT,'(/A,/A)') &
             ' There is an error in the Chemkin binary file...', &
             ' Check CHEMKIN INTERPRETER output for error conditions.'
     ENDIF
     IF (.NOT. VOK) THEN
        WRITE (LOUT,'(/A,A)') &
             ' Chemkin binary file is incompatible with Chemkin', &
             ' Library Version 4.9'
     ENDIF
     IF (.NOT. POK) THEN
        WRITE (LOUT,'(/A,A)') &
             ' Precision of Chemkin binary file does not agree with', &
             ' precision of Chemkin library'
     ENDIF
     STOP
  ENDIF
!
  READ (LINC, ERR=999) LENICK, LENRCK, LENCCK, MM, KK, II, &
       MAXSP, MAXTB, MAXTP, NTHCF, NIPAR, NITAR, &
       NIFAR, NRV, NFL, NTB, NLT, NRL, NW, NCHRG, &
       NSTO, NOR, MAXORD, CKMN
  REWIND LINC
!
  LENI = LENICK
  LENR = LENRCK
  LENC = LENCCK
  LI   = LENI
  LR   = LENR
  LC   = LENC
  RETURN
!
999 CONTINUE
  WRITE (LOUT, 50)
50 FORMAT (' Error reading Chemkin binary file.')
  STOP
END SUBROUTINE CKLEN
!----------------------------------------------------------------------C
!
!
!----------------------------------------------------------------------C
!
SUBROUTINE CKINDX (ICKWRK, RCKWRK, MM, KK, II, NFIT)
!
!  START PROLOGUE
!
!  SUBROUTINE CKINDX (ICKWRK, RCKWRK, MM, KK, II, NFIT)*
!     Returns a group of indices defining the size of the particular
!     reaction mechanism
!
!  INPUT
!     ICKWRK - Array of integer workspace.
!                   Data type - integer array
!                   Dimension ICKWRK(*) at least LENIWK.
!     RCKWRK - Array of real work space.
!                   Data type - real array
!                   Dimension RCKWRK(*) at least LENRWK.
!
!  OUTPUT
!     MM     - Total number of elements in mechanism.
!                   Data type - integer scalar
!     KK     - Total number of species in mechanism.
!                   Data type - integer scalar
!     II     - Total number of reactions in mechanism.
!                   Data type - integer scalar
!     NFIT   - number of coefficients in fits to thermodynamic data
!              for one temperature range; NFIT = number of
!              coefficients in polynomial fits to CP/R  +  2.
!                   Data type - integer scalar
!
!  END PROLOGUE
!
  IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
!
  INCLUDE 'ckstrt.fh'
!
  DIMENSION ICKWRK(*), RCKWRK(*)
!
  MM = NMM
  KK = NKK
  II = NII
  NFIT = NCP2
  RETURN
END SUBROUTINE CKINDX
!
!----------------------------------------------------------------------C
!
SUBROUTINE CKSYMS (CCKWRK, LOUT, KNAME, KERR)
!
!  START PROLOGUE
!
!  SUBROUTINE CKSYMS (CCKWRK, LOUT, KNAME, KERR)*
!     Returns the character strings of species names
!
!  INPUT
!     CCKWRK - Array of character work space.
!                   Data type - CHARACTER*16 array
!                   Dimension CCKWRK(*) at least LENCWK.
!     LOUT   - Output unit for printed diagnostics.
!                   Data type - integer scalar
!
!  OUTPUT
!     KNAME  - Species names.
!                   Data type - CHARACTER*(*) array
!                   Dimension KNAME(*) at least KK,
!                   the total number of species.
!     KERR   - Error flag; character length errors will result in
!              KERR = .TRUE.
!                   Data type - logical
!
!  END PROLOGUE
!
  IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
!
  CHARACTER*(*) CCKWRK(*), KNAME(*)
  LOGICAL KERR
!
  INCLUDE 'ckstrt.fh'
!
  KERR = .FALSE.
  ILEN = LEN(KNAME(1))
  DO  K = 1, NKK
     LT = ILASCH(CCKWRK(IcKK + K - 1))
     KNAME(K) = ' '
     IF (LT .LE. ILEN) THEN
        KNAME(K) = CCKWRK(IcKK+K-1)
     ELSE
        WRITE (LOUT,*) &
             ' Error in CKSYM...character string length too small '
        KERR = .TRUE.
     ENDIF
  end DO
  RETURN
END SUBROUTINE CKSYMS
!
!----------------------------------------------------------------------C
!
SUBROUTINE CKWT   (ICKWRK, RCKWRK, WT)
!
!  START PROLOGUE
!
!  SUBROUTINE CKWT   (ICKWRK, RCKWRK, WT)
!     Returns the molecular weights of the species
!
!  INPUT
!     ICKWRK - Array of integer workspace.
!                   Data type - integer array
!                   Dimension ICKWRK(*) at least LENIWK.
!     RCKWRK - Array of real work space.
!                   Data type - real array
!                   Dimension RCKWRK(*) at least LENRWK.
!
!  OUTPUT
!     WT     - Molecular weights of the species.
!                   cgs units - gm/mole
!                   Data type - real array
!                   Dimension WT(*) at least KK, the total number of
!                   species.
!
!  END PROLOGUE
!
  IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
!
  INCLUDE 'ckstrt.fh'
!
  DIMENSION ICKWRK(*), RCKWRK(*), WT(*)
!
  DO  K = 1, NKK
     WT(K) = RCKWRK(NcWT + K - 1)
  end DO
!
  RETURN
END SUBROUTINE CKWT
!
!----------------------------------------------------------------------C
!
SUBROUTINE CKMMWX (X, ICKWRK, RCKWRK, WTM)
!
!  START PROLOGUE
!
!  SUBROUTINE CKMMWX (X, ICKWRK, RCKWRK, WTM)
!     Returns the mean molecular weight of the gas mixture given the
!     mole fractions;  see Eq. (4).
!
!  INPUT
!     X      - Mole fractions of the species.
!                   cgs units - none
!                   Data type - real array
!                   Dimension X(*) at least KK, the total number of
!                   species.
!     ICKWRK - Array of integer workspace.
!                   Data type - integer array
!                   Dimension ICKWRK(*) at least LENIWK.
!     RCKWRK - Array of real work space.
!                   Data type - real array
!                   Dimension RCKWRK(*) at least LENRWK.
!
!  OUTPUT
!     WTM    - Mean molecular weight of the species mixture.
!                   cgs units - gm/mole
!                   Data type - real scalar
!
!  END PROLOGUE
!
  IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
!
  DIMENSION X(*), ICKWRK(*), RCKWRK(*)
!
  INCLUDE 'ckstrt.fh'
!
  WTM = 0.0d0
  DO  K = 1, NKK
     WTM = WTM + X(K)*RCKWRK(NcWT + K - 1)
  end DO
  RETURN
END SUBROUTINE CKMMWX
!
!----------------------------------------------------------------------C
!
!
!----------------------------------------------------------------------C
!
SUBROUTINE CKRP   (ICKWRK, RCKWRK, RU, RUC, PA)
!
!  START PROLOGUE
!
!  SUBROUTINE CKRP   (ICKWRK, RCKWRK, RU, RUC, PA)
!     Returns universal gas constants and the pressure of one standard
!     atmosphere
!
!  INPUT
!     ICKWRK - Array of integer workspace.
!                   Data type - integer array
!                   Dimension ICKWRK(*) at least LENIWK.
!     RCKWRK - Array of real work space.
!                   Data type - real array
!                   Dimension RCKWRK(*) at least LENRWK.
!
!  OUTPUT
!     RU     - Universal gas constant.
!                   cgs units - 8.314D7 ergs/(mole*K)
!                   Data type - real scalar
!     RUC    - Universal gas constant used only in conjuction with
!              activation energy.
!                   preferred units - 1.987 cal/(mole*K)
!                   Data type - real scalar
!     PA     - Pressure of one standard atmosphere.
!                   cgs units - 1.01325D6 dynes/cm**2
!                   Data type - real scalar
!
!  END PROLOGUE
!
  IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
!
  INCLUDE 'ckstrt.fh'
!
  DIMENSION ICKWRK(*), RCKWRK(*)
!
  RU  = RCKWRK(NcRU)
  RUC = RCKWRK(NcRC)
  PA  = RCKWRK(NcPA)
  RETURN
END SUBROUTINE CKRP
!
!----------------------------------------------------------------------C
!
SUBROUTINE MCINIT (LINKMC, LOUT, LENIMC, LENRMC, IMCWRK, RMCWRK)
!
  !Use MYMPI, only : myid
  IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  SUBROUTINE MCINIT (LINKMC, LOUT, LENIMC, LENRMC, IMCWRK, RMCWRK)
!
!    THIS SUBROUTINE SERVES TO READ THE LINKING FILE FROM THE FITTING
!    CODE AND TO CREATE THE INTERNAL STORAGE AND WORK ARRAYS, IMCWRK(*)
!    AND RMCWRK(*).  MCINIT MUST BE CALLED BEFORE ANY OTHER TRANSPORT
!    SUBROUTINE IS CALLED.  IT MUST BE CALLED AFTER THE CHEMKIN PACKAGE
!    IS INITIALIZED.
!
!  INPUT-
!    LINKMC  - LOGICAL UNIT NUMBER OF THE LINKING FILE.
!                  FITTING CODE WRITE TO DEFAULT UNIT 35
!    LOUT    - LOGICAL UNIT NUMBER FOR PRINTED OUTPUT.
!    LENIMC  - ACTUAL DIMENSION OF THE INTEGER STORAGE AND WORKING
!              SPACE, ARRAY IMCWRK(*).  LENIMC MUST BE AT LEAST:
!                LENIMC = 4*KK + NLITE
!                 WHERE, KK    = NUMBER OF SPECIES.
!                        NLITE = NUMBER OF SPECIES WITH MOLECULAR WEIGHT
!                                LESS THAN 5.
!    LENRMC  - ACTUAL DIMENSION OF THE FLOATING POINT STORAGE AND
!              WORKING SPACE, ARRAY RMCWRK(*).  LENRMC MUST BE AT LEAST:
!                LENRMC = KK*(19 + 2*NO + NO*NLITE) + (NO+15)*KK**2
!                 WHERE, KK    = NUMBER OF SPECIES.
!                        NO    = ORDER OF THE POLYNOMIAL FITS,
!                                DEFAULT, NO=4.
!                        NLITE = NUMBER OF SPECIES WITH MOLECULAR WEIGHT
!                                LESS THAN 5.
!
!  WORK-
!    IMCWRK  - ARRAY OF INTEGER STORAGE AND WORK SPACE.  THE STARTING
!              ADDRESSES FOR THE IMCWRK SPACE ARE STORED IN
!              COMMON /MCMCMC/.
!                  DIMENSION IMCWRK(*) AT LEAST LENIMC.
!    RMCWRK  - ARRAY OF FLOATING POINT STORAGE AND WORK SPACE.  THE
!              STARTING ADDRESSES FOR THE RMCWRK SPACE ARE STORED IN
!              COMMON /MCMCMC/.
!                  DIMENSION RMCWRK(*) AT LEAST LENRMC.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
  DIMENSION IMCWRK(*), RMCWRK(*)
  CHARACTER*16 VERS, PREC
  LOGICAL IOK, ROK, KERR
  COMMON /MCCONS/ VERS, PREC, KERR, LENI, LENR
!
  COMMON /MCMCMC/ RU, PATMOS, SMALL, NKK, NO, NLITE, INLIN, IKTDIF, &
       IPVT, NWT, NEPS, NSIG, NDIP, NPOL, NZROT, NLAM, &
       NETA, NDIF, NTDIF, NXX, NVIS, NXI, NCP, &
       NCROT, NCINT, NBIND, NEOK, NSGM, &
       NAST, NBST, NCST, NXL, NR, NWRK, K3
!
!
!         THE FOLLOWING NUMBER SMALL IS USED IN THE MIXTURE DIFFUSION
!        COEFFICIENT CALCULATION.  ITS USE ALLOWS A SMOOTH AND WELL
!        DEFINED DIFFUSION COEFFICIENT AS THE MIXTURE APPROACHES A
!        PURE SPECIES, EVEN THOUGH STRICTLY SPEAKING THERE DOES NOT
!        EXIST A DIFFUSION COEFFICIENT IN THIS CASE.  THE VALUE OF
!        "SMALL" SHOULD BE SMALL RELATIVE TO ANY SPECIES MOLE FRACTION
!        OF IMPORTANCE, BUT LARGE ENOUGH TO BE REPRESENTED ON THE
!        COMPUTER.
!
  SMALL = 1.0D-14
  RU    = 8.314E+07
  PATMOS= 1.01325E+06
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!          WRITE VERSION NUMBER
!
  !if(myid == 0)then
  !   WRITE (LOUT, 15)
  !end if
15 FORMAT(/' TRANLIB:  Multicomponent transport library,',&
       /'           CHEMKIN-II Version 3.1, March 1994',&
       /'           REAL*8')
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!        READ THE PROBLEM SIZE
!
  CALL MCLEN (LINKMC, LOUT, LI, LR)
  IOK = (LENIMC .GE. LI)
  ROK = (LENRMC .GE. LR)
!
  IF (.NOT.IOK .OR. .NOT.ROK) THEN
     IF (.NOT. IOK) WRITE (LOUT, 300) LI
     IF (.NOT. ROK) WRITE (LOUT, 350) LR
     STOP
  ENDIF
!
  REWIND LINKMC
  READ (LINKMC, ERR=999) VERS, PREC, KERR
  READ (LINKMC, ERR=999) LI, LR, NO, NKK, NLITE
!
  NK  = NO*NKK
  NK2 = NO*NKK*NKK
  K2  = NKK*NKK
  K3  = 3*NKK
  K32 = K3*K3
  NKT = NO*NKK*NLITE
!
!        APPORTION THE REAL WORK SPACE
!           THE POINTERS HAVE THE FOLLOWING MEANINGS:
!             NWT   - THE SPECIES MOLECULAR WEIGHTS.
!             NEPS  - THE EPSILON/K WELL DEPTH FOR THE SPECIES.
!             NSIG  - THE COLLISION DIAMETER FOR THE SPECIES.
!             NDIP  - THE DIPOLE MOMENTS FOR THE SPECIES.
!             NPOL  - THE POLARIZABILITIES FOR THE SPECIES.
!             NZROT - THE ROTATIONAL RELAXATION COLLISION NUMBERS.
!             NLAM  - THE COEFFICIENTS FOR THE CONDUCTIVITY FITS.
!             NETA  - THE COEFFICIENTS FOR THE VISCOSITY FITS.
!             NTDIF - THE COEFFICIENTS FOR THE THERMAL DIFFUSION
!                     RATIO FITS.
!             NXX   - THE MOLE FRACTIONS.
!             NVIS  - THE SPECIES VISCOSITIES.
!             NXI   - THE ROTATIONAL RELAXATION COLLISION NUMBERS BEFORE
!                     THE PARKER COFFECTION.
!             NCP   - THE SPECIES SPECIFIC HEATS.
!             NCROT - THE ROTATIONAL PARTS OF THE SPECIFIC HEATS.
!             NCINT - THE INTERNAL PARTS OF THE SPECIFIC HEATS.
!             NBIND - THE BINARY DIFFUSION COEFFICIENTS.
!             NEOK  - THE MATRIX OF REDUCED WELL DEPTHS.
!             NSGM  - THE MATRIX OF REDUCED COLLISION DIAMETERS.
!             NAST  - THE MATRIX OF A* COLLISION INTEGRALS FOR EACH
!                     SPECIES PAIR.
!             NBST  - THE MATRIX OF B* COLLISION INTEGRALS FOR EACH
!                     SPECIES PAIR.
!             NCST  - THE MATRIX OF C* COLLISION INTEGRALS FOR EACH
!                     SPECIES PAIR.
!             NXL   - THE "L" MATRIX.
!             NR    - THE RIGHT HAND SIDES OF THE LINEAR SYSTEM
!                     INVOLVING THE "L" MATRIX.
!             NWRK  - THE WORK SPACE NEEDED BY LINPACK TO SOLVE THE
!                     "L" MATRIX LINEAR SYSTEM.
!
  NWT  = 1
  NEPS = NWT + NKK
  NSIG = NEPS + NKK
  NDIP = NSIG + NKK
  NPOL = NDIP + NKK
  NZROT= NPOL + NKK
!
  NLAM = NZROT + NKK
  NETA = NLAM + NK
  NDIF = NETA + NK
  NTDIF= NDIF + NK2
!
  NXX  = NTDIF + NO*NKK*NLITE
  NVIS = NXX + NKK
  NXI  = NVIS + NKK
  NCP  = NXI + NKK
  NCROT= NCP + NKK
  NCINT= NCROT + NKK
!
  NBIND= NCINT + NKK
  NEOK = NBIND + K2
  NSGM = NEOK + K2
  NAST = NSGM + K2
  NBST = NAST + K2
  NCST = NBST + K2
!
  NXL  = NCST + K2
!
  NR   = NXL + K32
  NWRK = NR + K3
!      NTOT = NWRK + K3 - 1
!
!           APPORTION THE INTEGER WORK SPACE
!              THE POINTERS HAVE THE FOLLOWING MEANING:
!
!                INLIN - THE INDICATORS FOR THE MOLECULE LINEARITY.
!                IKTDIF- THE SPECIES INDICIES FOR THE "LIGHT" SPECIES.
!                IPVT  - THE PIVOT INDICIES FOR LINPACK CALLS.
!
  INLIN = 1
  IKTDIF= INLIN + NKK
  IPVT  = IKTDIF + NLITE
!      ITOT  = IPVT + K3 - 1
!
!        READ THE DATA FROM THE LINK FILE
!

  READ (LINKMC, ERR=999) PATMOS, (RMCWRK(NWT+N-1), &
       RMCWRK(NEPS+N-1), RMCWRK(NSIG+N-1), &
       RMCWRK(NDIP+N-1), RMCWRK(NPOL+N-1), RMCWRK(NZROT+N-1), &
       IMCWRK(INLIN+N-1), N=1,NKK), &
       (RMCWRK(NLAM+N-1), N=1,NK), (RMCWRK(NETA+N-1), N=1,NK), &
       (RMCWRK(NDIF+N-1), N=1,NK2), &
       (IMCWRK(IKTDIF+N-1), N=1,NLITE), (RMCWRK(NTDIF+N-1), N=1,NKT)
!
!        SET EPS/K AND SIG FOR ALL I,J PAIRS
!
  CALL MCEPSG (NKK, RMCWRK(NEPS), RMCWRK(NSIG), RMCWRK(NDIP), &
       RMCWRK(NPOL), RMCWRK(NEOK), RMCWRK(NSGM) )
!
300 FORMAT (10X,'IMCWRK MUST BE DIMENSIONED AT LEAST ', I5)
350 FORMAT (10X,'RMCWRK MUST BE DIMENSIONED AT LEAST ', I5)
  RETURN
999 CONTINUE
  WRITE (LOUT, *) ' Error reading Transport binary file...'
  STOP
END SUBROUTINE MCINIT
!
!---------------------------------------------------------------------
!
SUBROUTINE MCLEN (LINKMC, LOUT, LI, LR)
!
  IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
!
  PARAMETER (NLIST = 1)
  LOGICAL KERR, VOK, POK
  CHARACTER*16 LIST(NLIST), VERS, PREC
  COMMON /MCCONS/ VERS, PREC, KERR, LENI, LENR
!      DATA LIST/'1.7','1.8','1.9'/
  DATA LIST/'3.0'/
!
  VERS = ' '
  PREC = ' '
  LENI = 0
  LENR = 0
  LI   = LENI
  LR   = LENR
!
  REWIND (LINKMC)
  READ (LINKMC, ERR=999) VERS, PREC, KERR
!
  VOK = .FALSE.
  DO N = 1, NLIST
     IF (VERS .EQ. LIST(N)) VOK = .TRUE.
  end DO
!
  POK = .FALSE.
  IF (INDEX(PREC, 'DOUB') .GT. 0) POK = .TRUE.
!
  IF (KERR .OR. (.NOT.POK) .OR. (.NOT.VOK)) THEN
     IF (KERR) THEN
        WRITE (LOUT,'(/A,/A)') &
             ' There is an error in the transport binary file...', &
             ' Check TRANFIT output for error conditions.'
     ENDIF
     IF (.NOT. VOK) THEN
        WRITE (LOUT,'(/A,A)') &
             ' Transport binary file is incompatible with Transport', &
             ' Library Version 1.7'
     ENDIF
     IF (.NOT. POK) THEN
        WRITE (LOUT, '(/A,A)') &
             ' Precision of Transport binary file does not agree with', &
             ' precision of Transport Library'
     ENDIF
     STOP
  ENDIF
!
  READ (LINKMC, ERR=999) LENIMC, LENRMC, NO, NKK, NLITE
  REWIND (LINKMC)
  LENI = LENIMC
  LENR = LENRMC
  LI   = LENI
  LR   = LENR
  RETURN
!
999 CONTINUE
  WRITE (LOUT, 50)
50 FORMAT &
       (' Error reading Multi-component Transport binary file.')
  STOP
END SUBROUTINE MCLEN
!
!**************************************************
!
SUBROUTINE MCPREP (RMCWRK)
!
  !Use DRIVTP,ONLY : MWT1JK,MWT2JK
  IMPLICIT NONE
!
  REAL*8 ::  RMCWRK(*),SQRT8R,ZERO,ONE
  INTEGER :: J,K
!
  real*8 :: RU, PATMOS, SMALL
  INTEGER :: NKK, NO, NLITE, INLIN, IKTDIF,&
       IPVT, NWT, NEPS, NSIG, NDIP, NPOL, NZROT, NLAM,&
       NETA, NDIF, NTDIF, NXX, NVIS, NXI, NCP, NCROT,&
       NCINT, NBIND, NEOK, NSGM, NAST, NBST,&
       NCST, NXL, NR, NWRK, K3
  COMMON /MCMCMC/ RU, PATMOS, SMALL, NKK, NO, NLITE, INLIN, IKTDIF,&
       IPVT, NWT, NEPS, NSIG, NDIP, NPOL, NZROT, NLAM,&
       NETA, NDIF, NTDIF, NXX, NVIS, NXI, NCP, NCROT,&
       NCINT, NBIND, NEOK, NSGM, NAST, NBST,&
       NCST, NXL, NR, NWRK, K3
!
  ZERO   = 0.0d0
  ONE    = 1.0d0
  SQRT8R = ONE/SQRT(8d0)

  !DO  K = 1, NKK
  !   DO  J = 1, NKK
  !      MWT1JK(J,K) = (RMCWRK(NWT+J-1)/RMCWRK(NWT+K-1))**0.25D0
  !      MWT2JK(J,K) = SQRT8R/SQRT(ONE + RMCWRK(NWT+K-1) / RMCWRK(NWT+J-1))
  !   end DO
  !end DO
!
!
  RETURN
END SUBROUTINE MCPREP
!
!---------------------------------------------------------------------
!
SUBROUTINE MCEPSG (KK, EPS, SIG, DIP, POL, EOK, SGM)
!
  IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  SUBROUTINE MCEPSG (KK, EPS, SIG, DIP, POL, EOK, SGM)
!
!    THIS SUBROUTINE COMPUTES THE REDUCED WELL DEPTH EOK(I,J) AND
!    COLLISION DIAMETER SGM(I,J) FOR EACH I,J SPECIES PAIR.  THE
!    ROUTINE IS CALLED ONLY ONCE BY THE INITIALIZATION SUBROUTINE MCINIT
!    THIS ROUTINE IS NORMALLY NOT CALLED BY THE USER.
!
!  INPUT-
!    KK      - NUMBER OF SPECIES
!    EPS     - ARRAY OF LENNARD-JONES POTENTIAL WELL DEPTHS.
!                  CGS UNITS - K.
!                  DIMENSION EPS(*) AT LEAST KK
!    SIG     - ARRAY OF LENNARD-JONES COLLISION DIAMETERS.
!                  UNITS - ANGSTROMS.
!                  DIMENSION SIG(*) AT LEAST KK
!    DIP     - ARRAY OF DIPOLE MOMENTS
!                  UNITS - DEBYE
!                  DIMENSION DIP(*) AT LEAST KK
!    POL     - ARRAY OF POLARIZABILITIES.
!                  UNITS - ANGSTROMS**3.
!                  DIMENSION POL(*) AT LEAST KK
!
!  OUTPUT-
!    EOK     - MATRIX OF REDUCED WELL DEPTHS FOR EACH SPECIES PAIR.
!                  UNITS - K
!                  DIMENSION EOK(KDIM,*) EXACTLY KDIM FOR THE FIRST
!                    DIMENSION, AND AT LEAST KK FOR THE SECOND.
!    SGM     - MATRIX OF REDUCED COLLISION DIAMETERS FOR EACH SPECIES
!              PAIR.
!                  UNITS - ANGSTROMS.
!                  DIMENSION SGM(KDIM,*) EXACTLY KDIM FOR THE FIRST
!                    DIMENSION, AND AT LEAST KK FOR THE SECOND.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
  DIMENSION EPS(*), SIG(*), DIP(*), POL(*), EOK(KK,*), SGM(KK,*)
!
  DATA FDTCGS/1.0D-18/, FATCM/1.0D8/, &
       DIPMIN/1.0D-20/, BOLTZ/1.38056D-16/
!
  ONE = 1.0
!         COMPUTE AND STORE EPS/K AND SIGMA FOR ALL PAIRS
!
  DO J = 1, KK
!
     DO  K = 1, J
!
        IF((DIP(J).LT.DIPMIN .AND. DIP(K).GT.DIPMIN)) THEN
!
!                K IS POLAR, J IS NONPOLAR
!
           XI = ONE + 0.25*(POL(J)/SIG(J)**3) * &
                (FDTCGS**2*FATCM**3/BOLTZ) * &
                (DIP(K)**2/(EPS(K)*SIG(K)**3)) * &
                SQRT(EPS(K)/EPS(J))
           SGM(K,J) = 0.5 * (SIG(J)+SIG(K)) * XI**(-ONE/6.0)
           SGM(J,K) = SGM(K,J)
           EOK(K,J) = SQRT(EPS(J)*EPS(K)) * XI**2
           EOK(J,K) = EOK(K,J)
!
        ELSE IF((DIP(J).GT.DIPMIN .AND. DIP(K).LT.DIPMIN)) THEN
!
!             J IS POLAR, K IS NONPOLAR
!
           XI = ONE + 0.25*(POL(K)/SIG(K)**3) * &
                (FDTCGS**2*FATCM**3/BOLTZ) * &
                (DIP(J)**2/(EPS(J)*SIG(J)**3)) * &
                SQRT(EPS(J)/EPS(K))
           SGM(K,J) = 0.5 * (SIG(J)+SIG(K)) * XI**(-ONE / 6.0)
           SGM(J,K) = SGM(K,J)
           EOK(K,J) = SQRT(EPS(J)*EPS(K)) * XI**2
           EOK(J,K) = EOK(K,J)
!
        ELSE
!
!              NORMAL CASE, EITHER BOTH POLAR OR BOTH NONPOLAR
!
           SGM(K,J) = 0.5 * (SIG(J) + SIG(K))
           SGM(J,K) = SGM(K,J)
           EOK(K,J) = SQRT(EPS(J)*EPS(K))
           EOK(J,K) = EOK(K,J)
!
        ENDIF
     end DO
  end DO
!
  RETURN
END SUBROUTINE MCEPSG
!
!----------------------------------------------------------------------C
!
SUBROUTINE CKYTX  (Y, ICKWRK, RCKWRK, X)
!
!  START PROLOGUE
!
!  SUBROUTINE CKYTX  (Y, ICKWRK, RCKWRK, X)
!     Returns the mole fractions given the mass fractions;  see Eq. (6).
!
!  INPUT
!     Y      - Mass fractions of the species.
!                   cgs units - none
!                   Data type - real array
!                   Dimension Y(*) at least KK, the total number of
!                   species.
!     ICKWRK - Array of integer workspace.
!                   Data type - integer array
!                   Dimension ICKWRK(*) at least LENIWK.
!     RCKWRK - Array of real work space.
!                   Data type - real array
!                   Dimension RCKWRK(*) at least LENRWK.
!
!  OUTPUT
!     X      - Mole fractions of the species.
!                   cgs units - none
!                   Data type - real array
!                   Dimension X(*) at least KK, the total number of
!                   species.
!
!  END PROLOGUE
!
  IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
!
  INCLUDE 'ckstrt.fh'
!
  DIMENSION ICKWRK(*), RCKWRK(*), Y(*), X(*)
!
  SUMYOW = 0.0
  DO  K = 1, NKK
     SUMYOW = SUMYOW + Y(K)/RCKWRK(NcWT + K - 1)
  end DO
  DO  K = 1, NKK
     X(K) = Y(K)/(SUMYOW*RCKWRK(NcWT + K - 1))
  end DO
  RETURN
END SUBROUTINE CKYTX
!
!----------------------------------------------------------------------C
!
SUBROUTINE CKXTY  (X, ICKWRK, RCKWRK, Y)
!
!  START PROLOGUE
!
!  SUBROUTINE CKXTY  (X, ICKWRK, RCKWRK, Y)
!     Returns the mass fractions given the mole fractions;
!     see Eq. (9).
!
!  INPUT
!     X      - Mole fractions of the species.
!                   cgs units - none
!                   Data type - real array
!                   Dimension X(*) at least KK, the total number of
!                   species.
!     ICKWRK - Array of integer workspace.
!                   Data type - integer array
!                   Dimension ICKWRK(*) at least LENIWK.
!     RCKWRK - Array of real work space.
!                   Data type - real array
!                   Dimension RCKWRK(*) at least LENRWK.
!
!  OUTPUT
!     Y      - Mass fractions of the species.
!                   cgs units - none
!                   Data type - real array
!                   Dimension Y(*) at least KK, the total number of
!                   species.
!
!  END PROLOGUE
!
  IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
!
  INCLUDE 'ckstrt.fh'
!
  DIMENSION ICKWRK(*), RCKWRK(*), X(*), Y(*)
!
  SUM = 0.0d0
  DO K = 1, NKK
     SUM = SUM + X(K)*RCKWRK(NcWT + K - 1)
  end DO
!
  DO  K = 1, NKK
     Y(K) = X(K)*RCKWRK(NcWT + K - 1)/SUM
  end DO
  RETURN
END SUBROUTINE CKXTY

!
!----------------------------------------------------------------------C
!
SUBROUTINE CKXTCP (P, T, X, ICKWRK, RCKWRK, C)
!
!  START PROLOGUE
!
!  SUBROUTINE CKXTCP (P, T, X, ICKWRK, RCKWRK, C)
!     Returns the molar concentrations given the pressure,
!     temperature and mole fractions;  see Eq. (10).
!
!  INPUT
!     P      - Pressure.
!                   cgs units - dynes/cm**2
!                   Data type - real scalar
!     T      - Temperature.
!                   cgs units - K
!                   Data type - real scalar
!     X      - Mole fractions of the species.
!                   cgs units - none
!                   Data type - real array
!                   Dimension X(*) at least KK, the total number of
!                   species.
!     ICKWRK - Array of integer workspace.
!                   Data type - integer array
!                   Dimension ICKWRK(*) at least LENIWK.
!     RCKWRK - Array of real work space.
!                   Data type - real array
!                   Dimension RCKWRK(*) at least LENRWK.
!
!  OUTPUT
!     C      - Molar concentrations of the species.
!                   cgs units - mole/cm**3
!                   Data type - real array
!                   Dimension C(*) at least KK, the total number of
!                   species.
!
!  END PROLOGUE
!
  IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
!
  INCLUDE 'ckstrt.fh'
!
  DIMENSION ICKWRK(*), RCKWRK(*), X(*), C(*)
!
  PRUT = P/(RCKWRK(NcRU)*T)
  DO  K = 1, NKK
     C(K) = X(K)*PRUT
  END DO
  RETURN
END SUBROUTINE CKXTCP
!
!----------------------------------------------------------------------C
!
SUBROUTINE CKPC   (RHO, T, C, ICKWRK, RCKWRK, P)
!
!  START PROLOGUE
!
!  SUBROUTINE CKPC   (RHO, T, C, ICKWRK, RCKWRK, P)
!     Returns the pressure of the gas mixture given the mass density,
!     temperature and molar concentrations;  see Eq. (2).
!
!  INPUT
!     RHO    - Mass density.
!                   cgs units - gm/cm**3
!                   Data type - real scalar
!     T      - Temperature.
!                   cgs units - K
!                   Data type - real scalar
!     C      - Molar concentrations of the species.
!                   cgs unitsy - mole/cm**3
!                   Data type - real array
!                   Dimension C(*) at least KK, the total number of
!                   species.
!     ICKWRK - Array of integer workspace.
!                   Data type - integer array
!                   Dimension ICKWRK(*) at least LENIWK.
!     RCKWRK - Array of real work space.
!                   Data type - real array
!                   Dimension RCKWRK(*) at least LENRWK.
!
!  OUTPUT
!     P      - Pressure.
!                   cgs units - dynes/cm**2
!                   Data type - real scalar
!
!  END PROLOGUE
!
  IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
!
  INCLUDE 'ckstrt.fh'
  DIMENSION C(*), ICKWRK(*), RCKWRK(*)
!
  CTOT = 0d0
  SUM = 0d0
  DO  K = 1, NKK
     CTOT = CTOT + C(K)
     SUM  = SUM + C(K)*RCKWRK(NcWT + K - 1)
  END DO
  P    = RHO*RCKWRK(NcRU) * T * CTOT / SUM
  RETURN
END SUBROUTINE CKPC
!----------------------------------------------------------------------C
!
SUBROUTINE CKPY(RHO, T, Y, ICKWRK, RCKWRK, P)
!
!  START PROLOGUE
!
!  SUBROUTINE CKPY   (RHO, T, Y, ICKWRK, RCKWRK, P)
!     Returns the pressure of the gas mixture given the mass density,
!     temperature and mass fractions;  see Eq. (*).
!
!  INPUT
!     RHO    - Mass density.
!                   cgs units - gm/cm**3
!                   Data type - real scalar
!     T      - Temperature.
!                   cgs units - K
!                   Data type - real scalar
!     Y      - Mass fractions of the species.
!                   cgs units - none
!                   Data type - real array
!                   Dimension Y(*) at least KK, the total number of
!                   species.
!     ICKWRK - Array of integer workspace.
!                   Data type - integer array
!                   Dimension ICKWRK(*) at least LENIWK.
!     RCKWRK - Array of real work space.
!                   Data type - real array
!                   Dimension RCKWRK(*) at least LENRWK.
!
!  OUTPUT
!     P      - Pressure.
!                   cgs units - dynes/cm**2
!                   Data type - real scalar
!
!  END PROLOGUE
!
  IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
!
  DIMENSION Y(*), ICKWRK(*), RCKWRK(*)
!
  INCLUDE 'ckstrt.fh'
!
  SUMYOW = 0d0
  DO  K = 1, NKK
     SUMYOW = SUMYOW + Y(K)/RCKWRK(NcWT + K - 1)
  end DO
  P = RHO * RCKWRK(NcRU) * T * SUMYOW
  RETURN
END SUBROUTINE CKPY
!
!----------------------------------------------------------------------C
!
SUBROUTINE CKRHOY (P, T, Y, ICKWRK, RCKWRK, RHO)
!
!  START PROLOGUE
!
!  SUBROUTINE CKRHOY (P, T, Y, ICKWRK, RCKWRK, RHO)
!     Returns the mass density of the gas mixture given the pressure,
!     temperature and mass fractions;  see Eq. (2).
!
!  INPUT
!     P      - Pressure.
!                   cgs units - dynes/cm**2
!                   Data type - real scalar
!     T      - Temperature.
!                   cgs units - K
!                   Data type - real scalar
!     Y      - Mass fractions of the species.
!                   cgs units - none
!                   Data type - real array
!                   Dimension Y(*) at least KK, the total number of
!                   species.
!     ICKWRK - Array of integer workspace.
!                   Data type - integer array
!                   Dimension ICKWRK(*) at least LENIWK.
!     RCKWRK - Array of real work space.
!                   Data type - real array
!                   Dimension RCKWRK(*) at least LENRWK.
!
!  OUTPUT
!     RHO    - Mass density.
!                   cgs units - gm/cm**3
!                   Data type - real scalar
!
!  END PROLOGUE
!
  IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
!
  INCLUDE 'ckstrt.fh'
!
  DIMENSION ICKWRK(*), RCKWRK(*), Y(*)
!
  SUMYOW = 0.0
  DO  K = 1, NKK
     SUMYOW = SUMYOW + Y(K)/RCKWRK(NcWT + K - 1)
  END DO
  RHO = P/(SUMYOW*T*RCKWRK(NcRU))
  RETURN
END SUBROUTINE CKRHOY
!
!----------------------------------------------------------------------C
SUBROUTINE CKTEMPY (P, RHO, Y, ICKWRK, RCKWRK, T)
!
!  START PROLOGUE
! SUBROUTINE TO REUTRN TEMPERATURE FROM DENSITY, PRESSURE, MASS FRACTIONS

    IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
!
  INCLUDE 'ckstrt.fh'
!
  DIMENSION ICKWRK(*), RCKWRK(*), Y(*)
!
  SUMYOW = 0.0
  DO  K = 1, NKK
     SUMYOW = SUMYOW + Y(K)/RCKWRK(NcWT + K - 1)
  END DO
  T = P/(SUMYOW*RHO*RCKWRK(NcRU))
  RETURN
END SUBROUTINE CKTEMPY
!----------------------------------------------------------------------C
!
SUBROUTINE CKRgas (Y, ICKWRK, RCKWRK, Rgas)
!
!  START PROLOGUE
!
!  SUBROUTINE CKRgas (Y, ICKWRK, RCKWRK, Rgas)
!     Returns the Gas constant for the mixture
!
!  INPUT
!     Y      - Mass fractions of the species.
!                   cgs units - none
!                   Data type - real array
!                   Dimension Y(*) at least KK, the total number of
!                   species.
!     ICKWRK - Array of integer workspace.
!                   Data type - integer array
!                   Dimension ICKWRK(*) at least LENIWK.
!     RCKWRK - Array of real work space.
!                   Data type - real array
!                   Dimension RCKWRK(*) at least LENRWK.
!
!  OUTPUT
!     Rgas    - Gas constant
!                   cgs units - ergs/(gm*K)
!                   Data type - real scalar
!
!  END PROLOGUE
!
  IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
!
  INCLUDE 'ckstrt.fh'
!
  DIMENSION ICKWRK(*), RCKWRK(*), Y(*)
!
  SUMYOW = 0.0
  DO  K = 1, NKK
     SUMYOW = SUMYOW + Y(K)/RCKWRK(NcWT + K - 1)
  END DO
  RGas = SUMYOW*RCKWRK(NcRU)
!
  RETURN
END SUBROUTINE CKRgas

!
!----------------------------------------------------------------------C
!
SUBROUTINE CKRHOX (P, T, X, ICKWRK, RCKWRK, RHO)
!
!  START PROLOGUE
!
!  SUBROUTINE CKRHOX (P, T, X, ICKWRK, RCKWRK, RHO)
!     Returns the mass density of the gas mixture given the pressure,
!     temperature and mole fractions;  see Eq. (2).
!
!  INPUT
!     P      - Pressure.
!                   cgs units - dynes/cm**2
!                   Data type - real scalar
!     T      - Temperature.
!                   cgs units - K
!                   Data type - real scalar
!     X      - Mole fractions of the species.
!                   cgs units - none
!                   Data type - real array
!                   Dimension X(*) at least KK, the total number of
!                   species.
!     ICKWRK - Array of integer workspace.
!                   Data type - integer array
!                   Dimension ICKWRK(*) at least LENIWK.
!     RCKWRK - Array of real work space.
!                   Data type - real array
!                   Dimension RCKWRK(*) at least LENRWK.
!
!  OUTPUT
!     RHO    - Mass density.
!                   cgs units - gm/cm**3
!                   Data type - real scalar
!
!  END PROLOGUE
!
  IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
!
  INCLUDE 'ckstrt.fh'
  DIMENSION X(*), ICKWRK(*), RCKWRK(*)
!
!
  SUM = 0.0
  DO  K = 1, NKK
     SUM = SUM + X(K)*RCKWRK(NcWT + K - 1)
  end DO
!
  RHO = SUM * P / (RCKWRK(NcRU)*T)
  RETURN
END SUBROUTINE CKRHOX
!
!----------------------------------------------------------------------C
!
FUNCTION ILASCH   (STRING)
!   BEGIN PROLOGUE  ILASCH
!   DATE WRITTEN   850626
!   REVISION DATE  850626
!   CATEGORY NO.  M4.
!   KEYWORDS  CHARACTER STRINGS,SIGNIFICANT CHARACTERS
!   AUTHOR  CLARK,G.L.,GROUP C-3 LOS ALAMOS NATL LAB
!   PURPOSE  Determines last significant (non-blank) character
!            in character variable
!   DESCRIPTION
!
!-----------------------------------------------------------------------
!  IFIRCH locates the last non-blank character in a string of
!  arbitrary length.  If no characters are found, ILASCH is set = 0.
!  When used with the companion routine IFIRCH, the length of a string
!  can be determined, and/or a concatenated substring containing the
!  significant characters produced.
!  Note that the FORTRAN intrinsic function LEN returns the length
!  of a character string as declared, rather than as filled.  The
!  declared length includes leading and trailing blanks, and thus is
!  not useful in generating 'significant' substrings.
!-----------------------------------------------------------------------
!
!   REFERENCES  (NONE)
!   ROUTINES CALLED  (NONE)
!   END PROLOGUE IFIRCH
  IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
!
  CHARACTER*(*) STRING
!
!   FIRST EXECUTABLE STATEMENT ILASCH
  NLOOP = LEN(STRING)
  IF (NLOOP.EQ.0 .OR. STRING.EQ.' ') THEN
     ILASCH = 0
     RETURN
  ENDIF
!
  DO  I = NLOOP, 1, -1
     ILASCH = I
     IF (STRING(I:I) .NE. ' ') RETURN
  end DO
!
END FUNCTION ILASCH
!
!------------------------------------------------------------------
!
SUBROUTINE MCACON (T, X, RMCWRK, CONMIX)
!
  IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  SUBROUTINE MCACON (T, X, RMCWRK, CONMIX)
!
!    THIS SUBROUTINE COMPUTES THE MIXTURE THERMAL CONDUCTIVITY, GIVEN
!    THE TEMPERATURE AND THE SPECIES MOLE FRACTIONS.
!
!  INPUT-
!    T       - TEMPERATURE
!                  CGS UNITS - K.
!    X       - ARRAY OF MOLE FRACTIONS OF THE MIXTURE.
!                  DIMENSION X(*) AT LEAST KK.
!
!  WORK-
!    RMCWRK  - ARRAY OF FLOATING POINT STORAGE AND WORK SPACE.  THE
!              STARTING ADDRESSES FOR THE RMCWRK SPACE ARE STORED IN
!              COMMON /MCMCMC/.
!                  DIMENSION RMCWRK(*) AT LEAST LENRMC.
!
!  OUTPUT-
!    CONMIX  - MIXTURE THERMAL CONDUCTIVITY
!                  CGS UNITS - ERG/CM*K*S.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
  DIMENSION X(*), RMCWRK(*)
!
  COMMON /MCMCMC/ RU, PATMOS, SMALL, NKK, NO, NLITE, INLIN, IKTDIF,&
       IPVT, NWT, NEPS, NSIG, NDIP, NPOL, NZROT, NLAM,&
       NETA, NDIF, NTDIF, NXX, NVIS, NXI, NCP, NCROT,&
       NCINT, NBIND, NEOK, NSGM, NAST, NBST,&
       NCST, NXL, NR, NWRK, K3
!
  ONE  = 1.0d0
  ZERO = 0.0d0
!
!       IN THE FOLLOWING CALL:
!         THE PURE SPECIES CONDUCTIVITIES ARE IN RMCWRK(NXI)
!
  ALOGT = LOG(T)
  CALL MCEVAL (ALOGT, NKK, NO, RMCWRK(NLAM), RMCWRK(NXI))
!
  SUM = ZERO
  SUMR = ZERO
  DO K = 1, NKK
     RMCWRK(NXI+K-1) = EXP(RMCWRK(NXI+K-1))
     SUM =  SUM  + X(K)*RMCWRK(NXI+K-1)
     SUMR = SUMR + X(K)/RMCWRK(NXI+K-1)
  end DO
!
  CONMIX = 0.5d0 * (SUM + ONE/SUMR)
!
  RETURN
END SUBROUTINE MCACON
!
!---------------------------------------------------------------------
!
SUBROUTINE MCAVIS (T, X, RMCWRK, VISMIX)
!
  IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  SUBROUTINE MCAVIS (T, X, RMCWRK, VISMIX)
!
!    THIS SUBROUTINE COMPUTES THE MIXTURE VISCOSITY, GIVEN
!    THE TEMPERATURE AND THE SPECIES MOLE FRACTIONS.  IT USES MODIFICATI
!    OF THE WILKE SEMI-EMPIRICAL FORMULAS.
!
!  INPUT-
!    T       - TEMPERATURE
!                  CGS UNITS - K.
!    X       - ARRAY OF MOLE FRACTIONS OF THE MIXTURE.
!                  DIMENSION Y(*) AT LEAST KK.
!
!  WORK-
!    RMCWRK  - ARRAY OF FLOATING POINT STORAGE AND WORK SPACE.  THE
!              STARTING ADDRESSES FOR THE RMCWRK SPACE ARE STORED IN
!              COMMON /MCMCMC/.
!                  DIMENSION RMCWRK(*) AT LEAST LENRMC.
!
!  OUTPUT-
!    VISMIX  - MIXTURE VISCOSITY
!                  CGS UNITS - GM/CM*S.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
  DIMENSION X(*), RMCWRK(*)
!
  COMMON /MCMCMC/ RU, PATMOS, SMALL, NKK, NO, NLITE, INLIN, IKTDIF,&
       IPVT, NWT, NEPS, NSIG, NDIP, NPOL, NZROT, NLAM,&
       NETA, NDIF, NTDIF, NXX, NVIS, NXI, NCP, NCROT,&
       NCINT, NBIND, NEOK, NSGM, NAST, NBST,&
       NCST, NXL, NR, NWRK, K3
!
  SQRT8R = 0.3535534d0
  ZERO   = 0.0d0
  ONE    = 1.0d0
!
!       IN THE FOLLOWING CALL:
!         THE SPECIES MOLECULAR WEIGHTS ARE STORED IN RMCWRK(NWT)
!         THE PURE SPECIES VISCOSITIES ARE IN RMCWRK(NVIS)
!

  ALOGT = LOG(T)
  CALL MCEVAL (ALOGT, NKK, NO, RMCWRK(NETA), RMCWRK(NVIS))
  DO  K = 1, NKK
     RMCWRK(NVIS+K-1) = EXP(RMCWRK(NVIS+K-1))
  end DO
!
  SUMO = ZERO
  DO  K = 1, NKK
!
     SUMI = ZERO
     DO  J = 1, NKK
        TOP = (ONE + SQRT(RMCWRK(NVIS+K-1)/RMCWRK(NVIS+J-1)) * &
             (RMCWRK(NWT+J-1)/RMCWRK(NWT+K-1))**0.25D0 )**2
        BOT = SQRT(ONE + RMCWRK(NWT+K-1) / RMCWRK(NWT+J-1))
        PHIKJ = SQRT8R*TOP/BOT
        SUMI = SUMI + X(J)*PHIKJ
     end DO
     SUMO = SUMO + X(K)*RMCWRK(NVIS+K-1) / SUMI
  end DO
!
  VISMIX = SUMO
!
  RETURN
END SUBROUTINE MCAVIS
!
!---------------------------------------------------------------------
!
SUBROUTINE MCMODAVIS (T, X, RMCWRK, VISMIX)
!
  !Use DRIVTP,ONLY : MWT1JK,MWT2JK,SQRTVIS
  IMPLICIT NONE
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  SUBROUTINE MCMODAVIS (T, X, RMCWRK, VISMIX)
!
!    THIS SUBROUTINE COMPUTES THE MIXTURE VISCOSITY, GIVEN
!    THE TEMPERATURE AND THE SPECIES MOLE FRACTIONS.  IT USES MODIFICATI
!    OF THE WILKE SEMI-EMPIRICAL FORMULAS.
!
!  INPUT-
!    T       - TEMPERATURE
!                  CGS UNITS - K.
!    X       - ARRAY OF MOLE FRACTIONS OF THE MIXTURE.
!                  DIMENSION Y(*) AT LEAST KK.
!
!  WORK-
!    RMCWRK  - ARRAY OF FLOATING POINT STORAGE AND WORK SPACE.  THE
!              STARTING ADDRESSES FOR THE RMCWRK SPACE ARE STORED IN
!              COMMON /MCMCMC/.
!                  DIMENSION RMCWRK(*) AT LEAST LENRMC.
!
!  OUTPUT-
!    VISMIX  - MIXTURE VISCOSITY
!                  CGS UNITS - GM/CM*S.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
  INTEGER :: J,K
  REAL*8,INTENT(IN) :: T
  REAL*8,INTENT(OUT) :: VISMIX
  REAL*8 :: X(*), RMCWRK(*)
  REAL*8 :: ALOGT,SUMI,SUMO,ONE,PHIKJ
!
  real*8 :: RU, PATMOS, SMALL
  INTEGER :: NKK, NO, NLITE, INLIN, IKTDIF,&
       IPVT, NWT, NEPS, NSIG, NDIP, NPOL, NZROT, NLAM,&
       NETA, NDIF, NTDIF, NXX, NVIS, NXI, NCP, NCROT,&
       NCINT, NBIND, NEOK, NSGM, NAST, NBST,&
       NCST, NXL, NR, NWRK, K3
  COMMON /MCMCMC/ RU, PATMOS, SMALL, NKK, NO, NLITE, INLIN, IKTDIF,&
       IPVT, NWT, NEPS, NSIG, NDIP, NPOL, NZROT, NLAM,&
       NETA, NDIF, NTDIF, NXX, NVIS, NXI, NCP, NCROT,&
       NCINT, NBIND, NEOK, NSGM, NAST, NBST,&
       NCST, NXL, NR, NWRK, K3

!
!       IN THE FOLLOWING CALL:
!         THE SPECIES MOLECULAR WEIGHTS ARE STORED IN RMCWRK(NWT)
!         THE PURE SPECIES VISCOSITIES ARE IN RMCWRK(NVIS)
!
  ONE = 1d0
!
  !ALOGT = LOG(T)
  !CALL MCEVAL (ALOGT, NKK, NO, RMCWRK(NETA), RMCWRK(NVIS))
  !DO  K = 1, NKK
  !   RMCWRK(NVIS+K-1) = EXP(RMCWRK(NVIS+K-1))
  !   SQRTVIS(K) = SQRT(RMCWRK(NVIS+K-1))
  !end DO
!
  !SUMO = 0d0
  !DO  K = 1, NKK
!
   !  SUMI = 0d0
   !  DO  J = 1, NKK
   !     PHIKJ = (ONE +  SQRTVIS(K)/SQRTVIS(J)*MWT1JK(J,K) )**2
   !     SUMI = SUMI + X(J)*PHIKJ*MWT2JK(J,K)
   !  end DO
   !  SUMO = SUMO + X(K)*RMCWRK(NVIS+K-1) / SUMI
  !end DO
!
  VISMIX = SUMO
!
  RETURN
END SUBROUTINE MCMODAVIS
!
!---------------------------------------------------------------------
!
SUBROUTINE MCEVAL (TF, KK, NO, COF, VAL)
!
  IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  SUBROUTINE MCEVAL (TF, KK, NO, COF, VAL)
!
!    THIS SUBROUTINE USES HORNERS ALGORITHM TO EVALUATE A POLYNOMIAL
!    FIT.  THIS ROUTINE IS NOT NORMALLY CALLED BY THE PACKAGE USER.
!
!  INPUT-
!    TF      - INDEPENDENT VARIABLE OF FIT. EITHER TEMPERATURE
!              OR LOG TEMPERATURE.
!    KK      - NUMBER OF SPECIES.
!    NO      - ORDER OF FIT.
!    COF     - MATRIX OF FIT COEFFICIENTS.  COF(N,K) IS THE NTH
!              COEFFICIENT OF A FIT FOR KTH SPECIES PROPERTY.
!                 DIMENSION COF(NO,*) EXACTLY NO FOR THE FIRST
!                 DIMENSION AND AT LEAST KK FOR THE SECOND.
!
!  OUTPUT-
!    VAL     - ARRAY OF VALUES, EVALUATED FROM THE FIT AT TF.
!                 DIMENSION VAL(*) AT KEAST KK.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
  DIMENSION COF(NO,*), VAL(*)
!
  NOM1 = NO-1
!
  DO  K = 1, KK
     B = COF(NO,K)
     DO  I = 1, NOM1
        B = COF(NO-I,K) + B*TF
     END DO
     VAL(K) = B
  END DO
!
  RETURN
END SUBROUTINE MCEVAL

!
REAL*8 FUNCTION DELTA (I,K)
  IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
!
! STATEMENT FUNCTION FOR DELTA FUNCTION
!
  DELTA = FLOAT( MAX(-1, -IABS(I-K)) + 1 )
  RETURN
END FUNCTION DELTA

!REAL*8 FUNCTION DDOT(N,DX,INCX,DY,INCY)
!***BEGIN PROLOGUE  DDOT
!***DATE WRITTEN   791001   (YYMMDD)
!***REVISION DATE  861211   (YYMMDD)
!***CATEGORY NO.  D1A4
!***KEYWORDS  LIBRARY=SLATEC(BLAS),
!             TYPE=REAL*8(SDOT-S DDOT-D CDOTU-C),
!             INNER PRODUCT,LINEAR ALGEBRA,VECTOR
!***AUTHOR  LAWSON, C. L., (JPL)
!           HANSON, R. J., (SNLA)
!           KINCAID, D. R., (U. OF TEXAS)
!           KROGH, F. T., (JPL)
!***PURPOSE  D.P. inner product of d.p. vectors
!***DESCRIPTION
!
!                B L A S  Subprogram
!    Description of Parameters
!
!     --Input--
!        N  number of elements in input vector(s)
!       DX  real*8 vector with N elements
!     INCX  storage spacing between elements of DX
!       DY  real*8 vector with N elements
!     INCY  storage spacing between elements of DY
!
!     --Output--
!     DDOT  real*8 dot product (zero if N .LE. 0)
!
!     Returns the dot product of real*8 DX and DY.
!     DDOT = sum for I = 0 to N-1 of  DX(LX+I*INCX) * DY(LY+I*INCY)
!     where LX = 1 if INCX .GE. 0, else LX = (-INCX)*N, and LY is
!     defined in a similar way using INCY.
!***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
!                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
!                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
!                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  DDOT
!
!  REAL*8 :: DX(*),DY(*)
!***FIRST EXECUTABLE STATEMENT  DDOT
!  DDOT = 0.D0
!  IF(N.LE.0)RETURN
!  IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60
!5 CONTINUE
!
!         CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.
!
!  IX = 1
!  IY = 1
!  IF(INCX.LT.0)IX = (-N+1)*INCX + 1
!  IF(INCY.LT.0)IY = (-N+1)*INCY + 1
!  DO  I = 1,N
!     DDOT = DDOT + DX(IX)*DY(IY)
!     IX = IX + INCX
!     IY = IY + INCY
!  end DO
!  RETURN
!
!        CODE FOR BOTH INCREMENTS EQUAL TO 1.
!
!
!        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5.
!
!20 M = MOD(N,5)
!  IF( M .EQ. 0 ) GO TO 40
!  DO  I = 1,M
!     DDOT = DDOT + DX(I)*DY(I)
!  end DO
!  IF( N .LT. 5 ) RETURN
!40 MP1 = M + 1
!  DO  I = MP1,N,5
!     DDOT = DDOT + DX(I)*DY(I) + DX(I+1)*DY(I+1) +&
!          DX(I + 2)*DY(I + 2) + DX(I + 3)*DY(I + 3) + DX(I + 4)*DY(I + 4)
!  end DO
!  RETURN
!
!         CODE FOR POSITIVE EQUAL INCREMENTS .NE.1.
!
!60 CONTINUE
!  NS = N*INCX
!  DO  I=1,NS,INCX
!     DDOT = DDOT + DX(I)*DY(I)
!  end DO
!  RETURN
!END FUNCTION DDOT

!
!---------------------------------------------------------------------
!
SUBROUTINE MCMDIF (P, T, X, KDIM, IMCWRK, RMCWRK, D)
!
  IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  SUBROUTINE MCMDIF (P, T, X, KDIM, IMCWRK, RMCWRK, D)
!
!    THIS SUBROUTINE COMPUTES THE ORDINARY MULTICOMPONENT DIFFUSION
!    COEFFICIENTS, GIVEN THE PRESSURE, TEMPERATURE, AND MOLE FRACTIONS.
!
!  INPUT-
!    P       - PRESSURE
!                  CGS UNITS - DYNES/CM**2.
!    T       - TEMPERATURE
!                  CGS UNITS - K.
!    X       - ARRAY OF MOLE FRACTIONS OF THE MIXTURE.
!                  DIMENSION X(*) AT LEAST KK.
!    KDIM    - ACTUAL FIRST DIMENSION OF D(KDIM,KK).  KDIM MUST BE AT
!                LEAST THE NUMBER OF SPECIES, KK.
!
!  WORK-
!    IMCWRK  - ARRAY OF INTEGER STORAGE AND WORK SPACE.  THE
!              STARTING ADDRESSES FOR THE IMCWRK SPACE ARE STORED IN
!              COMMON /MCMCMC/.
!                  DIMENSION IMCWRK(*) AT LEAST LENIMC.
!    RMCWRK  - ARRAY OF FLOATING POINT STORAGE AND WORK SPACE.  THE
!              STARTING ADDRESSES FOR THE RMCWRK SPACE ARE STORED IN
!              COMMON /MCMCMC/.
!                  DIMENSION RMCWRK(*) AT LEAST LENRMC.
!
!  OUTPUT-
!    D       - MATRIX OF ORDINARY MULTICOMPONENT DIFFUSION COEFFICIENTS.
!                  CGS UNITS - CM**2/S
!                  DIMENSION DJK(KDIM,*) EXACTLY KDIM FOR THE FIRST
!                   DIMENSION AND AT LEAST KK FOR THE SECOND.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
  DIMENSION X(*), IMCWRK(*), RMCWRK(*), D(KDIM,*)
!
  COMMON /MCMCMC/ RU, PATMOS, SMALL, NKK, NO, NLITE, INLIN, IKTDIF,&
       IPVT, NWT, NEPS, NSIG, NDIP, NPOL, NZROT, NLAM,&
       NETA, NDIF, NTDIF, NXX, NVIS, NXI, NCP, NCROT,&
       NCINT, NBIND, NEOK, NSGM, NAST, NBST,&
       NCST, NXL, NR, NWRK, K3
!
  CALL MCORDF (P, T, X, NKK, KDIM, SMALL, RMCWRK(NWT), RMCWRK,&
       RMCWRK(NXX), RMCWRK(NBIND), RMCWRK(NXL),&
       RMCWRK(NWRK), IMCWRK(IPVT), D)
!
  RETURN
END SUBROUTINE MCMDIF
!
!---------------------------------------------------------------------
!
SUBROUTINE MCSDIF (P, T, KDIM, RMCWRK, DJK)
!
  IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  SUBROUTINE MCSDIF (P, T, KDIM, RMCWRK, DJK)
!
!    THIS SUBROUTINE COMPUTES THE BINARY DIFFUSION COEFFICIENTS, GIVEN
!    THE PRESSURE AND TEMPERATURE.
!
!  INPUT-
!    P       - PRESSURE
!                  CGS UNITS - DYNES/CM**2.
!    T       - TEMPERATURE
!                  CGS UNITS - K.
!    KDIM    - ACTUAL FIRST DIMENSION OF DJK(KDIM,KK)
!
!  WORK-
!    RMCWRK   - ARRAY OF FLOATING POINT STORAGE AND WORK SPACE.  THE
!              STARTING ADDRESSES FOR THE RMCWRK SPACE ARE STORED IN
!              COMMON /MCMCMC/.
!                  DIMENSION RMCWRK(*) AT LEAST LENRMC.
!
!  OUTPUT-
!    DJK     - MATRIX OF BINARY DIFFUSION COEFFICIENTS.  DJK(J,K) IS
!              DIFFUSION COEFFICIENT OF SPECIES J IN SPECIES K.
!                  CGS UNITS - CM**2/S
!                  DIMENSION DJK(KDIM,*) EXACTLY KDIM FOR THE FIRST
!                   DIMENSION AND AT LEAST KK FOR THE SECOND.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
  DIMENSION RMCWRK(*), DJK(KDIM,*)
!
  COMMON /MCMCMC/ RU, PATMOS, SMALL, NKK, NO, NLITE, INLIN, IKTDIF,&
       IPVT, NWT, NEPS, NSIG, NDIP, NPOL, NZROT, NLAM,&
       NETA, NDIF, NTDIF, NXX, NVIS, NXI, NCP, NCROT,&
       NCINT, NBIND, NEOK, NSGM, NAST, NBST,&
       NCST, NXL, NR, NWRK, K3
!
  PFAC = PATMOS/P
  ALOGT = LOG(T)
!
  DO  K = 1, NKK
     ISTART = NDIF + (K-1)*NO*NKK
     CALL MCEVAL (ALOGT, NKK, NO, RMCWRK(ISTART), DJK(1,K))
     DO  J = 1, NKK
        DJK(J,K) = EXP(DJK(J,K)) * PFAC
     end DO
  end DO
!
  RETURN
END SUBROUTINE MCSDIF
!
!---------------------------------------------------------------------
!
SUBROUTINE MCORDF (P, T, X, KK, KDIM, SMALL, WT, RMCWRK, XX,&
     BINDIF, XL0000, WORK, IPVT, D)
!
  IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  SUBROUTINE MCORDF (P, T, X, KK, KDIM, SMALL, WT, RMCWRK, XX,
! 1                   BINDIF, XL0000, WORK, IPVT, D)
!
!    THIS SUBROUTINE COMPUTES ORDINARY MULTICOMPONENT DIFFUSION COEFFICI
!    COEFFICIENT MATRIX.  IT DOES SO BY COMPUTING THE INVERSE OF THE
!    L00,00 MATRIX.  THIS ROUTINE IS NOT NORMALLY CALLED DIRECTLY BY THE
!    USER; THE USER CALLS MCMDIF, WHICH IN TURN CALLS MCORDF.
!
!  INPUT-
!    P       - PRESSURE
!                  CGS UNITS - DYNES/CM**2.
!    T       - TEMPERATURE
!                  CGS UNITS - K.
!    X       - ARRAY OF MOLE FRACTIONS OF THE MIXTURE.
!                  DIMENSION X(*) AT LEAST KK.
!    KK      - NUMBER OF SPECIES.
!    KDIM    - ACTUAL FIRST DIMENSION OF D(KDIM,KK).  KDIM MUST BE AT
!                LEAST THE NUMBER OF SPECIES, KK.
!    SMALL   - THE MOLE FRACTIONS USED IN THE TRANSPORT COMPUTATION
!              ARE GIVEN BY XX(K) = X(K) + SMALL.
!    WT      - THE ARRAY OF SPECIES MOLECULAR WEIGHTS.
!                   DIMENSION WT(*) AT LEAST KK.
!
!  WORK AND SCRATCH SPACE
!    RMCWRK  - ARRAY OF FLOATING POINT STORAGE AND WORK SPACE.  THE
!              STARTING ADDRESSES FOR THE RMCWRK SPACE ARE STORED IN
!              COMMON /MCMCMC/.
!                  DIMENSION RMCWRK(*) AT LEAST LENRMC.
!    XX      - THE SPECIES MOLE FRACTION ARRAY THAT IS USED IN THE
!              TRANSPORT COMPUTATION.  XX(K) = X(K) + SMALL.
!    BINDIF  - MATRIX OF BINARY DIFFUSION COEFFICIENTS.
!                  CGS UNITS - CM**2/S
!                  DIMENSION BINDIF(KDIM,*) EXACTLY KDIM FOR THE FIRST
!                   DIMENSION AND AT LEAST KK FOR THE SECOND.
!    XL0000  - THE L00,00 MATRIX.
!                  DIMENSION L0000(KDIM,*) EXACTLY KDIM FOR THE FIRST
!                   DIMENSION AND AT LEAST KK FOR THE SECOND.
!    WORK    - ARRAY OF WORK SPACE FOR THE INVERSION OF THE L00,00
!              MATRIX BY THE LINPACK ROUTINES SGEFA AND SGEDI.
!                  DIMENSION WORK(*) AT LEAST KK.
!    IPVT    - ARRAY OF PIVOT INDICES FOR THE INVERSION OF THE L00,00
!              MATRIX BY THE LINPACK ROUTINES SGEFA AND SGEDI.
!                  DIMENSION IPVT(*) AT LEAST KK.
!
!  OUTPUT-
!    D       - MATRIX OF ORDINARY MULTICOMPONENT DIFFUSION COEFFICIENTS.
!                  CGS UNITS - CM**2/S
!                  DIMENSION DJK(KDIM,*) EXACTLY KDIM FOR THE FIRST
!                   DIMENSION AND AT LEAST KK FOR THE SECOND.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
  DIMENSION X(*), WT(*), BINDIF(KK,*), XL0000(KK,*), IPVT(*),&
       DET(2), WORK(*), XX(*), RMCWRK(*), D(KDIM,*)
!
  DATA JOB/1/
!
  ONE = 1.0
  ZERO = 0.0
!
!         SET MINIMUM MOLE FRACTION TO SMALL
!
  DO  K = 1, KK
     XX(K) = X(K) + SMALL
  end DO
!
!        EVALUATE THE BINARY DIFFUSION COEFFICIENTS
!
  CALL MCSDIF (P, T, KK, RMCWRK, BINDIF)
!
!        ASSEMBLE L00,00
!
  DO  J = 1, KK
     DO I = 1, KK
        SUM = ZERO
        DO K = 1, KK
           SUM = SUM + XX(K) / (WT(I)*BINDIF(I,K)) *&
                ( WT(J)*XX(J) * (ONE - DELTA(I,K)) -&
                WT(I)*XX(I) * (DELTA(I,J) - DELTA(J,K)) )
        end DO
        XL0000(I,J) = 16.0 * T * SUM / (25.0 * P)
     end DO
  end DO
!
!        INVERT L00,00 USING LINPACK
  CALL DGEFA (XL0000, KK, KK, IPVT, INFO)
!
  IF (INFO .NE. 0) THEN
     WRITE (6, *) ' ERROR IN DGEFA, INFO = ', INFO
     STOP
  ENDIF
!
  CALL DGEDI (XL0000, KK, KK, IPVT, DET, WORK, JOB)
!
!
!        COMPUTE THE ORDINARY MULTICOMPONENT DIFFUSION COEFFICIENTS
!
  DO  J = 1, KK
     DO  I = 1, KK
        SUM = ZERO
        DO  K = 1, KK
           SUM = SUM + WT(K) * X(K)
        end DO
        D(I,J) = XX(I) * 16.0 * T * SUM&
             * (XL0000(I,J)-XL0000(I,I)) /(25.0 * P * WT(J) )
     end DO
  end DO
!
  RETURN
END SUBROUTINE MCORDF
!
!---------------------------------------------------------------------
!

SUBROUTINE CKDGECO(A,LDA,N,IPVT,RCOND,Z)
!***BEGIN PROLOGUE  DGECO
!***DATE WRITTEN   780814   (YYMMDD)
!***REVISION DATE  861211   (YYMMDD)
!***CATEGORY NO.  D2A1
!***KEYWORDS  LIBRARY=SLATEC(LINPACK),
!             TYPE=REAL*8 ::(SGECO-S DGECO-D CGECO-C),
!             CONDITION NUMBER,GENERAL MATRIX,LINEAR ALGEBRA,MATRIX,
!             MATRIX FACTORIZATION
!***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
!***PURPOSE  Factors a real*8 :: matrix by Gaussian elimination
!            and estimates the condition of the matrix.
!***DESCRIPTION
!
!     DGECO factors a real*8 :: matrix by Gaussian elimination
!     and estimates the condition of the matrix.
!
!     If  RCOND  is not needed, DGEFA is slightly faster.
!     To solve  A*X = B , follow DGECO by DGESL.
!     To compute  INVERSE(A)*C , follow DGECO by DGESL.
!     To compute  DETERMINANT(A) , follow DGECO by DGEDI.
!     To compute  INVERSE(A) , follow DGECO by DGEDI.
!
!     On Entry
!
!        A       REAL*8 ::(LDA, N)
!                the matrix to be factored.
!
!        LDA     INTEGER
!                the leading dimension of the array  A .
!
!        N       INTEGER
!                the order of the matrix  A .
!
!     On Return
!
!        A       an upper triangular matrix and the multipliers
!                which were used to obtain it.
!                The factorization can be written  A = L*U  where
!                L  is a product of permutation and unit lower
!                triangular matrices and  U  is upper triangular.
!
!        IPVT    INTEGER(N)
!                an INTEGER vector of pivot indices.
!
!        RCOND   REAL*8 ::
!                an estimate of the reciprocal condition of  A .
!                For the system  A*X = B , relative perturbations
!                in  A  and  B  of size  EPSILON  may cause
!                relative perturbations in  X  of size  EPSILON/RCOND .
!                If  RCOND  is so small that the logical expression
!                           1.0 + RCOND .EQ. 1.0
!                is true, then  A  may be singular to working
!                precision.  In particular,  RCOND  is zero  if
!                exact singularity is detected or the estimate
!                underflows.
!
!        Z       REAL*8 ::(N)
!                a work vector whose contents are usually unimportant.
!                If  A  is close to a singular matrix, then  Z  is
!                an approximate null vector in the sense that
!                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
!
!     LINPACK.  This version dated 08/14/78 .
!     Cleve Moler, University of New Mexico, Argonne National Lab.
!
!     Subroutines and Functions
!
!     LINPACK DGEFA
!     BLAS DAXPY,DDOT,DSCAL,DASUM
!     Fortran DABS,DMAX1,DSIGN
!***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
!                 *LINPACK USERS  GUIDE*, SIAM, 1979.
!***ROUTINES CALLED  DASUM,DAXPY,DDOT,DGEFA,DSCAL
!***END PROLOGUE  DGECO
  INTEGER LDA,N,IPVT(1)
  REAL*8 :: A(LDA,1),Z(1)
  REAL*8 :: RCOND
!
  REAL*8 :: DDOT,EK,T,WK,WKM
  REAL*8 :: ANORM,S,SM,YNORM,BNORM
  INTEGER INFO,J,K,KB,KP1,L
!
!     COMPUTE 1-NORM OF A
!
!***FIRST EXECUTABLE STATEMENT  DGECO
  ANORM = 0.0D0
  DO  J = 1, N
     BNORM = CKDASUM(N,A(1,J),1)
     ANORM = MAX(ANORM,BNORM)
  end DO
!
!     FACTOR
!
  CALL DGEFA(A,LDA,N,IPVT,INFO)
!
!     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
!     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .
!     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE
!     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE
!     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID
!     OVERFLOW.
!
!     SOLVE TRANS(U)*W = E
!
  EK = 1.0D0
  DO  J = 1, N
     Z(J) = 0.0D0
  end DO
  DO  K = 1, N
     IF (Z(K) .NE. 0.0D0) EK = DSIGN(EK,-Z(K))
     IF (DABS(EK-Z(K)) .LE. DABS(A(K,K))) GO TO 30
     S = DABS(A(K,K))/DABS(EK-Z(K))
     CALL DSCAL(N,S,Z,1)
     EK = S*EK
30   CONTINUE
     WK = EK - Z(K)
     WKM = -EK - Z(K)
     S = DABS(WK)
     SM = DABS(WKM)
     IF (A(K,K) .EQ. 0.0D0) GO TO 40
     WK = WK/A(K,K)
     WKM = WKM/A(K,K)
     GO TO 50
40   CONTINUE
     WK = 1.0D0
     WKM = 1.0D0
50   CONTINUE
     KP1 = K + 1
     IF (KP1 .GT. N) GO TO 90
     DO  J = KP1, N
        SM = SM + DABS(Z(J)+WKM*A(K,J))
        Z(J) = Z(J) + WK*A(K,J)
        S = S + DABS(Z(J))
     end DO
     IF (S .GE. SM) GO TO 80
     T = WKM - WK
     WK = WKM
     DO  J = KP1, N
        Z(J) = Z(J) + T*A(K,J)
     end DO
80   CONTINUE
90   CONTINUE
     Z(K) = WK
  end DO
  S = 1.0D0/CKDASUM(N,Z,1)
  CALL DSCAL(N,S,Z,1)
!
!     SOLVE TRANS(L)*Y = W
!
  DO  KB = 1, N
     K = N + 1 - KB
     IF (K .LT. N) Z(K) = Z(K) + DDOT(N-K,A(K+1,K),1,Z(K+1),1)
     IF (DABS(Z(K)) .LE. 1.0D0) GO TO 110
     S = 1.0D0/DABS(Z(K))
     CALL DSCAL(N,S,Z,1)
110  CONTINUE
     L = IPVT(K)
     T = Z(L)
     Z(L) = Z(K)
     Z(K) = T
  end DO
  S = 1.0D0/CKDASUM(N,Z,1)
  CALL DSCAL(N,S,Z,1)
!
  YNORM = 1.0D0
!
!     SOLVE L*V = Y
!
  DO  K = 1, N
     L = IPVT(K)
     T = Z(L)
     Z(L) = Z(K)
     Z(K) = T
     IF (K .LT. N) CALL DAXPY(N-K,T,A(K+1,K),1,Z(K+1),1)
     IF (DABS(Z(K)) .LE. 1.0D0) GO TO 130
     S = 1.0D0/DABS(Z(K))
     CALL DSCAL(N,S,Z,1)
     YNORM = S*YNORM
130  CONTINUE
  end DO
  S = 1.0D0/CKDASUM(N,Z,1)
  CALL DSCAL(N,S,Z,1)
  YNORM = S*YNORM
!
!     SOLVE  U*Z = V
!
  DO  KB = 1, N
     K = N + 1 - KB
     IF (DABS(Z(K)) .LE. DABS(A(K,K))) GO TO 150
     S = DABS(A(K,K))/DABS(Z(K))
     CALL DSCAL(N,S,Z,1)
     YNORM = S*YNORM
150  CONTINUE
     IF (A(K,K) .NE. 0.0D0) Z(K) = Z(K)/A(K,K)
     IF (A(K,K) .EQ. 0.0D0) Z(K) = 1.0D0
     T = -Z(K)
     CALL DAXPY(K-1,T,A(1,K),1,Z(1),1)
  end DO
!     MAKE ZNORM = 1.0
  S = 1.0D0/CKDASUM(N,Z,1)
  CALL DSCAL(N,S,Z,1)
  YNORM = S*YNORM
!
  IF (ANORM .NE. 0.0D0) RCOND = YNORM/ANORM
  IF (ANORM .EQ. 0.0D0) RCOND = 0.0D0
  RETURN
END SUBROUTINE CKDGECO

  REAL*8  FUNCTION CKDASUM(N,DX,INCX)
!***BEGIN PROLOGUE  CKDASUM
!***DATE WRITTEN   791001   (YYMMDD)
!***REVISION DATE  861211   (YYMMDD)
!***CATEGORY NO.  D1A3A
!***KEYWORDS  LIBRARY=SLATEC(BLAS),
!             TYPE=REAL*8 ::(SASUM-S DASUM-D SCASUM-C),ADD,
!             LINEAR ALGEBRA,MAGNITUDE,SUM,VECTOR
!***AUTHOR  LAWSON, C. L., (JPL)
!           HANSON, R. J., (SNLA)
!           KINCAID, D. R., (U. OF TEXAS)
!           KROGH, F. T., (JPL)
!***PURPOSE  Sum of magnitudes of d.p. vector components
!***DESCRIPTION
!
!                B L A S  Subprogram
!    Description of Parameters
!
!     --Input--
!        N  number of elements in input vector(s)
!       DX  real*8 :: vector with N elements
!     INCX  storage spacing between elements of DX
!
!     --Output--
!    DASUM  real*8 :: result (zero if N .LE. 0)
!
!     Returns sum of magnitudes of real*8 :: DX.
!     DASUM = sum from 0 to N-1 of DABS(DX(1+I*INCX))
!***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
!                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
!                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
!                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  DASUM
!
    REAL*8 :: DX(1)
!***FIRST EXECUTABLE STATEMENT  DASUM
    CKDASUM = 0.D0
    IF(N.LE.0)RETURN
    IF(INCX.EQ.1)GOTO 20
!
!        CODE FOR INCREMENTS NOT EQUAL TO 1.
!
    NS = N*INCX
    DO  I=1,NS,INCX
       CKDASUM = CKDASUM + DABS(DX(I))
    end DO
    RETURN
!
!        CODE FOR INCREMENTS EQUAL TO 1.
!
!
!        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 6.
!
20  M = MOD(N,6)
    IF( M .EQ. 0 ) GO TO 40
    DO  I = 1,M
       CKDASUM = CKDASUM + DABS(DX(I))
    end DO
    IF( N .LT. 6 ) RETURN
40  MP1 = M + 1
    DO  I = MP1,N,6
       CKCKDASUM = CKDASUM + DABS(DX(I)) + DABS(DX(I+1)) + DABS(DX(I+2)) &
            &+ DABS(DX(I+3)) + DABS(DX(I+4)) + DABS(DX(I+5))
    end DO
    RETURN
  END FUNCTION CKDASUM

!SUBROUTINE DGEFA(A,LDA,N,IPVT,INFO)
!***BEGIN PROLOGUE  DGEFA
!***DATE WRITTEN   780814   (YYMMDD)
!***REVISION DATE  861211   (YYMMDD)
!***CATEGORY NO.  D2A1
!***KEYWORDS  LIBRARY=SLATEC(LINPACK),
!             TYPE=REAL*8(SGEFA-S DGEFA-D CGEFA-C),
!             GENERAL MATRIX,LINEAR ALGEBRA,MATRIX,MATRIX FACTORIZATION
!***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
!***PURPOSE  Factors a real*8 matrix by Gaussian elimination.
!***DESCRIPTION
!
!     DGEFA factors a real*8 matrix by Gaussian elimination.
!
!     DGEFA is usually called by DGECO, but it can be called
!     directly with a saving in time if  RCOND  is not needed.
!     (Time for DGECO) = (1 + 9/N)*(Time for DGEFA) .
!
!     On Entry
!
!        A       REAL*8(LDA, N)
!                the matrix to be factored.
!
!        LDA     INTEGER
!                the leading dimension of the array  A .
!
!        N       INTEGER
!                the order of the matrix  A .
!
!     On Return
!
!        A       an upper triangular matrix and the multipliers
!                which were used to obtain it.
!                The factorization can be written  A = L*U  where
!                L  is a product of permutation and unit lower
!                triangular matrices and  U  is upper triangular.
!
!        IPVT    INTEGER(N)
!                an integer vector of pivot indices.
!
!        INFO    INTEGER
!                = 0  normal value.
!                = K  if  U(K,K) .EQ. 0.0 .  This is not an error
!                     condition for this subroutine, but it does
!                     indicate that DGESL or DGEDI will divide by zero
!                     if called.  Use  RCOND  in DGECO for a reliable
!                     indication of singularity.
!
!     LINPACK.  This version dated 08/14/78 .
!     Cleve Moler, University of New Mexico, Argonne National Lab.
!
!     Subroutines and Functions
!
!     BLAS DAXPY,DSCAL,IDAMAX
!***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
!                 *LINPACK USERS  GUIDE*, SIAM, 1979.
!***ROUTINES CALLED  DAXPY,DSCAL,IDAMAX
!***END PROLOGUE  DGEFA
!  INTEGER :: LDA,N,IPVT(*),INFO
!  REAL*8 :: A(LDA,*)
!
!  REAL*8 :: T
!  INTEGER ::IDAMAX,J,K,KP1,L,NM1
!
!     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
!
!***FIRST EXECUTABLE STATEMENT  DGEFA
!  INFO = 0
!  NM1 = N - 1
!  IF (NM1 .LT. 1) GO TO 70
!  DO  K = 1, NM1
!     KP1 = K + 1
!
!        FIND L = PIVOT INDEX
!
!     L = IDAMAX(N-K+1,A(K,K),1) + K - 1
!     IPVT(K) = L
!
!        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
!
!     IF (A(L,K) .EQ. 0.0D0) GO TO 40
!
!           INTERCHANGE IF NECESSARY
!
!     IF (L .EQ. K) GO TO 10
!     T = A(L,K)
!     A(L,K) = A(K,K)
!     A(K,K) = T
!10   CONTINUE
!
!           COMPUTE MULTIPLIERS
!
!     T = -1.0D0/A(K,K)
!     CALL DSCAL(N-K,T,A(K+1,K),1)
!
!           ROW ELIMINATION WITH COLUMN INDEXING
!
!     DO  J = KP1, N
!        T = A(L,J)
!        IF (L .EQ. K) GO TO 20
!        A(L,J) = A(K,J)
!        A(K,J) = T
!20      CONTINUE
!        CALL DAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)
!     end DO
!     GO TO 50
!40   CONTINUE
!     INFO = K
!50   CONTINUE
!  end DO
!70 CONTINUE
!  IPVT(N) = N
!  IF (A(N,N) .EQ. 0.0D0) INFO = N
!  RETURN
!END SUBROUTINE DGEFA
!
!-------------------------------------------------------------------
!

SUBROUTINE DGEDI(A,LDA,N,IPVT,DET,WORK,JOB)
!***BEGIN PROLOGUE  DGEDI
!***DATE WRITTEN   780814   (YYMMDD)
!***REVISION DATE  861211   (YYMMDD)
!***CATEGORY NO.  D3A1,D2A1
!***KEYWORDS  LIBRARY=SLATEC(LINPACK),
!             TYPE=REAL*8(SGEDI-S DGEDI-D CGEDI-C),
!             DETERMINANT,INVERSE,LINEAR ALGEBRA,MATRIX,
!             MATRIX FACTORIZATION
!***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
!***PURPOSE  Computes the determinant and inverse of a matrix using
!            factors computed by DGECO or DGEFA.
!***DESCRIPTION
!
!     DGEDI computes the determinant and inverse of a matrix
!     using the factors computed by DGECO or DGEFA.
!
!     On Entry
!
!        A       REAL*8(LDA, N)
!                the output from DGECO or DGEFA.
!
!        LDA     INTEGER
!                the leading dimension of the array  A .
!
!        N       INTEGER
!                the order of the matrix  A .
!
!        IPVT    INTEGER(N)
!                the pivot vector from DGECO or DGEFA.
!
!        WORK    REAL*8(N)
!                work vector.  Contents destroyed.
!
!        JOB     INTEGER
!                = 11   both determinant and inverse.
!                = 01   inverse only.
!                = 10   determinant only.
!
!     On Return
!
!        A       inverse of original matrix if requested.
!                Otherwise unchanged.
!
!        DET     REAL*8(2)
!                determinant of original matrix if requested.
!                Otherwise not referenced.
!                Determinant = DET(1) * 10.0**DET(2)
!                with  1.0 .LE. DABS(DET(1)) .LT. 10.0
!                or  DET(1) .EQ. 0.0 .
!
!     Error Condition
!
!        A division by zero will occur if the input factor contains
!        a zero on the diagonal and the inverse is requested.
!        It will not occur if the subroutines are called correctly
!        and if DGECO has set RCOND .GT. 0.0 or DGEFA has set
!        INFO .EQ. 0 .
!
!     LINPACK.  This version dated 08/14/78 .
!     Cleve Moler, University of New Mexico, Argonne National Lab.
!
!     Subroutines and Functions
!
!     BLAS DAXPY,DSCAL,DSWAP
!     Fortran DABS,MOD
!***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
!                 *LINPACK USERS  GUIDE*, SIAM, 1979.
!***ROUTINES CALLED  DAXPY,DSCAL,DSWAP
!***END PROLOGUE  DGEDI
  INTEGER :: LDA,N,IPVT(1),JOB
  REAL*8  :: A(LDA,1),DET(2),WORK(1)
!
  REAL*8  :: T
  REAL*8  :: TEN
  INTEGER  :: I,J,K,KB,KP1,L,NM1
!
!     COMPUTE DETERMINANT
!
!***FIRST EXECUTABLE STATEMENT  DGEDI
  IF (JOB/10 .EQ. 0) GO TO 70
  DET(1) = 1.0D0
  DET(2) = 0.0D0
  TEN = 10.0D0
  DO  I = 1, N
     IF (IPVT(I) .NE. I) DET(1) = -DET(1)
     DET(1) = A(I,I)*DET(1)
!        ...EXIT
     IF (DET(1) .EQ. 0.0D0) GO TO 60
10   IF (DABS(DET(1)) .GE. 1.0D0) GO TO 20
     DET(1) = TEN*DET(1)
     DET(2) = DET(2) - 1.0D0
     GO TO 10
20   CONTINUE
30   IF (DABS(DET(1)) .LT. TEN) GO TO 40
     DET(1) = DET(1)/TEN
     DET(2) = DET(2) + 1.0D0
     GO TO 30
40   CONTINUE
  end DO
60 CONTINUE
70 CONTINUE
!
!     COMPUTE INVERSE(U)
!
  IF (MOD(JOB,10) .EQ. 0) GO TO 150
  DO  K = 1, N
     A(K,K) = 1.0D0/A(K,K)
     T = -A(K,K)
     CALL DSCAL(K-1,T,A(1,K),1)
     KP1 = K + 1
     IF (N .LT. KP1) GO TO 90
     DO  J = KP1, N
        T = A(K,J)
        A(K,J) = 0.0D0
        CALL DAXPY(K,T,A(1,K),1,A(1,J),1)
     end DO
90   CONTINUE
  end DO
!
!        FORM INVERSE(U)*INVERSE(L)
!
  NM1 = N - 1
  IF (NM1 .LT. 1) GO TO 140
  DO KB = 1, NM1
     K = N - KB
     KP1 = K + 1
     DO  I = KP1, N
        WORK(I) = A(I,K)
        A(I,K) = 0.0D0
     end DO
     DO  J = KP1, N
        T = WORK(J)
        CALL DAXPY(N,T,A(1,J),1,A(1,K),1)
     end DO
     L = IPVT(K)
     IF (L .NE. K) CALL DSWAP(N,A(1,K),1,A(1,L),1)
  end DO
140 CONTINUE
150 CONTINUE
  RETURN
END SUBROUTINE DGEDI

!SUBROUTINE DGESL(A,LDA,N,IPVT,B,JOB)
!***BEGIN PROLOGUE  DGESL
!***DATE WRITTEN   780814   (YYMMDD)
!***REVISION DATE  861211   (YYMMDD)
!***CATEGORY NO.  D2A1
!***KEYWORDS  LIBRARY=SLATEC(LINPACK),
!             TYPE=REAL*8(SGESL-S DGESL-D CGESL-C),
!             LINEAR ALGEBRA,MATRIX,SOLVE
!***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
!***PURPOSE  Solves the Real*8 system  A*X=B or  TRANS(A)*X=B
!            using the factors computed by DGECO or DGEFA.
!***DESCRIPTION
!
!     DGESL solves the Real*8 system
!     A * X = B  or  TRANS(A) * X = B
!     using the factors computed by DGECO or DGEFA.
!
!     On Entry
!
!        A       REAL*8(LDA, N)
!                the output from DGECO or DGEFA.
!
!        LDA     INTEGER
!                the leading dimension of the array  A .
!
!        N       INTEGER
!                the order of the matrix  A .
!
!        IPVT    INTEGER(N)
!                the pivot vector from DGECO or DGEFA.
!
!        B       REAL*8(N)
!                the right hand side vector.
!
!        JOB     INTEGER
!                = 0         to solve  A*X = B ,
!                = nonzero   to solve  TRANS(A)*X = B  where
!                            TRANS(A)  is the transpose.
!
!     On Return
!
!        B       the solution vector  X .
!
!     Error Condition
!
!        A division by zero will occur if the input factor contains a
!        zero on the diagonal.  Technically this indicates singularity
!        but it is often caused by improper arguments or improper
!        setting of LDA .  It will not occur if the subroutines are
!        called correctly and if DGECO has set RCOND .GT. 0.0
!        or DGEFA has set INFO .EQ. 0 .
!
!     To compute  INVERSE(A) * C  where  C  is a matrix
!     with  P  columns
!           CALL DGECO(A,LDA,N,IPVT,RCOND,Z)
!           IF (RCOND is too small) GO TO ...
!           DO 10 J = 1, P
!              CALL DGESL(A,LDA,N,IPVT,C(1,J),0)
!        10 CONTINUE
!
!     LINPACK.  This version dated 08/14/78 .
!     Cleve Moler, University of New Mexico, Argonne National Lab.
!
!     Subroutines and Functions
!
!     BLAS DAXPY,DDOT
!***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
!                 *LINPACK USERS  GUIDE*, SIAM, 1979.
!***ROUTINES CALLED  DAXPY,DDOT
!***END PROLOGUE  DGESL
!  INTEGER :: LDA,N,IPVT(*),JOB
!  REAL*8 ::  A(LDA,*),B(*)
!
!  REAL*8 :: DDOT,T
!  INTEGER :: K,KB,L,NM1
!***FIRST EXECUTABLE STATEMENT  DGESL
!  NM1 = N - 1
!  IF (JOB .NE. 0) GO TO 50
!
!  JOB = 0 , SOLVE  A * X = B
!  FIRST SOLVE  L*Y = B
!
!  IF (NM1 .LT. 1) GO TO 30
!  Kloop:  DO  K = 1, NM1
!     L = IPVT(K)
!     T = B(L)
!     IF (L .EQ. K) GO TO 10
!     B(L) = B(K)
!     B(K) = T
!10   CONTINUE
!     CALL DAXPY(N-K,T,A(K+1,K),1,B(K+1),1)
!  end DO Kloop
!30 CONTINUE
!
!        NOW SOLVE  U*X = Y
!
!  DO  KB = 1, N
!     K = N + 1 - KB
!     B(K) = B(K)/A(K,K)
!     T = -B(K)
!     CALL DAXPY(K-1,T,A(1,K),1,B(1),1)
!  end DO
!  GO TO 100
!50 CONTINUE
!
!        JOB = NONZERO, SOLVE  TRANS(A) * X = B
!        FIRST SOLVE  TRANS(U)*Y = B
!
!  DO K = 1, N
!     T = DDOT(K-1,A(1,K),1,B(1),1)
!     B(K) = (B(K) - T)/A(K,K)
!  end DO
!
!        NOW SOLVE TRANS(L)*X = Y
!
!  IF (NM1 .LT. 1) GO TO 90
!  kbloop:   DO  KB = 1, NM1
!     K = N - KB
!     B(K) = B(K) + DDOT(N-K,A(K+1,K),1,B(K+1),1)
!     L = IPVT(K)
!     IF (L .EQ. K) GO TO 70
!     T = B(L)
!     B(L) = B(K)
!     B(K) = T
!70   CONTINUE
!  end DO kbloop
!90 CONTINUE
!100 CONTINUE
!  RETURN
!END SUBROUTINE DGESL
!
!-------------------------------------------------------------------
!
!SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
!***BEGIN PROLOGUE  DAXPY
!***DATE WRITTEN   791001   (YYMMDD)
!***REVISION DATE  861211   (YYMMDD)
!***CATEGORY NO.  D1A7
!***KEYWORDS  LIBRARY=SLATEC(BLAS),
!             TYPE=REAL*8(SAXPY-S DAXPY-D CAXPY-C),
!             LINEAR ALGEBRA,TRIAD,VECTOR
!***AUTHOR  LAWSON, C. L., (JPL)
!           HANSON, R. J., (SNLA)
!           KINCAID, D. R., (U. OF TEXAS)
!           KROGH, F. T., (JPL)
!***PURPOSE  D.P computation y = a*x + y
!***DESCRIPTION
!
!                B L A S  Subprogram
!    Description of Parameters
!
!     --Input--
!        N  number of elements in input vector(s)
!       DA  real*8 scalar multiplier
!       DX  real*8 vector with N elements
!     INCX  storage spacing between elements of DX
!       DY  real*8 vector with N elements
!     INCY  storage spacing between elements of DY
!
!     --Output--
!       DY  real*8 result (unchanged if N .LE. 0)
!
!     Overwrite real*8 DY with real*8 DA*DX + DY.
!     For I = 0 to N-1, replace  DY(LY+I*INCY) with DA*DX(LX+I*INCX) +
!       DY(LY+I*INCY), where LX = 1 if INCX .GE. 0, else LX = (-INCX)*N
!       and LY is defined in a similar way using INCY.
!***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
!                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
!                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
!                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  DAXPY
!
!  REAL*8 :: DX(*),DY(*),DA
!***FIRST EXECUTABLE STATEMENT  DAXPY
!  IF(N.LE.0.OR.DA.EQ.0.D0) RETURN
!  IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60
!5 CONTINUE
!
!        CODE FOR NONEQUAL OR NONPOSITIVE INCREMENTS.
!
!  IX = 1
!  IY = 1
!  IF(INCX.LT.0)IX = (-N+1)*INCX + 1
!  IF(INCY.LT.0)IY = (-N+1)*INCY + 1
!  DO  I = 1,N
!     DY(IY) = DY(IY) + DA*DX(IX)
!     IX = IX + INCX
!     IY = IY + INCY
!  end DO
!  RETURN
!
!        CODE FOR BOTH INCREMENTS EQUAL TO 1
!
!
!        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 4.
!
!20 M = MOD(N,4)
!  IF( M .EQ. 0 ) GO TO 40
!  DO  I = 1,M
!     DY(I) = DY(I) + DA*DX(I)
!  end DO
!  IF( N .LT. 4 ) RETURN
!40 MP1 = M + 1
!  DO I = MP1,N,4
!     DY(I) = DY(I) + DA*DX(I)
!     DY(I + 1) = DY(I + 1) + DA*DX(I + 1)
!     DY(I + 2) = DY(I + 2) + DA*DX(I + 2)
!     DY(I + 3) = DY(I + 3) + DA*DX(I + 3)
!  end DO
!  RETURN
!
!        CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.
!
!60 CONTINUE
!  NS = N*INCX
!  DO I=1,NS,INCX
!     DY(I) = DA*DX(I) + DY(I)
!  end DO
!  RETURN
!END SUBROUTINE DAXPY
!
!--------------------------------------------------------------------
!

!SUBROUTINE DSCAL(N,DA,DX,INCX)
!***BEGIN PROLOGUE  DSCAL
!***DATE WRITTEN   791001   (YYMMDD)
!***REVISION DATE  861211   (YYMMDD)
!***CATEGORY NO.  D1A6
!***KEYWORDS  LIBRARY=SLATEC(BLAS),
!             TYPE=REAL*8(SSCAL-S DSCAL-D CSCAL-C),
!             LINEAR ALGEBRA,SCALE,VECTOR
!***AUTHOR  LAWSON, C. L., (JPL)
!           HANSON, R. J., (SNLA)
!           KINCAID, D. R., (U. OF TEXAS)
!           KROGH, F. T., (JPL)
!***PURPOSE  D.P. vector scale x = a*x
!***DESCRIPTION
!
!                B L A S  Subprogram
!    Description of Parameters
!
!     --Input--
!        N  number of elements in input vector(s)
!       DA  real*8 scale factor
!       DX  real*8 vector with N elements
!     INCX  storage spacing between elements of DX
!
!     --Output--
!       DX  real*8 result (unchanged if N.LE.0)
!
!     Replace real*8 DX by real*8 DA*DX.
!     For I = 0 to N-1, replace DX(1+I*INCX) with  DA * DX(1+I*INCX)
!***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
!                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
!                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
!                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  DSCAL
!
!  REAL*8 :: DA,DX(*)
!***FIRST EXECUTABLE STATEMENT  DSCAL
!  IF(N.LE.0)RETURN
!  IF(INCX.EQ.1)GOTO 20
!
!        CODE FOR INCREMENTS NOT EQUAL TO 1.
!
!  NS = N*INCX
!  DO I = 1,NS,INCX
!     DX(I) = DA*DX(I)
!  end DO
!  RETURN
!
!        CODE FOR INCREMENTS EQUAL TO 1.
!
!
!        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5.
!
!20 M = MOD(N,5)
!  IF( M .EQ. 0 ) GO TO 40
!  DO  I = 1,M
!     DX(I) = DA*DX(I)
!  end DO
!  IF( N .LT. 5 ) RETURN
!40 MP1 = M + 1
!  DO I = MP1,N,5
!     DX(I) = DA*DX(I)
!     DX(I + 1) = DA*DX(I + 1)
!     DX(I + 2) = DA*DX(I + 2)
!     DX(I + 3) = DA*DX(I + 3)
!     DX(I + 4) = DA*DX(I + 4)
!  end DO
!  RETURN
!END SUBROUTINE DSCAL
!
!----------------------------------------------------------------
!
!INTEGER FUNCTION IDAMAX(N,DX,INCX)
!***BEGIN PROLOGUE  IDAMAX
!***DATE WRITTEN   791001   (YYMMDD)
!***REVISION DATE  861211   (YYMMDD)
!***CATEGORY NO.  D1A2
!***KEYWORDS  LIBRARY=SLATEC(BLAS),
!             TYPE=REAL*8(ISAMAX-S IDAMAX-D ICAMAX-C),
!             LINEAR ALGEBRA,MAXIMUM COMPONENT,VECTOR
!***AUTHOR  LAWSON, C. L., (JPL)
!           HANSON, R. J., (SNLA)
!           KINCAID, D. R., (U. OF TEXAS)
!           KROGH, F. T., (JPL)
!***PURPOSE  Find the smallest index of that component of a d.p. vector
!            having the maximum magnitude.
!***DESCRIPTION
!
!                B L A S  Subprogram
!    Description of Parameters
!
!     --Input--
!        N  number of elements in input vector(s)
!       DX  real*8 vector with N elements
!     INCX  storage spacing between elements of DX
!
!     --Output--
!   IDAMAX  smallest index (zero if N .LE. 0)
!
!     Find smallest index of maximum magnitude of real*8 DX.
!     IDAMAX =  first I, I = 1 to N, to minimize  ABS(DX(1-INCX+I*INCX)
!***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
!                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
!                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
!                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  IDAMAX
!
!  REAL*8 :: DX(*),DMAX,XMAG
!***FIRST EXECUTABLE STATEMENT  IDAMAX
!  IDAMAX = 0
!  IF(N.LE.0) RETURN
!  IDAMAX = 1
!  IF(N.LE.1)RETURN
!  IF(INCX.EQ.1)GOTO 20
!
!        CODE FOR INCREMENTS NOT EQUAL TO 1.
!
!  DMAX = DABS(DX(1))
!  NS = N*INCX
!  II = 1
!  DO I = 1,NS,INCX
!     XMAG = DABS(DX(I))
!     IF(XMAG.LE.DMAX) GO TO 5
!     IDAMAX = II
!     DMAX = XMAG
!5    II = II + 1
!  end DO
!  RETURN
!
!        CODE FOR INCREMENTS EQUAL TO 1.
!
!20 DMAX = DABS(DX(1))
!  DO I = 2,N
!     XMAG = DABS(DX(I))
!     IF(XMAG.LE.DMAX) CYCLE
!     IDAMAX = I
!     DMAX = XMAG
!  end DO
!  RETURN
!END FUNCTION IDAMAX
!
!-----------------------------------------------------------------
!
SUBROUTINE DSWAP(N,DX,INCX,DY,INCY)
!***BEGIN PROLOGUE  DSWAP
!***DATE WRITTEN   791001   (YYMMDD)
!***REVISION DATE  861211   (YYMMDD)
!***CATEGORY NO.  D1A5
!***KEYWORDS  LIBRARY=SLATEC(BLAS),
!             TYPE=REAL*8(SSWAP-S DSWAP-D CSWAP-C ISWAP-I),
!             INTERCHANGE,LINEAR ALGEBRA,VECTOR
!***AUTHOR  LAWSON, C. L., (JPL)
!           HANSON, R. J., (SNLA)
!           KINCAID, D. R., (U. OF TEXAS)
!           KROGH, F. T., (JPL)
!***PURPOSE  Interchange d.p. vectors
!***DESCRIPTION
!
!                B L A S  Subprogram
!    Description of Parameters
!
!     --Input--
!        N  number of elements in input vector(s)
!       DX  real*8 vector with N elements
!     INCX  storage spacing between elements of DX
!       DY  real*8 vector with N elements
!     INCY  storage spacing between elements of DY
!
!     --Output--
!       DX  input vector DY (unchanged if N .LE. 0)
!       DY  input vector DX (unchanged if N .LE. 0)
!
!     Interchange real*8 DX and real*8 DY.
!     For I = 0 to N-1, interchange  DX(LX+I*INCX) and DY(LY+I*INCY),
!     where LX = 1 if INCX .GE. 0, else LX = (-INCX)*N, and LY is
!     defined in a similar way using INCY.
!***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
!                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
!                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
!                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  DSWAP
!
  REAL*8 :: DX(*),DY(*),DTEMP1,DTEMP2,DTEMP3
!***FIRST EXECUTABLE STATEMENT  DSWAP
  IF(N.LE.0)RETURN
  IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60
5 CONTINUE
!
!       CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.
!
  IX = 1
  IY = 1
  IF(INCX.LT.0)IX = (-N+1)*INCX + 1
  IF(INCY.LT.0)IY = (-N+1)*INCY + 1
  DO  I = 1,N
     DTEMP1 = DX(IX)
     DX(IX) = DY(IY)
     DY(IY) = DTEMP1
     IX = IX + INCX
     IY = IY + INCY
  end DO
  RETURN
!
!       CODE FOR BOTH INCREMENTS EQUAL TO 1
!
!
!       CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 3.
!
20 M = MOD(N,3)
  IF( M .EQ. 0 ) GO TO 40
  DO I = 1,M
     DTEMP1 = DX(I)
     DX(I) = DY(I)
     DY(I) = DTEMP1
  end DO
  IF( N .LT. 3 ) RETURN
40 MP1 = M + 1
  DO  I = MP1,N,3
     DTEMP1 = DX(I)
     DTEMP2 = DX(I+1)
     DTEMP3 = DX(I+2)
     DX(I) = DY(I)
     DX(I+1) = DY(I+1)
     DX(I+2) = DY(I+2)
     DY(I) = DTEMP1
     DY(I+1) = DTEMP2
     DY(I+2) = DTEMP3
  end DO
  RETURN
60 CONTINUE
!
!     CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.
!
!      NS = N*INCX
  DO  I=1,NS,INCX
     DTEMP1 = DX(I)
     DX(I) = DY(I)
     DY(I) = DTEMP1
  end DO
  RETURN
END SUBROUTINE DSWAP

!
!----------------------------------------------------------------------C
!
SUBROUTINE CKCVCoeff (T, ICKWRK, RCKWRK, Cvs, IPolyOrder)
!
!  START PROLOGUE
!
!     Returns the specific heats at constant pressure in mass units;
!     see Eq. (26).
!
!  INPUT
!     T      - Temperature.
!                   cgs units - K
!                   Data type - real scalar
!     ICKWRK - Array of integer workspace.
!                   Data type - integer array
!                   Dimension ICKWRK(*) at least LENIWK.
!     RCKWRK - Array of real work space.
!                   Data type - real array
!                   Dimension RCKWRK(*) at least LENRWK.
!
!  OUTPUT
!     Cvs   - Specific heats at constant pressure in mass units
!              for the species.
!                   cgs units - ergs/(gm*K)
!                   Data type - real array
!                   Dimension CPMS(*) at least KK, the total number of
!                   species.
!    IPolyOrder - order of polynomila interpolation
!
!  END PROLOGUE
!
  IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
!
  INCLUDE 'ckstrt.fh'
!
  DIMENSION ICKWRK(*), RCKWRK(*), Cvs(NKK,*)
!
!
  IPolyOrder = NCP
  if(T < 0d0) RETURN

  species_loop : DO  K = 1, NKK
     L = 1
     range_loop : DO  N = 2, ICKWRK(IcNT + K - 1)-1
        TEMP = RCKWRK(NcTT + (K-1)*MXTP + N - 1)
        IF (T .GT. TEMP) L = L+1
     end DO range_loop
!
     termK = RCKWRK(NcRU)/RCKWRK(NcWT + K - 1)
     NA1 = NcAA + (L-1)*NCP2 + (K-1)*NCP2T
!
     poly_loop :DO N = 1, NCP
        CVS(K,N) = RCKWRK(NA1 + N - 1)*termK
     end DO poly_loop
     CVS(K,1) = CVS(K,1) - termK
     CVS(K,NCP+1) = RCKWRK(NA1 + NCP1 - 1)*termK
!
  end DO species_loop
  RETURN
END SUBROUTINE CKCVCoeff
!
!----------------------------------------------------------------------C
!
SUBROUTINE CKCPCoeff (T, ICKWRK, RCKWRK, Cps, IPolyOrder)
!
!  START PROLOGUE
!
!     Returns the specific heats at constant pressure in mass units;
!     see Eq. (26).
!
!  INPUT
!     T      - Temperature.
!                   cgs units - K
!                   Data type - real scalar
!     ICKWRK - Array of integer workspace.
!                   Data type - integer array
!                   Dimension ICKWRK(*) at least LENIWK.
!     RCKWRK - Array of real work space.
!                   Data type - real array
!                   Dimension RCKWRK(*) at least LENRWK.
!
!  OUTPUT
!     Cps   - Specific heats at constant pressure in mass units
!              for the species.
!                   cgs units - ergs/(gm*K)
!                   Data type - real array
!                   Dimension CPMS(*) at least KK, the total number of
!                   species.
!
!  END PROLOGUE
!
  IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
!
  INCLUDE 'ckstrt.fh'
!
  DIMENSION ICKWRK(*), RCKWRK(*), Cps(NKK,*)
!
!
  IPolyOrder = NCP
  if(T < 0d0) RETURN

  species_loop : DO  K = 1, NKK
     L = 1
     range_loop : DO  N = 2, ICKWRK(IcNT + K - 1)-1
        TEMP = RCKWRK(NcTT + (K-1)*MXTP + N - 1)
        IF (T .GT. TEMP) L = L+1
     end DO range_loop
!
     termK = RCKWRK(NcRU)/RCKWRK(NcWT + K - 1)
     NA1 = NcAA + (L-1)*NCP2 + (K-1)*NCP2T
!
     poly_loop :DO N = 1, NCP
        CPS(K,N) = RCKWRK(NA1 + N - 1)*termK
     end DO poly_loop
     CPS(K,NCP+1) = RCKWRK(NA1 + NCP1 - 1)*termK
!
  end DO species_loop
  RETURN
END SUBROUTINE CKCPCoeff
!
!----------------------------------------------------------------------C
!
SUBROUTINE CKCVCoeffML (T, ICKWRK, RCKWRK, Cvs, CVSInt, IPolyOrder)
!
!  START PROLOGUE
!
!     Returns the specific heats at constant volume in molar units;
!     see Eq. (26).
!
!  INPUT
!     T      - Temperature.
!                   cgs units - K
!                   Data type - real scalar
!     ICKWRK - Array of integer workspace.
!                   Data type - integer array
!                   Dimension ICKWRK(*) at least LENIWK.
!     RCKWRK - Array of real work space.
!                   Data type - real array
!                   Dimension RCKWRK(*) at least LENRWK.
!
!  OUTPUT
!     Cvs    - Specific heats at constant pressure in mass units
!              for the species.
!                   cgs units - ergs/(mole*K)
!                   Data type - real array
!                   Dimension CPMS(*) at least KK, the total number of
!                   species.
!
!  END PROLOGUE
!
  IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
!
  INCLUDE 'ckstrt.fh'
!
  DIMENSION ICKWRK(*), RCKWRK(*), Cvs(NKK,*), CvsInt(NKK,*)
!
!
  IPolyOrder = NCP
  if(T < 0d0) RETURN

  species_loop : DO  K = 1, NKK
     L = 1
     range_loop : DO  N = 2, ICKWRK(IcNT + K - 1)-1
        TEMP = RCKWRK(NcTT + (K-1)*MXTP + N - 1)
        IF (T .GT. TEMP) L = L+1
     end DO range_loop
!
     termK = RCKWRK(NcRU)
     NA1 = NcAA + (L-1)*NCP2 + (K-1)*NCP2T
!
     poly_loop :DO N = 1, NCP
        CVS(K,N) = RCKWRK(NA1 + N - 1)*termK
     end DO poly_loop
     CVS(K,1) = CVS(K,1) - termK
     CVS(K,NCP+1) = RCKWRK(NA1 + NCP1 - 1)*termK

     poly_loopII :DO N = 1, NCP
        CVSInt(K,N) =  CVS(K,N) / dble(N)
     end DO poly_loopII
!
  end DO species_loop
  RETURN
END SUBROUTINE CKCVCoeffML
!
!----------------------------------------------------------------------C
!
SUBROUTINE CKCVCoeffAvg(CY,Cvs,CvsAvg,CvsAvgInt)
!
!  START PROLOGUE
!
!  END PROLOGUE
!
  IMPLICIT NONE
!
  INCLUDE 'ckstrt.fh'
!
  REAL*8 :: CvsAvg(*),CvsAvgInt(*),Cvs(NKK,*),CY(*)
  INTEGER :: N,K
!
!
  poly_loop :DO N = 1, NCP+1

     CVSAVG(N) = Cvs(1,N)*CY(1)
     species_loop : DO  K = 2, NKK
!
        CVSAvg(N) = CVSAvg(N) + Cvs(K,N)*CY(K)
!
     END DO species_loop
!
  end DO poly_loop

  poly_loopII :DO N = 1, NCP
     CVSAvgInt(N) =  CVSAvg(N) / dble(N)
  end DO poly_loopII

  RETURN
END SUBROUTINE CKCVCoeffAvg
!
!-------------------------------------------------
!
SUBROUTINE CKMMWY (Y, ICKWRK, RCKWRK, WTM)
!
!  START PROLOGUE
!
!  SUBROUTINE CKMMWY (Y, ICKWRK, RCKWRK, WTM)
!     Returns the mean molecular weight of the gas mixture given the
!     mass fractions;  see Eq. (3).
!
!  INPUT
!     Y      - Mass fractions of the species.
!                   cgs units - none
!                   Data type - real array
!                   Dimension Y(*) at least KK, the total number of
!                   species.
!     ICKWRK - Array of integer workspace.
!                   Data type - integer array
!                   Dimension ICKWRK(*) at least LENIWK.
!     RCKWRK - Array of real work space.
!                   Data type - real array
!                   Dimension RCKWRK(*) at least LENRWK.
!
!  OUTPUT
!     WTM    - Mean molecular weight of the species mixture.
!                   cgs units - gm/mole
!                   Data type - real scalar
!
!  END PROLOGUE
!
!*****precision > double
  IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
!*****END precision > double
!*****precision > single
!        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
!*****END precision > single
!
  DIMENSION Y(*), ICKWRK(*), RCKWRK(*)
!
  INCLUDE 'ckstrt.fh'
!
  SUMYOW=0.0
  DO  K = 1, NKK
     SUMYOW = SUMYOW + Y(K)/RCKWRK(NcWT + K - 1)
  END DO
  WTM = 1.0/SUMYOW
  RETURN
END SUBROUTINE CKMMWY
!
!-----------------------------------------------------------------------
!
real*8 function solve_quintic_newton(U,a_298,b_298,Tref,success) Result(Tguess)
  implicit None
  INCLUDE 'ckstrt.fh'
  real*8,INTENT(IN) :: U,a_298(*),b_298(*),Tref
  Logical,INTENT(OUT) :: success
  real*8 :: funcT,deriv,Tinitguess
  integer :: k,iter,Maxiter

!Solve the quadratic as initial guess
!!>  Tguess = 0.5d0/b_298(2)*(-b_298(1)+sqrt(b_298(1)**2+4d0*b_298(2)*U))
  Tguess = Tref
  Tinitguess = Tguess

  funcT = b_298(1)*Tguess - U 
  do k = 2,NCP
     funcT = funcT + b_298(k)*Tguess**k
  enddo
  Maxiter = 100
  iter = 0
  do while (abs(funcT) > 1d-7*U .and. iter < Maxiter )
     deriv = a_298(1) 
     do k = 2,NCP
        deriv = deriv + a_298(k)*Tguess**(k-1)
     enddo
!
     Tguess = Tguess - funcT/deriv
!
     funcT = b_298(1)*Tguess - U 
     do k = 2,NCP
        funcT = funcT + b_298(k)*Tguess**k
     enddo
     iter = iter+1
  enddo

  success = iter<Maxiter-1
  if(.not.success) then
     write(*,'(i3,a,1p123e12.4)') iter,'CVP solve_quintic_newton, Tout,fT,U,Tguess',&
          &Tguess,funcT,U,Tinitguess
  end if

end function solve_quintic_newton
!
!----------------------------------------------------------------------C
!
SUBROUTINE CKHMS  (T, ICKWRK, RCKWRK, HMS)
!
!  START PROLOGUE
!
!  SUBROUTINE CKHMS  (T, ICKWRK, RCKWRK, HMS)
!     Returns the enthalpies in mass units;  see Eq. (27).
!
!  INPUT
!     T      - Temperature.
!                   cgs units - K
!                   Data type - real scalar
!     ICKWRK - Array of integer workspace.
!                   Data type - integer array
!                   Dimension ICKWRK(*) at least LENIWK.
!     RCKWRK - Array of real work space.
!                   Data type - real array
!                   Dimension RCKWRK(*) at least LENRWK.
!  OUTPUT
!     HMS    - Array of Enthalpies in mass units for the species.
!                   cgs units - ergs/gm
!                   Data type - real array
!                   Dimension HMS(*) at least KK, the total number of
!                   species.
!
!  END PROLOGUE
!
  IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
!
  INCLUDE 'ckstrt.fh'
!
  DIMENSION ICKWRK(*), RCKWRK(*), HMS(*), TN(10)
!
  RUT = T*RCKWRK(NcRU)
  TN(1)=1.0
  DO  N = 2, NCP
     TN(N) = T**(N-1)/N
  END DO
!
  DO  K = 1, NKK
     L = 1
     DO  N = 2, ICKWRK(IcNT + K - 1)-1
        TEMP = RCKWRK(NcTT + (K-1)*MXTP + N - 1)
        IF (T .GT. TEMP) L = L+1
     END DO
!
     NA1 = NcAA + (L-1)*NCP2 + (K-1)*NCP2T
     SUM = 0.0
     DO  N = 1, NCP
        SUM = SUM + TN(N)*RCKWRK(NA1 + N - 1)
     END DO
     HMS(K) = RUT * (SUM + RCKWRK(NA1 + NCP1 - 1)/T)&
          / RCKWRK(NcWT + K - 1)
  END DO
  RETURN
END SUBROUTINE CKHMS
!
!----------------------------------------------------------------------C
!
SUBROUTINE CKHBMS  (T, Y, ICKWRK, RCKWRK, HBMS)
!
!  START PROLOGUE
!
!  SUBROUTINE CKHBMS  (T, Y, ICKWRK, RCKWRK, HBMS)
!     Returns the enthalpies in mass units;  see Eq. (27).
!
!  INPUT
!     T      - Temperature.
!                   cgs units - K
!                   Data type - real scalar
!     ICKWRK - Array of integer workspace.
!                   Data type - integer array
!                   Dimension ICKWRK(*) at least LENIWK.
!     RCKWRK - Array of real work space.
!                   Data type - real array
!                   Dimension RCKWRK(*) at least LENRWK.
!  OUTPUT
!     HMS    - Array of Enthalpies in mass units for the species.
!                   cgs units - ergs/gm
!                   Data type - real array
!                   Dimension HMS(*) at least KK, the total number of
!                   species.
!
!  END PROLOGUE
!
  IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
!
  INCLUDE 'ckstrt.fh'
!
  DIMENSION ICKWRK(*), RCKWRK(*), Y(*), TN(10)
!
  RUT = T*RCKWRK(NcRU)
  TN(1)=1.0
  DO  N = 2, NCP
     TN(N) = T**(N-1)/N
  END DO
!
 HBMS = 0.0
  DO  K = 1, NKK
     L = 1
     DO  N = 2, ICKWRK(IcNT + K - 1)-1
        TEMP = RCKWRK(NcTT + (K-1)*MXTP + N - 1)
        IF (T .GT. TEMP) L = L+1
     END DO
!
     NA1 = NcAA + (L-1)*NCP2 + (K-1)*NCP2T
     SUM = 0.0
     DO  N = 1, NCP
        SUM = SUM + TN(N)*RCKWRK(NA1 + N - 1)
     END DO
     HBMS = HBMS+ (Y(K))*(RUT * (SUM + RCKWRK(NA1 + NCP1 - 1)/T)&
          / RCKWRK(NcWT + K - 1))
  END DO
  RETURN
END SUBROUTINE CKHBMS

!
!----------------------------------------------------------------------C
!
SUBROUTINE CKWmsYP  (P, T, Y, ICKWRK, RCKWRK, WDOT)
!
!  START PROLOGUE
!
!  SUBROUTINE CKWmsYP  (P, T, Y, ICKWRK, RCKWRK, WDOT)
!     Returns the mass production rates of the species given the
!     pressure, temperature and mass fractions;  see Eq. (49).
!
!  INPUT
!     P      - Pressure.
!                   cgs units - dynes/cm**2
!                   Data type - real scalar
!     T      - Temperature.
!                   cgs units - K
!                   Data type - real scalar
!     Y      - Mass fractions of the species.
!                   cgs units - none
!                   Data type - real array
!                   Dimension Y(*) at least KK, the total number of
!                   species.
!     ICKWRK - Array of integer workspace.
!                   Data type - integer array
!                   Dimension ICKWRK(*) at least LENIWK.
!     RCKWRK - Array of real work space.
!                   Data type - real array
!                   Dimension RCKWRK(*) at least LENRWK.
!
!  OUTPUT
!     WDOT   - Chemical mass production rates of the species.
!                   cgs units - gm/(cm**3*sec)
!                   Data type - real array
!                   Dimension WDOT(*) at least KK, the total number of
!                   species.
!                  
!                   Routines called CKRATT CKYTCP CKRATX
!
!  END PROLOGUE
!
  IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
!incf90
  INCLUDE 'ckstrt.fh'
!
  DIMENSION ICKWRK(*), RCKWRK(*), Y(*), WDOT(*)
!
  CALL CKRATT (RCKWRK, ICKWRK, NII, MXSP, RCKWRK(NcRU),&
       RCKWRK(NcPA), T, ICKWRK(IcNS), ICKWRK(IcNU),&
       ICKWRK(IcNK), NPAR+1, RCKWRK(NcCO), NREV,&
       ICKWRK(IcRV), RCKWRK(NcRV), NLAN, NLAR, ICKWRK(IcLT),&
       RCKWRK(NcLT), NRLT, ICKWRK(IcRL), RCKWRK(NcRL),&
       RCKWRK(NcK1), RCKWRK(NcKF), RCKWRK(NcKR),&
       RCKWRK(NcI1), NRNU, ICKWRK(IcRNU), RCKWRK(NcRNU))
!
  CALL CKYTCP (P, T, Y, ICKWRK, RCKWRK, RCKWRK(NcK1))
!
  CALL CKRATX (NII, NKK, MXSP, MXTB, T, RCKWRK(NcK1), ICKWRK(IcNS),&
       ICKWRK(IcNU), ICKWRK(IcNK), NPAR+1, RCKWRK(NcCO),&
       NFAL, ICKWRK(IcFL), ICKWRK(IcFO), ICKWRK(IcKF), NFAR, &
       RCKWRK(NcFL), NTHB, ICKWRK(IcTB), ICKWRK(IcKN), &
       RCKWRK(NcKT), ICKWRK(IcKT), RCKWRK(NcKF), &
       RCKWRK(NcKR), RCKWRK(NcI1), RCKWRK(NcI2), &
       RCKWRK(NcI3), NRNU, ICKWRK(IcRNU), RCKWRK(NcRNU),&
       NORD, ICKWRK(IcORD), MXORD, ICKWRK(IcKOR), &
       RCKWRK(NcKOR))
!
  DO  K = 1, NKK
     WDOT(K) = 0D0
  END DO
  N_loop: DO N = 1, MXSP
     I_loop: DO I = 1, NII
        K = ICKWRK(IcNK + (I-1)*MXSP + N - 1)
        ROP = RCKWRK(NcI1+I-1) - RCKWRK(NcI2+I-1)
        IF (K .NE. 0) THEN
           RNU = DBLE(ICKWRK(IcNU + (I-1)*MXSP + N - 1))
           WDOT(K) = WDOT(K) + RNU * ROP
        ENDIF
     END DO I_loop
  END DO N_loop
!
  IF (NRNU .LE. 0) GOTO 999
!
  L_loop: DO L = 1, NRNU
     I = ICKWRK(IcRNU + L - 1)
     N2_loop: DO  N = 1, MXSP
        K = ICKWRK(IcNK + (I-1)*MXSP + N - 1)
        ROP = RCKWRK(NcI1+I-1) - RCKWRK(NcI2+I-1)
        IF (K .NE. 0) THEN
           RNU = RCKWRK(NcRNU + (L-1)*MXSP + N - 1)
           WDOT(K) = WDOT(K) + RNU * ROP
        ENDIF
     END DO N2_loop
  END DO L_loop
!
999 CONTINUE
!
  molweight_loop: DO  K = 1, NKK
     WDOT(K) = WDOT(K)*RCKWRK(NcWT + K - 1)
  END DO  molweight_loop
!
  RETURN
END SUBROUTINE CKWmsYP
!
!----------------------------------------------------------------------C
!
!
!----------------------------------------------------------------------C
!
SUBROUTINE CKWmsYRHO  (RHO, T, Y, ICKWRK, RCKWRK, WDOT)
!
!  START PROLOGUE
!
!  SUBROUTINE CKWmsYP  (P, T, Y, ICKWRK, RCKWRK, WDOT)
!     Returns the mass production rates of the species given the
!     pressure, temperature and mass fractions;  see Eq. (49).
!
!  INPUT
!     RHO      - density
!                   cgs units - gm/cm^3
!                   Data type - real scalar
!     T      - Temperature.
!                   cgs units - K
!                   Data type - real scalar
!     Y      - Mass fractions of the species.
!                   cgs units - none
!                   Data type - real array
!                   Dimension Y(*) at least KK, the total number of
!                   species.
!     ICKWRK - Array of integer workspace.
!                   Data type - integer array
!                   Dimension ICKWRK(*) at least LENIWK.
!     RCKWRK - Array of real work space.
!                   Data type - real array
!                   Dimension RCKWRK(*) at least LENRWK.
!
!  OUTPUT
!     WDOT   - Chemical mass production rates of the species.
!                   cgs units - gm/(cm**3*sec)
!                   Data type - real array
!                   Dimension WDOT(*) at least KK, the total number of
!                   species.
!                  
!                   Routines called CKRATT CKYTCP CKRATX
!
!  END PROLOGUE
!
  IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
!incf90
  INCLUDE 'ckstrt.fh'
!
  DIMENSION ICKWRK(*), RCKWRK(*), Y(*), WDOT(*)
!
  CALL CKRATT (RCKWRK, ICKWRK, NII, MXSP, RCKWRK(NcRU),         &
       RCKWRK(NcPA), T, ICKWRK(IcNS), ICKWRK(IcNU),             &
       ICKWRK(IcNK), NPAR+1, RCKWRK(NcCO), NREV,                &
       ICKWRK(IcRV), RCKWRK(NcRV), NLAN, NLAR, ICKWRK(IcLT),    &
       RCKWRK(NcLT), NRLT, ICKWRK(IcRL), RCKWRK(NcRL),          &
       RCKWRK(NcK1), RCKWRK(NcKF), RCKWRK(NcKR),                &
       RCKWRK(NcI1), NRNU, ICKWRK(IcRNU), RCKWRK(NcRNU))
!
  DO K = 1,NKK
     RCKWRK(NcK1+K-1) = Y(K)*RHO/RCKWRK(NcWT + K - 1)
  END DO
!
  CALL CKRATX (NII, NKK, MXSP, MXTB, T, RCKWRK(NcK1), ICKWRK(IcNS), &
       ICKWRK(IcNU), ICKWRK(IcNK), NPAR+1, RCKWRK(NcCO),            &
       NFAL, ICKWRK(IcFL), ICKWRK(IcFO), ICKWRK(IcKF), NFAR,        &
       RCKWRK(NcFL), NTHB, ICKWRK(IcTB), ICKWRK(IcKN),              &
       RCKWRK(NcKT), ICKWRK(IcKT), RCKWRK(NcKF),                    &
       RCKWRK(NcKR), RCKWRK(NcI1), RCKWRK(NcI2),                    &
       RCKWRK(NcI3), NRNU, ICKWRK(IcRNU), RCKWRK(NcRNU),            &
       NORD, ICKWRK(IcORD), MXORD, ICKWRK(IcKOR),                   &
       RCKWRK(NcKOR))
!
  DO  K = 1, NKK
     WDOT(K) = 0D0
  END DO
  N_loop: DO N = 1, MXSP
     I_loop: DO I = 1, NII
        K = ICKWRK(IcNK + (I-1)*MXSP + N - 1)
        ROP = RCKWRK(NcI1+I-1) - RCKWRK(NcI2+I-1)
        IF (K .NE. 0) THEN
           RNU = DBLE(ICKWRK(IcNU + (I-1)*MXSP + N - 1))
           WDOT(K) = WDOT(K) + RNU * ROP
        ENDIF
     END DO I_loop
  END DO N_loop
!
  IF (NRNU .LE. 0) GOTO 999
!
  L_loop: DO L = 1, NRNU
     I = ICKWRK(IcRNU + L - 1)
     N2_loop: DO  N = 1, MXSP
        K = ICKWRK(IcNK + (I-1)*MXSP + N - 1)
        ROP = RCKWRK(NcI1+I-1) - RCKWRK(NcI2+I-1)
        IF (K .NE. 0) THEN
           RNU = RCKWRK(NcRNU + (L-1)*MXSP + N - 1)
           WDOT(K) = WDOT(K) + RNU * ROP
        ENDIF
     END DO N2_loop
  END DO L_loop
!
999 CONTINUE
!
  molweight_loop: DO  K = 1, NKK
     WDOT(K) = WDOT(K)*RCKWRK(NcWT + K - 1)
  END DO  molweight_loop
!
  RETURN
END SUBROUTINE CKWmsYRHO
!
!----------------------------------------------------------------------C
!
SUBROUTINE CKRATT (RCKWRK, ICKWRK, II, MAXSP, RU, PATM, T, NSPEC,&
     NU, NUNK, NPAR, PAR, NREV, IREV, RPAR, NLAN,&
     NLAR, ILAN, PLT, NRLT, IRLT, RPLT, SMH, RKFT,&
     RKRT, EQK, NRNU, IRNU, RNU)
!
!  START PROLOGUE
!
!  SUBROUTINE CKRATT (RCKWRK, ICKWRK, II, MAXSP, RU, PATM, T, NSPEC,
! 1                   NU, NUNK, NPAR, PAR, NREV, IREV, RPAR, NLAN,
! 2                   NLAR, ILAN, PLT, NRLT, IRLT, RPLT, SMH, RKFT,
! 3                   RKRT, EQK, NRNU, IRNU, RNU)
!
!  END PROLOGUE
!
  IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
!
  DIMENSION RCKWRK(*), ICKWRK(*), NSPEC(*), NU(MAXSP,*),&
       NUNK(MAXSP,*), PAR(NPAR,*), IREV(*), RPAR(NPAR,*),&
       ILAN(*), IRLT(*), PLT(NLAR,*), RPLT(NLAR,*), SMH(*),&
       RKFT(*), RKRT(*), EQK(*), IRNU(*), RNU(MAXSP,*)
!
  COMMON /MACH/ SMALL,BIG,EXPARG
!
  ALOGT = LOG(T)
!
  DO  I = 1, II
     RKFT(I) = PAR(1,I) * EXP(PAR(2,I)*ALOGT - PAR(3,I)/T)
  END DO
!
!Landau-Teller reactions
!
  DO  N = 1, NLAN
     I = ILAN(N)
     TFAC = PLT(1,N)/T**(1.0/3.0) + PLT(2,N)/T**(2.0/3.0)
     RKFT(I) = RKFT(I) * EXP(TFAC)
  END DO
!
!Reverse reaction 
!
  CALL CKSMH (T, ICKWRK, RCKWRK, SMH)
  DO  I = 1, II
     SUMSMH = 0.0
     DO N = 1, MAXSP
        IF (NUNK(N,I).NE.0) SUMSMH = SUMSMH + DBLE(NU(N,I))*SMH(NUNK(N,I))
     END DO
     IF (SUMSMH .NE. 0.0) EQK(I) = EXP(MIN(SUMSMH,EXPARG))
  END DO
!
  DO  N = 1, NRNU
     SUMSMH = 0.0
     I = IRNU(N)
     DO  L = 1, MAXSP
        IF (NUNK(L,I).NE.0) SUMSMH=SUMSMH+RNU(L,N)*SMH(NUNK(L,I))
     END DO
     IF (SUMSMH .NE. 0.0) EQK(I) = EXP(MIN(SUMSMH,EXPARG))
  END DO
!
  PFAC = PATM / (RU*T)
  DO  I = 1, II
     NUSUMK = NU(1,I)+NU(2,I)+NU(3,I)+NU(4,I)+NU(5,I)+NU(6,I)
     EQK(I) = EQK(I) * PFAC**NUSUMK
  END DO
  DO  N = 1, NRNU
     RNUSUM = RNU(1,N)+RNU(2,N)+RNU(3,N)+RNU(4,N)+RNU(5,N)+RNU(6,N)
     I = IRNU(N)
     PFR = PFAC ** RNUSUM
     EQK(I) = EQK(I) * PFR
  END DO
!
  DO  I = 1, II
!
!     RKR=0.0 for irreversible reactions, else RKR=RKF/MAX(EQK,SMALL)
!
     RKRT(I) = 0.0
     IF (NSPEC(I).GT.0) RKRT(I) = RKFT(I) / MAX(EQK(I),SMALL)
  END DO
!
!     if reverse parameters have been given:
!
  DO  N = 1, NREV
     I = IREV(N)
     RKRT(I) = RPAR(1,N) * EXP(RPAR(2,N)*ALOGT - RPAR(3,N)/T)
     EQK(I)  = RKFT(I)/RKRT(I)
  END DO
!
!     if reverse Landau-Teller parameters have been given:
!
  DO  N = 1, NRLT
     I = IRLT(N)
     TFAC = RPLT(1,N)/T**(1.0/3.0) + RPLT(2,N)/T**(2.0/3.0)
     RKRT(I) = RKRT(I) * EXP(TFAC)
     EQK(I) = RKFT(I)/RKRT(I)
  END DO
!
  RETURN
END SUBROUTINE CKRATT
!
!----------------------------------------------------------------------C
!
SUBROUTINE CKSMH  (T, ICKWRK, RCKWRK, SMH)
!
!  START PROLOGUE
!
!  SUBROUTINE CKSMH  (T, ICKWRK, RCKWRK, SMH)*
!     Returns the array of entropies minus enthalpies for the species.
!     It is normally not called directly by the user.
!
!  INPUT
!     T      - Temperature.
!                   cgs units - K
!                   Data type - real scalar
!     ICKWRK - Array of integer workspace.
!                   Data type - integer array
!                   Dimension ICKWRK(*) at least LENIWK.
!     RCKWRK - Array of real work space.
!                   Data type - real array
!                   Dimension RCKWRK(*) at least LENRWK.
!
!  OUTPUT
!     SMH    - Entropy minus enthalpy for the species,
!              SMH(K) = S(K)/R - H(K)/RT.
!                   cgs units - none
!                   Data type - real array
!                   Dimension SMH(*) at least KK, the total number of
!                   species.
!
!  END PROLOGUE
!
  IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
!
  INCLUDE 'ckstrt.fh'
!
  DIMENSION ICKWRK(*), RCKWRK(*), SMH(*), TN(10)
!
  TN(1) = LOG(T) - 1.0
  DO  N = 2, NCP
     TN(N) = T**(N-1)/((N-1)*N)
  END DO
!
  species_loop : DO  K = 1, NKK
     L = 1
     range_loop : DO  N = 2, ICKWRK(IcNT + K - 1)-1
        TEMP = RCKWRK(NcTT + (K-1)*MXTP + N - 1)
        IF (T .GT. TEMP) L = L+1
     end DO range_loop
!
     NA1 = NcAA + (L-1)*NCP2 + (K-1)*NCP2T
     SUM = 0.0
     poly_loop : DO  N = 1, NCP
        SUM = SUM + TN(N)*RCKWRK(NA1 + N - 1)
     end DO poly_loop
     SMH(K) = SUM + RCKWRK(NA1 + NCP2 - 1) - RCKWRK(NA1 + NCP1 - 1)/T
!
  END DO species_loop
  RETURN
END SUBROUTINE CKSMH
!
!----------------------------------------------------------------------C
!
SUBROUTINE JAC_CKSMH  (T, ICKWRK, RCKWRK, SMH, D_SMH)
!
!  START PROLOGUE
!
!  SUBROUTINE JAC_CKSMH  (T, ICKWRK, RCKWRK, SMH)*
!     Returns the array of entropies minus enthalpies for the species.
!     It is normally not called directly by the user.
!
!  INPUT
!     T      - Temperature.
!                   cgs units - K
!                   Data type - real scalar
!     ICKWRK - Array of integer workspace.
!                   Data type - integer array
!                   Dimension ICKWRK(*) at least LENIWK.
!     RCKWRK - Array of real work space.
!                   Data type - real array
!                   Dimension RCKWRK(*) at least LENRWK.
!
!  OUTPUT
!     SMH    - Entropy minus enthalpy for the species,
!              SMH(K) = S(K)/R - H(K)/RT.
!                   cgs units - none
!                   Data type - real array
!                   Dimension SMH(*) at least KK, the total number of
!                   species.
!
!  END PROLOGUE
!
  IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
!
  INCLUDE 'ckstrt.fh'
!
  DIMENSION ICKWRK(*), RCKWRK(*), SMH(*), D_SMH(*), TN(10),D_TN(10)
!
  TN(1) = LOG(T) - 1.0
  D_TN(1) = 1d0/T 
  DO  N = 2, NCP
     D_TN(N) = T**(N-2)/DBLE(N)
     TN(N)   = D_TN(N)*T/DBLE(N-1)
  END DO
  RINVTsq = 1d0/(T*T)
!
  species_loop : DO  K = 1, NKK
     L = 1
     range_loop : DO  N = 2, ICKWRK(IcNT + K - 1)-1
        TEMP = RCKWRK(NcTT + (K-1)*MXTP + N - 1)
        IF (T .GT. TEMP) L = L+1
     end DO range_loop
!
     NA1 = NcAA + (L-1)*NCP2 + (K-1)*NCP2T
     SUM = 0D0
     D_SUM = 0d0
     poly_loop : DO  N = 1, NCP
        SUM = SUM + TN(N)*RCKWRK(NA1 + N - 1)
        D_SUM = D_SUM + D_TN(N)*RCKWRK(NA1 + N - 1)
     end DO poly_loop
     SMH(K) = SUM + RCKWRK(NA1 + NCP2 - 1) - RCKWRK(NA1 + NCP1 - 1)/T
     D_SMH(K) = D_SUM + RCKWRK(NA1 + NCP1 - 1)*RINVTsq
!
  END DO species_loop
  RETURN
END SUBROUTINE JAC_CKSMH
!
!----------------------------------------------------------------------C
!
SUBROUTINE CKYTCP (P, T, Y, ICKWRK, RCKWRK, C)
!
!  START PROLOGUE
!
!  SUBROUTINE CKYTCP (P, T, Y, ICKWRK, RCKWRK, C)
!     Returns the molar concentrations given the pressure,
!     temperature and mass fractions;  see Eq. (7).
!
!  INPUT
!     P      - Pressure.
!                   cgs units - dynes/cm**2
!                   Data type - real scalar
!     T      - Temperature.
!                   cgs units - K
!                   Data type - real scalar
!     Y      - Mass fractions of the species.
!                   cgs units - none
!                   Data type - real array
!                   Dimension Y(*) at least KK, the total number of
!                   species.
!     ICKWRK - Array of integer workspace.
!                   Data type - integer array
!                   Dimension ICKWRK(*) at least LENIWK.
!     RCKWRK - Array of real work space.
!                   Data type - real array
!                   Dimension RCKWRK(*) at least LENRWK.
!
!  OUTPUT
!     C      - Molar concentrations of the species.
!                   cgs units - mole/cm**3
!                   Data type - real array
!                   Dimension C(*) at least KK, the total number of
!                   species.
!
!  END PROLOGUE
!
  IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
!
  INCLUDE 'ckstrt.fh'
!
  DIMENSION ICKWRK(*), RCKWRK(*), Y(*), C(*)
!
  SUMYOW = 0.0
  DO  K = 1, NKK
     SUMYOW = SUMYOW + Y(K) / RCKWRK(NcWT + K - 1)
  END DO
  SUMYOW = SUMYOW * T * RCKWRK(NcRU)
  DO  K = 1, NKK
     C(K) = P * Y(K) / (SUMYOW * RCKWRK(NcWT + K - 1))
  END DO
  RETURN
END SUBROUTINE CKYTCP
!
!----------------------------------------------------------------------C
!
SUBROUTINE CKRATX (II, KK, MAXSP, MAXTB, T, C, NSPEC, NU, NUNK,&
     NPAR, PAR, NFAL, IFAL, IFOP, KFAL, NFAR, FPAR, &
     NTHB, ITHB, NTBS, AIK, NKTB, RKFT, RKRT, RKF, &
     RKR, CTB, NRNU, IRNU, RNU, NORD, IORD, MXORD,&
     KORD, RORD)
!
!  START PROLOGUE
!
!  SUBROUTINE CKRATX (II, KK, MAXSP, MAXTB, T, C, NSPEC, NU, NUNK,
! 1                   NPAR, PAR, NFAL, IFAL, IFOP, KFAL, NFAR, FPAR, 
! 2                   NTHB, ITHB, NTBS, AIK, NKTB, RKFT, RKRT, RKF, 
! 3                   RKR, CTB, NRNU, IRNU, RNU, NORD, IORD, MXORD,
! 4                   KORD, RORD)
!
!  END PROLOGUE
!
  IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
!
  DIMENSION C(*), NSPEC(*), NU(MAXSP,*), NUNK(MAXSP,*), PAR(NPAR,*),&
       IFAL(*), IFOP(*), KFAL(*), FPAR(NFAR,*), ITHB(*),&
       NTBS(*), AIK(MAXTB,*), NKTB(MAXTB,*), RKFT(*),&
       RKRT(*), RKF(*), RKR(*), CTB(*), IRNU(*), RNU(MAXSP,*),&
       IORD(*), KORD(MXORD,*), RORD(MXORD,*)
!
  COMMON /MACH/ SMALL,BIG,EXPARG
!
  DO I = 1, II
     CTB(I) = 1D0
     RKF(I) = 0D0
     RKR(I) = 0D0
  END DO
!
!     third-body reactions
!
  IF (NTHB .GT. 0) THEN
     CTOT = 0.0
     DO  K = 1, KK
        CTOT = CTOT + C(K)
     END DO
     N_loop : DO  N = 1, NTHB
        CTB(ITHB(N)) = CTOT
        L_loop : DO  L = 1, NTBS(N)
           CTB(ITHB(N)) = CTB(ITHB(N)) + (AIK(L,N)-1.0)*C(NKTB(L,N))
        END DO L_loop
     END DO N_loop
  ENDIF
!
!     If fall-off (pressure correction):
!
  IF (NFAL .GT. 0) THEN
     ALOGT = LOG(T)
!
     DO N = 1, NFAL
!
        RKLOW = FPAR(1,N) * EXP(FPAR(2,N)*ALOGT - FPAR(3,N)/T)
!
!        CONCENTRATION OF THIRD BODY
!
        IF (KFAL(N) .EQ. 0) THEN
           PR = RKLOW * CTB(IFAL(N)) / RKFT(IFAL(N))
           CTB(IFAL(N)) = 1D0
        ELSE
           PR = RKLOW * C(KFAL(N)) / RKFT(IFAL(N))
        ENDIF
!
        PCOR = PR / (1D0 + PR)
!
        IF (IFOP(N) .GT. 1) THEN
           PRLOG = LOG10(MAX(PR,SMALL))
!
           IF (IFOP(N) .EQ. 2) THEN
!
!              8-PARAMETER SRI FORM
!
              XP = 1.0/(1.0 + PRLOG**2)
              FC = ((FPAR(4,N)*EXP(-FPAR(5,N)/T) &
                   + EXP(-T/FPAR(6,N))) **XP)&
                   * FPAR(7,N) * T**FPAR(8,N)
!
           ELSE
!
!              6-PARAMETER TROE FORM
!
              IF(FPAR(6,N) > SMALL) THEN
                 FCENT = FPAR(4,N) *  EXP(-T/FPAR(6,N))
              ELSE
                 FCENT = 0D0
              END IF
              IF(FPAR(5,N) > SMALL) FCENT = FCENT + (1.D0-FPAR(4,N)) * EXP(-T/FPAR(5,N))
!
!              7-PARAMETER TROE FORM
!
              IF (IFOP(N) .EQ. 4) FCENT = FCENT + EXP(-FPAR(7,N)/T)
!
              FCLOG = LOG10(MAX(FCENT,SMALL))
              XN    = 0.75 - 1.27*FCLOG
              CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
              FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
              FC = 10.0**FLOG
           ENDIF
           PCOR = FC * PCOR
        ENDIF
!
        RKFT(IFAL(N)) = RKFT(IFAL(N)) * PCOR
        RKRT(IFAL(N)) = RKRT(IFAL(N)) * PCOR
     END DO
  ENDIF
!
!     Multiply by the product of reactants and product of products
!     PAR(4,I) is a perturbation factor
!
  DO  I = 1, II
     RKFT(I) = RKFT(I) * CTB(I) * PAR(4,I)
     RKRT(I) = RKRT(I) * CTB(I) * PAR(4,I)
!
     IF (NU(1,I) .NE. 0) THEN
        RKF(I) = RKFT(I)*C(NUNK(1,I))**IABS(NU(1,I))
        RKR(I) = RKRT(I)*C(NUNK(4,I))**NU(4,I)
        IF (NUNK(2,I) .NE. 0) THEN
           RKF(I)= RKF(I) * C(NUNK(2,I))**IABS(NU(2,I))
           IF (NUNK(3,I) .NE. 0)&
                RKF(I) = RKF(I) * C(NUNK(3,I))**IABS(NU(3,I))
        ENDIF
        IF (NUNK(5,I) .NE. 0) THEN
           RKR(I) = RKR(I) * C(NUNK(5,I))**NU(5,I)
           IF (NUNK(6,I) .NE. 0) RKR(I)=RKR(I)*C(NUNK(6,I))**NU(6,I)
        ENDIF
     ENDIF
  END DO
!
  DO  N = 1, NRNU
     I = IRNU(N)
     C1 = C(NUNK(1,I)) ** ABS(RNU(1,N))
     C4 = C(NUNK(4,I)) ** RNU(4,N)
     RKF(I) = RKFT(I) * C1
     RKR(I) = RKRT(I) * C4
     IF (NUNK(2,I) .NE. 0) THEN
        C2 = C(NUNK(2,I)) ** ABS(RNU(2,N))
        RKF(I) = RKF(I) * C2
        IF (NUNK(3,I) .NE. 0) THEN
           C3 = C(NUNK(3,I)) ** ABS(RNU(3,N))
           RKF(I) = RKF(I) * C3
        ENDIF
     ENDIF
     IF (NUNK(5,I) .NE. 0) THEN
        C5 = C(NUNK(5,I)) ** RNU(5,N)
        RKR(I) = RKR(I) * C5
        IF (NUNK(6,I) .NE. 0) THEN
           C6 = C(NUNK(6,I)) ** RNU(6,N)
           RKR(I) = RKR(I) * C6
        ENDIF
     ENDIF
  END DO
!
  N2_loop: DO  N = 1, NORD
     I = IORD(N)
     RKF(I) = RKFT(I)
     RKR(I) = RKRT(I)
!
     L2_loop: DO L = 1, MXORD
        NK = KORD(L,N)
        IF (NK .LT. 0) THEN
           NK = IABS(NK)
           CNK = C(NK) ** RORD(L,N)
           RKF(I) = RKF(I) * CNK
        ELSEIF (NK .GT. 0) THEN
           CNK = C(NK) ** RORD(L,N)
           RKR(I) = RKR(I) * CNK
        ENDIF
     END DO L2_loop
  END DO N2_loop
!
  RETURN
END SUBROUTINE CKRATX
!
!----------------------------------------------------------------------C
!
SUBROUTINE CKGAM (T,X, ICKWRK, RCKWRK, GAM)
!
!  START PROLOGUE
!
!  SUBROUTINE CKCVBLOR (T,X, ICKWRK, RCKWRK, GAM)
!     Returns GAMMA
!
!  INPUT
!     T      - Temperature.
!                   cgs units - K
!                   Data type - real scalar
!     X      - mole fractions
!
!     ICKWRK - Array of integer workspace.
!                   Data type - integer array
!                   Dimension ICKWRK(*) at least LENIWK.
!     RCKWRK - Array of real work space.
!                   Data type - real array
!                   Dimension RCKWRK(*) at least LENRWK.
!
!  OUTPUT
!     GAM   - Gamma index for the speed of sound         
!                   cgs units - none
!                   Data type - scalar
!
!  END PROLOGUE
!
  IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
!
  INCLUDE 'ckstrt.fh'
!
  DIMENSION ICKWRK(*), RCKWRK(*), X(*), TN(10)
!
  TN(1) = 1.0
  DO  N = 2, NCP
     TN(N) = T**(N-1)
  END DO
!
  CVBLOR = -1d0
  species_loop: DO  K = 1, NKK
     L = 1
     range_loop: DO N = 2, ICKWRK(IcNT + K - 1)-1
        TEMP = RCKWRK(NcTT + (K-1)*MXTP + N - 1)
        IF (T .GT. TEMP) L = L+1
     end DO range_loop
!
     NA1 = NcAA + (L-1)*NCP2 + (K-1)*NCP2T
     sum = 0d0
     poly_loop: DO  N = 1, NCP
        sum = sum + TN(N)*RCKWRK(NA1 + N - 1)
     end DO poly_loop
     CVBLOR = CVBLOR + sum*X(K)
  end DO species_loop
  GAM = 1d0+1d0/CVBLOR
!
  RETURN
END SUBROUTINE CKGAM

!
!----------------------------------------------------------------------C
!
SUBROUTINE CKCPBMS (T, Y, ICKWRK, RCKWRK, CPBMS)
!
!  START PROLOGUE
!
!  SUBROUTINE CKCPBMS (T, Y, ICKWRK, RCKWRK, CPBMS)
!     Returns the mean specific heat at constant pressure;
!     see Eq. (34).
!
!  INPUT
!     T      - Temperature.
!                   cgs units - K
!                   Data type - real scalar
!     Y      - Mass fractions of the species.
!                   cgs units - none
!                   Data type - real array
!                   Dimension Y(*) at least KK, the total number of
!                   species.
!     ICKWRK - Array of integer workspace.
!                   Data type - integer array
!                   Dimension ICKWRK(*) at least LENIWK.
!     RCKWRK - Array of real work space.
!                   Data type - real array
!                   Dimension RCKWRK(*) at least LENRWK.
!
!  OUTPUT
!     CPBMS  - Mean specific heat at constant pressure in mass units.
!                   cgs units - ergs/(gm*K)
!                   Data type - real scalar
!
!  END PROLOGUE
!
  IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
!
  INCLUDE 'ckstrt.fh'
!
  DIMENSION ICKWRK(*), RCKWRK(*), Y(*),TN(10)
!
  TN(1) = 1.0
  DO  N = 2, NCP
     TN(N) = T**(N-1)
  END DO
!
  CPBMS = 0d0
  species_loop: DO  K = 1, NKK
     L = 1
     range_loop: DO N = 2, ICKWRK(IcNT + K - 1)-1
        TEMP = RCKWRK(NcTT + (K-1)*MXTP + N - 1)
        IF (T .GT. TEMP) L = L+1
     end DO range_loop
!
     NA1 = NcAA + (L-1)*NCP2 + (K-1)*NCP2T
     termK = RCKWRK(NcRU)/RCKWRK(NcWT + K - 1)
     sum = 0d0
     poly_loop: DO  N = 1, NCP
        sum = sum + TN(N)*RCKWRK(NA1 + N - 1)
     end DO poly_loop
     CPBMS = CPBMS + termK*sum*Y(K)
  end DO species_loop

  return
END SUBROUTINE CKCPBMS
!
!----------------------------------------------------------------------
!
SUBROUTINE JAC_CKCpMS (T,ICKWRK, RCKWRK, CPK,CPKDER,Hms)
!
!  START PROLOGUE
!
!  SUBROUTINE CKCPBMS (T, Y, ICKWRK, RCKWRK, CPBMS)
!     Returns the specific heat at constant pressure;
!     the derivative and the inegral for all species
!
!  INPUT
!     T      - Temperature.
!                   cgs units - K
!                   Data type - real scalar
!     ICKWRK - Array of integer workspace.
!                   Data type - integer array
!                   Dimension ICKWRK(*) at least LENIWK.
!     RCKWRK - Array of real work space.
!                   Data type - real array
!                   Dimension RCKWRK(*) at least LENRWK.
!
!  OUTPUT
!     CPK  -       specific heat at constant pressure
!                   cgs units - ergs/(gm*K)
!                   Data type - real size NKK
!
!  END PROLOGUE
!
  IMPLICIT NONE
!
  INCLUDE 'ckstrt.fh'
!
  INTEGER :: K,L,N,na1
  INTEGER :: ICKWRK(*)
  REAL*8 :: RCKWRK(*), CPK(*),CPKDER(*),hms(*)
  REAL*8 :: TN(10),TN_D(10),TN_I(10)
  REAL*8 :: T,temp,termK,sum,sum_I,sum_D
!
  TN(1) = 1D0
  TN_D(1) = 0d0
  TN_I(1) = T
  DO  N = 2, NCP
     TN(N) = T**(N-1)
     TN_D(N) = (N-1)*TN(N)/T
     TN_I(N) = TN(N)*T/N
  END DO
!
  species_loop: DO  K = 1, NKK
     L = 1
     range_loop: DO N = 2, ICKWRK(IcNT + K - 1)-1
        TEMP = RCKWRK(NcTT + (K-1)*MXTP + N - 1)
        IF (T .GT. TEMP) L = L+1
     end DO range_loop
!
     NA1 = NcAA + (L-1)*NCP2 + (K-1)*NCP2T
     termK = RCKWRK(NcRU)/RCKWRK(NcWT + K - 1)
     sum   = 0d0
     sum_I = 0d0
     sum_D = 0d0
     poly_loop: DO  N = 1, NCP
        sum   = sum   + TN(N)  *RCKWRK(NA1 + N - 1)
        sum_I = sum_I + TN_I(N)*RCKWRK(NA1 + N - 1)
        sum_D = sum_D + TN_D(N)*RCKWRK(NA1 + N - 1)
     end DO poly_loop
     CPK(k)   =  termK*sum
     Cpkder(k) =  termK*sum_D
     hms(k)    =  termK*(sum_I + RCKWRK(NA1 + NCP1 - 1))
  end DO species_loop

  return
END SUBROUTINE JAC_CKCPMS
!
!----------------------------------------------------------------------
!
SUBROUTINE MCADIF (P, T, X, RMCWRK, D)
!
  IMPLICIT real*8 (A-H, O-Z), INTEGER (I-N)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  SUBROUTINE MCADIF (P, T, X, RMCWRK, D)
!
!    THIS SUBROUTINE COMPUTES MIXTURE-AVERAGED DIFFUSION COEFFICIENTS
!    GIVEN THE PRESSURE, TEMPERATURE, AND SPECIES MASS FRACTIONS.
!
!  INPUT-
!    P       - PRESSURE
!                  CGS UNITS - DYNES/CM**2.
!    T       - TEMPERATURE
!                  CGS UNITS - K.
!    X       - ARRAY OF MOLE FRACTIONS OF THE MIXTURE.
!                  DIMENSION X(*) AT LEAST KK.
!
!  WORK-
!    RMCWRK  - ARRAY OF FLOATING POINT STORAGE AND WORK SPACE.  THE
!              STARTING ADDRESSES FOR THE RMCWRK SPACE ARE STORED IN
!              COMMON /MCMCMC/.
!                  DIMENSION RMCWRK(*) AT LEAST LENRMC.
!
!  OUTPUT-
!    D       - ARRAY OF MIXTURE DIFFUSION COEFFICIENTS
!                  CGS UNITS - CM**2/S.
!                  DIMENSION D(*) AT LEAST KK.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
  DIMENSION X(*), D(*), RMCWRK(*)
!
  COMMON /MCMCMC/ RU, PATMOS, SMALL, NKK, NO, NLITE, INLIN, IKTDIF,&
       IPVT, NWT, NEPS, NSIG, NDIP, NPOL, NZROT, NLAM,&
       NETA, NDIF, NTDIF, NXX, NVIS, NXI, NCP, NCROT,&
       NCINT, NBIND, NEOK, NSGM, NAST, NBST,&
       NCST, NXL, NR, NWRK, K3
!
  CALL MCEDIF (T, NO, NKK, X, RMCWRK(NDIF), RMCWRK(NWT), SMALL,&
       RMCWRK(NXX), RMCWRK(NBIND), D)
!
  DO  K = 1, NKK
     D(K) = D(K) * PATMOS/P
  end DO
!
  RETURN
END SUBROUTINE MCADIF
!
!----------------------------------------------------------------------
!
SUBROUTINE MCEDIF(T, NO, KK, X, COFD, WT, SMALL, XX, DJK, D)
!
  !USE MYMPI, ONLY : MYID
  IMPLICIT real*8 (A-H, O-Z), INTEGER (I-N)
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  SUBROUTINE MCEDIF(T, NO, KK, X, COFD, WT, SMALL, XX, DJK, D)
!
!    THIS SUBROUTINE IS USED INTERNALLY TO COMPUTE THE MIXTURE
!    DIFFUSION COEFFICIENTS.  NORMALLY NOT CALLED BY THE PACKAGE USER.
!
!  INPUT-
!    T       - TEMPERATURE
!                  CGS UNITS - K.
!    NO      - ORDER OF FIT.
!    KK      - NUMBER OF SPECIES.
!    X       - ARRAY OF MOLE FRACTIONS OF THE MIXTURE.
!                  DIMENSION X(*) AT LEAST KK.
!    COFD    - COEFFICIENTS OF THE FITS FOR THE BINARY DIFFUSION
!              COEFFICIENTS.
!                  DIMENSION COFD(NO,KK,*) EXACTLY NO FOR
!                   FIRST DIMENSION, KK FOR THE SECOND, AND AT
!                   LEAST KK FOR THE THIRD.
!    WT      - ARRAY OF SPECIES MOLECULAR WEIGHTS.
!                   DIMENSION WT(*) AT LEAST KK.
!    SMALL   - A SMALL NUMBER ADDED TO ALL MOLE FRACTIONS BEFORE
!                COMPUTING THE MIXTURE DIFFUSION COEFFICIENTS.
!                THIS PROCESS AVOIDS AN UNDEFINED SITUATION WHEN
!                A PURE SPECIES CONDITION IS APPROACHED.
!    XX      - ARRAY OF MOLE FRACTIONS PLUS "SMALL," TO AVOID THE
!              PROBLEM OF A PURE SPECIES.
!                  DIMENSION XX(*) AT LEAST KK.
!
!  WORK-
!    RMCWRK  - ARRAY OF FLOATING POINT STORAGE AND WORK SPACE.  THE
!              STARTING ADDRESSES FOR THE RMCWRK SPACE ARE STORED IN
!              COMMON /MCMCMC/.
!                  DIMENSION RMCWRK(*) AT LEAST LENRMC.
!
!  OUTPUT-
!    D       - ARRAY OF MIXTURE DIFFUSION COEFFICIENTS
!                  CGS UNITS - CM**2/S.
!                  DIMENSION D(*) AT LEAST KK.
!    DJK     - MATRIX OF BINARY DIFFUSION COEFFICIENTS.  DJK(J,K) IS
!              DIFFUSION COEFFICIENT OF SPECIES J IN SPECIES K.
!                  CGS UNITS - CM**2/S
!                  DIMENSION DJK(KDIM,*) EXACTLY KDIM FOR THE FIRST
!                   DIMENSION AND AT LEAST KK FOR THE SECOND.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
  DIMENSION X(*), COFD(NO,KK,*), WT(*), XX(*), DJK(KK,*), D(*)
!
  ZERO = 0.0d0
  ALOGT = LOG(T)
  DO  K = 1, KK
     CALL MCEVAL (ALOGT, KK, NO, COFD(1,1,K), DJK(1,K) )
  END DO

!Use a pade approximation of the exponential function to save time
  DO  K = 1, KK
     DO  J = 1, KK
        xpd = DJK(J,K)
        if(xpd > 3d0) then   !error of 0.1%%
           DJK(J,K) = EXP(DJK(J,K))
        else
           xpd2 = xpd*xpd
           xpd3 = xpd2*xpd
           xpd4 = xpd2*xpd2
           DJK(J,K) = (1680d0+840d0*xpd+180d0*xpd2+20d0*xpd3+xpd4)/&
                &     (1680d0-840d0*xpd+180d0*xpd2-20d0*xpd3+xpd4)
        end if
     end DO
  end DO

  if (KK == 1) then
      D(1) = DJK(1,1)
      return
   end if
!
  WTM = ZERO
  DO  K = 1, KK
     WTM = WTM + WT(K)*X(K)
     XX(K) = X(K) + SMALL
  end DO
  SUMALL = SUM(XX(1:KK)*WT(1:KK))
!
  spec_loop: DO  K = 1, KK
!
     SUMXW = SUMALL - XX(K)*WT(K)
     SUMXOD = ZERO
!
     sum_loop: DO  J = 1, KK
        IF (J .NE. K) THEN
           SUMXOD = SUMXOD + XX(J)/DJK(J,K)
        ENDIF
     end DO sum_loop
!
     D(K) = SUMXW/(WTM*SUMXOD)
!
  end DO spec_loop
!
  RETURN
END SUBROUTINE MCEDIF
!
!---------------------------------------------------------------------
!

SUBROUTINE MCMCDT (P, T, X, IMCWRK, RMCWRK, ICKWRK, CKWRK,&
     DT, COND)

  IMPLICIT Real*8 (A-H, O-Z), INTEGER (I-N)
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!   SUBROUTINE MCMCDT (P, T, X, KDIM, IMCWRK, RMCWRK, ICKWRK, CKWRK,
!  1                   DT, COND)
!
!    THIS SUBROUTINE COMPUTES THE THERMAL DIFFUSION COEFFICIENTS, AND
!    MIXTURE THERMAL CONDUCTIVITIES, GIVEN THE PRESSURE, TEMPERATURE,
!    AND MOLE FRACTIONS.
!
!  INPUT-
!    P       - PRESSURE
!                  CGS UNITS - DYNES/CM**2.
!    T       - TEMPERATURE
!                  CGS UNITS - K.
!    X       - ARRAY OF MOLE FRACTIONS OF THE MIXTURE.
!                  DIMENSION X(*) AT LEAST KK.
!
!  WORK-
!    IMCWRK  - ARRAY OF INTEGER STORAGE AND WORK SPACE.  THE
!              STARTING ADDRESSES FOR THE IMCWRK SPACE ARE STORED IN
!              COMMON /MCMCMC/.
!                  DIMENSION IMCWRK(*) AT LEAST LENIMC.
!    RMCWRK  - ARRAY OF FLOATING POINT STORAGE AND WORK SPACE.  THE
!              STARTING ADDRESSES FOR THE RMCWRK SPACE ARE STORED IN
!              COMMON /MCMCMC/.
!                  DIMENSION RMCWRK(*) AT LEAST LENRMC.
!    ICKWRK  - ARRAY OF INTEGER CHEMKIN STORAGE AND WORK SPACE.  SEE
!                  CHEMKIN DOCUMENTATION.
!                  DIMENSION ICKWRK(*) AT LEAST LENIWRK.
!    CKWRK   - ARRAY OF FLOATING POINT CHEMKIN STORAGE AND WORK SPACE.
!              SEE CHEMKIN DOCUMENTATION.
!                  DIMENSION CKWRK(*) AT LEAST LENWRK.
!
!  OUTPUT-
!    DT      - VECTOR OF THERMAL MULTICOMPONENT DIFFUSION COEFFICIENTS.
!                  CGS UNITS - GM/(CM*SEC)
!    COND    - MIXTURE THERMAL CONDUCTIVITY
!                  CGS UNITS - ERG/(CM*K*S).
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
  DIMENSION X(*), IMCWRK(*), RMCWRK(*), ICKWRK(*), CKWRK(*), DT(*)
!
  COMMON /MCMCMC/ RU, PATMOS, SMALL, NKK, NO, NLITE, INLIN, IKTDIF,&
       IPVT, NWT, NEPS, NSIG, NDIP, NPOL, NZROT, NLAM,&
       NETA, NDIF, NTDIF, NXX, NVIS, NXI, NCP, NCROT,&
       NCINT, NBIND, NEOK, NSGM, NAST, NBST,&
       NCST, NXL, NR, NWRK, K3
!
  CALL MCLMDT (P, T, X, NKK, NO, NETA, K3, SMALL, RMCWRK(NWT), RMCWRK(NEOK),&
       RMCWRK(NZROT), IMCWRK(INLIN), RMCWRK(NEPS),&
       ICKWRK, CKWRK, RMCWRK, RMCWRK(NXX), RMCWRK(NVIS),&
       RMCWRK(NAST), RMCWRK(NBST), RMCWRK(NCST),&
       RMCWRK(NXI),  RMCWRK(NCP), RMCWRK(NCROT),&
       RMCWRK(NCINT), RMCWRK(NXL), RMCWRK(NR), &
       RMCWRK(NBIND), IMCWRK(IPVT), DT, COND)
!
  RETURN
END SUBROUTINE MCMCDT

!
!---------------------------------------------------------------------
!
SUBROUTINE MCLMDT (P, T, X, KK, NO, NETA, KK3, SMALL, WT, EOK, ZROT, LIN,&
     EPS, ICKWRK, CKWRK, RMCWRK, XX, VIS, ASTAR,&
     BSTAR, CSTAR, XI, CPOR, CROTOR, CINTOR, XL,&
     R, BINDIF, IPVT, DT, COND)

  IMPLICIT Real*8 (A-H, O-Z), INTEGER (I-N)
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  SUBROUTINE MCLMDT (P, T, X, KK, KK3, SMALL, WT, EOK, ZROT,
! 1                   LIN, EPS, ICKWRK, CKWRK, RMCWRK, XX, VIS,
! 2                   ASTAR, BSTAR, CSTAR, XI, CPOR, CROTOR,
! 3                   CINTOR, XL, R, BINDIF, IPVT,
! 4                   DT, COND)
!
!
!    THIS SUBROUTINE COMPUTES THE THERMAL CONDUCTIVITY, AND THE THERMAL
!    DIFFUSION COEFFICIENT ARRAY.  IT DOES SO BY FIRST FORMING THE L
!    MATRIX, AND THEN SOLVING EQ. 24A.  THIS ROUTINE IS NOT NORMALLY CAL
!    DIRECTLY BY THE USER; THE USER CALLS MCMCDT, WHICH IN TURN CALLS
!
!  INPUT-
!    P       - PRESSURE
!                  CGS UNITS - DYNES/CM**2.
!    T       - TEMPERATURE
!                  CGS UNITS - K.
!    X       - ARRAY OF MOLE FRACTIONS OF THE MIXTURE.
!                  DIMENSION X(*) AT LEAST KK.
!    KK      - NUMBER OF SPECIES.
!    KK3     - THREE TIMES THE NUMBER OF SPECIES.  THE SIZE OF THE L
!              MATRIX IS KK3 * KK3.
!    SMALL   - THE MOLE FRACTIONS USED IN THE TRANSPORT COMPUTATION
!              ARE GIVEN BY XX(K) = X(K) + SMALL.
!    WT      - THE ARRAY OF SPECIES MOLECULAR WEIGHTS.
!                   DIMENSION WT(*) AT LEAST KK.
!    EOK     - MATRIX OF REDUCED WELL DEPTHS FOR EACH SPECIES PAIR.
!                  UNITS - K
!                  DIMENSION EOK(KK,*) EXACTLY KDIM FOR THE FIRST
!                    DIMENSION, AND AT LEAST KK FOR THE SECOND.
!    ZROT    - ARRAY OF ROTATIONAL COLLISION NUMBERS EVALUATED
!              AT 298K.
!                  UNITS - NONE
!                  DIMENSION ZROT(*) AT LEAST KK
!    LIN     - ARRAY OF FLAGS INDICATING WHETHER THE MOLECULE
!              LINEAR OR NOT.
!              NLIN=0, SINGLE ATOM.
!              NLIN=1, LINEAR MOLECULE.
!              NLIN=2, NONLINEAR MOLECULE.
!                  UNITS - NONE.
!                  DIMENSION NLIN(*) AT LEAST KK
!    EPS     - ARRAY OF LENNARD-JONES POTENTIAL WELL DEPTHS.
!                  CGS UNITS - K.
!                  DIMENSION EPS(*) AT LEAST KK
!
!  WORK AND SCRATCH SPACE
!    ICKWRK  - ARRAY OF CHEMKIN INTEGER WORK SPACE.
!    CKWRK   - ARRAY OF CHEMKIN REAL WORK SPACE.
!    RMCWRK  - ARRAY OF FLOATING POINT STORAGE AND WORK SPACE.  THE
!              STARTING ADDRESSES FOR THE RMCWRK SPACE ARE STORED IN
!              COMMON /MCMCMC/.
!                  DIMENSION RMCWRK(*) AT LEAST LENRMC.
!    XX      - THE SPECIES MOLE FRACTION ARRAY THAT IS USED IN THE
!              TRANSPORT COMPUTATION.  XX(K) = X(K) + SMALL.
!    VIS     - ARRAY OF SPECIES VISCOSITIES.  EVALUATED FROM MCSVIS.
!                  CGS UNITS - GM/CM-S
!                  DIMENSION VIS(*) AT LEAST KK.
!    ASTAR   - MATRIX OF COLLISION INTEGRALS A*, FOR EACH SPECIES PAIR.
!                  DIMENSION ASTAR(KDIM,*) EXACTLY KDIM FOR THE FIRST
!                   DIMENSION AND AT LEAST KK FOR THE SECOND.
!    BSTAR   - MATRIX OF COLLISION INTEGRALS B*, FOR EACH SPECIES PAIR.
!                  DIMENSION BSTAR(KDIM,*) EXACTLY KDIM FOR THE FIRST
!                   DIMENSION AND AT LEAST KK FOR THE SECOND.
!    CSTAR   - MATRIX OF COLLISION INTEGRALS C*, FOR EACH SPECIES PAIR.
!                  DIMENSION CSTAR(KDIM,*) EXACTLY KDIM FOR THE FIRST
!                   DIMENSION AND AT LEAST KK FOR THE SECOND.
!    XI      - ARRAY COLLISION NUMBERS FOR THE TRANSFER OF ROTATIONAL
!              ENERGY OF SPECIES I INTO TRANSLATIONAL ENERGY OF
!              SPECIES J (EQ. 42).  WE ASSUME THAT ALL XI(I,J) = XI(I,I)
!              SEE P. 132 FOR DISCUSSION.
!    CPOR    - ARRAY OF DIMENSIONLESS SPECIFIC HEATS, CP/R.  EVALUATED
!              FROM CKCPOR.
!                   DIMENSION CPOR(*) AT LEAST KK.
!    CROT    - ARRAY OF DIMENSIONLESS ROTATIONAL CONTRIBUTIONS TO THE
!              SPECIES SPECIFIC HEATS.
!                   DIMENSION CROT(*) AT LEAST KK.
!    CINT    - ARRAY OF DIMENSIONLESS INTERNAL CONTRIBUTIONS TO THE
!              SPECIES SPECIFIC HEATS.
!                   DIMENSION CINT(*) AT LEAST KK.
!    XL      - THE L MATRIX, EQ. 43 AND 49.
!                  DIMENSION XL(3*KK,*) EXACTLY KDIM FOR THE FIRST
!                   DIMENSION AND AT LEAST 3*KK FOR THE SECOND.
!    R       - ARRAY OF RIGHT HAND SIDES OF EQ. 24A.
!                   DIMENSION R(*) AT LEAST 3*KK.
!    BINDIF  - MATRIX OF BINARY DIFFUSION COEFFICIENTS.
!                  CGS UNITS - CM**2/S
!                  DIMENSION BINDIF(KK,*) EXACTLY KDIM FOR THE FIRST
!                   DIMENSION AND AT LEAST KK FOR THE SECOND.
!    IPVT    - ARRAY OF PIVOT INDICES FOR THE INVERSION OF THE XL
!              MATRIX BY THE LINPACK ROUTINES SGEFA AND SGEDI.
!                  DIMENSION IPVT(*) AT LEAST 3*KK.
!
!  OUTPUT-
!    DT      - VECTOR OF THERMAL MULTICOMPONENT DIFFUSION COEFFICIENTS.
!                  CGS UNITS - GM/(CM*SEC)
!                  DIMENSION DT(*) AT LEAST KK.
!    COND    - MIXTURE THERMAL CONDUCTIVITY
!                  CGS UNITS - ERG/(CM*K*S).
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
  DIMENSION X(*), ICKWRK(*), CKWRK(*), RMCWRK(*), WT(*), XX(*),&
       VIS(*), EOK(KK,*), ZROT(*), LIN(*), EPS(*),&
       ASTAR(KK,*), BSTAR(KK,*), CSTAR(KK,*), XI(*), CPOR(*),&
       CROTOR(*), CINTOR(*), XL(KK3,*), R(*), BINDIF(KK,*),&
       IPVT(*), DT(*), FITAST(7), FITBST(7), FITCST(7)
!
!            FITS OF A*, B*, AND C* AS FUNCTIONS OF LN(T*)
!
  DATA FITAST / .1106910525E+01, -.7065517161E-02,&
       -.1671975393E-01,  .1188708609E-01,&
       .7569367323E-03, -.1313998345E-02,&
       .1720853282E-03/
!
  DATA FITBST / .1199673577E+01, -.1140928763E+00,&
       -.2147636665E-02,  .2512965407E-01,&
       -.3030372973E-02, -.1445009039E-02,&
       .2492954809E-03/
!
  DATA FITCST / .8386993788E+00,  .4748325276E-01,&
       .3250097527E-01, -.1625859588E-01,&
       -.2260153363E-02,  .1844922811E-02,&
       -.2115417788E-03/
!
  DATA RU/8.314E+07/, PI/3.1415926535/ 
  DATA PI32O2/2.7842/, P2O4P2/4.4674/, PI32/5.5683/
!
  ONE = 1.0D0
  ZERO = 0.0D0
!
!         SET MINIMUM MOLE FRACTION TO SMALL
!
  DO  K = 1, KK
     XX(K) = X(K) +  SMALL
  end DO
!
!          DETERMINE A*, B*, AND C* FOR EACH SPECIES
!
  DO J = 1, KK
     DO  I = 1, KK
!
        TSLOG = LOG ( T/EOK(I,J) )
        T1 = TSLOG
        T2 = TSLOG*T1
        T3 = TSLOG*T2
        T4 = TSLOG*T3
        T5 = TSLOG*T4
        T6 = TSLOG*T5
        ASTAR(I,J) = FITAST(1)    + FITAST(2)*T1 + FITAST(3)*T2 +&
             FITAST(4)*T3 + FITAST(5)*T4 + FITAST(6)*T5 +&
             FITAST(7)*T6
        BSTAR(I,J) = FITBST(1)    + FITBST(2)*T1 + FITBST(3)*T2 +&
             FITBST(4)*T3 + FITBST(5)*T4 + FITBST(6)*T5 +&
             FITBST(7)*T6
        CSTAR(I,J) = FITCST(1)    + FITCST(2)*T1 + FITCST(3)*T2 +&
             FITCST(4)*T3 + FITCST(5)*T4 + FITCST(6)*T5 +&
             FITCST(7)*T6
     end DO
  end DO
!
!        EVALUATE THE BINARY DIFFUSION COEFFICIENTS AND VISCOSITY
!
  CALL MCSDIF (P, T, KK, RMCWRK, BINDIF)
!
  ALOGT = LOG(T)
  CALL MCEVAL (ALOGT, KK, NO, RMCWRK(NETA), VIS)
  Kvis: DO  K = 1, KK
     VIS(K) = EXP(VIS(K))
  end DO Kvis
!
  TRU = 6.0 * RU * T
  PFAC = 5.0 * P
  DO  K = 1, KK
!
!        EVALUATE BINARY SELF-DIFFUSION COEFFICIENTS FROM VISCOSITY
!
     BINDIF(K,K) = TRU * ASTAR(K,K) * VIS(K) / ( PFAC * WT(K) )
!
!         COMPUTE PARKER CORRECTION FOR ZROT
!
     DD = EPS(K) / T
     DR = EPS(K) / 298.0
     SQRTDD = SQRT(DD)
     SQRTDR = SQRT(DR)
     DD32 = SQRTDD*DD
     DR32 = SQRTDR*DR
     XI(K) = ( (ONE + PI32O2*SQRTDR + P2O4P2*DR + PI32*DR32) /&
          (ONE + PI32O2*SQRTDD + P2O4P2*DD + PI32*DD32)  )&
          * MAX(ONE, ZROT(K))
  end DO
!
!         ROTATIONAL AND INTERNAL PARTS OF SPECIFIC HEAT
!
  CALL CKCPOR (T, ICKWRK, CKWRK, CPOR)

  DO  K = 1, KK
     IF (LIN(K) .EQ. 0) THEN
        CROTOR(K) = ZERO
        CINTOR(K) = ZERO
     ELSEIF (LIN(K) .EQ. 1) THEN
        CROTOR(K) = ONE
        CINTOR(K) = CPOR(K) - 2.5
     ELSEIF (LIN(K) .EQ. 2) THEN
        CROTOR(K) = 3.0/2.0
        CINTOR(K) = CPOR(K) - 2.5
     ENDIF
  end DO
!
!        ASSEMBLE L00,00
!
  PFAC = 16.0 * T / (25.0 * P)
  jloop:  DO  J = 1, KK
     iloop:    DO  I = 1, KK
        SUM = ZERO
        kloop:     DO  K = 1, KK
           SUM = SUM + XX(K) / (WT(I)*BINDIF(I,K)) *&
                ( WT(J)*XX(J) * (ONE - DELTA(I,K)) -&
                WT(I)*XX(I) * (DELTA(I,J) - DELTA(J,K)) )
        end DO kloop
        XL(I,J) = PFAC * SUM
     end DO iloop
  end DO jloop
!
!         ASSEMBLE L00,10
!
  PFAC = 8.0 * T / (5.0 * P)
  jloop1: DO  J = 1, KK
     iloop1:    DO  I = 1, KK
        SUM = ZERO
        kloop1:    DO  K = 1, KK
           SUM = SUM + XX(J)*XX(K) * (DELTA(I,J)-DELTA(I,K)) *&
                WT(K) * (1.2*CSTAR(J,K) - ONE) /&
                ( (WT(J)+WT(K)) * BINDIF(J,K) )
        end DO kloop1
        XL(I, J+KK) = PFAC * SUM
     end DO iloop1
  end DO jloop1
!
!         ASSEMBLE L10,00
!
  DO J = 1, KK
     DO  I = 1, KK
        XL(I+KK, J) = XL(J, I+KK)
     end DO
  end DO
!
!         ASSEMBLE L01,00 AND L00,01
!
  DO  J = 1, KK
     DO  I = 1, KK
        XL(2*KK+I, J) = ZERO
        XL(I, 2*KK+J) = ZERO
     end DO
  end DO
!
!         ASSEMBLE L10,10
!
  PFAC = 16.0 * T / (25.0 * P)
  PIFAC = 5.0 / (3.0*PI)
  jloop2: DO  J = 1, KK
     iloop2: DO  I = 1, KK
        SUM = ZERO
        kloop2:  DO  K = 1, KK
           SUM = SUM + XX(I)*XX(K) /&
                ( (WT(I)+WT(K))**2 * BINDIF(I,K) ) *&
                (  (DELTA(J,K)-DELTA(I,J)) *&
                ( 7.5*WT(J)**2 + 6.25*WT(K)**2 -&
                3.0*WT(K)**2 * BSTAR(I,K) )  -&
                4.0*WT(J)*WT(K) * ASTAR(I,K) *&
                (DELTA(J,K)+DELTA(I,J)) *&
                ( ONE + PIFAC *&
                (CROTOR(I)/XI(I) + CROTOR(K)/XI(K)) ) )
        end DO kloop2
        XL(I+KK, J+KK) = PFAC * WT(I) * SUM / WT(J)
     end DO iloop2
  end DO jloop2
!
!         ASSEMBLE DIAGONAL OF L10,10 USING VISCOSITY, EQ. (49)
!
  PFAC = 16.0 * T / (25.0 * P)
  PIFAC = 5.0 / (3.0 * PI)
  iloop3: DO  I = 1, KK
     SUM = ZERO
     kloop3: DO  K = 1, KK
        IF (I .NE. K) THEN
           SUM = SUM + XX(I)*XX(K) /&
                ( (WT(I)+WT(K))**2 * BINDIF(I,K) ) *&
                ( ( 7.5*WT(I)**2 + 6.25*WT(K)**2 -&
                3.0*WT(K)**2 * BSTAR(I,K) )  +&
                4.0*WT(I)*WT(K) * ASTAR(I,K) *&
                ( ONE + PIFAC *&
                (CROTOR(I)/XI(I) + CROTOR(K)/XI(K)) ) )
        ENDIF
     end DO kloop3
     XL(I+KK, I+KK) = - (16.0/15.0) *&
          ( WT(I) * XX(I)**2 / (RU * VIS(I) ) ) *&
          ( ONE + 2.0 * PIFAC *&
          (CROTOR(I)/XI(I)) ) - PFAC*SUM
  end DO iloop3
!
!         ASSEMBLE L10,01 AND L01,10
!
  NN = 2*KK
  PFAC = 32.0d0 * T / (5.0d0 * PI * P)
  jloop4: DO  J = 1, KK
     IF (LIN(J) .NE. 0) THEN
        NN = NN + 1
        iloop4:   DO  I = 1, KK
           SUM = ZERO
           kloop4: DO  K = 1, KK
              SUM = SUM + WT(J) * ASTAR(J,K)&
                   * (DELTA(I,K)+DELTA(I,J)) *&
                   XX(J) * XX(K) * CROTOR(J) /&
                   ( (WT(J)+WT(K)) * BINDIF(J,K) * XI(J) )
           end DO kloop4
!              L10,01
           XL(I+KK, NN) = PFAC * SUM / CINTOR(J) 
!              L01,10
           XL(NN, I+KK) = XL(I+KK, NN)
        end DO iloop4
     ENDIF
  end DO jloop4
!
!         ASSEMBLE L10,01 AND L01,10 USING VISCOSITY, EQ. (49)
!
  NN = 2*KK
  PFAC = 32.0 * T / (5.0 * PI * P)
  PIFAC = 16.0 / (3.0 * PI)
  iloop5: DO  I = 1, KK
     IF (LIN(I) .NE. 0) THEN
        NN = NN + 1
        SUM = ZERO
        kloop5: DO K = 1, KK
           IF (I .NE. K) THEN
              SUM = SUM + WT(I) * ASTAR(I,K) *&
                   XX(I) * XX(K) * CROTOR(I) /&
                   ( (WT(I)+WT(K)) * BINDIF(I,K) * XI(I) )
           ENDIF
        end DO kloop5
!           L10,01
        XL(I+KK, NN) = PIFAC *&
             ( WT(I) * XX(I)**2 * CROTOR(I) /&
             (RU * VIS(I) * CINTOR(I) * XI(I) ) )+&
             PFAC * SUM / CINTOR(I) 
!           L01,10
        XL(NN, I+KK) = XL(I+KK, NN)
     ENDIF
  end DO iloop5
!
!        ASSEMBLE L01,01, USING VISCOSITY EQ. (49)
!
  DO  J = 1, KK
     DO  I = 1, KK
        XL(2*KK+I, 2*KK+J) = ZERO
     end DO
  end DO
!
  NN = 2*KK
  PFAC = 4.0 * T / P
  PIFAC = 5.0 * PI
  PIRU  = PI * RU
  iloop6: DO  I = 1, KK
     IF (LIN(I) .NE. 0) THEN
        NN = NN + 1
        SUM = ZERO
        kloop6: DO K = 1, KK
           SUM = SUM + XX(I)*XX(K) / BINDIF(I,K)
           IF (I .NE. K) THEN
              SUM = SUM +&
                   12.0 * XX(I)*XX(K) * WT(I) * ASTAR(I,K) *&
                   CROTOR(I) /&
                   ( PIFAC * CINTOR(I) * WT(K) * BINDIF(I,K)&
                   * XI(I)  )
           ENDIF
        end DO kloop6
        XL(NN, NN) = - 8.0 * WT(I) * XX(I)**2 * CROTOR(I) /&
             ( PIRU * CINTOR(I)**2 * VIS(I) * XI(I) ) -&
             PFAC * SUM / CINTOR(I)
     ENDIF
  end DO iloop6
!
!          ASSEMBLE THE RIGHT HAND SIDE FOR SOLVING EQ. (24)
!
  NN = 2*KK
  DO  I = 1, KK
     R(I)    = ZERO
     R(I+KK) = XX(I)
     IF (LIN(I) .NE. 0) THEN
        NN  = NN + 1
        R(NN) = XX(I)
     ENDIF
  end DO
!
!          FACTOR AND SOLVE EQ. (24)
!
  CALL DGEFA (XL, KK3, NN, IPVT, INFO)

  IF (INFO .NE. 0) THEN
     WRITE (6,*) ' ERROR IN SGEFA, INFO = ', INFO
     STOP
  ENDIF
!
! Solve the linear system of equations
!
  CALL DGESL (XL, KK3, NN, IPVT, R, 0)
!
!          FORM THERMAL DIFFUSION COEFFICIENTS
!
  RFAC = 8.0d0 / (5.0d0 * RU)
  DO K = 1, KK
     DT(K) = RFAC * WT(K) * XX(K) * R(K)
  end DO
!
!          FORM THE THERMAL CONDUCTIVITY
!
  CONDTR = ZERO
  DO  K = 1, KK
     CONDTR = CONDTR + X(K) * R(KK+K)
  end DO
  CONDTR = - 4.0 * CONDTR
!
  NN = 2*KK
  CONDIN = ZERO
  DO  K = 1, KK
     IF (LIN(K) .NE. 0) THEN
        NN = NN + 1
        CONDIN = CONDIN + X(K) * R(NN)
     ENDIF
  END DO

  CONDIN = - 4.0d0 * CONDIN
!
  COND = CONDTR + CONDIN
!
  RETURN
END SUBROUTINE MCLMDT
!
!----------------------------------------------------------------------
!

!
!----------------------------------------------------------------------C
!
SUBROUTINE CKCPOR (T, ICKWRK, RCKWRK, CPOR)
!
!  START PROLOGUE
!
!  SUBROUTINE CKCPOR (T, ICKWRK, RCKWRK, CPOR)
!     Returns the nondimensional specific heats at constant pressure;
!     see Eq. (19).
!
!  INPUT
!     T      - Temperature.
!                   cgs units - K
!                   Data type - real scalar
!     ICKWRK - Array of integer workspace.
!                   Data type - integer array
!                   Dimension ICKWRK(*) at least LENIWK.
!     RCKWRK - Array of real work space.
!                   Data type - real array
!                   Dimension RCKWRK(*) at least LENRWK.
!
!  OUTPUT
!     CPOR   - Nondimensional specific heats at constant pressure
!              for the species.
!                   cgs units - none
!                   Data type - real array
!                   Dimension CPOR(*) at least KK, the total number of
!                   species.
!
!  END PROLOGUE
!
  IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
!
  INCLUDE 'ckstrt.fh'
!
  DIMENSION ICKWRK(*), RCKWRK(*), TN(10), CPOR(*)
!
  TN(1) = 1.0
  DO  N = 2, NCP
     TN(N) = T**(N-1)
  end DO
!
  specloop: DO K = 1, NKK
     L = 1
     DO  N = 2, ICKWRK(IcNT + K - 1)-1
        TEMP = RCKWRK(NcTT + (K-1)*MXTP + N - 1)
        IF (T .GT. TEMP) L = L+1
     end DO
!
     NA1 = NcAA + (L-1)*NCP2 + (K-1)*NCP2T
     CPOR(K) = 0.0
     polyloop:   DO  N = 1, NCP
        CPOR(K) = CPOR(K) + TN(N)*RCKWRK(NA1 + N - 1)
     end DO polyloop
  end DO specloop
  RETURN
END SUBROUTINE CKCPOR
!
!----------------------------------------------------------------------C
!
SUBROUTINE CKSNUM (LINE, NEXP, LOUT, KRAY, NN, KNUM, NVAL,&
     RVAL, KERR)
!
!  START PROLOGUE
!
!  SUBROUTINE CKSNUM (LINE, NEXP, LOUT, KRAY, NN, KNUM, NVAL,
!                     RVAL, KERR)
!     This subroutine is called to parse a character string, LINE,
!     that is composed of several blank-delimited substrings.
!     It is expected that the first substring in LINE is also an
!     entry in a reference array of character strings, KRAY(*), in
!     which case the index position in KRAY(*) is returned as KNUM,
!     otherwise an error flag is returned.  The substrings following
!     the first are expected to represent numbers, and are converted
!     to elements of the array RVAL(*).  If NEXP substrings are not
!     found an error flag will be returned.  This allows format-free
!     input of combined alpha-numeric data.  For example, after
!     reading a line containing a species name followed by several
!     numerical values, the subroutine might be called to find
!     a Chemkin species index and convert the other substrings to
!     real values:
!
!     input:  LINE    = "N2  1.2"
!             NEXP    = 1, the number of values expected
!             LOUT    = 6, a logical unit number on which to write
!                       diagnostic messages.
!             KRAY(*) = "H2" "O2" "N2" "H" "O" "N" "OH" "H2O" "NO"
!             NN      = 9, the number of entries in KRAY(*)
!     output: KNUM    = 3, the index number of the substring in
!                       KRAY(*) which corresponds to the first
!                       substring in LINE
!             NVAL    = 1, the number of values found in LINE
!                       following the first substring
!             RVAL(*) = 1.200E+00, the substring converted to a number
!             KERR    = .FALSE.
!  INPUT
!     LINE   - A character string.
!                   Data type - CHARACTER*80
!     NEXP   - Number of real values to be found in character string.
!              If NEXP is negative, then IABS(NEXP) values are
!              expected.  However, it is not an error condition,
!              if less values are found.
!                   Data type - integer scalar
!     LOUT   - Output unit for printed diagnostics.
!                   Data type - integer scalar
!     KRAY   - Array of character strings.
!                   Data type - CHARACTER*(*)
!     NN     - Total number of character strings in KRAY.
!                   Data type - integer scalar
!
!  OUTPUT
!     KNUM   - Index number of character string in array which
!              corresponds to the first substring in LINE.
!                   Data type - integer scalar
!     NVAL   - Number of real values found in LINE.
!                   Data type - integer scalar
!     RVAL   - Array of real values found in LINE.
!                   Data type - real array
!                   Dimension RVAL(*) at least NEXP
!     KERR   - Error flag; syntax or dimensioning error,
!              corresponding string not found, or total of
!              values found is not the number of values expected,
!              will result in KERR = .TRUE.
!                   Data type - logical
!
!  END PROLOGUE
!     A '!' will comment out a line, or remainder of the line.
!
  IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
!
  CHARACTER LINE*(*), KRAY(*)*(*), ISTR*80
  DIMENSION RVAL(*)
  LOGICAL KERR, IERR
!
  NVAL = 0
  KERR = .FALSE.
  ILEN = MIN (IPPLEN(LINE), ILASCH(LINE))
  IF (ILEN .LE. 0) RETURN
!
  I1 = IFIRCH(LINE(:ILEN))
  I3 = INDEX(LINE(I1:ILEN),' ')
  IF (I3 .EQ. 0) I3 = ILEN - I1 + 1
  I2 = I1 + I3
  ISTR = ' '
  ISTR = LINE(I1:I2-1)
!
  CALL CKCOMP (ISTR, KRAY, NN, KNUM)
  IF (KNUM.EQ.0) THEN
     LT = MAX (ILASCH(ISTR), 1)
     WRITE (LOUT,'(A)')&
          ' Error in CKSNUM...'//ISTR(:LT)//' not found...'
     KERR = .TRUE.
  ENDIF
!
  ISTR = ' '
  ISTR = LINE(I2:ILEN)
  IF (NEXP .NE. 0)&
       CALL CKXNUM (ISTR, NEXP, LOUT, NVAL, RVAL, IERR)
!
  RETURN
END SUBROUTINE CKSNUM
!
!----------------------------------------------------------------------C
!
SUBROUTINE CKCOMP (IST, IRAY, II, I)
!
!  START PROLOGUE
!
!  SUBROUTINE CKCOMP (IST, IRAY, II, I)*
!     Returns the index of an element of a reference character
!     string array which corresponds to a character string;
!     leading and trailing blanks are ignored.
!
!
!  INPUT
!     IST   - A character string.
!                  Data type - CHARACTER*(*)
!     IRAY  - An array of character strings;
!                  Data type - CHARACTER*(*)
!                  Dimension IRAY(*) at least II
!     II    - The length of IRAY.
!                  Data type - integer scalar.
!
!  OUTPUT
!     I     - The first integer location in IRAY in which IST
!             corresponds to IRAY(I); if IST is not also an
!             entry in IRAY, I=0.
!
!  END PROLOGUE
!
!*****precision > double
  IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
!*****END precision > double
!*****precision > single
!      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
!*****END precision > single
!
  CHARACTER*(*) IST, IRAY(*)
!
  I = 0
  DO N = II, 1, -1
     IS1 = IFIRCH(IST)
     IS2 = ILASCH(IST)
     IR1 = IFIRCH(IRAY(N))
     IR2 = ILASCH(IRAY(N))
     IF ( IS2.GE.IS1 .AND. IS2.GT.0 .AND.&
          IR2.GE.IR1 .AND. IR2.GT.0 .AND.&
          IST(IS1:IS2).EQ.IRAY(N)(IR1:IR2) ) I = N
  end DO
  RETURN
END SUBROUTINE CKCOMP
!
!----------------------------------------------------------------------C
!
SUBROUTINE CKXNUM (LINE, NEXP, LOUT, NVAL, RVAL, KERR)
!
!  START PROLOGUE
!
!  SUBROUTINE CKXNUM (LINE, NEXP, LOUT, NVAL, RVAL, KERR)
!     This subroutine is called to parse a character string, LINE,
!     that is composed of several blank-delimited substrings.
!     Each substring is expected to represent a number, which
!     is converted to entries in the array of real numbers, RVAL(*).
!     NEXP is the number of values expected, and NVAL is the
!     number of values found.  This allows format-free input of
!     numerical data.  For example:
!
!     input:  LINE    = " 0.170E+14 0 47780.0"
!             NEXP    = 3, the number of values requested
!             LOUT    = 6, a logical unit number on which to write
!                       diagnostic messages.
!     output: NVAL    = 3, the number of values found
!             RVAL(*) = 1.700E+13, 0.000E+00, 4.778E+04
!             KERR    = .FALSE.
!
!  INPUT
!     LINE   - A character string.
!                   Data type - CHARACTER*80
!     NEXP   - Number of real values to be found in character string.
!              If NEXP is negative, then IABS(NEXP) values are
!              expected.  However, it is not an error condition,
!              if less values are found.
!                   Data type - integer scalar
!     LOUT   - Output unit for printed diagnostics.
!                   Data type - integer scalar
!
!  OUTPUT
!     NVAL   - Number of real values found in character string.
!                   Data type - integer scalar
!     RVAL   - Array of real values found.
!                   Data type - real array
!                   Dimension RVAL(*) at least NEXP
!     KERR   - Error flag;  syntax or dimensioning error results
!              in KERR = .TRUE.
!                   Data type - logical
!
!  END PROLOGUE
!
!     A '!' will comment out a line, or remainder of the line.
!
!*****precision > double
  IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
!*****END precision > double
!*****precision > single
!        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
!*****END precision > single
!
  CHARACTER LINE*(*), ITEMP*80
  DIMENSION RVAL(*), RTEMP(80)
  LOGICAL KERR
!
!----------Find Comment String (! signifies comment)
!
  ILEN = IPPLEN(LINE)
  NVAL = 0
  KERR = .FALSE.
!
  IF (ILEN .LE. 0) RETURN
  IF (ILEN .GT. 80) THEN
     WRITE (LOUT,*)     ' Error in CKXNUM...line length > 80 '
     WRITE (LOUT,'(A)') LINE
     KERR = .TRUE.
     RETURN
  ENDIF
!
  ITEMP = LINE(:ILEN)
  IF (NEXP .LT. 0) THEN
     CALL IPPARR (ITEMP, -1, NEXP, RTEMP, NVAL, IERR, LOUT)
  ELSE
     CALL IPPARR (ITEMP, -1, -NEXP, RTEMP, NVAL, IERR, LOUT)
     IF (IERR .EQ. 1) THEN
        WRITE (LOUT, *)    ' Syntax errors in CKXNUM...'
        WRITE (LOUT,'(A)') LINE
        KERR = .TRUE.
     ELSEIF (NVAL .NE. NEXP) THEN
        WRITE (LOUT,*) ' Error in CKXNUM...'
        WRITE (LOUT,'(A)') LINE
        KERR = .TRUE.
        WRITE (LOUT,*) NEXP,' values expected, ',&
             NVAL,' values found.'
     ENDIF
  ENDIF
  IF (NVAL .LE. IABS(NEXP)) THEN
     DO N = 1, NVAL
        RVAL(N) = RTEMP(N)
     end DO
  ENDIF
!
  RETURN
END SUBROUTINE CKXNUM
!
!----------------------------------------------------------------------C
!
SUBROUTINE CKELCF(MDIM, ICKWRK, RCKWRK, ELCF)
!
!  START PROLOGUE
!
!  SUBROUTINE CKELCF  (MDIM, ICKWRK, RCKWRK, ELCF)
!     Returns the elemental composition of the species
!
!  INPUT
!     MDIM   - First dimension of the two-dimensional array ELCF;
!              MDIM must be equal to or greater than the number of
!              elements, MM.
!                   Data type - integer scalar
!     ICKWRK - Array of integer workspace.
!                   Data type - integer array
!                   Dimension ICKWRK(*) at least LENIWK.
!     RCKWRK - Array of real work space.
!                   Data type - real array
!                   Dimension RCKWRK(*) at least LENRWK.
!
!  OUTPUT
!     ELCF    - Matrix of the elemental composition of the species;
!              ELCF(M,K) is the number of atoms of the Mth element
!              in the Kth species.
!                   Data type - integer array
!                   Dimension ELCF(MDIM,*) exactly MDIM (at least MM,
!                   the total number of elements in the problem) for
!                   the first dimension and at least KK, the total
!                   number of species, for the second.
!
!  END PROLOGUE
!
        IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
!
      DIMENSION ICKWRK(*), RCKWRK(*), ELCF(MDIM,*)
!
      INCLUDE 'ckstrt.fh'
!
      DO K = 1, NKK
         J = IcNC + (K-1)*NMM
         DO M = 1, NMM
            ELCF(M,K) = dble(ICKWRK(J + M - 1))
         END DO
      END DO
      RETURN
    END SUBROUTINE CKELCF
!
!----------------------------------------------------------------------C
!
SUBROUTINE IPPARR (STRING,ICARD,NEXPEC,RVAL,NFOUND,IERR,LOUT)
!   BEGIN PROLOGUE  IPPARR
!   REFER TO  IPGETR
!   DATE WRITTEN  850625   (YYMMDD)
!   REVISION DATE 851625   (YYMMDD)
!   CATEGORY NO.  J3.,J4.,M2.
!   KEYWORDS  PARSE
!   AUTHOR  CLARK,G.L.,GROUP C-3 LOS ALAMOS NAT'L LAB
!   PURPOSE  Parses real variables from a character variable.  Called
!            by IPGETR, the IOPAK routine used for interactive input.
!   DESCRIPTION
!
!-----------------------------------------------------------------------
!  IPPARR may be used for parsing an input record that contains real
!  values, but was read into a character variable instead of directly
!  into real variables.
!  The following benefits are gained by this approach:
!    - specification of only certain elements of the array is allowed,
!      thus letting the others retain default values
!    - variable numbers of values may be input in a record, up to a
!      specified maximum
!    - control remains with the calling program in case of an input
!      error
!    - diagnostics may be printed by IPPARR to indicate the nature
!      of input errors
!
!   The contents of STRING on input indicate which elements of RVAL
!   are to be changed from their entry values, and values to which
!   they should be changed on exit.  Commas and blanks serve as
!   delimiters, but multiple blanks are treated as a single delimeter.
!   Thus, an input record such as:
!     '   1.,   2,,4.e-5   , ,6.e-6'
!   is interpreted as the following set of instructions by IPGETR:
!
!     (1) set RVAL(1) = 1.0
!     (2) set RVAL(2) = 2.0
!     (3) leave RVAL(3) unchanged
!     (4) set RVAL(4) = 4.0E-05
!     (5) leave RVAL(5) unchanged
!     (6) set RVAL(6) = 6.0E-06
!
!   IPPARR will print diagnostics on the default output device, if
!   desired.
!
!   IPPARR is part of IOPAK, and is written in ANSI FORTRAN 77
!
!   Examples:
!
!      Assume RVAL = (0., 0., 0.) and NEXPEC = 3 on entry:
!
!   input string           RVAL on exit            IERR    NFOUND
!   -------------          ----------------------  ----    ------
!  '  2.34e-3,  3 45.1'    (2.34E-03, 3.0, 45.1)     0       3
!  '2,,3.-5'               (2.0, 0.0, 3.0E-05)       0       3
!  ',1.4,0.028E4'          (0.0, 1.4, 280.0)         0       3
!  '1.0, 2.a4, 3.0'        (1.0, 0.0, 0.0)           1       1
!  '1.0'                   (1.0, 0.0, 0.0)           2       1
!
!      Assume RVAL = (0.,0.,0.,0.) and NEXPEC = -4 on entry:
!
!   input string           RVAL on exit            IERR    NFOUND
!   -------------          ----------------------  ----    ------
!  '1.,2.'                 (1.0, 2.0)                0       2
!  ',,3  4.0'              (0.0, 0.0, 3.0, 4.0)      0       4
!  '1,,3,,5.0'             (0.0, 0.0, 3.0, 0.0)      3       4
!
!  arguments: (I=input,O=output)
!  -----------------------------
!  STRING (I) - the character string to be parsed.
!
!  ICARD  (I) - data statement number, and error processing flag
!         < 0 : no error messages printed
!         = 0 : print error messages, but not ICARD
!         > 0 : print error messages, and ICARD
!
!  NEXPEC (I) - number of real variables expected to be input.  If
!         < 0, the number is unknown, and any number of values
!         between 0 and abs(nexpec) may be input.  (see NFOUND)
!
!  PROMPT (I) - prompting string, character type.  A question
!         mark will be added to form the prompt at the screen.
!
!  RVAL (I,O) - the real value or values to be modified.  On entry,
!       the values are printed as defaults.  The formal parameter
!       corresponding to RVAL must be dimensioned at least NEXPEC
!       in the calling program if NEXPEC > 1.
!
!  NFOUND (O) - the number of real values represented in STRING,
!         only in the case that there were as many or less than
!         NEXPEC.
!
!  IERR (O) - error flag:
!       = 0 if no errors found
!       = 1 syntax errors or illegal values found
!       = 2 for too few values found (NFOUND < NEXPEC)
!       = 3 for too many values found (NFOUND > NEXPEC)
!-----------------------------------------------------------------------
!
!   REFERENCES  (NONE)
!   ROUTINES CALLED  IFIRCH,ILASCH
!   END PROLOGUE  IPPARR
!*****precision > double
  IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
!*****END precision > double
!*****precision > single
!      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
!*****END precision > single
!
  CHARACTER STRING*(*), ITEMP*80
  DIMENSION RVAL(*)
  CHARACTER *8 FMT(16)
  LOGICAL OKINCR
!
!   FIRST EXECUTABLE STATEMENT  IPPARR
  IERR   = 0
  NFOUND = 0
  NEXP = IABS(NEXPEC)
  IE = ILASCH(STRING)
  IF (IE .EQ. 0) GO TO 500
  NC = 1
!
!--- OKINCR is a flag that indicates it's OK to increment
!--- NFOUND, the index of the array into which the value
!--- should be read.  It is set negative when a space follows
!--- a real value substring, to keep incrementing from
!--- occurring if a comma should be encountered before the
!--- next value.
!
  OKINCR = .TRUE.
!
!--- begin overall loop on characters in string
!
100 CONTINUE
!
  IF (STRING(NC:NC) .EQ. ',') THEN
     IF (OKINCR) THEN
        NFOUND = NFOUND + 1
     ELSE
        OKINCR = .TRUE.
     ENDIF
!
     GO TO 450
  ENDIF
  IF (STRING(NC:NC) .EQ. ' ') GO TO 450
!
!--- first good character (non-delimeter) found - now find
!--- last good character
!
  IBS = NC
160 CONTINUE
  NC = NC + 1
  IF (NC .GT. IE) GO TO 180
  IF (STRING(NC:NC) .EQ. ' ')THEN
     OKINCR = .FALSE.
  ELSEIF (STRING(NC:NC) .EQ. ',')THEN
     OKINCR = .TRUE.
  ELSE
     GO TO 160
  ENDIF
!
!--- end of substring found - read value into real array
!
180 CONTINUE
  NFOUND = NFOUND + 1
  IF (NFOUND .GT. NEXP) THEN
     IERR = 3
     GO TO 500
  ENDIF
!
  DATA FMT/     ' (E1.0)', ' (E2.0)', ' (E3.0)', ' (E4.0)',&
       ' (E5.0)', ' (E6.0)', ' (E7.0)', ' (E8.0)', ' (E9.0)',&
       '(E10.0)', '(E11.0)', '(E12.0)', '(E13.0)', '(E14.0)',&
       '(E15.0)', '(E16.0)'/
  IES = NC - 1
  NCH = IES - IBS + 1
  ITEMP = ' '
  ITEMP = STRING(IBS:IES)
  READ (ITEMP(:NCH), FMT(NCH), ERR = 400) RVAL(NFOUND)
  GO TO 450
400 CONTINUE
  IERR = 1
  GO TO 510
450 CONTINUE
  NC = NC + 1
  IF (NC .LE. IE) GO TO 100
!
500 CONTINUE
  IF (NEXPEC .GT. 0 .AND. NFOUND .LT. NEXP) IERR = 2
510 CONTINUE
!
  IF (IERR .EQ. 0 .OR. ICARD .LT. 0) RETURN
  IF (ICARD .NE. 0) WRITE(LOUT,'(A,I3)')&
       '!! ERROR IN DATA STATEMENT NUMBER', ICARD
  IF (IERR .EQ. 1) WRITE(LOUT,'(A)')'SYNTAX ERROR, OR ILLEGAL VALUE'
  IF (IERR .EQ. 2) WRITE(LOUT,'(A,I2, A, I2)')&
       ' TOO FEW DATA ITEMS.  NUMBER FOUND = ' , NFOUND,&
       '  NUMBER EXPECTED = ', NEXPEC
  IF (IERR .EQ. 3) WRITE(LOUT,'(A,I2)')&
       ' TOO MANY DATA ITEMS.  NUMBER EXPECTED = ', NEXPEC
END SUBROUTINE IPPARR
!
!----------------------------------------------------------------------C
!
FUNCTION IPPLEN (LINE)
!
!  BEGIN PROLOGUE
!
!  FUNCTION IPPLEN (LINE)
!     Returns the effective length of a character string, i.e.,
!     the index of the last character before an exclamation mark (!)
!     indicating a comment.
!
!  INPUT
!     LINE  - A character string.
!                  Data type - CHARACTER*(*)
!
!  OUTPUT
!     IPPLEN - The effective length of the character string.
!                   Data type - integer scalar
!
!  END PROLOGUE
!
  IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
!
  CHARACTER LINE*(*)
!
  IN = IFIRCH(LINE)
  IF (IN.EQ.0 .OR. LINE(IN:IN) .EQ. '!') THEN
     IPPLEN = 0
  ELSE
     IN = INDEX(LINE,'!')
     IF (IN .EQ. 0) THEN
        IPPLEN = ILASCH(LINE)
     ELSE
        IPPLEN = ILASCH(LINE(:IN-1))
     ENDIF
  ENDIF
  RETURN
END FUNCTION IPPLEN
!
!----------------------------------------------------------------------C
!
FUNCTION IFIRCH   (STRING)
!   BEGIN PROLOGUE  IFIRCH
!   DATE WRITTEN   850626
!   REVISION DATE  850626
!   CATEGORY NO.  M4.
!   KEYWORDS  CHARACTER STRINGS,SIGNIFICANT CHARACTERS
!   AUTHOR  CLARK,G.L.,GROUP C-3 LOS ALAMOS NAT'L LAB
!   PURPOSE  Determines first significant (non-blank) character
!            in character variable
!   DESCRIPTION
!
!-----------------------------------------------------------------------
!  IFIRCH locates the first non-blank character in a string of
!  arbitrary length.  If no characters are found, IFIRCH is set = 0.
!  When used with the companion routine ILASCH, the length of a string
!  can be determined, and/or a concatenated substring containing the
!  significant characters produced.
!-----------------------------------------------------------------------
!
!   REFERENCES  (NONE)
!   ROUTINES CALLED  (NONE)
!   END PROLOGUE IFIRCH
  IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
!
  CHARACTER* (*)STRING
!
!   FIRST EXECUTABLE STATEMENT IFIRCH
  NLOOP = LEN(STRING)
!
  IF (NLOOP.EQ.0 .OR. STRING.EQ.' ') THEN
     IFIRCH = 0
     RETURN
  ENDIF
!
  DO I = 1, NLOOP
     IF (STRING(I:I) .NE. ' ') GO TO 120
  end DO
!
  IFIRCH = 0
  RETURN
120 CONTINUE
  IFIRCH = I
END FUNCTION IFIRCH
!
!
!----------------------------------------------------------------------C
!
subroutine readfile(fname,KK,ksym,Bout)
  IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
  character*(*) :: fname
  character :: ksym(KK)*16
  real*8 :: bout(*),val(10)
  character :: line*80
  logical kerr,ierr

  LIN = 599


  open (unit=LIN, STATUS='OLD', file=TRIM(fname))


40 CONTINUE
  READ  (LIN,7600,END=45)LINE

  ILEN = INDEX (LINE, '!')
  IF (ILEN .EQ. 1) GO TO 40

  ILEN = ILEN - 1
  IF (ILEN .LE. 0) ILEN = LEN(LINE)
  IF (INDEX(LINE(:ILEN), 'END') .EQ. 0) THEN
     IF (LINE(:ILEN) .NE. ' ') THEN
        CALL CKSNUM (LINE(:ILEN), 1, LOUT, KSYM, KK, KNUM,&
             NVAL, VAL, IERR)
        IF (IERR) THEN
           WRITE (LOUT,*) ' Error reading moles...'
           KERR = .TRUE.
        ELSE
           bout(KNUM) = VAL(1)
        ENDIF
     ENDIF
     GO TO 40
  ENDIF
!
  CLOSE(LIN)
!    
45 CONTINUE

7600 FORMAT (1A80)

  return
end subroutine readfile
!
!*************************************************************
!
!************************************************************************************
!
real*8 FUNCTION  FALLOFF_PCOR(T,PR,NFPAR,FPAR,ifop) Result(PCOR)
  implicit NONE
  integer,INTENT(IN) ::IFOP,NFPAR
  real*8,INTENT(IN) :: T,PR,Fpar(NFPAR)
  real*8 :: PRLOG,XP,FC,FCENT,FCLOG,XN,CPRLOG,FLOG
!
  real*8 :: SMALL,BIG,EXPARG
  COMMON /MACH/ SMALL,BIG,EXPARG
!
  PCOR = PR / (1D0 + PR)

  IF (IFOP .GT. 1) THEN
     PRLOG = LOG10(MAX(PR,SMALL))
!
     IF (IFOP .EQ. 2) THEN
!
!              8-PARAMETER SRI FORM
!
        XP = 1d0/(1d0 + PRLOG**2)
        FC = ((FPAR(4)*EXP(-FPAR(5)/T) &
             + EXP(-T/FPAR(6))) **XP)&
             * FPAR(7) * T**FPAR(8)
!
     ELSE
!
!              6-PARAMETER TROE FORM
!
        IF(FPAR(6) > SMALL) THEN
           FCENT = FPAR(4) *  EXP(-T/FPAR(6))
        ELSE
           FCENT = 0D0
        END IF
        IF(FPAR(5) > SMALL) FCENT = FCENT + (1.D0-FPAR(4)) * EXP(-T/FPAR(5))
!
!              7-PARAMETER TROE FORM
!
        IF (IFOP .EQ. 4) FCENT = FCENT + EXP(-FPAR(7)/T)
!
        FCLOG = LOG10(MAX(FCENT,SMALL))
        XN    = 0.75d0 - 1.27d0*FCLOG
        CPRLOG= PRLOG - (0.4d0 + 0.67d0*FCLOG)
        FLOG = FCLOG/(1d0 + (CPRLOG/(XN-0.14d0*CPRLOG))**2)
        FC = 10d0**FLOG
     ENDIF
     PCOR = FC * PCOR
  ENDIF

end FUNCTION FALLOFF_PCOR


!
!************************************************************
!
!----------------------------------------------------------------------C
real*8 function functionT(U,b_298,T) Result(fT)
  implicit none
  INCLUDE 'ckstrt.fh'
  real*8,INTENT(IN) :: U,b_298(*),T
  integer :: k

  fT = b_298(1)*T - U
  do k = 2,NCP
    fT = fT + b_298(k)*T**k
  enddo

end function functionT
!----------------------------------------------------------------------C

