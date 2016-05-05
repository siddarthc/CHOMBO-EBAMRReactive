#ifndef _EBPLANARREACTIVESHOCKF_F_H_
#define _EBPLANARREACTIVESHOCKF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSETPLANARSHOCK 
#define GUARDSETPLANARSHOCK 
// Prototype for Fortran procedure setplanarshock ...
//
void FORTRAN_NAME( SETPLANARSHOCK ,setplanarshock )(
      CHFp_CONST_REAL(aushockvel)
      ,CHFp_CONST_REAL(aucenter)
      ,CHFp_CONST_REAL(aupreshockpress)
      ,CHFp_CONST_REAL(aupreshocktemp)
      ,CHFp_CONST_VR(aumassfrac)
      ,CHFp_CONST_INT(aushocknorm)
      ,CHFp_CONST_INT(aushockback) );

#define FORT_SETPLANARSHOCK FORTRAN_NAME( inlineSETPLANARSHOCK, inlineSETPLANARSHOCK)
#define FORTNT_SETPLANARSHOCK FORTRAN_NAME( SETPLANARSHOCK, setplanarshock)

inline void FORTRAN_NAME(inlineSETPLANARSHOCK, inlineSETPLANARSHOCK)(
      CHFp_CONST_REAL(aushockvel)
      ,CHFp_CONST_REAL(aucenter)
      ,CHFp_CONST_REAL(aupreshockpress)
      ,CHFp_CONST_REAL(aupreshocktemp)
      ,CHFp_CONST_VR(aumassfrac)
      ,CHFp_CONST_INT(aushocknorm)
      ,CHFp_CONST_INT(aushockback) )
{
 CH_TIMELEAF("FORT_SETPLANARSHOCK");
 FORTRAN_NAME( SETPLANARSHOCK ,setplanarshock )(
      CHFt_CONST_REAL(aushockvel)
      ,CHFt_CONST_REAL(aucenter)
      ,CHFt_CONST_REAL(aupreshockpress)
      ,CHFt_CONST_REAL(aupreshocktemp)
      ,CHFt_CONST_VR(aumassfrac)
      ,CHFt_CONST_INT(aushocknorm)
      ,CHFt_CONST_INT(aushockback) );
}
#endif  // GUARDSETPLANARSHOCK 

#ifndef GUARDSETSINEWAVE 
#define GUARDSETSINEWAVE 
// Prototype for Fortran procedure setsinewave ...
//
void FORTRAN_NAME( SETSINEWAVE ,setsinewave )(
      CHFp_CONST_REALVECT(aucenter)
      ,CHFp_CONST_REAL(aurnot)
      ,CHFp_CONST_REAL(aupress)
      ,CHFp_CONST_REAL(autemp)
      ,CHFp_CONST_VR(aumassfrac)
      ,CHFp_CONST_REAL(aumach)
      ,CHFp_CONST_REAL(ausmallness) );

#define FORT_SETSINEWAVE FORTRAN_NAME( inlineSETSINEWAVE, inlineSETSINEWAVE)
#define FORTNT_SETSINEWAVE FORTRAN_NAME( SETSINEWAVE, setsinewave)

inline void FORTRAN_NAME(inlineSETSINEWAVE, inlineSETSINEWAVE)(
      CHFp_CONST_REALVECT(aucenter)
      ,CHFp_CONST_REAL(aurnot)
      ,CHFp_CONST_REAL(aupress)
      ,CHFp_CONST_REAL(autemp)
      ,CHFp_CONST_VR(aumassfrac)
      ,CHFp_CONST_REAL(aumach)
      ,CHFp_CONST_REAL(ausmallness) )
{
 CH_TIMELEAF("FORT_SETSINEWAVE");
 FORTRAN_NAME( SETSINEWAVE ,setsinewave )(
      CHFt_CONST_REALVECT(aucenter)
      ,CHFt_CONST_REAL(aurnot)
      ,CHFt_CONST_REAL(aupress)
      ,CHFt_CONST_REAL(autemp)
      ,CHFt_CONST_VR(aumassfrac)
      ,CHFt_CONST_REAL(aumach)
      ,CHFt_CONST_REAL(ausmallness) );
}
#endif  // GUARDSETSINEWAVE 

#ifndef GUARDGETPOSTSHOCKPROPS 
#define GUARDGETPOSTSHOCKPROPS 
// Prototype for Fortran procedure getpostshockprops ...
//
void FORTRAN_NAME( GETPOSTSHOCKPROPS ,getpostshockprops )( );

#define FORT_GETPOSTSHOCKPROPS FORTRAN_NAME( inlineGETPOSTSHOCKPROPS, inlineGETPOSTSHOCKPROPS)
#define FORTNT_GETPOSTSHOCKPROPS FORTRAN_NAME( GETPOSTSHOCKPROPS, getpostshockprops)

inline void FORTRAN_NAME(inlineGETPOSTSHOCKPROPS, inlineGETPOSTSHOCKPROPS)( )
{
 CH_TIMELEAF("FORT_GETPOSTSHOCKPROPS");
 FORTRAN_NAME( GETPOSTSHOCKPROPS ,getpostshockprops )( );
}
#endif  // GUARDGETPOSTSHOCKPROPS 

#ifndef GUARDPLANARSHOCKINIT 
#define GUARDPLANARSHOCKINIT 
// Prototype for Fortran procedure planarshockinit ...
//
void FORTRAN_NAME( PLANARSHOCKINIT ,planarshockinit )(
      CHFp_FRA(u)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) );

#define FORT_PLANARSHOCKINIT FORTRAN_NAME( inlinePLANARSHOCKINIT, inlinePLANARSHOCKINIT)
#define FORTNT_PLANARSHOCKINIT FORTRAN_NAME( PLANARSHOCKINIT, planarshockinit)

inline void FORTRAN_NAME(inlinePLANARSHOCKINIT, inlinePLANARSHOCKINIT)(
      CHFp_FRA(u)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_PLANARSHOCKINIT");
 FORTRAN_NAME( PLANARSHOCKINIT ,planarshockinit )(
      CHFt_FRA(u)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDPLANARSHOCKINIT 

#ifndef GUARDPOINTPLANARSHOCKINIT 
#define GUARDPOINTPLANARSHOCKINIT 
// Prototype for Fortran procedure pointplanarshockinit ...
//
void FORTRAN_NAME( POINTPLANARSHOCKINIT ,pointplanarshockinit )(
      CHFp_REAL(density)
      ,CHFp_REALVECT(momentum)
      ,CHFp_REAL(energy)
      ,CHFp_VR(specDense)
      ,CHFp_CONST_INTVECT(iv)
      ,CHFp_CONST_REAL(dx) );

#define FORT_POINTPLANARSHOCKINIT FORTRAN_NAME( inlinePOINTPLANARSHOCKINIT, inlinePOINTPLANARSHOCKINIT)
#define FORTNT_POINTPLANARSHOCKINIT FORTRAN_NAME( POINTPLANARSHOCKINIT, pointplanarshockinit)

inline void FORTRAN_NAME(inlinePOINTPLANARSHOCKINIT, inlinePOINTPLANARSHOCKINIT)(
      CHFp_REAL(density)
      ,CHFp_REALVECT(momentum)
      ,CHFp_REAL(energy)
      ,CHFp_VR(specDense)
      ,CHFp_CONST_INTVECT(iv)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_POINTPLANARSHOCKINIT");
 FORTRAN_NAME( POINTPLANARSHOCKINIT ,pointplanarshockinit )(
      CHFt_REAL(density)
      ,CHFt_REALVECT(momentum)
      ,CHFt_REAL(energy)
      ,CHFt_VR(specDense)
      ,CHFt_CONST_INTVECT(iv)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDPOINTPLANARSHOCKINIT 

#ifndef GUARDSINEWAVEINIT 
#define GUARDSINEWAVEINIT 
// Prototype for Fortran procedure sinewaveinit ...
//
void FORTRAN_NAME( SINEWAVEINIT ,sinewaveinit )(
      CHFp_FRA(u)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) );

#define FORT_SINEWAVEINIT FORTRAN_NAME( inlineSINEWAVEINIT, inlineSINEWAVEINIT)
#define FORTNT_SINEWAVEINIT FORTRAN_NAME( SINEWAVEINIT, sinewaveinit)

inline void FORTRAN_NAME(inlineSINEWAVEINIT, inlineSINEWAVEINIT)(
      CHFp_FRA(u)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SINEWAVEINIT");
 FORTRAN_NAME( SINEWAVEINIT ,sinewaveinit )(
      CHFt_FRA(u)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSINEWAVEINIT 

#ifndef GUARDPOINTSINEWAVEINIT 
#define GUARDPOINTSINEWAVEINIT 
// Prototype for Fortran procedure pointsinewaveinit ...
//
void FORTRAN_NAME( POINTSINEWAVEINIT ,pointsinewaveinit )(
      CHFp_REAL(density)
      ,CHFp_REALVECT(momentum)
      ,CHFp_REAL(energy)
      ,CHFp_VR(specDense)
      ,CHFp_CONST_INTVECT(iv)
      ,CHFp_CONST_REAL(dx) );

#define FORT_POINTSINEWAVEINIT FORTRAN_NAME( inlinePOINTSINEWAVEINIT, inlinePOINTSINEWAVEINIT)
#define FORTNT_POINTSINEWAVEINIT FORTRAN_NAME( POINTSINEWAVEINIT, pointsinewaveinit)

inline void FORTRAN_NAME(inlinePOINTSINEWAVEINIT, inlinePOINTSINEWAVEINIT)(
      CHFp_REAL(density)
      ,CHFp_REALVECT(momentum)
      ,CHFp_REAL(energy)
      ,CHFp_VR(specDense)
      ,CHFp_CONST_INTVECT(iv)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_POINTSINEWAVEINIT");
 FORTRAN_NAME( POINTSINEWAVEINIT ,pointsinewaveinit )(
      CHFt_REAL(density)
      ,CHFt_REALVECT(momentum)
      ,CHFt_REAL(energy)
      ,CHFt_VR(specDense)
      ,CHFt_CONST_INTVECT(iv)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDPOINTSINEWAVEINIT 

#ifndef GUARDSLOPEBCS 
#define GUARDSLOPEBCS 
// Prototype for Fortran procedure SLOPEBCS ...
//
void FORTRAN_NAME( SLOPEBCS ,slopebcs )(
      CHFp_FRA(dW)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi) );

#define FORT_SLOPEBCS FORTRAN_NAME( inlineSLOPEBCS, inlineSLOPEBCS)
#define FORTNT_SLOPEBCS FORTRAN_NAME( SLOPEBCS, slopebcs)

inline void FORTRAN_NAME(inlineSLOPEBCS, inlineSLOPEBCS)(
      CHFp_FRA(dW)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi) )
{
 CH_TIMELEAF("FORT_SLOPEBCS");
 FORTRAN_NAME( SLOPEBCS ,slopebcs )(
      CHFt_FRA(dW)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(dir)
      ,CHFt_BOX(loBox)
      ,CHFt_CONST_INT(hasLo)
      ,CHFt_BOX(hiBox)
      ,CHFt_CONST_INT(hasHi) );
}
#endif  // GUARDSLOPEBCS 

#ifndef GUARDEXTRAPBC 
#define GUARDEXTRAPBC 
// Prototype for Fortran procedure extrapbc ...
//
void FORTRAN_NAME( EXTRAPBC ,extrapbc )(
      CHFp_FRA(flux)
      ,CHFp_CONST_FRA(wextrap)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_INT(nspec)
      ,CHFp_BOX(box) );

#define FORT_EXTRAPBC FORTRAN_NAME( inlineEXTRAPBC, inlineEXTRAPBC)
#define FORTNT_EXTRAPBC FORTRAN_NAME( EXTRAPBC, extrapbc)

inline void FORTRAN_NAME(inlineEXTRAPBC, inlineEXTRAPBC)(
      CHFp_FRA(flux)
      ,CHFp_CONST_FRA(wextrap)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_INT(nspec)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_EXTRAPBC");
 FORTRAN_NAME( EXTRAPBC ,extrapbc )(
      CHFt_FRA(flux)
      ,CHFt_CONST_FRA(wextrap)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_CONST_INT(nspec)
      ,CHFt_BOX(box) );
}
#endif  // GUARDEXTRAPBC 

#ifndef GUARDSLIPWALLSOLIDBC 
#define GUARDSLIPWALLSOLIDBC 
// Prototype for Fortran procedure SLIPWALLSOLIDBC ...
//
void FORTRAN_NAME( SLIPWALLSOLIDBC ,slipwallsolidbc )(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(Wextrap)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(box) );

#define FORT_SLIPWALLSOLIDBC FORTRAN_NAME( inlineSLIPWALLSOLIDBC, inlineSLIPWALLSOLIDBC)
#define FORTNT_SLIPWALLSOLIDBC FORTRAN_NAME( SLIPWALLSOLIDBC, slipwallsolidbc)

inline void FORTRAN_NAME(inlineSLIPWALLSOLIDBC, inlineSLIPWALLSOLIDBC)(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(Wextrap)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SLIPWALLSOLIDBC");
 FORTRAN_NAME( SLIPWALLSOLIDBC ,slipwallsolidbc )(
      CHFt_FRA(F)
      ,CHFt_CONST_FRA(Wextrap)
      ,CHFt_CONST_INT(lohisign)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSLIPWALLSOLIDBC 

}

#endif
