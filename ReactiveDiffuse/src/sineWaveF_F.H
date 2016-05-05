#ifndef _SINEWAVEF_F_H_
#define _SINEWAVEF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSETSINEWAVE 
#define GUARDSETSINEWAVE 
// Prototype for Fortran procedure setsinewave ...
//
void FORTRAN_NAME( SETSINEWAVE ,setsinewave )(
      CHFp_CONST_REALVECT(aucenter)
      ,CHFp_CONST_REAL(aurnot)
      ,CHFp_CONST_REAL(aupinf)
      ,CHFp_CONST_REAL(autinf)
      ,CHFp_CONST_VR(auvinf)
      ,CHFp_CONST_VR(aumassfrac)
      ,CHFp_CONST_REAL(aumach)
      ,CHFp_CONST_REAL(ausmallness) );

#define FORT_SETSINEWAVE FORTRAN_NAME( inlineSETSINEWAVE, inlineSETSINEWAVE)
#define FORTNT_SETSINEWAVE FORTRAN_NAME( SETSINEWAVE, setsinewave)

inline void FORTRAN_NAME(inlineSETSINEWAVE, inlineSETSINEWAVE)(
      CHFp_CONST_REALVECT(aucenter)
      ,CHFp_CONST_REAL(aurnot)
      ,CHFp_CONST_REAL(aupinf)
      ,CHFp_CONST_REAL(autinf)
      ,CHFp_CONST_VR(auvinf)
      ,CHFp_CONST_VR(aumassfrac)
      ,CHFp_CONST_REAL(aumach)
      ,CHFp_CONST_REAL(ausmallness) )
{
 CH_TIMELEAF("FORT_SETSINEWAVE");
 FORTRAN_NAME( SETSINEWAVE ,setsinewave )(
      CHFt_CONST_REALVECT(aucenter)
      ,CHFt_CONST_REAL(aurnot)
      ,CHFt_CONST_REAL(aupinf)
      ,CHFt_CONST_REAL(autinf)
      ,CHFt_CONST_VR(auvinf)
      ,CHFt_CONST_VR(aumassfrac)
      ,CHFt_CONST_REAL(aumach)
      ,CHFt_CONST_REAL(ausmallness) );
}
#endif  // GUARDSETSINEWAVE 

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

}

#endif
