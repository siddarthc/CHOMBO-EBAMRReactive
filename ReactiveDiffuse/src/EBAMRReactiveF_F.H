#ifndef _EBAMRREACTIVEF_F_H_
#define _EBAMRREACTIVEF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDGETRELATIVEGRAD 
#define GUARDGETRELATIVEGRAD 
// Prototype for Fortran procedure getrelativegrad ...
//
void FORTRAN_NAME( GETRELATIVEGRAD ,getrelativegrad )(
      CHFp_FRA1(du)
      ,CHFp_CONST_FRA1(u)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(lobox)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_BOX(hibox)
      ,CHFp_CONST_INT(hashi)
      ,CHFp_BOX(centerbox) );

#define FORT_GETRELATIVEGRAD FORTRAN_NAME( inlineGETRELATIVEGRAD, inlineGETRELATIVEGRAD)
#define FORTNT_GETRELATIVEGRAD FORTRAN_NAME( GETRELATIVEGRAD, getrelativegrad)

inline void FORTRAN_NAME(inlineGETRELATIVEGRAD, inlineGETRELATIVEGRAD)(
      CHFp_FRA1(du)
      ,CHFp_CONST_FRA1(u)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(lobox)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_BOX(hibox)
      ,CHFp_CONST_INT(hashi)
      ,CHFp_BOX(centerbox) )
{
 CH_TIMELEAF("FORT_GETRELATIVEGRAD");
 FORTRAN_NAME( GETRELATIVEGRAD ,getrelativegrad )(
      CHFt_FRA1(du)
      ,CHFt_CONST_FRA1(u)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(lobox)
      ,CHFt_CONST_INT(haslo)
      ,CHFt_BOX(hibox)
      ,CHFt_CONST_INT(hashi)
      ,CHFt_BOX(centerbox) );
}
#endif  // GUARDGETRELATIVEGRAD 

#ifndef GUARDMAGNITUDE 
#define GUARDMAGNITUDE 
// Prototype for Fortran procedure magnitude ...
//
void FORTRAN_NAME( MAGNITUDE ,magnitude )(
      CHFp_FRA1(magdata)
      ,CHFp_CONST_FRA(data)
      ,CHFp_BOX(box) );

#define FORT_MAGNITUDE FORTRAN_NAME( inlineMAGNITUDE, inlineMAGNITUDE)
#define FORTNT_MAGNITUDE FORTRAN_NAME( MAGNITUDE, magnitude)

inline void FORTRAN_NAME(inlineMAGNITUDE, inlineMAGNITUDE)(
      CHFp_FRA1(magdata)
      ,CHFp_CONST_FRA(data)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_MAGNITUDE");
 FORTRAN_NAME( MAGNITUDE ,magnitude )(
      CHFt_FRA1(magdata)
      ,CHFt_CONST_FRA(data)
      ,CHFt_BOX(box) );
}
#endif  // GUARDMAGNITUDE 

}

#endif
