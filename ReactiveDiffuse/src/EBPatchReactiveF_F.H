#ifndef _EBPATCHREACTIVEF_F_H_
#define _EBPATCHREACTIVEF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDINITIALIZE_CHEMISTRY 
#define GUARDINITIALIZE_CHEMISTRY 
// Prototype for Fortran procedure INITIALIZE_CHEMISTRY ...
//
void FORTRAN_NAME( INITIALIZE_CHEMISTRY ,initialize_chemistry )( );

#define FORT_INITIALIZE_CHEMISTRY FORTRAN_NAME( inlineINITIALIZE_CHEMISTRY, inlineINITIALIZE_CHEMISTRY)
#define FORTNT_INITIALIZE_CHEMISTRY FORTRAN_NAME( INITIALIZE_CHEMISTRY, initialize_chemistry)

inline void FORTRAN_NAME(inlineINITIALIZE_CHEMISTRY, inlineINITIALIZE_CHEMISTRY)( )
{
 CH_TIMELEAF("FORT_INITIALIZE_CHEMISTRY");
 FORTRAN_NAME( INITIALIZE_CHEMISTRY ,initialize_chemistry )( );
}
#endif  // GUARDINITIALIZE_CHEMISTRY 

#ifndef GUARDGET_NSPECIES 
#define GUARDGET_NSPECIES 
// Prototype for Fortran procedure get_nspecies ...
//
void FORTRAN_NAME( GET_NSPECIES ,get_nspecies )(
      CHFp_INT(spec) );

#define FORT_GET_NSPECIES FORTRAN_NAME( inlineGET_NSPECIES, inlineGET_NSPECIES)
#define FORTNT_GET_NSPECIES FORTRAN_NAME( GET_NSPECIES, get_nspecies)

inline void FORTRAN_NAME(inlineGET_NSPECIES, inlineGET_NSPECIES)(
      CHFp_INT(spec) )
{
 CH_TIMELEAF("FORT_GET_NSPECIES");
 FORTRAN_NAME( GET_NSPECIES ,get_nspecies )(
      CHFt_INT(spec) );
}
#endif  // GUARDGET_NSPECIES 

#ifndef GUARDGETSMALL 
#define GUARDGETSMALL 
// Prototype for Fortran procedure getsmall ...
//
void FORTRAN_NAME( GETSMALL ,getsmall )(
      CHFp_REAL(ausmall)
      ,CHFp_REAL(ausmallp)
      ,CHFp_REAL(ausmallu)
      ,CHFp_REAL(ausmallr) );

#define FORT_GETSMALL FORTRAN_NAME( inlineGETSMALL, inlineGETSMALL)
#define FORTNT_GETSMALL FORTRAN_NAME( GETSMALL, getsmall)

inline void FORTRAN_NAME(inlineGETSMALL, inlineGETSMALL)(
      CHFp_REAL(ausmall)
      ,CHFp_REAL(ausmallp)
      ,CHFp_REAL(ausmallu)
      ,CHFp_REAL(ausmallr) )
{
 CH_TIMELEAF("FORT_GETSMALL");
 FORTRAN_NAME( GETSMALL ,getsmall )(
      CHFt_REAL(ausmall)
      ,CHFt_REAL(ausmallp)
      ,CHFt_REAL(ausmallu)
      ,CHFt_REAL(ausmallr) );
}
#endif  // GUARDGETSMALL 

#ifndef GUARDFLOORPRIM 
#define GUARDFLOORPRIM 
// Prototype for Fortran procedure floorprim ...
//
void FORTRAN_NAME( FLOORPRIM ,floorprim )(
      CHFp_BOX(dcalc)
      ,CHFp_FRA(state) );

#define FORT_FLOORPRIM FORTRAN_NAME( inlineFLOORPRIM, inlineFLOORPRIM)
#define FORTNT_FLOORPRIM FORTRAN_NAME( FLOORPRIM, floorprim)

inline void FORTRAN_NAME(inlineFLOORPRIM, inlineFLOORPRIM)(
      CHFp_BOX(dcalc)
      ,CHFp_FRA(state) )
{
 CH_TIMELEAF("FORT_FLOORPRIM");
 FORTRAN_NAME( FLOORPRIM ,floorprim )(
      CHFt_BOX(dcalc)
      ,CHFt_FRA(state) );
}
#endif  // GUARDFLOORPRIM 

#ifndef GUARDFLOORCONS 
#define GUARDFLOORCONS 
// Prototype for Fortran procedure floorcons ...
//
void FORTRAN_NAME( FLOORCONS ,floorcons )(
      CHFp_BOX(dcalc)
      ,CHFp_FRA(state) );

#define FORT_FLOORCONS FORTRAN_NAME( inlineFLOORCONS, inlineFLOORCONS)
#define FORTNT_FLOORCONS FORTRAN_NAME( FLOORCONS, floorcons)

inline void FORTRAN_NAME(inlineFLOORCONS, inlineFLOORCONS)(
      CHFp_BOX(dcalc)
      ,CHFp_FRA(state) )
{
 CH_TIMELEAF("FORT_FLOORCONS");
 FORTRAN_NAME( FLOORCONS ,floorcons )(
      CHFt_BOX(dcalc)
      ,CHFt_FRA(state) );
}
#endif  // GUARDFLOORCONS 

#ifndef GUARDUPDATE 
#define GUARDUPDATE 
// Prototype for Fortran procedure update ...
//
void FORTRAN_NAME( UPDATE ,update )(
      CHFp_BOX(dcalc)
      ,CHFp_FRA(state)
      ,CHFp_CONST_FRA(flux)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_INT(nconserved)
      ,CHFp_CONST_REAL(dtbydx) );

#define FORT_UPDATE FORTRAN_NAME( inlineUPDATE, inlineUPDATE)
#define FORTNT_UPDATE FORTRAN_NAME( UPDATE, update)

inline void FORTRAN_NAME(inlineUPDATE, inlineUPDATE)(
      CHFp_BOX(dcalc)
      ,CHFp_FRA(state)
      ,CHFp_CONST_FRA(flux)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_INT(nconserved)
      ,CHFp_CONST_REAL(dtbydx) )
{
 CH_TIMELEAF("FORT_UPDATE");
 FORTRAN_NAME( UPDATE ,update )(
      CHFt_BOX(dcalc)
      ,CHFt_FRA(state)
      ,CHFt_CONST_FRA(flux)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_CONST_INT(nconserved)
      ,CHFt_CONST_REAL(dtbydx) );
}
#endif  // GUARDUPDATE 

#ifndef GUARDREGUPDATE 
#define GUARDREGUPDATE 
// Prototype for Fortran procedure regupdate ...
//
void FORTRAN_NAME( REGUPDATE ,regupdate )(
      CHFp_BOX(dcalc)
      ,CHFp_FRA(consstate)
      ,CHFp_CONST_FRA(divf)
      ,CHFp_CONST_INT(nconserved)
      ,CHFp_CONST_REAL(dt) );

#define FORT_REGUPDATE FORTRAN_NAME( inlineREGUPDATE, inlineREGUPDATE)
#define FORTNT_REGUPDATE FORTRAN_NAME( REGUPDATE, regupdate)

inline void FORTRAN_NAME(inlineREGUPDATE, inlineREGUPDATE)(
      CHFp_BOX(dcalc)
      ,CHFp_FRA(consstate)
      ,CHFp_CONST_FRA(divf)
      ,CHFp_CONST_INT(nconserved)
      ,CHFp_CONST_REAL(dt) )
{
 CH_TIMELEAF("FORT_REGUPDATE");
 FORTRAN_NAME( REGUPDATE ,regupdate )(
      CHFt_BOX(dcalc)
      ,CHFt_FRA(consstate)
      ,CHFt_CONST_FRA(divf)
      ,CHFt_CONST_INT(nconserved)
      ,CHFt_CONST_REAL(dt) );
}
#endif  // GUARDREGUPDATE 

#ifndef GUARDGETGRAD 
#define GUARDGETGRAD 
// Prototype for Fortran procedure getgrad ...
//
void FORTRAN_NAME( GETGRAD ,getgrad )(
      CHFp_FRA1(du)
      ,CHFp_CONST_FRA1(u)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(lobox)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_BOX(hibox)
      ,CHFp_CONST_INT(hashi)
      ,CHFp_BOX(centerbox) );

#define FORT_GETGRAD FORTRAN_NAME( inlineGETGRAD, inlineGETGRAD)
#define FORTNT_GETGRAD FORTRAN_NAME( GETGRAD, getgrad)

inline void FORTRAN_NAME(inlineGETGRAD, inlineGETGRAD)(
      CHFp_FRA1(du)
      ,CHFp_CONST_FRA1(u)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(lobox)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_BOX(hibox)
      ,CHFp_CONST_INT(hashi)
      ,CHFp_BOX(centerbox) )
{
 CH_TIMELEAF("FORT_GETGRAD");
 FORTRAN_NAME( GETGRAD ,getgrad )(
      CHFt_FRA1(du)
      ,CHFt_CONST_FRA1(u)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(lobox)
      ,CHFt_CONST_INT(haslo)
      ,CHFt_BOX(hibox)
      ,CHFt_CONST_INT(hashi)
      ,CHFt_BOX(centerbox) );
}
#endif  // GUARDGETGRAD 

#ifndef GUARDGETDPTWO 
#define GUARDGETDPTWO 
// Prototype for Fortran procedure getdptwo ...
//
void FORTRAN_NAME( GETDPTWO ,getdptwo )(
      CHFp_FRA1(delta2p)
      ,CHFp_CONST_FRA1(delta1p)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(lobox)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_BOX(hibox)
      ,CHFp_CONST_INT(hashi)
      ,CHFp_BOX(centerbox) );

#define FORT_GETDPTWO FORTRAN_NAME( inlineGETDPTWO, inlineGETDPTWO)
#define FORTNT_GETDPTWO FORTRAN_NAME( GETDPTWO, getdptwo)

inline void FORTRAN_NAME(inlineGETDPTWO, inlineGETDPTWO)(
      CHFp_FRA1(delta2p)
      ,CHFp_CONST_FRA1(delta1p)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(lobox)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_BOX(hibox)
      ,CHFp_CONST_INT(hashi)
      ,CHFp_BOX(centerbox) )
{
 CH_TIMELEAF("FORT_GETDPTWO");
 FORTRAN_NAME( GETDPTWO ,getdptwo )(
      CHFt_FRA1(delta2p)
      ,CHFt_CONST_FRA1(delta1p)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(lobox)
      ,CHFt_CONST_INT(haslo)
      ,CHFt_BOX(hibox)
      ,CHFt_CONST_INT(hashi)
      ,CHFt_BOX(centerbox) );
}
#endif  // GUARDGETDPTWO 

#ifndef GUARDGETFLAT 
#define GUARDGETFLAT 
// Prototype for Fortran procedure getflat ...
//
void FORTRAN_NAME( GETFLAT ,getflat )(
      CHFp_FRA1(zetatwiddle)
      ,CHFp_CONST_FRA1(delta1p)
      ,CHFp_CONST_FRA1(delta2p)
      ,CHFp_CONST_FRA1(bulkmin)
      ,CHFp_BOX(box) );

#define FORT_GETFLAT FORTRAN_NAME( inlineGETFLAT, inlineGETFLAT)
#define FORTNT_GETFLAT FORTRAN_NAME( GETFLAT, getflat)

inline void FORTRAN_NAME(inlineGETFLAT, inlineGETFLAT)(
      CHFp_FRA1(zetatwiddle)
      ,CHFp_CONST_FRA1(delta1p)
      ,CHFp_CONST_FRA1(delta2p)
      ,CHFp_CONST_FRA1(bulkmin)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_GETFLAT");
 FORTRAN_NAME( GETFLAT ,getflat )(
      CHFt_FRA1(zetatwiddle)
      ,CHFt_CONST_FRA1(delta1p)
      ,CHFt_CONST_FRA1(delta2p)
      ,CHFt_CONST_FRA1(bulkmin)
      ,CHFt_BOX(box) );
}
#endif  // GUARDGETFLAT 

#ifndef GUARDMIN3PTS 
#define GUARDMIN3PTS 
// Prototype for Fortran procedure min3pts ...
//
void FORTRAN_NAME( MIN3PTS ,min3pts )(
      CHFp_FRA1(mindata)
      ,CHFp_CONST_FRA1(data)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(lobox)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_BOX(hibox)
      ,CHFp_CONST_INT(hashi)
      ,CHFp_BOX(centerbox) );

#define FORT_MIN3PTS FORTRAN_NAME( inlineMIN3PTS, inlineMIN3PTS)
#define FORTNT_MIN3PTS FORTRAN_NAME( MIN3PTS, min3pts)

inline void FORTRAN_NAME(inlineMIN3PTS, inlineMIN3PTS)(
      CHFp_FRA1(mindata)
      ,CHFp_CONST_FRA1(data)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(lobox)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_BOX(hibox)
      ,CHFp_CONST_INT(hashi)
      ,CHFp_BOX(centerbox) )
{
 CH_TIMELEAF("FORT_MIN3PTS");
 FORTRAN_NAME( MIN3PTS ,min3pts )(
      CHFt_FRA1(mindata)
      ,CHFt_CONST_FRA1(data)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(lobox)
      ,CHFt_CONST_INT(haslo)
      ,CHFt_BOX(hibox)
      ,CHFt_CONST_INT(hashi)
      ,CHFt_BOX(centerbox) );
}
#endif  // GUARDMIN3PTS 

#ifndef GUARDMINFLAT 
#define GUARDMINFLAT 
// Prototype for Fortran procedure minflat ...
//
void FORTRAN_NAME( MINFLAT ,minflat )(
      CHFp_FRA1(flattening)
      ,CHFp_CONST_FRA(zetadir)
      ,CHFp_CONST_FRA(du)
      ,CHFp_BOX(box) );

#define FORT_MINFLAT FORTRAN_NAME( inlineMINFLAT, inlineMINFLAT)
#define FORTNT_MINFLAT FORTRAN_NAME( MINFLAT, minflat)

inline void FORTRAN_NAME(inlineMINFLAT, inlineMINFLAT)(
      CHFp_FRA1(flattening)
      ,CHFp_CONST_FRA(zetadir)
      ,CHFp_CONST_FRA(du)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_MINFLAT");
 FORTRAN_NAME( MINFLAT ,minflat )(
      CHFt_FRA1(flattening)
      ,CHFt_CONST_FRA(zetadir)
      ,CHFt_CONST_FRA(du)
      ,CHFt_BOX(box) );
}
#endif  // GUARDMINFLAT 

#ifndef GUARDDIVUONED 
#define GUARDDIVUONED 
// Prototype for Fortran procedure divuoned ...
//
void FORTRAN_NAME( DIVUONED ,divuoned )(
      CHFp_FRA1(divu)
      ,CHFp_CONST_FRA1(velnorm)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(centerbox) );

#define FORT_DIVUONED FORTRAN_NAME( inlineDIVUONED, inlineDIVUONED)
#define FORTNT_DIVUONED FORTRAN_NAME( DIVUONED, divuoned)

inline void FORTRAN_NAME(inlineDIVUONED, inlineDIVUONED)(
      CHFp_FRA1(divu)
      ,CHFp_CONST_FRA1(velnorm)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(centerbox) )
{
 CH_TIMELEAF("FORT_DIVUONED");
 FORTRAN_NAME( DIVUONED ,divuoned )(
      CHFt_FRA1(divu)
      ,CHFt_CONST_FRA1(velnorm)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_BOX(centerbox) );
}
#endif  // GUARDDIVUONED 

#ifndef GUARDDIVUTRAN 
#define GUARDDIVUTRAN 
// Prototype for Fortran procedure divutran ...
//
void FORTRAN_NAME( DIVUTRAN ,divutran )(
      CHFp_FRA1(divu)
      ,CHFp_CONST_FRA1(slopevel)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(centerbox) );

#define FORT_DIVUTRAN FORTRAN_NAME( inlineDIVUTRAN, inlineDIVUTRAN)
#define FORTNT_DIVUTRAN FORTRAN_NAME( DIVUTRAN, divutran)

inline void FORTRAN_NAME(inlineDIVUTRAN, inlineDIVUTRAN)(
      CHFp_FRA1(divu)
      ,CHFp_CONST_FRA1(slopevel)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(centerbox) )
{
 CH_TIMELEAF("FORT_DIVUTRAN");
 FORTRAN_NAME( DIVUTRAN ,divutran )(
      CHFt_FRA1(divu)
      ,CHFt_CONST_FRA1(slopevel)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_BOX(centerbox) );
}
#endif  // GUARDDIVUTRAN 

#ifndef GUARDDIVUEDGE 
#define GUARDDIVUEDGE 
// Prototype for Fortran procedure divuedge ...
//
void FORTRAN_NAME( DIVUEDGE ,divuedge )(
      CHFp_FRA1(divu)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(lobox)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_BOX(hibox)
      ,CHFp_CONST_INT(hashi) );

#define FORT_DIVUEDGE FORTRAN_NAME( inlineDIVUEDGE, inlineDIVUEDGE)
#define FORTNT_DIVUEDGE FORTRAN_NAME( DIVUEDGE, divuedge)

inline void FORTRAN_NAME(inlineDIVUEDGE, inlineDIVUEDGE)(
      CHFp_FRA1(divu)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(lobox)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_BOX(hibox)
      ,CHFp_CONST_INT(hashi) )
{
 CH_TIMELEAF("FORT_DIVUEDGE");
 FORTRAN_NAME( DIVUEDGE ,divuedge )(
      CHFt_FRA1(divu)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_BOX(lobox)
      ,CHFt_CONST_INT(haslo)
      ,CHFt_BOX(hibox)
      ,CHFt_CONST_INT(hashi) );
}
#endif  // GUARDDIVUEDGE 

#ifndef GUARDARTVISC 
#define GUARDARTVISC 
// Prototype for Fortran procedure artvisc ...
//
void FORTRAN_NAME( ARTVISC ,artvisc )(
      CHFp_FRA(f)
      ,CHFp_CONST_FRA(u)
      ,CHFp_CONST_FRA1(divu)
      ,CHFp_CONST_REAL(coeff)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(numcons)
      ,CHFp_CONST_REAL(dx) );

#define FORT_ARTVISC FORTRAN_NAME( inlineARTVISC, inlineARTVISC)
#define FORTNT_ARTVISC FORTRAN_NAME( ARTVISC, artvisc)

inline void FORTRAN_NAME(inlineARTVISC, inlineARTVISC)(
      CHFp_FRA(f)
      ,CHFp_CONST_FRA(u)
      ,CHFp_CONST_FRA1(divu)
      ,CHFp_CONST_REAL(coeff)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(numcons)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_ARTVISC");
 FORTRAN_NAME( ARTVISC ,artvisc )(
      CHFt_FRA(f)
      ,CHFt_CONST_FRA(u)
      ,CHFt_CONST_FRA1(divu)
      ,CHFt_CONST_REAL(coeff)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box)
      ,CHFt_CONST_INT(numcons)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDARTVISC 

#ifndef GUARDDIVERGEF 
#define GUARDDIVERGEF 
// Prototype for Fortran procedure divergef ...
//
void FORTRAN_NAME( DIVERGEF ,divergef )(
      CHFp_BOX(dcalc)
      ,CHFp_FRA(divf)
      ,CHFp_CONST_FRA(flux)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_INT(nconserved)
      ,CHFp_CONST_REAL(dx) );

#define FORT_DIVERGEF FORTRAN_NAME( inlineDIVERGEF, inlineDIVERGEF)
#define FORTNT_DIVERGEF FORTRAN_NAME( DIVERGEF, divergef)

inline void FORTRAN_NAME(inlineDIVERGEF, inlineDIVERGEF)(
      CHFp_BOX(dcalc)
      ,CHFp_FRA(divf)
      ,CHFp_CONST_FRA(flux)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_INT(nconserved)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_DIVERGEF");
 FORTRAN_NAME( DIVERGEF ,divergef )(
      CHFt_BOX(dcalc)
      ,CHFt_FRA(divf)
      ,CHFt_CONST_FRA(flux)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_CONST_INT(nconserved)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDDIVERGEF 

#ifndef GUARDGETSOUNDSPEED 
#define GUARDGETSOUNDSPEED 
// Prototype for Fortran procedure getsoundspeed ...
//
void FORTRAN_NAME( GETSOUNDSPEED ,getsoundspeed )(
      CHFp_CONST_REAL(intenrg)
      ,CHFp_CONST_VR(massfrac)
      ,CHFp_REAL(soundspeed) );

#define FORT_GETSOUNDSPEED FORTRAN_NAME( inlineGETSOUNDSPEED, inlineGETSOUNDSPEED)
#define FORTNT_GETSOUNDSPEED FORTRAN_NAME( GETSOUNDSPEED, getsoundspeed)

inline void FORTRAN_NAME(inlineGETSOUNDSPEED, inlineGETSOUNDSPEED)(
      CHFp_CONST_REAL(intenrg)
      ,CHFp_CONST_VR(massfrac)
      ,CHFp_REAL(soundspeed) )
{
 CH_TIMELEAF("FORT_GETSOUNDSPEED");
 FORTRAN_NAME( GETSOUNDSPEED ,getsoundspeed )(
      CHFt_CONST_REAL(intenrg)
      ,CHFt_CONST_VR(massfrac)
      ,CHFt_REAL(soundspeed) );
}
#endif  // GUARDGETSOUNDSPEED 

#ifndef GUARDGETSOUNDSPEEDTEMP 
#define GUARDGETSOUNDSPEEDTEMP 
// Prototype for Fortran procedure getsoundspeedtemp ...
//
void FORTRAN_NAME( GETSOUNDSPEEDTEMP ,getsoundspeedtemp )(
      CHFp_CONST_REAL(temp)
      ,CHFp_CONST_VR(massfrac)
      ,CHFp_REAL(soundspeed) );

#define FORT_GETSOUNDSPEEDTEMP FORTRAN_NAME( inlineGETSOUNDSPEEDTEMP, inlineGETSOUNDSPEEDTEMP)
#define FORTNT_GETSOUNDSPEEDTEMP FORTRAN_NAME( GETSOUNDSPEEDTEMP, getsoundspeedtemp)

inline void FORTRAN_NAME(inlineGETSOUNDSPEEDTEMP, inlineGETSOUNDSPEEDTEMP)(
      CHFp_CONST_REAL(temp)
      ,CHFp_CONST_VR(massfrac)
      ,CHFp_REAL(soundspeed) )
{
 CH_TIMELEAF("FORT_GETSOUNDSPEEDTEMP");
 FORTRAN_NAME( GETSOUNDSPEEDTEMP ,getsoundspeedtemp )(
      CHFt_CONST_REAL(temp)
      ,CHFt_CONST_VR(massfrac)
      ,CHFt_REAL(soundspeed) );
}
#endif  // GUARDGETSOUNDSPEEDTEMP 

#ifndef GUARDCONS2PRM 
#define GUARDCONS2PRM 
// Prototype for Fortran procedure cons2prm ...
//
void FORTRAN_NAME( CONS2PRM ,cons2prm )(
      CHFp_BOX(dcalc)
      ,CHFp_CONST_FRA(u)
      ,CHFp_FRA(q)
      ,CHFp_CONST_INT(logflag)
      ,CHFp_CONST_INT(iverbose) );

#define FORT_CONS2PRM FORTRAN_NAME( inlineCONS2PRM, inlineCONS2PRM)
#define FORTNT_CONS2PRM FORTRAN_NAME( CONS2PRM, cons2prm)

inline void FORTRAN_NAME(inlineCONS2PRM, inlineCONS2PRM)(
      CHFp_BOX(dcalc)
      ,CHFp_CONST_FRA(u)
      ,CHFp_FRA(q)
      ,CHFp_CONST_INT(logflag)
      ,CHFp_CONST_INT(iverbose) )
{
 CH_TIMELEAF("FORT_CONS2PRM");
 FORTRAN_NAME( CONS2PRM ,cons2prm )(
      CHFt_BOX(dcalc)
      ,CHFt_CONST_FRA(u)
      ,CHFt_FRA(q)
      ,CHFt_CONST_INT(logflag)
      ,CHFt_CONST_INT(iverbose) );
}
#endif  // GUARDCONS2PRM 

#ifndef GUARDPOINTCONS2PRM 
#define GUARDPOINTCONS2PRM 
// Prototype for Fortran procedure pointcons2prm ...
//
void FORTRAN_NAME( POINTCONS2PRM ,pointcons2prm )(
      CHFp_VR(conserved)
      ,CHFp_VR(primitive)
      ,CHFp_INT(success)
      ,CHFp_CONST_INT(logflag) );

#define FORT_POINTCONS2PRM FORTRAN_NAME( inlinePOINTCONS2PRM, inlinePOINTCONS2PRM)
#define FORTNT_POINTCONS2PRM FORTRAN_NAME( POINTCONS2PRM, pointcons2prm)

inline void FORTRAN_NAME(inlinePOINTCONS2PRM, inlinePOINTCONS2PRM)(
      CHFp_VR(conserved)
      ,CHFp_VR(primitive)
      ,CHFp_INT(success)
      ,CHFp_CONST_INT(logflag) )
{
 CH_TIMELEAF("FORT_POINTCONS2PRM");
 FORTRAN_NAME( POINTCONS2PRM ,pointcons2prm )(
      CHFt_VR(conserved)
      ,CHFt_VR(primitive)
      ,CHFt_INT(success)
      ,CHFt_CONST_INT(logflag) );
}
#endif  // GUARDPOINTCONS2PRM 

#ifndef GUARDSECONDSLOPEDIFFS 
#define GUARDSECONDSLOPEDIFFS 
// Prototype for Fortran procedure secondslopediffs ...
//
void FORTRAN_NAME( SECONDSLOPEDIFFS ,secondslopediffs )(
      CHFp_FRA(deltawc)
      ,CHFp_FRA(deltawl)
      ,CHFp_FRA(deltawr)
      ,CHFp_CONST_FRA(w)
      ,CHFp_CONST_INT(numslopes)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(lobox)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_BOX(hibox)
      ,CHFp_CONST_INT(hashi)
      ,CHFp_BOX(centerbox) );

#define FORT_SECONDSLOPEDIFFS FORTRAN_NAME( inlineSECONDSLOPEDIFFS, inlineSECONDSLOPEDIFFS)
#define FORTNT_SECONDSLOPEDIFFS FORTRAN_NAME( SECONDSLOPEDIFFS, secondslopediffs)

inline void FORTRAN_NAME(inlineSECONDSLOPEDIFFS, inlineSECONDSLOPEDIFFS)(
      CHFp_FRA(deltawc)
      ,CHFp_FRA(deltawl)
      ,CHFp_FRA(deltawr)
      ,CHFp_CONST_FRA(w)
      ,CHFp_CONST_INT(numslopes)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(lobox)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_BOX(hibox)
      ,CHFp_CONST_INT(hashi)
      ,CHFp_BOX(centerbox) )
{
 CH_TIMELEAF("FORT_SECONDSLOPEDIFFS");
 FORTRAN_NAME( SECONDSLOPEDIFFS ,secondslopediffs )(
      CHFt_FRA(deltawc)
      ,CHFt_FRA(deltawl)
      ,CHFt_FRA(deltawr)
      ,CHFt_CONST_FRA(w)
      ,CHFt_CONST_INT(numslopes)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(lobox)
      ,CHFt_CONST_INT(haslo)
      ,CHFt_BOX(hibox)
      ,CHFt_CONST_INT(hashi)
      ,CHFt_BOX(centerbox) );
}
#endif  // GUARDSECONDSLOPEDIFFS 

#ifndef GUARDVLLIMITER 
#define GUARDVLLIMITER 
// Prototype for Fortran procedure vllimiter ...
//
void FORTRAN_NAME( VLLIMITER ,vllimiter )(
      CHFp_FRA(slopeprim)
      ,CHFp_CONST_FRA(slopeleft)
      ,CHFp_CONST_FRA(sloperigh)
      ,CHFp_BOX(dcalc) );

#define FORT_VLLIMITER FORTRAN_NAME( inlineVLLIMITER, inlineVLLIMITER)
#define FORTNT_VLLIMITER FORTRAN_NAME( VLLIMITER, vllimiter)

inline void FORTRAN_NAME(inlineVLLIMITER, inlineVLLIMITER)(
      CHFp_FRA(slopeprim)
      ,CHFp_CONST_FRA(slopeleft)
      ,CHFp_CONST_FRA(sloperigh)
      ,CHFp_BOX(dcalc) )
{
 CH_TIMELEAF("FORT_VLLIMITER");
 FORTRAN_NAME( VLLIMITER ,vllimiter )(
      CHFt_FRA(slopeprim)
      ,CHFt_CONST_FRA(slopeleft)
      ,CHFt_CONST_FRA(sloperigh)
      ,CHFt_BOX(dcalc) );
}
#endif  // GUARDVLLIMITER 

#ifndef GUARDPOINTVLLIMITER 
#define GUARDPOINTVLLIMITER 
// Prototype for Fortran procedure pointvllimiter ...
//
void FORTRAN_NAME( POINTVLLIMITER ,pointvllimiter )(
      CHFp_REAL(dqlim)
      ,CHFp_CONST_REAL(dql)
      ,CHFp_CONST_REAL(dqr) );

#define FORT_POINTVLLIMITER FORTRAN_NAME( inlinePOINTVLLIMITER, inlinePOINTVLLIMITER)
#define FORTNT_POINTVLLIMITER FORTRAN_NAME( POINTVLLIMITER, pointvllimiter)

inline void FORTRAN_NAME(inlinePOINTVLLIMITER, inlinePOINTVLLIMITER)(
      CHFp_REAL(dqlim)
      ,CHFp_CONST_REAL(dql)
      ,CHFp_CONST_REAL(dqr) )
{
 CH_TIMELEAF("FORT_POINTVLLIMITER");
 FORTRAN_NAME( POINTVLLIMITER ,pointvllimiter )(
      CHFt_REAL(dqlim)
      ,CHFt_CONST_REAL(dql)
      ,CHFt_CONST_REAL(dqr) );
}
#endif  // GUARDPOINTVLLIMITER 

#ifndef GUARDFORTHSLOPEDIFFS 
#define GUARDFORTHSLOPEDIFFS 
// Prototype for Fortran procedure forthslopediffs ...
//
void FORTRAN_NAME( FORTHSLOPEDIFFS ,forthslopediffs )(
      CHFp_FRA(delta4wc)
      ,CHFp_CONST_FRA(w)
      ,CHFp_CONST_FRA(delta2w)
      ,CHFp_CONST_INT(numslopes)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(lobox)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_BOX(hibox)
      ,CHFp_CONST_INT(hashi)
      ,CHFp_BOX(centerbox) );

#define FORT_FORTHSLOPEDIFFS FORTRAN_NAME( inlineFORTHSLOPEDIFFS, inlineFORTHSLOPEDIFFS)
#define FORTNT_FORTHSLOPEDIFFS FORTRAN_NAME( FORTHSLOPEDIFFS, forthslopediffs)

inline void FORTRAN_NAME(inlineFORTHSLOPEDIFFS, inlineFORTHSLOPEDIFFS)(
      CHFp_FRA(delta4wc)
      ,CHFp_CONST_FRA(w)
      ,CHFp_CONST_FRA(delta2w)
      ,CHFp_CONST_INT(numslopes)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(lobox)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_BOX(hibox)
      ,CHFp_CONST_INT(hashi)
      ,CHFp_BOX(centerbox) )
{
 CH_TIMELEAF("FORT_FORTHSLOPEDIFFS");
 FORTRAN_NAME( FORTHSLOPEDIFFS ,forthslopediffs )(
      CHFt_FRA(delta4wc)
      ,CHFt_CONST_FRA(w)
      ,CHFt_CONST_FRA(delta2w)
      ,CHFt_CONST_INT(numslopes)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(lobox)
      ,CHFt_CONST_INT(haslo)
      ,CHFt_BOX(hibox)
      ,CHFt_CONST_INT(hashi)
      ,CHFt_BOX(centerbox) );
}
#endif  // GUARDFORTHSLOPEDIFFS 

#ifndef GUARDAPPLYFLAT 
#define GUARDAPPLYFLAT 
// Prototype for Fortran procedure applyflat ...
//
void FORTRAN_NAME( APPLYFLAT ,applyflat )(
      CHFp_FRA(dw)
      ,CHFp_CONST_FRA1(flattening)
      ,CHFp_CONST_INT(numslopes)
      ,CHFp_BOX(box) );

#define FORT_APPLYFLAT FORTRAN_NAME( inlineAPPLYFLAT, inlineAPPLYFLAT)
#define FORTNT_APPLYFLAT FORTRAN_NAME( APPLYFLAT, applyflat)

inline void FORTRAN_NAME(inlineAPPLYFLAT, inlineAPPLYFLAT)(
      CHFp_FRA(dw)
      ,CHFp_CONST_FRA1(flattening)
      ,CHFp_CONST_INT(numslopes)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_APPLYFLAT");
 FORTRAN_NAME( APPLYFLAT ,applyflat )(
      CHFt_FRA(dw)
      ,CHFt_CONST_FRA1(flattening)
      ,CHFt_CONST_INT(numslopes)
      ,CHFt_BOX(box) );
}
#endif  // GUARDAPPLYFLAT 

#ifndef GUARDSETSOURCERZ 
#define GUARDSETSOURCERZ 
// Prototype for Fortran procedure setsourcerz ...
//
void FORTRAN_NAME( SETSOURCERZ ,setsourcerz )(
      CHFp_FRA(sourcerz)
      ,CHFp_CONST_FRA(state)
      ,CHFp_CONST_REAL(dr)
      ,CHFp_BOX(dcalc) );

#define FORT_SETSOURCERZ FORTRAN_NAME( inlineSETSOURCERZ, inlineSETSOURCERZ)
#define FORTNT_SETSOURCERZ FORTRAN_NAME( SETSOURCERZ, setsourcerz)

inline void FORTRAN_NAME(inlineSETSOURCERZ, inlineSETSOURCERZ)(
      CHFp_FRA(sourcerz)
      ,CHFp_CONST_FRA(state)
      ,CHFp_CONST_REAL(dr)
      ,CHFp_BOX(dcalc) )
{
 CH_TIMELEAF("FORT_SETSOURCERZ");
 FORTRAN_NAME( SETSOURCERZ ,setsourcerz )(
      CHFt_FRA(sourcerz)
      ,CHFt_CONST_FRA(state)
      ,CHFt_CONST_REAL(dr)
      ,CHFt_BOX(dcalc) );
}
#endif  // GUARDSETSOURCERZ 

#ifndef GUARDPOINTSETSOURCERZ 
#define GUARDPOINTSETSOURCERZ 
// Prototype for Fortran procedure pointsetsourcerz ...
//
void FORTRAN_NAME( POINTSETSOURCERZ ,pointsetsourcerz )(
      CHFp_REAL(densesource)
      ,CHFp_REAL(presssource)
      ,CHFp_CONST_REAL(dense)
      ,CHFp_CONST_REALVECT(momen)
      ,CHFp_CONST_REAL(energ)
      ,CHFp_CONST_REAL(radius) );

#define FORT_POINTSETSOURCERZ FORTRAN_NAME( inlinePOINTSETSOURCERZ, inlinePOINTSETSOURCERZ)
#define FORTNT_POINTSETSOURCERZ FORTRAN_NAME( POINTSETSOURCERZ, pointsetsourcerz)

inline void FORTRAN_NAME(inlinePOINTSETSOURCERZ, inlinePOINTSETSOURCERZ)(
      CHFp_REAL(densesource)
      ,CHFp_REAL(presssource)
      ,CHFp_CONST_REAL(dense)
      ,CHFp_CONST_REALVECT(momen)
      ,CHFp_CONST_REAL(energ)
      ,CHFp_CONST_REAL(radius) )
{
 CH_TIMELEAF("FORT_POINTSETSOURCERZ");
 FORTRAN_NAME( POINTSETSOURCERZ ,pointsetsourcerz )(
      CHFt_REAL(densesource)
      ,CHFt_REAL(presssource)
      ,CHFt_CONST_REAL(dense)
      ,CHFt_CONST_REALVECT(momen)
      ,CHFt_CONST_REAL(energ)
      ,CHFt_CONST_REAL(radius) );
}
#endif  // GUARDPOINTSETSOURCERZ 

#ifndef GUARDPRED 
#define GUARDPRED 
// Prototype for Fortran procedure pred ...
//
void FORTRAN_NAME( PRED ,pred )(
      CHFp_BOX(dcalc)
      ,CHFp_CONST_FRA(q)
      ,CHFp_CONST_FRA(dq)
      ,CHFp_FRA(qlo)
      ,CHFp_FRA(qhi)
      ,CHFp_CONST_INT(normdir)
      ,CHFp_CONST_REAL(dtbydx)
      ,CHFp_CONST_INT(useflat) );

#define FORT_PRED FORTRAN_NAME( inlinePRED, inlinePRED)
#define FORTNT_PRED FORTRAN_NAME( PRED, pred)

inline void FORTRAN_NAME(inlinePRED, inlinePRED)(
      CHFp_BOX(dcalc)
      ,CHFp_CONST_FRA(q)
      ,CHFp_CONST_FRA(dq)
      ,CHFp_FRA(qlo)
      ,CHFp_FRA(qhi)
      ,CHFp_CONST_INT(normdir)
      ,CHFp_CONST_REAL(dtbydx)
      ,CHFp_CONST_INT(useflat) )
{
 CH_TIMELEAF("FORT_PRED");
 FORTRAN_NAME( PRED ,pred )(
      CHFt_BOX(dcalc)
      ,CHFt_CONST_FRA(q)
      ,CHFt_CONST_FRA(dq)
      ,CHFt_FRA(qlo)
      ,CHFt_FRA(qhi)
      ,CHFt_CONST_INT(normdir)
      ,CHFt_CONST_REAL(dtbydx)
      ,CHFt_CONST_INT(useflat) );
}
#endif  // GUARDPRED 

#ifndef GUARDPOINTLIMIT 
#define GUARDPOINTLIMIT 
// Prototype for Fortran procedure pointlimit ...
//
void FORTRAN_NAME( POINTLIMIT ,pointlimit )(
      CHFp_REAL(dense)
      ,CHFp_REAL(denlo)
      ,CHFp_REAL(denhi)
      ,CHFp_REALVECT(veloc)
      ,CHFp_REALVECT(vello)
      ,CHFp_REALVECT(velhi)
      ,CHFp_REAL(press)
      ,CHFp_REAL(prelo)
      ,CHFp_REAL(prehi) );

#define FORT_POINTLIMIT FORTRAN_NAME( inlinePOINTLIMIT, inlinePOINTLIMIT)
#define FORTNT_POINTLIMIT FORTRAN_NAME( POINTLIMIT, pointlimit)

inline void FORTRAN_NAME(inlinePOINTLIMIT, inlinePOINTLIMIT)(
      CHFp_REAL(dense)
      ,CHFp_REAL(denlo)
      ,CHFp_REAL(denhi)
      ,CHFp_REALVECT(veloc)
      ,CHFp_REALVECT(vello)
      ,CHFp_REALVECT(velhi)
      ,CHFp_REAL(press)
      ,CHFp_REAL(prelo)
      ,CHFp_REAL(prehi) )
{
 CH_TIMELEAF("FORT_POINTLIMIT");
 FORTRAN_NAME( POINTLIMIT ,pointlimit )(
      CHFt_REAL(dense)
      ,CHFt_REAL(denlo)
      ,CHFt_REAL(denhi)
      ,CHFt_REALVECT(veloc)
      ,CHFt_REALVECT(vello)
      ,CHFt_REALVECT(velhi)
      ,CHFt_REAL(press)
      ,CHFt_REAL(prelo)
      ,CHFt_REAL(prehi) );
}
#endif  // GUARDPOINTLIMIT 

#ifndef GUARDPOINTPRED 
#define GUARDPOINTPRED 
// Prototype for Fortran procedure pointpred ...
//
void FORTRAN_NAME( POINTPRED ,pointpred )(
      CHFp_VR(primit)
      ,CHFp_VR(pslope)
      ,CHFp_VR(primlo)
      ,CHFp_VR(primhi)
      ,CHFp_CONST_INT(normdir)
      ,CHFp_CONST_REAL(dtbydx)
      ,CHFp_CONST_INT(useflat) );

#define FORT_POINTPRED FORTRAN_NAME( inlinePOINTPRED, inlinePOINTPRED)
#define FORTNT_POINTPRED FORTRAN_NAME( POINTPRED, pointpred)

inline void FORTRAN_NAME(inlinePOINTPRED, inlinePOINTPRED)(
      CHFp_VR(primit)
      ,CHFp_VR(pslope)
      ,CHFp_VR(primlo)
      ,CHFp_VR(primhi)
      ,CHFp_CONST_INT(normdir)
      ,CHFp_CONST_REAL(dtbydx)
      ,CHFp_CONST_INT(useflat) )
{
 CH_TIMELEAF("FORT_POINTPRED");
 FORTRAN_NAME( POINTPRED ,pointpred )(
      CHFt_VR(primit)
      ,CHFt_VR(pslope)
      ,CHFt_VR(primlo)
      ,CHFt_VR(primhi)
      ,CHFt_CONST_INT(normdir)
      ,CHFt_CONST_REAL(dtbydx)
      ,CHFt_CONST_INT(useflat) );
}
#endif  // GUARDPOINTPRED 

#ifndef GUARDINCSOURCE 
#define GUARDINCSOURCE 
// Prototype for Fortran procedure incsource ...
//
void FORTRAN_NAME( INCSOURCE ,incsource )(
      CHFp_FRA(prim)
      ,CHFp_CONST_FRA(source)
      ,CHFp_CONST_REAL(scale)
      ,CHFp_BOX(dcalc) );

#define FORT_INCSOURCE FORTRAN_NAME( inlineINCSOURCE, inlineINCSOURCE)
#define FORTNT_INCSOURCE FORTRAN_NAME( INCSOURCE, incsource)

inline void FORTRAN_NAME(inlineINCSOURCE, inlineINCSOURCE)(
      CHFp_FRA(prim)
      ,CHFp_CONST_FRA(source)
      ,CHFp_CONST_REAL(scale)
      ,CHFp_BOX(dcalc) )
{
 CH_TIMELEAF("FORT_INCSOURCE");
 FORTRAN_NAME( INCSOURCE ,incsource )(
      CHFt_FRA(prim)
      ,CHFt_CONST_FRA(source)
      ,CHFt_CONST_REAL(scale)
      ,CHFt_BOX(dcalc) );
}
#endif  // GUARDINCSOURCE 

#ifndef GUARDFLUXASSEMBLE 
#define GUARDFLUXASSEMBLE 
// Prototype for Fortran procedure fluxassemble ...
//
void FORTRAN_NAME( FLUXASSEMBLE ,fluxassemble )(
      CHFp_FRA(fluxregflux)
      ,CHFp_CONST_FRA(godunovflux)
      ,CHFp_CONST_REAL(dr)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(dcalc) );

#define FORT_FLUXASSEMBLE FORTRAN_NAME( inlineFLUXASSEMBLE, inlineFLUXASSEMBLE)
#define FORTNT_FLUXASSEMBLE FORTRAN_NAME( FLUXASSEMBLE, fluxassemble)

inline void FORTRAN_NAME(inlineFLUXASSEMBLE, inlineFLUXASSEMBLE)(
      CHFp_FRA(fluxregflux)
      ,CHFp_CONST_FRA(godunovflux)
      ,CHFp_CONST_REAL(dr)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(dcalc) )
{
 CH_TIMELEAF("FORT_FLUXASSEMBLE");
 FORTRAN_NAME( FLUXASSEMBLE ,fluxassemble )(
      CHFt_FRA(fluxregflux)
      ,CHFt_CONST_FRA(godunovflux)
      ,CHFt_CONST_REAL(dr)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(dcalc) );
}
#endif  // GUARDFLUXASSEMBLE 

#ifndef GUARDPRM2CONS 
#define GUARDPRM2CONS 
// Prototype for Fortran procedure prm2cons ...
//
void FORTRAN_NAME( PRM2CONS ,prm2cons )(
      CHFp_BOX(dcalc)
      ,CHFp_FRA(u)
      ,CHFp_CONST_FRA(q) );

#define FORT_PRM2CONS FORTRAN_NAME( inlinePRM2CONS, inlinePRM2CONS)
#define FORTNT_PRM2CONS FORTRAN_NAME( PRM2CONS, prm2cons)

inline void FORTRAN_NAME(inlinePRM2CONS, inlinePRM2CONS)(
      CHFp_BOX(dcalc)
      ,CHFp_FRA(u)
      ,CHFp_CONST_FRA(q) )
{
 CH_TIMELEAF("FORT_PRM2CONS");
 FORTRAN_NAME( PRM2CONS ,prm2cons )(
      CHFt_BOX(dcalc)
      ,CHFt_FRA(u)
      ,CHFt_CONST_FRA(q) );
}
#endif  // GUARDPRM2CONS 

#ifndef GUARDPOINTPRM2CONS 
#define GUARDPOINTPRM2CONS 
// Prototype for Fortran procedure pointprm2cons ...
//
void FORTRAN_NAME( POINTPRM2CONS ,pointprm2cons )(
      CHFp_VR(conserved)
      ,CHFp_VR(primitive) );

#define FORT_POINTPRM2CONS FORTRAN_NAME( inlinePOINTPRM2CONS, inlinePOINTPRM2CONS)
#define FORTNT_POINTPRM2CONS FORTRAN_NAME( POINTPRM2CONS, pointprm2cons)

inline void FORTRAN_NAME(inlinePOINTPRM2CONS, inlinePOINTPRM2CONS)(
      CHFp_VR(conserved)
      ,CHFp_VR(primitive) )
{
 CH_TIMELEAF("FORT_POINTPRM2CONS");
 FORTRAN_NAME( POINTPRM2CONS ,pointprm2cons )(
      CHFt_VR(conserved)
      ,CHFt_VR(primitive) );
}
#endif  // GUARDPOINTPRM2CONS 

#ifndef GUARDRIEMANN 
#define GUARDRIEMANN 
// Prototype for Fortran procedure riemann ...
//
void FORTRAN_NAME( RIEMANN ,riemann )(
      CHFp_BOX(d)
      ,CHFp_CONST_FRA(ql)
      ,CHFp_CONST_FRA(qr)
      ,CHFp_FRA(fluxgod)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_INT(nspec) );

#define FORT_RIEMANN FORTRAN_NAME( inlineRIEMANN, inlineRIEMANN)
#define FORTNT_RIEMANN FORTRAN_NAME( RIEMANN, riemann)

inline void FORTRAN_NAME(inlineRIEMANN, inlineRIEMANN)(
      CHFp_BOX(d)
      ,CHFp_CONST_FRA(ql)
      ,CHFp_CONST_FRA(qr)
      ,CHFp_FRA(fluxgod)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_INT(nspec) )
{
 CH_TIMELEAF("FORT_RIEMANN");
 FORTRAN_NAME( RIEMANN ,riemann )(
      CHFt_BOX(d)
      ,CHFt_CONST_FRA(ql)
      ,CHFt_CONST_FRA(qr)
      ,CHFt_FRA(fluxgod)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_INT(nspec) );
}
#endif  // GUARDRIEMANN 

#ifndef GUARDPOINTRIEMANN 
#define GUARDPOINTRIEMANN 
// Prototype for Fortran procedure pointriemann ...
//
void FORTRAN_NAME( POINTRIEMANN ,pointriemann )(
      CHFp_VR(priml)
      ,CHFp_VR(primr)
      ,CHFp_VR(fluxvec)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_INT(nspec) );

#define FORT_POINTRIEMANN FORTRAN_NAME( inlinePOINTRIEMANN, inlinePOINTRIEMANN)
#define FORTNT_POINTRIEMANN FORTRAN_NAME( POINTRIEMANN, pointriemann)

inline void FORTRAN_NAME(inlinePOINTRIEMANN, inlinePOINTRIEMANN)(
      CHFp_VR(priml)
      ,CHFp_VR(primr)
      ,CHFp_VR(fluxvec)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_INT(nspec) )
{
 CH_TIMELEAF("FORT_POINTRIEMANN");
 FORTRAN_NAME( POINTRIEMANN ,pointriemann )(
      CHFt_VR(priml)
      ,CHFt_VR(primr)
      ,CHFt_VR(fluxvec)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_INT(nspec) );
}
#endif  // GUARDPOINTRIEMANN 

#ifndef GUARDPOINTGETFLUX 
#define GUARDPOINTGETFLUX 
// Prototype for Fortran procedure pointgetflux ...
//
void FORTRAN_NAME( POINTGETFLUX ,pointgetflux )(
      CHFp_VR(flux)
      ,CHFp_VR(primitive)
      ,CHFp_CONST_INT(facedir) );

#define FORT_POINTGETFLUX FORTRAN_NAME( inlinePOINTGETFLUX, inlinePOINTGETFLUX)
#define FORTNT_POINTGETFLUX FORTRAN_NAME( POINTGETFLUX, pointgetflux)

inline void FORTRAN_NAME(inlinePOINTGETFLUX, inlinePOINTGETFLUX)(
      CHFp_VR(flux)
      ,CHFp_VR(primitive)
      ,CHFp_CONST_INT(facedir) )
{
 CH_TIMELEAF("FORT_POINTGETFLUX");
 FORTRAN_NAME( POINTGETFLUX ,pointgetflux )(
      CHFt_VR(flux)
      ,CHFt_VR(primitive)
      ,CHFt_CONST_INT(facedir) );
}
#endif  // GUARDPOINTGETFLUX 

#ifndef GUARDGETROEFLUX 
#define GUARDGETROEFLUX 
// Prototype for Fortran procedure getRoeflux ...
//
void FORTRAN_NAME( GETROEFLUX ,getroeflux )(
      CHFp_VR(flux)
      ,CHFp_VR(fluxl)
      ,CHFp_VR(fluxr)
      ,CHFp_VR(priml)
      ,CHFp_VR(primr)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_INT(nspec)
      ,CHFp_INT(i)
      ,CHFp_INT(j) );

#define FORT_GETROEFLUX FORTRAN_NAME( inlineGETROEFLUX, inlineGETROEFLUX)
#define FORTNT_GETROEFLUX FORTRAN_NAME( GETROEFLUX, getroeflux)

inline void FORTRAN_NAME(inlineGETROEFLUX, inlineGETROEFLUX)(
      CHFp_VR(flux)
      ,CHFp_VR(fluxl)
      ,CHFp_VR(fluxr)
      ,CHFp_VR(priml)
      ,CHFp_VR(primr)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_INT(nspec)
      ,CHFp_INT(i)
      ,CHFp_INT(j) )
{
 CH_TIMELEAF("FORT_GETROEFLUX");
 FORTRAN_NAME( GETROEFLUX ,getroeflux )(
      CHFt_VR(flux)
      ,CHFt_VR(fluxl)
      ,CHFt_VR(fluxr)
      ,CHFt_VR(priml)
      ,CHFt_VR(primr)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_INT(nspec)
      ,CHFt_INT(i)
      ,CHFt_INT(j) );
}
#endif  // GUARDGETROEFLUX 

#ifndef GUARDREACTIVESRC 
#define GUARDREACTIVESRC 
// Prototype for Fortran procedure reactivesrc ...
//
void FORTRAN_NAME( REACTIVESRC ,reactivesrc )(
      CHFp_BOX(dcalc)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_FRA(primitive) );

#define FORT_REACTIVESRC FORTRAN_NAME( inlineREACTIVESRC, inlineREACTIVESRC)
#define FORTNT_REACTIVESRC FORTRAN_NAME( REACTIVESRC, reactivesrc)

inline void FORTRAN_NAME(inlineREACTIVESRC, inlineREACTIVESRC)(
      CHFp_BOX(dcalc)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_FRA(primitive) )
{
 CH_TIMELEAF("FORT_REACTIVESRC");
 FORTRAN_NAME( REACTIVESRC ,reactivesrc )(
      CHFt_BOX(dcalc)
      ,CHFt_CONST_REAL(dt)
      ,CHFt_FRA(primitive) );
}
#endif  // GUARDREACTIVESRC 

#ifndef GUARDPOINTREACTIVESRC 
#define GUARDPOINTREACTIVESRC 
// Prototype for Fortran procedure pointreactivesrc ...
//
void FORTRAN_NAME( POINTREACTIVESRC ,pointreactivesrc )(
      CHFp_CONST_REAL(dt)
      ,CHFp_VR(primitive) );

#define FORT_POINTREACTIVESRC FORTRAN_NAME( inlinePOINTREACTIVESRC, inlinePOINTREACTIVESRC)
#define FORTNT_POINTREACTIVESRC FORTRAN_NAME( POINTREACTIVESRC, pointreactivesrc)

inline void FORTRAN_NAME(inlinePOINTREACTIVESRC, inlinePOINTREACTIVESRC)(
      CHFp_CONST_REAL(dt)
      ,CHFp_VR(primitive) )
{
 CH_TIMELEAF("FORT_POINTREACTIVESRC");
 FORTRAN_NAME( POINTREACTIVESRC ,pointreactivesrc )(
      CHFt_CONST_REAL(dt)
      ,CHFt_VR(primitive) );
}
#endif  // GUARDPOINTREACTIVESRC 

#ifndef GUARDDYDT 
#define GUARDDYDT 
// Prototype for Fortran procedure DYDT ...
//
void FORTRAN_NAME( DYDT ,dydt )(
      CHFp_REAL(p)
      ,CHFp_REAL(t)
      ,CHFp_REAL(rho)
      ,CHFp_VR(MassFrac)
      ,CHFp_VR(ydot)
      ,CHFp_INT(n) );

#define FORT_DYDT FORTRAN_NAME( inlineDYDT, inlineDYDT)
#define FORTNT_DYDT FORTRAN_NAME( DYDT, dydt)

inline void FORTRAN_NAME(inlineDYDT, inlineDYDT)(
      CHFp_REAL(p)
      ,CHFp_REAL(t)
      ,CHFp_REAL(rho)
      ,CHFp_VR(MassFrac)
      ,CHFp_VR(ydot)
      ,CHFp_INT(n) )
{
 CH_TIMELEAF("FORT_DYDT");
 FORTRAN_NAME( DYDT ,dydt )(
      CHFt_REAL(p)
      ,CHFt_REAL(t)
      ,CHFt_REAL(rho)
      ,CHFt_VR(MassFrac)
      ,CHFt_VR(ydot)
      ,CHFt_INT(n) );
}
#endif  // GUARDDYDT 

#ifndef GUARDGETRHOCV 
#define GUARDGETRHOCV 
// Prototype for Fortran procedure getrhocv ...
//
void FORTRAN_NAME( GETRHOCV ,getrhocv )(
      CHFp_BOX(dcalc)
      ,CHFp_CONST_FRA(dense)
      ,CHFp_CONST_FRA(temperature)
      ,CHFp_CONST_FRA(massFrac)
      ,CHFp_FRA(rhoCv) );

#define FORT_GETRHOCV FORTRAN_NAME( inlineGETRHOCV, inlineGETRHOCV)
#define FORTNT_GETRHOCV FORTRAN_NAME( GETRHOCV, getrhocv)

inline void FORTRAN_NAME(inlineGETRHOCV, inlineGETRHOCV)(
      CHFp_BOX(dcalc)
      ,CHFp_CONST_FRA(dense)
      ,CHFp_CONST_FRA(temperature)
      ,CHFp_CONST_FRA(massFrac)
      ,CHFp_FRA(rhoCv) )
{
 CH_TIMELEAF("FORT_GETRHOCV");
 FORTRAN_NAME( GETRHOCV ,getrhocv )(
      CHFt_BOX(dcalc)
      ,CHFt_CONST_FRA(dense)
      ,CHFt_CONST_FRA(temperature)
      ,CHFt_CONST_FRA(massFrac)
      ,CHFt_FRA(rhoCv) );
}
#endif  // GUARDGETRHOCV 

#ifndef GUARDPOINTGETRHOCV 
#define GUARDPOINTGETRHOCV 
// Prototype for Fortran procedure pointgetrhocv ...
//
void FORTRAN_NAME( POINTGETRHOCV ,pointgetrhocv )(
      CHFp_CONST_REAL(rho)
      ,CHFp_CONST_REAL(temp)
      ,CHFp_CONST_VR(massFrac)
      ,CHFp_REAL(rhoCv) );

#define FORT_POINTGETRHOCV FORTRAN_NAME( inlinePOINTGETRHOCV, inlinePOINTGETRHOCV)
#define FORTNT_POINTGETRHOCV FORTRAN_NAME( POINTGETRHOCV, pointgetrhocv)

inline void FORTRAN_NAME(inlinePOINTGETRHOCV, inlinePOINTGETRHOCV)(
      CHFp_CONST_REAL(rho)
      ,CHFp_CONST_REAL(temp)
      ,CHFp_CONST_VR(massFrac)
      ,CHFp_REAL(rhoCv) )
{
 CH_TIMELEAF("FORT_POINTGETRHOCV");
 FORTRAN_NAME( POINTGETRHOCV ,pointgetrhocv )(
      CHFt_CONST_REAL(rho)
      ,CHFt_CONST_REAL(temp)
      ,CHFt_CONST_VR(massFrac)
      ,CHFt_REAL(rhoCv) );
}
#endif  // GUARDPOINTGETRHOCV 

#ifndef GUARDGETSPECMASSDIFFCOEFF 
#define GUARDGETSPECMASSDIFFCOEFF 
// Prototype for Fortran procedure getspecmassdiffcoeff ...
//
void FORTRAN_NAME( GETSPECMASSDIFFCOEFF ,getspecmassdiffcoeff )(
      CHFp_FRA1(bcospec)
      ,CHFp_FRA(rhscospec)
      ,CHFp_CONST_FRA1(pressure)
      ,CHFp_CONST_FRA1(temperature)
      ,CHFp_CONST_FRA(specdense)
      ,CHFp_CONST_INT(ispec)
      ,CHFp_BOX(dcalc) );

#define FORT_GETSPECMASSDIFFCOEFF FORTRAN_NAME( inlineGETSPECMASSDIFFCOEFF, inlineGETSPECMASSDIFFCOEFF)
#define FORTNT_GETSPECMASSDIFFCOEFF FORTRAN_NAME( GETSPECMASSDIFFCOEFF, getspecmassdiffcoeff)

inline void FORTRAN_NAME(inlineGETSPECMASSDIFFCOEFF, inlineGETSPECMASSDIFFCOEFF)(
      CHFp_FRA1(bcospec)
      ,CHFp_FRA(rhscospec)
      ,CHFp_CONST_FRA1(pressure)
      ,CHFp_CONST_FRA1(temperature)
      ,CHFp_CONST_FRA(specdense)
      ,CHFp_CONST_INT(ispec)
      ,CHFp_BOX(dcalc) )
{
 CH_TIMELEAF("FORT_GETSPECMASSDIFFCOEFF");
 FORTRAN_NAME( GETSPECMASSDIFFCOEFF ,getspecmassdiffcoeff )(
      CHFt_FRA1(bcospec)
      ,CHFt_FRA(rhscospec)
      ,CHFt_CONST_FRA1(pressure)
      ,CHFt_CONST_FRA1(temperature)
      ,CHFt_CONST_FRA(specdense)
      ,CHFt_CONST_INT(ispec)
      ,CHFt_BOX(dcalc) );
}
#endif  // GUARDGETSPECMASSDIFFCOEFF 

#ifndef GUARDPOINTGETSPECMASSDIFFCOEFF 
#define GUARDPOINTGETSPECMASSDIFFCOEFF 
// Prototype for Fortran procedure pointgetspecmassdiffcoeff ...
//
void FORTRAN_NAME( POINTGETSPECMASSDIFFCOEFF ,pointgetspecmassdiffcoeff )(
      CHFp_REAL(bcospec)
      ,CHFp_VR(rhscospec)
      ,CHFp_CONST_REAL(pres)
      ,CHFp_CONST_REAL(temp)
      ,CHFp_CONST_VR(specdense)
      ,CHFp_CONST_INT(ispec) );

#define FORT_POINTGETSPECMASSDIFFCOEFF FORTRAN_NAME( inlinePOINTGETSPECMASSDIFFCOEFF, inlinePOINTGETSPECMASSDIFFCOEFF)
#define FORTNT_POINTGETSPECMASSDIFFCOEFF FORTRAN_NAME( POINTGETSPECMASSDIFFCOEFF, pointgetspecmassdiffcoeff)

inline void FORTRAN_NAME(inlinePOINTGETSPECMASSDIFFCOEFF, inlinePOINTGETSPECMASSDIFFCOEFF)(
      CHFp_REAL(bcospec)
      ,CHFp_VR(rhscospec)
      ,CHFp_CONST_REAL(pres)
      ,CHFp_CONST_REAL(temp)
      ,CHFp_CONST_VR(specdense)
      ,CHFp_CONST_INT(ispec) )
{
 CH_TIMELEAF("FORT_POINTGETSPECMASSDIFFCOEFF");
 FORTRAN_NAME( POINTGETSPECMASSDIFFCOEFF ,pointgetspecmassdiffcoeff )(
      CHFt_REAL(bcospec)
      ,CHFt_VR(rhscospec)
      ,CHFt_CONST_REAL(pres)
      ,CHFt_CONST_REAL(temp)
      ,CHFt_CONST_VR(specdense)
      ,CHFt_CONST_INT(ispec) );
}
#endif  // GUARDPOINTGETSPECMASSDIFFCOEFF 

#ifndef GUARDGETVISCOSITY 
#define GUARDGETVISCOSITY 
// Prototype for Fortran procedure getviscosity ...
//
void FORTRAN_NAME( GETVISCOSITY ,getviscosity )(
      CHFp_FRA1(eta)
      ,CHFp_FRA1(lambda)
      ,CHFp_CONST_FRA1(temperature)
      ,CHFp_CONST_FRA(specdense)
      ,CHFp_BOX(dcalc) );

#define FORT_GETVISCOSITY FORTRAN_NAME( inlineGETVISCOSITY, inlineGETVISCOSITY)
#define FORTNT_GETVISCOSITY FORTRAN_NAME( GETVISCOSITY, getviscosity)

inline void FORTRAN_NAME(inlineGETVISCOSITY, inlineGETVISCOSITY)(
      CHFp_FRA1(eta)
      ,CHFp_FRA1(lambda)
      ,CHFp_CONST_FRA1(temperature)
      ,CHFp_CONST_FRA(specdense)
      ,CHFp_BOX(dcalc) )
{
 CH_TIMELEAF("FORT_GETVISCOSITY");
 FORTRAN_NAME( GETVISCOSITY ,getviscosity )(
      CHFt_FRA1(eta)
      ,CHFt_FRA1(lambda)
      ,CHFt_CONST_FRA1(temperature)
      ,CHFt_CONST_FRA(specdense)
      ,CHFt_BOX(dcalc) );
}
#endif  // GUARDGETVISCOSITY 

#ifndef GUARDPOINTGETVISCOSITY 
#define GUARDPOINTGETVISCOSITY 
// Prototype for Fortran procedure pointgetviscosity ...
//
void FORTRAN_NAME( POINTGETVISCOSITY ,pointgetviscosity )(
      CHFp_REAL(eta)
      ,CHFp_REAL(lambda)
      ,CHFp_CONST_REAL(temp)
      ,CHFp_CONST_VR(specdense) );

#define FORT_POINTGETVISCOSITY FORTRAN_NAME( inlinePOINTGETVISCOSITY, inlinePOINTGETVISCOSITY)
#define FORTNT_POINTGETVISCOSITY FORTRAN_NAME( POINTGETVISCOSITY, pointgetviscosity)

inline void FORTRAN_NAME(inlinePOINTGETVISCOSITY, inlinePOINTGETVISCOSITY)(
      CHFp_REAL(eta)
      ,CHFp_REAL(lambda)
      ,CHFp_CONST_REAL(temp)
      ,CHFp_CONST_VR(specdense) )
{
 CH_TIMELEAF("FORT_POINTGETVISCOSITY");
 FORTRAN_NAME( POINTGETVISCOSITY ,pointgetviscosity )(
      CHFt_REAL(eta)
      ,CHFt_REAL(lambda)
      ,CHFt_CONST_REAL(temp)
      ,CHFt_CONST_VR(specdense) );
}
#endif  // GUARDPOINTGETVISCOSITY 

#ifndef GUARDGETCONDUCTIVITY 
#define GUARDGETCONDUCTIVITY 
// Prototype for Fortran procedure getconductivity ...
//
void FORTRAN_NAME( GETCONDUCTIVITY ,getconductivity )(
      CHFp_FRA1(kappa)
      ,CHFp_CONST_FRA1(temperature)
      ,CHFp_CONST_FRA(specdense)
      ,CHFp_BOX(dcalc) );

#define FORT_GETCONDUCTIVITY FORTRAN_NAME( inlineGETCONDUCTIVITY, inlineGETCONDUCTIVITY)
#define FORTNT_GETCONDUCTIVITY FORTRAN_NAME( GETCONDUCTIVITY, getconductivity)

inline void FORTRAN_NAME(inlineGETCONDUCTIVITY, inlineGETCONDUCTIVITY)(
      CHFp_FRA1(kappa)
      ,CHFp_CONST_FRA1(temperature)
      ,CHFp_CONST_FRA(specdense)
      ,CHFp_BOX(dcalc) )
{
 CH_TIMELEAF("FORT_GETCONDUCTIVITY");
 FORTRAN_NAME( GETCONDUCTIVITY ,getconductivity )(
      CHFt_FRA1(kappa)
      ,CHFt_CONST_FRA1(temperature)
      ,CHFt_CONST_FRA(specdense)
      ,CHFt_BOX(dcalc) );
}
#endif  // GUARDGETCONDUCTIVITY 

#ifndef GUARDPOINTGETCONDUCTIVITY 
#define GUARDPOINTGETCONDUCTIVITY 
// Prototype for Fortran procedure pointgetconductivity ...
//
void FORTRAN_NAME( POINTGETCONDUCTIVITY ,pointgetconductivity )(
      CHFp_REAL(kappa)
      ,CHFp_CONST_REAL(temp)
      ,CHFp_CONST_VR(specdense) );

#define FORT_POINTGETCONDUCTIVITY FORTRAN_NAME( inlinePOINTGETCONDUCTIVITY, inlinePOINTGETCONDUCTIVITY)
#define FORTNT_POINTGETCONDUCTIVITY FORTRAN_NAME( POINTGETCONDUCTIVITY, pointgetconductivity)

inline void FORTRAN_NAME(inlinePOINTGETCONDUCTIVITY, inlinePOINTGETCONDUCTIVITY)(
      CHFp_REAL(kappa)
      ,CHFp_CONST_REAL(temp)
      ,CHFp_CONST_VR(specdense) )
{
 CH_TIMELEAF("FORT_POINTGETCONDUCTIVITY");
 FORTRAN_NAME( POINTGETCONDUCTIVITY ,pointgetconductivity )(
      CHFt_REAL(kappa)
      ,CHFt_CONST_REAL(temp)
      ,CHFt_CONST_VR(specdense) );
}
#endif  // GUARDPOINTGETCONDUCTIVITY 

#ifndef GUARDFILLDKMATRIX 
#define GUARDFILLDKMATRIX 
// Prototype for Fortran procedure filldkmatrix ...
//
void FORTRAN_NAME( FILLDKMATRIX ,filldkmatrix )(
      CHFp_BOX(dcalc)
      ,CHFp_CONST_FRA(pressure)
      ,CHFp_CONST_FRA(temperature)
      ,CHFp_CONST_FRA(massFrac)
      ,CHFp_FRA(DiffCoeff) );

#define FORT_FILLDKMATRIX FORTRAN_NAME( inlineFILLDKMATRIX, inlineFILLDKMATRIX)
#define FORTNT_FILLDKMATRIX FORTRAN_NAME( FILLDKMATRIX, filldkmatrix)

inline void FORTRAN_NAME(inlineFILLDKMATRIX, inlineFILLDKMATRIX)(
      CHFp_BOX(dcalc)
      ,CHFp_CONST_FRA(pressure)
      ,CHFp_CONST_FRA(temperature)
      ,CHFp_CONST_FRA(massFrac)
      ,CHFp_FRA(DiffCoeff) )
{
 CH_TIMELEAF("FORT_FILLDKMATRIX");
 FORTRAN_NAME( FILLDKMATRIX ,filldkmatrix )(
      CHFt_BOX(dcalc)
      ,CHFt_CONST_FRA(pressure)
      ,CHFt_CONST_FRA(temperature)
      ,CHFt_CONST_FRA(massFrac)
      ,CHFt_FRA(DiffCoeff) );
}
#endif  // GUARDFILLDKMATRIX 

#ifndef GUARDPOINTFILLDKMATRIX 
#define GUARDPOINTFILLDKMATRIX 
// Prototype for Fortran procedure pointfilldkmatrix ...
//
void FORTRAN_NAME( POINTFILLDKMATRIX ,pointfilldkmatrix )(
      CHFp_CONST_REAL(pressure)
      ,CHFp_CONST_REAL(temperature)
      ,CHFp_CONST_VR(massFrac)
      ,CHFp_VR(DiffCoeff) );

#define FORT_POINTFILLDKMATRIX FORTRAN_NAME( inlinePOINTFILLDKMATRIX, inlinePOINTFILLDKMATRIX)
#define FORTNT_POINTFILLDKMATRIX FORTRAN_NAME( POINTFILLDKMATRIX, pointfilldkmatrix)

inline void FORTRAN_NAME(inlinePOINTFILLDKMATRIX, inlinePOINTFILLDKMATRIX)(
      CHFp_CONST_REAL(pressure)
      ,CHFp_CONST_REAL(temperature)
      ,CHFp_CONST_VR(massFrac)
      ,CHFp_VR(DiffCoeff) )
{
 CH_TIMELEAF("FORT_POINTFILLDKMATRIX");
 FORTRAN_NAME( POINTFILLDKMATRIX ,pointfilldkmatrix )(
      CHFt_CONST_REAL(pressure)
      ,CHFt_CONST_REAL(temperature)
      ,CHFt_CONST_VR(massFrac)
      ,CHFt_VR(DiffCoeff) );
}
#endif  // GUARDPOINTFILLDKMATRIX 

#ifndef GUARDFILLVISCANDCONDUCTIVECOEFF 
#define GUARDFILLVISCANDCONDUCTIVECOEFF 
// Prototype for Fortran procedure fillviscandconductivecoeff ...
//
void FORTRAN_NAME( FILLVISCANDCONDUCTIVECOEFF ,fillviscandconductivecoeff )(
      CHFp_BOX(dcalc)
      ,CHFp_CONST_FRA(massFrac)
      ,CHFp_CONST_FRA(dense)
      ,CHFp_CONST_FRA(temperature)
      ,CHFp_FRA(aco)
      ,CHFp_FRA(mu)
      ,CHFp_FRA(lambda)
      ,CHFp_FRA(kappa) );

#define FORT_FILLVISCANDCONDUCTIVECOEFF FORTRAN_NAME( inlineFILLVISCANDCONDUCTIVECOEFF, inlineFILLVISCANDCONDUCTIVECOEFF)
#define FORTNT_FILLVISCANDCONDUCTIVECOEFF FORTRAN_NAME( FILLVISCANDCONDUCTIVECOEFF, fillviscandconductivecoeff)

inline void FORTRAN_NAME(inlineFILLVISCANDCONDUCTIVECOEFF, inlineFILLVISCANDCONDUCTIVECOEFF)(
      CHFp_BOX(dcalc)
      ,CHFp_CONST_FRA(massFrac)
      ,CHFp_CONST_FRA(dense)
      ,CHFp_CONST_FRA(temperature)
      ,CHFp_FRA(aco)
      ,CHFp_FRA(mu)
      ,CHFp_FRA(lambda)
      ,CHFp_FRA(kappa) )
{
 CH_TIMELEAF("FORT_FILLVISCANDCONDUCTIVECOEFF");
 FORTRAN_NAME( FILLVISCANDCONDUCTIVECOEFF ,fillviscandconductivecoeff )(
      CHFt_BOX(dcalc)
      ,CHFt_CONST_FRA(massFrac)
      ,CHFt_CONST_FRA(dense)
      ,CHFt_CONST_FRA(temperature)
      ,CHFt_FRA(aco)
      ,CHFt_FRA(mu)
      ,CHFt_FRA(lambda)
      ,CHFt_FRA(kappa) );
}
#endif  // GUARDFILLVISCANDCONDUCTIVECOEFF 

#ifndef GUARDPOINTFILLVISCANDCONDUCTIVECOEFF 
#define GUARDPOINTFILLVISCANDCONDUCTIVECOEFF 
// Prototype for Fortran procedure pointfillviscandconductivecoeff ...
//
void FORTRAN_NAME( POINTFILLVISCANDCONDUCTIVECOEFF ,pointfillviscandconductivecoeff )(
      CHFp_CONST_VR(massFrac)
      ,CHFp_CONST_REAL(dense)
      ,CHFp_CONST_REAL(temperature)
      ,CHFp_REAL(aco)
      ,CHFp_REAL(mu)
      ,CHFp_REAL(lambda)
      ,CHFp_REAL(kappa) );

#define FORT_POINTFILLVISCANDCONDUCTIVECOEFF FORTRAN_NAME( inlinePOINTFILLVISCANDCONDUCTIVECOEFF, inlinePOINTFILLVISCANDCONDUCTIVECOEFF)
#define FORTNT_POINTFILLVISCANDCONDUCTIVECOEFF FORTRAN_NAME( POINTFILLVISCANDCONDUCTIVECOEFF, pointfillviscandconductivecoeff)

inline void FORTRAN_NAME(inlinePOINTFILLVISCANDCONDUCTIVECOEFF, inlinePOINTFILLVISCANDCONDUCTIVECOEFF)(
      CHFp_CONST_VR(massFrac)
      ,CHFp_CONST_REAL(dense)
      ,CHFp_CONST_REAL(temperature)
      ,CHFp_REAL(aco)
      ,CHFp_REAL(mu)
      ,CHFp_REAL(lambda)
      ,CHFp_REAL(kappa) )
{
 CH_TIMELEAF("FORT_POINTFILLVISCANDCONDUCTIVECOEFF");
 FORTRAN_NAME( POINTFILLVISCANDCONDUCTIVECOEFF ,pointfillviscandconductivecoeff )(
      CHFt_CONST_VR(massFrac)
      ,CHFt_CONST_REAL(dense)
      ,CHFt_CONST_REAL(temperature)
      ,CHFt_REAL(aco)
      ,CHFt_REAL(mu)
      ,CHFt_REAL(lambda)
      ,CHFt_REAL(kappa) );
}
#endif  // GUARDPOINTFILLVISCANDCONDUCTIVECOEFF 

#ifndef GUARDBRENT_METHOD 
#define GUARDBRENT_METHOD 
// Prototype for Fortran procedure brent_method ...
//
void FORTRAN_NAME( BRENT_METHOD ,brent_method )(
      CHFp_REAL(T)
      ,CHFp_REAL(a)
      ,CHFp_REAL(b)
      ,CHFp_REAL(cvtemp)
      ,CHFp_INT(success) );

#define FORT_BRENT_METHOD FORTRAN_NAME( inlineBRENT_METHOD, inlineBRENT_METHOD)
#define FORTNT_BRENT_METHOD FORTRAN_NAME( BRENT_METHOD, brent_method)

inline void FORTRAN_NAME(inlineBRENT_METHOD, inlineBRENT_METHOD)(
      CHFp_REAL(T)
      ,CHFp_REAL(a)
      ,CHFp_REAL(b)
      ,CHFp_REAL(cvtemp)
      ,CHFp_INT(success) )
{
 CH_TIMELEAF("FORT_BRENT_METHOD");
 FORTRAN_NAME( BRENT_METHOD ,brent_method )(
      CHFt_REAL(T)
      ,CHFt_REAL(a)
      ,CHFt_REAL(b)
      ,CHFt_REAL(cvtemp)
      ,CHFt_INT(success) );
}
#endif  // GUARDBRENT_METHOD 

}

#endif
