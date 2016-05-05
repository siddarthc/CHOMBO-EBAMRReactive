#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EBPatchReactive.H"
#include "EBPatchReactiveF_F.H"
#include "EBLGIntegrator.H"

#include "EBArith.H"
#include "PolyGeom.H"
#include "EBDebugOut.H"
#include "DebugOut.H"
#include "VoFIterator.H"
#include "FaceIterator.H"
#include "EBLoHiCenter.H"
#include "Stencils.H"
#include "EBArith.H"
#include "EBAMRIO.H"
#include "PolyGeom.H"
#include "ParmParse.H"
#include "parstream.H"
#include <iomanip>
#include <cmath>
#include <cstdio>
#include <string>
#include <sstream>
#include "CH_Timer.H"
#include "FabDataOps.H"
#include "TensorCFInterp.H" // just for gradIndex
#include "NamespaceHeader.H"

bool EBPatchReactive::s_conservativeSource  = true;
bool EBPatchReactive::s_verbose  = false;
int  EBPatchReactive::s_curLevel = -1;
int  EBPatchReactive::s_curComp  = -1;
int  EBPatchReactive::s_doingVel  = -1;
int  EBPatchReactive::s_doingAdvVel  = -1;
Real EBPatchReactive::s_maxWaveSpeed  = -1.0;
IntVect EBPatchReactive::s_maxWaveSpeedIV   = IntVect(D_DECL(-1,-1,-1));
IntVect EBPatchReactive::s_debugIV          = IntVect(D_DECL(16, 3, 0));
int     EBPatchReactive::s_whichLev=-1;

//-----------------------------------------------------------------------

EBPatchReactive::EBPatchReactive()
{
  m_isDefined  = false;
  m_isBCSet    = false;
  m_isBoxSet   = false;
  m_useAgg     = false;
  m_bc         = NULL;
}

EBPatchReactive::~EBPatchReactive()
{
  if (m_bc != NULL)
    {
      delete m_bc;
     }
}
/******/
void 
EBPatchReactive::
setGamma(const Real& a_gamma)
{
  m_gamma = a_gamma;
  m_isGammaSet = true;
}
/*****************************/
Real
EBPatchReactive::
getGamma() const
{
  CH_assert(m_isGammaSet);
  return m_gamma;
}
/******/
void
EBPatchReactive::
setSpecHeat(const Real& a_specHeat)
{
  m_specHeat = a_specHeat;
  m_isSpecHeatSet = true;
}
/*****************************/
Real
EBPatchReactive::
getSpecHeat() const
{
  CH_assert(m_isSpecHeatSet);
  return m_specHeat;
}
/******/
void 
EBPatchReactive::
setnSpecies(const int& a_nSpecies)
{
  m_nSpec = a_nSpecies;
  m_isnSpeciesSet = true;
}
/*****************************/
int 
EBPatchReactive::
getnSpecies() const
{
  CH_assert(m_isnSpeciesSet);
  return m_nSpec;
}
/******/
void 
EBPatchReactive::
setEBPhysIBC(const EBPhysIBCFactory& a_bcFact)
{
  // Delete old IBC object if any
  if (m_bc != NULL)
    {
      delete m_bc;
    }

  m_bc = a_bcFact.create();
  if (m_isDefined)
  m_bc->define(m_domain,m_dx);
  m_isBCSet = true;
}
/*****************************/
const EBPhysIBC*
EBPatchReactive::
getEBPhysIBC() const
{
  CH_assert(m_isBCSet);
  return m_bc;
}
/*****************************/
void
EBPatchReactive::
doRZCoords(bool a_doRZCoords)
{
  m_doRZCoords = a_doRZCoords;
}
/*****************************/
void
EBPatchReactive::
artificialViscosity(bool a_useArtificialVisc)
{
  m_useArtificialVisc = a_useArtificialVisc;
  m_isArtViscSet = true;
}
/*****************************/
Real
EBPatchReactive::
artificialViscosityCoefficient() const
{
  ParmParse pp;
  Real retval;
  pp.get("artificial_viscosity", retval);
  return retval;
}
/*****************************/
void
EBPatchReactive::
setSlopeParameters(bool a_useFourthOrderSlopes,
                   bool a_useZeroSlopes,
                   bool a_useFlattening,
                   bool a_useLimiting)
{
  // Slope flattening is only allowed with 4th order slopes
  CH_assert(a_useFourthOrderSlopes || !a_useFlattening);

  // Store the slope computation parameters
  m_useFourthOrderSlopes = a_useFourthOrderSlopes;
  m_useZeroSlopes        = a_useZeroSlopes;
  m_useFlattening        = a_useFlattening;
  m_useLimiting          = a_useLimiting;

  m_isSlopeSet = true;
}
/*****************************/
void
EBPatchReactive::
define(const ProblemDomain&  a_domain,
       const RealVect& a_dx,
       bool a_useAgg)
{
  // Store the domain and grid spacing
  m_isDefined = true;
  m_domain = a_domain;
  // Set the domain and grid spacing in the boundary condition object
  m_dx = a_dx;
  m_useAgg = a_useAgg;
  //figure out dxScale (used for unscaling the boundary area)
  {
    Real maxDx=0.0,dv=1.0;
    for (int idir = 0; idir < SpaceDim; idir++)
      {
        dv *= m_dx[idir];
        if ( m_dx[idir] > maxDx )
          {
            maxDx = m_dx[idir];
          }
      }
    m_dxScale = pow(maxDx,SpaceDim-1)/dv;
  }

  if (m_isBCSet)
    m_bc->define(a_domain, a_dx);
}
/******/
int 
EBPatchReactive::
numConserved() const 
{
   return CNUM+m_nSpec;
}
/******/
int
EBPatchReactive::
numFluxes() const
{
   return FNUM+m_nSpec;
}
/******/
int
EBPatchReactive::
numPrimitives() const
{
   return QNUM+m_nSpec;
}
/******/
int
EBPatchReactive::
numSlopes() const
{
  return QSLOPE+m_nSpec;
}
/******/
Vector<string>
EBPatchReactive::
stateNames()
{
  Vector<string> retval;
  
  retval.push_back("mixture mass-density");   
  if (m_doRZCoords)
    {
      retval.push_back("r-momentum");
      if (SpaceDim >= 2)
        {
         retval.push_back("z-momentum");
        }
      if (SpaceDim >= 3)
        {
          MayDay::Error();
        }
     }
  else
    {
       retval.push_back("x-momentum");
       if (SpaceDim >= 2)
         {
           retval.push_back("y-momentum");
         }
       if (SpaceDim >= 3)
         {
           retval.push_back("z-momentum");
         }
    }
  retval.push_back("energy_density");

  std::string str;
  std::string str1 ("mass-density of species ");
  for (int ivar = 0; ivar < m_nSpec; ivar++)
     {
       str = static_cast<ostringstream*>( &(ostringstream() << ivar) )->str();
       // std::stringstream ss;
       // ss << ivar;
       str = str1 + str;
       retval.push_back(str);
     }

  return retval;
}
/******/
Vector<string>
EBPatchReactive::
primNames()
{
  Vector<string> retval;

  ParmParse pp;
  int logflag = 0;
  if (pp.contains("logflag"))
    {
      pp.get("logflag", logflag);
    }
  if (logflag == 1)
    {
      retval.push_back("log10density");
    }
  else
    {
      retval.push_back("density");
    }
  if (m_doRZCoords)
    {
      retval.push_back("r-velocity");
      retval.push_back("z-velocity");
    }
  else
    {
      retval.push_back("x-velocity");
      retval.push_back("y-velocity");
    }

#if CH_SPACEDIM==3
  retval.push_back("z-velocity");
#endif

  if (logflag == 1)
    {
      retval.push_back("log10pressure");
      retval.push_back("log10entropy");
    }
  else
    {
      retval.push_back("pressure");
      retval.push_back("entropy");
    }
//  retval.push_back("internal_energy");
//  retval.push_back("cv_temperature");
  retval.push_back("soundspeed");
  retval.push_back("temperature");

//#ifdef MODIANO_PROBLEM
//  retval.push_back("modiano-velocity-axial");
//  retval.push_back("modiano-velocity-tangent0");

//#if CH_SPACEDIM==3
//  retval.push_back("modiano-velocity-tangent1");
//#endif
//#endif

  std::string str;
  std::string str1 ("mass-fraction of species ");
  for (int ivar = 0; ivar < m_nSpec; ivar++)
     {
       str = static_cast<ostringstream*>( &(ostringstream() << ivar) )->str();
       str = str1 + str;
       retval.push_back(str);
     }

  return retval;
}
/******/
int
EBPatchReactive::
densityIndex() const
{
  return QRHO;
}
/******/
/******/
int
EBPatchReactive::
densityIndexC() const
{
  return CRHO;
}
/******/
Interval
EBPatchReactive::
velocityInterval() const
{
#if CH_SPACEDIM==2
  Interval retval(QVELX, QVELY);
#elif CH_SPACEDIM==3
  Interval retval(QVELX, QVELZ);
#else
  bogus_spacedim();
#endif
  return retval;
}
/******/
Interval
EBPatchReactive::
momentumInterval() const
{
#if CH_SPACEDIM==2
  Interval retval(CMOMX, CMOMY);
#elif CH_SPACEDIM==3
  Interval retval(CMOMX, CMOMZ);
#else
  bogus_spacedim();
#endif
  return retval;
}
/******/
Interval
EBPatchReactive::
speciesMassFracInterval() const
{
  Interval retval(QSPEC1,QSPEC1+m_nSpec-1);
  return retval;
}
/******/
Interval
EBPatchReactive::
speciesDenseInterval() const
{
  Interval retval(CSPEC1,CSPEC1+m_nSpec-1);
  return retval;
}
/******/
int
EBPatchReactive::
pressureIndex() const
{
  return QPRES;
}
/******/
int
EBPatchReactive::
temperatureIndex() const
{
  return QTEMP;
}
/******/
/******/
int
EBPatchReactive::
energyIndexC() const
{
  return CENG;
}
/******/
int
EBPatchReactive::
spec1MassFracIndex() const
{
  return QSPEC1;
}
/******/
int
EBPatchReactive::
spec1DenseIndex() const
{
  return CSPEC1;
}
/******/
int
EBPatchReactive::
bulkModulusIndex() const
{
  //phil said this was OK 1-2-2002
  //MayDay::Warning("returning pressure for bulk modulus");
  return QPRES;
}
/******/
bool
EBPatchReactive::
usesFlattening() const
{
  CH_assert(m_isSlopeSet);

  return m_useFlattening;
}
/******/
bool
EBPatchReactive::
usesArtificialViscosity() const
{
  return m_useArtificialVisc;
}
/******/
bool
EBPatchReactive::
usesFourthOrderSlopes() const
{
  CH_assert(m_isSlopeSet);
  return m_useFourthOrderSlopes;
}
/******/
bool
EBPatchReactive::
usesZeroSlopes() const
{
  CH_assert(m_isSlopeSet);
  return m_useZeroSlopes;
}
/******/
bool
EBPatchReactive::
isDefined() const
{
  return m_isDefined && m_isBCSet && m_isBoxSet && m_isSlopeSet && m_isArtViscSet;
}
/******/
void
EBPatchReactive::
setValidBox(const Box&        a_validBox,
            const EBISBox&    a_ebisBox,
            const IntVectSet& a_coarseFineIVS,
            const Real&       a_cumulativeTime,
            const Real&       a_timestep)
{
  m_isBoxSet = true;
  m_dt       = a_timestep;
  m_time     = a_cumulativeTime;
  m_validBox = a_validBox;
  m_ebisBox  = a_ebisBox;

  //define the interpolation stencils
  IntVectSet ivsIrreg = m_ebisBox.getIrregIVS(m_validBox);
  FaceStop::WhichFaces facestop = FaceStop::SurroundingWithBoundary;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      BaseIFFAB<FaceStencil>& stenDir = m_interpStencils[idir];
      stenDir.define(ivsIrreg, m_ebisBox.getEBGraph(), idir, 1);
      for (FaceIterator faceit(ivsIrreg, m_ebisBox.getEBGraph(), idir, facestop);
          faceit.ok(); ++faceit)
        {
          const  FaceIndex& face = faceit();
          FaceStencil sten = EBArith::getInterpStencil(face, a_coarseFineIVS,
                                                       m_ebisBox, m_domain.domainBox());
          stenDir(face, 0) = sten;
        }
    }
  //I know that a lot of these objects are slightly larger than they
  //need to be.   Be warned, however, that it is
  //really important that all these things are the same size.
  //This allows me to make all kinds of wacky stenciling assumptions down
  //the road.  Mess with these definitions at thy peril.
  m_validBoxG4  = grow(m_validBox, 4);
  m_validBoxG4 &= m_domain;
  int numFlux = numFluxes();
  int numPrim = numPrimitives();
  m_primState.define(   m_ebisBox, m_validBoxG4, numPrim);
  m_primMinuTemp.define(m_ebisBox, m_validBoxG4, numPrim);
  m_primPlusTemp.define(m_ebisBox, m_validBoxG4, numPrim);
  m_primGdnv.define(m_ebisBox, m_validBoxG4, numPrim);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      m_primPlus[idir].define(m_ebisBox, m_validBoxG4,       numPrim);
      m_primMinu[idir].define(m_ebisBox, m_validBoxG4,       numPrim);
      m_fluxOne [idir].define(m_ebisBox, m_validBoxG4, idir, numFlux);

      IntVectSet irregIVSPlus, irregIVSMinu;
      computeCoveredFaces(m_coveredFacePlusG4[idir],
                          m_coveredSetsPlusG4[idir],
                          irregIVSPlus,
                          idir, Side::Hi, m_validBoxG4);
      computeCoveredFaces(m_coveredFaceMinuG4[idir],
                          m_coveredSetsMinuG4[idir],
                          irregIVSMinu,
                          idir, Side::Lo, m_validBoxG4);

      m_extendStatePlusG4  [idir].define(m_coveredSetsPlusG4[idir], m_ebisBox.getEBGraph(), numPrim);
      m_extendStateMinuG4  [idir].define(m_coveredSetsMinuG4[idir], m_ebisBox.getEBGraph(), numPrim);
      m_extendStateNormPlus[idir].define(m_coveredSetsPlusG4[idir], m_ebisBox.getEBGraph(), numPrim);
      m_extendStateNormMinu[idir].define(m_coveredSetsMinuG4[idir], m_ebisBox.getEBGraph(), numPrim);
      m_coveredFluxPlusG4  [idir].define(m_coveredSetsPlusG4[idir], m_ebisBox.getEBGraph(), numFlux);
      m_coveredFluxMinuG4  [idir].define(m_coveredSetsMinuG4[idir], m_ebisBox.getEBGraph(), numFlux);
      m_coveredFluxNormPlus[idir].define(m_coveredSetsPlusG4[idir], m_ebisBox.getEBGraph(), numFlux);
      m_coveredFluxNormMinu[idir].define(m_coveredSetsMinuG4[idir], m_ebisBox.getEBGraph(), numFlux);

      m_extendStatePlusG4[idir].setVal(0.);
      m_extendStateMinuG4[idir].setVal(0.);
      m_coveredFluxPlusG4[idir].setVal(0.);
      m_coveredFluxMinuG4[idir].setVal(0.);

      if (SpaceDim==3)
        {
          for (int jdir = 0; jdir < SpaceDim; jdir++)
            {
              m_fluxTwo[idir][jdir].define(m_ebisBox, m_validBoxG4, idir, numFlux);

              m_extendStatePlus3D[idir][jdir].define(m_coveredSetsPlusG4[idir], m_ebisBox.getEBGraph(), numPrim);
              m_extendStateMinu3D[idir][jdir].define(m_coveredSetsMinuG4[idir], m_ebisBox.getEBGraph(), numPrim);
              m_coveredFluxPlus3D[idir][jdir].define(m_coveredSetsPlusG4[idir], m_ebisBox.getEBGraph(), numFlux);
              m_coveredFluxMinu3D[idir][jdir].define(m_coveredSetsMinuG4[idir], m_ebisBox.getEBGraph(), numFlux);
              m_extendStatePlus3D[idir][jdir].setVal(0.);
              m_extendStateMinu3D[idir][jdir].setVal(0.);
              m_coveredFluxPlus3D[idir][jdir].setVal(0.);
              m_coveredFluxMinu3D[idir][jdir].setVal(0.);
            }
        }
    }
  IntVectSet        ivsIrregG4 =  m_ebisBox.getIrregIVS(m_validBoxG4);
  VoFIterator vofit(ivsIrregG4,   m_ebisBox.getEBGraph());
  m_irregVoFs = vofit.getVector();
  m_updateStencil.resize(m_irregVoFs.size());
  for (int ivof = 0; ivof<m_irregVoFs.size(); ivof++)
    {
      fillUpdateStencil(m_updateStencil[ivof], m_irregVoFs[ivof]);
    }

  if (m_useAgg)
    setSlopeStuff();

}
/******************/
void
EBPatchReactive::
fillUpdateStencil(EBPatchReactive::updateStencil_t& a_stencil, const VolIndex& a_vof)
{
  //i am going to use this assumption in many places
  int numVoFs = m_ebisBox.numVoFs(a_vof.gridIndex());
  if (numVoFs > 1)
    {
      a_stencil.m_vofOffset.m_multiValued = true;
      const BaseIVFAB<Real>& baseivfabPhi= m_primState.getMultiValuedFAB();
      a_stencil.m_vofOffset.m_offset = baseivfabPhi.getIndex(a_vof, 0) - baseivfabPhi.dataPtr(0);
    }
  else
    {
      a_stencil.m_vofOffset.m_multiValued = false;
      IntVect ncells  = m_validBoxG4.size();
      IntVect iv = a_vof.gridIndex() - m_validBoxG4.smallEnd();

      a_stencil.m_vofOffset.m_offset = iv[0] + iv[1]*ncells[0] ;
      if (SpaceDim==3)
        {
          a_stencil.m_vofOffset.m_offset +=    iv[2]*ncells[0]*ncells[1];
        }
    }
}
/****/
void
EBPatchReactive::
setSlopeStuff()
{
  //first get the elusive and irritating m_entirebox

  //this is the box sent to slope routines
  Box argBox[SpaceDim];
  getArgBox(argBox);

  //this is the raindance the slope routine uses to generate entirebox
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      Box box1, box2;
      Box loBox, hiBox, centerBox;
      int hasLo, hasHi;

      box1 = argBox[idir];
      box1.grow(idir, 1);
      box1 &= m_domain;

      box2 = argBox[idir];
      box2.grow(idir, 2);
      box2 &= m_domain;

      if (usesFourthOrderSlopes())
        {
          eblohicenter(loBox, hasLo, hiBox, hasHi, centerBox, m_entireBox[idir],
                       box2, m_domain, idir);
        }
      else
        {
          eblohicenter(loBox, hasLo, hiBox, hasHi, centerBox, m_entireBox[idir],
                       box1, m_domain, idir);
        }
    }

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      IntVectSet        ivsIrregG4 =  m_ebisBox.getIrregIVS(m_validBoxG4);
      ivsIrregG4 &= m_entireBox[idir] ;
      VoFIterator vofit(ivsIrregG4,   m_ebisBox.getEBGraph());
      const Vector<VolIndex>& vofs = vofit.getVector();
      m_slopVec[idir].resize(vofs.size());
      Vector<VoFStencil> slowStenLo(vofs.size());
      Vector<VoFStencil> slowStenHi(vofs.size());
      //for getting offsets
      EBCellFAB dumprim(m_ebisBox,      m_validBoxG4, numPrimitives());
      EBCellFAB dumslop(m_ebisBox, m_entireBox[idir], numSlopes());
      for (int ivof = 0; ivof < vofs.size(); ivof++)
        {
          const VolIndex& vof = vofs[ivof];
          const IntVect&  iv = vof.gridIndex();
          bool onLeftDomain = (iv[idir] == m_domain.domainBox().smallEnd(idir));
          bool onRighDomain = (iv[idir] == m_domain.domainBox().bigEnd(  idir));

          Vector<FaceIndex> facesLo = (m_ebisBox.getFaces(vof, idir, Side::Lo));
          Vector<FaceIndex> facesHi = (m_ebisBox.getFaces(vof, idir, Side::Hi));

          //the faces.size() == 1 thing should probably be done as >= 0
          //and we figure out what we want to do in the case of multiple faces.
          //as it is, I do not want to change the algorithm in the midst of trying
          //to do optimization.   Since pointgetslopes does it this way, this is the
          //way we shall do it.
          bool hasFacesLo = (facesLo.size()==1) && !onLeftDomain;
          bool hasFacesHi = (facesHi.size()==1) && !onRighDomain;

          m_slopVec[idir][ivof].hasLo = hasFacesLo;
          m_slopVec[idir][ivof].hasHi = hasFacesHi;
          m_slopVec[idir][ivof].slop_access.offset = dumslop.offset(  vof, 0);
          m_slopVec[idir][ivof].slop_access.dataID = dumslop.dataType(vof);
          if (hasFacesLo)
            {
              slowStenLo[ivof].add(facesLo[0].getVoF(Side::Lo), -1.0);
              slowStenLo[ivof].add(                        vof,  1.0);
            }
          if (hasFacesHi)
            {
              slowStenHi[ivof].add(facesHi[0].getVoF(Side::Hi),  1.0);
              slowStenHi[ivof].add(                        vof, -1.0);
            }
        }

      Vector<RefCountedPtr<BaseIndex   > > baseindice(vofs.size());
      Vector<RefCountedPtr<BaseStencil > > basestenlo(vofs.size());
      Vector<RefCountedPtr<BaseStencil > > basestenhi(vofs.size());
      for (int ivof= 0; ivof < vofs.size(); ivof++)
        {
          baseindice[ivof] =   RefCountedPtr<BaseIndex>(new   VolIndex(      vofs[ivof]));
          basestenlo[ivof] = RefCountedPtr<BaseStencil>(new VoFStencil(slowStenLo[ivof]));
          basestenhi[ivof] = RefCountedPtr<BaseStencil>(new VoFStencil(slowStenHi[ivof]));
        }

      m_slopStenLo[idir] = RefCountedPtr< AggStencil <EBCellFAB, EBCellFAB > >
        (new AggStencil <EBCellFAB, EBCellFAB >(baseindice, basestenlo, dumprim, dumslop));
      m_slopStenHi[idir] = RefCountedPtr< AggStencil <EBCellFAB, EBCellFAB > >
        (new AggStencil <EBCellFAB, EBCellFAB >(baseindice, basestenhi, dumprim, dumslop));

    }



}
/*****************************/
void
EBPatchReactive::
getArgBox(Box a_argBox[SpaceDim])
{
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (SpaceDim==2)
        {
          Box slopeBoxG2 = grow(m_validBox, 2);
          slopeBoxG2 &= m_domain;
          a_argBox[idir] = slopeBoxG2;
        }
      else
        {
          Box modBoxOpen = m_validBox;
          modBoxOpen.grow(3);
          modBoxOpen &= m_domain;
          a_argBox[idir] = modBoxOpen;
        }
    }
}
/*****************************/
Real
EBPatchReactive::
getMaxWaveSpeed(const EBCellFAB& a_consState,
                const Box& a_box)
{
  CH_TIME("EBPatchReactive::getMaxWaveSpeeed");
  CH_assert(m_isDefined && m_isBoxSet);
  Real speed = 0.0;
  const EBISBox& ebisBox = a_consState.getEBISBox();
  IntVectSet ivs(a_box);
  Real small, smallp, smallu, smallr;
  FORT_GETSMALL(CHF_REAL(small),
                CHF_REAL(smallp),
                CHF_REAL(smallu),
                CHF_REAL(smallr));
  //cannot just send to fortran because covered cell values
  //would get in and multiply-valued cells would use bogus
  //values that could also get in.  The could be avoided by
  //either masks or some clever manipulation of the underlying
  //basefab
  for (VoFIterator vofit(ivs, ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      //      const  IntVect& iv = vof.gridIndex();
      Real dense = Max(smallr, a_consState(vof, CRHO));
      Real eng   = Max(small,  a_consState(vof, CENG));
      Real vmax = 0.0;
      Real kinetic = 0.0;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          Real momen = a_consState(vof, CMOMX + idir);
          Real vel = Abs(momen/dense);
          kinetic += 0.5*vel*vel;
          vmax = Max(vmax, vel);
        }
      Real intEng = eng/dense - kinetic;
      intEng = Max(intEng, (Real)0.001*small);

      Vector<Real> MassFrac(m_nSpec);
      for (int ivar = 0; ivar < m_nSpec; ivar++)
        {
          MassFrac[ivar] = (a_consState(vof, CSPEC1+ivar))/dense;
        }

      Real cspeed;
      FORT_GETSOUNDSPEED(CHF_CONST_REAL(intEng),
                         CHF_CONST_VR(MassFrac),
                         CHF_REAL(cspeed));

      if (Abs(cspeed + vmax) > speed)
        {
          speed = cspeed+ vmax;
          if (speed > EBPatchReactive::getMaxWaveSpeed())
            {
              EBPatchReactive::setMaxWaveSpeed(speed);
              EBPatchReactive::setMaxWaveSpeedIV(vof.gridIndex());
             //s_maxWaveSpeed = speed;
             //s_maxWaveSpeedIV = vof.gridIndex();
            }
        }
    }
  return speed;

}
/******/
#ifdef CH_USE_HDF5

void EBPatchReactive::expressions(HDF5HeaderData& a_expressions)
{
  a_expressions.m_string["vector velocity"] = "momentum/<mass-density>";
  a_expressions.m_string["scalar kinetic_energy"] = "dot(velocity,velocity)/2";
  a_expressions.m_string["species mass fractions"] = "species density/<mass-density>";
  a_expressions.m_string["tempreature"] = "obtained from energ using numerical scheme";
  a_expressions.m_string["soundspeed"] = "sqrt((Cp/Cv)*Rgas*Temperature)";
  a_expressions.m_string["pressure"] = "<mass-density>*Rgas*Temperature";
}

#endif
/******/
void
EBPatchReactive::
consToPrim(EBCellFAB&       a_primState,
           const EBCellFAB& a_consState,
           const Box&       a_box,
           int              a_logflag,
           bool             a_verbose)
{
  CH_TIME("EBPatchReactive::consToPrim");
  CH_assert(m_isDefined && m_isBoxSet);
  CH_assert(m_isGammaSet);
  CH_assert(a_consState.getRegion().contains(a_box));
  //have to set the covered cell vals so need the cast
  //does not change real data
  //  EBCellFAB& conData = (EBCellFAB&) a_consState;
  //  setCoveredConsVals(conData);
  //debug
  if (a_verbose)
    {
      pout()  << "constoprim " << endl;
    }
  //end debug
  const BaseFab<Real>&   regCons = a_consState.getSingleValuedFAB();
  BaseFab<Real>&   regPrim = a_primState.getSingleValuedFAB();

  int iverbose = 0;
  if (a_verbose) iverbose = 1;

/*
  // debug
  if (a_logflag == 5)
    {
      pout() << "consState" << endl;
      FabDataOps::getFabData(regCons,a_consState.getEBISBox(),0);
    }
  // end debug
*/

  FORT_CONS2PRM(CHF_BOX(a_box),
                CHF_CONST_FRA(regCons),
                CHF_FRA(regPrim),
                CHF_CONST_INT(a_logflag),
                CHF_CONST_INT(iverbose));

  int nCons = numConserved();
  int nPrim = numPrimitives();
  IntVectSet ivsMulti = m_ebisBox.getMultiCells(a_box);
  for (VoFIterator vofit(ivsMulti, m_ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();

      Vector<Real> conserved(nCons);
      Vector<Real> primitive(nPrim);
      for (int ivar = 0; ivar < nCons; ivar++)
        {
          conserved[ivar] = a_consState(vof, ivar);
        }
       
      int sss = 0;
      FORT_POINTCONS2PRM(CHF_VR(conserved),
                         CHF_VR(primitive),
                         CHF_INT(sss),
                         CHF_CONST_INT(a_logflag));

      for (int ivar = 0; ivar < nPrim; ivar++)
        {
          a_primState(vof, ivar)  = primitive[ivar];
        }
      //  setCoveredPrimVals(a_primState);
    }
}
/******/
void
EBPatchReactive::
consToPrim(BaseIVFAB<Real>&  a_primState,
           const EBCellFAB&  a_consState,
           const IntVectSet& a_ivs)
{
  CH_TIME("EBPatchReactive::consToPrimIrrReg");
  CH_assert(isDefined());

  int nCons = numConserved();
  int nPrim = numPrimitives();
  for (VoFIterator vofit(a_ivs, m_ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();

      Vector<Real> conserved(nCons);
      Vector<Real> primitive(nPrim);
      for (int ivar = 0; ivar < nCons; ivar++)
        {
          conserved[ivar] = a_consState(vof, ivar);
        }

      int logflag = 0;
      int sss = 0;
      FORT_POINTCONS2PRM(CHF_VR(conserved),
                         CHF_VR(primitive),
                         CHF_INT(sss),
                         CHF_CONST_INT(logflag));

      for (int ivar = 0; ivar < nPrim; ivar++)
        {
          a_primState(vof, ivar)  = primitive[ivar];
        }
    }
}
/******/
void
EBPatchReactive::
consToPrim(BaseIVFAB<Real>&  a_primState,
           const BaseIVFAB<Real>&  a_consState,
           const IntVectSet& a_ivs)
{
  CH_TIME("EBPatchReactive::consToPrimIrr2");
  CH_assert(isDefined());
  int nCons = numConserved();
  int nPrim = numPrimitives();
  for (VoFIterator vofit(a_ivs, m_ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      Vector<Real> conserved(nCons);
      Vector<Real> primitive(nPrim);
      for (int ivar = 0; ivar < nCons; ivar++)
        {
          conserved[ivar] = a_consState(vof, ivar);
        }
      int logflag = 0;
      int sss = 0;
      FORT_POINTCONS2PRM(CHF_VR(conserved),
                         CHF_VR(primitive),
                         CHF_INT(sss),
                         CHF_CONST_INT(logflag));

      for (int ivar = 0; ivar < nPrim; ivar++)
        {
          a_primState(vof, ivar)  = primitive[ivar];
        }
    }
}
/******/
int
EBPatchReactive::
getCurLevel()
{
  return s_curLevel;
}
/*****/
void
EBPatchReactive::
setCurLevel(int a_curLevel)
{
  s_curLevel = a_curLevel;
}
/*****/
void
EBPatchReactive::
setCoveredConsVals(EBCellFAB& a_consState)
{
  CH_TIME("EBPatchReactive::setCoveredConsVals");
  a_consState.setInvalidData(0.0, CRHO);
  a_consState.setInvalidData(0.0, CMOMX);
  a_consState.setInvalidData(0.0, CMOMY);
#if CH_SPACEDIM==3
  a_consState.setInvalidData(0.0, CMOMZ);
#endif
  a_consState.setInvalidData(0E06, CENG);
  for (int ivar = 0; ivar < m_nSpec; ivar++)
   {
      a_consState.setInvalidData(4.0/m_nSpec, CSPEC1+ivar);
   }
}
/*****************************/
void EBPatchReactive::
computeFlattening(EBCellFAB&       a_flattening,
                  const EBCellFAB& a_primState,
                  const Box&       a_box)
{
  CH_TIME("EBPatchReactive::computeFlattening");
  CH_assert(isDefined());
  CH_assert(isDefined());
  CH_assert(a_primState.getRegion().contains(a_box));
  CH_assert(usesFourthOrderSlopes());
  CH_assert(a_flattening.nComp() == 1);
  CH_assert(a_flattening.getRegion().contains(a_box));

  BaseFab<Real>& regFlattening = a_flattening.getSingleValuedFAB();
  const BaseFab<Real>& regPrimState  = a_primState.getSingleValuedFAB();

  const Box& domainBox = m_domain.domainBox();
  EBCellFAB zetaDir(m_ebisBox, a_box, SpaceDim);
  EBCellFAB delta1U(m_ebisBox, a_box, SpaceDim);
  BaseFab<Real>& regZetaDir = zetaDir.getSingleValuedFAB();
  BaseFab<Real>& regDelta1U = delta1U.getSingleValuedFAB();

  Interval velInterval= velocityInterval();
  int v0index = velInterval.begin();

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      Box box1 = a_box;
      box1.grow(idir, 1);
      box1 &= m_domain;

      Box box2 = a_box;
      box2.grow(idir, 2);
      box2 &= m_domain;

      Box box3 = a_box;
      box3.grow(idir, 3);
      box3 &= m_domain;

      CH_assert(a_primState.getRegion().contains(box3));

      Box loBox, hiBox, centerBox, entireBox;
      int hasLo, hasHi;

      eblohicenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
                 box3, m_domain, idir);

      EBCellFAB delta1P(m_ebisBox, entireBox, 1);
      BaseFab<Real>& regDelta1P = delta1P.getSingleValuedFAB();

      int pressIndex = pressureIndex();
      //compute delta1P

      FORT_GETGRAD(CHF_FRA1(regDelta1P, 0),
                   CHF_CONST_FRA1(regPrimState, pressIndex),
                   CHF_CONST_INT(idir),
                   CHF_BOX(loBox),
                   CHF_CONST_INT(hasLo),
                   CHF_BOX(hiBox),
                   CHF_CONST_INT(hasHi),
                   CHF_BOX(centerBox));

      //update the irregular vofsn
      for (int ivof = 0; ivof<m_irregVoFs.size(); ivof++)
        {
          const VolIndex& vof = m_irregVoFs[ivof];
          if (entireBox.contains(vof.gridIndex()))
            {
              const IntVect&  iv = vof.gridIndex();
              //one-sided diffs on domain bndry
              bool onLeftDomain = iv[idir] == domainBox.smallEnd(idir);
              bool onRighDomain = iv[idir] == domainBox.bigEnd(idir);
              bool hasFacesLeft = (m_ebisBox.numFaces(vof, idir, Side::Lo) > 0) && !onLeftDomain;
              bool hasFacesRigh = (m_ebisBox.numFaces(vof, idir, Side::Hi) > 0) && !onRighDomain;

              Real valCent = a_primState(vof, pressIndex);
              Real dpl = 0;
              Real dpr = 0;
              Real dpc = 0;
              //compute one-sided diffs where you have them
              if (hasFacesLeft)
                {
                  Vector<FaceIndex> facesLeft =
                    m_ebisBox.getFaces(vof, idir, Side::Lo);
                  Real valLeft = 0.0;
                  for (int iface = 0; iface <facesLeft.size(); iface++)
                    {
                      VolIndex vofLeft = facesLeft[iface].getVoF(Side::Lo);
                      valLeft += a_primState(vofLeft, pressIndex);
                    }
                  valLeft /= Real(facesLeft.size());
                  dpl = valCent - valLeft;
                }
              if (hasFacesRigh)
                {
                  Vector<FaceIndex> facesRigh =
                    m_ebisBox.getFaces(vof, idir, Side::Hi);
                  Real valRigh = 0.0;
                  for (int iface = 0; iface <facesRigh.size(); iface++)
                    {
                      VolIndex vofRigh = facesRigh[iface].getVoF(Side::Hi);
                      valRigh += a_primState(vofRigh, pressIndex);
                    }
                  valRigh /= Real(facesRigh.size());
                  dpr = valRigh - valCent;
                }
              if (hasFacesLeft && hasFacesRigh)
                {
                  dpc = 0.5*(dpl+dpr);
                }
              else if (!hasFacesLeft && !hasFacesRigh)
                {
                  dpc = 0.0;
                }
              else if (hasFacesLeft && !hasFacesRigh)
                {
                  dpc = dpl;
                }
              else if (hasFacesRigh && !hasFacesLeft)
                {
                  dpc = dpr;
                }

              delta1P(vof, 0) = dpc;
            }
        }

      eblohicenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
                   box2, m_domain, idir);

      EBCellFAB delta2P(m_ebisBox, entireBox, 1);
      EBCellFAB bulkMin(m_ebisBox, entireBox, 1);
      EBCellFAB zetaTwi(m_ebisBox, entireBox, 1);
      BaseFab<Real>& regDelta2P = delta2P.getSingleValuedFAB();
      BaseFab<Real>& regBulkMin = bulkMin.getSingleValuedFAB();
      BaseFab<Real>& regZetaTwi = zetaTwi.getSingleValuedFAB();

      //compute delta2P

      FORT_GETDPTWO(CHF_FRA1(regDelta2P, 0),
                    CHF_CONST_FRA1(regDelta1P, 0),
                    CHF_CONST_INT(idir),
                    CHF_BOX(loBox),
                    CHF_CONST_INT(hasLo),
                    CHF_BOX(hiBox),
                    CHF_CONST_INT(hasHi),
                    CHF_BOX(centerBox));

      //compute min3(press, d) = bulkMin
      int bulkIndex = bulkModulusIndex();

      FORT_MIN3PTS(CHF_FRA1(regBulkMin, 0),
                   CHF_CONST_FRA1(regPrimState, bulkIndex),
                   CHF_CONST_INT(idir),
                   CHF_BOX(loBox),
                   CHF_CONST_INT(hasLo),
                   CHF_BOX(hiBox),
                   CHF_CONST_INT(hasHi),
                   CHF_BOX(centerBox));

      //compute preliminary flattening coeff

      FORT_GETFLAT(CHF_FRA1(regZetaTwi, 0),
                   CHF_CONST_FRA1(regDelta1P, 0),
                   CHF_CONST_FRA1(regDelta2P, 0),
                   CHF_CONST_FRA1(regBulkMin, 0),
                   CHF_BOX(entireBox));

      //update the irregular vofsn
      for (int ivof = 0; ivof<m_irregVoFs.size(); ivof++)
        {
          const VolIndex& vof = m_irregVoFs[ivof];
          if (entireBox.contains(vof.gridIndex()))
            {
              const IntVect&  iv = vof.gridIndex();
              //one-sided diffs on domain bndry
              bool onLeftDomain = iv[idir] == domainBox.smallEnd(idir);
              bool onRighDomain = iv[idir] == domainBox.bigEnd(idir);
              bool hasFacesLeft = (m_ebisBox.numFaces(vof, idir, Side::Lo) > 0) && !onLeftDomain;
              bool hasFacesRigh = (m_ebisBox.numFaces(vof, idir, Side::Hi) > 0) && !onRighDomain;

              Real dp2c;
              Real bulkMinL, bulkMinR, bulkMinC;

              Real dp1Cent = delta1P(vof, 0);
              Real bulkCent = a_primState(vof, bulkIndex);
              Real dp1Left=0;
              Real dp1Righ=0;
              Real bulkLeft=0;
              Real bulkRigh=0;
              //compute one-sided diffs where you have them
              if (hasFacesLeft)
                {
                  Vector<FaceIndex> facesLeft =
                    m_ebisBox.getFaces(vof, idir, Side::Lo);
                  dp1Left = 0.0;
                  bulkLeft = 0.0;
                  for (int iface = 0; iface <facesLeft.size(); iface++)
                    {
                      VolIndex vofLeft = facesLeft[iface].getVoF(Side::Lo);
                      dp1Left += delta1P(vofLeft, 0);
                      bulkLeft +=a_primState(vofLeft, bulkIndex);
                    }
                  dp1Left /= Real(facesLeft.size());
                  bulkLeft /= Real(facesLeft.size());

                  bulkMinL = Min(bulkCent, bulkLeft);
                }
              if (hasFacesRigh)
                {
                  Vector<FaceIndex> facesRigh =
                    m_ebisBox.getFaces(vof, idir, Side::Hi);
                  dp1Righ = 0.0;
                  bulkRigh = 0.0;
                  for (int iface = 0; iface <facesRigh.size(); iface++)
                    {
                      VolIndex vofRigh = facesRigh[iface].getVoF(Side::Hi);
                      dp1Righ += delta1P(vofRigh, 0);
                      bulkRigh += a_primState(vofRigh, bulkIndex);
                    }
                  dp1Righ /= Real(facesRigh.size());
                  bulkRigh /= Real(facesRigh.size());
                  bulkMinR = Min(bulkCent, bulkRigh);
                }

              if (hasFacesLeft && hasFacesRigh)
                {
                  dp2c = (dp1Left+dp1Righ);
                  bulkMinC = Min(bulkMinR, bulkMinL);
                }
              else if (!hasFacesLeft && !hasFacesRigh)
                {
                  dp2c = dp1Cent;
                  bulkMinC = bulkCent;
                }
              else if (hasFacesLeft && !hasFacesRigh)
                {
                  dp2c = dp1Cent + dp1Left;
                  bulkMinC = bulkMinL;
                }
              else if (hasFacesRigh && !hasFacesLeft)
                {
                  dp2c = dp1Cent + dp1Righ;
                  bulkMinC = bulkMinR;
                }

              delta2P(vof, 0) = dp2c;
              bulkMin(vof, 0) = bulkMinC;

              //this is that whole zeta(dp1, dp2, p0) func
              Real r0 = 0.75;  Real r1= 0.85; Real d = 0.33;
              Real d1pVoF = Abs(delta1P(vof, 0)) ;
              //
              //     bad idea among many
              Real smallp = 1.0e-2;
              Real d2pVoF = Max(Abs(delta2P(vof, 0)), smallp);

              Real strength = Abs(d1pVoF/bulkMinC);
              Real zetaFuncDP;
              if (strength >= d)
                {
                  Real ratio =  d1pVoF/d2pVoF;
                  if ( ratio <= r0)
                    {
                      zetaFuncDP = 1.0;
                    }
                  else if (ratio >= r1)
                    {
                      zetaFuncDP = 0.0;
                    }
                  else
                    {
                      zetaFuncDP= 1.0 - (ratio - r0)/(r1-r0);;
                    }
                }
              else //strength < d
                {
                  zetaFuncDP = 1.0;
                }
              zetaTwi(vof, 0) = zetaFuncDP;
            }
        }

      eblohicenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
                   box1, m_domain, idir);

      //compute zeta = min3(zetatwid) and delta1U

      FORT_MIN3PTS(CHF_FRA1(regZetaDir, idir),
                   CHF_CONST_FRA1(regZetaTwi, 0),
                   CHF_CONST_INT(idir),
                   CHF_BOX(loBox),
                   CHF_CONST_INT(hasLo),
                   CHF_BOX(hiBox),
                   CHF_CONST_INT(hasHi),
                   CHF_BOX(centerBox));

      int velIndex = v0index + idir;

      FORT_GETGRAD(CHF_FRA1(regDelta1U, idir),
                   CHF_CONST_FRA1(regPrimState, velIndex),
                   CHF_CONST_INT(idir),
                   CHF_BOX(loBox),
                   CHF_CONST_INT(hasLo),
                   CHF_BOX(hiBox),
                   CHF_CONST_INT(hasHi),
                   CHF_BOX(centerBox));

      for (int ivof = 0; ivof<m_irregVoFs.size(); ivof++)
        {
          const VolIndex& vof = m_irregVoFs[ivof];
          if (entireBox.contains(vof.gridIndex()))
            {
              const IntVect&  iv = vof.gridIndex();
              //one-sided diffs on domain bndry
              bool onLeftDomain = iv[idir] == domainBox.smallEnd(idir);
              bool onRighDomain = iv[idir] == domainBox.bigEnd(idir);
              bool hasFacesLeft = (m_ebisBox.numFaces(vof, idir, Side::Lo) > 0) && !onLeftDomain;
              bool hasFacesRigh = (m_ebisBox.numFaces(vof, idir, Side::Hi) > 0) && !onRighDomain;

              Real velLeft=0;
              Real velRigh=0;
              Real du1c=0;
              Real zetaMinL, zetaMinR, zetaMinC;
              Real zetaLeft, zetaRigh;
              Real velCent = a_primState(vof, velIndex);
              Real zetaCent = zetaTwi(vof, 0);
              //compute one-sided diffs where you have them
              if (hasFacesLeft)
                {
                  Vector<FaceIndex> facesLeft =
                    m_ebisBox.getFaces(vof, idir, Side::Lo);
                  velLeft  = 0.0;
                  zetaLeft = 0.0;
                  for (int iface = 0; iface <facesLeft.size(); iface++)
                    {
                      VolIndex vofLeft = facesLeft[iface].getVoF(Side::Lo);
                      velLeft += a_primState(vofLeft, velIndex);
                      zetaLeft +=    zetaTwi(vofLeft, 0);
                    }
                  velLeft /= Real(facesLeft.size());
                  zetaLeft /= Real(facesLeft.size());

                  zetaMinL = Min(zetaCent, zetaLeft);
                }
              if (hasFacesRigh)
                {
                  Vector<FaceIndex> facesRigh =
                    m_ebisBox.getFaces(vof, idir, Side::Hi);
                  velRigh = 0.0;
                  zetaRigh = 0.0;
                  for (int iface = 0; iface <facesRigh.size(); iface++)
                    {
                      VolIndex vofRigh = facesRigh[iface].getVoF(Side::Hi);
                      velRigh += a_primState(vofRigh, velIndex);
                      zetaRigh += zetaTwi(vofRigh, 0);
                    }
                  velRigh /= Real(facesRigh.size());
                  zetaRigh /= Real(facesRigh.size());
                  zetaMinR = Min(zetaCent, zetaRigh);
                }

              if (hasFacesLeft && hasFacesRigh)
                {
                  du1c = 0.5*(velRigh - velLeft);
                  zetaMinC = Min(zetaMinR, zetaMinL);
                }
              else if (!hasFacesLeft && !hasFacesRigh)
                {
                  du1c = 0.0;
                  zetaMinC = zetaCent;
                }
              else if (hasFacesLeft && !hasFacesRigh)
                {
                  du1c = velCent-velLeft;
                  zetaMinC = zetaMinL;
                }
              else if (hasFacesRigh && !hasFacesLeft)
                {
                  du1c = velRigh-velCent;
                  zetaMinC = zetaMinR;
                }

              delta1U(vof, idir) = du1c;
              zetaDir(vof, idir) = zetaMinC;
            }

        }
    }

  // take the minimum of all directions if the sum of velocity diffs
  // are negative.  unity (no flattening) otherwise.

  FORT_MINFLAT(CHF_FRA1(regFlattening, 0),
               CHF_CONST_FRA(regZetaDir),
               CHF_CONST_FRA(regDelta1U),
               CHF_BOX(a_box));

  for (int ivof = 0; ivof<m_irregVoFs.size(); ivof++)
    {
      const VolIndex& vof = m_irregVoFs[ivof];
      if (a_box.contains(vof.gridIndex()))
        {
          Real velDiffSum = 0.0;
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              velDiffSum += delta1U(vof, idir);
            }
          Real finalFlat;
          if (velDiffSum < 0)
            {
              finalFlat = zetaDir(vof, 0);
              for (int idir = 1; idir < SpaceDim; idir++)
                {
                  finalFlat = Min(finalFlat, zetaDir(vof, idir));
                }
              CH_assert(finalFlat >= 0.0);
            }
          else
            {
              finalFlat = 1.0;
            }
          a_flattening(vof, 0) = finalFlat;
        }
    }
}
/*****************************/
void
EBPatchReactive::
setSource(EBCellFAB&       a_source,
          const EBCellFAB& a_consState,
          const Box&       a_box)
{
  if (m_doRZCoords)
    {
      setRZSource(a_source,
                  a_consState,
                  a_box);
    }
}
/*****************************/
void
EBPatchReactive::
setRZSource(EBCellFAB&       a_source,
            const EBCellFAB& a_consState,
            const Box&       a_box)
{
  a_source.setVal(0.);
  const BaseFab<Real>& regConsState = a_consState.getSingleValuedFAB();
  BaseFab<Real>& regSource          =    a_source.getSingleValuedFAB();
  FORT_SETSOURCERZ(CHF_FRA(regSource),
                   CHF_CONST_FRA(regConsState),
                   CHF_CONST_REAL(m_dx[0]),
                   CHF_BOX(a_box));

  //pointwise operation so just have to do the multi valued cells
  //as an irregular iteration
  IntVectSet multiIVS = m_ebisBox.getMultiCells(a_box);
  for (VoFIterator vofit(multiIVS, m_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      const IntVect& iv = vof.gridIndex();
      Real dense = a_consState(vof, CRHO);
      Real energy = a_consState(vof, CENG);
      RealVect momen;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          momen[idir] = a_consState(vof, CMOMX + idir);
        }
      Real radius = (Real(iv[0]) + 0.5)*m_dx[0];
      Real densitySource, pressureSource;
      FORT_POINTSETSOURCERZ(CHF_REAL(densitySource),
                            CHF_REAL(pressureSource),
                            CHF_CONST_REAL(dense),
                            CHF_CONST_REALVECT(momen),
                            CHF_CONST_REAL(energy),
                            CHF_CONST_REAL(radius));
      //only the density and energy equations have non-zero source
      //so set the others to zero
      for (int ivar = 0; ivar < numConserved(); ivar++)
        {
          a_source(vof, ivar) =  0.0;
        }
      a_source(vof, QRHO) = densitySource;
      a_source(vof, QPRES) = pressureSource;
    }
}
/*****************************/
void EBPatchReactive::
primitivesAndDivergences(EBCellFAB&          a_nonConsDivF,
                         EBCellFAB&          a_consState,
                         BaseIVFAB<Real>     a_coveredPrimMinu[SpaceDim],
                         BaseIVFAB<Real>     a_coveredPrimPlus[SpaceDim],
                         Vector<VolIndex>    a_coveredFaceMinu[SpaceDim],
                         Vector<VolIndex>    a_coveredFacePlus[SpaceDim],
                         EBFluxFAB&          a_flux,
                         BaseIVFAB<Real>&    a_ebIrregFlux,
                         BaseIVFAB<Real>&    a_nonConservativeDivergence,
                         const EBCellFAB&    a_flattening,
                         const EBCellFAB&    a_source,
                         const Box&          a_box,
                         const IntVectSet&   a_ivsSmall,
                         const DataIndex&    a_dit,
                         bool                a_verbose)
{
  CH_TIME("EBPatchPolytropic::regularUpdate");
    CH_assert(isDefined());
  int numCons = numConserved();

  CH_assert(a_flux.nComp() >= numConserved());
  CH_assert(a_consState.nComp() >= numConserved());
  setCoveredConsVals(a_consState);

  EBCellFAB slopePrim[SpaceDim];
  EBCellFAB slopeNLim[SpaceDim];

  computeFluxes(a_flux,
                m_coveredFluxMinuG4, m_coveredFluxPlusG4,
                a_coveredFaceMinu,   a_coveredFacePlus,
                m_primState,    slopePrim, slopeNLim,
                a_flattening, a_consState, a_source,
                a_box, a_dit,a_verbose);

  //now I happen to know that m_extendedStateMinuG4 holds  coveredPrimMinu and so on
  Interval inter(0, numPrimitives()-1);

  for(int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      IntVectSet ivsMinu = m_coveredSetsMinuG4[faceDir] & a_box;
      IntVectSet ivsPlus = m_coveredSetsPlusG4[faceDir] & a_box;
      a_coveredPrimMinu[faceDir].define(ivsMinu, m_ebisBox.getEBGraph(), numPrimitives());
      a_coveredPrimPlus[faceDir].define(ivsPlus, m_ebisBox.getEBGraph(), numPrimitives());
      a_coveredPrimMinu[faceDir].copy(a_box, inter, a_box, m_extendStateMinuG4[faceDir], inter);
      a_coveredPrimPlus[faceDir].copy(a_box, inter, a_box, m_extendStatePlusG4[faceDir], inter);
    }

  //now that we have the fluxes, modify them with
  //artificial viscosity if that is called for.
  //the artificial viscosity correction to the
  //embedded boundary flux is computed inside
  //computeebirregflux
  if(usesArtificialViscosity())
    {
      EBFluxFAB openDivU(m_ebisBox, a_box, 1);

      getFaceDivergence(openDivU,
                        m_primState, slopeNLim,
                        a_box, a_ivsSmall);

      applyArtificialViscosity(a_flux,
                               m_coveredFluxMinuG4,
                               m_coveredFluxPlusG4,
                               m_coveredFaceMinuG4,
                               m_coveredFacePlusG4,
                               a_consState,
                               openDivU,
                               a_box,
                               a_ivsSmall);
    }

  //this is the stable, non-conservative estimate of the solution update
  nonconservativeDivergence(a_nonConsDivF, a_flux,
                            m_coveredFluxMinuG4,
                            m_coveredFluxPlusG4,
                            a_coveredFaceMinu,
                            a_coveredFacePlus,
                            a_box);

  //copy nonconservative div f into sparse output thingy
  IntVectSet ivsIrreg = a_ivsSmall;
  for(VoFIterator vofit(ivsIrreg, m_ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      for(int ivar = 0; ivar < numCons; ivar++)
        {
          a_nonConservativeDivergence(vof, ivar) = a_nonConsDivF(vof, ivar);
        }
    }

  //compute irregular boundary flux.  this includes
  //an artificial viscosity modification if appropriate
  computeEBIrregFlux( a_ebIrregFlux, m_primState,
                      slopeNLim, ivsIrreg, a_source);
}
/*****************************/
void
EBPatchReactive::
computeFluxes(EBFluxFAB&       a_flux,
              BaseIVFAB<Real>  a_coveredFluxMinu[SpaceDim],
              BaseIVFAB<Real>  a_coveredFluxPlus[SpaceDim],
              Vector<VolIndex> a_coveredFaceMinu[SpaceDim],
              Vector<VolIndex> a_coveredFacePlus[SpaceDim],
              EBCellFAB&       a_primState,
              EBCellFAB        a_slopePrim[SpaceDim],
              EBCellFAB        a_slopeNLim[SpaceDim],
              const EBCellFAB& a_flattening,
              const EBCellFAB& a_consState,
              const EBCellFAB& a_source,
              const Box&       a_box,
              const DataIndex& a_dit,
              bool             a_verbose)
{

  if (SpaceDim==2)
    {
      extrapolatePrim2D(m_primMinu, m_primPlus,
                        a_primState,    a_slopePrim, a_slopeNLim,
                        a_flattening, a_consState, a_source,
                        a_box, a_dit,a_verbose);
    }
  else if (SpaceDim==3)
    {

      extrapolatePrim3D(m_primMinu,   m_primPlus,
                        a_primState,  a_slopePrim, a_slopeNLim,
                        a_flattening, a_consState, a_source,
                        a_box, a_dit,a_verbose);
    }
  else
    {
      MayDay::Error("bogus SpaceDim");
    }

  IntVectSet        coveredSetsPlus[SpaceDim];
  IntVectSet        coveredSetsMinu[SpaceDim];
  Box faceBox[SpaceDim];
  Box bndryFaceBox[SpaceDim];
  //this keeps the fluxes from being calculated
  //on boundaries
  for (int dir1 = 0; dir1 < SpaceDim; ++dir1)
    {
      bndryFaceBox[dir1] = a_box;
      bndryFaceBox[dir1] &= m_domain;
      bndryFaceBox[dir1].surroundingNodes(dir1);

      faceBox[dir1] = a_box;
      faceBox[dir1].grow(dir1,1);
      faceBox[dir1] &= m_domain;
      faceBox[dir1].grow(dir1,-1);
      faceBox[dir1].surroundingNodes(dir1);
    }
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      IntVectSet irregIVSPlus, irregIVSMinu;
      computeCoveredFaces(a_coveredFacePlus[idir],
                          coveredSetsPlus[idir],
                          irregIVSPlus,
                          idir, Side::Hi, a_box);
      computeCoveredFaces(a_coveredFaceMinu[idir],
                          coveredSetsMinu[idir],
                          irregIVSMinu,
                          idir, Side::Lo, a_box);
    }
  for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      //extrapolate to covered faces using updated extrapolated state
      extrapToCoveredFaces(m_extendStateMinuG4[faceDir],
                           m_primMinu[faceDir],
                           m_primPlus[faceDir],
                           a_primState,
                           a_coveredFaceMinu[faceDir],
                           faceDir, Side::Lo, a_box);

      extrapToCoveredFaces(m_extendStatePlusG4[faceDir],
                           m_primMinu[faceDir],
                           m_primPlus[faceDir],
                           a_primState,
                           a_coveredFacePlus[faceDir],
                           faceDir, Side::Hi, a_box);

      //equation 1.19 top
      riemann(a_flux[faceDir], m_primPlus[faceDir], m_primMinu[faceDir],
              faceDir, faceBox[faceDir]);

      // Use the user supplied PhysBC object to obtain boundary fluxes
      m_bc->fluxBC(a_flux, a_primState, m_primMinu[faceDir],
                   Side::Lo,  m_time, m_ebisBox, a_dit,a_box, bndryFaceBox[faceDir], faceDir);
      m_bc->fluxBC(a_flux, a_primState, m_primPlus[faceDir],
                   Side::Hi,  m_time, m_ebisBox, a_dit,a_box, bndryFaceBox[faceDir], faceDir);

      //solve riemann problem between extended state and
      //value in vof to get covered flux.
      //Equation 1.19 bottom.
      riemann(a_coveredFluxMinu[faceDir],
              m_extendStateMinuG4[faceDir], m_primMinu[faceDir],
              a_coveredFaceMinu[faceDir], faceDir, Side::Lo, a_box);
      riemann(a_coveredFluxPlus[faceDir],
              m_extendStatePlusG4[faceDir], m_primPlus[faceDir],
              a_coveredFacePlus[faceDir], faceDir, Side::Hi, a_box);
    }
}
/******/
void
EBPatchReactive::
extrapolatePrim2D(EBCellFAB           a_primMinu[SpaceDim],
                  EBCellFAB           a_primPlus[SpaceDim],
                  EBCellFAB&          a_primState,
                  EBCellFAB           a_slopesPrim[SpaceDim],
                  EBCellFAB           a_slopesSeco[SpaceDim],
                  const EBCellFAB&    a_flattening,
                  const EBCellFAB&    a_consState,
                  const EBCellFAB&    a_source,
                  const Box&          a_box,
                  const DataIndex&    a_dit,
                  bool                a_verbose)
{
  CH_TIME("EBPatchReactive::extrapolatePrim2D");

  //now define the plethora of data holders that I need

  //now do the actual computation.
  //1. transform to primitive state
  int logflag = 0;
  Box primBox = a_consState.box();

  consToPrim(a_primState, a_consState, primBox, logflag); // debug change logflag

  doNormalDerivativeExtr2D(a_primMinu,
                           a_primPlus,
                           m_fluxOne,
                           m_coveredFluxNormMinu,
                           m_coveredFluxNormPlus,
                           m_coveredFaceMinuG4,
                           m_coveredFacePlusG4,
                           a_slopesPrim,
                           a_slopesSeco,
                           a_flattening,
                           a_primState,
                           a_source,
                           a_dit,
                           a_box);


  /**/
  // Do the final corrections to the fluxes
  //this only happens on the non-ghosted box
  finalExtrap2D(a_primMinu,
                a_primPlus,
                m_coveredFluxNormMinu,
                m_coveredFluxNormPlus,
                m_coveredFaceMinuG4,
                m_coveredFacePlusG4,
                m_fluxOne,
                a_primState,
                a_slopesPrim,
                a_slopesSeco,
                a_box
                );

  /**/



}                     
/*******************/
void EBPatchReactive::
doNormalDerivativeExtr2D(EBCellFAB              a_primMinu[SpaceDim],
                         EBCellFAB              a_primPlus[SpaceDim],
                         EBFaceFAB              a_fluxOne[SpaceDim],
                         BaseIVFAB<Real>        a_coveredFluxNormMinu[SpaceDim],
                         BaseIVFAB<Real>        a_coveredFluxNormPlus[SpaceDim],
                         Vector<VolIndex>       a_coveredFaceNormMinu[SpaceDim],
                         Vector<VolIndex>       a_coveredFaceNormPlus[SpaceDim],
                         EBCellFAB              a_slopesPrim[SpaceDim],
                         EBCellFAB              a_slopesSeco[SpaceDim],
                         const EBCellFAB&       a_flattening,
                         const EBCellFAB&       a_primState,
                         const EBCellFAB&       a_source,
                         const DataIndex&       a_dit,
                         const Box&             a_box)

{
  CH_TIME("EBPatchReactive::doNormalDerivative2D");
  Box slopeBoxG1 = grow(a_box, 1);
  Box slopeBoxG2 = grow(a_box, 2);
  slopeBoxG1 &= m_domain;
  slopeBoxG2 &= m_domain;
  Box cellBoxG2[SpaceDim];
  Box cellBoxG1[SpaceDim];
  Box faceBox[SpaceDim];
  Box bndryFaceBox[SpaceDim];
  int numSlop = numSlopes();
  for (int idir = 0; idir < SpaceDim; ++idir)
    {

      cellBoxG2[idir] = a_box;
      cellBoxG2[idir].grow(2);
      cellBoxG2[idir].grow(idir, -1);
      cellBoxG2[idir] &= m_domain;

      cellBoxG1[idir] = a_box;
      cellBoxG1[idir].grow(1);
      cellBoxG1[idir].grow(idir, -1);
      cellBoxG1[idir] &= m_domain;

      bndryFaceBox[idir] = a_box;
      bndryFaceBox[idir].grow(1);
      bndryFaceBox[idir] &= m_domain;
      bndryFaceBox[idir].surroundingNodes(idir);


      faceBox[idir] = a_box;
      faceBox[idir].grow(2);
      faceBox[idir] &= m_domain;
      faceBox[idir].grow(idir, -1);
      faceBox[idir].surroundingNodes(idir);
    }

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      a_slopesPrim[idir].define(m_ebisBox, slopeBoxG2, numSlop);
      a_slopesSeco[idir].define(m_ebisBox, slopeBoxG2, numSlop);

      m_extendStateNormPlus[idir].setVal(0.);
      m_extendStateNormMinu[idir].setVal(0.);
      a_coveredFluxNormPlus[idir].setVal(0.);
      a_coveredFluxNormMinu[idir].setVal(0.);
      a_primPlus[idir].setVal(0.);
      a_primMinu[idir].setVal(0.);
      a_fluxOne[idir].setVal(0.);
    }

  //check flattening coefficients have correct box
  if (usesFlattening())
    {
      if (!a_flattening.getRegion().contains(slopeBoxG2))
        {
          MayDay::Error("flattening defined over too small a region");
        }
    }

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      slope(a_slopesPrim[idir], a_slopesSeco[idir], a_primState,
            a_flattening, idir, slopeBoxG2, m_useAgg);
    }

  for (int idir = 0; idir < SpaceDim; idir++)
    {

      //3. Compute the effect of the normal derivative terms and the
      //source term on the extrapolation in space and time from cell
      //centers to faces. Equation 1.7 (top).
      normalPred(a_primMinu[idir], a_primPlus[idir],
                 a_primState, a_slopesPrim[idir],
                 m_dt/m_dx[idir], idir, slopeBoxG2);

      // If the source term is valid add it to the primitive quantities
      //if source term is not defined, assume it is zero.  lots of
      //problems will have no sources.  not the most elegant method
      //but it beats passing a null pointer.
      if (a_source.isDefined())
        {
          if (!s_conservativeSource)
            {
              incrementWithSource(a_primMinu[idir], a_source,
                                  0.5*m_dt, slopeBoxG2);
              incrementWithSource(a_primPlus[idir], a_source,
                                  0.5*m_dt, slopeBoxG2);
            }
          else
            {
              int logflag = 0;
              const Box& modBox = slopeBoxG2;
              EBCellFAB  consTemp(m_ebisBox, modBox, numConserved());

              primToCons(consTemp,  a_primMinu[idir],  modBox);
              
/*              // DEBUG:
              pout() << "timestep size, dt = " << m_dt << endl;
              pout() <<"consTemp:" << endl;
              FabDataOps::getFabData(consTemp,true);
              pout() << "source:" << endl;
              FabDataOps::getFabData(a_source, true);
              // end DEBUG
*/              
              incrementWithSource(consTemp,  a_source, 0.5*m_dt,  modBox); //
              consToPrim(a_primMinu[idir],  consTemp,    modBox, logflag);

              primToCons(consTemp,  a_primPlus[idir],  modBox);
              incrementWithSource(consTemp,  a_source,  0.5*m_dt, modBox); //
              consToPrim(a_primPlus[idir],  consTemp,   modBox, logflag);
            }
        }

      //4. compute estimates of the flux suitable for computing
      // 1D flux derivatives using a riemann solver for the interior R
      //and for the boundary RB.  Equation 1.8
      riemann(a_fluxOne[idir], a_primPlus[idir], a_primMinu[idir],
              idir, faceBox[idir]);

      //some wackiness to work around the fact that fluxfab did
      //not exist when a lot of this was written.
      Vector<EBFaceFAB*> bcfluxes(SpaceDim);
      for (int kdir = 0; kdir < SpaceDim; kdir++)
        {
          bcfluxes[kdir] = &(a_fluxOne[kdir]);
        }
      EBFluxFAB fluxAlias;
      fluxAlias.alias(bcfluxes);

      // Use the user supplied PhysBC object to obtain boundary fluxes
      m_bc->fluxBC(fluxAlias, a_primState, a_primMinu[idir],
                   Side::Lo, m_time, m_ebisBox, a_dit,a_box, bndryFaceBox[idir], idir);
      m_bc->fluxBC(fluxAlias, a_primState, a_primPlus[idir],
                   Side::Hi, m_time, m_ebisBox, a_dit,a_box, bndryFaceBox[idir], idir);

      //extrapolate to covered faces
      extrapToCoveredFaces(m_extendStateNormMinu[idir],
                           a_primMinu[idir],
                           a_primPlus[idir],
                           a_primState,
                           a_coveredFaceNormMinu[idir],
                           idir, Side::Lo, slopeBoxG1);

      extrapToCoveredFaces(m_extendStateNormPlus[idir],
                           a_primMinu[idir],
                           a_primPlus[idir],
                           a_primState,
                           a_coveredFaceNormPlus[idir],
                           idir, Side::Hi, slopeBoxG1);

      //solve riemann problem between extended state and
      //value in vof to get covered flux  Equation 1.9
      riemann(a_coveredFluxNormMinu[idir],
              m_extendStateNormMinu[idir], a_primMinu[idir],
              a_coveredFaceNormMinu[idir], idir, Side::Lo, slopeBoxG2);
      riemann(a_coveredFluxNormPlus[idir],
              m_extendStateNormPlus[idir], a_primPlus[idir],
              a_coveredFaceNormPlus[idir], idir, Side::Hi, slopeBoxG2);
    }
}
/*****************************/
void
EBPatchReactive::
slope(EBCellFAB&       a_slopePrim,
      EBCellFAB&       a_slopeNLim,
      const EBCellFAB& a_primState,
      const EBCellFAB& a_flattening,
      const int&       a_dir,
      const Box&       a_box,
      bool a_doAggregated)
{
  CH_TIME("EBPatchReactive::slope");
  int numSlope = numSlopes();
  CH_assert(a_slopePrim.nComp() == numSlopes());
  CH_assert(a_primState.nComp() >= numSlopes());
  EBCellFAB delta2W, deltaWL, deltaWR, deltaWC;

  delta2W.setVal(0.);
  deltaWC.setVal(0.);

  {
    if (!usesZeroSlopes())
     {
       CH_TIME("second_order_slopes");
       doSecondOrderSlopes(delta2W,
                           deltaWL,
                           deltaWR,
                           deltaWC,
                           a_primState,
                           a_dir,
                           a_box,
                           a_doAggregated);


       a_slopeNLim.copy(deltaWC);

       if (usesFourthOrderSlopes())
        {
          EBCellFAB delta4W;
          {
            CH_TIME("forth_order_slopes");
            doFourthOrderSlopes(delta4W,
                                deltaWC,
                                delta2W,
                                deltaWL,
                                deltaWR,
                                a_primState,
                                a_dir,
                                a_box);
          }

      //apply flattening coefficient
      //would be nice to just use *= operator here
      //but they have different numbers of variables
          if (usesFlattening())
            {
              CH_TIME("apply_flatttening");
              Box entireBox = delta4W.getRegion();
              const BaseFab<Real>& regFlattening = a_flattening.getSingleValuedFAB();
              BaseFab<Real>& regDelta4W = delta4W.getSingleValuedFAB();

              FORT_APPLYFLAT(CHF_FRA(regDelta4W),
                             CHF_CONST_FRA1(regFlattening, 0),
                             CHF_CONST_INT(numSlope),
                             CHF_BOX(entireBox));

              // all single valued cells multiplied in fortran
              IntVectSet multiIVS  = m_ebisBox.getMultiCells(entireBox);
              for (VoFIterator vofit(multiIVS, m_ebisBox.getEBGraph());
                  vofit.ok(); ++vofit)
                {
                  const VolIndex& vof = vofit();
                  for (int ivar = 0; ivar < numSlope; ivar++)
                    {
                      delta4W(vof, ivar) *= a_flattening(vof, 0);
                    }
                }
            }
            {
              CH_TIME("slope_copy1");
              a_slopePrim.copy(delta4W);
            }
        }//end if (usefourthorder slopes())
       else
        {
          {
            CH_TIME("slope_copy2");
            a_slopePrim.copy(delta2W);
          }
        }
     } // end if(!useZeroSlopes())
    
    else
     {
       a_slopePrim.setVal(0.);
     }

  }

  //debug
  //a_slopePrim.setVal(0.);
  //end debug
}
/*****************************/
void
EBPatchReactive::
doSecondOrderSlopes(EBCellFAB&       a_delta2W,
                    EBCellFAB&       a_deltaWL,
                    EBCellFAB&       a_deltaWR,
                    EBCellFAB&       a_deltaWC,
                    const EBCellFAB& a_primState,
                    const int&       a_dir,
                    const Box&       a_box,
                    bool a_doAggregated)
{
  CH_TIME("EBPR::doSecondOrderSlopes");

  int numSlope;
  Box box1, box2;
  Box loBox, hiBox, centerBox, entireBox;
  int hasLo, hasHi;
  {
    //    CH_TIME("doSecondOrderSlopes_prelims");
    CH_assert(m_isSlopeSet);
    numSlope = numSlopes();
    box1 = a_box;
    box1.grow(a_dir, 1);
    box1 &= m_domain;

    box2 = a_box;
    box2.grow(a_dir, 2);
    box2 &= m_domain;

    if (usesFourthOrderSlopes())
      {
        eblohicenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
                     box2, m_domain, a_dir);
      }
    else
      {
        eblohicenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
                     box1, m_domain, a_dir);
      }
  }
  {
    CH_TIME("EBPR::doSecondOrderSlopes_fab_defs");
    a_delta2W.define(m_ebisBox, entireBox, numSlope);
    a_deltaWL.define(m_ebisBox, entireBox, numSlope);
    a_deltaWR.define(m_ebisBox, entireBox, numSlope);
    a_deltaWC.define(m_ebisBox, entireBox, numSlope);
  }
  //  { not necessary
  //  CH_TIME("EBPG::doSecondOrderSlopes_zeros");
  //  a_delta2W.setVal(0.);
  //  a_deltaWL.setVal(0.);
  //  a_deltaWR.setVal(0.);
  //  a_deltaWC.setVal(0.);
  //  }

    BaseFab<Real>& regDelta2W = a_delta2W.getSingleValuedFAB();
    BaseFab<Real>& regDeltaWC = a_deltaWC.getSingleValuedFAB();
    BaseFab<Real>& regDeltaWL = a_deltaWL.getSingleValuedFAB();
    BaseFab<Real>& regDeltaWR = a_deltaWR.getSingleValuedFAB();
    const BaseFab<Real>& regPrimState = a_primState.getSingleValuedFAB();

    {
      CH_TIME("EBPR::doSecondOrderSlopes_regular_calc");
      FORT_SECONDSLOPEDIFFS(CHF_FRA(regDeltaWC),
                            CHF_FRA(regDeltaWL),
                            CHF_FRA(regDeltaWR),
                            CHF_CONST_FRA(regPrimState),
                            CHF_CONST_INT(numSlope),
                            CHF_CONST_INT(a_dir),
                            CHF_BOX(loBox),
                            CHF_CONST_INT(hasLo),
                            CHF_BOX(hiBox),
                            CHF_CONST_INT(hasHi),
                            CHF_BOX(centerBox));
    }
  //apply van leer limiter to regular cells.
  //limiting for irregular cells is more complicated
  //now that we are doing higher order limited
  //one-sided diffs
    {
      CH_TIME("copy_wc_to_2w");
      regDelta2W.copy(regDeltaWC);
    }
  if (m_useLimiting)
    {
      CH_TIME("limiting");
      FORT_VLLIMITER(CHF_FRA(regDelta2W),
                     CHF_CONST_FRA(regDeltaWL),
                     CHF_CONST_FRA(regDeltaWR),
                     CHF_BOX(centerBox));
    }

  if (a_doAggregated && m_useAgg)
    {
      if (usesFlattening())
        {
          MayDay::Warning("agg irreg slopes does not contain flattening limiter");
        }
      aggIrregSecondOrderSlopes(a_delta2W,
                                a_deltaWL,
                                a_deltaWR,
                                a_deltaWC,
                                a_primState,
                                a_dir,
                                entireBox);
    }
  else
    {
      irregSecondOrderSlopes(a_delta2W,
                             a_deltaWL,
                             a_deltaWR,
                             a_deltaWC,
                             a_primState,
                             a_dir,
                             a_box);
    }

  if (m_isBCSet)
    {
      //      CH_TIME("boundary_slopes");
      m_bc->setBndrySlopes(a_delta2W, a_primState, m_ebisBox, entireBox, a_dir);
    }
}

/******************/
void
EBPatchReactive::
aggIrregSecondOrderSlopes(EBCellFAB&       a_delta2W,
                          EBCellFAB&       a_deltaWL,
                          EBCellFAB&       a_deltaWH,
                          EBCellFAB&       a_deltaWC,
                          const EBCellFAB& a_primState,
                          const int&       a_dir,
                          const Box&       a_entireBox)
{
  CH_TIME("agg_irreg_second_order_slopes");
  CH_assert(m_useAgg);
  CH_assert(a_entireBox == m_entireBox[a_dir]);
  CH_assert(a_delta2W.getRegion() == m_entireBox[a_dir]);
  CH_assert(a_deltaWL.getRegion() == m_entireBox[a_dir]);
  CH_assert(a_deltaWH.getRegion() == m_entireBox[a_dir]);
  CH_assert(a_deltaWC.getRegion() == m_entireBox[a_dir]);
  CH_assert(a_primState.getRegion() == m_validBoxG4);
  //the slow version has some funky flattening logic that I do
  //not have the time to program into the fast version right now as I
  //am trying to make a deadline for a code that does not use flattening.

  //get low and high slopes (or left and right if you prefer)
  m_slopStenLo[a_dir]->apply(a_deltaWL, a_primState, 0, 0,  numSlopes(), false);
  m_slopStenHi[a_dir]->apply(a_deltaWH, a_primState, 0, 0,  numSlopes(), false);

  //at irregular cells, if one side does not exist,
  //set all to one-sided diffs.  try to do higher order
  //limited diffs if possible.

    for (int ivof = 0; ivof< m_slopVec[a_dir].size(); ivof++)
      {
        //slopes have all the same box so they will have the same offset and type
        const size_t& offsetS = m_slopVec[a_dir][ivof].slop_access.offset;
        const int   & dataIDS = m_slopVec[a_dir][ivof].slop_access.dataID;

        const bool& hasFacesLo = m_slopVec[a_dir][ivof].hasLo;
        const bool& hasFacesHi = m_slopVec[a_dir][ivof].hasHi;

        for (int ivar = 0; ivar < numSlopes(); ivar++)
          {
            Real*       ptrLo =   a_deltaWL.dataPtr(dataIDS,ivar);
            Real*       ptrHi =   a_deltaWH.dataPtr(dataIDS,ivar);
            Real*       ptrCe =   a_deltaWC.dataPtr(dataIDS,ivar);
            Real*       ptrSe =   a_delta2W.dataPtr(dataIDS,ivar);

            Real& dql = *(ptrLo + offsetS);
            Real& dqr = *(ptrHi + offsetS);
            Real& dqc = *(ptrCe + offsetS);
            Real& dqS = *(ptrSe + offsetS);
            //this logic from pointgetslopes to get dqc
            if (hasFacesLo && hasFacesHi)
              {
                dqc = 0.5*(dql+dqr);
                if (m_useLimiting)
                  {
                    Real dqmin = Min(Abs(dql)*2.0,Abs(dqr)*2.0);
                    if (dql*dqr < 0.)
                      {
                        dqmin = 0.;
                      }
                    if (dqc < 0.)
                      {
                        dqc = -Min(dqmin,Abs(dqc));
                      }
                    else
                      {
                        dqc = Min(dqmin,Abs(dqc));
                      }
                  }
              }
            else if (!hasFacesLo && !hasFacesHi)
              {
                dql = 0.0;
                dqr = 0.0;
                dqc = 0.0;
              }
            else if (hasFacesLo && !hasFacesHi)
              {
                dqr = dql;
                //no higher-order one sided diffs
                dqc = dql;
              }
            else if (hasFacesHi && !hasFacesLo)
              {
                dql = dqr;
                //no higher-order one sided diffs
                dqc = dqr;
              }
            else
              {
                MayDay::Error("EBPatchReactive::aggirregpointGetSlopes -- missed a case");
              }

            //now for the limiting bits from irregSlopes
            //looks like a lot of this is redundant
            if (!m_useLimiting)
              {
                dqS = dqc;
              }
            else
              {
                if (hasFacesLo && hasFacesHi)
                  {
                    Real dqlim = dqc;

                    FORTNT_POINTVLLIMITER(CHF_REAL(dqlim),
                                          CHF_CONST_REAL(dql),
                                          CHF_CONST_REAL(dqr));
                    dqS = dqlim;
                  }
                else if (!hasFacesLo && !hasFacesHi)
                  {
                    dqS = 0.0;
                  }
                else if (hasFacesLo && !hasFacesHi)
                  {
                    if (dqc*dql > 0.0)
                      {
                        Real rsign = 1.0;
                        if (dqc < 0.0)
                          {
                            rsign = -1.0;
                          }
                        dqS = rsign*Min(Abs(dqc), Abs(dql));
                      }
                    else
                      {
                        dqS = 0.0;
                      }

                  }
                else if (hasFacesHi && !hasFacesLo)
                  {
                    if (dqc*dqr > 0.0)
                      {
                        Real rsign = 1.0;
                        if (dqc < 0.0)
                          {
                            rsign = -1.0;
                          }
                        dqS = rsign*Min(Abs(dqc), Abs(dqr));
                      }
                    else
                      {
                        dqS = 0.0;
                      }
                  }
                else
                  {
                    MayDay::Error("EBPatchReactive::aggirreg2-- missed a case");
                  }
              }//end if using limiting

          } //end loop over variables
      }//end loop over vofs
}
       
/******************/
void
EBPatchReactive::
irregSecondOrderSlopes(EBCellFAB&       a_delta2W,
                       EBCellFAB&       a_deltaWL,
                       EBCellFAB&       a_deltaWR,
                       EBCellFAB&       a_deltaWC,
                       const EBCellFAB& a_primState,
                       const int&       a_dir,
                       const Box&       a_entireBox)
{
  //at irregular cells, if one side does not exist,
  //set all to one-sided diffs.  try to do higher order
  //limited diffs if possible.
  {
    CH_TIME("irreg_slopes");
    for (int ivof = 0; ivof<m_irregVoFs.size(); ivof++)
      {
        const VolIndex& vof = m_irregVoFs[ivof];
        if (a_entireBox.contains(vof.gridIndex()))
          {
            for (int ivar = 0; ivar < numSlopes(); ivar++)
              {
                bool hasFacesLeft, hasFacesRigh;
                Real dql, dqr, dqc;
                bool verbose =false;
                pointGetSlopes(dql, dqr,dqc,
                               hasFacesLeft,
                               hasFacesRigh,
                               vof, a_primState, a_dir, ivar, verbose);
                Real dqSec=0;
                if (!m_useLimiting)
                  {
                    dqSec = dqc;
                  }
                else
                  {
                    if (hasFacesLeft && hasFacesRigh)
                      {
                        Real dqlim = dqc;

                        FORT_POINTVLLIMITER(CHF_REAL(dqlim),
                                            CHF_CONST_REAL(dql),
                                            CHF_CONST_REAL(dqr));
                        dqSec = dqlim;
                      }
                    else if (!hasFacesLeft && !hasFacesRigh)
                      {
                        dqSec = 0.0;
                      }
                    else if (hasFacesLeft && !hasFacesRigh)
                      {
                        if (dqc*dql > 0.0)
                          {
                            Real rsign = 1.0;
                            if (dqc < 0.0)
                              {
                                rsign = -1.0;
                              }
                            dqSec = rsign*Min(Abs(dqc), Abs(dql));
                          }
                        else
                          {
                            dqSec = 0.0;
                          }

                      }
                    else if (hasFacesRigh && !hasFacesLeft)
                      {
                        if (dqc*dqr > 0.0)
                          {
                            Real rsign = 1.0;
                            if (dqc < 0.0)
                              {
                                rsign = -1.0;
                              }
                            dqSec = rsign*Min(Abs(dqc), Abs(dqr));
                          }
                        else
                          {
                            dqSec = 0.0;
                          }
                      }
                    else
                      {
                        MayDay::Error("EBPatchReactive::doSecondOrderSlopes -- missed a case");
                      }
                  }//end if using limiting

                if (usesFlattening())
                  {

                    int pindex = pressureIndex();
                    Real press = Max(a_primState(vof, pindex), Real(1.0e-10));
                    Real dqp, dqlp, dqrp;
                    verbose = false;
                    pointGetSlopes(dqlp, dqrp, dqp,
                                   hasFacesLeft,
                                   hasFacesRigh,
                                   vof, a_primState, a_dir, pindex,verbose);
                    //try to detect pressure discontinutity
                    dqp = Max(Abs(dqp), Abs(dqlp));
                    dqp = Max(Abs(dqp), Abs(dqrp));
                    if (Abs(dqp/press) > 0.1)
                      {
                        dqSec = 0.0;
                      }
                  }
                a_delta2W(vof, ivar) = dqSec;
                a_deltaWL(vof, ivar) = dql;
                a_deltaWR(vof, ivar) = dqr;
                a_deltaWC(vof, ivar) = dqc;
              } //end loop over variables
          }
      }//end loop over vofs
  }
}
/*****************************/
void EBPatchReactive::
pointGetSlopes(Real&            a_dql,
               Real&            a_dqr,
               Real&            a_dqc,
               bool&            a_hasFacesLeft,
               bool&            a_hasFacesRigh,
               const VolIndex&  a_vof,
               const EBCellFAB& a_primState,
               const int&       a_dir,
               const int&       a_ivar,
               const bool&      a_verbose)
{
 
   const IntVect&  iv = a_vof.gridIndex();
  //one-sided diffs on domain bndry
  const Box& domainBox = m_domain.domainBox();
  bool onLeftDomain = (iv[a_dir] == domainBox.smallEnd(a_dir));
  bool onRighDomain = (iv[a_dir] == domainBox.bigEnd(a_dir)  );
  VolIndex vofLeft;
  VolIndex vofRigh;
  a_hasFacesLeft = (m_ebisBox.numFaces(a_vof, a_dir, Side::Lo)==1) && !onLeftDomain;
  a_hasFacesRigh = (m_ebisBox.numFaces(a_vof, a_dir, Side::Hi)==1) && !onRighDomain;

  if (a_hasFacesLeft)
    {
      Vector<FaceIndex> facesLeft =
        m_ebisBox.getFaces(a_vof, a_dir, Side::Lo);
      vofLeft = facesLeft[0].getVoF(Side::Lo);
    }

  if (a_hasFacesRigh)
    {
      Vector<FaceIndex> facesRigh =
        m_ebisBox.getFaces(a_vof, a_dir, Side::Hi);
      vofRigh = facesRigh[0].getVoF(Side::Hi);
    }

  Real valCent = a_primState(a_vof, a_ivar);
  Real valLeft = 0.0;
  Real valRigh = 0.0;
  if (a_hasFacesLeft)
    {
      valLeft = a_primState(vofLeft, a_ivar);
    }
  if (a_hasFacesRigh)
    {
      valRigh = a_primState(vofRigh, a_ivar);
    }

  //at irregular cells, if one side does not exist,
  //set all to one-sided diffs.
  //if neither exists, set all slopes to zero
  //limiting sprinked in here to make higher order
  //one-sided diffs possible
  if (a_hasFacesLeft)
    {
      a_dql = valCent - valLeft;
    }
  if (a_hasFacesRigh)
    {
      a_dqr = valRigh - valCent;
    }
  if (a_hasFacesLeft && a_hasFacesRigh)
    {
      a_dqc = 0.5*(a_dql+a_dqr);
      if (m_useLimiting)
        {
          Real dqmin = Min(Abs(a_dql)*2.0,Abs(a_dqr)*2.0);
          if (a_dql*a_dqr < 0.)
            {
              dqmin = 0.;
            }
          if (a_dqc < 0.)
            {
              a_dqc = -Min(dqmin,Abs(a_dqc));
            }
          else
            {
              a_dqc = Min(dqmin,Abs(a_dqc));
            }
        }
    }
  else if (!a_hasFacesLeft && !a_hasFacesRigh)
    {
      a_dql = 0.0;
      a_dqr = 0.0;
      a_dqc = 0.0;
    }
  else if (a_hasFacesLeft && !a_hasFacesRigh)
    {
      a_dqr = a_dql;
      //no higher-order one sided diffs
      a_dqc = a_dql;
    }
  else if (a_hasFacesRigh && !a_hasFacesLeft)
    {
      a_dql = a_dqr;
      //no higher-order one sided diffs
      a_dqc = a_dqr;
    }
  else
    {
      MayDay::Error("EBPatchReactive::pointGetSlopes -- missed a case");
    }
  if (a_verbose)
    {
      pout() << "  a_vof="    << a_vof
             << ", vofLeft= " << vofLeft
             << ", vofRigh= " << vofRigh  << endl;
      pout() << "  hasFacesLeft=" << a_hasFacesLeft
             << ", hasFacesRigh=" << a_hasFacesRigh;
      pout() << ", valLeft=" << valLeft
             << ", valRigh=" << valRigh
             << ", valCent=" << valCent << endl;
      pout() << "  a_dql=" << a_dql
             << ", a_dqr=" << a_dqr
             << ", a_dqc=" << a_dqc << endl;
    }
  if (a_dqc != a_dqc)
    {
      // debug
      FabDataOps::getFabData(a_primState,1);
      // end debug
      MayDay::Error("EBPatchReactive::pointGetSlopes -- a_dqc != a_dqc");
    }
  //atttempts to detect discontinuties in constrained vars
  if (usesFlattening())
    {
      if ((Abs(a_dqc) > Abs(0.1*valCent)) && ((a_ivar == densityIndex()) || (a_ivar == pressureIndex())))
        {
          a_dql = 0.0;
          a_dqr = 0.0;
          a_dqc = 0.0;
        }
    }
  //
}

/*****************************/
void EBPatchReactive::
doFourthOrderSlopes(EBCellFAB&       a_delta4W,
                    EBCellFAB&       a_deltaWC,
                    const EBCellFAB& a_delta2W,
                    const EBCellFAB& a_deltaWL,
                    const EBCellFAB& a_deltaWR,
                    const EBCellFAB& a_primState,
                    const int&       a_dir,
                    const Box&       a_box)
{
  CH_TIME("dofourthorderslopes");
  Box loBox, hiBox, centerBox, entireBox;
  int hasLo, hasHi;
  int numSlope = numSlopes();
  const Box& domainBox = m_domain.domainBox();
  Box box1 = a_box;
  box1.grow(a_dir, 1);
  box1 &= m_domain;
  eblohicenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
               box1, m_domain, a_dir);

  a_delta4W.define(m_ebisBox, entireBox, numSlope);

  BaseFab<Real>& regDelta4W = a_delta4W.getSingleValuedFAB();

  const BaseFab<Real>& regDelta2W = a_delta2W.getSingleValuedFAB();
  //const BaseFab<Real>& regDeltaWC = a_deltaWC.getSingleValuedFAB();
  //const BaseFab<Real>& regDeltaWL = a_deltaWL.getSingleValuedFAB();
  //const BaseFab<Real>& regDeltaWR = a_deltaWR.getSingleValuedFAB();
  const BaseFab<Real>& regPrimState = a_primState.getSingleValuedFAB();
  {
    CH_TIME("reg4thOrderSlopes");
    FORT_FORTHSLOPEDIFFS(CHF_FRA(regDelta4W),
                         CHF_CONST_FRA(regPrimState),
                         CHF_CONST_FRA(regDelta2W),
                         CHF_CONST_INT(numSlope),
                         CHF_CONST_INT(a_dir),
                         CHF_BOX(loBox),
                         CHF_CONST_INT(hasLo),
                         CHF_BOX(hiBox),
                         CHF_CONST_INT(hasHi),
                         CHF_BOX(centerBox));
  }

  {
    CH_TIME("irreg4thOrderSlopes");
    for (int ivof = 0; ivof<m_irregVoFs.size(); ivof++)
      {
        const VolIndex& vof = m_irregVoFs[ivof];
        if (a_box.contains(vof.gridIndex()))
          {
            const IntVect&  iv = vof.gridIndex();
            //one-sided diffs on domain bndry
            bool onLeftDomain = iv[a_dir] == domainBox.smallEnd(a_dir);
            bool onRighDomain = iv[a_dir] == domainBox.bigEnd(a_dir);
            bool hasFacesLeft = (m_ebisBox.numFaces(vof, a_dir, Side::Lo) > 0) && !onLeftDomain;
            bool hasFacesRigh = (m_ebisBox.numFaces(vof, a_dir, Side::Hi) > 0) && !onRighDomain;

            for (int ivar = 0; ivar < numSlopes(); ivar++)
              {
                Real dq4;
                if (hasFacesRigh && hasFacesLeft)
                  {
                    Vector<FaceIndex> facesLeft =
                      m_ebisBox.getFaces(vof, a_dir, Side::Lo);
                    Real valLeft = 0.0;
                    Real sloLeft = 0.0;
                    for (int iface = 0; iface < facesLeft.size(); iface++)
                      {
                        VolIndex vofLeft = facesLeft[iface].getVoF(Side::Lo);
                        valLeft += a_primState(vofLeft, ivar);
                        sloLeft +=      a_delta2W(vofLeft, ivar);
                      }
                    valLeft /= Real(facesLeft.size());
                    sloLeft /= Real(facesLeft.size());

                    Vector<FaceIndex> facesRigh =
                      m_ebisBox.getFaces(vof, a_dir, Side::Hi);
                    Real valRigh = 0.0;
                    Real sloRigh = 0.0;
                    for (int iface = 0; iface <facesRigh.size(); iface++)
                      {
                        VolIndex vofRigh = facesRigh[iface].getVoF(Side::Hi);
                        valRigh += a_primState(vofRigh, ivar);
                        sloRigh +=     a_delta2W(vofRigh, ivar);
                      }
                    valRigh /= Real(facesRigh.size());
                    sloRigh /= Real(facesRigh.size());

                    Real dwr = valRigh - 0.25*sloRigh;
                    Real dwl = valLeft + 0.25*sloLeft;
                    dq4 = (2.0/3.0)*(dwr - dwl);

                  }
                else if (hasFacesRigh || hasFacesLeft)
                  {
                    dq4 = a_delta2W(vof, ivar);
                  }
                else
                  {
                    //no faces on either side.  set the slope to zero
                    dq4 = 0.0;
                  }
                a_delta4W(vof, ivar) = dq4;
              } //end loop over variables
          }
      }//end loop over vofs
    }

  if (m_useLimiting)
    {
      CH_TIME("4th order limiting");
      applyLimiter(a_delta4W, a_deltaWL, a_deltaWR,  a_dir, centerBox);
    }

}
/******/
void EBPatchReactive::
applyLimiter(EBCellFAB&       a_slopePrim,
             const EBCellFAB& a_slopeLeft,
             const EBCellFAB& a_slopeRigh,
             const int&       a_dir,
             const Box&       a_box)
{
  CH_TIME("EBPatchReactive::applyLimiter");
  CH_assert(isDefined());
  CH_assert(a_slopePrim.getRegion().contains(a_box));
  CH_assert(a_slopeLeft.getRegion().contains(a_box));
  CH_assert(a_slopeRigh.getRegion().contains(a_box));
  CH_assert((a_dir >= 0) && (a_dir < SpaceDim));
  CH_assert(a_slopePrim.nComp() == numSlopes());
  CH_assert(a_slopeLeft.nComp() == numSlopes());
  CH_assert(a_slopeRigh.nComp() == numSlopes());

  const BaseFab<Real>& regSlopeLeft = a_slopeLeft.getSingleValuedFAB();
  const BaseFab<Real>& regSlopeRigh = a_slopeRigh.getSingleValuedFAB();
  BaseFab<Real>& regSlopePrim       = a_slopePrim.getSingleValuedFAB();

  {
    CH_TIME("regular limiting");
    FORT_VLLIMITER(CHF_FRA(regSlopePrim),
                   CHF_CONST_FRA(regSlopeLeft),
                   CHF_CONST_FRA(regSlopeRigh),
                   CHF_BOX(a_box));
  }

  {
    CH_TIME("irregular limiting");
    IntVectSet ivsMulti = m_ebisBox.getMultiCells(a_box);
    for (VoFIterator vofit(ivsMulti, m_ebisBox.getEBGraph());
        vofit.ok(); ++vofit)
      {
        const VolIndex& vof = vofit();
        for (int ivar = 0; ivar < numSlopes(); ivar++)
          {
            //the output array arrives with the centered diff value
            const Real& dql = a_slopeLeft(vof, ivar);
            const Real& dqr = a_slopeRigh(vof, ivar);
            Real& dqlim = a_slopePrim(vof, ivar);

            FORT_POINTVLLIMITER(CHF_REAL(dqlim),
                                CHF_CONST_REAL(dql),
                                CHF_CONST_REAL(dqr));
          }
      }
  }
}
/*****************************/
void EBPatchReactive::
normalPred(EBCellFAB&       a_primLo,
           EBCellFAB&       a_primHi,
           const EBCellFAB& a_primState,
           const EBCellFAB& a_slopePrim,
           const Real&      a_dtbydx,
           const int&       a_dir,
           const Box&       a_box)
{
  CH_TIME("EBPatchReactive::normalPred");
  int nslope =  numSlopes();
  CH_assert(isDefined());
  CH_assert(a_primLo.nComp()    >= nslope);
  CH_assert(a_primHi.nComp()    >= nslope);
  CH_assert(a_primState.nComp() >= nslope);
  CH_assert(a_slopePrim.nComp() == nslope);

  Real dtbydx = a_dtbydx;
  const BaseFab<Real>& regState = a_primState.getSingleValuedFAB();
  const BaseFab<Real>& regSlope = a_slopePrim.getSingleValuedFAB();
  BaseFab<Real>& regPrimLo = a_primLo.getSingleValuedFAB();
  BaseFab<Real>& regPrimHi = a_primHi.getSingleValuedFAB();
  int useflat = 0;
  if (usesFlattening())
    useflat = 1;

  FORT_PRED(CHF_BOX(a_box),
            CHF_CONST_FRA(regState),
            CHF_CONST_FRA(regSlope),
            CHF_FRA(regPrimLo),
            CHF_FRA(regPrimHi),
            CHF_CONST_INT(a_dir),
            CHF_CONST_REAL(dtbydx),
            CHF_CONST_INT(useflat));

  IntVectSet ivsMulti = m_ebisBox.getMultiCells(a_box);
  for (VoFIterator vofit(ivsMulti, m_ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      Vector<Real> primlo(numPrimitives()), primhi(numPrimitives());
      Vector<Real> primit(numPrimitives()), pslope(numPrimitives());

            for (int ivar = 0; ivar < numPrimitives(); ivar++)
        {
          primit[ivar] = a_primState(vof, ivar);
          pslope[ivar] = a_slopePrim(vof, ivar);
        }

      FORT_POINTPRED(CHF_VR(primit),
                     CHF_VR(pslope),
                     CHF_VR(primlo),
                     CHF_VR(primhi),
                     CHF_CONST_INT(a_dir),
                     CHF_CONST_REAL(dtbydx),
                     CHF_CONST_INT(useflat));


      for (int ivar = 0; ivar < numPrimitives(); ivar++)
        {
          a_primLo(vof, ivar) = primlo[ivar];
          a_primHi(vof, ivar) = primhi[ivar];
        }
    }
}
/******/
void
EBPatchReactive::
incrementWithSource(EBCellFAB&       a_state,
                    const EBCellFAB& a_source,
                    const Real&      a_scale,
                    const Box&       a_box)
{
  //if conservative source is set to true,  state is conservative state
  //if conservative source is set to false, state is primitive    state
  CH_TIME("EBPatchReactive::incrmentWithSource");
  CH_assert(m_isBoxSet);
  CH_assert(a_source.nComp() == a_state.nComp());
  CH_assert(a_source.getRegion().contains(a_box));
  CH_assert(a_state.getRegion().contains(a_box));

  BaseFab<Real>&       primReg = a_state.getSingleValuedFAB();
  const BaseFab<Real>& sourReg =    a_source.getSingleValuedFAB();

  FORT_INCSOURCE(CHF_FRA(primReg),
                 CHF_CONST_FRA(sourReg),
                 CHF_CONST_REAL(a_scale),
                 CHF_BOX(a_box));

  IntVectSet ivsMulti = m_ebisBox.getMultiCells(a_box);
  for (VoFIterator vofit(ivsMulti, m_ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      for (int ivar = 0; ivar < a_source.nComp(); ivar++)
        {
          Real increm = a_scale*a_source(vof, ivar);
          a_state(vof, ivar) +=  increm;
        }
    }
}
/******/
void
EBPatchReactive::
primToCons(EBCellFAB&       a_consState,
           const EBCellFAB& a_primState,
           const Box&       a_box)
{
  CH_TIME("EBPatchReactive::primToCons");
  CH_assert(isDefined());
  const Box& regBox = a_box;
  CH_assert(a_primState.getRegion().contains(regBox));
  //set covered vals. does not change real data.
  //  EBCellFAB& primData = (EBCellFAB&) a_primState;
  //  setCoveredPrimVals(primData);

  BaseFab<Real>&   regCons = a_consState.getSingleValuedFAB();
  const BaseFab<Real>&  regPrim = a_primState.getSingleValuedFAB();


  FORT_PRM2CONS(CHF_BOX(regBox),
                CHF_FRA(regCons),
                CHF_CONST_FRA(regPrim));


  IntVectSet ivsMulti = m_ebisBox.getMultiCells(regBox);
  for (VoFIterator vofit(ivsMulti, m_ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();

      Vector<Real> conserved(numConserved());
      Vector<Real> primitive(numPrimitives());
      for (int ivar = 0; ivar < numPrimitives(); ivar++)
        {
          primitive[ivar] = a_primState(vof, ivar);
        }

      FORT_POINTPRM2CONS(CHF_VR(conserved),
                         CHF_VR(primitive));

      for (int ivar = 0; ivar < numConserved(); ivar++)
        {
          a_consState(vof, ivar)  = conserved[ivar];
        }
    }
  //  setCoveredConsVals(a_consState);
}
/******/
void
EBPatchReactive::
primToCons(BaseIVFAB<Real>&       a_consState,
           const BaseIVFAB<Real>& a_primState,
           const IntVectSet&      a_ivs)
{
  CH_TIME("EBPatchReactive::primToConsIrr");
  CH_assert(isDefined());
  for (VoFIterator vofit(a_ivs, m_ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();

      Vector<Real> conserved(numConserved());
      Vector<Real> primitive(numPrimitives());
      for (int ivar = 0; ivar < numPrimitives(); ivar++)
        {
          primitive[ivar] = a_primState(vof, ivar);
        }

      FORT_POINTPRM2CONS(CHF_VR(conserved),
                         CHF_VR(primitive));

      for (int ivar = 0; ivar < numConserved(); ivar++)
        {
          a_consState(vof, ivar)  = conserved[ivar];
        }
    }
}
/*****************************/
void
EBPatchReactive::
riemann(EBFaceFAB&       a_flux,
        const EBCellFAB& a_primLeft,
        const EBCellFAB& a_primRigh,
        const int&       a_dir,
        const Box&       a_box)
{
  CH_TIME("EBPatchReactive::riemann");
  //int nPrim = numPrimitives();
  Box cellBox = enclosedCells(a_box);

  //explicitly cast the left and right states to modifiable references
  //because we need to shift them to faces and then shift them back.
  BaseFab<Real>& regPrimRigh = (BaseFab<Real>&)a_primRigh.getSingleValuedFAB();
  BaseFab<Real>& regPrimLeft = (BaseFab<Real>&)a_primLeft.getSingleValuedFAB();
  //shift the cell centered stuff to edges
  regPrimRigh.shiftHalf(a_dir, -1);
  regPrimLeft.shiftHalf(a_dir,  1);

  CH_assert(regPrimRigh.box().contains(a_box));
  CH_assert(regPrimLeft.box().contains(a_box));
  CH_assert(a_flux.direction() == a_dir);
  CH_assert(a_flux.getRegion().contains(a_box));
  //find the fluxes at faces of singlevalued enclosed cells

  BaseFab<Real>& regFlux = a_flux.getSingleValuedFAB();

  CH_assert(regPrimLeft.nComp() == numPrimitives());
  CH_assert(regPrimRigh.nComp() == numPrimitives());

// regular fluxes
  FORT_RIEMANN(CHF_BOX(a_box),
                  CHF_CONST_FRA(regPrimLeft),
                  CHF_CONST_FRA(regPrimRigh),
                  CHF_FRA(regFlux),
                  CHF_CONST_INT(a_dir),
                  CHF_INT(m_nSpec));

  //shift the cell centered stuff back to cells
  regPrimRigh.shiftHalf(a_dir,  1);
  regPrimLeft.shiftHalf(a_dir, -1);

  //the box sent into this is face-centered.
  //we need to use the cell-centered one it surrounds.
  //this can be more than multivalued cells if there
  //are neighbors that are multivalued

  Box grownBox = cellBox;
  grownBox.grow(a_dir, 1);
  grownBox &= m_domain;
  IntVectSet ivsMulti = m_ebisBox.getIrregIVS(grownBox);

  FaceStop::WhichFaces stopCrit = FaceStop::SurroundingNoBoundary;
  for (FaceIterator faceit(ivsMulti, m_ebisBox.getEBGraph(), a_dir, stopCrit);
      faceit.ok(); ++faceit)
    {
      const FaceIndex& face = faceit();
      if (a_box.contains(face.gridIndex(Side::Hi)))
        {
          VolIndex vofl = face.getVoF(Side::Lo);
          VolIndex vofr = face.getVoF(Side::Hi);
          Vector<Real> ql(numPrimitives()), qr(numPrimitives()), fluxvec(numFluxes());
          for (int ivar = 0; ivar < numPrimitives(); ivar++)
            {
              ql[ivar] = a_primLeft(vofl, ivar);
              qr[ivar] = a_primRigh(vofr, ivar);
            }


          FORT_POINTRIEMANN(CHF_VR(ql), CHF_VR(qr), CHF_VR(fluxvec),
                            CHF_CONST_INT(a_dir), CHF_INT(m_nSpec));


          for (int ivar = 0; ivar < numFluxes(); ivar++)
            {
              a_flux(face, ivar) = fluxvec[ivar];
            }
        }
    }

  // fix boundary face
  for (SideIterator sit; sit.ok(); ++sit)
    {
      Box boundBox;
      if (sit() == Side::Lo)
        {
          boundBox = adjCellLo(m_domain, a_dir, 1);
        }
      else
        {
          boundBox = adjCellHi(m_domain, a_dir, 1);
        }
      boundBox.shift(a_dir, -sign(sit()));
      IntVectSet ivsBound(boundBox);
      ivsBound &= m_validBox;
      if (!ivsBound.isEmpty())
        {
          FaceStop::WhichFaces stopCrit = FaceStop::AllBoundaryOnly;
          int nPrim = numPrimitives();
          int nFlux = numFluxes();
          for (FaceIterator faceit(ivsBound, m_ebisBox.getEBGraph(), a_dir, stopCrit);
              faceit.ok(); ++faceit)
            {
              Vector<Real> prim(nPrim), flux(nFlux);
              const FaceIndex& face = faceit();

              for (int ivar = 0; ivar < nPrim; ivar++)
                {
                  Real primExtrap;
                  if (sit() == Side::Lo)
                    {
                      primExtrap = a_primRigh(faceit().getVoF(Side::Hi), ivar);
                    }
                  else
                    {
                      primExtrap = a_primLeft(faceit().getVoF(Side::Lo), ivar);
                    }
                  prim[ivar] = primExtrap;
                }
              FORT_POINTGETFLUX(CHF_VR(flux),
                                CHF_VR(prim),
                                CHF_CONST_INT(a_dir));

             for (int ivar = 0; ivar < numFluxes(); ivar++)
              {
              a_flux(face, ivar) = flux[ivar];
              }

            }
        }
    }
}
/*****************************/
void
EBPatchReactive::
riemann(BaseIVFAB<Real>&        a_coveredFlux,
        const BaseIVFAB<Real>&  a_exteState,
        const EBCellFAB&        a_primState,
        const Vector<VolIndex>& a_vofset,
        const int&              a_dir,
        const Side::LoHiSide&   a_sd,
        const Box&       a_box)
{
  CH_TIME("EBPatchReactive::riemannIrr");
  int nPrim = numPrimitives();
  int nFlux = numFluxes();
  for (int ivof = 0; ivof < a_vofset.size(); ivof++)
    {
      const VolIndex& vof = a_vofset[ivof];
      if (a_box.contains(vof.gridIndex()))
        {
          Vector<Real> ql(nPrim), qr(nPrim), flux(nFlux);

          if (a_sd == Side::Hi)
            {
              for (int ivar = 0; ivar < nPrim; ivar++)
                {
                  ql[ivar] = a_primState(vof, ivar);
                  qr[ivar] = a_exteState(vof, ivar);
                }
            }
          else
            {
              for (int ivar = 0; ivar < nPrim; ivar++)
                {
                  ql[ivar] = a_exteState(vof, ivar);
                  qr[ivar] = a_primState(vof, ivar);
                }
            }


          FORT_POINTRIEMANN(CHF_VR(ql), CHF_VR(qr), CHF_VR(flux),
                            CHF_CONST_INT(a_dir), CHF_INT(m_nSpec));

          for (int ivar = 0; ivar < nFlux; ivar++)
            {
              a_coveredFlux(vof, ivar) = flux[ivar];
            }
        }
    }
}
/******/
void
EBPatchReactive::
extrapToCoveredFaces(BaseIVFAB<Real>&        a_extendedPrim,
                     const EBCellFAB&        a_primMinu,
                     const EBCellFAB&        a_primPlus,
                     const EBCellFAB&        a_primState,
                     const Vector<VolIndex>& a_coveredFaces,
                     const int&              a_faceDir,
                     const Side::LoHiSide&   a_sd,
                     const Box&              a_box)
{
  CH_TIME("EBPatchReactive::extrapToCoveredFaces");
/*  
  //debug
  pout() << "primState in extrapToCoveredFaces" << endl;
  FabDataOps::getFabData(a_primState,1);

  pout() << "primPlus in extrapToCoveredFaces" << endl;
  FabDataOps::getFabData(a_primPlus,0);

  pout() << "primMinu in extrapToCoveredFaces" << endl;
  FabDataOps::getFabData(a_primMinu,0);
  // end debug 
*/
  for (int ivof = 0; ivof < a_coveredFaces.size(); ivof++)
    {
      const VolIndex& vof = a_coveredFaces[ivof];

      RealVect normal = m_ebisBox.normal(vof);
      if (a_box.contains(vof.gridIndex()))
        {
          const int numPrim = numPrimitives();
          Vector<Real> extPrim(numPrim, 0.0);

          if (SpaceDim== 2)
            {
              pointExtrapToCovered2D(extPrim,
                                     a_primMinu,
                                     a_primPlus,
                                     a_primState,
                                     a_faceDir,
                                     vof,
                                     normal,
                                     a_sd,
                                     numPrim);
            }
          else if (SpaceDim==3)
            {
              pointExtrapToCovered3D(extPrim,
                                     a_primMinu,
                                     a_primPlus,
                                     a_primState,
                                     a_faceDir,
                                     vof,
                                     normal,
                                     a_sd,
                                     numPrim);
            }
          else
            {
              MayDay::Error("Bogus SpaceDim");
            }

          for (int ivar = 0; ivar < numPrim; ivar++)
            {
              a_extendedPrim(vof, ivar) = extPrim[ivar];
            }
        }
    }

  const IntVectSet& ivs  = a_extendedPrim.getIVS();
  floorPrimitives(a_extendedPrim, ivs);
}
/********/
void
EBPatchReactive::
pointExtrapToCovered2D(Vector<Real>&           a_extrapVal,
                       const EBCellFAB&        a_primMinu,
                       const EBCellFAB&        a_primPlus,
                       const EBCellFAB&        a_primState,
                       const int&              a_faceDir,
                       const VolIndex&         a_vof,
                       const RealVect&         a_normal,
                       const Side::LoHiSide&   a_sd,
                       const int&              a_numPrim)
{
  CH_assert(SpaceDim==2);
  int tangenDir = 1 - a_faceDir;

  int signNorm = 1;
  int signTang = 1;
  if (a_normal[a_faceDir] < 0.0) signNorm = -1;
  if (Abs(a_normal[tangenDir]) == 0.0)
    {//since signTang is arbitrary for this case,
      //check if we are next to the high domain edge
      IntVect ivHi = a_vof.gridIndex();;
      ivHi[tangenDir] += 1;
      if (!m_domain.contains(ivHi))
        {
          signTang = -1;
        }
    }
  else if (a_normal[tangenDir] < 0.0)
    {
      signTang = -1;
    }

  //iv[0][0] is the vof. iv[1][0] is vofside.  iv[0][1] is vofup iv[1][1] is vofcorner
  IntVect   ivSten[2][2];
  bool      hasVoF[2][2];
  VolIndex     vof[2][2];
  Real         val[2][2];
  ivSten[0][0] = a_vof.gridIndex();
  ivSten[1][0] = ivSten[0][0] + signNorm*BASISV(a_faceDir);
  ivSten[0][1] = ivSten[0][0] - signNorm*BASISV(a_faceDir) + signTang*BASISV(tangenDir) ;;
  ivSten[1][1] = ivSten[0][0]                              + signTang*BASISV(tangenDir) ;

  int radius = 1;
  Vector<VolIndex> vofsStencil;
  EBArith::getAllVoFsInMonotonePath(vofsStencil, a_vof, m_ebisBox, radius);

  for (int ivar = 0; ivar < a_numPrim; ivar++)
    {
      for (int ix = 0; ix < 2; ix++)
        {
          for (int iy = 0; iy < 2; iy++)
            {
              hasVoF[ix][iy] = EBArith::isVoFHere(vof[ix][iy], vofsStencil, ivSten[ix][iy]);
              if (hasVoF[ix][iy])
                {
                  if (a_sd == Side::Hi)
                    {
                      val[ix][iy] = a_primMinu(vof[ix][iy], ivar);
                    }
                  else
                    {
                      val[ix][iy] = a_primPlus(vof[ix][iy], ivar);
                    }
                }
              else
                {
                  val[ix][iy] = 0.0;
                }
            }
        }

      int dirBigNorm, dirLitNorm;
      VolIndex vofBigNorm, vofLitNorm;
      bool hasVoFBigNorm,  hasVoFLitNorm;
      int  signBigNorm, signLitNorm;
      Real wBigNorm, wLitNorm;
      Real dWBigNorm[2],dWLitNorm[2];
      Real eps = 1.0e-12;
      if ((Abs(a_normal[a_faceDir]) < eps) || (Abs(a_normal[tangenDir]/a_normal[a_faceDir]) > m_dx[tangenDir]/m_dx[a_faceDir]))
        {
          dirBigNorm = tangenDir;
          dirLitNorm = a_faceDir;

          signBigNorm   = signTang;
          signLitNorm   = signNorm;
          hasVoFBigNorm = hasVoF[0][1];
          if (hasVoF[0][1])
            {
              wBigNorm      =    val[0][1];
              vofBigNorm    =    vof[0][1];
            }
          else
            {
              wBigNorm      =    val[0][0];
              vofBigNorm    =    vof[0][0];
            }
        }
      else
        {
          dirBigNorm = a_faceDir;
          dirLitNorm = tangenDir;

          signBigNorm   = signNorm;
          signLitNorm   = signTang;
          hasVoFBigNorm = hasVoF[1][0];

          if (hasVoF[1][0])
            {
              //              Real deltaW = a_slopesSeco[a_faceDir](vof[1][0]  , ivar);
              Real deltaW;
              coveredExtrapSlopes(deltaW, vof[1][0],a_primState, a_faceDir, ivar);
              wBigNorm      =    val[1][0] - signBigNorm*deltaW;
              vofBigNorm    =    vof[1][0];
            }
          else
            {
              wBigNorm    =    val[0][0];
              vofBigNorm  =    vof[0][0];
            }

        }
      hasVoFLitNorm = hasVoF[1][1];
      if (hasVoF[1][1])
        {
          wLitNorm      =    val[1][1];
          vofLitNorm    =    vof[1][1];
        }
      else
        {
          wLitNorm      =    val[0][0];
          vofLitNorm    =    vof[0][0];
        }

      if (hasVoFLitNorm && hasVoFBigNorm)
        {
          coveredExtrapSlopes(dWBigNorm[0],vofBigNorm,a_primState,0,ivar);
          coveredExtrapSlopes(dWBigNorm[1],vofBigNorm,a_primState,1,ivar);
          coveredExtrapSlopes(dWLitNorm[0],vofLitNorm,a_primState,0,ivar);
          coveredExtrapSlopes(dWLitNorm[1],vofLitNorm,a_primState,1,ivar);

          const Real xLit = Abs( (m_dx[dirBigNorm]*a_normal[dirLitNorm])/ (m_dx[dirLitNorm]*a_normal[dirBigNorm]) );
          const Real xBig = 1.0 - xLit;

          const Real dWBigExtrap = xLit*dWLitNorm[dirBigNorm] + xBig*dWBigNorm[dirBigNorm];
          const Real dWLitExtrap = xLit*dWLitNorm[dirLitNorm] + xBig*dWBigNorm[dirLitNorm];

          const Real wc = xLit*wLitNorm + xBig*wBigNorm;
          a_extrapVal[ivar] = wc - xBig*signBigNorm*dWBigExtrap - xLit*signLitNorm*dWLitExtrap;
          //debug turn off extrapolation bit
          //a_extrapVal[ivar] = wc;
          //end debug

          if (usesFlattening())
            {
              int pindex = pressureIndex();
              int dindex = densityIndex();

              if ((ivar == pindex) || (ivar == dindex))
                {
                  //if extrapolated to a negative density or pressure then turn
                  //off the slopes of the extrapolation
                  if (a_extrapVal[ivar] < 0.0)
                    {
                      a_extrapVal[ivar] = xLit*wLitNorm + xBig*wBigNorm;
                    }
                }
            }
        }
      else if (hasVoF[1][0])
        {//change methods, still second order
          Real deltaW;
          coveredExtrapSlopes(deltaW,vof[1][0],a_primState,a_faceDir,ivar);

          a_extrapVal[ivar] = val[1][0] - 2.0*signNorm*deltaW;
        }
      else if (hasVoFLitNorm)
        {//drop order
          a_extrapVal[ivar] = wLitNorm;
        }
      else if (hasVoFBigNorm)
        {
          a_extrapVal[ivar] = wBigNorm;
        }
      else
        {
          a_extrapVal[ivar] = a_primState(a_vof,ivar);
        }
    }
}

/*****************************/
void EBPatchReactive::
coveredExtrapSlopes(Real&            a_dq,
                    const VolIndex&  a_vof,
                    const EBCellFAB& a_primState,
                    const int&       a_dir,
                    const int&       a_ivar)
{
  CH_TIME("EBPatchReactive::coveredExtrapSlopes");
  bool hasFacesLeft, hasFacesRigh;
  Real dql, dqr;
  bool verbose = false;
  pointGetSlopes(dql, dqr,a_dq,
                 hasFacesLeft,
                 hasFacesRigh,
                 a_vof, a_primState, a_dir, a_ivar, verbose);
  if (usesFlattening())
    {

      int pindex = pressureIndex();
      int rindex = densityIndex();
      Real press = Max(a_primState(a_vof, pindex), Real(1.0e-10));
      Real dense = Max(a_primState(a_vof, rindex), Real(1.0e-10));
      Real dqp, dqd;

      verbose = false;
      pointGetSlopes(dql, dqr, dqp,
                     hasFacesLeft,
                     hasFacesRigh,
                     a_vof, a_primState, a_dir, pindex,verbose);

      //try to detect pressure discontinutity
      dqp = Max(Abs(dqp), Abs(dql));
      dqp = Max(Abs(dqp), Abs(dqr));
      pointGetSlopes(dql, dqr, dqd,
                     hasFacesLeft,
                     hasFacesRigh,
                     a_vof, a_primState, a_dir, rindex,verbose);

      //try to detect denisty discontinutity
      dqd = Max(Abs(dqd), Abs(dql));
      dqd = Max(Abs(dqd), Abs(dqr));
      if ((Abs(dqp/press) > 0.1) || (Abs(dqd/dense) > 0.1))
        {
          a_dq = 0.0;
        }
    }
  //debug! set slopes to zero!
  // a_dq = 0;
  //end debug
}
/*****************************/
void
EBPatchReactive::
pointExtrapToCovered3D(Vector<Real>&           a_extrapVal,
                       const EBCellFAB&        a_primMinu,
                       const EBCellFAB&        a_primPlus,
                       const EBCellFAB&        a_primState,
                       const int&              a_faceDir,
                       const VolIndex&         a_vof,
                       const RealVect&         a_normal,
                       const Side::LoHiSide&   a_sd,
                       const int&              a_numPrim)
{
  RealVect normal = a_normal;
  CH_assert(SpaceDim==3);
  bool dropExtrap = false;
  Tuple<int,CH_SPACEDIM-1> tangenDir = PolyGeom::computeTanDirs(a_faceDir);
  Real anormNorm = Abs(normal[a_faceDir]);
  Real anormTan[2];
  int signNorm = 1;
  if (normal[a_faceDir] < 0.0) signNorm = -1;
  int signTang[2];
  for (int itan = 0; itan < 2; itan++)
    {
      int tanDir = tangenDir[itan];
      anormTan[itan] = Abs(normal[tanDir]);
      signTang[itan] =  1;
      if (anormTan[itan] == 0.0)
        {//since signTang is arbitrary for this case,
          //check if we are next to the high domain edge
          IntVect ivHi = a_vof.gridIndex();;
          ivHi[tanDir] += 1;
          if (!m_domain.contains(ivHi))
            {
              signTang[itan] = -1;
            }
        }
      else if (normal[tanDir] < 0.0)
        {
          signTang[itan] = -1;
        }
    }

  const IntVect& ivVoF= a_vof.gridIndex();

  //whice one of the tangential directions has the largest normal
  int d1, d2;
  if (anormTan[0]/m_dx[tangenDir[0]] > anormTan[1]/m_dx[tangenDir[1]])
    {
      d1 = 0;
      d2 = 1;
    }
  else
    {
      d1 = 1;
      d2 = 0;
    }

  // figure out in which plane we are extrapolating
  bool faceDirOut = ((anormNorm/m_dx[a_faceDir] > anormTan[0]/m_dx[tangenDir[0]]) &&
                     (anormNorm/m_dx[a_faceDir] > anormTan[1]/m_dx[tangenDir[1]]));

  IntVect ivSten[2][2];
  bool hasVoF[2][2];
  VolIndex vofSten[2][2];

  if (faceDirOut)
    {// face direction has largest normal.

      ivSten[0][0] = ivVoF + signNorm*BASISV(a_faceDir);
      ivSten[1][0] = ivVoF + signTang[0]*BASISV(tangenDir[0]);
      ivSten[0][1] = ivVoF + signTang[1]*BASISV(tangenDir[1]);
      ivSten[1][1] = ivVoF + signTang[0]*BASISV(tangenDir[0]) + signTang[1]*BASISV(tangenDir[1]) ;
      d1 = 0;
      d2 = 1;

    }
  else
    { //tandir[d1] is the biggest

      ivSten[0][0] = ivVoF + signTang[d1]*BASISV(tangenDir[d1])                                      - signNorm*BASISV(a_faceDir);
      ivSten[1][0] = ivVoF + signTang[d1]*BASISV(tangenDir[d1])                                      ;
      ivSten[0][1] = ivVoF + signTang[d1]*BASISV(tangenDir[d1]) + signTang[d2]*BASISV(tangenDir[d2]) - signNorm*BASISV(a_faceDir);
      ivSten[1][1] = ivVoF + signTang[d1]*BASISV(tangenDir[d1]) + signTang[d2]*BASISV(tangenDir[d2]) ;

    }

  int radius = 1;
  Vector<VolIndex> vofsStencil;
  EBArith::getAllVoFsInMonotonePath(vofsStencil, a_vof, m_ebisBox, radius);

  for (int ix = 0; ix < 2; ix++)
    {
      for (int iy = 0; iy < 2; iy++)
        {
          hasVoF[ix][iy] =  EBArith::isVoFHere(vofSten[ix][iy], vofsStencil, ivSten[ix][iy]);
        }
    }
  bool hasAllVoFs = hasVoF[0][0] && hasVoF[1][0] && hasVoF[0][1] && hasVoF[1][1];

  //get the vof info for a 1D extrap stencil (if needed)
  VolIndex vofSten1D;
  const IntVect ivSten1D = ivVoF + signNorm*BASISV(a_faceDir);
  const bool has1DVoF =  EBArith::isVoFHere(vofSten1D, vofsStencil, ivSten1D);

  for (int ivar = 0; ivar < a_numPrim; ivar++)
    {
      if (hasAllVoFs)
        {
          Real WVal[2][2];
          for (int ix = 0; ix < 2; ix++)
            {
              for (int iy = 0; iy < 2; iy++)
                {
                  if (a_sd == Side::Hi)
                    {
                      WVal[ix][iy] = a_primMinu(vofSten[ix][iy], ivar);
                    }
                  else
                    {
                      WVal[ix][iy] = a_primPlus(vofSten[ix][iy], ivar);
                    }
                }

            }

          Real delWNorm[2][2];
          Real delWTan0[2][2];
          Real delWTan1[2][2];
          for (int ix = 0; ix < 2; ix++)
            {
              for (int iy = 0; iy < 2 ; iy++)
                {
                  coveredExtrapSlopes(delWNorm[ix][iy], vofSten[ix][iy], a_primState, a_faceDir   , ivar);
                  coveredExtrapSlopes(delWTan0[ix][iy], vofSten[ix][iy], a_primState, tangenDir[d1], ivar);
                  coveredExtrapSlopes(delWTan1[ix][iy], vofSten[ix][iy], a_primState, tangenDir[d2], ivar);
                }
            }

          if (faceDirOut)
            {
              // case 00 is special because we have to move it backwards in face direction to get it to plane
              Real deltaW;
              coveredExtrapSlopes(deltaW, vofSten[0][0],a_primState, a_faceDir, ivar);
              //Real slipSlope = a_slopesSeco[a_faceDir](vofSten[0][0], ivar);
              Real slipSlope = deltaW;
              WVal[0][0] -= signNorm*slipSlope;
              //debug to make compatible with 2d
              //WVal[0][1] = WVal[0][0];
              //end debug

              Real x0Val = (m_dx[a_faceDir]*anormTan[0])/(m_dx[tangenDir[0]]*anormNorm);
              Real x1Val = (m_dx[a_faceDir]*anormTan[1])/(m_dx[tangenDir[1]]*anormNorm);
              Real    funcVal = bilinearFunc(    WVal,x0Val,x1Val);
              Real deltaWNorm = bilinearFunc(delWNorm,x0Val,x1Val);
              Real deltaWTan0 = bilinearFunc(delWTan0,x0Val,x1Val);
              Real deltaWTan1 = bilinearFunc(delWTan1,x0Val,x1Val);
              Real dWdn =  -signNorm*deltaWNorm - signTang[d2]*x1Val*deltaWTan1 - signTang[d1]*x0Val*deltaWTan0;
              a_extrapVal[ivar] = funcVal + dWdn;

              if (usesFlattening())
                {
                  int pindex = pressureIndex();
                  int dindex = densityIndex();
                  if ((ivar == dindex) && (a_extrapVal[ivar] < 0.5*funcVal))
                    {
                      dropExtrap = true;
                    }

                  if ((ivar == pindex) || (ivar == dindex))
                    {
                      //if extrapolated to a negative density or pressure then turn
                      //off the slopes of the extrapolation
                      if (a_extrapVal[ivar] < 1.0e-8)
                        {
                          a_extrapVal[ivar] = funcVal;
                        }
                    }
                  if (dropExtrap)
                    {
                      a_extrapVal[ivar] = funcVal;
                    }
                }
            }
          else
            {
              //tanDir[d1] the biggestnormal
              Real x0Val = (m_dx[tangenDir[d1]]*anormNorm   )/(m_dx[    a_faceDir]*anormTan[d1]);
              Real x1Val = (m_dx[tangenDir[d1]]*anormTan[d2])/(m_dx[tangenDir[d2]]*anormTan[d1]);
              Real    funcVal = bilinearFunc(WVal,    x0Val,x1Val);
              Real deltaWNorm = bilinearFunc(delWNorm,x0Val,x1Val);
              Real deltaWTan0 = bilinearFunc(delWTan0,x0Val,x1Val);
              Real deltaWTan1 = bilinearFunc(delWTan1,x0Val,x1Val);

              Real dWdn =  -signNorm*x0Val*deltaWNorm - signTang[d2]*x1Val*deltaWTan1 - signTang[d1]*deltaWTan0;
              a_extrapVal[ivar] = funcVal + dWdn;

              if (usesFlattening())
                {
                  int pindex = pressureIndex();
                  int dindex = densityIndex();
                  if ((ivar == dindex) && (a_extrapVal[ivar] < 0.5*funcVal))
                    {
                      dropExtrap = true;
                    }
                  if ((ivar == pindex) || (ivar == dindex))
                    {
                      //if extrapolated to a negative density or pressure then turn
                      //off the slopes of the extrapolation
                      if (a_extrapVal[ivar] < 0.0)
                        {
                          a_extrapVal[ivar] = funcVal;
                        }
                    }
                  if (dropExtrap)
                    {
                      a_extrapVal[ivar] = funcVal;
                    }
                }
            }
        }
      else if (has1DVoF)
        {//try 1D extrap
          Real WVal1D;
          if (a_sd == Side::Hi)
            {
              WVal1D = a_primMinu(vofSten1D, ivar);
            }
          else
            {
              WVal1D = a_primPlus(vofSten1D, ivar);
            }
          Real deltaW;
          coveredExtrapSlopes(deltaW,vofSten1D,a_primState,a_faceDir,ivar);
          a_extrapVal[ivar] = WVal1D - 2.0*signNorm*deltaW;
        }
      else
        // At least one of the vofs in the stencil exists - we drop order
        // by using a weighted sum of the the wVal's.
        {
          a_extrapVal[ivar] = 0.;
          Real volTot = 0.;
          bool hasAVoF = false;
          for (int ix = 0; ix < 2; ix++)
            {
              for (int iy = 0; iy < 2; iy++)
                {
                  if (hasVoF[ix][iy])
                    {
                      hasAVoF = true;
                      Real wVal;
                      if (a_sd == Side::Hi)
                        {
                          wVal = a_primMinu(vofSten[ix][iy], ivar);
                        }
                      else
                        {
                          wVal = a_primPlus(vofSten[ix][iy], ivar);
                        }
                      Real volFrac = m_ebisBox.volFrac(vofSten[ix][iy]);
                      a_extrapVal[ivar] += wVal*volFrac;
                      volTot += volFrac;
                    }
                }
            }
          if (hasAVoF)
            {
              CH_assert(volTot > 0.0);
              a_extrapVal[ivar] = a_extrapVal[ivar]/volTot;
            }
          else
            {
              a_extrapVal[ivar] = a_primState(a_vof, ivar);
            }
        }
      // None of the vofs in the stencil exists. We use the value at the
      // cell center adjacent to the face.
    }
  //    if (!hasAllVoFs)
  //          cout << "Dropping order at covered face = " << a_vof <<
  //          "direction = " << a_faceDir <<  "\n";
}
/******/
void
EBPatchReactive::floorPrimitives(EBCellFAB&  a_primState,
                                   const Box&  a_box)
{
  CH_TIME("EBPatchReactive::floorPrimitives");
  BaseFab<Real>&       primReg = a_primState.getSingleValuedFAB();

  FORT_FLOORPRIM( CHF_BOX(a_box),
                  CHF_FRA(primReg));
  Real small, smallp, smallu, smallr;
  FORT_GETSMALL(CHF_REAL(small),
                CHF_REAL(smallp),
                CHF_REAL(smallu),
                CHF_REAL(smallr));
  for (int ivof = 0; ivof<m_irregVoFs.size(); ivof++)
    {
      const VolIndex& vof = m_irregVoFs[ivof];
      if (a_box.contains(vof.gridIndex()))
        {
          a_primState(vof, QRHO)  = Max(a_primState(vof, QRHO) ,  smallr);
          a_primState(vof, QPRES) = Max(a_primState(vof, QPRES),  smallp);
          a_primState(vof, QTEMP) = Max(a_primState(vof, QTEMP),  small);
         }
     }
}
/******/
void
EBPatchReactive::floorConserved(EBCellFAB&  a_consState,
                                  const Box&  a_box)
{
  CH_TIME("EBPatchReactive::floorConserved");
  BaseFab<Real>&       consReg = a_consState.getSingleValuedFAB();

  FORT_FLOORCONS( CHF_BOX(a_box),
                  CHF_FRA(consReg));

  Real small, smallp, smallu, smallr;
  FORT_GETSMALL(CHF_REAL(small),
                CHF_REAL(smallp),
                CHF_REAL(smallu),
                CHF_REAL(smallr));
  IntVectSet ivsIrreg = m_ebisBox.getIrregIVS(a_box);
  for (VoFIterator vofit(ivsIrreg, m_ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      a_consState(vof, CRHO) = Max(a_consState(vof, CRHO), smallr);
      a_consState(vof, CENG) = Max(a_consState(vof, CENG), small );
    }
}
/******/
void
EBPatchReactive::floorConserved(BaseIVFAB<Real>& a_consState,
                                  const IntVectSet& a_ivsIrreg)
{
  CH_TIME("EBPatchReactive::floorConservedIrr");
  Real small, smallp, smallu, smallr;
  FORT_GETSMALL(CHF_REAL(small),
                CHF_REAL(smallp),
                CHF_REAL(smallu),
                CHF_REAL(smallr));
  for (VoFIterator vofit(a_ivsIrreg, m_ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      a_consState(vof, CRHO) = Max(a_consState(vof, CRHO), smallr);
      a_consState(vof, CENG) = Max(a_consState(vof, CENG), small );
    }

}
/******/
void
EBPatchReactive::floorPrimitives(BaseIVFAB<Real>& a_primState,
                                   const IntVectSet& a_ivsIrreg)
{
  CH_TIME("EBPatchReactive::floorPrimitivesIrr");
  Real small, smallp, smallu, smallr;
  FORT_GETSMALL(CHF_REAL(small),
                CHF_REAL(smallp),
                CHF_REAL(smallu),
                CHF_REAL(smallr));
  for (VoFIterator vofit(a_ivsIrreg, m_ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      //floors
      a_primState(vof, QRHO)  = Max(a_primState(vof, QRHO) ,  smallr);
      a_primState(vof, QPRES) = Max(a_primState(vof, QPRES),  smallp);
    }
}
/*****************************/
void
EBPatchReactive::
finalExtrap2D(EBCellFAB              a_primMinu[SpaceDim],
              EBCellFAB              a_primPlus[SpaceDim],
              const BaseIVFAB<Real>  a_coveredFluxNormMinu[SpaceDim],
              const BaseIVFAB<Real>  a_coveredFluxNormPlus[SpaceDim],
              const Vector<VolIndex> a_coveredFaceNormMinu[SpaceDim],
              const Vector<VolIndex> a_coveredFaceNormPlus[SpaceDim],
              const EBFaceFAB        a_fluxOne[SpaceDim],
              const EBCellFAB&       a_primState,
              const EBCellFAB        a_slopesPrim[SpaceDim],
              const EBCellFAB        a_slopesSeco[SpaceDim],
              const Box&             a_box)

{
  CH_TIME("EBPatchReactive::finalExtrap2D");
  Box slopeBoxG1 = grow(a_box, 1);
  slopeBoxG1 &= m_domain;

  //7. compute fnal corrections to W due to the
  //final transverse derivatives.
  // Do the final corrections to the fluxes
  for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      // Correct the flux using fluxes in the remaining direction(s)
      for (int diffDir = 0; diffDir < SpaceDim; diffDir++)
        {
          // A different direction has been found
          if (diffDir != faceDir)
            {

              // In 2D,
              //the current primitive state is updated by a flux in
              // the other direction
              //a_coveredFlux holds the covered flux
              //due to the normal derivative
              //equation 1.18 line 3.
              updatePrim(a_primMinu[faceDir],
                         a_primPlus[faceDir], a_fluxOne[diffDir],
                         a_coveredFluxNormMinu[diffDir], a_coveredFluxNormPlus[diffDir],
                         a_coveredFaceNormMinu[diffDir], a_coveredFaceNormPlus[diffDir],
                         diffDir,  slopeBoxG1, (0.5)*m_dt/m_dx[diffDir]);

            } //end if dir2 != faceDir
        } //loop over dir2
    } //loop over facedir
}
/***********/
void
EBPatchReactive::
updatePrim(EBCellFAB&              a_primMinu,
           EBCellFAB&              a_primPlus,
           const EBFaceFAB&        a_flux,
           const BaseIVFAB<Real>&  a_coveredFluxMinu,
           const BaseIVFAB<Real>&  a_coveredFluxPlus,
           const Vector<VolIndex>& a_coveredFaceMinu,
           const Vector<VolIndex>& a_coveredFacePlus,
           const int&              a_dir,
           const Box&              a_box,
           const Real&             a_scale)
{
  CH_TIME("EBPatchReactive::updatePrim");
  CH_assert(isDefined());
  CH_assert(a_primPlus.nComp() == numPrimitives());
  CH_assert(a_primMinu.nComp() == numPrimitives());
  CH_assert(a_flux.nComp() == numFluxes());
  CH_assert((a_dir >= 0) && (a_dir < SpaceDim));
  CH_assert(a_flux.direction() == a_dir);
  CH_assert(a_flux.getRegion().contains(surroundingNodes(a_box, a_dir)));
  CH_assert(a_primPlus.getRegion().contains(a_box));
  CH_assert(a_primMinu.getRegion().contains(a_box));
  CH_assert(a_primPlus.getRegion() == a_primMinu.getRegion());

  const Box& primBox = a_primPlus.getRegion();
  int nCons =  numConserved();
  EBCellFAB consTemp(m_ebisBox, m_validBoxG4, nCons);
  primToCons(consTemp, a_primMinu, primBox);

  updateCons(consTemp, a_flux,
             a_coveredFluxMinu, a_coveredFluxPlus,
             a_coveredFaceMinu, a_coveredFacePlus,
             a_dir, a_box, a_scale);
  int logflag = 0;

//  pout() << "consTemp after updateCons from primMinu" << endl;
  consToPrim(a_primMinu, consTemp, m_validBoxG4, logflag);

  primToCons(consTemp, a_primPlus, primBox);
  updateCons(consTemp, a_flux,
             a_coveredFluxMinu, a_coveredFluxPlus,
             a_coveredFaceMinu, a_coveredFacePlus,
             a_dir, a_box, a_scale);
  
//  pout() << "consTemp after updateCons from primMinu" << endl;
  consToPrim(a_primPlus, consTemp, primBox, logflag);

}
/*****/
void
EBPatchReactive::
updateCons(EBCellFAB&              a_consState,
           const EBFaceFAB&        a_flux,
           const BaseIVFAB<Real>&  a_coveredFluxMinu,
           const BaseIVFAB<Real>&  a_coveredFluxPlus,
           const Vector<VolIndex>& a_coveredFaceMinu,
           const Vector<VolIndex>& a_coveredFacePlus,
           const int&              a_dir,
           const Box&              a_box,
           const Real&             a_scale)
{
  CH_TIME("EBPatchReactive::updateCons");
  CH_assert(isDefined());
  CH_assert(a_consState.nComp() <= a_flux.nComp());
  CH_assert(a_consState.nComp() == numConserved() );
  CH_assert((a_dir >= 0) && (a_dir < SpaceDim));
  CH_assert(a_flux.direction() == a_dir);
  CH_assert(a_flux.getRegion().contains(surroundingNodes(a_box, a_dir)));
  CH_assert(a_consState.getRegion().contains(a_box));

  Vector<Vector<Real> > cache;
  {
    CH_TIME("update cons cache");
    cacheEBCF(cache, a_consState);
  }
  {
    CH_TIME("EBPatchReactive::updateConsRegular");
    int ncons = numConserved();
    BaseFab<Real>&       consReg = a_consState.getSingleValuedFAB();
    const BaseFab<Real>& fluxReg = a_flux.getSingleValuedFAB();
/*
    // debug
    pout() << "Reg flux in updateCons" << endl;
    FabDataOps::getFabData(fluxReg);
*/
    FORT_UPDATE( CHF_BOX(a_box),
                 CHF_FRA(consReg),
                 CHF_CONST_FRA(fluxReg),
                 CHF_CONST_INT(a_dir),
                 CHF_CONST_INT(ncons),
                 CHF_CONST_REAL(a_scale));
  }

  {
    CH_TIME("update cons uncache");
    uncacheEBCF(a_consState, cache);
  }
  //update the irregular vofs
  {
    CH_TIME("EBPatchReactive::updateConsIrregular");
    for (int ivof = 0; ivof<m_irregVoFs.size(); ivof++)
      {
        const VolIndex& vof = m_irregVoFs[ivof];
        if (a_box.contains(vof.gridIndex()))
          {
            for (SideIterator sit; sit.ok(); ++sit)
              {
                int isign = sign(sit());
                Vector<FaceIndex> faces = m_ebisBox.getFaces(vof, a_dir, sit());
                for (int ivar = 0; ivar < a_consState.nComp(); ivar++)
                  {
                    //aggregate the fluxes over the faces for a_dir and sit()
                    Real aggFlux = 0.0;
                    for (int iface = 0; iface < faces.size(); iface++)
                      {
                        const FaceIndex& face = faces[iface];
                        Real faceFlux = a_flux(face, ivar);
                        aggFlux += faceFlux;
                      }

                    if (faces.size() > 1)
                      {
                        aggFlux /= Real(faces.size());
                      }

                    //dx is already divided out (part of scale)
                    aggFlux *= -isign*a_scale;

                     //a_consState modified in regular update, but restored in uncache
                    a_consState(vof, ivar) += aggFlux;
                  }
              }
          }
      }
  }
  {
    CH_TIME("EBPatchReactive::updateCons::Covered");
    //add in covered face fluxes
    for (SideIterator sit; sit.ok(); ++sit)
      {
        int isign = sign(sit());
        const BaseIVFAB<Real>*  coveredFluxPtr;
        const Vector<VolIndex>* coveredFacePtr;
        if (sit() == Side::Lo)
          {
            coveredFluxPtr = &(a_coveredFluxMinu);
            coveredFacePtr = &(a_coveredFaceMinu);
          }
        else
          {
            coveredFluxPtr = &(a_coveredFluxPlus);
            coveredFacePtr = &(a_coveredFacePlus);
          }
        const BaseIVFAB<Real>&  coveredFlux = *coveredFluxPtr;
        const Vector<VolIndex>& coveredFace = *coveredFacePtr;
        for (int ivof = 0; ivof < coveredFace.size(); ivof++)
          {
            const VolIndex& vof = coveredFace[ivof];
            //the input set can be bigger than the box over which we are integrating.
            //this was done for peformance reasonas so the sets to not get
            //calculated over and over.
            if (a_box.contains(vof.gridIndex()))
              {
                for (int ivar = 0; ivar < a_consState.nComp(); ivar++)
                  {

                    Real covFlux = coveredFlux(vof, ivar);
                    Real update =   isign*covFlux;

                    //dx is already divided out (part of scale)
                    update *= -a_scale;
                    Real state = a_consState(vof, ivar);
                    a_consState(vof, ivar) = state + update;

                  }
              }
          }
      }
  }

  {
    CH_TIME("EBPatchReactive::updateCons::floors");
    floorConserved(a_consState, a_box);
  }
}
/******/
void
EBPatchReactive::
cacheEBCF( Vector<Vector<Real> >& a_cache, const  EBCellFAB& a_input)
{
  CH_assert(a_input.getRegion() == m_validBoxG4);

  int nvar = a_input.nComp();
  a_cache.resize(nvar);
  for (int ivar = 0; ivar < nvar; ivar++)
    {
      a_cache[ivar].resize(m_irregVoFs.size());
    }
  for (int ivar = 0; ivar < nvar; ivar++)
    {
      const Real* inputSV = a_input.getSingleValuedFAB().dataPtr(ivar);
      const Real* inputMV =  a_input.getMultiValuedFAB().dataPtr(ivar);
      for (int ivof = 0; ivof<m_irregVoFs.size(); ivof++)
        {
          pointerOffset_t& vofOff = m_updateStencil[ivof].m_vofOffset;
          const Real* minuPtr;
          if (vofOff.m_multiValued)
            {
              minuPtr = inputMV + vofOff.m_offset;
            }
          else
            {
              minuPtr = inputSV + vofOff.m_offset;
            }

          a_cache[ivar][ivof] = *minuPtr;
        }
    }
}
/******************/
void
EBPatchReactive::
uncacheEBCF(EBCellFAB& a_output, const Vector<Vector<Real> >& a_cache)
{
  CH_assert(a_output.getRegion() == m_validBoxG4);
  CH_assert(a_output.nComp()     == a_cache.size());

  for (int ivar = 0; ivar < a_cache.size(); ivar++)
    {
      Real* outputSV = a_output.getSingleValuedFAB().dataPtr(ivar);
      Real* outputMV = a_output.getMultiValuedFAB().dataPtr(ivar);

      for (int ivof = 0; ivof<m_irregVoFs.size(); ivof++)
        {
          //const VolIndex& vof = m_irregVoFs[ivof];
          pointerOffset_t& vofOff = m_updateStencil[ivof].m_vofOffset;
          Real* minuPtr;
          if (vofOff.m_multiValued)
            {
              minuPtr = outputMV + vofOff.m_offset;
            }
          else
            {
              minuPtr = outputSV + vofOff.m_offset;
            }
          *minuPtr = a_cache[ivar][ivof];
        }
    }
}
/******************/
void
EBPatchReactive::
extrapolatePrim3D(EBCellFAB           a_primMinu[SpaceDim],
                  EBCellFAB           a_primPlus[SpaceDim],
                  EBCellFAB&          a_primState,
                  EBCellFAB           a_slopesPrim[SpaceDim],
                  EBCellFAB           a_slopesSeco[SpaceDim],
                  const EBCellFAB&    a_flattening,
                  const EBCellFAB&    a_consState,
                  const EBCellFAB&    a_source,
                  const Box&          a_box,
                  const DataIndex&    a_dit,
                  bool                a_verbose)
{
  CH_TIME("EBPatchReactive::extrapolatePrim3D");
  //define the plethora of data holders that I need

  //now do the actual computation.
  //1. transform to primitive state
  int numPrim = numPrimitives();
  Box primBox = a_consState.getRegion();
  CH_assert(m_domain.contains(primBox));
  a_primState.define(m_ebisBox, primBox, numPrim);
  int logflag = 0;
  consToPrim(a_primState, a_consState, primBox, logflag);

  //2.  compute slopes delta^d w
  //compute normal derivative stuff
  doNormalDerivativeExtr3D(a_primMinu,
                           a_primPlus,
                           m_fluxOne,
                           m_coveredFluxNormMinu,
                           m_coveredFluxNormPlus,
                           m_coveredFaceMinuG4,
                           m_coveredFacePlusG4,
                           a_slopesPrim,
                           a_slopesSeco,
                           a_flattening,
                           a_primState,
                           a_source,
                           a_dit,
                           a_box);

  //6. In 3D compute corrections to u corresponding to one set
  //   of transverse derivatives appropriate to obtain (1, 1, 1)
  //   coupling.
  // In 3D, compute some additional intermediate fluxes
  // NOTE:  The diagonal entries of this array of fluxes are not
  // used and will not be defined.

  do111coupling(m_fluxTwo,
                m_coveredFluxMinu3D,
                m_coveredFluxPlus3D,
                a_primMinu,
                a_primPlus,
                m_coveredFluxNormMinu,
                m_coveredFluxNormPlus,
                m_coveredFaceMinuG4,
                m_coveredFacePlusG4,
                m_fluxOne,
                a_primState,
                a_slopesPrim,
                a_slopesSeco,
                a_dit,
                a_box);

  finalExtrap3D(a_primMinu,
                a_primPlus,
                m_coveredFluxMinu3D,
                m_coveredFluxPlus3D,
                m_fluxTwo,
                a_primState,
                a_slopesPrim,
                a_slopesSeco,
                a_box);
}
/******************/
void EBPatchReactive::
doNormalDerivativeExtr3D(EBCellFAB              a_primMinu[SpaceDim],
                         EBCellFAB              a_primPlus[SpaceDim],
                         EBFaceFAB              a_fluxOne[SpaceDim],
                         BaseIVFAB<Real>        a_coveredFluxNormMinu[SpaceDim],
                         BaseIVFAB<Real>        a_coveredFluxNormPlus[SpaceDim],
                         Vector<VolIndex>       a_coveredFaceNormMinu[SpaceDim],
                         Vector<VolIndex>       a_coveredFaceNormPlus[SpaceDim],
                         EBCellFAB              a_slopesPrim[SpaceDim],
                         EBCellFAB              a_slopesSeco[SpaceDim],
                         const EBCellFAB&       a_flattening,
                         const EBCellFAB&       a_primState,
                         const EBCellFAB&       a_source,
                         const DataIndex&       a_dit,
                         const Box&             a_box)

{
  CH_TIME("EBPatchReactive::doNormalDerivative3D");
  Box modBoxCov[SpaceDim];
  Box modBoxOpen[SpaceDim];
  Box faceBox[SpaceDim];
  Box bndryFaceBox[SpaceDim];
  int numSlop = numSlopes();

  //set up the data structures
  for (int idir = 0; idir < SpaceDim; ++idir)
    {
      modBoxOpen[idir] = a_box;
      modBoxOpen[idir].grow(3);
      modBoxOpen[idir].grow(idir, -1);
      modBoxOpen[idir] &= m_domain;

      modBoxCov[idir] = a_box;
      modBoxCov[idir].grow(2);
      modBoxCov[idir].grow(idir, -1);
      modBoxCov[idir] &= m_domain;

      bndryFaceBox[idir] = modBoxCov[idir];
      bndryFaceBox[idir].surroundingNodes(idir);

      faceBox[idir] = modBoxCov[idir];
      faceBox[idir].grow(idir, 1);
      faceBox[idir] &= m_domain;
      faceBox[idir].grow(idir, -1);
      faceBox[idir].surroundingNodes(idir);

      a_slopesPrim[idir].define(m_ebisBox, modBoxOpen[idir], numSlop);
      a_slopesSeco[idir].define(m_ebisBox, modBoxOpen[idir], numSlop);

      m_extendStateNormPlus[idir].setVal(0.);
      m_extendStateNormMinu[idir].setVal(0.);
      a_coveredFluxNormPlus[idir].setVal(0.);
      a_coveredFluxNormMinu[idir].setVal(0.);
      a_primPlus[idir].setVal(0.);
      a_primMinu[idir].setVal(0.);
      a_fluxOne[idir].setVal(0.);
    }
  //all one big loop since no cross terms here
  for (int idir = 0; idir < SpaceDim; ++idir)
    {
      //check flattening coefficients have correct box and compute slopes
      if (usesFlattening())
        CH_assert(a_flattening.getRegion().contains(modBoxOpen[idir]));

      //compute slopes
      slope(a_slopesPrim[idir], a_slopesSeco[idir], a_primState,
            a_flattening, idir, modBoxOpen[idir]);

      //3. Compute the effect of the normal derivative terms and the
      //source term on the extrapolation in space and time from cell
      //centers to faces. Equation 1.7 (top).
      normalPred(a_primMinu[idir], a_primPlus[idir],
                 a_primState, a_slopesPrim[idir],
                 m_dt/m_dx[idir], idir, modBoxOpen[idir]);

      // If the source term is valid add it to the primitive quantities
      //if source term is not defined, assume it is zero.  lots of
      //problems will have no sources.  not the most elegant method
      //but it beats passing a null pointer.
      if (a_source.isDefined())
        {
          if (!s_conservativeSource)
            {
              incrementWithSource(a_primMinu[idir], a_source,
                                  0.5*m_dt, modBoxOpen[idir]);
              incrementWithSource(a_primPlus[idir], a_source,
                                  0.5*m_dt, modBoxOpen[idir]);
            }
          else
            {
              int logflag = 0;
              const Box& modBox = modBoxOpen[idir];
              EBCellFAB  consTemp(m_ebisBox, modBox, numConserved());

              primToCons(consTemp,  a_primMinu[idir],  modBox);
              incrementWithSource(consTemp,  a_source, 0.5*m_dt,  modBox); //
              consToPrim(a_primMinu[idir],  consTemp,    modBox, logflag);

              primToCons(consTemp,  a_primPlus[idir],  modBox);
              incrementWithSource(consTemp,  a_source,  0.5*m_dt, modBox); //
              consToPrim(a_primPlus[idir],  consTemp,   modBox, logflag);
            }

        }

      //4. compute estimates of the flux suitable for computing
      // 1D flux derivatives using a riemann solver for the interior R
      //and for the boundary RB.  Equation 1.8
/*
      // debug 
      pout() << "PrimPlus" << endl;
      FabDataOps::getFabData(a_primPlus[idir],1);
      pout() << "PrimMinu" << endl;
      FabDataOps::getFabData(a_primMinu[idir],0);
      // end debug
*/
      riemann(a_fluxOne[idir], a_primPlus[idir], a_primMinu[idir],
              idir, faceBox[idir]);

      //some wackiness to work around the fact that fluxfab did
      //not exist when a lot of this was written.
      Vector<EBFaceFAB*> bcfluxes(SpaceDim);
      for (int kdir = 0; kdir < SpaceDim; kdir++)
        {
          bcfluxes[kdir] = &(a_fluxOne[kdir]);
        }
      EBFluxFAB fluxAlias;
      fluxAlias.alias(bcfluxes);

      // Use the user supplied PhysBC object to obtain boundary fluxes
      m_bc->fluxBC(fluxAlias, a_primState, a_primMinu[idir],
                   Side::Lo, m_time, m_ebisBox, a_dit,a_box, bndryFaceBox[idir], idir);
      m_bc->fluxBC(fluxAlias, a_primState, a_primPlus[idir],
                   Side::Hi, m_time, m_ebisBox, a_dit,a_box, bndryFaceBox[idir], idir);

      //extrapolate to covered faces
      extrapToCoveredFaces(m_extendStateNormMinu[idir],
                           a_primMinu[idir],
                           a_primPlus[idir],
                           a_primState,
                           a_coveredFaceNormMinu[idir],
                           idir, Side::Lo, modBoxCov[idir]);

      extrapToCoveredFaces(m_extendStateNormPlus[idir],
                           a_primMinu[idir],
                           a_primPlus[idir],
                           a_primState,
                           a_coveredFaceNormPlus[idir],
                           idir, Side::Hi, modBoxCov[idir]);

      //solve riemann problem between extended state and
      //value in vof to get covered flux  Equation 1.9
      riemann(a_coveredFluxNormMinu[idir],
              m_extendStateNormMinu[idir], a_primMinu[idir],
              a_coveredFaceNormMinu[idir], idir, Side::Lo, modBoxOpen[idir]);
      riemann(a_coveredFluxNormPlus[idir],
              m_extendStateNormPlus[idir], a_primPlus[idir],
              a_coveredFaceNormPlus[idir], idir, Side::Hi, modBoxOpen[idir]);
    }
}

/*******************/
void EBPatchReactive::
do111coupling(EBFaceFAB              a_fluxTwo[SpaceDim][SpaceDim],
              BaseIVFAB<Real>        a_coveredFluxMinu3D[SpaceDim][SpaceDim],
              BaseIVFAB<Real>        a_coveredFluxPlus3D[SpaceDim][SpaceDim],
              const EBCellFAB        a_primMinu[SpaceDim],
              const EBCellFAB        a_primPlus[SpaceDim],
              const BaseIVFAB<Real>  a_coveredFluxNormMinu[SpaceDim],
              const BaseIVFAB<Real>  a_coveredFluxNormPlus[SpaceDim],
              const Vector<VolIndex> a_coveredFaceNormMinu[SpaceDim],
              const Vector<VolIndex> a_coveredFaceNormPlus[SpaceDim],
              const EBFaceFAB        a_fluxOne[SpaceDim],
              const EBCellFAB&       a_primState,
              const EBCellFAB        a_slopesPrim[SpaceDim],
              const EBCellFAB        a_slopesSeco[SpaceDim],
              const DataIndex&       a_dit,
              const Box&             a_box)
{
  CH_TIME("EBPatchReactive::do111coupling");
  Box slopeBoxG1 = grow(a_box, 1);
  Box slopeBoxG2 = grow(a_box, 2);
  slopeBoxG1 &= m_domain;
  slopeBoxG2 &= m_domain;
  Box faceBox[SpaceDim][SpaceDim];
  Box bndryFaceBox[SpaceDim][SpaceDim];
  Box modBoxCov[SpaceDim];
  Box modBoxOpen[SpaceDim][SpaceDim];
  for (int faceDir = 0; faceDir < SpaceDim; ++faceDir)
    {
      for (int diffDir = 0; diffDir < SpaceDim; ++diffDir)
        {

          modBoxOpen[faceDir][diffDir] = grow(a_box, 2);
          modBoxOpen[faceDir][diffDir].grow(diffDir, -1);
          modBoxOpen[faceDir][diffDir] &= m_domain;
        }
    }

  for (int modDir = 0; modDir < SpaceDim; ++modDir)
    {
      modBoxCov[modDir] = grow(a_box, 1);
      modBoxCov[modDir] &= m_domain;
      for (int faceDir = 0; faceDir < SpaceDim; ++faceDir)
        {
          bndryFaceBox[faceDir][modDir] = modBoxCov[modDir];
          bndryFaceBox[faceDir][modDir].surroundingNodes(faceDir);

          faceBox[faceDir][modDir] = modBoxCov[modDir];
          faceBox[faceDir][modDir].grow(faceDir, 1);
          faceBox[faceDir][modDir] &= m_domain;
          faceBox[faceDir][modDir].grow(faceDir, -1);
          faceBox[faceDir][modDir].surroundingNodes(faceDir);
        }
    }

  for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      for (int diffDir = 0; diffDir < SpaceDim; diffDir++)
        {
          if (diffDir != faceDir)
            {
              a_coveredFluxPlus3D[faceDir][diffDir].setVal(0.);
              a_coveredFluxMinu3D[faceDir][diffDir].setVal(0.);
              a_fluxTwo[faceDir][diffDir].setVal(0.);
              m_extendStateMinu3D[faceDir][diffDir].setVal(0.);
              m_extendStatePlus3D[faceDir][diffDir].setVal(0.);
            }
        }
    }
  for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      for (int diffDir = 0; diffDir < SpaceDim; diffDir++)
        {
          if (diffDir != faceDir)
            {
              int modDir = 3 - faceDir - diffDir;

              // Copy data for in place modification
              m_primMinuTemp.copy(a_primMinu[faceDir]);
              m_primPlusTemp.copy(a_primPlus[faceDir]);

              // Update the current, extrapolated primitive
              //variable using a flux
              // in a different direction.  uses covered fluxes calculated
              // above.  Equation 1.12
              updatePrim(m_primMinuTemp, m_primPlusTemp, a_fluxOne[diffDir],
                         a_coveredFluxNormMinu[diffDir], a_coveredFluxNormPlus[diffDir],
                         a_coveredFaceNormMinu[diffDir], a_coveredFaceNormPlus[diffDir],
                         diffDir,  modBoxOpen[faceDir][diffDir],  (1.0/3.0)*m_dt/m_dx[diffDir]);

              extrapToCoveredFaces(m_extendStateMinu3D[faceDir][diffDir],
                                   m_primMinuTemp,
                                   m_primPlusTemp,
                                   a_primState,
                                   m_coveredFaceMinuG4[faceDir],
                                   faceDir, Side::Lo, modBoxCov[modDir]);

              extrapToCoveredFaces(m_extendStatePlus3D[faceDir][diffDir],
                                   m_primMinuTemp,
                                   m_primPlusTemp,
                                   a_primState,
                                   m_coveredFacePlusG4[faceDir],
                                   faceDir, Side::Hi, modBoxCov[modDir]);

              // Solve the Riemann problem and get fluxes.  Eqution 1.15
              riemann(a_fluxTwo[faceDir][diffDir], m_primPlusTemp, m_primMinuTemp,
                      faceDir, faceBox[faceDir][modDir]);

              //some wackiness to work around the fact that fluxfab did
              //not exist when a lot of this was written.
              Vector<EBFaceFAB*> bcfluxes(SpaceDim);
              for (int idir = 0; idir < SpaceDim; idir++)
                {
                  bcfluxes[idir] = &(a_fluxTwo[idir][diffDir]);
                }
              EBFluxFAB fluxAlias;
              fluxAlias.alias(bcfluxes);

              m_bc->fluxBC(fluxAlias,  a_primState, m_primMinuTemp,
                           Side::Lo,  m_time, m_ebisBox, a_dit,a_box, bndryFaceBox[faceDir][diffDir], faceDir);
              m_bc->fluxBC(fluxAlias,  a_primState, m_primPlusTemp,
                           Side::Hi,  m_time, m_ebisBox, a_dit,a_box, bndryFaceBox[faceDir][diffDir], faceDir);

              //solve riemann problem between extended state and
              //value in vof to get covered flux.  Equation 1.14
              riemann(a_coveredFluxMinu3D[faceDir][diffDir],
                      m_extendStateMinu3D[faceDir][diffDir], m_primMinuTemp,
                      m_coveredFaceMinuG4[faceDir], faceDir, Side::Lo, modBoxOpen[faceDir][diffDir]);
              riemann(a_coveredFluxPlus3D[faceDir][diffDir],
                      m_extendStatePlus3D[faceDir][diffDir], m_primPlusTemp,
                      m_coveredFacePlusG4[faceDir], faceDir, Side::Hi, modBoxOpen[faceDir][diffDir]);

            }
        }
    }
}
/*******************/
void
EBPatchReactive::
finalExtrap3D(EBCellFAB              a_primMinu[SpaceDim],
              EBCellFAB              a_primPlus[SpaceDim],
              const BaseIVFAB<Real>  a_coveredFluxMinu3D[SpaceDim][SpaceDim],
              const BaseIVFAB<Real>  a_coveredFluxPlus3D[SpaceDim][SpaceDim],
              const EBFaceFAB        a_fluxTwo[SpaceDim][SpaceDim],
              const EBCellFAB&       a_primState,
              const EBCellFAB        a_slopesPrim[SpaceDim],
              const EBCellFAB        a_slopesSeco[SpaceDim],
              const Box&             a_box)

{
  CH_TIME("EBPatchReactive::finalExtrap3D");
  Box slopeBoxG1 = grow(a_box, 1);
  slopeBoxG1 &= m_domain;

  //7. compute fnal corrections to W due to the
  //final transverse derivatives.
  // Do the final corrections to the fluxes
  for (int dir1 = 0; dir1 < SpaceDim; dir1++)
    {
      // Correct the flux using fluxes in the remaining direction(s)
      for (int dir2 = 0; dir2 < SpaceDim; dir2++)
        {
          // A different direction has been found
          if (dir2 != dir1)
            {
              // In 3D,  find a direction different from the two above
              // here we need to use the coveredFlux1D for the covered flux

              int dir3 = 3 - dir2  - dir1;
              // A different direction has been found
              //const EBFluxFAB& fluxTwo = a_fluxTwo[dir3];
              // Update the conservative state
              //using both corrected fluxes in
              // the other two directions.  Equation 1.18 line 4
              updatePrim(a_primMinu[dir1], a_primPlus[dir1],
                         a_fluxTwo[dir2][dir3],
                         a_coveredFluxMinu3D[dir2][dir3],
                         a_coveredFluxPlus3D[dir2][dir3],
                         m_coveredFaceMinuG4[dir2],
                         m_coveredFacePlusG4[dir2],
                         dir2,  slopeBoxG1,  (0.5)*m_dt/m_dx[dir2]);

            }
        } //loop over dir2

    } //loop over dir1

}
/*******************/
void
EBPatchReactive::
computeCoveredFaces(Vector<VolIndex>&     a_coveredFace,
                    IntVectSet&           a_coveredSets,
                    IntVectSet&           a_irregIVS,
                    const int&            a_idir,
                    const Side::LoHiSide& a_sd,
                    const Box&            a_region)
{
  CH_TIME("EBPatchReactive::computeCoveredFaces");

  EBArith::computeCoveredFaces(a_coveredFace, a_coveredSets,
                               a_irregIVS, a_idir, a_sd,
                               m_ebisBox, a_region);
}
/*****************************/
/*
void
EBPatchReactive::
extrapToCoveredFaces(BaseIVFAB<Real>&        a_extendedPrim,
                     const EBCellFAB&        a_primMinu,
                     const EBCellFAB&        a_primPlus,
                     const EBCellFAB&        a_primState,
                     const Vector<VolIndex>& a_coveredFaces,
                     const int&              a_faceDir,
                     const Side::LoHiSide&   a_sd,
                     const Box&              a_box)
{
  CH_TIME("EBPatchReactive::extrapToCoveredFaces");
  for (int ivof = 0; ivof < a_coveredFaces.size(); ivof++)
    {
      const VolIndex& vof = a_coveredFaces[ivof];

      RealVect normal = m_ebisBox.normal(vof);
      if (a_box.contains(vof.gridIndex()))
        {
          const int numPrim = numPrimitives();
          Vector<Real> extPrim(numPrim, 0.0);
          if (SpaceDim== 2)
            {
              pointExtrapToCovered2D(extPrim,
                                     a_primMinu,
                                     a_primPlus,
                                     a_primState,
                                     a_faceDir,
                                     vof,
                                     normal,
                                     a_sd,
                                     numPrim);
            }
          else if (SpaceDim==3)
            {
              pointExtrapToCovered3D(extPrim,
                                     a_primMinu,
                                     a_primPlus,
                                     a_primState,
                                     a_faceDir,
                                     vof,
                                     normal,
                                     a_sd,
                                     numPrim);
            }
          else
            {
              MayDay::Error("Bogus SpaceDim");
            }

          for (int ivar = 0; ivar < numPrim; ivar++)
            {
              a_extendedPrim(vof, ivar) = extPrim[ivar];
            }
        }
    }

  const IntVectSet& ivs  = a_extendedPrim.getIVS();
  floorPrimitives(a_extendedPrim, ivs);
}
*/
/********/
void EBPatchReactive::
getFaceDivergence(EBFluxFAB&        a_openDivU,
                  const EBCellFAB&  a_primState,
                  const EBCellFAB   a_slopePrim[SpaceDim],
                  const Box&        a_box,
                  const IntVectSet& a_ivsIrreg)
{
  CH_TIME("EBPatchReactive::getFaceDivergence");
  //first set open divergence to zero because
  //we will compute the divergence additively
  CH_assert(a_openDivU.nComp() == 1);
  for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      a_openDivU[faceDir].setVal(0.);
    }

  //compute divergence on all cells as though they
  //were regular.  Because clarity is not the only virtue
  //in life, we will do this in fortran
  for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      //this box calc stuff poached directly from
      //the non-eb version.  The fortran used here
      //is slightly different since i am using the
      //existing slopes instead of recalculating the
      //gradients.

      // Now, we can calculate the divergence of the normal velocity
      // at the center normal-direction edges. To do this, we determine
      // which edges at which we have sufficient data to compute centered
      // estimates of h*(div(u)). At the remaining edges. i.e. those
      // corresponding to the physical boundaries, we use zeroth-order
      // extrapolation.
      EBFaceFAB& divUDir = a_openDivU[faceDir];

      Box divBox = divUDir.getRegion();
      divBox.enclosedCells(faceDir);
      divBox.grow(faceDir,1);
      divBox &= m_domain;

      Box loBox,hiBox,centerBox,entireBox;
      int hasLo,hasHi;

      eblohicenterFace(loBox,hasLo,hiBox,hasHi,centerBox,entireBox,
                     divBox,m_domain,faceDir);

      // All of the boxes computed above are shifted so as to be cell-centered,
      // with the index of the cell center being identified with the low edge.
      // We then shift a_divVel to be compatible with that convention on input,
      // and undo the shift on output.  Nutty though this may seem,
      // it still beats the crap out of rewriting it.

      loBox.shiftHalf(faceDir,1);
      centerBox.shiftHalf(faceDir,1);
      hiBox.shiftHalf(faceDir,1);

      BaseFab<Real>& regOpenDivU = divUDir.getSingleValuedFAB();
      const BaseFab<Real>& regPrimState = a_primState.getSingleValuedFAB();
      Interval velInt = velocityInterval();
      int inormVel = velInt.begin() + faceDir;

      //this does nothing
      regOpenDivU.shiftHalf(faceDir,1);

      //put in normal direction velocity comp
      FORT_DIVUONED(CHF_FRA1(regOpenDivU,0),
                    CHF_CONST_FRA1(regPrimState,inormVel),
                    CHF_CONST_INT(faceDir),
                    CHF_BOX(centerBox));

      //add in slopes
      for (int tranDir = 0; tranDir < SpaceDim; tranDir++)
        {
          if (tranDir != faceDir)
            {
              const BaseFab<Real>& regSlopeDir =  a_slopePrim[tranDir].getSingleValuedFAB();
              Interval velInt = velocityInterval();
              int itranVel = velInt.begin() + tranDir;
              FORT_DIVUTRAN(CHF_FRA1(regOpenDivU,0),
                            CHF_CONST_FRA1(regSlopeDir,itranVel),
                            CHF_CONST_INT(faceDir),
                            CHF_BOX(centerBox));
            }
        }

      //enforce bcs
      FORT_DIVUEDGE(CHF_FRA1(regOpenDivU,0),
                    CHF_CONST_INT(faceDir),
                    CHF_BOX(loBox),
                    CHF_CONST_INT(hasLo),
                    CHF_BOX(hiBox),
                    CHF_CONST_INT(hasHi));

      //this does nothing.
      regOpenDivU.shiftHalf(faceDir,-1);

    }

  //now for the irregular cells
  for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      EBFaceFAB& openDiv = a_openDivU[faceDir];
      FaceStop::WhichFaces stopCrit;
      if (m_domain.isPeriodic(faceDir))
        {
          stopCrit = FaceStop::SurroundingWithBoundary;
        }
      else
        {
          stopCrit = FaceStop::SurroundingNoBoundary;
        }
      FaceIterator faceit(a_ivsIrreg, m_ebisBox.getEBGraph(), faceDir,
                          stopCrit);

      //reset the divergence to zero
      for (faceit.reset(); faceit.ok(); ++faceit)
        {
          const FaceIndex& face = faceit();
          openDiv(face, 0) = 0.0;
        }
      //loop through divergence directions and add
      //in diffs
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          //use the normal velocity diff if the direction
          //is the same direction as the face normal.
          //use the ave of velocity slopes otherwise

          Interval velInt = velocityInterval();
          int velInd = velInt.begin() + idir;
          if (idir == faceDir) //use velocity diff
            {
              for (faceit.reset(); faceit.ok(); ++faceit)
                {
                  const FaceIndex& face = faceit();
                  Real normalDiv = 0.0;
                  for (SideIterator sit; sit.ok(); ++sit)
                    {
                      VolIndex vof = face.getVoF(sit());
                      Real velDir = a_primState(vof, velInd);
                      normalDiv += sign(sit())*velDir;
                    }
                  openDiv(face, 0) += normalDiv;
                }
            }
          else //use ave of slopes
            {
              const EBCellFAB& slopeDir = a_slopePrim[idir];
              for (faceit.reset(); faceit.ok(); ++faceit)
                {
                  const FaceIndex& face = faceit();
                  Real tranDiv = 0.0;
                  for (SideIterator sit; sit.ok(); ++sit)
                    {
                      VolIndex vof = face.getVoF(sit());
                      Real velSlop = slopeDir(vof, velInd);
                      //the 0.5 is there because we are doing an average
                      tranDiv += 0.5*velSlop;
                    }
                  openDiv(face, 0) += tranDiv;
                }
            }

        }
    }
}
/*****************************/
void EBPatchReactive::
applyArtificialViscosity(EBFluxFAB&             a_openFlux,
                         BaseIVFAB<Real>        a_coveredFluxMinu[SpaceDim],
                         BaseIVFAB<Real>        a_coveredFluxPlus[SpaceDim],
                         const Vector<VolIndex> a_coveredFaceMinu[SpaceDim],
                         const Vector<VolIndex> a_coveredFacePlus[SpaceDim],
                         const EBCellFAB&       a_consState,
                         const EBFluxFAB&       a_divVel,
                         const Box&             a_box,
                         const IntVectSet&      a_ivsIrreg)
{
  CH_TIME("EBPatchReactive::applyArtificialViscosity");
  // Get the artificial viscosity coefficient
  Real coeff = artificialViscosityCoefficient();

  Box faceBox[SpaceDim];
  for (int dir1 = 0; dir1 < SpaceDim; ++dir1)
    {
      faceBox[dir1] = a_box;
      faceBox[dir1].grow(dir1,1);
      faceBox[dir1] &= m_domain;
      faceBox[dir1].grow(dir1,-1);
      faceBox[dir1].surroundingNodes(dir1);
    }
  //copy the flux into a temp array so to
  //not screw things up at irregular cells
  //since this is an update in place
  int numCons = numConserved();
  
  // debug
  //CH_assert(a_openFlux.getRegion() == a_box);  // commented out!!!
  // end debug
  EBFluxFAB fluxSave(m_ebisBox, a_box, numCons);

  Interval interv(0, numCons-1);
  fluxSave.copy(a_box, interv, a_box, a_openFlux, interv);

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      const Box& faceBoxDir = faceBox[idir];
      BaseFab<Real>& regFluxDir = a_openFlux[idir].getSingleValuedFAB();
      const BaseFab<Real>& regConsState = a_consState.getSingleValuedFAB();
      const BaseFab<Real>& regDivVel = a_divVel[idir].getSingleValuedFAB();
      FORT_ARTVISC(CHF_FRA(regFluxDir),
                   CHF_CONST_FRA(regConsState),
                   CHF_CONST_FRA1(regDivVel,0),
                   CHF_CONST_REAL(coeff),
                   CHF_CONST_INT(idir),
                   CHF_BOX(faceBoxDir),
                   CHF_CONST_INT(numCons),
                   CHF_CONST_REAL(m_dx[idir]));
    }

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      EBFaceFAB& fluxDir = a_openFlux[idir];
      const EBFaceFAB& fluxSaveDir = fluxSave[idir];
      const EBFaceFAB& divergeUDir = a_divVel[idir];
      //save the change in the flux so we can
      //do 0th order extrapolation to the covered faces
      EBFaceFAB fluxChange(m_ebisBox, a_box, idir, numCons);
      fluxChange.setVal(0.);
      FaceStop::WhichFaces stopCrit;
      if (m_domain.isPeriodic(idir))
        {
          stopCrit = FaceStop::SurroundingWithBoundary;
        }
      else
        {
          stopCrit = FaceStop::SurroundingNoBoundary;
        }
      for (FaceIterator  faceit(a_ivsIrreg, m_ebisBox.getEBGraph(), idir,
                               stopCrit);
          faceit.ok(); ++faceit)
        {
          const FaceIndex& face = faceit();
          VolIndex voflo = face.getVoF(Side::Lo);
          VolIndex vofhi = face.getVoF(Side::Hi);
          for (int ivar = 0; ivar < numCons; ivar++)
            {
              Real dv = divergeUDir(face, 0);
              Real uhi= a_consState(vofhi, ivar);
              Real ulo= a_consState(voflo, ivar);

              Real fsave = fluxSaveDir(face, ivar);
              Real deltaF = -coeff*std::max(-dv, (Real)0.)*(uhi-ulo);

              fluxDir(face, ivar) = fsave + deltaF;
              fluxChange(face, ivar) = deltaF;

            }
        }
      //now do the covered face fluxes
      for (SideIterator sit; sit.ok(); ++sit)
        {
          BaseIVFAB<Real>*  coveredFluxPtr;
          const Vector<VolIndex>* coveredFacePtr;
          if (sit() == Side::Lo)
            {
              coveredFacePtr = &a_coveredFaceMinu[idir];
              coveredFluxPtr = &a_coveredFluxMinu[idir];
            }
          else
            {
              coveredFacePtr = &a_coveredFacePlus[idir];
              coveredFluxPtr = &a_coveredFluxPlus[idir];
            }
          BaseIVFAB<Real>&  coveredFlux = *coveredFluxPtr;
          const Vector<VolIndex>& coveredFace = *coveredFacePtr;
          for (int ivof = 0; ivof < coveredFace.size(); ivof++)
            {
              const VolIndex& vof = coveredFace[ivof];
              if (a_box.contains(vof.gridIndex()))
                {
                  Vector<FaceIndex> flipFaces =
                    m_ebisBox.getFaces(vof, idir, flip(sit()));
                  for (int ivar = 0; ivar < numCons; ivar++)
                    {
                      Real areaTot = 0.;
                      Real diffTot = 0.;
                      for (int iface = 0; iface < flipFaces.size(); iface++)
                        {
                          const FaceIndex& face = flipFaces[iface];
                          areaTot += m_ebisBox.areaFrac(face);
                          diffTot += areaTot*fluxChange(face, ivar);
                        }
                      if (areaTot > 0.0)
                        {
                          diffTot /= areaTot;
                        }
                      coveredFlux(vof, ivar) += diffTot;
                    }
                }
            }
        }
    }
}
/*****************************/
void
EBPatchReactive::
nonconservativeDivergence(EBCellFAB&             a_divF,
                          const EBFluxFAB&       a_flux,
                          const BaseIVFAB<Real>  a_coveredFluxMinu[SpaceDim],
                          const BaseIVFAB<Real>  a_coveredFluxPlus[SpaceDim],
                          const Vector<VolIndex> a_coveredFaceMinu[SpaceDim],
                          const Vector<VolIndex> a_coveredFacePlus[SpaceDim],
                          const Box&             a_box)
{
  CH_TIME("EBPatchReactive::nonconservativeDivergence");
  CH_assert(m_isDefined);
  CH_assert(m_isBoxSet);
  CH_assert(a_divF.nComp() >= numConserved());
  CH_assert(a_flux.nComp() >= numConserved());
  CH_assert(a_flux.getRegion().contains(a_box));
  CH_assert(a_divF.getRegion().contains(a_box));

  //set the divergence initially to zero
  //then loop through directions and increment the divergence
  //with each directions flux difference.
  int ncons = numConserved();
  a_divF.setVal(0.0);
  BaseFab<Real>&       regDivF = a_divF.getSingleValuedFAB();
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      //update for the regular vofs in the nonconservative
      //case  works for all single valued vofs.
      const EBFaceFAB& fluxDir = a_flux[idir];

      const BaseFab<Real>& regFluxDir = fluxDir.getSingleValuedFAB();

      /* do the regular vofs */

      FORT_DIVERGEF( CHF_BOX(a_box),
                     CHF_FRA(regDivF),
                     CHF_CONST_FRA(regFluxDir),
                     CHF_CONST_INT(idir),
                     CHF_CONST_INT(ncons),
                     CHF_CONST_REAL(m_dx[idir]));

    }
  //update the irregular vofsn
    for (int ivof = 0; ivof<m_irregVoFs.size(); ivof++)
      {
        const VolIndex& vof = m_irregVoFs[ivof];
        if (a_box.contains(vof.gridIndex()))
          {
            //divergence was set in regular update.  we reset it
            // to zero and recalc.
            for (int ivar = 0; ivar < ncons; ivar++)
              {
                Real irregDiv = 0.0;
                for (int idir = 0; idir < SpaceDim; idir++)
                  {
                    const EBFaceFAB& fluxDir = a_flux[idir];
                    for (SideIterator sit; sit.ok(); ++sit)
                      {
                        int isign = sign(sit());
                        Vector<FaceIndex> faces =
                          m_ebisBox.getFaces(vof, idir, sit());
                        Real update = 0.;
                        for (int iface = 0; iface < faces.size(); iface++)
                          {
                            const FaceIndex& face = faces[iface];
                            Real flux = fluxDir(face, ivar);
                            update += isign*flux;

                          }
                        if (faces.size() > 1)
                          update /= Real(faces.size());
                        irregDiv += update/m_dx[idir];
                      } //end loop over sides
                  }//end loop over directions
                a_divF(vof, ivar) = irregDiv;
              }//end loop over variables
          }
    }//end loop over irreg vofs

  //now correct for the covered fluxes
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      for (SideIterator sit; sit.ok(); ++sit)
        {
          int isign = sign(sit());
          const BaseIVFAB<Real>*  coveredFluxPtr;
          const Vector<VolIndex>* coveredFacePtr;
          if (sit() == Side::Lo)
            {
              coveredFluxPtr = &a_coveredFluxMinu[idir];
              coveredFacePtr = &a_coveredFaceMinu[idir];
            }
          else
            {
              coveredFluxPtr = &a_coveredFluxPlus[idir];
              coveredFacePtr = &a_coveredFacePlus[idir];
            }
          const BaseIVFAB<Real>&  coveredFlux = *coveredFluxPtr;
          const Vector<VolIndex>& coveredFace = *coveredFacePtr;
          for (int ivof = 0; ivof < coveredFace.size(); ivof++)
            {
              const VolIndex& vof = coveredFace[ivof];
              //the sets can be bigger than the box for performance reasons.
              if (a_box.contains(vof.gridIndex()))
                {
                  //face on this side is covered.  use covered flux.
                  for (int ivar = 0; ivar < ncons; ivar++)
                    {
                      Real flux = coveredFlux(vof, ivar);
                      Real update = isign*flux/m_dx[idir];
                      a_divF(vof, ivar) += update;
                    }
                }
            }
        }
    }
}
/*****************************/
void
EBPatchReactive::
computeEBIrregFlux(BaseIVFAB<Real>&  a_ebIrregFlux,
                   const EBCellFAB&  a_primState,
                   const EBCellFAB   a_slopePrim[SpaceDim],
                   const IntVectSet& a_irregIVS,
                   const EBCellFAB&  a_source)
{
  CH_TIME("EBPatchReactive::computeEBIrregFlux");
  EBCellFAB* primPtr = (EBCellFAB*)(&a_primState);
  EBCellFAB primT;
  if (a_source.isDefined() && (s_conservativeSource))
    {
      CH_assert(a_source.box().contains(m_validBox));
      primT.define(m_ebisBox, a_source.box(), numPrimitives());
      primT.copy(a_primState);
      EBCellFAB cons(m_ebisBox, a_source.box(), numConserved());
      primToCons(cons, primT, m_validBox);
      cons.plus(a_source, 0.5*m_dt);
      consToPrim(primT, cons, m_validBox, 0);

      primPtr = &primT;
    }
  EBCellFAB& prim = *primPtr;
  Real small, smallp, smallu, smallr;
  FORT_GETSMALL(CHF_REAL(small),
                CHF_REAL(smallp),
                CHF_REAL(smallu),
                CHF_REAL(smallr));
  int numCons = numConserved();
  int numPrim = numPrimitives();
  CH_assert(isDefined());
  CH_assert(a_ebIrregFlux.nComp() == numCons);
  CH_assert(prim.nComp() == numPrim);
  CH_assert(a_ebIrregFlux.getIVS().contains(a_irregIVS));

  for (VoFIterator vofit(a_irregIVS, m_ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
//      RealVect centroid = m_ebisBox.bndryCentroid(vof);
      CH_assert(prim(vof, QRHO) > 0.0);
      
      Real     dense  = prim(vof, QRHO);
      Real     press  = max(prim(vof, QPRES),smallp);      
      RealVect veloc;
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          veloc[idir] = prim(vof, QVELX+idir);
        }      

// NOTE: HARDWIRED!!! Irreg MomFlux is set to pressure

      RealVect normal = m_ebisBox.normal(vof);

      Real xstate[QNUM];
      xstate[QRHO]  = dense;
      xstate[QPRES] = press;
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          xstate[QVELX + idir] = veloc[idir];
        }

      if(a_source.isDefined() && (!s_conservativeSource))
        {
          for(int ivar = 0; ivar < QNUM; ivar++)
            {
              Real sourceval = a_source(vof, ivar);
              xstate[ivar] += 0.5*m_dt*sourceval;
            }
        }

      //enforce positivity
      Real pext = xstate[QPRES];
      Real rext = xstate[QRHO];
      if((pext < 0.0) || (rext < 0.0))
        {
          for(int ivar = 0; ivar < numPrim; ivar++)
            {
              xstate[ivar] = prim(vof, ivar);
            }
        }
      //now compute the boundary condition assuming
      //solid walls at the embedded boundary
      RealVect velocity;
      Real     density, pressure;
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          velocity[idir] = xstate[QVELX + idir];
        }

      pressure = xstate[QPRES];
      density  = xstate[QRHO];
      Real normalvel = PolyGeom::dot(normal, velocity);
      //HACK putin true normal for computing normal vel
      //RealVect truenorm = getTrueNorm(vof);
      //Real normalvel1 = PolyGeom::dot(truenorm, velocity);
      //normalvel = normalvel1;
      //END HACK
      if (density <= 0.0)
        {
          pout() << "Density is non-positive at " << vofit() << endl;
          // FIXME: Doesn't work in parallel.
          abort();
        }
      CH_assert(density > 0.0);
      pressure = Max(pressure, smallp);

      a_ebIrregFlux(vof, CRHO) = 0.0;
      a_ebIrregFlux(vof, CENG) = 0.0;
      for (int ivar=0; ivar<m_nSpec;ivar++)
        {
          a_ebIrregFlux(vof, CSPEC1+ivar) = 0.0;
        }

      for (int idir = 0; idir < SpaceDim; idir++)
        {
          int imom = CMOMX + idir;
          //negative because our normal points into the fluid
          //these are not separated out for RZ since the pressure flux
          //is the only flux at the EB
          a_ebIrregFlux(vof, imom) = -pressure*normal[idir];
        }
    }
}
/********/
void
EBPatchReactive::
assembleFluxReg(EBFaceFAB&       a_fluxRegFlux,
                const EBFaceFAB& a_godunovFlux,
                const int&       a_idir,
                const Box&       a_cellBox)
{
  CH_assert(a_idir==a_fluxRegFlux.direction());
  CH_assert(a_idir==a_godunovFlux.direction());
  CH_assert(a_fluxRegFlux.getCellRegion().contains(a_cellBox));
  CH_assert(a_godunovFlux.getCellRegion().contains(a_cellBox));
  CH_assert(a_godunovFlux.nComp() == numFluxes());
  CH_assert(a_fluxRegFlux.nComp() == numConserved());
  CH_assert(SpaceDim==2);

  int rdir = 0;
  BaseFab<Real>& regFluxRegFlux       = a_fluxRegFlux.getSingleValuedFAB();
  const BaseFab<Real>& regGodunovFlux = a_godunovFlux.getSingleValuedFAB();
  const Box& faceRegion = regGodunovFlux.box();
  CH_assert(regFluxRegFlux.box().contains(faceRegion));

  FORT_FLUXASSEMBLE(CHF_FRA(regFluxRegFlux),
                    CHF_CONST_FRA(regGodunovFlux),
                    CHF_CONST_REAL(m_dx[0]),
                    CHF_CONST_INT(a_idir),
                    CHF_BOX(faceRegion));

  IntVectSet ivsirreg = m_ebisBox.getIrregIVS(a_cellBox);
  FaceStop::WhichFaces facestop = FaceStop::SurroundingWithBoundary;
  for (FaceIterator faceit(ivsirreg, m_ebisBox.getEBGraph(), a_idir, facestop);
      faceit.ok(); ++faceit)
    {
      const FaceIndex& face = faceit();
      Real rind = face.gridIndex(Side::Hi)[rdir];
      Real faceRad;
      if (a_idir == rdir)
        {
          faceRad = m_dx[0]*rind;
        }
      else
        {
          Real rcent = m_ebisBox.centroid(face)[rdir];
          faceRad = m_dx[0]*(rind + 0.5 + rcent);
        }

      for (int ivar=0; ivar < numConserved(); ivar++)
        {
          Real flux = faceRad*a_godunovFlux(face, ivar);
          //add in pressure term
          if (((ivar == CMOMX) && (a_idir==0)) ||
             ((ivar == CMOMY) && (a_idir==1)))
            {
              flux += faceRad*a_godunovFlux(face, CPRES);
            }
          a_fluxRegFlux(face, ivar) = flux;
        }
    }
}
/*****************************/
void
EBPatchReactive::
assembleFluxIrr(BaseIFFAB<Real>&       a_fluxRegFlux,
                const BaseIFFAB<Real>& a_godunovFlux,
                const int&             a_idir,
                const Box&             a_cellBox,
                const IntVectSet&      a_set)
{
  int rdir = 0;
  IntVectSet ivsirreg = m_ebisBox.getIrregIVS(a_cellBox);
  FaceStop::WhichFaces facestop = FaceStop::SurroundingWithBoundary;
  for (FaceIterator faceit(ivsirreg, m_ebisBox.getEBGraph(), a_idir, facestop);
      faceit.ok(); ++faceit)
    {
      const FaceIndex& face = faceit();
      Real rind = face.gridIndex(Side::Hi)[rdir];
      Real faceRad;
      if (a_idir == rdir)
        {
         faceRad = m_dx[0]*rind;
        }
      else
        {
          Real rcent = m_ebisBox.centroid(face)[rdir];
          faceRad = m_dx[0]*(rind + 0.5 + rcent);
        }
      for (int ivar=0; ivar < numConserved(); ivar++)
        {
          Real flux = faceRad*a_godunovFlux(face, ivar);
          //add in pressure term
          if (((ivar == CMOMX) && (a_idir==0)) ||
             ((ivar == CMOMY) && (a_idir==1)))
            {
              flux += faceRad*a_godunovFlux(face, CPRES);
            }
          a_fluxRegFlux(face, ivar) = flux;
        }
    }
}
/*****************************/
void
EBPatchReactive::
interpolateFluxToCentroids(BaseIFFAB<Real>              a_centroidFlux[SpaceDim],
                           const BaseIFFAB<Real>* const a_fluxInterpolant[SpaceDim],
                           const IntVectSet&            a_irregIVS)
{
  CH_TIME("EBPatchReactive::interpolateFluxToCentroids");
  FaceStop::WhichFaces stopCrit = FaceStop::SurroundingWithBoundary;

  //could just use numFluxes but this allows unit testing
  //with a single variable interpolant
  int nflux = a_centroidFlux[0].nComp();
  //now loop through the irregular faces
  for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      BaseIFFAB<Real>& fluxDir = a_centroidFlux[faceDir];
      const BaseIFFAB<Real>& interpol = *(a_fluxInterpolant[faceDir]);
      const BaseIFFAB<FaceStencil>& stencils = m_interpStencils[faceDir];

      for (FaceIterator faceit(a_irregIVS, m_ebisBox.getEBGraph(), faceDir, stopCrit);
          faceit.ok(); ++faceit)
        {
          const FaceIndex& centFace = faceit();
          const FaceStencil& sten = stencils(centFace, 0);
          for (int ivar = 0; ivar < nflux; ivar++)
            {
              Real newFlux = 0.0;
              for (int isten = 0; isten < sten.size(); isten++)
                {
                  const FaceIndex& stenFace= sten.face(isten);
                  Real weight = sten.weight(isten);
                  Real interpFlux = interpol(stenFace, ivar);
                  newFlux += weight*interpFlux;
                }
              fluxDir(centFace, ivar) = newFlux;
            }

          //debug!  turn off interpolation
          // for (int ivar = 0; ivar < numFluxes(); ivar++)
          //   {
          //     fluxDir(centFace, ivar) = interpol(centFace, ivar);
          //   }
        }
    }
}
void 
EBPatchReactive::
hybridDivergence(EBCellFAB&             a_hybridDiv,
                 EBCellFAB&             a_consState,
                 BaseIVFAB<Real>&       a_massDiff,
                 const BaseIFFAB<Real>  a_centroidFlux[SpaceDim],
                 const BaseIVFAB<Real>& a_ebIrregFlux,
                 const BaseIVFAB<Real>& a_nonConsDivF,
                 const Box&             a_box,
                 const IntVectSet&      a_ivs)
{
  CH_TIME("EBPatchReactive::irregularUpdate");
  CH_assert(isDefined());

  //compute conservative divergences of fluxes
  BaseIVFAB<Real>         conservatDivF(a_ivs, m_ebisBox.getEBGraph(), numConserved());

  consUndividedDivergence(conservatDivF, a_centroidFlux, a_ebIrregFlux, a_ivs);


  //initialize mass difference to zero
  a_massDiff.setVal(0.0);
  int ncons = numConserved();
  //update the irreg vofs
  for (VoFIterator vofit(a_ivs, m_ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      //      const IntVect& iv = vof.gridIndex();
      Real volFrac = m_ebisBox.volFrac(vof);

      for (int ivar = 0; ivar < ncons ; ivar++)
        {
          //volume fraction not divided out when divergence
          //calculated.  so we do not multiply here.

          //conservative divergence is already
          //multiplied by volfrac
          Real kapConsDiv = conservatDivF(vof, ivar);
          Real ncDivergeF = a_nonConsDivF(vof, ivar);;

          //dm =-(1-kappa)*kappa*dt*divfcons- kappa*dt*(1-kappa)*divfnc
          //dm =-(1-kappa)*dt*(kappa*divfcons- kappa*divfnc)
          Real deltaM =
            -(1.0-volFrac)*m_dt*(kapConsDiv - volFrac*ncDivergeF);

          a_hybridDiv(vof, ivar) = kapConsDiv + (1.0-volFrac)*ncDivergeF;
          a_massDiff(vof, ivar) = deltaM;
        }
    }
}
/*****************************/
void
EBPatchReactive::
consUndividedDivergence(BaseIVFAB<Real>&       a_divF,
                        const BaseIFFAB<Real>  a_centroidFlux[SpaceDim],
                        const BaseIVFAB<Real>& a_ebIrregFlux,
                        const IntVectSet&      a_ivs)
{
  CH_TIME("EBPatchReactive::consUndividedDivergence");
  CH_assert(m_isDefined);
  CH_assert(m_isBoxSet);
  int ncons = a_divF.nComp();
  a_divF.setVal(0.0);
  for (VoFIterator vofit(a_ivs, m_ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();

      Real bndryArea = m_ebisBox.bndryArea(vof);
      for (int ivar = 0; ivar < ncons ; ivar++)
        {
          Real bndryFlux = a_ebIrregFlux(vof,ivar);
          Real update = 0.;
          for ( int idir = 0; idir < SpaceDim; idir++)
            {
              const BaseIFFAB<Real>& fluxDir = a_centroidFlux[idir];
              for (SideIterator sit; sit.ok(); ++sit)
                {
                  int isign = sign(sit());
                  Vector<FaceIndex> faces =
                    m_ebisBox.getFaces(vof, idir, sit());
                  for (int iface = 0; iface < faces.size(); iface++)
                    {
                      const FaceIndex& face = faces[iface];
                      Real areaFrac = m_ebisBox.areaFrac(face);
                      Real faceFlux = fluxDir(face, ivar);
                      update += isign*areaFrac*faceFlux/m_dx[idir];
                    }
                }
            }

          //add EB boundary conditions in divergence
          update += bndryFlux*bndryArea*m_dxScale;
          //note NOT divided by volfrac
          a_divF(vof, ivar) = update;
          a_divF(vof, ivar) = update;
        } //end loop over variables
    } //end loop over vofs
}
/*****************************/
void
EBPatchReactive::
finalUpdate(EBCellFAB&              a_consState,
            BaseIVFAB<Real>&        a_massDiff,
            const BaseIVFAB<Real>&  a_nonConsDivF,
            const BaseIVFAB<Real>&  a_conservDivF,
            const IntVectSet&       a_ivs)
{
  CH_TIME("EBPatchReactive::finalUpdate");
  CH_assert(a_nonConsDivF.nComp() == numConserved());
  CH_assert(a_conservDivF.nComp() == numConserved());
  CH_assert(a_consState.nComp()   == numConserved());
  CH_assert(a_massDiff.nComp()    == numConserved());

  //initialize mass difference to zero
  a_massDiff.setVal(0.0);
  int ncons = numConserved();
  //update the irreg vofs
  for (VoFIterator vofit(a_ivs, m_ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      //      const IntVect& iv = vof.gridIndex();
      Real volFrac = m_ebisBox.volFrac(vof);

      for (int ivar = 0; ivar < ncons ; ivar++)
        {
          //volume fraction not divided out when divergence
          //calculated.  so we do not multiply here.

          //conservative divergence is already
          //multiplied by volfrac
          Real kapConsDiv = a_conservDivF(vof, ivar);
          Real ncDivergeF = a_nonConsDivF(vof, ivar);

          //dt*DF^NC has already been added in regular update
          //subtract back off to get back to the old state
          Real updateIncrementNC =  -m_dt*ncDivergeF;
          Real oldState = a_consState(vof, ivar) - updateIncrementNC;

          //un+1 = un -dt(kappa*divfcons + (1-kapp)divfnc)
          Real newState = oldState
            - m_dt*kapConsDiv
            - m_dt*(1.0-volFrac)*ncDivergeF;
          //dm =-(1-kappa)*kappa*dt*divfcons- kappa*dt*(1-kappa)*divfnc
          //dm =-(1-kappa)*dt*(kappa*divfcons- kappa*divfnc)
          Real deltaM =
            -(1.0-volFrac)*m_dt*(kapConsDiv - volFrac*ncDivergeF);

          a_consState(vof, ivar) = newState;
          a_massDiff(vof, ivar) = deltaM;
        }
    }
}
/******/
Real
EBPatchReactive::
bilinearFunc(const Real  a_WVal[2][2],
             const Real& a_xd1,
             const Real& a_xd2)
{
  Real D =                a_WVal[0][0];
  Real A =(a_WVal[1][0] - D);
  Real B =(a_WVal[0][1] - D);
  Real C =(a_WVal[1][1] - D - A - B);
  Real retval = A*a_xd1 + B*a_xd2 + C*a_xd1*a_xd2 + D;
  return retval;
}
/******/
Real
EBPatchReactive::
getMaxWaveSpeed()
{
  return s_maxWaveSpeed;
}
/******/
IntVect
EBPatchReactive::
getMaxWaveSpeedIV()
{
  return s_maxWaveSpeedIV;
}
/******/
void 
EBPatchReactive::
setMaxWaveSpeed(Real a_maxWaveSpeed)
{
  s_maxWaveSpeed = a_maxWaveSpeed;
}
/******/
void
EBPatchReactive::
setMaxWaveSpeedIV(const IntVect& a_maxWaveSpeedIV)
{
  s_maxWaveSpeedIV = a_maxWaveSpeedIV;
}
/******/
void 
EBPatchReactive::
integrateReactiveSource(EBCellFAB&  a_consState,
                        const Box&  a_box,
                        const Real& a_dt)
{
  CH_assert(a_consState.box().contains(a_box));
 
  //NOTE :: a_conState.box() and a_box are different
 
  int logflag = 0;
  EBCellFAB primState(m_ebisBox, a_box, numPrimitives());
  consToPrim(primState, a_consState, a_box, logflag);

  BaseFab<Real>& regPrim = primState.getSingleValuedFAB();

  FORT_REACTIVESRC(CHF_BOX(a_box),
                   CHF_CONST_REAL(a_dt),
                   CHF_FRA(regPrim));

 IntVectSet ivsMulti = m_ebisBox.getMultiCells(a_box);
 for (VoFIterator vofit(ivsMulti, m_ebisBox.getEBGraph()); 
     vofit.ok(); ++vofit)
   {
     const VolIndex& vof = vofit();
     Vector<Real> primitive(numPrimitives());

     for (int ivar = 0; ivar < numPrimitives(); ivar++)
      {
        primitive[ivar] = primState(vof,ivar);
      }

     FORT_POINTREACTIVESRC(CHF_CONST_REAL(a_dt),
                           CHF_VR(primitive));

     for (int ivar = 0; ivar < numPrimitives(); ivar++)
      {
        primState(vof,ivar) = primitive[ivar];
      }
   }

  primToCons(a_consState,primState,a_box);
}
/******/
/******/
void EBPatchReactive::
getRhoCv(EBCellFAB&       a_rhoCv,
         const EBCellFAB& a_rho,
         const EBCellFAB& a_massFrac,
         const EBCellFAB& a_temperat,
         const Box&       a_box)
{
  CH_assert(m_isDefined && m_isBoxSet);
  CH_assert(a_rhoCv.getRegion().contains(a_box));
  CH_assert(a_rho.getRegion().contains(a_box));
  //have to set the covered cell vals so need the cast
  //does not change real data
  EBCellFAB& rho = (EBCellFAB&) a_rho;
  EBCellFAB& temperature = (EBCellFAB&) a_temperat;
  EBCellFAB& massFrac = (EBCellFAB&) a_massFrac;
   
  rho.setInvalidData(0e0,0);
  temperature.setInvalidData(0e0,0);

  for(int i=0; i<m_nSpec; i++)
   {  
     massFrac.setInvalidData(0e0,0);
   } 
 
  const BaseFab<Real>&   regRho = a_rho.getSingleValuedFAB();
  const BaseFab<Real>&   regTemp = a_temperat.getSingleValuedFAB();
  const BaseFab<Real>&   regMassFrac = a_massFrac.getSingleValuedFAB();
  BaseFab<Real>&         regRhoCv = a_rhoCv.getSingleValuedFAB();
  
  FORT_GETRHOCV(CHF_BOX(a_box),
                CHF_CONST_FRA(regRho),
                CHF_CONST_FRA(regTemp),
                CHF_CONST_FRA(regMassFrac),
                CHF_FRA(regRhoCv));
 
  IntVectSet ivsMulti = m_ebisBox.getMultiCells(a_box);
  for (VoFIterator vofit(ivsMulti, m_ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    { 
      const VolIndex& vof = vofit();
  
      const Real rho = a_rho(vof,0);
      const Real temp = a_temperat(vof,0);
      Vector<Real> massFrac(m_nSpec);
      for (int ivar=0; ivar<m_nSpec;ivar++)
       { 
         massFrac[ivar] = a_massFrac(vof,ivar);
       } 
  
      Real rhoCv;
  
      FORT_POINTGETRHOCV(CHF_CONST_REAL(rho),
                         CHF_CONST_REAL(temp),
                         CHF_CONST_VR(massFrac),
                         CHF_REAL(rhoCv));
    } 
  a_rhoCv.setInvalidData(0e0,0);
} 
/******/
void
EBPatchReactive::
fillDiffCoeffMatrix(EBCellFAB&       a_diffCoefMatrix,
                    const EBCellFAB& a_pressure,
                    const EBCellFAB& a_temperature,
                    const EBCellFAB& a_massFrac,
                    const Box&       a_box)
{
  CH_assert(m_isDefined && m_isBoxSet);
  CH_assert(a_diffCoefMatrix.getRegion().contains(a_box));
  CH_assert(a_pressure.getRegion().contains(a_box));
  //have to set the covered cell vals so need the cast
  //does not change real data
  EBCellFAB& pressure = (EBCellFAB&) a_pressure;
  EBCellFAB& temperature = (EBCellFAB&) a_temperature;
  EBCellFAB& massFrac = (EBCellFAB&) a_massFrac;

  pressure.setInvalidData(0e0,0);
  temperature.setInvalidData(0e0,0);

  for(int i=0; i<m_nSpec; i++)
   {
     massFrac.setInvalidData(0e0,0);
   }

  const BaseFab<Real>&   regPres = a_pressure.getSingleValuedFAB();
  const BaseFab<Real>&   regTemp = a_temperature.getSingleValuedFAB();
  const BaseFab<Real>&   regMassFrac = a_massFrac.getSingleValuedFAB();
  BaseFab<Real>&         regDkMatrix = a_diffCoefMatrix.getSingleValuedFAB();

  FORT_FILLDKMATRIX(CHF_BOX(a_box),
                   CHF_CONST_FRA(regPres),
                   CHF_CONST_FRA(regTemp),
                   CHF_CONST_FRA(regMassFrac),
                   CHF_FRA(regDkMatrix));


  IntVectSet ivsMulti = m_ebisBox.getMultiCells(a_box);
  for (VoFIterator vofit(ivsMulti, m_ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();

      const Real pres = a_pressure(vof,0);
      const Real temp = a_temperature(vof,0);
      Vector<Real> massFrac(m_nSpec);
      for (int ivar=0; ivar<m_nSpec;ivar++)
       {
         massFrac[ivar] = a_massFrac(vof,ivar);
       }

      Vector<Real> DkMatrix(m_nSpec*m_nSpec);

      FORT_POINTFILLDKMATRIX(CHF_CONST_REAL(pres),
                            CHF_CONST_REAL(temp),
                            CHF_CONST_VR(massFrac),
                            CHF_VR(DkMatrix));

      for (int ivar = 0; ivar < m_nSpec*m_nSpec; ivar++)
        {
          a_diffCoefMatrix(vof, ivar)  = DkMatrix[ivar];
        }
    }

  for (int ivar = 0; ivar < m_nSpec*m_nSpec; ivar++)
   {
     a_diffCoefMatrix.setInvalidData(0e0,ivar);
   }
}
/******/
void
EBPatchReactive::
fillViscousAndConductiveCoeff(EBCellFAB&       a_aco,
                              EBCellFAB&       a_mu,
                              EBCellFAB&       a_lambda,
                              EBCellFAB&       a_kappa,
                              const EBCellFAB& a_massFrac,
                              const EBCellFAB& a_dense,
                              const EBCellFAB& a_temperature,
                              const Box&       a_box)
{
  CH_assert(m_isDefined && m_isBoxSet);
  CH_assert(a_mu.getRegion().contains(a_box));
  CH_assert(a_dense.getRegion().contains(a_box));
  //have to set the covered cell vals so need the cast
  //does not change real data
  EBCellFAB& massFrac = (EBCellFAB&) a_massFrac;
  EBCellFAB& dense = (EBCellFAB&) a_dense;
  EBCellFAB& temperature = (EBCellFAB&) a_temperature;
  dense.setInvalidData(0e0,0);
  temperature.setInvalidData(0e0,0);
  for(int i=0; i<m_nSpec; i++)
   {
     massFrac.setInvalidData(0e0,0);
   }

  const BaseFab<Real>&   regMassFrac = a_massFrac.getSingleValuedFAB();
  const BaseFab<Real>&   regDense = a_dense.getSingleValuedFAB();
  const BaseFab<Real>&   regTemp = a_temperature.getSingleValuedFAB();

  BaseFab<Real>&         regAco = a_aco.getSingleValuedFAB();
  BaseFab<Real>&         regMu = a_mu.getSingleValuedFAB();
  BaseFab<Real>&         regLambda = a_lambda.getSingleValuedFAB();
  BaseFab<Real>&         regKappa = a_kappa.getSingleValuedFAB();

  FORT_FILLVISCANDCONDUCTIVECOEFF(CHF_BOX(a_box),
                                  CHF_CONST_FRA(regMassFrac),
                                  CHF_CONST_FRA(regDense),
                                  CHF_CONST_FRA(regTemp),
                                  CHF_FRA(regAco),
                                  CHF_FRA(regMu),
                                  CHF_FRA(regLambda),
                                  CHF_FRA(regKappa));

  IntVectSet ivsMulti = m_ebisBox.getMultiCells(a_box);
  for (VoFIterator vofit(ivsMulti, m_ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();

      Vector<Real> massFrac(m_nSpec);

      for (int ivar = 0; ivar < m_nSpec; ivar++)
        {
          massFrac[ivar] = a_massFrac(vof, ivar);
        }

      const Real dense = a_dense(vof,0);
      const Real temp = a_temperature(vof,0);

      Real aco;
      Real mu;
      Real lambda;
      Real kappa;

      FORT_POINTFILLVISCANDCONDUCTIVECOEFF(CHF_CONST_VR(massFrac),
                                           CHF_CONST_REAL(dense),
                                           CHF_CONST_REAL(temp),
                                           CHF_REAL(aco),
                                           CHF_REAL(mu),
                                           CHF_REAL(lambda),
                                           CHF_REAL(kappa));

      a_aco(vof,0)  = aco;
      a_mu(vof, 0)  = mu;
      a_lambda(vof,0) = lambda;
      a_kappa(vof,0) = kappa;
    }
  a_aco.setInvalidData(0e0,0);
  a_mu.setInvalidData(0e0,0);
  a_lambda.setInvalidData(0e0,0);
  a_kappa.setInvalidData(0e0,0);
}
/******/
#include "NamespaceFooter.H"
