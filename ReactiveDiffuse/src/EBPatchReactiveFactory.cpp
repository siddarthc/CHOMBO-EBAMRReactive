#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EBPatchReactiveFactory.H"
#include "EBPatchReactiveF_F.H"
#include "NamespaceHeader.H"

//void initChem();
int getnSpecies();
//-----------------------------------------------------------------------
EBPatchReactiveFactory::
EBPatchReactiveFactory(const EBPhysIBCFactory*       a_bcFactoryPtr,
                       const bool&                   a_useFourthOrderSlopes,
                       const bool&                   a_useZeroSlopes,
                       const bool&                   a_useFlattening,
                       const bool&                   a_useArtificialVisc,
                       const bool&                   a_useLimiting,
                       const bool&                   a_doRZCoords)
{
//  initChem();
  m_useLimiting          = a_useLimiting;
  m_useFourthOrderSlopes = a_useFourthOrderSlopes;
  m_useZeroSlopes        = a_useZeroSlopes;
  m_bcFactoryPtr         = a_bcFactoryPtr;
  m_useFlattening        = a_useFlattening;
  m_useArtificialVisc    = a_useArtificialVisc;
  m_doRZCoords           = a_doRZCoords;
  m_nSpec                = getnSpecies();
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
/******************/
EBPatchReactiveFactory::
~EBPatchReactiveFactory()
{
}
/******************/
/******************/
EBPatchReactive*
EBPatchReactiveFactory::
create() const 
{
  EBPatchReactive* retval = new EBPatchReactive();
  retval->setnSpecies(m_nSpec);
  retval->setSlopeParameters(m_useFourthOrderSlopes,
                             m_useZeroSlopes,
                             m_useFlattening,
                             m_useLimiting);
  retval->artificialViscosity(m_useArtificialVisc);
  retval->setEBPhysIBC(*m_bcFactoryPtr);
  retval->setGamma(m_gamma);
  retval->setSpecHeat(m_specHeat);
  retval->doRZCoords(m_doRZCoords);
  return retval;
}
/******************/
/******************/
/*
void initChem()
{
  FORT_INITCHEM();
}
*/
/******************/
/******************/
int getnSpecies()
{
  int nSpecies;
  FORT_GET_NSPECIES(CHF_INT(nSpecies));
  return nSpecies;
}
/******************/
/******************/
#include "NamespaceFooter.H"

