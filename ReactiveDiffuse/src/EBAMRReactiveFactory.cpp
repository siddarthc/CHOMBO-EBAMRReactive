#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "AMRLevel.H"
#include "EBAMRReactive.H"
#include "EBAMRReactiveFactory.H"
#include "NamespaceHeader.H"
/****************/
/****************/
AMRLevel*
EBAMRReactiveFactory::
new_amrlevel() const
{
  EBAMRReactive* amrr_ptr = new EBAMRReactive();

  amrr_ptr->CFL(m_cfl);
  amrr_ptr->patchReactive(m_patchReactive);
  amrr_ptr->refinementThreshold(m_refineThresh);
  amrr_ptr->tagBufferSize(m_tagBufferSize);
  amrr_ptr->initialDtMultiplier(m_initialDtMultiplier);
  amrr_ptr->domainLength(m_domainLength);
  amrr_ptr->redistRadius(m_redistRad);
  amrr_ptr->verbosity(m_verbosity);
  amrr_ptr->useMassRedistribution(m_useMassRedist);
  amrr_ptr->doRZCoords(m_doRZCoords);
  amrr_ptr->hasSourceTerm(m_hasSourceTerm);
  amrr_ptr->doSmushing(m_doSmushing);
  amrr_ptr->tagAll(m_tagAll);
  amrr_ptr->addReactionRates(m_addReactionRates);
  amrr_ptr->addDiffusion(m_addDiffusion);

  return (static_cast <AMRLevel*> (amrr_ptr));
}
/****************/
/****************/
EBAMRReactiveFactory::
~EBAMRReactiveFactory()
{
}
/****************/
/****************/
EBAMRReactiveFactory::
EBAMRReactiveFactory(const Real&                   a_initialDtMultiplier,
                     const Real&                   a_cfl,
                     const int &                   a_redistRad,
                     const RealVect&               a_domainLength,
                     const Real&                   a_refineThresh,
                     const int &                   a_tagBufferSize,
                     const int &                   a_verbosity,
                     const bool&                   a_useMassRedist,
                     const bool&                   a_doSmushing,
                     const bool&                   a_doRZCoords,
                     const bool&                   a_hasSourceTerm,
                     const bool&                   a_addReactionRates,
                     const bool&                   a_addDiffusion,
                     const EBPatchReactiveFactory* const a_patchReactive,
                     bool                          a_tagAll)
{
  m_tagAll = a_tagAll;

  m_initialDtMultiplier = a_initialDtMultiplier;
  m_cfl                 = a_cfl;
  m_useMassRedist       = a_useMassRedist;
  m_redistRad           = a_redistRad;
  m_domainLength        = a_domainLength;
  m_refineThresh        = a_refineThresh;
  m_tagBufferSize       = a_tagBufferSize;
  m_verbosity           = a_verbosity;
  m_doSmushing          = a_doSmushing;
  m_doRZCoords          = a_doRZCoords;
  m_hasSourceTerm       = a_hasSourceTerm;
  m_patchReactive       = a_patchReactive;
  m_addReactionRates    = a_addReactionRates;
  m_addDiffusion        = a_addDiffusion;
}
/****************/
/****************/
#include "NamespaceFooter.H"
                                                            
