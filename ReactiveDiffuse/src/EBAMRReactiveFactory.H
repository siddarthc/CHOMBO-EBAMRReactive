#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBAMRREACTIVEFACTORY_H_
#define _EBAMRREACTIVEFACTORY_H_

#include "EBAMRReactive.H"
#include "AMRLevelFactory.H"
#include "EBPatchReactiveFactory.H"
//#include "NamespaceHeader.H"

///
/**
 */
class EBAMRReactiveFactory : public AMRLevelFactory
{
public:
  ///
  AMRLevel* new_amrlevel() const;

  ///
  EBAMRReactiveFactory(const Real&                         a_initialDtMultiplier,
                       const Real&                         a_cfl,
                       const int &                         a_redistRad,
                       const RealVect&                     a_domainLength,
                       const Real&                         a_refineThresh,
                       const int &                         a_tagBufferSize,
                       const int &                         a_verbosity,
                       const bool&                         a_useMassRedist,
                       const bool&                         a_doSmushing,
                       const bool&                         a_doRZCoords,
                       const bool&                         a_hasSourceTerm,
                       const bool&                         a_addReactionRates,
                       const bool&                         a_addDiffusion,
                       const EBPatchReactiveFactory* const a_patchReactive,
                       bool                                a_tagAll = false);

  ///
  virtual ~EBAMRReactiveFactory();

protected:
  bool                   m_tagAll;
  Real                   m_initialDtMultiplier;
  Real                   m_cfl;
  int                    m_redistRad;
  RealVect               m_domainLength;
  Real                   m_refineThresh;
  int                    m_tagBufferSize;
  int                    m_verbosity;
  bool                   m_useMassRedist;
  bool                   m_doSmushing;
  bool                   m_doRZCoords;
  bool                   m_hasSourceTerm;
  bool                   m_addReactionRates;
  bool                   m_addDiffusion;
  const EBPatchReactiveFactory* m_patchReactive;


private:
  EBAMRReactiveFactory()
  {
    MayDay::Error("invalid operator");
  }

  //disallowed for all the usual reasons
  void operator=(const EBAMRReactiveFactory& a_input)
  {
    MayDay::Error("invalid operator");
  }
  EBAMRReactiveFactory(const EBAMRReactiveFactory& a_input)
  {
    MayDay::Error("invalid operator");
  }

};

//#include "NamespaceFooter.H"
#endif

