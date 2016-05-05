#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBPATCHREACTIVEFACTORY_H_
#define _EBPATCHREACTIVEFACTORY_H_

#include "EBPatchReactive.H"
//#include "EBPatchChemistryFactory.H"
//#include "NamespaceHeader.H"

///
/**
 */
class EBPatchReactiveFactory
{
public:
  //! Destructor.
  virtual ~EBPatchReactiveFactory();

  ///
  /**
   */
  EBPatchReactive* create() const;

  //! Create a new factory, specifying algorithmic parameters.
  //! \param a_useFourthOrderSlopes If set, the EBPatchReactive objects created
  //!                               by this factory will use fourth-order estimates
  //!                               for slopes.
  EBPatchReactiveFactory(const EBPhysIBCFactory* const a_bc,
                         const bool&                   a_useFourthOrderSlopes,
                         const bool&                   a_useZeroSlopes,
                         const bool&                   a_useFlattening,
                         const bool&                   a_useArtificialVisc,
                         const bool&                   a_useLimiting,
                         const bool&                   a_doRZCoords);


protected:
  const EBPhysIBCFactory* m_bcFactoryPtr;
  Real                    m_gamma;
  Real                    m_specHeat;
  bool                    m_useFourthOrderSlopes;
  bool                    m_useZeroSlopes;
  bool                    m_useFlattening;
  bool                    m_useArtificialVisc;
  bool                    m_doRZCoords;
  bool                    m_useLimiting;

  int                     m_nSpec;
private:

  //disallowed for all the usual reasons
  EBPatchReactiveFactory()
  {
    MayDay::Error("invalid operator");
  }
  void operator=(const EBPatchReactiveFactory& a_input)
  {
    MayDay::Error("invalid operator");
  }
  EBPatchReactiveFactory(const EBPatchReactiveFactory& a_input)
  {
    MayDay::Error("invalid operator");
  }
};
//#include "NamespaceFooter.H"
#endif
