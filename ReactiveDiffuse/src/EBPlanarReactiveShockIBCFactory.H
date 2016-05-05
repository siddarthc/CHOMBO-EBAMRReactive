#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif


#ifndef _EBPLANARREACTIVESHOCKIBCFACTORY_H_
#define _EBPLANARREACTIVESHOCKIBCFACTORY_H_

#include "EBCellFAB.H"
#include "EBFaceFAB.H"
#include "REAL.H"
#include "LevelData.H"
#include "ProblemDomain.H"
#include "EBPhysIBC.H"
#include "EBPhysIBCFactory.H"

#include "EBPlanarReactiveShockIBC.H"

//#include "EBPlanarReactiveShockF_F.H"

///
/**
 */
class EBPlanarReactiveShockIBCFactory: public EBPhysIBCFactory
{
public:

  ///
  ~EBPlanarReactiveShockIBCFactory();

  ///
  EBPlanarReactiveShockIBCFactory(const Real&     a_shockvel,
                                  const Real&     a_center,
                                  const int&      a_shocknorm,
                                  const bool&     a_shockbackward,
                                  const bool&     a_doRZCoords);

  ///
  EBPhysIBC* create() const;

protected:
  Real m_shockvel;
  Real m_center;
  int  m_shocknorm;
  bool m_shockbackward;
  bool m_doRZCoords;

  bool m_isChemistrySet;

private:
  EBPlanarReactiveShockIBCFactory()
  {
    MayDay::Error("Invalid operator");
  }

  void operator=(const EBPlanarReactiveShockIBCFactory& a_input)
  {
    MayDay::Error("Invalid operator");
  }

  EBPlanarReactiveShockIBCFactory(const EBPlanarReactiveShockIBCFactory& a_input)
  {
    MayDay::Error("Invalid operator");
  }


};

#endif

