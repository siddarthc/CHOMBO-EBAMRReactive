#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EBPlanarReactiveShockIBCFactory.H"
#include "EBPlanarReactiveShockIBC.H"
#include "EBPlanarReactiveShockF_F.H"
#include "EBPatchReactiveF_F.H"

void initChem();
/******************/
/******************/
EBPlanarReactiveShockIBCFactory::
EBPlanarReactiveShockIBCFactory(const Real&     a_shockvel,
                                const Real&     a_center,
                                const int&      a_shocknorm,
                                const bool&     a_shockbackward,
                                const bool&     a_doRZCoords)
  :EBPhysIBCFactory()
{
  initChem();
  m_shockvel      = a_shockvel;
  m_center        = a_center;
  m_shocknorm     = a_shocknorm;
  m_shockbackward = a_shockbackward;
  m_doRZCoords    = a_doRZCoords;
}
/******************/
/******************/

EBPlanarReactiveShockIBCFactory::
~EBPlanarReactiveShockIBCFactory()
{
}
/******************/
/******************/
EBPhysIBC*
EBPlanarReactiveShockIBCFactory::
create() const
{
  EBPlanarReactiveShockIBC* retval =
    new EBPlanarReactiveShockIBC(m_shockvel,
                                 m_center,
                                 m_shocknorm,
                                 m_shockbackward,
                                 m_doRZCoords);
  return static_cast<EBPhysIBC*>(retval);
}
/******************/
/******************/
void initChem()
{
  FORT_INITIALIZE_CHEMISTRY();
//  m_isChemistrySet = true;
}
/******************/
/******************/

