#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MIXEDPOISSONEBBC_H_
#define _MIXEDPOISSONEBBC_H_

#include "RefCountedPtr.H"
#include "Stencils.H"
#include "BaseEBBC.H"
#include "BaseMixBCValue.H"
#include "LayoutData.H"
#include "IntVectSet.H"
#include "EBStencil.H"
#include "NamespaceHeader.H"

///
/**
 */
class MixedPoissonEBBC: public BaseEBBC
{
public:

  ///
  MixedPoissonEBBC();

  MixedPoissonEBBC(const ProblemDomain& a_domain,
                       const EBISLayout&    a_layout,
                       const RealVect&      a_dx,
                       const IntVect*       a_ghostCellsPhi=0,
                       const IntVect*       a_ghostCellsRhs=0);

  void construct(const ProblemDomain& a_domain,
                 const EBISLayout&    a_layout,
                 const RealVect&      a_dx,
                 const IntVect*       a_ghostCellsPhi=0,
                 const IntVect*       a_ghostCellsRhs=0);


 ///
  /**
   */
  virtual ~MixedPoissonEBBC();
  ///
  /**
   */
  void setArguments(Real 							a_dvalue, 
					Real 							a_nvalue, 
					RefCountedPtr<BaseMixBCValue> 	a_func,
					int								a_ebType,
					RealVect					    a_dymin,
					Real 							a_epsLo,
					Real							a_epsHi);  
  ///
  /**
  */
  virtual void setOrder(int a_order);
  ///
  /**
   */
  virtual void define(const LayoutData<IntVectSet>& a_cfivs,
                      const Real&                   a_factor);
  ///
  /**
  //MD added for dielectric/gas interface at highest level of refinement
   */
  void defineCustom(const LayoutData<IntVectSet>& a_cfivs,
                    const Real&                   a_factor);
  ///
  /**
   */
  virtual LayoutData<BaseIVFAB<VoFStencil> >* getFluxStencil(int ivar);


  LayoutData<BaseIVFAB<Real> >& getFluxWeight()
  {
    return m_fluxWeight;
  }
  ///


   /**
   */
  virtual void applyEBFlux(EBCellFAB&                    a_lphi,
                           const EBCellFAB&              a_phi,
                           VoFIterator&                  a_vofit,
                           const LayoutData<IntVectSet>& a_cfivs,
                           const DataIndex&              a_dit,
                           const RealVect&               a_probLo,
                           const RealVect&               a_dx,
                           const Real&                   a_factor,
                           const bool&                   a_useHomogeneous,
                           const Real&                   a_time);
  ///

  /**  deprecated interface
   */
  /*
  virtual void getEBFlux(Real&                         a_flux,
                         const VolIndex&               a_vof,
                         const LevelData<EBCellFAB>&   a_phi,
                         const LayoutData<IntVectSet>& a_cfivs,
                         const DataIndex&              a_dit,
                         const RealVect&               a_probLo,
                         const RealVect&               a_dx,
                         const bool&                   a_useHomogeneous,
                         const Real&                   a_time,
                         const pair<int,Real>*         a_cacheHint=0 );

  */

  virtual void setEBBCSource(const LevelData<BaseIVFAB<Real> >&     a_rhos);

  bool getSecondOrderStencil(VoFStencil&          a_stencil,
                             Real&                a_weight,
                             Vector<VoFStencil>&  a_pointStencils,
                             Vector<Real>&        a_distanceAlongLine,
                             const VolIndex&      a_vof,
                             const EBISBox&       a_ebisBox,
                             const RealVect&      a_dx,
                             const IntVectSet&    a_cfivs);

  //LM added routine
  void setTime(Real* a_time)
  {
    if(m_timePtr == NULL) m_isDefined = false;
    m_timePtr = a_time;
  }

  //hack to let component of velocity get set
  static int s_velComp;
  static bool s_areaFracWeighted;

  static Real s_timePrint;

protected:

private:
  virtual void getFirstOrderStencil(VoFStencil&     a_stencil,
                                    Real&           a_weight,
                                    const VolIndex& a_vof,
                                    const EBISBox&  a_ebisBox,
                                    const RealVect& a_dx);

  virtual void getSecondOrderStencil(VoFStencil&     a_stencil,
                                     Real&           a_weight,
                                     const VolIndex& a_vof,
                                     const EBISBox&  a_ebisBox,
                                     const RealVect& a_dx,
                                     const IntVectSet& a_cfivs);


public:
  bool m_isDefined;
  Real m_dvalue;
  Real m_nvalue;
  RefCountedPtr<BaseMixBCValue> m_func;
  int m_ebType;
  bool m_onlyHomogeneous;
  int m_order;
 
  ProblemDomain m_domain;
  EBISLayout    m_layout;

  RealVect m_dx;
  RealVect m_dymin;  
  bool m_finest;
  Real m_epsLo;
  Real m_epsHi;
  Real* m_timePtr;

  //stencils for operator evaluation
  IntVect                   m_ghostCellsPhi;
  IntVect                   m_ghostCellsRHS;

  //the stuff below is to keep EBAMRPoissonOp working in the transition
/*
  virtual void define(const LayoutData<IntVectSet>& a_cfivs);
*/
  // m_fluxWeight is the weight of the Dirichlet BC value used when
  // calculating the EB flux as a linear combination of the stencil
  // cell values and the BC value.
  LayoutData<BaseIVFAB<Real> >       m_fluxWeight;
  LayoutData<BaseIVFAB<Real> >       m_custWeight;
  LayoutData<BaseIVFAB<VoFStencil> > m_fluxStencil;
  LayoutData<BaseIVFAB<VoFStencil> > m_custStencil;
  LayoutData<BaseIVFAB<Real> >       m_rhos;


};

///
/**
 */
class MixedPoissonEBBCFactory: public BaseEBBCFactory
{
public:
  ///
  /**
   */
  MixedPoissonEBBCFactory();

  ///
  /**
   */
  virtual ~MixedPoissonEBBCFactory();

  ///
  /**
   */
  virtual void setOrder(int a_order);

  ///
  /**
   */
  void setArguments(Real a_dvalue, Real a_nvalue, RefCountedPtr<BaseMixBCValue> a_func, int a_ebType, RealVect a_dymin, Real a_epsLo, Real a_epsHi);

  // LM added
  void setTime(Real* a_time);

  ///
  /**
   */
  virtual MixedPoissonEBBC* create(const ProblemDomain& a_domain,
                                       const EBISLayout&    a_layout,
                                       const RealVect&      a_dx,
                                       const IntVect*       a_ghostCellsPhi=0,
                                       const IntVect*       a_ghostCellsRhs=0);

private:
  Real m_dvalue;
  Real m_nvalue;
  RefCountedPtr<BaseMixBCValue> m_func;
  int m_ebType;
  bool m_onlyHomogeneous;
  int m_order;
  RealVect m_dymin;
  Real m_epsLo;
  Real m_epsHi;
  Real* m_timePtr;

};

#include "NamespaceFooter.H"
#endif
