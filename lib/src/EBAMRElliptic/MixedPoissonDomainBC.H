#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MIXEDPOISSONDOMAINBC_H_
#define _MIXEDPOISSONDOMAINBC_H_

#include "RefCountedPtr.H"
#include "BaseDomainBC.H"
#include "BaseMixBCValue.H"

#include "NamespaceHeader.H"

///
/**
 */
class MixedPoissonDomainBC: public BaseDomainBC
{
public:
  ///
  /**
   */
  MixedPoissonDomainBC();

  ///
  /**
   */
  virtual ~MixedPoissonDomainBC();

  ///
  /**
   */
  void setArguments(Vector<Vector<int> >&  a_domMixBc, Real a_dvalue, Real a_nvalue, RefCountedPtr<BaseMixBCValue> a_func);  
  
  ///
   /** use this for order of domain boundary
   */
  virtual void setEBOrder(int a_ebOrder);

  ///
  /**
     This is used by the projections to get
     velocity at domain faces where eb intersects.  Sets velocity to a dirichlet value.
   */
  virtual void getFaceVel(Real&                 a_faceFlux,
                          const FaceIndex&      a_face,
                          const EBFluxFAB&      a_vel,
                          const RealVect&       a_probLo,
                          const RealVect&       a_dx,
                          const int&            a_idir,
                          const int&            a_icomp,
                          const Real&           a_time,
                          const Side::LoHiSide& a_side,
                          const bool&           a_doDivFreeOutflow);

  ///
  /**
   */
  virtual void getFaceFlux(BaseFab<Real>&        a_faceFlux,
                           const BaseFab<Real>&  a_phi,
                           const RealVect&       a_probLo,
                           const RealVect&       a_dx,
                           const int&            a_idir,
                           const Side::LoHiSide& a_side,
                           const DataIndex&      a_dit,
                           const Real&           a_time,
                           const bool&           a_useHomogeneous);

  /**
     Elliptic solver flux.
   */
  virtual void getFaceFlux(Real&                 a_faceFlux,
                           const VolIndex&       a_vof,
                           const int&            a_comp,
                           const EBCellFAB&      a_phi,
                           const RealVect&       a_probLo,
                           const RealVect&       a_dx,
                           const int&            a_idir,
                           const Side::LoHiSide& a_side,
                           const DataIndex&      a_dit,
                           const Real&           a_time,
                           const bool&           a_useHomogeneous);

  ///
  /**
     Elliptic solver flux.
   */
  void getInhomFaceFlux(Real&                 a_faceFlux,
                        const VolIndex&       a_vof,
                        const int&            a_comp,
                        const EBCellFAB&      a_phi,
                        const RealVect&       a_probLo,
                        const RealVect&       a_dx,
                        const int&            a_idir,
                        const Side::LoHiSide& a_side,
                        const DataIndex&      a_dit,
                        const Real&           a_time);


///
/**
   A query of whether a_jvof is Dirichlet Domain boundary to a_ivof. a_jvof is the ghost vof; a_ivof is the valid computational vof.
   If the code runs into here it's got to be Dirichlet BC
   [MD] check if this true with mixed BC
 */
  bool isDirichletDom(const VolIndex&   a_ivof,
                    const VolIndex&   a_jvof,
                    const EBCellFAB&  a_phi) const
  {
   	return true;
  }

  ///
  /**
     This is used by the projections to get grad(phi) at domain faces (usually outflow).
     Called by InflowOutflowPoissonDomainBC::getFaceFlux.
   */
  virtual void getFaceGradPhi(Real&                 a_faceFlux,
                              const FaceIndex&      a_face,
                              const int&            a_comp,
                              const EBCellFAB&      a_phi,
                              const RealVect&       a_probLo,
                              const RealVect&       a_dx,
                              const int&            a_idir,
                              const Side::LoHiSide& a_side,
                              const DataIndex&      a_dit,
                              const Real&           a_time,
                              const bool&           a_useAreaFrac,
                              const RealVect&       a_centroid,
                              const bool&           a_useHomogeneous);

  /// 
  /**
     From Neumann:
	 This is used by the projections to get extrapolated
     grad(phi) at domain faces for cell-centered gradient.
     Called by getFaceFlux.
   */
  virtual void getFaceFluxGradPhi(Real&                 a_faceFlux,
                                  const FaceIndex&      a_face,
                                  const int&            a_comp,
                                  const EBCellFAB&      a_phi,
                                  const RealVect&       a_probLo,
                                  const RealVect&       a_dx,
                                  const int&            a_idir,
                                  const Side::LoHiSide& a_side,
                                  const DataIndex&      a_dit,
                                  const Real&           a_time,
                                  const bool&           a_useAreaFrac,
                                  const RealVect&       a_centroid,
                                  const bool&           a_useHomogeneous);

  ///
  /**
     This is used by the projections to get grad(phi) at domain faces (usually outflow).
     Called by InflowOutflowPoissonDomainBC::getFaceFlux.
   */
  void getInhomFaceGradPhi(Real&                 a_faceFlux,
                           const FaceIndex&      a_face,
                           const int&            a_comp,
                           const EBCellFAB&      a_phi,
                           const RealVect&       a_probLo,
                           const RealVect&       a_dx,
                           const int&            a_idir,
                           const Side::LoHiSide& a_side,
                           const DataIndex&      a_dit,
                           const Real&           a_time,
                           const bool&           a_useAreaFrac,
                           const RealVect&       a_centroid);

  ///

  /**
	From Neumann:
   */
  virtual void getInhomFaceFluxGradPhi(Real&                 a_faceFlux,
                                       const FaceIndex&      a_face,
                                       const int&            a_comp,
                                       const EBCellFAB&      a_phi,
                                       const RealVect&       a_probLo,
                                       const RealVect&       a_dx,
                                       const int&            a_idir,
                                       const Side::LoHiSide& a_side,
                                       const DataIndex&      a_dit,
                                       const Real&           a_time,
                                       const bool&           a_useAreaFrac,
                                       const RealVect&       a_centroid);
  /**
     Stencils
   */

  virtual void getFluxStencil(      VoFStencil&      a_stencil,
                              const VolIndex&        a_vof,
                              const int&             a_comp,
                              const RealVect&        a_dx,
                              const int&             a_idir,
                              const Side::LoHiSide&  a_side,
                              const EBISBox&         a_ebisBox);

  virtual void getFluxStencil(      VoFStencil&      a_stencil,
                              const FaceIndex&       a_face,
                              const int&             a_comp,
                              const RealVect&        a_dx,
                              const int&             a_idir,
                              const Side::LoHiSide&  a_side,
                              const EBISBox&         a_ebisBox);

  void getFirstOrderFluxStencil(      VoFStencil&      a_stencil,
                                const FaceIndex&       a_face,
                                const int&             a_comp,
                                const RealVect&        a_dx,
                                const int&             a_idir,
                                const Side::LoHiSide&  a_side,
                                const EBISBox&         a_ebisBox);

  void getSecondOrderFluxStencil(      VoFStencil&      a_stencil,
                                 const FaceIndex&       a_face,
                                 const int&             a_comp,
                                 const RealVect&        a_dx,
                                 const int&             a_idir,
                                 const Side::LoHiSide&  a_side,
                                 const EBISBox&         a_ebisBox);

  bool m_isDefined;

  void setTime(Real* a_time)
  {
    if(m_timePtr == NULL) m_isDefined = false;
    m_timePtr = a_time;
  }

private:
  Vector<Vector<int> >  m_domMixBc;
  Real m_dvalue;
  Real m_nvalue;
  RefCountedPtr<BaseMixBCValue> m_func;
  bool m_onlyHomogeneous;
  int m_ebOrder;
  Real* m_timePtr;
};

class MixedPoissonDomainBCFactory: public BaseDomainBCFactory
{
public:
  ///
  /**
   */
  MixedPoissonDomainBCFactory();

  ///
  /**
   */
  virtual ~MixedPoissonDomainBCFactory();

  ///
  /**
   */
  virtual MixedPoissonDomainBC* create(const ProblemDomain& a_domain,
                                           const EBISLayout&    a_layout,
                                           const RealVect&      a_dx);
  ///
  /**
   */
  void setArguments(Vector<Vector<int> >& a_domMixBc, Real a_dvalue, Real a_nvalue, RefCountedPtr<BaseMixBCValue> a_func);

  ///
  /**
   */
  virtual void setEBOrder(int a_ebOrder);

  void setTime(Real* a_time);

private:
  Vector<Vector<int> >   m_domMixBc;
  bool m_onlyHomogeneous;
  Real m_dvalue;
  Real m_nvalue;
  RefCountedPtr<BaseMixBCValue> m_func;
  int m_ebOrder;
  Real* m_timePtr;

};
#include "NamespaceFooter.H"
#endif
