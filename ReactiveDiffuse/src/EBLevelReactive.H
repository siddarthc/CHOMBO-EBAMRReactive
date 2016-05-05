#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBLEVELREACTIVE_H_
#define _EBLEVELREACTIVE_H_

#include "EBCellFAB.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "PiecewiseLinearFillPatch.H"
#include "EBPWLFillPatch.H"
#include "EBFluxRegister.H"
#include "ProblemDomain.H"
#include "EBLevelRedist.H"
//#include "NamespaceHeader.H"

#include "EBPatchReactiveFactory.H"

/// Level Reactive
/**
 */
class EBLevelReactive
{
public:
  /// Default constructor
  /**
     Object requires define(..) to be called before all other functions.
  */
  EBLevelReactive();

  /// Destructor
  /**
     destroys all objects created by define(..). Passed in
     data references of define(..) are left alone.
  */
  ~EBLevelReactive();

  /// Actual constructor.
  /**
     For the coarsest level, an empty DisjointBoxLayout
     is passed in for coaserDisjointBoxLayout.
     Inside the routine, we cast away const-ness on the data members
     for the assignment.  The arguments passed in are maintained const.
     (coding standards).  a_nRefine is the refinement ratio between
     this level and the next coarser level.  a_numGhosts is the number
     of ghost cells in each direction.
  */
  void define(const DisjointBoxLayout&            a_thisDBL,
              const DisjointBoxLayout&            a_coarDBL,
              const EBISLayout&                   a_thisEBISL,
              const EBISLayout&                   a_coarEBISL,
              const ProblemDomain&                a_DProblem,
              const int&                          a_nRefine,
              const RealVect&                     a_dx,
              const bool&                         a_useMassRedist,
              const bool&                         a_doSmushing,
              const bool&                         a_doRZCoords,
              const bool&                         a_hasSourceTerm,
              const EBPatchReactiveFactory* const a_patchReactive,
              const bool&                         a_hasCoarser,
              const bool&                         a_hasFiner);


  /// Get maximum wave speed
  Real getMaxWaveSpeed(const LevelData<EBCellFAB>& a_state);

  bool isDefined() const;

  /// outputs divergeF.  a_cons gets its ghost cells filled
  /**
     For the coarsest level, UcoarseOld and UCoarseNew are empty
     LevelData<EBCellFAB> objects.  Also, either UCoarseOld or
     UCoarseNew might be empty to indicate that t(nf) = t(nc) the
     one grid is at the current time and interpolation is not required
     for boundary condition generation.  U must be defined on the same
     DisjointBoxLayouts as were used in define(...)
     Coarse flux register is flux register with the coarse level.
     Fine flux register is the flux register with the fine level.
     To the finer level FR, this level is the coarse level.
     To the coarser level FR, this level is the fine level.
     outputs divergeF.  a_cons gets its ghost cells filled
  */
  void divergeF(LevelData<EBCellFAB>&         a_divergeF,
                LevelData<BaseIVFAB<Real> >&  a_massDiff,
                EBFluxRegister&               a_fineFluxRegister,
                EBFluxRegister&               a_coarFluxRegister,
                LevelData<EBCellFAB>&         a_consState,   //not really changed
                const LevelData<EBCellFAB>&   a_source,
                const LevelData<EBCellFAB>&   a_consStateCoarseOld,
                const LevelData<EBCellFAB>&   a_consStateCoarseNew,
                const Real&                   a_time,
                const Real&                   a_coarTimeOld,
                const Real&                   a_coarTimeNew,
                const Real&                   a_dt);


  ///
  void floorConserved(LevelData<EBCellFAB>&   a_consState,
                      Real                    a_time, 
                      Real                    a_dt);

  void setHyperbolicSrc(LevelData<EBCellFAB>* a_hyperbolicSrc);

  void integrateReactiveSource(LevelData<EBCellFAB>& a_consState,
                               const Box&            a_domain,  
                               const Real&           a_time,   
                               const Real&           a_dt);
                               

protected:
  void fillConsState(LevelData<EBCellFAB>&         a_consState,
                     const LevelData<EBCellFAB>&   a_consStateCoarseOld,
                     const LevelData<EBCellFAB>&   a_consStateCoarseNew,
                     const Real&                   a_time,
                     const Real&                   a_coarTimeOld,
                     const Real&                   a_coarTimeNew);

  void computeFlattening(Real a_time, Real a_dt,
                         LevelData<EBCellFAB>&         a_consState);

  void doRegularUpdate(LevelData<EBCellFAB>&         a_divergeF,
                       LevelData<EBCellFAB>&         a_cons,
                       EBFluxRegister&               a_fineFluxRegister,
                       EBFluxRegister&               a_coarFluxRegister,
                       const LevelData<EBCellFAB>&   a_source,
                       Real a_time, Real a_dt);

  void doIrregularUpdate(LevelData<EBCellFAB>&         a_divergeF,
                         LevelData<EBCellFAB>&         a_cons,
                         EBFluxRegister&               a_fineFluxRegister,
                         EBFluxRegister&               a_coarFluxRegister,
                         LevelData<BaseIVFAB<Real> >&  a_massDiff,
                         Real a_time, Real a_dt);

  //these are not grown by one.
  LayoutData<IntVectSet> m_irregSetsSmall;

  //these are grown by one in the directions != idir
  LayoutData<IntVectSet> m_irregSetsGrown[SpaceDim];

  LevelData<BaseIFFAB<Real> > m_fluxInterpolants[SpaceDim];
  LevelData<BaseIVFAB<Real> > m_nonConsDivergence;
  LevelData<BaseIVFAB<Real> > m_ebIrregFaceFlux;
  LayoutData<IntVectSet>      m_cfIVS;
  bool               m_hasCoarser;
  bool               m_hasFiner;
  DisjointBoxLayout  m_thisGrids;
  DisjointBoxLayout  m_coarGrids;
  EBISLayout         m_thisEBISL;
  EBISLayout         m_coarEBISL;
  bool               m_isDefined;
  EBPatchReactive*   m_ebPatchReactive;
  RealVect           m_dx;
  ProblemDomain      m_domain;
  int                m_refRatCrse;
  int                m_nCons;
  int                m_nFlux;
  int                m_nGhost;
  EBPWLFillPatch     m_patcher;
  bool               m_useMassRedist;
  EBLevelRedist      m_ebLevelRedist;
  static int         s_timestep;
  bool m_doSmushing;
  bool m_doRZCoords;
  bool m_hasSourceTerm;
  LevelData<EBCellFAB> m_flattening;
  LevelData<EBCellFAB>* m_qbar;
  bool m_SFD;
  Real m_chi;
  bool m_isHyperbolicSrcSet;
  LevelData<EBCellFAB>* m_hyperbolicSrc;
private:
  //disallowed for all the usual reasons
  void operator=(const EBLevelReactive& a_input)
  {
    MayDay::Error("invalid operator");
  }
  EBLevelReactive(const EBLevelReactive& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

//#include "NamespaceFooter.H"
#endif

