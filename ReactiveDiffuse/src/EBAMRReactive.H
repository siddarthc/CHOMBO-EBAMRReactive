#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBAMRREACTIVE_H_
#define _EBAMRREACTIVE_H_

#include "AMRLevel.H"
#include "EBCellFAB.H"
#include "BaseIVFAB.H"
#include "LevelData.H"
#include "EBLevelRedist.H"
#include "EBCoarToCoarRedist.H"
#include "EBCoarToFineRedist.H"
#include "EBFineToCoarRedist.H"
#include "EBCoarseAverage.H"
#include "EBPWLFineInterp.H"
#include "EBFluxRegister.H"
#include "Box.H"
#include "IntVectSet.H"
#include "Vector.H"
#include "DisjointBoxLayout.H"
//#include "NamespaceHeader.H"

#include "EBLevelReactive.H"

#include "EBConductivityOpFactory.H"
#include "EBViscousTensorOpFactory.H"
#include "EBLevelTGA.H"
#include "MomentumTGA.H"
#include "BiCGStabSolver.H"
#include "EBQuadCFInterp.H"
#include "BaseDomainBC.H"
#include "BaseEBBC.H"

///interface for load balancing
typedef int(*LoadBalanceFunc)(Vector<int>&         a_procs,
                               const Vector<Box>&   a_boxes,
                               const ProblemDomain& a_domain,
                               bool a_verbose);


///
/**
 */
class EBAMRReactive : public AMRLevel
{
public:

  static bool s_noEBCF;

  EBAMRReactive();

  virtual ~EBAMRReactive();

  virtual void clearSolvers();

 virtual void define(AMRLevel*  a_coarser_level_ptr,
              const Box& a_problem_domain,
              int        a_level,
              int        a_ref_ratio);

  virtual void define(AMRLevel*  a_coarser_level_ptr,
              const ProblemDomain& a_problem_domain,
              int        a_level,
              int        a_ref_ratio);


  /// object for physics integration
  virtual void patchReactive(const EBPatchReactiveFactory* const a_ebPatchReactiveFactory);

  /// CFL
  virtual void CFL(Real a_cfl);

  ///
  virtual void doSmushing(bool a_doSmushing);

  ///
  virtual   void doRZCoords(bool a_doRZCoords);

  ///
  virtual   void hasSourceTerm(bool a_hasSourceTerm);

  ///redistribution radius
  virtual   void redistRadius(int a_redistRad);

  /// longest length of the domain in real space
  virtual   void domainLength(RealVect a_domainLength);

  /// refinement threshold
  virtual   void refinementThreshold(Real a_refine_thresh);

  /// tag buffer size
  virtual   void tagBufferSize(int a_tag_buffer_size);

  virtual   void useMassRedistribution(bool a_useMassRedist);

  void tagAll(bool a_tagAll);

  void addReactionRates(bool a_addReactionRates);

  void addDiffusion(bool a_addDiffusion);

  static void setLoadBalance(LoadBalanceFunc a_func)
  {
    s_loadBalance      = a_func;
    s_isLoadBalanceSet = true;
  }

  /// initialize grids
  virtual void initialGrid(const Vector<Box>& a_new_grids);

  virtual void levelSetup();

//  EBAMRReactive* getCoarserLevel() const;
//  EBAMRReactive* getFinerLevel() const;

 /// initialize data
  virtual void initialData(); 

 /// fill all ghost cells on all levels accross all processors
  void fillGhosts();

  /// create tags at initialization
  virtual void tagCellsInit(IntVectSet& a_tags);

  /// create tags
  virtual void tagCells(IntVectSet& a_tags) ;
  
  /// things to do after initialization
  virtual void postInitialize(); 

  void postInitialGrid(const bool a_restart);

  void getHierarchyAndGrids(Vector<EBAMRReactive*>&        a_hierarchy,
                            Vector<DisjointBoxLayout>&   a_grids,
                            Vector<EBISLayout>&          a_ebisl,
                            Vector<EBLevelGrid>&         a_eblg,
                            Vector<int>&                 a_refRat,
                            ProblemDomain&               a_lev0Dom,
                            Real&                        a_lev0Dx,
                            Vector<ProblemDomain>&       a_domains,
                            Vector<Real>&                a_dxs);  

  /// setup fine->coarse stuff
  /**
     stuff that needs to be setup from the finer
     level.  A bunch of objects depend upon the layouts
     from both levels and the finer level changes more
     often from regrid so this needs to be called from the finer
     level
  */
  virtual void syncWithFineLevel();

  /// compute dt with initial data
  virtual Real computeInitialDt();

  /// regrid
  virtual void regrid(const Vector<Box>& a_new_grids);

  /// advance by one timestep
  virtual Real advance();

  void dumpDebug();
  
  void dumpDebug(const string& a_debstring);

  void sumConserved(Real& a_sumcons,
                    const int& a_ivar) const;

  /// things to do after a timestep
  virtual void postTimeStep();

  ///
  virtual Real getDt() const;

  //! get the LAST dt
  Real getDtOld()
  {
    return m_dtOld;
  }

  /// compute dt
  virtual Real computeDt();

  bool addReactionRates();
 
  bool addDiffusion();
 
  // state vector at new time
  LevelData<EBCellFAB>& getStateNew();
 
  // state vector at old time
  LevelData<EBCellFAB>& getStateOld();

  bool convergedToSteadyState();

  void defineSolvers();

  void getPrimState(LevelData<EBCellFAB>&       a_prim,
                    const LevelData<EBCellFAB>& a_cons);

  void setDiffusionCoefficients(const LevelData<EBCellFAB>& a_presCell,
                                const LevelData<EBCellFAB>& a_tempCell,
                                const LevelData<EBCellFAB>& a_specDensCell,
                                const LevelData<EBFluxFAB>& a_presFace,
                                const LevelData<EBFluxFAB>& a_tempFace,
                                const LevelData<EBFluxFAB>& a_specDensFace);

  //this one is explicit and also does reweighting of the 
  //redistribution object in the case of mass weighting
  void hyperbolicRedistribution(LevelData<EBCellFAB>& a_state);

  ///inviscid version of the algorithm
  void explicitAdvance(const LevelData<EBCellFAB>& a_divergeF);

  void hyperbolicSource(LevelData<EBCellFAB>& a_source);


  // advances the solution using the hyperbolic terms.
  //(\kappa \rho^\npo I - \frac{\kappa \dt}{2} L^m)\ubar^\npo = \kappa \rho^\npo \ubold^*
  void  getUStar(LevelData<EBCellFAB>      & a_UStar,
                 const LevelData<EBCellFAB>& a_UN,
                 const LevelData<EBCellFAB>& a_divergef);

  //if doNormalization = true, averages with neighboring cells
  //otherwise returns kappa * source
  void explicitHyperbolicSource(LevelData<EBCellFAB>&  a_specMassSrc,
                                LevelData<EBCellFAB>&  a_momenSrc,
                                LevelData<EBCellFAB>&  a_energSrc,
                                const LevelData<EBCellFAB>& a_state,
                                bool a_doNormalization);

  void kappaSpecMassDiffSrc(LevelData<EBCellFAB>& a_kappaSpecMassSrc,
                      const LevelData<EBCellFAB>& a_specMF,
                      const LevelData<EBCellFAB>* a_specMFCoar,
                      const LevelData<EBCellFAB>& a_state);

  void getMCDiffTerm(LevelData<EBCellFAB>& a_MCDiffTerm,
                     const LevelData<EBCellFAB>& a_massFrac, 
                     const LevelData<EBCellFAB>* a_massFracCoar,
                     const Real a_alpha, 
                     const Real a_beta,   
                     const int a_iSpec,
                     bool  a_explicitHyperbolicSrc, 
                     bool  a_applyBC);

  void kappaMomentumSrc(LevelData<EBCellFAB>& a_kappaDivSigma,
                        const LevelData<EBCellFAB>& a_velocity,
                        const LevelData<EBCellFAB>* a_veloCoar,
                        const LevelData<EBCellFAB>& a_state);

  void kappaEnergySrc(LevelData<EBCellFAB>& a_kappaEnergySrc,
                      const LevelData<EBCellFAB>& a_velocity,
                      const LevelData<EBCellFAB>* a_veloCoar,
                      const LevelData<EBCellFAB>& a_temperature,
                      const LevelData<EBCellFAB>* a_tempCoar,
                      const LevelData<EBCellFAB>& a_specMF,
                      const LevelData<EBCellFAB>& a_state);

  void getCoarserPrimitives(LevelData<EBCellFAB>* &  a_specMFCoar,
                            LevelData<EBCellFAB>* &  a_veloCoar, 
                            LevelData<EBCellFAB>* &  a_tempCoar);

  void addMCSpecDiff(LevelData<EBCellFAB>&   a_UStar); 

  void addViscosity(LevelData<EBCellFAB>&    a_UStar);

  void addConductivity(LevelData<EBCellFAB>& a_UStar);

  void
  getSingleLdOfU(LevelData<EBCellFAB>      & a_divSigmaU,
                 LevelData<EBCellFAB>      & UStar);

    ///
  /**
     (\rho E)^{**} = (\rho E)^* + \dt L^d(U^*)
     mass diff = kappa(1-kappa)*(kappaConsDissFcn - nonConsDissFcn)
  */
  void
  updateEnergyBySingleLdAndRedistribute(LevelData<EBCellFAB>      & a_divSigmaU,
                                        LevelData<EBCellFAB>      & a_kappaConsDivSigmaU,
                                        LevelData<EBCellFAB>      & a_nonConsDivSigmaU,
                                        LevelData<EBCellFAB>      & a_UStar);

  //
  /** (\kappa \rho^\npo C_v I - \frac{\kappa \dt} L^k)T^\npo = \kappa \rho^\npo C_v T^{**}
      (\rho E)^{n+1} = (\rho E)^{**} + \dt L^k(T**)
  **/
  void getDivKappaGradT(LevelData<EBCellFAB>& a_divKappaGradT,
                        LevelData<EBCellFAB>& a_Ustar);

  //put ustar into state
  void finalAdvance(LevelData<EBCellFAB>& a_Ustar);

  void setBCs();

#ifdef CH_USE_HDF5
  /// write checkpoint header
  virtual void writeCheckpointHeader(HDF5Handle& a_handle) const;
  /// write checkpoint data for this level
  virtual void writeCheckpointLevel(HDF5Handle& a_handle) const;
  /// read checkpoint header
  virtual void readCheckpointHeader(HDF5Handle& a_handle);
  /// read checkpoint data for this level
  virtual void readCheckpointLevel(HDF5Handle& a_handle);
  /// write plotfile header
  virtual void writePlotHeader(HDF5Handle& a_handle) const;
  void writePlotHeaderOld(HDF5Handle& a_handle) const;
  /// write plot file data for this level
  virtual void writePlotLevel(HDF5Handle& a_handle) const;
  void writePlotLevelOld(HDF5Handle& a_handle) const;
  static int s_NewPlotFile;
  static bool s_solversDefined;
#endif

protected:

  void resetWeights();

  void defineFactories(bool a_atHalfTime);

  void implicitReflux();

  void refluxRedistInteraction();  

  // This performs the redistribution of a quantity (mass, momentum, or energy) in m_stateNew. 
  // The quantity is determined by a_interv. If a_putIntoState is set, the redistribution is 
  // done EXPLICITLY, and the redistributed quantity is placed directly into m_stateNew.
  // Otherwise, the quantity is accumulated into m_redisRHS for the later implicit 
  // redistribution.
  void refluxRHSConserved();
  void explicitReflux(const Interval& a_interv);
  void coarseFineRedistribution(const Interval& a_interv);

  int getFinestLevel();

  void getHalfState(LevelData<EBCellFAB>& a_stateInt);

  static LoadBalanceFunc   s_loadBalance;
  static bool              s_isLoadBalanceSet;
  bool m_tagAll;
  bool m_useMassRedist;
  bool m_addReactionRates;
  bool m_addDiffusion;
  Box m_domainBox;
  Real calculateMass() const;

  EBAMRReactive* getCoarserLevel() const;
  EBAMRReactive* getFinerLevel() const;

  DisjointBoxLayout m_grids;
  EBISLayout m_ebisl;

  // grid spacing
  RealVect m_origin;
  RealVect m_dx;
  RealVect m_aspect;

  // domain length
  RealVect m_domainLength;

  // cfl number
  Real m_cfl;


  // refinement threshold for gradient
  Real m_refineThresh;

  // tag buffer size
  int m_tagBufferSize;
  int m_redistRad;

  bool m_doRZCoords;
  bool m_hasSourceTerm;
  bool m_doSmushing;
  bool m_hasCoarser;
  bool m_hasFiner;
  bool m_isDefined;

  // state vector at old time
  LevelData<EBCellFAB> m_stateOld;

  // state vector at new time
  LevelData<EBCellFAB> m_stateNew;


  Real m_dtOld;
  // number of components of m_state
  int m_nComp;

  // number of primitives
  int m_nPrim;

  // number of species
  int m_nSpec;

  int m_nGhost ;
  // names of components
  Vector<string> m_stateNames;
  Vector<string> m_primNames;

  // new time step
  Real m_dtNew;

  // averaging from coarse to fine level
  EBCoarseAverage m_ebCoarseAverage;

  // interpolation from fine to coarse level
  EBPWLFineInterp m_ebFineInterp;

 // level solver
  EBLevelReactive m_ebLevelReactive;

  // flux register
  EBFluxRegister m_ebFluxRegister;

  // patch integrator
  EBPatchReactive* m_ebPatchReactive;

  // patch integrator
  const EBPatchReactiveFactory* m_ebPatchReactiveFactory;

  //redistribution mass
  LevelData<BaseIVFAB<Real> > m_massDiff;
  LayoutData<IntVectSet> m_sets;

  EBLevelGrid    m_eblg;

  RedistStencil m_redStencil;
  EBCoarToFineRedist m_ebCoarToFineRedist;
  EBCoarToCoarRedist m_ebCoarToCoarRedist;
  EBFineToCoarRedist m_ebFineToCoarRedist;

  EBLevelRedist m_ebLevelRedist;

//  EBFluxRegister         m_divFFluxRegister;
  EBFluxRegister  m_specFluxRegister;
  EBFluxRegister  m_veloFluxRegister;
  EBFluxRegister  m_tempFluxRegister;

  Vector<RefCountedPtr<LevelData<EBCellFAB> > >         m_acoDiff; //better to have an aco for each species since we are 
                                                                   //dealing with ptrs
  Vector<RefCountedPtr<LevelData<EBFluxFAB> > >         m_bco;
  Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >  m_bcoIrreg;

  Vector<RefCountedPtr<LevelData<EBFluxFAB> > >         m_rhsco;
  Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >  m_rhscoIrreg;

  RefCountedPtr<LevelData<EBCellFAB> >                  m_acoVisc;
  RefCountedPtr<LevelData<EBFluxFAB> >                  m_eta;
  RefCountedPtr<LevelData<BaseIVFAB<Real> > >           m_etaIrreg;
  RefCountedPtr<LevelData<EBFluxFAB> >                  m_lambda;
  RefCountedPtr<LevelData<BaseIVFAB<Real> > >           m_lambdaIrreg;
  RefCountedPtr<LevelData<EBFluxFAB> >                  m_kappa;
  RefCountedPtr<LevelData<BaseIVFAB<Real> > >           m_kappaIrreg; 
  RefCountedPtr<LevelData<EBCellFAB> >                  m_acoCond;

  RefCountedPtr<EBQuadCFInterp>                         m_quadCFI;
  
  RefCountedPtr<BaseDomainBCFactory>                    m_specDomBC;
  RefCountedPtr<BaseDomainBCFactory>                    m_veloDomBC;
  RefCountedPtr<BaseDomainBCFactory>                    m_tempDomBC;
  RefCountedPtr<BaseDomainBCFactory>                    m_smoothDomBC;
  RefCountedPtr<BaseEBBCFactory>                        m_specEBBC;
  RefCountedPtr<BaseEBBCFactory>                        m_veloEBBC;
  RefCountedPtr<BaseEBBCFactory>                        m_tempEBBC;
  RefCountedPtr<BaseEBBCFactory>                        m_smoothEBBC; 

  bool                                                  m_doLazyRelax;  //viscTensorOp
   /*------------------------------Class Static Variables     ------------------------*/

//  static Real s_alpha, s_beta;
//  static Vector<RefCountedPtr<EBConductivityOpFactory> >                  s_diffuseOpFact;
  static Vector<RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > > > s_diffuseOpFact;
  static Vector<RefCountedPtr<EBLevelBackwardEuler> >                      s_diffuseLevBE;
  static Vector<RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB> > > >      s_diffuseAMRMG;
//  static Vector<Vector<EBConductivityOp *> >                              s_EBAMRDiffuseOps;

//  static RefCountedPtr<EBViscousTensorOpFactory>                         s_viscOpFact;
  static RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > >        s_viscOpFact;
  static RefCountedPtr<MomentumBackwardEuler>                            s_viscLevBE;
  static RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB> > >             s_viscAMRMG;
//  static Vector<EBViscousTensorOp *>                                     s_EBAMRviscOps;

//  static RefCountedPtr<EBConductivityOpFactory>                          s_condOpFact;
  static RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > >        s_condOpFact;
  static RefCountedPtr<EBLevelBackwardEuler>                             s_condLevBE;
  static RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB> > >             s_condAMRMG;
//  static Vector<EBConductivityOp *>                                      s_EBAMRcondOps;
 
  static BiCGStabSolver<LevelData<EBCellFAB> >                          s_botSolver;

  void fillCoefficients(const LevelData<EBCellFAB>& a_state);

  void setDiffuseCoeff(Vector<EBAMRReactive*>&     a_hierarchy,
                       Vector<DisjointBoxLayout>&  a_grids,
                       Vector<EBISLayout>&         a_ebisl,
                       Vector<EBLevelGrid>&        a_eblg,
                       Vector<int>&                a_refRat,
                       ProblemDomain&              a_lev0Dom,
                       Vector<Real>&               a_dxs,
                       Vector<ProblemDomain>&      a_domains);

  void fillDiffCoeffMatrix(LevelData<EBCellFAB>&       a_diffCoeffMatrix,
                           const LevelData<EBCellFAB>& a_pressure,
                           const LevelData<EBCellFAB>& a_temperature,
                           const LevelData<EBCellFAB>& a_specMassFrac,
                           const EBISLayout&           a_ebisl);

  void getEBDiffuseOpFactory(RefCountedPtr<EBConductivityOpFactory>&        a_factory,
                             const Vector<DisjointBoxLayout>&               a_grids,
                             const Vector<EBISLayout>&                      a_ebisl,
                             const Vector<EBLevelGrid>&                     a_eblg,
                             Vector<int>&                                   a_refRat,
                             Real&                                          a_lev0Dx,
                             Vector<ProblemDomain>&                         a_domains,
                             const int&                                     a_iSpec);

  void getSpeciesBCFactories(RefCountedPtr<BaseDomainBCFactory>&       a_domBC,
                           RefCountedPtr<BaseEBBCFactory>&             a_ebBC);

  void setViscousAndConductiveCoeff(Vector<EBAMRReactive*>&     a_hierarchy,
                                  Vector<DisjointBoxLayout>&    a_grids,
                                  Vector<EBISLayout>&           a_ebisl,
                                  Vector<EBLevelGrid>&          a_eblg,
                                  Vector<int>&                  a_refRat,
                                  ProblemDomain&                a_lev0Dom,
                                  Vector<Real>&                 a_dxs,
                                  Vector<ProblemDomain>&        a_domains);

  void fillViscousAndConductiveCoeff(LevelData<EBCellFAB>&         a_aco,
                                     LevelData<EBCellFAB>&         a_muCell,
                                     LevelData<EBCellFAB>&         a_lambdaCell,
                                     LevelData<EBCellFAB>&         a_kappaCell,
                                     const LevelData<EBCellFAB>&   a_massFrac,
                                     const LevelData<EBCellFAB>&   a_dense,
                                     const LevelData<EBCellFAB>&   a_temperature,
                                     const EBISLayout&             a_ebisl);

  // This performs the redistribution of a quantity (mass, momentum, or energy) in m_stateNew. 
  // The quantity is determined by a_interv. If a_putIntoState is set, the redistribution is 
  // done EXPLICITLY, and the redistributed quantity is placed directly into m_stateNew.
  // Otherwise, the quantity is accumulated into m_redisRHS for the later implicit 
  // redistribution.
  
  LevelData<EBCellFAB>       m_redisRHS; 

  void postTimeStepRefluxRedistDance();

  //(rho I - dt Ly) delta = dt*Dr*(Frho)
  void getRefluxDeltaMassFrac(Vector<LevelData<EBCellFAB>* >& a_deltaMassFrac,
                              Vector<LevelData<EBCellFAB>* >& a_dtReDivergeRHO,
                              int baseLev, int finestLev, Real baseDt);

  //(rho I - dt Lv) delta = dt*Dr(Fm)
  void getRefluxDeltaV(Vector<LevelData<EBCellFAB>* >&  a_deltaVelocity,
                       Vector<LevelData<EBCellFAB>* >&  a_dtRefluxDivergeM,
                       int baseLev, int finestLev, Real baseDt);


  //(rho Cv I - dt Lt) delta = dt*Dr(Fe)
  void getRefluxDeltaT(Vector<LevelData<EBCellFAB>* >& a_deltaTemperat,
                       Vector<LevelData<EBCellFAB>* >& a_dtRefluxDivergeE,
                       int baseLev, int finestLev, Real baseDt);

   //rhoi += dt*rho*(deltaY)
   void incrDenseByDeltaY(Vector<LevelData<EBCellFAB>* >& a_deltaMassFrac,
                          int baseLev, int finestLev);
 
   // rhoi += dt*rho*(deltaY)
   void incrDenseByDeltaY(LevelData<EBCellFAB>& a_state,
                          LevelData<EBCellFAB>& a_deltaMassFrac);

  //mom += dt*rho*(deltav)
  void incrMomentumByDeltaV(Vector<LevelData<EBCellFAB>* >& a_deltaVelocity,
                            int baseLev, int finestLev);

  //mom += dt*rho*(deltav)
  void incrMomentumByDeltaV(LevelData<EBCellFAB>& a_state,
                            LevelData<EBCellFAB>& a_deltaVelocity);

  //ene += dt*rho*Cv*(deltaT)
  void incrEnergyByDeltaT(Vector<LevelData<EBCellFAB>* >& a_deltaTemperat,
                          int baseLev, int finestLev);

  void incrEnergyByDeltaT(LevelData<EBCellFAB>& a_state,
                          LevelData<EBCellFAB>& a_deltaTemperat);

  void coarseFineIncrement();
 
private:
  //disallowed for all the usual reasons
  void operator=(const EBAMRReactive& a_input)
  {
    MayDay::Error("invalid operator");
  }
  EBAMRReactive(const EBAMRReactive& a_input)
  {
    MayDay::Error("invalid operator");
  }

};

//#include "NamespaceFooter.H"
#endif
                                 
