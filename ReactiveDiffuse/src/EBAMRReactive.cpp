#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cmath>

#include "parstream.H"
#include "ParmParse.H"
#include "PolyGeom.H"

#include "DebugOut.H"
#include "Box.H"
#include "Vector.H"
#include "IntVectSet.H"
#include "EBCellFAB.H"
#include "DisjointBoxLayout.H"
#include "LayoutIterator.H"
#include "CH_HDF5.H"
#include "SPMD.H"
#include "LoadBalance.H"
#include "ProblemDomain.H"
#include "BoxIterator.H"
#include "EBAMRIO.H"
#include "AMRLevel.H"
#include "EBCellFactory.H"
#include "BaseIVFactory.H"
#include "VoFIterator.H"
#include "EBIndexSpace.H"
#include "EBArith.H"
#include "EBFastFR.H"
#include "AMR.H"
#include "FabDataOps.H"
#include "EBLevelDataOps.H"
#include "KappaSquareNormal.H"
#include "FabDataOps.H"
#include "EBAMRDataOps.H"
#include "EBAMRReactive.H"
#include "EBPatchReactiveF_F.H"
#include "EBAMRReactiveF_F.H"
#include "DirichletPoissonEBBC.H"
#include   "NeumannPoissonEBBC.H"
#include "DirichletPoissonDomainBC.H"
#include   "NeumannPoissonDomainBC.H"
#include "DirichletViscousTensorEBBC.H"
#include   "NeumannViscousTensorEBBC.H"
#include "DirichletViscousTensorDomainBC.H"
#include   "NeumannViscousTensorDomainBC.H"
#include "DirichletConductivityDomainBC.H"
#include   "NeumannConductivityDomainBC.H"
#include "DirichletConductivityEBBC.H"
#include   "NeumannConductivityEBBC.H"
#include "EBLGIntegrator.H"

bool EBAMRReactive::s_isLoadBalanceSet = false;
LoadBalanceFunc EBAMRReactive::s_loadBalance = NULL;
IntVect ivdebamrg(D_DECL(16, 5, 0));
int EBAMRReactive::s_NewPlotFile = 0;
int debuglevel = 1;

Vector<RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > > > EBAMRReactive::s_diffuseOpFact;
RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > >          EBAMRReactive::s_viscOpFact;
RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > >          EBAMRReactive::s_condOpFact; 

Vector<RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB> > > >  EBAMRReactive::s_diffuseAMRMG;
RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB> > >           EBAMRReactive::s_viscAMRMG;
RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB> > >           EBAMRReactive::s_condAMRMG;

Vector<RefCountedPtr<EBLevelBackwardEuler> >  EBAMRReactive::s_diffuseLevBE;
RefCountedPtr<MomentumBackwardEuler>          EBAMRReactive::s_viscLevBE;
RefCountedPtr<EBLevelBackwardEuler>           EBAMRReactive::s_condLevBE;  

BiCGStabSolver<LevelData<EBCellFAB> >        EBAMRReactive::s_botSolver;
bool EBAMRReactive::s_noEBCF = false;
bool EBAMRReactive::s_solversDefined = false;


/***************************/
bool
EBAMRReactive::
convergedToSteadyState()
{
  EBCellFactory      cellFact(m_eblg.getEBISL());
  LevelData<EBCellFAB>  udiff(m_eblg.getDBL(), m_stateNew.nComp(), 4*IntVect::Unit, cellFact);
  EBLevelDataOps::setToZero(udiff);
  EBLevelDataOps::incr(udiff, m_stateOld, -1.0);
  EBLevelDataOps::incr(udiff, m_stateNew,  1.0);
  Real umax, umin, eps;
  Real dmax, dmin;
  int ivar;
  ParmParse pp;
  pp.get("convergence_metric", eps);
  pp.get("convergence_variable", ivar);
  EBLevelDataOps::getMaxMin(dmax, dmin, udiff, ivar);
  EBLevelDataOps::getMaxMin(umax, umin, m_stateNew, ivar);
  Real denom = 1;
  if(Abs(umax - umin) > eps)
    {
      denom = Abs(umax-umin);
    }
  Real maxdiff = Abs(dmax-dmin)/denom;
  pout() << "max difference in convergence variable = " << maxdiff << ", eps set to " << eps << endl;
  return (maxdiff < eps);
}
/***************************/
void EBAMRReactive::getPrimState(LevelData<EBCellFAB>&        a_prim,
                                 const LevelData<EBCellFAB>&  a_cons)
{
  EBLevelDataOps::setToZero(a_prim);
  for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      const Box& region = m_eblg.getDBL().get(dit());
      BaseFab<Real>& regPrim = a_prim[dit()].getSingleValuedFAB();
      const BaseFab<Real>& regCons = a_cons[dit()].getSingleValuedFAB();

      int logflag = 0;
      int verbose = 0;
      int sss = 1;
 
      FORT_CONS2PRM(CHF_BOX(region),
                    CHF_CONST_FRA(regCons),
                    CHF_FRA(regPrim),
                    CHF_CONST_INT(logflag),
                    CHF_CONST_INT(verbose));

      IntVectSet ivsMulti = m_eblg.getEBISL()[dit()].getMultiCells(region);
      for(VoFIterator vofit(ivsMulti, m_eblg.getEBISL()[dit()].getEBGraph()); vofit.ok(); ++vofit)
        {
          Vector<Real> cvec(m_nComp); Vector<Real> pvec(m_nPrim);
          for(int ivar = 0; ivar< m_nComp; ivar++)
            {
              cvec[ivar] = a_cons[dit()](vofit(), ivar);
            }

          FORT_POINTCONS2PRM(CHF_VR(cvec),
                             CHF_VR(pvec),
                             CHF_INT(sss),
                             CHF_INT(logflag));

          for(int ivar = 0; ivar < m_nPrim; ivar++) 
           {
             a_prim[dit()](vofit(), ivar) = pvec[ivar];
           }

        }
    }
}
/***************************/
 void 
EBAMRReactive::
addReactionRates(bool a_addReactionRates)
{
  m_addReactionRates = a_addReactionRates;
}
/***************************/ 
void 
EBAMRReactive::
addDiffusion(bool a_addDiffusion)
{
  m_addDiffusion = a_addDiffusion;
}
/***************************/
void
EBAMRReactive::
tagAll(bool a_tagAll)
{
  m_tagAll = a_tagAll;
}
/***************************/
void
EBAMRReactive::
doSmushing(bool a_doSmushing)
{
  m_doSmushing = a_doSmushing;
}
/***************************/
void
EBAMRReactive::
doRZCoords(bool a_doRZCoords)
{
  m_doRZCoords = a_doRZCoords;
}
/***************************/
void
EBAMRReactive::
hasSourceTerm(bool a_hasSourceTerm)
{
  m_hasSourceTerm = a_hasSourceTerm;
}
/***************************/
void
EBAMRReactive::
useMassRedistribution(bool a_useMassRedist)
{
  m_useMassRedist = a_useMassRedist;
}
/***************************/
EBAMRReactive::EBAMRReactive()
{
  m_cfl = 0.8;
  m_tagAll = false;
  m_useMassRedist = true;
  m_doRZCoords = false;
  m_hasSourceTerm = false;
  m_doSmushing = true;
  m_origin = RealVect::Zero;
  m_dx = RealVect::Unit;
  m_aspect = RealVect::Unit;
  m_domainLength = RealVect::Unit;
  m_refineThresh = 0.2;
  m_initial_dt_multiplier = 0.1;
  m_ebPatchReactive = NULL;
  m_redistRad = 1;
  m_isDefined = false;
  m_addReactionRates = false;
  m_addDiffusion = false;
}
/***************************/
void EBAMRReactive::redistRadius(int a_redistRad)
{
  m_redistRad = a_redistRad;
}
/***************************/
EBAMRReactive::~EBAMRReactive()
{
  if (m_ebPatchReactive != NULL)
    delete m_ebPatchReactive;

  if ((m_level == 0) && addDiffusion())
    { 
      clearSolvers();
    }
}
/***************************/
void EBAMRReactive::clearSolvers()
{
  for (int iSpec = 0; iSpec < m_nSpec; iSpec++)
   {
     s_diffuseOpFact[iSpec] = RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > >();
     s_diffuseAMRMG[iSpec] = RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB> > >();
     s_diffuseLevBE[iSpec] = RefCountedPtr<EBLevelBackwardEuler>();
   }
  s_viscOpFact = RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > >();
  s_viscAMRMG = RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB> > >();
  s_viscLevBE = RefCountedPtr<MomentumBackwardEuler>(); 

  s_condOpFact = RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > >();
  s_condAMRMG = RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB> > >();
  s_condLevBE = RefCountedPtr<EBLevelBackwardEuler>();
}
/***************************/
void EBAMRReactive::define(AMRLevel*  a_coarser_level_ptr,
                          const Box& a_problem_domain,
                          int        a_level,
                          int        a_ref_ratio)
{
  MayDay::Error("EBAMRReactive::define -\n\tShould never be called with a Box for a problem domain");
}
/***************************/
void EBAMRReactive::define(AMRLevel*  a_coarser_level_ptr,
                          const ProblemDomain& a_problem_domain,
                          int        a_level,
                          int        a_ref_ratio)
{
  CH_TIME("EBAMRReactive::define");
  if (s_verbosity >= 3)
    {
      pout() << "EBAMRReactive::define, level=" << a_level << endl;
    }

  m_isDefined = true;
  AMRLevel::define(a_coarser_level_ptr,
                   a_problem_domain,
                   a_level,
                   a_ref_ratio);
  m_domainBox = m_problem_domain.domainBox();
 
  if (a_coarser_level_ptr != NULL)
    {
      EBAMRReactive* amrr_ptr =
        dynamic_cast<EBAMRReactive*>(a_coarser_level_ptr);
      if (amrr_ptr == NULL)
        {
          pout() << "EBAMRR::define:cannot cast  to derived class"
                 << endl;
          MayDay::Error();
        }

      m_cfl = amrr_ptr->m_cfl;
      m_domainLength = amrr_ptr->m_domainLength;
      m_refineThresh = amrr_ptr->m_refineThresh;
      m_tagBufferSize = amrr_ptr->m_tagBufferSize;
    }
  for (int idir = 0; idir < SpaceDim; ++idir)
    {
      m_dx[idir] = m_domainLength[idir]/m_domainBox.size(idir);
    }
  m_nGhost = 4; 

  if (m_ebPatchReactive != NULL)
    delete m_ebPatchReactive;

  m_ebPatchReactive = m_ebPatchReactiveFactory->create();
  m_ebPatchReactive->define(m_problem_domain, m_dx);

  m_nComp = m_ebPatchReactive->numConserved();
  m_nSpec = m_ebPatchReactive->getnSpecies();
  m_nPrim = m_ebPatchReactive->numPrimitives();
  m_stateNames = m_ebPatchReactive->stateNames();
  m_primNames = m_ebPatchReactive->primNames();
}

/***************************/
void EBAMRReactive::patchReactive(const EBPatchReactiveFactory* const a_ebPatchReactiveFactory)
{
  m_ebPatchReactiveFactory = a_ebPatchReactiveFactory;
}
/***************************/
void EBAMRReactive::CFL(Real a_cfl)
{
  m_cfl = a_cfl;
}
/***************************/
void EBAMRReactive::domainLength(RealVect a_domainLength)
{
  m_domainLength = a_domainLength;
}
/***************************/
void EBAMRReactive::refinementThreshold(Real a_refineThresh)
{
  m_refineThresh = a_refineThresh;
}
/***************************/
void EBAMRReactive::tagBufferSize(int a_tagBufferSize)
{
  m_tagBufferSize = a_tagBufferSize;
}
/***************************/
bool EBAMRReactive::addReactionRates()
{
  CH_assert(m_isDefined);
  return m_addReactionRates;
}
/***************************/
bool EBAMRReactive::addDiffusion()
{
  CH_assert(m_isDefined);
  return m_addDiffusion;
} 
/***************************/
void EBAMRReactive::initialGrid(const Vector<Box>& a_new_grids)
{
  CH_TIME("EBAMRR::initialGrid");

  if (s_verbosity >= 3)
    {
      pout() << " in EBAMRReactive initialGrid for level " << m_level << endl;
    }
  Vector<Box>& newGrids =   (Vector<Box>&) a_new_grids;
  mortonOrdering(newGrids);
  const EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();
  m_level_grids = a_new_grids;

  // load balance and create boxlayout
  Vector<int> proc_map;
  if (s_isLoadBalanceSet)
    {
      s_loadBalance(proc_map,a_new_grids, m_domainBox, false);
    }
  else
    {
      LoadBalance(proc_map,a_new_grids);
    }
  if (s_verbosity >= 3)
    {
      pout() << " just loadbalanced " << m_level << endl;
    }

  m_grids = DisjointBoxLayout(a_new_grids,proc_map);
  if (s_verbosity >= 3)
    {
      pout() << "EBAMRReactive::initialgrid grids " << endl;
      dumpDBL(&m_grids);
    }

  //the ebisl has to know about the fact that we really have
  //four ghost cells.
  int nGhostEBISL = 6;
  if (s_verbosity >= 3)
    {
      pout() << " about to fill ebislayout  in EBAMRReactive initialGrid for level " << m_level << endl;
    }

  ebisPtr->fillEBISLayout(m_ebisl, m_grids,  m_domainBox, nGhostEBISL);

  if (s_verbosity >= 3)
    {
      pout() << " done with filling ebislayout  in EBAMRReactive initialGrid for level " << m_level << endl;
    }

  EBCellFactory factoryNew(m_ebisl);
  IntVect ivGhost = m_nGhost*IntVect::Unit;
  m_stateNew.define(m_grids,m_nComp, ivGhost, factoryNew);
  m_stateOld.define(m_grids,m_nComp, ivGhost, factoryNew);

  m_eblg.define(m_grids, m_problem_domain, m_nGhost, ebisPtr);

  // set up data structures
  levelSetup();
}
/***************************/
void
EBAMRReactive::levelSetup()
{
  CH_TIME("EBAMRG::levelSetup");
  if (s_verbosity >= 3)
    {
      pout() << " in EBAMRReactive levelSetup for level " << m_level << endl;
    }

  EBCellFactory factoryNew(m_eblg.getEBISL());
  IntVect ivGhost = m_nGhost*IntVect::Unit;
  m_stateNew.define(m_eblg.getDBL(),m_nComp, ivGhost, factoryNew);
  m_stateOld.define(m_eblg.getDBL(),m_nComp, ivGhost, factoryNew);
  m_redisRHS.define(m_eblg.getDBL(),m_nComp, ivGhost, factoryNew);
  EBLevelDataOps::setVal(m_redisRHS, 0.0);

  m_sets.define(m_eblg.getDBL());
  for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      m_sets[dit()] = m_eblg.getEBISL()[dit()].getIrregIVS(m_eblg.getDBL().get(dit()));
    }

  EBCellFactory       cellFact(m_eblg.getEBISL());
  EBFluxFactory       fluxFact(m_eblg.getEBISL());
  BaseIVFactory<Real> bivfFact(m_eblg.getEBISL(), m_sets);

  int nghost = 4;
  
  // Allocate the A coefficients 
  m_acoVisc = RefCountedPtr< LevelData<EBCellFAB> >       (new LevelData<EBCellFAB>       (m_eblg.getDBL(), 1, nghost*IntVect::Unit, cellFact));
  m_acoCond = RefCountedPtr< LevelData<EBCellFAB> >       (new LevelData<EBCellFAB>       (m_eblg.getDBL(), 1, nghost*IntVect::Unit, cellFact));

  // Allocate the species diffusive coefficients on regular and irregular cells
  m_bco.resize(m_nSpec);
  m_bcoIrreg.resize(m_nSpec);
  
  // Allocate the rhs here too
  m_acoDiff.resize(m_nSpec);
  m_rhsco.resize(m_nSpec);
  m_rhscoIrreg.resize(m_nSpec); 

  for (int iSpec = 0; iSpec <  m_nSpec; iSpec++)
   {
     m_acoDiff[iSpec]      = RefCountedPtr< LevelData<EBCellFAB> >       (new LevelData<EBCellFAB>       (m_eblg.getDBL(), 1, nghost*IntVect::Unit, cellFact));
     m_bco[iSpec]          = RefCountedPtr< LevelData<EBFluxFAB> >       (new LevelData<EBFluxFAB>       (m_eblg.getDBL(), 1, nghost*IntVect::Unit, fluxFact));
     m_bcoIrreg[iSpec]     = RefCountedPtr< LevelData<BaseIVFAB<Real> > >(new LevelData<BaseIVFAB<Real> >(m_eblg.getDBL(), 1, nghost*IntVect::Unit, bivfFact));
     m_rhsco[iSpec]      = RefCountedPtr< LevelData<EBFluxFAB> >       (new LevelData<EBFluxFAB>       (m_eblg.getDBL(), m_nSpec, nghost*IntVect::Unit, fluxFact));
     m_rhscoIrreg[iSpec] = RefCountedPtr< LevelData<BaseIVFAB<Real> > >(new LevelData<BaseIVFAB<Real> >(m_eblg.getDBL(), m_nSpec, nghost*IntVect::Unit, bivfFact)); 
   }  

  //Allocate the viscous tensor coefficients on regular and irregular cells
  m_eta         = RefCountedPtr< LevelData<EBFluxFAB> >       (new LevelData<EBFluxFAB>       (m_eblg.getDBL(), 1, nghost*IntVect::Unit, fluxFact));
  m_etaIrreg    = RefCountedPtr< LevelData<BaseIVFAB<Real> > >(new LevelData<BaseIVFAB<Real> >(m_eblg.getDBL(), 1, nghost*IntVect::Unit, bivfFact));
  m_lambda      = RefCountedPtr< LevelData<EBFluxFAB> >       (new LevelData<EBFluxFAB>       (m_eblg.getDBL(), 1, nghost*IntVect::Unit, fluxFact));
  m_lambdaIrreg = RefCountedPtr< LevelData<BaseIVFAB<Real> > >(new LevelData<BaseIVFAB<Real> >(m_eblg.getDBL(), 1, nghost*IntVect::Unit, bivfFact));

  // Allocate the thermal conductivity on regular and irregular cells
  m_kappa     = RefCountedPtr< LevelData<EBFluxFAB> >       (new LevelData<EBFluxFAB>       (m_eblg.getDBL(), 1, nghost*IntVect::Unit, fluxFact));
  m_kappaIrreg= RefCountedPtr< LevelData<BaseIVFAB<Real> > >(new LevelData<BaseIVFAB<Real> >(m_eblg.getDBL(), 1, nghost*IntVect::Unit, bivfFact));

  EBAMRReactive* coarPtr = getCoarserLevel();
  EBAMRReactive* finePtr = getFinerLevel();

  m_hasCoarser = (coarPtr != NULL);
  m_hasFiner   = (finePtr != NULL);

  //define redistribution object for this level
  //for now set to volume weighting
  m_ebLevelRedist.define(m_eblg.getDBL(),
                         m_eblg.getEBISL(),
                         m_eblg.getDomain(),
                         m_nComp);

  if (m_hasCoarser)
    {
      int nRefCrse = m_coarser_level_ptr->refRatio();
      const EBLevelGrid& coEBLG = coarPtr->m_eblg;
      {
        CH_TIME("ave_interp_defs");
        m_ebCoarseAverage.define(m_eblg.getDBL(),
                                 coEBLG.getDBL(),
                                 m_eblg.getEBISL(),
                                 coEBLG.getEBISL(),
                                 coEBLG.getDomain().domainBox(),
                                 nRefCrse,
                                 m_nComp, Chombo_EBIS::instance());
        m_ebFineInterp.define(m_eblg.getDBL(),
                              coEBLG.getDBL(),
                              m_eblg.getEBISL(),
                              coEBLG.getEBISL(),
                              coEBLG.getDomain().domainBox(),
                              nRefCrse,
                              m_nComp);
      }

      // maintain levelreactive
      m_ebLevelReactive.define(m_eblg.getDBL(),
                              coEBLG.getDBL(),
                              m_eblg.getEBISL(),
                              coEBLG.getEBISL(),
                              m_eblg.getDomain(),
                              nRefCrse,
                              m_dx,
                              m_useMassRedist,
                              m_doSmushing,
                              m_doRZCoords,
                              m_hasSourceTerm,
                              m_ebPatchReactiveFactory,
                              m_hasCoarser,
                              m_hasFiner);

      //define fine to coarse redistribution object
      //for now set to volume weighting
      if (!s_noEBCF)
       {
         CH_TIME("fineToCoar_defs");
         m_ebFineToCoarRedist.define(m_eblg.getDBL(), coEBLG.getDBL(),
                                     m_eblg.getEBISL(), coEBLG.getEBISL(),
                                     coEBLG.getDomain().domainBox(), nRefCrse, m_nComp);
         m_ebFineToCoarRedist.setToZero();
       }

      int nvarQuad = 1;  // species and temperature

      m_quadCFI = RefCountedPtr<EBQuadCFInterp>
        (new EBQuadCFInterp(m_eblg.getDBL(),
                            coEBLG.getDBL(),
                            m_eblg.getEBISL(),
                            coEBLG.getEBISL(),
                            coEBLG.getDomain(),
                            nRefCrse, nvarQuad,
                            (*m_eblg.getCFIVS()),
                            Chombo_EBIS::instance()));
                           

      coarPtr->syncWithFineLevel();
    }
  else
    {
      m_quadCFI = RefCountedPtr<EBQuadCFInterp>(new EBQuadCFInterp());
      m_ebLevelReactive.define(m_eblg.getDBL(),
                              DisjointBoxLayout(),
                              m_eblg.getEBISL(),
                              EBISLayout(),
                              m_eblg.getDomain(),
                              m_ref_ratio,
                              m_dx,
                              m_useMassRedist,
                              m_doSmushing,
                              m_doRZCoords,
                              m_hasSourceTerm,
                              m_ebPatchReactiveFactory,
                              m_hasCoarser,
                              m_hasFiner);
    }
  //set up mass redistribution array
  m_sets.define(m_grids);
  for (DataIterator dit = m_grids.dataIterator();
      dit.ok(); ++dit)
    {
      Box thisBox = m_grids.get(dit());
      m_sets[dit()] = m_ebisl[dit()].getIrregIVS(thisBox);
    }
  BaseIVFactory<Real> factory(m_ebisl, m_sets);
  //the ghost cells on the mass redistribution array
  //are tied to the redistribution radius
  nghost = 2*m_redistRad;
  ivGhost = nghost*IntVect::Unit;
  m_massDiff.define(m_grids, m_nComp, ivGhost, factory);
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      m_massDiff[dit()].setVal(0.0);
    }
}
/***************************/
EBAMRReactive*
EBAMRReactive::getCoarserLevel() const
{
  EBAMRReactive* retval = NULL;
  if (m_coarser_level_ptr != NULL)
    {
      retval = dynamic_cast <EBAMRReactive*> (m_coarser_level_ptr);

      if (retval == NULL)
        {
          pout() << "EBAMRR::getCoarserLevel: dynamic cast failed"
                 << endl;
          MayDay::Error();
        }
    }
  return retval;
}
/***************************/
EBAMRReactive*
EBAMRReactive::getFinerLevel() const
{
  EBAMRReactive* retval = NULL;
  if (m_finer_level_ptr != NULL)
    {
      retval = dynamic_cast <EBAMRReactive*> (m_finer_level_ptr);

      if (retval == NULL)
        {
          pout() << "EBAMRG::getFinerLevel: dynamic cast failed"
                 << endl;
          MayDay::Error();
        }
    }
  return retval;
}
/***************************/
void EBAMRReactive::initialData()
{
  CH_TIME("EBAMRR::initialData");

    if (s_verbosity >= 3)
    {
      pout() << " in EBAMRReactive initialData for level " << m_level << endl;
    }
 
  const EBPhysIBC* const ebphysIBCPtr =
  m_ebPatchReactive->getEBPhysIBC();

  int nComp = m_stateNew.nComp()-1;
  // initialize data on the levels to be 0.0 
  EBLevelDataOps::setVal(m_stateNew,0.0);
  EBLevelDataOps::setVal(m_stateOld,0.0);

  //initialize both new and old states to
  //be the same thing

  ebphysIBCPtr->initialize(m_stateNew, m_ebisl);
  ebphysIBCPtr->initialize(m_stateOld, m_ebisl);
}  
/***************************/
void 
EBAMRReactive::
defineSolvers()
{
  CH_TIME("EBAMRReactive::defineSolvers");
  ParmParse pp;
  bool tagAllIrregular = false;
  if(pp.contains("tag_all_irregular"))
    {
      pp.get("tag_all_irregular", tagAllIrregular);
    }
  if(tagAllIrregular) s_noEBCF = true;

  EBConductivityOp::setForceNoEBCF(s_noEBCF);
  EBViscousTensorOp::setForceNoEBCF(s_noEBCF);
 
  s_diffuseOpFact.resize(m_nSpec);
  
  defineFactories(true);  // NOTE: coefficients get filled here
  
  if (addDiffusion())
   {
     Vector<AMRLevel*> hierarchy = AMRLevel::getAMRLevelHierarchy();
     int nlevels = hierarchy.size();
      Vector<int>                                           refRat(nlevels);
      Vector<EBLevelGrid>                                   eblgs(nlevels);
      Vector<DisjointBoxLayout>                             grids(nlevels);

      EBAMRReactive* coarsestLevel = dynamic_cast<EBAMRReactive*>(hierarchy[0]);
      ProblemDomain lev0Dom        =  coarsestLevel->m_eblg.getDomain();
      for(int ilev = 0; ilev < nlevels; ilev++)
        {
          EBAMRReactive* reactiveLevel = (EBAMRReactive*)(hierarchy[ilev]);
          eblgs       [ilev]           = reactiveLevel->m_eblg;
          grids       [ilev]           = reactiveLevel->m_eblg.getDBL();
          refRat      [ilev]           = reactiveLevel->m_ref_ratio;
        }

      //alpha, beta get replaced in tga solves    
      int nghost = 4;
      IntVect giv = nghost*IntVect::Unit;
      RealVect origin =  RealVect::Zero;
      int numSmooth, numMG, maxIter, mgverb;
      Real tolerance, hang, normThresh;
      ParmParse pp("amrmultigrid");
      pp.get("num_smooth", numSmooth);
      pp.get("num_mg",     numMG);
      pp.get("hang_eps",   hang);
      pp.get("norm_thresh",normThresh);
      pp.get("tolerance",  tolerance);
      pp.get("max_iter",   maxIter);
      pp.get("verbosity",  mgverb);

      Real bottomCushion = 1.0;
      if (pp.contains("bottom_cushion"))
       {
         pp.get("bottom_cushion", bottomCushion);
       }

      s_diffuseAMRMG.resize(m_nSpec);
      s_diffuseLevBE.resize(m_nSpec);

      // species stuff
      for (int iSpec = 0; iSpec < m_nSpec; iSpec++)
       {
         // multigrid solver
         s_diffuseAMRMG[iSpec] = RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB> > >(new AMRMultiGrid<LevelData<EBCellFAB> >());
         s_diffuseAMRMG[iSpec]->define(lev0Dom, *s_diffuseOpFact[iSpec], &s_botSolver, nlevels);
         s_diffuseAMRMG[iSpec]->setSolverParameters(numSmooth, numSmooth, numSmooth, numMG, maxIter, tolerance, hang, normThresh);
         s_diffuseAMRMG[iSpec]->m_verbosity = mgverb;
         s_diffuseAMRMG[iSpec]->m_bottomSolverEpsCushion = bottomCushion;

         // BE Integrator
         s_diffuseLevBE[iSpec] = RefCountedPtr<EBLevelBackwardEuler>(new EBLevelBackwardEuler(grids, refRat, lev0Dom, s_diffuseOpFact[iSpec], s_diffuseAMRMG[iSpec]));
         s_diffuseLevBE[iSpec]->setEBLG(eblgs);
       } // end species loop

      // viscous stuff
      s_viscAMRMG = RefCountedPtr<AMRMultiGrid< LevelData<EBCellFAB> > >( new AMRMultiGrid< LevelData<EBCellFAB> > ());
      s_viscAMRMG->define(lev0Dom, *s_viscOpFact, &s_botSolver, nlevels);
      s_viscAMRMG->setSolverParameters(numSmooth, numSmooth, numSmooth, numMG, maxIter, tolerance, hang, normThresh);
      s_viscAMRMG->m_verbosity = mgverb;
      s_viscAMRMG->m_bottomSolverEpsCushion = bottomCushion;
      
      s_viscLevBE = RefCountedPtr<MomentumBackwardEuler>(new MomentumBackwardEuler(grids, refRat, lev0Dom, s_viscOpFact, s_viscAMRMG));
      s_viscLevBE->setEBLG(eblgs);

      // conductivity stuff
      s_condAMRMG = RefCountedPtr<AMRMultiGrid< LevelData<EBCellFAB> > >( new AMRMultiGrid< LevelData<EBCellFAB> > ());
      s_condAMRMG->define(lev0Dom, *s_condOpFact, &s_botSolver, nlevels);
      s_condAMRMG->setSolverParameters(numSmooth, numSmooth, numSmooth, numMG, maxIter, tolerance, hang, normThresh);
      s_condAMRMG->m_verbosity = mgverb;
      s_condAMRMG->m_bottomSolverEpsCushion = bottomCushion;
       
      s_condLevBE = RefCountedPtr<EBLevelBackwardEuler>(new EBLevelBackwardEuler(grids, refRat, lev0Dom, s_condOpFact, s_condAMRMG));
      s_condLevBE->setEBLG(eblgs);

   } // end addDiffusion  
}
/***************************/
void 
EBAMRReactive::
defineFactories(bool a_atHalfTime)
{
  CH_TIME("EBAMRReactive::defineFactories");
  if (addDiffusion())
    {
      Vector<AMRLevel*> hierarchy = AMRLevel::getAMRLevelHierarchy();
      int nlevels = hierarchy.size();
      Vector<int>                                           refRat(nlevels);
      Vector<EBLevelGrid>                                   eblgs(nlevels);
      Vector<DisjointBoxLayout>                             grids(nlevels);
      Vector<RefCountedPtr<LevelData<EBCellFAB> > >         acoVisc(nlevels);
      Vector<RefCountedPtr<LevelData<EBCellFAB> > >         acoCond(nlevels);
      Vector<RefCountedPtr<LevelData<EBFluxFAB> > >         eta(nlevels);
      Vector<RefCountedPtr<LevelData<EBFluxFAB> > >         lambda(nlevels);
      Vector<RefCountedPtr<LevelData<EBFluxFAB> > >         kappa(nlevels);
      Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >  etaIrreg(nlevels);
      Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >  lambdaIrreg(nlevels);
      Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >  kappaIrreg(nlevels);
      Vector<RefCountedPtr<EBQuadCFInterp> >                quadCFI(nlevels);

      Vector<Vector<RefCountedPtr<LevelData<EBCellFAB> > > >        acoDiff(m_nSpec);
      Vector<Vector<RefCountedPtr<LevelData<EBFluxFAB> > > >        bco(m_nSpec);
      Vector<Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > > > bcoIrreg(m_nSpec);

      EBAMRReactive* coarsestLevel = (EBAMRReactive*)(hierarchy[0]);
      Real           lev0Dx        = (coarsestLevel->m_dx[0]);
      ProblemDomain  lev0Dom       = coarsestLevel->m_eblg.getDomain();

      for (int iSpec = 0; iSpec < m_nSpec; iSpec++)
       {
         acoDiff[iSpec].resize(nlevels);
         bco[iSpec].resize(nlevels);
         bcoIrreg[iSpec].resize(nlevels);
       }

      for (int ilev = 0; ilev < nlevels; ilev++)
       {
         EBAMRReactive* reactiveLevel = dynamic_cast<EBAMRReactive*>(hierarchy[ilev]);
         EBCellFactory fact(reactiveLevel->m_eblg.getEBISL());
         LevelData<EBCellFAB> halfSt(reactiveLevel->m_eblg.getDBL(), m_nComp, 4*IntVect::Unit, fact);
         if (a_atHalfTime)
           {
             reactiveLevel->getHalfState(halfSt);
           }
         else
           {
             Interval interv(0,m_nComp-1);
             reactiveLevel->m_stateOld.copyTo(interv, halfSt, interv);
           }

         reactiveLevel->fillCoefficients(halfSt);

         eblgs      [ilev] = reactiveLevel->m_eblg;
         grids      [ilev] = reactiveLevel->m_eblg.getDBL();
         refRat     [ilev] = reactiveLevel->m_ref_ratio;
         acoVisc    [ilev] = reactiveLevel->m_acoVisc;
         acoCond    [ilev] = reactiveLevel->m_acoCond;
         eta        [ilev] = reactiveLevel->m_eta;
         etaIrreg   [ilev] = reactiveLevel->m_etaIrreg;
         lambda     [ilev] = reactiveLevel->m_lambda;
         lambdaIrreg[ilev] = reactiveLevel->m_lambdaIrreg;
         kappa      [ilev] = reactiveLevel->m_kappa;
         kappaIrreg [ilev] = reactiveLevel->m_kappaIrreg;
         quadCFI    [ilev] = reactiveLevel->m_quadCFI;
         
         for (int iSpec = 0; iSpec < m_nSpec; iSpec++)
          {
            acoDiff[iSpec][ilev] = reactiveLevel->m_acoDiff[iSpec];
            bco[iSpec][ilev] = reactiveLevel->m_bco[iSpec];
            bcoIrreg[iSpec][ilev] = reactiveLevel->m_bcoIrreg[iSpec];
          } // spec loop
       } // end level loop

      //alpha, beta get replaced in solvers
      Real alpha = 1;
      Real beta = 1;
      int nghost = 4;
      IntVect giv = nghost*IntVect::Unit;
 
      setBCs();  // sets BCs for diffusive solvers
 
      //species diffusion operator
      int relaxType = 0;  // not sure!
      for (int iSpec = 0; iSpec < m_nSpec; iSpec++)
       {
         s_diffuseOpFact[iSpec] =RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > > 
                                 (dynamic_cast<AMRLevelOpFactory<LevelData<EBCellFAB> >*>
                                  (new EBConductivityOpFactory(eblgs, quadCFI, alpha, beta,
                                  acoDiff[iSpec], bco[iSpec], bcoIrreg[iSpec], lev0Dx, refRat, 
                                  m_specDomBC, m_specEBBC, giv, giv, relaxType)));
       }  // end species loop

      //Viscous tensor operator
      bool noMG = true;
      m_doLazyRelax = true;
      s_viscOpFact = 
        RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > >
        (dynamic_cast<AMRLevelOpFactory<LevelData<EBCellFAB> >*>
         (new EBViscousTensorOpFactory(eblgs, alpha, beta, acoVisc, eta, 
                                       lambda, etaIrreg, lambdaIrreg,lev0Dx, refRat,
                                       m_veloDomBC, m_veloEBBC, giv, giv, -1, noMG))); 
       
      EBViscousTensorOp::doLazyRelax(m_doLazyRelax);

      //Thermal dissusion operator
      s_condOpFact = 
        RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > >
        (dynamic_cast<AMRLevelOpFactory<LevelData<EBCellFAB> >*>
         (new EBConductivityOpFactory(eblgs, quadCFI, alpha, beta, acoCond,
                                      kappa, kappaIrreg, lev0Dx,           
                                      refRat, m_tempDomBC, m_tempEBBC,                                                                                          
                                      giv, giv, relaxType)));

      m_eta->exchange(Interval(0,0));
      m_lambda->exchange(Interval(0,0));
      m_kappa->exchange(Interval(0,0));
      for (int iSpec = 0; iSpec < m_nSpec; iSpec++)
       {
         m_bco[iSpec]->exchange(Interval(0,0));
         m_rhsco[iSpec]->exchange(Interval(0,m_nSpec-1));
       }
    } // end addDiffusion
}
/***************************/
void EBAMRReactive::
setBCs()
{
  //Homogeneous Neumann boundary conditions for our smoothing operators
  NeumannPoissonEBBCFactory* neumbcSmooth = new NeumannPoissonEBBCFactory();
  neumbcSmooth->setValue(0.); 
  m_smoothEBBC = RefCountedPtr<BaseEBBCFactory>(neumbcSmooth);
  NeumannPoissonDomainBCFactory* neumdobcSmooth = new NeumannPoissonDomainBCFactory();
  neumdobcSmooth->setValue(0.);
  m_smoothDomBC = RefCountedPtr<BaseDomainBCFactory>(neumdobcSmooth);

  //no diffusion of species accross the embedded boundaries and domain boundaries
  NeumannConductivityEBBCFactory* neumbcspec =  new NeumannConductivityEBBCFactory();
  neumbcspec->setValue(0.);
  m_specEBBC = RefCountedPtr<BaseEBBCFactory>(neumbcspec);

  NeumannConductivityDomainBCFactory* neumdobcspec = new NeumannConductivityDomainBCFactory();
  neumdobcspec->setValue(0.);
  m_specDomBC = RefCountedPtr<BaseDomainBCFactory>(neumdobcspec);

  //no slip, no flow velocity accross the embedded boundaries and domain boundaries
  DirichletViscousTensorEBBCFactory* diribc = new DirichletViscousTensorEBBCFactory();
  diribc->setValue(0.);
  m_veloEBBC = RefCountedPtr<BaseEBBCFactory>(diribc);

  DirichletViscousTensorDomainBCFactory* diribcdom = new DirichletViscousTensorDomainBCFactory();
  diribcdom->setValue(0.);
  m_veloDomBC = RefCountedPtr<BaseDomainBCFactory>(diribcdom);

  // no heat flow accross the embedded boundaries and domain boundaries
  NeumannConductivityEBBCFactory* neumbctemp = new NeumannConductivityEBBCFactory();
  neumbctemp->setValue(0.); 
  m_tempEBBC = RefCountedPtr<BaseEBBCFactory>(neumbctemp);

  NeumannConductivityDomainBCFactory* neumdobc = new NeumannConductivityDomainBCFactory();
  neumdobc->setValue(0.);
  m_tempDomBC = RefCountedPtr<BaseDomainBCFactory>(neumdobc); 
} 
/***************************/
void EBAMRReactive::
getHalfState(LevelData<EBCellFAB>& a_stateInt)
{
  // interpolate state to n+1/2
  Real told = 0; Real tnew = 1; Real time = 0.5;
  EBArith::timeInterpolate(a_stateInt, m_stateOld, m_stateNew,
                           m_eblg.getDBL(), time, told, tnew);
}
/***************************/
void EBAMRReactive::
fillCoefficients(const LevelData<EBCellFAB>& a_state)
{
  CH_TIME("EBAMRReactive::fillCoefficients");
  
  EBCellFactory cellFact(m_eblg.getEBISL());
  EBFluxFactory fluxFact(m_eblg.getEBISL());
  int nghost = 4;  

  //finding the primitives at faces and using them to find coefficients
  LevelData<EBCellFAB> primCell    (m_eblg.getDBL(), m_nPrim, nghost*IntVect::Unit, cellFact);
  LevelData<EBCellFAB> presCell    (m_eblg.getDBL(),       1, nghost*IntVect::Unit, cellFact);
  LevelData<EBCellFAB> tempCell    (m_eblg.getDBL(),       1, nghost*IntVect::Unit, cellFact);
  LevelData<EBCellFAB> specDensCell(m_eblg.getDBL(), m_nSpec, nghost*IntVect::Unit, cellFact); 
  
  LevelData<EBFluxFAB> presFace    (m_eblg.getDBL(),       1, IntVect::Zero, fluxFact);
  LevelData<EBFluxFAB> tempFace    (m_eblg.getDBL(),       1, IntVect::Zero, fluxFact);  
  LevelData<EBFluxFAB> specDensFace(m_eblg.getDBL(), m_nSpec, IntVect::Zero, fluxFact);


  getPrimState(primCell, a_state);

  Interval consInterv(0,m_nComp-1);
  Interval dstInterv(0,0);
  Interval presInterv(QPRES,QPRES);
  Interval specDensInterv(CSPEC1,CSPEC1+m_nSpec-1);
  Interval tempInterv(QTEMP,QTEMP);
  Interval specInterv(0,m_nSpec-1);

  primCell.copyTo(presInterv, presCell, dstInterv);
  primCell.copyTo(tempInterv, tempCell, dstInterv);
  a_state.copyTo(specDensInterv, specDensCell, specInterv); 
  
  if(m_hasCoarser)
    {
      EBAMRReactive* coarPtr = getCoarserLevel();
      EBPWLFillPatch patcher(m_eblg.getDBL(), coarPtr->m_eblg.getDBL(), 
                             m_eblg.getEBISL(), coarPtr->m_eblg.getEBISL(),
                             coarPtr->m_eblg.getDomain(), m_ref_ratio, m_nComp, nghost);
            

      EBCellFactory coarFact(coarPtr->m_eblg.getEBISL());
      LevelData<EBCellFAB> coarPres(coarPtr->m_eblg.getDBL(), 1, nghost*IntVect::Unit, coarFact);
      LevelData<EBCellFAB> coarTemp(coarPtr->m_eblg.getDBL(), 1, nghost*IntVect::Unit, coarFact);
      LevelData<EBCellFAB> coarspecDens(coarPtr->m_eblg.getDBL(), m_nSpec, nghost*IntVect::Unit, coarFact);

      LevelData<EBCellFAB> coarPrim(coarPtr->m_eblg.getDBL(), m_nPrim, nghost*IntVect::Unit, coarFact);
      coarPtr->getPrimState(coarPrim, coarPtr->m_stateNew);
      coarPrim.copyTo(presInterv, coarPres, dstInterv);  
      coarPrim.copyTo(tempInterv, coarTemp, dstInterv);
      coarPtr->m_stateNew.copyTo(specDensInterv, coarspecDens, specInterv);

      patcher.interpolate(presCell, coarPres, coarPres, m_time, m_time, m_time, Interval(0,0));
      patcher.interpolate(tempCell, coarTemp, coarTemp, m_time, m_time, m_time, Interval(0,0));
      patcher.interpolate(specDensCell, coarspecDens, coarspecDens, m_time, m_time, m_time, Interval(0,m_nSpec-1));
    }

  EBLevelDataOps::averageCellToFaces(presFace, presCell, m_eblg.getDBL(), m_eblg.getEBISL(), m_eblg.getDomain(), 0);
  EBLevelDataOps::averageCellToFaces(tempFace, tempCell, m_eblg.getDBL(), m_eblg.getEBISL(), m_eblg.getDomain(), 0);

  for (int iSpec = 0; iSpec < m_nSpec; iSpec++)
   {
     EBLevelDataOps::averageCellToFaces(specDensFace, specDensCell, m_eblg.getDBL(), m_eblg.getEBISL(), m_eblg.getDomain(), iSpec);
   }
   
  setDiffusionCoefficients(presCell, tempCell, specDensCell, presFace, tempFace, specDensFace);
} 
/***************************/
void 
EBAMRReactive::
setDiffusionCoefficients(const LevelData<EBCellFAB>& a_presCell,
                         const LevelData<EBCellFAB>& a_tempCell,
                         const LevelData<EBCellFAB>& a_specDensCell,
                         const LevelData<EBFluxFAB>& a_presFace,
                         const LevelData<EBFluxFAB>& a_tempFace,
                         const LevelData<EBFluxFAB>& a_specDensFace)
{
  CH_TIME("EBAMRReactive::setDiffusionCoefficients");
  for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
   {
     Box cellBox = m_eblg.getDBL().get(dit());
     for (int idir = 0; idir < SpaceDim; idir++)
      {
        Box faceBox = surroundingNodes(cellBox, idir);
        // fill species by species even though this means repeating bunch of stuff in fort routines
        for (int iSpec = 0; iSpec < m_nSpec; iSpec++)
         {
           FORT_GETSPECMASSDIFFCOEFF(CHF_FRA1((*m_bco[iSpec])    [dit()][idir].getSingleValuedFAB(), 0),
                                     CHF_FRA((*m_rhsco[iSpec])   [dit()][idir].getSingleValuedFAB())   ,
                                     CHF_CONST_FRA1(a_presFace   [dit()][idir].getSingleValuedFAB(), 0),
                                     CHF_CONST_FRA1(a_tempFace   [dit()][idir].getSingleValuedFAB(), 0),
                                     CHF_CONST_FRA(a_specDensFace[dit()][idir].getSingleValuedFAB()),
                                     CHF_CONST_INT(iSpec),
                                     CHF_BOX(faceBox));          
          
           pout() << "bco for species" << iSpec << "for dir" << idir << endl;
           FabDataOps::getFabData((*m_bco[iSpec])[dit()][idir]);
         } // end species

        FORT_GETVISCOSITY(CHF_FRA1((*m_eta)           [dit()][idir].getSingleValuedFAB(), 0),
                          CHF_FRA1((*m_lambda)        [dit()][idir].getSingleValuedFAB(), 0),
                          CHF_CONST_FRA1(a_tempFace   [dit()][idir].getSingleValuedFAB(), 0),
                          CHF_CONST_FRA(a_specDensFace[dit][idir].getSingleValuedFAB())     ,
                          CHF_BOX(faceBox));

        pout() << "viscosity, eta for dir" << idir << endl;
        FabDataOps::getFabData((*m_eta)[dit()][idir]); 

        FORT_GETCONDUCTIVITY(CHF_FRA1((*m_kappa)   [dit()][idir].getSingleValuedFAB(), 0),
                             CHF_FRA1(a_tempFace   [dit()][idir].getSingleValuedFAB(), 0),
                             CHF_FRA(a_specDensFace[dit()][idir].getSingleValuedFAB())   ,
                             CHF_BOX(faceBox));
 
        pout() << "conductivity, kappa for dir" << idir << endl;
        FabDataOps::getFabData((*m_kappa)[dit()][idir]);

        Vector<FaceIndex> faces = m_eblg.getEBISL()[dit()].getEBGraph().getMultiValuedFaces(idir, cellBox);
        for (int iface = 0; iface < faces.size(); iface++)
         {
           Vector<Real> specDense(m_nSpec);
           for (int iSpec1 = 0; iSpec1 < m_nSpec; iSpec1++)
            {
              specDense[iSpec1] = a_specDensFace[dit()][idir](faces[iface], iSpec1);
            }

           for (int iSpec = 0; iSpec < m_nSpec; iSpec++)
            {
              Vector<Real> rhsco(m_nSpec); 
              
              FORT_POINTGETSPECMASSDIFFCOEFF(CHF_REAL((*m_bco[iSpec])   [dit()][idir](faces[iface], 0)),
                                             CHF_VR(rhsco),
                                             CHF_CONST_REAL(a_presFace  [dit()][idir](faces[iface], 0)),
                                             CHF_CONST_REAL(a_tempFace  [dit()][idir](faces[iface], 0)),
                                             CHF_CONST_VR(specDense),
                                             CHF_CONST_INT(iSpec));

              for (int iSpec2 = 0; iSpec2 < m_nSpec; iSpec2++)
               {
                  (*m_rhsco[iSpec]) [dit()][idir](faces[iface],iSpec2) = rhsco[iSpec2];
               } // end species 2 
            } // end species

           FORT_POINTGETVISCOSITY(CHF_REAL((*m_eta)          [dit()][idir](faces[iface], 0)),
                                  CHF_REAL((*m_lambda)       [dit()][idir](faces[iface], 0)),
                                  CHF_CONST_REAL(a_tempFace  [dit()][idir](faces[iface], 0)),
                                  CHF_CONST_VR(specDense));

           FORT_POINTGETCONDUCTIVITY(CHF_REAL((*m_kappa)        [dit()][idir](faces[iface], 0)),
                                     CHF_CONST_REAL(a_tempFace  [dit()][idir](faces[iface], 0)),
                                     CHF_CONST_VR(specDense));

         } // end faces loop
      } // end dir loop
    
     IntVectSet ivs = m_eblg.getEBISL()[dit()].getIrregIVS(cellBox); //irreg because we are also setting the irregular coeffs here
     for(VoFIterator vofit(ivs, m_eblg.getEBISL()[dit()].getEBGraph()); vofit.ok(); ++vofit)
      {
        Vector<Real> specDense(m_nSpec);
        for (int iSpec1 = 0; iSpec1 < m_nSpec; iSpec1++)
         {
           specDense[iSpec1] = a_specDensCell[dit()](vofit(), iSpec1);
         }

        for (int iSpec = 0; iSpec < m_nSpec; iSpec++)
         {
           Vector<Real> rhsco(m_nSpec);

           FORT_POINTGETSPECMASSDIFFCOEFF(CHF_REAL((*m_bcoIrreg[iSpec]) [dit()](vofit(), 0)),
                                          CHF_VR(rhsco)   ,
                                          CHF_CONST_REAL(a_presCell     [dit()](vofit(), 0)), 
                                          CHF_CONST_REAL(a_tempCell     [dit()](vofit(), 0)),
                                          CHF_CONST_VR(specDense)   ,
                                          CHF_CONST_INT(iSpec));

           for (int iSpec2 = 0; iSpec2 < m_nSpec; iSpec2++)
            {
              (*m_rhscoIrreg[iSpec])[dit()](vofit(), iSpec2) = rhsco[iSpec2];
            } // end species 2
         } // end species 

        FORT_POINTGETVISCOSITY(CHF_REAL((*m_etaIrreg)     [dit()](vofit(), 0)),
                               CHF_REAL((*m_lambdaIrreg)  [dit()](vofit(), 0)),
                               CHF_CONST_REAL(a_tempCell  [dit()](vofit(), 0)),
                               CHF_CONST_VR(specDense));

        FORT_POINTGETCONDUCTIVITY(CHF_REAL((*m_kappaIrreg)   [dit()](vofit(), 0)),
                                  CHF_CONST_REAL(a_tempCell  [dit()](vofit(), 0)),
                                  CHF_CONST_VR(specDense));
      } // end vofit  
   } // end dit
}
/***************************/
void EBAMRReactive::tagCellsInit(IntVectSet& a_tags)
{
  if (s_verbosity >= 3)
    {
      pout() << "EBAMRReactive::tagCellsInit for level " << m_level << endl;
    }

  // Since tags are calculated using only current time step data, use
  // the same tagging function for initialization and for regridding.
  tagCells(a_tags);
}
/***************************/
void EBAMRReactive::tagCells(IntVectSet& a_tags)
{
  CH_TIME("EBAMRReactive::tagCells");
  if (s_verbosity >= 3)
    {
      pout() << "EBAMRReactive::tagCells for level " << m_level << endl;
    }

  // Create tags based on undivided gradient of density
  IntVectSet localTags;
  if (!m_tagAll)
    {
      ParmParse pp;
      int itagAllIrreg;
      int queryVal = pp.query("tag_all_irreg", itagAllIrreg);
      bool tagAllIrregCells = ((queryVal != 0) && (itagAllIrreg == 1));

      Real xmaxTag, ymaxTag, xminTag, yminTag;
      vector<Real> problo(SpaceDim);
      bool tagTheLot = false;

      // compute gradients and all that to compute tags
      // If there is a coarser level interpolate undefined ghost cells
      //only interpolate the density
      int densityIndex = m_ebPatchReactive->densityIndexC();
      Interval intervDensity(densityIndex, densityIndex);
      EBCellFactory factory(m_ebisl);
      int nCons = m_ebPatchReactive->numConserved();
      LevelData<EBCellFAB> consTemp(m_grids, nCons, IntVect::Unit, factory);
      Interval consInterv(0, nCons-1);
      m_stateNew.copyTo(consInterv, consTemp, consInterv);
      if (m_hasCoarser)
        {
          const EBAMRReactive* amrReacCoarserPtr = getCoarserLevel();
          int refRatCrse = amrReacCoarserPtr->refRatio();
          int nghost = 1;
          EBPWLFillPatch patcher(m_grids,
                                 amrReacCoarserPtr->m_grids,
                                 m_ebisl,
                                 amrReacCoarserPtr->m_ebisl,
                                 amrReacCoarserPtr->m_domainBox,
                                 refRatCrse, m_nComp, nghost);
 
          Real coarTimeOld = 0.0;
          Real coarTimeNew = 1.0;
          Real fineTime    = 0.0;
          patcher.interpolate(consTemp,
                              amrReacCoarserPtr->m_stateOld,
                              amrReacCoarserPtr->m_stateNew,
                              coarTimeOld,
                              coarTimeNew,
                              fineTime,
                              intervDensity);
        }
      consTemp.exchange(intervDensity);

      // Compute undivided gradient
      for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
        {
          const Box& b = m_grids.get(dit());
          const EBISBox& ebisBox = m_ebisl[dit()];
          EBCellFAB gradFab(ebisBox, b, SpaceDim);
          const EBCellFAB& stateFab = consTemp[dit()];
          BaseFab<Real>& regGradFab = gradFab.getSingleValuedFAB();
          const BaseFab<Real>& regStateFab = stateFab.getSingleValuedFAB();

          for (int idir = 0; idir < SpaceDim; ++idir)
            {
              const Box bCenter = b & grow(m_problem_domain,-BASISV(idir));
              const Box bLo     = b & adjCellLo(bCenter,idir);
              const int hasLo = ! bLo.isEmpty();
              const Box bHi     = b & adjCellHi(bCenter,idir);
              const int hasHi = ! bHi.isEmpty();

              FORT_GETRELATIVEGRAD(CHF_FRA1(regGradFab,idir),
                                   CHF_CONST_FRA1(regStateFab,0),
                                   CHF_CONST_INT(idir),
                                   CHF_BOX(bLo),
                                   CHF_CONST_INT(hasLo),
                                   CHF_BOX(bHi),
                                   CHF_CONST_INT(hasHi),
                                   CHF_BOX(bCenter));
 
             //do one-sided diffs where necessary at irregular cells.
             IntVectSet ivsIrreg = ebisBox.getIrregIVS(b);
             for (VoFIterator vofit(ivsIrreg, ebisBox.getEBGraph());
                  vofit.ok(); ++vofit)
               {
                 const VolIndex& vof = vofit();
                 const IntVect&  iv = vof.gridIndex();
                 //one-sided diffs on domain bndry
                 bool onLeftDomain = iv[idir] == m_domainBox.smallEnd(idir);
                 bool onRighDomain = iv[idir] == m_domainBox.bigEnd(idir);
                 bool hasFacesLeft = (ebisBox.numFaces(vof, idir, Side::Lo) > 0) && !onLeftDomain;
                 bool hasFacesRigh = (ebisBox.numFaces(vof, idir, Side::Hi) > 0) && !onRighDomain;

                 Real valCent = stateFab(vof, densityIndex);
                 Real dpl=0;
                 Real dpr=0;
                 Real dpc=0;
                 //compute one-sided diffs where you have them
                 if (hasFacesLeft)
                   {
                     Vector<FaceIndex> facesLeft = ebisBox.getFaces(vof, idir, Side::Lo);
                     Real valLeft = 0.0;
                     for (int iface = 0; iface <facesLeft.size(); iface++)
                       {
                         VolIndex vofLeft = facesLeft[iface].getVoF(Side::Lo);
                         valLeft += stateFab(vofLeft, densityIndex);
                       }
                     valLeft /= Real(facesLeft.size());
                     dpl = valCent - valLeft;
                   } 
                 if (hasFacesRigh)
                   {
                     Vector<FaceIndex> facesRigh = ebisBox.getFaces(vof, idir, Side::Hi);
                     Real valRigh = 0.0;
                     for (int iface = 0; iface <facesRigh.size(); iface++)
                      {
                        VolIndex vofRigh = facesRigh[iface].getVoF(Side::Hi);
                        valRigh += stateFab(vofRigh, densityIndex);
                      }
                    valRigh /= Real(facesRigh.size());
                    dpr = valRigh - valCent;
                   }
                 if (hasFacesLeft && hasFacesRigh) 
                    {
                      dpc = 0.5*(dpl+dpr);
                    }
                 else if (!hasFacesLeft && !hasFacesRigh)
                    {
                      dpc = 0.0;
                    }
                 else if (hasFacesLeft && !hasFacesRigh)
                   {
                     dpc = dpl;
                   }
                 else if (hasFacesRigh && !hasFacesLeft) 
                   {
                     dpc = dpr;
                   }

                 gradFab(vof,idir) = dpc/valCent;
               }
            }

          EBCellFAB gradMagFab(ebisBox, b, 1);
          BaseFab<Real>& regGradMagFab = gradMagFab.getSingleValuedFAB();
          FORT_MAGNITUDE(CHF_FRA1(regGradMagFab,0),
                         CHF_CONST_FRA(regGradFab),
                         CHF_BOX(b));

          //pointwise op so just have to iterate over multivalued cells 
          IntVectSet ivsMulti = ebisBox.getMultiCells(b);
          for (VoFIterator vofit(ivsMulti, ebisBox.getEBGraph());
               vofit.ok(); ++vofit)
            {
              const VolIndex& vof = vofit();
              Real mag = 0.0;
              for (int idir = 0; idir < SpaceDim; idir++)
                {
                  Real graddir = gradFab(vof, idir);
                  mag += graddir*graddir;
                }
              mag = sqrt(mag);
              gradMagFab(vof, 0) = mag;
            }

          // Tag where gradient exceeds threshold

          IntVectSet ivsTot(b);
          for (VoFIterator vofit(ivsTot, ebisBox.getEBGraph());
              vofit.ok(); ++vofit)
            {
              const VolIndex& vof = vofit();
              const IntVect& iv = vof.gridIndex();
              Real xval = 0;
              Real yval = 0;
              if (((xval < xmaxTag) && (yval < ymaxTag)))
                {
                  if (gradMagFab(vof, 0) >= m_refineThresh)
                    {
                      localTags |= iv;
                    }
                }
            }
        
          //refine all irregular cells
          //this is probably not ideal but is used for replicate
          //scaling tests.
          if (tagAllIrregCells)
            {
              IntVectSet irregIVS = ebisBox.getIrregIVS(b);
              localTags |= irregIVS;
            }
        }

      localTags.grow(m_tagBufferSize);

      // Need to do this in two steps unless a IntVectSet::operator &=
      //  (ProblemDomain) operator is defined
      Box localTagsBox = localTags.minBox();
      localTagsBox &= m_problem_domain; 
      localTags &= localTagsBox;
    }//if !tagall
  else
    {
      localTags.makeEmpty();
      for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
        {
          const Box& b = m_grids.get(dit());
          const EBISBox& ebisBox = m_ebisl[dit()];
          IntVectSet ivsTot(b);
          for (VoFIterator vofit(ivsTot, ebisBox.getEBGraph());
              vofit.ok(); ++vofit)
            {
              localTags |= vofit().gridIndex();
            }
        }
    }

  a_tags = localTags;
}
/***************************/
void EBAMRReactive::postInitialize()
{
}
/***************************/
LevelData<EBCellFAB>&
EBAMRReactive::getStateNew()
{
  return m_stateNew;
}
/***************************/
LevelData<EBCellFAB>&
EBAMRReactive::getStateOld()
{
  return m_stateOld;
}
/***************************/
void EBAMRReactive::fillGhosts()
{
  // fill all the ghost cells distributed accross all processors
  Vector<EBAMRReactive*>       hierarchy;
  Vector<DisjointBoxLayout>    grids;
  Vector<EBISLayout>           ebisl;
  Vector<EBLevelGrid>          eblg;
  Vector<int>                  refRat;
  ProblemDomain                lev0Dom;
  Real                         lev0Dx;
  Vector<ProblemDomain>        domains;
  Vector<Real>                 dxs;
  getHierarchyAndGrids(hierarchy, grids, ebisl, eblg, refRat, lev0Dom, lev0Dx, domains, dxs);

  int nlevels = hierarchy.size();

  Vector<LevelData<EBCellFAB>* > data(nlevels);
  for(int i=0; i< nlevels; i++) data[i] = &(hierarchy[i]->m_stateNew);
  EBAMRDataOps::pwlFillPatchAll(data,eblg,refRat);
  EBAMRDataOps::quadCFInterpAll(data,eblg,refRat);
  EBAMRDataOps::averageDown(data,eblg,refRat);

  // debug:
   int a = EBAMRDataOps::countVoF(grids,ebisl,domains);
  // end debug

}
/***************************/
void
EBAMRReactive::
getHierarchyAndGrids(Vector<EBAMRReactive*>&      a_hierarchy,
                     Vector<DisjointBoxLayout>&   a_grids,
                     Vector<EBISLayout>&          a_ebisl,
                     Vector<EBLevelGrid>&         a_eblg,
                     Vector<int>&                 a_refRat,
                     ProblemDomain&               a_lev0Dom,
                     Real&                        a_lev0Dx,
                     Vector<ProblemDomain>&       a_domains,
                     Vector<Real>&                a_dxs)
{
  // retrieve an array of all the AMRLevel objects in the entire hierarchy
  Vector<AMRLevel*> hierarchy = AMRLevel::getAMRLevelHierarchy();
  int nlevels = hierarchy.size();

  a_hierarchy.resize(nlevels);
  a_refRat.resize(nlevels);
  a_grids.resize(nlevels);
  a_ebisl.resize(nlevels);
  a_eblg.resize(nlevels);
  a_domains.resize(nlevels);
  a_dxs.resize(nlevels);

  EBAMRReactive* coarsestLevel = (EBAMRReactive*)(hierarchy[0]);
  a_lev0Dx     = coarsestLevel->m_dx[0];
  a_lev0Dom    = coarsestLevel->m_problem_domain;

  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();

  for (int ilev = 0; ilev < nlevels; ilev++)
    {
      EBAMRReactive* adLevel = (EBAMRReactive*)(hierarchy[ilev]);

      a_hierarchy[ilev] = adLevel;
      a_grids[ilev] = adLevel->m_grids;
      a_refRat[ilev] = adLevel->m_ref_ratio;
      const ProblemDomain& domainLev = adLevel->problemDomain();
      ebisPtr->fillEBISLayout(a_ebisl[ilev], a_grids[ilev], domainLev, m_nGhost);
      a_eblg[ilev] = EBLevelGrid(a_grids[ilev], a_ebisl[ilev], domainLev);
      a_domains[ilev] = domainLev;
      a_dxs[ilev] = adLevel->m_dx[0];
    }
}
/***************************/
void EBAMRReactive::
postInitialGrid(const bool a_restart)
{
}

/***************************/
void EBAMRReactive::syncWithFineLevel()
{
  CH_TIME("EBAMRR::syncWithFineLevel");
  if (s_verbosity >= 3)
    {
      pout() << " in EBAMRReactive syncWithFineLevel for level " << m_level << endl;
    }
  //stuff that needs to be setup from the finer
  //level.  A bunch of objects depend upon the layouts
  //from both levels and the finer level changes more
  //often from regrid so this needs to be called from the finer
  //level
  CH_assert(m_hasFiner);
  if (m_hasFiner)
    {
      EBAMRReactive* finePtr = getFinerLevel();
      int nRefFine = refRatio();
      const EBLevelGrid& finer_eblg = finePtr->m_eblg;
      const DisjointBoxLayout& finer_dbl = finer_eblg.getDBL();
      const EBISLayout& finer_ebisl = finer_eblg.getEBISL();

      // maintain flux registers
      m_ebFluxRegister.define(finer_dbl,
                                m_eblg.getDBL(),
                                finer_ebisl,
                                m_eblg.getEBISL(),
                                m_eblg.getDomain().domainBox(),
                                nRefFine,
                                m_nComp, Chombo_EBIS::instance(), s_noEBCF);

      m_veloFluxRegister.define(finer_dbl,
                                m_eblg.getDBL(),
                                finer_ebisl,
                                m_eblg.getEBISL(),
                                m_eblg.getDomain().domainBox(),
                                nRefFine,
                                SpaceDim, Chombo_EBIS::instance(), s_noEBCF);

      m_tempFluxRegister.define(finer_dbl,
                                m_eblg.getDBL(),
                                finer_ebisl,
                                m_eblg.getEBISL(),
                                m_eblg.getDomain().domainBox(),
                                nRefFine,
                                1, Chombo_EBIS::instance(), s_noEBCF);

      // set all registers to zero
      m_ebFluxRegister.setToZero();
      m_veloFluxRegister.setToZero();
      m_tempFluxRegister.setToZero();       
 
      m_specFluxRegister.define(finer_dbl,
                                m_eblg.getDBL(),
                                finer_ebisl,
                                m_eblg.getEBISL(),
                                m_eblg.getDomain().domainBox(),
                                nRefFine,
                                m_nSpec, Chombo_EBIS::instance(), s_noEBCF);

      m_specFluxRegister.setToZero();

      //define fine to coarse redistribution object
      //for now set to volume weighting
      if (!s_noEBCF)
       {
         m_ebCoarToFineRedist.define(finer_dbl, m_eblg.getDBL(),
                                     finer_ebisl,  m_eblg.getEBISL(),
                                     m_eblg.getDomain().domainBox(), nRefFine , m_nComp, 
                                     1, Chombo_EBIS::instance());

         //define coarse to coarse redistribution object
         m_ebCoarToCoarRedist.define(finer_dbl, m_eblg.getDBL(),
                                     finer_ebisl,  m_eblg.getEBISL(),
                                     m_eblg.getDomain().domainBox(), nRefFine , m_nComp);
       }

      // set all the registers to zero
      if (!s_noEBCF)
       {
         m_ebCoarToFineRedist.setToZero();
         m_ebCoarToCoarRedist.setToZero();
       }
    }
}
/***************************/
Real EBAMRReactive::computeInitialDt()
{
  Real maxwavespeed =  m_ebLevelReactive.getMaxWaveSpeed(m_stateNew);
  CH_assert(Abs(m_dx[0]-m_dx[1])<1.e-9);
  Real newDT = m_initial_dt_multiplier * m_dx[0] /maxwavespeed;

  return newDT;
}
/***************************/
void EBAMRReactive::regrid(const Vector<Box>& a_new_grids)
{
  CH_TIME("EBAMRReactive::regrid");
  if (s_verbosity >= 3)
    {
      pout() << " in EBAMRReactive regrid for level " << m_level << endl;
    }
  Vector<Box>& newGrids =   (Vector<Box>&) a_new_grids;
  {
    CH_TIME("mortonordering");
    mortonOrdering(newGrids);
  }

  const EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();
  // save data for later copy
  //not using m_ebisl because it gets wiped later
  //the ebisl has to know about the fact that we really have
  //four ghost cells.
  int nGhostEBISL = 6;
  Interval interv(0,m_nComp-1);
  LevelData<EBCellFAB> stateSaved;
  IntVect ivGhost = m_nGhost*IntVect::Unit;
  {
    CH_TIME("defines and copies and ebisls");
    EBISLayout ebislOld;
    ebisPtr->fillEBISLayout(ebislOld, m_grids, m_domainBox, nGhostEBISL);
    EBCellFactory factoryOld(ebislOld);
    stateSaved.define(m_grids, m_nComp, ivGhost, factoryOld);
    m_stateNew.copyTo(interv, stateSaved, interv);
  }
  //create grids and ebis layouts
  //the ebisl has to know about the fact that we really have
  //four ghost cells.
  m_level_grids = a_new_grids;
  Vector<int> proc_map;

  if (s_isLoadBalanceSet)
    {
      s_loadBalance(proc_map,a_new_grids, m_domainBox, false);
    }
  else
    {
      LoadBalance(proc_map,a_new_grids);
    }

  m_grids= DisjointBoxLayout(a_new_grids,proc_map);
  m_eblg.define(m_grids, m_problem_domain, 6, ebisPtr);

  ebisPtr->fillEBISLayout(m_ebisl, m_grids, m_domainBox, nGhostEBISL);

  EBCellFactory factoryNew(m_ebisl);
  // reshape state with new grids
  m_stateNew.define(m_grids,m_nComp,ivGhost, factoryNew);
  m_stateOld.define(m_grids,m_nComp,ivGhost, factoryNew);

  // set up data structures
  levelSetup();

  // interpolate to coarser level
  if (m_hasCoarser)
    {
      EBAMRReactive* coarPtr = getCoarserLevel();
      m_ebFineInterp.interpolate(m_stateNew,
                                 coarPtr->m_stateNew,
                                 interv);
    }

  // copy from old state
  stateSaved.copyTo(interv,m_stateNew, interv);
}
/***************************/
Real EBAMRReactive::advance()
{
  // advance the conservative state by one time step
  CH_TIME("EBAMRReactive::advance");
  EBPatchReactive::s_whichLev = m_level;
  EBViscousTensorOp::s_step = AMR::s_step;
  //this is so I can output a rhs
  EBViscousTensorOp::s_whichLev = m_level;
  m_dtOld = m_dt;

  dumpDebug(string("going into advance"));

  if (addDiffusion() && (m_level == 0)) defineSolvers();

  if (s_verbosity >= 3)
    {
      pout() << " in EBAMRReactive advance for level =" << m_level << ", with dt = " << m_dt << endl;
    }
  {
    CH_TIME("copy new to old");
    m_stateNew.copyTo(m_stateNew.interval(),
                      m_stateOld,
                      m_stateOld.interval());
  }

  //set up arguments to step
  //undefined lfr in case we need it
  EBFluxRegister lfr;
  //undefined leveldata in case we need it
  const LevelData<EBCellFAB> ld;

  //set arguments to dummy arguments and
  //then fix if they are available
  EBFluxRegister* coarFR = &lfr;
  EBFluxRegister* fineFR = &lfr;
  const LevelData<EBCellFAB>* coarDataOld = &ld;
  const LevelData<EBCellFAB>* coarDataNew = &ld;

  Real told = 0.0;
  Real tnew = 0.0;
  if (m_hasCoarser)
    {
      EBAMRReactive* coarPtr = getCoarserLevel();
      coarFR = &coarPtr->m_ebFluxRegister;
      coarDataOld = &coarPtr->m_stateOld;
      coarDataNew = &coarPtr->m_stateNew;
      tnew = coarPtr->m_time;
      told = tnew - coarPtr->m_dt;
      //time should never be greater than the newest coarse
      //time.  time might be very slightly smaller than
      //told because of the above subtraction.
      Real eps = 1.0e-10;
      if ((m_time > tnew) || (m_time < (told - eps)))
        {
          MayDay::Error("out of bounds time input to AMRReactive");
        }
      //correct for said floating-point nastiness
      m_time = Max(m_time, told);
    }
  if (m_hasFiner)
    {
      //recall that my flux register goes between my
      //level and the next finer level
      fineFR = &m_ebFluxRegister;
    }

#ifndef NDEBUG
  if (!m_hasCoarser && (s_verbosity > 1))
    {
      Real summass;
      int densityIndex = m_ebPatchReactive->densityIndex();
      sumConserved(summass,  densityIndex);

      pout() << "sum mass = " << summass << endl;
    }
#endif

  int nghost = 4;
  EBCellFactory fact(m_eblg.getEBISL());
  LevelData<EBCellFAB> source(m_eblg.getDBL(), m_nComp, nghost*IntVect::Unit, fact);
  LevelData<EBCellFAB> divergeF(m_eblg.getDBL(), m_nComp, nghost*IntVect::Unit, fact);

  hyperbolicSource(source);

  pout() << "getting diverence of flux" << endl;

  //PDivU and rhouedele are not multiplied by the volume fraction (that happens in the integrator.
  //finding the divergence of Flux
  bool onlyDiffusion = false;
  ParmParse pp;

  if (pp.contains("do_only_diffusion"))
      {
        pp.get("do_only_diffusion", onlyDiffusion);
      }
  
  
  if (onlyDiffusion)
    {
      EBLevelDataOps::setToZero(divergeF);
    }
  else
   {
     m_ebLevelReactive.divergeF(divergeF,
                                m_massDiff,
                                *fineFR,
                                *coarFR,
                                m_stateOld,
                                source,
                                *coarDataOld,
                                *coarDataNew,
                                m_time,
                                told,
                                tnew,
                                m_dt); 
    }

  coarseFineIncrement();

  if (!(addDiffusion()))
    {
      pout() << "explicit Advance" << endl;
      explicitAdvance(divergeF);
    }
  else
    {
      CH_TIME("diffusion dance");
      //step 1) get u*
      pout() << "step 1: getting Ustar (U^n + explicit contribution)" << endl;

      //U* = Un - dt*LE(U^n)
      LevelData<EBCellFAB>       UStar(m_eblg.getDBL(), m_nComp, nghost*IntVect::Unit, fact);
      getUStar(UStar, m_stateOld, divergeF);

      pout() << "step 1.5: doing explicit redistribution into ustar of hyperbolic mass difference" <<endl;
      //redistribute stuff in hyperbolic registers.
      hyperbolicRedistribution(UStar);

      //step 2) update species mass densities and add them to UStar
      pout() << "step 2: updating species multicomponent diffusion" << endl;
      addMCSpecDiff(UStar); 

      //step 3) update momentum with viscous terms
      pout() << "step 3: updating momentum with viscocity" << endl;
      addViscosity(UStar);

      //step 4) update energy with conductive terms
      pout() << "step 4: updating energy with conductivity " << endl;
      addConductivity(UStar);
 
      pout() << "step 5: putting state into m_statenew and flooring" << endl;
      finalAdvance(UStar);
    }

  Real new_dt;
  {
    CH_TIME("time_step_calc");
    //deal with time and time step
    Real maxWaveSpeed = m_ebLevelReactive.getMaxWaveSpeed(m_stateOld);
    new_dt = m_cfl*m_dx[0]/maxWaveSpeed;
    pout() << "max wave speed = " << maxWaveSpeed << ", dx = " << m_dx[0] <<", new_dt = "  << new_dt << endl;
    m_time += m_dt;
  }

  //save stable timestep to save computational effort
  m_dtNew = new_dt;

  if (addReactionRates())
   {
     // step 6:
     pout() << "step 6: integrating reaction source terms" << endl;
     // computations are not done on ghost cells. Ghost cells are later filled in posTimeStep
     {
       m_ebLevelReactive.integrateReactiveSource(m_stateNew,m_domainBox,m_time,new_dt);
     }
   } 
  return new_dt; 
}
/***************************/ 
void EBAMRReactive::
getUStar(LevelData<EBCellFAB>      & a_UStar,
         const LevelData<EBCellFAB>& a_UN,
         const LevelData<EBCellFAB>& a_divergeF)
{
  CH_TIME("EBAMRReactive::getUStar");

  EBCellFactory fact(m_eblg.getEBISL());
  //advance everything explicitly
  for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      EBCellFAB  dtLcU(m_eblg.getEBISL()[dit()], m_eblg.getDBL().get(dit()), m_nComp);
      dtLcU.setVal(0.);
      dtLcU += a_divergeF[dit()];
      dtLcU *= m_dt;

      a_UStar[dit()].setVal(0.);
      a_UStar[dit()] += a_UN[dit()];
      a_UStar[dit()] -= dtLcU;
    }

  m_ebLevelReactive.floorConserved(a_UStar, m_time, m_dt);
}
/***************************/
void EBAMRReactive::
explicitAdvance(const LevelData<EBCellFAB>& a_divergeF)
{
  CH_TIME("EBAMRCNS::explicitAdvance");
  CH_assert(!(addDiffusion()));
  //advance everything explicitly
  for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      EBCellFAB dtDivergeF(m_eblg.getEBISL()[dit()], m_eblg.getDBL().get(dit()), a_divergeF.nComp());
      dtDivergeF.setVal(0.);
      dtDivergeF += a_divergeF[dit()];
      dtDivergeF *= m_dt;
      m_stateNew[dit()] -= dtDivergeF;
    }
  hyperbolicRedistribution(m_stateNew);

  m_ebLevelReactive.floorConserved(m_stateNew, m_time, m_dt);
}
/***************************/
void EBAMRReactive::
hyperbolicRedistribution(LevelData<EBCellFAB>& a_state)
{
  CH_TIME("EBAMRReactive::hyperbolic redistribution");
  //explicit redistribution here because these are the hyperbolic terms
  m_ebLevelRedist.setToZero();
  if(m_useMassRedist)
    {
      //if use mass weighting, need to
      //fix weights of redistribution object
      int densityIndex = m_ebPatchReactive->densityIndex();
      a_state.exchange(Interval(0, m_nComp-1));
      m_ebLevelRedist.resetWeights(a_state, densityIndex);
    }

  Interval consInterv(0, m_nComp-1);
  for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      m_ebLevelRedist.increment(m_massDiff[dit()], dit(), consInterv);
    }
  //
  m_ebLevelRedist.redistribute(a_state, consInterv);

  m_ebLevelRedist.setToZero();
}
/***************************/
void EBAMRReactive::
hyperbolicSource(LevelData<EBCellFAB>&  a_source)
{
  CH_TIME("EBAMRReactive::hyperbolicSource");

  EBPatchReactive::useConservativeSource(true);
  //this is important because sometimes there is no diffusion
  EBLevelDataOps::setVal(a_source, 0.0);
  if(addDiffusion())
   {
     EBCellFactory fact(m_eblg.getEBISL());
     int nghost = 4;
     LevelData<EBCellFAB> specMassSrc(m_eblg.getDBL(), m_nSpec, nghost*IntVect::Unit, fact);
     LevelData<EBCellFAB> momenSrc(m_eblg.getDBL(),   SpaceDim, nghost*IntVect::Unit, fact);
     LevelData<EBCellFAB> energSrc(m_eblg.getDBL(),          1, nghost*IntVect::Unit, fact);
     explicitHyperbolicSource(specMassSrc, momenSrc, energSrc, m_stateOld, true);

     Interval specSrcInterv(CSPEC1, CSPEC1+m_nSpec-1);
     Interval momeSrcInterv(CMOMX, CMOMX+SpaceDim-1);
     Interval enerSrcInterv(CENG, CENG);

     Interval specDstInterv(0, m_nSpec-1);
     Interval momeDstInterv(0, SpaceDim-1);
     Interval enerDstInterv(0, 0);

     for (DataIterator dit = a_source.dataIterator(); dit.ok(); ++dit)
      {
        const Box& region = m_eblg.getDBL().get(dit());
        a_source[dit()].copy(region, specSrcInterv, region, specMassSrc[dit()], specDstInterv);
        a_source[dit()].copy(region, momeSrcInterv, region, momenSrc[dit()], momeDstInterv);
        a_source[dit()].copy(region, enerSrcInterv, region, energSrc[dit()], enerDstInterv);

      } // end dit

     a_source.exchange();
   } // end addDiffusion()
}
/***************************/
void EBAMRReactive::
explicitHyperbolicSource(LevelData<EBCellFAB>&  a_specMassSrc,
                         LevelData<EBCellFAB>&  a_momenSrc,
                         LevelData<EBCellFAB>&  a_energSrc,
                         const LevelData<EBCellFAB>& a_state,
                         bool a_doNormalization)
{
  CH_TIME("EBAMRReactive::explicitHyperbolicSource");
  pout()<< "computing explicit hyperbolic source" << endl;
  EBCellFactory fact(m_eblg.getEBISL());
  int nghost = 4;
  LevelData<EBCellFAB> primState(m_eblg.getDBL(), m_nPrim, nghost*IntVect::Unit, fact);
  LevelData<EBCellFAB> specMF(m_eblg.getDBL(), m_nSpec, nghost*IntVect::Unit, fact);
  LevelData<EBCellFAB> velocity(m_eblg.getDBL(), SpaceDim, nghost*IntVect::Unit, fact);
  LevelData<EBCellFAB> temperature(m_eblg.getDBL(), 1, nghost*IntVect::Unit, fact);

  LevelData<EBCellFAB>* specMFCoar = NULL;
  LevelData<EBCellFAB>* velocityCoar = NULL;
  LevelData<EBCellFAB>* temperatureCoar = NULL;

  getPrimState(primState, a_state);

  Interval specSrcInterv(QSPEC1, QSPEC1+m_nSpec-1);
  Interval veloSrcInterv(QVELX, QVELX+SpaceDim-1);
  Interval tempSrcInterv(QTEMP, QTEMP);

  Interval specDstInterv(0, m_nSpec-1);
  Interval veloDstInterv(0, SpaceDim-1);
  Interval tempDstInterv(0, 0);

  primState.copyTo(specSrcInterv, specMF, specDstInterv);
  primState.copyTo(veloSrcInterv, velocity, veloDstInterv);
  primState.copyTo(tempSrcInterv, temperature, tempDstInterv);

  getCoarserPrimitives(specMFCoar, velocityCoar, temperatureCoar);

  //get the source terms to be added explicitly to hyperbolic advance
  kappaSpecMassDiffSrc(a_specMassSrc, specMF, specMFCoar, a_state);
  kappaMomentumSrc(a_momenSrc, velocity, velocityCoar, a_state);
  kappaEnergySrc(a_energSrc, velocity, velocityCoar, temperature, temperatureCoar, specMF, a_state);

  //the getcoarser stuff did news
  if (m_hasCoarser)
   {
     delete specMFCoar;
     delete velocityCoar;
     delete temperatureCoar;
   }

  if (a_doNormalization)
   {
     //finally all these things are multiplied by kappa so we need to make 
     //them normalized 
     KappaSquareNormal normalizinator(m_eblg);
     normalizinator(a_specMassSrc);
     normalizinator(a_momenSrc);
     normalizinator(a_energSrc);
   }
}
/***************************/
void EBAMRReactive::
addMCSpecDiff(LevelData<EBCellFAB>&       a_UStar)
{
  CH_TIME("EBAMRReactive::addSpecMFDiff");

  int nghost = 4;
  EBCellFactory fact(m_eblg.getEBISL());
  LevelData<EBCellFAB> primOld(m_eblg.getDBL(), m_nPrim, nghost*IntVect::Unit, fact);
  LevelData<EBCellFAB> MFnew(m_eblg.getDBL(), m_nSpec, nghost*IntVect::Unit, fact);
  LevelData<EBCellFAB> MFold(m_eblg.getDBL(), m_nSpec, nghost*IntVect::Unit, fact);
  LevelData<EBCellFAB> sumMF(m_eblg.getDBL(), 1, nghost*IntVect::Unit, fact);

  EBLevelDataOps::setToZero(sumMF);

  getPrimState(primOld, a_UStar);

  Interval specSrcInterv(QSPEC1, QSPEC1+m_nSpec-1);
  Interval specDstInterv(0, m_nSpec-1);

  primOld.copyTo(specSrcInterv, MFold, specDstInterv);

  LevelData<EBCellFAB>* MFCoarPtr   = NULL;
  LevelData<EBCellFAB>* veloCoarPtr = NULL;
  LevelData<EBCellFAB>* tempCoarPtr = NULL;

  getCoarserPrimitives(MFCoarPtr, veloCoarPtr, tempCoarPtr);
 
  Interval srcInt(CRHO,CRHO);
  Interval dstInt(0, 0);
  for (int iSpec = 0; iSpec < m_nSpec; iSpec++)
   {
     for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
      {
        const Box& region = m_eblg.getDBL().get(dit());
        (*m_acoDiff[iSpec])[dit()].copy(region, dstInt, region, a_UStar[dit()], srcInt);
      }
    }  
  
  Real tCoarOld = 0.0; 
  Real tCoarNew = 0.0;

  //species by species
  for (int iSpec = 0; iSpec < m_nSpec; iSpec++)
   {
     LevelData<EBCellFAB> specMFold(m_eblg.getDBL(), 1, nghost*IntVect::Unit, fact);
     LevelData<EBCellFAB> specMFnew(m_eblg.getDBL(), 1, nghost*IntVect::Unit, fact);
     LevelData<EBCellFAB> rhs(m_eblg.getDBL(), 1, nghost*IntVect::Unit, fact);
     EBLevelDataOps::setToZero(rhs);
  
     Interval srcInt1(iSpec, iSpec);
     Interval srcInt2(QSPEC1+iSpec, QSPEC1+iSpec);
     Interval dstInt1(0, 0);
     MFold.copyTo(srcInt1, specMFold, dstInt1);

     EBFluxRegister*       coarMFFRPtr  = NULL;
     EBFluxRegister*       fineMFFRPtr  = NULL;
     LevelData<EBCellFAB>* specMFCoarOldPtr = NULL;
     LevelData<EBCellFAB>* specMFCoarNewPtr = NULL;

     if (m_hasCoarser)
      {

        // get coarser data if we have it
        EBAMRReactive* coarPtr = getCoarserLevel(); 
        tCoarNew = coarPtr->m_time; 
        tCoarOld = tCoarNew - coarPtr->m_dt;
        const EBLevelGrid& ceblg = coarPtr->m_eblg;
        EBCellFactory cfact(ceblg.getEBISL());
        LevelData<EBCellFAB> primCoarOld(ceblg.getDBL(), m_nPrim, nghost*IntVect::Unit, cfact);
        LevelData<EBCellFAB> primCoarNew(ceblg.getDBL(), m_nPrim, nghost*IntVect::Unit, cfact);
        coarPtr->getPrimState(primCoarOld, coarPtr->m_stateOld);
        coarPtr->getPrimState(primCoarNew, coarPtr->m_stateNew);

        coarMFFRPtr = &coarPtr->m_specFluxRegister;
        
        specMFCoarOldPtr = new LevelData<EBCellFAB>(ceblg.getDBL(), 1, nghost*IntVect::Unit, cfact);
        specMFCoarNewPtr = new LevelData<EBCellFAB>(ceblg.getDBL(), 1, nghost*IntVect::Unit, cfact);
        
        primCoarOld.copyTo(srcInt2, *specMFCoarOldPtr, dstInt1);
        primCoarNew.copyTo(srcInt2, *specMFCoarNewPtr, dstInt1);
      }

     if (m_hasFiner)
      {
        fineMFFRPtr = &m_specFluxRegister; 
      }      

      // set rhs:
      Real alpha = 0; Real beta = 1;
      getMCDiffTerm(rhs, MFold, MFCoarPtr, alpha, beta, iSpec, false, true);       

     //Backward Euler Integrator:
     s_diffuseLevBE[iSpec]->updateSoln(specMFnew, specMFold, rhs, fineMFFRPtr, coarMFFRPtr, 
                                       specMFCoarOldPtr, specMFCoarNewPtr, m_time, tCoarOld, tCoarNew, m_dt,
                                       m_level, true, true, iSpec);

      specMFnew.copyTo(dstInt1, MFnew, srcInt1);

     for(DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
      {
        sumMF[dit()] += specMFnew[dit()];  
      }
 
     if (m_hasCoarser)
       {
         delete specMFCoarOldPtr;
         delete specMFCoarNewPtr;
       } 
   } //end iSpec

  //enforce sum of mass fractions = 1 : Yi* = Yi/sum(Yi)
  for (int iSpec = 0; iSpec < m_nSpec; iSpec++)
    {
      for(DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
       {
         int isrc = 0; int idst = iSpec; int inco = 1;
         MFnew[dit()].divide(sumMF[dit()], isrc, idst, inco);
       } 
    }        

  EBLevelDataOps::scale(MFold, -1./m_dt);
  EBLevelDataOps::scale(MFnew,  1./m_dt);
  LevelData<EBCellFAB> divMF(m_eblg.getDBL(), m_nSpec, nghost*IntVect::Unit, fact);
  for(DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      divMF[dit()].setVal(0.);
      //sets divMF = (MFnew - MFold)/dt
      divMF[dit()] += MFnew[dit()];
      divMF[dit()] += MFold[dit()]; //see scale above
      //sets divMF = rho(MFnew - MFold)/dt
      for (int iSpec = 0; iSpec < m_nSpec; iSpec++)
       {
         int isrc = 0; int idst = iSpec; int inco = 1;
         divMF[dit()].mult((*m_acoDiff[iSpec])[dit()], isrc, idst, inco);
       }

      //now add dt*divMF into ustar
      EBCellFAB dtDivMF(m_eblg.getEBISL()[dit()], m_eblg.getDBL().get(dit()), m_nSpec);
      dtDivMF.setVal(0.);
      dtDivMF += divMF[dit()];
      dtDivMF *= m_dt;

      int isrc = 0; int idst = CSPEC1; int inco = m_nSpec;
      a_UStar[dit()].plus(dtDivMF, isrc, idst, inco);
    }

  if (m_hasCoarser)
   {
     delete MFCoarPtr;
     delete veloCoarPtr;
     delete tempCoarPtr;
   }
}
/***************************/
void EBAMRReactive::
addViscosity(LevelData<EBCellFAB>& a_UStar)
{
  CH_TIME("EBAMRReactive::addViscosity");

  int nghost = 4;
  EBCellFactory fact(m_eblg.getEBISL());
  LevelData<EBCellFAB> divSigma(m_eblg.getDBL(), SpaceDim, nghost*IntVect::Unit, fact);
  LevelData<EBCellFAB> velold(m_eblg.getDBL(), SpaceDim, nghost*IntVect::Unit, fact);
  LevelData<EBCellFAB> velnew(m_eblg.getDBL(), SpaceDim, nghost*IntVect::Unit, fact);
  LevelData<EBCellFAB> rhsZero(m_eblg.getDBL(), SpaceDim, nghost*IntVect::Unit, fact);
  LevelData<EBCellFAB> primOld(m_eblg.getDBL(), m_nPrim, nghost*IntVect::Unit, fact); 
  EBLevelDataOps::setToZero(rhsZero);

  getPrimState(primOld, a_UStar);

  Interval velInt(0,SpaceDim-1);
  Interval srcInt(QVELX, QVELX+SpaceDim-1);

  primOld.copyTo(srcInt, velold, velInt);

  EBFluxRegister*       coarVelFRPtr = NULL;
  EBFluxRegister*       fineVelFRPtr = NULL;
  LevelData<EBCellFAB>* vCoarOldPtr = NULL;
  LevelData<EBCellFAB>* vCoarNewPtr = NULL;
  int ncomp = SpaceDim;
  Real tCoarOld = 0.0;
  Real tCoarNew = 0.0;

  if(m_hasCoarser)
    {
      // Get coarser data if we have it.
      EBAMRReactive* coarPtr = getCoarserLevel();
      coarVelFRPtr = &coarPtr->m_veloFluxRegister;
      tCoarNew = coarPtr->m_time;
      tCoarOld = tCoarNew - coarPtr->m_dt;
      const EBLevelGrid& ceblg = coarPtr->m_eblg;
      EBCellFactory cfact(ceblg.getEBISL());
      vCoarNewPtr = new LevelData<EBCellFAB>(ceblg.getDBL(), ncomp, nghost*IntVect::Unit, cfact);
      vCoarOldPtr = new LevelData<EBCellFAB>(ceblg.getDBL(), ncomp, nghost*IntVect::Unit, cfact);

      LevelData<EBCellFAB> primCoarOld(ceblg.getDBL(), m_nPrim, nghost*IntVect::Unit, cfact);
      LevelData<EBCellFAB> primCoarNew(ceblg.getDBL(), m_nPrim, nghost*IntVect::Unit, cfact);

      coarPtr->getPrimState(primCoarOld, coarPtr->m_stateOld);
      coarPtr->getPrimState(primCoarNew, coarPtr->m_stateNew);

      primCoarOld.copyTo(srcInt, *vCoarOldPtr, velInt);
      primCoarNew.copyTo(srcInt, *vCoarNewPtr, velInt);
    }

  if(m_hasFiner)
    {
      // Get finer flux registers if we have them.
      fineVelFRPtr = &m_veloFluxRegister;
    }  

  Interval srcInt1(CRHO, CRHO);
  Interval dstInt(0, 0);
  for(DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      const Box& region = m_eblg.getDBL().get(dit());
      (*m_acoVisc)[dit()].copy(region, dstInt, region, a_UStar[dit()], srcInt1);
    }

  //backward Euler Integrator:
  s_viscLevBE->updateSoln(velnew, velold, rhsZero,  fineVelFRPtr, coarVelFRPtr,
                          vCoarOldPtr, vCoarNewPtr, m_time, tCoarOld, tCoarNew, m_dt,
                          m_level, true);

  EBLevelDataOps::scale(velold, -1./m_dt);
  EBLevelDataOps::scale(velnew,  1./m_dt);
  for(DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      divSigma[dit()].setVal(0.);
      //sets divsigma = (velold - velold)/dt
      divSigma[dit()] += velnew[dit()];
      divSigma[dit()] += velold[dit()]; //see the scale above
      //sets divsigma = rho(velold - velold)/dt
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          int isrc = 0; int idst = idir; int inco = 1;
          divSigma[dit()].mult((*m_acoVisc)[dit()], isrc, idst, inco);
        }

      //now add divsigma*dt into ustar
      EBCellFAB dtDivSigma(m_eblg.getEBISL()[dit()], m_eblg.getDBL().get(dit()), SpaceDim);
      dtDivSigma.setVal(0.);
      dtDivSigma += divSigma[dit()];
      dtDivSigma *= m_dt;

      int isrc = 0; int idst = CMOMX; int inco = SpaceDim;
      a_UStar[dit()].plus(dtDivSigma, isrc, idst, inco);
    }  

  if(m_hasCoarser)
    {
      delete vCoarOldPtr;
      delete vCoarNewPtr;
    }
}
/***************************/
void EBAMRReactive::
addConductivity(LevelData<EBCellFAB>& a_UStar)
{
  CH_TIME("EBAMRReactive::addConductivity");
  
  // neglecting time variation of rho*Cv
  int nghost = 4;
  EBCellFactory fact(m_eblg.getEBISL());
  LevelData<EBCellFAB>   divSigmaU(m_eblg.getDBL(), 1, nghost*IntVect::Unit, fact);
  LevelData<EBCellFAB> dtDivKGradT(m_eblg.getDBL(), 1, nghost*IntVect::Unit, fact);

  getSingleLdOfU(divSigmaU,  a_UStar);  //UStar gets updated in here
  getDivKappaGradT(dtDivKGradT, a_UStar); //UStar gets updated in here
}
/***************************/
void EBAMRReactive::
getSingleLdOfU(LevelData<EBCellFAB>      & a_divSigmaU,
               LevelData<EBCellFAB>      & a_UStar)
{
  CH_TIME("EBAMRReactive::getSingleLdOfU");

  int nghost = 4;
  EBCellFactory fact(m_eblg.getEBISL());
  LevelData<EBCellFAB> velStar(m_eblg.getDBL(), SpaceDim, nghost*IntVect::Unit, fact);
  LevelData<EBCellFAB> velHalf(m_eblg.getDBL(), SpaceDim, nghost*IntVect::Unit, fact);
  LevelData<EBCellFAB> primStar(m_eblg.getDBL(), m_nPrim, nghost*IntVect::Unit, fact);
  
  getPrimState(primStar, a_UStar);

  Interval velInt(0,SpaceDim-1);
  Interval srcInt(QVELX, QVELX+SpaceDim-1);

  primStar.copyTo(srcInt, velStar, velInt);
  LevelData<EBCellFAB>  kappaConsDivSigmaU(m_eblg.getDBL(), 1,nghost*IntVect::Unit, fact);
  LevelData<EBCellFAB>    nonConsDivSigmaU(m_eblg.getDBL(), 1,nghost*IntVect::Unit, fact);
  
  LevelData<EBCellFAB>* MFCoar = NULL;
  LevelData<EBCellFAB>* veloCoar = NULL;
  LevelData<EBCellFAB>* tempCoar = NULL;

  getCoarserPrimitives(MFCoar, veloCoar, tempCoar);
  s_viscLevBE->getKappaDivSigmaU(kappaConsDivSigmaU,
                                 velStar,
                                 veloCoar,
                                 m_level);

  //now get the non conservative version
  EBLevelDataOps::setToZero(nonConsDivSigmaU);
  EBLevelDataOps::incr     (nonConsDivSigmaU, kappaConsDivSigmaU, 1.0);
  KappaSquareNormal normalizinator(m_eblg);
  normalizinator(nonConsDivSigmaU);

  //(\rho E)^{**} = (\rho E)^* + \dt L^d(U^*)
  //     mass diff = kappa(1-kappa)*(kappaConsDissFcn - nonConsDissFcn)
  updateEnergyBySingleLdAndRedistribute(a_divSigmaU,
                                        kappaConsDivSigmaU,
                                        nonConsDivSigmaU,
                                        a_UStar);

  if(m_hasCoarser)
    {
      delete MFCoar;
      delete veloCoar;
      delete tempCoar;
    }
}
/***************************/
void EBAMRReactive::
updateEnergyBySingleLdAndRedistribute(LevelData<EBCellFAB>      & a_divSigmaU,
                                      LevelData<EBCellFAB>      & a_kappaConsDivSigmaU,
                                      LevelData<EBCellFAB>      & a_nonConsDivSigmaU,
                                      LevelData<EBCellFAB>      & a_UStar)
{
  CH_TIME("EBAMRReactive::updateEnergyBySingleLdAndRedistribute");
  //make the dissipation function = kappa*cons + (1-kappa)*noncons
    EBLevelDataOps::setToZero(a_divSigmaU);
  EBLevelDataOps::incr(a_divSigmaU, a_kappaConsDivSigmaU, 1.0);
  for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      const Box& region = m_eblg.getDBL().get(dit());
      IntVectSet ivs = m_eblg.getEBISL()[dit()].getIrregIVS(region);
      for(VoFIterator vofit(ivs, m_eblg.getEBISL()[dit()].getEBGraph()); vofit.ok(); ++vofit)
        {
          Real kappa   =  m_eblg.getEBISL()[dit()].volFrac(vofit());   // kappa here is volfrac, 
          Real kapCons =  a_kappaConsDivSigmaU[dit()](vofit(), 0);    //     not to be confused with conductivity
          Real nonCons =    a_nonConsDivSigmaU[dit()](vofit(), 0);
          a_divSigmaU[dit()](vofit(), 0) = kapCons + (1.-kappa)*nonCons;
          m_massDiff[dit()](vofit(), CENG) = m_dt*((1.-kappa)*kapCons - kappa*(1.-kappa)*nonCons);
        }
    }
  Interval enerint(CENG, CENG);
  for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      EBCellFAB incr(m_eblg.getEBISL()[dit()], m_eblg.getDBL().get(dit()), 1);
      incr.setVal(0.);
      //this makes incr = dissipation function = divsigma u
      incr += a_divSigmaU[dit()];
      //this makes incr = dt*(divsigmau)
      incr *= m_dt;
      //this adds effects of dissipation into the energy
      int isrc = 0; int idst = CENG; int inco = 1;
      a_UStar[dit()].plus(incr, isrc, idst, inco);

      m_ebLevelRedist.increment(m_massDiff[dit()], dit(), enerint);
      m_massDiff[dit()].setVal(0.);
    }

  //this smushes the energy difference into the state

  a_UStar.exchange(Interval(0, m_nComp-1));
  m_ebLevelRedist.resetWeights(a_UStar, CENG);

  m_ebLevelRedist.redistribute(m_redisRHS, enerint);
  m_ebLevelRedist.setToZero();
}
/***************************/
void EBAMRReactive::
getDivKappaGradT(LevelData<EBCellFAB>& a_dtDivKappaGradT,
                 LevelData<EBCellFAB>& a_UStar)
{
  CH_TIME("EBAMRReactive::getDivKappaGradT");

  int nghost = 4;
  EBCellFactory fact(m_eblg.getEBISL());
  LevelData<EBCellFAB>  Told(m_eblg.getDBL(), 1,nghost*IntVect::Unit, fact);
  LevelData<EBCellFAB>  Tnew(m_eblg.getDBL(), 1,nghost*IntVect::Unit, fact);
  LevelData<EBCellFAB> primOld(m_eblg.getDBL(), m_nPrim, nghost*IntVect::Unit, fact);
  LevelData<EBCellFAB> densOld(m_eblg.getDBL(), 1,nghost*IntVect::Unit, fact);
  LevelData<EBCellFAB> MFold(m_eblg.getDBL(), m_nSpec,nghost*IntVect::Unit, fact);

  getPrimState(primOld, a_UStar);

  Interval tempInt(0,0);
  Interval specInt(0, m_nSpec-1);
  Interval srcInt1(QTEMP, QTEMP);
  Interval srcInt2(QRHO,QRHO);
  Interval srcInt3(QSPEC1, QSPEC1+m_nSpec-1);

  primOld.copyTo(srcInt1, Told, tempInt);
  primOld.copyTo(srcInt2, densOld, tempInt);
  primOld.copyTo(srcInt3, MFold, specInt);

  EBFluxRegister*       coarTempFRPtr = NULL;
  EBFluxRegister*       fineTempFRPtr = NULL;
  LevelData<EBCellFAB>* TCoarOldPtr = NULL;
  LevelData<EBCellFAB>* TCoarNewPtr = NULL;
  Real tCoarOld= 0;
  Real tCoarNew= 0;
  if(m_hasCoarser)
    {
      // Get coarser data if we have it.
      EBAMRReactive* coarPtr = getCoarserLevel();
      coarTempFRPtr = &coarPtr->m_tempFluxRegister;
      tCoarNew = coarPtr->m_time;
      tCoarOld = tCoarNew - coarPtr->m_dt;

      const EBLevelGrid& ceblg = coarPtr->m_eblg;
      EBCellFactory cfact(ceblg.getEBISL());
      TCoarNewPtr = new LevelData<EBCellFAB>(ceblg.getDBL(), 1, nghost*IntVect::Unit, cfact);
      TCoarOldPtr = new LevelData<EBCellFAB>(ceblg.getDBL(), 1, nghost*IntVect::Unit, cfact);

       LevelData<EBCellFAB> primCoarOld(ceblg.getDBL(), m_nPrim, nghost*IntVect::Unit, cfact);
       LevelData<EBCellFAB> primCoarNew(ceblg.getDBL(), m_nPrim, nghost*IntVect::Unit, cfact);

       coarPtr->getPrimState(primCoarOld, coarPtr->m_stateOld);
       coarPtr->getPrimState(primCoarNew, coarPtr->m_stateNew);

       primCoarOld.copyTo(srcInt1, *TCoarOldPtr, tempInt);
       primCoarNew.copyTo(srcInt1, *TCoarNewPtr, tempInt);
    }
  if(m_hasFiner)
    {
      // Get finer flux registers if we have them.
      fineTempFRPtr = &m_tempFluxRegister;
    }
  // Set the a coefficients for the thermal conduction equation to Cv * rho, 
  // specifying both the old and new values so that the density gets linearly
  // interpolated.
  for(DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
   {
     const Box& region = m_eblg.getDBL().get(dit());
     const EBISBox& ebisbox = m_eblg.getEBISL()[dit()];
     EBCellFAB& rho = densOld[dit()];
     const EBCellFAB& massFrac = MFold[dit()];
     const EBCellFAB& temperat = Told[dit()];

     // cfivs, time, timestep fake and not used here
     Real faket = 1.0;
     int logflag = 0;
     IntVectSet emptyivs;
     m_ebPatchReactive->setValidBox(region, ebisbox, emptyivs, faket, faket);
     m_ebPatchReactive->getRhoCv((*m_acoCond)[dit()], rho, massFrac, temperat, region);
   }

  LevelData<EBCellFAB> rhsTemp(m_eblg.getDBL(), 1,  nghost*IntVect::Unit, fact);
  EBLevelDataOps::setToZero(rhsTemp);
  Interval srcComp(CENG, CENG);
  Interval dstComp(0,0);

  m_redisRHS.copyTo(srcComp, rhsTemp, dstComp);
  EBLevelDataOps::setToZero(m_redisRHS);

  //defineSolvers();

  //backward Euler integrator
  s_condLevBE->updateSoln(Tnew, Told, rhsTemp,  fineTempFRPtr, coarTempFRPtr,
                          TCoarOldPtr, TCoarNewPtr, m_time, tCoarOld, tCoarNew, m_dt/2.,
                          m_level, true);

  EBLevelDataOps::scale(Told, -1.0);
  EBLevelDataOps::scale(Tnew,  1.0);

  for(DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      a_dtDivKappaGradT[dit()].setVal(0.);
      //sets divkgradt = ((tnew-told)/dt)*dt
      a_dtDivKappaGradT[dit()] += Told[dit()]; //neg sign is in the scale call above
      a_dtDivKappaGradT[dit()] += Tnew[dit()];

      //this makes divkgradt = (rho cv (tnew-told)/dt)*dt
      a_dtDivKappaGradT[dit()] *= (*m_acoCond)[dit()];

      int isrc = 0; int idst = CENG; int inco = 1;
      a_UStar[dit()].plus(a_dtDivKappaGradT[dit()], isrc, idst, inco);
    }

  if(m_hasCoarser)
    {
      delete TCoarOldPtr;
      delete TCoarNewPtr;
    }
}
/***************************/
void EBAMRReactive::
finalAdvance(LevelData<EBCellFAB>& a_Ustar)
{
  CH_TIME("EBAMRReactive::finalAdvance");
  //set stateNew to udbst
  Interval comps(0, m_nComp-1);
  a_Ustar.copyTo(comps, m_stateNew, comps);

  m_ebLevelReactive.floorConserved(m_stateNew, m_time, m_dt);
}
/***************************/
void EBAMRReactive::
kappaSpecMassDiffSrc(LevelData<EBCellFAB>& a_kappaSpecMassSrc,
                     const LevelData<EBCellFAB>& a_specMF,
                     const LevelData<EBCellFAB>* a_specMFCoar,
                     const LevelData<EBCellFAB>& a_state)
{
  CH_TIME("EBAMRReactive::kappaSpecMassDiffSrc");
  //set the a coefficient to time n
  Interval srcInt(CRHO,CRHO);
  Interval dstInt(0, 0);
  for (int iSpec = 0; iSpec < m_nSpec; iSpec++)
   {
     for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
      { 
        const Box& region = m_eblg.getDBL().get(dit());
        (*m_acoDiff[iSpec])[dit()].copy(region, dstInt, region, a_state[dit()], srcInt);
      }
    }

  //species by species
  for (int iSpec = 0; iSpec < m_nSpec; iSpec++)
   {
     Real alpha = 0; Real beta = 1; // want just the div(flux) part of the operator
     // Compute the mass diffusion term.  coefficient is unity because we want the straight operator
     int nghost = 4;
     EBCellFactory fact(m_eblg.getEBISL());
     LevelData<EBCellFAB> kappaSpecSrc(m_eblg.getDBL(),  1, nghost*IntVect::Unit, fact);
     LevelData<EBCellFAB> massFrac(m_eblg.getDBL(), 1, nghost*IntVect::Unit, fact);
     LevelData<EBCellFAB>* massFracCoar = NULL;
     Interval specInt(iSpec, iSpec);
     Interval thisInt(0, 0);
     a_specMF.copyTo(specInt, massFrac, thisInt);

     if(m_hasCoarser)
      {
        massFracCoar = new LevelData<EBCellFAB>(m_eblg.getDBL(),  1, nghost*IntVect::Unit, fact);
        a_specMFCoar->copyTo(specInt, *massFracCoar, thisInt);
      }

     bool applyBC = true;
     bool explicitHyperbolicSrc = true;  
     s_diffuseLevBE[iSpec]->applyOperator(kappaSpecSrc,
                                          massFrac,
                                          massFracCoar,
                                          m_level, alpha, beta,applyBC);

     bool addMCDiff = true;  //MCDiff for MultiComponent Diffusion

     ParmParse pp;
     if (pp.contains("addMCDiffTerm"))
      {
        pp.get("addMCDiffTerm", addMCDiff);
      }

     if (addMCDiff)
      {
        bool explicitHyperbolicSrc = true;
        LevelData<EBCellFAB> MCDiffTerm(m_eblg.getDBL(), 1, nghost*IntVect::Unit, fact);
        getMCDiffTerm(MCDiffTerm, a_specMF, a_specMFCoar, alpha, beta, iSpec, explicitHyperbolicSrc, applyBC);
        for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
         {
           kappaSpecSrc[dit()] += MCDiffTerm[dit()]; 
         }
      }

     if (m_hasCoarser)
      {
        delete massFracCoar;
      }
 
     kappaSpecSrc.copyTo(thisInt, a_kappaSpecMassSrc, specInt);
   } // end iSpec
}
/***************************/
void EBAMRReactive::
getMCDiffTerm(LevelData<EBCellFAB>& a_MCDiffTerm,
              const LevelData<EBCellFAB>& a_massFrac,
              const LevelData<EBCellFAB>* a_massFracCoar,
              const Real a_alpha,
              const Real a_beta,  
              const int a_iSpec,
              bool a_explicitHyperbolicSrc,
              bool  a_applyBC)
{
  CH_TIME("EBAMRReactive::getMCDiffTerm");
  EBLevelDataOps::setToZero(a_MCDiffTerm);
  int nghost = 4;
  EBCellFactory fact(m_eblg.getEBISL());
  EBFluxFactory fluxFact(m_eblg.getEBISL());
  BaseIVFactory<Real> bivFact(m_eblg.getEBISL(), m_sets);
  Interval thisInterv(0, 0);

  LevelData<EBFluxFAB> tempData(m_eblg.getDBL(), 1, nghost*IntVect::Unit, fluxFact);
  LevelData<BaseIVFAB<Real> > tempDataIrreg(m_eblg.getDBL(), 1, nghost*IntVect::Unit, bivFact);

  //swap pointers, set bco to rhs and set bco back to its original value. So life gets easy
  // not the best way, I know. 
  (*m_bco[a_iSpec]).copyTo(thisInterv, tempData, thisInterv); 
  (*m_bcoIrreg[a_iSpec]).copyTo(thisInterv, tempDataIrreg, thisInterv);

  for (int iSpec = 0; iSpec < m_nSpec; iSpec++)
   {
     Interval specInterv(iSpec,iSpec);
     LevelData<EBCellFAB> kappaSpecMCDiffSrc(m_eblg.getDBL(),  1, nghost*IntVect::Unit, fact);
     LevelData<EBCellFAB> massFrac(m_eblg.getDBL(), 1, nghost*IntVect::Unit, fact);
     LevelData<EBCellFAB>* massFracCoar = NULL; 
     a_massFrac.copyTo(specInterv, massFrac, thisInterv); 
 
     if(m_hasCoarser)
      { 
        massFracCoar = new LevelData<EBCellFAB>(m_eblg.getDBL(),  1, nghost*IntVect::Unit, fact);
        a_massFracCoar->copyTo(specInterv, *massFracCoar, thisInterv); 
      }
     (*m_rhsco[a_iSpec]).copyTo(specInterv, *m_bco[a_iSpec], thisInterv);
     (*m_rhscoIrreg[a_iSpec]).copyTo(specInterv, *m_bcoIrreg[a_iSpec], thisInterv);
     
     //get the kappa*div(flux)
     s_diffuseLevBE[iSpec]->applyOperator(kappaSpecMCDiffSrc,
                                          massFrac,
                                          massFracCoar,
                                          m_level, a_alpha, a_beta, a_applyBC);

     for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
       { 
         a_MCDiffTerm[dit()] += kappaSpecMCDiffSrc[dit()]; 
       } 

     if (m_hasCoarser)
      {
        delete massFracCoar;
      }
   } // end iSpec

  // put the original value of bco in its place
  tempData.copyTo(thisInterv, *m_bco[a_iSpec], thisInterv);
  tempDataIrreg.copyTo(thisInterv, *m_bcoIrreg[a_iSpec], thisInterv);
}
/***************************/
void EBAMRReactive::
kappaMomentumSrc(LevelData<EBCellFAB>& a_kappaDivSigma,
                 const LevelData<EBCellFAB>& a_velocity,
                 const LevelData<EBCellFAB>* a_veloCoar,
                 const LevelData<EBCellFAB>& a_state)
{
  CH_TIME("EBAMRReactive::kappaMomentumSrc");
  //set the a coefficients to time n
  Interval srcInt(CRHO, CRHO);
  Interval dstInt(0, 0);
  for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
   {
     const Box& region = m_eblg.getDBL().get(dit());
      (*m_acoVisc)[dit()].copy(region, dstInt, region, a_state[dit()], srcInt);
    }

  Real alpha = 0;   Real beta = 1; //want just the div(flux) part of the operator
  // Compute the viscous diffusion term.  coefficient is unity because we want the straight operator
  s_viscLevBE->applyOperator(a_kappaDivSigma,
                             a_velocity,
                             a_veloCoar,
                             m_level, alpha, beta, true);
}
/***************************/
void EBAMRReactive::
kappaEnergySrc(LevelData<EBCellFAB>& a_kappaEnergySrc,
               const LevelData<EBCellFAB>& a_velocity,
               const LevelData<EBCellFAB>* a_veloCoar,
               const LevelData<EBCellFAB>& a_temperature,
               const LevelData<EBCellFAB>* a_tempCoar,
               const LevelData<EBCellFAB>& a_specMF,
               const LevelData<EBCellFAB>& a_state)
{
  CH_TIME("EBAMRReactive::kappaEnergySrc");
  //set the rhoCv to time n
  Interval srcInt(CRHO, CRHO);
  Interval dstInt(0,0);
  int nghost = 4;
  EBCellFactory fact(m_eblg.getEBISL());

  LevelData<EBCellFAB> dens(m_eblg.getDBL(), 1, nghost*(IntVect::Unit), fact);
  a_state.copyTo(srcInt, dens, dstInt);

  for(DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
   {
     const Box& region = m_eblg.getDBL().get(dit());
     const EBISBox& ebisbox = m_eblg.getEBISL()[dit()];
     EBCellFAB& rho = dens[dit()];
     const EBCellFAB& massFrac = a_specMF[dit()];
     const EBCellFAB& temperat = a_temperature[dit()];

     // cfivs, time, timestep fake and not used here
     Real faket = 1.0;
     int logflag = 0;
     IntVectSet emptyivs;
     m_ebPatchReactive->setValidBox(region, ebisbox, emptyivs, faket, faket);
     m_ebPatchReactive->getRhoCv((*m_acoCond)[dit()], rho, massFrac, temperat, region);
   }

  // Compute the heat conduction term = volfrac(div kappa grad T)
  Real alpha = 0;   Real beta = 1; //want just the div(flux) part of the operator
  s_condLevBE->applyOperator(a_kappaEnergySrc, a_temperature,  a_tempCoar,
                             m_level, alpha, beta,  true);
 
  /**/
  /// add volfrac*div( sigma u) to a_kappaEnergySource 
  LevelData<EBCellFAB>  kappaDivSigmaU(m_eblg.getDBL(), 1,nghost*IntVect::Unit, fact);
  s_viscLevBE->getKappaDivSigmaU(kappaDivSigmaU,
                                 a_velocity,
                                 a_veloCoar,
                                 m_level);


  for(DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      a_kappaEnergySrc[dit()] += kappaDivSigmaU[dit()];
    }
  /**/
}
/***************************/
void EBAMRReactive::
getCoarserPrimitives(LevelData<EBCellFAB>* &  a_specMFCoar,
                     LevelData<EBCellFAB>* &  a_veloCoar,
                     LevelData<EBCellFAB>* &  a_tempCoar)
{
  CH_TIME("EBAMRReactive::getCoarserPrimitives");
  LevelData<EBCellFAB>* specMFCoar = NULL;
  LevelData<EBCellFAB>* veloCoar = NULL;
  LevelData<EBCellFAB>* tempCoar = NULL;

  if (m_hasCoarser)
   {
     EBAMRReactive* coarPtr = getCoarserLevel();
     EBCellFactory factCoar(coarPtr->m_eblg.getEBISL());
     int nghost = 4;
     specMFCoar = new LevelData<EBCellFAB>(coarPtr->m_eblg.getDBL(), m_nSpec, nghost*IntVect::Unit, factCoar);
     veloCoar = new LevelData<EBCellFAB>(coarPtr->m_eblg.getDBL(), SpaceDim, nghost*IntVect::Unit, factCoar);
     tempCoar = new LevelData<EBCellFAB>(coarPtr->m_eblg.getDBL(),        1, nghost*IntVect::Unit, factCoar);

      LevelData<EBCellFAB> consCoar(coarPtr->m_eblg.getDBL(), m_nComp, nghost*IntVect::Unit, factCoar);
      LevelData<EBCellFAB> primCoar(coarPtr->m_eblg.getDBL(), m_nPrim, nghost*IntVect::Unit, factCoar);
  
      Real tCoarNew = coarPtr->m_time;
      Real tCoarOld = tCoarNew - coarPtr->m_dt;

      //interpolate coarse solution to fine time
      EBArith::timeInterpolate(consCoar, coarPtr->m_stateOld, coarPtr->m_stateNew, coarPtr->m_eblg.getDBL(), m_time, tCoarOld, tCoarNew);

      getPrimState(primCoar, consCoar);

      Interval specSrcInterv(QSPEC1, QSPEC1+m_nSpec-1);
      Interval veloSrcInterv(QVELX, QVELX+SpaceDim-1);
      Interval tempSrcInterv(QTEMP, QTEMP);
      Interval specDstInterv(0, m_nSpec-1);
      Interval veloDstInterv(0, SpaceDim-1);
      Interval tempDstInterv(0, 0);
   
      primCoar.copyTo(specSrcInterv, *specMFCoar, specDstInterv);
      primCoar.copyTo(veloSrcInterv, *veloCoar, veloDstInterv);
      primCoar.copyTo(tempSrcInterv, *tempCoar, tempDstInterv);
   } // end hasCoarser

  a_specMFCoar = specMFCoar;
  a_veloCoar = veloCoar;
  a_tempCoar = tempCoar; 
}
/***************************/
void EBAMRReactive::coarseFineIncrement()
{
  Interval interv(0, m_nComp-1);
  // increment redistribution register between this level
  // and next coarser level
  {
    CH_TIME("coarse-fine guano");
    if (m_hasCoarser)
      {
        for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
          {
            m_ebFineToCoarRedist.increment(m_massDiff[dit()], dit(), interv);
          }
      }
 
    //initialize redistirbution register between this level
    // and next finer level
    //this includes re-redistribution registers
    if (m_hasFiner)
      {
 
        m_ebCoarToFineRedist.setToZero();
        m_ebCoarToCoarRedist.setToZero();
        for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
          {
            BaseIVFAB<Real>& massDiffFAB = m_massDiff[dit()];
            m_ebCoarToFineRedist.increment(massDiffFAB, dit(), interv);
            m_ebCoarToCoarRedist.increment(massDiffFAB, dit(), interv);
          }
      }
  }
}
/***************************/
void
EBAMRReactive::dumpDebug()
{
  dumpDebug("arg");
}
/***************************/
void
EBAMRReactive::dumpDebug(const string& a_debstring)
{
  int ilev = m_level;
  //  if (ilev == debuglevel)
  if (0)
    {
      pout()    << setprecision(10)
                << setiosflags(ios::showpoint)
                << setiosflags(ios::scientific);
      for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
        {
          const Box& grid =m_grids.get(dit());
          if (grid.contains(ivdebamrg))
            {
              pout() << a_debstring << ", lev = " << ilev << ", iv=" << ivdebamrg<<  ":";
              for (int ivar = 0; ivar < m_stateNew.nComp(); ivar++)
                {
                  VolIndex vof(ivdebamrg, 0);
                  Real sold = m_stateOld[dit()](vof, ivar);
                  Real snew = m_stateNew[dit()](vof, ivar);
                  pout() << " ( " << snew << ", "  << sold << " ) ";

                }
              pout() << endl;
            }
        }
    }
}
/***************************/
void
EBAMRReactive::sumConserved(Real& a_sumcons,
                           const int& a_ivar) const
{
  CH_TIME("EBAMRR::sumConserved");
  Real sumgrid = 0;
  for (DataIterator dit= m_grids.dataIterator(); dit.ok(); ++dit)
    {
      const EBCellFAB& state = m_stateNew[dit()];
      Box thisBox = m_grids.get(dit());
      IntVectSet uberIVS(thisBox);
      const EBISBox& ebisBox = m_ebisl[dit()];
      for (VoFIterator vofit(uberIVS, ebisBox.getEBGraph());
          vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          Real consVal = state(vof, a_ivar);
          Real volFrac = ebisBox.volFrac(vof);
          Real volume = volFrac;
          if (m_doRZCoords)
            {
              Real cellVol, kvol;
              EBArith::getKVolRZ(kvol, cellVol, ebisBox, m_dx[0], vof);
              volume = cellVol*kvol;
            }
          sumgrid += consVal*volume;
        }
    }

  Vector<Real> all_sum;
  gather(all_sum,sumgrid,uniqueProc(SerialTask::compute));
  Real sumallgrid = 0.;
  if (procID() == uniqueProc(SerialTask::compute))
    {
      for (int i = 0; i < all_sum.size(); ++i)
        {
          sumallgrid += all_sum[i];
        }
    }
  broadcast(sumallgrid,uniqueProc(SerialTask::compute));
  a_sumcons = sumallgrid;
}
/***************************/
void EBAMRReactive::postTimeStep()
{
  CH_TIME("EBAMRReactive::postTimeStep");
  if (s_verbosity >= 3)
    {
      pout() << " in EBAMRReactive postTimeStep for level " << m_level << endl;
    }

  // Average the finer grid level to be consistent with this one

  if (m_hasFiner)
    {
      Interval interv(0, m_nComp-1);
      EBAMRReactive* finePtr = getFinerLevel();
      finePtr->m_ebCoarseAverage.average(m_stateNew,
                                         finePtr->m_stateNew,
                                         interv);
    }

  // this does the refluxing and redistribution evil dance
  postTimeStepRefluxRedistDance();
  m_ebLevelReactive.floorConserved(m_stateNew, m_time, m_dt);

  for(DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      m_massDiff[dit()].setVal(0.0);
      m_redisRHS[dit()].setVal(0.0);
    } 
}
/****************************/
void EBAMRReactive::
postTimeStepRefluxRedistDance()
{
  CH_TIME("EBAMRReactive::evilRedistDance");
  resetWeights();
  refluxRedistInteraction();
  coarseFineRedistribution(Interval(0,m_nComp-1));
  explicitReflux(Interval(0, m_nComp-1));
  ParmParse pp;
  bool turn_off_reflux = false;
  pp.query("turn_off_implicit_reflux", turn_off_reflux);
  if (turn_off_reflux)
    {
      pout() << "implicit reflux turned off" << endl;
    }
  if (addDiffusion() && !turn_off_reflux)
    {
      refluxRHSConserved();
      // defineSolvers();
      implicitReflux();
    }
}
/****************************/
void EBAMRReactive::
implicitReflux()
{
  CH_TIME("EBAMRReactive::implicitReflux");
  pout() << "EBAMRReactive::implicit redist/reflux for level" << m_level << endl;
  Real crseTime = -1.0;
  Real timeEps = 1.0e-2*m_dt;
  if (m_level > 0) crseTime = m_coarser_level_ptr->time();

  if (m_level == 0 || (abs(crseTime - m_time) < timeEps))
    {
      int finestLev = getFinestLevel();
      //if I have not caught up to the next coarser level,
      //then I need to do implicit{reflux,redist} for my level and all
      //finer levels
      Vector<AMRLevel*> hierarchy = AMRLevel::getAMRLevelHierarchy();
      Vector<LevelData<EBCellFAB>* >    deltaVelocity(finestLev+1, NULL);
      Vector<LevelData<EBCellFAB>* >    deltaTemperat(finestLev+1, NULL);
      Vector<LevelData<EBCellFAB>* > dtRefluxDivergeM(finestLev+1, NULL);
      Vector<LevelData<EBCellFAB>* > dtRefluxDivergeE(finestLev+1, NULL);
   
      // species stuff
      Vector<LevelData<EBCellFAB>* >  deltaMassFrac(finestLev+1, NULL);
      Vector<LevelData<EBCellFAB>* >  dtRefluxDivergeRHO(finestLev+1, NULL);
 
      for (int ilev = 0; ilev <= finestLev; ilev++)
        {
          EBAMRReactive* reactiveLevel = dynamic_cast<EBAMRReactive*>(hierarchy[ilev]);

          EBCellFactory fact(reactiveLevel->m_eblg.getEBISL());
          deltaVelocity[ilev]      = new LevelData<EBCellFAB>(reactiveLevel->m_eblg.getDBL(), SpaceDim, 4*IntVect::Unit, fact);
          deltaTemperat[ilev]      = new LevelData<EBCellFAB>(reactiveLevel->m_eblg.getDBL(),        1, 4*IntVect::Unit, fact);
          dtRefluxDivergeE[ilev]   = new LevelData<EBCellFAB>(reactiveLevel->m_eblg.getDBL(),        1, 4*IntVect::Unit, fact);
          dtRefluxDivergeM[ilev]   = new LevelData<EBCellFAB>(reactiveLevel->m_eblg.getDBL(), SpaceDim, 4*IntVect::Unit, fact);
          deltaMassFrac[ilev]      = new LevelData<EBCellFAB>(reactiveLevel->m_eblg.getDBL(),  m_nSpec, 4*IntVect::Unit, fact);
          dtRefluxDivergeRHO[ilev] = new LevelData<EBCellFAB>(reactiveLevel->m_eblg.getDBL(),  m_nSpec, 4*IntVect::Unit, fact);

          //copy redistRHS to reflux divergence holders
          Interval srcComp, dstComp;
          srcComp = Interval(CMOMX, CMOMX+SpaceDim-1);
          dstComp = Interval(0, SpaceDim-1);
          reactiveLevel->m_redisRHS.copyTo(srcComp, *dtRefluxDivergeM[ilev], dstComp);
          srcComp = Interval(CENG, CENG);
          dstComp = Interval(0, 0);
          reactiveLevel->m_redisRHS.copyTo(srcComp, *dtRefluxDivergeE[ilev], dstComp);
          srcComp = Interval(CSPEC1, CSPEC1+m_nSpec-1);
          dstComp = Interval(0, m_nSpec-1);
          reactiveLevel->m_redisRHS.copyTo(srcComp, *dtRefluxDivergeRHO[ilev], dstComp);

          EBLevelDataOps::scale(*dtRefluxDivergeM[ilev], reactiveLevel->m_dt);
          EBLevelDataOps::scale(*dtRefluxDivergeE[ilev], reactiveLevel->m_dt);
          EBLevelDataOps::scale(*dtRefluxDivergeRHO[ilev], reactiveLevel->m_dt);

          //these MATTER especially with subcycling
          EBLevelDataOps::setVal(*deltaVelocity[ilev], 0.0);
          EBLevelDataOps::setVal(*deltaTemperat[ilev], 0.0);
          EBLevelDataOps::setVal(*deltaMassFrac[ilev], 0.0); 
        }  // end for level loop
  
      pout() << "getting spec continuity reflux/redist increment" << endl;
      int baseLev = Max(m_level-1, 0);
      Real baseDt = ((EBAMRReactive*)hierarchy[baseLev])->m_dt;

      //(rho I - dt Ly) delta = dt* Dr(Frho) 
      getRefluxDeltaMassFrac(deltaMassFrac, dtRefluxDivergeRHO, baseLev, finestLev, baseDt);

      // rhoi += dt* Ly(deltaY)
      incrDenseByDeltaY(deltaMassFrac, baseLev, finestLev);

      pout() << "getting momentum reflux/redist increment"  << endl;
      //(rho I - dt Lv) delta = dt* Dr(Fm)
      getRefluxDeltaV(  deltaVelocity, dtRefluxDivergeM, baseLev, finestLev, baseDt);

      //mom += dt* Lv(deltav)
      incrMomentumByDeltaV(deltaVelocity, baseLev, finestLev);

      pout() << "getting energy reflux/redist increment"  << endl;
      //(rho I - dt Lt) delta =dt* Dr(Fe)
      getRefluxDeltaT(deltaTemperat, dtRefluxDivergeE, baseLev, finestLev, baseDt);

      //ene += dt* Lt(deltaT)
      incrEnergyByDeltaT(deltaTemperat, baseLev, finestLev);

      for(int ilev = 0; ilev <= finestLev; ilev++)
        {
          delete    deltaVelocity[ilev];
          delete    deltaTemperat[ilev];
          delete    deltaMassFrac[ilev];
          delete dtRefluxDivergeM[ilev];
          delete dtRefluxDivergeE[ilev];
          delete dtRefluxDivergeRHO[ilev];
        }

      //reset redistribution RHS to avoid double counting
      for(int ilev = m_level; ilev <= finestLev; ilev++)
        {
          EBAMRReactive* reactiveLevel = dynamic_cast<EBAMRReactive*>(hierarchy[ilev]);
          EBLevelDataOps::setToZero(reactiveLevel->m_redisRHS);
        }
    }// end if m_level == 0        
}
/****************************/
// (rho I - dt Ly) delta = dt*Dr(Frho)
void EBAMRReactive::
getRefluxDeltaMassFrac(Vector<LevelData<EBCellFAB>* >& a_deltaMassFrac,
                       Vector<LevelData<EBCellFAB>* >& a_dtRefluxDivergeRHO,
                       int a_baseLev, int a_finestLev, Real a_dtBase)
{
  CH_TIME("EBAMRReactive::getRefluxDeltaMassFrac");
  Vector<AMRLevel*> hierarchy = AMRLevel::getAMRLevelHierarchy();
  //set identity coefficients
  for (int ilev = 0; ilev <= a_finestLev; ilev++)
    {
      EBAMRReactive* reactiveLevel = dynamic_cast<EBAMRReactive*>(hierarchy[ilev]);
      LevelData<EBCellFAB>& state = reactiveLevel->m_stateNew;
      Vector<RefCountedPtr<LevelData<EBCellFAB> > >& acoef = reactiveLevel->m_acoDiff;
      Interval srcInt(CRHO, CRHO);
      Interval dstInt(0, 0);

      
      for (int iSpec = 0; iSpec < m_nSpec; iSpec++)
       {
         for (DataIterator dit = reactiveLevel->m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
          {
             const Box& region = reactiveLevel->m_eblg.getDBL().get(dit());
             (*acoef[iSpec])[dit()].copy(region, dstInt, region, state[dit()], srcInt);
          }
       }
    }

  for (int iSpec = 0; iSpec < m_nSpec; iSpec++)
   {
     //set alpha to 1 and beta = -dtbase
     s_diffuseLevBE[iSpec]->resetSolverAlphaAndBeta(1.0, -a_dtBase);
     
     //the sad way of doing it :-(
     Vector<LevelData<EBCellFAB>* > deltaSpecMF(a_finestLev+1, NULL);
     Vector<LevelData<EBCellFAB>* > dtRefluxDivergeRHO(a_finestLev+1, NULL);
      
     for (int ilev = 0; ilev <= a_finestLev; ilev++)
      {
        EBAMRReactive* reactiveLevel = dynamic_cast<EBAMRReactive*>(hierarchy[ilev]);
        EBCellFactory fact(reactiveLevel->m_eblg.getEBISL());
        deltaSpecMF[ilev]        = new LevelData<EBCellFAB>(reactiveLevel->m_eblg.getDBL(), 1, 4*IntVect::Unit, fact);  
        dtRefluxDivergeRHO[ilev] = new LevelData<EBCellFAB>(reactiveLevel->m_eblg.getDBL(), 1, 4*IntVect::Unit, fact);
        
        Interval srcComp, dstComp;
        srcComp = Interval(iSpec, iSpec);
        dstComp = Interval(0, 0);
        a_dtRefluxDivergeRHO[ilev]->copyTo(srcComp, *dtRefluxDivergeRHO[ilev], dstComp); 
      }
 
     //solve equation (rho I - dt Ly) delta = dt*Dr(Frho)
     //rhs already multiplied by dt
     //first true = zero phi
     //second true = force homogeneous bcs (solving for delta Y)
     s_diffuseAMRMG[iSpec]->solve(deltaSpecMF, dtRefluxDivergeRHO, a_finestLev, a_baseLev, true, true);

     for (int ilev = 0; ilev <= a_finestLev; ilev++)
      {
        Interval srcComp, dstComp; 
        dstComp = Interval(iSpec, iSpec); 
        srcComp = Interval(0, 0); 
        deltaSpecMF[ilev]->copyTo(srcComp, *a_deltaMassFrac[ilev], dstComp);
        
        delete deltaSpecMF[ilev];
        delete dtRefluxDivergeRHO[ilev];
      }
   }  
}
/****************************/
// (rho I - dt Lv) delta = dt*Dr(Fm) 
void EBAMRReactive::
getRefluxDeltaV(Vector<LevelData<EBCellFAB>* >&  a_deltaVelocity,
                Vector<LevelData<EBCellFAB>* >&  a_dtRefluxDivergeM,
                int  a_baseLev, int a_finestLev, Real a_dtBase)
{
  CH_TIME("EBAMRReactive::getRefluxDeltaV");
  Vector<AMRLevel*> hierarchy = AMRLevel::getAMRLevelHierarchy();
  //set identity coefficients
  for(int ilev = 0; ilev <= a_finestLev; ilev++)
    {
      EBAMRReactive* reactiveLevel = dynamic_cast<EBAMRReactive*>(hierarchy[ilev]);
      LevelData<EBCellFAB>& state =  reactiveLevel->m_stateNew;
      LevelData<EBCellFAB>& acoef = *reactiveLevel->m_acoVisc;
      Interval srcInt(CRHO, CRHO);
      Interval dstInt(0, 0);
      for(DataIterator dit = reactiveLevel->m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
        {
          const Box& region = reactiveLevel->m_eblg.getDBL().get(dit());
          acoef[dit()].copy(region, dstInt, region, state[dit()], srcInt);
        }
    }

  //set alpha to 1, beta = -dtbase
  s_viscLevBE->resetSolverAlphaAndBeta(1.0, -a_dtBase);

  //solve equation (rho I - dt Lv) delta = dt*Dr(Fm)
  //rhs already multiplied by dt
  //first true = zero phi
  //second true = force homogeneous bcs (solving for a delta v)
  s_viscAMRMG->solve(a_deltaVelocity, a_dtRefluxDivergeM, a_finestLev, a_baseLev, true, true);
} 
/****************************/
//(rho Cv I - dt Lt) delta = dt*Dr(Fe)
void EBAMRReactive::
getRefluxDeltaT(Vector<LevelData<EBCellFAB>* >& a_deltaTemperat,
                      Vector<LevelData<EBCellFAB>* >& a_dtRefluxDivergeE,
                      int a_baseLev, int a_finestLev, Real a_dtBase)
{
  CH_TIME("EBAMRReactive::getRefluxDeltaT"); 
  Vector<AMRLevel*> hierarchy = AMRLevel::getAMRLevelHierarchy();
  //set identity coefficients
  for(int ilev = 0; ilev <= a_finestLev; ilev++)
    {
      EBAMRReactive* reactiveLevel = dynamic_cast<EBAMRReactive*>(hierarchy[ilev]);
      Interval srcInt(CRHO, CRHO);
      Interval dstInt(0, 0);
      LevelData<EBCellFAB>& state =  reactiveLevel->m_stateNew;
      LevelData<EBCellFAB>& acoef = *reactiveLevel->m_acoCond;

      const DisjointBoxLayout& dbl = reactiveLevel->m_eblg.getDBL();

      EBCellFactory fact(reactiveLevel->m_eblg.getEBISL());
      int nghost = 4;

      LevelData<EBCellFAB> primState(dbl, m_nPrim, (IntVect::Unit)*nghost, fact);
      LevelData<EBCellFAB> levelMassFrac(dbl, m_nSpec, (IntVect::Unit)*nghost, fact);
      LevelData<EBCellFAB> levelDense(dbl, 1, (IntVect::Unit)*nghost, fact);
      LevelData<EBCellFAB> levelTemperature(dbl, 1, (IntVect::Unit)*nghost, fact);    

      getPrimState(primState, state);

      EBLevelDataOps::setVal(levelMassFrac,0e0);
      EBLevelDataOps::setVal(levelTemperature,0e0);
      EBLevelDataOps::setVal(levelDense,0e0);
 
      EBPatchReactive* patch = reactiveLevel->m_ebPatchReactive;
      int massFracIndex = patch->spec1MassFracIndex();
      int denseIndex = patch->densityIndexC();
      int tempIndex = patch->temperatureIndex();
      Interval specInterv(massFracIndex,massFracIndex+m_nSpec-1);
      Interval denseInterv(denseIndex,denseIndex);
      Interval tempInterv(tempIndex,tempIndex);
      Interval thisInterv1(0,m_nSpec-1);
      Interval thisInterv2(0,0);

      primState.copyTo(denseInterv,levelDense,thisInterv2);
      primState.copyTo(tempInterv,levelTemperature,thisInterv2);
      primState.copyTo(specInterv,levelMassFrac,thisInterv1);     

      for(DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
        {
          const Box& region   = dbl.get(dit());
          const EBISBox& ebisbox = reactiveLevel->m_eblg.getEBISL()[dit()];
//          EBCellFAB& rhoCv    = acoef[dit()];
          EBCellFAB& rho      = levelDense[dit()];
          EBCellFAB& massFrac = levelMassFrac[dit()];
          EBCellFAB& temperat = levelTemperature[dit()];

          // cfivs, time, timestep fake and not used here
          Real faket = 1.0;
          int logflag = 0;
          IntVectSet emptyivs;
          patch->setValidBox(region, ebisbox, emptyivs, faket, faket);
          patch->getRhoCv(acoef[dit()], rho, massFrac, temperat, region);     
        }
    }
  //set alpha to 1, beta = -dtbase
  s_condLevBE->resetSolverAlphaAndBeta(1.0, -a_dtBase);

  //solve equation (rho C_v I - dt Lv) delta = dt*Dr(Fm)
  //rhs already multiplied by dt
  //first true = zero phi
  //second true = force homogeneous bcs (solving for a delta t)
  s_condAMRMG->solve(a_deltaTemperat, a_dtRefluxDivergeE, a_finestLev, a_baseLev, true, true);
}
/****************************/
//rhoi += dt* Ly(deltaY)
// deltaRHOi gets multiplied by rho
void EBAMRReactive::
incrDenseByDeltaY(Vector<LevelData<EBCellFAB>* >& a_deltaMassFrac, 
                  int a_baseLev, int a_finestLev)
{
  CH_TIME("EBAMRReactive::incrDenseByDeltaY");
  Vector<AMRLevel*> hierarchy = AMRLevel::getAMRLevelHierarchy();
  for(int ilev = a_baseLev; ilev <= a_finestLev; ilev++)
    {
      EBAMRReactive* reactiveLevel   = dynamic_cast<EBAMRReactive*>(hierarchy[ilev]);
      reactiveLevel->incrDenseByDeltaY(reactiveLevel->m_stateNew, *a_deltaMassFrac[ilev]);
    }
}
/****************************/
void EBAMRReactive::
incrDenseByDeltaY(LevelData<EBCellFAB>& a_state,
                  LevelData<EBCellFAB>& a_deltaMassFrac)
{
  CH_TIME("EBAMRReactive::incrDenseByDeltaYLevel");
  for(DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      const Box& b = m_eblg.getDBL()[dit()];
      const EBISBox& ebisBox = m_eblg.getEBISL()[dit()];
      EBCellFAB deltaRhoSpec(ebisBox, b, m_nSpec);
      deltaRhoSpec.setVal(0.);
      deltaRhoSpec += a_deltaMassFrac[dit()];

      EBCellFAB& stateNew  = a_state[dit()];
      //delta rhoi now = dYi 
      //make delta = rho dYi
      for(int iSpec = 0; iSpec <m_nSpec; iSpec++)
        {
          int isrc = CRHO; int idst = iSpec; int inco = 1;
          deltaRhoSpec.mult(stateNew, isrc, idst, inco);
        }
      //add change in species density to the state
      {
        int isrc = 0; int idst = CSPEC1; int inco = m_nSpec;
        stateNew.plus(deltaRhoSpec, isrc, idst, inco);
      }
    }
}
/****************************/
//mom += dt* Lv(deltav)
// deltaM gets multiplied by rho
void EBAMRReactive::
incrMomentumByDeltaV(Vector<LevelData<EBCellFAB>* >& a_deltaVelocity,
                     int a_baseLev, int a_finestLev)
{
  CH_TIME("EBAMRReactive::incrMomentumByDeltaV");
  Vector<AMRLevel*> hierarchy = AMRLevel::getAMRLevelHierarchy();
  for(int ilev = a_baseLev; ilev <= a_finestLev; ilev++)
    {
      EBAMRReactive* reactiveLevel   = dynamic_cast<EBAMRReactive*>(hierarchy[ilev]);
      reactiveLevel->incrMomentumByDeltaV(reactiveLevel->m_stateNew, *a_deltaVelocity[ilev]);
    }
}
/****************************/
void EBAMRReactive::
incrMomentumByDeltaV(LevelData<EBCellFAB>& a_state, 
                     LevelData<EBCellFAB>& a_deltaVelocity)
  
{ 
  CH_TIME("EBAMRReactive::incrMomentumByDeltaVLevel");
  for(DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      const Box& b = m_eblg.getDBL()[dit()];
      const EBISBox& ebisBox = m_eblg.getEBISL()[dit()];
      EBCellFAB deltaMom(ebisBox, b, SpaceDim);
      deltaMom.setVal(0.);
      deltaMom += a_deltaVelocity[dit()];

      EBCellFAB& stateNew  = a_state[dit()];
      //delta mom now = dV 
      //make delta = rho dV
      for(int idir = 0; idir <SpaceDim; idir++)
        {
          int isrc = CRHO; int idst = idir; int inco = 1;
          deltaMom.mult(stateNew, isrc, idst, inco);
        }
      //add change in momentum to the state
      {
        int isrc = 0; int idst = CMOMX; int inco = SpaceDim;
        stateNew.plus(deltaMom, isrc, idst, inco);
      }
    }
}
/****************************/
//ene += dt* Lt(deltaT)
// deltaT gets multiplied by rho cv
void EBAMRReactive::
incrEnergyByDeltaT(Vector<LevelData<EBCellFAB>* >& a_deltaTemperat,
                   int a_baseLev, int a_finestLev) 
{
  CH_TIME("EBAMRReactive::incrEnergyByDeltaT");
  Vector<AMRLevel*> hierarchy = AMRLevel::getAMRLevelHierarchy();
  for(int ilev = a_baseLev; ilev <= a_finestLev; ilev++)
    {
      EBAMRReactive* reactiveLevel   = dynamic_cast<EBAMRReactive*>(hierarchy[ilev]);
      reactiveLevel->incrEnergyByDeltaT(reactiveLevel->m_stateNew, *a_deltaTemperat[ilev]);
    }
}
/****************************/
//ene += dt* Lt(deltaT)
// deltaT gets multiplied by rho cv
void EBAMRReactive::
incrEnergyByDeltaT(LevelData<EBCellFAB>& a_state,
                   LevelData<EBCellFAB>& a_deltaTemperat)
{ 
  CH_TIME("EBAMRReactive::incrEnergyByDeltaTLevel");

  const DisjointBoxLayout& dbl = m_eblg.getDBL();
  EBCellFactory fact(m_eblg.getEBISL());
  int nghost = 4; 

  LevelData<EBCellFAB> primState(dbl, m_nPrim, (IntVect::Unit)*nghost, fact);
  LevelData<EBCellFAB> levelSpecMassFrac(dbl, m_nSpec, (IntVect::Unit)*nghost, fact);
  LevelData<EBCellFAB> levelDense(dbl, 1, (IntVect::Unit)*nghost, fact);
  LevelData<EBCellFAB> levelTemperature(dbl, 1, (IntVect::Unit)*nghost, fact);
  LevelData<EBCellFAB> levelrhoCv(dbl, 1, (IntVect::Unit)*nghost, fact);
 
  getPrimState(primState, a_state);

  EBLevelDataOps::setVal(levelSpecMassFrac,0e0);
  EBLevelDataOps::setVal(levelTemperature,0e0);
  EBLevelDataOps::setVal(levelDense,0e0);
  
  EBPatchReactive* patch = m_ebPatchReactive;
  int massFracIndex = patch->spec1MassFracIndex();
  int denseIndex = patch->densityIndexC();
  int tempIndex = patch->temperatureIndex();
  Interval specInterv(massFracIndex,massFracIndex+m_nSpec-1);
  Interval denseInterv(denseIndex,denseIndex);
  Interval tempInterv(tempIndex,tempIndex);
  Interval thisInterv1(0,m_nSpec-1);
  Interval thisInterv2(0,0);

  primState.copyTo(denseInterv,levelDense,thisInterv2);
  primState.copyTo(tempInterv,levelTemperature,thisInterv2);
  primState.copyTo(specInterv,levelSpecMassFrac,thisInterv1);

  for(DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
    {
      const Box& b = dbl[dit()];
      const EBISBox& ebisBox = m_eblg.getEBISL()[dit()];
      EBCellFAB deltaEnergy(ebisBox, b, 1);
      deltaEnergy.setVal(0.);
      deltaEnergy += a_deltaTemperat[dit()];

      //deltaEnergy now = dT 
      //make deltaEnergy = rho  Cv dT
      {
        EBCellFAB& rhoCv    = levelrhoCv[dit()];
        EBCellFAB& rho      = levelDense[dit()];
        EBCellFAB& massFrac = levelSpecMassFrac[dit()]; 
        EBCellFAB& temperat = levelTemperature[dit()];

        // cfivs, time, timestep fake and not used here
        Real faket = 1.0; 
        int logflag = 0;
        IntVectSet emptyivs; 
        patch->setValidBox(b, ebisBox, emptyivs, faket, faket); 
        patch->getRhoCv(rhoCv, rho, massFrac, temperat, b);
 
        int isrc = 0; int idst = 0; int inco = 1;
        deltaEnergy.mult(rhoCv, isrc, idst, inco);       
      }
      //add change in energy to the state
      {
        int isrc = 0; int idst = CENG; int inco = 1;
        a_state[dit()].plus(deltaEnergy, isrc, idst, inco);
      }
    }
}
/****************************/
int EBAMRReactive::
getFinestLevel()
{
  int imax = 0;
  Vector<AMRLevel*> hierarchy = AMRLevel::getAMRLevelHierarchy();
  for(int ilev = 0; ilev < hierarchy.size(); ilev++)
    {
      EBAMRReactive* reactiveLevel = dynamic_cast<EBAMRReactive*>(hierarchy[ilev]);
      int numBoxes = 0;
      if(reactiveLevel->m_eblg.getDBL().isClosed()) numBoxes = reactiveLevel->m_eblg.getDBL().size();
      if(numBoxes > 0)
        imax++;
      else
        break;
    }
  return imax-1;
}
/****************************/
void EBAMRReactive::
refluxRHSConserved()
{
  CH_TIME("EBAMRReactive::refluxRHSConserved");
  // this does the refluxing and redistribution evil dance
  Interval interm(CMOMX, CMOMX+SpaceDim-1);
  Interval intere(CENG, CENG);
  Interval specInterv(CSPEC1, CSPEC1+m_nSpec-1);
  if (m_hasFiner)
    {
      // reflux from finer level solution
      CH_assert(Abs(m_dx[0]-m_dx[1])<1.e-9);
      Real scale = -1.0/m_dx[0];
      Interval interfluxm(0, SpaceDim-1);
      Interval interfluxe(0, 0);
      m_veloFluxRegister.reflux(m_redisRHS, interm, interfluxm, scale);
      m_tempFluxRegister.reflux(m_redisRHS, intere, interfluxe, scale);

      m_veloFluxRegister.setToZero();
      m_tempFluxRegister.setToZero();

      Interval interfluxspec(0, m_nSpec-1);
      m_specFluxRegister.reflux(m_redisRHS, specInterv, interfluxspec, scale);
      m_specFluxRegister.setToZero();
    }
}
/****************************/
void EBAMRReactive::
explicitReflux(const Interval& a_interv)
{
  CH_TIME("EBAMRReactive::explicitReflux");
  if (m_hasFiner)
   {
     Real scale = -1.0/m_dx[0];
     m_ebFluxRegister.reflux(m_stateNew, a_interv, scale);
     m_ebFluxRegister.setToZero();
   }
}
/****************************/
void EBAMRReactive::
coarseFineRedistribution(const Interval& a_interv)
{
  CH_TIME("EBAMRReactive::coarseFineRedistribution");
  if (m_hasCoarser)
   {
      if (m_doSmushing && (!s_noEBCF))
        {
          // redistribute to coarser level
          EBAMRReactive* coarPtr = getCoarserLevel();
          // put mass directly into state
          m_ebFineToCoarRedist.redistribute(coarPtr->m_stateNew, a_interv);
          m_ebFineToCoarRedist.setToZero();
         }
    }
  if (m_hasFiner)
    {
      EBAMRReactive* finePtr = getFinerLevel();
      if(m_doSmushing && (!s_noEBCF))
        {
          //redistibute to finer level
          m_ebCoarToFineRedist.redistribute(finePtr->m_stateNew, a_interv);
          //do the re redistirubtion
          m_ebCoarToCoarRedist.redistribute(         m_stateNew, a_interv);
          m_ebCoarToFineRedist.setToZero();
          m_ebCoarToCoarRedist.setToZero();
        }
    }
}
/****************************/
void EBAMRReactive::
resetWeights()
{
  CH_TIME("EBAMRReactive::resetWeights");
  // if use mass weighting, need to fix weights of redistribution objects
  if (m_useMassRedist)
   {
     int densevar = m_ebPatchReactive->densityIndex();
     m_stateNew.exchange(Interval(0,m_nComp-1));
     m_ebLevelRedist.resetWeights(m_stateNew, densevar);

     if (m_hasCoarser && (!s_noEBCF))
      {
        EBAMRReactive* coarPtr = getCoarserLevel();
        coarPtr->m_stateNew.exchange(Interval(0,m_nComp-1));
        m_ebFineToCoarRedist.resetWeights(coarPtr->m_stateNew, densevar);
      }

    if (m_hasFiner && (!s_noEBCF))
      {
        m_ebCoarToFineRedist.resetWeights(m_stateNew, densevar);
        m_ebCoarToFineRedist.resetWeights(m_stateNew, densevar);
      }
    }
}
/****************************/
void EBAMRReactive::
refluxRedistInteraction()
{
  // the flux register must modify the redistribution registers
  CH_TIME("EBAMRReactive::refluxRedistInteraction");
  if (m_hasFiner && (!s_noEBCF))
    {
      Real scale = -1.0/m_dx[0];
      Interval interv(0,m_nComp-1);
      m_ebFluxRegister.incrementRedistRegister(m_ebCoarToFineRedist,
                                                 interv, scale);
  
      m_ebFluxRegister.incrementRedistRegister(m_ebCoarToCoarRedist,
                                                 interv, scale);
    }
}  
/****************************/
Real EBAMRReactive::computeDt()
{
  Real newDt;
  newDt = m_dtNew;

  return newDt;
}
/***************************/
Real
EBAMRReactive::getDt() const
{
  return m_dt;
}
/***************************/
#ifdef CH_USE_HDF5
/***************************/
void EBAMRReactive::writeCheckpointHeader(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
    {
      pout() << "EBAMRReactive::writeCheckpointHeader" << endl;
    }

  //stuff in non-eb checkpoint header is already
  //set in the define function, such as
  // Setup the number of components
  // Setup the component names
  //so i will eliminate the middleman
}
/***************************/
void EBAMRReactive::readCheckpointHeader(HDF5Handle& a_handle)
{
  if (s_verbosity >= 3)
    {
      pout() << "EBAMRReactive::readCheckpointHeader" << endl;
    }

  //stuff in non-eb checkpoint header is already
  //set in the define function, such as
  // Setup the number of components
  // Setup the component names
  // So i will eliminate the middleman.
}
/***************************/
void EBAMRReactive::writeCheckpointLevel(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
    {
      pout() << "EBAMRReactive::writeCheckpointLevel" << endl;
    }

  //reasons for deviations from non-eb stuff
  // the tag buffer size is set by the factory.
  // dx is set in the define function
  // the problem domain is set in the define function

  // Setup the level string
  char levelStr[20];
  sprintf(levelStr,"%d",m_level);
  const std::string label = std::string("level_") + levelStr;

  a_handle.setGroup(label);

  // Setup the level header information
  HDF5HeaderData header;

  header.m_int ["ref_ratio"]       = m_ref_ratio;
  header.m_real["dt"]              = m_dt;

  if (s_verbosity >= 3)
    {
      pout() << "EBAMRReactive::writeCheckpointLevel 2" << endl;
    }

  // Write the header for this level
  header.writeToFile(a_handle);

  // Write the data for this level
  if (s_verbosity >= 3)
    {
      pout() << "EBAMRReactive::writeCheckpointLevel 3" << endl;
    }
  write(a_handle,m_grids);
  if (s_verbosity >= 3)
    {
      pout() << "EBAMRReactive::writeCheckpointLevel 4" << endl;
    }
  write(a_handle,m_stateOld,"dataOld");
  if (s_verbosity >= 3)
    {
      pout() << "EBAMRReactive::writeCheckpointLevel 5" << endl;
    }
  write(a_handle,m_stateNew,"dataNew");
  if (s_verbosity >= 3)
    {
      pout() << "EBAMRReactive::writeCheckpointLevel 6" << endl;
    }
}
/***************************/
/***************************/
void EBAMRReactive::readCheckpointLevel(HDF5Handle& a_handle)
{
  CH_assert(m_isDefined);
  if (s_verbosity >= 3)
    {
      pout() << "EBAMRReactive::readCheckpointLevel" << endl;
    }

  // Setup the level string
  char levelStr[20];
  sprintf(levelStr,"%d",m_level);
  const std::string label = std::string("level_") + levelStr;

  // Read the header for this level
  a_handle.setGroup(label);

  HDF5HeaderData header;
  header.readFromFile(a_handle);

  // Get the refinement ratio
  if (header.m_int.find("ref_ratio") == header.m_int.end())
    {
      MayDay::Error("::readCheckpointLevel: file does not contain ref_ratio");
    }
  m_ref_ratio = header.m_int["ref_ratio"];

  //reasons for deviations from non-eb stuff
  // the tag buffer size is set by the factory.
  // dx is set in the define function
  // the problem domain is set in the define function

  // Get dt
  if (header.m_real.find("dt") == header.m_real.end())
    {
      MayDay::Error("readCheckpointLevel: file does not contain dt");
    }
  m_dt = header.m_real["dt"];

  // Get the grids
  Vector<Box> vboxGrids;
  const int gridStatus = read(a_handle, vboxGrids);
  if (gridStatus != 0)
    {
      MayDay::Error("readCheckpointLevel: file has no grids");
    }

  Vector<int> proc_map;
  if (s_isLoadBalanceSet)
    {
      s_loadBalance(proc_map,vboxGrids, m_domainBox, false);
    }
  else
    {
      LoadBalance(proc_map,vboxGrids);
    }
  broadcast(proc_map, uniqueProc(SerialTask::compute));

  m_grids= DisjointBoxLayout(vboxGrids,proc_map);
  //this keeps the order of the AMRLevel m_level_grids
  //consistent with m_grids
  LayoutIterator lit = m_grids.layoutIterator();
  for (lit.begin(); lit.ok(); ++lit)
    {
      Box b = m_grids.get(lit());
      m_level_grids.push_back(b);
    }

  //the ebisl has to know about the fact that we really have
  //four ghost cells.
  int nGhostEBISL = 6;
  const EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();
  ebisPtr->fillEBISLayout(m_ebisl, m_grids,
                          m_domainBox, nGhostEBISL);
  EBCellFactory factoryNew(m_ebisl);
  //m_nghost is set in define function
  IntVect ivGhost = m_nGhost*IntVect::Unit;
  m_stateNew.define(m_grids,m_nComp, ivGhost, factoryNew);
  m_stateOld.define(m_grids,m_nComp, ivGhost, factoryNew);
  
  m_eblg.define(m_grids, m_problem_domain, m_nGhost, ebisPtr);

  //  Interval vars(0, m_nComp-1);
  //the false says to not redefine the data
  int dataStatusNew = read<EBCellFAB>(a_handle,
                                      m_stateNew,
                                      "dataNew",
                                      m_grids,
                                      Interval(),
                                      false);

  int dataStatusOld = read<EBCellFAB>(a_handle,
                                      m_stateOld,
                                      "dataOld",
                                      m_grids,
                                      Interval(),
                                      false);

  if ((dataStatusNew != 0) || (dataStatusOld != 0))
    {
      MayDay::Error("file does not contain state data");
    }
  // Set up data structures

  levelSetup();
}
/***************************/
void EBAMRReactive::writePlotHeader(HDF5Handle& a_handle) const
{
    if (s_NewPlotFile == 0)
    {
      writePlotHeaderOld(a_handle);
      return;
    }
  Vector<int> refRatios;
  const EBAMRReactive* current = this;
  int nlevs = 0;
  while (current != NULL)
  {
    refRatios.push_back(current->refRatio());
    nlevs++;
    current = (const EBAMRReactive*)(current-> m_finer_level_ptr);
  }

  headerEBFile(a_handle, nlevs, refRatios,
               m_problem_domain.domainBox(),
               m_origin, m_dx, m_aspect,
               m_stateNew.ghostVect());

  writeCellCenteredNames(a_handle, m_stateNames);
  a_handle.setGroup("/Expressions");
  HDF5HeaderData expressions;
  m_ebPatchReactive->expressions(expressions);
  expressions.writeToFile(a_handle);

}
/***************************/
void EBAMRReactive::writePlotHeaderOld(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
    {
      pout() << "EBAMRReactive::writePlotHeader" << endl;
    }

  HDF5HeaderData header;
  // Setup the number of components
  //have to add in a lot of geometric crap.
  // 3 norms + 6 area fracs + 1 distance + 1 volFrac
  // =  11 extra components
  //forces 3d
  int nCons = m_ebPatchReactive->numConserved();
  int nPrim = m_ebPatchReactive->numPrimitives() ;
  int consAndPrim = nCons + nPrim;

  int indexVolFrac = consAndPrim;
  int indexAreaFrac = indexVolFrac+1;
  int indexNormal = indexAreaFrac+ 2*SpaceDim;
  int indexDist = indexNormal+SpaceDim;
  int nCompTotal = indexDist+1;
  CH_assert(nCompTotal == consAndPrim + 3*SpaceDim+2);

  Vector<string> names(nCompTotal);

  for (int i = 0; i < nCons; i++)
    {
      names[i] = m_stateNames[i];
    }
  for (int i = 0; i < nPrim; i++)
    {
      names[nCons + i] = m_primNames[i];
    }

  string volFracName("fraction-0");
  Vector<string> normName(3);
  Vector<string> areaName(6);
  string distName("distance-0");

  normName[0] = "xnormal-0";
  normName[1] = "ynormal-0";
  normName[2] = "znormal-0";

  areaName[0] = "xAreafractionLo-0";
  areaName[1] = "xAreafractionHi-0";
  areaName[2] = "yAreafractionLo-0";
  areaName[3] = "yAreafractionHi-0";
  areaName[4] = "zAreafractionLo-0";
  areaName[5] = "zAreafractionHi-0";

  names[indexVolFrac] = volFracName;

  for (int i=0; i < 2*SpaceDim; i++)
    {
      names[indexAreaFrac+i] = areaName[i];
    }

  for (int i=0; i < SpaceDim; i++)
    {
      names[indexNormal+i] = normName[i];
    }

  names[indexDist] = distName;
 
  //now output this into the hdf5 handle
  header.m_int["num_components"] = nCompTotal;
  // Setup the component names
  char compStr[30];
  for (int comp = 0; comp < nCompTotal; ++comp)
    {
      sprintf(compStr,"component_%d",comp);
      header.m_string[compStr] = names[comp];
    }

  // Write the header to the file
  header.writeToFile(a_handle);

  if (s_verbosity >= 4)
    {
      pout() << header << endl;
    }
}
/***************************/
void EBAMRReactive::writePlotLevel(HDF5Handle& a_handle) const
{
  if (s_NewPlotFile == 0)
    {
      writePlotLevelOld(a_handle);
      return;
    }
  if (s_verbosity >= 3)
    {
      pout() << "EBAMRReactive::writePlotLevel " << m_level<< endl;
    }
  writeCellCentered(a_handle, m_level, &m_stateNew);
}
/***************************/
void EBAMRReactive::writePlotLevelOld(HDF5Handle& a_handle) const
{

  if (s_verbosity >= 3)
    {
      pout() << "EBAMRReactive::writePlotLevel" << endl;
    }
  //fill output data including a bunch of geometric crud
  int nCons = m_ebPatchReactive->numConserved();
  int nPrim = m_ebPatchReactive->numPrimitives() ;
  int consAndPrim = nCons + nPrim;
  int indexVolFrac = consAndPrim;
  int indexAreaFrac = indexVolFrac+1;
  int indexNormal = indexAreaFrac+ 2*SpaceDim;
  int indexDist = indexNormal+SpaceDim;
  int nCompTotal = indexDist+1;
  CH_assert(nCompTotal == consAndPrim + 3*SpaceDim+2);

  Vector<Real> coveredValuesCons(nCons, -10.0);
  Vector<Real> coveredValuesPrim(nPrim, -10.0);

  Vector<Real> coveredValues;
  coveredValues.append(coveredValuesCons);
  coveredValues.append(coveredValuesPrim);

  LevelData<FArrayBox> fabData(m_grids, nCompTotal, IntVect::Zero);

#ifdef CH_MPI
    MPI_Barrier(Chombo_MPI::comm);
#endif
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      const EBISBox& ebisbox = m_ebisl[dit()];
      const Box& grid = m_grids.get(dit());
      const EBCellFAB& consfab = m_stateNew[dit()];
      EBCellFAB primfab(ebisbox, grid, nPrim);
      //cfivs, time and timestep fake and not used here
      Real faket = 1.0;
      IntVectSet emptyivs;
      ParmParse pp;
      int logflag = 1;
      if (pp.contains("logflag"))
        {
          pp.get("logflag", logflag);
        }
      m_ebPatchReactive->setValidBox(grid, ebisbox, emptyivs, faket, faket);
      m_ebPatchReactive->consToPrim(primfab, consfab, grid, logflag);

      FArrayBox& currentFab = fabData[dit()];

      // copy regular data
      currentFab.copy(consfab.getSingleValuedFAB(),0,0,nCons);
      currentFab.copy(primfab.getSingleValuedFAB(),0,nCons,nPrim);

      // set default volume fraction
      currentFab.setVal(1.0,indexVolFrac);

      // set default area fractions
      for (int i=0; i < 2*SpaceDim; i++)
        {
          currentFab.setVal(1.0,indexAreaFrac+i);
        }

      // set default normal
      for (int i=0; i < SpaceDim; i++)
        {
          currentFab.setVal(0.0,indexNormal+i);
        }

      // set default distance of EB from corner
      currentFab.setVal(0.0,indexDist);

      // set special values
      // iterate through the current grid
      // NOTE:  this is probably an inefficient way to do this
      for (BoxIterator bit(grid); bit.ok(); ++bit)
        {
          const IntVect& iv = bit();
          // set special values for covered cells
          if (ebisbox.isCovered(iv))
            {
              for (int icomp = 0; icomp < consAndPrim; icomp++)
                {
                  Real cval = coveredValues[icomp];

                  currentFab(iv,icomp) = cval;
                }
              // volume fraction is zero
              currentFab(iv,indexVolFrac) = 0.0;

              // area fractions are zero
              for (int i=0; i < 2*SpaceDim; i++)
                {
                  currentFab(iv,indexAreaFrac+i) = 0.0;
                }
            }

          // set special values for irregular cells
          if (ebisbox.isIrregular(iv))
            {
              Vector<VolIndex> vofs = ebisbox.getVoFs(iv);
              Real volFrac = ebisbox.volFrac(vofs[0]);
              RealVect normal = ebisbox.normal(vofs[0]);

              // set volume fraction
              currentFab(iv,indexVolFrac) = volFrac;

              // set area fractions--use only the first face you find
              for (int i=0; i < SpaceDim; i++)
                {
                  Vector<FaceIndex> faces;

                  faces = ebisbox.getFaces(vofs[0],i,Side::Lo);
                  if (faces.size() == 0)
                    {
                      currentFab(iv,indexAreaFrac+2*i) = 0.0;
                    }
                  else
                    {
                      currentFab(iv,indexAreaFrac+2*i) =
                        ebisbox.areaFrac(faces[0]);
                    }

                  faces = ebisbox.getFaces(vofs[0],i,Side::Hi);
                  if (faces.size() == 0)
                    {
                      currentFab(iv,indexAreaFrac+2*i+1) = 0.0;
                    }
                  else
                    {
                      currentFab(iv,indexAreaFrac+2*i+1) =
                        ebisbox.areaFrac(faces[0]);
                    }
                }

              // set normal
              for (int i=0; i < SpaceDim; i++)
                {
                  currentFab(iv,indexNormal+i) = normal[i];
                }

              // set distance unless the length of the normal is zero
              Real length = PolyGeom::dot(normal,normal);

              if (length > 0)
                {
                  Real dist = PolyGeom::computeAlpha(volFrac,normal)*m_dx[0];
                  currentFab(iv,indexDist) = -dist;
                }
            } //end if (isIrregular)
        }//end loop over cells
    }//end loop over grids

#ifdef CH_MPI
    MPI_Barrier(Chombo_MPI::comm);
#endif
  // Setup the level string
  char levelStr[20];
  sprintf(levelStr,"%d",m_level);
  const std::string label = std::string("level_") + levelStr;

  a_handle.setGroup(label);

  // Setup the level header information
  HDF5HeaderData header;

  header.m_int ["ref_ratio"]   = m_ref_ratio;
  header.m_real["dx"]          = m_dx[0];
  header.m_real["dt"]          = m_dt;
  header.m_real["time"]        = m_time;
  header.m_box ["prob_domain"] = m_problem_domain.domainBox();

  // Write the header for this level
  header.writeToFile(a_handle);

  if (s_verbosity >= 4)
    {
      pout() << header << endl;
    }

  // Write the data for this level
  write(a_handle,fabData.boxLayout());
  write(a_handle,fabData,"data");
}

#endif
