#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBPATCHREACTIVE_H_
#define _EBPATCHREACTIVE_H_

#include "Box.H"
#include "IntVectSet.H"
#include "Vector.H"
#include "CH_HDF5.H"
#include "EBCellFAB.H"
#include "EBFluxFAB.H"
#include "EBISBox.H"
#include "BaseIVFAB.H"
#include "EBPhysIBC.H"
#include "EBPhysIBCFactory.H"
#include "Stencils.H"
#include "AggStencil.H"
//#include "NamespaceHeader.H"

#include "Stencils.H"
#include "CoveredFaceFAB.H"
#include "BaseIFFAB.H"

#include "NamespaceHeader.H"

///
/**
 */
class EBPatchReactive
{ 
public:
  ///
  /**
   */
  EBPatchReactive();
  
  ///
  /**
   */
  virtual ~EBPatchReactive();

  ///
  /**
   */
  void 
  setGamma(const Real& a_gamma);

  ///
  /**
   */
  Real getGamma() const;

  ///
  /**
   */
  void 
  setSpecHeat(const Real& a_specHeat);

  ///
  /**
   */
  Real getSpecHeat() const;

  ///
  /**
   */
  void 
  setnSpecies(const int& a_nSpec);

  ///
  /**
   */
  int getnSpecies() const;

  ///
  /**
   */
  void
  setEBPhysIBC(const EBPhysIBCFactory& a_bc);

  ///
  /**
   */
  const EBPhysIBC*  getEBPhysIBC() const;

  /// Set parameters for slope computations
  /**
   */
  void setSlopeParameters(bool a_fourthOrderSlopes,
                          bool a_zeroSlopes,
                          bool a_flattening,
                          bool a_useLimiting);
                       

  void
  primitivesAndDivergences(EBCellFAB&          a_nonConsDivF,
                           EBCellFAB&          a_consState,
                           BaseIVFAB<Real>     a_coveredPrimMinu[SpaceDim],
                           BaseIVFAB<Real>     a_coveredPrimPlus[SpaceDim],
                           Vector<VolIndex>    a_coveredFaceMinu[SpaceDim],
                           Vector<VolIndex>    a_coveredFacePlus[SpaceDim],
                           EBFluxFAB&          a_flux,
                           BaseIVFAB<Real>&    a_ebIrregFlux,
                           BaseIVFAB<Real>&    a_nonConservativeDivergence,
                           const EBCellFAB&    a_flattening,
                           const EBCellFAB&    a_source,
                           const Box&          a_box,
                           const IntVectSet&   a_ivsSmall,
                           const DataIndex&    a_dit,
                           bool                a_verbose);

  /**
   */
  void fillDiffCoeffMatrix(EBCellFAB& a_diffCoefMatrix,
                           const EBCellFAB& a_pressure,
                           const EBCellFAB& a_temperature,
                           const EBCellFAB& a_massFrac,
                           const Box&       a_box);

  void fillViscousAndConductiveCoeff(EBCellFAB&       a_aco,
                                     EBCellFAB&       a_mu,
                                     EBCellFAB&       a_lambda,
                                     EBCellFAB&       a_kappa,
                                     const EBCellFAB& a_massFrac,
                                     const EBCellFAB& a_dense,
                                     const EBCellFAB& a_temperature,
                                     const Box&       a_box);

  void getRhoCv(EBCellFAB&       a_rhoCv,
                const EBCellFAB& a_rho,
                const EBCellFAB& a_massFrac,
                const EBCellFAB& a_temperat,
                const Box&       a_box);

  ///
  /**
   */
  void doRZCoords(bool a_doRZCoords);

  ///
  /**
   */
  void artificialViscosity(bool a_useArtificialVisc);

  ///
  /**
   */
  void
  define(const ProblemDomain& a_domain,
         const RealVect& a_dx, bool a_useAgg = false);


  ///deprecated interface
  void
  define(const ProblemDomain& a_domain,
         const Real& a_dx)
  {
    m_useAgg = false;
    RealVect dxVect = a_dx*RealVect::Unit;
    define(a_domain, dxVect);
  }

  ///
  /**
     Returns number of components for conserved variables.
  */
  int numConserved() const;

  ///
  /**
     Returns number of components for flux variables.
  */
  int numFluxes() const;

  ///
  /**
     Returns number of components for primitive variables.
  */
  int numPrimitives() const;

  ///
  /**
   */
  int numSlopes() const;

  ///
  /**
     Return the names of the variables.  A default
     implementation is provided that puts in generic names.
  */
  Vector<string> stateNames();

  ///
  /**
     Return the names of the variables.  A default
     implementation is provided that puts in generic names.
  */
  Vector<string> primNames();

  ///
  /**
     needs coarse-fine IVS to know where to drop order for interpolation
     virtual in case you need to add anything to definition
   */
  virtual void
  setValidBox(const Box& a_validBox,
              const EBISBox& a_ebisBox,
              const IntVectSet& a_coarseFineIVS,
              const Real& a_time,
              const Real& a_dt);


  //for testing
  void getArgBox(Box a_argBox[SpaceDim]);

  //stuff for aggregated slopes
  void setSlopeStuff();

 ///
  /**
   */
  Real
  getMaxWaveSpeed(const EBCellFAB& a_consState,
                  const Box& a_box);


  ///
  /**
     Compute the primitive state given the conserved state.
     W = W(U).
  */
  void
  consToPrim(EBCellFAB&       a_primState,
             const EBCellFAB& a_conState,
             const Box&       a_box,
             int              a_logflag,
             bool             a_verbose=false);

  void
  consToPrim(BaseIVFAB<Real>&       a_primState,
             const BaseIVFAB<Real>& a_conState,
             const IntVectSet&      a_ivs);

  void
  consToPrim(BaseIVFAB<Real>&  a_primState,
             const EBCellFAB&  a_conState,
             const IntVectSet& a_ivs);

  static IntVect s_debugIV;
  static int     s_whichLev;

    ///
  /**
   */
  void
  setCoveredConsVals(EBCellFAB& a_consState);
 
    void
  computeFlattening(EBCellFAB&       a_flattening,
                    const EBCellFAB& a_primState,
                    const Box&       a_box);

  ///
  /**
     Returns the interval of component indices in the primitive variable
     EBCellFAB for the velocities.
  */
  Interval velocityInterval() const;

  ///
  /**
     Returns the component index for the density.
  */
  virtual int densityIndex() const;

  ///
  /**
     Returns the component index for the pressure. Called only if flattening is used.
  */
  int pressureIndex() const;
  ///
  /**
     Returns the component index for the bulk modulus, used as a
     normalization to measure shock strength in flattening.
     Called only if flattening is used.
  */
  int bulkModulusIndex() const;

  virtual int densityIndexC() const; // conserved vector

  ///
  /**
     Returns the component index for the energy in C vector.
  */
  virtual int energyIndexC() const;

  int temperatureIndex() const;

  ///
  //Returns species 1 mass frac and dense index
  int spec1MassFracIndex() const;
  int spec1DenseIndex() const;

  Interval momentumInterval() const;

  Interval speciesMassFracInterval() const;
  Interval speciesDenseInterval() const;

 //RZ func 
  void
  setSource(EBCellFAB&       a_source,
              const EBCellFAB& a_consState,
              const Box&       a_box);

  void
  floorPrimitives(EBCellFAB&  a_primState,
                  const Box&  a_box);

  void
  floorConserved(EBCellFAB&  a_consState,
                 const Box&  a_box);


  void
  floorConserved(BaseIVFAB<Real>&  a_consState,
                 const IntVectSet&  a_set);

  void
  floorPrimitives(BaseIVFAB<Real>&  a_primState,
                  const IntVectSet&  a_set);

   ///
  /**
     Return true if the application is using flattening.
  */
  bool usesFlattening() const;

  ///
  /**
     Return true if the application is using artificial viscosity.
  */
  bool usesArtificialViscosity() const;

  ///
  /**
     Return true if you are using fourth-order slopes.
     Return false if you are using second-order slopes.
  */
  bool usesFourthOrderSlopes() const;

  bool usesZeroSlopes() const; 

  bool isDefined() const;

  ///
  /**
   */
  void
  computeFluxes(EBFluxFAB&          a_flux,
                BaseIVFAB<Real>     a_coveredFluxLo[SpaceDim],
                BaseIVFAB<Real>     a_coveredFluxHi[SpaceDim],
                Vector<VolIndex>    a_coveredFaceLo[SpaceDim],
                Vector<VolIndex>    a_coveredFaceHi[SpaceDim],
                EBCellFAB&          a_primState,
                EBCellFAB           a_slopesPrim[SpaceDim],
                EBCellFAB           a_slopesSeco[SpaceDim],
                const EBCellFAB&    a_flattening,
                const EBCellFAB&    a_consState,
                const EBCellFAB&    a_source,
                const Box&          a_box,
                const DataIndex&    a_dit,
                bool                a_verbose);

  ///
  /**
   */
  void
  extrapolatePrim2D(EBCellFAB           a_primMinu[SpaceDim],
                    EBCellFAB           a_primPlus[SpaceDim],
                    EBCellFAB&          a_primState,
                    EBCellFAB           a_slopesPrim[SpaceDim],
                    EBCellFAB           a_slopesSeco[SpaceDim],
                    const EBCellFAB&    a_flattening,
                    const EBCellFAB&    a_consState,
                    const EBCellFAB&    a_source,
                    const Box&          a_box,
                    const DataIndex&    a_dit,
                    bool                a_verbose);

  ///
  /**
   */
  void
  doNormalDerivativeExtr2D(EBCellFAB              a_primMinu[SpaceDim],
                           EBCellFAB              a_primPlus[SpaceDim],
                           EBFaceFAB              a_fluxOne[SpaceDim],
                           BaseIVFAB<Real>        a_coveredFluxNormMinu[SpaceDim],
                           BaseIVFAB<Real>        a_coveredFluxNormPlus[SpaceDim],
                           Vector<VolIndex>       a_coveredFaceNormMinu[SpaceDim],
                           Vector<VolIndex>       a_coveredFaceNormPlus[SpaceDim],
                           EBCellFAB              a_slopesPrim[SpaceDim],
                           EBCellFAB              a_slopesSeco[SpaceDim],
                           const EBCellFAB&       a_flattening,
                           const EBCellFAB&       a_primState,
                           const EBCellFAB&       a_source,
                           const DataIndex&       a_dit,
                           const Box&             a_box);

  void
  extrapolatePrim3D(EBCellFAB           a_primMinu[SpaceDim],
                    EBCellFAB           a_primPlus[SpaceDim],
                    EBCellFAB&          a_primState,
                    EBCellFAB           a_slopesPrim[SpaceDim],
                    EBCellFAB           a_slopesSeco[SpaceDim],
                    const EBCellFAB&    a_flattening,
                    const EBCellFAB&    a_consState,
                    const EBCellFAB&    a_source,
                    const Box&          a_box,
                    const DataIndex&    a_dit,
                    bool                a_verbose);

  void
  doNormalDerivativeExtr3D(EBCellFAB              a_primMinu[SpaceDim],
                           EBCellFAB              a_primPlus[SpaceDim],
                           EBFaceFAB              a_fluxOne[SpaceDim],
                           BaseIVFAB<Real>        a_coveredFluxNormMinu[SpaceDim],
                           BaseIVFAB<Real>        a_coveredFluxNormPlus[SpaceDim],
                           Vector<VolIndex>       a_coveredFaceNormMinu[SpaceDim],
                           Vector<VolIndex>       a_coveredFaceNormPlus[SpaceDim],
                           EBCellFAB              a_slopesPrim[SpaceDim],
                           EBCellFAB              a_slopesSeco[SpaceDim],
                           const EBCellFAB&       a_flattening,
                           const EBCellFAB&       a_primState,
                           const EBCellFAB&       a_source,
                           const DataIndex&       a_dit,
                           const Box&             a_box);

  ///
  /**
     Compute the limited slope /a_dq/ of the primitive variables
     /a_q/ for the components in the interval /a_interval/,
     Calls applyLimiter.
  */
  void
  slope(EBCellFAB&       a_slopePrim,
        EBCellFAB&       a_slopeNLim,
        const EBCellFAB& a_primState,
        const EBCellFAB& a_flattening,
        const int&       a_dir,
        const Box&       a_box,
        bool a_doAggregated = false);

  void 
  doSecondOrderSlopes(EBCellFAB&       a_delta2W,
                      EBCellFAB&       a_deltaWL,
                      EBCellFAB&       a_deltaWR,
                      EBCellFAB&       a_deltaWC,
                      const EBCellFAB& a_primState,
                      const int&       a_dir,
                      const Box&       a_box,
                      bool a_doAggregated = false);

  virtual void
  aggIrregSecondOrderSlopes(EBCellFAB&       a_delta2W,
                            EBCellFAB&       a_deltaWL,
                            EBCellFAB&       a_deltaWR,
                            EBCellFAB&       a_deltaWC,
                            const EBCellFAB& a_primState,
                            const int&       a_dir,
                            const Box&       a_box);


  virtual void
  irregSecondOrderSlopes(EBCellFAB&       a_delta2W,
                         EBCellFAB&       a_deltaWL,
                         EBCellFAB&       a_deltaWR,
                         EBCellFAB&       a_deltaWC,
                         const EBCellFAB& a_primState,
                         const int&       a_dir,
                         const Box&       a_box);

  void
  pointGetSlopes(Real&            a_dql,
                 Real&            a_dqr,
                 Real&            a_dqc,
                 bool&            a_hasFacesLeft,
                 bool&            a_hasFacesRigh,
                 const VolIndex&  a_vof,
                 const EBCellFAB& a_primState,
                 const int&       a_dir,
                 const int&       a_ivar,
                 const bool&      a_verbose);

  void
  doFourthOrderSlopes(EBCellFAB& a_delta4W,
                      EBCellFAB& a_deltaWC,
                      const EBCellFAB& a_delta2W,
                      const EBCellFAB& a_deltaWL,
                      const EBCellFAB& a_deltaWR,
                      const EBCellFAB& a_primState,
                      const int& a_dir,
                      const Box& a_box);

  ///
  /**
     Given left and right one-sided undivided differences /a_dql,a_dqr/,
     apply van Leer limiter $vL$ defined in section
     to a_dq.
  */
  void
  applyLimiter(EBCellFAB&       a_slopePrim,
               const EBCellFAB& a_slopePrimLeft,
               const EBCellFAB& a_slopePrimRigh,
               const int&       a_dir,
               const Box&       a_box);

  virtual void
  normalPred(EBCellFAB&       a_primLo,
             EBCellFAB&       a_primHi,
             const EBCellFAB& a_primState,
             const EBCellFAB& a_slopePrim,
             const Real&      a_scale,
             const int&       a_dir,
             const Box&       a_box);

  ///
  /**
   */
  void
  primToCons(EBCellFAB&       a_consState,
             const EBCellFAB& a_primState,
             const Box&       a_box);

  ///
  /**
   */
  void
  primToCons(BaseIVFAB<Real>&       a_consState,
             const BaseIVFAB<Real>& a_primState,
             const IntVectSet&      a_ivs);

  ///

  /**
     Given input left and right states, compute a suitably-upwinded
     flux (e.g. by solving a Riemann problem), as in
  */
  void
  riemann(EBFaceFAB&       a_flux,
          const EBCellFAB& a_primLeft,
          const EBCellFAB& a_primRight,
          const int&       a_dir,
          const Box&       a_box);

  ///
  /**
     Given input left and right states, compute a suitably-upwinded
     flux (e.g. by solving a Riemann problem).
  */
  void
  riemann(BaseIVFAB<Real>&        a_coveredFlux,
          const BaseIVFAB<Real>&  a_extendedState,
          const EBCellFAB&        a_primState,
          const Vector<VolIndex>& a_region,
          const int&              a_dir,
          const Side::LoHiSide&   a_sd,
          const Box&       a_box);

  void
  do111coupling(EBFaceFAB              a_fluxTwoVec[SpaceDim][SpaceDim],
                BaseIVFAB<Real>        a_coveredFlux3DMinu[SpaceDim][SpaceDim],
                BaseIVFAB<Real>        a_coveredFlux3DPlus[SpaceDim][SpaceDim],
                const EBCellFAB        a_primMinu[SpaceDim],
                const EBCellFAB        a_primPlus[SpaceDim],
                const BaseIVFAB<Real>  a_coveredFluxNormMinu[SpaceDim],
                const BaseIVFAB<Real>  a_coveredFluxNormPlus[SpaceDim],
                const Vector<VolIndex> a_coveredFaceNormMinu[SpaceDim],
                const Vector<VolIndex> a_coveredFaceNormPlus[SpaceDim],
                const EBFaceFAB        a_fluxOne[SpaceDim],
                const EBCellFAB&       a_primState,
                const EBCellFAB        a_slopesPrim[SpaceDim],
                const EBCellFAB        a_slopesSeco[SpaceDim],
                const DataIndex&       a_dit,
                const Box&             a_box);

  Real bilinearFunc(const Real  a_WVal[2][2],
                    const Real& a_xd1,
                    const Real& a_xd2);


   void integrateReactiveSource(EBCellFAB&  a_consState,
                                const Box&  a_box,
                                const Real& a_dt);

  ///
  /**
     Returns value of artificial viscosity. Called only if artificial
     viscosity is being used.
  */
  Real artificialViscosityCoefficient() const;

  //these are the objects I need to make updatePrim less of a dog
  struct
  {
    int     m_offset;
    bool    m_multiValued;
  } typedef pointerOffset_t;

  struct
  {
    pointerOffset_t m_vofOffset;
  } typedef updateStencil_t;

  //gets offsets
  void fillUpdateStencil(EBPatchReactive::updateStencil_t& a_sten,
                         const VolIndex&                   a_vof); 

#ifdef CH_USE_HDF5
  virtual void expressions(HDF5HeaderData& a_expressions);
#endif

  //debugging hooks.  not for the faint of heart
  static void setVerbose(bool a_verbose);
  static void setCurLevel(int a_curLevel);
  static void setCurComp(int a_curComp);
  static void setDoingVel(int a_yesorno);
  static void setDoingAdvVel(int a_yesorno);
  static int getDoingAdvVel();
  static int getDoingVel();
  static int getCurLevel();
  static int getCurComp();
  static bool getVerbose();
  static Real getMaxWaveSpeed();
  static void setMaxWaveSpeed(Real a_maxWaveSpeedIV);
  static IntVect getMaxWaveSpeedIV();
  static void    setMaxWaveSpeedIV(const IntVect& a_maxWaveSpeedIV);

  EBCellFAB& getPrimState()
  {
    return m_primState;
  }
  BaseIVFAB<Real>* getCoveredFluxPlus()
  {
    return m_coveredFluxPlusG4;
  }
  BaseIVFAB<Real>* getCoveredFluxMinu()
  {
    return m_coveredFluxMinuG4;
  }

  Vector<VolIndex>* getCoveredFaceMinu()
  {
    return m_coveredFaceMinuG4;
  }

  Vector<VolIndex>* getCoveredFacePlus()
  {
    return m_coveredFacePlusG4;
  }

  ///
  /**
     set to true if the source you will provide
     is in conservative variables.   Default is false
  */
  static void useConservativeSource(bool a_conservativeSource)
  {
    s_conservativeSource = a_conservativeSource;
  }

  ///non-virtual stuff
  void
  incrementWithSource(EBCellFAB&       a_primState,
                      const EBCellFAB& a_source,
                      const Real&      a_scale,
                      const Box&       a_box) ;

  void
  extrapToCoveredFaces(BaseIVFAB<Real>&        a_extendedPrim,
                       const EBCellFAB&        a_primMinu,
                       const EBCellFAB&        a_primPlus,
                       const EBCellFAB&        a_primState,
                       const Vector<VolIndex>& a_coveredFaces,
                       const int&              a_faceDir,
                       const Side::LoHiSide&   a_sd,
                       const Box&              a_box);

  void
  pointExtrapToCovered2D(Vector<Real>&           a_extrapVal,
                         const EBCellFAB&        a_primMinu,
                         const EBCellFAB&        a_primPlus,
                         const EBCellFAB&        a_primState,
                         const int&              a_faceDir,
                         const VolIndex&         a_vof,
                         const RealVect&         a_normal,
                         const Side::LoHiSide&   a_sd,
                         const int&              a_numPrim);

  ///
  void
  pointExtrapToCovered3D(Vector<Real>&           a_extrapVal,
                         const EBCellFAB&        a_primMinu,
                         const EBCellFAB&        a_primPlus,
                         const EBCellFAB&        a_primState,
                         const int&              a_faceDir,
                         const VolIndex&         a_vof,
                         const RealVect&         a_normal,
                         const Side::LoHiSide&   a_sd,
                         const int&              a_numPrim);


  void
  coveredExtrapSlopes(Real&            a_dqc,
                      const VolIndex&  a_vof,
                      const EBCellFAB& a_primState,
                      const int&       a_dir,
                      const int&       a_ivar);

  void
  finalExtrap2D(EBCellFAB              a_primMinu[SpaceDim],
                EBCellFAB              a_primPlus[SpaceDim],
                const BaseIVFAB<Real>  a_coveredFluxNormMinu[SpaceDim],
                const BaseIVFAB<Real>  a_coveredFluxNormPlus[SpaceDim],
                const Vector<VolIndex> a_coveredFaceNormMinu[SpaceDim],
                const Vector<VolIndex> a_coveredFaceNormPlus[SpaceDim],
                const EBFaceFAB        a_fluxOne[SpaceDim],
                const EBCellFAB&       a_primState,
                const EBCellFAB        a_slopesPrim[SpaceDim],
                const EBCellFAB        a_slopesSeco[SpaceDim],
                const Box&             a_box);

  void
  finalExtrap3D(EBCellFAB              a_primMinu[SpaceDim],
                EBCellFAB              a_primPlus[SpaceDim],
                const BaseIVFAB<Real>  a_coveredFlux3DMinu[SpaceDim][SpaceDim],
                const BaseIVFAB<Real>  a_coveredFlux3DPlus[SpaceDim][SpaceDim],
                const EBFaceFAB        a_fluxTwoVec[SpaceDim][SpaceDim],
                const EBCellFAB&       a_primState,
                const EBCellFAB        a_slopesPrim[SpaceDim],
                const EBCellFAB        a_slopesSeco[SpaceDim],
                const Box&             a_box); 

  void
  updatePrim(EBCellFAB&              a_primMinu,
             EBCellFAB&              a_primPlus,
             const EBFaceFAB&        a_flux,
             const BaseIVFAB<Real>&  a_coveredFluxMinu,
             const BaseIVFAB<Real>&  a_coveredFluxPlus,
             const Vector<VolIndex>& a_coveredFaceMinu,
             const Vector<VolIndex>& a_coveredFacePlus,
             const int&              a_dir,
             const Box&              a_box,
             const Real&             a_scale) ;

  void
  updateCons(EBCellFAB&              a_consState,
             const EBFaceFAB&        a_flux,
             const BaseIVFAB<Real>&  a_coveredFluxMinu,
             const BaseIVFAB<Real>&  a_coveredFluxPlus,
             const Vector<VolIndex>& a_coveredFaceMinu,
             const Vector<VolIndex>& a_coveredFacePlus,
             const int&              a_dir,
             const Box&              a_box,
             const Real&             a_scale);

  //defines and fills cache
  void   cacheEBCF(Vector<Vector<Real> >&       a_cache,
                   const EBCellFAB&             a_input);

  //puts values of cache into output
  void uncacheEBCF(EBCellFAB&                   a_output,
                   const Vector<Vector<Real> >& a_cache);

  void
  computeCoveredFaces(Vector<VolIndex>&     a_coveredFace,
                      IntVectSet&           a_coveredIVS,
                      IntVectSet&           a_irregIVS,
                      const int&            a_idir,
                      const Side::LoHiSide& a_sd,
                      const Box&            a_region);

  void
  getFaceDivergence(EBFluxFAB&        a_openDivU,
                    const EBCellFAB&  a_primState,
                    const EBCellFAB   a_slopePrim[SpaceDim],
                    const Box&        a_box,
                    const IntVectSet& a_ivsIrreg);

  void
  applyArtificialViscosity(EBFluxFAB&        a_openFlux,
                           BaseIVFAB<Real>  a_coveredFluxMinu[SpaceDim],
                           BaseIVFAB<Real>  a_coveredFluxPlus[SpaceDim],
                           const Vector<VolIndex> a_coveredFaceMinu[SpaceDim],
                           const Vector<VolIndex> a_coveredFacePlus[SpaceDim],
                           const EBCellFAB&  a_consState,
                           const EBFluxFAB&  a_divVel,
                           const Box&        a_box,
                           const IntVectSet& a_ivsIrreg);

 void 
  nonconservativeDivergence(EBCellFAB&             a_divF,
                            const EBFluxFAB&       a_flux,
                            const BaseIVFAB<Real>  a_coveredFluxMinu[SpaceDim],
                            const BaseIVFAB<Real>  a_coveredFluxPlus[SpaceDim],
                            const Vector<VolIndex> a_coveredFaceMinu[SpaceDim],
                            const Vector<VolIndex> a_coveredFacePlus[SpaceDim],
                            const Box&             a_box);

 void
  computeEBIrregFlux(BaseIVFAB<Real>&  a_ebIrregFlux,
                     const EBCellFAB&  a_primState,
                     const EBCellFAB   a_slopePrim[SpaceDim],
                     const IntVectSet& a_irregIVS,
                     const EBCellFAB&  a_flux);

  ///rz func.
  virtual void
  assembleFluxReg(EBFaceFAB&       a_fluxRegFlux,
                  const EBFaceFAB& a_godunovFlux,
                  const int&       a_idir,
                  const Box&       a_cellBox);


  ///rz func.
  virtual void
  assembleFluxIrr(BaseIFFAB<Real>&       a_fluxRegFlux,
                  const BaseIFFAB<Real>& a_godunovFlux,
                  const int&             a_idir,
                  const Box&             a_cellBox,
                  const IntVectSet&      a_set);


  virtual void
  interpolateFluxToCentroids(BaseIFFAB<Real>              a_centroidFlux[SpaceDim],
                             const BaseIFFAB<Real>* const a_fluxInterpolant[SpaceDim],
                             const IntVectSet&            a_irregIVS);


  virtual void
  hybridDivergence(EBCellFAB&             a_hybridDiv,
                   EBCellFAB&             a_consState,
                   BaseIVFAB<Real>&       a_massDiff,
                   const BaseIFFAB<Real>  a_centroidFlux[SpaceDim],
                   const BaseIVFAB<Real>& a_ebIrregFlux,
                   const BaseIVFAB<Real>& a_nonConservativeDivergence,
                   const Box&             a_box,
                   const IntVectSet&      a_ivs);

  void
  consUndividedDivergence(BaseIVFAB<Real>&       a_divF,
                          const BaseIFFAB<Real>  a_centroidFlux[SpaceDim],
                          const BaseIVFAB<Real>& a_ebIrregFlux,
                          const IntVectSet&      a_ivs);

  void
  finalUpdate(EBCellFAB&             a_consState,
              BaseIVFAB<Real>&       a_massDiff,
              const BaseIVFAB<Real>& a_nonConsDivF,
              const BaseIVFAB<Real>& a_conservDivF,
              const IntVectSet&      a_ivs);

  void
  setRZSource(EBCellFAB&       a_source,
              const EBCellFAB& a_consState,
              const Box&       a_box);

protected:

  //the wild variety of crap I need to stencilize slopes
  Box m_modBoxOpen[SpaceDim];

  static bool s_conservativeSource;
  //stuff for debugging output
  static int  s_curLevel;
  static int  s_curComp;
  static int  s_doingVel;
  static int  s_doingAdvVel;
  static bool s_verbose;
  static Real s_maxWaveSpeed;
  static IntVect s_maxWaveSpeedIV;
  //set in define()
  ProblemDomain m_domain;
  RealVect m_dx;
  Real m_dxScale;
  bool m_useAgg;
  BaseIFFAB<FaceStencil> m_interpStencils[SpaceDim];
  Box m_validBox;

  EBISBox m_ebisBox;
  Real m_time;
  Real m_dt;

  //bool m_useAgg;

  bool m_isDefined;
  bool m_isBCSet;
  bool m_isBoxSet;
  bool m_isSlopeSet;
  bool m_isArtViscSet;
  bool m_useFourthOrderSlopes;
  bool m_useZeroSlopes;
  bool m_useFlattening;
  bool m_useLimiting;
  bool m_useArtificialVisc;
  Vector<VolIndex> m_irregVoFs;

  bool m_isnSpeciesSet;

  bool m_isBCNull;

  Vector<updateStencil_t>  m_updateStencil;
  //set by factory
  EBPhysIBC* m_bc;

  //all stuff sent to UpdatePrim, etc needs to be built with this
  //set of info if any of these stenciling optimizations are to have any hope
  Box                     m_validBoxG4;
  IntVectSet       m_coveredSetsPlusG4[SpaceDim];
  IntVectSet       m_coveredSetsMinuG4[SpaceDim];
  Vector<VolIndex> m_coveredFacePlusG4[SpaceDim];
  Vector<VolIndex> m_coveredFaceMinuG4[SpaceDim];

  BaseIVFAB<Real>  m_extendStatePlusG4[SpaceDim];
  BaseIVFAB<Real>  m_extendStateMinuG4[SpaceDim];
  BaseIVFAB<Real>  m_coveredFluxPlusG4[SpaceDim];
  BaseIVFAB<Real>  m_coveredFluxMinuG4[SpaceDim];

  BaseIVFAB<Real>  m_extendStatePlus3D[SpaceDim][SpaceDim];
  BaseIVFAB<Real>  m_extendStateMinu3D[SpaceDim][SpaceDim];
  BaseIVFAB<Real>  m_coveredFluxPlus3D[SpaceDim][SpaceDim];
  BaseIVFAB<Real>  m_coveredFluxMinu3D[SpaceDim][SpaceDim];

  EBCellFAB        m_primPlus[SpaceDim];
  EBCellFAB        m_primMinu[SpaceDim];
  EBCellFAB        m_primState;
  EBCellFAB        m_primMinuTemp;
  EBCellFAB        m_primPlusTemp;

  EBFluxFAB        m_primGdnv;
  EBFaceFAB        m_fluxOne[SpaceDim];
  EBFaceFAB        m_fluxTwo[SpaceDim][SpaceDim];

  BaseIVFAB<Real>  m_coveredFluxNormMinu[SpaceDim];
  BaseIVFAB<Real>  m_coveredFluxNormPlus[SpaceDim];
  BaseIVFAB<Real>  m_extendStateNormMinu[SpaceDim];
  BaseIVFAB<Real>  m_extendStateNormPlus[SpaceDim];

  RefCountedPtr< AggStencil<EBCellFAB, EBCellFAB> > m_slopStenLo[SpaceDim];
  RefCountedPtr< AggStencil<EBCellFAB, EBCellFAB> > m_slopStenHi[SpaceDim];
  Box m_entireBox[SpaceDim];

 struct
  {
    size_t offset;
    int    dataID;
  } typedef access_t;

  struct
  {
    //slopes have all the same box so they will have the same offset and type
    access_t slop_access; //offsets for slopes
    bool hasLo, hasHi;
  } typedef slop_logic_t;

  Vector<slop_logic_t> m_slopVec[SpaceDim];
   
private:

  bool m_isGammaSet;
  Real m_gamma;
  bool m_isSpecHeatSet;
  Real m_specHeat;
  bool m_doRZCoords;
  bool m_isnSpecSet;
  int m_nSpec;

private:
  //disallowed for all the usual reasons
  void operator=(const EBPatchReactive& a_input)
  {
    MayDay::Error("invalid operator");
  }
  EBPatchReactive(const EBPatchReactive& a_input)
  {
    MayDay::Error("invalid operator");
  }

};  
#include "NamespaceFooter.H"
#endif
