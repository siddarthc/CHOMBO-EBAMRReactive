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

#include "BoxIterator.H"
#include "VoFIterator.H"
#include "Stencils.H"
#include "EBArith.H"

#include "MixedPoissonEBBC.H"
#include "CH_Timer.H"
#include "NamespaceHeader.H"

// This file only works for mixed BC applied along an horizontal EB along a grid line
// medium 1 (e.g. gas) over the EB and medium 2 (e.g. dielectric) under the EB
// BC identifiers: 0=Neumann, 1=Dirichlet, 10=Mixed

int MixedPoissonEBBC::s_velComp = 0;

bool MixedPoissonEBBC::s_areaFracWeighted = false;
Real MixedPoissonEBBC::s_timePrint = 0;

MixedPoissonEBBC::MixedPoissonEBBC()
{
  m_dataBased = false;
}

MixedPoissonEBBC::MixedPoissonEBBC(const ProblemDomain& a_domain,
				   const EBISLayout&    a_layout,
				   const RealVect&      a_dx,
				   const IntVect*       a_ghostCellsPhi /*=0*/,
				   const IntVect*       a_ghostCellsRhs /*=0*/)
{
  construct(a_domain, a_layout, a_dx, a_ghostCellsPhi, a_ghostCellsRhs);
}

void
MixedPoissonEBBC::construct(const ProblemDomain& a_domain,
			    const EBISLayout&    a_layout,
			    const RealVect&      a_dx,
			    const IntVect*       a_ghostCellsPhi /*=0*/,
			    const IntVect*       a_ghostCellsRhs /*=0*/)
{
  m_ghostCellsPhi = (*a_ghostCellsPhi) ;
  m_ghostCellsRHS = (*a_ghostCellsRhs) ;

  CH_assert( bool(a_ghostCellsPhi) == bool(a_ghostCellsRhs) ); // !xor

  m_dvalue = 12345.6789;
  m_nvalue = 12345.6789;
  m_func = RefCountedPtr<BaseMixBCValue>();
  m_ebType = 1;
  m_order = 1;
  m_dymin = 12345.6789*RealVect::Unit;
  m_epsLo = 12345.6789;
  m_epsHi = 12345.6789;

  m_onlyHomogeneous = true;
  m_domain = a_domain;
  m_layout = a_layout;

  m_dx = a_dx;
  m_dataBased = false;
  m_isDefined = false;
}

MixedPoissonEBBC::~MixedPoissonEBBC()
{
}

void MixedPoissonEBBC::setArguments(Real 							a_dvalue, 
				    Real 							a_nvalue, 
				    RefCountedPtr<BaseMixBCValue> 	a_func,
				    int								a_ebType,
				    RealVect						a_dymin,
				    Real							a_epsLo,
				    Real							a_epsHi)
{
  m_dvalue = a_dvalue;
  m_nvalue = a_nvalue;
  m_func = a_func;
  m_ebType = a_ebType;
  m_dymin = a_dymin;
  m_epsLo = a_epsLo;
  m_epsHi = a_epsHi;
  m_onlyHomogeneous = false;
}

void MixedPoissonEBBC::setOrder(int a_order)
{
  CH_assert(a_order >= 1 && a_order <= 2);

  if (m_order != a_order)
    {
      m_isDefined = false;
    }
  m_order = a_order;
}

//------------------------------------------------------------------------------------------------------------------------------------

void MixedPoissonEBBC::define(const LayoutData<IntVectSet>& a_cfivs,
                              const Real&                   a_factor)
{
  bool isDirBC;
  if ((m_ebType == 1) || (m_ebType == 10))
    {	
      const DisjointBoxLayout& dbl = m_layout.getDisjointLayout();

      LayoutData<VoFIterator > vofItIrreg;
      vofItIrreg.define(dbl); // vofiterator cache

      m_fluxStencil.define(dbl);
      m_fluxWeight.define(dbl);
      m_custStencil.define(dbl);
      m_custWeight.define(dbl);
      m_rhos.define(dbl);

      //make the Dirichlet stencils
      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
	{
	  const Box& curBox = dbl[dit()];
	  const EBISBox& curEBISBox = m_layout[dit()];
	  const EBGraph& curEBGraph = curEBISBox.getEBGraph();
	  const IntVectSet& cfivsThisBox = a_cfivs[dit()];

	  IntVectSet notRegular;
	  int nComps = 1;

	  notRegular |= curEBISBox.getIrregIVS  (curBox);
	  notRegular |= curEBISBox.getMultiCells(curBox);

	  vofItIrreg[dit()].define(notRegular,curEBISBox.getEBGraph());

	  BaseIVFAB<VoFStencil>& curStencilBaseIVFAB = m_fluxStencil[dit()];     
	  BaseIVFAB<Real>&       curWeightBaseIVFAB  = m_fluxWeight[dit()];
 	  BaseIVFAB<VoFStencil>& custStencilBaseIVFAB = m_custStencil[dit()];
 	  BaseIVFAB<Real>&       custWeightBaseIVFAB  = m_custWeight[dit()];

	  curStencilBaseIVFAB.define(notRegular,curEBGraph,nComps);     
	  curWeightBaseIVFAB.define(notRegular,curEBGraph,nComps);
	  custStencilBaseIVFAB.define(notRegular,curEBGraph,nComps);
	  custWeightBaseIVFAB.define(notRegular,curEBGraph,nComps);
	  m_rhos[dit()].define(notRegular,curEBGraph,nComps);
	  m_rhos[dit()].setVal(0e0);

	  for (VoFIterator vofit(notRegular,curEBGraph); vofit.ok(); ++vofit)
	    {
	      const VolIndex& vof = vofit();
	      VoFStencil& curStencil = curStencilBaseIVFAB(vof,0);
	      Real&  curWeight  = curWeightBaseIVFAB(vof,0);
          
	      Real areaFrac = curEBISBox.bndryArea(vof);
	      Real&  custWeight  = custWeightBaseIVFAB(vof,0);
	      VoFStencil& custStencil = custStencilBaseIVFAB(vof,0);

	      if (m_order == 1)
		{
		  getFirstOrderStencil(curStencil,curWeight,vof,curEBISBox,m_dx);
		}
	      else if (m_order == 2)
		{
		  getSecondOrderStencil(curStencil,curWeight,vof,curEBISBox,m_dx,cfivsThisBox);
		}
	      else
		{
		  MayDay::Error("MixedPoissonEBBC::define stencil order not 1 or 2");
		}
	      //Need to magnify weight in fluxStencil with areafrac*factor, factor = 1/dx;
	      //Pass the factor in because m_dx[0] in here may not be same as in EBAMRPoissonOp
	      curStencil *= areaFrac*a_factor;
	      if (s_areaFracWeighted)
		{
		  curStencil *= curEBISBox.areaFracScaling(vof);
		}
	      custStencil = curStencil;	

	    }   
	}

      // Add correction to stencil for surface potential at Dielectric/gas boundary if required,
      defineCustom(a_cfivs, a_factor);
      m_isDefined = true;

    }

  else
    {
      //no flux stencil for Neumann
      ;
    }
}

//--------------------------------------------------------------------------------------------------------------
// called by define above, accounts for Dielectric/Air EB BC, only at the highest level of refinement
// Assumes that the EB is a horizontal line separating 2 rows of irregular cells

void MixedPoissonEBBC::defineCustom(const LayoutData<IntVectSet>& a_cfivs,
				    const Real&                   a_factor)
{
  bool m_finest = (m_dx[1] == m_dymin[1]); 
  pout() << "MixedPoissonEBBC::defineCustom, bool m_finest " << m_finest <<  ", " << m_dx[1] << ", " << m_dymin[1] << endl;	
  if (m_finest && (m_ebType == 10))
    {	

      int comp = 0;
      bool isDirBC;
      const DisjointBoxLayout& dbl = m_layout.getDisjointLayout();

      LayoutData<VoFIterator > vofItIrreg;
      vofItIrreg.define(dbl); // vofiterator cache

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
	{
	  const Box& curBox = dbl[dit()];
	  const EBISBox& curEBISBox = m_layout[dit()];
	  const EBGraph& curEBGraph = curEBISBox.getEBGraph();
	  const IntVectSet& cfivsThisBox = a_cfivs[dit()];

	  IntVectSet notRegular;
	  int nComps = 1;

	  notRegular |= curEBISBox.getIrregIVS  (curBox);
	  //assuming no multi-cells here
	  //notRegular |= curEBISBox.getMultiCells(curBox);

	  vofItIrreg[dit()].define(notRegular,curEBISBox.getEBGraph());

	  const BaseIVFAB<Real>& curWeightBaseIVFAB = m_fluxWeight[dit()];
	  const BaseIVFAB<VoFStencil>& custStencilBaseIVFAB = m_custStencil[dit()];
	  BaseIVFAB<VoFStencil>& curStencilBaseIVFAB = m_fluxStencil[dit()];
	  BaseIVFAB<Real>&       custWeightBaseIVFAB  = m_custWeight[dit()];

	  for (VoFIterator vofit(notRegular,curEBGraph); vofit.ok(); ++vofit)
	    { 
	      const VolIndex& vof = vofit();

	      // Compute the bndryCentroid location in physical coordinates		  
	      const IntVect& iv = vof.gridIndex();
	      RealVect point = RealVect::Unit;
	      point *= 0.5;
	      point += iv;
	      point += curEBISBox.bndryCentroid(vof);
	      point *= m_dx;
	      point += RealVect::Zero; //a_probLo assumed 0;
		  	
	      isDirBC = m_func->isDirichlet(point);
	      if (!isDirBC)
		{
		  //Stencil info for current volume
		  const Real& curWeight = curWeightBaseIVFAB(vof,comp);
		  const VoFStencil& custSten = custStencilBaseIVFAB(vof,0);
		  const Real& areaFrac = curEBISBox.bndryArea(vof);
		  VoFStencil& curSten = curStencilBaseIVFAB(vof,0);		      		
		  Real&  custWeight  = custWeightBaseIVFAB(vof,0);
	
		  //Get stencil information for the vof on the opposite side of the EB
		  //hardwired source parameter for now, to be passed from higher level array later on

		  Real denom, eps, epsFlip;				
		  int signVof = 1;								
		  const IntVect& ivLo = iv - signVof*BASISV(1);		
		  if (curBox.contains(ivLo))
		    {
		      if (notRegular.contains(ivLo))
			{
			  eps = m_epsHi;
			  epsFlip = m_epsLo;
			}
		      else
			{
			  const IntVect& ivHi = iv + signVof*BASISV(1);
			  if (curBox.contains(ivHi))
			    {
			      if (notRegular.contains(ivHi))
				{
				  signVof = -1;	
				  eps = m_epsLo;
				  epsFlip= m_epsHi;				
				}
			      else
				{
				  pout() << "MixedPoissonEBBC::defineCustom, vofHi/Lo indexes outside of irregular IVS " << ivLo << ivHi << endl;
				  MayDay::Error("EB BC defineCustom vof's out of irregular IVS");
				}
			    }
			}
							
		    }
		  const IntVect& ivFlip = iv - signVof*BASISV(1);					 
		  Vector<VolIndex> vofsFlip = curEBISBox.getVoFs(ivFlip);	
		  CH_assert(vofsFlip.size() == 1);
		  const VolIndex& vofFlip =  vofsFlip[0];				

		  const Real& curWeightFlip = curWeightBaseIVFAB(vofFlip,comp);	
		  const VoFStencil& custStenFlip = custStencilBaseIVFAB(vofFlip,comp);	
		  const Real& areaFracFlip = curEBISBox.bndryArea(vofFlip);
			
		  //Stencil correction for phi_s
		  //add assert denom non-zero 
		      
		  denom = (eps*curWeight*areaFrac-epsFlip*curWeightFlip*areaFracFlip);
		  Real coef = -eps*curWeight*areaFrac/denom;	
		  Real coefFlip = epsFlip*curWeight*areaFrac/denom;	
			  			
		  curSten *= (1.0 + coef)/coefFlip;
		  curSten += custStenFlip;
		  curSten *= coefFlip;

		  if (s_areaFracWeighted)
		    {
		      denom *= curEBISBox.areaFracScaling(vof);
		    }				
		  custWeight = signVof/denom;
		}
	    }			
	}
    }
}

//--------------------------------------------------------------------------------------------------------------

LayoutData<BaseIVFAB<VoFStencil> >* MixedPoissonEBBC::getFluxStencil(int ivar)
{
  if ((m_ebType == 1) || (m_ebType == 10))
    {
      return &m_fluxStencil;
    }
  else
    {
      //Neumann
      return NULL;
    }
}

//--------------------------------------------------------------------------------------------------------------
/*
// combined Dirichlet & Neumann
*/

void MixedPoissonEBBC::applyEBFlux(EBCellFAB&                    a_lphi,
				   const EBCellFAB&              a_phi,
				   VoFIterator&                  a_vofit,
				   const LayoutData<IntVectSet>& a_cfivs,
				   const DataIndex&              a_dit,
				   const RealVect&               a_probLo,
				   const RealVect&               a_dx,
				   const Real&                   a_factor,
				   const bool&                   a_useHomogeneous,
				   const Real&                   a_time)
{
  CH_TIME("MixedPoissonEBBC::applyEBFlux");
  CH_assert(a_lphi.nComp() == 1 );
  CH_assert(a_phi.nComp() == 1);
  const EBISBox& ebisBox = a_phi.getEBISBox();
  int comp = 0;
  bool isDirBC;
  Real flux=0;
  const Real& time = *m_timePtr;
  bool m_finest = (m_dx[1] == m_dymin[1]); 
  //MDCancel

  Real buffer=0;

  for (a_vofit.reset(); a_vofit.ok(); ++a_vofit)
    {
      const VolIndex& vof = a_vofit();     
      isDirBC = false;

	
      if (m_dataBased)
        {
	  if (m_ebType == 1)
	    {
	      Real value;
	      const BaseIVFAB<Real>& curWeightBaseIVFAB = m_fluxWeight[a_dit];
	      const Real& curWeight = curWeightBaseIVFAB(vof,comp);			      
	      value = (*m_data)[a_dit](vof, 0);
	      flux = curWeight * value;
	    }
	  else
	    {
	      //m_ebType = 0
	      flux = (*m_data)[a_dit](vof, 0);
	    }
        }
      else if (m_ebType == 10)
        {
          // Compute the bndryCentroid location in physical coordinates	
          const IntVect& iv = vof.gridIndex();	  
	  RealVect point = RealVect::Unit;
          point *= 0.5;
          point += iv;
          point += ebisBox.bndryCentroid(vof);
          point *= m_dx;
          point += a_probLo;
	  RealVect normal = ebisBox.normal(vof);
		  	    
	  isDirBC = m_func->isDirichlet(point);
	  if (isDirBC)
	    {
	      const BaseIVFAB<Real>& curWeightBaseIVFAB = m_fluxWeight[a_dit];
	      const Real& curWeight = curWeightBaseIVFAB(vof,comp);

	      Real value;				
	      value = m_func->valueDir(point, normal, time, s_velComp);
	      flux = curWeight * value;
	      buffer = value;
	    }
	  else
	    // Dielectric-Gas EB interface condition handled through below source term and stencil	        
	    {

	      if (m_finest)
		{ 
		  const BaseIVFAB<Real>& curWeightBaseIVFAB = m_fluxWeight[a_dit];
		  const Real& curWeight = curWeightBaseIVFAB(vof,comp);
		        	
		  //hardwired source parameter from poisson.inputs now, to be passed from higher level array later on				
		  
		  //BaseIVFAB<Real>& rhosBaseIVFAB = (*m_rhos[curlev])[a_dit];
		  //Real& curRhos = rhosBaseIVFAB(vof,comp);
		  // LM integrating div(eps grad phi) = -e \sum_k Z_k n_k, I get [eps grad phi]_g - [eps grad phi]_s = -rhos, where rhos is the surface charge
		  // Thus, I get a (-) negative sign here, whereby MD had a positive one
		  Real src = -m_rhos[a_dit](vof,comp);
		  const BaseIVFAB<Real>& custWeightBaseIVFAB = m_custWeight[a_dit];
		  const Real& custWeight = custWeightBaseIVFAB(vof,comp);	

		  flux = curWeight * src * custWeight;				    
		}
	    }
        }
      else
        {
          if (m_onlyHomogeneous)
            {
              MayDay::Error("MixedPoissonEBBC::applyEBFlux called with undefined inhomogeneous BC");
            }
	  else if (m_ebType == 1)
	    {
	      Real value;
	      const BaseIVFAB<Real>& curWeightBaseIVFAB = m_fluxWeight[a_dit];
	      const Real& curWeight = curWeightBaseIVFAB(vof,comp);
	      value = m_dvalue;
	      flux = curWeight * value;
	    }
	  else
	    {
	      flux = m_nvalue;
	    }
        }

      const Real& areaFrac = ebisBox.bndryArea(vof);
      flux *= areaFrac;
      Real* lphiPtr = NULL;

      int offset;
      //get ghosted box of lphi
      Box boxLph = a_lphi.getSingleValuedFAB().box();
      //check to see if multi-valued cell or not
      if (ebisBox.numVoFs(vof.gridIndex()) > 1)
        {
          const IntVectSet& ivsLph = ebisBox.getMultiCells(boxLph);
          const EBGraph& ebgraph = ebisBox.getEBGraph();
          BaseIVFAB<Real> baseivfabLph(ivsLph, ebgraph, 1);

          offset = baseivfabLph.getIndex(vof, 0) - baseivfabLph.dataPtr(0);
          Real* multiValuedPtrLph = a_lphi.getMultiValuedFAB().dataPtr(0);
          lphiPtr  = multiValuedPtrLph + offset;
        }
      else
	{
          IntVect ncellsLph = boxLph.size();
	  const IntVect& smallendLph = boxLph.smallEnd();
	  IntVect ivLph = vof.gridIndex()  - smallendLph;
	  offset = ivLph[0] + ivLph[1]*ncellsLph[0] ;
#if CH_SPACEDIM==3
	  offset +=  ivLph[2]*ncellsLph[0]*ncellsLph[1];
#endif
	  Real* singleValuedPtrLph = a_lphi.getSingleValuedFAB().dataPtr();		  
	  lphiPtr  = singleValuedPtrLph + offset;     
        }
      Real& lphi = *lphiPtr;
      lphi += flux * a_factor; 
    }
  /*if(abs(buffer)>0 && time > s_timePrint){
    pout() << "MixedPoissonEBBC::applyEBFlux, m_timePtr=" << time << ", dirichlet value on EB electrode =" << buffer << endl;s_timePrint=time;}*/
	  
	
}

//--------------------------------------------------------------------------------------------------------------
/*
// Stencil functions - below only second order used for mixed BC
*/

bool
MixedPoissonEBBC::
getSecondOrderStencil(VoFStencil&          a_stencil,
                      Real&                a_weight,
                      Vector<VoFStencil>&  a_pointStencils,
                      Vector<Real>&        a_distanceAlongLine,
                      const VolIndex&      a_vof,
                      const EBISBox&       a_ebisBox,
                      const RealVect&      a_dx,
                      const IntVectSet&    a_cfivs)
{

  a_stencil.clear();
  bool dropOrder = false;
  EBArith::johanStencil(dropOrder, a_pointStencils, a_distanceAlongLine,
                        a_vof, a_ebisBox, a_dx, a_cfivs);
  if (dropOrder)
    {
      return true;
    }

  //if we got this far, sizes should be at least 2
  CH_assert(a_distanceAlongLine.size() >= 2);
  CH_assert(a_pointStencils.size() >= 2);
  Real x1 = a_distanceAlongLine[0];
  Real x2 = a_distanceAlongLine[1];
  //fit quadratic function to line and find the gradient at the origin
  //grad = (x2*x2*(phi1-phi0)-x1*x1(phi2-phi0))/(x2*x2*x1 - x1*x1*x2);
  Real denom = x2*x2*x1 - x1*x1*x2;
  //not done by reference because i want point stencils checkable externally.
  VoFStencil phi1Sten = a_pointStencils[0];
  VoFStencil phi2Sten = a_pointStencils[1];
  phi1Sten *=-x2*x2/denom;
  phi2Sten *= x1*x1/denom;
  //weight is the multiplier of the inhomogeneous value (phi0)
  a_weight =-x1*x1/denom + x2*x2/denom;
  a_stencil += phi1Sten;
  a_stencil += phi2Sten;
  //if we got this far, we have a second order stencil;
  return false;
}

void MixedPoissonEBBC::getFirstOrderStencil(VoFStencil&     a_stencil,
					    Real&           a_weight,
					    const VolIndex& a_vof,
					    const EBISBox&  a_ebisBox,
					    const RealVect& a_dx)
{
  CH_TIME("MixedPoissonEBBC::getFirstOrderStencil1");
  EBArith::getLeastSquaresGradSten(a_stencil, a_weight, a_vof, a_ebisBox, a_dx, m_domain);
}

void MixedPoissonEBBC::getSecondOrderStencil(VoFStencil&       a_stencil,
					     Real&             a_weight,
					     const VolIndex&   a_vof,
					     const EBISBox&    a_ebisBox,
					     const RealVect&   a_dx,
					     const IntVectSet& a_cfivs)
{
  Vector<VoFStencil>  pointStencils;
  Vector<Real>        distanceAlongLine;
  bool needToDropOrder = getSecondOrderStencil(a_stencil, a_weight,
                                               pointStencils, distanceAlongLine,
                                               a_vof, a_ebisBox, a_dx, a_cfivs);
  if (needToDropOrder)
    {
      getFirstOrderStencil(a_stencil,a_weight,a_vof,a_ebisBox,a_dx);
    }
}

//--------------------------------------------------------------------------------------------------------------
/*
// Factories:
*/

MixedPoissonEBBCFactory::MixedPoissonEBBCFactory()
{
  m_dvalue = 12345.6789;
  m_nvalue = 12345.6789;
  m_func = RefCountedPtr<BaseMixBCValue>();
  m_ebType = 1;
  m_onlyHomogeneous = false;
  m_order = 1;
  m_dymin = 12345.6789*RealVect::Unit;
  m_epsLo = 12345.6789;
  m_epsHi = 12345.6789;
}

MixedPoissonEBBCFactory::~MixedPoissonEBBCFactory()
{
}

void MixedPoissonEBBCFactory::setOrder(int a_order)
{
  CH_assert(a_order >= 1 && a_order <= 2);
  m_order = a_order;
}

void MixedPoissonEBBCFactory::setArguments(Real a_dvalue, Real a_nvalue, RefCountedPtr<BaseMixBCValue> a_func, int a_ebType, RealVect a_dymin, Real a_epsLo, Real a_epsHi)
{
  m_dvalue = a_dvalue;
  m_nvalue = a_nvalue;
  m_func = a_func;
  m_ebType = a_ebType;
  m_dymin = a_dymin;
  m_epsLo = a_epsLo;
  m_epsHi = a_epsHi;
  m_onlyHomogeneous = false;
}

// LM added time routine
void MixedPoissonEBBCFactory::setTime(Real* a_time)
{
  m_timePtr = a_time;
}

MixedPoissonEBBC* MixedPoissonEBBCFactory::create(const ProblemDomain& a_domain,
						  const EBISLayout&    a_layout,
						  const RealVect&      a_dx,
						  const IntVect*       a_ghostCellsPhi /*=0*/,
						  const IntVect*       a_ghostCellsRhs /*=0*/)
{
  CH_TIME("MixedPoissonEBBC::create");

  MixedPoissonEBBC* fresh = new MixedPoissonEBBC(a_domain,a_layout,a_dx,a_ghostCellsPhi,
						 a_ghostCellsRhs);
  fresh->setArguments(m_dvalue, m_nvalue, m_func, m_ebType, m_dymin, m_epsLo, m_epsHi);
  fresh->setOrder(m_order);
  // LM added call
  fresh->setTime(m_timePtr);
  return fresh;
}
/******************************************/
void MixedPoissonEBBC::setEBBCSource(const LevelData<BaseIVFAB<Real> >&     a_rhos)
{
  //LM setting the rhos. Using the MD routines defineCustom & applyEBFlux
  bool m_finest = (m_dx[1] == m_dymin[1]); 
  if (m_finest && (m_ebType == 10))
    {	

      int comp = 0;
      bool isDirBC;
      const DisjointBoxLayout& dbl = m_layout.getDisjointLayout();

      LayoutData<VoFIterator > vofItIrreg;
      vofItIrreg.define(dbl); // vofiterator cache

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
	{
	  const Box& curBox = dbl[dit()];
	  const EBISBox& curEBISBox = m_layout[dit()];
	  const EBGraph& curEBGraph = curEBISBox.getEBGraph();

	  IntVectSet notRegular;
	  notRegular |= curEBISBox.getIrregIVS  (curBox);
	  //assuming no multi-cells here
	  //notRegular |= curEBISBox.getMultiCells(curBox);
	  vofItIrreg[dit()].define(notRegular,curEBISBox.getEBGraph());

	  BaseIVFAB<Real>& LocalIVFAB = m_rhos[dit()];
	  const BaseIVFAB<Real>& ExternalIVFAB = a_rhos[dit()];

	  for (VoFIterator vofit(notRegular,curEBGraph); vofit.ok(); ++vofit)
	    { 
	      const VolIndex& vof = vofit();

	      // Compute the bndryCentroid location in physical coordinates		  
	      const IntVect& iv = vof.gridIndex();RealVect point = RealVect::Unit;point *= 0.5; point += iv;point += curEBISBox.bndryCentroid(vof);point *= m_dx;point += RealVect::Zero; //a_probLo assumed 0;
		  	
	      isDirBC = m_func->isDirichlet(point);
	      if (!isDirBC)
		{			
		  int signVof = 1;								
		  const IntVect& ivLo = iv - signVof*BASISV(1);		
		  if (curBox.contains(ivLo))
		    {
		      if (notRegular.contains(ivLo)) // means it is on the high side, duh!
			LocalIVFAB(vof,comp)	= ExternalIVFAB(vof,comp);
		      else
			{//Low side; find the Hi Side
			  const IntVect& ivHi = iv + signVof*BASISV(1);
			  if (!notRegular.contains(ivHi)) pout() << vof <<  " MixedPoissonEBBC::setEBBCSource (LM) messed up (it will probably crash) " << endl;
			  Vector<VolIndex> vofsFlip = curEBISBox.getVoFs(ivHi);	
			  CH_assert(vofsFlip.size() == 1);
			  const VolIndex& vofFlip =  vofsFlip[0];
			  LocalIVFAB(vof,comp)	= ExternalIVFAB(vofFlip,comp);
			}
							
		    }//curBox.contains(ivLo)
		}//!isDirBC
	    }//vofit
	}//dit
    }//m_finest && (m_ebType == 10)
}
#include "NamespaceFooter.H"
