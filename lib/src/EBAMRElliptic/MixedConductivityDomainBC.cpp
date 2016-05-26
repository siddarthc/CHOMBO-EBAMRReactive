#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "BoxIterator.H"
#include "EBArith.H"
#include "EBArithF_F.H"
#include "PolyGeom.H"

#include "MixedConductivityDomainBC.H"
#include "DirichletPoissonDomainBCF_F.H"
#include "NamespaceHeader.H"


MixedConductivityDomainBC::MixedConductivityDomainBC()
  : m_onlyHomogeneous(true),
    m_dvalue(12345.6789),
    m_nvalue(12345.6789),
    m_func(RefCountedPtr<BaseMixBCValue>()),
    m_ebOrder(2)
{
}

MixedConductivityDomainBC::~MixedConductivityDomainBC()
{
}

/*************/
// BC identifiers: 0=Neumann, 1=Dirichlet, 10=Mixed
void MixedConductivityDomainBC::setArguments(Vector<Vector<int> >&  a_domMixBc, Real a_dvalue, Real a_nvalue, RefCountedPtr<BaseMixBCValue> a_func)
{
  m_domMixBc = a_domMixBc;
  m_dvalue = a_dvalue;
  m_nvalue = a_nvalue;
  m_func = a_func;
  m_onlyHomogeneous = false; 
}

void MixedConductivityDomainBC::setEBOrder(int a_ebOrder)
{
  CH_assert(a_ebOrder >= 1);
  CH_assert(a_ebOrder <= 2);
  m_ebOrder = a_ebOrder;
}
 
///
/****************
   This never gets called.  InflowOutflowPoissonDomainBC::getFaceVel takes care of it.
*/
void
MixedConductivityDomainBC::
getFaceVel(Real&                 a_faceFlux,
           const FaceIndex&      a_face,
           const EBFluxFAB&      a_vel,
           const RealVect&       a_probLo,
           const RealVect&       a_dx,
           const int&            a_idir,
           const int&            a_icomp,
           const Real&           a_time,
           const Side::LoHiSide& a_side,
           const bool&           a_doDivFreeOutflow)
{ 
  //int iside = -sign(a_side);
  //int jside = (iside+1)/2;
  int iside;
  if (a_side == Side::Lo)
    {
      iside = 1;
    }
  else
    {
      iside = -1;
    }
  int jside = (iside+1)/2;
  const Real& time = *m_timePtr;
  if ((m_domMixBc[jside][a_idir]==1) || (m_domMixBc[jside][a_idir]==10))
    {
      CH_assert(a_idir == a_face.direction());
      Real value;
      if (m_domMixBc[jside][a_idir]==10)
        //Mixed BC function of space
	{
	  RealVect pt;
	  IntVect iv = a_face.gridIndex(Side::Hi);
	  for (int idir = 0; idir < SpaceDim; idir++)
	    {
	      if (idir != a_face.direction())
		{
		  Real ptval = a_dx[a_idir]*(Real(iv[idir]) + 0.5);
		  pt[idir] = ptval;
		}
	      else
		{
		  pt[idir] = a_dx[a_idir]*(Real(iv[idir]));
		}
	    }
	  bool isDirBC;
	  isDirBC = m_func->isDirichlet(pt);
	  if (isDirBC)
	    {
	      RealVect normal = EBArith::getDomainNormal(a_idir, a_side);
	      value = m_func->valueDir(pt, normal, time,a_icomp);
	      a_faceFlux = value;
	    }
	  else
	    //Neumann
	    {
	      a_faceFlux =  a_vel[a_idir](a_face, a_icomp);
	    }
	}
      else
	//D not function of space
	{
	  a_faceFlux = m_dvalue;
	} 		
    }
  else
    //N not function of space
    {
      a_faceFlux =  a_vel[a_idir](a_face, a_icomp);
    }
}

void MixedConductivityDomainBC::getFaceFlux(BaseFab<Real>&   a_faceFlux,
				       const BaseFab<Real>&  a_phi,
				       const RealVect&       a_probLo,
				       const RealVect&       a_dx,
				       const int&            a_idir,
				       const Side::LoHiSide& a_side,
				       const DataIndex&      a_dit,
				       const Real&           a_time,
				       const bool&           a_useHomogeneous)
{
  CH_TIME("MixedConductivityDomainBC::getFaceFlux");
  CH_assert(a_phi.nComp() == 1);
  const Real& time = *m_timePtr;
  for (int comp=0; comp<a_phi.nComp(); comp++)
    {
      const Box& box = a_faceFlux.box();

      int iside;
      if (a_side == Side::Lo)
        {
          iside = 1;
        }
      else
        {
          iside = -1;
        }
      int jside = (iside+1)/2;

      if (a_useHomogeneous)
        {
          Real value = 0.0;

	  if (m_domMixBc[jside][a_idir]==1)
	    //Dirichlet not function of position
	    {
	      FORT_SETDIRICHLETFACEFLUX(CHF_FRA1(a_faceFlux,comp),
					CHF_CONST_FRA1(a_phi,comp),
					CHF_CONST_REAL(value),
					CHF_CONST_REALVECT(a_dx),
					CHF_CONST_INT(a_idir),
					CHF_CONST_INT(iside),
					CHF_BOX(box));
	    }
		   
	  else if (m_domMixBc[jside][a_idir]==10)
	    //Dirichlet function of position 
	    {
              Real ihdx;
              ihdx = 2.0 / a_dx[a_idir];
              BoxIterator bit(box);

              for (bit.begin(); bit.ok(); ++bit)
                {
		  IntVect iv = bit();
		  //Dirichlet (x)
		  IntVect ivNeigh = iv;
		  ivNeigh[a_idir] += sign(a_side);
		  const VolIndex vof      = VolIndex(iv,     0);
		  const VolIndex vofNeigh = VolIndex(ivNeigh,0);
		  const FaceIndex face = FaceIndex(vof,vofNeigh,a_idir);
		  const RealVect  point = EBArith::getFaceLocation(face,a_dx,a_probLo);
		  bool isDirBC;
		  isDirBC = m_func->isDirichlet(point);
		  if (isDirBC)
		    {
		      FORT_SETDIRICHLETFACEFLUX(CHF_FRA1(a_faceFlux,comp),
						CHF_CONST_FRA1(a_phi,comp),
						CHF_CONST_REAL(value),
						CHF_CONST_REALVECT(a_dx),
						CHF_CONST_INT(a_idir),
						CHF_CONST_INT(iside),
						CHF_BOX(box));
		    }
		  else
		    //Neumann (x)
		    {
		      a_faceFlux(iv,comp) = iside * value;
		    }
		}
            }					   	
	  else
	    //Neumann
	    {			
	      a_faceFlux.setVal(iside * value);
	    }
        }
      else
	//a_useHomogeneous=false
        {
          if (m_domMixBc[jside][a_idir]==10)
	    //fonctional
            {
              Real ihdx;
              ihdx = 2.0 / a_dx[a_idir];
              BoxIterator bit(box);

              for (bit.begin(); bit.ok(); ++bit)
                {
                  IntVect iv = bit();
		  //Dirichlet
		  IntVect ivNeigh = iv;
		  ivNeigh[a_idir] += sign(a_side);
		  const VolIndex vof      = VolIndex(iv,     0);
		  const VolIndex vofNeigh = VolIndex(ivNeigh,0);
		  const FaceIndex face = FaceIndex(vof,vofNeigh,a_idir);
		  const RealVect  point = EBArith::getFaceLocation(face,a_dx,a_probLo);
		  const RealVect normal = EBArith::getDomainNormal(a_idir,a_side);
		  bool isDirBC;
		  isDirBC = m_func->isDirichlet(point);
		  if (isDirBC)
		    {
		      Real value = m_func->valueDir(face,a_side,a_dit,point,normal,time,comp);
		      Real phiVal = a_phi(iv,comp);
		      a_faceFlux(iv,comp) = iside * ihdx * (phiVal - value);
		    }
		  else
		    //Neumann
		    {
		      RealVect point = EBArith::getIVLocation(iv,a_dx,a_probLo);
		      point[a_idir] -= iside * 0.5 * a_dx[a_idir];//point is now at the face center
		      RealVect normal = RealVect::Zero;
		      normal[a_idir] = iside;
		      a_faceFlux(iv,comp) = iside * m_func->valueNeu(iv,a_dit,point,normal,time,comp);
		    }
		}
	    }
          else
	    //non-functional D or N
            {
              if (m_onlyHomogeneous)
                {
                  MayDay::Error("MixedConductivityDomainBC::getFaceFlux called with undefined inhomogeneous BC");
                }

              Real value = m_dvalue;			  
	      if (m_domMixBc[jside][a_idir]==1)
		//Dirichlet
              	{
		  FORT_SETDIRICHLETFACEFLUX(CHF_FRA1(a_faceFlux,comp),
					    CHF_CONST_FRA1(a_phi,comp),
					    CHF_CONST_REAL(value),
					    CHF_CONST_REALVECT(a_dx),
					    CHF_CONST_INT(a_idir),
					    CHF_CONST_INT(iside),
					    CHF_BOX(box));
		}
	      else
		//Neumann - was * value = m_value
		{
		  a_faceFlux.setVal(iside * m_nvalue);
		}
            }
        }
    }
   //again, following the odd convention of EBAMRPoissonOp
  //(because I am reusing its BC classes),
  //the input flux here is CELL centered and the input box
  //is the box adjacent to the domain boundary on the valid side.
  //because I am not insane (yet) I will just shift the flux's box
  //over and multiply by the appropriate coefficient
  a_faceFlux.shiftHalf(a_idir, -sign(a_side));
  const Box& faceBox = a_faceFlux.box();
  const BaseFab<Real>&   regCoef = (*m_bcoef)[a_dit][a_idir].getSingleValuedFAB();
  int  isrc = 0;
  int  idst = 0;
  int  inum = 1;
  FORT_MULTIPLYTWOFAB(CHF_FRA(a_faceFlux),
                      CHF_CONST_FRA(regCoef),
                      CHF_BOX(faceBox),
                      CHF_INT(isrc),CHF_INT(idst),CHF_INT(inum));

  //shift flux back to cell centered land
  a_faceFlux.shiftHalf(a_idir,  sign(a_side));

} 

void MixedConductivityDomainBC::getFaceFlux(Real&                     a_faceFlux,
				       const VolIndex&       a_vof,
				       const int&            a_comp,
				       const EBCellFAB&      a_phi,
				       const RealVect&       a_probLo,
				       const RealVect&       a_dx,
				       const int&            a_idir,
				       const Side::LoHiSide& a_side,
				       const DataIndex&      a_dit,
				       const Real&           a_time,
				       const bool&           a_useHomogeneous)
{
  //int iside = -sign(a_side);
  //int jside = (iside+1)/2;
  int iside;
  if (a_side == Side::Lo)
    {
      iside = 1;
    }
  else
    {
      iside = -1;
    }
  int jside = (iside+1)/2; 
  const Real& time = *m_timePtr;	
  const EBISBox& ebisBox = a_phi.getEBISBox();
  const ProblemDomain& domainBox = ebisBox.getDomain();
  a_faceFlux = 0.0;
  Vector<FaceIndex> faces = ebisBox.getFaces(a_vof,a_idir,a_side);
  if (faces.size() > 0)
    {
      if (faces.size()==1)
        {
          IntVectSet cfivs;
          FaceStencil faceSten = EBArith::getInterpStencil(faces[0],
                                                           cfivs,
                                                           ebisBox,
                                                           domainBox);
          for (int isten=0; isten < faceSten.size(); isten++)
            {
              const Real& weight = faceSten.weight(isten);
              const FaceIndex& face = faceSten.face(isten);
              Real thisFaceFlux;
              const RealVect centroid = ebisBox.centroid(face);
	      if ((m_domMixBc[jside][a_idir]==1) || (m_domMixBc[jside][a_idir]==10))
		//note getFaceGradPhi calls getFaceFluxGradPhi for the Neumann part of mixed BC 10
		{
		  getFaceGradPhi(thisFaceFlux,face,a_comp,a_phi,a_probLo,a_dx,a_idir,
				 a_side,a_dit,time,false,centroid,a_useHomogeneous);
		}
	      else
		{
		  getFaceFluxGradPhi(thisFaceFlux,face,a_comp,a_phi,a_probLo,a_dx,a_idir,
				     a_side,a_dit,time,false,centroid,a_useHomogeneous);
		}
			 
              a_faceFlux += thisFaceFlux*weight;
            }
          a_faceFlux *= ebisBox.areaFrac(faces[0]);
        }
      else
        {
          MayDay::Error("MixedConductivityDomainBC::getHigherOrderFaceFlux has multi-valued faces (or could be 0 faces)");
        }
    }
  Real bcoave = 0;
  Real areaTot = 0;
  for (int iface = 0; iface < faces.size(); iface++)
    {
      Real areaFrac = ebisBox.areaFrac(faces[iface]);
      Real bcoFace  = (*m_bcoef)[a_dit][a_idir](faces[iface], 0);
      bcoave += areaFrac*bcoFace;
      areaTot += areaFrac;
    }
  if (areaTot > 1.0e-8)
    {
      bcoave /= areaTot;
    }
  a_faceFlux *= bcoave;
}

void MixedConductivityDomainBC::getFaceGradPhi(Real&            a_faceFlux,
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
					  const bool&           a_useHomogeneous)
{
  //int iside = -sign(a_side);
  //int jside = (iside+1)/2;
  int iside;
  if (a_side == Side::Lo)
    {
      iside = 1;
    }
  else
    {
      iside = -1;
    }
  int jside = (iside+1)/2;
  RealVect point = EBArith::getFaceLocation(a_face,a_dx,a_probLo);
  bool isDirBC;
  isDirBC = m_func->isDirichlet(point);
  const EBISBox& ebisBox = a_phi.getEBISBox();
	
  if ((m_domMixBc[jside][a_idir]==1) || ((m_domMixBc[jside][a_idir]==10) && (isDirBC)))
    {	
      const Real ihdx = 2.0 / a_dx[a_idir];
      Real value = -1.e99;
      if (a_useHomogeneous)
	{
	  value = 0.0;
	}
      else if (m_domMixBc[jside][a_idir]==10)
	{   
	  const RealVect normal = EBArith::getDomainNormal(a_idir,a_side);     			
	  value = m_func->valueDir(a_face,a_side,a_dit,point,normal,a_time,a_comp);
	}
    		
      else
	{
	  if (m_onlyHomogeneous)
	    {
	      MayDay::Error("MixedConductivityDomainBC::getFaceFlux called with undefined inhomogeneous BC");
	    }
	  value = m_dvalue;
	}

      const VolIndex& vof = a_face.getVoF(flip(a_side));
      a_faceFlux = iside * ihdx * (a_phi(vof,a_comp) - value);
    }
  else 
    // Neumann (0 & mixed case 10 with isDirBC = false)
    {
      getFaceFluxGradPhi(a_faceFlux,a_face,a_comp,a_phi,a_probLo,
			 a_dx,a_idir,a_side,a_dit,a_time,a_useAreaFrac,a_centroid,a_useHomogeneous);	
    }

  if (a_useAreaFrac)
    {
      MayDay::Error("DirichletPoissonDomainBC::getFaceFlux -- useAreaFrac=TRUE");
      a_faceFlux *= ebisBox.areaFrac(a_face);
    }
}
///
/****************
From Neumann, used in getFaceFlux and getFaceGradPhi, not used by Dirichlet
*/
void MixedConductivityDomainBC::getFaceFluxGradPhi(Real&                       a_faceFlux,
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
					      const bool&           a_useHomogeneous)
{
  const int iside = -sign(a_side);
  int jside = (iside+1)/2;  
  Real flux = -1.e99;
  if (a_useHomogeneous)
    {
      flux = 0.0;
    }
  else if (m_domMixBc[jside][a_idir]==10)
    {
      const RealVect normal = EBArith::getDomainNormal(a_idir,a_side);
      RealVect point = EBArith::getFaceLocation(a_face,a_dx,a_probLo);

      bool isDirBC;
      isDirBC = m_func->isDirichlet(point);
      if (!isDirBC)
	{
	  flux = m_func->valueNeu(a_face,a_side,a_dit,point,normal,a_time,a_comp);
	}
    }
  else
    {
      if (m_onlyHomogeneous)
	{
	  MayDay::Error("NeumannConductivityDomainBC::getFaceFlux called with undefined inhomogeneous BC");
	}
      flux = m_nvalue;
    }
  a_faceFlux = iside*flux;
}


/****************
 // From Dirichlet:
 */
void MixedConductivityDomainBC::getFluxStencil(VoFStencil&            a_stencil,
					  const VolIndex&        a_vof,
					  const int&             a_comp,
					  const RealVect&        a_dx,
					  const int&             a_idir,
					  const Side::LoHiSide&  a_side,
					  const EBISBox&         a_ebisBox)
{

  const ProblemDomain& domainBox = a_ebisBox.getDomain();
  Vector<FaceIndex> faces = a_ebisBox.getFaces(a_vof, a_idir, a_side);
  if (faces.size() > 0)
    {
      if (faces.size()==1)
        {
          CH_assert(faces[0].isBoundary());
	  // LM changes. If neumann BC: no changes to the stencils
	  const RealVect  point = EBArith::getFaceLocation(faces[0],a_dx,RealVect::Zero);
	  bool isDirBC = m_func->isDirichlet(point);
	  if(!isDirBC)
	    {
	      a_stencil.clear();
	      return;
	    }
	  //end of LM changes
          IntVectSet cfivs;
          FaceStencil faceSten = EBArith::getInterpStencil(faces[0],
                                                           cfivs,
                                                           a_ebisBox,
                                                           domainBox);
          for (int isten=0; isten < faceSten.size(); isten++)
            {
              const Real& weight = faceSten.weight(isten);
              const FaceIndex& face = faceSten.face(isten);

              VoFStencil thisStencil;
              getFluxStencil(thisStencil, face, a_comp, a_dx, a_idir, a_side, a_ebisBox);
              thisStencil *= weight;
              a_stencil += thisStencil;
            }
          a_stencil *= a_ebisBox.areaFrac(faces[0]);
        }
      else
        {
          MayDay::Error("MixedConductivityDomainBC::getFluxStencil has multi-valued faces");
        }
    }
  a_stencil *= 1./a_dx[a_idir];
  a_stencil *= Real(sign(a_side));//this sign is for which side is differenced
}

void MixedConductivityDomainBC::getFluxStencil(VoFStencil&            a_stencil,
					  const FaceIndex&       a_face,
					  const int&             a_comp,
					  const RealVect&        a_dx,
					  const int&             a_idir,
					  const Side::LoHiSide&  a_side,
					  const EBISBox&         a_ebisBox)
{
  // LM changes
  // if neumann BC: no changes to the stencils
  const RealVect  point = EBArith::getFaceLocation(a_face,a_dx,RealVect::Zero);
  bool isDirBC = m_func->isDirichlet(point);
  if(!isDirBC)
    {
      a_stencil.clear();
      return;
    }
  //end of LM changes

  if (m_ebOrder == 1)//pressure
    {
      getFirstOrderFluxStencil(a_stencil, a_face, a_comp,
                               a_dx, a_idir, a_side,
                               a_ebisBox);
    }
  else if (m_ebOrder == 2)//velocity
    {
      getSecondOrderFluxStencil(a_stencil, a_face, a_comp,
                                a_dx, a_idir, a_side,
                                a_ebisBox);
    }
  else
    {
      MayDay::Error("MixedConductivityDomainBC::getFluxStencil -- bad BC order");
    }
}

void MixedConductivityDomainBC::getFirstOrderFluxStencil(      VoFStencil&      a_stencil,
							  const FaceIndex&       a_face,
							  const int&             a_comp,
							  const RealVect&        a_dx,
							  const int&             a_idir,
							  const Side::LoHiSide&  a_side,
							  const EBISBox&         a_ebisBox)
{
  const VolIndex& vof = a_face.getVoF(flip(a_side));
  const Real isign = Real(sign(a_side));//this sign is for the extrapolation direction
  const Real weight = -isign*2.0/a_dx[a_idir];
  a_stencil.add(vof, weight, a_comp);
}

void MixedConductivityDomainBC::getSecondOrderFluxStencil(      VoFStencil&      a_stencil,
							   const FaceIndex&       a_face,
							   const int&             a_comp,
							   const RealVect&        a_dx,
							   const int&             a_idir,
							   const Side::LoHiSide&  a_side,
							   const EBISBox&         a_ebisBox)
{
  const VolIndex& vof = a_face.getVoF(flip(a_side));
  const Real isign = Real(sign(a_side));//this sign is for the extrapolation direction
  const Real weight = isign*1.0/(3.0*a_dx[a_idir]);
  Vector<FaceIndex> facesInsideDomain = a_ebisBox.getFaces(vof,a_idir,flip(a_side));
  if (facesInsideDomain.size() == 1)
    {
      const VolIndex& vofNextInsideDomain = facesInsideDomain[0].getVoF(flip(a_side));
      a_stencil.add(vofNextInsideDomain,      weight, a_comp);
      a_stencil.add(                vof, -9.0*weight, a_comp);
    }
  else
    {
      a_stencil.add(vof, -6.0*weight, a_comp);
    }
}

/****************
 // Factory Routines:
 */
MixedConductivityDomainBCFactory::MixedConductivityDomainBCFactory()
  : m_dvalue(12345.6789),
    m_nvalue(12345.6789),
    m_func(RefCountedPtr<BaseMixBCValue>()),
    m_ebOrder(2),
    m_onlyHomogeneous(true)
{
  m_domMixBc.resize(2); 
  for (int i = 0; i < SpaceDim; i++) {m_domMixBc[i].resize(SpaceDim,0);}
}

MixedConductivityDomainBCFactory::~MixedConductivityDomainBCFactory()
{
}

void MixedConductivityDomainBCFactory::setArguments(Vector<Vector<int> >& a_domMixBc, Real a_dvalue, Real a_nvalue, RefCountedPtr<BaseMixBCValue> a_func)
{
  m_domMixBc = a_domMixBc;
  m_dvalue = a_dvalue;
  m_nvalue = a_nvalue;
  m_func = a_func;
  m_onlyHomogeneous = false;
}

//use this for order of domain boundary
void MixedConductivityDomainBCFactory::setEBOrder(int a_ebOrder)
{
  CH_assert(a_ebOrder >= 1);
  CH_assert(a_ebOrder <= 2);

  m_ebOrder = a_ebOrder;
}

void MixedConductivityDomainBCFactory::setTime(Real* a_time)
{
  m_timePtr = a_time;
}

MixedConductivityDomainBC*
MixedConductivityDomainBCFactory::
create(const ProblemDomain& a_domain,
       const EBISLayout&    a_layout,
       const RealVect&      a_dx)
{
  MixedConductivityDomainBC* newBC = new MixedConductivityDomainBC();
  newBC->setEBOrder(m_ebOrder);  
  //  if (!m_onlyHomogeneous)
  //    {
  newBC->setArguments(m_domMixBc, m_dvalue, m_nvalue, m_func);
  newBC->setTime(m_timePtr);
  //	}
  return newBC;
}
#include "NamespaceFooter.H"
