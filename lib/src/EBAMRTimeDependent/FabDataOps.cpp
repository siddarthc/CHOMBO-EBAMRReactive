#include <cmath>

#include "FabDataOps.H"
#include "BoxIterator.H"
#include "VoFIterator.H"
#include "FaceIterator.H"
#include "IntVect.H"
#include "IntVectSet.H"
#include "EBISBox.H"
#include "parstream.H"
#include "Vector.H"
//#include "REAL.H"

void FabDataOps::getFabData(BaseFab<Real>&   a_Fab)
{
    int nComp = a_Fab.nComp();
    Box domain = a_Fab.box();
   // EBISBox ebisBox;
    BoxIterator bit(domain);
    for (bit.begin(); bit.ok(); ++bit)
    {
      const IntVect& iv = bit();
//       vector<Real> data(nComp);
//       Real data[nComp];
//       a_Fab.getVal(&data,iv);
       for (int ivar = 0; ivar < nComp; ivar++)
        { 
          Real val = a_Fab(iv,ivar); 
          pout() << "    icomp = " << ivar << " vof = " << iv << " val = " << val << std::endl;
        }
    }
}

void FabDataOps::getFabData(const BaseFab<Real>&   a_Fab)
{
    int nComp = a_Fab.nComp();
    Box domain = a_Fab.box();
    //EBISBox ebisBox;
    BoxIterator bit(domain);
    for (bit.begin(); bit.ok(); ++bit)
    {
       const IntVect& iv = bit();
//       vector<Real> data(nComp);
//       Real data[nComp];
//       a_Fab.getVal(&data,iv);

       for (int ivar = 0; ivar < nComp; ivar++)
        {
          Real val = a_Fab(iv,ivar);
          pout() << "    icomp = " << ivar << " vof = " << iv << " val = " << val << std::endl;
        }
    }
}

void FabDataOps::getFabData(BaseFab<Real>&   a_Fab, const EBISBox& a_ebisBox, bool verbose)
{
    int nComp = a_Fab.nComp();
    const EBISBox& ebisBox = a_ebisBox;
    ProblemDomain probDomain = ebisBox.getDomain();
    Box domain = a_Fab.box() & probDomain;
    BoxIterator bit(domain);
    for (bit.begin(); bit.ok(); ++bit)
    {
       const IntVect& iv = bit();
//       vector<Real> data(nComp);
//       Real data[nComp];
//       a_Fab.getVal(&data,iv);

     if (verbose){
      pout() << " vof = " << iv;
      if (probDomain.contains(iv))
      {       
       if (ebisBox.isRegular(iv))
        {
          pout() << "regular cell" << std::endl;
        }
       else if (ebisBox.isIrregular(iv))
        {
          pout() << "irregular cell" << std::endl;
        }
       else if (ebisBox.isCovered(iv))
        {
          pout() << "covered cell" << std::endl;
        }
       else
        {
          pout() << "bogus cell type" << std::endl;
        }
      }
     else
      {
        pout() << "ghost cell" << std::endl;
      }

     }

       for (int ivar = 0; ivar < nComp; ivar++)
        {
          Real val = a_Fab(iv,ivar);
          pout() << "    icomp(mod) = " << ivar << " vof = " << iv << " val = " << val << std::endl;
        }
    }
}


void FabDataOps::getFabData(const BaseFab<Real>&   a_Fab, const EBISBox& a_ebisBox, bool verbose)
{
    int nComp = a_Fab.nComp();
    const EBISBox& ebisBox = a_ebisBox;
    ProblemDomain probDomain = ebisBox.getDomain();
    Box domain = a_Fab.box() & probDomain;
    BoxIterator bit(domain);
    for (bit.begin(); bit.ok(); ++bit)
    {
       const IntVect& iv = bit();
//       vector<Real> data(nComp);
//       Real data[nComp];
//       a_Fab.getVal(&data,iv);

      if (verbose){
      pout() << " vof = " << iv;
      if (probDomain.contains(iv))
      { 
       if (ebisBox.isRegular(iv))
        {
          pout() << "regular cell" << std::endl;
        }
       else if (ebisBox.isIrregular(iv))
        {
          pout() << "irregular cell" << std::endl;
        }
       else if (ebisBox.isCovered(iv))
        {
          pout() << "covered cell" << std::endl;
        }
       else
        {
          pout() << "bogus cell type" << std::endl;
        }
      }
     else
      {
        pout() << "ghost cell" << std::endl;
      }

     }

     for (int ivar = 0; ivar < nComp; ivar++)
      {
        Real val = a_Fab(iv,ivar);
        pout() << "    icomp = " << ivar << " vof = " << iv << " val = " << val << std::endl;
      }
    }
}


/*
void FabDataOps::getFabData(BaseIVFAB<Real>&   a_Fab, bool verbose)
{
    int nComp = a_Fab.nComp();
    Box domain = a_Fab.box();
    EBISBox ebisBox;
    BoxIterator bit(domain);
    for (bit.begin(); bit.ok(); ++bit)
    {
       const IntVect& iv = bit();
//       vector<Real> data(nComp);
//       Real data[nComp];
//       a_Fab.getVal(&data,iv);
      if(verbose){
       pout() << " vof = " << iv;
       if (ebisBox.isRegular(iv))
        {
          pout() << "regular cell" << std::endl;
        }
       else if (ebisBox.isIrregular(iv))
        {
          pout() << "irregular cell" << std::endl;
        }
       else if (ebisBox.isCovered(iv))
        {
          pout() << "covered cell" << std::endl;
        }
       else
        {
          pout() << "bogus cell type" << std::endl;
        }
       }

       for (int ivar = 0; ivar < nComp; ivar++)
        {
          Real val = a_Fab(iv,ivar);
          pout() << "    icomp = " << ivar << " vof = " << iv << " val = " << val << std::endl;
        }
    }
}


void FabDataOps::getFabData(const BaseIVFAB<Real>&   a_Fab, bool verbose)
{
    int nComp = a_Fab.nComp();
    Box domain = a_Fab.box();
    EBISBox ebisBox;
    BoxIterator bit(domain);
    for (bit.begin(); bit.ok(); ++bit)
    {
       const IntVect& iv = bit();
//       vector<Real> data(nComp);
//       Real data[nComp];
//       a_Fab.getVal(&data,iv);
      if(verbose){
       pout() << " vof = " << iv;
       if (ebisBox.isRegular(iv))
        {
          pout() << "regular cell" << std::endl;
        }
       else if (ebisBox.isIrregular(iv))
        {
          pout() << "irregular cell" << std::endl;
        }
       else if (ebisBox.isCovered(iv))
        {
          pout() << "covered cell" << std::endl;
        }
       else
        {
          pout() << "bogus cell type" << std::endl;
        }
       }
       for (int ivar = 0; ivar < nComp; ivar++)
        {
          Real val = a_Fab(iv,ivar);
          pout() << "    icomp = " << ivar << " vof = " << iv << " val = " << val << std::endl;
        }
    }
}

*/
void FabDataOps::getFabData(EBCellFAB&   a_Fab, bool verbose)
{
   Box fabBox = a_Fab.box();
   BoxIterator bit(fabBox);
   const EBISBox& ebisBox = a_Fab.getEBISBox();
   ProblemDomain probDomain = ebisBox.getDomain();  

   if (verbose){
    for (bit.begin(); bit.ok(); ++bit)
    {
     const IntVect& iv = bit();
     pout() << " vof = " << iv;
     if (probDomain.contains(iv))
      {

       if (ebisBox.isRegular(iv))
        {
          pout() << "regular cell" << std::endl;
        }
       else if (ebisBox.isIrregular(iv))
        {
          pout() << "irregular cell" << std::endl;
        }
       else if (ebisBox.isCovered(iv))
        {
          pout() << "covered cell" << std::endl;
        }
       else
        {
          pout() << "bogus cell type" << std::endl;
        }
      }
     else
      {
        pout() << "ghost cell" << std::endl;
      }
    }
  }

   BaseFab<Real>& regFab = a_Fab.getSingleValuedFAB();
 //  pout() << "regular FAB values" << std::endl;
   getFabData(regFab,ebisBox,0);
/*
   BaseIVFAB<Real>& irregFab = a_Fab.getMultiValuedFAB();
   pout() << "irregular FAB values" << std::endl;
   getFabData(irregFab,0);
*/

   int nComp = regFab.nComp();

   IntVectSet ivsMulti = ebisBox.getMultiCells(fabBox);
//   pout() << "irregular FAB values" << std::endl;
   for (VoFIterator vofit(ivsMulti, ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      IntVect iv = vof.gridIndex();
      for (int ivar = 0; ivar < nComp; ivar++)
        {
          Real val = a_Fab(vof, ivar);
          pout() << "    icomp = " << ivar << " vof = " << iv << " val = " << val << std::endl;
        }
     }

}

void FabDataOps::getFabData(const EBCellFAB&   a_Fab, bool verbose)
{ 
   Box fabBox = a_Fab.box();
   BoxIterator bit(fabBox);
   const EBISBox& ebisBox = a_Fab.getEBISBox();
   ProblemDomain probDomain = ebisBox.getDomain();

   if (verbose){
    for (bit.begin(); bit.ok(); ++bit)
    {
     const IntVect& iv = bit();
     pout() << " vof = " << iv;
     if (probDomain.contains(iv))
      { 
       if (ebisBox.isRegular(iv))
        {
          pout() << "regular cell" << std::endl;
        }
       else if (ebisBox.isIrregular(iv))
        {
          pout() << "irregular cell" << std::endl;
        }
       else if (ebisBox.isCovered(iv))
        {
          pout() << "covered cell" << std::endl;
        }
       else
        {
          pout() << "bogus cell type" << std::endl;
        }
      }
     else
      {
        pout() << "ghost cell" << std::endl;
      }
     
    }
  }

   BaseFab<Real>& regFab = (BaseFab<Real>&)a_Fab.getSingleValuedFAB();
   pout() << "regular FAB values" << std::endl;
   getFabData(regFab,ebisBox,0);
/*
   BaseIFFAB<Real>& irregFab = a_Fab.getMultiValuedFAB();
   pout() << "irregular FAB values" << std::endl;
   getFabData(irregFab,1);
*/
   
   int nComp = regFab.nComp();

   IntVectSet ivsMulti = ebisBox.getMultiCells(fabBox);
   pout() << "irregular FAB values" << std::endl;
   for (VoFIterator vofit(ivsMulti, ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      IntVect iv = vof.gridIndex();
      for (int ivar = 0; ivar < nComp; ivar++)
        {
          Real val = a_Fab(vof, ivar);
          pout() << "    icomp = " << ivar << " vof = " << iv << " val = " << val << std::endl;
        }
     }

}

void FabDataOps::getFabData(const EBCellFAB&   a_Fab, VolIndex a_vof)
{
   Box fabBox = a_Fab.box();
//   BoxIterator bit(fabBox);
   const EBISBox& ebisBox = a_Fab.getEBISBox();
//   ProblemDomain probDomain = ebisBox.getDomain();
   IntVect iv = a_vof.gridIndex();
   int nComp = a_Fab.nComp();
   for (int ivar = 0; ivar < nComp; ivar++)
        {
          Real val = a_Fab(a_vof, ivar);
          pout() << "    icomp = " << ivar << " vof = " << iv << " val = " << val << std::endl;
        }

}

void FabDataOps::getFabData(EBFaceFAB&   a_Fab)
{
  int nComp = a_Fab.nComp();
  //  Box cellBox = enclosedCells(a_box);
  int dir = a_Fab.direction();
 //  CH_assert(a_Fab.direction() == a_dir);
 //  CH_assert(a_Fab.getRegion().contains(a_box));

   BaseFab<Real>& regFab = a_Fab.getSingleValuedFAB();
   getFabData(regFab); 
  //check regular stuff 
 //  dataIsNANINF = checkFabData(regFab, a_box);
 
   //check irregular stuff 
  
 //  Box grownBox = cellBox; 
 //  grownBox.grow(a_dir,1);
/* 
   const EBISBox& ebisBox = a_Fab.getEBISBox();
 //  IntVectSet ivsMulti =    
   FaceStop::WhichFaces stopCritGrid = FaceStop::SurroundingWithBoundary;
   IntVectSet ivs(a_Fab.getRegion());
   for (FaceIterator faceit(ivs,ebisBox.getEBGraph(),dir,stopCritGrid); faceit.ok(); ++faceit)
    {
      const FaceIndex& face = faceit();
      IntVect vof = face.gridIndex(Side::Lo);
      for (int ivar = 0; ivar < nComp; ivar++)
       { 
         Real val = a_Fab(face,ivar);
         pout() << "    icomp = " << ivar << " vof = " << vof << " val = " << val << std::endl; 
       } 
     }
*/ 
}

void FabDataOps::getFabData(const EBFaceFAB&   a_Fab)
{
   int nComp = a_Fab.nComp();
  //  Box cellBox = enclosedCells(a_box);
   int dir = a_Fab.direction(); 
 //  CH_assert(a_Fab.direction() == a_dir);
 //  CH_assert(a_Fab.getRegion().contains(a_box));
 
 //  BaseFab<Real>& regFab = a_Fab.getSingleValuedFAB();
   //check regular stuff  
 //  dataIsNANINF = checkFabData(regFab, a_box);
  
   //check irregular stuff 
   
 //  Box grownBox = cellBox; 
 //  grownBox.grow(a_dir,1); 
   const EBISBox& ebisBox = a_Fab.getEBISBox();
 //  IntVectSet ivsMulti =    
   FaceStop::WhichFaces stopCritGrid = FaceStop::SurroundingWithBoundary;
   IntVectSet ivs(a_Fab.getRegion());
   for (FaceIterator faceit(ivs,ebisBox.getEBGraph(),dir,stopCritGrid); faceit.ok(); ++faceit)
    { 
      const FaceIndex& face = faceit();
      IntVect vof = face.gridIndex(Side::Lo);
      for (int ivar = 0; ivar < nComp; ivar++)
       {
         Real val = a_Fab(face,ivar); 
         pout() << "    icomp = " << ivar << " vof = " << vof << " val = " << val << std::endl; 
       }  
     }  
}

void FabDataOps::getFabData(BaseIVFAB<Real>& a_coveredFab, const Vector<VolIndex>& a_vofset, const Box& a_box)
{
   int nComp = a_coveredFab.nComp();
   for (int ivof = 0; ivof < a_vofset.size(); ivof++)
    {
      const VolIndex& vof = a_vofset[ivof];
      if (a_box.contains(vof.gridIndex()))
       {
         for (int ivar =0; ivar < nComp; ivar++)
          {
            Real val = a_coveredFab(vof,ivar);
            pout() << "    icomp = " << ivar << " vof = " << vof << " val = " << val << std::endl;
          
       }
     }
    }
}
 

void FabDataOps::getLevelData(LevelData<EBCellFAB>& a_data, bool verbose)
{
  int nComp = a_data.nComp();
  for (DataIterator dit=a_data.dataIterator();dit.ok();++dit)
   {
     const EBCellFAB& dataEBFAB = a_data[dit()];
     const Box& region = dataEBFAB.getRegion();
     IntVectSet ivsBox(region);
     const EBISBox& ebisBox = dataEBFAB.getEBISBox();
     for (VoFIterator vofit(ivsBox, ebisBox.getEBGraph());vofit.ok();++vofit)
      {
        const VolIndex& vof = vofit();
        IntVect iv = vof.gridIndex();
        if (verbose)
         {
          pout() << " vof = " << iv;
          if (ebisBox.isRegular(iv))
           {
             pout() << "regular cell" << std::endl;
           }
          else if (ebisBox.isIrregular(iv))
           {
             pout() << "irregular cell" << std::endl;
           }
          else if (ebisBox.isCovered(iv))
           {
             pout() << "covered cell" << std::endl;
           }
          else
           {
             pout() << "bogus cell type" << std::endl;
           }
         }

       for (int ivar = 0; ivar < nComp; ivar++)
        {
          Real val = dataEBFAB(vof,ivar);
          pout() << "    icomp = " << ivar << " vof = " << iv << " val = " << val << std::endl;
        }
    }
  }
}

void FabDataOps::getLevelData(const LevelData<EBCellFAB>& a_data, bool verbose)
{
  int nComp = a_data.nComp();
  for (DataIterator dit=a_data.dataIterator();dit.ok();++dit)
   {
     const EBCellFAB& dataEBFAB = a_data[dit()];
     const Box& region = dataEBFAB.getRegion();
     IntVectSet ivsBox(region);
     const EBISBox& ebisBox = dataEBFAB.getEBISBox();
     for (VoFIterator vofit(ivsBox, ebisBox.getEBGraph());vofit.ok();++vofit)
      {
        const VolIndex& vof = vofit();
        IntVect iv = vof.gridIndex();
        if (verbose)
         {
          pout() << " vof = " << iv;
          if (ebisBox.isRegular(iv))
           {
             pout() << "regular cell" << std::endl;
           }
          else if (ebisBox.isIrregular(iv))
           {
             pout() << "irregular cell" << std::endl;
           }
          else if (ebisBox.isCovered(iv))
           {
             pout() << "covered cell" << std::endl;
           }
          else
           {
             pout() << "bogus cell type" << std::endl;
           }
         }
       
       for (int ivar = 0; ivar < nComp; ivar++)
        {
          Real val = dataEBFAB(vof,ivar);
          pout() << "    icomp = " << ivar << " vof = " << iv << " val = " << val << std::endl;
        }
    }
   }
}

bool FabDataOps::checkFabData(BaseFab<Real>& a_Fab, const Box& a_box)
{
    bool dataIsNANINF = false;

    int nComp = a_Fab.nComp();
    Box domain = a_Fab.box();
   
   // CH_assert(domain.contains(a_box));
   // EBISBox ebisBox;
    BoxIterator bit(a_box);
    for (bit.begin(); bit.ok(); ++bit)
    {
      const IntVect& iv = bit();
//       vector<Real> data(nComp);
//       Real data[nComp];
//       a_Fab.getVal(&data,iv);
       for (int ivar = 0; ivar < nComp; ivar++)
        {
          Real val = a_Fab(iv,ivar);
          if (isnan(val) || isinf(val) || Abs(val)>1.e40)
           {
             pout() << "    icomp = " << ivar << " vof = " << iv << " val = " << val << std::endl;
             dataIsNANINF = true;
           }
        }
    }
    if (dataIsNANINF)
      {
        MayDay::Warning("Found a NAN or Infinity.");
      }
    return dataIsNANINF;
}
 

bool FabDataOps::checkFabData(EBFaceFAB& a_Fab, const int& a_dir, const Box& a_box)
{
   bool dataIsNANINF = false;
   int nComp = a_Fab.nComp();
  // Box cellBox = enclosedCells(a_box);
   const Box fabBox = a_Fab.getRegion();
 
   CH_assert(a_Fab.direction() == a_dir);
  // CH_assert(a_Fab.getRegion().contains(a_box));

   BaseFab<Real>& regFab = a_Fab.getSingleValuedFAB();
   //check regular stuff
   dataIsNANINF = checkFabData(regFab, a_box);

   //check irregular stuff
 
 //  Box grownBox = cellBox;
 //  grownBox.grow(a_dir,1);
/*
   const EBISBox& ebisBox = a_Fab.getEBISBox();
   //IntVectSet ivsMulti = ebisBox.getMultiCells(cellBox); 
   //IntVectSet ivs(fabBox);  
   FaceStop::WhichFaces stopCritGrid = FaceStop::SurroundingNoBoundary;
   IntVectSet ivs(a_Fab.getRegion());
   for (FaceIterator faceit(ivs,ebisBox.getEBGraph(),a_dir,stopCritGrid); faceit.ok(); ++faceit)
    {
      const FaceIndex& face = faceit();
      //IntVect vof = face.gridIndex(Side::Lo);
      for (int ivar = 0; ivar < nComp; ivar++)
       {
         Real val = a_Fab(face,ivar);
         if (isnan(val) || isinf(val) || Abs(val)>1.e40)
          {
            pout() << "    icomp = " << ivar << "val = " << val << std::endl;
            dataIsNANINF = true;
            MayDay::Warning("Found a NAN or Infinity.");
          }
       }
     }
*/
/*   if(dataIsNANINF)
     {
       MayDay::Warning("Found a NAN or Infinity.");
     }*/
    return dataIsNANINF;
}

bool FabDataOps::checkFabData(BaseIVFAB<Real>& a_coveredFab, const Vector<VolIndex>& a_vofset, const Box& a_box)
{
   bool dataIsNANINF = false;
   int nComp = a_coveredFab.nComp();
   for (int ivof = 0; ivof < a_vofset.size(); ivof++)
    {
      const VolIndex& vof = a_vofset[ivof];
      if (a_box.contains(vof.gridIndex()))
       {
         for (int ivar =0; ivar < nComp; ivar++)
          {
            Real val = a_coveredFab(vof,ivar);
            if (isnan(val) || isinf(val) || Abs(val)>1.e40)
          {
            pout() << "    icomp = " << ivar << " vof = " << vof << " val = " << val << std::endl;
            dataIsNANINF = true;
          }
       }
     }
    }
   if(dataIsNANINF)
     {
       MayDay::Warning("Found a NAN or Infinity.");
     }
    return dataIsNANINF;
}

bool FabDataOps::checkFabData(EBCellFAB& a_Fab)
{
   bool dataIsNANINF = false;
   int nComp = a_Fab.nComp();
   const Box& region = a_Fab.getRegion();
   IntVectSet ivs(region);
   const EBISBox& ebisBox = a_Fab.getEBISBox();
   for (VoFIterator vofit(ivs, ebisBox.getEBGraph());vofit.ok();++vofit)
     {
       const VolIndex& vof = vofit();
       for (int ivar =0; ivar < nComp; ivar++)
          {
            Real val = a_Fab(vof,ivar);
            if (isnan(val) || isinf(val) || Abs(val)>1.e40)
             {
              pout() << "    icomp = " << ivar << " vof = " << vof << " val = " << val << std::endl;
              dataIsNANINF = true;
             }
         } 
     }
   if(dataIsNANINF)
     {
       MayDay::Warning("Found a NAN or Infinity.");
     }
    return dataIsNANINF;
} 

bool FabDataOps::checkFabData(const EBCellFAB& a_Fab)
{
   bool dataIsNANINF = false;
   int nComp = a_Fab.nComp(); 
   const Box& region = a_Fab.getRegion(); 
   IntVectSet ivs(region); 
   const EBISBox& ebisBox = a_Fab.getEBISBox();
   for (VoFIterator vofit(ivs, ebisBox.getEBGraph());vofit.ok();++vofit)
     { 
       const VolIndex& vof = vofit(); 
       for (int ivar =0; ivar < nComp; ivar++)
          {
            Real val = a_Fab(vof,ivar);
            if (isnan(val) || isinf(val) || Abs(val)>1.e40)
          { 
            pout() << "    icomp = " << ivar << " vof = " << vof << " val = " << val << std::endl;
            dataIsNANINF = true;
          } 
       } 
     } 
   if(dataIsNANINF) 
     { 
       MayDay::Warning("Found a NAN or Infinity.");
     } 
    return dataIsNANINF; 
}  

bool FabDataOps::checkFabData(EBFluxFAB& a_Flux, const Box& a_box)
{
   bool dataIsNANINF = false;
   bool a = false;
   for (int idir = 0; idir < CH_SPACEDIM; idir++)
    {
      EBFaceFAB& fab = a_Flux[idir];
      dataIsNANINF = checkFabData(fab,idir,a_box);
      if(dataIsNANINF)
       {
         //MayDay::Warning("Found a NAN or Infinity.");
         a = true;
       }
     }
    
   return a;
}

   
    
   
