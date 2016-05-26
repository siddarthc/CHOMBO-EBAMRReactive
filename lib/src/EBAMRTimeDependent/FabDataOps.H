#ifndef _FABDATAOPS_H_
#define _FABDATAOPS_H_

#include "BaseFab.H"
#include "BaseIFFAB.H"
#include "BaseIVFAB.H"
#include "REAL.H"
#include "EBCellFAB.H"
#include "EBFaceFAB.H"
#include "EBFluxFAB.H"
#include "NamespaceHeader.H"
#include "EBISBox.H"

class FabDataOps
{
public:

   FabDataOps()
   {
   }

   ~FabDataOps()
   {
   }


  static void getFabData(BaseFab<Real>&   a_Fab);
  
  static void getFabData(const BaseFab<Real>&   a_Fab);
 
  static void getFabData(BaseFab<Real>&   a_Fab, const EBISBox& a_ebisBox, bool verbose);

  static void getFabData(const BaseFab<Real>&   a_Fab, const EBISBox& a_ebisBox, bool verbose);

 
//  static void getFabData(BaseIVFAB<Real>&   a_Fab, bool verbose);

//  static void getFabData(const BaseIVFAB<Real>&   a_Fab, bool verbose);
 
  static void getFabData(EBCellFAB&   a_Fab, bool verbose);

  static void getFabData(const EBCellFAB& a_Fab, bool verbose);

  static void getFabData(const EBCellFAB& a_Fab, VolIndex vof);

  static void getFabData(EBFaceFAB&   a_Fab);

  static void getFabData(const EBFaceFAB& a_Fab);

  static void getFabData(BaseIVFAB<Real>& a_coveredFab, const Vector<VolIndex>& a_vofset, const Box& a_box);

  static void getLevelData(LevelData<EBCellFAB>& a_data, bool verbose);
  
  static void getLevelData(const LevelData<EBCellFAB>& a_data, bool verbose);
  
//  static void getFabData(EBFluxFAB& a_Fab, bool verbose);

//  static void getFabData(const EBFluxFAB& a_Fab, bool verbose);

  static bool checkFabData(BaseFab<Real>& a_Fab, const Box& a_box);

  static bool checkFabData(EBFaceFAB& a_Fab, const int& a_dir, const Box& a_box);

  static bool checkFabData(BaseIVFAB<Real>& a_coveredFab, const Vector<VolIndex>& a_vofset, const Box& a_box);

  static bool checkFabData(EBCellFAB& a_Fab);

  static bool checkFabData(const EBCellFAB& a_Fab);

  static bool checkFabData(EBFluxFAB& a_Fab, const Box& a_box);

//  static bool checkFabData(const EBFluxFAB& a_Fab);

private:

};
#include "NamespaceFooter.H"
#endif
