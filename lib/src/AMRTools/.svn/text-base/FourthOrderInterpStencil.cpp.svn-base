#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// #include <cstdio>

#include "FourthOrderInterpStencil.H"
#include "BoxIterator.H"
#include "CFMatrixF_F.H"

//////////////////////////////////////////////////////////////////////////////
// Constructor - set up some defaults
FourthOrderInterpStencil::FourthOrderInterpStencil()
{
  m_defined = false;
}

//////////////////////////////////////////////////////////////////////////////
// Destructor - free up storage
FourthOrderInterpStencil::~FourthOrderInterpStencil()
{
}


//////////////////////////////////////////////////////////////////////////////
// Define the object
void FourthOrderInterpStencil::define(
                                      const IntVect&         a_bdryOffset,
                                      const int&             a_refineCoarse,
                                      const int&             a_degree)
{
  m_bdryOffset = a_bdryOffset;
  m_refineCoarse = a_refineCoarse;
  m_degree = a_degree;

  m_baseFineBox.define(IntVect::Zero, (m_refineCoarse - 1) * IntVect::Unit);

  { // Find m_coarseBaseIndices.
    m_coarseBaseIndices.clear();
    // INNER SET:
    // First find the center of the inner set.
    IntVect center = IntVect::Zero;
    for (int idir = 0; idir < SpaceDim; idir++)
      {
        switch (m_bdryOffset[idir])
          {
          case -1:
            { // boundary is immediately to left
              center[idir] = 1;
              break;
            }
          case 1:
            { // boundary is immediately to right
              center[idir] = -1;
              break;
            }
          }
      }
    Box innerBox = Box(center - IntVect::Unit,
                       center + IntVect::Unit);
    for (BoxIterator bit(innerBox); bit.ok(); ++bit)
      {
        IntVect iv = bit();
        for (int idirComp = 0; idirComp < SpaceDim; idirComp++)
          m_coarseBaseIndices.push_back(iv[idirComp]);
      }

    // OUTER SET:
    for (int idir = 0; idir < SpaceDim; idir++)
      { // Add cells beyond innerBox, if possible, in dimension idir
        int bdryOffsetDir = m_bdryOffset[idir];
        if (bdryOffsetDir <= 0)
          { // point available on right
            IntVect iv = IntVect::Zero;
            iv[idir] = center[idir] + 2;
            for (int idirComp = 0; idirComp < SpaceDim; idirComp++)
              m_coarseBaseIndices.push_back(iv[idirComp]);
          }
        if (bdryOffsetDir >= 0)
          { // point available on left
            IntVect iv = IntVect::Zero;
            iv[idir] = center[idir] - 2;
            for (int idirComp = 0; idirComp < SpaceDim; idirComp++)
              m_coarseBaseIndices.push_back(iv[idirComp]);
          }
      }
  }
  m_stencilSize = m_coarseBaseIndices.size() / SpaceDim;

  m_coarseToFineFab.define(m_baseFineBox, m_stencilSize);

  // number of nonzero powers:  sum > 0 and sum <= degree
  int numPowers = 0;
  for (int isum = 1; isum <= m_degree; isum++)
    {
      // number of ways to put
      // isum indistinguishable balls into
      // SpaceDim distinguishable boxes
      // is equal to C(isum + SpaceDim - 1, SpaceDim - 1)
      int icount = 1;
      for (int ib = 1; ib <= SpaceDim-1; ib++)
        {
          icount *= (isum + ib);
        }
      for (int ib = 1; ib <= SpaceDim-1; ib++)
        {
          icount /= ib;
        }
      numPowers += icount;
    }

  int numFinePoints = m_baseFineBox.numPts();

  Box degreeBox(IntVect::Zero, m_degree * IntVect::Unit);

  FORT_GETCOARSEFINEINTERPMATRIX(CHF_FRA(m_coarseToFineFab),
                                 CHF_CONST_INT(m_refineCoarse),
                                 CHF_CONST_INT(m_stencilSize),
                                 CHF_CONST_INT(numPowers),
                                 CHF_CONST_INT(m_degree),
                                 CHF_CONST_INT(numFinePoints),
                                 CHF_CONST_VI(m_coarseBaseIndices),
                                 CHF_BOX(m_baseFineBox),
                                 CHF_BOX(degreeBox));

  // The hard part is filling in m_coarseToFineFab.
  // Maybe use Vector<IntVect> m_powers :  see powers2d(degree).
  // Look at coarsetofineconstrainedmatrix
  // and findconstrainedmatrix2d
  // and xpower2dfine0avgall
  // and xpower2dcoarseind0avg
  // and LAPACK for matrix multiplication, transpose?,
  // and of course inverse:
  // DGETRF to get LU decomp, then DGETRI to get inverse.
  // Then convert the results to the closest rational numbers,
  // knowing the factors 256 and 1428|512|212|70|32|10,
  // and check that the rationalized matrix is the inverse,
  // by multiplying integers.

  // Everything is defined now.
  m_defined = true;
}

//////////////////////////////////////////////////////////////////////////////
void FourthOrderInterpStencil::fillFine(
                                        FArrayBox&         a_fineFab,
                                        const FArrayBox&   a_coarseFab,
                                        const IntVect&     a_coarseDataCell,
                                        const IntVect&     a_coarseToFineOffset) const
{
  CH_assert(m_defined);
  IntVect fineBase = m_refineCoarse * a_coarseDataCell + a_coarseToFineOffset;
  // added by petermc, 20 Aug 2009:  fill in a_fineFab only where you can.
  // (previously, had set shiftedFineBox = m_baseFineBox)
  Box shiftedFineBox(a_fineFab.box());
  shiftedFineBox.shift(-fineBase);
  shiftedFineBox &= m_baseFineBox; // intersect with [0:nref-1]^D
  // Fill in a_fineFab at fineBaseBox + fineBase.
  FORT_APPLYCOARSEFINEINTERP(CHF_FRA(a_fineFab),
                             CHF_CONST_FRA(a_coarseFab),
                             CHF_CONST_FRA(m_coarseToFineFab),
                             CHF_CONST_INTVECT(fineBase),
                             CHF_CONST_INTVECT(a_coarseDataCell),
                             CHF_CONST_VI(m_coarseBaseIndices),
                             CHF_BOX(shiftedFineBox));

  // dummy statement in order to get around gdb bug
  int dummy_unused = 0; dummy_unused = 0;
}
