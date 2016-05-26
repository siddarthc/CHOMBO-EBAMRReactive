#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _FOURTHORDERINTERPSTENCIL_H_
#define _FOURTHORDERINTERPSTENCIL_H_

#include "FArrayBox.H"
#include "IntVect.H"
#include "Vector.H"
#include "UsingNamespace.H"

/// Time interpolator class using 4th-order Runge-Kutta

/**
 */
class FourthOrderInterpStencil
{
public:
  /// Default constructor
  /**
     Object requires define() to be called before all other functions.
   */
  FourthOrderInterpStencil();

  /// Full constructor
  FourthOrderInterpStencil(/// in range -(a_degree-1):(a_degree-1), specifies displacement to boundary in each dimension
                           const IntVect&         a_bdryOffset,
                           /// refinement ratio
                           const int&             a_refineCoarse,
                           /// max degree of polynomials
                           const int&             a_degree = 3)
  {
    define(a_bdryOffset, a_refineCoarse, a_degree);
  }

  /// Destructor
  /**
     Destroys all objects created by define(). Passed in data references
     of define() are left alone.
   */
  ~FourthOrderInterpStencil();

  /// Actual constructor.
  /**
     Set up object.
   */
  void define(/// in range -(a_degree-1):(a_degree-1), specifies displacement to boundary in each dimension
              const IntVect&         a_bdryOffset,
              /// refinement ratio
              const int&             a_refineCoarse,
              /// max degree of polynomials
              const int&             a_degree = 3);

  /// Interpolate to all the fine cells within one coarse cell.
  /**
     Fill the fine cells inside a_coarseDataCell with interpolated data
     from a_coarseFab.

     The fine cells to be filled in are at
     m_baseFineBox + m_refineCoarse * (a_coarseDataCell + a_coarseToFineOffset).
     The coarse cells from which we draw data are
     (a_coarseDataCell + m_coarseBaseIndices[i*D+[0:2]]:  0 <= i < m_stencilSize).
   */
  void fillFine(/// to fill at m_baseFineBox + m_refineCoarse*(a_coarseBase + a_coarseToFineOffset)
                FArrayBox&         a_fineFab,
                /// coarse data
                const FArrayBox&   a_coarseFab,
                /// coarse cell with fine subcells to fill in
                const IntVect&     a_coarseDataCell,
                /// offset from base coarse cell to coarsened fine cells
                const IntVect&     a_coarseToFineOffset = IntVect::Zero) const;

protected:

  /// refinement ratio
  int m_refineCoarse;

  /// max degree of polynomials
  int m_degree;

  /// displacement from domain boundary in each dimension, in range -(a_degree-1):(a_degree-1)
  IntVect m_bdryOffset;

  /// 0:(m_refineCoarse-1) in every dimension
  Box m_baseFineBox;

  /// number of coarse cells in the stencil
  int m_stencilSize;

  /// length m_stencilSize*SpaceDim, lists indices of coarse cell corresponding to each component of m_coarseToFineFab
  Vector<int> m_coarseBaseIndices;

  /// m_stencilSize coefficients of coarse cells for each cell in m_baseFineBox
  FArrayBox m_coarseToFineFab;

  /// whether define() has been called
  bool m_defined;

private:

  // Disallowed for all the usual reasons
  void operator=(const FourthOrderInterpStencil& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  FourthOrderInterpStencil(const FourthOrderInterpStencil& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif
