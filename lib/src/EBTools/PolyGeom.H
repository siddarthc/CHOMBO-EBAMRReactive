#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//  ANAG, LBNL

#ifndef _POLYGEOM_H_
#define _POLYGEOM_H_

#include "REAL.H"
#include "RealVect.H"
#include "Tuple.H"
#include "Vector.H"

#include "VolIndex.H"

#include "NamespaceHeader.H"

class EBISBox;
///
/**
   PolyGeom encapsulates the functions
   to generate centroids, normals, etc given
   the geometric information.
   It is meant to be used simply as a static class.
   Its only data member (a static) is its tolerance.
*/
class PolyGeom
{
public:

  ///
  /**
     return the inverse of the input matrix.   done using cramers rule
   */
  static void invertMatrix(Real      a_AInverse[SpaceDim][SpaceDim],
                           const Real       a_A[SpaceDim][SpaceDim], bool a_test);

  ///
  /**
   // calculate the cofactor of element (row,col)
   //this is the matrix which excludes row and column row, col
   */
  static void getMinor(Real      a_Aminor[SpaceDim-1][SpaceDim-1],
                       const Real     a_A[SpaceDim  ][SpaceDim  ],
                       int a_row, int acol);

  ///
  static Real determinantSD(const Real a_A[SpaceDim][SpaceDim]);

  ///
  static Real determinantSDM1(const Real a_A[SpaceDim-1][SpaceDim-1]);

  ///
  /**
     find equation of line
     between a line formed by a_pointOnLine = p0 and a_direction=v
     l(t) = p0 + v t
     and a point in space a_point.
     returns closest point on line and the normal between the
     closest point and the point in space
   */
  static void
  pointToLine(RealVect& a_closestPt,
              RealVect& a_normal,
              const RealVect& a_point,
              const RealVect& a_pointOnLine,
              const RealVect& a_direction);

  ///return distance from point to a plane
  static Real
  distancePointToPlane(const RealVect& a_point,
                       const RealVect& a_normal,
                       const RealVect& a_pointOnLine);

  ///
  /**
   */
  static void setTolerance(const Real& a_tolerance);

  ///
  /**
   */
  static void setVectDx(const RealVect& a_vectDx);

  ///
  /**
   */
  static void setVolumeTolerance(const Real& a_tolerance);

  ///
  /**
   */
  static void setAreaTolerance(const Real& a_tolerance);

  ///
  /**
   */
  static void setLengthTolerance(const Real& a_tolerance);

  ///
  /**
     compute the cross product between xvec0 and xvec1
     (returns xvec1 x xvec0)
  */
  static RealVect cross(const RealVect& a_xvec1,
                        const RealVect& a_xvec0);

  ///
  /**
     compute the dot product between xvec0 and xvec1
     (returns xvec1 dot xvec0)
  */
  static Real dot(const RealVect& a_xvec1,
                  const RealVect& a_xvec0);

  ///
  /**
   */
  static const Real& getTolerance();

  ///
  /**
   */
  static const RealVect& getVectDx();

  ///
  /**
   */
  static const Real& getVolumeTolerance();

  ///
  /**
   */
  static const Real& getAreaTolerance();

  ///
  /**
   */
  static const Real& getLengthTolerance();

  ///
  /**
     Return the normal to the boundary at the
     input VoF.  If said normal is undefined,
     returns the zero vector.
  */
  static RealVect normal(const VolIndex& a_vof,
                         const EBISBox& a_ebisBox,
                         const Real& a_bndryArea);

  ///
  /**
   */
  static Real bndryArea(const VolIndex& a_vof,
                        const EBISBox& a_ebisBox);

  ///
  /**
     Return the inhomogeneous component of the plane
     eqution ni xi = alpha.
     Volume fraction must be positive and <= 1.0
     and the normals can be in any order.
  */
  static Real computeAlpha(const Real& a_volFrac,
                           const RealVect& a_normal);

  ///
  /**
     Compute the volume fraction in the cell given
     alpha and the normal.
  */
  static Real computeVolume(const Real& a_alpha,
                            const RealVect& a_normal);

  ///
  /**
     Get the normal and alpha given the points in a polygon.
     Solves the linear system \\
     ni xi - alpha = 0\\
     \sum ni*ni = 1\\
     for ni and alpha.
     Always returns nromal[updir] >= 0
  */
  static void computeNormalAndAlpha(Real& a_alpha,
                                    RealVect& a_normal,
                                    const int& a_upDir,
                                    const Tuple<RealVect, CH_SPACEDIM>& a_poly);

  // solve ax = b (gives icomp component of x using cramers rule).
  static Real
  matrixSolveComp(const Vector<Vector<Real> >& a_A,
                  const  Vector<Real>& a_rhs,
                  const int& a_icomp);

  //return the determinant of a
  static Real
  determinant(const Vector< Vector<Real> >& a_A);

  /**
     andzvolume computes the volume under
     the curve (covered volume).
     Does so using Scardovelli
     and Zaleski's direct formulae.
     The normal has to be all positive numbers.
     The normals also must be ordered from lowest
     to highest.
  */
  static Real sAndZVolume(const Real& a_alpha,
                          const RealVect& a_normal);

  ///
  PolyGeom()
  {;}
  ///
  ~PolyGeom()
  {;}

  ///
  /**
     Sort vector from low to high.  All comps must be positive.
     Return reverse mapping
  */
  static void sortVector(RealVect& vect, IntVect& ivmap);

  ///make vector all pos and return the signs
  static void posifyVector(RealVect& vect, IntVect& signs);

  ///make vector into a unit vector and return the sum of squares
  static void unifyVector(RealVect& normal, Real& sumSquare);

  ///
  static Tuple<int, CH_SPACEDIM-1> computeTanDirs(int upDir);

protected:

  static Real   twoDFunc(const Real& arg);
  static Real threeDFunc(const Real& arg);

  static Real computeAnyVolume(const Real& a_alpha,
                               const Real& a_norm0,
                               const Real& a_norm1,
                               const Real& a_norm2);

  static RealVect s_vectDx;

  static Real s_tolerance;

  static Real s_lengthtolerance;

  static Real s_areatolerance;

  static RealVect tetCentroid(const RealVect& normal,
                              const Real& alpha);

  static Real tetVolume(const RealVect& normal,
                        const Real& alpha);


private:
};

#include "NamespaceFooter.H"
#endif
