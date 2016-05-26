#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _BASEMIXBCVALUE_H_
#define _BASEMIXBCVALUE_H_

#include "RealVect.H"
#include "EBArith.H"
#include "NamespaceHeader.H"

class BaseMixBCValue
{
public:
  BaseMixBCValue()
  {
  }

  virtual ~BaseMixBCValue()
  {
  }

  virtual bool isDirichlet(const RealVect& a_point) const = 0;
  virtual bool isDirichlet(const RealVect& a_point, const Real& a_dx) const 
  {
    return isDirichlet(a_point);
  };

  virtual Real valueDir(const RealVect& a_point,
                     const RealVect& a_normal,
                     const Real&     a_time,
                     const int&      a_comp) const = 0;

  virtual Real valueDir(const FaceIndex&      a_face,
                     const Side::LoHiSide& a_side,
                     const DataIndex&      a_dit,
                     const RealVect&       a_point,
                     const RealVect&       a_normal,
                     const Real&           a_time,
                     const int&            a_comp) const
  {
    return valueDir(a_point,a_normal,a_time,a_comp);
  }

  virtual Real valueNeu(const RealVect& a_point,
                     const RealVect& a_normal,
                     const Real&     a_time,
                     const int&      a_comp) const = 0;

  virtual Real valueNeu(const FaceIndex&      a_face,
                     const Side::LoHiSide& a_side,
                     const DataIndex&      a_dit,
                     const RealVect&       a_point,
                     const RealVect&       a_normal,
                     const Real&           a_time,
                     const int&            a_comp) const
  {
    return valueNeu(a_point,a_normal,a_time,a_comp);
  }

   virtual Real valueNeu(const IntVect&   a_iv,
                     const DataIndex& a_dit,
                     const RealVect&  a_point,
                     const RealVect&  a_normal,
                     const Real&      a_time,
                     const int&       a_comp) const
  {
    return valueNeu(a_point,a_normal,a_time,a_comp);
  }

  virtual Real valueNeu(const VolIndex&  a_vof,
                     const RealVect&  a_centroid,
                     const RealVect&  a_normal,
                     const RealVect&  a_dx,
                     const RealVect&  a_probLo,
                     const DataIndex& a_dit,
                     const Real&      a_time,
                     const int&       a_comp) const
  {
    RealVect point = a_centroid;
    point *= a_dx;
    point += EBArith::getVofLocation(a_vof,a_dx,a_probLo);
    return valueNeu(point,a_normal,a_time,a_comp);
  }


};

#include "NamespaceFooter.H"
#endif
