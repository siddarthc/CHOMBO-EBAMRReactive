#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBLGINTEGRATOR_H_
#define _EBLGINTEGRATOR_H_
      
#if (CH_SPACEDIM == 2)

#define QRHO  0
#define QVELX 1
#define QVELY 2
#define QPRES 3
#define QENTR 4
#define QC 5
#define QTEMP 6
#ifdef MODIANO_PROBLEM
#define QMVAX 7
#define QMVTN 8
#define QNUM  9
#define QSLOPE 9
#define QSPEC1 9
#else
#define QNUM 7
#define QSLOPE 7
#define QSPEC1 7
#endif
#define CRHO  0
#define CMOMX 1
#define CMOMY 2
#define CENG  3
#define CNUM  4
#define CSPEC1 4

#ifdef RZ_COORDS
#define FNUM  5
#define FSPEC1 5
#define CPRES 4
#else
#define FNUM  4
#define FSPEC1 4
#define CPRES 0
#endif

#endif

#if (CH_SPACEDIM == 3)
#define CPRES 0

#define QRHO  0
#define QVELX 1
#define QVELY 2
#define QVELZ 3
#define QPRES 4
#define QENTR 5
#define QC    6
#define QTEMP 7
#ifdef MODIANO_PROBLEM
#define QMVAX 8
#define QMVTN1 9
#define QMVTN2 10
#define QNUM  11
#define QSLOPE 11
#define QSPEC1 11
#else
#define QNUM 8
#define QSLOPE 8
#define QSPEC1 8
#endif

#define CRHO  0
#define CMOMX 1
#define CMOMY 2
#define CMOMZ 3
#define CENG  4
#define CNUM  5
#define FNUM  5
#define CSPEC1 5
#define FSPEC1 5
#endif

#endif

