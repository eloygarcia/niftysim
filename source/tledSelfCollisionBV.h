// =========================================================================
// File:       tledSelfCollisionBV.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    May 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledSelfCollisionBV_H
#define tledSelfCollisionBV_H

#include "tledBV.h"
#include "tledAABB.h"
#include "tledOBB.h"
#ifdef _GPU_
#include "tledCUDAHelpers.h"
#endif

/**
 * \brief Adds self-collision detection related members to any type of BV 
 * \ingroup contact
 */
template <class TBaseBV>
struct tledSelfCollisionBV : public TBaseBV {    
  /**
   * \name Types
   * @{
   */
public:
  typedef TBaseBV BaseBV;
  typedef BaseBV Superclass;
  /** @} */

  /**
   * \name Self-collision detection and lazy updating facilities
   * @{
   */
public:
  int UpdateCounter;

  /** VolinoAngle: Surface cone opening angle. If value = -1 self-collision criterion not applicable (probably because bounded geometry not connected) */
  float VolinoAngle;

  /** Surface cone axis */
  float VolinoAxis[3];

  /** Minimum element diameter in bounded geometry */
  float SubtreeMinH;
  /** @} */

#ifdef _GPU_
  /**
   * \name GPU API
   * @{
   */
public:
  typedef typename TBaseBV::GPUBV GPUBV;
#endif

  /**
   * \name VMT Angle Bounds
   * @{
   */
public:  
  static float GetVolinoThresholdAngle(void) { return 0.95f*tledPi; }
  static float GetConeMinAngle(void) { return 0.01f*tledPi/180.0f; }
  /** @} */

  /**
   * \name Computation
   * @{
   */
public:
  static void CopyBoundsFromBV(tledSelfCollisionBV &r_dstBV, const tledSelfCollisionBV &srcBV) { TBaseBV::CopyBoundsFromBV(r_dstBV, srcBV); }
  static void Merge(tledSelfCollisionBV &r_bv0, const tledSelfCollisionBV &bv1) { TBaseBV::Merge(r_bv0, bv1); }
  /** @} */

  /**
   * \name Geometry Queries
   * @{
   */
public:
  void SetDisconnectedGeometryFlag(void) { VolinoAngle = -1.f; }

  /** If true, the surface cone self-collision criterion is not applicable. */
  bool HasDisconnectedGeometry(void) const { return VolinoAngle == -1.f; }
  /** @} */

  tledSelfCollisionBV() {
#ifndef NDEBUG
    std::fill(VolinoAxis, VolinoAxis + 3, std::numeric_limits<float>::quiet_NaN());
#endif
    VolinoAngle = SubtreeMinH = std::numeric_limits<float>::quiet_NaN();
    UpdateCounter = -1;
  }
};

typedef tledSelfCollisionBV<tledAABB<2> > tledBinarySelfCollisionAABB;
typedef tledSelfCollisionBV<tledOBB<2> > tledBinarySelfCollisionOBB;

typedef tledSelfCollisionBV<tledAABB<4> > tledQuarternarySelfCollisionAABB;
typedef tledSelfCollisionBV<tledOBB<4> > tledQuarternarySelfCollisionOBB;

#endif
