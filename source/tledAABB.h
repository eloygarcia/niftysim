// =========================================================================
// File:       tledAABB.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    February 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledAABB_H
#define tledAABB_H

#include "tledBV.h"
#include "tledVectorArithmetic.h"
#ifdef _GPU_
#include "tledCUDAHelpers.h"
#endif

#include <limits>
#include <algorithm>
#include <cassert>

/**
 * \brief AABB data-structure for AABB hierarchies
 * \ingroup contact
 */
template <const int t_numChildBVs>
class tledAABB : public tledBV<t_numChildBVs> {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledBV<t_numChildBVs> Superclass;
  /** @} */

  /**
   * \name Traits
   * @{
   */
public:
  static const bool CanRotate = false;

  /** BV geometry type ID */
  static const int BVGeometryID = 0;
  /** @} */

  /**
   * \name Bounds
   * @{
   */
public:
  /**
   * AABB bounds:
   * First index
   * 0 -> x
   * 1 -> y
   * 2 -> z
   * 2nd index:
   * 0 -> lower bound
   * 1 -> upper bound
   */
  float Bounds[3][2];
  /** @} */

  /**
   * \name Misc. Queries
   * @{
   */
public:
  float ComputeVolume(void) const { return (this->Bounds[2][1] - this->Bounds[2][0])*(this->Bounds[1][1] - this->Bounds[1][0])*(this->Bounds[0][1] - this->Bounds[0][0]); }
  /** @} */

  /**
   * \name Simple Collision Queries
   * @{
   */
public:
  virtual bool IsInside(const float pt[]) const {
    return pt[0] > this->Bounds[0][0] && pt[0] < this->Bounds[0][1] 
      && pt[1] > this->Bounds[1][0] && pt[1] < this->Bounds[1][1]
      && pt[2] > this->Bounds[2][0] && pt[2] < this->Bounds[2][1];
  }

  virtual bool DoesIntersect(const tledBV<t_numChildBVs> &bv) const;
  virtual bool DoesIntersect(const float a[], const float b[]) const;
  /** @} */

  /**
   * \name Computation
   * @{
   */
public:
  virtual void ComputeFromNodeList(const int *nodeListStart, const int *nodeListEnd, const float nodes[]);
  virtual void ExpandWithNodeList(const int *nodeListStart, const int *nodeListEnd, const float nodes[]);
  virtual void AddMargin(const float margin);
  virtual void Translate(const float t[]);

  static void CopyBoundsFromBV(tledAABB &r_dstBV, const tledAABB &srcBV);
  static void Merge(tledAABB &r_bv0, const tledAABB &bv1);
  /** @} */

#ifdef _GPU_
  /**
   * \name GPU API
   * @{
   */
public:
  struct GPUBV : public Superclass::GPUBV {
    float2 Bounds[3];
  };

public:
  virtual void InitGPU(typename Superclass::GPUBV &r_dst);
  /** @} */
#endif
}; 

typedef tledAABB<2> tledBinaryAABB;

#include "tledAABB.tpp"

#endif
