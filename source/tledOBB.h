// =========================================================================
// File:       tledOBB.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    March 2014
// 
// Copyright (c) 2014, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledOBB_H
#define tledOBB_H

#include "tledBV.h"
#include "tledVectorArithmetic.h"
#include "tledMatrixFunctions.h"
#ifdef _GPU_
#include "tledCUDAHelpers.h"
#endif

#include <limits>
#include <algorithm>
#include <cassert>

/**
 * \brief Oriented bounding boxes
 * \ingroup contact
 */
template <const int t_numChildBVs>
class tledOBB : public tledBV<t_numChildBVs> {
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
  static const bool CanRotate = true;

  /** BV geometry type ID */
  static const int BVGeometryID = 1;
  /** @} */

  /**
   * \name Bounds
   * @{
   */
public:
  float Axes[3][3];
  float Extents[3];
  float Centroid[3];
  /** @} */

  /**
   * \name Misc. Queries
   * @{
   */
public:
  float ComputeVolume(void) const { return 8*this->Extents[0]*this->Extents[1]*this->Extents[2]; }
  /** @} */

  /**
   * \name Simple Collision Queries
   * @{
   */
private:
  bool _DoOverlapOnAxis(const float axis[], const tledOBB &obb) const;
  const float* _ComputeSATAxis(float* computeBuffer, bool &r_isSameAxis, const int axisIndex, const tledOBB &obb) const;
  float _ComputeBoundsOnAxis(const float axis[]) const;

public:
  virtual bool IsInside(const float pt[]) const;
  virtual bool DoesIntersect(const tledBV<t_numChildBVs> &bv) const;
  virtual bool DoesIntersect(const float a[], const float b[]) const;
  /** @} */

  /**
   * \name Computation
   * @{
   */
private:
  void _UpdateExtent(const float x[]);
  void _ResetExtent(void);
  void _ComputeVertices(float (*p_vDst)[3]) const;
  void _CorrectAxisInversion(void);
  void _SwapAxes(const int a0Ind, const int a1Ind);
  void _SortAxes(void);

public:
  virtual void ComputeFromNodeList(const int *nodeListStart, const int *nodeListEnd, const float nodes[]);
  virtual void ExpandWithNodeList(const int *nodeListStart, const int *nodeListEnd, const float nodes[]);
  virtual void AddMargin(const float margin);
  virtual void Translate(const float t[]);

  static void CopyBoundsFromBV(tledOBB &r_dstBV, const tledOBB &srcBV);
  static void Merge(tledOBB &r_bv0, const tledOBB &bv1);
  /** @} */

#ifdef _GPU_
  /**
   * \name GPU API
   * @{
   */
public:
  struct GPUBV : public Superclass::GPUBV {
    float3 Extents;
    float3 Centroid;    
    float3 Axes[3];
  };

public:
  virtual void InitGPU(typename Superclass::GPUBV &r_dst);
  /** @} */
#endif
};

typedef tledOBB<2> tledBinaryOBB;

#include "tledOBB.tpp"
#endif
