// =========================================================================
// File:       tledBV.h
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
#ifndef tledBV_H
#define tledBV_H

#include "tledHelper.h"

#include <algorithm>

/**
 * \brief Basic BV Interface
 */
template <const int t_NumChildBVs>
class tledBV {
  /**
   * \name Hierarchy Information
   * @{
   */
public:
  /** Number of children per BV */
  static const int NumberOfChildBVs = t_NumChildBVs;

  /** BV geometry type ID */
  static const int BVGeometryID = -1;

  /** BV suitable for rotational transforms? */
  static const bool CanRotate = false;

  /**
   * Invalid index entries (e.g. no parent for root BV) are indicated with a value -1.
   */
  int PrimitiveIndex, ParentIndex;

  /**
   * \brief Indices of child BVs. 
   *
   * In higher order hierarchies (> 2), invalid entries (i.e. unoccupied entries in BVs that have less than NumberOfChildBVs children) hold a -1 value.
   */
  int ChildIndices[NumberOfChildBVs];
  /** @} */

  /**
   * \name Computation
   * @{
   */
public:
  virtual void ComputeFromNodeList(const int *nodeListStart, const int *nodeListEnd, const float nodes[]) = 0;
  virtual void ExpandWithNodeList(const int *nodeListStart, const int *nodeListEnd, const float nodes[]) = 0;
  virtual void AddMargin(const float margin) = 0;
  virtual void Translate(const float t[]) = 0;

  /** Copies BV geometries. Should be implemented by every BV type for its respective set of compatible BV classes */
  static void CopyBoundsFromBV(tledBV &r_dstBV, const tledBV &srcBV) { std::abort(); }
  /** Merges bv0 and bv1, storing the result in bv0. Should be implemented by every BV type for its respective set of compatible BV classes */
  static void Merge(tledBV &r_bv0, const tledBV &bv) { std::abort(); }    

  /** Rotates the BV about COR */
  virtual void Rotate(const float m[][3], const float cor[]) {;}
  /** @} */

  /**
   * \name Simple Collision Queries
   * @{
   */
public:
  /** Point-inside-volume test */
  virtual bool IsInside(const float pt[]) const = 0;

  /** Volume-edge intersection test */
  virtual bool DoesIntersect(const float a[], const float b[]) const = 0;

  /** Volume-volume intersection test */
  virtual bool DoesIntersect(const tledBV<t_NumChildBVs> &bv) const = 0;
  /** @} */   

#ifdef _GPU_
  /**
   * \name GPU API
   * @{
   */
public:
  struct GPUBV {
    static const int NumberOfChildBVs = t_NumChildBVs;

    int PrimitiveIndex, ParentIndex;
    int ChildIndices[t_NumChildBVs];    
  };

public:
  virtual void InitGPU(GPUBV &r_dst);
  /** @} */
#endif
};

#include "tledBV.tpp"
#endif
