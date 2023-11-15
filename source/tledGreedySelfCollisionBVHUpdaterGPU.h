// =========================================================================
// File:       tledGreedySelfCollisionBVHUpdaterGPU.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    April 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledGreedySelfCollisionBVHUpdaterGPU_H
#define tledGreedySelfCollisionBVHUpdaterGPU_H

#include "tledBottomUpBVHUpdaterGPU.h"

template <class TBVH>
class tledGreedySelfCollisionBVHUpdaterGPU : public tledBottomUpBVHUpdaterGPU<TBVH> {
  /**
   * \name Types
   * @{
   */
public:
  typedef TBVH BVH;
  typedef typename TBVH::BoundingVolume BoundingVolume;
  typedef typename TBVH::ContactMesh ContactMesh;
  typedef typename ContactMesh::Facet Facet;
  typedef tledBottomUpBVHUpdaterGPU<TBVH> Superclass;
  /** @} */

  /**
   * \name Initialisation
   * @{
   */
public:
  virtual void Init(void);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledGreedySelfCollisionBVHUpdaterGPU(BVH &r_bvh) : Superclass(r_bvh) {}
  tledGreedySelfCollisionBVHUpdaterGPU(void) {}
  virtual ~tledGreedySelfCollisionBVHUpdaterGPU(void) {}
  /** @} */
};

#if defined __CUDACC__ && !defined __GPU_TEST_LINK_CONFLICT_NO_INCLUDE
#include "tledGreedySelfCollisionBVHUpdaterGPU.tpp"
#endif

#endif
