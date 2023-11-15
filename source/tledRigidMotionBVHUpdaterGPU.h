// =========================================================================
// File:       tledRigidMotionBVHUpdaterGPU.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   CUDA
// Created:    September 2014
// 
// Copyright (c) 2014, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledRigidMotionBVHUpdaterGPU_H
#define tledRigidMotionBVHUpdaterGPU_H

#include "tledRigidMotionBVHUpdater.h"
#include "tledGreedySelfCollisionBVHUpdaterGPU.h"

template <class TBVH>
class tledRigidMotionBVHUpdaterGPU : public tledRigidMotionBVHUpdater<TBVH> {
  /**
   * \name Types
   * @{
   */
public:
  typedef TBVH BVH;
  typedef typename BVH::ContactMesh ContactMesh;
  typedef tledRigidMotionBVHUpdater<TBVH> Superclass;
  /** @} */

  /**
   * \name Updating
   * @{
   */
private:
  void _TranslateBVH(const float tInc[]);

public:
  virtual void UpdateBVH(void);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledRigidMotionBVHUpdaterGPU(BVH &r_bvh) : Superclass(r_bvh) {}
  virtual ~tledRigidMotionBVHUpdaterGPU(void) {}
  /** @} */
};

#include "tledRigidMotionBVHUpdaterGPU.tpp"
#endif
