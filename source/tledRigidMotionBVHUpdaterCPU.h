// =========================================================================
// File:       tledRigidMotionBVHUpdaterCPU.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    September 2014
// 
// Copyright (c) 2014, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledRigidMotionBVHUpdaterCPU_H
#define tledRigidMotionBVHUpdaterCPU_H

#include "tledRigidMotionBVHUpdater.h"

template <class TBVH>
class tledRigidMotionBVHUpdaterCPU : public tledRigidMotionBVHUpdater<TBVH> {
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
public:
  virtual void UpdateBVH(void);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledRigidMotionBVHUpdaterCPU(BVH &r_bvh) : Superclass(r_bvh) {}
  tledRigidMotionBVHUpdaterCPU(void) {}
  virtual ~tledRigidMotionBVHUpdaterCPU(void) {}
  /** @} */
};

#include "tledRigidMotionBVHUpdaterCPU.tpp"

#endif
