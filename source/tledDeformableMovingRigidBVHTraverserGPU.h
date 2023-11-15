// =========================================================================
// File:       tledDeformableMovingRigidBVHTraverserGPU.h
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
#ifndef tledDeformableMovingRigidBVHTraverserGPU_H
#define tledDeformableMovingRigidBVHTraverserGPU_H

#ifdef _GPU_
#include "tledDeformableRigidBVHTraverserGPU.h"
#include "tledSelfCollisionBVHGPU.h"
#include "tledDynamicBVHGPU.h"

/**
 * \name Interface for deformable/moving rigid body contact search (GPU version)
 * \ingroup contact
 * @{
 */
class tledDeformableMovingRigidBVHTraverserGPU : public tledDeformableRigidBVHTraverserGPU {
  /**
   * \name Construction, Destruction
   * @{
   */
private:
  template <class TBV>
  static tledDeformableMovingRigidBVHTraverserGPU* _CreateTraverserForBVType(tledSelfCollisionBVH &r_defBVH, const tledDynamicBVH &rigBVH);

public:
  static tledDeformableMovingRigidBVHTraverserGPU* CreateTraverser(tledSelfCollisionBVH &r_defBVH, const tledDynamicBVH &rigBVH);

  virtual ~tledDeformableMovingRigidBVHTraverserGPU(void) {}
  /** @} */
};
#endif 
#endif
