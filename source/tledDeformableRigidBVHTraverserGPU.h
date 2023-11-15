// =========================================================================
// File:       tledDeformableRigidBVHTraverserGPU.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   CUDA
// Created:    July 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledDeformableRigidBVHTraverserGPU_H
#define tledDeformableRigidBVHTraverserGPU_H

#ifdef _GPU_
#include "tledBVHTraverserGPU.h"
#include "tledSelfCollisionBVHGPU.h"
#include "tledStaticBVHGPU.h"

/**
 * \brief BVH traverser for detecting collisions between deformable and rigid geometry.
 * \ingroup contact
 */
class tledDeformableRigidBVHTraverserGPU : public tledBVHTraverserGPU {
  /**
   * \name Construction, Destruction
   * @{
   */
public:
  static tledDeformableRigidBVHTraverserGPU* CreateTraverser(tledSelfCollisionBVHGPU &r_bvh, const tledBVH &rigBVH);
  
  virtual ~tledDeformableRigidBVHTraverserGPU(void) {}
  /** @} */
};

#endif
#endif
