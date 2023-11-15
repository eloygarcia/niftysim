// =========================================================================
// File:       tledDeformableDeformableBVHTraverserGPU.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    November 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledDeformableDeformableBVHTraverserGPU_H
#define tledDeformableDeformableBVHTraverserGPU_H
#ifndef _GPU_

#error "CUDA file included despite CUDA-support not being enabled at compile time."

#else
#include "tledBVHTraverserGPU.h"
#include "tledCUDAMemoryBlock.h"

/**
 * \brief BVH traverser for detecting collisions between deformable geometry with no topological connection.
 * \ingroup contact
 */
class tledDeformableDeformableBVHTraverserGPU : public tledBVHTraverserGPU {
  /**
   * \name Construction, Destruction
   * @{
   */
public:
  static tledDeformableDeformableBVHTraverserGPU* CreateTraverser(tledSelfCollisionBVH &r_bvh);
  
  virtual ~tledDeformableDeformableBVHTraverserGPU(void) {}
  /** @} */
};

#endif
#endif
