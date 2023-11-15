// =========================================================================
// File:       tledTrackingDeformableDeformableBVHTraverserGPU.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    December 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledTrackingDeformableDeformableBVHTraverserGPU_H
#define tledTrackingDeformableDeformableBVHTraverserGPU_H

#include "tledTrackingBVHTraverserGPU.h"
#include "tledDeformableDeformableBVHTraverserGPU.h"

/**
 * \ingroup contact
 */
class tledTrackingDeformableDeformableBVHTraverserGPU : public tledDeformableDeformableBVHTraverserGPU, public tledTrackingBVHTraverser {
  /**
   * \name Construction, Destruction
   * @{
   */
public:
  static tledDeformableDeformableBVHTraverserGPU* CreateTraverser(tledSelfCollisionBVH &r_bvh);

  virtual ~tledTrackingDeformableDeformableBVHTraverserGPU(void) {}
  /** @} */
};

/**
 * \brief Deformable-deformable collision detection with tracking.
 * \ingroup contact
 */
template <class TBaseDeformableDeformableBVHTraverserGPU>
class tledTrackingDeformableDeformableBVHTraverserImplGPU : public tledTrackingBVHTraverserImplGPU<TBaseDeformableDeformableBVHTraverserGPU> {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledTrackingBVHTraverserImplGPU<TBaseDeformableDeformableBVHTraverserGPU> Superclass;
  typedef typename TBaseDeformableDeformableBVHTraverserGPU::MasterBVH MasterBVH;
  typedef typename TBaseDeformableDeformableBVHTraverserGPU::SlaveBVH SlaveBVH;
  typedef typename MasterBVH::BoundingVolume BoundingVolume;
  typedef typename MasterBVH::ContactMesh MasterMesh;
  typedef typename MasterBVH::ContactMesh SlaveMesh;
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  virtual void Init(tledUnstructuredContactManager &r_manager);

  tledTrackingDeformableDeformableBVHTraverserImplGPU(SlaveBVH &r_slaveBVH, const MasterBVH &masterBVH) : Superclass(r_slaveBVH, masterBVH) {}
  tledTrackingDeformableDeformableBVHTraverserImplGPU(SlaveBVH &r_bvh) : Superclass(r_bvh, r_bvh) {}

  virtual ~tledTrackingDeformableDeformableBVHTraverserImplGPU(void);
  /** @} */
};

#endif
