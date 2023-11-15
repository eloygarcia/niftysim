// =========================================================================
// File:       tledDeformableMovingRigidBVHTraverserCPU.h
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
#ifndef tledDeformableMovingRigidBVHTraverserCPU_H
#define tledDeformableMovingRigidBVHTraverserCPU_H

#include "tledDeformableRigidBVHTraverserCPU.h"
#include "tledDynamicBVH.h"

/**
 * \name Interface for deformable/moving rigid body contact search (CPU version)
 * \ingroup contact
 * @{
 */
class tledDeformableMovingRigidBVHTraverserCPU : public tledDeformableRigidBVHTraverserCPU {
  /**
   * \name Construction, Destruction
   * @{
   */
private:
  template <class TBV>
  static tledDeformableMovingRigidBVHTraverserCPU* _CreateTraverserForBVType(tledSelfCollisionBVH &r_defBVH, const tledDynamicBVH &rigBVH);

public:
  static tledDeformableMovingRigidBVHTraverserCPU* CreateTraverser(tledSelfCollisionBVH &r_defBVH, const tledDynamicBVH &rigBVH);

  virtual ~tledDeformableMovingRigidBVHTraverserCPU(void) {}
  /** @} */
};

template <class TDeformableBVH, class TRigidBVH, class TAPI = tledDeformableMovingRigidBVHTraverserCPU>
class tledDeformableMovingRigidBVHTraverserImplCPU : public tledDeformableRigidBVHTraverserImplCPU<TDeformableBVH, TRigidBVH, TAPI> {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledDeformableRigidBVHTraverserImplCPU<TDeformableBVH, TRigidBVH, TAPI> Superclass;
  typedef TRigidBVH MasterBVH;
  typedef TDeformableBVH SlaveBVH;
  typedef typename MasterBVH::ContactMesh MasterMesh;
  typedef typename SlaveBVH::ContactMesh SlaveMesh;
  /** @} */

  /**
   * \name Detection
   * @{
   */
protected:
  virtual bool DoNarrowPhaseNodeFacetDefMaster(const int defPrimIndex, const int rigNodeIndex);
  virtual bool DoNarrowPhaseEdgeEdgeDefMaster(const int defEdgeIndex, const int rigEdgeIndex);
  virtual bool DoNarrowPhaseNodeFacetRigMaster(const int defNodeIndex, const int rigFacetIndex);
  virtual bool DoNarrowPhaseEdgeEdgeRigMaster(const int defEdgeIndex, const int rigEdgeIndex);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledDeformableMovingRigidBVHTraverserImplCPU(SlaveBVH &r_slaveBVH, const MasterBVH &masterBVH) : Superclass(r_slaveBVH, masterBVH) {}
  virtual ~tledDeformableMovingRigidBVHTraverserImplCPU(void) {}
  /** @} */  
};

#include "tledDeformableMovingRigidBVHTraverserCPU.tpp"

#endif
