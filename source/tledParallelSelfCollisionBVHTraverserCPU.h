// =========================================================================
// File:       tledParallelSelfCollisionBVHTraverser.h
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
#ifndef tledParallelSelfCollisionBVHTraverserCPU_H
#define tledParallelSelfCollisionBVHTraverserCPU_H
#include "tledParallelDeformableDeformableBVHTraverserCPU.h"

/**
 * \brief Adds checks for primitive adjacency.
 * \ingroup contact
 */
template <class TBaseTraverser>
class tledParallelSelfCollisionBVHTraverserCPU : public tledParallelDeformableDeformableBVHTraverserCPU<TBaseTraverser> {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledParallelDeformableDeformableBVHTraverserCPU<TBaseTraverser> Superclass;
  typedef typename Superclass::MasterBVH MasterBVH;
  typedef typename Superclass::SlaveBVH SlaveBVH;
  typedef typename Superclass::MasterMesh MasterMesh;
  typedef typename Superclass::SlaveMesh SlaveMesh;
  typedef typename Superclass::NodeFacetNarrowPhaseSubTraverser NodeFacetNarrowPhaseSubTraverser;
  typedef typename Superclass::EdgeEdgeNarrowPhaseSubTraverser EdgeEdgeNarrowPhaseSubTraverser;
  /** @} */

  /**
   * \name Detection
   * @{
   */
private:
  class _BroadPhaseSubTraverser;

protected:
  virtual typename tledParallelBVHTraverserImplCPU<TBaseTraverser>::BroadPhaseSubTraverser* LaunchSubTraverser(const int slaveStartInd, const int masterStartInd);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledParallelSelfCollisionBVHTraverserCPU(SlaveBVH &r_bvh) : Superclass(r_bvh) {}
  virtual ~tledParallelSelfCollisionBVHTraverserCPU(void) {}
  /** @} */
};

#include "tledParallelSelfCollisionBVHTraverserCPU.tpp"
#endif
