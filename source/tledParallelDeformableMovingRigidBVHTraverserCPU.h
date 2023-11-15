// =========================================================================
// File:       tledParallelDeformableMovingRigidBVHTraverserCPU.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    March 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledParallelDeformableMovingRigidBVHTraverserCPU_H
#define tledParallelDeformableMovingRigidBVHTraverserCPU_H

#include "tledParallelDeformableRigidBVHTraverserCPU.h"

#ifdef _USE_BOOST_
template <class TBaseTraverser>
class tledParallelDeformableMovingRigidBVHTraverserCPU : public tledParallelDeformableRigidBVHTraverserCPU<TBaseTraverser> {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledParallelDeformableRigidBVHTraverserCPU<TBaseTraverser> Superclass;
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
  class _NodeFacetNarrowPhaseDefMasterSubTraverser;
  class _NodeFacetNarrowPhaseRigMasterSubTraverser;
  class _EdgeEdgeNarrowPhaseDefMasterSubTraverser;
  class _EdgeEdgeNarrowPhaseRigMasterSubTraverser;

protected:  
  virtual NodeFacetNarrowPhaseSubTraverser* LaunchNodeFacetNarrowPhaseWorker(const std::vector<std::pair<int, int> >::const_iterator itemsBegin, const std::vector<std::pair<int, int> >::const_iterator itemsEnd, const int threadInd);
  virtual EdgeEdgeNarrowPhaseSubTraverser* LaunchEdgeEdgeNarrowPhaseWorker(const std::vector<std::pair<int, int> >::const_iterator itemsBegin, const std::vector<std::pair<int, int> >::const_iterator itemsEnd, const int threadInd);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledParallelDeformableMovingRigidBVHTraverserCPU(SlaveBVH &r_bvh, const MasterBVH &masterBVH) : Superclass(r_bvh, masterBVH) {}
  virtual ~tledParallelDeformableMovingRigidBVHTraverserCPU(void) {}
  /** @} */
};

#include "tledParallelDeformableMovingRigidBVHTraverserCPU.tpp"
#endif
#endif
