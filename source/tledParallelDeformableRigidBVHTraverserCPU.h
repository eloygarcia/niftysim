// =========================================================================
// File:       tledParallelDeformableRigidBVHTraverserCPU.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    February 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledParallelDeformableRigidBVHTraverserCPU_H
#define tledParallelDeformableRigidBVHTraverserCPU_H

#include "tledParallelBVHTraverserCPU.h"
#include "tledDeformableRigidBVHTraverserCPU.h"

#include <vector>
#include <algorithm>
#include <cassert>

#ifdef _USE_BOOST_

/**
 * \brief CPU-parallel BVH traversal for rigid-deformable contact search
 * \ingroup contact
 */
template <class TBaseTraverser>
class tledParallelDeformableRigidBVHTraverserCPU : public tledParallelBVHTraverserImplCPU<TBaseTraverser> {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledParallelBVHTraverserImplCPU<TBaseTraverser> Superclass;
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
  class BroadPhaseSubTraverser;

protected:  
  virtual void RunBroadPhase(void);

  virtual typename Superclass::BroadPhaseSubTraverser* LaunchSubTraverser(const int slaveStartInd, const int masterStartInd);
  virtual NodeFacetNarrowPhaseSubTraverser* LaunchNodeFacetNarrowPhaseWorker(const std::vector<std::pair<int, int> >::const_iterator itemsBegin, const std::vector<std::pair<int, int> >::const_iterator itemsEnd, const int threadInd);
  virtual EdgeEdgeNarrowPhaseSubTraverser* LaunchEdgeEdgeNarrowPhaseWorker(const std::vector<std::pair<int, int> >::const_iterator itemsBegin, const std::vector<std::pair<int, int> >::const_iterator itemsEnd, const int threadInd);

  virtual void PreNodeFacetNarrowPhaseHook(void);
  virtual void PreEdgeEdgeNarrowPhaseHook(void);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledParallelDeformableRigidBVHTraverserCPU(SlaveBVH &r_bvh, const MasterBVH &masterBVH) : Superclass(r_bvh, masterBVH) {}
  virtual ~tledParallelDeformableRigidBVHTraverserCPU(void) {}
  /** @} */
};

#include "tledParallelDeformableRigidBVHTraverserCPU.tpp"
#endif
#endif
