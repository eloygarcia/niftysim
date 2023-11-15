// =========================================================================
// File:       tledParallelDeformableDeformableBVHTraverserCPU.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    January 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledParallelDeformableDeformableBVHTraverserCPU_H
#define tledParallelDeformableDeformableBVHTraverserCPU_H

#include "tledParallelBVHTraverserCPU.h"
#include "tledDeformableDeformableBVHTraverserCPU.h"

#include <algorithm>

#ifdef _USE_BOOST_

/**
 * \ingroup contact
 */
template <class TBaseTraverser>
class tledParallelDeformableDeformableBVHTraverserCPU : public tledParallelBVHTraverserImplCPU<TBaseTraverser> {
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
  class _NodeFacetNarrowPhaseSubTraverser;
  class _EdgeEdgeNarrowPhaseSubTraverser;

protected:  
  virtual void RunBroadPhase(void);

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
  tledParallelDeformableDeformableBVHTraverserCPU(SlaveBVH &r_bvh) : Superclass(r_bvh, r_bvh) {}
  virtual ~tledParallelDeformableDeformableBVHTraverserCPU(void) {}
  /** @} */
};

#include "tledParallelDeformableDeformableBVHTraverserCPU.tpp"
#endif
#endif
