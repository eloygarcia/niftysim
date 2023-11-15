// =========================================================================
// File:       tledDeformableDeformableBVHTraverserCPU.h
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
#ifndef tledDeformableDeformableBVHTraverserCPU_H
#define tledDeformableDeformableBVHTraverserCPU_H

#include "tledSelfCollisionBVH.h"
#include "tledBVHTraverserCPU.h"

/**
 * \brief BVH traverser for detecting deformable-deformable contacts (not self-collisions, though).
 * \ingroup contact
 */
class tledDeformableDeformableBVHTraverserCPU : public tledBVHTraverserCPU {  
  /**
   * \name Construction, Destruction
   * @{
   */
public:
  static tledDeformableDeformableBVHTraverserCPU* CreateTraverser(tledSelfCollisionBVH &r_bvh);

  virtual ~tledDeformableDeformableBVHTraverserCPU(void) {}
  /** @} */
};

/**
 * \ingroup contact
 */
template <class TBVH, class TAPI = tledDeformableDeformableBVHTraverserCPU>
class tledDeformableDeformableBVHTraverserImplCPU : public tledBVHTraverserImplCPU<TBVH, TBVH, TAPI> {  
  /**
   * \name Types
   * @{
   */
public:
  typedef tledBVHTraverserImplCPU<TBVH, TBVH, TAPI> Superclass;
  typedef TBVH MasterBVH;
  typedef TBVH SlaveBVH;
  typedef typename TBVH::ContactMesh MasterMesh;
  typedef typename TBVH::ContactMesh SlaveMesh;
  /** @} */

  /**
   * \name Detection
   * @{
   */
protected:
  virtual bool DoNarrowPhaseNodeFacetDetection(const int nodeInd, const int primInd);
  virtual bool DoNarrowPhaseEdgeDetection(const int slaveEdgeInd, const int masterEdgeInd);

  bool DoNarrowPhaseNodeFacetIterativeProjection(float *p_xi, const int nodeInd, const int primInd);

  virtual void ProcessNodeFacetNarrowPhaseItems(const std::vector<std::pair<int, int> >::const_iterator itemsBegin, const std::vector<std::pair<int, int> >::const_iterator itemsEnd);
  virtual void ProcessEdgeEdgeNarrowPhaseItems(const std::vector<std::pair<int, int> >::const_iterator itemsBegin, const std::vector<std::pair<int, int> >::const_iterator itemsEnd);

  virtual void RunBroadPhase(void);

  /** Returns a list of subtree-root BV indices that need to be checked for collisions. */
  virtual const std::vector<int>& GetCollisionCandidateBVs(void) const { return this->GetMasterBVH().GetNonAdjacentGeometryNodes(); }
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
protected:
  tledDeformableDeformableBVHTraverserImplCPU(TBVH &r_slaveBVH, const TBVH &masterBVH) : Superclass(r_slaveBVH, masterBVH) {}

public:
  tledDeformableDeformableBVHTraverserImplCPU(TBVH &r_bvh) : Superclass(r_bvh, r_bvh) {}
  virtual ~tledDeformableDeformableBVHTraverserImplCPU(void) {}
  /** @} */
};

#include "tledDeformableDeformableBVHTraverserCPU.tpp"
#endif
