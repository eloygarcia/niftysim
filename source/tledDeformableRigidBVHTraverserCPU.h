// =========================================================================
// File:       tledDeformableRigidBVHTraverserCPU.h
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
#ifndef tledDeformableRigidBVHTraverserCPU_H
#define tledDeformableRigidBVHTraverserCPU_H

#include "tledHelper.h"
#include "tledSelfCollisionBVH.h"
#include "tledVectorArithmetic.h"
#include "tledBVHTraverserCPU.h"

/**
 * \brief BVH traversal for rigid-deformable contact search
 * \ingroup contact
 */
class tledDeformableRigidBVHTraverserCPU : public tledBVHTraverserCPU {  
  /**
   * \name Detection
   * @{
   */
protected:
  static bool FindRootsFixedAndMovingEdge(float &r_t, float &r_r, float &r_q, const float fixedA[], const float fixedB[], const float movC0[], const float movD0[], const float movC1[], const float movD1[]);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
private:
  template <class TBV>
  static tledDeformableRigidBVHTraverserCPU* _CreateTraverserForBVType(tledSelfCollisionBVH &r_defBVH, const tledBVH &rigBVH);

public:
  static tledDeformableRigidBVHTraverserCPU* CreateTraverser(tledSelfCollisionBVH &r_defBVH, const tledBVH &rigBVH);

  virtual ~tledDeformableRigidBVHTraverserCPU(void) {}
  /** @} */
};

template <class TDeformableBVH, class TRigidBVH, class TAPI = tledDeformableRigidBVHTraverserCPU>
class tledDeformableRigidBVHTraverserImplCPU : public tledBVHTraverserImplCPU<TRigidBVH, TDeformableBVH, TAPI> {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledBVHTraverserImplCPU<TRigidBVH, TDeformableBVH, TAPI> Superclass;
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
  virtual void AddNarrowPhaseTests(const int sFacetInd, const int mFacetInd) { this->AddNarrowPhaseTestsSwitching(sFacetInd, mFacetInd); }
  
  virtual bool DoNarrowPhaseNodeFacetDefMaster(const int defPrimIndex, const int rigNodeIndex);
  virtual bool DoNarrowPhaseEdgeEdgeDefMaster(const int defEdgeIndex, const int rigEdgeIndex);
  virtual bool DoNarrowPhaseNodeFacetRigMaster(const int defNodeIndex, const int rigFacetIndex);
  virtual bool DoNarrowPhaseEdgeEdgeRigMaster(const int defEdgeIndex, const int rigEdgeIndex);

  virtual void ProcessNodeFacetNarrowPhaseItems(const std::vector<std::pair<int, int> >::const_iterator itemsBegin, const std::vector<std::pair<int, int> >::const_iterator itemsEnd);
  virtual void ProcessEdgeEdgeNarrowPhaseItems(const std::vector<std::pair<int, int> >::const_iterator itemsBegin, const std::vector<std::pair<int, int> >::const_iterator itemsEnd);

  virtual void RunBroadPhase(void);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledDeformableRigidBVHTraverserImplCPU(SlaveBVH &r_slaveBVH, const MasterBVH &masterBVH) : Superclass(r_slaveBVH, masterBVH) {}
  virtual ~tledDeformableRigidBVHTraverserImplCPU(void) {}
  /** @} */
};

#include "tledDeformableRigidBVHTraverserCPU.tpp"
#endif
