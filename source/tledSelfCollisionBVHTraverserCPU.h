// =========================================================================
// File:       tledSelfCollisionBVHTraverserCPU.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    July 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledSelfCollisionBVHTraverserCPU_H
#define tledSelfCollisionBVHTraverserCPU_H

#include "tledDeformableDeformableBVHTraverserCPU.h"

/**
 * \brief BVH traverser for detecting deformable-deformable contacts including self-collisions.
 * \ingroup contact
 */
class tledSelfCollisionBVHTraverserCPU : public tledDeformableDeformableBVHTraverserCPU {
  /**
   * \name Construction, Destruction
   * @{
   */
public:
  static tledSelfCollisionBVHTraverserCPU* CreateTraverser(tledSelfCollisionBVH &r_bvh);

  virtual ~tledSelfCollisionBVHTraverserCPU(void) {}
  /** @} */
};

/**
 * \ingroup contact
 */
template <class TBVH, class TAPI = tledSelfCollisionBVHTraverserCPU>
class tledSelfCollisionBVHTraverserImplCPU : public tledDeformableDeformableBVHTraverserImplCPU<TBVH, TAPI> {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledDeformableDeformableBVHTraverserImplCPU<TBVH, TAPI> Superclass;
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
  virtual void AddNarrowPhaseTests(const int sFacetInd, const int mFacetInd); 

  /** Returns a list of subtree-root BV indices that need to be checked for collisions. */
  virtual const std::vector<int>& GetCollisionCandidateBVs(void) const { return this->GetMasterBVH().GetSelfCollisionCandidates(); }
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
protected:
  tledSelfCollisionBVHTraverserImplCPU(TBVH &r_slaveBVH, const TBVH &masterBVH) : Superclass(r_slaveBVH, masterBVH) {}

public:
  tledSelfCollisionBVHTraverserImplCPU(TBVH &r_bvh) : Superclass(r_bvh) {}
  virtual ~tledSelfCollisionBVHTraverserImplCPU(void) {}
  /** @} */
};

#include "tledSelfCollisionBVHTraverserCPU.tpp"
#endif
