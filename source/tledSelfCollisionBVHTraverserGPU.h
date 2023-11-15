// =========================================================================
// File:       tledSelfCollisionBVHTraverserGPU.h
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
#ifndef tledSelfCollisionBVHTraverserGPU_H
#define tledSelfCollisionBVHTraverserGPU_H

#ifdef _GPU_

/**
 * \brief Self-collision detection 
 * \ingroup contact
 */
class tledSelfCollisionBVHTraverserGPU : public tledDeformableDeformableBVHTraverserGPU {
  /**
   * \name Construction, Destruction
   * @{
   */
public:
  static tledDeformableDeformableBVHTraverserGPU* CreateTraverser(tledSelfCollisionBVH &r_bvh);
  
  virtual ~tledDeformableDeformableBVHTraverserGPU(void) {}
  /** @} */  
};

template <class TBVH, class TAPI = tledSelfCollisionBVHTraverserGPU>
class tledSelfCollisionBVHTraverserImplGPU : public tledDeformableDeformableBVHTraverserImplGPU<TBVH, TAPI> {
    /**
   * \name Types
   * @{
   */
public:
  typedef tledDeformableDeformableBVHTraverserImplGPU<TBVH, TAPI> Superclass;
  typedef TBVH MasterBVH;
  typedef TBVH SlaveBVH;
  typedef typename MasterBVH::BoundingVolume BoundingVolume;
  typedef typename TBVH::ContactMesh MasterMesh;
  typedef typename TBVH::ContactMesh SlaveMesh;

  typedef tledDeformableDeformableBVHTraverserGPU::NodeFacetNarrowPhaseResult<MasterMesh::Facet::NumberOfVertices> NodeFacetNarrowPhaseResult;
  typedef tledDeformableDeformableBVHTraverserGPU::EdgeEdgeNarrowPhaseResult EdgeEdgeNarrowPhaseResult;
  /** @} */

  /**
   * \name Detection
   * @{
   */
protected:
  void InitStartBVPairs(void);
  /** @} */

  
};

#endif
#endif
