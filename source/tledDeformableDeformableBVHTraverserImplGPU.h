// =========================================================================
// File:       tledDeformableDeformableBVHTraverserImplGPU.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    January 2015
// 
// Copyright (c) 2015, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledDeformableDeformableBVHTraverserImplGPU_H
#define tledDeformableDeformableBVHTraverserImplGPU_H

#include "tledDeformableDeformableBVHTraverserGPU.h"
#include "tledBVHTraverserImplGPU.h"

/**
 * \brief BVH traverser implementation for deformable-deformable collision detection.
 * \ingroup contact
 */
template <class TBVH, class TAPI = tledDeformableDeformableBVHTraverserGPU>
class tledDeformableDeformableBVHTraverserImplGPU : public tledBVHTraverserImplGPU<TBVH, TBVH, TAPI> {  
  /**
   * \name Types
   * @{
   */
public:
  typedef tledBVHTraverserImplGPU<TBVH, TBVH, TAPI> Superclass;
  typedef TBVH MasterBVH;
  typedef TBVH SlaveBVH;
  typedef typename MasterBVH::BoundingVolume BoundingVolume;
  typedef typename TBVH::ContactMesh MasterMesh;
  typedef typename TBVH::ContactMesh SlaveMesh;
  typedef typename MasterMesh::GPUSurface GPUSurface;
  /** @} */

  /**
   * \name Detection
   * @{
   */
public:
  typedef typename Superclass::NodeFacetInitialProjection NodeFacetInitialProjection;
  typedef typename Superclass::NodeFacetNarrowPhaseResult NodeFacetNarrowPhaseResult;
  typedef typename Superclass::EdgeEdgeInitialProjection EdgeEdgeInitialProjection;
  typedef typename Superclass::EdgeEdgeNarrowPhaseResult EdgeEdgeNarrowPhaseResult;

private:
  const GPUSurface* _GetGPUSurface(void) const { return static_cast<const GPUSurface*>(this->GetMasterMesh().GetDeviceGPUSurface()); }

protected:
  virtual int GetCurrentMasterBVHOrder(void) const { return this->GetMasterBVH().GetBVHOrder(); }
  virtual void InitStartBVPairs(void);

  virtual const std::vector<int>& GetHostStartBVs(void) const;

  virtual tledCUDADeviceMemoryBlock* InvokeNodeFacetNarrowPhaseStage1(int *dp_outItemCounter, cudaStream_t stream, const int2 *inItems, const int numNodeFacet);
  virtual tledCUDADeviceMemoryBlock* InvokeEdgeEdgeNarrowPhaseStage1(int *dp_outItemCounter, cudaStream_t stream, const int2 *inItems, const int numEdgeEdge);
  virtual tledCUDADeviceMemoryBlock* InvokeNodeFacetNarrowPhaseStage2(int *dp_outItemCounter, cudaStream_t stream, const NodeFacetInitialProjection *inItems, const int numNodeFacet);
  virtual tledCUDADeviceMemoryBlock* InvokeEdgeEdgeNarrowPhaseStage2(int *dp_outItemCounter, cudaStream_t stream, const EdgeEdgeInitialProjection *inItems, const int numEdgeEdge);
  /** @} */

  /**
   * \name Construction, Destruction, Initialisation
   * @{
   */
protected:
  /** For compatibility with template wrappers */
  tledDeformableDeformableBVHTraverserImplGPU(TBVH &r_slaveBVH, const TBVH &masterBVH) : Superclass(r_slaveBVH, masterBVH) {}

public:
  tledDeformableDeformableBVHTraverserImplGPU(TBVH &r_bvh) : Superclass(r_bvh, r_bvh) {}

  virtual ~tledDeformableDeformableBVHTraverserImplGPU(void);
  /** @} */
};

#if defined __CUDACC__ && !defined  __GPU_TEST_LINK_CONFLICT_NO_INCLUDE
#include "tledDeformableDeformableBVHTraverserImplGPU.tpp"
#endif

#endif
