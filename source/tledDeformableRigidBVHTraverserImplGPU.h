// =========================================================================
// File:       tledDeformableRigidBVHTraverserImplGPU.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   CUDA
// Created:    January 2015
// 
// Copyright (c) 2015, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledDeformableRigidBVHTraverserImplGPU_H
#define tledDeformableRigidBVHTraverserImplGPU_H

#include "tledDeformableRigidBVHTraverserGPU.h"
#include "tledBVHTraverserImplGPU.h"
#include "tledCUDAMemoryBlock.h"

/**
 * \brief BVH traverser implementation for deformable-deformable collision detection.
 * \ingroup contact
 */
template <class TDeformableBVH, class TRigidBVH, class TAPI = tledDeformableRigidBVHTraverserGPU>
class tledDeformableRigidBVHTraverserImplGPU : public tledBVHTraverserImplGPU<TRigidBVH, TDeformableBVH, TAPI> {  
  /**
   * \name Types
   * @{
   */
public:
  typedef tledBVHTraverserImplGPU<TRigidBVH, TDeformableBVH, TAPI> Superclass;
  typedef TRigidBVH MasterBVH;
  typedef TDeformableBVH SlaveBVH;
  typedef typename MasterBVH::BoundingVolume BoundingVolume;
  typedef typename MasterBVH::ContactMesh MasterMesh;
  typedef typename SlaveBVH::ContactMesh SlaveMesh;
  typedef typename MasterMesh::GPUSurface MasterGPUSurface;
  typedef typename SlaveMesh::GPUSurface SlaveGPUSurface;
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
  const MasterGPUSurface* _GetMasterGPUSurface(void) const { return static_cast<const MasterGPUSurface*>(this->GetMasterMesh().GetDeviceGPUSurface()); }
  const SlaveGPUSurface* _GetSlaveGPUSurface(void) const { return static_cast<const SlaveGPUSurface*>(this->GetSlaveMesh().GetDeviceGPUSurface()); }

protected:
  virtual void InitStartBVPairs(void);
  
  virtual tledCUDADeviceMemoryBlock* InvokeNodeFacetNarrowPhaseStage1(int *dp_outItemCounter, cudaStream_t stream, const int2 *inItems, const int numNodeFacet);
  virtual tledCUDADeviceMemoryBlock* InvokeEdgeEdgeNarrowPhaseStage1(int *dp_outItemCounter, cudaStream_t stream, const int2 *inItems, const int numEdgeEdge);
  virtual tledCUDADeviceMemoryBlock* InvokeNodeFacetNarrowPhaseStage2(int *dp_outItemCounter, cudaStream_t stream, const NodeFacetInitialProjection *inItems, const int numNodeFacet);
  virtual tledCUDADeviceMemoryBlock* InvokeEdgeEdgeNarrowPhaseStage2(int *dp_outItemCounter, cudaStream_t stream, const EdgeEdgeInitialProjection *inItems, const int numEdgeEdge);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledDeformableRigidBVHTraverserImplGPU(SlaveBVH &r_slaveBVH, const MasterBVH &masterBVH) : Superclass(r_slaveBVH, masterBVH) {}
  virtual ~tledDeformableRigidBVHTraverserImplGPU(void) {}
  /** @} */
};

#ifdef __CUDACC__
#include "tledDeformableRigidBVHTraverserImplGPU.tpp"
#endif

#endif
