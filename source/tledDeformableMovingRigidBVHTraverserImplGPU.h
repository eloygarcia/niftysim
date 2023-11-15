// =========================================================================
// File:       tledDeformableMovingRigidBVHTraverserImplGPU.h
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
#ifndef tledDeformableMovingRigidBVHTraverserImplGPU_H
#define tledDeformableMovingRigidBVHTraverserImplGPU_H

/** 
 * Deformable/dynamic rigid body contact search implementation. 
 * \ingroup contact
 */
template <class TDeformableBVH, class TRigidBVH, class TAPI = tledDeformableMovingRigidBVHTraverserGPU>
class tledDeformableMovingRigidBVHTraverserImplGPU : public tledDeformableRigidBVHTraverserImplGPU<TDeformableBVH, TRigidBVH, TAPI> {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledDeformableRigidBVHTraverserImplGPU<TDeformableBVH, TRigidBVH, TAPI> Superclass;
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
  /* 
   * MSVC++ 2010/CUDA 5.5 compiler bug (2014-10-09): typedef for NodeFacetInitialProjection, EdgeEdgeInitialProjection in class def. not recognised in cu-file. 
   * The following def's for {EdgeEdge,NodeFacet}InitialResult are compiler bug workarounds.
   */
  struct NodeFacetInitialProjection : public Superclass::NodeFacetInitialProjection {};
  typedef typename Superclass::NodeFacetNarrowPhaseResult NodeFacetNarrowPhaseResult;
  struct EdgeEdgeInitialProjection : public Superclass::EdgeEdgeInitialProjection {};
  typedef typename Superclass::EdgeEdgeNarrowPhaseResult EdgeEdgeNarrowPhaseResult;

private:
  const MasterGPUSurface* _GetMasterGPUSurface(void) const { return static_cast<const MasterGPUSurface*>(this->GetMasterMesh().GetDeviceGPUSurface()); }
  const SlaveGPUSurface* _GetSlaveGPUSurface(void) const { return static_cast<const SlaveGPUSurface*>(this->GetSlaveMesh().GetDeviceGPUSurface()); }

protected:
  virtual tledCUDADeviceMemoryBlock* InvokeNodeFacetNarrowPhaseStage1(int *dp_outItemCounter, cudaStream_t stream, const int2 *inItems, const int numNodeFacet);
  virtual tledCUDADeviceMemoryBlock* InvokeEdgeEdgeNarrowPhaseStage2(int *dp_outItemCounter, cudaStream_t stream, const typename Superclass::EdgeEdgeInitialProjection *inItems, const int numEdgeEdge);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledDeformableMovingRigidBVHTraverserImplGPU(SlaveBVH &r_slaveBVH, const MasterBVH &masterBVH) : Superclass(r_slaveBVH, masterBVH) {}
  virtual ~tledDeformableMovingRigidBVHTraverserImplGPU(void) {}
  /** @} */  
};

#ifdef __CUDACC__
#include "tledDeformableMovingRigidBVHTraverserImplGPU.tpp"
#endif

#endif
