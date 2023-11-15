// =========================================================================
// File:       tledDeformableMovingRigidContactSolverImplGPU.h
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
#ifndef tledDeformableMovingRigidContactSolverImplGPU_H
#define tledDeformableMovingRigidContactSolverImplGPU_H

#include "tledDeformableMovingRigidContactSolverGPU.h"
#include "tledDeformableRigidContactSolverImplGPU.h"

template <class TDeformableContactMesh, class TRigidContactMesh, class TAPI = tledDeformableMovingRigidContactSolverGPU>
class tledDeformableMovingRigidContactSolverImplGPU : public tledDeformableRigidContactSolverImplGPU<TDeformableContactMesh, TRigidContactMesh, TAPI> {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledDeformableRigidContactSolverImplGPU<TDeformableContactMesh, TRigidContactMesh, TAPI> Superclass;
  typedef TDeformableContactMesh ContactMesh;
  typedef TRigidContactMesh RigidMesh;
  typedef typename ContactMesh::Facet Facet;
  typedef typename ContactMesh::GPUSurface GPUSurface;
  typedef tledDeformableRigidBVHTraverserGPU::NodeFacetNarrowPhaseResult<ContactMesh::Facet::NumberOfVertices> NodeProjection;
  typedef tledDeformableRigidBVHTraverserGPU::EdgeEdgeNarrowPhaseResult EdgeProjection;
  /** @} */

  /**
   * \name Contact Search
   * @{
   */
protected:
  virtual tledBVHTraverserGPU* InstantiateBVHTraverser(void);
  /** @} */

  /**
   * \name Responses
   * @{
   */
public:
  typedef typename Superclass::ContactResponse ContactResponse;

  class NodeProjectionToResponseRigidMasterTransform;
  class NodeProjectionToResponseDeformableMasterTransform;

  template <const bool t_isDeformableMaster>
  class EdgeProjectionToResponseTransform;

private:
  class _MovingRigidTransformBase;

private:
  static float3 _MakeTimeStepIncrement(const float xform[], const int numSteps);

protected:
  virtual void ComputeNodeFacetResponsesRigidMaster(ContactResponse *dp_dst, const float4 *dpc_uNexts, const float4 *dpc_uCurrs, const cudaStream_t stream);
  virtual void ComputeNodeFacetResponsesDeformableMaster(ContactResponse *dp_dst, const float4 *dpc_uNexts, const float4 *dpc_uCurrs, const cudaStream_t stream);

  virtual void ComputeEdgeEdgeResponsesRigidMaster(ContactResponse *dp_dst, const float4 *dpc_uNexts, const float4 *dpc_uCurrs, const cudaStream_t stream);
  virtual void ComputeEdgeEdgeResponsesDeformableMaster(ContactResponse *dp_dst, const float4 *dpc_uNexts, const float4 *dpc_uCurrs, const cudaStream_t stream);
  /** @} */

  /**
   * \name Construction and Destruction
   * @{
   */
public:
  tledDeformableMovingRigidContactSolverImplGPU(tledUnstructuredContactManager &r_contactManager, const int rigidSurfaceIndex) : Superclass(r_contactManager, rigidSurfaceIndex) {}
  virtual ~tledDeformableMovingRigidContactSolverImplGPU(void) {}
  /** @} */
};

#ifdef __CUDACC__
#include "tledDeformableMovingRigidContactSolverImplGPU.tpp"
#endif

#endif
