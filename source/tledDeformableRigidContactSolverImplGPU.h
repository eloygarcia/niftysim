// =========================================================================
// File:       tledDeformableRigidContactSolverImplGPU.h
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
#ifndef tledDeformableRigidContactSolverImplGPU_H
#define tledDeformableRigidContactSolverImplGPU_H

#include "tledDeformableRigidContactSolverGPU.h"
#include "tledDeformableRigidBVHTraverserGPU.h"
#include "tledRigidContactSurfaceGPU.h"

/**
 * \brief Deformable-rigid contact handler implementation.
 * \ingroup contact
 */
template <class TDeformableContactMesh, class TRigidContactMesh, class TAPI = tledDeformableRigidContactSolverGPU>
class tledDeformableRigidContactSolverImplGPU : public tledContactSolverImplGPU<TDeformableContactMesh, TAPI> {
  /** 
   * \name Types
   * @{
   */
public:
  typedef TDeformableContactMesh ContactMesh;
  typedef TRigidContactMesh RigidMesh;
  typedef typename ContactMesh::Facet Facet;
  typedef tledContactSolverImplGPU<TDeformableContactMesh, TAPI> Superclass;
  typedef tledDeformableRigidBVHTraverserGPU::NodeFacetNarrowPhaseResult<ContactMesh::Facet::NumberOfVertices> NodeProjection;
  typedef tledDeformableRigidBVHTraverserGPU::EdgeEdgeNarrowPhaseResult EdgeProjection;
  typedef tledContactSolverGPU::ContactResponse ContactResponse;
  typedef typename ContactMesh::GPUSurface GPUSurface;
  /** @} */

  /**
   * \name Responses
   * @{
   */
public:
  class NodeProjectionToResponseRigidMasterTransform;
  class NodeProjectionToResponseDeformableMasterTransform;

  template <const bool t_isDeformableMaster>
  class EdgeProjectionToResponseTransform;

protected:
  virtual void ComputeNodeFacetResponsesRigidMaster(ContactResponse *dp_dst, const float4 *dpc_uNexts, const float4 *dpc_uCurrs, const cudaStream_t stream);
  virtual void ComputeNodeFacetResponsesDeformableMaster(ContactResponse *dp_dst, const float4 *dpc_uNexts, const float4 *dpc_uCurrs, const cudaStream_t stream);

  virtual void ComputeEdgeEdgeResponsesRigidMaster(ContactResponse *dp_dst, const float4 *dpc_uNexts, const float4 *dpc_uCurrs, const cudaStream_t stream);
  virtual void ComputeEdgeEdgeResponsesDeformableMaster(ContactResponse *dp_dst, const float4 *dpc_uNexts, const float4 *dpc_uCurrs, const cudaStream_t stream);

  virtual void ComputeNodeFacetResponses(ContactResponse *dp_dst, const float4 *dpc_uNexts, const float4 *dpc_uCurrs, const cudaStream_t nodeFacetStream);
  virtual void ReleaseNodeFacetResponseComputationMemory(void);

  virtual void ComputeEdgeEdgeResponses(ContactResponse *dp_dst, const float4 *dpc_uNexts, const float4 *dpc_uCurrs, const cudaStream_t edgeEdgeStream);
  virtual void ReleaseEdgeEdgeResponseComputationMemory(void);

  virtual int GetNumberOfNodeFacetResponses(void) const { return (this->DoMaster()? 1 : ContactMesh::Facet::NumberOfVertices)*this->GetBVHTraverser().GetNumberOfNodeFacetContacts(); }
  virtual int GetNumberOfEdgeEdgeResponses(void) const { return 2*this->GetBVHTraverser().GetNumberOfEdgeEdgeContacts(); }

public:
  virtual bool DoFriction(void) const;
  /** @} */

  /**
   * \name Meshes
   * @{
   */
private:
  const RigidMesh *mpc_RigidSurface;

public:
  const RigidMesh& GetRigidSurface(void) const { return *mpc_RigidSurface; }
  ContactMesh& GetDeformableSurface(void) { return this->GetMesh(); }
  const ContactMesh& GetDeformableSurface(void) const { return this->GetMesh(); }
  /** @} */

  /**
   * \name Construction, Destruction, Initialisation
   * @{
   */
protected:
  virtual tledBVHTraverserGPU* InstantiateBVHTraverser(void);

public:
  tledDeformableRigidContactSolverImplGPU(tledUnstructuredContactManager &r_contactRes, const int surfaceIndex);
  virtual ~tledDeformableRigidContactSolverImplGPU(void) {}  
  /** @} */
};

#ifdef __CUDACC__
#include "tledDeformableRigidContactSolverImplGPU.tpp"
#endif

#endif
