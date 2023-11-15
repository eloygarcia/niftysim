// =========================================================================
// File:       tledDeformableDeformableContactSolverImplGPU.h
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
#ifndef tledDeformableDeformableContactSolverImplGPU_H
#define tledDeformableDeformableContactSolverImplGPU_H

#include "tledDeformableDeformableContactSolverGPU.h"
#include "tledContactSolverImplGPU.h"
#include "tledDeformableDeformableBVHTraverserGPU.h"

/**
 * \brief Deformable-deformable contact handler implementation.
 * \ingroup contact
 */
template <class TContactMesh, class TAPI = tledDeformableDeformableContactSolverGPU>
class tledDeformableDeformableContactSolverImplGPU : public tledContactSolverImplGPU<TContactMesh, TAPI> {
  /** 
   * \name Types
   * @{
   */
public:
  typedef TContactMesh ContactMesh;
  typedef typename ContactMesh::Facet Facet;
  typedef tledContactSolverImplGPU<TContactMesh, TAPI> Superclass;
  typedef tledContactSolverGPU::ContactResponse ContactResponse;
  typedef typename ContactMesh::GPUSurface GPUSurface;
  typedef typename Superclass::NodeProjection NodeProjection;
  typedef typename Superclass::EdgeProjection EdgeProjection;
  /** @} */

  /**
   * \name Responses
   * @{
   */
public:
  class NodeProjectionToResponseTransform;
  class EdgeProjectionToResponseTransform;

protected:
  virtual void ComputeNodeFacetResponses(ContactResponse *dp_dst, const float4 *dpc_uNexts, const float4 *dpc_uCurrs, const cudaStream_t nodeFacetStream);
  virtual void ReleaseNodeFacetResponseComputationMemory(void);  

  virtual void ComputeEdgeEdgeResponses(ContactResponse *dp_dst, const float4 *dpc_uNexts, const float4 *dpc_uCurrs, const cudaStream_t edgeEdgeStream);
  virtual void ReleaseEdgeEdgeResponseComputationMemory(void);

  virtual int GetNumberOfNodeFacetResponses(void) const { return (1 + ContactMesh::Facet::NumberOfVertices)*this->GetBVHTraverser().GetNumberOfNodeFacetContacts(); }
  virtual int GetNumberOfEdgeEdgeResponses(void) const { return 4*this->GetBVHTraverser().GetNumberOfEdgeEdgeContacts(); }

public:
  virtual bool DoFriction(void) const { return this->GetMesh().GetFrictionCoefficient() > 0; }
  /** @} */

  /**
   * \name Construction, Destruction, Initialisation
   * @{
   */
protected:
  virtual tledBVHTraverserGPU* InstantiateBVHTraverser(void);

public:
  tledDeformableDeformableContactSolverImplGPU(tledUnstructuredContactManager &r_contactRes) : Superclass(r_contactRes) {}
  virtual ~tledDeformableDeformableContactSolverImplGPU(void) {}  
  /** @} */
};

#ifdef __CUDACC__
#include "tledDeformableDeformableContactSolverImplGPU.tpp"
#endif

#endif
