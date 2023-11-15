// =========================================================================
// File:       tledDeformableRigidContactSolverImplGPU.tpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   CUDA
// Created:    July 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#include "tledCUDAMemoryBlock.h"
#include "tledCUDAHelpers.h"
#include "tledContactSolverGPU_kernels.h"

#include <cassert>

template <class TDeformableContactMesh, class TRigidContactMesh, class TAPI>
tledBVHTraverserGPU* tledDeformableRigidContactSolverImplGPU<TDeformableContactMesh, TRigidContactMesh, TAPI>::InstantiateBVHTraverser() {
  return tledDeformableRigidBVHTraverserGPU::CreateTraverser(static_cast<tledSelfCollisionBVHGPU&>(this->GetManager().GetDeformableBVH()), static_cast<const tledStaticBVHGPU&>(this->GetManager().GetRigidBVH(this->GetRigidSurfaceIndex())));
}

template <class TDeformableContactMesh, class TRigidContactMesh, class TAPI>
class tledDeformableRigidContactSolverImplGPU<TDeformableContactMesh, TRigidContactMesh, TAPI>::NodeProjectionToResponseRigidMasterTransform : public Superclass::NodeProjectionToResponseTransform {
protected:
  __device__ float3 ComputeRelativeNodeFacetVelocity(const NodeProjection &proj) const {
    return this->GetNodeVelocity(proj.ContactNodeIndices[0]);
  }

public:
  __device__ ContactResponse operator()(const NodeProjection &proj, const int nodeInd) const {
    float effDepth = this->ComputeEffectiveDepth(proj), lambda;
    ContactResponse response;

    tledCudaAssert(nodeInd == 0);
    response.NodeIndex = proj.ContactNodeIndices[0];
    if (effDepth < 0) {
      float beta;

      lambda = this->ComputeBaseLambdaFromDepth(effDepth);
      beta = this->template ComputeSlaveBeta<true, false>(proj);
      lambda *= this->GetNodeMass(proj.ContactNodeIndices[0])*beta;
      tledCudaPrintf(0, "Slave penetration (g = %f): lambda(%d) = %f\n", effDepth); 
    } else {
      float vMag = dot(this->ComputeRelativeNodeFacetVelocity(proj), proj.Normal);

      lambda = this->ComputeRateResponseWeight(effDepth)*this->template ComputeStoppingLambdaClamped<true, false>(vMag, proj)/this->GetTimeStep();
      tledCudaPrintf(0, "Slave rate-response (g = %f): lambda(%d) = %f\n", effDepth); 
    }

    response.Force = lambda*proj.Normal;
    tledCudaPrintf(0, "Applying node-facet contact force (slave) to %d: %f %f %f\n", response.NodeIndex, response.Force.x, response.Force.y, response.Force.z);

    return response;
  }

public:
  NodeProjectionToResponseRigidMasterTransform(const float *dpc_masses, const float3 *dpc_v, const float dt, const float off, const float rateDelta) : Superclass::NodeProjectionToResponseTransform(dpc_masses, dpc_v, dt, off, rateDelta) {}
};

template <class TDeformableContactMesh, class TRigidContactMesh, class TAPI>
class tledDeformableRigidContactSolverImplGPU<TDeformableContactMesh, TRigidContactMesh, TAPI>::NodeProjectionToResponseDeformableMasterTransform : public Superclass::NodeProjectionToResponseTransform {
protected:
  __device__ float3 ComputeRelativeNodeFacetVelocity(const NodeProjection &proj) const {
    float3 v = proj.ShapeValues[0]*this->GetNodeVelocity(proj.ContactNodeIndices[1]);

    for (int n = 1; n < ContactMesh::Facet::NumberOfVertices; n++) {
      v -= proj.ShapeValues[n]*this->GetNodeVelocity(proj.ContactNodeIndices[1+n]);
    }

    return v;
  }

public:
  __device__ ContactResponse operator()(const NodeProjection &proj, const int nodeInd) const {
    float effDepth = this->ComputeEffectiveDepth(proj), lambda;
    ContactResponse response;

    tledCudaAssert(nodeInd < ContactMesh::Facet::NumberOfVertices);
    response.NodeIndex = proj.ContactNodeIndices[1+nodeInd];
    if (effDepth < 0) {
      float beta;

      lambda = this->ComputeBaseLambdaFromDepth(effDepth);
      beta = 1 - this->template ComputeSlaveBeta<false, true>(proj);
      lambda *= -this->GetNodeMass(proj.ContactNodeIndices[1+nodeInd])*beta;
    } else {
      float vMag = dot(this->ComputeRelativeNodeFacetVelocity(proj), proj.Normal);

      /* No division by dt in v! */
      lambda = this->ComputeRateResponseWeight(effDepth)*this->template ComputeStoppingLambdaClamped<false, true>(vMag, proj)/this->GetTimeStep(); 
    }

    lambda *= proj.ShapeValues[nodeInd]/this->ComputeShapeSquareSum(proj);    
    response.Force = lambda*proj.Normal;
    tledCudaPrintf(0, "Applying node-facet contact force (master) to %d: %f %f %f\n", response.NodeIndex, response.Force.x, response.Force.y, response.Force.z);

    return response;
  }

public:
  NodeProjectionToResponseDeformableMasterTransform(const float *dpc_masses, const float3 *dpc_v, const float dt, const float off, const float rateDelta) : Superclass::NodeProjectionToResponseTransform(dpc_masses, dpc_v, dt, off, rateDelta) {}
};

template <class TDeformableContactMesh, class TRigidContactMesh, class TAPI>
template <const bool t_isDeformableMaster>
class tledDeformableRigidContactSolverImplGPU<TDeformableContactMesh, TRigidContactMesh, TAPI>::EdgeProjectionToResponseTransform : public Superclass::EdgeProjectionToResponseTransform {
protected:
  __device__ float ComputeRelativeEdgeEdgeVelocityMagnitude(const float4 &shapeVals, const EdgeProjection &proj) const {
    if (t_isDeformableMaster) return -dot(shapeVals.z*this->GetNodeVelocity(proj.MasterEdge.x) + shapeVals.w*this->GetNodeVelocity(proj.MasterEdge.y), proj.Normal)/this->GetTimeStep();
    else return dot(shapeVals.x*this->GetNodeVelocity(proj.SlaveEdge.x) + shapeVals.y*this->GetNodeVelocity(proj.SlaveEdge.y), proj.Normal)/this->GetTimeStep();
  }

  __device__ int GetNodeIndex(const EdgeProjection &proj, const int nodeInd) const {
    if (t_isDeformableMaster) return nodeInd == 0? proj.MasterEdge.x : proj.MasterEdge.y;
    else return nodeInd == 0? proj.SlaveEdge.x : proj.SlaveEdge.y;
  }

  __device__ const float& GetShapeFunctionValue(const float4 &shapeVals, const int nodeInd) const {
    if (t_isDeformableMaster) return nodeInd == 0? shapeVals.z : shapeVals.w;
    else return nodeInd == 0? shapeVals.x : shapeVals.y;
  }

public:
  __device__ ContactResponse operator()(const EdgeProjection &proj, const int nodeInd) const {
    const bool doSlave = t_isDeformableMaster? false : true;
    const bool doMaster = t_isDeformableMaster? true : false;

    float effDepth = this->ComputeEffectiveDepth(proj), lambda;
    float4 shapeVals = this->template ComputeShapeValues<doSlave, doMaster>(proj); 
    ContactResponse response;

    tledCudaAssert(nodeInd >= 0 && nodeInd < 2);
    response.NodeIndex = this->GetNodeIndex(proj, nodeInd);
    if (effDepth > 0) {
      float vMag = this->ComputeRelativeEdgeEdgeVelocityMagnitude(shapeVals, proj);

      lambda = -this->ComputeRateResponseWeight(effDepth)*this->template ComputeStoppingLambdaClamped<doSlave, doMaster>(vMag, shapeVals, proj);
    } else lambda = this->GetNodeMass(response.NodeIndex)*this->ComputeBaseLambdaFromDepth(effDepth);
    lambda = this->GetShapeFunctionValue(shapeVals, nodeInd)*lambda/this->template ComputeShapeSquareSum<doSlave, doMaster>(shapeVals);
    if (t_isDeformableMaster) lambda *= -1.f;
    response.Force = lambda*proj.Normal;
    tledCudaPrintf(0, "Applying edge-edge contact force (%s) to %d: %f %f %f\n", (t_isDeformableMaster? "master" : "slave"), response.NodeIndex, response.Force.x, response.Force.y, response.Force.z);
    
    return response;
  }

public:
  EdgeProjectionToResponseTransform(const float *dpc_masses, const float3 *dpc_v, const float dt, const float off, const float rateThrsh) : Superclass::EdgeProjectionToResponseTransform(dpc_masses, dpc_v, dt, off, rateThrsh) {}
};

template <class TDeformableContactMesh, class TRigidContactMesh, class TAPI>
tledDeformableRigidContactSolverImplGPU<TDeformableContactMesh, TRigidContactMesh, TAPI>::tledDeformableRigidContactSolverImplGPU(tledUnstructuredContactManager &r_contactRes, const int surfaceIndex) : Superclass(r_contactRes) { 
  this->SetRigidSurfaceIndex(surfaceIndex); 
  mpc_RigidSurface = &this->GetManager().template GetRigidSurface<TRigidContactMesh>(this->GetRigidSurfaceIndex());
}

template <class TDeformableContactMesh, class TRigidContactMesh, class TAPI>
void tledDeformableRigidContactSolverImplGPU<TDeformableContactMesh, TRigidContactMesh, TAPI>::ComputeNodeFacetResponsesRigidMaster(ContactResponse *dp_dst, const float4 *dpc_uNexts, const float4 *dpc_uCurrs, const cudaStream_t stream) {
  const tledCUDADeviceMemoryBlock &projections = this->GetBVHTraverser().GetNodeFacetResults();
  const float *dpc_nodeMasses = this->GetHostGPUSurface().NodeMasses;
  const int blockSize = 128;
  const int numContacts = this->GetBVHTraverser().GetNumberOfNodeFacetContacts();
  const int numBlocks = tledCUDAHelpers::GetNumberOfBlocks(numContacts, blockSize);
  const NodeProjectionToResponseRigidMasterTransform xForm(dpc_nodeMasses, this->GetNodeVelocities(), this->GetDt(), this->GetManager().GetSafetyMargin(), this->GetNodeCloseDistance());

  tledContactSolverGPU_kernels::ProjectionToResponseTransformKernel<NodeProjection, NodeProjectionToResponseRigidMasterTransform> <<<numBlocks, blockSize, 0, stream>>> (dp_dst, projections.GetBuffer<NodeProjection>(), numContacts, xForm, 1);
  tledDeviceSyncDebug;
}

template <class TDeformableContactMesh, class TRigidContactMesh, class TAPI>
void tledDeformableRigidContactSolverImplGPU<TDeformableContactMesh, TRigidContactMesh, TAPI>::ComputeNodeFacetResponsesDeformableMaster(ContactResponse *dp_dst, const float4 *dpc_uNexts, const float4 *dpc_uCurrs, const cudaStream_t stream) {
  const tledCUDADeviceMemoryBlock &projections = this->GetBVHTraverser().GetNodeFacetResults();
  const float *dpc_nodeMasses = this->GetHostGPUSurface().NodeMasses;
  const int numResponsesPerContact = TDeformableContactMesh::Facet::NumberOfVertices;
  const int blockSize = 128;
  const int numContacts = this->GetBVHTraverser().GetNumberOfNodeFacetContacts();
  const int numBlocks = tledCUDAHelpers::GetNumberOfBlocks(numContacts*numResponsesPerContact, blockSize);
  const NodeProjectionToResponseDeformableMasterTransform xForm(dpc_nodeMasses, this->GetNodeVelocities(), this->GetDt(), this->GetManager().GetSafetyMargin(), this->GetNodeCloseDistance());

  tledContactSolverGPU_kernels::ProjectionToResponseTransformKernel<NodeProjection, NodeProjectionToResponseDeformableMasterTransform> <<<blockSize, numBlocks, 0, stream>>> (dp_dst, projections.GetBuffer<NodeProjection>(), numContacts, xForm, numResponsesPerContact);
  tledDeviceSyncDebug;
}

template <class TDeformableContactMesh, class TRigidContactMesh, class TAPI>
void tledDeformableRigidContactSolverImplGPU<TDeformableContactMesh, TRigidContactMesh, TAPI>::ReleaseNodeFacetResponseComputationMemory() {
}

template <class TDeformableContactMesh, class TRigidContactMesh, class TAPI>
void tledDeformableRigidContactSolverImplGPU<TDeformableContactMesh, TRigidContactMesh, TAPI>::ReleaseEdgeEdgeResponseComputationMemory() {
}

template <class TDeformableContactMesh, class TRigidContactMesh, class TAPI>
void tledDeformableRigidContactSolverImplGPU<TDeformableContactMesh, TRigidContactMesh, TAPI>::ComputeEdgeEdgeResponsesRigidMaster(ContactResponse *dp_dst, const float4 *dpc_uNexts, const float4 *dpc_uCurrs, const cudaStream_t stream) {
  typedef EdgeProjectionToResponseTransform<false> __XForm;

  const tledCUDADeviceMemoryBlock &projections = this->GetBVHTraverser().GetEdgeEdgeResults();
  const float *dpc_nodeMasses = this->GetHostGPUSurface().NodeMasses;
  const int blockSize = 128;
  const int numBlocks = tledCUDAHelpers::GetNumberOfBlocks(this->GetBVHTraverser().GetNumberOfEdgeEdgeContacts()*2, blockSize);	
  const __XForm xForm(dpc_nodeMasses, this->GetNodeVelocities(), this->GetDt(), this->GetManager().GetSafetyMargin(), this->GetNodeCloseDistance());

  assert(projections.IsActive());
  tledContactSolverGPU_kernels::ProjectionToResponseTransformKernel<EdgeProjection, __XForm> <<<numBlocks, blockSize, 0, stream>>> (dp_dst, projections.GetBuffer<EdgeProjection>(), this->GetBVHTraverser().GetNumberOfEdgeEdgeContacts(), xForm, 2);
}

template <class TDeformableContactMesh, class TRigidContactMesh, class TAPI>
void tledDeformableRigidContactSolverImplGPU<TDeformableContactMesh, TRigidContactMesh, TAPI>::ComputeEdgeEdgeResponsesDeformableMaster(ContactResponse *dp_dst, const float4 *dpc_uNexts, const float4 *dpc_uCurrs, const cudaStream_t stream) {
  typedef EdgeProjectionToResponseTransform<true> __XForm;

  const tledCUDADeviceMemoryBlock &projections = this->GetBVHTraverser().GetEdgeEdgeResults();
  const int blockSize = 128;
  const int numBlocks = tledCUDAHelpers::GetNumberOfBlocks(this->GetBVHTraverser().GetNumberOfEdgeEdgeContacts()*2, blockSize);	
  const float *dpc_nodeMasses = this->GetHostGPUSurface().NodeMasses;
  const __XForm xForm(dpc_nodeMasses, this->GetNodeVelocities(), this->GetDt(), this->GetManager().GetSafetyMargin(), this->GetNodeCloseDistance());

  assert(projections.IsActive());
  tledContactSolverGPU_kernels::ProjectionToResponseTransformKernel<EdgeProjection, __XForm> <<<numBlocks, blockSize, 0, stream>>> (dp_dst, projections.GetBuffer<EdgeProjection>(), this->GetBVHTraverser().GetNumberOfEdgeEdgeContacts(), xForm, 2);
}

template <class TDeformableContactMesh, class TRigidContactMesh, class TAPI>
bool tledDeformableRigidContactSolverImplGPU<TDeformableContactMesh, TRigidContactMesh, TAPI>::DoFriction() const { 
  return this->GetRigidSurface().GetFrictionCoefficient() > 0; 
}

template <class TDeformableContactMesh, class TRigidContactMesh, class TAPI>
void tledDeformableRigidContactSolverImplGPU<TDeformableContactMesh, TRigidContactMesh, TAPI>::ComputeNodeFacetResponses(ContactResponse *dp_dst, const float4 *dpc_uNexts, const float4 *dpc_uCurrs, const cudaStream_t nodeFacetStream) {
  if (this->DoMaster()) this->ComputeNodeFacetResponsesRigidMaster(dp_dst, dpc_uNexts, dpc_uCurrs, nodeFacetStream);
  else this->ComputeNodeFacetResponsesDeformableMaster(dp_dst, dpc_uNexts, dpc_uCurrs, nodeFacetStream);
}

template <class TDeformableContactMesh, class TRigidContactMesh, class TAPI>
void tledDeformableRigidContactSolverImplGPU<TDeformableContactMesh, TRigidContactMesh, TAPI>::ComputeEdgeEdgeResponses(ContactResponse *dp_dst, const float4 *dpc_uNexts, const float4 *dpc_uCurrs, const cudaStream_t edgeEdgeStream) {
  if (this->DoMaster()) this->ComputeEdgeEdgeResponsesRigidMaster(dp_dst, dpc_uNexts, dpc_uCurrs, edgeEdgeStream);
  else this->ComputeEdgeEdgeResponsesDeformableMaster(dp_dst, dpc_uNexts, dpc_uCurrs, edgeEdgeStream);
}
