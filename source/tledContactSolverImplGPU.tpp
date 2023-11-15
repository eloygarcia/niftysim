// =========================================================================
// File:       tledContactSolverImplGPU.tpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   CUDA
// Created:    May 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledContactSolverImplGPU_TPP
#define tledContactSolverImplGPU_TPP

#include "tledCUDAMemoryBlock.h"

#include <thrust/sort.h>
#include <thrust/unique.h>
#include <thrust/reduce.h>
#include <thrust/iterator/discard_iterator.h>
#include <thrust/iterator/counting_iterator.h>

#include "tledContactSolverGPU_kernels.tpp"

template <class TContactMesh, class TAPI>
__device__ float tledContactSolverImplGPU<TContactMesh, TAPI>::ProjectionToResponseTransform::ComputeRateResponseWeight(const float effectiveDepth) const {
  const float w = fmaxf(0.f, 1 - effectiveDepth/this->GetRateResponseMaxDistance());
  
  tledCudaAssert(w >= -1e-6f && w < 1 + 1e-6f);
  tledCudaPrintf(2, "Rate response weight g = %f (tau = %f): w = %f\n", effectiveDepth, this->GetRateResponseMaxDistance(), w);

  return w;
}

template <class TContactMesh, class TAPI>
__device__ float tledContactSolverImplGPU<TContactMesh, TAPI>::ProjectionToResponseTransform::ComputeBaseLambdaFromDepth(const float effectiveDepth) const {
  return effectiveDepth/(this->GetTimeStep()*this->GetTimeStep());
}

template <class TContactMesh, class TAPI>
template <const bool t_doSlave, const bool t_doMaster>
__device__ float tledContactSolverImplGPU<TContactMesh, TAPI>::NodeProjectionToResponseTransform::ComputeStoppingLambda(const float vMag, const NodeProjection &proj) const {
  float lambda = 0;

  if (t_doMaster) {
    for (int n = 0; n < ContactMesh::Facet::NumberOfVertices; n++) {
      lambda += proj.ShapeValues[n]*proj.ShapeValues[n]/this->GetNodeMass(proj.ContactNodeIndices[1+n]);      
    }
    lambda /= this->ComputeShapeSquareSum(proj);
  }

  if (t_doSlave) {
    lambda += 1/this->GetNodeMass(proj.ContactNodeIndices[0]);
  }
    
  return fabsf(vMag/(this->GetTimeStep()*lambda));
}

template <class TContactMesh, class TAPI>
template <const bool t_doSlave, const bool t_doMaster>
__device__ float tledContactSolverImplGPU<TContactMesh, TAPI>::NodeProjectionToResponseTransform::ComputeStoppingLambdaClamped(const float vMag, const NodeProjection &proj) const {
  return vMag > 0? 0.f : this->template ComputeStoppingLambda<t_doSlave, t_doMaster>(vMag, proj); 
}

template <class TContactMesh, class TAPI>
__device__ float tledContactSolverImplGPU<TContactMesh, TAPI>::NodeProjectionToResponseTransform::ComputeShapeSquareSum(const NodeProjection &proj) const {
  float shapeSqrSum = 0;

  for (int n = 0; n < ContactMesh::Facet::NumberOfVertices; n++) shapeSqrSum += proj.ShapeValues[n]*proj.ShapeValues[n];
    
  return shapeSqrSum;
}

template <class TContactMesh, class TAPI>
template <const bool t_doSlave, const bool t_doMaster>
__device__ float tledContactSolverImplGPU<TContactMesh, TAPI>::NodeProjectionToResponseTransform::ComputeSlaveBeta(const NodeProjection &proj) const {
  if (t_doMaster != t_doSlave) return t_doSlave? 1 : 0;
  else {
    float beta = 0;
    
    for (int n = 0; n < ContactMesh::Facet::NumberOfVertices; n++) {
      beta += proj.ShapeValues[n]*this->GetNodeMass(proj.ContactNodeIndices[1+n]);
    }
    beta = beta/(beta + this->GetNodeMass(proj.ContactNodeIndices[0]));
      
    return beta;
  }
}

template <class TContactMesh, class TAPI>
template <const bool t_doSlave, const bool t_doMaster>
__device__ float tledContactSolverImplGPU<TContactMesh, TAPI>::EdgeProjectionToResponseTransform::ComputeShapeSquareSum(const float4 &shapeVals) const {
  float shapeSqr = 0.f;

  if (t_doSlave) shapeSqr += shapeVals.x*shapeVals.x + shapeVals.y*shapeVals.y;
  if (t_doMaster) shapeSqr += shapeVals.z*shapeVals.z + shapeVals.w*shapeVals.w;

  return shapeSqr;
}

template <class TContactMesh, class TAPI>
template <const bool t_doSlave, const bool t_doMaster>
__device__ float tledContactSolverImplGPU<TContactMesh, TAPI>::EdgeProjectionToResponseTransform::ComputeStoppingLambda(const float vMag, const float4 &shapeVals, const EdgeProjection &proj) const {
  float lambda = 0.f;

  if (t_doMaster) {
    float masterDenom = 0.f;
    float sqrSum = 0.f;

    masterDenom += shapeVals.z*shapeVals.z/this->GetNodeMass(proj.MasterEdge.x);
    sqrSum += shapeVals.z*shapeVals.z;

    masterDenom += shapeVals.w*shapeVals.w/this->GetNodeMass(proj.MasterEdge.y);
    sqrSum += shapeVals.w*shapeVals.w;
    
    lambda += masterDenom/sqrSum;
  } 

  if (t_doSlave) {    
    float slaveDenom = 0.f;
    float sqrSum = 0.f;
    
    slaveDenom += shapeVals.x*shapeVals.x/this->GetNodeMass(proj.SlaveEdge.x);
    sqrSum += shapeVals.x*shapeVals.x;

    slaveDenom += shapeVals.y*shapeVals.y/this->GetNodeMass(proj.SlaveEdge.y);
    sqrSum += shapeVals.y*shapeVals.y;

    lambda += slaveDenom/sqrSum;
  }

  lambda = std::fabs(vMag/(this->GetTimeStep()*lambda));

  return lambda;
}

template <class TContactMesh, class TAPI>
template <const bool t_doSlave, const bool t_doMaster>
__device__ float tledContactSolverImplGPU<TContactMesh, TAPI>::EdgeProjectionToResponseTransform::ComputeStoppingLambdaClamped(const float vMag, const float4 &shapeVals, const EdgeProjection &proj) const {
  return vMag > 0? 0.f : this->template ComputeStoppingLambda<t_doSlave, t_doMaster>(vMag, shapeVals, proj);
}

template <class TContactMesh, class TAPI>
template <const bool t_doSlave, const bool t_doMaster>
__device__ float4 tledContactSolverImplGPU<TContactMesh, TAPI>::EdgeProjectionToResponseTransform::ComputeShapeValues(const EdgeProjection &proj) const {
  float4 shapeVals;

  if (t_doSlave) {
    shapeVals.x = (1 - proj.Xi.y);
    shapeVals.y = proj.Xi.y;
  } 

  if (t_doMaster) {
    shapeVals.z = (1 - proj.Xi.z);
    shapeVals.w = proj.Xi.z;
  }

  return shapeVals;
}

template <class TContactMesh, class TAPI>
void tledContactSolverImplGPU<TContactMesh, TAPI>::InitBVHTraverser() {
  mp_BVHTraverser->SetNarrowPhaseMaxDistance(this->GetNodeCloseDistance());
  mp_BVHTraverser->Init(this->GetManager());
}

template <class TContactMesh, class TAPI>
void tledContactSolverImplGPU<TContactMesh, TAPI>::ConsolidateResponses(float4 *dp_R, tledCUDADeviceMemoryBlock &r_responses, const int numResponses) const {
  const int blockSize = 256;

  tledLogDebugStream(tledHelper::Info() << "Have " << numResponses << " response forces.");
  thrust::sort(r_responses.GetThrustBuffer<ContactResponse>(), r_responses.GetThrustBuffer<ContactResponse>() + numResponses, tledContactSolverGPU_kernels::ContactResponseOrdering());
  tledDeviceSyncDebug;
  tledContactSolverGPU_kernels::ConvertResponsesToExternalForce <<<tledCUDAHelpers::GetNumberOfBlocks(numResponses, blockSize), blockSize>>> (dp_R, r_responses.GetBuffer<ContactResponse>(), numResponses, this->GetHostGPUSurface().SurfaceToVolumeNodeIndexMap);
  tledDeviceSyncDebug;
  r_responses.ToggleActive();  
  assert(!r_responses.IsActive());
}

template <class TContactMesh, class TAPI>
tledContactSolverImplGPU<TContactMesh, TAPI>::~tledContactSolverImplGPU() {
  if (mp_BVHTraverser != NULL) delete mp_BVHTraverser;
}

template <class TContactMesh, class TAPI>
void tledContactSolverImplGPU<TContactMesh, TAPI>::Init() {
  mp_BVHTraverser = this->InstantiateBVHTraverser();
  this->InitBVHTraverser();
}

template <class TContactMesh, class TAPI>
bool tledContactSolverImplGPU<TContactMesh, TAPI>::FindCollisions() {
  this->GetBVHTraverser().SetDoMaster(this->DoMaster());  
  this->GetBVHTraverser().FindCollisions();
  assert(this->GetBVHTraverser().GetNumberOfNodeFacetContacts() == 0 || this->GetBVHTraverser().GetNodeFacetResults().IsActive());
  assert(this->GetBVHTraverser().GetNumberOfEdgeEdgeContacts() == 0 || this->GetBVHTraverser().GetEdgeEdgeResults().IsActive());

  return this->GetBVHTraverser().GetNumberOfNodeFacetContacts() > 0 || this->GetBVHTraverser().GetNumberOfEdgeEdgeContacts() > 0;
}

template <class TContactMesh, class TAPI>
bool tledContactSolverImplGPU<TContactMesh, TAPI>::ComputeContactResponses(float4 *dp_R, const float4 *dpc_uNexts, const float4 *dpc_uCurrs) {
  tledLogDebugStream(tledHelper::Info() << "Computing contact responses.");
  this->ToggleDoMaster();  

  if (this->FindCollisions()) {    
    const int numNodeFacetContacts = this->GetBVHTraverser().GetNumberOfNodeFacetContacts();
    const int numEdgeEdgeContacts = this->GetBVHTraverser().GetNumberOfEdgeEdgeContacts();       

    cudaStream_t nfStream, eeStream;
    tledCUDADeviceMemoryBlock &r_responses = tledCUDADeviceMemoryBlock::GetNextFreeBufferWithSize<ContactResponse>(this->GetNumberOfNodeFacetResponses() + this->GetNumberOfEdgeEdgeResponses());

    tledCUDADeviceMemoryBlock::SaveAllocationCounter(-(numEdgeEdgeContacts > 0) - (numNodeFacetContacts > 0));

    this->ComputeNodeVelocities(dpc_uNexts, dpc_uCurrs);
    assert(numNodeFacetContacts > 0 || this->GetNumberOfNodeFacetResponses() == 0);
    assert(numEdgeEdgeContacts > 0 || this->GetNumberOfEdgeEdgeResponses() == 0);

    if (numNodeFacetContacts > 0) {
      tledCheckCUDAErrors(cudaStreamCreate(&nfStream));
      tledLogDebugStream(tledHelper::Info() << "Processing " << numNodeFacetContacts << " node-facet contacts, expecting " << this->GetNumberOfNodeFacetResponses() << " responses.");
      this->ComputeNodeFacetResponses(r_responses.GetBuffer<ContactResponse>(), dpc_uNexts, dpc_uCurrs, nfStream);
    }

    if (numEdgeEdgeContacts > 0) {
      tledCheckCUDAErrors(cudaStreamCreate(&eeStream));
      tledLogDebugStream(tledHelper::Info() << "Processing " << numEdgeEdgeContacts << " edge-edge contacts, expecting " << this->GetNumberOfEdgeEdgeResponses() << " responses.");
      this->ComputeEdgeEdgeResponses(r_responses.GetBuffer<ContactResponse>() + this->GetNumberOfNodeFacetResponses(), dpc_uNexts, dpc_uCurrs, eeStream);
    }

    if (numEdgeEdgeContacts > 0) {
      tledCheckCUDAErrors(cudaStreamSynchronize(eeStream));
      this->GetBVHTraverser().GetEdgeEdgeResults().ToggleActive();
      assert(!this->GetBVHTraverser().GetEdgeEdgeResults().IsActive());
      this->ReleaseEdgeEdgeResponseComputationMemory();
      tledCheckCUDAErrors(cudaStreamDestroy(eeStream));
    }

    if (numNodeFacetContacts > 0) {
      tledCheckCUDAErrors(cudaStreamSynchronize(nfStream));
      this->GetBVHTraverser().GetNodeFacetResults().ToggleActive();
      assert(!this->GetBVHTraverser().GetNodeFacetResults().IsActive());
      this->ReleaseNodeFacetResponseComputationMemory();
      tledCheckCUDAErrors(cudaStreamDestroy(nfStream));
    }

    this->ReleaseNodeVelocities();

    tledCUDADeviceMemoryBlock::CheckAllocationCounter();
    assert(r_responses.IsActive());
    this->ConsolidateResponses(dp_R, r_responses, this->GetNumberOfNodeFacetResponses() + this->GetNumberOfEdgeEdgeResponses());
    assert(!r_responses.IsActive());

    return true;
  }

  return false;
}

template <class TContactMesh, class TAPI>
void tledContactSolverImplGPU<TContactMesh, TAPI>::ComputeNodeVelocities(const float4 *dpc_uNexts, const float4 *dpc_uCurrs) {
  const int blockSize = 256;
  const int numBlocks = tledCUDAHelpers::GetNumberOfBlocks(this->GetMesh().GetNumberOfNodes(), blockSize);

  assert(mp_NodeVelocityBuffer == NULL);
  mp_NodeVelocityBuffer = &tledCUDADeviceMemoryBlock::GetNextFreeBufferWithSize<float3>(this->GetMesh().GetNumberOfNodes());

  tledContactSolverGPU_kernels::ComputeNodeVelocityKernel <<<numBlocks, blockSize>>> (mp_NodeVelocityBuffer->GetBuffer<float3>(), dpc_uNexts, dpc_uCurrs, this->GetHostGPUSurface().SurfaceToVolumeNodeIndexMap, this->GetMesh().GetNumberOfNodes());
  
  assert(mp_NodeVelocityBuffer->IsActive());
}

template <class TContactMesh, class TAPI>
void tledContactSolverImplGPU<TContactMesh, TAPI>::ReleaseNodeVelocities() {
  assert(mp_NodeVelocityBuffer->IsActive());
  mp_NodeVelocityBuffer->ToggleActive();
  mp_NodeVelocityBuffer = NULL;
}

#endif
