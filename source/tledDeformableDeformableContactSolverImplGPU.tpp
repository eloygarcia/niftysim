// =========================================================================
// File:       tledDeformableDeformableContactSolverImplGPU.tpp
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
#include "tledCUDAHelpers.h"
#include "tledDeformableDeformableBVHTraverserGPU.cu"

#ifdef _TRACKING_CONTACTS
#include "tledTrackingDeformableDeformableBVHTraverserGPU.h"
#include "tledTrackingBVHTraverserGPU.h"

#include "tledTrackingDeformableDeformableBVHTraverserGPU.cu"
#include "tledTrackingBVHTraverserGPU.cu"
#endif

#include <thrust/transform.h>
#include <thrust/iterator/counting_iterator.h>

#include <cassert>

template <class TContactMesh, class TAPI>
class tledDeformableDeformableContactSolverImplGPU<TContactMesh, TAPI>::NodeProjectionToResponseTransform : public Superclass::NodeProjectionToResponseTransform {
protected:
  __device__ float3 ComputeRelativeNodeFacetVelocity(const NodeProjection &proj) const {
    float3 v = this->GetNodeVelocity(proj.ContactNodeIndices[0]);

    for (int n = 0; n < ContactMesh::Facet::NumberOfVertices; n++) {
      v -= proj.ShapeValues[n]*this->GetNodeVelocity(proj.ContactNodeIndices[1+n]);
    }

    return v/this->GetTimeStep();
  }

public:
  __device__ ContactResponse operator()(const NodeProjection &proj, const int nodeInd) const {
    const float effDepth = this->ComputeEffectiveDepth(proj);

    ContactResponse response;

    tledCudaAssert(nodeInd < ContactMesh::Facet::NumberOfVertices + 1 && nodeInd >= 0);
    response.NodeIndex = proj.ContactNodeIndices[nodeInd];
    tledCudaPrintf(0, "node-facet contact %d-(%d, %d, %d): depth = %f (thrsh = %f)\n", proj.ContactNodeIndices[0], proj.ContactNodeIndices[1], proj.ContactNodeIndices[2], proj.ContactNodeIndices[3], effDepth, this->GetRateResponseMaxDistance());
    if (effDepth < this->GetRateResponseMaxDistance()) {
      float lambda;

      if (effDepth < 0) {
	float beta;

	lambda = this->ComputeBaseLambdaFromDepth(effDepth);
	beta = this->template ComputeSlaveBeta<true, true>(proj);
	if (nodeInd > 0) beta = 1 - beta;
	lambda *= this->GetNodeMass(proj.ContactNodeIndices[nodeInd])*beta;
      } else {
	float vMag = dot(this->ComputeRelativeNodeFacetVelocity(proj), proj.Normal);

	lambda = -this->ComputeRateResponseWeight(effDepth)*this->template ComputeStoppingLambdaClamped<true, true>(vMag, proj); 
      }
      if (nodeInd > 0) {
	lambda *= -proj.ShapeValues[nodeInd-1]/this->ComputeShapeSquareSum(proj);
      }
      response.Force = lambda*proj.Normal;
    } else {
      response.Force = make_float3(0.f, 0.f, 0.f);
    }

    tledCudaPrintf(0, "Node-facet response (%s) %d (%d): %f, %f, %f\n", (nodeInd > 0? "master" : "slave"), response.NodeIndex, nodeInd, response.Force.x, response.Force.y, response.Force.z);

    return response;
  }

public:
  NodeProjectionToResponseTransform(const float *dpc_masses, const float3 *dpc_v, const float dt, const float off, const float rateDelta) : Superclass::NodeProjectionToResponseTransform(dpc_masses, dpc_v, dt, off, rateDelta) {}
};

template <class TContactMesh, class TAPI>
class tledDeformableDeformableContactSolverImplGPU<TContactMesh, TAPI>::EdgeProjectionToResponseTransform : public Superclass::EdgeProjectionToResponseTransform {
public:
  __device__ ContactResponse operator()(const EdgeProjection &proj, const int nodeInd) const {
    const float effDepth = proj.Xi.x - this->GetSafetyOffset();

    ContactResponse response;

    switch (nodeInd) {
    case 0:
      response.NodeIndex = proj.SlaveEdge.x;
      break;

    case 1:
      response.NodeIndex = proj.SlaveEdge.y;
      break;

    case 2:
      response.NodeIndex = proj.MasterEdge.x;
      break;

    default:
      response.NodeIndex = proj.MasterEdge.y;
    }

    tledCudaPrintf(2, "edge-edge contact (%d, %d)-(%d, %d): depth = %f (thrsh = %f)\n", proj.SlaveEdge.x, proj.SlaveEdge.y, proj.MasterEdge.x, proj.MasterEdge.y, proj.Xi.x, this->GetRateResponseMaxDistance());
    if (effDepth < this->GetRateResponseMaxDistance()) {
      const float4 shapeVals = this->template ComputeShapeValues<true, true>(proj); 

      if (effDepth >= 0) {
	float3 relV;
	float lambda, gRate;	

	relV = this->GetNodeVelocity(proj.SlaveEdge.x)*shapeVals.x + this->GetNodeVelocity(proj.SlaveEdge.y)*shapeVals.y;
	relV = relV - (this->GetNodeVelocity(proj.MasterEdge.x)*shapeVals.z + this->GetNodeVelocity(proj.MasterEdge.y)*shapeVals.w);
	gRate = dot(proj.Normal, relV);
	lambda = this->template ComputeStoppingLambdaClamped<true, true>(gRate/this->GetTimeStep(), shapeVals, proj);	  
	
	switch (nodeInd) {
	case 0:
	  lambda = -lambda*shapeVals.x/this->template ComputeShapeSquareSum<true, false>(shapeVals);
	  break;

	case 1:
	  lambda = -lambda*shapeVals.y/this->template ComputeShapeSquareSum<true, false>(shapeVals);
	  break;

	case 2:
	  lambda = lambda*shapeVals.z/this->template ComputeShapeSquareSum<false, true>(shapeVals);
	  break;

	default:
	  lambda = lambda*shapeVals.w/this->template ComputeShapeSquareSum<false, true>(shapeVals);
	  break;	   
	}
	response.Force = this->ComputeRateResponseWeight(effDepth)*lambda*proj.Normal;

	tledCudaPrintf(0, "Applying edge-edge rate response (%s) to %d: %f %f %f\n", (nodeInd < 2? "master" : "slave"), response.NodeIndex, response.Force.x, response.Force.y, response.Force.z);
      } else {
	float beta, shapeSqrSum, lambda;
	float2 masses;

	lambda = this->ComputeBaseLambdaFromDepth(effDepth);
	masses.x = shapeVals.x*this->GetNodeMass(proj.SlaveEdge.x) + shapeVals.y*this->GetNodeMass(proj.SlaveEdge.y);
	masses.y = shapeVals.z*this->GetNodeMass(proj.MasterEdge.x) + shapeVals.w*this->GetNodeMass(proj.MasterEdge.y);
	beta = masses.y/(masses.x + masses.y);
	if (nodeInd < 2) {
	  shapeSqrSum = this->template ComputeShapeSquareSum<true, false>(shapeVals);
	  lambda = (nodeInd == 0? shapeVals.x : shapeVals.y)*lambda;
	} else {
	  beta = 1 - beta;
	  shapeSqrSum = this->template ComputeShapeSquareSum<false, true>(shapeVals);
	  lambda = -(nodeInd == 2? shapeVals.z : shapeVals.w)*lambda;
	}
   
	lambda = this->GetNodeMass(response.NodeIndex)*beta*lambda/shapeSqrSum;
	response.Force = lambda*proj.Normal;
	tledCudaPrintf(0, "Edge-edge penetration response (%s) %d: %f, %f, %f\n", (nodeInd < 2? "master" : "slave"), response.NodeIndex, response.Force.x, response.Force.y, response.Force.z);
      }
    } else {
      response.Force = make_float3(0, 0, 0);
    }
    
    return response;
  }

public:
  EdgeProjectionToResponseTransform(const float *dpc_masses, const float3 *dpc_v, const float dt, const float off, const float rateThrsh) : Superclass::EdgeProjectionToResponseTransform(dpc_masses, dpc_v, dt, off, rateThrsh) {}
};

template <class TContactMesh, class TAPI>
tledBVHTraverserGPU* tledDeformableDeformableContactSolverImplGPU<TContactMesh, TAPI>::InstantiateBVHTraverser() {
  if (this->DoSelfCollision()) {
    tledFatalNotYetImplementedError;

    return NULL;
  } else {
#ifdef _TRACKING_CONTACTS
    return tledTrackingDeformableDeformableBVHTraverserGPU::CreateTraverser(this->GetManager().GetDeformableBVH());
#else
    return tledDeformableDeformableBVHTraverserGPU::CreateTraverser(this->GetManager().GetDeformableBVH());
#endif
  }
}

template <class TContactMesh, class TAPI>
void tledDeformableDeformableContactSolverImplGPU<TContactMesh, TAPI>::ComputeNodeFacetResponses(ContactResponse *dp_dst, const float4 *dpc_uNexts, const float4 *dpc_uCurrs, const cudaStream_t nodeFacetStream) {
  const tledDeformableDeformableBVHTraverserGPU &traverser = static_cast<const tledDeformableDeformableBVHTraverserGPU&>(this->GetBVHTraverser());
  const int blockSize = 128;
  const int numContacts = traverser.GetNumberOfNodeFacetContacts();
  const int numResponsesPerProjection = 1 + ContactMesh::Facet::NumberOfVertices;
  const int numNodeFacetResponses = numResponsesPerProjection*numContacts;
  const int numBlocks = tledCUDAHelpers::GetNumberOfBlocks(numNodeFacetResponses, blockSize);	
  const float *dpc_nodeMasses = this->GetHostGPUSurface().NodeMasses;
  const tledCUDADeviceMemoryBlock &projections = traverser.GetNodeFacetResults();
  const NodeProjectionToResponseTransform xForm(dpc_nodeMasses, this->GetNodeVelocities(), this->GetDt(), this->GetManager().GetSafetyMargin(), this->GetNodeCloseDistance());

  assert(projections.IsActive());
  tledContactSolverGPU_kernels::ProjectionToResponseTransformKernel<NodeProjection, NodeProjectionToResponseTransform> <<<numBlocks, blockSize, 0, nodeFacetStream>>> (dp_dst, projections.GetBuffer<NodeProjection>(), numContacts, xForm, numResponsesPerProjection);
}

template <class TContactMesh, class TAPI>
void tledDeformableDeformableContactSolverImplGPU<TContactMesh, TAPI>::ComputeEdgeEdgeResponses(ContactResponse *dp_dst, const float4 *dpc_uNexts, const float4 *dpc_uCurrs, const cudaStream_t stream) {
  const tledDeformableDeformableBVHTraverserGPU &traverser = static_cast<const tledDeformableDeformableBVHTraverserGPU&>(this->GetBVHTraverser());
  const tledCUDADeviceMemoryBlock &projections = traverser.GetEdgeEdgeResults();
  const int blockSize = 128;
  const int numContacts = traverser.GetNumberOfEdgeEdgeContacts();
  const int numBlocks = tledCUDAHelpers::GetNumberOfBlocks(numContacts*4, blockSize);	
  const float *dpc_nodeMasses = this->GetHostGPUSurface().NodeMasses;
  const EdgeProjectionToResponseTransform xForm(dpc_nodeMasses, this->GetNodeVelocities(), this->GetDt(), this->GetManager().GetSafetyMargin(), this->GetNodeCloseDistance());
	
  tledContactSolverGPU_kernels::ProjectionToResponseTransformKernel<EdgeProjection, EdgeProjectionToResponseTransform> <<<numBlocks, blockSize, 0, stream>>> (dp_dst, projections.GetBuffer<EdgeProjection>(), numContacts, xForm, 4);
}

template <class TContactMesh, class TAPI>
void tledDeformableDeformableContactSolverImplGPU<TContactMesh, TAPI>::ReleaseNodeFacetResponseComputationMemory() {
}

template <class TContactMesh, class TAPI>
void tledDeformableDeformableContactSolverImplGPU<TContactMesh, TAPI>::ReleaseEdgeEdgeResponseComputationMemory() {
}
