// =========================================================================
// File:       tledDeformableMovingRigidContactSolverGPU.tpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   CUDA
// Created:    September 2014
// 
// Copyright (c) 2014, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

#include "tledCUDAHelpers.h"

template <class TDeformableContactMesh, class TRigidContactMesh, class TAPI>
class tledDeformableMovingRigidContactSolverImplGPU<TDeformableContactMesh, TRigidContactMesh, TAPI>::_MovingRigidTransformBase {
private:
  const float3 *mdpc_RigXCurr;
  const float3 mc_Translation, mc_Rotations, mc_COR;

protected:
  __device__ const float3& GetNodeCoordinates(const int nodeInd) const { return mdpc_RigXCurr[nodeInd]; }
  __device__ const float3& GetTranslation(void) const { return mc_Translation; }
  __device__ const float3& GetCOR(void) const { return mc_COR; }
  __device__ bool HasRotation(void) const { return fabs(mc_Rotations.x) + fabs(mc_Rotations.y) + fabs(mc_Rotations.z) > 0.f; }
  __device__ const float3& GetRotations(void) const { return mc_Rotations; }

public:
  /* Pass per time-step rotation, translation for direct application */
  _MovingRigidTransformBase(const float3 *dpc_rigXCurr, const float3 &t, const float3 &rots, const float3 &cor) : mdpc_RigXCurr(dpc_rigXCurr), mc_Translation(t), mc_Rotations(rots), mc_COR(cor) {}
};

template <class TDeformableContactMesh, class TRigidContactMesh, class TAPI>
class tledDeformableMovingRigidContactSolverImplGPU<TDeformableContactMesh, TRigidContactMesh, TAPI>::NodeProjectionToResponseRigidMasterTransform : public Superclass::NodeProjectionToResponseTransform, public _MovingRigidTransformBase {
protected:
  __device__ float3 ComputeRelativeNodeFacetVelocity(const NodeProjection &proj) const {
    float3 v = this->GetNodeVelocity(proj.ContactNodeIndices[0]);

    /*
     * Only per time-step translation/rotation assumed.
     * Translation uniform over whole surface.
     */
    if (this->HasRotation()) {
      tledCudaAssert(false);
    }
    v += this->GetTranslation();

    return v;
  }

public:
  __device__ ContactResponse operator()(const NodeProjection &proj, const int nodeInd) const {
    float effDepth = this->ComputeEffectiveDepth(proj), lambda;
    ContactResponse response;

    tledCudaAssert(nodeInd == 0);
    response.NodeIndex = proj.ContactNodeIndices[0];
    if (effDepth < 0) {
      lambda = this->ComputeBaseLambdaFromDepth(effDepth);
      lambda *= this->GetNodeMass(proj.ContactNodeIndices[0]);
      tledCudaPrintf(1, "Applying node-facet penetration response (slave %d, g = %f): base mag = %f\n", response.NodeIndex, effDepth, lambda);
    } else {
      float vMag = dot(this->ComputeRelativeNodeFacetVelocity(proj), proj.Normal);

      lambda = this->ComputeRateResponseWeight(effDepth)*this->template ComputeStoppingLambdaClamped<true, false>(vMag, proj)/this->GetTimeStep();       
      tledCudaPrintf(1, "Applying node-facet rate response (slave %d, g = %f): base mag = %f\n", response.NodeIndex, effDepth, lambda);
    }

    response.Force = lambda*proj.Normal;
    tledCudaPrintf(1, "Applying node-facet contact force (slave) to %d: %f %f %f\n", response.NodeIndex, response.Force.x, response.Force.y, response.Force.z);

    return response;
  }

public:
  NodeProjectionToResponseRigidMasterTransform(const float *dpc_masses, const float3 *dpc_v, const float3 *dpc_rigXCurr, const float3 &t, const float3 &rots, const float3 &cor, const float dt, const float off, const float rateDelta) 
    : Superclass::NodeProjectionToResponseTransform(dpc_masses, dpc_v, dt, off, rateDelta), _MovingRigidTransformBase(dpc_rigXCurr, t, rots, cor) {}
};

template <class TDeformableContactMesh, class TRigidContactMesh, class TAPI>
class tledDeformableMovingRigidContactSolverImplGPU<TDeformableContactMesh, TRigidContactMesh, TAPI>::NodeProjectionToResponseDeformableMasterTransform : public Superclass::NodeProjectionToResponseTransform, public _MovingRigidTransformBase {
protected:
  __device__ float3 ComputeRelativeNodeFacetVelocity(const NodeProjection &proj) const {
    float3 v;

    if (this->HasRotation()) {
      tledCudaAssert(false);
    } else {
      v = this->GetTranslation();
    }

    for (int n = 0; n < ContactMesh::Facet::NumberOfVertices; n++) {
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

      lambda = this->ComputeRateResponseWeight(effDepth)*this->template ComputeStoppingLambdaClamped<false, true>(vMag, proj)/this->GetTimeStep(); 
    }

    lambda *= proj.ShapeValues[nodeInd]/this->ComputeShapeSquareSum(proj);    
    response.Force = lambda*proj.Normal;
    tledCudaPrintf(1, "Applying node-facet contact force (master) to %d: %f %f %f\n", response.NodeIndex, response.Force.x, response.Force.y, response.Force.z);

    return response;
  }

public:
  NodeProjectionToResponseDeformableMasterTransform(const float *dpc_masses, const float3 *dpc_v, const float3 *dpc_rigXCurr, const float3 &t, const float3 &rots, const float3 &cor, const float dt, const float off, const float rateDelta) 
    : Superclass::NodeProjectionToResponseTransform(dpc_masses, dpc_v, dt, off, rateDelta), _MovingRigidTransformBase(dpc_rigXCurr, t, rots, cor) {}
};

template <class TDeformableContactMesh, class TRigidContactMesh, class TAPI>
template <const bool t_isDeformableMaster>
class tledDeformableMovingRigidContactSolverImplGPU<TDeformableContactMesh, TRigidContactMesh, TAPI>::EdgeProjectionToResponseTransform : public Superclass::template EdgeProjectionToResponseTransform<t_isDeformableMaster>, public _MovingRigidTransformBase {
protected:
  __device__ float ComputeRelativeEdgeEdgeVelocityMagnitude(const float4 &shapeVals, const EdgeProjection &proj) const {
    float3 v;

    if (t_isDeformableMaster) {
      v = this->GetTranslation();
      v -= shapeVals.z*this->GetNodeVelocity(proj.MasterEdge.x) + shapeVals.w*this->GetNodeVelocity(proj.MasterEdge.y);
    } else {
      v = shapeVals.x*this->GetNodeVelocity(proj.SlaveEdge.x) + shapeVals.y*this->GetNodeVelocity(proj.SlaveEdge.y);
      v -= this->GetTranslation();
    }

    return dot(v, proj.Normal)/this->GetTimeStep();
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
    tledCudaPrintf(1, "Applying edge-edge contact force (%s) to %d: %f %f %f\n", (t_isDeformableMaster? "master" : "slave"), response.NodeIndex, response.Force.x, response.Force.y, response.Force.z);
    
    return response;
  }

public:
  EdgeProjectionToResponseTransform(const float *dpc_masses, const float3 *dpc_v, const float3 *dpc_rigXCurr, const float3 &t, const float3 &rots, const float3 &cor, const float dt, const float off, const float rateThrsh) 
    : Superclass::template EdgeProjectionToResponseTransform<t_isDeformableMaster>(dpc_masses, dpc_v, dt, off, rateThrsh), _MovingRigidTransformBase(dpc_rigXCurr, t, rots, cor) {}
};

template <class TDeformableContactMesh, class TRigidContactMesh, class TAPI>
float3 tledDeformableMovingRigidContactSolverImplGPU<TDeformableContactMesh, TRigidContactMesh, TAPI>::_MakeTimeStepIncrement(const float xform[], const int numSteps) {
  float sXform[3];

  return tledCUDAHelpers::ConvertToFloat3(tledVectorArithmetic::ScalarDiv(sXform, xform, numSteps));
}

template <class TDeformableContactMesh, class TRigidContactMesh, class TAPI>
tledBVHTraverserGPU* tledDeformableMovingRigidContactSolverImplGPU<TDeformableContactMesh, TRigidContactMesh, TAPI>::InstantiateBVHTraverser() {
  return tledDeformableMovingRigidBVHTraverserGPU::CreateTraverser(this->GetManager().GetDeformableBVH(), static_cast<const tledDynamicBVH&>(this->GetManager().GetRigidBVH(this->GetRigidSurfaceIndex())));
}

template <class TDeformableContactMesh, class TRigidContactMesh, class TAPI>
void tledDeformableMovingRigidContactSolverImplGPU<TDeformableContactMesh, TRigidContactMesh, TAPI>::ComputeNodeFacetResponsesRigidMaster(ContactResponse *dp_dst, const float4 *dpc_uNexts, const float4 *dpc_uCurrs, const cudaStream_t stream) {
  const tledCUDADeviceMemoryBlock &projections = this->GetBVHTraverser().GetNodeFacetResults();
  const float *dpc_nodeMasses = this->GetHostGPUSurface().NodeMasses;
  const int blockSize = 128;
  const int numContacts = this->GetBVHTraverser().GetNumberOfNodeFacetContacts();
  const int numBlocks = tledCUDAHelpers::GetNumberOfBlocks(numContacts, blockSize);  
  const NodeProjectionToResponseRigidMasterTransform xForm(dpc_nodeMasses, this->GetNodeVelocities(), this->GetRigidSurface().GetAllOnDeviceNodeCoordinates(), _MakeTimeStepIncrement(this->GetRigidSurface().GetTotalTranslation(), this->GetRigidSurface().GetTotalNumberOfSteps()), _MakeTimeStepIncrement(this->GetRigidSurface().GetAllTotalRotationAngles(), this->GetRigidSurface().GetTotalNumberOfSteps()), tledCUDAHelpers::ConvertToFloat3(this->GetRigidSurface().GetRotationCentre()), this->GetDt(), this->GetManager().GetSafetyMargin(), this->GetNodeCloseDistance());

  tledContactSolverGPU_kernels::ProjectionToResponseTransformKernel<NodeProjection, NodeProjectionToResponseRigidMasterTransform> <<<numBlocks, blockSize, 0, stream>>> (dp_dst, projections.GetBuffer<NodeProjection>(), numContacts, xForm, 1);
  tledDeviceSyncDebug;
}

template <class TDeformableContactMesh, class TRigidContactMesh, class TAPI>
void tledDeformableMovingRigidContactSolverImplGPU<TDeformableContactMesh, TRigidContactMesh, TAPI>::ComputeNodeFacetResponsesDeformableMaster(ContactResponse *dp_dst, const float4 *dpc_uNexts, const float4 *dpc_uCurrs, const cudaStream_t stream) {
  const tledCUDADeviceMemoryBlock &projections = this->GetBVHTraverser().GetNodeFacetResults();
  const float *dpc_nodeMasses = this->GetHostGPUSurface().NodeMasses;
  const int numResponsesPerContact = TDeformableContactMesh::Facet::NumberOfVertices;
  const int blockSize = 128;
  const int numContacts = this->GetBVHTraverser().GetNumberOfNodeFacetContacts();
  const int numBlocks = tledCUDAHelpers::GetNumberOfBlocks(numContacts*numResponsesPerContact, blockSize);
  const NodeProjectionToResponseDeformableMasterTransform xForm(dpc_nodeMasses, this->GetNodeVelocities(), this->GetRigidSurface().GetAllOnDeviceNodeCoordinates(), _MakeTimeStepIncrement(this->GetRigidSurface().GetTotalTranslation(), this->GetRigidSurface().GetTotalNumberOfSteps()), _MakeTimeStepIncrement(this->GetRigidSurface().GetAllTotalRotationAngles(), this->GetRigidSurface().GetTotalNumberOfSteps()), tledCUDAHelpers::ConvertToFloat3(this->GetRigidSurface().GetRotationCentre()), this->GetDt(), this->GetManager().GetSafetyMargin(), this->GetNodeCloseDistance());

  tledContactSolverGPU_kernels::ProjectionToResponseTransformKernel<NodeProjection, NodeProjectionToResponseDeformableMasterTransform> <<<blockSize, numBlocks, 0, stream>>> (dp_dst, projections.GetBuffer<NodeProjection>(), numContacts, xForm, numResponsesPerContact);
  tledDeviceSyncDebug;
}

template <class TDeformableContactMesh, class TRigidContactMesh, class TAPI>
void tledDeformableMovingRigidContactSolverImplGPU<TDeformableContactMesh, TRigidContactMesh, TAPI>::ComputeEdgeEdgeResponsesDeformableMaster(ContactResponse *dp_dst, const float4 *dpc_uNexts, const float4 *dpc_uCurrs, const cudaStream_t stream) {
  typedef EdgeProjectionToResponseTransform<true> __XForm;

  const tledCUDADeviceMemoryBlock &projections = this->GetBVHTraverser().GetEdgeEdgeResults();
  const int blockSize = 128;
  const int numBlocks = tledCUDAHelpers::GetNumberOfBlocks(this->GetBVHTraverser().GetNumberOfEdgeEdgeContacts()*2, blockSize);	
  const float *dpc_nodeMasses = this->GetHostGPUSurface().NodeMasses;
  const __XForm xForm(dpc_nodeMasses, this->GetNodeVelocities(), this->GetRigidSurface().GetAllOnDeviceNodeCoordinates(), _MakeTimeStepIncrement(this->GetRigidSurface().GetTotalTranslation(), this->GetRigidSurface().GetTotalNumberOfSteps()), _MakeTimeStepIncrement(this->GetRigidSurface().GetAllTotalRotationAngles(), this->GetRigidSurface().GetTotalNumberOfSteps()), tledCUDAHelpers::ConvertToFloat3(this->GetRigidSurface().GetRotationCentre()), this->GetDt(), this->GetManager().GetSafetyMargin(), this->GetNodeCloseDistance());

  assert(projections.IsActive());
  tledContactSolverGPU_kernels::ProjectionToResponseTransformKernel<EdgeProjection, __XForm> <<<numBlocks, blockSize, 0, stream>>> (dp_dst, projections.GetBuffer<EdgeProjection>(), this->GetBVHTraverser().GetNumberOfEdgeEdgeContacts(), xForm, 2);
}

template <class TDeformableContactMesh, class TRigidContactMesh, class TAPI>
void tledDeformableMovingRigidContactSolverImplGPU<TDeformableContactMesh, TRigidContactMesh, TAPI>::ComputeEdgeEdgeResponsesRigidMaster(ContactResponse *dp_dst, const float4 *dpc_uNexts, const float4 *dpc_uCurrs, const cudaStream_t stream) {
  typedef EdgeProjectionToResponseTransform<false> __XForm;

  const tledCUDADeviceMemoryBlock &projections = this->GetBVHTraverser().GetEdgeEdgeResults();
  const float *dpc_nodeMasses = this->GetHostGPUSurface().NodeMasses;
  const int blockSize = 128;
  const int numBlocks = tledCUDAHelpers::GetNumberOfBlocks(this->GetBVHTraverser().GetNumberOfEdgeEdgeContacts()*2, blockSize);	
  const __XForm xForm(dpc_nodeMasses, this->GetNodeVelocities(), this->GetRigidSurface().GetAllOnDeviceNodeCoordinates(), _MakeTimeStepIncrement(this->GetRigidSurface().GetTotalTranslation(), this->GetRigidSurface().GetTotalNumberOfSteps()), _MakeTimeStepIncrement(this->GetRigidSurface().GetAllTotalRotationAngles(), this->GetRigidSurface().GetTotalNumberOfSteps()), tledCUDAHelpers::ConvertToFloat3(this->GetRigidSurface().GetRotationCentre()), this->GetDt(), this->GetManager().GetSafetyMargin(), this->GetNodeCloseDistance());

  assert(projections.IsActive());
  tledContactSolverGPU_kernels::ProjectionToResponseTransformKernel<EdgeProjection, __XForm> <<<numBlocks, blockSize, 0, stream>>> (dp_dst, projections.GetBuffer<EdgeProjection>(), this->GetBVHTraverser().GetNumberOfEdgeEdgeContacts(), xForm, 2);
}
