// =========================================================================
// File:       tledDeformableMovingRigidContactSolverCPU.tpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    December 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

template <class TDeformableContactMesh, class TRigidContactMesh, class TAPI>
tledBVHTraverserCPU* tledDeformableMovingRigidContactSolverImplCPU<TDeformableContactMesh, TRigidContactMesh, TAPI>::InstantiateBVHTraverser() {
  return tledDeformableMovingRigidBVHTraverserCPU::CreateTraverser(this->GetManager().GetDeformableBVH(), static_cast<const tledDynamicBVH&>(this->GetManager().GetRigidBVH(this->GetRigidSurfaceIndex())));
}

template <class TDeformableContactMesh, class TRigidContactMesh, class TAPI>
float* tledDeformableMovingRigidContactSolverImplCPU<TDeformableContactMesh, TRigidContactMesh, TAPI>::ComputeRigidSurfaceNodeVelocity(float *p_v, const int nodeIndex) const {
  if (this->GetRigidSurface().HasRotation()) {
    std::abort();
  } else {
    float tPrev[3];

    tledVectorArithmetic::ScalarDiv(tledVectorArithmetic::Sub(p_v, this->GetRigidSurface().GetTranslation(p_v), this->GetRigidSurface().GetTranslation(tPrev, this->GetRigidSurface().GetCurrentStep() - 1)), this->GetDt());
  }

  return p_v;
}

template <class TDeformableContactMesh, class TRigidContactMesh, class TAPI>
void tledDeformableMovingRigidContactSolverImplCPU<TDeformableContactMesh, TRigidContactMesh, TAPI>::ComputeRelativeNodeFacetVelocityMaster(float *p_v, const NodeFacetConstraintItem &ci, const float uNexts[], const float uCurrs[]) const {
  Superclass::ComputeRelativeNodeFacetVelocityMaster(p_v, ci, uNexts, uCurrs);
  for (int v = 0; v < RigidMesh::Facet::NumberOfVertices; v++) {
    float rigV[3];

    assert(ci.ContactNodeIndices[1+v] >= 0 && ci.ContactNodeIndices[1+v] < this->GetRigidSurface().GetNumberOfNodes());
    tledVectorArithmetic::Sub(p_v, p_v, tledVectorArithmetic::ScalarMul(this->ComputeRigidSurfaceNodeVelocity(rigV, ci.ContactNodeIndices[1+v]), ci.ShapeValues[v]));
  }
}

template <class TDeformableContactMesh, class TRigidContactMesh, class TAPI>
void tledDeformableMovingRigidContactSolverImplCPU<TDeformableContactMesh, TRigidContactMesh, TAPI>::ComputeRelativeNodeFacetVelocitySlave(float *p_v, const NodeFacetConstraintItem &ci, const float uNexts[], const float uCurrs[]) const {
  Superclass::ComputeRelativeNodeFacetVelocitySlave(p_v, ci, uNexts, uCurrs);

  {
    float rigV[3];

    assert(ci.ContactNodeIndices[0] >= 0 && ci.ContactNodeIndices[0] < this->GetRigidSurface().GetNumberOfNodes());
    tledVectorArithmetic::Add(p_v, p_v, this->ComputeRigidSurfaceNodeVelocity(rigV, ci.ContactNodeIndices[0]));
  }
}

template <class TDeformableContactMesh, class TRigidContactMesh, class TAPI>
void tledDeformableMovingRigidContactSolverImplCPU<TDeformableContactMesh, TRigidContactMesh, TAPI>::ComputeRelativeEdgeEdgeVelocityMaster(float *p_v, const EdgeEdgeConstraintItem &ci, const float uNexts[], const float uCurrs[]) const {
  Superclass::ComputeRelativeEdgeEdgeVelocityMaster(p_v, ci, uNexts, uCurrs);

  for (int i = 0; i < 2; i++) {
    float rigV[3];

    assert(ci.MasterNodeIndices[i] >= 0 && ci.MasterNodeIndices[i] < this->GetRigidSurface().GetNumberOfNodes());
    tledVectorArithmetic::ScalarMul(this->ComputeRigidSurfaceNodeVelocity(rigV, ci.MasterNodeIndices[i]), ci.MasterShapeValues[i]);
    tledVectorArithmetic::Sub(p_v, p_v, rigV);
  }
}

template <class TDeformableContactMesh, class TRigidContactMesh, class TAPI>
void tledDeformableMovingRigidContactSolverImplCPU<TDeformableContactMesh, TRigidContactMesh, TAPI>::ComputeRelativeEdgeEdgeVelocitySlave(float *p_v, const EdgeEdgeConstraintItem &ci, const float uNexts[], const float uCurrs[]) const {
  Superclass::ComputeRelativeEdgeEdgeVelocitySlave(p_v, ci, uNexts, uCurrs);

  for (int i = 0; i < 2; i++) {
    float rigV[3];

    assert(ci.SlaveNodeIndices[i] >= 0 && ci.SlaveNodeIndices[i] < this->GetRigidSurface().GetNumberOfNodes());
    tledVectorArithmetic::ScalarMul(this->ComputeRigidSurfaceNodeVelocity(rigV, ci.SlaveNodeIndices[i]), ci.SlaveShapeValues[i]);
    tledVectorArithmetic::Add(p_v, p_v, rigV);
  }
}
