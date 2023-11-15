// =========================================================================
// File:       tledDeformableRigidContactSolverCPU.tpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    October 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

template <class TDeformableContactMesh, class TRigidContactMesh, class TAPI>
tledDeformableRigidContactSolverImplCPU<TDeformableContactMesh, TRigidContactMesh, TAPI>::tledDeformableRigidContactSolverImplCPU(tledUnstructuredContactManager &r_contactManager, const int rigidSurfaceIndex) : Superclass(r_contactManager) {  
  this->SetRigidSurfaceIndex(rigidSurfaceIndex);
  mpc_RigidSurface = &r_contactManager.GetRigidSurface<TRigidContactMesh>(rigidSurfaceIndex);    
}
      
template <class TDeformableContactMesh, class TRigidContactMesh, class TAPI>
tledBVHTraverserCPU* tledDeformableRigidContactSolverImplCPU<TDeformableContactMesh, TRigidContactMesh, TAPI>::InstantiateBVHTraverser() {
  return tledDeformableRigidBVHTraverserCPU::CreateTraverser(this->GetManager().GetDeformableBVH(), this->GetManager().GetRigidBVH(this->GetRigidSurfaceIndex()));
}

template <class TDeformableContactMesh, class TRigidContactMesh, class TAPI>
bool tledDeformableRigidContactSolverImplCPU<TDeformableContactMesh, TRigidContactMesh, TAPI>::ComputeContactForces(float *p_f, const float uNexts[], const float uCurrs[]) {
  if (this->DoMaster()) return this->ComputeSlaveResponses(p_f, uNexts, uCurrs);    
  else return this->ComputeMasterResponses(p_f, uNexts, uCurrs);
}
 
template <class TDeformableContactMesh, class TRigidContactMesh, class TAPI>
bool tledDeformableRigidContactSolverImplCPU<TDeformableContactMesh, TRigidContactMesh, TAPI>::ComputeMasterResponses(float *p_fs, const float uNexts[], const float uCurrs[]) {
  using namespace tledVectorArithmetic;

  std::vector<NodeFacetConstraintItem> nodeConstraints;
  std::vector<EdgeEdgeConstraintItem> edgeConstraints;

  for (std::vector<int>::const_iterator ic_slaveNodeInd = this->GetBVHTraverser().GetSlaveNodeIndices().begin(); ic_slaveNodeInd < this->GetBVHTraverser().GetSlaveNodeIndices().end(); ic_slaveNodeInd++) {	
    const Facet &mFacet = this->GetDeformableSurface().GetFacet(this->GetBVHTraverser().GetNodeProjection(*ic_slaveNodeInd).FacetIndex);
    const float *minXi = this->GetBVHTraverser().GetNodeProjection(*ic_slaveNodeInd).CollisionCoords;

    NodeFacetConstraintItem ci;

    assert(this->GetBVHTraverser().GetNodeProjection(*ic_slaveNodeInd).FacetIndex >= 0 && this->GetBVHTraverser().GetNodeProjection(*ic_slaveNodeInd).FacetIndex < GetDeformableSurface().GetNumberOfFacets());

    ci.ContactNodeIndices[0] = *ic_slaveNodeInd;
    std::copy(mFacet.NodeIndices, mFacet.NodeIndices + Facet::NumberOfVertices, ci.ContactNodeIndices + 1);

    assert(Facet::NumberOfVertices != 3 || (minXi[0] >= 0 && minXi[1] >= 0 && minXi[0] + minXi[1] < 1 + 1e-5));
    this->GetDeformableSurface().ComputeShapeValues(ci.ShapeValues, minXi[0], minXi[1]);
    std::copy(this->GetBVHTraverser().GetNodeProjection(*ic_slaveNodeInd).Normal, this->GetBVHTraverser().GetNodeProjection(*ic_slaveNodeInd).Normal + 3, ci.Normal);

    ci.GapValue = minXi[2] - this->GetManager().GetSafetyMargin();

    this->GetBVHTraverser().ResetNodeProjection(*ic_slaveNodeInd);
    if (this->template ComputeNodeFacetResponse<true, false>(ci) || this->template ComputeNodeFacetRateResponse<true, false>(ci, uNexts, uCurrs)) {
      nodeConstraints.push_back(ci);
    } 
  } /* for active slave nodes */

  for (std::vector<int>::const_iterator ic_slaveEdgeInd = this->GetBVHTraverser().GetSlaveEdgeIndices().begin(); ic_slaveEdgeInd < this->GetBVHTraverser().GetSlaveEdgeIndices().end(); ic_slaveEdgeInd++) {
    const std::pair<int, int> &masterEdge = this->GetDeformableSurface().GetEdge(this->GetBVHTraverser().GetEdgeProjection(*ic_slaveEdgeInd).EdgeIndex);
    const std::pair<int, int> &slaveEdge = this->GetRigidSurface().GetEdge(*ic_slaveEdgeInd);
    const float r = this->GetBVHTraverser().GetEdgeProjection(*ic_slaveEdgeInd).CollisionCoords[1];
    const float q = this->GetBVHTraverser().GetEdgeProjection(*ic_slaveEdgeInd).CollisionCoords[2];

    EdgeEdgeConstraintItem ci;

    ci.SlaveNodeIndices[0] = slaveEdge.first, ci.SlaveNodeIndices[1] = slaveEdge.second;
    ci.MasterNodeIndices[0] = masterEdge.first, ci.MasterNodeIndices[1] = masterEdge.second;

    ci.SlaveShapeValues[0] = (1 - r), ci.SlaveShapeValues[1] = r;
    ci.MasterShapeValues[0] = (1 - q), ci.MasterShapeValues[1] = q;

    ci.GapValue = this->GetBVHTraverser().GetEdgeProjection(*ic_slaveEdgeInd).CollisionCoords[0] - this->GetManager().GetSafetyMargin();
    std::copy(this->GetBVHTraverser().GetEdgeProjection(*ic_slaveEdgeInd).Normal, this->GetBVHTraverser().GetEdgeProjection(*ic_slaveEdgeInd).Normal + 3, ci.Normal);

    this->GetBVHTraverser().ResetEdgeProjection(*ic_slaveEdgeInd);

    if (this->template ComputeEdgeEdgeResponse<true, false>(ci) || this->template ComputeEdgeEdgeRateResponse<true, false>(ci, uNexts, uCurrs)) {
      edgeConstraints.push_back(ci);
    }
  }

  this->template ApplyContactForces<true, false>(p_fs, nodeConstraints, edgeConstraints);

  if (this->DoFriction()) {
    this->template ComputeNodeFacetFrictionResponse<true, false>(nodeConstraints, uNexts, uCurrs, this->GetRigidSurface().GetFrictionCoefficient());
    this->template ComputeEdgeEdgeFrictionResponse<true, false>(edgeConstraints, uNexts, uCurrs, this->GetRigidSurface().GetFrictionCoefficient());
    this->template ApplyFrictionForces<true, false>(p_fs, nodeConstraints, edgeConstraints);
  }

  return nodeConstraints.size() > 0 || edgeConstraints.size() > 0;
} /* ComputeMasterResponses */

template <class TDeformableContactMesh, class TRigidContactMesh, class TAPI>
bool tledDeformableRigidContactSolverImplCPU<TDeformableContactMesh, TRigidContactMesh, TAPI>::ComputeSlaveResponses(float *p_fs, const float uNexts[], const float uCurrs[]) {
  using namespace tledVectorArithmetic;

  std::vector<NodeFacetConstraintItem> nodeConstraints;
  std::vector<EdgeEdgeConstraintItem> edgeConstraints;

  for (std::vector<int>::const_iterator ic_slaveInd = this->GetBVHTraverser().GetSlaveNodeIndices().begin(); ic_slaveInd < this->GetBVHTraverser().GetSlaveNodeIndices().end(); ic_slaveInd++) {
    const int minMFacetInd = this->GetBVHTraverser().GetNodeProjection(*ic_slaveInd).FacetIndex;
    const float minXi = this->GetBVHTraverser().GetNodeProjection(*ic_slaveInd).CollisionCoords[0];
    const float minEta = this->GetBVHTraverser().GetNodeProjection(*ic_slaveInd).CollisionCoords[1];
    const float minDepth = this->GetBVHTraverser().GetNodeProjection(*ic_slaveInd).CollisionCoords[2];

    NodeFacetConstraintItem ci;

    assert(minMFacetInd >= 0 && minMFacetInd < this->GetRigidSurface().GetNumberOfFacets());
    assert(*ic_slaveInd >= 0 && *ic_slaveInd < this->GetDeformableSurface().GetNumberOfNodes());
    assert(this->GetRigidSurface().GetSlaveNodeMask()[*ic_slaveInd]);

    ci.ContactNodeIndices[0] = *ic_slaveInd;      
    std::copy(this->GetRigidSurface().GetFacet(minMFacetInd).NodeIndices, this->GetRigidSurface().GetFacet(minMFacetInd).NodeIndices + RigidMesh::Facet::NumberOfVertices, ci.ContactNodeIndices + 1);
    this->GetMesh().ComputeShapeValues(ci.ShapeValues, minXi, minEta);
    this->GetRigidSurface().ComputeC0Normal(ci.Normal, this->GetRigidSurface().GetFacet(minMFacetInd).NodeIndices, ci.ShapeValues);
      
    ci.GapValue = minDepth - this->GetManager().GetSafetyMargin();
    this->GetBVHTraverser().ResetNodeProjection(*ic_slaveInd);
    if (this->template ComputeNodeFacetResponse<false, true>(ci) || this->template ComputeNodeFacetRateResponse<false, true>(ci, uNexts, uCurrs)) {
      nodeConstraints.push_back(ci);
    } 
  }

  for (std::vector<int>::const_iterator ic_slaveEdgeInd = this->GetBVHTraverser().GetSlaveEdgeIndices().begin(); ic_slaveEdgeInd < this->GetBVHTraverser().GetSlaveEdgeIndices().end(); ic_slaveEdgeInd++) {
    const std::pair<int, int> &slaveEdge = this->GetDeformableSurface().GetEdge(*ic_slaveEdgeInd);
    const std::pair<int, int> &masterEdge = this->GetRigidSurface().GetEdge(this->GetBVHTraverser().GetEdgeProjection(*ic_slaveEdgeInd).EdgeIndex);
    const float r = this->GetBVHTraverser().GetEdgeProjection(*ic_slaveEdgeInd).CollisionCoords[1];
    const float q = this->GetBVHTraverser().GetEdgeProjection(*ic_slaveEdgeInd).CollisionCoords[2];

    EdgeEdgeConstraintItem ci;

    assert(this->GetBVHTraverser().GetEdgeProjection(*ic_slaveEdgeInd).EdgeIndex >= 0 && this->GetBVHTraverser().GetEdgeProjection(*ic_slaveEdgeInd).EdgeIndex < GetRigidSurface().GetNumberOfEdges());
    ci.SlaveNodeIndices[0] = slaveEdge.first, ci.SlaveNodeIndices[1] = slaveEdge.second;
    ci.MasterNodeIndices[0] = masterEdge.first, ci.MasterNodeIndices[1] = masterEdge.second;

    ci.SlaveShapeValues[0] = (1 - r), ci.SlaveShapeValues[1] = r;
    ci.MasterShapeValues[0] = (1 - q), ci.MasterShapeValues[1] = q;
    
    std::copy(this->GetBVHTraverser().GetEdgeProjection(*ic_slaveEdgeInd).Normal, this->GetBVHTraverser().GetEdgeProjection(*ic_slaveEdgeInd).Normal + 3, ci.Normal);
    ci.GapValue = this->GetBVHTraverser().GetEdgeProjection(*ic_slaveEdgeInd).CollisionCoords[0] - this->GetManager().GetSafetyMargin();

    this->GetBVHTraverser().ResetEdgeProjection(*ic_slaveEdgeInd);
    if (this->template ComputeEdgeEdgeResponse<false, true>(ci) || this->template ComputeEdgeEdgeRateResponse<false, true>(ci, uNexts, uCurrs)) {
      edgeConstraints.push_back(ci);
    }
  }

  this->template ApplyContactForces<false, true>(p_fs, nodeConstraints, edgeConstraints);

  if (this->DoFriction()) {
    this->template ComputeNodeFacetFrictionResponse<false, true>(nodeConstraints, uNexts, uCurrs, this->GetRigidSurface().GetFrictionCoefficient());
    this->template ComputeEdgeEdgeFrictionResponse<false, true>(edgeConstraints, uNexts, uCurrs, this->GetRigidSurface().GetFrictionCoefficient());
    this->template ApplyFrictionForces<false, true>(p_fs, nodeConstraints, edgeConstraints);
  }

  return nodeConstraints.size() > 0 || edgeConstraints.size() > 0;
}
