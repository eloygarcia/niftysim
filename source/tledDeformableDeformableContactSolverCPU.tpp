// =========================================================================
// File:       tledDeformableDeformableContactSolverCPU.tpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    October 2012
// 
// Copyright (c) 2011, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
template <class TContactMesh, class TAPI>
tledDeformableDeformableContactSolverImplCPU<TContactMesh, TAPI>::tledDeformableDeformableContactSolverImplCPU(tledUnstructuredContactManager &r_contactRes) : Superclass(r_contactRes) {
}

template <class TContactMesh, class TAPI>
tledBVHTraverserCPU* tledDeformableDeformableContactSolverImplCPU<TContactMesh, TAPI>::InstantiateBVHTraverser() {
  if (this->DoSelfCollision()) return tledSelfCollisionBVHTraverserCPU::CreateTraverser(this->GetManager().GetDeformableBVH());
  else return tledDeformableDeformableBVHTraverserCPU::CreateTraverser(this->GetManager().GetDeformableBVH());
}

template <class TContactMesh, class TAPI>
bool tledDeformableDeformableContactSolverImplCPU<TContactMesh, TAPI>::ComputeContactForces(float *p_f, const float uNexts[], const float uCurrs[]) {
  using namespace tledVectorArithmetic;

  std::vector<EdgeEdgeConstraintItem> edgeConstraints;
  std::vector<NodeFacetConstraintItem> nodeConstraints;
  
  for (std::vector<int>::const_iterator ic_slaveInd = this->GetBVHTraverser().GetSlaveEdgeIndices().begin(); ic_slaveInd < this->GetBVHTraverser().GetSlaveEdgeIndices().end(); ic_slaveInd++) {
    const tledBVHTraverserCPU::MasterEdgeContactData &edgeMaster = this->GetBVHTraverser().GetEdgeProjection(*ic_slaveInd);
    const std::pair<int, int> &slaveEdge = this->GetMesh().GetAllEdges()[*ic_slaveInd];
    const std::pair<int, int> &masterEdge = this->GetMesh().GetAllEdges()[edgeMaster.EdgeIndex];
    const float r = edgeMaster.CollisionCoords[1];
    const float q = edgeMaster.CollisionCoords[2];

    EdgeEdgeConstraintItem ci;

    assert(edgeMaster.CollisionCoords[0] < this->GetNodeCloseDistance() && r >= 0 && r <= 1 && q >= 0 && q <= 1);

    ci.SlaveNodeIndices[0] = slaveEdge.first, ci.SlaveNodeIndices[1] = slaveEdge.second;
    ci.MasterNodeIndices[0] = masterEdge.first, ci.MasterNodeIndices[1] = masterEdge.second;
    
    ci.SlaveShapeValues[0] = (1 - r), ci.SlaveShapeValues[1] = r;
    ci.MasterShapeValues[0] = (1 - q), ci.MasterShapeValues[1] = q;

    ci.GapValue = edgeMaster.CollisionCoords[0] - this->GetManager().GetSafetyMargin();
    std::copy(edgeMaster.Normal, edgeMaster.Normal + 3, ci.Normal);
    
    this->GetBVHTraverser().ResetEdgeProjection(*ic_slaveInd);    
    if (Superclass::template ComputeEdgeEdgeResponse<true, true>(ci) || Superclass::template ComputeEdgeEdgeRateResponse<true, true>(ci, uNexts, uCurrs)) {
      edgeConstraints.push_back(ci);
    } 
  }

  for (std::vector<int>::const_iterator ic_slaveInd = this->GetBVHTraverser().GetSlaveNodeIndices().begin(); ic_slaveInd < this->GetBVHTraverser().GetSlaveNodeIndices().end(); ic_slaveInd++) {
    const int contactNodeInd = *ic_slaveInd;    
    const tledBVHTraverserCPU::MasterFacetContactData &collPt = this->GetBVHTraverser().GetNodeProjection(*ic_slaveInd);
    const int minDistFacetInd = collPt.FacetIndex;
    const float minDist = collPt.CollisionCoords[2];
    const float minXi = collPt.CollisionCoords[0];
    const float minEta = collPt.CollisionCoords[1];

    NodeFacetConstraintItem ci;

    assert(minDist < this->GetNodeCloseDistance());
    assert(minDistFacetInd >= 0 && minDistFacetInd < this->GetMesh().GetNumberOfFacets());

    ci.ContactNodeIndices[0] = contactNodeInd;
    std::copy(this->GetMesh().GetFacet(minDistFacetInd).NodeIndices, this->GetMesh().GetFacet(minDistFacetInd).NodeIndices + Facet::NumberOfVertices, ci.ContactNodeIndices + 1);

    tledLogDebugStream(tledHelper::Info() << "Facet collision " << contactNodeInd << " - " << minDistFacetInd << " @ (" << minXi << ", " << minEta << ")\n"
		       << "Penetration depth: " << minDist);

    this->GetMesh().ComputeShapeValues(ci.ShapeValues, minXi, minEta);
    std::copy(collPt.Normal, collPt.Normal + 3, ci.Normal);

    this->GetBVHTraverser().ResetNodeProjection(*ic_slaveInd);

    ci.GapValue = minDist - this->GetManager().GetSafetyMargin();
    if (this->template ComputeNodeFacetResponse<true, true>(ci) || this->template ComputeNodeFacetRateResponse<true, true>(ci, uNexts, uCurrs)) {
      nodeConstraints.push_back(ci);
    } 
  } /* for surface nodes */

  this->template ApplyContactForces<true, true>(p_f, nodeConstraints, edgeConstraints);

  if (this->DoFriction()) {
    this->template ComputeNodeFacetFrictionResponse<true, true>(nodeConstraints, uNexts, uCurrs, this->GetMesh().GetFrictionCoefficient());
    this->template ComputeEdgeEdgeFrictionResponse<true, true>(edgeConstraints, uNexts, uCurrs, this->GetMesh().GetFrictionCoefficient());
    this->template ApplyFrictionForces<true, true>(p_f, nodeConstraints, edgeConstraints);
  }

  return nodeConstraints.size() > 0 || edgeConstraints.size() > 0;
}

