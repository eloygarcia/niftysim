// =========================================================================
// File:       tledDeformableDeformableBVHTraverserCPU.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    November 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

template <class TBVH, class TAPI>
bool tledDeformableDeformableBVHTraverserImplCPU<TBVH, TAPI>::DoNarrowPhaseNodeFacetDetection(const int nodeInd, const int primInd) {
  assert(std::find(this->GetMasterMesh().GetFacet(primInd).NodeIndices, this->GetMasterMesh().GetFacet(primInd).NodeIndices + MasterMesh::Facet::NumberOfVertices, nodeInd) == this->GetMasterMesh().GetFacet(primInd).NodeIndices + MasterMesh::Facet::NumberOfVertices);
  assert(&this->GetMasterMesh() == &this->GetSlaveMesh());
  assert(this->GetNodeProjectionBuffer()[nodeInd].FacetIndex != primInd);
  if (this->DoBehindFacetTest(this->GetSlaveMesh().GetOldNodeCoordinates(nodeInd), this->GetSlaveMesh().GetNormalisedOldFacetNormalCached(primInd), this->GetMasterMesh().GetOldNodeCoordinates(this->GetMasterMesh().GetFacet(primInd).NodeIndices[0]))) {
    float xi[3];

    if (this->DoNodeFacetInitialProjection(xi, this->GetSlaveMesh().GetNodeCoordinates(nodeInd), this->GetSlaveMesh().GetFacetProjectionOperatorCached(primInd), this->GetNodeProjectionBuffer()[nodeInd].CollisionCoords)) {
      return this->DoNarrowPhaseNodeFacetIterativeProjection(xi, nodeInd, primInd);
    } /* if not was behind facet */
  } /* if not incident facet */

  return false;
}

template <class TBVH, class TAPI>
bool tledDeformableDeformableBVHTraverserImplCPU<TBVH, TAPI>::DoNarrowPhaseNodeFacetIterativeProjection(float *p_xi, const int nodeInd, const int primInd) {
  const int *facet = this->GetMasterMesh().GetFacet(primInd).NodeIndices;
  const int numFacetVtcs = SlaveMesh::Facet::NumberOfVertices;

  float facetNodes[numFacetVtcs][3], facetNodeNormals[numFacetVtcs][3];
  float n[3];

  for (int v = 0; v < numFacetVtcs; v++) {
    const float *nx = this->GetSlaveMesh().GetNodeCoordinates(facet[v]), *nn = this->GetSlaveMesh().GetNodeNormalCached(facet[v]);

    std::copy(nx, nx + 3, facetNodes[v]);
    std::copy(nn, nn + 3, facetNodeNormals[v]);
  }

  this->ProjectOntoFacetC0Iterative(p_xi, n, this->GetSlaveMesh().GetNodeCoordinates(nodeInd), facetNodes, facetNodeNormals);
  if (this->IsBetterNodeProjection(p_xi, this->GetNodeProjectionBuffer()[nodeInd].CollisionCoords)) {
    this->GetNodeProjectionBuffer()[nodeInd].FacetIndex = primInd;
    std::copy(p_xi, p_xi + 3, this->GetNodeProjectionBuffer()[nodeInd].CollisionCoords);
    std::copy(n, n + 3, this->GetNodeProjectionBuffer()[nodeInd].Normal);
	
    return true;
  }    

  return false;
}

template <class TBVH, class TAPI>
bool tledDeformableDeformableBVHTraverserImplCPU<TBVH, TAPI>::DoNarrowPhaseEdgeDetection(const int slaveEdgeInd, const int masterEdgeInd) {
  const std::pair<int, int> &sEdge = this->GetSlaveMesh().GetAllEdges()[slaveEdgeInd];
  const std::pair<int, int> &mEdge = this->GetMasterMesh().GetAllEdges()[masterEdgeInd];
  const float *A1 = this->GetSlaveMesh().GetNodeCoordinates(sEdge.first);
  const float *B1 = this->GetSlaveMesh().GetNodeCoordinates(sEdge.second);
  const float *C1 = this->GetMasterMesh().GetNodeCoordinates(mEdge.first);
  const float *D1 = this->GetMasterMesh().GetNodeCoordinates(mEdge.second);      
  const float *A0 = this->GetSlaveMesh().GetOldNodeCoordinates(sEdge.first);
  const float *B0 = this->GetSlaveMesh().GetOldNodeCoordinates(sEdge.second);
  const float *C0 = this->GetMasterMesh().GetOldNodeCoordinates(mEdge.first);
  const float *D0 = this->GetMasterMesh().GetOldNodeCoordinates(mEdge.second);

  float xi[3], masterN[3];
  
  assert(mEdge.first != sEdge.first && mEdge.first != sEdge.second && mEdge.second != sEdge.first && mEdge.second != sEdge.second);
  if (this->ComputeEdgeEdgeClosestPointParameters(xi[1], xi[2], A1, B1, C1, D1)) {   
    if (this->template ComputeEdgeEdgePenetrationDepth<true, true> (xi[0], masterN, A0, B0, A1, B1, this->GetSlaveMesh().GetNodeNormalCached(sEdge.first), this->GetSlaveMesh().GetNodeNormalCached(sEdge.second), xi[1], C0, D0, C1, D1, this->GetSlaveMesh().GetNodeNormalCached(mEdge.first), this->GetSlaveMesh().GetNodeNormalCached(mEdge.second), xi[2]) && this->IsBetterEdgeProjection(xi, this->GetEdgeProjectionBuffer()[slaveEdgeInd].CollisionCoords)) {
      std::copy(xi, xi + 3, this->GetEdgeProjectionBuffer()[slaveEdgeInd].CollisionCoords);
      std::copy(masterN, masterN + 3, this->GetEdgeProjectionBuffer()[slaveEdgeInd].Normal);
      this->GetEdgeProjectionBuffer()[slaveEdgeInd].EdgeIndex = masterEdgeInd;
	
      return true;
    }
  }

  return false;
}

template <class TBVH, class TAPI>
void tledDeformableDeformableBVHTraverserImplCPU<TBVH, TAPI>::RunBroadPhase() {
  this->GetContactNodeBuffer().clear();
  this->GetContactEdgeBuffer().clear();

  for (std::vector<int>::const_iterator ic_candInd = this->GetCollisionCandidateBVs().begin(); ic_candInd < this->GetCollisionCandidateBVs().end(); ic_candInd++) {
    const typename MasterBVH::BoundingVolume &collParent = this->GetMasterBVH().GetBV(*ic_candInd);
    
    assert(collParent.PrimitiveIndex == -1);      
    if (MasterBVH::BVHOrder == 2) {
      if (this->DoMaster()) this->DoBroadPhaseDetectionRecursive(collParent.ChildIndices[0], collParent.ChildIndices[1]);
      else this->DoBroadPhaseDetectionRecursive(collParent.ChildIndices[1], collParent.ChildIndices[0]);
    } else {
      for (int c0 = 0; c0 < MasterBVH::BVHOrder - 1; c0++) if (collParent.ChildIndices[c0] >= 0) {
	  for (int c1 = c0 + 1; c1 < MasterBVH::BVHOrder; c1++) if (collParent.ChildIndices[c1] >= 0) {
	      if (this->DoMaster()) this->DoBroadPhaseDetectionRecursive(collParent.ChildIndices[c0], collParent.ChildIndices[c1]);
	      else this->DoBroadPhaseDetectionRecursive(collParent.ChildIndices[c1], collParent.ChildIndices[c0]);
	    }
	}
    }
  } /* for candidate BVs */  
  this->FinishBroadPhase();
}

template <class TBVH, class TAPI>
void tledDeformableDeformableBVHTraverserImplCPU<TBVH, TAPI>::ProcessNodeFacetNarrowPhaseItems(const std::vector<std::pair<int, int> >::const_iterator itemsBegin, const std::vector<std::pair<int, int> >::const_iterator itemsEnd) {
  int lastAddedSlave = -1;

  assert(this->GetContactNodeBuffer().size() == 0);
  for (std::vector<std::pair<int, int> >::const_iterator ic_i = itemsBegin; ic_i < itemsEnd; ic_i++) {    
    assert(ic_i + 1 == itemsEnd || ic_i->first <= (ic_i + 1)->first);
    if (this->DoNarrowPhaseNodeFacetDetection(ic_i->first, ic_i->second) && lastAddedSlave != ic_i->first) {
      this->GetContactNodeBuffer().push_back(ic_i->first);
      lastAddedSlave = ic_i->first;
    }    
    assert(lastAddedSlave <= ic_i->first);
    assert(lastAddedSlave == ic_i->first || this->GetNodeProjection(ic_i->first).FacetIndex < 0);
    assert(lastAddedSlave != ic_i->first || this->GetNodeProjection(ic_i->first).FacetIndex >= 0);
  } /* for broad-phase contact nodes */
  assert(this->GetContactNodeBuffer().size() == tledHelper::MakeSortedUnique(this->GetContactNodeBuffer()).size());  
}

template <class TBVH, class TAPI>
void tledDeformableDeformableBVHTraverserImplCPU<TBVH, TAPI>::ProcessEdgeEdgeNarrowPhaseItems(const std::vector<std::pair<int, int> >::const_iterator itemsBegin, const std::vector<std::pair<int, int> >::const_iterator itemsEnd) {
  int lastAddedSlave = -1;

  assert(this->GetContactEdgeBuffer().size() == 0);
  for (std::vector<std::pair<int, int> >::const_iterator ic_i = itemsBegin; ic_i < itemsEnd; ic_i++) {
    assert(ic_i + 1 == itemsEnd || ic_i->first <= (ic_i + 1)->first);
    if (this->DoNarrowPhaseEdgeDetection(ic_i->first, ic_i->second) && lastAddedSlave != ic_i->first) {
      this->GetContactEdgeBuffer().push_back(ic_i->first);
      lastAddedSlave = ic_i->first;
    }	
    assert(lastAddedSlave == ic_i->first || this->GetEdgeProjection(ic_i->first).EdgeIndex < 0);
    assert(lastAddedSlave != ic_i->first || this->GetEdgeProjection(ic_i->first).EdgeIndex >= 0);
  }
  assert(this->GetContactEdgeBuffer().size() == tledHelper::MakeSortedUnique(this->GetContactEdgeBuffer()).size());
}

