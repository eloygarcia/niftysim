// =========================================================================
// File:       tledDeformableRigidBVHTraverserCPU.tpp
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

template <class TDeformableBVH, class TRigidBVH, class TAPI>
bool tledDeformableRigidBVHTraverserImplCPU<TDeformableBVH, TRigidBVH, TAPI>::DoNarrowPhaseNodeFacetDefMaster(const int defPrimIndex, const int rigNodeIndex) {
  const float *cNode = this->GetMasterMesh().GetNodeCoordinates(rigNodeIndex);
  const int *mFacet = this->GetSlaveMesh().GetFacet(defPrimIndex).NodeIndices;

  float xi[3], n[3];

  assert(rigNodeIndex >= 0 && rigNodeIndex < this->GetMasterMesh().GetNumberOfNodes());
  assert(defPrimIndex >= 0 && defPrimIndex < this->GetSlaveMesh().GetNumberOfFacets());

  if (this->DoBehindFacetTest(cNode, this->GetSlaveMesh().GetNormalisedOldFacetNormalCached(defPrimIndex), this->GetSlaveMesh().GetOldNodeCoordinates(mFacet[0]))
      && this->DoNodeFacetInitialProjection(xi, cNode, this->GetSlaveMesh().GetFacetProjectionOperatorCached(defPrimIndex), this->GetNodeProjectionBuffer()[rigNodeIndex].CollisionCoords)) {
    const float facetNodes[][3] = {{this->GetSlaveMesh().GetNodeCoordinates(mFacet[0])[0], this->GetSlaveMesh().GetNodeCoordinates(mFacet[0])[1], this->GetSlaveMesh().GetNodeCoordinates(mFacet[0])[2]},
				   {this->GetSlaveMesh().GetNodeCoordinates(mFacet[1])[0], this->GetSlaveMesh().GetNodeCoordinates(mFacet[1])[1], this->GetSlaveMesh().GetNodeCoordinates(mFacet[1])[2]},
				   {this->GetSlaveMesh().GetNodeCoordinates(mFacet[2])[0], this->GetSlaveMesh().GetNodeCoordinates(mFacet[2])[1], this->GetSlaveMesh().GetNodeCoordinates(mFacet[2])[2]}};
    const float facetNodeNormals[][3] = {{this->GetSlaveMesh().GetNodeNormalCached(mFacet[0])[0], this->GetSlaveMesh().GetNodeNormalCached(mFacet[0])[1], this->GetSlaveMesh().GetNodeNormalCached(mFacet[0])[2]},
					 {this->GetSlaveMesh().GetNodeNormalCached(mFacet[1])[0], this->GetSlaveMesh().GetNodeNormalCached(mFacet[1])[1], this->GetSlaveMesh().GetNodeNormalCached(mFacet[1])[2]},
					 {this->GetSlaveMesh().GetNodeNormalCached(mFacet[2])[0], this->GetSlaveMesh().GetNodeNormalCached(mFacet[2])[1], this->GetSlaveMesh().GetNodeNormalCached(mFacet[2])[2]}};
    
    this->ProjectOntoFacetC0Iterative(xi, n, cNode, facetNodes, facetNodeNormals);    
    if (this->IsBetterNodeProjection(xi, this->GetNodeProjectionBuffer()[rigNodeIndex].CollisionCoords)) {
      this->GetNodeProjectionBuffer()[rigNodeIndex].FacetIndex = defPrimIndex;
      std::copy(xi, xi + 3, this->GetNodeProjectionBuffer()[rigNodeIndex].CollisionCoords);
      std::copy(n, n + 3, this->GetNodeProjectionBuffer()[rigNodeIndex].Normal);
      
      return true;
    }
  } /* if is potential collision */

  return false;
}

template <class TDeformableBVH, class TRigidBVH, class TAPI>
bool tledDeformableRigidBVHTraverserImplCPU<TDeformableBVH, TRigidBVH, TAPI>::DoNarrowPhaseEdgeEdgeDefMaster(const int defEdgeIndex, const int rigEdgeIndex) {
  const std::pair<int, int> &sEdge = this->GetMasterMesh().GetAllEdges()[rigEdgeIndex];
  const std::pair<int, int> &mEdge = this->GetSlaveMesh().GetAllEdges()[defEdgeIndex];
  const float *slave0 = this->GetMasterMesh().GetNodeCoordinates(sEdge.first), *slave1 = this->GetMasterMesh().GetNodeCoordinates(sEdge.second);
  const float *master0T1 = this->GetSlaveMesh().GetNodeCoordinates(mEdge.first), *master1T1 = this->GetSlaveMesh().GetNodeCoordinates(mEdge.second);
  
  float xi[3], n[3];

  tledBVHTraverserCPU::MasterEdgeContactData &r_projection = this->GetEdgeProjectionBuffer()[rigEdgeIndex];

  assert(rigEdgeIndex >= 0 && rigEdgeIndex < (int)this->GetMasterMesh().GetAllEdges().size());
  assert(defEdgeIndex >= 0 && defEdgeIndex < (int)this->GetSlaveMesh().GetAllEdges().size());
  
  if (this->ComputeEdgeEdgeClosestPointParameters(xi[1], xi[2], slave0, slave1, master0T1, master1T1)) {
    const float *master0T0 = this->GetSlaveMesh().GetOldNodeCoordinates(mEdge.first), *master1T0 = this->GetSlaveMesh().GetOldNodeCoordinates(mEdge.second);
    const float *slaveN0 = this->GetMasterMesh().GetNodeNormal(sEdge.first), *slaveN1 = this->GetMasterMesh().GetNodeNormal(sEdge.second);
    const float *masterN0 = this->GetSlaveMesh().GetNodeNormalCached(mEdge.first), *masterN1 = this->GetSlaveMesh().GetNodeNormalCached(mEdge.second);

    if (this->template ComputeEdgeEdgePenetrationDepth<false, true> (xi[0], n, slave0, slave1, slave0, slave1, slaveN0, slaveN1, xi[1], master0T0, master1T0, master0T1, master1T1, masterN0, masterN1, xi[2]) && this->IsBetterEdgeProjection(xi, r_projection.CollisionCoords)) {
      std::copy(xi, xi + 3, r_projection.CollisionCoords);
      std::copy(n, n + 3, r_projection.Normal);
      r_projection.EdgeIndex = defEdgeIndex;

      return true;
    }
  }

  return false;
}

template <class TDeformableBVH, class TRigidBVH, class TAPI>
bool tledDeformableRigidBVHTraverserImplCPU<TDeformableBVH, TRigidBVH, TAPI>::DoNarrowPhaseNodeFacetRigMaster(const int defNodeInd, const int rigFacetInd) {
  const float *oldP = this->GetSlaveMesh().GetOldNodeCoordinates(defNodeInd);
  const float *currP = this->GetSlaveMesh().GetNodeCoordinates(defNodeInd);
  const int *mFacet = this->GetMasterMesh().GetFacet(rigFacetInd).NodeIndices;
  const float *mV0 = this->GetMasterMesh().GetNodeCoordinates(mFacet[0]);

  typename Superclass::MasterFacetContactData &r_projection = this->GetNodeProjectionBuffer()[defNodeInd];
  float xi[3], mNormal[3];
  
  assert(this->DoMaster());
  this->GetMasterMesh().ComputeNormalisedFacetNormal(mNormal, rigFacetInd);
  if (this->DoBehindFacetTest(oldP, mNormal, mV0) && this->DoNodeFacetInitialProjection(xi, currP, this->GetMasterMesh().GetFacetProjectionOperator(rigFacetInd), this->GetNodeProjectionBuffer()[defNodeInd].CollisionCoords)) {
    const float facetNodes[][3] = {{this->GetMasterMesh().GetNodeCoordinates(mFacet[0])[0], this->GetMasterMesh().GetNodeCoordinates(mFacet[0])[1], this->GetMasterMesh().GetNodeCoordinates(mFacet[0])[2]},
				   {this->GetMasterMesh().GetNodeCoordinates(mFacet[1])[0], this->GetMasterMesh().GetNodeCoordinates(mFacet[1])[1], this->GetMasterMesh().GetNodeCoordinates(mFacet[1])[2]},
				   {this->GetMasterMesh().GetNodeCoordinates(mFacet[2])[0], this->GetMasterMesh().GetNodeCoordinates(mFacet[2])[1], this->GetMasterMesh().GetNodeCoordinates(mFacet[2])[2]}};
    const float facetNodeNormals[][3] = {{this->GetMasterMesh().GetNodeNormal(mFacet[0])[0], this->GetMasterMesh().GetNodeNormal(mFacet[0])[1], this->GetMasterMesh().GetNodeNormal(mFacet[0])[2]},
					 {this->GetMasterMesh().GetNodeNormal(mFacet[1])[0], this->GetMasterMesh().GetNodeNormal(mFacet[1])[1], this->GetMasterMesh().GetNodeNormal(mFacet[1])[2]},
					 {this->GetMasterMesh().GetNodeNormal(mFacet[2])[0], this->GetMasterMesh().GetNodeNormal(mFacet[2])[1], this->GetMasterMesh().GetNodeNormal(mFacet[2])[2]}};
    
    this->ProjectOntoFacetC0Iterative(xi, mNormal, currP, facetNodes, facetNodeNormals);
    if (this->IsBetterNodeProjection(xi, this->GetNodeProjectionBuffer()[defNodeInd].CollisionCoords)) {
      r_projection.FacetIndex = rigFacetInd;
      std::copy(xi, xi + 3, r_projection.CollisionCoords);
      std::copy(mNormal, mNormal + 3, r_projection.Normal);
      
      return true;
    }
  }

  return false;
}

template <class TDeformableBVH, class TRigidBVH, class TAPI>
bool tledDeformableRigidBVHTraverserImplCPU<TDeformableBVH, TRigidBVH, TAPI>::DoNarrowPhaseEdgeEdgeRigMaster(const int defEdgeIndex, const int rigEdgeIndex) {
  const std::pair<int, int> &sEdge = this->GetSlaveMesh().GetAllEdges()[defEdgeIndex];
  const std::pair<int, int> &mEdge = this->GetMasterMesh().GetAllEdges()[rigEdgeIndex];
  const float *master0 = this->GetMasterMesh().GetNodeCoordinates(mEdge.first), *master1 = this->GetMasterMesh().GetNodeCoordinates(mEdge.second);
  const float *slave0T1 = this->GetSlaveMesh().GetNodeCoordinates(sEdge.first), *slave1T1 = this->GetSlaveMesh().GetNodeCoordinates(sEdge.second);

  float xi[3], n[3];
  
  assert(rigEdgeIndex >= 0 && rigEdgeIndex < (int)this->GetMasterMesh().GetAllEdges().size());
  assert(defEdgeIndex >= 0 && defEdgeIndex < (int)this->GetSlaveMesh().GetAllEdges().size());
  
  assert(this->DoMaster());
  if (this->ComputeEdgeEdgeClosestPointParameters(xi[1], xi[2], slave0T1, slave1T1, master0, master1)) {
    const float *slave0T0 = this->GetSlaveMesh().GetOldNodeCoordinates(sEdge.first), *slave1T0 = this->GetSlaveMesh().GetOldNodeCoordinates(sEdge.second);
    const float *masterN0 = this->GetMasterMesh().GetNodeNormal(mEdge.first), *masterN1 = this->GetMasterMesh().GetNodeNormal(mEdge.second);
    const float *slaveN0 = this->GetSlaveMesh().GetNodeNormalCached(sEdge.first), *slaveN1 = this->GetSlaveMesh().GetNodeNormalCached(sEdge.second);

    if (this->template ComputeEdgeEdgePenetrationDepth<true, false> (xi[0], n, slave0T0, slave1T0, slave0T1, slave1T1, slaveN0, slaveN1, xi[1], master0, master1, master0, master1, masterN0, masterN1, xi[2]) && this->IsBetterEdgeProjection(xi, this->GetEdgeProjectionBuffer()[defEdgeIndex].CollisionCoords)) {
      this->GetEdgeProjectionBuffer()[defEdgeIndex].EdgeIndex = rigEdgeIndex;
      std::copy(xi, xi + 3, this->GetEdgeProjectionBuffer()[defEdgeIndex].CollisionCoords);
      std::copy(n, n + 3, this->GetEdgeProjectionBuffer()[defEdgeIndex].Normal);
      assert(this->GetEdgeProjectionBuffer()[defEdgeIndex].CollisionCoords[2] >= 0 && this->GetEdgeProjectionBuffer()[defEdgeIndex].CollisionCoords[2] <= 1 && this->GetEdgeProjectionBuffer()[defEdgeIndex].CollisionCoords[1] >= 0 && this->GetEdgeProjectionBuffer()[defEdgeIndex].CollisionCoords[1] <= 1);

      return true;
    }
  }

  return false;
}

template <class TDeformableBVH, class TRigidBVH, class TAPI>
void tledDeformableRigidBVHTraverserImplCPU<TDeformableBVH, TRigidBVH, TAPI>::RunBroadPhase() {
  this->GetContactNodeBuffer().clear();
  this->GetContactEdgeBuffer().clear();

  this->DoBroadPhaseDetectionRecursive(0, 0);
  this->FinishBroadPhase();
}

template <class TDeformableBVH, class TRigidBVH, class TAPI>
void tledDeformableRigidBVHTraverserImplCPU<TDeformableBVH, TRigidBVH, TAPI>::ProcessNodeFacetNarrowPhaseItems(const std::vector<std::pair<int, int> >::const_iterator itemsBegin, const std::vector<std::pair<int, int> >::const_iterator itemsEnd) {
  using namespace tledHelper;

  int lastAddedSlave = -1;

  if (!this->DoMaster()) {
    for (std::vector<std::pair<int, int> >::const_iterator ic_i = itemsBegin; ic_i < itemsEnd; ic_i++) {    
      if (DoNarrowPhaseNodeFacetDefMaster(ic_i->second, ic_i->first) && lastAddedSlave != ic_i->first) {
	this->GetContactNodeBuffer().push_back(ic_i->first);
	lastAddedSlave = ic_i->first;
      }
      assert(ic_i->first == lastAddedSlave || this->GetNodeProjection(ic_i->first).FacetIndex < 0);
      assert(lastAddedSlave != ic_i->first || this->GetNodeProjection(ic_i->first).FacetIndex >= 0);
    }    
  } else {
    for (std::vector<std::pair<int, int> >::const_iterator ic_i = itemsBegin; ic_i < itemsEnd; ic_i++) {    
      if (DoNarrowPhaseNodeFacetRigMaster(ic_i->first, ic_i->second) && lastAddedSlave != ic_i->first) {
	this->GetContactNodeBuffer().push_back(ic_i->first);
	lastAddedSlave = ic_i->first;
      }
      assert(lastAddedSlave == ic_i->first || this->GetNodeProjection(ic_i->first).FacetIndex < 0);
      assert(lastAddedSlave != ic_i->first || this->GetNodeProjection(ic_i->first).FacetIndex >= 0);
    }
  }
}

template <class TDeformableBVH, class TRigidBVH, class TAPI>
void tledDeformableRigidBVHTraverserImplCPU<TDeformableBVH, TRigidBVH, TAPI>::ProcessEdgeEdgeNarrowPhaseItems(const std::vector<std::pair<int, int> >::const_iterator itemsBegin, const std::vector<std::pair<int, int> >::const_iterator itemsEnd) {
  using namespace tledHelper;

  int lastAddedSlave = -1;

  if (!this->DoMaster()) {
    for (std::vector<std::pair<int, int> >::const_iterator ic_i = itemsBegin; ic_i < itemsEnd; ic_i++) {
      if (DoNarrowPhaseEdgeEdgeDefMaster(ic_i->second, ic_i->first) && lastAddedSlave != ic_i->first) {
	lastAddedSlave = ic_i->first;
	this->GetContactEdgeBuffer().push_back(ic_i->first);
      }
      assert(lastAddedSlave == ic_i->first || this->GetEdgeProjection(ic_i->first).EdgeIndex < 0);
      assert(lastAddedSlave != ic_i->first || this->GetEdgeProjection(ic_i->first).EdgeIndex >= 0);
    }      
  } else {
    for (std::vector<std::pair<int, int> >::const_iterator ic_i = itemsBegin; ic_i < itemsEnd; ic_i++) {
      if (DoNarrowPhaseEdgeEdgeRigMaster(ic_i->first, ic_i->second) && lastAddedSlave != ic_i->first) {
	this->GetContactEdgeBuffer().push_back(ic_i->first);
	lastAddedSlave = ic_i->first;
      }
      assert(lastAddedSlave == ic_i->first || this->GetEdgeProjection(ic_i->first).EdgeIndex < 0);
      assert(lastAddedSlave != ic_i->first || this->GetEdgeProjection(ic_i->first).EdgeIndex >= 0);
    }      
  } /* if master else ... */  
}
