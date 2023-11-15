// =========================================================================
// File:       tledDeformableMovingRigidBVHTraverserCPU.tpp
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
bool tledDeformableMovingRigidBVHTraverserImplCPU<TDeformableBVH, TRigidBVH, TAPI>::DoNarrowPhaseNodeFacetDefMaster(const int defPrimIndex, const int rigNodeIndex) {
  const int *facet = this->GetSlaveMesh().GetFacet(defPrimIndex).NodeIndices;
  const float *oV0 = this->GetSlaveMesh().GetOldNodeCoordinates(facet[0]);
  const float *oldN = this->GetSlaveMesh().GetNormalisedOldFacetNormalCached(defPrimIndex);
  const float *currCNode = this->GetMasterMesh().GetNodeCoordinates(rigNodeIndex);
  const float *oldP = this->GetMasterMesh().GetOldNodeCoordinates(rigNodeIndex);
  
  assert(defPrimIndex >= 0 && defPrimIndex < this->GetSlaveMesh().GetNumberOfFacets());
  assert(rigNodeIndex >= 0 && rigNodeIndex < this->GetMasterMesh().GetNumberOfNodes());
  assert(this->GetNodeProjectionBuffer()[rigNodeIndex].FacetIndex != defPrimIndex);
  if (this->DoBehindFacetTest(oldP, oldN, oV0)) {
    float xi[3], n[3];

    /*
     * Cull situations where we cannot go back through facet, as useless to contact solver.
     */
    if (this->DoNodeFacetInitialProjection(xi, currCNode, this->GetSlaveMesh().GetFacetProjectionOperatorCached(defPrimIndex), this->GetNodeProjectionBuffer()[rigNodeIndex].CollisionCoords)) {
      const float facetNodes[][3] = {{this->GetSlaveMesh().GetNodeCoordinates(facet[0])[0], this->GetSlaveMesh().GetNodeCoordinates(facet[0])[1], this->GetSlaveMesh().GetNodeCoordinates(facet[0])[2]},
				     {this->GetSlaveMesh().GetNodeCoordinates(facet[1])[0], this->GetSlaveMesh().GetNodeCoordinates(facet[1])[1], this->GetSlaveMesh().GetNodeCoordinates(facet[1])[2]},
				     {this->GetSlaveMesh().GetNodeCoordinates(facet[2])[0], this->GetSlaveMesh().GetNodeCoordinates(facet[2])[1], this->GetSlaveMesh().GetNodeCoordinates(facet[2])[2]}};
      const float facetNodeNormals[][3] = {{this->GetSlaveMesh().GetNodeNormalCached(facet[0])[0], this->GetSlaveMesh().GetNodeNormalCached(facet[0])[1], this->GetSlaveMesh().GetNodeNormalCached(facet[0])[2]},
					   {this->GetSlaveMesh().GetNodeNormalCached(facet[1])[0], this->GetSlaveMesh().GetNodeNormalCached(facet[1])[1], this->GetSlaveMesh().GetNodeNormalCached(facet[1])[2]},
					   {this->GetSlaveMesh().GetNodeNormalCached(facet[2])[0], this->GetSlaveMesh().GetNodeNormalCached(facet[2])[1], this->GetSlaveMesh().GetNodeNormalCached(facet[2])[2]}};

      this->ProjectOntoFacetC0Iterative(xi, n, currCNode, facetNodes, facetNodeNormals);
      assert(xi[0] >= 0 && xi[1] >= 0 && xi[0] + xi[1] <= 1 + 1e-6);

      if (this->IsBetterNodeProjection(xi, this->GetNodeProjectionBuffer()[rigNodeIndex].CollisionCoords)) {
	this->GetNodeProjectionBuffer()[rigNodeIndex].FacetIndex = defPrimIndex;
	std::copy(xi, xi + 3, this->GetNodeProjectionBuffer()[rigNodeIndex].CollisionCoords);
	std::copy(n, n + 3, this->GetNodeProjectionBuffer()[rigNodeIndex].Normal);

	return true;
      }
    }
  } /* if not was behind facet */

  return false;
}

template <class TDeformableBVH, class TRigidBVH, class TAPI>
bool tledDeformableMovingRigidBVHTraverserImplCPU<TDeformableBVH, TRigidBVH, TAPI>::DoNarrowPhaseEdgeEdgeDefMaster(const int defEdgeIndex, const int rigEdgeIndex) {
  const std::pair<int, int> &sEdge = this->GetMasterMesh().GetAllEdges()[rigEdgeIndex];
  const std::pair<int, int> &mEdge = this->GetSlaveMesh().GetAllEdges()[defEdgeIndex];
  const float *master0T1 = this->GetSlaveMesh().GetNodeCoordinates(mEdge.first), *master1T1 = this->GetSlaveMesh().GetNodeCoordinates(mEdge.second);
  const float *slave0T1 = this->GetMasterMesh().GetNodeCoordinates(sEdge.first), *slave1T1 = this->GetMasterMesh().GetNodeCoordinates(sEdge.second);      

  float n[3], xi[3];

  assert(defEdgeIndex >= 0 && defEdgeIndex < this->GetSlaveMesh().GetNumberOfEdges());
  assert(rigEdgeIndex >= 0 && rigEdgeIndex < this->GetMasterMesh().GetNumberOfEdges());
  if (this->ComputeEdgeEdgeClosestPointParameters(xi[1], xi[2], slave0T1, slave1T1, master0T1, master1T1)) {
    const float *master0T0 = this->GetSlaveMesh().GetOldNodeCoordinates(mEdge.first), *master1T0 = this->GetSlaveMesh().GetOldNodeCoordinates(mEdge.second);
    const float *masterN0 = this->GetSlaveMesh().GetNodeNormalCached(mEdge.first), *masterN1 = this->GetSlaveMesh().GetNodeNormalCached(mEdge.second);
    const float *slave0T0 = this->GetMasterMesh().GetOldNodeCoordinates(sEdge.first), *slave1T0 = this->GetMasterMesh().GetOldNodeCoordinates(sEdge.second);
    const float *slaveN0 = this->GetMasterMesh().GetNodeNormal(sEdge.first), *slaveN1 = this->GetMasterMesh().GetNodeNormal(sEdge.second);
    
    if (this->template ComputeEdgeEdgePenetrationDepth<true, true> (xi[0], n, slave0T0, slave1T0, slave0T1, slave1T1, slaveN0, slaveN1, xi[1], master0T0, master1T0, master0T1, master1T1, masterN0, masterN1, xi[2]) && this->IsBetterEdgeProjection(xi, this->GetEdgeProjectionBuffer()[rigEdgeIndex].CollisionCoords)) {
      std::copy(xi, xi + 3, this->GetEdgeProjectionBuffer()[rigEdgeIndex].CollisionCoords);
      std::copy(n, n + 3, this->GetEdgeProjectionBuffer()[rigEdgeIndex].Normal);
      this->GetEdgeProjectionBuffer()[rigEdgeIndex].EdgeIndex = defEdgeIndex;
      
      return true;
    }
  }

  return false;
}

template <class TDeformableBVH, class TRigidBVH, class TAPI>
bool tledDeformableMovingRigidBVHTraverserImplCPU<TDeformableBVH, TRigidBVH, TAPI>::DoNarrowPhaseNodeFacetRigMaster(const int defNodeIndex, const int rigFacetIndex) {
  const int *facet = this->GetMasterMesh().GetFacet(rigFacetIndex).NodeIndices;
  const float *oV0 = this->GetMasterMesh().GetOldNodeCoordinates(facet[0]);
  const float *oldN = this->GetMasterMesh().GetOldFacetNormal(rigFacetIndex);
  const float *currCNode = this->GetSlaveMesh().GetNodeCoordinates(defNodeIndex);
  const float *oldP = this->GetSlaveMesh().GetOldNodeCoordinates(defNodeIndex);
  
  assert(rigFacetIndex >= 0 && rigFacetIndex < this->GetMasterMesh().GetNumberOfFacets());
  assert(defNodeIndex >= 0 && defNodeIndex < this->GetSlaveMesh().GetNumberOfNodes());
  assert(this->GetNodeProjectionBuffer()[defNodeIndex].FacetIndex != rigFacetIndex);
  if (this->DoBehindFacetTest(oldP, oldN, oV0)) {
    float xi[3], n[3];
    
    /*
     * Cull situations where we cannot go back through facet, as useless to contact solver.
     */
    if (this->DoNodeFacetInitialProjection(xi, currCNode, this->GetMasterMesh().GetFacetProjectionOperator(rigFacetIndex), this->GetNodeProjectionBuffer()[defNodeIndex].CollisionCoords)) {
      const float facetNodes[][3] = {{this->GetMasterMesh().GetNodeCoordinates(facet[0])[0], this->GetMasterMesh().GetNodeCoordinates(facet[0])[1], this->GetMasterMesh().GetNodeCoordinates(facet[0])[2]},
				     {this->GetMasterMesh().GetNodeCoordinates(facet[1])[0], this->GetMasterMesh().GetNodeCoordinates(facet[1])[1], this->GetMasterMesh().GetNodeCoordinates(facet[1])[2]},
				     {this->GetMasterMesh().GetNodeCoordinates(facet[2])[0], this->GetMasterMesh().GetNodeCoordinates(facet[2])[1], this->GetMasterMesh().GetNodeCoordinates(facet[2])[2]}};
      const float facetNodeNormals[][3] = {{this->GetMasterMesh().GetNodeNormal(facet[0])[0], this->GetMasterMesh().GetNodeNormal(facet[0])[1], this->GetMasterMesh().GetNodeNormal(facet[0])[2]},
					   {this->GetMasterMesh().GetNodeNormal(facet[1])[0], this->GetMasterMesh().GetNodeNormal(facet[1])[1], this->GetMasterMesh().GetNodeNormal(facet[1])[2]},
					   {this->GetMasterMesh().GetNodeNormal(facet[2])[0], this->GetMasterMesh().GetNodeNormal(facet[2])[1], this->GetMasterMesh().GetNodeNormal(facet[2])[2]}};

      this->ProjectOntoFacetC0Iterative(xi, n, currCNode, facetNodes, facetNodeNormals);
      if (this->IsBetterNodeProjection(xi, this->GetNodeProjectionBuffer()[defNodeIndex].CollisionCoords)) {
	this->GetNodeProjectionBuffer()[defNodeIndex].FacetIndex = rigFacetIndex;
	std::copy(xi, xi + 3, this->GetNodeProjectionBuffer()[defNodeIndex].CollisionCoords);
	std::copy(n, n + 3, this->GetNodeProjectionBuffer()[defNodeIndex].Normal);
	
	return true;
      }
    }
  } /* if not was behind facet */

  return false;  
}

template <class TDeformableBVH, class TRigidBVH, class TAPI>
bool tledDeformableMovingRigidBVHTraverserImplCPU<TDeformableBVH, TRigidBVH, TAPI>::DoNarrowPhaseEdgeEdgeRigMaster(const int defEdgeIndex, const int rigEdgeIndex) {
  const std::pair<int, int> &sEdge = this->GetSlaveMesh().GetAllEdges()[defEdgeIndex];
  const std::pair<int, int> &mEdge = this->GetMasterMesh().GetAllEdges()[rigEdgeIndex];  
  const float *slave0T1 = this->GetSlaveMesh().GetNodeCoordinates(sEdge.first), *slave1T1 = this->GetSlaveMesh().GetNodeCoordinates(sEdge.second);
  const float *master0T1 = this->GetMasterMesh().GetNodeCoordinates(mEdge.first), *master1T1 = this->GetMasterMesh().GetNodeCoordinates(mEdge.second);      

  float xi[3], n[3];

  assert(defEdgeIndex >= 0 && defEdgeIndex < this->GetSlaveMesh().GetNumberOfEdges());
  assert(rigEdgeIndex >= 0 && rigEdgeIndex < this->GetMasterMesh().GetNumberOfEdges());
  if (this->ComputeEdgeEdgeClosestPointParameters(xi[1], xi[2], slave0T1, slave1T1, master0T1, master1T1)) {
    const float *slave0T0 = this->GetSlaveMesh().GetOldNodeCoordinates(sEdge.first), *slave1T0 = this->GetSlaveMesh().GetOldNodeCoordinates(sEdge.second);
    const float *slaveN0 = this->GetSlaveMesh().GetNodeNormalCached(sEdge.first), *slaveN1 = this->GetSlaveMesh().GetNodeNormalCached(sEdge.second);    
    const float *master0T0 = this->GetMasterMesh().GetOldNodeCoordinates(mEdge.first), *master1T0 = this->GetMasterMesh().GetOldNodeCoordinates(mEdge.second);
    const float *masterN0 = this->GetMasterMesh().GetNodeNormal(mEdge.first), *masterN1 = this->GetMasterMesh().GetNodeNormal(mEdge.second);

    if (this->template ComputeEdgeEdgePenetrationDepth<true, true> (xi[0], n, slave0T0, slave1T0, slave0T1, slave1T1, slaveN0, slaveN1, xi[1], master0T0, master1T0, master0T1, master1T1, masterN0, masterN1, xi[2]) && this->IsBetterEdgeProjection(xi, this->GetEdgeProjectionBuffer()[defEdgeIndex].CollisionCoords)) {
      std::copy(xi, xi + 3, this->GetEdgeProjectionBuffer()[defEdgeIndex].CollisionCoords);
      std::copy(n, n + 3, this->GetEdgeProjectionBuffer()[defEdgeIndex].Normal);
      this->GetEdgeProjectionBuffer()[defEdgeIndex].EdgeIndex = rigEdgeIndex;

      return true;
    }
  }

  return false;
}
