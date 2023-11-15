// =========================================================================
// File:       tledParallelDeformableMovingRigidBVHTraverserCPU.tpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    March 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

template <class TBaseTraverser>
class tledParallelDeformableMovingRigidBVHTraverserCPU<TBaseTraverser>::_NodeFacetNarrowPhaseDefMasterSubTraverser : public NodeFacetNarrowPhaseSubTraverser {
protected:
  virtual bool DoNarrowPhaseNodeFacetDefMaster(const int defPrimIndex, const int rigNodeIndex);
  virtual bool DoNarrowPhaseEdgeEdgeDefMaster(const int defEdgeIndex, const int rigEdgeIndex) { std::abort(); return false; }
  virtual bool DoNarrowPhaseNodeFacetRigMaster(const int defNodeIndex, const int rigFacetIndex) { std::abort(); return false; }
  virtual bool DoNarrowPhaseEdgeEdgeRigMaster(const int defEdgeIndex, const int rigEdgeIndex) { std::abort(); return false; }
  
public:  
  _NodeFacetNarrowPhaseDefMasterSubTraverser(std::vector<int> &r_slaves, std::vector<tledBVHTraverserCPU::MasterFacetContactData> &r_masters, SlaveBVH &r_slaveBVH, const MasterBVH &masterBVH) : NodeFacetNarrowPhaseSubTraverser(r_slaves, r_masters, r_slaveBVH, masterBVH) {}
  virtual ~_NodeFacetNarrowPhaseDefMasterSubTraverser(void) {}
};

template <class TBaseTraverser>
class tledParallelDeformableMovingRigidBVHTraverserCPU<TBaseTraverser>::_NodeFacetNarrowPhaseRigMasterSubTraverser : public NodeFacetNarrowPhaseSubTraverser {
protected:
  virtual bool DoNarrowPhaseNodeFacetDefMaster(const int defPrimIndex, const int rigNodeIndex) { std::abort(); return false; }
  virtual bool DoNarrowPhaseEdgeEdgeDefMaster(const int defEdgeIndex, const int rigEdgeIndex) { std::abort(); return false; }
  virtual bool DoNarrowPhaseNodeFacetRigMaster(const int defNodeIndex, const int rigFacetIndex);
  virtual bool DoNarrowPhaseEdgeEdgeRigMaster(const int defEdgeIndex, const int rigEdgeIndex) { std::abort(); return false; }
  
public:  
  _NodeFacetNarrowPhaseRigMasterSubTraverser(std::vector<int> &r_slaves, std::vector<tledBVHTraverserCPU::MasterFacetContactData> &r_masters, SlaveBVH &r_slaveBVH, const MasterBVH &masterBVH) : NodeFacetNarrowPhaseSubTraverser(r_slaves, r_masters, r_slaveBVH, masterBVH) {}
  virtual ~_NodeFacetNarrowPhaseRigMasterSubTraverser(void) {}
};

template <class TBaseTraverser>
class tledParallelDeformableMovingRigidBVHTraverserCPU<TBaseTraverser>::_EdgeEdgeNarrowPhaseDefMasterSubTraverser : public EdgeEdgeNarrowPhaseSubTraverser {
protected:
  virtual bool DoNarrowPhaseNodeFacetDefMaster(const int defPrimIndex, const int rigNodeIndex) { std::abort(); return false; }
  virtual bool DoNarrowPhaseEdgeEdgeDefMaster(const int defEdgeIndex, const int rigEdgeIndex);
  virtual bool DoNarrowPhaseNodeFacetRigMaster(const int defNodeIndex, const int rigFacetIndex) { std::abort(); return false; }
  virtual bool DoNarrowPhaseEdgeEdgeRigMaster(const int defEdgeIndex, const int rigEdgeIndex) { std::abort(); return false; }

public:  
  _EdgeEdgeNarrowPhaseDefMasterSubTraverser(std::vector<int> &r_slaves, std::vector<tledBVHTraverserCPU::MasterEdgeContactData> &r_masters, SlaveBVH &r_slaveBVH, const MasterBVH &masterBVH) : EdgeEdgeNarrowPhaseSubTraverser(r_slaves, r_masters, r_slaveBVH, masterBVH) {}
  virtual ~_EdgeEdgeNarrowPhaseDefMasterSubTraverser(void) {}
};

template <class TBaseTraverser>
class tledParallelDeformableMovingRigidBVHTraverserCPU<TBaseTraverser>::_EdgeEdgeNarrowPhaseRigMasterSubTraverser : public EdgeEdgeNarrowPhaseSubTraverser {
protected:
  virtual bool DoNarrowPhaseNodeFacetDefMaster(const int defPrimIndex, const int rigNodeIndex) { std::abort(); return false; }
  virtual bool DoNarrowPhaseEdgeEdgeDefMaster(const int defEdgeIndex, const int rigEdgeIndex) { std::abort(); return false; }
  virtual bool DoNarrowPhaseNodeFacetRigMaster(const int defNodeIndex, const int rigFacetIndex) { std::abort(); return false; }
  virtual bool DoNarrowPhaseEdgeEdgeRigMaster(const int defEdgeIndex, const int rigEdgeIndex);

public:  
  _EdgeEdgeNarrowPhaseRigMasterSubTraverser(std::vector<int> &r_slaves, std::vector<tledBVHTraverserCPU::MasterEdgeContactData> &r_masters, SlaveBVH &r_slaveBVH, const MasterBVH &masterBVH) : EdgeEdgeNarrowPhaseSubTraverser(r_slaves, r_masters, r_slaveBVH, masterBVH) {}
  virtual ~_EdgeEdgeNarrowPhaseRigMasterSubTraverser(void) {}
};

template <class TBaseTraverser>
bool tledParallelDeformableMovingRigidBVHTraverserCPU<TBaseTraverser>::_NodeFacetNarrowPhaseDefMasterSubTraverser::DoNarrowPhaseNodeFacetDefMaster(const int defPrimIndex, const int rigNodeIndex) {
  const int *facet = this->GetSlaveMesh().GetFacet(defPrimIndex).NodeIndices;
  const float *oV0 = this->GetSlaveMesh().GetOldNodeCoordinates(facet[0]);
  const float *oldN = this->GetSlaveMesh().GetNormalisedOldFacetNormalCached(defPrimIndex);
  const float *currCNode = this->GetMasterMesh().GetNodeCoordinates(rigNodeIndex);
  const float *oldP = this->GetMasterMesh().GetOldNodeCoordinates(rigNodeIndex);
  
  typename Superclass::MasterFacetContactData &r_projection = this->GetNodeProjectionBuffer()[rigNodeIndex];
  
  assert(defPrimIndex >= 0 && defPrimIndex < this->GetSlaveMesh().GetNumberOfFacets());
  assert(rigNodeIndex >= 0 && rigNodeIndex < this->GetMasterMesh().GetNumberOfNodes());
  assert(r_projection.FacetIndex != defPrimIndex);
  if (this->DoBehindFacetTest(oldP, oldN, oV0)) {
    float xi[3], n[3];

    /*
     * Cull situations where we cannot go back through facet, as useless to contact solver.
     */
    if (this->DoNodeFacetInitialProjection(xi, currCNode, this->GetSlaveMesh().GetFacetProjectionOperator(defPrimIndex), r_projection.CollisionCoords)) {
      float facetNodes[SlaveMesh::Facet::NumberOfVertices][3], facetNodeNormals[SlaveMesh::Facet::NumberOfVertices][3];

      for (int vInd = 0; vInd < SlaveMesh::Facet::NumberOfVertices; vInd++) {
	float const *pc_x, *pc_n;

	pc_x = this->GetSlaveMesh().GetNodeCoordinates(facet[vInd]);
	std::copy(pc_x, pc_x + 3, facetNodes[vInd]);
	pc_n = this->GetSlaveMesh().GetNodeNormalCached(facet[vInd]);
	std::copy(pc_n, pc_n + 3, facetNodeNormals[vInd]);
      }

      this->ProjectOntoFacetC0Iterative(xi, n, currCNode, facetNodes, facetNodeNormals);
      assert(xi[0] >= 0 && xi[1] >= 0 && xi[0] + xi[1] <= 1 + 1e-6);

      if (this->IsBetterNodeProjection(xi, r_projection.CollisionCoords)) {
	r_projection.FacetIndex = defPrimIndex;
	std::copy(xi, xi + 3, r_projection.CollisionCoords);
	std::copy(n, n + 3, r_projection.Normal);

	return true;
      }
    }
  } /* if not was behind facet */

  return false;
}

template <class TBaseTraverser>
bool tledParallelDeformableMovingRigidBVHTraverserCPU<TBaseTraverser>::_NodeFacetNarrowPhaseRigMasterSubTraverser::DoNarrowPhaseNodeFacetRigMaster(const int defNodeIndex, const int rigFacetIndex) {
  const int *facet = this->GetMasterMesh().GetFacet(rigFacetIndex).NodeIndices;
  const float *oV0 = this->GetMasterMesh().GetOldNodeCoordinates(facet[0]);
  const float *oldN = this->GetMasterMesh().GetOldFacetNormal(rigFacetIndex);
  const float *currCNode = this->GetSlaveMesh().GetNodeCoordinates(defNodeIndex);
  const float *oldP = this->GetSlaveMesh().GetOldNodeCoordinates(defNodeIndex);
  
  typename Superclass::MasterFacetContactData &r_projection = this->GetNodeProjectionBuffer()[defNodeIndex];

  assert(rigFacetIndex >= 0 && rigFacetIndex < this->GetMasterMesh().GetNumberOfFacets());
  assert(defNodeIndex >= 0 && defNodeIndex < this->GetSlaveMesh().GetNumberOfNodes());
  assert(r_projection.FacetIndex != rigFacetIndex);
  if (this->DoBehindFacetTest(oldP, oldN, oV0)) {
    float xi[3], n[3];
    
    /*
     * Cull situations where we cannot go back through facet, as useless to contact solver.
     */
    if (this->DoNodeFacetInitialProjection(xi, currCNode, this->GetMasterMesh().GetFacetProjectionOperator(rigFacetIndex), r_projection.CollisionCoords)) {
      const float facetNodes[][3] = {{this->GetMasterMesh().GetNodeCoordinates(facet[0])[0], this->GetMasterMesh().GetNodeCoordinates(facet[0])[1], this->GetMasterMesh().GetNodeCoordinates(facet[0])[2]},
				     {this->GetMasterMesh().GetNodeCoordinates(facet[1])[0], this->GetMasterMesh().GetNodeCoordinates(facet[1])[1], this->GetMasterMesh().GetNodeCoordinates(facet[1])[2]},
				     {this->GetMasterMesh().GetNodeCoordinates(facet[2])[0], this->GetMasterMesh().GetNodeCoordinates(facet[2])[1], this->GetMasterMesh().GetNodeCoordinates(facet[2])[2]}};
      const float facetNodeNormals[][3] = {{this->GetMasterMesh().GetNodeNormal(facet[0])[0], this->GetMasterMesh().GetNodeNormal(facet[0])[1], this->GetMasterMesh().GetNodeNormal(facet[0])[2]},
					   {this->GetMasterMesh().GetNodeNormal(facet[1])[0], this->GetMasterMesh().GetNodeNormal(facet[1])[1], this->GetMasterMesh().GetNodeNormal(facet[1])[2]},
					   {this->GetMasterMesh().GetNodeNormal(facet[2])[0], this->GetMasterMesh().GetNodeNormal(facet[2])[1], this->GetMasterMesh().GetNodeNormal(facet[2])[2]}};

      this->ProjectOntoFacetC0Iterative(xi, n, currCNode, facetNodes, facetNodeNormals);
      if (this->IsBetterNodeProjection(xi, r_projection.CollisionCoords)) {
	r_projection.FacetIndex = rigFacetIndex;
	std::copy(xi, xi + 3, r_projection.CollisionCoords);
	std::copy(n, n + 3, r_projection.Normal);

	return true;
      }
    }
  } /* if not was behind facet */

  return false;
}

template <class TBaseTraverser>
bool tledParallelDeformableMovingRigidBVHTraverserCPU<TBaseTraverser>::_EdgeEdgeNarrowPhaseDefMasterSubTraverser::DoNarrowPhaseEdgeEdgeDefMaster(const int defEdgeIndex, const int rigEdgeIndex) {
  const std::pair<int, int> &sEdge = this->GetMasterMesh().GetAllEdges()[rigEdgeIndex];
  const std::pair<int, int> &mEdge = this->GetSlaveMesh().GetAllEdges()[defEdgeIndex];
  const float *master0T1 = this->GetSlaveMesh().GetNodeCoordinates(mEdge.first), *master1T1 = this->GetSlaveMesh().GetNodeCoordinates(mEdge.second);
  const float *slave0T1 = this->GetMasterMesh().GetNodeCoordinates(sEdge.first), *slave1T1 = this->GetMasterMesh().GetNodeCoordinates(sEdge.second);      

  tledBVHTraverserCPU::MasterEdgeContactData &r_projection = this->GetEdgeProjectionBuffer()[rigEdgeIndex];
  float n[3], xi[3];

  assert(defEdgeIndex >= 0 && defEdgeIndex < this->GetSlaveMesh().GetNumberOfEdges());
  assert(rigEdgeIndex >= 0 && rigEdgeIndex < this->GetMasterMesh().GetNumberOfEdges());
  assert(mEdge.first >= 0 && mEdge.first < this->GetSlaveMesh().GetNumberOfNodes());
  assert(mEdge.second >= 0 && mEdge.second < this->GetSlaveMesh().GetNumberOfNodes());      
  if (this->ComputeEdgeEdgeClosestPointParameters(xi[1], xi[2], slave0T1, slave1T1, master0T1, master1T1)) {
    const float *master0T0 = this->GetMasterMesh().GetOldNodeCoordinates(mEdge.first), *master1T0 = this->GetSlaveMesh().GetOldNodeCoordinates(mEdge.second);
    const float *masterN0 = this->GetSlaveMesh().GetNodeNormal(mEdge.first), *masterN1 = this->GetSlaveMesh().GetNodeNormal(mEdge.second);
    const float *slave0T0 = this->GetMasterMesh().GetOldNodeCoordinates(sEdge.first), *slave1T0 = this->GetMasterMesh().GetOldNodeCoordinates(sEdge.second);
    const float *slaveN0 = this->GetMasterMesh().GetNodeNormal(sEdge.first), *slaveN1 = this->GetMasterMesh().GetNodeNormal(sEdge.second);
    
    if (this->template ComputeEdgeEdgePenetrationDepth<true, true> (xi[0], n, slave0T0, slave1T0, slave0T1, slave1T1, slaveN0, slaveN1, xi[1], master0T0, master1T0, master0T1, master1T1, masterN0, masterN1, xi[2]) && this->IsBetterEdgeProjection(xi, r_projection.CollisionCoords)) {
      std::copy(xi, xi + 3, r_projection.CollisionCoords);
      std::copy(n, n + 3, r_projection.Normal);
      r_projection.EdgeIndex = defEdgeIndex;

      return true;
    }
  }
  
  return false;
}

template <class TBaseTraverser>
bool tledParallelDeformableMovingRigidBVHTraverserCPU<TBaseTraverser>::_EdgeEdgeNarrowPhaseRigMasterSubTraverser::DoNarrowPhaseEdgeEdgeRigMaster(const int defEdgeIndex, const int rigEdgeIndex) {
  const std::pair<int, int> &sEdge = this->GetSlaveMesh().GetAllEdges()[defEdgeIndex];
  const std::pair<int, int> &mEdge = this->GetMasterMesh().GetAllEdges()[rigEdgeIndex];
  const float *slave0T1 = this->GetSlaveMesh().GetNodeCoordinates(sEdge.first), *slave1T1 = this->GetSlaveMesh().GetNodeCoordinates(sEdge.second);
  const float *master0T1 = this->GetMasterMesh().GetNodeCoordinates(mEdge.first), *master1T1 = this->GetMasterMesh().GetNodeCoordinates(mEdge.second);      

  tledBVHTraverserCPU::MasterEdgeContactData &r_projection = this->GetEdgeProjectionBuffer()[defEdgeIndex];
  float xi[3], n[3];

  assert(defEdgeIndex >= 0 && defEdgeIndex < this->GetSlaveMesh().GetNumberOfEdges());
  assert(rigEdgeIndex >= 0 && rigEdgeIndex < this->GetMasterMesh().GetNumberOfEdges());
  if (this->ComputeEdgeEdgeClosestPointParameters(xi[1], xi[2], slave0T1, slave1T1, master0T1, master1T1)) {
    const float *slave0T0 = this->GetSlaveMesh().GetOldNodeCoordinates(sEdge.first), *slave1T0 = this->GetSlaveMesh().GetOldNodeCoordinates(sEdge.second);
    const float *slaveN0 = this->GetSlaveMesh().GetNodeNormal(sEdge.first), *slaveN1 = this->GetSlaveMesh().GetNodeNormal(sEdge.second);    
    const float *master0T0 = this->GetMasterMesh().GetOldNodeCoordinates(mEdge.first), *master1T0 = this->GetMasterMesh().GetOldNodeCoordinates(mEdge.second);
    const float *masterN0 = this->GetMasterMesh().GetNodeNormal(mEdge.first), *masterN1 = this->GetMasterMesh().GetNodeNormal(mEdge.second);

    if (this->template ComputeEdgeEdgePenetrationDepth<true, true> (xi[0], n, slave0T0, slave1T0, slave0T1, slave1T1, slaveN0, slaveN1, xi[1], master0T0, master1T0, master0T1, master1T1, masterN0, masterN1, xi[2]) && this->IsBetterEdgeProjection(xi, r_projection.CollisionCoords)) {
      std::copy(xi, xi + 3, r_projection.CollisionCoords);
      std::copy(n, n + 3, r_projection.Normal);
      r_projection.EdgeIndex = rigEdgeIndex;
      
      return true;
    }
  }

  return false;
}

template <class TBaseTraverser>
typename tledParallelDeformableMovingRigidBVHTraverserCPU<TBaseTraverser>::NodeFacetNarrowPhaseSubTraverser* tledParallelDeformableMovingRigidBVHTraverserCPU<TBaseTraverser>::LaunchNodeFacetNarrowPhaseWorker(const std::vector<std::pair<int, int> >::const_iterator itemsBegin, const std::vector<std::pair<int, int> >::const_iterator itemsEnd, const int threadInd) {
  NodeFacetNarrowPhaseSubTraverser *p_sub = this->DoMaster()? static_cast<NodeFacetNarrowPhaseSubTraverser*>(new _NodeFacetNarrowPhaseRigMasterSubTraverser(this->GetThreadContactNodeBuffer(threadInd), this->GetThreadNodeProjectionBuffer(threadInd), this->GetSlaveBVH(), this->GetMasterBVH()))
    : static_cast<NodeFacetNarrowPhaseSubTraverser*>(new _NodeFacetNarrowPhaseDefMasterSubTraverser(this->GetThreadContactNodeBuffer(threadInd), this->GetThreadNodeProjectionBuffer(threadInd), this->GetSlaveBVH(), this->GetMasterBVH()));

  this->InitNarrowPhaseSubTraverser(*p_sub);
  p_sub->SetBoostThread(*new boost::thread(NodeFacetNarrowPhaseSubTraverser::ItemProcessorWorker, p_sub, itemsBegin, itemsEnd));
 
  return p_sub;
}

template <class TBaseTraverser>
typename tledParallelDeformableMovingRigidBVHTraverserCPU<TBaseTraverser>::EdgeEdgeNarrowPhaseSubTraverser* tledParallelDeformableMovingRigidBVHTraverserCPU<TBaseTraverser>::LaunchEdgeEdgeNarrowPhaseWorker(const std::vector<std::pair<int, int> >::const_iterator itemsBegin, const std::vector<std::pair<int, int> >::const_iterator itemsEnd, const int threadInd) {
  EdgeEdgeNarrowPhaseSubTraverser *p_sub = this->DoMaster()? static_cast<EdgeEdgeNarrowPhaseSubTraverser*>(new _EdgeEdgeNarrowPhaseRigMasterSubTraverser(this->GetThreadContactEdgeBuffer(threadInd), this->GetThreadEdgeProjectionBuffer(threadInd), this->GetSlaveBVH(), this->GetMasterBVH()))
    : static_cast<EdgeEdgeNarrowPhaseSubTraverser*>(new _EdgeEdgeNarrowPhaseDefMasterSubTraverser(this->GetThreadContactEdgeBuffer(threadInd), this->GetThreadEdgeProjectionBuffer(threadInd), this->GetSlaveBVH(), this->GetMasterBVH()));

  this->InitNarrowPhaseSubTraverser(*p_sub);
  p_sub->SetBoostThread(*new boost::thread(EdgeEdgeNarrowPhaseSubTraverser::ItemProcessorWorker, p_sub, itemsBegin, itemsEnd));

  return p_sub;
}
