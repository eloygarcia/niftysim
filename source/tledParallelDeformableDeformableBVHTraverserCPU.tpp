// =========================================================================
// File:       tledParallelDeformableDeformableBVHTraverserCPU.tpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    January 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

template <class TBaseTraverser>
class tledParallelDeformableDeformableBVHTraverserCPU<TBaseTraverser>::_NodeFacetNarrowPhaseSubTraverser : public Superclass::NodeFacetNarrowPhaseSubTraverser {
protected:
  virtual bool DoNarrowPhaseNodeFacetDetection(const int nodeInd, const int primInd);
  virtual bool DoNarrowPhaseEdgeDetection(const int slaveEdgeInd, const int masterEdgeInd) { std::abort(); return false; }

public:
  _NodeFacetNarrowPhaseSubTraverser(std::vector<int> &r_slaves, std::vector<tledBVHTraverserCPU::MasterFacetContactData> &r_masters, SlaveBVH &r_slaveBVH, const MasterBVH &masterBVH) : NodeFacetNarrowPhaseSubTraverser(r_slaves, r_masters, r_slaveBVH, masterBVH) {}
  virtual ~_NodeFacetNarrowPhaseSubTraverser(void) {}
};

template <class TBaseTraverser>
class tledParallelDeformableDeformableBVHTraverserCPU<TBaseTraverser>::_EdgeEdgeNarrowPhaseSubTraverser : public Superclass::EdgeEdgeNarrowPhaseSubTraverser {
protected:
  virtual bool DoNarrowPhaseNodeFacetDetection(const int nodeInd, const int primInd) { std::abort(); return false; }
  virtual bool DoNarrowPhaseEdgeDetection(const int slaveEdgeInd, const int masterEdgeInd);

public:
  _EdgeEdgeNarrowPhaseSubTraverser(std::vector<int> &r_slaves, std::vector<tledBVHTraverserCPU::MasterEdgeContactData> &r_masters, SlaveBVH &r_slaveBVH, const MasterBVH &masterBVH) : EdgeEdgeNarrowPhaseSubTraverser(r_slaves, r_masters, r_slaveBVH, masterBVH) {}
  virtual ~_EdgeEdgeNarrowPhaseSubTraverser(void) {}
};

template <class TBaseTraverser>
void tledParallelDeformableDeformableBVHTraverserCPU<TBaseTraverser>::RunBroadPhase() {
  const int baseNumThreads = this->GetNumberOfThreads();

  this->GetContactNodeBuffer().clear();
  this->GetContactEdgeBuffer().clear();
  for (int t = 0; t < (int)this->GetCollisionCandidateBVs().size();) {
    int blockSize = std::min((int)this->GetCollisionCandidateBVs().size() - t, this->GetNumberOfThreads());
    std::vector<typename Superclass::BroadPhaseSubTraverser*> vp_subs;

    assert(this->GetNumberOfBroadPhaseThreads() == 0);
    if (blockSize >= this->GetNumberOfThreads()/2) {
      while ((int)vp_subs.size() < blockSize) {
	const int candInd = this->GetCollisionCandidateBVs()[t+vp_subs.size()];
	const typename MasterBVH::BoundingVolume &collParent = this->GetMasterBVH().GetBV(candInd);
      
	if (this->DoMaster()) vp_subs.push_back(this->LaunchSubTraverser(collParent.ChildIndices[0], collParent.ChildIndices[1]));
	else vp_subs.push_back(this->LaunchSubTraverser(collParent.ChildIndices[1], collParent.ChildIndices[0]));
      }
    } else {
      const int oldNumThreads = this->GetNumberOfThreads();

      assert(blockSize < oldNumThreads/2 && blockSize > 0);
      this->SetNumberOfThreads(this->GetNumberOfThreads()/blockSize);
      assert(this->GetNumberOfThreads() > 0);
      for (std::vector<int>::const_iterator ic_candInd = this->GetCollisionCandidateBVs().begin() + t; ic_candInd < this->GetCollisionCandidateBVs().end(); ic_candInd++) {
	const typename MasterBVH::BoundingVolume &collParent = this->GetMasterBVH().GetBV(*ic_candInd);

	std::vector<typename Superclass::BroadPhaseSubTraverser*> vp_candSubs;
       
	assert(collParent.PrimitiveIndex == -1);          
	if (this->DoMaster())  vp_candSubs = this->DescendUntilFullyLoaded(collParent.ChildIndices[0], collParent.ChildIndices[1]);
	else vp_candSubs = this->DescendUntilFullyLoaded(collParent.ChildIndices[1], collParent.ChildIndices[0]);

	vp_subs.insert(vp_subs.end(), vp_candSubs.begin(), vp_candSubs.end());
      }
      this->SetNumberOfThreads(oldNumThreads);
    } /* Many candidates else */

    assert((int)vp_subs.size() <= this->GetNumberOfThreads() + (this->GetNumberOfThreads() - 1)*MasterBVH::BoundingVolume::NumberOfChildBVs*SlaveBVH::BoundingVolume::NumberOfChildBVs);

    this->CollectBroadPhaseSubResults(vp_subs);    
    t += blockSize;
  } /* for candidate blocks */

  this->SetNumberOfThreads(baseNumThreads);
  this->FinishBroadPhase();
}

template <class TBaseTraverser>
bool tledParallelDeformableDeformableBVHTraverserCPU<TBaseTraverser>::_EdgeEdgeNarrowPhaseSubTraverser::DoNarrowPhaseEdgeDetection(const int slaveEdgeInd, const int masterEdgeInd) {
  const std::pair<int, int> &sEdge = this->GetSlaveMesh().GetAllEdges()[slaveEdgeInd];
  const float *slave0T0 = this->GetSlaveMesh().GetOldNodeCoordinates(sEdge.first), *slave1T0 = this->GetSlaveMesh().GetOldNodeCoordinates(sEdge.second);
  const float *slave0T1 = this->GetSlaveMesh().GetNodeCoordinates(sEdge.first), *slave1T1 = this->GetSlaveMesh().GetNodeCoordinates(sEdge.second);
  const std::pair<int, int> &mEdge = this->GetMasterMesh().GetAllEdges()[masterEdgeInd];
  const float *master0T0 = this->GetMasterMesh().GetOldNodeCoordinates(mEdge.first),*master1T0 = this->GetMasterMesh().GetOldNodeCoordinates(mEdge.second);
  const float *master0T1 = this->GetMasterMesh().GetNodeCoordinates(mEdge.first), *master1T1 = this->GetMasterMesh().GetNodeCoordinates(mEdge.second);      

  float xi[3], n[3];

  if (this->ComputeEdgeEdgeClosestPointParameters(xi[1], xi[2], slave0T1, slave1T1, master0T1, master1T1)) {   
    const float *masterN0 = this->GetSlaveMesh().GetNodeNormal(mEdge.first), *masterN1 = this->GetSlaveMesh().GetNodeNormal(mEdge.second);
    const float *slaveN0 = this->GetSlaveMesh().GetNodeNormal(sEdge.first), *slaveN1 = this->GetSlaveMesh().GetNodeNormal(sEdge.second);

    if (this->template ComputeEdgeEdgePenetrationDepth<true, true> (xi[0], n, slave0T0, slave1T0, slave0T1, slave1T1, slaveN0, slaveN1, xi[1], master0T0, master1T0, master0T1, master1T1, masterN0, masterN1, xi[2]) && this->IsBetterEdgeProjection(xi, this->GetEdgeProjectionBuffer()[slaveEdgeInd].CollisionCoords)) {
      std::copy(xi, xi + 3, this->GetEdgeProjectionBuffer()[slaveEdgeInd].CollisionCoords);
      std::copy(n, n + 3, this->GetEdgeProjectionBuffer()[slaveEdgeInd].Normal);
      this->GetEdgeProjectionBuffer()[slaveEdgeInd].EdgeIndex = masterEdgeInd;	

      return true;
    }
  }    

  return false;
}

template <class TBaseTraverser>
bool tledParallelDeformableDeformableBVHTraverserCPU<TBaseTraverser>::_NodeFacetNarrowPhaseSubTraverser::DoNarrowPhaseNodeFacetDetection(const int nodeInd, const int primInd) {
  float xi[3], n[3];

  if (this->DoBehindFacetTest(this->GetSlaveMesh().GetOldNodeCoordinates(nodeInd), this->GetMasterMesh().GetNormalisedOldFacetNormal(primInd), this->GetMasterMesh().GetOldNodeCoordinates(this->GetMasterMesh().GetFacet(primInd).NodeIndices[0]))
      && this->DoNodeFacetInitialProjection(xi, this->GetSlaveMesh().GetNodeCoordinates(nodeInd), this->GetMasterMesh().GetFacetProjectionOperator(primInd), this->GetNodeProjectionBuffer()[nodeInd].CollisionCoords)) {
    const int *facet = this->GetMasterMesh().GetFacet(primInd).NodeIndices;

    float facetNodes[MasterMesh::Facet::NumberOfVertices][3], facetNodeNormals[MasterMesh::Facet::NumberOfVertices][3];
    
    for (int vInd = 0; vInd < MasterMesh::Facet::NumberOfVertices; vInd++) {
      const float *x = this->GetSlaveMesh().GetNodeCoordinates(facet[vInd]), *n = this->GetSlaveMesh().GetNodeNormalCached(facet[vInd]);
      
      std::copy(x, x + 3, facetNodes[vInd]);
      std::copy(n, n + 3, facetNodeNormals[vInd]);
    }

    this->ProjectOntoFacetC0Iterative(xi, n, this->GetSlaveMesh().GetNodeCoordinates(nodeInd), facetNodes, facetNodeNormals);
    if (this->IsBetterNodeProjection(xi, this->GetNodeProjectionBuffer()[nodeInd].CollisionCoords)) {
      this->GetNodeProjectionBuffer()[nodeInd].FacetIndex = primInd;
      std::copy(xi, xi + 3, this->GetNodeProjectionBuffer()[nodeInd].CollisionCoords);	
      std::copy(n, n + 3, this->GetNodeProjectionBuffer()[nodeInd].Normal);

      return true;
    }    
  }

  return false;
}

template <class TBaseTraverser>
void tledParallelDeformableDeformableBVHTraverserCPU<TBaseTraverser>::PreNodeFacetNarrowPhaseHook() {
  Superclass::PreNodeFacetNarrowPhaseHook();
  this->UpdateNodeFacetDeformableMasterGeometry();
}

template <class TBaseTraverser>
typename tledParallelDeformableDeformableBVHTraverserCPU<TBaseTraverser>::NodeFacetNarrowPhaseSubTraverser* tledParallelDeformableDeformableBVHTraverserCPU<TBaseTraverser>::LaunchNodeFacetNarrowPhaseWorker(const std::vector<std::pair<int, int> >::const_iterator itemsBegin, const std::vector<std::pair<int, int> >::const_iterator itemsEnd, const int threadInd) {
  NodeFacetNarrowPhaseSubTraverser *p_sub = new _NodeFacetNarrowPhaseSubTraverser(this->GetThreadContactNodeBuffer(threadInd), this->GetThreadNodeProjectionBuffer(threadInd), this->GetSlaveBVH(), this->GetMasterBVH());

  this->InitNarrowPhaseSubTraverser(*p_sub);
  p_sub->SetBoostThread(*new boost::thread(NodeFacetNarrowPhaseSubTraverser::ItemProcessorWorker, p_sub, itemsBegin, itemsEnd));

  return p_sub;
}

template <class TBaseTraverser>
void tledParallelDeformableDeformableBVHTraverserCPU<TBaseTraverser>::PreEdgeEdgeNarrowPhaseHook() {
  Superclass::PreEdgeEdgeNarrowPhaseHook();
  this->UpdateEdgeEdgeDeformableMasterSlaveGeometry();
}

template <class TBaseTraverser>
typename tledParallelDeformableDeformableBVHTraverserCPU<TBaseTraverser>::EdgeEdgeNarrowPhaseSubTraverser* tledParallelDeformableDeformableBVHTraverserCPU<TBaseTraverser>::LaunchEdgeEdgeNarrowPhaseWorker(const std::vector<std::pair<int, int> >::const_iterator itemsBegin, const std::vector<std::pair<int, int> >::const_iterator itemsEnd, const int threadInd) {
  EdgeEdgeNarrowPhaseSubTraverser *p_sub = new _EdgeEdgeNarrowPhaseSubTraverser(this->GetThreadContactEdgeBuffer(threadInd), this->GetThreadEdgeProjectionBuffer(threadInd), this->GetSlaveBVH(), this->GetMasterBVH());

  this->InitNarrowPhaseSubTraverser(*p_sub);
  p_sub->SetBoostThread(*new boost::thread(EdgeEdgeNarrowPhaseSubTraverser::ItemProcessorWorker, p_sub, itemsBegin, itemsEnd));

  return p_sub;
}
