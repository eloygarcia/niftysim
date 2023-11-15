// =========================================================================
// File:       tledParallelDeformableRigidBVHTraverserCPU.tpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    February 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

template <class TBaseTraverser>
class tledParallelDeformableRigidBVHTraverserCPU<TBaseTraverser>::BroadPhaseSubTraverser : public tledParallelBVHTraverserImplCPU<TBaseTraverser>::BroadPhaseSubTraverser {
public:
  typedef typename tledParallelBVHTraverserImplCPU<TBaseTraverser>::BroadPhaseSubTraverser Superclass;

protected:
  virtual void AddNarrowPhaseTests(const int sFacetInd, const int mFacetInd) { this->AddNarrowPhaseTestsSwitching(sFacetInd, mFacetInd); }

public:
  BroadPhaseSubTraverser(std::vector<std::pair<int, int> > &r_nodeFacetPairs, std::vector<std::pair<int, int> > &r_edgeEdgePairs, SlaveBVH &r_slaveBVH, const MasterBVH &masterBVH, const bool doMaster) : Superclass(r_nodeFacetPairs, r_edgeEdgePairs, r_slaveBVH, masterBVH) { 
    this->tledBVHTraverser::SetDoMaster(doMaster);
  }
  virtual ~BroadPhaseSubTraverser(void) {}
};

template <class TBaseTraverser>
class tledParallelDeformableRigidBVHTraverserCPU<TBaseTraverser>::_NodeFacetNarrowPhaseDefMasterSubTraverser : public NodeFacetNarrowPhaseSubTraverser {
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
class tledParallelDeformableRigidBVHTraverserCPU<TBaseTraverser>::_NodeFacetNarrowPhaseRigMasterSubTraverser : public NodeFacetNarrowPhaseSubTraverser {
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
class tledParallelDeformableRigidBVHTraverserCPU<TBaseTraverser>::_EdgeEdgeNarrowPhaseDefMasterSubTraverser : public EdgeEdgeNarrowPhaseSubTraverser {
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
class tledParallelDeformableRigidBVHTraverserCPU<TBaseTraverser>::_EdgeEdgeNarrowPhaseRigMasterSubTraverser : public EdgeEdgeNarrowPhaseSubTraverser {
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
typename tledParallelBVHTraverserImplCPU<TBaseTraverser>::BroadPhaseSubTraverser* tledParallelDeformableRigidBVHTraverserCPU<TBaseTraverser>::LaunchSubTraverser(const int slaveStartInd, const int masterStartInd) {
  BroadPhaseSubTraverser *p_sub = new BroadPhaseSubTraverser(this->GetThreadNodeFacetNarrowPhasePairBuffer(this->GetNumberOfBroadPhaseThreads()), this->GetThreadEdgeEdgeNarrowPhasePairBuffer(this->GetNumberOfBroadPhaseThreads()),
							     this->GetSlaveBVH(), this->GetMasterBVH(), this->DoMaster());
  boost::thread *p_thread = new boost::thread(tledParallelBVHTraverserImplCPU<TBaseTraverser>::BroadPhaseSubTraverser::BroadPhaseWorker, p_sub, slaveStartInd, masterStartInd);

#ifndef NDEBUG
  std::cout << "Launching sub-thread " << this->GetNumberOfBroadPhaseThreads() << std::endl;
#endif

  this->IncNumberOfBroadPhaseThreads();
  p_sub->SetBoostThread(*p_thread);

  return p_sub;
}

template <class TBaseTraverser>
void tledParallelDeformableRigidBVHTraverserCPU<TBaseTraverser>::RunBroadPhase() {
  std::vector<typename Superclass::BroadPhaseSubTraverser*> vp_subs;

  this->GetContactNodeBuffer().clear();
  this->GetContactEdgeBuffer().clear();

  vp_subs = this->DescendUntilFullyLoaded(0, 0);
  this->CollectBroadPhaseSubResults(vp_subs);
  this->FinishBroadPhase();
}

template <class TBaseTraverser>
bool tledParallelDeformableRigidBVHTraverserCPU<TBaseTraverser>::_NodeFacetNarrowPhaseDefMasterSubTraverser::DoNarrowPhaseNodeFacetDefMaster(const int defPrimIndex, const int rigNodeIndex) {
  const float *cNode = this->GetMasterMesh().GetNodeCoordinates(rigNodeIndex);
  const int *mFacet = this->GetSlaveMesh().GetFacet(defPrimIndex).NodeIndices;

  float xi[3], n[3];

  assert(rigNodeIndex >= 0 && rigNodeIndex < this->GetMasterMesh().GetNumberOfNodes());
  assert(defPrimIndex >= 0 && defPrimIndex < this->GetSlaveMesh().GetNumberOfFacets());

  if (this->DoBehindFacetTest(cNode, this->GetSlaveMesh().GetNormalisedOldFacetNormal(defPrimIndex), this->GetSlaveMesh().GetOldNodeCoordinates(mFacet[0]))
      && this->DoNodeFacetInitialProjection(xi, cNode, this->GetSlaveMesh().GetFacetProjectionOperator(defPrimIndex), this->GetNodeProjectionBuffer()[rigNodeIndex].CollisionCoords)) {
    float facetNodes[SlaveMesh::Facet::NumberOfVertices][3], facetNodeNormals[SlaveMesh::Facet::NumberOfVertices][3];
    
    for (int vInd = 0; vInd < SlaveMesh::Facet::NumberOfVertices; vInd++) {
      float const *pc_x, *pc_n;
      
      assert(mFacet[vInd] >= 0 && mFacet[vInd] < this->GetSlaveMesh().GetNumberOfNodes());
      
      pc_x = this->GetSlaveMesh().GetNodeCoordinates(mFacet[vInd]);
      std::copy(pc_x, pc_x + 3, facetNodes[vInd]);
      
      pc_n = this->GetSlaveMesh().GetNodeNormalCached(mFacet[vInd]);
      std::copy(pc_n, pc_n + 3, facetNodeNormals[vInd]);
    }
    
    this->ProjectOntoFacetC0Iterative(xi, n, cNode, facetNodes, facetNodeNormals);    
    if (this->IsBetterNodeProjection(xi, this->GetNodeProjection(rigNodeIndex).CollisionCoords)) {
      this->GetNodeProjectionBuffer()[rigNodeIndex].FacetIndex = defPrimIndex;
      std::copy(xi, xi + 3, this->GetNodeProjectionBuffer()[rigNodeIndex].CollisionCoords);
      std::copy(n, n + 3, this->GetNodeProjectionBuffer()[rigNodeIndex].Normal);

      return true;
    }
  } /* if is potential collision */

  return false;
}

template <class TBaseTraverser>
bool tledParallelDeformableRigidBVHTraverserCPU<TBaseTraverser>::_NodeFacetNarrowPhaseRigMasterSubTraverser::DoNarrowPhaseNodeFacetRigMaster(const int defNodeInd, const int rigFacetInd) {
  const float *oldP = this->GetSlaveMesh().GetOldNodeCoordinates(defNodeInd);
  const float *currP = this->GetSlaveMesh().GetNodeCoordinates(defNodeInd);
  const int *mFacet = this->GetMasterMesh().GetFacet(rigFacetInd).NodeIndices;
  const float *mV0 = this->GetMasterMesh().GetNodeCoordinates(mFacet[0]);

  typename Superclass::MasterFacetContactData &r_projection = this->GetNodeProjectionBuffer()[defNodeInd];
  float xi[3], mNormal[3];
  
  assert(this->DoMaster());
  this->GetMasterMesh().ComputeNormalisedFacetNormal(mNormal, rigFacetInd);
  if (this->DoBehindFacetTest(oldP, mNormal, mV0) && this->DoNodeFacetInitialProjection(xi, currP, this->GetMasterMesh().GetFacetProjectionOperator(rigFacetInd), this->GetNodeProjectionBuffer()[defNodeInd].CollisionCoords)) {
    float facetNodes[MasterMesh::Facet::NumberOfVertices][3], facetNodeNormals[MasterMesh::Facet::NumberOfVertices][3];
    
    for (int vInd = 0; vInd < MasterMesh::Facet::NumberOfVertices; vInd++) {
      const float *x = this->GetMasterMesh().GetNodeCoordinates(mFacet[vInd]), *n = this->GetMasterMesh().GetNodeNormal(mFacet[vInd]);

      std::copy(x, x + 3, facetNodes[vInd]);
      std::copy(n, n + 3, facetNodeNormals[vInd]);
    }

    this->ProjectOntoFacetC0Iterative(xi, mNormal, currP, facetNodes, facetNodeNormals);
    if (this->IsBetterNodeProjection(xi, r_projection.CollisionCoords)) {
      r_projection.FacetIndex = rigFacetInd;
      std::copy(xi, xi + 3, r_projection.CollisionCoords);
      std::copy(mNormal, mNormal + 3, r_projection.Normal);

      return true;
    }
  }

  return false;
}

template <class TBaseTraverser>
bool tledParallelDeformableRigidBVHTraverserCPU<TBaseTraverser>::_EdgeEdgeNarrowPhaseDefMasterSubTraverser::DoNarrowPhaseEdgeEdgeDefMaster(const int defEdgeIndex, const int rigEdgeIndex) {
  const std::pair<int, int> &sEdge = this->GetMasterMesh().GetAllEdges()[rigEdgeIndex];
  const std::pair<int, int> &mEdge = this->GetSlaveMesh().GetAllEdges()[defEdgeIndex];
  const float *slave0 = this->GetMasterMesh().GetNodeCoordinates(sEdge.first), *slave1 = this->GetMasterMesh().GetNodeCoordinates(sEdge.second);
  const float *master0T1 = this->GetSlaveMesh().GetNodeCoordinates(mEdge.first), *master1T1 = this->GetSlaveMesh().GetNodeCoordinates(mEdge.second);

  float xi[3], n[3];
  tledBVHTraverserCPU::MasterEdgeContactData &r_projection = this->GetEdgeProjectionBuffer()[rigEdgeIndex];

  assert(rigEdgeIndex >= 0 && rigEdgeIndex < (int)this->GetMasterMesh().GetAllEdges().size());
  assert(defEdgeIndex >= 0 && defEdgeIndex < (int)this->GetSlaveMesh().GetAllEdges().size());
  assert(mEdge.first >= 0 && mEdge.first < this->GetSlaveMesh().GetNumberOfNodes());
  assert(mEdge.second >= 0 && mEdge.second < this->GetSlaveMesh().GetNumberOfNodes());      
  if (this->ComputeEdgeEdgeClosestPointParameters(xi[1], xi[2], slave0, slave1, master0T1, master1T1)) {
    const float *master0T0 = this->GetSlaveMesh().GetOldNodeCoordinates(mEdge.first), *master1T0 = this->GetSlaveMesh().GetOldNodeCoordinates(mEdge.second);
    const float *slaveN0 = this->GetMasterMesh().GetNodeNormal(sEdge.first), *slaveN1 = this->GetMasterMesh().GetNodeNormal(sEdge.second);
    const float *masterN0 = this->GetSlaveMesh().GetNodeNormal(mEdge.first), *masterN1 = this->GetSlaveMesh().GetNodeNormal(mEdge.second);

    if (this->template ComputeEdgeEdgePenetrationDepth<false, true> (xi[0], n, slave0, slave1, slave0, slave1, slaveN0, slaveN1, xi[1], master0T0, master1T0, master0T1, master1T1, masterN0, masterN1, xi[2]) && this->IsBetterEdgeProjection(xi, r_projection.CollisionCoords)) {
      std::copy(xi, xi + 3, r_projection.CollisionCoords);
      std::copy(n, n + 3, r_projection.Normal);
      r_projection.EdgeIndex = defEdgeIndex;

      return true;
    }
  }

  return false;
}

template <class TBaseTraverser>
bool tledParallelDeformableRigidBVHTraverserCPU<TBaseTraverser>::_EdgeEdgeNarrowPhaseRigMasterSubTraverser::DoNarrowPhaseEdgeEdgeRigMaster(const int defEdgeIndex, const int rigEdgeIndex) { 
  const std::pair<int, int> &sEdge = this->GetSlaveMesh().GetAllEdges()[defEdgeIndex];
  const std::pair<int, int> &mEdge = this->GetMasterMesh().GetAllEdges()[rigEdgeIndex];
  const float *slave0T1 = this->GetSlaveMesh().GetNodeCoordinates(sEdge.first), *slave1T1 = this->GetSlaveMesh().GetNodeCoordinates(sEdge.second);
  const float *master0 = this->GetMasterMesh().GetNodeCoordinates(mEdge.first), *master1 = this->GetMasterMesh().GetNodeCoordinates(mEdge.second);

  tledBVHTraverserCPU::MasterEdgeContactData &r_projection = this->GetEdgeProjectionBuffer()[defEdgeIndex];
  float xi[3], n[3];
  
  assert(rigEdgeIndex >= 0 && rigEdgeIndex < (int)this->GetMasterMesh().GetAllEdges().size());
  assert(defEdgeIndex >= 0 && defEdgeIndex < (int)this->GetSlaveMesh().GetAllEdges().size());
  
  assert(this->DoMaster());
  if (this->ComputeEdgeEdgeClosestPointParameters(xi[1], xi[2], slave0T1, slave1T1, master0, master1)) {
    const float *slave0T0 = this->GetSlaveMesh().GetOldNodeCoordinates(sEdge.first), *slave1T0 = this->GetSlaveMesh().GetOldNodeCoordinates(sEdge.second);
    const float *masterN0 = this->GetMasterMesh().GetNodeNormal(mEdge.first), *masterN1 = this->GetMasterMesh().GetNodeNormal(mEdge.second);
    const float *slaveN0 = this->GetSlaveMesh().GetNodeNormal(sEdge.first), *slaveN1 = this->GetSlaveMesh().GetNodeNormal(sEdge.second);

    if (this->template ComputeEdgeEdgePenetrationDepth<true, false> (xi[0], n, slave0T0, slave1T0, slave0T1, slave1T1, slaveN0, slaveN1, xi[1], master0, master1, master0, master1, masterN0, masterN1, xi[2]) && this->IsBetterEdgeProjection(xi, r_projection.CollisionCoords)) {
      r_projection.EdgeIndex = rigEdgeIndex;
      std::copy(xi, xi + 3, r_projection.CollisionCoords);
      std::copy(n, n + 3, r_projection.Normal);
      assert(r_projection.CollisionCoords[2] >= 0 && r_projection.CollisionCoords[2] <= 1 && r_projection.CollisionCoords[1] >= 0 && r_projection.CollisionCoords[1] <= 1);

      return true;
    }
  }

  return false;
}

template <class TBaseTraverser>
typename tledParallelDeformableRigidBVHTraverserCPU<TBaseTraverser>::NodeFacetNarrowPhaseSubTraverser* tledParallelDeformableRigidBVHTraverserCPU<TBaseTraverser>::LaunchNodeFacetNarrowPhaseWorker(const std::vector<std::pair<int, int> >::const_iterator itemsBegin, const std::vector<std::pair<int, int> >::const_iterator itemsEnd, const int threadInd) {
  NodeFacetNarrowPhaseSubTraverser *p_sub = this->DoMaster()? static_cast<NodeFacetNarrowPhaseSubTraverser*>(new _NodeFacetNarrowPhaseRigMasterSubTraverser(this->GetThreadContactNodeBuffer(threadInd), this->GetThreadNodeProjectionBuffer(threadInd), this->GetSlaveBVH(), this->GetMasterBVH()))
    : static_cast<NodeFacetNarrowPhaseSubTraverser*>(new _NodeFacetNarrowPhaseDefMasterSubTraverser(this->GetThreadContactNodeBuffer(threadInd), this->GetThreadNodeProjectionBuffer(threadInd), this->GetSlaveBVH(), this->GetMasterBVH()));

  this->InitNarrowPhaseSubTraverser(*p_sub);
  p_sub->SetBoostThread(*new boost::thread(NodeFacetNarrowPhaseSubTraverser::ItemProcessorWorker, p_sub, itemsBegin, itemsEnd));
 
  return p_sub;
}

template <class TBaseTraverser>
typename tledParallelDeformableRigidBVHTraverserCPU<TBaseTraverser>::EdgeEdgeNarrowPhaseSubTraverser* tledParallelDeformableRigidBVHTraverserCPU<TBaseTraverser>::LaunchEdgeEdgeNarrowPhaseWorker(const std::vector<std::pair<int, int> >::const_iterator itemsBegin, const std::vector<std::pair<int, int> >::const_iterator itemsEnd, const int threadInd) {
  EdgeEdgeNarrowPhaseSubTraverser *p_sub = this->DoMaster()? static_cast<EdgeEdgeNarrowPhaseSubTraverser*>(new _EdgeEdgeNarrowPhaseRigMasterSubTraverser(this->GetThreadContactEdgeBuffer(threadInd), this->GetThreadEdgeProjectionBuffer(threadInd), this->GetSlaveBVH(), this->GetMasterBVH()))
    : static_cast<EdgeEdgeNarrowPhaseSubTraverser*>(new _EdgeEdgeNarrowPhaseDefMasterSubTraverser(this->GetThreadContactEdgeBuffer(threadInd), this->GetThreadEdgeProjectionBuffer(threadInd), this->GetSlaveBVH(), this->GetMasterBVH()));

  this->InitNarrowPhaseSubTraverser(*p_sub);
  p_sub->SetBoostThread(*new boost::thread(EdgeEdgeNarrowPhaseSubTraverser::ItemProcessorWorker, p_sub, itemsBegin, itemsEnd));

  return p_sub;
}

template <class TBaseTraverser>
void tledParallelDeformableRigidBVHTraverserCPU<TBaseTraverser>::PreNodeFacetNarrowPhaseHook() {
  Superclass::PreNodeFacetNarrowPhaseHook();
  if (!this->DoMaster()) this->UpdateNodeFacetDeformableMasterGeometry();
}

template <class TBaseTraverser>
void tledParallelDeformableRigidBVHTraverserCPU<TBaseTraverser>::PreEdgeEdgeNarrowPhaseHook() {
  Superclass::PreEdgeEdgeNarrowPhaseHook();
  if (!this->DoMaster()) this->UpdateEdgeEdgeDeformableMasterGeometry();
  else this->UpdateEdgeEdgeDeformableSlaveGeometry();
}
