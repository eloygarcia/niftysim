// =========================================================================
// File:       tledParallelSelfCollisionBVHTraverserCPU.tpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    April 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

template <class TBaseTraverser>
class tledParallelSelfCollisionBVHTraverserCPU<TBaseTraverser>::_BroadPhaseSubTraverser : public tledParallelBVHTraverserImplCPU<TBaseTraverser>::BroadPhaseSubTraverser {
public:
  typedef typename tledParallelBVHTraverserImplCPU<TBaseTraverser>::BroadPhaseSubTraverser Superclass;

protected:
  virtual void AddNarrowPhaseTests(const int sFacetInd, const int mFacetInd) { 
    typedef typename SlaveMesh::Facet __SlaveFacet;
    typedef typename MasterMesh::Facet __MasterFacet;

    const __SlaveFacet &sFacet = this->GetSlaveMesh().GetFacet(sFacetInd);
    const __MasterFacet &mFacet = this->GetMasterMesh().GetFacet(mFacetInd);

    for (int const *pc_sn = sFacet.NodeIndices; pc_sn < sFacet.NodeIndices + __SlaveFacet::NumberOfVertices; pc_sn++) {
      for (int const *pc_mn = mFacet.NodeIndices; pc_mn < mFacet.NodeIndices + __MasterFacet::NumberOfVertices; pc_mn++) {
	if (*pc_mn == *pc_sn) return;
      }
    }
    
    Superclass::AddNarrowPhaseTests(sFacetInd, mFacetInd);
  }

public:
  _BroadPhaseSubTraverser(std::vector<std::pair<int, int> > &r_nodeFacetPairs, std::vector<std::pair<int, int> > &r_edgeEdgePairs, SlaveBVH &r_slaveBVH, const MasterBVH &masterBVH) : Superclass(r_nodeFacetPairs, r_edgeEdgePairs, r_slaveBVH, masterBVH) {}
  virtual ~_BroadPhaseSubTraverser(void) {}
};

template <class TBaseTraverser>
typename tledParallelBVHTraverserImplCPU<TBaseTraverser>::BroadPhaseSubTraverser* tledParallelSelfCollisionBVHTraverserCPU<TBaseTraverser>::LaunchSubTraverser(const int slaveStartInd, const int masterStartInd) {
  _BroadPhaseSubTraverser *p_sub = new _BroadPhaseSubTraverser(this->GetThreadNodeFacetNarrowPhasePairBuffer(this->GetNumberOfBroadPhaseThreads()), this->GetThreadEdgeEdgeNarrowPhasePairBuffer(this->GetNumberOfBroadPhaseThreads()), this->GetSlaveBVH(), this->GetMasterBVH());
  boost::thread *p_thread = new boost::thread(tledParallelBVHTraverserImplCPU<TBaseTraverser>::BroadPhaseSubTraverser::BroadPhaseWorker, p_sub, slaveStartInd, masterStartInd);

  this->IncNumberOfBroadPhaseThreads();
  p_sub->SetBoostThread(*p_thread);

  return p_sub;
}
