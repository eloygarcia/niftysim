// =========================================================================
// File:       tledParallelBVHTraverserCPU.tpp
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
class tledParallelBVHTraverserImplCPU<TBaseTraverser>::SubTraverserAPI : public tledBVHTraverserCPU {
};

template <class TBaseTraverser>
class tledParallelBVHTraverserImplCPU<TBaseTraverser>::BroadPhaseSubTraverser : public tledBVHTraverserImplCPU<MasterBVH, SlaveBVH, SubTraverserAPI>, public TraverserWorker {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledBVHTraverserImplCPU<MasterBVH, SlaveBVH, SubTraverserAPI> Superclass;
  /** @} */

  /**
   * \name Detection
   * @{
   */
private:
  std::vector<std::pair<int, int> > &mr_NodeFacetPairs, &mr_EdgeEdgePairs;

protected:
  virtual void RunBroadPhase(void) { std::abort(); }
  virtual void RunNarrowPhase(void) { std::abort(); }

  virtual void ProcessNodeFacetNarrowPhaseItems(const std::vector<std::pair<int, int> >::const_iterator itemsBegin, const std::vector<std::pair<int, int> >::const_iterator itemsEnd) { std::abort(); }
  virtual void ProcessEdgeEdgeNarrowPhaseItems(const std::vector<std::pair<int, int> >::const_iterator itemsBegin, const std::vector<std::pair<int, int> >::const_iterator itemsEnd) { std::abort(); }

public:
  virtual std::vector<std::pair<int, int> >& GetNodeFacetNarrowPhasePairs(void) { return mr_NodeFacetPairs; }
  virtual std::vector<std::pair<int, int> >& GetEdgeEdgeNarrowPhasePairs(void) { return mr_EdgeEdgePairs; }
  virtual const std::vector<std::pair<int, int> >& GetNodeFacetNarrowPhasePairs(void) const { return mr_NodeFacetPairs; }
  virtual const std::vector<std::pair<int, int> >& GetEdgeEdgeNarrowPhasePairs(void) const { return mr_EdgeEdgePairs; }

  virtual std::vector<int>& GetContactNodeBuffer(void) { std::abort(); return Superclass::GetContactNodeBuffer(); }
  virtual std::vector<int>& GetContactEdgeBuffer(void) { std::abort(); return Superclass::GetContactEdgeBuffer(); }

  virtual const std::vector<int>& GetSlaveNodeIndices(void) const { std::abort(); return Superclass::GetSlaveNodeIndices(); }
  virtual const std::vector<int>& GetSlaveEdgeIndices(void) const { std::abort(); return Superclass::GetSlaveEdgeIndices(); }

  static void BroadPhaseWorker(BroadPhaseSubTraverser *p_traverser, const int slaveStartInd, const int masterStartInd);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  BroadPhaseSubTraverser(std::vector<std::pair<int, int> > &r_nodeFacetBuffer, std::vector<std::pair<int, int> > &r_edgeEdgeBuffer, SlaveBVH &r_slaveBVH, const MasterBVH &masterBVH) : Superclass(r_slaveBVH, masterBVH), mr_NodeFacetPairs(r_nodeFacetBuffer), mr_EdgeEdgePairs(r_edgeEdgeBuffer) { 
    assert(r_nodeFacetBuffer.size() == 0);
    assert(r_edgeEdgeBuffer.size() == 0);
  }

  virtual ~BroadPhaseSubTraverser(void) {}
  /** @} */
};

template <class TBaseTraverser>
void tledParallelBVHTraverserImplCPU<TBaseTraverser>::BroadPhaseSubTraverser::BroadPhaseWorker(BroadPhaseSubTraverser *p_traverser, const int slaveStartInd, const int masterStartInd) {
  assert(p_traverser->GetNodeFacetNarrowPhasePairs().size() == 0 && p_traverser->GetEdgeEdgeNarrowPhasePairs().size() == 0);
  p_traverser->DoBroadPhaseDetectionRecursive(masterStartInd, slaveStartInd);
}

template <class TBaseTraverser>
void tledParallelBVHTraverserImplCPU<TBaseTraverser>::TraverserWorker::Join() { 
  if (mp_Thread != NULL) {
    mp_Thread->join(); 
    delete mp_Thread;
  }
}

template <class TBaseTraverser>
typename tledParallelBVHTraverserImplCPU<TBaseTraverser>::BroadPhaseSubTraverser* tledParallelBVHTraverserImplCPU<TBaseTraverser>::LaunchSubTraverser(const int slaveStartInd, const int masterStartInd) {
  BroadPhaseSubTraverser *p_sub = new BroadPhaseSubTraverser(this->GetThreadNodeFacetNarrowPhasePairBuffer(this->GetNumberOfBroadPhaseThreads()), this->GetThreadEdgeEdgeNarrowPhasePairBuffer(this->GetNumberOfBroadPhaseThreads()),
							     this->GetSlaveBVH(), this->GetMasterBVH());
  boost::thread *p_thread = new boost::thread(BroadPhaseSubTraverser::BroadPhaseWorker, p_sub, slaveStartInd, masterStartInd);

  tledLogDebugStream(tledHelper::Info() << "Launching sub-thread " << this->GetNumberOfBroadPhaseThreads());

  this->IncNumberOfBroadPhaseThreads();
  p_sub->SetBoostThread(*p_thread);

  return p_sub;
}

template <class TBaseTraverser>
std::vector<typename tledParallelBVHTraverserImplCPU<TBaseTraverser>::BroadPhaseSubTraverser*> tledParallelBVHTraverserImplCPU<TBaseTraverser>::DescendUntilFullyLoaded(const int slaveStartInd, const int masterStartInd) {
  std::vector<BroadPhaseSubTraverser*> vp_subs;
  std::deque<std::pair<int, int> > nextLevel;

  nextLevel.push_back(std::pair<int, int>(slaveStartInd, masterStartInd));
  while ((int)nextLevel.size() < this->GetNumberOfThreads() && nextLevel.size() > 0) {
    const std::pair<int, int> item = nextLevel.front();
    const typename MasterBVH::BoundingVolume &bv0 = this->GetMasterBVH().GetBV(item.second);
    const typename SlaveBVH::BoundingVolume &bv1 = this->GetSlaveBVH().GetBV(item.first);

    nextLevel.pop_front();
    
    if (bv0.PrimitiveIndex >= 0 || bv1.PrimitiveIndex >= 0) {
      this->DoBroadPhaseDetectionRecursive(masterStartInd, slaveStartInd);
    } else {
      for (int const *pc_masterChild = bv0.ChildIndices; pc_masterChild < bv0.ChildIndices + MasterBVH::BoundingVolume::NumberOfChildBVs; pc_masterChild++) for (int const *pc_slaveChild = bv1.ChildIndices; pc_slaveChild < bv1.ChildIndices + SlaveBVH::BoundingVolume::NumberOfChildBVs; pc_slaveChild++) {
	  const typename MasterBVH::BoundingVolume &child0 = this->GetMasterBVH().GetBV(*pc_masterChild);
	  const typename SlaveBVH::BoundingVolume &child1 = this->GetSlaveBVH().GetBV(*pc_slaveChild);
	  
	  if (child0.DoesIntersect(child1)) {
	    std::pair<int, int> newPair(*pc_slaveChild, *pc_masterChild);
	    
	    nextLevel.push_back(newPair);
	  }
	}
    }
  }
  
  tledLogDebugStream(tledHelper::Info() << "Finished collecting potential thread starting points, having found " << (int)nextLevel.size());

  while (nextLevel.size() > 0) {
    const std::pair<int, int> item = nextLevel.front();

    vp_subs.push_back(this->LaunchSubTraverser(item.first, item.second));
    nextLevel.pop_front();
    assert((int)vp_subs.size() <= this->GetNumberOfThreads() + MasterBVH::BoundingVolume::NumberOfChildBVs*SlaveBVH::BoundingVolume::NumberOfChildBVs - 1);
  }

  return vp_subs;
}

template <class TBaseTraverser>
void tledParallelBVHTraverserImplCPU<TBaseTraverser>::CollectBroadPhaseSubResults(std::vector<BroadPhaseSubTraverser*> &rvp_subTraversers) {
  for (typename std::vector<BroadPhaseSubTraverser*>::iterator ip_sub = rvp_subTraversers.begin(); ip_sub < rvp_subTraversers.end(); ip_sub++) {
    (*ip_sub)->Join();
    this->GetEdgeEdgeNarrowPhasePairs().insert(this->GetEdgeEdgeNarrowPhasePairs().end(), (*ip_sub)->GetEdgeEdgeNarrowPhasePairs().begin(), (*ip_sub)->GetEdgeEdgeNarrowPhasePairs().end());
    (*ip_sub)->GetEdgeEdgeNarrowPhasePairs().clear();
    this->GetNodeFacetNarrowPhasePairs().insert(this->GetNodeFacetNarrowPhasePairs().end(), (*ip_sub)->GetNodeFacetNarrowPhasePairs().begin(), (*ip_sub)->GetNodeFacetNarrowPhasePairs().end());
    (*ip_sub)->GetNodeFacetNarrowPhasePairs().clear();
    delete *ip_sub;
    *ip_sub = NULL;
  }  

#ifndef NDEBUG
  for (typename std::vector<BroadPhaseSubTraverser*>::iterator ip_sub = rvp_subTraversers.begin(); ip_sub < rvp_subTraversers.end(); ip_sub++) {
    assert(*ip_sub == NULL);
  }
  for (int t = 0; t < this->GetNumberOfBroadPhaseThreads(); t++) {
    assert(this->GetThreadNodeFacetNarrowPhasePairBuffer(t).size() == 0 && this->GetThreadEdgeEdgeNarrowPhasePairBuffer(t).size() == 0);
  }
#endif

  this->ResetBroadPhaseThreadCounter();  
}

template <class TBaseTraverser>
void tledParallelBVHTraverserImplCPU<TBaseTraverser>::MergeNodeFacetNarrowPhaseResults(NodeFacetNarrowPhaseSubTraverser *p_worker) {
  std::vector<int> &r_threadSlaveNodeList = p_worker->GetContactNodeBuffer();

  p_worker->Join();
  for (std::vector<int>::const_iterator ic_nInd = r_threadSlaveNodeList.begin(); ic_nInd < r_threadSlaveNodeList.end(); ic_nInd++) {
    typename Superclass::MasterFacetContactData &r_threadProjection = p_worker->GetNodeProjectionBuffer()[*ic_nInd];
    
    assert(r_threadProjection.FacetIndex >= 0);
    if (this->IsBetterNodeProjection(r_threadProjection.CollisionCoords, this->GetNodeProjectionBuffer()[*ic_nInd].CollisionCoords)) {
      assert(this->GetNodeProjectionBuffer()[*ic_nInd].FacetIndex < 0 || std::find(this->GetContactNodeBuffer().begin(), this->GetContactNodeBuffer().end(), *ic_nInd) != this->GetContactNodeBuffer().end());
      assert(this->GetNodeProjectionBuffer()[*ic_nInd].FacetIndex >= 0 || std::find(this->GetContactNodeBuffer().begin(), this->GetContactNodeBuffer().end(), *ic_nInd) == this->GetContactNodeBuffer().end());
      if (this->GetNodeProjection(*ic_nInd).FacetIndex < 0) this->GetContactNodeBuffer().push_back(*ic_nInd);
      this->GetNodeProjectionBuffer()[*ic_nInd] = r_threadProjection;
    }
    r_threadProjection.Reset();
  }	
  
  r_threadSlaveNodeList.clear();
  delete p_worker;
  assert(this->GetContactNodeBuffer().size() == tledHelper::MakeSortedUnique(this->GetContactNodeBuffer()).size());      
}

template <class TBaseTraverser>
void tledParallelBVHTraverserImplCPU<TBaseTraverser>::MergeEdgeEdgeNarrowPhaseResults(EdgeEdgeNarrowPhaseSubTraverser *p_worker) {
  std::vector<int> &r_threadSlaveEdgeList = p_worker->GetContactEdgeBuffer();

  p_worker->Join();
  for (std::vector<int>::const_iterator ic_eInd = r_threadSlaveEdgeList.begin(); ic_eInd < r_threadSlaveEdgeList.end(); ic_eInd++) {
    typename Superclass::MasterEdgeContactData &r_threadProjection = p_worker->GetEdgeProjectionBuffer()[*ic_eInd];

    if (this->GetEdgeProjection(*ic_eInd).EdgeIndex < 0 || this->IsBetterEdgeProjection(r_threadProjection.CollisionCoords, this->GetEdgeProjection(*ic_eInd).CollisionCoords)) {
      assert(this->GetEdgeProjectionBuffer()[*ic_eInd].EdgeIndex < 0 || std::find(this->GetContactEdgeBuffer().begin(), this->GetContactEdgeBuffer().end(), *ic_eInd) != this->GetContactEdgeBuffer().end());
      assert(this->GetEdgeProjectionBuffer()[*ic_eInd].EdgeIndex >= 0 || std::find(this->GetContactEdgeBuffer().begin(), this->GetContactEdgeBuffer().end(), *ic_eInd) == this->GetContactEdgeBuffer().end());

      if (this->GetEdgeProjectionBuffer()[*ic_eInd].EdgeIndex < 0) this->GetContactEdgeBuffer().push_back(*ic_eInd);
      this->GetEdgeProjectionBuffer()[*ic_eInd] = r_threadProjection;
    }
    r_threadProjection.Reset();
  }	

  r_threadSlaveEdgeList.clear();
  delete p_worker;
}

template <class TBaseTraverser>
void tledParallelBVHTraverserImplCPU<TBaseTraverser>::_ComputeNumberOfNodeFacetNarrowPhaseTests(int *p_dst, const tledParallelBVHTraverserImplCPU *pc_bt, const std::vector<int>::const_iterator ic_nodeIndBegin, const std::vector<int>::const_iterator ic_nodeIndEnd) {
  *p_dst = 0;
  for (std::vector<int>::const_iterator ic_n = ic_nodeIndBegin; ic_n < ic_nodeIndEnd; ic_n++) *p_dst += pc_bt->GetMasterFacetIndexBuffer()[*ic_n].size();
}

template <class TBaseTraverser>
void tledParallelBVHTraverserImplCPU<TBaseTraverser>::_ComputeNumberOfEdgeEdgeNarrowPhaseTests(int *p_dst, const tledParallelBVHTraverserImplCPU *pc_bt, const std::vector<int>::const_iterator ic_edgeIndBegin, const std::vector<int>::const_iterator ic_edgeIndEnd) {
  *p_dst = 0;
  for (std::vector<int>::const_iterator ic_e = ic_edgeIndBegin; ic_e < ic_edgeIndEnd; ic_e++) *p_dst += pc_bt->GetMasterEdgeIndexBuffer()[*ic_e].size();
}

template <class TBaseTraverser>
void tledParallelBVHTraverserImplCPU<TBaseTraverser>::RunNarrowPhase(void) {
  assert(this->GetContactNodeBuffer().size() == 0 && this->GetContactEdgeBuffer().size() == 0);
  assert(this->GetNumberOfBroadPhaseThreads() == 0);

  if ((int)(this->GetNodeFacetNarrowPhasePairs().size() + this->GetEdgeEdgeNarrowPhasePairs().size()) < 64*this->GetNumberOfThreads()) {
    Superclass::RunNarrowPhase();
  } else {
    using namespace tledHelper;  

    this->PreNodeFacetNarrowPhaseHook();
    if (this->GetNodeFacetNarrowPhasePairs().size() > 0) {
      const int szThreadBlock = this->GetThreadBlockSize(this->GetNodeFacetNarrowPhasePairs().size());

      std::vector<std::pair<int, int> >::const_iterator ic_bStart;
      std::vector<NodeFacetNarrowPhaseSubTraverser*> vp_workers;
    
      for (ic_bStart = this->GetNodeFacetNarrowPhasePairs().begin(); ic_bStart < this->GetNodeFacetNarrowPhasePairs().end() - szThreadBlock; ic_bStart += szThreadBlock) {
	vp_workers.push_back(this->LaunchNodeFacetNarrowPhaseWorker(ic_bStart, ic_bStart + szThreadBlock, vp_workers.size()));
      }
      assert((int)vp_workers.size() <= this->GetNumberOfThreads() - 1);
      this->ProcessNodeFacetNarrowPhaseItems(ic_bStart, this->GetNodeFacetNarrowPhasePairs().end());

      for (typename std::vector<NodeFacetNarrowPhaseSubTraverser*>::iterator ip_worker = vp_workers.begin(); ip_worker < vp_workers.end(); ip_worker++)  this->MergeNodeFacetNarrowPhaseResults(*ip_worker);
    }
    this->GetNodeFacetNarrowPhasePairs().clear();

    this->PreEdgeEdgeNarrowPhaseHook();
    if (this->GetEdgeEdgeNarrowPhasePairs().size() > 0) {      
      const int szThreadBlock = this->GetThreadBlockSize(this->GetEdgeEdgeNarrowPhasePairs().size());

      std::vector<std::pair<int, int> >::const_iterator ic_cEdgeInd;
      std::vector<EdgeEdgeNarrowPhaseSubTraverser*> vp_workers;

      for (ic_cEdgeInd = this->GetEdgeEdgeNarrowPhasePairs().begin(); ic_cEdgeInd < this->GetEdgeEdgeNarrowPhasePairs().end() - szThreadBlock; ic_cEdgeInd += szThreadBlock) {
	vp_workers.push_back(this->LaunchEdgeEdgeNarrowPhaseWorker(ic_cEdgeInd, ic_cEdgeInd + szThreadBlock, vp_workers.size()));
      }
      assert((int)vp_workers.size() <= this->GetNumberOfThreads() - 1);
      this->ProcessEdgeEdgeNarrowPhaseItems(ic_cEdgeInd, this->GetEdgeEdgeNarrowPhasePairs().end());
      for (typename std::vector<EdgeEdgeNarrowPhaseSubTraverser*>::iterator ip_worker = vp_workers.begin(); ip_worker < vp_workers.end(); ip_worker++) this->MergeEdgeEdgeNarrowPhaseResults(*ip_worker);
    }
    this->GetEdgeEdgeNarrowPhasePairs().clear();
  }
}

template <class TBaseTraverser>
void tledParallelBVHTraverserImplCPU<TBaseTraverser>::_NodeFacetStage1GeometryUpdateWorker(SlaveMesh *p_surface, const std::vector<int>::const_iterator ic_facetIndBegin, const std::vector<int>::const_iterator ic_facetIndEnd) {
  for (std::vector<int>::const_iterator ic_f = ic_facetIndBegin; ic_f < ic_facetIndEnd; ic_f++) {
    p_surface->GetNormalisedOldFacetNormalCached(*ic_f);
    p_surface->GetFacetProjectionOperatorCached(*ic_f);
  }
}

template <class TBaseTraverser>
void tledParallelBVHTraverserImplCPU<TBaseTraverser>::_NodeNormalUpdateWorker(SlaveMesh *p_surface, const std::vector<int>::const_iterator ic_nodeIndBegin, const std::vector<int>::const_iterator ic_nodeIndEnd) {
  for (std::vector<int>::const_iterator ic_n = ic_nodeIndBegin; ic_n < ic_nodeIndEnd; ic_n++) {
    p_surface->GetNodeNormalCached(*ic_n);
  }
}

template <class TBaseTraverser>
void tledParallelBVHTraverserImplCPU<TBaseTraverser>::UpdateNodeFacetDeformableMasterGeometry() {
  std::vector<int> masterFacetIndices;
  std::vector<boost::thread*> vp_threads;

  masterFacetIndices.reserve(this->GetNodeFacetNarrowPhasePairs().size());
  for (std::vector<std::pair<int, int> >::const_iterator ic_i = this->GetNodeFacetNarrowPhasePairs().begin(); ic_i < this->GetNodeFacetNarrowPhasePairs().end(); ic_i++) {
    masterFacetIndices.push_back(ic_i->second);
  }
  masterFacetIndices = tledHelper::MakeSortedUnique(masterFacetIndices);

  {
    const int szThreadBlock = this->GetThreadBlockSize(masterFacetIndices.size());

    std::vector<int>::const_iterator ic_n;

    for (ic_n = masterFacetIndices.begin(); ic_n < masterFacetIndices.end() - szThreadBlock; ic_n += szThreadBlock) {
      vp_threads.push_back(new boost::thread(_NodeFacetStage1GeometryUpdateWorker, &this->GetSlaveMesh(), ic_n, ic_n + szThreadBlock));
    }
    _NodeFacetStage1GeometryUpdateWorker(&this->GetSlaveMesh(), ic_n, masterFacetIndices.end());

    this->JoinThreads(vp_threads);
  }    
}

template <class TBaseTraverser>
void tledParallelBVHTraverserImplCPU<TBaseTraverser>::_LaunchNodeUpdateWorker(const std::vector<int> &nodeIndices) {
  const int szThreadBlock = this->GetThreadBlockSize(nodeIndices.size());

  std::vector<boost::thread*> vp_threads;
  std::vector<int>::const_iterator ic_n;

  for (ic_n = nodeIndices.begin(); ic_n < nodeIndices.end() - szThreadBlock; ic_n += szThreadBlock) {
    vp_threads.push_back(new boost::thread(_NodeNormalUpdateWorker, &this->GetSlaveMesh(), ic_n, ic_n + szThreadBlock));
  }
  _NodeNormalUpdateWorker(&this->GetSlaveMesh(), ic_n, nodeIndices.end());

  this->JoinThreads(vp_threads);
}

template <class TBaseTraverser>
void tledParallelBVHTraverserImplCPU<TBaseTraverser>::UpdateEdgeEdgeDeformableSlaveGeometry() {
  std::vector<int> slaveIndices;

  slaveIndices.clear();
  slaveIndices.reserve(this->GetEdgeEdgeNarrowPhasePairs().size()*2);
  for (std::vector<std::pair<int, int> >::const_iterator ic_i = this->GetEdgeEdgeNarrowPhasePairs().begin(); ic_i < this->GetEdgeEdgeNarrowPhasePairs().end(); ic_i++) {
      slaveIndices.push_back(this->GetSlaveMesh().GetEdge(ic_i->first).first);
      slaveIndices.push_back(this->GetSlaveMesh().GetEdge(ic_i->first).second);
  }

  _LaunchNodeUpdateWorker(tledHelper::MakeSortedUnique(slaveIndices));
}

template <class TBaseTraverser>
void tledParallelBVHTraverserImplCPU<TBaseTraverser>::UpdateEdgeEdgeDeformableMasterGeometry() {
  std::vector<int> masterIndices;

  masterIndices.clear();
  masterIndices.reserve(this->GetEdgeEdgeNarrowPhasePairs().size()*2);
  for (std::vector<std::pair<int, int> >::const_iterator ic_i = this->GetEdgeEdgeNarrowPhasePairs().begin(); ic_i < this->GetEdgeEdgeNarrowPhasePairs().end(); ic_i++) {
      masterIndices.push_back(this->GetSlaveMesh().GetEdge(ic_i->second).first);
      masterIndices.push_back(this->GetSlaveMesh().GetEdge(ic_i->second).second);
  }
  _LaunchNodeUpdateWorker(tledHelper::MakeSortedUnique(masterIndices));
}

template <class TBaseTraverser>
void tledParallelBVHTraverserImplCPU<TBaseTraverser>::UpdateEdgeEdgeDeformableMasterSlaveGeometry() {
  std::vector<int> nodeIndices;

  nodeIndices.clear();
  nodeIndices.reserve(this->GetEdgeEdgeNarrowPhasePairs().size()*4);
  for (std::vector<std::pair<int, int> >::const_iterator ic_i = this->GetEdgeEdgeNarrowPhasePairs().begin(); ic_i < this->GetEdgeEdgeNarrowPhasePairs().end(); ic_i++) {
    nodeIndices.push_back(this->GetSlaveMesh().GetEdge(ic_i->second).first);
    nodeIndices.push_back(this->GetSlaveMesh().GetEdge(ic_i->second).second);
    nodeIndices.push_back(this->GetSlaveMesh().GetEdge(ic_i->first).first);
    nodeIndices.push_back(this->GetSlaveMesh().GetEdge(ic_i->first).second);
  }
  _LaunchNodeUpdateWorker(tledHelper::MakeSortedUnique(nodeIndices));
}

template <class TBaseTraverser>
void tledParallelBVHTraverserImplCPU<TBaseTraverser>::Init(tledUnstructuredContactManager &r_manager) {
  Superclass::Init(r_manager);
}

template <class TBaseTraverser>
void tledParallelBVHTraverserImplCPU<TBaseTraverser>::InitNarrowPhaseSubTraverser(NarrowPhaseSubTraverser &r_sub) const {
  r_sub.SetNarrowPhaseMaxDistance(this->GetNarrowPhaseMaxDistance());
  r_sub.SetDoMaster(this->DoMaster());
}

template <class TBaseTraverser>
void tledParallelBVHTraverserImplCPU<TBaseTraverser>::_SortUniqueWorker(std::vector<std::pair<int, int> >::iterator *pi_dstBegin, std::vector<std::pair<int, int> >::iterator *pi_dstEnd) {
  std::sort(*pi_dstBegin, *pi_dstEnd);
  *pi_dstEnd = std::unique(*pi_dstBegin, *pi_dstEnd);
}

template <class TBaseTraverser>
void tledParallelBVHTraverserImplCPU<TBaseTraverser>::_UniqueRangeMergeWorker(std::vector<std::pair<int, int> >::iterator *pi_dstBegin, std::vector<std::pair<int, int> >::iterator *pi_dstEnd, const std::vector<std::pair<int, int> >::iterator src0Begin, const std::vector<std::pair<int, int> >::iterator src0End, const std::vector<std::pair<int, int> >::iterator src1Begin, const std::vector<std::pair<int, int> >::iterator src1End) {
  assert(std::unique(src0Begin, src0End) == src0End && std::unique(src1Begin, src1End) == src1End);
  *pi_dstEnd = std::set_union(src0Begin, src0End, src1Begin, src1End, *pi_dstBegin, tledBVHTraverserCPU::NarrowPhaseOrdering());  
  assert(*pi_dstEnd - *pi_dstBegin >= std::max(src1End - src1Begin, src0End - src0Begin));
}

template <class TBaseTraverser>
void tledParallelBVHTraverserImplCPU<TBaseTraverser>::_MakeUniqueParallel(std::vector<std::pair<int, int> > &r_range) const {
  typedef std::vector<std::pair<int, int> >::iterator __RangeIterator;
  typedef std::pair<__RangeIterator, __RangeIterator> __Range;

  std::vector<boost::thread*> vp_threads;
  std::vector<__Range> subRanges(this->GetNumberOfThreads());
  std::vector<std::pair<int, int> > outAlt(r_range.size());
  int numRanges = 0;
#ifndef NDEBUG
  std::vector<std::pair<int, int> > fullRange = tledHelper::MakeSortedUnique(r_range, std::equal_to<std::pair<int, int> >(),  tledBVHTraverserCPU::NarrowPhaseOrdering());
#endif
  
  {
    const int blockSize = this->GetThreadBlockSize(r_range.size());

    __RangeIterator i_start;
      
    assert(blockSize > 1);
    for (i_start = r_range.begin(); i_start < r_range.end() - blockSize; i_start += blockSize, numRanges++) {
      subRanges[numRanges].first = i_start;
      subRanges[numRanges].second = i_start + blockSize;
      vp_threads.push_back(new boost::thread(_SortUniqueWorker, &subRanges[numRanges].first, &subRanges[numRanges].second)); 
    }
    
    subRanges[numRanges].first = i_start;
    subRanges[numRanges].second = r_range.end();
    _SortUniqueWorker(&subRanges[numRanges].first, &subRanges[numRanges].second);
    numRanges++;

    this->JoinThreads(vp_threads);
#ifndef NDEBUG
    for (std::vector<__Range>::const_iterator ic_r = subRanges.begin(); ic_r < subRanges.begin() + numRanges; ic_r++) {
      assert(ic_r->second > ic_r->first);
    }
#endif
  }

  while (numRanges > 1) {    
    std::vector<std::pair<int, int> > &r_dstVector = subRanges[0].first == r_range.begin()? outAlt : r_range;
    __RangeIterator i_out = r_dstVector.begin();
    std::vector<__Range> newRanges(subRanges.size()/2 + 1);
    int numNewRanges = 0;
    
    assert(vp_threads.size() == 0);
    for (int s = 0; s + 1 < numRanges; s += 2) {
      const __Range src0 = subRanges[s], src1 = subRanges[s+1];      

      newRanges[numNewRanges].first = i_out;
      assert(src0.second > src0.first && src1.second > src1.first);
      vp_threads.push_back(new boost::thread(_UniqueRangeMergeWorker, &newRanges[numNewRanges].first, &newRanges[numNewRanges].second, src0.first, src0.second, src1.first, src1.second));
      i_out += (src0.second - src0.first) + (src1.second - src1.first);
      numNewRanges += 1;
    }

    if (numRanges%2 == 1) {
      newRanges[numNewRanges].first = i_out;
      newRanges[numNewRanges].second = std::copy(subRanges[numRanges-1].first, subRanges[numRanges-1].second, newRanges[numNewRanges].first);
      numNewRanges += 1;
    }
    
    assert(numNewRanges <= (int)newRanges.size() && numNewRanges >= numRanges/2);
    this->JoinThreads(vp_threads);
#ifndef NDEBUG
    for (std::vector<__Range>::const_iterator ic_r = newRanges.begin(); ic_r < newRanges.begin() + numNewRanges; ic_r++) {
      assert(ic_r->second > ic_r->first);
    }
#endif
    subRanges = newRanges;
    numRanges = numNewRanges;
  }
  assert(numRanges == 1);
  if (r_range.begin() != subRanges[0].first) {
    std::copy(subRanges[0].first, subRanges[0].second, r_range.begin());
    subRanges[0].second = r_range.begin() + (subRanges[0].second - subRanges[0].first);
  }
  r_range.erase(subRanges[0].second, r_range.end());
  assert(std::unique(r_range.begin(), r_range.end()) == r_range.end());
  assert(r_range.size() == fullRange.size() && std::equal(fullRange.begin(), fullRange.end(), r_range.begin()));
}

template <class TBaseTraverser>
void tledParallelBVHTraverserImplCPU<TBaseTraverser>::FinishBroadPhase() {
  if ((int)this->GetNodeFacetNarrowPhasePairs().size() > 128*this->GetNumberOfThreads() && (int)this->GetEdgeEdgeNarrowPhasePairs().size() > 128*this->GetNumberOfThreads()) {
    _MakeUniqueParallel(this->GetNodeFacetNarrowPhasePairs());
    _MakeUniqueParallel(this->GetEdgeEdgeNarrowPhasePairs());
  } else Superclass::FinishBroadPhase();
}
