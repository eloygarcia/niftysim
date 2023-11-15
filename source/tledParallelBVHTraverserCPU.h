// =========================================================================
// File:       tledParallelBVHTraverserCPU.h
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
#ifndef tledParallelBVHTraverserCPU_H
#define tledParallelBVHTraverserCPU_H

#include "tledBVHTraverserCPU.h"
#include "tledParallel.h"
#include "tledUnstructuredContactManager.h"

#include <utility>
#include <cassert>
#include <deque>
#include <algorithm>
#include <vector>

#ifdef _USE_BOOST_

/**
 * \brief Shared resources for CPU-parallel BVH traversers
 * \ingroup contact
 */
class tledParallelBVHTraverserCPU : public tledParallel {
  /**
   * \name Memory Management
   * @{
   */
private:
  static std::vector<std::vector<int>*> svpv_ThreadSlaveNodeIndices, svpv_ThreadSlaveEdgeIndices;
  static std::vector<std::vector<tledBVHTraverserCPU::MasterFacetContactData>*> svpv_ThreadNodeProjectionBuffer;
  static std::vector<std::vector<tledBVHTraverserCPU::MasterEdgeContactData>*> svpv_ThreadEdgeProjectionBuffer;
  static std::vector<std::vector<std::pair<int, int> >*> svpv_ThreadNodeFacetNarrowPhasePairs, svpv_ThreadEdgeEdgeNarrowPhasePairs;

public:
  /** \brief Node-facet projection buffers for multi-threaded traversers */
  std::vector<tledBVHTraverserCPU::MasterFacetContactData>& GetThreadNodeProjectionBuffer(const int threadInd);
  /** \brief Edge-edge projection buffers for multi-threaded traversers */
  std::vector<tledBVHTraverserCPU::MasterEdgeContactData>& GetThreadEdgeProjectionBuffer(const int threadInd);

  /** \brief Slave-node buffers for multi-threaded traversers */
  std::vector<int>& GetThreadContactNodeBuffer(const int threadInd);
  /** \brief Slave-edge buffers for multi-threaded traversers */
  std::vector<int>& GetThreadContactEdgeBuffer(const int threadInd);

  std::vector<std::pair<int, int> >& GetThreadNodeFacetNarrowPhasePairBuffer(const int threadInd);
  std::vector<std::pair<int, int> >& GetThreadEdgeEdgeNarrowPhasePairBuffer(const int threadInd);

  virtual int GetMaxNumberOfSlaveNodes(void) const = 0;
  virtual int GetMaxNumberOfSlaveEdges(void) const = 0;

public:
  /** Explicitly deallocate shared resources */
  static void DeallocateSharedBuffers(void);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  virtual ~tledParallelBVHTraverserCPU(void) {}
  /** @} */
};

/**
 * \brief Adds CPU parallelism to a BVH traverser
 * \ingroup contact
 */
template <class TBaseTraverser>
class tledParallelBVHTraverserImplCPU : public TBaseTraverser, public tledParallelBVHTraverserCPU {
  /**
   * \name Types
   * @{
   */
public:
  typedef TBaseTraverser Superclass;
  typedef typename Superclass::MasterBVH MasterBVH;
  typedef typename Superclass::SlaveBVH SlaveBVH;
  typedef typename Superclass::MasterMesh MasterMesh;
  typedef typename Superclass::SlaveMesh SlaveMesh;
  typedef tledBVHTraverserCPU::MasterEdgeContactData MasterEdgeContactData;
  typedef tledBVHTraverserCPU::MasterFacetContactData MasterFacetContactData;

protected:
  class TraverserWorker;

  class SubTraverserAPI;
  class BroadPhaseSubTraverser;

  class NarrowPhaseSubTraverser;
  class NodeFacetNarrowPhaseSubTraverser;
  class EdgeEdgeNarrowPhaseSubTraverser;
  /** @} */

  /**
   * \name Detection
   * @{
   */
private:
  int m_NumBPSubThreads;

private:
  static void _NodeFacetStage1GeometryUpdateWorker(SlaveMesh *p_surface, const std::vector<int>::const_iterator ic_facetIndBegin, const std::vector<int>::const_iterator ic_facetIndEnd);
  static void _NodeNormalUpdateWorker(SlaveMesh *p_surface, const std::vector<int>::const_iterator ic_nodeIndBegin, const std::vector<int>::const_iterator ic_nodeIndEnd);
  static void _ComputeNumberOfNodeFacetNarrowPhaseTests(int *p_dst, const tledParallelBVHTraverserImplCPU *pc_bt, const std::vector<int>::const_iterator ic_nodeIndBegin, const std::vector<int>::const_iterator ic_nodeIndEnd);
  static void _ComputeNumberOfEdgeEdgeNarrowPhaseTests(int *p_dst, const tledParallelBVHTraverserImplCPU *pc_bt, const std::vector<int>::const_iterator ic_edgeIndBegin, const std::vector<int>::const_iterator ic_edgeIndEnd);
  
  static void _UniqueRangeMergeWorker(std::vector<std::pair<int, int> >::iterator *pi_dstBegin, std::vector<std::pair<int, int> >::iterator *pi_dstEnd, const std::vector<std::pair<int, int> >::iterator src0Begin, const std::vector<std::pair<int, int> >::iterator src0End, const std::vector<std::pair<int, int> >::iterator src1Begin, const std::vector<std::pair<int, int> >::iterator src1End);
  static void _SortUniqueWorker(std::vector<std::pair<int, int> >::iterator *pi_dstBegin, std::vector<std::pair<int, int> >::iterator *pi_dstEnd);
  void _MakeUniqueParallel(std::vector<std::pair<int, int> > &r_range) const;

  void _LaunchNodeUpdateWorker(const std::vector<int> &nodeIndices);

protected: 
  int GetNumberOfBroadPhaseThreads(void) const { return m_NumBPSubThreads; }
  void IncNumberOfBroadPhaseThreads(void) { m_NumBPSubThreads += 1; }
  void ResetBroadPhaseThreadCounter(void) { m_NumBPSubThreads = 0; }  

  virtual void RunBroadPhase(void) = 0;
  virtual void RunNarrowPhase(void);

  virtual void FinishBroadPhase(void);

  /** 
   * Updates the geometry required for node-facet narrow-phase tests with deformable, lazy-updating masters. 
   * Attention: Since this requires read-write access it is assumed that the object's slave mesh is the master for this particular narrow-phase.
   */
  void UpdateNodeFacetDeformableMasterGeometry(void);

  /** 
   * Updates the geometry required for edge-edge narrow-phase tests with deformable, lazy-updating masters. 
   * Attention: Since this requires read-write access it is assumed that the object's slave mesh is the master for this particular narrow-phase.
   */
  void UpdateEdgeEdgeDeformableMasterGeometry(void);

  /** 
   * Updates the geometry required for edge-edge narrow-phase tests with deformable, lazy-updating slaves. 
   */
  void UpdateEdgeEdgeDeformableSlaveGeometry(void);

  /** 
   * Updates the geometry required for edge-edge narrow-phase tests with deformable geometry only. 
   */
  void UpdateEdgeEdgeDeformableMasterSlaveGeometry(void);

  /** Performs a pair-wise, breadth-first descent in the BVHs until enough BV pairs have been collected to proceed in parallel with the desired number of threads. */
  virtual std::vector<BroadPhaseSubTraverser*> DescendUntilFullyLoaded(const int slaveStartInd, const int masterStartInd);

  /** Converts the output of the BroadPhaseSubTraverser objects to something that can be used in narrow-phase tests, cleans up the memory used in the broad-phase */
  virtual void CollectBroadPhaseSubResults(std::vector<BroadPhaseSubTraverser*> &rvp_subTraversers);

  virtual NodeFacetNarrowPhaseSubTraverser* LaunchNodeFacetNarrowPhaseWorker(const std::vector<std::pair<int, int> >::const_iterator itemsBegin, const std::vector<std::pair<int, int> >::const_iterator itemsEnd, const int threadInd) = 0;
  virtual EdgeEdgeNarrowPhaseSubTraverser* LaunchEdgeEdgeNarrowPhaseWorker(const std::vector<std::pair<int, int> >::const_iterator itemsBegin, const std::vector<std::pair<int, int> >::const_iterator itemsEnd, const int threadInd) = 0;
  
  /** Initialises the search parameters on a sub-traverser object */
  void InitNarrowPhaseSubTraverser(NarrowPhaseSubTraverser &r_sub) const;

  /** Called before the node-facet narrow-phase is started, can be used for carrying out any necessary geometry updates */
  virtual void PreNodeFacetNarrowPhaseHook(void) {}

  /** Called before the edge-edge narrow-phase is started, can be used for carrying out any necessary geometry updates */
  virtual void PreEdgeEdgeNarrowPhaseHook(void) {}

  /** Merges the node-facet narrow-phase results of sub-object with those of the object on which the function is called, also deletes the p_worker. */
  void MergeNodeFacetNarrowPhaseResults(NodeFacetNarrowPhaseSubTraverser *p_worker);

  /** Merges the node-facet narrow-phase results of sub-object with those of the object on which the function is called, also deletes the p_worker. */
  void MergeEdgeEdgeNarrowPhaseResults(EdgeEdgeNarrowPhaseSubTraverser *p_worker);

  /** 
   * \brief Launches a new broad-phase worker thread.
   *
   * This member function itself cannot be assumed to be thread-safe! 
   */
  virtual BroadPhaseSubTraverser* LaunchSubTraverser(const int slaveStartInd, const int masterStartInd);

  virtual int GetMaxNumberOfSlaveNodes(void) const { return std::max(this->GetMasterMesh().GetNumberOfNodes(), this->GetSlaveMesh().GetNumberOfNodes()); }
  virtual int GetMaxNumberOfSlaveEdges(void) const { return std::max(this->GetMasterMesh().GetNumberOfEdges(), this->GetSlaveMesh().GetNumberOfEdges()); }
  /** @} */

  /**
   * \name Construction, Destruction, Initialisation
   * @{
   */
public:
  virtual void Init(tledUnstructuredContactManager &r_manager);

  tledParallelBVHTraverserImplCPU(SlaveBVH &r_slaveBVH, const MasterBVH &masterBVH) : Superclass(r_slaveBVH, masterBVH), m_NumBPSubThreads(0) {}
  virtual ~tledParallelBVHTraverserImplCPU(void) {}
  /** @} */
};

template <class TBaseTraverser>
class tledParallelBVHTraverserImplCPU<TBaseTraverser>::TraverserWorker {
  /**
   * \name Thread Control
   * @{
   */
private:
  boost::thread *mp_Thread;

public:
  void SetBoostThread(boost::thread &r_thread) { mp_Thread = &r_thread; }
  void Join(void);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  TraverserWorker(void) : mp_Thread(NULL) {}
  virtual ~TraverserWorker(void) {}
  /** @} */
};

template <class TBaseTraverser>
class tledParallelBVHTraverserImplCPU<TBaseTraverser>::NarrowPhaseSubTraverser : public TBaseTraverser, public TraverserWorker {
  /**
   * \name Memory Management
   * @{
   */
public:
  virtual std::vector<std::pair<int, int> >& GetNodeFacetNarrowPhasePairs(void) { std::abort(); return TBaseTraverser::GetNodeFacetNarrowPhasePairs(); }
  virtual std::vector<std::pair<int, int> >& GetEdgeEdgeNarrowPhasePairs(void) { std::abort(); return TBaseTraverser::GetEdgeEdgeNarrowPhasePairs(); }
  virtual const std::vector<std::pair<int, int> >& GetNodeFacetNarrowPhasePairs(void) const { std::abort(); return TBaseTraverser::GetNodeFacetNarrowPhasePairs(); }
  virtual const std::vector<std::pair<int, int> >& GetEdgeEdgeNarrowPhasePairs(void) const { std::abort(); return TBaseTraverser::GetNodeFacetNarrowPhasePairs(); }

  virtual std::vector<int>& GetContactNodeBuffer(void) { std::abort(); return TBaseTraverser::GetContactNodeBuffer(); }
  virtual std::vector<int>& GetContactEdgeBuffer(void) { std::abort(); return TBaseTraverser::GetContactEdgeBuffer(); }

  virtual std::vector<tledBVHTraverserCPU::MasterEdgeContactData>& GetEdgeProjectionBuffer(void) { std::abort(); return TBaseTraverser::GetEdgeProjectionBuffer(); }
  virtual std::vector<tledBVHTraverserCPU::MasterFacetContactData>& GetNodeProjectionBuffer(void) { std::abort(); return TBaseTraverser::GetNodeProjectionBuffer(); }
  virtual const std::vector<MasterEdgeContactData>& GetEdgeProjectionBuffer(void) const { std::abort(); return TBaseTraverser::GetEdgeProjectionBuffer(); }
  virtual const std::vector<MasterFacetContactData>& GetNodeProjectionBuffer(void) const { std::abort(); return TBaseTraverser::GetNodeProjectionBuffer(); }

public:
  virtual const std::vector<int>& GetSlaveNodeIndices(void) const { std::abort(); return TBaseTraverser::GetSlaveNodeIndices(); }
  virtual const std::vector<int>& GetSlaveEdgeIndices(void) const { std::abort(); return TBaseTraverser::GetSlaveEdgeIndices(); }
  /** @} */

public:
  NarrowPhaseSubTraverser(SlaveBVH &r_slaveBVH, const MasterBVH &masterBVH) : Superclass(r_slaveBVH, masterBVH) {}
  virtual ~NarrowPhaseSubTraverser(void) {}
};

/**
 * \brief Node-facet narrow-phase worker base class 
 */
template <class TBaseTraverser>
class tledParallelBVHTraverserImplCPU<TBaseTraverser>::NodeFacetNarrowPhaseSubTraverser : public NarrowPhaseSubTraverser {
  /**
   * \name Memory Management
   * @{
   */
private:
  std::vector<int> &mr_ContactNodeIndices;
  std::vector<MasterFacetContactData> &mr_MasterFacets;  

protected:
  virtual std::vector<int>& GetContactNodeBuffer(void) { return mr_ContactNodeIndices; }
  virtual std::vector<tledBVHTraverserCPU::MasterFacetContactData>& GetNodeProjectionBuffer(void) { return mr_MasterFacets; }
  virtual const std::vector<MasterFacetContactData>& GetNodeProjectionBuffer(void) const { return mr_MasterFacets; }

public:
  virtual const std::vector<int>& GetSlaveNodeIndices(void) const { return mr_ContactNodeIndices; }
  /** @} */

  /**
   * \name Processing
   * @{
   */
public:
  static void ItemProcessorWorker(NodeFacetNarrowPhaseSubTraverser *p_traverser, const std::vector<std::pair<int, int> >::const_iterator itemsBegin, const std::vector<std::pair<int, int> >::const_iterator itemsEnd) { p_traverser->ProcessNodeFacetNarrowPhaseItems(itemsBegin, itemsEnd); }
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:  
  NodeFacetNarrowPhaseSubTraverser(std::vector<int> &r_slaves, std::vector<MasterFacetContactData> &r_masters, SlaveBVH &r_slaveBVH, const MasterBVH &masterBVH) : NarrowPhaseSubTraverser(r_slaveBVH, masterBVH), mr_ContactNodeIndices(r_slaves), mr_MasterFacets(r_masters) {}
  virtual ~NodeFacetNarrowPhaseSubTraverser(void) {}
  /** @} */

  friend class tledParallelBVHTraverserImplCPU;
};

/**
 * \brief Edge-edge narrow-phase worker base class 
 */
template <class TBaseTraverser>
class tledParallelBVHTraverserImplCPU<TBaseTraverser>::EdgeEdgeNarrowPhaseSubTraverser : public NarrowPhaseSubTraverser {
  /**
   * \name Memory Management
   * @{
   */
private:
  std::vector<int> &mr_ContactEdgeIndices;
  std::vector<MasterEdgeContactData> &mr_MasterEdges;

protected:
  virtual std::vector<int>& GetContactEdgeBuffer(void) { return mr_ContactEdgeIndices; }
  virtual std::vector<tledBVHTraverserCPU::MasterEdgeContactData>& GetEdgeProjectionBuffer(void) { return mr_MasterEdges; }
  virtual const std::vector<MasterEdgeContactData>& GetEdgeProjectionBuffer(void) const { return mr_MasterEdges; }

public:
  virtual const std::vector<int>& GetSlaveEdgeIndices(void) const { return mr_ContactEdgeIndices; }
  /** @} */

  /**
   * \name Processing
   * @{
   */
public:
  static void ItemProcessorWorker(EdgeEdgeNarrowPhaseSubTraverser *p_traverser, const std::vector<std::pair<int, int> >::const_iterator itemsBegin, const std::vector<std::pair<int, int> >::const_iterator itemsEnd) { p_traverser->ProcessEdgeEdgeNarrowPhaseItems(itemsBegin, itemsEnd); }
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:  
  EdgeEdgeNarrowPhaseSubTraverser(std::vector<int> &r_slaves, std::vector<MasterEdgeContactData> &r_masters, SlaveBVH &r_slaveBVH, const MasterBVH &masterBVH) : NarrowPhaseSubTraverser(r_slaveBVH, masterBVH), mr_ContactEdgeIndices(r_slaves), mr_MasterEdges(r_masters) {}
  virtual ~EdgeEdgeNarrowPhaseSubTraverser(void) {}
  /** @} */

  friend class tledParallelBVHTraverserImplCPU;
};

#include "tledParallelBVHTraverserCPU.tpp"
#endif
#endif
