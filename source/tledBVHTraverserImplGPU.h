// =========================================================================
// File:       tledBVHTraverserImplGPU.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    January 2015
// 
// Copyright (c) 2015, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledBVHTraverserImplGPU_H
#define tledBVHTraverserImplGPU_H

#include "tledBVHTraverserGPU.h"

#include <algorithm>
#include <cassert>

/**
 * \brief Base class for GPU contact search classes.
 * \ingroup contact
 * 
 * Public API, template parameter TAPI, must be covariant with tledBVHTraverserGPU and tledBVHTraverser
 */
template <class TMasterBVH, class TSlaveBVH, class TAPI>
class tledBVHTraverserImplGPU : public tledBVHTraverserImpl<TMasterBVH, TSlaveBVH, TAPI> {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledBVHTraverserImpl<TMasterBVH, TSlaveBVH, TAPI> Superclass;
  typedef TMasterBVH MasterBVH;
  typedef TSlaveBVH SlaveBVH;
  typedef typename TMasterBVH::ContactMesh MasterMesh;
  typedef typename TSlaveBVH::ContactMesh SlaveMesh;
  /** @} */

  /**
   * \name Meshes
   * @{
   */
protected:
  typename SlaveMesh::GPUSurface* GetSlaveGPUSurface(void) { return static_cast<typename SlaveMesh::GPUSurface*>(this->GetSlaveMesh().GetDeviceGPUSurface()); }
  const typename SlaveMesh::GPUSurface* GetSlaveGPUSurface(void) const { return static_cast<const typename SlaveMesh::GPUSurface*>(this->GetSlaveMesh().GetDeviceGPUSurface()); }
  const typename MasterMesh::GPUSurface* GetMasterGPUSurface(void) const { return static_cast<const typename MasterMesh::GPUSurface*>(this->GetMasterMesh().GetDeviceGPUSurface()); }
  /** @} */

  /**
   * \name Results
   * @{ 
   */
private:
  int *mdp_NumNodeFacetItems, *mdp_NumEdgeEdgeItems, *mhp_NumResults;

protected:
  /** General-purpose pinned copy buffer for result counters */
  int* GetHostItemCounter(void) { return mhp_NumResults; }

  /** Fetches a counter value from the device */
  int RetrieveCounterValue(const int *dpc_counter);

  /** Fetches a counter value from the device (will have synchronise the given stream, but other streams proceed unaffected) */
  int RetrieveCounterValue(const int *dpc_counter, const cudaStream_t stream);

  /** 
   * \brief Pointer to counter used in node-facet narrow-phase.
   */
  int* GetOnDeviceNodeFacetItemCounter(void) { return mdp_NumNodeFacetItems; }
  
  /** 
   * Pointer to counter used in edge-edge narrow-phase.
   */
  int* GetOnDeviceEdgeEdgeItemCounter(void) { return mdp_NumEdgeEdgeItems; }
  /** @} */

  /**
   * \name Broad-Phase
   * @{
   */
private:
  int m_LastMaxBPBufferCapacity;
  int *mdpp_BPPairCounters[2];
  int2 *mdp_BPSearchPairs;
  int m_NumSearchPairs, m_NumSearchLevels;

protected:
  /** Converts a list of master/slave leaf BV pairs to a list of corresponding surface primitives. */
  virtual void ConvertBVToPrimitivePairs(int2 *dp_list, const int *dpc_numItems, const int numItmesUB);

  /** Setter for descendant provided search settings, has to be called at least once during initialisation. */
  virtual void SetBroadPhaseSearchSettings(const int2 *searchBVPairs, const int numSearchBVPairs, const int numSearchLevels);

  /** Hook provided for descendant setup of BP search settings. */
  virtual void InitStartBVPairs(void) = 0;

  /** \brief Capacity required for initial list of broad-phase contact search BV pairs. */
  int GetNumberOfInitialBroadPhaseSearchPairs(void) const { return m_NumSearchPairs; }
  
  /** \brief Returns max. BVH depth number of levels that have to be visited in current master BVH */
  int GetNumberOfSearchBVHLevels(void) const { return m_NumSearchLevels; }

  /** \brief Generates the initial list BV pairs to be tested by descending to root level in the slave BVH */
  static std::vector<int2> GenerateBroadPhaseCheckPairs(const std::vector<int> &slaveLeaves, const int masterRootInd);

  /** \brief Starts a broad-phase contact search. Can be used as hook for any start-up tasks, however parent function must be called by sub-classes. */
  virtual void InitialiseBroadPhaseSearch(int2 *dp_dst, int *dp_counter);

  /** \brief Filters out non-intersecting bv pairs */
  virtual void FilterNonIntersectingBroadPhaseBVPairs(int2 *dp_outList, int *dp_outCounter, const int2 *dpc_inList, const int *dpc_inCounter, const int maxNumItems, const cudaStream_t stream);

  /** \brief Descends one level in the master BVH during broad-phase traversal */
  void DescendInMasterBVH(int2 *dp_outList, int *dp_outCounter, const int2 *dpc_inList, const int *dpc_inCounter, const int maxNumItems, const cudaStream_t stream);

  /** Post-processing of narrow-phase candidate index pairs */
  virtual void PrepareNarrowPhase(const int2 *dpc_inBPList, const int *dpc_actualNumResultItems, const int maxBPResultItems);

  /** Invokes the node-facet broad-phase result post-processing kernel. Asynchronous operation preferred. */
  virtual void InvokeNodeFacetBroadPhaseQueueKernel(int2 *dp_out, int *dp_counter, cudaStream_t stream, const int2 *dpc_inItems, const int *dpc_numItems, const int maxNumItems);

  /** Invokes the edge-edge broad-phase result post-processing kernel. Asynchronous operation preferred. */
  virtual void InvokeEdgeEdgeBroadPhaseQueueKernel(int2 *dp_out, int *dp_counter, cudaStream_t stream, const int2 *dpc_inItems, const int *dpc_numItems, const int maxNumItems);

  virtual void RunBroadPhase(void);
  /** @} */

  /**
   * \name Narrow Phase
   * @{
   */
public:
  /** First intermediate result in node to facet projection */
  struct NodeFacetInitialProjection {
    int ContactNodeIndices[SlaveMesh::Facet::NumberOfVertices+1];
    float3 Xi;
  };

  /** First intermediate result in edge-edge narrow phase */
  struct EdgeEdgeInitialProjection {
    int2 SlaveEdge, MasterEdge;
    float3 A1, B1;
    float3 C1, D1;
    float2 Xi;
  };

  typedef tledBVHTraverserGPU::NodeFacetNarrowPhaseResult<MasterMesh::Facet::NumberOfVertices> NodeFacetNarrowPhaseResult;
  typedef tledBVHTraverserGPU::EdgeEdgeNarrowPhaseResult EdgeEdgeNarrowPhaseResult;

protected:
  /** Launches the kernel(s) for the first stage of node-facet narrow-phase testing */
  virtual tledCUDADeviceMemoryBlock* InvokeNodeFacetNarrowPhaseStage1(int *dp_outItemCounter, cudaStream_t stream, const int2 *inItems, const int numNodeFacet) = 0;

  /** Launches the kernel(s) for the first stage of edge-edge narrow-phase testing */
  virtual tledCUDADeviceMemoryBlock* InvokeEdgeEdgeNarrowPhaseStage1(int *dp_outItemCounter, cudaStream_t stream, const int2 *inItems, const int numEdgeEdge) = 0;

  /** Launches the kernel(s) for the second stage of node-facet narrow-phase testing */
  virtual tledCUDADeviceMemoryBlock* InvokeNodeFacetNarrowPhaseStage2(int *dp_outItemCounter, cudaStream_t stream, const NodeFacetInitialProjection *inItems, const int numNodeFacet) = 0;

  /** Launches the kernel(s) for the second stage of edge-edge narrow-phase testing */
  virtual tledCUDADeviceMemoryBlock* InvokeEdgeEdgeNarrowPhaseStage2(int *dp_outItemCounter, cudaStream_t stream, const EdgeEdgeInitialProjection *inItems, const int numEdgeEdge) = 0;

  virtual void RunNarrowPhase(void);
  /** @} */

  /**
   * \name Main Functions
   * @{
   */
public:
  virtual void FindCollisions(void);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  /** Final initialisation after BVH has been constructed etc. */
  virtual void Init(tledUnstructuredContactManager &r_manager);

  tledBVHTraverserImplGPU(TSlaveBVH &r_slaveBVH, const TMasterBVH &masterBVH);
  virtual ~tledBVHTraverserImplGPU(void);
  /** @} */
};

#if defined __CUDACC__ && !defined  __GPU_TEST_LINK_CONFLICT_NO_INCLUDE
#include "tledBVHTraverserImplGPU.tpp"
#endif

#endif
