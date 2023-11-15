// =========================================================================
// File:       tledBVHTraverserImplGPU.tpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   CUDA
// Created:    January 2015
// 
// Copyright (c) 2015, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

#include "tledHelper.h"
#include "tledCUDAHelpers.h"
#include "tledComplexCUDAHelpers.h"

#include "tledBVHTraverserGPU_kernels.tpp"

template <class TMasterBVH, class TSlaveBVH, class TAPI>
void tledBVHTraverserImplGPU<TMasterBVH, TSlaveBVH, TAPI>::InvokeNodeFacetBroadPhaseQueueKernel(int2 *dp_out, int *dp_counter, cudaStream_t stream, const int2 *inItems, const int *dpc_numItems, const int maxNumItems) {
  const int blockSize = 128;
  const int numBlks = tledCUDAHelpers::GetNumberOfBlocks(maxNumItems, blockSize);

  tledLogDebugStream(tledHelper::Info() << "Post-processing <= "<< maxNumItems << " node-facet broad-phase results in " << (this->DoMaster()? "master" : "slave") << " mode");
  assert(maxNumItems >= this->RetrieveCounterValue(dpc_numItems));
  if (this->DoMaster()) {
    tledBVHTraverserGPU_kernels::AddNarrowPhaseNodeFacetTests<SlaveMesh::Facet::NumberOfVertices, blockSize, true> <<<numBlks, blockSize, 0, stream>>> (dp_out, dp_counter, this->GetSlaveMesh().GetOnDeviceFacetVertexIndices(), inItems, dpc_numItems);    
  } else {
    tledBVHTraverserGPU_kernels::AddNarrowPhaseNodeFacetTests<MasterMesh::Facet::NumberOfVertices, blockSize, false> <<<numBlks, blockSize, 0, stream>>> (dp_out, dp_counter, this->GetMasterMesh().GetOnDeviceFacetVertexIndices(), inItems, dpc_numItems);    
  }
}

template <class TMasterBVH, class TSlaveBVH, class TAPI>
void tledBVHTraverserImplGPU<TMasterBVH, TSlaveBVH, TAPI>::InvokeEdgeEdgeBroadPhaseQueueKernel(int2 *dp_out, int *dp_counter, cudaStream_t stream, const int2 *dpc_inItems, const int *dpc_numItems, const int maxNumItems) {
  const int blockSize = 128;
  const int numBlks = tledCUDAHelpers::GetNumberOfBlocks(this->DoMaster()? maxNumItems*SlaveMesh::Facet::NumberOfVertices : maxNumItems*MasterMesh::Facet::NumberOfVertices, blockSize);

  tledLogDebugStream(tledHelper::Info() << "Post-processing <= "<< maxNumItems << " edge-edge broad-phase results in " << (this->DoMaster()? "master" : "slave") << " mode");
  assert(maxNumItems >= this->RetrieveCounterValue(dpc_numItems));
  if (this->DoMaster()) {
    tledBVHTraverserGPU_kernels::AddNarrowPhaseEdgeEdgeTests<SlaveMesh::Facet::NumberOfVertices, blockSize, true> <<<numBlks, blockSize, 0, stream>>> (dp_out, dp_counter, this->GetSlaveMesh().GetOnDeviceFacetEdgeIndices(), this->GetMasterMesh().GetOnDeviceFacetEdgeIndices(), dpc_inItems, dpc_numItems);    
  } else {
    tledBVHTraverserGPU_kernels::AddNarrowPhaseEdgeEdgeTests<MasterMesh::Facet::NumberOfVertices, blockSize, false> <<<numBlks, blockSize, 0, stream>>> (dp_out, dp_counter, this->GetMasterMesh().GetOnDeviceFacetEdgeIndices(), this->GetSlaveMesh().GetOnDeviceFacetEdgeIndices(), dpc_inItems, dpc_numItems);    
  }
}

template <class TMasterBVH, class TSlaveBVH, class TAPI>
tledBVHTraverserImplGPU<TMasterBVH, TSlaveBVH, TAPI>::tledBVHTraverserImplGPU(TSlaveBVH &r_slaveBVH, const TMasterBVH &masterBVH) : Superclass(r_slaveBVH, masterBVH) {
  m_LastMaxBPBufferCapacity = 0;
  mdp_BPSearchPairs = NULL;
  m_NumSearchPairs = m_NumSearchLevels = 0;

  tledCUDAHelpers::AllocateHostMemory(mhp_NumResults);
  tledCUDAHelpers::AllocateDeviceMemory(mdp_NumNodeFacetItems, 2);
  tledCUDAHelpers::AllocateDeviceMemory(mdp_NumEdgeEdgeItems, 2);
  tledCUDAHelpers::AllocateDeviceMemory(mdpp_BPPairCounters[0]);
  tledCUDAHelpers::AllocateDeviceMemory(mdpp_BPPairCounters[1]);
}

template <class TMasterBVH, class TSlaveBVH, class TAPI>
tledBVHTraverserImplGPU<TMasterBVH, TSlaveBVH, TAPI>::~tledBVHTraverserImplGPU() {
  this->ResetResults();

  tledCheckCUDAErrors(cudaFree(mdp_NumNodeFacetItems));
  tledCheckCUDAErrors(cudaFree(mdp_NumEdgeEdgeItems));
  tledCheckCUDAErrors(cudaFreeHost(mhp_NumResults));

  tledCheckCUDAErrors(cudaFree(mdpp_BPPairCounters[0]));
  tledCheckCUDAErrors(cudaFree(mdpp_BPPairCounters[1]));

  if (mdp_BPSearchPairs != NULL) {
    tledCheckCUDAErrors(cudaFree(mdp_BPSearchPairs));
  }
}

template <class TMasterBVH, class TSlaveBVH, class TAPI>
int tledBVHTraverserImplGPU<TMasterBVH, TSlaveBVH, TAPI>::RetrieveCounterValue(const int *dpc_counter) {
  tledCUDAHelpers::CopyFromDevice(this->GetHostItemCounter(), dpc_counter);

  return *this->GetHostItemCounter();
}

template <class TMasterBVH, class TSlaveBVH, class TAPI>
int tledBVHTraverserImplGPU<TMasterBVH, TSlaveBVH, TAPI>::RetrieveCounterValue(const int *dpc_counter, const cudaStream_t stream) {
  tledCUDAHelpers::CopyFromDevice(this->GetHostItemCounter(), dpc_counter, 1, stream);
  tledCheckCUDAErrors(cudaStreamSynchronize(stream));

  return *this->GetHostItemCounter();
}

template <class TMasterBVH, class TSlaveBVH, class TAPI>
std::vector<int2> tledBVHTraverserImplGPU<TMasterBVH, TSlaveBVH, TAPI>::GenerateBroadPhaseCheckPairs(const std::vector<int> &slaveLeaves, const int masterRootInd) {
  std::vector<int2> copyBuffer;
  
  copyBuffer.reserve(slaveLeaves.size());
  for (std::vector<int>::const_iterator ic_s = slaveLeaves.begin(); ic_s < slaveLeaves.end(); ic_s++) {
    int2 item;

    item.x = masterRootInd;
    item.y = *ic_s;
    copyBuffer.push_back(item);
  }

  return copyBuffer;
}

template <class TMasterBVH, class TSlaveBVH, class TAPI>
void tledBVHTraverserImplGPU<TMasterBVH, TSlaveBVH, TAPI>::PrepareNarrowPhase(const int2 *dpc_inBPList, const int *dpc_bpItemCounter, const int maxBPResultItems) {
  const int maxNumNodeFacet = this->DoMaster()? SlaveMesh::Facet::NumberOfVertices*maxBPResultItems : MasterMesh::Facet::NumberOfVertices*maxBPResultItems;
  const int maxNumEdgeEdge = SlaveMesh::Facet::NumberOfVertices*MasterMesh::Facet::NumberOfVertices*maxBPResultItems;
  
  cudaStream_t nfStream, eeStream;
  tledCUDADeviceMemoryBlock &r_nfBuf = tledCUDADeviceMemoryBlock::GetNextFreeBufferWithSize<int2>(maxNumNodeFacet);
  tledCUDADeviceMemoryBlock &r_eeBuf = tledCUDADeviceMemoryBlock::GetNextFreeBufferWithSize<int2>(maxNumEdgeEdge);

  tledCheckCUDAErrors(cudaStreamCreate(&nfStream));
  tledCheckCUDAErrors(cudaStreamCreate(&eeStream));
  tledLogDebugStream(tledHelper::Info() << "Post-processing " << maxBPResultItems << " items in broad-phase output queue.");

  tledLogDebugStream(tledHelper::Info() << "Starting node-facet post-processing on <= " << maxBPResultItems << " unprocessed BP items");
  tledCUDAHelpers::SetZero(this->GetOnDeviceNodeFacetItemCounter(), 1, nfStream);
  this->InvokeNodeFacetBroadPhaseQueueKernel(r_nfBuf.GetBuffer<int2>(), this->GetOnDeviceNodeFacetItemCounter(), nfStream, dpc_inBPList, dpc_bpItemCounter, maxBPResultItems);  
  assert(this->RetrieveCounterValue(this->GetOnDeviceNodeFacetItemCounter(), nfStream) <= r_nfBuf.GetMaxNumberOfItems<int2>());
  tledCUDAHelpers::SetZero(this->GetOnDeviceEdgeEdgeItemCounter(), 1, eeStream);
  tledLogDebugStream(tledHelper::Info() << "Starting edge-edge post-processing on <= " << maxBPResultItems << " unprocessed BP items");
  this->InvokeEdgeEdgeBroadPhaseQueueKernel(r_eeBuf.GetBuffer<int2>(), this->GetOnDeviceEdgeEdgeItemCounter(), eeStream, dpc_inBPList, dpc_bpItemCounter, maxBPResultItems);
  assert(this->RetrieveCounterValue(this->GetOnDeviceEdgeEdgeItemCounter(), eeStream) <= r_eeBuf.GetMaxNumberOfItems<int2>());
  
  this->SetNumberOfNodeFacetContacts(this->RetrieveCounterValue(this->GetOnDeviceNodeFacetItemCounter(), nfStream));
  tledCheckCUDAErrors(cudaStreamDestroy(nfStream));
  this->SetNumberOfEdgeEdgeContacts(this->RetrieveCounterValue(this->GetOnDeviceEdgeEdgeItemCounter(), eeStream));
  tledCheckCUDAErrors(cudaStreamDestroy(eeStream));
  
  this->SetEdgeEdgeResultBuffer(r_eeBuf);
  this->SetNodeFacetResultBuffer(r_nfBuf);
}

template <class TMasterBVH, class TSlaveBVH, class TAPI>
void tledBVHTraverserImplGPU<TMasterBVH, TSlaveBVH, TAPI>::FindCollisions() {
  this->ResetResults();
  Superclass::FindCollisions();
}

template <class TMasterBVH, class TSlaveBVH, class TAPI>
void tledBVHTraverserImplGPU<TMasterBVH, TSlaveBVH, TAPI>::RunBroadPhase() {
  const int bvhOrder = TMasterBVH::BVHOrder;

  cudaStream_t bpStream, readBackStream;
  tledCUDADeviceMemoryBlock *pp_searchPairBuffers[2];
  int currentNumPairs;

  pp_searchPairBuffers[0] = &tledCUDADeviceMemoryBlock::GetNextFreeBufferWithSize<int2>(m_LastMaxBPBufferCapacity = std::max(m_LastMaxBPBufferCapacity, std::max(128, this->GetNumberOfInitialBroadPhaseSearchPairs()*bvhOrder)));
  pp_searchPairBuffers[1] = &tledCUDADeviceMemoryBlock::GetNextFreeBufferWithSize<int2>(m_LastMaxBPBufferCapacity);

  this->InitialiseBroadPhaseSearch(pp_searchPairBuffers[1]->GetBuffer<int2>(), mdpp_BPPairCounters[1]);  
  assert(this->GetNumberOfInitialBroadPhaseSearchPairs() == this->RetrieveCounterValue(mdpp_BPPairCounters[1]));
  
  tledCUDAHelpers::SetZero(mdpp_BPPairCounters[0]);
  this->FilterNonIntersectingBroadPhaseBVPairs(pp_searchPairBuffers[0]->GetBuffer<int2>(), mdpp_BPPairCounters[0], pp_searchPairBuffers[1]->GetBuffer<int2>(), mdpp_BPPairCounters[1], this->GetNumberOfInitialBroadPhaseSearchPairs(), 0);
  currentNumPairs = this->GetNumberOfInitialBroadPhaseSearchPairs();

  tledCheckCUDAErrors(cudaStreamCreate(&bpStream));
  tledCheckCUDAErrors(cudaStreamCreate(&readBackStream));
  for (int l = 1; l <= this->GetNumberOfSearchBVHLevels(); l++) {    
    assert(this->RetrieveCounterValue(mdpp_BPPairCounters[0])*bvhOrder <= m_LastMaxBPBufferCapacity);
    assert(pp_searchPairBuffers[0]->IsActive() && pp_searchPairBuffers[1]->IsActive());
#ifndef NDEBUG
    tledCheckCUDAErrors(cudaMemset(pp_searchPairBuffers[1]->GetBuffer<int2>(), -1, sizeof(int2)*m_LastMaxBPBufferCapacity));
#endif

    tledCUDAHelpers::SetZero(mdpp_BPPairCounters[1], 1, bpStream);
    this->DescendInMasterBVH(pp_searchPairBuffers[1]->GetBuffer<int2>(), mdpp_BPPairCounters[1], pp_searchPairBuffers[0]->GetBuffer<int2>(), mdpp_BPPairCounters[0], currentNumPairs, bpStream);
    currentNumPairs = this->RetrieveCounterValue(mdpp_BPPairCounters[0], readBackStream)*bvhOrder;
    tledCheckCUDAErrors(cudaStreamSynchronize(bpStream));

#ifndef NDEBUG
    tledCheckCUDAErrors(cudaMemset(pp_searchPairBuffers[0]->GetBuffer<int2>(), -1, sizeof(int2)*m_LastMaxBPBufferCapacity));
#endif
    tledCUDAHelpers::SetZero(mdpp_BPPairCounters[0], 1);
    this->FilterNonIntersectingBroadPhaseBVPairs(pp_searchPairBuffers[0]->GetBuffer<int2>(), mdpp_BPPairCounters[0], pp_searchPairBuffers[1]->GetBuffer<int2>(), mdpp_BPPairCounters[1], currentNumPairs, 0);

    if (currentNumPairs == 0) {
      tledLogDebugStream(tledHelper::Info() << "Early broad-phase exit after descending " << l << " levels in master BVH, no more candidate pairs left..");
      break;
    }

    if (currentNumPairs*bvhOrder > m_LastMaxBPBufferCapacity) {      
      m_LastMaxBPBufferCapacity = int(1.2f*currentNumPairs*bvhOrder);
      tledLogDebugStream(tledHelper::Info() << "Increasing broad-phase buffer capacity to " << m_LastMaxBPBufferCapacity);
      for (int b = 0; b < 2; b++) pp_searchPairBuffers[b]->Grow<int2>(m_LastMaxBPBufferCapacity);	
    }    
  }
  tledLogDebugStream(tledHelper::Info() << "Exiting broad-phase search with " << currentNumPairs << " BV-pairs remaining in list.");
  tledCheckCUDAErrors(cudaStreamDestroy(bpStream));
  tledCheckCUDAErrors(cudaStreamDestroy(readBackStream));

  pp_searchPairBuffers[1]->ToggleActive();
  assert(!pp_searchPairBuffers[1]->IsActive());
  pp_searchPairBuffers[1] = NULL;

  assert(this->GetNumberOfNodeFacetContacts() == 0 && this->GetNumberOfEdgeEdgeContacts() == 0);
  if (currentNumPairs > 0) {
    this->ConvertBVToPrimitivePairs(pp_searchPairBuffers[0]->GetBuffer<int2>(), mdpp_BPPairCounters[0], currentNumPairs);
    this->PrepareNarrowPhase(pp_searchPairBuffers[0]->GetBuffer<int2>(), mdpp_BPPairCounters[0], currentNumPairs);
  } 
  pp_searchPairBuffers[0]->ToggleActive();
  assert(!pp_searchPairBuffers[0]->IsActive());
}

template <class TMasterBVH, class TSlaveBVH, class TAPI>
void tledBVHTraverserImplGPU<TMasterBVH, TSlaveBVH, TAPI>::ConvertBVToPrimitivePairs(int2 *dp_list, const int *dpc_inCounter, const int maxNumItems) {
  const int blockSize = 256;
  const int numBlks = tledCUDAHelpers::GetNumberOfBlocks(maxNumItems, blockSize);

  tledLogDebugStream(tledHelper::Info() << "Invoking BV-list conversion with " << numBlks << " blocks");
  tledBVHTraverserGPU_kernels::ConvertBVToPrimitivePairsKernel <<<numBlks, blockSize>>> (dp_list, dpc_inCounter, this->GetMasterBVH().GetOnDeviceBVs(), this->GetSlaveBVH().GetOnDeviceBVs());
  tledDeviceSyncDebug;
}

template <class TMasterBVH, class TSlaveBVH, class TAPI>
void tledBVHTraverserImplGPU<TMasterBVH, TSlaveBVH, TAPI>::FilterNonIntersectingBroadPhaseBVPairs(int2 *dp_outList, int *dp_outCounter, const int2 *dpc_inList, const int *dpc_inCounter, const int maxNumItems, const cudaStream_t stream) {
  const int blockSize = 128;
  const int numBlks = tledCUDAHelpers::GetNumberOfBlocks(maxNumItems, blockSize);

  tledLogDebugStream(tledHelper::Info() << "Invoking broad-phase search filter kernel for non-intersecting BV-pairs with " << numBlks << " blocks");
  tledBVHTraverserGPU_kernels::FilterNonIntersectingBroadPhaseBVPairsKernel<blockSize> <<<numBlks, blockSize, 0, stream>>> (dp_outList, dp_outCounter, dpc_inList, dpc_inCounter, this->GetMasterBVH().GetOnDeviceBVs(), this->GetSlaveBVH().GetOnDeviceBVs());
  tledDeviceSyncDebug;
}

template <class TMasterBVH, class TSlaveBVH, class TAPI>
void tledBVHTraverserImplGPU<TMasterBVH, TSlaveBVH, TAPI>::DescendInMasterBVH(int2 *dp_outList, int *dp_outCounter, const int2 *dpc_inList, const int *dpc_inCounter, const int maxNumItems, const cudaStream_t stream) {
  const int blockSize = 128;
  const int numBlks = tledCUDAHelpers::GetNumberOfBlocks(maxNumItems, blockSize);

  tledLogDebugStream(tledHelper::Info() << "Descending in master BVH, updating search BV-pairs with " << numBlks << " blocks");
  tledBVHTraverserGPU_kernels::DescendInMasterBVHKernel<blockSize> <<<numBlks, blockSize, 0, stream>>> (dp_outList, dp_outCounter, dpc_inList, dpc_inCounter, this->GetMasterBVH().GetOnDeviceBVs());
  tledDeviceSyncDebug;
}

template <class TMasterBVH, class TSlaveBVH, class TAPI>
void tledBVHTraverserImplGPU<TMasterBVH, TSlaveBVH, TAPI>::Init(tledUnstructuredContactManager &r_manager) {
  this->ResetResults();
  this->InitStartBVPairs();
}

template <class TMasterBVH, class TSlaveBVH, class TAPI>
void tledBVHTraverserImplGPU<TMasterBVH, TSlaveBVH, TAPI>::InitialiseBroadPhaseSearch(int2 *dp_dst, int *dp_counter) {
  cudaStream_t devCopyStream, hostCopyStream;

  tledLogDebugStream(tledHelper::Info() << "Building initial list of " << m_NumSearchPairs << " BV-pairs for broad-phase search.");
  tledCheckCUDAErrors(cudaStreamCreate(&devCopyStream));
  tledCheckCUDAErrors(cudaMemcpyAsync(dp_dst, mdp_BPSearchPairs, m_NumSearchPairs*sizeof(int2), cudaMemcpyDeviceToDevice, devCopyStream));
  tledCheckCUDAErrors(cudaStreamCreate(&hostCopyStream));
  tledCUDAHelpers::CopyToDevice(dp_counter, &m_NumSearchPairs, 1, hostCopyStream);
  tledCheckCUDAErrors(cudaStreamSynchronize(devCopyStream));
  tledCheckCUDAErrors(cudaStreamDestroy(devCopyStream));
  tledCheckCUDAErrors(cudaStreamSynchronize(hostCopyStream));
  tledCheckCUDAErrors(cudaStreamDestroy(hostCopyStream));
}

template <class TMasterBVH, class TSlaveBVH, class TAPI>
void tledBVHTraverserImplGPU<TMasterBVH, TSlaveBVH, TAPI>::SetBroadPhaseSearchSettings(const int2 *searchBVPairs, const int numSearchBVPairs, const int numSearchLevels) {
  if (mdp_BPSearchPairs != NULL) {
    tledCheckCUDAErrors(cudaFree(mdp_BPSearchPairs));
  }

  tledCUDAHelpers::AllocateDeviceMemory(mdp_BPSearchPairs, numSearchBVPairs);
  tledCUDAHelpers::CopyToDevice(mdp_BPSearchPairs, searchBVPairs, numSearchBVPairs);

  m_NumSearchPairs = numSearchBVPairs;
  m_NumSearchLevels = numSearchLevels;
}

template <class TMasterBVH, class TSlaveBVH, class TAPI>
void tledBVHTraverserImplGPU<TMasterBVH, TSlaveBVH, TAPI>::RunNarrowPhase() {
  tledDeviceSyncDebug;
  tledLogDebugStream(tledHelper::Info() << "Starting narrow-phase with " <<  this->GetNumberOfNodeFacetContacts() << " node-facet and " <<  this->GetNumberOfEdgeEdgeContacts() << " edge-edge pairs.");
  if (this->GetNumberOfNodeFacetContacts() > 0 || this->GetNumberOfEdgeEdgeContacts() > 0) {
    typedef typename MasterMesh::GPUSurface __GPUSurface;

    tledCUDADeviceMemoryBlock *p_nfStage1Buffer = NULL, *p_nfStage2Buffer = NULL, *p_eeStage1Buffer = NULL, *p_eeStage2Buffer = NULL;
    int numEdgeEdge = this->GetNumberOfEdgeEdgeContacts(), numNodeFacet = this->GetNumberOfNodeFacetContacts();
    cudaStream_t nfStream, eeStream;

    tledCheckCUDAErrors(cudaStreamCreate(&nfStream));
    tledCheckCUDAErrors(cudaStreamCreate(&eeStream));
    this->GetSlaveMesh().UpdateAllNormals();

    if (numNodeFacet > 0) {
      const int2 *items = this->GetNodeFacetResults().template GetBuffer<int2>();

      tledCUDAHelpers::SetZero(this->GetOnDeviceNodeFacetItemCounter(), 1, nfStream);      
      p_nfStage1Buffer = this->InvokeNodeFacetNarrowPhaseStage1(this->GetOnDeviceNodeFacetItemCounter(), nfStream, items, numNodeFacet);
    }

    if (numEdgeEdge > 0) {      
      const int2 *items = this->GetEdgeEdgeResults().template GetBuffer<int2>();

      tledCUDAHelpers::SetZero(this->GetOnDeviceEdgeEdgeItemCounter(), 1, eeStream);      
      p_eeStage1Buffer = this->InvokeEdgeEdgeNarrowPhaseStage1(this->GetOnDeviceEdgeEdgeItemCounter(), eeStream, items, numEdgeEdge);
    }

    tledDeviceSyncDebug;
    if (numNodeFacet > 0) {
      numNodeFacet = this->RetrieveCounterValue(this->GetOnDeviceNodeFacetItemCounter(), nfStream);
      this->GetNodeFacetResults().ToggleActive();
    }
    assert(!this->GetNodeFacetResults().IsActive());

    if (numEdgeEdge > 0) {
      numEdgeEdge = this->RetrieveCounterValue(this->GetOnDeviceEdgeEdgeItemCounter(), eeStream);
      this->GetEdgeEdgeResults().ToggleActive();
    }
    assert(!this->GetEdgeEdgeResults().IsActive());

    if (numNodeFacet > 0) {
      typedef tledBVHTraverserGPU_kernels::NodeFacetProjectionOrdering<NodeFacetNarrowPhaseResult> __Ordering;

      tledCUDAHelpers::SetZero(this->GetOnDeviceNodeFacetItemCounter(), 1, nfStream);      
      tledLogDebugStream(tledHelper::Info() << "Launching 2nd stage node-facet narrow-phase kernel with " << numNodeFacet << " input items.");
      p_nfStage2Buffer = this->InvokeNodeFacetNarrowPhaseStage2(this->GetOnDeviceNodeFacetItemCounter(), nfStream, p_nfStage1Buffer->GetBuffer<NodeFacetInitialProjection>(), numNodeFacet);
      tledComplexCUDAHelpers::MakeSortedUnique<NodeFacetNarrowPhaseResult, __Ordering>(*p_nfStage2Buffer, this->GetOnDeviceNodeFacetItemCounter(), __Ordering(), nfStream);
    } else if (p_nfStage1Buffer != NULL) {
      p_nfStage1Buffer->ToggleActive();
      assert(!p_nfStage1Buffer->IsActive());
    }

    if (numEdgeEdge > 0) {
      typedef tledBVHTraverserGPU_kernels::EdgeEdgeProjectionOrdering<EdgeEdgeNarrowPhaseResult> __Ordering;

      tledCUDAHelpers::SetZero(this->GetOnDeviceEdgeEdgeItemCounter(), 1, eeStream);
      tledLogDebugStream(tledHelper::Info() << "Launching 2nd stage edge-edge narrow-phase kernel with " << numEdgeEdge << " input items.");
      p_eeStage2Buffer = this->InvokeEdgeEdgeNarrowPhaseStage2(this->GetOnDeviceEdgeEdgeItemCounter(), eeStream, p_eeStage1Buffer->GetBuffer<EdgeEdgeInitialProjection>(), numEdgeEdge);
      tledComplexCUDAHelpers::MakeSortedUnique<EdgeEdgeNarrowPhaseResult, __Ordering>(*p_eeStage2Buffer, this->GetOnDeviceEdgeEdgeItemCounter(), __Ordering(), eeStream);
    } else if (p_eeStage1Buffer != NULL) {
      p_eeStage1Buffer->ToggleActive();
      assert(!p_eeStage1Buffer->IsActive());
    }

    if (numNodeFacet > 0) {
      assert(p_nfStage2Buffer->IsActive());

      numNodeFacet = this->RetrieveCounterValue(this->GetOnDeviceNodeFacetItemCounter(), nfStream);
      p_nfStage1Buffer->ToggleActive();      	

      if (numNodeFacet > 0) this->SetNodeFacetResultBuffer(*p_nfStage2Buffer);      
      else p_nfStage2Buffer->ToggleActive();

      tledLogDebugStream(tledHelper::Info() << "Have " << numNodeFacet << " unique final node-facet projections.");
      assert(!p_nfStage1Buffer->IsActive());
    }

    if (numEdgeEdge > 0) {
      assert(p_eeStage2Buffer->IsActive());

      numEdgeEdge = this->RetrieveCounterValue(this->GetOnDeviceEdgeEdgeItemCounter(), eeStream);
      p_eeStage1Buffer->ToggleActive();

      if (numEdgeEdge > 0) this->SetEdgeEdgeResultBuffer(*p_eeStage2Buffer);
      else p_eeStage2Buffer->ToggleActive();

      tledLogDebugStream(tledHelper::Info() << "Have " << numEdgeEdge << " unique final edge-edge projections.");
    } 
    assert((numEdgeEdge == 0 && (p_eeStage1Buffer == NULL || !p_eeStage1Buffer->IsActive())) || (numEdgeEdge > 0 && p_eeStage2Buffer->IsActive()));

    tledCheckCUDAErrors(cudaStreamDestroy(nfStream));
    tledCheckCUDAErrors(cudaStreamDestroy(eeStream));

    this->SetNumberOfEdgeEdgeContacts(numEdgeEdge);
    this->SetNumberOfNodeFacetContacts(numNodeFacet);
  } else {
    this->ResetResults();
  } /* if have narrow-phase items else .. */

  tledLogDebugStream(tledHelper::Info() << "On traverser exit: " << this->GetNumberOfNodeFacetContacts() << " node-facet contacts, " << this->GetNumberOfEdgeEdgeContacts() << " edge-edge contacts; "
		     << tledCUDADeviceMemoryBlock::GetNumberOfActiveBuffers() << " active device memory blocks");
}
