// =========================================================================
// File:       tledTrackingBVHTraverserGPU.cu
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   CUDA
// Created:    December 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledTrackingBVHTraverserGPU_CU
#define tledTrackingBVHTraverserGPU_CU

#include "tledTrackingBVHTraverserGPU.h"
#include "tledCUDAHelpers.h"

#include "tledTrackingBVHTraverserGPU_kernels.cu"

template <class TBaseTraverserGPU>
void tledTrackingBVHTraverserImplGPU<TBaseTraverserGPU>::ExtractTrackingInformationFromContacts() {
  cudaStream_t nfStream, eeStream;
  tledCUDADeviceMemoryBlock *p_edgeEdgeBuffer = NULL, *p_nodeFacetBuffer = NULL;
  
  tledDeviceSyncDebug;
  tledCheckCUDAErrors(cudaStreamCreate(&nfStream));
  tledCheckCUDAErrors(cudaStreamCreate(&eeStream));
  
  tledLogDebugStream(tledHelper::Info() << "Extracting tracking information (" << (this->DoMaster()? "master" : "slave") << ") from " << this->GetNumberOfNodeFacetContacts() << " node-facet, " << this->GetNumberOfEdgeEdgeContacts() << " edge-edge contacts.");
  if (this->DoMaster()) {
    if (this->GetNumberOfNodeFacetContacts() > 0) {
      assert(this->GetMasterNodeFacetTrackingBuffer() == NULL);
      _UpdateTrackedNodeFacetContacts<MasterMesh::Facet::NumberOfVertices>(p_nodeFacetBuffer, this->GetOnDeviceNodeFacetItemCounter(), nfStream, this->GetNodeFacetResults(), this->GetNumberOfNodeFacetContacts(), this->GetMasterNodeFacetNeighbours(), this->GetMasterNodeFacetNeighbourRanges(), this->GetMaxMasterNodeFacetNeighbourRangeSize());
      this->SetMasterNodeFacetTrackingBuffer(p_nodeFacetBuffer);
    }

    if (this->GetNumberOfEdgeEdgeContacts() > 0) {
      assert(this->GetMasterEdgeEdgeTrackingBuffer() == NULL);
      _UpdateTrackedEdgeEdgeContacts(p_edgeEdgeBuffer, this->GetOnDeviceEdgeEdgeItemCounter(), eeStream, this->GetEdgeEdgeResults(), this->GetNumberOfEdgeEdgeContacts(), this->GetMasterNodeEdgeNeighbours(), this->GetMasterNodeEdgeNeighbourRanges(), this->GetSlaveNodeEdgeNeighbours(), this->GetSlaveNodeEdgeNeighbourRanges(), this->GetMaxSlaveNodeEdgeNeighbourRangeSize(), this->GetMaxMasterNodeEdgeNeighbourRangeSize());
      this->SetMasterEdgeEdgeTrackingBuffer(p_edgeEdgeBuffer);      
    }

    if (this->GetNumberOfNodeFacetContacts() > 0) m_NumberOfTrackedMasterNodeFacetContacts = this->RetrieveCounterValue(this->GetOnDeviceNodeFacetItemCounter(), nfStream);
    else m_NumberOfTrackedMasterNodeFacetContacts = 0;

    if (this->GetNumberOfEdgeEdgeContacts() > 0) m_NumberOfTrackedMasterEdgeEdgeContacts = this->RetrieveCounterValue(this->GetOnDeviceEdgeEdgeItemCounter(), eeStream);
    else m_NumberOfTrackedMasterEdgeEdgeContacts = 0;
  } else {
    if (this->GetNumberOfNodeFacetContacts() > 0) {
      assert(this->GetSlaveNodeFacetTrackingBuffer() == NULL);
      _UpdateTrackedNodeFacetContacts<SlaveMesh::Facet::NumberOfVertices>(p_nodeFacetBuffer, this->GetOnDeviceNodeFacetItemCounter(), nfStream, this->GetNodeFacetResults(), this->GetNumberOfNodeFacetContacts(), this->GetSlaveNodeFacetNeighbours(), this->GetSlaveNodeFacetNeighbourRanges(), this->GetMaxSlaveNodeFacetNeighbourRangeSize());
      this->SetSlaveNodeFacetTrackingBuffer(p_nodeFacetBuffer);
    }

    if (this->GetNumberOfEdgeEdgeContacts() > 0) {
      assert(this->GetSlaveEdgeEdgeTrackingBuffer() == NULL);
      _UpdateTrackedEdgeEdgeContacts(p_edgeEdgeBuffer, this->GetOnDeviceEdgeEdgeItemCounter(), eeStream, this->GetEdgeEdgeResults(), this->GetNumberOfEdgeEdgeContacts(), this->GetSlaveNodeEdgeNeighbours(), this->GetSlaveNodeEdgeNeighbourRanges(), this->GetMasterNodeEdgeNeighbours(), this->GetMasterNodeEdgeNeighbourRanges(), this->GetMaxMasterNodeEdgeNeighbourRangeSize(), this->GetMaxSlaveNodeEdgeNeighbourRangeSize());
      this->SetSlaveEdgeEdgeTrackingBuffer(p_edgeEdgeBuffer);
    }

    if (this->GetNumberOfNodeFacetContacts() > 0) m_NumberOfTrackedSlaveNodeFacetContacts = this->RetrieveCounterValue(this->GetOnDeviceNodeFacetItemCounter(), nfStream);
    else m_NumberOfTrackedSlaveNodeFacetContacts = 0;

    if (this->GetNumberOfEdgeEdgeContacts() > 0) m_NumberOfTrackedSlaveEdgeEdgeContacts = this->RetrieveCounterValue(this->GetOnDeviceEdgeEdgeItemCounter(), eeStream);
    else m_NumberOfTrackedSlaveEdgeEdgeContacts = 0;
  }

  tledCheckCUDAErrors(cudaStreamDestroy(nfStream));
  tledCheckCUDAErrors(cudaStreamDestroy(eeStream));
  tledLogDebugStream(tledHelper::Info() << "Exiting contact tracking with " << (this->DoMaster()? m_NumberOfTrackedMasterNodeFacetContacts : m_NumberOfTrackedSlaveNodeFacetContacts) << "node-facet, "
		     << (this->DoMaster()? m_NumberOfTrackedMasterEdgeEdgeContacts : m_NumberOfTrackedSlaveEdgeEdgeContacts) << " edge-edge contact candidates found.");
}

template <class TBaseTraverserGPU>
void tledTrackingBVHTraverserImplGPU<TBaseTraverserGPU>::PreFullBVHTraversalHook() {
  Superclass::PreFullBVHTraversalHook();
  if (this->DoMaster()) {
    if (this->GetMasterNodeFacetTrackingBuffer() != NULL) {
      assert(this->GetMasterNodeFacetTrackingBuffer()->IsActive());
      this->GetMasterNodeFacetTrackingBuffer()->ToggleActive();
      this->SetMasterNodeFacetTrackingBuffer(NULL);
    }

    if (this->GetMasterEdgeEdgeTrackingBuffer() != NULL) {
      assert(this->GetMasterEdgeEdgeTrackingBuffer()->IsActive());
      this->GetMasterEdgeEdgeTrackingBuffer()->ToggleActive();
      this->SetMasterEdgeEdgeTrackingBuffer(NULL);
    }    
  } else {
    if (this->GetSlaveNodeFacetTrackingBuffer() != NULL) {
      assert(this->GetSlaveNodeFacetTrackingBuffer()->IsActive());
      this->GetSlaveNodeFacetTrackingBuffer()->ToggleActive();
      this->SetSlaveNodeFacetTrackingBuffer(NULL);
    }

    if (this->GetSlaveEdgeEdgeTrackingBuffer() != NULL) {
      assert(this->GetSlaveEdgeEdgeTrackingBuffer()->IsActive());
      this->GetSlaveEdgeEdgeTrackingBuffer()->ToggleActive();
      this->SetSlaveEdgeEdgeTrackingBuffer(NULL);
    }    
  }
}

template <class TBaseTraverserGPU>
void tledTrackingBVHTraverserImplGPU<TBaseTraverserGPU>::ConvertHostToDeviceNeighbourLists(int* &rdp_list, int2* &rdp_ranges, const std::vector<int> &list, const std::vector<std::pair<int, int> > &ranges) {
  std::vector<int2> hostRanges;

  tledCUDAHelpers::AllocateDeviceMemory(rdp_list, list.size());
  tledCUDAHelpers::CopyToDevice(rdp_list, list);
  hostRanges.resize(ranges.size());
  for (size_t i = 0; i < ranges.size(); i++) {
    hostRanges[i].x = ranges[i].first;
    hostRanges[i].y = ranges[i].second;
  }
  tledCUDAHelpers::AllocateDeviceMemory(rdp_ranges, hostRanges.size());
  tledCUDAHelpers::CopyToDevice(rdp_ranges, hostRanges);
}

template <class TBaseTraverserGPU>
template <const int t_numFacetVertices>
void tledTrackingBVHTraverserImplGPU<TBaseTraverserGPU>::_UpdateTrackedNodeFacetContacts(tledCUDADeviceMemoryBlock* &rp_nodeFacetCandidatePairs, int *dp_counter, cudaStream_t stream, const tledCUDADeviceMemoryBlock &nodeFacetResultBuffer, const int numContacts, const int *dpc_masterNodeFacetNeighbours, const int2 *dpc_masterNodeFacetNeighbourRanges, const int maxMasterRangeSize) {
  const int blockSize = 256;
  const int maxNumOutItems = maxMasterRangeSize*t_numFacetVertices*numContacts;
  const int numBlocks = tledCUDAHelpers::GetNumberOfBlocks(maxNumOutItems, blockSize);
  const NodeFacetNarrowPhaseResult *results = nodeFacetResultBuffer.GetBuffer<NodeFacetNarrowPhaseResult>();

  int2 *dp_dst;

  assert(nodeFacetResultBuffer.IsActive());
  tledCUDAHelpers::SetZero(dp_counter, 1, stream);
  rp_nodeFacetCandidatePairs = &tledCUDADeviceMemoryBlock::GetNextFreeBufferWithSize<int2>(maxNumOutItems);
  dp_dst = rp_nodeFacetCandidatePairs->GetBuffer<int2>();
  tledLogDebugStream(tledHelper::Info() << "Updating node-facet tracking data from " << numContacts << " previously detected contacts, yielding up to " << maxNumOutItems << " candidates.");
  tledDeviceSyncDebug;
  tledTrackingBVHTraverserGPU_kernels::UpdateTrackedNodeFacetContacts<NodeFacetNarrowPhaseResult, blockSize, t_numFacetVertices> <<< numBlocks, blockSize, 0, stream >>> (dp_dst, dp_counter, results, numContacts, dpc_masterNodeFacetNeighbours, dpc_masterNodeFacetNeighbourRanges);
  tledDeviceSyncDebug;
}

template <class TBaseTraverserGPU>
void tledTrackingBVHTraverserImplGPU<TBaseTraverserGPU>::_UpdateTrackedEdgeEdgeContacts(tledCUDADeviceMemoryBlock* &rp_edgeEdgeCandidatePairs, int *dp_counter, cudaStream_t stream, const tledCUDADeviceMemoryBlock &edgeEdgeResultBuffer, const int numContacts, const int *dpc_masterNodeEdgeNeigbours, const int2 *dpc_masterNodeEdgeNeigbourRanges, const int *dpc_slaveNodeEdgeNeigbours, const int2 *dpc_slaveNodeEdgeNeigbourRanges, const int maxMasterRangeSize, const int maxSlaveRangeSize) {
  const int blockSize = 256;
  const int maxNumOutItems = maxMasterRangeSize*maxSlaveRangeSize*2*2*numContacts;
  const int numBlocks = tledCUDAHelpers::GetNumberOfBlocks(maxNumOutItems, blockSize);
  const EdgeEdgeNarrowPhaseResult *results = edgeEdgeResultBuffer.GetBuffer<EdgeEdgeNarrowPhaseResult>();

  int2 *dp_dst;

  assert(edgeEdgeResultBuffer.IsActive());
  tledCUDAHelpers::SetZero(dp_counter, 1, stream);
  rp_edgeEdgeCandidatePairs = &tledCUDADeviceMemoryBlock::GetNextFreeBufferWithSize<int2>(maxNumOutItems);
  dp_dst = rp_edgeEdgeCandidatePairs->GetBuffer<int2>();
  tledLogDebugStream(tledHelper::Info() << "Updating edge-edge tracking data from " << numContacts << " previously detected contacts.");
  tledTrackingBVHTraverserGPU_kernels::UpdateTrackedEdgeEdgeContacts<EdgeEdgeNarrowPhaseResult, blockSize> <<< numBlocks, blockSize, 0, stream >>> (dp_dst, dp_counter, results, numContacts, dpc_masterNodeEdgeNeigbours, dpc_masterNodeEdgeNeigbourRanges, dpc_slaveNodeEdgeNeigbours, dpc_slaveNodeEdgeNeigbourRanges);
  tledDeviceSyncDebug;
}

template <class TBaseTraverserGPU>
void tledTrackingBVHTraverserImplGPU<TBaseTraverserGPU>::RunTrackingBroadPhase() {
  if (this->DoMaster()) {
    tledLogDebugStream(tledHelper::Info() << "Returning previously computed master contacs: " << this->GetNumberOfTrackedMasterNodeFacetContacts() << " node-facet, " << this->GetNumberOfTrackedMasterEdgeEdgeContacts() << "edge-edge");
    this->SetNodeFacetResultBuffer(*this->GetMasterNodeFacetTrackingBuffer());
    this->SetSlaveNodeFacetTrackingBuffer(NULL);
    this->SetNumberOfNodeFacetContacts(this->GetNumberOfTrackedMasterNodeFacetContacts());

    this->SetEdgeEdgeResultBuffer(*this->GetMasterEdgeEdgeTrackingBuffer());
    this->SetMasterEdgeEdgeTrackingBuffer(NULL);
    this->SetNumberOfEdgeEdgeContacts(this->GetNumberOfTrackedMasterEdgeEdgeContacts());
  } else {
    tledLogDebugStream(tledHelper::Info() << "Returning previously computed slave contacs: " << this->GetNumberOfTrackedSlaveNodeFacetContacts() << " node-facet, " << this->GetNumberOfTrackedSlaveEdgeEdgeContacts() << "edge-edge");
    this->SetNodeFacetResultBuffer(*this->GetSlaveNodeFacetTrackingBuffer());
    this->SetSlaveNodeFacetTrackingBuffer(NULL);
    this->SetNumberOfNodeFacetContacts(this->GetNumberOfTrackedSlaveNodeFacetContacts());

    this->SetEdgeEdgeResultBuffer(*this->GetSlaveEdgeEdgeTrackingBuffer());
    this->SetSlaveEdgeEdgeTrackingBuffer(NULL);
    this->SetNumberOfEdgeEdgeContacts(this->GetNumberOfTrackedSlaveEdgeEdgeContacts());
  }
}

template <class TBaseTraverserGPU>
tledTrackingBVHTraverserImplGPU<TBaseTraverserGPU>::tledTrackingBVHTraverserImplGPU(SlaveBVH &r_slaveBVH, const MasterBVH &masterBVH) : Superclass(r_slaveBVH, masterBVH) {
  mdp_SlaveNodeEdgeNeighbourRanges = mdp_SlaveNodeFacetNeighbourRanges = NULL;
  mdp_SlaveNodeEdgeNeighbours = mdp_SlaveNodeFacetNeighbours = NULL;
  mdp_MasterNodeEdgeNeighbourRanges = mdp_MasterNodeFacetNeighbourRanges = NULL;
  mdp_MasterNodeEdgeNeighbours = mdp_MasterNodeFacetNeighbours = NULL;

  mp_TrackedSlaveNodeFacetContacts = mp_TrackedSlaveEdgeEdgeContacts = NULL;
  mp_TrackedMasterNodeFacetContacts = mp_TrackedMasterEdgeEdgeContacts = NULL;
}

template <class TBaseTraverserGPU>
void tledTrackingBVHTraverserImplGPU<TBaseTraverserGPU>::DeallocateSlaveNeighbourData() {
  if (mdp_SlaveNodeFacetNeighbours != NULL) {
    tledCheckCUDAErrors(cudaFree(mdp_SlaveNodeFacetNeighbours));
    tledCheckCUDAErrors(cudaFree(mdp_SlaveNodeFacetNeighbourRanges));
  }

  if (mdp_SlaveNodeEdgeNeighbours != NULL) {
    tledCheckCUDAErrors(cudaFree(mdp_SlaveNodeEdgeNeighbours));
    tledCheckCUDAErrors(cudaFree(mdp_SlaveNodeEdgeNeighbourRanges));
  }
}

template <class TBaseTraverserGPU>
void tledTrackingBVHTraverserImplGPU<TBaseTraverserGPU>::DeallocateMasterNeighbourData() {
  if (mdp_MasterNodeFacetNeighbours != NULL) {
    tledCheckCUDAErrors(cudaFree(mdp_MasterNodeFacetNeighbours));
    tledCheckCUDAErrors(cudaFree(mdp_MasterNodeFacetNeighbourRanges));
  }

  if (mdp_MasterNodeEdgeNeighbours != NULL) {
    tledCheckCUDAErrors(cudaFree(mdp_MasterNodeEdgeNeighbours));
    tledCheckCUDAErrors(cudaFree(mdp_MasterNodeEdgeNeighbourRanges));
  }
}

template <class TBaseTraverserGPU>
void tledTrackingBVHTraverserImplGPU<TBaseTraverserGPU>::SetSlaveNodeFacetNeighbours(int *dp_neighbours, int2 *dp_ranges, const int maxRangeSize) {
  mdp_SlaveNodeFacetNeighbours = dp_neighbours;
  mdp_SlaveNodeFacetNeighbourRanges = dp_ranges;
  this->SetMaxSlaveNodeFacetNeighbourRangeSize(maxRangeSize);
}

template <class TBaseTraverserGPU>
void tledTrackingBVHTraverserImplGPU<TBaseTraverserGPU>::SetSlaveNodeEdgeNeighbours(int *dp_neighbours, int2 *dp_ranges, const int maxRangeSize) {
  mdp_SlaveNodeEdgeNeighbours = dp_neighbours;
  mdp_SlaveNodeEdgeNeighbourRanges = dp_ranges;
  this->SetMaxSlaveNodeEdgeNeighbourRangeSize(maxRangeSize);
}

template <class TBaseTraverserGPU>
void tledTrackingBVHTraverserImplGPU<TBaseTraverserGPU>::SetMasterNodeFacetNeighbours(int *dp_neighbours, int2 *dp_ranges, const int maxRangeSize) {
  mdp_MasterNodeFacetNeighbours = dp_neighbours;
  mdp_MasterNodeFacetNeighbourRanges = dp_ranges;
  this->SetMaxMasterNodeFacetNeighbourRangeSize(maxRangeSize);
}

template <class TBaseTraverserGPU>
void tledTrackingBVHTraverserImplGPU<TBaseTraverserGPU>::SetMasterNodeEdgeNeighbours(int *dp_neighbours, int2 *dp_ranges, const int maxRangeSize) {
  mdp_MasterNodeEdgeNeighbours = dp_neighbours;
  mdp_MasterNodeEdgeNeighbourRanges = dp_ranges;
  this->SetMaxMasterNodeEdgeNeighbourRangeSize(maxRangeSize);
}

#endif
