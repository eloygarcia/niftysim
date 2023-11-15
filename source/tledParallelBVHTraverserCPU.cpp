// =========================================================================
// File:       tledParallelBVHTraverserCPU.h
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
#include "tledParallelBVHTraverserCPU.h"

#ifdef _USE_BOOST_
std::vector<std::vector<int>*> tledParallelBVHTraverserCPU::svpv_ThreadSlaveNodeIndices(0), tledParallelBVHTraverserCPU::svpv_ThreadSlaveEdgeIndices(0);
std::vector<std::vector<tledBVHTraverserCPU::MasterFacetContactData>*> tledParallelBVHTraverserCPU::svpv_ThreadNodeProjectionBuffer(0);
std::vector<std::vector<tledBVHTraverserCPU::MasterEdgeContactData>*> tledParallelBVHTraverserCPU::svpv_ThreadEdgeProjectionBuffer(0);
std::vector<std::vector<std::pair<int, int> >*> tledParallelBVHTraverserCPU::svpv_ThreadNodeFacetNarrowPhasePairs(0), tledParallelBVHTraverserCPU::svpv_ThreadEdgeEdgeNarrowPhasePairs(0);

template <class TItem>
static std::vector<TItem>& _ResizeThreadBufferIfNecessary(std::vector<std::vector<TItem>*> &rvpv_buffer, const int threadInd, const size_t szBuffer) {
  while (threadInd >= (int)rvpv_buffer.size()) rvpv_buffer.push_back(new std::vector<TItem>(szBuffer));
  if (rvpv_buffer[threadInd]->size() < szBuffer) rvpv_buffer[threadInd]->resize(szBuffer);
  assert((int)rvpv_buffer.size() > threadInd);

  return *rvpv_buffer[threadInd];
}

std::vector<tledBVHTraverserCPU::MasterFacetContactData>& tledParallelBVHTraverserCPU::GetThreadNodeProjectionBuffer(const int threadInd) {
  return _ResizeThreadBufferIfNecessary(svpv_ThreadNodeProjectionBuffer, threadInd, this->GetMaxNumberOfSlaveNodes());
}

std::vector<tledBVHTraverserCPU::MasterEdgeContactData>& tledParallelBVHTraverserCPU::GetThreadEdgeProjectionBuffer(const int threadInd) {
  return _ResizeThreadBufferIfNecessary(svpv_ThreadEdgeProjectionBuffer, threadInd, this->GetMaxNumberOfSlaveEdges());
}

std::vector<int>& tledParallelBVHTraverserCPU::GetThreadContactNodeBuffer(const int threadInd) {
  return _ResizeThreadBufferIfNecessary(svpv_ThreadSlaveNodeIndices, threadInd, 0);
}

std::vector<int>& tledParallelBVHTraverserCPU::GetThreadContactEdgeBuffer(const int threadInd) {
  return _ResizeThreadBufferIfNecessary(svpv_ThreadSlaveEdgeIndices, threadInd, 0);
}

std::vector<std::pair<int, int> >& tledParallelBVHTraverserCPU::GetThreadNodeFacetNarrowPhasePairBuffer(const int threadInd) {
  return _ResizeThreadBufferIfNecessary(svpv_ThreadNodeFacetNarrowPhasePairs, threadInd, 0);
}

std::vector<std::pair<int, int> >& tledParallelBVHTraverserCPU::GetThreadEdgeEdgeNarrowPhasePairBuffer(const int threadInd) {
  return _ResizeThreadBufferIfNecessary(svpv_ThreadEdgeEdgeNarrowPhasePairs, threadInd, 0);
}

template <typename TItem>
static void _DeallocateThreadBuffers(std::vector<std::vector<TItem>*> &rvpv_buffer) {
  for (typename std::vector<std::vector<TItem>*>::iterator ip_tb = rvpv_buffer.begin(); ip_tb < rvpv_buffer.end(); ip_tb++) delete *ip_tb;
  rvpv_buffer.clear();
}

void tledParallelBVHTraverserCPU::DeallocateSharedBuffers() {
  _DeallocateThreadBuffers(svpv_ThreadSlaveNodeIndices);
  _DeallocateThreadBuffers(svpv_ThreadSlaveEdgeIndices);
  _DeallocateThreadBuffers(svpv_ThreadNodeProjectionBuffer);
  _DeallocateThreadBuffers(svpv_ThreadEdgeProjectionBuffer);
  _DeallocateThreadBuffers(svpv_ThreadNodeFacetNarrowPhasePairs);
  _DeallocateThreadBuffers(svpv_ThreadEdgeEdgeNarrowPhasePairs);
}

#endif
