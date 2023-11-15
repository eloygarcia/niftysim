// =========================================================================
// File:       tledBVHTraverserGPU.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    May 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledBVHTraverserGPU_H
#define tledBVHTraverserGPU_H

#include "tledHelper.h"
#include "tledCUDAHelpers.h"
#include "tledBVHTraverser.h"
#include "tledCUDAMemoryBlock.h"

/**
 * \brief Facilities for BVH traversal on the GPU
 * \ingroup contact
 */
class tledBVHTraverserGPU : public tledBVHTraverser {
  /**
   * \name Result Types
   * @{
   */
public:
  /** Node-facet projection base class */
  template <const int t_numberOfFacetVertices>
  struct NodeFacetNarrowPhaseResult {
    int ContactNodeIndices[t_numberOfFacetVertices+1];
    float ShapeValues[t_numberOfFacetVertices];
    float GapValue;
    float3 Normal;
  };

  /** Edge-edge collision */
  struct EdgeEdgeNarrowPhaseResult {
    int2 SlaveEdge, MasterEdge;
    float3 Xi;
    float3 Normal;
  };
  /** @} */

  /**
   * \name Results
   * @{ 
   */
private:
  int m_NumEdgeEdge, m_NumNodeFacet;
  tledCUDADeviceMemoryBlock *mp_EdgeEdgeResultBuffer, *mp_NodeFacetResulBuffer;

protected:
  /** Sets number of results to 0 for all types of contacts */
  virtual void ResetResults(void);

  void SetNumberOfEdgeEdgeContacts(const int num) { m_NumEdgeEdge = num; }
  void SetNumberOfNodeFacetContacts(const int num) { m_NumNodeFacet = num; }

  void SetEdgeEdgeResultBuffer(tledCUDADeviceMemoryBlock &r_buffer) { mp_EdgeEdgeResultBuffer = &r_buffer; }
  void SetNodeFacetResultBuffer(tledCUDADeviceMemoryBlock &r_buffer) { mp_NodeFacetResulBuffer = &r_buffer; }

public:
  int GetNumberOfEdgeEdgeContacts(void) const { return m_NumEdgeEdge; }
  int GetNumberOfNodeFacetContacts(void) const { return m_NumNodeFacet; }

  /** Edge-edge collision events, client is expected to do memory cleanup, i.e. release buffer, when finished. */
  tledCUDADeviceMemoryBlock& GetEdgeEdgeResults(void) { return *mp_EdgeEdgeResultBuffer; }

  /** Edge-edge collision events, client is expected to do memory cleanup, i.e. release buffer, when finished. */
  const tledCUDADeviceMemoryBlock& GetEdgeEdgeResults(void) const { return *mp_EdgeEdgeResultBuffer; }

  /** Node-facet collision events, client is expected to do memory cleanup, i.e. release buffer, when finished. */
  tledCUDADeviceMemoryBlock& GetNodeFacetResults(void) { return *mp_NodeFacetResulBuffer; }

  /** Node-facet collision events, client is expected to do memory cleanup, i.e. release buffer, when finished. */
  const tledCUDADeviceMemoryBlock& GetNodeFacetResults(void) const { return *mp_NodeFacetResulBuffer; }
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledBVHTraverserGPU(void);
  virtual ~tledBVHTraverserGPU(void);
  /** @} */
};

#endif
