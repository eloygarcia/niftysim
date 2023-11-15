// =========================================================================
// File:       tledTrackingBVHTraverserGPU.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    December 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledTrackingBVHTraverserGPU_H
#define tledTrackingBVHTraverserGPU_H

#include "tledCUDAHelpers.h"
#include "tledBVHTraverserGPU.h"
#include "tledTrackingBVHTraverser.h"
#include "tledCUDAMemoryBlock.h"

template <class TBaseTraverserGPU>
class tledTrackingBVHTraverserImplGPU : public tledTrackingBVHTraverserImpl<TBaseTraverserGPU> {
  /**
   * \name Imported Types
   * @{
   */
public:
  typedef tledTrackingBVHTraverserImpl<TBaseTraverserGPU> Superclass;
  typedef typename Superclass::MasterBVH MasterBVH;
  typedef typename Superclass::SlaveBVH SlaveBVH;
  typedef typename MasterBVH::ContactMesh MasterMesh;
  typedef typename SlaveBVH::ContactMesh SlaveMesh;

  typedef tledBVHTraverserGPU::NodeFacetNarrowPhaseResult<MasterMesh::Facet::NumberOfVertices> NodeFacetNarrowPhaseResult;
  typedef tledBVHTraverserGPU::EdgeEdgeNarrowPhaseResult EdgeEdgeNarrowPhaseResult;
  /** @} */

  /**
   * \name Tracking
   * @{
   */
private:
  int2 *mdp_SlaveNodeEdgeNeighbourRanges, *mdp_SlaveNodeFacetNeighbourRanges;
  int *mdp_SlaveNodeEdgeNeighbours, *mdp_SlaveNodeFacetNeighbours;

  int2 *mdp_MasterNodeEdgeNeighbourRanges, *mdp_MasterNodeFacetNeighbourRanges;
  int *mdp_MasterNodeEdgeNeighbours, *mdp_MasterNodeFacetNeighbours;

  int m_NumberOfTrackedSlaveNodeFacetContacts, m_NumberOfTrackedSlaveEdgeEdgeContacts;
  tledCUDADeviceMemoryBlock *mp_TrackedSlaveNodeFacetContacts, *mp_TrackedSlaveEdgeEdgeContacts;

  int m_NumberOfTrackedMasterNodeFacetContacts, m_NumberOfTrackedMasterEdgeEdgeContacts;
  tledCUDADeviceMemoryBlock *mp_TrackedMasterNodeFacetContacts, *mp_TrackedMasterEdgeEdgeContacts;

private:
  template <const int t_numFacetVertices>
  static void _UpdateTrackedNodeFacetContacts(tledCUDADeviceMemoryBlock* &rp_nodeFacetCandidatePairs, int *dp_counter, cudaStream_t stream, const tledCUDADeviceMemoryBlock &resultBuffer, const int numContacts, const int *dpc_masterNodeFacetNeighbours, const int2 *dpc_masterNodeFacetNeighbourRanges, const int maxRangeSize);
  static void _UpdateTrackedEdgeEdgeContacts(tledCUDADeviceMemoryBlock* &rp_edgeEdgeCandidatePairs, int *dp_counter, cudaStream_t stream, const tledCUDADeviceMemoryBlock &resultBuffer, const int numContacts, const int *dpc_masterNodeEdgeNeigbours, const int2 *dpc_masterNodeEdgeNeigbourRanges, const int *dpc_slaveNodeEdgeNeigbours, const int2 *dpc_slaveNodeEdgeNeigbourRanges, const int maxMasterRangeSize, const int maxSlaveRangeSize);

protected:
  static void ConvertHostToDeviceNeighbourLists(int* &rdp_list, int2* &rdp_ranges, const std::vector<int> &list, const std::vector<std::pair<int, int> > &ranges);

  void SetSlaveNodeEdgeNeighbours(int *dp_neighbours, int2 *dp_ranges, const int maxRangeSize);
  void SetSlaveNodeFacetNeighbours(int *dp_neighbours, int2 *dp_ranges, const int maxRangeSize);
  void DeallocateSlaveNeighbourData(void);

  const int* GetSlaveNodeEdgeNeighbours(void) const { return mdp_SlaveNodeEdgeNeighbours; }
  const int2* GetSlaveNodeEdgeNeighbourRanges(void) const { return mdp_SlaveNodeEdgeNeighbourRanges; }
  const int* GetSlaveNodeFacetNeighbours(void) const { return mdp_SlaveNodeFacetNeighbours; }
  const int2* GetSlaveNodeFacetNeighbourRanges(void) const { return mdp_SlaveNodeFacetNeighbourRanges; }

  void SetMasterNodeEdgeNeighbours(int *dp_neighbours, int2 *dp_ranges, const int maxRangeSize);
  void SetMasterNodeFacetNeighbours(int *dp_neighbours, int2 *dp_ranges, const int maxRangeSize);
  void DeallocateMasterNeighbourData(void);

  const int* GetMasterNodeEdgeNeighbours(void) const { return mdp_MasterNodeEdgeNeighbours; }
  const int2* GetMasterNodeEdgeNeighbourRanges(void) const { return mdp_MasterNodeEdgeNeighbourRanges; }
  const int* GetMasterNodeFacetNeighbours(void) const { return mdp_MasterNodeFacetNeighbours; }
  const int2* GetMasterNodeFacetNeighbourRanges(void) const { return mdp_MasterNodeFacetNeighbourRanges; }

  int GetNumberOfTrackedSlaveNodeFacetContacts(void) const { return m_NumberOfTrackedSlaveNodeFacetContacts; }
  tledCUDADeviceMemoryBlock* GetSlaveNodeFacetTrackingBuffer(void) { return mp_TrackedSlaveNodeFacetContacts; }
  void SetSlaveNodeFacetTrackingBuffer(tledCUDADeviceMemoryBlock *p_buffer) { mp_TrackedSlaveNodeFacetContacts = p_buffer; }

  int GetNumberOfTrackedSlaveEdgeEdgeContacts(void) const { return m_NumberOfTrackedSlaveEdgeEdgeContacts; }
  tledCUDADeviceMemoryBlock* GetSlaveEdgeEdgeTrackingBuffer(void) { return mp_TrackedSlaveEdgeEdgeContacts; }
  void SetSlaveEdgeEdgeTrackingBuffer(tledCUDADeviceMemoryBlock *p_buffer) { mp_TrackedSlaveEdgeEdgeContacts = p_buffer; }

  int GetNumberOfTrackedMasterNodeFacetContacts(void) const { return m_NumberOfTrackedMasterNodeFacetContacts; }
  tledCUDADeviceMemoryBlock* GetMasterNodeFacetTrackingBuffer(void) { return mp_TrackedMasterNodeFacetContacts; }
  void SetMasterNodeFacetTrackingBuffer(tledCUDADeviceMemoryBlock *p_buffer) { mp_TrackedMasterNodeFacetContacts = p_buffer; }

  int GetNumberOfTrackedMasterEdgeEdgeContacts(void) const { return m_NumberOfTrackedMasterEdgeEdgeContacts; }
  tledCUDADeviceMemoryBlock* GetMasterEdgeEdgeTrackingBuffer(void) { return mp_TrackedMasterEdgeEdgeContacts; }
  void SetMasterEdgeEdgeTrackingBuffer(tledCUDADeviceMemoryBlock *p_buffer) { mp_TrackedMasterEdgeEdgeContacts = p_buffer; }

  virtual void ExtractTrackingInformationFromContacts(void);  
  /** @} */  

  /**
   * \name Contact Search
   * @{
   */
protected:
  virtual void PreFullBVHTraversalHook(void);
  virtual void RunTrackingBroadPhase(void);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:  
  tledTrackingBVHTraverserImplGPU(SlaveBVH &r_slaveBVH, const MasterBVH &masterBVH);
  virtual ~tledTrackingBVHTraverserImplGPU(void) {}
  /** @} */
};

#endif
