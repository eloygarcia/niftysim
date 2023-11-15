// =========================================================================
// File:       tledTrackingBVHTraverser.h
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
#ifndef tledTrackingBVHTraverser_H
#define tledTrackingBVHTraverser_H

#include "tledBVHTraverser.h"

/**
 * \brief Tracking-specific API for tracking BVH traversers
 * \ingroup contact
 */
class tledTrackingBVHTraverser {
  /**
   * \name Tracking Control
   * @{
   */
private:
  int m_NumTrackingIts;

public:
  /**
   * \brief Number of iterations between full BVH traversals when only tracking is employed.
   * @{
   */
  int GetNumberOfTrackingIterations(void) const { return m_NumTrackingIts; }
  virtual void SetNumberOfTrackingIterations(const int numIterations) { m_NumTrackingIts = numIterations; }
  /** @} */
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledTrackingBVHTraverser(void) : m_NumTrackingIts(-1) {}
  virtual ~tledTrackingBVHTraverser(void) {}
  /** @} */  
};

/**
 * \brief Base class for BVH traverser that can track contacts
 * \ingroup contact
 */
template <class TBaseTraverser>
class tledTrackingBVHTraverserImpl : public TBaseTraverser {
  /**
   * \name Imported Types
   * @{
   */
public:
  typedef TBaseTraverser Superclass;
  typedef typename Superclass::MasterBVH MasterBVH;
  typedef typename Superclass::SlaveBVH SlaveBVH;
  typedef typename MasterBVH::ContactMesh MasterMesh;
  typedef typename SlaveBVH::ContactMesh SlaveMesh;
  /** @} */

  /**
   * \name Tracking Topological Information
   * @{
   */
private:
  int m_MaxSlaveNodeNodeNeighbourRangeSize, m_MaxSlaveNodeEdgeNeighbourRangeSize, m_MaxSlaveNodeFacetNeighbourRangeSize;
  int m_MaxMasterNodeNodeNeighbourRangeSize, m_MaxMasterNodeEdgeNeighbourRangeSize, m_MaxMasterNodeFacetNeighbourRangeSize; 

private:
  static void _ConvertNodeListsToLinearList(std::vector<int> &r_list, std::vector<std::pair<int, int> > &r_ranges, int &r_maxRangeSize, const std::vector<std::vector<int> > &nodeLists);

protected:
  void SetMaxSlaveNodeEdgeNeighbourRangeSize(const int size) { m_MaxSlaveNodeEdgeNeighbourRangeSize = size; }
  void SetMaxSlaveNodeFacetNeighbourRangeSize(const int size) { m_MaxSlaveNodeFacetNeighbourRangeSize = size; }

  void SetMaxMasterNodeEdgeNeighbourRangeSize(const int size) { m_MaxMasterNodeEdgeNeighbourRangeSize = size; }
  void SetMaxMasterNodeFacetNeighbourRangeSize(const int size) { m_MaxMasterNodeFacetNeighbourRangeSize = size; }

  int GetMaxSlaveNodeEdgeNeighbourRangeSize(void) const { return m_MaxSlaveNodeEdgeNeighbourRangeSize; }
  int GetMaxSlaveNodeFacetNeighbourRangeSize(void) const { return m_MaxSlaveNodeFacetNeighbourRangeSize; }

  int GetMaxMasterNodeEdgeNeighbourRangeSize(void) const { return m_MaxMasterNodeEdgeNeighbourRangeSize; }
  int GetMaxMasterNodeFacetNeighbourRangeSize(void) const { return m_MaxMasterNodeFacetNeighbourRangeSize; }

  template <class TSurface>
  static void ExtractNodeFacetNeighbourHood(std::vector<int> &r_list, std::vector<std::pair<int, int> > &r_ranges, int &r_maxRangeSize, const TSurface &surface);
  template <class TSurface>
  static void ExtractNodeEdgeNeighbourHood(std::vector<int> &r_list, std::vector<std::pair<int, int> > &r_ranges, int &r_maxRangeSize, const TSurface &surface);
  /** @} */

  /**
   * \name Collision Detection
   * @{
   */
private:
  int m_TrackingStep;

protected:
  /** \brief Steps since last full search */
  int GetTrackingStep(void) const { return m_TrackingStep; }
  
  /** \brief Simplified, tracking-based broad-phase contact search. */
  virtual void RunTrackingBroadPhase(void) = 0;      

  /** \brief Hook called at end of narrow-phase, extracts the information necessary for tracking in the next steps */
  virtual void ExtractTrackingInformationFromContacts(void) = 0;

  /** Allows for clean up of unused memory allocated for tracking prior to a full BVH traversal */
  virtual void PreFullBVHTraversalHook(void) {}

  virtual void RunBroadPhase(void);
  virtual void RunNarrowPhase(void);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  virtual void Init(tledUnstructuredContactManager &r_manager);

  tledTrackingBVHTraverserImpl(SlaveBVH &r_slaveBVH, const MasterBVH &masterBVH);
  virtual ~tledTrackingBVHTraverserImpl(void) {}
  /** @} */  
};

#include "tledTrackingBVHTraverser.tpp"

#endif
