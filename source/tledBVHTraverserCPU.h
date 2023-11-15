// =========================================================================
// File:       tledBVHTraverserCPU.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    November 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledBVHTraverserCPU_H
#define tledBVHTraverserCPU_H

#include "tledBVHTraverser.h"
#include "tledVectorArithmetic.h"
#include "tledHelper.h"
#include "tledContactSurface.h"

#include <limits>
#include <vector>
#include <algorithm>
#include <utility>

class tledUnstructuredContactManager;

/**
 * \brief Contact search API for CPU contact modelling
 * \ingroup contact
 *
 * Contains only CPU specific functions for output of contact events, for general contact search API see tledBVHTraverser
 */
class tledBVHTraverserCPU : public tledBVHTraverser {
  /**
   * \name Results
   * @{
   */
public:
  /**
   * \brief Basic data-structure holding information about master surface components colliding with the slave.
   */
  struct MasterContactData {
    /** Coordinates of collision, interpretation depends on type of collision */
    float CollisionCoords[3];

    /** Master-surface normal used in projection (always normalised). */
    float Normal[3];

  public:
    /** Resets the projection to a safe default value indicating the absence of a collision */
    virtual void Reset(void) { std::fill(CollisionCoords, CollisionCoords + 3, 1e10f); }
  };

  /** 
   * \brief Holds the information required for computing a node-facet collision response for a known slave node.
   */
  struct MasterFacetContactData : public MasterContactData {
    int FacetIndex;    

  public:
    virtual void Reset(void) { FacetIndex = -1, MasterContactData::Reset(); }

  public:
    /** Initialises the object with the &quot;no-collision&quot; default value */
    MasterFacetContactData(void) { Reset(); }
  };

  /**
   * \brief Edge-edge collision master data
   */
  struct MasterEdgeContactData : public MasterContactData {
    int EdgeIndex;

  public:
    virtual void Reset(void) { EdgeIndex = -1, MasterContactData::Reset(); }

  public:
    /** Initialises the object with the &quot;no-collision&quot; default value */
    MasterEdgeContactData(void) { Reset(); }
  };
  /** @} */

  /**
   * \name Result Buffers
   * @{
   */
protected:
  /** Ordering operator for narrow-phase items */
  struct NarrowPhaseOrdering;

private:
  static std::vector<int> s_ContactNodeIndices, s_ContactEdgeIndices;
  static std::vector<std::pair<int, int> > svv_NodeFacetNarrowPhasePairs, svv_EdgeEdgeNarrowPhasePairs;
  static std::vector<MasterEdgeContactData> s_MasterEdges;
  static std::vector<MasterFacetContactData> s_MasterFacets;  

protected:
  virtual std::vector<std::pair<int, int> >& GetNodeFacetNarrowPhasePairs(void) { return svv_NodeFacetNarrowPhasePairs; }
  virtual std::vector<std::pair<int, int> >& GetEdgeEdgeNarrowPhasePairs(void) { return svv_EdgeEdgeNarrowPhasePairs; }
  virtual const std::vector<std::pair<int, int> >& GetNodeFacetNarrowPhasePairs(void) const { return svv_NodeFacetNarrowPhasePairs; }
  virtual const std::vector<std::pair<int, int> >& GetEdgeEdgeNarrowPhasePairs(void) const { return svv_EdgeEdgeNarrowPhasePairs; }

  virtual std::vector<int>& GetContactNodeBuffer(void) { return s_ContactNodeIndices; }
  virtual std::vector<int>& GetContactEdgeBuffer(void) { return s_ContactEdgeIndices; }

  virtual std::vector<MasterEdgeContactData>& GetEdgeProjectionBuffer(void) { return s_MasterEdges; }
  virtual std::vector<MasterFacetContactData>& GetNodeProjectionBuffer(void) { return s_MasterFacets; }

  virtual const std::vector<MasterEdgeContactData>& GetEdgeProjectionBuffer(void) const { return s_MasterEdges; }
  virtual const std::vector<MasterFacetContactData>& GetNodeProjectionBuffer(void) const { return s_MasterFacets; }

public:
  const tledBVHTraverserCPU::MasterEdgeContactData& GetEdgeProjection(const int slaveEdgeIndex) { return this->GetEdgeProjectionBuffer()[slaveEdgeIndex]; }
  const tledBVHTraverserCPU::MasterFacetContactData& GetNodeProjection(const int slaveNodeIndex) { return this->GetNodeProjectionBuffer()[slaveNodeIndex]; }

  virtual const std::vector<int>& GetSlaveNodeIndices(void) const { return s_ContactNodeIndices; }
  virtual const std::vector<int>& GetSlaveEdgeIndices(void) const { return s_ContactEdgeIndices; }

  /** Resets a collision information bucket after the information has been retrieved by the solver. */
  void ResetEdgeProjection(const int slaveEdgeIndex) { s_MasterEdges[slaveEdgeIndex].Reset(); }

  /** Resets a collision information bucket after the information has been retrieved by the solver. */
  void ResetNodeProjection(const int slaveNodeIndex) { s_MasterFacets[slaveNodeIndex].Reset(); }

  /** True if there are any active slaves in the buffers */
  static bool HaveActiveContacts(void) { return (s_ContactNodeIndices.size() > 0 || s_ContactEdgeIndices.size() > 0); }
  /** @} */

  /**
   * \name Search
   * @{
   */
protected:
  /** Sorts narrow-phase items and removes duplicates, needs to be called before entering narrow-phase search. */
  virtual void FinishBroadPhase(void);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  virtual ~tledBVHTraverserCPU(void) {}
  /** @} */
};
  
/**
 * \brief Base class for CPU contact search classes.
 * \ingroup contact
 * 
 * Public API TAPI must be covariant with tledBVHTraverserCPU and tledBVHTraverser
 */
template <class TMasterBVH, class TSlaveBVH, class TAPI>
class tledBVHTraverserImplCPU : public tledBVHTraverserImpl<TMasterBVH, TSlaveBVH, TAPI> {
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
   * \name Detection
   * @{
   */  
protected: 
  void AddNodeFacetPairs(const int *nodesBegin, const int *nodesEnd, const int masterFacetInd);
  void AddEdgeEdgePairs(const int *slaveEdgesBegin, const int *slaveEdgesEnd, const int *masterEdgesBegin, const int *masterEdgesEnd);

  virtual void AddNarrowPhaseTests(const int sFacetInd, const int mFacetInd); 

  /** Swaps master/slave if requested by client. Non-default behaviour, has to be used to overload AddNarrowPhaseTests. */
  void AddNarrowPhaseTestsSwitching(const int sFacetInd, const int mFacetInd);

  void DoBroadPhaseDetectionRecursive(const int masterStartInd, const int slaveStartInd);

  /** Checks that the node to be tested was in front of the facet before */
  bool DoBehindFacetTest(const float oldSlaveNodePosition[], const float oldFacetNormal[], const float oldFacetV0[]) const;

  /** Tests if there is any mesh intersection in the final configuration, hence testing if there's any need to check the co-planarity criterion */
  bool DoEdgeBehindEdgeTest(const float slave0[], const float slave1[], const float master0[], const float master1[], const float masterNormal0[], const float masterNormal1[]) const;

  /** Projection onto master surface with continuous normals */
  bool ProjectOntoFacetC0Iterative(float *p_xi, float *p_n, const float x[], const float facetNodes[][3], const float nodeNormals[][3]) const;

  /** 
   * \brief Computes for a segment AB and a segment CD the respective parameter (r <- AB; q <- CD) of the point of lease distance.
   * \return true iff both parameters are in [0, 1].
   */
  static bool ComputeEdgeEdgeClosestPointParameters(float &r_r, float &r_q, const float A[], const float B[], const float C[], const float D[]);

  /** 
   * \brief Computes a penetration depth wrt. a normal, given by two node normals and a parameter, for a given pair of edges &amp; parameters. Also performs some checks on the type of motion the geometry underwent during the time period and the surface normals to reject false positives.
   *
   * If t_isMoving* is false the T0 and T1 value for the respective node positions of the edge are assumed to be identical (T0 value used, T1 can be NULL).
   * \return true if the gap value is negative, ie. a penetration has occurred and the motion of the edges can explain this penetration.
   */
  template <const bool t_isMovingSlave, const bool t_isMovingMaster>
  bool ComputeEdgeEdgePenetrationDepth(float &r_g, float *p_n, 
				       const float slave0T0[], const float slave1T0[], const float slave0T1[], const float slave1T1[], const float slaveN0[], const float slaveN1[], const float r, 
				       const float master0T0[], const float master1T0[], const float master0T1[], const float master1T1[], const float masterN0[], const float masterN1[], const float q);

  /** 
   * \brief Checks that the slave node's initial projection lie within the bounds of the facet, and compares the new projection to the current closest point on the master surface.
   * 
   * Even if the test fails, the initial guess buffer is still modified.
   */
  bool DoNodeFacetInitialProjection(float *p_xiInitial, const float slaveNodePosition[], const float facetProjectionOperator[], const float currentBestProjection[]) const;

  /** Compares two projections of the same node, returns true if the first one's the better. */
  bool IsBetterNodeProjection(const float xi0[], const float xi1[]) const;

  /** Compares two projections of the same edge, returns true if the first one's the better. */
  bool IsBetterEdgeProjection(const float xi0[], const float xi1[]) const;

  virtual void ProcessNodeFacetNarrowPhaseItems(const std::vector<std::pair<int, int> >::const_iterator itemsBegin, const std::vector<std::pair<int, int> >::const_iterator itemsEnd) = 0;
  virtual void ProcessEdgeEdgeNarrowPhaseItems(const std::vector<std::pair<int, int> >::const_iterator itemsBegin, const std::vector<std::pair<int, int> >::const_iterator itemsEnd) = 0;

  /** Proecesses all narrow-phase candidates found in broad-phase. */
  virtual void RunNarrowPhase(void);

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

  tledBVHTraverserImplCPU(TSlaveBVH &r_slaveBVH, const TMasterBVH &masterBVH) : Superclass(r_slaveBVH, masterBVH) {}
  virtual ~tledBVHTraverserImplCPU(void) {}
  /** @} */
};

#include "tledBVHTraverserCPU.tpp"
#endif
