// =========================================================================
// File:       tledGeometricSelfCollisionBVHUpdater.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    April 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledGeometricSelfCollisionBVHUpdater_H
#define tledGeometricSelfCollisionBVHUpdater_H

#include "tledDynamicBVHUpdater.h"
#include "tledVectorArithmetic.h"

#include <algorithm>
#include <vector>
#include <cassert>

/**
 * \brief Update strategy for self-collision BVHs based on geometric criteria. 
 * \ingroup contact
 */
template <class TBVH>
class tledGeometricSelfCollisionBVHUpdater : public tledDynamicBVHUpdaterImpl<TBVH> {
  /**
   * \name Types
   * @{
   */
public:
  typedef TBVH BVH;
  typedef typename TBVH::BoundingVolume BoundingVolume;
  typedef typename TBVH::ContactMesh ContactMesh;
  typedef typename ContactMesh::Facet Facet;
  typedef tledDynamicBVHUpdaterImpl<TBVH> Superclass;

  /**
   * \brief Holds the information needed for the default update.
   */
  struct UpdateNodeInfo;
  /** @} */

  /**
   * \name Geometry 
   * @{
   */
protected:
  std::vector<float> m_UpdateNodeLastUpdateNodePositions, m_UpdateNodeCurrPositions;

public:
  /** 
   * \brief Computes the controid of nodes bounded by a BVH update node.
   */
  void ComputeBoundedNodeCentroid(UpdateNodeInfo &r_updateInfo);

  std::vector<int>& GetUpdateNodeIndices(void) { return m_UpdateLevelNodeNodeIndices; }
  const std::vector<int>& GetUpdateNodeIndices(void) const { return m_UpdateLevelNodeNodeIndices; }

  const std::vector<float>& GetUpdateNodePreviousPositions(void) const { return m_UpdateNodeLastUpdateNodePositions; }
  std::vector<float>& GetUpdateNodePreviousPositions(void) { return m_UpdateNodeLastUpdateNodePositions; }

  const std::vector<float>& GetUpdateNodeCurrentPositions(void) const { return m_UpdateNodeCurrPositions; }
  std::vector<float>& GetUpdateNodeCurrentPositions(void) { return m_UpdateNodeCurrPositions; }
  /** @} */

  /**
   * \name Updating
   * @{
   */
protected:
  std::vector<UpdateNodeInfo> m_DefaultUpdateInfo;
  std::vector<int> m_UpdateLevelNodePrimitiveIndices, m_UpdateLevelNodeNodeIndices;

public:
  /** Possible self-coll.: Nodes with this status need to be checked for any type of collision and updated accordingly */
  static const int ACTIVE = 0x1;

  /** (Near-) Zero velocity, no updated needed under any circumstances */
  static const int INACTIVE = (1 << 1); 

  /**
   * Non-adjacent: no Volino, checks assume potentially self-colliding by default, unless inactive.
   */
  static const int NONADJACENT = (1 << 2);

  /**
   * If this flag is set, no self-collision is possible regardless of surface configuration. Sub-tree can be updated with a rigid transform.
   */
  static const int RIGID_MOTION = (1 << 3);    

  /**
   * Flag indicating a translated tree (no rotation unlike with RIGID_MOTION)
   */
  static const int TRANSLATION = (1 << 4);

public:
#ifndef NDEBUG
  bool IsUpdateNode(const int bvInd) const;
#endif

  /** Saves the state of an update node for gauging the deformation in the next time step. */
  void StoreUpdateNodePositions(UpdateNodeInfo &r_updateInfo);

  /** Returns the vector of Update start nodes */
  std::vector<UpdateNodeInfo>& GetUpdateNodes(void) { return m_DefaultUpdateInfo; }
  const std::vector<UpdateNodeInfo>& GetUpdateNodes(void) const { return m_DefaultUpdateInfo; }
  /** @} */

  /**
   * \name Initialisation
   * @{
   */
protected:
  float m_MaxConeAngle;
  int m_MaxNumUN, m_MinNumUN;

private:
#ifndef NDEBUG
  bool _CheckPathHasUpdateNode(const int bvInd) const;
#endif

protected:
  void InitialiseUpdateNodePositionBuffers(void);
#ifndef NDEBUG
  void CheckContainmentRecursive(std::vector<int> &r_nodeInds, const int bvInd, const float trans[]) const;
#endif

public:
  /** Setter for max. number of update start nodes */
  void SetMaxNumberOfUpdateNodes(const int numUN) { m_MaxNumUN = numUN; }

  /** Setter for min. number of update start nodes */
  void SetMinNumberOfUpdateNodes(const int numUN) { m_MinNumUN = numUN; }

  /** Setter for surface cone threshold angle for update-start node initialisation */
  void SetMaxUpdateNodeTargetConeAngle(const float coneAngle) { m_MaxConeAngle = coneAngle; }

  virtual void Init(void);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
private:
  void _SetDefaults();

public:
  tledGeometricSelfCollisionBVHUpdater(BVH &r_bvh) : Superclass(r_bvh) { _SetDefaults(); }
  tledGeometricSelfCollisionBVHUpdater(void) { _SetDefaults(); }
  virtual ~tledGeometricSelfCollisionBVHUpdater(void) {}
  /** @} */
};

/**
 * \brief Holds the information about/for deciding the "up-to-dateness" of a BVH subtree
 */
template <class TBVH>
struct tledGeometricSelfCollisionBVHUpdater<TBVH>::UpdateNodeInfo {
  int BVIndex;
  int Status;
  int NodeStartIndex, NodeEndIndex;
  float LastUpdateCentroid[3], CurrentCentroid[3];
  float Rotation[3][3], Translation[3];
  float LastUpdateConeAngle, MaxNonRigidDisplacement;
  int RigidUpdateCounter;
  
#ifndef NDEBUG
  UpdateNodeInfo(void) {
    BVIndex = NodeStartIndex = NodeEndIndex = -1;
    std::fill(&Rotation[0][0], &Rotation[2][3], std::numeric_limits<float>::quiet_NaN());
    std::fill(Translation, Translation + 3, std::numeric_limits<float>::quiet_NaN());
    std::fill(LastUpdateCentroid, LastUpdateCentroid + 3, std::numeric_limits<float>::quiet_NaN());
    std::fill(CurrentCentroid, CurrentCentroid + 3, std::numeric_limits<float>::quiet_NaN());
    LastUpdateConeAngle = MaxNonRigidDisplacement = std::numeric_limits<float>::quiet_NaN();
    RigidUpdateCounter = -1;
    Status = -1;
  }
#endif
};

#include "tledGeometricSelfCollisionBVHUpdater.tpp"
#endif
