// =========================================================================
// File:       tledSelfCollisionBVH.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    May 2011
// 
// Copyright (c) 2011, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

#ifndef tledSelfCollisionBVH_H
#define tledSelfCollisionBVH_H

#include "tledHelper.h"
#include "tledDynamicBVH.h"
#include "tledDynamicBVHUpdater.h"
#include "tledSelfCollisionBVHXMLExporter.h"
#include "tledDeformableContactSurface.h"
#include "xmlParser.h"

#include <utility>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <limits>
#include <string>

/**
 * \brief BVH capable of detecting deformable-deformable and self-collisions.  
 * \ingroup contact
 */
class tledSelfCollisionBVH : public tledDynamicBVH {
  /**
   * \name Construction, Destruction
   * @{
   */
public:
  /** Factory function for CPU/GPU (if enabled) deformable geometry BVHs */
  static tledSelfCollisionBVH* CreateBVH(tledDeformableContactSurface &r_mesh, const std::string &bvType, const float margin, const float maxDisplacement, const bool useGPU);

  /** Factory function for loading from XML exports. */
  static tledSelfCollisionBVH* CreateBVH(tledDeformableContactSurface &r_mesh, const XMLNode root, const bool useGPU);

  tledSelfCollisionBVH(void) {}
  virtual ~tledSelfCollisionBVH(void) {}
  /** @} */
};

/**
 * \brief BVH capable of detecting deformable-deformable and self-collisions.  
 * \ingroup contact
 *
 *
 * The template argument TSurface (i.e. the type of the contact surface being bounded) has to be API-compatible (not necessarily via class polymorphism) with tledDeformableContactSurfaceImpl.
 */
template <class TSurface, class TBV, class TAPI = tledSelfCollisionBVH>
class tledSelfCollisionBVHImpl : public tledDynamicBVHImpl<TSurface, tledSelfCollisionBV<TBV>, TAPI> {
  /**
   * \name Types
   * @{
   */
public:
  typedef TSurface ContactMesh;
  typedef typename ContactMesh::Facet Facet;
  typedef tledSelfCollisionBV<TBV> BoundingVolume;
  typedef tledDynamicBVHImpl<ContactMesh, BoundingVolume, TAPI> Superclass;
  typedef typename Superclass::BVHBuilder BVHBuilder;

  static const int BVHOrder = BoundingVolume::NumberOfChildBVs;
  /** @} */

  /**
   * \name Collision Detection
   * @{
   */
private:
  /** \internal Only contains roots of subtrees that might be in collision, no redundancies */
  std::vector<int> m_SelfCollisionCandidateNodes, m_NewSelfCollisionCandidateNodes, m_NonAdjacentChildrenNodes;

public:
  /** List of BVs that contain some geometry that is not topologically connected. */
  const std::vector<int>& GetNonAdjacentGeometryNodes(void) const { return m_NonAdjacentChildrenNodes; }

  /** List of BVs flagged by the self-collision criterion as potentially containing self-collisions. Includes GetNonAdjacentGeometryNodes. */
  const std::vector<int>& GetSelfCollisionCandidates(void) const { return m_SelfCollisionCandidateNodes; }

  std::vector<int>& GetNonAdjacentGeometryNodes(void) { return m_NonAdjacentChildrenNodes; }
  std::vector<int>& GetSelfCollisionCandidates(void) { return m_SelfCollisionCandidateNodes; }

  /** Adds a self-collision candidate node to the list of self-collision candidates, mainly intended for use by BVH builders. */
  void AddSelfCollisionCandidate(const int bvInd);
  /** @} */

  /**
   * \name Initialisation
   * @{
   */
public:
  virtual void Init(tledBVHCreator &r_bvhBuilder);
  /** @} */

  /**
   * \name Update
   * @{
   */
private:
  /* Updates the node information */
  void _UpdateNodePositions(void);

#ifndef NDEBUG
  bool _CheckVolinoComprehensive(const int bvInd) const;
  bool _CheckIntersectionRecursive(const int subtree0Ind, const int subtree1Ind) const;
  bool _CheckCompletenessRecursive(const int subtreeInd) const;
#endif  

protected:
  /* Computes an interior node surface cone from two child cones */
  virtual void ComputeSurfaceConePairwise(float &r_angle, float *p_axis, const float angle0, const float axis0[], const float angle1, const float axis1[]);
  virtual void RefitBVCommon(const int bvInd);

public:
  virtual void RefitLeafBV(const int bvInd);
  virtual void RefitInteriorBV(const int bvInd);

  virtual bool DoesNeedUpdate(const int bvInd) const { return this->GetBV(bvInd).UpdateCounter < this->GetUpdateCounter(); }

  virtual void TranslateSubTree(const int rootBVInd, const float t[]);
  virtual void TransformSubTree(const int rootBVInd, const float m[][3], const float cor[], const float t[]);

  virtual void UpdateBottomUpRecursive(const int bvInd);  

  /** Checks if one surface cone is contained in another. */
  static bool IsConeContainedInCone(const float axis[], const float angle, const float testAxis[], const float testAngle);

  virtual void Update(void);
  /** @} */

  /**
   * \name Hierarchy Topology Queries
   * @{
   */
private:
  std::vector<int> m_GeometryClusterSubtrees;

public:
  /** Returns the BV indices associated with geometry clusters (groups of connected primitives). */
  const std::vector<int>& GetGeometryClusterSubtreeRootIndices(void) const { return m_GeometryClusterSubtrees; }
  std::vector<int>& GetGeometryClusterSubtreeRootIndices(void) { return m_GeometryClusterSubtrees; }
  /** @} */

  /**
   * \name XML Export
   * @{
   */
public:
  virtual XMLNode ExportToXML(void) const;
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  /** Creates a new BVH with a given update algorithm */
  tledSelfCollisionBVHImpl(ContactMesh &r_mesh);
  virtual ~tledSelfCollisionBVHImpl(void) {}
  /** @} */
};

#include "tledSelfCollisionBVH.tpp"

#endif
