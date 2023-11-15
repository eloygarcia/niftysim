// =========================================================================
// File:       tledSelfCollisionBVHCreator.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    August 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledSelfCollisionBVHCreator_H
#define tledSelfCollisionBVHCreator_H

#include "tledBVHCreator.h"
#include "tledSelfCollisionBVH.h"

#include <set>

/**
 * \brief Generates BVHs where every BV holds only connected primitives
 * \ingroup contact
 */
template <class TBVH>
class tledSelfCollisionBVHCreator : public tledBVHCreatorImpl<TBVH, typename TBVH::ContactMesh> {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledBVHCreatorImpl<TBVH, typename TBVH::ContactMesh> Superclass;
  typedef TBVH BVH;
  typedef typename BVH::BoundingVolume BoundingVolume;
  typedef typename BVH::ContactMesh::Facet Facet;

  static const int BVHOrder = BoundingVolume::NumberOfChildBVs;
  /** @} */

  /**
   * \name Builder Data
   * @{
   */
protected:
  class BVPrimitiveSet;

protected:
  std::vector<int> m_ClusterStartInds;
  std::vector<BoundingVolume> m_PrimitiveBVs;

protected:
  virtual void InitBVHCreationData(void);
  virtual void FreeBVHCreationData(void);

  const std::vector<BoundingVolume>& GetPrimitiveBVs(void) const { return m_PrimitiveBVs; }
  const BoundingVolume& GetPrimitiveBV(const int pInd) const { return GetPrimitiveBVs()[pInd]; }
  BoundingVolume& GetPrimitiveBV(const int pInd) { return m_PrimitiveBVs[pInd]; }
  /** @} */

  /**
   * \name BVH Building Routines
   * @{
   */
private:
  void _ComputeClusterVolinoAxis(float *p_vAxis, const std::vector<int>::const_iterator pIndsBegin, const std::vector<int>::const_iterator pIndsEnd);

  int _ReinsertBVsRecursive(const std::vector<BoundingVolume> &src, const int bvInd, const int parentInd);

  /* Collapses a lower-order hierarchy, eliminating every other level in the hierarchy */
  void _CollapseSubTree(const int subtreeRoot, const int subTreeOrder);  

  /* Starts a BV split by determining the two primitives that are the farthest apart wrt. an axis (can be spatial or object-oriented) */
  void _InitialiseChildPrimitiveSets(const std::vector<int>::iterator &pIndsBegin, const std::vector<int>::iterator &pIndsEnd, const BoundingVolume &bv);
  std::vector<int>::iterator _FindOptimumAdjacentPrimitive(std::vector<int>::iterator pIndsBegin, std::vector<int>::iterator pIndsEnd, const BVPrimitiveSet &cluster, const BoundingVolume &clusterBV);
  /* Recursive BV split */
  void _SplitBV2(std::vector<int>::iterator ppi_childIndexBounds[][2], const BoundingVolume &bv, const std::vector<int>::iterator pIndsBegin, const std::vector<int>::iterator pIndsEnd, const int child1Offset);
  
protected:  
  /** Inserts a BVH node with 2 children (regradless of BV type), only to be used before top-down stage. */
  int InsertParentBinary(const int child0Ind, const int child1Ind);

  const std::vector<int>& GetClusterStartIndices(void) const { return m_ClusterStartInds; }

  std::vector<int> FindClusters(std::vector<int>::iterator i_pIndsBegin, std::vector<int>::iterator i_pIndsEnd) const;
  virtual void SplitBV(std::vector<int>::iterator p_childIndexBounds[][2], const int BVIndex, const std::vector<int>::iterator &pIndicesBegin, const std::vector<int>::iterator &pIndicesEnd);

  /** Performs the bottom-up stage, creating the BVs for disconnected parts of the geometry. */
  virtual void InitialiseClusters(void);

  /** 
   * Hook that is called before bottom-up grouping on BVs listed in r_activeBVs is performed, this can be used to perform grouping of BVs bounding special geometries.
   * The order of the BVs existing in the BVH at this time may (and likely will) subsequently change.<br />
   * At this time in the BVH construction process it can be assumed that there is a one-to-one correspondence between BVs and clusters (i.e. i-th cluster start index <=> i-th BV).
   */
  virtual void PreBottomUpHook(std::vector<int> &r_activeBVs) {}

  virtual void GenerateMain(void);
  /** @} */
};

/** Helper class used in partitioning primitive sets */
template <class TBVH>
class tledSelfCollisionBVHCreator<TBVH>::BVPrimitiveSet {
public:
  typedef TBVH BVH;
  typedef typename BVH::ContactMesh::Facet Facet;

  /**
   * \name Set Operations
   * @{
   */
private:
  std::set<int> m_NodeSet;

public:
  void Insert(const Facet &facet);
  bool IsAdjacent(const Facet &facet) const;
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  BVPrimitiveSet(const Facet &facet) { Insert(facet); }
  BVPrimitiveSet(void) {}
  /** @} */
};

#include "tledSelfCollisionBVHCreator.tpp"
#endif
