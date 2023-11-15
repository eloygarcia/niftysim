// =========================================================================
// File:       tledBVH.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    January 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledBVH_H
#define tledBVH_H

#include "tledBVHCreator.h"
#include "tledContactSurface.h"
#include "tledBVHXMLExporter.h"

#include <vector>
#include <limits>

/**
 * \brief Basic bounding-volume hierarchy API
 * \ingroup contact
 */
class tledBVH {
  /**
   * \name Mesh
   * @{
   */
public:
  /** Access to underlying mesh through most general type. */
  virtual const tledContactSurface& GetUnspecifiedMesh(void) const = 0;
  /** @} */

  /**
   * \name Hierarchy Queries
   * @{
   */
public:
  /** Number of leaf BVs = number of mesh primitives. Mainly useful for consistency checks. */
  virtual int GetNumberOfLeafs(void) const = 0;

  /** Maximum depth of a subtree given by root index */
  virtual int GetSubtreeMaxDepth(const int rootIndex) const = 0;

  /** Returns the leafs of a given subtree, reserve vector of suitable size for performance. */
  virtual void GetSubtreeLeafs(std::vector<int> &r_leafIndexBuffer, const int rootIndex) const = 0;

  /** \brief Total hierarchy depth */
  int GetMaxDepth(void) const { return GetSubtreeMaxDepth(0); }
  /** @} */

  /**
   * \name Parameters
   * @{
   */
private:
  float m_BVMargin;

public:
  /** (Max.) Number of children per BVH node */
  virtual int GetBVHOrder(void) const = 0;

  /** Geometry ID of BV type */
  virtual int GetBVGeometryID(void) const = 0;

  /** Safety margin around bounded geometry */
  float GetMargin(void) const { return m_BVMargin; }
  virtual void SetMargin(const float bvMargin) { m_BVMargin = bvMargin; }
  /** @} */

  /**
   * \name XML Export
   * @{
   */
public:
  virtual XMLNode ExportToXML(void) const = 0;
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  virtual void Init(tledBVHCreator &r_builder) = 0;

  /** Called at the end of loading when reading a previous XML export. */
  virtual void LoadFromXMLPostloadHook(void) {}

  tledBVH(void) : m_BVMargin(std::numeric_limits<float>::quiet_NaN()) {}
  virtual ~tledBVH(void) {}
  /** @} */
};

/**
 * Basis of bounding-volume hierarchy implementations
 * \ingroup contact
 */
template <class TContactMesh, class TBV, class TAPI = tledBVH>
class tledBVHImpl : public TAPI {
  /**
   * \name Types
   * @{
   */
public:
  typedef TContactMesh ContactMesh;
  typedef typename ContactMesh::Facet Facet;
  typedef TBV BoundingVolume;
  typedef tledBVHCreatorImpl<tledBVHImpl, ContactMesh> BVHBuilder;
  /** @} */

  /**
   * \name BVH Traits
   * @{
   */
public:
  /** Number of children per node */
  static const int BVHOrder = BoundingVolume::NumberOfChildBVs;
  /** @} */

  /**
   * \name Mesh
   * @{
   */
private:
  const ContactMesh *mpc_Mesh;
  
public:
  const ContactMesh& GetMesh(void) const { return *mpc_Mesh; }
  virtual const tledContactSurface& GetUnspecifiedMesh(void) const { return GetMesh(); }
  virtual int GetNumberOfLeafs(void) const { return this->GetLeafBVIndices().size(); }
  /** @} */

  /**
   * \name Index Queries
   * @{
   */
private:
  std::vector<int> m_LeafBVInds, m_PrimAABBInds;

public:
  std::vector<int>& GetLeafBVIndices(void) { return m_LeafBVInds; }
  const std::vector<int>& GetLeafBVIndices(void) const { return m_LeafBVInds; }

  /** Returns for a given primitive the corresponding BV */
  int GetPrimitiveBVIndex(const int pInd) const { return m_PrimAABBInds[pInd]; }

  /** Returns for a given primitive the corresponding BV */
  std::vector<int>& GetAllPrimitiveBVIndices(void) { return m_PrimAABBInds; }

  /** Returns for a given primitive the corresponding BV */
  const std::vector<int>& GetAllPrimitiveBVIndices(void) const { return m_PrimAABBInds; }

  /**
   * Recursively compiles a list of primitives contained in the subtree whose root is given by subtreeRootAABBInd
   */
  void CompilePrimitiveListRecursive(std::vector<int> &r_dstList, const int subtreeRootAABBInd) const;

  bool IsLeaf(const int bvInd) const { return GetBV(bvInd).PrimitiveIndex >= 0; }
  /** @} */

  /**
   * \name BV Access
   * @{
   */
private:
  std::vector<TBV> m_BVs;

public:
  virtual int GetBVHOrder(void) const { return BVHOrder; }
  virtual int GetBVGeometryID(void) const { return BoundingVolume::BVGeometryID; }

  int GetNumberOfBVs(void) const { return m_BVs.size(); }
  std::vector<BoundingVolume>& GetBVs(void) { return m_BVs; }
  const std::vector<BoundingVolume>& GetBVs(void) const { return m_BVs; }

  BoundingVolume& GetBV(const int bvInd) { return GetBVs()[bvInd]; }
  const BoundingVolume& GetBV(const int bvInd) const { return GetBVs()[bvInd]; }

  /** Maximum depth of a subtree given by root index */
  virtual int GetSubtreeMaxDepth(const int rootIndex) const;

  /** Returns the leafs of a given subtree, reserve vector of suitable size for performance. */
  virtual void GetSubtreeLeafs(std::vector<int> &r_leafIndexBuffer, const int rootIndex) const;
  /** @} */

  /**
   * \name BV Refitting
   *
   *
   * ComputeBounds*: Only computation of bounds, no additional (e.g. self-collision information) is updated<br>
   * Refit*: Actual refitting and comprehensive updating of BV information.
   * @{
   */
protected:
  /** Implements update operations common to all types of BVs */
  virtual void RefitBVCommon(const int bvInd) {}

public:
  virtual void ComputeBoundsFromNodes(BoundingVolume &r_bv, const int *nodeListStart, const int *nodeListEnd) const;
  virtual void ComputeBoundsFromPrimitives(BoundingVolume &r_bv, const int *primitiveListStart, const int *primitiveListEnd) const;
  virtual void ComputeBoundsFromChildren(BoundingVolume &r_bv) const;

  virtual void RefitLeafBV(const int bvInd);
  virtual void RefitInteriorBV(const int bvInd);
  /** @} */

  /**
   * \name Initialisation
   * @{
   */
public:
  /** Initialises the BVH before the first detection pass. */
  virtual void Init(tledBVHCreator &r_bvhBuilder);
  /** @} */

  /**
   * \name Collision Query API
   * @{
   */
public:
  virtual bool DoIntersect(const int bvInd0, const int bvInd1) const { return GetBV(bvInd0).DoesIntersect(GetBV(bvInd1)); }
  virtual bool IsPointInside(const float x[], const int bvInd) const { return GetBV(bvInd).IsInside(x); }
  virtual bool DoesRayIntersect(const float a[], const float b[], const int bvInd) const { return GetBV(bvInd).DoesIntersect(a, b); }
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
  tledBVHImpl(const ContactMesh &mesh) : mpc_Mesh(&mesh) {}
  virtual ~tledBVHImpl(void) {}
  /** @} */
}; /* tledBVH */

#include "tledBVH.tpp"
#endif
