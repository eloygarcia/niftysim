// =========================================================================
// File:       tledBVHCreator.h
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
#ifndef tledBVHCreator_H
#define tledBVHCreator_H

#include "tledHelper.h"
#include "tledVectorArithmetic.h"

#include <algorithm>
#include <vector>
#include <cassert>
#include <iterator>
#ifndef NDEBUG
#include <iostream>
#endif

class tledBVH;

/**
 * \brief Basic BVH generator API
 */
class tledBVHCreator {
  /**
   * \name BVH Generation
   * @{
   */
public:
  virtual void Generate(void) = 0;
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  virtual ~tledBVHCreator(void) {}
  /** @} */
};

/**
 * \brief Base class for BVH generator implementations 
 * \ingroup contact
 */
template <class TBVH, class TMesh>
class tledBVHCreatorImpl : public tledBVHCreator {
  /**
   * \name Types
   * @{
   */
public:
  typedef TMesh ContactMesh;
  typedef typename TBVH::BoundingVolume BoundingVolume;
  typedef typename ContactMesh::Facet Facet;
  /** @} */

  /**
   * \name Mesh 
   * @{
   */
private:
  const TMesh *mpc_Mesh;

public:
  void SetMesh(const ContactMesh &mesh) { mpc_Mesh = &mesh; }
  const ContactMesh& GetMesh(void) const { return *mpc_Mesh; }
  /** @} */

  /**
   * \name Generation and Output
   * @{
   */
protected:
  TBVH *mp_BVH;
  std::vector<float> m_PrimitiveCentroids;
  std::vector<int> m_BVPrimitiveIndices;

protected:
  /** Initialise data required for BVH construction */
  virtual void InitBVHCreationData(void);

  /** Free data used in BVH construction */
  virtual void FreeBVHCreationData(void);

  /** Initialise the output data structure */
  virtual void InitBVHStruct(void);

  void InitLeafData(void);

  /** Returns the vector of primitive centroids, only available between _InitBVHCreationData and _FreeBVHCreationData */
  const std::vector<float>& GetPrimitiveCentroids(void) const { return m_PrimitiveCentroids; }
  std::vector<float>& GetPrimitiveCentroids(void) { return m_PrimitiveCentroids; }

  /** Returns a reference to the vector of BV primitives used in initialisation */
  std::vector<int>& GetBVPrimitiveIndices(void) { return m_BVPrimitiveIndices; }

  /** Splits the BV into two or more child BVs */
  virtual void SplitBV(std::vector<int>::iterator p_childIndexBounds[][2], const int BVIndex, const std::vector<int>::iterator &pIndicesBegin, const std::vector<int>::iterator &pIndicesEnd) = 0;

  virtual void InitialiseTopDownRecursive(const int aabbInd, const std::vector<int>::iterator pIndsBegin, const std::vector<int>::iterator pIndsEnd);

  virtual void GenerateMain(void) = 0;

  void ComputeCentroids(float *p_dst) const;

public:
  /** 
   * Inserts a new BV, sets parent index to the user specified value (children, primitive indices initialised with -1 if debug).
   *
   *
   * \return Index of newly inserted BV.
   */
  int AddBV(const int parentIndex);

  /**
   * \brief Create BVH in user specified buffer
   */
  virtual void SetBVHBuffer(TBVH &r_dst);

  /**
   * \brief Dynamically allocate output buffer (memory clean-up is responsibility of client)
   *
   *
   * Mesh must be set prior to allocation.
   */
  void AllocateBVHBuffer(void) { SetBVHBuffer(*(new TBVH(GetMesh()))); }

  TBVH& GetOutput(void) { return *mp_BVH; }
  const TBVH& GetOutput(void) const { return *mp_BVH; }

  virtual void Generate(void);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  virtual ~tledBVHCreatorImpl(void) {}
  /** @} */
}; 

#ifndef NDEBUG
template <class TMesh>
bool CheckConnectivity(const std::vector<int> &primIndList, const TMesh &mesh);
#endif

#include "tledBVHCreator.tpp"
#endif
