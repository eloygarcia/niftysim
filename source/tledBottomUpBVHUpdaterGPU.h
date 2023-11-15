// =========================================================================
// File:       tledBottomUpBVHUpdaterGPU.h
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
#ifndef tledBottomUpBVHUpdaterGPU_H
#define tledBottomUpBVHUpdaterGPU_H

#include "tledDynamicBVHUpdater.h"
#include "tledCUDAHelpers.h"
#include "tledCUDAMemoryBlock.h"

#include <iostream>
#include <cstdlib>
#include <algorithm>

/**
 * \brief Updates a BVH subtree, defined by a set of leafs and a root, bottom up.
 * \ingroup contact
 */
template <class TBVH>
class tledBottomUpBVHUpdaterGPU : public tledDynamicBVHUpdaterImpl<TBVH> {
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
  /** @} */

  /**
   * \name Subtree Definition
   * @{
   */
private:
  const int *mpc_StartIndices;
  int *mdp_Levels;
  int2 *mdp_LevelBounds;
  int m_NumStartIndices;
  int m_SubtreeParentIndex;
  int m_NumLevels;
  std::vector<int2> m_LevelBounds;

protected:
  const int* GetOnDeviceLevels(void) const { return mdp_Levels; }
  const int2* GetOnDeviceLevelIndexBounds(void) const { return mdp_LevelBounds; }
  const int2& GetLevelBounds(const int levelInd) const { return m_LevelBounds[levelInd]; }

public:
  /** @{ */
  /** Start indices of the update */
  void SetStartIndices(const int *startIndices, const int numIndices);
  const int* GetStartIndices(void) const { return mpc_StartIndices; }
  int GetNumberOfStartIndices(void) const { return m_NumStartIndices; }
  /** @} */

  /** @{ */
  /** \brief Node at which updating is stopped, defaults to -1 (i.e. parent of the BVH root) */
  void SetSubtreeRootParent(const int rootInd) { m_SubtreeParentIndex = rootInd; }
  int GetSubtreeRootParent(void) const { return m_SubtreeParentIndex; }
  /** @} */
  /** @} */ 

  /**
   * \name Updating
   * @{
   */
private:
  void _RefitInteriorLevel(const int levelIndex);
  void _RefitLeafLevel(void);

public:
  virtual void UpdateBVH(void);
  /** @} */

  /**
   * \name Initialisation
   * @{
   */
private:
  void _ProgressLevel(std::vector<bool> &r_isUpdated, std::vector<int> &r_nextLevel, const std::vector<int> &levels, const int2 currBounds);
  std::vector<int> _FilterNonUpdateable(const std::vector<int> &potNextLevel, const std::vector<bool> &isUpdated);

public:
  virtual void Init(void);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledBottomUpBVHUpdaterGPU(BVH &r_bvh) : Superclass(r_bvh), mpc_StartIndices(NULL), m_SubtreeParentIndex(-1), mdp_Levels(NULL), mdp_LevelBounds(NULL) {}
  tledBottomUpBVHUpdaterGPU(void) : mpc_StartIndices(NULL), m_SubtreeParentIndex(-1), mdp_Levels(NULL), mdp_LevelBounds(NULL) {}
  virtual ~tledBottomUpBVHUpdaterGPU(void);
  /** @} */
};

#if defined __CUDACC__ && !defined  __GPU_TEST_LINK_CONFLICT_NO_INCLUDE
#include "tledBottomUpBVHUpdaterGPU.tpp"
#endif

#endif
