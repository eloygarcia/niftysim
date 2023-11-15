// =========================================================================
// File:       tledGreedySelfCollisionBVHUpdater.h
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
#ifndef tledGreedySelfCollisionBVHUpdater_H
#define tledGreedySelfCollisionBVHUpdater_H

/**
 * \brief Greedy update strategy for self-collision BVHs
 * \ingroup contact
 */
template <class TBVH>
class tledGreedySelfCollisionBVHUpdater : public tledDynamicBVHUpdaterImpl<TBVH> {
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
   * \name Updating
   * @{
   */
public:
  virtual void UpdateBVH(void);
  /** @} */

  /**
   * \name Initialisation
   * @{
   */
public:
  virtual void Init(void) {}
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledGreedySelfCollisionBVHUpdater(BVH &r_bvh) : Superclass(r_bvh) {}
  tledGreedySelfCollisionBVHUpdater(void) {}
  virtual ~tledGreedySelfCollisionBVHUpdater(void) {}
  /** @} */
};

#include "tledGreedySelfCollisionBVHUpdater.tpp"
#endif
