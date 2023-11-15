// =========================================================================
// File:       tledParallelGeometricSelfCollisionBVHUpdaterCPU.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    February 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledParallelGeometricSelfCollisionBVHUpdaterCPU_H
#define tledParallelGeometricSelfCollisionBVHUpdaterCPU_H
#ifdef BOOST_DISABLE_THREADS
#undef _USE_BOOST_
#endif

#ifdef _USE_BOOST_
#include "tledGeometricSelfCollisionBVHUpdaterCPU.h"

#include <boost/thread.hpp>

/**
 * \brief CPU-parallel updater for self-collision BVHs
 * \ingroup contact
 */
template <class TBVH>
class tledParallelGeometricSelfCollisionBVHUpdaterCPU : public tledGeometricSelfCollisionBVHUpdaterCPU<TBVH> {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledGeometricSelfCollisionBVHUpdaterCPU<TBVH> Superclass;
  typedef TBVH BVH;
  typedef typename Superclass::BoundingVolume BoundingVolume;
  typedef typename Superclass::ContactMesh ContactMesh;
  typedef typename Superclass::Facet Facet;
  typedef typename Superclass::UpdateNodeInfo UpdateNodeInfo;
  /** @} */

  /**
   * \name Parallelism
   * @{
   */
private:
  int m_NumThreads;

public:
  int GetNumberOfThreads(void) { return m_NumThreads; }
  void SetNumberOfThreads(int numThreads) { m_NumThreads = numThreads; }
  /** @} */

  /**
   * \name Updating
   * @{
   */
protected:
  static void UpdateUpdateStatusWorker(tledParallelGeometricSelfCollisionBVHUpdaterCPU *p_updater, UpdateNodeInfo *p_updateNode);

  virtual void UpdateBVH(void);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledParallelGeometricSelfCollisionBVHUpdaterCPU(BVH &r_bvh) : Superclass(r_bvh), m_NumThreads(boost::thread::hardware_concurrency()) {}
  tledParallelGeometricSelfCollisionBVHUpdaterCPU(void) : m_NumThreads(boost::thread::hardware_concurrency()) {}
  virtual ~tledParallelGeometricSelfCollisionBVHUpdaterCPU(void) {}
  /** @} */
};

#include "tledParallelGeometricSelfCollisionBVHUpdaterCPU.tpp"
#endif
#endif
