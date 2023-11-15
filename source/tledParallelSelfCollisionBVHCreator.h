// =========================================================================
// File:       tledParallelSelfCollisionBVHCreator.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    December 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledParallelSelfCollisionBVHCreator_H
#define tledParallelSelfCollisionBVHCreator_H

#ifdef BOOST_DISABLE_THREADS
#undef _USE_BOOST_
#endif

#ifdef _USE_BOOST_
#include "tledSelfCollisionBVHCreator.h"

#include <boost/thread.hpp>
#include <cassert>

/**
 * \brief CPU-parallel BVH builder for general deformable-deformable and self-collsion detection
 * \ingroup contact
 */
template <class TBVH>
class tledParallelSelfCollisionBVHCreator : public tledSelfCollisionBVHCreator<TBVH> {
  /**
   * \name Types
   * @{
   */
public:  
  typedef tledSelfCollisionBVHCreator<TBVH> Superclass;
  typedef TBVH BVH;
  typedef typename TBVH::ContactMesh ContactMesh;
  typedef typename ContactMesh::Facet Facet;
  typedef typename TBVH::BoundingVolume BoundingVolume;
  /** @} */

  /**
   * \name Mesh
   * @{
   */
public:
  const ContactMesh& GetMesh(void) const { return Superclass::GetMesh(); }
  ContactMesh& GetMesh(void) { return const_cast<ContactMesh&>(Superclass::GetMesh()); }
  /** @} */

  /**
   * \name BVH Generation
   * @{
   */
private:
  void _InsertSubBVH(const BVH &subBVH, const int rootBVIndex, const int parentIndex);
  void _CopyAttributes(const BVH &subBVH);
  void _InitThread(BVH &r_bvh, tledParallelSelfCollisionBVHCreator &r_builder);

protected:
  static void BuilderWorker(tledParallelSelfCollisionBVHCreator &r_builder, std::vector<int> &r_primitiveIndices, const int threadSetTargetSize);
  void SplitSetInitialiseTopDownRecursive(std::vector<int> &r_primitiveIndices, const int threadSetTargetSize);

public:
  virtual void GenerateMain(void);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledParallelSelfCollisionBVHCreator(void) {}
  virtual ~tledParallelSelfCollisionBVHCreator(void) {}
  /** @} */
};

#include "tledParallelSelfCollisionBVHCreator.tpp"
#endif 
#endif
