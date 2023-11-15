// =========================================================================
// File:       tledNarrowConeSelfCollisionBVHGPU.h
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
#ifndef tledNarrowConeSelfCollisionBVHGPU_H
#define tledNarrowConeSelfCollisionBVHGPU_H

#include "tledNarrowConeSelfCollisionBVH.h"
#include "tledSelfCollisionBVHGPU.h"

#ifdef _GPU_
#include <vector_functions.h>
#endif

template <class TSurface, class TBV, class TAPI = tledSelfCollisionBVHGPU>
class tledNarrowConeSelfCollisionBVHImplGPU : public tledBVHImplGPU<tledNarrowConeSelfCollisionBVH<TSurface, TBV, TAPI> > {    
  /**
   * \name Types
   * @{
   */
public:
  typedef TSurface ContactMesh;
  typedef typename ContactMesh::Facet Facet;
  typedef tledSelfCollisionBV<TBV> BoundingVolume;
  typedef tledBVHImplGPU<tledNarrowConeSelfCollisionBVH<TSurface, TBV, tledSelfCollisionBVHGPU> > Superclass;
  typedef typename Superclass::BVHBuilder BVHBuilder;
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledNarrowConeSelfCollisionBVHImplGPU(ContactMesh &r_mesh) : Superclass(r_mesh) {}
  virtual ~tledNarrowConeSelfCollisionBVHImplGPU(void) {}
  /** @} */
};
#endif
