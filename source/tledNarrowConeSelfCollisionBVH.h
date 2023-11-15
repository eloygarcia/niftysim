// =========================================================================
// File:       tledNarrowConeSelfCollisionBVH.h
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
#ifndef tledNarrowConeSelfCollisionBVH_H
#define tledNarrowConeSelfCollisionBVH_H

#include "tledSelfCollisionBVH.h"

template <class TSurface, class TBV, class TAPI = tledSelfCollisionBVH>
class tledNarrowConeSelfCollisionBVH : public tledSelfCollisionBVHImpl<TSurface, TBV, TAPI> {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledSelfCollisionBVHImpl<TSurface, TBV, TAPI> Superclass;
  typedef typename Superclass::ContactMesh ContactMesh;
  typedef typename ContactMesh::Facet Facet;
  typedef tledSelfCollisionBV<TBV> BoundingVolume;
  typedef typename Superclass::BVHBuilder BVHBuilder;
  /** @} */

  /**
   * \name Update
   * @{
   */
protected:
  virtual void ComputeSurfaceConePairwise(float &r_angle, float *p_axis, const float angle0, const float axis0[], const float angle1, const float axis1[]);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledNarrowConeSelfCollisionBVH(ContactMesh &r_mesh) : Superclass(r_mesh) {}
  virtual ~tledNarrowConeSelfCollisionBVH(void) {}
  /** @} */
};

#include "tledNarrowConeSelfCollisionBVH.tpp"
#endif
