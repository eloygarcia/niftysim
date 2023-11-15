// =========================================================================
// File:       tledGeometricSelfCollisionBVHUpdaterCPU.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    September 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledGeometricSelfCollisionBVHUpdaterCPU_H
#define tledGeometricSelfCollisionBVHUpdaterCPU_H
#include "tledGeometricSelfCollisionBVHUpdater.h"
#include "tledVectorArithmetic.h"
#include "tledMatrixFunctions.h"

template <class TBVH>
class tledGeometricSelfCollisionBVHUpdaterCPU : public tledGeometricSelfCollisionBVHUpdater<TBVH> {
  /**
   * \name Types
   * @{
   */
public:
  typedef TBVH BVH;
  typedef typename TBVH::BoundingVolume BoundingVolume;
  typedef typename TBVH::ContactMesh ContactMesh;
  typedef typename ContactMesh::Facet Facet;
  typedef tledGeometricSelfCollisionBVHUpdater<TBVH> Superclass;
  typedef typename Superclass::UpdateNodeInfo UpdateNodeInfo;
  /** @} */

  /**
   * \name Updating
   * @{
   */
protected:
  /** 
   * \brief Determines whether a hierarchy subtree requires updating by checking the displacements since the last fixed position.
   *
   * Return codes:<br>
   * 0: No significant motion: No update needed what so ever -> INACTIVE
   * 1: Immediate update required, since self-collision possible -> ACTIVE
   * 2: No immediate update required, but have significant non-rigid deformation exceeding BV margin -> NO STATUS FLAG
   * 3: Have significant rigid motion, but bounds are still valid -> RIGID_MOTION
   */
  int ComputeUpdateStatus(UpdateNodeInfo &r_updateInfo);

  /** Updates the status flag of an UpdateNodeInfo object with the value returned by ComputeUpdateStatus. */
  void UpdateSubtreeUpdateStatus(UpdateNodeInfo &r_updateInfo);

  /** Performs the update of a subtree based on the information held in its corresponding UpdateNodeInfo object */
  void PerformSubtreeUpdate(UpdateNodeInfo &r_updateInfo);

  virtual void UpdateBVH(void);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledGeometricSelfCollisionBVHUpdaterCPU(BVH &r_bvh) : Superclass(r_bvh) {}
  tledGeometricSelfCollisionBVHUpdaterCPU(void) {}
  virtual ~tledGeometricSelfCollisionBVHUpdaterCPU(void) {}
  /** @} */
};

#include "tledGeometricSelfCollisionBVHUpdaterCPU.tpp"

#endif
