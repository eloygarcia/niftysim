// =========================================================================
// File:       tledDeformableMovingRigidContactSolverCPU.h
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
#ifndef tledDeformableMovingRigidContactSolverCPU_H
#define tledDeformableMovingRigidContactSolverCPU_H

#include "tledMovingRigidContactSurfaceCPU.h"
#include "tledDeformableRigidContactSolverCPU.h"
#include "tledDeformableMovingRigidBVHTraverserCPU.h"

/**
 * \name Interface for the modelling of contact between deformable and moving rigid bodies (CPU version).
 * \ingroup contact
 */
class tledDeformableMovingRigidContactSolverCPU : public tledDeformableRigidContactSolverCPU {  
  /**
   * \name Construction and Destruction
   * @{
   */
public:
  static tledDeformableRigidContactSolverCPU* CreateContactSolver(tledUnstructuredContactManager &r_manager, const int surfaceIndex);

  tledDeformableMovingRigidContactSolverCPU(void) {}
  virtual ~tledDeformableMovingRigidContactSolverCPU(void) {}
  /** @} */
};

template <class TDeformableContactMesh, class TRigidContactMesh, class TAPI = tledDeformableMovingRigidContactSolverCPU>
class tledDeformableMovingRigidContactSolverImplCPU : public tledDeformableRigidContactSolverImplCPU<TDeformableContactMesh, TRigidContactMesh, TAPI> {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledDeformableRigidContactSolverImplCPU<TDeformableContactMesh, TRigidContactMesh, TAPI> Superclass;
  typedef TDeformableContactMesh ContactMesh;
  typedef TRigidContactMesh RigidMesh;
  typedef typename ContactMesh::Facet Facet;
  typedef typename Superclass::NodeFacetConstraintItem NodeFacetConstraintItem;
  typedef typename Superclass::EdgeEdgeConstraintItem EdgeEdgeConstraintItem;
  /** @} */

  /**
   * \name Responses
   * @{ 
   */
protected:
  float* ComputeRigidSurfaceNodeVelocity(float *p_v, const int nodeIndex) const;

  virtual void ComputeRelativeNodeFacetVelocityMaster(float *p_v, const NodeFacetConstraintItem &ci, const float uNexts[], const float uCurrs[]) const;
  virtual void ComputeRelativeNodeFacetVelocitySlave(float *p_v, const NodeFacetConstraintItem &ci, const float uNexts[], const float uCurrs[]) const;

  virtual void ComputeRelativeEdgeEdgeVelocityMaster(float *p_v, const EdgeEdgeConstraintItem &ci, const float uNexts[], const float uCurrs[]) const;
  virtual void ComputeRelativeEdgeEdgeVelocitySlave(float *p_v, const EdgeEdgeConstraintItem &ci, const float uNexts[], const float uCurrs[]) const;
  /** @} */

protected:
  virtual tledBVHTraverserCPU* InstantiateBVHTraverser(void);

  /**
   * \name Construction and Destruction
   * @{
   */
public:
  tledDeformableMovingRigidContactSolverImplCPU(tledUnstructuredContactManager &r_contactManager, const int rigidSurfaceIndex) : Superclass(r_contactManager, rigidSurfaceIndex) {}
  virtual ~tledDeformableMovingRigidContactSolverImplCPU(void) {}
  /** @} */
};

#include "tledDeformableMovingRigidContactSolverCPU.tpp"
#endif
