// =========================================================================
// File:       tledDeformableRigidContactSolverCPU.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    May 2011
// 
// Copyright (c) 2011, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

#ifndef tledDeformableRigidContactSolverCPU_H
#define tledDeformableRigidContactSolverCPU_H

#include "tledModel.h"
#include "tledContactSolver.h"
#include "tledContactSolverCPU.h"
#include "tledUnstructuredContactManager.h"
#include "tledDeformableRigidBVHTraverserCPU.h"
#include "tledDeformableRigidContactSolver.h"
#include "tledVectorArithmetic.h"
#include "tledHelper.h"

class tledDeformableRigidContactSolverCPU : public tledDeformableRigidContactSolver, public tledContactSolverCPU {
  /**
   * \name Construction and Destruction
   * @{
   */
public:
  static tledDeformableRigidContactSolverCPU* CreateContactSolver(tledUnstructuredContactManager &r_manager, const int surfaceIndex);

  tledDeformableRigidContactSolverCPU(void) {}
  virtual ~tledDeformableRigidContactSolverCPU(void) {}
  /** @} */
};

template <class TDeformableContactMesh, class TRigidContactMesh, class TAPI = tledDeformableRigidContactSolverCPU>
class tledDeformableRigidContactSolverImplCPU : public tledContactSolverImplCPU<TDeformableContactMesh, TAPI> {
public:
  typedef tledContactSolverImplCPU<TDeformableContactMesh, TAPI> Superclass;
  typedef TDeformableContactMesh ContactMesh;
  typedef TRigidContactMesh RigidMesh;
  typedef typename ContactMesh::Facet Facet;

protected:
  typedef typename Superclass::EdgeEdgeConstraintItem EdgeEdgeConstraintItem;
  typedef typename Superclass::NodeFacetConstraintItem NodeFacetConstraintItem;

  /**
   * \name Responses
   * @{ 
   */
protected:
  bool ComputeMasterResponses(float *p_fs, const float uNexts[], const float uCurrs[]);
  bool ComputeSlaveResponses(float *p_fs, const float uNexts[], const float uCurrs[]);
  virtual bool ComputeContactForces(float *p_f, const float uNexts[], const float uCurrs[]);

  virtual bool DoFriction(void) const { return this->GetRigidSurface().GetFrictionCoefficient() > 0; }
  /** @} */

protected:
  virtual tledBVHTraverserCPU* InstantiateBVHTraverser(void);

  /**
   * \name Meshes
   * @{
   */
private:
  const RigidMesh *mpc_RigidSurface;

public:
  const RigidMesh& GetRigidSurface(void) const { return *mpc_RigidSurface; }
  ContactMesh& GetDeformableSurface(void) { return this->GetMesh(); }
  const ContactMesh& GetDeformableSurface(void) const { return this->GetMesh(); }
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledDeformableRigidContactSolverImplCPU(tledUnstructuredContactManager &r_contactManager, const int rigidSurfaceIndex);
  virtual ~tledDeformableRigidContactSolverImplCPU(void) {}
  /** @} */
};

#include "tledDeformableRigidContactSolverCPU.tpp"

#endif
