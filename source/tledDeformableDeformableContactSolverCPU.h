// =========================================================================
// File:       tledDeformableDeformableContactSolverCPU.h
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

#ifndef tledDeformableDeformableContactSolverCPU_H
#define tledDeformableDeformableContactSolverCPU_H

#include "tledContactSolver.h"
#include "tledVectorArithmetic.h"
#include "tledDeformableDeformableContactSolver.h"
#include "tledContactSolverCPU.h"
#include "tledDeformableDeformableBVHTraverserCPU.h"
#include "tledSelfCollisionBVHTraverserCPU.h"

#include <string>
#include <iostream>
#include <algorithm>
#include <cstdlib>

/**
 * \brief CPU solver for deformable-deformable contacts (including self-collisions)
 * \ingroup contact
 */
class tledDeformableDeformableContactSolverCPU : public tledDeformableDeformableContactSolver, public tledContactSolverCPU {
  /**
   * \name Construction, Destruction
   * @{
   */
public:
  static tledDeformableDeformableContactSolverCPU* CreateContactSolver(tledUnstructuredContactManager &r_manager);
  virtual ~tledDeformableDeformableContactSolverCPU(void) {}
  /** @} */  
}; 

/**
 * \brief CPU implementation of deformable-deformable contact solver 
 */
template <class TContactMesh, class TAPI = tledDeformableDeformableContactSolverCPU>
class tledDeformableDeformableContactSolverImplCPU : public tledContactSolverImplCPU<TContactMesh, TAPI> {
public:
  typedef TContactMesh ContactMesh;
  typedef typename ContactMesh::Facet Facet;
  typedef tledContactSolverImplCPU<TContactMesh, TAPI> Superclass;

protected:
  typedef typename Superclass::EdgeEdgeConstraintItem EdgeEdgeConstraintItem;
  typedef typename Superclass::NodeFacetConstraintItem NodeFacetConstraintItem;

protected:
  virtual tledBVHTraverserCPU* InstantiateBVHTraverser(void);
  virtual bool ComputeContactForces(float *p_f, const float uNexts[], const float uCurrs[]);

public:
  virtual bool DoFriction(void) const { return this->GetMesh().GetFrictionCoefficient() > 0; }

public:
  tledDeformableDeformableContactSolverImplCPU(tledUnstructuredContactManager &r_contactRes);
  virtual ~tledDeformableDeformableContactSolverImplCPU(void) {}
};

#include "tledDeformableDeformableContactSolverCPU.tpp"

#endif
