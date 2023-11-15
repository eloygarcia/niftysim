// =========================================================================
// File:       tledContactSolver.h
// Purpose:    General contact solver interface
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

#ifndef tledContactSolver_H
#define tledContactSolver_H

#include "tledDeformableContactSurface.h"
#include "tledUnstructuredContactManager.h"
#include "tledDynamicBVH.h"

#include <cassert>
#include <cstdlib>
#include <vector>

class tledSolver;

/**
 * \brief General contact solver front-end API
 */
class tledContactSolver {
  /**
   * \name Responses
   * @{
   */
private:
  bool m_DoMaster;

public:
  /** Master/slave flag, interpretation depends on type of collisions to be resolved */
  virtual bool DoMaster(void) const { return m_DoMaster; }

  /** Toggle master/slave switch */
  virtual void ToggleDoMaster(void) { m_DoMaster = !m_DoMaster; }

  /** Tells if friction responses are to be computed. */
  virtual bool DoFriction(void) const = 0;
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  /** Has to be called before use, but after all contact model parameters have been set. */
  virtual void Init(void) = 0;

  tledContactSolver(void) : m_DoMaster(true) {}
  virtual ~tledContactSolver(void) {}
  /** @} */
}; 

/**
 * \brief Base class for contact solver implementations, provides some members generally useful in solving contact problems
 */
template <class TContactMesh, class TContactSolverInterface>
class tledContactSolverImpl : public TContactSolverInterface {
  /**
   * \name Types
   * @{
   */
public:
  typedef TContactSolverInterface Superclass;
  typedef TContactMesh ContactMesh;
  typedef typename ContactMesh::Facet Facet;
  /** @} */

  /**
   * \name Misc. Getters
   * @{
   */
private:
  ContactMesh &mr_Mesh;
  tledUnstructuredContactManager &mr_Manager;

public:
  /** Contact-surface mesh getter */
  const ContactMesh& GetMesh(void) const { return mr_Mesh; }
  /** Contact-surface mesh getter */
  ContactMesh& GetMesh(void) { return mr_Mesh; }

  tledUnstructuredContactManager& GetManager(void) { return mr_Manager; }
  const tledUnstructuredContactManager& GetManager(void) const { return mr_Manager; }
  /** @} */

  /**
   * \name Response Computation
   * @{
   */
protected:
  /** Distance at which two nodes are considered close, and hence need to be slowed down */
  float GetNodeCloseDistance(void) const { return this->GetManager().GetCloseDistance(); }
  
  /** Simulation time step */
  float GetDt(void) const { return this->GetManager().GetDt(); }
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  virtual void Init(void);

  tledContactSolverImpl(tledUnstructuredContactManager &r_manager);
  virtual ~tledContactSolverImpl(void) {}
  /** @} */
};

#include "tledContactSolver.tpp"

#endif
