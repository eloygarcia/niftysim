// =========================================================================
// File:       tledShellSolverCPU.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    March 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledShellSolverCPU_H
#define tledShellSolverCPU_H

#include "tledShellSolver.h"
#include "tledSolverCPU.h"

/**
 * \brief tledSolverCPU backend for shell/membrane problems
 * \ingroup shell
 * \ingroup solver
 */
class tledShellSolverCPU : public tledShellSolver {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledShellSolver Superclass;
  /** @} */

  /**
   * \name Solution Computation
   * @{
   */
public:
  /** CPU force computation function */
  virtual void ComputeNewForces(float *p_F, const float u[]);

  /** Computes the element bending moments for all elements in a given element set (3 components/element, format: x, y, xy) */
  virtual void ComputeElementSetBendingMoments(float *p_M, const int elSetIndex, const float uCurr[]) const;
  /** @} */

  /**
   * \name Setup 
   * @{
   */
public:
  class ElementSet : public tledShellSolver::ElementSet {
  public:
    typedef tledShellSolver::ElementSet Superclass;

  public:
    virtual void ComputeForces(float *p_f, const float uCurr[]) = 0;
    /** 
     * \brief Computes the three bending moments (Mx, My, Mxy) for each element in the element set (local numbering).
     *
     * Not defined for all element types, default is an all-zero field 
     * \param p_M, 3*#elements bending moment buffer.
     */
    virtual void ComputeElementBendingMoments(float *p_M, const float uCurr[]) const = 0;

  public:
    ElementSet(tledShellMaterial &r_mat, const std::vector<int> &elInds) : Superclass(r_mat, elInds) {}
    virtual ~ElementSet(void) {}
  };

protected:
  template <class TShellElement>
  class ElementSetImpl : public tledShellSolver::ElementSetImpl<TShellElement, ElementSet> {
  public:
    typedef tledShellSolver::ElementSetImpl<TShellElement, ElementSet> Superclass;

  public:
    virtual void ComputeForces(float *p_f, const float uCurr[]);

    virtual void ComputeElementBendingMoments(float *p_M, const float uCurr[]) const {
      std::fill(p_M, p_M + 3*this->m_Elements.size(), 0.0f);
    }

  protected:
    ElementSetImpl(tledShellMaterial &r_mat) : Superclass(r_mat) {}

  public:
    ElementSetImpl(tledShellMaterial &r_mat, const typename TShellElement::Surface &surface, const std::vector<int> &elInds) : Superclass(r_mat, surface, elInds) {}

    virtual ~ElementSetImpl(void) {}
  };

protected:
  virtual tledShellSolver::ElementSet* CreateElementSetFromIndices3(tledShellMaterial &r_material, const tledSurface &surface, const std::vector<int> &elInds);

public:
  virtual void Init(tledSolver &solver, const tledModel &model) { Superclass::Init(solver, model); }
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledShellSolverCPU(void) {}
  virtual ~tledShellSolverCPU(void) {}
  /** @} */
};
#endif
