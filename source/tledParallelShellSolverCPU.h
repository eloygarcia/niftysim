// =========================================================================
// File:       tledParallelShellSolverCPU.h
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
#ifndef tledParallelShellSolverCPU_H
#define tledParallelShellSolverCPU_H
#include "tledParallel.h"
#include "tledShellSolverCPU.h"

#ifdef _USE_BOOST_
class tledParallelShellSolverCPU : public tledShellSolverCPU, public tledParallel {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledShellSolverCPU Superclass;
  /** @} */

  /**
   * \name Computation
   * @{
   */
private:
  float *mp_ThreadResultBuffer;
  int m_NumNodes;
  /** @} */

  /**
   * \name Element Sets
   * @{
   */
protected:
  template <class TShellElement>
  class ElementSetImpl : public tledShellSolverCPU::ElementSetImpl<TShellElement>, public tledParallel {
  private:
    int m_NumNodes;
    float *mp_ThreadResultBuffer;

  public:
    virtual void ComputeForces(float *p_f, const float uCurr[]);
    virtual void ComputeElementThicknesses(float *p_dst, const float U[]) const;

  public:
    ElementSetImpl(tledShellMaterial &r_mat, float *p_resultBuffer, const typename TShellElement::Surface &surface, const std::vector<int> &elInds, const int numNodes) : tledShellSolverCPU::ElementSetImpl<TShellElement>(r_mat, surface, elInds), m_NumNodes(numNodes), mp_ThreadResultBuffer(p_resultBuffer) {}
    virtual ~ElementSetImpl(void) {}
  };

protected:
  virtual tledShellSolver::ElementSet* CreateElementSetFromIndices3(tledShellMaterial &r_material, const tledSurface &surface, const std::vector<int> &elInds);
  /** @} */

  /**
   * \name Construction, Destruction, Initialisation
   * @{
   */
public:
  virtual void Init(tledSolver &solver, const tledModel &model);

  tledParallelShellSolverCPU(void) : mp_ThreadResultBuffer(NULL) {}
  virtual ~tledParallelShellSolverCPU(void);
  /** @} */
};
#endif
#endif
