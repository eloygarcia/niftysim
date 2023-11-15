// =========================================================================
// File:       tledShellSolverGPU.h
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
#ifndef tledShellSolverGPU_H
#define tledShellSolverGPU_H

#include "tledShellSolver.h"
#include "tledCUDAHelpers.h"

/**
 * \brief tledSolverGPU backend for shell/membrane problems
 * \ingroup shell
 * \ingroup solver
 */
class tledShellSolverGPU : public tledShellSolver {
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
  /** GPU force computation function */
  void ComputeNewForces(void);

  /** Accumulates the effective nodal forces (-shell internal forces) */
  void ComputeNodalForces(float4 *dp_f);
  /** @} */

  /**
   * \name Setup 
   * @{
   */
public:
  class ElementSet : public tledShellSolver::ElementSet {
  public:
    typedef tledShellSolver::ElementSet Superclass;

    /**
     * \name Computation API
     * @{
     */
  public:
    virtual void ComputeForces(void) = 0;
    virtual void ComputeNodalForces(float4 *dp_f) = 0;
    /** @} */

    /**
     * \name Memory Management
     * @{
     */
  public:
    /** Allocates GPU memory, binds texture names, etc. */
    virtual void InitGPU(const tledSolver &mainSolver) = 0;

    /** Deallocates GPU memory */
    virtual void FreeGPU(void) = 0;
    /** @} */

  public:
    ElementSet(tledShellMaterial &r_mat, const std::vector<int> &elInds) : Superclass(r_mat, elInds) {}
    virtual ~ElementSet(void) {}
  };

private:
  template <class TShellElement>
  class _ElementSetImpl : public tledShellSolver::ElementSetImpl<TShellElement, ElementSet> {
  public:
    typedef tledShellSolver::ElementSetImpl<TShellElement, ElementSet> Superclass;

  private:
    typename TShellElement::GPUElement *mdp_GPUElements;
    typename tledShellMaterial::GPUMaterial *mdp_GPUMaterial;
    float4 *mdp_ShellElementForces;
    int *mdp_NodeIndexLookupTable;
    int *mdp_NodeElementVertexLookupTable;
    int2 *mdp_NodeElementVertexLookupBaseIndex;
    int m_NumberOfNodes;
    int m_NodeElementVertexLookupTableSize;

  private:
    void _InitGPU(const tledSolver &mainSolver);
    void _InitDisplacementTextures(void);

  public:
    virtual void InitGPU(const tledSolver &mainSolver);
    virtual void FreeGPU(void);
    virtual void ComputeForces(void);
    virtual void ComputeNodalForces(float4 *dp_u);

  public:
    _ElementSetImpl(tledShellMaterial &r_mat, const typename TShellElement::Surface &surface, const std::vector<int> &elInds) : Superclass(r_mat, surface, elInds) {}
    virtual ~_ElementSetImpl(void) { FreeGPU(); }
  };

protected:
  tledShellSolver::ElementSet* CreateElementSetFromIndices3(tledShellMaterial &r_material, const tledSurface &surface, const std::vector<int> &elInds);

public:
  virtual void Init(tledSolver &r_solver, const tledModel &model);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  virtual ~tledShellSolverGPU(void) {}
  /** @} */
};

#endif
