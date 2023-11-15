// =========================================================================
// File:       tledTractionConstraint.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    March 2014
// 
// Copyright (c) 2014, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledTractionConstraint_H
#define tledTractionConstraint_H

#include "tledSurfaceConstraint.h"

#include <algorithm>
#include <vector>

/**
 * \brief Surface traction boundary condition
 * \ingroup constraints
 */
class tledTractionConstraint : public tledSurfaceConstraint {  
  /**
   * \name Types
   * @{
   */
public:
  typedef tledSurfaceConstraint Superclass;
  /** @} */

  /**
   * \name Traction 
   * @{
   */
private:
  std::vector<float> m_FaceTractions, m_CurrentNodeTractions;

private:
  template <class TSurface>
  void _UpdateNodeForces(const int step, const double dt, const double T);

public:
  virtual void SetFaces(const int type, const std::vector<int> &faces);
  void SetFaceTractions(const std::vector<float> &ts);

  /** Peak (amplitude = 1) traction for given facet index. */
  const float* GetFacetPeakTraction(const int fInd) const { return &m_FaceTractions[3*fInd]; }
  const std::vector<float>& GetFaceTractions(void) const { return m_FaceTractions; }
  /** @} */

  /**
   * \name Constraint Manager Interface
   * @{
   */
public:
  virtual std::vector<float>* GetForceVal(int dof, int step, double dt, double T);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  /** See SetSolver, SetFaces, SetFaceTractions, SetLoadShape */
  tledTractionConstraint(tledSolver &r_solver, const int faceType, const std::vector<int> &faces, const std::vector<float> &tractions, const loadShape loadShape);
  virtual ~tledTractionConstraint(void) {}
  /** @} */
};


#endif
