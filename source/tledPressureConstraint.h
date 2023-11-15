// =========================================================================
// File:       tledPressureConstraint.h
// Purpose:    Pressure constraint class. Used for applying uniform surface
//             pressure loads.
// Package:    NiftySim: Nonlinear finite element program
// Author:     Zeike Taylor
// Language:   C++
// Created:    March 2011
// 
// Copyright (c) 2011, University of Queensland. All rights reserved.
// MedTeQ Centre
// See the LICENSE.txt file in the root folder
// 
// ztaylor@itee.uq.edu.au
// =========================================================================


#ifndef tledPressureConstraint_H
#define tledPressureConstraint_H

#include "tledSurfaceConstraint.h"

/**
 * \brief Pressure boundary condition
 * \ingroup constraints
 */
class tledPressureConstraint : public tledSurfaceConstraint {   
  /**
   * \name Types
   * @{
   */
public:
  typedef tledSurfaceConstraint Superclass;
  /** @} */   

  /**
   * \name Pressure Force
   * @{
   */
private:
  float m_Mag;
  std::vector<float> m_CurrentNodePressureForces;

private:
  template <class TSurface>
  void _UpdateNodePressureForces(const int step, const double dt, const double T);

public:
  /** (Peak) pressure magnitude */
  float GetMag(void) const { return m_Mag; }  
  void SetMag(const float mag) { m_Mag = mag;}
  /** @} */

  /**
   * \name Geometry
   * @{
   */
public:
  virtual void SetFaces(const int type, const std::vector<int> &faces);
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
  tledPressureConstraint(void);
  virtual ~tledPressureConstraint();
  
  /** Legacy constructor signature */
  tledPressureConstraint(std::vector<int> faces, int facetype, tledSolver* solver, float mag, enum loadShape ls);

  /** See SetSolver, SetFaces, SetMag, SetLoadShape */
  tledPressureConstraint(tledSolver &r_solver, const int faceType, const std::vector<int> &faces, const float mag, const loadShape loadShape);      
  /** @} */
};

#endif // tledPressureConstraint_H
