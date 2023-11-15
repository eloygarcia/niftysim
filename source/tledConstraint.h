// =========================================================================
// File:       tledConstraint.h
// Purpose:    Constraint base class
// Package:    NiftySim: Nonlinear finite element program
// Author:     Zeike Taylor
// Language:   C++
// Created:    June 2008
// 
// Copyright (c) 2010, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// ztaylor@itee.uq.edu.au
// =========================================================================


#ifndef tledConstraint_H
#define tledConstraint_H

#include "tledHelper.h"

#include <vector>
#include <iostream>

/**
 * \defgroup constraints Constraints
 *
 * Boundary conditions and loading
 */

/**
 * \brief Time profiles for constraint amplitudes
 * \ingroup constraints
 */
enum loadShape
{
   POLY345,
   RAMP,
   STEP,
   HILLY,
   RAMPHOLD
};

enum loadShape atols(const char* str);

/**
 * \brief Constraint base class
 * \ingroup constraints
 */
class tledConstraint
{
protected:
  static std::vector<int> emptyIntVec;
  static std::vector<float> emptyFloatVec;

  /**
   * \name Load-Shape
   *
   * Constraint time profile.
   * @{
   */
private:
  enum loadShape m_LS;

public:
  virtual void SetLoadShape(const enum loadShape ls) { m_LS = ls; }
  enum loadShape GetLoadShape(void) const {return m_LS;}
  /** @} */

public:
   tledConstraint() {;}
   virtual ~tledConstraint() {;}

  virtual std::vector<int>& GetInd(void) = 0;
  virtual const std::vector<int>& GetInd(void) const = 0;
   
  /**
   * \name Constraint queries
   *
   * Default action: return an empty std::vector for affected DOFs, values.
   * @{
   */
public:
  virtual std::vector<int>* GetDispInd(int dof) { return &emptyIntVec; }
  virtual std::vector<float>* GetDispVal(int dof, int step, double dt, double T) { return &emptyFloatVec; }
  virtual std::vector<int>* GetForceInd(int dof) { return &emptyIntVec; }
  virtual std::vector<float>* GetForceVal(int dof, int step, double dt, double T) { return &emptyFloatVec; }

protected:
  /** Returns for a given point in time the amplitude of the constraint (in [0, 1]). */
  virtual float ComputeAmplitude(const double TR, const enum loadShape ls) const;
  /** @} */
};


#endif // tledConstraint_H
