// =========================================================================
// File:       tledDispConstraint.h
// Purpose:    Displacement constraint class. Used for applying nodal displacement loads.
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


#ifndef tledDispConstraint_H
#define tledDispConstraint_H

#include "tledConstraint.h"

/**
 * \brief Non-zero displacement boundary conditions
 * \ingroup constraints
 */
class tledDispConstraint : public tledConstraint
{
public:
   tledDispConstraint();
   tledDispConstraint(std::vector<int> indices, std::vector<float> magnitudes, enum loadShape ls, int dof);
   virtual ~tledDispConstraint();

  /**
   * \name Get params
   * @{
   */
public:
  /** Affected degree of freedom (0,1,2) */
  int GetDOF() {return DOF;}
  /** Node indices */
  const std::vector<int>& GetInd() const {return Ind;}
  /** Node indices */
  std::vector<int>& GetInd() {return Ind;}
  std::vector<float> GetMag() {return Mag;}
  /** @} */
   
  /**
   * \name Set params
   * @{
   */
public:
  void SetDOF(int dof) {DOF = dof;}   
  void SetInd(std::vector<int> ind) {Ind = ind;}
  void SetMag(std::vector<float> mag) {Mag = mag;}
  /** @} */
   
  /**
   * \name Methods required by tledConstraintManager
   * @{
   */
public:
  virtual std::vector<int>* GetDispInd(int dof);
  /** Displacements for current step */
  virtual std::vector<float>* GetDispVal(int dof, int step, double dt, double T);
  /** @} */

private:
   int DOF; 
   std::vector<int> Ind;  
   std::vector<float> Mag;
   std::vector<float> U;  // Displacements for current step
};


#endif // tledDispConstraint_H
