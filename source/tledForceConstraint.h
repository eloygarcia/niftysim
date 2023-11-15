// =========================================================================
// File:       tledForceConstraint.h
// Purpose:    Force constraint class. Used for applying nodal force loads.
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


#include "tledConstraint.h"

#ifndef tledForceConstraint_H
#define tledForceConstraint_H

/**
 * \brief Force constraints: can be used for natural BC as well as body forces
 * \ingroup constraints
 */
class tledForceConstraint : public tledConstraint
{
public:
   tledForceConstraint();
   tledForceConstraint(std::vector<int> indices, std::vector<float> magnitudes, enum loadShape ls, int dof);
   virtual ~tledForceConstraint();

   // Get params
   int GetDOF() const {return DOF;}
  virtual std::vector<int>& GetInd() {return Ind;}
  virtual const std::vector<int>& GetInd() const {return Ind;}
  std::vector<float>& GetMag() {return Mag;}
  const std::vector<float>& GetMag() const {return Mag;}
   
   // Set params
   void SetDOF(int dof) {DOF = dof;}
   void SetInd(std::vector<int> ind) {Ind = ind;}
   void SetMag(std::vector<float> mag) {Mag = mag;}
   
   // Methods required by tledConstraintManager
   std::vector<int>* GetForceInd(int dof);
   std::vector<float>* GetForceVal(int dof, int step, double dt, double T);

private:
   int DOF; // Affected degree of freedom (0,1,2)
   std::vector<int> Ind;  // Node indices
   std::vector<float> Mag;
   std::vector<float> R;  // Forces for current step   
};


#endif // tledForceConstraint_H
