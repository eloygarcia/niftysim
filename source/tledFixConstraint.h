// =========================================================================
// File:       tledFixConstraint.h
// Purpose:    Fix constraint class. Used for fixing nodal DOFs.
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

#ifndef tledFixConstraint_H
#define tledFixConstraint_H

/**
 * \brief Zero-displacement boundary condition
 * \ingroup constraints
 */
class tledFixConstraint : public tledConstraint
{
public:
   tledFixConstraint();
   tledFixConstraint(std::vector<int> indices, int dof);
   virtual ~tledFixConstraint();

   // Get params
   int GetDOF() {return DOF;}
   virtual std::vector<int>& GetInd() {return Ind;}
   virtual const std::vector<int>& GetInd() const {return Ind;}
   
   // Set params
   void SetDOF(int dof) {DOF = dof;}
   void SetInd(std::vector<int> ind) {Ind = ind;}
   
   // Methods required by tledConstraintManager
  virtual std::vector<int>* GetDispInd(int dof);
  std::vector<int>* GetFixInd(int dof);

private:
   int DOF;	// Affected degree of freedom (0,1,2)
   std::vector<int> Ind;	// Node indices
};


#endif // tledFixConstraint_H
