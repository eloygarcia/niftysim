// =========================================================================
// File:       tledConstraintManager.h
// Purpose:    Class for managing all constraints in an analysis
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


#ifndef tledConstraintManager_H
#define tledConstraintManager_H

#include "tledConstraint.h"
#include "tledDispConstraint.h"
#include "tledForceConstraint.h"
#include "tledFixConstraint.h"
#include "tledGravityConstraint.h"
#include "tledPressureConstraint.h"
#include "tledTractionConstraint.h"
#include "tledModel.h"
#include "tledSolver.h"
#include "xmlParser.h"

#include <fstream>

/**
 * \brief Boundary condition manager
 * \ingroup solver
 * \ingroup constraints
 */
class tledConstraintManager
{
  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledConstraintManager();
  /** Construction and initialisation */
  tledConstraintManager(tledModel* Model, tledSolver* Solver);  
  virtual ~tledConstraintManager();
  /** @} */

  /**
   * \name Constraint Queries
   * @{
   */
public:
   // Methods required by Solver
   std::vector<int>* GetDispInd(int dof);   // Displaced node indices
   /** Cannot be called before an initialising call to tledConstraintManager::GetDispInd. */
   std::vector<float>* GetDispVal(int dof, int step, double dt, double T); // Displacements
   std::vector<int>* GetForceInd(int dof);  // Forced node indices
   /** Cannot be called before an initialising call to tledConstraintManager::GetForceInd. */
   std::vector<float>* GetForceVal(int dof, int step, double dt, double T);  // Forces
   std::vector<int>* GetFixInd(int dof);    // Fixed node indices
   
   // Other utilities
   int GetNumConstraints() const {return Constraints.size() + FixConstraints.size();}
   int GetNumDispConstraints() const;
   int GetNumForceConstraints() const;
   int GetNumFixConstraints() const {return FixConstraints.size();}
   tledConstraint* GetConstraint(int num);
   const tledConstraint* GetConstraint(int num) const;
   tledFixConstraint* GetFixConstraint(int num);
   const tledFixConstraint* GetFixConstraint(int num) const;
  /** @} */

  /**
   * \name Setup
   * @{
   */
protected:
  /** Create constraint of type cType corresponding the constraint with index cInd in the model. */
  virtual tledConstraint* CreateConstraint(const std::string &cType, const int cInd, const tledModel &model);

public:
   void AddConstraint(tledConstraint* constraint);
   void AddFixConstraint(tledFixConstraint* constraint);

  /** Assembles simulation constraints by querying "model" */
  virtual void Init(tledSolver &r_solver, const tledModel &model);
  /** @} */

private:
   std::vector<tledConstraint*> Constraints;
   std::vector<tledFixConstraint*> FixConstraints;
   std::vector<int> IndDisp[3];
   std::vector<int> IndForce[3];
   std::vector<int> IndFix[3];
   std::vector<float> ValDisp[3];
   std::vector<float> ValForce[3];
};


#endif // tledConstraintManager_H
