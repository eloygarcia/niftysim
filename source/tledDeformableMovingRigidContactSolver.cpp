// =========================================================================
// File:       tledDeformableMovingRigidContactSolver.cpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    December 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#include "tledDeformableMovingRigidContactSolver.h"
#include "tledDeformableMovingRigidContactSolverCPU.h"

#include <iostream>
#include <cstdlib>

tledDeformableRigidContactSolver* tledDeformableMovingRigidContactSolver::CreateContactSolver(tledUnstructuredContactManager &r_manager, const int surfaceIndex) {
  if (r_manager.UseGPU()) {
    return NULL;
  } else {
    return tledDeformableMovingRigidContactSolverCPU::CreateContactSolver(r_manager, surfaceIndex);
  }
}
