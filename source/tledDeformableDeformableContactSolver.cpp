// =========================================================================
// File:       tledDeformableDeformableContactSolver.cpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    May 2011
// 
// Copyright (c) 2011, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#include "tledDeformableContactSurface.h"
#include "tledDeformableDeformableContactSolver.h"
#include "tledDeformableDeformableContactSolverCPU.h"
#ifdef GPU_GP_CONTACT
#include "tledDeformableDeformableContactSolverGPU.h"
#endif

#include <cstdlib>
#include <iostream>

tledDeformableDeformableContactSolver* tledDeformableDeformableContactSolver::CreateContactSolver(tledUnstructuredContactManager &r_manager) {
  if (r_manager.UseGPU()) {
    tledDeformableDeformableContactSolver *p_solver = NULL;

#ifdef GPU_GP_CONTACT
    p_solver = tledDeformableDeformableContactSolverGPU::CreateContactSolver(r_manager);
#else
    std::cerr << "In tledDeformableDeformableContactSolver::CreateContactSolver: Experimental contact-modelling features are not enabled in this build.\n";
    std::abort();
#endif

    return p_solver;
  } else {
    return tledDeformableDeformableContactSolverCPU::CreateContactSolver(r_manager);
  }
}
