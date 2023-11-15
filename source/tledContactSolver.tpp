// =========================================================================
// File:       tledContactSolver.tpp
// Purpose:    Contact solver interface and general purpose functions
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

template <class TContactMesh, class TContactSolverInterface>
tledContactSolverImpl<TContactMesh, TContactSolverInterface>::tledContactSolverImpl(tledUnstructuredContactManager &r_manager) : mr_Mesh(r_manager.GetDeformableSurface<TContactMesh>()), mr_Manager(r_manager) {
  assert(this->GetNodeCloseDistance() >= 0);
}

template <class TContactMesh, class TContactSolverInterface>
void tledContactSolverImpl<TContactMesh, TContactSolverInterface>::Init() {  
}

