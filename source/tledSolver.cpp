// =========================================================================
// File:       tledSolver.cpp
// Purpose:    Main finite element object base class
// Package:    NiftySim: Nonlinear finite element program
// Author:     Zeike Taylor
// Language:   C++
// Created:    July 2007
// 
// Copyright (c) 2010, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// ztaylor@itee.uq.edu.au
// =========================================================================


#include "tledSolver.h"
#include "tledShellSolver.h"

using namespace std;

tledSolver::~tledSolver() {
  if (this->HasMembrane()) delete mp_ShellSolver;
}

void tledSolver::PrintKineticEnergy() {
  cout.setf(ios::fixed, ios::floatfield);
  cout << "E kinetic:       " << GetKineticEnergy() << endl;
}

void tledSolver::PrintStrainEnergy() {
  cout.setf(ios::fixed, ios::floatfield);
  cout << "E strain:        " << GetStrainEnergy() << endl;
}
