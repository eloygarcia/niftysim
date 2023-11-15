// =========================================================================
// File:       tledFixConstraint.cpp
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


#include "tledFixConstraint.h"

using namespace std;

tledFixConstraint::tledFixConstraint()
{
}

tledFixConstraint::tledFixConstraint(vector<int> indices, int dof)
{
   // Assign inputs
   DOF = dof;
   Ind = indices;
}

tledFixConstraint::~tledFixConstraint()
{
   Ind.resize(0);
}

vector<int>* tledFixConstraint::GetDispInd(int dof)
{
   if (dof != DOF)
     return &tledConstraint::emptyIntVec;
   
   return &Ind;
}

vector<int>* tledFixConstraint::GetFixInd(int dof)
{
   if (dof != DOF)
     return &tledConstraint::emptyIntVec;

   return &Ind;
}
