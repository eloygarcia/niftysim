// =========================================================================
// File:       tledDispConstraint.cpp
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


#include "tledDispConstraint.h"

using namespace std;

tledDispConstraint::tledDispConstraint()
{
}

tledDispConstraint::tledDispConstraint(vector<int> indices, vector<float> magnitudes, enum loadShape ls, int dof)
{
   // Check inputs
  if (indices.size() != magnitudes.size()) {
     tledLogErrorStream(tledHelper::Warning() << "Vectors 'indices' and 'magnitudes' must be of equal length");
  }
   
   // Assign inputs
   DOF = dof;
   this->SetLoadShape(ls);
   Ind = indices;
   Mag = magnitudes;
}

tledDispConstraint::~tledDispConstraint()
{
   Ind.resize(0);
   Mag.resize(0);
   U.resize(0);
}

vector<int>* tledDispConstraint::GetDispInd(int dof)
{
   if (dof != DOF)
      return &tledConstraint::emptyIntVec;
   
   return &Ind;
}

vector<float>* tledDispConstraint::GetDispVal(int dof, int step, double dt, double T)
{
   if (dof != DOF)
      return &tledConstraint::emptyFloatVec;
   
   // Elapsed relative (i.e. in [0,1]) simulation time
   double TR = dt*(step+1)/T;
   // Compute current vals
   float Amp = this->ComputeAmplitude(TR, this->GetLoadShape());
   if (U.size() != Ind.size())
      U.resize(Ind.size());
   for (int node = 0; node < (int)Ind.size(); node++)
      U[node] = Amp*Mag[node];
   
   return &U;
}

