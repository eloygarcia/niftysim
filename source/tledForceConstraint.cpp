// =========================================================================
// File:       tledForceConstraint.cpp
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


#include "tledForceConstraint.h"

using namespace std;

tledForceConstraint::tledForceConstraint()
{
}

tledForceConstraint::tledForceConstraint(vector<int> indices, vector<float> magnitudes, enum loadShape ls, int dof)
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

tledForceConstraint::~tledForceConstraint()
{
   Ind.resize(0);
   Mag.resize(0);
   R.resize(0);
}

vector<int>* tledForceConstraint::GetForceInd(int dof)
{
   if (dof != DOF)
      return &tledConstraint::emptyIntVec;
   
   return &Ind;
}

vector<float>* tledForceConstraint::GetForceVal(int dof, int step, double dt, double T)
{
   if (dof != DOF)
     return &tledConstraint::emptyFloatVec;
   
   // Elapsed relative (i.e. in [0,1]) simulation time
   double TR = dt*(step+1)/T;
   // Compute current vals
   if (R.size() != Ind.size())
      R.resize(Ind.size());
   float Amp = this->ComputeAmplitude(TR, this->GetLoadShape());
   for (int node = 0; node < (int)Ind.size(); node++)
      R[node] = Amp*Mag[node];
   
   return &R;
}
