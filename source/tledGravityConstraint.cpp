// =========================================================================
// File:       .cpp
// Purpose:    Gravity constraint class. Used for applying gravitational
//             loads.
// Package:    NiftySim: Nonlinear finite element program
// Author:     Zeike Taylor
// Language:   C++
// Created:    March 2011
// 
// Copyright (c) 2011, University of Queensland. All rights reserved.
// MedTeQ Centre
// See the LICENSE.txt file in the root folder
// 
// ztaylor@itee.uq.edu.au
// =========================================================================


#include "tledGravityConstraint.h"

using namespace std;

tledGravityConstraint::tledGravityConstraint()
{
}

tledGravityConstraint::tledGravityConstraint(vector<int> ind, vector<float> mass, float mag, vector<float> dir, enum loadShape ls)
{
   // Check inputs
  if (dir.size() != 3) {
     tledLogErrorStream(tledHelper::FatalError() << "Vector 'dir' must have length = 3");
  }

  if (ind.size() != mass.size()) {
    tledLogErrorStream(tledHelper::FatalError() << "Vectors 'ind' and 'mass' must be same length");
  }
   
   // Assign inputs
   this->SetLoadShape(ls);
   Ind = ind;
   Mass = mass;
   AccelMag = mag;
   AccelDir = dir;   
}

tledGravityConstraint::~tledGravityConstraint()
{
   Ind.resize(0);
   R.resize(0);
}

void tledGravityConstraint::SetAccelDir(vector<float> dir)
{
   if (dir.size() != 3) {
     tledLogErrorStream(tledHelper::FatalError() << "Vector of invalid length passed");

     return;
   }
   AccelDir = dir;
}

vector<int>* tledGravityConstraint::GetForceInd(int dof)
{
   return &Ind;
}

vector<float>* tledGravityConstraint::GetForceVal(int dof, int step, double dt, double T)
{
   // Elapsed relative (i.e. in [0,1]) simulation time
   double TR = dt*(step+1)/T;
   // Compute current vals
   if (R.size() != Ind.size())
      R.resize(Ind.size());
   float Amp = this->ComputeAmplitude(TR, this->GetLoadShape());
   for (int node = 0; node < (int)Ind.size(); node++)
      R[node] = Amp*AccelDir[dof]*Mass[node]*AccelMag;
   
   return &R;
}
