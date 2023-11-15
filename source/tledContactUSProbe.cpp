// =========================================================================
// File:       tledContactUSProbe.cpp
// Purpose:    Create instance of a rigid ultrasound probe for contact (CPU
//             version)
// Package:    NiftySim: Nonlinear finite element program
// Author:     Zeike Taylor
// Language:   C++
// Created:    April 2011
// 
// Copyright (c) 2011, University of Queensland. All rights reserved.
// MedTeQ Centre
// See the LICENSE.txt file in the root folder
// 
// ztaylor@itee.uq.edu.au
// =========================================================================

#ifndef _GPU_

#include "tledContactUSProbe.h"

using namespace std;

tledContactUSProbe::tledContactUSProbe(vector<float> orig, vector<float> axis, float R, float L, vector<int> slvs, vector<float> origdisp, float radchng, int NumNodes)
{
   // CPU version not implemented
}

void tledContactUSProbe::Update(double TR)
{
   // CPU version not implemented
}

vector<float> tledContactUSProbe::GetStartOriginV()
{
   // CPU version not implemented
   vector<float> emptyVector;
   return emptyVector;
}

vector<float> tledContactUSProbe::GetStartAxisV()
{
   // CPU version not implemented
   vector<float> emptyVector;
   return emptyVector;
}

vector<float> tledContactUSProbe::GetOriginDispV()
{
   // CPU version not implemented
   vector<float> emptyVector;
   return emptyVector;
}

void tledContactUSProbe::SetOriginDisp(vector<float> disp)
{
   // CPU version not implemented
}

void tledContactUSProbe::SetStartOrigin(vector<float> orig)
{
   // CPU version not implemented
}

void tledContactUSProbe::SetStartAxis(vector<float> axis)
{
   // CPU version not implemented
}

#endif // _GPU_

