// =========================================================================
// File:       tledContactCylinder.cpp
// Purpose:    Create instance of a rigid cylinder for contact (CPU version)
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

#include "tledContactCylinder.h"

using namespace std;

tledContactCylinder::tledContactCylinder(vector<float> orig, vector<float> axis, float R, float L, vector<int> slvs, vector<float> origdisp, float radchng, int NumNodes)
{
   // CPU version not implemented
}

void tledContactCylinder::Update(double TR)
{
   // CPU version not implemented
}

vector<float> tledContactCylinder::GetStartOriginV()
{
   // CPU version not implemented
   vector<float> emptyVector;
   return emptyVector;
}

vector<float> tledContactCylinder::GetStartAxisV()
{
   // CPU version not implemented
   vector<float> emptyVector;
   return emptyVector;
}

vector<float> tledContactCylinder::GetOriginDispV()
{
   // CPU version not implemented
   vector<float> emptyVector;
   return emptyVector;
}

void tledContactCylinder::SetOriginDisp(vector<float> disp)
{
   // CPU version not implemented
}

void tledContactCylinder::SetStartOrigin(vector<float> orig)
{
   // CPU version not implemented
}

void tledContactCylinder::SetStartAxis(vector<float> axis)
{
   // CPU version not implemented
}


#endif // _GPU_
