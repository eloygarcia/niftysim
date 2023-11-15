// =========================================================================
// File:       tledContactPlate.cpp
// Purpose:    Create instance of a rigid rectangular sheet for contact modelling
//             (CPU version)
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

#include "tledContactPlate.h"

using namespace std;

tledContactPlate::tledContactPlate(vector<float> A, vector<float> B, vector<float> C, vector<int> slvs, vector<float> disp, int NumNodes)
{
   // CPU version not implemented
}

void tledContactPlate::Update(double TR)
{
   // CPU version not implemented
}

vector<float> tledContactPlate::GetStartCrnrAV()
{
   // CPU version not implemented
   vector<float> emptyVector;
   return emptyVector;
}

vector<float> tledContactPlate::GetStartCrnrBV()
{
   // CPU version not implemented
   vector<float> emptyVector;
   return emptyVector;
}

vector<float> tledContactPlate::GetStartCrnrCV()
{
   // CPU version not implemented
   vector<float> emptyVector;
   return emptyVector;
}

vector<float> tledContactPlate::GetDispV()
{
   // CPU version not implemented
   vector<float> emptyVector;
   return emptyVector;
}

void tledContactPlate::SetDisp(vector<float> disp)
{
   // CPU version not implemented
}

void tledContactPlate::SetStartCrnrA(vector<float> A)
{
   // CPU version not implemented
}

void tledContactPlate::SetStartCrnrB(vector<float> B)
{
   // CPU version not implemented
}

void tledContactPlate::SetStartCrnrC(vector<float> C)
{
   // CPU version not implemented
}

#endif // _GPU_

