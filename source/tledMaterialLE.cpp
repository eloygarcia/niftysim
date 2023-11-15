// =========================================================================
// File:       tledMaterialLE.cpp
// Purpose:    Material class - linear elastic
// Package:    NiftySim: Nonlinear finite element program
// Author:     Zeike Taylor
// Language:   C++
// Created:    March 2007
// 
// Copyright (c) 2010, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// ztaylor@itee.uq.edu.au
// =========================================================================


#include "tledMaterialLE.h"
#include "tledMatrixFunctions.h"
#include <iostream>

using namespace std;

tledMaterialLE::tledMaterialLE()
{
}

tledMaterialLE::tledMaterialLE(float* MatParams, float Density) : tledMaterial(Density)
{
   E = MatParams[0];
   P = MatParams[1];

   c1 = E*(1 - P)/(2*(1 + P)*(1 - 2*P));
   c2 = c1*P/(1 - P);
   c3 = c1 + 2*c2;
   c4 = E/2/(1+P);
}

tledMaterialLE::~tledMaterialLE()
{
}

void tledMaterialLE::GetHGLameParams(float* Lambda, float* Mu)
{
   *Lambda = E*P/((1+P)*(1-2*P));
   *Mu = E/(2*(1+P));
}

void tledMaterialLE::SetParams(vector<float> params)
{
   E = params[0];
   P = params[1];

   c1 = E*(1 - P)/(2*(1 + P)*(1 - 2*P));
   c2 = c1*P/(1 - P);
   c3 = c1 + 2*c2;
   c4 = E/2/(1+P);
}

void tledMaterialLE::ComputeSPKStress(float X[3][3], float *SPK)
{
   float XT11,XT12,XT13,XT21,XT22,XT23,XT31,XT32,XT33;
   float C11,C12,C13,C22,C23,C33;

   // Transpose of deformation gradient
   XT11 = X[0][0]; XT12 = X[1][0]; XT13 = X[2][0];
   XT21 = X[0][1]; XT22 = X[1][1]; XT23 = X[2][1];
   XT31 = X[0][2]; XT32 = X[1][2]; XT33 = X[2][2];

   // Right Cauchy-Green deformation tensor
   C11 = XT11*XT11 + XT12*XT12 + XT13*XT13;
   C12 = XT11*XT21 + XT12*XT22 + XT13*XT23;
   C13 = XT11*XT31 + XT12*XT32 + XT13*XT33;
   C22 = XT21*XT21 + XT22*XT22 + XT23*XT23;
   C23 = XT21*XT31 + XT22*XT32 + XT23*XT33;
   C33 = XT31*XT31 + XT32*XT32 + XT33*XT33;

   //// 2nd Piola-Kirchhoff stress
   SPK[0] = c1*C11 + c2*(C22 + C33) - c3;
   SPK[1] = c1*C22 + c2*(C11 + C33) - c3;
   SPK[2] = c1*C33 + c2*(C11 + C22) - c3;
   SPK[3] = c4*C12;
   SPK[4] = c4*C23;
   SPK[5] = c4*C13;
   
//    cout << "SPK = " << SPK[0] << " " << SPK[1] << " " << SPK[2] << " " << SPK[3] << " " << SPK[4] << " " << SPK[5] << endl;
}

void tledMaterialLE::ComputeSPKStress(float X[3][3], float SPK[3][3])
{
   float XT11,XT12,XT13,XT21,XT22,XT23,XT31,XT32,XT33;
   float C11,C12,C13,C22,C23,C33;

   // Transpose of deformation gradient
   XT11 = X[0][0]; XT12 = X[1][0]; XT13 = X[2][0];
   XT21 = X[0][1]; XT22 = X[1][1]; XT23 = X[2][1];
   XT31 = X[0][2]; XT32 = X[1][2]; XT33 = X[2][2];

   // Right Cauchy-Green deformation tensor
   C11 = XT11*XT11 + XT12*XT12 + XT13*XT13;
   C12 = XT11*XT21 + XT12*XT22 + XT13*XT23;
   C13 = XT11*XT31 + XT12*XT32 + XT13*XT33;
   C22 = XT21*XT21 + XT22*XT22 + XT23*XT23;
   C23 = XT21*XT31 + XT22*XT32 + XT23*XT33;
   C33 = XT31*XT31 + XT32*XT32 + XT33*XT33;

   SPK[0][0] = c1*C11 + c2*(C22 + C33) - c3; // S00
   SPK[1][1] = c1*C22 + c2*(C11 + C33) - c3; // S11
   SPK[2][2] = c1*C33 + c2*(C11 + C22) - c3; // S22
   SPK[0][1] = c4*C12; // S01
   SPK[1][0] = SPK[0][1];
   SPK[1][2] = c4*C23; // S12
   SPK[2][1] = SPK[1][2];
   SPK[0][2] = c4*C13; // S02
   SPK[2][0] = SPK[0][2];
   
//    cout << "SPK = " << SPK[0][0] << " " << SPK[1][1] << " " << SPK[2][2] << " " << SPK[0][1] << " " << SPK[1][2] << " " << SPK[0][2] << endl;
}
