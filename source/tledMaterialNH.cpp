// =========================================================================
// File:       tledMaterialNH.cpp
// Purpose:    Material class - Neo-Hookean
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

#include "tledMaterialNH.h"
#include "tledMatrixFunctions.h"

#include <iostream>
#include <cmath>

using namespace std;

tledMaterialNH::tledMaterialNH()
{
}

tledMaterialNH::tledMaterialNH(float* MatParams, float Density) : tledMaterial(Density)
{
   // MatParams[0] = M, MatParams[1] = K
   M = MatParams[0];
   K = MatParams[1];
}

tledMaterialNH::~tledMaterialNH()
{
}

void tledMaterialNH::GetHGLameParams(float* Lambda, float* Mu)
{
   *Lambda = K - 2*M/3;
   *Mu = M;
}

void tledMaterialNH::SetParams(vector<float> params)
{
   M = params[0];
   K = params[1];
}

void tledMaterialNH::ComputeSPKStress(float X[3][3], float *SPK)
{
   float J,detC,x1,x2;
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

   // Determinant of deformation gradient
   J = XT11*(XT22*XT33 - XT32*XT23) + XT12*(XT23*XT31 - XT21*XT33) + XT13*(XT21*XT32 - XT22*XT31);
   // Determinant of C
   detC = C11*(C22*C33 - C23*C23) - C22*C13*C13 - C12*(C33*C12 - 2*C13*C23);

   x1 = (float)pow((double)J,-(double)2/(double)3)*M;
   x2 = (K*J*(J-1) - x1*(C11+C22+C33)/3)/detC;

   SPK[0] = (C22*C33 - C23*C23)*x2 + x1;
   SPK[1] = (C11*C33 - C13*C13)*x2 + x1;
   SPK[2] = (C11*C22 - C12*C12)*x2 + x1;
   SPK[3] = (C13*C23 - C12*C33)*x2;
   SPK[4] = (C12*C13 - C23*C11)*x2;
   SPK[5] = (C12*C23 - C13*C22)*x2;
}

void tledMaterialNH::ComputeSPKStress(float X[3][3], float SPK[3][3])
{
   float J,detC,x1,x2;
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

   // Determinant of deformation gradient
   J = XT11*(XT22*XT33 - XT32*XT23) + XT12*(XT23*XT31 - XT21*XT33) + XT13*(XT21*XT32 - XT22*XT31);
   // Determinant of C
   detC = C11*(C22*C33 - C23*C23) - C22*C13*C13 - C12*(C33*C12 - 2*C13*C23);

   x1 = (float)pow((double)J,-(double)2/(double)3)*M;
   x2 = (K*J*(J-1) - x1*(C11+C22+C33)/3)/detC;

   // 2nd Piola-Kirchhoff stress
   SPK[0][0] = (C22*C33 - C23*C23)*x2 + x1; // S00
   SPK[1][1] = (C11*C33 - C13*C13)*x2 + x1; // S11
   SPK[2][2] = (C11*C22 - C12*C12)*x2 + x1; // S22
   SPK[0][1] = (C13*C23 - C12*C33)*x2; // S01
   SPK[1][0] = SPK[0][1];
   SPK[1][2] = (C12*C13 - C23*C11)*x2; // S12
   SPK[2][1] = SPK[1][2];
   SPK[0][2] = (C12*C23 - C13*C22)*x2; // S02
   SPK[2][0] = SPK[0][2];
}
