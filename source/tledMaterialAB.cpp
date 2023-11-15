// =========================================================================
// File:       tledMaterialAB.cpp
// Purpose:    Material class - Arruda-Boyce
// Package:    NiftySim: Nonlinear finite element program
// Author:     Zeike Taylor
// Language:   C++
// Created:    October 2010
// 
// Copyright (c) 2010, University of Queensland. All rights reserved.
// MedTeQ Centre
// See the LICENSE.txt file in the root folder
// 
// ztaylor@itee.uq.edu.au
// =========================================================================

#include "tledMaterialAB.h"
#include "tledMatrixFunctions.h"

#include <iostream>
#include <cmath>

using namespace std;

tledMaterialAB::tledMaterialAB()
{
}

tledMaterialAB::tledMaterialAB(float* MatParams, float Density) : tledMaterial(Density)
{      
   M = MatParams[0]; // Initial shear modulus
   Lm = MatParams[1]; // Locking stretch
   K = MatParams[2]; // Initial bulk modulus
}

tledMaterialAB::~tledMaterialAB()
{
}

void tledMaterialAB::GetHGLameParams(float* Lambda, float* Mu)
{
   *Lambda = K - 2*M/3;
   *Mu = M;
}

void tledMaterialAB::SetParams(vector<float> params)
{
   M = params[0]; // Initial shear modulus
   Lm = params[1]; // Locking stretch
   K = params[2]; // Initial bulk modulus
}

void tledMaterialAB::ComputeSPKStress(float X[3][3], float *SPK)
{
   // Transpose of deformation gradient
   float XT11 = X[0][0]; float XT12 = X[1][0]; float XT13 = X[2][0];
   float XT21 = X[0][1]; float XT22 = X[1][1]; float XT23 = X[2][1];
   float XT31 = X[0][2]; float XT32 = X[1][2]; float XT33 = X[2][2];
   
   // Right Cauchy-Green deformation tensor
   float C11 = XT11*XT11 + XT12*XT12 + XT13*XT13;
   float C12 = XT11*XT21 + XT12*XT22 + XT13*XT23;
   float C13 = XT11*XT31 + XT12*XT32 + XT13*XT33;
   float C22 = XT21*XT21 + XT22*XT22 + XT23*XT23;
   float C23 = XT21*XT31 + XT22*XT32 + XT23*XT33;
   float C33 = XT31*XT31 + XT32*XT32 + XT33*XT33;
   
   // Determinant of deformation gradient
   float J = XT11*(XT22*XT33 - XT32*XT23) + XT12*(XT23*XT31 - XT21*XT33) + XT13*(XT21*XT32 - XT22*XT31);
   // Determinant of C
   float detC = C11*(C22*C33 - C23*C23) - C22*C13*C13 - C12*(C33*C12 - 2*C13*C23);
   
   float I1b = (float)pow((double)J,-(double)2/(double)3)*(C11+C22+C33);
   // c1 = 1/2; c2 = 1/20; c3 = 11/1050; c4 = 19/7000; c5 = 519/673750;
   float Lm2 = Lm*Lm;
   float c1 = 0.5f;
   float c2 = 0.05f; c2 *= 2*I1b/Lm2;
   float I1b2 = I1b*I1b;
   Lm2 *= Lm*Lm; // Lm2 = Lm^4
   float c3 = 0.01047619047619f; c3 *= 3*I1b2/Lm2;
   I1b2 *= I1b; // I1b2 = I1b^3
   Lm2 *= Lm*Lm; // Lm2 = Lm^6
   float c4 = 0.002714285714286f; c4 *= 4*I1b2/Lm2;
   I1b2 *= I1b; // I1b2 = I1b^4
   Lm2 *= Lm*Lm; // Lm2 = Lm^8
   float c5 = 7.703153988868275e-4f; c5 *= 5*I1b2/Lm2;
   float x1 = (float)pow((double)J,-(double)2/(double)3)*2*M*(c1 + c2 + c3 + c4 + c5);
   float x2 = (float)((0.5*K*J*(J - 1/J) - x1*(C11+C22+C33)/3)/detC);
   
   SPK[0] = (C22*C33 - C23*C23)*x2 + x1;
   SPK[1] = (C11*C33 - C13*C13)*x2 + x1;
   SPK[2] = (C11*C22 - C12*C12)*x2 + x1;
   SPK[3] = (C13*C23 - C12*C33)*x2;
   SPK[4] = (C12*C13 - C23*C11)*x2;
   SPK[5] = (C12*C23 - C13*C22)*x2;
}

void tledMaterialAB::ComputeSPKStress(float X[3][3], float SPK[3][3])
{
   // Transpose of deformation gradient
   float XT11 = X[0][0]; float XT12 = X[1][0]; float XT13 = X[2][0];
   float XT21 = X[0][1]; float XT22 = X[1][1]; float XT23 = X[2][1];
   float XT31 = X[0][2]; float XT32 = X[1][2]; float XT33 = X[2][2];
   
   // Right Cauchy-Green deformation tensor
   float C11 = XT11*XT11 + XT12*XT12 + XT13*XT13;
   float C12 = XT11*XT21 + XT12*XT22 + XT13*XT23;
   float C13 = XT11*XT31 + XT12*XT32 + XT13*XT33;
   float C22 = XT21*XT21 + XT22*XT22 + XT23*XT23;
   float C23 = XT21*XT31 + XT22*XT32 + XT23*XT33;
   float C33 = XT31*XT31 + XT32*XT32 + XT33*XT33;
   
   // Determinant of deformation gradient
   float J = XT11*(XT22*XT33 - XT32*XT23) + XT12*(XT23*XT31 - XT21*XT33) + XT13*(XT21*XT32 - XT22*XT31);
   // Determinant of C
   float detC = C11*(C22*C33 - C23*C23) - C22*C13*C13 - C12*(C33*C12 - 2*C13*C23);
   
   float I1b = (float)pow((double)J,-(double)2/(double)3)*(C11+C22+C33);
   // c1 = 1/2; c2 = 1/20; c3 = 11/1050; c4 = 19/7000' c5 = 519/673750
   float Lm2 = Lm*Lm;
   float c1 = 0.5f;
   float c2 = 0.05f; c2 *= 2*I1b/Lm2;
   float I1b2 = I1b*I1b;
   Lm2 *= Lm*Lm; // Lm2 = Lm^4
   float c3 = 0.01047619047619f; c3 *= 3*I1b2/Lm2;
   I1b2 *= I1b; // I1b2 = I1b^3
   Lm2 *= Lm*Lm; // Lm2 = Lm^6
   float c4 = 0.002714285714286f; c4 *= 4*I1b2/Lm2;
   I1b2 *= I1b; // I1b2 = I1b^4
   Lm2 *= Lm*Lm; // Lm2 = Lm^8
   float c5 = 7.703153988868275e-4f; c5 *= 5*I1b2/Lm2;
   float x1 = (float)pow((double)J,-(double)2/(double)3)*2*M*(c1 + c2 + c3 + c4 + c5);
   float x2 = (float)((0.5*K*J*(J - 1/J) - x1*(C11+C22+C33)/3)/detC);
   
   SPK[0][0] = (C22*C33 - C23*C23)*x2 + x1;
   SPK[1][1] = (C11*C33 - C13*C13)*x2 + x1;
   SPK[2][2] = (C11*C22 - C12*C12)*x2 + x1;
   SPK[0][1] = (C13*C23 - C12*C33)*x2;
   SPK[1][0] = SPK[0][1];
   SPK[1][2] = (C12*C13 - C23*C11)*x2;
   SPK[2][1] = SPK[1][2];
   SPK[0][2] = (C12*C23 - C13*C22)*x2;
   SPK[2][0] = SPK[0][2];
   
//   cout << "SPK = " << SPK[0][0] << " " << SPK[1][1] << " " << SPK[2][2] << endl;
}
