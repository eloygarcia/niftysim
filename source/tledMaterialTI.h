// =========================================================================
// File:       tledMaterialTI.h
// Purpose:    Material class - Transversely isotropic
// Package:    NiftySim: Nonlinear finite element program
// Author:     Zeike Taylor
// Language:   C++
// Created:    January 2008
// 
// Copyright (c) 2010, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// ztaylor@itee.uq.edu.au
// =========================================================================


#ifndef tledMaterialTI_H
#define tledMaterialTI_H

#include "tledMaterial.h"


class tledMaterialTI : public tledMaterial
{
public:
   tledMaterialTI();
   tledMaterialTI(float* MatParams, float Density);
   virtual ~tledMaterialTI();

   void ComputeSPKStress(float X[3][3], float* SPK);
   void ComputeSPKStress(float X[3][3], float SPK[3][3]);
   void GetHGLameParams(float* Lambda, float* Mu);
   void GetKappa(float* kappa){*kappa = K;}
   void SetTimeStep(float dt) {;}
   void SetParams(std::vector<float> params);

private:
   float K;	// Bulk modulus
   float M;
   float E;	// Material parameter eta
   float A11;	// Structure tensor components (symmetric, so only need 6 components)
   float A22;
   float A33;
   float A12;
   float A23;
   float A13;
};

#endif // tledMaterialTI_H
