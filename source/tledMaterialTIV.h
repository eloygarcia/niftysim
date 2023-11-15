// =========================================================================
// File:       tledMaterialTIV.h
// Purpose:    Material class - Transversely isotropic neo-Hookean viscoelastic
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


#ifndef tledMaterialTIV_H
#define tledMaterialTIV_H

#include "tledMaterial.h"

/**
 * \brief Transversely isotropic neo-Hookean viscoelastic constitutive model
 * \ingroup constitutive
 */
class tledMaterialTIV : public tledMaterial
{
public:
   tledMaterialTIV();
   tledMaterialTIV(float* MatParams, float Density);
   virtual ~tledMaterialTIV();

   void ComputeSPKStress(float X[3][3], float* SPK);
   void ComputeSPKStress(float X[3][3], float SPK[3][3]);
   void GetHGLameParams(float* Lambda, float* Mu);
   void GetKappa(float* kappa){*kappa = K;}
   void SetTimeStep(float dt);
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
   int Ni;	// Number of Prony terms
   int Nv;
   float* Ai;	// Backward Euler params: Ai1,Bi1,Ai2,Bi2,...,AiN,BiN --> length = 2Ni
   float* Av;
   float* Di;	// State variables: ([Di1,1-6],[Di2,1-6],...,[DiN,1-6]) --> length = 6Ni
   float* Dv;
   float Dt; // Time step
};

#endif // tledMaterialTIV_H
