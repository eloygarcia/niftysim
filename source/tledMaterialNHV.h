// =========================================================================
// File:       tledMaterialNHV.h
// Purpose:    Material class - Neo-Hookean viscoelastic
// Package:    NiftySim: Nonlinear finite element program
// Author:     Zeike Taylor
// Language:   C++
// Created:    September 2007
// 
// Copyright (c) 2010, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// ztaylor@itee.uq.edu.au
// =========================================================================


#ifndef tledMaterialNHV_H
#define tledMaterialNHV_H

#include "tledMaterial.h"

/**
 * \brief Neo-Hookean viscoelastic constitutive model
 * \ingroup constitutive
 */
class tledMaterialNHV : public tledMaterial
{
public:
   tledMaterialNHV();
   tledMaterialNHV(float* MatParams, float Density);
   virtual ~tledMaterialNHV();

   void ComputeSPKStress(float X[3][3], float* SPK);
   void ComputeSPKStress(float X[3][3], float SPK[3][3]);
   void GetHGLameParams(float* Lambda, float* Mu);
   void GetKappa(float* kappa){*kappa = K;}
   void SetTimeStep(float dt);
   void SetParams(std::vector<float> params);

private:
   float K;	// Bulk modulus
   float M;
   int Ni;	// Number of Prony terms
   int Nv;
   float* Ai;	// Backward Euler params: Ai1,Bi1,Ai2,Bi2,...,AiN,BiN --> length = 2Ni
   float* Av;
   float* Di;	// State variables: ([Di1,1-6],[Di2,1-6],...,[DiN,1-6]) --> length = 6Ni
   float* Dv;
   float Dt; // Time step
};

#endif // tledMaterialNHV_H
