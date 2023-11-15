// =========================================================================
// File:       tledMaterialPY.h
// Purpose:    Material class - Polynomial (N=2)
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

#ifndef tledMaterialPY_H
#define tledMaterialPY_H

#include "tledMaterial.h"

/**
 * \brief Polynomial (N=2) constitutive model
 * \ingroup constitutive
 */
class tledMaterialPY : public tledMaterial
{
public:
   tledMaterialPY();
   tledMaterialPY(float* MatParams, float Density);
   virtual ~tledMaterialPY();
   
   void ComputeSPKStress(float X[3][3], float* SPK);
   void ComputeSPKStress(float X[3][3], float SPK[3][3]);
   void GetHGLameParams(float* Lambda, float* Mu);
   void GetKappa(float* kappa){*kappa = K;}
   void SetTimeStep(float dt) {;}
   void SetParams(std::vector<float> params);
   
private:
   float c10, c01, c20, c02, c11; // Params
   float K;	// Initial bulk modulus
   float M; // Initial shear modulus
};

#endif // tledMaterialPY_H
