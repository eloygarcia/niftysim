// =========================================================================
// File:       tledMaterialAB.h
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

#ifndef tledMaterialAB_H
#define tledMaterialAB_H

#include "tledMaterial.h"

/**
 * \brief Arruda-Boyce constitutive model
 * \ingroup constitutive
 */
class tledMaterialAB : public tledMaterial
{
public:
   tledMaterialAB();
   tledMaterialAB(float* MatParams, float Density);
   virtual ~tledMaterialAB();
   
   void ComputeSPKStress(float X[3][3], float* SPK);
   void ComputeSPKStress(float X[3][3], float SPK[3][3]);
   void GetHGLameParams(float* Lambda, float* Mu);
   void GetKappa(float* kappa) {*kappa = K;}
   void SetTimeStep(float dt) {;}
   void SetParams(std::vector<float> params);

private:
   float K;	// Bulk modulus
   float M;
   float Lm;
};

#endif // tledMaterialAB_H
