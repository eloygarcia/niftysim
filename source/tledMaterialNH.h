// =========================================================================
// File:       tledMaterialNH.h
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


#ifndef tledMaterialNH_H
#define tledMaterialNH_H

#include "tledMaterial.h"

/**
 * \brief Neo-Hookean constitutive model
 * \ingroup constitutive
 */
class tledMaterialNH : public tledMaterial
{
public:
   tledMaterialNH();
   tledMaterialNH(float* MatParams, float Density);
   virtual ~tledMaterialNH();

   void ComputeSPKStress(float X[3][3], float* SPK);
   void ComputeSPKStress(float X[3][3], float SPK[3][3]);
   void GetHGLameParams(float* Lambda, float* Mu);
   void GetKappa(float* kappa){*kappa = K;}
   void SetTimeStep(float dt) {;}
   void SetParams(std::vector<float> params);

private:
   float K;	// Bulk modulus
   float M;
};

#endif // tledMaterialNH_H
