// =========================================================================
// File:       tledMaterialLE.h
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


#ifndef tledMaterialLE_H
#define tledMaterialLE_H

#include "tledMaterial.h"

/**
 * \brief Linear elastic constitutive model
 * \ingroup constitutive
 */
class tledMaterialLE : public tledMaterial
{
public:
   tledMaterialLE();
   tledMaterialLE(float* MatParams, float Density);
   virtual ~tledMaterialLE();

   void ComputeSPKStress(float X[3][3], float* SPK);
   void ComputeSPKStress(float X[3][3], float SPK[3][3]);
   void GetHGLameParams(float* Lambda, float* Mu);
   void GetKappa(float* kappa){*kappa = E/(3*(1-2*P));}
   void SetTimeStep(float dt) {;}
   void SetParams(std::vector<float> params);

private:
   float E;   // Young's modulus
   float P;   // Poisson's ratio
   float c1;  // Elastic coefficients
   float c2;
   float c3;
   float c4;
};

#endif // tledMaterialLE_H
