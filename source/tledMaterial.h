// =========================================================================
// File:       tledMaterial.h
// Purpose:    Material base class
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


#ifndef tledMaterial_H
#define tledMaterial_H

#include <vector>

/**
 * \defgroup constitutive Constitutive Models
 */

/**
 * \brief Constitutive model base class
 * \ingroup constitutive
 */
class tledMaterial
{
  /**
   * \name Construction, Destruction
   * @{
   */
public:
   tledMaterial();
   tledMaterial(const float rho);

   virtual ~tledMaterial() {;}
   /** @} */

   /**
    * \name Mass Density
    * @{
    */
private:
   float m_D;	// Density

public:
   float GetDensity(void) const { return m_D; }
   void SetDensity(const float rho) { m_D = rho; }
   /** @} */

   virtual void ComputeSPKStress(float X[3][3], float* SPK) = 0;
   virtual void ComputeSPKStress(float X[3][3], float SPK[3][3]) = 0;
   virtual void GetHGLameParams(float* Lambda, float* Mu) = 0;
   virtual void GetKappa(float* kappa) = 0;
   virtual void SetTimeStep(float dt) = 0;
   virtual void SetParams(std::vector<float> params) = 0;   
};

#endif // tledMaterial_H
