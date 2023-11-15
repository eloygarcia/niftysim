// =========================================================================
// File:       tledShellMaterialNonLinear.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    March 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledShellMaterialNonLinear_H
#define tledShellMaterialNonLinear_H
#include "tledShellMaterial.h"

/**
 * \brief Base class for non-linear shell/membrane consistitutive models.
 * \ingroup shell
 */
class tledShellMaterialNonLinear : public tledShellMaterial {
  /**
   * \name Traits
   * @{
   */
public:
  static const int NumberOfStrainComponents = 8;
  static const int NumberOfStressComponents = 4;
  /** @} */

  /**
   * \name Computation
   * @{
   */
public:
  virtual void ComputeStress(float *p_stressOut, const float strainIn[]) const = 0;
  /** @} */  

  /**
   * \name Parameters
   * @{
   */
public:
  virtual bool IsNonLinear(void) const { return true; }
  /** @} */
};

#endif
