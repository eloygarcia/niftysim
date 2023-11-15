// =========================================================================
// File:       tledShellMaterialLinearThickPlateDecorator.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    July 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledShellMaterialLinearThickPlateDecorator_H
#define tledShellMaterialLinearThickPlateDecorator_H

#include "tledShellMaterialLinearPlateDecorator.h"

/**
 * \name Thick shell material with shear stresses
 * \ingroup shell
 */
template <class TMembraneMaterial>
class tledShellMaterialLinearThickPlateDecorator : public tledShellMaterialLinearPlateDecorator<TMembraneMaterial> {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledShellMaterialLinearPlateDecorator<TMembraneMaterial> Superclass;
  /** @} */

  /**
   * \name Traits
   * @{
   */
public:
  static const int NumberOfStrainComponents = Superclass::NumberOfStrainComponents;
  static const int NumberOfStressComponents = Superclass::NumberOfStressComponents;
  static const int NumberOfParameters = Superclass::NumberOfParameters;
  /** @} */

  /** 
   * \name Parameters
   * @{
   */
public:
  float GetPlateShearG(void) const { return this->GetBendingE()/(2 + 2*this->GetBendingNu()); }
  /** @} */

  /**
   * \name Computation
   * @{
   */
public:
  virtual void ComputeShearStress(float *p_stressOut, const float shearStrain[]) const;
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  virtual ~tledShellMaterialLinearThickPlateDecorator(void) {}
  /** @} */
};

#include "tledShellMaterialLinearThickPlateDecorator.tpp"
#endif
