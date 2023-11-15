// =========================================================================
// File:       tledShellMaterialLinearPlateDecorator.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    May 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledShellMaterialLinearPlateDecorator_H
#define tledShellMaterialLinearPlateDecorator_H

#include "tledShellMaterial.h"
#ifdef _GPU_
#include "tledCUDAHelpers.h"
#endif

#include <cassert>


/**
 * \name Adds bending stiffness to properties to membrane materials.
 * \ingroup shell
 */
template <class TMembraneMaterial>
class tledShellMaterialLinearPlateDecorator : public TMembraneMaterial {
  /**
   * \name Traits
   * @{
   */
public:
  static const int NumberOfStrainComponents = TMembraneMaterial::NumberOfStrainComponents + 3;
  static const int NumberOfStressComponents = TMembraneMaterial::NumberOfStressComponents + 3;
  static const int NumberOfParameters = TMembraneMaterial::NumberOfParameters + 2;
  /** @} */

  /**
   * \name Types
   * @{
   */
public:
  typedef TMembraneMaterial Superclass;
  /** @} */

  /**
   * \name Parameters
   * @{
   */
private:
  float m_E, m_Nu;

public:
  /** Expects a two component vector containing Young's modulus and Poisson ratio, in that order, following the membrane material parameters. */
  virtual void SetParameters(const std::vector<float> &elasticParameters);

  /** Poisson ratio */
  float GetBendingNu(void) const { return m_Nu; }  
  /** Young's modulus */
  float GetBendingE(void) const { return m_E; }

  virtual bool HasBendingStiffness(void) const { return true; }
  /** @} */

  /**
   * \name Computation
   * @{
   */
public:
  virtual void ComputeBendingResponse(float *p_stress, const float curvature[]) const;

  /** Computes the stress for a given strain. */
  virtual void ComputeStress(float *p_stressOut, const float strainCurvatures[]) const;
  /** @} */

#ifdef _GPU_
  /**
   * \name GPU API
   * @{
   */
public:
  struct GPUMaterial : public Superclass::GPUMaterial {
    float BendingE, BendingNu;
  };

protected:
  void InitHostGPURepresentation(tledShellMaterial::GPUMaterial &r_sdst) const;

public:
  virtual tledShellMaterial::GPUMaterial* InitGPU() const;
  /** @} */  
#endif 

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  virtual ~tledShellMaterialLinearPlateDecorator(void) {}
  /** @} */
};

#include "tledShellMaterialLinearPlateDecorator.tpp"
#endif
