// =========================================================================
// File:       tledMembraneMaterialNeoHookean.h
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
#ifndef tledMembraneMaterialNeoHookean_H
#define tledMembraneMaterialNeoHookean_H
#include "tledShellMaterialNonLinear.h"

/**
 * \brief Incompressible neo-hookean material model based on Bonet et al.: &quot;Finite element analysis of air supported membrane structures&quot; (2000)
 * \ingroup shell
 */
class tledMembraneMaterialNeoHookean : public tledShellMaterialNonLinear {
  /**
   * \name Traits
   * @{
   */
public:
  static const int NumberOfStrainComponents = 8;
  static const int NumberOfStressComponents = 4;
  static const int NumberOfParameters = 1;
  /** @} */

  /**
   * \name Computation
   * @{
   */
public:
  /** 
   * Computes Piola-Kirchhoff stress.<br>
   * Input deformation data format (all in row-major format):
   * <ul>
   * <li>\f$C_0\f$ (initial Cauchy-Green def. tensor): 4 components</li>
   * <li>\f$C_n\f$ (current Cauchy-Green def. tensor): 4 components</li>
   * </ul>
   */
  virtual void ComputeStress(float *p_stressOut, const float defIn[]) const;

  void ComputeStress(float *p_stressOut, const float C0[], const float Cn[]) const;
  /** @} */  

  /**
   * \name Parameters
   * @{
   */
private:
  float m_Mu;

public:
  virtual void SetParameters(const std::vector<float> &elasticParameters);
  virtual bool HasBendingStiffness(void) const { return false; }  

  /** Shear modulus */
  float GetMu(void) const { return m_Mu; }
  /** @} */

#ifdef _GPU_
  /**
   * \name GPU API
   * @{
   */
public:
  struct GPUMaterial : public tledShellMaterial::GPUMaterial {
    float Mu;
  };

public:
  virtual tledShellMaterial::GPUMaterial* InitGPU() const;
  /** @} */  
#endif 
};
#endif
