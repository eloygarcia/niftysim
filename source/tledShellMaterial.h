// =========================================================================
// File:       tledShellMaterial.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    February 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledShellMaterial_H
#define tledShellMaterial_H

#include <vector>
#include <limits>

/**
 * \brief Constitutive models for shell/membrane modelling
 * \ingroup shell
 */
class tledShellMaterial {
  /**
   * \name Computation
   * @{
   */
public:
  /** Computes the stress for a given strain. Expected strain representation is material dependent. */
  virtual void ComputeStress(float *p_stressOut, const float strainIn[]) const = 0;
  /** @} */

  /**
   * \name Parameters
   * @{
   */
private:
  float m_Rho, m_Thickness;

public:
  void SetDensity(const float rho) { m_Rho = rho; }
  float GetDensity(void) const { return m_Rho; }

  void SetThickness(const float t) { m_Thickness = t; }
  float GetThickness(void) const { return m_Thickness; }

  /** Setter for elastic parameters. Density, thickness have their own setters. */
  virtual void SetParameters(const std::vector<float> &elasticParameters) = 0;

  /** Query type of constitutive equation (linear/non-linear). */
  virtual bool IsNonLinear(void) const = 0;

  /** Bending stiffness query (has plate/membrane-only) */
  virtual bool HasBendingStiffness(void) const = 0;
  /** @} */

#ifdef _GPU_
  /**
   * \name GPU API
   * @{
   */
public:
  struct GPUMaterial {
    float Rho, Thickness;
  };

protected:
  virtual void InitHostGPURepresentation(tledShellMaterial::GPUMaterial &r_dst) const { r_dst.Rho = m_Rho, r_dst.Thickness = m_Thickness; }

public:
  virtual GPUMaterial* InitGPU() const = 0;
  /** @} */  
#endif 

  /**
   * \name Constructors, Destructor
   * @{
   */
public:
  tledShellMaterial(void) {
#ifndef NDEBUG
    m_Rho = std::numeric_limits<float>::quiet_NaN();
#endif
  }

  virtual ~tledShellMaterial(void) {}
  /** @} */
};
#endif
