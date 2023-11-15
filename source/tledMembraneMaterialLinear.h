// =========================================================================
// File:       tledMembraneMaterialLinear.h
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
#ifndef tledMembraneMaterialLinear_H
#define tledMembraneMaterialLinear_H
#include <vector>
#include <limits>

#include "tledShellMaterial.h"

/**
 * \brief Essentially plane-stress linear elasticity (membrane stiffness only)
 * \ingroup shell
 */
class tledMembraneMaterialLinear : public tledShellMaterial {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledShellMaterial Superclass;
  /** @} */

  /**
   * \name Traits
   * @{
   */
public:
  static const int NumberOfStrainComponents = 3;
  static const int NumberOfStressComponents = 3;
  static const int NumberOfParameters = 2;
  /** @} */

  /**
   * \name Stress Computation
   * @{
   */
public:
  virtual void ComputeStress(float *p_stressOut, const float strainIn[]) const;
  /** @} */

  /**
   * \name Parameters
   * @{
   */
private:
  float m_E, m_Nu;

public:
  /** Expects a two component vector containing Young's modulus and Poisson ratio, in that order. */
  virtual void SetParameters(const std::vector<float> &elasticParameters);
  virtual bool IsNonLinear(void) const { return false; }
  virtual bool HasBendingStiffness(void) const { return false; }

  /** Poisson ratio */
  float GetNu(void) const { return m_Nu; }  
  /** Young's modulus */
  float GetE(void) const { return m_E; }
  /** @} */

#ifdef _GPU_
  /**
   * \name GPU API
   * @{
   */
public:
  struct GPUMaterial : public tledShellMaterial::GPUMaterial {
    float E, Nu;
  };

protected:
  virtual void InitHostGPURepresentation(tledShellMaterial::GPUMaterial &r_dst) const;

public:
  virtual tledShellMaterial::GPUMaterial* InitGPU() const;
  /** @} */  
#endif 

  /**
   * \name Constructors, Destructor
   * @{
   */
public:
  tledMembraneMaterialLinear(void) {
#ifndef NDEBUG
    m_E = m_Nu = std::numeric_limits<float>::quiet_NaN();
#endif
  }

  virtual ~tledMembraneMaterialLinear(void) {}
  /** @} */
};

#endif
