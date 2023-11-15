// =========================================================================
// File:       tledElementMembrane.h
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
#ifndef tledElementMembrane_H
#define tledElementMembrane_H

#ifdef _GPU_
#include <vector_types.h>
#endif

#include <cstdlib>

#include "tledShellMesh.h"
#include "tledShellMaterial.h"

/**
 * \brief Membrane element interface 
 * \ingroup shell
 */
class tledElementMembrane {
  /**
   * \name Forces, Materials, etc.
   * @{
   */
private:
  const tledShellMaterial *mpc_Material;
  
public:
  virtual void ComputeElementMass(float *p_dst) const = 0;
  virtual void ComputeElementForces(float *p_dst, const float U[]) = 0;  

  void SetMaterial(const tledShellMaterial &mat) { mpc_Material = &mat; }
  const tledShellMaterial& GetMaterial(void) const { return *mpc_Material; }

  /** Returns current thickness at element centre */
  virtual float ComputeCurrentThickness(const float U[]) const = 0;

  /** Computes current normal */
  virtual float* ComputeCurrentNormal(float *p_n, const float U[]) const = 0;
  /** @} */

#ifdef _GPU_
  /**
   * \name Static API and Types Associated with GPU On-Device Representation
   * @{
   */
public:
  struct GPUElement {
    float Area, Thickness;
    int4 ElementNodeIndices;
  };

public:
  virtual void InitGPU(GPUElement &r_dst) = 0;
  /** @} */
#endif

  /**
   * \name Element Geometry
   * @{
   */
private:
  float m_Area;

protected:
  void SetArea(const float area) { m_Area = area; }

public:
  float GetThickness(void) const { return GetMaterial().GetThickness(); }
  float GetArea(void) const { return m_Area; }
  float GetDensity(void) const { return GetMaterial().GetDensity(); }
  /** @} */

public:
  virtual ~tledElementMembrane() {}
};

#include "tledElementMembrane.tpp"
#endif
