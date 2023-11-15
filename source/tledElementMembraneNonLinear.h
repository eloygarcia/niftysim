// =========================================================================
// File:       tledElementMembraneNonLinear.h
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
#ifndef tledElementMembraneNonLinear_H
#define tledElementMembraneNonLinear_H

#include "tledElementMembraneImpl.h"

#ifdef _GPU_
#include <vector_functions.h>
#endif

/**
 * Large-deformation membrane based on Bonet, J., Wood, R., Mahaney, J., & Heywood, P. (2000). Finite element analysis of air supported membrane structures. Computer Methods in Applied Mechanics and Engineering, 190(5-7), 579–595. doi:10.1016/S0045-7825(99)00428-4
 * \ingroup shell
 */
template <const int t_numFacetVertices>
class tledElementMembraneNonLinear : public tledElementMembraneImpl<t_numFacetVertices> {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledElementMembraneImpl<t_numFacetVertices> Superclass;
  typedef tledShellMesh<t_numFacetVertices> Surface;
  typedef typename Surface::Facet Facet;
  /** @} */

  /**
   * \name Element Geometry
   * @{
   */
public:
  virtual void InitialiseElement(const Surface &surf, const int facetIndex);
  /** @} */

  /**
   * \name Forces, Materials, etc.
   * @{
   */
private:
  float m_InitialCauchyGreenTensor[2][2], m_InitialJacobian[2][3], m_J[2][3], m_LastStrain[8];

private:
  static void _ComputeCauchyGreenTensor(float *p_tensor, const float J[][3]);

protected:
  virtual float* ComputeStrain(float *p_dst, const float U[]);
  virtual void ComputeForces(float *p_dst, const float stress[]);

public:
  virtual void ComputeElementForces(float *p_dst, const float U[]) { this->template ComputeElementForcesTemplate<8, 4>(p_dst, U); }
  virtual float ComputeCurrentThickness(const float U[]) const;
  virtual float* ComputeCurrentNormal(float *p_n, const float U[]) const { return tledVectorArithmetic::Cross(p_n, m_J[0], m_J[1]); }
  /** @} */

#ifdef _GPU_
  /**
   * \name GPU On-Device Representation
   * @{
   */
public:
  struct GPUElement : public Superclass::GPUElement {
    float3 InitialJacobian[2];
    float2 InitialCauchyGreenTensor[2];
    float3 J[2];
  };

public:
  virtual void InitGPU(tledElementMembrane::GPUElement &r_dst);
  /** @} */
#endif

public:
  tledElementMembraneNonLinear(void);
  virtual ~tledElementMembraneNonLinear() {}
};

#include "tledElementMembraneNonLinear.tpp"
#endif
