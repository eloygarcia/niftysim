// =========================================================================
// File:       tledElementMembraneSimpleLinear.h
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
#ifndef tledElementMembraneSimpleLinear_H
#define tledElementMembraneSimpleLinear_H

#include "tledElementMembraneImpl.h"

#ifdef _GPU_
#include <vector_functions.h>
#endif

#include <vector>
#include <limits>
#include <algorithm>

template <const int t_numFacetVertices>
class tledElementMembraneSimpleLinear : public tledElementMembraneImpl<t_numFacetVertices> {
  /**
   * \name Traits
   * @{
   */
public:
  static const int NumberOfStrainComponents = 3;
  static const int NumberOfStressComponents = 3;
  /** @} */

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
protected:
  //float m_B[3][2*t_numFacetVertices];
  float m_PhiInvT[2][2];
  float m_ElementBasis[3][3];
  float m_X0[t_numFacetVertices][3];

#ifndef NDEBUG
  int m_FacetIndex;
#endif

public:
  /** Updates the corotational basis */
  void UpdateBasis(const float x0[], const float u[]);

  virtual void InitialiseElement(const typename Superclass::Surface &surf, const int facetIndex);
  /** @} */

#ifdef _GPU_
  /**
   * \name GPU On-Device Representation
   * @{
   */
public:
  struct GPUElement : public Superclass::GPUElement {
    float3 ElementBasis[3];
    float2 PhiInvT[2];
  };

public:
  virtual void InitGPU(tledElementMembrane::GPUElement &r_dst);
  /** @} */
#endif

  /**
   * \name Forces, Materials, etc.
   * @{
   */
protected:
  virtual void ComputeForces(float *p_dst, const float planeStress[]);
  virtual float* ComputeStrain(float *p_dst, const float U[]);  

public:
  virtual void ComputeElementForces(float *p_dst, const float U[]) { this->template ComputeElementForcesTemplate<3, 3>(p_dst, U); }
  virtual float ComputeCurrentThickness(const float U[]) const;

  /** Returns initial normal, no updating as only admissible for very small displacements! */
  virtual float* ComputeCurrentNormal(float *p_n, const float U[]) const { std::copy(m_ElementBasis[2], m_ElementBasis[2] + 3, p_n); return p_n; }
  /** @} */

public:
  virtual ~tledElementMembraneSimpleLinear() {}
};

#include "tledElementMembraneSimpleLinear.tpp"
#endif
