// =========================================================================
// File:       tledElementMembraneImpl.h
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

#ifndef tledElementMembraneImpl_H
#define tledElementMembraneImpl_H

#include <cmath>
#include <cassert>

#include "tledElementMembrane.h"

/** Internal use only */
template <const int t_numFacetVertices>
class tledElementMembraneImpl : public tledElementMembrane {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledElementMembrane Superclass;
  typedef tledShellMesh<t_numFacetVertices> Surface;
  typedef typename Surface::Facet Facet;
  /** @} */

  /**
   * \name CPU Force Computation
   * @{
   */
protected:
  virtual float* ComputeStrain(float *p_dst, const float U[]) = 0;
  virtual void ComputeForces(float *p_dst, const float stress[]) = 0;

  /** Basic element force computation member function can be used by child classes for overloading tledElementMembrane::ComputeElementForces */
  template <const int t_numStrainComponents, const int t_numStressComponents>
  void ComputeElementForcesTemplate(float *p_dst, const float U[]);
  /** @} */

  /**
   * \name Element Geometry
   * @{
   */
private:
  int m_FacetVertexIndices[t_numFacetVertices];

protected:
  /** Computes the element basis from the node positions held in X, as a byproduct the element area is computed, too. */
  void ComputeElementBasis(float (*p_basis)[3], float &r_area, const float X[]) const;

public:
  const int* GetFacetVertexIndices(void) const { return m_FacetVertexIndices; }
  virtual void InitialiseElement(const Surface &surf, const int facetIndex);
  virtual void ComputeElementMass(float *p_dst) const;
  /** @} */

#ifdef _GPU_
  /**
   * \name GPU On-Device Representation
   * @{
   */
public:
  virtual void InitGPU(tledElementMembrane::GPUElement &r_dst);
  /** @} */
#endif
};

#include "tledElementMembraneImpl.tpp"

#endif
