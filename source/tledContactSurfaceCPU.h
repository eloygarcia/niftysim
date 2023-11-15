// =========================================================================
// File:       tledContactSurfaceCPU.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    July 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledContactSurfaceCPU_H
#define tledContactSurfaceCPU_H

#include "tledContactSurface.h"

typedef tledContactSurface tledContactSurfaceCPU;

/**
 * \brief CPU contact surface base class implementation
 * \ingroup contact 	 
 */
template <class TBaseSurface>
class tledContactSurfaceImplCPU : public TBaseSurface {
  /**
   * \name Types
   * @{
   */
public:
  typedef TBaseSurface Superclass;
  typedef typename Superclass::Facet Facet;
  /** @} */

  /**
   * \name Interpolation
   * @{
   */
public:
  static float* ComputeShapeValues(float *p_buffer, const float xi, const float eta);
  /** @} */

  /**
   * \name Surface Projection
   * @{
   */
public:
  /**
   * \brief Projects coordinate x onto a facet
   *
   * @param projOp projection operator computed with ComputeFacetProjectionOperator
   * \return true iff coordinate is inside facet, otherwise the local coordinates returned in p_xi are invalid.
   */
  static bool ProjectOntoFacet(float *p_xi, const float x[], const float projOp[]);
  /** @} */
};

#include "tledContactSurfaceCPU.tpp"
#endif
