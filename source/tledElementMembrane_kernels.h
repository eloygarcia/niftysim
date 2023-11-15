// =========================================================================
// File:       tledElementMembrane_kernels.h
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
#ifndef tledElementMembrane_kernels_H
#define tledElementMembrane_kernels_H

#ifdef _GPU_
/**
 * \brief Namespace of GPU shell/membrane element functions
 * \ingroup shell
 */
namespace tledElementMembrane_kernels {
  /** Computes the strain from the current displacement */
  template <class TElement>
  __device__ void ComputeStrain(float *p_strain, typename TElement::GPUElement *gp_elements);

  /** Computes the vertex nodal forces from the given stress */
  template <class TElement>
  __device__ void ComputeForces(float4 *p_fDst, typename TElement::GPUElement *gp_elements, const float stress[]);
}

#endif
#endif
