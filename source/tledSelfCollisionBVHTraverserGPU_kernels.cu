// =========================================================================
// File:       tledSelfCollisionBVHTraverserGPU_kernels.cu
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   CUDA
// Created:    July 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledSelfCollisionBVHTraverserGPU_kernels_CU
#define tledSelfCollisionBVHTraverserGPU_kernels_CU

namespace tledSelfCollisionBVHTraverserGPU_kernels {
  template <class TSurface>
  struct _PrimitiveNonAdjacencyPredicate {
  private:
    const TSurface &mc_Surface;

  public:
    __device__ bool operator()(const int2 &p) const {
      return !tledDeformableContactSurfaceGPU_kernels::IsAdjacent<TSurface>(mc_Surface, p.x, p.y);
    }

  public:
    __device__ _PrimitiveNonAdjacencyPredicate(const TSurface &mesh) : mc_Surface(mesh) {}
  };

  template <const int t_blockSize, class TSurface>
  __device__ void RemoveAdjacent(int2 *p_list, unsigned short &r_numItems, const typename TSurface::GPUSurface &mesh) {
    tledCUDAHelpers::KeepIf<t_blockSize, int2, _PrimitiveNonAdjacencyPredicate<TSurface> > (p_list, r_numItems, _PrimitiveNonAdjacencyPredicate<TSurface>(mesh));
  }
}
#endif
