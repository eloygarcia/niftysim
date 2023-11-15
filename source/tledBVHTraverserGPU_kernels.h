// =========================================================================
// File:       tledBVHTraverserGPU_kernels.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    December 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledBVHTraverserGPU_kernels_H
#define tledBVHTraverserGPU_kernels_H

#include "tledCUDAHelpers.h"
#include "tledCUDA_operators.h"

namespace tledBVHTraverserGPU_kernels {
  template <class TProjection>
  struct NodeFacetProjectionSameSlavePredicate {
    __device__ bool operator()(const TProjection &p0, const TProjection &p1) const {
      return p0.ContactNodeIndices[0] == p1.ContactNodeIndices[0];
    }
  };

  template <class TProjection>
  struct NodeFacetProjectionOrdering {
    __device__ bool operator()(const TProjection &p0, const TProjection &p1) const {
      return p0.ContactNodeIndices[0] < p1.ContactNodeIndices[0] || (p0.ContactNodeIndices[0] == p1.ContactNodeIndices[0] && fabsf(p0.GapValue) < fabsf(p1.GapValue));
    }
  };

  template <class TProjection>
  struct EdgeEdgeProjectionSameSlavePredicate {
    __device__ __host__ bool operator()(const TProjection &p0, const TProjection &p1) const {
      return p0.SlaveEdge == p1.SlaveEdge;
    }
  };

  template <class TProjection>
  struct EdgeEdgeProjectionOrdering {
    __device__ __host__ bool operator()(const TProjection &p0, const TProjection &p1) const {
      const int2 &slave0 = p0.SlaveEdge;
      const int2 &slave1 = p1.SlaveEdge;

      return slave0.x < slave1.x || (slave0.x == slave1.x && slave0.y < slave1.y) || (slave0 == slave1 && p0.Xi.x < p1.Xi.x);
    }
  };

  /** Ordering suitable for removing dupes, only compares indices */
  template <class TProjection, const int t_numMasterIndices>
  struct NodeFacetProjectionIndexOrdering {
    __device__ bool operator()(const TProjection &a, const TProjection &b) const {
      for (int i = 0; i < t_numMasterIndices + 1; i++) {
	if (a.ContactNodeIndices[i] < b.ContactNodeIndices[i]) return true;
	if (a.ContactNodeIndices[i] > b.ContactNodeIndices[i]) return false;
      }

      return false;
    }
  };

  /** Ordering suitable for removing dupes, only compares indices */
  template <class TProjection>
  class EdgeEdgeProjectionIndexOrdering {
  private:
    __device__ bool _Lower(const int2 &a, const int2 &b) const {
      return a.x < b.x || (a.x == b.x && a.y < b.y);
    }

  public:
    __device__ bool operator()(const TProjection &a, const TProjection &b) const {    
      return _Lower(a.SlaveEdge, b.SlaveEdge) || (a.SlaveEdge == b.SlaveEdge && _Lower(a.MasterEdge, b.MasterEdge));
    }
  };

  struct BroadPhaseResultOrdering {
    __device__ bool operator()(const int2 &a, const int2 &b) const {
      return a.x < b.x || (a.x == b.x && a.y < b.y);
    }
  };

  template <class TProjection, class TSurface>
  __device__ bool ProjectOntoSurfaceC0(TProjection &r_proj, const float3 &x, const TSurface &mesh, const float3 nodeNormals[]);

  template <class TBVHTraverser, class TMasterBV, class TSlaveBV, class TNarrowPhaseCommitFunction, const int t_blockSize>
  __global__ void DoBroadPhaseSelfCollisionSearch(TBVHTraverser *p_traverser, const TMasterBV *masterBVs, const TSlaveBV *slaveBVs, const TNarrowPhaseCommitFunction commitFunct);

  /** 1-2-1 CUDA port of tledBVHTraverserCPU::ComputeEdgeEdgeClosestPointParameters */
  __device__ bool ComputeEdgeEdgeClosestPointParameters(float &r_r, float &r_q, const float3 &A, const float3 &B, const float3 &C, const float3 &D);

  /** 1-2-1 CUDA port of tledBVHTraverserCPU::ComputeEdgeEdgePenetrationDepth, including checks. */
  template <const bool t_isMovingSlave, const bool t_isMovingMaster>
  __device__ bool ComputeEdgeEdgePenetrationDepth(float &r_g, float3 &r_n, 
						  const float3 &slave0T0, const float3 &slave1T0, const float3& slave0T1, const float3 &slave1T1, const float3 &slaveN0, const float3 &slaveN1, const float r, 
						  const float3 &master0T0, const float3 &master1T0, const float3 &master0T1, const float3 &master1T1, const float3 &masterN0, const float3 &masterN1, const float q,
						  const float safetyMargin);

  /** Performs the first stage of edge-edge collision detection: determining the local coordinates of the contact points. */
  template <const int t_blockSize, class TMasterSurface, class TSlaveSurface, class TEdgeEdgeResult>
  __global__ void EdgeEdgeNarrowPhaseStage1(TEdgeEdgeResult *p_edgeList, int *p_numEdges, const TMasterSurface *pc_mMesh, const TSlaveSurface *pc_sMesh, const int2 *narrowPhasePairs, const int numPairs);

  /** Performs the iterative projection of a node onto a surface based on intiatial guesses obtained in the first stage of the node-facet narrow-phase. */
  template <const int t_blockSize, class TFinalProjection, class TInitialProjection, class TMasterMesh, class TSlaveMesh>
  __global__ void NodeFacetNarrowPhaseStage2(TFinalProjection *p_projections, int *p_numProjections, const TInitialProjection *initialProjections, const int numInitial, const TMasterMesh *pc_mMesh, const TSlaveMesh *pc_sMesh);
}

#endif
