// =========================================================================
// File:       tledElementShellBSTP1_kernels.cu
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   Cuda
// Created:    August 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledElementShellBSTP1_kernels_CU
#define tledElementShellBSTP1_kernels_CU

#include "tledElementShellBSTP1.h"
#include "tledElementShellBSTP1Edge.h"
#include "tledShellSolver_kernels.h"
#include "tledSolverGPU_kernels.h"

namespace tledElementShellBSTP1_kernels {
  __device__ bool _IsComplete(const tledElementShellBSTP1Edge::GPUElement &edgeElement) {
    return edgeElement.NeighbourNodeIndex >= 0;
  }

  __device__ bool _IsClamped(const tledElementShellBSTP1Edge::GPUElement &edgeElement) {
    return !_IsComplete(edgeElement) && edgeElement.EdgeNormalAndLength0.z == edgeElement.EdgeNormalAndLength0.z;
  }

  __device__ int _GetBufferIndexForEdgePatchNode(const int n, const int edgeIndex) {
    tledCudaAssert(n < 4 && n >= 0);

    if (n == 3) return 3 + edgeIndex;
    else return (edgeIndex + n + 1)%3;
  }

  __device__ void _ComputeDeformationGradient(tledElementShellBSTP1Edge::GPUElement &r_edgeElement, const float3 triaX[], const int elementEdgeIndex) {    
    tledCudaAssert(r_edgeElement.NeighbourNodeIndex >= 0);

    {
      const float3 xNeighbour = r_edgeElement.NeighbourX0 + tledSolverGPU_kernels::GetCurrentDisplacement(r_edgeElement.NeighbourNodeIndex);    
      
      r_edgeElement.DeformationGradient[0] = xNeighbour*r_edgeElement.NeighbourPatchShapeGradient.x;
      r_edgeElement.DeformationGradient[1] = xNeighbour*r_edgeElement.NeighbourPatchShapeGradient.y;
    }

    for (int n = 0; n < 3; n++) {
      r_edgeElement.DeformationGradient[0] += triaX[(elementEdgeIndex+n+1)%3]*r_edgeElement.PatchShapeGradients[n].x;
      r_edgeElement.DeformationGradient[1] += triaX[(elementEdgeIndex+n+1)%3]*r_edgeElement.PatchShapeGradients[n].y;    
    }
  }

  __device__ void _ComputeH(float3 *p_hs, tledElementShellBSTP1Edge::GPUElement &r_edgeElement, const tledElementShellBSTP1::GPUElement &tria, const int elementEdgeIndex) {
    const float2 linearGradient = tria.LinearShapeGradients[elementEdgeIndex];
   
    p_hs[0] += 2*linearGradient.x*r_edgeElement.DeformationGradient[0];
    p_hs[1] += 2*linearGradient.y*r_edgeElement.DeformationGradient[1];
    p_hs[2] += 2*(linearGradient.y*r_edgeElement.DeformationGradient[0] + linearGradient.x*r_edgeElement.DeformationGradient[1]);
    tledCudaAssert(_IsComplete(r_edgeElement) || (std::fabs(dot(linearGradient.x*r_edgeElement.DeformationGradient[0], tria.CurrentNormal)) < 1e-3f && std::fabs(dot(linearGradient.y*r_edgeElement.DeformationGradient[1], tria.CurrentNormal)) < 1e-3f));
  }

  __device__ void _ComputeClampedOrSymmetryH(float3 *p_hs, tledElementShellBSTP1Edge::GPUElement &r_edgeElement, const tledElementShellBSTP1::GPUElement &tria, const float3 triaX[], const float thicknessRatio, const int elementEdgeIndex) {
    const float2 linearGradient = tria.LinearShapeGradients[elementEdgeIndex];

    float3 phiS, phiN;
    
    phiS = (triaX[(elementEdgeIndex+1+1)%3] - triaX[(elementEdgeIndex+1)%3])*2/r_edgeElement.EdgeNormalAndLength0.z;
    phiN = r_edgeElement.PhiN0/(norm(phiS)*thicknessRatio/4);

    p_hs[0] -= linearGradient.x*r_edgeElement.EdgeNormalAndLength0.y*phiS;
    p_hs[0] += linearGradient.x*r_edgeElement.EdgeNormalAndLength0.x*phiN;

    p_hs[1] += linearGradient.y*r_edgeElement.EdgeNormalAndLength0.x*phiS;
    p_hs[1] += linearGradient.y*r_edgeElement.EdgeNormalAndLength0.y*phiN;

    p_hs[2] += (linearGradient.x*r_edgeElement.EdgeNormalAndLength0.x - linearGradient.y*r_edgeElement.EdgeNormalAndLength0.y)*phiS;
    p_hs[2] += (linearGradient.y*r_edgeElement.EdgeNormalAndLength0.x + linearGradient.x*r_edgeElement.EdgeNormalAndLength0.y)*phiN;
  }

  __device__ void _ComputeMembraneForces(float4 *p_fDst, const float sigma[], const tledElementShellBSTP1Edge::GPUElement &edgeElement, const tledElementShellBSTP1::GPUElement &tria, const int elementEdgeIndex) {
    float3 tmpF[2], tmp;
    float3 sigmaEdge;    

    sigmaEdge = make_float3(sigma[0]/3, sigma[1]/3, sigma[2]/3);

    tmpF[0] = edgeElement.DeformationGradient[0]*sigmaEdge.x + edgeElement.DeformationGradient[1]*sigmaEdge.z;   
    tmpF[1] = edgeElement.DeformationGradient[1]*sigmaEdge.y + edgeElement.DeformationGradient[0]*sigmaEdge.z;

    for (int n = 0; n < 3; n++) {
      tmp = (tmpF[0]*edgeElement.PatchShapeGradients[n].x + tmpF[1]*edgeElement.PatchShapeGradients[n].y)*tria.Area*tria.Thickness;
      p_fDst[_GetBufferIndexForEdgePatchNode(n, elementEdgeIndex)] += make_float4(tmp.x, tmp.y, tmp.z, 0);
      tledCudaAssert(_GetBufferIndexForEdgePatchNode(n, elementEdgeIndex) < 3);
      tledCudaAssert(p_fDst[_GetBufferIndexForEdgePatchNode(n, elementEdgeIndex)].w == 0);
    }

    if (_IsComplete(edgeElement)) {
      tmp = (tmpF[0]*edgeElement.NeighbourPatchShapeGradient.x + tmpF[1]*edgeElement.NeighbourPatchShapeGradient.y)*tria.Area*tria.Thickness;
      p_fDst[_GetBufferIndexForEdgePatchNode(3, elementEdgeIndex)] = make_float4(tmp.x, tmp.y, tmp.z, 0);
      tledCudaAssert(_GetBufferIndexForEdgePatchNode(3, elementEdgeIndex) >= 3 && _GetBufferIndexForEdgePatchNode(3, elementEdgeIndex) < 6);
      tledCudaAssert(p_fDst[_GetBufferIndexForEdgePatchNode(3, elementEdgeIndex)].w == 0);
    }
  }

  __device__ void _ComputeBendingForces(float4 *p_fDst, const float M[], const tledElementShellBSTP1Edge::GPUElement &edgeElement, const float3 hs[], const tledElementShellBSTP1::GPUElement &tria, const int elementEdgeIndex) {
    const float2 linearGradient = tria.LinearShapeGradients[elementEdgeIndex];

    {
      float3 gradientRho;
      float2 rho;

      rho.x = dot(tria.ContravariantBasis[0], hs[0]);
      rho.y = dot(tria.ContravariantBasis[1], hs[0]);

      gradientRho.x = dot(linearGradient, rho);
    
      rho.x = dot(tria.ContravariantBasis[0], hs[1]);
      rho.y = dot(tria.ContravariantBasis[1], hs[1]);
    
      gradientRho.y = dot(linearGradient, rho);
    
      rho.x = dot(tria.ContravariantBasis[0], hs[2]);
      rho.y = dot(tria.ContravariantBasis[1], hs[2]);
    
      gradientRho.z = dot(linearGradient, rho);      
    
      p_fDst[_GetBufferIndexForEdgePatchNode(2, elementEdgeIndex)] -= 2*tria.Area*dot(gradientRho, make_float3(M[0], M[1], M[2]))*tria.CurrentNormal;
    }

    if (_IsClamped(edgeElement)) {
      float tmpHM;

      tmpHM = -M[0]*linearGradient.x*edgeElement.EdgeNormalAndLength0.y + M[1]*linearGradient.y*edgeElement.EdgeNormalAndLength0.x + M[2]*(linearGradient.x*edgeElement.EdgeNormalAndLength0.x + linearGradient.y*edgeElement.EdgeNormalAndLength0.y);
      tmpHM *= 2*tria.Area/edgeElement.EdgeNormalAndLength0.z;
      p_fDst[_GetBufferIndexForEdgePatchNode(0, elementEdgeIndex)] -= tria.CurrentNormal*tmpHM;
      p_fDst[_GetBufferIndexForEdgePatchNode(1, elementEdgeIndex)] += tria.CurrentNormal*tmpHM;
    } else {
      float2 tmp;

      tmp.x = linearGradient.x*M[0] + linearGradient.y*M[2];
      tmp.y = linearGradient.y*M[1] + linearGradient.x*M[2];
      tmp *= 2*tria.Area;

      for (int n = 0; n < 3; n++) p_fDst[_GetBufferIndexForEdgePatchNode(n, elementEdgeIndex)] += tria.CurrentNormal*dot(edgeElement.PatchShapeGradients[n], tmp);
      if (_IsComplete(edgeElement)) p_fDst[_GetBufferIndexForEdgePatchNode(3, elementEdgeIndex)] += tria.CurrentNormal*dot(edgeElement.NeighbourPatchShapeGradient, tmp);
    }    
  }
}

namespace tledElementMembrane_kernels {  
  template <>
  __device__ void ComputeStrain<tledElementShellBSTP1>(float *p_strainDst, tledElementShellBSTP1::GPUElement *gp_elements) {
    using namespace tledElementShellBSTP1_kernels;

    const int elInd = blockIdx.x*blockDim.x + threadIdx.x;

    tledElementShellBSTP1::GPUElement &r_element = gp_elements[elInd];
    float3 centralDefGrad[2], xCurr[3];
    float thicknessRatio;

    for (int c = 0; c < 2; c++) centralDefGrad[c] = make_float3(0, 0, 0);

    for (int n = 0; n < 3; n++) {
      xCurr[n] = r_element.X0s[n] + tledSolverGPU_kernels::GetCurrentDisplacement(n == 0? r_element.ElementNodeIndices.x : (n == 1? r_element.ElementNodeIndices.y : r_element.ElementNodeIndices.z));
      centralDefGrad[0] += xCurr[n]*r_element.LinearShapeGradients[n].x;
      centralDefGrad[1] += xCurr[n]*r_element.LinearShapeGradients[n].y;
    }

    r_element.CurrentNormal = cross(centralDefGrad[0], centralDefGrad[1]);
    thicknessRatio = 1/norm(r_element.CurrentNormal);
    r_element.CurrentNormal *= thicknessRatio;

    r_element.ContravariantBasis[0] = cross(centralDefGrad[1], r_element.CurrentNormal)*thicknessRatio;
    r_element.ContravariantBasis[1] = cross(r_element.CurrentNormal, centralDefGrad[0])*thicknessRatio;    

    for (int c = 0; c < 6; c++) p_strainDst[c] = 0;
    for (int c = 0; c < 3; c++) r_element.H[c] = make_float3(0, 0, 0);

    for (int e = 0; e < 3; e++) {
      if (_IsComplete(r_element.EdgeElements[e])) {
	_ComputeDeformationGradient(r_element.EdgeElements[e], xCurr, (e + 2)%3);
      } else {
	r_element.EdgeElements[e].DeformationGradient[0] = centralDefGrad[0], r_element.EdgeElements[e].DeformationGradient[1] = centralDefGrad[1];
	tledCudaAssert(std::fabs(dot(r_element.EdgeElements[e].DeformationGradient[0], r_element.CurrentNormal)) < 1e-3f && std::fabs(dot(r_element.EdgeElements[e].DeformationGradient[1], r_element.CurrentNormal)) < 1e-3f);	
      }

      if (_IsClamped(r_element.EdgeElements[e])) _ComputeClampedOrSymmetryH(r_element.H, r_element.EdgeElements[e], r_element, xCurr, thicknessRatio, (e + 2)%3);
      else _ComputeH(r_element.H, r_element.EdgeElements[e], r_element, (e + 2)%3);
      
      for (int c = 0; c < 2; c++) p_strainDst[c] += dot(r_element.EdgeElements[e].DeformationGradient[c], r_element.EdgeElements[e].DeformationGradient[c]);
      p_strainDst[2] += dot(r_element.EdgeElements[e].DeformationGradient[0], r_element.EdgeElements[e].DeformationGradient[1]);      
    }

    for (int c = 0; c < 2; c++) p_strainDst[c] = 0.5*(p_strainDst[c]/3 - 1);
    p_strainDst[2] /= 3;
    for (int c = 0; c < 3; c++) p_strainDst[3+c] = dot(r_element.CurrentNormal, r_element.H[c]);
  }

  template <>
  __device__ void ComputeForces<tledElementShellBSTP1>(float4 *p_fDst, tledElementShellBSTP1::GPUElement *gp_elements, const float stress[]) {
    using namespace tledElementShellBSTP1_kernels;

    const int elInd = blockIdx.x*blockDim.x + threadIdx.x;
    const tledElementShellBSTP1::GPUElement &element = gp_elements[elInd];

    for (int n = 0; n < 3; n++) *(p_fDst + elInd*3*2 + n) = make_float4(0, 0, 0, 0);

    for (int e = 0; e < 3; e++) {
      _ComputeMembraneForces(p_fDst + elInd*3*2, stress, element.EdgeElements[e], element, (e + 2)%3);
      _ComputeBendingForces(p_fDst + elInd*3*2, stress + 3, element.EdgeElements[e], element.H, element, (e + 2)%3);
    }
  }
}
#endif
