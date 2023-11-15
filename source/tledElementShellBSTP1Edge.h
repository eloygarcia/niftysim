// =========================================================================
// File:       tledElementShellBSTP1Edge.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    July 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledElementShellBSTP1Edge_H
#define tledElementShellBSTP1Edge_H

#include "tledElementMembraneImpl.h"
#include "tledSurfaceTopology.h"

#ifdef _GPU_
#include <vector_functions.h>
#endif

class tledElementShellBSTP1;

class tledElementShellBSTP1Edge {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledShellMesh<3> Surface;
  friend class tledElementShellBSTP1;
  /** @} */

  /**
   * \name Computation
   * @{
   */
private:
  float m_DeformationGradient[2][3];
  float m_PatchShapeGradients[4][2];
  static const float s_PatchShapeGradients[3][4][2];

public:
  const float* GetDeformationGradient(const int component) const { return m_DeformationGradient[component]; }
  float* GetDeformationGradient(const int component) { return m_DeformationGradient[component]; }

  void ComputeDeformationGradient(const float triangleXs[][3], const float U[]);
  void ComputeH(float (*p_hs)[3]);
  void ComputeClampedOrSymmetryH(float (*p_hs)[3], const float U[], const float tr);

  void ComputeMembraneForces(float *p_fU, const float sigma[]) const;
  void ComputeBendingForces(float *p_fU, const float M[], const float hs[][3]) const;
  void ComputeClampedBendingForces(float *p_fU, const float M[], const float hs[][3]);
  /** @} */

  /**
   * \name Element Geometry and Topology
   * @{
   */
private:
  const tledElementShellBSTP1 *mpc_Owner;
  const float *mc_X0s;

  int m_EdgePatchNodeIndices[4];
  int m_NumPatchNodes;
  int m_ElementEdgeIndex;
  /* Only used with clamped edges */
  float m_L0, m_EdgeNormal[2];  
  float m_PhiN0[3];

protected:
  /** Returns the local index assigned to the triangle element being registered. */
  void RegisterPrimaryTriangle(const tledElementShellBSTP1 &te, const int localEdgeIndex);

  /** Returns the local index assigned to the triangle element being registered. */
  void RegisterSecondaryTriangle(const tledElementShellBSTP1 &te, const int localEdgeIndex);

  bool IsComplete(void) const { return m_NumPatchNodes == 4; }
  bool IsClamped(void) const { return m_L0 == m_L0; }  

  /** Clamping requires that the nodes are also fixed through standard essential boundary conditions. Symmetry edges are marked with this procedure, too. */
  void ClampEdge(void);

public:
  void ComputeShearAngleMass(float *p_mass);
  void InitPatchTopologyAndGeometry(const tledSurfaceTopology<Surface> &surfTopo, const int edgeIndex);

  /** Last stage of initialisation, initialises Jacobians, etc. after geometry and topology has been determined. */
  void InitTensors(const float elementBasis[][3]);
  /** @} */

#ifdef _GPU_
  /**
   * \name GPU Execution
   * @{
   */
public:
  struct GPUElement {
    /** 
     * For complete elements NeighbourNodeIndex holds the node index of the node of the adjacent element associated with this edge. 
     * If the element lies at a boundary this field is set to a negative value and can be used for identify the element as one lying at a boundary.
     */
    int NeighbourNodeIndex;
    float3 DeformationGradient[2];
    float2 PatchShapeGradients[3];

    union {
      float2 NeighbourPatchShapeGradient;
      /** PhiN0: Boundary normal. Only used in connjunction with symmetry/clamped boundaries. */
      float3 PhiN0;
    };

    union {
      /* Clamped edges don't have neighbours...*/
      float3 NeighbourX0;
      float3 EdgeNormalAndLength0;
    };
  };

public:
  void InitGPU(GPUElement &r_dst);
  /** @} */
#endif

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledElementShellBSTP1Edge(void);
  /** @} */
};
#endif
