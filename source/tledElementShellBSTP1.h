// =========================================================================
// File:       tledElementShellBSTP1.h
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
#ifndef tledElementShellBSTP1_H
#define tledElementShellBSTP1_H

#include "tledSurfaceTopology.h"
#include "tledElementMembraneImpl.h"
#include "tledElementShellBSTP1Edge.h"

#include <limits>
#include <algorithm>

/**
 * \brief Shell element implementation based on Flores, Onate: "Improvements in the membrane 
 * behaviour of the three node rotation-free BST shell triangle using an assumed strain approach" (2005). 
 * \ingroup shell
 */
class tledElementShellBSTP1 : public tledElementMembraneImpl<3> {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledElementMembraneImpl<3> Superclass;
  friend class tledElementShellBSTP1Edge;
  /** @} */

  /**
   * \name Computation
   * @{
   */
private:
  float m_LinearShapeGradients[3][2];
  float m_CurrentNormal[3];
  float m_ContravariantBasis[2][3];
  float m_Hs[3][3], m_ThicknessRatio;

private:
  void _ComputeCurrentConfiguration(float (*p_deformationGrad)[3], float (*p_nodePos)[3], const float U[]) const;

protected:
  virtual float* ComputeStrain(float *p_dst, const float U[]);
  virtual void ComputeForces(float *p_dst, const float stress[]);
  float GetThicknessRatio(void) const { return m_ThicknessRatio; }

public:
  /** Current configuration element normal: Used by edge elements for curvature computation */
  const float* GetCurrentNormal(void) const { return m_CurrentNormal; }

  /** Returns unnormalised current normal */
  virtual float* ComputeCurrentNormal(float *p_n, const float U[]) const { return tledVectorArithmetic::ScalarDiv(p_n, m_CurrentNormal, m_ThicknessRatio); }

  /** Standard linear shape function gradients */
  const float* GetLinearGradient(const int nodeIndex) const { return m_LinearShapeGradients[nodeIndex]; }

  virtual void ComputeElementForces(float *p_dst, const float U[]);

  virtual float ComputeCurrentThickness(const float U[]) const;
  
  /** See Flores 2005 */
  const float* GetContravariantAxis(const int axisIndex) const { return m_ContravariantBasis[axisIndex]; }
  /** @} */

  /**
   * \name Element Geometry
   * @{
   */
private:
  tledElementShellBSTP1Edge m_Edges[3];
  float m_ElementBasis[3][3];
  const float *mc_X0;

protected:
  void RegisterSecondaryTriangle(const int localEdgeIndex, const int neighbourLocalEdgeIndex, const tledElementShellBSTP1 &neighbour);

public:
  const float* GetElementAxis(const int aInd) const { return m_ElementBasis[aInd]; }

  /** Do not use! This element requires a surface topology object for init. */
  virtual void InitialiseElement(const Surface &surf, const int facetIndex) { tledFatalError("Requires a surface topology to initialise"); }
  /** Element initialisation function, requires a surface topology object unlike other 2D elements. */
  virtual void InitialiseElement(const tledSurfaceTopology<Surface> &surfTopo, const int facetIndex);  
  void InitialiseEdgeElements(std::vector<tledElementShellBSTP1> &r_elements, const tledSurfaceTopology<Surface> &surfTopo, const int facetIndex, const std::vector<int> &surface2elSetIndexMap);
  void ClampNodes(const std::vector<bool> &clampedNodeMask);
  void FinaliseInitialisation(void);
  /** @} */

#ifdef _GPU_
  /**
   * \name GPU On-Device Representation
   * @{
   */
public:
  struct GPUElement : public Superclass::GPUElement {
    tledElementShellBSTP1Edge::GPUElement EdgeElements[3];
    float2 LinearShapeGradients[3];
    float3 H[3];    
    float3 X0s[3];
    float3 ContravariantBasis[2], CurrentNormal;
  };

public:
  virtual void InitGPU(tledElementMembrane::GPUElement &r_dst);
  /** @} */
#endif

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledElementShellBSTP1(void);
  virtual ~tledElementShellBSTP1(void) {}
  /** @} */
};

#endif
