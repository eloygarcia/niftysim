// =========================================================================
// File:       tledRigidContactSurfaceCPU.h
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
#ifndef tledRigidContactSurfaceCPU_H
#define tledRigidContactSurfaceCPU_H

#include "tledRigidContactSurface.h"
#include "tledContactSurfaceCPU.h"

class tledRigidContactSurfaceCPU : public tledRigidContactSurface {
  /**
   * \name Construction, Destruction
   * @{
   */
public:
  static tledRigidContactSurfaceCPU* CreateSurface(const XMLNode &meshSpec);
  virtual ~tledRigidContactSurfaceCPU(void) {}
  /** @} */
};

template <class TBaseSurface>
class tledRigidContactSurfaceImplCPU : public tledContactSurfaceImplCPU<TBaseSurface> {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledContactSurfaceImplCPU<TBaseSurface> Superclass;
  typedef typename Superclass::Facet Facet;
  /** @} */

  /**
   * \name Facet Information
   * @{
   */
protected:
  struct FacetData {
    float Normal[3];
    float ProjectionOperator[12];
  };
  
private:
  std::vector<FacetData> m_FacetData;

protected:
  FacetData& GetFacetData(const int fInd) { return m_FacetData[fInd]; }
  const FacetData& GetFacetData(const int fInd) const { return m_FacetData[fInd]; }

public:
  virtual void SetNumberOfFacets(const int numFacets);
  const float* GetFacetProjectionOperator(const int facetIndex) const { return m_FacetData[facetIndex].ProjectionOperator; }
  bool ProjectOntoFacet(float *p_xi, const float x[], const int facetIndex) const { return Superclass::ProjectOntoFacet(p_xi, x, GetFacetProjectionOperator(facetIndex)); }  
  static bool ProjectOntoFacet(float *p_xi, const float x[], const float projOp[]) { return Superclass::ProjectOntoFacet(p_xi, x, projOp); }  
  /** @} */
  
  /**
   * \name Normals 
   * @{
   */
private:
  std::vector<float> m_NodeNormals;

protected:
  void SetNodeNormal(const int nInd, const float n[]) { std::copy(n, n + 3, m_NodeNormals.begin() + 3*nInd); }

public:
  /** Normalised facet normals. */
  const float* GetFacetNormal(const int facetIndex) const { return m_FacetData[facetIndex].Normal; }

  /** Normalised node normals */
  const float* GetNodeNormal(const int nodeIndex) const { return &m_NodeNormals.front() + 3*nodeIndex; }

  const float* GetAllNodeNormals(void) const { return &m_NodeNormals.front(); }
  float* GetAllNodeNormals(void) { return &m_NodeNormals.front(); }
  void SetAllNodeNormals(const float normals[]);
  
  /** Computes the continuous normal at an arbitrary position (given as shape function values) inside a facet. */
  float* ComputeC0Normal(float *p_dstNormal, const int facet[], const float shapeValues[]) const;  
  /** @} */

  /**
   * \name Nodes
   * @{
   */
public:
  virtual void SetNumberOfNodes(const int numNodes);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
protected:
  virtual void InitNormals(void);

public:
  /** Preferably instantiated through tledRigidContactSurface::CreateSurface */
  tledRigidContactSurfaceImplCPU(void) {}

  virtual ~tledRigidContactSurfaceImplCPU(void) {}
  /** @} */
};

/**
 * \brief Shorthand for rigid triangle surface meshes.
 * \ingroup contact
 */
typedef tledRigidContactSurfaceImplCPU<tledRigidContactSurfaceImpl<3, tledRigidContactSurfaceCPU> > tledRigidContactSurfaceT3CPU;

/**
 * \brief Shorthand for rigid quad surface meshes.
 * \ingroup contact
 */
typedef tledRigidContactSurfaceImplCPU<tledRigidContactSurfaceImpl<4, tledRigidContactSurfaceCPU> > tledRigidContactSurfaceQ4CPU;

#include "tledRigidContactSurfaceCPU.tpp"
#endif
