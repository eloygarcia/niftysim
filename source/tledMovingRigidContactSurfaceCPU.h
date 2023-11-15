// =========================================================================
// File:       tledMovingRigidContactSurfaceCPU.h
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
#ifndef tledMovingRigidContactSurfaceCPU_H
#define tledMovingRigidContactSurfaceCPU_H

#include "tledMovingRigidContactSurface.h"
#include "tledRigidContactSurfaceCPU.h"
#include "tledContactSurfaceCPU.h"

#include <vector>
#include <algorithm>
#include <cmath>

class tledMovingRigidContactSurfaceCPU : public tledMovingRigidContactSurface, public tledRigidContactSurfaceCPU {
  /**
   * \name Transform
   * @{
   */
protected:
  static void RotateVector(float *p_n, const float rc[], const float rs[]);
  void ApplyRotationsAndTranslationToNodeList(float *p_nodes, const float rc[], const float rs[], const float tPCor[]);
  void ApplyTranslationToNodeList(float *p_nodes, const float t[]);
  void ApplyRotationsToNormals(float *p_normals, const float n0s[], const int numNorms, const float rc[], const float rs[]);  
  /** @} */

  /**
   * \name Nodes
   * @{
   */
private:
  std::vector<float> m_OldNodeCoordinates, m_Nodes0;
  std::vector<float> m_NodeNormals0;

public:
  void SetNumberOfNodes(const int numNodes);

  virtual void SetAllNodeCoordinates0(const float nodes[]);
  const float* GetAllNodeCoordinates0(void) const { return &m_Nodes0.front(); } 
  float* GetAllNodeCoordinates0(void) { return &m_Nodes0.front(); } 

  virtual void SetAllOldNodeCoordinates(const float nodes[]);
  const float* GetAllOldNodeCoordinates(void) const { return &m_OldNodeCoordinates.front(); }  
  float* GetAllOldNodeCoordinates(void) { return &m_OldNodeCoordinates.front(); }  
  const float* GetOldNodeCoordinates(const int nIndex) const { return &m_OldNodeCoordinates[3*nIndex]; }

  float* GetAllNodeNormals0(void) { return &m_NodeNormals0.front(); }
  const float* GetAllNodeNormals0(void) const { return &m_NodeNormals0.front(); }
  const float* GetNodeNormal0(const int nInd) const { return this->GetAllNodeNormals0() + 3*nInd; }
  /** @} */

  /**
   * \name Facets
   * @{
   */
private:
  std::vector<float> m_OldFacetNormals, m_FacetNormals0;

public:
  virtual void SetNumberOfFacets(const int numFacets);

  float* GetAllOldFacetNormals(void) { return &m_OldFacetNormals.front(); }
  const float* GetAllOldFacetNormals(void) const { return &m_OldFacetNormals.front(); }

  float* GetAllFacetNormals0(void) { return &m_FacetNormals0.front(); }
  const float* GetAllFacetNormals0(void) const { return &m_FacetNormals0.front(); }

  /** Facet normals corresponding to the configuration given by the old node positions */
  const float* GetOldFacetNormal(const int fInd) const { return this->GetAllOldFacetNormals() + 3*fInd; }

  /** Reference configuration facet normals */
  const float* GetFacetNormal0(const int fInd) const { return this->GetAllFacetNormals0() + 3*fInd; }
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
protected:
  tledMovingRigidContactSurfaceCPU(void) {}

public:
  static tledMovingRigidContactSurfaceCPU* CreateSurface(const std::string &type);
  virtual ~tledMovingRigidContactSurfaceCPU(void) {}
  /** @} */
};

/**
 * \brief Moving rigid surface CPU implementation.
 * \ingroup contact
 */
template <class TBaseSurface> 
class tledMovingRigidContactSurfaceImplCPU : public tledRigidContactSurfaceImplCPU<TBaseSurface> {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledRigidContactSurfaceImplCPU<TBaseSurface> Superclass;
  /** @} */

  /**
   * \name Transform
   * @{
   */
protected:
  /** Requires prior update of node positions */
  void TranslateProjectionOperators(void);

public:
  virtual void Update(void);
  /** @} */

  /**
   * \name Output and Visualisation Untemplated Interface
   * @{
   */
public:
  virtual void ResetNodes(void);
  /** @} */

  /**
   * \name Construction, Destruction
   */
public:
  class MovingSurfaceXMLImporter;

protected:
  virtual void InitNormals(void);

public:
  /** Should be instantiated through tledRigidContactSurface::CreateSurface */
  tledMovingRigidContactSurfaceImplCPU(void) {}

  virtual ~tledMovingRigidContactSurfaceImplCPU(void) {}
  /** @} */
};

template <class TSurface> 
class tledMovingRigidContactSurfaceImplCPU<TSurface>::MovingSurfaceXMLImporter : public Superclass::MovingSurfaceXMLImporter {
protected:
  virtual void ReadBasicMeshSpec(void);

public:
  virtual ~MovingSurfaceXMLImporter(void) {}
};

/**
 * \brief Shorthand for moving rigid triangle surface meshes.
 * \ingroup contact
 */
typedef tledMovingRigidContactSurfaceImplCPU<tledMovingRigidContactSurfaceImpl<3, tledMovingRigidContactSurfaceCPU> > tledMovingRigidContactSurfaceT3CPU;

/**
 * \brief Shorthand for moving rigid quad surface meshes.
 * \ingroup contact
 */
typedef tledMovingRigidContactSurfaceImplCPU<tledMovingRigidContactSurfaceImpl<4, tledMovingRigidContactSurfaceCPU> > tledMovingRigidContactSurfaceQ4CPU;

#include "tledMovingRigidContactSurfaceCPU.tpp"
#endif
