// =========================================================================
// File:       tledRigidContactSurface.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    April 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledRigidContactSurface_H
#define tledRigidContactSurface_H
#include "xmlParser.h"
#include "tledContactSurface.h"
#include "tledVTKSurfaceLoader.h"
#include "tledXMLSurfaceCreator.h"
#include "tledContactSurfaceCreator.h"
#include "tledSurfaceTopology.h"
#include "tledVectorArithmetic.h"

#include <vector>
#include <string>
#include <cassert>

/**
 * \brief Spatially fixed, rigid surfaces
 * \ingroup contact
 */
class tledRigidContactSurface : public tledContactSurface {
  /**
   * \name Slave Nodes
   * @{
   */
private: 
  std::vector<unsigned char> m_SlaveMask;

public:
  const bool* GetSlaveNodeMask(void) const { return (const bool*)&m_SlaveMask.front(); }
  void SetSlaveNodeIndices(const std::vector<int> &slaveNodeIndices, const int numDeformableNodes);
  /** @} */

  /**
   * \name Surface Type Queries
   * @{
   */
public:
  virtual bool IsMoving(void) const { return false; }
  /** @} */

  /**
   * \name Output and Visualisation Untemplated Interface
   * @{
   */
public:
  /** 
   * Returns the final displacement of moving surfaces (IsMoving = true), does nothing in the case of static ones. 
   * May require reset of the surface prior to use (ResetNodes).
   */
  virtual float* GetFinalDisplacement(float *p_dst) const;

  /** Returns the surface to its initial configuration (nodes only) */
  virtual void ResetNodes(void) {}
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
protected:
  tledRigidContactSurface(void) {}

  static bool IsVTKSurfaceXMLRepresentation(const XMLNode rootNode);
  static bool IsMovingSurfaceXMLRepresentation(const XMLNode rootNode);
  static std::string ExtractSurfaceType(const XMLNode rootNode);

public:
  /** To be called at the end of the surface definition process, precomputes secondary information such as surface normals. */
  virtual void Init(void) = 0;

  static tledRigidContactSurface* CreateSurface(const XMLNode &meshSpec, const bool useGPU);
  virtual ~tledRigidContactSurface(void) {}
  /** @} */
};

/**
 * \brief Rigid surface implementation, templated wrt. facet type
 * \ingroup contact
 */
template <const int t_numFacetVertices, class TAPI = tledRigidContactSurface>
class tledRigidContactSurfaceImpl : public tledContactSurfaceImpl<t_numFacetVertices, TAPI> {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledContactSurfaceImpl<t_numFacetVertices, TAPI> Superclass;
  typedef typename Superclass::Facet Facet;
  /** @} */

  /**
   * \name Facet Information
   * @{
   */
public:
  virtual void SetNumberOfFacets(const int numFacets);
  /** @} */
  
  /**
   * \name Construction, Destruction
   * @{
   */
public:
  class BasicSurfaceXMLImporter;

protected:
  void ComputeNodeNormals(float *p_dst) const;
  virtual void InitNormals(void) = 0;

public:
  virtual void Init(void) { this->InitNormals(); }

  /** Preferably instantiated through tledRigidContactSurface::CreateSurface */
  tledRigidContactSurfaceImpl(void) {}

  virtual ~tledRigidContactSurfaceImpl(void) {}
  /** @} */
};

template <const int t_numFacetVertices, class TAPI>
class tledRigidContactSurfaceImpl<t_numFacetVertices, TAPI>::BasicSurfaceXMLImporter : public tledXMLImporter<tledRigidContactSurfaceImpl> {
protected:
  virtual void ReadBasicMeshSpec(void);

public:  
  virtual void Import(void);

public:
  virtual ~BasicSurfaceXMLImporter(void) {}
};

#include "tledRigidContactSurface.tpp"
#endif
