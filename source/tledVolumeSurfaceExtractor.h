// =========================================================================
// File:       tledVolumeSurfaceExtractor.h
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
#ifndef tledVolumeSurfaceExtractor_H
#define tledVolumeSurfaceExtractor_H

#include "tledBasicSurfaceCreator.h"
#include "tledMeshTopology.h"
#include "tledVectorArithmetic.h"

#include <algorithm>
#include <vector>

/**
 * \brief Extracts the surface of a solid mesh
 * \ingroup surface
 */
template <class TSurface>
class tledVolumeSurfaceExtractor : public tledBasicSurfaceCreator<TSurface> {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledBasicSurfaceCreator<TSurface> Superclass;
  /** @} */

  /**
   * \name Topology
   * @{
   */
public:
  static const int NumberOfElementVertices = TSurface::Facet::NumberOfVertices == 3? 4 : 8;

private:
  const tledMesh *mpc_Mesh;
  tledMeshTopology<NumberOfElementVertices> m_Topology;
  std::vector<int> m_SurfaceFacetIndices;

public:
  const tledMesh& GetMesh(void) const { return *mpc_Mesh; }

  /** 
   * \brief Returns the topology using which surface was created. 
   *
   * Only accessible after Create.
   */
  const tledMeshTopology<NumberOfElementVertices>& GetMeshTopology(void) const { return m_Topology; }
  
  /** 
   * \brief Returns the topology using which surface was created (R/W). 
   *
   * Only accessible after Create.
   */
  tledMeshTopology<NumberOfElementVertices>& GetMeshTopology(void) { return m_Topology; }

  /**
   * \brief List of facets constituting the mesh's surface.
   */
  const std::vector<int>& GetSurfaceFacetIndices(void) const { return m_SurfaceFacetIndices; }
  /** @} */

  /**
   * \name Main 
   * @{
   */
protected:
  /** Copies the facet definitions from the mesh topology object */
  void InitFacets();  
  
  /** Sets the node container on the output mesh, is surface type-specific. */
  virtual void InitNodes() = 0;

  /** Orients the output surface's facets such that normals point outward wrt. the solid mesh. */
  void OrientFacets();

public:
  virtual void Create();
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  /** The mesh isn't read until Create is called. */
  tledVolumeSurfaceExtractor(const tledMesh &mesh) : mpc_Mesh(&mesh), m_Topology(mesh) {}

  virtual ~tledVolumeSurfaceExtractor(void) {}
  /** @} */
};

#include "tledVolumeSurfaceExtractor.tpp"

#endif
