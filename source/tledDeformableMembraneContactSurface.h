// =========================================================================
// File:       tledDeformableMembraneContactSurface.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    October 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledDeformableMembraneContactSurface_H
#define tledDeformableMembraneContactSurface_H

#include "tledDeformableContactSurface.h"
#include "tledDeformableMembraneContactSurfaceXMLExporter.h"
#include "tledShellMesh.h"

#include <limits>
#include <cassert>
#include <algorithm>
#include <vector>

/**
 * \brief Untemplated API for contact surfaces comprising a membrane component
 * \ingroup contact
 */
class tledDeformableMembraneContactSurface {
  /**
   * \name Membrane Specific API
   * @{
   */
private:
  int m_MembraneNodeMinIndex, m_MembraneFacetMinIndex, m_MembraneEdgeMinIndex;
  int m_NumMembraneNodes, m_NumMembraneFacets, m_NumMembraneEdges;

public:
  /** Index of first membrane facet in contact surface */
  int GetMembraneFacetBaseIndex(void) const { return m_MembraneFacetMinIndex; }
  void SetMembraneFacetBaseIndex(const int fInd) { m_MembraneFacetMinIndex = fInd; }

  /** Number of membrane elements comprised in mesh. The number of facet objects is twice this number as there's one facet for each side of the membrane element */
  int GetNumberOfMembraneElements(void) const { return m_NumMembraneFacets; }
  void SetNumberOfMembraneElements(const int numEls) { m_NumMembraneFacets = numEls; }

  /** Index of first membrane edge */
  int GetMembraneEdgeBaseIndex(void) const { return m_MembraneEdgeMinIndex; }
  void SetMembraneEdgeBaseIndex(const int eInd) { m_MembraneEdgeMinIndex = eInd; }
  
  /** Number of membrane edges (single-sided) */
  int GetNumberOfMembraneEdges(void) const { return m_NumMembraneEdges; }
  void SetNumberOfMembraneEdges(const int numEdges) { m_NumMembraneEdges = numEdges; }

  /** Index of first membrane node */
  int GetMembraneNodeBaseIndex(void) const { return m_MembraneNodeMinIndex; }
  void SetMembraneNodeBaseIndex(const int nInd) { m_MembraneNodeMinIndex = nInd; }

  /** Number of membrane nodes (single-sided) */
  int GetNumberOfMembraneNodes(void) const { return m_NumMembraneNodes; }
  void SetNumberOfMembraneNodes(const int numNodes) { m_NumMembraneNodes = numNodes; }
  /** @} */

  /**
   * \name Construction, Destruction, Instantiation
   * @{
   */
protected:
  tledDeformableMembraneContactSurface(void) {}

public:
  static tledDeformableContactSurface* CreateSurface(const tledMesh &volumeMesh, const tledSurface &membraneMesh, const bool useGPU);
  static tledDeformableContactSurface* CreateSurface(const XMLNode xmlRep, const bool useGPU);
  static tledDeformableContactSurface* CreateSurface(const std::string &type, const bool useGPU);

  virtual ~tledDeformableMembraneContactSurface(void) {}
  /** @} */
};

/**
 * \brief Membrane contact surface implementation
 * \ingroup contact
 */
template <const int t_numFacetVertices, class TAPI = tledDeformableMembraneContactSurface>
class tledDeformableMembraneContactSurfaceImpl : public tledDeformableContactSurfaceImpl<t_numFacetVertices, TAPI> {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledDeformableContactSurfaceImpl<t_numFacetVertices, TAPI> Superclass;
  typedef typename Superclass::Facet Facet;
  /** @} */
  
  /**
   * \name XML Export
   * @{
   */
public:
  virtual XMLNode ExportToXML(void) const;
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
protected:
  void ConstructFromSolidMeshAndMembrane(const tledMesh &volMesh, const tledSurface &membrane);

public:
  tledDeformableMembraneContactSurfaceImpl(void) : Superclass() {}
  virtual ~tledDeformableMembraneContactSurfaceImpl(void) {}
  /** @} */
};

#include "tledDeformableMembraneContactSurface.tpp"
#endif
